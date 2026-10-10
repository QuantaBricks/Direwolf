! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Fixed-density ("direct"/Pulay) second derivative of the semilocal XC

module mod_xc_hess
use MOL_info
use GRID_info
use mod_integrals, only: engine_puream
use mod_xc, only: xc_eval, xc_eval_fxc, xc_eval_fxc_mgga, xc_eval_fxc_mgga_pol, xc_uses_tau
use mod_fixq, only: to_q, from_q
implicit none
INCLUDE 'parameter.h'
private
public :: xc_hess_direct
public :: xc_fxc_apply
public :: xc_fxc_project_vo
public :: xc_fxc_project, xc_fxc_c, xc_fxc_set_loc, xc_fxc_clear_loc
public :: xc_vxc_deriv_atom, xc_vxc_deriv_occ
public :: xc_hess_direct_u, xc_vxc_deriv_occ_u, xc_fxc_apply_u, xc_fxc_c_u
public :: xc_resp_grid_build, xc_resp_grid_clear
public :: eval_bf_derivs

integer,parameter :: MAX_TERM = 6
integer,parameter :: MAX_LPOW = 7
integer :: nbf
integer,allocatable,save :: bf_atom(:),bf_shell(:),bf_nterm(:),bf_idim(:,:)
integer,allocatable,save :: bf_ta(:,:),bf_tb(:,:),bf_tc(:,:)
integer,allocatable,save :: bf_start(:),bf_len(:)
real(8),allocatable,save :: bf_coef(:,:)
logical,save :: tables_built = .false.
integer,save :: nrg = 0
real(8),save,public :: fxu_t(3) = 0.0d0
logical,save :: loc_ready = .false.
integer,save :: loc_n = 0, loc_nocc = 0
real(8),allocatable,save :: loc_Co(:,:), loc_Cl(:,:), loc_Q(:,:)
logical,save :: sh_built = .false.
integer,save :: nsh = 0
integer,allocatable,save :: sh_atom(:), sh_loc(:), sh_first(:), sh_dim(:)
real(8),allocatable,save :: sh_rc(:)
real(8),allocatable,save :: rg_coor(:,:), rg_w(:)

integer,save :: vxc_cache_nConts = -1, vxc_cache_natoms = -1
real(8),allocatable,save :: vxc_cache_P(:,:),vxc_cache(:,:,:,:)

interface
   subroutine GTOeval_shell_batch(i,j,coor_batch,nb,off,ld,val0_out,val1_out,val2_out)
   implicit none
   integer,intent(in) :: i,j,nb,off,ld
   real(8),intent(in) :: coor_batch(3,nb)
   real(8),intent(inout) :: val0_out(ld,nb), val1_out(ld,3,nb)
   real(8),intent(inout),optional :: val2_out(ld,3,3,nb)
   end subroutine GTOeval_shell_batch
end interface

contains

subroutine xc_kern_c(mg, rho, sig, tau, o)
implicit none
logical,intent(in) :: mg
real(8),intent(in) :: rho, sig, tau
real(8),intent(out) :: o(9)
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1),ta(1),tb(1),vta(1),vtb(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)
real(8) :: g1(1),g2(1),g3(1),g4(1),g5(1),g6(1),r1(1),s1(1),t1(1)
o = 0.0d0
ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
saa(1)=sig; sbb(1)=sig; sab(1)=sig; st(1)=sig
if (mg) then
   ta(1) = 0.5d0*tau ; tb(1) = 0.5d0*tau
   call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                vra,vrb,vsaa,vsbb,vsab,ta,tb,vta,vtb)
   r1(1) = rho ; s1(1) = sig ; t1(1) = tau
   call xc_eval_fxc_mgga(1,r1,s1,t1,g1,g2,g3,g4,g5,g6)
   o(1) = d1r(1) ; o(2) = d1s(1) ; o(3) = g1(1) ; o(4) = g2(1) ; o(5) = g3(1)
   o(6) = vta(1) ; o(7) = g4(1) ; o(8) = g5(1) ; o(9) = g6(1)
else
   call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                vra,vrb,vsaa,vsbb,vsab)
   call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
   o(1) = d1r(1) ; o(2) = d1s(1) ; o(3) = v2r2(1) ; o(4) = v2rs(1) ; o(5) = v2s2(1)
endif
end subroutine xc_kern_c

subroutine build_tables()
implicit none
integer :: i,j,l,idim0,a,b,c
if (tables_built) return
nbf = nconts
allocate(bf_atom(nbf),bf_shell(nbf),bf_nterm(nbf))
allocate(bf_idim(MAX_TERM,nbf),bf_ta(MAX_TERM,nbf),bf_tb(MAX_TERM,nbf),bf_tc(MAX_TERM,nbf))
allocate(bf_coef(MAX_TERM,nbf))
allocate(bf_start(natoms),bf_len(natoms))
l = 0
do i = 1,natoms
   do j = 1,size(atoms(i)%shell)
      if (engine_puream .and. atoms(i)%shell(j)%angMoment .ge. 5) then
         print *, 'ERROR: xc_hess puream Lt=',atoms(i)%shell(j)%angMoment,' not implemented'
         stop 1
      else if (engine_puream .and. atoms(i)%shell(j)%angMoment .eq. 4) then
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_ta(1,l)=3;bf_tb(1,l)=1;bf_tc(1,l)=0;bf_idim(1,l)=2;bf_coef(1,l)=sqrt(5.0d0)/2.0d0
         bf_ta(2,l)=1;bf_tb(2,l)=3;bf_tc(2,l)=0;bf_idim(2,l)=7;bf_coef(2,l)=-sqrt(5.0d0)/2.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_ta(1,l)=2;bf_tb(1,l)=1;bf_tc(1,l)=1;bf_idim(1,l)=5;bf_coef(1,l)=3.0d0*sqrt(2.0d0)/4.0d0
         bf_ta(2,l)=0;bf_tb(2,l)=3;bf_tc(2,l)=1;bf_idim(2,l)=12;bf_coef(2,l)=-sqrt(10.0d0)/4.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=1;bf_tb(1,l)=1;bf_tc(1,l)=2;bf_idim(1,l)=9;bf_coef(1,l)=3.0d0*sqrt(7.0d0)/7.0d0
         bf_ta(2,l)=3;bf_tb(2,l)=1;bf_tc(2,l)=0;bf_idim(2,l)=2;bf_coef(2,l)=-sqrt(35.0d0)/14.0d0
         bf_ta(3,l)=1;bf_tb(3,l)=3;bf_tc(3,l)=0;bf_idim(3,l)=7;bf_coef(3,l)=-sqrt(35.0d0)/14.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=0;bf_tb(1,l)=1;bf_tc(1,l)=3;bf_idim(1,l)=14;bf_coef(1,l)=sqrt(70.0d0)/7.0d0
         bf_ta(2,l)=2;bf_tb(2,l)=1;bf_tc(2,l)=1;bf_idim(2,l)=5;bf_coef(2,l)=-3.0d0*sqrt(14.0d0)/28.0d0
         bf_ta(3,l)=0;bf_tb(3,l)=3;bf_tc(3,l)=1;bf_idim(3,l)=12;bf_coef(3,l)=-3.0d0*sqrt(70.0d0)/28.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=6
         bf_ta(1,l)=0;bf_tb(1,l)=0;bf_tc(1,l)=4;bf_idim(1,l)=15;bf_coef(1,l)=1.0d0
         bf_ta(2,l)=2;bf_tb(2,l)=0;bf_tc(2,l)=2;bf_idim(2,l)=6;bf_coef(2,l)=-3.0d0*sqrt(105.0d0)/35.0d0
         bf_ta(3,l)=0;bf_tb(3,l)=2;bf_tc(3,l)=2;bf_idim(3,l)=13;bf_coef(3,l)=-3.0d0*sqrt(105.0d0)/35.0d0
         bf_ta(4,l)=4;bf_tb(4,l)=0;bf_tc(4,l)=0;bf_idim(4,l)=1;bf_coef(4,l)=3.0d0/8.0d0
         bf_ta(5,l)=2;bf_tb(5,l)=2;bf_tc(5,l)=0;bf_idim(5,l)=4;bf_coef(5,l)=3.0d0*sqrt(105.0d0)/140.0d0
         bf_ta(6,l)=0;bf_tb(6,l)=4;bf_tc(6,l)=0;bf_idim(6,l)=11;bf_coef(6,l)=3.0d0/8.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=1;bf_tb(1,l)=0;bf_tc(1,l)=3;bf_idim(1,l)=10;bf_coef(1,l)=sqrt(70.0d0)/7.0d0
         bf_ta(2,l)=3;bf_tb(2,l)=0;bf_tc(2,l)=1;bf_idim(2,l)=3;bf_coef(2,l)=-3.0d0*sqrt(70.0d0)/28.0d0
         bf_ta(3,l)=1;bf_tb(3,l)=2;bf_tc(3,l)=1;bf_idim(3,l)=8;bf_coef(3,l)=-3.0d0*sqrt(14.0d0)/28.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=4
         bf_ta(1,l)=2;bf_tb(1,l)=0;bf_tc(1,l)=2;bf_idim(1,l)=6;bf_coef(1,l)=3.0d0*sqrt(21.0d0)/14.0d0
         bf_ta(2,l)=0;bf_tb(2,l)=2;bf_tc(2,l)=2;bf_idim(2,l)=13;bf_coef(2,l)=-3.0d0*sqrt(21.0d0)/14.0d0
         bf_ta(3,l)=4;bf_tb(3,l)=0;bf_tc(3,l)=0;bf_idim(3,l)=1;bf_coef(3,l)=-sqrt(5.0d0)/4.0d0
         bf_ta(4,l)=0;bf_tb(4,l)=4;bf_tc(4,l)=0;bf_idim(4,l)=11;bf_coef(4,l)=sqrt(5.0d0)/4.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_ta(1,l)=3;bf_tb(1,l)=0;bf_tc(1,l)=1;bf_idim(1,l)=3;bf_coef(1,l)=sqrt(10.0d0)/4.0d0
         bf_ta(2,l)=1;bf_tb(2,l)=2;bf_tc(2,l)=1;bf_idim(2,l)=8;bf_coef(2,l)=-3.0d0*sqrt(2.0d0)/4.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=4;bf_tb(1,l)=0;bf_tc(1,l)=0;bf_idim(1,l)=1;bf_coef(1,l)=sqrt(35.0d0)/8.0d0
         bf_ta(2,l)=2;bf_tb(2,l)=2;bf_tc(2,l)=0;bf_idim(2,l)=4;bf_coef(2,l)=-3.0d0*sqrt(3.0d0)/4.0d0
         bf_ta(3,l)=0;bf_tb(3,l)=4;bf_tc(3,l)=0;bf_idim(3,l)=11;bf_coef(3,l)=sqrt(35.0d0)/8.0d0
      else if (engine_puream .and. atoms(i)%shell(j)%angMoment .eq. 3) then
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_ta(1,l)=2;bf_tb(1,l)=1;bf_tc(1,l)=0;bf_idim(1,l)=2;bf_coef(1,l)=3.0d0*sqrt(2.0d0)/4.0d0
         bf_ta(2,l)=0;bf_tb(2,l)=3;bf_tc(2,l)=0;bf_idim(2,l)=7;bf_coef(2,l)=-sqrt(10.0d0)/4.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_ta(1,l)=1;bf_tb(1,l)=1;bf_tc(1,l)=1;bf_idim(1,l)=5;bf_coef(1,l)=1.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=0;bf_tb(1,l)=1;bf_tc(1,l)=2;bf_idim(1,l)=9;bf_coef(1,l)=sqrt(30.0d0)/5.0d0
         bf_ta(2,l)=2;bf_tb(2,l)=1;bf_tc(2,l)=0;bf_idim(2,l)=2;bf_coef(2,l)=-sqrt(30.0d0)/20.0d0
         bf_ta(3,l)=0;bf_tb(3,l)=3;bf_tc(3,l)=0;bf_idim(3,l)=7;bf_coef(3,l)=-sqrt(6.0d0)/4.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=0;bf_tb(1,l)=0;bf_tc(1,l)=3;bf_idim(1,l)=10;bf_coef(1,l)=1.0d0
         bf_ta(2,l)=2;bf_tb(2,l)=0;bf_tc(2,l)=1;bf_idim(2,l)=3;bf_coef(2,l)=-3.0d0*sqrt(5.0d0)/10.0d0
         bf_ta(3,l)=0;bf_tb(3,l)=2;bf_tc(3,l)=1;bf_idim(3,l)=8;bf_coef(3,l)=-3.0d0*sqrt(5.0d0)/10.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=1;bf_tb(1,l)=0;bf_tc(1,l)=2;bf_idim(1,l)=6;bf_coef(1,l)=sqrt(30.0d0)/5.0d0
         bf_ta(2,l)=3;bf_tb(2,l)=0;bf_tc(2,l)=0;bf_idim(2,l)=1;bf_coef(2,l)=-sqrt(6.0d0)/4.0d0
         bf_ta(3,l)=1;bf_tb(3,l)=2;bf_tc(3,l)=0;bf_idim(3,l)=4;bf_coef(3,l)=-sqrt(30.0d0)/20.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_ta(1,l)=2;bf_tb(1,l)=0;bf_tc(1,l)=1;bf_idim(1,l)=3;bf_coef(1,l)=sqrt(3.0d0)/2.0d0
         bf_ta(2,l)=0;bf_tb(2,l)=2;bf_tc(2,l)=1;bf_idim(2,l)=8;bf_coef(2,l)=-sqrt(3.0d0)/2.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_ta(1,l)=3;bf_tb(1,l)=0;bf_tc(1,l)=0;bf_idim(1,l)=1;bf_coef(1,l)=sqrt(10.0d0)/4.0d0
         bf_ta(2,l)=1;bf_tb(2,l)=2;bf_tc(2,l)=0;bf_idim(2,l)=4;bf_coef(2,l)=-3.0d0*sqrt(2.0d0)/4.0d0
      else if (engine_puream .and. atoms(i)%shell(j)%angMoment .eq. 2) then
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_ta(1,l)=1;bf_tb(1,l)=1;bf_tc(1,l)=0;bf_idim(1,l)=2;bf_coef(1,l)=1.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_ta(1,l)=0;bf_tb(1,l)=1;bf_tc(1,l)=1;bf_idim(1,l)=5;bf_coef(1,l)=1.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=3
         bf_ta(1,l)=0;bf_tb(1,l)=0;bf_tc(1,l)=2;bf_idim(1,l)=6;bf_coef(1,l)=1.0d0
         bf_ta(2,l)=2;bf_tb(2,l)=0;bf_tc(2,l)=0;bf_idim(2,l)=1;bf_coef(2,l)=-0.5d0
         bf_ta(3,l)=0;bf_tb(3,l)=2;bf_tc(3,l)=0;bf_idim(3,l)=4;bf_coef(3,l)=-0.5d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
         bf_ta(1,l)=1;bf_tb(1,l)=0;bf_tc(1,l)=1;bf_idim(1,l)=3;bf_coef(1,l)=1.0d0
         l=l+1; bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=2
         bf_ta(1,l)=2;bf_tb(1,l)=0;bf_tc(1,l)=0;bf_idim(1,l)=1;bf_coef(1,l)=sqrt(3.0d0)/2.0d0
         bf_ta(2,l)=0;bf_tb(2,l)=2;bf_tc(2,l)=0;bf_idim(2,l)=4;bf_coef(2,l)=-sqrt(3.0d0)/2.0d0
      else
         idim0 = 0
         do a = atoms(i)%shell(j)%angMoment,0,-1
            do b = atoms(i)%shell(j)%angMoment-a,0,-1
               c = atoms(i)%shell(j)%angMoment-a-b
               l = l+1; idim0 = idim0+1
               bf_atom(l)=i; bf_shell(l)=j; bf_nterm(l)=1
               bf_ta(1,l)=a; bf_tb(1,l)=b; bf_tc(1,l)=c
               bf_idim(1,l)=idim0; bf_coef(1,l)=1.0d0
            enddo
         enddo
       endif
    enddo
 enddo
 do i = 1,natoms
    bf_start(i) = 0
    bf_len(i) = 0
 enddo
 do l = 1,nbf
    i = bf_atom(l)
    if (bf_start(i) .eq. 0) bf_start(i) = l
    bf_len(i) = bf_len(i) + 1
 enddo
 tables_built = .true.
 end subroutine build_tables

real(8) function mono_deriv(a,b,c, nx,ny,nz, S, px,py,pz, L)
implicit none
integer,intent(in) :: a,b,c,nx,ny,nz,L
real(8),intent(in) :: S(0:3), px(0:L),py(0:L),pz(0:L)
integer,parameter :: MAXT=16
real(8) :: tc(MAXT)
integer :: ta(MAXT),tb(MAXT),tcc(MAXT),tl(MAXT)
integer :: n,nt,i,d
real(8) :: acc
tc(1)=1.0d0; ta(1)=a; tb(1)=b; tcc(1)=c; tl(1)=0; n=1
do d = 1,nx
   do i = 1,n
      tc(n+i)=tc(i)*ta(i); ta(n+i)=ta(i)-1; tb(n+i)=tb(i); tcc(n+i)=tcc(i); tl(n+i)=tl(i)
      tc(i)=-2.0d0*tc(i); ta(i)=ta(i)+1; tl(i)=tl(i)+1
   enddo
   n=2*n
enddo
do d = 1,ny
   do i = 1,n
      tc(n+i)=tc(i)*tb(i); tb(n+i)=tb(i)-1; ta(n+i)=ta(i); tcc(n+i)=tcc(i); tl(n+i)=tl(i)
      tc(i)=-2.0d0*tc(i); tb(i)=tb(i)+1; tl(i)=tl(i)+1
   enddo
   n=2*n
enddo
do d = 1,nz
   do i = 1,n
      tc(n+i)=tc(i)*tcc(i); tcc(n+i)=tcc(i)-1; ta(n+i)=ta(i); tb(n+i)=tb(i); tl(n+i)=tl(i)
      tc(i)=-2.0d0*tc(i); tcc(i)=tcc(i)+1; tl(i)=tl(i)+1
   enddo
   n=2*n
enddo
acc = 0.0d0
do i = 1,n
   if (tl(i) .le. 3) acc = acc + tc(i)*px(max(ta(i),0))*py(max(tb(i),0))*pz(max(tcc(i),0))*S(tl(i))
enddo
mono_deriv = acc
end function mono_deriv

subroutine eval_bf_derivs(coor, v0,v1,v2,v3, maxder)
implicit none
real(8),intent(in) :: coor(3)
real(8),intent(out) :: v0(nbf), v1(nbf,3), v2(nbf,3,3), v3(nbf,3,3,3)
integer,intent(in),optional :: maxder
integer,parameter :: NGMAX = 256
integer :: i,j,l,kk,gi,nG,idim,a,b,c,md
real(8) :: CR(3),dist2,E(NGMAX),S(0:3)
real(8) :: px(0:MAX_LPOW),py(0:MAX_LPOW),pz(0:MAX_LPOW)
real(8),parameter :: Z0=0.0d0
call build_tables()
md = 3
if (present(maxder)) md = maxder
v0=0; v1=0; v2=0; v3=0
do i = 1,natoms
   CR = coor - atoms(i)%coor*ans2bohr
   dist2 = sum(CR**2)
   px(0)=1.0d0; py(0)=1.0d0; pz(0)=1.0d0
   do l = 1,MAX_LPOW
      px(l)=px(l-1)*CR(1); py(l)=py(l-1)*CR(2); pz(l)=pz(l-1)*CR(3)
   enddo
   do l = 1,nbf
      if (bf_atom(l) .ne. i) cycle
      j = bf_shell(l)
      nG = atoms(i)%shell(j)%nGauss
      if (nG .gt. NGMAX) then
         print *,'ERROR: xc_hess nGauss > NGMAX'; stop 1
      endif
      do gi = 1,nG
         E(gi) = DEXP(-atoms(i)%shell(j)%exponents(gi)*dist2)
      enddo
      do kk = 1,bf_nterm(l)
         idim = bf_idim(kk,l)
         S = 0.0d0
         do gi = 1,nG
            S(0) = S(0) + atoms(i)%shell(j)%cnVal(gi,idim)*E(gi)
            S(1) = S(1) + atoms(i)%shell(j)%exponents(gi)*atoms(i)%shell(j)%cnVal(gi,idim)*E(gi)
            S(2) = S(2) + atoms(i)%shell(j)%exponents(gi)**2*atoms(i)%shell(j)%cnVal(gi,idim)*E(gi)
            S(3) = S(3) + atoms(i)%shell(j)%exponents(gi)**3*atoms(i)%shell(j)%cnVal(gi,idim)*E(gi)
         enddo
          a=bf_ta(kk,l); b=bf_tb(kk,l); c=bf_tc(kk,l)
          v0(l) = v0(l) + bf_coef(kk,l)*mono_deriv(a,b,c,0,0,0,S,px,py,pz,MAX_LPOW)
          v1(l,1) = v1(l,1) + bf_coef(kk,l)*mono_deriv(a,b,c,1,0,0,S,px,py,pz,MAX_LPOW)
          v1(l,2) = v1(l,2) + bf_coef(kk,l)*mono_deriv(a,b,c,0,1,0,S,px,py,pz,MAX_LPOW)
          v1(l,3) = v1(l,3) + bf_coef(kk,l)*mono_deriv(a,b,c,0,0,1,S,px,py,pz,MAX_LPOW)
          if (md .ge. 2) then
          v2(l,1,1) = v2(l,1,1) + bf_coef(kk,l)*mono_deriv(a,b,c,2,0,0,S,px,py,pz,MAX_LPOW)
          v2(l,2,2) = v2(l,2,2) + bf_coef(kk,l)*mono_deriv(a,b,c,0,2,0,S,px,py,pz,MAX_LPOW)
          v2(l,3,3) = v2(l,3,3) + bf_coef(kk,l)*mono_deriv(a,b,c,0,0,2,S,px,py,pz,MAX_LPOW)
          v2(l,1,2) = v2(l,1,2) + bf_coef(kk,l)*mono_deriv(a,b,c,1,1,0,S,px,py,pz,MAX_LPOW)
          v2(l,1,3) = v2(l,1,3) + bf_coef(kk,l)*mono_deriv(a,b,c,1,0,1,S,px,py,pz,MAX_LPOW)
          v2(l,2,3) = v2(l,2,3) + bf_coef(kk,l)*mono_deriv(a,b,c,0,1,1,S,px,py,pz,MAX_LPOW)
          endif
          if (md .ge. 3) then
          v3(l,1,1,1) = v3(l,1,1,1) + bf_coef(kk,l)*mono_deriv(a,b,c,3,0,0,S,px,py,pz,MAX_LPOW)
          v3(l,2,2,2) = v3(l,2,2,2) + bf_coef(kk,l)*mono_deriv(a,b,c,0,3,0,S,px,py,pz,MAX_LPOW)
          v3(l,3,3,3) = v3(l,3,3,3) + bf_coef(kk,l)*mono_deriv(a,b,c,0,0,3,S,px,py,pz,MAX_LPOW)
          v3(l,1,1,2) = v3(l,1,1,2) + bf_coef(kk,l)*mono_deriv(a,b,c,2,1,0,S,px,py,pz,MAX_LPOW)
          v3(l,1,1,3) = v3(l,1,1,3) + bf_coef(kk,l)*mono_deriv(a,b,c,2,0,1,S,px,py,pz,MAX_LPOW)
          v3(l,1,2,2) = v3(l,1,2,2) + bf_coef(kk,l)*mono_deriv(a,b,c,1,2,0,S,px,py,pz,MAX_LPOW)
          v3(l,2,2,3) = v3(l,2,2,3) + bf_coef(kk,l)*mono_deriv(a,b,c,0,2,1,S,px,py,pz,MAX_LPOW)
          v3(l,1,3,3) = v3(l,1,3,3) + bf_coef(kk,l)*mono_deriv(a,b,c,1,0,2,S,px,py,pz,MAX_LPOW)
          v3(l,2,3,3) = v3(l,2,3,3) + bf_coef(kk,l)*mono_deriv(a,b,c,0,1,2,S,px,py,pz,MAX_LPOW)
          v3(l,1,2,3) = v3(l,1,2,3) + bf_coef(kk,l)*mono_deriv(a,b,c,1,1,1,S,px,py,pz,MAX_LPOW)
          endif
       enddo
    enddo
 enddo
 do l = 1,nbf
    if (md .ge. 2) then
       v2(l,2,1)=v2(l,1,2); v2(l,3,1)=v2(l,1,3); v2(l,3,2)=v2(l,2,3)
    endif
    if (md .ge. 3) call fill_v3(v3(l,:,:,:))
 enddo
 end subroutine eval_bf_derivs

subroutine fill_v3(t)
implicit none
real(8),intent(inout) :: t(3,3,3)
real(8) :: u(3,3,3)
integer :: i,j,k,a,b,c
u = t
do i = 1,3
   do j = 1,3
      do k = 1,3
         a = min(i,j,k)
         c = max(i,j,k)
         b = i+j+k-a-c
         t(i,j,k) = u(a,b,c)
      enddo
   enddo
enddo
end subroutine fill_v3

subroutine xc_hess_direct_pt(nConts, natoms, Ptot, d2E_xc)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
external :: dgemv,dgemm
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: d2E_xc(3,natoms,3,natoms)

integer :: ig,ia,ib,m,n,dd,x,y
real(8) :: coor(3),w,rho,sig
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: A0(:),A1(:,:),U(:,:,:)
real(8),allocatable :: B1m(:,:),B2m(:,:,:),C2(:,:,:,:),C2b(:,:,:)
real(8),allocatable :: B1m_all(:,:,:),B2m_all(:,:,:,:)
real(8),allocatable :: rhoA(:,:),rhoAB(:,:,:,:)
real(8),allocatable :: gammaA(:,:,:),gammaAB(:,:,:,:,:),sigmaA(:,:),sigmaAB(:,:,:,:)
real(8),allocatable :: d2E_xc_loc(:,:,:,:)
real(8) :: gamma(3)
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vr(1),vs(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)
real(8) :: vrho,vsigma
logical :: ti
character(len=8) :: tienv

call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_hess nConts/=nbf',nConts,nbf; stop 1
endif
d2E_xc = 0.0d0
call get_environment_variable("ENGINE_XC_HESS_TI", tienv)
ti = (trim(tienv) .ne. "0")
!$omp parallel default(shared) &
!$omp&   private(ig,ia,ib,m,n,dd,x,y,coor,w,rho,sig, &
!$omp&           v0,v1,v2,v3,A0,A1,U,B1m,B2m,B1m_all,B2m_all,C2,C2b,rhoA,rhoAB, &
!$omp&           gammaA,gammaAB,sigmaA,sigmaAB,d2E_xc_loc,gamma, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vr,vs,vra,vrb, &
!$omp&           vsaa,vsbb,vsab,v2r2,v2rs,v2s2,vrho,vsigma)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(A0(nConts),A1(nConts,3),U(3,natoms,nConts))
allocate(B1m(nConts,3),B2m(nConts,3,3),C2(nConts,3,3,3),C2b(nConts,3,3))
allocate(B1m_all(nConts,3,natoms),B2m_all(nConts,3,3,natoms))
allocate(rhoA(3,natoms),rhoAB(3,natoms,3,natoms))
allocate(gammaA(3,3,natoms),gammaAB(3,3,3,natoms,natoms))
allocate(sigmaA(3,natoms),sigmaAB(3,natoms,3,natoms))
allocate(d2E_xc_loc(3,natoms,3,natoms))
d2E_xc_loc = 0.0d0

!$omp do schedule(static)
do ig = 1,ngrids
   w = Grids(ig)%weight
   if (w .eq. 0.0d0) cycle
   coor = Grids(ig)%coor
   if (ti) then
      call eval_bf_derivs(coor, v0,v1,v2,v3, 2)
   else
      call eval_bf_derivs(coor, v0,v1,v2,v3)
   endif
   call dgemv('N',nConts,nConts,1.0d0,Ptot,nConts,v0,1,0.0d0,A0,1)
   call dgemm('N','N',nConts,3,nConts,1.0d0,Ptot,nConts,v1,nConts,0.0d0,A1,nConts)
   rho = sum(v0*A0)
   do dd = 1,3
      gamma(dd) = 2.0d0*sum(v1(:,dd)*A0)
   enddo
   sig = sum(gamma**2)

   ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
   saa(1)=sig; sbb(1)=sig; sab(1)=sig; st(1)=sig
   call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                vra,vrb,vsaa,vsbb,vsab)
   call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
   vrho = d1r(1); vsigma = d1s(1)

   U = 0.0d0
   do ia = 1,natoms
      do m = bf_start(ia),bf_start(ia)+bf_len(ia)-1
         do x = 1,3
            U(x,ia,:) = U(x,ia,:) + Ptot(m,:)*v1(m,x)
         enddo
      enddo
   enddo

   do ia = 1,natoms
      do x = 1,3
         rhoA(x,ia) = 0.0d0
         do m = bf_start(ia),bf_start(ia)+bf_len(ia)-1
            rhoA(x,ia) = rhoA(x,ia) - 2.0d0*v1(m,x)*A0(m)
         enddo
      enddo
      do dd = 1,3
         do x = 1,3
            gammaA(dd,x,ia) = 0.0d0
            do m = bf_start(ia),bf_start(ia)+bf_len(ia)-1
               gammaA(dd,x,ia) = gammaA(dd,x,ia) &
                  - 2.0d0*( v2(m,dd,x)*A0(m) + v1(m,x)*A1(m,dd) )
            enddo
         enddo
      enddo
      do x = 1,3
         sigmaA(x,ia) = 2.0d0*sum(gamma*gammaA(:,x,ia))
      enddo
   enddo

   do ib = 1,natoms
      do m = 1,nConts
         do y = 1,3
            B1m_all(m,y,ib) = 0.0d0
         enddo
         do dd = 1,3
            do y = 1,3
               B2m_all(m,dd,y,ib) = 0.0d0
            enddo
         enddo
      enddo
      do n = bf_start(ib),bf_start(ib)+bf_len(ib)-1
         do y = 1,3
            do m = 1,nConts
               B1m_all(m,y,ib) = B1m_all(m,y,ib) + Ptot(m,n)*v1(n,y)
            enddo
         enddo
         do dd = 1,3
            do y = 1,3
               do m = 1,nConts
                  B2m_all(m,dd,y,ib) = B2m_all(m,dd,y,ib) + Ptot(m,n)*v2(n,dd,y)
               enddo
            enddo
         enddo
      enddo
   enddo

   do ia = 1,natoms
      do ib = 1,natoms
         if (ti .and. ia .eq. ib) cycle
         do x = 1,3
            do y = 1,3
               rhoAB(x,ia,y,ib) = 0.0d0
               if (ia .eq. ib) then
                  do m = bf_start(ia),bf_start(ia)+bf_len(ia)-1
                     rhoAB(x,ia,y,ib) = rhoAB(x,ia,y,ib) + 2.0d0*v2(m,x,y)*A0(m)
                  enddo
               endif
               do m = bf_start(ia),bf_start(ia)+bf_len(ia)-1
                  rhoAB(x,ia,y,ib) = rhoAB(x,ia,y,ib) + v1(m,x)*U(y,ib,m)
               enddo
               do n = bf_start(ib),bf_start(ib)+bf_len(ib)-1
                  rhoAB(x,ia,y,ib) = rhoAB(x,ia,y,ib) + v1(n,y)*U(x,ia,n)
               enddo
            enddo
         enddo
         if (ia .ne. ib) then
            do dd = 1,3
               do x = 1,3
                  do y = 1,3
                     gammaAB(dd,x,y,ia,ib) = 0.0d0
                     do m = bf_start(ia),bf_start(ia)+bf_len(ia)-1
                        gammaAB(dd,x,y,ia,ib) = gammaAB(dd,x,y,ia,ib) + 2.0d0*( &
                             v2(m,dd,x)*B1m_all(m,y,ib) &
                           + v1(m,x)*B2m_all(m,dd,y,ib) )
                     enddo
                  enddo
               enddo
            enddo
         else
            C2 = 0.0d0; C2b = 0.0d0
            do n = bf_start(ia),bf_start(ia)+bf_len(ia)-1
               do m = 1,nConts
                  do dd = 1,3
                     do y = 1,3
                        C2b(m,dd,y) = C2b(m,dd,y) + Ptot(m,n)*v2(n,dd,y)
                     enddo
                     do x = 1,3
                        do y = 1,3
                           C2(m,dd,x,y) = C2(m,dd,x,y) + Ptot(m,n)*v3(n,dd,x,y)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
            do dd = 1,3
               do x = 1,3
                  do y = 1,3
                     gammaAB(dd,x,y,ia,ib) = 0.0d0
                     do m = bf_start(ia),bf_start(ia)+bf_len(ia)-1
                        gammaAB(dd,x,y,ia,ib) = gammaAB(dd,x,y,ia,ib) + 2.0d0*( &
                             v3(m,dd,x,y)*A0(m) &
                           + v2(m,dd,x)*B1m_all(m,y,ib) &
                           + v2(m,x,y)*A1(m,dd) &
                           + v1(m,x)*C2b(m,dd,y) )
                     enddo
                  enddo
               enddo
            enddo
         endif
         do x = 1,3
            do y = 1,3
               sigmaAB(x,ia,y,ib) = 0.0d0
               do dd = 1,3
                  sigmaAB(x,ia,y,ib) = sigmaAB(x,ia,y,ib) + 2.0d0*( &
                       gammaA(dd,x,ia)*gammaA(dd,y,ib) &
                     + gamma(dd)*gammaAB(dd,x,y,ia,ib) )
               enddo
            enddo
         enddo
          do x = 1,3
             do y = 1,3
                d2E_xc_loc(x,ia,y,ib) = d2E_xc_loc(x,ia,y,ib) + w*( &
                     vrho*rhoAB(x,ia,y,ib) + vsigma*sigmaAB(x,ia,y,ib) &
                   + v2r2(1)*rhoA(x,ia)*rhoA(y,ib) &
                   + v2rs(1)*(rhoA(x,ia)*sigmaA(y,ib) + sigmaA(x,ia)*rhoA(y,ib)) &
                   + v2s2(1)*sigmaA(x,ia)*sigmaA(y,ib) )
             enddo
          enddo
       enddo
    enddo
 enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) d2E_xc = d2E_xc + d2E_xc_loc
   enddo
end block

deallocate(v0,v1,v2,v3,A0,A1,U,B1m,B2m,B1m_all,B2m_all,C2,C2b,rhoA,rhoAB,gammaA,gammaAB,sigmaA,sigmaAB,d2E_xc_loc)
!$omp end parallel
block
   character(len=8) :: e
   integer :: ja
   if (ti) then
      do ia = 1,natoms
         d2E_xc(:,ia,:,ia) = 0.0d0
         do ja = 1,natoms
            if (ja .ne. ia) d2E_xc(:,ia,:,ia) = d2E_xc(:,ia,:,ia) - d2E_xc(:,ia,:,ja)
         enddo
      enddo
   endif
end block
end subroutine xc_hess_direct_pt

subroutine build_shell_tables()
implicit none
integer :: l, s, i, j
if (sh_built) return
call build_tables()
nsh = 0
do l = 1,nbf
   if (l .eq. 1) then
      nsh = 1
   else if (bf_atom(l) .ne. bf_atom(l-1) .or. bf_shell(l) .ne. bf_shell(l-1)) then
      nsh = nsh + 1
   endif
enddo
allocate(sh_atom(nsh), sh_loc(nsh), sh_first(nsh), sh_dim(nsh), sh_rc(nsh))
sh_first = 0 ; sh_dim = 0
s = 0
do l = 1,nbf
   if (l .eq. 1) then
      s = 1
   else if (bf_atom(l) .ne. bf_atom(l-1) .or. bf_shell(l) .ne. bf_shell(l-1)) then
      s = s + 1
   endif
   if (sh_dim(s) .eq. 0) then
      sh_atom(s) = bf_atom(l) ; sh_loc(s) = bf_shell(l) ; sh_first(s) = l
   endif
   sh_dim(s) = sh_dim(s) + 1
enddo
do s = 1,nsh
   i = sh_atom(s); j = sh_loc(s)
   sh_rc(s) = sqrt(-log(1.0d-7)/minval(atoms(i)%shell(j)%exponents)) + 2.0d0
enddo
sh_built = .true.
end subroutine build_shell_tables

subroutine xc_blk_shells(ng, coor, nsig, gidx, nsl, sl)
implicit none
integer,intent(in) :: ng
real(8),intent(in) :: coor(3,ng)
integer,intent(out) :: nsig, gidx(nbf), nsl, sl(nsh)
real(8) :: cen(3), rad, d
integer :: k, s, a
cen = 0.0d0
do k = 1,ng
   cen = cen + coor(:,k)
enddo
cen = cen/dble(ng)
rad = 0.0d0
do k = 1,ng
   rad = max(rad, sum((coor(:,k)-cen)**2))
enddo
rad = sqrt(rad)
nsig = 0 ; nsl = 0
do s = 1,nsh
   d = sqrt(sum((atoms(sh_atom(s))%coor*ans2bohr - cen)**2))
   if (d - rad .gt. sh_rc(s)) cycle
   nsl = nsl + 1
   sl(nsl) = s
   do a = 0,sh_dim(s)-1
      gidx(nsig+a+1) = sh_first(s) + a
   enddo
   nsig = nsig + sh_dim(s)
enddo
end subroutine xc_blk_shells

subroutine xc_blk_ao(ng, coor, nsig, nsl, sl, X0, X1, X2)
implicit none
integer,intent(in) :: ng, nsig, nsl, sl(nsl)
real(8),intent(in) :: coor(3,ng)
real(8),intent(out) :: X0(nsig,ng), X1(nsig,3,ng)
real(8),intent(out),optional :: X2(nsig,3,3,ng)
integer :: s, off
X0 = 0.0d0 ; X1 = 0.0d0
if (present(X2)) X2 = 0.0d0
off = 0
do s = 1,nsl
   if (present(X2)) then
      call GTOeval_shell_batch(sh_atom(sl(s)), sh_loc(sl(s)), coor, ng, off, nsig, X0, X1, X2)
   else
      call GTOeval_shell_batch(sh_atom(sl(s)), sh_loc(sl(s)), coor, ng, off, nsig, X0, X1)
   endif
   off = off + sh_dim(sl(s))
enddo
end subroutine xc_blk_ao

subroutine xc_blk_compress(ng, ns, X0r, X1r, X2r, has2, gidx, vthr, ns2, X0, X1, X2)
implicit none
integer,intent(in) :: ng, ns
real(8),intent(in) :: X0r(ns,ng), X1r(ns,3,ng), X2r(ns,9,*), vthr
logical,intent(in) :: has2
integer,intent(inout) :: gidx(ns)
integer,intent(out) :: ns2
real(8),intent(out) :: X0(*), X1(*), X2(*)
logical :: keep(ns)
integer :: m, k, d, m2
keep = .false.
do k = 1,ng
   do m = 1,ns
      if (abs(X0r(m,k)) .ge. vthr .or. abs(X1r(m,1,k)) .ge. vthr .or. &
          abs(X1r(m,2,k)) .ge. vthr .or. abs(X1r(m,3,k)) .ge. vthr) keep(m) = .true.
   enddo
   if (has2) then
      do d = 1,9
         do m = 1,ns
            if (abs(X2r(m,d,k)) .ge. vthr) keep(m) = .true.
         enddo
      enddo
   endif
enddo
ns2 = 0
do m = 1,ns
   if (keep(m)) then
      ns2 = ns2 + 1
      gidx(ns2) = gidx(m)
   endif
enddo
do k = 1,ng
   m2 = 0
   do m = 1,ns
      if (.not. keep(m)) cycle
      m2 = m2 + 1
      X0(m2 + (k-1)*ns2) = X0r(m,k)
      do d = 1,3
         X1(m2 + (d-1)*ns2 + (k-1)*3*ns2) = X1r(m,d,k)
      enddo
      if (has2) then
         do d = 1,9
            X2(m2 + (d-1)*ns2 + (k-1)*9*ns2) = X2r(m,d,k)
         enddo
      endif
   enddo
enddo
end subroutine xc_blk_compress

subroutine xc_hess_direct(nConts, natoms, Ptot, d2E_xc)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
external :: dgemm
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: d2E_xc(3,natoms,3,natoms)
integer,parameter :: NGC = 128
integer :: ic0, ng, k, ig, nsig, nsl, dd, x, y, m, n, ia, ib, j, la, nl
integer,allocatable :: gidx(:), sl(:), mat(:), latom(:), aloc(:)
real(8),allocatable :: cc(:,:), wg(:), X0(:), X1(:), X2(:), A0(:), A1(:), Pl(:)
real(8),allocatable :: Lm(:), Rm(:), F(:), Lf(:), Rf(:), Df(:), d2loc(:,:,:,:)
real(8),allocatable :: kv(:,:), Zs(:)
character(len=8) :: envs
integer :: hd_nchunk, nq
logical :: mg
real(8) :: hd_s1, hd_s2

call get_environment_variable("ENGINE_XC_HESS_TI", envs)
if (trim(envs) .eq. "0") then
   call xc_hess_direct_pt(nConts, natoms, Ptot, d2E_xc)
   return
endif
call get_environment_variable("ENGINE_XC_HESS_PT", envs)
if (trim(envs) .eq. "1") then
   call xc_hess_direct_pt(nConts, natoms, Ptot, d2E_xc)
   return
endif
call build_shell_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_hess nConts/=nbf',nConts,nbf; stop 1
endif
d2E_xc = 0.0d0
mg = xc_uses_tau()
nq = merge(6,5,mg)
hd_nchunk = 0 ; hd_s1 = 0.0d0 ; hd_s2 = 0.0d0
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,ig,nsig,nsl,dd,x,y,m,n,ia,ib,j,la,nl,gidx,sl,mat,latom,aloc, &
!$omp&           cc,wg,X0,X1,X2,A0,A1,Pl,Lm,Rm,F,Lf,Rf,Df,d2loc,kv,Zs)
allocate(gidx(nbf), sl(nsh), mat(nbf), latom(natoms), aloc(natoms))
allocate(cc(3,NGC), wg(NGC), kv(9,NGC))
allocate(X0(nbf*NGC), X1(3*nbf*NGC), X2(9*nbf*NGC), A0(nbf*NGC), A1(3*nbf*NGC), Pl(nbf*nbf))
allocate(Lm(3*nbf*NGC), Rm(3*nbf*NGC), F(9*nbf*nbf))
allocate(Zs(merge(9*nbf*NGC,1,mg)))
allocate(Lf(3*natoms*nq*NGC), Rf(3*natoms*nq*NGC), Df(9*natoms*natoms))
allocate(d2loc(3,natoms,3,natoms))
d2loc = 0.0d0
aloc = 0
!$omp do schedule(static,1)
do ic0 = 1, ngrids, NGC
   ng = min(NGC, ngrids-ic0+1)
   do k = 1,ng
      cc(:,k) = Grids(ic0+k-1)%coor
      wg(k) = Grids(ic0+k-1)%weight
   enddo
   if (all(wg(1:ng) .eq. 0.0d0)) cycle
   call xc_blk_shells(ng, cc, nsig, gidx, nsl, sl)
   if (nsig .eq. 0) cycle
   !$omp atomic
   hd_nchunk = hd_nchunk + 1
   !$omp atomic
   hd_s1 = hd_s1 + dble(nsig)
   !$omp atomic
   hd_s2 = hd_s2 + dble(nsig)**2
   nl = 0
   do m = 1,nsig
      ia = bf_atom(gidx(m))
      if (aloc(ia) .eq. 0) then
         nl = nl + 1
         latom(nl) = ia
         aloc(ia) = nl
      endif
      mat(m) = aloc(ia)
   enddo
   call hd_chunk(ng, nsig, nl, nsl, sl, gidx, mat, latom, cc, wg, kv, d2loc, &
                 X0, X1, X2, A0, A1, Pl, Lm, Rm, F, Lf, Rf, Df, Zs)
   do la = 1,nl
      aloc(latom(la)) = 0
   enddo
enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) d2E_xc = d2E_xc + d2loc
   enddo
end block
deallocate(gidx,sl,mat,latom,aloc,cc,wg,kv,X0,X1,X2,A0,A1,Pl,Lm,Rm,F,Lf,Rf,Df,d2loc,Zs)
!$omp end parallel
print '(A,I0,A,F8.1,A,F8.1,A,I0)', '  XC Hessian chunks ', hd_nchunk, '  mean nsig ', hd_s1/max(1,hd_nchunk), &
      '  rms nsig ', sqrt(hd_s2/max(1,hd_nchunk)), ' of ', nbf
do ia = 1,natoms
   d2E_xc(:,ia,:,ia) = 0.0d0
   do ib = 1,natoms
      if (ib .ne. ia) d2E_xc(:,ia,:,ia) = d2E_xc(:,ia,:,ia) - d2E_xc(:,ia,:,ib)
   enddo
enddo

contains

subroutine hd_chunk(ng, ns, nl, nsl, sl, gidx, mat, latom, cc, wg, kv, d2loc, &
                    X0, X1, X2, A0, A1, Pl, Lm, Rm, F, Lf, Rf, Df, Zs)
integer,intent(in) :: ng, ns, nl, nsl, sl(nsl), gidx(ns), mat(ns), latom(nl)
real(8),intent(in) :: cc(3,ng), wg(ng)
real(8),intent(inout) :: kv(9,ng), d2loc(3,natoms,3,natoms)
real(8),intent(inout) :: X0(ns,ng), X1(ns,3,ng), X2(ns,3,3,ng), A0(ns,ng), A1(ns,3,ng), Pl(ns,ns)
real(8),intent(inout) :: Lm(ns,3,ng), Rm(ns,3,ng), F(ns,3,ns,3)
real(8),intent(inout) :: Lf(3,nl,nq,ng), Rf(3,nl,nq,ng), Df(3,nl,3,nl)
real(8),intent(inout) :: Zs(*)
real(8) :: t, fr, fs, ft, rho, sig, tau, gam(3), kk(9)
integer :: m, n, x, y, dd, k, a, b
call xc_blk_ao(ng, cc, ns, nsl, sl, X0, X1, X2)
do n = 1,ns
   do m = 1,ns
      Pl(m,n) = Ptot(gidx(m),gidx(n))
   enddo
enddo
call dgemm('N','N',ns,ng,ns,1.0d0,Pl,ns,X0,ns,0.0d0,A0,ns)
call dgemm('N','N',ns,3*ng,ns,1.0d0,Pl,ns,X1,ns,0.0d0,A1,ns)
Lf = 0.0d0 ; Rf = 0.0d0
do k = 1,ng
   kv(:,k) = 0.0d0
   Lm(:,:,k) = 0.0d0 ; Rm(:,:,k) = 0.0d0
   if (wg(k) .eq. 0.0d0) cycle
   rho = sum(X0(:,k)*A0(:,k))
   do dd = 1,3
      gam(dd) = 2.0d0*sum(X1(:,dd,k)*A0(:,k))
   enddo
   sig = sum(gam**2)
   tau = 0.0d0
   if (mg) tau = 0.5d0*sum(X1(:,:,k)*A1(:,:,k))
   call xc_kern_c(mg, rho, sig, tau, kk)
   kv(:,k) = wg(k)*kk
   do x = 1,3
      do m = 1,ns
         t = gam(1)*X2(m,1,x,k) + gam(2)*X2(m,2,x,k) + gam(3)*X2(m,3,x,k)
         Lm(m,x,k) = kv(1,k)*X1(m,x,k) + 2.0d0*kv(2,k)*t
         Rm(m,x,k) = 2.0d0*kv(2,k)*t
      enddo
   enddo
   do m = 1,ns
      a = mat(m)
      do x = 1,3
         Lf(x,a,1,k) = Lf(x,a,1,k) - 2.0d0*X1(m,x,k)*A0(m,k)
         do dd = 1,3
            Lf(x,a,2+dd,k) = Lf(x,a,2+dd,k) - 2.0d0*(X2(m,dd,x,k)*A0(m,k) + X1(m,x,k)*A1(m,dd,k))
         enddo
         if (mg) Lf(x,a,6,k) = Lf(x,a,6,k) - (X2(m,1,x,k)*A1(m,1,k) + X2(m,2,x,k)*A1(m,2,k) + X2(m,3,x,k)*A1(m,3,k))
      enddo
   enddo
   do a = 1,nl
      do x = 1,3
         Lf(x,a,2,k) = 2.0d0*(gam(1)*Lf(x,a,3,k) + gam(2)*Lf(x,a,4,k) + gam(3)*Lf(x,a,5,k))
         fr = Lf(x,a,1,k) ; fs = Lf(x,a,2,k)
         Rf(x,a,1,k) = kv(3,k)*fr + kv(4,k)*fs
         Rf(x,a,2,k) = kv(4,k)*fr + kv(5,k)*fs
         Rf(x,a,3,k) = 2.0d0*kv(2,k)*Lf(x,a,3,k)
         Rf(x,a,4,k) = 2.0d0*kv(2,k)*Lf(x,a,4,k)
         Rf(x,a,5,k) = 2.0d0*kv(2,k)*Lf(x,a,5,k)
         if (mg) then
            ft = Lf(x,a,6,k)
            Rf(x,a,1,k) = Rf(x,a,1,k) + kv(7,k)*ft
            Rf(x,a,2,k) = Rf(x,a,2,k) + kv(8,k)*ft
            Rf(x,a,6,k) = kv(7,k)*fr + kv(8,k)*fs + kv(9,k)*ft
         endif
      enddo
   enddo
enddo
call dgemm('N','T',3*ns,3*ns,ng,1.0d0,Lm,3*ns,X1,3*ns,0.0d0,F,3*ns)
call dgemm('N','T',3*ns,3*ns,ng,1.0d0,X1,3*ns,Rm,3*ns,1.0d0,F,3*ns)
if (mg) then
   do k = 1,ng
      Zs(1+9*ns*(k-1):9*ns*k) = 0.5d0*kv(6,k)*reshape(X2(:,:,:,k),[9*ns])
   enddo
   call dgemm('N','T',3*ns,3*ns,3*ng,1.0d0,Zs,3*ns,X2,3*ns,1.0d0,F,3*ns)
endif
do n = 1,ns
   b = latom(mat(n))
   do y = 1,3
      do m = 1,ns
         a = latom(mat(m))
         t = 2.0d0*Pl(m,n)
         d2loc(1,a,y,b) = d2loc(1,a,y,b) + t*F(m,1,n,y)
         d2loc(2,a,y,b) = d2loc(2,a,y,b) + t*F(m,2,n,y)
         d2loc(3,a,y,b) = d2loc(3,a,y,b) + t*F(m,3,n,y)
      enddo
   enddo
enddo
call dgemm('N','T',3*nl,3*nl,nq*ng,1.0d0,Lf,3*nl,Rf,3*nl,0.0d0,Df,3*nl)
do b = 1,nl
   do a = 1,nl
      d2loc(:,latom(a),:,latom(b)) = d2loc(:,latom(a),:,latom(b)) + Df(:,a,:,b)
   enddo
enddo
end subroutine hd_chunk

end subroutine xc_hess_direct

subroutine xc_fxc_apply(nConts, nbatch, Ptot, Dtrial, Vout)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, nbatch
real(8),intent(in) :: Ptot(nConts,nConts), Dtrial(nConts,nConts,nbatch)
real(8),intent(inout) :: Vout(nConts,nConts,nbatch)
external :: dgemv
integer :: ig,m,n,dd,ib
real(8) :: coor(3),w,rho,sig,drho,dsig,Acoef
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: A0(:),Dv0(:,:),cb(:)
real(8),allocatable :: Vout_loc(:,:,:)
real(8) :: gam(3),dgam(3),Bcoef(3),bm,cm
real(8) :: e_rr,e_rs,e_ss,e_s
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vr(1),vs(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)

call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_fxc_apply nConts/=nbf'; stop 1
endif
!$omp parallel default(shared) &
!$omp&   private(ig,m,n,dd,ib,coor,w,rho,sig,drho,dsig,Acoef, &
!$omp&           v0,v1,v2,v3,A0,Dv0,Vout_loc,gam,dgam,Bcoef,e_rr,e_rs,e_ss,e_s, &
!$omp&           bm,cm,cb, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vr,vs,vra,vrb, &
!$omp&           vsaa,vsbb,vsab,v2r2,v2rs,v2s2)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(A0(nConts),Dv0(nConts,nbatch),cb(nConts))
allocate(Vout_loc(nConts,nConts,nbatch))
Vout_loc = 0.0d0
!$omp do schedule(static)
do ig = 1,ngrids
   w = Grids(ig)%weight
   if (w .eq. 0.0d0) cycle
   coor = Grids(ig)%coor
   call eval_bf_derivs(coor, v0,v1,v2,v3, 1)
   call dgemv('N',nConts,nConts,1.0d0,Ptot,nConts,v0,1,0.0d0,A0,1)
   rho = sum(v0*A0)
   do dd = 1,3
      gam(dd) = 2.0d0*sum(v1(:,dd)*A0)
   enddo
   sig = gam(1)**2+gam(2)**2+gam(3)**2
   ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
   saa(1)=sig; sbb(1)=sig; sab(1)=sig; st(1)=sig
   call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                vra,vrb,vsaa,vsbb,vsab)
   call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
   e_rr=v2r2(1); e_rs=v2rs(1); e_ss=v2s2(1); e_s=d1s(1)
   do ib = 1,nbatch
      call dgemv('N',nConts,nConts,1.0d0,Dtrial(:,:,ib),nConts,v0,1,0.0d0,Dv0(:,ib),1)
   enddo
   do ib = 1,nbatch
      drho = dot_product(v0, Dv0(:,ib))
      do dd = 1,3
         dgam(dd) = 2.0d0*dot_product(v1(:,dd), Dv0(:,ib))
      enddo
      dsig = 2.0d0*(gam(1)*dgam(1)+gam(2)*dgam(2)+gam(3)*dgam(3))
      Acoef = e_rr*drho + e_rs*dsig
      do dd = 1,3
         Bcoef(dd) = 2.0d0*(e_rs*drho + e_ss*dsig)*gam(dd) + 2.0d0*e_s*dgam(dd)
      enddo
      do m = 1,nConts
         cb(m) = Bcoef(1)*v1(m,1)+Bcoef(2)*v1(m,2)+Bcoef(3)*v1(m,3)
      enddo
      do n = 1,nConts
         bm = w*(cb(n)+Acoef*v0(n))
         cm = w*v0(n)
         do m = 1,n
            Vout_loc(m,n,ib) = Vout_loc(m,n,ib) + v0(m)*bm + cb(m)*cm
         enddo
      enddo
   enddo
enddo
!$omp end do
do ib = 1,nbatch
   do m = 1,nConts
      do n = 1,m-1
         Vout_loc(m,n,ib) = Vout_loc(n,m,ib)
      enddo
   enddo
enddo
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) Vout = Vout + Vout_loc
   enddo
end block
deallocate(v0,v1,v2,v3,A0,Dv0,cb,Vout_loc)
!$omp end parallel
end subroutine xc_fxc_apply

subroutine xc_vxc_deriv_all_pt(nConts, natoms, Ptot, dVxc_all)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: dVxc_all(nConts,nConts,3,natoms)
external :: dgemv,dgemm
integer :: ig,m,n,dd,ia,ic
real(8) :: coor(3),w,rho
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: A0(:),A1(:,:)
real(8),allocatable :: cb(:),f1(:),g2m(:),dvloc(:,:,:,:)
real(8) :: gam(3),gamA(3),rhoA(3),sigA(3),er,ers,ess,es
real(8) :: de_r(3),de_g(3,3),eg(3),bm,cm
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vr(1),vs(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)

call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_vxc_deriv_all nConts/=nbf'; stop 1
endif
dVxc_all = 0.0d0
!$omp parallel default(shared) &
!$omp&   private(ig,m,n,dd,ia,ic,coor,w,rho,v0,v1,v2,v3,A0,A1,cb,f1,g2m,dvloc, &
!$omp&           gam,gamA,rhoA,sigA,er,ers,ess,es,de_r,de_g,eg,bm,cm, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vr,vs,vra,vrb, &
!$omp&           vsaa,vsbb,vsab,v2r2,v2rs,v2s2)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(A0(nConts),A1(nConts,3))
allocate(cb(nConts),f1(nConts),g2m(nConts))
allocate(dvloc(nConts,nConts,3,natoms))
dvloc = 0.0d0
!$omp do schedule(static)
do ig = 1,ngrids
   w = Grids(ig)%weight
   if (w .eq. 0.0d0) cycle
   coor = Grids(ig)%coor
   call eval_bf_derivs(coor, v0,v1,v2,v3, 2)
   call dgemv('N',nConts,nConts,1.0d0,Ptot,nConts,v0,1,0.0d0,A0,1)
   call dgemm('N','N',nConts,3,nConts,1.0d0,Ptot,nConts,v1,nConts,0.0d0,A1,nConts)
   rho = sum(v0*A0)
   do dd = 1,3
      gam(dd) = 2.0d0*sum(v1(:,dd)*A0)
   enddo
   block
      real(8) :: g2
      g2 = gam(1)**2+gam(2)**2+gam(3)**2
      ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
      saa(1)=g2; sbb(1)=g2; sab(1)=g2; st(1)=g2
   end block
   call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                vra,vrb,vsaa,vsbb,vsab)
   call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
   er=v2r2(1); ers=v2rs(1); ess=v2s2(1); es=d1s(1)
   do dd = 1,3
      eg(dd) = 2.0d0*gam(dd)*es
   enddo
   do n = 1,nConts
      f1(n) = w*( d1r(1)*v0(n) + eg(1)*v1(n,1)+eg(2)*v1(n,2)+eg(3)*v1(n,3) )
   enddo
   do ia = 1,natoms
      do ic = 1,3
         rhoA(ic) = 0.0d0
         do m = 1,nConts
            if (bf_atom(m) .ne. ia) cycle
            rhoA(ic) = rhoA(ic) - 2.0d0*v1(m,ic)*A0(m)
         enddo
         sigA(ic) = 0.0d0
         do dd = 1,3
            gamA(dd) = 0.0d0
            do m = 1,nConts
               if (bf_atom(m) .ne. ia) cycle
               gamA(dd) = gamA(dd) - 2.0d0*( v2(m,dd,ic)*A0(m) + v1(m,ic)*A1(m,dd) )
            enddo
            sigA(ic) = sigA(ic) + 2.0d0*gam(dd)*gamA(dd)
         enddo
         de_r(ic) = er*rhoA(ic) + ers*sigA(ic)
         do dd = 1,3
            de_g(dd,ic) = 2.0d0*gamA(dd)*es + 2.0d0*gam(dd)*(ers*rhoA(ic)+ess*sigA(ic))
         enddo
      enddo
      do ic = 1,3
         do m = 1,nConts
            cb(m) = de_g(1,ic)*v1(m,1)+de_g(2,ic)*v1(m,2)+de_g(3,ic)*v1(m,3)
         enddo
         do m = 1,nConts
            if (bf_atom(m) .ne. ia) cycle
            g2m(m) = w*( eg(1)*v2(m,1,ic)+eg(2)*v2(m,2,ic)+eg(3)*v2(m,3,ic) )
         enddo
         do n = 1,nConts
            bm = w*(cb(n)+de_r(ic)*v0(n))
            cm = w*v0(n)
            do m = 1,nConts
               dvloc(m,n,ic,ia) = dvloc(m,n,ic,ia) + v0(m)*bm + cb(m)*cm
            enddo
            do m = 1,nConts
               if (bf_atom(m) .ne. ia) cycle
               dvloc(m,n,ic,ia) = dvloc(m,n,ic,ia) - v1(m,ic)*f1(n) - v0(n)*g2m(m)
               dvloc(n,m,ic,ia) = dvloc(n,m,ic,ia) - v1(m,ic)*f1(n) - v0(n)*g2m(m)
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) dVxc_all = dVxc_all + dvloc
   enddo
end block
deallocate(v0,v1,v2,v3,A0,A1,cb,f1,g2m,dvloc)
!$omp end parallel
end subroutine xc_vxc_deriv_all_pt

subroutine xc_vxc_deriv_all(nConts, natoms, Ptot, dVxc_all)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: dVxc_all(nConts,nConts,3,natoms)
external :: dgemm
integer,parameter :: NGC = 64
integer :: ic0, ng, k, ig, m, n, dd, ia, ic, y, N3
real(8) :: w, rho, gsq
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: V0c(:,:), V1c(:,:,:), V2c(:,:,:,:), A0c(:,:), A1c(:,:,:)
real(8),allocatable :: wg(:), gam(:,:), er(:), ers(:), ess(:), es(:), d1(:)
real(8),allocatable :: der(:,:), deg(:,:,:), Z(:,:), F1(:,:), G2(:,:)
real(8),allocatable :: Mloc(:,:,:), Qloc(:,:,:), Mtot(:,:,:), Qtot(:,:,:)
real(8) :: rhoA(3), gamA(3), sigA(3), eg(3)
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)
character(len=8) :: envs

call get_environment_variable("ENGINE_XC_VXCDERIV_PT", envs)
if (trim(envs) .eq. "1") then
   call xc_vxc_deriv_all_pt(nConts, natoms, Ptot, dVxc_all)
   return
endif
call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_vxc_deriv_all nConts/=nbf'; stop 1
endif
N3 = 3*natoms
allocate(Mtot(nConts,nConts,N3), Qtot(nConts,nConts,3))
Mtot = 0.0d0 ; Qtot = 0.0d0
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,ig,m,n,dd,ia,ic,y,w,rho,gsq,v0,v1,v2,v3,V0c,V1c,V2c,A0c,A1c, &
!$omp&           wg,gam,er,ers,ess,es,d1,der,deg,Z,F1,G2,Mloc,Qloc,rhoA,gamA,sigA,eg, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vra,vrb,vsaa,vsbb,vsab,v2r2,v2rs,v2s2)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(V0c(nConts,NGC), V1c(nConts,NGC,3), V2c(nConts,NGC,3,3), A0c(nConts,NGC), A1c(nConts,NGC,3))
allocate(wg(NGC), gam(3,NGC), er(NGC), ers(NGC), ess(NGC), es(NGC), d1(NGC))
allocate(der(NGC,N3), deg(3,NGC,N3), Z(nConts,NGC), F1(nConts,NGC), G2(nConts,NGC))
allocate(Mloc(nConts,nConts,N3), Qloc(nConts,nConts,3))
Mloc = 0.0d0 ; Qloc = 0.0d0
!$omp do schedule(static)
do ic0 = 1, ngrids, NGC
   ng = min(NGC, ngrids-ic0+1)
   do k = 1,ng
      ig = ic0 + k - 1
      wg(k) = Grids(ig)%weight
      if (wg(k) .eq. 0.0d0) then
         V0c(:,k) = 0.0d0 ; V1c(:,k,:) = 0.0d0 ; V2c(:,k,:,:) = 0.0d0
         cycle
      endif
      call eval_bf_derivs(Grids(ig)%coor, v0,v1,v2,v3, 2)
      V0c(:,k) = v0
      do dd = 1,3
         V1c(:,k,dd) = v1(:,dd)
         do ic = 1,3
            V2c(:,k,dd,ic) = v2(:,dd,ic)
         enddo
      enddo
   enddo
   call dgemm('N','N',nConts,ng,nConts,1.0d0,Ptot,nConts,V0c,nConts,0.0d0,A0c,nConts)
   do dd = 1,3
      call dgemm('N','N',nConts,ng,nConts,1.0d0,Ptot,nConts,V1c(1,1,dd),nConts,0.0d0,A1c(1,1,dd),nConts)
   enddo
   do k = 1,ng
      er(k)=0.0d0; ers(k)=0.0d0; ess(k)=0.0d0; es(k)=0.0d0; d1(k)=0.0d0; gam(:,k)=0.0d0
      if (wg(k) .eq. 0.0d0) cycle
      rho = sum(V0c(:,k)*A0c(:,k))
      do dd = 1,3
         gam(dd,k) = 2.0d0*sum(V1c(:,k,dd)*A0c(:,k))
      enddo
      gsq = gam(1,k)**2+gam(2,k)**2+gam(3,k)**2
      ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
      saa(1)=gsq; sbb(1)=gsq; sab(1)=gsq; st(1)=gsq
      call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                   vra,vrb,vsaa,vsbb,vsab)
      call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
      er(k)=v2r2(1); ers(k)=v2rs(1); ess(k)=v2s2(1); es(k)=d1s(1); d1(k)=d1r(1)
   enddo
   do k = 1,ng
      if (wg(k) .eq. 0.0d0) then
         der(k,:) = 0.0d0 ; deg(:,k,:) = 0.0d0
         cycle
      endif
      do ia = 1,natoms
         do ic = 1,3
            rhoA(ic) = 0.0d0
            sigA(ic) = 0.0d0
            do dd = 1,3
               gamA(dd) = 0.0d0
            enddo
            do m = bf_start(ia), bf_start(ia)+bf_len(ia)-1
               rhoA(ic) = rhoA(ic) - 2.0d0*V1c(m,k,ic)*A0c(m,k)
               do dd = 1,3
                  gamA(dd) = gamA(dd) - 2.0d0*( V2c(m,k,dd,ic)*A0c(m,k) + V1c(m,k,ic)*A1c(m,k,dd) )
               enddo
            enddo
            do dd = 1,3
               sigA(ic) = sigA(ic) + 2.0d0*gam(dd,k)*gamA(dd)
            enddo
            y = ic + 3*(ia-1)
            der(k,y) = er(k)*rhoA(ic) + ers(k)*sigA(ic)
            do dd = 1,3
               deg(dd,k,y) = 2.0d0*gamA(dd)*es(k) + 2.0d0*gam(dd,k)*(ers(k)*rhoA(ic)+ess(k)*sigA(ic))
            enddo
         enddo
      enddo
   enddo
   do y = 1,N3
      do k = 1,ng
         do m = 1,nConts
            Z(m,k) = wg(k)*( deg(1,k,y)*V1c(m,k,1) + deg(2,k,y)*V1c(m,k,2) + deg(3,k,y)*V1c(m,k,3) &
                             + 0.5d0*der(k,y)*V0c(m,k) )
         enddo
      enddo
      call dgemm('N','T',nConts,nConts,ng,1.0d0,V0c,nConts,Z,nConts,1.0d0,Mloc(1,1,y),nConts)
   enddo
   do k = 1,ng
      do dd = 1,3
         eg(dd) = 2.0d0*gam(dd,k)*es(k)
      enddo
      do m = 1,nConts
         F1(m,k) = wg(k)*( d1(k)*V0c(m,k) + eg(1)*V1c(m,k,1)+eg(2)*V1c(m,k,2)+eg(3)*V1c(m,k,3) )
      enddo
   enddo
   do ic = 1,3
      do k = 1,ng
         do dd = 1,3
            eg(dd) = 2.0d0*gam(dd,k)*es(k)
         enddo
         do m = 1,nConts
            G2(m,k) = wg(k)*( eg(1)*V2c(m,k,1,ic)+eg(2)*V2c(m,k,2,ic)+eg(3)*V2c(m,k,3,ic) )
         enddo
      enddo
      call dgemm('N','T',nConts,nConts,ng,1.0d0,V1c(1,1,ic),nConts,F1,nConts,1.0d0,Qloc(1,1,ic),nConts)
      call dgemm('N','T',nConts,nConts,ng,1.0d0,G2,nConts,V0c,nConts,1.0d0,Qloc(1,1,ic),nConts)
   enddo
enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) then
         Mtot = Mtot + Mloc
         Qtot = Qtot + Qloc
      endif
   enddo
end block
deallocate(v0,v1,v2,v3,V0c,V1c,V2c,A0c,A1c,wg,gam,er,ers,ess,es,d1,der,deg,Z,F1,G2,Mloc,Qloc)
!$omp end parallel
do ia = 1,natoms
   do ic = 1,3
      y = ic + 3*(ia-1)
      dVxc_all(:,:,ic,ia) = Mtot(:,:,y) + transpose(Mtot(:,:,y))
      do m = bf_start(ia), bf_start(ia)+bf_len(ia)-1
         dVxc_all(m,:,ic,ia) = dVxc_all(m,:,ic,ia) - Qtot(m,:,ic)
         dVxc_all(:,m,ic,ia) = dVxc_all(:,m,ic,ia) - Qtot(m,:,ic)
      enddo
   enddo
enddo
deallocate(Mtot, Qtot)
end subroutine xc_vxc_deriv_all

subroutine xc_vxc_deriv_occ(nConts, natoms, nocc, Ptot, Co, VC)
use omp_lib
implicit none
external :: dgemm
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: Ptot(nConts,nConts), Co(nConts,nocc)
real(8),intent(out) :: VC(nConts,nocc,3,natoms)
integer,parameter :: NGC = 128
integer :: ic0, ng, k, nsig, nsl, m, ia, la, nl, c, mt, nt
integer,allocatable :: gidx(:), sl(:), mat(:), latom(:), aloc(:), ar0(:), ar1(:)
real(8),allocatable :: cc(:,:), wg(:), X0(:), X1(:), X2(:), A0(:), A1(:), Pl(:), Col(:)
real(8),allocatable :: kv(:,:), gm(:,:), W0(:), FC1(:), F1(:), G2(:), Z(:), WZ(:), T1(:), OL(:)
real(8),allocatable :: QCl(:,:,:), QC(:,:,:), X0r(:), X1r(:), X2r(:), S3(:), G3(:), R3(:)
real(8) :: vthr
integer :: ns2
logical :: mg
character(len=16) :: envs
integer(omp_lock_kind),allocatable :: alk(:)
integer(16),allocatable :: VCq(:,:,:,:)

mg = xc_uses_tau()
vthr = 1.0d-7
call get_environment_variable("ENGINE_XC_PROJ_VTHR", envs)
if (len_trim(envs) .gt. 0) read(envs,*) vthr
call build_shell_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_vxc_deriv_occ nConts/=nbf'; stop 1
endif
VC = 0.0d0
allocate(VCq(nConts,nocc,3,natoms))
VCq = 0_16
allocate(QC(nConts,nocc,3), alk(natoms))
QC = 0.0d0
do ia = 1,natoms
   call omp_init_lock(alk(ia))
enddo
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,nsig,nsl,m,ia,la,nl,c,mt,nt,gidx,sl,mat,latom,aloc,ar0,ar1, &
!$omp&           cc,wg,X0,X1,X2,A0,A1,Pl,Col,kv,gm,W0,FC1,F1,G2,Z,WZ,T1,OL,QCl,X0r,X1r,X2r,ns2,S3,G3,R3)
allocate(gidx(nbf), sl(nsh), mat(nbf), latom(natoms), aloc(natoms), ar0(natoms), ar1(natoms))
allocate(cc(3,NGC), wg(NGC), kv(9,NGC), gm(3,NGC))
allocate(S3(merge(3*NGC*nocc,1,mg)), R3(merge(3*NGC*nocc,1,mg)), G3(merge(9*nbf*NGC,1,mg)))
allocate(X0(nbf*NGC), X1(3*nbf*NGC), X2(9*nbf*NGC), A0(nbf*NGC), A1(3*nbf*NGC), Pl(nbf*nbf))
allocate(Col(nbf*nocc), W0(NGC*nocc), FC1(NGC*nocc), F1(nbf*NGC), G2(3*nbf*NGC), Z(nbf*NGC))
allocate(WZ(NGC*nocc), T1(NGC*nocc), OL(nbf*nocc), QCl(nConts,nocc,3))
allocate(X0r(nbf*NGC), X1r(3*nbf*NGC), X2r(9*nbf*NGC))
QCl = 0.0d0
aloc = 0
!$omp do schedule(static,1)
do ic0 = 1, ngrids, NGC
   ng = min(NGC, ngrids-ic0+1)
   do k = 1,ng
      cc(:,k) = Grids(ic0+k-1)%coor
      wg(k) = Grids(ic0+k-1)%weight
   enddo
   if (all(wg(1:ng) .eq. 0.0d0)) cycle
   call xc_blk_shells(ng, cc, nsig, gidx, nsl, sl)
   if (nsig .eq. 0) cycle
   call xc_blk_ao(ng, cc, nsig, nsl, sl, X0r, X1r, X2r)
   call xc_blk_compress(ng, nsig, X0r, X1r, X2r, .true., gidx, vthr, ns2, X0, X1, X2)
   if (ns2 .eq. 0) cycle
   nsig = ns2
   nl = 0
   do m = 1,nsig
      ia = bf_atom(gidx(m))
      if (aloc(ia) .eq. 0) then
         nl = nl + 1
         latom(nl) = ia
         aloc(ia) = nl
         ar0(nl) = m
      endif
      mat(m) = aloc(ia)
      ar1(aloc(ia)) = m
   enddo
   call vo_chunk(ng, nsig, nl, nsl, sl, gidx, mat, latom, ar0, ar1, cc, wg, kv, gm, &
                 X0, X1, X2, A0, A1, Pl, Col, W0, FC1, F1, G2, Z, WZ, T1, OL, QCl, S3, G3, R3)
   do la = 1,nl
      aloc(latom(la)) = 0
   enddo
enddo
!$omp end do
nt = omp_get_num_threads()
do mt = 0,nt-1
   !$omp barrier
   if (omp_get_thread_num() .eq. mt) QC = QC + QCl
enddo
deallocate(gidx,sl,mat,latom,aloc,ar0,ar1,cc,wg,kv,gm,X0,X1,X2,A0,A1,Pl,Col,W0,FC1,F1,G2,Z,WZ,T1,OL,QCl,X0r,X1r,X2r,S3,G3,R3)
!$omp end parallel
VC = from_q(VCq)
deallocate(VCq)
do ia = 1,natoms
   do c = 1,3
      do m = bf_start(ia), bf_start(ia)+bf_len(ia)-1
         VC(m,:,c,ia) = VC(m,:,c,ia) - QC(m,:,c)
      enddo
   enddo
   call omp_destroy_lock(alk(ia))
enddo
deallocate(QC, alk)

contains

subroutine vo_chunk(ng, ns, nl, nsl, sl, gidx, mat, latom, ar0, ar1, cc, wg, kv, gm, &
                    X0, X1, X2, A0, A1, Pl, Col, W0, FC1, F1, G2, Z, WZ, T1, OL, QCl, S3, G3, R3)
integer,intent(in) :: ng, ns, nl, nsl, sl(nsl), gidx(ns), mat(ns), latom(nl), ar0(nl), ar1(nl)
real(8),intent(in) :: cc(3,ng), wg(ng)
real(8),intent(inout) :: kv(9,ng), gm(3,ng)
real(8),intent(inout) :: X0(ns,ng), X1(ns,3,ng), X2(ns,3,3,ng), A0(ns,ng), A1(ns,3,ng), Pl(ns,ns)
real(8),intent(inout) :: Col(ns,nocc), W0(ng,nocc), FC1(ng,nocc), F1(ns,ng), G2(ns,ng,3), Z(ns,ng)
real(8),intent(inout) :: WZ(ng,nocc), T1(ng,nocc), OL(ns,nocc), QCl(nConts,nocc,3)
real(8),intent(inout) :: S3(3,ng,*), G3(ns,3,ng,*), R3(3,ng,*)
real(8) :: rho, sig, tau, eg(3), rA, sA, tA, gA(3), der, deg(3), dvt(ng)
integer :: m, n, i, k, dd, c, la, a, r0, r1, nA
do n = 1,ns
   do m = 1,ns
      Pl(m,n) = Ptot(gidx(m),gidx(n))
   enddo
enddo
do i = 1,nocc
   do m = 1,ns
      Col(m,i) = Co(gidx(m),i)
   enddo
enddo
call dgemm('N','N',ns,ng,ns,1.0d0,Pl,ns,X0,ns,0.0d0,A0,ns)
call dgemm('N','N',ns,3*ng,ns,1.0d0,Pl,ns,X1,ns,0.0d0,A1,ns)
do k = 1,ng
   kv(:,k) = 0.0d0 ; gm(:,k) = 0.0d0
   if (wg(k) .eq. 0.0d0) cycle
   rho = sum(X0(:,k)*A0(:,k))
   do dd = 1,3
      gm(dd,k) = 2.0d0*sum(X1(:,dd,k)*A0(:,k))
   enddo
   sig = sum(gm(:,k)**2)
   tau = 0.0d0
   if (mg) tau = 0.5d0*sum(X1(:,:,k)*A1(:,:,k))
   call xc_kern_c(mg, rho, sig, tau, kv(:,k))
enddo
do k = 1,ng
   eg = 2.0d0*gm(:,k)*kv(2,k)
   do m = 1,ns
      F1(m,k) = wg(k)*( kv(1,k)*X0(m,k) + eg(1)*X1(m,1,k) + eg(2)*X1(m,2,k) + eg(3)*X1(m,3,k) )
   enddo
   do c = 1,3
      do m = 1,ns
         G2(m,k,c) = wg(k)*( eg(1)*X2(m,1,c,k) + eg(2)*X2(m,2,c,k) + eg(3)*X2(m,3,c,k) )
      enddo
   enddo
enddo
call dgemm('T','N',ng,nocc,ns,1.0d0,X0,ns,Col,ns,0.0d0,W0,ng)
call dgemm('T','N',ng,nocc,ns,1.0d0,F1,ns,Col,ns,0.0d0,FC1,ng)
do c = 1,3
   call dgemm('N','N',ns,nocc,ng,1.0d0,X1(1,c,1),3*ns,FC1,ng,0.0d0,OL,ns)
   call dgemm('N','N',ns,nocc,ng,1.0d0,G2(1,1,c),ns,W0,ng,1.0d0,OL,ns)
   do i = 1,nocc
      do m = 1,ns
         QCl(gidx(m),i,c) = QCl(gidx(m),i,c) + OL(m,i)
      enddo
   enddo
enddo
if (mg) then
   call dgemm('T','N',3*ng,nocc,ns,1.0d0,X1,ns,Col,ns,0.0d0,S3,3*ng)
   do c = 1,3
      do k = 1,ng
         do dd = 1,3
            G3(1:ns,dd,k,c) = 0.5d0*wg(k)*kv(6,k)*X2(1:ns,dd,c,k)
         enddo
      enddo
      call dgemm('N','N',ns,nocc,3*ng,1.0d0,G3(1,1,1,c),ns,S3,3*ng,0.0d0,OL,ns)
      do i = 1,nocc
         do m = 1,ns
            QCl(gidx(m),i,c) = QCl(gidx(m),i,c) + OL(m,i)
         enddo
      enddo
   enddo
endif
do la = 1,nl
   a = latom(la) ; r0 = ar0(la) ; r1 = ar1(la) ; nA = r1 - r0 + 1
   do c = 1,3
      do k = 1,ng
         dvt(k) = 0.0d0
         if (wg(k) .eq. 0.0d0) then
            Z(:,k) = 0.0d0
            cycle
         endif
         rA = 0.0d0 ; gA = 0.0d0 ; tA = 0.0d0
         do m = r0,r1
            rA = rA - 2.0d0*X1(m,c,k)*A0(m,k)
            do dd = 1,3
               gA(dd) = gA(dd) - 2.0d0*( X2(m,dd,c,k)*A0(m,k) + X1(m,c,k)*A1(m,dd,k) )
            enddo
            if (mg) tA = tA - (X2(m,1,c,k)*A1(m,1,k) + X2(m,2,c,k)*A1(m,2,k) + X2(m,3,c,k)*A1(m,3,k))
         enddo
         sA = 2.0d0*sum(gm(:,k)*gA)
         der = kv(3,k)*rA + kv(4,k)*sA
         deg = 2.0d0*gA*kv(2,k) + 2.0d0*gm(:,k)*(kv(4,k)*rA + kv(5,k)*sA)
         if (mg) then
            der = der + kv(7,k)*tA
            deg = deg + 2.0d0*gm(:,k)*kv(8,k)*tA
            dvt(k) = 0.5d0*wg(k)*(kv(7,k)*rA + kv(8,k)*sA + kv(9,k)*tA)
         endif
         do m = 1,ns
            Z(m,k) = wg(k)*( deg(1)*X1(m,1,k) + deg(2)*X1(m,2,k) + deg(3)*X1(m,3,k) &
                             + 0.5d0*der*X0(m,k) )
         enddo
      enddo
      call dgemm('T','N',ng,nocc,ns,1.0d0,Z,ns,Col,ns,0.0d0,WZ,ng)
      call dgemm('T','N',ng,nocc,nA,-1.0d0,G2(r0,1,c),ns,Col(r0,1),ns,1.0d0,WZ,ng)
      call dgemm('T','N',ng,nocc,nA,1.0d0,X1(r0,c,1),3*ns,Col(r0,1),ns,0.0d0,T1,ng)
      call dgemm('N','N',ns,nocc,ng,1.0d0,X0,ns,WZ,ng,0.0d0,OL,ns)
      call dgemm('N','N',ns,nocc,ng,1.0d0,Z,ns,W0,ng,1.0d0,OL,ns)
      call dgemm('N','N',ns,nocc,ng,-1.0d0,F1,ns,T1,ng,1.0d0,OL,ns)
      if (mg) then
         call dgemm('T','N',3*ng,nocc,nA,-1.0d0,G3(r0,1,1,c),ns,Col(r0,1),ns,0.0d0,R3,3*ng)
         do i = 1,nocc
            do k = 1,ng
               R3(1:3,k,i) = R3(1:3,k,i) + dvt(k)*S3(1:3,k,i)
            enddo
         enddo
         call dgemm('N','N',ns,nocc,3*ng,1.0d0,X1,ns,R3,3*ng,1.0d0,OL,ns)
      endif
      call omp_set_lock(alk(a))
      do i = 1,nocc
         do m = 1,ns
            VCq(gidx(m),i,c,a) = VCq(gidx(m),i,c,a) + to_q(OL(m,i))
         enddo
      enddo
      call omp_unset_lock(alk(a))
   enddo
enddo
end subroutine vo_chunk

end subroutine xc_vxc_deriv_occ

subroutine xc_vxc_deriv_atom(nConts, natoms, Ptot, iatom, ic, dVxc)
implicit none
integer,intent(in) :: nConts, natoms, iatom, ic
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: dVxc(nConts,nConts)
logical :: valid
valid = allocated(vxc_cache) .and. (vxc_cache_nConts == nConts) &
        .and. (vxc_cache_natoms == natoms)
if (valid) valid = all(vxc_cache_P == Ptot)
if (.not. valid) then
   if (allocated(vxc_cache)) deallocate(vxc_cache,vxc_cache_P)
   allocate(vxc_cache_P(nConts,nConts),vxc_cache(nConts,nConts,3,natoms))
   vxc_cache_P = Ptot
   call xc_vxc_deriv_all(nConts, natoms, Ptot, vxc_cache)
   vxc_cache_nConts = nConts
   vxc_cache_natoms = natoms
endif
dVxc = vxc_cache(:,:,ic,iatom)
end subroutine xc_vxc_deriv_atom

subroutine xc_fxc_project_vo_pt(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, nocc, nvir, nbatch
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: U(nvir,nocc,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xproj(nvir,nocc,nbatch)
external :: dgemv,dgemm
integer :: ig,a,i,dd,ib
real(8) :: coor(3),w,rho,sig,drho,dsig,Acoef
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: A0(:),Pv(:),Po(:),dPv(:,:),dPo(:,:)
real(8),allocatable :: Yl(:),Zr(:),cfa(:),Xloc(:,:,:)
real(8) :: gam(3),dgam(3),Bcoef(3),cfi
real(8) :: e_rr,e_rs,e_ss,e_s
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vr(1),vs(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)

call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_fxc_project_vo nConts/=nbf'; stop 1
endif
Xproj = 0.0d0
!$omp parallel default(shared) &
!$omp&   private(ig,a,i,dd,ib,coor,w,rho,sig,drho,dsig,Acoef, &
!$omp&           v0,v1,v2,v3,A0,Pv,Po,dPv,dPo,Yl,Zr,cfa,Xloc, &
!$omp&           gam,dgam,Bcoef,cfi,e_rr,e_rs,e_ss,e_s, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vr,vs,vra,vrb, &
!$omp&           vsaa,vsbb,vsab,v2r2,v2rs,v2s2)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(A0(nConts))
allocate(Pv(nvir),Po(nocc),dPv(nvir,3),dPo(nocc,3))
allocate(Yl(nvir),Zr(nocc),cfa(nvir))
allocate(Xloc(nvir,nocc,nbatch))
Xloc = 0.0d0
!$omp do schedule(static)
do ig = 1,ngrids
   w = Grids(ig)%weight
   if (w .eq. 0.0d0) cycle
   coor = Grids(ig)%coor
   call eval_bf_derivs(coor, v0,v1,v2,v3, 1)
   call dgemv('N',nConts,nConts,1.0d0,Ptot,nConts,v0,1,0.0d0,A0,1)
   rho = sum(v0*A0)
   do dd = 1,3
      gam(dd) = 2.0d0*sum(v1(:,dd)*A0)
   enddo
   sig = gam(1)**2+gam(2)**2+gam(3)**2
   ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
   saa(1)=sig; sbb(1)=sig; sab(1)=sig; st(1)=sig
   call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                vra,vrb,vsaa,vsbb,vsab)
   call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
   e_rr=v2r2(1); e_rs=v2rs(1); e_ss=v2s2(1); e_s=d1s(1)
   call dgemv('T',nConts,nvir,1.0d0,C_a(1,nocc+1),nConts,v0,1,0.0d0,Pv,1)
   call dgemm('T','N',nvir,3,nConts,1.0d0,C_a(1,nocc+1),nConts,v1,nConts,0.0d0,dPv,nvir)
   call dgemv('T',nConts,nocc,1.0d0,C_a,nConts,v0,1,0.0d0,Po,1)
   call dgemm('T','N',nocc,3,nConts,1.0d0,C_a,nConts,v1,nConts,0.0d0,dPo,nocc)
   do ib = 1,nbatch
      call dgemv('N',nvir,nocc,2.0d0,U(1,1,ib),nvir,Po,1,0.0d0,Yl,1)
      call dgemv('T',nvir,nocc,2.0d0,U(1,1,ib),nvir,Pv,1,0.0d0,Zr,1)
      drho = 2.0d0*dot_product(Pv,Yl)
      do dd = 1,3
         dgam(dd) = 2.0d0*( dot_product(dPv(:,dd),Yl) + dot_product(dPo(:,dd),Zr) )
      enddo
      dsig = 2.0d0*(gam(1)*dgam(1)+gam(2)*dgam(2)+gam(3)*dgam(3))
      Acoef = e_rr*drho + e_rs*dsig
      do dd = 1,3
         Bcoef(dd) = 2.0d0*(e_rs*drho + e_ss*dsig)*gam(dd) + 2.0d0*e_s*dgam(dd)
      enddo
      do a = 1,nvir
         cfa(a) = w*( Bcoef(1)*dPv(a,1)+Bcoef(2)*dPv(a,2)+Bcoef(3)*dPv(a,3) )
      enddo
      do i = 1,nocc
         cfi = w*( Acoef*Po(i) + Bcoef(1)*dPo(i,1)+Bcoef(2)*dPo(i,2)+Bcoef(3)*dPo(i,3) )
         do a = 1,nvir
            Xloc(a,i,ib) = Xloc(a,i,ib) + Pv(a)*cfi + cfa(a)*Po(i)
         enddo
      enddo
   enddo
enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) Xproj = Xproj + Xloc
   enddo
end block
deallocate(v0,v1,v2,v3,A0,Pv,Po,dPv,dPo,Yl,Zr,cfa,Xloc)
!$omp end parallel
end subroutine xc_fxc_project_vo_pt

subroutine xc_fxc_project_vo_chk(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, nocc, nvir, nbatch
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: U(nvir,nocc,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xproj(nvir,nocc,nbatch)
external :: dgemm
integer,parameter :: NGC = 128
integer :: ic0, ng, k, ig, dd, ib, a, i
real(8) :: rho, sig
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: V0c(:,:), V1c(:,:,:), A0c(:,:)
real(8),allocatable :: Pv(:,:), Po(:,:), dPv(:,:,:), dPo(:,:,:)
real(8),allocatable :: wg(:), gam(:,:), err(:), ers(:), ess(:), es(:)
real(8),allocatable :: Y(:,:), Z(:,:), cfi(:,:), cfa(:,:), Xloc(:,:,:)
real(8) :: drho, dgam(3), dsig, Acoef, Bcoef(3)
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)
character(len=8) :: envs
integer :: npg

call get_environment_variable("ENGINE_XC_PROJ_PT", envs)
if (trim(envs) .eq. "1") then
   call xc_fxc_project_vo_pt(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
   return
endif
call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_fxc_project_vo nConts/=nbf'; stop 1
endif
Xproj = 0.0d0
npg = ngrids
if (nrg .gt. 0) npg = nrg
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,ig,dd,ib,a,i,rho,sig,v0,v1,v2,v3,V0c,V1c,A0c,Pv,Po,dPv,dPo, &
!$omp&           wg,gam,err,ers,ess,es,Y,Z,cfi,cfa,Xloc,drho,dgam,dsig,Acoef,Bcoef, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vra,vrb,vsaa,vsbb,vsab,v2r2,v2rs,v2s2)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(V0c(nConts,NGC), V1c(nConts,NGC,3), A0c(nConts,NGC))
allocate(Pv(NGC,nvir), Po(NGC,nocc), dPv(NGC,nvir,3), dPo(NGC,nocc,3))
allocate(wg(NGC), gam(NGC,3), err(NGC), ers(NGC), ess(NGC), es(NGC))
allocate(Y(NGC,nvir), Z(NGC,nocc), cfi(NGC,nocc), cfa(NGC,nvir))
allocate(Xloc(nvir,nocc,nbatch))
Xloc = 0.0d0
!$omp do schedule(static)
do ic0 = 1, npg, NGC
   ng = min(NGC, npg-ic0+1)
   do k = 1,ng
      ig = ic0 + k - 1
      if (nrg .gt. 0) then
         wg(k) = rg_w(ig)
      else
         wg(k) = Grids(ig)%weight
      endif
      if (wg(k) .eq. 0.0d0) then
         V0c(:,k) = 0.0d0 ; V1c(:,k,:) = 0.0d0
         cycle
      endif
      if (nrg .gt. 0) then
         call eval_bf_derivs(rg_coor(:,ig), v0,v1,v2,v3, 1)
      else
         call eval_bf_derivs(Grids(ig)%coor, v0,v1,v2,v3, 1)
      endif
      V0c(:,k) = v0
      do dd = 1,3
         V1c(:,k,dd) = v1(:,dd)
      enddo
   enddo
   call dgemm('N','N',nConts,ng,nConts,1.0d0,Ptot,nConts,V0c,nConts,0.0d0,A0c,nConts)
   do k = 1,ng
      err(k) = 0.0d0 ; ers(k) = 0.0d0 ; ess(k) = 0.0d0 ; es(k) = 0.0d0 ; gam(k,:) = 0.0d0
      if (wg(k) .eq. 0.0d0) cycle
      rho = sum(V0c(:,k)*A0c(:,k))
      do dd = 1,3
         gam(k,dd) = 2.0d0*sum(V1c(:,k,dd)*A0c(:,k))
      enddo
      sig = gam(k,1)**2 + gam(k,2)**2 + gam(k,3)**2
      ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
      saa(1)=sig; sbb(1)=sig; sab(1)=sig; st(1)=sig
      call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                   vra,vrb,vsaa,vsbb,vsab)
      call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
      err(k)=v2r2(1); ers(k)=v2rs(1); ess(k)=v2s2(1); es(k)=d1s(1)
   enddo
   call dgemm('T','N',ng,nvir,nConts,1.0d0,V0c,nConts,C_a(1,nocc+1),nConts,0.0d0,Pv,NGC)
   call dgemm('T','N',ng,nocc,nConts,1.0d0,V0c,nConts,C_a,nConts,0.0d0,Po,NGC)
   do dd = 1,3
      call dgemm('T','N',ng,nvir,nConts,1.0d0,V1c(1,1,dd),nConts,C_a(1,nocc+1),nConts,0.0d0,dPv(1,1,dd),NGC)
      call dgemm('T','N',ng,nocc,nConts,1.0d0,V1c(1,1,dd),nConts,C_a,nConts,0.0d0,dPo(1,1,dd),NGC)
   enddo
   do ib = 1,nbatch
      call dgemm('N','T',ng,nvir,nocc,2.0d0,Po,NGC,U(1,1,ib),nvir,0.0d0,Y,NGC)
      call dgemm('N','N',ng,nocc,nvir,2.0d0,Pv,NGC,U(1,1,ib),nvir,0.0d0,Z,NGC)
      do k = 1,ng
         if (wg(k) .eq. 0.0d0) then
            cfi(k,:) = 0.0d0 ; cfa(k,:) = 0.0d0
            cycle
         endif
         drho = 0.0d0
         dgam = 0.0d0
         do a = 1,nvir
            drho = drho + Pv(k,a)*Y(k,a)
            dgam(1) = dgam(1) + dPv(k,a,1)*Y(k,a)
            dgam(2) = dgam(2) + dPv(k,a,2)*Y(k,a)
            dgam(3) = dgam(3) + dPv(k,a,3)*Y(k,a)
         enddo
         do i = 1,nocc
            dgam(1) = dgam(1) + dPo(k,i,1)*Z(k,i)
            dgam(2) = dgam(2) + dPo(k,i,2)*Z(k,i)
            dgam(3) = dgam(3) + dPo(k,i,3)*Z(k,i)
         enddo
         drho = 2.0d0*drho
         dgam = 2.0d0*dgam
         dsig = 2.0d0*(gam(k,1)*dgam(1)+gam(k,2)*dgam(2)+gam(k,3)*dgam(3))
         Acoef = err(k)*drho + ers(k)*dsig
         do dd = 1,3
            Bcoef(dd) = 2.0d0*(ers(k)*drho + ess(k)*dsig)*gam(k,dd) + 2.0d0*es(k)*dgam(dd)
         enddo
         do i = 1,nocc
            cfi(k,i) = wg(k)*( Acoef*Po(k,i) + Bcoef(1)*dPo(k,i,1)+Bcoef(2)*dPo(k,i,2)+Bcoef(3)*dPo(k,i,3) )
         enddo
         do a = 1,nvir
            cfa(k,a) = wg(k)*( Bcoef(1)*dPv(k,a,1)+Bcoef(2)*dPv(k,a,2)+Bcoef(3)*dPv(k,a,3) )
         enddo
      enddo
      call dgemm('T','N',nvir,nocc,ng,1.0d0,Pv,NGC,cfi,NGC,1.0d0,Xloc(1,1,ib),nvir)
      call dgemm('T','N',nvir,nocc,ng,1.0d0,cfa,NGC,Po,NGC,1.0d0,Xloc(1,1,ib),nvir)
   enddo
enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) Xproj = Xproj + Xloc
   enddo
end block
deallocate(v0,v1,v2,v3,V0c,V1c,A0c,Pv,Po,dPv,dPo,wg,gam,err,ers,ess,es,Y,Z,cfi,cfa,Xloc)
!$omp end parallel
end subroutine xc_fxc_project_vo_chk

subroutine xc_fxc_project_vo_grp(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
use omp_lib
implicit none
integer,intent(in) :: nConts, nocc, nvir, nbatch
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: U(nvir,nocc,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xproj(nvir,nocc,nbatch)
external :: dgemm
integer :: NGC, NBGRP
integer :: ic0, ng, k, ig, dd, a, i, b, g, gg, g0, nbk, ngrp, npg, it
real(8) :: rho, sig, dsig
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: V0c(:,:), V1c(:,:,:), A0c(:,:)
real(8),allocatable :: Pv(:,:), Po(:,:), dPv(:,:,:), dPo(:,:,:)
real(8),allocatable :: wg(:), gam(:,:), err(:), ers(:), ess(:), es(:)
real(8),allocatable :: Ut(:,:,:), Y(:,:,:), Z(:,:,:), cfi(:,:,:), cfa(:,:,:)
real(8),allocatable :: drho(:), dgam(:,:), Ac(:), Bc(:,:), Xt(:,:,:), W(:,:)
integer(omp_lock_kind),allocatable :: glk(:)
integer(16),allocatable :: Xq(:,:,:)
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)
character(len=8) :: envs

NGC = 128 ; NBGRP = 8
call get_environment_variable("ENGINE_XC_PROJ_NGC", envs)
if (len_trim(envs) .gt. 0) read(envs,*) NGC
call get_environment_variable("ENGINE_XC_PROJ_NBG", envs)
if (len_trim(envs) .gt. 0) read(envs,*) NBGRP
call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_fxc_project_vo nConts/=nbf'; stop 1
endif
Xproj = 0.0d0
allocate(Xq(nvir,nocc,nbatch))
Xq = 0_16
npg = ngrids
if (nrg .gt. 0) npg = nrg
ngrp = (nbatch + NBGRP - 1)/NBGRP
allocate(Ut(nocc,nvir,nbatch), glk(ngrp))
do b = 1,nbatch
   Ut(:,:,b) = transpose(U(:,:,b))
enddo
do g = 1,ngrp
   call omp_init_lock(glk(g))
enddo
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,ig,dd,a,i,b,g,gg,g0,nbk,it,rho,sig,dsig,v0,v1,v2,v3,V0c,V1c,A0c, &
!$omp&           Pv,Po,dPv,dPo,wg,gam,err,ers,ess,es,Y,Z,cfi,cfa,drho,dgam,Ac,Bc,Xt,W, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vra,vrb,vsaa,vsbb,vsab,v2r2,v2rs,v2s2)
it = omp_get_thread_num()
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(V0c(nConts,NGC), V1c(nConts,NGC,3), A0c(nConts,NGC))
allocate(Pv(NGC,nvir), Po(NGC,nocc), dPv(NGC,nvir,3), dPo(NGC,nocc,3))
allocate(wg(NGC), gam(NGC,3), err(NGC), ers(NGC), ess(NGC), es(NGC))
allocate(Y(NGC,nvir,NBGRP), Z(NGC,nocc,NBGRP), cfi(NGC,nocc,NBGRP), cfa(NGC,nvir,NBGRP))
allocate(drho(NGC), dgam(NGC,3), Ac(NGC), Bc(NGC,3), Xt(nvir,nocc,NBGRP), W(nvir*NBGRP,nocc))
!$omp do schedule(static)
do ic0 = 1, npg, NGC
   ng = min(NGC, npg-ic0+1)
   do k = 1,ng
      ig = ic0 + k - 1
      if (nrg .gt. 0) then
         wg(k) = rg_w(ig)
      else
         wg(k) = Grids(ig)%weight
      endif
      if (wg(k) .eq. 0.0d0) then
         V0c(:,k) = 0.0d0 ; V1c(:,k,:) = 0.0d0
         cycle
      endif
      if (nrg .gt. 0) then
         call eval_bf_derivs(rg_coor(:,ig), v0,v1,v2,v3, 1)
      else
         call eval_bf_derivs(Grids(ig)%coor, v0,v1,v2,v3, 1)
      endif
      V0c(:,k) = v0
      do dd = 1,3
         V1c(:,k,dd) = v1(:,dd)
      enddo
   enddo
   call dgemm('N','N',nConts,ng,nConts,1.0d0,Ptot,nConts,V0c,nConts,0.0d0,A0c,nConts)
   do k = 1,ng
      err(k) = 0.0d0 ; ers(k) = 0.0d0 ; ess(k) = 0.0d0 ; es(k) = 0.0d0 ; gam(k,:) = 0.0d0
      if (wg(k) .eq. 0.0d0) cycle
      rho = sum(V0c(:,k)*A0c(:,k))
      do dd = 1,3
         gam(k,dd) = 2.0d0*sum(V1c(:,k,dd)*A0c(:,k))
      enddo
      sig = gam(k,1)**2 + gam(k,2)**2 + gam(k,3)**2
      ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
      saa(1)=sig; sbb(1)=sig; sab(1)=sig; st(1)=sig
      call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                   vra,vrb,vsaa,vsbb,vsab)
      call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
      err(k)=wg(k)*v2r2(1); ers(k)=wg(k)*v2rs(1); ess(k)=wg(k)*v2s2(1); es(k)=wg(k)*d1s(1)
   enddo
   call dgemm('T','N',ng,nvir,nConts,1.0d0,V0c,nConts,C_a(1,nocc+1),nConts,0.0d0,Pv,NGC)
   call dgemm('T','N',ng,nocc,nConts,1.0d0,V0c,nConts,C_a,nConts,0.0d0,Po,NGC)
   do dd = 1,3
      call dgemm('T','N',ng,nvir,nConts,1.0d0,V1c(1,1,dd),nConts,C_a(1,nocc+1),nConts,0.0d0,dPv(1,1,dd),NGC)
      call dgemm('T','N',ng,nocc,nConts,1.0d0,V1c(1,1,dd),nConts,C_a,nConts,0.0d0,dPo(1,1,dd),NGC)
   enddo
   do gg = 0,ngrp-1
      g = mod(gg + it, ngrp) + 1
      g0 = (g-1)*NBGRP + 1
      nbk = min(NBGRP, nbatch-g0+1)
      call dgemm('N','N',ng,nvir*nbk,nocc,2.0d0,Po,NGC,Ut(1,1,g0),nocc,0.0d0,Y,NGC)
      call dgemm('N','N',ng,nocc*nbk,nvir,2.0d0,Pv,NGC,U(1,1,g0),nvir,0.0d0,Z,NGC)
      do b = 1,nbk
         drho(1:ng) = 0.0d0 ; dgam(1:ng,:) = 0.0d0
         do a = 1,nvir
            do k = 1,ng
               drho(k)   = drho(k)   + Pv(k,a)*Y(k,a,b)
               dgam(k,1) = dgam(k,1) + dPv(k,a,1)*Y(k,a,b)
               dgam(k,2) = dgam(k,2) + dPv(k,a,2)*Y(k,a,b)
               dgam(k,3) = dgam(k,3) + dPv(k,a,3)*Y(k,a,b)
            enddo
         enddo
         do i = 1,nocc
            do k = 1,ng
               dgam(k,1) = dgam(k,1) + dPo(k,i,1)*Z(k,i,b)
               dgam(k,2) = dgam(k,2) + dPo(k,i,2)*Z(k,i,b)
               dgam(k,3) = dgam(k,3) + dPo(k,i,3)*Z(k,i,b)
            enddo
         enddo
         do k = 1,ng
            drho(k) = 2.0d0*drho(k)
            dgam(k,:) = 2.0d0*dgam(k,:)
            dsig = 2.0d0*(gam(k,1)*dgam(k,1)+gam(k,2)*dgam(k,2)+gam(k,3)*dgam(k,3))
            Ac(k) = err(k)*drho(k) + ers(k)*dsig
            do dd = 1,3
               Bc(k,dd) = 2.0d0*(ers(k)*drho(k) + ess(k)*dsig)*gam(k,dd) + 2.0d0*es(k)*dgam(k,dd)
            enddo
         enddo
         do i = 1,nocc
            do k = 1,ng
               cfi(k,i,b) = Ac(k)*Po(k,i) + Bc(k,1)*dPo(k,i,1) + Bc(k,2)*dPo(k,i,2) + Bc(k,3)*dPo(k,i,3)
            enddo
         enddo
         do a = 1,nvir
            do k = 1,ng
               cfa(k,a,b) = Bc(k,1)*dPv(k,a,1) + Bc(k,2)*dPv(k,a,2) + Bc(k,3)*dPv(k,a,3)
            enddo
         enddo
      enddo
      call dgemm('T','N',nvir,nocc*nbk,ng,1.0d0,Pv,NGC,cfi,NGC,0.0d0,Xt,nvir)
      call dgemm('T','N',nvir*nbk,nocc,ng,1.0d0,cfa,NGC,Po,NGC,0.0d0,W,nvir*NBGRP)
      call omp_set_lock(glk(g))
      do b = 1,nbk
         do i = 1,nocc
            Xq(:,i,g0+b-1) = Xq(:,i,g0+b-1) + to_q(Xt(:,i,b) + W((b-1)*nvir+1:b*nvir,i))
         enddo
      enddo
      call omp_unset_lock(glk(g))
   enddo
enddo
!$omp end do
deallocate(v0,v1,v2,v3,V0c,V1c,A0c,Pv,Po,dPv,dPo,wg,gam,err,ers,ess,es,Y,Z,cfi,cfa,drho,dgam,Ac,Bc,Xt,W)
!$omp end parallel
Xproj = Xproj + from_q(Xq)
deallocate(Xq)
do g = 1,ngrp
   call omp_destroy_lock(glk(g))
enddo
deallocate(Ut, glk)
end subroutine xc_fxc_project_vo_grp

subroutine xc_fxc_project_vo(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
implicit none
integer,intent(in) :: nConts, nocc, nvir, nbatch
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: U(nvir,nocc,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xproj(nvir,nocc,nbatch)
external :: dgemm
real(8),allocatable :: L(:,:,:), Xt(:,:,:)
integer :: b
character(len=8) :: envs
call get_environment_variable("ENGINE_XC_PROJ_GRP", envs)
if (trim(envs) .eq. "1") then
   call xc_fxc_project_vo_grp(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
   return
endif
call get_environment_variable("ENGINE_XC_PROJ_CHK", envs)
if (trim(envs) .eq. "1") then
   call xc_fxc_project_vo_chk(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
   return
endif
call get_environment_variable("ENGINE_XC_PROJ_PT", envs)
if (trim(envs) .eq. "1") then
   call xc_fxc_project_vo_pt(nConts, nocc, nvir, nbatch, C_a, U, Ptot, Xproj)
   return
endif
allocate(L(nConts,nocc,nbatch), Xt(nConts,nocc,nbatch))
!$omp parallel do schedule(static)
do b = 1,nbatch
   call dgemm('N','N',nConts,nocc,nvir,2.0d0,C_a(1,nocc+1),nConts,U(1,1,b),nvir,0.0d0,L(1,1,b),nConts)
enddo
!$omp end parallel do
call xc_fxc_c(nConts, nocc, C_a, nbatch, L, Ptot, Xt)
!$omp parallel do schedule(static)
do b = 1,nbatch
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,Xt(1,1,b),nConts,0.0d0,Xproj(1,1,b),nvir)
enddo
!$omp end parallel do
deallocate(L, Xt)
end subroutine xc_fxc_project_vo

subroutine xc_fxc_c(nConts, nocc, Co, nbatch, L, Ptot, Xt)
use omp_lib
implicit none
integer,intent(in) :: nConts, nocc, nbatch
real(8),intent(in) :: Co(nConts,nocc), L(nConts,nocc,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xt(nConts,nocc,nbatch)
external :: dgemm
integer :: NGC
integer :: NBG0
character(len=8) :: nbge, ngce
integer :: ic0, ng, k, nsig, nsl, ns2, g, gg, ngrp, npg, it, nbmax, nbk, g0
integer,allocatable :: gidx(:), gidx2(:), sl(:), gst(:)
real(8),allocatable :: cc(:,:), wg(:), X0r(:), X1r(:), X0(:), X1(:), A0(:), Pl(:), Col(:), Po(:), dPo(:)
real(8),allocatable :: kv(:,:), gm(:,:), Ul(:), Z(:), Yt(:), cfi(:), H(:), XL(:), sc(:), A1(:), dZ(:)
integer(omp_lock_kind),allocatable :: glk(:)
integer(16),allocatable :: Xq(:,:,:)
real(8) :: vthr, othr
real(8),allocatable :: Cu(:,:), Lu(:,:,:), Tq(:,:), dPoc(:)
integer,allocatable :: olist(:)
integer :: no2, bb
logical :: use_loc, mg
character(len=16) :: envs

mg = xc_uses_tau()

vthr = 1.0d-7
call get_environment_variable("ENGINE_XC_PROJ_VTHR", envs)
if (len_trim(envs) .gt. 0) read(envs,*) vthr
othr = vthr
call get_environment_variable("ENGINE_XC_PROJ_OTHR", envs)
if (len_trim(envs) .gt. 0) read(envs,*) othr
call build_shell_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_fxc_c nConts/=nbf'; stop 1
endif
npg = ngrids
if (nrg .gt. 0) npg = nrg
NGC = 128
call get_environment_variable("ENGINE_XC_FXC_NGC", ngce)
if (len_trim(ngce) .gt. 0) read(ngce,*) NGC
NBG0 = 8
call get_environment_variable("ENGINE_XC_FXC_NBG", nbge)
if (len_trim(nbge) .gt. 0) read(nbge,*) NBG0
ngrp = (nbatch + NBG0 - 1)/NBG0
allocate(gst(ngrp+1), glk(ngrp))
do g = 1,ngrp+1
   gst(g) = 1 + ((g-1)*nbatch)/ngrp
enddo
nbmax = 0
do g = 1,ngrp
   nbmax = max(nbmax, gst(g+1)-gst(g))
enddo
use_loc = .false.
if (loc_ready .and. loc_nocc .eq. nocc .and. loc_n .eq. nConts) use_loc = all(Co .eq. loc_Co)
allocate(Cu(nConts,nocc), Lu(nConts,nocc,nbatch))
if (use_loc) then
   Cu = loc_Cl
   !$omp parallel do schedule(static)
   do bb = 1,nbatch
      call dgemm('N','N',nConts,nocc,nocc,1.0d0,L(1,1,bb),nConts,loc_Q,nocc,0.0d0,Lu(1,1,bb),nConts)
   enddo
   !$omp end parallel do
else
   Cu = Co
   Lu = L
endif
Xt = 0.0d0
allocate(Xq(nConts,nocc,nbatch))
Xq = 0_16
do g = 1,ngrp
   call omp_init_lock(glk(g))
enddo
!$omp parallel default(shared) &
!$omp&   private(olist,no2,dPoc,ic0,ng,k,nsig,nsl,ns2,g,gg,it,nbk,g0,gidx,gidx2,sl,cc,wg,X0r,X1r,X0,X1,A0,Pl,Col,Po,dPo, &
!$omp&           kv,gm,Ul,Z,Yt,cfi,H,XL,sc,A1,dZ)
it = omp_get_thread_num()
allocate(gidx(nbf), gidx2(nbf), sl(nsh), cc(3,NGC), wg(NGC), kv(7,NGC), gm(3,NGC))
allocate(A1(merge(3*nbf*NGC,1,mg)), dZ(merge(3*NGC*nbmax*nocc,1,mg)))
allocate(X0r(nbf*NGC), X1r(3*nbf*NGC), X0(nbf*NGC), X1(3*nbf*NGC), A0(nbf*NGC), Pl(nbf*nbf), Col(nbf*nocc))
allocate(olist(nocc), dPoc(3*NGC*nocc))
allocate(Po(NGC*nocc), dPo(3*NGC*nocc), Ul(nbf*nbmax*nocc), Z(NGC*nbmax*nocc))
allocate(Yt(nbf*nbmax*NGC), cfi(NGC*nbmax*nocc), H(nbf*nbmax*NGC), XL(nbf*nbmax*nocc), sc(9*NGC*nbmax))
!$omp do schedule(dynamic)
do ic0 = 1, npg, NGC
   ng = min(NGC, npg-ic0+1)
   do k = 1,ng
      if (nrg .gt. 0) then
         cc(:,k) = rg_coor(:,ic0+k-1) ; wg(k) = rg_w(ic0+k-1)
      else
         cc(:,k) = Grids(ic0+k-1)%coor ; wg(k) = Grids(ic0+k-1)%weight
      endif
   enddo
   if (all(wg(1:ng) .eq. 0.0d0)) cycle
   call xc_blk_shells(ng, cc, nsig, gidx, nsl, sl)
   if (nsig .eq. 0) cycle
   call xc_blk_ao(ng, cc, nsig, nsl, sl, X0r, X1r)
   call pj_compress(ng, nsig, X0r, X1r, gidx, vthr, ns2, X0, X1, gidx2)
   if (ns2 .eq. 0) cycle
   call pj_inv(ng, ns2, gidx2, wg, kv, gm, X0, X1, A0, A1, Pl, Col, Po, dPo)
   call pj_occ(ng, Po, dPo, othr, no2, olist, dPoc)
   if (no2 .eq. 0) cycle
   do gg = 0,ngrp-1
      g = mod(gg + it, ngrp) + 1
      g0 = gst(g) ; nbk = gst(g+1) - g0
      call pj_group(ng, ns2, nbk, no2, g0, gidx2, olist, X0, X1, Po, dPoc, kv, gm, Ul, Z, Yt, cfi, H, XL, sc, dZ)
      call omp_set_lock(glk(g))
      call pj_add(ns2, nbk, no2, g0, gidx2, olist, XL)
      call omp_unset_lock(glk(g))
   enddo
enddo
!$omp end do
deallocate(olist,dPoc)
deallocate(gidx,gidx2,sl,cc,wg,kv,gm,X0r,X1r,X0,X1,A0,Pl,Col,Po,dPo,Ul,Z,Yt,cfi,H,XL,sc,A1,dZ)
!$omp end parallel
Xt = from_q(Xq)
deallocate(Xq)
do g = 1,ngrp
   call omp_destroy_lock(glk(g))
enddo
deallocate(gst, glk)
if (use_loc) then
   allocate(Tq(nConts,nocc))
   do bb = 1,nbatch
      call dgemm('N','T',nConts,nocc,nocc,1.0d0,Xt(1,1,bb),nConts,loc_Q,nocc,0.0d0,Tq,nConts)
      Xt(:,:,bb) = Tq
   enddo
   deallocate(Tq)
endif
deallocate(Cu, Lu)

contains

subroutine pj_compress(ng, ns, X0r, X1r, gidx, vthr, ns2, X0, X1, gidx2)
integer,intent(in) :: ng, ns, gidx(ns)
real(8),intent(in) :: X0r(ns,ng), X1r(ns,3,ng), vthr
integer,intent(out) :: ns2, gidx2(ns)
real(8),intent(out) :: X0(*), X1(*)
logical :: keep(ns)
integer :: m, k, d, m2
keep = .false.
do k = 1,ng
   do m = 1,ns
      if (abs(X0r(m,k)) .ge. vthr .or. abs(X1r(m,1,k)) .ge. vthr .or. &
          abs(X1r(m,2,k)) .ge. vthr .or. abs(X1r(m,3,k)) .ge. vthr) keep(m) = .true.
   enddo
enddo
ns2 = 0
do m = 1,ns
   if (keep(m)) then
      ns2 = ns2 + 1
      gidx2(ns2) = gidx(m)
   endif
enddo
do k = 1,ng
   m2 = 0
   do m = 1,ns
      if (.not. keep(m)) cycle
      m2 = m2 + 1
      X0(m2 + (k-1)*ns2) = X0r(m,k)
      do d = 1,3
         X1(m2 + (d-1)*ns2 + (k-1)*3*ns2) = X1r(m,d,k)
      enddo
   enddo
enddo
end subroutine pj_compress

subroutine pj_inv(ng, ns, gidx, wg, kv, gm, X0, X1, A0, A1, Pl, Col, Po, dPo)
integer,intent(in) :: ng, ns, gidx(ns)
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng)
real(8),intent(inout) :: kv(7,ng), gm(3,ng), A0(ns,ng), A1(ns,3,ng), Pl(ns,ns), Col(ns,nocc), Po(ng,nocc), dPo(ng,nocc,3)
real(8) :: rho, sig, tau, kk(9)
integer :: m, n, i, k, dd
do n = 1,ns
   do m = 1,ns
      Pl(m,n) = Ptot(gidx(m),gidx(n))
   enddo
enddo
do i = 1,nocc
   do m = 1,ns
      Col(m,i) = Cu(gidx(m),i)
   enddo
enddo
call dgemm('N','N',ns,ng,ns,1.0d0,Pl,ns,X0,ns,0.0d0,A0,ns)
if (mg) call dgemm('N','N',ns,3*ng,ns,1.0d0,Pl,ns,X1,ns,0.0d0,A1,ns)
do k = 1,ng
   kv(:,k) = 0.0d0 ; gm(:,k) = 0.0d0
   if (wg(k) .eq. 0.0d0) cycle
   rho = sum(X0(:,k)*A0(:,k))
   do dd = 1,3
      gm(dd,k) = 2.0d0*sum(X1(:,dd,k)*A0(:,k))
   enddo
   sig = sum(gm(:,k)**2)
   tau = 0.0d0
   if (mg) tau = 0.5d0*sum(X1(:,:,k)*A1(:,:,k))
   call xc_kern_c(mg, rho, sig, tau, kk)
   kv(1,k) = wg(k)*kk(3) ; kv(2,k) = wg(k)*kk(4) ; kv(3,k) = wg(k)*kk(5) ; kv(4,k) = wg(k)*kk(2)
   kv(5,k) = wg(k)*kk(7) ; kv(6,k) = wg(k)*kk(8) ; kv(7,k) = wg(k)*kk(9)
enddo
call dgemm('T','N',ng,nocc,ns,1.0d0,X0,ns,Col,ns,0.0d0,Po,ng)
do dd = 1,3
   call dgemm('T','N',ng,nocc,ns,1.0d0,X1(1,dd,1),3*ns,Col,ns,0.0d0,dPo(1,1,dd),ng)
enddo
end subroutine pj_inv

subroutine pj_occ(ng, Po, dPo, vthr, no2, olist, dPoc)
integer,intent(in) :: ng
real(8),intent(inout) :: Po(ng,nocc)
real(8),intent(in) :: dPo(ng,nocc,3), vthr
integer,intent(out) :: no2, olist(nocc)
real(8),intent(out) :: dPoc(ng,*)
integer :: i, k, d
logical :: keep
no2 = 0
do i = 1,nocc
   keep = .false.
   do k = 1,ng
      if (abs(Po(k,i)) .ge. vthr .or. abs(dPo(k,i,1)) .ge. vthr .or. &
          abs(dPo(k,i,2)) .ge. vthr .or. abs(dPo(k,i,3)) .ge. vthr) then
         keep = .true.
         exit
      endif
   enddo
   if (.not. keep) cycle
   no2 = no2 + 1
   olist(no2) = i
   if (no2 .ne. i) Po(:,no2) = Po(:,i)
enddo
do d = 1,3
   do i = 1,no2
      dPoc(:,i + (d-1)*no2) = dPo(:,olist(i),d)
   enddo
enddo
end subroutine pj_occ

subroutine pj_group(ng, ns, nbk, no, g0, gidx, olist, X0, X1, Po, dPo, kv, gm, Ul, Z, Yt, cfi, H, XL, sc, dZ)
integer,intent(in) :: ng, ns, nbk, no, g0, gidx(ns), olist(no)
real(8),intent(in) :: X0(ns,ng), X1(ns,3,ng), Po(ng,no), dPo(ng,no,3), kv(7,ng), gm(3,ng)
real(8),intent(inout) :: Ul(ns,nbk,no), Z(ng,nbk,no), Yt(ns,nbk,ng), cfi(ng,nbk,no)
real(8),intent(inout) :: H(ns,nbk,ng), XL(ns,nbk,no), sc(ng,nbk,9), dZ(ng,nbk,no,*)
real(8) :: z1, dg(3), dsig, dtau, Ac, Bc(3), Tc
integer :: m, i, k, b, d
do i = 1,no
   do b = 1,nbk
      do m = 1,ns
         Ul(m,b,i) = Lu(gidx(m),olist(i),g0+b-1)
      enddo
   enddo
enddo
call dgemm('T','N',ng,nbk*no,ns,1.0d0,X0,ns,Ul,ns,0.0d0,Z,ng)
call dgemm('N','T',ns*nbk,ng,no,1.0d0,Ul,ns*nbk,Po,ng,0.0d0,Yt,ns*nbk)
sc(:,:,1:4) = 0.0d0
do i = 1,no
   do b = 1,nbk
      do k = 1,ng
         z1 = Z(k,b,i)
         sc(k,b,1) = sc(k,b,1) + Po(k,i)*z1
         sc(k,b,2) = sc(k,b,2) + dPo(k,i,1)*z1
         sc(k,b,3) = sc(k,b,3) + dPo(k,i,2)*z1
         sc(k,b,4) = sc(k,b,4) + dPo(k,i,3)*z1
      enddo
   enddo
enddo
if (mg) then
   do d = 1,3
      call dgemm('T','N',ng,nbk*no,ns,1.0d0,X1(1,d,1),3*ns,Ul,ns,0.0d0,dZ(1,1,1,d),ng)
   enddo
   sc(:,:,9) = 0.0d0
   do d = 1,3
      do i = 1,no
         do b = 1,nbk
            do k = 1,ng
               sc(k,b,9) = sc(k,b,9) + dPo(k,i,d)*dZ(k,b,i,d)
            enddo
         enddo
      enddo
   enddo
endif
do k = 1,ng
   do b = 1,nbk
      dg(1) = sc(k,b,2) ; dg(2) = sc(k,b,3) ; dg(3) = sc(k,b,4)
      do m = 1,ns
         dg(1) = dg(1) + X1(m,1,k)*Yt(m,b,k)
         dg(2) = dg(2) + X1(m,2,k)*Yt(m,b,k)
         dg(3) = dg(3) + X1(m,3,k)*Yt(m,b,k)
      enddo
      z1 = 2.0d0*sc(k,b,1)
      dg = 2.0d0*dg
      dsig = 2.0d0*(gm(1,k)*dg(1)+gm(2,k)*dg(2)+gm(3,k)*dg(3))
      Ac = kv(1,k)*z1 + kv(2,k)*dsig
      Bc = 2.0d0*(kv(2,k)*z1 + kv(3,k)*dsig)*gm(:,k) + 2.0d0*kv(4,k)*dg
      if (mg) then
         dtau = sc(k,b,9)
         Ac = Ac + kv(5,k)*dtau
         Bc = Bc + 2.0d0*kv(6,k)*dtau*gm(:,k)
         Tc = kv(5,k)*z1 + kv(6,k)*dsig + kv(7,k)*dtau
         sc(k,b,9) = 0.5d0*Tc
      endif
      sc(k,b,5) = Ac ; sc(k,b,6) = Bc(1) ; sc(k,b,7) = Bc(2) ; sc(k,b,8) = Bc(3)
      do m = 1,ns
         H(m,b,k) = Bc(1)*X1(m,1,k) + Bc(2)*X1(m,2,k) + Bc(3)*X1(m,3,k)
      enddo
   enddo
enddo
do i = 1,no
   do b = 1,nbk
      do k = 1,ng
         cfi(k,b,i) = sc(k,b,5)*Po(k,i) + sc(k,b,6)*dPo(k,i,1) + sc(k,b,7)*dPo(k,i,2) + sc(k,b,8)*dPo(k,i,3)
      enddo
   enddo
enddo
call dgemm('N','N',ns,nbk*no,ng,1.0d0,X0,ns,cfi,ng,0.0d0,XL,ns)
call dgemm('N','N',ns*nbk,no,ng,1.0d0,H,ns*nbk,Po,ng,1.0d0,XL,ns*nbk)
if (mg) then
   do d = 1,3
      do k = 1,ng
         do b = 1,nbk
            H(1:ns,b,k) = sc(k,b,9)*X1(1:ns,d,k)
         enddo
      enddo
      call dgemm('N','N',ns*nbk,no,ng,1.0d0,H,ns*nbk,dPo(1,1,d),ng,1.0d0,XL,ns*nbk)
   enddo
endif
end subroutine pj_group

subroutine pj_add(ns, nbk, no, g0, gidx, olist, XL)
integer,intent(in) :: ns, nbk, no, g0, gidx(ns), olist(no)
real(8),intent(in) :: XL(ns,nbk,no)
integer :: m, i, b
do b = 1,nbk
   do i = 1,no
      do m = 1,ns
         Xq(gidx(m),olist(i),g0+b-1) = Xq(gidx(m),olist(i),g0+b-1) + to_q(XL(m,b,i))
      enddo
   enddo
enddo
end subroutine pj_add

end subroutine xc_fxc_c

subroutine xc_fxc_set_loc(nConts, nocc, C)
use mod_localization, only: pipek_mezey_localize
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C(nConts,nConts)
real(8),allocatable :: S(:,:), SCl(:,:)
integer :: info, i, j, nsw
real(8) :: t0
external :: dgemm, dpotrf, dpotri
t0 = omp_get_wtime()
call xc_fxc_clear_loc()
allocate(S(nConts,nConts))
call dgemm('N','T',nConts,nConts,nConts,1.0d0,C,nConts,C,nConts,0.0d0,S,nConts)
call dpotrf('U',nConts,S,nConts,info)
if (info .eq. 0) call dpotri('U',nConts,S,nConts,info)
if (info .ne. 0) then
   deallocate(S)
   return
endif
do j = 1,nConts
   do i = j+1,nConts
      S(i,j) = S(j,i)
   enddo
enddo
allocate(loc_Co(nConts,nocc), loc_Cl(nConts,nocc), loc_Q(nocc,nocc), SCl(nConts,nocc))
loc_Co = C(:,1:nocc)
call pipek_mezey_localize(nConts, S, loc_Co, nocc, loc_Cl, nsw)
call dgemm('N','N',nConts,nocc,nConts,1.0d0,S,nConts,loc_Cl,nConts,0.0d0,SCl,nConts)
call dgemm('T','N',nocc,nocc,nConts,1.0d0,loc_Co,nConts,SCl,nConts,0.0d0,loc_Q,nocc)
deallocate(S, SCl)
loc_n = nConts ; loc_nocc = nocc ; loc_ready = .true.
print '(A,I0,A,F7.2,A)', '  CPSCF XC: occupied orbitals localized (Pipek-Mezey, ', nsw, ' sweeps, ', &
      omp_get_wtime()-t0, ' s)'
end subroutine xc_fxc_set_loc

subroutine xc_fxc_clear_loc()
implicit none
if (allocated(loc_Co)) deallocate(loc_Co, loc_Cl, loc_Q)
loc_ready = .false.
end subroutine xc_fxc_clear_loc

subroutine xc_kern_u(mg, ra, rb, ga, gb, ta, tb, fr, fs, ft, K)
implicit none
logical,intent(in) :: mg
real(8),intent(in) :: ra, rb, ga(3), gb(3), ta, tb
real(8),intent(out) :: fr(2), fs(3), ft(2), K(7,7)
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1),ta1(1),tb1(1),vta(1),vtb(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6),Km(7,7,1)
ra1(1) = ra ; rb1(1) = rb
saa(1) = sum(ga*ga) ; sab(1) = sum(ga*gb) ; sbb(1) = sum(gb*gb)
st(1) = saa(1) + 2.0d0*sab(1) + sbb(1)
ft = 0.0d0
K = 0.0d0
if (mg) then
   ta1(1) = ta ; tb1(1) = tb
   call xc_eval(2,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vra,vrb,vsaa,vsbb,vsab,ta1,tb1,vta,vtb)
   call xc_eval_fxc_mgga_pol(1,ra1,rb1,saa,sab,sbb,ta1,tb1,Km)
   K = Km(:,:,1)
   ft(1) = vta(1) ; ft(2) = vtb(1)
else
   call xc_eval(2,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vra,vrb,vsaa,vsbb,vsab)
   call xc_eval_fxc(2,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
   K(1,1) = v2r2(1) ; K(1,2) = v2r2(2) ; K(2,2) = v2r2(3)
   K(1,3) = v2rs(1) ; K(1,4) = v2rs(2) ; K(1,5) = v2rs(3)
   K(2,3) = v2rs(4) ; K(2,4) = v2rs(5) ; K(2,5) = v2rs(6)
   K(3,3) = v2s2(1) ; K(3,4) = v2s2(2) ; K(3,5) = v2s2(3)
   K(4,4) = v2s2(4) ; K(4,5) = v2s2(5) ; K(5,5) = v2s2(6)
   K(2,1) = K(1,2)
   K(3,1) = K(1,3) ; K(4,1) = K(1,4) ; K(5,1) = K(1,5)
   K(3,2) = K(2,3) ; K(4,2) = K(2,4) ; K(5,2) = K(2,5)
   K(4,3) = K(3,4) ; K(5,3) = K(3,5) ; K(5,4) = K(4,5)
endif
fr(1) = vra(1) ; fr(2) = vrb(1)
fs(1) = vsaa(1) ; fs(2) = vsab(1) ; fs(3) = vsbb(1)
end subroutine xc_kern_u

pure subroutine xc_dy_u(ga, gb, dra, drb, dga, dgb, dta, dtb, dY)
implicit none
real(8),intent(in) :: ga(3), gb(3), dra, drb, dga(3), dgb(3), dta, dtb
real(8),intent(out) :: dY(7)
dY(1) = dra ; dY(2) = drb
dY(3) = 2.0d0*sum(ga*dga)
dY(4) = sum(ga*dgb) + sum(gb*dga)
dY(5) = 2.0d0*sum(gb*dgb)
dY(6) = dta ; dY(7) = dtb
end subroutine xc_dy_u

pure subroutine xc_dv_u(fs, ga, gb, dF, dga, dgb, Aa, Ab, Ba, Bb)
implicit none
real(8),intent(in) :: fs(3), ga(3), gb(3), dF(7), dga(3), dgb(3)
real(8),intent(out) :: Aa, Ab, Ba(3), Bb(3)
Aa = dF(1) ; Ab = dF(2)
Ba = 2.0d0*dF(3)*ga + 2.0d0*fs(1)*dga + dF(4)*gb + fs(2)*dgb
Bb = 2.0d0*dF(5)*gb + 2.0d0*fs(3)*dgb + dF(4)*ga + fs(2)*dga
end subroutine xc_dv_u

subroutine xc_hess_direct_u(nConts, natoms, Pa, Pb, d2E_xc)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
external :: dgemm
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: d2E_xc(3,natoms,3,natoms)
integer,parameter :: NGC = 128
integer :: ncmp
logical :: mg
integer :: ic0, ng, k, nsig, nsl, m, ia, ib, la, nl, mt, nt
integer,allocatable :: gidx(:), sl(:), mat(:), latom(:), aloc(:)
real(8),allocatable :: cc(:,:), wg(:), X0(:), X1(:), X2(:), W(:), d2loc(:,:,:,:)

call build_shell_tables()
mg = xc_uses_tau()
ncmp = merge(13,11,mg)
d2E_xc = 0.0d0
!$omp parallel default(shared) private(ic0,ng,k,nsig,nsl,m,ia,ib,la,nl,mt,nt,gidx,sl,mat,latom,aloc, &
!$omp&   cc,wg,X0,X1,X2,W,d2loc)
allocate(gidx(nbf), sl(nsh), mat(nbf), latom(natoms), aloc(natoms), cc(3,NGC), wg(NGC))
allocate(X0(nbf*NGC), X1(3*nbf*NGC), X2(9*nbf*NGC), d2loc(3,natoms,3,natoms))
allocate(W(11*nbf*nbf + 14*nbf*NGC + 6*natoms*ncmp*NGC + 9*natoms*natoms + 9*nbf*NGC))
d2loc = 0.0d0
aloc = 0
!$omp do schedule(static,1)
do ic0 = 1, ngrids, NGC
   ng = min(NGC, ngrids-ic0+1)
   do k = 1,ng
      cc(:,k) = Grids(ic0+k-1)%coor
      wg(k) = Grids(ic0+k-1)%weight
   enddo
   if (all(wg(1:ng) .eq. 0.0d0)) cycle
   call xc_blk_shells(ng, cc, nsig, gidx, nsl, sl)
   if (nsig .eq. 0) cycle
   nl = 0
   do m = 1,nsig
      ia = bf_atom(gidx(m))
      if (aloc(ia) .eq. 0) then
         nl = nl + 1 ; latom(nl) = ia ; aloc(ia) = nl
      endif
      mat(m) = aloc(ia)
   enddo
   call xc_blk_ao(ng, cc, nsig, nsl, sl, X0, X1, X2)
   call hdu_chunk(ng, nsig, nl, gidx, mat, latom, wg, X0, X1, X2, W, d2loc)
   do la = 1,nl
      aloc(latom(la)) = 0
   enddo
enddo
!$omp end do
nt = omp_get_num_threads()
do mt = 0,nt-1
   !$omp barrier
   if (omp_get_thread_num() .eq. mt) d2E_xc = d2E_xc + d2loc
enddo
deallocate(gidx,sl,mat,latom,aloc,cc,wg,X0,X1,X2,W,d2loc)
!$omp end parallel
do ia = 1,natoms
   d2E_xc(:,ia,:,ia) = 0.0d0
   do ib = 1,natoms
      if (ib .ne. ia) d2E_xc(:,ia,:,ia) = d2E_xc(:,ia,:,ia) - d2E_xc(:,ia,:,ib)
   enddo
enddo

contains

subroutine hdu_chunk(ng, ns, nl, gidx, mat, latom, wg, X0, X1, X2, W, d2loc)
integer,intent(in) :: ng, ns, nl, gidx(ns), mat(ns), latom(nl)
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng), X2(ns,3,3,ng)
real(8),intent(inout) :: W(*), d2loc(3,natoms,3,natoms)
integer :: o
o = 1
call hdu_core(ng, ns, nl, gidx, mat, latom, wg, X0, X1, X2, &
   W(o), W(o+ns*ns), W(o+2*ns*ns), W(o+2*ns*ns+ns*ng), W(o+2*ns*ns+2*ns*ng), W(o+2*ns*ns+5*ns*ng), &
   W(o+2*ns*ns+8*ns*ng), W(o+2*ns*ns+11*ns*ng), W(o+2*ns*ns+14*ns*ng), &
   W(o+2*ns*ns+14*ns*ng+9*ns*ns), W(o+2*ns*ns+14*ns*ng+9*ns*ns+3*nl*ncmp*ng), &
   W(o+2*ns*ns+14*ns*ng+9*ns*ns+6*nl*ncmp*ng), d2loc, &
   W(o+2*ns*ns+14*ns*ng+9*ns*ns+6*nl*ncmp*ng+9*nl*nl))
end subroutine hdu_chunk

subroutine hdu_core(ng, ns, nl, gidx, mat, latom, wg, X0, X1, X2, Pla, Plb, A0a, A0b, A1a, A1b, &
                    Lm, Rm, F, Lf, Rf, Df, d2loc, Zs)
integer,intent(in) :: ng, ns, nl, gidx(ns), mat(ns), latom(nl)
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng), X2(ns,3,3,ng)
real(8),intent(inout) :: Pla(ns,ns), Plb(ns,ns), A0a(ns,ng), A0b(ns,ng), A1a(ns,3,ng), A1b(ns,3,ng)
real(8),intent(inout) :: Lm(ns,3,ng), Rm(ns,3,ng), F(ns,3,ns,3)
real(8),intent(inout) :: Lf(3,nl,ncmp,ng), Rf(3,nl,ncmp,ng), Df(3,nl,3,nl), d2loc(3,natoms,3,natoms)
real(8),intent(inout) :: Zs(*)
real(8) :: ra, rb, ga(3), gb(3), fr(2), fs(3), ft(2), K(7,7), g(3,2), t, dYv(7), dFv(7), w, ta, tb
integer :: nv
integer :: m, n, x, y, d, kk, a, b, s
do n = 1,ns
   do m = 1,ns
      Pla(m,n) = Pa(gidx(m),gidx(n)) ; Plb(m,n) = Pb(gidx(m),gidx(n))
   enddo
enddo
call dgemm('N','N',ns,ng,ns,1.0d0,Pla,ns,X0,ns,0.0d0,A0a,ns)
call dgemm('N','N',ns,ng,ns,1.0d0,Plb,ns,X0,ns,0.0d0,A0b,ns)
call dgemm('N','N',ns,3*ng,ns,1.0d0,Pla,ns,X1,ns,0.0d0,A1a,ns)
call dgemm('N','N',ns,3*ng,ns,1.0d0,Plb,ns,X1,ns,0.0d0,A1b,ns)
Lf = 0.0d0 ; Rf = 0.0d0
nv = merge(7,5,mg)
do s = 1,2
   do kk = 1,ng
      Lm(:,:,kk) = 0.0d0 ; Rm(:,:,kk) = 0.0d0
      if (mg) Zs(1+9*ns*(kk-1):9*ns*kk) = 0.0d0
      w = wg(kk)
      if (w .eq. 0.0d0) cycle
      ra = sum(X0(:,kk)*A0a(:,kk)) ; rb = sum(X0(:,kk)*A0b(:,kk))
      do d = 1,3
         ga(d) = 2.0d0*sum(X1(:,d,kk)*A0a(:,kk)) ; gb(d) = 2.0d0*sum(X1(:,d,kk)*A0b(:,kk))
      enddo
      ta = 0.0d0 ; tb = 0.0d0
      if (mg) then
         ta = 0.5d0*sum(X1(:,:,kk)*A1a(:,:,kk)) ; tb = 0.5d0*sum(X1(:,:,kk)*A1b(:,:,kk))
      endif
      call xc_kern_u(mg, ra, rb, ga, gb, ta, tb, fr, fs, ft, K)
      if (mg) Zs(1+9*ns*(kk-1):9*ns*kk) = 0.5d0*w*ft(s)*reshape(X2(:,:,:,kk),[9*ns])
      g(:,1) = 2.0d0*fs(1)*ga + fs(2)*gb
      g(:,2) = 2.0d0*fs(3)*gb + fs(2)*ga
      do x = 1,3
         do m = 1,ns
            t = g(1,s)*X2(m,1,x,kk) + g(2,s)*X2(m,2,x,kk) + g(3,s)*X2(m,3,x,kk)
            Lm(m,x,kk) = w*(fr(s)*X1(m,x,kk) + t)
            Rm(m,x,kk) = w*t
         enddo
      enddo
      if (s .eq. 2) cycle
      do m = 1,ns
         a = mat(m)
         do x = 1,3
            Lf(x,a,1,kk) = Lf(x,a,1,kk) - 2.0d0*X1(m,x,kk)*A0a(m,kk)
            Lf(x,a,2,kk) = Lf(x,a,2,kk) - 2.0d0*X1(m,x,kk)*A0b(m,kk)
            do d = 1,3
               Lf(x,a,5+d,kk) = Lf(x,a,5+d,kk) - 2.0d0*(X2(m,d,x,kk)*A0a(m,kk) + X1(m,x,kk)*A1a(m,d,kk))
               Lf(x,a,8+d,kk) = Lf(x,a,8+d,kk) - 2.0d0*(X2(m,d,x,kk)*A0b(m,kk) + X1(m,x,kk)*A1b(m,d,kk))
            enddo
            if (mg) then
               Lf(x,a,12,kk) = Lf(x,a,12,kk) - (X2(m,1,x,kk)*A1a(m,1,kk) + X2(m,2,x,kk)*A1a(m,2,kk) + X2(m,3,x,kk)*A1a(m,3,kk))
               Lf(x,a,13,kk) = Lf(x,a,13,kk) - (X2(m,1,x,kk)*A1b(m,1,kk) + X2(m,2,x,kk)*A1b(m,2,kk) + X2(m,3,x,kk)*A1b(m,3,kk))
            endif
         enddo
      enddo
      do a = 1,nl
         do x = 1,3
            if (mg) then
               call xc_dy_u(ga, gb, Lf(x,a,1,kk), Lf(x,a,2,kk), Lf(x,a,6:8,kk), Lf(x,a,9:11,kk), &
                            Lf(x,a,12,kk), Lf(x,a,13,kk), dYv)
            else
               call xc_dy_u(ga, gb, Lf(x,a,1,kk), Lf(x,a,2,kk), Lf(x,a,6:8,kk), Lf(x,a,9:11,kk), 0.0d0, 0.0d0, dYv)
            endif
            Lf(x,a,3:5,kk) = dYv(3:5)
            dFv(1:nv) = matmul(K(1:nv,1:nv), dYv(1:nv))
            Rf(x,a,1:5,kk) = w*dFv(1:5)
            if (mg) then
               Rf(x,a,12,kk) = w*dFv(6) ; Rf(x,a,13,kk) = w*dFv(7)
            endif
            Rf(x,a,6:8,kk) = w*(2.0d0*fs(1)*Lf(x,a,6:8,kk) + fs(2)*Lf(x,a,9:11,kk))
            Rf(x,a,9:11,kk) = w*(fs(2)*Lf(x,a,6:8,kk) + 2.0d0*fs(3)*Lf(x,a,9:11,kk))
         enddo
      enddo
   enddo
   call dgemm('N','T',3*ns,3*ns,ng,1.0d0,Lm,3*ns,X1,3*ns,0.0d0,F,3*ns)
   call dgemm('N','T',3*ns,3*ns,ng,1.0d0,X1,3*ns,Rm,3*ns,1.0d0,F,3*ns)
   if (mg) call dgemm('N','T',3*ns,3*ns,3*ng,1.0d0,Zs,3*ns,X2,3*ns,1.0d0,F,3*ns)
   do n = 1,ns
      b = latom(mat(n))
      do y = 1,3
         do m = 1,ns
            a = latom(mat(m))
            if (s .eq. 1) then
               t = 2.0d0*Pla(m,n)
            else
               t = 2.0d0*Plb(m,n)
            endif
            d2loc(1,a,y,b) = d2loc(1,a,y,b) + t*F(m,1,n,y)
            d2loc(2,a,y,b) = d2loc(2,a,y,b) + t*F(m,2,n,y)
            d2loc(3,a,y,b) = d2loc(3,a,y,b) + t*F(m,3,n,y)
         enddo
      enddo
   enddo
enddo
call dgemm('N','T',3*nl,3*nl,ncmp*ng,1.0d0,Lf,3*nl,Rf,3*nl,0.0d0,Df,3*nl)
do b = 1,nl
   do a = 1,nl
      d2loc(:,latom(a),:,latom(b)) = d2loc(:,latom(a),:,latom(b)) + Df(:,a,:,b)
   enddo
enddo
end subroutine hdu_core

end subroutine xc_hess_direct_u

subroutine xc_vxc_deriv_occ_u(nConts, natoms, noa, nob, Pa, Pb, Ca, Cb, VCa, VCb)
use omp_lib
implicit none
external :: dgemm
integer,intent(in) :: nConts, natoms, noa, nob
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts), Ca(nConts,noa), Cb(nConts,nob)
real(8),intent(out) :: VCa(nConts,noa,3,natoms), VCb(nConts,nob,3,natoms)
integer,parameter :: NGC = 128
integer :: ic0, ng, k, nsig, nsl, m, ia, la, nl, c, mt, nt, ns2
integer,allocatable :: gidx(:), sl(:), mat(:), latom(:), aloc(:), ar0(:), ar1(:)
real(8),allocatable :: cc(:,:), wg(:), X0r(:), X1r(:), X2r(:), X0(:), X1(:), X2(:), W(:)
real(8),allocatable :: QCa(:,:,:), QCb(:,:,:), QCla(:,:,:), QClb(:,:,:)
integer(omp_lock_kind),allocatable :: alk(:)
integer(16),allocatable :: VCqa(:,:,:,:), VCqb(:,:,:,:)
real(8) :: vthr
character(len=16) :: envs
integer :: nomx
logical :: mg

mg = xc_uses_tau()
vthr = 1.0d-7
call get_environment_variable("ENGINE_XC_PROJ_VTHR", envs)
if (len_trim(envs) .gt. 0) read(envs,*) vthr
call build_shell_tables()
nomx = max(noa, nob)
VCa = 0.0d0 ; VCb = 0.0d0
allocate(VCqa(nConts,noa,3,natoms), VCqb(nConts,nob,3,natoms))
VCqa = 0_16 ; VCqb = 0_16
allocate(QCa(nConts,noa,3), QCb(nConts,nob,3), alk(natoms))
QCa = 0.0d0 ; QCb = 0.0d0
do ia = 1,natoms
   call omp_init_lock(alk(ia))
enddo
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,nsig,nsl,m,ia,la,nl,c,mt,nt,ns2,gidx,sl,mat,latom,aloc,ar0,ar1, &
!$omp&           cc,wg,X0r,X1r,X2r,X0,X1,X2,W,QCla,QClb)
allocate(gidx(nbf), sl(nsh), mat(nbf), latom(natoms), aloc(natoms), ar0(natoms), ar1(natoms))
allocate(cc(3,NGC), wg(NGC), X0r(nbf*NGC), X1r(3*nbf*NGC), X2r(9*nbf*NGC))
allocate(X0(nbf*NGC), X1(3*nbf*NGC), X2(9*nbf*NGC))
allocate(W(2*nbf*nbf + 8*nbf*NGC + 2*nbf*nomx + 6*NGC*nomx + 9*nbf*NGC + 2*nbf*NGC + nbf*nomx &
           + 64*NGC + 6*NGC*nomx + 9*nbf*NGC + NGC))
allocate(QCla(nConts,noa,3), QClb(nConts,nob,3))
QCla = 0.0d0 ; QClb = 0.0d0
aloc = 0
!$omp do schedule(static,1)
do ic0 = 1, ngrids, NGC
   ng = min(NGC, ngrids-ic0+1)
   do k = 1,ng
      cc(:,k) = Grids(ic0+k-1)%coor
      wg(k) = Grids(ic0+k-1)%weight
   enddo
   if (all(wg(1:ng) .eq. 0.0d0)) cycle
   call xc_blk_shells(ng, cc, nsig, gidx, nsl, sl)
   if (nsig .eq. 0) cycle
   call xc_blk_ao(ng, cc, nsig, nsl, sl, X0r, X1r, X2r)
   call xc_blk_compress(ng, nsig, X0r, X1r, X2r, .true., gidx, vthr, ns2, X0, X1, X2)
   if (ns2 .eq. 0) cycle
   nsig = ns2
   nl = 0
   do m = 1,nsig
      ia = bf_atom(gidx(m))
      if (aloc(ia) .eq. 0) then
         nl = nl + 1 ; latom(nl) = ia ; aloc(ia) = nl ; ar0(nl) = m
      endif
      mat(m) = aloc(ia)
      ar1(aloc(ia)) = m
   enddo
   call vou_chunk(ng, nsig, nl, gidx, latom, ar0, ar1, wg, X0, X1, X2, W, QCla, QClb)
   do la = 1,nl
      aloc(latom(la)) = 0
   enddo
enddo
!$omp end do
nt = omp_get_num_threads()
do mt = 0,nt-1
   !$omp barrier
   if (omp_get_thread_num() .eq. mt) then
      QCa = QCa + QCla ; QCb = QCb + QClb
   endif
enddo
deallocate(gidx,sl,mat,latom,aloc,ar0,ar1,cc,wg,X0r,X1r,X2r,X0,X1,X2,W,QCla,QClb)
!$omp end parallel
VCa = from_q(VCqa) ; VCb = from_q(VCqb)
deallocate(VCqa, VCqb)
do ia = 1,natoms
   do c = 1,3
      do m = bf_start(ia), bf_start(ia)+bf_len(ia)-1
         VCa(m,:,c,ia) = VCa(m,:,c,ia) - QCa(m,:,c)
         VCb(m,:,c,ia) = VCb(m,:,c,ia) - QCb(m,:,c)
      enddo
   enddo
   call omp_destroy_lock(alk(ia))
enddo
deallocate(QCa, QCb, alk)

contains

subroutine vou_chunk(ng, ns, nl, gidx, latom, ar0, ar1, wg, X0, X1, X2, W, QCla, QClb)
integer,intent(in) :: ng, ns, nl, gidx(ns), latom(nl), ar0(nl), ar1(nl)
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng), X2(ns,3,3,ng)
real(8),intent(inout) :: W(*), QCla(nConts,noa,3), QClb(nConts,nob,3)
integer :: o(16)
o(1) = 1
o(2) = o(1) + ns*ns
o(3) = o(2) + ns*ns
o(4) = o(3) + ns*ng
o(5) = o(4) + ns*ng
o(6) = o(5) + 3*ns*ng
o(7) = o(6) + 3*ns*ng
o(8) = o(7) + ns*nomx
o(9) = o(8) + 6*ng*nomx
o(10) = o(9) + 9*ns*ng
o(11) = o(10) + 2*ns*ng
o(12) = o(11) + ns*nomx
o(13) = o(12) + 64*ng
o(14) = o(13) + 3*ng*nomx
o(15) = o(14) + 3*ng*nomx
o(16) = o(15) + 9*ns*ng
call vou_core(ng, ns, nl, gidx, latom, ar0, ar1, wg, X0, X1, X2, W(o(1)), W(o(2)), W(o(3)), W(o(4)), &
              W(o(5)), W(o(6)), W(o(7)), W(o(8)), W(o(9)), W(o(10)), W(o(11)), W(o(12)), &
              W(o(13)), W(o(14)), W(o(15)), W(o(16)), QCla, QClb)
end subroutine vou_chunk

subroutine vou_core(ng, ns, nl, gidx, latom, ar0, ar1, wg, X0, X1, X2, Pla, Plb, A0a, A0b, A1a, A1b, &
                    Col, WS, G2, FZ, OL, KV, S3, R3, G3, dvt, QCla, QClb)
integer,intent(in) :: ng, ns, nl, gidx(ns), latom(nl), ar0(nl), ar1(nl)
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng), X2(ns,3,3,ng)
real(8),intent(inout) :: Pla(ns,ns), Plb(ns,ns), A0a(ns,ng), A0b(ns,ng), A1a(ns,3,ng), A1b(ns,3,ng)
real(8),intent(inout) :: Col(ns,nomx), WS(ng,nomx,6), G2(ns,ng,3), FZ(ns,ng,2), OL(ns,nomx), KV(64,ng)
real(8),intent(inout) :: S3(3,ng,nomx), R3(3,ng,nomx), G3(ns,3,ng,3), dvt(ng)
real(8),intent(inout) :: QCla(nConts,noa,3), QClb(nConts,nob,3)
real(8) :: fr(2), fs(3), ft(2), K(7,7), ga(3), gb(3), g(3,2), rhoA(2), gAt(3,2), tAt(2), dYv(7), dF(7)
real(8) :: Aa, Ab, Bv(3,2), eg(3), ta, tb
integer :: m, n, i, kk, d, c, la, a, r0, r1, nA, s, no, nv
nv = merge(7,5,mg)
do n = 1,ns
   do m = 1,ns
      Pla(m,n) = Pa(gidx(m),gidx(n)) ; Plb(m,n) = Pb(gidx(m),gidx(n))
   enddo
enddo
call dgemm('N','N',ns,ng,ns,1.0d0,Pla,ns,X0,ns,0.0d0,A0a,ns)
call dgemm('N','N',ns,ng,ns,1.0d0,Plb,ns,X0,ns,0.0d0,A0b,ns)
call dgemm('N','N',ns,3*ng,ns,1.0d0,Pla,ns,X1,ns,0.0d0,A1a,ns)
call dgemm('N','N',ns,3*ng,ns,1.0d0,Plb,ns,X1,ns,0.0d0,A1b,ns)
do kk = 1,ng
   KV(:,kk) = 0.0d0
   if (wg(kk) .eq. 0.0d0) cycle
   do d = 1,3
      ga(d) = 2.0d0*sum(X1(:,d,kk)*A0a(:,kk)) ; gb(d) = 2.0d0*sum(X1(:,d,kk)*A0b(:,kk))
   enddo
   ta = 0.0d0 ; tb = 0.0d0
   if (mg) then
      ta = 0.5d0*sum(X1(:,:,kk)*A1a(:,:,kk)) ; tb = 0.5d0*sum(X1(:,:,kk)*A1b(:,:,kk))
   endif
   call xc_kern_u(mg, sum(X0(:,kk)*A0a(:,kk)), sum(X0(:,kk)*A0b(:,kk)), ga, gb, ta, tb, fr, fs, ft, K)
   KV(1:2,kk) = fr ; KV(3:5,kk) = fs ; KV(6:8,kk) = ga ; KV(9:11,kk) = gb ; KV(12:13,kk) = ft
   KV(14:62,kk) = reshape(K,(/49/))
enddo
do s = 1,2
   if (s .eq. 1) then
      no = noa
      do i = 1,no
         do m = 1,ns
            Col(m,i) = Ca(gidx(m),i)
         enddo
      enddo
   else
      no = nob
      do i = 1,no
         do m = 1,ns
            Col(m,i) = Cb(gidx(m),i)
         enddo
      enddo
   endif
   do kk = 1,ng
      ga = KV(6:8,kk) ; gb = KV(9:11,kk) ; fs = KV(3:5,kk)
      if (s .eq. 1) then
         eg = 2.0d0*fs(1)*ga + fs(2)*gb
      else
         eg = 2.0d0*fs(3)*gb + fs(2)*ga
      endif
      do m = 1,ns
         FZ(m,kk,1) = wg(kk)*( KV(s,kk)*X0(m,kk) + eg(1)*X1(m,1,kk) + eg(2)*X1(m,2,kk) + eg(3)*X1(m,3,kk) )
      enddo
      do c = 1,3
         do m = 1,ns
            G2(m,kk,c) = wg(kk)*( eg(1)*X2(m,1,c,kk) + eg(2)*X2(m,2,c,kk) + eg(3)*X2(m,3,c,kk) )
         enddo
      enddo
   enddo
   call dgemm('T','N',ng,no,ns,1.0d0,X0,ns,Col,ns,0.0d0,WS(1,1,1),ng)
   call dgemm('T','N',ng,no,ns,1.0d0,FZ(1,1,1),ns,Col,ns,0.0d0,WS(1,1,2),ng)
   do c = 1,3
      call dgemm('N','N',ns,no,ng,1.0d0,X1(1,c,1),3*ns,WS(1,1,2),ng,0.0d0,OL,ns)
      call dgemm('N','N',ns,no,ng,1.0d0,G2(1,1,c),ns,WS(1,1,1),ng,1.0d0,OL,ns)
      do i = 1,no
         do m = 1,ns
            if (s .eq. 1) then
               QCla(gidx(m),i,c) = QCla(gidx(m),i,c) + OL(m,i)
            else
               QClb(gidx(m),i,c) = QClb(gidx(m),i,c) + OL(m,i)
            endif
         enddo
      enddo
   enddo
   if (mg) then
      call dgemm('T','N',3*ng,no,ns,1.0d0,X1,ns,Col,ns,0.0d0,S3,3*ng)
      do c = 1,3
         do kk = 1,ng
            do d = 1,3
               G3(1:ns,d,kk,c) = 0.5d0*wg(kk)*KV(11+s,kk)*X2(1:ns,d,c,kk)
            enddo
         enddo
         call dgemm('N','N',ns,no,3*ng,1.0d0,G3(1,1,1,c),ns,S3,3*ng,0.0d0,OL,ns)
         do i = 1,no
            do m = 1,ns
               if (s .eq. 1) then
                  QCla(gidx(m),i,c) = QCla(gidx(m),i,c) + OL(m,i)
               else
                  QClb(gidx(m),i,c) = QClb(gidx(m),i,c) + OL(m,i)
               endif
            enddo
         enddo
      enddo
   endif
   do la = 1,nl
      a = latom(la) ; r0 = ar0(la) ; r1 = ar1(la) ; nA = r1 - r0 + 1
      do c = 1,3
         do kk = 1,ng
            dvt(kk) = 0.0d0
            if (wg(kk) .eq. 0.0d0) then
               FZ(:,kk,2) = 0.0d0
               cycle
            endif
            rhoA = 0.0d0 ; gAt = 0.0d0 ; tAt = 0.0d0
            do m = r0,r1
               rhoA(1) = rhoA(1) - 2.0d0*X1(m,c,kk)*A0a(m,kk)
               rhoA(2) = rhoA(2) - 2.0d0*X1(m,c,kk)*A0b(m,kk)
               do d = 1,3
                  gAt(d,1) = gAt(d,1) - 2.0d0*( X2(m,d,c,kk)*A0a(m,kk) + X1(m,c,kk)*A1a(m,d,kk) )
                  gAt(d,2) = gAt(d,2) - 2.0d0*( X2(m,d,c,kk)*A0b(m,kk) + X1(m,c,kk)*A1b(m,d,kk) )
               enddo
               if (mg) then
                  tAt(1) = tAt(1) - (X2(m,1,c,kk)*A1a(m,1,kk) + X2(m,2,c,kk)*A1a(m,2,kk) + X2(m,3,c,kk)*A1a(m,3,kk))
                  tAt(2) = tAt(2) - (X2(m,1,c,kk)*A1b(m,1,kk) + X2(m,2,c,kk)*A1b(m,2,kk) + X2(m,3,c,kk)*A1b(m,3,kk))
               endif
            enddo
            ga = KV(6:8,kk) ; gb = KV(9:11,kk) ; fs = KV(3:5,kk)
            K = reshape(KV(14:62,kk),(/7,7/))
            call xc_dy_u(ga, gb, rhoA(1), rhoA(2), gAt(:,1), gAt(:,2), tAt(1), tAt(2), dYv)
            dF(1:nv) = matmul(K(1:nv,1:nv), dYv(1:nv))
            if (mg) dvt(kk) = 0.5d0*wg(kk)*dF(5+s)
            call xc_dv_u(fs, ga, gb, dF, gAt(:,1), gAt(:,2), Aa, Ab, Bv(:,1), Bv(:,2))
            if (s .eq. 1) then
               dF(1) = Aa
            else
               dF(1) = Ab
            endif
            do m = 1,ns
               FZ(m,kk,2) = wg(kk)*( Bv(1,s)*X1(m,1,kk) + Bv(2,s)*X1(m,2,kk) + Bv(3,s)*X1(m,3,kk) &
                                     + 0.5d0*dF(1)*X0(m,kk) )
            enddo
         enddo
         call dgemm('T','N',ng,no,ns,1.0d0,FZ(1,1,2),ns,Col,ns,0.0d0,WS(1,1,3),ng)
         call dgemm('T','N',ng,no,nA,-1.0d0,G2(r0,1,c),ns,Col(r0,1),ns,1.0d0,WS(1,1,3),ng)
         call dgemm('T','N',ng,no,nA,1.0d0,X1(r0,c,1),3*ns,Col(r0,1),ns,0.0d0,WS(1,1,4),ng)
         call dgemm('N','N',ns,no,ng,1.0d0,X0,ns,WS(1,1,3),ng,0.0d0,OL,ns)
         call dgemm('N','N',ns,no,ng,1.0d0,FZ(1,1,2),ns,WS(1,1,1),ng,1.0d0,OL,ns)
         call dgemm('N','N',ns,no,ng,-1.0d0,FZ(1,1,1),ns,WS(1,1,4),ng,1.0d0,OL,ns)
         if (mg) then
            call dgemm('T','N',3*ng,no,nA,-1.0d0,G3(r0,1,1,c),ns,Col(r0,1),ns,0.0d0,R3,3*ng)
            do i = 1,no
               do kk = 1,ng
                  R3(1:3,kk,i) = R3(1:3,kk,i) + dvt(kk)*S3(1:3,kk,i)
               enddo
            enddo
            call dgemm('N','N',ns,no,3*ng,1.0d0,X1,ns,R3,3*ng,1.0d0,OL,ns)
         endif
         call omp_set_lock(alk(a))
         do i = 1,no
            do m = 1,ns
               if (s .eq. 1) then
                  VCqa(gidx(m),i,c,a) = VCqa(gidx(m),i,c,a) + to_q(OL(m,i))
               else
                  VCqb(gidx(m),i,c,a) = VCqb(gidx(m),i,c,a) + to_q(OL(m,i))
               endif
            enddo
         enddo
         call omp_unset_lock(alk(a))
      enddo
   enddo
enddo
end subroutine vou_core

end subroutine xc_vxc_deriv_occ_u

subroutine xc_fxc_apply_u(nConts, nb, Pa, Pb, Da, Db, Va, Vb)
use omp_lib
implicit none
external :: dgemm
integer,intent(in) :: nConts, nb
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: Da(nConts,nConts,nb), Db(nConts,nConts,nb)
real(8),intent(inout) :: Va(nConts,nConts,nb), Vb(nConts,nConts,nb)
integer,parameter :: NGC = 128, NBG = 8
integer :: ic0, ng, k, nsig, nsl, ns2, g, gg, ngrp, it, npg
integer,allocatable :: gidx(:), sl(:)
real(8),allocatable :: cc(:,:), wg(:), X0r(:), X1r(:), X0(:), X1(:), W(:)
integer(omp_lock_kind),allocatable :: glk(:)
integer(16),allocatable :: Vqa(:,:,:), Vqb(:,:,:)
real(8) :: vthr
character(len=16) :: envs
logical :: mg

mg = xc_uses_tau()
vthr = 1.0d-7
call get_environment_variable("ENGINE_XC_PROJ_VTHR", envs)
if (len_trim(envs) .gt. 0) read(envs,*) vthr
call build_shell_tables()
npg = ngrids
if (nrg .gt. 0) npg = nrg
ngrp = (nb + NBG - 1)/NBG
allocate(glk(ngrp))
do g = 1,ngrp
   call omp_init_lock(glk(g))
enddo
allocate(Vqa(nConts,nConts,nb), Vqb(nConts,nConts,nb))
Vqa = 0_16 ; Vqb = 0_16
!$omp parallel default(shared) private(ic0,ng,k,nsig,nsl,ns2,g,gg,it,gidx,sl,cc,wg,X0r,X1r,X0,X1,W)
it = omp_get_thread_num()
allocate(gidx(nbf), sl(nsh), cc(3,NGC), wg(NGC), X0r(nbf*NGC), X1r(3*nbf*NGC), X0(nbf*NGC), X1(3*nbf*NGC))
allocate(W(6*nbf*nbf + 12*nbf*NGC + 64*NGC + merge(18*nbf*NGC,0,mg)))
!$omp do schedule(dynamic)
do ic0 = 1, npg, NGC
   ng = min(NGC, npg-ic0+1)
   do k = 1,ng
      if (nrg .gt. 0) then
         cc(:,k) = rg_coor(:,ic0+k-1) ; wg(k) = rg_w(ic0+k-1)
      else
         cc(:,k) = Grids(ic0+k-1)%coor ; wg(k) = Grids(ic0+k-1)%weight
      endif
   enddo
   if (all(wg(1:ng) .eq. 0.0d0)) cycle
   call xc_blk_shells(ng, cc, nsig, gidx, nsl, sl)
   if (nsig .eq. 0) cycle
   call xc_blk_ao(ng, cc, nsig, nsl, sl, X0r, X1r)
   call xc_blk_compress(ng, nsig, X0r, X1r, X0r, .false., gidx, vthr, ns2, X0, X1, X0r)
   if (ns2 .eq. 0) cycle
   call fau_chunk(ng, ns2, gidx, wg, X0, X1, W, it, ngrp)
enddo
!$omp end do
deallocate(gidx,sl,cc,wg,X0r,X1r,X0,X1,W)
!$omp end parallel
Va = Va + from_q(Vqa) ; Vb = Vb + from_q(Vqb)
deallocate(Vqa, Vqb)
do g = 1,ngrp
   call omp_destroy_lock(glk(g))
enddo
deallocate(glk)

contains

subroutine fau_chunk(ng, ns, gidx, wg, X0, X1, W, it, ngrp)
integer,intent(in) :: ng, ns, gidx(ns), it, ngrp
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng)
real(8),intent(inout) :: W(*)
integer :: o(10)
o(1) = 1
o(2) = o(1) + 2*ns*ns
o(3) = o(2) + 2*ns*ns
o(4) = o(3) + 4*ns*ng
o(5) = o(4) + 6*ns*ng
o(6) = o(5) + 2*ns*ng
o(7) = o(6) + 2*ns*ns
o(8) = o(7) + 64*ng
o(9) = o(8) + 6*ns*ng
o(10) = o(9) + 6*ns*ng
call fau_core(ng, ns, gidx, wg, X0, X1, W(o(1)), W(o(2)), W(o(3)), W(o(4)), W(o(5)), W(o(6)), W(o(7)), &
              W(o(8)), W(o(9)), W(o(10)), it, ngrp)
end subroutine fau_chunk

subroutine fau_core(ng, ns, gidx, wg, X0, X1, Pl, Dl, A0, A1, H, Vl, KV, A1r, DX, Ht, it, ngrp)
integer,intent(in) :: ng, ns, gidx(ns), it, ngrp
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng)
real(8),intent(inout) :: Pl(ns,ns,2), Dl(ns,ns,2), A0(ns,ng,4), A1(ns,3,ng,2), H(ns,ng,2), Vl(ns,ns,2), KV(64,ng)
real(8),intent(inout) :: A1r(ns,3,ng,2), DX(ns,3,ng,2), Ht(ns,3,ng,2)
real(8) :: fr(2), fs(3), ft(2), K(7,7), ga(3), gb(3), dra, drb, dga(3), dgb(3), dta, dtb, Y(7), dF(7), Aa, Ab, Ba(3), Bb(3), ta, tb
integer :: m, n, kk, d, b, g, gg, g0, g1, s, nv
nv = merge(7,5,mg)
do n = 1,ns
   do m = 1,ns
      Pl(m,n,1) = Pa(gidx(m),gidx(n)) ; Pl(m,n,2) = Pb(gidx(m),gidx(n))
   enddo
enddo
call dgemm('N','N',ns,ng,ns,1.0d0,Pl(1,1,1),ns,X0,ns,0.0d0,A0(1,1,1),ns)
call dgemm('N','N',ns,ng,ns,1.0d0,Pl(1,1,2),ns,X0,ns,0.0d0,A0(1,1,2),ns)
if (mg) then
   call dgemm('N','N',ns,3*ng,ns,1.0d0,Pl(1,1,1),ns,X1,ns,0.0d0,A1r(1,1,1,1),ns)
   call dgemm('N','N',ns,3*ng,ns,1.0d0,Pl(1,1,2),ns,X1,ns,0.0d0,A1r(1,1,1,2),ns)
endif
do kk = 1,ng
   KV(:,kk) = 0.0d0
   if (wg(kk) .eq. 0.0d0) cycle
   do d = 1,3
      ga(d) = 2.0d0*sum(X1(:,d,kk)*A0(:,kk,1)) ; gb(d) = 2.0d0*sum(X1(:,d,kk)*A0(:,kk,2))
   enddo
   ta = 0.0d0 ; tb = 0.0d0
   if (mg) then
      ta = 0.5d0*sum(X1(:,:,kk)*A1r(:,:,kk,1)) ; tb = 0.5d0*sum(X1(:,:,kk)*A1r(:,:,kk,2))
   endif
   call xc_kern_u(mg, sum(X0(:,kk)*A0(:,kk,1)), sum(X0(:,kk)*A0(:,kk,2)), ga, gb, ta, tb, fr, fs, ft, K)
   KV(1:2,kk) = fr ; KV(3:5,kk) = fs ; KV(6:8,kk) = ga ; KV(9:11,kk) = gb ; KV(12:13,kk) = ft
   KV(14:62,kk) = reshape(K,(/49/))
enddo
do gg = 0,ngrp-1
   g = mod(gg + it, ngrp) + 1
   g0 = (g-1)*NBG + 1 ; g1 = min(nb, g*NBG)
   do b = g0,g1
      do n = 1,ns
         do m = 1,ns
            Dl(m,n,1) = Da(gidx(m),gidx(n),b) ; Dl(m,n,2) = Db(gidx(m),gidx(n),b)
         enddo
      enddo
      call dgemm('N','N',ns,ng,ns,1.0d0,Dl(1,1,1),ns,X0,ns,0.0d0,A0(1,1,3),ns)
      call dgemm('N','N',ns,ng,ns,1.0d0,Dl(1,1,2),ns,X0,ns,0.0d0,A0(1,1,4),ns)
      if (mg) then
         call dgemm('N','N',ns,3*ng,ns,1.0d0,Dl(1,1,1),ns,X1,ns,0.0d0,DX(1,1,1,1),ns)
         call dgemm('N','N',ns,3*ng,ns,1.0d0,Dl(1,1,2),ns,X1,ns,0.0d0,DX(1,1,1,2),ns)
      endif
      do kk = 1,ng
         if (wg(kk) .eq. 0.0d0) then
            H(:,kk,:) = 0.0d0
            if (mg) Ht(:,:,kk,:) = 0.0d0
            cycle
         endif
         dra = sum(X0(:,kk)*A0(:,kk,3)) ; drb = sum(X0(:,kk)*A0(:,kk,4))
         do d = 1,3
            dga(d) = 2.0d0*sum(X1(:,d,kk)*A0(:,kk,3)) ; dgb(d) = 2.0d0*sum(X1(:,d,kk)*A0(:,kk,4))
         enddo
         fs = KV(3:5,kk) ; ga = KV(6:8,kk) ; gb = KV(9:11,kk)
         K = reshape(KV(14:62,kk),(/7,7/))
         dta = 0.0d0 ; dtb = 0.0d0
         if (mg) then
            dta = 0.5d0*sum(X1(:,:,kk)*DX(:,:,kk,1)) ; dtb = 0.5d0*sum(X1(:,:,kk)*DX(:,:,kk,2))
         endif
         call xc_dy_u(ga, gb, dra, drb, dga, dgb, dta, dtb, Y)
         dF(1:nv) = matmul(K(1:nv,1:nv), Y(1:nv))
         call xc_dv_u(fs, ga, gb, dF, dga, dgb, Aa, Ab, Ba, Bb)
         if (mg) then
            Ht(:,:,kk,1) = 0.25d0*wg(kk)*dF(6)*X1(:,:,kk)
            Ht(:,:,kk,2) = 0.25d0*wg(kk)*dF(7)*X1(:,:,kk)
         endif
         do m = 1,ns
            H(m,kk,1) = wg(kk)*( 0.5d0*Aa*X0(m,kk) + Ba(1)*X1(m,1,kk) + Ba(2)*X1(m,2,kk) + Ba(3)*X1(m,3,kk) )
            H(m,kk,2) = wg(kk)*( 0.5d0*Ab*X0(m,kk) + Bb(1)*X1(m,1,kk) + Bb(2)*X1(m,2,kk) + Bb(3)*X1(m,3,kk) )
         enddo
      enddo
      do s = 1,2
         call dgemm('N','T',ns,ns,ng,1.0d0,X0,ns,H(1,1,s),ns,0.0d0,Vl(1,1,s),ns)
         if (mg) call dgemm('N','T',ns,ns,3*ng,1.0d0,Ht(1,1,1,s),ns,X1,ns,1.0d0,Vl(1,1,s),ns)
      enddo
      call omp_set_lock(glk(g))
      do n = 1,ns
         do m = 1,ns
            Vqa(gidx(m),gidx(n),b) = Vqa(gidx(m),gidx(n),b) + to_q(Vl(m,n,1) + Vl(n,m,1))
            Vqb(gidx(m),gidx(n),b) = Vqb(gidx(m),gidx(n),b) + to_q(Vl(m,n,2) + Vl(n,m,2))
         enddo
      enddo
      call omp_unset_lock(glk(g))
   enddo
enddo
end subroutine fau_core

end subroutine xc_fxc_apply_u

subroutine xc_fxc_c_u(nConts, noa, nob, Ca, Cb, nbatch, La, Lb, Pa, Pb, Xa, Xb)
use omp_lib
implicit none
integer,intent(in) :: nConts, noa, nob, nbatch
real(8),intent(in) :: Ca(nConts,noa), Cb(nConts,nob), La(nConts,noa,nbatch), Lb(nConts,nob,nbatch)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: Xa(nConts,noa,nbatch), Xb(nConts,nob,nbatch)
external :: dgemm
integer,parameter :: NGC = 128, NBG0 = 8
integer :: ic0, ng, k, nsig, nsl, ns2, g, gg, ngrp, npg, it, nbmax, nbk, g0, nom
integer,allocatable :: gidx(:), gidx2(:), sl(:), gst(:)
real(8),allocatable :: cc(:,:), wg(:), X0r(:), X1r(:), X0(:), X1(:), W(:)
integer(omp_lock_kind),allocatable :: glk(:)
integer(16),allocatable :: Xqa(:,:,:), Xqb(:,:,:)
real(8) :: vthr
character(len=16) :: envs
logical :: mg

mg = xc_uses_tau()
vthr = 1.0d-7
call get_environment_variable("ENGINE_XC_PROJ_VTHR", envs)
if (len_trim(envs) .gt. 0) read(envs,*) vthr
call build_shell_tables()
npg = ngrids
if (nrg .gt. 0) npg = nrg
nom = max(noa, nob)
ngrp = (nbatch + NBG0 - 1)/NBG0
allocate(gst(ngrp+1), glk(ngrp))
do g = 1,ngrp+1
   gst(g) = 1 + ((g-1)*nbatch)/ngrp
enddo
nbmax = 0
do g = 1,ngrp
   nbmax = max(nbmax, gst(g+1)-gst(g))
enddo
Xa = 0.0d0 ; Xb = 0.0d0
allocate(Xqa(nConts,noa,nbatch), Xqb(nConts,nob,nbatch))
Xqa = 0_16 ; Xqb = 0_16
do g = 1,ngrp
   call omp_init_lock(glk(g))
enddo
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,nsig,nsl,ns2,g,gg,it,nbk,g0,gidx,gidx2,sl,cc,wg,X0r,X1r,X0,X1,W)
it = omp_get_thread_num()
allocate(gidx(nbf), gidx2(nbf), sl(nsh), cc(3,NGC), wg(NGC))
allocate(X0r(nbf*NGC), X1r(3*nbf*NGC), X0(nbf*NGC), X1(3*nbf*NGC))
allocate(W(2*nbf*nbf + 2*nbf*NGC + 2*nbf*nom + 8*NGC*nom + 64*NGC + merge(6*nbf*NGC,0,mg) &
           + 2*(nbf*nbmax*nom + NGC*nbmax*nom + nbf*nbmax*NGC + NGC*nbmax*nom + nbf*nbmax*NGC &
                + nbf*nbmax*nom + 9*NGC*nbmax + merge(3*NGC*nbmax*nom,0,mg))))
!$omp do schedule(dynamic)
do ic0 = 1, npg, NGC
   ng = min(NGC, npg-ic0+1)
   do k = 1,ng
      if (nrg .gt. 0) then
         cc(:,k) = rg_coor(:,ic0+k-1) ; wg(k) = rg_w(ic0+k-1)
      else
         cc(:,k) = Grids(ic0+k-1)%coor ; wg(k) = Grids(ic0+k-1)%weight
      endif
   enddo
   if (all(wg(1:ng) .eq. 0.0d0)) cycle
   call xc_blk_shells(ng, cc, nsig, gidx, nsl, sl)
   if (nsig .eq. 0) cycle
   call xc_blk_ao(ng, cc, nsig, nsl, sl, X0r, X1r)
   gidx2(1:nsig) = gidx(1:nsig)
   call xc_blk_compress(ng, nsig, X0r, X1r, X0r, .false., gidx2, vthr, ns2, X0, X1, X0r)
   if (ns2 .eq. 0) cycle
   call pu_chunk(ng, ns2, gidx2, wg, X0, X1, W, it)
enddo
!$omp end do
deallocate(gidx,gidx2,sl,cc,wg,X0r,X1r,X0,X1,W)
!$omp end parallel
Xa = from_q(Xqa) ; Xb = from_q(Xqb)
deallocate(Xqa, Xqb)
do g = 1,ngrp
   call omp_destroy_lock(glk(g))
enddo
deallocate(gst, glk)

contains

subroutine pu_chunk(ng, ns, gidx, wg, X0, X1, W, it)
integer,intent(in) :: ng, ns, gidx(ns), it
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng)
real(8),intent(inout) :: W(*)
integer :: o(20), gsz(2), s
o(1) = 1
o(2) = o(1) + 2*ns*ns
o(3) = o(2) + 2*ns*ng
o(4) = o(3) + 2*ns*nom
o(5) = o(4) + 8*ng*nom
o(6) = o(5) + 64*ng
o(7) = o(6) + merge(6*ns*ng,0,mg)
gsz(1) = ns*nbmax*nom + NGC*nbmax*nom + ns*nbmax*NGC + NGC*nbmax*nom + ns*nbmax*NGC + ns*nbmax*nom + 9*NGC*nbmax &
         + merge(3*NGC*nbmax*nom,0,mg)
call pu_core(ng, ns, gidx, wg, X0, X1, W(o(1)), W(o(2)), W(o(3)), W(o(4)), W(o(4)+2*ng*nom), W(o(5)), &
             W(o(6)), W(o(7)), W(o(7)+gsz(1)), it)
end subroutine pu_chunk

subroutine pu_core(ng, ns, gidx, wg, X0, X1, Pl, A0, Col, Po, dPo, KV, A1r, WA, WB, it)
integer,intent(in) :: ng, ns, gidx(ns), it
real(8),intent(in) :: wg(ng), X0(ns,ng), X1(ns,3,ng)
real(8),intent(inout) :: Pl(ns,ns,2), A0(ns,ng,2), Col(ns,nom,2), Po(ng,nom,2), dPo(ng,nom,3,2), KV(64,ng)
real(8),intent(inout) :: A1r(ns,3,ng,2), WA(*), WB(*)
real(8) :: ga(3), gb(3), fr(2), fs(3), ft(2), K(7,7), t0, t1, ta, tb
integer :: m, n, i, kk, d, s, no, g, gg, g0, nbk
t0 = omp_get_wtime()
do n = 1,ns
   do m = 1,ns
      Pl(m,n,1) = Pa(gidx(m),gidx(n)) ; Pl(m,n,2) = Pb(gidx(m),gidx(n))
   enddo
enddo
do s = 1,2
   call dgemm('N','N',ns,ng,ns,1.0d0,Pl(1,1,s),ns,X0,ns,0.0d0,A0(1,1,s),ns)
   if (mg) call dgemm('N','N',ns,3*ng,ns,1.0d0,Pl(1,1,s),ns,X1,ns,0.0d0,A1r(1,1,1,s),ns)
   no = merge(noa, nob, s .eq. 1)
   do i = 1,no
      do m = 1,ns
         if (s .eq. 1) then
            Col(m,i,s) = Ca(gidx(m),i)
         else
            Col(m,i,s) = Cb(gidx(m),i)
         endif
      enddo
   enddo
   call dgemm('T','N',ng,no,ns,1.0d0,X0,ns,Col(1,1,s),ns,0.0d0,Po(1,1,s),ng)
   do d = 1,3
      call dgemm('T','N',ng,no,ns,1.0d0,X1(1,d,1),3*ns,Col(1,1,s),ns,0.0d0,dPo(1,1,d,s),ng)
   enddo
enddo
t1 = omp_get_wtime()
!$omp atomic
fxu_t(1) = fxu_t(1) + (t1-t0)
t0 = t1
do kk = 1,ng
   KV(:,kk) = 0.0d0
   if (wg(kk) .eq. 0.0d0) cycle
   do d = 1,3
      ga(d) = 2.0d0*sum(X1(:,d,kk)*A0(:,kk,1)) ; gb(d) = 2.0d0*sum(X1(:,d,kk)*A0(:,kk,2))
   enddo
   ta = 0.0d0 ; tb = 0.0d0
   if (mg) then
      ta = 0.5d0*sum(X1(:,:,kk)*A1r(:,:,kk,1)) ; tb = 0.5d0*sum(X1(:,:,kk)*A1r(:,:,kk,2))
   endif
   call xc_kern_u(mg, sum(X0(:,kk)*A0(:,kk,1)), sum(X0(:,kk)*A0(:,kk,2)), ga, gb, ta, tb, fr, fs, ft, K)
   KV(1:3,kk) = fs ; KV(4:6,kk) = ga ; KV(7:9,kk) = gb
   KV(10:58,kk) = reshape(K,(/49/))
   KV(59,kk) = wg(kk)
enddo
t1 = omp_get_wtime()
!$omp atomic
fxu_t(2) = fxu_t(2) + (t1-t0)
t0 = t1
do gg = 0,ngrp-1
   g = mod(gg + it, ngrp) + 1
   g0 = gst(g) ; nbk = gst(g+1) - g0
   call pu_group(ng, ns, nbk, g0, gidx, X0, X1, Po, dPo, KV, WA, WB)
   call omp_set_lock(glk(g))
   call pu_add(ns, nbk, g0, gidx, WA, WB)
   call omp_unset_lock(glk(g))
enddo
t1 = omp_get_wtime()
!$omp atomic
fxu_t(3) = fxu_t(3) + (t1-t0)
end subroutine pu_core

subroutine pu_group(ng, ns, nbk, g0, gidx, X0, X1, Po, dPo, KV, WA, WB)
integer,intent(in) :: ng, ns, nbk, g0, gidx(ns)
real(8),intent(in) :: X0(ns,ng), X1(ns,3,ng), Po(ng,nom,2), dPo(ng,nom,3,2), KV(64,ng)
real(8),intent(inout) :: WA(*), WB(*)
integer :: o(8)
o(1) = 1
o(2) = o(1) + ns*nbmax*nom
o(3) = o(2) + NGC*nbmax*nom
o(4) = o(3) + ns*nbmax*NGC
o(5) = o(4) + NGC*nbmax*nom
o(6) = o(5) + ns*nbmax*NGC
o(7) = o(6) + ns*nbmax*nom
o(8) = o(7) + 9*NGC*nbmax
call pu_spin_in(ng, ns, nbk, noa, g0, gidx, X0, X1, Po(1,1,1), dPo(1,1,1,1), La, &
                WA(o(1)), WA(o(2)), WA(o(3)), WA(o(7)), WA(o(8)))
call pu_spin_in(ng, ns, nbk, nob, g0, gidx, X0, X1, Po(1,1,2), dPo(1,1,1,2), Lb, &
                WB(o(1)), WB(o(2)), WB(o(3)), WB(o(7)), WB(o(8)))
call pu_kernel(ng, ns, nbk, X1, KV, WA(o(3)), WB(o(3)), WA(o(7)), WB(o(7)), WA(o(5)), WB(o(5)))
call pu_spin_out(ng, ns, nbk, noa, X0, X1, Po(1,1,1), dPo(1,1,1,1), WA(o(4)), WA(o(5)), WA(o(6)), WA(o(7)))
call pu_spin_out(ng, ns, nbk, nob, X0, X1, Po(1,1,2), dPo(1,1,1,2), WB(o(4)), WB(o(5)), WB(o(6)), WB(o(7)))
end subroutine pu_group

subroutine pu_spin_in(ng, ns, nbk, no, g0, gidx, X0, X1, Po, dPo, Ls, Ul, Z, Yt, sc, dZ)
integer,intent(in) :: ng, ns, nbk, no, g0, gidx(ns)
real(8),intent(in) :: X0(ns,ng), X1(ns,3,ng), Po(ng,nom), dPo(ng,nom,3), Ls(nConts,no,nbatch)
real(8),intent(inout) :: Ul(ns,nbk,no), Z(ng,nbk,no), Yt(ns,nbk,ng), sc(ng,nbk,9), dZ(ng,nbk,no,*)
integer :: m, i, b, k, d
real(8) :: z1
do i = 1,no
   do b = 1,nbk
      do m = 1,ns
         Ul(m,b,i) = Ls(gidx(m),i,g0+b-1)
      enddo
   enddo
enddo
call dgemm('T','N',ng,nbk*no,ns,1.0d0,X0,ns,Ul,ns,0.0d0,Z,ng)
call dgemm('N','T',ns*nbk,ng,no,1.0d0,Ul,ns*nbk,Po,ng,0.0d0,Yt,ns*nbk)
sc(:,:,1:4) = 0.0d0
do i = 1,no
   do b = 1,nbk
      do k = 1,ng
         z1 = Z(k,b,i)
         sc(k,b,1) = sc(k,b,1) + Po(k,i)*z1
         sc(k,b,2) = sc(k,b,2) + dPo(k,i,1)*z1
         sc(k,b,3) = sc(k,b,3) + dPo(k,i,2)*z1
         sc(k,b,4) = sc(k,b,4) + dPo(k,i,3)*z1
      enddo
   enddo
enddo
if (mg) then
   do d = 1,3
      call dgemm('T','N',ng,nbk*no,ns,1.0d0,X1(1,d,1),3*ns,Ul,ns,0.0d0,dZ(1,1,1,d),ng)
   enddo
   sc(:,:,9) = 0.0d0
   do d = 1,3
      do i = 1,no
         do b = 1,nbk
            do k = 1,ng
               sc(k,b,9) = sc(k,b,9) + dPo(k,i,d)*dZ(k,b,i,d)
            enddo
         enddo
      enddo
   enddo
endif
end subroutine pu_spin_in

subroutine pu_kernel(ng, ns, nbk, X1, KV, Yta, Ytb, sca, scb, Ha, Hb)
integer,intent(in) :: ng, ns, nbk
real(8),intent(in) :: X1(ns,3,ng), KV(64,ng), Yta(ns,nbk,ng), Ytb(ns,nbk,ng)
real(8),intent(inout) :: sca(ng,nbk,9), scb(ng,nbk,9), Ha(ns,nbk,ng), Hb(ns,nbk,ng)
real(8) :: dga(3), dgb(3), dY(7), dF(7), Aa, Ab, Ba(3), Bb(3), Km(7,7), w, dta, dtb
integer :: m, k, b, d, nv
nv = merge(7,5,mg)
do k = 1,ng
   w = KV(59,k)
   Km = reshape(KV(10:58,k),(/7,7/))
   do b = 1,nbk
      if (w .eq. 0.0d0) then
         sca(k,b,5:9) = 0.0d0 ; scb(k,b,5:9) = 0.0d0
         Ha(:,b,k) = 0.0d0 ; Hb(:,b,k) = 0.0d0
         cycle
      endif
      do d = 1,3
         dga(d) = 2.0d0*(sca(k,b,1+d) + sum(X1(:,d,k)*Yta(:,b,k)))
         dgb(d) = 2.0d0*(scb(k,b,1+d) + sum(X1(:,d,k)*Ytb(:,b,k)))
      enddo
      dta = 0.0d0 ; dtb = 0.0d0
      if (mg) then
         dta = sca(k,b,9) ; dtb = scb(k,b,9)
      endif
      call xc_dy_u(KV(4:6,k), KV(7:9,k), 2.0d0*sca(k,b,1), 2.0d0*scb(k,b,1), dga, dgb, dta, dtb, dY)
      dF(1:nv) = matmul(Km(1:nv,1:nv), dY(1:nv))
      call xc_dv_u(KV(1:3,k), KV(4:6,k), KV(7:9,k), dF, dga, dgb, Aa, Ab, Ba, Bb)
      sca(k,b,5) = w*Aa ; sca(k,b,6:8) = w*Ba
      scb(k,b,5) = w*Ab ; scb(k,b,6:8) = w*Bb
      if (mg) then
         sca(k,b,9) = 0.5d0*w*dF(6) ; scb(k,b,9) = 0.5d0*w*dF(7)
      endif
      do m = 1,ns
         Ha(m,b,k) = w*(Ba(1)*X1(m,1,k) + Ba(2)*X1(m,2,k) + Ba(3)*X1(m,3,k))
         Hb(m,b,k) = w*(Bb(1)*X1(m,1,k) + Bb(2)*X1(m,2,k) + Bb(3)*X1(m,3,k))
      enddo
   enddo
enddo
end subroutine pu_kernel

subroutine pu_spin_out(ng, ns, nbk, no, X0, X1, Po, dPo, cfi, H, XL, sc)
integer,intent(in) :: ng, ns, nbk, no
real(8),intent(in) :: X0(ns,ng), X1(ns,3,ng), Po(ng,nom), dPo(ng,nom,3), sc(ng,nbk,9)
real(8),intent(inout) :: H(ns,nbk,ng), cfi(ng,nbk,no), XL(ns,nbk,no)
integer :: i, b, k, d
do i = 1,no
   do b = 1,nbk
      do k = 1,ng
         cfi(k,b,i) = sc(k,b,5)*Po(k,i) + sc(k,b,6)*dPo(k,i,1) + sc(k,b,7)*dPo(k,i,2) + sc(k,b,8)*dPo(k,i,3)
      enddo
   enddo
enddo
call dgemm('N','N',ns,nbk*no,ng,1.0d0,X0,ns,cfi,ng,0.0d0,XL,ns)
call dgemm('N','N',ns*nbk,no,ng,1.0d0,H,ns*nbk,Po,ng,1.0d0,XL,ns*nbk)
if (mg) then
   do d = 1,3
      do k = 1,ng
         do b = 1,nbk
            H(1:ns,b,k) = sc(k,b,9)*X1(1:ns,d,k)
         enddo
      enddo
      call dgemm('N','N',ns*nbk,no,ng,1.0d0,H,ns*nbk,dPo(1,1,d),ng,1.0d0,XL,ns*nbk)
   enddo
endif
end subroutine pu_spin_out

subroutine pu_add(ns, nbk, g0, gidx, WA, WB)
integer,intent(in) :: ns, nbk, g0, gidx(ns)
real(8),intent(in) :: WA(*), WB(*)
call pu_add1(ns, nbk, noa, g0, gidx, WA(xl_off(ns)), Xqa)
call pu_add1(ns, nbk, nob, g0, gidx, WB(xl_off(ns)), Xqb)
end subroutine pu_add

integer function xl_off(ns)
integer,intent(in) :: ns
xl_off = 1 + ns*nbmax*nom + NGC*nbmax*nom + ns*nbmax*NGC + NGC*nbmax*nom + ns*nbmax*NGC
end function xl_off

subroutine pu_add1(ns, nbk, no, g0, gidx, XL, Xs)
integer,intent(in) :: ns, nbk, no, g0, gidx(ns)
real(8),intent(in) :: XL(ns,nbk,no)
integer(16),intent(inout) :: Xs(nConts,no,nbatch)
integer :: m, i, b
do b = 1,nbk
   do i = 1,no
      do m = 1,ns
         Xs(gidx(m),i,g0+b-1) = Xs(gidx(m),i,g0+b-1) + to_q(XL(m,b,i))
      enddo
   enddo
enddo
end subroutine pu_add1

end subroutine xc_fxc_c_u

subroutine xc_fxc_project_pt(nConts, nDL, nDR, CDL, CDR, nPL, nPR, CPL, CPR, nbatch, Xd, Ptot, Xproj)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, nDL, nDR, nPL, nPR, nbatch
real(8),intent(in) :: CDL(nConts,nDL), CDR(nConts,nDR)
real(8),intent(in) :: CPL(nConts,nPL), CPR(nConts,nPR)
real(8),intent(in) :: Xd(nDL,nDR,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xproj(nPL,nPR,nbatch)
external :: dgemv,dgemm
integer :: ig,a,i,dd,ib
real(8) :: coor(3),w,rho,sig,drho,dsig,Acoef
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: A0(:)
real(8),allocatable :: PLd(:),PRd(:),dPLd(:,:),dPRd(:,:)
real(8),allocatable :: PLp(:),PRp(:),dPLp(:,:),dPRp(:,:)
real(8),allocatable :: Yl(:),Zr(:),cfa(:),Xloc(:,:,:)
real(8) :: gam(3),dgam(3),Bcoef(3),cfi
real(8) :: e_rr,e_rs,e_ss,e_s
real(8) :: ra1(1),rb1(1),saa(1),sbb(1),sab(1),st(1)
real(8) :: ex1(1),ec1(1),d1r(1),d1s(1),vr(1),vs(1),vra(1),vrb(1),vsaa(1),vsbb(1),vsab(1)
real(8) :: v2r2(3),v2rs(6),v2s2(6)

call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_fxc_project nConts/=nbf'; stop 1
endif
Xproj = 0.0d0
!$omp parallel default(shared) &
!$omp&   private(ig,a,i,dd,ib,coor,w,rho,sig,drho,dsig,Acoef, &
!$omp&           v0,v1,v2,v3,A0,PLd,PRd,dPLd,dPRd,PLp,PRp,dPLp,dPRp,Yl,Zr,cfa,Xloc, &
!$omp&           gam,dgam,Bcoef,cfi,e_rr,e_rs,e_ss,e_s, &
!$omp&           ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s,vr,vs,vra,vrb, &
!$omp&           vsaa,vsbb,vsab,v2r2,v2rs,v2s2)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(A0(nConts))
allocate(PLd(nDL),PRd(nDR),dPLd(nDL,3),dPRd(nDR,3))
allocate(PLp(nPL),PRp(nPR),dPLp(nPL,3),dPRp(nPR,3))
allocate(Yl(nDL),Zr(nDR),cfa(nPL))
allocate(Xloc(nPL,nPR,nbatch))
Xloc = 0.0d0
!$omp do schedule(static)
do ig = 1,ngrids
   w = Grids(ig)%weight
   if (w .eq. 0.0d0) cycle
   coor = Grids(ig)%coor
   call eval_bf_derivs(coor, v0,v1,v2,v3, 1)
   call dgemv('N',nConts,nConts,1.0d0,Ptot,nConts,v0,1,0.0d0,A0,1)
   rho = sum(v0*A0)
   do dd = 1,3
      gam(dd) = 2.0d0*sum(v1(:,dd)*A0)
   enddo
   sig = gam(1)**2+gam(2)**2+gam(3)**2
   ra1(1)=0.5d0*rho; rb1(1)=0.5d0*rho
   saa(1)=sig; sbb(1)=sig; sab(1)=sig; st(1)=sig
   call xc_eval(1,0,1,ra1,rb1,saa,sbb,sab,st,ex1,ec1,d1r,d1s, &
                vra,vrb,vsaa,vsbb,vsab)
   call xc_eval_fxc(1,1,ra1,rb1,saa,sbb,sab,st,v2r2,v2rs,v2s2)
   e_rr=v2r2(1); e_rs=v2rs(1); e_ss=v2s2(1); e_s=d1s(1)
   call dgemv('T',nConts,nDL,1.0d0,CDL,nConts,v0,1,0.0d0,PLd,1)
   call dgemm('T','N',nDL,3,nConts,1.0d0,CDL,nConts,v1,nConts,0.0d0,dPLd,nDL)
   call dgemv('T',nConts,nDR,1.0d0,CDR,nConts,v0,1,0.0d0,PRd,1)
   call dgemm('T','N',nDR,3,nConts,1.0d0,CDR,nConts,v1,nConts,0.0d0,dPRd,nDR)
   call dgemv('T',nConts,nPL,1.0d0,CPL,nConts,v0,1,0.0d0,PLp,1)
   call dgemm('T','N',nPL,3,nConts,1.0d0,CPL,nConts,v1,nConts,0.0d0,dPLp,nPL)
   call dgemv('T',nConts,nPR,1.0d0,CPR,nConts,v0,1,0.0d0,PRp,1)
   call dgemm('T','N',nPR,3,nConts,1.0d0,CPR,nConts,v1,nConts,0.0d0,dPRp,nPR)
   do ib = 1,nbatch
      call dgemv('N',nDL,nDR,1.0d0,Xd(1,1,ib),nDL,PRd,1,0.0d0,Yl,1)
      call dgemv('T',nDL,nDR,1.0d0,Xd(1,1,ib),nDL,PLd,1,0.0d0,Zr,1)
      drho = 2.0d0*dot_product(PLd,Yl)
      do dd = 1,3
         dgam(dd) = 2.0d0*( dot_product(dPLd(:,dd),Yl) + dot_product(dPRd(:,dd),Zr) )
      enddo
      dsig = 2.0d0*(gam(1)*dgam(1)+gam(2)*dgam(2)+gam(3)*dgam(3))
      Acoef = e_rr*drho + e_rs*dsig
      do dd = 1,3
         Bcoef(dd) = 2.0d0*(e_rs*drho + e_ss*dsig)*gam(dd) + 2.0d0*e_s*dgam(dd)
      enddo
      do a = 1,nPL
         cfa(a) = w*( Bcoef(1)*dPLp(a,1)+Bcoef(2)*dPLp(a,2)+Bcoef(3)*dPLp(a,3) )
      enddo
      do i = 1,nPR
         cfi = w*( Acoef*PRp(i) + Bcoef(1)*dPRp(i,1)+Bcoef(2)*dPRp(i,2)+Bcoef(3)*dPRp(i,3) )
         do a = 1,nPL
            Xloc(a,i,ib) = Xloc(a,i,ib) + PLp(a)*cfi + cfa(a)*PRp(i)
         enddo
      enddo
   enddo
enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) Xproj = Xproj + Xloc
   enddo
end block
deallocate(v0,v1,v2,v3,A0,PLd,PRd,dPLd,dPRd,PLp,PRp,dPLp,dPRp,Yl,Zr,cfa,Xloc)
!$omp end parallel
end subroutine xc_fxc_project_pt

subroutine xc_fxc_project_chk(nConts, nDL, nDR, CDL, CDR, nPL, nPR, CPL, CPR, nbatch, Xd, Ptot, Xproj)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, nDL, nDR, nPL, nPR, nbatch
real(8),intent(in) :: CDL(nConts,nDL), CDR(nConts,nDR)
real(8),intent(in) :: CPL(nConts,nPL), CPR(nConts,nPR)
real(8),intent(in) :: Xd(nDL,nDR,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xproj(nPL,nPR,nbatch)
external :: dgemm
integer,parameter :: NGC = 128
integer :: ic0, ng, k, ig, dd, ib, a, i
real(8) :: rho, sig, drho, dgam(3), dsig, Acoef, Bcoef(3), tau, kk(9), dtau, Tc
logical :: mg
real(8),allocatable :: A1c(:,:,:), dY(:,:,:), T2(:,:), ert(:), est(:), ett(:), tcs(:)
real(8),allocatable :: v0(:),v1(:,:),v2(:,:,:),v3(:,:,:,:)
real(8),allocatable :: V0c(:,:), V1c(:,:,:), A0c(:,:), wg(:), gam(:,:), err(:), ers(:), ess(:), es(:)
real(8),allocatable :: PLd(:,:), PRd(:,:), dPLd(:,:,:), dPRd(:,:,:), PLp(:,:), PRp(:,:), dPLp(:,:,:), dPRp(:,:,:)
real(8),allocatable :: Y(:,:), Z(:,:), cfi(:,:), cfa(:,:), Xloc(:,:,:)
character(len=8) :: envs
integer :: npg

call get_environment_variable("ENGINE_XC_PROJ_PT", envs)
if (trim(envs) .eq. "1") then
   call xc_fxc_project_pt(nConts, nDL, nDR, CDL, CDR, nPL, nPR, CPL, CPR, nbatch, Xd, Ptot, Xproj)
   return
endif
mg = xc_uses_tau()
call build_tables()
if (nConts .ne. nbf) then
   print *,'ERROR: xc_fxc_project nConts/=nbf'; stop 1
endif
Xproj = 0.0d0
npg = ngrids
if (nrg .gt. 0) npg = nrg
!$omp parallel default(shared) &
!$omp&   private(ic0,ng,k,ig,dd,ib,a,i,rho,sig,drho,dgam,dsig,Acoef,Bcoef,v0,v1,v2,v3,V0c,V1c,A0c, &
!$omp&           wg,gam,err,ers,ess,es,PLd,PRd,dPLd,dPRd,PLp,PRp,dPLp,dPRp,Y,Z,cfi,cfa,Xloc, &
!$omp&           kk,tau,dtau,Tc,A1c,dY,T2,ert,est,ett,tcs)
allocate(v0(nConts),v1(nConts,3),v2(nConts,3,3),v3(nConts,3,3,3))
allocate(V0c(nConts,NGC), V1c(nConts,NGC,3), A0c(nConts,NGC))
allocate(wg(NGC), gam(NGC,3), err(NGC), ers(NGC), ess(NGC), es(NGC))
allocate(PLd(NGC,nDL), PRd(NGC,nDR), dPLd(NGC,nDL,3), dPRd(NGC,nDR,3))
allocate(PLp(NGC,nPL), PRp(NGC,nPR), dPLp(NGC,nPL,3), dPRp(NGC,nPR,3))
allocate(Y(NGC,nDL), Z(NGC,nDR), cfi(NGC,nPR), cfa(NGC,nPL), Xloc(nPL,nPR,nbatch))
allocate(ert(NGC), est(NGC), ett(NGC), tcs(NGC))
allocate(A1c(nConts,merge(NGC,1,mg),merge(3,1,mg)), dY(NGC,merge(nDL,1,mg),merge(3,1,mg)), T2(NGC,merge(nPR,1,mg)))
ert = 0.0d0 ; est = 0.0d0 ; ett = 0.0d0
Xloc = 0.0d0
!$omp do schedule(static)
do ic0 = 1, npg, NGC
   ng = min(NGC, npg-ic0+1)
   do k = 1,ng
      ig = ic0 + k - 1
      if (nrg .gt. 0) then
         wg(k) = rg_w(ig)
      else
         wg(k) = Grids(ig)%weight
      endif
      if (wg(k) .eq. 0.0d0) then
         V0c(:,k) = 0.0d0 ; V1c(:,k,:) = 0.0d0
         cycle
      endif
      if (nrg .gt. 0) then
         call eval_bf_derivs(rg_coor(:,ig), v0,v1,v2,v3, 1)
      else
         call eval_bf_derivs(Grids(ig)%coor, v0,v1,v2,v3, 1)
      endif
      V0c(:,k) = v0
      do dd = 1,3
         V1c(:,k,dd) = v1(:,dd)
      enddo
   enddo
   call dgemm('N','N',nConts,ng,nConts,1.0d0,Ptot,nConts,V0c,nConts,0.0d0,A0c,nConts)
   if (mg) then
      do dd = 1,3
         call dgemm('N','N',nConts,ng,nConts,1.0d0,Ptot,nConts,V1c(1,1,dd),nConts,0.0d0,A1c(1,1,dd),nConts)
      enddo
   endif
   do k = 1,ng
      err(k) = 0.0d0 ; ers(k) = 0.0d0 ; ess(k) = 0.0d0 ; es(k) = 0.0d0 ; gam(k,:) = 0.0d0
      ert(k) = 0.0d0 ; est(k) = 0.0d0 ; ett(k) = 0.0d0
      if (wg(k) .eq. 0.0d0) cycle
      rho = sum(V0c(:,k)*A0c(:,k))
      do dd = 1,3
         gam(k,dd) = 2.0d0*sum(V1c(:,k,dd)*A0c(:,k))
      enddo
      sig = gam(k,1)**2 + gam(k,2)**2 + gam(k,3)**2
      tau = 0.0d0
      if (mg) tau = 0.5d0*sum(V1c(:,k,:)*A1c(:,k,:))
      call xc_kern_c(mg, rho, sig, tau, kk)
      err(k)=kk(3); ers(k)=kk(4); ess(k)=kk(5); es(k)=kk(2)
      ert(k)=kk(7); est(k)=kk(8); ett(k)=kk(9)
   enddo
   call dgemm('T','N',ng,nDL,nConts,1.0d0,V0c,nConts,CDL,nConts,0.0d0,PLd,NGC)
   call dgemm('T','N',ng,nDR,nConts,1.0d0,V0c,nConts,CDR,nConts,0.0d0,PRd,NGC)
   call dgemm('T','N',ng,nPL,nConts,1.0d0,V0c,nConts,CPL,nConts,0.0d0,PLp,NGC)
   call dgemm('T','N',ng,nPR,nConts,1.0d0,V0c,nConts,CPR,nConts,0.0d0,PRp,NGC)
   do dd = 1,3
      call dgemm('T','N',ng,nDL,nConts,1.0d0,V1c(1,1,dd),nConts,CDL,nConts,0.0d0,dPLd(1,1,dd),NGC)
      call dgemm('T','N',ng,nDR,nConts,1.0d0,V1c(1,1,dd),nConts,CDR,nConts,0.0d0,dPRd(1,1,dd),NGC)
      call dgemm('T','N',ng,nPL,nConts,1.0d0,V1c(1,1,dd),nConts,CPL,nConts,0.0d0,dPLp(1,1,dd),NGC)
      call dgemm('T','N',ng,nPR,nConts,1.0d0,V1c(1,1,dd),nConts,CPR,nConts,0.0d0,dPRp(1,1,dd),NGC)
   enddo
   do ib = 1,nbatch
      if (mg) then
         do dd = 1,3
            call dgemm('N','T',ng,nDL,nDR,1.0d0,dPRd(1,1,dd),NGC,Xd(1,1,ib),nDL,0.0d0,dY(1,1,dd),NGC)
         enddo
      endif
      call dgemm('N','T',ng,nDL,nDR,1.0d0,PRd,NGC,Xd(1,1,ib),nDL,0.0d0,Y,NGC)
      call dgemm('N','N',ng,nDR,nDL,1.0d0,PLd,NGC,Xd(1,1,ib),nDL,0.0d0,Z,NGC)
      do k = 1,ng
         if (wg(k) .eq. 0.0d0) then
            cfi(k,:) = 0.0d0 ; cfa(k,:) = 0.0d0 ; tcs(k) = 0.0d0
            cycle
         endif
         drho = 0.0d0
         dgam = 0.0d0
         do a = 1,nDL
            drho = drho + PLd(k,a)*Y(k,a)
            dgam(1) = dgam(1) + dPLd(k,a,1)*Y(k,a)
            dgam(2) = dgam(2) + dPLd(k,a,2)*Y(k,a)
            dgam(3) = dgam(3) + dPLd(k,a,3)*Y(k,a)
         enddo
         do i = 1,nDR
            dgam(1) = dgam(1) + dPRd(k,i,1)*Z(k,i)
            dgam(2) = dgam(2) + dPRd(k,i,2)*Z(k,i)
            dgam(3) = dgam(3) + dPRd(k,i,3)*Z(k,i)
         enddo
         drho = 2.0d0*drho
         dgam = 2.0d0*dgam
         dsig = 2.0d0*(gam(k,1)*dgam(1)+gam(k,2)*dgam(2)+gam(k,3)*dgam(3))
         Acoef = err(k)*drho + ers(k)*dsig
         do dd = 1,3
            Bcoef(dd) = 2.0d0*(ers(k)*drho + ess(k)*dsig)*gam(k,dd) + 2.0d0*es(k)*dgam(dd)
         enddo
         if (mg) then
            dtau = 0.0d0
            do dd = 1,3
               do a = 1,nDL
                  dtau = dtau + dPLd(k,a,dd)*dY(k,a,dd)
               enddo
            enddo
            Acoef = Acoef + ert(k)*dtau
            Bcoef = Bcoef + 2.0d0*est(k)*dtau*gam(k,:)
            Tc = ert(k)*drho + est(k)*dsig + ett(k)*dtau
            tcs(k) = 0.5d0*wg(k)*Tc
         endif
         do i = 1,nPR
            cfi(k,i) = wg(k)*( Acoef*PRp(k,i) + Bcoef(1)*dPRp(k,i,1)+Bcoef(2)*dPRp(k,i,2)+Bcoef(3)*dPRp(k,i,3) )
         enddo
         do a = 1,nPL
            cfa(k,a) = wg(k)*( Bcoef(1)*dPLp(k,a,1)+Bcoef(2)*dPLp(k,a,2)+Bcoef(3)*dPLp(k,a,3) )
         enddo
      enddo
      call dgemm('T','N',nPL,nPR,ng,1.0d0,PLp,NGC,cfi,NGC,1.0d0,Xloc(1,1,ib),nPL)
      call dgemm('T','N',nPL,nPR,ng,1.0d0,cfa,NGC,PRp,NGC,1.0d0,Xloc(1,1,ib),nPL)
      if (mg) then
         do dd = 1,3
            do i = 1,nPR
               T2(1:ng,i) = tcs(1:ng)*dPRp(1:ng,i,dd)
            enddo
            call dgemm('T','N',nPL,nPR,ng,1.0d0,dPLp(1,1,dd),NGC,T2,NGC,1.0d0,Xloc(1,1,ib),nPL)
         enddo
      endif
   enddo
enddo
!$omp end do
block
   integer :: mt, nt
   nt = omp_get_num_threads()
   do mt = 0,nt-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) Xproj = Xproj + Xloc
   enddo
end block
deallocate(v0,v1,v2,v3,V0c,V1c,A0c,wg,gam,err,ers,ess,es,PLd,PRd,dPLd,dPRd,PLp,PRp,dPLp,dPRp,Y,Z,cfi,cfa,Xloc, &
           ert,est,ett,tcs,A1c,dY,T2)
!$omp end parallel
end subroutine xc_fxc_project_chk

subroutine xc_fxc_project(nConts, nDL, nDR, CDL, CDR, nPL, nPR, CPL, CPR, nbatch, Xd, Ptot, Xproj)
implicit none
integer,intent(in) :: nConts, nDL, nDR, nPL, nPR, nbatch
real(8),intent(in) :: CDL(nConts,nDL), CDR(nConts,nDR)
real(8),intent(in) :: CPL(nConts,nPL), CPR(nConts,nPR)
real(8),intent(in) :: Xd(nDL,nDR,nbatch)
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Xproj(nPL,nPR,nbatch)
external :: dgemm
real(8),allocatable :: L(:,:,:), Xt(:,:,:)
integer :: b
character(len=8) :: e1, e2
call get_environment_variable("ENGINE_XC_PROJ_GRP", e1)
call get_environment_variable("ENGINE_XC_PROJ_CHK", e2)
if (trim(e1) .eq. "1" .or. trim(e2) .eq. "1" .or. nDR .ne. nPR .or. nPL .le. nPR) then
   call xc_fxc_project_chk(nConts, nDL, nDR, CDL, CDR, nPL, nPR, CPL, CPR, nbatch, Xd, Ptot, Xproj)
   return
endif
if (any(CDR .ne. CPR)) then
   call xc_fxc_project_chk(nConts, nDL, nDR, CDL, CDR, nPL, nPR, CPL, CPR, nbatch, Xd, Ptot, Xproj)
   return
endif
allocate(L(nConts,nDR,nbatch), Xt(nConts,nDR,nbatch))
!$omp parallel do schedule(static)
do b = 1,nbatch
   call dgemm('N','N',nConts,nDR,nDL,1.0d0,CDL,nConts,Xd(1,1,b),nDL,0.0d0,L(1,1,b),nConts)
enddo
!$omp end parallel do
call xc_fxc_c(nConts, nDR, CDR, nbatch, L, Ptot, Xt)
!$omp parallel do schedule(static)
do b = 1,nbatch
   call dgemm('T','N',nPL,nPR,nConts,1.0d0,CPL,nConts,Xt(1,1,b),nConts,0.0d0,Xproj(1,1,b),nPL)
enddo
!$omp end parallel do
deallocate(L, Xt)
end subroutine xc_fxc_project

subroutine xc_resp_grid_build(tier)
integer,intent(in) :: tier
integer :: info, save_tier, ig
if (.not. xcgrid_use_theory_table) return
save_tier = xcgrid_theory_tier
if (tier .eq. save_tier) return
call xcgrid_free_derived()
xcgrid_theory_tier = tier
call xcgrid_theory_apply()
call gridgen(info)
nrg = ngrids
if (allocated(rg_coor)) deallocate(rg_coor, rg_w)
allocate(rg_coor(3,nrg), rg_w(nrg))
do ig = 1,nrg
   rg_coor(:,ig) = Grids(ig)%coor
   rg_w(ig) = Grids(ig)%weight
enddo
call xcgrid_free_derived()
xcgrid_theory_tier = save_tier
call xcgrid_theory_apply()
call gridgen(info)
call GTOeval(info)
print '(A,I0,A,I0,A)', "  CPSCF XC response grid: ", nrg, " points (SCF grid ", ngrids, ")"
end subroutine xc_resp_grid_build

subroutine xc_resp_grid_clear()
nrg = 0
if (allocated(rg_coor)) deallocate(rg_coor, rg_w)
end subroutine xc_resp_grid_clear

end module mod_xc_hess
