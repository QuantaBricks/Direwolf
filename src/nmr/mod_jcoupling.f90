! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_jcoupling: indirect nuclear spin-spin coupling, Fermi-contact term, closed-shell RHF / RKS.

module mod_jcoupling
implicit none
private
public :: jcoupling_fc, jcoupling_pso_sd, jcoupling_dso, gyromagnetic_ratio, jcoupling_to_hz, jcoupling_to_hz_all

contains

function gyromagnetic_ratio(z) result(g)
integer,intent(in) :: z
real(8) :: g
select case (z)
case (1);  g =  2.675221874d8
case (3);  g =  1.039790d8
case (6);  g =  6.728284d7
case (7);  g = -2.7126180d7
case (8);  g = -3.62808d7
case (9);  g =  2.518148d8
case (14); g = -5.3190d7
case (15); g =  1.08394d8
case (16); g =  2.0517d7
case (17); g =  2.6242d7
case default; g = 0.0d0
end select
end function gyromagnetic_ratio

subroutine jc_solve(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, mode, nrhs, H, X, niter, meta)
use mod_cpscf, only: cpscf_apb_rhf, cpscf_apb_amb_rhf, cpscf_response_begin, cpscf_response_end, cpscf_triplet
use mod_exchange, only: RS_omega
implicit none
integer,intent(in) :: nConts, nocc, mode, nrhs
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts), Pa(nConts,nConts), Pb(nConts,nConts), HF_exchange_frac
real(8),intent(in) :: H(nConts-nocc,nocc,nrhs)
real(8),intent(out) :: X(nConts-nocc,nocc,nrhs)
integer,intent(out) :: niter
logical,intent(in),optional :: meta
logical :: use_meta
integer :: nvir, nt, k, a, i, it, m, m0, mmax, info, lw
real(8),allocatable :: R(:,:), T(:,:), dl(:), V(:,:), AV(:,:), Hf(:,:), Xf(:,:)
real(8),allocatable :: Mr(:,:), Ms(:,:), Bred(:,:), Y(:,:), ev(:), wk(:), bn(:)
integer,allocatable :: ipv(:)
logical,allocatable :: act(:)
logical :: exx, plain
real(8),parameter :: TOL = 1.0d-7

nvir = nConts - nocc
nt = nvir*nocc
exx = (HF_exchange_frac /= 0.0d0 .or. RS_omega > 0.0d0)
use_meta = .false.
if (present(meta)) use_meta = meta .and. mode == 2
plain = (mode == 2 .and. .not. exx .and. .not. use_meta)
allocate(dl(nt), Hf(nt,nrhs), Xf(nt,nrhs))
do i = 1,nocc
   do a = 1,nvir
      dl(a+(i-1)*nvir) = eLev_a(nocc+a) - eLev_a(i)
   enddo
enddo
Hf = reshape(H, (/nt,nrhs/))
niter = 0
if (plain) then
   do k = 1,nrhs
      Xf(:,k) = Hf(:,k)/dl
   enddo
   X = reshape(Xf, (/nvir,nocc,nrhs/))
   return
endif

mmax = min(nt, max(60, 12*nrhs))
allocate(V(nt,mmax), AV(nt,mmax), R(nt,nrhs), T(nt,nrhs), bn(nrhs), act(nrhs))
do k = 1,nrhs
   bn(k) = sqrt(sum(Hf(:,k)**2))
   T(:,k) = Hf(:,k)/dl
enddo
act = bn > 1.0d-14
m = 0 ; m0 = 0
call add_directions()
cpscf_triplet = (mode == 1)
call cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac, scf_grid=.true.)
Xf = 0.0d0
do it = 1,60
   if (m == m0) exit
   niter = it
   call apply_new(m0+1, m)
   m0 = m
   allocate(Mr(m,m), Bred(m,nrhs), Y(m,nrhs), ipv(m), Ms(m,m), ev(m), wk(1))
   call dgemm('T','N',m,m,nt,1.0d0,V,nt,AV,nt,0.0d0,Mr,m)
   call dgemm('T','N',m,nrhs,nt,1.0d0,V,nt,Hf,nt,0.0d0,Bred,m)
   Ms = 0.5d0*(Mr + transpose(Mr))
   call dsyev('N','U',m,Ms,m,ev,wk,-1,info)
   lw = int(wk(1)) ; deallocate(wk) ; allocate(wk(lw))
   call dsyev('N','U',m,Ms,m,ev,wk,lw,info)
   if (info /= 0) error stop 'J coupling: diagonalization of the reduced Hessian failed'
   if (ev(1) <= 0.0d0) then
      call cpscf_response_end()
      cpscf_triplet = .false.
      if (mode == 1) error stop 'J coupling: triplet instability of the reference (Hartree-Fock on a multiple bond): '// &
                                'use a hybrid functional'
      error stop 'J coupling: singlet instability of the reference (negative eigenvalue of A-B)'
   endif
   Y = Bred
   call dgesv(m,nrhs,Mr,m,ipv,Y,m,info)
   if (info /= 0) error stop 'J coupling: the reduced response system is singular'
   call dgemm('N','N',nt,nrhs,m,1.0d0,V,nt,Y,m,0.0d0,Xf,nt)
   call dgemm('N','N',nt,nrhs,m,1.0d0,AV,nt,Y,m,0.0d0,R,nt)
   R = Hf - R
   deallocate(Mr, Bred, Y, ipv, Ms, ev, wk)
   do k = 1,nrhs
      act(k) = bn(k) > 1.0d-14
      if (act(k)) act(k) = sqrt(sum(R(:,k)**2)) > TOL*bn(k)
      T(:,k) = R(:,k)/dl
   enddo
   if (.not. any(act)) exit
   call add_directions()
   if (m == m0) exit
enddo
call cpscf_response_end()
cpscf_triplet = .false.
if (any(act)) error stop 'J coupling: the response did not converge'
X = reshape(Xf, (/nvir,nocc,nrhs/))

contains

subroutine add_directions()
integer :: kk, jj, q
real(8) :: nrm
do kk = 1,nrhs
   if (.not. act(kk)) cycle
   nrm = sqrt(sum(T(:,kk)**2))
   if (nrm == 0.0d0) cycle
   T(:,kk) = T(:,kk)/nrm
   do q = 1,2
      do jj = 1,m
         T(:,kk) = T(:,kk) - sum(V(:,jj)*T(:,kk))*V(:,jj)
      enddo
   enddo
   nrm = sqrt(sum(T(:,kk)**2))
   if (nrm < 1.0d-6) cycle
   if (m >= mmax) error stop 'J coupling: the response subspace is full'
   m = m + 1
   V(:,m) = T(:,kk)/nrm
enddo
end subroutine add_directions

subroutine add_current_kernel(nb, Xr, Qr)
use mod_nmr, only: nmr_xc_resp
integer,intent(in) :: nb
real(8),intent(in) :: Xr(nvir,nocc,nb)
real(8),intent(inout) :: Qr(nvir,nocc,nb)
real(8),allocatable :: dd(:,:,:), kl(:,:,:), t1(:,:), t2(:,:), Cvm(:,:), Com(:,:)
integer :: ix, p, q
allocate(dd(nConts,nConts,nb), kl(nConts,nConts,nb), t1(nConts,nocc), t2(nConts,nConts), Cvm(nConts,nvir), Com(nConts,nocc))
Cvm = C_a(:,nocc+1:nConts) ; Com = C_a(:,1:nocc)
do ix = 1,nb
   call dgemm('N','N',nConts,nocc,nvir,1.0d0,Cvm,nConts,Xr(1,1,ix),nvir,0.0d0,t1,nConts)
   call dgemm('N','T',nConts,nConts,nocc,1.0d0,t1,nConts,Com,nConts,0.0d0,t2,nConts)
   do q = 1,nConts
      do p = 1,nConts
         dd(p,q,ix) = 2.0d0*(t2(p,q) - t2(q,p))
      enddo
   enddo
enddo
call nmr_xc_resp(nConts, nb, dd, kl, 1)
do ix = 1,nb
   call dgemm('N','N',nConts,nocc,nConts,1.0d0,kl(1,1,ix),nConts,Com,nConts,0.0d0,t1,nConts)
   call dgemm('T','N',nvir,nocc,nConts,-1.0d0,Cvm,nConts,t1,nConts,1.0d0,Qr(1,1,ix),nvir)
enddo
deallocate(dd, kl, t1, t2, Cvm, Com)
end subroutine add_current_kernel

subroutine apply_new(j1, j2)
integer,intent(in) :: j1, j2
real(8),allocatable :: Xr(:,:,:), Qr(:,:,:), Mm(:,:,:)
integer :: nb, jx
nb = j2 - j1 + 1
allocate(Xr(nvir,nocc,nb), Qr(nvir,nocc,nb), Mm(nvir,nocc,1))
do jx = 1,nb
   Xr(:,:,jx) = reshape(V(:,j1+jx-1), (/nvir,nocc/))
enddo
if (mode == 2) then
   call cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, 0, Mm, Mm, nb, Xr, Qr)
   if (use_meta) call add_current_kernel(nb, Xr, Qr)
elseif (exx) then
   call cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, nb, Xr, Qr, 0, Mm, Mm)
else
   call cpscf_apb_rhf(nConts, nocc, nvir, nb, C_a, eLev_a, HF_exchange_frac, Xr, Qr)
endif
do jx = 1,nb
   AV(:,j1+jx-1) = reshape(Qr(:,:,jx), (/nt/))
enddo
deallocate(Xr, Qr, Mm)
end subroutine apply_new

end subroutine jc_solve

subroutine jcoupling_fc(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nat, coor_bohr, S, niter)
use mod_cube, only: cube_screen_begin, cube_screen_end, cube_ao_values
implicit none
integer,intent(in) :: nConts, nocc, nat
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts), Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac, coor_bohr(3,nat)
real(8),intent(out) :: S(nat,nat)
integer,intent(out) :: niter
integer :: nvir, k, a, i
real(8),allocatable :: phi(:,:), vmo(:,:), H(:,:,:), X(:,:,:), Phit(:,:)
logical :: own_screen

nvir = nConts - nocc
allocate(phi(nConts,nat), vmo(nConts,nat), Phit(nat,nConts), H(nvir,nocc,nat), X(nvir,nocc,nat))
call cube_screen_begin(own_screen)
call cube_ao_values(nat, coor_bohr, Phit)
call cube_screen_end(own_screen)
phi = transpose(Phit)
call dgemm('T','N',nConts,nat,nConts,1.0d0,C_a,nConts,phi,nConts,0.0d0,vmo,nConts)
do k = 1,nat
   do i = 1,nocc
      do a = 1,nvir
         H(a,i,k) = vmo(nocc+a,k)*vmo(i,k)
      enddo
   enddo
enddo
call jc_solve(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, 1, nat, H, X, niter)
call dgemm('T','N',nat,nat,nvir*nocc,1.0d0,H,nvir*nocc,X,nvir*nocc,0.0d0,S,nat)
S = 0.5d0*(S + transpose(S))
end subroutine jcoupling_fc

subroutine jcoupling_pso_sd(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nat, coor_bohr, &
                            want_pso, want_sd, Spso, Ssd, niter)
use mod_cube, only: cube_screen_begin, cube_screen_end, cube_ao_values
use mod_integrals, only: integrals_nmr_nucleus, integrals_hessian_nuc_atom
use mod_xc, only: xc_uses_tau
use mod_nmr, only: nmr_xc_giao
implicit none
integer,intent(in) :: nConts, nocc, nat
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts), Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac, coor_bohr(3,nat)
logical :: meta
logical,intent(in) :: want_pso, want_sd
real(8),intent(out) :: Spso(nat,nat), Ssd(nat,nat)
integer,intent(out) :: niter
integer,parameter :: AB(2,6) = reshape((/1,1, 2,2, 3,3, 1,2, 1,3, 2,3/), (/2,6/))
real(8),parameter :: WT(6) = (/1.0d0,1.0d0,1.0d0,2.0d0,2.0d0,2.0d0/), pi = 3.14159265358979d0
integer :: nvir, nt, k, l, c, i, a, nj
real(8),allocatable :: PSO(:,:,:), A11(:,:,:), A01(:,:,:), CC(:,:,:), CB(:,:,:), CK(:,:,:), Hp(:,:,:), Hs(:,:,:)
real(8),allocatable :: Xp(:,:,:), Xs(:,:,:), Mao(:,:), Tm(:,:), Phit(:,:), phi(:,:), Cv(:,:), Co(:,:)
logical :: own_screen

meta = want_pso .and. xc_uses_tau()
if (meta) then
   block
      real(8),allocatable :: Vx(:,:,:,:), Dt(:,:)
      allocate(Vx(nConts,nConts,3,1), Dt(nConts,nConts))
      Dt = Pa + Pb
      call nmr_xc_giao(nConts, 1, Dt, 1, coor_bohr(:,1), Vx)
   end block
endif
nvir = nConts - nocc
nt = nvir*nocc
niter = 0
Spso = 0.0d0 ; Ssd = 0.0d0
allocate(Mao(nConts,nConts), Tm(nvir,nConts), Cv(nConts,nvir), Co(nConts,nocc), Phit(1,nConts), phi(nConts,1))
Cv = C_a(:,nocc+1:nConts) ; Co = C_a(:,1:nocc)
if (want_pso) allocate(Hp(nvir,nocc,3*nat), Xp(nvir,nocc,3*nat))
if (want_sd) allocate(Hs(nvir,nocc,6*nat), Xs(nvir,nocc,6*nat))
if (want_pso) allocate(PSO(nConts,nConts,3), A11(nConts,nConts,9), A01(nConts,nConts,9))
if (want_sd) allocate(CC(nConts,nConts,9), CB(nConts,nConts,9), CK(nConts,nConts,9))
do k = 1,nat
   if (want_pso) then
      call integrals_nmr_nucleus(nConts, coor_bohr(:,k), PSO, A11, A01)
      do c = 1,3
         call mo_vo(PSO(:,:,c), Hp(:,:,3*(k-1)+c))
      enddo
   endif
   if (want_sd) then
      call integrals_hessian_nuc_atom(coor_bohr(:,k), nConts, CC, CB, CK)
      call cube_screen_begin(own_screen)
      call cube_ao_values(1, coor_bohr(:,k), Phit)
      call cube_screen_end(own_screen)
      phi(:,1) = Phit(1,:)
      do c = 1,6
         Mao = CC(:,:,3*(AB(1,c)-1)+AB(2,c))
         Mao = 0.5d0*(Mao + transpose(Mao))
         if (AB(1,c) == AB(2,c)) then
            do l = 1,nConts
               do i = 1,nConts
                  Mao(i,l) = Mao(i,l) + 4.0d0*pi/3.0d0*phi(i,1)*phi(l,1)
               enddo
            enddo
         endif
         call mo_vo(Mao, Hs(:,:,6*(k-1)+c))
      enddo
   endif
enddo
if (want_pso) then
   call jc_solve(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, 2, 3*nat, Hp, Xp, nj, meta)
   niter = max(niter, nj)
   do k = 1,nat
      do l = 1,nat
         do c = 1,3
            Spso(l,k) = Spso(l,k) + sum(Hp(:,:,3*(l-1)+c)*Xp(:,:,3*(k-1)+c))
         enddo
      enddo
   enddo
   Spso = 0.5d0*(Spso + transpose(Spso))
endif
if (want_sd) then
   call jc_solve(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, 1, 6*nat, Hs, Xs, nj)
   niter = max(niter, nj)
   do k = 1,nat
      do l = 1,nat
         do c = 1,6
            Ssd(l,k) = Ssd(l,k) + WT(c)*sum(Hs(:,:,6*(l-1)+c)*Xs(:,:,6*(k-1)+c))
         enddo
      enddo
   enddo
   Ssd = 0.5d0*(Ssd + transpose(Ssd))
endif

contains

subroutine mo_vo(M, out)
real(8),intent(in) :: M(nConts,nConts)
real(8),intent(out) :: out(nvir,nocc)
call dgemm('T','N',nvir,nConts,nConts,1.0d0,Cv,nConts,M,nConts,0.0d0,Tm,nvir)
call dgemm('N','N',nvir,nocc,nConts,1.0d0,Tm,nvir,Co,nConts,0.0d0,out,nvir)
end subroutine mo_vo

end subroutine jcoupling_pso_sd

subroutine jcoupling_dso(nConts, Pa, Pb, nat, coor_bohr, Dso)
use GRID_info, only: Grids, ngrids
use mod_cube, only: cube_screen_begin, cube_screen_end
implicit none
integer,intent(in) :: nConts, nat
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts), coor_bohr(3,nat)
real(8),intent(out) :: Dso(nat,nat)
interface
   subroutine GTOeval_point(coor_grid, val0, val1)
   use MOL_info, only: nconts
   real(8),intent(in) :: coor_grid(3)
   real(8),intent(out) :: val0(nconts), val1(nconts,3)
   end subroutine GTOeval_point
   subroutine gridgen(info)
   integer :: info
   end subroutine gridgen
end interface
integer,parameter :: NBT = 128
real(8),allocatable :: P(:,:), Phi(:,:), T(:,:), v0(:), v1(:,:), Dl(:,:), rho(:), wt(:), rr(:,:), rk(:,:), ri(:)
integer :: nblk, ib, ip0, nbn, k, kk, l, info_g
real(8) :: dd
logical :: own

if (.not. allocated(Grids) .or. ngrids <= 0) then
   call gridgen(info_g)
   if (.not. allocated(Grids) .or. ngrids <= 0) error stop 'J coupling (DSO): no integration grid'
endif
allocate(P(nConts,nConts))
P = Pa + Pb
Dso = 0.0d0
call cube_screen_begin(own)
nblk = (ngrids + NBT - 1)/NBT
!$omp parallel private(Phi, T, v0, v1, Dl, rho, wt, rr, rk, ri, ib, ip0, nbn, k, kk, l, dd)
allocate(Phi(NBT,nConts), T(NBT,nConts), v0(nConts), v1(nConts,3), Dl(nat,nat), rho(NBT), wt(NBT), rr(NBT,3), rk(3,nat), ri(nat))
Dl = 0.0d0
!$omp do schedule(static,1)
do ib = 1,nblk
   ip0 = 1 + (ib-1)*NBT
   nbn = min(NBT, ngrids - ip0 + 1)
   do k = 1,nbn
      call GTOeval_point(Grids(ip0+k-1)%coor, v0, v1)
      Phi(k,:) = v0
      wt(k) = Grids(ip0+k-1)%weight
      rr(k,:) = Grids(ip0+k-1)%coor
   enddo
   call dgemm('N','N',nbn,nConts,nConts,1.0d0,Phi,NBT,P,nConts,0.0d0,T,NBT)
   do k = 1,nbn
      rho(k) = sum(T(k,:)*Phi(k,:))
      do kk = 1,nat
         rk(:,kk) = rr(k,:) - coor_bohr(:,kk)
         ri(kk) = sqrt(sum(rk(:,kk)**2))
      enddo
      do kk = 1,nat
         if (ri(kk) < 1.0d-8) cycle
         do l = kk+1,nat
            if (ri(l) < 1.0d-8) cycle
            dd = sum(rk(:,kk)*rk(:,l))/(ri(kk)**3*ri(l)**3)
            Dl(kk,l) = Dl(kk,l) + wt(k)*rho(k)*dd
         enddo
      enddo
   enddo
enddo
!$omp end do
!$omp critical
Dso = Dso + Dl
!$omp end critical
deallocate(Phi, T, v0, v1, Dl, rho, wt, rr, rk, ri)
!$omp end parallel
call cube_screen_end(own)
Dso = Dso + transpose(Dso)
end subroutine jcoupling_dso

subroutine jcoupling_to_hz_all(nat, znum, Sfc, Ssd, Spso, Dso, Jfc, Jsd, Jpso, Jdso, Jtot, Kred)
integer,intent(in) :: nat, znum(nat)
real(8),intent(in) :: Sfc(nat,nat), Ssd(nat,nat), Spso(nat,nat), Dso(nat,nat)
real(8),intent(out) :: Jfc(nat,nat), Jsd(nat,nat), Jpso(nat,nat), Jdso(nat,nat), Jtot(nat,nat), Kred(nat,nat)
real(8),parameter :: mu0 = 1.25663706212d-6, ge = 2.0d0, muB = 9.2740100783d-24
real(8),parameter :: hbar = 1.054571817d-34, hpl = 6.62607015d-34, Eh = 4.3597447222071d-18, a0 = 5.29177210903d-11
real(8),parameter :: qe = 1.602176634d-19, me = 9.1093837015d-31, pi = 3.14159265358979d0
real(8) :: Pi_(nat), Pd(nat), Cp(nat), g(nat)
integer :: k, l
do k = 1,nat
   g(k) = gyromagnetic_ratio(znum(k))
   Pi_(k) = 2.0d0*mu0/3.0d0*ge*muB*g(k)*hbar
   Pd(k) = mu0/(4.0d0*pi)*ge*muB*g(k)*hbar
   Cp(k) = qe*mu0*g(k)*hbar**2/(4.0d0*pi*me*a0**3)
enddo
Jfc = 0.0d0 ; Jsd = 0.0d0 ; Jpso = 0.0d0 ; Jdso = 0.0d0 ; Jtot = 0.0d0 ; Kred = 0.0d0
do k = 1,nat
   do l = 1,nat
      if (k == l .or. g(k) == 0.0d0 .or. g(l) == 0.0d0) cycle
      Jfc(k,l) = -(Pi_(k)*Pi_(l)/hpl)*Sfc(k,l)/(Eh*a0**6)
      Jsd(k,l) = -(Pd(k)*Pd(l)/(3.0d0*hpl))*Ssd(k,l)/(Eh*a0**6)
      Jpso(k,l) = -4.0d0*Cp(k)*Cp(l)*Spso(k,l)/(3.0d0*hpl*Eh)
      Jdso(k,l) = 2.0d0/(3.0d0*hpl)*(qe**2/me)*(mu0/(4.0d0*pi))**2*g(k)*g(l)*hbar**2*Dso(k,l)/a0**4
      Jtot(k,l) = Jfc(k,l) + Jsd(k,l) + Jpso(k,l) + Jdso(k,l)
      Kred(k,l) = 4.0d0*pi*pi*Jtot(k,l)/(hpl*g(k)*g(l))*1.0d-19
   enddo
enddo
end subroutine jcoupling_to_hz_all

subroutine jcoupling_to_hz(nat, znum, S, J, Kred)
integer,intent(in) :: nat, znum(nat)
real(8),intent(in) :: S(nat,nat)
real(8),intent(out) :: J(nat,nat), Kred(nat,nat)
real(8),parameter :: mu0 = 1.25663706212d-6, ge = 2.0d0, muB = 9.2740100783d-24
real(8),parameter :: hbar = 1.054571817d-34, hpl = 6.62607015d-34, Eh = 4.3597447222071d-18, a0 = 5.29177210903d-11
real(8),parameter :: pi = 3.14159265358979d0
real(8) :: Pi_(nat), g(nat)
integer :: k, l
do k = 1,nat
   g(k) = gyromagnetic_ratio(znum(k))
   Pi_(k) = 2.0d0*mu0/3.0d0*ge*muB*g(k)*hbar
enddo
J = 0.0d0 ; Kred = 0.0d0
do k = 1,nat
   do l = 1,nat
      if (k == l .or. g(k) == 0.0d0 .or. g(l) == 0.0d0) cycle
      J(k,l) = -(Pi_(k)*Pi_(l)/hpl)*S(k,l)/(Eh*a0**6)
      Kred(k,l) = 4.0d0*pi*pi*J(k,l)/(hpl*g(k)*g(l))*1.0d-19
   enddo
enddo
end subroutine jcoupling_to_hz

end module mod_jcoupling
