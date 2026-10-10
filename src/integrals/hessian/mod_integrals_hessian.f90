! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! 1-electron analytic second-derivative (Hessian) integral kernels.

submodule (mod_integrals) hessian_impl
implicit none
contains

module subroutine integrals_hessian_1e_shell(nConts, s1aa, s1ab, h1aa, h1ab)
implicit none
integer,intent(in) :: nConts
real(8),intent(out) :: s1aa(nConts,nConts,9), s1ab(nConts,nConts,9)
real(8),intent(out) :: h1aa(nConts,nConts,9), h1ab(nConts,nConts,9)

integer :: i,j,di,dj,e1,e2,a,b
integer :: shls(2)
real(8),allocatable :: buf_ovlp_aa(:,:,:), buf_ovlp_ab(:,:,:)
real(8),allocatable :: buf_kin_aa(:,:,:), buf_kin_ab(:,:,:)
real(8),allocatable :: buf_nuc_aa(:,:,:), buf_nuc_ab(:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: si_off, num_off

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

!$omp parallel private(i,j,shls,di,dj,e1,e2,a,b) &
!$omp&   private(buf_ovlp_aa,buf_ovlp_ab,buf_kin_aa,buf_kin_ab,buf_nuc_aa,buf_nuc_ab)
!$omp do schedule(dynamic) collapse(2)
do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i
      shls(2) = j
      di = shell_dim(i)
      dj = shell_dim(j)
      allocate(buf_ovlp_aa(di,dj,9), buf_ovlp_ab(di,dj,9))
      allocate(buf_kin_aa(di,dj,9), buf_kin_ab(di,dj,9))
      allocate(buf_nuc_aa(di,dj,9), buf_nuc_ab(di,dj,9))
      call ovlp1e_ipip_engine(buf_ovlp_aa, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      call ovlp1e_ipovlpip_engine(buf_ovlp_ab, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      call kin1e_ipip_engine(buf_kin_aa, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      call kin1e_ipkinip_engine(buf_kin_ab, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      call nuc1e_ipip_engine(buf_nuc_aa, shls, atm_nuc, nAtoms_nuc, bas, nBases, env_nuc, 0_8)
      call nuc1e_ipnucip_engine(buf_nuc_ab, shls, atm_nuc, nAtoms_nuc, bas, nBases, env_nuc, 0_8)

      if (necpbas .gt. 0) then
         block
            real(8),allocatable :: be(:,:,:)
            integer :: dims(2)
            allocate(be(di,dj,9))
            dims(1) = di ; dims(2) = dj
            call ecp1e_ipipnuc_engine(be, dims, shls, atm, size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
            buf_nuc_aa = buf_nuc_aa + be
            call ecp1e_ipnucip_engine(be, dims, shls, atm, size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
            buf_nuc_ab = buf_nuc_ab + be
            deallocate(be)
         end block
      endif
      do a = 1,di
         do b = 1,dj
            e1 = ao_offset(i)+a
            e2 = ao_offset(j)+b
            s1aa(e1,e2,:) = buf_ovlp_aa(a,b,:)*NorVEC(e1)*NorVEC(e2)
            s1ab(e1,e2,:) = buf_ovlp_ab(a,b,:)*NorVEC(e1)*NorVEC(e2)
            h1aa(e1,e2,:) = (buf_kin_aa(a,b,:)+buf_nuc_aa(a,b,:))*NorVEC(e1)*NorVEC(e2)
            h1ab(e1,e2,:) = (buf_kin_ab(a,b,:)+buf_nuc_ab(a,b,:))*NorVEC(e1)*NorVEC(e2)
         enddo
      enddo
      deallocate(buf_ovlp_aa, buf_ovlp_ab, buf_kin_aa, buf_kin_ab, buf_nuc_aa, buf_nuc_ab)
   enddo
enddo
!$omp end do
!$omp end parallel

deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_1e_shell

module subroutine integrals_hessian_nuc_atom(coor_bohr, nConts, CC, CB, CK)
implicit none
real(8),intent(in) :: coor_bohr(3)
integer,intent(in) :: nConts
real(8),intent(out) :: CC(nConts,nConts,9), CB(nConts,nConts,9), CK(nConts,nConts,9)

integer :: i,j,di,dj,e1,e2,a,b,c
integer :: shls(2)
real(8),allocatable :: buf_aa(:,:,:), buf_ab(:,:,:), buf_kk(:,:,:)
real(8),allocatable :: AA(:,:,:), AB(:,:,:), KK(:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: si_off, num_off

allocate(AA(nConts,nConts,9), AB(nConts,nConts,9), KK(nConts,nConts,9))
call integrals_hessian_nuc_atom_raw(coor_bohr, nConts, AA, AB)

do e1 = 1,nConts
   do e2 = 1,nConts
      KK(e1,e2,:) = AA(e2,e1,:)
   enddo
enddo

do c = 1,9
   CB(:,:,c) = -(AA(:,:,c)+AB(:,:,c))
   CK(:,:,c) = -(AB(:,:,c)+KK(:,:,c))
   CC(:,:,c) = AA(:,:,c) + 2.0d0*AB(:,:,c) + KK(:,:,c)
enddo

deallocate(AA, AB, KK)
end subroutine integrals_hessian_nuc_atom

module subroutine integrals_hessian_ecp_atom_raw(iatm, nConts, AAe, ABe, has)
implicit none
integer,intent(in) :: iatm, nConts
real(8),intent(out) :: AAe(nConts,nConts,9), ABe(nConts,nConts,9)
logical,intent(out) :: has
integer :: i, j, di, dj, e1, e2, a, b, dims(2), shls(2), num_off
real(8),allocatable :: ba(:,:,:), bb(:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
AAe = 0.0d0 ; ABe = 0.0d0
has = .false.
do i = nBases+1, nBases+necpbas
   if (basECP(1,i)+1 .eq. iatm) has = .true.
enddo
if (.not. has) return
envECP(18) = dble(iatm-1)
allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do i = 0,nBases-1
   ao_offset(i) = num_off
   shell_dim(i) = cgto_engine(i, bas)
   num_off = num_off + shell_dim(i)
enddo
!$omp parallel do schedule(dynamic) private(i,j,di,dj,e1,e2,a,b,dims,shls,ba,bb)
do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i ; shls(2) = j
      di = shell_dim(i) ; dj = shell_dim(j)
      dims(1) = di ; dims(2) = dj
      allocate(ba(di,dj,9), bb(di,dj,9))
      call ecp1e_ipiprinv_engine(ba, dims, shls, atm, size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
      call ecp1e_iprinvip_engine(bb, dims, shls, atm, size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
      do a = 1,di
         do b = 1,dj
            e1 = ao_offset(i)+a ; e2 = ao_offset(j)+b
            AAe(e1,e2,:) = ba(a,b,:)*NorVEC(e1)*NorVEC(e2)
            ABe(e1,e2,:) = bb(a,b,:)*NorVEC(e1)*NorVEC(e2)
         enddo
      enddo
      deallocate(ba, bb)
   enddo
enddo
!$omp end parallel do
deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_ecp_atom_raw

module subroutine integrals_ecp_iprinv_atom(iatm, nConts, Ecp, has)
implicit none
integer,intent(in) :: iatm, nConts
real(8),intent(out) :: Ecp(nConts,nConts,3)
logical,intent(out) :: has
integer :: i, j, di, dj, e1, e2, a, b, dims(2), shls(2), num_off
real(8),allocatable :: ba(:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
Ecp = 0.0d0
has = .false.
do i = nBases+1, nBases+necpbas
   if (basECP(1,i)+1 .eq. iatm) has = .true.
enddo
if (.not. has) return
envECP(18) = dble(iatm-1)
allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do i = 0,nBases-1
   ao_offset(i) = num_off
   shell_dim(i) = cgto_engine(i, bas)
   num_off = num_off + shell_dim(i)
enddo
!$omp parallel do schedule(dynamic) private(i,j,di,dj,e1,e2,a,b,dims,shls,ba)
do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i ; shls(2) = j
      di = shell_dim(i) ; dj = shell_dim(j)
      dims(1) = di ; dims(2) = dj
      allocate(ba(di,dj,3))
      call ecp1e_iprinv_engine(ba, dims, shls, atm, size(atm,2), basECP, nBases+necpbas, envECP, 0_8, 0_8)
      do a = 1,di
         do b = 1,dj
            e1 = ao_offset(i)+a ; e2 = ao_offset(j)+b
            Ecp(e1,e2,:) = ba(a,b,:)*NorVEC(e1)*NorVEC(e2)
         enddo
      enddo
      deallocate(ba)
   enddo
enddo
!$omp end parallel do
deallocate(ao_offset, shell_dim)
end subroutine integrals_ecp_iprinv_atom

module subroutine integrals_hessian_nuc_atom_raw(coor_bohr, nConts, AA, AB)
implicit none
real(8),intent(in) :: coor_bohr(3)
integer,intent(in) :: nConts
real(8),intent(out) :: AA(nConts,nConts,9), AB(nConts,nConts,9)

integer :: i,j,di,dj,e1,e2,a,b
integer :: shls(2)
real(8),allocatable :: buf_aa(:,:,:), buf_ab(:,:,:), env_t(:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: si_off, num_off

env(5) = coor_bohr(1)
env(6) = coor_bohr(2)
env(7) = coor_bohr(3)

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

!$omp parallel default(shared) private(i,j,di,dj,e1,e2,a,b,shls,buf_aa,buf_ab,env_t)
allocate(env_t(size(env)))
env_t = env
!$omp do schedule(dynamic) collapse(2)
do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i
      shls(2) = j
      di = shell_dim(i)
      dj = shell_dim(j)
      allocate(buf_aa(di,dj,9), buf_ab(di,dj,9))
      call rinv1e_ipiprinv_engine(buf_aa, shls, atm, size(atm,2), bas, nBases, env_t, 0_8)
      call rinv1e_iprinvip_engine(buf_ab, shls, atm, size(atm,2), bas, nBases, env_t, 0_8)
      do a = 1,di
         do b = 1,dj
            e1 = ao_offset(i)+a
            e2 = ao_offset(j)+b
            AA(e1,e2,:) = buf_aa(a,b,:)*NorVEC(e1)*NorVEC(e2)
            AB(e1,e2,:) = buf_ab(a,b,:)*NorVEC(e1)*NorVEC(e2)
         enddo
      enddo
      deallocate(buf_aa, buf_ab)
   enddo
enddo
!$omp end do
deallocate(env_t)
!$omp end parallel

deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_nuc_atom_raw

module subroutine integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e, omega)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: d2E_2e(3,natoms,3,natoms)
real(8),intent(in),optional :: omega
real(8) :: jscale

integer :: ish,jsh,ksh,lsh,di,dj,dk,dl
integer :: atI,atJ,atK,atL
integer :: a,b,c,d,e1,e2,e3,e4,nq
integer :: shlsA(4), shlsB(4), shlsC(4)
real(8),allocatable :: buf1(:), buf2(:), buf3(:)
real(8),allocatable :: coefA(:), coefB(:), coefC(:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: si_off, num_off, maxdim
real(8) :: coefv, pair_bound
real(8) :: tk0, t_kern, t_coef, t_contr
integer(8) :: opt_ipip1, opt_ipvip1, opt_ip1ip2
integer(8) :: nquad_done
real(8) :: Mii(3,3), Mij(3,3), Mik(3,3), Mil(3,3)
real(8) :: Mji(3,3), Mjj(3,3), Mjk(3,3), Mjl(3,3)
real(8) :: Mki(3,3), Mkj(3,3), Mkk(3,3), Mkl(3,3)
real(8) :: Mli(3,3), Mlj(3,3), Mlk(3,3), Mll(3,3)
real(8) :: omp_get_wtime
external :: omp_get_wtime
real(8) :: tm0, tm_pre, tm_opt, tm_par, tm_post
logical :: tm_on
character(len=4) :: tm_env
real(8) :: sum_ii(9), sum_ij(9), sum_ik(9), sum_il(9)
real(8) :: sum_jj(9), sum_ji(9), sum_jk(9), sum_jl(9)
real(8) :: sum_kk(9), sum_kl(9), sum_ki(9), sum_kj(9)
real(8) :: sum_ll(9), sum_lk(9), sum_li(9), sum_lj(9)
real(8),allocatable :: dmax_shell(:,:)
real(8),allocatable :: Ptot(:,:)
real(8) :: hess_2e_cutoff
integer :: ij, kl, npair, nperm, iperm, jperm
integer :: sh(4), ao(4), shtup(4), useperm(4,8)
integer,allocatable :: pair_i(:), pair_j(:)
logical :: dup
integer,parameter :: PERMMAP(4,8) = reshape( (/ &
   1,2,3,4,   2,1,3,4,   1,2,4,3,   2,1,4,3,  &
   3,4,1,2,   4,3,1,2,   3,4,2,1,   4,3,2,1 /), (/4,8/) )

call get_environment_variable("ENGINE_HESS2E_TIMING", tm_env)
tm_on = (trim(tm_env) .eq. "1")
tm0 = omp_get_wtime()
allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo
maxdim = maxval(shell_dim)

npair = nBases*(nBases+1)/2
allocate(pair_i(0:npair-1), pair_j(0:npair-1))
ij = 0
do ish = 0,nBases-1
   do jsh = 0,ish
      pair_i(ij) = ish
      pair_j(ij) = jsh
      ij = ij + 1
   enddo
enddo

allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts, abs(Pa)+abs(Pb), dmax_shell)
hess_2e_cutoff = 1.0d-12

allocate(Ptot(nConts,nConts))
Ptot = Pa + Pb

d2E_2e = 0.0d0
t_kern = 0.0d0; t_coef = 0.0d0; t_contr = 0.0d0; nquad_done = 0

jscale = 1.0d0
if (present(omega)) then
   if (omega > 0.0d0) then
      jscale = 0.0d0
      env(9) = omega
   endif
endif
tm_pre = omp_get_wtime() - tm0 ; tm0 = omp_get_wtime()
call twoe_ipip1_optimizer_engine(opt_ipip1, atm, size(atm,2), bas, nBases, env)
call twoe_ipvip1_optimizer_engine(opt_ipvip1, atm, size(atm,2), bas, nBases, env)
call twoe_ip1ip2_optimizer_engine(opt_ip1ip2, atm, size(atm,2), bas, nBases, env)

tm_opt = omp_get_wtime() - tm0 ; tm0 = omp_get_wtime()
block
   real(8),allocatable :: d2E_local(:,:,:,:)
   !$omp parallel default(shared) &
   !$omp&   private(ish,jsh,ksh,lsh,di,dj,dk,dl,atI,atJ,atK,atL) &
   !$omp&   private(a,b,c,d,e1,e2,e3,e4,nq,coefv,pair_bound) &
   !$omp&   private(ij,kl,nperm,iperm,jperm,sh,ao,shtup,useperm,dup) &
   !$omp&   private(shlsA,shlsB,shlsC) &
   !$omp&   private(buf1,buf2,buf3,coefA,coefB,coefC) &
   !$omp&   private(sum_ii,sum_ij,sum_ik,sum_il,sum_jj,sum_ji,sum_jk,sum_jl) &
   !$omp&   private(sum_kk,sum_kl,sum_ki,sum_kj,sum_ll,sum_lk,sum_li,sum_lj) &
   !$omp&   private(d2E_local,tk0) &
   !$omp&   private(Mii,Mij,Mik,Mil,Mji,Mjj,Mjk,Mjl) &
   !$omp&   private(Mki,Mkj,Mkk,Mkl,Mli,Mlj,Mlk,Mll) &
   !$omp&   reduction(+:t_kern,t_coef,t_contr,nquad_done)
   allocate(buf1(maxdim**4*9), buf2(maxdim**4*9), buf3(maxdim**4*9))
   allocate(coefA(maxdim**4), coefB(maxdim**4), coefC(maxdim**4))
   allocate(d2E_local(3,natoms,3,natoms))
   d2E_local = 0.0d0
   !$omp do schedule(static,1)
   do ij = 0,npair-1
      ish = pair_i(ij)
      jsh = pair_j(ij)
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      atI = bas(1,ish+1)+1
      atJ = bas(1,jsh+1)+1
      di = shell_dim(ish); dj = shell_dim(jsh)
      do kl = 0,ij
         ksh = pair_i(kl)
         lsh = pair_j(kl)
               pair_bound = schwarz_bound(ish,jsh)*schwarz_bound(ksh,lsh)
               if (pair_bound .lt. SCHWARZ_CUTOFF) cycle
               if (pair_bound*max(dmax_shell(ish,jsh)*dmax_shell(ksh,lsh), &
                                  dmax_shell(ish,ksh)*dmax_shell(jsh,lsh), &
                                  dmax_shell(ish,lsh)*dmax_shell(jsh,ksh)) &
                   .lt. hess_2e_cutoff) cycle

               atK = bas(1,ksh+1)+1
               atL = bas(1,lsh+1)+1
               dk = shell_dim(ksh); dl = shell_dim(lsh)
               nq = di*dj*dk*dl

               tk0 = omp_get_wtime()
               sh = (/ish,jsh,ksh,lsh/)
               nperm = 0
               do iperm = 1,8
                  shtup = (/ sh(PERMMAP(1,iperm)), sh(PERMMAP(2,iperm)), &
                             sh(PERMMAP(3,iperm)), sh(PERMMAP(4,iperm)) /)
                  dup = .false.
                  do jperm = 1,nperm
                     if ( shtup(1) .eq. sh(useperm(1,jperm)) .and. &
                          shtup(2) .eq. sh(useperm(2,jperm)) .and. &
                          shtup(3) .eq. sh(useperm(3,jperm)) .and. &
                          shtup(4) .eq. sh(useperm(4,jperm)) ) then
                        dup = .true.
                        exit
                     endif
                  enddo
                  if (.not. dup) then
                     nperm = nperm + 1
                     useperm(:,nperm) = PERMMAP(:,iperm)
                  endif
               enddo

               do a = 1,di
                  e1 = ao_offset(ish)+a
                  do b = 1,dj
                     e2 = ao_offset(jsh)+b
                     do c = 1,dk
                        e3 = ao_offset(ksh)+c
                        do d = 1,dl
                           e4 = ao_offset(lsh)+d
                           ao = (/e1,e2,e3,e4/)
                           coefv = 0.0d0
                           do iperm = 1,nperm
                              coefv = coefv &
                                 + jscale*0.5d0*Ptot(ao(useperm(1,iperm)),ao(useperm(2,iperm))) &
                                        *Ptot(ao(useperm(3,iperm)),ao(useperm(4,iperm))) &
                                 - 0.5d0*HF_exchange_frac &
                                   *( Pa(ao(useperm(1,iperm)),ao(useperm(3,iperm))) &
                                     *Pa(ao(useperm(2,iperm)),ao(useperm(4,iperm))) &
                                    + Pb(ao(useperm(1,iperm)),ao(useperm(3,iperm))) &
                                     *Pb(ao(useperm(2,iperm)),ao(useperm(4,iperm))) )
                           enddo
                           coefv = coefv * NorVEC(e1)*NorVEC(e2)*NorVEC(e3)*NorVEC(e4)
                           coefA(a + (b-1)*di + (c-1)*di*dj + (d-1)*di*dj*dk) = coefv
                           coefB(b + (a-1)*dj + (c-1)*dj*di + (d-1)*dj*di*dk) = coefv
                           coefC(c + (d-1)*dk + (a-1)*dk*dl + (b-1)*dk*dl*di) = coefv
                        enddo
                     enddo
                  enddo
               enddo

               t_coef = t_coef + (omp_get_wtime()-tk0)

               shlsA = (/ish,jsh,ksh,lsh/)
               shlsB = (/jsh,ish,ksh,lsh/)
               shlsC = (/ksh,lsh,ish,jsh/)

               tk0 = omp_get_wtime()
               call twoe_ipip1_engine(buf1, shlsA, atm, size(atm,2), bas, nBases, env, opt_ipip1)
               call twoe_ipvip1_engine(buf2, shlsA, atm, size(atm,2), bas, nBases, env, opt_ipvip1)
               call twoe_ip1ip2_engine(buf3, shlsA, atm, size(atm,2), bas, nBases, env, opt_ip1ip2)
               t_kern = t_kern + (omp_get_wtime()-tk0)
               tk0 = omp_get_wtime()
               call hessian_contract_comp9(buf1, nq, coefA, sum_ii)
               call hessian_contract_comp9(buf2, nq, coefA, sum_ij)
               call hessian_contract_comp9(buf3, nq, coefA, sum_ik)
               t_contr = t_contr + (omp_get_wtime()-tk0)

               tk0 = omp_get_wtime()
               call twoe_ipip1_engine(buf1, shlsB, atm, size(atm,2), bas, nBases, env, opt_ipip1)
               call twoe_ip1ip2_engine(buf3, shlsB, atm, size(atm,2), bas, nBases, env, opt_ip1ip2)
               t_kern = t_kern + (omp_get_wtime()-tk0)
               tk0 = omp_get_wtime()
               call hessian_contract_comp9(buf1, nq, coefB, sum_jj)
               call hessian_contract_comp9(buf3, nq, coefB, sum_jk)
               t_contr = t_contr + (omp_get_wtime()-tk0)

               tk0 = omp_get_wtime()
               call twoe_ipip1_engine(buf1, shlsC, atm, size(atm,2), bas, nBases, env, opt_ipip1)
               t_kern = t_kern + (omp_get_wtime()-tk0)
               tk0 = omp_get_wtime()
               call hessian_contract_comp9(buf1, nq, coefC, sum_kk)
               t_contr = t_contr + (omp_get_wtime()-tk0)

               Mii = transpose(reshape(sum_ii,(/3,3/)))
               Mij = transpose(reshape(sum_ij,(/3,3/)))
               Mik = transpose(reshape(sum_ik,(/3,3/)))
               Mjj = transpose(reshape(sum_jj,(/3,3/)))
               Mjk = transpose(reshape(sum_jk,(/3,3/)))
               Mkk = transpose(reshape(sum_kk,(/3,3/)))

               Mji = transpose(Mij)
               Mki = transpose(Mik)
               Mkj = transpose(Mjk)
               Mli = -(Mii + Mji + Mki)
               Mlj = -(Mij + Mjj + Mkj)
               Mlk = -(Mik + Mjk + Mkk)
               Mil = transpose(Mli)
               Mjl = transpose(Mlj)
               Mkl = transpose(Mlk)
               Mll = -(Mil + Mjl + Mkl)

               d2E_local(:,atI,:,atI) = d2E_local(:,atI,:,atI) + Mii
               d2E_local(:,atI,:,atJ) = d2E_local(:,atI,:,atJ) + Mij
               d2E_local(:,atI,:,atK) = d2E_local(:,atI,:,atK) + Mik
               d2E_local(:,atI,:,atL) = d2E_local(:,atI,:,atL) + Mil

               d2E_local(:,atJ,:,atI) = d2E_local(:,atJ,:,atI) + Mji
               d2E_local(:,atJ,:,atJ) = d2E_local(:,atJ,:,atJ) + Mjj
               d2E_local(:,atJ,:,atK) = d2E_local(:,atJ,:,atK) + Mjk
               d2E_local(:,atJ,:,atL) = d2E_local(:,atJ,:,atL) + Mjl

               d2E_local(:,atK,:,atI) = d2E_local(:,atK,:,atI) + Mki
               d2E_local(:,atK,:,atJ) = d2E_local(:,atK,:,atJ) + Mkj
               d2E_local(:,atK,:,atK) = d2E_local(:,atK,:,atK) + Mkk
               d2E_local(:,atK,:,atL) = d2E_local(:,atK,:,atL) + Mkl

               d2E_local(:,atL,:,atI) = d2E_local(:,atL,:,atI) + Mli
               d2E_local(:,atL,:,atJ) = d2E_local(:,atL,:,atJ) + Mlj
               d2E_local(:,atL,:,atK) = d2E_local(:,atL,:,atK) + Mlk
               d2E_local(:,atL,:,atL) = d2E_local(:,atL,:,atL) + Mll
               nquad_done = nquad_done + 1
      enddo
   enddo
   !$omp end do
   block
      integer :: merge_tid, merge_nthreads
      merge_nthreads = omp_get_num_threads()
      do merge_tid = 0, merge_nthreads-1
         !$omp barrier
         if (omp_get_thread_num() .eq. merge_tid) then
            d2E_2e = d2E_2e + d2E_local
         endif
      enddo
   end block
   deallocate(buf1, buf2, buf3, coefA, coefB, coefC, d2E_local)
   !$omp end parallel
end block
tm_par = omp_get_wtime() - tm0 ; tm0 = omp_get_wtime()

block
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_HESSIAN_2E_TIMING", envval, status=envstat)
   if (envstat .eq. 0) then
      print *, '   2e-direct split (thread-summed s): kernels=', t_kern, ' coef=', t_coef, ' contract=', t_contr
      print *, '   2e-direct quadruples kept =', nquad_done, '  kernel calls =', nquad_done*6
   endif
end block

if (tm_on) print '(A,4F8.3)', ' HESS2E wall: prep,optimizers,parallel,rest:', &
   tm_pre, tm_opt, tm_par, omp_get_wtime() - tm0
call cintdel_optimizer(opt_ipip1)
call cintdel_optimizer(opt_ipvip1)
call cintdel_optimizer(opt_ip1ip2)
env(9) = 0.0d0
deallocate(ao_offset, shell_dim, dmax_shell, Ptot, pair_i, pair_j)
end subroutine integrals_hessian_2e_direct

subroutine hessian_contract_comp9(buf, nq, coef, s9)
implicit none
integer,intent(in) :: nq
real(8),intent(in) :: buf(nq,9), coef(nq)
real(8),intent(out) :: s9(9)
integer :: ic
do ic = 1,9
   s9(ic) = dot_product(buf(:,ic), coef)
enddo
end subroutine hessian_contract_comp9

module subroutine integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Eout)
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: Eout

integer :: ish,jsh,ksh,lsh,di,dj,dk,dl
integer :: a,b,c,d,e1,e2,e3,e4
integer :: shls(4)
real(8),allocatable :: buf_val(:,:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: si_off, num_off
real(8) :: coefv

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

Eout = 0.0d0

do ish = 0,nBases-1
   do jsh = 0,nBases-1
      do ksh = 0,nBases-1
         do lsh = 0,nBases-1
            di = shell_dim(ish); dj = shell_dim(jsh)
            dk = shell_dim(ksh); dl = shell_dim(lsh)
            allocate(buf_val(di,dj,dk,dl))
            shls = (/ish,jsh,ksh,lsh/)
            call twoe_engine(buf_val, shls, atm, size(atm,2), bas, nBases, env, 0_8)
            do a = 1,di
               e1 = ao_offset(ish)+a
               do b = 1,dj
                  e2 = ao_offset(jsh)+b
                  do c = 1,dk
                     e3 = ao_offset(ksh)+c
                     do d = 1,dl
                        e4 = ao_offset(lsh)+d
                        coefv = ( 0.5d0*(Pa(e1,e2)+Pb(e1,e2))*(Pa(e3,e4)+Pb(e3,e4)) &
                                - 0.5d0*HF_exchange_frac*(Pa(e1,e3)*Pa(e2,e4)+Pb(e1,e3)*Pb(e2,e4)) ) &
                                * NorVEC(e1)*NorVEC(e2)*NorVEC(e3)*NorVEC(e4)
                        Eout = Eout + buf_val(a,b,c,d)*coefv
                     enddo
                  enddo
               enddo
            enddo
            deallocate(buf_val)
         enddo
      enddo
   enddo
enddo

deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_e2e_fixed_density

module subroutine integrals_hessian_debug_check2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac)
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

integer :: atomA, atomB, c1, c2
real(8) :: h, Ep, Em, E0, Epp, Epm, Imp, Imm
real(8),allocatable :: d2E_2e(:,:,:,:)
real(8) :: d2_fd(3,3)
real(8) :: save3A(3), save3B(3)
integer :: envA, envB

h = 1.0d-3
atomA = 1
atomB = min(2,natoms)
envA = atm(2,atomA)+1
envB = atm(2,atomB)+1
save3A = env(envA:envA+2)
save3B = env(envB:envB+2)

allocate(d2E_2e(3,natoms,3,natoms))
call integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e)
call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, E0)
print *, 'HESSIAN DEBUG 2e-direct: E2e(fixed P) = ', E0

do c1 = 1,3
   do c2 = 1,3
      if (c1 == c2) then
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) + h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Ep)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) - h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Em)
         d2_fd(c1,c2) = (Ep - 2.0d0*E0 + Em)/(h*h)
      else
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) + h
         env(envA+c2-1) = save3A(c2) + h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Epp)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) + h
         env(envA+c2-1) = save3A(c2) - h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Epm)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) - h
         env(envA+c2-1) = save3A(c2) + h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Imp)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) - h
         env(envA+c2-1) = save3A(c2) - h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Imm)
         d2_fd(c1,c2) = (Epp - Epm - Imp + Imm)/(4.0d0*h*h)
      endif
   enddo
enddo
env(envA:envA+2) = save3A

print *, 'FD  d2E/dRA2   (row1) = ', d2_fd(1,:)
print *, 'FD  d2E/dRA2   (row2) = ', d2_fd(2,:)
print *, 'FD  d2E/dRA2   (row3) = ', d2_fd(3,:)
print *, 'analytic (row1)       = ', d2E_2e(1,atomA,:,atomA)
print *, 'analytic (row2)       = ', d2E_2e(2,atomA,:,atomA)
print *, 'analytic (row3)       = ', d2E_2e(3,atomA,:,atomA)

if (atomB /= atomA) then
   do c1 = 1,3
      do c2 = 1,3
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) + h
         env(envB+c2-1) = save3B(c2) + h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Epp)
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) + h
         env(envB+c2-1) = save3B(c2) - h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Epm)
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) - h
         env(envB+c2-1) = save3B(c2) + h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Imp)
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) - h
         env(envB+c2-1) = save3B(c2) - h
         call integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Imm)
         d2_fd(c1,c2) = (Epp - Epm - Imp + Imm)/(4.0d0*h*h)
      enddo
   enddo
   env(envA:envA+2) = save3A
   env(envB:envB+2) = save3B
   print *, 'FD  d2E/dRA.dRB (row1) = ', d2_fd(1,:)
   print *, 'FD  d2E/dRA.dRB (row2) = ', d2_fd(2,:)
   print *, 'FD  d2E/dRA.dRB (row3) = ', d2_fd(3,:)
   print *, 'analytic (row1)        = ', d2E_2e(1,atomA,:,atomB)
   print *, 'analytic (row2)        = ', d2E_2e(2,atomA,:,atomB)
   print *, 'analytic (row3)        = ', d2E_2e(3,atomA,:,atomB)
endif

deallocate(d2E_2e)
end subroutine integrals_hessian_debug_check2e_direct

module subroutine integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Eout)
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Ptot(nConts,nConts), Wtot(nConts,nConts)
real(8),intent(out) :: Eout

integer :: ish,jsh,di,dj,a,b,e1,e2
integer :: shls(4)
real(8),allocatable :: buf_kin(:,:), buf_nuc(:,:), buf_s(:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: si_off, num_off

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

Eout = 0.0d0
do ish = 0,nBases-1
   do jsh = 0,nBases-1
      shls(1) = ish
      shls(2) = jsh
      di = shell_dim(ish); dj = shell_dim(jsh)
      allocate(buf_kin(di,dj), buf_nuc(di,dj), buf_s(di,dj))
      call kin1e_engine(buf_kin, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      call nuc1e_engine(buf_nuc, shls, atm_nuc, nAtoms_nuc, bas, nBases, env_nuc, 0_8)
      call ovlp1e_engine(buf_s, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      do a = 1,di
         e1 = ao_offset(ish)+a
         do b = 1,dj
            e2 = ao_offset(jsh)+b
            Eout = Eout + ( Ptot(e1,e2)*(buf_kin(a,b)+buf_nuc(a,b)) - Wtot(e1,e2)*buf_s(a,b) ) &
                        * NorVEC(e1)*NorVEC(e2)
         enddo
      enddo
      deallocate(buf_kin, buf_nuc, buf_s)
   enddo
enddo
deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_e1e_fixed_density

module subroutine integrals_hessian_debug_check1e(nConts, natoms, Ptot, Wtot)
use MOL_info, only: atoms
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Ptot(nConts,nConts), Wtot(nConts,nConts)
INCLUDE 'parameter.h'

integer :: atomA, atomB, c1, c2, i, mu
real(8) :: h, Ep, Em, E0, Epp, Epm, Imp, Imm
real(8),allocatable :: d2E_1e(:,:,:,:)
real(8) :: d2_fd(3,3)
real(8) :: save3A(3), save3B(3)
integer :: envA, envB
integer,allocatable :: ao_atom(:), znum(:), ecpce(:)
real(8),allocatable :: coor_bohr(:,:)

allocate(ao_atom(nConts), znum(natoms), ecpce(natoms), coor_bohr(3,natoms))
mu = 0
do i = 1,natoms
   znum(i) = atoms(i)%charge
   ecpce(i) = atoms(i)%ecpCoreElec
   coor_bohr(:,i) = atoms(i)%coor*ans2bohr
   block
      integer :: j
      do j = 1,atoms(i)%nconts
         mu = mu+1
         ao_atom(mu) = i
      enddo
   end block
enddo

allocate(d2E_1e(3,natoms,3,natoms))
call calc_hessian_1e(nConts, natoms, ao_atom, Ptot, Wtot, znum, ecpce, coor_bohr, d2E_1e)

h = 1.0d-3
atomA = 1
atomB = min(2,natoms)
envA = atm(2,atomA)+1
envB = atm(2,atomB)+1
save3A = env(envA:envA+2)
save3B = env(envB:envB+2)

call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, E0)
print *, 'HESSIAN DEBUG 1e: E1e(fixed P,W) = ', E0

do c1 = 1,3
   do c2 = 1,3
      if (c1 == c2) then
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) + h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Ep)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) - h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Em)
         d2_fd(c1,c2) = (Ep - 2.0d0*E0 + Em)/(h*h)
      else
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) + h
         env(envA+c2-1) = save3A(c2) + h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Epp)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) + h
         env(envA+c2-1) = save3A(c2) - h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Epm)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) - h
         env(envA+c2-1) = save3A(c2) + h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Imp)
         env(envA:envA+2) = save3A
         env(envA+c1-1) = save3A(c1) - h
         env(envA+c2-1) = save3A(c2) - h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Imm)
         d2_fd(c1,c2) = (Epp - Epm - Imp + Imm)/(4.0d0*h*h)
      endif
   enddo
enddo
env(envA:envA+2) = save3A

print *, 'FD  d2E/dRA2   (row1) = ', d2_fd(1,:)
print *, 'FD  d2E/dRA2   (row2) = ', d2_fd(2,:)
print *, 'FD  d2E/dRA2   (row3) = ', d2_fd(3,:)
print *, 'analytic (row1)       = ', d2E_1e(1,atomA,:,atomA)
print *, 'analytic (row2)       = ', d2E_1e(2,atomA,:,atomA)
print *, 'analytic (row3)       = ', d2E_1e(3,atomA,:,atomA)

if (atomB /= atomA) then
   do c1 = 1,3
      do c2 = 1,3
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) + h
         env(envB+c2-1) = save3B(c2) + h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Epp)
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) + h
         env(envB+c2-1) = save3B(c2) - h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Epm)
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) - h
         env(envB+c2-1) = save3B(c2) + h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Imp)
         env(envA:envA+2) = save3A
         env(envB:envB+2) = save3B
         env(envA+c1-1) = save3A(c1) - h
         env(envB+c2-1) = save3B(c2) - h
         call integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Imm)
         d2_fd(c1,c2) = (Epp - Epm - Imp + Imm)/(4.0d0*h*h)
      enddo
   enddo
   env(envA:envA+2) = save3A
   env(envB:envB+2) = save3B
   print *, 'FD  d2E/dRA.dRB (row1) = ', d2_fd(1,:)
   print *, 'FD  d2E/dRA.dRB (row2) = ', d2_fd(2,:)
   print *, 'FD  d2E/dRA.dRB (row3) = ', d2_fd(3,:)
   print *, 'analytic (row1)        = ', d2E_1e(1,atomA,:,atomB)
   print *, 'analytic (row2)        = ', d2E_1e(2,atomA,:,atomB)
   print *, 'analytic (row3)        = ', d2E_1e(3,atomA,:,atomB)
endif

deallocate(d2E_1e, ao_atom, znum, ecpce, coor_bohr)
end subroutine integrals_hessian_debug_check1e

module subroutine integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2, dK2, dK2b, omega)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, iatom
real(8),intent(in),optional :: omega
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: dJ2(nConts,nConts,3), dK2(nConts,nConts,3)
real(8),intent(out),optional :: dK2b(nConts,nConts,3)

integer :: ish,jsh,ksh,lsh,di,dj,dk,dl
integer :: shls2e(4)
real(8),allocatable :: buf(:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer(8) :: opt
integer :: a,b,c,d,e1,e2,e3,e4
real(8) :: intVal(3)
real(8),allocatable :: dmax_shell(:,:)
real(8) :: density_role_cutoff
integer :: iw, nwork, maxdim
integer,allocatable :: work_i(:), work_j(:)

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
block
   integer :: si_off, num_off
   num_off = 0
   do si_off = 0,nBases-1
      ao_offset(si_off) = num_off
      shell_dim(si_off) = cgto_engine(si_off, bas)
      num_off = num_off + shell_dim(si_off)
   enddo
end block

if (present(omega)) then
   if (omega > 0.0d0) env(9) = omega
endif
call twoe_ip1_optimizer_engine(opt, atm, size(atm,2), bas, nBases, env)

allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts, abs(Pa)+abs(Pb), dmax_shell)
density_role_cutoff = 1.0d-9
maxdim = maxval(shell_dim)
allocate(work_i(nBases*nBases), work_j(nBases*nBases))

nwork = 0
do ish = 0,nBases-1
   if (bas(1,ish+1)+1 .ne. iatom) cycle
   do jsh = 0,nBases-1
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      nwork = nwork + 1
      work_i(nwork) = ish
      work_j(nwork) = jsh
   enddo
enddo

dJ2 = 0.0d0
dK2 = 0.0d0
if (present(dK2b)) dK2b = 0.0d0
   block
      real(8),allocatable :: dJ2_local(:,:,:), dK2_local(:,:,:), dK2b_local(:,:,:)
      !$omp parallel default(shared) &
      !$omp&   private(iw,ish,di,jsh,dj,ksh,dk,lsh,dl,shls2e,a,b,c,d,e1,e2,e3,e4,intVal) &
      !$omp&   private(buf,dJ2_local,dK2_local,dK2b_local)
      allocate(dJ2_local(nConts,nConts,3), dK2_local(nConts,nConts,3))
      allocate(dK2b_local(nConts,nConts,3))
      allocate(buf(maxdim**4*3))
      dJ2_local = 0.0d0
      dK2_local = 0.0d0
      dK2b_local = 0.0d0
      !$omp do schedule(static,1)
      do iw = 1,nwork
         ish = work_i(iw)
         jsh = work_j(iw)
         di = shell_dim(ish)
         dj = shell_dim(jsh)
         do ksh = 0,nBases-1
            dk = shell_dim(ksh)
            do lsh = 0,ksh
               dl = shell_dim(lsh)
               if (schwarz_bound(ish,jsh)*schwarz_bound(ksh,lsh) .lt. SCHWARZ_CUTOFF) cycle
               if (schwarz_bound(ish,jsh)*schwarz_bound(ksh,lsh) &
                   *max(dmax_shell(ish,jsh), dmax_shell(lsh,ish), dmax_shell(ksh,ish)) &
                   .lt. density_role_cutoff) cycle
               shls2e = (/ish,jsh,ksh,lsh/)
                call cphf_density_role_quad(buf, shls2e, di, dj, dk, dl, opt, &
                                             nConts, ao_offset(ish), ao_offset(jsh), &
                                             ao_offset(ksh), ao_offset(lsh), &
                                             Pa, Pb, dJ2_local, dK2_local, dK2b_local, ksh .ne. lsh)
            enddo
         enddo
      enddo
      !$omp end do
      block
         integer :: merge_tid, merge_nthreads
         merge_nthreads = omp_get_num_threads()
         do merge_tid = 0, merge_nthreads-1
            !$omp barrier
            if (omp_get_thread_num() .eq. merge_tid) then
               dJ2 = dJ2 + dJ2_local
               dK2 = dK2 + dK2_local
               if (present(dK2b)) dK2b = dK2b + dK2b_local
            endif
         enddo
      end block
      deallocate(dJ2_local, dK2_local, dK2b_local, buf)
      !$omp end parallel
   end block

call cintdel_optimizer(opt)
env(9) = 0.0d0
deallocate(ao_offset, shell_dim, dmax_shell, work_i, work_j)
end subroutine integrals_cphf_2e_density_role

subroutine cphf_density_role_quad(buf, shls2e, di, dj, dk, dl, opt, nConts, &
                                   off_i, off_j, off_k, off_l, Pa, Pb, dJ2_local, dK2_local, dK2b_local, add_swap)
implicit none
integer,intent(in) :: di,dj,dk,dl,nConts,off_i,off_j,off_k,off_l
integer,intent(in) :: shls2e(4)
integer(8),intent(in) :: opt
real(8),intent(out) :: buf(di,dj,dk,dl,3)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(inout) :: dJ2_local(nConts,nConts,3), dK2_local(nConts,nConts,3)
real(8),intent(inout) :: dK2b_local(nConts,nConts,3)
logical,intent(in) :: add_swap
integer :: a,b,c,d,e1,e2,e3,e4
real(8) :: intVal(3)

call twoe_ip1_engine(buf, shls2e, atm, size(atm,2), bas, nBases, env, opt)
do a = 1,di
   e1 = off_i+a
   do b = 1,dj
      e2 = off_j+b
      do c = 1,dk
         e3 = off_k+c
         do d = 1,dl
            e4 = off_l+d
            intVal = buf(a,b,c,d,:)*NorVEC(e1)*NorVEC(e2)*NorVEC(e3)*NorVEC(e4)
            dJ2_local(e3,e4,:) = dJ2_local(e3,e4,:) + intVal*(Pa(e1,e2)+Pb(e1,e2))
            dK2_local(e3,e2,:) = dK2_local(e3,e2,:) + intVal*Pa(e4,e1)
            dK2b_local(e3,e2,:) = dK2b_local(e3,e2,:) + intVal*Pb(e4,e1)
            if (add_swap) then
               dJ2_local(e4,e3,:) = dJ2_local(e4,e3,:) + intVal*(Pa(e1,e2)+Pb(e1,e2))
               dK2_local(e4,e2,:) = dK2_local(e4,e2,:) + intVal*Pa(e3,e1)
               dK2b_local(e4,e2,:) = dK2b_local(e4,e2,:) + intVal*Pb(e3,e1)
            endif
         enddo
      enddo
   enddo
enddo
end subroutine cphf_density_role_quad

end submodule hessian_impl
