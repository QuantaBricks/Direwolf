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

allocate(AA(nConts,nConts,9), AB(nConts,nConts,9), KK(nConts,nConts,9))
do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i
      shls(2) = j
      di = shell_dim(i)
      dj = shell_dim(j)
      allocate(buf_aa(di,dj,9), buf_ab(di,dj,9))
      call rinv1e_ipiprinv_engine(buf_aa, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      call rinv1e_iprinvip_engine(buf_ab, shls, atm, size(atm,2), bas, nBases, env, 0_8)
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
deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_nuc_atom

module subroutine integrals_hessian_nuc_atom_raw(coor_bohr, nConts, AA, AB)
implicit none
real(8),intent(in) :: coor_bohr(3)
integer,intent(in) :: nConts
real(8),intent(out) :: AA(nConts,nConts,9), AB(nConts,nConts,9)

integer :: i,j,di,dj,e1,e2,a,b
integer :: shls(2)
real(8),allocatable :: buf_aa(:,:,:), buf_ab(:,:,:)
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

do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i
      shls(2) = j
      di = shell_dim(i)
      dj = shell_dim(j)
      allocate(buf_aa(di,dj,9), buf_ab(di,dj,9))
      call rinv1e_ipiprinv_engine(buf_aa, shls, atm, size(atm,2), bas, nBases, env, 0_8)
      call rinv1e_iprinvip_engine(buf_ab, shls, atm, size(atm,2), bas, nBases, env, 0_8)
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

deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_nuc_atom_raw

module subroutine integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e)
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: d2E_2e(3,natoms,3,natoms)

integer :: ish,jsh,ksh,lsh,di,dj,dk,dl
integer :: atI,atJ,atK,atL
integer :: a,b,c,d,e1,e2,e3,e4
integer :: shlsA(4), shlsB(4), shlsC(4), shlsD(4)
real(8),allocatable :: bA_ipip1(:,:,:,:,:), bA_ipvip1(:,:,:,:,:), bA_ip1ip2(:,:,:,:,:)
real(8),allocatable :: bB_ipip1(:,:,:,:,:), bB_ipvip1(:,:,:,:,:), bB_ip1ip2(:,:,:,:,:)
real(8),allocatable :: bC_ipip1(:,:,:,:,:), bC_ipvip1(:,:,:,:,:), bC_ip1ip2(:,:,:,:,:)
real(8),allocatable :: bD_ipip1(:,:,:,:,:), bD_ipvip1(:,:,:,:,:), bD_ip1ip2(:,:,:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: si_off, num_off
real(8) :: coefv
real(8) :: sum_ii(9), sum_ij(9), sum_ik(9), sum_il(9)
real(8) :: sum_jj(9), sum_ji(9), sum_jk(9), sum_jl(9)
real(8) :: sum_kk(9), sum_kl(9), sum_ki(9), sum_kj(9)
real(8) :: sum_ll(9), sum_lk(9), sum_li(9), sum_lj(9)

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

d2E_2e = 0.0d0

do ish = 0,nBases-1
   atI = bas(1,ish+1)+1
   do jsh = 0,nBases-1
      atJ = bas(1,jsh+1)+1
      do ksh = 0,nBases-1
         atK = bas(1,ksh+1)+1
         do lsh = 0,nBases-1
            atL = bas(1,lsh+1)+1

            di = shell_dim(ish); dj = shell_dim(jsh)
            dk = shell_dim(ksh); dl = shell_dim(lsh)

            allocate(bA_ipip1(di,dj,dk,dl,9), bA_ipvip1(di,dj,dk,dl,9), bA_ip1ip2(di,dj,dk,dl,9))
            allocate(bB_ipip1(dj,di,dk,dl,9), bB_ipvip1(dj,di,dk,dl,9), bB_ip1ip2(dj,di,dk,dl,9))
            allocate(bC_ipip1(dk,dl,di,dj,9), bC_ipvip1(dk,dl,di,dj,9), bC_ip1ip2(dk,dl,di,dj,9))
            allocate(bD_ipip1(dl,dk,di,dj,9), bD_ipvip1(dl,dk,di,dj,9), bD_ip1ip2(dl,dk,di,dj,9))

            shlsA = (/ish,jsh,ksh,lsh/)
            shlsB = (/jsh,ish,ksh,lsh/)
            shlsC = (/ksh,lsh,ish,jsh/)
            shlsD = (/lsh,ksh,ish,jsh/)

            call twoe_ipip1_engine(bA_ipip1, shlsA, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ipvip1_engine(bA_ipvip1, shlsA, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ip1ip2_engine(bA_ip1ip2, shlsA, atm, size(atm,2), bas, nBases, env, 0_8)

            call twoe_ipip1_engine(bB_ipip1, shlsB, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ipvip1_engine(bB_ipvip1, shlsB, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ip1ip2_engine(bB_ip1ip2, shlsB, atm, size(atm,2), bas, nBases, env, 0_8)

            call twoe_ipip1_engine(bC_ipip1, shlsC, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ipvip1_engine(bC_ipvip1, shlsC, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ip1ip2_engine(bC_ip1ip2, shlsC, atm, size(atm,2), bas, nBases, env, 0_8)

            call twoe_ipip1_engine(bD_ipip1, shlsD, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ipvip1_engine(bD_ipvip1, shlsD, atm, size(atm,2), bas, nBases, env, 0_8)
            call twoe_ip1ip2_engine(bD_ip1ip2, shlsD, atm, size(atm,2), bas, nBases, env, 0_8)

            sum_ii=0d0; sum_ij=0d0; sum_ik=0d0; sum_il=0d0
            sum_jj=0d0; sum_ji=0d0; sum_jk=0d0; sum_jl=0d0
            sum_kk=0d0; sum_kl=0d0; sum_ki=0d0; sum_kj=0d0
            sum_ll=0d0; sum_lk=0d0; sum_li=0d0; sum_lj=0d0

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

                        sum_ii = sum_ii + bA_ipip1(a,b,c,d,:) *coefv
                        sum_ij = sum_ij + bA_ipvip1(a,b,c,d,:)*coefv
                        sum_ik = sum_ik + bA_ip1ip2(a,b,c,d,:)*coefv

                        sum_jj = sum_jj + bB_ipip1(b,a,c,d,:) *coefv
                        sum_ji = sum_ji + bB_ipvip1(b,a,c,d,:)*coefv
                        sum_jk = sum_jk + bB_ip1ip2(b,a,c,d,:)*coefv

                        sum_kk = sum_kk + bC_ipip1(c,d,a,b,:) *coefv
                        sum_kl = sum_kl + bC_ipvip1(c,d,a,b,:)*coefv
                        sum_ki = sum_ki + bC_ip1ip2(c,d,a,b,:)*coefv

                        sum_ll = sum_ll + bD_ipip1(d,c,a,b,:) *coefv
                        sum_lk = sum_lk + bD_ipvip1(d,c,a,b,:)*coefv
                        sum_li = sum_li + bD_ip1ip2(d,c,a,b,:)*coefv
                     enddo
                  enddo
               enddo
            enddo

            sum_il = -(sum_ii + sum_ij + sum_ik)
            sum_jl = -(sum_jj + sum_ji + sum_jk)
            sum_kj = -(sum_kk + sum_kl + sum_ki)
            sum_lj = -(sum_ll + sum_lk + sum_li)

            d2E_2e(:,atI,:,atI) = d2E_2e(:,atI,:,atI) + transpose(reshape(sum_ii,(/3,3/)))
            d2E_2e(:,atI,:,atJ) = d2E_2e(:,atI,:,atJ) + transpose(reshape(sum_ij,(/3,3/)))
            d2E_2e(:,atI,:,atK) = d2E_2e(:,atI,:,atK) + transpose(reshape(sum_ik,(/3,3/)))
            d2E_2e(:,atI,:,atL) = d2E_2e(:,atI,:,atL) + transpose(reshape(sum_il,(/3,3/)))

            d2E_2e(:,atJ,:,atJ) = d2E_2e(:,atJ,:,atJ) + transpose(reshape(sum_jj,(/3,3/)))
            d2E_2e(:,atJ,:,atI) = d2E_2e(:,atJ,:,atI) + transpose(reshape(sum_ji,(/3,3/)))
            d2E_2e(:,atJ,:,atK) = d2E_2e(:,atJ,:,atK) + transpose(reshape(sum_jk,(/3,3/)))
            d2E_2e(:,atJ,:,atL) = d2E_2e(:,atJ,:,atL) + transpose(reshape(sum_jl,(/3,3/)))

            d2E_2e(:,atK,:,atK) = d2E_2e(:,atK,:,atK) + transpose(reshape(sum_kk,(/3,3/)))
            d2E_2e(:,atK,:,atL) = d2E_2e(:,atK,:,atL) + transpose(reshape(sum_kl,(/3,3/)))
            d2E_2e(:,atK,:,atI) = d2E_2e(:,atK,:,atI) + transpose(reshape(sum_ki,(/3,3/)))
            d2E_2e(:,atK,:,atJ) = d2E_2e(:,atK,:,atJ) + transpose(reshape(sum_kj,(/3,3/)))

            d2E_2e(:,atL,:,atL) = d2E_2e(:,atL,:,atL) + transpose(reshape(sum_ll,(/3,3/)))
            d2E_2e(:,atL,:,atK) = d2E_2e(:,atL,:,atK) + transpose(reshape(sum_lk,(/3,3/)))
            d2E_2e(:,atL,:,atI) = d2E_2e(:,atL,:,atI) + transpose(reshape(sum_li,(/3,3/)))
            d2E_2e(:,atL,:,atJ) = d2E_2e(:,atL,:,atJ) + transpose(reshape(sum_lj,(/3,3/)))

            deallocate(bA_ipip1, bA_ipvip1, bA_ip1ip2)
            deallocate(bB_ipip1, bB_ipvip1, bB_ip1ip2)
            deallocate(bC_ipip1, bC_ipvip1, bC_ip1ip2)
            deallocate(bD_ipip1, bD_ipvip1, bD_ip1ip2)
         enddo
      enddo
   enddo
enddo

deallocate(ao_offset, shell_dim)
end subroutine integrals_hessian_2e_direct

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

module subroutine integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2, dK2)
implicit none
integer,intent(in) :: nConts, iatom
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: dJ2(nConts,nConts,3), dK2(nConts,nConts,3)

integer :: ish,jsh,ksh,lsh,di,dj,dk,dl
integer :: shls2e(4)
real(8),allocatable :: buf(:,:,:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer(8) :: opt
integer :: a,b,c,d,e1,e2,e3,e4
real(8) :: intVal(3)

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

call twoe_ip1_optimizer_engine(opt, atm, size(atm,2), bas, nBases, env)

dJ2 = 0.0d0
dK2 = 0.0d0
do ish = 0,nBases-1
   if (bas(1,ish+1)+1 .ne. iatom) cycle
   di = shell_dim(ish)
   do jsh = 0,nBases-1
      dj = shell_dim(jsh)
      do ksh = 0,nBases-1
         dk = shell_dim(ksh)
         do lsh = 0,nBases-1
            dl = shell_dim(lsh)
            shls2e = (/ish,jsh,ksh,lsh/)
            allocate(buf(di,dj,dk,dl,3))
            call twoe_ip1_engine(buf, shls2e, atm, size(atm,2), bas, nBases, env, opt)
            do a = 1,di
               e1 = ao_offset(ish)+a
               do b = 1,dj
                  e2 = ao_offset(jsh)+b
                  do c = 1,dk
                     e3 = ao_offset(ksh)+c
                     do d = 1,dl
                        e4 = ao_offset(lsh)+d
                        intVal = buf(a,b,c,d,:)*NorVEC(e1)*NorVEC(e2)*NorVEC(e3)*NorVEC(e4)
                        dJ2(e3,e4,:) = dJ2(e3,e4,:) + intVal*(Pa(e1,e2)+Pb(e1,e2))
                        dK2(e3,e2,:) = dK2(e3,e2,:) + intVal*Pa(e4,e1)
                     enddo
                  enddo
               enddo
            enddo
            deallocate(buf)
         enddo
      enddo
   enddo
enddo

call cintdel_optimizer(opt)
deallocate(ao_offset, shell_dim)
end subroutine integrals_cphf_2e_density_role

end submodule hessian_impl
