! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: the closed-shell orbital Hessian applied to vir-occ vectors,

submodule (mod_cpscf) cpscf_apb_impl
implicit none
contains

module subroutine cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac, scf_grid)
use mod_integrals, only: engine_use_df_j
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
logical,intent(in),optional :: scf_grid
logical :: on_scf
on_scf = .false.
if (present(scf_grid)) on_scf = scf_grid
cpscf_xc_active = (HF_exchange_frac < 1.0d0)
if (cpscf_xc_active) then
   if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
   allocate(cpscf_xc_Ptot(nConts,nConts))
   cpscf_xc_Ptot = Pa + Pb
endif
cpscf_force_dfj = engine_use_df_j
block
   use mod_integrals, only: integrals_cosx_resp_cache_clear
   call integrals_cosx_resp_cache_clear()
end block
if (cpscf_xc_active .and. .not. on_scf) then
   block
      use mod_xc_hess, only: xc_resp_grid_build, xc_resp_grid_build_td
      use GRID_info, only: XCGRID_TIER_LO
      character(8) :: e
      call get_environment_variable("ENGINE_CPSCF_XC_GRID", e)
      if (trim(e) == "lo") then
         call xc_resp_grid_build(XCGRID_TIER_LO)
      else if (trim(e) .ne. "full") then
         call xc_resp_grid_build_td()
      endif
   end block
endif
cphf_voo_track = .false.
cphf_omega2 = 0.0d0
cphf_t_vj = 0.0d0 ; cphf_t_vk = 0.0d0 ; cphf_t_vxc = 0.0d0
cphf_n_vk = 0 ; cphf_nd_vk = 0
end subroutine cpscf_response_begin

module subroutine cpscf_response_end()
use mod_xc_hess, only: xc_resp_grid_clear
use mod_integrals, only: integrals_cosx_resp_cache_clear
implicit none
call xc_resp_grid_clear()
call integrals_cosx_resp_cache_clear()
cpscf_force_dfj = .false.
cpscf_xc_active = .false.
if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
cpscf_xc_uks = .false.
if (allocated(cpscf_xc_Pa)) deallocate(cpscf_xc_Pa, cpscf_xc_Pb)
end subroutine cpscf_response_end

module subroutine cpscf_apb_rhf(nConts, nocc, nvir, nvec, C_a, eLev_a, HF_exchange_frac, X, Q)
use mod_xc_hess, only: xc_fxc_project_vo
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc, nvir, nvec
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: X(nvir,nocc,nvec)
real(8),intent(out) :: Q(nvir,nocc,nvec)
real(8),allocatable :: Dtr(:,:,:), Jtr(:,:,:), Katr(:,:,:), Xp(:,:,:), Tmp(:,:), T2(:,:), CU(:,:)
real(8),allocatable :: VC(:,:,:), LU(:,:,:)
logical :: vcm
integer :: c, i, a
real(8) :: tv0

if (nvec < 1) return
vcm = cpscf_use_vind_c(HF_exchange_frac)
if (vcm) then
   allocate(Dtr(1,1,1), Jtr(1,1,1), Katr(1,1,1), LU(nConts,nocc,nvec))
else
   allocate(Dtr(nConts,nConts,nvec), Jtr(nConts,nConts,nvec), Katr(nConts,nConts,nvec))
endif
allocate(Tmp(nConts,nocc), T2(nConts,nConts), CU(nConts,nocc))
do c = 1,nvec
   call dgemm('N','N',nConts,nocc,nvir,1.0d0,C_a(1,nocc+1),nConts,X(1,1,c),nvir,0.0d0,CU,nConts)
   if (vcm) then
      LU(:,:,c) = 2.0d0*CU
   else
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,CU,nConts,C_a,nConts,0.0d0,T2,nConts)
      Dtr(:,:,c) = 2.0d0*(T2 + transpose(T2))
   endif
enddo
if (cpscf_triplet .and. HF_exchange_frac /= 0.0d0) &
   error stop 'cpscf_apb_rhf: triplet with exact exchange goes through cpscf_apb_amb_rhf'
if (cpscf_triplet) then
   if (vcm) then
      allocate(VC(nConts,nocc,nvec))
      VC = 0.0d0
   else
      Jtr = 0.0d0
   endif
else if (vcm) then
   allocate(VC(nConts,nocc,nvec))
   call cpscf_vind_c(nConts, nvec, HF_exchange_frac, nocc, C_a, LU, VC)
else
   call cpscf_vind(nConts, nvec, HF_exchange_frac, Dtr, Jtr, Katr, skip_xc=.true.)
endif
if (cpscf_xc_active) then
   allocate(Xp(nvir,nocc,nvec))
   tv0 = omp_get_wtime()
   call apb_xc_vo(nConts, nocc, nvir, nvec, C_a, X, Xp)
   cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
endif
do c = 1,nvec
   if (vcm) then
      Tmp = VC(:,:,c)
   else
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Jtr(1,1,c),nConts,C_a,nConts,0.0d0,Tmp,nConts)
   endif
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,Tmp,nConts,0.0d0,Q(1,1,c),nvir)
   if (cpscf_xc_active) Q(:,:,c) = Q(:,:,c) + Xp(:,:,c)
   do i = 1,nocc
      do a = 1,nvir
         Q(a,i,c) = Q(a,i,c) + (eLev_a(nocc+a) - eLev_a(i))*X(a,i,c)
      enddo
   enddo
enddo
deallocate(Dtr, Jtr, Katr, Tmp, T2, CU)
if (allocated(Xp)) deallocate(Xp)
if (allocated(VC)) deallocate(VC, LU)
end subroutine cpscf_apb_rhf

module subroutine cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, np, Xp, Qp, nm, Xm, Qm, Qoo)
use mod_xc_hess, only: xc_fxc_project_vo
use mod_integrals, only: integrals_build_exchange_cosx_batch_c2g, integrals_build_exchange_cosx_batch_raw, &
                         integrals_build_coulomb_batch, &
                         integrals_build_coulomb_df_batch, integrals_build_exchange_batch, direct_mode
use mod_exchange, only: cosx_enabled, RS_omega, RS_beta
use mod_integrals, only: integrals_build_exchange_lr_batch
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc, nvir, np, nm
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Xp(nvir,nocc,np), Xm(nvir,nocc,nm)
real(8),intent(out) :: Qp(nvir,nocc,np), Qm(nvir,nocc,nm)
real(8),intent(out),optional :: Qoo(nocc,nocc,np)
real(8),allocatable :: L1(:,:,:), L2(:,:,:), Ls(:,:,:), KC(:,:,:), KC2(:,:,:), VC(:,:,:), Xk(:,:,:), &
                       Gd(:,:,:), Dp(:,:,:), Jb(:,:,:), Kb(:,:,:), Kd(:,:,:), T(:,:), Sv(:,:), Dv(:,:)
integer :: c, i, a, npr
real(8) :: tv0, hx
logical :: cosx_ao, dmode0

npr = max(np, nm)
if (npr < 1) return
block
   character(8) :: eao
   call get_environment_variable('ENGINE_TDDFT_COSX_AO', eao)
   cosx_ao = (trim(eao) == '1')
end block
hx = 0.5d0*HF_exchange_frac
if (HF_exchange_frac == 0.0d0 .and. RS_omega > 0.0d0) hx = 0.5d0*RS_beta
allocate(T(nConts,nocc), Sv(nvir,nocc), Dv(nvir,nocc))
allocate(L1(nConts,nocc,npr), L2(nConts,nocc,npr), Ls(nConts,nocc,max(np,1)))
do c = 1,npr
   Sv = 0.0d0 ; Dv = 0.0d0
   if (c <= np) Sv = Xp(:,:,c)
   if (c <= nm) Dv = Xm(:,:,c)
   call dgemm('N','N',nConts,nocc,nvir,2.0d0,C_a(1,nocc+1),nConts,Sv+Dv,nvir,0.0d0,L1(1,1,c),nConts)
   call dgemm('N','N',nConts,nocc,nvir,2.0d0,C_a(1,nocc+1),nConts,Sv-Dv,nvir,0.0d0,L2(1,1,c),nConts)
   if (c <= np) call dgemm('N','N',nConts,nocc,nvir,2.0d0,C_a(1,nocc+1),nConts,Sv,nvir,0.0d0,Ls(1,1,c),nConts)
enddo
allocate(KC(nConts,nocc,npr), KC2(nConts,nocc,npr))
if (cosx_enabled .and. .not. cosx_ao) then
   tv0 = omp_get_wtime()
   call integrals_build_exchange_cosx_batch_c2g(nConts, npr, L1, L2, nocc, C_a, KC, KC2)
   cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   do c = 1,npr
      T = 0.5d0*(KC(:,:,c) + KC2(:,:,c))
      KC2(:,:,c) = 0.5d0*(KC2(:,:,c) - KC(:,:,c))
      KC(:,:,c) = T
   enddo
else if (cosx_enabled) then
   allocate(Gd(nConts,nConts,npr), Kb(nConts,nConts,npr), Kd(nConts,nConts,1))
   do c = 1,npr
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,L1(1,1,c),nConts,C_a,nConts,0.0d0,Gd(1,1,c),nConts)
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,C_a,nConts,L2(1,1,c),nConts,1.0d0,Gd(1,1,c),nConts)
   enddo
   tv0 = omp_get_wtime()
   call integrals_build_exchange_cosx_batch_raw(nConts, npr, Gd, Kb)
   cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   do c = 1,npr
      Kd(:,:,1) = 0.5d0*(Kb(:,:,c) + transpose(Kb(:,:,c)))
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Kd(1,1,1),nConts,C_a,nConts,0.0d0,KC(1,1,c),nConts)
      Kd(:,:,1) = 0.5d0*(transpose(Kb(:,:,c)) - Kb(:,:,c))
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Kd(1,1,1),nConts,C_a,nConts,0.0d0,KC2(1,1,c),nConts)
   enddo
   deallocate(Gd, Kb, Kd)
else
   dmode0 = direct_mode
   direct_mode = .true.
   allocate(Gd(nConts,nConts,npr), Kb(nConts,nConts,npr), Kd(nConts,nConts,npr))
   do c = 1,npr
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,L1(1,1,c),nConts,C_a,nConts,0.0d0,Gd(1,1,c),nConts)
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,C_a,nConts,L2(1,1,c),nConts,1.0d0,Gd(1,1,c),nConts)
   enddo
   tv0 = omp_get_wtime()
   if (HF_exchange_frac /= 0.0d0) then
      call integrals_build_exchange_batch(nConts, npr, Gd, Gd, Kb, Kd)
      Kb = HF_exchange_frac*Kb
   else
      Kb = 0.0d0
   endif
   if (RS_omega > 0.0d0) then
      block
         real(8),allocatable :: Kl(:,:,:)
         allocate(Kl(nConts,nConts,npr))
         call integrals_build_exchange_lr_batch(nConts, npr, Gd, Gd, Kl, Kd)
         Kb = Kb + RS_beta*Kl
         deallocate(Kl)
      end block
   endif
   direct_mode = dmode0
   hx = 0.5d0
   cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   do c = 1,npr
      Kd(:,:,1) = 0.5d0*(Kb(:,:,c) + transpose(Kb(:,:,c)))
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Kd(1,1,1),nConts,C_a,nConts,0.0d0,KC(1,1,c),nConts)
      Kd(:,:,1) = 0.5d0*(Kb(:,:,c) - transpose(Kb(:,:,c)))
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Kd(1,1,1),nConts,C_a,nConts,0.0d0,KC2(1,1,c),nConts)
   enddo
   deallocate(Gd, Kb, Kd)
endif
if (np > 0) then
   allocate(VC(nConts,nocc,np))
   if (cpscf_triplet) then
      VC = 0.0d0
   else
      call cpscf_vind_c(nConts, np, 0.0d0, nocc, C_a, Ls(:,:,1:np), VC)
   endif
   do c = 1,np
      VC(:,:,c) = VC(:,:,c) - hx*KC(:,:,c)
      call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,VC(1,1,c),nConts,0.0d0,Qp(1,1,c),nvir)
      if (present(Qoo)) call dgemm('T','N',nocc,nocc,nConts,1.0d0,C_a,nConts,VC(1,1,c),nConts,0.0d0,Qoo(1,1,c),nocc)
   enddo
   deallocate(VC)
endif
do c = 1,nm
   call dgemm('T','N',nvir,nocc,nConts,-hx,C_a(1,nocc+1),nConts,KC2(1,1,c),nConts,0.0d0,Qm(1,1,c),nvir)
enddo
if (cpscf_xc_active .and. np > 0) then
   allocate(Xk(nvir,nocc,np))
   tv0 = omp_get_wtime()
   call apb_xc_vo(nConts, nocc, nvir, np, C_a, Xp, Xk)
   cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
   Qp = Qp + Xk
   deallocate(Xk)
endif
do i = 1,nocc
   do a = 1,nvir
      Qp(a,i,:) = Qp(a,i,:) + (eLev_a(nocc+a) - eLev_a(i))*Xp(a,i,:)
      Qm(a,i,:) = Qm(a,i,:) + (eLev_a(nocc+a) - eLev_a(i))*Xm(a,i,:)
   enddo
enddo
deallocate(T, Sv, Dv, L1, L2, Ls, KC, KC2)
end subroutine cpscf_apb_amb_rhf

module subroutine cpscf_apb_amb1_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, nv, X, Qp, Qm)
use mod_xc_hess, only: xc_fxc_project_vo
use mod_integrals, only: integrals_build_exchange_cosx_batch_c2, integrals_build_coulomb_batch, &
                         integrals_build_coulomb_df_batch, integrals_build_exchange_batch, direct_mode
use mod_exchange, only: cosx_enabled, RS_omega
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc, nvir, nv
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: X(nvir,nocc,nv)
real(8),intent(out) :: Qp(nvir,nocc,nv), Qm(nvir,nocc,nv)
real(8),allocatable :: L2(:,:,:), KC(:,:,:), KC2(:,:,:), VC(:,:,:), Xk(:,:,:), Dn(:,:,:), Dp(:,:,:), Jb(:,:,:), &
                       Kb(:,:,:), Kd(:,:,:), T(:,:), CU(:,:)
integer :: c, i, a
real(8) :: tv0, hx
logical :: dmode0

if (nv < 1) return
hx = 0.5d0*HF_exchange_frac
allocate(CU(nConts,nocc), T(nConts,nocc))
if (cosx_enabled) then
   allocate(L2(nConts,nocc,nv), KC(nConts,nocc,nv), KC2(nConts,nocc,nv), VC(nConts,nocc,nv))
   do c = 1,nv
      call dgemm('N','N',nConts,nocc,nvir,2.0d0,C_a(1,nocc+1),nConts,X(1,1,c),nvir,0.0d0,L2(1,1,c),nConts)
   enddo
   tv0 = omp_get_wtime()
   call integrals_build_exchange_cosx_batch_c2(nConts, nv, L2, nocc, C_a, KC, KC2)
   cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   if (cpscf_triplet) then
      VC = 0.0d0
   else
      call cpscf_vind_c(nConts, nv, 0.0d0, nocc, C_a, L2, VC)
   endif
   do c = 1,nv
      VC(:,:,c) = VC(:,:,c) - hx*(KC(:,:,c) + KC2(:,:,c))
      call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,VC(1,1,c),nConts,0.0d0,Qp(1,1,c),nvir)
      KC2(:,:,c) = KC2(:,:,c) - KC(:,:,c)
      call dgemm('T','N',nvir,nocc,nConts,-hx,C_a(1,nocc+1),nConts,KC2(1,1,c),nConts,0.0d0,Qm(1,1,c),nvir)
   enddo
   deallocate(L2, KC, KC2, VC)
else
   if (RS_omega > 0.0d0) error stop 'TDDFT: range-separated hybrids need K=cosx'
   dmode0 = direct_mode
   direct_mode = .true.
   allocate(Dn(nConts,nConts,nv), Dp(nConts,nConts,nv), Jb(nConts,nConts,nv), Kb(nConts,nConts,nv), &
            Kd(nConts,nConts,nv))
   do c = 1,nv
      call dgemm('N','N',nConts,nocc,nvir,1.0d0,C_a(1,nocc+1),nConts,X(1,1,c),nvir,0.0d0,CU,nConts)
      call dgemm('N','T',nConts,nConts,nocc,2.0d0,CU,nConts,C_a,nConts,0.0d0,Dn(1,1,c),nConts)
      Dp(:,:,c) = Dn(:,:,c) + transpose(Dn(:,:,c))
   enddo
   tv0 = omp_get_wtime()
   call integrals_build_exchange_batch(nConts, nv, Dn, Dn, Kb, Kd)
   direct_mode = dmode0
   cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   tv0 = omp_get_wtime()
   if (cpscf_dfj_response()) then
      call integrals_build_coulomb_df_batch(nConts, nv, Dp, Jb)
   else
      call integrals_build_coulomb_batch(nConts, nv, Dp, Jb)
   endif
   cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
   do c = 1,nv
      Jb(:,:,c) = Jb(:,:,c) - hx*(Kb(:,:,c) + transpose(Kb(:,:,c)))
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Jb(1,1,c),nConts,C_a,nConts,0.0d0,T,nConts)
      call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,T,nConts,0.0d0,Qp(1,1,c),nvir)
      Kd(:,:,c) = Kb(:,:,c) - transpose(Kb(:,:,c))
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Kd(1,1,c),nConts,C_a,nConts,0.0d0,T,nConts)
      call dgemm('T','N',nvir,nocc,nConts,-hx,C_a(1,nocc+1),nConts,T,nConts,0.0d0,Qm(1,1,c),nvir)
   enddo
   deallocate(Dn, Dp, Jb, Kb, Kd)
endif
if (cpscf_xc_active) then
   allocate(Xk(nvir,nocc,nv))
   tv0 = omp_get_wtime()
   call apb_xc_vo(nConts, nocc, nvir, nv, C_a, X, Xk)
   cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
   Qp = Qp + Xk
   deallocate(Xk)
endif
do i = 1,nocc
   do a = 1,nvir
      Qp(a,i,:) = Qp(a,i,:) + (eLev_a(nocc+a) - eLev_a(i))*X(a,i,:)
      Qm(a,i,:) = Qm(a,i,:) + (eLev_a(nocc+a) - eLev_a(i))*X(a,i,:)
   enddo
enddo
deallocate(CU, T)
end subroutine cpscf_apb_amb1_rhf

subroutine apb_xc_vo(nConts, nocc, nvir, nv, C_a, X, Xk)
use mod_xc_hess, only: xc_fxc_project_vo, xc_fxc_project_vo_t
integer,intent(in) :: nConts, nocc, nvir, nv
real(8),intent(in) :: C_a(nConts,nConts), X(nvir,nocc,nv)
real(8),intent(out) :: Xk(nvir,nocc,nv)
if (cpscf_triplet) then
   call xc_fxc_project_vo_t(nConts, nocc, nvir, nv, C_a, X, cpscf_xc_Ptot, Xk)
else
   call xc_fxc_project_vo(nConts, nocc, nvir, nv, C_a, X, cpscf_xc_Ptot, Xk)
endif
end subroutine apb_xc_vo

module subroutine cpscf_apb_amb_uhf(nConts, noa, nob, Ca, ea, Cb, eb, HF_exchange_frac, np, Xp, Qp, nm, Xm, Qm)
use mod_integrals, only: integrals_build_exchange_cosx_batch_c2g, integrals_build_exchange_batch, &
                         integrals_build_exchange_lr_batch, direct_mode
use mod_exchange, only: cosx_enabled, RS_omega, RS_beta
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, noa, nob, np, nm
real(8),intent(in) :: Ca(nConts,nConts), ea(nConts), Cb(nConts,nConts), eb(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Xp((nConts-noa)*noa+(nConts-nob)*nob,np), Xm((nConts-noa)*noa+(nConts-nob)*nob,nm)
real(8),intent(out) :: Qp((nConts-noa)*noa+(nConts-nob)*nob,np), Qm((nConts-noa)*noa+(nConts-nob)*nob,nm)
integer :: nva, nvb, na, npr, c, i, a, k, s, no, nv, off
real(8) :: hx, tv0
logical :: exx, dmode0
real(8),allocatable :: Lsa(:,:,:), Lsb(:,:,:), VCa(:,:,:), VCb(:,:,:), T(:,:)
real(8),allocatable :: L1(:,:,:), L2(:,:,:), KC(:,:,:), KC2(:,:,:), KCa(:,:,:), KC2a(:,:,:), KCb(:,:,:), KC2b(:,:,:)
real(8),allocatable :: Ga(:,:,:), Gb(:,:,:), Ka(:,:,:), Kb(:,:,:), Kl(:,:,:), Klb(:,:,:), Kd(:,:)

nva = nConts - noa ; nvb = nConts - nob
na = nva*noa
npr = max(np, nm)
if (npr < 1) return
exx = (HF_exchange_frac /= 0.0d0 .or. RS_omega > 0.0d0)
Qp = 0.0d0 ; Qm = 0.0d0
if (exx) then
   allocate(KCa(nConts,noa,npr), KC2a(nConts,noa,npr), KCb(nConts,nob,npr), KC2b(nConts,nob,npr))
   if (cosx_enabled) then
      hx = HF_exchange_frac
      if (HF_exchange_frac == 0.0d0) hx = RS_beta
      do s = 1,2
         if (s == 1) then ; no = noa ; nv = nva ; off = 0 ; else ; no = nob ; nv = nvb ; off = na ; endif
         if (no == 0) cycle
         allocate(L1(nConts,no,npr), L2(nConts,no,npr), KC(nConts,no,npr), KC2(nConts,no,npr), T(nConts,no))
         do c = 1,npr
            call pair_factor(s, c, no, nv, off, L1(1,1,c), L2(1,1,c))
         enddo
         tv0 = omp_get_wtime()
         if (s == 1) then
            call integrals_build_exchange_cosx_batch_c2g(nConts, npr, L1, L2, no, Ca, KC, KC2)
         else
            call integrals_build_exchange_cosx_batch_c2g(nConts, npr, L1, L2, no, Cb, KC, KC2)
         endif
         cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
         do c = 1,npr
            T = 0.5d0*(KC(:,:,c) + KC2(:,:,c))
            KC2(:,:,c) = 0.5d0*(KC2(:,:,c) - KC(:,:,c))
            KC(:,:,c) = T
         enddo
         if (s == 1) then ; KCa = KC ; KC2a = KC2 ; else ; KCb = KC ; KC2b = KC2 ; endif
         deallocate(L1, L2, KC, KC2, T)
      enddo
   else
      hx = 1.0d0
      dmode0 = direct_mode
      direct_mode = .true.
      allocate(Ga(nConts,nConts,npr), Gb(nConts,nConts,npr), Ka(nConts,nConts,npr), Kb(nConts,nConts,npr))
      allocate(L1(nConts,max(noa,nob),npr), L2(nConts,max(noa,nob),npr))
      do s = 1,2
         if (s == 1) then ; no = noa ; nv = nva ; off = 0 ; else ; no = nob ; nv = nvb ; off = na ; endif
         do c = 1,npr
            call pair_factor(s, c, no, nv, off, L1(1,1,c), L2(1,1,c))
            if (s == 1) then
               call dgemm('N','T',nConts,nConts,no,1.0d0,L1(1,1,c),nConts,Ca,nConts,0.0d0,Ga(1,1,c),nConts)
               call dgemm('N','T',nConts,nConts,no,1.0d0,Ca,nConts,L2(1,1,c),nConts,1.0d0,Ga(1,1,c),nConts)
            else
               call dgemm('N','T',nConts,nConts,no,1.0d0,L1(1,1,c),nConts,Cb,nConts,0.0d0,Gb(1,1,c),nConts)
               call dgemm('N','T',nConts,nConts,no,1.0d0,Cb,nConts,L2(1,1,c),nConts,1.0d0,Gb(1,1,c),nConts)
            endif
         enddo
      enddo
      deallocate(L1, L2)
      tv0 = omp_get_wtime()
      if (HF_exchange_frac /= 0.0d0) then
         call integrals_build_exchange_batch(nConts, npr, Ga, Gb, Ka, Kb)
         Ka = HF_exchange_frac*Ka ; Kb = HF_exchange_frac*Kb
      else
         Ka = 0.0d0 ; Kb = 0.0d0
      endif
      if (RS_omega > 0.0d0) then
         allocate(Kl(nConts,nConts,npr), Klb(nConts,nConts,npr))
         call integrals_build_exchange_lr_batch(nConts, npr, Ga, Gb, Kl, Klb)
         Ka = Ka + RS_beta*Kl ; Kb = Kb + RS_beta*Klb
         deallocate(Kl, Klb)
      endif
      direct_mode = dmode0
      cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
      allocate(Kd(nConts,nConts))
      do c = 1,npr
         Kd = 0.5d0*(Ka(:,:,c) + transpose(Ka(:,:,c)))
         call dgemm('N','N',nConts,noa,nConts,1.0d0,Kd,nConts,Ca,nConts,0.0d0,KCa(1,1,c),nConts)
         Kd = 0.5d0*(Ka(:,:,c) - transpose(Ka(:,:,c)))
         call dgemm('N','N',nConts,noa,nConts,1.0d0,Kd,nConts,Ca,nConts,0.0d0,KC2a(1,1,c),nConts)
         Kd = 0.5d0*(Kb(:,:,c) + transpose(Kb(:,:,c)))
         call dgemm('N','N',nConts,nob,nConts,1.0d0,Kd,nConts,Cb,nConts,0.0d0,KCb(1,1,c),nConts)
         Kd = 0.5d0*(Kb(:,:,c) - transpose(Kb(:,:,c)))
         call dgemm('N','N',nConts,nob,nConts,1.0d0,Kd,nConts,Cb,nConts,0.0d0,KC2b(1,1,c),nConts)
      enddo
      deallocate(Ga, Gb, Ka, Kb, Kd)
   endif
endif
if (np > 0) then
   allocate(Lsa(nConts,noa,np), Lsb(nConts,nob,np), VCa(nConts,noa,np), VCb(nConts,nob,np))
   do c = 1,np
      call dgemm('N','N',nConts,noa,nva,1.0d0,Ca(1,noa+1),nConts,Xp(1,c),nva,0.0d0,Lsa(1,1,c),nConts)
      if (nob > 0) call dgemm('N','N',nConts,nob,nvb,1.0d0,Cb(1,nob+1),nConts,Xp(na+1,c),nvb,0.0d0,Lsb(1,1,c),nConts)
   enddo
   call cpscf_vind_spin_c(nConts, np, 0.0d0, noa, nob, Ca(:,1:noa), Cb(:,1:nob), Lsa, Lsb, VCa, VCb)
   do c = 1,np
      if (exx) then
         VCa(:,:,c) = VCa(:,:,c) - hx*KCa(:,:,c)
         if (nob > 0) VCb(:,:,c) = VCb(:,:,c) - hx*KCb(:,:,c)
      endif
      call dgemm('T','N',nva,noa,nConts,1.0d0,Ca(1,noa+1),nConts,VCa(1,1,c),nConts,0.0d0,Qp(1,c),nva)
      if (nob > 0) call dgemm('T','N',nvb,nob,nConts,1.0d0,Cb(1,nob+1),nConts,VCb(1,1,c),nConts,0.0d0,Qp(na+1,c),nvb)
   enddo
   deallocate(Lsa, Lsb, VCa, VCb)
endif
if (exx) then
   do c = 1,nm
      call dgemm('T','N',nva,noa,nConts,-hx,Ca(1,noa+1),nConts,KC2a(1,1,c),nConts,0.0d0,Qm(1,c),nva)
      if (nob > 0) call dgemm('T','N',nvb,nob,nConts,-hx,Cb(1,nob+1),nConts,KC2b(1,1,c),nConts,0.0d0,Qm(na+1,c),nvb)
   enddo
   deallocate(KCa, KC2a, KCb, KC2b)
endif
do c = 1,np
   do i = 1,noa
      do a = 1,nva
         k = a + (i-1)*nva
         Qp(k,c) = Qp(k,c) + (ea(noa+a) - ea(i))*Xp(k,c)
      enddo
   enddo
   do i = 1,nob
      do a = 1,nvb
         k = na + a + (i-1)*nvb
         Qp(k,c) = Qp(k,c) + (eb(nob+a) - eb(i))*Xp(k,c)
      enddo
   enddo
enddo
do c = 1,nm
   do i = 1,noa
      do a = 1,nva
         k = a + (i-1)*nva
         Qm(k,c) = Qm(k,c) + (ea(noa+a) - ea(i))*Xm(k,c)
      enddo
   enddo
   do i = 1,nob
      do a = 1,nvb
         k = na + a + (i-1)*nvb
         Qm(k,c) = Qm(k,c) + (eb(nob+a) - eb(i))*Xm(k,c)
      enddo
   enddo
enddo

contains

subroutine pair_factor(sp, cc, nos, nvs, ofs, F1, F2)
integer,intent(in) :: sp, cc, nos, nvs, ofs
real(8),intent(out) :: F1(nConts,nos), F2(nConts,nos)
real(8),allocatable :: Sv(:,:), Dv(:,:)
integer :: ii, aa
allocate(Sv(nvs,nos), Dv(nvs,nos))
Sv = 0.0d0 ; Dv = 0.0d0
do ii = 1,nos
   do aa = 1,nvs
      if (cc <= np) Sv(aa,ii) = Xp(ofs + aa + (ii-1)*nvs, cc)
      if (cc <= nm) Dv(aa,ii) = Xm(ofs + aa + (ii-1)*nvs, cc)
   enddo
enddo
if (sp == 1) then
   call dgemm('N','N',nConts,nos,nvs,1.0d0,Ca(1,noa+1),nConts,Sv+Dv,nvs,0.0d0,F1,nConts)
   call dgemm('N','N',nConts,nos,nvs,1.0d0,Ca(1,noa+1),nConts,Sv-Dv,nvs,0.0d0,F2,nConts)
else
   call dgemm('N','N',nConts,nos,nvs,1.0d0,Cb(1,nob+1),nConts,Sv+Dv,nvs,0.0d0,F1,nConts)
   call dgemm('N','N',nConts,nos,nvs,1.0d0,Cb(1,nob+1),nConts,Sv-Dv,nvs,0.0d0,F2,nConts)
endif
deallocate(Sv, Dv)
end subroutine pair_factor

end subroutine cpscf_apb_amb_uhf

module subroutine cpscf_response_begin_uks(nConts, Pa, Pb, HF_exchange_frac, scf_grid)
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
logical,intent(in),optional :: scf_grid
call cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac, scf_grid)
cpscf_xc_uks = (HF_exchange_frac < 1.0d0)
if (allocated(cpscf_xc_Pa)) deallocate(cpscf_xc_Pa, cpscf_xc_Pb)
if (cpscf_xc_uks) then
   allocate(cpscf_xc_Pa(nConts,nConts), cpscf_xc_Pb(nConts,nConts))
   cpscf_xc_Pa = Pa ; cpscf_xc_Pb = Pb
endif
end subroutine cpscf_response_begin_uks

end submodule cpscf_apb_impl
