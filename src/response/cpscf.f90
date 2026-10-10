! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: the coupled-perturbed self-consistent-field response

module mod_cpscf
implicit none
private
public :: cpscf_vind
public :: cpscf_vind_spin
public :: cpscf_debug_amat_cosx_symmetry
public :: cphf_solve_rhf, cpscf_solve_uhf, cphf_solve_efield
public :: cphf_debug_check_residual, cphf_debug_dump_dP, cphf_debug_grad_check
public :: cphf_hess_response, cphf_hess_response_spin, cphf_hess_nuc
public :: cphf_debug_hessian_row, cphf_debug_hessian_row_uhf, cphf_debug_e1_fixed
public :: cphf_debug_ccbk_raw, cphf_debug_dF_raw
public :: cphf_fock_deriv_atom, cpscf_match_scf_requested, cpscf_rij_mode, cpscf_warn_mixed_direct
public :: cphf_debug_amat_batch_check
public :: cphf_debug_batch_scaling
real(8),public :: cphf_t_amat = 0.0d0
real(8),public :: cphf_t_fockderiv = 0.0d0
real(8),public :: cphf_t_iter = 0.0d0
real(8),save :: cphf_c_iter = 0.0d0, cphf_c_pass1 = 0.0d0, cphf_c_pass3 = 0.0d0, cphf_c_fockderiv = 0.0d0, tstage_c = 0.0d0, ctmp_c
real(8),public :: cphf_t_vj = 0.0d0, cphf_t_vk = 0.0d0, cphf_t_vxc = 0.0d0
integer,public :: cphf_n_vk = 0, cphf_nd_vk = 0
integer,public :: cphf_n_iter_total = 0
real(8),public :: cphf_t_forceint = 0.0d0
real(8),public :: cphf_t_pass1 = 0.0d0
real(8),public :: cphf_t_pass3 = 0.0d0
real(8),allocatable,save :: cphf_voo_conn(:,:,:), cphf_voo_u(:,:,:)
logical,save :: cphf_voo_track = .false., cphf_voo_ready = .false.
real(8),public :: cphf_t_mo = 0.0d0

real(8),allocatable,save :: cpscf_xc_Ptot(:,:)
logical,save :: cpscf_xc_active = .false.
real(8),allocatable,save :: cpscf_cosx_kxa(:,:,:,:), cpscf_cosx_kxb(:,:,:,:)
logical,save :: cpscf_cosx_kx_proj = .false.
real(8),allocatable,save :: cpscf_lr_dKa(:,:,:), cpscf_lr_dKb(:,:,:)
real(8),allocatable,save :: cpscf_xc_vxco(:,:,:,:)
logical,save :: cpscf_xc_uks = .false.
logical,save :: cpscf_ecp_has = .false.
integer,save :: cpscf_ecp_atom = -1
real(8),allocatable,save :: cpscf_ecp_E(:,:,:)
real(8),allocatable,save :: cpscf_xc_Pa(:,:), cpscf_xc_Pb(:,:), cpscf_xc_vxca(:,:,:,:), cpscf_xc_vxcb(:,:,:,:)
real(8),allocatable,save :: cpscf_dfj_jx(:,:,:,:)
logical,save :: cpscf_force_dfj = .false.
real(8),save :: cphf_conv_override = -1.0d0

interface
module subroutine cpscf_hess_2e(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e, Cocc_a, Cocc_b)
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts), HF_exchange_frac
real(8),intent(out) :: d2E_2e(3,natoms,3,natoms)
real(8),intent(in),optional :: Cocc_a(:,:), Cocc_b(:,:)
end subroutine cpscf_hess_2e

module subroutine cpscf_hess_2e_release()
implicit none
end subroutine cpscf_hess_2e_release

module subroutine cpscf_cosx_hess_xchk(nConts, natoms, Pa, Pb, HF_exchange_frac, d2K_cosx, Cocc)
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts), HF_exchange_frac
real(8),intent(in) :: d2K_cosx(3,natoms,3,natoms)
real(8),intent(in) :: Cocc(:,:)
end subroutine cpscf_cosx_hess_xchk

module subroutine cphf_debug_hessian_row(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
implicit none
integer,intent(in) :: nConts, natoms, nocc, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_hessian_row

module subroutine cphf_debug_hessian_row_uhf(nConts, natoms, nocc_a, nocc_b, C_a, eLev_a, C_b, eLev_b, &
                                       Pa, Pb, HF_exchange_frac, iatom, ic)
implicit none
integer,intent(in) :: nConts, natoms, nocc_a, nocc_b, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: C_b(nConts,nConts), eLev_b(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_hessian_row_uhf

module subroutine cphf_build_amat(nConts, nocc, C_a, HF_exchange_frac, Amat)
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: Amat(nConts-nocc,nocc,nConts-nocc,nocc)
end subroutine cphf_build_amat

module subroutine cpscf_debug_amat_cosx_symmetry(nConts, nocc, C_a, HF_exchange_frac)
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cpscf_debug_amat_cosx_symmetry

module subroutine cphf_debug_amat_batch_check(nConts, nocc, C_a, HF_exchange_frac)
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_amat_batch_check

module subroutine cphf_debug_batch_scaling(nConts, nocc, C_a)
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
end subroutine cphf_debug_batch_scaling

module subroutine cphf_debug_check_residual(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_check_residual

module subroutine cphf_debug_dump_dP(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
implicit none
integer,intent(in) :: nConts, natoms, nocc, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_dump_dP

module subroutine cphf_debug_grad_check(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_grad_check

module subroutine cphf_debug_e1_fixed(nConts, Pa, Pb)
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
end subroutine cphf_debug_e1_fixed

module subroutine cphf_debug_ccbk_raw(nConts, natoms, iatom_c, mu, nu)
implicit none
integer,intent(in) :: nConts, natoms, iatom_c, mu, nu
end subroutine cphf_debug_ccbk_raw

module subroutine cphf_debug_dF_raw(nConts, natoms, iatom, ic, Pa, Pb, HF_exchange_frac)
implicit none
integer,intent(in) :: nConts, natoms, iatom, ic
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_dF_raw

end interface

contains

subroutine cphf_fock_deriv_atom(nConts, ao_atom, iatom, ic, &
                                 dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, &
                                 dF_full, dS_full)
use MOL_info, only: atoms, natoms
use mod_integrals, only: integrals_nuc_attraction_deriv, integrals_cphf_2e_density_role
use mod_xc_hess, only: xc_vxc_deriv_atom
use mod_profile, only: prof_start, prof_stop
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_Fx
implicit none
integer,intent(in) :: nConts, ao_atom(nConts), iatom, ic
real(8),intent(in) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)
real(8),intent(in) :: dJi(nConts,nConts,3), dKa(nConts,nConts,3)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: dF_full(nConts,nConts), dS_full(nConts,nConts)
INCLUDE 'parameter.h'

integer :: mu,nu
real(8),allocatable :: Dr(:,:,:)
real(8) :: Zeff, kfrac, jfac

kfrac = HF_exchange_frac
if (allocated(cpscf_cosx_kxa)) kfrac = 0.0d0
jfac = 1.0d0
if (allocated(cpscf_dfj_jx)) jfac = 0.0d0

dF_full = 0.0d0
dS_full = 0.0d0
do nu = 1,nConts
   do mu = 1,nConts
      if (ao_atom(mu) == iatom) then
         dF_full(mu,nu) = dF_full(mu,nu) - (dHcore(mu,nu,ic) + jfac*dJi(mu,nu,ic) - kfrac*dKa(mu,nu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(mu,nu,ic)
      endif
      if (ao_atom(nu) == iatom) then
         dF_full(mu,nu) = dF_full(mu,nu) - (dHcore(nu,mu,ic) + jfac*dJi(nu,mu,ic) - kfrac*dKa(nu,mu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(nu,mu,ic)
      endif
   enddo
enddo

allocate(Dr(nConts,nConts,3))
call prof_start("fd_nuc_attr_deriv")
call integrals_nuc_attraction_deriv(atoms(iatom)%coor*ans2bohr, nConts, Dr)
call prof_stop("fd_nuc_attr_deriv")
Zeff = atoms(iatom)%charge - atoms(iatom)%ecpCoreElec
dF_full = dF_full - Zeff*(Dr(:,:,ic) + transpose(Dr(:,:,ic)))
call cpscf_ecp_centre(nConts, iatom, Dr)
if (cpscf_ecp_has) dF_full = dF_full + Dr(:,:,ic) + transpose(Dr(:,:,ic))
deallocate(Dr)

if (.not. (jfac == 0.0d0 .and. kfrac == 0.0d0)) then
block
   integer, save :: cache_iatom = -1, cache_nConts = -1
   real(8), allocatable, save :: dJ2_cache(:,:,:), dK2_cache(:,:,:)
   real(8), allocatable, save :: Pa_cache(:,:), Pb_cache(:,:)
   logical :: cache_valid

   cache_valid = allocated(dJ2_cache) .and. (cache_iatom == iatom) .and. (cache_nConts == nConts)
   if (cache_valid) cache_valid = all(Pa_cache == Pa) .and. all(Pb_cache == Pb)

   if (.not. cache_valid) then
      if (allocated(dJ2_cache)) deallocate(dJ2_cache, dK2_cache, Pa_cache, Pb_cache)
      allocate(dJ2_cache(nConts,nConts,3), dK2_cache(nConts,nConts,3))
      allocate(Pa_cache(nConts,nConts), Pb_cache(nConts,nConts))
      call integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2_cache, dK2_cache)
      cache_iatom = iatom
      cache_nConts = nConts
      Pa_cache = Pa
      Pb_cache = Pb
   endif

   do nu = 1,nConts
      do mu = 1,nConts
         dF_full(mu,nu) = dF_full(mu,nu) &
            - (jfac*dJ2_cache(mu,nu,ic) - kfrac*dK2_cache(mu,nu,ic)) &
            - (jfac*dJ2_cache(nu,mu,ic) - kfrac*dK2_cache(nu,mu,ic))
      enddo
   enddo
end block
endif
if (allocated(cpscf_lr_dKa)) then
   block
      use mod_exchange, only: RS_beta
      real(8),allocatable :: dK2a(:,:,:), dK2b(:,:,:)
      allocate(dK2a(nConts,nConts,3), dK2b(nConts,nConts,3))
      call cpscf_lr_density_role(nConts, iatom, Pa, Pb, dK2a, dK2b)
      do nu = 1,nConts
         do mu = 1,nConts
            if (ao_atom(mu) == iatom) dF_full(mu,nu) = dF_full(mu,nu) + RS_beta*cpscf_lr_dKa(mu,nu,ic)
            if (ao_atom(nu) == iatom) dF_full(mu,nu) = dF_full(mu,nu) + RS_beta*cpscf_lr_dKa(nu,mu,ic)
            dF_full(mu,nu) = dF_full(mu,nu) + RS_beta*(dK2a(mu,nu,ic) + dK2a(nu,mu,ic))
         enddo
      enddo
      deallocate(dK2a, dK2b)
   end block
endif
if (allocated(cpscf_dfj_jx)) dF_full = dF_full + cpscf_dfj_jx(:,:,ic,iatom)
if (allocated(cpscf_cosx_kxa) .and. .not. cpscf_cosx_kx_proj) &
   dF_full = dF_full - HF_exchange_frac*cpscf_cosx_kxa(:,:,ic,iatom)

if (HF_exchange_frac < 1.0d0 .and. .not. allocated(cpscf_xc_vxco)) then
   block
      real(8),allocatable :: Ptot(:,:), dVxc(:,:)
      allocate(Ptot(nConts,nConts), dVxc(nConts,nConts))
      Ptot = Pa + Pb
      call prof_start("fd_xc_vxc_deriv")
      call xc_vxc_deriv_atom(nConts, natoms, Ptot, iatom, ic, dVxc)
      call prof_stop("fd_xc_vxc_deriv")
      dF_full = dF_full + dVxc
      deallocate(Ptot, dVxc)
   end block
endif
if (cosmo_hess_active) dF_full = dF_full + cosmo_hess_Fx(:,:,ic,iatom)
end subroutine cphf_fock_deriv_atom

subroutine cpscf_lr_density_role(nConts, iatom, Pa, Pb, dK2a, dK2b)
use mod_integrals, only: integrals_cphf_2e_density_role
use mod_exchange, only: RS_omega
implicit none
integer,intent(in) :: nConts, iatom
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: dK2a(nConts,nConts,3), dK2b(nConts,nConts,3)
integer,save :: c_iatom = -1, c_n = -1
real(8),allocatable,save :: c_a(:,:,:), c_b(:,:,:), c_Pa(:,:), c_Pb(:,:)
real(8),allocatable :: dJ2(:,:,:)
logical :: valid
valid = allocated(c_a) .and. (c_iatom == iatom) .and. (c_n == nConts)
if (valid) valid = all(c_Pa == Pa) .and. all(c_Pb == Pb)
if (.not. valid) then
   if (allocated(c_a)) deallocate(c_a, c_b, c_Pa, c_Pb)
   allocate(c_a(nConts,nConts,3), c_b(nConts,nConts,3), c_Pa(nConts,nConts), c_Pb(nConts,nConts))
   allocate(dJ2(nConts,nConts,3))
   call integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2, c_a, c_b, omega=RS_omega)
   deallocate(dJ2)
   c_iatom = iatom ; c_n = nConts ; c_Pa = Pa ; c_Pb = Pb
endif
dK2a = c_a ; dK2b = c_b
end subroutine cpscf_lr_density_role

subroutine cpscf_ecp_centre(nConts, iatom, Dr)
use mod_integrals, only: integrals_ecp_iprinv_atom
integer,intent(in) :: nConts, iatom
real(8),intent(inout) :: Dr(nConts,nConts,3)
if (iatom .ne. cpscf_ecp_atom .or. .not. allocated(cpscf_ecp_E)) then
   if (allocated(cpscf_ecp_E)) deallocate(cpscf_ecp_E)
   allocate(cpscf_ecp_E(nConts,nConts,3))
   call integrals_ecp_iprinv_atom(iatom, nConts, cpscf_ecp_E, cpscf_ecp_has)
   cpscf_ecp_atom = iatom
endif
if (cpscf_ecp_has) Dr = cpscf_ecp_E
end subroutine cpscf_ecp_centre

subroutine cphf_fock_deriv_atom_spin(nConts, ao_atom, iatom, ic, &
                                      dS, dHcore, dJi, dKa, dKb, HF_exchange_frac, Pa, Pb, &
                                      dF_a, dF_b, dS_full)
use MOL_info, only: atoms
use mod_integrals, only: integrals_nuc_attraction_deriv, integrals_cphf_2e_density_role
use mod_profile, only: prof_start, prof_stop
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_Fx
implicit none
integer,intent(in) :: nConts, ao_atom(nConts), iatom, ic
real(8),intent(in) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)
real(8),intent(in) :: dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: dF_a(nConts,nConts), dF_b(nConts,nConts), dS_full(nConts,nConts)
INCLUDE 'parameter.h'

integer :: mu,nu
real(8),allocatable :: Dr(:,:,:)
real(8) :: Zeff, kfrac, jfac

kfrac = HF_exchange_frac
if (allocated(cpscf_cosx_kxa)) kfrac = 0.0d0
jfac = 1.0d0
if (allocated(cpscf_dfj_jx)) jfac = 0.0d0

dF_a = 0.0d0
dF_b = 0.0d0
dS_full = 0.0d0
do nu = 1,nConts
   do mu = 1,nConts
      if (ao_atom(mu) == iatom) then
         dF_a(mu,nu) = dF_a(mu,nu) - (dHcore(mu,nu,ic) + jfac*dJi(mu,nu,ic) - kfrac*dKa(mu,nu,ic))
         dF_b(mu,nu) = dF_b(mu,nu) - (dHcore(mu,nu,ic) + jfac*dJi(mu,nu,ic) - kfrac*dKb(mu,nu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(mu,nu,ic)
      endif
      if (ao_atom(nu) == iatom) then
         dF_a(mu,nu) = dF_a(mu,nu) - (dHcore(nu,mu,ic) + jfac*dJi(nu,mu,ic) - kfrac*dKa(nu,mu,ic))
         dF_b(mu,nu) = dF_b(mu,nu) - (dHcore(nu,mu,ic) + jfac*dJi(nu,mu,ic) - kfrac*dKb(nu,mu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(nu,mu,ic)
      endif
   enddo
enddo

allocate(Dr(nConts,nConts,3))
call prof_start("fd_nuc_attr_deriv")
call integrals_nuc_attraction_deriv(atoms(iatom)%coor*ans2bohr, nConts, Dr)
call prof_stop("fd_nuc_attr_deriv")
Zeff = atoms(iatom)%charge - atoms(iatom)%ecpCoreElec
dF_a = dF_a - Zeff*(Dr(:,:,ic) + transpose(Dr(:,:,ic)))
dF_b = dF_b - Zeff*(Dr(:,:,ic) + transpose(Dr(:,:,ic)))
call cpscf_ecp_centre(nConts, iatom, Dr)
if (cpscf_ecp_has) then
   dF_a = dF_a + Dr(:,:,ic) + transpose(Dr(:,:,ic))
   dF_b = dF_b + Dr(:,:,ic) + transpose(Dr(:,:,ic))
endif
deallocate(Dr)

if (.not. (jfac == 0.0d0 .and. kfrac == 0.0d0)) then
block
   integer, save :: cache_iatom = -1, cache_nConts = -1
   real(8), allocatable, save :: dJ2_cache(:,:,:), dK2a_cache(:,:,:), dK2b_cache(:,:,:)
   real(8), allocatable, save :: Pa_cache(:,:), Pb_cache(:,:)
   logical :: cache_valid

   cache_valid = allocated(dJ2_cache) .and. (cache_iatom == iatom) .and. (cache_nConts == nConts)
   if (cache_valid) cache_valid = all(Pa_cache == Pa) .and. all(Pb_cache == Pb)

   if (.not. cache_valid) then
      if (allocated(dJ2_cache)) deallocate(dJ2_cache, dK2a_cache, dK2b_cache, Pa_cache, Pb_cache)
      allocate(dJ2_cache(nConts,nConts,3), dK2a_cache(nConts,nConts,3), dK2b_cache(nConts,nConts,3))
      allocate(Pa_cache(nConts,nConts), Pb_cache(nConts,nConts))
      call integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2_cache, dK2a_cache, dK2b_cache)
      cache_iatom = iatom
      cache_nConts = nConts
      Pa_cache = Pa
      Pb_cache = Pb
   endif

   do nu = 1,nConts
      do mu = 1,nConts
         dF_a(mu,nu) = dF_a(mu,nu) &
            - (jfac*dJ2_cache(mu,nu,ic) - kfrac*dK2a_cache(mu,nu,ic)) &
            - (jfac*dJ2_cache(nu,mu,ic) - kfrac*dK2a_cache(nu,mu,ic))
         dF_b(mu,nu) = dF_b(mu,nu) &
            - (jfac*dJ2_cache(mu,nu,ic) - kfrac*dK2b_cache(mu,nu,ic)) &
            - (jfac*dJ2_cache(nu,mu,ic) - kfrac*dK2b_cache(nu,mu,ic))
      enddo
   enddo
end block
endif
if (allocated(cpscf_lr_dKa)) then
   block
      use mod_exchange, only: RS_beta
      real(8),allocatable :: dK2a(:,:,:), dK2b(:,:,:)
      allocate(dK2a(nConts,nConts,3), dK2b(nConts,nConts,3))
      call cpscf_lr_density_role(nConts, iatom, Pa, Pb, dK2a, dK2b)
      do nu = 1,nConts
         do mu = 1,nConts
            if (ao_atom(mu) == iatom) then
               dF_a(mu,nu) = dF_a(mu,nu) + RS_beta*cpscf_lr_dKa(mu,nu,ic)
               dF_b(mu,nu) = dF_b(mu,nu) + RS_beta*cpscf_lr_dKb(mu,nu,ic)
            endif
            if (ao_atom(nu) == iatom) then
               dF_a(mu,nu) = dF_a(mu,nu) + RS_beta*cpscf_lr_dKa(nu,mu,ic)
               dF_b(mu,nu) = dF_b(mu,nu) + RS_beta*cpscf_lr_dKb(nu,mu,ic)
            endif
            dF_a(mu,nu) = dF_a(mu,nu) + RS_beta*(dK2a(mu,nu,ic) + dK2a(nu,mu,ic))
            dF_b(mu,nu) = dF_b(mu,nu) + RS_beta*(dK2b(mu,nu,ic) + dK2b(nu,mu,ic))
         enddo
      enddo
      deallocate(dK2a, dK2b)
   end block
endif
if (allocated(cpscf_dfj_jx)) then
   dF_a = dF_a + cpscf_dfj_jx(:,:,ic,iatom)
   dF_b = dF_b + cpscf_dfj_jx(:,:,ic,iatom)
endif
if (allocated(cpscf_cosx_kxa) .and. .not. cpscf_cosx_kx_proj) then
   dF_a = dF_a - HF_exchange_frac*cpscf_cosx_kxa(:,:,ic,iatom)
   dF_b = dF_b - HF_exchange_frac*cpscf_cosx_kxb(:,:,ic,iatom)
endif
if (cosmo_hess_active) then
   dF_a = dF_a + cosmo_hess_Fx(:,:,ic,iatom)
   dF_b = dF_b + cosmo_hess_Fx(:,:,ic,iatom)
endif
end subroutine cphf_fock_deriv_atom_spin

logical function cpscf_match_scf_requested() result(req)
logical, save :: checked = .false., val = .false.
character(32) :: envv
integer :: st
if (.not. checked) then
   call get_environment_variable("ENGINE_CPSCF_MATCH_SCF", envv, status=st)
   val = (st == 0)
   checked = .true.
endif
req = val
end function cpscf_match_scf_requested

logical function cpscf_dfj_response() result(r)
use mod_integrals, only: engine_use_df_j
r = engine_use_df_j .and. (cpscf_match_scf_requested() .or. allocated(cpscf_dfj_jx) .or. cpscf_force_dfj)
end function cpscf_dfj_response

logical function cpscf_rij_mode(HF_exchange_frac) result(r)
use mod_integrals, only: engine_use_df_j
use mod_exchange, only: cosx_enabled
real(8),intent(in) :: HF_exchange_frac
character(8) :: e
call get_environment_variable("ENGINE_CPSCF_EXACT_J", e)
r = engine_use_df_j .and. (cosx_enabled .or. HF_exchange_frac == 0.0d0) .and. (trim(e) .ne. "1")
end function cpscf_rij_mode

subroutine cpscf_warn_mixed_direct()
use mod_integrals, only: engine_use_df_j
if (.not. cpscf_match_scf_requested()) return
if (.not. engine_use_df_j) return
print *, 'WARNING: RI-J response requested (ENGINE_CPSCF_MATCH_SCF), but the CPSCF'
print *, '         right-hand side still differentiates EXACT J: this Hessian'
print *, '         mixes Coulomb approximations.'
end subroutine cpscf_warn_mixed_direct

subroutine cpscf_cosx_k(nConts, nd, D, K)
use mod_integrals, only: integrals_build_exchange_cosx, integrals_build_exchange_cosx_batch, &
                         cosx_set_response_mode
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nd
real(8),intent(in) :: D(nConts,nConts,nd)
real(8),intent(out) :: K(nConts,nConts,nd)
real(8),allocatable :: Kdisc(:,:), Kl(:,:,:)
character(len=8) :: e1, e2
integer :: b
real(8) :: tb
call get_environment_variable("ENGINE_CPSCF_COSX_LOOP", e1)
call get_environment_variable("ENGINE_COSX_BATCH_CHECK", e2)
if (trim(e1) .eq. "1" .or. trim(e2) .eq. "1") then
   allocate(Kdisc(nConts,nConts), Kl(nConts,nConts,nd))
   call cosx_set_response_mode(.true.)
   do b = 1,nd
      call integrals_build_exchange_cosx(nConts, D(:,:,b), D(:,:,b), Kl(:,:,b), Kdisc, .false.)
   enddo
   call cosx_set_response_mode(.false.)
   deallocate(Kdisc)
   if (trim(e2) .ne. "1") then
      K = Kl
      deallocate(Kl)
      return
   endif
endif
call get_environment_variable("ENGINE_COSX_BATCH_TIME", e1)
if (trim(e1) .eq. "1") then
   tb = omp_get_wtime()
   call integrals_build_exchange_cosx_batch(nConts, nd, D, K)
   print '(A,I0,A,F9.3,A)', " COSX BATCH nd=", nd, "  ", omp_get_wtime()-tb, " s"
   call flush(6)
else
   call integrals_build_exchange_cosx_batch(nConts, nd, D, K)
endif
if (trim(e2) .eq. "1") then
   print '(A,I0,A,ES10.3,A,ES10.3)', " COSX BATCH CHECK nd=", nd, "  max|K_batch - K_loop| = ", &
         maxval(abs(K - Kl)), "  max|K| = ", maxval(abs(Kl))
   deallocate(Kl)
endif
end subroutine cpscf_cosx_k

logical function cpscf_use_vind_c(HF_exchange_frac)
use mod_exchange, only: cosx_enabled
real(8),intent(in) :: HF_exchange_frac
character(len=8) :: e
call get_environment_variable("ENGINE_CPSCF_COSX_MIRROR", e)
cpscf_use_vind_c = cosx_enabled .and. HF_exchange_frac /= 0.0d0 .and. trim(e) .ne. "1"
call get_environment_variable("ENGINE_CPSCF_COSX_NONVAR", e)
if (trim(e) .eq. "0") cpscf_use_vind_c = .false.
if (HF_exchange_frac == 0.0d0 .and. cpscf_dfj_response()) then
   call get_environment_variable("ENGINE_CPSCF_DENSE_J", e)
   cpscf_use_vind_c = (trim(e) .ne. "1")
endif
end function cpscf_use_vind_c

subroutine cpscf_vind_c(nConts, nbatch, HF_exchange_frac, nocc, Co, Lf, VC)
use mod_integrals, only: integrals_build_coulomb_batch, integrals_build_coulomb_df_batch_c, &
                         integrals_build_exchange_cosx_batch_c
use omp_lib, only: omp_get_wtime
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_vind_c
implicit none
integer,intent(in) :: nConts, nbatch, nocc
real(8),intent(in) :: HF_exchange_frac, Co(nConts,nocc)
real(8),intent(in) :: Lf(nConts,nocc,nbatch)
real(8),intent(out) :: VC(nConts,nocc,nbatch)
real(8),allocatable :: Dch(:,:,:), Jch(:,:,:), KC(:,:,:)
integer :: b, c0, nch
real(8) :: tv0
character(len=8) :: e1
allocate(KC(nConts,nocc,nbatch))
tv0 = omp_get_wtime()
if (cpscf_dfj_response()) then
   call integrals_build_coulomb_df_batch_c(nConts, nbatch, Lf, nocc, Co, VC)
else
   allocate(Dch(nConts,nConts,16), Jch(nConts,nConts,16))
   do c0 = 0, nbatch-1, 16
      nch = min(16, nbatch-c0)
      do b = 1,nch
         call dgemm('N','T',nConts,nConts,nocc,1.0d0,Lf(1,1,c0+b),nConts,Co,nConts,0.0d0,Dch(1,1,b),nConts)
         Dch(:,:,b) = Dch(:,:,b) + transpose(Dch(:,:,b))
      enddo
      call integrals_build_coulomb_batch(nConts, nch, Dch(:,:,1:nch), Jch(:,:,1:nch))
      do b = 1,nch
         call dgemm('N','N',nConts,nocc,nConts,1.0d0,Jch(1,1,b),nConts,Co,nConts,0.0d0,VC(1,1,c0+b),nConts)
      enddo
   enddo
   deallocate(Dch, Jch)
endif
cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
if (HF_exchange_frac == 0.0d0) then
   deallocate(KC)
if (cosmo_hess_active) call cosmo_hess_vind_c(nConts, nbatch, nocc, Co, Lf, VC)
   return
endif
tv0 = omp_get_wtime()
call integrals_build_exchange_cosx_batch_c(nConts, nbatch, Lf, nocc, Co, KC)
call get_environment_variable("ENGINE_COSX_BATCH_TIME", e1)
if (trim(e1) .eq. "1") print '(A,I0,A,F9.3,A)', " COSX BATCH_C nd=", nbatch, "  ", omp_get_wtime()-tv0, " s"
cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
cphf_n_vk = cphf_n_vk + 1
cphf_nd_vk = cphf_nd_vk + nbatch
!$omp parallel do private(b)
do b = 1,nbatch
   VC(:,:,b) = VC(:,:,b) - 0.5d0*HF_exchange_frac*KC(:,:,b)
enddo
!$omp end parallel do
deallocate(KC)
if (cosmo_hess_active) call cosmo_hess_vind_c(nConts, nbatch, nocc, Co, Lf, VC)
end subroutine cpscf_vind_c

subroutine cpscf_vind(nConts, nbatch, HF_exchange_frac, Dtrial, Jbuf, Kbuf, skip_xc)
use mod_integrals, only: integrals_build_exchange_lr_batch, integrals_build_jk_batch, integrals_build_coulomb_batch, &
                         integrals_build_coulomb_df_batch, integrals_build_exchange_batch, &
                         integrals_build_exchange_cosx, cosx_set_response_mode, engine_use_df_j
use mod_exchange, only: cosx_enabled
use mod_xc_hess, only: xc_fxc_apply
use omp_lib, only: omp_get_wtime
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_add
implicit none
integer,intent(in) :: nConts, nbatch
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Dtrial(nConts,nConts,nbatch)
real(8),intent(inout) :: Jbuf(nConts,nConts,nbatch), Kbuf(nConts,nConts,nbatch)
logical,intent(in),optional :: skip_xc
integer :: b
real(8),allocatable :: Kdisc(:,:,:)
logical :: want_xc
real(8) :: tv0

if ((.not. cosx_enabled .or. HF_exchange_frac == 0.0d0) .and. &
    (.not. cpscf_dfj_response())) then
   if (HF_exchange_frac == 0.0d0) then
      call integrals_build_coulomb_batch(nConts, nbatch, Dtrial, Jbuf)
   else
      call integrals_build_jk_batch(nConts, nbatch, Dtrial, Jbuf, Kbuf)
      Jbuf = Jbuf - 0.5d0*HF_exchange_frac*Kbuf
   endif
   if (allocated(cpscf_lr_dKa)) then
      block
         use mod_exchange, only: RS_beta
         real(8),allocatable :: Kl(:,:,:), Kl2(:,:,:)
         allocate(Kl(nConts,nConts,nbatch), Kl2(nConts,nConts,nbatch))
         call integrals_build_exchange_lr_batch(nConts, nbatch, Dtrial(:,:,1:nbatch), Dtrial(:,:,1:nbatch), Kl, Kl2)
         Jbuf = Jbuf - 0.5d0*RS_beta*Kl
         deallocate(Kl, Kl2)
      end block
   endif
else
   if (cpscf_dfj_response()) then
      tv0 = omp_get_wtime()
      call integrals_build_coulomb_df_batch(nConts, nbatch, Dtrial, Jbuf)
      cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
   else
      tv0 = omp_get_wtime()
      call integrals_build_coulomb_batch(nConts, nbatch, Dtrial, Jbuf)
      cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
   endif

   if (cosx_enabled) then
      if (HF_exchange_frac /= 0.0d0) then
         tv0 = omp_get_wtime()
         call cpscf_cosx_k(nConts, nbatch, Dtrial, Kbuf)
         cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
         cphf_n_vk = cphf_n_vk + 1
         cphf_nd_vk = cphf_nd_vk + nbatch
         Jbuf = Jbuf - 0.5d0*HF_exchange_frac*Kbuf
      endif
   else
      if (HF_exchange_frac /= 0.0d0) then
         allocate(Kdisc(nConts,nConts,nbatch))
         call integrals_build_exchange_batch(nConts, nbatch, Dtrial, Dtrial, Kbuf, Kdisc)
         deallocate(Kdisc)
         Jbuf = Jbuf - 0.5d0*HF_exchange_frac*Kbuf
      endif
   endif
endif

want_xc = cpscf_xc_active
if (present(skip_xc)) then
   if (skip_xc) want_xc = .false.
endif
if (want_xc) then
   if (.not. allocated(cpscf_xc_Ptot)) then
      print *,'ERROR: cpscf_xc_active but Ptot not set'; stop 1
   endif
   call xc_fxc_apply(nConts, nbatch, cpscf_xc_Ptot, Dtrial, Jbuf)
endif
if (cosmo_hess_active) call cosmo_hess_add(nConts, nbatch, Dtrial, Jbuf)
end subroutine cpscf_vind

subroutine cpscf_vind_spin(nConts, nbatch, nspin, HF_exchange_frac, Dtrial, Dtot, Vout, Kscr)
use mod_integrals, only: integrals_build_exchange_lr_batch, integrals_build_jk_batch_uhf, integrals_build_coulomb_df_batch, &
                         integrals_build_coulomb_batch, integrals_build_exchange_batch, &
                         integrals_build_exchange_cosx, cosx_set_response_mode, engine_use_df_j
use mod_exchange, only: cosx_enabled
use mod_xc_hess, only: xc_fxc_apply_u
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_add_spin
implicit none
integer,intent(in) :: nConts, nbatch, nspin
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Dtrial(nConts,nConts,nbatch,nspin)
real(8),intent(in) :: Dtot(nConts,nConts,nbatch)
real(8),intent(inout) :: Vout(nConts,nConts,nbatch,nspin)
real(8),intent(inout) :: Kscr(nConts,nConts,nbatch,nspin)
integer :: b

if (nspin .ne. 2) then
   print *, 'ERROR: cpscf_vind_spin only implements nspin=2 (UHF); got', nspin
   stop 1
endif

if ((.not. cosx_enabled .or. HF_exchange_frac == 0.0d0) .and. &
    (.not. cpscf_dfj_response())) then
   call integrals_build_jk_batch_uhf(nConts, nbatch, Dtot, &
             Dtrial(:,:,:,1), Dtrial(:,:,:,2), &
             Vout(:,:,:,1), Kscr(:,:,:,1), Kscr(:,:,:,2))
   Vout(:,:,:,2) = Vout(:,:,:,1)
   Vout(:,:,:,1) = Vout(:,:,:,1) - HF_exchange_frac*Kscr(:,:,:,1)
   Vout(:,:,:,2) = Vout(:,:,:,2) - HF_exchange_frac*Kscr(:,:,:,2)
   if (allocated(cpscf_lr_dKa)) then
      block
         use mod_exchange, only: RS_beta
         real(8),allocatable :: Kla(:,:,:), Klb(:,:,:)
         allocate(Kla(nConts,nConts,nbatch), Klb(nConts,nConts,nbatch))
         call integrals_build_exchange_lr_batch(nConts, nbatch, Dtrial(:,:,1:nbatch,1), Dtrial(:,:,1:nbatch,2), Kla, Klb)
         Vout(:,:,1:nbatch,1) = Vout(:,:,1:nbatch,1) - RS_beta*Kla
         Vout(:,:,1:nbatch,2) = Vout(:,:,1:nbatch,2) - RS_beta*Klb
         deallocate(Kla, Klb)
      end block
   endif
   if (cpscf_xc_uks) call xc_fxc_apply_u(nConts, nbatch, cpscf_xc_Pa, cpscf_xc_Pb, &
                                         Dtrial(:,:,:,1), Dtrial(:,:,:,2), Vout(:,:,:,1), Vout(:,:,:,2))
if (cosmo_hess_active) call cosmo_hess_add_spin(nConts, nbatch, Dtot, Vout)
   return
endif

if (cpscf_dfj_response()) then
   call integrals_build_coulomb_df_batch(nConts, nbatch, Dtot, Vout(:,:,:,1))
else
   call integrals_build_coulomb_batch(nConts, nbatch, Dtot, Vout(:,:,:,1))
endif

if (cosx_enabled) then
   if (HF_exchange_frac /= 0.0d0) then
      call cpscf_cosx_k(nConts, 2*nbatch, Dtrial, Kscr)
   else
      Kscr = 0.0d0
   endif
else
   if (HF_exchange_frac /= 0.0d0) then
      call integrals_build_exchange_batch(nConts, nbatch, Dtrial(:,:,:,1), Dtrial(:,:,:,2), &
                                           Kscr(:,:,:,1), Kscr(:,:,:,2))
   else
      Kscr = 0.0d0
   endif
endif
Vout(:,:,:,2) = Vout(:,:,:,1)
Vout(:,:,:,1) = Vout(:,:,:,1) - HF_exchange_frac*Kscr(:,:,:,1)
Vout(:,:,:,2) = Vout(:,:,:,2) - HF_exchange_frac*Kscr(:,:,:,2)
if (cpscf_xc_uks) call xc_fxc_apply_u(nConts, nbatch, cpscf_xc_Pa, cpscf_xc_Pb, &
                                      Dtrial(:,:,:,1), Dtrial(:,:,:,2), Vout(:,:,:,1), Vout(:,:,:,2))
if (cosmo_hess_active) call cosmo_hess_add_spin(nConts, nbatch, Dtot, Vout)
end subroutine cpscf_vind_spin

subroutine cphf_solve_vo_pcg(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                             RHSfix, Uvo_all, niter_total, conv_thr)
use mod_xc_hess, only: xc_fxc_project_vo
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc, nvir, npert
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac, conv_thr
real(8),intent(in) :: RHSfix(nvir,nocc,npert)
real(8),intent(out) :: Uvo_all(nvir,nocc,npert)
integer,intent(out) :: niter_total
real(8),allocatable :: dlt(:,:), R(:,:,:), Z(:,:,:), Pd(:,:,:), Q(:,:,:), rz(:), Qoo(:,:,:)
logical,allocatable :: conv(:)
integer,allocatable :: act(:)
integer :: p, it, nact, k, a, i
real(8) :: alpha, beta, rznew, pq

allocate(dlt(nvir,nocc), R(nvir,nocc,npert), Z(nvir,nocc,npert), Pd(nvir,nocc,npert))
allocate(Q(nvir,nocc,npert), rz(npert), conv(npert), act(npert))
if (cphf_voo_track) allocate(Qoo(nocc,nocc,npert))
do i = 1,nocc
   do a = 1,nvir
      dlt(a,i) = eLev_a(nocc+a) - eLev_a(i)
   enddo
enddo
block
   character(len=16) :: solv
   call get_environment_variable("ENGINE_CPHF_SOLVER", solv)
   if (trim(solv) .ne. "pcg") then
      call bcg_rhf()
      deallocate(dlt, R, Z, Pd, Q, rz, conv, act)
      if (cphf_voo_track) then
         deallocate(Qoo)
         cphf_voo_ready = .true.
      endif
      return
   endif
end block
do p = 1,npert
   Uvo_all(:,:,p) = -RHSfix(:,:,p)/dlt
   act(p) = p
enddo
call apply_A(npert, act, Uvo_all, Q)
if (cphf_voo_track) cphf_voo_u = Qoo
niter_total = npert
do p = 1,npert
   R(:,:,p) = -RHSfix(:,:,p) - Q(:,:,p)
   Z(:,:,p) = R(:,:,p)/dlt
   Pd(:,:,p) = Z(:,:,p)
   rz(p) = sum(R(:,:,p)*Z(:,:,p))
   conv(p) = maxval(abs(Z(:,:,p))) < conv_thr
enddo
do it = 1,200
   nact = 0
   do p = 1,npert
      if (.not. conv(p)) then
         nact = nact + 1
         act(nact) = p
      endif
   enddo
   if (nact .eq. 0) exit
   niter_total = niter_total + nact
   call apply_A(nact, act, Pd, Q)
   do k = 1,nact
      p = act(k)
      pq = sum(Pd(:,:,p)*Q(:,:,p))
      alpha = rz(p)/pq
      Uvo_all(:,:,p) = Uvo_all(:,:,p) + alpha*Pd(:,:,p)
      if (cphf_voo_track) cphf_voo_u(:,:,p) = cphf_voo_u(:,:,p) + alpha*Qoo(:,:,p)
      R(:,:,p) = R(:,:,p) - alpha*Q(:,:,p)
      Z(:,:,p) = R(:,:,p)/dlt
      rznew = sum(R(:,:,p)*Z(:,:,p))
      beta = rznew/rz(p)
      rz(p) = rznew
      Pd(:,:,p) = Z(:,:,p) + beta*Pd(:,:,p)
      conv(p) = maxval(abs(Z(:,:,p))) < conv_thr
   enddo
enddo
deallocate(dlt, R, Z, Pd, Q, rz, conv, act)
if (cphf_voo_track) then
   deallocate(Qoo)
   cphf_voo_ready = .true.
endif

contains

subroutine bcg_rhf()
integer :: nt_, it_, k_, p_, q_, nk_, info_, lwork_
real(8),allocatable :: G_(:,:), ev_(:), work_(:), T_(:,:), W_(:,:), Wo(:,:)
logical :: verb_
character(len=8) :: e_
call get_environment_variable("ENGINE_CPHF_VERBOSE", e_)
verb_ = (trim(e_) .eq. "1")
nt_ = nvir*nocc
do p_ = 1,npert
   Uvo_all(:,:,p_) = -RHSfix(:,:,p_)/dlt
   act(p_) = p_
enddo
call apply_A(npert, act, Uvo_all, Q)
if (cphf_voo_track) cphf_voo_u = Qoo
niter_total = npert
R = -RHSfix - Q
nk_ = npert
do it_ = 1,200
   do p_ = 1,npert
      Z(:,:,p_) = R(:,:,p_)/dlt
      conv(p_) = maxval(abs(Z(:,:,p_))) < conv_thr
      if (conv(p_)) Z(:,:,p_) = 0.0d0
   enddo
   if (verb_) print '(A,I4,A,I5)', '  RHF BCG it', it_, '  unconverged', count(.not. conv)
   if (all(conv)) exit
   if (it_ .eq. 1) then
      Pd = Z
   else
      allocate(T_(nk_,npert), W_(nt_,npert))
      call dgemm('T','N',nk_,npert,nt_,1.0d0,Q,nt_,Z,nt_,0.0d0,T_,nk_)
      W_ = reshape(Z,(/nt_,npert/))
      call dgemm('N','N',nt_,npert,nk_,-1.0d0,Pd,nt_,T_,nk_,1.0d0,W_,nt_)
      Pd = reshape(W_,(/nvir,nocc,npert/))
      deallocate(T_, W_)
   endif
   q_ = 0
   do p_ = 1,npert
      if (.not. conv(p_)) then
         q_ = q_ + 1
         if (q_ .ne. p_) Pd(:,:,q_) = Pd(:,:,p_)
      endif
   enddo
   nk_ = q_
   do k_ = 1,nk_
      act(k_) = k_
   enddo
   call apply_A(nk_, act, Pd, Q)
   niter_total = niter_total + nk_
   allocate(G_(nk_,nk_), ev_(nk_), work_(1))
   call dgemm('T','N',nk_,nk_,nt_,1.0d0,Pd,nt_,Q,nt_,0.0d0,G_,nk_)
   G_ = 0.5d0*(G_ + transpose(G_))
   call dsyev('V','U',nk_,G_,nk_,ev_,work_,-1,info_)
   lwork_ = int(work_(1)) ; deallocate(work_) ; allocate(work_(lwork_))
   call dsyev('V','U',nk_,G_,nk_,ev_,work_,lwork_,info_)
   q_ = 0
   do k_ = 1,nk_
      if (ev_(k_) .gt. 1.0d-12*ev_(nk_)) then
         q_ = q_ + 1
         G_(:,q_) = G_(:,k_)/sqrt(ev_(k_))
      endif
   enddo
   allocate(W_(nt_,q_))
   call dgemm('N','N',nt_,q_,nk_,1.0d0,Pd,nt_,G_,nk_,0.0d0,W_,nt_)
   Pd(:,:,1:q_) = reshape(W_,(/nvir,nocc,q_/))
   call dgemm('N','N',nt_,q_,nk_,1.0d0,Q,nt_,G_,nk_,0.0d0,W_,nt_)
   Q(:,:,1:q_) = reshape(W_,(/nvir,nocc,q_/))
   deallocate(W_)
   if (cphf_voo_track) then
      allocate(Wo(nocc*nocc,q_))
      call dgemm('N','N',nocc*nocc,q_,nk_,1.0d0,Qoo,nocc*nocc,G_,nk_,0.0d0,Wo,nocc*nocc)
      Qoo(:,:,1:q_) = reshape(Wo,(/nocc,nocc,q_/))
      deallocate(Wo)
   endif
   deallocate(G_, ev_, work_)
   nk_ = q_
   allocate(T_(nk_,npert))
   call dgemm('T','N',nk_,npert,nt_,1.0d0,Pd,nt_,R,nt_,0.0d0,T_,nk_)
   call dgemm('N','N',nt_,npert,nk_,1.0d0,Pd,nt_,T_,nk_,1.0d0,Uvo_all,nt_)
   call dgemm('N','N',nt_,npert,nk_,-1.0d0,Q,nt_,T_,nk_,1.0d0,R,nt_)
   if (cphf_voo_track) call dgemm('N','N',nocc*nocc,npert,nk_,1.0d0,Qoo,nocc*nocc,T_,nk_,1.0d0,cphf_voo_u,nocc*nocc)
   deallocate(T_)
enddo
end subroutine bcg_rhf

subroutine apply_A(n_act, alist, X, Qout)
integer,intent(in) :: n_act, alist(n_act)
real(8),intent(in) :: X(nvir,nocc,npert)
real(8),intent(inout) :: Qout(nvir,nocc,npert)
real(8),allocatable :: Dtr(:,:,:), Jtr(:,:,:), Katr(:,:,:), Xp(:,:,:), Ua(:,:,:), Tmp(:,:), T2(:,:), CU(:,:)
real(8),allocatable :: VC(:,:,:), LU(:,:,:)
logical :: vcm
integer :: c, pp
real(8) :: tv0
vcm = cpscf_use_vind_c(HF_exchange_frac)
if (vcm) then
   allocate(Dtr(1,1,1), Jtr(1,1,1), Katr(1,1,1), LU(nConts,nocc,n_act))
else
   allocate(Dtr(nConts,nConts,n_act), Jtr(nConts,nConts,n_act), Katr(nConts,nConts,n_act))
endif
allocate(Tmp(nConts,nocc), T2(nConts,nConts), CU(nConts,nocc))
do c = 1,n_act
   pp = alist(c)
   call dgemm('N','N',nConts,nocc,nvir,1.0d0,C_a(1,nocc+1),nConts,X(1,1,pp),nvir,0.0d0,CU,nConts)
   if (vcm) then
      LU(:,:,c) = 2.0d0*CU
   else
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,CU,nConts,C_a,nConts,0.0d0,T2,nConts)
      Dtr(:,:,c) = 2.0d0*(T2 + transpose(T2))
   endif
enddo
if (vcm) then
   allocate(VC(nConts,nocc,n_act))
   call cpscf_vind_c(nConts, n_act, HF_exchange_frac, nocc, C_a, LU, VC)
else
   call cpscf_vind(nConts, n_act, HF_exchange_frac, Dtr, Jtr, Katr, skip_xc=.true.)
endif
if (cpscf_xc_active) then
   allocate(Xp(nvir,nocc,n_act), Ua(nvir,nocc,n_act))
   do c = 1,n_act
      Ua(:,:,c) = X(:,:,alist(c))
   enddo
   tv0 = omp_get_wtime()
   call xc_fxc_project_vo(nConts, nocc, nvir, n_act, C_a, Ua, cpscf_xc_Ptot, Xp)
   cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
endif
do c = 1,n_act
   pp = alist(c)
   if (vcm) then
      Tmp = VC(:,:,c)
   else
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Jtr(1,1,c),nConts,C_a,nConts,0.0d0,Tmp,nConts)
   endif
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,Tmp,nConts,0.0d0,Qout(1,1,pp),nvir)
   if (cphf_voo_track) call dgemm('T','N',nocc,nocc,nConts,1.0d0,C_a,nConts,Tmp,nConts,0.0d0,Qoo(1,1,pp),nocc)
   if (cpscf_xc_active) Qout(:,:,pp) = Qout(:,:,pp) + Xp(:,:,c)
   Qout(:,:,pp) = dlt*X(:,:,pp) + Qout(:,:,pp)
enddo
deallocate(Dtr, Jtr, Katr, Tmp, T2, CU)
if (allocated(Xp)) deallocate(Xp, Ua)
if (allocated(VC)) deallocate(VC, LU)
end subroutine apply_A
end subroutine cphf_solve_vo_pcg

subroutine cphf_solve_vo_direct(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                                 RHSfix, Uvo_all, niter_total)
use mod_integrals, only: direct_mode
use mod_xc_hess, only: xc_fxc_project_vo
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc, nvir, npert
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: RHSfix(nvir,nocc,npert)
real(8),intent(out) :: Uvo_all(nvir,nocc,npert)
integer,intent(out) :: niter_total

integer(8), parameter :: cache_budget_bytes = 2_8*1024_8*1024_8
integer(8), parameter :: diis_budget_bytes  = 256_8*1024_8*1024_8
integer :: bchunk, ndiis
integer :: iter, p, i, a, cstart, cend, nb, cidx, nact
integer :: info, nuse, ihist, r, s
real(8) :: dmax_p, tv0
real(8) :: conv_thr
integer, parameter :: NTHR = 5
real(8), parameter :: thr_list(NTHR) = (/1.0d-4,1.0d-6,1.0d-8,1.0d-10,1.0d-12/)
integer :: thr_passes(NTHR), it
logical,allocatable :: thr_hit(:,:)
logical,allocatable :: converged(:)
integer,allocatable :: actlist(:), nhist(:), ipiv(:)
real(8),allocatable :: Dtr(:,:,:), Jtr(:,:,:), Katr(:,:,:)
real(8),allocatable :: Tmp(:,:), Unew(:,:,:)
real(8),allocatable :: dsU(:,:,:,:), dsE(:,:,:,:)
real(8),allocatable :: Bmat(:,:), rhs(:)
real(8),allocatable :: Xproj(:,:,:), Uact(:,:,:)

if (direct_mode) then
   bchunk = npert
else
   bchunk = int(max(1_8, cache_budget_bytes / (2_8*int(nConts,8)*int(nConts,8)*8_8)))
   bchunk = max(1, min(npert, bchunk))
endif

ndiis = int(max(1_8, diis_budget_bytes / &
         (2_8*int(nvir,8)*int(nocc,8)*int(npert,8)*8_8)))
ndiis = max(2, min(10, ndiis))

allocate(Dtr(nConts,nConts,bchunk), Jtr(nConts,nConts,bchunk))
allocate(Katr(nConts,nConts,bchunk))
allocate(Tmp(nConts,nocc))
allocate(Unew(nvir,nocc,npert), converged(npert), actlist(npert), nhist(npert))
allocate(Xproj(nvir,nocc,bchunk), Uact(nvir,nocc,bchunk))
allocate(dsU(nvir,nocc,ndiis,npert), dsE(nvir,nocc,ndiis,npert))

Uvo_all = 0.0d0
converged = .false.
nhist = 0
niter_total = 0
allocate(thr_hit(NTHR,npert))
thr_hit = .false.
thr_passes = 0

conv_thr = 5.0d-5
block
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_CPHF_CONV", envval, status=envstat)
   if (envstat .eq. 0) read(envval,*) conv_thr
   if (cphf_conv_override > 0.0d0) conv_thr = cphf_conv_override
end block
block
   character(16) :: solv
   call get_environment_variable("ENGINE_CPHF_SOLVER", solv)
   if (trim(solv) .ne. "diis") then
      call cphf_solve_vo_pcg(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                             RHSfix, Uvo_all, niter_total, conv_thr)
      return
   endif
end block

do iter = 1,200
   nact = 0
   do p = 1,npert
      if (.not. converged(p)) then
         nact = nact + 1
         actlist(nact) = p
      endif
   enddo
   if (nact .eq. 0) exit
   niter_total = niter_total + nact

   cstart = 1
   do while (cstart .le. nact)
      cend = min(nact, cstart+bchunk-1)
      nb = cend - cstart + 1
      do cidx = 1,nb
         p = actlist(cstart+cidx-1)
         Dtr(:,:,cidx) = 2.0d0*( matmul(C_a(:,nocc+1:nConts), matmul(Uvo_all(:,:,p), transpose(C_a(:,1:nocc)))) &
                       + matmul(C_a(:,1:nocc), matmul(transpose(Uvo_all(:,:,p)), transpose(C_a(:,nocc+1:nConts)))) )
      enddo
      call cpscf_vind(nConts, nb, HF_exchange_frac, Dtr(:,:,1:nb), Jtr(:,:,1:nb), Katr(:,:,1:nb), skip_xc=.true.)
      if (cpscf_xc_active) then
         do cidx = 1,nb
            p = actlist(cstart+cidx-1)
            Uact(:,:,cidx) = Uvo_all(:,:,p)
         enddo
         tv0 = omp_get_wtime()
         call xc_fxc_project_vo(nConts, nocc, nvir, nb, C_a, Uact(:,:,1:nb), cpscf_xc_Ptot, Xproj(:,:,1:nb))
         cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
      endif
      do cidx = 1,nb
         p = actlist(cstart+cidx-1)
         Tmp = matmul(Jtr(:,:,cidx), C_a(:,1:nocc))
         do i = 1,nocc
            do a = 1,nvir
               if (cpscf_xc_active) then
                  Unew(a,i,p) = ( RHSfix(a,i,p) + dot_product(C_a(:,nocc+a), Tmp(:,i)) &
                                  + Xproj(a,i,cidx) ) / (eLev_a(i)-eLev_a(nocc+a))
               else
                  Unew(a,i,p) = ( RHSfix(a,i,p) + dot_product(C_a(:,nocc+a), Tmp(:,i)) ) &
                                / (eLev_a(i)-eLev_a(nocc+a))
               endif
            enddo
         enddo
      enddo
      cstart = cend + 1
   enddo

   do cidx = 1,nact
      p = actlist(cidx)
      dmax_p = maxval(abs(Unew(:,:,p)-Uvo_all(:,:,p)))

      ihist = mod(nhist(p), ndiis) + 1
      dsU(:,:,ihist,p) = Unew(:,:,p)
      dsE(:,:,ihist,p) = Unew(:,:,p) - Uvo_all(:,:,p)
      nhist(p) = nhist(p) + 1
      Uvo_all(:,:,p) = Unew(:,:,p)

      nuse = min(nhist(p), ndiis)
      if (nuse .ge. 2) then
         allocate(Bmat(nuse+1,nuse+1), rhs(nuse+1), ipiv(nuse+1))
         do r = 1,nuse
            do s = 1,nuse
               Bmat(r,s) = sum(dsE(:,:,r,p)*dsE(:,:,s,p))
            enddo
         enddo
         Bmat(1:nuse,nuse+1) = -1.0d0
         Bmat(nuse+1,1:nuse) = -1.0d0
         Bmat(nuse+1,nuse+1) = 0.0d0
         rhs = 0.0d0
         rhs(nuse+1) = -1.0d0
         call dgesv(nuse+1, 1, Bmat, nuse+1, ipiv, rhs, nuse+1, info)
         if (info .eq. 0) then
            Uvo_all(:,:,p) = 0.0d0
            do r = 1,nuse
               Uvo_all(:,:,p) = Uvo_all(:,:,p) + rhs(r)*dsU(:,:,r,p)
            enddo
         endif
         deallocate(Bmat, rhs, ipiv)
      endif

      do it = 1,NTHR
         if (.not. thr_hit(it,p)) then
            thr_passes(it) = thr_passes(it) + 1
            if (dmax_p < thr_list(it)) thr_hit(it,p) = .true.
         endif
      enddo

      if (dmax_p < conv_thr) converged(p) = .true.
   enddo
enddo

block
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_CPHF_THRESHOLD_SCAN", envval, status=envstat)
   if (envstat .eq. 0) then
      print *, ' CPHF passes vs convergence threshold (total over', npert, 'perturbations):'
      do it = 1,NTHR
         print *, '   thr=', thr_list(it), ' passes=', thr_passes(it)
      enddo
   endif
end block
deallocate(thr_hit)

deallocate(Dtr,Jtr,Katr,Tmp,Unew,converged,actlist,nhist,dsU,dsE)
deallocate(Xproj,Uact)
end subroutine cphf_solve_vo_direct

subroutine cpscf_solve_vo_direct_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, &
                                       C_a, eLev_a, C_b, eLev_b, HF_exchange_frac, &
                                       RHSfix_a, RHSfix_b, Uvo_a, Uvo_b, niter_total)
use mod_integrals, only: direct_mode
implicit none
integer,intent(in) :: nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: C_b(nConts,nConts), eLev_b(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: RHSfix_a(nvir_a,nocc_a,npert), RHSfix_b(nvir_b,nocc_b,npert)
real(8),intent(out) :: Uvo_a(nvir_a,nocc_a,npert), Uvo_b(nvir_b,nocc_b,npert)
integer,intent(out) :: niter_total

integer(8), parameter :: cache_budget_bytes = 2_8*1024_8*1024_8
integer(8), parameter :: diis_budget_bytes  = 256_8*1024_8*1024_8
integer :: bchunk, ndiis, nuvo
integer :: iter, p, i, a, cstart, cend, nb, cidx, nact
integer :: info, nuse, ihist, r, s, off
real(8) :: dmax_p
real(8) :: conv_thr
integer, parameter :: NTHR = 5
real(8), parameter :: thr_list(NTHR) = (/1.0d-4,1.0d-6,1.0d-8,1.0d-10,1.0d-12/)
integer :: thr_passes(NTHR), it
logical,allocatable :: thr_hit(:,:)
logical,allocatable :: converged(:)
integer,allocatable :: actlist(:), nhist(:), ipiv(:)
real(8),allocatable :: Dtr(:,:,:,:), Dtot(:,:,:), Vtr(:,:,:,:), Kscr(:,:,:,:)
real(8),allocatable :: Tmpa(:,:), Tmpb(:,:), Unew_a(:,:,:), Unew_b(:,:,:)
real(8),allocatable :: dsU(:,:,:), dsE(:,:,:), flat_u(:), flat_unew(:), flat_err(:)
real(8),allocatable :: Bmat(:,:), rhs(:)

nuvo = nvir_a*nocc_a + nvir_b*nocc_b

if (direct_mode) then
   bchunk = npert
else
   bchunk = int(max(1_8, cache_budget_bytes / (7_8*int(nConts,8)*int(nConts,8)*8_8)))
   bchunk = max(1, min(npert, bchunk))
endif

ndiis = int(max(1_8, diis_budget_bytes / (2_8*int(nuvo,8)*int(npert,8)*8_8)))
ndiis = max(2, min(10, ndiis))

allocate(Dtr(nConts,nConts,bchunk,2), Dtot(nConts,nConts,bchunk))
allocate(Vtr(nConts,nConts,bchunk,2), Kscr(nConts,nConts,bchunk,2))
allocate(Tmpa(nConts,nocc_a), Tmpb(nConts,nocc_b))
allocate(Unew_a(nvir_a,nocc_a,npert), Unew_b(nvir_b,nocc_b,npert))
allocate(converged(npert), actlist(npert), nhist(npert))
allocate(dsU(nuvo,ndiis,npert), dsE(nuvo,ndiis,npert))
allocate(flat_u(nuvo), flat_unew(nuvo), flat_err(nuvo))

Uvo_a = 0.0d0
Uvo_b = 0.0d0
converged = .false.
nhist = 0
niter_total = 0
allocate(thr_hit(NTHR,npert))
thr_hit = .false.
thr_passes = 0

conv_thr = 5.0d-5
block
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_CPHF_CONV", envval, status=envstat)
   if (envstat .eq. 0) read(envval,*) conv_thr
end block

do iter = 1,200
   nact = 0
   do p = 1,npert
      if (.not. converged(p)) then
         nact = nact + 1
         actlist(nact) = p
      endif
   enddo
   if (nact .eq. 0) exit
   niter_total = niter_total + nact

   cstart = 1
   do while (cstart .le. nact)
      cend = min(nact, cstart+bchunk-1)
      nb = cend - cstart + 1
      do cidx = 1,nb
         p = actlist(cstart+cidx-1)
         Dtr(:,:,cidx,1) = matmul(C_a(:,nocc_a+1:nConts), matmul(Uvo_a(:,:,p), transpose(C_a(:,1:nocc_a)))) &
                         + matmul(C_a(:,1:nocc_a), matmul(transpose(Uvo_a(:,:,p)), transpose(C_a(:,nocc_a+1:nConts))))
         Dtr(:,:,cidx,2) = matmul(C_b(:,nocc_b+1:nConts), matmul(Uvo_b(:,:,p), transpose(C_b(:,1:nocc_b)))) &
                         + matmul(C_b(:,1:nocc_b), matmul(transpose(Uvo_b(:,:,p)), transpose(C_b(:,nocc_b+1:nConts))))
         Dtot(:,:,cidx) = Dtr(:,:,cidx,1) + Dtr(:,:,cidx,2)
      enddo
      call cpscf_vind_spin(nConts, nb, 2, HF_exchange_frac, Dtr(:,:,1:nb,:), &
                           Dtot(:,:,1:nb), Vtr(:,:,1:nb,:), Kscr(:,:,1:nb,:))
      do cidx = 1,nb
         p = actlist(cstart+cidx-1)
         Tmpa = matmul(Vtr(:,:,cidx,1), C_a(:,1:nocc_a))
         Tmpb = matmul(Vtr(:,:,cidx,2), C_b(:,1:nocc_b))
         do i = 1,nocc_a
            do a = 1,nvir_a
               Unew_a(a,i,p) = ( RHSfix_a(a,i,p) + dot_product(C_a(:,nocc_a+a), Tmpa(:,i)) ) &
                               / (eLev_a(i)-eLev_a(nocc_a+a))
            enddo
         enddo
         do i = 1,nocc_b
            do a = 1,nvir_b
               Unew_b(a,i,p) = ( RHSfix_b(a,i,p) + dot_product(C_b(:,nocc_b+a), Tmpb(:,i)) ) &
                               / (eLev_b(i)-eLev_b(nocc_b+a))
            enddo
         enddo
      enddo
      cstart = cend + 1
   enddo

   do cidx = 1,nact
      p = actlist(cidx)
      dmax_p = max( maxval(abs(Unew_a(:,:,p)-Uvo_a(:,:,p))), &
                    maxval(abs(Unew_b(:,:,p)-Uvo_b(:,:,p))) )

      flat_unew(1:nvir_a*nocc_a) = reshape(Unew_a(:,:,p), (/nvir_a*nocc_a/))
      flat_u(1:nvir_a*nocc_a)    = reshape(Uvo_a(:,:,p),  (/nvir_a*nocc_a/))
      off = nvir_a*nocc_a
      flat_unew(off+1:nuvo) = reshape(Unew_b(:,:,p), (/nvir_b*nocc_b/))
      flat_u(off+1:nuvo)    = reshape(Uvo_b(:,:,p),  (/nvir_b*nocc_b/))
      flat_err = flat_unew - flat_u

      ihist = mod(nhist(p), ndiis) + 1
      dsU(:,ihist,p) = flat_unew
      dsE(:,ihist,p) = flat_err
      nhist(p) = nhist(p) + 1
      Uvo_a(:,:,p) = Unew_a(:,:,p)
      Uvo_b(:,:,p) = Unew_b(:,:,p)

      nuse = min(nhist(p), ndiis)
      if (nuse .ge. 2) then
         allocate(Bmat(nuse+1,nuse+1), rhs(nuse+1), ipiv(nuse+1))
         do r = 1,nuse
            do s = 1,nuse
               Bmat(r,s) = dot_product(dsE(:,r,p), dsE(:,s,p))
            enddo
         enddo
         Bmat(1:nuse,nuse+1) = -1.0d0
         Bmat(nuse+1,1:nuse) = -1.0d0
         Bmat(nuse+1,nuse+1) = 0.0d0
         rhs = 0.0d0
         rhs(nuse+1) = -1.0d0
         call dgesv(nuse+1, 1, Bmat, nuse+1, ipiv, rhs, nuse+1, info)
         if (info .eq. 0) then
            flat_u = 0.0d0
            do r = 1,nuse
               flat_u = flat_u + rhs(r)*dsU(:,r,p)
            enddo
            Uvo_a(:,:,p) = reshape(flat_u(1:nvir_a*nocc_a), (/nvir_a,nocc_a/))
            Uvo_b(:,:,p) = reshape(flat_u(off+1:nuvo), (/nvir_b,nocc_b/))
         endif
         deallocate(Bmat, rhs, ipiv)
      endif

      do it = 1,NTHR
         if (.not. thr_hit(it,p)) then
            thr_passes(it) = thr_passes(it) + 1
            if (dmax_p < thr_list(it)) thr_hit(it,p) = .true.
         endif
      enddo

      if (dmax_p < conv_thr) converged(p) = .true.
   enddo
enddo

block
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_CPHF_THRESHOLD_SCAN", envval, status=envstat)
   if (envstat .eq. 0) then
      print *, ' CPHF passes vs convergence threshold (total over', npert, 'perturbations):'
      do it = 1,NTHR
         print *, '   thr=', thr_list(it), ' passes=', thr_passes(it)
      enddo
   endif
end block
deallocate(thr_hit)

deallocate(Dtr,Dtot,Vtr,Kscr,Tmpa,Tmpb,Unew_a,Unew_b,converged,actlist,nhist,dsU,dsE)
deallocate(flat_u,flat_unew,flat_err)
end subroutine cpscf_solve_vo_direct_spin

subroutine cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, &
                          U_full, h1_MO_all, s1_MO_all, mo_e1_all)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force, integrals_compute_force_1e, integrals_compute_force_exchange_lr
use mod_exchange, only: exchange_build, RS_omega, cosx_enabled
use mod_integrals, only: direct_mode
use mod_profile, only: prof_start, prof_stop, prof_add, prof_add_cpu
use mod_xc_hess, only: xc_fxc_project, xc_vxc_deriv_occ, xc_fxc_c, xc_fxc_set_loc, xc_fxc_clear_loc
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: U_full(nConts,nocc,3,natoms)
real(8),intent(out) :: h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms)
real(8),intent(out) :: mo_e1_all(nocc,nocc,3,natoms)

integer :: nvir
integer :: npert, ipert
integer :: pchunk, pstart, pend, pn, pidx
integer(8), parameter :: pchunk_cache_bytes = 2_8*1024_8*1024_8
real(8),allocatable :: s1oo_all(:,:,:)
real(8),allocatable :: Dbat(:,:,:), Jbat(:,:,:), Kabat(:,:,:), VCb(:,:,:), Lb(:,:,:)
logical :: vcm1
real(8) :: tstage
real(8) :: omp_get_wtime
external :: omp_get_wtime
real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: RHSfix(:,:,:), Uvo_all(:,:,:)
real(8),allocatable :: dF_full(:,:), dS_full(:,:), h1_MO(:,:), s1_MO(:,:)
real(8),allocatable :: s1oo(:,:), Pconn(:,:), Kbc(:,:), Tmp(:,:)
real(8),allocatable :: Uvo(:,:), Pvo(:,:), Ptot2(:,:), v1oo(:,:)
real(8),allocatable :: Din1(:,:,:), Jout1(:,:,:)
real(8),allocatable :: Xd_oo(:,:,:), Xd_vo(:,:,:), Xp_vo(:,:,:), Xo1(:,:,:), Xo2(:,:,:)
logical :: xc_proj_mo
integer :: iat,ic,i,j,a,b,mu

nvir = nConts - nocc

cpscf_xc_active = (HF_exchange_frac < 1.0d0)
if (cpscf_xc_active) then
   if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
   allocate(cpscf_xc_Ptot(nConts,nConts))
   cpscf_xc_Ptot = Pa + Pb
endif
if (cpscf_xc_active) then
   block
      use mod_xc_hess, only: xc_resp_grid_build
      use GRID_info, only: XCGRID_TIER_LO
      character(8) :: e
      call get_environment_variable("ENGINE_CPSCF_XC_GRID", e)
      if (trim(e) .ne. "full") call xc_resp_grid_build(XCGRID_TIER_LO)
   end block
endif

cpscf_ecp_atom = -1
allocate(dS(nConts,nConts,3), dHcore(nConts,nConts,3))
allocate(dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3))
tstage = omp_get_wtime()
call cpu_time(tstage_c)
call prof_start("hess_force_int")
if (allocated(cpscf_dfj_jx) .and. (allocated(cpscf_cosx_kxa) .or. HF_exchange_frac == 0.0d0)) then
   call integrals_compute_force_1e(nConts, dS, dHcore)
   dJi = 0.0d0 ; dKa = 0.0d0 ; dKb = 0.0d0
else
   call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
endif
call prof_stop("hess_force_int")
if (RS_omega > 0.0d0 .and. .not. cosx_enabled) then
   if (allocated(cpscf_lr_dKa)) deallocate(cpscf_lr_dKa, cpscf_lr_dKb)
   allocate(cpscf_lr_dKa(nConts,nConts,3), cpscf_lr_dKb(nConts,nConts,3))
   call integrals_compute_force_exchange_lr(nConts, Pa, Pb, cpscf_lr_dKa, cpscf_lr_dKb)
endif
cphf_t_forceint = omp_get_wtime() - tstage

allocate(ao_atom(nConts))
mu = 0
do iat = 1,natoms
   do i = 1,atoms(iat)%nconts
      mu = mu+1
      ao_atom(mu) = iat
   enddo
enddo

cphf_t_amat = 0.0d0

allocate(dF_full(nConts,nConts), dS_full(nConts,nConts))
allocate(h1_MO(nConts,nConts), s1_MO(nConts,nConts))
allocate(s1oo(nocc,nocc), Pconn(nConts,nConts))
allocate(Kbc(nConts,nConts), Tmp(nConts,nocc))
allocate(Uvo(nvir,nocc), Pvo(nConts,nConts), Ptot2(nConts,nConts), v1oo(nocc,nocc))
allocate(Din1(nConts,nConts,1), Jout1(nConts,nConts,1))

U_full = 0.0d0
cphf_t_fockderiv = 0.0d0 ; cphf_c_fockderiv = 0.0d0
cphf_t_iter = 0.0d0
cphf_t_pass1 = 0.0d0
cphf_t_pass3 = 0.0d0
cphf_t_mo = 0.0d0
cphf_n_iter_total = 0

npert = 3*natoms
allocate(RHSfix(nvir,nocc,npert), Uvo_all(nvir,nocc,npert))
allocate(s1oo_all(nocc,nocc,npert))
block
   character(16) :: solv0
   call get_environment_variable("ENGINE_CPHF_SOLVER", solv0)
   cphf_voo_track = (trim(solv0) .ne. "diis") .and. .not. cpscf_xc_active
end block
vcm1 = cpscf_use_vind_c(HF_exchange_frac)
cphf_voo_ready = .false.
if (allocated(cphf_voo_conn)) deallocate(cphf_voo_conn)
if (allocated(cphf_voo_u)) deallocate(cphf_voo_u)
if (cphf_voo_track) allocate(cphf_voo_conn(nocc,nocc,npert), cphf_voo_u(nocc,nocc,npert))
if (direct_mode) then
   pchunk = npert
else
   pchunk = int(max(1_8, pchunk_cache_bytes / (2_8*int(nConts,8)*int(nConts,8)*8_8)))
   pchunk = max(1, min(npert, pchunk))
endif
if (vcm1) then
   allocate(Dbat(1,1,1), Jbat(1,1,1), Kabat(1,1,1))
else
   allocate(Dbat(nConts,nConts,pchunk), Jbat(nConts,nConts,pchunk))
   allocate(Kabat(nConts,nConts,pchunk))
endif
allocate(Xd_oo(nocc,nocc,pchunk), Xd_vo(nvir,nocc,pchunk), Xp_vo(nvir,nocc,pchunk))
allocate(Xo1(nocc,nocc,pchunk), Xo2(nocc,nocc,pchunk))

block
   character(len=8) :: e1, e2
   call get_environment_variable("ENGINE_XC_PROJ_GRP", e1)
   call get_environment_variable("ENGINE_XC_PROJ_CHK", e2)
   xc_proj_mo = (trim(e1) .eq. "1" .or. trim(e2) .eq. "1")
   call get_environment_variable("ENGINE_XC_PROJ_LOC", e1)
   if (HF_exchange_frac < 1.0d0 .and. trim(e1) .eq. "1") call xc_fxc_set_loc(nConts, nocc, C_a)
end block
ipert = 0
if (HF_exchange_frac < 1.0d0) then
   block
      character(len=8) :: e
      call get_environment_variable("ENGINE_XC_VXCDERIV_FULL", e)
      if (trim(e) .ne. "1") then
         allocate(cpscf_xc_vxco(nConts,nocc,3,natoms))
         tstage = omp_get_wtime()
         call cpu_time(tstage_c)
         call prof_start("fd_xc_vxc_deriv")
         call xc_vxc_deriv_occ(nConts, natoms, nocc, Pa+Pb, C_a(:,1:nocc), cpscf_xc_vxco)
         call prof_stop("fd_xc_vxc_deriv")
         cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
         call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
      endif
   end block
endif
do iat = 1,natoms
   do ic = 1,3
      ipert = ipert + 1
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, dF_full, dS_full)
      cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
      call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      block
         real(8) :: FC(nConts,nocc)
         call dgemm('N','N',nConts,nocc,nConts,1.0d0,dF_full,nConts,C_a,nConts,0.0d0,FC,nConts)
         if (cpscf_cosx_kx_proj) FC = FC - HF_exchange_frac*cpscf_cosx_kxa(:,1:nocc,ic,iat)
         if (allocated(cpscf_xc_vxco)) FC = FC + cpscf_xc_vxco(:,:,ic,iat)
         call dgemm('T','N',nConts,nocc,nConts,1.0d0,C_a,nConts,FC,nConts,0.0d0,h1_MO_all(1,1,ic,iat),nConts)
         call dgemm('N','N',nConts,nocc,nConts,1.0d0,dS_full,nConts,C_a,nConts,0.0d0,FC,nConts)
         call dgemm('T','N',nConts,nocc,nConts,1.0d0,C_a,nConts,FC,nConts,0.0d0,s1_MO_all(1,1,ic,iat),nConts)
      end block
      s1_MO(:,1:nocc) = s1_MO_all(:,:,ic,iat)
      cphf_t_mo = cphf_t_mo + (omp_get_wtime() - tstage)

      s1oo_all(:,:,ipert) = s1_MO(1:nocc,1:nocc)
   enddo
enddo

tstage = omp_get_wtime()
call cpu_time(tstage_c)
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      if (.not. vcm1) Dbat(:,:,pidx) = -2.0d0*matmul(C_a(:,1:nocc), &
                        matmul(s1oo_all(:,:,ipert), transpose(C_a(:,1:nocc))))
      Xd_oo(:,:,pidx) = -s1oo_all(:,:,ipert)
   enddo
   if (vcm1) then
      allocate(VCb(nConts,nocc,pn), Lb(nConts,nocc,pn))
      do pidx = 1,pn
         Lb(:,:,pidx) = -matmul(C_a(:,1:nocc), s1oo_all(:,:,pstart+pidx-1))
      enddo
      call cpscf_vind_c(nConts, pn, HF_exchange_frac, nocc, C_a(:,1:nocc), Lb, VCb)
      deallocate(Lb)
   else
      call cpscf_vind(nConts, pn, HF_exchange_frac, Dbat(:,:,1:pn), Jbat(:,:,1:pn), &
                      Kabat(:,:,1:pn), skip_xc=.true.)
   endif
   if (cpscf_xc_active) then
      call xc_fxc_project(nConts, nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                          nvir, nocc, C_a(:,nocc+1:nConts), C_a(:,1:nocc), &
                          pn, Xd_oo(:,:,1:pn), cpscf_xc_Ptot, Xp_vo(:,:,1:pn))
   endif
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO(:,1:nocc) = h1_MO_all(:,:,ic,iat)
      s1_MO(:,1:nocc) = s1_MO_all(:,:,ic,iat)
      if (vcm1) then
         Tmp = VCb(:,:,pidx)
      else
         Tmp = matmul(Jbat(:,:,pidx), C_a(:,1:nocc))
      endif
      if (cphf_voo_track) cphf_voo_conn(:,:,ipert) = matmul(transpose(C_a(:,1:nocc)), Tmp)

      do i = 1,nocc
         do a = 1,nvir
            RHSfix(a,i,ipert) = h1_MO(nocc+a,i) - s1_MO(nocc+a,i)*eLev_a(i) &
                                + dot_product(C_a(:,nocc+a), Tmp(:,i))
            if (cpscf_xc_active) RHSfix(a,i,ipert) = RHSfix(a,i,ipert) + Xp_vo(a,i,pidx)
         enddo
      enddo
   enddo
   if (allocated(VCb)) deallocate(VCb)
   pstart = pend + 1
enddo
cphf_t_pass1 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass1 = ctmp_c - tstage_c
call prof_add("hess_cphf_rhs", cphf_t_pass1)

call prof_add_cpu("hess_cphf_rhs", cphf_c_pass1)
tstage = omp_get_wtime()
call cpu_time(tstage_c)
call cphf_solve_vo_direct(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                          RHSfix, Uvo_all, cphf_n_iter_total)
cphf_t_iter = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_iter = ctmp_c - tstage_c
call prof_add("hess_cphf_solve", cphf_t_iter)

call prof_add_cpu("hess_cphf_solve", cphf_c_iter)
tstage = omp_get_wtime()
call cpu_time(tstage_c)
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      Uvo = Uvo_all(:,:,ipert)
      s1oo = s1oo_all(:,:,ipert)
      Pconn = -2.0d0*matmul(C_a(:,1:nocc), matmul(s1oo, transpose(C_a(:,1:nocc))))
      Pvo = 2.0d0*( matmul(C_a(:,nocc+1:nConts), matmul(Uvo, transpose(C_a(:,1:nocc)))) &
                  + matmul(C_a(:,1:nocc), matmul(transpose(Uvo), transpose(C_a(:,nocc+1:nConts)))) )
      if (.not. vcm1) Dbat(:,:,pidx) = Pconn + Pvo
      Xd_oo(:,:,pidx) = -s1oo
      Xd_vo(:,:,pidx) = 2.0d0*Uvo

      do i = 1,nocc
         do j = 1,nocc
            U_full(j,i,ic,iat) = -0.5d0*s1oo(j,i)
         enddo
         do a = 1,nvir
            U_full(nocc+a,i,ic,iat) = Uvo(a,i)
         enddo
      enddo
   enddo
   if (.not. cphf_voo_ready) then
      if (vcm1) then
         allocate(VCb(nConts,nocc,pn), Lb(nConts,nocc,pn))
         do pidx = 1,pn
            ipert = pstart + pidx - 1
            Lb(:,:,pidx) = -matmul(C_a(:,1:nocc), s1oo_all(:,:,ipert)) &
                           + 2.0d0*matmul(C_a(:,nocc+1:nConts), Uvo_all(:,:,ipert))
         enddo
         call cpscf_vind_c(nConts, pn, HF_exchange_frac, nocc, C_a(:,1:nocc), Lb, VCb)
         deallocate(Lb)
      else
         call cpscf_vind(nConts, pn, HF_exchange_frac, Dbat(:,:,1:pn), Jbat(:,:,1:pn), &
                         Kabat(:,:,1:pn), skip_xc=.true.)
      endif
   endif
   if (cpscf_xc_active) then
      if (xc_proj_mo) then
         call xc_fxc_project(nConts, nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                             nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                             pn, Xd_oo(:,:,1:pn), cpscf_xc_Ptot, Xo1(:,:,1:pn))
         call xc_fxc_project(nConts, nvir, nocc, C_a(:,nocc+1:nConts), C_a(:,1:nocc), &
                             nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                             pn, Xd_vo(:,:,1:pn), cpscf_xc_Ptot, Xo2(:,:,1:pn))
      else
         block
            real(8),allocatable :: Lx(:,:,:), Xt(:,:,:)
            allocate(Lx(nConts,nocc,pn), Xt(nConts,nocc,pn))
            do pidx = 1,pn
               call dgemm('N','N',nConts,nocc,nocc,1.0d0,C_a,nConts,Xd_oo(1,1,pidx),nocc,0.0d0,Lx(1,1,pidx),nConts)
               call dgemm('N','N',nConts,nocc,nvir,1.0d0,C_a(1,nocc+1),nConts,Xd_vo(1,1,pidx),nvir,1.0d0,Lx(1,1,pidx),nConts)
            enddo
            call xc_fxc_c(nConts, nocc, C_a(:,1:nocc), pn, Lx, cpscf_xc_Ptot, Xt)
            do pidx = 1,pn
               call dgemm('T','N',nocc,nocc,nConts,1.0d0,C_a,nConts,Xt(1,1,pidx),nConts,0.0d0,Xo1(1,1,pidx),nocc)
            enddo
            Xo2(:,:,1:pn) = 0.0d0
            deallocate(Lx, Xt)
         end block
      endif
   endif
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO(:,1:nocc) = h1_MO_all(:,:,ic,iat)
      s1_MO(:,1:nocc) = s1_MO_all(:,:,ic,iat)
      if (cphf_voo_ready) then
         v1oo = cphf_voo_conn(:,:,ipert) + cphf_voo_u(:,:,ipert)
      else if (vcm1) then
         v1oo = matmul(transpose(C_a(:,1:nocc)), VCb(:,:,pidx))
      else
         v1oo = matmul(transpose(C_a(:,1:nocc)), matmul(Jbat(:,:,pidx), C_a(:,1:nocc)))
      endif
      if (cpscf_xc_active) v1oo = v1oo + Xo1(:,:,pidx) + Xo2(:,:,pidx)
      do i = 1,nocc
         do j = 1,nocc
            mo_e1_all(i,j,ic,iat) = h1_MO(i,j) - 0.5d0*(eLev_a(i)+eLev_a(j))*s1_MO(i,j) + v1oo(i,j)
         enddo
      enddo
   enddo
   if (allocated(VCb)) deallocate(VCb)
   pstart = pend + 1
enddo
cphf_t_pass3 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass3 = ctmp_c - tstage_c
call prof_add("hess_cphf_obs", cphf_t_pass3)
call prof_add_cpu("hess_cphf_obs", cphf_c_pass3)
block
   use mod_xc_hess, only: xc_resp_grid_clear
   call xc_resp_grid_clear()
end block

call prof_add("hess_fock_deriv", cphf_t_fockderiv)

call prof_add_cpu("hess_fock_deriv", cphf_c_fockderiv)
deallocate(dS,dHcore,dJi,dKa,dKb,ao_atom,RHSfix,Uvo_all,s1oo_all)
if (allocated(cpscf_lr_dKa)) deallocate(cpscf_lr_dKa, cpscf_lr_dKb)
deallocate(Dbat,Jbat,Kabat)
deallocate(Xd_oo,Xd_vo,Xp_vo,Xo1,Xo2)
deallocate(dF_full,dS_full,h1_MO,s1_MO,s1oo,Pconn,Kbc,Tmp)
if (allocated(cpscf_xc_vxco)) deallocate(cpscf_xc_vxco)
call xc_fxc_clear_loc()
deallocate(Uvo,Pvo,Ptot2,v1oo)
deallocate(Din1,Jout1)
if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
cpscf_xc_active = .false.
end subroutine cphf_solve_rhf

subroutine cpscf_vind_spin_c(nConts, nb, HF_exchange_frac, noa, nob, Ca, Cb, La, Lb, VCa, VCb)
use mod_exchange, only: RS_beta
use mod_integrals, only: integrals_build_exchange_lr_batch, &
   integrals_build_coulomb_df_batch_c, integrals_build_exchange_cosx_batch_c, &
                         integrals_build_jk_batch_uhf, integrals_build_coulomb_batch, integrals_build_exchange_batch, &
                         integrals_build_coulomb_df_batch
use mod_exchange, only: cosx_enabled
use mod_xc_hess, only: xc_fxc_c_u
use omp_lib, only: omp_get_wtime
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_vind_spin_c
implicit none
integer,intent(in) :: nConts, nb, noa, nob
real(8),intent(in) :: HF_exchange_frac, Ca(nConts,noa), Cb(nConts,nob)
real(8),intent(in) :: La(nConts,noa,nb), Lb(nConts,nob,nb)
real(8),intent(out) :: VCa(nConts,noa,nb), VCb(nConts,nob,nb)
real(8),allocatable :: Lc(:,:,:), Cc(:,:), JC(:,:,:), KC(:,:,:), Xa(:,:,:), Xb(:,:,:)
real(8),allocatable :: Da(:,:,:), Db(:,:,:), Dt(:,:,:), Jt(:,:,:), Ka(:,:,:), Kb(:,:,:), Kla(:,:,:), Klb(:,:,:)
integer :: b, c0, nch, ib2
logical :: rij, kcosx
real(8) :: tv0
rij = cpscf_dfj_response()
kcosx = cosx_enabled .and. HF_exchange_frac /= 0.0d0
if (rij .and. (kcosx .or. HF_exchange_frac == 0.0d0)) then
tv0 = omp_get_wtime()
allocate(Lc(nConts,noa+nob,nb), Cc(nConts,noa+nob), JC(nConts,noa+nob,nb))
Cc(:,1:noa) = Ca ; Cc(:,noa+1:noa+nob) = Cb
do b = 1,nb
   Lc(:,1:noa,b) = La(:,:,b) ; Lc(:,noa+1:noa+nob,b) = Lb(:,:,b)
enddo
call integrals_build_coulomb_df_batch_c(nConts, nb, Lc, noa+nob, Cc, JC)
do b = 1,nb
   VCa(:,:,b) = JC(:,1:noa,b) ; VCb(:,:,b) = JC(:,noa+1:noa+nob,b)
enddo
deallocate(Lc, Cc, JC)
cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
if (HF_exchange_frac /= 0.0d0) then
   tv0 = omp_get_wtime()
   allocate(KC(nConts,max(noa,nob),nb))
   call integrals_build_exchange_cosx_batch_c(nConts, nb, La, noa, Ca, KC(:,1:noa,:))
   VCa = VCa - HF_exchange_frac*KC(:,1:noa,:)
   call integrals_build_exchange_cosx_batch_c(nConts, nb, Lb, nob, Cb, KC(:,1:nob,:))
   VCb = VCb - HF_exchange_frac*KC(:,1:nob,:)
   deallocate(KC)
   cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   cphf_n_vk = cphf_n_vk + 2
   cphf_nd_vk = cphf_nd_vk + 2*nb
endif
else
   tv0 = omp_get_wtime()
   allocate(Da(nConts,nConts,16), Db(nConts,nConts,16), Dt(nConts,nConts,16), Jt(nConts,nConts,16))
   allocate(Ka(nConts,nConts,16), Kb(nConts,nConts,16), Kla(nConts,nConts,16), Klb(nConts,nConts,16))
   do c0 = 0, nb-1, 16
      nch = min(16, nb-c0)
      do ib2 = 1,nch
         call dgemm('N','T',nConts,nConts,noa,1.0d0,La(1,1,c0+ib2),nConts,Ca,nConts,0.0d0,Da(1,1,ib2),nConts)
         Da(:,:,ib2) = Da(:,:,ib2) + transpose(Da(:,:,ib2))
         call dgemm('N','T',nConts,nConts,nob,1.0d0,Lb(1,1,c0+ib2),nConts,Cb,nConts,0.0d0,Db(1,1,ib2),nConts)
         Db(:,:,ib2) = Db(:,:,ib2) + transpose(Db(:,:,ib2))
         Dt(:,:,ib2) = Da(:,:,ib2) + Db(:,:,ib2)
      enddo
      if (.not. rij .and. .not. kcosx) then
         call integrals_build_jk_batch_uhf(nConts, nch, Dt(:,:,1:nch), Da(:,:,1:nch), Db(:,:,1:nch), &
                                           Jt(:,:,1:nch), Ka(:,:,1:nch), Kb(:,:,1:nch))
      else
         if (rij) then
            call integrals_build_coulomb_df_batch(nConts, nch, Dt(:,:,1:nch), Jt(:,:,1:nch))
         else
            call integrals_build_coulomb_batch(nConts, nch, Dt(:,:,1:nch), Jt(:,:,1:nch))
         endif
         Ka(:,:,1:nch) = 0.0d0 ; Kb(:,:,1:nch) = 0.0d0
         if (HF_exchange_frac /= 0.0d0 .and. .not. kcosx) &
            call integrals_build_exchange_batch(nConts, nch, Da(:,:,1:nch), Db(:,:,1:nch), Ka(:,:,1:nch), Kb(:,:,1:nch))
      endif
      if (allocated(cpscf_lr_dKa)) &
         call integrals_build_exchange_lr_batch(nConts, nch, Da(:,:,1:nch), Db(:,:,1:nch), Kla(:,:,1:nch), Klb(:,:,1:nch))
      do ib2 = 1,nch
         Jt(:,:,ib2) = Jt(:,:,ib2) - HF_exchange_frac*Ka(:,:,ib2)
         if (allocated(cpscf_lr_dKa)) Jt(:,:,ib2) = Jt(:,:,ib2) - RS_beta*Kla(:,:,ib2)
         call dgemm('N','N',nConts,noa,nConts,1.0d0,Jt(1,1,ib2),nConts,Ca,nConts,0.0d0,VCa(1,1,c0+ib2),nConts)
         Jt(:,:,ib2) = Jt(:,:,ib2) + HF_exchange_frac*(Ka(:,:,ib2) - Kb(:,:,ib2))
         if (allocated(cpscf_lr_dKa)) Jt(:,:,ib2) = Jt(:,:,ib2) + RS_beta*(Kla(:,:,ib2) - Klb(:,:,ib2))
         call dgemm('N','N',nConts,nob,nConts,1.0d0,Jt(1,1,ib2),nConts,Cb,nConts,0.0d0,VCb(1,1,c0+ib2),nConts)
      enddo
   enddo
   deallocate(Da, Db, Dt, Jt, Ka, Kb, Kla, Klb)
   cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
   if (kcosx) then
      tv0 = omp_get_wtime()
      allocate(KC(nConts,max(noa,nob),nb))
      call integrals_build_exchange_cosx_batch_c(nConts, nb, La, noa, Ca, KC(:,1:noa,:))
      VCa = VCa - HF_exchange_frac*KC(:,1:noa,:)
      call integrals_build_exchange_cosx_batch_c(nConts, nb, Lb, nob, Cb, KC(:,1:nob,:))
      VCb = VCb - HF_exchange_frac*KC(:,1:nob,:)
      deallocate(KC)
      cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   endif
endif
if (cpscf_xc_uks) then
   tv0 = omp_get_wtime()
   allocate(Xa(nConts,noa,nb), Xb(nConts,nob,nb))
   call xc_fxc_c_u(nConts, noa, nob, Ca, Cb, nb, La, Lb, cpscf_xc_Pa, cpscf_xc_Pb, Xa, Xb)
   VCa = VCa + Xa ; VCb = VCb + Xb
   deallocate(Xa, Xb)
   cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
endif
if (cosmo_hess_active) call cosmo_hess_vind_spin_c(nConts, nb, noa, nob, Ca, Cb, La, Lb, VCa, VCb)
end subroutine cpscf_vind_spin_c

subroutine cpscf_solve_vo_pcg_spin(nConts, noa, nva, nob, nvb, npert, Ca, ea, Cb, eb, &
                                   HF_exchange_frac, RHSa, RHSb, Ua, Ub, niter_total, conv_thr)
implicit none
integer,intent(in) :: nConts, noa, nva, nob, nvb, npert
real(8),intent(in) :: Ca(nConts,nConts), ea(nConts), Cb(nConts,nConts), eb(nConts)
real(8),intent(in) :: HF_exchange_frac, conv_thr
real(8),intent(in) :: RHSa(nva,noa,npert), RHSb(nvb,nob,npert)
real(8),intent(out) :: Ua(nva,noa,npert), Ub(nvb,nob,npert)
integer,intent(out) :: niter_total
real(8),allocatable :: ma(:,:), mb(:,:), da(:,:), db(:,:), Ra(:,:,:), Rb(:,:,:), Za(:,:,:), Zb(:,:,:)
real(8),allocatable :: Pa_(:,:,:), Pb_(:,:,:), Qa(:,:,:), Qb(:,:,:), rz(:)
logical,allocatable :: conv(:)
integer,allocatable :: act(:)
integer :: p, it, nact, k, a, i
real(8) :: alpha, beta, rznew, pq
logical :: verb
character(len=8) :: ev
call get_environment_variable("ENGINE_CPHF_VERBOSE", ev)
verb = (trim(ev) .eq. "1")
allocate(da(nva,noa), db(nvb,nob), rz(npert), conv(npert), act(npert))
allocate(Ra(nva,noa,npert), Za(nva,noa,npert), Pa_(nva,noa,npert), Qa(nva,noa,npert))
allocate(Rb(nvb,nob,npert), Zb(nvb,nob,npert), Pb_(nvb,nob,npert), Qb(nvb,nob,npert))
do i = 1,noa
   do a = 1,nva
      da(a,i) = ea(noa+a) - ea(i)
   enddo
enddo
do i = 1,nob
   do a = 1,nvb
      db(a,i) = eb(nob+a) - eb(i)
   enddo
enddo
if (verb) print '(A,2ES11.3)', '  UHF PCG min orbital-energy gaps a / b:', minval(da), minval(db)
block
   real(8) :: pfl
   character(len=16) :: e
   pfl = 0.0d0
   call get_environment_variable("ENGINE_CPHF_PFLOOR", e)
   if (len_trim(e) .gt. 0) read(e,*) pfl
   allocate(ma(nva,noa), mb(nvb,nob))
   ma = max(da, pfl) ; mb = max(db, pfl)
end block
block
   real(8) :: gthr
   character(len=16) :: e
   integer,parameter :: NMAXD = 64
   integer :: nd, ks, ka_(NMAXD), ki_(NMAXD), kspin(NMAXD), q
   real(8),allocatable :: La1(:,:,:), Lb1(:,:,:), V1a(:,:,:), V1b(:,:,:)
   gthr = 0.05d0
   call get_environment_variable("ENGINE_CPHF_DIAGGAP", e)
   if (len_trim(e) .gt. 0) read(e,*) gthr
   nd = 0
   do i = 1,noa
      do a = 1,nva
         if (da(a,i) .lt. gthr .and. nd .lt. NMAXD) then
            nd = nd + 1 ; ka_(nd) = a ; ki_(nd) = i ; kspin(nd) = 1
         endif
      enddo
   enddo
   do i = 1,nob
      do a = 1,nvb
         if (db(a,i) .lt. gthr .and. nd .lt. NMAXD) then
            nd = nd + 1 ; ka_(nd) = a ; ki_(nd) = i ; kspin(nd) = 2
         endif
      enddo
   enddo
   if (nd .gt. 0) then
      allocate(La1(nConts,noa,nd), Lb1(nConts,nob,nd), V1a(nConts,noa,nd), V1b(nConts,nob,nd))
      La1 = 0.0d0 ; Lb1 = 0.0d0
      do q = 1,nd
         if (kspin(q) .eq. 1) then
            La1(:,ki_(q),q) = Ca(:,noa+ka_(q))
         else
            Lb1(:,ki_(q),q) = Cb(:,nob+ka_(q))
         endif
      enddo
      call cpscf_vind_spin_c(nConts, nd, HF_exchange_frac, noa, nob, Ca(:,1:noa), Cb(:,1:nob), La1, Lb1, V1a, V1b)
      do q = 1,nd
         if (kspin(q) .eq. 1) then
            ks = ka_(q)
            ma(ks,ki_(q)) = max(ma(ks,ki_(q)), da(ks,ki_(q)) + sum(Ca(:,noa+ks)*V1a(:,ki_(q),q)))
         else
            ks = ka_(q)
            mb(ks,ki_(q)) = max(mb(ks,ki_(q)), db(ks,ki_(q)) + sum(Cb(:,nob+ks)*V1b(:,ki_(q),q)))
         endif
      enddo
      if (verb) print '(A,I0,A)', '  UHF PCG exact diagonal for ', nd, ' small-gap pairs'
      deallocate(La1, Lb1, V1a, V1b)
   endif
end block
do p = 1,npert
   Ua(:,:,p) = -RHSa(:,:,p)/ma
   Ub(:,:,p) = -RHSb(:,:,p)/mb
   act(p) = p
enddo
call apply_A_spin(npert, act, Ua, Ub, Qa, Qb)
niter_total = npert
do p = 1,npert
   Ra(:,:,p) = -RHSa(:,:,p) - Qa(:,:,p) ; Rb(:,:,p) = -RHSb(:,:,p) - Qb(:,:,p)
   Za(:,:,p) = Ra(:,:,p)/ma ; Zb(:,:,p) = Rb(:,:,p)/mb
   Pa_(:,:,p) = Za(:,:,p) ; Pb_(:,:,p) = Zb(:,:,p)
   rz(p) = sum(Ra(:,:,p)*Za(:,:,p)) + sum(Rb(:,:,p)*Zb(:,:,p))
   conv(p) = max(maxval(abs(Za(:,:,p))), maxval(abs(Zb(:,:,p)))) < conv_thr
enddo
do it = 1,200
   nact = 0
   do p = 1,npert
      if (.not. conv(p)) then
         nact = nact + 1
         act(nact) = p
      endif
   enddo
   if (nact .eq. 0) exit
   if (verb) print '(A,I4,A,I5,A,ES10.3)', '  UHF PCG it', it, '  active', nact, '  max|z|', &
      maxval((/ (max(maxval(abs(Za(:,:,act(k)))), maxval(abs(Zb(:,:,act(k))))), k=1,nact) /))
   niter_total = niter_total + nact
   call apply_A_spin(nact, act, Pa_, Pb_, Qa, Qb)
   do k = 1,nact
      p = act(k)
      pq = sum(Pa_(:,:,p)*Qa(:,:,p)) + sum(Pb_(:,:,p)*Qb(:,:,p))
      alpha = rz(p)/pq
      Ua(:,:,p) = Ua(:,:,p) + alpha*Pa_(:,:,p) ; Ub(:,:,p) = Ub(:,:,p) + alpha*Pb_(:,:,p)
      Ra(:,:,p) = Ra(:,:,p) - alpha*Qa(:,:,p) ; Rb(:,:,p) = Rb(:,:,p) - alpha*Qb(:,:,p)
      Za(:,:,p) = Ra(:,:,p)/ma ; Zb(:,:,p) = Rb(:,:,p)/mb
      rznew = sum(Ra(:,:,p)*Za(:,:,p)) + sum(Rb(:,:,p)*Zb(:,:,p))
      beta = rznew/rz(p)
      rz(p) = rznew
      Pa_(:,:,p) = Za(:,:,p) + beta*Pa_(:,:,p) ; Pb_(:,:,p) = Zb(:,:,p) + beta*Pb_(:,:,p)
      conv(p) = max(maxval(abs(Za(:,:,p))), maxval(abs(Zb(:,:,p)))) < conv_thr
   enddo
enddo
deallocate(ma, mb, da, db, Ra, Rb, Za, Zb, Pa_, Pb_, Qa, Qb, rz, conv, act)

contains

subroutine apply_A_spin(n_act, alist, Xa, Xb, Qa_, Qb_)
integer,intent(in) :: n_act, alist(n_act)
real(8),intent(in) :: Xa(nva,noa,npert), Xb(nvb,nob,npert)
real(8),intent(inout) :: Qa_(nva,noa,npert), Qb_(nvb,nob,npert)
real(8),allocatable :: La(:,:,:), Lb(:,:,:), VCa(:,:,:), VCb(:,:,:)
integer :: c, pp
allocate(La(nConts,noa,n_act), Lb(nConts,nob,n_act), VCa(nConts,noa,n_act), VCb(nConts,nob,n_act))
!$omp parallel do private(c,pp)
do c = 1,n_act
   pp = alist(c)
   call dgemm('N','N',nConts,noa,nva,1.0d0,Ca(1,noa+1),nConts,Xa(1,1,pp),nva,0.0d0,La(1,1,c),nConts)
   call dgemm('N','N',nConts,nob,nvb,1.0d0,Cb(1,nob+1),nConts,Xb(1,1,pp),nvb,0.0d0,Lb(1,1,c),nConts)
enddo
!$omp end parallel do
call cpscf_vind_spin_c(nConts, n_act, HF_exchange_frac, noa, nob, Ca(:,1:noa), Cb(:,1:nob), La, Lb, VCa, VCb)
!$omp parallel do private(c,pp)
do c = 1,n_act
   pp = alist(c)
   call dgemm('T','N',nva,noa,nConts,1.0d0,Ca(1,noa+1),nConts,VCa(1,1,c),nConts,0.0d0,Qa_(1,1,pp),nva)
   call dgemm('T','N',nvb,nob,nConts,1.0d0,Cb(1,nob+1),nConts,VCb(1,1,c),nConts,0.0d0,Qb_(1,1,pp),nvb)
   Qa_(:,:,pp) = Qa_(:,:,pp) + da*Xa(:,:,pp)
   Qb_(:,:,pp) = Qb_(:,:,pp) + db*Xb(:,:,pp)
enddo
!$omp end parallel do
deallocate(La, Lb, VCa, VCb)
end subroutine apply_A_spin
end subroutine cpscf_solve_vo_pcg_spin

subroutine cpscf_solve_vo_bcg_spin(nConts, noa, nva, nob, nvb, npert, Ca, ea, Cb, eb, &
                                   HF_exchange_frac, RHSa, RHSb, Ua, Ub, niter_total, conv_thr)
implicit none
integer,intent(in) :: nConts, noa, nva, nob, nvb, npert
real(8),intent(in) :: Ca(nConts,nConts), ea(nConts), Cb(nConts,nConts), eb(nConts)
real(8),intent(in) :: HF_exchange_frac, conv_thr
real(8),intent(in) :: RHSa(nva,noa,npert), RHSb(nvb,nob,npert)
real(8),intent(out) :: Ua(nva,noa,npert), Ub(nvb,nob,npert)
integer,intent(out) :: niter_total
integer :: na, nt, it, k, p, i, a, nk, info, lwork, q
real(8),allocatable :: dm(:), X(:,:), R(:,:), Z(:,:), P_(:,:), Q_(:,:), G(:,:), ev(:), work(:), T(:,:), W(:,:)
logical,allocatable :: conv(:)
logical :: verb
character(len=8) :: e
external :: dgemm, dsyev
call get_environment_variable("ENGINE_CPHF_VERBOSE", e)
verb = (trim(e) .eq. "1")
na = nva*noa ; nt = na + nvb*nob
allocate(dm(nt), X(nt,npert), R(nt,npert), Z(nt,npert), P_(nt,npert), Q_(nt,npert), conv(npert))
do i = 1,noa
   do a = 1,nva
      dm(a + (i-1)*nva) = ea(noa+a) - ea(i)
   enddo
enddo
do i = 1,nob
   do a = 1,nvb
      dm(na + a + (i-1)*nvb) = eb(nob+a) - eb(i)
   enddo
enddo
do p = 1,npert
   R(1:na,p) = -reshape(RHSa(:,:,p),(/na/)) ; R(na+1:nt,p) = -reshape(RHSb(:,:,p),(/nt-na/))
   X(:,p) = R(:,p)/dm
enddo
call applyA(npert, X, Q_)
niter_total = npert
R = R - Q_
nk = npert
do it = 1,200
   do p = 1,npert
      Z(:,p) = R(:,p)/dm
      conv(p) = maxval(abs(Z(:,p))) < conv_thr
      if (conv(p)) Z(:,p) = 0.0d0
   enddo
   if (verb) print '(A,I4,A,I5,A,ES10.3)', '  UHF BCG it', it, '  unconverged', count(.not. conv), &
      '  max|z|', maxval(abs(R/spread(dm,2,npert)))
   if (all(conv)) exit
   if (it .eq. 1) then
      P_(:,1:npert) = Z ; nk = npert
   else
      allocate(T(nk,npert))
      call dgemm('T','N',nk,npert,nt,1.0d0,Q_,nt,Z,nt,0.0d0,T,nk)
      allocate(W(nt,npert))
      W = Z
      call dgemm('N','N',nt,npert,nk,-1.0d0,P_,nt,T,nk,1.0d0,W,nt)
      P_(:,1:npert) = W ; nk = npert
      deallocate(T, W)
   endif
   q = 0
   do p = 1,npert
      if (.not. conv(p)) then
         q = q + 1
         if (q .ne. p) P_(:,q) = P_(:,p)
      endif
   enddo
   nk = q
   call applyA(nk, P_, Q_)
   niter_total = niter_total + nk
   allocate(G(nk,nk), ev(nk), work(1))
   call dgemm('T','N',nk,nk,nt,1.0d0,P_,nt,Q_,nt,0.0d0,G,nk)
   G = 0.5d0*(G + transpose(G))
   call dsyev('V','U',nk,G,nk,ev,work,-1,info)
   lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
   call dsyev('V','U',nk,G,nk,ev,work,lwork,info)
   q = 0
   do k = 1,nk
      if (ev(k) .gt. 1.0d-12*ev(nk)) then
         q = q + 1
         G(:,q) = G(:,k)/sqrt(ev(k))
      endif
   enddo
   allocate(W(nt,q))
   call dgemm('N','N',nt,q,nk,1.0d0,P_,nt,G,nk,0.0d0,W,nt)
   P_(:,1:q) = W
   call dgemm('N','N',nt,q,nk,1.0d0,Q_,nt,G,nk,0.0d0,W,nt)
   Q_(:,1:q) = W
   deallocate(W, G, ev, work)
   nk = q
   allocate(T(nk,npert))
   call dgemm('T','N',nk,npert,nt,1.0d0,P_,nt,R,nt,0.0d0,T,nk)
   call dgemm('N','N',nt,npert,nk,1.0d0,P_,nt,T,nk,1.0d0,X,nt)
   call dgemm('N','N',nt,npert,nk,-1.0d0,Q_,nt,T,nk,1.0d0,R,nt)
   deallocate(T)
enddo
do p = 1,npert
   Ua(:,:,p) = reshape(X(1:na,p),(/nva,noa/)) ; Ub(:,:,p) = reshape(X(na+1:nt,p),(/nvb,nob/))
enddo
deallocate(dm, X, R, Z, P_, Q_, conv)

contains

subroutine applyA(n, Pin, Qo)
integer,intent(in) :: n
real(8),intent(in) :: Pin(nt,*)
real(8),intent(inout) :: Qo(nt,*)
real(8),allocatable :: La(:,:,:), Lb(:,:,:), VCa(:,:,:), VCb(:,:,:), Ta(:,:), Tb(:,:)
integer :: c
allocate(La(nConts,noa,n), Lb(nConts,nob,n), VCa(nConts,noa,n), VCb(nConts,nob,n))
!$omp parallel do private(c)
do c = 1,n
   call dgemm('N','N',nConts,noa,nva,1.0d0,Ca(1,noa+1),nConts,Pin(1,c),nva,0.0d0,La(1,1,c),nConts)
   call dgemm('N','N',nConts,nob,nvb,1.0d0,Cb(1,nob+1),nConts,Pin(na+1,c),nvb,0.0d0,Lb(1,1,c),nConts)
enddo
!$omp end parallel do
call cpscf_vind_spin_c(nConts, n, HF_exchange_frac, noa, nob, Ca(:,1:noa), Cb(:,1:nob), La, Lb, VCa, VCb)
!$omp parallel do private(c)
do c = 1,n
   call dgemm('T','N',nva,noa,nConts,1.0d0,Ca(1,noa+1),nConts,VCa(1,1,c),nConts,0.0d0,Qo(1,c),nva)
   call dgemm('T','N',nvb,nob,nConts,1.0d0,Cb(1,nob+1),nConts,VCb(1,1,c),nConts,0.0d0,Qo(na+1,c),nvb)
   Qo(1:nt,c) = Qo(1:nt,c) + dm*Pin(1:nt,c)
enddo
!$omp end parallel do
deallocate(La, Lb, VCa, VCb)
end subroutine applyA
end subroutine cpscf_solve_vo_bcg_spin

subroutine cpscf_solve_uhf(nConts, natoms, nocc_a, nocc_b, C_a, eLev_a, C_b, eLev_b, &
                           Pa, Pb, HF_exchange_frac, &
                           U_full_a, U_full_b, h1_MO_all_a, s1_MO_all_a, mo_e1_all_a, &
                           h1_MO_all_b, s1_MO_all_b, mo_e1_all_b)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force, integrals_compute_force_1e, integrals_compute_force_exchange_lr
use mod_integrals, only: direct_mode
use mod_profile, only: prof_start, prof_stop, prof_add, prof_add_cpu
use mod_xc_hess, only: xc_vxc_deriv_occ_u
use mod_exchange, only: cosx_enabled, RS_omega
implicit none
integer,intent(in) :: nConts, natoms, nocc_a, nocc_b
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: C_b(nConts,nConts), eLev_b(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: U_full_a(nConts,nocc_a,3,natoms), U_full_b(nConts,nocc_b,3,natoms)
real(8),intent(out) :: h1_MO_all_a(nConts,nocc_a,3,natoms), s1_MO_all_a(nConts,nocc_a,3,natoms)
real(8),intent(out) :: mo_e1_all_a(nocc_a,nocc_a,3,natoms)
real(8),intent(out) :: h1_MO_all_b(nConts,nocc_b,3,natoms), s1_MO_all_b(nConts,nocc_b,3,natoms)
real(8),intent(out) :: mo_e1_all_b(nocc_b,nocc_b,3,natoms)

integer :: nvir_a, nvir_b
integer :: npert, ipert
integer :: pchunk, pstart, pend, pn, pidx
integer(8), parameter :: pchunk_cache_bytes = 2_8*1024_8*1024_8
real(8),allocatable :: s1oo_all_a(:,:,:), s1oo_all_b(:,:,:)
real(8) :: tstage
real(8) :: omp_get_wtime
external :: omp_get_wtime
real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: RHSfix_a(:,:,:), RHSfix_b(:,:,:), Uvo_all_a(:,:,:), Uvo_all_b(:,:,:)
real(8),allocatable :: dF_a(:,:), dF_b(:,:), dS_full(:,:)
real(8),allocatable :: h1_MO_a(:,:), s1_MO_a(:,:), h1_MO_b(:,:), s1_MO_b(:,:)
real(8),allocatable :: s1oo_a(:,:), s1oo_b(:,:), Tmpa(:,:), Tmpb(:,:)
real(8),allocatable :: Uvo_a(:,:), Uvo_b(:,:), Pvo_a(:,:), Pvo_b(:,:), v1oo_a(:,:), v1oo_b(:,:)
real(8),allocatable :: Dtr(:,:,:,:), Dtot(:,:,:), Vtr(:,:,:,:), Kscr(:,:,:,:)
integer :: iat,ic,i,j,a,mu
logical :: fm
real(8),allocatable :: Lfa(:,:,:), Lfb(:,:,:), VCa(:,:,:), VCb(:,:,:)
real(8) :: conv_thr

nvir_a = nConts - nocc_a
nvir_b = nConts - nocc_b

cpscf_ecp_atom = -1
allocate(dS(nConts,nConts,3), dHcore(nConts,nConts,3))
allocate(dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3))
tstage = omp_get_wtime()
call cpu_time(tstage_c)
call prof_start("hess_force_int")
if (allocated(cpscf_dfj_jx) .and. (allocated(cpscf_cosx_kxa) .or. HF_exchange_frac == 0.0d0)) then
   call integrals_compute_force_1e(nConts, dS, dHcore)
   dJi = 0.0d0 ; dKa = 0.0d0 ; dKb = 0.0d0
else
   call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
endif
call prof_stop("hess_force_int")
if (RS_omega > 0.0d0 .and. .not. cosx_enabled) then
   if (allocated(cpscf_lr_dKa)) deallocate(cpscf_lr_dKa, cpscf_lr_dKb)
   allocate(cpscf_lr_dKa(nConts,nConts,3), cpscf_lr_dKb(nConts,nConts,3))
   call integrals_compute_force_exchange_lr(nConts, Pa, Pb, cpscf_lr_dKa, cpscf_lr_dKb)
endif
cphf_t_forceint = omp_get_wtime() - tstage

allocate(ao_atom(nConts))
mu = 0
do iat = 1,natoms
   do i = 1,atoms(iat)%nconts
      mu = mu+1
      ao_atom(mu) = iat
   enddo
enddo

cphf_t_amat = 0.0d0
allocate(dF_a(nConts,nConts), dF_b(nConts,nConts), dS_full(nConts,nConts))
allocate(h1_MO_a(nConts,nConts), s1_MO_a(nConts,nConts))
allocate(h1_MO_b(nConts,nConts), s1_MO_b(nConts,nConts))
allocate(s1oo_a(nocc_a,nocc_a), s1oo_b(nocc_b,nocc_b))
allocate(Tmpa(nConts,nocc_a), Tmpb(nConts,nocc_b))

U_full_a = 0.0d0
U_full_b = 0.0d0
cphf_t_fockderiv = 0.0d0 ; cphf_c_fockderiv = 0.0d0
cphf_t_iter = 0.0d0
if (HF_exchange_frac < 1.0d0) then
   cpscf_xc_uks = .true.
   if (allocated(cpscf_xc_Pa)) deallocate(cpscf_xc_Pa, cpscf_xc_Pb)
   allocate(cpscf_xc_Pa(nConts,nConts), cpscf_xc_Pb(nConts,nConts))
   cpscf_xc_Pa = Pa ; cpscf_xc_Pb = Pb
   if (allocated(cpscf_xc_vxca)) deallocate(cpscf_xc_vxca, cpscf_xc_vxcb)
   allocate(cpscf_xc_vxca(nConts,nocc_a,3,natoms), cpscf_xc_vxcb(nConts,nocc_b,3,natoms))
   tstage = omp_get_wtime()
   call cpu_time(tstage_c)
   call prof_start("fd_xc_vxc_deriv")
   call xc_vxc_deriv_occ_u(nConts, natoms, nocc_a, nocc_b, Pa, Pb, C_a(:,1:nocc_a), C_b(:,1:nocc_b), &
                           cpscf_xc_vxca, cpscf_xc_vxcb)
   call prof_stop("fd_xc_vxc_deriv")
   cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
   call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
endif
if (cpscf_xc_uks) then
   block
      use mod_xc_hess, only: xc_resp_grid_build
      use GRID_info, only: XCGRID_TIER_LO
      character(8) :: e
      call get_environment_variable("ENGINE_CPSCF_XC_GRID", e)
      if (trim(e) .ne. "full") call xc_resp_grid_build(XCGRID_TIER_LO)
   end block
endif
cphf_t_pass1 = 0.0d0
cphf_t_pass3 = 0.0d0
cphf_t_mo = 0.0d0
cphf_n_iter_total = 0

npert = 3*natoms
allocate(RHSfix_a(nvir_a,nocc_a,npert), RHSfix_b(nvir_b,nocc_b,npert))
allocate(Uvo_all_a(nvir_a,nocc_a,npert), Uvo_all_b(nvir_b,nocc_b,npert))
allocate(s1oo_all_a(nocc_a,nocc_a,npert), s1oo_all_b(nocc_b,nocc_b,npert))

if (direct_mode) then
   pchunk = npert
else
   pchunk = int(max(1_8, pchunk_cache_bytes / (7_8*int(nConts,8)*int(nConts,8)*8_8)))
   pchunk = max(1, min(npert, pchunk))
endif
block
   character(len=8) :: e
   call get_environment_variable("ENGINE_CPSCF_UHF_DENSE", e)
   fm = (trim(e) .ne. "1")
end block
if (fm) then
   allocate(Dtr(1,1,1,2), Dtot(1,1,1), Vtr(1,1,1,2), Kscr(1,1,1,2))
else
   allocate(Dtr(nConts,nConts,pchunk,2), Dtot(nConts,nConts,pchunk))
   allocate(Vtr(nConts,nConts,pchunk,2), Kscr(nConts,nConts,pchunk,2))
endif

ipert = 0
do iat = 1,natoms
   do ic = 1,3
      ipert = ipert + 1
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      call cphf_fock_deriv_atom_spin(nConts, ao_atom, iat, ic, dS, dHcore, dJi, dKa, dKb, &
                                     HF_exchange_frac, Pa, Pb, dF_a, dF_b, dS_full)
      cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
      call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      h1_MO_a = matmul(transpose(C_a), matmul(dF_a, C_a))
      s1_MO_a = matmul(transpose(C_a), matmul(dS_full, C_a))
      h1_MO_b = matmul(transpose(C_b), matmul(dF_b, C_b))
      if (cpscf_xc_uks) then
         h1_MO_a(:,1:nocc_a) = h1_MO_a(:,1:nocc_a) + matmul(transpose(C_a), cpscf_xc_vxca(:,:,ic,iat))
         h1_MO_b(:,1:nocc_b) = h1_MO_b(:,1:nocc_b) + matmul(transpose(C_b), cpscf_xc_vxcb(:,:,ic,iat))
      endif
      if (cpscf_cosx_kx_proj) then
         h1_MO_a(:,1:nocc_a) = h1_MO_a(:,1:nocc_a) &
            - HF_exchange_frac*matmul(transpose(C_a), cpscf_cosx_kxa(:,1:nocc_a,ic,iat))
         h1_MO_b(:,1:nocc_b) = h1_MO_b(:,1:nocc_b) &
            - HF_exchange_frac*matmul(transpose(C_b), cpscf_cosx_kxb(:,1:nocc_b,ic,iat))
      endif
      s1_MO_b = matmul(transpose(C_b), matmul(dS_full, C_b))
      h1_MO_all_a(:,:,ic,iat) = h1_MO_a(:,1:nocc_a)
      s1_MO_all_a(:,:,ic,iat) = s1_MO_a(:,1:nocc_a)
      h1_MO_all_b(:,:,ic,iat) = h1_MO_b(:,1:nocc_b)
      s1_MO_all_b(:,:,ic,iat) = s1_MO_b(:,1:nocc_b)
      cphf_t_mo = cphf_t_mo + (omp_get_wtime() - tstage)

      s1oo_all_a(:,:,ipert) = s1_MO_a(1:nocc_a,1:nocc_a)
      s1oo_all_b(:,:,ipert) = s1_MO_b(1:nocc_b,1:nocc_b)
   enddo
enddo

tstage = omp_get_wtime()

call cpu_time(tstage_c)
if (fm) then
   allocate(Lfa(nConts,nocc_a,npert), Lfb(nConts,nocc_b,npert), VCa(nConts,nocc_a,npert), VCb(nConts,nocc_b,npert))
   do ipert = 1,npert
      call dgemm('N','N',nConts,nocc_a,nocc_a,-0.5d0,C_a,nConts,s1oo_all_a(1,1,ipert),nocc_a,0.0d0,Lfa(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_b,nocc_b,-0.5d0,C_b,nConts,s1oo_all_b(1,1,ipert),nocc_b,0.0d0,Lfb(1,1,ipert),nConts)
   enddo
   call cpscf_vind_spin_c(nConts, npert, HF_exchange_frac, nocc_a, nocc_b, C_a(:,1:nocc_a), C_b(:,1:nocc_b), &
                          Lfa, Lfb, VCa, VCb)
   do ipert = 1,npert
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      call dgemm('T','N',nvir_a,nocc_a,nConts,1.0d0,C_a(1,nocc_a+1),nConts,VCa(1,1,ipert),nConts,0.0d0,RHSfix_a(1,1,ipert),nvir_a)
      call dgemm('T','N',nvir_b,nocc_b,nConts,1.0d0,C_b(1,nocc_b+1),nConts,VCb(1,1,ipert),nConts,0.0d0,RHSfix_b(1,1,ipert),nvir_b)
      do i = 1,nocc_a
         do a = 1,nvir_a
            RHSfix_a(a,i,ipert) = RHSfix_a(a,i,ipert) + h1_MO_all_a(nocc_a+a,i,ic,iat) &
                                  - s1_MO_all_a(nocc_a+a,i,ic,iat)*eLev_a(i)
         enddo
      enddo
      do i = 1,nocc_b
         do a = 1,nvir_b
            RHSfix_b(a,i,ipert) = RHSfix_b(a,i,ipert) + h1_MO_all_b(nocc_b+a,i,ic,iat) &
                                  - s1_MO_all_b(nocc_b+a,i,ic,iat)*eLev_b(i)
         enddo
      enddo
   enddo
else
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      Dtr(:,:,pidx,1) = -matmul(C_a(:,1:nocc_a), &
                        matmul(s1oo_all_a(:,:,ipert), transpose(C_a(:,1:nocc_a))))
      Dtr(:,:,pidx,2) = -matmul(C_b(:,1:nocc_b), &
                        matmul(s1oo_all_b(:,:,ipert), transpose(C_b(:,1:nocc_b))))
      Dtot(:,:,pidx) = Dtr(:,:,pidx,1) + Dtr(:,:,pidx,2)
   enddo
   call cpscf_vind_spin(nConts, pn, 2, HF_exchange_frac, Dtr(:,:,1:pn,:), &
                        Dtot(:,:,1:pn), Vtr(:,:,1:pn,:), Kscr(:,:,1:pn,:))
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO_a(:,1:nocc_a) = h1_MO_all_a(:,:,ic,iat)
      s1_MO_a(:,1:nocc_a) = s1_MO_all_a(:,:,ic,iat)
      h1_MO_b(:,1:nocc_b) = h1_MO_all_b(:,:,ic,iat)
      s1_MO_b(:,1:nocc_b) = s1_MO_all_b(:,:,ic,iat)
      Tmpa = matmul(Vtr(:,:,pidx,1), C_a(:,1:nocc_a))
      Tmpb = matmul(Vtr(:,:,pidx,2), C_b(:,1:nocc_b))
      do i = 1,nocc_a
         do a = 1,nvir_a
            RHSfix_a(a,i,ipert) = h1_MO_a(nocc_a+a,i) - s1_MO_a(nocc_a+a,i)*eLev_a(i) &
                                  + dot_product(C_a(:,nocc_a+a), Tmpa(:,i))
         enddo
      enddo
      do i = 1,nocc_b
         do a = 1,nvir_b
            RHSfix_b(a,i,ipert) = h1_MO_b(nocc_b+a,i) - s1_MO_b(nocc_b+a,i)*eLev_b(i) &
                                  + dot_product(C_b(:,nocc_b+a), Tmpb(:,i))
         enddo
      enddo
   enddo
   pstart = pend + 1
enddo
endif
cphf_t_pass1 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass1 = ctmp_c - tstage_c
call prof_add("hess_cphf_rhs", cphf_t_pass1)

call prof_add_cpu("hess_cphf_rhs", cphf_c_pass1)
tstage = omp_get_wtime()

call cpu_time(tstage_c)
if (fm) then
   conv_thr = 5.0d-5
   block
      character(32) :: envval
      integer :: envstat
      call get_environment_variable("ENGINE_CPHF_CONV", envval, status=envstat)
      if (envstat .eq. 0) read(envval,*) conv_thr
   end block
   block
      character(len=16) :: solv
      call get_environment_variable("ENGINE_CPHF_SOLVER", solv)
      if (trim(solv) .eq. "pcg") then
         call cpscf_solve_vo_pcg_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, C_a, eLev_a, C_b, eLev_b, &
                                      HF_exchange_frac, RHSfix_a, RHSfix_b, Uvo_all_a, Uvo_all_b, &
                                      cphf_n_iter_total, conv_thr)
      else
         call cpscf_solve_vo_bcg_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, C_a, eLev_a, C_b, eLev_b, &
                                      HF_exchange_frac, RHSfix_a, RHSfix_b, Uvo_all_a, Uvo_all_b, &
                                      cphf_n_iter_total, conv_thr)
      endif
   end block
else
   call cpscf_solve_vo_direct_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, &
                                   C_a, eLev_a, C_b, eLev_b, HF_exchange_frac, &
                                   RHSfix_a, RHSfix_b, Uvo_all_a, Uvo_all_b, cphf_n_iter_total)
endif
cphf_t_iter = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_iter = ctmp_c - tstage_c
call prof_add("hess_cphf_solve", cphf_t_iter)

call prof_add_cpu("hess_cphf_solve", cphf_c_iter)
allocate(Uvo_a(nvir_a,nocc_a), Uvo_b(nvir_b,nocc_b))
allocate(Pvo_a(nConts,nConts), Pvo_b(nConts,nConts))
allocate(v1oo_a(nocc_a,nocc_a), v1oo_b(nocc_b,nocc_b))

tstage = omp_get_wtime()

call cpu_time(tstage_c)
if (fm) then
   do ipert = 1,npert
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      call dgemm('N','N',nConts,nocc_a,nocc_a,-0.5d0,C_a,nConts,s1oo_all_a(1,1,ipert),nocc_a,0.0d0,Lfa(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_a,nvir_a,1.0d0,C_a(1,nocc_a+1),nConts,Uvo_all_a(1,1,ipert),nvir_a,1.0d0,Lfa(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_b,nocc_b,-0.5d0,C_b,nConts,s1oo_all_b(1,1,ipert),nocc_b,0.0d0,Lfb(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_b,nvir_b,1.0d0,C_b(1,nocc_b+1),nConts,Uvo_all_b(1,1,ipert),nvir_b,1.0d0,Lfb(1,1,ipert),nConts)
      do i = 1,nocc_a
         do j = 1,nocc_a
            U_full_a(j,i,ic,iat) = -0.5d0*s1oo_all_a(j,i,ipert)
         enddo
         do a = 1,nvir_a
            U_full_a(nocc_a+a,i,ic,iat) = Uvo_all_a(a,i,ipert)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            U_full_b(j,i,ic,iat) = -0.5d0*s1oo_all_b(j,i,ipert)
         enddo
         do a = 1,nvir_b
            U_full_b(nocc_b+a,i,ic,iat) = Uvo_all_b(a,i,ipert)
         enddo
      enddo
   enddo
   call cpscf_vind_spin_c(nConts, npert, HF_exchange_frac, nocc_a, nocc_b, C_a(:,1:nocc_a), C_b(:,1:nocc_b), &
                          Lfa, Lfb, VCa, VCb)
   do ipert = 1,npert
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      call dgemm('T','N',nocc_a,nocc_a,nConts,1.0d0,C_a,nConts,VCa(1,1,ipert),nConts,0.0d0,v1oo_a,nocc_a)
      call dgemm('T','N',nocc_b,nocc_b,nConts,1.0d0,C_b,nConts,VCb(1,1,ipert),nConts,0.0d0,v1oo_b,nocc_b)
      do i = 1,nocc_a
         do j = 1,nocc_a
            mo_e1_all_a(i,j,ic,iat) = h1_MO_all_a(i,j,ic,iat) - 0.5d0*(eLev_a(i)+eLev_a(j))*s1_MO_all_a(i,j,ic,iat) &
                                      + v1oo_a(i,j)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            mo_e1_all_b(i,j,ic,iat) = h1_MO_all_b(i,j,ic,iat) - 0.5d0*(eLev_b(i)+eLev_b(j))*s1_MO_all_b(i,j,ic,iat) &
                                      + v1oo_b(i,j)
         enddo
      enddo
   enddo
   deallocate(Lfa, Lfb, VCa, VCb)
else
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      Uvo_a = Uvo_all_a(:,:,ipert)
      Uvo_b = Uvo_all_b(:,:,ipert)
      Pvo_a = matmul(C_a(:,nocc_a+1:nConts), matmul(Uvo_a, transpose(C_a(:,1:nocc_a)))) &
            + matmul(C_a(:,1:nocc_a), matmul(transpose(Uvo_a), transpose(C_a(:,nocc_a+1:nConts))))
      Pvo_b = matmul(C_b(:,nocc_b+1:nConts), matmul(Uvo_b, transpose(C_b(:,1:nocc_b)))) &
            + matmul(C_b(:,1:nocc_b), matmul(transpose(Uvo_b), transpose(C_b(:,nocc_b+1:nConts))))
      Dtr(:,:,pidx,1) = -matmul(C_a(:,1:nocc_a), &
                        matmul(s1oo_all_a(:,:,ipert), transpose(C_a(:,1:nocc_a)))) + Pvo_a
      Dtr(:,:,pidx,2) = -matmul(C_b(:,1:nocc_b), &
                        matmul(s1oo_all_b(:,:,ipert), transpose(C_b(:,1:nocc_b)))) + Pvo_b
      Dtot(:,:,pidx) = Dtr(:,:,pidx,1) + Dtr(:,:,pidx,2)

      do i = 1,nocc_a
         do j = 1,nocc_a
            U_full_a(j,i,ic,iat) = -0.5d0*s1oo_all_a(j,i,ipert)
         enddo
         do a = 1,nvir_a
            U_full_a(nocc_a+a,i,ic,iat) = Uvo_a(a,i)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            U_full_b(j,i,ic,iat) = -0.5d0*s1oo_all_b(j,i,ipert)
         enddo
         do a = 1,nvir_b
            U_full_b(nocc_b+a,i,ic,iat) = Uvo_b(a,i)
         enddo
      enddo
   enddo
   call cpscf_vind_spin(nConts, pn, 2, HF_exchange_frac, Dtr(:,:,1:pn,:), &
                        Dtot(:,:,1:pn), Vtr(:,:,1:pn,:), Kscr(:,:,1:pn,:))
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO_a(:,1:nocc_a) = h1_MO_all_a(:,:,ic,iat)
      s1_MO_a(:,1:nocc_a) = s1_MO_all_a(:,:,ic,iat)
      h1_MO_b(:,1:nocc_b) = h1_MO_all_b(:,:,ic,iat)
      s1_MO_b(:,1:nocc_b) = s1_MO_all_b(:,:,ic,iat)
      v1oo_a = matmul(transpose(C_a(:,1:nocc_a)), matmul(Vtr(:,:,pidx,1), C_a(:,1:nocc_a)))
      v1oo_b = matmul(transpose(C_b(:,1:nocc_b)), matmul(Vtr(:,:,pidx,2), C_b(:,1:nocc_b)))
      do i = 1,nocc_a
         do j = 1,nocc_a
            mo_e1_all_a(i,j,ic,iat) = h1_MO_a(i,j) - 0.5d0*(eLev_a(i)+eLev_a(j))*s1_MO_a(i,j) + v1oo_a(i,j)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            mo_e1_all_b(i,j,ic,iat) = h1_MO_b(i,j) - 0.5d0*(eLev_b(i)+eLev_b(j))*s1_MO_b(i,j) + v1oo_b(i,j)
         enddo
      enddo
   enddo
   pstart = pend + 1
enddo
endif
cphf_t_pass3 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass3 = ctmp_c - tstage_c
call prof_add("hess_cphf_obs", cphf_t_pass3)

call prof_add_cpu("hess_cphf_obs", cphf_c_pass3)
call prof_add("hess_fock_deriv", cphf_t_fockderiv)
call prof_add_cpu("hess_fock_deriv", cphf_c_fockderiv)
deallocate(dS,dHcore,dJi,dKa,dKb,ao_atom)
if (allocated(cpscf_lr_dKa)) deallocate(cpscf_lr_dKa, cpscf_lr_dKb)
deallocate(RHSfix_a,RHSfix_b,Uvo_all_a,Uvo_all_b,s1oo_all_a,s1oo_all_b)
deallocate(dF_a,dF_b,dS_full,h1_MO_a,s1_MO_a,h1_MO_b,s1_MO_b,s1oo_a,s1oo_b,Tmpa,Tmpb)
deallocate(Uvo_a,Uvo_b,Pvo_a,Pvo_b,v1oo_a,v1oo_b)
deallocate(Dtr,Dtot,Vtr,Kscr)
if (allocated(cpscf_xc_vxca)) deallocate(cpscf_xc_vxca, cpscf_xc_vxcb)
if (allocated(cpscf_xc_Pa)) deallocate(cpscf_xc_Pa, cpscf_xc_Pb)
if (cpscf_xc_uks) then
   block
      use mod_xc_hess, only: xc_resp_grid_clear
      call xc_resp_grid_clear()
   end block
endif
cpscf_xc_uks = .false.
end subroutine cpscf_solve_uhf

subroutine cphf_hess_response(nConts, natoms, nocc, U_full, h1_MO_all, s1_MO_all, mo_e1_all, eLev_a, d2E_response)
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: U_full(nConts,nocc,3,natoms)
real(8),intent(in) :: h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms)
real(8),intent(in) :: mo_e1_all(nocc,nocc,3,natoms)
real(8),intent(in) :: eLev_a(nConts)
real(8),intent(out) :: d2E_response(3,natoms,3,natoms)

integer :: ia,ja,x,y,i
real(8),allocatable :: hs(:,:)
real(8) :: t1,t3

allocate(hs(nConts,nocc))
do ia = 1,natoms
   do x = 1,3
      do i = 1,nocc
         hs(:,i) = h1_MO_all(:,i,x,ia) - eLev_a(i)*s1_MO_all(:,i,x,ia)
      enddo
      do ja = 1,natoms
         do y = 1,3
            t1 = sum(hs*U_full(:,:,y,ja))
            t3 = sum(s1_MO_all(1:nocc,:,x,ia)*mo_e1_all(:,:,y,ja))
            d2E_response(x,ia,y,ja) = 4.0d0*t1 - 2.0d0*t3
         enddo
      enddo
   enddo
enddo
deallocate(hs)
end subroutine cphf_hess_response

subroutine cphf_hess_response_spin(nConts, natoms, nocc_a, nocc_b, &
                                    U_full_a, U_full_b, &
                                    h1_MO_all_a, s1_MO_all_a, mo_e1_all_a, eLev_a, &
                                    h1_MO_all_b, s1_MO_all_b, mo_e1_all_b, eLev_b, &
                                    d2E_response)
implicit none
integer,intent(in) :: nConts, natoms, nocc_a, nocc_b
real(8),intent(in) :: U_full_a(nConts,nocc_a,3,natoms), U_full_b(nConts,nocc_b,3,natoms)
real(8),intent(in) :: h1_MO_all_a(nConts,nocc_a,3,natoms), s1_MO_all_a(nConts,nocc_a,3,natoms)
real(8),intent(in) :: mo_e1_all_a(nocc_a,nocc_a,3,natoms)
real(8),intent(in) :: h1_MO_all_b(nConts,nocc_b,3,natoms), s1_MO_all_b(nConts,nocc_b,3,natoms)
real(8),intent(in) :: mo_e1_all_b(nocc_b,nocc_b,3,natoms)
real(8),intent(in) :: eLev_a(nConts), eLev_b(nConts)
real(8),intent(out) :: d2E_response(3,natoms,3,natoms)

integer :: ia,ja,x,y,i
real(8),allocatable :: hs_a(:,:), hs_b(:,:)
real(8) :: t1a,t1b,t3a,t3b

allocate(hs_a(nConts,nocc_a), hs_b(nConts,nocc_b))
do ia = 1,natoms
   do x = 1,3
      do i = 1,nocc_a
         hs_a(:,i) = h1_MO_all_a(:,i,x,ia) - eLev_a(i)*s1_MO_all_a(:,i,x,ia)
      enddo
      do i = 1,nocc_b
         hs_b(:,i) = h1_MO_all_b(:,i,x,ia) - eLev_b(i)*s1_MO_all_b(:,i,x,ia)
      enddo
      do ja = 1,natoms
         do y = 1,3
            t1a = sum(hs_a*U_full_a(:,:,y,ja))
            t1b = sum(hs_b*U_full_b(:,:,y,ja))
            t3a = sum(s1_MO_all_a(1:nocc_a,:,x,ia)*mo_e1_all_a(:,:,y,ja))
            t3b = sum(s1_MO_all_b(1:nocc_b,:,x,ia)*mo_e1_all_b(:,:,y,ja))
            d2E_response(x,ia,y,ja) = 2.0d0*t1a + 2.0d0*t1b - t3a - t3b
         enddo
      enddo
   enddo
enddo
deallocate(hs_a, hs_b)
end subroutine cphf_hess_response_spin

subroutine cphf_hess_nuc(natoms, znum, ecpce, coor_bohr, d2Enn)
implicit none
integer,intent(in) :: natoms, znum(natoms), ecpce(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(out) :: d2Enn(3,natoms,3,natoms)

integer :: i,k,x,y
real(8) :: r(3), r1, qi, qk, blk(3,3)

d2Enn = 0.0d0
do i = 1,natoms
   qi = real(znum(i)-ecpce(i),8)
   do k = 1,natoms
      if (k == i) cycle
      qk = real(znum(k)-ecpce(k),8)
      r = coor_bohr(:,i)-coor_bohr(:,k)
      r1 = sqrt(sum(r**2))
      do y = 1,3
         do x = 1,3
            blk(x,y) = qi*qk*3.0d0*r(x)*r(y)/r1**5
         enddo
         blk(y,y) = blk(y,y) - qi*qk/r1**3
      enddo
      d2Enn(:,i,:,i) = d2Enn(:,i,:,i) + blk
      d2Enn(:,i,:,k) = d2Enn(:,i,:,k) - blk
   enddo
enddo
end subroutine cphf_hess_nuc

subroutine cphf_solve_efield(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, alpha, niter)
use mod_integrals, only: integrals_dipole_ao, engine_use_df_j
use mod_exchange, only: RS_omega, cosx_enabled
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: alpha(3,3)
integer,intent(out) :: niter
integer :: nvir, x, y
real(8),allocatable :: R(:,:,:), T(:,:), RHS(:,:,:), U(:,:,:)

if (RS_omega > 0.0d0 .and. .not. cosx_enabled) &
   error stop 'polarizability: range-separated hybrids with exact exchange are not implemented'
nvir = nConts - nocc
cpscf_xc_active = (HF_exchange_frac < 1.0d0)
if (cpscf_xc_active) then
   if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
   allocate(cpscf_xc_Ptot(nConts,nConts))
   cpscf_xc_Ptot = Pa + Pb
endif
cpscf_force_dfj = engine_use_df_j
cphf_conv_override = 1.0d-8
cphf_voo_track = .false.
cphf_t_vj = 0.0d0 ; cphf_t_vk = 0.0d0 ; cphf_t_vxc = 0.0d0
cphf_n_vk = 0 ; cphf_nd_vk = 0

allocate(R(nConts,nConts,3), T(nConts,nocc), RHS(nvir,nocc,3), U(nvir,nocc,3))
call integrals_dipole_ao(nConts, R)
do x = 1,3
   call dgemm('N','N',nConts,nocc,nConts,1.0d0,R(1,1,x),nConts,C_a,nConts,0.0d0,T,nConts)
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,T,nConts,0.0d0,RHS(1,1,x),nvir)
enddo
call cphf_solve_vo_direct(nConts, nocc, nvir, 3, C_a, eLev_a, HF_exchange_frac, RHS, U, niter)
do x = 1,3
   do y = 1,3
      alpha(x,y) = -4.0d0*sum(RHS(:,:,x)*U(:,:,y))
   enddo
enddo
alpha = 0.5d0*(alpha + transpose(alpha))
cpscf_force_dfj = .false.
cphf_conv_override = -1.0d0
cpscf_xc_active = .false.
deallocate(R, T, RHS, U)
end subroutine cphf_solve_efield

end module mod_cpscf
