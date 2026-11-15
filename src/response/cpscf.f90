! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: the coupled-perturbed self-consistent-field response

module mod_cpscf
use MOL_info, only: engine_verbose
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
public :: cpscf_use_vind_c, cpscf_vind_c, cpscf_dfj_response, cpscf_cosx_k
public :: cpscf_ecp_centre, cpscf_lr_density_role, cphf_solve_vo_direct
public :: cpscf_apb_rhf, cpscf_apb_amb_rhf, cpscf_apb_amb1_rhf, cpscf_response_begin, cpscf_response_end
public :: cpscf_apb_amb_uhf, cpscf_response_begin_uks
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
logical,save,public :: cpscf_triplet = .false.
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
real(8),save :: cphf_omega2 = 0.0d0

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

module subroutine cphf_debug_hessian_row(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic &
                                         )
implicit none
integer,intent(in) :: nConts, natoms, nocc, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
end subroutine cphf_debug_hessian_row

module subroutine cphf_debug_hessian_row_uhf(nConts, natoms, nocc_a, nocc_b, C_a, eLev_a, C_b, eLev_b, &
                                       Pa, Pb, HF_exchange_frac, iatom, ic &
                                       )
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

module subroutine cpscf_vind_spin(nConts, nbatch, nspin, HF_exchange_frac, Dtrial, Dtot, Vout, Kscr)
implicit none
integer,intent(in) :: nConts, nbatch, nspin
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Dtrial(nConts,nConts,nbatch,nspin)
real(8),intent(in) :: Dtot(nConts,nConts,nbatch)
real(8),intent(inout) :: Vout(nConts,nConts,nbatch,nspin)
real(8),intent(inout) :: Kscr(nConts,nConts,nbatch,nspin)
end subroutine cpscf_vind_spin

module subroutine cpscf_solve_uhf(nConts, natoms, nocc_a, nocc_b, C_a, eLev_a, C_b, eLev_b, &
                                  Pa, Pb, HF_exchange_frac, &
                                  U_full_a, U_full_b, h1_MO_all_a, s1_MO_all_a, mo_e1_all_a, &
                                  h1_MO_all_b, s1_MO_all_b, mo_e1_all_b)
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
end subroutine cpscf_solve_uhf

module subroutine cphf_solve_efield(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, alpha, niter, omega)
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: alpha(3,3)
integer,intent(out) :: niter
real(8),intent(in),optional :: omega
end subroutine cphf_solve_efield
module subroutine cphf_solve_vo_direct(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                                       RHSfix, Uvo_all, niter_total)
implicit none
integer,intent(in) :: nConts, nocc, nvir, npert
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: RHSfix(nvir,nocc,npert)
real(8),intent(out) :: Uvo_all(nvir,nocc,npert)
integer,intent(out) :: niter_total
end subroutine cphf_solve_vo_direct

module subroutine cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, &
                                 U_full, h1_MO_all, s1_MO_all, mo_e1_all)
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: U_full(nConts,nocc,3,natoms)
real(8),intent(out) :: h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms)
real(8),intent(out) :: mo_e1_all(nocc,nocc,3,natoms)
end subroutine cphf_solve_rhf
module subroutine cpscf_response_begin(nConts, Pa, Pb, HF_exchange_frac, scf_grid)
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
logical,intent(in),optional :: scf_grid
end subroutine cpscf_response_begin

module subroutine cpscf_response_end()
implicit none
end subroutine cpscf_response_end

module subroutine cpscf_apb_rhf(nConts, nocc, nvir, nvec, C_a, eLev_a, HF_exchange_frac, X, Q)
implicit none
integer,intent(in) :: nConts, nocc, nvir, nvec
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: X(nvir,nocc,nvec)
real(8),intent(out) :: Q(nvir,nocc,nvec)
end subroutine cpscf_apb_rhf

module subroutine cpscf_apb_amb_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, np, Xp, Qp, nm, Xm, Qm, Qoo)
implicit none
integer,intent(in) :: nConts, nocc, nvir, np, nm
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Xp(nvir,nocc,np), Xm(nvir,nocc,nm)
real(8),intent(out) :: Qp(nvir,nocc,np), Qm(nvir,nocc,nm)
real(8),intent(out),optional :: Qoo(nocc,nocc,np)
end subroutine cpscf_apb_amb_rhf

module subroutine cpscf_apb_amb_uhf(nConts, noa, nob, Ca, ea, Cb, eb, HF_exchange_frac, np, Xp, Qp, nm, Xm, Qm)
implicit none
integer,intent(in) :: nConts, noa, nob, np, nm
real(8),intent(in) :: Ca(nConts,nConts), ea(nConts), Cb(nConts,nConts), eb(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Xp((nConts-noa)*noa+(nConts-nob)*nob,np), Xm((nConts-noa)*noa+(nConts-nob)*nob,nm)
real(8),intent(out) :: Qp((nConts-noa)*noa+(nConts-nob)*nob,np), Qm((nConts-noa)*noa+(nConts-nob)*nob,nm)
end subroutine cpscf_apb_amb_uhf

module subroutine cpscf_response_begin_uks(nConts, Pa, Pb, HF_exchange_frac, scf_grid)
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
logical,intent(in),optional :: scf_grid
end subroutine cpscf_response_begin_uks

module subroutine cpscf_vind_spin_c(nConts, nb, HF_exchange_frac, noa, nob, Ca, Cb, La, Lb, VCa, VCb)
implicit none
integer,intent(in) :: nConts, nb, noa, nob
real(8),intent(in) :: HF_exchange_frac, Ca(nConts,noa), Cb(nConts,nob)
real(8),intent(in) :: La(nConts,noa,nb), Lb(nConts,nob,nb)
real(8),intent(out) :: VCa(nConts,noa,nb), VCb(nConts,nob,nb)
end subroutine cpscf_vind_spin_c

module subroutine cpscf_apb_amb1_rhf(nConts, nocc, nvir, C_a, eLev_a, HF_exchange_frac, nv, X, Qp, Qm)
implicit none
integer,intent(in) :: nConts, nocc, nvir, nv
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: X(nvir,nocc,nv)
real(8),intent(out) :: Qp(nvir,nocc,nv), Qm(nvir,nocc,nv)
end subroutine cpscf_apb_amb1_rhf
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
if (engine_verbose .ge. 2) then
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
if (engine_verbose .ge. 2) print '(A,I0,A,F9.3,A)', " COSX BATCH_C nd=", nbatch, "  ", omp_get_wtime()-tv0, " s"
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
use MOL_info, only: nPointCharges, pointcharge_q, pointcharge_coor
implicit none
INCLUDE 'parameter.h'
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
   do k = 1,nPointCharges
      qk = pointcharge_q(k)
      r = coor_bohr(:,i) - pointcharge_coor(:,k)*ans2bohr
      r1 = sqrt(sum(r**2))
      do y = 1,3
         do x = 1,3
            blk(x,y) = qi*qk*3.0d0*r(x)*r(y)/r1**5
         enddo
         blk(y,y) = blk(y,y) - qi*qk/r1**3
      enddo
      d2Enn(:,i,:,i) = d2Enn(:,i,:,i) + blk
   enddo
enddo
end subroutine cphf_hess_nuc

end module mod_cpscf
