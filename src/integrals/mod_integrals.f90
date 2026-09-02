! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_integrals: the sole boundary between the rest of Engine and the

module mod_integrals
implicit none
private
public :: integrals_init
public :: eri_get
public :: integrals_build_coulomb
public :: integrals_build_jk_fused
public :: jk_fusion_active
public :: integrals_compute_force
public :: integrals_compute_force_exchange_lr
public :: integrals_compute_force_1e
public :: integrals_hessian_1e_shell
public :: integrals_hessian_nuc_atom
public :: integrals_hessian_nuc_atom_raw
public :: integrals_cphf_2e_density_role
public :: integrals_hessian_2e_direct
public :: integrals_hessian_e2e_fixed_density
public :: integrals_hessian_debug_check2e_direct
public :: integrals_hessian_e1e_fixed_density
public :: integrals_hessian_debug_check1e
public :: integrals_nuc_attraction_deriv, integrals_nuc_deriv_allatoms
public :: integrals_ecp_force
public :: integrals_finalize
public :: integrals_set_accuracy
public :: integrals_build_exchange
public :: integrals_build_exchange_cosx
public :: integrals_build_exchange_cosx_lr
public :: integrals_build_exchange_cosx_sr
public :: integrals_force_exchange_cosx
public :: integrals_build_exchange_lr
public :: integrals_build_coulomb_df
public :: df_build_coulomb_core
public :: integrals_build_exchange_df
public :: integrals_build_exchange_df_lr
public :: integrals_force_df
public :: integrals_force_df_lr
public :: decide_df_force_mode
public :: df_direct_mode
public :: direct_mode
public :: df_force_direct_mode
public :: df_energy_need_gb
public :: df_energy_avail_gb
public :: df_force_need_gb
public :: df_force_avail_gb
public :: engine_use_df, engine_use_df_j, engine_use_df_k, engine_do_force
public :: cosx_promote_md
public :: cosx_force_hi_grid
public :: near_singular_overlap
public :: engine_puream
public :: norvec_is_identity
public :: engine_df_aux_basis
public :: engine_estimate_only
public :: integrals_estimate_aux_size
public :: nRec
public :: nContsAux
public :: compute_atomic_density
public :: calc_properties
public :: resp_charges
public :: engine_harris_guess
public :: cosmo_build_one_tess_matrix
public :: cosmo_build_one_tess_matrix_shellderiv
public :: esp_at_grid_batch
public :: esp_at_grid_batch_df

logical :: engine_use_df = .false.
logical :: engine_use_df_j = .false.
logical :: engine_use_df_k = .false.
logical :: engine_do_force = .false.
logical :: cosx_promote_md = .false.
logical :: cosx_force_hi_grid = .false.
logical :: engine_puream = .false.
logical :: norvec_is_identity = .false.
logical :: engine_harris_guess = .true.
character(len=64) :: engine_df_aux_basis = ""
logical :: engine_estimate_only = .false.

integer,allocatable :: atm(:,:)
integer,allocatable :: bas(:,:)
real(8),allocatable :: env(:)
integer             :: nBases = 0
integer,allocatable :: atm_nuc(:,:)
real(8),allocatable :: env_nuc(:)
integer             :: nAtoms_nuc = 0

integer,allocatable :: basECP(:,:)
real(8),allocatable :: envECP(:)
integer             :: necpbas = 0
real(8),allocatable :: NorVEC(:)
real(8),allocatable :: NorVECsph(:)
real(8),allocatable :: TWOEI(:)
integer(8)          :: nRec = 0

logical              :: direct_mode = .false.

real(8),allocatable :: schwarz_bound(:,:)
real(8),parameter   :: SCHWARZ_CUTOFF = 1.0d-8
logical :: near_singular_overlap = .false.

real(8) :: density_screen_cutoff = SCHWARZ_CUTOFF
real(8),parameter :: DENSITY_SCREEN_SCALE = 1.0d-5
real(8),parameter :: DENSITY_SCREEN_MAX = 1.0d-7

real(8) :: density_screen_cutoff_incr = SCHWARZ_CUTOFF
real(8),parameter :: DENSITY_SCREEN_INCR_FLOOR = 1.0d-16
real(8),parameter :: DENSITY_SCREEN_INCR_SCALE = 1.0d-6

real(8),parameter :: DENSITY_SCREEN_INCR_SCALE_J = 1.0d-7
real(8) :: density_screen_cutoff_incr_j = SCHWARZ_CUTOFF

integer(8) :: int2e_opt = 0_8

integer(8) :: int2e_opt_lr = 0_8

integer,allocatable :: ao_shell(:)

logical :: coulomb_list_built = .false.
integer,allocatable :: cPK_pq(:),cPK_rs(:)
real(8),allocatable :: cPK_val(:)
integer(8) :: cPK_n = 0
integer :: coulomb_npairs = 0

logical :: sig_pairs_built = .false.
integer,allocatable :: sig_pair_ishl(:),sig_pair_jshl(:)
integer :: n_sig_pairs = 0

logical :: exchange_list_built = .false.
integer,allocatable :: cK_i(:),cK_j(:),cK_k(:),cK_l(:)
real(8),allocatable :: cK_val(:)
integer(8) :: cK_n = 0

logical :: incr_fock_enabled = .true.
logical :: incr_fock_env_checked = .false.
integer :: incr_full_period = 20
integer :: incr_since_full = 0
logical :: incr_has_prev = .false.
real(8),allocatable :: incr_Ptot_prev(:,:), incr_J_prev(:,:)
real(8),allocatable :: incrf_Da_prev(:,:), incrf_Db_prev(:,:)
real(8),allocatable :: incrf_J_prev(:,:), incrf_Ka_prev(:,:), incrf_Kb_prev(:,:)
logical :: incrf_has_prev = .false.
integer :: incrf_since_full = 0

integer,parameter :: DF_NAUX_PER_L = 10
real(8),parameter :: DF_EXP_LO = 0.5d0, DF_EXP_HI = 2.5d0
logical :: df_built = .false.
integer :: nBasesAux = 0, nContsAux = 0
integer,allocatable :: basDF(:,:)
real(8),allocatable :: envDF(:)
real(8),allocatable :: NorVECAux(:)
real(8),allocatable :: dfB_compact(:,:)
real(8),allocatable :: dfB_compact_LR(:,:)
real(8),allocatable :: df_evec(:,:)
real(8),allocatable :: df_evalinv(:)
real(8),parameter :: DF_EVAL_FLOOR = 1.0d-8
integer(8) :: int3c2e_opt = 0_8, int2c2e_opt = 0_8
integer(8) :: int3c2e_opt_lr = 0_8
real(8),allocatable :: aux_shell_bound(:)
logical :: df_lr_optimizer_built = .false.
logical :: df_direct_mode = .false.
logical :: df_force_direct_mode = .false.
real(8) :: df_energy_need_gb = 0.0d0, df_energy_avail_gb = 0.0d0
real(8) :: df_force_need_gb = 0.0d0, df_force_avail_gb = 0.0d0
real(8),allocatable :: df_direct_Jlocal_pool(:,:,:)
real(8),allocatable :: df_direct_bveclocal_pool(:,:)
logical :: df_direct_pool_built = .false.
real(8),allocatable :: df_direct_cvec_cache(:)
real(8),allocatable :: df_direct_cvec_cache_Ptot(:,:)
logical :: df_direct_cvec_cache_valid = .false.
integer :: df_n_triples = 0
integer(2),allocatable :: df_triple_ish(:), df_triple_jsh(:), df_triple_ksh(:)
integer,allocatable :: df_shell_ncgto(:), df_shell_offset(:)

logical :: dfK_built = .false.
logical :: df_full_lr_metric_built = .false.
real(8),allocatable :: dfK(:,:,:)
real(8),allocatable :: dfK_compact(:,:)
real(8),allocatable :: dfC_full_compact(:,:)
integer,allocatable :: df_pair_row(:), df_pair_col(:)
integer,allocatable :: df_row_start(:), df_row_count(:)
integer :: df_total_pairs = 0, df_max_row_count = 0
integer,allocatable :: df_pair_mirror(:)
logical :: df_Minvhalf_built = .false.
real(8),allocatable :: df_Minvhalf(:,:)
logical :: df_Minv_full_built = .false.
real(8),allocatable :: df_Minv_full(:,:)

interface
   subroutine openblas_set_num_threads(num_threads) bind(c, name="openblas_set_num_threads")
   use iso_c_binding, only: c_int
   integer(c_int),value :: num_threads
   end subroutine openblas_set_num_threads
end interface

interface
module subroutine integrals_init(info)
use MOL_info
use mod_meminfo, only: get_memory_budget_bytes, get_store_threshold_bytes
    implicit none
    integer :: info
end subroutine integrals_init

module subroutine ensure_int2e_optimizer()
    implicit none
end subroutine ensure_int2e_optimizer

module subroutine integrals_set_accuracy(prms)
    implicit none
    real(8),intent(in) :: prms
end subroutine integrals_set_accuracy

module subroutine compute_schwarz_bounds()
    implicit none
end subroutine compute_schwarz_bounds

module subroutine build_significant_pairs()
    implicit none
end subroutine build_significant_pairs

module real(8) function eri_get(i,j,k,l) result(val)
    implicit none
    integer,intent(in) :: i,j,k,l
end function eri_get

module subroutine integrals_finalize()
    implicit none
end subroutine integrals_finalize

module subroutine Normal(iL,di,bas,nBasis,buf1e,Norfac,nConts)
    implicit none
    integer :: iL,di,nConts,nBasis
    real(8) :: buf1e(di,di)
    integer :: bas(8,nbasis)
    real(8) :: Norfac(nConts)
end subroutine Normal

module integer(8) function INDEX_2E(i,j,k,l)
    implicit none
    integer  ::   i,j,k,l
end function INDEX_2E

module subroutine  store2e(shls,bas,buf2e,TWOEI,di,dj,dk,dl,nConts,nRec,nBases,Norfac)
implicit none
integer :: di, dj, dk, dl, nConts, nBases
integer(8) :: nRec
integer :: shls(4)
integer :: bas(8,nBases)
real(8) :: buf2e(di,dj,dk,dl)
real(8) :: TWOEI(nRec)
real(8) :: Norfac(nConts)
end subroutine store2e

module subroutine  store1e(shls,di,dj,bas,nBases,buf1e,Norfac,nConts,S)
implicit none
integer :: di, dj, nConts,  nBases
integer :: shls(4)
integer :: bas(8,nBases)
real(8) :: S(nConts,nConts)
real(8) :: Norfac(nConts)
real(8) :: buf1e(di,dj)
end subroutine store1e

module subroutine  store1edrv(shls,di,dj,ao_offset,nBases,buf1e,Norfac,nConts,S)
implicit none
integer :: di, dj, nConts,  nBases
integer :: shls(2)
integer :: ao_offset(0:nBases-1)
real(8) :: S(nConts,nConts,3)
real(8) :: Norfac(nConts)
real(8) :: buf1e(di,dj,3)
end subroutine store1edrv

module subroutine  store2edrv(shls,ao_offset,buf2e,dJi,dKa,dKb,Da,Db,di,dj,dk,dl,nConts,nBases,Norfac)
implicit none
integer :: di, dj, dk, dl, nConts, nBases
integer :: shls(4)
integer :: ao_offset(0:nBases-1)
real(8) :: buf2e(di,dj,dk,dl,3)
real(8) :: dJi(nConts,nConts,3),dKa(nConts,nConts,3),dKb(nConts,nConts,3)
real(8) :: Da(nConts,nConts),Db(nConts,nConts)
real(8) :: Norfac(nConts)
end subroutine store2edrv

module subroutine compute_shell_density_bound(nConts, Ptot, dmax_shell)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: dmax_shell(0:nBases-1,0:nBases-1)
end subroutine compute_shell_density_bound

module subroutine integrals_build_coulomb(nConts, Ptot, J)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
end subroutine integrals_build_coulomb

module subroutine integrals_build_jk_fused(nConts, Ptot, Da, Db, J, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine integrals_build_jk_fused

module logical function jk_fusion_active() result(ok)
    implicit none
end function jk_fusion_active

module logical function is_canonical_quartet(e1,e2,e3,e4) result(ok)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4
end function is_canonical_quartet

module subroutine accumulate_quartet_J(e1,e2,e3,e4,val,nConts,Ptot,J)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4,nConts
    real(8),intent(in) :: val
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(inout) :: J(nConts,nConts)
end subroutine accumulate_quartet_J

module subroutine expand_quartet_tuples(e1,e2,e3,e4,n_out,idx_out)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4
    integer,intent(out) :: n_out
    integer,intent(out) :: idx_out(4,8)
end subroutine expand_quartet_tuples

module subroutine build_coulomb_pk_list(nConts)
    implicit none
    integer,intent(in) :: nConts
end subroutine build_coulomb_pk_list

module subroutine build_coulomb_store(nConts, Ptot, J)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
end subroutine build_coulomb_store

module subroutine build_coulomb_direct(nConts, Ptot, J, is_incremental)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    logical,intent(in),optional :: is_incremental
end subroutine build_coulomb_direct

module subroutine accumulate_quartet_K(e1,e2,e3,e4,val,nConts,D,K)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4,nConts
    real(8),intent(in) :: val
    real(8),intent(in) :: D(nConts,nConts)
    real(8),intent(inout) :: K(nConts,nConts)
end subroutine accumulate_quartet_K

module subroutine integrals_build_exchange(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine integrals_build_exchange

module subroutine integrals_build_exchange_cosx(nConts, Da, Db, Ka, Kb, need_force)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    logical,intent(in) :: need_force
end subroutine integrals_build_exchange_cosx

module subroutine integrals_build_exchange_cosx_lr(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine integrals_build_exchange_cosx_lr

module subroutine integrals_build_exchange_cosx_sr(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine integrals_build_exchange_cosx_sr

module subroutine integrals_force_exchange_cosx(nConts, natoms, Da, Db, RSomega, forceKa, forceKb)
    implicit none
    integer,intent(in) :: nConts, natoms
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(in) :: RSomega
    real(8),intent(out) :: forceKa(natoms,3), forceKb(natoms,3)
end subroutine integrals_force_exchange_cosx

module subroutine integrals_build_exchange_lr(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine integrals_build_exchange_lr

module subroutine build_exchange_pk_list(nConts)
    implicit none
    integer,intent(in) :: nConts
end subroutine build_exchange_pk_list

module subroutine unpack_pair(pq,p,q)
    implicit none
    integer,intent(in) :: pq
    integer,intent(out) :: p,q
end subroutine unpack_pair

module subroutine build_exchange_store(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine build_exchange_store

module subroutine build_exchange_direct(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine build_exchange_direct

module subroutine ensure_int2e_optimizer_lr()
    implicit none
end subroutine ensure_int2e_optimizer_lr

module subroutine build_exchange_lr_direct(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
end subroutine build_exchange_lr_direct

module subroutine integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(out) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)
    real(8),intent(out) :: dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3)
end subroutine integrals_compute_force

module subroutine integrals_compute_force_exchange_lr(nConts, Pa, Pb, dKa_lr, dKb_lr)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(out) :: dKa_lr(nConts,nConts,3), dKb_lr(nConts,nConts,3)
end subroutine integrals_compute_force_exchange_lr

module subroutine integrals_compute_force_1e(nConts, dS, dHcore)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(out) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)
end subroutine integrals_compute_force_1e

module subroutine integrals_nuc_deriv_allatoms(nConts, natoms_in, Ptot, force_out)
implicit none
integer,intent(in) :: nConts, natoms_in
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: force_out(natoms_in,3)
end subroutine integrals_nuc_deriv_allatoms

module subroutine integrals_nuc_attraction_deriv(coor_bohr, nConts, Dr)
    implicit none
    real(8),intent(in) :: coor_bohr(3)
    integer,intent(in) :: nConts
    real(8),intent(out) :: Dr(nConts,nConts,3)
end subroutine integrals_nuc_attraction_deriv

module subroutine integrals_ecp_force(nConts_in, natoms_in, Ptot, ecpForce)
use MOL_info
implicit none
    integer,intent(in) :: nConts_in, natoms_in
    real(8),intent(in) :: Ptot(nConts_in,nConts_in)
    real(8),intent(out) :: ecpForce(natoms_in,3)
end subroutine integrals_ecp_force

module subroutine integrals_hessian_1e_shell(nConts, s1aa, s1ab, h1aa, h1ab)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(out) :: s1aa(nConts,nConts,9), s1ab(nConts,nConts,9)
    real(8),intent(out) :: h1aa(nConts,nConts,9), h1ab(nConts,nConts,9)
end subroutine integrals_hessian_1e_shell

module subroutine integrals_hessian_nuc_atom(coor_bohr, nConts, CC, CB, CK)
    implicit none
    real(8),intent(in) :: coor_bohr(3)
    integer,intent(in) :: nConts
    real(8),intent(out) :: CC(nConts,nConts,9), CB(nConts,nConts,9), CK(nConts,nConts,9)
end subroutine integrals_hessian_nuc_atom

module subroutine integrals_hessian_nuc_atom_raw(coor_bohr, nConts, AA, AB)
    implicit none
    real(8),intent(in) :: coor_bohr(3)
    integer,intent(in) :: nConts
    real(8),intent(out) :: AA(nConts,nConts,9), AB(nConts,nConts,9)
end subroutine integrals_hessian_nuc_atom_raw

module subroutine integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2, dK2)
    implicit none
    integer,intent(in) :: nConts, iatom
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(out) :: dJ2(nConts,nConts,3), dK2(nConts,nConts,3)
end subroutine integrals_cphf_2e_density_role

module subroutine integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e)
    implicit none
    integer,intent(in) :: nConts, natoms
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(in) :: HF_exchange_frac
    real(8),intent(out) :: d2E_2e(3,natoms,3,natoms)
end subroutine integrals_hessian_2e_direct

module subroutine integrals_hessian_e2e_fixed_density(nConts, Pa, Pb, HF_exchange_frac, Eout)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(in) :: HF_exchange_frac
    real(8),intent(out) :: Eout
end subroutine integrals_hessian_e2e_fixed_density

module subroutine integrals_hessian_debug_check2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac)
    implicit none
    integer,intent(in) :: nConts, natoms
    real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
    real(8),intent(in) :: HF_exchange_frac
end subroutine integrals_hessian_debug_check2e_direct

module subroutine integrals_hessian_e1e_fixed_density(nConts, Ptot, Wtot, Eout)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts), Wtot(nConts,nConts)
    real(8),intent(out) :: Eout
end subroutine integrals_hessian_e1e_fixed_density

module subroutine integrals_hessian_debug_check1e(nConts, natoms, Ptot, Wtot)
    implicit none
    integer,intent(in) :: nConts, natoms
    real(8),intent(in) :: Ptot(nConts,nConts), Wtot(nConts,nConts)
end subroutine integrals_hessian_debug_check1e

module subroutine build_df_aux_basis()
use MOL_info
implicit none
end subroutine build_df_aux_basis

module subroutine build_df_aux_basis_fromfile(auxname)
use MOL_info
implicit none
character(len=*),intent(in) :: auxname
end subroutine build_df_aux_basis_fromfile

module subroutine build_df_integrals()
use MOL_info
use mod_meminfo, only: get_memory_budget_bytes, get_store_threshold_bytes
implicit none
end subroutine build_df_integrals

module subroutine build_df_integrals_lr()
use MOL_info
implicit none
end subroutine build_df_integrals_lr

module subroutine integrals_estimate_aux_size()
implicit none
end subroutine integrals_estimate_aux_size

module subroutine integrals_build_coulomb_df(nConts_in, Ptot, J)
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: J(nConts_in,nConts_in)
end subroutine integrals_build_coulomb_df

module subroutine df_build_coulomb_core(nConts_in, Ptot, J, incremental_call)
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: J(nConts_in,nConts_in)
logical,intent(in),optional :: incremental_call
end subroutine df_build_coulomb_core

module subroutine build_df_triple_list()
implicit none
end subroutine build_df_triple_list

module subroutine df_build_coulomb_direct(nConts_in, Ptot, J, incremental_call)
use MOL_info
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: J(nConts_in,nConts_in)
logical,intent(in),optional :: incremental_call
end subroutine df_build_coulomb_direct

module subroutine ensure_df_built()
implicit none
end subroutine ensure_df_built

module subroutine build_df_metric_invhalf()
implicit none
end subroutine build_df_metric_invhalf

module subroutine build_df_exchange_metric()
implicit none
end subroutine build_df_exchange_metric

module subroutine build_df_exchange_metric_full_lr()
implicit none
end subroutine build_df_exchange_metric_full_lr

module subroutine integrals_build_exchange_df(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
end subroutine integrals_build_exchange_df

module subroutine integrals_build_exchange_df_lr(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
end subroutine integrals_build_exchange_df_lr

module subroutine df_build_exchange_direct(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
use MOL_info
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
end subroutine df_build_exchange_direct

module subroutine df_build_exchange_direct_lr(nConts_in, Ca, nOccA, Cb, nOccB, Ka, Kb)
use MOL_info
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: Ka(nConts_in,nConts_in), Kb(nConts_in,nConts_in)
end subroutine df_build_exchange_direct_lr

module subroutine build_df_minv_full()
implicit none
end subroutine build_df_minv_full

module subroutine ensure_df_lr_optimizer_built()
implicit none
end subroutine ensure_df_lr_optimizer_built

module subroutine integrals_force_df(nConts_in, natoms_in, Pa, Pb, Ca, nOccA, Cb, nOccB, &
                                      force_J, force_Ka, force_Kb, need_k)
implicit none
integer,intent(in) :: nConts_in, natoms_in, nOccA, nOccB
real(8),intent(in) :: Pa(nConts_in,nConts_in), Pb(nConts_in,nConts_in)
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: force_J(natoms_in,3), force_Ka(natoms_in,3), force_Kb(natoms_in,3)
logical,intent(in) :: need_k
end subroutine integrals_force_df

module subroutine integrals_force_df_lr(nConts_in, natoms_in, Ca, nOccA, Cb, nOccB, &
                                         force_Ka_lr, force_Kb_lr)
implicit none
integer,intent(in) :: nConts_in, natoms_in, nOccA, nOccB
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: force_Ka_lr(natoms_in,3), force_Kb_lr(natoms_in,3)
end subroutine integrals_force_df_lr

module subroutine df_force_2c2e_contract(nT, W, natoms_in, force_out)
implicit none
integer,intent(in) :: natoms_in, nT
real(8),intent(in) :: W(nContsAux,nContsAux,nT)
real(8),intent(inout) :: force_out(natoms_in,3,nT)
end subroutine df_force_2c2e_contract

module subroutine integrals_force_df_direct_J(nConts_in, natoms_in, Ptot, force_J)
implicit none
integer,intent(in) :: nConts_in, natoms_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: force_J(natoms_in,3)
end subroutine integrals_force_df_direct_J

module subroutine integrals_force_df_direct_JK(nConts_in, natoms_in, Pa, Pb, Ca, nOccA, Cb, nOccB, &
                                                force_J, force_Ka, force_Kb)
implicit none
integer,intent(in) :: nConts_in, natoms_in, nOccA, nOccB
real(8),intent(in) :: Pa(nConts_in,nConts_in), Pb(nConts_in,nConts_in)
real(8),intent(in) :: Ca(nConts_in,nConts_in), Cb(nConts_in,nConts_in)
real(8),intent(out) :: force_J(natoms_in,3), force_Ka(natoms_in,3), force_Kb(natoms_in,3)
end subroutine integrals_force_df_direct_JK

module subroutine decide_df_force_mode()
implicit none
end subroutine decide_df_force_mode

module subroutine compute_atomic_density(ia, closed_pairing, Pa_atom, Pb_atom, ok)
implicit none
integer,intent(in) :: ia
logical,intent(in) :: closed_pairing
real(8),allocatable,intent(out) :: Pa_atom(:,:), Pb_atom(:,:)
logical,intent(out) :: ok
end subroutine compute_atomic_density

module subroutine calc_properties(Natoms_in, MLcharge_out, RESPcharge_in)
implicit none
integer,intent(in) :: Natoms_in
real(8),intent(out) :: MLcharge_out(Natoms_in)
real(8),intent(in),optional :: RESPcharge_in(Natoms_in)
end subroutine calc_properties

module subroutine resp_charges(Natoms_in, RESPcharge_out)
implicit none
integer,intent(in) :: Natoms_in
real(8),intent(out) :: RESPcharge_out(Natoms_in)
end subroutine resp_charges

module subroutine cosmo_build_one_tess_matrix(tess_pos_bohr, Bk)
implicit none
real(8),intent(in) :: tess_pos_bohr(3)
real(8),intent(out) :: Bk(:,:)
end subroutine cosmo_build_one_tess_matrix

module subroutine esp_at_grid_batch(npts, grid_bohr, Ptot, V_elec)
implicit none
integer,intent(in) :: npts
real(8),intent(in) :: grid_bohr(3,npts)
real(8),intent(in) :: Ptot(:,:)
real(8),intent(out) :: V_elec(npts)
end subroutine esp_at_grid_batch

module subroutine esp_at_grid_batch_df(npts, grid_bohr, Ptot, V_elec, ok)
implicit none
integer,intent(in) :: npts
real(8),intent(in) :: grid_bohr(3,npts)
real(8),intent(in) :: Ptot(:,:)
real(8),intent(out) :: V_elec(npts)
logical,intent(out) :: ok
end subroutine esp_at_grid_batch_df

module subroutine cosmo_build_one_tess_matrix_shellderiv(tess_pos_bohr, dBk)
implicit none
real(8),intent(in) :: tess_pos_bohr(3)
real(8),intent(out) :: dBk(:,:,:)
end subroutine cosmo_build_one_tess_matrix_shellderiv

module integer function cgto_engine(ish, bas) result(n)
implicit none
integer,intent(in) :: ish
integer,intent(in) :: bas(*)
end function cgto_engine

module subroutine ovlp1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine ovlp1e_engine

module subroutine kin1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine kin1e_engine

module subroutine nuc1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine nuc1e_engine

module subroutine r1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine r1e_engine

module subroutine ecp1e_engine(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
implicit none
real(8) :: buf(*)
integer :: dims(2)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
integer(8),intent(in) :: cache
end subroutine ecp1e_engine

module subroutine ovlp1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine ovlp1e_ip_engine

module subroutine kin1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine kin1e_ip_engine

module subroutine nuc1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine nuc1e_ip_engine

module subroutine grids1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine grids1e_ip_engine

module subroutine grids1e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine grids1e_engine

module subroutine grids1e_engine_cached(buf, shls, atm, natm, bas, nbas, env, opt, cache)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
real(8) :: cache(*)
end subroutine grids1e_engine_cached

module subroutine rinv1e_ip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine rinv1e_ip_engine

module subroutine twoe_ip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoe_ip1_engine

module subroutine ovlp1e_ipip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine ovlp1e_ipip_engine

module subroutine ovlp1e_ipovlpip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine ovlp1e_ipovlpip_engine

module subroutine kin1e_ipip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine kin1e_ipip_engine

module subroutine kin1e_ipkinip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine kin1e_ipkinip_engine

module subroutine nuc1e_ipip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine nuc1e_ipip_engine

module subroutine nuc1e_ipnucip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine nuc1e_ipnucip_engine

module subroutine twoe_ipip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoe_ipip1_engine

module subroutine twoe_ip1ip2_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoe_ip1ip2_engine

module subroutine twoe_ipvip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoe_ipvip1_engine

module subroutine twoe_ipip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine twoe_ipip1_optimizer_engine

module subroutine rinv1e_ipiprinv_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine rinv1e_ipiprinv_engine

module subroutine rinv1e_iprinvip_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine rinv1e_iprinvip_engine

module subroutine twoe_ip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine twoe_ip1_optimizer_engine

module subroutine threec2e_ip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine threec2e_ip1_engine

module subroutine threec2e_ip2_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine threec2e_ip2_engine

module subroutine threec2e_ip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine threec2e_ip1_optimizer_engine

module subroutine threec2e_ip2_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine threec2e_ip2_optimizer_engine

module subroutine twoc2e_ip1_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoc2e_ip1_engine

module subroutine twoc2e_ip2_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoc2e_ip2_engine

module subroutine twoc2e_ip1_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine twoc2e_ip1_optimizer_engine

module subroutine twoc2e_ip2_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine twoc2e_ip2_optimizer_engine

module subroutine ecp1e_ipnuc_engine(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
implicit none
real(8) :: buf(*)
integer :: dims(2)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
integer(8),intent(in) :: cache
end subroutine ecp1e_ipnuc_engine

module subroutine ecp1e_iprinv_engine(buf, dims, shls, atm, natm, bas, nbas, env, opt, cache)
implicit none
real(8) :: buf(*)
integer :: dims(2)
integer :: shls(2)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
integer(8),intent(in) :: cache
end subroutine ecp1e_iprinv_engine

module subroutine twoe_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoe_engine

module subroutine threec2e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine threec2e_engine

module subroutine twoc2e_engine(buf, shls, atm, natm, bas, nbas, env, opt)
implicit none
real(8) :: buf(*)
integer :: shls(4)
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
integer(8),intent(in) :: opt
end subroutine twoc2e_engine

module subroutine twoe_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine twoe_optimizer_engine

module subroutine threec2e_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine threec2e_optimizer_engine

module subroutine twoc2e_optimizer_engine(opt, atm, natm, bas, nbas, env)
implicit none
integer(8),intent(out) :: opt
integer :: atm(*)
integer,intent(in) :: natm
integer :: bas(*)
integer,intent(in) :: nbas
real(8) :: env(*)
end subroutine twoc2e_optimizer_engine

end interface

end module mod_integrals
