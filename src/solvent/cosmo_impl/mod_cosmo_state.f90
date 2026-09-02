! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_state: the module-level state shared by cosmo_init/

module mod_cosmo_state
implicit none

logical :: cosmo_enabled = .false.

integer :: n_tess = 0
real(8),allocatable :: tess_coor(:,:)
real(8),allocatable :: tess_area(:)
integer,allocatable :: tess_atom(:)
real(8),allocatable :: V_nuc(:)
real(8) :: f_eps = 1.0d0

integer :: saved_natoms = 0
integer,allocatable :: saved_Z(:), saved_ecpCoreElec(:)
real(8) :: saved_radii_scale = 0.0d0, saved_avg_area_ang2 = 0.0d0
real(8) :: saved_Qm = 0.0d0
real(8),allocatable :: saved_Ainv1(:)
real(8) :: saved_sum_Ainv1 = 0.0d0
character(len=16) :: saved_cavity_type = ''
real(8) :: debug_sigma_rav = 0.0d0
character(len=256) :: debug_sigma_profile_file = ''
real(8),allocatable :: debug_prev_hist(:)
logical :: use_smd_radii = .false.
real(8) :: saved_smd_alpha = 0.0d0

real(8),allocatable :: A_fac(:,:)
integer,allocatable :: A_ipiv(:)

logical :: store_mode = .false.
real(8),allocatable :: B_store(:,:,:)
integer :: nConts_saved = 0

real(8),allocatable :: q_last(:)

real(8),allocatable :: Vtot_last(:)

logical :: cosmors_request = .false.
character(256) :: cosmors_cosmo_file_req = ''
character(64) :: cosmors_compound_name_req = ''

end module mod_cosmo_state
