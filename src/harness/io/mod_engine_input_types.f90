! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Plain data carrier for a fully-parsed Engine input: everything a

module mod_engine_input_types
implicit none
private
public :: engine_input_t, EI_MAXATOM

integer, parameter :: EI_MAXATOM = 500

type :: engine_input_t
   integer       :: ncenters = 0
   integer       :: imult = 1
   integer       :: icharge = 0
   character(30) :: functional = 'PBE_PBE'
   character(30) :: baselabel = 'def2svp'
   character(30) :: ecplabel = ''
   character(256):: basedir = ''
   character(16) :: J = 'RI'
   character(16) :: K = 'cosx'
   character(64) :: ri_aux_basis = ''
   logical       :: spherical = .true.
   logical       :: harris_guess = .true.
   logical       :: calc_force = .true.
   logical       :: vv10_nonself = .false.
   real(8)       :: mem_cap_gb = 0.0d0
   logical       :: estimate_only = .false.
   integer       :: n_threads = 0
   integer       :: verbose = 1
   integer       :: scf_conv_level = 0

   logical       :: opt_run = .false.
   integer       :: opt_maxcyc = 100
   integer       :: opt_conv_level = 0
   real(8)       :: opt_trust = 0.3d0
   character(16) :: opt_coord = 'ric'
   logical       :: opt_restart = .true.

   integer, allocatable       :: atomchg(:)
   real(8), allocatable       :: x(:), y(:), z(:)
   character(30), allocatable :: atom_basis(:), atom_ecp(:)

   integer :: npc = 0
   real(8), allocatable :: pc_q(:), pc_x(:), pc_y(:), pc_z(:)

   logical        :: chk_read = .false.
   logical        :: chk_write = .false.
   character(256) :: chk_file = ''

   logical        :: molden_write = .false.
   character(256) :: molden_file = ''
   logical        :: molden_read = .false.
   character(256) :: molden_read_file = ''

   logical       :: resp_charges_on = .true.

   logical       :: cosmo_on = .false.
   real(8)       :: cosmo_epsilon = 78.4d0
   real(8)       :: cosmo_radii_scale = 1.2d0
   real(8)       :: cosmo_avg_area = 0.3d0
   real(8)       :: cosmo_sigma_rav = 0.5d0
   character(16) :: cosmo_cavity_type = 'gepol'
   real(8)       :: cosmo_rsolv = 0.0d0
   integer       :: cosmo_ks_nseg = 92
   integer       :: cosmo_ks_nface = 1082
   character(256):: cosmo_sigma_profile_file = ''
   character(32) :: cosmo_solvent = ''
   logical       :: cosmo_smd = .false.

end type engine_input_t

end module mod_engine_input_types
