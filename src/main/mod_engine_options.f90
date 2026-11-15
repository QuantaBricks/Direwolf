! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_engine_options: everything a caller asks of EngineUp beyond the molecule, the method and the

module mod_engine_options
use mod_engine_report, only: engine_report_t, report_hook_i
implicit none
private
public :: engine_options_t

type :: engine_options_t
   integer :: verbose = 1
   logical :: quiet = .false.
   logical :: force_dense_allow = .true.
   logical :: cosx_force_loose = .false.
   logical :: cosx_hi_reconverge = .false.
   logical :: resp_charges = .false.
   logical :: esp_charges = .false.
   real(8),allocatable :: esp_q(:)
   logical :: pc_grad_request = .false.
   real(8),allocatable :: pc_grad(:,:)
   logical :: cube_den = .false.
   logical :: cube_esp = .false.
   character(len=256) :: nics = ''
   character(len=64) :: nmr_j_atoms = 'all'
   character(len=32) :: nmr_j_terms = 'all'
   logical :: nmr_j = .false.
   logical :: nmr = .false.
   logical :: cube_every = .false.
   logical :: cube_esp_exact = .false.
   logical :: cosmors_on = .false.
   logical :: cosmo_cpcm = .false.
   character(len=256) :: cosmors_file = ''
   character(len=64)  :: cosmors_compound = ''
   integer :: tddft_nroots = 0
   logical :: tddft_robust = .false.
   character(len=64) :: tddft_cube = ''
   logical :: tddft_triplets = .false.
   logical :: tddft_tda = .false.
   logical :: tddft_print = .true.
   real(8) :: tddft_eps_opt = 1.78d0
   logical :: polar = .false.
   character(len=256) :: polar_omega = ''
   logical :: uks = .false.
   logical :: freq = .false.
   logical :: hessian_return = .false.
   logical :: freq_thermo = .true.
   integer :: freq_unit = 6
   integer :: freq_mult = 1
   real(8) :: freq_temp = 298.15d0
   real(8) :: freq_pres = 101325.0d0
   character(len=256) :: freq_base = 'direwolf'
   logical :: mayer = .false.
   real(8),allocatable :: mayer_bo(:,:)
   integer :: scf_niter = 0
   real(8) :: scf_dP = 0.0d0, scf_dE = 0.0d0
   real(8) :: e_hcore = 0.0d0, e_coul = 0.0d0, e_exx = 0.0d0, e_xc = 0.0d0, e_nuc = 0.0d0
   logical :: vv10_used = .false.
   real(8) :: e_vv10 = 0.0d0
   procedure(report_hook_i), pointer, nopass :: progress => null()
   type(engine_report_t) :: rep
   logical :: done = .false.
   logical :: uks_restart = .false.
   real(8) :: uks_lambda = 0.0d0
end type engine_options_t

public :: engine_output_stem
contains

function engine_output_stem(outfile) result(stem)
character(len=*),intent(in) :: outfile
character(len=256) :: stem
integer :: lo
lo = len_trim(outfile)
stem = outfile
if (lo > 4) then
   if (outfile(lo-3:lo) == '.out') stem = outfile(1:lo-4)
endif
end function engine_output_stem

end module mod_engine_options
