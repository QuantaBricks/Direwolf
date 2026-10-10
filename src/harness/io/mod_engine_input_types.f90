! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Plain data carrier for a fully-parsed Engine input: everything a

module mod_engine_input_types
implicit none
private
public :: engine_input_t, EI_MAXATOM
public :: parse_runtype, runtype_string

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
   logical       :: freq = .false.
   logical       :: polar = .false.
   logical       :: vv10_nonself = .false.
   real(8)       :: mem_cap_gb = 0.0d0
   logical       :: estimate_only = .false.
   integer       :: n_threads = 0
   integer       :: verbose = 1
   integer       :: scf_conv_level = 0
   integer       :: newton_cg_iter = -1
   integer       :: stability_mode = -1

   logical       :: irc_run = .false.
   character(16) :: irc_direction = 'both'
   real(8)       :: irc_step = 0.10d0
   integer       :: irc_maxpoints = 50
   integer       :: irc_hess_every = 5
   real(8)       :: irc_stop(3) = [3.0d-4, 1.0d-2, 1.0d0]

   logical       :: opt_run = .false.
   integer       :: opt_maxcyc = 100
   integer       :: opt_conv_level = 0
   real(8)       :: opt_gmax = 0.0d0, opt_grms = 0.0d0, opt_smax = 0.0d0, opt_srms = 0.0d0, opt_de = 0.0d0
   real(8)       :: opt_trust = 0.3d0
   character(16) :: opt_coord = 'ric'
   logical       :: opt_restart = .true.
   logical       :: opt_write_ric = .false.
   character(256):: opt_hessian_file = ''
   logical       :: opt_ts = .false.
   integer       :: opt_ts_hess_every = 0
   character(8)  :: opt_hess = 'model'
   integer       :: opt_hess_every = 0
   character(8)  :: opt_hess_recalc = ''
   character(8)  :: opt_xtb_method = 'gfn2'
   character(256):: opt_xtb_cmd = 'xtb'
   character(256):: opt_xtb_args = ''
   character(256):: opt_freeze = ''
   character(512):: opt_constraints = ''

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

contains

subroutine parse_runtype(str, do_opt, do_freq, do_polar, ok, do_irc)
character(len=*), intent(in)  :: str
logical, intent(out)          :: do_opt, do_freq, do_polar, ok
logical, intent(out), optional :: do_irc
logical :: irc
character(len=64) :: w
integer :: i, n, c, nw
logical :: sp
do_opt = .false.; do_freq = .false.; do_polar = .false.; ok = .true.; sp = .false.; nw = 0; irc = .false.
n = len_trim(str); w = ''; i = 1
do while (i <= n + 1)
   c = 0
   if (i <= n) c = ichar(str(i:i))
   if ((c >= ichar('a') .and. c <= ichar('z')) .or. (c >= ichar('A') .and. c <= ichar('Z'))) then
      if (len_trim(w) < len(w)) then
         if (c >= ichar('A') .and. c <= ichar('Z')) c = c + 32
         w(len_trim(w)+1:len_trim(w)+1) = char(c)
      endif
   else
      if (i <= n .and. c /= ichar(',') .and. c /= ichar(' ') .and. c /= 9) ok = .false.
      if (len_trim(w) > 0) then
         nw = nw + 1
         select case (trim(w))
         case ('energy', 'sp')
            sp = .true.
         case ('opt', 'optimize')
            do_opt = .true.
         case ('freq', 'frequency')
            do_freq = .true.
         case ('polar', 'polarizability')
            do_polar = .true.
         case ('irc')
            irc = .true.
         case default
            ok = .false.
         end select
         w = ''
      endif
   endif
   i = i + 1
enddo
if (nw == 0) sp = .true.
if (sp .and. (do_opt .or. do_freq .or. do_polar .or. irc)) ok = .false.
if (do_opt .and. do_polar) ok = .false.
if (irc .and. (do_opt .or. do_freq .or. do_polar)) ok = .false.
if (present(do_irc)) then
   do_irc = irc
else if (irc) then
   ok = .false.
endif
end subroutine parse_runtype

function runtype_string(do_opt, do_freq, do_polar, do_irc) result(str)
logical, intent(in) :: do_opt, do_freq, do_polar
logical, intent(in), optional :: do_irc
character(len=24) :: str
str = ''
if (present(do_irc)) then
   if (do_irc) then
      str = 'irc'
      return
   endif
endif
if (do_opt) str = 'opt'
if (do_freq) then
   if (len_trim(str) > 0) then
      str = trim(str)//',freq'
   else
      str = 'freq'
   endif
endif
if (do_polar) then
   if (len_trim(str) > 0) then
      str = trim(str)//',polar'
   else
      str = 'polar'
   endif
endif
if (len_trim(str) == 0) str = 'energy'
end function runtype_string

end module mod_engine_input_types
