! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Plain data carrier for a fully-parsed Engine input: everything a

module mod_engine_input_types
implicit none
private
public :: engine_input_t, EI_MAXATOM, spec_freq_pressure_pa, spec_tddft_cube
public :: parse_runtype, runtype_string, parse_polar_omega, hess_every_code, hess_every_str, nto_cube_items
public :: resolve_cosmo_avg_area

integer, parameter :: EI_MAXATOM = 1000

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
   character(128):: polar_omega = ''
   integer       :: tddft_nroots = 0
   logical       :: tddft_robust = .false.
   logical       :: tddft_triplets = .false.
   logical       :: tddft_tda = .false.
   real(8)       :: tddft_eps_opt(1) = [1.78d0]
   logical       :: uks = .false.
   logical       :: jcoupling = .false.
   character(32) :: jcoupling_terms = 'all'
   character(64) :: jcoupling_atoms = 'all'
   logical       :: vv10_nonself = .false.
   real(8)       :: mem_cap_gb = 0.0d0
   logical       :: estimate_only = .false.
   integer       :: n_threads = 0
   integer       :: verbose = 1
   integer       :: scf_conv_level = 0
   integer       :: newton_cg_iter = -1
   integer       :: stability_mode = -1

   logical       :: irc_run = .false.
   logical       :: nmr_run
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
   logical       :: opt_fragment = .false.
   character(256):: opt_hessian_file = ''
   logical       :: opt_ts = .false.
   character(8)  :: opt_hess = 'model'
   integer       :: opt_hess_every = 0
   character(8)  :: opt_hess_recalc = ''
   character(8)  :: opt_xtb_method = 'gfn2'
   character(256):: opt_xtb_cmd = 'xtb'
   character(256):: opt_xtb_args = ''
   logical       :: opt_fragopts(2) = [.false., .true.]
   character(256):: opt_frag_charges = '', opt_frag_mults = ''
   character(512):: opt_constraints = ''
   character(256):: opt_ts_follow = ''

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
   real(8)       :: freq_cond(3) = [298.15d0, 1.0d0, 0.0d0]
   character(96) :: cube = ''
   character(256) :: nics = ''

   logical       :: cosmo_on = .false.
   real(8)       :: cosmo_epsilon = 78.4d0
   real(8)       :: cosmo_radii_scale = 1.2d0
   real(8)       :: cosmo_avg_area = -1.0d0
   real(8)       :: cosmo_sigma_rav = 0.5d0
   character(16) :: cosmo_cavity_type = 'gepol'
   real(8)       :: cosmo_rsolv = 0.0d0
   integer       :: cosmo_ks_nseg = 92
   integer       :: cosmo_ks_nface = 1082
   character(256):: cosmo_sigma_profile_file = ''
   character(32) :: cosmo_solvent = ''
   logical       :: cosmo_smd = .false.
   character(8)  :: cosmo_convention = 'cosmo'

end type engine_input_t

contains

function spec_tddft_cube(spec) result(t)
type(engine_input_t), intent(in) :: spec
character(len=64) :: t
integer :: i, j
t = ''
i = 1
do while (i <= len_trim(spec%cube))
   j = index(spec%cube(i:), ',')
   if (j == 0) j = len_trim(spec%cube) - i + 2
   if (index(spec%cube(i:i+j-2), ':') > 0 .or. trim(adjustl(spec%cube(i:i+j-2))) == 'none' .or. &
       trim(adjustl(spec%cube(i:i+j-2))) == 'ct') then
      if (len_trim(t) > 0) t = trim(t)//','
      t = trim(t)//trim(adjustl(spec%cube(i:i+j-2)))
   endif
   i = i + j
enddo
end function spec_tddft_cube

pure function spec_freq_pressure_pa(spec) result(p)
type(engine_input_t), intent(in) :: spec
real(8) :: p
p = spec%freq_cond(2)*101325.0d0
if (spec%freq_cond(3) > 0.0d0) p = spec%freq_cond(3)*1000.0d0*8.314462618d0*spec%freq_cond(1)
end function spec_freq_pressure_pa

subroutine parse_runtype(str, do_opt, do_freq, do_polar, ok, do_irc, do_nmr)
character(len=*), intent(in)  :: str
logical, intent(out)          :: do_opt, do_freq, do_polar, ok
logical, intent(out), optional :: do_irc, do_nmr
logical :: irc, nmr
character(len=64) :: w
integer :: i, n, c, nw
logical :: sp
do_opt = .false.; do_freq = .false.; do_polar = .false.; ok = .true.; sp = .false.; nw = 0; irc = .false.; nmr = .false.
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
         case ('nmr')
            nmr = .true.
         case default
            ok = .false.
         end select
         w = ''
      endif
   endif
   i = i + 1
enddo
if (nw == 0) sp = .true.
if (sp .and. (do_opt .or. do_freq .or. do_polar .or. irc .or. nmr)) ok = .false.
if (nmr .and. (do_opt .or. do_freq .or. do_polar .or. irc)) ok = .false.
if (do_opt .and. do_polar) ok = .false.
if (irc .and. (do_opt .or. do_freq .or. do_polar)) ok = .false.
if (present(do_nmr)) then
   do_nmr = nmr
else if (nmr) then
   ok = .false.
endif
if (present(do_irc)) then
   do_irc = irc
else if (irc) then
   ok = .false.
endif
end subroutine parse_runtype

function runtype_string(do_opt, do_freq, do_polar, do_irc, do_nmr) result(str)
logical, intent(in) :: do_opt, do_freq, do_polar
logical, intent(in), optional :: do_irc, do_nmr
character(len=24) :: str
str = ''
if (present(do_nmr)) then
   if (do_nmr) then
      str = 'nmr'
      return
   endif
endif
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

subroutine parse_polar_omega(str, omega, n, ok)
character(len=*), intent(in) :: str
real(8), intent(out) :: omega(16)
integer, intent(out) :: n
logical, intent(out) :: ok
real(8), parameter :: HARTREE_EV = 27.211386245988d0, HARTREE_NM = 45.56335252760d0
character(len=len(str)) :: work
character(len=64) :: tok
integer :: i, l, ios, k
real(8) :: v
work = str
do i = 1, len(work)
   if (work(i:i) == ',') work(i:i) = ' '
enddo
n = 0; ok = .true.; omega = 0.0d0
work = adjustl(work)
do while (len_trim(work) > 0)
   k = index(trim(work), ' ')
   if (k == 0) k = len_trim(work) + 1
   tok = work(1:k-1)
   work = adjustl(work(k:))
   l = len_trim(tok)
   if (n >= 16) then
      ok = .false. ; return
   endif
   if (l > 2 .and. (tok(l-1:l) == 'nm' .or. tok(l-1:l) == 'NM')) then
      read(tok(1:l-2), *, iostat=ios) v
      if (ios /= 0 .or. v <= 0.0d0) then
         ok = .false. ; return
      endif
      v = HARTREE_NM / v
   else if (l > 2 .and. (tok(l-1:l) == 'ev' .or. tok(l-1:l) == 'EV' .or. tok(l-1:l) == 'eV')) then
      read(tok(1:l-2), *, iostat=ios) v
      if (ios /= 0) then
         ok = .false. ; return
      endif
      v = v / HARTREE_EV
   else
      read(tok(1:l), *, iostat=ios) v
      if (ios /= 0) then
         ok = .false. ; return
      endif
   endif
   if (v < 0.0d0) then
      ok = .false. ; return
   endif
   n = n + 1
   omega(n) = v
enddo
end subroutine parse_polar_omega

function nto_cube_items(str, with_cube) result(t)
character(len=*), intent(in) :: str
logical, intent(in), optional :: with_cube
character(len=64) :: t
character(len=len(str)) :: w
integer :: i, k, ios, st
t = ''
w = str
do i = 1, len(w)
   if (w(i:i) == ',') w(i:i) = ' '
enddo
i = 1
do while (i <= len_trim(w))
   if (w(i:i) == ' ') then
      i = i + 1 ; cycle
   endif
   k = i
   do while (k <= len_trim(w))
      if (w(k:k) == ' ') exit
      k = k + 1
   enddo
   read(w(i:k-1), *, iostat=ios) st
   if (ios /= 0 .or. st < 1) then
      print *, "Input error: nto = '"//trim(str)//"' (expected state numbers such as '1, 2')"
      stop 1
   endif
   if (len_trim(t) > 0) t = trim(t)//', '
   if (present(with_cube)) then
      if (.not. with_cube) then
         t = trim(t)//'ntom:'//trim(adjustl(itoa(st)))
         i = k
         cycle
      endif
   endif
   t = trim(t)//'nto:'//trim(adjustl(itoa(st)))
   i = k
enddo
end function nto_cube_items

function itoa(i) result(s)
integer, intent(in) :: i
character(len=12) :: s
write(s,'(I0)') i
end function itoa

function hess_every_code(str) result(code)
character(len=*), intent(in) :: str
integer :: code, ios
character(len=len(str)) :: t
t = adjustl(str)
code = 0
if (len_trim(t) == 0) return
select case (t(1:len_trim(t)))
case ('first', 'First', 'FIRST'); code = -1
case ('all', 'All', 'ALL');       code = 1
case default
   read(t, *, iostat=ios) code
   if (ios /= 0 .or. code < 0) then
      print *, "Input error: hess_every = '"//trim(t)//"' (expected first, all or a number of cycles)"
      stop 1
   endif
end select
end function hess_every_code

function hess_every_str(code) result(str)
integer, intent(in) :: code
character(len=8) :: str
if (code < 0) then
   str = 'first'
else if (code == 0) then
   str = '0'
else
   write(str,'(I0)') code
endif
end function hess_every_str

subroutine resolve_cosmo_avg_area(spec)
type(engine_input_t), intent(inout) :: spec
if (spec%cosmo_avg_area <= 0.0d0) then
   spec%cosmo_avg_area = 0.3d0
   if (spec%nmr_run) spec%cosmo_avg_area = 0.1d0
endif
end subroutine resolve_cosmo_avg_area

end module mod_engine_input_types
