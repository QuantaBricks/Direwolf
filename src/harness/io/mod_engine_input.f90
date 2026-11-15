! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Top-level Engine input-file parser: reads an .inp file's namelist

module mod_engine_input
use mod_engine_input_types, only: engine_input_t, EI_MAXATOM, parse_runtype, hess_every_code, nto_cube_items, resolve_cosmo_avg_area
implicit none
private
public :: engine_read_input

integer, parameter :: MAXLINES = 10000

contains

subroutine engine_read_input(infile, spec)
use mod_engine_input_block, only: find_namelist_block, read_input_lines, strip_key_prefix, rename_keys, quote_value
use mod_engine_input_atoms_namelist, only: parse_atoms_namelist_mode
use mod_engine_input_atoms_xyz, only: read_xyz
use mod_engine_input_atoms_column, only: is_old_namelist_atoms, parse_atoms_column, frozen_ranges, append_freeze
use mod_cosmo, only: cosmo_solvent_epsilon
character(len=*), intent(in)      :: infile
type(engine_input_t), intent(out) :: spec

integer       :: ncenters, imult, icharge
character(30) :: functional, baselabel, ecplabel
character(256):: basedir
character(10) :: unit
character(256):: xyzfile
character(16) :: J, K
character(64) :: ri_aux_basis
character(128):: polar_omega
integer       :: tddft_nroots, n_roots_before
logical       :: tddft_robust, tddft_triplets, tddft_tda, uks, jcoupling
real(8)       :: tddft_eps_opt
integer       :: stability, newton_cg
logical       :: spherical, harris_guess, calc_force, vv10_nonself, freq
real(8)       :: mem_cap_gb
logical       :: estimate_only
integer       :: n_threads, verbose
character(16) :: scf_conv
character(32) :: runtype
logical       :: opt_block, freq_block
integer       :: i_hdr

integer       :: opt_maxcyc
character(16) :: opt_conv, opt_coord
real(8)       :: opt_trust
real(8)       :: opt_gmax, opt_grms, opt_smax, opt_srms, opt_de
logical       :: opt_restart, opt_write_ric, opt_fragment, fragment
character(256):: opt_hessian_file
logical       :: opt_ts
integer       :: opt_ts_hess_every
character(8)  :: opt_hess, opt_hess_recalc, opt_xtb_method
character(8)  :: opt_hess_every
character(256):: opt_xtb_cmd, opt_xtb_args
character(256):: opt_freeze, opt_frag_charges, opt_frag_mults
logical       :: opt_embed, opt_prerelax, frg_embed, frg_prerelax
character(256):: frag_charges, frag_mults
character(512):: opt_constraints
character(256):: opt_ts_follow

character(16) :: irc_direction
real(8)       :: irc_step, irc_gmax, irc_dbo
logical       :: irc_corrector
integer       :: irc_maxpoints, irc_hess_every

integer       :: atomchg(EI_MAXATOM)
real(8)       :: x(EI_MAXATOM), y(EI_MAXATOM), z(EI_MAXATOM)
character(30) :: atom_basis(EI_MAXATOM), atom_ecp(EI_MAXATOM)
logical :: atom_frozen(EI_MAXATOM)

integer :: npc
real(8) :: pc_q(EI_MAXATOM), pc_x(EI_MAXATOM), pc_y(EI_MAXATOM), pc_z(EI_MAXATOM)

logical        :: chk_read, chk_write
character(256) :: chk_file
logical        :: molden_write, molden_read
character(256) :: molden_file, molden_read_file

logical        :: resp_charges_on
character(32)  :: cube
character(256) :: nics
character(32) :: jcoupling_terms
character(64) :: jcoupling_atoms
character(64)  :: tddft_cube, tddft_nto
logical        :: tddft_nto_cube
real(8)        :: freq_temp, freq_pres, freq_conc
logical        :: on
real(8)        :: epsilon, radii_scale, sigma_rav, avg_area
character(16)  :: cavity_type
real(8)        :: rsolv
integer        :: ks_nseg, ks_nface
character(256) :: sigma_profile_file
character(32)  :: solvent
logical        :: smd
character(8)   :: convention

integer,parameter :: UNSET_I = -999999
real(8),parameter :: UNSET_R = -9.999d20
integer :: prof_grid_cache
integer :: prof_xc_shellpair, prof_xc_incr, prof_incremental_fock, prof_diis_window
integer :: prof_force_df_direct, prof_force_df_store, prof_force_df_force_direct
integer :: prof_vv10_radial, prof_vv10_angular
real(8) :: prof_xc_pscreen, prof_level_shift
character(16) :: prof_xc_grid_tier, prof_opt_hess
integer :: prof_p86_vwn
real(8) :: prof_vv10_rcut, prof_vv10_rhocut
integer :: prof_cosx_hi_level, prof_cosx_md_level, prof_cosx_lo_level

character(4096), allocatable :: lines(:)
integer :: nlines, i, na, start, finish
logical :: found
character(len=10000) :: buf
integer :: nrad, nsph

interface
   function c_setenv(name, value, overwrite) bind(C, name="setenv") result(r)
   use iso_c_binding, only: c_char, c_int
   character(kind=c_char), intent(in) :: name(*)
   character(kind=c_char), intent(in) :: value(*)
   integer(c_int), value :: overwrite
   integer(c_int) :: r
   end function
end interface

namelist /molecule/ ncenters, imult, icharge, functional, baselabel, ecplabel, basedir, unit, xyzfile, &
                     J, K, ri_aux_basis, polar_omega, tddft_nroots, tddft_robust, tddft_triplets, tddft_tda, &
                     tddft_eps_opt, uks, jcoupling, &
                     spherical, harris_guess, calc_force, freq, vv10_nonself, &
                     mem_cap_gb, estimate_only, n_threads, verbose, scf_conv, resp_charges_on, cube, tddft_cube, nics, &
                     jcoupling_terms, jcoupling_atoms, &
                     freq_temp, freq_pres, freq_conc, runtype, stability, newton_cg
namelist /opt/ opt_maxcyc, opt_conv, opt_trust, opt_coord, opt_restart, opt_write_ric, &
                opt_hessian_file, opt_ts, opt_ts_hess_every, opt_hess, opt_hess_every, &
                opt_hess_recalc, opt_xtb_method, opt_xtb_cmd, opt_xtb_args, opt_freeze, opt_constraints, opt_ts_follow, &
                opt_embed, opt_frag_charges, opt_frag_mults, opt_prerelax, &
                frg_embed, frg_prerelax, frag_charges, frag_mults, &
                opt_gmax, opt_grms, opt_smax, opt_srms, opt_de, opt_fragment, fragment
namelist /tddft/ tddft_nroots, tddft_robust, tddft_triplets, tddft_tda, tddft_eps_opt, &
                tddft_cube, tddft_nto, tddft_nto_cube
namelist /freq_blk/ freq_temp, freq_pres, &
                freq_conc
namelist /irc/ irc_direction, irc_step, irc_maxpoints, irc_hess_every, irc_gmax, irc_dbo, irc_corrector
namelist /pointcharges/ npc, pc_q, pc_x, pc_y, pc_z
namelist /checkpoint/ chk_read, chk_write, chk_file
namelist /molden/ molden_write, molden_file, molden_read, molden_read_file
namelist /cosmo/ on, epsilon, radii_scale, avg_area, sigma_rav, &
                 cavity_type, rsolv, ks_nseg, ks_nface, sigma_profile_file, &
                 solvent, smd, convention
namelist /professional/ prof_grid_cache, &
                 prof_xc_shellpair, prof_xc_pscreen, prof_xc_incr, &
                 prof_incremental_fock, prof_diis_window, &
                 prof_force_df_direct, prof_force_df_store, prof_force_df_force_direct, &
                 prof_vv10_radial, prof_vv10_angular, prof_vv10_rcut, prof_vv10_rhocut, &
                 prof_cosx_hi_level, prof_cosx_md_level, prof_cosx_lo_level, &
                 prof_xc_grid_tier, prof_level_shift, prof_opt_hess, prof_p86_vwn

call read_input_lines(infile, lines, nlines, MAXLINES)

if (is_gaussian_extension(infile)) then
   print *, 'Error: Gaussian input files (.com/.gjf) are not supported - write the input as a ', &
            'native namelist (&molecule/&atoms) or a .toml: ', trim(infile)
   stop 1
endif

call find_namelist_block(lines, nlines, 'molecule', buf, start, finish, found)
if (.not. found) then
   print *, 'Error: &molecule namelist not found in ', trim(infile), &
            ' (foreign input files are not supported; use the native namelist or a .toml)'
   stop 1
endif

ncenters = 0; imult = 1; icharge = 0
functional = 'PBE_PBE'; baselabel = 'def2svp'; ecplabel = ''; basedir = ''
unit = 'angstrom'; xyzfile = ''
J = 'RI'; K = 'cosx'; ri_aux_basis = ''; polar_omega = ''; tddft_nroots = 0; tddft_robust = .false.
tddft_triplets = .false. ; tddft_tda = .false. ; tddft_eps_opt = 1.78d0 ; uks = .false. ; jcoupling = .false. ; stability = -1
newton_cg = -1
spherical = .true.; harris_guess = .true.; calc_force = .true.; freq = .false.; vv10_nonself = .false.
mem_cap_gb = 0.0d0; estimate_only = .false.; n_threads = 0; verbose = 1; scf_conv = 'regular'
runtype = 'energy'; opt_block = .false.
opt_maxcyc = 100; opt_conv = 'normal'; opt_trust = 0.3d0; opt_coord = 'ric'; opt_restart = .true.
opt_write_ric = .false.; opt_hessian_file = ''; opt_fragment = .false.; fragment = .false.
opt_ts = .false.; opt_ts_hess_every = 0
opt_hess = 'model'; opt_hess_every = ''; opt_hess_recalc = ''; opt_xtb_method = 'gfn2'
opt_xtb_cmd = 'xtb'; opt_xtb_args = ''
opt_freeze = ''; opt_constraints = ''; opt_ts_follow = ''
opt_embed = .false.; opt_frag_charges = ''; opt_frag_mults = ''; opt_prerelax = .true.
frg_embed = .false.; frag_charges = ''; frag_mults = ''; frg_prerelax = .true.
opt_gmax = 0.0d0; opt_grms = 0.0d0; opt_smax = 0.0d0; opt_srms = 0.0d0; opt_de = 0.0d0
irc_direction = 'both'; irc_step = 0.10d0; irc_maxpoints = 50; irc_hess_every = 5; irc_gmax = 3.0d-4
irc_dbo = 1.0d-2; irc_corrector = .true.
atomchg = 0; x = 0.0d0; y = 0.0d0; z = 0.0d0
atom_basis = ''; atom_ecp = ''; atom_frozen = .false.
npc = 0; pc_q = 0.0d0; pc_x = 0.0d0; pc_y = 0.0d0; pc_z = 0.0d0
chk_read = .false.; chk_write = .false.; chk_file = ''
molden_write = .false.; molden_file = ''
molden_read = .false.; molden_read_file = ''
resp_charges_on = .true.
cube = ''; tddft_cube = ''; nics = ''; jcoupling_terms = 'all'; jcoupling_atoms = 'all'
freq_temp = 298.15d0; freq_pres = 1.0d0; freq_conc = 0.0d0
on = .false.; epsilon = 78.4d0; radii_scale = 1.2d0
avg_area = -1.0d0; sigma_rav = 0.5d0
sigma_profile_file = ''
cavity_type = 'gepol'; rsolv = 0.0d0
ks_nseg = 92; ks_nface = 1082
solvent = ''; smd = .false.; convention = 'cosmo'
prof_grid_cache = UNSET_I
prof_xc_shellpair = UNSET_I; prof_xc_incr = UNSET_I; prof_incremental_fock = UNSET_I; prof_diis_window = UNSET_I
prof_force_df_direct = UNSET_I; prof_force_df_store = UNSET_I; prof_force_df_force_direct = UNSET_I
prof_vv10_radial = UNSET_I; prof_vv10_angular = UNSET_I
prof_xc_pscreen = UNSET_R; prof_level_shift = UNSET_R
prof_xc_grid_tier = ''; prof_opt_hess = ''; prof_p86_vwn = UNSET_I
prof_vv10_rcut = UNSET_R; prof_vv10_rhocut = UNSET_R
prof_cosx_hi_level = UNSET_I; prof_cosx_md_level = UNSET_I; prof_cosx_lo_level = UNSET_I

read(buf, nml=molecule)

call find_namelist_block(lines, nlines, 'atoms', buf, start, finish, found)
if (found) then
   if (is_old_namelist_atoms(lines, start, finish)) then
      call parse_atoms_namelist_mode(buf, ncenters, atomchg, x, y, z)
   else
      call parse_atoms_column(lines, start, finish, &
              icharge, imult, unit, baselabel, ecplabel, functional, &
              atomchg, x, y, z, atom_basis, atom_ecp, atom_frozen, na)
      if (ncenters == 0) then
         ncenters = na
      elseif (ncenters /= na) then
         print *, 'Error: ncenters mismatch: ', ncenters, ' vs atoms: ', na
         stop 1
      endif
   endif
endif

call find_namelist_block(lines, nlines, 'opt', buf, start, finish, found)
if (found) then
   call rename_keys(buf, &
     [character(len=24) :: 'max_steps', 'maxcyc', 'convergence', 'conv', 'max_step_size', 'trust', 'coord', 'restart', &
        'write_ric', 'hessian_file', 'ts', 'ts_hess_every', 'hessian', 'hess', 'recomp_hess', 'recompute_hessian', 'hess_every', &
        'hess_recalc', 'xtb_method', 'xtb_cmd', 'xtb_args', 'constraints', 'ts_follow', 'embed', 'prerelax', &
        'max_force', 'gmax', 'rms_force', 'grms', 'max_disp', 'max_displacement', 'smax', 'rms_disp', 'rms_displacement', 'srms', &
        'energy_change', 'de'], &
     [character(len=24) :: 'opt_maxcyc', 'opt_maxcyc', 'opt_conv', 'opt_conv', 'opt_trust', 'opt_trust', 'opt_coord', &
        'opt_restart', 'opt_write_ric', 'opt_hessian_file', 'opt_ts', 'opt_ts_hess_every', 'opt_hess', 'opt_hess', &
        'opt_hess_every', 'opt_hess_every', 'opt_hess_every', 'opt_hess_recalc', 'opt_xtb_method', 'opt_xtb_cmd', 'opt_xtb_args', &
        'opt_constraints', 'opt_ts_follow', 'opt_embed', 'opt_prerelax', 'opt_gmax', 'opt_gmax', 'opt_grms', 'opt_grms', &
        'opt_smax', 'opt_smax', 'opt_smax', 'opt_srms', 'opt_srms', 'opt_srms', 'opt_de', 'opt_de'])
   call quote_value(buf, 'opt_hess_every')
   read(buf, nml=opt)
   opt_block = .true.
endif

call find_namelist_block(lines, nlines, 'tddft', buf, start, finish, found)
if (found) then
   call rename_keys(buf, [character(len=24) :: 'nroots', 'triplets', 'tda', 'robust', 'eps_opt', 'cube', 'nto', &
        'nto_cube'], &
        [character(len=24) :: 'tddft_nroots', 'tddft_triplets', 'tddft_tda', 'tddft_robust', 'tddft_eps_opt', &
        'tddft_cube', 'tddft_nto', 'tddft_nto_cube'])
   tddft_nto = '' ; tddft_nto_cube = .true.
   n_roots_before = tddft_nroots
   tddft_nroots = -1
   read(buf, nml=tddft)
   if (tddft_nroots < 0) tddft_nroots = merge(n_roots_before, 5, n_roots_before > 0)
   if (len_trim(tddft_nto) > 0) then
      if (len_trim(tddft_cube) > 0) tddft_cube = trim(tddft_cube)//', '
      tddft_cube = trim(tddft_cube)//nto_cube_items(tddft_nto, tddft_nto_cube)
   endif
endif

call find_namelist_block(lines, nlines, 'freq', buf, start, finish, found)
freq_block = found
if (found) then
   i_hdr = index(buf, '&freq')
   buf = buf(:i_hdr+4)//'_blk'//buf(i_hdr+5:)
   call rename_keys(buf, [character(len=24) :: 'temperature', 'pressure', 'concentration'], &
        [character(len=24) :: 'freq_temp', 'freq_pres', 'freq_conc'])
   read(buf, nml=freq_blk)
endif

call find_namelist_block(lines, nlines, 'irc', buf, start, finish, found)
if (found) then
   call rename_keys(buf, [character(len=24) :: 'direction', 'step', 'max_points', 'recomp_hess', 'recompute_hessian', 'max_force', &
        'bond_order_change', 'corrector'], &
        [character(len=24) :: 'irc_direction', 'irc_step', 'irc_maxpoints', 'irc_hess_every', 'irc_hess_every', &
        'irc_gmax', 'irc_dbo', &
        'irc_corrector'])
   read(buf, nml=irc)
endif

call find_namelist_block(lines, nlines, 'pointcharges', buf, start, finish, found)
if (found) read(buf, nml=pointcharges)

call find_namelist_block(lines, nlines, 'checkpoint', buf, start, finish, found)
if (found) read(buf, nml=checkpoint)
if (chk_write .and. trim(chk_file) == '') chk_file = trim(make_outfile_name(infile))
if (chk_write) then
   i = len_trim(chk_file)
   if (i > 4) then
      if (chk_file(i-3:i) == '.out') chk_file = chk_file(1:i-4) // '.ckp'
   endif
endif

call find_namelist_block(lines, nlines, 'cosmo', buf, start, finish, found)
if (found) then
   on = .true.
   call strip_key_prefix(buf, 'cosmo_', [character(len=24) :: 'on', 'epsilon', 'radii_scale', 'avg_area', 'sigma_rav', &
        'cavity_type', 'rsolv', 'ks_nseg', 'ks_nface', 'sigma_profile_file', 'solvent', 'smd'])
   read(buf, nml=cosmo)
endif
if (trim(solvent) /= '') epsilon = cosmo_solvent_epsilon(trim(solvent))

call find_namelist_block(lines, nlines, 'molden', buf, start, finish, found)
if (found) read(buf, nml=molden)
if (molden_write .and. trim(molden_file) == '') molden_file = trim(make_outfile_name(infile))
if (molden_write) then
   i = len_trim(molden_file)
   if (i > 4) then
      if (molden_file(i-3:i) == '.out') molden_file = molden_file(1:i-4) // '.molden'
   endif
endif

call find_namelist_block(lines, nlines, 'professional', buf, start, finish, found)
if (found) read(buf, nml=professional)
if (prof_grid_cache /= UNSET_I)         call setenv_int("ENGINE_GRID_CACHE", prof_grid_cache)
if (prof_xc_shellpair /= UNSET_I)       call setenv_int("ENGINE_XC_SHELLPAIR", prof_xc_shellpair)
if (prof_xc_pscreen /= UNSET_R)         call setenv_real("ENGINE_XC_PSCREEN", prof_xc_pscreen)
if (prof_xc_incr /= UNSET_I)            call setenv_int("ENGINE_XC_INCR", prof_xc_incr)
if (len_trim(prof_xc_grid_tier) > 0)    call setenv_str("ENGINE_XC_THEORY_TIER", prof_xc_grid_tier)
if (prof_level_shift /= UNSET_R)        call setenv_real("ENGINE_LEVEL_SHIFT", prof_level_shift)
if (len_trim(prof_opt_hess) > 0)        call setenv_str("ENGINE_OPT_HESS", prof_opt_hess)
if (prof_p86_vwn /= UNSET_I)            call setenv_int("ENGINE_P86_VWN", prof_p86_vwn)
if (prof_incremental_fock /= UNSET_I)   call setenv_int("ENGINE_INCREMENTAL_FOCK", prof_incremental_fock)
if (prof_diis_window /= UNSET_I)        call setenv_int("ENGINE_DIIS_WINDOW", prof_diis_window)
if (prof_force_df_direct /= UNSET_I)    call setenv_int("ENGINE_FORCE_DF_DIRECT", prof_force_df_direct)
if (prof_force_df_store /= UNSET_I)     call setenv_int("ENGINE_FORCE_DF_STORE", prof_force_df_store)
if (prof_force_df_force_direct /= UNSET_I) call setenv_int("ENGINE_FORCE_DF_FORCE_DIRECT", prof_force_df_force_direct)
if (prof_vv10_radial /= UNSET_I)        call setenv_int("ENGINE_VV10_RADIAL", prof_vv10_radial)
if (prof_vv10_angular /= UNSET_I)       call setenv_int("ENGINE_VV10_ANGULAR", prof_vv10_angular)
if (prof_vv10_rcut /= UNSET_R)          call setenv_real("ENGINE_VV10_RCUT", prof_vv10_rcut)
if (prof_vv10_rhocut /= UNSET_R)        call setenv_real("ENGINE_VV10_RHOCUT", prof_vv10_rhocut)
if (prof_cosx_hi_level /= UNSET_I .and. prof_cosx_hi_level /= 2) then
   call cosx_level_to_nrad_nsph(prof_cosx_hi_level, 18, 194, nrad, nsph)
   call setenv_int("ENGINE_COSX_NRAD", nrad)
   call setenv_int("ENGINE_COSX_NSPH", nsph)
endif
if (prof_cosx_md_level /= UNSET_I .and. prof_cosx_md_level /= 2) then
   call cosx_level_to_nrad_nsph(prof_cosx_md_level, 14, 110, nrad, nsph)
   call setenv_int("ENGINE_COSX_NRAD_MD", nrad)
   call setenv_int("ENGINE_COSX_NSPH_MD", nsph)
endif
if (prof_cosx_lo_level /= UNSET_I .and. prof_cosx_lo_level /= 2) then
   call cosx_level_to_nrad_nsph(prof_cosx_lo_level, 14, 50, nrad, nsph)
   call setenv_int("ENGINE_COSX_NRAD_LO", nrad)
   call setenv_int("ENGINE_COSX_NSPH_LO", nsph)
endif

if (xyzfile /= '') then
   block
      logical :: xex
      integer :: ls
      inquire(file=trim(xyzfile), exist=xex)
      if (.not. xex .and. xyzfile(1:1) /= '/') then
         ls = index(infile, '/', back=.true.)
         if (ls > 0) then
            inquire(file=infile(1:ls)//trim(xyzfile), exist=xex)
            if (xex) xyzfile = infile(1:ls)//trim(xyzfile)
         endif
      endif
   end block
   call read_xyz(xyzfile, atomchg, x, y, z, na, atom_frozen)
   if (ncenters == 0) then
      ncenters = na
   elseif (ncenters /= na) then
      print *, 'Error: ncenters mismatch: ', ncenters, ' vs XYZ file: ', na
      stop 1
   endif
endif

if (ncenters <= 0 .or. ncenters > EI_MAXATOM) then
   print *, 'Invalid ncenters: ', ncenters
   stop 1
endif

call convert_unit(atomchg, x, y, z, ncenters, unit)

spec%ncenters = ncenters
spec%imult = imult
spec%icharge = icharge
spec%functional = functional
spec%baselabel = baselabel
block
   character(len=30) :: fu
   integer :: kk, cc, np
   fu = '' ; np = 0
   do kk = 1, len_trim(functional)
      cc = ichar(functional(kk:kk))
      if (functional(kk:kk) == '-' .or. functional(kk:kk) == '_' .or. &
          functional(kk:kk) == ' ') cycle
      if (cc >= ichar('a') .and. cc <= ichar('z')) cc = cc - 32
      np = np + 1
      fu(np:np) = char(cc)
   enddo
   select case (trim(fu))
   case ('B973C');    spec%baselabel = 'mtzvp'
   case ('R2SCAN3C'); spec%baselabel = 'mtzvpp'
   case ('WB97X3C');  spec%baselabel = 'vdzp'
   end select
   if (trim(spec%baselabel) /= trim(baselabel)) &
      write(*,'(A)') ' [Notice] '//trim(functional)//' is a fixed-basis composite'// &
         ' method - using its required basis "'//trim(spec%baselabel)//'"'
end block
spec%ecplabel = ecplabel
spec%basedir = basedir
spec%J = J
spec%K = K
spec%ri_aux_basis = ri_aux_basis
spec%polar_omega = polar_omega
spec%tddft_nroots = tddft_nroots
spec%tddft_robust = tddft_robust
spec%tddft_triplets = tddft_triplets
spec%tddft_tda = tddft_tda
spec%tddft_eps_opt(1) = tddft_eps_opt
spec%uks = uks
spec%jcoupling = jcoupling
spec%stability_mode = stability
spec%newton_cg_iter = newton_cg
spec%spherical = spherical
spec%harris_guess = harris_guess
spec%calc_force = calc_force
spec%vv10_nonself = vv10_nonself
spec%mem_cap_gb = mem_cap_gb
spec%estimate_only = estimate_only
spec%n_threads = n_threads
spec%verbose = verbose
block
   character(16) :: sc
   sc = adjustl(scf_conv)
   select case (trim(sc))
   case ('regular', 'Regular', 'REGULAR', '')
      spec%scf_conv_level = 0
   case ('fine', 'Fine', 'FINE')
      spec%scf_conv_level = 1
   case ('tight', 'Tight', 'TIGHT')
      spec%scf_conv_level = 2
   case default
      write(*,'(A)') "Input error: scf_conv = '"//trim(sc)// &
           "' not recognized (use regular / fine / tight)"
      stop 1
   end select
end block

block
   character(16) :: oc
   logical :: rt_opt, rt_freq, rt_polar, rt_ok, rt_irc, rt_nmr
   call parse_runtype(runtype, rt_opt, rt_freq, rt_polar, rt_ok, rt_irc, rt_nmr)
   if (.not. rt_ok .and. rt_opt .and. rt_polar) then
      write(*,'(A)') "Input error: runtype = '"//trim(runtype)//"': polar cannot be combined with opt "// &
                     "(optimize first, then run polar)"
      stop 1
   endif
   if (.not. rt_ok) then
      write(*,'(A)') "Input error: runtype = '"//trim(runtype)// &
           "' not recognized (use energy, opt, freq, polar, or a comma-separated list such as 'opt, freq')"
      stop 1
   endif
   spec%opt_run = (rt_opt .or. opt_block) .and. .not. rt_irc
   spec%irc_run = rt_irc
   spec%nmr_run = rt_nmr
   if (opt_block .and. .not. rt_opt) &
      write(*,'(A)') " [Notice] an &opt block is given: optimizing the geometry " // &
                     "(add opt to runtype to say so explicitly)"
   spec%freq = freq .or. rt_freq .or. (freq_block .and. .not. rt_irc .and. .not. rt_nmr)
   spec%polar = rt_polar
   oc = adjustl(opt_conv)
   select case (trim(oc))
   case ('normal', 'Normal', 'NORMAL', 'regular', '')
      spec%opt_conv_level = 0
   case ('tight', 'Tight', 'TIGHT')
      spec%opt_conv_level = 1
   case ('loose', 'Loose', 'LOOSE')
      spec%opt_conv_level = -1
   case default
      write(*,'(A)') "Input error: opt_conv = '"//trim(oc)// &
           "' not recognized (use loose / normal / tight)"
      stop 1
   end select
end block
spec%opt_maxcyc  = opt_maxcyc
spec%opt_trust   = opt_trust
spec%opt_coord   = opt_coord
spec%opt_restart = opt_restart
spec%opt_write_ric = opt_write_ric
spec%opt_fragment = opt_fragment .or. fragment
spec%opt_hessian_file = opt_hessian_file
spec%opt_ts = opt_ts
spec%opt_hess = opt_hess
spec%opt_hess_every = hess_every_code(opt_hess_every)
if (spec%opt_hess_every == 0) spec%opt_hess_every = opt_ts_hess_every
spec%opt_hess_recalc = opt_hess_recalc
spec%opt_xtb_method = opt_xtb_method
spec%opt_xtb_cmd = opt_xtb_cmd
spec%opt_xtb_args = opt_xtb_args
spec%opt_fragopts(1) = opt_embed .or. frg_embed
spec%opt_fragopts(2) = opt_prerelax .and. frg_prerelax
spec%opt_frag_charges = opt_frag_charges
if (len_trim(frag_charges) > 0) spec%opt_frag_charges = frag_charges
spec%opt_frag_mults = opt_frag_mults
if (len_trim(frag_mults) > 0) spec%opt_frag_mults = frag_mults
spec%opt_gmax = opt_gmax; spec%opt_grms = opt_grms; spec%opt_smax = opt_smax
spec%opt_srms = opt_srms; spec%opt_de = opt_de
spec%opt_constraints = opt_constraints
if (len_trim(opt_freeze) > 0) call append_freeze(spec%opt_constraints, opt_freeze)
if (any(atom_frozen(1:ncenters))) call append_freeze(spec%opt_constraints, frozen_ranges(atom_frozen, ncenters))
spec%opt_ts_follow = opt_ts_follow
spec%irc_direction = irc_direction; spec%irc_step = irc_step; spec%irc_maxpoints = irc_maxpoints
spec%irc_hess_every = irc_hess_every; spec%irc_stop = [irc_gmax, irc_dbo, merge(1.0d0, 0.0d0, irc_corrector)]

allocate(spec%atomchg(ncenters), spec%x(ncenters), spec%y(ncenters), spec%z(ncenters))
allocate(spec%atom_basis(ncenters), spec%atom_ecp(ncenters))
spec%atomchg = atomchg(1:ncenters)
spec%x = x(1:ncenters); spec%y = y(1:ncenters); spec%z = z(1:ncenters)
spec%atom_basis = atom_basis(1:ncenters)
spec%atom_ecp = atom_ecp(1:ncenters)

spec%npc = npc
allocate(spec%pc_q(max(npc,1)), spec%pc_x(max(npc,1)), spec%pc_y(max(npc,1)), spec%pc_z(max(npc,1)))
spec%pc_q(1:max(npc,1)) = pc_q(1:max(npc,1))
spec%pc_x(1:max(npc,1)) = pc_x(1:max(npc,1))
spec%pc_y(1:max(npc,1)) = pc_y(1:max(npc,1))
spec%pc_z(1:max(npc,1)) = pc_z(1:max(npc,1))

spec%chk_read = chk_read; spec%chk_write = chk_write; spec%chk_file = chk_file
spec%molden_write = molden_write; spec%molden_file = molden_file
spec%molden_read = molden_read; spec%molden_read_file = molden_read_file

spec%resp_charges_on = resp_charges_on
spec%cube = cube
spec%nics = nics
spec%jcoupling_terms = jcoupling_terms
spec%jcoupling_atoms = jcoupling_atoms
if (len_trim(tddft_cube) > 0) then
   spec%cube = trim(cube)
   if (len_trim(spec%cube) > 0) spec%cube = trim(spec%cube)//', '
   spec%cube = trim(spec%cube)//trim(tddft_cube)
endif
spec%freq_cond(1) = freq_temp; spec%freq_cond(2) = freq_pres; spec%freq_cond(3) = freq_conc
spec%cosmo_on = on; spec%cosmo_epsilon = epsilon
spec%cosmo_radii_scale = radii_scale; spec%cosmo_avg_area = avg_area
call resolve_cosmo_avg_area(spec)
spec%cosmo_sigma_rav = sigma_rav; spec%cosmo_cavity_type = cavity_type
spec%cosmo_rsolv = rsolv; spec%cosmo_ks_nseg = ks_nseg
spec%cosmo_ks_nface = ks_nface; spec%cosmo_sigma_profile_file = sigma_profile_file
spec%cosmo_solvent = solvent; spec%cosmo_smd = smd; spec%cosmo_convention = convention

call print_input_summary(infile, spec%ncenters, spec%atomchg, spec%x, spec%y, spec%z, spec%atom_basis, spec%atom_ecp, &
     spec%baselabel, spec%ecplabel, spec%functional, spec%J, spec%K, spec%icharge, spec%imult &
     )

contains

   subroutine setenv_str(name, val)
   use iso_c_binding, only: c_null_char, c_int
   character(len=*), intent(in) :: name, val
   integer(c_int) :: rc
   rc = c_setenv(trim(name)//c_null_char, trim(val)//c_null_char, 1_c_int)
   end subroutine

   subroutine setenv_int(name, val)
   character(len=*), intent(in) :: name
   integer, intent(in) :: val
   character(len=32) :: valstr
   write(valstr,'(I0)') val
   call setenv_str(name, trim(adjustl(valstr)))
   end subroutine

   subroutine setenv_real(name, val)
   character(len=*), intent(in) :: name
   real(8), intent(in) :: val
   character(len=32) :: valstr
   write(valstr,'(ES16.8)') val
   call setenv_str(name, trim(adjustl(valstr)))
   end subroutine

   subroutine cosx_level_to_nrad_nsph(level, def_nrad, def_nsph, out_nrad, out_nsph)
   integer, intent(in)  :: level, def_nrad, def_nsph
   integer, intent(out) :: out_nrad, out_nsph
   integer :: sph_ladder(11), pos, i
   sph_ladder = [26, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302]
   pos = 0
   do i = 1, 11
      if (sph_ladder(i) == def_nsph) pos = i
   enddo
   if (pos == 0) pos = 5
   select case (level)
   case (1)
      out_nrad = max(6, def_nrad - 4)
      out_nsph = sph_ladder(max(1, pos-1))
   case (3)
      out_nrad = def_nrad + 4
      out_nsph = sph_ladder(min(11, pos+1))
   case default
      out_nrad = def_nrad
      out_nsph = def_nsph
   end select
   end subroutine

end subroutine engine_read_input

subroutine print_input_summary(infile, ncenters, atomchg, x, y, z, atom_basis, atom_ecp, baselabel, ecplabel, &
                               functional, jmode, kmode, icharge, imult &
                               )
use mod_engine_input_elements, only: elem_to_name
character(len=*), intent(in) :: infile, baselabel, ecplabel, functional, jmode, kmode
integer, intent(in) :: ncenters, atomchg(*), icharge, imult
real(8), intent(in) :: x(*), y(*), z(*)
character(len=*), intent(in) :: atom_basis(*), atom_ecp(*)
integer :: i

print *, '[Input] ', trim(infile)
print *, '[HEADEREND]'
print *
print *, '[SYSTEM]'
print '("  Atoms: ", I0, "   (coordinates in Angstrom)")', ncenters
do i = 1, ncenters
   print '("    ", A2, 3F14.6, "   basis=", A, "   ecp=", A)', elem_to_name(atomchg(i)), &
      x(i), y(i), z(i), trim(atom_label(atom_basis(i), baselabel)), &
      trim(atom_ecp_label(atom_ecp(i), ecplabel))
enddo
print *, '  Functional: ', trim(functional), '  Basis: ', trim(baselabel)
if (ecplabel /= '') print *, '  ECP: ', trim(ecplabel)
print '(A)', '  Method: '//trim(jk_method_label(jmode, kmode))// &
             ' (J='//trim(jmode)//', K='//trim(kmode)//')'
print *, '  Charge: ', icharge, '  Multiplicity: ', imult

contains

   character(len=32) function jk_method_label(Jm, Km)
   character(len=*), intent(in) :: Jm, Km
   character(len=16) :: ju, ku
   ju = upper(Jm); ku = upper(Km)
   if (trim(ju) == 'RI' .and. trim(ku) == 'COSX') then
      jk_method_label = 'RIJCOSX'
   else if (trim(ju) == 'RI' .and. trim(ku) == 'RI') then
      jk_method_label = 'RIJK (density fitting)'
   else if (trim(ju) == 'RI' .and. trim(ku) == 'EXACT') then
      jk_method_label = 'RIJONX (RI-J + exact K)'
   else if (trim(ju) == 'EXACT' .and. trim(ku) == 'EXACT') then
      jk_method_label = 'conventional (exact J/K)'
   else
      jk_method_label = 'J='//trim(ju)//' K='//trim(ku)
   endif
   end function jk_method_label

   character(len=16) function upper(s)
   character(len=*), intent(in) :: s
   integer :: k, ic
   upper = s
   do k = 1, len_trim(s)
      ic = ichar(upper(k:k))
      if (ic >= ichar('a') .and. ic <= ichar('z')) upper(k:k) = char(ic - 32)
   enddo
   end function upper

   character(len=40) function atom_label(override, global)
   character(len=*), intent(in) :: override, global
   if (override /= '') then
      atom_label = trim(override)
   else
      atom_label = '(global: '//trim(global)//')'
   endif
   end function

   character(len=40) function atom_ecp_label(override, global)
   character(len=*), intent(in) :: override, global
   if (override /= '') then
      atom_ecp_label = trim(override)
   elseif (global /= '') then
      atom_ecp_label = '(global: '//trim(global)//')'
   else
      atom_ecp_label = '(= basis, if any)'
   endif
   end function

end subroutine print_input_summary

logical function is_gaussian_extension(name)
character(len=*), intent(in) :: name
character(len=8) :: ext
integer :: l, i, ic
l = len_trim(name)
is_gaussian_extension = .false.
if (l > 4) then
   ext = name(l-3:l)
   do i = 1, 4
      ic = ichar(ext(i:i))
      if (ic >= ichar('A') .and. ic <= ichar('Z')) ext(i:i) = char(ic + 32)
   enddo
   is_gaussian_extension = (ext == '.com' .or. ext == '.gjf')
endif
end function

logical function has_extension(name, ext)
character(len=*), intent(in) :: name, ext
integer :: l, le, i, ic
character(len=16) :: tail
has_extension = .false.
l = len_trim(name); le = len_trim(ext)
if (le < 1 .or. le > 16 .or. l <= le) return
tail = name(l-le+1:l)
do i = 1, le
   ic = ichar(tail(i:i))
   if (ic >= ichar('A') .and. ic <= ichar('Z')) tail(i:i) = char(ic + 32)
enddo
has_extension = (tail(1:le) == ext(1:le))
end function

function make_outfile_name(name) result(outname)
character(len=*), intent(in) :: name
character(len=256) :: outname
integer :: l
l = len_trim(name)
if (l > 4) then
   if (name(l-3:l) == '.inp') then
      outname = name(1:l-4) // '.out'
      return
   endif
endif
outname = trim(name) // '.out'
end function

subroutine convert_unit(atomchg, x, y, z, ncenters, unit)
integer, intent(inout) :: atomchg(:)
real(8), intent(inout) :: x(:), y(:), z(:)
integer, intent(in) :: ncenters
character(len=*), intent(in) :: unit

real(8), parameter :: ans2bohr = 1.88972612456506d0
integer :: i
if (unit == 'bohr' .or. unit == 'Bohr' .or. unit == 'BOHR') then
   do i = 1, ncenters
      x(i) = x(i) / ans2bohr
      y(i) = y(i) / ans2bohr
      z(i) = z(i) / ans2bohr
   enddo
   print *, '[Input] Coordinates converted from Bohr to Angstrom'
endif
end subroutine

end module mod_engine_input
