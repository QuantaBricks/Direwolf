! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Top-level Engine input-file parser: reads an .inp file's namelist

module mod_engine_input
use mod_engine_input_types, only: engine_input_t, EI_MAXATOM
implicit none
private
public :: engine_read_input

integer, parameter :: MAXLINES = 10000

contains

subroutine engine_read_input(infile, spec)
use mod_engine_input_block, only: find_namelist_block, read_input_lines
use mod_engine_input_atoms_namelist, only: parse_atoms_namelist_mode
use mod_engine_input_atoms_xyz, only: read_xyz
use mod_engine_input_atoms_column, only: is_old_namelist_atoms, parse_atoms_column
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
logical       :: spherical, harris_guess, calc_force, vv10_nonself
real(8)       :: mem_cap_gb
logical       :: estimate_only
integer       :: n_threads, verbose
character(16) :: scf_conv
character(16) :: runtype

integer       :: opt_maxcyc
character(16) :: opt_conv, opt_coord
real(8)       :: opt_trust
logical       :: opt_restart

integer       :: atomchg(EI_MAXATOM)
real(8)       :: x(EI_MAXATOM), y(EI_MAXATOM), z(EI_MAXATOM)
character(30) :: atom_basis(EI_MAXATOM), atom_ecp(EI_MAXATOM)

integer :: npc
real(8) :: pc_q(EI_MAXATOM), pc_x(EI_MAXATOM), pc_y(EI_MAXATOM), pc_z(EI_MAXATOM)

logical        :: chk_read, chk_write
character(256) :: chk_file
logical        :: molden_write, molden_read
character(256) :: molden_file, molden_read_file

logical        :: resp_charges_on
logical        :: cosmo_on
real(8)        :: cosmo_epsilon, cosmo_radii_scale, cosmo_sigma_rav, cosmo_avg_area
character(16)  :: cosmo_cavity_type
real(8)        :: cosmo_rsolv
integer        :: cosmo_ks_nseg, cosmo_ks_nface
character(256) :: cosmo_sigma_profile_file
character(32)  :: cosmo_solvent
logical        :: cosmo_smd

integer,parameter :: UNSET_I = -999999
real(8),parameter :: UNSET_R = -9.999d20
integer :: prof_grid_cache, prof_xc_dyngrid, prof_xc_dyngrid_sph
integer :: prof_xc_shellpair, prof_xc_incr, prof_incremental_fock, prof_diis_window
integer :: prof_force_df_direct, prof_force_df_store, prof_force_df_force_direct
integer :: prof_vv10_radial, prof_vv10_angular
integer :: prof_xc_grid_level
integer :: prof_xc_grid_level_coarse, prof_xc_grid_level_fine
real(8) :: prof_xc_dyngrid_prms, prof_xc_dyngrid_rscale, prof_xc_pscreen
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
                     J, K, ri_aux_basis, spherical, harris_guess, calc_force, vv10_nonself, &
                     mem_cap_gb, estimate_only, n_threads, verbose, scf_conv, resp_charges_on, runtype
namelist /opt/ opt_maxcyc, opt_conv, opt_trust, opt_coord, opt_restart
namelist /pointcharges/ npc, pc_q, pc_x, pc_y, pc_z
namelist /checkpoint/ chk_read, chk_write, chk_file
namelist /molden/ molden_write, molden_file, molden_read, molden_read_file
namelist /cosmo/ cosmo_on, cosmo_epsilon, cosmo_radii_scale, cosmo_avg_area, cosmo_sigma_rav, &
                 cosmo_cavity_type, cosmo_rsolv, cosmo_ks_nseg, cosmo_ks_nface, cosmo_sigma_profile_file, &
                 cosmo_solvent, cosmo_smd
namelist /professional/ prof_grid_cache, prof_xc_dyngrid, prof_xc_dyngrid_prms, prof_xc_dyngrid_rscale, &
                 prof_xc_dyngrid_sph, prof_xc_shellpair, prof_xc_pscreen, prof_xc_incr, &
                 prof_incremental_fock, prof_diis_window, &
                 prof_force_df_direct, prof_force_df_store, prof_force_df_force_direct, &
                 prof_vv10_radial, prof_vv10_angular, prof_vv10_rcut, prof_vv10_rhocut, &
                 prof_xc_grid_level, prof_xc_grid_level_coarse, prof_xc_grid_level_fine, &
                 prof_cosx_hi_level, prof_cosx_md_level, prof_cosx_lo_level

call read_input_lines(infile, lines, nlines, MAXLINES)

if (is_gaussian_extension(infile)) then
endif

call find_namelist_block(lines, nlines, 'molecule', buf, start, finish, found)
if (.not. found) then
   print *, 'Error: &molecule namelist not found in ', trim(infile)
   stop 1
endif

ncenters = 0; imult = 1; icharge = 0
functional = 'PBE_PBE'; baselabel = 'def2svp'; ecplabel = ''; basedir = ''
unit = 'angstrom'; xyzfile = ''
J = 'RI'; K = 'cosx'; ri_aux_basis = ''
spherical = .true.; harris_guess = .true.; calc_force = .true.; vv10_nonself = .false.
mem_cap_gb = 0.0d0; estimate_only = .false.; n_threads = 0; verbose = 1; scf_conv = 'regular'
runtype = 'energy'
opt_maxcyc = 100; opt_conv = 'normal'; opt_trust = 0.3d0; opt_coord = 'ric'; opt_restart = .true.
atomchg = 0; x = 0.0d0; y = 0.0d0; z = 0.0d0
atom_basis = ''; atom_ecp = ''
npc = 0; pc_q = 0.0d0; pc_x = 0.0d0; pc_y = 0.0d0; pc_z = 0.0d0
chk_read = .false.; chk_write = .false.; chk_file = ''
molden_write = .false.; molden_file = ''
molden_read = .false.; molden_read_file = ''
resp_charges_on = .true.
cosmo_on = .false.; cosmo_epsilon = 78.4d0; cosmo_radii_scale = 1.2d0
cosmo_avg_area = 0.3d0; cosmo_sigma_rav = 0.5d0
cosmo_sigma_profile_file = ''
cosmo_cavity_type = 'gepol'; cosmo_rsolv = 0.0d0
cosmo_ks_nseg = 92; cosmo_ks_nface = 1082
cosmo_solvent = ''; cosmo_smd = .false.
prof_grid_cache = UNSET_I; prof_xc_dyngrid = UNSET_I; prof_xc_dyngrid_sph = UNSET_I
prof_xc_shellpair = UNSET_I; prof_xc_incr = UNSET_I; prof_incremental_fock = UNSET_I; prof_diis_window = UNSET_I
prof_force_df_direct = UNSET_I; prof_force_df_store = UNSET_I; prof_force_df_force_direct = UNSET_I
prof_vv10_radial = UNSET_I; prof_vv10_angular = UNSET_I
prof_xc_dyngrid_prms = UNSET_R; prof_xc_dyngrid_rscale = UNSET_R; prof_xc_pscreen = UNSET_R
prof_vv10_rcut = UNSET_R; prof_vv10_rhocut = UNSET_R
prof_xc_grid_level = UNSET_I
prof_xc_grid_level_coarse = UNSET_I; prof_xc_grid_level_fine = UNSET_I
prof_cosx_hi_level = UNSET_I; prof_cosx_md_level = UNSET_I; prof_cosx_lo_level = UNSET_I

read(buf, nml=molecule)

call find_namelist_block(lines, nlines, 'atoms', buf, start, finish, found)
if (found) then
   if (is_old_namelist_atoms(lines, start, finish)) then
      call parse_atoms_namelist_mode(buf, ncenters, atomchg, x, y, z)
   else
      call parse_atoms_column(lines, start, finish, &
              icharge, imult, unit, baselabel, ecplabel, functional, &
              atomchg, x, y, z, atom_basis, atom_ecp, na)
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
   read(buf, nml=opt)
   runtype = 'opt'
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
if (found) read(buf, nml=cosmo)
if (trim(cosmo_solvent) /= '') cosmo_epsilon = cosmo_solvent_epsilon(trim(cosmo_solvent))

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
if (prof_xc_dyngrid /= UNSET_I)         call setenv_int("ENGINE_XC_DYNGRID", prof_xc_dyngrid)
if (prof_xc_dyngrid_prms /= UNSET_R)    call setenv_real("ENGINE_XC_DYNGRID_PRMS", prof_xc_dyngrid_prms)
if (prof_xc_dyngrid_rscale /= UNSET_R)  call setenv_real("ENGINE_XC_DYNGRID_RSCALE", prof_xc_dyngrid_rscale)
if (prof_xc_dyngrid_sph /= UNSET_I)     call setenv_int("ENGINE_XC_DYNGRID_SPH", prof_xc_dyngrid_sph)
if (prof_xc_shellpair /= UNSET_I)       call setenv_int("ENGINE_XC_SHELLPAIR", prof_xc_shellpair)
if (prof_xc_pscreen /= UNSET_R)         call setenv_real("ENGINE_XC_PSCREEN", prof_xc_pscreen)
if (prof_xc_incr /= UNSET_I)            call setenv_int("ENGINE_XC_INCR", prof_xc_incr)
if (prof_incremental_fock /= UNSET_I)   call setenv_int("ENGINE_INCREMENTAL_FOCK", prof_incremental_fock)
if (prof_diis_window /= UNSET_I)        call setenv_int("ENGINE_DIIS_WINDOW", prof_diis_window)
if (prof_force_df_direct /= UNSET_I)    call setenv_int("ENGINE_FORCE_DF_DIRECT", prof_force_df_direct)
if (prof_force_df_store /= UNSET_I)     call setenv_int("ENGINE_FORCE_DF_STORE", prof_force_df_store)
if (prof_force_df_force_direct /= UNSET_I) call setenv_int("ENGINE_FORCE_DF_FORCE_DIRECT", prof_force_df_force_direct)
if (prof_vv10_radial /= UNSET_I)        call setenv_int("ENGINE_VV10_RADIAL", prof_vv10_radial)
if (prof_vv10_angular /= UNSET_I)       call setenv_int("ENGINE_VV10_ANGULAR", prof_vv10_angular)
if (prof_vv10_rcut /= UNSET_R)          call setenv_real("ENGINE_VV10_RCUT", prof_vv10_rcut)
if (prof_vv10_rhocut /= UNSET_R)        call setenv_real("ENGINE_VV10_RHOCUT", prof_vv10_rhocut)
if (prof_xc_grid_level /= UNSET_I)        call setenv_int("ENGINE_XC_GRID_LEVEL", prof_xc_grid_level)
if (prof_xc_grid_level_coarse /= UNSET_I) call setenv_int("ENGINE_XC_GRID_LEVEL_COARSE", prof_xc_grid_level_coarse)
if (prof_xc_grid_level_fine /= UNSET_I)   call setenv_int("ENGINE_XC_GRID_LEVEL_FINE", prof_xc_grid_level_fine)
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
   call read_xyz(xyzfile, atomchg, x, y, z, na)
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
spec%ecplabel = ecplabel
spec%basedir = basedir
spec%J = J
spec%K = K
spec%ri_aux_basis = ri_aux_basis
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
   character(16) :: rt, oc
   rt = adjustl(runtype)
   select case (trim(rt))
   case ('energy', 'Energy', 'ENERGY', 'sp', 'SP', '')
      spec%opt_run = .false.
   case ('opt', 'Opt', 'OPT', 'optimize', 'optimization')
      spec%opt_run = .true.
   case default
      write(*,'(A)') "Input error: runtype = '"//trim(rt)// &
           "' not recognized (use energy / opt)"
      stop 1
   end select
   oc = adjustl(opt_conv)
   select case (trim(oc))
   case ('normal', 'Normal', 'NORMAL', 'regular', '')
      spec%opt_conv_level = 0
   case ('tight', 'Tight', 'TIGHT')
      spec%opt_conv_level = 1
   case default
      write(*,'(A)') "Input error: opt_conv = '"//trim(oc)// &
           "' not recognized (use normal / tight)"
      stop 1
   end select
end block
spec%opt_maxcyc  = opt_maxcyc
spec%opt_trust   = opt_trust
spec%opt_coord   = opt_coord
spec%opt_restart = opt_restart

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
spec%cosmo_on = cosmo_on; spec%cosmo_epsilon = cosmo_epsilon
spec%cosmo_radii_scale = cosmo_radii_scale; spec%cosmo_avg_area = cosmo_avg_area
spec%cosmo_sigma_rav = cosmo_sigma_rav; spec%cosmo_cavity_type = cosmo_cavity_type
spec%cosmo_rsolv = cosmo_rsolv; spec%cosmo_ks_nseg = cosmo_ks_nseg
spec%cosmo_ks_nface = cosmo_ks_nface; spec%cosmo_sigma_profile_file = cosmo_sigma_profile_file
spec%cosmo_solvent = cosmo_solvent; spec%cosmo_smd = cosmo_smd

call print_input_summary(infile, spec)

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

subroutine print_input_summary(infile, spec)
use mod_engine_input_types, only: engine_input_t
use mod_engine_input_elements, only: elem_to_name
character(len=*), intent(in)     :: infile
type(engine_input_t), intent(in) :: spec
integer :: i

print *, '[Input] ', trim(infile)
print *, '[HEADEREND]'
print *
print *, '[SYSTEM]'
print '("  Atoms: ", I0, "   (coordinates in Angstrom)")', spec%ncenters
do i = 1, spec%ncenters
   print '("    ", A2, 3F14.6, "   basis=", A, "   ecp=", A)', elem_to_name(spec%atomchg(i)), &
      spec%x(i), spec%y(i), spec%z(i), trim(atom_label(spec%atom_basis(i), spec%baselabel)), &
      trim(atom_ecp_label(spec%atom_ecp(i), spec%ecplabel))
enddo
print *, '  Functional: ', trim(spec%functional), '  Basis: ', trim(spec%baselabel)
if (spec%ecplabel /= '') print *, '  ECP: ', trim(spec%ecplabel)
print '(A)', '  Method: '//trim(jk_method_label(spec%J, spec%K))// &
             ' (J='//trim(spec%J)//', K='//trim(spec%K)//')'
print *, '  Charge: ', spec%icharge, '  Multiplicity: ', spec%imult

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
