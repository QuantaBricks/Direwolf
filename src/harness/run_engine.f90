! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! run_engine: a thin, file-driven example front end for the Engine

program run_engine
use mod_engine_input, only: engine_read_input
use mod_engine_input_types, only: engine_input_t, spec_freq_pressure_pa, spec_tddft_cube
use mod_engine_options, only: engine_options_t
use mod_profile, only: prof_report
use mod_output_format, only: set_output_path, finalize_output_file
use mod_version, only: engine_version, engine_git_version
use mod_checkpoint, only: checkpoint_read, checkpoint_write
use mod_engineup_interface, only: EngineUp
use mod_engine_results, only: write_results_block, finish_engine_call, report_progress
use mod_opt_scan, only: scan_or_opt_run
use mod_json_out, only: json_open, json_set_meta, json_set_geometry, json_flush
use mod_engine_input_elements, only: elem_to_name
use mod_engine_input_types, only: runtype_string, parse_polar_omega
use mod_irc, only: irc_run
use mod_citations, only: write_citations
implicit none

type(engine_input_t) :: spec
type(engine_options_t) :: opts

integer,allocatable  :: atomchg_a(:)
real(8),allocatable  :: coord_a(:,:)
real(8),allocatable  :: force_out(:,:), MLcharge_out(:)
real(8),allocatable  :: pc_charge_a(:), pc_coord_a(:,:)
real(8),allocatable  :: dens_in_a(:,:), dens_in_b(:,:)
real(8),allocatable  :: dens_out_a(:,:), dens_out_b(:,:)

integer       :: iconv
real(8)       :: energy_out, econv
real(8)       :: mem_grid_gb, mem_2e_gb

character(256) :: infile, outfile
integer        :: nargs, i
character(256) :: argv(2)
integer        :: nargs_eff
integer(8)     :: out_size0
logical        :: render_only = .false.
logical        :: spec_td_state = .false.
integer        :: out_unit

nargs = command_argument_count()
argv = ''
nargs_eff = 0
do i = 1, nargs
   block
      character(256) :: a
      call get_command_argument(i, a)
      nargs_eff = nargs_eff + 1
      if (nargs_eff <= 2) argv(nargs_eff) = a
   end block
enddo
if (nargs_eff >= 1) then
   infile = argv(1)
else
   infile = 'input.inp'
endif

if (nargs_eff >= 2) then
   outfile = argv(2)
else
   outfile = make_outfile_name(infile)
endif

block
   integer :: verbose_peek
   call peek_verbose(infile, verbose_peek)
   opts%verbose = verbose_peek
   if (render_only) then
      verbose_peek = 2
      opts%verbose = 2
      out_unit = 11
   endif
   if (verbose_peek <= 1) then
      close(6)
      if (verbose_peek == 1) then
         open(unit=6, file=trim(outfile), status='replace', action='write')
         call set_output_path(trim(outfile))
         write(6,'(A)') '[HEADER]'
         write(6,'(A,A)') 'Direwolf version: ', engine_version
         write(6,'(A,A)') 'Build: ', engine_git_version
         out_unit = 6
      else
         open(unit=6, file='/dev/null', status='old', action='write')
         out_unit = 11
      endif
   else
      out_unit = 11
   endif
end block

call engine_read_input(infile, spec)

block
   integer :: lo
   lo = len_trim(outfile)
   opts%freq_base = outfile
   if (lo > 4) then
      if (outfile(lo-3:lo) == '.out') opts%freq_base = outfile(1:lo-4)
   endif
end block

allocate(atomchg_a(spec%ncenters), coord_a(spec%ncenters, 3))
allocate(force_out(spec%ncenters, 3), MLcharge_out(spec%ncenters))
do i = 1, spec%ncenters
   atomchg_a(i) = spec%atomchg(i)
   coord_a(i, 1) = spec%x(i)
   coord_a(i, 2) = spec%y(i)
   coord_a(i, 3) = spec%z(i)
enddo

block
   character(len=256) :: jname
   character(len=3), allocatable :: jsym(:)
   real(8), allocatable :: jxyz(:,:)
   integer :: lj, ja
   lj = len_trim(outfile)
   jname = trim(outfile) // '.json'
   if (lj > 4) then
      if (outfile(lj-3:lj) == '.out') jname = outfile(1:lj-4) // '.json'
   endif
   call json_open(trim(jname))
   allocate(jsym(spec%ncenters), jxyz(3, spec%ncenters))
   do ja = 1, spec%ncenters
      jsym(ja) = elem_to_name(spec%atomchg(ja))
      jxyz(:, ja) = [spec%x(ja), spec%y(ja), spec%z(ja)]
   enddo
   call json_set_meta(trim(engine_version), trim(engine_git_version), trim(spec%functional), trim(spec%baselabel), &
                      spec%icharge, spec%imult, &
                      trim(runtype_string(spec%opt_run, spec%freq, spec%polar, spec%irc_run, spec%nmr_run)), &
                      jsym, spec%atomchg(1:spec%ncenters), jxyz)
   call json_set_geometry(jxyz)
   call json_flush()
end block

allocate(pc_charge_a(max(spec%npc, 1)), pc_coord_a(max(spec%npc, 1), 3))
do i = 1, spec%npc
   pc_charge_a(i) = spec%pc_q(i)
   pc_coord_a(i, 1) = spec%pc_x(i)
   pc_coord_a(i, 2) = spec%pc_y(i)
   pc_coord_a(i, 3) = spec%pc_z(i)
enddo

opts%resp_charges = spec%resp_charges_on
if (.not. (spec%opt_run .or. spec%irc_run)) then
   opts%cube_den = index(spec%cube, 'den') > 0 .or. index(spec%cube, 'DEN') > 0
   opts%cube_esp = index(spec%cube, 'esp') > 0 .or. index(spec%cube, 'ESP') > 0
   opts%cube_esp_exact = index(spec%cube, 'esp_exact') > 0 .or. index(spec%cube, 'ESP_EXACT') > 0
endif
opts%progress => report_progress

allocate(dens_in_a(0,0), dens_in_b(0,0))
if (spec%chk_read) then
   block
      real(8),allocatable :: tmp_a(:,:), tmp_b(:,:)
      logical :: chk_ok
      call checkpoint_read(trim(spec%chk_file), tmp_a, tmp_b, chk_ok)
      if (chk_ok) then
         deallocate(dens_in_a, dens_in_b)
         call move_alloc(tmp_a, dens_in_a)
         call move_alloc(tmp_b, dens_in_b)
      endif
   end block
endif

if (spec%stability_mode >= 0) then
   block
      use iso_c_binding, only: c_char, c_int, c_null_char
      interface
         function c_setenv_st(name, value, overwrite) bind(C, name="setenv") result(r)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
         end function
      end interface
      character(len=16) :: st_str
      integer(c_int) :: st_rc
      write(st_str,'(I0)') spec%stability_mode
      st_rc = c_setenv_st("ENGINE_STABILITY"//c_null_char, trim(st_str)//c_null_char, 1_c_int)
   end block
endif

if (spec%irc_run) then
   if (out_unit == 11) then
      open(unit=11, file=trim(outfile), status='replace', action='write')
      write(11, '(A,A)') 'Direwolf version     = ', engine_version
      write(11, '(A,A)') 'Build                = ', engine_git_version
   endif
   call irc_run(spec, outfile, out_unit)
   call write_citations(spec, out_unit, bib_name(outfile))
   if (out_unit == 11) close(11)
   if (out_unit == 6) then
      call finalize_output_file()
   endif
   print *, 'Done. Output written to ', trim(outfile)
   stop
endif

if (spec%opt_run) then
   if (out_unit == 11) then
      open(unit=11, file=trim(outfile), status='replace', action='write')
      write(11, '(A,A)') 'Direwolf version     = ', engine_version
      write(11, '(A,A)') 'Build                = ', engine_git_version
   endif
   call scan_or_opt_run(spec, outfile, out_unit)
   call write_citations(spec, out_unit, bib_name(outfile))
   if (out_unit == 11) close(11)
   if (out_unit == 6) then
      call finalize_output_file()
   endif
   print *, 'Done. Output written to ', trim(outfile)
   stop
endif

if (spec%newton_cg_iter >= 0) then
   block
      use iso_c_binding, only: c_char, c_int, c_null_char
      interface
         function c_setenv_nc(name, value, overwrite) bind(C, name="setenv") result(r)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: name(*), value(*)
         integer(c_int), value :: overwrite
         integer(c_int) :: r
         end function
      end interface
      character(len=16) :: nc_str
      integer(c_int) :: nc_rc
      write(nc_str,'(I0)') spec%newton_cg_iter
      nc_rc = c_setenv_nc("ENGINE_NEWTON_ITER"//c_null_char, trim(nc_str)//c_null_char, 1_c_int)
   end block
endif
opts%polar = spec%polar
opts%nmr = spec%nmr_run
if (spec%nmr_run) then
   if (trim(spec%J) /= 'RI' .and. trim(spec%J) /= 'exact') &
      write(*,'(A)') " [Notice] runtype nmr: J = 'exact' (the Coulomb part of the shielding is exact or RI)"
   if (trim(spec%K) /= 'exact' .and. trim(spec%K) /= 'cosx') &
      write(*,'(A)') " [Notice] runtype nmr: K = 'exact' (the exchange of the shielding is exact or COSX)"
   if (trim(spec%J) /= 'RI') spec%J = 'exact'
   if (trim(spec%K) /= 'cosx') spec%K = 'exact'
   spec%calc_force = .false.
   spec%scf_conv_level = max(spec%scf_conv_level, 1)
endif
opts%nmr_j = spec%jcoupling
opts%nmr_j_terms = spec%jcoupling_terms
opts%nmr_j_atoms = spec%jcoupling_atoms
if (spec%jcoupling) then
   block
      character(len=32) :: tl
      tl = spec%jcoupling_terms
      if (trim(tl) /= 'all') then
         do while (len_trim(tl) > 0)
            if (index(tl, ',') > 0) then
               if (.not. any(tl(1:index(tl, ',')-1) == (/'fc ', 'sd ', 'pso', 'dso'/))) &
                  error stop "jcoupling_terms: 'all' or fc, sd, pso, dso separated by commas"
               tl = tl(index(tl, ',')+1:)
            else
               if (.not. any(trim(tl) == (/'fc ', 'sd ', 'pso', 'dso'/))) &
                  error stop "jcoupling_terms: 'all' or fc, sd, pso, dso separated by commas"
               exit
            endif
         enddo
      endif
   end block
endif
opts%nics = spec%nics
if (len_trim(spec%nics) > 0 .and. .not. spec%nmr_run) error stop 'nics needs runtype nmr'
if (spec%jcoupling) then
   if (.not. spec%nmr_run) error stop 'jcoupling needs runtype nmr'
   if (spec%uks .or. spec%imult /= 1) error stop 'jcoupling: closed-shell molecules only'
endif
opts%polar_omega = spec%polar_omega
if (len_trim(spec%polar_omega) > 0) then
   block
      real(8) :: om(16)
      integer :: nom
      logical :: omok
      call parse_polar_omega(spec%polar_omega, om, nom, omok)
      if (.not. omok) error stop 'polar_omega: use numbers in a.u., or with a unit (nm, ev), at most 16'
      if (.not. spec%polar) error stop 'polar_omega is only used with runtype polar'
   end block
endif
if (spec%polar) then
   if (spec%opt_run) error stop 'polar with opt is not implemented: optimize first, then run polar on the result'
   spec%scf_conv_level = max(spec%scf_conv_level, 1)
endif
opts%tddft_nroots = spec%tddft_nroots
opts%tddft_cube = spec_tddft_cube(spec)
if (spec%opt_run .or. spec%irc_run) opts%tddft_nroots = 0
opts%tddft_robust = spec%tddft_robust
opts%tddft_triplets = spec%tddft_triplets
opts%tddft_tda = spec%tddft_tda
opts%tddft_eps_opt = spec%tddft_eps_opt(1)
opts%uks = spec%uks
opts%uks_restart = .false.
if (spec%tddft_nroots < 0) error stop 'tddft_nroots must be 0 (no TDDFT) or a positive number of states'
if (spec%tddft_nroots > 0) then
   spec%scf_conv_level = max(spec%scf_conv_level, 1)
   spec_td_state = .false.
   if (.not. spec_td_state .and. .not. (spec%opt_run .or. spec%irc_run .or. spec%freq)) &
      spec%calc_force = .false.
endif
opts%freq = spec%freq
if (spec%freq) then
   spec%scf_conv_level = max(spec%scf_conv_level, 1)
   if (out_unit == 11) then
      open(unit=11, file=trim(outfile), status='replace', action='write')
      write(11, '(A,A)') 'Direwolf version     = ', engine_version
      write(11, '(A,A)') 'Build                = ', engine_git_version
   endif
   opts%freq_unit = out_unit
   opts%freq_mult = spec%imult
endif
if (spec%freq_cond(1) <= 0.0d0 .or. spec%freq_cond(2) <= 0.0d0 .or. spec%freq_cond(3) < 0.0d0) &
   error stop 'freq_temp, freq_pres must be positive and freq_conc not negative'
opts%freq_temp = spec%freq_cond(1)
opts%freq_pres = spec_freq_pressure_pa(spec)
select case (trim(spec%cosmo_convention))
case ('cosmo', '')
   opts%cosmo_cpcm = .false.
case ('cpcm')
   opts%cosmo_cpcm = .true.
case default
   error stop '&cosmo convention must be cosmo or cpcm'
end select
out_size0 = -1
if (out_unit == 6) then
   flush(6)
   inquire(unit=6, size=out_size0)
endif
do
call EngineUp(spec%ncenters, spec%imult, spec%icharge, spec%functional, &
              coord_a, atomchg_a, spec%baselabel, spec%ecplabel, &
              spec%atom_basis, spec%atom_ecp, &
              trim(spec%J), trim(spec%K), trim(spec%ri_aux_basis), spec%spherical, spec%harris_guess, &
              spec%calc_force, spec%vv10_nonself, spec%mem_cap_gb, &
              spec%estimate_only, spec%n_threads, &
              spec%npc, pc_charge_a, pc_coord_a, &
              spec%molden_write, trim(spec%molden_file), spec%molden_read, trim(spec%molden_read_file), &
              spec%cosmo_on, spec%cosmo_epsilon, spec%cosmo_radii_scale, spec%cosmo_avg_area, spec%cosmo_sigma_rav, &
              trim(spec%cosmo_cavity_type), spec%cosmo_rsolv, spec%cosmo_ks_nseg, spec%cosmo_ks_nface, &
              trim(spec%cosmo_sigma_profile_file), &
              spec%cosmo_smd, trim(spec%cosmo_solvent), &
              force_out, energy_out, MLcharge_out, iconv, econv, &
              mem_grid_gb, mem_2e_gb, spec%scf_conv_level, trim(spec%basedir), &
              dens_in_a, dens_in_b, dens_out_a, dens_out_b, opts)
if (.not. opts%uks_restart) exit
opts%uks_restart = .false.
opts%uks = .true.
call move_alloc(dens_out_a, dens_in_a)
call move_alloc(dens_out_b, dens_in_b)
if (out_unit == 6 .and. out_size0 >= 0) then
   block
      character(len=32) :: szs
      flush(6)
      close(6)
      write(szs,'(I0)') out_size0
      call execute_command_line('truncate -s '//trim(szs)//' "'//trim(outfile)//'"')
      open(unit=6, file=trim(outfile), status='old', position='append', action='write')
   end block
endif
print '(A,ES10.3,A)', 'Restricted solution unstable toward an unrestricted one (lowest triplet-type eigenvalue of A+B =', &
      opts%uks_lambda, '): restarting as an unrestricted (UKS) calculation'
enddo

call finish_engine_call(opts, spec, outfile, out_unit)

if (spec%chk_write .and. iconv .eq. 1) call checkpoint_write(trim(spec%chk_file), dens_out_a, dens_out_b)

if (out_unit == 11) then
   open(unit=11, file=trim(outfile), status='replace', action='write')
   write(11, '(A,A)')        'Direwolf version     = ', engine_version
   write(11, '(A,A)')        'Build                = ', engine_git_version
endif
if (.not. spec%estimate_only) then
   call write_results_block(out_unit, spec%ncenters, energy_out, force_out, iconv, spec%calc_force, opts)
   if (out_unit == 6) then
      write(out_unit, '(A)') '[RESULTSEND]'
      write(out_unit, '(A)') ''
      call prof_report(out_unit)
   else
      call prof_report()
   endif
   call write_citations(spec, out_unit, bib_name(outfile))
endif
if (out_unit == 11) close(11)

if (out_unit == 6) then
   call finalize_output_file()
endif

print *, 'Done. Output written to ', trim(outfile)

contains

function bib_name(outfile) result(f)
character(len=*), intent(in) :: outfile
character(len=512) :: f
integer :: l
l = len_trim(outfile)
if (l > 4 .and. outfile(max(1,l-3):l) == '.out') then
   f = outfile(1:l-4)//'.bib'
else
   f = trim(outfile)//'.bib'
endif
end function bib_name

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
if (l > 5) then
   if (name(l-4:l) == '.toml') then
      outname = name(1:l-5) // '.out'
      return
   endif
endif
outname = trim(name) // '.out'
end function

subroutine peek_verbose(infile, verbose_peek)
use mod_engine_input_block, only: find_namelist_block, read_input_lines
character(len=*), intent(in) :: infile
integer, intent(out) :: verbose_peek
character(4096), allocatable :: lines(:)
character(len=10000) :: buf
integer :: nlines, start, finish, vpos, eqpos, ios
logical :: found
integer, parameter :: MAXLINES = 10000

verbose_peek = 1
call read_input_lines(infile, lines, nlines, MAXLINES)
call find_namelist_block(lines, nlines, 'molecule', buf, start, finish, found)
if (.not. found) return

vpos = index(buf, 'verbose')
if (vpos == 0) return
eqpos = index(buf(vpos:), '=')
if (eqpos == 0) return
read(buf(vpos+eqpos:), *, iostat=ios) verbose_peek
if (ios /= 0) verbose_peek = 1
end subroutine

end program run_engine
