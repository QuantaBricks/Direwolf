! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! run_engine: a thin, file-driven example front end for the Engine

program run_engine
use mod_engine_input, only: engine_read_input
use mod_engine_input_types, only: engine_input_t
use MOL_info, only: engine_verbose
use mod_integrals, only: df_direct_mode, df_force_direct_mode, &
                          df_energy_need_gb, df_energy_avail_gb, &
                          df_force_need_gb, df_force_avail_gb
use mod_profile, only: prof_report
use mod_version, only: engine_version, engine_git_version
use mod_vv10, only: vv10_report
use mod_checkpoint, only: checkpoint_read, checkpoint_write
use mod_engineup_interface, only: EngineUp
use mod_engine_results, only: write_results_block
use mod_geomopt, only: geomopt_run
implicit none

type(engine_input_t) :: spec

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
integer        :: out_unit

nargs = command_argument_count()
if (nargs >= 1) then
   call get_command_argument(1, infile)
else
   infile = 'input.inp'
endif

if (nargs >= 2) then
   call get_command_argument(2, outfile)
else
   outfile = make_outfile_name(infile)
endif

block
   integer :: verbose_peek
   call peek_verbose(infile, verbose_peek)
   engine_verbose = verbose_peek
   if (verbose_peek <= 1) then
      close(6)
      if (verbose_peek == 1) then
         open(unit=6, file=trim(outfile), status='replace', action='write')
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

allocate(atomchg_a(spec%ncenters), coord_a(spec%ncenters, 3))
allocate(force_out(spec%ncenters, 3), MLcharge_out(spec%ncenters))
do i = 1, spec%ncenters
   atomchg_a(i) = spec%atomchg(i)
   coord_a(i, 1) = spec%x(i)
   coord_a(i, 2) = spec%y(i)
   coord_a(i, 3) = spec%z(i)
enddo

allocate(pc_charge_a(max(spec%npc, 1)), pc_coord_a(max(spec%npc, 1), 3))
do i = 1, spec%npc
   pc_charge_a(i) = spec%pc_q(i)
   pc_coord_a(i, 1) = spec%pc_x(i)
   pc_coord_a(i, 2) = spec%pc_y(i)
   pc_coord_a(i, 3) = spec%pc_z(i)
enddo

block
use MOL_info, only: resp_charges_on
resp_charges_on = spec%resp_charges_on
end block

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

if (spec%opt_run) then
   if (out_unit == 11) then
      open(unit=11, file=trim(outfile), status='replace', action='write')
      write(11, '(A,A)') 'Direwolf version     = ', engine_version
      write(11, '(A,A)') 'Build                = ', engine_git_version
   endif
   call geomopt_run(spec, outfile, out_unit)
   if (out_unit == 11) close(11)
   if (out_unit == 6) then
      close(6)
      call normalize_output_indentation(outfile)
   endif
   print *, 'Done. Output written to ', trim(outfile)
   stop
endif

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
              dens_in_a, dens_in_b, dens_out_a, dens_out_b)

if (spec%chk_write .and. iconv .eq. 1) call checkpoint_write(trim(spec%chk_file), dens_out_a, dens_out_b)

if (out_unit == 11) then
   open(unit=11, file=trim(outfile), status='replace', action='write')
   write(11, '(A,A)')        'Direwolf version     = ', engine_version
   write(11, '(A,A)')        'Build                = ', engine_git_version
endif
if (.not. spec%estimate_only) then
   call write_results_block(out_unit, spec%ncenters, energy_out, force_out, iconv, spec%calc_force)
   if (out_unit == 6) then
      write(out_unit, '(A)') '[RESULTSEND]'
      write(out_unit, '(A)') ''
      call prof_report(out_unit)
   endif
endif
if (out_unit == 11) close(11)

if (out_unit == 6) then
   close(6)
   call normalize_output_indentation(outfile)
endif

print *, 'Done. Output written to ', trim(outfile)

contains

subroutine normalize_output_indentation(filename)
character(len=*), intent(in) :: filename
character(len=512), allocatable :: lines(:)
character(len=512) :: raw, trimmed
integer :: nlines, i, unit_no, ios

open(newunit=unit_no, file=trim(filename), status='old', action='read', iostat=ios)
if (ios /= 0) return
nlines = 0
do
   read(unit_no, '(A)', iostat=ios) raw
   if (ios /= 0) exit
   nlines = nlines + 1
enddo
close(unit_no)
if (nlines == 0) return

allocate(lines(nlines))
open(newunit=unit_no, file=trim(filename), status='old', action='read', iostat=ios)
if (ios /= 0) return
do i = 1, nlines
   read(unit_no, '(A)', iostat=ios) raw
   if (ios /= 0) exit
   lines(i) = raw
enddo
close(unit_no)

open(newunit=unit_no, file=trim(filename), status='replace', action='write')
do i = 1, nlines
   trimmed = adjustl(lines(i))
   if (len_trim(trimmed) == 0) then
      write(unit_no, '(A)') ''
   else if (trimmed(1:1) == '[') then
      write(unit_no, '(A)') trim(trimmed)
   else if ((trimmed(1:1) >= '0' .and. trimmed(1:1) <= '9') .or. trimmed(1:1) == '-') then
      write(unit_no, '(A)') '  '//trim(lines(i))
   else
      write(unit_no, '(A)') '  '//trim(trimmed)
   endif
enddo
close(unit_no)
end subroutine normalize_output_indentation

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
