! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! The single-point results block written to a Direwolf .out file:

module mod_engine_results
implicit none
private
public :: write_results_block, finish_engine_call, report_progress

contains

subroutine finish_engine_call(opts, spec, outfile, out_unit)
use mod_engine_options, only: engine_options_t
use mod_engine_input_types, only: engine_input_t
use mod_citations, only: write_citations, citations_bib_name
use mod_json_out, only: json_write_report
use mod_output_format, only: finalize_output_file
type(engine_options_t), intent(in) :: opts
type(engine_input_t), intent(in) :: spec
character(len=*), intent(in) :: outfile
integer, intent(in) :: out_unit
call json_write_report(opts%rep)
if (opts%done) then
   call write_citations(spec, out_unit, citations_bib_name(outfile))
   call finalize_output_file()
   stop
endif
end subroutine finish_engine_call

subroutine report_progress(rep)
use mod_engine_report, only: engine_report_t
use mod_json_out, only: json_write_report
type(engine_report_t), intent(in) :: rep
call json_write_report(rep)
end subroutine report_progress

subroutine write_results_block(out_unit, ncenters, energy_out, force_out, iconv, want_forces, opts)
use mod_engine_options, only: engine_options_t
use mod_json_out, only: json_result
type(engine_options_t), intent(in) :: opts
integer, intent(in) :: out_unit, ncenters, iconv
real(8), intent(in) :: energy_out
real(8), intent(in) :: force_out(ncenters, 3)
logical, intent(in) :: want_forces

integer :: i

write(out_unit, '(A, A)')     'converged            = ', merge('Yes', 'No ', iconv == 1)
write(out_unit, '(A, F20.9)') 'Total Energy (Hartree)  = ', energy_out
block
   character(len=8) :: hiprec_env
   call get_environment_variable('ENGINE_ENERGY_HIPREC', hiprec_env)
   if (trim(hiprec_env) .eq. '1') &
      write(out_unit, '(A, ES28.18)') 'Total Energy HIPREC     = ', energy_out
end block
write(out_unit, '(A, F20.9)') 'Hcore trace (Hartree)   = ', opts%e_hcore
write(out_unit, '(A, F20.9)') 'Ecoul (Hartree)         = ', opts%e_coul
write(out_unit, '(A, F20.9)') 'Exact exchange (Hartree)= ', opts%e_exx
write(out_unit, '(A, F20.9)') 'Exc (Hartree)           = ', opts%e_xc
write(out_unit, '(A, F20.9)') 'Nuclear repulsion (Ha)  = ', opts%e_nuc

if (opts%vv10_used) write(out_unit, '(A, F20.9)') 'VV10 NLC Energy (Hartree) = ', opts%e_vv10

call json_result(iconv == 1, energy_out, opts%e_hcore, opts%e_coul, &
     opts%e_exx, opts%e_xc, opts%e_nuc, &
     opts%vv10_used, opts%e_vv10, want_forces, force_out)

if (want_forces) then
   write(out_unit, '(A)') 'forces(hartree/bohr):'
   do i = 1, ncenters
      write(out_unit, '(3F18.9)') force_out(i, :)
   enddo
endif
end subroutine write_results_block

end module mod_engine_results
