! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! The single-point results block written to a Direwolf .out file:

module mod_engine_results
implicit none
private
public :: write_results_block

contains

subroutine write_results_block(out_unit, ncenters, energy_out, force_out, iconv, want_forces)
use mod_scf_history, only: scf_hist_last_exc, scf_hist_last_ecoul, &
     scf_hist_last_exact_exchange, scf_hist_last_hcore, scf_hist_last_erep
use mod_vv10, only: vv10_report
integer, intent(in) :: out_unit, ncenters, iconv
real(8), intent(in) :: energy_out
real(8), intent(in) :: force_out(ncenters, 3)
logical, intent(in) :: want_forces

integer :: i
logical :: vv10_used
real(8) :: vv10_enl

write(out_unit, '(A, A)')     'converged            = ', merge('Yes', 'No ', iconv == 1)
write(out_unit, '(A, F20.9)') 'Total Energy (Hartree)  = ', energy_out
write(out_unit, '(A, F20.9)') 'Hcore trace (Hartree)   = ', scf_hist_last_hcore()
write(out_unit, '(A, F20.9)') 'Ecoul (Hartree)         = ', scf_hist_last_ecoul()
write(out_unit, '(A, F20.9)') 'Exact exchange (Hartree)= ', scf_hist_last_exact_exchange()
write(out_unit, '(A, F20.9)') 'Exc (Hartree)           = ', scf_hist_last_exc()
write(out_unit, '(A, F20.9)') 'Nuclear repulsion (Ha)  = ', scf_hist_last_erep()

call vv10_report(vv10_enl, vv10_used)
if (vv10_used) write(out_unit, '(A, F20.9)') 'VV10 NLC Energy (Hartree) = ', vv10_enl

if (want_forces) then
   write(out_unit, '(A)') 'forces(hartree/bohr):'
   do i = 1, ncenters
      write(out_unit, '(3F18.9)') force_out(i, :)
   enddo
endif
end subroutine write_results_block

end module mod_engine_results
