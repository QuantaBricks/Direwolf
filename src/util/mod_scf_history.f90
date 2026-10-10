! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_scf_history: records the per-iteration SCF convergence trace

module mod_scf_history
use mod_meminfo, only: get_process_memory_bytes
implicit none
private
public :: scf_hist_reset, scf_hist_record, scf_hist_report, scf_hist_last_exc
public :: scf_hist_niter, scf_hist_last_prms, scf_hist_last_de
public :: scf_hist_set_ecomponents, scf_hist_last_ecoul, scf_hist_last_exact_exchange
public :: scf_hist_last_hcore, scf_hist_last_erep

integer,parameter :: MAX_ITER = 500

integer  :: hist_iter(MAX_ITER)
real(8)  :: hist_prms(MAX_ITER), hist_dE(MAX_ITER), hist_Etot(MAX_ITER)
real(8)  :: hist_Exc(MAX_ITER)
real(8)  :: hist_diag_t(MAX_ITER), hist_dft_t(MAX_ITER), hist_total_t(MAX_ITER)
real(8)  :: hist_rss_gb(MAX_ITER), hist_peak_rss_gb(MAX_ITER)
integer  :: n_hist = 0

real(8),save :: last_ecoul = 0.0d0
real(8),save :: last_exact_exchange = 0.0d0
real(8),save :: last_hcore = 0.0d0
real(8),save :: last_erep = 0.0d0

contains

integer function scf_hist_niter()
    implicit none
    scf_hist_niter = n_hist
end function scf_hist_niter

real(8) function scf_hist_last_prms()
    implicit none
    if (n_hist .ge. 1) then
       scf_hist_last_prms = hist_prms(n_hist)
    else
       scf_hist_last_prms = 0.0d0
    endif
end function scf_hist_last_prms

real(8) function scf_hist_last_de()
    implicit none
    if (n_hist .ge. 1) then
       scf_hist_last_de = hist_dE(n_hist)
    else
       scf_hist_last_de = 0.0d0
    endif
end function scf_hist_last_de

subroutine scf_hist_reset()
    implicit none
    n_hist = 0
end subroutine scf_hist_reset

real(8) function scf_hist_last_exc()
    implicit none
    if (n_hist .ge. 1) then
       scf_hist_last_exc = hist_Exc(n_hist)
    else
       scf_hist_last_exc = 0.0d0
    endif
end function scf_hist_last_exc

subroutine scf_hist_set_ecomponents(ecoul, exact_exchange, hcore_trace, erep)
    implicit none
    real(8),intent(in) :: ecoul, exact_exchange, hcore_trace, erep
    last_ecoul = ecoul
    last_exact_exchange = exact_exchange
    last_hcore = hcore_trace
    last_erep = erep
end subroutine scf_hist_set_ecomponents

real(8) function scf_hist_last_ecoul()
    implicit none
    scf_hist_last_ecoul = last_ecoul
end function scf_hist_last_ecoul

real(8) function scf_hist_last_exact_exchange()
    implicit none
    scf_hist_last_exact_exchange = last_exact_exchange
end function scf_hist_last_exact_exchange

real(8) function scf_hist_last_hcore()
    implicit none
    scf_hist_last_hcore = last_hcore
end function scf_hist_last_hcore

real(8) function scf_hist_last_erep()
    implicit none
    scf_hist_last_erep = last_erep
end function scf_hist_last_erep

subroutine scf_hist_record(iter, prms, dE, Etot, Exc, diag_t, dft_t, total_t)
    implicit none
    integer,intent(in) :: iter
    real(8),intent(in) :: prms, dE, Etot, Exc, diag_t, dft_t, total_t
    integer(8) :: rss_bytes, hwm_bytes
    if (n_hist .ge. MAX_ITER) return
    n_hist = n_hist + 1
    hist_iter(n_hist) = iter
    hist_prms(n_hist) = prms
    hist_dE(n_hist) = dE
    hist_Etot(n_hist) = Etot
    hist_Exc(n_hist) = Exc
    hist_diag_t(n_hist) = diag_t
    hist_dft_t(n_hist) = dft_t
    hist_total_t(n_hist) = total_t
    call get_process_memory_bytes(rss_bytes, hwm_bytes)
    hist_rss_gb(n_hist) = real(rss_bytes,8)/1024.0d0**3
    hist_peak_rss_gb(n_hist) = real(hwm_bytes,8)/1024.0d0**3
end subroutine scf_hist_record

subroutine scf_hist_report(unit)
    use iso_fortran_env, only: output_unit
    implicit none
    integer,intent(in),optional :: unit
    integer :: i, out_unit
    out_unit = output_unit
    if (present(unit)) out_unit = unit
    write(out_unit,*) "===== SCF iteration history ====="
    write(out_unit,'(A4,2X,A14,2X,A14,2X,A16,2X,A14,5(2X,A10))') &
          "iter","density_change","E_change","total_energy","Exc", &
          "diag_s","dft_s","total_s","rss_gb","peak_gb"
    do i = 1,n_hist
       write(out_unit,'(I4,2X,F14.9,2X,F14.9,2X,F16.9,2X,F14.9,5(2X,F10.3))') &
             hist_iter(i), hist_prms(i), hist_dE(i), hist_Etot(i), hist_Exc(i), &
             hist_diag_t(i), hist_dft_t(i), hist_total_t(i), &
             hist_rss_gb(i), hist_peak_rss_gb(i)
    enddo
    write(out_unit,*) "==================================="
end subroutine scf_hist_report

end module mod_scf_history
