! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_profile: minimal wall-clock section timer, used to profile

module mod_profile
use mod_meminfo, only: get_process_memory_bytes
implicit none
private
public :: prof_reset, prof_start, prof_stop, prof_report, fmt_gb, itoa

interface itoa
   module procedure itoa_i4, itoa_i8
end interface itoa

integer,parameter :: MAX_SECTIONS = 32
integer,parameter :: NAME_LEN = 32

character(NAME_LEN) :: sec_name(MAX_SECTIONS)
integer(8)          :: sec_start_count(MAX_SECTIONS)
integer(8)          :: sec_total_count(MAX_SECTIONS)
integer             :: sec_call_count(MAX_SECTIONS)
real(8)             :: sec_rss_gb(MAX_SECTIONS)
real(8)             :: sec_hwm_gb(MAX_SECTIONS)
integer             :: n_sections = 0
integer(8)          :: count_rate = 0

contains

integer function find_or_add(name) result(idx)
    implicit none
    character(*),intent(in) :: name
    integer :: i
    do i = 1,n_sections
       if (trim(sec_name(i)) == trim(name)) then
          idx = i
          return
       endif
    enddo
    if (n_sections >= MAX_SECTIONS) then
       print *,"mod_profile: MAX_SECTIONS exceeded, dropping section ",trim(name)
       idx = 1
       return
    endif
    n_sections = n_sections + 1
    idx = n_sections
    sec_name(idx) = name
    sec_total_count(idx) = 0
    sec_call_count(idx) = 0
end function find_or_add

subroutine prof_reset()
    implicit none
    n_sections = 0
    sec_total_count = 0
    sec_call_count = 0
    call system_clock(count_rate=count_rate)
end subroutine prof_reset

subroutine prof_start(name)
    implicit none
    character(*),intent(in) :: name
    integer :: idx
    idx = find_or_add(name)
    call system_clock(sec_start_count(idx))
end subroutine prof_start

subroutine prof_stop(name)
    implicit none
    character(*),intent(in) :: name
    integer :: idx
    integer(8) :: now, rss_bytes, hwm_bytes
    idx = find_or_add(name)
    call system_clock(now)
    sec_total_count(idx) = sec_total_count(idx) + (now - sec_start_count(idx))
    sec_call_count(idx) = sec_call_count(idx) + 1
    call get_process_memory_bytes(rss_bytes, hwm_bytes)
    sec_rss_gb(idx) = real(rss_bytes,8)/1024.0d0**3
    sec_hwm_gb(idx) = real(hwm_bytes,8)/1024.0d0**3
end subroutine prof_stop

subroutine prof_report(unit)
    use iso_fortran_env, only: output_unit
    implicit none
    integer,intent(in),optional :: unit
    integer :: i, out_unit
    real(8) :: secs
    if (count_rate == 0) return
    out_unit = output_unit
    if (present(unit)) out_unit = unit
    write(out_unit,'(A)') '[PROFILE]'
    write(out_unit,*) "===== PROFILE report ====="
    do i = 1,n_sections
       secs = real(sec_total_count(i),8) / real(count_rate,8)
       write(out_unit,'("PROFILE  ",A20,2X,F12.4," s",2X,I8," calls",2X, &
             &"rss=",F8.3,"GB",2X,"peak_rss=",F8.3,"GB")') &
             trim(sec_name(i)), secs, sec_call_count(i), sec_rss_gb(i), sec_hwm_gb(i)
    enddo
    write(out_unit,*) "==========================="
    write(out_unit,'(A)') '[PROFILEEND]'
    write(out_unit,'(A)') ''
end subroutine prof_report

function fmt_gb(val) result(str)
real(8), intent(in) :: val
character(len=16) :: str
write(str, '(F0.3)') val
str = adjustl(str)
if (str(1:1) == '.') then
   str = '0'//str
else if (len_trim(str) >= 2) then
   if (str(1:2) == '-.') str = '-0'//str(2:)
endif
end function fmt_gb

function itoa_i4(n) result(str)
integer, intent(in) :: n
character(len=24) :: str
write(str, '(I0)') n
end function itoa_i4

function itoa_i8(n) result(str)
integer(8), intent(in) :: n
character(len=24) :: str
write(str, '(I0)') n
end function itoa_i8

end module mod_profile
