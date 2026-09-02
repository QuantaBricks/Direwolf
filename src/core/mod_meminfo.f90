! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_meminfo: reports currently-available system memory, used by

module mod_meminfo
implicit none
private
public :: get_available_memory_bytes
public :: get_memory_budget_bytes
public :: get_store_threshold_bytes
public :: engine_mem_cap_bytes
public :: get_process_memory_bytes
public :: engine_avail_at_start_bytes

integer(8),parameter :: FALLBACK_BYTES = 8_8*1024_8*1024_8*1024_8

integer(8) :: engine_mem_cap_bytes = 0_8

integer(8) :: engine_avail_at_start_bytes = 0_8

contains

integer(8) function get_available_memory_bytes() result(bytes)
    implicit none
    integer :: unit,ios
    character(256) :: line
    integer(8) :: kb
    logical :: found

    bytes = FALLBACK_BYTES
    found = .false.

    open(newunit=unit,file="/proc/meminfo",status="old",action="read",iostat=ios)
    if (ios /= 0) return

    do
       read(unit,"(A)",iostat=ios) line
       if (ios /= 0) exit
       if (index(line,"MemAvailable:") == 1) then
          read(line(14:),*,iostat=ios) kb
          if (ios == 0) then
             bytes = kb * 1024_8
             found = .true.
          endif
          exit
       endif
    enddo
    close(unit)

    if (.not. found) bytes = FALLBACK_BYTES

end function get_available_memory_bytes

integer(8) function get_memory_budget_bytes() result(bytes)
    implicit none
    integer(8) :: available_now
    available_now = get_available_memory_bytes()
    if (engine_mem_cap_bytes .gt. 0_8) then
       bytes = min(engine_mem_cap_bytes, available_now)
    else
       bytes = available_now
    endif
end function get_memory_budget_bytes

subroutine get_process_memory_bytes(rss_bytes, hwm_bytes)
    implicit none
    integer(8),intent(out) :: rss_bytes, hwm_bytes
    integer :: unit, ios
    character(256) :: line
    integer(8) :: val_kb
    rss_bytes = 0_8
    hwm_bytes = 0_8
    open(newunit=unit,file="/proc/self/status",status="old",action="read",iostat=ios)
    if (ios /= 0) return
    do
       read(unit,"(A)",iostat=ios) line
       if (ios /= 0) exit
       if (index(line,"VmRSS:") == 1) then
          read(line(7:),*,iostat=ios) val_kb
          if (ios == 0) rss_bytes = val_kb * 1024_8
       else if (index(line,"VmHWM:") == 1) then
          read(line(7:),*,iostat=ios) val_kb
          if (ios == 0) hwm_bytes = val_kb * 1024_8
       endif
    enddo
    close(unit)
end subroutine get_process_memory_bytes

integer(8) function get_store_threshold_bytes() result(bytes)
    implicit none
    integer(8) :: budget, reserve
    real(8),parameter :: safety_fraction = 0.5d0, reserve_fraction = 0.1d0
    integer(8),parameter :: reserve_cap_bytes = 4_8*1024_8*1024_8*1024_8
    budget = get_memory_budget_bytes()
    reserve = min(reserve_cap_bytes, int(budget*reserve_fraction, 8))
    bytes = max(0_8, int(budget*safety_fraction, 8) - reserve)
end function get_store_threshold_bytes

end module mod_meminfo
