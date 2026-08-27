! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Reads/writes SCF density-matrix checkpoint files for restart.

module mod_checkpoint
implicit none
private
public :: checkpoint_write, checkpoint_try_read

contains

subroutine checkpoint_write(filename)
use MOL_info, only: nConts, Natoms, Pa, Pb
implicit none
character(len=*),intent(in) :: filename
integer :: unit_chk

open(newunit=unit_chk, file=trim(filename), form='unformatted', &
     access='stream', status='replace', action='write')
write(unit_chk) nConts, Natoms
write(unit_chk) Pa
write(unit_chk) Pb
close(unit_chk)

end subroutine checkpoint_write

subroutine checkpoint_try_read(filename, Pa_out, Pb_out, success)
use MOL_info, only: nConts, Natoms
implicit none
character(len=*),intent(in) :: filename
real(8),intent(out) :: Pa_out(nConts,nConts),Pb_out(nConts,nConts)
logical,intent(out) :: success
integer :: unit_chk, nConts_chk, Natoms_chk, ios
logical :: file_exists

success = .false.
inquire(file=trim(filename), exist=file_exists)
if (.not. file_exists) return

open(newunit=unit_chk, file=trim(filename), form='unformatted', &
     access='stream', status='old', action='read', iostat=ios)
if (ios .ne. 0) return

read(unit_chk, iostat=ios) nConts_chk, Natoms_chk
if (ios .ne. 0 .or. nConts_chk .ne. nConts .or. Natoms_chk .ne. Natoms) then
   close(unit_chk)
   return
endif

read(unit_chk, iostat=ios) Pa_out
if (ios .ne. 0) then
   close(unit_chk)
   return
endif
read(unit_chk, iostat=ios) Pb_out
close(unit_chk)
if (ios .ne. 0) return

success = .true.

end subroutine checkpoint_try_read

end module mod_checkpoint
