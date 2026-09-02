! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Reads/writes SCF density-matrix checkpoint files, for warm-starting a

module mod_checkpoint
implicit none
private
public :: checkpoint_write, checkpoint_read

contains

subroutine checkpoint_write(filename, Pa, Pb)
implicit none
character(len=*),intent(in) :: filename
real(8),intent(in) :: Pa(:,:), Pb(:,:)
integer :: unit_chk, nConts

nConts = size(Pa,1)
open(newunit=unit_chk, file=trim(filename), form='unformatted', &
     access='stream', status='replace', action='write')
write(unit_chk) nConts
write(unit_chk) Pa
write(unit_chk) Pb
close(unit_chk)

end subroutine checkpoint_write

subroutine checkpoint_read(filename, Pa_out, Pb_out, success)
implicit none
character(len=*),intent(in) :: filename
real(8),allocatable,intent(out) :: Pa_out(:,:), Pb_out(:,:)
logical,intent(out) :: success
integer :: unit_chk, nConts_chk, ios
logical :: file_exists

success = .false.
inquire(file=trim(filename), exist=file_exists)
if (.not. file_exists) return

open(newunit=unit_chk, file=trim(filename), form='unformatted', &
     access='stream', status='old', action='read', iostat=ios)
if (ios .ne. 0) return

read(unit_chk, iostat=ios) nConts_chk
if (ios .ne. 0 .or. nConts_chk .le. 0) then
   close(unit_chk)
   return
endif

allocate(Pa_out(nConts_chk,nConts_chk), Pb_out(nConts_chk,nConts_chk))
read(unit_chk, iostat=ios) Pa_out
if (ios .eq. 0) read(unit_chk, iostat=ios) Pb_out
close(unit_chk)
if (ios .ne. 0) then
   deallocate(Pa_out, Pb_out)
   return
endif

success = .true.

end subroutine checkpoint_read

end module mod_checkpoint
