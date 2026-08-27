! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Generic "&name ... &end" block finder, shared by every namelist group

module mod_engine_input_block
implicit none
private
public :: find_namelist_block, read_input_lines

contains

subroutine read_input_lines(filename, lines, nlines, maxlines)
character(len=*), intent(in)               :: filename
character(len=4096), allocatable, intent(out) :: lines(:)
integer, intent(out)                       :: nlines
integer, intent(in)                        :: maxlines

integer :: unit_no, ios

allocate(lines(maxlines))
open(newunit=unit_no, file=trim(filename), status='old', action='read', iostat=ios)
if (ios /= 0) then
   print *, 'Cannot open input file: ', trim(filename)
   stop 1
endif
nlines = 0
do
   if (nlines >= maxlines) exit
   read(unit_no, '(A)', iostat=ios) lines(nlines + 1)
   if (ios /= 0) exit
   nlines = nlines + 1
enddo
close(unit_no)
end subroutine

subroutine find_namelist_block(lines, nlines, nname, buf, start, finish, found)
character(len=*), intent(in)  :: lines(:)
integer, intent(in)           :: nlines
character(len=*), intent(in)  :: nname
character(len=*), intent(out) :: buf
integer, intent(out)          :: start, finish
logical, intent(out)          :: found

integer :: i, lpad
character(len=4096) :: line

found = .false.
start = 0
finish = 0
buf = ''

do i = 1, nlines
   line = trim(adjustl(lines(i)))
   if (line == '&' // trim(nname)) then
      start = i
      exit
   endif
enddo
if (start == 0) return

do i = start + 1, nlines
   if (trim(adjustl(lines(i))) == '&end') then
      finish = i
      exit
   endif
enddo
if (finish == 0) then
   print *, 'Error: &end not found for &', trim(nname)
   stop 1
endif

do i = start, finish
   line = lines(i)
   lpad = len_trim(line)
   if (lpad > 0) then
      buf = trim(buf) // new_line('a') // trim(line)
   endif
enddo

found = .true.
end subroutine

end module mod_engine_input_block
