! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_output_format: one layout rule for the verbose=1 output file.

module mod_output_format
implicit none
private
public :: out_path, set_output_path, finalize_output_file, normalize_output_file
public :: set_buildmode, print_buildmode

character(len=512), save :: out_path = ''
character(len=64), save :: bm_j = '', bm_k = '', bm_vxc = ''
character(len=*), parameter :: TABLE_BLOCKS = '|NORMALMODES|DIPOLEDERIV|IR|THERMO|POLARIZABILITY|IRCWALK|'

contains

subroutine set_buildmode(j, k, vxc)
character(len=*), intent(in) :: j, k, vxc
bm_j = j; bm_k = k; bm_vxc = vxc
end subroutine set_buildmode

subroutine print_buildmode()
if (len_trim(bm_j) == 0) return
write(*,'(A)') '[BUILDMODE]'
write(*,'(A)') '  J:   '//trim(bm_j)
write(*,'(A)') '  K:   '//trim(bm_k)
write(*,'(A)') '  Vxc: '//trim(bm_vxc)
write(*,'(A)') '[BUILDMODEEND]'
bm_j = ''
end subroutine print_buildmode

subroutine set_output_path(name)
character(len=*), intent(in) :: name
out_path = name
end subroutine set_output_path

subroutine finalize_output_file()
logical :: opened
if (len_trim(out_path) == 0) return
inquire(unit=6, opened=opened)
if (opened) close(6)
call normalize_output_file(trim(out_path))
out_path = ''
open(unit=6, file='/dev/stdout', status='old', action='write')
end subroutine finalize_output_file

logical function is_marker(s)
character(len=*), intent(in) :: s
is_marker = .false.
if (len_trim(s) == 0) return
is_marker = s(1:1) == '['
end function is_marker

logical function is_end_marker(s)
character(len=*), intent(in) :: s
integer :: l
is_end_marker = .false.
if (.not. is_marker(s)) return
l = len_trim(s)
if (l < 6) return
is_end_marker = s(l-3:l) == 'END]' .and. index(s(1:l), ' ') == 0
end function is_end_marker

function marker_name(s) result(name)
character(len=*), intent(in) :: s
character(len=64) :: name
integer :: e
name = ''
e = index(s, ']')
if (e > 2) name = s(2:min(e-1, 65))
end function marker_name

integer function lead_blanks(s)
character(len=*), intent(in) :: s
lead_blanks = verify(s, ' ') - 1
if (lead_blanks < 0) lead_blanks = 0
end function lead_blanks

subroutine normalize_output_file(filename)
character(len=*), intent(in) :: filename
integer, parameter :: LL = 1024
character(len=LL), allocatable :: lines(:)
character(len=LL) :: raw
integer :: nlines, i, j, k, kk, u, ios, m, nblank
logical :: after_end, wrote_any, table_block, loose_text
character(len=64) :: cur_block
character(len=64) :: stack(16), name
integer :: depth

open(newunit=u, file=trim(filename), status='old', action='read', iostat=ios)
if (ios /= 0) return
nlines = 0
do
   read(u, '(A)', iostat=ios) raw
   if (ios /= 0) exit
   nlines = nlines + 1
enddo
if (nlines == 0) then
   close(u); return
endif
allocate(lines(nlines))
rewind(u)
do i = 1, nlines
   read(u, '(A)', iostat=ios) lines(i)
   if (ios /= 0) lines(i) = ''
enddo
close(u)

open(newunit=u, file=trim(filename), status='replace', action='write')
after_end = .false.
wrote_any = .false.
nblank = 0
depth = 0
cur_block = ''
loose_text = .false.
i = 1
do while (i <= nlines)
   if (len_trim(lines(i)) == 0) then
      nblank = nblank + 1
      i = i + 1
      cycle
   endif
   if (nblank > 0 .and. wrote_any .and. .not. after_end) write(u, '(A)') ''
   nblank = 0
   after_end = .false.
   if (is_marker(adjustl(lines(i)))) then
      if (loose_text .and. depth == 0 .and. .not. is_end_marker(adjustl(lines(i)))) then
         write(u, '(A)') ''
         write(u, '(A)') ''
      endif
      loose_text = .false.
      write(u, '(A)') trim(adjustl(lines(i)))
      wrote_any = .true.
      if (is_end_marker(adjustl(lines(i)))) then
         name = marker_name(adjustl(lines(i)))
         name = name(1:max(len_trim(name)-3,0))
         do kk = depth, 1, -1
            if (trim(stack(kk)) == trim(name)) then
               depth = kk - 1
               exit
            endif
         enddo
         cur_block = ''
         if (depth > 0) cur_block = stack(depth)
         k = i + 1
         do while (k <= nlines)
            if (len_trim(lines(k)) > 0) exit
            k = k + 1
         enddo
         if (k <= nlines) then
            write(u, '(A)') ''
            write(u, '(A)') ''
         endif
         after_end = .true.
      else if (depth < 16) then
         depth = depth + 1
         stack(depth) = marker_name(adjustl(lines(i)))
         cur_block = stack(depth)
      endif
      i = i + 1
      cycle
   endif
   j = i
   m = huge(0)
   do while (j <= nlines)
      if (len_trim(lines(j)) == 0) exit
      if (is_marker(adjustl(lines(j)))) exit
      m = min(m, lead_blanks(lines(j)))
      j = j + 1
   enddo
   table_block = len_trim(cur_block) > 0 .and. index(TABLE_BLOCKS, '|'//trim(cur_block)//'|') > 0
   do k = i, j-1
      if (table_block) then
         write(u, '(A)') '  '//trim(lines(k)(m+1:))
      else if (verify(lines(k)(lead_blanks(lines(k))+1:lead_blanks(lines(k))+1), '0123456789-') == 0) then
         write(u, '(A)') '  '//trim(lines(k))
      else
         write(u, '(A)') '  '//trim(adjustl(lines(k)))
      endif
   enddo
   wrote_any = .true.
   loose_text = (depth == 0)
   i = j
enddo
close(u)
deallocate(lines)
end subroutine normalize_output_file

end module mod_output_format
