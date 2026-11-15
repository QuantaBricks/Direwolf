! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Generic "&name ... &end" block finder, shared by every namelist group

module mod_engine_input_block
implicit none
private
public :: find_namelist_block, read_input_lines, strip_key_prefix, rename_keys, quote_value

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

subroutine strip_key_prefix(buf, prefix, names)
character(len=*), intent(inout) :: buf
character(len=*), intent(in)    :: prefix, names(:)
character(len=len(buf)) :: low
integer :: i, n, lp, k, ln, j
character :: q
logical :: hit

lp = len_trim(prefix)
n = len_trim(buf)
low = buf
do i = 1, n
   if (low(i:i) >= 'A' .and. low(i:i) <= 'Z') low(i:i) = char(ichar(low(i:i)) + 32)
enddo
q = ' '
i = 1
do while (i <= n)
   if (q /= ' ') then
      if (buf(i:i) == q) q = ' '
      i = i + 1
      cycle
   endif
   if (buf(i:i) == "'" .or. buf(i:i) == '"') then
      q = buf(i:i)
      i = i + 1
      cycle
   endif
   hit = .false.
   if (i + lp - 1 <= n) then
      if (low(i:i+lp-1) == prefix) then
         if (i == 1) then
            hit = .true.
         else
            hit = .not. (is_word(low(i-1:i-1)))
         endif
      endif
   endif
   if (hit) then
      hit = .false.
      do k = 1, size(names)
         ln = len_trim(names(k))
         if (i + lp + ln - 1 > n) cycle
         if (low(i+lp:i+lp+ln-1) /= trim(names(k))) cycle
         if (i + lp + ln <= n) then
            if (is_word(low(i+lp+ln:i+lp+ln))) cycle
         endif
         hit = .true.
         exit
      enddo
   endif
   if (hit) then
      do j = i, n - lp
         buf(j:j) = buf(j+lp:j+lp)
         low(j:j) = low(j+lp:j+lp)
      enddo
      buf(n-lp+1:n) = ' '
      low(n-lp+1:n) = ' '
      n = n - lp
   else
      i = i + 1
   endif
enddo

contains
   logical function is_word(c)
   character, intent(in) :: c
   is_word = (c >= 'a' .and. c <= 'z') .or. (c >= '0' .and. c <= '9') .or. c == '_'
   end function is_word
end subroutine strip_key_prefix

subroutine rename_keys(buf, alias, target)
character(len=*), intent(inout) :: buf
character(len=*), intent(in)    :: alias(:), target(:)
character(len=len(buf)) :: low, out
integer :: i, n, k, ln, m, j, kh
character :: q
logical :: hit

n = len_trim(buf)
low = buf
do i = 1, n
   if (low(i:i) >= 'A' .and. low(i:i) <= 'Z') low(i:i) = char(ichar(low(i:i)) + 32)
enddo
out = ' '
m = 0
q = ' '
i = 1
do while (i <= n)
   if (q /= ' ') then
      if (buf(i:i) == q) q = ' '
      m = m + 1 ; out(m:m) = buf(i:i)
      i = i + 1
      cycle
   endif
   if (buf(i:i) == '!') then
      do while (i <= n)
         if (buf(i:i) == new_line('a')) exit
         m = m + 1 ; out(m:m) = buf(i:i)
         i = i + 1
      enddo
      cycle
   endif
   if (buf(i:i) == "'" .or. buf(i:i) == '"') q = buf(i:i)
   hit = .false.
   kh = 0
   if (q == ' ') then
      if (i == 1) then
         hit = .true.
      else
         hit = .not. is_word(low(i-1:i-1))
      endif
   endif
   if (hit) then
      hit = .false.
      do k = 1, size(alias)
         ln = len_trim(alias(k))
         if (i + ln - 1 > n) cycle
         if (low(i:i+ln-1) /= trim(alias(k))) cycle
         if (i + ln <= n) then
            if (is_word(low(i+ln:i+ln))) cycle
         endif
         j = i + ln
         do while (j <= n)
            if (low(j:j) /= ' ') exit
            j = j + 1
         enddo
         if (j > n) cycle
         if (low(j:j) /= '=' .and. low(j:j) /= '(') cycle
         hit = .true.
         kh = k
         exit
      enddo
   endif
   if (hit) then
      ln = len_trim(target(kh))
      if (m + ln + 1 > len(out)) exit
      out(m+1:m+ln) = trim(target(kh))
      m = m + ln
      i = i + len_trim(alias(kh))
      cycle
   endif
   m = m + 1 ; out(m:m) = buf(i:i)
   i = i + 1
enddo
if (i > n) buf = out

contains
   logical function is_word(c)
   character, intent(in) :: c
   is_word = (c >= 'a' .and. c <= 'z') .or. (c >= '0' .and. c <= '9') .or. c == '_'
   end function is_word
end subroutine rename_keys

subroutine quote_value(buf, key)
character(len=*), intent(inout) :: buf
character(len=*), intent(in)    :: key
character(len=len(buf)) :: low
integer :: i, n, lk, j, e
n = len_trim(buf)
low = buf
do i = 1, n
   if (low(i:i) >= 'A' .and. low(i:i) <= 'Z') low(i:i) = char(ichar(low(i:i)) + 32)
enddo
lk = len_trim(key)
i = 1
do while (i + lk - 1 <= n)
   if (low(i:i) == '!') then
      do while (i <= n)
         if (buf(i:i) == new_line('a')) exit
         i = i + 1
      enddo
      cycle
   endif
   if (low(i:i+lk-1) == key) then
      if (i == 1 .or. .not. word_char(low(i-1:i-1))) then
         j = i + lk
         do while (j <= n)
            if (buf(j:j) /= ' ') exit
            j = j + 1
         enddo
         if (j < n .and. buf(j:j) == '=') then
            j = j + 1
            do while (j <= n)
               if (buf(j:j) /= ' ') exit
               j = j + 1
            enddo
            if (j <= n .and. buf(j:j) /= "'" .and. buf(j:j) /= '"') then
               e = j
               do while (e <= n)
                  if (buf(e:e) == ' ' .or. buf(e:e) == ',' .or. buf(e:e) == '/' .or. buf(e:e) == '!' &
                      .or. buf(e:e) == new_line('a')) exit
                  e = e + 1
               enddo
               if (n + 2 <= len(buf)) then
                  buf(e+2:n+2) = buf(e:n)
                  buf(e+1:e+1) = "'"
                  buf(j+1:e) = buf(j:e-1)
                  buf(j:j) = "'"
               endif
               return
            endif
         endif
      endif
   endif
   i = i + 1
enddo
contains
   logical function word_char(c)
   character, intent(in) :: c
   word_char = (c >= 'a' .and. c <= 'z') .or. (c >= '0' .and. c <= '9') .or. c == '_'
   end function word_char
end subroutine quote_value

end module mod_engine_input_block
