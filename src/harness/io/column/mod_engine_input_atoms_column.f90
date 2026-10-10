! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! &atoms Mode C: column format (new, recommended) - mixes optional

module mod_engine_input_atoms_column
implicit none
private
public :: is_old_namelist_atoms, parse_atoms_column, frozen_ranges, append_constraint

contains

logical function is_old_namelist_atoms(lines, start, finish)
character(len=*), intent(in) :: lines(:)
integer, intent(in) :: start, finish
character(len=256) :: line, key
integer :: i, eq
is_old_namelist_atoms = .false.
do i = start + 1, finish - 1
   line = trim(adjustl(lines(i)))
   if (line == '' .or. line(1:1) == '!' .or. line(1:1) == '#') cycle
   eq = index(line, '=')
   if (eq == 0) cycle
   key = trim(adjustl(line(:eq-1)))
   if (key == 'element' .or. key == 'atomchg' .or. &
       key == 'x' .or. key == 'y' .or. key == 'z') then
      is_old_namelist_atoms = .true.
      return
   endif
enddo
end function

subroutine parse_atoms_column(lines, start, finish, &
      icharge, imult, unit, baselabel, ecplabel, functional, &
      atomchg, x, y, z, atom_basis, atom_ecp, atom_frozen, na)
use mod_engine_input_elements, only: elem_normalize, elem_sym_to_z
character(len=*), intent(in)    :: lines(:)
integer, intent(in)             :: start, finish
integer, intent(inout)          :: icharge, imult
character(len=*), intent(inout) :: unit, baselabel, ecplabel, functional
integer, intent(out)            :: atomchg(:)
real(8), intent(out)            :: x(:), y(:), z(:)
character(len=*), intent(out)   :: atom_basis(:), atom_ecp(:)
logical, intent(out)            :: atom_frozen(:)
integer, intent(out)            :: na

character(len=256) :: line
character(len=64)  :: toks(8)
integer :: i, k, nopt, zval, ntok, ios
real(8) :: xv, yv, zv

na = 0
atom_frozen = .false.
do i = start + 1, finish - 1
   line = trim(adjustl(lines(i)))
   if (line == '' .or. line(1:1) == '!' .or. line(1:1) == '#') cycle
   if (index(line, '=') > 0) then
      call parse_atoms_keyval(line)
   else
      call split_tokens(line, toks, ntok)
      if (ntok < 4) then
         print *, 'Error parsing atom line in &atoms: ', trim(line)
         stop 1
      endif
      read(toks(2), *, iostat=ios) xv
      if (ios == 0) read(toks(3), *, iostat=ios) yv
      if (ios == 0) read(toks(4), *, iostat=ios) zv
      if (ios /= 0) then
         print *, 'Error parsing atom line in &atoms: ', trim(line)
         stop 1
      endif
      zval = elem_sym_to_z(elem_normalize(toks(1)))
      if (zval == 0) then
         print *, 'Unknown element: ', trim(toks(1))
         stop 1
      endif
      na = na + 1
      if (na > size(atomchg)) then
         print *, 'Too many atoms (max ', size(atomchg), ')'
         stop 1
      endif
      atomchg(na) = zval
      x(na) = xv; y(na) = yv; z(na) = zv
      atom_basis(na) = ''
      atom_ecp(na) = ''
      nopt = 0
      do k = 5, min(ntok, size(toks))
         if (trim(toks(k)) == 'F' .or. trim(toks(k)) == 'f') then
            atom_frozen(na) = .true.
            cycle
         endif
         nopt = nopt + 1
         if (trim(toks(k)) == '-') cycle
         if (nopt == 1) atom_basis(na) = toks(k)
         if (nopt == 2) atom_ecp(na) = toks(k)
      enddo
   endif
enddo

contains

   subroutine parse_atoms_keyval(line)
   character(len=*), intent(in) :: line
   character(len=256) :: key, val
   integer :: eq
   eq = index(line, '=')
   if (eq == 0) return
   key = trim(adjustl(line(:eq-1)))
   val = trim(adjustl(line(eq+1:)))
   if (len_trim(val) >= 2) then
      if ((val(1:1) == "'" .and. val(len_trim(val):len_trim(val)) == "'") .or. &
          (val(1:1) == '"' .and. val(len_trim(val):len_trim(val)) == '"')) then
         val = val(2:len_trim(val)-1)
      endif
   endif
   select case (key)
   case ('charge');       read(val, *) icharge
   case ('mult');         read(val, *) imult
   case ('unit');         unit = val
   case ('basis');        baselabel = val
   case ('ecp');          ecplabel = val
   case ('functional');   functional = val
   case default
      print *, 'Unknown key in &atoms: ', trim(key)
      stop 1
   end select
   end subroutine

end subroutine

function frozen_ranges(flags, n) result(str)
logical, intent(in) :: flags(:)
integer, intent(in) :: n
character(len=256) :: str
character(len=24) :: piece
integer :: i, j
str = ''
i = 1
do while (i <= n)
   if (.not. flags(i)) then
      i = i + 1
      cycle
   endif
   j = i
   do while (j < n)
      if (.not. flags(j+1)) exit
      j = j + 1
   enddo
   if (j > i) then
      write(piece,'(I0,A,I0)') i, '-', j
   else
      write(piece,'(I0)') i
   endif
   if (len_trim(str) + len_trim(piece) + 1 > len(str)) then
      print *, 'Too many separate frozen-atom ranges for opt_freeze (use opt_freeze = ... directly)'
      stop 1
   endif
   str = trim(str)//' '//trim(piece)
   i = j + 1
enddo
str = adjustl(str)
end function frozen_ranges

subroutine append_constraint(cons, nat_in_item, ia, has_val, val)
character(len=*), intent(inout) :: cons
integer, intent(in) :: nat_in_item, ia(4)
logical, intent(in) :: has_val
real(8), intent(in) :: val
character(len=64) :: piece
character(len=24) :: tmp
integer :: k
select case (nat_in_item)
case (2); piece = 'bond'
case (3); piece = 'angle'
case default; piece = 'dihedral'
end select
do k = 1, nat_in_item
   write(tmp,'(I0)') ia(k)
   piece = trim(piece)//' '//trim(tmp)
enddo
if (has_val) then
   write(tmp,'(F0.6)') val
   piece = trim(piece)//' '//trim(tmp)
endif
if (len_trim(cons) + len_trim(piece) + 2 > len(cons)) then
   print *, 'Too many constraints for opt_constraints'
   stop 1
endif
if (len_trim(cons) > 0) then
   cons = trim(cons)//'; '//trim(piece)
else
   cons = trim(piece)
endif
end subroutine append_constraint

subroutine split_tokens(line, tok, ntok)
character(len=*), intent(in) :: line
character(len=*), intent(out) :: tok(:)
integer, intent(out) :: ntok
integer :: i, len_s, tstart
logical :: in_tok

ntok = 0
in_tok = .false.
tstart = 1
len_s = len_trim(line)
do i = 1, len_s
   if (line(i:i) /= ' ' .and. line(i:i) /= char(9)) then
      if (.not. in_tok) then
         in_tok = .true.
         tstart = i
      endif
   else
      if (in_tok) then
         ntok = ntok + 1
         if (ntok <= size(tok)) tok(ntok) = line(tstart:i-1)
         in_tok = .false.
      endif
   endif
enddo
if (in_tok) then
   ntok = ntok + 1
   if (ntok <= size(tok)) tok(ntok) = line(tstart:len_s)
endif
end subroutine

end module mod_engine_input_atoms_column
