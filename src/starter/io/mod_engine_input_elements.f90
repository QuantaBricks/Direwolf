! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Element symbol <-> atomic number, shared by all three &atoms input

module mod_engine_input_elements
implicit none
private
public :: elem_normalize, elem_sym_to_z, elem_to_name

character(len=2), parameter :: elem_syms(86) = [ &
  'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne', &
  'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca', &
  'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
  'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr', &
  'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
  'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd', &
  'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
  'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
  'Tl','Pb','Bi','Po','At','Rn']

contains

character(len=2) function elem_normalize(symbol)
character(len=*), intent(in) :: symbol
character(len=30) :: s
integer :: i, l
s = adjustl(symbol)
l = len_trim(s)
if (l >= 1 .and. s(1:1) >= 'a' .and. s(1:1) <= 'z') &
   s(1:1) = char(ichar(s(1:1)) - 32)
do i = 2, l
   if (s(i:i) >= 'A' .and. s(i:i) <= 'Z') &
      s(i:i) = char(ichar(s(i:i)) + 32)
enddo
elem_normalize = trim(s)
end function

integer function elem_sym_to_z(symbol)
character(len=*), intent(in) :: symbol
integer :: i
elem_sym_to_z = 0
do i = 1, 86
   if (trim(symbol) == trim(elem_syms(i))) then
      elem_sym_to_z = i
      return
   endif
enddo
end function

character(len=2) function elem_to_name(z)
integer, intent(in) :: z
if (z >= 1 .and. z <= 86) then
   elem_to_name = elem_syms(z)
else
   elem_to_name = '? '
endif
end function

end module mod_engine_input_elements
