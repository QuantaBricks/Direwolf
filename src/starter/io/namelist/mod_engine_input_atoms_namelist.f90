! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! &atoms Mode A: classic namelist format (backward-compatible).

module mod_engine_input_atoms_namelist
implicit none
private
public :: parse_atoms_namelist_mode

contains

subroutine parse_atoms_namelist_mode(buf, ncenters, atomchg, x, y, z)
use mod_engine_input_elements, only: elem_normalize, elem_sym_to_z
character(len=*), intent(in)    :: buf
integer, intent(in)             :: ncenters
integer, intent(inout)          :: atomchg(:)
real(8), intent(inout)          :: x(:), y(:), z(:)

character(2) :: element(size(atomchg))
integer      :: i

namelist /atoms/ element, atomchg, x, y, z

element = ''
read(buf, nml=atoms)

if (element(1) /= '') then
   do i = 1, ncenters
      atomchg(i) = elem_sym_to_z(elem_normalize(element(i)))
      if (atomchg(i) == 0) then
         print *, 'Unknown element: ', trim(element(i))
         stop 1
      endif
   enddo
endif
end subroutine

end module mod_engine_input_atoms_namelist
