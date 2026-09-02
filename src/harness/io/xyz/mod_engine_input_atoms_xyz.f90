! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! &atoms Mode B: external XYZ file, triggered by &molecule's xyzfile

module mod_engine_input_atoms_xyz
implicit none
private
public :: read_xyz

contains

subroutine read_xyz(filename, atomchg, x, y, z, na)
use mod_engine_input_elements, only: elem_normalize, elem_sym_to_z
character(len=*), intent(in) :: filename
integer, intent(out) :: atomchg(:)
real(8), intent(out) :: x(:), y(:), z(:)
integer, intent(out) :: na

character(len=256) :: elem_str
integer :: i, ios, zval
real(8) :: xv, yv, zv

open(unit=12, file=trim(filename), status='old', action='read', iostat=ios)
if (ios /= 0) then
   print *, 'Error: cannot open XYZ file: ', trim(filename)
   stop 1
endif

read(12, *, iostat=ios) na
if (ios /= 0 .or. na <= 0) then
   print *, 'Error reading natoms from XYZ file: ', trim(filename)
   stop 1
endif

read(12, '(A)', iostat=ios)
if (ios /= 0) then
   print *, 'Error reading comment line from XYZ file: ', trim(filename)
   stop 1
endif

do i = 1, na
   read(12, *, iostat=ios) elem_str, xv, yv, zv
   if (ios /= 0) then
      print *, 'Error reading atom', i, ' from XYZ file: ', trim(filename)
      stop 1
   endif
   zval = elem_sym_to_z(elem_normalize(elem_str))
   if (zval == 0) then
      print *, 'Unknown element in XYZ file: ', trim(elem_str)
      stop 1
   endif
   atomchg(i) = zval
   x(i) = xv; y(i) = yv; z(i) = zv
enddo
close(12)

print *, '[XYZ] Read ', na, ' atoms from ', trim(filename)
end subroutine

end module mod_engine_input_atoms_xyz
