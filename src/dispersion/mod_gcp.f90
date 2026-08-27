! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_gcp: Grimme's geometrical counterpoise (gCP) correction (H. Kruse,

module mod_gcp
use gcp, only: gcp_call
implicit none
private
public :: gcp_correction

contains

subroutine gcp_correction(natoms, atomic_number, coor_bohr, method, Egcp, dEgcp)
implicit none
integer,intent(in) :: natoms
integer,intent(in) :: atomic_number(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
character(*),intent(in) :: method
real(8),intent(out) :: Egcp
real(8),intent(out) :: dEgcp(3,natoms)

real(8) :: lat(3,3), glat(3,3)
character(len=20) :: method20
integer :: iz_work(natoms)

iz_work = atomic_number
lat = 0.0d0
method20 = trim(method)
call gcp_call(natoms, coor_bohr, lat, iz_work, Egcp, dEgcp, glat, &
              .true., .false., .false., method20, .false., .false.)

end subroutine gcp_correction

end module mod_gcp
