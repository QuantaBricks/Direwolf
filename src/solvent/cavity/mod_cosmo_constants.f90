! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_constants: shared unit-conversion/geometry constants used by

module mod_cosmo_constants
implicit none

real(8),parameter :: bohr2ang = 0.52917721067d0
real(8),parameter :: PI_cosmo = 3.14159265358979324d0

integer,parameter :: MAXV = 12
real(8),parameter :: COSMO_EDGE_TOL = 1.0d-9
real(8),parameter :: COSMO_INTER_TOL = 1.0d-12
integer,parameter :: COSMO_INTER_MAXIT = 300

integer,parameter :: COSMO_MAXSLOT = 2*MAXV+1

end module mod_cosmo_constants
