! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_radii: per-element cavity radius tables (plain Bondi*scale

module mod_cosmo_radii
use mod_cosmo_constants, only: bohr2ang
implicit none

real(8),parameter :: COSMO_DEFAULT_RADIUS_ANG = 2.00d0
integer,parameter :: N_BONDI = 54
real(8),parameter :: bondi_radius_ang(N_BONDI) = (/ &
   1.20d0, 1.40d0, 1.82d0, 1.53d0, 1.92d0, 1.70d0, 1.55d0, 1.52d0, 1.47d0, 1.54d0, &
   2.27d0, 1.73d0, 1.84d0, 2.10d0, 1.80d0, 1.80d0, 1.75d0, 1.88d0, 2.75d0, 2.31d0, &
   2.15d0, 2.11d0, 2.07d0, 2.06d0, 2.05d0, 2.04d0, 2.00d0, 1.97d0, 1.96d0, 2.01d0, &
   1.87d0, 2.11d0, 1.85d0, 1.90d0, 1.85d0, 2.02d0, 3.03d0, 2.49d0, 2.32d0, 2.23d0, &
   2.18d0, 2.17d0, 2.16d0, 2.13d0, 2.10d0, 1.63d0, 1.72d0, 1.58d0, 1.93d0, 2.17d0, &
   2.06d0, 2.06d0, 1.98d0, 2.16d0 /)

logical :: cosmo_radius_warned(200) = .false.

integer,parameter :: N_KSOPT = 53
real(8),parameter :: ks_optimized_radius_ang(N_KSOPT) = (/ &
   1.30d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.00d0, 1.83d0, 1.72d0, 1.72d0, 0.0d0, &
   0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.16d0, 2.05d0, 0.0d0, 0.0d0, 0.0d0, &
   0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
   0.0d0, 0.0d0, 0.0d0, 0.0d0, 2.16d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
   0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
   0.0d0, 0.0d0, 2.32d0 /)

contains

function cosmo_radius_bohr(Z, radii_scale) result(R_bohr)
integer,intent(in) :: Z
real(8),intent(in) :: radii_scale
real(8) :: R_bohr, R_ang
if (Z >= 1 .and. Z <= N_BONDI) then
   R_ang = bondi_radius_ang(Z)
else
   R_ang = COSMO_DEFAULT_RADIUS_ANG
   if (Z >= 1 .and. Z <= 200) then
      if (.not. cosmo_radius_warned(Z)) then
         print '("COSMO: no tabulated Bondi radius for Z=",I0,&
                 &" - using default ",F5.2," Ang (may reduce cavity accuracy)")', &
                 Z, COSMO_DEFAULT_RADIUS_ANG
         cosmo_radius_warned(Z) = .true.
      endif
   endif
endif
R_bohr = R_ang * radii_scale / bohr2ang
end function cosmo_radius_bohr

function cosmo_radius_ks_bohr(Z) result(R_bohr)
integer,intent(in) :: Z
real(8) :: R_bohr, R_ang
if (Z >= 1 .and. Z <= N_KSOPT) then
   if (ks_optimized_radius_ang(Z) > 0.0d0) then
      R_ang = ks_optimized_radius_ang(Z)
   else
      R_ang = 1.17d0*cosmo_radius_bohr(Z, 1.0d0)*bohr2ang
   endif
else
   R_ang = 1.17d0*cosmo_radius_bohr(Z, 1.0d0)*bohr2ang
endif
R_bohr = R_ang/bohr2ang
end function cosmo_radius_ks_bohr

function cosmo_smd_radius_bohr(Z, alpha) result(R_bohr)
integer,intent(in) :: Z
real(8),intent(in) :: alpha
real(8) :: R_bohr, R_ang
select case (Z)
case (1);  R_ang = 1.20d0
case (6);  R_ang = 1.85d0
case (7);  R_ang = 1.89d0
case (8)
   if (alpha >= 0.43d0) then
      R_ang = 1.52d0
   else
      R_ang = 1.52d0 + 1.8d0*(0.43d0 - alpha)
   endif
case (9);  R_ang = 1.73d0
case (14); R_ang = 2.47d0
case (15); R_ang = 2.12d0
case (16); R_ang = 2.49d0
case (17); R_ang = 2.38d0
case (35); R_ang = 2.60d0
case (53); R_ang = 2.74d0
case default
   R_ang = cosmo_radius_bohr(Z, 1.0d0)*bohr2ang
end select
R_bohr = R_ang/bohr2ang
end function cosmo_smd_radius_bohr

end module mod_cosmo_radii
