! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_smd_tables: SMD CDS constants and per-solvent surface-tension

module mod_smd_tables
use mod_cosmo_cavity, only: cosmo_radius_bohr, bohr2ang
implicit none

real(8),parameter :: hartree2kcal = 627.509451d0
real(8),parameter :: PI_smd = 3.14159265358979d0

integer,parameter :: NATCNV(100) = (/ &
   1,0,0,0,0, 2,3,4,5,0, 0,0,0,11,9, 6,7,0,0,0, 0,0,0,0,0, &
   0,0,0,0,0, 0,0,0,0,8, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, &
   0,0,10,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, &
   0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0, 0,0,0,0,0 /)

contains

function smd_rkkval() result(RKKVAL)
real(8) :: RKKVAL(15,15)
RKKVAL = 0.0d0
RKKVAL(1,1:5)  = (/0.00d0,1.55d0,1.55d0,1.55d0,0.00d0/)
RKKVAL(1,6:10) = (/2.14d0,0.00d0,0.00d0,0.00d0,0.00d0/)
RKKVAL(2,1:5)  = (/1.55d0,1.84d0,1.84d0,1.84d0,1.84d0/)
RKKVAL(2,6:10) = (/2.20d0,2.10d0,2.30d0,2.20d0,2.60d0/)
RKKVAL(3,1:5)  = (/1.55d0,1.84d0,1.85d0,1.50d0,0.00d0/)
RKKVAL(4,1:5)  = (/1.55d0,1.84d0,1.50d0,2.75d0,0.00d0/)
RKKVAL(4,6:10) = (/1.71d0,0.00d0,0.00d0,2.10d0,0.00d0/)
RKKVAL(4,11)   = 2.10d0
RKKVAL(5,1:5)  = (/0.00d0,1.84d0,0.00d0,0.00d0,0.00d0/)
RKKVAL(6,1:5)  = (/2.14d0,2.20d0,0.00d0,1.71d0,0.00d0/)
RKKVAL(6,6:10) = (/2.75d0,0.00d0,0.00d0,2.50d0,0.00d0/)
RKKVAL(7,1:5)  = (/0.00d0,2.10d0,0.00d0,0.00d0,0.00d0/)
RKKVAL(8,1:5)  = (/0.00d0,2.30d0,0.00d0,0.00d0,0.00d0/)
RKKVAL(9,1:5)  = (/0.00d0,2.20d0,0.00d0,2.10d0,0.00d0/)
RKKVAL(9,6:10) = (/2.50d0,0.00d0,0.00d0,0.00d0,0.00d0/)
RKKVAL(10,1:5) = (/0.00d0,2.60d0,0.00d0,0.00d0,0.00d0/)
RKKVAL(11,4)   = 2.10d0
end function smd_rkkval

subroutine smd_cds_aq(SIGMA, HSIGMA, CSSIGM)
real(8),intent(out) :: SIGMA(150), HSIGMA(150), CSSIGM
SIGMA = 0.0d0
HSIGMA = 0.0d0
SIGMA(1) = 48.69d0
SIGMA(6) = 129.74d0
SIGMA(9) = 38.18d0
SIGMA(16) = -9.10d0
SIGMA(17) = 9.82d0
SIGMA(35) = -8.72d0
SIGMA(101) = -72.95d0
SIGMA(103) = 68.69d0
SIGMA(105) = -48.22d0
SIGMA(106) = 121.98d0
SIGMA(114) = 68.85d0
SIGMA(116) = 84.10d0
HSIGMA(6) = -60.77d0
CSSIGM = 0.0d0
end subroutine smd_cds_aq

subroutine smd_cds_naq(SIGMA, HSIGMA, CSSIGM, sol_alpha, sol_beta, sol_phi, sol_gamma, sol_psi, sol_n)
real(8),intent(out) :: SIGMA(150), HSIGMA(150), CSSIGM
real(8),intent(in) :: sol_alpha, sol_beta, sol_phi, sol_gamma, sol_psi, sol_n
real(8) :: SIGMA_N(150), SIGMA_A(150), SIGMA_B(150)
real(8) :: HSIGMA_N(150), HSIGMA_A(150), HSIGMA_B(150)
real(8) :: SIGMA_MOL(4)
integer :: i

SIGMA_N = 0.0d0; SIGMA_A = 0.0d0; SIGMA_B = 0.0d0
HSIGMA_N = 0.0d0; HSIGMA_A = 0.0d0; HSIGMA_B = 0.0d0

SIGMA_N(6) = 58.10d0;  SIGMA_N(7) = 32.62d0;  SIGMA_N(8) = -17.56d0
SIGMA_N(14) = -18.04d0
SIGMA_N(16) = -33.17d0; SIGMA_N(17) = -24.31d0
SIGMA_N(35) = -35.42d0
SIGMA_N(101) = -62.05d0; SIGMA_N(103) = -15.70d0
SIGMA_N(110) = -99.76d0

SIGMA_A(6) = 48.10d0; SIGMA_A(8) = 193.06d0
SIGMA_A(103) = 95.99d0; SIGMA_A(105) = -41.00d0
SIGMA_A(110) = 152.20d0

SIGMA_B(6) = 32.87d0; SIGMA_B(8) = -43.79d0
SIGMA_B(104) = -128.16d0; SIGMA_B(106) = 79.13d0

HSIGMA_N(6) = -36.37d0; HSIGMA_N(8) = -19.39d0

SIGMA_MOL = (/0.35d0, 0.00d0, -4.19d0, -6.68d0/)

do i = 1,150
   SIGMA(i)  = SIGMA_N(i)*sol_n  + SIGMA_A(i)*sol_alpha  + SIGMA_B(i)*sol_beta
   HSIGMA(i) = HSIGMA_N(i)*sol_n + HSIGMA_A(i)*sol_alpha + HSIGMA_B(i)*sol_beta
enddo

CSSIGM = SIGMA_MOL(1)*sol_gamma + SIGMA_MOL(2)*sol_beta**2 &
       + SIGMA_MOL(3)*sol_phi**2 + SIGMA_MOL(4)*sol_psi**2
end subroutine smd_cds_naq

function smd_sasa_radius_ang(Z) result(R_ang)
integer,intent(in) :: Z
real(8) :: R_ang
R_ang = cosmo_radius_bohr(Z, 1.0d0)*bohr2ang + 0.4d0
end function smd_sasa_radius_ang

end module mod_smd_tables
