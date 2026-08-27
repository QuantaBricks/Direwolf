! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_smd_cds: SMD's "CDS" (Cavitation, Dispersion, Solvent-structure)

module mod_smd_cds
use mod_smd_tables, only: hartree2kcal, smd_sasa_radius_ang, smd_cds_aq, smd_cds_naq
use mod_smd_sts, only: smd_sts_and_grad
use mod_smd_daareal, only: smd_daareal
use mod_cosmo_cavity, only: bohr2ang
implicit none
private
public :: smd_cds_energy_force

contains

subroutine smd_cds_energy_force(natoms, Z, coor_bohr, &
                                 sol_n, sol_alpha, sol_beta, sol_gamma, sol_phi, sol_psi, is_water, &
                                 E_cds, force_cds)
integer,intent(in) :: natoms, Z(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(in) :: sol_n, sol_alpha, sol_beta, sol_gamma, sol_phi, sol_psi
logical,intent(in) :: is_water
real(8),intent(out) :: E_cds
real(8),intent(out) :: force_cds(3,natoms)

real(8) :: coor_ang(3,natoms), RAD(natoms), RIJ(natoms,natoms), URIJ(3,natoms,natoms)
real(8) :: SIGMA(150), HSIGMA(150), CSSIGM
real(8) :: STS(natoms), DSTS(3,natoms,natoms)
real(8) :: AREA(natoms), DATAR(3,natoms,natoms)
real(8) :: AREA0, RAS2, GCDS
integer :: ia, k1, k2, NCROSS
integer,allocatable :: NC(:)
real(8),allocatable :: DAREA(:,:)
real(8) :: dcdsx, dcdsy, dcdsz

coor_ang = coor_bohr*bohr2ang

do ia = 1,natoms
   RAD(ia) = smd_sasa_radius_ang(Z(ia))
enddo

RIJ = 0.0d0; URIJ = 0.0d0
do k1 = 1,natoms
   do k2 = 1,natoms
      if (k1 == k2) cycle
      RIJ(k1,k2) = sqrt(sum((coor_ang(:,k1)-coor_ang(:,k2))**2))
      URIJ(:,k1,k2) = (coor_ang(:,k2)-coor_ang(:,k1))/RIJ(k1,k2)
   enddo
enddo

if (is_water) then
   call smd_cds_aq(SIGMA, HSIGMA, CSSIGM)
else
   call smd_cds_naq(SIGMA, HSIGMA, CSSIGM, sol_alpha, sol_beta, sol_phi, sol_gamma, sol_psi, sol_n)
endif

call smd_sts_and_grad(natoms, Z, RIJ, URIJ, SIGMA, HSIGMA, STS, DSTS)

allocate(NC(0:natoms), DAREA(3,0:natoms))
DATAR = 0.0d0
GCDS = 0.0d0
do k1 = 1,natoms
   call smd_daareal(natoms, RAD, RIJ, URIJ, k1, AREA0, NCROSS, NC, DAREA)
   RAS2 = RAD(k1)**2
   AREA(k1) = AREA0*RAS2
   GCDS = GCDS + AREA(k1)*(STS(k1)+CSSIGM)*0.001d0
   do ia = 0,NCROSS
      DATAR(:,NC(ia),k1) = DAREA(:,ia)*RAS2
   enddo
enddo
deallocate(NC, DAREA)

E_cds = GCDS/hartree2kcal

do ia = 1,natoms
   dcdsx = 0.0d0; dcdsy = 0.0d0; dcdsz = 0.0d0
   do k1 = 1,natoms
      dcdsx = dcdsx + (STS(k1)+CSSIGM)*DATAR(1,ia,k1)*0.001d0 + DSTS(1,ia,k1)*AREA(k1)*0.001d0
      dcdsy = dcdsy + (STS(k1)+CSSIGM)*DATAR(2,ia,k1)*0.001d0 + DSTS(2,ia,k1)*AREA(k1)*0.001d0
      dcdsz = dcdsz + (STS(k1)+CSSIGM)*DATAR(3,ia,k1)*0.001d0 + DSTS(3,ia,k1)*AREA(k1)*0.001d0
   enddo
   force_cds(1,ia) = dcdsx/hartree2kcal*bohr2ang
   force_cds(2,ia) = dcdsy/hartree2kcal*bohr2ang
   force_cds(3,ia) = dcdsz/hartree2kcal*bohr2ang
enddo
end subroutine smd_cds_energy_force

end module mod_smd_cds
