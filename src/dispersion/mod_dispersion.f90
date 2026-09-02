! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_dispersion: empirical (Grimme-style) dispersion correction, added

module mod_dispersion
implicit none
private
public :: dispersion_d2
public :: dispersion_chg

real(8),parameter :: RvdW_D2(55) = (/ &
   2.07869858790000d0, 1.89161571498900d0, 1.91240270086800d0, 1.55902394092500d0, &
   2.66073419251200d0, 2.80624309366500d0, 2.74388213602800d0, 2.63994720663300d0, &
   2.53601227723800d0, 2.43207734784300d0, 2.34892940432700d0, 2.16184653141600d0, &
   2.57758624899600d0, 3.09726089597100d0, 3.24276979712400d0, 3.22198281124500d0, &
   3.18040883948700d0, 3.09726089597100d0, 3.01411295245500d0, 2.80624309366500d0, &
   2.78545610778600d0, 2.95175199481800d0, 2.95175199481800d0, 2.95175199481800d0, &
   2.95175199481800d0, 2.95175199481800d0, 2.95175199481800d0, 2.95175199481800d0, &
   2.95175199481800d0, 2.95175199481800d0, 2.95175199481800d0, 3.11804788185000d0, &
   3.26355678300300d0, 3.32591774064000d0, 3.34670472651900d0, 3.30513075476100d0, &
   3.26355678300300d0, 3.07647391009200d0, 3.03489993833400d0, 3.09726089597100d0, &
   3.09726089597100d0, 3.09726089597100d0, 3.09726089597100d0, 3.09726089597100d0, &
   3.09726089597100d0, 3.09726089597100d0, 3.09726089597100d0, 3.09726089597100d0, &
   3.09726089597100d0, 3.15962185360800d0, 3.40906568415600d0, 3.55457458530900d0, &
   3.57536157118800d0, 3.57536157118800d0, 3.55457458530900d0 /)

real(8),parameter :: C6_D2(55) = (/ &
   0.0000000000000000d00, 2.4283353778422600d00, 1.3876202159098600d00, 2.7925856845186000d01, &
   2.7925856845186000d01, 5.4290640947473300d01, 3.0354192223028200d01, 2.1334660819614100d01, &
   1.2141676889211300d01, 1.3008939524155000d01, 1.0927509200290200d01, 9.9041392910566400d01, &
   9.9041392910566400d01, 1.8715527662084300d02, 1.6009668241060000d02, 1.3598678115916600d02, &
   9.6613057532724100d01, 8.7940431183287500d01, 7.9961614941805800d01, 1.8732872914783100d02, &
   1.8732872914783100d02, 1.8732872914783100d02, 1.8732872914783100d02, 1.8732872914783100d02, &
   1.8732872914783100d02, 1.8732872914783100d02, 1.8732872914783100d02, 1.8732872914783100d02, &
   1.8732872914783100d02, 1.8732872914783100d02, 1.8732872914783100d02, 2.9469584335385700d02, &
   2.9660382115073300d02, 2.8394178668055500d02, 2.1924399411375800d02, 2.1629530115495000d02, &
   2.0831648491346800d02, 4.2790738408120400d02, 4.2790738408120400d02, 4.2790738408120400d02, &
   4.2790738408120400d02, 4.2790738408120400d02, 4.2790738408120400d02, 4.2790738408120400d02, &
   4.2790738408120400d02, 4.2790738408120400d02, 4.2790738408120400d02, 4.2790738408120400d02, &
   4.2790738408120400d02, 6.4732483072195000d02, 6.7143473197338400d02, 6.6675151374468900d02, &
   5.5053832066223800d02, 5.4637546001450800d02, 5.2018412843920900d02 /)

real(8),parameter :: D2_ALPHA = 20.0d0

real(8),parameter :: CHG_ALPHA = 6.0d0

contains

subroutine dispersion_d2(natoms, atomic_number, coor_bohr, s6, Edisp, dEdisp)
implicit none
integer,intent(in) :: natoms
integer,intent(in) :: atomic_number(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(in) :: s6
real(8),intent(out) :: Edisp
real(8),intent(out) :: dEdisp(3,natoms)

integer :: i,j,zi,zj
real(8) :: dxyz(3), R2, R, R6, Rm6, R0, expo, fdamp, C6ij
real(8) :: dfdamp_dR, dRm6_dR, dE_dR

Edisp = 0.0d0
dEdisp = 0.0d0

do i = 1,natoms
   zi = atomic_number(i)
   if (zi .lt. 1 .or. zi .gt. 54) cycle
   do j = 1,i-1
      zj = atomic_number(j)
      if (zj .lt. 1 .or. zj .gt. 54) cycle

      dxyz = coor_bohr(:,i) - coor_bohr(:,j)
      R2 = sum(dxyz*dxyz)
      R = sqrt(R2)
      R6 = R2*R2*R2
      Rm6 = 1.0d0/R6

      C6ij = sqrt(C6_D2(zi+1)*C6_D2(zj+1))
      R0 = RvdW_D2(zi+1)+RvdW_D2(zj+1)
      expo = exp(-D2_ALPHA*(R/R0 - 1.0d0))
      fdamp = 1.0d0/(1.0d0+expo)

      Edisp = Edisp - s6*C6ij*Rm6*fdamp

      dRm6_dR = -6.0d0*Rm6/R
      dfdamp_dR = fdamp*fdamp*expo*(D2_ALPHA/R0)
      dE_dR = -s6*C6ij*(dRm6_dR*fdamp + Rm6*dfdamp_dR)

      dEdisp(:,i) = dEdisp(:,i) + dE_dR*(dxyz/R)
      dEdisp(:,j) = dEdisp(:,j) - dE_dR*(dxyz/R)
   enddo
enddo

end subroutine dispersion_d2

subroutine dispersion_chg(natoms, atomic_number, coor_bohr, Edisp, dEdisp)
implicit none
integer,intent(in) :: natoms
integer,intent(in) :: atomic_number(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(out) :: Edisp
real(8),intent(out) :: dEdisp(3,natoms)

integer :: i,j,zi,zj
real(8) :: dxyz(3), R2, R, R6, Rm6, R0, RR0m12, fdamp, C6ij
real(8) :: dfdamp_dR, dRm6_dR, dE_dR

Edisp = 0.0d0
dEdisp = 0.0d0

do i = 1,natoms
   zi = atomic_number(i)
   if (zi .lt. 1 .or. zi .gt. 54) cycle
   do j = 1,i-1
      zj = atomic_number(j)
      if (zj .lt. 1 .or. zj .gt. 54) cycle

      dxyz = coor_bohr(:,i) - coor_bohr(:,j)
      R2 = sum(dxyz*dxyz)
      R = sqrt(R2)
      R6 = R2*R2*R2
      Rm6 = 1.0d0/R6

      C6ij = sqrt(C6_D2(zi+1)*C6_D2(zj+1))
      R0 = RvdW_D2(zi+1)+RvdW_D2(zj+1)
      RR0m12 = (R/R0)**(-12.0d0)
      fdamp = 1.0d0/(1.0d0+CHG_ALPHA*RR0m12)

      Edisp = Edisp - C6ij*Rm6*fdamp

      dRm6_dR = -6.0d0*Rm6/R
      dfdamp_dR = fdamp*fdamp*12.0d0*CHG_ALPHA*RR0m12/R
      dE_dR = -C6ij*(dRm6_dR*fdamp + Rm6*dfdamp_dR)

      dEdisp(:,i) = dEdisp(:,i) + dE_dR*(dxyz/R)
      dEdisp(:,j) = dEdisp(:,j) - dE_dR*(dxyz/R)
   enddo
enddo

end subroutine dispersion_chg

end module mod_dispersion
