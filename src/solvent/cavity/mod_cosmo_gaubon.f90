! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_gaubon: exact Gauss-Bonnet spherical-polygon area/centroid

module mod_cosmo_gaubon
use mod_cosmo_constants, only: PI_cosmo, MAXV, COSMO_MAXSLOT
use mod_cosmo_dual
implicit none

contains

subroutine cosmo_gaubon_area(verts, ccc, intsph, nv, i_atom, Rvdw, coor_bohr, area, centroid)
real(8),intent(in) :: verts(3,MAXV), ccc(3,MAXV)
integer,intent(in) :: intsph(MAXV)
integer,intent(in) :: nv, i_atom
real(8),intent(in) :: Rvdw(:)
real(8),intent(in) :: coor_bohr(:,:)
real(8),intent(out) :: area
real(8),intent(out) :: centroid(3)

real(8) :: phin(MAXV), beta(MAXV), weight(0:MAXV)
real(8) :: sum1, sum2, v1(3), v2(3), d1(3), d2(3), costn, cosphin
real(8) :: p1(3), p2(3), p3(3), p3b(3), u1(3), u2(3), dn1, dn2, dnorm
integer :: n, n1, n0, n2, nother
real(8) :: pp(3)

sum1 = 0.0d0
do n = 1,nv
   n2 = n+1; if (n == nv) n2 = 1
   v1 = verts(:,n) - ccc(:,n)
   v2 = verts(:,n2) - ccc(:,n)
   cosphin = sum(v1*v2)/(sqrt(sum(v1*v1))*sqrt(sum(v2*v2)))
   if (cosphin > 1.0d0) cosphin = 1.0d0
   if (cosphin < -1.0d0) cosphin = -1.0d0
   phin(n) = acos(cosphin)

   nother = intsph(n)
   d1 = coor_bohr(:,nother) - coor_bohr(:,i_atom)
   dn1 = sqrt(sum(d1*d1))
   if (dn1 <= 1.0d-14) dn1 = 1.0d0
   d2 = verts(:,n) - coor_bohr(:,i_atom)
   dn2 = sqrt(sum(d2*d2))
   costn = sum(d1*d2)/(dn1*dn2)
   sum1 = sum1 + phin(n)*costn
enddo

weight(1:nv) = phin(1:nv)
weight(0) = weight(nv)
do n = nv,1,-1
   weight(n) = weight(n) + weight(n-1)
enddo

sum2 = 0.0d0
do n1 = 1,nv
   n0 = n1-1; if (n0 == 0) n0 = nv
   n2 = n1+1; if (n2 > nv) n2 = 1

   p1 = verts(:,n1) - ccc(:,n0)
   p2 = verts(:,n0) - ccc(:,n0)
   call cosmo_cross(p1, p2, p3)
   call cosmo_cross(p1, p3, p3b)
   u1 = p3b/sqrt(sum(p3b*p3b))

   p1 = verts(:,n1) - ccc(:,n1)
   p2 = verts(:,n2) - ccc(:,n1)
   call cosmo_cross(p1, p2, p3)
   call cosmo_cross(p1, p3, p3b)
   u2 = p3b/sqrt(sum(p3b*p3b))

   beta(n1) = acos(max(-1.0d0,min(1.0d0,sum(u1*u2))))
   sum2 = sum2 + beta(n1)
enddo

area = Rvdw(i_atom)*Rvdw(i_atom)*(real(2-nv,8)*PI_cosmo + sum1 + sum2)
if (area < 0.0d0) area = 0.0d0

pp = 0.0d0
do n = 1,nv
   pp = pp + (verts(:,n)-coor_bohr(:,i_atom))*weight(n)
enddo
dnorm = sqrt(sum(pp*pp))
centroid = coor_bohr(:,i_atom) + pp*Rvdw(i_atom)/dnorm
end subroutine cosmo_gaubon_area

subroutine cosmo_gaubon_area_deriv(verts, ccc, intsph, nv, i_atom, Rvdw, coor_bohr, &
                                    vJacI, vAtomA, vJacA, vAtomB, vJacB, cJacI, cJacGov, &
                                    area, centroid, nslots, atom_of_slot, dArea, dCentroid)
real(8),intent(in) :: verts(3,MAXV), ccc(3,MAXV)
integer,intent(in) :: intsph(MAXV)
integer,intent(in) :: nv, i_atom
real(8),intent(in) :: Rvdw(:)
real(8),intent(in) :: coor_bohr(:,:)
real(8),intent(in) :: vJacI(3,3,MAXV), vJacA(3,3,MAXV), vJacB(3,3,MAXV)
integer,intent(in) :: vAtomA(MAXV), vAtomB(MAXV)
real(8),intent(in) :: cJacI(3,3,MAXV), cJacGov(3,3,MAXV)
real(8),intent(out) :: area
real(8),intent(out) :: centroid(3)
integer,intent(out) :: nslots
integer,intent(out) :: atom_of_slot(COSMO_MAXSLOT)
real(8),intent(out) :: dArea(3,COSMO_MAXSLOT)
real(8),intent(out) :: dCentroid(3,3,COSMO_MAXSLOT)

type(dvec3) :: vD(MAXV), cD(MAXV), CiD
type(dvec3) :: v1D, v2D, d1D, d2D, p1D, p2D, p3D, p3bD, u1D, u2D, ppD, centroidD
type(dscal) :: phinD(MAXV), betaD(MAXV), weightD(0:MAXV)
type(dscal) :: sum1D, sum2D, cosphinD, costnD, dn1D, dn2D, dnormD, areaD
integer :: n, n1, n0, n2, nother, slotI, slotA, slotB, slotGov, slotOther

nslots = 0
atom_of_slot = 0
slotI = cosmo_get_slot(atom_of_slot, nslots, i_atom)

do n = 1,nv
   slotA = 0
   if (vAtomA(n) /= 0) slotA = cosmo_get_slot(atom_of_slot, nslots, vAtomA(n))
   slotB = 0
   if (vAtomB(n) /= 0) slotB = cosmo_get_slot(atom_of_slot, nslots, vAtomB(n))
   vD(n) = const_dvec3(verts(:,n))
   vD(n)%g(:,:,slotI) = vJacI(:,:,n)
   if (slotA > 0) vD(n)%g(:,:,slotA) = vD(n)%g(:,:,slotA) + vJacA(:,:,n)
   if (slotB > 0) vD(n)%g(:,:,slotB) = vD(n)%g(:,:,slotB) + vJacB(:,:,n)

   slotGov = 0
   if (intsph(n) /= i_atom) slotGov = cosmo_get_slot(atom_of_slot, nslots, intsph(n))
   cD(n) = const_dvec3(ccc(:,n))
   cD(n)%g(:,:,slotI) = cJacI(:,:,n)
   if (slotGov > 0) cD(n)%g(:,:,slotGov) = cD(n)%g(:,:,slotGov) + cJacGov(:,:,n)
enddo

CiD = const_dvec3(coor_bohr(:,i_atom))
CiD%g(:,:,slotI) = cosmo_idn3()

sum1D = const_dscal(0.0d0)
do n = 1,nv
   n2 = n+1; if (n == nv) n2 = 1
   v1D = vD(n) - cD(n)
   v2D = vD(n2) - cD(n)
   cosphinD = ddot(v1D,v2D) / (dsqrt_d(ddot(v1D,v1D)) * dsqrt_d(ddot(v2D,v2D)))
   phinD(n) = dacos_d(cosphinD)

   nother = intsph(n)
   slotOther = slotI
   if (nother /= i_atom) slotOther = cosmo_get_slot(atom_of_slot, nslots, nother)
   d1D = const_dvec3(coor_bohr(:,nother) - coor_bohr(:,i_atom))
   d1D%g(:,:,slotOther) = d1D%g(:,:,slotOther) + cosmo_idn3()
   d1D%g(:,:,slotI) = d1D%g(:,:,slotI) - cosmo_idn3()
   dn1D = dsqrt_d(ddot(d1D,d1D))
   if (dn1D%v <= 1.0d-14) dn1D = const_dscal(1.0d0)
   d2D = vD(n) - CiD
   dn2D = dsqrt_d(ddot(d2D,d2D))
   costnD = ddot(d1D,d2D) / (dn1D*dn2D)
   sum1D = sum1D + phinD(n)*costnD
enddo

weightD(1:nv) = phinD(1:nv)
weightD(0) = weightD(nv)
do n = nv,1,-1
   weightD(n) = weightD(n) + weightD(n-1)
enddo

sum2D = const_dscal(0.0d0)
do n1 = 1,nv
   n0 = n1-1; if (n0 == 0) n0 = nv
   n2 = n1+1; if (n2 > nv) n2 = 1

   p1D = vD(n1) - cD(n0)
   p2D = vD(n0) - cD(n0)
   p3D = dcross(p1D,p2D)
   p3bD = dcross(p1D,p3D)
   u1D = p3bD / dsqrt_d(ddot(p3bD,p3bD))

   p1D = vD(n1) - cD(n1)
   p2D = vD(n2) - cD(n1)
   p3D = dcross(p1D,p2D)
   p3bD = dcross(p1D,p3D)
   u2D = p3bD / dsqrt_d(ddot(p3bD,p3bD))

   betaD(n1) = dacos_d(ddot(u1D,u2D))
   sum2D = sum2D + betaD(n1)
enddo

areaD = (Rvdw(i_atom)*Rvdw(i_atom)) * (const_dscal(real(2-nv,8)*PI_cosmo) + sum1D + sum2D)
area = areaD%v
if (area < 0.0d0) area = 0.0d0

ppD = const_dvec3((/0.0d0,0.0d0,0.0d0/))
do n = 1,nv
   ppD = ppD + (vD(n)-CiD)*weightD(n)
enddo
dnormD = dsqrt_d(ddot(ppD,ppD))
centroidD = CiD + ppD*(Rvdw(i_atom)/dnormD)
centroid = centroidD%v

dArea(:,1:nslots) = areaD%g(:,1:nslots)
dCentroid(:,:,1:nslots) = centroidD%g(:,:,1:nslots)
end subroutine cosmo_gaubon_area_deriv

end module mod_cosmo_gaubon
