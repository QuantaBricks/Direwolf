! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_cosmo_cavity: cavity-build ORCHESTRATION only - assembles the

module mod_cosmo_cavity
use mod_cosmo_constants
use mod_cosmo_dual
use mod_cosmo_radii
use mod_cosmo_geodesic
use mod_cosmo_clip
use mod_cosmo_gaubon
implicit none

contains

subroutine cosmo_build_cavity(natoms, Z, coor_bohr, radii_scale, avg_area_ang2, &
                               tess_coor, tess_area, tess_atom, n_tess, Rvdw_override)
integer,intent(in) :: natoms
integer,intent(in) :: Z(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(in) :: radii_scale
real(8),intent(in) :: avg_area_ang2
real(8),allocatable,intent(out) :: tess_coor(:,:), tess_area(:)
integer,allocatable,intent(out) :: tess_atom(:)
integer,intent(out) :: n_tess
real(8),intent(in),optional :: Rvdw_override(natoms)

real(8),allocatable :: Rvdw(:)
real(8),allocatable :: tri(:,:,:)
real(8) :: verts(3,MAXV), ccc(3,MAXV), area, centroid(3)
integer :: intsph(MAXV), nv
integer :: ia, it, nf, n_tri, count_tess
logical :: discarded
real(8) :: target_area_bohr2

target_area_bohr2 = avg_area_ang2/(bohr2ang*bohr2ang)

allocate(Rvdw(natoms))
if (present(Rvdw_override)) then
   Rvdw = Rvdw_override
else
   do ia = 1,natoms
      Rvdw(ia) = cosmo_radius_bohr(Z(ia), radii_scale)
   enddo
endif

count_tess = 0
do ia = 1,natoms
   nf = cosmo_choose_nf(Rvdw(ia), target_area_bohr2)
   call cosmo_geodesic_unit_triangles(nf, tri, n_tri)
   do it = 1,n_tri
      block
         real(8) :: tri_scaled(3,3)
         integer :: k
         do k = 1,3
            tri_scaled(:,k) = coor_bohr(:,ia) + Rvdw(ia)*tri(:,k,it)
         enddo
         call cosmo_clip_tessera(tri_scaled, ia, natoms, Rvdw, coor_bohr, &
                                  verts, ccc, intsph, nv, discarded)
         if (.not. discarded) count_tess = count_tess+1
      end block
   enddo
   deallocate(tri)
enddo

n_tess = count_tess
allocate(tess_coor(3,n_tess), tess_area(n_tess), tess_atom(n_tess))

count_tess = 0
do ia = 1,natoms
   nf = cosmo_choose_nf(Rvdw(ia), target_area_bohr2)
   call cosmo_geodesic_unit_triangles(nf, tri, n_tri)
   do it = 1,n_tri
      block
         real(8) :: tri_scaled(3,3)
         integer :: k
         do k = 1,3
            tri_scaled(:,k) = coor_bohr(:,ia) + Rvdw(ia)*tri(:,k,it)
         enddo
         call cosmo_clip_tessera(tri_scaled, ia, natoms, Rvdw, coor_bohr, &
                                  verts, ccc, intsph, nv, discarded)
         if (.not. discarded) then
            call cosmo_gaubon_area(verts, ccc, intsph, nv, ia, Rvdw, coor_bohr, area, centroid)
            if (area > 0.0d0) then
               count_tess = count_tess+1
               tess_coor(:,count_tess) = centroid
               tess_area(count_tess) = area
               tess_atom(count_tess) = ia
            endif
         endif
      end block
   enddo
   deallocate(tri)
enddo
n_tess = count_tess

deallocate(Rvdw)
end subroutine cosmo_build_cavity
subroutine cosmo_build_cavity_deriv(natoms, Z, coor_bohr, radii_scale, avg_area_ang2, &
                                     tess_coor, tess_area, tess_atom, n_tess, &
                                     tess_nslot, tess_atom_of_slot, tess_dArea, tess_dCentroid, &
                                     Rvdw_override)
integer,intent(in) :: natoms
integer,intent(in) :: Z(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(in) :: radii_scale
real(8),intent(in) :: avg_area_ang2
real(8),allocatable,intent(out) :: tess_coor(:,:), tess_area(:)
integer,allocatable,intent(out) :: tess_atom(:)
integer,intent(out) :: n_tess
integer,allocatable,intent(out) :: tess_nslot(:)
integer,allocatable,intent(out) :: tess_atom_of_slot(:,:)
real(8),allocatable,intent(out) :: tess_dArea(:,:,:)
real(8),allocatable,intent(out) :: tess_dCentroid(:,:,:,:)
real(8),intent(in),optional :: Rvdw_override(natoms)

real(8),allocatable :: Rvdw(:)
real(8),allocatable :: tri(:,:,:)
real(8) :: verts(3,MAXV), ccc(3,MAXV), area, centroid(3)
real(8) :: vJacI(3,3,MAXV), vJacA(3,3,MAXV), vJacB(3,3,MAXV)
integer :: vAtomA(MAXV), vAtomB(MAXV)
real(8) :: cJacI(3,3,MAXV), cJacGov(3,3,MAXV)
integer :: intsph(MAXV), nv
integer :: ia, it, nf, n_tri, count_tess, n_max, nslot
logical :: discarded
real(8) :: target_area_bohr2
integer :: atom_of_slot(COSMO_MAXSLOT)
real(8) :: dArea_l(3,COSMO_MAXSLOT), dCentroid_l(3,3,COSMO_MAXSLOT)

target_area_bohr2 = avg_area_ang2/(bohr2ang*bohr2ang)

allocate(Rvdw(natoms))
if (present(Rvdw_override)) then
   Rvdw = Rvdw_override
else
   do ia = 1,natoms
      Rvdw(ia) = cosmo_radius_bohr(Z(ia), radii_scale)
   enddo
endif

n_max = 0
do ia = 1,natoms
   nf = cosmo_choose_nf(Rvdw(ia), target_area_bohr2)
   call cosmo_geodesic_unit_triangles(nf, tri, n_tri)
   n_max = n_max + n_tri
   deallocate(tri)
enddo

allocate(tess_coor(3,n_max), tess_area(n_max), tess_atom(n_max))
allocate(tess_nslot(n_max))
allocate(tess_atom_of_slot(COSMO_MAXSLOT,n_max))
allocate(tess_dArea(3,COSMO_MAXSLOT,n_max))
allocate(tess_dCentroid(3,3,COSMO_MAXSLOT,n_max))

count_tess = 0
do ia = 1,natoms
   nf = cosmo_choose_nf(Rvdw(ia), target_area_bohr2)
   call cosmo_geodesic_unit_triangles(nf, tri, n_tri)
   do it = 1,n_tri
      block
         real(8) :: tri_scaled(3,3)
         integer :: k
         do k = 1,3
            tri_scaled(:,k) = coor_bohr(:,ia) + Rvdw(ia)*tri(:,k,it)
         enddo
         call cosmo_clip_tessera_deriv(tri_scaled, ia, natoms, Rvdw, coor_bohr, &
                                        verts, ccc, intsph, nv, discarded, &
                                        vJacI, vAtomA, vJacA, vAtomB, vJacB, cJacI, cJacGov)
         if (.not. discarded) then
            call cosmo_gaubon_area_deriv(verts, ccc, intsph, nv, ia, Rvdw, coor_bohr, &
                                          vJacI, vAtomA, vJacA, vAtomB, vJacB, cJacI, cJacGov, &
                                          area, centroid, nslot, atom_of_slot, dArea_l, dCentroid_l)
            if (area > 0.0d0) then
               count_tess = count_tess+1
               tess_coor(:,count_tess) = centroid
               tess_area(count_tess) = area
               tess_atom(count_tess) = ia
               tess_nslot(count_tess) = nslot
               tess_atom_of_slot(:,count_tess) = atom_of_slot
               tess_dArea(:,:,count_tess) = dArea_l
               tess_dCentroid(:,:,:,count_tess) = dCentroid_l
            endif
         endif
      end block
   enddo
   deallocate(tri)
enddo
n_tess = count_tess

block
   real(8),allocatable :: tc(:,:), ta(:)
   integer,allocatable :: tat(:), tns(:), tas(:,:)
   real(8),allocatable :: tda(:,:,:), tdc(:,:,:,:)
   allocate(tc(3,n_tess), ta(n_tess), tat(n_tess), tns(n_tess))
   allocate(tas(COSMO_MAXSLOT,n_tess), tda(3,COSMO_MAXSLOT,n_tess), tdc(3,3,COSMO_MAXSLOT,n_tess))
   tc = tess_coor(:,1:n_tess)
   ta = tess_area(1:n_tess)
   tat = tess_atom(1:n_tess)
   tns = tess_nslot(1:n_tess)
   tas = tess_atom_of_slot(:,1:n_tess)
   tda = tess_dArea(:,:,1:n_tess)
   tdc = tess_dCentroid(:,:,:,1:n_tess)
   call move_alloc(tc, tess_coor)
   call move_alloc(ta, tess_area)
   call move_alloc(tat, tess_atom)
   call move_alloc(tns, tess_nslot)
   call move_alloc(tas, tess_atom_of_slot)
   call move_alloc(tda, tess_dArea)
   call move_alloc(tdc, tess_dCentroid)
end block

deallocate(Rvdw)
end subroutine cosmo_build_cavity_deriv
subroutine cosmo_build_cavity_ks(natoms, Z, coor_bohr, rsolv_ang, n_seg_target, n_face_target, &
                                  tess_coor, tess_area, tess_atom, n_tess)
integer,intent(in) :: natoms
integer,intent(in) :: Z(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(in) :: rsolv_ang
integer,intent(in) :: n_seg_target, n_face_target
real(8),allocatable,intent(out) :: tess_coor(:,:), tess_area(:)
integer,allocatable,intent(out) :: tess_atom(:)
integer,intent(out) :: n_tess

real(8),allocatable :: Rvdw(:), Rtest(:)
real(8),allocatable :: seg_unit(:,:), face_unit(:,:)
integer,allocatable :: face_to_seg(:), surv_count(:)
integer :: n_seg, n_face, nf_seg, nf_face
real(8) :: rsolv_bohr, face_pos(3), area_per_face, dbest, d2
integer :: ia, ja, iseg, ifac, ibest, count_tess
logical :: eliminated

rsolv_bohr = rsolv_ang/bohr2ang

nf_seg = cosmo_choose_nf_points(n_seg_target)
call cosmo_geodesic_unit_points(nf_seg, seg_unit, n_seg)
nf_face = cosmo_choose_nf_points(n_face_target)
call cosmo_geodesic_unit_points(nf_face, face_unit, n_face)

allocate(face_to_seg(n_face))
do ifac = 1,n_face
   ibest = 1
   dbest = sum((face_unit(:,ifac)-seg_unit(:,1))**2)
   do iseg = 2,n_seg
      d2 = sum((face_unit(:,ifac)-seg_unit(:,iseg))**2)
      if (d2 < dbest) then
         dbest = d2
         ibest = iseg
      endif
   enddo
   face_to_seg(ifac) = ibest
enddo

allocate(Rvdw(natoms), Rtest(natoms))
do ia = 1,natoms
   Rvdw(ia) = cosmo_radius_ks_bohr(Z(ia))
   Rtest(ia) = Rvdw(ia) + rsolv_bohr
enddo

allocate(tess_coor(3,natoms*n_seg), tess_area(natoms*n_seg), tess_atom(natoms*n_seg))
allocate(surv_count(n_seg))
count_tess = 0
do ia = 1,natoms
   surv_count = 0
   area_per_face = 4.0d0*PI_cosmo*Rvdw(ia)*Rvdw(ia)/real(n_face,8)
   do ifac = 1,n_face
      face_pos = coor_bohr(:,ia) + Rtest(ia)*face_unit(:,ifac)
      eliminated = .false.
      do ja = 1,natoms
         if (ja == ia) cycle
         if (sqrt(sum((face_pos-coor_bohr(:,ja))**2)) < Rtest(ja)) then
            eliminated = .true.
            exit
         endif
      enddo
      if (.not. eliminated) surv_count(face_to_seg(ifac)) = surv_count(face_to_seg(ifac)) + 1
   enddo
   do iseg = 1,n_seg
      if (surv_count(iseg) > 0) then
         count_tess = count_tess+1
         tess_coor(:,count_tess) = coor_bohr(:,ia) + Rvdw(ia)*seg_unit(:,iseg)
         tess_area(count_tess) = real(surv_count(iseg),8)*area_per_face
         tess_atom(count_tess) = ia
      endif
   enddo
enddo
n_tess = count_tess

block
   real(8),allocatable :: tc(:,:), ta(:)
   integer,allocatable :: tat(:)
   allocate(tc(3,n_tess), ta(n_tess), tat(n_tess))
   tc = tess_coor(:,1:n_tess)
   ta = tess_area(1:n_tess)
   tat = tess_atom(1:n_tess)
   call move_alloc(tc, tess_coor)
   call move_alloc(ta, tess_area)
   call move_alloc(tat, tess_atom)
end block

deallocate(Rvdw, Rtest, seg_unit, face_unit, face_to_seg, surv_count)
end subroutine cosmo_build_cavity_ks
subroutine cosmo_build_cavity_yk(natoms, Z, coor_bohr, n_seg_target, &
                                  tess_coor, tess_area, tess_atom, n_tess)
integer,intent(in) :: natoms
integer,intent(in) :: Z(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
integer,intent(in) :: n_seg_target
real(8),allocatable,intent(out) :: tess_coor(:,:), tess_area(:)
integer,allocatable,intent(out) :: tess_atom(:)
integer,intent(out) :: n_tess

real(8),allocatable :: Rvdw(:), rin(:), rout(:)
real(8),allocatable :: seg_unit(:,:)
integer :: n_seg, nf_seg
real(8) :: area_per_seg, seg_pos(3), xi, dist, r, window_frac
integer :: ia, ja, iseg, count_tess
real(8),parameter :: ALPHAI = 0.5d0
real(8),parameter :: GAMMAS = 1.0d0
integer,parameter :: MINBEM = 2
real(8),parameter :: XI_TOL = 1.0d-10

nf_seg = cosmo_choose_nf_points(n_seg_target)
call cosmo_geodesic_unit_points(nf_seg, seg_unit, n_seg)

allocate(Rvdw(natoms), rin(natoms), rout(natoms))
window_frac = GAMMAS*sqrt(0.25d0**MINBEM)
do ia = 1,natoms
   Rvdw(ia) = cosmo_radius_ks_bohr(Z(ia))
   rin(ia) = Rvdw(ia)*(1.0d0 - ALPHAI*window_frac)
   rout(ia) = Rvdw(ia)*(1.0d0 + (1.0d0-ALPHAI)*window_frac)
enddo

allocate(tess_coor(3,natoms*n_seg), tess_area(natoms*n_seg), tess_atom(natoms*n_seg))
count_tess = 0
do ia = 1,natoms
   area_per_seg = 4.0d0*PI_cosmo*Rvdw(ia)*Rvdw(ia)/real(n_seg,8)
   do iseg = 1,n_seg
      seg_pos = coor_bohr(:,ia) + Rvdw(ia)*seg_unit(:,iseg)
      xi = 1.0d0
      do ja = 1,natoms
         if (ja == ia) cycle
         dist = sqrt(sum((seg_pos-coor_bohr(:,ja))**2))
         if (dist <= rin(ja)) then
            xi = 0.0d0
            exit
         else if (dist < rout(ja)) then
            r = (dist-rin(ja))/(rout(ja)-rin(ja))
            xi = xi * r*r*r*(10.0d0-15.0d0*r+6.0d0*r*r)
         endif
      enddo
      if (xi > XI_TOL) then
         count_tess = count_tess+1
         tess_coor(:,count_tess) = seg_pos
         tess_area(count_tess) = xi*area_per_seg
         tess_atom(count_tess) = ia
      endif
   enddo
enddo
n_tess = count_tess

block
   real(8),allocatable :: tc(:,:), ta(:)
   integer,allocatable :: tat(:)
   allocate(tc(3,n_tess), ta(n_tess), tat(n_tess))
   tc = tess_coor(:,1:n_tess)
   ta = tess_area(1:n_tess)
   tat = tess_atom(1:n_tess)
   call move_alloc(tc, tess_coor)
   call move_alloc(ta, tess_area)
   call move_alloc(tat, tess_atom)
end block

deallocate(Rvdw, rin, rout, seg_unit)
end subroutine cosmo_build_cavity_yk
integer function cosmo_iswig_leb_order(Zi, Rvdw_ang) result(ord)
integer,intent(in) :: Zi
real(8),intent(in) :: Rvdw_ang
integer,parameter :: LEB_ORDERS_TAB(19) = &
   (/6,14,26,38,50,74,86,110,146,170,194,230,266,302,350,434,590,770,974/)
integer :: n_want, i
select case (Zi)
case (1); ord = 110
case (6); ord = 302
case (7); ord = 194
case (8); ord = 194
case default
   n_want = nint(5.0d0*4.0d0*PI_cosmo*Rvdw_ang**2)
   ord = LEB_ORDERS_TAB(size(LEB_ORDERS_TAB))
   do i = 1,size(LEB_ORDERS_TAB)
      if (LEB_ORDERS_TAB(i) >= n_want) then
         ord = LEB_ORDERS_TAB(i)
         exit
      endif
   enddo
end select
end function cosmo_iswig_leb_order

subroutine cosmo_build_cavity_iswig(natoms, Z, coor_bohr, charge_density_target, &
                                     tess_coor, tess_area, tess_atom, n_tess)
integer,intent(in) :: natoms
integer,intent(in) :: Z(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(in) :: charge_density_target
real(8),allocatable,intent(out) :: tess_coor(:,:), tess_area(:)
integer,allocatable,intent(out) :: tess_atom(:)
integer,intent(out) :: n_tess

real(8),parameter :: ZETA_GLOBAL = 4.902d0
real(8),allocatable :: Rvdw(:)
real(8),allocatable :: lx(:), ly(:), lz(:), lw(:)
integer :: n_seg, n_max
real(8) :: seg_pos(3), gtot, rij, garg, zeta_seg
integer :: ia, ja, iseg, count_tess, max_tess
real(8),parameter :: G_TOL = 1.0d-7

allocate(Rvdw(natoms))
do ia = 1,natoms
   if (Z(ia) == 1) then
      Rvdw(ia) = 1.10d0*1.2d0/bohr2ang
   else
      Rvdw(ia) = cosmo_radius_bohr(Z(ia), 1.2d0)
   endif
enddo

n_max = 0
do ia = 1,natoms
   n_max = n_max + cosmo_iswig_leb_order(Z(ia), Rvdw(ia)*bohr2ang)
enddo
allocate(tess_coor(3,n_max), tess_area(n_max), tess_atom(n_max))
allocate(lx(974), ly(974), lz(974), lw(974))

count_tess = 0
do ia = 1,natoms
   call cosmo_lebedev_points(cosmo_iswig_leb_order(Z(ia), Rvdw(ia)*bohr2ang), lx, ly, lz, lw, n_seg)
   do iseg = 1,n_seg
      seg_pos = coor_bohr(:,ia) + Rvdw(ia)*(/lx(iseg),ly(iseg),lz(iseg)/)
      zeta_seg = ZETA_GLOBAL/(Rvdw(ia)*sqrt(lw(iseg)*4.0d0*PI_cosmo))
      gtot = 1.0d0
      do ja = 1,natoms
         if (ja == ia) cycle
         rij = sqrt(sum((seg_pos-coor_bohr(:,ja))**2))
         garg = 1.0d0 - 0.5d0*(erf(zeta_seg*(Rvdw(ja)-rij)) + erf(zeta_seg*(Rvdw(ja)+rij)))
         if (garg < G_TOL) garg = 0.0d0
         gtot = gtot*garg
         if (gtot < G_TOL) exit
      enddo
      if (gtot > G_TOL) then
         count_tess = count_tess+1
         tess_coor(:,count_tess) = seg_pos
         tess_area(count_tess) = gtot*lw(iseg)*4.0d0*PI_cosmo*Rvdw(ia)*Rvdw(ia)
         tess_atom(count_tess) = ia
      endif
   enddo
enddo
n_tess = count_tess

block
   real(8),allocatable :: tc(:,:), ta(:)
   integer,allocatable :: tat(:)
   allocate(tc(3,n_tess), ta(n_tess), tat(n_tess))
   tc = tess_coor(:,1:n_tess)
   ta = tess_area(1:n_tess)
   tat = tess_atom(1:n_tess)
   call move_alloc(tc, tess_coor)
   call move_alloc(ta, tess_area)
   call move_alloc(tat, tess_atom)
end block

deallocate(Rvdw, lx, ly, lz, lw)
end subroutine cosmo_build_cavity_iswig

end module mod_cosmo_cavity
