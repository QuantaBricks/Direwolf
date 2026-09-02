! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_geodesic: sphere discretization utilities shared by the

module mod_cosmo_geodesic
use mod_cosmo_constants, only: PI_cosmo
implicit none

contains

subroutine cosmo_subdivide_face(v1, v2, v3, nf, tri)
real(8),intent(in) :: v1(3), v2(3), v3(3)
integer,intent(in) :: nf
real(8),intent(out) :: tri(3,3,nf*nf)
real(8) :: grid(3,0:nf,0:nf)
integer :: a, b, c, itri
real(8) :: g(3), gnorm

do a = 0,nf
   do b = 0,nf-a
      c = nf-a-b
      g = real(a,8)*v1 + real(b,8)*v2 + real(c,8)*v3
      gnorm = sqrt(sum(g*g))
      grid(:,a,b) = g/gnorm
   enddo
enddo

itri = 0
do a = 0,nf-1
   do b = 0,nf-1-a
      itri = itri+1
      tri(:,1,itri) = grid(:,a,b)
      tri(:,2,itri) = grid(:,a+1,b)
      tri(:,3,itri) = grid(:,a,b+1)
      if (a+b < nf-1) then
         itri = itri+1
         tri(:,1,itri) = grid(:,a+1,b)
         tri(:,2,itri) = grid(:,a+1,b+1)
         tri(:,3,itri) = grid(:,a,b+1)
      endif
   enddo
enddo
end subroutine cosmo_subdivide_face

subroutine cosmo_geodesic_unit_triangles(nf, tri, n_tri)
integer,intent(in) :: nf
real(8),allocatable,intent(out) :: tri(:,:,:)
integer,intent(out) :: n_tri
real(8) :: v1(3), v2(3), v3(3)
integer :: sx, sy, sz, iface
n_tri = 8*nf*nf
allocate(tri(3,3,n_tri))
iface = 0
do sx = -1,1,2
   do sy = -1,1,2
      do sz = -1,1,2
         v1 = (/real(sx,8), 0.0d0, 0.0d0/)
         v2 = (/0.0d0, real(sy,8), 0.0d0/)
         v3 = (/0.0d0, 0.0d0, real(sz,8)/)
         call cosmo_subdivide_face(v1, v2, v3, nf, tri(:,:,iface*nf*nf+1:(iface+1)*nf*nf))
         iface = iface+1
      enddo
   enddo
enddo
end subroutine cosmo_geodesic_unit_triangles

function cosmo_choose_nf(R_bohr, target_area_bohr2) result(nf)
real(8),intent(in) :: R_bohr, target_area_bohr2
integer :: nf, itsnum
itsnum = nint(4.0d0*PI_cosmo*R_bohr*R_bohr/target_area_bohr2)
nf = nint(sqrt(real(itsnum,8)/8.0d0))
if (nf < 2) nf = 2
end function cosmo_choose_nf

subroutine cosmo_subdivide_face_points(v1, v2, v3, nf, pts)
real(8),intent(in) :: v1(3), v2(3), v3(3)
integer,intent(in) :: nf
real(8),intent(out) :: pts(3,(nf+1)*(nf+2)/2)
integer :: a, b, c, ip
real(8) :: g(3), gnorm
ip = 0
do a = 0,nf
   do b = 0,nf-a
      c = nf-a-b
      g = real(a,8)*v1 + real(b,8)*v2 + real(c,8)*v3
      gnorm = sqrt(sum(g*g))
      ip = ip+1
      pts(:,ip) = g/gnorm
   enddo
enddo
end subroutine cosmo_subdivide_face_points

subroutine cosmo_geodesic_unit_points(nf, pts, n_pts)
integer,intent(in) :: nf
real(8),allocatable,intent(out) :: pts(:,:)
integer,intent(out) :: n_pts
real(8),allocatable :: raw(:,:), face_pts(:,:)
integer :: n_raw_max, n_raw, sx, sy, sz, i, j
real(8) :: v1(3), v2(3), v3(3)
logical :: dup
real(8),parameter :: DEDUP_TOL = 1.0d-10

n_raw_max = 8*(nf+1)*(nf+2)/2
allocate(raw(3,n_raw_max))
allocate(face_pts(3,(nf+1)*(nf+2)/2))
n_raw = 0
do sx = -1,1,2
   do sy = -1,1,2
      do sz = -1,1,2
         v1 = (/real(sx,8), 0.0d0, 0.0d0/)
         v2 = (/0.0d0, real(sy,8), 0.0d0/)
         v3 = (/0.0d0, 0.0d0, real(sz,8)/)
         call cosmo_subdivide_face_points(v1, v2, v3, nf, face_pts)
         do i = 1,size(face_pts,2)
            raw(:,n_raw+1) = face_pts(:,i)
            n_raw = n_raw+1
         enddo
      enddo
   enddo
enddo
deallocate(face_pts)

allocate(pts(3,n_raw))
n_pts = 0
do i = 1,n_raw
   dup = .false.
   do j = 1,n_pts
      if (sum((raw(:,i)-pts(:,j))**2) < DEDUP_TOL) then
         dup = .true.
         exit
      endif
   enddo
   if (.not. dup) then
      n_pts = n_pts+1
      pts(:,n_pts) = raw(:,i)
   endif
enddo
deallocate(raw)
end subroutine cosmo_geodesic_unit_points

function cosmo_choose_nf_points(n_target) result(nf)
integer,intent(in) :: n_target
integer :: nf
nf = nint(sqrt(max(1.0d0, real(n_target-2,8)/4.0d0)))
if (nf < 1) nf = 1
end function cosmo_choose_nf_points
subroutine cosmo_lebedev_points(n_requested, x, y, z, w, n_actual)
integer,intent(in) :: n_requested
real(8),intent(out) :: x(:), y(:), z(:), w(:)
integer,intent(out) :: n_actual
select case (n_requested)
case (6);   call LD0006(x,y,z,w,n_actual)
case (14);  call LD0014(x,y,z,w,n_actual)
case (26);  call LD0026(x,y,z,w,n_actual)
case (38);  call LD0038(x,y,z,w,n_actual)
case (50);  call LD0050(x,y,z,w,n_actual)
case (74);  call LD0074(x,y,z,w,n_actual)
case (86);  call LD0086(x,y,z,w,n_actual)
case (110); call LD0110(x,y,z,w,n_actual)
case (146); call LD0146(x,y,z,w,n_actual)
case (170); call LD0170(x,y,z,w,n_actual)
case (194); call LD0194(x,y,z,w,n_actual)
case (230); call LD0230(x,y,z,w,n_actual)
case (266); call LD0266(x,y,z,w,n_actual)
case (302); call LD0302(x,y,z,w,n_actual)
case (350); call LD0350(x,y,z,w,n_actual)
case (434); call LD0434(x,y,z,w,n_actual)
case (590); call LD0590(x,y,z,w,n_actual)
case (770); call LD0770(x,y,z,w,n_actual)
case (974); call LD0974(x,y,z,w,n_actual)
case default
   print '("COSMO: n_lebedev=",I0," not in the supported set - using 110")', n_requested
   call LD0110(x,y,z,w,n_actual)
end select
end subroutine cosmo_lebedev_points

end module mod_cosmo_geodesic
