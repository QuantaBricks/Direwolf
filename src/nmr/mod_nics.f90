! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_nics: nucleus-independent chemical shifts. A probe is a point in space with no nucleus and no basis function; its GIAO

module mod_nics
implicit none
private
public :: nics_setup

integer,parameter :: MAXRING = 8

contains

function rcov(z) result(r)
integer,intent(in) :: z
real(8) :: r
real(8),parameter :: tab(36) = [0.31d0, 0.28d0, 1.28d0, 0.96d0, 0.84d0, 0.76d0, 0.71d0, 0.66d0, 0.57d0, 0.58d0, &
   1.66d0, 1.41d0, 1.21d0, 1.11d0, 1.07d0, 1.05d0, 1.02d0, 1.06d0, 2.03d0, 1.76d0, 1.70d0, 1.60d0, 1.53d0, 1.39d0, &
   1.39d0, 1.32d0, 1.26d0, 1.24d0, 1.32d0, 1.22d0, 1.22d0, 1.20d0, 1.19d0, 1.20d0, 1.20d0, 1.16d0]
r = 1.5d0
if (z >= 1 .and. z <= 36) r = tab(z)
end function rcov

subroutine nics_setup(req, nat, znum, coor_bohr, np, pos, label, normal, ring, ring_atoms, ring_size, nring, ok, msg)
character(len=*),intent(in) :: req
integer,intent(in) :: nat, znum(nat)
real(8),intent(in) :: coor_bohr(3,nat)
integer,intent(out) :: np, nring
real(8),allocatable,intent(out) :: pos(:,:), normal(:,:)
character(len=40),allocatable,intent(out) :: label(:)
integer,allocatable,intent(out) :: ring(:), ring_atoms(:,:), ring_size(:)
logical,intent(out) :: ok
character(len=*),intent(out) :: msg
real(8),parameter :: ANG = 1.0d0/0.529177210903d0
character(len=len(req)) :: r
ok = .true. ; msg = '' ; np = 0 ; nring = 0
r = adjustl(req)
if (trim(r) == 'rings') then
   call find_rings(nat, znum, coor_bohr, nring, ring_atoms, ring_size)
   if (nring == 0) then
      ok = .false. ; msg = "nics = 'rings': the molecule has no ring (3 to 8 members)" ; return
   endif
   np = 3*nring
   allocate(pos(3,np), normal(3,np), label(np), ring(np))
   call ring_probes()
else
   call parse_points()
endif

contains

subroutine ring_probes()
integer :: k, j, a, info
real(8) :: c(3), cov(3,3), ev(3), wk(20), n(3), d(3)
do k = 1,nring
   c = 0.0d0
   do j = 1,ring_size(k)
      c = c + coor_bohr(:,ring_atoms(j,k))
   enddo
   c = c/ring_size(k)
   cov = 0.0d0
   do j = 1,ring_size(k)
      d = coor_bohr(:,ring_atoms(j,k)) - c
      do a = 1,3
         cov(:,a) = cov(:,a) + d*d(a)
      enddo
   enddo
   call dsyev('V','U',3,cov,3,ev,wk,20,info)
   n = cov(:,1)
   n = n/sqrt(sum(n*n))
   pos(:,3*k-2) = c ; pos(:,3*k-1) = c + ANG*n ; pos(:,3*k) = c - ANG*n
   normal(:,3*k-2:3*k) = spread(n, 2, 3)
   ring(3*k-2:3*k) = k
   write(label(3*k-2),'(A,I0,A)') 'ring ', k, ' NICS(0)'
   write(label(3*k-1),'(A,I0,A)') 'ring ', k, ' NICS(+1)'
   write(label(3*k),'(A,I0,A)') 'ring ', k, ' NICS(-1)'
enddo
end subroutine ring_probes

subroutine parse_points()
integer :: i, j, n, ios
character(len=len(req)) :: s, item
real(8) :: v(3)
n = 1
do i = 1,len_trim(r)
   if (r(i:i) == ';') n = n + 1
enddo
np = n
allocate(pos(3,np), normal(3,np), label(np), ring(np))
s = r
do i = 1,np
   j = index(s, ';')
   if (j == 0) then
      item = s ; s = ''
   else
      item = s(1:j-1) ; s = s(j+1:)
   endif
   read(item, *, iostat=ios) v
   if (ios /= 0) then
      ok = .false. ; msg = "nics: use 'rings' or points 'x y z; x y z' in Angstrom" ; return
   endif
   pos(:,i) = v*ANG
   normal(:,i) = [0.0d0, 0.0d0, 1.0d0]
   ring(i) = 0
   write(label(i),'(A,I0)') 'point ', i
enddo
nring = 0
allocate(ring_atoms(MAXRING,1), ring_size(1))
ring_atoms = 0 ; ring_size = 0
end subroutine parse_points

end subroutine nics_setup

subroutine find_rings(nat, znum, coor_bohr, nring, ring_atoms, ring_size)
integer,intent(in) :: nat, znum(nat)
real(8),intent(in) :: coor_bohr(3,nat)
integer,intent(out) :: nring
integer,allocatable,intent(out) :: ring_atoms(:,:), ring_size(:)
real(8),parameter :: ANG = 1.0d0/0.529177210903d0
logical :: bond(nat,nat), dup
integer :: u, v, i, j, k, head, tail, cur, len_, qu(nat), prev(nat), dist(nat), path(MAXRING), cand(MAXRING,nat*4), csz(nat*4), nc
integer :: key(MAXRING), t
nc = 0 ; bond = .false.
do i = 1,nat
   if (znum(i) <= 1) cycle
   do j = i+1,nat
      if (znum(j) <= 1) cycle
      if (sqrt(sum((coor_bohr(:,i) - coor_bohr(:,j))**2))/ANG < 1.25d0*(rcov(znum(i)) + rcov(znum(j)))) then
         bond(i,j) = .true. ; bond(j,i) = .true.
      endif
   enddo
enddo
do u = 1,nat
   do v = u+1,nat
      if (.not. bond(u,v)) cycle
      dist = -1 ; prev = 0 ; head = 1 ; tail = 1 ; qu(1) = u ; dist(u) = 0
      do while (head <= tail)
         cur = qu(head) ; head = head + 1
         if (cur == v) exit
         if (dist(cur) >= MAXRING - 1) cycle
         do k = 1,nat
            if (.not. bond(cur,k) .or. dist(k) >= 0) cycle
            if (cur == u .and. k == v) cycle
            dist(k) = dist(cur) + 1 ; prev(k) = cur
            tail = tail + 1 ; qu(tail) = k
         enddo
      enddo
      if (dist(v) < 0 .or. dist(v) + 1 < 3 .or. dist(v) + 1 > MAXRING) cycle
      len_ = 0 ; cur = v
      do while (cur /= u)
         len_ = len_ + 1 ; path(len_) = cur ; cur = prev(cur)
      enddo
      len_ = len_ + 1 ; path(len_) = u
      key = 0 ; key(1:len_) = path(1:len_)
      do i = 2,len_
         t = key(i) ; j = i - 1
         do while (j >= 1)
            if (key(j) <= t) exit
            key(j+1) = key(j) ; j = j - 1
         enddo
         key(j+1) = t
      enddo
      dup = .false.
      do i = 1,nc
         if (csz(i) == len_) then
            if (all(cand(1:len_,i) == key(1:len_))) dup = .true.
         endif
      enddo
      if (dup) cycle
      nc = nc + 1 ; cand(:,nc) = key ; csz(nc) = len_
   enddo
enddo
nring = nc
allocate(ring_atoms(MAXRING,max(nring,1)), ring_size(max(nring,1)))
ring_atoms = 0 ; ring_size = 0
do i = 1,nring
   ring_atoms(:,i) = cand(:,i) ; ring_size(i) = csz(i)
enddo
end subroutine find_rings

end module mod_nics
