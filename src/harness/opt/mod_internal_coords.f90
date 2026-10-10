! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Redundant internal coordinates (RIC) for geometry optimization.

module mod_internal_coords
use mod_wilson_bvec, only: bvec_stretch, bvec_bend, bvec_torsion, bvec_linbend, &
                            val_stretch, val_bend, val_torsion, val_linbend, &
                            wrap_pi, perp_pair
use mod_lindh_hessian, only: lindh_rho, lindh_row, &
                             K_STRETCH, K_BEND, K_TORSION
implicit none
private
public :: ic_set_t, ic_build, ic_rebuild_if_changed, ic_nprim
public :: ic_values, ic_bmatrix, ic_ginv, ic_delta_q
public :: ic_grad_to_internal, ic_back_transform, ic_diag_hessian
public :: ic_lindh_hessian, ic_rank, ic_write

integer, parameter :: KIND_STRETCH = 1, KIND_BEND = 2, KIND_TORSION = 3
integer, parameter :: KIND_LINBEND = 4
integer, parameter :: KIND_IMPROPER = 5
real(8), parameter :: BOND_FAC   = 1.30d0
real(8), parameter :: AUX_BRIDGE_FAC = 2.20d0
real(8), parameter :: LINEAR_TOL = 2.9670d0
real(8), parameter :: LINEAR_HYST = 0.0873d0
real(8), parameter :: ANG2BOHR   = 1.88972612456506d0
integer, parameter :: MAX_PRIM   = 40000
integer, public :: ic_n_dropped = 0

type :: ic_set_t
   integer :: nprim = 0
   integer :: nat   = 0
   integer, allocatable :: kind(:)
   integer, allocatable :: at(:,:)
   integer, allocatable :: Z(:)
   real(8), allocatable :: uref(:,:)
   integer :: nbond = 0
   integer, allocatable :: bond(:,:)
end type ic_set_t

contains

integer function ic_nprim(ic)
type(ic_set_t), intent(in) :: ic
ic_nprim = ic%nprim
end function ic_nprim

real(8) function rcov_ang(Z)
integer, intent(in) :: Z
real(8), parameter :: r(86) = [ &
  0.31d0,0.28d0,1.28d0,0.96d0,0.84d0,0.76d0,0.71d0,0.66d0,0.57d0,0.58d0, &
  1.66d0,1.41d0,1.21d0,1.11d0,1.07d0,1.05d0,1.02d0,1.06d0,2.03d0,1.76d0, &
  1.70d0,1.60d0,1.53d0,1.39d0,1.50d0,1.42d0,1.38d0,1.24d0,1.32d0,1.22d0, &
  1.22d0,1.20d0,1.19d0,1.20d0,1.20d0,1.16d0,2.20d0,1.95d0,1.90d0,1.75d0, &
  1.64d0,1.54d0,1.47d0,1.46d0,1.42d0,1.39d0,1.45d0,1.44d0,1.42d0,1.39d0, &
  1.39d0,1.38d0,1.39d0,1.40d0,2.44d0,2.15d0,2.07d0,2.04d0,2.03d0,2.01d0, &
  1.99d0,1.98d0,1.98d0,1.96d0,1.94d0,1.92d0,1.92d0,1.89d0,1.90d0,1.87d0, &
  1.87d0,1.75d0,1.70d0,1.62d0,1.51d0,1.44d0,1.41d0,1.36d0,1.36d0,1.32d0, &
  1.45d0,1.46d0,1.48d0,1.40d0,1.50d0,1.50d0 ]
if (Z >= 1 .and. Z <= 86) then
   rcov_ang = r(Z)
else
   rcov_ang = 1.50d0
endif
end function rcov_ang

subroutine build_bonds(nat, Z, coord, bond, nbond)
integer, intent(in)  :: nat, Z(nat)
real(8), intent(in)  :: coord(3,nat)
integer, allocatable, intent(out) :: bond(:,:)
integer, intent(out) :: nbond
integer :: i, j, cap, comp(nat), comp0(nat), nc, ci, cj, bi, bj
integer :: nextra
real(8) :: d, cut, dmin
integer, allocatable :: tmp(:,:)

cap = 32*nat + 64
allocate(tmp(2,cap))
nbond = 0
do i = 1, nat
   do j = i+1, nat
      d   = sqrt(sum((coord(:,i) - coord(:,j))**2))
      cut = BOND_FAC * (rcov_ang(Z(i)) + rcov_ang(Z(j))) * ANG2BOHR
      if (d < cut) then
         if (nbond < cap) then
            nbond = nbond + 1
            tmp(1,nbond) = i ; tmp(2,nbond) = j
         else
            ic_n_dropped = ic_n_dropped + 1
         endif
      endif
   enddo
enddo

do i = 1, nat
   comp(i) = i
enddo
do i = 1, nbond
   call uf_union(comp, nat, tmp(1,i), tmp(2,i))
enddo

comp0 = comp
nextra = 0
do i = 1, nat
   do j = i+1, nat
      if (uf_find(comp0, nat, i) == uf_find(comp0, nat, j)) cycle
      d   = sqrt(sum((coord(:,i) - coord(:,j))**2))
      cut = AUX_BRIDGE_FAC * (rcov_ang(Z(i)) + rcov_ang(Z(j))) * ANG2BOHR
      if (d < cut .and. nbond < cap) then
         nbond = nbond + 1
         tmp(1,nbond) = i ; tmp(2,nbond) = j
         nextra = nextra + 1
         call uf_union(comp, nat, i, j)
      endif
   enddo
enddo

do
   nc = uf_ncomp(comp, nat)
   if (nc <= 1) exit
   dmin = huge(1.0d0) ; bi = 0 ; bj = 0
   do i = 1, nat
      do j = i+1, nat
         ci = uf_find(comp, nat, i) ; cj = uf_find(comp, nat, j)
         if (ci == cj) cycle
         d = sqrt(sum((coord(:,i) - coord(:,j))**2))
         if (d < dmin) then
            dmin = d ; bi = i ; bj = j
         endif
      enddo
   enddo
   if (bi == 0) exit
   if (nbond < cap) then
      nbond = nbond + 1
      tmp(1,nbond) = bi ; tmp(2,nbond) = bj
   endif
   call uf_union(comp, nat, bi, bj)
enddo

allocate(bond(2,max(nbond,1)))
if (nbond > 0) bond(:,1:nbond) = tmp(:,1:nbond)
end subroutine build_bonds

recursive integer function uf_find(p, n, a) result(r)
integer, intent(inout) :: p(n)
integer, intent(in)    :: n, a
if (p(a) /= a) p(a) = uf_find(p, n, p(a))
r = p(a)
end function uf_find

subroutine uf_union(p, n, a, b)
integer, intent(inout) :: p(n)
integer, intent(in)    :: n, a, b
integer :: ra, rb
ra = uf_find(p, n, a) ; rb = uf_find(p, n, b)
if (ra /= rb) p(ra) = rb
end subroutine uf_union

integer function uf_ncomp(p, n)
integer, intent(inout) :: p(n)
integer, intent(in)    :: n
integer :: i
uf_ncomp = 0
do i = 1, n
   if (uf_find(p, n, i) == i) uf_ncomp = uf_ncomp + 1
enddo
end function uf_ncomp

subroutine ic_build(ic, nat, Z, coord)
type(ic_set_t), intent(out) :: ic
integer, intent(in) :: nat, Z(nat)
real(8), intent(in) :: coord(3,nat)

integer, allocatable :: bond(:,:)
integer :: nbond, i, j, k, l, b, nb, np
integer :: deg(nat)
integer, allocatable :: nbr(:,:)
integer, allocatable :: pk(:), pa(:,:)
real(8), allocatable :: pu(:,:)
real(8) :: ang, u1(3), u2(3)

ic_n_dropped = 0
call build_bonds(nat, Z, coord, bond, nbond)

deg = 0
do b = 1, nbond
   deg(bond(1,b)) = deg(bond(1,b)) + 1
   deg(bond(2,b)) = deg(bond(2,b)) + 1
enddo
nb = max(1, maxval(deg))
allocate(nbr(nb,nat))
deg = 0
do b = 1, nbond
   i = bond(1,b) ; j = bond(2,b)
   deg(i) = deg(i) + 1 ; nbr(deg(i),i) = j
   deg(j) = deg(j) + 1 ; nbr(deg(j),j) = i
enddo

allocate(pk(MAX_PRIM), pa(4,MAX_PRIM), pu(3,MAX_PRIM))
pu = 0.0d0
np = 0

do b = 1, nbond
   call add_prim(pk, pa, np, KIND_STRETCH, bond(1,b), bond(2,b), 0, 0)
enddo

do j = 1, nat
   do i = 1, deg(j)
      do k = i+1, deg(j)
         ang = val_bend(coord(:,nbr(i,j)), coord(:,j), coord(:,nbr(k,j)))
         if (ang > LINEAR_TOL) then
            call perp_pair(coord(:,nbr(k,j)) - coord(:,nbr(i,j)), u1, u2)
            call add_prim(pk, pa, np, KIND_LINBEND, nbr(i,j), j, nbr(k,j), 0, pu, u1)
            call add_prim(pk, pa, np, KIND_LINBEND, nbr(i,j), j, nbr(k,j), 0, pu, u2)
            cycle
         endif
         call add_prim(pk, pa, np, KIND_BEND, nbr(i,j), j, nbr(k,j), 0)
      enddo
   enddo
enddo

do j = 1, nat
   if (deg(j) /= 3) cycle
   call add_prim(pk, pa, np, KIND_IMPROPER, nbr(1,j), nbr(2,j), j, nbr(3,j))
enddo

do b = 1, nbond
   j = bond(1,b) ; k = bond(2,b)
   do i = 1, deg(j)
      if (nbr(i,j) == k) cycle
      if (val_bend(coord(:,nbr(i,j)), coord(:,j), coord(:,k)) > LINEAR_TOL) cycle
      do l = 1, deg(k)
         if (nbr(l,k) == j .or. nbr(l,k) == nbr(i,j)) cycle
         if (val_bend(coord(:,j), coord(:,k), coord(:,nbr(l,k))) > LINEAR_TOL) cycle
         call add_prim(pk, pa, np, KIND_TORSION, nbr(i,j), j, k, nbr(l,k))
      enddo
   enddo
enddo

ic%nat   = nat
ic%nprim = np
allocate(ic%Z(nat))
ic%Z = Z
allocate(ic%kind(np), ic%at(4,np), ic%uref(3,np))
ic%kind      = pk(1:np)
ic%at(:,1:np) = pa(:,1:np)
ic%uref(:,1:np) = pu(:,1:np)
ic%nbond = nbond
allocate(ic%bond(2,max(nbond,1)))
if (nbond > 0) ic%bond(:,1:nbond) = bond(:,1:nbond)
end subroutine ic_build

subroutine add_prim(pk, pa, np, kd, a1, a2, a3, a4, pu, uvec)
integer, intent(inout) :: pk(:), pa(:,:), np
integer, intent(in)    :: kd, a1, a2, a3, a4
real(8), intent(inout), optional :: pu(:,:)
real(8), intent(in),    optional :: uvec(3)
if (np >= size(pk)) then
   ic_n_dropped = ic_n_dropped + 1
   return
endif
np = np + 1
pk(np) = kd
pa(1,np) = a1 ; pa(2,np) = a2 ; pa(3,np) = a3 ; pa(4,np) = a4
if (present(pu) .and. present(uvec)) pu(:,np) = uvec
end subroutine add_prim

subroutine ic_rebuild_if_changed(ic, nat, Z, coord, changed)
type(ic_set_t), intent(inout) :: ic
integer, intent(in) :: nat, Z(nat)
real(8), intent(in) :: coord(3,nat)
logical, intent(out) :: changed
integer, allocatable :: bond(:,:)
integer :: nbond, b, p
logical :: same
real(8) :: ang
call build_bonds(nat, Z, coord, bond, nbond)
same = (nbond == ic%nbond)
if (same) then
   do b = 1, nbond
      if (bond(1,b) /= ic%bond(1,b) .or. bond(2,b) /= ic%bond(2,b)) then
         same = .false. ; exit
      endif
   enddo
endif

if (same) then
   do p = 1, ic%nprim
      if (ic%kind(p) == KIND_BEND) then
         ang = val_bend(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), coord(:,ic%at(3,p)))
         if (ang > LINEAR_TOL) then
            same = .false. ; exit
         endif
      else if (ic%kind(p) == KIND_LINBEND) then
         ang = val_bend(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), coord(:,ic%at(3,p)))
         if (ang < LINEAR_TOL - LINEAR_HYST) then
            same = .false. ; exit
         endif
      endif
   enddo
endif

changed = .not. same
if (changed) then
   block
      type(ic_set_t) :: newic
      call ic_build(newic, nat, Z, coord)
      ic = newic
   end block
else
   call ic_refresh_linbend_axes(ic, coord)
endif
end subroutine ic_rebuild_if_changed

subroutine ic_refresh_linbend_axes(ic, coord)
type(ic_set_t), intent(inout) :: ic
real(8), intent(in) :: coord(3,ic%nat)
integer :: p
real(8) :: ax(3), nax, u1(3), u2(3)
real(8), parameter :: PARALLEL_TOL = 0.5d0
do p = 1, ic%nprim
   if (ic%kind(p) /= KIND_LINBEND) cycle
   ax  = coord(:,ic%at(3,p)) - coord(:,ic%at(1,p))
   nax = sqrt(sum(ax*ax))
   if (nax < 1.0d-6) cycle
   ax = ax / nax
   if (abs(dot_product(ic%uref(:,p), ax)) < PARALLEL_TOL) cycle
   call perp_pair(ax, u1, u2)
   if (p > 1) then
      if (ic%kind(p-1) == KIND_LINBEND .and. &
          ic%at(1,p-1) == ic%at(1,p) .and. ic%at(3,p-1) == ic%at(3,p)) then
         ic%uref(:,p-1) = u1 ; ic%uref(:,p) = u2 ; cycle
      endif
   endif
   if (p < ic%nprim) then
      if (ic%kind(p+1) == KIND_LINBEND .and. &
          ic%at(1,p+1) == ic%at(1,p) .and. ic%at(3,p+1) == ic%at(3,p)) then
         ic%uref(:,p) = u1 ; ic%uref(:,p+1) = u2 ; cycle
      endif
   endif
   ic%uref(:,p) = u1
enddo
end subroutine ic_refresh_linbend_axes

subroutine ic_values(ic, coord, q)
type(ic_set_t), intent(in) :: ic
real(8), intent(in)  :: coord(3,ic%nat)
real(8), intent(out) :: q(ic%nprim)
integer :: p
do p = 1, ic%nprim
   select case (ic%kind(p))
   case (KIND_STRETCH)
      q(p) = val_stretch(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)))
   case (KIND_BEND)
      q(p) = val_bend(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), coord(:,ic%at(3,p)))
   case (KIND_TORSION, KIND_IMPROPER)
      q(p) = val_torsion(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), &
                         coord(:,ic%at(3,p)), coord(:,ic%at(4,p)))
   case (KIND_LINBEND)
      q(p) = val_linbend(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), &
                         coord(:,ic%at(3,p)), ic%uref(:,p))
   end select
enddo
end subroutine ic_values

subroutine ic_delta_q(ic, qa, qb, dq)
type(ic_set_t), intent(in) :: ic
real(8), intent(in)  :: qa(ic%nprim), qb(ic%nprim)
real(8), intent(out) :: dq(ic%nprim)
integer :: p
do p = 1, ic%nprim
   if (ic%kind(p) == KIND_TORSION .or. ic%kind(p) == KIND_IMPROPER) then
      dq(p) = wrap_pi(qa(p) - qb(p))
   else
      dq(p) = qa(p) - qb(p)
   endif
enddo
end subroutine ic_delta_q

subroutine ic_bmatrix(ic, coord, B)
type(ic_set_t), intent(in) :: ic
real(8), intent(in)  :: coord(3,ic%nat)
real(8), intent(out) :: B(ic%nprim, 3*ic%nat)
integer :: p, a
real(8) :: s1(3), s2(3), s3(3), s4(3)
logical :: ok
B = 0.0d0
do p = 1, ic%nprim
   select case (ic%kind(p))
   case (KIND_STRETCH)
      call bvec_stretch(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), s1, s2, ok)
      if (.not. ok) cycle
      call put(B, p, ic%at(1,p), s1) ; call put(B, p, ic%at(2,p), s2)
   case (KIND_BEND)
      call bvec_bend(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), coord(:,ic%at(3,p)), &
                     s1, s2, s3, ok)
      if (.not. ok) cycle
      call put(B, p, ic%at(1,p), s1) ; call put(B, p, ic%at(2,p), s2)
      call put(B, p, ic%at(3,p), s3)
   case (KIND_TORSION, KIND_IMPROPER)
      call bvec_torsion(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), &
                        coord(:,ic%at(3,p)), coord(:,ic%at(4,p)), s1, s2, s3, s4, ok)
      if (.not. ok) cycle
      call put(B, p, ic%at(1,p), s1) ; call put(B, p, ic%at(2,p), s2)
      call put(B, p, ic%at(3,p), s3) ; call put(B, p, ic%at(4,p), s4)
   case (KIND_LINBEND)
      call bvec_linbend(coord(:,ic%at(1,p)), coord(:,ic%at(2,p)), &
                        coord(:,ic%at(3,p)), ic%uref(:,p), s1, s2, s3, ok)
      if (.not. ok) cycle
      call put(B, p, ic%at(1,p), s1) ; call put(B, p, ic%at(2,p), s2)
      call put(B, p, ic%at(3,p), s3)
   end select
enddo
contains
   subroutine put(BB, ip, iat, sv)
   real(8), intent(inout) :: BB(ic%nprim, 3*ic%nat)
   integer, intent(in) :: ip, iat
   real(8), intent(in) :: sv(3)
   BB(ip, 3*(iat-1)+1) = BB(ip, 3*(iat-1)+1) + sv(1)
   BB(ip, 3*(iat-1)+2) = BB(ip, 3*(iat-1)+2) + sv(2)
   BB(ip, 3*(iat-1)+3) = BB(ip, 3*(iat-1)+3) + sv(3)
   end subroutine put
end subroutine ic_bmatrix

subroutine ic_write(ic, coord, unit_out, tag)
type(ic_set_t), intent(in) :: ic
real(8), intent(in) :: coord(3,ic%nat)
integer, intent(in) :: unit_out
character(len=*), intent(in) :: tag
real(8), allocatable :: q(:)
integer :: p
allocate(q(ic%nprim))
call ic_values(ic, coord, q)
write(unit_out,'(A)') '# Direwolf internal coordinate set: '//trim(tag)
write(unit_out,'(A,I0,A,I0)') '# natom = ', ic%nat, '   nprim = ', ic%nprim
write(unit_out,'(A)') '# geometry (bohr):'
do p = 1, ic%nat
   write(unit_out,'(A,I5,3F16.9)') '# ', p, coord(:,p)
enddo
write(unit_out,'(A)') '# kind  atoms                    value        [ref axis]'
do p = 1, ic%nprim
   select case (ic%kind(p))
   case (KIND_STRETCH)
      write(unit_out,'(A,2I5,20X,F16.9)') 'B ', ic%at(1,p), ic%at(2,p), q(p)
   case (KIND_BEND)
      write(unit_out,'(A,3I5,15X,F16.9)') 'A ', ic%at(1,p), ic%at(2,p), ic%at(3,p), q(p)
   case (KIND_LINBEND)
      write(unit_out,'(A,3I5,15X,F16.9,3F14.9)') 'L ', ic%at(1,p), ic%at(2,p), &
           ic%at(3,p), q(p), ic%uref(:,p)
   case (KIND_TORSION)
      write(unit_out,'(A,4I5,10X,F16.9)') 'D ', ic%at(1,p), ic%at(2,p), ic%at(3,p), &
           ic%at(4,p), q(p)
   case (KIND_IMPROPER)
      write(unit_out,'(A,4I5,10X,F16.9)') 'O ', ic%at(1,p), ic%at(2,p), ic%at(3,p), &
           ic%at(4,p), q(p)
   end select
enddo
end subroutine ic_write

integer function ic_rank(ic, coord)
type(ic_set_t), intent(in) :: ic
real(8), intent(in) :: coord(3,ic%nat)
real(8), allocatable :: B(:,:), G(:,:)
integer :: np, ndof
np = ic%nprim ; ndof = 3*ic%nat
ic_rank = 0
if (np < 1) return
allocate(B(np,ndof), G(np,np))
call ic_bmatrix(ic, coord, B)
call ic_ginv(B, np, ndof, G, ic_rank)
end function ic_rank

subroutine ic_ginv(B, np, ndof, Ginv, rank)
integer, intent(in)  :: np, ndof
real(8), intent(in)  :: B(np, ndof)
real(8), intent(out) :: Ginv(np, np)
integer, intent(out), optional :: rank
real(8), allocatable :: G(:,:), w(:), work(:), V(:,:)
integer :: i, j, lwork, info
real(8) :: tol
allocate(G(np,np), w(np), V(np,np))
call dgemm('N','T', np, np, ndof, 1.0d0, B, np, B, np, 0.0d0, G, np)
V = G
lwork = max(1, 3*np + 2*np*np)
allocate(work(lwork))
call dsyev('V', 'U', np, V, np, w, work, lwork, info)
Ginv = 0.0d0
if (present(rank)) rank = 0
if (info /= 0) return
tol = max(1.0d-10, 1.0d-6 * maxval(w))
if (present(rank)) rank = count(w > tol)
do i = 1, np
   if (w(i) > tol) then
      do j = 1, np
         Ginv(:,j) = Ginv(:,j) + (V(:,i) * V(j,i)) / w(i)
      enddo
   endif
enddo
end subroutine ic_ginv

subroutine ic_grad_to_internal(ic, coord, gx, gq)
type(ic_set_t), intent(in) :: ic
real(8), intent(in)  :: coord(3,ic%nat)
real(8), intent(in)  :: gx(3*ic%nat)
real(8), intent(out) :: gq(ic%nprim)
real(8), allocatable :: B(:,:), Ginv(:,:), bg(:)
integer :: np, ndof
np = ic%nprim ; ndof = 3*ic%nat
allocate(B(np,ndof), Ginv(np,np), bg(np))
call ic_bmatrix(ic, coord, B)
call ic_ginv(B, np, ndof, Ginv)
call dgemv('N', np, ndof, 1.0d0, B, np, gx, 1, 0.0d0, bg, 1)
call dgemv('N', np, np, 1.0d0, Ginv, np, bg, 1, 0.0d0, gq, 1)
end subroutine ic_grad_to_internal

subroutine ic_back_transform(ic, coord0, dq, dx, ok)
type(ic_set_t), intent(in) :: ic
real(8), intent(in)  :: coord0(3,ic%nat)
real(8), intent(in)  :: dq(ic%nprim)
real(8), intent(out) :: dx(3*ic%nat)
logical, intent(out) :: ok
real(8), allocatable :: B(:,:), Ginv(:,:), BtGi(:,:), q0(:), qc(:), dcur(:), resid(:), step(:), xit(:)
real(8) :: rmax, rfirst
integer :: np, ndof, it, i
integer, parameter :: MAXIT = 25

np = ic%nprim ; ndof = 3*ic%nat
allocate(B(np,ndof), Ginv(np,np), BtGi(ndof,np))
allocate(q0(np), qc(np), dcur(np), resid(np), step(ndof), xit(ndof))

call ic_bmatrix(ic, coord0, B)
call ic_ginv(B, np, ndof, Ginv)
call dgemm('T','N', ndof, np, np, 1.0d0, B, np, Ginv, np, 0.0d0, BtGi, ndof)

call ic_values(ic, coord0, q0)
xit = reshape(coord0, [ndof])
dx  = 0.0d0
ok  = .false.
rfirst = huge(1.0d0)

call dgemv('N', ndof, np, 1.0d0, BtGi, ndof, dq, 1, 0.0d0, dx, 1)

do it = 1, MAXIT
   call ic_values(ic, reshape(xit, [3,ic%nat]), qc)
   call ic_delta_q(ic, qc, q0, dcur)
   resid = dq - dcur
   rmax = maxval(abs(resid))
   if (it == 1) rfirst = rmax
   if (rmax < 1.0d-6) then
      ok = .true. ; exit
   endif
   if (it >= 5 .and. rmax > 2.0d0*rfirst) exit
   call dgemv('N', ndof, np, 1.0d0, BtGi, ndof, resid, 1, 0.0d0, step, 1)
   xit = xit + step
   if (maxval(abs(step)) < 1.0d-9) then
      ok = .true. ; exit
   endif
enddo

if (ok) then
   dx = xit - reshape(coord0, [ndof])
else if (rmax < 1.0d-4) then
   dx = xit - reshape(coord0, [ndof])
   ok = .true.
endif
end subroutine ic_back_transform

subroutine ic_diag_hessian(ic, Hq)
type(ic_set_t), intent(in) :: ic
real(8), intent(out) :: Hq(ic%nprim, ic%nprim)
integer :: p
Hq = 0.0d0
do p = 1, ic%nprim
   select case (ic%kind(p))
   case (KIND_STRETCH) ; Hq(p,p) = 0.50d0
   case (KIND_BEND)    ; Hq(p,p) = 0.20d0
   case (KIND_TORSION) ; Hq(p,p) = 0.10d0
   case (KIND_IMPROPER); Hq(p,p) = 0.10d0
   case (KIND_LINBEND) ; Hq(p,p) = 0.20d0
   end select
enddo
end subroutine ic_diag_hessian

subroutine ic_lindh_hessian(ic, coord, Hq)
type(ic_set_t), intent(in) :: ic
real(8), intent(in)  :: coord(3,ic%nat)
real(8), intent(out) :: Hq(ic%nprim, ic%nprim)
integer :: p, a, b, c, d, ra, rb, rc, rd
real(8) :: v
real(8), parameter :: HMIN_STRETCH = 1.0d-2
real(8), parameter :: HMIN_BEND    = 5.0d-3
real(8), parameter :: HMIN_TORSION = 1.0d-3

Hq = 0.0d0
do p = 1, ic%nprim
   a = ic%at(1,p) ; b = ic%at(2,p) ; c = ic%at(3,p) ; d = ic%at(4,p)
   select case (ic%kind(p))
   case (KIND_STRETCH)
      ra = lindh_row(ic%Z(a)) ; rb = lindh_row(ic%Z(b))
      v = K_STRETCH * lindh_rho(coord(:,a), coord(:,b), ra, rb)
      Hq(p,p) = max(v, HMIN_STRETCH)
   case (KIND_BEND, KIND_LINBEND)
      ra = lindh_row(ic%Z(a)) ; rb = lindh_row(ic%Z(b)) ; rc = lindh_row(ic%Z(c))
      v = K_BEND * lindh_rho(coord(:,a), coord(:,b), ra, rb) &
                 * lindh_rho(coord(:,b), coord(:,c), rb, rc)
      Hq(p,p) = max(v, HMIN_BEND)
   case (KIND_TORSION)
      ra = lindh_row(ic%Z(a)) ; rb = lindh_row(ic%Z(b))
      rc = lindh_row(ic%Z(c)) ; rd = lindh_row(ic%Z(d))
      v = K_TORSION * lindh_rho(coord(:,a), coord(:,b), ra, rb) &
                    * lindh_rho(coord(:,b), coord(:,c), rb, rc) &
                    * lindh_rho(coord(:,c), coord(:,d), rc, rd)
      Hq(p,p) = max(v, HMIN_TORSION)
   case (KIND_IMPROPER)
      ra = lindh_row(ic%Z(a)) ; rb = lindh_row(ic%Z(b))
      rc = lindh_row(ic%Z(c)) ; rd = lindh_row(ic%Z(d))
      v = K_BEND * lindh_rho(coord(:,c), coord(:,a), rc, ra) &
                 * lindh_rho(coord(:,c), coord(:,b), rc, rb) &
                 * lindh_rho(coord(:,c), coord(:,d), rc, rd)
      Hq(p,p) = max(v, HMIN_BEND)
   end select
enddo
end subroutine ic_lindh_hessian

end module mod_internal_coords
