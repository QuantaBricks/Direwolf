! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Lindh model Hessian - a cheap analytic Cartesian Hessian built from

module mod_lindh_hessian
use mod_wilson_bvec, only: bvec_stretch, bvec_bend, bvec_torsion
implicit none
private
public :: lindh_cart_hessian
public :: lindh_rho, lindh_row
public :: K_STRETCH, K_BEND, K_TORSION

real(8), parameter :: L_ALPHA(3,3) = reshape([ &
   1.0000d0, 0.3949d0, 0.3949d0, &
   0.3949d0, 0.2800d0, 0.2800d0, &
   0.3949d0, 0.2800d0, 0.2800d0 ], [3,3])
real(8), parameter :: L_RREF(3,3) = reshape([ &
   1.35d0, 2.10d0, 2.53d0, &
   2.10d0, 2.87d0, 3.40d0, &
   2.53d0, 3.40d0, 3.40d0 ], [3,3])

real(8), parameter :: K_STRETCH = 0.45d0
real(8), parameter :: K_BEND    = 0.15d0
real(8), parameter :: K_TORSION = 0.005d0

real(8), parameter :: RHO_TAU = 1.0d-6

contains

subroutine lindh_cart_hessian(natom, Z, coord, H)
integer, intent(in)  :: natom
integer, intent(in)  :: Z(natom)
real(8), intent(in)  :: coord(3,natom)
real(8), intent(out) :: H(3*natom, 3*natom)

integer :: row(natom)
integer :: ia, ib, ic, id
real(8) :: rho_ab, rho_bc, rho_cd, coef
real(8) :: sA(3), sB(3), sC(3), sD(3)
logical :: ok

do ia = 1, natom
   row(ia) = lindh_row(Z(ia))
enddo

H = 0.0d0

do ia = 1, natom
   do ib = ia+1, natom
      rho_ab = lindh_rho(coord(:,ia), coord(:,ib), row(ia), row(ib))
      if (rho_ab < RHO_TAU) cycle
      call bvec_stretch(coord(:,ia), coord(:,ib), sA, sB, ok)
      if (.not. ok) cycle
      coef = K_STRETCH * rho_ab
      call accum_pair(H, natom, ia, sA, ib, sB, coef)
   enddo
enddo

do ib = 1, natom
   do ia = 1, natom
      if (ia == ib) cycle
      rho_ab = lindh_rho(coord(:,ia), coord(:,ib), row(ia), row(ib))
      if (rho_ab < RHO_TAU) cycle
      do ic = ia+1, natom
         if (ic == ib) cycle
         rho_bc = lindh_rho(coord(:,ib), coord(:,ic), row(ib), row(ic))
         if (rho_ab * rho_bc < RHO_TAU) cycle
         call bvec_bend(coord(:,ia), coord(:,ib), coord(:,ic), sA, sB, sC, ok)
         if (.not. ok) cycle
         coef = K_BEND * rho_ab * rho_bc
         call accum_triple(H, natom, ia, sA, ib, sB, ic, sC, coef)
      enddo
   enddo
enddo

do ib = 1, natom
   do ic = ib+1, natom
      rho_bc = lindh_rho(coord(:,ib), coord(:,ic), row(ib), row(ic))
      if (rho_bc < RHO_TAU) cycle
      do ia = 1, natom
         if (ia == ib .or. ia == ic) cycle
         rho_ab = lindh_rho(coord(:,ia), coord(:,ib), row(ia), row(ib))
         if (rho_ab * rho_bc < RHO_TAU) cycle
         do id = 1, natom
            if (id == ia .or. id == ib .or. id == ic) cycle
            rho_cd = lindh_rho(coord(:,ic), coord(:,id), row(ic), row(id))
            if (rho_ab * rho_bc * rho_cd < RHO_TAU) cycle
            call bvec_torsion(coord(:,ia), coord(:,ib), coord(:,ic), coord(:,id), &
                              sA, sB, sC, sD, ok)
            if (.not. ok) cycle
            coef = K_TORSION * rho_ab * rho_bc * rho_cd
            call accum_quad(H, natom, ia, sA, ib, sB, ic, sC, id, sD, coef)
         enddo
      enddo
   enddo
enddo

end subroutine lindh_cart_hessian

integer function lindh_row(Zi)
integer, intent(in) :: Zi
if (Zi <= 2) then
   lindh_row = 1
else if (Zi <= 10) then
   lindh_row = 2
else
   lindh_row = 3
endif
end function lindh_row

real(8) function lindh_rho(ra, rb, ri, rj)
real(8), intent(in) :: ra(3), rb(3)
integer, intent(in) :: ri, rj
real(8) :: r2, a, rref
r2   = sum((ra - rb)**2)
a    = L_ALPHA(ri, rj)
rref = L_RREF(ri, rj)
lindh_rho = exp(a * (rref*rref - r2))
end function lindh_rho

subroutine accum_pair(H, natom, i1, s1, i2, s2, coef)
integer, intent(in)    :: natom, i1, i2
real(8), intent(in)    :: s1(3), s2(3), coef
real(8), intent(inout) :: H(3*natom, 3*natom)
integer :: at(2), a, b, p, q
real(8) :: sv(3,2)
at = [i1, i2] ; sv(:,1) = s1 ; sv(:,2) = s2
do a = 1, 2
   do b = 1, 2
      do p = 1, 3
         do q = 1, 3
            H(3*(at(a)-1)+p, 3*(at(b)-1)+q) = H(3*(at(a)-1)+p, 3*(at(b)-1)+q) &
                 + coef * sv(p,a) * sv(q,b)
         enddo
      enddo
   enddo
enddo
end subroutine accum_pair

subroutine accum_triple(H, natom, i1, s1, i2, s2, i3, s3, coef)
integer, intent(in)    :: natom, i1, i2, i3
real(8), intent(in)    :: s1(3), s2(3), s3(3), coef
real(8), intent(inout) :: H(3*natom, 3*natom)
integer :: at(3), a, b, p, q
real(8) :: sv(3,3)
at = [i1, i2, i3] ; sv(:,1) = s1 ; sv(:,2) = s2 ; sv(:,3) = s3
do a = 1, 3
   do b = 1, 3
      do p = 1, 3
         do q = 1, 3
            H(3*(at(a)-1)+p, 3*(at(b)-1)+q) = H(3*(at(a)-1)+p, 3*(at(b)-1)+q) &
                 + coef * sv(p,a) * sv(q,b)
         enddo
      enddo
   enddo
enddo
end subroutine accum_triple

subroutine accum_quad(H, natom, i1, s1, i2, s2, i3, s3, i4, s4, coef)
integer, intent(in)    :: natom, i1, i2, i3, i4
real(8), intent(in)    :: s1(3), s2(3), s3(3), s4(3), coef
real(8), intent(inout) :: H(3*natom, 3*natom)
integer :: at(4), a, b, p, q
real(8) :: sv(3,4)
at = [i1, i2, i3, i4]
sv(:,1) = s1 ; sv(:,2) = s2 ; sv(:,3) = s3 ; sv(:,4) = s4
do a = 1, 4
   do b = 1, 4
      do p = 1, 3
         do q = 1, 3
            H(3*(at(a)-1)+p, 3*(at(b)-1)+q) = H(3*(at(a)-1)+p, 3*(at(b)-1)+q) &
                 + coef * sv(p,a) * sv(q,b)
         enddo
      enddo
   enddo
enddo
end subroutine accum_quad

end module mod_lindh_hessian
