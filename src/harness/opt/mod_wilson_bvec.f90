! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Wilson B-vectors and values for the four primitive internal

module mod_wilson_bvec
implicit none
private
public :: bvec_stretch, bvec_bend, bvec_torsion
public :: val_stretch, val_bend, val_torsion
public :: cross3, wrap_pi

contains

pure function cross3(a, b) result(c)
real(8), intent(in) :: a(3), b(3)
real(8) :: c(3)
c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = a(3)*b(1) - a(1)*b(3)
c(3) = a(1)*b(2) - a(2)*b(1)
end function cross3

elemental function wrap_pi(x) result(y)
real(8), intent(in) :: x
real(8) :: y
real(8), parameter :: PI = 3.14159265358979323846d0
y = x - 2.0d0*PI * nint(x / (2.0d0*PI))
end function wrap_pi

real(8) function val_stretch(ra, rb)
real(8), intent(in) :: ra(3), rb(3)
val_stretch = sqrt(sum((ra - rb)**2))
end function val_stretch

real(8) function val_bend(ra, rb, rc)
real(8), intent(in) :: ra(3), rb(3), rc(3)
real(8) :: u(3), v(3), cph
u = ra - rb ; u = u / max(sqrt(sum(u*u)), 1.0d-300)
v = rc - rb ; v = v / max(sqrt(sum(v*v)), 1.0d-300)
cph = max(-1.0d0, min(1.0d0, dot_product(u, v)))
val_bend = acos(cph)
end function val_bend

real(8) function val_torsion(ra, rb, rc, rd)
real(8), intent(in) :: ra(3), rb(3), rc(3), rd(3)
real(8) :: b1(3), b2(3), b3(3), n1(3), n2(3), m1(3), x, y
b1 = rb - ra
b2 = rc - rb
b3 = rd - rc
n1 = cross3(b1, b2)
n2 = cross3(b2, b3)
m1 = cross3(n1, b2 / max(sqrt(sum(b2*b2)), 1.0d-300))
x = dot_product(n1, n2)
y = dot_product(m1, n2)
val_torsion = atan2(y, x)
end function val_torsion

subroutine bvec_stretch(ra, rb, sa, sb, ok)
real(8), intent(in)  :: ra(3), rb(3)
real(8), intent(out) :: sa(3), sb(3)
logical, intent(out) :: ok
real(8) :: d(3), r
d = ra - rb
r = sqrt(sum(d*d))
if (r < 1.0d-6) then
   ok = .false.; sa = 0.0d0; sb = 0.0d0; return
endif
sa = d / r
sb = -sa
ok = .true.
end subroutine bvec_stretch

subroutine bvec_bend(ra, rb, rc, sa, sb, sc, ok)
real(8), intent(in)  :: ra(3), rb(3), rc(3)
real(8), intent(out) :: sa(3), sb(3), sc(3)
logical, intent(out) :: ok
real(8) :: u(3), v(3), ru, rv, cph, sph
u = ra - rb ; ru = sqrt(sum(u*u))
v = rc - rb ; rv = sqrt(sum(v*v))
if (ru < 1.0d-6 .or. rv < 1.0d-6) then
   ok = .false.; sa = 0.0d0; sb = 0.0d0; sc = 0.0d0; return
endif
u = u / ru ; v = v / rv
cph = dot_product(u, v)
sph = sqrt(max(1.0d0 - cph*cph, 0.0d0))
if (sph < 1.0d-4) then
   ok = .false.; sa = 0.0d0; sb = 0.0d0; sc = 0.0d0; return
endif
sa = (cph*u - v) / (ru*sph)
sc = (cph*v - u) / (rv*sph)
sb = -(sa + sc)
ok = .true.
end subroutine bvec_bend

subroutine bvec_torsion(ra, rb, rc, rd, sa, sb, sc, sd, ok)
real(8), intent(in)  :: ra(3), rb(3), rc(3), rd(3)
real(8), intent(out) :: sa(3), sb(3), sc(3), sd(3)
logical, intent(out) :: ok
real(8) :: b1(3), b2(3), b3(3), n1(3), n2(3), nb2, n1sq, n2sq, c1, c2
b1 = rb - ra
b2 = rc - rb
b3 = rd - rc
n1 = cross3(b1, b2)
n2 = cross3(b2, b3)
nb2  = sqrt(sum(b2*b2))
n1sq = sum(n1*n1)
n2sq = sum(n2*n2)
if (nb2 < 1.0d-6 .or. n1sq < 1.0d-8 .or. n2sq < 1.0d-8) then
   ok = .false.; sa = 0.0d0; sb = 0.0d0; sc = 0.0d0; sd = 0.0d0; return
endif
sa =  (nb2 / n1sq) * n1
sd = -(nb2 / n2sq) * n2
c1 = dot_product(b1, b2) / (nb2*nb2)
c2 = dot_product(b3, b2) / (nb2*nb2)
sb = -(1.0d0 + c1)*sa + c2*sd
sc = -(1.0d0 + c2)*sd + c1*sa
ok = .true.
end subroutine bvec_torsion

end module mod_wilson_bvec
