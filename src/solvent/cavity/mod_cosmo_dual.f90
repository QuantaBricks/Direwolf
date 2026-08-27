! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_cosmo_dual: forward-mode dual-number scaffolding (dscal/dvec3)

module mod_cosmo_dual
use mod_cosmo_constants, only: COSMO_MAXSLOT
implicit none

type dscal
   real(8) :: v
   real(8) :: g(3,COSMO_MAXSLOT)
end type dscal

type dvec3
   real(8) :: v(3)
   real(8) :: g(3,3,COSMO_MAXSLOT)
end type dvec3

interface operator(+)
   module procedure dscal_add, dvec3_add
end interface
interface operator(-)
   module procedure dscal_sub, dvec3_sub, dscal_neg, dvec3_neg
end interface
interface operator(*)
   module procedure dscal_mul, real_dscal_mul, dscal_real_mul, real_dvec3_mul, dvec3_real_mul, &
                     dvec3_dscal_mul, dscal_dvec3_mul
end interface
interface operator(/)
   module procedure dscal_div, dscal_real_div, real_dscal_div, dvec3_dscal_div
end interface

contains

subroutine cosmo_cross(a, b, c)
real(8),intent(in) :: a(3), b(3)
real(8),intent(out) :: c(3)
c(1) = a(2)*b(3) - a(3)*b(2)
c(2) = a(3)*b(1) - a(1)*b(3)
c(3) = a(1)*b(2) - a(2)*b(1)
end subroutine cosmo_cross

function const_dscal(x) result(r)
real(8),intent(in) :: x
type(dscal) :: r
r%v = x; r%g = 0.0d0
end function const_dscal

function const_dvec3(x) result(r)
real(8),intent(in) :: x(3)
type(dvec3) :: r
r%v = x; r%g = 0.0d0
end function const_dvec3

function dvec3_sparse(x, slotA, JA, slotB, JB) result(r)
real(8),intent(in) :: x(3)
integer,intent(in) :: slotA, slotB
real(8),intent(in) :: JA(3,3), JB(3,3)
type(dvec3) :: r
r%v = x; r%g = 0.0d0
if (slotA > 0) r%g(:,:,slotA) = r%g(:,:,slotA) + JA
if (slotB > 0) r%g(:,:,slotB) = r%g(:,:,slotB) + JB
end function dvec3_sparse

function dscal_add(a,b) result(r)
type(dscal),intent(in) :: a,b
type(dscal) :: r
r%v = a%v+b%v; r%g = a%g+b%g
end function dscal_add
function dvec3_add(a,b) result(r)
type(dvec3),intent(in) :: a,b
type(dvec3) :: r
r%v = a%v+b%v; r%g = a%g+b%g
end function dvec3_add
function dscal_sub(a,b) result(r)
type(dscal),intent(in) :: a,b
type(dscal) :: r
r%v = a%v-b%v; r%g = a%g-b%g
end function dscal_sub
function dvec3_sub(a,b) result(r)
type(dvec3),intent(in) :: a,b
type(dvec3) :: r
r%v = a%v-b%v; r%g = a%g-b%g
end function dvec3_sub
function dscal_neg(a) result(r)
type(dscal),intent(in) :: a
type(dscal) :: r
r%v = -a%v; r%g = -a%g
end function dscal_neg
function dvec3_neg(a) result(r)
type(dvec3),intent(in) :: a
type(dvec3) :: r
r%v = -a%v; r%g = -a%g
end function dvec3_neg
function dscal_mul(a,b) result(r)
type(dscal),intent(in) :: a,b
type(dscal) :: r
r%v = a%v*b%v; r%g = a%g*b%v + a%v*b%g
end function dscal_mul
function real_dscal_mul(x,b) result(r)
real(8),intent(in) :: x
type(dscal),intent(in) :: b
type(dscal) :: r
r%v = x*b%v; r%g = x*b%g
end function real_dscal_mul
function dscal_real_mul(a,x) result(r)
type(dscal),intent(in) :: a
real(8),intent(in) :: x
type(dscal) :: r
r%v = a%v*x; r%g = a%g*x
end function dscal_real_mul
function real_dvec3_mul(x,b) result(r)
real(8),intent(in) :: x
type(dvec3),intent(in) :: b
type(dvec3) :: r
r%v = x*b%v; r%g = x*b%g
end function real_dvec3_mul
function dvec3_real_mul(a,x) result(r)
type(dvec3),intent(in) :: a
real(8),intent(in) :: x
type(dvec3) :: r
r%v = a%v*x; r%g = a%g*x
end function dvec3_real_mul
function dvec3_dscal_mul(a,b) result(r)
type(dvec3),intent(in) :: a
type(dscal),intent(in) :: b
type(dvec3) :: r
integer :: i
r%v = a%v*b%v
do i = 1,3
   r%g(i,:,:) = a%g(i,:,:)*b%v + a%v(i)*b%g(:,:)
enddo
end function dvec3_dscal_mul
function dscal_dvec3_mul(a,b) result(r)
type(dscal),intent(in) :: a
type(dvec3),intent(in) :: b
type(dvec3) :: r
r = dvec3_dscal_mul(b,a)
end function dscal_dvec3_mul
function dscal_div(a,b) result(r)
type(dscal),intent(in) :: a,b
type(dscal) :: r
r%v = a%v/b%v; r%g = (a%g*b%v - a%v*b%g)/(b%v*b%v)
end function dscal_div
function dscal_real_div(a,x) result(r)
type(dscal),intent(in) :: a
real(8),intent(in) :: x
type(dscal) :: r
r%v = a%v/x; r%g = a%g/x
end function dscal_real_div
function real_dscal_div(x,b) result(r)
real(8),intent(in) :: x
type(dscal),intent(in) :: b
type(dscal) :: r
r%v = x/b%v; r%g = -x*b%g/(b%v*b%v)
end function real_dscal_div

function cosmo_idn3() result(I3)
real(8) :: I3(3,3)
integer :: a
I3 = 0.0d0
do a = 1,3
   I3(a,a) = 1.0d0
enddo
end function cosmo_idn3

function cosmo_get_slot(atom_of_slot, nslots, atom_idx) result(islot)
integer,intent(inout) :: atom_of_slot(COSMO_MAXSLOT)
integer,intent(inout) :: nslots
integer,intent(in) :: atom_idx
integer :: islot, s
do s = 1,nslots
   if (atom_of_slot(s) == atom_idx) then
      islot = s
      return
   endif
enddo
nslots = nslots+1
atom_of_slot(nslots) = atom_idx
islot = nslots
end function cosmo_get_slot
function dvec3_dscal_div(a,b) result(r)
type(dvec3),intent(in) :: a
type(dscal),intent(in) :: b
type(dvec3) :: r
integer :: i
r%v = a%v/b%v
do i = 1,3
   r%g(i,:,:) = (a%g(i,:,:)*b%v - a%v(i)*b%g(:,:))/(b%v*b%v)
enddo
end function dvec3_dscal_div

function ddot(a,b) result(r)
type(dvec3),intent(in) :: a,b
type(dscal) :: r
integer :: i
r%v = sum(a%v*b%v)
r%g = 0.0d0
do i = 1,3
   r%g = r%g + a%g(i,:,:)*b%v(i) + a%v(i)*b%g(i,:,:)
enddo
end function ddot

function dcross(a,b) result(r)
type(dvec3),intent(in) :: a,b
type(dvec3) :: r
integer :: s,c
r%v(1) = a%v(2)*b%v(3)-a%v(3)*b%v(2)
r%v(2) = a%v(3)*b%v(1)-a%v(1)*b%v(3)
r%v(3) = a%v(1)*b%v(2)-a%v(2)*b%v(1)
do s = 1,COSMO_MAXSLOT
   do c = 1,3
      r%g(1,c,s) = a%g(2,c,s)*b%v(3)-a%g(3,c,s)*b%v(2) + a%v(2)*b%g(3,c,s)-a%v(3)*b%g(2,c,s)
      r%g(2,c,s) = a%g(3,c,s)*b%v(1)-a%g(1,c,s)*b%v(3) + a%v(3)*b%g(1,c,s)-a%v(1)*b%g(3,c,s)
      r%g(3,c,s) = a%g(1,c,s)*b%v(2)-a%g(2,c,s)*b%v(1) + a%v(1)*b%g(2,c,s)-a%v(2)*b%g(1,c,s)
   enddo
enddo
end function dcross

function dsqrt_d(a) result(r)
type(dscal),intent(in) :: a
type(dscal) :: r
r%v = sqrt(a%v)
r%g = a%g/(2.0d0*r%v)
end function dsqrt_d

function dacos_d(a) result(r)
type(dscal),intent(in) :: a
type(dscal) :: r
real(8) :: x, s
x = a%v
if (x > 1.0d0) x = 1.0d0
if (x < -1.0d0) x = -1.0d0
r%v = acos(x)
s = 1.0d0 - x*x
if (s < 1.0d-14) s = 1.0d-14
r%g = -a%g/sqrt(s)
end function dacos_d

end module mod_cosmo_dual
