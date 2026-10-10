! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Minimization pieces of the geometry optimizer: the RFO step, the BFGS

module mod_opt_min
use mod_internal_coords, only: ic_set_t, ic_bmatrix, ic_ginv
implicit none
private
public :: regularize_hessian, ric_floor_hessian, bfgs_update, solve_rfo

contains

subroutine regularize_hessian(H, n, out_unit)
integer, intent(in)    :: n, out_unit
real(8), intent(inout) :: H(n,n)
real(8), allocatable :: V(:,:), w(:), work(:), T(:,:)
integer :: i, j, lwork, info, nneg, nlift
real(8), parameter :: HMIN = 5.0d-3
allocate(V(n,n), w(n))
V = H
lwork = max(1, 3*n + 2*n*n)
allocate(work(lwork))
call dsyev('V', 'U', n, V, n, w, work, lwork, info)
if (info /= 0) then
   write(out_unit,'(A)') '[OPT] could not diagonalize the supplied Hessian - using it as given'
   return
endif
nneg  = count(w < -1.0d-6)
nlift = count(w < HMIN)
do i = 1, n
   w(i) = max(w(i), HMIN)
enddo
allocate(T(n,n))
do j = 1, n
   do i = 1, n
      T(i,j) = V(i,j) * w(j)
   enddo
enddo
call dgemm('N','T', n, n, n, 1.0d0, T, n, V, n, 0.0d0, H, n)
if (nneg > 0) then
   write(out_unit,'(A,I0,A,I0,A,ES9.2)') &
        '[OPT]  supplied Hessian had ', nneg, ' negative eigenvalue(s); ', nlift, &
        ' eigenvalue(s) raised to ', HMIN
else
   write(out_unit,'(A,I0,A,ES9.2)') '[OPT]  supplied Hessian: ', nlift, &
        ' eigenvalue(s) raised to ', HMIN
endif
end subroutine regularize_hessian

subroutine ric_floor_hessian(ic, coord3, nat, ni, Hq, out_unit)
type(ic_set_t), intent(in) :: ic
integer, intent(in)    :: nat, ni, out_unit
real(8), intent(in)    :: coord3(3,nat)
real(8), intent(inout) :: Hq(ni,ni)
real(8), parameter :: HQ_FLOOR = 1.0d-4
real(8), allocatable :: B(:,:), Ginv(:,:), G(:,:), P(:,:), T(:,:), V(:,:), w(:), work(:), Pv(:)
integer :: ndof, k, nneg, nlift, lwork, info
ndof = 3*nat
allocate(B(ni,ndof), Ginv(ni,ni), G(ni,ni), P(ni,ni), T(ni,ni), V(ni,ni), w(ni), Pv(ni))
call ic_bmatrix(ic, coord3, B)
call ic_ginv(B, ni, ndof, Ginv)
call dgemm('N','T', ni, ni, ndof, 1.0d0, B, ni, B, ni, 0.0d0, G, ni)
call dgemm('N','N', ni, ni, ni, 1.0d0, Ginv, ni, G, ni, 0.0d0, P, ni)
call dgemm('N','N', ni, ni, ni, 1.0d0, Hq, ni, P, ni, 0.0d0, T, ni)
call dgemm('T','N', ni, ni, ni, 1.0d0, P, ni, T, ni, 0.0d0, V, ni)
V = 0.5d0*(V + transpose(V))
lwork = max(1, 3*ni + ni*ni)
allocate(work(lwork))
call dsyev('V', 'U', ni, V, ni, w, work, lwork, info)
if (info /= 0) return
nneg = 0 ; nlift = 0
do k = 1, ni
   call dgemv('N', ni, ni, 1.0d0, P, ni, V(:,k), 1, 0.0d0, Pv, 1)
   if (dot_product(Pv, V(:,k)) < 0.5d0) then
      w(k) = 0.0d0
      cycle
   endif
   if (w(k) < -1.0d-6) nneg = nneg + 1
   if (abs(w(k)) < HQ_FLOOR .or. w(k) < 0.0d0) then
      nlift = nlift + 1
      w(k) = max(abs(w(k)), HQ_FLOOR)
   endif
enddo
do k = 1, ni
   T(:,k) = V(:,k)*w(k)
enddo
call dgemm('N','T', ni, ni, ni, 1.0d0, T, ni, V, ni, 0.0d0, Hq, ni)
if (nlift > 0) write(out_unit,'(A,I0,A,I0,A,ES9.2)') '  [opt] internal Hessian: ', nneg, &
     ' negative eigenvalue(s) flipped, ', nlift, ' raised to at least ', HQ_FLOOR
end subroutine ric_floor_hessian

subroutine bfgs_update(n, H, s, y)
integer, intent(in)    :: n
real(8), intent(inout) :: H(n,n)
real(8), intent(in)    :: s(n), y(n)
real(8), allocatable   :: Hs(:)
real(8) :: sy, sHs
integer :: i, j
sy = dot_product(s, y)
if (sy <= 1.0d-10) return
allocate(Hs(n))
call dgemv('N', n, n, 1.0d0, H, n, s, 1, 0.0d0, Hs, 1)
sHs = dot_product(s, Hs)
if (sHs <= 1.0d-12) return
do j = 1, n
   do i = 1, n
      H(i,j) = H(i,j) - Hs(i)*Hs(j)/sHs + y(i)*y(j)/sy
   enddo
enddo
end subroutine bfgs_update

subroutine solve_rfo(n, H, g, step, out_unit)
integer, intent(in)  :: n, out_unit
real(8), intent(in)  :: H(n,n), g(n)
real(8), intent(out) :: step(n)
real(8), allocatable :: aug(:,:), w(:), work(:)
integer :: lwork, info
real(8) :: vlast
allocate(aug(n+1,n+1), w(n+1))
aug = 0.0d0
aug(1:n,1:n) = H
aug(1:n,n+1) = g
aug(n+1,1:n) = g
lwork = max(1, 3*(n+1) + (n+1)*(n+1))
allocate(work(lwork))
call dsyev('V', 'U', n+1, aug, n+1, w, work, lwork, info)
if (info /= 0) then
   write(out_unit,'(A,I0)') '  [rfo] dsyev info=', info
   step = -g
   return
endif
vlast = aug(n+1,1)
if (abs(vlast) > 1.0d-8) then
   step = aug(1:n,1) / vlast
else
   step = -g
endif
end subroutine solve_rfo

end module mod_opt_min
