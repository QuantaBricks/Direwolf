! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! GEDIIS (Li & Frisch, J. Chem. Theory Comput. 2, 835 (2006)) for the

module mod_opt_diis
implicit none
private
public :: gediis_solve

contains

subroutine gediis_solve(m, nq, E, dQ, G, cmax, c, qoff, gstar, Estar, ok)
integer, intent(in) :: m, nq
real(8), intent(in) :: E(m), dQ(nq,m), G(nq,m), cmax
real(8), intent(out) :: c(m), qoff(nq), gstar(nq), Estar
logical, intent(out) :: ok
real(8), allocatable :: A(:,:), K(:,:), rhs(:)
integer :: i, j, info
integer, allocatable :: ipiv(:)
ok = .false.
c = 0.0d0 ; qoff = 0.0d0 ; gstar = 0.0d0 ; Estar = 0.0d0
if (m < 2) return
allocate(A(m,m), K(m+1,m+1), rhs(m+1), ipiv(m+1))
do j = 1, m
   do i = 1, m
      A(i,j) = dot_product(G(:,i) - G(:,j), dQ(:,i) - dQ(:,j))
   enddo
enddo
K = 0.0d0
K(1:m,1:m) = -0.5d0 * A
K(1:m,m+1) = -1.0d0
K(m+1,1:m) = 1.0d0
rhs(1:m) = -E
rhs(m+1) = 1.0d0
call dgesv(m+1, 1, K, m+1, ipiv, rhs, m+1, info)
if (info /= 0) return
c = rhs(1:m)
if (maxval(abs(c)) > cmax) return
do i = 1, m
   qoff  = qoff  + c(i) * dQ(:,i)
   gstar = gstar + c(i) * G(:,i)
enddo
Estar = dot_product(c, E) - 0.25d0 * dot_product(c, matmul(A, c))
ok = .true.
end subroutine gediis_solve

end module mod_opt_diis
