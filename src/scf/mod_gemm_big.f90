! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

module mod_gemm_big
implicit none
private
public :: gemm_big, level_shift_big
contains
subroutine gemm_big(ta, tb, n, A, B, C)
character(len=1),intent(in) :: ta, tb
integer,intent(in) :: n
real(8),intent(in) :: A(n,n), B(n,n)
real(8),intent(out) :: C(n,n)
logical :: done
done = .false.
if (.not. done) call dgemm(ta, tb, n, n, n, 1.0d0, A, n, B, n, 0.0d0, C, n)
end subroutine gemm_big

function level_shift_big(n, nvirt, shift, S, Cv, F, W1, W2) result(done)
integer,intent(in) :: n, nvirt
real(8),intent(in) :: shift, S(n,n), Cv(n,nvirt)
real(8),intent(inout) :: F(n,n)
real(8),intent(out) :: W1(n,n), W2(n,n)
logical :: done
done = .false.
end function level_shift_big
end module mod_gemm_big
