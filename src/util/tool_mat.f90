! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Small dense-matrix utility routines (e.g. upper-triangular-to-full expansion).

subroutine  UPtri2Full(A,n,B)
     implicit none
     integer :: n
     real(8) :: A((n-1)*n/2+1)
     integer :: i,j
     real(8) :: B(n,n)

do i = 1,n
  do j = 1,n
     if (j .ge. i) then
        B(i,j)= A(j*(j-1)/2+i)
     else
        B(i,j)= A(i*(i-1)/2+j)
     endif
  enddo
enddo

end

integer function factorial(n)
integer i,n

factorial = 1
do i = 1,n
   factorial = factorial*i
enddo

end function

integer function factorial2(n)
integer i,n

factorial2 = 1
do i = 1,n
   factorial2 = factorial2*(2*i-1)
enddo

end function

