! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Matrix pretty-printing utilities for log output.

subroutine prtLowMat(A,n,C)

integer :: n
real(8) :: A(n*(n+1)/2+1)
character(20) :: C

write(*,"(5X,A30)") C
do i = 1,n
  do j = 1,n
     if (j .ge. i) then
        write(*,"(F10.6)",advance='no') A(j*(j-1)/2+i)
     else
        write(*,"(F10.6)",advance='no') A(i*(i-1)/2+j)
     endif
  enddo
  write(*,*)
enddo
end subroutine

subroutine prtMat(A,n,C)

integer :: n
real(8) :: A(n,n)
character(20) :: C

write(*,"(5X,A30)") C
do i = 1,n
  do j = 1,n
        write(*,"(F11.7)",advance='no') A(i,j)
  enddo
  write(*,*)
enddo
end subroutine

