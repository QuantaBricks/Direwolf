! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! 128-bit fixed-point accumulation (value * 2^64 held in an integer(16))

module mod_fixq
implicit none
private
public :: to_q, from_q

real(8), parameter :: S32 = 4294967296.0d0

contains

elemental integer(16) function to_q(x)
real(8), intent(in) :: x
real(8) :: t
integer(8) :: hi, lo
t = x*S32
hi = int(t, 8)
lo = int((t - real(hi,8))*S32, 8)
to_q = ishft(int(hi,16), 32) + int(lo,16)
end function to_q

elemental real(8) function from_q(q)
integer(16), intent(in) :: q
integer(8) :: hi, lo
hi = int(shifta(q, 32), 8)
lo = int(iand(q, int(z'FFFFFFFF',16)), 8)
from_q = (real(hi,8) + real(lo,8)/S32)/S32
end function from_q

end module mod_fixq
