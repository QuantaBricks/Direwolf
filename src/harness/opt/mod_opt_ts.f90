! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Transition-state search pieces of the geometry optimizer (opt_ts):

module mod_opt_ts
implicit none
private
public :: ts_on, ts_nneg, ts_vprev, TS_MAX_RECALC
public :: ts_start, ts_forget_mode, ts_endgame_refresh_due, ts_lost_mode_recalc_due
public :: solve_prfo, bofill_update

logical, save :: ts_on = .false.
real(8), allocatable, save :: ts_vprev(:)
integer, save :: ts_nneg = 1

integer, parameter :: TS_MAX_RECALC = 5
real(8), parameter :: TS_REFRESH_GMAX = 2.0d-3

contains

subroutine ts_start(on)
logical, intent(in) :: on
ts_on = on
if (allocated(ts_vprev)) deallocate(ts_vprev)
ts_nneg = 1
end subroutine ts_start

subroutine ts_forget_mode()
if (allocated(ts_vprev)) deallocate(ts_vprev)
end subroutine ts_forget_mode

logical function ts_endgame_refresh_due(cyc, gmax, nrecalc, multi_done)
integer, intent(in) :: cyc, nrecalc
real(8), intent(in) :: gmax
logical, intent(in) :: multi_done
ts_endgame_refresh_due = ts_on .and. cyc > 1 .and. ts_nneg > 1 .and. gmax < TS_REFRESH_GMAX &
                         .and. .not. multi_done .and. nrecalc < TS_MAX_RECALC
end function ts_endgame_refresh_due

logical function ts_lost_mode_recalc_due(hess_fresh, nrecalc)
logical, intent(in) :: hess_fresh
integer, intent(in) :: nrecalc
ts_lost_mode_recalc_due = ts_on .and. ts_nneg == 0 .and. .not. hess_fresh .and. nrecalc < TS_MAX_RECALC
end function ts_lost_mode_recalc_due

subroutine solve_prfo(n, H, g, step, out_unit)
integer, intent(in)  :: n, out_unit
real(8), intent(in)  :: H(n,n), g(n)
real(8), intent(out) :: step(n)
real(8), allocatable :: V(:,:), b(:), F(:), s(:), work(:), aug(:,:), w(:)
integer :: lwork, info, k, kf, m, nneg
real(8) :: ov, ovbest, lam_p, vlast
allocate(V(n,n), b(n), F(n), s(n))
V = H
lwork = max(1, 3*n + n*n)
allocate(work(lwork))
call dsyev('V', 'U', n, V, n, b, work, lwork, info)
if (info /= 0) then
   write(out_unit,'(A,I0)') '  [ts] dsyev info=', info
   step = -g
   return
endif
call dgemv('T', n, n, 1.0d0, V, n, g, 1, 0.0d0, F, 1)

kf = 1
ovbest = 1.0d0
if (allocated(ts_vprev)) then
   if (size(ts_vprev) == n) then
      ovbest = -1.0d0
      do k = 1, n
         ov = abs(dot_product(V(:,k), ts_vprev))
         if (ov > ovbest) then
            ovbest = ov
            kf = k
         endif
      enddo
      if (ovbest < 0.5d0) then
         write(out_unit,'(A,F6.3,A)') '  [ts] followed mode lost (overlap ', ovbest, ') - using the lowest'
         kf = 1
      endif
   endif
endif
nneg = count(b < -1.0d-6)
ts_nneg = nneg
write(out_unit,'(A,I0,A,ES10.3,A,ES10.3,A,I0,A,F6.3)') '  [ts] following mode ', kf, &
     '  eigenvalue ', b(kf), '  lowest ', b(1), '  negative ', nneg, '  overlap ', ovbest
if (nneg == 0) write(out_unit,'(A)') '  [ts] WARNING: the Hessian has no negative eigenvalue'

s = 0.0d0
if (abs(b(kf)) + abs(F(kf)) > 1.0d-14) then
   lam_p = 0.5d0*b(kf) + 0.5d0*sqrt(b(kf)**2 + 4.0d0*F(kf)**2)
   if (abs(b(kf) - lam_p) > 1.0d-12) s(kf) = -F(kf) / (b(kf) - lam_p)
endif
m = n - 1
if (m >= 1) then
   allocate(aug(m+1,m+1), w(m+1))
   aug = 0.0d0
   m = 0
   do k = 1, n
      if (k == kf) cycle
      m = m + 1
      aug(m,m) = b(k)
      aug(m,n) = F(k)
      aug(n,m) = F(k)
   enddo
   deallocate(work)
   lwork = max(1, 3*n + n*n)
   allocate(work(lwork))
   call dsyev('V', 'U', n, aug, n, w, work, lwork, info)
   if (info == 0) then
      vlast = aug(n,1)
      m = 0
      do k = 1, n
         if (k == kf) cycle
         m = m + 1
         if (abs(vlast) > 1.0d-8) then
            s(k) = aug(m,1) / vlast
         else
            s(k) = -F(k) / max(b(k), 1.0d-3)
         endif
      enddo
   else
      write(out_unit,'(A,I0)') '  [ts] dsyev (minimisation block) info=', info
      m = 0
      do k = 1, n
         if (k == kf) cycle
         s(k) = -F(k) / max(b(k), 1.0d-3)
      enddo
   endif
endif
call dgemv('N', n, n, 1.0d0, V, n, s, 1, 0.0d0, step, 1)

if (allocated(ts_vprev)) then
   if (size(ts_vprev) /= n) deallocate(ts_vprev)
endif
if (.not. allocated(ts_vprev)) allocate(ts_vprev(n))
ts_vprev = V(:,kf)
end subroutine solve_prfo

subroutine bofill_update(n, H, s, y)
integer, intent(in)    :: n
real(8), intent(inout) :: H(n,n)
real(8), intent(in)    :: s(n), y(n)
real(8) :: r(n), ss, rr, rs, phi
integer :: i, j
call dgemv('N', n, n, 1.0d0, H, n, s, 1, 0.0d0, r, 1)
r = y - r
ss = dot_product(s, s)
rr = dot_product(r, r)
rs = dot_product(r, s)
if (ss < 1.0d-20 .or. rr < 1.0d-20) return
phi = rs*rs / (rr*ss)
do j = 1, n
   do i = 1, n
      H(i,j) = H(i,j) + (1.0d0 - phi) * ((r(i)*s(j) + s(i)*r(j))/ss - rs*s(i)*s(j)/(ss*ss))
      if (abs(rs) > 1.0d-14) H(i,j) = H(i,j) + phi * r(i)*r(j)/rs
   enddo
enddo
end subroutine bofill_update

end module mod_opt_ts
