! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Transition-state search pieces of the geometry optimizer (opt_ts):

module mod_opt_ts
use mod_wilson_bvec, only: bvec_stretch, bvec_bend, bvec_torsion
implicit none
private
public :: ts_on, ts_nneg, ts_vprev, TS_MAX_RECALC
public :: ts_start, ts_forget_mode, ts_endgame_refresh_due, ts_lost_mode_recalc_due
public :: solve_prfo, bofill_update
public :: ts_follow_init, ts_follow_context, ts_abort

logical, save :: ts_on = .false.
real(8), allocatable, save :: ts_vprev(:)
integer, save :: ts_nneg = 1
logical, save :: ts_fresh = .true.
logical, save :: ts_abort = .false.

integer, parameter :: MAXFOL = 16
integer, save :: nfol = 0, fkind(MAXFOL) = 0, fat(4,MAXFOL) = 0
logical, allocatable, save :: fmask(:)
logical, save :: fol_atoms = .false.
integer, save :: ctx_nat = 0
real(8), allocatable, save :: ctx_g3(:,:), ctx_M(:,:)

integer, parameter :: TS_MAX_RECALC = 5
real(8), parameter :: TS_REFRESH_GMAX = 2.0d-3

contains

subroutine ts_start(on)
logical, intent(in) :: on
ts_on = on
if (allocated(ts_vprev)) deallocate(ts_vprev)
ts_nneg = 1
ts_fresh = .true.
ts_abort = .false.
end subroutine ts_start

subroutine ts_forget_mode()
if (allocated(ts_vprev)) deallocate(ts_vprev)
ts_fresh = .true.
end subroutine ts_forget_mode

logical function ts_endgame_refresh_due(cyc, gmax, nrecalc, multi_done)
integer, intent(in) :: cyc, nrecalc
real(8), intent(in) :: gmax
logical, intent(in) :: multi_done
ts_endgame_refresh_due = ts_on .and. .not. ts_abort .and. cyc > 1 .and. ts_nneg > 1 .and. gmax < TS_REFRESH_GMAX &
                         .and. .not. multi_done .and. nrecalc < TS_MAX_RECALC
end function ts_endgame_refresh_due

logical function ts_lost_mode_recalc_due(hess_fresh, nrecalc)
logical, intent(in) :: hess_fresh
integer, intent(in) :: nrecalc
ts_lost_mode_recalc_due = ts_on .and. .not. ts_abort .and. ts_nneg == 0 .and. .not. hess_fresh .and. nrecalc < TS_MAX_RECALC
end function ts_lost_mode_recalc_due

subroutine ts_follow_init(str, nat, out_unit, ok)
character(len=*), intent(in) :: str
integer, intent(in) :: nat, out_unit
logical, intent(out) :: ok
character(len=len(str)+1) :: rest, item
integer :: p, ios, nw, k, a, b, na, ia(4)
character(len=32) :: w(8)
ok = .true.
nfol = 0 ; fol_atoms = .false.
if (allocated(fmask)) deallocate(fmask)
if (len_trim(str) == 0) return
rest = trim(str)//';'
do
   p = index(rest, ';')
   if (p == 0) exit
   item = adjustl(rest(1:p-1))
   rest = rest(p+1:)
   if (len_trim(item) == 0) cycle
   if (item(1:1) == '-' .or. item(1:1) == '+') then
      write(out_unit,'(A)') '  [ts] opt_ts_follow: no signs - list the coordinates, "'//trim(item)//'"'
      ok = .false. ; return
   endif
   call words(item, w, nw)
   select case (lower(trim(w(1))))
   case ('bond', 'angle', 'dihedral')
      k = merge(2, merge(3, 4, trim(lower(w(1))) == 'angle'), trim(lower(w(1))) == 'bond')
      if (nw /= k+1) then
         write(out_unit,'(A)') '  [ts] opt_ts_follow: "'//trim(item)//'" needs '//merge('2', merge('3','4',k==3), k==2)//' atoms'
         ok = .false. ; return
      endif
      do a = 1, k
         read(w(a+1), *, iostat=ios) ia(a)
         if (ios /= 0 .or. ia(a) < 1 .or. ia(a) > nat) then
            write(out_unit,'(A)') '  [ts] opt_ts_follow: bad atom number in "'//trim(item)//'"'
            ok = .false. ; return
         endif
      enddo
      if (nfol >= MAXFOL) then
         write(out_unit,'(A)') '  [ts] opt_ts_follow: too many items' ; ok = .false. ; return
      endif
      nfol = nfol + 1
      fkind(nfol) = k ; fat(:,nfol) = 0 ; fat(1:k,nfol) = ia(1:k)
   case ('atoms')
      if (.not. allocated(fmask)) then
         allocate(fmask(nat)) ; fmask = .false.
      endif
      do a = 2, nw
         b = index(w(a), '-')
         if (b > 1) then
            read(w(a)(1:b-1), *, iostat=ios) ia(1)
            if (ios == 0) read(w(a)(b+1:), *, iostat=ios) ia(2)
            if (ios /= 0 .or. ia(1) < 1 .or. ia(2) > nat .or. ia(1) > ia(2)) then
               write(out_unit,'(A)') '  [ts] opt_ts_follow: bad atom range "'//trim(w(a))//'"'
               ok = .false. ; return
            endif
            fmask(ia(1):ia(2)) = .true.
         else
            read(w(a), *, iostat=ios) ia(1)
            if (ios /= 0 .or. ia(1) < 1 .or. ia(1) > nat) then
               write(out_unit,'(A)') '  [ts] opt_ts_follow: bad atom number "'//trim(w(a))//'"'
               ok = .false. ; return
            endif
            fmask(ia(1)) = .true.
         endif
      enddo
      na = count(fmask)
      if (na == 0) then
         write(out_unit,'(A)') '  [ts] opt_ts_follow: empty atoms list' ; ok = .false. ; return
      endif
      fol_atoms = .true.
   case default
      write(out_unit,'(A)') '  [ts] opt_ts_follow: unknown item "'//trim(item)//'" (bond, angle, dihedral, atoms)'
      ok = .false. ; return
   end select
enddo
write(out_unit,'(A)') '  [ts] mode to follow: '//trim(str)
end subroutine ts_follow_init

subroutine ts_follow_context(nat, g3, ni, B, Ginv)
integer, intent(in) :: nat
real(8), intent(in) :: g3(3,nat)
integer, intent(in), optional :: ni
real(8), intent(in), optional :: B(:,:), Ginv(:,:)
if (nfol == 0 .and. .not. fol_atoms) return
ctx_nat = nat
if (allocated(ctx_g3)) deallocate(ctx_g3)
if (allocated(ctx_M)) deallocate(ctx_M)
allocate(ctx_g3(3,nat), source=g3)
if (present(ni) .and. present(B) .and. present(Ginv)) then
   allocate(ctx_M(3*nat,ni))
   call dgemm('T', 'N', 3*nat, ni, ni, 1.0d0, B, ni, Ginv, ni, 0.0d0, ctx_M, 3*nat)
endif
end subroutine ts_follow_context

integer function pick_follow_mode(n, V, b, out_unit)
integer, intent(in) :: n, out_unit
real(8), intent(in) :: V(n,n), b(n)
real(8) :: tc(3*ctx_nat), dx(3*ctx_nat), Q(3*ctx_nat,MAXFOL), s(3,4), sc, best, nt, nd, wf
integer :: i, j, k, a, nq, nnegc
logical :: ok
pick_follow_mode = 0
if (nfol == 0 .and. .not. fol_atoms) return
if (ctx_nat == 0 .or. .not. allocated(ctx_g3)) return
if (n /= 3*ctx_nat) then
   if (.not. allocated(ctx_M)) return
   if (size(ctx_M,2) /= n) return
endif
nq = 0
do i = 1, nfol
   s = 0.0d0
   k = fkind(i)
   select case (k)
   case (2); call bvec_stretch(ctx_g3(:,fat(1,i)), ctx_g3(:,fat(2,i)), s(:,1), s(:,2), ok)
   case (3); call bvec_bend(ctx_g3(:,fat(1,i)), ctx_g3(:,fat(2,i)), ctx_g3(:,fat(3,i)), s(:,1), s(:,2), s(:,3), ok)
   case default
      call bvec_torsion(ctx_g3(:,fat(1,i)), ctx_g3(:,fat(2,i)), ctx_g3(:,fat(3,i)), ctx_g3(:,fat(4,i)), &
                        s(:,1), s(:,2), s(:,3), s(:,4), ok)
   end select
   tc = 0.0d0
   do j = 1, k
      a = fat(j,i)
      tc(3*a-2:3*a) = s(:,j)
   enddo
   do j = 1, nq
      tc = tc - dot_product(Q(:,j), tc)*Q(:,j)
   enddo
   nt = sqrt(sum(tc*tc))
   if (nt > 1.0d-8) then
      nq = nq + 1
      Q(:,nq) = tc/nt
   endif
enddo
if (nfol > 0 .and. nq == 0) return
best = -1.0d0
nnegc = count(b < -1.0d-6)
do k = 1, n
   if (nnegc > 0 .and. b(k) >= -1.0d-6) cycle
   if (allocated(ctx_M) .and. n /= 3*ctx_nat) then
      call dgemv('N', 3*ctx_nat, n, 1.0d0, ctx_M, 3*ctx_nat, V(:,k), 1, 0.0d0, dx, 1)
   else
      dx = V(:,k)
   endif
   nd = sqrt(sum(dx*dx))
   if (nd < 1.0d-6) cycle
   sc = 1.0d0
   if (nfol > 0) then
      sc = 0.0d0
      do j = 1, nq
         sc = sc + dot_product(Q(:,j), dx)**2
      enddo
      sc = sqrt(sc) / nd
   endif
   if (fol_atoms) then
      wf = 0.0d0
      do a = 1, ctx_nat
         if (fmask(a)) wf = wf + sum(dx(3*a-2:3*a)**2)
      enddo
      sc = sc * wf / (nd*nd)
   endif
   if (sc > best) then
      best = sc ; pick_follow_mode = k
   endif
enddo
if (pick_follow_mode > 0) &
   write(out_unit,'(A,I0,A,F6.3,A,ES10.3)') '  [ts] requested motion matches mode ', pick_follow_mode, &
        '  (score ', best, ')  eigenvalue ', b(pick_follow_mode)
end function pick_follow_mode

subroutine words(sin, w, n)
character(len=*), intent(in) :: sin
character(len=*), intent(out) :: w(:)
integer, intent(out) :: n
integer :: i, j, L
n = 0 ; i = 1 ; L = len_trim(sin)
w = ''
do while (i <= L)
   if (sin(i:i) == ' ' .or. sin(i:i) == ',') then
      i = i + 1 ; cycle
   endif
   j = i
   do while (j <= L)
      if (sin(j:j) == ' ' .or. sin(j:j) == ',') exit
      j = j + 1
   enddo
   if (n < size(w)) then
      n = n + 1 ; w(n) = sin(i:j-1)
   endif
   i = j
enddo
end subroutine words

function lower(sv) result(t)
character(len=*), intent(in) :: sv
character(len=len(sv)) :: t
integer :: i, c
t = sv
do i = 1, len(sv)
   c = iachar(sv(i:i))
   if (c >= 65 .and. c <= 90) t(i:i) = achar(c+32)
enddo
end function lower

subroutine solve_prfo(n, H, g, step, out_unit)
integer, intent(in)  :: n, out_unit
real(8), intent(in)  :: H(n,n), g(n)
real(8), intent(out) :: step(n)
real(8), allocatable :: V(:,:), b(:), F(:), s(:), work(:), aug(:,:), w(:)
integer :: lwork, info, k, kf, m, nneg
real(8) :: ov, ovbest, lam_p, vlast
logical :: following
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
following = (nfol > 0 .or. fol_atoms)
if (following) then
   if (.not. allocated(ts_vprev)) then
      k = pick_follow_mode(n, V, b, out_unit)
      if (k > 0) kf = k
   endif
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
endif
nneg = count(b < -1.0d-6)
ts_nneg = nneg
if (following) then
   write(out_unit,'(A,I0,A,ES10.3,A,ES10.3,A,I0,A,F6.3)') '  [ts] following mode ', kf, &
        '  eigenvalue ', b(kf), '  lowest ', b(1), '  negative ', nneg, '  overlap  ', ovbest
else
   write(out_unit,'(A,ES10.3,A,I0)') '  [ts] climbing the most negative mode  eigenvalue ', b(1), '  negative ', nneg
endif
if (nneg == 0) write(out_unit,'(A)') '  [ts] WARNING: the Hessian has no negative eigenvalue'
if (.not. following) then
   if (nneg == 0 .and. ts_fresh) then
      ts_abort = .true.
      call ts_stop_message(out_unit, 'The Hessian has no negative eigenvalue: this geometry is not near a transition state. ' // &
           'Give opt_ts_follow (bond / angle / dihedral / atoms) to choose the motion to climb, or start nearer the saddle point.')
      step = 0.0d0
      return
   elseif (nneg >= 2 .and. ts_fresh) then
      write(out_unit,'(A,I0,A,ES10.3,A)') '  [ts] ', nneg, ' negative eigenvalues; climbing the most negative one (', b(1), &
           '). opt_ts_follow chooses another.'
   endif
endif
ts_fresh = .false.

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

if (following) then
   if (allocated(ts_vprev)) then
      if (size(ts_vprev) /= n) deallocate(ts_vprev)
   endif
   if (.not. allocated(ts_vprev)) allocate(ts_vprev(n))
   ts_vprev = V(:,kf)
endif
end subroutine solve_prfo

subroutine ts_stop_message(out_unit, msg)
integer, intent(in) :: out_unit
character(len=*), intent(in) :: msg
write(out_unit,'(A)') '  [ts] ERROR: '//msg
write(0,'(A)') 'Direwolf: transition-state search stopped - '//msg
end subroutine ts_stop_message

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
