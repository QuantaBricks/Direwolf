! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Rigid-body coordinates of covalent fragments, decoupled from the

module mod_opt_fragment
use mod_opt_min, only: bfgs_update, solve_rfo
implicit none
private
public :: frag_init, frag_reset, frag_is_active, frag_nfrag, frag_expected_intra
public :: frag_rigid_step, frag_update, frag_get_fid, frag_trust_update, frag_trust_shrink, frag_trust_info
public :: frag_set_model, lj_params, frag_classical_opt

integer, parameter :: NC = 6
real(8) :: KT0 = 2.0d-2
real(8) :: KR0 = 2.0d-2
real(8), parameter :: ROT_ARM = 2.0d0

real(8), parameter :: TR0 = 1.0d0, TR_MIN = 0.1d0, TR_MAX = 4.0d0
real(8) :: trust_r = TR0, last_snr = 0.0d0
logical :: last_hit = .false.

logical :: active = .false.
integer :: nat_f = 0, nfrag = 0, nr = 0
integer, allocatable :: fid(:)
real(8), allocatable :: Hr(:,:)
real(8), allocatable :: u_pend(:), gu_at(:)
logical :: have_pend = .false.

contains

subroutine frag_set_model(nat, Z, coord)
integer, intent(in) :: nat, Z(nat)
real(8), intent(in) :: coord(3,nat)
real(8), parameter :: STEP_T = 1.0d-2, STEP_R = 5.0d-3
real(8) :: Hw(nr,nr), hstep(nr), up(nr), Epp, Epm, Emp, Emm, w(nr), work(8*nr)
real(8), allocatable :: V(:,:)
integer :: j, k, info, i
if (.not. active) return
do i = 1, nfrag
   hstep(NC*(i-1)+1:NC*(i-1)+3) = STEP_T
   hstep(NC*(i-1)+4:NC*(i-1)+6) = STEP_R
enddo
do j = 1, nr
   do k = j, nr
      up = 0.0d0 ; up(j) = up(j) + hstep(j) ; up(k) = up(k) + hstep(k)
      Epp = wall_energy(Z, rigid_move(coord, up))
      up = 0.0d0 ; up(j) = up(j) + hstep(j) ; up(k) = up(k) - hstep(k)
      Epm = wall_energy(Z, rigid_move(coord, up))
      up = 0.0d0 ; up(j) = up(j) - hstep(j) ; up(k) = up(k) + hstep(k)
      Emp = wall_energy(Z, rigid_move(coord, up))
      up = 0.0d0 ; up(j) = up(j) - hstep(j) ; up(k) = up(k) - hstep(k)
      Emm = wall_energy(Z, rigid_move(coord, up))
      Hw(j,k) = (Epp - Epm - Emp + Emm) / (4.0d0*hstep(j)*hstep(k))
      Hw(k,j) = Hw(j,k)
   enddo
enddo
allocate(V(nr,nr))
V = Hw
call dsyev('V', 'U', nr, V, nr, w, work, 8*nr, info)
if (info /= 0) return
Hw = 0.0d0
do k = 1, nr
   if (w(k) > 0.0d0) then
      do j = 1, nr
         do i = 1, nr
            Hw(i,j) = Hw(i,j) + w(k)*V(i,k)*V(j,k)
         enddo
      enddo
   endif
enddo
Hr = Hr + Hw
end subroutine frag_set_model

subroutine frag_classical_opt(nat, Z, q, coord, maxit, E_start, E_end, nit, dmax)
integer, intent(in)    :: nat, Z(nat), maxit
real(8), intent(in)    :: q(nat)
real(8), intent(inout) :: coord(3,nat)
real(8), intent(out)   :: E_start, E_end, dmax
integer, intent(out)   :: nit
real(8) :: Hc(nr,nr), g(3,nat), gn(3,nat), gu(nr), gun(nr), u(nr), Hu(nr), cnew(3,nat), c0(3,nat)
real(8) :: E, En, pred, snr, tr, ratio
integer :: it, k, i
E_start = 0.0d0 ; E_end = 0.0d0 ; nit = 0 ; dmax = 0.0d0
if (.not. active) return
Hc = 0.0d0
do k = 1, nr
   Hc(k,k) = KT0
enddo
c0 = coord
tr = 0.5d0
call cl_eg(Z, q, coord, E, g)
E_start = E
do it = 1, maxit
   nit = it
   call rigid_gradient(coord, reshape(g, [3*nat]), gu)
   if (maxval(abs(gu)) < 2.0d-5) exit
   call solve_rfo(nr, Hc, gu, u, 6)
   snr = 0.0d0
   do k = 1, nfrag
      snr = snr + sum(u(NC*(k-1)+1:NC*(k-1)+3)**2) + ROT_ARM**2 * sum(u(NC*(k-1)+4:NC*(k-1)+6)**2)
   enddo
   snr = sqrt(snr)
   if (snr > tr .and. snr > 0.0d0) then
      u = u * (tr / snr) ; snr = tr
   endif
   cnew = rigid_move(coord, u)
   call cl_eg(Z, q, cnew, En, gn)
   call dgemv('N', nr, nr, 1.0d0, Hc, nr, u, 1, 0.0d0, Hu, 1)
   pred = dot_product(gu, u) + 0.5d0 * dot_product(u, Hu)
   if (En < E) then
      call rigid_gradient(cnew, reshape(gn, [3*nat]), gun)
      call bfgs_update(nr, Hc, u, gun - gu)
      ratio = (En - E) / min(pred, -1.0d-14)
      coord = cnew ; E = En ; g = gn
      if (ratio > 0.75d0) then
         tr = min(2.0d0*tr, 5.0d0)
      else if (ratio < 0.25d0) then
         tr = max(0.5d0*snr, 1.0d-4)
      endif
   else
      tr = 0.25d0 * snr
      if (tr < 1.0d-5) exit
   endif
enddo
E_end = E
do i = 1, nat
   dmax = max(dmax, sqrt(sum((coord(:,i) - c0(:,i))**2)))
enddo
end subroutine frag_classical_opt

subroutine cl_eg(Z, q, coord, E, g)
integer, intent(in)  :: Z(nat_f)
real(8), intent(in)  :: q(nat_f), coord(3,nat_f)
real(8), intent(out) :: E, g(3,nat_f)
real(8), parameter :: ANG2BOHR = 1.0d0/0.529177210903d0, KCAL2HA = 1.0d0/627.5094740631d0
real(8), parameter :: R_HX = 1.55d0*ANG2BOHR, R_HH = 1.80d0*ANG2BOHR, K_WALL = 0.2d0
integer :: i, j
real(8) :: si, ei, sj, ej, sg, ep, d(3), r, s6, dEdr, r0
logical :: hi, hj
E = 0.0d0 ; g = 0.0d0
do i = 1, nat_f
   call lj_params(Z(i), si, ei)
   hi = (Z(i) == 1)
   do j = i+1, nat_f
      if (fid(i) == fid(j)) cycle
      call lj_params(Z(j), sj, ej)
      hj = (Z(j) == 1)
      d = coord(:,i) - coord(:,j)
      r = sqrt(sum(d*d))
      E = E + q(i)*q(j)/r
      dEdr = -q(i)*q(j)/r**2
      g(:,i) = g(:,i) + dEdr*d/r ; g(:,j) = g(:,j) - dEdr*d/r
      if (ei > 0.0d0 .and. ej > 0.0d0) then
         sg = 0.5d0*(si + sj)*ANG2BOHR ; ep = sqrt(ei*ej)*KCAL2HA
         s6 = (sg/r)**6
         E = E + 4.0d0*ep*(s6*s6 - s6)
         dEdr = 4.0d0*ep*(-12.0d0*s6*s6 + 6.0d0*s6)/r
         g(:,i) = g(:,i) + dEdr*d/r ; g(:,j) = g(:,j) - dEdr*d/r
      endif
      if (hi .or. hj) then
         r0 = merge(R_HH, R_HX, hi .and. hj)
         if (r < r0) then
            E = E + K_WALL*(r0 - r)**2
            dEdr = -2.0d0*K_WALL*(r0 - r)
            g(:,i) = g(:,i) + dEdr*d/r ; g(:,j) = g(:,j) - dEdr*d/r
         endif
      endif
   enddo
enddo
end subroutine cl_eg

function rigid_move(coord, u) result(cn)
real(8), intent(in) :: coord(3,nat_f), u(nr)
real(8) :: cn(3,nat_f), cen(3,nfrag), r(3), th(3), ax(3), ang, c, s
integer :: i, k, cnt(nfrag)
cen = 0.0d0 ; cnt = 0
do i = 1, nat_f
   cen(:,fid(i)) = cen(:,fid(i)) + coord(:,i) ; cnt(fid(i)) = cnt(fid(i)) + 1
enddo
do k = 1, nfrag
   cen(:,k) = cen(:,k) / real(cnt(k),8)
enddo
do i = 1, nat_f
   k = fid(i)
   th = u(NC*(k-1)+4:NC*(k-1)+6) ; ang = sqrt(sum(th*th))
   r = coord(:,i) - cen(:,k)
   if (ang > 1.0d-12) then
      ax = th/ang ; c = cos(ang) ; s = sin(ang)
      r = r*c + [ax(2)*r(3)-ax(3)*r(2), ax(3)*r(1)-ax(1)*r(3), ax(1)*r(2)-ax(2)*r(1)]*s + ax*dot_product(ax,r)*(1.0d0-c)
   endif
   cn(:,i) = cen(:,k) + r + u(NC*(k-1)+1:NC*(k-1)+3)
enddo
end function rigid_move

real(8) function wall_energy(Z, coord)
integer, intent(in) :: Z(nat_f)
real(8), intent(in) :: coord(3,nat_f)
real(8), parameter :: ANG2BOHR = 1.0d0/0.529177210903d0, KCAL2HA = 1.0d0/627.5094740631d0
integer :: i, j
real(8) :: si, ei, sj, ej, sg, ep, r
wall_energy = 0.0d0
do i = 1, nat_f
   call lj_params(Z(i), si, ei)
   if (ei <= 0.0d0) cycle
   do j = i+1, nat_f
      if (fid(i) == fid(j)) cycle
      call lj_params(Z(j), sj, ej)
      if (ej <= 0.0d0) cycle
      sg = 0.5d0*(si + sj)*ANG2BOHR ; ep = sqrt(ei*ej)*KCAL2HA
      r = sqrt(sum((coord(:,i) - coord(:,j))**2))
      wall_energy = wall_energy + 4.0d0*ep*(sg/r)**12
   enddo
enddo
end function wall_energy

subroutine lj_params(Z, sigma, eps)
integer, intent(in)  :: Z
real(8), intent(out) :: sigma, eps
select case (Z)
case (1)  ; sigma = 0.0d0 ; eps = 0.0d0
case (6)  ; sigma = 3.50d0 ; eps = 0.066d0
case (7)  ; sigma = 3.25d0 ; eps = 0.170d0
case (8)  ; sigma = 3.07d0 ; eps = 0.170d0
case (9)  ; sigma = 2.95d0 ; eps = 0.061d0
case (15) ; sigma = 3.74d0 ; eps = 0.200d0
case (16) ; sigma = 3.55d0 ; eps = 0.250d0
case (17) ; sigma = 3.50d0 ; eps = 0.265d0
case default ; sigma = 3.50d0 ; eps = 0.100d0
end select
end subroutine lj_params

subroutine frag_get_fid(nat, fidout)
integer, intent(in)  :: nat
integer, intent(out) :: fidout(nat)
fidout = 1
if (active .and. allocated(fid)) fidout = fid
end subroutine frag_get_fid

logical function frag_is_active()
frag_is_active = active
end function frag_is_active

integer function frag_nfrag()
frag_nfrag = nfrag
end function frag_nfrag

subroutine frag_reset()
active = .false. ; nfrag = 0 ; nr = 0 ; have_pend = .false.
if (allocated(fid)) deallocate(fid)
if (allocated(Hr)) deallocate(Hr)
if (allocated(u_pend)) deallocate(u_pend)
if (allocated(gu_at)) deallocate(gu_at)
end subroutine frag_reset

subroutine frag_init(nat, nbond, bond)
integer, intent(in) :: nat, nbond, bond(2,*)
integer :: par(nat), i, j, b, ra, rb, lab(nat)
character(len=32) :: env
call frag_reset()
call get_environment_variable('TMP_HR', env)
if (len_trim(env) > 0) read(env,*) KT0, KR0
do i = 1, nat
   par(i) = i
enddo
do b = 1, nbond
   ra = root(bond(1,b)) ; rb = root(bond(2,b))
   if (ra /= rb) par(ra) = rb
enddo
lab = 0 ; nfrag = 0
allocate(fid(nat))
do i = 1, nat
   j = root(i)
   if (lab(j) == 0) then
      nfrag = nfrag + 1 ; lab(j) = nfrag
   endif
   fid(i) = lab(j)
enddo
if (nfrag < 2) then
   call frag_reset() ; return
endif
nat_f = nat
nr = NC*nfrag
allocate(Hr(nr,nr), u_pend(nr), gu_at(nr))
Hr = 0.0d0
do i = 1, nfrag
   do j = 1, 3
      Hr(NC*(i-1)+j, NC*(i-1)+j) = KT0
      Hr(NC*(i-1)+3+j, NC*(i-1)+3+j) = KR0
   enddo
enddo
u_pend = 0.0d0 ; gu_at = 0.0d0
trust_r = TR0 ; last_snr = 0.0d0 ; last_hit = .false.
active = .true.
contains
   integer function root(a)
   integer, intent(in) :: a
   integer :: x
   x = a
   do while (par(x) /= x)
      par(x) = par(par(x))
      x = par(x)
   enddo
   root = x
   end function root
end subroutine frag_init

integer function frag_expected_intra(nat, coord)
integer, intent(in) :: nat
real(8), intent(in) :: coord(3,nat)
integer :: k, i, n, first, second
real(8) :: u(3), w(3), cr(3), nu
logical :: linear
frag_expected_intra = 0
if (.not. active) return
do k = 1, nfrag
   n = count(fid == k)
   if (n <= 1) cycle
   first = 0 ; second = 0 ; linear = .true.
   do i = 1, nat
      if (fid(i) /= k) cycle
      if (first == 0) then
         first = i
      else if (second == 0) then
         second = i
      else
         u = coord(:,second) - coord(:,first) ; w = coord(:,i) - coord(:,first)
         cr(1) = u(2)*w(3)-u(3)*w(2) ; cr(2) = u(3)*w(1)-u(1)*w(3) ; cr(3) = u(1)*w(2)-u(2)*w(1)
         nu = sqrt(sum(u*u))
         if (sqrt(sum(cr*cr)) > 1.0d-3 * max(nu, 1.0d0)**2) linear = .false.
      endif
   enddo
   frag_expected_intra = frag_expected_intra + 3*n - merge(5, 6, linear)
enddo
end function frag_expected_intra

subroutine rigid_gradient(coord3, grad, gu)
real(8), intent(in)  :: coord3(3,nat_f), grad(3*nat_f)
real(8), intent(out) :: gu(nr)
real(8) :: cen(3,nfrag), r(3), g(3)
integer :: i, k, cnt(nfrag)
cen = 0.0d0 ; cnt = 0
do i = 1, nat_f
   cen(:,fid(i)) = cen(:,fid(i)) + coord3(:,i) ; cnt(fid(i)) = cnt(fid(i)) + 1
enddo
do k = 1, nfrag
   cen(:,k) = cen(:,k) / real(cnt(k),8)
enddo
gu = 0.0d0
do i = 1, nat_f
   k = fid(i)
   g = grad(3*(i-1)+1:3*(i-1)+3)
   r = coord3(:,i) - cen(:,k)
   gu(NC*(k-1)+1:NC*(k-1)+3) = gu(NC*(k-1)+1:NC*(k-1)+3) + g
   gu(NC*(k-1)+4) = gu(NC*(k-1)+4) + r(2)*g(3) - r(3)*g(2)
   gu(NC*(k-1)+5) = gu(NC*(k-1)+5) + r(3)*g(1) - r(1)*g(3)
   gu(NC*(k-1)+6) = gu(NC*(k-1)+6) + r(1)*g(2) - r(2)*g(1)
enddo
end subroutine rigid_gradient

subroutine frag_rigid_step(coord_g, grad, coord3, trust_unused, dxr, pred_r, snr)
real(8), intent(in)  :: coord_g(3,nat_f), grad(3*nat_f), coord3(3,nat_f), trust_unused
real(8), intent(out) :: dxr(3*nat_f), pred_r, snr
real(8) :: gu(nr), u(nr), Hu(nr), cen(3,nfrag), x(3), th(3), ang, ax(3), c, s, r(3), dm, d(3), fac
integer :: i, k, cnt(nfrag), q
dxr = 0.0d0 ; pred_r = 0.0d0 ; snr = 0.0d0
if (.not. active) return
call rigid_gradient(coord_g, grad, gu)
call solve_rfo(nr, Hr, gu, u, 6)
snr = 0.0d0
do k = 1, nfrag
   snr = snr + sum(u(NC*(k-1)+1:NC*(k-1)+3)**2) + ROT_ARM**2 * sum(u(NC*(k-1)+4:NC*(k-1)+6)**2)
enddo
snr = sqrt(snr)
last_hit = .false.
if (snr > trust_r .and. snr > 0.0d0) then
   u = u * (trust_r / snr) ; snr = trust_r ; last_hit = .true.
endif
cen = 0.0d0 ; cnt = 0
do i = 1, nat_f
   cen(:,fid(i)) = cen(:,fid(i)) + coord3(:,i) ; cnt(fid(i)) = cnt(fid(i)) + 1
enddo
do k = 1, nfrag
   cen(:,k) = cen(:,k) / real(cnt(k),8)
enddo
call displace(u, dxr)
dm = 0.0d0
do i = 1, nat_f
   dm = max(dm, sqrt(sum(dxr(3*(i-1)+1:3*(i-1)+3)**2)))
enddo
if (dm > 2.0d0*trust_r .and. dm > 0.0d0) then
   fac = 2.0d0*trust_r / dm
   u = u * fac ; snr = snr * fac ; last_hit = .true.
   call displace(u, dxr)
endif
call dgemv('N', nr, nr, 1.0d0, Hr, nr, u, 1, 0.0d0, Hu, 1)
pred_r = dot_product(gu, u) + 0.5d0 * dot_product(u, Hu)
u_pend = u ; gu_at = gu ; have_pend = .true. ; last_snr = snr

contains
   subroutine displace(uu, dd)
   real(8), intent(in)  :: uu(nr)
   real(8), intent(out) :: dd(3*nat_f)
   integer :: ia, kk
   do ia = 1, nat_f
      kk = fid(ia)
      th = uu(NC*(kk-1)+4:NC*(kk-1)+6)
      ang = sqrt(sum(th*th))
      r = coord3(:,ia) - cen(:,kk)
      if (ang > 1.0d-12) then
         ax = th / ang ; c = cos(ang) ; s = sin(ang)
         x = r*c + cross(ax, r)*s + ax*dot_product(ax, r)*(1.0d0 - c)
      else
         x = r
      endif
      d = x - r + uu(NC*(kk-1)+1:NC*(kk-1)+3)
      dd(3*(ia-1)+1:3*(ia-1)+3) = d
   enddo
   end subroutine displace
   function cross(a, b) result(cr)
   real(8), intent(in) :: a(3), b(3)
   real(8) :: cr(3)
   cr(1) = a(2)*b(3)-a(3)*b(2) ; cr(2) = a(3)*b(1)-a(1)*b(3) ; cr(3) = a(1)*b(2)-a(2)*b(1)
   end function cross
end subroutine frag_rigid_step

subroutine frag_trust_update(ratio)
real(8), intent(in) :: ratio
if (.not. active) return
if (ratio < 0.25d0) then
   trust_r = max(0.5d0*last_snr, TR_MIN)
else if (ratio > 0.75d0 .and. last_hit) then
   trust_r = min(2.0d0*trust_r, TR_MAX)
endif
end subroutine frag_trust_update

subroutine frag_trust_info(tr, snr, hit)
real(8), intent(out) :: tr, snr
logical, intent(out) :: hit
tr = trust_r ; snr = last_snr ; hit = last_hit
end subroutine frag_trust_info

subroutine frag_trust_shrink()
if (.not. active) return
trust_r = max(0.25d0*last_snr, TR_MIN)
end subroutine frag_trust_shrink

subroutine frag_update(coord_g, grad)
real(8), intent(in) :: coord_g(3,nat_f), grad(3*nat_f)
real(8) :: gu(nr)
if (.not. active .or. .not. have_pend) return
call rigid_gradient(coord_g, grad, gu)
call bfgs_update(nr, Hr, u_pend, gu - gu_at)
have_pend = .false.
end subroutine frag_update

end module mod_opt_fragment
