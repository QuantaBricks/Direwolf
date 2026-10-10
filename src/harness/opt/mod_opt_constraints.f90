! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Constraints for the geometry optimizer: frozen atoms and frozen (or

module mod_opt_constraints
use mod_wilson_bvec, only: bvec_stretch, bvec_bend, bvec_torsion, &
                            val_stretch, val_bend, val_torsion, wrap_pi
implicit none
private
public :: cons_init, cons_active, cons_has_frozen, cons_update, cons_basis, cons_project_grad, cons_correct
public :: cons_max_dev
public :: orth_cols

integer, parameter :: MAXCON = 64
real(8), parameter :: ANG2BOHR = 1.88972612456506d0
real(8), parameter :: PI = 3.14159265358979323846d0
real(8), parameter :: CONS_LIM_BOND = 0.10d0
real(8), parameter :: CONS_LIM_ANGLE = 0.10d0

integer, save :: nfrozen = 0, ncon = 0
integer, allocatable, save :: frozen(:)
real(8), allocatable, save :: x0(:,:)
integer, save :: ckind(MAXCON) = 0
integer, save :: cat(4,MAXCON) = 0
real(8), save :: cfinal(MAXCON) = 0.0d0, cref(MAXCON) = 0.0d0
real(8), allocatable, save :: V(:,:)
integer, save :: nv = 0

contains

logical function cons_has_frozen()
cons_has_frozen = (nfrozen > 0)
end function cons_has_frozen

logical function cons_active()
cons_active = (nfrozen > 0 .or. ncon > 0)
end function cons_active

subroutine cons_init(freeze_str, cons_str, nat, g3, out_unit, ok)
character(len=*), intent(in) :: freeze_str, cons_str
integer, intent(in) :: nat, out_unit
real(8), intent(in) :: g3(3,nat)
logical, intent(out) :: ok
integer :: i, p, q, a, b, ios, nf, ntok, k, ia(4)
character(len=len(cons_str)) :: item
logical :: isfrozen(nat)
character(len=32) :: tok(8)
character(len=24) :: atoms_str
real(8) :: val
ok = .true.
nfrozen = 0 ; ncon = 0
if (allocated(frozen)) deallocate(frozen, x0)
if (allocated(V)) deallocate(V)
nv = 0

isfrozen = .false.
item = freeze_str
do i = 1, len_trim(item)
   if (item(i:i) == ',') item(i:i) = ' '
enddo
p = 1
do while (p <= len_trim(item))
   do while (p <= len_trim(item))
      if (item(p:p) /= ' ') exit
      p = p + 1
   enddo
   if (p > len_trim(item)) exit
   q = p
   do while (q <= len_trim(item))
      if (item(q:q) == ' ') exit
      q = q + 1
   enddo
   k = index(item(p:q-1), '-')
   if (k > 0) then
      read(item(p:p+k-2), *, iostat=ios) a
      if (ios == 0) read(item(p+k:q-1), *, iostat=ios) b
   else
      read(item(p:q-1), *, iostat=ios) a
      b = a
   endif
   if (ios /= 0 .or. a < 1 .or. b > nat .or. a > b) then
      write(out_unit,'(A)') '[OPT] opt_freeze: cannot read "'//item(p:q-1)//'"'
      ok = .false. ; return
   endif
   isfrozen(a:b) = .true.
   p = q
enddo
nf = count(isfrozen)
if (nf > 0) then
   allocate(frozen(nf), x0(3,nf))
   k = 0
   do i = 1, nat
      if (isfrozen(i)) then
         k = k + 1 ; frozen(k) = i ; x0(:,k) = g3(:,i)
      endif
   enddo
   nfrozen = nf
endif

p = 1
do while (p <= len_trim(cons_str))
   q = index(cons_str(p:), ';')
   if (q == 0) then
      item = cons_str(p:)
      p = len_trim(cons_str) + 1
   else
      item = cons_str(p:p+q-2)
      p = p + q
   endif
   if (len_trim(item) == 0) cycle
   tok = ''
   ntok = 0
   call split_words(item, tok, ntok)
   select case (trim(lower(tok(1))))
   case ('bond');     k = 2
   case ('angle');    k = 3
   case ('dihedral'); k = 4
   case default
      write(out_unit,'(A)') '[OPT] opt_constraints: unknown item "'//trim(item)//'"'
      ok = .false. ; return
   end select
   if (ntok < 1 + k) then
      write(out_unit,'(A)') '[OPT] opt_constraints: too few atoms in "'//trim(item)//'"'
      ok = .false. ; return
   endif
   ia = 0
   do i = 1, k
      read(tok(1+i), *, iostat=ios) ia(i)
      if (ios /= 0 .or. ia(i) < 1 .or. ia(i) > nat) then
         write(out_unit,'(A)') '[OPT] opt_constraints: bad atom number in "'//trim(item)//'"'
         ok = .false. ; return
      endif
   enddo
   if (ncon >= MAXCON) then
      write(out_unit,'(A,I0)') '[OPT] opt_constraints: at most ', MAXCON
      ok = .false. ; return
   endif
   ncon = ncon + 1
   ckind(ncon) = k
   cat(:,ncon) = ia
   cref(ncon) = value_of(ncon, g3, nat)
   cfinal(ncon) = cref(ncon)
   if (ntok >= 2 + k) then
      read(tok(2+k), *, iostat=ios) val
      if (ios /= 0) then
         write(out_unit,'(A)') '[OPT] opt_constraints: bad value in "'//trim(item)//'"'
         ok = .false. ; return
      endif
      if (k == 2) then
         cfinal(ncon) = val * ANG2BOHR
      else
         cfinal(ncon) = val * PI / 180.0d0
      endif
   endif
enddo
if (nfrozen > 0) write(out_unit,'(A,I0,A)') '[OPT]  frozen atoms: ', nfrozen, ' (opt_freeze)'
do i = 1, ncon
   write(atoms_str,'(4(I0,1X))') cat(1:ckind(i),i)
   write(out_unit,'(A,I0,A,A,A,F10.4,A,F10.4)') '[OPT]  constraint ', i, ': ', &
        trim(kname(ckind(i)))//' '//trim(atoms_str), '   start ', shown(i, cref(i)), '   target ', shown(i, cfinal(i))
enddo
end subroutine cons_init

character(len=8) function kname(k)
integer, intent(in) :: k
select case (k)
case (2); kname = 'bond'
case (3); kname = 'angle'
case default; kname = 'dihedral'
end select
end function kname

real(8) function shown(i, v)
integer, intent(in) :: i
real(8), intent(in) :: v
if (ckind(i) == 2) then
   shown = v / ANG2BOHR
else
   shown = v * 180.0d0 / PI
endif
end function shown

function lower(s) result(t)
character(len=*), intent(in) :: s
character(len=len(s)) :: t
integer :: i, c
t = s
do i = 1, len(s)
   c = iachar(s(i:i))
   if (c >= 65 .and. c <= 90) t(i:i) = achar(c+32)
enddo
end function lower

subroutine split_words(s, w, n)
character(len=*), intent(in) :: s
character(len=*), intent(out) :: w(:)
integer, intent(out) :: n
integer :: p, q
n = 0 ; w = ''
p = 1
do while (p <= len_trim(s) .and. n < size(w))
   if (s(p:p) == ' ') then
      p = p + 1 ; cycle
   endif
   q = p
   do while (q <= len_trim(s))
      if (s(q:q) == ' ') exit
      q = q + 1
   enddo
   n = n + 1 ; w(n) = s(p:q-1)
   p = q
enddo
end subroutine split_words

real(8) function value_of(i, g3, nat)
integer, intent(in) :: i, nat
real(8), intent(in) :: g3(3,nat)
select case (ckind(i))
case (2); value_of = val_stretch(g3(:,cat(1,i)), g3(:,cat(2,i)))
case (3); value_of = val_bend(g3(:,cat(1,i)), g3(:,cat(2,i)), g3(:,cat(3,i)))
case default
   value_of = val_torsion(g3(:,cat(1,i)), g3(:,cat(2,i)), g3(:,cat(3,i)), g3(:,cat(4,i)))
end select
end function value_of

subroutine jac_row(i, g3, nat, row)
integer, intent(in) :: i, nat
real(8), intent(in) :: g3(3,nat)
real(8), intent(out) :: row(3*nat)
real(8) :: s(3,4)
logical :: ok
integer :: j, a
row = 0.0d0
s = 0.0d0
select case (ckind(i))
case (2); call bvec_stretch(g3(:,cat(1,i)), g3(:,cat(2,i)), s(:,1), s(:,2), ok)
case (3); call bvec_bend(g3(:,cat(1,i)), g3(:,cat(2,i)), g3(:,cat(3,i)), s(:,1), s(:,2), s(:,3), ok)
case default
   call bvec_torsion(g3(:,cat(1,i)), g3(:,cat(2,i)), g3(:,cat(3,i)), g3(:,cat(4,i)), &
                     s(:,1), s(:,2), s(:,3), s(:,4), ok)
end select
do j = 1, ckind(i)
   a = cat(j,i)
   row(3*a-2:3*a) = s(:,j)
enddo
end subroutine jac_row

subroutine orth_cols(A, n, m, Q, nq)
integer, intent(in) :: n, m
real(8), intent(in) :: A(n,m)
real(8), intent(out) :: Q(n,m)
integer, intent(out) :: nq
real(8) :: w(n), nrm, nrm0
integer :: k, j, pass
nq = 0
do k = 1, m
   w = A(:,k)
   nrm0 = sqrt(sum(w*w))
   if (nrm0 < 1.0d-12) cycle
   do pass = 1, 2
      do j = 1, nq
         w = w - dot_product(Q(:,j), w) * Q(:,j)
      enddo
   enddo
   nrm = sqrt(sum(w*w))
   if (nrm > 1.0d-8 * nrm0) then
      nq = nq + 1
      Q(:,nq) = w / nrm
   endif
enddo
end subroutine orth_cols

subroutine cons_update(nat, g3, ric_mode)
integer, intent(in) :: nat
real(8), intent(in) :: g3(3,nat)
logical, intent(in) :: ric_mode
real(8), allocatable :: A(:,:), T(:,:), Tq(:,:), Nn(:,:), Wn(:,:)
real(8) :: Tfull(3*nat,6)
integer :: i, m, k, nf3, ntr, nn_, ax, j
m = 3*nfrozen + ncon
if (allocated(V)) deallocate(V)
if (m == 0) then
   nv = 0
   return
endif
allocate(A(3*nat,m), V(3*nat,m))
A = 0.0d0
k = 0
if (ric_mode .and. nfrozen > 0) then
   nf3 = 3*nfrozen
   call tr_vectors(nat, g3, Tfull)
   allocate(T(nf3,6), Tq(nf3,6), Nn(nf3,nf3), Wn(nf3,nf3))
   do i = 1, nfrozen
      T(3*i-2:3*i,:) = Tfull(3*frozen(i)-2:3*frozen(i),:)
   enddo
   call orth_cols(T, nf3, 6, Tq, ntr)
   Nn = 0.0d0
   do j = 1, nf3
      Nn(j,j) = 1.0d0
   enddo
   do j = 1, ntr
      do i = 1, nf3
         do ax = 1, nf3
            Nn(ax,i) = Nn(ax,i) - Tq(ax,j)*Tq(i,j)
         enddo
      enddo
   enddo
   call orth_cols(Nn, nf3, nf3, Wn, nn_)
   do j = 1, nn_
      k = k + 1
      do i = 1, nfrozen
         A(3*frozen(i)-2:3*frozen(i), k) = Wn(3*i-2:3*i, j)
      enddo
   enddo
else
   do i = 1, nfrozen
      A(3*frozen(i)-2, k+1) = 1.0d0
      A(3*frozen(i)-1, k+2) = 1.0d0
      A(3*frozen(i),   k+3) = 1.0d0
      k = k + 3
   enddo
endif
do i = 1, ncon
   k = k + 1
   call jac_row(i, g3, nat, A(:,k))
enddo
call orth_cols(A, 3*nat, m, V, nv)
end subroutine cons_update

subroutine cons_basis(nat, g3, ric_mode, Vout, nvout)
integer, intent(in) :: nat
real(8), intent(in) :: g3(3,nat)
logical, intent(in) :: ric_mode
real(8), allocatable, intent(out) :: Vout(:,:)
integer, intent(out) :: nvout
call cons_update(nat, g3, ric_mode)
nvout = nv
allocate(Vout(3*nat, max(nv,1)))
if (nv > 0) Vout(:,1:nv) = V(:,1:nv)
end subroutine cons_basis

subroutine cons_project_grad(nat, g3, grad)
integer, intent(in) :: nat
real(8), intent(in) :: g3(3,nat)
real(8), intent(inout) :: grad(3*nat)
integer :: j
if (.not. cons_active()) return
call cons_update(nat, g3, .false.)
do j = 1, nv
   grad = grad - dot_product(V(:,j), grad) * V(:,j)
enddo
end subroutine cons_project_grad

subroutine tr_vectors(nat, g3, Tm)
integer, intent(in) :: nat
real(8), intent(in) :: g3(3,nat)
real(8), intent(out) :: Tm(3*nat,6)
real(8) :: com(3), d(3), raw(3*nat,6)
integer :: i, ax, nq
com = sum(g3, dim=2) / dble(nat)
raw = 0.0d0
do ax = 1, 3
   do i = 1, nat
      raw(3*(i-1)+ax, ax) = 1.0d0
      d = g3(:,i) - com
      select case (ax)
      case (1); raw(3*(i-1)+2, 3+ax) = -d(3); raw(3*(i-1)+3, 3+ax) =  d(2)
      case (2); raw(3*(i-1)+1, 3+ax) =  d(3); raw(3*(i-1)+3, 3+ax) = -d(1)
      case (3); raw(3*(i-1)+1, 3+ax) = -d(2); raw(3*(i-1)+2, 3+ax) =  d(1)
      end select
   enddo
enddo
Tm = raw
end subroutine tr_vectors

subroutine align_frozen(nat, g3)
integer, intent(in) :: nat
real(8), intent(inout) :: g3(3,nat)
real(8) :: Tm(3*nat,6), TF(3*nfrozen,6), rF(3*nfrozen), M(6,6), c(6), dxv(3*nat)
integer :: it, i, info, ipiv(6)
if (nfrozen == 0) return
do it = 1, 6
   call tr_vectors(nat, g3, Tm)
   do i = 1, nfrozen
      TF(3*i-2:3*i,:) = Tm(3*frozen(i)-2:3*frozen(i),:)
      rF(3*i-2:3*i) = g3(:,frozen(i)) - x0(:,i)
   enddo
   if (maxval(abs(rF)) < 1.0d-12) exit
   M = matmul(transpose(TF), TF)
   do i = 1, 6
      M(i,i) = M(i,i) + 1.0d-9
   enddo
   c = -matmul(transpose(TF), rF)
   call dgesv(6, 1, M, 6, ipiv, c, 6, info)
   if (info /= 0) exit
   dxv = matmul(Tm, c)
   g3 = g3 + reshape(dxv, [3,nat])
enddo
end subroutine align_frozen

real(8) function cons_max_dev(nat, g3)
integer, intent(in) :: nat
real(8), intent(in) :: g3(3,nat)
real(8) :: d
integer :: i
cons_max_dev = 0.0d0
do i = 1, ncon
   d = value_of(i, g3, nat) - cfinal(i)
   if (ckind(i) == 4) d = wrap_pi(d)
   cons_max_dev = max(cons_max_dev, abs(d))
enddo
end function cons_max_dev

subroutine cons_correct(nat, g3, out_unit)
integer, intent(in) :: nat, out_unit
real(8), intent(inout) :: g3(3,nat)
real(8), allocatable :: C(:,:), r(:), CC(:,:), lam(:)
real(8) :: tgt(MAXCON), d, lim, rmax
integer :: i, k, m, it, info, ax
integer, allocatable :: ipiv(:)
if (.not. cons_active()) return
call align_frozen(nat, g3)
do i = 1, ncon
   d = cfinal(i) - cref(i)
   if (ckind(i) == 4) d = wrap_pi(d)
   lim = merge(CONS_LIM_BOND, CONS_LIM_ANGLE, ckind(i) == 2)
   tgt(i) = cref(i) + max(-lim, min(lim, d))
enddo
m = 3*nfrozen + ncon
allocate(C(m,3*nat), r(m), CC(m,m), lam(m), ipiv(m))
do it = 1, 12
   C = 0.0d0
   k = 0
   do i = 1, nfrozen
      do ax = 1, 3
         k = k + 1
         C(k, 3*frozen(i)-3+ax) = 1.0d0
         r(k) = g3(ax, frozen(i)) - x0(ax, i)
      enddo
   enddo
   do i = 1, ncon
      k = k + 1
      call jac_row(i, g3, nat, C(k,:))
      d = value_of(i, g3, nat) - tgt(i)
      if (ckind(i) == 4) d = wrap_pi(d)
      r(k) = d
   enddo
   rmax = maxval(abs(r))
   if (rmax < 1.0d-11) exit
   CC = matmul(C, transpose(C))
   do i = 1, m
      CC(i,i) = CC(i,i) + 1.0d-10
   enddo
   lam = r
   call dgesv(m, 1, CC, m, ipiv, lam, m, info)
   if (info /= 0) exit
   do i = 1, nat
      do k = 1, 3
         g3(k,i) = g3(k,i) - dot_product(C(:,3*(i-1)+k), lam)
      enddo
   enddo
enddo
do i = 1, ncon
   cref(i) = tgt(i)
   write(out_unit,'(A,I0,A,A,F10.5,A,F10.5,A)') '  [cons] ', i, ' ', trim(kname(ckind(i)))//' ', &
        shown(i, value_of(i, g3, nat)), '  (final ', shown(i, cfinal(i)), ')'
enddo
end subroutine cons_correct

end module mod_opt_constraints
