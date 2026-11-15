! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Relaxed scan: constraints = 'bond 1 2 1.60 1.80 10' (also angle / dihedral; Angstrom or

module mod_opt_scan
use mod_engine_input_types, only: engine_input_t
use mod_geomopt, only: geomopt_run, geomopt_last_E, geomopt_last_ncyc, geomopt_last_ok
use mod_wilson_bvec, only: val_stretch, val_bend, val_torsion
use mod_engine_input_elements, only: elem_to_name
use mod_json_out, only: json_scan
implicit none
private
public :: scan_or_opt_run

real(8), parameter :: PI = 3.14159265358979323846d0
real(8), parameter :: HARTREE_TO_KCAL = 627.5094740631d0
real(8), parameter :: ANG2BOHR = 1.88972612456506d0

contains

subroutine scan_or_opt_run(spec, outfile, out_unit)
type(engine_input_t), intent(inout) :: spec
character(len=*), intent(in) :: outfile
integer, intent(in) :: out_unit

integer, parameter :: MAXD = 3
character(len=len(spec%opt_constraints)) :: rest, cons
character(len=16) :: kind(MAXD)
character(len=64) :: names(MAXD)
integer :: ia(4,MAXD), kat(MAXD), np(MAXD), nd, ntot, t, j, i, nat, u, ios, nfail, npt, k, base, q, par, d(MAXD)
real(8) :: v0(MAXD), v1(MAXD), eref
logical :: found, restart0
real(8), allocatable :: E(:), erel(:), tgt(:,:), val(:,:), Vx(:,:,:)
integer, allocatable :: ncyc(:), gidx(:,:)
logical, allocatable :: ok(:)
character(len=300) :: scanfile
character(len=96) :: piece
character(len=24) :: atoms_str
character(len=400) :: label, hdr

call find_scan(spec%opt_constraints, rest, kind, ia, kat, v0, v1, np, nd, found, out_unit)
if (.not. found) then
   call geomopt_run(spec, outfile, out_unit)
   return
endif
if (spec%opt_ts) then
   write(out_unit,'(A)') '[SCAN] ERROR: a scan cannot be combined with a transition-state search'
   return
endif
if (any(np(1:nd) < 1)) then
   write(out_unit,'(A)') '[SCAN] ERROR: the number of points must be at least 1'
   return
endif
ntot = product(np(1:nd))
nat = spec%ncenters
allocate(E(ntot), erel(ntot), ok(ntot), Vx(3,nat,ntot), tgt(nd,ntot), val(nd,ntot), ncyc(ntot), gidx(nd,ntot))
E = 0.0d0 ; ok = .false. ; tgt = 0.0d0 ; val = 0.0d0 ; ncyc = 0 ; gidx = 0
do j = 1, nd
   write(atoms_str,'(4(I0,1X))') ia(1:kat(j),j)
   names(j) = trim(kind(j))//' '//trim(atoms_str)
enddo
write(out_unit,'(A)') ''
if (nd == 1) then
   write(out_unit,'(A,I0,A,A,F10.4,A,F10.4)') '[SCAN] ', ntot, ' points: '//trim(names(1)), &
        '   from ', v0(1), '  to ', v1(1)
else
   write(out_unit,'(A,I0,A,I0,A)') '[SCAN] ', ntot, ' points on a ', nd, &
        '-dimensional grid (snake order, each point starts from the previous one):'
   do j = 1, nd
      write(out_unit,'(A,I0,A,A,F10.4,A,F10.4,A,I0,A)') '  q', j, ' = '//trim(names(j)), '   from ', v0(j), '  to ', v1(j), &
           '   ', np(j), ' points'
   enddo
endif
nfail = 0
npt = 0
restart0 = spec%opt_restart
do t = 0, ntot - 1
   do j = 1, nd
      base = 1
      do k = j + 1, nd
         base = base * np(k)
      enddo
      q = mod(t / base, np(j))
      par = mod(t / (base * np(j)), 2)
      d(j) = merge(np(j) - 1 - q, q, par == 1)
   enddo
   i = t + 1
   cons = rest
   label = ''
   do j = 1, nd
      if (np(j) > 1) then
         tgt(j,i) = v0(j) + (v1(j) - v0(j)) * dble(d(j)) / dble(np(j) - 1)
      else
         tgt(j,i) = v0(j)
      endif
      gidx(j,i) = d(j) + 1
      write(piece,'(F0.6)') tgt(j,i)
      write(atoms_str,'(4(I0,1X))') ia(1:kat(j),j)
      if (len_trim(cons) > 0) cons = trim(cons)//'; '
      cons = trim(cons)//trim(kind(j))//' '//trim(atoms_str)//' '//trim(piece)
      write(piece,'(F10.4)') tgt(j,i)
      label = trim(label)//'   '//trim(names(j))//' = '//trim(adjustl(piece))
   enddo
   spec%opt_constraints = cons
   spec%opt_restart = (t > 0)
   write(out_unit,'(A)') ''
   write(out_unit,'(A,I0,A,I0,A)') '[SCAN] point ', i, ' of ', ntot, trim(label)
   call geomopt_run(spec, outfile, out_unit)
   if (geomopt_last_ncyc == 0 .and. .not. geomopt_last_ok) then
      write(out_unit,'(A)') '[SCAN] the optimization could not start (see above); scan stopped'
      exit
   endif
   npt = i
   E(i) = geomopt_last_E ; ok(i) = geomopt_last_ok ; ncyc(i) = geomopt_last_ncyc
   if (.not. ok(i)) nfail = nfail + 1
   Vx(1,:,i) = spec%x(1:nat) ; Vx(2,:,i) = spec%y(1:nat) ; Vx(3,:,i) = spec%z(1:nat)
   do j = 1, nd
      val(j,i) = coord_value(kind(j), ia(:,j), kat(j), Vx(:,:,i), nat)
   enddo
enddo
spec%opt_constraints = trim(rest)
spec%opt_restart = restart0

if (npt < 1) return

scanfile = outfile
k = index(scanfile, '.', back=.true.)
if (k > 0) scanfile = scanfile(:k-1)
scanfile = trim(scanfile)//'_scan.xyz'
open(newunit=u, file=trim(scanfile), status='replace', action='write', iostat=ios)
if (ios == 0) then
   do i = 1, npt
      label = ''
      do j = 1, nd
         write(piece,'(F12.5)') val(j,i)
         if (j > 1) then
            label = trim(label)//'   '//trim(names(j))//' = '//trim(adjustl(piece))
         else
            label = trim(names(j))//' = '//trim(adjustl(piece))
         endif
      enddo
      write(u,'(I0)') nat
      write(u,'(A,F18.10,A)') trim(label)//'   E = ', E(i), merge('                 ', '  (not converged)', ok(i))
      do k = 1, nat
         write(u,'(A2,3F16.9)') elem_to_name(spec%atomchg(k)), Vx(:,k,i)
      enddo
   enddo
   close(u)
endif

eref = minval(E(1:npt))
if (any(ok(1:npt))) eref = minval(E(1:npt), mask=ok(1:npt))
erel(1:npt) = (E(1:npt) - eref) * HARTREE_TO_KCAL
write(out_unit,'(A)') ''
write(out_unit,'(A)') '[SCAN] result'
if (nd == 1) then
   hdr = '  point   '//trim(kind(1))//'         E (Ha)           E rel (kcal/mol)   converged'
   write(out_unit,'(A)') trim(hdr)
   do i = 1, npt
      write(out_unit,'(2X,I4,3X,F10.4,2X,F18.10,2X,F14.4,5X,A)') i, val(1,i), E(i), erel(i), merge('yes', 'NO ', ok(i))
   enddo
else
   do j = 1, nd
      write(out_unit,'(A,I0,A)') '  q', j, ' = '//trim(names(j))
   enddo
   hdr = '  point'
   do j = 1, nd
      write(piece,'(A,I0)') 'q', j
      hdr = trim(hdr)//'  '//repeat(' ', 8)//trim(piece)
   enddo
   write(out_unit,'(A)') trim(hdr)//'         E (Ha)           E rel (kcal/mol)   converged'
   do i = 1, npt
      write(label,'(2X,I5)') i
      do j = 1, nd
         write(piece,'(F11.4)') val(j,i)
         label = trim(label)//' '//piece(1:11)
      enddo
      write(piece,'(F18.10)') E(i)
      label = trim(label)//'  '//trim(piece)
      write(piece,'(F14.4)') erel(i)
      label = trim(label)//'  '//piece(1:14)//'     '//merge('yes', 'NO ', ok(i))
      write(out_unit,'(A)') trim(label)
   enddo
endif
if (nfail > 0) write(out_unit,'(A,I0,A)') '  ', nfail, ' point(s) did not converge (their energies are those of the last step)'
if (nd == 2) call write_matrix(out_unit, np(1), np(2), npt, gidx, erel, tgt, ok)
call json_scan(names, nd, npt, tgt, val, E, erel, ok, ncyc, Vx)
if (ios == 0) write(out_unit,'(A)') '  geometries: '//trim(scanfile)
end subroutine scan_or_opt_run

subroutine write_matrix(u, n1, n2, npt, gidx, erel, tgt, ok)
integer, intent(in) :: u, n1, n2, npt, gidx(:,:)
real(8), intent(in) :: erel(:), tgt(:,:)
logical, intent(in) :: ok(:)
real(8) :: m(n1,n2), c1(n1), c2(n2)
logical :: have(n1,n2)
integer :: i, a, b
character(len=2000) :: line
character(len=16) :: piece
m = 0.0d0 ; have = .false. ; c1 = 0.0d0 ; c2 = 0.0d0
do i = 1, npt
   a = gidx(1,i) ; b = gidx(2,i)
   m(a,b) = erel(i) ; have(a,b) = .true. ; c1(a) = tgt(1,i) ; c2(b) = tgt(2,i)
enddo
write(u,'(A)') ''
write(u,'(A)') '[SCAN] E rel (kcal/mol), rows q1, columns q2'
line = '            '
do b = 1, n2
   write(piece,'(F10.4)') c2(b)
   line = trim(line)//' '//piece(1:10)
enddo
write(u,'(A)') trim(line)
do a = 1, n1
   write(line,'(F10.4,2X)') c1(a)
   do b = 1, n2
      if (have(a,b)) then
         write(piece,'(F10.4)') m(a,b)
      else
         piece = '        --'
      endif
      line = trim(line)//' '//piece(1:10)
   enddo
   write(u,'(A)') trim(line)
enddo
end subroutine write_matrix

subroutine find_scan(cons, rest, kind, ia, kat, v0, v1, np, nd, found, out_unit)
character(len=*), intent(in) :: cons
character(len=*), intent(out) :: rest
character(len=*), intent(out) :: kind(:)
integer, intent(out) :: ia(:,:), kat(:), np(:), nd
real(8), intent(out) :: v0(:), v1(:)
logical, intent(out) :: found
integer, intent(in) :: out_unit
character(len=len(cons)) :: item
character(len=32) :: tok(10)
integer :: p, q, ntok, i, ios, k, ib(4)
real(8) :: x0, x1, rn

rest = '' ; found = .false. ; ia = 0 ; kat = 0 ; np = 0 ; v0 = 0.0d0 ; v1 = 0.0d0 ; kind = '' ; nd = 0
p = 1
do while (p <= len_trim(cons))
   q = index(cons(p:), ';')
   if (q == 0) then
      item = cons(p:)
      p = len_trim(cons) + 1
   else
      item = cons(p:p+q-2)
      p = p + q
   endif
   if (len_trim(item) == 0) cycle
   call split_words(item, tok, ntok)
   select case (lower(trim(tok(1))))
   case ('bond');     k = 2
   case ('angle');    k = 3
   case ('dihedral'); k = 4
   case default
      k = 0
   end select
   if (k > 0 .and. ntok == 1 + k + 3) then
      ib = 0
      ios = 0
      do i = 1, k
         read(tok(1+i), *, iostat=ios) ib(i)
         if (ios /= 0) exit
      enddo
      if (ios == 0) read(tok(2+k), *, iostat=ios) x0
      if (ios == 0) read(tok(3+k), *, iostat=ios) x1
      if (ios == 0) read(tok(4+k), *, iostat=ios) rn
      if (ios == 0) then
         if (nd >= size(kind)) then
            write(out_unit,'(A,I0,A)') '[SCAN] ERROR: at most ', size(kind), ' scan coordinates'
            stop 1
         endif
         nd = nd + 1
         found = .true.
         kind(nd) = lower(trim(tok(1)))
         ia(:,nd) = ib ; kat(nd) = k
         v0(nd) = x0 ; v1(nd) = x1 ; np(nd) = nint(rn)
         cycle
      endif
   endif
   if (len_trim(rest) > 0) rest = trim(rest)//'; '
   rest = trim(rest)//trim(adjustl(item))
enddo
end subroutine find_scan

real(8) function coord_value(kind, ia, k, xyz, nat)
character(len=*), intent(in) :: kind
integer, intent(in) :: ia(4), k, nat
real(8), intent(in) :: xyz(3,nat)
real(8) :: g(3,4)
integer :: j
do j = 1, k
   g(:,j) = xyz(:,ia(j)) * ANG2BOHR
enddo
select case (kind)
case ('bond');  coord_value = val_stretch(g(:,1), g(:,2)) / ANG2BOHR
case ('angle'); coord_value = val_bend(g(:,1), g(:,2), g(:,3)) * 180.0d0 / PI
case default;   coord_value = val_torsion(g(:,1), g(:,2), g(:,3), g(:,4)) * 180.0d0 / PI
end select
end function coord_value

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

end module mod_opt_scan
