! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Natural transition orbitals of TDDFT states.

module mod_tddft_nto
use MOL_info, only: nConts
use mod_molden, only: write_molden_mos
use mod_cube, only: cube_grid, cube_screen_begin, cube_screen_end, cube_ao_values, cube_chunk_points, cube_write_values
implicit none
private
public :: tddft_nto_rks, tddft_nto_uks

integer, parameter :: MAXPAIR = 6
real(8), parameter :: WMIN = 1.0d-3
real(8), parameter :: WCUBE = 1.0d-2

contains

subroutine tddft_nto_rks(nc, nocc, C, nroots, XpY, spec, stem)
integer,intent(in) :: nc, nocc, nroots
real(8),intent(in) :: C(nc,nc), XpY(nc-nocc,nocc,nroots)
character(len=*),intent(in) :: spec, stem
integer :: nst, st(64), is, s, nv, m, k, nk, no
logical :: cub(64)
real(8), allocatable :: U(:,:), VT(:,:), S1(:), A(:,:), coef(:,:), w(:), ene(:), occ(:)
character(len=8), allocatable :: lab(:), spn(:)
character(len=300) :: fn
character(len=2000) :: line
character(len=16) :: piece
call parse_nto(spec, nroots, nst, st, cub)
if (nst == 0) return
nv = nc - nocc ; no = nocc ; m = min(nv, no)
print '(A)', '[NTO]'
print '(A)', 'Natural transition orbitals: SVD of X+Y, weight = squared singular value (normalized to 1)'
do is = 1, nst
   s = st(is)
   allocate(A(nv,no), U(nv,m), VT(m,no), S1(m))
   A = XpY(:,:,s)
   call svd(A, nv, no, U, S1, VT, m)
   S1 = S1**2 / sum(S1**2)
   nk = npair(S1, m)
   allocate(coef(nc,2*nk), w(2*nk), ene(2*nk), occ(2*nk), lab(2*nk), spn(2*nk))
   do k = 1, nk
      call dgemv('N', nc, no, 1.0d0, C(1,1), nc, VT(k,1), m, 0.0d0, coef(1,k), 1)
      call dgemv('N', nc, nv, 1.0d0, C(1,nocc+1), nc, U(1,k), 1, 0.0d0, coef(1,nk+k), 1)
      write(lab(k),'(A,I0)') 'H', k
      write(lab(nk+k),'(A,I0)') 'P', k
      ene(k) = -S1(k) ; ene(nk+k) = S1(k)
      occ(k) = S1(k) ; occ(nk+k) = 0.0d0
      w(k) = S1(k) ; w(nk+k) = S1(k)
   enddo
   spn = 'Alpha'
   write(fn,'(A,I0,A)') trim(stem)//'.s', s, '.nto.molden'
   call write_molden_mos(trim(fn), 2*nk, coef, lab, ene, occ, spn)
   if (cub(is)) call write_cubes(stem, s, nc, 2*nk, coef, lab, w)
   line = ''
   do k = 1, nk
      write(piece,'(F8.4)') S1(k)
      line = trim(line)//'  '//trim(adjustl(piece))
   enddo
   print '(A,I0,A,A)', 'state ', s, '  weights of the NTO pairs:', trim(line)
   print '(A)', '  file: '//trim(fn)//'  (hole orbitals H1.., particle orbitals P1..)'
   deallocate(A, U, VT, S1, coef, w, ene, occ, lab, spn)
enddo
print '(A)', '[NTOEND]'
call flush(6)
end subroutine tddft_nto_rks

subroutine tddft_nto_uks(nc, noa, nob, Ca, Cb, nroots, XpY, spec, stem)
integer,intent(in) :: nc, noa, nob, nroots
real(8),intent(in) :: Ca(nc,nc), Cb(nc,nc), XpY(:,:)
character(len=*),intent(in) :: spec, stem
integer :: nst, st(64), is, s, nva, nvb, ma, mb, k, nka, nkb, nk, o
logical :: cub(64)
real(8), allocatable :: Ua(:,:), VTa(:,:), Sa(:), Ub(:,:), VTb(:,:), Sb(:), Aa(:,:), Ab(:,:), coef(:,:), ene(:), occ(:), wt(:)
character(len=8), allocatable :: lab(:), spn(:)
character(len=300) :: fn
character(len=2000) :: line
character(len=16) :: piece
real(8) :: tot
call parse_nto(spec, nroots, nst, st, cub)
if (nst == 0) return
nva = nc - noa ; nvb = nc - nob ; ma = min(nva, noa) ; mb = min(nvb, nob)
print '(A)', '[NTO]'
print '(A)', 'Natural transition orbitals: SVD of X+Y per spin, weight = squared singular value (normalized over both spins)'
do is = 1, nst
   s = st(is)
   allocate(Aa(nva,noa), Ab(nvb,nob), Ua(nva,ma), VTa(ma,noa), Sa(ma), Ub(nvb,mb), VTb(mb,nob), Sb(mb))
   Aa = reshape(XpY(1:nva*noa, s), [nva, noa])
   Ab = reshape(XpY(nva*noa+1:nva*noa+nvb*nob, s), [nvb, nob])
   if (ma > 0) call svd(Aa, nva, noa, Ua, Sa, VTa, ma)
   if (mb > 0) call svd(Ab, nvb, nob, Ub, Sb, VTb, mb)
   tot = sum(Sa**2) + sum(Sb**2)
   Sa = Sa**2 / tot ; Sb = Sb**2 / tot
   nka = npair(Sa, ma) ; nkb = npair(Sb, mb)
   if (ma == 0) nka = 0
   if (mb == 0) nkb = 0
   nk = nka + nkb
   allocate(coef(nc,2*nk), ene(2*nk), occ(2*nk), lab(2*nk), spn(2*nk), wt(2*nk))
   o = 0
   do k = 1, nka
      o = o + 1
      call dgemv('N', nc, noa, 1.0d0, Ca(1,1), nc, VTa(k,1), ma, 0.0d0, coef(1,o), 1)
      write(lab(o),'(A,I0)') 'Ha', k ; ene(o) = -Sa(k) ; occ(o) = Sa(k) ; spn(o) = 'Alpha' ; wt(o) = Sa(k)
      o = o + 1
      call dgemv('N', nc, nva, 1.0d0, Ca(1,noa+1), nc, Ua(1,k), 1, 0.0d0, coef(1,o), 1)
      write(lab(o),'(A,I0)') 'Pa', k ; ene(o) = Sa(k) ; occ(o) = 0.0d0 ; spn(o) = 'Alpha' ; wt(o) = Sa(k)
   enddo
   do k = 1, nkb
      o = o + 1
      call dgemv('N', nc, nob, 1.0d0, Cb(1,1), nc, VTb(k,1), mb, 0.0d0, coef(1,o), 1)
      write(lab(o),'(A,I0)') 'Hb', k ; ene(o) = -Sb(k) ; occ(o) = Sb(k) ; spn(o) = 'Beta ' ; wt(o) = Sb(k)
      o = o + 1
      call dgemv('N', nc, nvb, 1.0d0, Cb(1,nob+1), nc, Ub(1,k), 1, 0.0d0, coef(1,o), 1)
      write(lab(o),'(A,I0)') 'Pb', k ; ene(o) = Sb(k) ; occ(o) = 0.0d0 ; spn(o) = 'Beta ' ; wt(o) = Sb(k)
   enddo
   write(fn,'(A,I0,A)') trim(stem)//'.s', s, '.nto.molden'
   call write_molden_mos(trim(fn), 2*nk, coef, lab, ene, occ, spn)
   if (cub(is)) call write_cubes(stem, s, nc, 2*nk, coef, lab, wt)
   line = ''
   do k = 1, nka
      write(piece,'(F8.4)') Sa(k)
      line = trim(line)//'  '//trim(adjustl(piece))
   enddo
   print '(A,I0,A,A)', 'state ', s, '  alpha weights:', trim(line)
   line = ''
   do k = 1, nkb
      write(piece,'(F8.4)') Sb(k)
      line = trim(line)//'  '//trim(adjustl(piece))
   enddo
   print '(A,A)', '          beta weights: ', trim(line)
   print '(A)', '  file: '//trim(fn)//'  (hole orbitals Ha/Hb, particle orbitals Pa/Pb)'
   deallocate(Aa, Ab, Ua, VTa, Sa, Ub, VTb, Sb, coef, ene, occ, lab, spn, wt)
enddo
print '(A)', '[NTOEND]'
call flush(6)
end subroutine tddft_nto_uks

subroutine write_cubes(stem, state, nc, nmo, coef, lab, w)
character(len=*),intent(in) :: stem
integer,intent(in) :: state, nc, nmo
real(8),intent(in) :: coef(nc,nmo), w(nmo)
character(len=*),intent(in) :: lab(nmo)
integer,parameter :: NB = 2000
integer :: sel(nmo), nsel, o, n(3), ntot, ip0, nb_now, ix
real(8) :: org(3), h
real(8),allocatable :: csel(:,:), vals(:,:), Phi(:,:), xyz(:,:)
logical :: own
character(len=300) :: fname
nsel = 0
do o = 1, nmo
   if (w(o) >= WCUBE .or. lab(o)(2:) == '1' .or. lab(o)(3:) == '1') then
      nsel = nsel + 1 ; sel(nsel) = o
   endif
enddo
if (nsel == 0) return
call cube_grid(6.0d0, 0.2d0, 2000000, n, org, h)
ntot = n(1)*n(2)*n(3)
allocate(csel(nc,nsel), vals(ntot,nsel), Phi(NB,nc), xyz(3,NB))
do o = 1, nsel
   csel(:,o) = coef(:,sel(o))
enddo
call cube_screen_begin(own)
do ip0 = 1, ntot, NB
   nb_now = min(NB, ntot - ip0 + 1)
   call cube_chunk_points(ip0, nb_now, n, org, h, xyz)
   call cube_ao_values(nb_now, xyz, Phi)
   call dgemm('N', 'N', nb_now, nsel, nc, 1.0d0, Phi, NB, csel, nc, 0.0d0, vals(ip0,1), ntot)
enddo
call cube_screen_end(own)
do o = 1, nsel
   write(fname,'(A,I0,A)') trim(stem)//'.s', state, '.nto.'//trim(lower(lab(sel(o))))//'.cube'
   call cube_write_values(trim(fname), 'natural transition orbital '//trim(lab(sel(o)))//', state '//trim(itoa(state)), &
        vals(:,o), n, org, h)
   print '(A,A)', '  cube file written: ', trim(fname)
enddo
deallocate(csel, vals, Phi, xyz)
end subroutine write_cubes

function lower(t) result(u)
character(len=*),intent(in) :: t
character(len=len(t)) :: u
integer :: i, c
u = t
do i = 1, len(t)
   c = iachar(t(i:i))
   if (c >= 65 .and. c <= 90) u(i:i) = achar(c + 32)
enddo
end function lower

function itoa(i) result(t)
integer,intent(in) :: i
character(len=12) :: t
write(t,'(I0)') i
end function itoa

integer function npair(w, m)
real(8),intent(in) :: w(:)
integer,intent(in) :: m
integer :: k
npair = 0
do k = 1, min(m, MAXPAIR)
   if (w(k) >= WMIN) npair = k
enddo
if (m > 0) npair = max(npair, 1)
end function npair

subroutine svd(A, n, p, U, S, VT, m)
integer,intent(in) :: n, p, m
real(8),intent(inout) :: A(n,p)
real(8),intent(out) :: U(n,m), S(m), VT(m,p)
real(8), allocatable :: work(:)
real(8) :: wq(1)
integer :: info, lwork
call dgesvd('S', 'S', n, p, A, n, S, U, n, VT, m, wq, -1, info)
lwork = max(1, int(wq(1)))
allocate(work(lwork))
call dgesvd('S', 'S', n, p, A, n, S, U, n, VT, m, work, lwork, info)
if (info /= 0) error stop 'tddft nto: singular value decomposition failed'
end subroutine svd

subroutine parse_nto(spec, nroots, nst, st, cub)
character(len=*),intent(in) :: spec
integer,intent(in) :: nroots
integer,intent(out) :: nst, st(64)
logical,intent(out) :: cub(64)
character(len=len(spec)) :: low, tok
integer :: i, j, ic, ios, k
nst = 0
low = spec
do i = 1, len(low)
   if (low(i:i) >= 'A' .and. low(i:i) <= 'Z') low(i:i) = char(ichar(low(i:i)) + 32)
enddo
i = 1
do while (i <= len_trim(low))
   j = index(low(i:), ',')
   if (j == 0) then
      tok = adjustl(low(i:)); i = len_trim(low) + 1
   else
      tok = adjustl(low(i:i+j-2)); i = i + j
   endif
   ic = index(tok, ':')
   if (ic == 0) cycle
   if (trim(tok(1:ic-1)) /= 'nto' .and. trim(tok(1:ic-1)) /= 'ntom') cycle
   read(tok(ic+1:), *, iostat=ios) k
   if (ios /= 0 .or. k < 1 .or. k > nroots) error stop 'nto: the state is not among the computed states'
   if (any(st(1:nst) == k)) cycle
   nst = nst + 1
   if (nst > 64) error stop 'nto: at most 64 states'
   st(nst) = k
   cub(nst) = (trim(tok(1:ic-1)) == 'nto')
enddo
end subroutine parse_nto

end module mod_tddft_nto
