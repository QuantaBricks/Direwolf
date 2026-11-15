! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Hole and electron densities of the TDDFT states, their charge-transfer measures and cube files.

module mod_tddft_density
use MOL_info, only: Natoms, atoms
use mod_cube, only: cube_grid, cube_screen_begin, cube_screen_end, cube_ao_values, cube_chunk_points, cube_write_values
implicit none
private
public :: tddft_hole_electron, tddft_hole_electron_uks, tddft_he_wanted

contains

function tddft_he_wanted(spec) result(want)
character(len=*),intent(in) :: spec
logical :: want
character(len=len(spec)) :: low, tok
integer :: i, j, ic
want = .false.
low = spec
do i = 1,len(low)
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
   if (trim(tok) == 'ct') then
      want = .true.
   else
      ic = index(tok, ':')
      if (ic > 0) then
         select case (trim(tok(1:ic-1)))
         case ('hole', 'elec', 'diff', 'trans'); want = .true.
         end select
      endif
   endif
enddo
end function tddft_he_wanted

subroutine tddft_hole_electron(nConts, nocc, C_a, nroots, XpY, XmY, spec, stem, ct)
implicit none
integer,intent(in) :: nConts, nocc, nroots
real(8),intent(in) :: C_a(nConts,nConts), XpY(nConts-nocc,nocc,nroots), XmY(nConts-nocc,nocc,nroots)
character(len=*),intent(in) :: spec, stem
real(8),intent(out) :: ct(5,nroots)
real(8) :: z(1,1,1)
call hole_electron_core(nConts, 1, nocc, 0, C_a, C_a, XpY, XmY, z, z, nroots, sqrt(2.0d0), spec, stem, ct)
end subroutine tddft_hole_electron

subroutine tddft_hole_electron_uks(nConts, noa, nob, Ca, Cb, nroots, XpY, XmY, spec, stem, ct)
implicit none
integer,intent(in) :: nConts, noa, nob, nroots
real(8),intent(in) :: Ca(nConts,nConts), Cb(nConts,nConts)
real(8),intent(in) :: XpY((nConts-noa)*noa+(nConts-nob)*nob,nroots), XmY((nConts-noa)*noa+(nConts-nob)*nob,nroots)
character(len=*),intent(in) :: spec, stem
real(8),intent(out) :: ct(5,nroots)
integer :: nva, nvb, na, r
real(8),allocatable :: xa(:,:,:), ma(:,:,:), xb(:,:,:), mb(:,:,:)
nva = nConts - noa ; nvb = nConts - nob ; na = nva*noa
allocate(xa(nva,noa,nroots), ma(nva,noa,nroots), xb(max(nvb,1),max(nob,1),nroots), mb(max(nvb,1),max(nob,1),nroots))
xb = 0.0d0 ; mb = 0.0d0
do r = 1,nroots
   xa(:,:,r) = reshape(XpY(1:na,r), [nva,noa])
   ma(:,:,r) = reshape(XmY(1:na,r), [nva,noa])
   if (nob > 0) then
      xb(1:nvb,1:nob,r) = reshape(XpY(na+1:na+nvb*nob,r), [nvb,nob])
      mb(1:nvb,1:nob,r) = reshape(XmY(na+1:na+nvb*nob,r), [nvb,nob])
   endif
enddo
call hole_electron_core(nConts, 2, noa, nob, Ca, Cb, xa, ma, xb, mb, nroots, 1.0d0, spec, stem, ct)
end subroutine tddft_hole_electron_uks

subroutine hole_electron_core(nConts, nspin, nocc1, nocc2, C1, C2, X1, M1, X2, M2, nroots, tfac, spec, stem, ct)
use omp_lib
use MOL_info, only: Natoms, atoms
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: nConts, nspin, nocc1, nocc2, nroots
real(8),intent(in) :: C1(nConts,nConts), C2(nConts,nConts), tfac
real(8),intent(in) :: X1(:,:,:), M1(:,:,:), X2(:,:,:), M2(:,:,:)
character(len=*),intent(in) :: spec, stem
real(8),intent(out) :: ct(5,nroots)

integer,parameter :: NB = 2048, GRP = 5
integer :: nreq, ireq, r, r0, r1, ng, k, ip0, nb_now, n(3), ntot, nkeep, ia, o, ix, iy, iz, ip, n23, sp
integer :: nocc, nvir, nom
real(8) :: org(3), h, w, rc2, d2
character(len=8) :: kind_req(32)
integer :: state_req(32)
logical :: own, any_req
integer,allocatable :: keep(:)
real(8),allocatable :: Phi(:,:), psi(:,:), chi(:,:), AB(:,:), PH(:,:), PS(:,:), xyz(:,:), XS(:,:)
real(8),allocatable :: rh(:,:), re(:,:), rt(:,:), vals(:,:), rhs(:), res(:), rts(:)
real(8),allocatable :: Nh(:), Ne(:), mh(:,:), me(:,:), qp(:), ov(:), Cs(:,:), xw(:,:), mw(:,:)
real(8) :: x(3), dh, de, cen_h(3), cen_e(3), sig2

call parse_spec(spec, nroots, nreq, kind_req, state_req)
any_req = nreq > 0
if (any_req) then
   call cube_grid(6.0d0, 0.2d0, 3000000, n, org, h)
else
   call cube_grid(6.0d0, hmin_metrics(), ngrid_metrics(), n, org, h)
endif
ntot = n(1)*n(2)*n(3)
n23 = n(2)*n(3)
w = h**3
rc2 = 36.0d0
allocate(keep(ntot))
nkeep = 0
!$omp parallel do private(ip, ix, iy, iz, x, d2, ia) schedule(static) reduction(+:nkeep)
do ip = 1,ntot
   iz = mod(ip-1, n(3)); iy = mod((ip-1)/n(3), n(2)); ix = (ip-1)/n23
   x = org + h*[real(ix,8), real(iy,8), real(iz,8)]
   keep(ip) = 0
   do ia = 1,Natoms
      d2 = sum((x - atoms(ia)%coor*ans2bohr)**2)
      if (d2 < rc2) then
         keep(ip) = 1
         nkeep = nkeep + 1
         exit
      endif
   enddo
enddo
!$omp end parallel do
call compact_index(ntot, keep, nkeep)
allocate(Nh(nroots), Ne(nroots), mh(4,nroots), me(4,nroots), qp(nroots), ov(nroots))
Nh = 0.0d0; Ne = 0.0d0; mh = 0.0d0; me = 0.0d0; qp = 0.0d0; ov = 0.0d0
if (any_req) then
   allocate(vals(ntot,nreq)); vals = 0.0d0
endif
nom = max(nocc1, nocc2)
allocate(Phi(NB,nConts), psi(NB,nom), chi(NB,nConts), xyz(3,NB), Cs(nConts,nConts), xw(nConts,nom), mw(nConts,nom))
allocate(AB(NB,2*nom*min(GRP,nroots)), PH(nom,nom*min(GRP,nroots)), PS(NB,nom))
allocate(XS(nConts,2*nom*min(GRP,nroots)))
allocate(rh(NB,GRP), re(NB,GRP), rt(NB,GRP), rhs(NB), res(NB), rts(NB))
call cube_screen_begin(own)

do r0 = 1,nroots,GRP
   r1 = min(r0 + GRP - 1, nroots)
   ng = r1 - r0 + 1
   do ip0 = 1,nkeep,NB
      nb_now = min(NB, nkeep - ip0 + 1)
      do k = 1,nb_now
         ip = keep(ip0 + k - 1)
         iz = mod(ip-1, n(3)); iy = mod((ip-1)/n(3), n(2)); ix = (ip-1)/n23
         xyz(:,k) = org + h*[real(ix,8), real(iy,8), real(iz,8)]
      enddo
      call cube_ao_values(nb_now, xyz, Phi)
      rh = 0.0d0; re = 0.0d0; rt = 0.0d0
      do sp = 1,nspin
         if (sp == 1) then
            nocc = nocc1; Cs = C1
         else
            nocc = nocc2; Cs = C2
         endif
         if (nocc == 0) cycle
         nvir = nConts - nocc
         call dgemm('N','N',nb_now,nocc,nConts,1.0d0,Phi,NB,Cs,nConts,0.0d0,psi,NB)
         call dgemm('N','N',nb_now,nvir,nConts,1.0d0,Phi,NB,Cs(1,nocc+1),nConts,0.0d0,chi,NB)
         do r = r0,r1
            o = (r - r0)*2*nocc
            if (sp == 1) then
               XS(1:nvir,o+1:o+nocc) = X1(1:nvir,1:nocc,r)
               XS(1:nvir,o+nocc+1:o+2*nocc) = M1(1:nvir,1:nocc,r)
            else
               XS(1:nvir,o+1:o+nocc) = X2(1:nvir,1:nocc,r)
               XS(1:nvir,o+nocc+1:o+2*nocc) = M2(1:nvir,1:nocc,r)
            endif
            o = (r - r0)*nocc
            if (sp == 1) then
               xw(1:nvir,1:nocc) = X1(1:nvir,1:nocc,r); mw(1:nvir,1:nocc) = M1(1:nvir,1:nocc,r)
            else
               xw(1:nvir,1:nocc) = X2(1:nvir,1:nocc,r); mw(1:nvir,1:nocc) = M2(1:nvir,1:nocc,r)
            endif
            call dgemm('T','N',nocc,nocc,nvir,0.5d0,xw,nConts,xw,nConts,0.0d0,PH(1,o+1),nom)
            call dgemm('T','N',nocc,nocc,nvir,0.5d0,mw,nConts,mw,nConts,1.0d0,PH(1,o+1),nom)
         enddo
         call dgemm('N','N',nb_now,2*nocc*ng,nvir,1.0d0,chi,NB,XS,nConts,0.0d0,AB,NB)
         do r = r0,r1
            o = (r - r0)*2*nocc
            call dgemm('N','N',nb_now,nocc,nocc,1.0d0,psi,NB,PH(1,(r-r0)*nocc+1),nom,0.0d0,PS,NB)
            do k = 1,nb_now
               res(k) = 0.5d0*sum(AB(k,o+1:o+2*nocc)**2)
               rhs(k) = sum(PS(k,1:nocc)*psi(k,1:nocc))
               rts(k) = tfac*sum(psi(k,1:nocc)*AB(k,o+1:o+nocc))
            enddo
            rh(1:nb_now,r-r0+1) = rh(1:nb_now,r-r0+1) + rhs(1:nb_now)
            re(1:nb_now,r-r0+1) = re(1:nb_now,r-r0+1) + res(1:nb_now)
            rt(1:nb_now,r-r0+1) = rt(1:nb_now,r-r0+1) + rts(1:nb_now)
         enddo
      enddo
      do r = r0,r1
         call add_state(r, r - r0 + 1, nb_now)
      enddo
   enddo
enddo
call cube_screen_end(own)

do r = 1,nroots
   dh = Nh(r)*w; de = Ne(r)*w
   cen_h = mh(1:3,r)/Nh(r); cen_e = me(1:3,r)/Ne(r)
   ct(1,r) = sqrt(sum((cen_e - cen_h)**2))/ans2bohr
   ct(2,r) = qp(r)*w
   sig2 = mh(4,r)/Nh(r) - sum(cen_h**2); ct(3,r) = sqrt(max(sig2, 0.0d0))/ans2bohr
   sig2 = me(4,r)/Ne(r) - sum(cen_e**2); ct(4,r) = sqrt(max(sig2, 0.0d0))/ans2bohr
   ct(5,r) = ov(r)*w/sqrt(dh*de)
enddo

do ireq = 1,nreq
   block
      character(len=256) :: fname
      character(len=64) :: what
      write(fname,'(A,A,I0,A,A,A)') trim(stem), '.s', state_req(ireq), '.', trim(kind_req(ireq)), '.cube'
      select case (trim(kind_req(ireq)))
      case ('hole');  what = 'hole density (e/bohr^3)'
      case ('elec');  what = 'electron density (e/bohr^3)'
      case ('diff');  what = 'electron minus hole density (e/bohr^3)'
      case default;   what = 'transition density (1/bohr^3)'
      end select
      call cube_write_values(trim(fname), trim(what)//', state '//itoa(state_req(ireq)), vals(:,ireq), n, org, h)
      print '(A,A)', '  cube file written: ', trim(fname)
   end block
enddo
deallocate(Phi, psi, chi, AB, PH, PS, XS, xyz, rh, re, rt, rhs, res, rts, keep, Cs, xw, mw)
if (allocated(vals)) deallocate(vals)

contains
   function itoa(m) result(t)
   integer,intent(in) :: m
   character(len=12) :: t
   write(t,'(I0)') m
   end function itoa

   subroutine add_state(rr, cc, nbn)
   integer,intent(in) :: rr, cc, nbn
   integer :: kk, jj
   real(8) :: xx(3)
   do kk = 1,nbn
      xx = xyz(:,kk)
      Nh(rr) = Nh(rr) + rh(kk,cc); Ne(rr) = Ne(rr) + re(kk,cc)
      mh(1:3,rr) = mh(1:3,rr) + rh(kk,cc)*xx; mh(4,rr) = mh(4,rr) + rh(kk,cc)*sum(xx*xx)
      me(1:3,rr) = me(1:3,rr) + re(kk,cc)*xx; me(4,rr) = me(4,rr) + re(kk,cc)*sum(xx*xx)
      qp(rr) = qp(rr) + max(re(kk,cc) - rh(kk,cc), 0.0d0)
      ov(rr) = ov(rr) + sqrt(rh(kk,cc)*re(kk,cc))
   enddo
   do ireq = 1,nreq
      if (state_req(ireq) /= rr) cycle
      do kk = 1,nbn
         jj = keep(ip0 + kk - 1)
         select case (trim(kind_req(ireq)))
         case ('hole');  vals(jj,ireq) = rh(kk,cc)
         case ('elec');  vals(jj,ireq) = re(kk,cc)
         case ('diff');  vals(jj,ireq) = re(kk,cc) - rh(kk,cc)
         case ('trans'); vals(jj,ireq) = rt(kk,cc)
         end select
      enddo
   enddo
   end subroutine add_state
end subroutine hole_electron_core

function hmin_metrics() result(x)
real(8) :: x
integer :: ios
character(len=32) :: e
x = 0.35d0
call get_environment_variable("ENGINE_TDDFT_CT_H", e)
if (len_trim(e) > 0) then
   read(e, *, iostat=ios) x
   if (ios /= 0 .or. x < 0.05d0) x = 0.35d0
endif
end function hmin_metrics

function ngrid_metrics() result(m)
integer :: m, ios
character(len=32) :: e
m = 600000
call get_environment_variable("ENGINE_TDDFT_CT_POINTS", e)
if (len_trim(e) > 0) then
   read(e, *, iostat=ios) m
   if (ios /= 0 .or. m < 1000) m = 600000
endif
end function ngrid_metrics

subroutine compact_index(ntot, keep, nkeep)
integer,intent(in) :: ntot, nkeep
integer,intent(inout) :: keep(ntot)
integer :: ip, m
m = 0
do ip = 1,ntot
   if (keep(ip) == 1) then
      m = m + 1
      keep(m) = ip
   endif
enddo
end subroutine compact_index

subroutine parse_spec(spec, nroots, nreq, kinds, states)
character(len=*),intent(in) :: spec
integer,intent(in) :: nroots
integer,intent(out) :: nreq, states(32)
character(len=8),intent(out) :: kinds(32)
character(len=len(spec)) :: low, tok
integer :: i, j, p, ic, ios, st
nreq = 0
low = spec
do i = 1,len(low)
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
   if (len_trim(tok) == 0) cycle
   if (trim(tok) == 'ct') cycle
   ic = index(tok, ':')
   if (ic == 0) error stop 'tddft_cube: write kind:state, e.g. diff:1, hole:1, elec:1, trans:2'
   read(tok(ic+1:), *, iostat=ios) st
   if (ios /= 0 .or. st < 1 .or. st > nroots) error stop 'tddft_cube: the state is not among the computed states'
   select case (trim(tok(1:ic-1)))
   case ('hole', 'elec', 'diff', 'trans')
   case ('nto', 'ntom')
      cycle
   case default
      error stop 'tddft_cube: kind must be hole, elec, diff, trans or nto'
   end select
   nreq = nreq + 1
   if (nreq > 32) error stop 'tddft_cube: at most 32 files'
   kinds(nreq) = tok(1:ic-1)
   states(nreq) = st
enddo
end subroutine parse_spec

end module mod_tddft_density
