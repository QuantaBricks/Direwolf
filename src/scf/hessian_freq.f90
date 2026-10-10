! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Harmonic frequencies from a Cartesian Hessian: mass-weight, project out

module mod_hessian_freq
implicit none
private
public :: hessian_write_full, hessian_frequencies, atomic_mass, hess_classical, &
          print_normal_modes, ir_intensities, print_ir, write_molden_freq, write_modes_csv, freq_unit, &
          print_thermo, freq_energy, freq_mult, freq_temp, freq_pres, freq_base, freq_thermo

real(8),allocatable,save :: hess_classical(:,:,:,:)
integer,save :: freq_unit = 6
character(len=256),save :: freq_base = 'direwolf'
logical,save :: freq_thermo = .true.
real(8),save :: freq_energy = 0d0
integer,save :: freq_mult = 1
real(8),save :: freq_temp = 298.15d0
real(8),save :: freq_pres = 101325d0

contains

function atomic_mass(z) result(m)
implicit none
integer,intent(in) :: z
real(8) :: m
real(8),parameter :: mass_table(86) = (/ &
   1.008d0,      4.002602d0,   6.94d0,       9.0121831d0,  10.81d0,      &
   12.011d0,     14.007d0,     15.999d0,     18.998403163d0, 20.1797d0,  &
   22.98976928d0, 24.305d0,    26.9815385d0, 28.085d0,     30.973761998d0, &
   32.06d0,      35.45d0,      39.948d0,     39.0983d0,    40.078d0,     &
   44.955908d0,  47.867d0,     50.9415d0,    51.9961d0,    54.938044d0,  &
   55.845d0,     58.933194d0,  58.6934d0,    63.546d0,     65.38d0,      &
   69.723d0,     72.630d0,     74.921595d0,  78.971d0,     79.904d0,     &
   83.798d0,     85.4678d0,    87.62d0,      88.90584d0,   91.224d0,     &
   92.90637d0,   95.95d0,      97.90721d0,   101.07d0,     102.90550d0,  &
   106.42d0,     107.8682d0,   112.414d0,    114.818d0,    118.710d0,    &
   121.760d0,    127.60d0,     126.90447d0,  131.293d0,    132.90545196d0, &
   137.327d0,    138.90547d0,  140.116d0,    140.90766d0,  144.242d0,    &
   144.91276d0,  150.36d0,     151.964d0,    157.25d0,     158.92535d0,  &
   162.500d0,    164.93033d0,  167.259d0,    168.93422d0,  173.045d0,    &
   174.9668d0,   178.49d0,     180.94788d0,  183.84d0,     186.207d0,    &
   190.23d0,     192.217d0,    195.084d0,    196.966569d0, 200.592d0,    &
   204.38d0,     207.2d0,      208.98040d0,  208.98243d0,  209.98715d0,  &
   222.01758d0 /)
if (z .ge. 1 .and. z .le. 86) then
   m = mass_table(z)
else
   m = 2.5d0*real(z,8)
endif
end function atomic_mass

subroutine hessian_write_full(unit_out, natoms, d2E)
implicit none
integer,intent(in) :: unit_out, natoms
real(8),intent(in) :: d2E(3,natoms,3,natoms)
integer :: n, i, j, ia, ic, ja, jc, cstart, cend, c
real(8),allocatable :: H(:,:)

n = 3*natoms
allocate(H(n,n))
do ja = 1,natoms
   do jc = 1,3
      j = 3*(ja-1)+jc
      do ia = 1,natoms
         do ic = 1,3
            i = 3*(ia-1)+ic
            H(i,j) = d2E(ic,ia,jc,ja)
         enddo
      enddo
   enddo
enddo

write(unit_out,'(A)') ''
write(unit_out,'(A)') '$hessian'
write(unit_out,'(I0)') n
cstart = 1
do while (cstart .le. n)
   cend = min(n, cstart+4)
   write(unit_out,'(5(I18))') (c-1, c = cstart,cend)
   do i = 1,n
      write(unit_out,'(I6,5(1X,ES18.10E2))') i-1, (H(i,c), c = cstart,cend)
   enddo
   cstart = cend + 1
enddo
deallocate(H)
end subroutine hessian_write_full

subroutine hessian_frequencies(natoms, znum, coor_bohr, d2E, freq_cm, nfreq, disp, nzero)
implicit none
integer,intent(in) :: natoms, znum(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(in) :: d2E(3,natoms,3,natoms)
real(8),intent(out) :: freq_cm(3*natoms)
integer,intent(out) :: nfreq
real(8),intent(out),optional :: disp(3*natoms,3*natoms)
integer,intent(out),optional :: nzero

real(8),parameter :: AU2CM = 5140.4981d0
real(8),parameter :: PROJ_TOL = 1.0d-8

integer :: n, i, j, ia, ic, ja, jc, k, nproj, info, lwork
real(8),allocatable :: Hm(:,:), P(:,:), V(:,:), eval(:), work(:)
real(8),allocatable :: mass(:), sqm(:)
real(8) :: com(3), totm, r(3), nrm, dotv
real(8) :: rotv(3*natoms,3), trav(3*natoms,3)
real(8),allocatable :: basis(:,:)

n = 3*natoms
allocate(mass(natoms), sqm(n))
totm = 0.0d0
com = 0.0d0
do ia = 1,natoms
   mass(ia) = atomic_mass(znum(ia))
   totm = totm + mass(ia)
   com = com + mass(ia)*coor_bohr(:,ia)
enddo
com = com/totm
do ia = 1,natoms
   do ic = 1,3
      sqm(3*(ia-1)+ic) = sqrt(mass(ia))
   enddo
enddo

allocate(Hm(n,n))
do ja = 1,natoms
   do jc = 1,3
      j = 3*(ja-1)+jc
      do ia = 1,natoms
         do ic = 1,3
            i = 3*(ia-1)+ic
            Hm(i,j) = d2E(ic,ia,jc,ja)/(sqm(i)*sqm(j))
         enddo
      enddo
   enddo
enddo

trav = 0.0d0
rotv = 0.0d0
do ia = 1,natoms
   r = coor_bohr(:,ia) - com
   do ic = 1,3
      trav(3*(ia-1)+ic,ic) = sqm(3*(ia-1)+ic)
   enddo
   rotv(3*(ia-1)+2,1) = -r(3)*sqm(3*(ia-1)+2)
   rotv(3*(ia-1)+3,1) =  r(2)*sqm(3*(ia-1)+3)
   rotv(3*(ia-1)+1,2) =  r(3)*sqm(3*(ia-1)+1)
   rotv(3*(ia-1)+3,2) = -r(1)*sqm(3*(ia-1)+3)
   rotv(3*(ia-1)+1,3) = -r(2)*sqm(3*(ia-1)+1)
   rotv(3*(ia-1)+2,3) =  r(1)*sqm(3*(ia-1)+2)
enddo

allocate(basis(n,6))
nproj = 0
do k = 1,6
   if (k .le. 3) then
      basis(:,nproj+1) = trav(:,k)
   else
      basis(:,nproj+1) = rotv(:,k-3)
   endif
   do j = 1,nproj
      dotv = dot_product(basis(:,nproj+1), basis(:,j))
      basis(:,nproj+1) = basis(:,nproj+1) - dotv*basis(:,j)
   enddo
   nrm = sqrt(dot_product(basis(:,nproj+1), basis(:,nproj+1)))
   if (nrm .gt. PROJ_TOL) then
      basis(:,nproj+1) = basis(:,nproj+1)/nrm
      nproj = nproj + 1
   endif
enddo

allocate(P(n,n))
P = 0.0d0
do i = 1,n
   P(i,i) = 1.0d0
enddo
do k = 1,nproj
   do j = 1,n
      do i = 1,n
         P(i,j) = P(i,j) - basis(i,k)*basis(j,k)
      enddo
   enddo
enddo
allocate(V(n,n))
V = matmul(P, matmul(Hm, P))

allocate(eval(n))
lwork = 1 + 6*n + 2*n*n
allocate(work(lwork))
block
   integer,allocatable :: iwork(:)
   integer :: liwork
   liwork = 3 + 5*n
   allocate(iwork(liwork))
   call dsyevd('V', 'U', n, V, n, eval, work, lwork, iwork, liwork, info)
   deallocate(iwork)
end block

block
   integer :: perm(n), nz, m
   real(8) :: wgt
   nz = 0
   do k = 1,n
      wgt = 0.0d0
      do j = 1,nproj
         wgt = wgt + dot_product(basis(:,j), V(:,k))**2
      enddo
      if (wgt .gt. 0.5d0) then
         nz = nz + 1
         perm(nz) = k
      endif
   enddo
   m = nz
   do k = 1,n
      wgt = 0.0d0
      do j = 1,nproj
         wgt = wgt + dot_product(basis(:,j), V(:,k))**2
      enddo
      if (wgt .le. 0.5d0) then
         m = m + 1
         perm(m) = k
      endif
   enddo
   if (present(nzero)) nzero = nz
   nfreq = n
   do i = 1,n
      k = perm(i)
      if (eval(k) .ge. 0.0d0) then
         freq_cm(i) =  AU2CM*sqrt(eval(k))
      else
         freq_cm(i) = -AU2CM*sqrt(-eval(k))
      endif
      if (present(disp)) disp(:,i) = V(:,k)/sqm
   enddo
end block

deallocate(Hm,P,V,eval,work,mass,sqm,basis)
end subroutine hessian_frequencies

subroutine print_normal_modes(natoms, znum, freq_cm, disp, nskip, redmass, inten)
implicit none
integer,intent(in) :: natoms, znum(natoms), nskip
real(8),intent(in) :: freq_cm(3*natoms), disp(3*natoms,3*natoms)
real(8),intent(out) :: redmass(3*natoms)
real(8),intent(in) :: inten(3*natoms)
integer :: n, k, ia, c0, k0, k1, nb, j
character(len=26) :: lab
n = 3*natoms
redmass = 0.0d0
do k = nskip+1,n
   redmass(k) = 1.0d0/sum(disp(:,k)**2)
enddo
write(freq_unit,'(A)') ' [NORMALMODES]'
write(freq_unit,'(A)') ' Cartesian displacements of the vibrations (unit norm), three modes per group.'
write(freq_unit,'(A)') ' Imaginary frequencies are negative; translations and rotations are not listed.'
do k0 = nskip+1,n,3
   k1 = min(k0+2,n)
   write(freq_unit,'(A)') ''
   write(freq_unit,'(A26,3I32)') ' Mode', (k-nskip, k = k0,k1)
   write(freq_unit,'(A26,3F32.4)') ' Frequency (cm-1)', freq_cm(k0:k1)
   write(freq_unit,'(A26,3F32.4)') ' IR intensity (km/mol)', inten(k0:k1)
   write(freq_unit,'(A26,3F32.4)') ' Reduced mass (amu)', redmass(k0:k1)
   write(freq_unit,'(A26,3F32.4)') ' Force constant (mDyne/A)', (fc(k), k = k0,k1)
   write(freq_unit,'(A26,3(2X,3A10))') ' Atom', ('        dx', '        dy', '        dz', k = k0,k1)
   do ia = 1,natoms
      lab = ''
      write(lab,'(1X,I4,2X,A)') ia, trim(elem_sym(znum(ia)))
      c0 = 3*(ia-1)
      write(freq_unit,'(A26,3(2X,3F10.4))') lab, ((disp(c0+j,k)*sqrt(redmass(k)), j = 1,3), k = k0,k1)
   enddo
enddo
write(freq_unit,'(A)') ' [NORMALMODESEND]'
contains
real(8) function fc(k)
integer,intent(in) :: k
real(8),parameter :: PI = 3.14159265358979d0, CL = 2.99792458d8, AMU = 1.66053906660d-27
fc = 4.0d0*PI**2*CL**2*(freq_cm(k)*100.0d0)**2*redmass(k)*AMU/100.0d0
if (freq_cm(k) < 0.0d0) fc = -fc
end function fc
end subroutine print_normal_modes

subroutine ir_intensities(natoms, dmu, disp, inten)
implicit none
integer,intent(in) :: natoms
real(8),intent(in) :: dmu(3,3*natoms), disp(3*natoms,3*natoms)
real(8),intent(out) :: inten(3*natoms)
real(8),parameter :: NA = 6.02214076d23, EPS0 = 8.8541878128d-12, CL = 2.99792458d8
real(8),parameter :: ECH = 1.602176634d-19, AMU = 1.66053906660d-27
real(8) :: fac, dq(3)
integer :: k, n
n = 3*natoms
fac = NA/(12.0d0*EPS0*CL**2)*ECH**2/AMU/1.0d3
do k = 1,n
   dq = matmul(dmu, disp(:,k))
   inten(k) = fac*sum(dq**2)
enddo
end subroutine ir_intensities

subroutine print_ir(n, nskip, freq_cm, inten)
implicit none
integer,intent(in) :: n, nskip
real(8),intent(in) :: freq_cm(n), inten(n)
integer :: k
write(freq_unit,'(A)') ' [IR] IR spectrum (harmonic)'
write(freq_unit,'(A)') '   mode    freq(cm-1)   intensity(km/mol)'
do k = nskip+1,n
   write(freq_unit,'(I7,F14.4,F17.4)') k-nskip, freq_cm(k), inten(k)
enddo
write(freq_unit,'(A)') ' [IREND]'
end subroutine print_ir

subroutine write_molden_freq(fname, natoms, znum, coor_bohr, freq_cm, disp, redmass, nskip, inten)
implicit none
character(len=*),intent(in) :: fname
integer,intent(in) :: natoms, znum(natoms), nskip
real(8),intent(in) :: coor_bohr(3,natoms), freq_cm(3*natoms), disp(3*natoms,3*natoms), redmass(3*natoms)
real(8),intent(in),optional :: inten(3*natoms)
integer :: u, k, ia, c0, n
n = 3*natoms
open(newunit=u, file=trim(fname), status='replace', action='write')
write(u,'(A)') '[Molden Format]'
write(u,'(A)') '[FREQ]'
do k = nskip+1,n
   write(u,'(F12.4)') freq_cm(k)
enddo
write(u,'(A)') '[FR-COORD]'
do ia = 1,natoms
   write(u,'(A2,3F16.8)') elem_sym(znum(ia)), coor_bohr(:,ia)
enddo
write(u,'(A)') '[FR-NORM-COORD]'
do k = nskip+1,n
   write(u,'(A,I0)') 'vibration ', k-nskip
   do ia = 1,natoms
      c0 = 3*(ia-1)
      write(u,'(3F14.8)') disp(c0+1:c0+3,k)*sqrt(redmass(k))
   enddo
enddo
if (present(inten)) then
   write(u,'(A)') '[INT]'
   do k = nskip+1,n
      write(u,'(F14.6)') inten(k)
   enddo
endif
close(u)
end subroutine write_molden_freq

subroutine write_modes_csv(fname, natoms, znum, freq_cm, disp, redmass, nskip, inten)
implicit none
character(len=*),intent(in) :: fname
integer,intent(in) :: natoms, znum(natoms), nskip
real(8),intent(in) :: freq_cm(3*natoms), disp(3*natoms,3*natoms), redmass(3*natoms), inten(3*natoms)
integer :: u, k, ia, c0
real(8) :: d(3)
open(newunit=u, file=trim(fname), status='replace', action='write')
write(u,'(A)') 'mode,freq_cm1,ir_kmmol,redmass_amu,atom,element,dx,dy,dz,amplitude'
do k = nskip+1,3*natoms
   do ia = 1,natoms
      c0 = 3*(ia-1)
      d = disp(c0+1:c0+3,k)*sqrt(redmass(k))
      write(u,'(I0,",",F0.3,",",F0.5,",",F0.5,",",I0,",",A,",",3(ES14.6,","),ES14.6)') &
         k-nskip, freq_cm(k), inten(k), redmass(k), ia, trim(elem_sym(znum(ia))), d, sqrt(sum(d**2))
   enddo
enddo
close(u)
end subroutine write_modes_csv

function elem_sym(z) result(s)
implicit none
integer,intent(in) :: z
character(len=2) :: s
character(len=2),parameter :: t(86) = [character(len=2) :: &
 'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca', &
 'Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr', &
 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg', &
 'Tl','Pb','Bi','Po','At','Rn']
if (z >= 1 .and. z <= 86) then
   s = t(z)
else
   s = 'X '
endif
end function elem_sym

function rot_symmetry_number(natoms, znum, r, linear) result(sigma)
implicit none
integer,intent(in) :: natoms, znum(natoms)
real(8),intent(in) :: r(3,natoms)
logical,intent(in) :: linear
integer :: sigma
real(8),parameter :: tol = 0.2d0, pi = 3.14159265358979d0
real(8),allocatable :: ax(:,:), rots(:,:,:)
integer :: nax, i, j, n, k, nrot, l, nax0
real(8) :: v(3), nv, Rm(3,3), c, sn, th, u(3), cen(3)
real(8),allocatable :: rc(:,:)
logical :: ok, dup

if (linear) then
   sigma = 2
   do i = 1,natoms
      ok = .false.
      do j = 1,natoms
         if (znum(j) == znum(i) .and. sum((r(:,j)+r(:,i))**2) < tol**2) ok = .true.
      enddo
      if (.not. ok) sigma = 1
   enddo
   return
endif
allocate(ax(3, natoms + natoms*natoms + 3 + 20000), rots(3,3,1024), rc(3,natoms))
cen = sum(r, dim=2) / dble(natoms)
do i = 1,natoms
   rc(:,i) = r(:,i) - cen
enddo
nax = 0
do i = 1,natoms
   call addax(rc(:,i))
   do j = i+1,natoms
      call addax(rc(:,i)+rc(:,j))
      call addax([rc(2,i)*rc(3,j)-rc(3,i)*rc(2,j), rc(3,i)*rc(1,j)-rc(1,i)*rc(3,j), rc(1,i)*rc(2,j)-rc(2,i)*rc(1,j)])
   enddo
enddo
if (natoms <= 80) then
   nax0 = nax
   do i = 1,natoms
      do j = i+1,natoms
         if (znum(j) /= znum(i) .or. abs(norm2(rc(:,j))-norm2(rc(:,i))) > 0.1d0) cycle
         do k = j+1,natoms
            if (znum(k) /= znum(i) .or. abs(norm2(rc(:,k))-norm2(rc(:,i))) > 0.1d0) cycle
            v = rc(:,i)+rc(:,j)+rc(:,k)
            nv = norm2(v)
            if (nv < 1d-3 .or. nax >= nax0 + 20000) cycle
            v = v/nv
            dup = .false.
            do l = nax0+1,nax
               if (abs(abs(sum(ax(:,l)*v)) - 1d0) < 1d-3) then; dup = .true.; exit; endif
            enddo
            if (.not. dup) then; nax = nax + 1; ax(:,nax) = v; endif
         enddo
      enddo
   enddo
endif
nrot = 1
rots(:,:,1) = 0d0
do i = 1,3; rots(i,i,1) = 1d0; enddo
do l = 1,nax
   u = ax(:,l)
   do n = 2,8
      do k = 1,n-1
         th = 2d0*pi*k/n
         c = cos(th); sn = sin(th)
         Rm(1,1) = c+u(1)*u(1)*(1-c);      Rm(1,2) = u(1)*u(2)*(1-c)-u(3)*sn; Rm(1,3) = u(1)*u(3)*(1-c)+u(2)*sn
         Rm(2,1) = u(2)*u(1)*(1-c)+u(3)*sn; Rm(2,2) = c+u(2)*u(2)*(1-c);      Rm(2,3) = u(2)*u(3)*(1-c)-u(1)*sn
         Rm(3,1) = u(3)*u(1)*(1-c)-u(2)*sn; Rm(3,2) = u(3)*u(2)*(1-c)+u(1)*sn; Rm(3,3) = c+u(3)*u(3)*(1-c)
         ok = .true.
         do i = 1,natoms
            v = matmul(Rm, rc(:,i))
            c = 1d9
            do j = 1,natoms
               if (znum(j) == znum(i)) c = min(c, sum((v-rc(:,j))**2))
            enddo
            if (c > tol**2) then; ok = .false.; exit; endif
         enddo
         if (.not. ok) cycle
         dup = .false.
         do j = 1,nrot
            if (maxval(abs(rots(:,:,j)-Rm)) < 0.1d0) then; dup = .true.; exit; endif
         enddo
         if (.not. dup .and. nrot < 1024) then
            nrot = nrot + 1
            rots(:,:,nrot) = Rm
         endif
      enddo
   enddo
enddo
sigma = nrot
contains
   subroutine addax(w)
   real(8),intent(in) :: w(3)
   real(8) :: nw
   nw = sqrt(sum(w**2))
   if (nw > 1d-3) then
      nax = nax + 1
      ax(:,nax) = w/nw
   endif
   end subroutine addax
end function rot_symmetry_number

subroutine print_thermo(natoms, znum, coor_bohr, freq_cm, nzero, mult)
implicit none
integer,intent(in) :: natoms, znum(natoms), nzero, mult
real(8),intent(in) :: coor_bohr(3,natoms), freq_cm(3*natoms)
real(8),parameter :: kB = 1.380649d-23, hP = 6.62607015d-34, NA = 6.02214076d23
real(8),parameter :: amu = 1.66053906660d-27, bohr = 0.529177210903d-10
real(8),parameter :: Rg = 8.314462618d0, Eh = 4.3597447222071d-18
real(8),parameter :: kcal = 4184d0, cm2K = 1.438776877d0, pi = 3.14159265358979d0
real(8) :: T, P, mtot, com(3), r(3,natoms), Iner(3,3), mom(3), work(64), mkg
real(8) :: qt, St, Et, Sr, Er, Sv, Ev, Cv, zpe, Se, th, x, qr, thr(3), Cvv, Crot
real(8) :: srrho, nuhz, mu, mup, sfree, wq
real(8) :: Etot, Htot, Gtot, Stot, e0
integer :: i, j, info, sigma, nv
logical :: linear, atom

T = freq_temp; P = freq_pres
atom = (natoms == 1)
mtot = 0d0; com = 0d0
do i = 1,natoms
   mtot = mtot + atomic_mass(znum(i))
   com = com + atomic_mass(znum(i))*coor_bohr(:,i)
enddo
com = com/mtot
do i = 1,natoms; r(:,i) = coor_bohr(:,i) - com; enddo
Iner = 0d0
do i = 1,natoms
   do j = 1,3
      Iner(j,j) = Iner(j,j) + atomic_mass(znum(i))*sum(r(:,i)**2)
   enddo
   Iner = Iner - atomic_mass(znum(i))*spread(r(:,i),2,3)*spread(r(:,i),1,3)
enddo
call dsyev('N','U',3,Iner,3,mom,work,64,info)
linear = (.not. atom) .and. mom(1) < 1d-4*max(mom(3),1d-12)
sigma = 1
if (.not. atom) sigma = rot_symmetry_number(natoms, znum, r, linear)

mkg = mtot*amu
qt = (2d0*pi*mkg*kB*T/hP**2)**1.5d0 * kB*T/P
St = Rg*(log(qt) + 2.5d0)
Et = 1.5d0*Rg*T
Cv = 1.5d0*Rg
Sr = 0d0; Er = 0d0; Crot = 0d0
if (.not. atom) then
   thr = 0d0
   do i = 1,3
      if (mom(i) > 1d-8) thr(i) = (hP/(2d0*pi))**2/(2d0*mom(i)*amu*bohr**2*kB)
   enddo
   if (linear) then
      qr = T/(sigma*thr(3))
      Sr = Rg*(log(qr) + 1d0); Er = Rg*T; Crot = Rg
   else
      qr = sqrt(pi)/sigma * T**1.5d0/sqrt(thr(1)*thr(2)*thr(3))
      Sr = Rg*(log(qr) + 1.5d0); Er = 1.5d0*Rg*T; Crot = 1.5d0*Rg
   endif
endif
Sv = 0d0; Ev = 0d0; zpe = 0d0; nv = 0; Cvv = 0d0
do i = nzero+1, 3*natoms
   if (freq_cm(i) <= 1d0) cycle
   nv = nv + 1
   th = cm2K*freq_cm(i); x = th/T
   zpe = zpe + 0.5d0*Rg*th
   Ev = Ev + Rg*th*(0.5d0 + 1d0/(exp(x)-1d0))
   srrho = Rg*(x/(exp(x)-1d0) - log(1d0-exp(-x)))
   nuhz = freq_cm(i)*2.99792458d10
   mu = hP/(8d0*pi**2*nuhz)
   mup = mu*1d-44/(mu + 1d-44)
   sfree = Rg*(0.5d0 + log(sqrt(8d0*pi**3*mup*kB*T/hP**2)))
   wq = 1d0/(1d0 + (100d0/freq_cm(i))**4)
   Sv = Sv + wq*srrho + (1d0-wq)*sfree
   Cvv = Cvv + Rg*x*x*exp(x)/(exp(x)-1d0)**2
enddo
Cv = Cv + Crot + Cvv
Se = Rg*log(dble(mult))
Etot = Et + Er + Ev
Stot = St + Sr + Sv + Se
Htot = Etot + Rg*T
Gtot = Htot - T*Stot
e0 = freq_energy
write(freq_unit,'(A)') ' [THERMO]'
write(freq_unit,'(A,F9.2,A,F10.1,A)') ' Ideal gas, rigid rotor, harmonic oscillator;  T =', T, ' K,  P =', P, ' Pa'
write(freq_unit,'(A,I3,A,I3)') ' Rotational symmetry number =', sigma, '   multiplicity =', mult
write(freq_unit,'(A,F14.4,A)') ' Total mass            =', mtot, ' amu'
write(freq_unit,'(A,F14.6,A,F10.3,A)') ' Zero-point correction =', zpe/NA/Eh, ' Eh  (', zpe/kcal, ' kcal/mol)'
write(freq_unit,'(A,F14.6,A,F10.3,A)') ' Thermal correction to E =', (Etot)/NA/Eh, ' Eh  (', Etot/kcal, ' kcal/mol)'
write(freq_unit,'(A,F14.6,A,F10.3,A)') ' Thermal correction to H =', Htot/NA/Eh, ' Eh  (', Htot/kcal, ' kcal/mol)'
write(freq_unit,'(A,F14.6,A,F10.3,A)') ' Thermal correction to G =', Gtot/NA/Eh, ' Eh  (', Gtot/kcal, ' kcal/mol)'
write(freq_unit,'(A)') '                 E (kcal/mol)   Cv (cal/mol-K)   S (cal/mol-K)'
write(freq_unit,'(A,3F16.3)') ' Total         ', Etot/kcal, Cv/4.184d0, Stot/4.184d0
write(freq_unit,'(A,3F16.3)') ' Electronic    ', 0d0, 0d0, Se/4.184d0
write(freq_unit,'(A,3F16.3)') ' Translational ', Et/kcal, 1.5d0*Rg/4.184d0, St/4.184d0
write(freq_unit,'(A,3F16.3)') ' Rotational    ', Er/kcal, Crot/4.184d0, Sr/4.184d0
write(freq_unit,'(A,3F16.3)') ' Vibrational   ', Ev/kcal, Cvv/4.184d0, Sv/4.184d0
if (e0 /= 0d0) then
   write(freq_unit,'(A,F18.8)') ' Sum of electronic and zero-point energies =', e0 + zpe/NA/Eh
   write(freq_unit,'(A,F18.8)') ' Sum of electronic and thermal energies    =', e0 + Etot/NA/Eh
   write(freq_unit,'(A,F18.8)') ' Sum of electronic and thermal enthalpies  =', e0 + Htot/NA/Eh
   write(freq_unit,'(A,F18.8)') ' Sum of electronic and thermal free energies =', e0 + Gtot/NA/Eh
endif
write(freq_unit,'(A)') ' [THERMOEND]'
end subroutine print_thermo

end module mod_hessian_freq
