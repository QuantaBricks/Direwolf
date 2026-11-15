! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Gaussian-format cube files of the converged electron density and electrostatic potential.

module mod_cube
implicit none
private
public :: write_cube_files
public :: cube_grid, cube_screen_begin, cube_screen_end, cube_ao_values, cube_chunk_points, cube_write_values

interface
   subroutine GTOeval_point(coor_grid, val0, val1)
   use MOL_info, only: nconts
   real(8),intent(in) :: coor_grid(3)
   real(8),intent(out) :: val0(nconts), val1(nconts,3)
   end subroutine GTOeval_point
end interface

contains

subroutine write_cube_files(stem, do_den, do_esp, esp_exact)
use MOL_info
use GRID_info, only: rcut2_shared, maxShell_shared, maxGauss_shared
use mod_integrals, only: esp_at_grid_batch, esp_at_grid_batch_df, df_build_coulomb_core, engine_use_df_j
implicit none
INCLUDE 'parameter.h'
character(len=*),intent(in) :: stem
logical,intent(in) :: do_den, do_esp, esp_exact
real(8),parameter :: MARGIN = 6.0d0, HMIN = 0.2d0
integer,parameter :: NMAX = 3000000, NB = 4096
real(8) :: lo(3), hi(3), org(3), h, vol, t1
integer(8) :: c0, c1, crate
integer :: n(3), ia, i, j, ntot, nocc, ip0, nb_now, k, ix, iy, iz, ip, info, lwork, liwork
real(8),allocatable :: Ptot(:,:), L(:,:), ev(:), work(:), Phi(:,:), T(:,:), v0(:,:), v1(:,:,:), xyz(:,:), vals(:), Vc(:)
integer,allocatable :: iwork(:)
logical :: own_rcut, use_fit, ok_fit
real(8),allocatable :: Jtmp(:,:), xs(:,:), Vf(:), Ve(:)
integer,parameter :: NSAMP = 512
real(8) :: dist

call system_clock(c0, crate)
do ia = 1,Natoms
   if (ia .eq. 1) then
      lo = atoms(1)%coor*ans2bohr; hi = lo
   else
      lo = min(lo, atoms(ia)%coor*ans2bohr); hi = max(hi, atoms(ia)%coor*ans2bohr)
   endif
enddo
lo = lo - MARGIN; hi = hi + MARGIN
vol = product(hi - lo)
h = max(HMIN, (vol/real(NMAX,8))**(1.0d0/3.0d0))
n = int((hi - lo)/h) + 1
org = lo
ntot = n(1)*n(2)*n(3)

own_rcut = .not. allocated(rcut2_shared)
if (own_rcut) then
   maxShell_shared = 0; maxGauss_shared = 0
   do ia = 1,Natoms
      maxShell_shared = max(maxShell_shared, size(atoms(ia)%shell))
      do j = 1,size(atoms(ia)%shell)
         maxGauss_shared = max(maxGauss_shared, atoms(ia)%shell(j)%nGauss)
      enddo
   enddo
   allocate(rcut2_shared(maxShell_shared,Natoms)); rcut2_shared = 0
   do ia = 1,Natoms
      do j = 1,size(atoms(ia)%shell)
         rcut2_shared(j,ia) = (sqrt(-log(1.0d-9)/minval(atoms(ia)%shell(j)%exponents)) + 2.0d0)**2
      enddo
   enddo
endif

allocate(Ptot(nConts,nConts))
Ptot = Pa + Pb

nocc = 0
if (do_den) then
   allocate(L(nConts,nConts), ev(nConts))
   L = Ptot
   lwork = 1 + 6*nConts + 2*nConts**2; liwork = 3 + 5*nConts
   allocate(work(lwork), iwork(liwork))
   call dsyevd('V', 'U', nConts, L, nConts, ev, work, lwork, iwork, liwork, info)
   deallocate(work, iwork)
   do k = 1,nConts
      if (ev(k) > 1.0d-10) then
         nocc = nocc + 1
         L(:,nocc) = L(:,k)*sqrt(ev(k))
      endif
   enddo
endif

allocate(vals(ntot))
if (do_den) then
   allocate(Phi(NB,nConts), T(NB,max(nocc,1)), xyz(3,NB))
   do ip0 = 1,ntot,NB
      nb_now = min(NB, ntot - ip0 + 1)
      call chunk_points(ip0, nb_now, n, org, h, xyz)
      !$omp parallel private(v0, v1, k)
      allocate(v0(nConts,1), v1(nConts,3,1))
      !$omp do schedule(dynamic,16)
      do k = 1,nb_now
         call GTOeval_point(xyz(:,k), v0(:,1), v1(:,:,1))
         Phi(k,:) = v0(:,1)
      enddo
      !$omp end do
      deallocate(v0, v1)
      !$omp end parallel
      call dgemm('N', 'N', nb_now, nocc, nConts, 1.0d0, Phi, NB, L, nConts, 0.0d0, T, NB)
      do k = 1,nb_now
         vals(ip0+k-1) = sum(T(k,1:nocc)**2)
      enddo
   enddo
   deallocate(Phi, T, xyz)
   call system_clock(c1); print '(A,F8.2,A)', '  cube density:  ', real(c1 - c0,8)/real(crate,8), ' s'
   call write_one(trim(stem)//'.den.cube', 'electron density (e/bohr^3)')
endif

if (do_esp) then
   allocate(xyz(3,NB), Vc(NB))
   use_fit = engine_use_df_j .and. .not. esp_exact
   if (use_fit) then
      allocate(Jtmp(nConts,nConts))
      call df_build_coulomb_core(nConts, Ptot, Jtmp, .false.)
      deallocate(Jtmp)
      allocate(xs(3,NSAMP), Vf(NSAMP), Ve(NSAMP))
      do k = 1,NSAMP
         call chunk_points(1 + int((int(k-1,8)*int(ntot-1,8))/int(max(NSAMP-1,1),8)), 1, n, org, h, xs(:,k:k))
      enddo
      call esp_at_grid_batch_df(NSAMP, xs, Ptot, Vf, ok_fit)
      if (ok_fit) then
         call esp_at_grid_batch(NSAMP, xs, Ptot, Ve)
         print '(A)', '  cube potential from the RI-J fitted density; against the exact density at 512 grid points:'
         print '(A,ES9.2,A,ES9.2,A)', '    max |dV| ', maxval(abs(Vf - Ve)), ' rms ', sqrt(sum((Vf - Ve)**2)/NSAMP), ' Eh/e'
      else
         use_fit = .false.
      endif
      deallocate(xs, Vf, Ve)
   endif
   do ip0 = 1,ntot,NB
      nb_now = min(NB, ntot - ip0 + 1)
      call chunk_points(ip0, nb_now, n, org, h, xyz)
      if (use_fit) then
         call esp_at_grid_batch_df(nb_now, xyz(:,1:nb_now), Ptot, Vc(1:nb_now), ok_fit)
      else
         call esp_at_grid_batch(nb_now, xyz(:,1:nb_now), Ptot, Vc(1:nb_now))
      endif
      !$omp parallel do private(ia, dist) schedule(static)
      do k = 1,nb_now
         do ia = 1,Natoms
            dist = sqrt(sum((xyz(:,k) - atoms(ia)%coor*ans2bohr)**2))
            Vc(k) = Vc(k) + real(atoms(ia)%charge - atoms(ia)%ecpCoreElec,8)/max(dist, 0.05d0)
         enddo
      enddo
      !$omp end parallel do
      vals(ip0:ip0+nb_now-1) = Vc(1:nb_now)
   enddo
   deallocate(xyz, Vc)
   call system_clock(c1); print '(A,F8.2,A)', '  cube potential:', real(c1 - c0,8)/real(crate,8), ' s (cumulative)'
   call write_one(trim(stem)//'.esp.cube', 'electrostatic potential (Eh/e)')
endif

if (own_rcut) deallocate(rcut2_shared)
call system_clock(c1)
t1 = real(c1 - c0,8)/real(crate,8)
print '(A,F8.2,A)', '  cube files: ', t1, ' s'

contains

subroutine write_one(fname, what)
character(len=*),intent(in) :: fname, what
integer :: iu, ia2, ix2, iy2, iz2, c
open(newunit=iu, file=fname, status='replace', action='write')
write(iu,'(A)') 'Direwolf cube: '//what
write(iu,'(A,I0,A,3(I0,1X),A,F6.3,A)') 'grid ', ntot, ' points: ', n, 'spacing ', h, ' bohr'
write(iu,'(I5,3F12.6)') Natoms, org
write(iu,'(I5,3F12.6)') n(1), h, 0.0d0, 0.0d0
write(iu,'(I5,3F12.6)') n(2), 0.0d0, h, 0.0d0
write(iu,'(I5,3F12.6)') n(3), 0.0d0, 0.0d0, h
do ia2 = 1,Natoms
   write(iu,'(I5,4F12.6)') atoms(ia2)%charge, real(atoms(ia2)%charge - atoms(ia2)%ecpCoreElec,8), &
        atoms(ia2)%coor*ans2bohr
enddo
c = 0
do ix2 = 1,n(1)
   do iy2 = 1,n(2)
      do iz2 = 1,n(3)
         c = c + 1
         write(iu,'(ES14.5)', advance='no') vals(((ix2-1)*n(2) + (iy2-1))*n(3) + iz2)
         if (mod(iz2,6) .eq. 0 .or. iz2 .eq. n(3)) write(iu,*)
      enddo
   enddo
enddo
close(iu)
print '(3A,I0,A,F5.3,A)', '  cube file written: ', trim(fname), ', ', ntot, ' points, spacing ', h, ' bohr'
end subroutine write_one

end subroutine write_cube_files

subroutine cube_grid(margin, hmin, nmax, n, org, h)
use MOL_info, only: Natoms, atoms
implicit none
INCLUDE 'parameter.h'
real(8),intent(in) :: margin, hmin
integer,intent(in) :: nmax
integer,intent(out) :: n(3)
real(8),intent(out) :: org(3), h
real(8) :: lo(3), hi(3)
integer :: ia
do ia = 1,Natoms
   if (ia .eq. 1) then
      lo = atoms(1)%coor*ans2bohr; hi = lo
   else
      lo = min(lo, atoms(ia)%coor*ans2bohr); hi = max(hi, atoms(ia)%coor*ans2bohr)
   endif
enddo
lo = lo - margin; hi = hi + margin
h = max(hmin, (product(hi - lo)/real(nmax,8))**(1.0d0/3.0d0))
n = int((hi - lo)/h) + 1
org = lo
end subroutine cube_grid

subroutine cube_screen_begin(own)
use MOL_info, only: Natoms, atoms
use GRID_info, only: rcut2_shared, maxShell_shared, maxGauss_shared
implicit none
logical,intent(out) :: own
integer :: ia, j
own = .not. allocated(rcut2_shared)
if (own) then
   maxShell_shared = 0; maxGauss_shared = 0
   do ia = 1,Natoms
      maxShell_shared = max(maxShell_shared, size(atoms(ia)%shell))
      do j = 1,size(atoms(ia)%shell)
         maxGauss_shared = max(maxGauss_shared, atoms(ia)%shell(j)%nGauss)
      enddo
   enddo
   allocate(rcut2_shared(maxShell_shared,Natoms)); rcut2_shared = 0
   do ia = 1,Natoms
      do j = 1,size(atoms(ia)%shell)
         rcut2_shared(j,ia) = (sqrt(-log(1.0d-9)/minval(atoms(ia)%shell(j)%exponents)) + 2.0d0)**2
      enddo
   enddo
endif
end subroutine cube_screen_begin

subroutine cube_screen_end(own)
use GRID_info, only: rcut2_shared
implicit none
logical,intent(in) :: own
if (own) deallocate(rcut2_shared)
end subroutine cube_screen_end

subroutine cube_ao_values(nb, xyz, Phi)
use MOL_info, only: nConts
implicit none
integer,intent(in) :: nb
real(8),intent(in) :: xyz(3,nb)
real(8),intent(out) :: Phi(:,:)
real(8),allocatable :: v0(:,:), v1(:,:,:)
integer :: k
!$omp parallel private(v0, v1, k)
allocate(v0(nConts,1), v1(nConts,3,1))
!$omp do schedule(dynamic,16)
do k = 1,nb
   call GTOeval_point(xyz(:,k), v0(:,1), v1(:,:,1))
   Phi(k,1:nConts) = v0(:,1)
enddo
!$omp end do
deallocate(v0, v1)
!$omp end parallel
end subroutine cube_ao_values

subroutine cube_chunk_points(ip0, nb, n, org, h, xyz)
integer,intent(in) :: ip0, nb, n(3)
real(8),intent(in) :: org(3), h
real(8),intent(out) :: xyz(3,*)
call chunk_points(ip0, nb, n, org, h, xyz)
end subroutine cube_chunk_points

subroutine cube_write_values(fname, what, vals, n, org, h)
use MOL_info, only: Natoms, atoms
implicit none
INCLUDE 'parameter.h'
character(len=*),intent(in) :: fname, what
integer,intent(in) :: n(3)
real(8),intent(in) :: vals(:), org(3), h
integer :: iu, ia2, ix2, iy2, iz2
open(newunit=iu, file=fname, status='replace', action='write')
write(iu,'(A)') 'Direwolf cube: '//what
write(iu,'(A,I0,A,3(I0,1X),A,F6.3,A)') 'grid ', n(1)*n(2)*n(3), ' points: ', n, 'spacing ', h, ' bohr'
write(iu,'(I5,3F12.6)') Natoms, org
write(iu,'(I5,3F12.6)') n(1), h, 0.0d0, 0.0d0
write(iu,'(I5,3F12.6)') n(2), 0.0d0, h, 0.0d0
write(iu,'(I5,3F12.6)') n(3), 0.0d0, 0.0d0, h
do ia2 = 1,Natoms
   write(iu,'(I5,4F12.6)') atoms(ia2)%charge, real(atoms(ia2)%charge - atoms(ia2)%ecpCoreElec,8), &
        atoms(ia2)%coor*ans2bohr
enddo
do ix2 = 1,n(1)
   do iy2 = 1,n(2)
      do iz2 = 1,n(3)
         write(iu,'(ES14.5)', advance='no') vals(((ix2-1)*n(2) + (iy2-1))*n(3) + iz2)
         if (mod(iz2,6) .eq. 0 .or. iz2 .eq. n(3)) write(iu,*)
      enddo
   enddo
enddo
close(iu)
end subroutine cube_write_values

subroutine chunk_points(ip0, nb, n, org, h, xyz)
integer,intent(in) :: ip0, nb, n(3)
real(8),intent(in) :: org(3), h
real(8),intent(out) :: xyz(3,*)
integer :: k, ip, ix, iy, iz
do k = 1,nb
   ip = ip0 + k - 2
   iz = mod(ip, n(3)); ip = ip/n(3)
   iy = mod(ip, n(2)); ix = ip/n(2)
   xyz(:,k) = org + h*[real(ix,8), real(iy,8), real(iz,8)]
enddo
end subroutine chunk_points

end module mod_cube
