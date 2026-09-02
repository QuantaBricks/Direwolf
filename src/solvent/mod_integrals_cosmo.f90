! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_integrals_cosmo: the ONE integral routine COSMO needs from the

submodule (mod_integrals) cosmo_integrals_impl
implicit none
contains

module subroutine cosmo_build_one_tess_matrix(tess_pos_bohr, Bk)
use MOL_info, only: Natoms
implicit none
real(8),intent(in) :: tess_pos_bohr(3)
real(8),intent(out) :: Bk(:,:)

integer,allocatable :: atm_k(:,:)
real(8),allocatable :: env_k(:)
real(8),allocatable :: buf1eV(:,:)
integer :: nAtoms_k, offpoint_k
integer :: i,j,di,dj,n,shls(4)

nAtoms_k = Natoms + 1
allocate(atm_k(6,nAtoms_k))
atm_k = 0
atm_k(:,1:Natoms) = atm(:,1:Natoms)
atm_k(1,1:Natoms) = 0

allocate(env_k(size(env)+4))
env_k = 0.0d0
env_k(1:size(env)) = env(1:size(env))
offpoint_k = size(env)
atm_k(2,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = tess_pos_bohr(1)
env_k(offpoint_k+2) = tess_pos_bohr(2)
env_k(offpoint_k+3) = tess_pos_bohr(3)
offpoint_k = offpoint_k + 3
atm_k(3,nAtoms_k) = 3
atm_k(5,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = 1.0d0
offpoint_k = offpoint_k + 1

n = 0
do i = 0,nBases-1
   n = max(n, cgto_engine(i,bas))
enddo
allocate(buf1eV(n,n))
Bk = 0.0d0
do i = 0,nBases-1
   do j = i,nBases-1
      shls(1) = i
      shls(2) = j
      di = cgto_engine(i,bas)
      dj = cgto_engine(j,bas)
      call nuc1e_engine(buf1eV(1:di,1:dj), shls, atm_k, nAtoms_k, bas, nBases, env_k, 0_8)
      call store1e(shls, di, dj, bas, nBases, buf1eV(1:di,1:dj), NorVEC, size(Bk,1), Bk)
   enddo
enddo

deallocate(buf1eV, atm_k, env_k)
end subroutine cosmo_build_one_tess_matrix

module subroutine esp_at_grid_batch(npts, grid_bohr, Ptot, V_elec)
use MOL_info, only: Natoms
implicit none
integer,intent(in) :: npts
real(8),intent(in) :: grid_bohr(3,npts)
real(8),intent(in) :: Ptot(:,:)
real(8),intent(out) :: V_elec(npts)

integer,parameter :: CHUNK = 256
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),allocatable :: env_g(:), buf(:), int_cache(:)
integer,allocatable :: ao_offset(:)
integer :: nenv, off_grids, ip, ish, jsh, di, dj, is, js, nmax
integer :: g0, g1, ng, g, p, q, e1, e2, shls(4)
real(8) :: w, nij

nenv = size(env)
off_grids = nenv
allocate(env_g(nenv + 3*npts))
env_g = 0.0d0
env_g(1:nenv) = env(1:nenv)
do ip = 1,npts
   env_g(off_grids + 3*(ip-1) + 1) = grid_bohr(1,ip)
   env_g(off_grids + 3*(ip-1) + 2) = grid_bohr(2,ip)
   env_g(off_grids + 3*(ip-1) + 3) = grid_bohr(3,ip)
enddo
env_g(12) = dble(npts)
env_g(13) = dble(off_grids)

allocate(ao_offset(0:nBases-1))
ip = 0
nmax = 0
do ish = 0,nBases-1
   ao_offset(ish) = ip
   ip = ip + cgto_engine(ish,bas)
   nmax = max(nmax, cgto_engine(ish,bas))
enddo

!$omp parallel private(buf,int_cache,g0,g1,ng,ish,jsh,di,dj,is,js,shls,g,p,q,e1,e2,w,nij)
allocate(buf(nmax*nmax*CHUNK), int_cache(INT_CACHE_SIZE))
!$omp do schedule(dynamic)
do g0 = 1,npts,CHUNK
   g1 = min(g0+CHUNK-1, npts)
   ng = g1 - g0 + 1
   V_elec(g0:g1) = 0.0d0
   do ish = 0,nBases-1
      di = cgto_engine(ish,bas)
      is = ao_offset(ish)
      do jsh = ish,nBases-1
         dj = cgto_engine(jsh,bas)
         js = ao_offset(jsh)
         shls(1) = ish
         shls(2) = jsh
         shls(3) = g0-1
         shls(4) = g1
         call grids1e_engine_cached(buf, shls, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         do q = 1,dj
            e2 = js + q
            do p = 1,di
               e1 = is + p
               if (ish .eq. jsh .and. e2 .lt. e1) cycle
               nij = NorVEC(e1)*NorVEC(e2)
               if (nij .eq. 0.0d0) cycle
               w = Ptot(e1,e2)*nij
               if (e2 .gt. e1) w = 2.0d0*w
               do g = 1,ng
                  V_elec(g0+g-1) = V_elec(g0+g-1) - w*buf(g + ng*((p-1) + di*(q-1)))
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf, int_cache)
!$omp end parallel

deallocate(env_g, ao_offset)
end subroutine esp_at_grid_batch

module subroutine esp_at_grid_batch_df(npts, grid_bohr, Ptot, V_elec, ok)
use MOL_info, only: Natoms
implicit none
integer,intent(in) :: npts
real(8),intent(in) :: grid_bohr(3,npts)
real(8),intent(in) :: Ptot(:,:)
real(8),intent(out) :: V_elec(npts)
logical,intent(out) :: ok

integer,parameter :: CHUNK = 256
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),parameter :: UNIT_S_COEFF = 2.0d0*sqrt(acos(-1.0d0))
real(8),allocatable :: envB(:), buf(:), int_cache(:)
integer,allocatable :: basB(:,:)
integer :: nenvB, off_grids, off_dummy, ip, ksh, dk, offP, nmax
integer :: g0, g1, ng, g, r, shls(4), nshB, idummy

ok = .false.
if (.not. df_built) return
if (.not. df_direct_cvec_cache_valid) return
if (.not. allocated(df_direct_cvec_cache)) return
if (.not. allocated(df_direct_cvec_cache_Ptot)) return
if (size(df_direct_cvec_cache_Ptot,1) .ne. size(Ptot,1)) return
if (any(df_direct_cvec_cache_Ptot .ne. Ptot)) return
ok = .true.

nshB = nBases + nBasesAux + 1
idummy = nBases + nBasesAux
nenvB = size(envDF)
off_dummy = nenvB
off_grids = nenvB + 2
allocate(basB(8,nshB), envB(nenvB + 2 + 3*npts))
basB = 0
basB(1:8,1:nBases+nBasesAux) = basDF(1:8,1:nBases+nBasesAux)
envB = 0.0d0
envB(1:nenvB) = envDF(1:nenvB)

basB(1,nshB) = 0
basB(2,nshB) = 0
basB(3,nshB) = 1
basB(4,nshB) = 1
basB(6,nshB) = off_dummy
basB(7,nshB) = off_dummy + 1
envB(off_dummy+1) = 0.0d0
envB(off_dummy+2) = UNIT_S_COEFF

do ip = 1,npts
   envB(off_grids + 3*(ip-1) + 1) = grid_bohr(1,ip)
   envB(off_grids + 3*(ip-1) + 2) = grid_bohr(2,ip)
   envB(off_grids + 3*(ip-1) + 3) = grid_bohr(3,ip)
enddo
envB(12) = dble(npts)
envB(13) = dble(off_grids)

nmax = 0
do ksh = nBases,nBases+nBasesAux-1
   nmax = max(nmax, df_shell_ncgto(ksh))
enddo

!$omp parallel private(buf,int_cache,g0,g1,ng,ksh,dk,offP,shls,g,r)
allocate(buf(nmax*CHUNK), int_cache(INT_CACHE_SIZE))
!$omp do schedule(dynamic)
do g0 = 1,npts,CHUNK
   g1 = min(g0+CHUNK-1, npts)
   ng = g1 - g0 + 1
   V_elec(g0:g1) = 0.0d0
   do ksh = nBases,nBases+nBasesAux-1
      dk = df_shell_ncgto(ksh)
      offP = df_shell_offset(ksh)
      shls(1) = ksh
      shls(2) = idummy
      shls(3) = g0-1
      shls(4) = g1
      call grids1e_engine_cached(buf, shls, atm, size(atm,2), basB, nshB, envB, 0_8, int_cache)
      do r = 1,dk
         do g = 1,ng
            V_elec(g0+g-1) = V_elec(g0+g-1) &
               - df_direct_cvec_cache(offP+r)*NorVECAux(offP+r)*buf(g + ng*(r-1))
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf, int_cache)
!$omp end parallel

deallocate(basB, envB)
end subroutine esp_at_grid_batch_df

module subroutine cosmo_build_one_tess_matrix_shellderiv(tess_pos_bohr, dBk)
use MOL_info, only: Natoms
implicit none
real(8),intent(in) :: tess_pos_bohr(3)
real(8),intent(out) :: dBk(:,:,:)

integer,allocatable :: atm_k(:,:)
real(8),allocatable :: env_k(:)
real(8),allocatable :: buf1eV(:,:,:)
integer,allocatable :: ao_offset(:), shell_dim(:)
integer :: nAtoms_k, offpoint_k
integer :: i,j,di,dj,shls(2)
integer :: si_off, num_off

nAtoms_k = Natoms + 1
allocate(atm_k(6,nAtoms_k))
atm_k = 0
atm_k(:,1:Natoms) = atm(:,1:Natoms)
atm_k(1,1:Natoms) = 0

allocate(env_k(size(env)+4))
env_k = 0.0d0
env_k(1:size(env)) = env(1:size(env))
offpoint_k = size(env)
atm_k(2,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = tess_pos_bohr(1)
env_k(offpoint_k+2) = tess_pos_bohr(2)
env_k(offpoint_k+3) = tess_pos_bohr(3)
offpoint_k = offpoint_k + 3
atm_k(3,nAtoms_k) = 3
atm_k(5,nAtoms_k) = offpoint_k
env_k(offpoint_k+1) = 1.0d0
offpoint_k = offpoint_k + 1

allocate(ao_offset(0:nBases-1), shell_dim(0:nBases-1))
num_off = 0
do si_off = 0,nBases-1
   ao_offset(si_off) = num_off
   shell_dim(si_off) = cgto_engine(si_off, bas)
   num_off = num_off + shell_dim(si_off)
enddo

dBk = 0.0d0
do i = 0,nBases-1
   do j = 0,nBases-1
      shls(1) = i
      shls(2) = j
      di = shell_dim(i)
      dj = shell_dim(j)
      allocate(buf1eV(di,dj,3))
      call nuc1e_ip_engine(buf1eV, shls, atm_k, nAtoms_k, bas, nBases, env_k, 0_8)
      call store1edrv(shls, di, dj, ao_offset, nBases, buf1eV, NorVEC, size(dBk,1), dBk)
      deallocate(buf1eV)
   enddo
enddo

deallocate(atm_k, env_k, ao_offset, shell_dim)
end subroutine cosmo_build_one_tess_matrix_shellderiv

end submodule cosmo_integrals_impl
