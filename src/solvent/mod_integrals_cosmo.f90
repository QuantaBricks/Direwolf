! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

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
