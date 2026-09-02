! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_init: cosmo_init - builds the cavity, the (static) conductor

module mod_cosmo_init
use mod_cosmo_cavity, only: cosmo_build_cavity, bohr2ang, PI_cosmo
use mod_cosmo_state
implicit none

contains

subroutine cosmo_init(natoms, Z, ecpCoreElec, coor_ang, epsilon, radii_scale, avg_area_ang2, &
                       cavity_type, rsolv_ang, ks_nseg, ks_nface, nConts_in, molecular_charge, smd_alpha)
use mod_cosmo_cavity, only: cosmo_build_cavity_ks, cosmo_build_cavity_yk, cosmo_build_cavity_iswig, &
                             cosmo_smd_radius_bohr
use mod_integrals, only: cosmo_build_one_tess_matrix
use mod_meminfo, only: get_memory_budget_bytes, get_store_threshold_bytes
implicit none
integer,intent(in) :: natoms
integer,intent(in) :: Z(natoms)
integer,intent(in) :: ecpCoreElec(natoms)
real(8),intent(in) :: coor_ang(3,natoms)
real(8),intent(in) :: epsilon, radii_scale
real(8),intent(in) :: avg_area_ang2
real(8),intent(in),optional :: smd_alpha
character(len=*),intent(in) :: cavity_type
real(8),intent(in) :: rsolv_ang
integer,intent(in) :: ks_nseg, ks_nface
integer,intent(in) :: nConts_in
real(8),intent(in) :: molecular_charge

real(8),allocatable :: coor_bohr(:,:), Rvdw_smd(:)
real(8) :: Zeff, d
integer :: ia, k, l, LWORK, info
real(8),allocatable :: WORK(:)
real(8) :: wsize(1)

allocate(coor_bohr(3,natoms))
coor_bohr = coor_ang * (1.0d0/bohr2ang)

saved_natoms = natoms
allocate(saved_Z(natoms), saved_ecpCoreElec(natoms))
saved_Z = Z
saved_ecpCoreElec = ecpCoreElec
saved_radii_scale = radii_scale
saved_avg_area_ang2 = avg_area_ang2
saved_Qm = molecular_charge
saved_cavity_type = trim(cavity_type)
use_smd_radii = present(smd_alpha)
if (use_smd_radii) saved_smd_alpha = smd_alpha

if (trim(cavity_type) == 'ks1993') then
   call cosmo_build_cavity_ks(natoms, Z, coor_bohr, rsolv_ang, ks_nseg, ks_nface, &
                              tess_coor, tess_area, tess_atom, n_tess)
else if (trim(cavity_type) == 'yk1999') then
   call cosmo_build_cavity_yk(natoms, Z, coor_bohr, ks_nseg, &
                              tess_coor, tess_area, tess_atom, n_tess)
else if (trim(cavity_type) == 'iswig') then
   call cosmo_build_cavity_iswig(natoms, Z, coor_bohr, rsolv_ang, &
                                 tess_coor, tess_area, tess_atom, n_tess)
else if (use_smd_radii) then
   allocate(Rvdw_smd(natoms))
   do ia = 1,natoms
      Rvdw_smd(ia) = cosmo_smd_radius_bohr(Z(ia), smd_alpha)
   enddo
   call cosmo_build_cavity(natoms, Z, coor_bohr, radii_scale, avg_area_ang2, &
                            tess_coor, tess_area, tess_atom, n_tess, Rvdw_smd)
   deallocate(Rvdw_smd)
else
   call cosmo_build_cavity(natoms, Z, coor_bohr, radii_scale, avg_area_ang2, &
                            tess_coor, tess_area, tess_atom, n_tess)
endif

block
integer :: kk, n_kept
real(8) :: area_ang2
real(8),allocatable :: tc2(:,:), ta2(:)
integer,allocatable :: tat2(:)
real(8),parameter :: MIN_TESS_AREA_ANG2 = 1.0d-4
n_kept = 0
do kk = 1,n_tess
   area_ang2 = tess_area(kk)*bohr2ang*bohr2ang
   if (area_ang2 >= MIN_TESS_AREA_ANG2) n_kept = n_kept+1
enddo
if (n_kept < n_tess) then
   allocate(tc2(3,n_kept), ta2(n_kept), tat2(n_kept))
   n_kept = 0
   do kk = 1,n_tess
      area_ang2 = tess_area(kk)*bohr2ang*bohr2ang
      if (area_ang2 >= MIN_TESS_AREA_ANG2) then
         n_kept = n_kept+1
         tc2(:,n_kept) = tess_coor(:,kk)
         ta2(n_kept) = tess_area(kk)
         tat2(n_kept) = tess_atom(kk)
      endif
   enddo
   print '("COSMO: dropped ",I0," near-degenerate tesserae (area<",ES9.2," Ang^2) before charge solve")', &
         n_tess-n_kept, MIN_TESS_AREA_ANG2
   call move_alloc(tc2, tess_coor)
   call move_alloc(ta2, tess_area)
   call move_alloc(tat2, tess_atom)
   n_tess = n_kept
endif
end block

print '("COSMO: cavity built (",A,") - ",I0," atoms, ",I0," surface tesserae, eps=",F8.2)', &
      trim(cavity_type), natoms, n_tess, epsilon
f_eps = (epsilon-1.0d0)/(epsilon+0.5d0)

allocate(V_nuc(n_tess))
do k = 1,n_tess
   V_nuc(k) = 0.0d0
   do ia = 1,natoms
      Zeff = real(Z(ia) - ecpCoreElec(ia),8)
      d = sqrt(sum((tess_coor(:,k)-coor_bohr(:,ia))**2))
      V_nuc(k) = V_nuc(k) + Zeff/d
   enddo
enddo

allocate(A_fac(n_tess,n_tess))
!$omp parallel do private(l) schedule(static)
do k = 1,n_tess
   do l = 1,n_tess
      if (k == l) then
         A_fac(k,l) = 1.07d0*sqrt(4.0d0*PI_cosmo/tess_area(k))
      else
         A_fac(k,l) = 1.0d0/sqrt(sum((tess_coor(:,k)-tess_coor(:,l))**2))
      endif
   enddo
enddo
!$omp end parallel do

allocate(A_ipiv(n_tess))
call dsytrf('U', n_tess, A_fac, n_tess, A_ipiv, wsize, -1, info)
LWORK = int(wsize(1))
allocate(WORK(LWORK))
call dsytrf('U', n_tess, A_fac, n_tess, A_ipiv, WORK, LWORK, info)
deallocate(WORK)
if (info /= 0) then
   print '("COSMO: WARNING - dsytrf factorization info=",I0," (cavity matrix may be singular)")', info
endif

allocate(saved_Ainv1(n_tess))
saved_Ainv1 = 1.0d0
call dsytrs('U', n_tess, 1, A_fac, n_tess, A_ipiv, saved_Ainv1, n_tess, info)
saved_sum_Ainv1 = sum(saved_Ainv1)

nConts_saved = nConts_in
store_mode = (int(n_tess,8)*int(nConts_in,8)*int(nConts_in,8)*8_8 <= get_store_threshold_bytes())
print '("COSMO: tessera-ESP-matrix mode = ",A)', merge("STORE ","DIRECT",store_mode)
if (store_mode) then
   allocate(B_store(nConts_in,nConts_in,n_tess))
   !$omp parallel do schedule(dynamic)
   do k = 1,n_tess
      call cosmo_build_one_tess_matrix(tess_coor(:,k), B_store(:,:,k))
   enddo
   !$omp end parallel do
endif

allocate(q_last(n_tess))
q_last = 0.0d0
cosmo_enabled = .true.
deallocate(coor_bohr)
end subroutine cosmo_init

end module mod_cosmo_init
