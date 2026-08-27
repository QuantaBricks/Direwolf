! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! COSX (chain-of-spheres) exchange matrix build: numerical-integration grid setup plus Fock exchange contraction.

submodule (mod_integrals) exchange_cosx_impl
implicit none

integer,parameter :: BATCH = 200

interface
   subroutine gridgen_nlc(nrad, nsph, npts, coor_out, weight_out, coor_in, per_atom_period_scale, intacc_eps)
   integer,intent(in) :: nrad, nsph
   integer,intent(out) :: npts
   real(8),allocatable,intent(out) :: coor_out(:,:), weight_out(:)
   real(8),optional,intent(in) :: coor_in(3,*)
   logical,optional,intent(in) :: per_atom_period_scale
   real(8),optional,intent(in) :: intacc_eps
   end subroutine gridgen_nlc
end interface

interface
   recursive subroutine mergesort_morton(key, idx, key_tmp, idx_tmp, lo, hi)
   integer(8),intent(inout) :: key(:)
   integer,intent(inout) :: idx(:)
   integer(8),intent(inout) :: key_tmp(:)
   integer,intent(inout) :: idx_tmp(:)
   integer,intent(in) :: lo,hi
   end subroutine mergesort_morton
end interface

logical,save :: cosx_ready = .false.
integer(8),save :: dbg_full_calls = 0, dbg_lr_calls = 0
integer(8),save :: dbg_full_ticks = 0_8, dbg_lr_ticks = 0_8
logical,save :: cosx_overlap_fit = .true.
integer,save :: cosx_nBases_cached = -1
real(8),save :: cosx_blkrad_max = 1.0d30
real(8),save :: cosx_kscreen = 1.0d-7
real(8),save :: cosx_kscreen_incr = 1.0d-9
logical,save :: cosx_no_incremental = .false.
logical,save :: cosx_no_incremental_lr = .false.
logical,save :: cosx_no_incremental_sr = .false.
logical,save :: cosx_peratom_lomd = .false.
real(8),save :: cosx_ext_threshold = 1.0d-6
real(8),save :: cosx_dp_skip = 1.0d-6
integer,allocatable,save :: cosx_shell_atom(:), cosx_shell_local(:)
integer,allocatable,save :: cosx_ao_offset(:), cosx_shell_dim(:)
real(8),allocatable,save,target :: cosx_coor_lo(:,:), cosx_w_lo(:), cosx_Qfull_lo(:,:)
real(8),allocatable,save,target :: cosx_coor_md(:,:), cosx_w_md(:), cosx_Qfull_md(:,:)
real(8),allocatable,save,target :: cosx_coor_hi(:,:), cosx_w_hi(:), cosx_Qfull_hi(:,:)
real(8),allocatable,save,target :: cosx_QfullT_lo(:,:), cosx_QfullT_md(:,:), cosx_QfullT_hi(:,:)
integer,save :: cosx_npts_lo = 0, cosx_npts_md = 0, cosx_npts_hi = 0
integer,allocatable,save,target :: cosx_blkptr_lo(:), cosx_blkidx_lo(:)
integer,allocatable,save,target :: cosx_blkptr_md(:), cosx_blkidx_md(:)
integer,allocatable,save,target :: cosx_blkptr_hi(:), cosx_blkidx_hi(:)
integer,allocatable,save,target :: cosx_blkstart_lo(:), cosx_blkstart_md(:), cosx_blkstart_hi(:)
integer,save :: cosx_nblk_lo = 0, cosx_nblk_md = 0, cosx_nblk_hi = 0
real(8),allocatable,save,target :: cosx_blkcen_lo(:,:), cosx_blkrad_lo(:)
real(8),allocatable,save,target :: cosx_blkcen_md(:,:), cosx_blkrad_md(:)
real(8),allocatable,save,target :: cosx_blkcen_hi(:,:), cosx_blkrad_hi(:)
real(8),allocatable,save :: cosx_shell_cen(:,:), cosx_shell_extent(:)
integer,save :: cosx_call_count = 0
logical,save :: cosx_incr_valid = .false.
integer,save :: cosx_stage_prev = 0
integer,save :: cosx_stage_settle = 0
integer,parameter :: COSX_REBASE_PERIOD = 3
integer,save :: cosx_calls_since_full = 0
real(8),allocatable,save :: cosx_D_incr_a(:,:), cosx_D_incr_b(:,:)
real(8),allocatable,save :: cosx_Ka_prev(:,:), cosx_Kb_prev(:,:)
real(8),allocatable,save :: cosx_D_incr_a_lr(:,:), cosx_D_incr_b_lr(:,:)
real(8),allocatable,save :: cosx_Ka_prev_lr(:,:), cosx_Kb_prev_lr(:,:)
logical,save :: cosx_incr_valid_lr = .false.
integer,save :: cosx_stage_prev_lr = 0
integer,save :: cosx_stage_settle_lr = 0
real(8),allocatable,save :: cosx_D_incr_a_sr(:,:), cosx_D_incr_b_sr(:,:)
real(8),allocatable,save :: cosx_Ka_prev_sr(:,:), cosx_Kb_prev_sr(:,:)
logical,save :: cosx_incr_valid_sr = .false.
integer,save :: cosx_stage_prev_sr = 0
integer,save :: cosx_stage_settle_sr = 0
real(8),allocatable,save,target :: cosx_esp_bound(:,:)
real(8),save :: cosx_esp_bound_max = 0.0d0
real(8),allocatable,save,target :: cosx_esp_mono(:,:)
real(8),allocatable,save,target :: cosx_esp_bound_lr(:,:)
logical,save :: cosx_esp_bound_lr_ready = .false.
real(8),save :: cosx_esp_bound_lr_omega = -1.0d0
real(8),allocatable,save,target :: cosx_esp_bound_sr(:,:)
logical,save :: cosx_esp_bound_sr_ready = .false.
real(8),save :: cosx_esp_bound_sr_omega = -1.0d0
real(8),save :: cosx_esp_bound_lr_max = 0.0d0
real(8),save :: cosx_kscreen_lr = 1.0d-7
real(8),save :: cosx_esp_bound_sr_max = 0.0d0
real(8),save :: cosx_kscreen_sr = 1.0d-7
real(8),save :: cosx_kscreen_incr_lr = 1.0d-9
real(8),save :: cosx_kscreen_incr_sr = 1.0d-9
real(8),allocatable,save,target :: cosx_Xblock_lo(:), cosx_Xblock_md(:), cosx_Xblock_hi(:)
integer,allocatable,save,target :: cosx_ext_ptr(:), cosx_ext_idx(:)
integer,allocatable,save :: cosx_ext_ptr_sym(:), cosx_ext_idx_sym(:)
integer,allocatable,save,target :: cosx_ext_ptr_lr(:), cosx_ext_idx_lr(:)
integer,allocatable,save :: cosx_ext_ptr_lr_sym(:), cosx_ext_idx_lr_sym(:)
integer,allocatable,save,target :: cosx_ext_ptr_sr(:), cosx_ext_idx_sr(:)
integer,allocatable,save :: cosx_ext_ptr_sr_sym(:), cosx_ext_idx_sr_sym(:)
real(8),allocatable,save :: cosx_G_hi_a(:,:), cosx_G_hi_b(:,:)
logical,save :: cosx_G_hi_valid = .false.
real(8),allocatable,save :: cosx_G_hi_a_lr(:,:), cosx_G_hi_b_lr(:,:)
logical,save :: cosx_G_hi_valid_lr = .false.

contains

subroutine cosx_sort_grid_morton(n, coor, w)
implicit none
integer,intent(in) :: n
real(8),intent(inout) :: coor(3,n), w(n)
integer(8),allocatable :: morton(:), morton_tmp(:)
integer,allocatable :: idx(:), idx_tmp(:)
real(8),allocatable :: coor_sorted(:,:), w_sorted(:)
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
integer(8) :: qx,qy,qz
integer,parameter :: QBITS = 20
integer(8),parameter :: QMAX = 2_8**QBITS - 1
integer :: i,b

allocate(morton(n),idx(n))
xmin = minval(coor(1,1:n)); xmax = maxval(coor(1,1:n))
ymin = minval(coor(2,1:n)); ymax = maxval(coor(2,1:n))
zmin = minval(coor(3,1:n)); zmax = maxval(coor(3,1:n))
!$omp parallel do private(i,qx,qy,qz,b)
do i = 1,n
   if (xmax .gt. xmin) then
      qx = int((coor(1,i)-xmin)/(xmax-xmin) * QMAX, 8)
   else
      qx = 0_8
   endif
   if (ymax .gt. ymin) then
      qy = int((coor(2,i)-ymin)/(ymax-ymin) * QMAX, 8)
   else
      qy = 0_8
   endif
   if (zmax .gt. zmin) then
      qz = int((coor(3,i)-zmin)/(zmax-zmin) * QMAX, 8)
   else
      qz = 0_8
   endif
   morton(i) = 0_8
   do b = 0,QBITS-1
      if (btest(qx,b)) morton(i) = ibset(morton(i), 3*b)
      if (btest(qy,b)) morton(i) = ibset(morton(i), 3*b+1)
      if (btest(qz,b)) morton(i) = ibset(morton(i), 3*b+2)
   enddo
   idx(i) = i
enddo
!$omp end parallel do

allocate(morton_tmp(n),idx_tmp(n))
call mergesort_morton(morton, idx, morton_tmp, idx_tmp, 1, n)
deallocate(morton_tmp,idx_tmp,morton)

allocate(coor_sorted(3,n), w_sorted(n))
do i = 1,n
   coor_sorted(:,i) = coor(:,idx(i))
   w_sorted(i) = w(idx(i))
enddo
coor(:,1:n) = coor_sorted(:,1:n)
w(1:n) = w_sorted(1:n)
deallocate(coor_sorted,w_sorted,idx)

end subroutine cosx_sort_grid_morton

subroutine cosx_dgemm(ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
implicit none
character,intent(in) :: ta, tb
integer,intent(in) :: m, n, k, lda, ldb, ldc
real(8),intent(in) :: alpha, beta
real(8),intent(in) :: A(lda,*), B(ldb,*)
real(8),intent(inout) :: C(ldc,*)
call dgemm(ta,tb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
end subroutine cosx_dgemm

subroutine cosx_build_ext_lr(nBases_in)
implicit none
integer,intent(in) :: nBases_in
integer :: si_e, sj_e, nnz_e, idx_e

if (allocated(cosx_ext_ptr_lr)) deallocate(cosx_ext_ptr_lr)
if (allocated(cosx_ext_idx_lr)) deallocate(cosx_ext_idx_lr)
allocate(cosx_ext_ptr_lr(0:nBases_in))
cosx_ext_ptr_lr(0) = 0
do si_e = 0,nBases_in-1
   nnz_e = 0
   do sj_e = si_e,nBases_in-1
      if (cosx_esp_bound_lr(si_e,sj_e) .gt. cosx_ext_threshold) nnz_e = nnz_e + 1
   enddo
   cosx_ext_ptr_lr(si_e+1) = cosx_ext_ptr_lr(si_e) + nnz_e
enddo
allocate(cosx_ext_idx_lr(cosx_ext_ptr_lr(nBases_in)))
idx_e = 0
do si_e = 0,nBases_in-1
   do sj_e = si_e,nBases_in-1
      if (cosx_esp_bound_lr(si_e,sj_e) .gt. cosx_ext_threshold) then
         idx_e = idx_e + 1
         cosx_ext_idx_lr(idx_e) = sj_e
      endif
   enddo
enddo

if (allocated(cosx_ext_ptr_lr_sym)) deallocate(cosx_ext_ptr_lr_sym)
if (allocated(cosx_ext_idx_lr_sym)) deallocate(cosx_ext_idx_lr_sym)
allocate(cosx_ext_ptr_lr_sym(0:nBases_in))
cosx_ext_ptr_lr_sym(0) = 0
do si_e = 0,nBases_in-1
   nnz_e = 0
   do sj_e = 0,nBases_in-1
      if (cosx_esp_bound_lr(si_e,sj_e) .gt. cosx_ext_threshold) nnz_e = nnz_e + 1
   enddo
   cosx_ext_ptr_lr_sym(si_e+1) = cosx_ext_ptr_lr_sym(si_e) + nnz_e
enddo
allocate(cosx_ext_idx_lr_sym(cosx_ext_ptr_lr_sym(nBases_in)))
idx_e = 0
do si_e = 0,nBases_in-1
   do sj_e = 0,nBases_in-1
      if (cosx_esp_bound_lr(si_e,sj_e) .gt. cosx_ext_threshold) then
         idx_e = idx_e + 1
         cosx_ext_idx_lr_sym(idx_e) = sj_e
      endif
   enddo
enddo
end subroutine cosx_build_ext_lr

subroutine cosx_build_ext_sr(nBases_in)
implicit none
integer,intent(in) :: nBases_in
integer :: si_e, sj_e, nnz_e, idx_e

if (allocated(cosx_ext_ptr_sr)) deallocate(cosx_ext_ptr_sr)
if (allocated(cosx_ext_idx_sr)) deallocate(cosx_ext_idx_sr)
allocate(cosx_ext_ptr_sr(0:nBases_in))
cosx_ext_ptr_sr(0) = 0
do si_e = 0,nBases_in-1
   nnz_e = 0
   do sj_e = si_e,nBases_in-1
      if (cosx_esp_bound_sr(si_e,sj_e) .gt. cosx_ext_threshold) nnz_e = nnz_e + 1
   enddo
   cosx_ext_ptr_sr(si_e+1) = cosx_ext_ptr_sr(si_e) + nnz_e
enddo
allocate(cosx_ext_idx_sr(cosx_ext_ptr_sr(nBases_in)))
idx_e = 0
do si_e = 0,nBases_in-1
   do sj_e = si_e,nBases_in-1
      if (cosx_esp_bound_sr(si_e,sj_e) .gt. cosx_ext_threshold) then
         idx_e = idx_e + 1
         cosx_ext_idx_sr(idx_e) = sj_e
      endif
   enddo
enddo

if (allocated(cosx_ext_ptr_sr_sym)) deallocate(cosx_ext_ptr_sr_sym)
if (allocated(cosx_ext_idx_sr_sym)) deallocate(cosx_ext_idx_sr_sym)
allocate(cosx_ext_ptr_sr_sym(0:nBases_in))
cosx_ext_ptr_sr_sym(0) = 0
do si_e = 0,nBases_in-1
   nnz_e = 0
   do sj_e = 0,nBases_in-1
      if (cosx_esp_bound_sr(si_e,sj_e) .gt. cosx_ext_threshold) nnz_e = nnz_e + 1
   enddo
   cosx_ext_ptr_sr_sym(si_e+1) = cosx_ext_ptr_sr_sym(si_e) + nnz_e
enddo
allocate(cosx_ext_idx_sr_sym(cosx_ext_ptr_sr_sym(nBases_in)))
idx_e = 0
do si_e = 0,nBases_in-1
   do sj_e = 0,nBases_in-1
      if (cosx_esp_bound_sr(si_e,sj_e) .gt. cosx_ext_threshold) then
         idx_e = idx_e + 1
         cosx_ext_idx_sr_sym(idx_e) = sj_e
      endif
   enddo
enddo
end subroutine cosx_build_ext_sr

subroutine cosx_build_shared(nConts)
use MOL_info, only: atoms, nAtoms
use GRID_info, only: rcut2_shared, maxShell_shared, maxGauss_shared
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: nConts
integer :: si, sj, i, j
character(len=16) :: envchar
integer :: envstat
real(8),parameter :: SCREEN_EPS = 1.0d-7, SCREEN_MARGIN = 2.0d0
integer :: pi_,pj_,ni,nj,offi,offj,nnz,idx
real(8) :: ei,ej,ci,cj,r2,dist
real(8) :: psum,kab

if (.not. allocated(rcut2_shared)) then
   maxShell_shared = 0
   maxGauss_shared = 0
   do i = 1,nAtoms
      maxShell_shared = max(maxShell_shared,size(atoms(i)%shell))
      do j = 1,size(atoms(i)%shell)
         maxGauss_shared = max(maxGauss_shared,atoms(i)%shell(j)%nGauss)
      enddo
   enddo
   allocate(rcut2_shared(maxShell_shared,nAtoms))
   rcut2_shared = 0
   do i = 1,nAtoms
      do j = 1,size(atoms(i)%shell)
         rcut2_shared(j,i) = (sqrt(-log(SCREEN_EPS)/minval(atoms(i)%shell(j)%exponents)) + SCREEN_MARGIN)**2
      enddo
   enddo
endif

call get_environment_variable("ENGINE_COSX_OVERLAP_FIT", envchar, status=envstat)
if (envstat .eq. 0) then
   cosx_overlap_fit = .not. (trim(envchar).eq.'0' .or. trim(envchar).eq.'false')
endif
if (allocated(cosx_shell_atom)) deallocate(cosx_shell_atom)
if (allocated(cosx_shell_local)) deallocate(cosx_shell_local)
if (allocated(cosx_ao_offset)) deallocate(cosx_ao_offset)
if (allocated(cosx_shell_dim)) deallocate(cosx_shell_dim)
allocate(cosx_shell_atom(0:nBases-1), cosx_shell_local(0:nBases-1))
allocate(cosx_ao_offset(0:nBases-1), cosx_shell_dim(0:nBases-1))

si = 0
do i = 1,nAtoms
   do j = 1,size(atoms(i)%shell)
      cosx_shell_atom(si) = i
      cosx_shell_local(si) = j
      si = si + 1
   enddo
enddo
if (si .ne. nBases) then
   print *, "COSX: shell traversal mismatch (atoms/shell count ",si, &
            ") vs nBases (",nBases,") - grid<->AO mapping would be wrong"
   call flush(6)
   stop 1
endif

cosx_ao_offset(0) = 0
do si = 0,nBases-1
   cosx_shell_dim(si) = cgto_engine(si, bas)
   if (si .lt. nBases-1) cosx_ao_offset(si+1) = cosx_ao_offset(si) + cosx_shell_dim(si)
enddo

cosx_nBases_cached = nBases

call get_environment_variable("ENGINE_COSX_BLKRAD", envchar, status=envstat)
if (envstat .eq. 0) read(envchar,*) cosx_blkrad_max
call get_environment_variable("ENGINE_COSX_KSCREEN", envchar, status=envstat)
if (envstat .eq. 0) read(envchar,*) cosx_kscreen
call get_environment_variable("ENGINE_COSX_KSCREEN_INCR", envchar, status=envstat)
if (envstat .eq. 0) read(envchar,*) cosx_kscreen_incr
call get_environment_variable("ENGINE_COSX_NO_INCREMENTAL", envchar, status=envstat)
if (envstat .eq. 0) cosx_no_incremental = (trim(envchar).eq.'1' .or. trim(envchar).eq.'true')
call get_environment_variable("ENGINE_COSX_NO_INCREMENTAL_LR", envchar, status=envstat)
if (envstat .eq. 0) cosx_no_incremental_lr = (trim(envchar).eq.'1' .or. trim(envchar).eq.'true')
call get_environment_variable("ENGINE_COSX_PERATOM_LOMD", envchar, status=envstat)
if (envstat .eq. 0) cosx_peratom_lomd = (trim(envchar).eq.'1' .or. trim(envchar).eq.'true')
call get_environment_variable("ENGINE_COSX_NO_INCREMENTAL_SR", envchar, status=envstat)
if (envstat .eq. 0) cosx_no_incremental_sr = (trim(envchar).eq.'1' .or. trim(envchar).eq.'true')
call get_environment_variable("ENGINE_COSX_DP_SKIP", envchar, status=envstat)
if (envstat .eq. 0) read(envchar,*) cosx_dp_skip
call get_environment_variable("ENGINE_COSX_EXT_THRESHOLD", envchar, status=envstat)
if (envstat .eq. 0) read(envchar,*) cosx_ext_threshold

if (allocated(cosx_esp_bound)) deallocate(cosx_esp_bound)
if (allocated(cosx_esp_mono)) deallocate(cosx_esp_mono)
if (allocated(cosx_ext_ptr)) deallocate(cosx_ext_ptr)
if (allocated(cosx_ext_idx)) deallocate(cosx_ext_idx)
allocate(cosx_esp_bound(0:nBases-1,0:nBases-1))
allocate(cosx_esp_mono(0:nBases-1,0:nBases-1))
if (allocated(cosx_shell_cen)) deallocate(cosx_shell_cen,cosx_shell_extent)
allocate(cosx_shell_cen(3,0:nBases-1), cosx_shell_extent(0:nBases-1))
do si = 0,nBases-1
   cosx_shell_cen(:,si) = atoms(cosx_shell_atom(si))%coor*ans2bohr
   cosx_shell_extent(si) = sqrt(rcut2_shared(cosx_shell_local(si),cosx_shell_atom(si)))
enddo

cosx_esp_bound = 0.0d0
cosx_esp_mono = 0.0d0
!$omp parallel do default(shared) private(si,sj,ni,nj,offi,offj,r2,pi_,pj_,ei,ci,ej,cj,psum,kab) schedule(dynamic)
do si = 0,nBases-1
   ni = bas(3,si+1)
   offi = bas(6,si+1)
   do sj = si,nBases-1
      nj = bas(3,sj+1)
      offj = bas(6,sj+1)
      r2 = sum((cosx_shell_cen(:,si)-cosx_shell_cen(:,sj))**2)
      do pi_ = 1,ni
         ei = env(offi+pi_)
         ci = env(bas(7,si+1)+pi_)
         do pj_ = 1,nj
            ej = env(offj+pj_)
            cj = env(bas(7,sj+1)+pj_)
            psum = ei + ej
            kab = abs(ci*cj)*exp(-r2*ei*ej/psum)
            cosx_esp_bound(si,sj) = cosx_esp_bound(si,sj) + &
               kab*2.0d0*acos(-1.0d0)/psum
            cosx_esp_mono(si,sj) = cosx_esp_mono(si,sj) + &
               kab*(acos(-1.0d0)/psum)**1.5d0
         enddo
      enddo
      cosx_esp_bound(sj,si) = cosx_esp_bound(si,sj)
      cosx_esp_mono(sj,si) = cosx_esp_mono(si,sj)
   enddo
enddo
!$omp end parallel do
cosx_esp_bound_max = maxval(cosx_esp_bound)

allocate(cosx_ext_ptr(0:nBases))
cosx_ext_ptr(0) = 0
do si = 0,nBases-1
   nnz = 0
   do sj = si,nBases-1
      if (cosx_esp_bound(si,sj) .gt. cosx_ext_threshold) nnz = nnz + 1
   enddo
   cosx_ext_ptr(si+1) = cosx_ext_ptr(si) + nnz
enddo
allocate(cosx_ext_idx(cosx_ext_ptr(nBases)))
idx = 0
do si = 0,nBases-1
   do sj = si,nBases-1
      if (cosx_esp_bound(si,sj) .gt. cosx_ext_threshold) then
         idx = idx + 1
         cosx_ext_idx(idx) = sj
      endif
   enddo
enddo

if (allocated(cosx_ext_ptr_sym)) deallocate(cosx_ext_ptr_sym)
if (allocated(cosx_ext_idx_sym)) deallocate(cosx_ext_idx_sym)
allocate(cosx_ext_ptr_sym(0:nBases))
cosx_ext_ptr_sym(0) = 0
do si = 0,nBases-1
   nnz = 0
   do sj = 0,nBases-1
      if (cosx_esp_bound(si,sj) .gt. cosx_ext_threshold) nnz = nnz + 1
   enddo
   cosx_ext_ptr_sym(si+1) = cosx_ext_ptr_sym(si) + nnz
enddo
allocate(cosx_ext_idx_sym(cosx_ext_ptr_sym(nBases)))
idx = 0
do si = 0,nBases-1
   do sj = 0,nBases-1
      if (cosx_esp_bound(si,sj) .gt. cosx_ext_threshold) then
         idx = idx + 1
         cosx_ext_idx_sym(idx) = sj
      endif
   enddo
enddo

end subroutine cosx_build_shared

subroutine cosx_build_grid_distprune(nrad, nsph, intacc_eps, coor_out, w_out, npts_out)
use MOL_info, only: atoms, nAtoms, engine_verbose
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: nrad, nsph
real(8),intent(in) :: intacc_eps
real(8),allocatable,intent(out) :: coor_out(:,:), w_out(:)
integer,intent(out) :: npts_out
integer :: nrad_atom(nAtoms), row_i, iatm, i, j, Ntemp, ii, jj, label
integer :: atom_pts_i(nAtoms), atom_off(nAtoms), cursphpot
real(8) :: Rij(nAtoms,nAtoms), aij_mat(nAtoms,nAtoms), acoor(3,nAtoms)
real(8) :: parm, ratm, radx, radr, radw, ratio, chi, uij, aij
real(8) :: potx(974),poty(974),potz(974),potw(974)
real(8) :: rdist(nAtoms), Pvec(nAtoms), rmiu, tmps, sij

do iatm = 1,nAtoms
   acoor(:,iatm) = atoms(iatm)%coor*ans2bohr
   if (atoms(iatm)%charge .le. 2) then
      row_i = 1
   else if (atoms(iatm)%charge .le. 10) then
      row_i = 2
   else if (atoms(iatm)%charge .le. 18) then
      row_i = 3
   else if (atoms(iatm)%charge .le. 36) then
      row_i = 4
   else if (atoms(iatm)%charge .le. 54) then
      row_i = 5
   else if (atoms(iatm)%charge .le. 86) then
      row_i = 6
   else
      row_i = 7
   endif
   nrad_atom(iatm) = max(1, int(intacc_eps*15.0d0 - 40.0d0 + 5.0d0*row_i))
   if (cosx_peratom_lomd .and. row_i .ge. 3) nrad_atom(iatm) = nrad_atom(iatm)*3
enddo

npts_out = 0
do iatm = 1,nAtoms
   atom_pts_i(iatm) = 0
   parm = covrad(atoms(iatm)%charge)/2*ans2bohr
   ratm = covrad(atoms(iatm)%charge)*ans2bohr
   do i = 1,nrad_atom(iatm)
      radx = cos(i*PI/(nrad_atom(iatm)+1))
      radr = (1+radx)/(1-radx)*parm
      ratio = radr/ratm
      cursphpot = cosx_sgx_prune(nsph, atoms(iatm)%charge, ratio)
      atom_pts_i(iatm) = atom_pts_i(iatm) + cursphpot
   enddo
   atom_off(iatm) = npts_out
   npts_out = npts_out + atom_pts_i(iatm)
enddo
if (engine_verbose .ge. 2) &
   print '(A,F6.3,A,I0)', "  COSX distprune gen nrad(row-based, eps=",intacc_eps,") npts=",npts_out
call flush(6)
allocate(coor_out(3,npts_out))
allocate(w_out(npts_out))

Rij = 0.0d0
aij_mat = 0.0d0
do ii = 1,nAtoms
   do jj = 1,nAtoms
      if (ii .eq. jj) cycle
      Rij(ii,jj) = dsqrt(sum((acoor(:,ii) - acoor(:,jj))**2))
      chi = covrad(atoms(ii)%charge)/covrad(atoms(jj)%charge)
      uij = (chi-1)/(chi+1)
      aij = uij/(uij**2 -1)
      if (aij .gt. 0.5d0) aij = 0.5d0
      if (aij .lt. -0.5d0) aij = -0.5d0
      aij_mat(ii,jj) = aij
   enddo
enddo

!$omp parallel do default(shared) private(iatm,parm,ratm,i,cursphpot,potx,poty,potz,potw,Ntemp) &
!$omp&   private(radx,radr,radw,ratio,j,label,ii,jj,rdist,Pvec,rmiu,tmps,sij) schedule(dynamic)
do iatm = 1,nAtoms
   parm = covrad(atoms(iatm)%charge)/2*ans2bohr
   ratm = covrad(atoms(iatm)%charge)*ans2bohr
   label = atom_off(iatm)
   do i = 1,nrad_atom(iatm)
      radx = cos(i*PI/(nrad_atom(iatm)+1))
      radr = (1+radx)/(1-radx)*parm
      ratio = radr/ratm
      cursphpot = cosx_sgx_prune(nsph, atoms(iatm)%charge, ratio)
      call cosx_ld_dispatch(cursphpot, potx,poty,potz,potw, Ntemp)
      radw = 2*PI/(nrad_atom(iatm)+1)*parm**3*(1+radx)**2.5D0/(1-radx)**3.5D0*4*PI
      do j = 1,Ntemp
         label = label + 1
         coor_out(1,label) = radr*potx(j) + acoor(1,iatm)
         coor_out(2,label) = radr*poty(j) + acoor(2,iatm)
         coor_out(3,label) = radr*potz(j) + acoor(3,iatm)
         do ii = 1,nAtoms
            rdist(ii) = dsqrt(sum((coor_out(:,label) - acoor(:,ii))**2))
         enddo
         Pvec = 1.0d0
         do ii = 1,nAtoms-1
            do jj = ii+1,nAtoms
               rmiu = (rdist(ii)-rdist(jj))/Rij(ii,jj)
               rmiu = rmiu + aij_mat(ii,jj)*(1-rmiu**2)
               tmps = 1.5d0*rmiu - 0.5d0*rmiu**3
               tmps = 1.5d0*tmps  - 0.5d0*tmps**3
               tmps = 1.5d0*tmps  - 0.5d0*tmps**3
               sij = 0.5d0*(1-tmps)
               Pvec(ii) = Pvec(ii)*sij
               Pvec(jj) = Pvec(jj)*(1.0d0-sij)
            enddo
         enddo
         w_out(label) = radw*potw(j)*Pvec(iatm)/sum(Pvec)
      enddo
   enddo
enddo
!$omp end parallel do

end subroutine cosx_build_grid_distprune

integer function cosx_sgx_prune(nsph, charge, ratio) result(cursphpot)
implicit none
integer,intent(in) :: nsph, charge
real(8),intent(in) :: ratio
integer :: level, place, k
integer :: lvls(8), leb(8,5)
real(8) :: alphas(4)
lvls = (/50,110,194,302,434,590,770,770/)
leb(1,:) = (/14, 26, 50, 50, 26/)
leb(2,:) = (/14, 26, 50,110, 50/)
leb(3,:) = (/26, 50,110,194,110/)
leb(4,:) = (/26,110,194,302,194/)
leb(5,:) = (/26,194,302,434,302/)
leb(6,:) = (/50,302,434,590,434/)
leb(7,:) = (/110,434,590,770,590/)
leb(8,:) = (/770,770,770,770,770/)
level = 1
do k = 1,8
   if (nsph .gt. lvls(k)) level = k+1
enddo
level = min(level,8)
if (charge .le. 2) then
   alphas = (/0.25d0, 0.5d0, 1.0d0, 4.5d0/)
else if (charge .le. 10) then
   alphas = (/0.1667d0, 0.5d0, 0.9d0, 3.5d0/)
else
   alphas = (/0.1d0, 0.4d0, 0.8d0, 2.5d0/)
endif
place = 1
do k = 1,4
   if (ratio .gt. alphas(k)) place = k+1
enddo
cursphpot = leb(level,place)
end function cosx_sgx_prune

subroutine cosx_ld_dispatch(nsph_in, px,py,pz,pw, nt)
implicit none
integer,intent(in) :: nsph_in
real(8),intent(out) :: px(974),py(974),pz(974),pw(974)
integer,intent(out) :: nt
select case (nsph_in)
case (6);   call LD0006(px,py,pz,pw,nt)
case (14);  call LD0014(px,py,pz,pw,nt)
case (26);  call LD0026(px,py,pz,pw,nt)
case (38);  call LD0038(px,py,pz,pw,nt)
case (50);  call LD0050(px,py,pz,pw,nt)
case (74);  call LD0074(px,py,pz,pw,nt)
case (86);  call LD0086(px,py,pz,pw,nt)
case (110); call LD0110(px,py,pz,pw,nt)
case (146); call LD0146(px,py,pz,pw,nt)
case (170); call LD0170(px,py,pz,pw,nt)
case (194); call LD0194(px,py,pz,pw,nt)
case (230); call LD0230(px,py,pz,pw,nt)
case (266); call LD0266(px,py,pz,pw,nt)
case (302); call LD0302(px,py,pz,pw,nt)
case (350); call LD0350(px,py,pz,pw,nt)
case (434); call LD0434(px,py,pz,pw,nt)
case (590); call LD0590(px,py,pz,pw,nt)
case (770); call LD0770(px,py,pz,pw,nt)
case (974); call LD0974(px,py,pz,pw,nt)
case default
   print *, "cosx_ld_dispatch: unsupported spherical order", nsph_in, "- using 110"
   call LD0110(px,py,pz,pw,nt)
end select
end subroutine cosx_ld_dispatch

subroutine cosx_malloc_trim()
use iso_c_binding, only: c_int
implicit none
interface
   function c_malloc_trim(pad) bind(C, name="malloc_trim") result(ret)
   import :: c_int
   integer(c_int),value :: pad
   integer(c_int) :: ret
   end function c_malloc_trim
end interface
integer(c_int) :: ret
ret = c_malloc_trim(0_c_int)
end subroutine cosx_malloc_trim

subroutine cosx_build_one_grid(nConts, nrad, nsph, coor_out, w_out, npts_out, Qfull_out, &
                                blkptr_out, blkidx_out, blkcen_out, blkrad_out, &
                                QfullT_out, blkstart_out, nblk_out, per_atom_period_scale, intacc_eps)
use MOL_info, only: S, engine_verbose
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, nrad, nsph
real(8),allocatable,intent(out) :: coor_out(:,:), w_out(:)
integer,intent(out) :: npts_out
real(8),allocatable,intent(out) :: Qfull_out(:,:)
integer,allocatable,intent(out) :: blkptr_out(:), blkidx_out(:)
real(8),allocatable,intent(out) :: blkcen_out(:,:), blkrad_out(:)
real(8),allocatable,intent(out) :: QfullT_out(:,:)
integer,allocatable,intent(out) :: blkstart_out(:)
integer,intent(out) :: nblk_out
logical,optional,intent(in) :: per_atom_period_scale
real(8),optional,intent(in) :: intacc_eps
integer :: ib, nblocks, bstart, nb, g, si, kk
real(8),allocatable :: coor_batch(:,:), w_batch(:), Xb(:,:), Xbs(:,:), val1_dummy(:,:,:)
real(8),allocatable :: S_num(:,:), S_num_local(:,:)
integer,allocatable :: ipiv(:)
integer :: info
real(8) :: blk_cen(3), blk_rad, dist
integer :: nnz, idx
real(8) :: bbmin(3), bbmax(3), cand_min(3), cand_max(3), halfdiag
integer :: ip, cur_start, cur_n
logical,save :: cosx_distprune = .true.
logical,save :: cosx_distprune_checked = .false.
character(len=8) :: cosx_distprune_env

if (.not. cosx_distprune_checked) then
   cosx_distprune_checked = .true.
   call get_environment_variable("ENGINE_COSX_DISTPRUNE", cosx_distprune_env)
   if (len_trim(cosx_distprune_env) .gt. 0) cosx_distprune = (trim(cosx_distprune_env) .eq. "1")
   if (cosx_distprune .and. engine_verbose .ge. 2) &
      print *,"COSX: using distance-ratio angular pruning (isolated, COSX-only)"
endif

if (cosx_distprune .and. present(intacc_eps)) then
   call cosx_build_grid_distprune(nrad, nsph, intacc_eps, coor_out, w_out, npts_out)
else
   call gridgen_nlc(nrad, nsph, npts_out, coor_out, w_out, per_atom_period_scale=per_atom_period_scale, intacc_eps=intacc_eps)
endif
call cosx_sort_grid_morton(npts_out, coor_out, w_out)

allocate(blkstart_out(0:npts_out))
nblk_out = 0
cur_start = 1
bbmin = coor_out(:,1); bbmax = coor_out(:,1)
cur_n = 1
do ip = 2,npts_out
   cand_min = min(bbmin, coor_out(:,ip))
   cand_max = max(bbmax, coor_out(:,ip))
   halfdiag = 0.5d0*sqrt(sum((cand_max-cand_min)**2))
   if (cur_n .ge. BATCH .or. halfdiag .gt. cosx_blkrad_max) then
      blkstart_out(nblk_out) = cur_start - 1
      nblk_out = nblk_out + 1
      cur_start = ip
      bbmin = coor_out(:,ip); bbmax = coor_out(:,ip)
      cur_n = 1
   else
      bbmin = cand_min; bbmax = cand_max
      cur_n = cur_n + 1
   endif
enddo
blkstart_out(nblk_out) = cur_start - 1
nblk_out = nblk_out + 1
blkstart_out(nblk_out) = npts_out
nblocks = nblk_out

allocate(blkptr_out(0:nblocks))
allocate(blkcen_out(3,0:nblocks-1), blkrad_out(0:nblocks-1))
blkptr_out(0) = 0
do ib = 0,nblocks-1
   bstart = blkstart_out(ib) + 1
   nb = blkstart_out(ib+1) - blkstart_out(ib)
   blk_cen = sum(coor_out(:,bstart:bstart+nb-1), dim=2) / dble(nb)
   blk_rad = 0.0d0
   do g = bstart,bstart+nb-1
      blk_rad = max(blk_rad, sqrt(sum((coor_out(:,g)-blk_cen)**2)))
   enddo
   blkcen_out(:,ib) = blk_cen
   blkrad_out(ib) = blk_rad
   nnz = 0
   do si = 0,nBases-1
      dist = sqrt(sum((cosx_shell_cen(:,si)-blk_cen)**2))
      if (dist .le. cosx_shell_extent(si)+blk_rad) nnz = nnz + 1
   enddo
   blkptr_out(ib+1) = blkptr_out(ib) + nnz
enddo
allocate(blkidx_out(blkptr_out(nblocks)))
idx = 0
do ib = 0,nblocks-1
   blk_cen = blkcen_out(:,ib)
   blk_rad = blkrad_out(ib)
   do si = 0,nBases-1
      dist = sqrt(sum((cosx_shell_cen(:,si)-blk_cen)**2))
      if (dist .le. cosx_shell_extent(si)+blk_rad) then
         idx = idx + 1
         blkidx_out(idx) = si
      endif
   enddo
enddo

if (.not. cosx_overlap_fit) then
   allocate(Qfull_out(0,0))
   allocate(QfullT_out(0,0))
   return
endif
allocate(S_num(nConts,nConts))
S_num = 0.0d0
nblocks = nblk_out
!$omp parallel default(shared) private(ib,bstart,nb,g,si,coor_batch,w_batch,Xb,Xbs,val1_dummy,kk) &
!$omp&   private(S_num_local)
allocate(S_num_local(nConts,nConts))
S_num_local = 0.0d0
!$omp do schedule(static)
do ib = 0,nblocks-1
   bstart = blkstart_out(ib) + 1
   nb = blkstart_out(ib+1) - blkstart_out(ib)
   allocate(coor_batch(3,nb), w_batch(nb))
   coor_batch = coor_out(:,bstart:bstart+nb-1)
   w_batch = w_out(bstart:bstart+nb-1)
   allocate(Xb(nConts,nb), val1_dummy(nConts,3,nb))
   Xb = 0.0d0
   do kk = blkptr_out(ib)+1, blkptr_out(ib+1)
      si = blkidx_out(kk)
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                                coor_batch, nb, cosx_ao_offset(si), nConts, Xb, val1_dummy)
   enddo
   deallocate(val1_dummy)
   do g = 1,nb
      Xb(:,g) = Xb(:,g) * sqrt(abs(w_batch(g)))
   enddo
   allocate(Xbs(nConts,nb))
   do g = 1,nb
      Xbs(:,g) = sign(1.0d0,w_batch(g)) * Xb(:,g)
   enddo
   call dgemm('N','T',nConts,nConts,nb,1.0d0,Xbs,nConts,Xb,nConts,1.0d0,S_num_local,nConts)
   deallocate(coor_batch,w_batch,Xb,Xbs)
enddo
!$omp end do
block
   integer :: merge_tid, merge_nthreads
   merge_nthreads = omp_get_num_threads()
   do merge_tid = 0, merge_nthreads-1
      !$omp barrier
      if (omp_get_thread_num() .eq. merge_tid) then
         S_num = S_num + S_num_local
      endif
   enddo
end block
deallocate(S_num_local)
!$omp end parallel

allocate(Qfull_out(nConts,nConts))
Qfull_out = S(1:nConts,1:nConts)
allocate(ipiv(nConts))
call dgesv(nConts, nConts, S_num, nConts, ipiv, Qfull_out, nConts, info)
if (info .ne. 0) then
   print *, "COSX: overlap-fitting dgesv failed (info=",info,") - S_num", &
            " singular/ill-conditioned; falling back to unfitted X (Qfull=I)"
   call flush(6)
   Qfull_out = 0.0d0
   do si = 1,nConts
      Qfull_out(si,si) = 1.0d0
   enddo
endif
deallocate(S_num, ipiv)

allocate(QfullT_out(nConts,nConts))
QfullT_out = transpose(Qfull_out)

end subroutine cosx_build_one_grid

end submodule exchange_cosx_impl
