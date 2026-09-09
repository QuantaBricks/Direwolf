! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Density-fitting setup: auxiliary basis construction and 3-center/2-center integral assembly.

submodule (mod_integrals) df_setup_impl
implicit none
contains

module subroutine build_df_aux_basis()
use MOL_info
implicit none
integer :: i,j,k,Lmax_atom,Lcap,L,ii
real(8) :: alpha_min,alpha_max,ratio,alpha
integer :: nEnvOld, offpoint, shidx, cartdim
real(8),external :: CINTgto_norm

nBasesAux = 0
nContsAux = 0
do i = 1,nAtoms
   Lmax_atom = 0
   do j = 1,size(atoms(i)%shell)
      Lmax_atom = max(Lmax_atom, atoms(i)%shell(j)%angMoment)
   enddo
   Lcap = min(2*Lmax_atom,4)
   do L = 0,Lcap
      nBasesAux = nBasesAux + DF_NAUX_PER_L
      cartdim = (L+1)*(L+2)/2
      nContsAux = nContsAux + DF_NAUX_PER_L*cartdim
   enddo
enddo

allocate(basDF(8,nBases+nBasesAux))
nEnvOld = size(env)
allocate(envDF(nEnvOld + 2*nBasesAux))
basDF = 0
envDF = 0.0d0
basDF(1:8,1:nBases) = bas(1:8,1:nBases)
envDF(1:nEnvOld) = env(1:nEnvOld)

offpoint = nEnvOld
shidx = nBases
do i = 1,nAtoms
   Lmax_atom = 0
   alpha_min = huge(1.0d0)
   alpha_max = 0.0d0
   do j = 1,size(atoms(i)%shell)
      Lmax_atom = max(Lmax_atom, atoms(i)%shell(j)%angMoment)
      alpha_min = min(alpha_min, minval(atoms(i)%shell(j)%exponents))
      alpha_max = max(alpha_max, maxval(atoms(i)%shell(j)%exponents))
   enddo
   Lcap = min(2*Lmax_atom,4)
   do L = 0,Lcap
      if (DF_NAUX_PER_L .gt. 1) then
         ratio = (alpha_max*DF_EXP_HI/(alpha_min*DF_EXP_LO))**(1.0d0/(DF_NAUX_PER_L-1))
      else
         ratio = 1.0d0
      endif
      do k = 1,DF_NAUX_PER_L
         alpha = alpha_min*DF_EXP_LO*ratio**(k-1)
         shidx = shidx + 1
         basDF(1,shidx) = i-1
         basDF(2,shidx) = L
         basDF(3,shidx) = 1
         basDF(4,shidx) = 1
         basDF(6,shidx) = offpoint
         offpoint = offpoint + 1
         envDF(offpoint) = alpha
         basDF(7,shidx) = offpoint
         offpoint = offpoint + 1
         envDF(offpoint) = CINTgto_norm(L,alpha)
      enddo
   enddo
enddo
nBasesAux = shidx - nBases

allocate(NorVECAux(nContsAux))
NorVECAux = 0.0d0
block
   real(8),allocatable :: buf1e(:,:)
   integer :: shls(4), di
   real(8),allocatable :: NorTmp(:)
   allocate(NorTmp(nConts+nContsAux))
   NorTmp = 0.0d0
   NorTmp(1:nConts) = NorVEC(1:nConts)
   do ii = nBases,nBases+nBasesAux-1
      shls(1) = ii
      shls(2) = ii
      di = cgto_engine(ii, basDF)
      allocate(buf1e(di,di))
      call ovlp1e_engine(buf1e, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, 0_8)
      call Normal(ii, di, basDF, nBases+nBasesAux, buf1e, NorTmp, nConts+nContsAux)
      deallocate(buf1e)
   enddo
   NorVECAux(1:nContsAux) = NorTmp(nConts+1:nConts+nContsAux)
   deallocate(NorTmp)
end block
end subroutine build_df_aux_basis

module subroutine build_df_aux_basis_fromfile(auxname)
use MOL_info
use mod_basis_files, only: basis_file_for_label
implicit none
character(len=*),intent(in) :: auxname
integer :: i,k,jsh,shidx,nGauss,offpoint,nEnvOld,ios,L,numShellAtom
character(len=30) :: base_temp, basis_char
character(len=200) :: bpath
real(8) :: expo, coeff
real(8) :: coeffs(50)
real(8),external :: CINTgto_norm
integer,allocatable :: nShellAtom(:)

write(base_temp,"(I2.2,A1,A)") atoms(1)%charge, "-", auxname
bpath = basis_file_for_label(base_temp)
if (len_trim(bpath) == 0) stop "DF aux basis label not found in data/bases"
open(unit=202, file=trim(bpath))

allocate(nShellAtom(nAtoms))
nBasesAux = 0
nContsAux = 0
do i = 1,nAtoms
   write(base_temp,"(I2.2,A1,A)") atoms(i)%charge, "-", auxname
   rewind(202)
   do
      read(202,*,iostat=ios) basis_char
      if (ios .ne. 0) stop "DF aux basis label not found in data/bases"
      if (trim(basis_char) .eq. trim(base_temp)) exit
   enddo
   numShellAtom = 0
   do
      read(202,*) basis_char, basis_char
      if (trim(basis_char) .eq. "END") exit
      select case(trim(basis_char))
      case("s"); L = 0
      case("p"); L = 1
      case("d"); L = 2
      case("f"); L = 3
      case("g"); L = 4
      case("h"); L = 5
      case("i"); L = 6
      case default
         cycle
      end select
      numShellAtom = numShellAtom + 1
      nBasesAux = nBasesAux + 1
      if (engine_puream) then
         nContsAux = nContsAux + (2*L+1)
      else
         nContsAux = nContsAux + (L+1)*(L+2)/2
      endif
   enddo
   nShellAtom(i) = numShellAtom
enddo
rewind(202)

allocate(basDF(8,nBases+nBasesAux))
nEnvOld = size(env)
allocate(envDF(nEnvOld + 20*nBasesAux))
basDF = 0
envDF = 0.0d0
basDF(1:8,1:nBases) = bas(1:8,1:nBases)
envDF(1:nEnvOld) = env(1:nEnvOld)

offpoint = nEnvOld
shidx = nBases
do i = 1,nAtoms
   write(base_temp,"(I2.2,A1,A)") atoms(i)%charge, "-", auxname
   rewind(202)
   do
      read(202,*,iostat=ios) basis_char
      if (ios .ne. 0) stop "DF aux basis label not found in data/bases"
      if (trim(basis_char) .eq. trim(base_temp)) exit
   enddo
   do jsh = 1,nShellAtom(i)
      read(202,*) nGauss, basis_char
      select case(trim(basis_char))
      case("s"); L = 0
      case("p"); L = 1
      case("d"); L = 2
      case("f"); L = 3
      case("g"); L = 4
      case("h"); L = 5
      case("i"); L = 6
      end select
      shidx = shidx + 1
      basDF(1,shidx) = i-1
      basDF(2,shidx) = L
      basDF(3,shidx) = nGauss
      basDF(4,shidx) = 1
      basDF(6,shidx) = offpoint
      do k = 1,nGauss
         read(202,*) expo, coeffs(k)
         envDF(offpoint+k) = expo
      enddo
      offpoint = offpoint + nGauss
      basDF(7,shidx) = offpoint
      do k = 1,nGauss
         envDF(offpoint+k) = coeffs(k)*CINTgto_norm(L, envDF(basDF(6,shidx)+k))
      enddo
      offpoint = offpoint + nGauss
   enddo
enddo
nBasesAux = shidx - nBases
close(202)
deallocate(nShellAtom)

allocate(NorVECAux(nContsAux))
NorVECAux = 0.0d0
block
   real(8),allocatable :: buf1e(:,:)
   integer :: shls(4), di, ii
   real(8),allocatable :: NorTmp(:)
   allocate(NorTmp(nConts+nContsAux))
   NorTmp = 0.0d0
   NorTmp(1:nConts) = NorVEC(1:nConts)
   do ii = nBases,nBases+nBasesAux-1
      shls(1) = ii
      shls(2) = ii
      di = cgto_engine(ii, basDF)
      allocate(buf1e(di,di))
      call ovlp1e_engine(buf1e, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, 0_8)
      call Normal(ii, di, basDF, nBases+nBasesAux, buf1e, NorTmp, nConts+nContsAux)
      deallocate(buf1e)
   enddo
   NorVECAux(1:nContsAux) = NorTmp(nConts+1:nConts+nContsAux)
   deallocate(NorTmp)
end block
end subroutine build_df_aux_basis_fromfile

module subroutine build_df_integrals()
use MOL_info
use mod_meminfo, only: engine_avail_at_start_bytes, get_process_memory_bytes
use mod_exchange, only: HF_exchange_frac, RS_omega
use mod_profile, only: prof_start, prof_stop, fmt_gb, itoa
use omp_lib, only: omp_get_max_threads
implicit none
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offP,offQ,offi,offj,info
integer :: shls(4)
real(8),allocatable :: buf2c(:,:), buf3c(:,:,:)
real(8) :: maxdiag
real(8),allocatable :: metric(:,:), eval(:), work(:)
integer,allocatable :: iwork(:)
integer :: lwork, liwork, i
integer,allocatable :: aux_offset(:)
integer,allocatable :: pair_slot_2d(:,:)
real(8) :: val3c
integer :: slot_pq, slot_qp
integer(8) :: df_total_pairs_est
integer(8) :: diag_rss, diag_hwm

call threec2e_optimizer_engine(int3c2e_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)
call twoc2e_optimizer_engine(int2c2e_opt, atm, nAtoms, basDF, nBases+nBasesAux, envDF)

if (engine_verbose .ge. 2) then
   call get_process_memory_bytes(diag_rss, diag_hwm)
   print '(A)', "  DIAG before metric alloc: nBases="//trim(itoa(nBases))// &
           "  nBasesAux="//trim(itoa(nBasesAux))//"  nContsAux="//trim(itoa(nContsAux))// &
           "  RSS="//trim(fmt_gb(real(diag_rss,8)/1024**3))//"GB"
endif
allocate(metric(nContsAux,nContsAux))
allocate(aux_shell_bound(nBasesAux))
metric = 0.0d0
aux_shell_bound = 0.0d0
allocate(aux_offset(nBases:nBases+nBasesAux-1))
offP = 0
do ish = nBases,nBases+nBasesAux-1
   aux_offset(ish) = offP
   offP = offP + cgto_engine(ish, basDF)
enddo
!$omp parallel do private(ish,jsh,di,dj,offP,offQ,shls,buf2c,p,q,maxdiag) schedule(dynamic)
do ish = nBases,nBases+nBasesAux-1
   di = cgto_engine(ish, basDF)
   offP = aux_offset(ish)
   do jsh = nBases,ish
      dj = cgto_engine(jsh, basDF)
      offQ = aux_offset(jsh)
      shls(1) = ish
      shls(2) = jsh
      allocate(buf2c(di,dj))
      call twoc2e_engine(buf2c, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int2c2e_opt)
      do q = 1,dj
         do p = 1,di
            metric(offP+p,offQ+q) = buf2c(p,q)*NorVECAux(offP+p)*NorVECAux(offQ+q)
            metric(offQ+q,offP+p) = metric(offP+p,offQ+q)
         enddo
      enddo
      if (ish .eq. jsh) then
         maxdiag = 0.0d0
         do p = 1,di
            maxdiag = max(maxdiag, metric(offP+p,offP+p))
         enddo
         aux_shell_bound(ish-nBases+1) = sqrt(max(maxdiag,0.0d0))
      endif
      deallocate(buf2c)
   enddo
enddo
!$omp end parallel do
deallocate(aux_offset)

maxdiag = maxval(aux_shell_bound)

block
   integer,allocatable :: fn_shell_est(:)
   integer :: m_est, n_est, ish_est, jsh_est, off_est, di_est
   allocate(fn_shell_est(nConts))
   off_est = 0
   do ish_est = 0,nBases-1
      di_est = cgto_engine(ish_est, basDF)
      fn_shell_est(off_est+1:off_est+di_est) = ish_est
      off_est = off_est + di_est
   enddo
   df_total_pairs_est = 0_8
   do m_est = 1,nConts
      ish_est = fn_shell_est(m_est)
      do n_est = 1,nConts
         jsh_est = fn_shell_est(n_est)
         if (schwarz_bound(ish_est,jsh_est)*maxdiag .ge. SCHWARZ_CUTOFF) &
            df_total_pairs_est = df_total_pairs_est + 1_8
      enddo
   enddo
   deallocate(fn_shell_est)
end block

block
   use GRID_info, only: grid_cache_capacity, ngrids, max_batch_nsig
   use mod_mem_predict, only: predict_energy_peak_bytes, GRID_RAGGED_FRACTION_EST
   integer(8) :: df_need_bytes, df_avail_bytes, grid_cache_bytes
   character(len=8) :: force_direct_env
   logical :: need_k_est
   need_k_est = (HF_exchange_frac .gt. 1.0d-12)
   grid_cache_bytes = int(int(max_batch_nsig,8)*int(grid_cache_capacity,8)*4_8*4_8 &
                          *GRID_RAGGED_FRACTION_EST, 8) &
                     + 8_8*int(ngrids,8)*8_8
   df_need_bytes = predict_energy_peak_bytes(nConts, nAtoms, nContsAux, need_k_est, &
                                              .true., grid_cache_bytes, &
                                              rs_omega_active=(RS_omega .gt. 0.0d0 .and. engine_use_df_k), &
                                              df_total_pairs_actual=df_total_pairs_est)
   df_avail_bytes = engine_avail_at_start_bytes
   df_direct_mode = .true.
   if (engine_use_df_k) df_direct_mode = (df_need_bytes > df_avail_bytes)
   if (.not. need_k_est .and. RS_omega .le. 0.0d0) &
      df_direct_mode = (df_need_bytes > df_avail_bytes)
   force_direct_env = ""
   call get_environment_variable("ENGINE_FORCE_DF_DIRECT", force_direct_env)
   if (trim(force_direct_env) .eq. "1") df_direct_mode = .true.
   force_direct_env = ""
   call get_environment_variable("ENGINE_FORCE_DF_STORE", force_direct_env)
   if (trim(force_direct_env) .eq. "1") df_direct_mode = .false.
   df_energy_need_gb = real(df_need_bytes,8)/1024**3
   df_energy_avail_gb = real(df_avail_bytes,8)/1024**3
   if (engine_verbose .ge. 2) &
      print '(A)', "  [round 1/energy] predicted total peak="//trim(fmt_gb(df_energy_need_gb))// &
              "GB  available="//trim(fmt_gb(df_energy_avail_gb))//"GB"
   print '(A)', "  Build mode: "//trim(merge("DIRECT (in-core)       ","STORE (density fitting)",df_direct_mode))
end block
if (df_direct_mode) then
   call build_df_triple_list()
endif
if (.not. df_direct_mode) then
block
   integer,allocatable :: fn_shell(:)
   integer :: m,n,ish3,jsh3,slot
   allocate(fn_shell(nConts))
   offi = 0
   do ish = 0,nBases-1
      di = cgto_engine(ish, basDF)
      do p = 1,di
         fn_shell(offi+p) = ish
      enddo
      offi = offi + di
   enddo
   allocate(df_row_start(nConts), df_row_count(nConts))
   slot = 0
   do m = 1,nConts
      ish3 = fn_shell(m)
      df_row_start(m) = slot+1
      do n = 1,nConts
         jsh3 = fn_shell(n)
         if (schwarz_bound(ish3,jsh3)*maxdiag .ge. SCHWARZ_CUTOFF) slot = slot+1
      enddo
      df_row_count(m) = slot - df_row_start(m) + 1
   enddo
   df_total_pairs = slot
   df_max_row_count = maxval(df_row_count)
   allocate(df_pair_row(df_total_pairs), df_pair_col(df_total_pairs))
   allocate(pair_slot_2d(nConts,nConts))
   pair_slot_2d = 0
   slot = 0
   do m = 1,nConts
      ish3 = fn_shell(m)
      do n = 1,nConts
         jsh3 = fn_shell(n)
         if (schwarz_bound(ish3,jsh3)*maxdiag .ge. SCHWARZ_CUTOFF) then
            slot = slot+1
            df_pair_row(slot) = m
            df_pair_col(slot) = n
            pair_slot_2d(m,n) = slot
         endif
      enddo
   enddo
   deallocate(fn_shell)
   allocate(df_pair_mirror(df_total_pairs))
   !$omp parallel do
   do slot = 1,df_total_pairs
      df_pair_mirror(slot) = pair_slot_2d(df_pair_col(slot), df_pair_row(slot))
   enddo
   !$omp end parallel do
end block

call prof_start("df_compact_fill")
allocate(dfB_compact(df_total_pairs,nContsAux))
dfB_compact = 0.0d0
allocate(aux_offset(0:nBases-1))
offi = 0
do ish = 0,nBases-1
   aux_offset(ish) = offi
   offi = offi + cgto_engine(ish, basDF)
enddo
!$omp parallel do private(ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,shls,buf3c,val3c,slot_pq,slot_qp) &
!$omp&            schedule(dynamic)
do ish = 0,nBases-1
   di = cgto_engine(ish, basDF)
   offi = aux_offset(ish)
   do jsh = 0,ish
      dj = cgto_engine(jsh, basDF)
      offj = aux_offset(jsh)
      if (schwarz_bound(ish,jsh)*maxdiag .lt. SCHWARZ_CUTOFF) then
         cycle
      endif
      offP = 0
      do ksh = nBases,nBases+nBasesAux-1
         dk = cgto_engine(ksh, basDF)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) then
            offP = offP + dk
            cycle
         endif
         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         allocate(buf3c(di,dj,dk))
         call threec2e_engine(buf3c, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val3c = buf3c(p,q,r) * NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  slot_pq = pair_slot_2d(offi+p, offj+q)
                  dfB_compact(slot_pq, offP+r) = val3c
                  if (offi+p .ne. offj+q) then
                     slot_qp = pair_slot_2d(offj+q, offi+p)
                     dfB_compact(slot_qp, offP+r) = val3c
                  endif
               enddo
            enddo
         enddo
         deallocate(buf3c)
         offP = offP + dk
      enddo
   enddo
enddo
!$omp end parallel do
deallocate(aux_offset)
deallocate(pair_slot_2d)
call prof_stop("df_compact_fill")
endif

allocate(eval(nContsAux))
allocate(work(1))
allocate(iwork(1))
call dsyevd('V', 'L', nContsAux, metric, nContsAux, eval, work, -1, iwork, -1, info)
lwork = int(work(1))
liwork = iwork(1)
deallocate(work,iwork)
allocate(work(lwork))
allocate(iwork(liwork))
call openblas_set_num_threads(omp_get_max_threads())
call dsyevd('V', 'L', nContsAux, metric, nContsAux, eval, work, lwork, iwork, liwork, info)
call openblas_set_num_threads(1)
deallocate(work,iwork)
if (info .ne. 0) then
   print *, "DF metric eigendecomposition failed, info=", info
   stop 1
endif
allocate(df_evec(nContsAux,nContsAux))
allocate(df_evalinv(nContsAux))
df_evec = metric
do i = 1,nContsAux
   if (eval(i) .gt. DF_EVAL_FLOOR*eval(nContsAux)) then
      df_evalinv(i) = 1.0d0/eval(i)
   else
      df_evalinv(i) = 0.0d0
   endif
enddo
deallocate(eval)
deallocate(metric)
end subroutine build_df_integrals

module subroutine build_df_integrals_lr()
use MOL_info
use mod_exchange, only: RS_omega
use mod_profile, only: prof_start, prof_stop, fmt_gb, itoa
implicit none
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP
integer :: shls(4)
real(8),allocatable :: buf3c(:,:,:)
real(8) :: maxdiag, val3c
integer,allocatable :: aux_offset(:)
integer,allocatable :: pair_slot_2d(:,:)
integer :: slot, slot_pq, slot_qp

if (.not. (RS_omega .gt. 0.0d0)) return

maxdiag = maxval(aux_shell_bound)

allocate(pair_slot_2d(nConts,nConts))
pair_slot_2d = 0
do slot = 1,df_total_pairs
   pair_slot_2d(df_pair_row(slot), df_pair_col(slot)) = slot
enddo

envDF(9) = RS_omega
call threec2e_optimizer_engine(int3c2e_opt_lr, atm, nAtoms, basDF, nBases+nBasesAux, envDF)

call prof_start("df_compact_fill_lr")
allocate(dfB_compact_LR(df_total_pairs,nContsAux))
dfB_compact_LR = 0.0d0
allocate(aux_offset(0:nBases-1))
offi = 0
do ish = 0,nBases-1
   aux_offset(ish) = offi
   offi = offi + cgto_engine(ish, basDF)
enddo
!$omp parallel do private(ish,jsh,ksh,p,q,r,di,dj,dk,offi,offj,offP,shls,buf3c,val3c,slot_pq,slot_qp) &
!$omp&            schedule(dynamic)
do ish = 0,nBases-1
   di = cgto_engine(ish, basDF)
   offi = aux_offset(ish)
   do jsh = 0,ish
      dj = cgto_engine(jsh, basDF)
      offj = aux_offset(jsh)
      if (schwarz_bound(ish,jsh)*maxdiag .lt. SCHWARZ_CUTOFF) then
         cycle
      endif
      offP = 0
      do ksh = nBases,nBases+nBasesAux-1
         dk = cgto_engine(ksh, basDF)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) then
            offP = offP + dk
            cycle
         endif
         shls(1) = ish
         shls(2) = jsh
         shls(3) = ksh
         allocate(buf3c(di,dj,dk))
         call threec2e_engine(buf3c, shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt_lr)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  val3c = buf3c(p,q,r) * NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
                  slot_pq = pair_slot_2d(offi+p, offj+q)
                  dfB_compact_LR(slot_pq, offP+r) = val3c
                  if (offi+p .ne. offj+q) then
                     slot_qp = pair_slot_2d(offj+q, offi+p)
                     dfB_compact_LR(slot_qp, offP+r) = val3c
                  endif
               enddo
            enddo
         enddo
         deallocate(buf3c)
         offP = offP + dk
      enddo
   enddo
enddo
!$omp end parallel do
deallocate(aux_offset)
deallocate(pair_slot_2d)
envDF(9) = 0.0d0
call prof_stop("df_compact_fill_lr")
end subroutine build_df_integrals_lr

module subroutine decide_df_force_mode()
use MOL_info, only: nConts, nAtoms, engine_verbose
use GRID_info, only: grid_cache_capacity, ngrids, max_batch_nsig
use mod_meminfo, only: engine_avail_at_start_bytes
use mod_exchange, only: HF_exchange_frac, RS_omega
use mod_mem_predict, only: predict_energy_peak_bytes, predict_force_extra_bytes, &
                            GRID_RAGGED_FRACTION_EST
use mod_profile, only: fmt_gb
implicit none
integer(8) :: total_need_bytes, df_avail_bytes, grid_cache_bytes
character(len=8) :: force_direct_env
logical :: need_k_est

if (df_direct_mode) then
   df_force_direct_mode = .true.
   return
endif

need_k_est = (HF_exchange_frac .gt. 1.0d-12)
grid_cache_bytes = int(int(max_batch_nsig,8)*int(grid_cache_capacity,8)*4_8*4_8 &
                       *GRID_RAGGED_FRACTION_EST, 8) &
                  + 8_8*int(ngrids,8)*8_8
total_need_bytes = predict_energy_peak_bytes(nConts, nAtoms, nContsAux, need_k_est, &
                                              .true., grid_cache_bytes, &
                                              rs_omega_active=(RS_omega .gt. 0.0d0 .and. engine_use_df_k), &
                                              df_total_pairs_actual=int(df_total_pairs,8)) &
                  + predict_force_extra_bytes(nConts, nContsAux, need_k_est, .true.)
df_avail_bytes = engine_avail_at_start_bytes
df_force_direct_mode = (total_need_bytes > df_avail_bytes)
df_force_need_gb = real(total_need_bytes,8)/1024**3
df_force_avail_gb = real(df_avail_bytes,8)/1024**3

force_direct_env = ""
call get_environment_variable("ENGINE_FORCE_DF_FORCE_DIRECT", force_direct_env)
if (trim(force_direct_env) .eq. "1") df_force_direct_mode = .true.

if (df_force_direct_mode) call build_df_triple_list()

if (engine_verbose .ge. 2) &
   print '(A)', "  [round 2/force] predicted total peak w/ Tstack="//trim(fmt_gb(df_force_need_gb))// &
           "GB  available="//trim(fmt_gb(df_force_avail_gb))// &
           "GB  force DF mode="//trim(merge("DIRECT","STORE ",df_force_direct_mode))
end subroutine decide_df_force_mode

module subroutine ensure_df_built()
use mod_exchange, only: RS_omega, cosx_enabled
implicit none
if (.not. df_built) then
   if (len_trim(engine_df_aux_basis) .gt. 0) then
      call build_df_aux_basis_fromfile(trim(engine_df_aux_basis))
   else
      call build_df_aux_basis()
   endif
   call build_df_integrals()
   if (RS_omega .gt. 0.0d0 .and. engine_use_df_k) then
      if (df_direct_mode) then
         call ensure_df_lr_optimizer_built()
      else
         call build_df_integrals_lr()
      endif
   endif
   df_built = .true.
endif
end subroutine ensure_df_built

module subroutine build_df_metric_invhalf()
use omp_lib, only: omp_get_max_threads
implicit none
real(8),allocatable :: Vscaled(:,:)
integer :: i

allocate(Vscaled(nContsAux,nContsAux))
do i = 1,nContsAux
   Vscaled(:,i) = df_evec(:,i)*sqrt(df_evalinv(i))
enddo
allocate(df_Minvhalf(nContsAux,nContsAux))
call openblas_set_num_threads(omp_get_max_threads())
call dgemm('N','T', nContsAux, nContsAux, nContsAux, 1.0d0, &
           Vscaled, nContsAux, df_evec, nContsAux, 0.0d0, df_Minvhalf, nContsAux)
call openblas_set_num_threads(1)
deallocate(Vscaled)
df_Minvhalf_built = .true.
end subroutine build_df_metric_invhalf

module subroutine build_df_minv_full()
use omp_lib, only: omp_get_max_threads
implicit none

if (.not. df_Minvhalf_built) call build_df_metric_invhalf()
allocate(df_Minv_full(nContsAux,nContsAux))
call openblas_set_num_threads(omp_get_max_threads())
call dgemm('N','N', nContsAux, nContsAux, nContsAux, 1.0d0, &
           df_Minvhalf, nContsAux, df_Minvhalf, nContsAux, 0.0d0, df_Minv_full, nContsAux)
call openblas_set_num_threads(1)
df_Minv_full_built = .true.
end subroutine build_df_minv_full

module subroutine ensure_df_lr_optimizer_built()
use MOL_info
use mod_exchange, only: RS_omega
implicit none

if (df_lr_optimizer_built) return
envDF(9) = RS_omega
call threec2e_optimizer_engine(int3c2e_opt_lr, atm, nAtoms, basDF, nBases+nBasesAux, envDF)
envDF(9) = 0.0d0
df_lr_optimizer_built = .true.
end subroutine ensure_df_lr_optimizer_built

end submodule df_setup_impl
