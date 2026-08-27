! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! Density-fitting (RI-J) Coulomb matrix build.

submodule (mod_integrals) df_coulomb_impl
implicit none
contains

module subroutine integrals_build_coulomb_df(nConts_in, Ptot, J)
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: J(nConts_in,nConts_in)
real(8),allocatable :: dP(:,:), dJ(:,:)
character(len=8) :: envval
logical :: do_full

if (.not. incr_fock_env_checked) then
   envval = ""
   call get_environment_variable("ENGINE_INCREMENTAL_FOCK", envval)
   if (trim(envval) .eq. "0") incr_fock_enabled = .false.
   incr_fock_env_checked = .true.
endif

do_full = (.not. incr_fock_enabled) .or. (.not. incr_has_prev) &
          .or. (incr_since_full .ge. incr_full_period)

if (do_full) then
   call df_build_coulomb_core(nConts_in, Ptot, J)
   incr_since_full = 0
else
   allocate(dP(nConts_in,nConts_in), dJ(nConts_in,nConts_in))
   dP = Ptot - incr_Ptot_prev
   call df_build_coulomb_core(nConts_in, dP, dJ, incremental_call=.true.)
   J = incr_J_prev + dJ
   deallocate(dP, dJ)
   incr_since_full = incr_since_full + 1
endif

if (incr_fock_enabled) then
   if (.not. allocated(incr_Ptot_prev)) then
      allocate(incr_Ptot_prev(nConts_in,nConts_in), incr_J_prev(nConts_in,nConts_in))
   endif
   incr_Ptot_prev = Ptot
   incr_J_prev = J
   incr_has_prev = .true.
endif
end subroutine integrals_build_coulomb_df

module subroutine df_build_coulomb_core(nConts_in, Ptot, J, incremental_call)
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: J(nConts_in,nConts_in)
logical,intent(in),optional :: incremental_call
real(8),allocatable :: bvec(:), tvec(:), cvec(:), bprime(:), Ptot_compact(:), J_compact(:)
integer :: slot

call ensure_df_built()

if (df_direct_mode) then
   call df_build_coulomb_direct(nConts_in, Ptot, J, incremental_call)
   return
endif

allocate(Ptot_compact(df_total_pairs), J_compact(df_total_pairs))
!$omp parallel do private(slot)
do slot = 1,df_total_pairs
   Ptot_compact(slot) = Ptot(df_pair_row(slot), df_pair_col(slot))
enddo
!$omp end parallel do
call openblas_set_num_threads(omp_get_max_threads())
if (dfK_built) then
   allocate(bprime(nContsAux))
   call dgemv('T', df_total_pairs, nContsAux, 1.0d0, dfK_compact, df_total_pairs, &
              Ptot_compact, 1, 0.0d0, bprime, 1)
   call dgemv('N', df_total_pairs, nContsAux, 1.0d0, dfK_compact, df_total_pairs, &
              bprime, 1, 0.0d0, J_compact, 1)
   deallocate(bprime)
else
   allocate(bvec(nContsAux), tvec(nContsAux), cvec(nContsAux))
   call dgemv('T', df_total_pairs, nContsAux, 1.0d0, dfB_compact, df_total_pairs, &
              Ptot_compact, 1, 0.0d0, bvec, 1)
   call dgemv('T', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, bvec, 1, 0.0d0, tvec, 1)
   tvec = tvec * df_evalinv
   call dgemv('N', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, tvec, 1, 0.0d0, cvec, 1)
   if (.not. allocated(df_direct_cvec_cache)) then
      allocate(df_direct_cvec_cache(nContsAux))
      allocate(df_direct_cvec_cache_Ptot(nConts_in,nConts_in))
   endif
   df_direct_cvec_cache = cvec
   df_direct_cvec_cache_Ptot = Ptot
   df_direct_cvec_cache_valid = .true.
   call dgemv('N', df_total_pairs, nContsAux, 1.0d0, dfB_compact, df_total_pairs, &
              cvec, 1, 0.0d0, J_compact, 1)
   deallocate(bvec,tvec,cvec)
endif
call openblas_set_num_threads(1)
J = 0.0d0
!$omp parallel do private(slot)
do slot = 1,df_total_pairs
   J(df_pair_row(slot), df_pair_col(slot)) = J_compact(slot)
enddo
!$omp end parallel do
deallocate(Ptot_compact,J_compact)
end subroutine df_build_coulomb_core

module subroutine build_df_triple_list()
implicit none
integer :: ish,jsh,ksh,di,dj,dk,offP,t
real(8) :: maxdiag
integer,allocatable :: count_per_ish(:), offset_per_ish(:)

maxdiag = maxval(aux_shell_bound)

if (allocated(df_shell_ncgto)) deallocate(df_shell_ncgto, df_shell_offset)
allocate(df_shell_ncgto(0:nBases+nBasesAux-1), df_shell_offset(0:nBases+nBasesAux-1))
t = 0
do ish = 0,nBases-1
   df_shell_ncgto(ish) = cgto_engine(ish, basDF)
   df_shell_offset(ish) = t
   t = t + df_shell_ncgto(ish)
enddo
t = 0
do ksh = nBases,nBases+nBasesAux-1
   df_shell_ncgto(ksh) = cgto_engine(ksh, basDF)
   df_shell_offset(ksh) = t
   t = t + df_shell_ncgto(ksh)
enddo

allocate(count_per_ish(0:nBases-1))
!$omp parallel do private(ish,jsh,ksh) schedule(dynamic)
do ish = 0,nBases-1
   count_per_ish(ish) = 0
   do jsh = 0,ish
      if (schwarz_bound(ish,jsh)*maxdiag .ge. SCHWARZ_CUTOFF) then
         do ksh = nBases,nBases+nBasesAux-1
            if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .ge. SCHWARZ_CUTOFF) then
               count_per_ish(ish) = count_per_ish(ish) + 1
            endif
         enddo
      endif
   enddo
enddo
!$omp end parallel do

allocate(offset_per_ish(0:nBases-1))
t = 0
do ish = 0,nBases-1
   offset_per_ish(ish) = t
   t = t + count_per_ish(ish)
enddo
df_n_triples = t

if (nBases+nBasesAux .gt. 32767) then
   print *, "ERROR: nBases+nBasesAux =",nBases+nBasesAux," exceeds int16 range (32767) - ", &
            "df_triple_ish/jsh/ksh cannot safely index shells this large. ", &
            "See df_triple_ish's own module-level comment (mod_integrals.f90)."
   stop 1
endif

if (allocated(df_triple_ish)) deallocate(df_triple_ish,df_triple_jsh,df_triple_ksh)
allocate(df_triple_ish(df_n_triples),df_triple_jsh(df_n_triples),df_triple_ksh(df_n_triples))

!$omp parallel do private(ish,jsh,ksh,t) schedule(dynamic)
do ish = 0,nBases-1
   t = offset_per_ish(ish)
   do jsh = 0,ish
      if (schwarz_bound(ish,jsh)*maxdiag .ge. SCHWARZ_CUTOFF) then
         do ksh = nBases,nBases+nBasesAux-1
            if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .ge. SCHWARZ_CUTOFF) then
               t = t + 1
               df_triple_ish(t) = ish; df_triple_jsh(t) = jsh; df_triple_ksh(t) = ksh
            endif
         enddo
      endif
   enddo
enddo
!$omp end parallel do
deallocate(count_per_ish, offset_per_ish)
end subroutine build_df_triple_list

module subroutine df_build_coulomb_direct(nConts_in, Ptot, J, incremental_call)
use MOL_info
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: J(nConts_in,nConts_in)
logical,intent(in),optional :: incremental_call
real(8),allocatable :: bvec(:), tvec(:), cvec(:)
real(8),allocatable :: buf3c_flat(:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offP,offi,offj,t
integer :: tid
integer :: shls(4)
real(8) :: val
real(8),allocatable :: dmax_shell(:,:), cvec_shell_bound(:)
real(8) :: maxdiag, this_cutoff
logical :: is_incremental
integer :: maxdi, maxdj, maxdk

is_incremental = present(incremental_call) .and. incremental_call
this_cutoff = merge(density_screen_cutoff_incr, density_screen_cutoff, is_incremental)
allocate(bvec(nContsAux), tvec(nContsAux), cvec(nContsAux))
bvec = 0.0d0
J = 0.0d0
maxdiag = maxval(aux_shell_bound)

if (.not. df_direct_pool_built) then
   allocate(df_direct_Jlocal_pool(nConts_in,nConts_in,omp_get_max_threads()))
   allocate(df_direct_bveclocal_pool(nContsAux,omp_get_max_threads()))
   df_direct_pool_built = .true.
endif

allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts_in, Ptot, dmax_shell)

maxdi = maxval(df_shell_ncgto(0:nBases-1))
maxdj = maxdi
maxdk = maxval(df_shell_ncgto(nBases:nBases+nBasesAux-1))

call prof_start("df_direct_bpass")
!$omp parallel private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf3c_flat,p,q,r,val,t) &
!$omp&   private(tid)
tid = omp_get_thread_num() + 1
df_direct_bveclocal_pool(:,tid) = 0.0d0
allocate(buf3c_flat(maxdi*maxdj*maxdk))
!$omp do schedule(static,1)
do t = 1,df_n_triples
   ish = df_triple_ish(t); jsh = df_triple_jsh(t); ksh = df_triple_ksh(t)
   di = df_shell_ncgto(ish); dj = df_shell_ncgto(jsh); dk = df_shell_ncgto(ksh)
   offi = df_shell_offset(ish); offj = df_shell_offset(jsh); offP = df_shell_offset(ksh)
   if (schwarz_bound(ish,jsh)*maxdiag*dmax_shell(ish,jsh) .lt. this_cutoff) cycle
   shls(1) = ish
   shls(2) = jsh
   shls(3) = ksh
   call threec2e_engine(buf3c_flat(1:di*dj*dk), shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
   do r = 1,dk
      do q = 1,dj
         do p = 1,di
            val = buf3c_flat((r-1)*di*dj+(q-1)*di+p)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            if (ish .ne. jsh) then
               df_direct_bveclocal_pool(offP+r,tid) = df_direct_bveclocal_pool(offP+r,tid) &
                                                      + val*(Ptot(offi+p,offj+q)+Ptot(offj+q,offi+p))
            else
               df_direct_bveclocal_pool(offP+r,tid) = df_direct_bveclocal_pool(offP+r,tid) &
                                                      + val*Ptot(offi+p,offj+q)
            endif
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf3c_flat)
block
   integer :: merge_tid, merge_nthreads
   merge_nthreads = omp_get_num_threads()
   do merge_tid = 0, merge_nthreads-1
      !$omp barrier
      if (omp_get_thread_num() .eq. merge_tid) then
         bvec = bvec + df_direct_bveclocal_pool(:,tid)
      endif
   enddo
end block
!$omp end parallel
call prof_stop("df_direct_bpass")

call prof_start("df_direct_solve")
call dgemv('T', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, bvec, 1, 0.0d0, tvec, 1)
tvec = tvec * df_evalinv
call dgemv('N', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, tvec, 1, 0.0d0, cvec, 1)
if (.not. allocated(df_direct_cvec_cache)) then
   allocate(df_direct_cvec_cache(nContsAux))
   allocate(df_direct_cvec_cache_Ptot(nConts_in,nConts_in))
endif
df_direct_cvec_cache = cvec
df_direct_cvec_cache_Ptot = Ptot
df_direct_cvec_cache_valid = .true.
allocate(cvec_shell_bound(nBasesAux))
offP = 0
do ksh = nBases,nBases+nBasesAux-1
   dk = cgto_engine(ksh, basDF)
   cvec_shell_bound(ksh-nBases+1) = maxval(abs(cvec(offP+1:offP+dk)))
   offP = offP + dk
enddo
call prof_stop("df_direct_solve")

call prof_start("df_direct_jpass")
!$omp parallel private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf3c_flat,p,q,r,val,t) &
!$omp&   private(tid)
tid = omp_get_thread_num() + 1
df_direct_Jlocal_pool(:,:,tid) = 0.0d0
allocate(buf3c_flat(maxdi*maxdj*maxdk))
!$omp do schedule(static,1)
do t = 1,df_n_triples
   ish = df_triple_ish(t); jsh = df_triple_jsh(t); ksh = df_triple_ksh(t)
   di = df_shell_ncgto(ish); dj = df_shell_ncgto(jsh); dk = df_shell_ncgto(ksh)
   offi = df_shell_offset(ish); offj = df_shell_offset(jsh); offP = df_shell_offset(ksh)
   if (schwarz_bound(ish,jsh)*cvec_shell_bound(ksh-nBases+1) .lt. this_cutoff) cycle
   shls(1) = ish
   shls(2) = jsh
   shls(3) = ksh
   call threec2e_engine(buf3c_flat(1:di*dj*dk), shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
   do r = 1,dk
      do q = 1,dj
         do p = 1,di
            val = buf3c_flat((r-1)*di*dj+(q-1)*di+p)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)*cvec(offP+r)
            df_direct_Jlocal_pool(offi+p,offj+q,tid) = df_direct_Jlocal_pool(offi+p,offj+q,tid) + val
            if (ish .ne. jsh) then
               df_direct_Jlocal_pool(offj+q,offi+p,tid) = df_direct_Jlocal_pool(offj+q,offi+p,tid) + val
            endif
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf3c_flat)
block
   integer :: merge_tid, merge_nthreads
   merge_nthreads = omp_get_num_threads()
   do merge_tid = 0, merge_nthreads-1
      !$omp barrier
      if (omp_get_thread_num() .eq. merge_tid) then
         J = J + df_direct_Jlocal_pool(:,:,tid)
      endif
   enddo
end block
!$omp end parallel
call prof_stop("df_direct_jpass")

deallocate(bvec,tvec,cvec)
deallocate(dmax_shell,cvec_shell_bound)
end subroutine df_build_coulomb_direct

end submodule df_coulomb_impl
