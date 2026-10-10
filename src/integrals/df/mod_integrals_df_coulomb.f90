! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

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

module subroutine integrals_build_coulomb_df_batch(nConts, nbatch, Ptot, J)
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts, nbatch
real(8),intent(in) :: Ptot(nConts,nConts,nbatch)
real(8),intent(out) :: J(nConts,nConts,nbatch)
integer :: b, slot
real(8),allocatable :: Ptot_compact(:,:), J_compact(:,:)
real(8),allocatable :: Bvec(:,:), Tvec(:,:), Cvec(:,:)
external :: dgemm

call ensure_df_built()

if (df_direct_mode) then
   block
      character(len=8) :: e
      call get_environment_variable("ENGINE_DFJ_BATCH_LOOP", e)
      if (trim(e) .eq. "1") then
         do b = 1,nbatch
            call df_build_coulomb_direct(nConts, Ptot(:,:,b), J(:,:,b))
         enddo
      else
         call df_build_coulomb_direct_batch(nConts, nbatch, Ptot, J)
      endif
   end block
   return
endif

allocate(Ptot_compact(df_total_pairs,nbatch), J_compact(df_total_pairs,nbatch))
!$omp parallel do private(slot,b)
do slot = 1,df_total_pairs
   do b = 1,nbatch
      Ptot_compact(slot,b) = Ptot(df_pair_row(slot), df_pair_col(slot), b)
   enddo
enddo
!$omp end parallel do
call openblas_set_num_threads(omp_get_max_threads())
if (dfK_built) then
   allocate(Bvec(nContsAux,nbatch))
   call dgemm('T','N', nContsAux, nbatch, df_total_pairs, 1.0d0, dfK_compact, df_total_pairs, &
              Ptot_compact, df_total_pairs, 0.0d0, Bvec, nContsAux)
   call dgemm('N','N', df_total_pairs, nbatch, nContsAux, 1.0d0, dfK_compact, df_total_pairs, &
              Bvec, nContsAux, 0.0d0, J_compact, df_total_pairs)
   deallocate(Bvec)
else
   allocate(Bvec(nContsAux,nbatch), Tvec(nContsAux,nbatch), Cvec(nContsAux,nbatch))
   call dgemm('T','N', nContsAux, nbatch, df_total_pairs, 1.0d0, dfB_compact, df_total_pairs, &
              Ptot_compact, df_total_pairs, 0.0d0, Bvec, nContsAux)
   call dgemm('T','N', nContsAux, nbatch, nContsAux, 1.0d0, df_evec, nContsAux, &
              Bvec, nContsAux, 0.0d0, Tvec, nContsAux)
   do b = 1,nbatch
      do slot = 1,nContsAux
         Tvec(slot,b) = Tvec(slot,b)*df_evalinv(slot)
      enddo
   enddo
   call dgemm('N','N', nContsAux, nbatch, nContsAux, 1.0d0, df_evec, nContsAux, &
              Tvec, nContsAux, 0.0d0, Cvec, nContsAux)
   if (.not. allocated(df_direct_cvec_cache)) then
      allocate(df_direct_cvec_cache(nContsAux))
      allocate(df_direct_cvec_cache_Ptot(nConts,nConts))
   endif
   df_direct_cvec_cache = Cvec(:,nbatch)
   df_direct_cvec_cache_Ptot = Ptot(:,:,nbatch)
   df_direct_cvec_cache_valid = .true.
   call dgemm('N','N', df_total_pairs, nbatch, nContsAux, 1.0d0, dfB_compact, df_total_pairs, &
              Cvec, nContsAux, 0.0d0, J_compact, df_total_pairs)
   deallocate(Bvec,Tvec,Cvec)
endif
call openblas_set_num_threads(1)
J = 0.0d0
!$omp parallel do private(slot,b)
do slot = 1,df_total_pairs
   do b = 1,nbatch
      J(df_pair_row(slot), df_pair_col(slot), b) = J_compact(slot,b)
   enddo
enddo
!$omp end parallel do
deallocate(Ptot_compact, J_compact)
end subroutine integrals_build_coulomb_df_batch

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

module subroutine integrals_build_coulomb_df_batch_c(nConts, nb, Lf, nop, Co, JC)
use omp_lib, only: omp_get_max_threads
implicit none
integer,intent(in) :: nConts, nb, nop
real(8),intent(in) :: Lf(nConts,nop,nb), Co(nConts,nop)
real(8),intent(out) :: JC(nConts,nop,nb)
real(8),allocatable :: Pc(:,:), Jc2(:,:), Bvec(:,:), Tvec(:,:), Cvec(:,:), Cg(:,:), Lm(:,:), Tm(:,:)
real(8),allocatable :: Dch(:,:,:), Jch(:,:,:)
integer :: m, s0, cnt, b, k, sl, sm, o, c0, nch
real(8) :: t
external :: dgemm

call ensure_df_built()
if (df_direct_mode .and. .not. dfK_built) then
   block
      character(len=8) :: e
      call get_environment_variable("ENGINE_DFJ_BATCHC_DENSE", e)
      if (trim(e) .ne. "1") then
         call df_build_coulomb_direct_batch_c(nConts, nb, Lf, nop, Co, JC)
         return
      endif
   end block
endif
if (df_direct_mode .or. dfK_built) then
   allocate(Dch(nConts,nConts,16), Jch(nConts,nConts,16))
   do c0 = 0, nb-1, 16
      nch = min(16, nb-c0)
      do b = 1,nch
         call dgemm('N','T',nConts,nConts,nop,1.0d0,Lf(1,1,c0+b),nConts,Co,nConts,0.0d0,Dch(1,1,b),nConts)
         Dch(:,:,b) = Dch(:,:,b) + transpose(Dch(:,:,b))
      enddo
      call integrals_build_coulomb_df_batch(nConts, nch, Dch(:,:,1:nch), Jch(:,:,1:nch))
      do b = 1,nch
         call dgemm('N','N',nConts,nop,nConts,1.0d0,Jch(1,1,b),nConts,Co,nConts,0.0d0,JC(1,1,c0+b),nConts)
      enddo
   enddo
   deallocate(Dch, Jch)
   return
endif

allocate(Pc(df_total_pairs,nb), Jc2(df_total_pairs,nb))
!$omp parallel private(m,s0,cnt,k,Cg,Lm,b)
allocate(Cg(df_max_row_count,nop), Lm(nop,nb))
!$omp do schedule(dynamic,8)
do m = 1,nConts
   s0 = df_row_start(m) ; cnt = df_row_count(m)
   if (cnt .eq. 0) cycle
   do k = 1,cnt
      Cg(k,:) = Co(df_pair_col(s0+k-1),:)
   enddo
   do b = 1,nb
      Lm(:,b) = Lf(m,:,b)
   enddo
   call dgemm('N','N',cnt,nb,nop,1.0d0,Cg,df_max_row_count,Lm,nop,0.0d0,Pc(s0,1),df_total_pairs)
enddo
!$omp end do
deallocate(Cg, Lm)
!$omp end parallel
!$omp parallel do private(sl,sm,t,b) schedule(static)
do sl = 1,df_total_pairs
   sm = df_pair_mirror(sl)
   if (sm .lt. sl) cycle
   do b = 1,nb
      t = Pc(sl,b) + Pc(sm,b)
      Pc(sl,b) = t
      Pc(sm,b) = t
   enddo
enddo
!$omp end parallel do
call openblas_set_num_threads(omp_get_max_threads())
allocate(Bvec(nContsAux,nb), Tvec(nContsAux,nb), Cvec(nContsAux,nb))
call dgemm('T','N', nContsAux, nb, df_total_pairs, 1.0d0, dfB_compact, df_total_pairs, &
           Pc, df_total_pairs, 0.0d0, Bvec, nContsAux)
call dgemm('T','N', nContsAux, nb, nContsAux, 1.0d0, df_evec, nContsAux, &
           Bvec, nContsAux, 0.0d0, Tvec, nContsAux)
do b = 1,nb
   Tvec(:,b) = Tvec(:,b)*df_evalinv
enddo
call dgemm('N','N', nContsAux, nb, nContsAux, 1.0d0, df_evec, nContsAux, &
           Tvec, nContsAux, 0.0d0, Cvec, nContsAux)
call dgemm('N','N', df_total_pairs, nb, nContsAux, 1.0d0, dfB_compact, df_total_pairs, &
           Cvec, nContsAux, 0.0d0, Jc2, df_total_pairs)
call openblas_set_num_threads(1)
deallocate(Bvec, Tvec, Cvec, Pc)
!$omp parallel private(m,s0,cnt,k,Cg,Tm,o)
allocate(Cg(df_max_row_count,nop), Tm(nb,nop))
!$omp do schedule(dynamic,8)
do m = 1,nConts
   s0 = df_row_start(m) ; cnt = df_row_count(m)
   if (cnt .eq. 0) then
      JC(m,:,:) = 0.0d0
      cycle
   endif
   do k = 1,cnt
      Cg(k,:) = Co(df_pair_col(s0+k-1),:)
   enddo
   call dgemm('T','N',nb,nop,cnt,1.0d0,Jc2(s0,1),df_total_pairs,Cg,df_max_row_count,0.0d0,Tm,nb)
   do o = 1,nop
      JC(m,o,:) = Tm(:,o)
   enddo
enddo
!$omp end do
deallocate(Cg, Tm)
!$omp end parallel
deallocate(Jc2)
end subroutine integrals_build_coulomb_df_batch_c

module subroutine df_build_coulomb_direct_batch(nConts_in, nb, Ptot, J)
use MOL_info
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads, omp_get_wtime
implicit none
integer,intent(in) :: nConts_in, nb
real(8),intent(in) :: Ptot(nConts_in,nConts_in,nb)
real(8),intent(out) :: J(nConts_in,nConts_in,nb)
real(8),allocatable :: Bm(:,:), Tm(:,:), Cm(:,:), Bloc(:,:), buf(:), bn(:,:), Dt(:,:), Jt(:,:)
real(8),allocatable :: dmax_shell(:,:), cbound(:), Pmax(:,:)
real(8),allocatable :: bnP(:,:), Tb(:,:), Cg(:,:)
integer,allocatable :: poff(:), pdk(:)
integer :: np, ncol, it0, it1, s0
integer,allocatable :: iptr(:), ilist(:), cnt(:), aosh(:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offP,offi,offj,t,it,b,maxdi,maxdk,mu,nu,pj,pdj,poffj
integer :: shls(4)
real(8) :: maxdiag, ti0, tint, tall0
integer(8) :: ncalls
character(len=8) :: venv
external :: dgemm

maxdi = maxval(df_shell_ncgto(0:nBases-1))
maxdk = maxval(df_shell_ncgto(nBases:nBases+nBasesAux-1))
maxdiag = maxval(aux_shell_bound)
allocate(Pmax(nConts_in,nConts_in), dmax_shell(0:nBases-1,0:nBases-1))
Pmax = 0.0d0
do b = 1,nb
   Pmax = max(Pmax, abs(Ptot(:,:,b)))
enddo
call compute_shell_density_bound(nConts_in, Pmax, dmax_shell)
deallocate(Pmax)

allocate(cnt(0:nBases-1), iptr(0:nBases), ilist(df_n_triples))
cnt = 0
do t = 1,df_n_triples
   cnt(df_triple_ish(t)) = cnt(df_triple_ish(t)) + 1
enddo
iptr(0) = 0
do ish = 0,nBases-1
   iptr(ish+1) = iptr(ish) + cnt(ish)
enddo
cnt = 0
do t = 1,df_n_triples
   ish = df_triple_ish(t)
   cnt(ish) = cnt(ish) + 1
   ilist(iptr(ish)+cnt(ish)) = t
enddo

call prof_start("df_direct_bpass")
tint = 0.0d0 ; ncalls = 0 ; tall0 = omp_get_wtime()
allocate(Bm(nContsAux,nb))
Bm = 0.0d0
!$omp parallel private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf,bn,Dt,Bloc,p,q,r,t,it,b,pj,ti0) &
!$omp&  private(bnP,Tb,poff,pdk,np,ncol,it0,it1,s0) reduction(+:tint,ncalls)
allocate(buf(maxdi*maxdi*maxdk), bn(maxdi*maxdi,maxdk), Dt(maxdi*maxdi,nb), Bloc(nContsAux,nb))
allocate(bnP(maxdi*maxdi,nContsAux), Tb(nContsAux,nb), poff(nBasesAux), pdk(nBasesAux))
Bloc = 0.0d0
!$omp do schedule(static,1)
do ish = 0,nBases-1
   di = df_shell_ncgto(ish); offi = df_shell_offset(ish)
   it0 = iptr(ish)+1
   do while (it0 .le. iptr(ish+1))
      jsh = df_triple_jsh(ilist(it0))
      it1 = it0
      do while (it1+1 .le. iptr(ish+1))
         if (df_triple_jsh(ilist(it1+1)) .ne. jsh) exit
         it1 = it1 + 1
      enddo
      if (schwarz_bound(ish,jsh)*maxdiag*dmax_shell(ish,jsh) .lt. density_screen_cutoff) then
         it0 = it1 + 1
         cycle
      endif
      dj = df_shell_ncgto(jsh); offj = df_shell_offset(jsh)
      do b = 1,nb
         do q = 1,dj
            do p = 1,di
               if (ish .ne. jsh) then
                  Dt((q-1)*di+p,b) = Ptot(offi+p,offj+q,b) + Ptot(offj+q,offi+p,b)
               else
                  Dt((q-1)*di+p,b) = Ptot(offi+p,offj+q,b)
               endif
            enddo
         enddo
      enddo
      np = 0 ; ncol = 0
      do it = it0,it1
         ksh = df_triple_ksh(ilist(it))
         dk = df_shell_ncgto(ksh); offP = df_shell_offset(ksh)
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         ti0 = omp_get_wtime()
         call threec2e_engine(buf(1:di*dj*dk), shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
         tint = tint + (omp_get_wtime() - ti0)
         ncalls = ncalls + 1
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  bnP((q-1)*di+p,ncol+r) = buf((r-1)*di*dj+(q-1)*di+p)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
               enddo
            enddo
         enddo
         np = np + 1
         poff(np) = offP ; pdk(np) = dk
         ncol = ncol + dk
      enddo
      call dgemm('T','N',ncol,nb,di*dj,1.0d0,bnP,maxdi*maxdi,Dt,maxdi*maxdi,0.0d0,Tb,nContsAux)
      do b = 1,nb
         s0 = 0
         do r = 1,np
            Bloc(poff(r)+1:poff(r)+pdk(r),b) = Bloc(poff(r)+1:poff(r)+pdk(r),b) + Tb(s0+1:s0+pdk(r),b)
            s0 = s0 + pdk(r)
         enddo
      enddo
      it0 = it1 + 1
   enddo
enddo
!$omp end do
block
   integer :: mt
   do mt = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) Bm = Bm + Bloc
   enddo
end block
deallocate(buf, bn, Dt, Bloc, bnP, Tb, poff, pdk)
!$omp end parallel
call prof_stop("df_direct_bpass")
call get_environment_variable("ENGINE_DFJ_BATCH_VERBOSE", venv)
if (trim(venv) .eq. "1") then
   print '(A,I0,A,I0,A,I0,A,F8.3,A,F9.3)', " DFJ batch nb=", nb, "  triples ", df_n_triples, " computed ", ncalls, &
         "  b-pass wall ", omp_get_wtime()-tall0, "  int thread-s ", tint
   call flush(6)
endif

call prof_start("df_direct_solve")
allocate(Tm(nContsAux,nb), Cm(nContsAux,nb))
call dgemm('T','N',nContsAux,nb,nContsAux,1.0d0,df_evec,nContsAux,Bm,nContsAux,0.0d0,Tm,nContsAux)
do b = 1,nb
   Tm(:,b) = Tm(:,b)*df_evalinv
enddo
call dgemm('N','N',nContsAux,nb,nContsAux,1.0d0,df_evec,nContsAux,Tm,nContsAux,0.0d0,Cm,nContsAux)
allocate(cbound(nBasesAux))
offP = 0
do ksh = nBases,nBases+nBasesAux-1
   dk = cgto_engine(ksh, basDF)
   cbound(ksh-nBases+1) = maxval(abs(Cm(offP+1:offP+dk,:)))
   offP = offP + dk
enddo
call prof_stop("df_direct_solve")

call prof_start("df_direct_jpass")
J = 0.0d0
!$omp parallel private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf,bn,Jt,p,q,r,t,it,b,pj,pdj,poffj) &
!$omp&  private(bnP,Cg,np,ncol,it0,it1)
allocate(buf(maxdi*maxdi*maxdk), bn(maxdi*maxdi,maxdk), Jt(maxdi*maxdi,nb))
allocate(bnP(maxdi*maxdi,nContsAux), Cg(nContsAux,nb))
!$omp do schedule(dynamic)
do ish = 0,nBases-1
   di = df_shell_ncgto(ish); offi = df_shell_offset(ish)
   it0 = iptr(ish)+1
   do while (it0 .le. iptr(ish+1))
      jsh = df_triple_jsh(ilist(it0))
      it1 = it0
      do while (it1+1 .le. iptr(ish+1))
         if (df_triple_jsh(ilist(it1+1)) .ne. jsh) exit
         it1 = it1 + 1
      enddo
      dj = df_shell_ncgto(jsh); offj = df_shell_offset(jsh)
      ncol = 0
      do it = it0,it1
         ksh = df_triple_ksh(ilist(it))
         if (schwarz_bound(ish,jsh)*cbound(ksh-nBases+1) .lt. density_screen_cutoff) cycle
         dk = df_shell_ncgto(ksh); offP = df_shell_offset(ksh)
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         call threec2e_engine(buf(1:di*dj*dk), shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  bnP((q-1)*di+p,ncol+r) = buf((r-1)*di*dj+(q-1)*di+p)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
               enddo
            enddo
         enddo
         do b = 1,nb
            Cg(ncol+1:ncol+dk,b) = Cm(offP+1:offP+dk,b)
         enddo
         ncol = ncol + dk
      enddo
      if (ncol .gt. 0) then
         call dgemm('N','N',di*dj,nb,ncol,1.0d0,bnP,maxdi*maxdi,Cg,nContsAux,0.0d0,Jt,maxdi*maxdi)
         do b = 1,nb
            do q = 1,dj
               do p = 1,di
                  J(offi+p,offj+q,b) = J(offi+p,offj+q,b) + Jt((q-1)*di+p,b)
               enddo
            enddo
         enddo
      endif
      it0 = it1 + 1
   enddo
enddo
!$omp end do
deallocate(buf, bn, Jt, bnP, Cg)
!$omp end parallel
allocate(aosh(nConts_in))
do ish = 0,nBases-1
   do p = 1,df_shell_ncgto(ish)
      aosh(df_shell_offset(ish)+p) = ish
   enddo
enddo
!$omp parallel do private(b,mu,nu,ish,jsh)
do b = 1,nb
   do nu = 1,nConts_in
      do mu = 1,nConts_in
         if (aosh(mu) .lt. aosh(nu)) J(mu,nu,b) = J(nu,mu,b)
      enddo
   enddo
enddo
!$omp end parallel do
call prof_stop("df_direct_jpass")
deallocate(Bm, Tm, Cm, cbound, dmax_shell, cnt, iptr, ilist, aosh)
end subroutine df_build_coulomb_direct_batch

module subroutine df_build_coulomb_direct_batch_c(nConts_in, nb, Lf, nop, Co, JC)
use MOL_info
use mod_profile, only: prof_start, prof_stop
use mod_fixq, only: to_q, from_q
use omp_lib, only: omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, omp_destroy_lock, &
                   omp_get_num_threads, omp_get_thread_num, omp_get_wtime
implicit none
integer,intent(in) :: nConts_in, nb, nop
real(8),intent(in) :: Lf(nConts_in,nop,nb), Co(nConts_in,nop)
real(8),intent(out) :: JC(nConts_in,nop,nb)
real(8),allocatable :: LB(:,:,:), CoT(:,:), Lnrm(:), Cnrm(:), smax(:,:)
integer(16),allocatable :: JB(:,:,:)
real(8),allocatable :: Bm(:,:), Tm(:,:), Cm(:,:), Bloc(:,:), buf(:), Dt(:,:), Tb(:,:), bnP(:,:)
real(8),allocatable :: Cg(:,:), Jt(:,:), Wi(:,:,:), Wj(:,:,:), cbound(:)
integer,allocatable :: iptr(:), ilist(:), cnt(:), poff(:), pdk(:)
integer(omp_lock_kind),allocatable :: slk(:)
integer :: ish,jsh,ksh,p,q,r,di,dj,dk,offP,offi,offj,t,it,b,maxdi,maxdk,m,o
integer :: np, ncol, it0, it1, s0, ii, r0
integer,allocatable :: ord(:)
real(8),allocatable :: wrk(:)
integer :: shls(4)
real(8) :: maxdiag, fd, sd, v, tq(10), tq0, tqa, tqb
character(len=8) :: venv
external :: dgemm

tq = 0.0d0 ; tq0 = omp_get_wtime()
maxdi = maxval(df_shell_ncgto(0:nBases-1))
maxdk = maxval(df_shell_ncgto(nBases:nBases+nBasesAux-1))
maxdiag = maxval(aux_shell_bound)
allocate(LB(nop,nb,nConts_in), CoT(nop,nConts_in), JB(nop,nb,nConts_in))
!$omp parallel do private(m,b)
do m = 1,nConts_in
   CoT(:,m) = Co(m,:)
   do b = 1,nb
      LB(:,b,m) = Lf(m,:,b)
   enddo
   JB(:,:,m) = 0_16
enddo
!$omp end parallel do
allocate(Lnrm(nConts_in), Cnrm(nConts_in), smax(0:nBases-1,0:nBases-1))
do m = 1,nConts_in
   Cnrm(m) = sqrt(sum(CoT(:,m)**2))
   Lnrm(m) = 0.0d0
   do b = 1,nb
      Lnrm(m) = max(Lnrm(m), sqrt(sum(LB(:,b,m)**2)))
   enddo
enddo
do jsh = 0,nBases-1
   do ish = 0,nBases-1
      sd = 0.0d0
      do q = df_shell_offset(jsh)+1, df_shell_offset(jsh)+df_shell_ncgto(jsh)
         do p = df_shell_offset(ish)+1, df_shell_offset(ish)+df_shell_ncgto(ish)
            sd = max(sd, Lnrm(p)*Cnrm(q) + Cnrm(p)*Lnrm(q))
         enddo
      enddo
      smax(ish,jsh) = sd
   enddo
enddo

allocate(cnt(0:nBases-1), iptr(0:nBases), ilist(df_n_triples))
cnt = 0
do t = 1,df_n_triples
   cnt(df_triple_ish(t)) = cnt(df_triple_ish(t)) + 1
enddo
iptr(0) = 0
do ish = 0,nBases-1
   iptr(ish+1) = iptr(ish) + cnt(ish)
enddo
cnt = 0
do t = 1,df_n_triples
   ish = df_triple_ish(t)
   cnt(ish) = cnt(ish) + 1
   ilist(iptr(ish)+cnt(ish)) = t
enddo

tq(1) = omp_get_wtime()-tq0 ; tq0 = omp_get_wtime()
allocate(ord(0:nBases-1), wrk(0:nBases-1))
do ish = 0,nBases-1
   ord(ish) = ish
   wrk(ish) = dble(iptr(ish+1)-iptr(ish))*dble(df_shell_ncgto(ish))
enddo
call sort_desc(nBases, wrk, ord)
call prof_start("df_direct_bpass")
allocate(Bm(nContsAux,nb))
Bm = 0.0d0
!$omp parallel private(ii,ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf,Dt,Bloc,p,q,r,it,b,fd,sd,tqa,r0) &
!$omp&  private(bnP,Tb,poff,pdk,np,ncol,it0,it1,s0)
allocate(buf(maxdi*maxdi*maxdk), Dt(maxdi*maxdi,nb), Bloc(nContsAux,nb))
allocate(bnP(maxdi*maxdi,nContsAux), Tb(nContsAux,nb), poff(nBasesAux), pdk(nBasesAux))
Bloc = 0.0d0
!$omp do schedule(static,1)
do ii = 0,nBases-1
   ish = ord(ii)
   di = df_shell_ncgto(ish); offi = df_shell_offset(ish)
   it0 = iptr(ish)+1
   do while (it0 .le. iptr(ish+1))
      jsh = df_triple_jsh(ilist(it0))
      it1 = it0
      do while (it1+1 .le. iptr(ish+1))
         if (df_triple_jsh(ilist(it1+1)) .ne. jsh) exit
         it1 = it1 + 1
      enddo
      if (schwarz_bound(ish,jsh)*maxdiag*smax(ish,jsh) .lt. density_screen_cutoff) then
         it0 = it1 + 1
         cycle
      endif
      dj = df_shell_ncgto(jsh); offj = df_shell_offset(jsh)
      fd = merge(2.0d0, 1.0d0, ish .ne. jsh)
      tqa = omp_get_wtime()
      do b = 1,nb
         do q = 1,dj
            do p = 1,di
               sd = 0.0d0
               do o = 1,nop
                  sd = sd + LB(o,b,offi+p)*CoT(o,offj+q) + CoT(o,offi+p)*LB(o,b,offj+q)
               enddo
               Dt((q-1)*di+p,b) = fd*sd
            enddo
         enddo
      enddo
      !$omp atomic
      tq(5) = tq(5) + (omp_get_wtime()-tqa)
      tqa = omp_get_wtime()
      np = 0 ; ncol = 0
      do it = it0,it1
         ksh = df_triple_ksh(ilist(it))
         dk = df_shell_ncgto(ksh); offP = df_shell_offset(ksh)
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         call threec2e_engine(buf(1:di*dj*dk), shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  bnP((q-1)*di+p,ncol+r) = buf((r-1)*di*dj+(q-1)*di+p)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
               enddo
            enddo
         enddo
         np = np + 1
         poff(np) = offP ; pdk(np) = dk
         ncol = ncol + dk
      enddo
      s0 = 0 ; r0 = 1
      do r = 1,np
         if (r .lt. np) then
            if (poff(r+1) .eq. poff(r)+pdk(r)) cycle
         endif
         call dgemm('T','N',poff(r)+pdk(r)-poff(r0),nb,di*dj,1.0d0,bnP(1,s0+1),maxdi*maxdi, &
                    Dt,maxdi*maxdi,1.0d0,Bloc(poff(r0)+1,1),nContsAux)
         s0 = s0 + poff(r)+pdk(r)-poff(r0)
         r0 = r + 1
      enddo
      !$omp atomic
      tq(6) = tq(6) + (omp_get_wtime()-tqa)
      it0 = it1 + 1
   enddo
enddo
!$omp end do
block
   integer :: mt
   do mt = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() .eq. mt) Bm = Bm + Bloc
   enddo
end block
deallocate(buf, Dt, Bloc, bnP, Tb, poff, pdk)
!$omp end parallel
call prof_stop("df_direct_bpass")
tq(2) = omp_get_wtime()-tq0 ; tq0 = omp_get_wtime()

call prof_start("df_direct_solve")
allocate(Tm(nContsAux,nb), Cm(nContsAux,nb))
call dgemm('T','N',nContsAux,nb,nContsAux,1.0d0,df_evec,nContsAux,Bm,nContsAux,0.0d0,Tm,nContsAux)
do b = 1,nb
   Tm(:,b) = Tm(:,b)*df_evalinv
enddo
call dgemm('N','N',nContsAux,nb,nContsAux,1.0d0,df_evec,nContsAux,Tm,nContsAux,0.0d0,Cm,nContsAux)
allocate(cbound(nBasesAux))
offP = 0
do ksh = nBases,nBases+nBasesAux-1
   dk = cgto_engine(ksh, basDF)
   cbound(ksh-nBases+1) = maxval(abs(Cm(offP+1:offP+dk,:)))
   offP = offP + dk
enddo
call prof_stop("df_direct_solve")

tq(3) = omp_get_wtime()-tq0 ; tq0 = omp_get_wtime()
call prof_start("df_direct_jpass")
allocate(slk(0:nBases-1))
do ish = 0,nBases-1
   call omp_init_lock(slk(ish))
enddo
!$omp parallel private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf,Jt,p,q,r,it,b,v) &
!$omp&  private(bnP,Cg,ncol,it0,it1,Wi,Wj,tqb)
allocate(buf(maxdi*maxdi*maxdk), Jt(maxdi*maxdi,nb))
allocate(bnP(maxdi*maxdi,nContsAux), Cg(nContsAux,nb), Wi(nop,nb,maxdi), Wj(nop,nb,maxdi))
!$omp do schedule(dynamic)
do ish = 0,nBases-1
   di = df_shell_ncgto(ish); offi = df_shell_offset(ish)
   Wi(:,:,1:di) = 0.0d0
   it0 = iptr(ish)+1
   do while (it0 .le. iptr(ish+1))
      jsh = df_triple_jsh(ilist(it0))
      it1 = it0
      do while (it1+1 .le. iptr(ish+1))
         if (df_triple_jsh(ilist(it1+1)) .ne. jsh) exit
         it1 = it1 + 1
      enddo
      dj = df_shell_ncgto(jsh); offj = df_shell_offset(jsh)
      tqb = omp_get_wtime()
      ncol = 0
      do it = it0,it1
         ksh = df_triple_ksh(ilist(it))
         if (schwarz_bound(ish,jsh)*cbound(ksh-nBases+1) .lt. density_screen_cutoff) cycle
         dk = df_shell_ncgto(ksh); offP = df_shell_offset(ksh)
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         call threec2e_engine(buf(1:di*dj*dk), shls, atm, nAtoms, basDF, nBases+nBasesAux, envDF, int3c2e_opt)
         do r = 1,dk
            do q = 1,dj
               do p = 1,di
                  bnP((q-1)*di+p,ncol+r) = buf((r-1)*di*dj+(q-1)*di+p)*NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
               enddo
            enddo
         enddo
         do b = 1,nb
            Cg(ncol+1:ncol+dk,b) = Cm(offP+1:offP+dk,b)
         enddo
         ncol = ncol + dk
      enddo
      !$omp atomic
      tq(7) = tq(7) + (omp_get_wtime()-tqb)
      tqb = omp_get_wtime()
      if (ncol .gt. 0) then
         call dgemm('N','N',di*dj,nb,ncol,1.0d0,bnP,maxdi*maxdi,Cg,nContsAux,0.0d0,Jt,maxdi*maxdi)
         !$omp atomic
         tq(8) = tq(8) + (omp_get_wtime()-tqb)
         tqb = omp_get_wtime()
         do p = 1,di
            do b = 1,nb
               do q = 1,dj
                  v = Jt((q-1)*di+p,b)
                  Wi(:,b,p) = Wi(:,b,p) + v*CoT(:,offj+q)
               enddo
            enddo
         enddo
         if (jsh .ne. ish) then
            do q = 1,dj
               do b = 1,nb
                  Wj(:,b,q) = 0.0d0
                  do p = 1,di
                     v = Jt((q-1)*di+p,b)
                     Wj(:,b,q) = Wj(:,b,q) + v*CoT(:,offi+p)
                  enddo
               enddo
            enddo
            !$omp atomic
            tq(9) = tq(9) + (omp_get_wtime()-tqb)
            tqb = omp_get_wtime()
            call omp_set_lock(slk(jsh))
            do q = 1,dj
               JB(:,:,offj+q) = JB(:,:,offj+q) + to_q(Wj(:,:,q))
            enddo
            call omp_unset_lock(slk(jsh))
            !$omp atomic
            tq(10) = tq(10) + (omp_get_wtime()-tqb)
         endif
      endif
      it0 = it1 + 1
   enddo
   call omp_set_lock(slk(ish))
   do p = 1,di
      JB(:,:,offi+p) = JB(:,:,offi+p) + to_q(Wi(:,:,p))
   enddo
   call omp_unset_lock(slk(ish))
enddo
!$omp end do
deallocate(buf, Jt, bnP, Cg, Wi, Wj)
!$omp end parallel
do ish = 0,nBases-1
   call omp_destroy_lock(slk(ish))
enddo
!$omp parallel do private(m,b)
do m = 1,nConts_in
   do b = 1,nb
      JC(m,:,b) = from_q(JB(:,b,m))
   enddo
enddo
!$omp end parallel do
call prof_stop("df_direct_jpass")
tq(4) = omp_get_wtime()-tq0
call get_environment_variable("ENGINE_DFJ_BATCH_VERBOSE", venv)
if (trim(venv) .eq. "1") print '(A,I0,A,4F8.3,A,2F9.2)', " DFJC direct nb=", nb, "  setup/b-pass/solve/j-pass ", tq(1:4), &
   "  thread-s Dt / int+dgemm ", tq(5), tq(6)
if (trim(venv) .eq. "1") print '(A,4F9.2)', "   j-pass thread-s int / dgemm / W loops / lock+add ", tq(7:10)
if (trim(venv) .eq. "1") call flush(6)
deallocate(ord, wrk)
deallocate(LB, CoT, JB, Lnrm, Cnrm, smax, Bm, Tm, Cm, cbound, cnt, iptr, ilist, slk)
contains
subroutine sort_desc(n, w, idx)
integer,intent(in) :: n
real(8),intent(inout) :: w(0:n-1)
integer,intent(inout) :: idx(0:n-1)
integer :: a, c, ti
real(8) :: tw
do a = 1,n-1
   tw = w(a) ; ti = idx(a) ; c = a - 1
   do while (c .ge. 0)
      if (w(c) .ge. tw) exit
      w(c+1) = w(c) ; idx(c+1) = idx(c) ; c = c - 1
   enddo
   w(c+1) = tw ; idx(c+1) = ti
enddo
end subroutine sort_desc
end subroutine df_build_coulomb_direct_batch_c

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
