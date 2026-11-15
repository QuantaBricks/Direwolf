! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Density-fitting (RI-J) Coulomb matrix build.

submodule (mod_integrals) df_coulomb_impl
use iso_c_binding, only: c_int, c_long, c_double
implicit none
real(8),allocatable :: dfj_cmp_J(:,:), dfj_cmp_b(:)
integer(8) :: mdj_ghs = -1
integer :: mdj_npair = -1
integer,allocatable :: mdj_pish(:), mdj_pjsh(:)
interface
   function mdj_init(atm_, bas_, nbas_, norb_, env_, pure_) bind(C, name='mdj_init')
      import :: c_int, c_long, c_double
      integer(c_int) :: atm_(*), bas_(*), nbas_, norb_, pure_
      real(c_double) :: env_(*)
      integer(c_long) :: mdj_init
   end function
   subroutine mdj_bfinal(G_, b_) bind(C, name='mdj_bfinal')
      import :: c_double
      real(c_double) :: G_(*), b_(*)
   end subroutine
   subroutine mdj_jprep(c_, H_) bind(C, name='mdj_jprep')
      import :: c_double
      real(c_double) :: c_(*), H_(*)
   end subroutine
   subroutine mdj_grad_pair(ish_, jsh_, D_, H_, schw_, auxb_, cut0_, gA_, gB_, gC_) bind(C, name='mdj_grad_pair')
      import :: c_int, c_double
      integer(c_int) :: ish_, jsh_
      real(c_double) :: D_(*), H_(*), schw_, auxb_(*), cut0_, gA_(3), gB_(3), gC_(*)
   end subroutine
   subroutine mdj_bpass_pair_g(ish_, jsh_, D_, schw_, auxb_, cut_, G_) bind(C, name='mdj_bpass_pair')
      import :: c_int, c_double
      integer(c_int) :: ish_, jsh_
      real(c_double) :: D_(*), schw_, auxb_(*), cut_, G_(*)
   end subroutine
   subroutine mdj_b1_pair(ish_, jsh_, D_, schw_, auxb_, cut_, ldG_, G0_, GI_, GJ_, GK_) bind(C, name='mdj_b1_pair')
      import :: c_int, c_long, c_double
      integer(c_int) :: ish_, jsh_
      integer(c_long) :: ldG_
      real(c_double) :: D_(*), schw_, auxb_(*), cut_, G0_(*), GI_(*), GJ_(*), GK_(*)
   end subroutine
   subroutine mdj_hess_pair(ish_, jsh_, D_, H_, schw_, auxb_, cut_, nat_, aI_, aJ_, d2_) bind(C, name='mdj_hess_pair')
      import :: c_int, c_double
      integer(c_int) :: ish_, jsh_, nat_, aI_, aJ_
      real(c_double) :: D_(*), H_(*), schw_, auxb_(*), cut_, d2_(*)
   end subroutine
   subroutine mdj_jx_pair(ish_, jsh_, H_, schw_, auxb_, cut_, nat_, SA_, SB_, nK_, Klist_, SK_) bind(C, name='mdj_jx_pair')
      import :: c_int, c_double
      integer(c_int) :: ish_, jsh_, nat_, nK_, Klist_(*)
      real(c_double) :: H_(*), schw_, auxb_(*), cut_, SA_(*), SB_(*), SK_(*)
   end subroutine
   subroutine mdj_bpass_multi(ish_, jsh_, D_, ldD_, nb_, schw_, auxb_, cut_, G_, ldG_) bind(C, name='mdj_bpass_multi')
      import :: c_int, c_long, c_double
      integer(c_int) :: ish_, jsh_, ldD_, nb_
      integer(c_long) :: ldG_
      real(c_double) :: D_(*), schw_, auxb_(*), cut_, G_(*)
   end subroutine
   subroutine mdj_jpass_multi(ish_, jsh_, H_, ldH_, nb_, schw_, auxb_, cut0_, cb_, cut_, J_, ldJ_) &
         bind(C, name='mdj_jpass_multi')
      import :: c_int, c_long, c_double
      integer(c_int) :: ish_, jsh_, nb_, ldJ_
      integer(c_long) :: ldH_
      real(c_double) :: H_(*), schw_, auxb_(*), cut0_, cb_(*), cut_, J_(*)
   end subroutine
end interface
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
   block
      logical :: ext_done
      ext_done = .false.
      if (.not. ext_done) call df_build_coulomb_core(nConts_in, Ptot, J)
   end block
   incr_since_full = 0
else
   allocate(dP(nConts_in,nConts_in), dJ(nConts_in,nConts_in))
   dP = Ptot - incr_Ptot_prev
   block
      logical :: ext_done
      ext_done = .false.
      if (.not. ext_done) call df_build_coulomb_core(nConts_in, dP, dJ, incremental_call=.true.)
   end block
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
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts, nbatch
real(8),intent(in) :: Ptot(nConts,nConts,nbatch)
real(8),intent(out) :: J(nConts,nConts,nbatch)
integer :: b, slot
real(8),allocatable :: Ptot_compact(:,:), J_compact(:,:)
real(8),allocatable :: Bvec(:,:), Tvec(:,:), Cvec(:,:)
external :: dgemm

call prof_start("df_setup")
call ensure_df_built()
call prof_stop("df_setup")

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
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(out) :: J(nConts_in,nConts_in)
logical,intent(in),optional :: incremental_call
real(8),allocatable :: bvec(:), tvec(:), cvec(:), bprime(:), Ptot_compact(:), J_compact(:)
integer :: slot

call prof_start("df_setup")
call ensure_df_built()
call prof_stop("df_setup")

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

if (df_triple_defer) return

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
use mod_profile, only: prof_start, prof_stop
implicit none
integer,intent(in) :: nConts, nb, nop
real(8),intent(in) :: Lf(nConts,nop,nb), Co(nConts,nop)
real(8),intent(out) :: JC(nConts,nop,nb)
real(8),allocatable :: Pc(:,:), Jc2(:,:), Bvec(:,:), Tvec(:,:), Cvec(:,:), Cg(:,:), Lm(:,:), Tm(:,:)
real(8),allocatable :: Dch(:,:,:), Jch(:,:,:)
integer :: m, s0, cnt, b, k, sl, sm, o, c0, nch
real(8) :: t
external :: dgemm

call prof_start("df_setup")
call ensure_df_built()
call prof_stop("df_setup")
if (df_direct_mode .and. .not. dfK_built) then
   block
      character(len=8) :: e
      logical :: dense
      call get_environment_variable("ENGINE_DFJ_BATCHC_DENSE", e)
      dense = (trim(e) .eq. "1")
      if (.not. dense) then
         call df_build_coulomb_direct_batch_c(nConts, nb, Lf, nop, Co, JC)
         return
      endif
   end block
endif
if (df_direct_mode .or. dfK_built) then
   allocate(Dch(nConts,nConts,16), Jch(nConts,nConts,16))
   do c0 = 0, nb-1, 16
      nch = min(16, nb-c0)
      !$omp parallel do private(b) schedule(static)
      do b = 1,nch
         call dgemm('N','T',nConts,nConts,nop,1.0d0,Lf(1,1,c0+b),nConts,Co,nConts,0.0d0,Dch(1,1,b),nConts)
         Dch(:,:,b) = Dch(:,:,b) + transpose(Dch(:,:,b))
      enddo
      !$omp end parallel do
      call integrals_build_coulomb_df_batch(nConts, nch, Dch(:,:,1:nch), Jch(:,:,1:nch))
      !$omp parallel do private(b) schedule(static)
      do b = 1,nch
         call dgemm('N','N',nConts,nop,nConts,1.0d0,Jch(1,1,b),nConts,Co,nConts,0.0d0,JC(1,1,c0+b),nConts)
      enddo
      !$omp end parallel do
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

if (.not. allocated(df_triple_ish)) call build_df_triple_list()
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
if (engine_verbose .ge. 2) then
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
logical :: md
integer(c_int) :: ishc, jshc, nbc, ldc
real(8),allocatable :: Gm(:,:), Gloc(:,:), Hm(:,:), Dn(:,:), cNa(:)
real(8) :: cut0

tq = 0.0d0 ; tq0 = omp_get_wtime()
maxdi = maxval(df_shell_ncgto(0:nBases-1))
maxdk = maxval(df_shell_ncgto(nBases:nBases+nBasesAux-1))
maxdiag = maxval(aux_shell_bound)
md = (dfj_md_mode() == 1)
if (md) then
   block
      character(len=8) :: nbe
      integer :: nbmax_md, nbios
      nbmax_md = 40
      call get_environment_variable("ENGINE_DFJ_MD_NBMAX", nbe)
      if (len_trim(nbe) > 0) then
         read(nbe, *, iostat=nbios) nbmax_md
      endif
      if (nb > nbmax_md) md = .false.
   end block
endif
cut0 = SCHWARZ_CUTOFF
nbc = nb ; ldc = maxdi*maxdi
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

if (.not. allocated(df_triple_ish)) call build_df_triple_list()
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
if (md) then
   allocate(Gm(mdj_ghs,nb))
   Gm = 0.0d0
endif
!$omp parallel private(ii,ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf,Dt,Bloc,p,q,r,it,b,fd,sd,tqa,r0) &
!$omp&  private(bnP,Tb,poff,pdk,np,ncol,it0,it1,s0,Gloc,Dn,ishc,jshc)
allocate(buf(maxdi*maxdi*maxdk), Dt(maxdi*maxdi,nb), Bloc(nContsAux,nb))
allocate(bnP(maxdi*maxdi,nContsAux), Tb(nContsAux,nb), poff(nBasesAux), pdk(nBasesAux))
if (md) then
   allocate(Gloc(mdj_ghs,nb), Dn(maxdi*maxdi,nb))
   Gloc = 0.0d0
endif
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
      if (md) then
         do b = 1,nb
            do q = 1,dj
               do p = 1,di
                  Dn((q-1)*di+p,b) = Dt((q-1)*di+p,b)*NorVEC(offi+p)*NorVEC(offj+q)
               enddo
            enddo
         enddo
         ishc = ish ; jshc = jsh
         call mdj_bpass_multi(ishc, jshc, Dn, ldc, nbc, schwarz_bound(ish,jsh), aux_shell_bound, cut0, Gloc, mdj_ghs)
         !$omp atomic
         tq(6) = tq(6) + (omp_get_wtime()-tqa)
         it0 = it1 + 1
         cycle
      endif
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
      if (omp_get_thread_num() .eq. mt) then
         if (md) then
            Gm = Gm + Gloc
         else
            Bm = Bm + Bloc
         endif
      endif
   enddo
end block
deallocate(buf, Dt, Bloc, bnP, Tb, poff, pdk)
if (md) deallocate(Gloc, Dn)
!$omp end parallel
if (md) then
   do b = 1,nb
      call mdj_bfinal(Gm(1,b), Bm(1,b))
      Bm(:,b) = Bm(:,b)*NorVECAux
   enddo
   deallocate(Gm)
endif
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
if (md) then
   allocate(Hm(mdj_ghs,nb), cNa(nContsAux))
   do b = 1,nb
      cNa = Cm(:,b)*NorVECAux
      call mdj_jprep(cNa, Hm(1,b))
   enddo
endif
!$omp parallel private(ish,jsh,ksh,di,dj,dk,offi,offj,offP,shls,buf,Jt,p,q,r,it,b,v) &
!$omp&  private(bnP,Cg,ncol,it0,it1,Wi,Wj,tqb,ishc,jshc)
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
      if (md) then
         ishc = ish ; jshc = jsh
         call mdj_jpass_multi(ishc, jshc, Hm, mdj_ghs, nbc, schwarz_bound(ish,jsh), aux_shell_bound, cut0, &
                              cbound, density_screen_cutoff, Jt, ldc)
         do b = 1,nb
            do q = 1,dj
               do p = 1,di
                  Jt((q-1)*di+p,b) = Jt((q-1)*di+p,b)*NorVEC(offi+p)*NorVEC(offj+q)
               enddo
            enddo
         enddo
         ncol = -1
         it = it1 + 1
      endif
      do it = it0,it1
         if (md) exit
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
      if (ncol .ne. 0) then
         if (ncol .gt. 0) call dgemm('N','N',di*dj,nb,ncol,1.0d0,bnP,maxdi*maxdi,Cg,nContsAux,0.0d0,Jt,maxdi*maxdi)
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
if (engine_verbose .ge. 2) print '(A,I0,A,4F8.3,A,2F9.2)', " DFJC direct nb=", nb, "  setup/b-pass/solve/j-pass ", tq(1:4), &
   "  thread-s Dt / int+dgemm ", tq(5), tq(6)
if (engine_verbose .ge. 2) print '(A,4F9.2)', "   j-pass thread-s int / dgemm / W loops / lock+add ", tq(7:10)
if (engine_verbose .ge. 2) call flush(6)
deallocate(ord, wrk)
if (md) deallocate(Hm, cNa)
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

if (dfj_md_mode() > 0) then
   block
      real(8),allocatable :: Jm(:,:), bm(:)
      if (dfj_md_mode() == 1) then
         call df_coulomb_direct_md(nConts_in, Ptot, J, this_cutoff, dmax_shell, maxdiag)
         deallocate(bvec, tvec, cvec, dmax_shell)
         return
      endif
      allocate(Jm(nConts_in,nConts_in))
      call df_coulomb_direct_md(nConts_in, Ptot, Jm, this_cutoff, dmax_shell, maxdiag, bm)
      dfj_cmp_J = Jm ; dfj_cmp_b = bm
   end block
endif
if (.not. allocated(df_triple_ish)) call build_df_triple_list()
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

if (dfj_md_mode() == 2) then
   print '(A,2ES12.3,A,2ES12.3)', '[dfj md] max|db|, max|b| ', maxval(abs(dfj_cmp_b - bvec)), maxval(abs(bvec)), &
         '  max|dJ|, max|J| ', maxval(abs(dfj_cmp_J - J)), maxval(abs(J))
   deallocate(dfj_cmp_J, dfj_cmp_b)
endif
deallocate(bvec,tvec,cvec)
deallocate(dmax_shell,cvec_shell_bound)
end subroutine df_build_coulomb_direct

integer function dfj_md_mode()
use iso_c_binding, only: c_int
integer, save :: mode = -1
character(len=8) :: e
integer(c_int) :: nb_all, norb_c, pure_c
if (mode < 0) then
   call get_environment_variable("ENGINE_DFJ_MD", e)
   mode = 1
   if (trim(e) == "0") mode = 0
   if (trim(e) == "2") mode = 2
endif
if (mode > 0 .and. mdj_stale) then
   nb_all = nBases + nBasesAux ; norb_c = nBases ; pure_c = merge(1, 0, engine_puream)
   mdj_ghs = mdj_init(atm, basDF, nb_all, norb_c, envDF, pure_c)
   mdj_npair = -1
   mdj_stale = .false.
endif
dfj_md_mode = mode
if (mode > 0 .and. mdj_ghs < 0) dfj_md_mode = 0
end function dfj_md_mode

subroutine df_coulomb_direct_md(nConts_in, Ptot, J, this_cutoff, dmax_shell, maxdiag, bout)
use MOL_info
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in), this_cutoff, dmax_shell(0:nBases-1,0:nBases-1), maxdiag
real(8),intent(out) :: J(nConts_in,nConts_in)
real(8),allocatable,intent(out),optional :: bout(:)
interface
   subroutine mdj_bpass_pair(ish_, jsh_, D_, schw_, auxb_, cut_, G_) bind(C, name='mdj_bpass_pair')
      import :: c_int, c_double
      integer(c_int) :: ish_, jsh_
      real(c_double) :: D_(*), schw_, auxb_(*), cut_, G_(*)
   end subroutine
   subroutine mdj_jpass_pair(ish_, jsh_, H_, schw_, auxb_, cut0_, cb_, cut_, Js_) bind(C, name='mdj_jpass_pair')
      import :: c_int, c_double
      integer(c_int) :: ish_, jsh_
      real(c_double) :: H_(*), schw_, auxb_(*), cut0_, cb_(*), cut_, Js_(*)
   end subroutine
end interface
real(8),allocatable :: G(:), Gl(:,:), bvec(:), tvec(:), cvec(:), cN(:), H(:), cb(:), Jl(:,:,:), Dblk(:,:), Jblk(:,:)
integer(c_int) :: ish, jsh
integer :: ip, di, dj, offi, offj, p, q, tid, ksh, offP, dk, mt
real(8) :: schw, cut0

cut0 = SCHWARZ_CUTOFF
call mdj_pairs_ready(maxdiag)
mt = omp_get_max_threads()
allocate(G(mdj_ghs), Gl(mdj_ghs, mt), bvec(nContsAux), tvec(nContsAux), cvec(nContsAux))
G = 0.0d0

call prof_start("df_direct_bpass")
!$omp parallel private(ip, ish, jsh, di, dj, offi, offj, p, q, tid, schw, Dblk)
tid = omp_get_thread_num() + 1
Gl(:,tid) = 0.0d0
allocate(Dblk(maxval(df_shell_ncgto(0:nBases-1)), maxval(df_shell_ncgto(0:nBases-1))))
!$omp do schedule(static,1)
do ip = 1, mdj_npair
   ish = mdj_pish(ip) ; jsh = mdj_pjsh(ip)
   schw = schwarz_bound(ish,jsh)
   if (schw*maxdiag*dmax_shell(ish,jsh) < this_cutoff) cycle
   di = df_shell_ncgto(ish) ; dj = df_shell_ncgto(jsh)
   offi = df_shell_offset(ish) ; offj = df_shell_offset(jsh)
   do q = 1, dj
      do p = 1, di
         if (ish /= jsh) then
            Dblk(p+(q-1)*di,1) = (Ptot(offi+p,offj+q) + Ptot(offj+q,offi+p))*NorVEC(offi+p)*NorVEC(offj+q)
         else
            Dblk(p+(q-1)*di,1) = Ptot(offi+p,offj+q)*NorVEC(offi+p)*NorVEC(offj+q)
         endif
      enddo
   enddo
   call mdj_bpass_pair(ish, jsh, Dblk, schw, aux_shell_bound, cut0, Gl(1,tid))
enddo
!$omp end do
deallocate(Dblk)
block
   integer :: m
   do m = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() == m) G = G + Gl(:,tid)
   enddo
end block
!$omp end parallel
bvec = 0.0d0
call mdj_bfinal(G, bvec)
bvec = bvec*NorVECAux
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
allocate(cb(nBasesAux), cN(nContsAux), H(mdj_ghs))
offP = 0
do ksh = nBases, nBases+nBasesAux-1
   dk = df_shell_ncgto(ksh)
   cb(ksh-nBases+1) = maxval(abs(cvec(offP+1:offP+dk)))
   offP = offP + dk
enddo
cN = cvec*NorVECAux
call mdj_jprep(cN, H)
call prof_stop("df_direct_solve")

call prof_start("df_direct_jpass")
J = 0.0d0
!$omp parallel private(ip, ish, jsh, di, dj, offi, offj, p, q, schw, Jblk)
allocate(Jblk(maxval(df_shell_ncgto(0:nBases-1)), maxval(df_shell_ncgto(0:nBases-1))))
!$omp do schedule(static,1)
do ip = 1, mdj_npair
   ish = mdj_pish(ip) ; jsh = mdj_pjsh(ip)
   schw = schwarz_bound(ish,jsh)
   di = df_shell_ncgto(ish) ; dj = df_shell_ncgto(jsh)
   offi = df_shell_offset(ish) ; offj = df_shell_offset(jsh)
   call mdj_jpass_pair(ish, jsh, H, schw, aux_shell_bound, cut0, cb, this_cutoff, Jblk)
   do q = 1, dj
      do p = 1, di
         J(offi+p,offj+q) = Jblk(p+(q-1)*di,1)*NorVEC(offi+p)*NorVEC(offj+q)
      enddo
   enddo
   if (ish /= jsh) then
      do q = 1, dj
         do p = 1, di
            J(offj+q,offi+p) = Jblk(p+(q-1)*di,1)*NorVEC(offi+p)*NorVEC(offj+q)
         enddo
      enddo
   endif
enddo
!$omp end do
deallocate(Jblk)
!$omp end parallel
call prof_stop("df_direct_jpass")
if (present(bout)) then
   allocate(bout(nContsAux)) ; bout = bvec
endif
deallocate(G, Gl, bvec, tvec, cvec, cN, H, cb)

end subroutine df_coulomb_direct_md

module subroutine df_mdj_force_j(nConts_in, Ptot, natoms_in, force3c, cvec, ok)
use MOL_info
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts_in, natoms_in
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
real(8),intent(inout) :: force3c(natoms_in,3)
real(8),intent(out) :: cvec(*)
logical,intent(out) :: ok
real(8),allocatable :: G(:), Gl(:,:), bvec(:), tvec(:), cN(:), H(:), gCl(:,:,:), gC(:,:), Dblk(:)
real(8),allocatable :: fl(:,:,:)
integer(c_int) :: ish, jsh
integer :: ip, di, dj, offi, offj, p, q, tid, mt, ia, ja, m
real(8) :: schw, cut0, maxdiag, gA(3), gB(3)

ok = (dfj_md_mode() >= 1)
if (.not. ok) return
maxdiag = maxval(aux_shell_bound)
call mdj_pairs_ready(maxdiag)
cut0 = SCHWARZ_CUTOFF
mt = omp_get_max_threads()
allocate(cN(nContsAux), H(mdj_ghs))
if (df_direct_cvec_cache_valid .and. all(df_direct_cvec_cache_Ptot .eq. Ptot)) then
   cvec(1:nContsAux) = df_direct_cvec_cache
else
   call prof_start("df_force_direct_bpass")
   allocate(G(mdj_ghs), Gl(mdj_ghs,mt), bvec(nContsAux), tvec(nContsAux))
   G = 0.0d0
!$omp parallel private(ip, ish, jsh, di, dj, offi, offj, p, q, tid, schw, Dblk)
   tid = omp_get_thread_num() + 1
   Gl(:,tid) = 0.0d0
   allocate(Dblk(maxval(df_shell_ncgto(0:nBases-1))**2))
!$omp do schedule(static,1)
   do ip = 1, mdj_npair
      ish = mdj_pish(ip) ; jsh = mdj_pjsh(ip)
      schw = schwarz_bound(ish,jsh)
      di = df_shell_ncgto(ish) ; dj = df_shell_ncgto(jsh)
      offi = df_shell_offset(ish) ; offj = df_shell_offset(jsh)
      call pair_dens(Ptot, ish, jsh, di, dj, offi, offj, Dblk)
      call mdj_bpass_pair_g(ish, jsh, Dblk, schw, aux_shell_bound, cut0, Gl(1,tid))
   enddo
!$omp end do
   deallocate(Dblk)
   do m = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() == m) G = G + Gl(:,tid)
   enddo
!$omp end parallel
   bvec = 0.0d0
   call mdj_bfinal(G, bvec)
   bvec = bvec*NorVECAux
   call dgemv('T', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, bvec, 1, 0.0d0, tvec, 1)
   tvec = tvec * df_evalinv
   call dgemv('N', nContsAux, nContsAux, 1.0d0, df_evec, nContsAux, tvec, 1, 0.0d0, cvec, 1)
   deallocate(G, Gl, bvec, tvec)
   call prof_stop("df_force_direct_bpass")
endif
cN = cvec(1:nContsAux)*NorVECAux
call mdj_jprep(cN, H)

call prof_start("df_force_direct_contract")
allocate(gC(3,nAtoms), gCl(3,nAtoms,mt), fl(natoms_in,3,mt))
gC = 0.0d0
!$omp parallel private(ip, ish, jsh, di, dj, offi, offj, p, q, tid, schw, Dblk, gA, gB, ia, ja)
tid = omp_get_thread_num() + 1
gCl(:,:,tid) = 0.0d0
fl(:,:,tid) = 0.0d0
allocate(Dblk(maxval(df_shell_ncgto(0:nBases-1))**2))
!$omp do schedule(dynamic)
do ip = 1, mdj_npair
   ish = mdj_pish(ip) ; jsh = mdj_pjsh(ip)
   schw = schwarz_bound(ish,jsh)
   di = df_shell_ncgto(ish) ; dj = df_shell_ncgto(jsh)
   offi = df_shell_offset(ish) ; offj = df_shell_offset(jsh)
   call pair_dens(Ptot, ish, jsh, di, dj, offi, offj, Dblk)
   call mdj_grad_pair(ish, jsh, Dblk, H, schw, aux_shell_bound, cut0, gA, gB, gCl(1,1,tid))
   ia = basDF(1,ish+1) + 1 ; ja = basDF(1,jsh+1) + 1
   fl(ia,:,tid) = fl(ia,:,tid) - gA
   fl(ja,:,tid) = fl(ja,:,tid) - gB
enddo
!$omp end do
deallocate(Dblk)
!$omp end parallel
do m = 1, mt
   force3c = force3c + fl(:,:,m)
   gC = gC + gCl(:,:,m)
enddo
do ia = 1, min(nAtoms, natoms_in)
   force3c(ia,:) = force3c(ia,:) - gC(:,ia)
enddo
call prof_stop("df_force_direct_contract")
deallocate(cN, H, gC, gCl, fl)

contains

subroutine pair_dens(Pt, is, js, ni, nj, oi, oj, Db)
real(8),intent(in) :: Pt(nConts_in,nConts_in)
integer(c_int),intent(in) :: is, js
integer,intent(in) :: ni, nj, oi, oj
real(8),intent(out) :: Db(*)
integer :: a, b
do b = 1, nj
   do a = 1, ni
      if (is /= js) then
         Db(a+(b-1)*ni) = (Pt(oi+a,oj+b) + Pt(oj+b,oi+a))*NorVEC(oi+a)*NorVEC(oj+b)
      else
         Db(a+(b-1)*ni) = Pt(oi+a,oj+b)*NorVEC(oi+a)*NorVEC(oj+b)
      endif
   enddo
enddo
end subroutine pair_dens

end subroutine df_mdj_force_j

subroutine mdj_pairs_ready(maxdiag)
real(8),intent(in) :: maxdiag
integer(c_int) :: ish, jsh
if (mdj_npair >= 0) return
if (allocated(mdj_pish)) deallocate(mdj_pish, mdj_pjsh)
mdj_npair = 0
do ish = 0, nBases-1
   do jsh = 0, ish
      if (schwarz_bound(ish,jsh)*maxdiag >= SCHWARZ_CUTOFF) mdj_npair = mdj_npair + 1
   enddo
enddo
allocate(mdj_pish(mdj_npair), mdj_pjsh(mdj_npair))
mdj_npair = 0
do ish = 0, nBases-1
   do jsh = 0, ish
      if (schwarz_bound(ish,jsh)*maxdiag >= SCHWARZ_CUTOFF) then
         mdj_npair = mdj_npair + 1 ; mdj_pish(mdj_npair) = ish ; mdj_pjsh(mdj_npair) = jsh
      endif
   enddo
enddo
end subroutine mdj_pairs_ready

module subroutine df_mdj_hess_b1(nConts_in, natoms_in, Ptot, own, nth, bvec, bx, ok)
use MOL_info
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
use iso_c_binding, only: c_int, c_long, c_double
implicit none
integer,intent(in) :: nConts_in, natoms_in, nth
real(8),intent(in) :: Ptot(nConts_in,nConts_in)
integer,intent(in) :: own(0:nBases-1,0:nBases-1)
real(8),intent(out) :: bvec(nContsAux), bx(natoms_in,3,nContsAux)
logical,intent(out) :: ok
real(8),allocatable :: G0t(:,:), Gat(:,:,:,:), GKt(:,:,:), G0(:), Ga(:,:,:), GK(:,:), Dblk(:), cb(:), tmp(:)
real(8),allocatable :: dmax_shell(:,:)
integer(c_int) :: ish, jsh
integer :: ip, di, dj, offi, offj, tid, atI, atJ, it, ia, ix, ik, nshp, io
integer(c_long) :: ldG
real(8) :: schw, cutv, maxdiag, gbmax
character(len=16) :: ge

ok = (dfj_md_mode() >= 1)
if (.not. ok) return
gbmax = 3.0d0
call get_environment_variable("ENGINE_DFJ_MD_HESS_GB", ge)
if (len_trim(ge) > 0) read(ge, *, iostat=io) gbmax
if (8.0d0*(3.0d0*natoms_in + 4.0d0)*mdj_ghs*nth > gbmax*1.0d9) then
   ok = .false.
   return
endif
maxdiag = maxval(aux_shell_bound)
call mdj_pairs_ready(maxdiag)
ldG = mdj_ghs
allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts_in, Ptot, dmax_shell)
allocate(G0t(mdj_ghs,nth), Gat(mdj_ghs,3,natoms_in,nth), GKt(mdj_ghs,3,nth))
G0t = 0.0d0 ; Gat = 0.0d0 ; GKt = 0.0d0
!$omp parallel default(shared) num_threads(nth) private(ip, ish, jsh, di, dj, offi, offj, tid, atI, atJ, schw, cutv, Dblk)
tid = omp_get_thread_num() + 1
if (omp_get_num_threads() .ne. nth) error stop 'df_mdj_hess_b1: thread team smaller than planned'
allocate(Dblk(maxval(df_shell_ncgto(0:nBases-1))**2))
do ip = 1, mdj_npair
   ish = mdj_pish(ip) ; jsh = mdj_pjsh(ip)
   if (own(ish,jsh) .ne. tid) cycle
   schw = schwarz_bound(ish,jsh)
   if (dmax_shell(ish,jsh) <= 0.0d0) cycle
   cutv = max(SCHWARZ_CUTOFF, density_screen_cutoff/dmax_shell(ish,jsh))
   di = df_shell_ncgto(ish) ; dj = df_shell_ncgto(jsh)
   offi = df_shell_offset(ish) ; offj = df_shell_offset(jsh)
   atI = basDF(1,ish+1) + 1 ; atJ = basDF(1,jsh+1) + 1
   call mdj_pair_dens(Ptot, ish, jsh, di, dj, offi, offj, Dblk, nConts_in)
   call mdj_b1_pair(ish, jsh, Dblk, schw, aux_shell_bound, cutv, ldG, G0t(1,tid), Gat(1,1,atI,tid), &
                    Gat(1,1,atJ,tid), GKt(1,1,tid))
enddo
deallocate(Dblk)
!$omp end parallel
allocate(G0(mdj_ghs), Ga(mdj_ghs,3,natoms_in), GK(mdj_ghs,3), cb(nContsAux), tmp(nContsAux))
G0 = 0.0d0 ; Ga = 0.0d0 ; GK = 0.0d0
do tid = 1, nth
   G0 = G0 + G0t(:,tid)
   Ga = Ga + Gat(:,:,:,tid)
   GK = GK + GKt(:,:,tid)
enddo
deallocate(G0t, Gat, GKt)
bvec = 0.0d0
call mdj_bfinal(G0, bvec)
bvec = bvec*NorVECAux
bx = 0.0d0
do ia = 1, natoms_in
   do ix = 1, 3
      tmp = 0.0d0
      call mdj_bfinal(Ga(1,ix,ia), tmp)
      bx(ia,ix,:) = tmp*NorVECAux
   enddo
enddo
do ix = 1, 3
   tmp = 0.0d0
   call mdj_bfinal(GK(1,ix), tmp)
   tmp = tmp*NorVECAux
   io = 0
   do ik = nBases, nBases+nBasesAux-1
      nshp = df_shell_ncgto(ik)
      ia = basDF(1,ik+1) + 1
      do it = 1, nshp
         bx(ia,ix,io+it) = bx(ia,ix,io+it) + tmp(io+it)
      enddo
      io = io + nshp
   enddo
enddo
deallocate(G0, Ga, GK, cb, tmp, dmax_shell)
end subroutine df_mdj_hess_b1

module subroutine df_mdj_hess_pass3(nConts_in, natoms_in, Ptot, cvec, own, nth, d2E_t, ok)
use MOL_info
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
use iso_c_binding, only: c_int, c_long, c_double
implicit none
integer,intent(in) :: nConts_in, natoms_in, nth
real(8),intent(in) :: Ptot(nConts_in,nConts_in), cvec(nContsAux)
integer,intent(in) :: own(0:nBases-1,0:nBases-1)
real(8),intent(inout) :: d2E_t(3,natoms_in,3,natoms_in,nth)
logical,intent(out) :: ok
real(8),allocatable :: H(:), cN(:), Dblk(:), dmax_shell(:,:)
integer(c_int) :: ish, jsh, nat_c, aI, aJ
integer :: ip, di, dj, offi, offj, tid
real(8) :: schw, cutv, maxdiag

ok = (dfj_md_mode() >= 1)
if (.not. ok) return
maxdiag = maxval(aux_shell_bound)
call mdj_pairs_ready(maxdiag)
allocate(H(mdj_ghs), cN(nContsAux), dmax_shell(0:nBases-1,0:nBases-1))
cN = cvec*NorVECAux
call mdj_jprep(cN, H)
call compute_shell_density_bound(nConts_in, Ptot, dmax_shell)
nat_c = natoms_in
!$omp parallel default(shared) num_threads(nth) private(ip, ish, jsh, di, dj, offi, offj, tid, aI, aJ, schw, cutv, Dblk)
tid = omp_get_thread_num() + 1
if (omp_get_num_threads() .ne. nth) error stop 'df_mdj_hess_pass3: thread team smaller than planned'
allocate(Dblk(maxval(df_shell_ncgto(0:nBases-1))**2))
do ip = 1, mdj_npair
   ish = mdj_pish(ip) ; jsh = mdj_pjsh(ip)
   if (own(ish,jsh) .ne. tid) cycle
   schw = schwarz_bound(ish,jsh)
   if (dmax_shell(ish,jsh) <= 0.0d0) cycle
   cutv = max(SCHWARZ_CUTOFF, density_screen_cutoff/dmax_shell(ish,jsh))
   di = df_shell_ncgto(ish) ; dj = df_shell_ncgto(jsh)
   offi = df_shell_offset(ish) ; offj = df_shell_offset(jsh)
   aI = basDF(1,ish+1) ; aJ = basDF(1,jsh+1)
   call mdj_pair_dens(Ptot, ish, jsh, di, dj, offi, offj, Dblk, nConts_in)
   call mdj_hess_pair(ish, jsh, Dblk, H, schw, aux_shell_bound, cutv, nat_c, aI, aJ, d2E_t(1,1,1,1,tid))
enddo
deallocate(Dblk)
!$omp end parallel
deallocate(H, cN, dmax_shell)
end subroutine df_mdj_hess_pass3

subroutine mdj_pair_dens(Pt, is, js, ni, nj, oi, oj, Db, n)
integer,intent(in) :: n
real(8),intent(in) :: Pt(n,n)
integer(c_int),intent(in) :: is, js
integer,intent(in) :: ni, nj, oi, oj
real(8),intent(out) :: Db(*)
integer :: a, b
do b = 1, nj
   do a = 1, ni
      if (is /= js) then
         Db(a+(b-1)*ni) = (Pt(oi+a,oj+b) + Pt(oj+b,oi+a))*NorVEC(oi+a)*NorVEC(oj+b)
      else
         Db(a+(b-1)*ni) = Pt(oi+a,oj+b)*NorVEC(oi+a)*NorVEC(oj+b)
      endif
   enddo
enddo
end subroutine mdj_pair_dens

module subroutine df_mdj_jx_prepare(cvec, H, ok)
implicit none
real(8),intent(in) :: cvec(nContsAux)
real(8),allocatable,intent(out) :: H(:)
logical,intent(out) :: ok
real(8),allocatable :: cN(:)
ok = (dfj_md_mode() >= 1)
block
   character(len=4) :: je
   call get_environment_variable("ENGINE_DFJ_MD_JX", je)
   if (trim(je) == "0") ok = .false.
end block
if (.not. ok) return
allocate(H(mdj_ghs), cN(nContsAux))
cN = cvec*NorVECAux
call mdj_jprep(cN, H)
end subroutine df_mdj_jx_prepare

module subroutine df_mdj_jx_pair(ish, jsh, H, nat, SA, SB, nK, Klist, SK)
implicit none
integer,intent(in) :: ish, jsh, nat
real(8),intent(in) :: H(*)
real(8),intent(out) :: SA(*), SB(*), SK(*)
integer,intent(out) :: nK, Klist(nat)
integer(c_int) :: ic, jc, nc, nkc
real(8) :: schw, cut
ic = ish ; jc = jsh ; nc = nat
schw = schwarz_bound(ish,jsh)
cut = SCHWARZ_CUTOFF
call mdj_jx_pair(ic, jc, H, schw, aux_shell_bound, cut, nc, SA, SB, nkc, Klist, SK)
nK = nkc
end subroutine df_mdj_jx_pair

end submodule df_coulomb_impl
