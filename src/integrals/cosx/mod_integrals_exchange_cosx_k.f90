! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! COSX full-range K builder - split out of mod_integrals_exchange_cosx.f90

submodule (mod_integrals:exchange_cosx_impl) exchange_cosx_k_impl
contains
module subroutine integrals_build_exchange_cosx(nConts, Da, Db, Ka, Kb, need_force)
use omp_lib, only: omp_get_wtime, omp_get_thread_num, omp_get_num_threads
use MOL_info, only: engine_verbose
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
logical,intent(in) :: need_force

integer :: ib, nblocks, bstart, nb, g
integer :: si, sj, di, dj, i, j, ao_i, ao_j, maxd, sj_idx
integer :: shls1e(4), off_grids
real(8) :: nij, wij, x_block_max, k_bound
real(8),allocatable :: coor_batch(:,:), w_batch(:)
real(8),allocatable :: Xb(:,:), val1_dummy(:,:,:)
real(8),allocatable :: Fa(:,:), Fb(:,:), Ga(:,:), Gb(:,:)
real(8),allocatable :: env_g(:), buf(:)
real(8),allocatable :: int_cache(:)
integer,parameter :: INT_CACHE_SIZE = 4000000
real(8),allocatable :: Ka_local(:,:), Kb_local(:,:)
real(8),allocatable :: Qb(:,:)
real(8),allocatable :: F_gmax(:)
integer :: nrad_lo, nsph_lo, nrad_md, nsph_md, nrad_hi, nsph_hi, npts, envstat
logical :: nrad_lo_explicit, nrad_md_explicit, nrad_hi_explicit
character(len=16) :: envchar
real(8),pointer :: coor_p(:,:), w_p(:), Qfull_p(:,:)
integer,pointer :: cosx_blkptr_p(:), cosx_blkidx_p(:)
integer,pointer :: blkstart_p(:)
real(8),allocatable :: D_pair_max(:,:)
real(8),allocatable :: D_weighted_max(:)
real(8) :: pair_bound
integer,allocatable :: ao2shell(:)
integer :: nu, mu, ii
real(8) :: dval
real(8) :: local_D_max
integer :: kk
integer :: stage_now
logical :: stage_changed
real(8),allocatable :: Da_eff(:,:), Db_eff(:,:)
logical :: do_incremental
real(8) :: kscreen_use
integer(8) :: dbg_t0, dbg_t1, dbg_tick_local
real(8) :: xmax_this_call, xmax_local
integer :: skip_local, skip_total
integer(8) :: n_geo_local, n_schwarz_local, n_calls_local
integer(8) :: n_geo_total, n_schwarz_total, n_calls_total
integer :: skip_pts_local, skip_pts_total
real(8) :: global_D_max, block_bound
real(8),pointer :: cosx_Xblock_p(:)
logical :: is_full_hi_build
real(8),pointer :: blkcen_p(:,:), blkrad_p(:)
real(8) :: esp_eff, r_far, seg_ab(3), seg_len2, seg_t, seg_cl(3)
integer,allocatable :: loc_ao(:)
integer,allocatable :: run_a(:), run_o(:), run_len(:)
logical,allocatable :: tau_hit(:)
integer,allocatable :: trun_a(:), trun_len(:)
integer :: ntrun, it
integer :: nloc, cofs, maxnloc, iao, nrun, ir, a0, dsh, prev_end
real(8) :: bet
real(8),allocatable :: Xc(:,:), Kc(:,:)
real(8),pointer,contiguous :: QfullT_p(:,:)

if (.not. cosx_ready .or. cosx_nBases_cached .ne. nBases) then
   call cosx_build_shared(nConts)

   nrad_hi = 18
   nsph_hi = 194
   call get_environment_variable("ENGINE_COSX_NRAD", envchar, status=envstat)
   nrad_hi_explicit = (envstat .eq. 0)
   if (envstat .eq. 0) read(envchar,*) nrad_hi
   call get_environment_variable("ENGINE_COSX_NSPH", envchar, status=envstat)
   if (envstat .eq. 0) read(envchar,*) nsph_hi

   nrad_md = 14
   nsph_md = 110
   call get_environment_variable("ENGINE_COSX_NRAD_MD", envchar, status=envstat)
   nrad_md_explicit = (envstat .eq. 0)
   if (envstat .eq. 0) read(envchar,*) nrad_md
   call get_environment_variable("ENGINE_COSX_NSPH_MD", envchar, status=envstat)
   if (envstat .eq. 0) read(envchar,*) nsph_md

   nrad_lo = 14
   nsph_lo = 50
   call get_environment_variable("ENGINE_COSX_NRAD_LO", envchar, status=envstat)
   nrad_lo_explicit = (envstat .eq. 0)
   if (envstat .eq. 0) read(envchar,*) nrad_lo
   call get_environment_variable("ENGINE_COSX_NSPH_LO", envchar, status=envstat)
   if (envstat .eq. 0) read(envchar,*) nsph_lo

   if (allocated(cosx_coor_hi)) deallocate(cosx_coor_hi,cosx_w_hi,cosx_Qfull_hi)
   if (allocated(cosx_QfullT_hi)) deallocate(cosx_QfullT_hi)
   if (allocated(cosx_blkptr_hi)) deallocate(cosx_blkptr_hi,cosx_blkidx_hi)
   if (allocated(cosx_blkcen_hi)) deallocate(cosx_blkcen_hi,cosx_blkrad_hi)
   if (allocated(cosx_blkstart_hi)) deallocate(cosx_blkstart_hi)
   if (nrad_hi_explicit) then
      call cosx_build_one_grid(nConts, nrad_hi, nsph_hi, cosx_coor_hi, cosx_w_hi, cosx_npts_hi, cosx_Qfull_hi, &
                                cosx_blkptr_hi, cosx_blkidx_hi, cosx_blkcen_hi, cosx_blkrad_hi, &
                                cosx_QfullT_hi, cosx_blkstart_hi, cosx_nblk_hi, &
                                per_atom_period_scale=.true.)
   else
      call cosx_build_one_grid(nConts, nrad_hi, nsph_hi, cosx_coor_hi, cosx_w_hi, cosx_npts_hi, cosx_Qfull_hi, &
                                cosx_blkptr_hi, cosx_blkidx_hi, cosx_blkcen_hi, cosx_blkrad_hi, &
                                cosx_QfullT_hi, cosx_blkstart_hi, cosx_nblk_hi, &
                                per_atom_period_scale=.true., intacc_eps=4.338d0)
   endif
   if (allocated(cosx_coor_md)) deallocate(cosx_coor_md,cosx_w_md,cosx_Qfull_md)
   if (allocated(cosx_QfullT_md)) deallocate(cosx_QfullT_md)
   if (allocated(cosx_blkptr_md)) deallocate(cosx_blkptr_md,cosx_blkidx_md)
   if (allocated(cosx_blkcen_md)) deallocate(cosx_blkcen_md,cosx_blkrad_md)
   if (allocated(cosx_blkstart_md)) deallocate(cosx_blkstart_md)
   if (nrad_md_explicit) then
      call cosx_build_one_grid(nConts, nrad_md, nsph_md, cosx_coor_md, cosx_w_md, cosx_npts_md, cosx_Qfull_md, &
                                cosx_blkptr_md, cosx_blkidx_md, cosx_blkcen_md, cosx_blkrad_md, &
                                cosx_QfullT_md, cosx_blkstart_md, cosx_nblk_md, &
                                per_atom_period_scale=cosx_peratom_lomd)
   else
      call cosx_build_one_grid(nConts, nrad_md, nsph_md, cosx_coor_md, cosx_w_md, cosx_npts_md, cosx_Qfull_md, &
                                cosx_blkptr_md, cosx_blkidx_md, cosx_blkcen_md, cosx_blkrad_md, &
                                cosx_QfullT_md, cosx_blkstart_md, cosx_nblk_md, intacc_eps=4.020d0, &
                                per_atom_period_scale=cosx_peratom_lomd)
   endif
   if (allocated(cosx_coor_lo)) deallocate(cosx_coor_lo,cosx_w_lo,cosx_Qfull_lo)
   if (allocated(cosx_QfullT_lo)) deallocate(cosx_QfullT_lo)
   if (allocated(cosx_blkptr_lo)) deallocate(cosx_blkptr_lo,cosx_blkidx_lo)
   if (allocated(cosx_blkcen_lo)) deallocate(cosx_blkcen_lo,cosx_blkrad_lo)
   if (allocated(cosx_blkstart_lo)) deallocate(cosx_blkstart_lo)
   if (nrad_lo_explicit) then
      call cosx_build_one_grid(nConts, nrad_lo, nsph_lo, cosx_coor_lo, cosx_w_lo, cosx_npts_lo, cosx_Qfull_lo, &
                                cosx_blkptr_lo, cosx_blkidx_lo, cosx_blkcen_lo, cosx_blkrad_lo, &
                                cosx_QfullT_lo, cosx_blkstart_lo, cosx_nblk_lo, &
                                per_atom_period_scale=cosx_peratom_lomd)
   else
      call cosx_build_one_grid(nConts, nrad_lo, nsph_lo, cosx_coor_lo, cosx_w_lo, cosx_npts_lo, cosx_Qfull_lo, &
                                cosx_blkptr_lo, cosx_blkidx_lo, cosx_blkcen_lo, cosx_blkrad_lo, &
                                cosx_QfullT_lo, cosx_blkstart_lo, cosx_nblk_lo, intacc_eps=3.816d0, &
                                per_atom_period_scale=cosx_peratom_lomd)
   endif

   print '(A,I0,A,I0,A,I0)', "  COSX grid points: hi=",cosx_npts_hi," md=",cosx_npts_md," lo=",cosx_npts_lo

   if (allocated(cosx_Xblock_hi)) deallocate(cosx_Xblock_hi)
   allocate(cosx_Xblock_hi(cosx_nblk_hi))
   cosx_Xblock_hi = 0.0d0
   if (allocated(cosx_Xblock_md)) deallocate(cosx_Xblock_md)
   allocate(cosx_Xblock_md(cosx_nblk_md))
   cosx_Xblock_md = 0.0d0
   if (allocated(cosx_Xblock_lo)) deallocate(cosx_Xblock_lo)
   allocate(cosx_Xblock_lo(cosx_nblk_lo))
   cosx_Xblock_lo = 0.0d0

   cosx_call_count = 0
   cosx_promote_md = .false.
   cosx_force_hi_grid = .false.
   cosx_incr_valid = .false.
   cosx_stage_prev = 0
   cosx_stage_settle = 0
   cosx_calls_since_full = 0
   cosx_incr_valid_lr = .false.
   cosx_stage_prev_lr = 0
   cosx_stage_settle_lr = 0
   cosx_incr_valid_sr = .false.
   cosx_stage_prev_sr = 0
   cosx_stage_settle_sr = 0
   cosx_G_hi_valid = .false.
   cosx_G_hi_valid_lr = .false.
   cosx_ready = .true.
endif

cosx_call_count = cosx_call_count + 1

if (cosx_force_hi_grid) then
   coor_p => cosx_coor_hi;  w_p => cosx_w_hi;  npts = cosx_npts_hi;  Qfull_p => cosx_Qfull_hi
   cosx_blkptr_p => cosx_blkptr_hi;  cosx_blkidx_p => cosx_blkidx_hi
   blkstart_p => cosx_blkstart_hi;  nblocks = cosx_nblk_hi
   blkcen_p => cosx_blkcen_hi;  blkrad_p => cosx_blkrad_hi;  QfullT_p => cosx_QfullT_hi
   stage_now = 3
else if (cosx_promote_md) then
   coor_p => cosx_coor_md;  w_p => cosx_w_md;  npts = cosx_npts_md;  Qfull_p => cosx_Qfull_md
   cosx_blkptr_p => cosx_blkptr_md;  cosx_blkidx_p => cosx_blkidx_md
   blkstart_p => cosx_blkstart_md;  nblocks = cosx_nblk_md
   blkcen_p => cosx_blkcen_md;  blkrad_p => cosx_blkrad_md;  QfullT_p => cosx_QfullT_md
   stage_now = 2
else
   coor_p => cosx_coor_lo;  w_p => cosx_w_lo;  npts = cosx_npts_lo;  Qfull_p => cosx_Qfull_lo
   cosx_blkptr_p => cosx_blkptr_lo;  cosx_blkidx_p => cosx_blkidx_lo
   blkstart_p => cosx_blkstart_lo;  nblocks = cosx_nblk_lo
   blkcen_p => cosx_blkcen_lo;  blkrad_p => cosx_blkrad_lo;  QfullT_p => cosx_QfullT_lo
   stage_now = 1
endif

stage_changed = (stage_now .ne. cosx_stage_prev)
do_incremental = cosx_incr_valid .and. (.not. stage_changed) .and. (cosx_stage_settle .eq. 0) &
                  .and. (.not. cosx_no_incremental)
kscreen_use = merge(cosx_kscreen_incr, cosx_kscreen, do_incremental)
cosx_G_hi_valid = .false.
is_full_hi_build = need_force .and. (.not. do_incremental) .and. (stage_now .eq. 3)
if (is_full_hi_build) then
   if (allocated(cosx_G_hi_a)) then
      if (size(cosx_G_hi_a,1) .ne. npts .or. size(cosx_G_hi_a,2) .ne. nConts) then
         deallocate(cosx_G_hi_a,cosx_G_hi_b)
      endif
   endif
   if (.not. allocated(cosx_G_hi_a)) allocate(cosx_G_hi_a(npts,nConts), cosx_G_hi_b(npts,nConts))
endif
if (do_incremental) then
   allocate(Da_eff(nConts,nConts), Db_eff(nConts,nConts))
   Da_eff = Da - cosx_D_incr_a
   Db_eff = Db - cosx_D_incr_b
else
   allocate(Da_eff(nConts,nConts), Db_eff(nConts,nConts))
   Da_eff = Da
   Db_eff = Db
endif

Ka = 0.0d0
Kb = 0.0d0

maxd = 0
do si = 0,nBases-1
   maxd = max(maxd, cosx_shell_dim(si))
enddo

select case (stage_now)
case (1); cosx_Xblock_p => cosx_Xblock_lo
case (2); cosx_Xblock_p => cosx_Xblock_md
case default; cosx_Xblock_p => cosx_Xblock_hi
end select
if (do_incremental) then
   global_D_max = max(maxval(abs(Da_eff)), maxval(abs(Db_eff)))
   allocate(D_pair_max(0:nBases-1,0:nBases-1))
   D_pair_max = 0.0d0
   allocate(ao2shell(nConts))
   do si = 0,nBases-1
      do ii = cosx_ao_offset(si)+1, cosx_ao_offset(si)+cosx_shell_dim(si)
         ao2shell(ii) = si
      enddo
   enddo
   !$omp parallel do default(shared) private(sj,nu,mu,si,dval) schedule(dynamic)
   do sj = 0,nBases-1
      do nu = cosx_ao_offset(sj)+1, cosx_ao_offset(sj)+cosx_shell_dim(sj)
         do mu = 1,nConts
            si = ao2shell(mu)
            dval = max(abs(Da_eff(mu,nu)), abs(Db_eff(mu,nu)))
            if (dval .gt. D_pair_max(si,sj)) D_pair_max(si,sj) = dval
         enddo
      enddo
   enddo
   !$omp end parallel do
   deallocate(ao2shell)
   allocate(D_weighted_max(0:nBases-1))
   !$omp parallel do default(shared) private(si) schedule(dynamic)
   do si = 0,nBases-1
      D_weighted_max(si) = maxval(cosx_esp_bound(si,:) * D_pair_max(si,:))
   enddo
   !$omp end parallel do
else
   global_D_max = 0.0d0
endif

maxnloc = 0
do ib = 0,nblocks-1
   nloc = 0
   do kk = cosx_blkptr_p(ib)+1, cosx_blkptr_p(ib+1)
      nloc = nloc + cosx_shell_dim(cosx_blkidx_p(kk))
   enddo
   maxnloc = max(maxnloc, nloc)
enddo
xmax_this_call = 0.0d0
skip_total = 0
n_geo_total = 0; n_schwarz_total = 0; n_calls_total = 0
skip_pts_total = 0

!$omp parallel default(shared) &
!$omp&   private(ib,bstart,nb,g,si,sj,di,dj,i,j,ao_i,ao_j,shls1e,off_grids,nij,wij) &
!$omp&   private(coor_batch,w_batch,Xb,val1_dummy,Fa,Fb,Ga,Gb,env_g,buf,int_cache,Ka_local,Kb_local,Qb) &
!$omp&   private(F_gmax,x_block_max,k_bound,pair_bound,sj_idx,xmax_local,skip_local,local_D_max,kk,skip_pts_local) &
!$omp&   private(n_geo_local,n_schwarz_local,n_calls_local) &
!$omp&   private(esp_eff,r_far,seg_ab,seg_len2,seg_t,seg_cl) &
!$omp&   private(loc_ao,nloc,cofs,iao,Xc,Kc) &
!$omp&   private(run_a,run_o,run_len,nrun,ir,a0,dsh,prev_end,bet) &
!$omp&   private(tau_hit,trun_a,trun_len,ntrun,it) &
!$omp&   private(dbg_t0,dbg_t1,dbg_tick_local)
dbg_tick_local = 0_8
allocate(buf(BATCH*maxd*maxd))
allocate(int_cache(INT_CACHE_SIZE))
allocate(loc_ao(maxnloc))
allocate(run_a(nBases), run_o(nBases), run_len(nBases))
allocate(tau_hit(0:nBases-1), trun_a(nBases), trun_len(nBases))
allocate(Ka_local(nConts,nConts), Kb_local(nConts,nConts))
Ka_local = 0.0d0
Kb_local = 0.0d0
xmax_local = 0.0d0
skip_local = 0
n_geo_local = 0; n_schwarz_local = 0; n_calls_local = 0
skip_pts_local = 0
allocate(env_g(size(env) + 3*BATCH))
env_g(1:size(env)) = env(1:size(env))
off_grids = size(env)

!$omp do schedule(static,1)
do ib = 0,nblocks-1
   bstart = blkstart_p(ib) + 1
   nb = blkstart_p(ib+1) - blkstart_p(ib)

   if (do_incremental) then
      local_D_max = 0.0d0
      do kk = cosx_blkptr_p(ib)+1, cosx_blkptr_p(ib+1)
         si = cosx_blkidx_p(kk)
         do sj_idx = cosx_ext_ptr_sym(si)+1, cosx_ext_ptr_sym(si+1)
            sj = cosx_ext_idx_sym(sj_idx)
            if (D_pair_max(si,sj) .le. 0.0d0) cycle
            seg_ab = cosx_shell_cen(:,sj) - cosx_shell_cen(:,si)
            seg_len2 = sum(seg_ab*seg_ab)
            if (seg_len2 .gt. 0.0d0) then
               seg_t = sum((blkcen_p(:,ib) - cosx_shell_cen(:,si))*seg_ab)/seg_len2
               seg_t = max(0.0d0, min(1.0d0, seg_t))
               seg_cl = cosx_shell_cen(:,si) + seg_t*seg_ab
            else
               seg_cl = cosx_shell_cen(:,si)
            endif
            r_far = sqrt(sum((blkcen_p(:,ib) - seg_cl)**2)) - blkrad_p(ib)
            if (r_far .gt. 0.0d0) then
               esp_eff = min(cosx_esp_bound(si,sj), cosx_esp_mono(si,sj)/r_far)
            else
               esp_eff = cosx_esp_bound(si,sj)
            endif
            local_D_max = max(local_D_max, esp_eff*D_pair_max(si,sj))
         enddo
      enddo
      local_D_max = max(local_D_max, cosx_ext_threshold*global_D_max)
      block_bound = (cosx_Xblock_p(ib+1)**2) * dble(maxd) * local_D_max
      if (block_bound .lt. kscreen_use) then
         skip_local = skip_local + 1
         skip_pts_local = skip_pts_local + nb
         cycle
      endif
   endif

   allocate(coor_batch(3,nb), w_batch(nb))
   coor_batch = coor_p(:,bstart:bstart+nb-1)
   w_batch = w_p(bstart:bstart+nb-1)

   nloc = 0
   nrun = 0
   prev_end = -1
   do kk = cosx_blkptr_p(ib)+1, cosx_blkptr_p(ib+1)
      si = cosx_blkidx_p(kk)
      a0 = cosx_ao_offset(si) + 1
      dsh = cosx_shell_dim(si)
      if (nrun .gt. 0 .and. a0 .eq. prev_end+1) then
         run_len(nrun) = run_len(nrun) + dsh
      else
         nrun = nrun + 1
         run_a(nrun) = a0
         run_o(nrun) = nloc
         run_len(nrun) = dsh
      endif
      prev_end = a0 + dsh - 1
      do iao = 1, dsh
         nloc = nloc + 1
         loc_ao(nloc) = cosx_ao_offset(si) + iao
      enddo
   enddo

   if (nloc .eq. 0) then
      deallocate(coor_batch,w_batch)
      skip_local = skip_local + 1
      skip_pts_local = skip_pts_local + nb
      cycle
   endif

   allocate(Xc(nloc,nb), val1_dummy(nloc,3,nb))
   Xc = 0.0d0
   cofs = 0
   do kk = cosx_blkptr_p(ib)+1, cosx_blkptr_p(ib+1)
      si = cosx_blkidx_p(kk)
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
                                coor_batch, nb, cofs, nloc, Xc, val1_dummy)
      cofs = cofs + cosx_shell_dim(si)
   enddo
   deallocate(val1_dummy)
   do g = 1,nb
      Xc(:,g) = Xc(:,g) * sqrt(abs(w_batch(g)))
   enddo

   allocate(Fa(nb,nConts), Fb(nb,nConts))
   do ir = 1,nrun
      bet = merge(1.0d0, 0.0d0, ir .gt. 1)
      call cosx_dgemm('T','N',nb,nConts,run_len(ir),1.0d0, &
                 Xc(run_o(ir)+1,1),nloc, Da_eff(run_a(ir),1),nConts, bet,Fa,nb)
      call cosx_dgemm('T','N',nb,nConts,run_len(ir),1.0d0, &
                 Xc(run_o(ir)+1,1),nloc, Db_eff(run_a(ir),1),nConts, bet,Fb,nb)
   enddo

   x_block_max = maxval(abs(Xc))
   xmax_local = max(xmax_local, x_block_max)
   if (.not. do_incremental) cosx_Xblock_p(ib+1) = x_block_max
   allocate(F_gmax(0:nBases-1))
   do si = 0,nBases-1
      F_gmax(si) = max(maxval(abs(Fa(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)))), &
                        maxval(abs(Fb(:,cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)))))
   enddo

   allocate(Ga(nb,nConts), Gb(nb,nConts))
   Ga = 0.0d0
   Gb = 0.0d0
   tau_hit = .false.

   do g = 1,nb
      env_g(off_grids+3*(g-1)+1) = coor_batch(1,g)
      env_g(off_grids+3*(g-1)+2) = coor_batch(2,g)
      env_g(off_grids+3*(g-1)+3) = coor_batch(3,g)
   enddo
   env_g(12) = dble(nb)
   env_g(13) = dble(off_grids)

   do si = 0,nBases-1
      di = cosx_shell_dim(si)
      do sj_idx = cosx_ext_ptr(si)+1, cosx_ext_ptr(si+1)
         sj = cosx_ext_idx(sj_idx)
         n_geo_local = n_geo_local + 1
         if (schwarz_bound(si,sj) .lt. SCHWARZ_CUTOFF) cycle
         n_schwarz_local = n_schwarz_local + 1
         seg_ab = cosx_shell_cen(:,sj) - cosx_shell_cen(:,si)
         seg_len2 = sum(seg_ab*seg_ab)
         if (seg_len2 .gt. 0.0d0) then
            seg_t = sum((blkcen_p(:,ib) - cosx_shell_cen(:,si))*seg_ab)/seg_len2
            seg_t = max(0.0d0, min(1.0d0, seg_t))
            seg_cl = cosx_shell_cen(:,si) + seg_t*seg_ab
         else
            seg_cl = cosx_shell_cen(:,si)
         endif
         r_far = sqrt(sum((blkcen_p(:,ib) - seg_cl)**2)) - blkrad_p(ib)
         if (r_far .gt. 0.0d0) then
            esp_eff = min(cosx_esp_bound(si,sj), cosx_esp_mono(si,sj)/r_far)
         else
            esp_eff = cosx_esp_bound(si,sj)
         endif
         if (do_incremental) then
            pair_bound = x_block_max*esp_eff*D_pair_max(si,sj)
            if (pair_bound .lt. kscreen_use) cycle
         endif
         k_bound = x_block_max*esp_eff*max(F_gmax(si),F_gmax(sj))
         if (k_bound .lt. kscreen_use) cycle
         n_calls_local = n_calls_local + 1
         tau_hit(si) = .true.
         tau_hit(sj) = .true.
         dj = cosx_shell_dim(sj)
         shls1e(1) = si
         shls1e(2) = sj
         shls1e(3) = 0
         shls1e(4) = nb
         call system_clock(dbg_t0)
         call grids1e_engine_cached(buf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         call system_clock(dbg_t1)
         dbg_tick_local = dbg_tick_local + (dbg_t1-dbg_t0)
         do j = 1,dj
            ao_j = cosx_ao_offset(sj) + j
            do i = 1,di
               ao_i = cosx_ao_offset(si) + i
               nij = NorVEC(ao_i)*NorVEC(ao_j)
               if (nij .eq. 0.0d0) cycle
               do g = 1,nb
                  wij = sign(1.0d0,w_batch(g)) * nij * buf(g + nb*((i-1)+di*(j-1)))
                  Ga(g,ao_i) = Ga(g,ao_i) + wij*Fa(g,ao_j)
                  Gb(g,ao_i) = Gb(g,ao_i) + wij*Fb(g,ao_j)
                  if (sj .ne. si) then
                     Ga(g,ao_j) = Ga(g,ao_j) + wij*Fa(g,ao_i)
                     Gb(g,ao_j) = Gb(g,ao_j) + wij*Fb(g,ao_i)
                  endif
               enddo
            enddo
         enddo
      enddo
   enddo

   if (is_full_hi_build) then
      cosx_G_hi_a(bstart:bstart+nb-1,:) = Ga
      cosx_G_hi_b(bstart:bstart+nb-1,:) = Gb
   endif

   ntrun = 0
   do it = 0,nBases-1
      if (.not. tau_hit(it)) cycle
      if (ntrun .gt. 0 .and. cosx_ao_offset(it)+1 .eq. trun_a(ntrun)+trun_len(ntrun)) then
         trun_len(ntrun) = trun_len(ntrun) + cosx_shell_dim(it)
      else
         ntrun = ntrun + 1
         trun_a(ntrun) = cosx_ao_offset(it) + 1
         trun_len(ntrun) = cosx_shell_dim(it)
      endif
   enddo

   do it = 1,ntrun
   enddo

   if (cosx_overlap_fit) then
      allocate(Qb(nConts,nb))
      do ir = 1,nrun
         bet = merge(1.0d0, 0.0d0, ir .gt. 1)
         call cosx_dgemm('N','N',nConts,nb,run_len(ir),1.0d0, &
                    QfullT_p(:,run_a(ir):run_a(ir)+run_len(ir)-1),nConts, &
                    Xc(run_o(ir)+1,1),nloc, bet,Qb,nConts)
      enddo
      do it = 1,ntrun
         call cosx_dgemm('N','N',nConts,trun_len(it),nb,1.0d0, &
                    Qb,nConts, Ga(1,trun_a(it)),nb, 1.0d0,Ka_local(1,trun_a(it)),nConts)
         call cosx_dgemm('N','N',nConts,trun_len(it),nb,1.0d0, &
                    Qb,nConts, Gb(1,trun_a(it)),nb, 1.0d0,Kb_local(1,trun_a(it)),nConts)
      enddo
      deallocate(Qb)
   else
      allocate(Kc(nloc,nConts))
      do it = 1,ntrun
         call cosx_dgemm('N','N',nloc,trun_len(it),nb,1.0d0, &
                    Xc,nloc, Ga(1,trun_a(it)),nb, 0.0d0,Kc(1,trun_a(it)),nloc)
      enddo
      do it = 1,ntrun
         do iao = 1,nloc
            Ka_local(loc_ao(iao),trun_a(it):trun_a(it)+trun_len(it)-1) = &
               Ka_local(loc_ao(iao),trun_a(it):trun_a(it)+trun_len(it)-1) &
               + Kc(iao,trun_a(it):trun_a(it)+trun_len(it)-1)
         enddo
      enddo
      do it = 1,ntrun
         call cosx_dgemm('N','N',nloc,trun_len(it),nb,1.0d0, &
                    Xc,nloc, Gb(1,trun_a(it)),nb, 0.0d0,Kc(1,trun_a(it)),nloc)
      enddo
      do it = 1,ntrun
         do iao = 1,nloc
            Kb_local(loc_ao(iao),trun_a(it):trun_a(it)+trun_len(it)-1) = &
               Kb_local(loc_ao(iao),trun_a(it):trun_a(it)+trun_len(it)-1) &
               + Kc(iao,trun_a(it):trun_a(it)+trun_len(it)-1)
         enddo
      enddo
      deallocate(Kc)
   endif

   deallocate(coor_batch,w_batch,Xc,Fa,Fb,Ga,Gb,F_gmax)
enddo
!$omp end do

block
   integer :: merge_tid, merge_nthreads
   merge_nthreads = omp_get_num_threads()
   do merge_tid = 0, merge_nthreads-1
      !$omp barrier
      if (omp_get_thread_num() .eq. merge_tid) then
      Ka = Ka + Ka_local
      Kb = Kb + Kb_local
      xmax_this_call = max(xmax_this_call, xmax_local)
      skip_total = skip_total + skip_local
      n_geo_total = n_geo_total + n_geo_local
      n_schwarz_total = n_schwarz_total + n_schwarz_local
      n_calls_total = n_calls_total + n_calls_local
      skip_pts_total = skip_pts_total + skip_pts_local
      dbg_full_ticks = dbg_full_ticks + dbg_tick_local
      dbg_full_calls = dbg_full_calls + n_calls_local
      endif
   enddo
end block
deallocate(buf,int_cache,Ka_local,Kb_local,env_g)
deallocate(loc_ao,run_a,run_o,run_len,tau_hit,trun_a,trun_len)
!$omp end parallel

if (is_full_hi_build) cosx_G_hi_valid = .true.

if (engine_verbose .ge. 2) then
   print '(A,I1,A,ES9.2,A,I0,A,I0,A,I0,A,I0,A)', "COSX: stage=",stage_now, &
         " global_D_max=",global_D_max," skipped ",skip_total,"/",nblocks," blocks (", &
         skip_pts_total,"/",npts," grid points)"
   print '(A,I0,A,I0,A,I0,A)', "COSX:   geo-pairs entered=",n_geo_total, &
         " schwarz-survived=",n_schwarz_total," grids1e_engine calls=",n_calls_total
   call flush(6)
endif

Ka = 0.5d0*(Ka + transpose(Ka))
Kb = 0.5d0*(Kb + transpose(Kb))

if (do_incremental) then
   Ka = Ka + cosx_Ka_prev
   Kb = Kb + cosx_Kb_prev
endif

if (allocated(cosx_Ka_prev)) deallocate(cosx_Ka_prev,cosx_Kb_prev)
allocate(cosx_Ka_prev(nConts,nConts), cosx_Kb_prev(nConts,nConts))
cosx_Ka_prev = Ka
cosx_Kb_prev = Kb
if (allocated(cosx_D_incr_a)) deallocate(cosx_D_incr_a,cosx_D_incr_b)
allocate(cosx_D_incr_a(nConts,nConts), cosx_D_incr_b(nConts,nConts))
cosx_D_incr_a = Da
cosx_D_incr_b = Db
cosx_incr_valid = .true.
if (stage_changed) then
   cosx_stage_settle = 1
else if (cosx_stage_settle .gt. 0) then
   cosx_stage_settle = cosx_stage_settle - 1
endif
if (do_incremental) then
   cosx_calls_since_full = cosx_calls_since_full + 1
else
   cosx_calls_since_full = 0
endif
cosx_stage_prev = stage_now
deallocate(Da_eff, Db_eff)
if (allocated(D_pair_max)) deallocate(D_pair_max)
if (allocated(D_weighted_max)) deallocate(D_weighted_max)
call cosx_malloc_trim()

end subroutine integrals_build_exchange_cosx
end submodule exchange_cosx_k_impl
