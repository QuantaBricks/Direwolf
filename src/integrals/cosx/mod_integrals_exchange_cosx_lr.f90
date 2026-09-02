! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX long-range K_LR (erf(omega*r)/r) builder - split out of mod_integrals_exchange_cosx.f90

submodule (mod_integrals:exchange_cosx_impl) exchange_cosx_lr_impl
contains
module subroutine integrals_build_exchange_cosx_lr(nConts, Da, Db, Ka, Kb)
use mod_exchange, only: RS_omega
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
use MOL_info, only: engine_verbose
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)

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
integer :: npts
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
real(8) :: xmax_this_call, xmax_local
integer :: skip_local, skip_total
integer(8) :: n_geo_local, n_schwarz_local, n_calls_local
integer(8) :: n_geo_total, n_schwarz_total, n_calls_total
integer :: skip_pts_local, skip_pts_total
real(8) :: global_D_max, block_bound
real(8) :: kscreen_use_lr
integer(8) :: dbg_t0, dbg_t1, dbg_tick_local
real(8),pointer :: cosx_Xblock_p(:)
integer :: ni_b, nj_b, offi_b, offj_b, pi_b, pj_b
real(8) :: r2_b, ei_b, ci_b, ej_b, cj_b, p_b
real(8),pointer :: blkcen_p(:,:), blkrad_p(:)
real(8) :: esp_eff, r_far, seg_ab(3), seg_len2, seg_t, seg_cl(3)
integer,allocatable :: loc_ao(:)
integer,allocatable :: run_a(:), run_o(:), run_len(:)
integer :: nloc, cofs, maxnloc, iao, nrun, ir, a0, dsh, prev_end
real(8) :: bet
real(8),allocatable :: Xc(:,:), Kc(:,:)
logical,allocatable :: tau_hit(:)
integer,allocatable :: trun_a(:), trun_len(:)
integer :: ntrun, it
real(8),pointer,contiguous :: QfullT_p(:,:)

if (.not. cosx_ready .or. cosx_nBases_cached .ne. nBases) then
   print *, "COSX-LR: called before the full-range COSX grid/screening ", &
            "cache was built - scf_fock.f90 must call exchange_build_cosx ", &
            "before exchange_build_cosx_lr in the same Fock build"
   call flush(6)
   stop 1
endif

if ((.not. cosx_esp_bound_lr_ready) .or. (cosx_esp_bound_lr_omega .ne. RS_omega)) then
   if (allocated(cosx_esp_bound_lr)) deallocate(cosx_esp_bound_lr)
   allocate(cosx_esp_bound_lr(0:nBases-1,0:nBases-1))
   cosx_esp_bound_lr = 0.0d0
   !$omp parallel do default(shared) &
   !$omp&   private(si,sj,ni_b,nj_b,offi_b,offj_b,r2_b,pi_b,pj_b,ei_b,ci_b,ej_b,cj_b,p_b) schedule(dynamic)
   do si = 0,nBases-1
      ni_b = bas(3,si+1)
      offi_b = bas(6,si+1)
      do sj = si,nBases-1
         nj_b = bas(3,sj+1)
         offj_b = bas(6,sj+1)
         r2_b = sum((cosx_shell_cen(:,si)-cosx_shell_cen(:,sj))**2)
         do pi_b = 1,ni_b
            ei_b = env(offi_b+pi_b)
            ci_b = env(bas(7,si+1)+pi_b)
            do pj_b = 1,nj_b
               ej_b = env(offj_b+pj_b)
               cj_b = env(bas(7,sj+1)+pj_b)
               p_b = ei_b + ej_b
               cosx_esp_bound_lr(si,sj) = cosx_esp_bound_lr(si,sj) + &
                  abs(ci_b*cj_b)*exp(-r2_b*ei_b*ej_b/p_b) * &
                  2.0d0*acos(-1.0d0)*RS_omega/(p_b*sqrt(p_b+RS_omega*RS_omega))
            enddo
         enddo
         cosx_esp_bound_lr(sj,si) = cosx_esp_bound_lr(si,sj)
      enddo
   enddo
   !$omp end parallel do
   cosx_esp_bound_lr_max = maxval(cosx_esp_bound_lr)
   call cosx_build_ext_lr(nBases)
   cosx_kscreen_lr = cosx_kscreen
   cosx_kscreen_incr_lr = cosx_kscreen_incr * (cosx_esp_bound_lr_max/cosx_esp_bound_max)
   cosx_esp_bound_lr_ready = .true.
   cosx_esp_bound_lr_omega = RS_omega
endif

if (cosx_force_hi_grid) then
   coor_p => cosx_coor_hi;  w_p => cosx_w_hi;  npts = cosx_npts_hi;  Qfull_p => cosx_Qfull_hi
   cosx_blkptr_p => cosx_blkptr_hi;  cosx_blkidx_p => cosx_blkidx_hi
   blkstart_p => cosx_blkstart_hi;  nblocks = cosx_nblk_hi
   blkcen_p => cosx_blkcen_hi;  blkrad_p => cosx_blkrad_hi
   QfullT_p => cosx_QfullT_hi
   stage_now = 3
else if (cosx_promote_md) then
   coor_p => cosx_coor_md;  w_p => cosx_w_md;  npts = cosx_npts_md;  Qfull_p => cosx_Qfull_md
   cosx_blkptr_p => cosx_blkptr_md;  cosx_blkidx_p => cosx_blkidx_md
   blkstart_p => cosx_blkstart_md;  nblocks = cosx_nblk_md
   blkcen_p => cosx_blkcen_md;  blkrad_p => cosx_blkrad_md
   QfullT_p => cosx_QfullT_md
   stage_now = 2
else
   coor_p => cosx_coor_lo;  w_p => cosx_w_lo;  npts = cosx_npts_lo;  Qfull_p => cosx_Qfull_lo
   cosx_blkptr_p => cosx_blkptr_lo;  cosx_blkidx_p => cosx_blkidx_lo
   blkstart_p => cosx_blkstart_lo;  nblocks = cosx_nblk_lo
   blkcen_p => cosx_blkcen_lo;  blkrad_p => cosx_blkrad_lo
   QfullT_p => cosx_QfullT_lo
   stage_now = 1
endif

stage_changed = (stage_now .ne. cosx_stage_prev_lr)
do_incremental = cosx_incr_valid_lr .and. (.not. stage_changed) .and. (cosx_stage_settle_lr .eq. 0) &
                  .and. (.not. cosx_no_incremental) .and. (.not. cosx_no_incremental_lr)
kscreen_use_lr = merge(cosx_kscreen_incr_lr, cosx_kscreen_lr, do_incremental)
if (do_incremental) then
   allocate(Da_eff(nConts,nConts), Db_eff(nConts,nConts))
   Da_eff = Da - cosx_D_incr_a_lr
   Db_eff = Db - cosx_D_incr_b_lr
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
      D_weighted_max(si) = maxval(cosx_esp_bound_lr(si,:) * D_pair_max(si,:))
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
!$omp&   private(loc_ao,nloc,cofs,iao,Xc,Kc,run_a,run_o,run_len,nrun,ir,a0,dsh,prev_end,bet) &
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
env_g(9) = RS_omega

!$omp do schedule(static)
do ib = 0,nblocks-1
   bstart = blkstart_p(ib) + 1
   nb = blkstart_p(ib+1) - blkstart_p(ib)

   if (do_incremental) then
      local_D_max = 0.0d0
      do kk = cosx_blkptr_p(ib)+1, cosx_blkptr_p(ib+1)
         si = cosx_blkidx_p(kk)
         do sj_idx = cosx_ext_ptr_lr_sym(si)+1, cosx_ext_ptr_lr_sym(si+1)
            sj = cosx_ext_idx_lr_sym(sj_idx)
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
               esp_eff = min(cosx_esp_bound_lr(si,sj), cosx_esp_mono(si,sj)/r_far)
            else
               esp_eff = cosx_esp_bound_lr(si,sj)
            endif
            local_D_max = max(local_D_max, esp_eff*D_pair_max(si,sj))
         enddo
      enddo
      local_D_max = max(local_D_max, cosx_ext_threshold*global_D_max)
      block_bound = (cosx_Xblock_p(ib+1)**2) * dble(maxd) * local_D_max
      if (block_bound .lt. kscreen_use_lr) then
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
      do sj_idx = cosx_ext_ptr_lr(si)+1, cosx_ext_ptr_lr(si+1)
         sj = cosx_ext_idx_lr(sj_idx)
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
            esp_eff = min(cosx_esp_bound_lr(si,sj), cosx_esp_mono(si,sj)/r_far)
         else
            esp_eff = cosx_esp_bound_lr(si,sj)
         endif
         if (do_incremental) then
            pair_bound = x_block_max*esp_eff*D_pair_max(si,sj)
            if (pair_bound .lt. kscreen_use_lr) cycle
         endif
         k_bound = x_block_max*esp_eff*max(F_gmax(si),F_gmax(sj))
         if (k_bound .lt. kscreen_use_lr) cycle
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
      dbg_lr_ticks = dbg_lr_ticks + dbg_tick_local
      dbg_lr_calls = dbg_lr_calls + n_calls_local
      endif
   enddo
end block
deallocate(buf,int_cache,Ka_local,Kb_local,env_g)
deallocate(loc_ao,run_a,run_o,run_len,tau_hit,trun_a,trun_len)
!$omp end parallel

if (engine_verbose .ge. 2) then
   print '(A,I1,A,ES9.2,A,I0,A,I0,A,I0,A,I0,A)', "COSX-LR: stage=",stage_now, &
         " global_D_max=",global_D_max," skipped ",skip_total,"/",nblocks," blocks (", &
         skip_pts_total,"/",npts," grid points)"
   print '(A,I0,A,I0,A,I0,A)', "COSX-LR:   geo-pairs entered=",n_geo_total, &
         " schwarz-survived=",n_schwarz_total," grids1e_engine calls=",n_calls_total
   block
      integer(8) :: dbg_rate
      call system_clock(count_rate=dbg_rate)
      print '(A,ES12.4,A,I0,A,ES12.4,A,I0,A,F8.2,A)', &
         "DBG: full ticks=",dble(dbg_full_ticks)," calls=",dbg_full_calls, &
         "  LR ticks=",dble(dbg_lr_ticks)," calls=",dbg_lr_calls, &
         "  ratio(LR_ns_per_call/full_ns_per_call)=", &
         (dble(dbg_lr_ticks)/max(1_8,dbg_lr_calls))/max(1.0d-30,(dble(dbg_full_ticks)/max(1_8,dbg_full_calls)))
   end block
   call flush(6)
endif

Ka = 0.5d0*(Ka + transpose(Ka))
Kb = 0.5d0*(Kb + transpose(Kb))

if (do_incremental) then
   Ka = Ka + cosx_Ka_prev_lr
   Kb = Kb + cosx_Kb_prev_lr
endif

if (allocated(cosx_Ka_prev_lr)) deallocate(cosx_Ka_prev_lr,cosx_Kb_prev_lr)
allocate(cosx_Ka_prev_lr(nConts,nConts), cosx_Kb_prev_lr(nConts,nConts))
cosx_Ka_prev_lr = Ka
cosx_Kb_prev_lr = Kb
if (allocated(cosx_D_incr_a_lr)) deallocate(cosx_D_incr_a_lr,cosx_D_incr_b_lr)
allocate(cosx_D_incr_a_lr(nConts,nConts), cosx_D_incr_b_lr(nConts,nConts))
cosx_D_incr_a_lr = Da
cosx_D_incr_b_lr = Db
cosx_incr_valid_lr = .true.
if (stage_changed) then
   cosx_stage_settle_lr = 1
else if (cosx_stage_settle_lr .gt. 0) then
   cosx_stage_settle_lr = cosx_stage_settle_lr - 1
endif
cosx_stage_prev_lr = stage_now
deallocate(Da_eff, Db_eff)
if (allocated(D_pair_max)) deallocate(D_pair_max)
if (allocated(D_weighted_max)) deallocate(D_weighted_max)
call cosx_malloc_trim()

end subroutine integrals_build_exchange_cosx_lr
end submodule exchange_cosx_lr_impl
