! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Do DFT calculation

subroutine  DFT_calc(Exc,Fxc_a,Fxc_b,info,need_deriv)

use MOL_info
use GRID_INFO
use mod_xc, only: xc_eval, xc_uses_tau, xc_uses_vv10
use mod_vv10, only: vv10_evaluate, vv10_nlc_grid_build, vv10_nlc_basis_build, &
                     NLC_CHUNK, vv10_nlc_chunk_info, vv10_nlc_chunk_val, engine_vv10_nonself, &
                     vv10_active_now
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
INCLUDE 'parameter.h'
    integer,parameter :: BATCH = 64
    real(8),parameter :: SIG_CUTOFF = 1.0d-7
    integer    :: info
    logical    :: need_deriv
    integer    :: igrid,ip,nb,batch_start
    integer    :: ii,jj,n_sig,n_sig_wide,dd
    logical    :: mgga_active
    logical    :: vv10_active
    real(8)    :: Fxc_a_s_val,Fxc_b_s_val
    integer,allocatable :: sig_idx(:),narrow_idx(:)
    integer,allocatable :: wide_idx(:)
    integer    :: wide_ld
    integer    :: ib,nbatches,ib_cache,local_start
    logical,save :: sp_mode = .false., sp_env_checked = .false.
    character(len=8) :: sp_envval
    real(8)    :: Exc
    real(8)    :: Fxc_a(nconts,nconts),Fxc_b(nconts,nconts)
    real(8)    :: Temp_B(3),Temp_Ba(3),Temp_Bb(3)
    real(8),allocatable :: ex_b(:),ec_b(:)
    real(8),allocatable :: sigma_aa_b(:),sigma_ab_b(:),sigma_bb_b(:),sigma_t_b(:)

    real(8),allocatable :: ValBatch_own(:,:),Val1Batch_own(:,:,:)
    real(8),allocatable :: WeightBatch(:)
    real(8),allocatable :: Wa(:,:),Wb(:,:)
    real(8),allocatable :: rho_a_b(:),rho_b_b(:)
    real(8),allocatable :: drv_a_b(:,:),drv_b_b(:,:)
    real(8),allocatable :: d1rho_b(:),d1sig_b(:)
    real(8),allocatable :: vrhoa_b(:),vrhob_b(:),vsigmaaa_b(:),vsigmabb_b(:),vsigmaab_b(:)
    real(8),allocatable :: tau_a_b(:),tau_b_b(:),vtaua_b(:),vtaub_b(:)
    real(8),allocatable :: Val1d_tmp(:,:),Wa1_tmp(:,:),ValTauW(:,:)
    real(8),allocatable :: Et_b(:)
    real(8),allocatable :: Ta(:,:),Tb(:,:)
    real(8),allocatable :: ValW1a(:,:),ValW1b(:,:)
    real(8),allocatable :: Pa_s(:,:),Pb_s(:,:),Fxc_a_s(:,:),Fxc_b_s(:,:)
    real(8),allocatable :: ValBatch_sig(:,:)
    real(8),allocatable :: val1_sig(:,:,:)
    real(8),allocatable :: ValBatch_wide(:,:),Val1Batch_wide(:,:,:)
    real(8),allocatable :: maxval_bf(:)
    real(8),allocatable :: catPab(:,:),catWab(:,:)
    real(8),allocatable :: Fxc_a_local(:,:),Fxc_b_local(:,:)
    real(8) :: Exc_local

    Fxc_a = 0
    Fxc_b = 0
    Exc = 0
    mgga_active = xc_uses_tau()
    vv10_active = xc_uses_vv10()

    if (.not. sp_env_checked) then
       sp_envval = ""
       call get_environment_variable("ENGINE_XC_SHELLPAIR", sp_envval)
       sp_mode = (trim(sp_envval) .eq. "1")
       sp_env_checked = .true.
       if (sp_mode) print *,"XC build = shell-pair driven (ENGINE_XC_SHELLPAIR=1)"
    endif
    if (sp_mode .and. Multi .eq. 1 .and. .not. mgga_active .and. .not. vv10_active &
        .and. .not. need_deriv) then
       call DFT_calc_shellpair(Exc,Fxc_a,Fxc_b,info)
       return
    endif

    if (need_deriv .and. .not. allocated(TempD_all)) allocate(TempD_all(10,ngrids))

    if (vv10_active .and. .not. engine_vv10_nonself .and. vv10_active_now) call DFT_vv10_nlc_pass(Fxc_a, Fxc_b, Exc)

    !$omp parallel &
    !$omp&   private(ib,batch_start,nb,ib_cache,local_start,igrid,ip,ii,jj,n_sig,n_sig_wide,dd,sig_idx,narrow_idx) &
    !$omp&   private(Fxc_a_s_val,Fxc_b_s_val) &
    !$omp&   private(Temp_B,Temp_Ba,Temp_Bb) &
    !$omp&   private(ValBatch_own,Val1Batch_own,ValBatch_wide,Val1Batch_wide,maxval_bf,wide_idx) &
    !$omp&   private(WeightBatch,Wa,Wb,rho_a_b,rho_b_b,drv_a_b,drv_b_b) &
    !$omp&   private(d1rho_b,d1sig_b,Et_b,Ta,Tb,ValW1a,ValW1b) &
    !$omp&   private(Pa_s,Pb_s,Fxc_a_s,Fxc_b_s,ValBatch_sig,val1_sig) &
    !$omp&   private(catPab,catWab,ex_b,ec_b,sigma_aa_b,sigma_ab_b,sigma_bb_b,sigma_t_b) &
    !$omp&   private(vrhoa_b,vrhob_b,vsigmaaa_b,vsigmabb_b,vsigmaab_b) &
    !$omp&   private(tau_a_b,tau_b_b,vtaua_b,vtaub_b,Val1d_tmp,Wa1_tmp,ValTauW) &
    !$omp&   private(Fxc_a_local,Fxc_b_local,Exc_local)
    allocate(WeightBatch(BATCH))
    allocate(Wa(nconts,BATCH))
    allocate(rho_a_b(BATCH),drv_a_b(BATCH,3))
    allocate(d1rho_b(BATCH),d1sig_b(BATCH))
    allocate(Et_b(BATCH))
    allocate(Ta(nconts,BATCH))
    allocate(ValW1a(nconts,BATCH))
    allocate(Pa_s(nconts,nconts),Fxc_a_s(nconts,nconts))
    allocate(ValBatch_sig(nconts,BATCH))
    allocate(val1_sig(nconts,3,BATCH))
    allocate(ex_b(BATCH),ec_b(BATCH))
    allocate(sigma_t_b(BATCH))
    allocate(sigma_aa_b(BATCH),sigma_ab_b(BATCH),sigma_bb_b(BATCH))
    allocate(vrhoa_b(BATCH),vrhob_b(BATCH),vsigmaaa_b(BATCH),vsigmabb_b(BATCH),vsigmaab_b(BATCH))
    if (Multi .ne. 1) then
       allocate(Wb(nconts,BATCH))
       allocate(rho_b_b(BATCH),drv_b_b(BATCH,3))
       allocate(Tb(nconts,BATCH))
       allocate(ValW1b(nconts,BATCH))
       allocate(Pb_s(nconts,nconts),Fxc_b_s(nconts,nconts))
       allocate(catPab(2*nconts,nconts),catWab(2*nconts,BATCH))
    endif
    if (mgga_active) then
       allocate(tau_a_b(BATCH),vtaua_b(BATCH),vtaub_b(BATCH))
       allocate(Val1d_tmp(nconts,BATCH),Wa1_tmp(nconts,BATCH),ValTauW(nconts,BATCH))
       if (Multi .ne. 1) allocate(tau_b_b(BATCH))
    endif
    allocate(ValBatch_own(nconts,BATCH))
    allocate(Val1Batch_own(nconts,3,BATCH))
    wide_ld = max(max_batch_nsig, max_dft_batch_nsig)
    allocate(ValBatch_wide(wide_ld,BATCH),Val1Batch_wide(wide_ld,3,BATCH))
    allocate(maxval_bf(wide_ld))
    allocate(sig_idx(nconts),narrow_idx(wide_ld),wide_idx(wide_ld))
    allocate(Fxc_a_local(nconts,nconts),Fxc_b_local(nconts,nconts))
    Fxc_a_local = 0
    Fxc_b_local = 0
    Exc_local = 0

    nbatches = (ngrids + BATCH - 1) / BATCH
    !$omp do schedule(static,16)
    do ib = 0, nbatches-1
       batch_start = ib*BATCH + 1
       nb = min(BATCH, ngrids-batch_start+1)

       ib_cache = (batch_start-1)/CACHE_BATCH + 1
       n_sig_wide = batch_nsig(ib_cache)
       wide_idx(1:n_sig_wide) = batch_sig_idx(1:n_sig_wide,ib_cache)

       if (grid_cache_mode .and. batch_start+nb-1 .le. grid_cache_capacity) then
          local_start = batch_start - (ib_cache-1)*CACHE_BATCH
          ValBatch_wide(1:n_sig_wide,1:nb) = val_blocks(ib_cache)%val0(1:n_sig_wide,local_start:local_start+nb-1)
          Val1Batch_wide(1:n_sig_wide,1:3,1:nb) = val_blocks(ib_cache)%val1(1:n_sig_wide,1:3,local_start:local_start+nb-1)
       else
          if (dft_loose_built .and. xc_geo_level .lt. XC_NLEVEL) then
             n_sig_wide = dft_batch_shells_lvl(ib+1,xc_geo_level)%n_sig
             wide_idx(1:n_sig_wide) = dft_batch_shells_lvl(ib+1,xc_geo_level)%sig_idx(1:n_sig_wide)
          else
             n_sig_wide = dft_batch_shells(ib+1)%n_sig
             wide_idx(1:n_sig_wide) = dft_batch_shells(ib+1)%sig_idx(1:n_sig_wide)
          endif
          call GTOeval_batch_compact(ib+1,batch_start,nb,wide_ld,ValBatch_wide,Val1Batch_wide)
       endif

       maxval_bf(1:n_sig_wide) = 0.0d0
       do ip = 1,nb
          maxval_bf(1:n_sig_wide) = max(maxval_bf(1:n_sig_wide), abs(ValBatch_wide(1:n_sig_wide,ip)))
       enddo
       n_sig = 0
       do ii = 1,n_sig_wide
          if (maxval_bf(ii) .gt. xc_sig_now) then
             n_sig = n_sig + 1
             sig_idx(n_sig) = wide_idx(ii)
             narrow_idx(n_sig) = ii
          endif
       enddo
       ValBatch_sig(1:n_sig,1:nb) = ValBatch_wide(narrow_idx(1:n_sig),1:nb)
       val1_sig(1:n_sig,1:3,1:nb) = Val1Batch_wide(narrow_idx(1:n_sig),1:3,1:nb)
       do ip = 1,nb
          igrid = batch_start+ip-1
          WeightBatch(ip) = grids(igrid)%weight
          if (need_deriv .and. .not. associated(grids(igrid)%TempD)) grids(igrid)%TempD => TempD_all(:,igrid)
       enddo
       do jj = 1,n_sig
          do ii = 1,jj
             Pa_s(ii,jj) = Pa(sig_idx(ii),sig_idx(jj))
             if (ii .ne. jj) Pa_s(jj,ii) = Pa_s(ii,jj)
          enddo
       enddo
       if (Multi .ne. 1) then
          do jj = 1,n_sig
             do ii = 1,jj
                Pb_s(ii,jj) = Pb(sig_idx(ii),sig_idx(jj))
                if (ii .ne. jj) Pb_s(jj,ii) = Pb_s(ii,jj)
             enddo
          enddo
          catPab(1:n_sig,1:n_sig) = Pa_s(1:n_sig,1:n_sig)
          catPab(n_sig+1:2*n_sig,1:n_sig) = Pb_s(1:n_sig,1:n_sig)
          call dgemm('N','N',2*n_sig,nb,n_sig,1.0d0,catPab,2*nconts,ValBatch_sig,nconts,0.0d0,catWab,2*nconts)
          Wa(1:n_sig,1:nb) = catWab(1:n_sig,1:nb)
          Wb(1:n_sig,1:nb) = catWab(n_sig+1:2*n_sig,1:nb)
       else
          call dgemm('N','N',n_sig,nb,n_sig,1.0d0,Pa_s,nconts,ValBatch_sig,nconts,0.0d0,Wa,nconts)
       endif

       if (mgga_active) then
          tau_a_b(1:nb) = 0.0d0
          do dd = 1,3
             Val1d_tmp(1:n_sig,1:nb) = val1_sig(1:n_sig,dd,1:nb)
             call dgemm('N','N',n_sig,nb,n_sig,1.0d0,Pa_s,nconts,Val1d_tmp,nconts,0.0d0,Wa1_tmp,nconts)
             tau_a_b(1:nb) = tau_a_b(1:nb) + sum(Val1d_tmp(1:n_sig,1:nb)*Wa1_tmp(1:n_sig,1:nb),dim=1)
          enddo
          tau_a_b(1:nb) = 0.5d0*tau_a_b(1:nb)
          if (Multi .ne. 1) then
             tau_b_b(1:nb) = 0.0d0
             do dd = 1,3
                Val1d_tmp(1:n_sig,1:nb) = val1_sig(1:n_sig,dd,1:nb)
                call dgemm('N','N',n_sig,nb,n_sig,1.0d0,Pb_s,nconts,Val1d_tmp,nconts,0.0d0,Wa1_tmp,nconts)
                tau_b_b(1:nb) = tau_b_b(1:nb) + sum(Val1d_tmp(1:n_sig,1:nb)*Wa1_tmp(1:n_sig,1:nb),dim=1)
             enddo
             tau_b_b(1:nb) = 0.5d0*tau_b_b(1:nb)
          endif
       endif

       rho_a_b(1:nb) = sum(ValBatch_sig(1:n_sig,1:nb)*Wa(1:n_sig,1:nb),dim=1)
       drv_a_b(1:nb,1) = 2*sum(Wa(1:n_sig,1:nb)*val1_sig(1:n_sig,1,1:nb),dim=1)
       drv_a_b(1:nb,2) = 2*sum(Wa(1:n_sig,1:nb)*val1_sig(1:n_sig,2,1:nb),dim=1)
       drv_a_b(1:nb,3) = 2*sum(Wa(1:n_sig,1:nb)*val1_sig(1:n_sig,3,1:nb),dim=1)
       if (Multi .eq. 1) then
          sigma_t_b(1:nb) = (2*drv_a_b(1:nb,1))**2 + (2*drv_a_b(1:nb,2))**2 + (2*drv_a_b(1:nb,3))**2
          if (mgga_active) then
             call xc_eval(Multi,1,nb,rho_a_b(1:nb),rho_a_b(1:nb), &
                          sigma_aa_b(1:nb),sigma_bb_b(1:nb),sigma_ab_b(1:nb),sigma_t_b(1:nb), &
                          ex_b(1:nb),ec_b(1:nb),d1rho_b(1:nb),d1sig_b(1:nb), &
                          vrhoa_b(1:nb),vrhob_b(1:nb),vsigmaaa_b(1:nb),vsigmabb_b(1:nb),vsigmaab_b(1:nb), &
                          tau_a_b(1:nb),tau_a_b(1:nb),vtaua_b(1:nb),vtaub_b(1:nb))
          else
             call xc_eval(Multi,1,nb,rho_a_b(1:nb),rho_a_b(1:nb), &
                          sigma_aa_b(1:nb),sigma_bb_b(1:nb),sigma_ab_b(1:nb),sigma_t_b(1:nb), &
                          ex_b(1:nb),ec_b(1:nb),d1rho_b(1:nb),d1sig_b(1:nb), &
                          vrhoa_b(1:nb),vrhob_b(1:nb),vsigmaaa_b(1:nb),vsigmabb_b(1:nb),vsigmaab_b(1:nb))
          endif
          do ip = 1,nb
             igrid = batch_start+ip-1
             Et_b(ip) = ex_b(ip)+ec_b(ip)
             Temp_B = 4*d1sig_b(ip)*drv_a_b(ip,:)
             if (need_deriv) then
                TempD_all(1,igrid) = d1rho_b(ip)
                TempD_all(2:4,igrid) = Temp_B
                if (mgga_active) TempD_all(9,igrid) = vtaua_b(ip)
             endif
             Ta(1:n_sig,ip) = val1_sig(1:n_sig,1,ip)*Temp_B(1) + val1_sig(1:n_sig,2,ip)*Temp_B(2) &
                            + val1_sig(1:n_sig,3,ip)*Temp_B(3)
          enddo
       else
          rho_b_b(1:nb) = sum(ValBatch_sig(1:n_sig,1:nb)*Wb(1:n_sig,1:nb),dim=1)
          drv_b_b(1:nb,1) = 2*sum(Wb(1:n_sig,1:nb)*val1_sig(1:n_sig,1,1:nb),dim=1)
          drv_b_b(1:nb,2) = 2*sum(Wb(1:n_sig,1:nb)*val1_sig(1:n_sig,2,1:nb),dim=1)
          drv_b_b(1:nb,3) = 2*sum(Wb(1:n_sig,1:nb)*val1_sig(1:n_sig,3,1:nb),dim=1)
          do ip = 1,nb
             sigma_aa_b(ip) = sum(drv_a_b(ip,:)**2)
             sigma_bb_b(ip) = sum(drv_b_b(ip,:)**2)
             sigma_ab_b(ip) = sum(drv_a_b(ip,:)*drv_b_b(ip,:))
          enddo
          if (mgga_active) then
             call xc_eval(Multi,1,nb,rho_a_b(1:nb),rho_b_b(1:nb), &
                          sigma_aa_b(1:nb),sigma_bb_b(1:nb),sigma_ab_b(1:nb),sigma_t_b(1:nb), &
                          ex_b(1:nb),ec_b(1:nb),d1rho_b(1:nb),d1sig_b(1:nb), &
                          vrhoa_b(1:nb),vrhob_b(1:nb),vsigmaaa_b(1:nb),vsigmabb_b(1:nb),vsigmaab_b(1:nb), &
                          tau_a_b(1:nb),tau_b_b(1:nb),vtaua_b(1:nb),vtaub_b(1:nb))
          else
             call xc_eval(Multi,1,nb,rho_a_b(1:nb),rho_b_b(1:nb), &
                          sigma_aa_b(1:nb),sigma_bb_b(1:nb),sigma_ab_b(1:nb),sigma_t_b(1:nb), &
                          ex_b(1:nb),ec_b(1:nb),d1rho_b(1:nb),d1sig_b(1:nb), &
                          vrhoa_b(1:nb),vrhob_b(1:nb),vsigmaaa_b(1:nb),vsigmabb_b(1:nb),vsigmaab_b(1:nb))
          endif
          do ip = 1,nb
             igrid = batch_start+ip-1
             Et_b(ip) = ex_b(ip)+ec_b(ip)
             Temp_Ba = 2*vsigmaaa_b(ip)*drv_a_b(ip,:) + vsigmaab_b(ip)*drv_b_b(ip,:)
             Temp_Bb = 2*vsigmabb_b(ip)*drv_b_b(ip,:) + vsigmaab_b(ip)*drv_a_b(ip,:)
             if (need_deriv) then
                TempD_all(1,igrid) = vrhoa_b(ip)
                TempD_all(2:4,igrid) = Temp_Ba
                TempD_all(5,igrid) = vrhob_b(ip)
                TempD_all(6:8,igrid) = Temp_Bb
                if (mgga_active) then
                   TempD_all(9,igrid) = vtaua_b(ip)
                   TempD_all(10,igrid) = vtaub_b(ip)
                endif
             endif
             Ta(1:n_sig,ip) = val1_sig(1:n_sig,1,ip)*Temp_Ba(1) + val1_sig(1:n_sig,2,ip)*Temp_Ba(2) &
                            + val1_sig(1:n_sig,3,ip)*Temp_Ba(3)
             Tb(1:n_sig,ip) = val1_sig(1:n_sig,1,ip)*Temp_Bb(1) + val1_sig(1:n_sig,2,ip)*Temp_Bb(2) &
                            + val1_sig(1:n_sig,3,ip)*Temp_Bb(3)
          enddo
       endif
       do ip = 1,nb
          Exc_local = Exc_local + Et_b(ip)*WeightBatch(ip)
       enddo

       if (Multi .eq. 1) then
          do ip = 1,nb
             ValW1a(1:n_sig,ip) = 0.5d0*ValBatch_sig(1:n_sig,ip)*(d1rho_b(ip)*WeightBatch(ip)) &
                                + Ta(1:n_sig,ip)*WeightBatch(ip)
          enddo
          call dgemm('N','T',n_sig,n_sig,nb,1.0d0,ValW1a,nconts,ValBatch_sig,nconts,0.0d0,Fxc_a_s,nconts)
          if (mgga_active) then
             do dd = 1,3
                Val1d_tmp(1:n_sig,1:nb) = val1_sig(1:n_sig,dd,1:nb)
                do ip = 1,nb
                   ValTauW(1:n_sig,ip) = 0.25d0*Val1d_tmp(1:n_sig,ip)*(vtaua_b(ip)*WeightBatch(ip))
                enddo
                call dgemm('N','T',n_sig,n_sig,nb,1.0d0,ValTauW,nconts,Val1d_tmp,nconts,1.0d0,Fxc_a_s,nconts)
             enddo
          endif
          do jj = 1,n_sig
             do ii = 1,jj
                Fxc_a_s_val = Fxc_a_s(ii,jj) + Fxc_a_s(jj,ii)
                Fxc_a_local(sig_idx(ii),sig_idx(jj)) = Fxc_a_local(sig_idx(ii),sig_idx(jj)) + Fxc_a_s_val
                if (ii .ne. jj) Fxc_a_local(sig_idx(jj),sig_idx(ii)) = Fxc_a_local(sig_idx(jj),sig_idx(ii)) + Fxc_a_s_val
             enddo
          enddo
       else
          do ip = 1,nb
             ValW1a(1:n_sig,ip) = 0.5d0*ValBatch_sig(1:n_sig,ip)*(vrhoa_b(ip)*WeightBatch(ip)) &
                                + Ta(1:n_sig,ip)*WeightBatch(ip)
             ValW1b(1:n_sig,ip) = 0.5d0*ValBatch_sig(1:n_sig,ip)*(vrhob_b(ip)*WeightBatch(ip)) &
                                + Tb(1:n_sig,ip)*WeightBatch(ip)
          enddo
          call dgemm('N','T',n_sig,n_sig,nb,1.0d0,ValW1a,nconts,ValBatch_sig,nconts,0.0d0,Fxc_a_s,nconts)
          call dgemm('N','T',n_sig,n_sig,nb,1.0d0,ValW1b,nconts,ValBatch_sig,nconts,0.0d0,Fxc_b_s,nconts)
          if (mgga_active) then
             do dd = 1,3
                Val1d_tmp(1:n_sig,1:nb) = val1_sig(1:n_sig,dd,1:nb)
                do ip = 1,nb
                   ValTauW(1:n_sig,ip) = 0.25d0*Val1d_tmp(1:n_sig,ip)*(vtaua_b(ip)*WeightBatch(ip))
                enddo
                call dgemm('N','T',n_sig,n_sig,nb,1.0d0,ValTauW,nconts,Val1d_tmp,nconts,1.0d0,Fxc_a_s,nconts)
                do ip = 1,nb
                   ValTauW(1:n_sig,ip) = 0.25d0*Val1d_tmp(1:n_sig,ip)*(vtaub_b(ip)*WeightBatch(ip))
                enddo
                call dgemm('N','T',n_sig,n_sig,nb,1.0d0,ValTauW,nconts,Val1d_tmp,nconts,1.0d0,Fxc_b_s,nconts)
             enddo
          endif
          do jj = 1,n_sig
             do ii = 1,jj
                Fxc_a_s_val = Fxc_a_s(ii,jj) + Fxc_a_s(jj,ii)
                Fxc_b_s_val = Fxc_b_s(ii,jj) + Fxc_b_s(jj,ii)
                Fxc_a_local(sig_idx(ii),sig_idx(jj)) = Fxc_a_local(sig_idx(ii),sig_idx(jj)) + Fxc_a_s_val
                Fxc_b_local(sig_idx(ii),sig_idx(jj)) = Fxc_b_local(sig_idx(ii),sig_idx(jj)) + Fxc_b_s_val
                if (ii .ne. jj) then
                   Fxc_a_local(sig_idx(jj),sig_idx(ii)) = Fxc_a_local(sig_idx(jj),sig_idx(ii)) + Fxc_a_s_val
                   Fxc_b_local(sig_idx(jj),sig_idx(ii)) = Fxc_b_local(sig_idx(jj),sig_idx(ii)) + Fxc_b_s_val
                endif
             enddo
          enddo
       endif
    enddo
    !$omp end do

    deallocate(WeightBatch)
    deallocate(ValBatch_sig)
    deallocate(val1_sig)
    deallocate(Wa,rho_a_b,drv_a_b,d1rho_b,d1sig_b,Et_b,Ta)
    deallocate(ValW1a)
    deallocate(ex_b,ec_b,sigma_t_b,sigma_aa_b,sigma_ab_b,sigma_bb_b)
    deallocate(vrhoa_b,vrhob_b,vsigmaaa_b,vsigmabb_b,vsigmaab_b)
    deallocate(Pa_s,Fxc_a_s)
    if (allocated(Wb)) deallocate(Wb)
    if (allocated(rho_b_b)) deallocate(rho_b_b,drv_b_b)
    if (allocated(Tb)) deallocate(Tb)
    if (allocated(ValW1b)) deallocate(ValW1b)
    if (allocated(catPab)) deallocate(catPab,catWab)
    if (allocated(Pb_s)) deallocate(Pb_s,Fxc_b_s)
    if (allocated(tau_a_b)) deallocate(tau_a_b,vtaua_b,vtaub_b)
    if (allocated(tau_b_b)) deallocate(tau_b_b)
    if (allocated(Val1d_tmp)) deallocate(Val1d_tmp,Wa1_tmp,ValTauW)
    deallocate(ValBatch_own,Val1Batch_own)
    deallocate(ValBatch_wide,Val1Batch_wide,maxval_bf)
    deallocate(sig_idx,narrow_idx,wide_idx)

    block
      integer :: merge_tid, merge_nthreads
      merge_nthreads = omp_get_num_threads()
      do merge_tid = 0, merge_nthreads-1
         !$omp barrier
         if (omp_get_thread_num() .eq. merge_tid) then
            Exc = Exc + Exc_local
            Fxc_a = Fxc_a + Fxc_a_local
            Fxc_b = Fxc_b + Fxc_b_local
         endif
      enddo
    end block
    deallocate(Fxc_a_local,Fxc_b_local)
    !$omp end parallel

    if (Multi .eq. 1) Fxc_b = Fxc_a
end subroutine
