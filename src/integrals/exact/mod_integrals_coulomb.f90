! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Exact (non-fitted) 4-center Coulomb matrix build via direct ERI quartets.

submodule (mod_integrals) coulomb_impl
implicit none
contains

module subroutine integrals_build_coulomb(nConts, Ptot, J)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    real(8),allocatable :: dP(:,:), dJ(:,:)
    character(len=8) :: envval
    logical :: do_full
    integer :: ios

    if (.not. incr_fock_env_checked) then
       envval = ""
       call get_environment_variable("ENGINE_INCREMENTAL_FOCK", envval)
       if (trim(envval) .eq. "1") incr_fock_enabled = .true.
       if (trim(envval) .eq. "0") incr_fock_enabled = .false.
       envval = ""
       call get_environment_variable("ENGINE_INCR_FULL_PERIOD", envval)
       if (len_trim(envval) .gt. 0) read(envval,*,iostat=ios) incr_full_period
       incr_fock_env_checked = .true.
    endif

    do_full = (.not. incr_fock_enabled) .or. (.not. incr_has_prev) &
              .or. (incr_since_full .ge. incr_full_period)

    if (do_full) then
       if (direct_mode) then
          call build_coulomb_direct(nConts, Ptot, J)
       else
          call build_coulomb_store(nConts, Ptot, J)
       endif
       incr_since_full = 0
    else
       allocate(dP(nConts,nConts), dJ(nConts,nConts))
       dP = Ptot - incr_Ptot_prev
       if (direct_mode) then
          call build_coulomb_direct(nConts, dP, dJ, is_incremental=.true.)
       else
          call build_coulomb_store(nConts, dP, dJ)
       endif
       J = incr_J_prev + dJ
       deallocate(dP, dJ)
       incr_since_full = incr_since_full + 1
    endif

    if (incr_fock_enabled) then
       if (.not. allocated(incr_Ptot_prev)) then
          allocate(incr_Ptot_prev(nConts,nConts), incr_J_prev(nConts,nConts))
       endif
       incr_Ptot_prev = Ptot
       incr_J_prev = J
       incr_has_prev = .true.
    endif
end subroutine integrals_build_coulomb

module logical function is_canonical_quartet(e1,e2,e3,e4) result(ok)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4
    integer :: pq,rs
    if (e1 < e2 .or. e3 < e4) then
       ok = .false.
       return
    endif
    pq = e1*(e1+1)/2 + e2
    rs = e3*(e3+1)/2 + e4
    ok = (pq >= rs)
end function is_canonical_quartet

module subroutine accumulate_quartet_J(e1,e2,e3,e4,val,nConts,Ptot,J)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4,nConts
    real(8),intent(in) :: val
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(inout) :: J(nConts,nConts)
    logical :: s12,s34,sbk
    integer :: a1,a2,a3,a4

    a1 = e1+1; a2 = e2+1; a3 = e3+1; a4 = e4+1
    s12 = (e1 .eq. e2)
    s34 = (e3 .eq. e4)
    sbk = (e1 .eq. e3) .and. (e2 .eq. e4)

    J(a1,a2) = J(a1,a2) + val*Ptot(a3,a4)
    if (.not. s12) J(a2,a1) = J(a2,a1) + val*Ptot(a3,a4)
    if (.not. s34) J(a1,a2) = J(a1,a2) + val*Ptot(a4,a3)
    if ((.not. s12) .and. (.not. s34)) J(a2,a1) = J(a2,a1) + val*Ptot(a4,a3)
    if (.not. sbk) then
       J(a3,a4) = J(a3,a4) + val*Ptot(a1,a2)
       if (.not. s34) J(a4,a3) = J(a4,a3) + val*Ptot(a1,a2)
       if (.not. s12) J(a3,a4) = J(a3,a4) + val*Ptot(a2,a1)
       if ((.not. s12) .and. (.not. s34)) J(a4,a3) = J(a4,a3) + val*Ptot(a2,a1)
    endif
end subroutine accumulate_quartet_J

module subroutine expand_quartet_tuples(e1,e2,e3,e4,n_out,idx_out)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4
    integer,intent(out) :: n_out
    integer,intent(out) :: idx_out(4,8)
    integer :: perm(4,8)
    logical :: seen(8)
    integer :: p,q

    perm(:,1) = (/ e1,e2,e3,e4 /)
    perm(:,2) = (/ e2,e1,e3,e4 /)
    perm(:,3) = (/ e1,e2,e4,e3 /)
    perm(:,4) = (/ e2,e1,e4,e3 /)
    perm(:,5) = (/ e3,e4,e1,e2 /)
    perm(:,6) = (/ e4,e3,e1,e2 /)
    perm(:,7) = (/ e3,e4,e2,e1 /)
    perm(:,8) = (/ e4,e3,e2,e1 /)

    seen = .false.
    n_out = 0
    do p = 1,8
       if (seen(p)) cycle
       do q = p+1,8
          if (.not.seen(q)) then
             if (all(perm(:,q) == perm(:,p))) seen(q) = .true.
          endif
       enddo
       n_out = n_out+1
       idx_out(:,n_out) = perm(:,p)
    enddo
end subroutine expand_quartet_tuples

module subroutine build_coulomb_pk_list(nConts)
    implicit none
    integer,intent(in) :: nConts
    integer :: p,q,r,s,pq,rs
    real(8) :: val
    integer(8) :: dbg_total, dbg_skipped

    coulomb_npairs = nConts*(nConts+1)/2
    dbg_total = 0
    dbg_skipped = 0
    cPK_n = 0
    do p = 1,nConts
       do q = 1,p
          pq = p*(p-1)/2 + q
          do r = 1,p
             do s = 1,r
                rs = r*(r-1)/2 + s
                if (rs > pq) cycle
                if (schwarz_bound(ao_shell(p),ao_shell(q))*schwarz_bound(ao_shell(r),ao_shell(s)) &
                    .lt. SCHWARZ_CUTOFF) cycle
                cPK_n = cPK_n + 1
             enddo
          enddo
       enddo
    enddo

    allocate(cPK_pq(cPK_n),cPK_rs(cPK_n),cPK_val(cPK_n))

    cPK_n = 0
    do p = 1,nConts
       do q = 1,p
          pq = p*(p-1)/2 + q
          do r = 1,p
             do s = 1,r
                rs = r*(r-1)/2 + s
                if (rs > pq) cycle
                dbg_total = dbg_total + 1
                if (schwarz_bound(ao_shell(p),ao_shell(q))*schwarz_bound(ao_shell(r),ao_shell(s)) &
                    .lt. SCHWARZ_CUTOFF) then
                   dbg_skipped = dbg_skipped + 1
                   cycle
                endif
                val = eri_get(p,q,r,s)
                cPK_n = cPK_n + 1
                cPK_pq(cPK_n) = pq
                cPK_rs(cPK_n) = rs
                cPK_val(cPK_n) = val
             enddo
          enddo
       enddo
    enddo
    print *,"DIAG: build_coulomb_store(packed) visited=",dbg_total," skipped=",dbg_skipped, &
            " terms=",cPK_n," (vs nConts^2/2=",coulomb_npairs,")"
end subroutine build_coulomb_pk_list

module subroutine build_coulomb_store(nConts, Ptot, J)
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    integer(8) :: t
    integer :: p,q,pq
    real(8),allocatable :: D_vec(:),J_vec(:)
    real(8),allocatable :: J_vec_local(:)

    if (.not. coulomb_list_built) then
       call build_coulomb_pk_list(nConts)
       coulomb_list_built = .true.
    endif

    allocate(D_vec(coulomb_npairs),J_vec(coulomb_npairs))

    do p = 1,nConts
       do q = 1,p
          pq = p*(p-1)/2 + q
          if (p .eq. q) then
             D_vec(pq) = Ptot(p,q)
          else
             D_vec(pq) = 2.0d0*Ptot(p,q)
          endif
       enddo
    enddo

    J_vec = 0
    !$omp parallel default(shared) private(J_vec_local)
    allocate(J_vec_local(coulomb_npairs))
    J_vec_local = 0.0d0
    !$omp do schedule(static)
    do t = 1,cPK_n
       J_vec_local(cPK_pq(t)) = J_vec_local(cPK_pq(t)) + cPK_val(t)*D_vec(cPK_rs(t))
       if (cPK_rs(t) .ne. cPK_pq(t)) then
          J_vec_local(cPK_rs(t)) = J_vec_local(cPK_rs(t)) + cPK_val(t)*D_vec(cPK_pq(t))
       endif
    enddo
    !$omp end do
    block
       integer :: merge_tid, merge_nthreads
       merge_nthreads = omp_get_num_threads()
       do merge_tid = 0, merge_nthreads-1
          !$omp barrier
          if (omp_get_thread_num() .eq. merge_tid) then
             J_vec = J_vec + J_vec_local
          endif
       enddo
    end block
    deallocate(J_vec_local)
    !$omp end parallel

    do p = 1,nConts
       do q = 1,p
          pq = p*(p-1)/2 + q
          J(p,q) = J_vec(pq)
          J(q,p) = J_vec(pq)
       enddo
    enddo

    deallocate(D_vec,J_vec)
end subroutine build_coulomb_store

module subroutine build_coulomb_direct(nConts, Ptot, J, is_incremental)
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    logical,intent(in),optional :: is_incremental
    logical :: incremental_feed
    real(8) :: cut_now

    integer :: ijsh,klsh,si,sj,sk,sl,di,dj,dk,dl
    integer :: is,js,ks,ls,num,ii
    integer :: shls(4)
    real(8),allocatable :: buf2e(:,:,:,:)
    real(8),allocatable :: J_local(:,:)
    integer :: p,q,r,s,e1,e2,e3,e4
    integer :: a1,a2,a3,a4
    real(8) :: val
    real(8),allocatable :: dmax_shell(:,:)
    integer,allocatable :: ao_offset(:),pq_min_shell(:),pq_max_shell(:)
    integer :: e1t,e2t,pqv
    integer(8) :: n_visit,n_cut_sw,n_cut_dn,n_call
    integer(8) :: n_tight(0:15)
    logical :: tightness_stats
    character(len=8) :: ts_env
    integer :: ts_stat, ib

    n_visit = 0; n_cut_sw = 0; n_cut_dn = 0; n_call = 0
    n_tight = 0
    ts_env = ""
    call get_environment_variable('ENGINE_JK_STATS', ts_env, status=ts_stat)
    tightness_stats = (ts_stat .eq. 0 .and. trim(ts_env) .eq. '2')
    J = 0
    incremental_feed = present(is_incremental)
    if (incremental_feed) incremental_feed = is_incremental
    cut_now = merge(density_screen_cutoff_incr_j, density_screen_cutoff, incremental_feed)
    block
      real(8) :: t0,t1
      call cpu_time(t0)
    allocate(dmax_shell(0:nBases-1,0:nBases-1))
    call compute_shell_density_bound(nConts, Ptot, dmax_shell)
    call build_significant_pairs()
      call cpu_time(t1)
      if (tightness_stats) write(6,'(a,f8.4,a)') ' SERIAL1 (density_bound+sig_pairs): ',t1-t0,' s'
    end block

    allocate(ao_offset(0:nBases-1))
    num = 0
    do si = 0,nBases-1
       ao_offset(si) = num
       num = num + cgto_engine(si,bas)
    enddo
    block
      real(8) :: t0,t1
      call cpu_time(t0)
    allocate(pq_min_shell(0:n_sig_pairs-1),pq_max_shell(0:n_sig_pairs-1))
    do ijsh = 0,n_sig_pairs-1
       si = sig_pair_ishl(ijsh); sj = sig_pair_jshl(ijsh)
       di = cgto_engine(si,bas); dj = cgto_engine(sj,bas)
       pq_min_shell(ijsh) = huge(1)
       pq_max_shell(ijsh) = -1
       do e1t = ao_offset(si), ao_offset(si)+di-1
          do e2t = ao_offset(sj), ao_offset(sj)+dj-1
             if (e1t .ge. e2t) then
                pqv = e1t*(e1t+1)/2 + e2t
                pq_min_shell(ijsh) = min(pq_min_shell(ijsh),pqv)
                pq_max_shell(ijsh) = max(pq_max_shell(ijsh),pqv)
             endif
          enddo
       enddo
    enddo
      call cpu_time(t1)
      if (tightness_stats) write(6,'(a,f8.4,a)') ' SERIAL2 (pq_min/max_shell): ',t1-t0,' s'
    end block

    !$omp parallel default(shared) private(si,sj,di,dj,klsh,sk,sl,dk,dl,shls,buf2e) &
    !$omp&   private(is,js,ks,ls,p,q,r,s,e1,e2,e3,e4,val,a1,a2,a3,a4) &
    !$omp&   private(J_local) reduction(+:n_visit,n_cut_sw,n_cut_dn,n_call,n_tight)
    allocate(J_local(nConts,nConts))
    J_local = 0.0d0
    !$omp do schedule(static,1)
    do ijsh = 0, n_sig_pairs-1
       si = sig_pair_ishl(ijsh)
       sj = sig_pair_jshl(ijsh)
       di = cgto_engine(si, bas)
       dj = cgto_engine(sj, bas)
       do klsh = 0, ijsh
          sk = sig_pair_ishl(klsh)
          sl = sig_pair_jshl(klsh)
          n_visit = n_visit + 1
          if (schwarz_bound(si,sj)*schwarz_bound(sk,sl) .lt. SCHWARZ_CUTOFF) then
             n_cut_sw = n_cut_sw + 1
             cycle
          endif
          if (schwarz_bound(si,sj)*schwarz_bound(sk,sl) &
              *max(dmax_shell(si,sj),dmax_shell(sk,sl)) .lt. cut_now) then
             n_cut_dn = n_cut_dn + 1
             cycle
          endif
          n_call = n_call + 1
          dk = cgto_engine(sk, bas)
          dl = cgto_engine(sl, bas)
          shls(1) = si; shls(2) = sj; shls(3) = sk; shls(4) = sl
          allocate(buf2e(di,dj,dk,dl))
          call twoe_engine(buf2e, shls, atm, size(atm,2), bas, nBases, env, int2e_opt)

          if (tightness_stats) then
             block
                real(8) :: amax, ratio
                integer :: bkt
                amax = maxval(abs(buf2e))
                if (amax .gt. 0.0d0) then
                   ratio = schwarz_bound(si,sj)*schwarz_bound(sk,sl)/amax
                   bkt = min(15, max(0, int(log10(max(ratio,1.0d0)))))
                else
                   bkt = 15
                endif
                n_tight(bkt) = n_tight(bkt) + 1
             end block
          endif

          is = ao_offset(si); js = ao_offset(sj); ks = ao_offset(sk); ls = ao_offset(sl)

          if (norvec_is_identity .and. klsh .ne. ijsh &
              .and. is .gt. js+dj-1 .and. ks .gt. ls+dl-1 &
              .and. pq_max_shell(klsh) .lt. pq_min_shell(ijsh)) then
             do p = 1,di
                a1 = is+p
                do q = 1,dj
                   a2 = js+q
                   do r = 1,dk
                      a3 = ks+r
                      do s = 1,dl
                         a4 = ls+s
                         val = buf2e(p,q,r,s)
                         J_local(a1,a2) = J_local(a1,a2) + val*Ptot(a3,a4)
                         J_local(a2,a1) = J_local(a2,a1) + val*Ptot(a3,a4)
                         J_local(a1,a2) = J_local(a1,a2) + val*Ptot(a4,a3)
                         J_local(a2,a1) = J_local(a2,a1) + val*Ptot(a4,a3)
                         J_local(a3,a4) = J_local(a3,a4) + val*Ptot(a1,a2)
                         J_local(a4,a3) = J_local(a4,a3) + val*Ptot(a1,a2)
                         J_local(a3,a4) = J_local(a3,a4) + val*Ptot(a2,a1)
                         J_local(a4,a3) = J_local(a4,a3) + val*Ptot(a2,a1)
                      enddo
                   enddo
                enddo
             enddo
          else if (norvec_is_identity) then
             do p = 1,di
                do q = 1,dj
                   do r = 1,dk
                      do s = 1,dl
                         e1 = is+p-1; e2 = js+q-1; e3 = ks+r-1; e4 = ls+s-1
                         val = buf2e(p,q,r,s)
                         if (is_canonical_quartet(e1,e2,e3,e4)) then
                            call accumulate_quartet_J(e1,e2,e3,e4,val,nConts,Ptot,J_local)
                         endif
                         if (klsh .ne. ijsh) then
                            if (is_canonical_quartet(e3,e4,e1,e2)) then
                               call accumulate_quartet_J(e3,e4,e1,e2,val,nConts,Ptot,J_local)
                            endif
                         endif
                      enddo
                   enddo
                enddo
             enddo
          else
          do p = 1,di
             do q = 1,dj
                do r = 1,dk
                   do s = 1,dl
                      e1 = is+p-1; e2 = js+q-1; e3 = ks+r-1; e4 = ls+s-1
                      val = buf2e(p,q,r,s)*NorVEC(e1+1)*NorVEC(e2+1)*NorVEC(e3+1)*NorVEC(e4+1)
                      if (is_canonical_quartet(e1,e2,e3,e4)) then
                         call accumulate_quartet_J(e1,e2,e3,e4,val,nConts,Ptot,J_local)
                      endif
                      if (klsh .ne. ijsh) then
                         if (is_canonical_quartet(e3,e4,e1,e2)) then
                            call accumulate_quartet_J(e3,e4,e1,e2,val,nConts,Ptot,J_local)
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
          endif
          deallocate(buf2e)
       enddo
    enddo
    !$omp end do
    block
       integer :: merge_tid, merge_nthreads
       merge_nthreads = omp_get_num_threads()
       do merge_tid = 0, merge_nthreads-1
          !$omp barrier
          if (omp_get_thread_num() .eq. merge_tid) then
             J = J + J_local
          endif
       enddo
    end block
    deallocate(J_local)
    !$omp end parallel
    block
       character(len=8) :: stats_env
       integer :: stats_stat
       call get_environment_variable('ENGINE_JK_STATS', stats_env, status=stats_stat)
       if (tightness_stats) then
          write(6,'(a)') ' JTIGHT  Schwarz bound / actual max|(ij|kl)|, per decade:'
          do ib = 0,15
             if (n_tight(ib) .gt. 0) write(6,'(a,i2,a,i12)') &
                '   >=1e',ib,' :',n_tight(ib)
          enddo
          flush(6)
       endif
       if (stats_stat .eq. 0 .and. (trim(stats_env) .eq. '1' .or. trim(stats_env) .eq. '2')) then
          write(6,'(a,i4,a,i12,a,i12,a,i12,a,i12,a,f6.2,a)') &
             ' JSTATS pairs=',n_sig_pairs,' visited=',n_visit, &
             ' cut_schwarz=',n_cut_sw,' cut_density=',n_cut_dn, &
             ' cint_calls=',n_call,' (',100.0d0*real(n_call,8)/max(1.0d0,real(n_visit,8)),'% survive)'
          flush(6)
       endif
    end block
    deallocate(dmax_shell)
    deallocate(ao_offset,pq_min_shell,pq_max_shell)
end subroutine build_coulomb_direct

module logical function jk_fusion_active() result(ok)
    implicit none
    ok = direct_mode .and. (.not. engine_use_df_j) .and. (.not. engine_use_df_k)
end function jk_fusion_active

subroutine jk_fused_core(nConts, Ptot, Da, Db, J, Ka, Kb, is_incremental)
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    logical,intent(in) :: is_incremental

    integer :: ijsh,klsh,si,sj,sk,sl,di,dj,dk,dl
    integer :: is,js,ks,ls,num
    integer :: shls(4)
    real(8),allocatable :: buf2e(:,:,:,:)
    real(8),allocatable :: J_local(:,:), Ka_local(:,:), Kb_local(:,:)
    integer :: p,q,r,s,e1,e2,e3,e4
    real(8) :: val, sw, bj, bk, cut_now
    real(8) :: nrm1, nrm12, nrm123
    integer :: pq_idx, rs_idx
    logical :: bra_ok, ket_ok
    real(8),allocatable :: dmaxJ(:,:), dmaxK(:,:), Dsum(:,:)
    integer,allocatable :: ao_offset(:)
    logical :: keepJ, keepK
    integer,allocatable :: pq_min_shell(:),pq_max_shell(:)
    integer :: e1t,e2t,pqv,a1,a2,a3,a4

    J = 0.0d0
    Ka = 0.0d0
    Kb = 0.0d0

    call build_significant_pairs()

    allocate(dmaxJ(0:nBases-1,0:nBases-1), dmaxK(0:nBases-1,0:nBases-1))
    allocate(Dsum(nConts,nConts))
    Dsum = abs(Da) + abs(Db)
    call compute_shell_density_bound(nConts, Ptot, dmaxJ)
    call compute_shell_density_bound(nConts, Dsum, dmaxK)
    deallocate(Dsum)

    allocate(ao_offset(0:nBases-1))
    num = 0
    do si = 0,nBases-1
       ao_offset(si) = num
       num = num + cgto_engine(si,bas)
    enddo
    allocate(pq_min_shell(0:n_sig_pairs-1),pq_max_shell(0:n_sig_pairs-1))
    do ijsh = 0,n_sig_pairs-1
       si = sig_pair_ishl(ijsh); sj = sig_pair_jshl(ijsh)
       di = cgto_engine(si,bas); dj = cgto_engine(sj,bas)
       pq_min_shell(ijsh) = huge(1)
       pq_max_shell(ijsh) = -1
       do e1t = ao_offset(si), ao_offset(si)+di-1
          do e2t = ao_offset(sj), ao_offset(sj)+dj-1
             if (e1t .ge. e2t) then
                pqv = e1t*(e1t+1)/2 + e2t
                pq_min_shell(ijsh) = min(pq_min_shell(ijsh),pqv)
                pq_max_shell(ijsh) = max(pq_max_shell(ijsh),pqv)
             endif
          enddo
       enddo
    enddo

    !$omp parallel default(shared) private(si,sj,di,dj,klsh,sk,sl,dk,dl,shls,buf2e) &
    !$omp&   private(is,js,ks,ls,p,q,r,s,e1,e2,e3,e4,val,sw,bj,bk,cut_now,keepJ,keepK) &
    !$omp&   private(nrm1,nrm12,nrm123,pq_idx,rs_idx,bra_ok,ket_ok,a1,a2,a3,a4) &
    !$omp&   private(J_local,Ka_local,Kb_local)
    allocate(J_local(nConts,nConts), Ka_local(nConts,nConts), Kb_local(nConts,nConts))
    J_local = 0.0d0
    Ka_local = 0.0d0
    Kb_local = 0.0d0
    !$omp do schedule(static,1)
    do ijsh = 0, n_sig_pairs-1
       si = sig_pair_ishl(ijsh)
       sj = sig_pair_jshl(ijsh)
       di = cgto_engine(si, bas)
       dj = cgto_engine(sj, bas)
       do klsh = 0, ijsh
          sk = sig_pair_ishl(klsh)
          sl = sig_pair_jshl(klsh)
          sw = schwarz_bound(si,sj)*schwarz_bound(sk,sl)
          if (sw .lt. SCHWARZ_CUTOFF) cycle
          bj = sw*max(dmaxJ(si,sj),dmaxJ(sk,sl))
          bk = sw*max(max(dmaxK(si,sk),dmaxK(sj,sl)), &
                      max(dmaxK(si,sl),dmaxK(sj,sk)))
          cut_now = merge(density_screen_cutoff_incr, density_screen_cutoff, is_incremental)
          keepJ = (bj .ge. cut_now)
          keepK = (bk .ge. cut_now)
          if (.not. (keepJ .or. keepK)) cycle
          dk = cgto_engine(sk, bas)
          dl = cgto_engine(sl, bas)
          shls(1) = si; shls(2) = sj; shls(3) = sk; shls(4) = sl
          allocate(buf2e(di,dj,dk,dl))
          call twoe_engine(buf2e, shls, atm, size(atm,2), bas, nBases, env, int2e_opt)
          is = ao_offset(si); js = ao_offset(sj); ks = ao_offset(sk); ls = ao_offset(sl)
          if (norvec_is_identity .and. klsh .ne. ijsh &
              .and. is .gt. js+dj-1 .and. ks .gt. ls+dl-1 &
              .and. pq_max_shell(klsh) .lt. pq_min_shell(ijsh)) then
             do p = 1,di
                a1 = is+p
                do q = 1,dj
                   a2 = js+q
                   do r = 1,dk
                      a3 = ks+r
                      do s = 1,dl
                         a4 = ls+s
                         val = buf2e(p,q,r,s)
                         if (keepJ) then
                            J_local(a1,a2) = J_local(a1,a2) + val*Ptot(a3,a4)
                            J_local(a2,a1) = J_local(a2,a1) + val*Ptot(a3,a4)
                            J_local(a1,a2) = J_local(a1,a2) + val*Ptot(a4,a3)
                            J_local(a2,a1) = J_local(a2,a1) + val*Ptot(a4,a3)
                            J_local(a3,a4) = J_local(a3,a4) + val*Ptot(a1,a2)
                            J_local(a4,a3) = J_local(a4,a3) + val*Ptot(a1,a2)
                            J_local(a3,a4) = J_local(a3,a4) + val*Ptot(a2,a1)
                            J_local(a4,a3) = J_local(a4,a3) + val*Ptot(a2,a1)
                         endif
                         if (keepK) then
                            Ka_local(a1,a3) = Ka_local(a1,a3) + val*Da(a2,a4)
                            Kb_local(a1,a3) = Kb_local(a1,a3) + val*Db(a2,a4)
                            Ka_local(a2,a3) = Ka_local(a2,a3) + val*Da(a1,a4)
                            Kb_local(a2,a3) = Kb_local(a2,a3) + val*Db(a1,a4)
                            Ka_local(a1,a4) = Ka_local(a1,a4) + val*Da(a2,a3)
                            Kb_local(a1,a4) = Kb_local(a1,a4) + val*Db(a2,a3)
                            Ka_local(a2,a4) = Ka_local(a2,a4) + val*Da(a1,a3)
                            Kb_local(a2,a4) = Kb_local(a2,a4) + val*Db(a1,a3)
                            Ka_local(a3,a1) = Ka_local(a3,a1) + val*Da(a4,a2)
                            Kb_local(a3,a1) = Kb_local(a3,a1) + val*Db(a4,a2)
                            Ka_local(a4,a1) = Ka_local(a4,a1) + val*Da(a3,a2)
                            Kb_local(a4,a1) = Kb_local(a4,a1) + val*Db(a3,a2)
                            Ka_local(a3,a2) = Ka_local(a3,a2) + val*Da(a4,a1)
                            Kb_local(a3,a2) = Kb_local(a3,a2) + val*Db(a4,a1)
                            Ka_local(a4,a2) = Ka_local(a4,a2) + val*Da(a3,a1)
                            Kb_local(a4,a2) = Kb_local(a4,a2) + val*Db(a3,a1)
                         endif
                      enddo
                   enddo
                enddo
             enddo
          else
          do p = 1,di
             e1 = is+p-1
             nrm1 = NorVEC(e1+1)
             do q = 1,dj
                e2 = js+q-1
                nrm12 = nrm1*NorVEC(e2+1)
                pq_idx = e1*(e1+1)/2 + e2
                bra_ok = (e1 .ge. e2)
                do r = 1,dk
                   e3 = ks+r-1
                   nrm123 = nrm12*NorVEC(e3+1)
                   do s = 1,dl
                      e4 = ls+s-1
                      rs_idx = e3*(e3+1)/2 + e4
                      ket_ok = (e3 .ge. e4)
                      val = buf2e(p,q,r,s)*nrm123*NorVEC(e4+1)
                      if (bra_ok .and. ket_ok .and. pq_idx .ge. rs_idx) then
                         call accumulate_quartet_jk(e1,e2,e3,e4,val,nConts,Ptot,Da,Db, &
                                                    J_local,Ka_local,Kb_local,keepJ,keepK)
                      endif
                      if (klsh .ne. ijsh) then
                         if (ket_ok .and. bra_ok .and. rs_idx .ge. pq_idx) then
                            call accumulate_quartet_jk(e3,e4,e1,e2,val,nConts,Ptot,Da,Db, &
                                                       J_local,Ka_local,Kb_local,keepJ,keepK)
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
          endif
          deallocate(buf2e)
       enddo
    enddo
    !$omp end do
    block
       integer :: merge_tid, merge_nthreads
       merge_nthreads = omp_get_num_threads()
       do merge_tid = 0, merge_nthreads-1
          !$omp barrier
          if (omp_get_thread_num() .eq. merge_tid) then
             J = J + J_local
             Ka = Ka + Ka_local
             Kb = Kb + Kb_local
          endif
       enddo
    end block
    deallocate(J_local,Ka_local,Kb_local)
    !$omp end parallel

    deallocate(dmaxJ,dmaxK,ao_offset,pq_min_shell,pq_max_shell)
end subroutine jk_fused_core

module subroutine integrals_build_jk_fused(nConts, Ptot, Da, Db, J, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    real(8),allocatable :: dP(:,:), dDa(:,:), dDb(:,:)
    real(8),allocatable :: dJ(:,:), dKa(:,:), dKb(:,:)
    logical :: do_full

    do_full = (.not. incr_fock_enabled) .or. (.not. incrf_has_prev) &
              .or. (incrf_since_full .ge. incr_full_period)

    if (do_full) then
       call jk_fused_core(nConts, Ptot, Da, Db, J, Ka, Kb, .false.)
       incrf_since_full = 0
    else
       allocate(dP(nConts,nConts), dDa(nConts,nConts), dDb(nConts,nConts))
       allocate(dJ(nConts,nConts), dKa(nConts,nConts), dKb(nConts,nConts))
       dP  = Ptot - (incrf_Da_prev + incrf_Db_prev)
       dDa = Da - incrf_Da_prev
       dDb = Db - incrf_Db_prev
       call jk_fused_core(nConts, dP, dDa, dDb, dJ, dKa, dKb, .true.)
       J  = incrf_J_prev  + dJ
       Ka = incrf_Ka_prev + dKa
       Kb = incrf_Kb_prev + dKb
       deallocate(dP,dDa,dDb,dJ,dKa,dKb)
       incrf_since_full = incrf_since_full + 1
    endif

    if (incr_fock_enabled) then
       if (.not. allocated(incrf_Da_prev)) then
          allocate(incrf_Da_prev(nConts,nConts), incrf_Db_prev(nConts,nConts))
          allocate(incrf_J_prev(nConts,nConts))
          allocate(incrf_Ka_prev(nConts,nConts), incrf_Kb_prev(nConts,nConts))
       endif
       incrf_Da_prev = Da
       incrf_Db_prev = Db
       incrf_J_prev  = J
       incrf_Ka_prev = Ka
       incrf_Kb_prev = Kb
       incrf_has_prev = .true.
    endif
end subroutine integrals_build_jk_fused

subroutine accumulate_quartet_jk(e1,e2,e3,e4,val,nConts,Ptot,Da,Db,J,Ka,Kb,doJ,doK)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4,nConts
    real(8),intent(in) :: val
    real(8),intent(in) :: Ptot(nConts,nConts), Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(inout) :: J(nConts,nConts), Ka(nConts,nConts), Kb(nConts,nConts)
    logical,intent(in) :: doJ, doK
    logical :: s12,s34,sbk,n12,n34
    integer :: a1,a2,a3,a4

    a1 = e1+1; a2 = e2+1; a3 = e3+1; a4 = e4+1
    s12 = (e1 .eq. e2);  n12 = .not. s12
    s34 = (e3 .eq. e4);  n34 = .not. s34
    sbk = (e1 .eq. e3) .and. (e2 .eq. e4)

    if (doJ) then
       J(a1,a2) = J(a1,a2) + val*Ptot(a3,a4)
       if (n12) J(a2,a1) = J(a2,a1) + val*Ptot(a3,a4)
       if (n34) J(a1,a2) = J(a1,a2) + val*Ptot(a4,a3)
       if (n12 .and. n34) J(a2,a1) = J(a2,a1) + val*Ptot(a4,a3)
       if (.not. sbk) then
          J(a3,a4) = J(a3,a4) + val*Ptot(a1,a2)
          if (n34) J(a4,a3) = J(a4,a3) + val*Ptot(a1,a2)
          if (n12) J(a3,a4) = J(a3,a4) + val*Ptot(a2,a1)
          if (n12 .and. n34) J(a4,a3) = J(a4,a3) + val*Ptot(a2,a1)
       endif
    endif

    if (doK) then
       Ka(a1,a3) = Ka(a1,a3) + val*Da(a2,a4)
       Kb(a1,a3) = Kb(a1,a3) + val*Db(a2,a4)
       if (n12) then
          Ka(a2,a3) = Ka(a2,a3) + val*Da(a1,a4)
          Kb(a2,a3) = Kb(a2,a3) + val*Db(a1,a4)
       endif
       if (n34) then
          Ka(a1,a4) = Ka(a1,a4) + val*Da(a2,a3)
          Kb(a1,a4) = Kb(a1,a4) + val*Db(a2,a3)
       endif
       if (n12 .and. n34) then
          Ka(a2,a4) = Ka(a2,a4) + val*Da(a1,a3)
          Kb(a2,a4) = Kb(a2,a4) + val*Db(a1,a3)
       endif
       if (.not. sbk) then
          Ka(a3,a1) = Ka(a3,a1) + val*Da(a4,a2)
          Kb(a3,a1) = Kb(a3,a1) + val*Db(a4,a2)
          if (n34) then
             Ka(a4,a1) = Ka(a4,a1) + val*Da(a3,a2)
             Kb(a4,a1) = Kb(a4,a1) + val*Db(a3,a2)
          endif
          if (n12) then
             Ka(a3,a2) = Ka(a3,a2) + val*Da(a4,a1)
             Kb(a3,a2) = Kb(a3,a2) + val*Db(a4,a1)
          endif
          if (n12 .and. n34) then
             Ka(a4,a2) = Ka(a4,a2) + val*Da(a3,a1)
             Kb(a4,a2) = Kb(a4,a2) + val*Db(a3,a1)
          endif
       endif
    endif
end subroutine accumulate_quartet_jk

end submodule coulomb_impl
