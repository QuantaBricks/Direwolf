! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Exact (non-fitted) 4-center exchange matrix build, including the long-range (range-separated) variant.

submodule (mod_integrals) exchange_impl
implicit none
contains

module subroutine accumulate_quartet_K(e1,e2,e3,e4,val,nConts,D,K)
    implicit none
    integer,intent(in) :: e1,e2,e3,e4,nConts
    real(8),intent(in) :: val
    real(8),intent(in) :: D(nConts,nConts)
    real(8),intent(inout) :: K(nConts,nConts)
    logical :: s12,s34,sbk
    integer :: a1,a2,a3,a4

    a1 = e1+1; a2 = e2+1; a3 = e3+1; a4 = e4+1
    s12 = (e1 .eq. e2)
    s34 = (e3 .eq. e4)
    sbk = (e1 .eq. e3) .and. (e2 .eq. e4)

    K(a1,a3) = K(a1,a3) + val*D(a2,a4)
    if (.not. s12) K(a2,a3) = K(a2,a3) + val*D(a1,a4)
    if (.not. s34) K(a1,a4) = K(a1,a4) + val*D(a2,a3)
    if ((.not. s12) .and. (.not. s34)) K(a2,a4) = K(a2,a4) + val*D(a1,a3)
    if (.not. sbk) then
       K(a3,a1) = K(a3,a1) + val*D(a4,a2)
       if (.not. s34) K(a4,a1) = K(a4,a1) + val*D(a3,a2)
       if (.not. s12) K(a3,a2) = K(a3,a2) + val*D(a4,a1)
       if ((.not. s12) .and. (.not. s34)) K(a4,a2) = K(a4,a2) + val*D(a3,a1)
    endif
end subroutine accumulate_quartet_K

module subroutine integrals_build_exchange(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    if (direct_mode) then
       call build_exchange_direct(nConts, Da, Db, Ka, Kb)
    else
       call build_exchange_store(nConts, Da, Db, Ka, Kb)
    endif
end subroutine integrals_build_exchange

module subroutine integrals_build_exchange_lr(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    call build_exchange_lr_direct(nConts, Da, Db, Ka, Kb)
end subroutine integrals_build_exchange_lr

module subroutine build_exchange_lr_direct(nConts, Da, Db, Ka, Kb)
    use mod_exchange, only: RS_omega
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)

    integer :: ijsh,klsh,si,sj,sk,sl,di,dj,dk,dl
    integer :: is,js,ks,ls,num,ii
    integer :: shls(4)
    real(8),allocatable :: buf2e(:,:,:,:)
    real(8),allocatable :: Ka_local(:,:), Kb_local(:,:)
    integer :: p,q,r,s,e1,e2,e3,e4
    real(8) :: val
    real(8),allocatable :: dmax_shell(:,:)
    integer,allocatable :: ao_offset(:),pq_min_shell(:),pq_max_shell(:)
    integer :: e1t,e2t,pqv

    call ensure_int2e_optimizer_lr()

    Ka = 0
    Kb = 0
    allocate(dmax_shell(0:nBases-1,0:nBases-1))
    block
       real(8),allocatable :: Dsum(:,:)
       allocate(Dsum(nConts,nConts))
       Dsum = abs(Da) + abs(Db)
       call compute_shell_density_bound(nConts, Dsum, dmax_shell)
       deallocate(Dsum)
    end block
    call build_significant_pairs()

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

    env(9) = RS_omega

    !$omp parallel default(shared) private(si,sj,di,dj,klsh,sk,sl,dk,dl,shls,buf2e) &
    !$omp&   private(is,js,ks,ls,p,q,r,s,e1,e2,e3,e4,val) &
    !$omp&   private(Ka_local,Kb_local)
    allocate(Ka_local(nConts,nConts), Kb_local(nConts,nConts))
    Ka_local = 0.0d0
    Kb_local = 0.0d0
    !$omp do schedule(static,1)   ! see build_coulomb_direct: ordered merge alone is not enough
    do ijsh = 0, n_sig_pairs-1
       si = sig_pair_ishl(ijsh)
       sj = sig_pair_jshl(ijsh)
       di = cgto_engine(si, bas)
       dj = cgto_engine(sj, bas)
       do klsh = 0, ijsh
          sk = sig_pair_ishl(klsh)
          sl = sig_pair_jshl(klsh)
          if (schwarz_bound(si,sj)*schwarz_bound(sk,sl) .lt. SCHWARZ_CUTOFF) cycle
          if (schwarz_bound(si,sj)*schwarz_bound(sk,sl) &
              *max(max(dmax_shell(si,sk),dmax_shell(sj,sl)), &
                   max(dmax_shell(si,sl),dmax_shell(sj,sk))) &
              .lt. density_screen_cutoff) cycle
          dk = cgto_engine(sk, bas)
          dl = cgto_engine(sl, bas)
          shls(1) = si; shls(2) = sj; shls(3) = sk; shls(4) = sl
          allocate(buf2e(di,dj,dk,dl))
          call twoe_engine(buf2e, shls, atm, size(atm,2), bas, nBases, env, int2e_opt_lr)

          is = ao_offset(si); js = ao_offset(sj); ks = ao_offset(sk); ls = ao_offset(sl)

          do p = 1,di
             do q = 1,dj
                do r = 1,dk
                   do s = 1,dl
                      e1 = is+p-1; e2 = js+q-1; e3 = ks+r-1; e4 = ls+s-1
                      val = buf2e(p,q,r,s)*NorVEC(e1+1)*NorVEC(e2+1)*NorVEC(e3+1)*NorVEC(e4+1)
                      if (is_canonical_quartet(e1,e2,e3,e4)) then
                         call accumulate_quartet_K(e1,e2,e3,e4,val,nConts,Da,Ka_local)
                         call accumulate_quartet_K(e1,e2,e3,e4,val,nConts,Db,Kb_local)
                      endif
                      if (klsh .ne. ijsh) then
                         if (is_canonical_quartet(e3,e4,e1,e2)) then
                            call accumulate_quartet_K(e3,e4,e1,e2,val,nConts,Da,Ka_local)
                            call accumulate_quartet_K(e3,e4,e1,e2,val,nConts,Db,Kb_local)
                         endif
                      endif
                   enddo
                enddo
             enddo
          enddo
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
             Ka = Ka + Ka_local
             Kb = Kb + Kb_local
          endif
       enddo
    end block
    deallocate(Ka_local,Kb_local)
    !$omp end parallel

    env(9) = 0.0d0

    deallocate(dmax_shell)
    deallocate(ao_offset,pq_min_shell,pq_max_shell)
end subroutine build_exchange_lr_direct

module subroutine build_exchange_pk_list(nConts)
    implicit none
    integer,intent(in) :: nConts
    integer(8) :: t
    integer :: p,q,r,s,pq,rs,n_out,idx_out(4,8),m

    if (.not. coulomb_list_built) then
       call build_coulomb_pk_list(nConts)
       coulomb_list_built = .true.
    endif

    cK_n = 0
    do t = 1,cPK_n
       pq = cPK_pq(t); rs = cPK_rs(t)
       call unpack_pair(pq,p,q)
       call unpack_pair(rs,r,s)
       call expand_quartet_tuples(p-1,q-1,r-1,s-1,n_out,idx_out)
       do m = 1,n_out
          if (idx_out(1,m) .ge. idx_out(3,m)) cK_n = cK_n + 1
       enddo
    enddo

    allocate(cK_i(cK_n),cK_j(cK_n),cK_k(cK_n),cK_l(cK_n),cK_val(cK_n))

    cK_n = 0
    do t = 1,cPK_n
       pq = cPK_pq(t); rs = cPK_rs(t)
       call unpack_pair(pq,p,q)
       call unpack_pair(rs,r,s)
       call expand_quartet_tuples(p-1,q-1,r-1,s-1,n_out,idx_out)
       do m = 1,n_out
          if (idx_out(1,m) .lt. idx_out(3,m)) cycle
          cK_n = cK_n + 1
          cK_i(cK_n) = idx_out(1,m)+1
          cK_j(cK_n) = idx_out(3,m)+1
          cK_k(cK_n) = idx_out(2,m)+1
          cK_l(cK_n) = idx_out(4,m)+1
          cK_val(cK_n) = cPK_val(t)
       enddo
    enddo
end subroutine build_exchange_pk_list

module subroutine unpack_pair(pq,p,q)
    implicit none
    integer,intent(in) :: pq
    integer,intent(out) :: p,q
    p = int(0.5d0*(1.0d0+sqrt(8.0d0*dble(pq)-7.0d0)))
    do while (p*(p-1)/2 .ge. pq)
       p = p - 1
    enddo
    do while ((p+1)*p/2 .lt. pq)
       p = p + 1
    enddo
    q = pq - p*(p-1)/2
end subroutine unpack_pair

module subroutine build_exchange_store(nConts, Da, Db, Ka, Kb)
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    integer(8) :: t
    real(8),allocatable :: Ka_local(:,:), Kb_local(:,:)

    if (.not. exchange_list_built) then
       call build_exchange_pk_list(nConts)
       exchange_list_built = .true.
    endif

    Ka = 0
    Kb = 0
    !$omp parallel default(shared) private(Ka_local,Kb_local)
    allocate(Ka_local(nConts,nConts), Kb_local(nConts,nConts))
    Ka_local = 0.0d0
    Kb_local = 0.0d0
    !$omp do schedule(static)
    do t = 1,cK_n
       Ka_local(cK_i(t),cK_j(t)) = Ka_local(cK_i(t),cK_j(t)) + cK_val(t)*Da(cK_k(t),cK_l(t))
       Kb_local(cK_i(t),cK_j(t)) = Kb_local(cK_i(t),cK_j(t)) + cK_val(t)*Db(cK_k(t),cK_l(t))
       if (cK_i(t) .ne. cK_j(t)) then
          Ka_local(cK_j(t),cK_i(t)) = Ka_local(cK_j(t),cK_i(t)) + cK_val(t)*Da(cK_k(t),cK_l(t))
          Kb_local(cK_j(t),cK_i(t)) = Kb_local(cK_j(t),cK_i(t)) + cK_val(t)*Db(cK_k(t),cK_l(t))
       endif
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
          endif
       enddo
    end block
    deallocate(Ka_local,Kb_local)
    !$omp end parallel
end subroutine build_exchange_store

module subroutine build_exchange_direct(nConts, Da, Db, Ka, Kb)
    use omp_lib, only: omp_get_thread_num, omp_get_num_threads
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)

    integer :: ijsh,klsh,si,sj,sk,sl,di,dj,dk,dl
    integer :: is,js,ks,ls,num,ii
    integer :: shls(4)
    real(8),allocatable :: buf2e(:,:,:,:)
    real(8),allocatable :: Ka_local(:,:), Kb_local(:,:)
    integer :: p,q,r,s,e1,e2,e3,e4
    integer :: a1,a2,a3,a4
    real(8) :: val
    real(8),allocatable :: dmax_shell(:,:)
    integer,allocatable :: ao_offset(:),pq_min_shell(:),pq_max_shell(:)
    integer :: e1t,e2t,pqv

    Ka = 0
    Kb = 0
    allocate(dmax_shell(0:nBases-1,0:nBases-1))
    block
       real(8),allocatable :: Dsum(:,:)
       allocate(Dsum(nConts,nConts))
       Dsum = abs(Da) + abs(Db)
       call compute_shell_density_bound(nConts, Dsum, dmax_shell)
       deallocate(Dsum)
    end block
    call build_significant_pairs()

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
    !$omp&   private(is,js,ks,ls,p,q,r,s,e1,e2,e3,e4,val,a1,a2,a3,a4) &
    !$omp&   private(Ka_local,Kb_local)
    allocate(Ka_local(nConts,nConts), Kb_local(nConts,nConts))
    Ka_local = 0.0d0
    Kb_local = 0.0d0
    !$omp do schedule(static,1)   ! see build_coulomb_direct: ordered merge alone is not enough
    do ijsh = 0, n_sig_pairs-1
       si = sig_pair_ishl(ijsh)
       sj = sig_pair_jshl(ijsh)
       di = cgto_engine(si, bas)
       dj = cgto_engine(sj, bas)
       do klsh = 0, ijsh
          sk = sig_pair_ishl(klsh)
          sl = sig_pair_jshl(klsh)
          if (schwarz_bound(si,sj)*schwarz_bound(sk,sl) .lt. SCHWARZ_CUTOFF) cycle
          if (schwarz_bound(si,sj)*schwarz_bound(sk,sl) &
              *max(max(dmax_shell(si,sk),dmax_shell(sj,sl)), &
                   max(dmax_shell(si,sl),dmax_shell(sj,sk))) &
              .lt. density_screen_cutoff) cycle
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
                         Ka_local(a1,a3) = Ka_local(a1,a3) + val*Da(a2,a4)
                         Ka_local(a2,a3) = Ka_local(a2,a3) + val*Da(a1,a4)
                         Ka_local(a1,a4) = Ka_local(a1,a4) + val*Da(a2,a3)
                         Ka_local(a2,a4) = Ka_local(a2,a4) + val*Da(a1,a3)
                         Ka_local(a3,a1) = Ka_local(a3,a1) + val*Da(a4,a2)
                         Ka_local(a4,a1) = Ka_local(a4,a1) + val*Da(a3,a2)
                         Ka_local(a3,a2) = Ka_local(a3,a2) + val*Da(a4,a1)
                         Ka_local(a4,a2) = Ka_local(a4,a2) + val*Da(a3,a1)
                         Kb_local(a1,a3) = Kb_local(a1,a3) + val*Db(a2,a4)
                         Kb_local(a2,a3) = Kb_local(a2,a3) + val*Db(a1,a4)
                         Kb_local(a1,a4) = Kb_local(a1,a4) + val*Db(a2,a3)
                         Kb_local(a2,a4) = Kb_local(a2,a4) + val*Db(a1,a3)
                         Kb_local(a3,a1) = Kb_local(a3,a1) + val*Db(a4,a2)
                         Kb_local(a4,a1) = Kb_local(a4,a1) + val*Db(a3,a2)
                         Kb_local(a3,a2) = Kb_local(a3,a2) + val*Db(a4,a1)
                         Kb_local(a4,a2) = Kb_local(a4,a2) + val*Db(a3,a1)
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
                            call accumulate_quartet_K(e1,e2,e3,e4,val,nConts,Da,Ka_local)
                            call accumulate_quartet_K(e1,e2,e3,e4,val,nConts,Db,Kb_local)
                         endif
                         if (klsh .ne. ijsh) then
                            if (is_canonical_quartet(e3,e4,e1,e2)) then
                               call accumulate_quartet_K(e3,e4,e1,e2,val,nConts,Da,Ka_local)
                               call accumulate_quartet_K(e3,e4,e1,e2,val,nConts,Db,Kb_local)
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
                         call accumulate_quartet_K(e1,e2,e3,e4,val,nConts,Da,Ka_local)
                         call accumulate_quartet_K(e1,e2,e3,e4,val,nConts,Db,Kb_local)
                      endif
                      if (klsh .ne. ijsh) then
                         if (is_canonical_quartet(e3,e4,e1,e2)) then
                            call accumulate_quartet_K(e3,e4,e1,e2,val,nConts,Da,Ka_local)
                            call accumulate_quartet_K(e3,e4,e1,e2,val,nConts,Db,Kb_local)
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
             Ka = Ka + Ka_local
             Kb = Kb + Kb_local
          endif
       enddo
    end block
    deallocate(Ka_local,Kb_local)
    !$omp end parallel
    deallocate(dmax_shell)
    deallocate(ao_offset,pq_min_shell,pq_max_shell)
end subroutine build_exchange_direct

end submodule exchange_impl
