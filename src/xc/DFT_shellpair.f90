! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Shell-pair-driven XC build (ENGINE_XC_SHELLPAIR=1).

subroutine DFT_calc_shellpair(Exc,Fxc_a,Fxc_b,info)
use MOL_info
use GRID_INFO
use mod_xc, only: xc_eval
    implicit none
INCLUDE 'parameter.h'
    integer,parameter :: BATCH = DFT_BATCH
    integer,parameter :: XC_GEMM_MIN_K = 32
    integer,parameter :: ROWGRP = 48
    real(8),parameter :: SIG_CUT = 1.0d-7
    integer :: info
    real(8) :: Exc
    real(8) :: Fxc_a(nconts,nconts),Fxc_b(nconts,nconts)

    integer :: A,B,ka,ib,nb,bstart,ip,ii,jj,d,igrid,n_sig,maxnd,offA,offB
    integer :: ndA,ndB,ao0A,ao0B,iAt,iBt,nat_b
    integer,allocatable :: at_off(:),at_nd(:),at_atom(:)
    integer,allocatable :: rg_off(:),rg_nd(:),rg_a0(:),rg_a1(:)
    integer :: nrg,irg,k0,k1,nkeep0
    logical :: do_incr,any_live
    real(8),allocatable :: at_phimx(:)
    integer :: ndone_loc,nskip_loc,n_keep
    integer,allocatable :: keep(:)
    real(8),allocatable :: phimx(:),Pk(:,:),v0k(:,:),Wfull(:,:)
    real(8),save :: xc_pscreen = 0.0d0
    integer,save :: npair_done = 0, npair_skip = 0
    character(len=16) :: ps_env
    integer :: ps_stat
    real(8) :: pv,acc,accd(3),w,Fval
    real(8) :: TempB(3)
    real(8),allocatable :: rho_g(:),drv_g(:,:)
    real(8),allocatable :: phiA(:,:),phiB(:,:),d1A(:,:,:),d1B(:,:,:),coorb(:,:)
    integer,allocatable :: sig_idx(:)
    real(8),allocatable :: val0(:,:),val1(:,:,:)
    real(8),allocatable :: rho_b_(:),drv_b_(:,:),sig_b(:)
    real(8),allocatable :: rho_r_(:),drv_r_(:,:),sig_r(:)
    real(8),allocatable :: ex_r(:),ec_r(:),d1r_r(:),d1s_r(:)
    real(8),allocatable :: ex_b(:),ec_b(:),d1r(:),d1s(:)
    real(8),allocatable :: si1(:),si2(:),si3(:),vo1(:),vo2(:),vo3(:),vo4(:),vo5(:)
    real(8),allocatable :: WeightB(:),Ta(:,:),ValW(:,:),Fxc_s(:,:),Fxc_local(:,:)
    real(8),allocatable :: PAblk(:,:),WAblk(:,:),Ps(:,:)
    real(8) :: Exc_local

    call GTOeval_build_shell_tables()
    ps_env = ""
    call get_environment_variable("ENGINE_XC_PSCREEN", ps_env)
    if (len_trim(ps_env) .gt. 0) read(ps_env,*,iostat=ps_stat) xc_pscreen
    npair_done = 0
    npair_skip = 0

    ps_env = ""
    call get_environment_variable("ENGINE_XC_INCR", ps_env)
    xc_incr_on = (trim(ps_env) .eq. "1")
    if (xc_incr_on) then
       if (.not. allocated(xcr_Pa)) then
          allocate(xcr_Pa(nconts,nconts),xcr_rho(ngrids),xcr_grad(3,ngrids))
          allocate(xcr_Fxc(nconts,nconts),atpair_ddp(natoms,natoms))
       endif
    endif
    if (.not. allocated(atpair_dp)) allocate(atpair_dp(natoms,natoms))
    atpair_dp = 0.0d0
    do jj = 1,nconts
       do ii = 1,jj
          pv = abs(Pa(ii,jj))
          if (pv .gt. atpair_dp(ao_atom(ii),ao_atom(jj))) then
             atpair_dp(ao_atom(ii),ao_atom(jj)) = pv
             atpair_dp(ao_atom(jj),ao_atom(ii)) = pv
          endif
       enddo
    enddo

    do_incr = xc_incr_on .and. xc_incr_primed
    if (do_incr) then
       atpair_ddp = 0.0d0
       do jj = 1,nconts
          do ii = 1,jj
             pv = abs(Pa(ii,jj)-xcr_Pa(ii,jj))
             if (pv .gt. atpair_ddp(ao_atom(ii),ao_atom(jj))) then
                atpair_ddp(ao_atom(ii),ao_atom(jj)) = pv
                atpair_ddp(ao_atom(jj),ao_atom(ii)) = pv
             endif
          enddo
       enddo
    endif

    Fxc_a = 0
    Fxc_b = 0
    Exc = 0
    maxnd = max(max_dft_batch_nsig,1)

    !$omp parallel &
    !$omp&  private(ib,bstart,nb,ip,ii,jj,d,igrid,n_sig,w,TempB,Fval) &
    !$omp&  private(A,B,ndA,ndB,ao0A,ao0B,offA,offB,pv,acc,accd,PAblk,WAblk,Ps) &
    !$omp&  private(iAt,iBt,nat_b,at_off,at_nd,at_atom,at_phimx,ndone_loc,nskip_loc) &
    !$omp&  private(nrg,irg,rg_off,rg_nd,rg_a0,rg_a1,k0,k1,nkeep0,any_live) &
    !$omp&  private(n_keep,keep,phimx,Pk,v0k,Wfull) &
    !$omp&  private(sig_idx,val0,val1,rho_b_,drv_b_,sig_b,ex_b,ec_b,d1r,d1s) &
    !$omp&  private(rho_r_,drv_r_,sig_r,ex_r,ec_r,d1r_r,d1s_r) &
    !$omp&  private(si1,si2,si3,vo1,vo2,vo3,vo4,vo5,WeightB,Ta,ValW,Fxc_s) &
    !$omp&  private(Fxc_local,Exc_local)
    allocate(sig_idx(max(max_dft_batch_nsig,1)))
    allocate(val0(max(max_dft_batch_nsig,1),BATCH),val1(max(max_dft_batch_nsig,1),3,BATCH))
    allocate(rho_b_(BATCH),drv_b_(BATCH,3),sig_b(BATCH))
    allocate(rho_r_(BATCH),drv_r_(BATCH,3),sig_r(BATCH))
    allocate(ex_r(BATCH),ec_r(BATCH),d1r_r(BATCH),d1s_r(BATCH))
    allocate(ex_b(BATCH),ec_b(BATCH),d1r(BATCH),d1s(BATCH))
    allocate(si1(BATCH),si2(BATCH),si3(BATCH))
    allocate(vo1(BATCH),vo2(BATCH),vo3(BATCH),vo4(BATCH),vo5(BATCH))
    si1 = 0.0d0; si2 = 0.0d0; si3 = 0.0d0
    allocate(WeightB(BATCH),Ta(max(max_dft_batch_nsig,1),BATCH))
    allocate(ValW(max(max_dft_batch_nsig,1),BATCH))
    allocate(Fxc_s(max(max_dft_batch_nsig,1),max(max_dft_batch_nsig,1)))
    allocate(PAblk(maxnd,max(max_dft_batch_nsig,1)),WAblk(maxnd,BATCH))
    allocate(Ps(max(max_dft_batch_nsig,1),max(max_dft_batch_nsig,1)))
    allocate(at_off(max(max_dft_batch_nsig,1)),at_nd(max(max_dft_batch_nsig,1)))
    allocate(at_atom(max(max_dft_batch_nsig,1)),at_phimx(max(max_dft_batch_nsig,1)))
    allocate(rg_off(max(max_dft_batch_nsig,1)),rg_nd(max(max_dft_batch_nsig,1)))
    allocate(rg_a0(max(max_dft_batch_nsig,1)),rg_a1(max(max_dft_batch_nsig,1)))
    allocate(keep(max(max_dft_batch_nsig,1)),phimx(max(max_dft_batch_nsig,1)))
    allocate(Pk(max(max_dft_batch_nsig,1),max(max_dft_batch_nsig,1)))
    allocate(v0k(max(max_dft_batch_nsig,1),BATCH))
    allocate(Wfull(max(max_dft_batch_nsig,1),max(max_dft_batch_nsig,BATCH)))
    allocate(Fxc_local(nconts,nconts))
    Fxc_local = 0.0d0
    Exc_local = 0.0d0
    ndone_loc = 0
    nskip_loc = 0

    !$omp do schedule(dynamic,16)
    do ib = 1,n_dft_batches
       n_sig = dft_batch_shells(ib)%n_sig
       if (n_sig .eq. 0) cycle
       bstart = (ib-1)*BATCH + 1
       nb = min(BATCH, ngrids-bstart+1)
       sig_idx(1:n_sig) = dft_batch_shells(ib)%sig_idx(1:n_sig)
       if (do_incr) then
          any_live = .false.
          do iAt = 1,batch_natom(ib)
             do iBt = 1,batch_natom(ib)
                if (atpair_ddp(batch_atoms(iAt,ib),batch_atoms(iBt,ib)) .ge. xc_pscreen) then
                   any_live = .true.
                   exit
                endif
             enddo
             if (any_live) exit
          enddo
          if (.not. any_live) then
             nskip_loc = nskip_loc + 1
             cycle
          endif
          ndone_loc = ndone_loc + 1
       endif
       do ip = 1,nb
          WeightB(ip) = Grids(bstart+ip-1)%weight
       enddo
       call GTOeval_batch_compact(ib,bstart,nb,max_dft_batch_nsig,val0,val1)

       do ii = 1,n_sig
          phimx(ii) = 0.0d0
          do ip = 1,nb
             phimx(ii) = max(phimx(ii),abs(val0(ii,ip)))
          enddo
       enddo
       n_keep = 0
       do ii = 1,n_sig
          if (phimx(ii) .gt. SIG_CUT) then
             n_keep = n_keep + 1
             keep(n_keep) = ii
          endif
       enddo
       if (n_keep .lt. n_sig) then
          do ip = 1,nb
             do ii = 1,n_keep
                v0k(ii,ip) = val0(keep(ii),ip)
             enddo
          enddo
          do ip = 1,nb
             do ii = 1,n_keep
                val0(ii,ip) = v0k(ii,ip)
             enddo
          enddo
          do d = 1,3
             do ip = 1,nb
                do ii = 1,n_keep
                   v0k(ii,ip) = val1(keep(ii),d,ip)
                enddo
             enddo
             do ip = 1,nb
                do ii = 1,n_keep
                   val1(ii,d,ip) = v0k(ii,ip)
                enddo
             enddo
          enddo
          do ii = 1,n_keep
             sig_idx(ii) = sig_idx(keep(ii))
             phimx(ii) = phimx(keep(ii))
          enddo
          n_sig = n_keep
       endif
       if (n_sig .eq. 0) cycle

       rho_b_(1:nb) = 0.0d0
       drv_b_(1:nb,1) = 0.0d0
       drv_b_(1:nb,2) = 0.0d0
       drv_b_(1:nb,3) = 0.0d0
       if (do_incr) then
          do jj = 1,n_sig
             do ii = 1,jj
                Ps(ii,jj) = Pa(sig_idx(ii),sig_idx(jj)) - xcr_Pa(sig_idx(ii),sig_idx(jj))
                if (ii .ne. jj) Ps(jj,ii) = Ps(ii,jj)
             enddo
          enddo
       else
          do jj = 1,n_sig
             do ii = 1,jj
                Ps(ii,jj) = Pa(sig_idx(ii),sig_idx(jj))
                if (ii .ne. jj) Ps(jj,ii) = Ps(ii,jj)
             enddo
          enddo
       endif
       if (do_incr) then
          do jj = 1,n_sig
             do ii = 1,jj
                Ps(ii,jj) = Pa(sig_idx(ii),sig_idx(jj)) - xcr_Pa(sig_idx(ii),sig_idx(jj))
                if (ii .ne. jj) Ps(jj,ii) = Ps(ii,jj)
             enddo
          enddo
       else
          do jj = 1,n_sig
             do ii = 1,jj
                Ps(ii,jj) = Pa(sig_idx(ii),sig_idx(jj))
                if (ii .ne. jj) Ps(jj,ii) = Ps(ii,jj)
             enddo
          enddo
       endif

       nat_b = 0
       ii = 1
       do while (ii .le. n_sig)
          nat_b = nat_b + 1
          at_off(nat_b)  = ii - 1
          at_atom(nat_b) = ao_atom(sig_idx(ii))
          at_nd(nat_b)   = 0
          do while (ii .le. n_sig)
             if (ao_atom(sig_idx(ii)) .ne. at_atom(nat_b)) exit
             at_nd(nat_b) = at_nd(nat_b) + 1
             ii = ii + 1
          enddo
       enddo
       do iAt = 1,nat_b
          at_phimx(iAt) = 0.0d0
          do ii = 1,at_nd(iAt)
             at_phimx(iAt) = max(at_phimx(iAt), phimx(at_off(iAt)+ii))
          enddo
       enddo

       rho_b_(1:nb) = 0.0d0
       drv_b_(1:nb,1) = 0.0d0
       drv_b_(1:nb,2) = 0.0d0
       drv_b_(1:nb,3) = 0.0d0
       nrg = 0
       iAt = 1
       do while (iAt .le. nat_b)
          nrg = nrg + 1
          rg_off(nrg) = at_off(iAt)
          rg_nd(nrg)  = 0
          rg_a0(nrg)  = iAt
          do while (iAt .le. nat_b)
             rg_nd(nrg) = rg_nd(nrg) + at_nd(iAt)
             iAt = iAt + 1
             if (rg_nd(nrg) .ge. ROWGRP) exit
          enddo
          rg_a1(nrg) = iAt - 1
       enddo

       do irg = 1,nrg
          ndA  = rg_nd(irg)
          offA = rg_off(irg)
          if (ndA .eq. 0) cycle
          n_keep = 0
          nkeep0 = 0
          do iBt = 1,nat_b
             ndB = at_nd(iBt)
             if (ndB .eq. 0) cycle
             offB = at_off(iBt)
             pv = 0.0d0
             do iAt = rg_a0(irg),rg_a1(irg)
                if (do_incr) then
                   pv = max(pv, atpair_ddp(at_atom(iAt),at_atom(iBt))*at_phimx(iAt))
                else
                   pv = max(pv, atpair_dp(at_atom(iAt),at_atom(iBt))*at_phimx(iAt))
                endif
             enddo
             if (pv*at_phimx(iBt) .lt. xc_pscreen) then
                nskip_loc = nskip_loc + 1
                cycle
             endif
             ndone_loc = ndone_loc + 1
             do jj = 1,ndB
                n_keep = n_keep + 1
                keep(n_keep) = offB + jj
             enddo
          enddo
          if (n_keep .eq. 0) then
             do ip = 1,nb
                do ii = 1,ndA
                   Wfull(offA+ii,ip) = 0.0d0
                enddo
             enddo
             cycle
          endif
          k0 = 1
          do while (k0 .le. n_keep)
             k1 = k0
             do while (k1 .lt. n_keep)
                if (keep(k1+1) .ne. keep(k1)+1) exit
                k1 = k1 + 1
             enddo
             call dgemm('N','N',ndA,nb,k1-k0+1,1.0d0,Ps(offA+1,keep(k0)),max_dft_batch_nsig, &
                        val0(keep(k0),1),max_dft_batch_nsig, &
                        merge(0.0d0,1.0d0,k0 .eq. 1),Wfull(offA+1,1),max_dft_batch_nsig)
             k0 = k1 + 1
          enddo
       enddo
       do ip = 1,nb
          acc = 0.0d0
          accd = 0.0d0
          do ii = 1,n_sig
             acc = acc + val0(ii,ip)*Wfull(ii,ip)
             accd(1) = accd(1) + val1(ii,1,ip)*Wfull(ii,ip)
             accd(2) = accd(2) + val1(ii,2,ip)*Wfull(ii,ip)
             accd(3) = accd(3) + val1(ii,3,ip)*Wfull(ii,ip)
          enddo
          rho_b_(ip) = acc
          drv_b_(ip,1) = 2.0d0*accd(1)
          drv_b_(ip,2) = 2.0d0*accd(2)
          drv_b_(ip,3) = 2.0d0*accd(3)
       enddo
       if (do_incr) then
          do ip = 1,nb
             igrid = bstart+ip-1
             rho_r_(ip)   = xcr_rho(igrid)
             drv_r_(ip,1) = xcr_grad(1,igrid)
             drv_r_(ip,2) = xcr_grad(2,igrid)
             drv_r_(ip,3) = xcr_grad(3,igrid)
             rho_b_(ip)   = rho_r_(ip)   + rho_b_(ip)
             drv_b_(ip,1) = drv_r_(ip,1) + drv_b_(ip,1)
             drv_b_(ip,2) = drv_r_(ip,2) + drv_b_(ip,2)
             drv_b_(ip,3) = drv_r_(ip,3) + drv_b_(ip,3)
             xcr_rho(igrid)    = rho_b_(ip)
             xcr_grad(1,igrid) = drv_b_(ip,1)
             xcr_grad(2,igrid) = drv_b_(ip,2)
             xcr_grad(3,igrid) = drv_b_(ip,3)
             sig_r(ip) = (2*drv_r_(ip,1))**2 + (2*drv_r_(ip,2))**2 + (2*drv_r_(ip,3))**2
          enddo
       endif
       if (xc_incr_on .and. .not. do_incr) then
          do ip = 1,nb
             igrid = bstart+ip-1
             xcr_rho(igrid)    = rho_b_(ip)
             xcr_grad(1,igrid) = drv_b_(ip,1)
             xcr_grad(2,igrid) = drv_b_(ip,2)
             xcr_grad(3,igrid) = drv_b_(ip,3)
          enddo
       endif
       do ip = 1,nb
          sig_b(ip) = (2*drv_b_(ip,1))**2 + (2*drv_b_(ip,2))**2 + (2*drv_b_(ip,3))**2
       enddo

       call xc_eval(Multi,1,nb,rho_b_(1:nb),rho_b_(1:nb),si1(1:nb),si2(1:nb),si3(1:nb), &
                    sig_b(1:nb),ex_b(1:nb),ec_b(1:nb),d1r(1:nb),d1s(1:nb), &
                    vo1(1:nb),vo2(1:nb),vo3(1:nb),vo4(1:nb),vo5(1:nb))
       if (do_incr) then
          call xc_eval(Multi,1,nb,rho_r_(1:nb),rho_r_(1:nb),si1(1:nb),si2(1:nb),si3(1:nb), &
                       sig_r(1:nb),ex_r(1:nb),ec_r(1:nb),d1r_r(1:nb),d1s_r(1:nb), &
                       vo1(1:nb),vo2(1:nb),vo3(1:nb),vo4(1:nb),vo5(1:nb))
       endif
       do ip = 1,nb
          w = WeightB(ip)
          if (do_incr) then
             Exc_local = Exc_local + ((ex_b(ip)+ec_b(ip))-(ex_r(ip)+ec_r(ip)))*w
             TempB(1) = 4*d1s(ip)*drv_b_(ip,1) - 4*d1s_r(ip)*drv_r_(ip,1)
             TempB(2) = 4*d1s(ip)*drv_b_(ip,2) - 4*d1s_r(ip)*drv_r_(ip,2)
             TempB(3) = 4*d1s(ip)*drv_b_(ip,3) - 4*d1s_r(ip)*drv_r_(ip,3)
             d1r(ip) = d1r(ip) - d1r_r(ip)
             Ta(1:n_sig,ip) = val1(1:n_sig,1,ip)*TempB(1) + val1(1:n_sig,2,ip)*TempB(2) &
                            + val1(1:n_sig,3,ip)*TempB(3)
             ValW(1:n_sig,ip) = 0.5d0*val0(1:n_sig,ip)*(d1r(ip)*w) + Ta(1:n_sig,ip)*w
             cycle
          endif
          Exc_local = Exc_local + (ex_b(ip)+ec_b(ip))*w
          TempB(1) = 4*d1s(ip)*drv_b_(ip,1)
          TempB(2) = 4*d1s(ip)*drv_b_(ip,2)
          TempB(3) = 4*d1s(ip)*drv_b_(ip,3)
          Ta(1:n_sig,ip) = val1(1:n_sig,1,ip)*TempB(1) + val1(1:n_sig,2,ip)*TempB(2) &
                         + val1(1:n_sig,3,ip)*TempB(3)
          ValW(1:n_sig,ip) = 0.5d0*val0(1:n_sig,ip)*(d1r(ip)*w) + Ta(1:n_sig,ip)*w
       enddo
       call dgemm('N','T',n_sig,n_sig,nb,1.0d0,ValW,max_dft_batch_nsig, &
                  val0,max_dft_batch_nsig,0.0d0,Fxc_s,max_dft_batch_nsig)
       do jj = 1,n_sig
          do ii = 1,jj
             Fval = Fxc_s(ii,jj) + Fxc_s(jj,ii)
             Fxc_local(sig_idx(ii),sig_idx(jj)) = Fxc_local(sig_idx(ii),sig_idx(jj)) + Fval
             if (ii .ne. jj) Fxc_local(sig_idx(jj),sig_idx(ii)) = Fxc_local(sig_idx(jj),sig_idx(ii)) + Fval
          enddo
       enddo
    enddo
    !$omp end do

    !$omp critical
    npair_done = npair_done + ndone_loc
    npair_skip = npair_skip + nskip_loc
    Exc = Exc + Exc_local
    Fxc_a = Fxc_a + Fxc_local
    !$omp end critical
    deallocate(sig_idx,val0,val1,rho_b_,drv_b_,sig_b,ex_b,ec_b,d1r,d1s)
    deallocate(rho_r_,drv_r_,sig_r,ex_r,ec_r,d1r_r,d1s_r)
    deallocate(si1,si2,si3,vo1,vo2,vo3,vo4,vo5,WeightB,Ta,ValW,Fxc_s)
    deallocate(Fxc_local,PAblk,WAblk,Ps,at_off,at_nd,at_atom,at_phimx)
    deallocate(rg_off,rg_nd,rg_a0,rg_a1)
    deallocate(keep,phimx,Pk,v0k,Wfull)
    !$omp end parallel

    if (xc_pscreen .gt. 0.0d0) print '(" XC atom pairs: ",I0," done, ",I0," skipped (",F5.1,"%)")', &
          npair_done, npair_skip, 100.0*real(npair_skip)/real(max(npair_done+npair_skip,1))
    if (xc_incr_on) then
       if (do_incr) then
          Fxc_a = xcr_Fxc + Fxc_a
          Exc   = xcr_Exc + Exc
       endif
       xcr_Fxc = Fxc_a
       xcr_Exc = Exc
       xcr_Pa  = Pa
       xc_incr_primed = .true.
    endif
    Fxc_b = Fxc_a
end subroutine
