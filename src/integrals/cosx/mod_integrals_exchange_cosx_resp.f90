! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX response K for the CPSCF (many trial densities at once), split out

submodule (mod_integrals:exchange_cosx_impl) exchange_cosx_resp_impl
real(8),save :: cosx_batch_omega = 0.0d0
contains
module subroutine integrals_build_exchange_cosx_batch(nConts, nd, D, K)
use mod_exchange, only: RS_omega, RS_beta, HF_exchange_frac
integer,intent(in) :: nConts, nd
real(8),intent(in) :: D(nConts,nConts,nd)
real(8),intent(out) :: K(nConts,nConts,nd)
real(8) :: dum(1,1,1)
real(8),allocatable :: Klr(:,:,:)
call cosx_batch_core(nConts, nd, D, K, nd, 1, dum, dum, 0, dum)
if (RS_omega .gt. 0.0d0 .and. HF_exchange_frac .gt. 0.0d0) then
   allocate(Klr(nConts,nConts,nd))
   cosx_batch_omega = RS_omega
   call cosx_batch_core(nConts, nd, D, Klr, nd, 1, dum, dum, 0, dum)
   cosx_batch_omega = 0.0d0
   K(:,:,1:nd) = K(:,:,1:nd) + (RS_beta/HF_exchange_frac)*Klr
   deallocate(Klr)
endif
end subroutine integrals_build_exchange_cosx_batch

module subroutine integrals_build_exchange_cosx_batch_c(nConts, nd, Lf, nop, Co, KC)
use mod_exchange, only: RS_omega, RS_beta, HF_exchange_frac
integer,intent(in) :: nConts, nd, nop
real(8),intent(in) :: Lf(nConts,nop,nd), Co(nConts,nop)
real(8),intent(out) :: KC(nConts,nop,nd)
real(8) :: dum(1,1,1)
real(8),allocatable :: KClr(:,:,:)
call cosx_batch_core(nConts, nd, dum, dum, 0, nop, Co, KC, nd, Lf)
if (RS_omega .gt. 0.0d0 .and. HF_exchange_frac .gt. 0.0d0) then
   allocate(KClr(nConts,nop,nd))
   cosx_batch_omega = RS_omega
   call cosx_batch_core(nConts, nd, dum, dum, 0, nop, Co, KClr, nd, Lf)
   cosx_batch_omega = 0.0d0
   KC(:,:,1:nd) = KC(:,:,1:nd) + (RS_beta/HF_exchange_frac)*KClr
   deallocate(KClr)
endif
end subroutine integrals_build_exchange_cosx_batch_c

subroutine cosx_batch_core(nConts, nd, D, K, nkd, nop, Co, KC, nkc, Lf)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_wtime, &
                   omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, omp_destroy_lock
use mod_fixq, only: to_q, from_q
implicit none
integer,intent(in) :: nConts, nd
integer,intent(in) :: nkd, nop, nkc
real(8),intent(in) :: D(nConts,nConts,*), Co(nConts,*), Lf(nConts,nop,*)
real(8),intent(inout) :: K(nConts,nConts,*), KC(nConts,nop,*)
logical :: proj
real(8),allocatable :: ACo(:,:), Hc(:,:,:), Lt(:,:), CoT(:,:), Yb(:,:), Zb(:,:)
integer :: nrow
integer,parameter :: INT_CACHE_SIZE = 4000000
integer :: n, ib, bstart, nb, g, si, sj, sj_idx, di, dj, i, j, ao_i, ao_j, m, kk, cofs
integer :: g0, ng, nbs, ncol, ndc, c0, nc, shls1e(4), off_grids, maxd, nloc
logical :: do_var, fit
real(8) :: kscr_loc
real(8) :: nij, esp_eff, r_far, seg_ab(3), seg_cl(3), seg_len2, seg_t, xbm, k_bound, sg
real(8),allocatable :: Xd(:,:), Qb(:,:), Xc(:,:), vdum(:,:,:), RR(:,:,:), OO(:,:,:), Asub(:,:,:)
real(8),allocatable :: Kloc(:,:), Fgm(:), buf(:), env_g(:), int_cache(:)
integer,allocatable :: run_a(:), run_len(:)
integer :: nrun, ir, a0
character(len=32) :: envs
real(8),pointer :: coor_p(:,:), w_p(:), blkcen_p(:,:), blkrad_p(:)
real(8),pointer,contiguous :: QfullT_p(:,:)
integer,pointer :: blkptr_p(:), blkidx_p(:), blkstart_p(:)
integer :: nblk_p
real(8),allocatable :: Dp(:,:,:), KV(:,:), KVs(:,:), KT(:,:)
integer(16),allocatable :: KTq(:,:)
real(8) :: tb0, tb(12), tw(4), tw0, tbusy, bmax, bmin, tmg
integer,allocatable :: border(:)
real(8),allocatable :: bcost(:)
integer :: ibo
logical :: btime
logical :: rsub, kshr
integer,parameter :: KCOL = 32
integer(omp_lock_kind),allocatable :: klk(:)
real(8),allocatable :: Tk(:,:)
integer :: rg0, cc0, ncc, cend, kshm, nkg, mgs, kg, mg0, mgn
real(8),allocatable :: LfP(:,:,:)
integer :: nrp, rb0, nrg
logical :: coop
integer :: ibS, nbB, bstartS, nlocS, nrunS, ntl, offgS, rb0S, nrgS, rgS, jr, r0, r1, nthr, mytid
integer,allocatable :: shofs(:), runA(:), runL(:), tl0(:), tln(:)
real(8),allocatable :: XdS(:,:), XcS(:,:), vdS(:,:,:), QbS(:,:), AsubS(:,:,:), RRS(:,:,:), HcS(:,:,:)
real(8),allocatable :: YbS(:,:), ZbS(:,:), FgmS(:), envGS(:), ACoP(:,:), bufP(:), icP(:)
real(8) :: xbmS, tc(8), tc0

n = nConts
proj = (nkc .gt. 0)
kscr_loc = cosx_kscreen
if (cosx_resp_follow_scf) then
   block
      real(8) :: dmx, fac
      character(len=16) :: ef
      dmx = maxval(abs(D(:,:,1:nd)))
      fac = 1.0d0
      call get_environment_variable("ENGINE_NEWTON_KFAC", ef)
      if (len_trim(ef) .gt. 0) read(ef,*) fac
      kscr_loc = max(1.0d-14, cosx_kscreen*fac*min(1.0d0, dmx/1.0d-2))
   end block
endif
call get_environment_variable("ENGINE_CPSCF_COSX_GRID", envs)
if (len_trim(envs) .eq. 0) envs = "lo"
if (cosx_resp_follow_scf) then
   if (cosx_force_hi_grid) then
      envs = "hi"
   else if (cosx_promote_md) then
      envs = "md"
   else
      envs = "lo"
   endif
endif
select case (trim(envs))
case ("md")
   coor_p => cosx_coor_md; w_p => cosx_w_md; QfullT_p => cosx_QfullT_md
   blkptr_p => cosx_blkptr_md; blkidx_p => cosx_blkidx_md; blkstart_p => cosx_blkstart_md
   blkcen_p => cosx_blkcen_md; blkrad_p => cosx_blkrad_md; nblk_p = cosx_nblk_md
case ("lo")
   coor_p => cosx_coor_lo; w_p => cosx_w_lo; QfullT_p => cosx_QfullT_lo
   blkptr_p => cosx_blkptr_lo; blkidx_p => cosx_blkidx_lo; blkstart_p => cosx_blkstart_lo
   blkcen_p => cosx_blkcen_lo; blkrad_p => cosx_blkrad_lo; nblk_p = cosx_nblk_lo
case default
   coor_p => cosx_coor_hi; w_p => cosx_w_hi; QfullT_p => cosx_QfullT_hi
   blkptr_p => cosx_blkptr_hi; blkidx_p => cosx_blkidx_hi; blkstart_p => cosx_blkstart_hi
   blkcen_p => cosx_blkcen_hi; blkrad_p => cosx_blkrad_hi; nblk_p = cosx_nblk_hi
end select
fit = cosx_overlap_fit
do_var = fit
call get_environment_variable("ENGINE_COSX_K_NONVAR", envs)
if (trim(envs) .eq. "1") do_var = .false.
call get_environment_variable("ENGINE_CPSCF_COSX_NONVAR", envs)
if (trim(envs) .ne. "0" .and. .not. cosx_resp_follow_scf) do_var = .false.
if (proj .and. do_var) then
   print *, "cosx_batch_core: the C-projected response K is single-leg only"
   stop 1
endif
ncol = merge(2, 1, do_var)
nrow = merge(nop, n, proj)
call get_environment_variable("ENGINE_CPSCF_COSX_RSUB", envs)
rsub = proj .and. trim(envs) .ne. "0"
nrp = 64
call get_environment_variable("ENGINE_CPSCF_COSX_RPTS", envs)
if (len_trim(envs) .gt. 0) read(envs,*) nrp
call get_environment_variable("ENGINE_CPSCF_COSX_KSHARED", envs)
kshr = proj .and. trim(envs) .ne. "0"
call get_environment_variable("ENGINE_CPSCF_COSX_COOP", envs)
coop = proj .and. rsub .and. .not. do_var .and. trim(envs) .eq. "1"
kshm = 1
if (trim(envs) .eq. "2") kshm = 2
mgs = max(1, (nd+15)/16)
nkg = (nd+mgs-1)/mgs
if (kshr) then
   allocate(klk(max((n+KCOL-1)/KCOL, nkg)))
   do kk = 1,size(klk)
      call omp_init_lock(klk(kk))
   enddo
endif
maxd = maxval(cosx_shell_dim(0:nBases-1))
nbs = max(1, min(BATCH, 4000000/(n*n)))
if (coop) nbs = max(1, min(BATCH, 16000000/(n*n)))
call get_environment_variable("ENGINE_CPSCF_COSX_NBS", envs)
if (coop .and. len_trim(envs) .gt. 0) read(envs,*) nbs
nrp = max(nbs, (nrp/nbs)*nbs)
ndc = max(1, min(nd, 32000000/(n*n)))
if (proj) ndc = nd
if (.not. proj) K(:,:,1:nd) = 0.0d0
allocate(border(nblk_p), bcost(nblk_p))
do ib = 0, nblk_p-1
   nloc = 0
   do kk = blkptr_p(ib)+1, blkptr_p(ib+1)
      nloc = nloc + cosx_shell_dim(blkidx_p(kk))
   enddo
   bcost(ib+1) = dble(blkstart_p(ib+1) - blkstart_p(ib))*dble(nloc)
   border(ib+1) = ib
enddo
do ib = 2, nblk_p
   ibo = border(ib)
   tb0 = bcost(ib)
   kk = ib - 1
   do while (kk .ge. 1)
      if (bcost(kk) .ge. tb0) exit
      border(kk+1) = border(kk) ; bcost(kk+1) = bcost(kk)
      kk = kk - 1
   enddo
   border(kk+1) = ibo ; bcost(kk+1) = tb0
enddo
call get_environment_variable("ENGINE_COSX_BATCH_TIME", envs)
btime = (trim(envs) .eq. "1")
tb = 0.0d0 ; tw = 0.0d0 ; tw0 = omp_get_wtime()

do c0 = 0, nd-1, ndc
nc = min(ndc, nd-c0)
if (proj) then
   allocate(Dp(1,1,1), KT(nrow*nc,n), Lt(nop*nc,n), CoT(nop,n))
   CoT = transpose(Co(:,1:nop))
   !$omp parallel do private(kk,m)
   do kk = 1,n
      do m = 1,nc
         Lt(nop*(m-1)+1:nop*m,kk) = Lf(kk,:,c0+m)
      enddo
   enddo
   !$omp end parallel do
   if (rsub) then
      allocate(LfP(n,nc,nop))
      !$omp parallel do private(kk,m)
      do kk = 1,nop
         do m = 1,nc
            LfP(:,m,kk) = Lf(:,kk,c0+m)
         enddo
      enddo
      !$omp end parallel do
   endif
else
   allocate(Dp(n,nc,n), KT(nrow*nc,n))
endif
if (kshr) allocate(KTq(size(KT,1),n))
!$omp parallel do private(kk)
do kk = 1,n
   KT(:,kk) = 0.0d0
   if (kshr) KTq(:,kk) = 0_16
enddo
!$omp end parallel do
tw(1) = tw(1) + (omp_get_wtime()-tw0) ; tw0 = omp_get_wtime()
if (do_var) then
   allocate(KVs(n*nc,n))
   KVs = 0.0d0
endif
if (.not. proj) then
!$omp parallel do private(kk,m)
do kk = 1,n
   do m = 1,nc
      Dp(:,m,kk) = D(:,kk,c0+m)
   enddo
enddo
!$omp end parallel do
endif
tw(1) = tw(1) + (omp_get_wtime()-tw0) ; tw0 = omp_get_wtime()
bmax = 0.0d0 ; bmin = 1.0d30 ; tmg = 0.0d0
if (coop) then
allocate(XdS(n,BATCH), RRS(n,nc,nrp), AsubS(n,n,nbs), HcS(nop,nc,BATCH), YbS(nop,nrp), ZbS(nop*nc,nrp))
allocate(FgmS(0:nBases-1), shofs(nBases), runA(nBases), runL(nBases), tl0(n), tln(n))
allocate(envGS(size(env) + 3*BATCH))
envGS(1:size(env)) = env(1:size(env))
envGS(9) = cosx_batch_omega
offgS = size(env)
if (fit) allocate(QbS(n,BATCH))
!$omp parallel default(shared) &
!$omp&  private(ibo,g,si,sj,sj_idx,di,dj,i,j,ao_i,ao_j,m,kk,g0,ng,shls1e,ir,jr,r0,r1,nthr,mytid) &
!$omp&  private(nij,esp_eff,r_far,seg_ab,seg_cl,seg_len2,seg_t,k_bound,sg,ACoP,bufP,icP)
allocate(bufP(BATCH*maxd*maxd), icP(INT_CACHE_SIZE), ACoP(n,nop))
!$omp master
tc = 0.0d0 ; tc0 = omp_get_wtime()
!$omp end master
nthr = omp_get_num_threads()
mytid = omp_get_thread_num()
do ibo = 1, nblk_p
   !$omp single
   ibS = border(ibo)
   bstartS = blkstart_p(ibS) + 1
   nbB = blkstart_p(ibS+1) - blkstart_p(ibS)
   nlocS = 0
   nrunS = 0
   do kk = blkptr_p(ibS)+1, blkptr_p(ibS+1)
      si = blkidx_p(kk)
      shofs(kk-blkptr_p(ibS)) = nlocS
      nlocS = nlocS + cosx_shell_dim(si)
      if (nrunS .gt. 0) then
         if (cosx_ao_offset(si)+1 .eq. runA(nrunS)+runL(nrunS)) then
            runL(nrunS) = runL(nrunS) + cosx_shell_dim(si)
            cycle
         endif
      endif
      nrunS = nrunS + 1
      runA(nrunS) = cosx_ao_offset(si) + 1
      runL(nrunS) = cosx_shell_dim(si)
   enddo
   ntl = 0
   do ir = 1,nrunS
      cc0 = runA(ir)
      cend = runA(ir) + runL(ir) - 1
      do while (cc0 .le. cend)
         ntl = ntl + 1
         tl0(ntl) = cc0
         tln(ntl) = min(cend, KCOL*((cc0-1)/KCOL+1)) - cc0 + 1
         cc0 = cc0 + tln(ntl)
      enddo
   enddo
   if (nlocS .gt. 0) then
      allocate(XcS(nlocS,nbB), vdS(nlocS,3,nbB))
      XcS = 0.0d0
      do g = 1,nbB
         envGS(offgS+3*(g-1)+1:offgS+3*g) = coor_p(:,bstartS+g-1)
      enddo
      envGS(12) = dble(nbB)
      envGS(13) = dble(offgS)
   endif
   !$omp end single
   if (nlocS .eq. 0) cycle
   !$omp do schedule(dynamic)
   do kk = blkptr_p(ibS)+1, blkptr_p(ibS+1)
      si = blkidx_p(kk)
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
           coor_p(:,bstartS:bstartS+nbB-1), nbB, shofs(kk-blkptr_p(ibS)), nlocS, XcS, vdS)
   enddo
   !$omp end do
   !$omp do schedule(static)
   do g = 1,nbB
      XdS(:,g) = 0.0d0
      do kk = blkptr_p(ibS)+1, blkptr_p(ibS+1)
         si = blkidx_p(kk)
         do i = 1,cosx_shell_dim(si)
            XdS(cosx_ao_offset(si)+i,g) = XcS(shofs(kk-blkptr_p(ibS))+i,g)*sqrt(abs(w_p(bstartS+g-1)))
         enddo
      enddo
   enddo
   !$omp end do
   if (fit) then
      r0 = (nbB*mytid)/nthr + 1
      r1 = (nbB*(mytid+1))/nthr
      if (r1 .ge. r0) then
         do ir = 1,nrunS
            call dgemm('N','N',n,r1-r0+1,runL(ir),1.0d0,QfullT_p(1,runA(ir)),n,XdS(runA(ir),r0),n, &
                       merge(1.0d0,0.0d0,ir.gt.1),QbS(1,r0),n)
         enddo
      endif
   endif
   !$omp barrier
   !$omp master
   tc(1) = tc(1) + (omp_get_wtime()-tc0) ; tc0 = omp_get_wtime()
   !$omp end master
   do g0 = 0, nbB-1, nbs
      ng = min(nbs, nbB-g0)
      !$omp single
      rb0S = (g0/nrp)*nrp
      rgS = g0 - rb0S
      nrgS = min(nrp, nbB-g0)
      !$omp end single
      if (g0 .eq. rb0S) then
         r0 = (nop*mytid)/nthr + 1
         r1 = (nop*(mytid+1))/nthr
         if (r1 .ge. r0) then
            do ir = 1,nrunS
               call dgemm('N','N',r1-r0+1,nrgS,runL(ir),1.0d0,CoT(r0,runA(ir)),nop,XdS(runA(ir),g0+1),n, &
                          merge(1.0d0,0.0d0,ir.gt.1),YbS(r0,1),nop)
            enddo
         endif
         r0 = (nop*nc*mytid)/nthr + 1
         r1 = (nop*nc*(mytid+1))/nthr
         if (r1 .ge. r0) then
            do ir = 1,nrunS
               call dgemm('N','N',r1-r0+1,nrgS,runL(ir),1.0d0,Lt(r0,runA(ir)),nop*nc,XdS(runA(ir),g0+1),n, &
                          merge(1.0d0,0.0d0,ir.gt.1),ZbS(r0,1),nop*nc)
            enddo
         endif
         !$omp barrier
         r0 = (nc*nrgS*mytid)/nthr + 1
         r1 = (nc*nrgS*(mytid+1))/nthr
         if (r1 .ge. r0) call dgemm('N','N',n,r1-r0+1,nop,1.0d0,Co(1,1),n,ZbS(1+nop*(r0-1),1),nop,0.0d0, &
                                    RRS(1+n*(r0-1),1,1),n)
         !$omp barrier
         r0 = (n*nc*mytid)/nthr + 1
         r1 = (n*nc*(mytid+1))/nthr
         if (r1 .ge. r0) call dgemm('N','N',r1-r0+1,nrgS,nop,1.0d0,LfP(r0,1,1),n*nc,YbS(1,1),nop,1.0d0, &
                                    RRS(r0,1,1),n*nc)
         !$omp barrier
         !$omp master
         tc(2) = tc(2) + (omp_get_wtime()-tc0) ; tc0 = omp_get_wtime()
         !$omp end master
      endif
      !$omp do schedule(static)
      do si = 0,nBases-1
         FgmS(si) = maxval(abs(RRS(cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si),:,rgS+1:rgS+ng)))
      enddo
      !$omp end do nowait
      !$omp single
      xbmS = maxval(abs(XdS(:,g0+1:g0+ng)))
      if (fit) xbmS = max(xbmS, maxval(abs(QbS(:,g0+1:g0+ng))))
      !$omp end single
      !$omp master
      tc(7) = tc(7) + (omp_get_wtime()-tc0) ; tc0 = omp_get_wtime()
      !$omp end master
      !$omp do schedule(static)
      do g = 1,ng
         AsubS(:,:,g) = 0.0d0
      enddo
      !$omp end do
      !$omp master
      tc(3) = tc(3) + (omp_get_wtime()-tc0) ; tc0 = omp_get_wtime()
      !$omp end master
      !$omp do schedule(dynamic)
      do si = 0,nBases-1
         di = cosx_shell_dim(si)
         do sj_idx = cosx_ext_ptr(si)+1, cosx_ext_ptr(si+1)
            sj = cosx_ext_idx(sj_idx)
            if (schwarz_bound(si,sj) .lt. SCHWARZ_CUTOFF) cycle
            seg_ab = cosx_shell_cen(:,sj) - cosx_shell_cen(:,si)
            seg_len2 = sum(seg_ab*seg_ab)
            if (seg_len2 .gt. 0.0d0) then
               seg_t = sum((blkcen_p(:,ibS) - cosx_shell_cen(:,si))*seg_ab)/seg_len2
               seg_t = max(0.0d0, min(1.0d0, seg_t))
               seg_cl = cosx_shell_cen(:,si) + seg_t*seg_ab
            else
               seg_cl = cosx_shell_cen(:,si)
            endif
            r_far = sqrt(sum((blkcen_p(:,ibS) - seg_cl)**2)) - blkrad_p(ibS)
            if (r_far .gt. 0.0d0) then
               esp_eff = min(cosx_esp_bound(si,sj), cosx_esp_mono(si,sj)/r_far)
            else
               esp_eff = cosx_esp_bound(si,sj)
            endif
            k_bound = xbmS*esp_eff*max(FgmS(si),FgmS(sj))
            if (k_bound .lt. kscr_loc) cycle
            dj = cosx_shell_dim(sj)
            shls1e(1) = si ; shls1e(2) = sj ; shls1e(3) = g0 ; shls1e(4) = g0+ng
            call grids1e_engine_cached(bufP, shls1e, atm, size(atm,2), bas, nBases, envGS, 0_8, icP)
            do j = 1,dj
               ao_j = cosx_ao_offset(sj) + j
               do i = 1,di
                  ao_i = cosx_ao_offset(si) + i
                  nij = NorVEC(ao_i)*NorVEC(ao_j)
                  if (nij .eq. 0.0d0) cycle
                  do g = 1,ng
                     sg = sign(1.0d0,w_p(bstartS+g0+g-1))*nij*bufP(g + ng*((i-1)+di*(j-1)))
                     AsubS(ao_i,ao_j,g) = AsubS(ao_i,ao_j,g) + sg
                     if (sj .ne. si) AsubS(ao_j,ao_i,g) = AsubS(ao_j,ao_i,g) + sg
                  enddo
               enddo
            enddo
         enddo
      enddo
      !$omp end do
      !$omp master
      tc(4) = tc(4) + (omp_get_wtime()-tc0) ; tc0 = omp_get_wtime()
      !$omp end master
      !$omp do schedule(static)
      do g = 1,ng
         call dgemm('N','N',n,nop,n,1.0d0,AsubS(1,1,g),n,Co(1,1),n,0.0d0,ACoP(1,1),n)
         call dgemm('T','N',nop,nc,n,1.0d0,ACoP(1,1),n,RRS(1,1,rgS+g),n,0.0d0,HcS(1,1,g0+g),nop)
      enddo
      !$omp end do
      !$omp master
      tc(5) = tc(5) + (omp_get_wtime()-tc0) ; tc0 = omp_get_wtime()
      !$omp end master
   enddo
   !$omp do schedule(dynamic)
   do jr = 1,ntl
      call dgemm('N','T',nop*nc,tln(jr),nbB,1.0d0,HcS(1,1,1),nop*nc,XdS(tl0(jr),1),n, &
                 1.0d0,KT(1,tl0(jr)),nop*nc)
   enddo
   !$omp end do
   !$omp master
   tc(6) = tc(6) + (omp_get_wtime()-tc0) ; tc0 = omp_get_wtime()
   !$omp end master
   !$omp single
   deallocate(XcS, vdS)
   !$omp end single
enddo
deallocate(bufP, icP, ACoP)
!$omp end parallel
if (btime) print '(A,7F8.3)', " COSX COOP wall: AO+Q / R / zero A / ints / A.R / K acc / screen:", tc(1:7)
deallocate(XdS, RRS, AsubS, HcS, YbS, ZbS, FgmS, shofs, runA, runL, tl0, tln, envGS)
if (allocated(QbS)) deallocate(QbS)
else
!$omp parallel default(shared) reduction(+:tb) private(tb0,tbusy) &
!$omp&  private(ib,ibo,bstart,nb,g,si,sj,sj_idx,di,dj,i,j,ao_i,ao_j,m,kk,cofs,g0,ng,shls1e,off_grids,nloc) &
!$omp&  private(nij,esp_eff,r_far,seg_ab,seg_cl,seg_len2,seg_t,xbm,k_bound,sg) &
!$omp&  private(Xd,Qb,Xc,vdum,RR,OO,Asub,Kloc,KV,Fgm,buf,env_g,int_cache,run_a,run_len,nrun,ir,a0,ACo,Hc,Yb,Zb) &
!$omp&  private(Tk,rg0,cc0,ncc,cend,kg,mg0,mgn,rb0,nrg)
if (do_var) then
   allocate(KV(n*nc,n))
   KV = 0.0d0
endif
if (kshr .and. kshm .eq. 2) then
   allocate(Kloc(1,1), Tk(nrow*mgs,n))
else if (kshr) then
   allocate(Kloc(1,1), Tk(nrow*nc,KCOL))
else
   allocate(Kloc(nrow*nc,n))
endif
if (rsub) then
   allocate(RR(n,ncol*nc,nrp))
else
   allocate(RR(n,ncol*nc,BATCH))
endif
allocate(Asub(n,n,nbs), Fgm(0:nBases-1))
if (proj .and. rsub) then
   allocate(ACo(n,nop), Hc(nop,nc,BATCH), OO(1,1,1), Yb(nop,nrp), Zb(nop*nc,nrp))
else if (proj) then
   allocate(ACo(n,nop), Hc(nop,nc,BATCH), OO(1,1,1), Yb(nop,BATCH), Zb(nop*nc,BATCH))
else
   allocate(OO(n,ncol*nc,BATCH))
endif
allocate(buf(BATCH*maxd*maxd), int_cache(INT_CACHE_SIZE), env_g(size(env) + 3*BATCH))
allocate(run_a(nBases), run_len(nBases))
Kloc = 0.0d0
env_g(1:size(env)) = env(1:size(env))
env_g(9) = cosx_batch_omega
off_grids = size(env)
tbusy = omp_get_wtime()
!$omp do schedule(static,1)
do ibo = 1, nblk_p
   ib = border(ibo)
   bstart = blkstart_p(ib) + 1
   nb = blkstart_p(ib+1) - blkstart_p(ib)
   nloc = 0
   do kk = blkptr_p(ib)+1, blkptr_p(ib+1)
      nloc = nloc + cosx_shell_dim(blkidx_p(kk))
   enddo
   if (nloc .eq. 0) cycle
   nrun = 0
   do kk = blkptr_p(ib)+1, blkptr_p(ib+1)
      si = blkidx_p(kk)
      a0 = cosx_ao_offset(si) + 1
      if (nrun .gt. 0) then
         if (a0 .eq. run_a(nrun)+run_len(nrun)) then
            run_len(nrun) = run_len(nrun) + cosx_shell_dim(si)
            cycle
         endif
      endif
      nrun = nrun + 1
      run_a(nrun) = a0
      run_len(nrun) = cosx_shell_dim(si)
   enddo
   tb0 = omp_get_wtime()
   allocate(Xd(n,nb), Xc(nloc,nb), vdum(nloc,3,nb))
   Xc = 0.0d0
   cofs = 0
   do kk = blkptr_p(ib)+1, blkptr_p(ib+1)
      si = blkidx_p(kk)
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
           coor_p(:,bstart:bstart+nb-1), nb, cofs, nloc, Xc, vdum)
      cofs = cofs + cosx_shell_dim(si)
   enddo
   Xd = 0.0d0
   cofs = 0
   do kk = blkptr_p(ib)+1, blkptr_p(ib+1)
      si = blkidx_p(kk)
      do i = 1, cosx_shell_dim(si)
         Xd(cosx_ao_offset(si)+i,:) = Xc(cofs+i,:)*sqrt(abs(w_p(bstart:bstart+nb-1)))
      enddo
      cofs = cofs + cosx_shell_dim(si)
   enddo
   deallocate(Xc, vdum)
   if (fit) then
      allocate(Qb(n,nb))
      do ir = 1,nrun
         call dgemm('N','N',n,nb,run_len(ir),1.0d0,QfullT_p(1,run_a(ir)),n,Xd(run_a(ir),1),n, &
                    merge(1.0d0,0.0d0,ir.gt.1),Qb,n)
      enddo
   endif
   tb(6) = tb(6) + (omp_get_wtime()-tb0)
   do g = 1,nb
      env_g(off_grids+3*(g-1)+1:off_grids+3*g) = coor_p(:,bstart+g-1)
   enddo
   env_g(12) = dble(nb)
   env_g(13) = dble(off_grids)

   tb0 = omp_get_wtime()
   if (rsub) then
      continue
   else if (proj) then
      do ir = 1,nrun
         call dgemm('N','N',nop,nb,run_len(ir),1.0d0,CoT(1,run_a(ir)),nop,Xd(run_a(ir),1),n, &
                    merge(1.0d0,0.0d0,ir.gt.1),Yb,nop)
         call dgemm('N','N',nop*nc,nb,run_len(ir),1.0d0,Lt(1,run_a(ir)),nop*nc,Xd(run_a(ir),1),n, &
                    merge(1.0d0,0.0d0,ir.gt.1),Zb,nop*nc)
      enddo
      do m = 1,nc
         call dgemm('N','N',n,nb,nop,1.0d0,Lf(1,1,c0+m),n,Yb,nop,0.0d0,RR(1,m,1),n*ncol*nc)
         call dgemm('N','N',n,nb,nop,1.0d0,Co(1,1),n,Zb(nop*(m-1)+1,1),nop*nc,1.0d0,RR(1,m,1),n*ncol*nc)
      enddo
   else
   do ir = 1,nrun
      call dgemm('N','N',n*nc,nb,run_len(ir),1.0d0,Dp(1,1,run_a(ir)),n*nc,Xd(run_a(ir),1),n, &
                 merge(1.0d0,0.0d0,ir.gt.1),RR(1,1,1),n*ncol*nc)
   enddo
   endif
   if (do_var) call dgemm('N','N',n*nc,nb,n,1.0d0,Dp(1,1,1),n*nc,Qb,n,0.0d0,RR(1,nc+1,1),n*ncol*nc)
   tb(1) = tb(1) + (omp_get_wtime()-tb0)
   if (btime .and. .not. rsub) then
      do g = 1,nb
         do m = 1,ncol*nc
            sg = maxval(abs(RR(:,m,g)))
            tb(7) = tb(7) + 1.0d0
            do kk = 1,4
               if (sg .ge. 10.0d0**(-5-kk)) tb(7+kk) = tb(7+kk) + 1.0d0
            enddo
         enddo
      enddo
   endif
   do g0 = 0, nb-1, nbs
      ng = min(nbs, nb-g0)
      rg0 = g0
      if (rsub) then
         rb0 = (g0/nrp)*nrp
         rg0 = g0 - rb0
      endif
      if (rsub .and. g0 .eq. rb0) then
         tb0 = omp_get_wtime()
         nrg = min(nrp, nb-g0)
         do ir = 1,nrun
            call dgemm('N','N',nop,nrg,run_len(ir),1.0d0,CoT(1,run_a(ir)),nop,Xd(run_a(ir),g0+1),n, &
                       merge(1.0d0,0.0d0,ir.gt.1),Yb,nop)
            call dgemm('N','N',nop*nc,nrg,run_len(ir),1.0d0,Lt(1,run_a(ir)),nop*nc,Xd(run_a(ir),g0+1),n, &
                       merge(1.0d0,0.0d0,ir.gt.1),Zb,nop*nc)
         enddo
         call dgemm('N','N',n,nc*nrg,nop,1.0d0,Co(1,1),n,Zb(1,1),nop,0.0d0,RR(1,1,1),n)
         call dgemm('N','N',n*nc,nrg,nop,1.0d0,LfP(1,1,1),n*nc,Yb(1,1),nop,1.0d0,RR(1,1,1),n*nc)
         tb(1) = tb(1) + (omp_get_wtime()-tb0)
      endif
      tb0 = omp_get_wtime()
      xbm = maxval(abs(Xd(:,g0+1:g0+ng)))
      if (fit) xbm = max(xbm, maxval(abs(Qb(:,g0+1:g0+ng))))
      do si = 0,nBases-1
         Fgm(si) = maxval(abs(RR(cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si),:,rg0+1:rg0+ng)))
      enddo
      Asub(:,:,1:ng) = 0.0d0
      tb(2) = tb(2) + (omp_get_wtime()-tb0) ; tb0 = omp_get_wtime()
      do si = 0,nBases-1
         di = cosx_shell_dim(si)
         do sj_idx = cosx_ext_ptr(si)+1, cosx_ext_ptr(si+1)
            sj = cosx_ext_idx(sj_idx)
            if (schwarz_bound(si,sj) .lt. SCHWARZ_CUTOFF) cycle
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
            k_bound = xbm*esp_eff*max(Fgm(si),Fgm(sj))
            if (k_bound .lt. kscr_loc) cycle
            dj = cosx_shell_dim(sj)
            shls1e(1) = si ; shls1e(2) = sj ; shls1e(3) = g0 ; shls1e(4) = g0+ng
            call grids1e_engine_cached(buf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            do j = 1,dj
               ao_j = cosx_ao_offset(sj) + j
               do i = 1,di
                  ao_i = cosx_ao_offset(si) + i
                  nij = NorVEC(ao_i)*NorVEC(ao_j)
                  if (nij .eq. 0.0d0) cycle
                  do g = 1,ng
                     sg = sign(1.0d0,w_p(bstart+g0+g-1))*nij*buf(g + ng*((i-1)+di*(j-1)))
                     Asub(ao_i,ao_j,g) = Asub(ao_i,ao_j,g) + sg
                     if (sj .ne. si) Asub(ao_j,ao_i,g) = Asub(ao_j,ao_i,g) + sg
                  enddo
               enddo
            enddo
         enddo
      enddo
      tb(3) = tb(3) + (omp_get_wtime()-tb0) ; tb0 = omp_get_wtime()
      if (proj) then
         do g = 1,ng
            call dgemm('N','N',n,nop,n,1.0d0,Asub(1,1,g),n,Co(1,1),n,0.0d0,ACo(1,1),n)
            call dgemm('T','N',nop,nc,n,1.0d0,ACo(1,1),n,RR(1,1,rg0+g),n,0.0d0,Hc(1,1,g0+g),nop)
         enddo
      else
         do g = 1,ng
            call dgemm('N','N',n,ncol*nc,n,1.0d0,Asub(1,1,g),n,RR(1,1,g0+g),n,0.0d0,OO(1,1,g0+g),n)
         enddo
      endif
      tb(4) = tb(4) + (omp_get_wtime()-tb0)
   enddo
   tb0 = omp_get_wtime()
   if (proj .and. kshr .and. kshm .eq. 2) then
      do kg = 0, (nc+mgs-1)/mgs-1
         mg0 = mgs*mod(kg + omp_get_thread_num(), (nc+mgs-1)/mgs)
         mgn = min(mgs, nc-mg0)
         do ir = 1,nrun
            call dgemm('N','T',nop*mgn,run_len(ir),nb,1.0d0,Hc(1,mg0+1,1),nop*nc,Xd(run_a(ir),1),n, &
                       0.0d0,Tk(1,run_a(ir)),nop*mgs)
         enddo
         call omp_set_lock(klk(mg0/mgs+1))
         do ir = 1,nrun
            do kk = run_a(ir), run_a(ir)+run_len(ir)-1
               KTq(nop*mg0+1:nop*(mg0+mgn),kk) = KTq(nop*mg0+1:nop*(mg0+mgn),kk) + to_q(Tk(1:nop*mgn,kk))
            enddo
         enddo
         call omp_unset_lock(klk(mg0/mgs+1))
      enddo
   endif
   do ir = 1,nrun
      if (proj .and. kshr .and. kshm .eq. 2) then
         exit
      else if (proj .and. kshr) then
         cc0 = run_a(ir)
         cend = run_a(ir) + run_len(ir) - 1
         do while (cc0 .le. cend)
            ncc = min(cend, KCOL*((cc0-1)/KCOL+1)) - cc0 + 1
            call dgemm('N','T',nop*nc,ncc,nb,1.0d0,Hc(1,1,1),nop*nc,Xd(cc0,1),n, &
                       0.0d0,Tk,nop*nc)
            call omp_set_lock(klk((cc0-1)/KCOL+1))
            KTq(:,cc0:cc0+ncc-1) = KTq(:,cc0:cc0+ncc-1) + to_q(Tk(:,1:ncc))
            call omp_unset_lock(klk((cc0-1)/KCOL+1))
            cc0 = cc0 + ncc
         enddo
         cycle
      else if (proj) then
         call dgemm('N','T',nop*nc,run_len(ir),nb,1.0d0,Hc(1,1,1),nop*nc,Xd(run_a(ir),1),n, &
                    1.0d0,Kloc(1,run_a(ir)),nop*nc)
         cycle
      endif
      call dgemm('N','T',n*nc,run_len(ir),nb,1.0d0,OO(1,1,1),n*ncol*nc,Xd(run_a(ir),1),n, &
                 1.0d0,Kloc(1,run_a(ir)),n*nc)
      if (do_var) call dgemm('N','T',n*nc,run_len(ir),nb,1.0d0,OO(1,nc+1,1),n*ncol*nc,Xd(run_a(ir),1),n, &
                 1.0d0,KV(1,run_a(ir)),n*nc)
   enddo
   tb(5) = tb(5) + (omp_get_wtime()-tb0)
   deallocate(Xd)
   if (allocated(Qb)) deallocate(Qb)
enddo
!$omp end do nowait
tbusy = omp_get_wtime() - tbusy
!$omp critical
bmax = max(bmax, tbusy) ; bmin = min(bmin, tbusy)
!$omp end critical
!$omp barrier
tbusy = omp_get_wtime()
block
   integer :: tid
   do tid = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() .eq. tid) then
         if (.not. kshr) KT = KT + Kloc
         if (do_var) KVs = KVs + KV
      endif
   enddo
end block
!$omp master
tmg = omp_get_wtime() - tbusy
!$omp end master
deallocate(Kloc, RR, OO, Asub, Fgm, buf, int_cache, env_g, run_a, run_len)
if (allocated(ACo)) deallocate(ACo, Hc, Yb, Zb)
if (allocated(Tk)) deallocate(Tk)
if (do_var) deallocate(KV)
!$omp end parallel
endif
if (allocated(KTq)) then
   !$omp parallel do private(kk)
   do kk = 1,n
      KT(:,kk) = KT(:,kk) + from_q(KTq(:,kk))
   enddo
   !$omp end parallel do
   deallocate(KTq)
endif
tw(2) = tw(2) + (omp_get_wtime()-tw0) ; tw0 = omp_get_wtime()
deallocate(Dp)
if (allocated(Lt)) deallocate(Lt, CoT)
if (allocated(LfP)) deallocate(LfP)
if (proj) then
   !$omp parallel do private(m,kk)
   do m = 1,nc
      if (fit) then
         call dgemm('N','T',n,nop,n,1.0d0,QfullT_p(1,1),n,KT(nop*(m-1)+1,1),nop*nc,0.0d0,KC(1,1,c0+m),n)
      else
         do kk = 1,n
            KC(kk,:,c0+m) = KT(nop*(m-1)+1:nop*m,kk)
         enddo
      endif
   enddo
   !$omp end parallel do
else
!$omp parallel do private(m,kk)
do m = 1,nc
   if (fit) then
      call dgemm('N','T',n,n,n,1.0d0,KT(n*(m-1)+1,1),n*nc,QfullT_p(1,1),n,0.0d0,K(1,1,c0+m),n)
   else
      do kk = 1,n
         K(:,kk,c0+m) = KT(n*(m-1)+1:n*m,kk)
      enddo
   endif
   if (do_var) then
      do kk = 1,n
         K(:,kk,c0+m) = K(:,kk,c0+m) + KVs(n*(m-1)+1:n*m,kk)
      enddo
   endif
enddo
!$omp end parallel do
endif
deallocate(KT)
if (do_var) deallocate(KVs)
tw(3) = tw(3) + (omp_get_wtime()-tw0) ; tw0 = omp_get_wtime()
enddo
if (btime) print '(A,5F9.2)', " COSX BATCH (thread-s) D.x / screen+zero / ints+scatter / A.R / K acc:", tb(1:5)
if (btime) print '(A,4F7.3)', " COSX BATCH (g,m) columns max|D_m x_g| >= 1e-6/-7/-8/-9:", tb(8:11)/max(tb(7),1.0d0)

call get_environment_variable("ENGINE_CPSCF_COSX_UNSYM", envs)
if (proj) then
   continue
else if (trim(envs) .eq. "1" .and. .not. do_var) then
   !$omp parallel do private(m)
   do m = 1,nd
      K(:,:,m) = transpose(K(:,:,m))
   enddo
   !$omp end parallel do
else
!$omp parallel do private(m)
do m = 1,nd
   if (do_var) K(:,:,m) = 0.5d0*K(:,:,m)
   K(:,:,m) = 0.5d0*(K(:,:,m) + transpose(K(:,:,m)))
enddo
!$omp end parallel do
endif
deallocate(border, bcost)
if (kshr) then
   do kk = 1,size(klk)
      call omp_destroy_lock(klk(kk))
   enddo
   deallocate(klk)
endif
call cosx_malloc_trim()
tw(4) = tw(4) + (omp_get_wtime()-tw0)
if (btime) print '(A,4F9.3)', " COSX BATCH wall: prep / parallel pass / K finalize / sym+trim:", tw
if (btime) print '(A,F8.3,A,2F8.3,A,F8.3)', " COSX BATCH AO+Qb thread-s", tb(6), "  thread busy max/min", bmax, bmin, &
   "  ordered merge", tmg
end subroutine cosx_batch_core
end submodule exchange_cosx_resp_impl
