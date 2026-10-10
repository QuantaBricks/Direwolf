! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! RI-J (density-fitted Coulomb) direct ("Pulay") second-derivative Hessian

submodule (mod_integrals) hessian_df_j_impl
implicit none
contains

module subroutine integrals_hessian_2e_df_j(nConts, natoms, Ptot, d2E)
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads
implicit none
real(8) :: tq0, tq(5), tb0
integer,allocatable :: own(:,:), uord(:), ui(:), uj(:)
real(8),allocatable :: ucost(:), tload(:)
integer :: nunits, u, tmin, kk0
real(8) :: ck
real(8),allocatable :: tbusy(:)
logical :: tq_on
character(len=4) :: tq_env
real(8), external :: omp_get_wtime
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: d2E(3,natoms,3,natoms)

integer :: ish,jsh,ksh,di,dj,dk,p,q,r,a,b,tid,nth
integer :: atomI,atomJ,atomK,offi,offj,offP
integer :: nshb, shls(4), shlsw(4)
integer,allocatable :: ao_offset0(:), aux_offset0(:)
real(8),allocatable :: bvec(:), cvec(:)
real(8),allocatable :: bx(:,:,:), bc(:,:,:), wmat(:,:,:), uvec(:,:,:)
real(8),allocatable :: bvec_t(:,:), bx_t(:,:,:,:), bc_t(:,:,:,:), d2E_t(:,:,:,:,:)
real(8),allocatable :: buf(:), coef(:,:,:), dmax_shell(:,:)
real(8) :: Ppq, ns, s3(3), s9(9), crow, fac, v3
real(8) :: M3(3,3), M3t(3,3), MII(3,3), MKK(3,3), MIK(3,3)
logical :: do_direct, do_cross, do_metric, do_3c
character(32) :: envx
integer :: stx
integer(8) :: opt_base, opt_ip1, opt_ip2, opt_ipip1, opt_ipip2, opt_ipvip1, opt_ip1ip2

do_direct = .true.; do_cross = .true.
call get_environment_variable("ENGINE_DFJ_NODIRECT", envx, status=stx); if (stx == 0) do_direct = .false.
call get_environment_variable("ENGINE_DFJ_NOCROSS", envx, status=stx); if (stx == 0) do_cross = .false.
do_metric = do_direct; do_3c = do_direct
call get_environment_variable("ENGINE_DFJ_NOMETRIC", envx, status=stx); if (stx == 0) do_metric = .false.
call get_environment_variable("ENGINE_DFJ_NO3C", envx, status=stx); if (stx == 0) do_3c = .false.

call prof_start("hess_2e_df_j")
call ensure_df_built()
if (.not. allocated(df_evec) .or. .not. allocated(df_evalinv)) then
   print *, 'ERROR: RI-J Hessian needs the df_evec/df_evalinv metric factorization.'
   stop 1
endif

nshb = nBases + nBasesAux
allocate(ao_offset0(0:nBases-1))
offi = 0
do ish = 0,nBases-1
   ao_offset0(ish) = offi
   offi = offi + cgto_engine(ish, basDF)
enddo
allocate(aux_offset0(nBases:nshb-1))
offP = 0
do ksh = nBases,nshb-1
   aux_offset0(ksh) = offP
   offP = offP + cgto_engine(ksh, basDF)
enddo

allocate(bvec(nContsAux), cvec(nContsAux))
allocate(bx(natoms,3,nContsAux), bc(natoms,3,nContsAux))
allocate(wmat(natoms,3,nContsAux), uvec(natoms,3,nContsAux))
allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts, Ptot, dmax_shell)

nth = omp_get_max_threads()
allocate(bvec_t(nContsAux,nth), bx_t(natoms,3,nContsAux,nth))
allocate(bc_t(natoms,3,nContsAux,nth), d2E_t(3,natoms,3,natoms,nth))

bvec = 0.0d0; bx = 0.0d0; bc = 0.0d0; d2E = 0.0d0

call get_environment_variable("ENGINE_HESS2E_TIMING", tq_env)
tq_on = (trim(tq_env) .eq. "1")
tq = 0.0d0
tq0 = omp_get_wtime()
call threec2e_optimizer_engine(opt_base, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ip1_optimizer_engine(opt_ip1, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ip2_optimizer_engine(opt_ip2, atm, size(atm,2), basDF, nshb, envDF)
bvec_t = 0.0d0; bx_t = 0.0d0
allocate(own(0:nBases-1,0:nBases-1)) ; own = 0
nunits = 0
do ish = 0,nBases-1
   do jsh = 0,ish
      if (schwarz_bound(ish,jsh) .ge. SCHWARZ_CUTOFF) nunits = nunits + 1
   enddo
enddo
allocate(ui(nunits), uj(nunits), ucost(nunits), uord(nunits), tload(nth))
u = 0
do ish = 0,nBases-1
   do jsh = 0,ish
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      u = u + 1
      ui(u) = ish ; uj(u) = jsh
      ck = 0.0d0
      do ksh = nBases,nshb-1
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) &
             *dmax_shell(ish,jsh) .lt. density_screen_cutoff) cycle
         ck = ck + real(cgto_engine(ksh, basDF),8)
      enddo
      ucost(u) = real(cgto_engine(ish, basDF),8)*real(cgto_engine(jsh, basDF),8)*ck + 1.0d-6
      uord(u) = u
   enddo
enddo
call sort_desc_idx(nunits, ucost, uord)
tload = 0.0d0
do kk0 = 1,nunits
   u = uord(kk0)
   tmin = minloc(tload, dim=1)
   own(ui(u),uj(u)) = tmin
   tload(tmin) = tload(tmin) + ucost(u)
enddo
deallocate(ui, uj, ucost, uord)
!$omp parallel default(shared) num_threads(nth) &
!$omp&   private(ish,jsh,ksh,di,dj,dk,p,q,r,a,atomI,atomJ,atomK,offi,offj,offP, &
!$omp&           shls,shlsw,buf,Ppq,ns,tid,fac,v3)
tid = omp_get_thread_num() + 1
if (omp_get_num_threads() .ne. nth) error stop 'hessian_2e_df_j: thread team smaller than planned'
do ish = 0,nBases-1
   di = cgto_engine(ish, basDF); atomI = basDF(1,ish+1) + 1; offi = ao_offset0(ish)
   do jsh = 0,ish
      fac = merge(2.0d0, 1.0d0, jsh .ne. ish)
      dj = cgto_engine(jsh, basDF); atomJ = basDF(1,jsh+1) + 1; offj = ao_offset0(jsh)
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      if (own(ish,jsh) .ne. tid) cycle
      do ksh = nBases,nshb-1
         dk = cgto_engine(ksh, basDF); atomK = basDF(1,ksh+1) + 1; offP = aux_offset0(ksh)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) &
             *dmax_shell(ish,jsh) .lt. density_screen_cutoff) cycle
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         shlsw(1) = jsh; shlsw(2) = ish; shlsw(3) = ksh

         allocate(buf(di*dj*dk))
         call threec2e_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_base)
         do r = 1,dk; do q = 1,dj; do p = 1,di
            Ppq = Ptot(offi+p,offj+q)
            if (Ppq .eq. 0.0d0) cycle
            ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            bvec_t(offP+r,tid) = bvec_t(offP+r,tid) + buf(((r-1)*dj+(q-1))*di+p)*Ppq*ns*fac
         enddo; enddo; enddo
         deallocate(buf)

         allocate(buf(di*dj*dk*3))
         call threec2e_ip1_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip1)
         do r = 1,dk; do q = 1,dj; do p = 1,di
            Ppq = Ptot(offi+p,offj+q)
            if (Ppq .eq. 0.0d0) cycle
            ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            do a = 1,3
               v3 = buf(((r-1)*dj+(q-1))*di+p+(a-1)*di*dj*dk)*Ppq*ns*fac
               bx_t(atomI,a,offP+r,tid) = bx_t(atomI,a,offP+r,tid) + v3
               bx_t(atomJ,a,offP+r,tid) = bx_t(atomJ,a,offP+r,tid) - v3
            enddo
         enddo; enddo; enddo
         deallocate(buf)

         allocate(buf(di*dj*dk*3))
         call threec2e_ip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip2)
         do r = 1,dk; do q = 1,dj; do p = 1,di
            Ppq = Ptot(offi+p,offj+q)
            if (Ppq .eq. 0.0d0) cycle
            ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            do a = 1,3
               v3 = buf(((r-1)*dj+(q-1))*di+p+(a-1)*di*dj*dk)*Ppq*ns*fac
               bx_t(atomK,a,offP+r,tid) = bx_t(atomK,a,offP+r,tid) + v3
               bx_t(atomJ,a,offP+r,tid) = bx_t(atomJ,a,offP+r,tid) - v3
            enddo
         enddo; enddo; enddo
         deallocate(buf)
      enddo
   enddo
enddo
!$omp end parallel
do tid = 1,nth
   bvec = bvec + bvec_t(:,tid)
   bx = bx + bx_t(:,:,:,tid)
enddo
call cintdel_optimizer(opt_base)
call cintdel_optimizer(opt_ip1)
call cintdel_optimizer(opt_ip2)

tq(1) = omp_get_wtime() - tq0 ; tq0 = omp_get_wtime()
call df_metric_inverse(df_evec, df_evalinv, nContsAux, bvec, cvec)

tq(2) = omp_get_wtime() - tq0 ; tq0 = omp_get_wtime()
call twoc2e_ip1_optimizer_engine(opt_ip1, atm, size(atm,2), basDF, nshb, envDF)
call twoc2e_ip2_optimizer_engine(opt_ip2, atm, size(atm,2), basDF, nshb, envDF)
call twoc2e_ipip1_optimizer_engine(opt_ipip1, atm, size(atm,2), basDF, nshb, envDF)
call twoc2e_ip1ip2_optimizer_engine(opt_ip1ip2, atm, size(atm,2), basDF, nshb, envDF)
bc_t = 0.0d0; d2E_t = 0.0d0
!$omp parallel default(shared) &
!$omp&   private(ish,jsh,di,dj,p,q,a,atomI,atomJ,offP,offj,shls,shlsw,buf,coef, &
!$omp&           ns,crow,s9,M3,M3t,tid)
tid = omp_get_thread_num() + 1
!$omp do schedule(static,1)
do ish = nBases,nshb-1
   di = cgto_engine(ish, basDF); atomI = basDF(1,ish+1) + 1; offP = aux_offset0(ish)
   do jsh = nBases,nshb-1
      dj = cgto_engine(jsh, basDF); atomJ = basDF(1,jsh+1) + 1; offj = aux_offset0(jsh)
      if (aux_shell_bound(ish-nBases+1)*aux_shell_bound(jsh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
      shls(1) = ish; shls(2) = jsh
      shlsw(1) = jsh; shlsw(2) = ish

      allocate(buf(di*dj*3))
      call twoc2e_ip1_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip1)
      do q = 1,dj; do p = 1,di
         ns = NorVECAux(offP+p)*NorVECAux(offj+q)
         crow = cvec(offj+q)*ns
         do a = 1,3
            bc_t(atomI,a,offP+p,tid) = bc_t(atomI,a,offP+p,tid) + buf(((q-1)*di+p)+(a-1)*di*dj)*crow
         enddo
      enddo; enddo
      deallocate(buf)
      allocate(buf(di*dj*3))
      call twoc2e_ip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip2)
      do q = 1,dj; do p = 1,di
         ns = NorVECAux(offP+p)*NorVECAux(offj+q)
         crow = cvec(offj+q)*ns
         do a = 1,3
            bc_t(atomJ,a,offP+p,tid) = bc_t(atomJ,a,offP+p,tid) + buf(((q-1)*di+p)+(a-1)*di*dj)*crow
         enddo
      enddo; enddo
      deallocate(buf)

      allocate(coef(di,dj,1))
      do q = 1,dj; do p = 1,di
         coef(p,q,1) = cvec(offP+p)*cvec(offj+q)*NorVECAux(offP+p)*NorVECAux(offj+q)
      enddo; enddo
      allocate(buf(di*dj*9))
      call twoc2e_ipip1_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ipip1)
      call contract9_2c(buf, di, dj, coef, s9)
      if (do_metric) then
         M3 = -0.5d0*transpose(reshape(s9,(/3,3/)))
         d2E_t(:,atomI,:,atomI,tid) = d2E_t(:,atomI,:,atomI,tid) + M3
      endif
      deallocate(buf)
      allocate(buf(di*dj*9))
      call twoc2e_ipip1_engine(buf, shlsw, atm, size(atm,2), basDF, nshb, envDF, opt_ipip1)
      call contract9_2c_sw(buf, di, dj, coef, s9)
      if (do_metric) then
         M3 = -0.5d0*transpose(reshape(s9,(/3,3/)))
         d2E_t(:,atomJ,:,atomJ,tid) = d2E_t(:,atomJ,:,atomJ,tid) + M3
      endif
      deallocate(buf)
      allocate(buf(di*dj*9))
      call twoc2e_ip1ip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip1ip2)
      call contract9_2c(buf, di, dj, coef, s9)
      if (do_metric) then
         M3 = -0.5d0*transpose(reshape(s9,(/3,3/)))
         M3t = transpose(M3)
         d2E_t(:,atomI,:,atomJ,tid) = d2E_t(:,atomI,:,atomJ,tid) + M3
         d2E_t(:,atomJ,:,atomI,tid) = d2E_t(:,atomJ,:,atomI,tid) + M3t
      endif
      deallocate(buf)
      deallocate(coef)
   enddo
enddo
!$omp end do
!$omp end parallel
do tid = 1,nth
   bc = bc + bc_t(:,:,:,tid)
   d2E = d2E + d2E_t(:,:,:,:,tid)
enddo
call cintdel_optimizer(opt_ip1)
call cintdel_optimizer(opt_ip2)
call cintdel_optimizer(opt_ipip1)
call cintdel_optimizer(opt_ip1ip2)

do a = 1,3
   do ish = 1,natoms
      do r = 1,nContsAux
         wmat(ish,a,r) = bx(ish,a,r) - bc(ish,a,r)
      enddo
   enddo
enddo
do a = 1,3
   do ish = 1,natoms
      call df_metric_inverse(df_evec, df_evalinv, nContsAux, wmat(ish,a,:), uvec(ish,a,:))
   enddo
enddo

tq(3) = omp_get_wtime() - tq0 ; tq0 = omp_get_wtime()
call threec2e_ipip1_optimizer_engine(opt_ipip1, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ipvip1_optimizer_engine(opt_ipvip1, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ip1ip2_optimizer_engine(opt_ip1ip2, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ipip2_optimizer_engine(opt_ipip2, atm, size(atm,2), basDF, nshb, envDF)
d2E_t = 0.0d0
allocate(tbusy(size(d2E_t,5))) ; tbusy = 0.0d0
!$omp parallel default(shared) num_threads(size(d2E_t,5)) &
!$omp&   private(ish,jsh,ksh,di,dj,dk,p,q,r,atomI,atomJ,atomK,offi,offj,offP, &
!$omp&           shls,shlsw,buf,coef,s9,M3,M3t,MII,MKK,MIK,tid,fac,tb0)
tid = omp_get_thread_num() + 1
if (omp_get_num_threads() .ne. size(d2E_t,5)) error stop 'hessian_2e_df_j: thread team smaller than planned'
tb0 = omp_get_wtime()
do ish = 0,nBases-1
   di = cgto_engine(ish, basDF); atomI = basDF(1,ish+1) + 1; offi = ao_offset0(ish)
   do jsh = 0,ish
      fac = merge(2.0d0, 1.0d0, jsh .ne. ish)
      dj = cgto_engine(jsh, basDF); atomJ = basDF(1,jsh+1) + 1; offj = ao_offset0(jsh)
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      if (own(ish,jsh) .ne. tid) cycle
      do ksh = nBases,nshb-1
         dk = cgto_engine(ksh, basDF); atomK = basDF(1,ksh+1) + 1; offP = aux_offset0(ksh)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) &
             *dmax_shell(ish,jsh) .lt. density_screen_cutoff) cycle
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         shlsw(1) = jsh; shlsw(2) = ish; shlsw(3) = ksh

         allocate(coef(di,dj,dk))
         do r = 1,dk; do q = 1,dj; do p = 1,di
            coef(p,q,r) = fac*cvec(offP+r)*Ptot(offi+p,offj+q) &
                          *NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
         enddo; enddo; enddo

         if (do_3c) then
            allocate(buf(di*dj*dk*9))
            call threec2e_ipip1_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ipip1)
            call contract9_3c(buf, di, dj, dk, coef, s9)
            MII = transpose(reshape(s9,(/3,3/)))
            call threec2e_ipip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ipip2)
            call contract9_3c(buf, di, dj, dk, coef, s9)
            MKK = transpose(reshape(s9,(/3,3/)))
            call threec2e_ip1ip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip1ip2)
            call contract9_3c(buf, di, dj, dk, coef, s9)
            MIK = transpose(reshape(s9,(/3,3/)))
            deallocate(buf)
            M3 = -(MII + MIK)
            M3t = -(MIK + MKK)
            d2E_t(:,atomI,:,atomI,tid) = d2E_t(:,atomI,:,atomI,tid) + MII
            d2E_t(:,atomK,:,atomK,tid) = d2E_t(:,atomK,:,atomK,tid) + MKK
            d2E_t(:,atomJ,:,atomJ,tid) = d2E_t(:,atomJ,:,atomJ,tid) + MII + MIK + transpose(MIK) + MKK
            d2E_t(:,atomI,:,atomK,tid) = d2E_t(:,atomI,:,atomK,tid) + MIK
            d2E_t(:,atomK,:,atomI,tid) = d2E_t(:,atomK,:,atomI,tid) + transpose(MIK)
            d2E_t(:,atomI,:,atomJ,tid) = d2E_t(:,atomI,:,atomJ,tid) + M3
            d2E_t(:,atomJ,:,atomI,tid) = d2E_t(:,atomJ,:,atomI,tid) + transpose(M3)
            d2E_t(:,atomJ,:,atomK,tid) = d2E_t(:,atomJ,:,atomK,tid) + M3t
            d2E_t(:,atomK,:,atomJ,tid) = d2E_t(:,atomK,:,atomJ,tid) + transpose(M3t)
         endif
         deallocate(coef)
      enddo
   enddo
enddo
tbusy(tid) = omp_get_wtime() - tb0
!$omp end parallel
deallocate(own, tload)
if (tq_on) print '(A,3F8.3)', ' HESS2EDFJ pass3 thread busy min,mean,max:', minval(tbusy), sum(tbusy)/size(tbusy), maxval(tbusy)
deallocate(tbusy)
do tid = 1,nth
   d2E = d2E + d2E_t(:,:,:,:,tid)
enddo
call cintdel_optimizer(opt_ipip1)
call cintdel_optimizer(opt_ipvip1)
call cintdel_optimizer(opt_ip1ip2)
call cintdel_optimizer(opt_ipip2)

if (do_cross) then
   do ish = 1,natoms
      do a = 1,3
         do jsh = 1,natoms
            do b = 1,3
               crow = 0.0d0
               do r = 1,nContsAux
                  crow = crow + wmat(jsh,b,r)*uvec(ish,a,r)
               enddo
               d2E(a,ish,b,jsh) = d2E(a,ish,b,jsh) + crow
            enddo
         enddo
      enddo
   enddo
endif

block
   character(32) :: dbgenv
   integer :: dbgst
   call get_environment_variable("ENGINE_DFJ_DEBUG", dbgenv, status=dbgst)
   if (dbgst == 0) then
      print *, 'DFJ DEBUG max|b|=',maxval(abs(bvec)),' max|c|=',maxval(abs(cvec))
      print *, 'DFJ DEBUG max|bx|=',maxval(abs(bx)),' max|bc|=',maxval(abs(bc))
      print *, 'DFJ DEBUG max|w|=',maxval(abs(wmat)),' max|u|=',maxval(abs(uvec))
      print *, 'DFJ DEBUG metric cond(M) =',maxval(df_evalinv)/max(minval(df_evalinv),1d-300)
      block
         real(8) :: sbx, sbc
         integer :: aa, bb, rr
         sbx = 0.0d0; sbc = 0.0d0
         do aa = 1,3
            do rr = 1,nContsAux
               crow = 0.0d0
               do bb = 1,natoms
                  crow = crow + bx(bb,aa,rr)
               enddo
               sbx = max(sbx, abs(crow))
               crow = 0.0d0
               do bb = 1,natoms
                  crow = crow + bc(bb,aa,rr)
               enddo
               sbc = max(sbc, abs(crow))
            enddo
         enddo
         print *, 'DFJ DEBUG max|sum_B bx| =',sbx,' max|sum_B bc| =',sbc
      end block
   endif
end block

deallocate(bvec,cvec,bx,bc,wmat,uvec,bvec_t,bx_t,bc_t,d2E_t)
deallocate(ao_offset0,aux_offset0,dmax_shell)
call prof_stop("hess_2e_df_j")
tq(4) = omp_get_wtime() - tq0
if (tq_on) print '(A,4F8.3)', ' HESS2EDFJ pass1,Minv,pass2,pass3:', tq(1), tq(2), tq(3), tq(4)
end subroutine integrals_hessian_2e_df_j

module subroutine integrals_df_j_fock_deriv(nConts, natoms, Ptot, Jx)
use mod_profile, only: prof_start, prof_stop
use omp_lib, only: omp_get_max_threads, omp_get_thread_num, omp_get_wtime, omp_get_num_threads
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Jx(nConts,nConts,3,natoms)

integer :: ish,jsh,ksh,di,dj,dk,p,q,r,a,tid,nth,mu,nu,N3,maxdo,maxda,ncol
integer :: atomI,atomJ,atomK,offi,offj,offP
integer :: nshb, shls(4), shlsw(4)
integer,allocatable :: ao_offset0(:), aux_offset0(:), ao_atom(:), ao_shell0(:)
real(8),allocatable :: bvec(:), cvec(:), bx(:,:,:), bc(:,:,:), wtmp(:), uT(:,:)
real(8),allocatable :: bvec_t(:,:), bx_t(:,:,:,:), bc_t(:,:,:,:), buf(:), S12(:,:,:)
real(8),allocatable :: dmax_shell(:,:), T(:,:), bnP(:,:), uTg(:,:)
real(8) :: Ppq, ns, crow, v, fac
integer(8) :: opt_base, opt_ip1, opt_ip2
real(8) :: tj0, tj(5)
integer,allocatable :: own(:,:), uord(:), ui(:), uj(:)
real(8),allocatable :: ucost(:), tload(:)
integer :: nunits, u, tmin, kk0
real(8) :: ck

call prof_start("hess_df_j_fock")
tj0 = omp_get_wtime()
call ensure_df_built()
N3 = 3*natoms
nshb = nBases + nBasesAux
allocate(ao_offset0(0:nBases-1), ao_atom(nConts), ao_shell0(nConts))
offi = 0
do ish = 0,nBases-1
   ao_offset0(ish) = offi
   do p = 1,cgto_engine(ish, basDF)
      ao_atom(offi+p) = basDF(1,ish+1) + 1
      ao_shell0(offi+p) = ish
   enddo
   offi = offi + cgto_engine(ish, basDF)
enddo
allocate(aux_offset0(nBases:nshb-1))
offP = 0
do ksh = nBases,nshb-1
   aux_offset0(ksh) = offP
   offP = offP + cgto_engine(ksh, basDF)
enddo
allocate(bvec(nContsAux), cvec(nContsAux), bx(natoms,3,nContsAux), bc(natoms,3,nContsAux))
allocate(dmax_shell(0:nBases-1,0:nBases-1))
call compute_shell_density_bound(nConts, Ptot, dmax_shell)
nth = omp_get_max_threads()
allocate(bvec_t(nContsAux,nth), bx_t(natoms,3,nContsAux,nth), bc_t(natoms,3,nContsAux,nth))

call threec2e_optimizer_engine(opt_base, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ip1_optimizer_engine(opt_ip1, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ip2_optimizer_engine(opt_ip2, atm, size(atm,2), basDF, nshb, envDF)
bvec_t = 0.0d0; bx_t = 0.0d0
allocate(own(0:nBases-1,0:nBases-1)) ; own = 0
nunits = 0
do ish = 0,nBases-1
   do jsh = 0,ish
      if (schwarz_bound(ish,jsh) .ge. SCHWARZ_CUTOFF) nunits = nunits + 1
   enddo
enddo
allocate(ui(nunits), uj(nunits), ucost(nunits), uord(nunits), tload(nth))
u = 0
do ish = 0,nBases-1
   do jsh = 0,ish
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      u = u + 1
      ui(u) = ish ; uj(u) = jsh
      ck = 0.0d0
      do ksh = nBases,nshb-1
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) &
             *dmax_shell(ish,jsh) .lt. density_screen_cutoff) cycle
         ck = ck + real(cgto_engine(ksh, basDF),8)
      enddo
      ucost(u) = real(cgto_engine(ish, basDF),8)*real(cgto_engine(jsh, basDF),8)*ck + 1.0d-6
      uord(u) = u
   enddo
enddo
call sort_desc_idx(nunits, ucost, uord)
tload = 0.0d0
do kk0 = 1,nunits
   u = uord(kk0)
   tmin = minloc(tload, dim=1)
   own(ui(u),uj(u)) = tmin
   tload(tmin) = tload(tmin) + ucost(u)
enddo
deallocate(ui, uj, ucost, uord)
!$omp parallel default(shared) num_threads(nth) &
!$omp&   private(ish,jsh,ksh,di,dj,dk,p,q,r,a,atomI,atomJ,atomK,offi,offj,offP, &
!$omp&           shls,shlsw,buf,Ppq,ns,tid,fac,v)
tid = omp_get_thread_num() + 1
if (omp_get_num_threads() .ne. nth) error stop 'integrals_df_j_fock_deriv: thread team smaller than planned'
do ish = 0,nBases-1
   di = cgto_engine(ish, basDF); atomI = basDF(1,ish+1) + 1; offi = ao_offset0(ish)
   do jsh = 0,ish
      fac = merge(2.0d0, 1.0d0, jsh .ne. ish)
      dj = cgto_engine(jsh, basDF); atomJ = basDF(1,jsh+1) + 1; offj = ao_offset0(jsh)
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      if (own(ish,jsh) .ne. tid) cycle
      do ksh = nBases,nshb-1
         dk = cgto_engine(ksh, basDF); atomK = basDF(1,ksh+1) + 1; offP = aux_offset0(ksh)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) &
             *dmax_shell(ish,jsh) .lt. density_screen_cutoff) cycle
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         shlsw(1) = jsh; shlsw(2) = ish; shlsw(3) = ksh
         allocate(buf(di*dj*dk*3))
         call threec2e_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_base)
         do r = 1,dk; do q = 1,dj; do p = 1,di
            Ppq = Ptot(offi+p,offj+q)
            ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            bvec_t(offP+r,tid) = bvec_t(offP+r,tid) + buf(((r-1)*dj+(q-1))*di+p)*Ppq*ns*fac
         enddo; enddo; enddo
         call threec2e_ip1_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip1)
         do r = 1,dk; do q = 1,dj; do p = 1,di
            Ppq = Ptot(offi+p,offj+q)
            ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            do a = 1,3
               v = buf(((r-1)*dj+(q-1))*di+p+(a-1)*di*dj*dk)*Ppq*ns*fac
               bx_t(atomI,a,offP+r,tid) = bx_t(atomI,a,offP+r,tid) + v
               bx_t(atomJ,a,offP+r,tid) = bx_t(atomJ,a,offP+r,tid) - v
            enddo
         enddo; enddo; enddo
         call threec2e_ip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip2)
         do r = 1,dk; do q = 1,dj; do p = 1,di
            Ppq = Ptot(offi+p,offj+q)
            ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            do a = 1,3
               v = buf(((r-1)*dj+(q-1))*di+p+(a-1)*di*dj*dk)*Ppq*ns*fac
               bx_t(atomK,a,offP+r,tid) = bx_t(atomK,a,offP+r,tid) + v
               bx_t(atomJ,a,offP+r,tid) = bx_t(atomJ,a,offP+r,tid) - v
            enddo
         enddo; enddo; enddo
         deallocate(buf)
      enddo
   enddo
enddo
!$omp end parallel
deallocate(own, tload)
bvec = 0.0d0; bx = 0.0d0
do tid = 1,nth
   bvec = bvec + bvec_t(:,tid)
   bx = bx + bx_t(:,:,:,tid)
enddo
call df_metric_inverse(df_evec, df_evalinv, nContsAux, bvec, cvec)
tj(1) = omp_get_wtime()-tj0 ; tj0 = omp_get_wtime()

call cintdel_optimizer(opt_ip1)
call cintdel_optimizer(opt_ip2)
call twoc2e_ip1_optimizer_engine(opt_ip1, atm, size(atm,2), basDF, nshb, envDF)
call twoc2e_ip2_optimizer_engine(opt_ip2, atm, size(atm,2), basDF, nshb, envDF)
bc_t = 0.0d0
!$omp parallel default(shared) private(ish,jsh,di,dj,p,q,a,atomI,atomJ,offP,offj,shls,buf,ns,crow,tid)
tid = omp_get_thread_num() + 1
!$omp do schedule(static,1)
do ish = nBases,nshb-1
   di = cgto_engine(ish, basDF); atomI = basDF(1,ish+1) + 1; offP = aux_offset0(ish)
   do jsh = nBases,nshb-1
      dj = cgto_engine(jsh, basDF); atomJ = basDF(1,jsh+1) + 1; offj = aux_offset0(jsh)
      if (aux_shell_bound(ish-nBases+1)*aux_shell_bound(jsh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
      shls(1) = ish; shls(2) = jsh
      allocate(buf(di*dj*3))
      call twoc2e_ip1_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip1)
      do q = 1,dj; do p = 1,di
         ns = NorVECAux(offP+p)*NorVECAux(offj+q)
         crow = cvec(offj+q)*ns
         do a = 1,3
            bc_t(atomI,a,offP+p,tid) = bc_t(atomI,a,offP+p,tid) + buf(((q-1)*di+p)+(a-1)*di*dj)*crow
         enddo
      enddo; enddo
      call twoc2e_ip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip2)
      do q = 1,dj; do p = 1,di
         ns = NorVECAux(offP+p)*NorVECAux(offj+q)
         crow = cvec(offj+q)*ns
         do a = 1,3
            bc_t(atomJ,a,offP+p,tid) = bc_t(atomJ,a,offP+p,tid) + buf(((q-1)*di+p)+(a-1)*di*dj)*crow
         enddo
      enddo; enddo
      deallocate(buf)
   enddo
enddo
!$omp end do
!$omp end parallel
bc = 0.0d0
do tid = 1,nth
   bc = bc + bc_t(:,:,:,tid)
enddo
call cintdel_optimizer(opt_ip1)
call cintdel_optimizer(opt_ip2)
deallocate(bvec_t, bx_t, bc_t)

allocate(uT(N3,nContsAux), wtmp(nContsAux))
do atomI = 1,natoms
   do a = 1,3
      wtmp = bx(atomI,a,:) - bc(atomI,a,:)
      call df_metric_inverse(df_evec, df_evalinv, nContsAux, wtmp, bvec)
      uT(a+3*(atomI-1),:) = bvec
   enddo
enddo
deallocate(wtmp)

tj(2) = omp_get_wtime()-tj0 ; tj0 = omp_get_wtime()
maxdo = 0
do ish = 0,nBases-1
   maxdo = max(maxdo, cgto_engine(ish, basDF))
enddo
maxda = 0
do ksh = nBases,nshb-1
   maxda = max(maxda, cgto_engine(ksh, basDF))
enddo
Jx = 0.0d0
call threec2e_ip1_optimizer_engine(opt_ip1, atm, size(atm,2), basDF, nshb, envDF)
call threec2e_ip2_optimizer_engine(opt_ip2, atm, size(atm,2), basDF, nshb, envDF)
!$omp parallel default(shared) &
!$omp&   private(ish,jsh,ksh,di,dj,dk,p,q,r,a,atomK,offi,offj,offP,shls,buf,ns,v,mu,nu,T,S12,atomI,atomJ, &
!$omp&           bnP,uTg,ncol)
allocate(buf(maxdo*maxdo*maxda*3), T(maxdo*maxdo,N3), S12(maxdo*maxdo,3,2))
allocate(bnP(maxdo*maxdo,nContsAux), uTg(N3,nContsAux))
!$omp do schedule(dynamic)
do ish = nBases-1,0,-1
   di = cgto_engine(ish, basDF); offi = ao_offset0(ish); atomI = basDF(1,ish+1) + 1
   do jsh = 0,ish
      dj = cgto_engine(jsh, basDF); offj = ao_offset0(jsh); atomJ = basDF(1,jsh+1) + 1
      if (schwarz_bound(ish,jsh) .lt. SCHWARZ_CUTOFF) cycle
      S12(1:di*dj,:,:) = 0.0d0
      ncol = 0
      do ksh = nBases,nshb-1
         dk = cgto_engine(ksh, basDF); atomK = basDF(1,ksh+1) + 1; offP = aux_offset0(ksh)
         if (schwarz_bound(ish,jsh)*aux_shell_bound(ksh-nBases+1) .lt. SCHWARZ_CUTOFF) cycle
         shls(1) = ish; shls(2) = jsh; shls(3) = ksh
         call threec2e_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_base)
         do r = 1,dk; do q = 1,dj; do p = 1,di
            ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
            bnP((q-1)*di+p,ncol+r) = buf(((r-1)*dj+(q-1))*di+p)*ns
         enddo; enddo; enddo
         uTg(:,ncol+1:ncol+dk) = uT(:,offP+1:offP+dk)
         ncol = ncol + dk
         call threec2e_ip1_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip1)
         do a = 1,3
            do r = 1,dk; do q = 1,dj; do p = 1,di
               ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
               S12((q-1)*di+p,a,1) = S12((q-1)*di+p,a,1) &
                  + buf(((r-1)*dj+(q-1))*di+p+(a-1)*di*dj*dk)*ns*cvec(offP+r)
            enddo; enddo; enddo
         enddo
         call threec2e_ip2_engine(buf, shls, atm, size(atm,2), basDF, nshb, envDF, opt_ip2)
         do a = 1,3
            do r = 1,dk; do q = 1,dj; do p = 1,di
               ns = NorVEC(offi+p)*NorVEC(offj+q)*NorVECAux(offP+r)
               v = buf(((r-1)*dj+(q-1))*di+p+(a-1)*di*dj*dk)*ns*cvec(offP+r)
               Jx(offi+p,offj+q,a,atomK) = Jx(offi+p,offj+q,a,atomK) - v
               S12((q-1)*di+p,a,2) = S12((q-1)*di+p,a,2) + v
            enddo; enddo; enddo
         enddo
      enddo
      if (ncol .gt. 0) then
         call dgemm('N','T',di*dj,N3,ncol,1.0d0,bnP,maxdo*maxdo,uTg,N3,0.0d0,T,maxdo*maxdo)
      else
         T(1:di*dj,:) = 0.0d0
      endif
      do a = 1,3
         do q = 1,dj
            do p = 1,di
               Jx(offi+p,offj+q,a,atomI) = Jx(offi+p,offj+q,a,atomI) - S12((q-1)*di+p,a,1)
               Jx(offi+p,offj+q,a,atomJ) = Jx(offi+p,offj+q,a,atomJ) &
                  + S12((q-1)*di+p,a,1) + S12((q-1)*di+p,a,2)
            enddo
         enddo
      enddo
      do a = 1,natoms
         do r = 1,3
            do q = 1,dj
               do p = 1,di
                  Jx(offi+p,offj+q,r,a) = Jx(offi+p,offj+q,r,a) - T((q-1)*di+p, r+3*(a-1))
               enddo
            enddo
         enddo
      enddo
   enddo
enddo
!$omp end do
deallocate(buf, T, S12, bnP, uTg)
!$omp end parallel
tj(3) = omp_get_wtime()-tj0 ; tj0 = omp_get_wtime()
!$omp parallel do private(a,r,nu,mu) collapse(2)
do a = 1,natoms
   do r = 1,3
      do nu = 1,nConts
         do mu = 1,nConts
            if (ao_shell0(mu) .lt. ao_shell0(nu)) Jx(mu,nu,r,a) = Jx(nu,mu,r,a)
         enddo
      enddo
   enddo
enddo
!$omp end parallel do
call cintdel_optimizer(opt_base)
call cintdel_optimizer(opt_ip1)
call cintdel_optimizer(opt_ip2)
deallocate(uT, bvec, cvec, bx, bc, dmax_shell, ao_offset0, aux_offset0, ao_atom, ao_shell0)
tj(4) = omp_get_wtime()-tj0
block
   character(8) :: ev
   call get_environment_variable("ENGINE_DFJ_VERBOSE", ev)
   if (trim(ev) .eq. "1") print '(A,4F8.2)', ' DFJ FOCK timing: pass1 / metric+u / Jx pass / mirror:', tj(1:4)
end block
call prof_stop("hess_df_j_fock")
end subroutine integrals_df_j_fock_deriv

module subroutine hessian_dfj_read_density(nConts, fname, P)
implicit none
integer,intent(in) :: nConts
character(len=*),intent(in) :: fname
real(8),intent(out) :: P(nConts,nConts)
integer :: u, i, j
open(newunit=u, file=fname, status='old', action='read')
do i = 1,nConts
   read(u,*) (P(i,j), j = 1,nConts)
enddo
close(u)
end subroutine hessian_dfj_read_density

subroutine sort_desc_idx(n, key, idx)
implicit none
integer,intent(in) :: n
real(8),intent(in) :: key(n)
integer,intent(inout) :: idx(n)
integer :: gap, i, j, t
gap = n/2
do while (gap .gt. 0)
   do i = gap+1, n
      t = idx(i)
      j = i
      do while (j .gt. gap)
         if (key(idx(j-gap)) .gt. key(t) .or. (key(idx(j-gap)) .eq. key(t) .and. idx(j-gap) .lt. t)) exit
         idx(j) = idx(j-gap)
         j = j - gap
      enddo
      idx(j) = t
   enddo
   gap = gap/2
enddo
end subroutine sort_desc_idx

subroutine contract9_3c(buf, di, dj, dk, coef, s9)
implicit none
integer,intent(in) :: di,dj,dk
real(8),intent(in) :: buf(di*dj*dk*9)
real(8),intent(in) :: coef(di,dj,dk)
real(8),intent(out) :: s9(9)
integer :: ic,p,q,r,nq,cr
real(8) :: s
nq = di*dj*dk
do ic = 1,9
   s = 0.0d0
   do r = 1,dk
      do q = 1,dj
         do p = 1,di
            cr = ((r-1)*dj+(q-1))*di + p
            s = s + buf(cr + (ic-1)*nq)*coef(p,q,r)
         enddo
      enddo
   enddo
   s9(ic) = s
enddo
end subroutine contract9_3c

subroutine contract9_3c_sw(buf, di, dj, dk, coef, s9)
implicit none
integer,intent(in) :: di,dj,dk
real(8),intent(in) :: buf(di*dj*dk*9)
real(8),intent(in) :: coef(di,dj,dk)
real(8),intent(out) :: s9(9)
integer :: ic,p,q,r,nq,cr
real(8) :: s
nq = di*dj*dk
do ic = 1,9
   s = 0.0d0
   do r = 1,dk
      do p = 1,di
         do q = 1,dj
            cr = ((r-1)*di+(p-1))*dj + q
            s = s + buf(cr + (ic-1)*nq)*coef(p,q,r)
         enddo
      enddo
   enddo
   s9(ic) = s
enddo
end subroutine contract9_3c_sw

subroutine contract9_2c(buf, di, dj, coef, s9)
implicit none
integer,intent(in) :: di,dj
real(8),intent(in) :: buf(di*dj*9)
real(8),intent(in) :: coef(di,dj,1)
real(8),intent(out) :: s9(9)
integer :: ic,p,q,nq,cr
real(8) :: s
nq = di*dj
do ic = 1,9
   s = 0.0d0
   do q = 1,dj
      do p = 1,di
         cr = (q-1)*di + p
         s = s + buf(cr + (ic-1)*nq)*coef(p,q,1)
      enddo
   enddo
   s9(ic) = s
enddo
end subroutine contract9_2c

subroutine contract9_2c_sw(buf, di, dj, coef, s9)
implicit none
integer,intent(in) :: di,dj
real(8),intent(in) :: buf(di*dj*9)
real(8),intent(in) :: coef(di,dj,1)
real(8),intent(out) :: s9(9)
integer :: ic,p,q,nq,cr
real(8) :: s
nq = di*dj
do ic = 1,9
   s = 0.0d0
   do p = 1,di
      do q = 1,dj
         cr = (p-1)*dj + q
         s = s + buf(cr + (ic-1)*nq)*coef(p,q,1)
      enddo
   enddo
   s9(ic) = s
enddo
end subroutine contract9_2c_sw

subroutine df_metric_inverse(evec, evalinv, naux, vin, vout)
implicit none
integer,intent(in) :: naux
real(8),intent(in) :: evec(naux,naux), evalinv(naux)
real(8),intent(in) :: vin(naux)
real(8),intent(out) :: vout(naux)
real(8),allocatable :: w(:)
integer :: i
allocate(w(naux))
w = matmul(transpose(evec), vin)
do i = 1,naux
   w(i) = w(i)*evalinv(i)
enddo
vout = matmul(evec, w)
deallocate(w)
end subroutine df_metric_inverse

end submodule hessian_df_j_impl
