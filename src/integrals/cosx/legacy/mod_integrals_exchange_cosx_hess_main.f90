! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic Hessian - the main grid pass.

submodule (mod_integrals:exchange_cosx_hess_impl) exchange_cosx_hess_main_impl
implicit none
contains

module subroutine hess_main_pass()
real(8),allocatable :: xb(:,:), x1b(:,:,:), x2b(:,:,:,:), wv(:)
real(8),allocatable :: bv(:,:,:), av(:,:,:), Gv(:,:,:), Gtv(:,:,:), hv(:,:,:), htv(:,:,:), Psi(:,:,:)
real(8),allocatable :: b1(:,:,:,:), a1(:,:,:,:), kv(:,:,:,:), mv(:,:,:,:), kp(:,:,:,:), mp(:,:,:,:)
real(8),allocatable :: epsy(:,:), wa(:,:,:), Hp(:,:,:,:), Hpp(:,:,:,:), wh(:,:)
real(8),allocatable :: Hloc(:,:,:,:), gloc(:,:), Vloc(:,:,:), Zloc(:,:,:,:)
real(8),allocatable :: env_g(:), buf(:), dbuf(:), dbuf2(:), dbpp(:), dbvp(:), dbpp2(:), int_cache(:)
real(8),allocatable :: Aval(:), Abra(:,:), Aket(:,:), th(:)
real(8) :: x_block_max, k_bound, esp_eff, r_far, seg_ab(3), seg_cl(3), seg_len2, seg_t, nij
real(8) :: sbb(3,3), sbk(3,3), skk(3,3), ipp, ivp, ivpt, ikk, tt
integer :: ib, nb, bstart, g, si, sj, sj_idx, di, dj, i, j, ao_i, ao_j, atom_i, atom_j
integer :: shls1e(4), off_grids, nloc, loc(natoms), yy, yi, yj, dd, ee, ob, okk, pg
real(8),allocatable :: Fmax(:), Tmax(:), sqa(:), Fu(:), Tu(:), x1t(:,:,:)
logical :: do_d2
real(8) :: dfac, Lw, kx_bound, m1, sth, hdiff, hmag
real(8),allocatable :: h34(:), cbuf(:)
integer :: nq, ii, jj
logical :: fused, hchk, hdiag
integer :: d2gate
character(len=32) :: envs2
real(8),allocatable :: Rg(:,:,:), Og(:,:,:)
real(8),allocatable,target :: Abuf(:)
real(8),pointer :: Asub(:,:,:)
integer :: g0, ng, nbs, mrhs, gb, gp, ngp
integer,parameter :: NGRP = 8
integer :: ip_
real(8),allocatable :: qv(:,:), q1(:,:,:), Qx(:,:), wq(:,:), Kloc(:,:,:,:)
real(8),allocatable :: GC(:,:), WC(:,:), QC(:,:), KC(:,:)
real(8),allocatable :: QVc(:,:), Q1C(:,:,:), KCa(:,:,:), RR(:,:,:,:), RS(:,:), wGC(:,:), XP(:), XWloc(:,:,:)
real(8) :: xmx
real(8),allocatable :: kav(:,:,:,:), xtg(:,:)
real(8),allocatable :: A1C(:,:,:,:), KCs(:,:,:,:), KT(:,:,:), GCt(:,:), anrm(:,:,:), XKloc(:,:,:,:)
real(8),allocatable :: aou(:,:), aow(:,:), abu(:), abw(:)
real(8) :: Lpt(BATCH)
real(8) :: tpre, tint, tloop, tpost, tq1, tq2, tdns, tq3, ts(10), tsq, sloc, spar, snb, ssig(4), te1, te2, cn1, cn2, tei
real(8) :: afill, nsubb, arow, tkx(6), sact(4), rsig(4,4), rblk
logical :: scrd
real(8),allocatable :: Vp(:,:), Fp(:,:), Up(:,:), Fup(:,:), Lp(:), amx(:), a1mx(:)
real(8) :: scr(24), S2(3,3,3), colthr, tua, ncol_ar, nptall, P2, bpt, bh, bk, bhs, wgt
logical,allocatable :: csig(:,:)
integer :: clist(2*N3), nb0, na0, c, oo
real(8),allocatable :: cmb(:,:), cma(:,:), sw(:)
real(8) :: fw, tw, fuu, tuu
real(8),allocatable :: XPa(:,:,:)
integer :: nbh, nch, ich, cb
integer :: nseg, s0, nsg
logical :: acmp
integer :: acap, nsa, nbsc, ic, jc
logical :: kcomp
integer :: kvcap, kg0
logical :: bcomp
real(8),allocatable :: B1C(:,:,:,:), Tb(:,:,:)
logical :: bfly
logical :: psub
real(8),allocatable :: Tsub(:,:,:,:)
real(8),allocatable :: Csl(:,:,:), Bg(:,:), bnm(:,:,:)
real(8) :: cms(2)
integer :: ncd, cand(2*N3)
real(8),allocatable :: KCc(:,:,:,:), MCc(:,:,:,:), XAv(:,:,:), Tkc(:,:)
logical,allocatable :: pok(:), pd2(:), shon(:)
integer,allocatable :: aoc(:), slist(:)
integer,allocatable :: ch_ib(:), ch_st(:), ch_nb(:)

tpre = 0.0d0 ; tint = 0.0d0 ; tloop = 0.0d0 ; tpost = 0.0d0 ; tdns = 0.0d0 ; ts = 0.0d0
sloc = 0.0d0 ; spar = 0.0d0 ; snb = 0.0d0 ; ssig = 0.0d0
tkx = 0.0d0 ; sact = 0.0d0
te1 = 0.0d0 ; te2 = 0.0d0 ; cn1 = 0.0d0 ; cn2 = 0.0d0 ; afill = 0.0d0 ; nsubb = 0.0d0 ; arow = 0.0d0
fused = .false.
call get_environment_variable("ENGINE_COSX_HESS_FUSED", envs2)
if (trim(envs2) .eq. "1") fused = .true.
nseg = max(1, min(BATCH, 4000000/(n*n)))
call get_environment_variable("ENGINE_COSX_HESS_D2SEG", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) nseg
if (fused) nseg = 0
if (nseg .gt. 0) then
   nbs = max(1, min(BATCH, 16000000/(n*n)))
   nbs = max(nseg, (nbs/nseg)*nseg)
else
   nbs = max(1, min(BATCH, 4000000/(n*n)))
endif
call get_environment_variable("ENGINE_COSX_HESS34_CHECK", envs2)
hchk = (trim(envs2) .eq. "1")
call get_environment_variable("ENGINE_COSX_HESS_DIAG", envs2)
hdiag = (trim(envs2) .eq. "1")
hdiff = 0.0d0 ; hmag = 0.0d0
scr = 0.0d0 ; ncol_ar = 0.0d0 ; nptall = 0.0d0 ; rsig = 0.0d0 ; rblk = 0.0d0
colthr = 1.0d-7
call get_environment_variable("ENGINE_COSX_HESS_COLTHR", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) colthr
t3lr = lowrank .and. kproj
call get_environment_variable("ENGINE_COSX_HESS_T3LR", envs2)
if (trim(envs2) .eq. "0") t3lr = .false.
noa1 = t3lr .and. .not. var_k .and. .not. hdiag .and. .not. scrd
call get_environment_variable("ENGINE_COSX_HESS_NOA1", envs2)
if (trim(envs2) .eq. "0") noa1 = .false.
kxx = kproj .and. .not. kxmir
abnd = .false.
call get_environment_variable("ENGINE_COSX_HESS_ABOUND", envs2)
if (trim(envs2) .eq. "1") abnd = .true.
apr = noa1 .and. .not. abnd .and. want_kx
call get_environment_variable("ENGINE_COSX_HESS_APRUNE", envs2)
if (trim(envs2) .eq. "0") apr = .false.
call get_environment_variable("ENGINE_COSX_HESS_KXX", envs2)
if (trim(envs2) .eq. "0") kxx = .false.
if (noa1) then
   allocate(cnrm(n,ns), ao_sh(n))
   do si = 0,nBases-1
      ao_sh(cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)) = si
   enddo
   do is = 1,ns
      do i = 1,n
         cnrm(i,is) = sqrt(sum(Cp(i,1:nop,is)**2))
      enddo
   enddo
endif
call get_environment_variable("ENGINE_COSX_HESS_SCRDIAG", envs2)
scrd = (trim(envs2) .eq. "1")
d2gate = 1
call get_environment_variable("ENGINE_COSX_HESS_D2GATE", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) d2gate
call get_environment_variable("ENGINE_COSX_HESS_SUBPTS", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) nbs
acmp = (nseg .gt. 0) .and. .not. scrd .and. .not. hdiag
call get_environment_variable("ENGINE_COSX_HESS_ACMP", envs2)
if (trim(envs2) .eq. "0") acmp = .false.
acap = 8000000
call get_environment_variable("ENGINE_COSX_HESS_ACAP", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) acap
kcomp = lowrank .and. kproj .and. want_kx .and. t3lr .and. (nseg .gt. 0) .and. &
        .not. (var_k .or. kxmir .or. kxd .or. hdiag .or. scrd)
call get_environment_variable("ENGINE_COSX_HESS_KCOMP", envs2)
if (trim(envs2) .eq. "0") kcomp = .false.
kvcap = 8000000
call get_environment_variable("ENGINE_COSX_HESS_KVCAP", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) kvcap
bcomp = .false.
call get_environment_variable("ENGINE_COSX_HESS_BCOMP", envs2)
if (trim(envs2) .eq. "1") bcomp = kcomp .and. apr
bfly = kcomp .and. apr .and. acmp .and. .not. bcomp
call get_environment_variable("ENGINE_COSX_HESS_BFLY", envs2)
if (trim(envs2) .eq. "0") bfly = .false.
allocate(sqa(0:nBases-1))
do si = 0,nBases-1
   sqa(si) = 0.0d0
   do ip_ = 1,bas(3,si+1)
      sqa(si) = max(sqa(si), sqrt(env(bas(6,si+1)+ip_)))
   enddo
enddo
nbh = BATCH
call get_environment_variable("ENGINE_COSX_HESS_BLKPTS", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) nbh
nch = 0
do ib = 0,nblocks-1
   nch = nch + (cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib) + nbh - 1)/nbh
enddo
allocate(ch_ib(nch), ch_st(nch), ch_nb(nch))
nch = 0
do ib = 0,nblocks-1
   do cb = cosx_blkstart_hi(ib), cosx_blkstart_hi(ib+1)-1, nbh
      nch = nch + 1
      ch_ib(nch) = ib
      ch_st(nch) = cb + 1
      ch_nb(nch) = min(nbh, cosx_blkstart_hi(ib+1) - cb)
   enddo
enddo
kxonly = kxx .and. .not. var_k
xkshr = kxx
call get_environment_variable("ENGINE_COSX_HESS_XKSHARED", envs2)
if (trim(envs2) .eq. "0") xkshr = .false.
psub = .false.
call get_environment_variable("ENGINE_COSX_HESS_PSUB", envs2)
if (trim(envs2) .eq. "1") psub = bfly .and. (lrz .or. lrfit) .and. xkshr .and. fit .and. .not. var_k
if (kxx) then
   allocate(XKsum(n,nop,N3,ns))
   !$omp parallel do private(yy) schedule(static)
   do yy = 1,N3
      XKsum(:,:,yy,:) = 0.0d0
   enddo
   !$omp end parallel do
endif
if (xkshr) then
   allocate(xklk((N3+7)/8))
   do i = 1,(N3+7)/8
      call omp_init_lock(xklk(i))
   enddo
endif
if (fit .and. zmode .eq. 1) then
   allocate(zlock(natoms))
   do i = 1,natoms
      call omp_init_lock(zlock(i))
   enddo
endif
!$omp parallel default(shared) &
!$omp&  private(ich) &
!$omp&  private(xb,x1b,x2b,wv,bv,av,Gv,Gtv,hv,htv,Psi,b1,a1,kv,mv,kp,mp,epsy,wa,Hp,Hpp,wh) &
!$omp&  private(Hloc,gloc,Vloc,Zloc,env_g,buf,dbuf,dbuf2,dbpp,dbvp,dbpp2,int_cache,Aval,Abra,Aket,th) &
!$omp&  private(x_block_max,k_bound,esp_eff,r_far,seg_ab,seg_cl,seg_len2,seg_t,nij,sbb,sbk,skk) &
!$omp&  private(ipp,ivp,ivpt,ikk,tt,ib,nb,bstart,g,si,sj,sj_idx,di,dj,i,j,ao_i,ao_j,atom_i,atom_j) &
!$omp&  private(shls1e,off_grids,nloc,loc,yy,yi,yj,dd,ee,ob,okk,pg,is,Fmax,qv,q1,Qx,wq,Kloc,kav,xtg,Tmax,do_d2) &
!$omp&  private(GC,WC,QC,KC,QVc,Q1C,KCa,RR,RS,wGC,XP,XWloc,xmx,A1C,KCs,KT,GCt,anrm,XKloc,aou,aow,abu,abw) &
!$omp&  private(tq1,tq2,dfac,Fu,Tu,x1t,Lw,kx_bound,m1,sth,Asub,Rg,Og,g0,ng,mrhs,tq3,gb,gp,ngp) &
!$omp&  private(tsq,h34,cbuf,nq,ii,jj,tei,s0,nsg,nsa,nbsc,ic,jc,pok,pd2,shon,aoc,slist,Abuf,kg0,KCc,MCc,XAv,Tkc,B1C,Tb) &
!$omp&  private(Csl,Bg,bnm,cms,ncd,cand) &
!$omp&  private(Vp,Fp,Up,Fup,Lp,amx,a1mx,S2,P2,bpt,bh,bk,bhs,wgt,Lpt,csig,clist,nb0,na0,c,tua,cmb,cma,sw,fw,tw,fuu,tuu,XPa,oo) &
!$omp&  reduction(+:tpre,tint,tloop,tpost,tdns,ts,sloc,spar,snb,ssig,te1,te2,cn1,cn2,afill,nsubb,arow,tkx,sact,scr,ncol_ar,nptall) &
!$omp&  reduction(+:rsig,rblk) &
!$omp&  reduction(max:hdiff,hmag)
allocate(Hloc(N3,N3,NT,ns), gloc(N3,ns))
if (want_kx) then
   allocate(Kloc(n,nop,N3,ns))
   if (.not. kxonly) Kloc = 0.0d0
endif
if (kxx .and. .not. xkshr) then
   allocate(XKloc(n,nop,N3,ns))
   XKloc = 0.0d0
else
   allocate(XKloc(1,1,1,1))
endif
if (kproj .and. fit) then
   allocate(XWloc(n,nop,ns))
else
   allocate(XWloc(1,1,1))
endif
XWloc = 0.0d0
Hloc = 0.0d0
gloc = 0.0d0
if (fit) then
   if (lrfit .or. lrz) then
      allocate(Vloc(1,1,1), Zloc(n,nop,3,ns))
   else if (zt7) then
      allocate(Vloc(n,n,ns), Zloc(1,1,1,ns))
   else if (zmode .eq. 1) then
      allocate(Vloc(n,n,ns), Zloc(n,n,3,ns))
   else
      allocate(Vloc(n,n,ns), Zloc(n,n,N3,ns))
   endif
   Vloc = 0.0d0
   Zloc = 0.0d0
endif
allocate(buf(BATCH*maxd*maxd), dbuf(BATCH*maxd*maxd*3), dbuf2(BATCH*maxd*maxd*3))
allocate(dbpp(BATCH*maxd*maxd*9), dbvp(BATCH*maxd*maxd*9), dbpp2(BATCH*maxd*maxd*9))
allocate(int_cache(INT_CACHE_SIZE), Fmax(0:nBases-1), Tmax(0:nBases-1))
allocate(Fu(0:nBases-1), Tu(0:nBases-1))
mrhs = 2 + 2*N3
if (acmp) then
   allocate(Abuf(max(acap, n*n*nseg)), pok(cosx_ext_ptr(nBases)), pd2(cosx_ext_ptr(nBases)), &
            shon(0:nBases-1), aoc(n), slist(n))
else
   allocate(Abuf(n*n*nbs))
   Asub(1:n,1:n,1:nbs) => Abuf
endif
allocate(Rg(n,mrhs,NGRP), Og(n,mrhs,NGRP))
if (bfly) allocate(Csl(n,nop,ns), Bg(nop,N3), bnm(BATCH,N3,ns))
allocate(h34(BATCH*maxd*maxd*34), cbuf(BATCH*maxd*maxd*9))
allocate(env_g(size(env) + 3*BATCH))
env_g(1:size(env)) = env(1:size(env))
off_grids = size(env)
env_g(9) = hs_rsomega

!$omp do schedule(dynamic)
do ich = 1,nch
   tq1 = omp_get_wtime()
   ib = ch_ib(ich)
   bstart = ch_st(ich)
   nb = ch_nb(ich)
   allocate(xb(n,nb), x1b(n,3,nb), x2b(n,3,3,nb), wv(nb))
   tsq = omp_get_wtime()
   call hess_block_ao(ib, nb, bstart, xb, x1b, x2b)
   ts(1) = ts(1) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
   call hess_block_atoms(ib, nloc, loc)
   if (hdiag) then
   block
      logical :: inl(natoms), pseen(natoms)
      integer :: gq
      inl = .false. ; pseen = .false.
      inl(loc(1:nloc)) = .true.
      do gq = 1,nb
         if (.not. inl(cosx_atom_hi(bstart+gq-1)) .and. .not. pseen(cosx_atom_hi(bstart+gq-1))) then
            pseen(cosx_atom_hi(bstart+gq-1)) = .true.
            spar = spar + nb
         endif
      enddo
      sloc = sloc + dble(nloc)*nb
      snb = snb + nb
      block
         integer :: kq, it4
         real(8) :: vm, thr4(4)
         thr4 = (/1.0d-6, 1.0d-8, 1.0d-10, 1.0d-12/)
         block
            integer :: aa, it5
            real(8) :: am
            do aa = 1,natoms
               am = 0.0d0
               do gq = 1,nb
                  am = max(am, maxval(abs(x1b(at_off(aa)+1:at_off(aa)+at_n(aa),:,gq))), &
                               maxval(abs(xb(at_off(aa)+1:at_off(aa)+at_n(aa),gq))))
               enddo
               do it5 = 1,4
                  if (am .ge. 10.0d0**(-5-2*(it5-1))) sact(it5) = sact(it5) + nb
               enddo
            enddo
         end block
         do kq = 1,n
            vm = max(maxval(abs(xb(kq,:))), maxval(abs(x1b(kq,:,:))), maxval(abs(x2b(kq,:,:,:))))
            do it4 = 1,4
               if (vm .ge. thr4(it4)) ssig(it4) = ssig(it4) + nb
            enddo
         enddo
      end block
   end block
   endif
   wv = cosx_w_hi(bstart:bstart+nb-1)
   nptall = nptall + nb

   allocate(bv(nb,n,ns), av(nb,n,ns), Gv(nb,n,ns), Gtv(nb,n,ns))
   if (apr) then
      allocate(a1(1,1,1,ns))
   else
      allocate(a1(nb,N3,n,ns))
   endif
   if (bcomp) then
      allocate(B1C(nb,nop,N3,ns))
   else if (bfly) then
      allocate(B1C(nb,nop,N3,ns), b1(1,1,1,1))
   else
      allocate(b1(nb,N3,n,ns))
   endif
   if (kcomp .and. psub) then
      allocate(XAv(nb,nop,ns))
   else if (kcomp) then
      allocate(KCc(nb,N3,nop,ns), MCc(nb,N3,nop,ns), XAv(nb,nop,ns))
   else
      allocate(kv(nb,N3,n,ns), mv(nb,N3,n,ns))
   endif
   allocate(kp(nb,n,3,ns), mp(nb,n,3,ns))
   if (t3lr) allocate(A1C(nb,nop,N3,ns))
   do is = 1,ns
      if (lowrank) then
         allocate(XPa(nb,nop,1))
         call dgemm('T','N',nb,nop,n,1.0d0,xb,n,Cp(1,1,is),n,0.0d0,XPa,nb)
         call dgemm('N','T',nb,n,nop,1.0d0,XPa,nb,Cp(1,1,is),n,0.0d0,bv(:,:,is),nb)
         call dgemm('T','N',nb,nop,n,1.0d0,xb,n,PtCp(1,1,is),n,0.0d0,XPa,nb)
         call dgemm('N','T',nb,n,nop,1.0d0,XPa,nb,Cp(1,1,is),n,0.0d0,av(:,:,is),nb)
         if (kcomp) XAv(:,:,is) = XPa(:,:,1)
         deallocate(XPa)
         if (bcomp .or. bfly) then
            B1C(:,:,:,is) = 0.0d0
            call hess_tangent_occ(nb, bstart, nloc, loc, x1b, Cp(:,:,is), B1C(:,:,:,is))
         else
            call hess_tangent_lr(nb, bstart, nloc, loc, x1b, Dsp(:,:,is), Cp(:,:,is), Cp(:,:,is), b1(:,:,:,is))
         endif
         if (.not. noa1) call hess_tangent_lr(nb, bstart, nloc, loc, x1b, DPt(:,:,is), PtCp(:,:,is), Cp(:,:,is), a1(:,:,:,is))
      else
         call dgemm('T','N',nb,n,n,1.0d0,xb,n,Dsp(:,:,is),n,0.0d0,bv(:,:,is),nb)
         call dgemm('T','N',nb,n,n,1.0d0,xb,n,DPt(:,:,is),n,0.0d0,av(:,:,is),nb)
         call hess_tangent(nb, bstart, nloc, loc, x1b, Dsp(:,:,is), b1(:,:,:,is))
         call hess_tangent(nb, bstart, nloc, loc, x1b, DPt(:,:,is), a1(:,:,:,is))
      endif
      ts(2) = ts(2) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
      if (fit .and. lowrank) then
         allocate(XPa(nb,nop,N3))
         call dgemm('T','N',nb,nop*N3,n,1.0d0,xb,n,PtC(1,1,1,is),n,0.0d0,XPa,nb)
         if (.not. noa1) then
            do yy = 1,N3
               call dgemm('N','T',nb,n,nop,1.0d0,XPa(1,1,yy),nb,Cp(1,1,is),n,1.0d0,a1(1,yy,1,is),nb*N3)
            enddo
         endif
         if (t3lr) A1C(:,:,:,is) = XPa
         deallocate(XPa)
      else if (fit) then
         do yy = 1,N3
            call dgemm('T','T',nb,n,n,1.0d0,xb,n,DPy(:,:,yy,is),n,1.0d0,a1(1,yy,1,is),nb*N3)
         enddo
      endif
      if (t3lr) then
         if (.not. fit) A1C(:,:,:,is) = 0.0d0
         call hess_tangent_occ(nb, bstart, nloc, loc, x1b, PtCp(:,:,is), A1C(:,:,:,is))
         if (noa1 .and. .not. abnd .and. .not. apr) then
            allocate(XPa(nb,N3,nop))
            do oo = 1,nop
               do yy = 1,N3
                  XPa(:,yy,oo) = A1C(:,oo,yy,is)
               enddo
            enddo
            call dgemm('N','T',nb*N3,n,nop,1.0d0,XPa,nb*N3,Cp(1,1,is),n,0.0d0,a1(1,1,1,is),nb*N3)
            deallocate(XPa)
         endif
      endif
      ts(3) = ts(3) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
   enddo
   Gv = 0.0d0 ; Gtv = 0.0d0 ; kp = 0.0d0 ; mp = 0.0d0
   if (.not. kcomp) then
      kv = 0.0d0 ; mv = 0.0d0
   endif
   if (want_kx .and. var_k) then
      allocate(kav(nb,N3,n,ns))
      kav = 0.0d0
   endif
   allocate(Hp(nb,N3,3,ns), Hpp(nb,3,3,ns), Aval(nb), Abra(nb,3), Aket(nb,3), th(nb))
   Hp = 0.0d0 ; Hpp = 0.0d0

   allocate(csig(2*N3,nb), cmb(nb,N3), cma(nb,N3), sw(nb))
   sw = sqrt(abs(wv))
   cmb = 0.0d0 ; cma = 0.0d0
   if (noa1 .and. abnd) then
      allocate(anrm(nb,N3,ns), aou(nop,ns), aow(nop,ns), abu(ns), abw(ns))
      aou = 0.0d0 ; aow = 0.0d0 ; abu = 0.0d0 ; abw = 0.0d0
      do is = 1,ns
         do yy = 1,N3
            do oo = 1,nop
               do g = 1,nb
                  aou(oo,is) = max(aou(oo,is), abs(A1C(g,oo,yy,is)))
                  aow(oo,is) = max(aow(oo,is), sw(g)*abs(A1C(g,oo,yy,is)))
               enddo
            enddo
            do g = 1,nb
               anrm(g,yy,is) = sqrt(sum(A1C(g,1:nop,yy,is)**2))
               abu(is) = max(abu(is), anrm(g,yy,is))
               abw(is) = max(abw(is), sw(g)*anrm(g,yy,is))
            enddo
         enddo
      enddo
   endif
   do si = 0,nBases-1
      fw = 0.0d0 ; tw = 0.0d0 ; fuu = 0.0d0 ; tuu = 0.0d0
      do is = 1,ns
         do i = cosx_ao_offset(si)+1, cosx_ao_offset(si)+cosx_shell_dim(si)
            do g = 1,nb
               tt = max(abs(bv(g,i,is)), abs(av(g,i,is)))
               fuu = max(fuu, tt)
               fw = max(fw, sw(g)*tt)
            enddo
            if (bcomp .or. bfly) then
               continue
            else if ((noa1 .and. abnd) .or. apr) then
            do yy = 1,N3
               do g = 1,nb
                  cmb(g,yy) = max(cmb(g,yy), abs(b1(g,yy,i,is)))
                  tt = abs(b1(g,yy,i,is))
                  tuu = max(tuu, tt)
                  tw = max(tw, sw(g)*tt)
               enddo
            enddo
            if (.not. apr) then
            tuu = max(tuu, min(abu(is)*cnrm(i,is), sum(aou(:,is)*abs(Cp(i,1:nop,is)))))
            tw = max(tw, min(abw(is)*cnrm(i,is), sum(aow(:,is)*abs(Cp(i,1:nop,is)))))
            endif
            else
            do yy = 1,N3
               do g = 1,nb
                  cmb(g,yy) = max(cmb(g,yy), abs(b1(g,yy,i,is)))
                  cma(g,yy) = max(cma(g,yy), abs(a1(g,yy,i,is)))
                  tt = max(abs(b1(g,yy,i,is)), abs(a1(g,yy,i,is)))
                  tuu = max(tuu, tt)
                  tw = max(tw, sw(g)*tt)
               enddo
            enddo
            endif
         enddo
      enddo
      Fmax(si) = fw
      Tmax(si) = max(fw, tw)
      Fu(si) = fuu
      Tu(si) = max(fuu, tuu)
   enddo
   if (apr) then
      do is = 1,ns
         if (bcomp .or. bfly) call a1_shell_max(nb, is, sw, B1C(:,:,:,is), Tu, Tmax)
         call a1_shell_max(nb, is, sw, A1C(:,:,:,is), Tu, Tmax)
      enddo
   endif
   do g = 1,nb
      do yy = 1,N3
         csig(yy,g) = (sw(g)*cmb(g,yy) .ge. colthr)
         csig(N3+yy,g) = (sw(g)*cma(g,yy) .ge. colthr)
      enddo
   enddo
   deallocate(cmb, cma, sw)
   if (allocated(anrm)) deallocate(anrm, aou, aow, abu, abw)
   tsq = omp_get_wtime()
   Lw = 0.0d0
   if (want_kx) then
      allocate(qv(nb,n))
      call dgemm('T','N',nb,n,n,1.0d0,xb,n,Pt,n,0.0d0,qv,nb)
      if (.not. kproj) then
         allocate(q1(nb,N3,n))
         call hess_tangent(nb, bstart, nloc, loc, x1b, Pt, q1)
         if (fit) then
            do yy = 1,N3
               call dgemm('T','T',nb,n,n,1.0d0,xb,n,Py(:,:,yy),n,1.0d0,q1(1,yy,1),nb*N3)
            enddo
         endif
         if (var_k) then
            allocate(x1t(nb,N3,n))
            call hess_tangent(nb, bstart, nloc, loc, x1b, Idn, x1t)
         endif
      endif
      call kx_free_leg_max(nb, bstart, xb, x1b, wv, qv, Lpt)
      Lw = maxval(Lpt)
   endif
   if (scrd) then
      allocate(Vp(nb,0:nBases-1), Fp(nb,0:nBases-1), Up(nb,0:nBases-1), Fup(nb,0:nBases-1), Lp(nb))
      allocate(amx(nb), a1mx(nb))
      Vp = 0.0d0 ; Fp = 0.0d0 ; Up = 0.0d0 ; Fup = 0.0d0 ; Lp = 0.0d0
      do si = 0,nBases-1
         do is = 1,ns
            do i = cosx_ao_offset(si)+1, cosx_ao_offset(si)+cosx_shell_dim(si)
               do g = 1,nb
                  tt = max(abs(bv(g,i,is)), abs(av(g,i,is)))
                  Fup(g,si) = max(Fup(g,si), tt)
                  Fp(g,si) = max(Fp(g,si), sqrt(abs(wv(g)))*tt)
                  tt = max(tt, maxval(abs(b1(g,:,i,is))), maxval(abs(a1(g,:,i,is))))
                  Up(g,si) = max(Up(g,si), tt)
                  Vp(g,si) = max(Vp(g,si), sqrt(abs(wv(g)))*tt)
               enddo
            enddo
         enddo
      enddo
      scr(12) = scr(12) + nb
      scr(13) = scr(13) + nb*cosx_blkrad_hi(ib)
      scr(14) = scr(14) + nb*dble(nloc)
      do g = 1,nb
         do yy = 1,N3
            tt = 0.0d0
            do is = 1,ns
               tt = max(tt, maxval(abs(b1(g,yy,:,is))), maxval(abs(a1(g,yy,:,is))))
            enddo
            tt = sqrt(abs(wv(g)))*tt
            do i = 1,5
               if (tt .ge. 10.0d0**(-3-i)) scr(14+i) = scr(14+i) + 1.0d0
            enddo
         enddo
      enddo
      if (want_kx) Lp(1:nb) = Lpt(1:nb)
   endif
   do g = 1,nb
      env_g(off_grids+3*(g-1)+1:off_grids+3*g) = cosx_coor_hi(:,bstart+g-1)
   enddo
   env_g(12) = dble(nb)
   env_g(13) = dble(off_grids)
   ts(4) = ts(4) + (omp_get_wtime() - tsq)
   tq2 = omp_get_wtime() ; tpre = tpre + (tq2 - tq1) ; tq1 = tq2

   nbsc = nbs
   nsa = n
   if (acmp) then
      shon = .false.
      do si = 0,nBases-1
         do sj_idx = cosx_ext_ptr(si)+1, cosx_ext_ptr(si+1)
            sj = cosx_ext_idx(sj_idx)
            pok(sj_idx) = .false.
            if (schwarz_bound(si,sj) .lt. SCHWARZ_CUTOFF) cycle
            seg_ab = cosx_shell_cen(:,sj) - cosx_shell_cen(:,si)
            seg_len2 = sum(seg_ab*seg_ab)
            if (seg_len2 .gt. 0.0d0) then
               seg_t = sum((cosx_blkcen_hi(:,ib) - cosx_shell_cen(:,si))*seg_ab)/seg_len2
               seg_t = max(0.0d0, min(1.0d0, seg_t))
               seg_cl = cosx_shell_cen(:,si) + seg_t*seg_ab
            else
               seg_cl = cosx_shell_cen(:,si)
            endif
            r_far = sqrt(sum((cosx_blkcen_hi(:,ib) - seg_cl)**2)) - cosx_blkrad_hi(ib)
            if (r_far .gt. 0.0d0) then
               esp_eff = min(cosx_esp_bound(si,sj), cosx_esp_mono(si,sj)/r_far)
            else
               esp_eff = cosx_esp_bound(si,sj)
            endif
            dfac = 2.0d0*(sqa(si)+sqa(sj))
            k_bound = esp_eff*max(Tmax(si)*Tmax(sj), dfac*Fmax(si)*Fmax(sj))
            if (want_kx) then
               kx_bound = esp_eff*Lw*max(Tu(sj)*cmx(si), Tu(si)*cmx(sj), &
                                         dfac*Fu(sj)*cmx(si), dfac*Fu(si)*cmx(sj))
               if (k_bound .lt. kscreen_h .and. kx_bound .lt. kscreen_hk) cycle
            else
               if (k_bound .lt. kscreen_h) cycle
            endif
            pok(sj_idx) = .true.
            pd2(sj_idx) = (esp_eff*dfac*dfac*Fmax(si)*Fmax(sj) .ge. kscreen_h2)
            shon(si) = .true.
            shon(sj) = .true.
         enddo
      enddo
      nsa = 0
      aoc = 0
      do si = 0,nBases-1
         if (.not. shon(si)) cycle
         do i = cosx_ao_offset(si)+1, cosx_ao_offset(si)+cosx_shell_dim(si)
            nsa = nsa + 1
            aoc(i) = nsa
            slist(nsa) = i
         enddo
      enddo
      nbsc = min(BATCH, acap/max(1,nsa*nsa))
      nbsc = max(nseg, (nbsc/nseg)*nseg)
      Asub(1:max(nsa,1),1:max(nsa,1),1:nbsc) => Abuf
   endif
   if (bfly) then
      tq3 = omp_get_wtime()
      do is = 1,ns
         cms(is) = 0.0d0
         do jc = 1,nsa
            Csl(jc,:,is) = Cp(slist(jc),1:nop,is)
            cms(is) = max(cms(is), cnrm(slist(jc),is))
         enddo
         do yy = 1,N3
            do g = 1,nb
               bnm(g,yy,is) = sqrt(abs(wv(g)))*sqrt(sum(B1C(g,:,yy,is)**2))*cms(is)
            enddo
         enddo
      enddo
      ts(2) = ts(2) + (omp_get_wtime() - tq3)
   endif
   if (kcomp) then
      if (bcomp) then
         nbsc = min(nbsc, max(nseg, ((kvcap/max(1,3*N3*n*ns))/nseg)*nseg))
         allocate(b1(nbsc,N3,n,ns))
      else
         nbsc = min(nbsc, max(nseg, ((kvcap/max(1,2*N3*n*ns))/nseg)*nseg))
      endif
      allocate(kv(nbsc,N3,n,ns), mv(nbsc,N3,n,ns))
      if (psub) allocate(KCc(nbsc,N3,nop,ns), MCc(nbsc,N3,nop,ns))
   endif
   do g0 = 0, nb-1, nbsc
   ng = min(nbsc, nb-g0)
   Asub(:,:,1:ng) = 0.0d0
   kg0 = g0
   if (kcomp) then
      kg0 = 0
      kv = 0.0d0 ; mv = 0.0d0
   endif
   if (bcomp) then
      tq3 = omp_get_wtime()
      allocate(Tb(ng,N3,n))
      do is = 1,ns
         call hess_tangent_lr(ng, bstart+g0, nloc, loc, x1b(:,:,g0+1:g0+ng), Dsp(:,:,is), Cp(:,:,is), &
                              Cp(:,:,is), Tb)
         b1(1:ng,:,:,is) = Tb
      enddo
      deallocate(Tb)
      do g = g0+1,g0+ng
         tt = sqrt(abs(wv(g)))
         do yy = 1,N3
            csig(yy,g) = (tt*maxval(abs(b1(g-g0,yy,:,:))) .ge. colthr)
            csig(N3+yy,g) = .false.
         enddo
      enddo
      ts(2) = ts(2) + (omp_get_wtime() - tq3)
   endif
   do si = 0,nBases-1
      di = cosx_shell_dim(si)
      atom_i = cosx_shell_atom(si)
      do sj_idx = cosx_ext_ptr(si)+1, cosx_ext_ptr(si+1)
         sj = cosx_ext_idx(sj_idx)
         if (acmp) then
            if (.not. pok(sj_idx)) cycle
            dfac = 2.0d0*(sqa(si)+sqa(sj))
            do_d2 = pd2(sj_idx)
         else
         if (schwarz_bound(si,sj) .lt. SCHWARZ_CUTOFF) cycle
         seg_ab = cosx_shell_cen(:,sj) - cosx_shell_cen(:,si)
         seg_len2 = sum(seg_ab*seg_ab)
         if (seg_len2 .gt. 0.0d0) then
            seg_t = sum((cosx_blkcen_hi(:,ib) - cosx_shell_cen(:,si))*seg_ab)/seg_len2
            seg_t = max(0.0d0, min(1.0d0, seg_t))
            seg_cl = cosx_shell_cen(:,si) + seg_t*seg_ab
         else
            seg_cl = cosx_shell_cen(:,si)
         endif
         r_far = sqrt(sum((cosx_blkcen_hi(:,ib) - seg_cl)**2)) - cosx_blkrad_hi(ib)
         if (r_far .gt. 0.0d0) then
            esp_eff = min(cosx_esp_bound(si,sj), cosx_esp_mono(si,sj)/r_far)
         else
            esp_eff = cosx_esp_bound(si,sj)
         endif
         dfac = 2.0d0*(sqa(si)+sqa(sj))
         k_bound = esp_eff*max(Tmax(si)*Tmax(sj), dfac*Fmax(si)*Fmax(sj))
         if (want_kx) then
            kx_bound = esp_eff*Lw*max(Tu(sj)*cmx(si), Tu(si)*cmx(sj), &
                                      dfac*Fu(sj)*cmx(si), dfac*Fu(si)*cmx(sj))
            if (k_bound .lt. kscreen_h .and. kx_bound .lt. kscreen_hk) cycle
         else
            if (k_bound .lt. kscreen_h) cycle
         endif
         do_d2 = (esp_eff*dfac*dfac*Fmax(si)*Fmax(sj) .ge. kscreen_h2)
         endif
         dj = cosx_shell_dim(sj)
         if (scrd) then
            wgt = dble(ng*di*dj)
            scr(1) = scr(1) + wgt
            if (k_bound .lt. kscreen_h) scr(11) = scr(11) + wgt
            bh = 0.0d0 ; bk = 0.0d0
            do g = g0+1,g0+ng
               bh = max(bh, Vp(g,si)*Vp(g,sj), dfac*max(Vp(g,si)*Fp(g,sj), Fp(g,si)*Vp(g,sj)))
               if (want_kx) bk = max(bk, Lp(g)*max(Up(g,si), Up(g,sj), dfac*Fup(g,si), dfac*Fup(g,sj)))
            enddo
            if (esp_eff*bh .ge. kscreen_h .or. esp_eff*bk .ge. kscreen_hk) scr(2) = scr(2) + wgt
            if (esp_eff*bh .ge. kscreen_h) scr(10) = scr(10) + wgt
         endif
         atom_j = cosx_shell_atom(sj)
         shls1e(1) = si ; shls1e(2) = sj ; shls1e(3) = g0 ; shls1e(4) = g0+ng
         tq2 = omp_get_wtime()
         if (fused) then
            call grids1e_hess34_engine_cached(h34, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            nq = ng*di*dj
            buf(1:nq) = h34(1:nq)
            dbuf(1:3*nq) = h34(nq+1:4*nq)
            dbpp(1:9*nq) = h34(7*nq+1:16*nq)
            dbvp(1:9*nq) = h34(16*nq+1:25*nq)
            do dd = 1,3
               do jj = 1,dj
                  do ii = 1,di
                     do g = 1,ng
                        dbuf2(g+ng*((jj-1)+dj*((ii-1)+di*(dd-1)))) = h34(g+ng*((ii-1)+di*((jj-1)+dj*(3+dd))))
                     enddo
                  enddo
               enddo
            enddo
            do dd = 1,9
               do jj = 1,dj
                  do ii = 1,di
                     do g = 1,ng
                        dbpp2(g+ng*((jj-1)+dj*((ii-1)+di*(dd-1)))) = h34(g+ng*((ii-1)+di*((jj-1)+dj*(24+dd))))
                     enddo
                  enddo
               enddo
            enddo
            if (do_d2 .and. d2gate .eq. 1) then
               m1 = 0.0d0
               do ob = 1, ng*di*dj*3
                  m1 = max(m1, abs(dbuf(ob)))
               enddo
               if (sj .ne. si) then
                  do ob = 1, ng*di*dj*3
                     m1 = max(m1, abs(dbuf2(ob)))
                  enddo
               endif
               sth = 0.0d0
               do is = 1,ns
                  do j = 1,dj
                     ao_j = cosx_ao_offset(sj) + j
                     do i = 1,di
                        ao_i = cosx_ao_offset(si) + i
                        nij = NorVEC(ao_i)*NorVEC(ao_j)
                        tt = 0.0d0
                        do g = 1,ng
                           tt = tt + abs(wv(g0+g)*(av(g0+g,ao_i,is)*bv(g0+g,ao_j,is) &
                                                 + av(g0+g,ao_j,is)*bv(g0+g,ao_i,is)))
                        enddo
                        sth = max(sth, nij*tt)
                     enddo
                  enddo
               enddo
               do_d2 = (sth*m1*dfac .ge. kscreen_h2)
            endif
            if (hchk) then
               call grids1e_engine_cached(cbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               hdiff = max(hdiff, maxval(abs(cbuf(1:nq)-buf(1:nq))))
               call grids1e_ip_engine_cached(cbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               hdiff = max(hdiff, maxval(abs(cbuf(1:3*nq)-dbuf(1:3*nq))))
               call grids1e_ipip_engine_cached(cbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               hdiff = max(hdiff, maxval(abs(cbuf(1:9*nq)-dbpp(1:9*nq))))
               call grids1e_ipvip_engine_cached(cbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               hdiff = max(hdiff, maxval(abs(cbuf(1:9*nq)-dbvp(1:9*nq))))
               if (sj .ne. si) then
                  shls1e(1) = sj ; shls1e(2) = si
                  call grids1e_ip_engine_cached(cbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
                  hdiff = max(hdiff, maxval(abs(cbuf(1:3*nq)-dbuf2(1:3*nq))))
                  call grids1e_ipip_engine_cached(cbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
                  hdiff = max(hdiff, maxval(abs(cbuf(1:9*nq)-dbpp2(1:9*nq))))
                  shls1e(1) = si ; shls1e(2) = sj
               endif
               hmag = max(hmag, maxval(abs(h34(1:34*nq))))
            endif
         endif
         if (.not. fused) then
         tei = omp_get_wtime() ; cn1 = cn1 + 1.0d0
         call grids1e_engine_cached(buf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            call grids1e_ip_engine_cached(dbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            if (sj .ne. si) then
               shls1e(1) = sj ; shls1e(2) = si
               call grids1e_ip_engine_cached(dbuf2, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               shls1e(1) = si ; shls1e(2) = sj
            endif
            te1 = te1 + (omp_get_wtime() - tei)
            if (do_d2 .and. d2gate .eq. 1 .and. nseg .eq. 0) then
               m1 = 0.0d0
               do ob = 1, ng*di*dj*3
                  m1 = max(m1, abs(dbuf(ob)))
               enddo
               if (sj .ne. si) then
                  do ob = 1, ng*di*dj*3
                     m1 = max(m1, abs(dbuf2(ob)))
                  enddo
               endif
               sth = 0.0d0
               do is = 1,ns
                  do j = 1,dj
                     ao_j = cosx_ao_offset(sj) + j
                     do i = 1,di
                        ao_i = cosx_ao_offset(si) + i
                        nij = NorVEC(ao_i)*NorVEC(ao_j)
                        tt = 0.0d0
                        do g = 1,ng
                           tt = tt + abs(wv(g0+g)*(av(g0+g,ao_i,is)*bv(g0+g,ao_j,is) &
                                                 + av(g0+g,ao_j,is)*bv(g0+g,ao_i,is)))
                        enddo
                        sth = max(sth, nij*tt)
                     enddo
                  enddo
               enddo
               do_d2 = (sth*m1*dfac .ge. kscreen_h2)
            endif
            tei = omp_get_wtime()
            if (do_d2 .and. nseg .eq. 0) then
               call grids1e_ipip_engine_cached(dbpp, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               call grids1e_ipvip_engine_cached(dbvp, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               if (sj .ne. si) then
                  shls1e(1) = sj ; shls1e(2) = si
                  call grids1e_ipip_engine_cached(dbpp2, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
                  shls1e(1) = si ; shls1e(2) = sj
               endif
            endif
            te2 = te2 + (omp_get_wtime() - tei)
            if (do_d2 .and. nseg .eq. 0) cn2 = cn2 + 1.0d0
         endif
         tint = tint + (omp_get_wtime() - tq2)
         if (scrd) then
            amx = 0.0d0 ; a1mx = 0.0d0
            do j = 1,dj
               ao_j = cosx_ao_offset(sj) + j
               do i = 1,di
                  ao_i = cosx_ao_offset(si) + i
                  nij = NorVEC(ao_i)*NorVEC(ao_j)
                  do g = 1,ng
                     amx(g) = max(amx(g), abs(nij*buf(g + ng*((i-1)+di*(j-1)))))
                     do dd = 1,3
                        a1mx(g) = max(a1mx(g), abs(nij*dbuf(g + ng*((i-1)+di*((j-1)+dj*(dd-1))))))
                        if (sj .ne. si) a1mx(g) = max(a1mx(g), abs(nij*dbuf2(g + ng*((j-1)+dj*((i-1)+di*(dd-1))))))
                     enddo
                  enddo
               enddo
            enddo
            bh = 0.0d0 ; bk = 0.0d0 ; bhs = 0.0d0
            do g = 1,ng
               bpt = amx(g)*Vp(g0+g,si)*Vp(g0+g,sj) + a1mx(g)*max(Vp(g0+g,si)*Fp(g0+g,sj), Fp(g0+g,si)*Vp(g0+g,sj))
               bh = max(bh, bpt)
               bhs = bhs + bpt
               if (want_kx) bk = max(bk, Lp(g0+g)*(amx(g)*max(Up(g0+g,si),Up(g0+g,sj)) &
                                      + a1mx(g)*max(Fup(g0+g,si),Fup(g0+g,sj))))
            enddo
            if (bh .ge. kscreen_h .or. bk .ge. kscreen_hk) scr(3) = scr(3) + wgt
            if (bhs .ge. kscreen_h .or. bk*ng .ge. kscreen_hk) scr(4) = scr(4) + wgt
            if (bh .lt. kscreen_h .and. bk .ge. kscreen_hk) scr(9) = scr(9) + wgt
            if (do_d2) scr(5) = scr(5) + wgt
            S2 = 0.0d0 ; P2 = 0.0d0
         endif

         do j = 1,dj
            ao_j = cosx_ao_offset(sj) + j
            do i = 1,di
               ao_i = cosx_ao_offset(si) + i
               nij = NorVEC(ao_i)*NorVEC(ao_j)
               if (nij .eq. 0.0d0) cycle
               do g = 1,ng
                  Aval(g) = nij*buf(g + ng*((i-1)+di*(j-1)))
               enddo
               if (acmp) then
                  ic = aoc(ao_i) ; jc = aoc(ao_j)
               else
                  ic = ao_i ; jc = ao_j
               endif
               do g = 1,ng
                  Asub(ic,jc,g) = Asub(ic,jc,g) + Aval(g)
               enddo
               if (sj .ne. si) then
                  do g = 1,ng
                     Asub(jc,ic,g) = Asub(jc,ic,g) + Aval(g)
                  enddo
               endif
               do dd = 1,3
                  ob = ng*((i-1) + di*((j-1) + dj*(dd-1)))
                  if (sj .ne. si) then
                     okk = ng*((j-1) + dj*((i-1) + di*(dd-1)))
                     do g = 1,ng
                        Abra(g,dd) = -nij*dbuf(ob+g)
                        Aket(g,dd) = -nij*dbuf2(okk+g)
                     enddo
                  else
                     okk = ng*((j-1) + di*((i-1) + dj*(dd-1)))
                     do g = 1,ng
                        Abra(g,dd) = -nij*dbuf(ob+g)
                        Aket(g,dd) = -nij*dbuf(okk+g)
                     enddo
                  endif
               enddo
               do is = 1,ns
                  do dd = 1,3
                     yi = dd + 3*(atom_i-1)
                     yj = dd + 3*(atom_j-1)
                     do g = 1,ng
                        kv(kg0+g,yi,ao_i,is) = kv(kg0+g,yi,ao_i,is) + Abra(g,dd)*bv(g0+g,ao_j,is)
                        kv(kg0+g,yj,ao_i,is) = kv(kg0+g,yj,ao_i,is) + Aket(g,dd)*bv(g0+g,ao_j,is)
                        kp(g0+g,ao_i,dd,is) = kp(g0+g,ao_i,dd,is) - (Abra(g,dd)+Aket(g,dd))*bv(g0+g,ao_j,is)
                        mv(kg0+g,yi,ao_i,is) = mv(kg0+g,yi,ao_i,is) + Abra(g,dd)*av(g0+g,ao_j,is)
                        mv(kg0+g,yj,ao_i,is) = mv(kg0+g,yj,ao_i,is) + Aket(g,dd)*av(g0+g,ao_j,is)
                        mp(g0+g,ao_i,dd,is) = mp(g0+g,ao_i,dd,is) - (Abra(g,dd)+Aket(g,dd))*av(g0+g,ao_j,is)
                     enddo
                  enddo
                  if (sj .ne. si) then
                     do dd = 1,3
                        yi = dd + 3*(atom_i-1)
                        yj = dd + 3*(atom_j-1)
                        do g = 1,ng
                           kv(kg0+g,yi,ao_j,is) = kv(kg0+g,yi,ao_j,is) + Abra(g,dd)*bv(g0+g,ao_i,is)
                           kv(kg0+g,yj,ao_j,is) = kv(kg0+g,yj,ao_j,is) + Aket(g,dd)*bv(g0+g,ao_i,is)
                           kp(g0+g,ao_j,dd,is) = kp(g0+g,ao_j,dd,is) - (Abra(g,dd)+Aket(g,dd))*bv(g0+g,ao_i,is)
                           mv(kg0+g,yi,ao_j,is) = mv(kg0+g,yi,ao_j,is) + Abra(g,dd)*av(g0+g,ao_i,is)
                           mv(kg0+g,yj,ao_j,is) = mv(kg0+g,yj,ao_j,is) + Aket(g,dd)*av(g0+g,ao_i,is)
                           mp(g0+g,ao_j,dd,is) = mp(g0+g,ao_j,dd,is) - (Abra(g,dd)+Aket(g,dd))*av(g0+g,ao_i,is)
                        enddo
                     enddo
                  endif
                  if (.not. do_d2 .or. nseg .gt. 0) cycle
                  if (sj .ne. si) then
                     do g = 1,ng
                        th(g) = wv(g0+g)*(av(g0+g,ao_i,is)*bv(g0+g,ao_j,is) + av(g0+g,ao_j,is)*bv(g0+g,ao_i,is))
                     enddo
                  else
                     do g = 1,ng
                        th(g) = wv(g0+g)*av(g0+g,ao_i,is)*bv(g0+g,ao_j,is)
                     enddo
                  endif
                  sbb = 0.0d0 ; sbk = 0.0d0 ; skk = 0.0d0
                  do dd = 1,3
                     do ee = 1,3
                        do g = 1,ng
                           ipp  = nij*dbpp(g + ng*((i-1)+di*((j-1)+dj*((dd-1)*3+(ee-1)))))
                           ivp  = nij*dbvp(g + ng*((i-1)+di*((j-1)+dj*((dd-1)*3+(ee-1)))))
                           ivpt = nij*dbvp(g + ng*((i-1)+di*((j-1)+dj*((ee-1)*3+(dd-1)))))
                           if (sj .ne. si) then
                              ikk = nij*dbpp2(g + ng*((j-1)+dj*((i-1)+di*((dd-1)*3+(ee-1)))))
                           else
                              ikk = nij*dbpp(g + ng*((j-1)+di*((i-1)+dj*((dd-1)*3+(ee-1)))))
                           endif
                           sbb(dd,ee) = sbb(dd,ee) + th(g)*ipp
                           sbk(dd,ee) = sbk(dd,ee) + th(g)*ivp
                           skk(dd,ee) = skk(dd,ee) + th(g)*ikk
                           if (scrd) P2 = max(P2, abs(th(g)*(ipp+ivp)), abs(th(g)*(ivpt+ikk)))
                           Hp(g0+g,dd+3*(atom_i-1),ee,is) = Hp(g0+g,dd+3*(atom_i-1),ee,is) - th(g)*(ipp+ivp)
                           Hp(g0+g,dd+3*(atom_j-1),ee,is) = Hp(g0+g,dd+3*(atom_j-1),ee,is) - th(g)*(ivpt+ikk)
                           Hpp(g0+g,dd,ee,is) = Hpp(g0+g,dd,ee,is) + th(g)*(ipp+ivp+ivpt+ikk)
                        enddo
                     enddo
                  enddo
                  if (scrd) then
                     S2(:,:,1) = S2(:,:,1) + sbb
                     S2(:,:,2) = S2(:,:,2) + sbk
                     S2(:,:,3) = S2(:,:,3) + skk
                  endif
                  do dd = 1,3
                     do ee = 1,3
                        yi = dd + 3*(atom_i-1)
                        Hloc(yi, ee+3*(atom_i-1), 5, is) = Hloc(yi, ee+3*(atom_i-1), 5, is) + sbb(dd,ee)
                        Hloc(yi, ee+3*(atom_j-1), 5, is) = Hloc(yi, ee+3*(atom_j-1), 5, is) + sbk(dd,ee)
                        Hloc(ee+3*(atom_j-1), yi, 5, is) = Hloc(ee+3*(atom_j-1), yi, 5, is) + sbk(dd,ee)
                        Hloc(dd+3*(atom_j-1), ee+3*(atom_j-1), 5, is) = &
                           Hloc(dd+3*(atom_j-1), ee+3*(atom_j-1), 5, is) + skk(dd,ee)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         if (nseg .gt. 0 .and. do_d2) then
            do s0 = 0, ng-1, nseg
               nsg = min(nseg, ng-s0)
               if (d2gate .eq. 1) then
                  m1 = 0.0d0
                  do ob = 0, di*dj*3-1
                     do g = 1,nsg
                        m1 = max(m1, abs(dbuf(s0+g+ng*ob)))
                     enddo
                  enddo
                  if (sj .ne. si) then
                     do ob = 0, di*dj*3-1
                        do g = 1,nsg
                           m1 = max(m1, abs(dbuf2(s0+g+ng*ob)))
                        enddo
                     enddo
                  endif
                  sth = 0.0d0
                  do is = 1,ns
                     do j = 1,dj
                        ao_j = cosx_ao_offset(sj) + j
                        do i = 1,di
                           ao_i = cosx_ao_offset(si) + i
                           nij = NorVEC(ao_i)*NorVEC(ao_j)
                           tt = 0.0d0
                           do g = g0+s0+1,g0+s0+nsg
                              tt = tt + abs(wv(g)*(av(g,ao_i,is)*bv(g,ao_j,is) + av(g,ao_j,is)*bv(g,ao_i,is)))
                           enddo
                           sth = max(sth, nij*tt)
                        enddo
                     enddo
                  enddo
                  if (sth*m1*dfac .lt. kscreen_h2) cycle
               endif
               tei = omp_get_wtime()
               shls1e(3) = g0+s0 ; shls1e(4) = g0+s0+nsg
               call grids1e_ipip_engine_cached(dbpp, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               call grids1e_ipvip_engine_cached(dbvp, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
               if (sj .ne. si) then
                  shls1e(1) = sj ; shls1e(2) = si
                  call grids1e_ipip_engine_cached(dbpp2, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
                  shls1e(1) = si ; shls1e(2) = sj
               endif
               shls1e(3) = g0 ; shls1e(4) = g0+ng
               tt = omp_get_wtime() - tei
               te2 = te2 + tt ; tint = tint + tt
               cn2 = cn2 + 1.0d0
               do is = 1,ns
                  do j = 1,dj
                     ao_j = cosx_ao_offset(sj) + j
                     atom_j = cosx_shell_atom(sj)
                     do i = 1,di
                        ao_i = cosx_ao_offset(si) + i
                        nij = NorVEC(ao_i)*NorVEC(ao_j)
                        if (nij .eq. 0.0d0) cycle
                        if (sj .ne. si) then
                           do g = 1,nsg
                              th(g) = wv(g0+s0+g)*(av(g0+s0+g,ao_i,is)*bv(g0+s0+g,ao_j,is) &
                                                  + av(g0+s0+g,ao_j,is)*bv(g0+s0+g,ao_i,is))
                           enddo
                        else
                           do g = 1,nsg
                              th(g) = wv(g0+s0+g)*av(g0+s0+g,ao_i,is)*bv(g0+s0+g,ao_j,is)
                           enddo
                        endif
                        sbb = 0.0d0 ; sbk = 0.0d0 ; skk = 0.0d0
                        do dd = 1,3
                           do ee = 1,3
                              do g = 1,nsg
                                 ipp  = nij*dbpp(g + nsg*((i-1)+di*((j-1)+dj*((dd-1)*3+(ee-1)))))
                                 ivp  = nij*dbvp(g + nsg*((i-1)+di*((j-1)+dj*((dd-1)*3+(ee-1)))))
                                 ivpt = nij*dbvp(g + nsg*((i-1)+di*((j-1)+dj*((ee-1)*3+(dd-1)))))
                                 if (sj .ne. si) then
                                    ikk = nij*dbpp2(g + nsg*((j-1)+dj*((i-1)+di*((dd-1)*3+(ee-1)))))
                                 else
                                    ikk = nij*dbpp(g + nsg*((j-1)+di*((i-1)+dj*((dd-1)*3+(ee-1)))))
                                 endif
                                 sbb(dd,ee) = sbb(dd,ee) + th(g)*ipp
                                 sbk(dd,ee) = sbk(dd,ee) + th(g)*ivp
                                 skk(dd,ee) = skk(dd,ee) + th(g)*ikk
                                 Hp(g0+s0+g,dd+3*(atom_i-1),ee,is) = Hp(g0+s0+g,dd+3*(atom_i-1),ee,is) - th(g)*(ipp+ivp)
                                 Hp(g0+s0+g,dd+3*(atom_j-1),ee,is) = Hp(g0+s0+g,dd+3*(atom_j-1),ee,is) - th(g)*(ivpt+ikk)
                                 Hpp(g0+s0+g,dd,ee,is) = Hpp(g0+s0+g,dd,ee,is) + th(g)*(ipp+ivp+ivpt+ikk)
                              enddo
                           enddo
                        enddo
                        do dd = 1,3
                           do ee = 1,3
                              yi = dd + 3*(atom_i-1)
                              Hloc(yi, ee+3*(atom_i-1), 5, is) = Hloc(yi, ee+3*(atom_i-1), 5, is) + sbb(dd,ee)
                              Hloc(yi, ee+3*(atom_j-1), 5, is) = Hloc(yi, ee+3*(atom_j-1), 5, is) + sbk(dd,ee)
                              Hloc(ee+3*(atom_j-1), yi, 5, is) = Hloc(ee+3*(atom_j-1), yi, 5, is) + sbk(dd,ee)
                              Hloc(dd+3*(atom_j-1), ee+3*(atom_j-1), 5, is) = &
                                 Hloc(dd+3*(atom_j-1), ee+3*(atom_j-1), 5, is) + skk(dd,ee)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         endif
         if (scrd .and. do_d2) then
            P2 = max(P2*ng, maxval(abs(S2)))
            if (P2 .ge. kscreen_h2) scr(6) = scr(6) + wgt
            if (P2 .ge. 1.0d1*kscreen_h2) scr(7) = scr(7) + wgt
            if (P2 .ge. 1.0d2*kscreen_h2) scr(8) = scr(8) + wgt
         endif
      enddo
   enddo
   if (hdiag) then
   block
      integer :: kq
      real(8) :: nzc
      nzc = 0.0d0
      do kq = 1,n
         if (any(Asub(:,kq,1) .ne. 0.0d0)) nzc = nzc + 1.0d0
      enddo
      afill = afill + dble(count(Asub(:,:,1) .ne. 0.0d0))/dble(n)**2
      arow = arow + nzc/dble(n)
      nsubb = nsubb + 1.0d0
   end block
   endif
   tq3 = omp_get_wtime()
   do is = 1,ns
      do gb = 1,ng,NGRP
         ngp = min(NGRP, ng-gb+1)
         nb0 = 0
         if (bfly) then
            ncd = 0
            do yy = 1,N3
               if (any(bnm(g0+gb:g0+gb+ngp-1,yy,is) .ge. colthr)) then
                  ncd = ncd + 1
                  cand(ncd) = yy
               endif
            enddo
            if (ncd .gt. 0) then
               do gp = 1,ngp
                  do c = 1,ncd
                     Bg(:,c) = B1C(g0+gb+gp-1,:,cand(c),is)
                  enddo
                  call dgemm('N','N',nsa,ncd,nop,1.0d0,Csl(1,1,is),n,Bg,nop,0.0d0,Rg(1,3,gp),n)
               enddo
               do c = 1,ncd
                  tt = 0.0d0
                  do gp = 1,ngp
                     tt = max(tt, sqrt(abs(wv(g0+gb+gp-1)))*maxval(abs(Rg(1:nsa,2+c,gp))))
                  enddo
                  if (tt .ge. colthr) then
                     nb0 = nb0 + 1
                     clist(nb0) = cand(c)
                     if (nb0 .ne. c) then
                        do gp = 1,ngp
                           Rg(1:nsa,2+nb0,gp) = Rg(1:nsa,2+c,gp)
                        enddo
                     endif
                  endif
               enddo
            endif
         else
         do yy = 1,N3
            if (any(csig(yy,g0+gb:g0+gb+ngp-1))) then
               nb0 = nb0 + 1
               clist(nb0) = yy
            endif
         enddo
         endif
         na0 = 0
         if (allocated(kav)) then
            do yy = 1,N3
               if (any(csig(N3+yy,g0+gb:g0+gb+ngp-1))) then
                  na0 = na0 + 1
                  clist(nb0+na0) = yy
               endif
            enddo
         endif
         do jc = 1,nsa
            j = jc
            if (acmp) j = slist(jc)
            do gp = 1,ngp
               Rg(jc,1,gp) = bv(g0+gb+gp-1,j,is)
               Rg(jc,2,gp) = av(g0+gb+gp-1,j,is)
            enddo
            if (.not. bfly) then
            do c = 1,nb0
               yy = clist(c)
               do gp = 1,ngp
                  Rg(jc,2+c,gp) = b1(merge(0,g0,bcomp)+gb+gp-1,yy,j,is)
               enddo
            enddo
            endif
            do c = nb0+1,nb0+na0
               yy = clist(c)
               do gp = 1,ngp
                  Rg(jc,2+c,gp) = a1(g0+gb+gp-1,yy,j,is)
               enddo
            enddo
         enddo
         do gp = 1,ngp
            call dgemm('N','N',nsa,2+nb0+na0,nsa,1.0d0,Abuf(1+size(Asub,1)**2*(gb+gp-2)),size(Asub,1),Rg(1,1,gp),n, &
                       0.0d0,Og(1,1,gp),n)
         enddo
         do ic = 1,nsa
            i = ic
            if (acmp) i = slist(ic)
            do gp = 1,ngp
               Gv(g0+gb+gp-1,i,is)  = Gv(g0+gb+gp-1,i,is)  + Og(ic,1,gp)
               Gtv(g0+gb+gp-1,i,is) = Gtv(g0+gb+gp-1,i,is) + Og(ic,2,gp)
            enddo
            do c = 1,nb0
               yy = clist(c)
               do gp = 1,ngp
                  kv(kg0+gb+gp-1,yy,i,is) = kv(kg0+gb+gp-1,yy,i,is) + Og(ic,2+c,gp)
               enddo
            enddo
            do c = nb0+1,nb0+na0
               yy = clist(c)
               do gp = 1,ngp
                  kav(g0+gb+gp-1,yy,i,is) = kav(g0+gb+gp-1,yy,i,is) + Og(ic,2+c,gp)
               enddo
            enddo
         enddo
         ncol_ar = ncol_ar + dble(ngp*(nb0+na0))
      enddo
   enddo
   tdns = tdns + (omp_get_wtime() - tq3)
   if (kcomp) then
      tq3 = omp_get_wtime()
      allocate(Tkc(nbsc*N3,nop))
      do is = 1,ns
         do ee = 1,3
            do g = g0+1,g0+ng
               pg = cosx_atom_hi(bstart+g-1)
               kv(g-g0,ee+3*(pg-1),:,is) = kv(g-g0,ee+3*(pg-1),:,is) + kp(g,:,ee,is)
               mv(g-g0,ee+3*(pg-1),:,is) = mv(g-g0,ee+3*(pg-1),:,is) + mp(g,:,ee,is)
            enddo
         enddo
         call dgemm('N','N',nbsc*N3,nop,n,1.0d0,kv(1,1,1,is),nbsc*N3,Cp(1,1,is),n,0.0d0,Tkc,nbsc*N3)
         do oo = 1,nop
            do yy = 1,N3
               KCc(merge(0,g0,psub)+1:merge(0,g0,psub)+ng,yy,oo,is) = Tkc(nbsc*(yy-1)+1:nbsc*(yy-1)+ng,oo)
            enddo
         enddo
         call dgemm('N','N',nbsc*N3,nop,n,1.0d0,mv(1,1,1,is),nbsc*N3,Cp(1,1,is),n,0.0d0,Tkc,nbsc*N3)
         do oo = 1,nop
            do yy = 1,N3
               MCc(merge(0,g0,psub)+1:merge(0,g0,psub)+ng,yy,oo,is) = Tkc(nbsc*(yy-1)+1:nbsc*(yy-1)+ng,oo)
            enddo
         enddo
      enddo
      deallocate(Tkc)
      tkx(3) = tkx(3) + (omp_get_wtime() - tq3)
   endif
   if (psub) then
      tq3 = omp_get_wtime()
      do is = 1,ns
         call hess_post_sub(ng, bstart+g0, ib, is, nloc, loc, xb(:,g0+1:g0+ng), x1b(:,:,g0+1:g0+ng), &
                            x2b(:,:,:,g0+1:g0+ng), wv(g0+1:g0+ng), av(g0+1:g0+ng,:,is), Gv(g0+1:g0+ng,:,is), &
                            Gtv(g0+1:g0+ng,:,is), qv(g0+1:g0+ng,:), KCc(1:ng,:,:,is), MCc(1:ng,:,:,is), &
                            XAv(g0+1:g0+ng,:,is), A1C(g0+1:g0+ng,:,:,is), B1C(g0+1:g0+ng,:,:,is), &
                            Hloc(:,:,:,is), gloc(:,is), Zloc(:,:,:,is), XWloc, Kloc, XKloc, tkx, ts)
      enddo
      tpost = tpost + (omp_get_wtime() - tq3)
   endif
   enddo

   do is = 1,ns
      call hess_scatter_par(nb, bstart, Hp(:,:,:,is), Hpp(:,:,:,is), Hloc(:,:,5,is))
   enddo

   tq2 = omp_get_wtime() ; tloop = tloop + (tq2 - tq1) ; tq1 = tq2
   if (.not. kcomp) then
   do is = 1,ns
      do ee = 1,3
         do g = 1,nb
            pg = cosx_atom_hi(bstart+g-1)
            kv(g,ee+3*(pg-1),:,is) = kv(g,ee+3*(pg-1),:,is) + kp(g,:,ee,is)
            mv(g,ee+3*(pg-1),:,is) = mv(g,ee+3*(pg-1),:,is) + mp(g,:,ee,is)
         enddo
      enddo
   enddo
   endif

   if (psub) then
      deallocate(qv)
   else
   allocate(hv(nb,n,ns), htv(nb,n,ns), Psi(nb,n,ns), epsy(nb,N3), wa(nb,N3,1), wh(nb,n))
   tsq = omp_get_wtime()
   if (want_kx .and. kproj) then
      allocate(QVc(nb,nop), GC(nb,nop), WC(nb,nop), wGC(nb,nop), Q1C(nb,N3,nop), KCa(nb,N3,nop))
      allocate(RR(nb,nop,N3,3), RS(N3*n,nop), wq(nb,n))
      if (fit) then
         allocate(XP(nb*nop*N3))
      else
         allocate(XP(1))
      endif
      do is = 1,ns
         do i = 1,n
            wq(:,i) = wv(:)*qv(:,i)
         enddo
         if (kcomp) then
            call kx_proj_leg(nb, bstart, is, xb, wv, wq, qv, x1b, Pt, PtCp(:,:,is), nloc, loc, &
                             Gv(:,:,is), kv(:,:,:,is), .true., &
                             QVc, GC, WC, wGC, Q1C, KCc(:,:,:,is), RR, RS, XP, Kloc, XWloc, tkx, kxx, XKloc, .true.)
         else
         call kx_proj_leg(nb, bstart, is, xb, wv, wq, qv, x1b, Pt, PtCp(:,:,is), nloc, loc, &
                          Gv(:,:,is), kv(:,:,:,is), .true., &
                          QVc, GC, WC, wGC, Q1C, KCa, RR, RS, XP, Kloc, XWloc, tkx, kxx, XKloc, .false.)
         if (t3lr) then
            if (.not. allocated(KCs)) allocate(KCs(nb,N3,nop,ns))
            KCs(:,:,:,is) = KCa
         endif
         endif
      enddo
      if (var_k) then
         do i = 1,n
            qv(:,i) = xb(i,:)
            wq(:,i) = wv(:)*xb(i,:)
         enddo
         do is = 1,ns
            kav(:,:,:,is) = kav(:,:,:,is) + mv(:,:,:,is)
            call kx_proj_leg(nb, bstart, is, xb, wv, wq, qv, x1b, Idn, Cp(:,:,is), nloc, loc, &
                             Gtv(:,:,is), kav(:,:,:,is), .false., &
                             QVc, GC, WC, wGC, Q1C, KCa, RR, RS, XP, Kloc, XWloc, tkx, .false., XKloc, .false.)
         enddo
      endif
      deallocate(QVc, GC, WC, wGC, Q1C, KCa, RR, RS, wq, qv)
      if (kcomp) call move_alloc(KCc, KCs)
      if (allocated(XP)) deallocate(XP)
   else if (want_kx) then
      allocate(Qx(nb,n), wq(nb,n))
      if (kproj) allocate(GC(nb,nop), WC(nb,nop), QC(nb,nop), KC(nb,nop))
      do i = 1,n
         wq(:,i) = wv(:)*qv(:,i)
      enddo
      do is = 1,ns
         if (kproj) then
            call dgemm('N','N',nb,nop,n,1.0d0,Gv(:,:,is),nb,Cp(1,1,is),n,0.0d0,GC,nb)
            call dgemm('N','N',nb,nop,n,1.0d0,wq,nb,Cp(1,1,is),n,0.0d0,WC,nb)
         endif
         do yy = 1,N3
            do i = 1,n
               Qx(:,i) = dwdr(mod(yy-1,3)+1,(yy-1)/3+1,bstart:bstart+nb-1)*qv(:,i) + wv(:)*q1(:,yy,i)
            enddo
            if (kproj) then
               call dgemm('N','N',nb,nop,n,1.0d0,Qx,nb,Cp(1,1,is),n,0.0d0,QC,nb)
               call dgemm('N','N',nb,nop,n,1.0d0,kv(1,yy,1,is),nb*N3,Cp(1,1,is),n,0.0d0,KC,nb)
               call dgemm('T','N',n,nop,nb,1.0d0,Qx,nb,GC,nb,1.0d0,Kloc(1,1,yy,is),n)
               call dgemm('T','N',n,nop,nb,1.0d0,Gv(:,:,is),nb,QC,nb,1.0d0,Kloc(1,1,yy,is),n)
               call dgemm('T','N',n,nop,nb,1.0d0,wq,nb,KC,nb,1.0d0,Kloc(1,1,yy,is),n)
               call dgemm('T','N',n,nop,nb,1.0d0,kv(1,yy,1,is),nb*N3,WC,nb,1.0d0,Kloc(1,1,yy,is),n)
            else
               call dgemm('T','N',n,n,nb,1.0d0,Qx,nb,Gv(:,:,is),nb,1.0d0,Kloc(1,1,yy,is),n)
               call dgemm('T','N',n,n,nb,1.0d0,wq,nb,kv(1,yy,1,is),nb*N3,1.0d0,Kloc(1,1,yy,is),n)
            endif
         enddo
      enddo
      if (var_k) then
         call move_alloc(x1t, q1)
         do i = 1,n
            qv(:,i) = xb(i,:)
            wq(:,i) = wv(:)*xb(i,:)
         enddo
         do is = 1,ns
            if (kproj) then
               call dgemm('N','N',nb,nop,n,1.0d0,Gtv(:,:,is),nb,Cp(1,1,is),n,0.0d0,GC,nb)
               call dgemm('N','N',nb,nop,n,1.0d0,wq,nb,Cp(1,1,is),n,0.0d0,WC,nb)
            endif
            do yy = 1,N3
               do i = 1,n
                  Qx(:,i) = dwdr(mod(yy-1,3)+1,(yy-1)/3+1,bstart:bstart+nb-1)*qv(:,i) + wv(:)*q1(:,yy,i)
               enddo
               kav(:,yy,:,is) = kav(:,yy,:,is) + mv(:,yy,:,is)
               if (kproj) then
                  call dgemm('N','N',nb,nop,n,1.0d0,Qx,nb,Cp(1,1,is),n,0.0d0,QC,nb)
                  call dgemm('N','N',nb,nop,n,1.0d0,kav(1,yy,1,is),nb*N3,Cp(1,1,is),n,0.0d0,KC,nb)
                  call dgemm('T','N',n,nop,nb,1.0d0,Qx,nb,GC,nb,1.0d0,Kloc(1,1,yy,is),n)
                  call dgemm('T','N',n,nop,nb,1.0d0,Gtv(:,:,is),nb,QC,nb,1.0d0,Kloc(1,1,yy,is),n)
                  call dgemm('T','N',n,nop,nb,1.0d0,wq,nb,KC,nb,1.0d0,Kloc(1,1,yy,is),n)
                  call dgemm('T','N',n,nop,nb,1.0d0,kav(1,yy,1,is),nb*N3,WC,nb,1.0d0,Kloc(1,1,yy,is),n)
               else
                  call dgemm('T','N',n,n,nb,1.0d0,Qx,nb,Gtv(:,:,is),nb,1.0d0,Kloc(1,1,yy,is),n)
                  call dgemm('T','N',n,n,nb,1.0d0,wq,nb,kav(1,yy,1,is),nb*N3,1.0d0,Kloc(1,1,yy,is),n)
               endif
            enddo
         enddo
      endif
      deallocate(qv, q1, Qx, wq)
      if (kproj) deallocate(GC, WC, QC, KC)
   endif
   ts(5) = ts(5) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
   do is = 1,ns
      call dgemm('N','N',nb,n,n,1.0d0,Gv(:,:,is),nb,Dsp(:,:,is),n,0.0d0,hv(:,:,is),nb)
      call dgemm('N','N',nb,n,n,1.0d0,Gtv(:,:,is),nb,Dsp(:,:,is),n,0.0d0,htv(:,:,is),nb)
      if (fit) then
         Psi(:,:,is) = htv(:,:,is)
         call dgemm('N','N',nb,n,n,1.0d0,hv(:,:,is),nb,Pf,n,1.0d0,Psi(:,:,is),nb)
      else
         Psi(:,:,is) = hv(:,:,is) + htv(:,:,is)
      endif
      do g = 1,nb
         epsall(bstart+g-1,is) = dot_product(av(g,:,is), Gv(g,:,is))
      enddo
      epsy = 0.0d0
      if (t3lr .and. allocated(KCs)) then
         allocate(GCt(nb,nop))
         call dgemm('N','N',nb,nop,n,1.0d0,Gv(:,:,is),nb,Cp(1,1,is),n,0.0d0,GCt,nb)
         do yy = 1,N3
            do oo = 1,nop
               do g = 1,nb
                  epsy(g,yy) = epsy(g,yy) + A1C(g,oo,yy,is)*GCt(g,oo)
               enddo
            enddo
         enddo
         deallocate(GCt)
         if (kcomp) then
            do oo = 1,nop
               do yy = 1,N3
                  do g = 1,nb
                     epsy(g,yy) = epsy(g,yy) + XAv(g,oo,is)*KCs(g,yy,oo,is)
                  enddo
               enddo
            enddo
         else
         do i = 1,n
            do yy = 1,N3
               do g = 1,nb
                  epsy(g,yy) = epsy(g,yy) + av(g,i,is)*kv(g,yy,i,is)
               enddo
            enddo
         enddo
         endif
      else
      do i = 1,n
         do yy = 1,N3
            do g = 1,nb
               epsy(g,yy) = epsy(g,yy) + a1(g,yy,i,is)*Gv(g,i,is) + av(g,i,is)*kv(g,yy,i,is)
            enddo
         enddo
      enddo
      endif
      do yy = 1,N3
         do g = 1,nb
            gloc(yy,is) = gloc(yy,is) + dwdr(mod(yy-1,3)+1,(yy-1)/3+1,bstart+g-1)*epsall(bstart+g-1,is) &
                          + wv(g)*epsy(g,yy)
         enddo
      enddo
      call dgemm('N','N',N3,N3,nb,1.0d0,dwdr(1,1,bstart),N3,epsy,nb,1.0d0,Hloc(1,1,2,is),N3)
      ts(6) = ts(6) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
      if (hdiag) then
      block
         integer :: yq, iq, t4
         real(8) :: rm(4)
         rblk = rblk + 1.0d0
         do yq = 1,N3
            rm = 0.0d0
            do iq = 1,n
               do g = 1,nb
                  rm(1) = max(rm(1), sqrt(abs(wv(g)))*abs(a1(g,yq,iq,is)))
                  rm(2) = max(rm(2), sqrt(abs(wv(g)))*abs(kv(g,yq,iq,is)))
                  rm(3) = max(rm(3), sqrt(abs(wv(g)))*abs(mv(g,yq,iq,is)))
                  rm(4) = max(rm(4), sqrt(abs(wv(g)))*abs(b1(g,yq,iq,is)))
               enddo
            enddo
            do t4 = 1,4
               where (rm .ge. 10.0d0**(-4-2*t4)) rsig(:,t4) = rsig(:,t4) + 1.0d0
            enddo
         enddo
      end block
      endif
      if (t3lr .and. allocated(KCs)) then
         allocate(KT(nb,nop,N3))
         do yy = 1,N3
            do oo = 1,nop
               do g = 1,nb
                  KT(g,oo,yy) = KCs(g,yy,oo,is)
                  A1C(g,oo,yy,is) = wv(g)*A1C(g,oo,yy,is)
               enddo
            enddo
         enddo
         call dgemm('T','N',N3,N3,nb*nop,1.0d0,A1C(1,1,1,is),nb*nop,KT,nb*nop,1.0d0,Hloc(1,1,3,is),N3)
         if (kcomp) then
            do yy = 1,N3
               do oo = 1,nop
                  do g = 1,nb
                     KT(g,oo,yy) = wv(g)*MCc(g,yy,oo,is)
                  enddo
               enddo
            enddo
         else
         allocate(GCt(nb*N3,nop))
         call dgemm('N','N',nb*N3,nop,n,1.0d0,mv(1,1,1,is),nb*N3,Cp(1,1,is),n,0.0d0,GCt,nb*N3)
         do yy = 1,N3
            do oo = 1,nop
               do g = 1,nb
                  KT(g,oo,yy) = wv(g)*GCt(g+nb*(yy-1),oo)
               enddo
            enddo
         enddo
         deallocate(GCt)
         endif
         if (bcomp .or. bfly) then
            call dgemm('T','N',N3,N3,nb*nop,1.0d0,KT,nb*nop,B1C(1,1,1,is),nb*nop,1.0d0,Hloc(1,1,3,is),N3)
         else
         A1C(:,:,:,is) = 0.0d0
         call hess_tangent_occ(nb, bstart, nloc, loc, x1b, Cp(:,:,is), A1C(:,:,:,is))
         call dgemm('T','N',N3,N3,nb*nop,1.0d0,KT,nb*nop,A1C(1,1,1,is),nb*nop,1.0d0,Hloc(1,1,3,is),N3)
         endif
         deallocate(KT)
      else
      do i = 1,n
         do yy = 1,N3
            wa(:,yy,1) = wv(:)*a1(:,yy,i,is)
         enddo
         call dgemm('T','N',N3,N3,nb,1.0d0,wa,nb,kv(1,1,i,is),nb,1.0d0,Hloc(1,1,3,is),N3)
         do yy = 1,N3
            wa(:,yy,1) = wv(:)*mv(:,yy,i,is)
         enddo
         call dgemm('T','N',N3,N3,nb,1.0d0,wa,nb,b1(1,1,i,is),nb,1.0d0,Hloc(1,1,3,is),N3)
      enddo
      endif
      ts(7) = ts(7) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
      Hp(:,:,:,1) = 0.0d0 ; Hpp(:,:,:,1) = 0.0d0
      call hess_xxy_contract(nb, bstart, ib, x2b, wv, Psi(:,:,is), Hloc(:,:,4,is), Hp(:,:,:,1), Hpp(:,:,:,1))
      call hess_scatter_par(nb, bstart, Hp(:,:,:,1), Hpp(:,:,:,1), Hloc(:,:,4,is))
      if (lrfit .or. lrz) then
         block
            real(8),allocatable :: wgc(:,:)
            allocate(wgc(nb,nop))
            call dgemm('N','N',nb,nop,n,1.0d0,Gv(:,:,is),nb,Cp(1,1,is),n,0.0d0,wgc,nb)
            do oo = 1,nop
               wgc(:,oo) = wv(:)*wgc(:,oo)
            enddo
            call hess_rows_shared_c(nb, bstart, is, x1b, wgc, Zloc(:,:,:,is))
            deallocate(wgc)
         end block
      else if (fit) then
         do g = 1,nb
            wh(g,:) = wv(g)*hv(g,:,is)
         enddo
         call dgemm('N','N',n,n,nb,1.0d0,xb,n,wh,nb,1.0d0,Vloc(1,1,is),n)
         if (zt7) then
            call hess_rows_t7(nb, bstart, nloc, loc, x1b, wh, Hloc(:,:,7,is))
         else if (zmode .eq. 1) then
            call hess_rows_shared(nb, bstart, is, x1b, wh, Zloc(:,:,:,is))
         else
            call hess_rows_par(nb, bstart, nloc, loc, x1b, wh, Zloc(:,:,:,is))
         endif
      endif
      ts(8) = ts(8) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
   enddo
   deallocate(hv, htv, Psi, epsy, wa, wh)
   endif
   deallocate(xb, x1b, x2b, wv, bv, av, Gv, Gtv, b1, a1, kv, mv, kp, mp)
   deallocate(Hp, Hpp, Aval, Abra, Aket, th)
   if (allocated(Vp)) deallocate(Vp, Fp, Up, Fup, Lp, amx, a1mx)
   if (allocated(kav)) deallocate(kav)
   if (allocated(A1C)) deallocate(A1C)
   if (allocated(KCs)) deallocate(KCs)
   if (allocated(KCc)) deallocate(KCc)
   if (allocated(B1C)) deallocate(B1C)
   if (allocated(MCc)) deallocate(MCc, XAv)
   deallocate(csig)
   tpost = tpost + (omp_get_wtime() - tq1)
enddo
!$omp end do
block
   integer :: tid
   do tid = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() .eq. tid) then
         Hterm = Hterm + Hloc
         if (kxx .and. .not. xkshr) XKsum = XKsum + XKloc
         if (want_kx .and. .not. kxonly) Kx = Kx + Kloc
         if (kproj .and. fit) XWg = XWg + XWloc
         grad1 = grad1 + gloc
         if (lrfit .or. lrz) then
            do is = 1,ns
               do ee = 1,3
                  do i = 1,natoms
                     Zt(at_off(i)+1:at_off(i)+at_n(i),:,ee+3*(i-1),is) = &
                        Zt(at_off(i)+1:at_off(i)+at_n(i),:,ee+3*(i-1),is) &
                        - Zloc(at_off(i)+1:at_off(i)+at_n(i),:,ee,is)
                  enddo
               enddo
            enddo
         else if (fit) then
            Vm = Vm + Vloc
            if (zmode .eq. 1) then
               do is = 1,ns
                  do ee = 1,3
                     do i = 1,natoms
                        Zy(at_off(i)+1:at_off(i)+at_n(i),:,ee+3*(i-1),is) = &
                           Zy(at_off(i)+1:at_off(i)+at_n(i),:,ee+3*(i-1),is) &
                           - Zloc(at_off(i)+1:at_off(i)+at_n(i),:,ee,is)
                     enddo
                  enddo
               enddo
            else if (.not. zt7) then
               Zy = Zy + Zloc
            endif
         endif
      endif
   enddo
end block
deallocate(XKloc)
deallocate(Hloc, gloc, buf, dbuf, dbuf2, dbpp, dbvp, dbpp2, int_cache, env_g, Fmax, Tmax, Fu, Tu, Abuf, Rg, Og, h34, cbuf)
nullify(Asub)
if (acmp) deallocate(pok, pd2, shon, aoc, slist)
if (fit) deallocate(Vloc, Zloc)
if (want_kx) deallocate(Kloc)
if (allocated(XWloc)) deallocate(XWloc)
!$omp end parallel
deallocate(sqa, ch_ib, ch_st, ch_nb)
if (kxx) then
   do is = 1,ns
      call dgemm('T','N',n,nop*N3,n,1.0d0,Pt,n,XKsum(1,1,1,is),n,1.0d0,Kx(1,1,1,is),n)
   enddo
   deallocate(XKsum)
endif
if (xkshr) then
   do i = 1,size(xklk)
      call omp_destroy_lock(xklk(i))
   enddo
   deallocate(xklk)
endif
if (allocated(zlock)) then
   do i = 1,natoms
      call omp_destroy_lock(zlock(i))
   enddo
   deallocate(zlock)
endif
if (verbose) then
   print '(A,5F10.2)', " COSX HESS main (thread-summed s): setup+tangents / ESP int / pair acc / dense A.R / post:", &
         tpre, tint, tloop-tint-tdns, tdns, tpost
   if (hchk) print '(A,ES10.3,A,ES10.3)', " COSX HESS34 CHECK max|fused - separate| =", hdiff, "  max|int| =", hmag
   print '(A,F7.2,A,I0,A,F7.3)', " COSX HESS block atoms: <nloc>_pts =", sloc/max(snb,1.0d0), " of ", natoms, &
         "   <parents outside loc>_pts =", spar/max(snb,1.0d0)
   print '(A,4F8.1,A,I0)', " COSX HESS <significant AOs>_pts at max|x,x',x''| >= 1e-6/-8/-10/-12:", &
         ssig/max(snb,1.0d0), "  of n=", n
   print '(A,2F10.2,A,2ES11.3)', " COSX HESS ESP (thread-s) val+ip+ip2 / ipip+ipvip+ipip2:", te1, te2, &
         "   calls:", cn1, cn2
   print '(A,F7.3,A,F7.3)', " COSX HESS A fill: <nonzero entries/n^2> =", afill/max(nsubb,1.0d0), &
         "   <nonzero columns/n> =", arow/max(nsubb,1.0d0)
   if (hdiag) print '(A,4F8.2,A,I0)', " COSX HESS <active atoms>_pts (unweighted max|x|,|x'| >= 1e-5/-7/-9/-11):", &
         sact/max(snb,1.0d0), "  of ", natoms
   if (hdiag) print '(A,16F7.1,A,I0)', " COSX HESS rows/block a1,kv,mv,b1 @1e-6 | 1e-8 | 1e-10 | 1e-12:", &
         rsig/max(rblk,1.0d0), " of ", N3
   print '(A,F8.2,A,I0)', " COSX HESS A.R tangent columns per point:", ncol_ar/max(nptall,1.0d0), " of ", 2*N3
   print '(A,6F9.2)', " COSX HESS Kx (thread-s): small dgemms / Q1C / KCa / RR build / 3 big dgemms / qT wGC+RS+scatter:", tkx
   print '(A,8F9.2)', " COSX HESS sub (thread-s): AOeval/tangent/fitDPy/Kxbound | Kx/T2eps/T3/T4+fitVZ:", ts(1:8)
   if (scrd) then
      print '(A,ES11.3)', " COSX SCRDIAG entered (sum ng*di*dj):", scr(1)
      print '(A,4F8.4)', " COSX SCRDIAG frac: kx-only-entered / pointwise-pre / measured-max / measured-sum:", &
            scr(11)/scr(1), scr(2)/scr(1), scr(3)/scr(1), scr(4)/scr(1)
      print '(A,2F8.4)', " COSX SCRDIAG frac: pointwise-pre Hess-part / measured kx-only:", scr(10)/scr(1), scr(9)/scr(1)
      print '(A,F8.3,A,F8.2,A,I0)', " COSX SCRDIAG <block radius>_pts (bohr):", scr(13)/scr(12), &
            "  <nloc>_pts:", scr(14)/scr(12), " of ", natoms
      print '(A,5F8.2,A,I0)', " COSX SCRDIAG <signif. tangent cols>_pts at sqrt|w|max|b^y,a^y| >= 1e-4..1e-8:", &
            scr(15:19)/scr(12), "  of 3N=", N3
      print '(A,4F8.4)', " COSX SCRDIAG d2: entered frac / measured >= h2, 10h2, 100h2 (of d2 entered):", &
            scr(5)/scr(1), scr(6)/max(scr(5),1.0d0), scr(7)/max(scr(5),1.0d0), scr(8)/max(scr(5),1.0d0)
   endif
   call flush(6)
endif

end subroutine hess_main_pass

end submodule exchange_cosx_hess_main_impl
