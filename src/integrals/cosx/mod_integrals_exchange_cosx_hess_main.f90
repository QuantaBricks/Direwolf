! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic Hessian - the main grid pass.

submodule (mod_integrals:exchange_cosx_hess_impl) exchange_cosx_hess_main_impl
implicit none
contains

module subroutine hess_main_pass()
real(8),allocatable :: xb(:,:), x1b(:,:,:), x2b(:,:,:,:), wv(:)
real(8),allocatable :: bv(:,:,:), av(:,:,:), Gv(:,:,:), Gtv(:,:,:), hv(:,:,:), htv(:,:,:), Psi(:,:,:)
real(8),allocatable :: kv(:,:,:,:), mv(:,:,:,:), kp(:,:,:,:), mp(:,:,:,:)
real(8),allocatable :: epsy(:,:), Hp(:,:,:,:), Hpp(:,:,:,:)
real(8),allocatable :: Hloc(:,:,:,:), gloc(:,:), Zloc(:,:,:,:)
real(8),allocatable :: env_g(:), buf(:), dbuf(:), dbuf2(:), dbpp(:), dbvp(:), dbpp2(:), int_cache(:)
real(8),allocatable :: Aval(:), Abra(:,:), Aket(:,:), th(:)
real(8) :: k_bound, esp_eff, r_far, seg_ab(3), seg_cl(3), seg_len2, seg_t, nij
real(8) :: sbb(3,3), sbk(3,3), skk(3,3), ipp, ivp, ivpt, ikk, tt
integer :: ib, nb, bstart, g, si, sj, sj_idx, di, dj, i, j, ao_i, ao_j, atom_i, atom_j
integer :: shls1e(4), off_grids, nloc, loc(natoms), yy, yi, yj, dd, ee, ob, okk, pg
real(8),allocatable :: Fmax(:), Tmax(:), sqa(:), Fu(:), Tu(:)
logical :: do_d2
real(8) :: dfac, Lw, kx_bound, m1, sth
character(len=32) :: envs2
real(8),allocatable :: Rg(:,:,:), Og(:,:,:)
integer,parameter :: NGRP = 8
integer :: mrhs, gb, gp, ngp, nb0, c, oo, ncd, cand(2*N3), clist(2*N3)
real(8),allocatable,target :: Abuf(:)
real(8),pointer :: Asub(:,:,:)
integer :: g0, ng, ip_
real(8),allocatable :: qv(:,:), XWloc(:,:,:)
real(8),allocatable :: A1C(:,:,:,:), B1C(:,:,:,:), KT(:,:,:), GCt(:,:), XPa(:,:,:)
real(8),allocatable :: KCc(:,:,:,:), MCc(:,:,:,:), XAv(:,:,:), Tkc(:,:)
real(8) :: Lpt(BATCH)
real(8) :: tpre, tint, tloop, tpost, tq1, tq2, tdns, tq3, ts(10), tsq, te1, te2, cn1, cn2, tei
real(8) :: tkx(6), ncol_ar, nptall
real(8),allocatable :: sw(:)
integer :: nbh, nch, ich, cb
integer :: nseg, s0, nsg
integer :: acap, nsa, nbsc, ic, jc
integer :: kvcap
real(8),allocatable :: Csl(:,:,:), Bg(:,:), bnm(:,:,:)
real(8) :: cms(2), colthr
logical :: psub
logical,allocatable :: pok(:), pd2(:), shon(:)
integer,allocatable :: aoc(:), slist(:)
integer,allocatable :: ch_ib(:), ch_st(:), ch_nb(:)

if (.not. (want_kx .and. kproj .and. lowrank)) then
   print *, "cosx hess_main_pass: needs dK/dR.C with the occupied coefficients (D = C C^T),"
   print *, "the dense / full-Kx paths are archived"
   print *, "in src/integrals/cosx/legacy (commit ffdd90d1)."
   stop 1
endif
tpre = 0.0d0 ; tint = 0.0d0 ; tloop = 0.0d0 ; tpost = 0.0d0 ; tdns = 0.0d0 ; ts = 0.0d0
tkx = 0.0d0 ; te1 = 0.0d0 ; te2 = 0.0d0 ; cn1 = 0.0d0 ; cn2 = 0.0d0
ncol_ar = 0.0d0 ; nptall = 0.0d0
nseg = max(1, min(BATCH, 4000000/(n*n)))
call get_environment_variable("ENGINE_COSX_HESS_D2SEG", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) nseg
nseg = max(1, nseg)
colthr = 1.0d-7
call get_environment_variable("ENGINE_COSX_HESS_COLTHR", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) colthr
allocate(cnrm(n,ns), ao_sh(n))
do si = 0,nBases-1
   ao_sh(cosx_ao_offset(si)+1:cosx_ao_offset(si)+cosx_shell_dim(si)) = si
enddo
do is = 1,ns
   do i = 1,n
      cnrm(i,is) = sqrt(sum(Cp(i,1:nop,is)**2))
   enddo
enddo
acap = 8000000
call get_environment_variable("ENGINE_COSX_HESS_ACAP", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) acap
kvcap = 8000000
call get_environment_variable("ENGINE_COSX_HESS_KVCAP", envs2)
if (len_trim(envs2) .gt. 0) read(envs2,*) kvcap
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
psub = .false.
call get_environment_variable("ENGINE_COSX_HESS_PSUB", envs2)
if (trim(envs2) .eq. "1") psub = fit
allocate(XKsum(n,nop,N3,ns))
!$omp parallel do private(yy) schedule(static)
do yy = 1,N3
   XKsum(:,:,yy,:) = 0_16
enddo
!$omp end parallel do
allocate(xklk((N3+7)/8))
do i = 1,(N3+7)/8
   call omp_init_lock(xklk(i))
enddo
if (fit) then
   allocate(zlock(natoms))
   do i = 1,natoms
      call omp_init_lock(zlock(i))
   enddo
endif
!$omp parallel default(shared) &
!$omp&  private(ich,xb,x1b,x2b,wv,bv,av,Gv,Gtv,hv,htv,Psi,kv,mv,kp,mp,epsy,Hp,Hpp) &
!$omp&  private(Hloc,gloc,Zloc,env_g,buf,dbuf,dbuf2,dbpp,dbvp,dbpp2,int_cache,Aval,Abra,Aket,th) &
!$omp&  private(k_bound,esp_eff,r_far,seg_ab,seg_cl,seg_len2,seg_t,nij,sbb,sbk,skk) &
!$omp&  private(ipp,ivp,ivpt,ikk,tt,ib,nb,bstart,g,si,sj,sj_idx,di,dj,i,j,ao_i,ao_j,atom_i,atom_j) &
!$omp&  private(shls1e,off_grids,nloc,loc,yy,yi,yj,dd,ee,ob,okk,pg,is,Fmax,Tmax,Fu,Tu,do_d2) &
!$omp&  private(qv,XWloc,A1C,B1C,KT,GCt,XPa) &
!$omp&  private(tq1,tq2,tq3,tsq,tei,dfac,Lw,kx_bound,m1,sth,Asub,Abuf,Rg,Og,g0,ng,mrhs,gb,gp,ngp) &
!$omp&  private(s0,nsg,nsa,nbsc,ic,jc,pok,pd2,shon,aoc,slist,KCc,MCc,XAv,Tkc) &
!$omp&  private(Csl,Bg,bnm,cms,ncd,cand,Lpt,clist,nb0,c,oo,sw) &
!$omp&  reduction(+:tpre,tint,tloop,tpost,tdns,ts,te1,te2,cn1,cn2,tkx,ncol_ar,nptall)
allocate(Hloc(N3,N3,NT,ns), gloc(N3,ns))
if (fit) then
   allocate(XWloc(n,nop,ns), Zloc(n,nop,3,ns))
   Zloc = 0.0d0
else
   allocate(XWloc(1,1,1))
endif
XWloc = 0.0d0
Hloc = 0.0d0
gloc = 0.0d0
allocate(buf(BATCH*maxd*maxd), dbuf(BATCH*maxd*maxd*3), dbuf2(BATCH*maxd*maxd*3))
allocate(dbpp(BATCH*maxd*maxd*9), dbvp(BATCH*maxd*maxd*9), dbpp2(BATCH*maxd*maxd*9))
allocate(int_cache(INT_CACHE_SIZE), Fmax(0:nBases-1), Tmax(0:nBases-1))
allocate(Fu(0:nBases-1), Tu(0:nBases-1))
mrhs = 2 + 2*N3
allocate(Abuf(max(acap, n*n*nseg)), pok(cosx_ext_ptr(nBases)), pd2(cosx_ext_ptr(nBases)), &
         shon(0:nBases-1), aoc(n), slist(n))
allocate(Rg(n,mrhs,NGRP), Og(n,mrhs,NGRP))
allocate(Csl(n,nop,ns), Bg(nop,N3), bnm(BATCH,N3,ns))
allocate(env_g(size(env) + 3*BATCH))
env_g(1:size(env)) = env(1:size(env))
off_grids = size(env)
env_g(9) = hs_rsomega

!$omp do schedule(static,1)
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
   wv = cosx_w_hi(bstart:bstart+nb-1)
   nptall = nptall + nb
   allocate(bv(nb,n,ns), av(nb,n,ns), Gv(nb,n,ns), Gtv(nb,n,ns))
   allocate(B1C(nb,nop,N3,ns), A1C(nb,nop,N3,ns), XAv(nb,nop,ns))
   if (.not. psub) allocate(KCc(nb,N3,nop,ns), MCc(nb,N3,nop,ns))
   allocate(kp(nb,n,3,ns), mp(nb,n,3,ns))
   do is = 1,ns
      allocate(XPa(nb,nop,1))
      call dgemm('T','N',nb,nop,n,1.0d0,xb,n,Cp(1,1,is),n,0.0d0,XPa,nb)
      call dgemm('N','T',nb,n,nop,1.0d0,XPa,nb,Cp(1,1,is),n,0.0d0,bv(:,:,is),nb)
      call dgemm('T','N',nb,nop,n,1.0d0,xb,n,PtCp(1,1,is),n,0.0d0,XPa,nb)
      call dgemm('N','T',nb,n,nop,1.0d0,XPa,nb,Cp(1,1,is),n,0.0d0,av(:,:,is),nb)
      XAv(:,:,is) = XPa(:,:,1)
      deallocate(XPa)
      B1C(:,:,:,is) = 0.0d0
      call hess_tangent_occ(nb, bstart, nloc, loc, x1b, Cp(:,:,is), B1C(:,:,:,is))
      ts(2) = ts(2) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
      if (fit) then
         call dgemm('T','N',nb,nop*N3,n,1.0d0,xb,n,PtC(1,1,1,is),n,0.0d0,A1C(1,1,1,is),nb)
      else
         A1C(:,:,:,is) = 0.0d0
      endif
      call hess_tangent_occ(nb, bstart, nloc, loc, x1b, PtCp(:,:,is), A1C(:,:,:,is))
      ts(3) = ts(3) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
   enddo
   Gv = 0.0d0 ; Gtv = 0.0d0 ; kp = 0.0d0 ; mp = 0.0d0
   allocate(Hp(nb,N3,3,ns), Hpp(nb,3,3,ns), Aval(nb), Abra(nb,3), Aket(nb,3), th(nb))
   Hp = 0.0d0 ; Hpp = 0.0d0

   allocate(sw(nb))
   sw = sqrt(abs(wv))
   do si = 0,nBases-1
      Fmax(si) = 0.0d0 ; Fu(si) = 0.0d0
      do is = 1,ns
         do i = cosx_ao_offset(si)+1, cosx_ao_offset(si)+cosx_shell_dim(si)
            do g = 1,nb
               tt = max(abs(bv(g,i,is)), abs(av(g,i,is)))
               Fu(si) = max(Fu(si), tt)
               Fmax(si) = max(Fmax(si), sw(g)*tt)
            enddo
         enddo
      enddo
      Tmax(si) = Fmax(si)
      Tu(si) = Fu(si)
   enddo
   do is = 1,ns
      call a1_shell_max(nb, is, sw, B1C(:,:,:,is), Tu, Tmax)
      call a1_shell_max(nb, is, sw, A1C(:,:,:,is), Tu, Tmax)
   enddo
   deallocate(sw)
   tsq = omp_get_wtime()
   allocate(qv(nb,n))
   call dgemm('T','N',nb,n,n,1.0d0,xb,n,Pt,n,0.0d0,qv,nb)
   call kx_free_leg_max(nb, bstart, xb, x1b, wv, qv, Lpt)
   Lw = maxval(Lpt)
   do g = 1,nb
      env_g(off_grids+3*(g-1)+1:off_grids+3*g) = cosx_coor_hi(:,bstart+g-1)
   enddo
   env_g(12) = dble(nb)
   env_g(13) = dble(off_grids)
   ts(4) = ts(4) + (omp_get_wtime() - tsq)
   tq2 = omp_get_wtime() ; tpre = tpre + (tq2 - tq1) ; tq1 = tq2
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
         kx_bound = esp_eff*Lw*max(Tu(sj)*cmx(si), Tu(si)*cmx(sj), &
                                   dfac*Fu(sj)*cmx(si), dfac*Fu(si)*cmx(sj))
         if (k_bound .lt. kscreen_h .and. kx_bound .lt. kscreen_hk) cycle
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
   nbsc = min(nbsc, max(nseg, ((kvcap/max(1,2*N3*n*ns))/nseg)*nseg))
   Asub(1:max(nsa,1),1:max(nsa,1),1:nbsc) => Abuf
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
   allocate(kv(nbsc,N3,n,ns), mv(nbsc,N3,n,ns))
   if (psub) allocate(KCc(nbsc,N3,nop,ns), MCc(nbsc,N3,nop,ns))

   do g0 = 0, nb-1, nbsc
   ng = min(nbsc, nb-g0)
   Asub(:,:,1:ng) = 0.0d0
   kv = 0.0d0 ; mv = 0.0d0
   do si = 0,nBases-1
      di = cosx_shell_dim(si)
      atom_i = cosx_shell_atom(si)
      do sj_idx = cosx_ext_ptr(si)+1, cosx_ext_ptr(si+1)
         sj = cosx_ext_idx(sj_idx)
         if (.not. pok(sj_idx)) cycle
         dfac = 2.0d0*(sqa(si)+sqa(sj))
         do_d2 = pd2(sj_idx)
         dj = cosx_shell_dim(sj)
         atom_j = cosx_shell_atom(sj)
         shls1e(1) = si ; shls1e(2) = sj ; shls1e(3) = g0 ; shls1e(4) = g0+ng
         tq2 = omp_get_wtime()
         tei = omp_get_wtime() ; cn1 = cn1 + 1.0d0
         call grids1e_engine_cached(buf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         call grids1e_ip_engine_cached(dbuf, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
         if (sj .ne. si) then
            shls1e(1) = sj ; shls1e(2) = si
            call grids1e_ip_engine_cached(dbuf2, shls1e, atm, size(atm,2), bas, nBases, env_g, 0_8, int_cache)
            shls1e(1) = si ; shls1e(2) = sj
         endif
         te1 = te1 + (omp_get_wtime() - tei)
         tint = tint + (omp_get_wtime() - tq2)

         do j = 1,dj
            ao_j = cosx_ao_offset(sj) + j
            do i = 1,di
               ao_i = cosx_ao_offset(si) + i
               nij = NorVEC(ao_i)*NorVEC(ao_j)
               if (nij .eq. 0.0d0) cycle
               do g = 1,ng
                  Aval(g) = nij*buf(g + ng*((i-1)+di*(j-1)))
               enddo
               ic = aoc(ao_i) ; jc = aoc(ao_j)
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
                        kv(g,yi,ao_i,is) = kv(g,yi,ao_i,is) + Abra(g,dd)*bv(g0+g,ao_j,is)
                        kv(g,yj,ao_i,is) = kv(g,yj,ao_i,is) + Aket(g,dd)*bv(g0+g,ao_j,is)
                        kp(g0+g,ao_i,dd,is) = kp(g0+g,ao_i,dd,is) - (Abra(g,dd)+Aket(g,dd))*bv(g0+g,ao_j,is)
                        mv(g,yi,ao_i,is) = mv(g,yi,ao_i,is) + Abra(g,dd)*av(g0+g,ao_j,is)
                        mv(g,yj,ao_i,is) = mv(g,yj,ao_i,is) + Aket(g,dd)*av(g0+g,ao_j,is)
                        mp(g0+g,ao_i,dd,is) = mp(g0+g,ao_i,dd,is) - (Abra(g,dd)+Aket(g,dd))*av(g0+g,ao_j,is)
                     enddo
                  enddo
                  if (sj .ne. si) then
                     do dd = 1,3
                        yi = dd + 3*(atom_i-1)
                        yj = dd + 3*(atom_j-1)
                        do g = 1,ng
                           kv(g,yi,ao_j,is) = kv(g,yi,ao_j,is) + Abra(g,dd)*bv(g0+g,ao_i,is)
                           kv(g,yj,ao_j,is) = kv(g,yj,ao_j,is) + Aket(g,dd)*bv(g0+g,ao_i,is)
                           kp(g0+g,ao_j,dd,is) = kp(g0+g,ao_j,dd,is) - (Abra(g,dd)+Aket(g,dd))*bv(g0+g,ao_i,is)
                           mv(g,yi,ao_j,is) = mv(g,yi,ao_j,is) + Abra(g,dd)*av(g0+g,ao_i,is)
                           mv(g,yj,ao_j,is) = mv(g,yj,ao_j,is) + Aket(g,dd)*av(g0+g,ao_i,is)
                           mp(g0+g,ao_j,dd,is) = mp(g0+g,ao_j,dd,is) - (Abra(g,dd)+Aket(g,dd))*av(g0+g,ao_i,is)
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo

         if (do_d2) then
            do s0 = 0, ng-1, nseg
               nsg = min(nseg, ng-s0)
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
      enddo
   enddo

   tq3 = omp_get_wtime()
   do is = 1,ns
      do gb = 1,ng,NGRP
         ngp = min(NGRP, ng-gb+1)
         nb0 = 0
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
         do jc = 1,nsa
            j = slist(jc)
            do gp = 1,ngp
               Rg(jc,1,gp) = bv(g0+gb+gp-1,j,is)
               Rg(jc,2,gp) = av(g0+gb+gp-1,j,is)
            enddo
         enddo
         do gp = 1,ngp
            call dgemm('N','N',nsa,2+nb0,nsa,1.0d0,Abuf(1+size(Asub,1)**2*(gb+gp-2)),size(Asub,1),Rg(1,1,gp),n, &
                       0.0d0,Og(1,1,gp),n)
         enddo
         do ic = 1,nsa
            i = slist(ic)
            do gp = 1,ngp
               Gv(g0+gb+gp-1,i,is)  = Gv(g0+gb+gp-1,i,is)  + Og(ic,1,gp)
               Gtv(g0+gb+gp-1,i,is) = Gtv(g0+gb+gp-1,i,is) + Og(ic,2,gp)
            enddo
            do c = 1,nb0
               yy = clist(c)
               do gp = 1,ngp
                  kv(gb+gp-1,yy,i,is) = kv(gb+gp-1,yy,i,is) + Og(ic,2+c,gp)
               enddo
            enddo
         enddo
         ncol_ar = ncol_ar + dble(ngp*nb0)
      enddo
   enddo
   tdns = tdns + (omp_get_wtime() - tq3)

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
   if (psub) then
      tq3 = omp_get_wtime()
      do is = 1,ns
         call hess_post_sub(ng, bstart+g0, ib, is, nloc, loc, xb(:,g0+1:g0+ng), x1b(:,:,g0+1:g0+ng), &
                            x2b(:,:,:,g0+1:g0+ng), wv(g0+1:g0+ng), av(g0+1:g0+ng,:,is), Gv(g0+1:g0+ng,:,is), &
                            Gtv(g0+1:g0+ng,:,is), KCc(1:ng,:,:,is), MCc(1:ng,:,:,is), &
                            XAv(g0+1:g0+ng,:,is), A1C(g0+1:g0+ng,:,:,is), B1C(g0+1:g0+ng,:,:,is), &
                            Hloc(:,:,:,is), gloc(:,is), Zloc(:,:,:,is), XWloc, tkx, ts)
      enddo
      tpost = tpost + (omp_get_wtime() - tq3)
   endif
   enddo

   do is = 1,ns
      call hess_scatter_par(nb, bstart, Hp(:,:,:,is), Hpp(:,:,:,is), Hloc(:,:,5,is))
   enddo
   tq2 = omp_get_wtime() ; tloop = tloop + (tq2 - tq1) ; tq1 = tq2

   if (psub) then
      deallocate(qv)
   else
   allocate(hv(nb,n,ns), htv(nb,n,ns), Psi(nb,n,ns), epsy(nb,N3))
   tsq = omp_get_wtime()
   do is = 1,ns
      call kx_proj_leg(nb, bstart, is, xb, wv, x1b, nloc, loc, Gv(:,:,is), KCc(:,:,:,is), XWloc, tkx)
   enddo
   deallocate(qv)
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
      do oo = 1,nop
         do yy = 1,N3
            do g = 1,nb
               epsy(g,yy) = epsy(g,yy) + XAv(g,oo,is)*KCc(g,yy,oo,is)
            enddo
         enddo
      enddo
      do yy = 1,N3
         do g = 1,nb
            gloc(yy,is) = gloc(yy,is) + dwdr(mod(yy-1,3)+1,(yy-1)/3+1,bstart+g-1)*epsall(bstart+g-1,is) &
                          + wv(g)*epsy(g,yy)
         enddo
      enddo
      call dgemm('N','N',N3,N3,nb,1.0d0,dwdr(1,1,bstart),N3,epsy,nb,1.0d0,Hloc(1,1,2,is),N3)
      ts(6) = ts(6) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
      allocate(KT(nb,nop,N3))
      do yy = 1,N3
         do oo = 1,nop
            do g = 1,nb
               KT(g,oo,yy) = KCc(g,yy,oo,is)
               A1C(g,oo,yy,is) = wv(g)*A1C(g,oo,yy,is)
            enddo
         enddo
      enddo
      call dgemm('T','N',N3,N3,nb*nop,1.0d0,A1C(1,1,1,is),nb*nop,KT,nb*nop,1.0d0,Hloc(1,1,3,is),N3)
      do yy = 1,N3
         do oo = 1,nop
            do g = 1,nb
               KT(g,oo,yy) = wv(g)*MCc(g,yy,oo,is)
            enddo
         enddo
      enddo
      call dgemm('T','N',N3,N3,nb*nop,1.0d0,KT,nb*nop,B1C(1,1,1,is),nb*nop,1.0d0,Hloc(1,1,3,is),N3)
      deallocate(KT)
      ts(7) = ts(7) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
      Hp(:,:,:,1) = 0.0d0 ; Hpp(:,:,:,1) = 0.0d0
      call hess_xxy_contract(nb, bstart, ib, x2b, wv, Psi(:,:,is), Hloc(:,:,4,is), Hp(:,:,:,1), Hpp(:,:,:,1))
      call hess_scatter_par(nb, bstart, Hp(:,:,:,1), Hpp(:,:,:,1), Hloc(:,:,4,is))
      if (fit) then
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
      endif
      ts(8) = ts(8) + (omp_get_wtime() - tsq) ; tsq = omp_get_wtime()
   enddo
   deallocate(hv, htv, Psi, epsy)
   endif
   deallocate(xb, x1b, x2b, wv, bv, av, Gv, Gtv, kv, mv, kp, mp)
   deallocate(Hp, Hpp, Aval, Abra, Aket, th)
   deallocate(A1C, B1C, KCc, MCc, XAv)
   tpost = tpost + (omp_get_wtime() - tq1)
enddo
!$omp end do
block
   integer :: tid
   do tid = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() .eq. tid) then
         Hterm = Hterm + Hloc
         if (fit) then
            XWg = XWg + XWloc
            do is = 1,ns
               do ee = 1,3
                  do i = 1,natoms
                     Zt(at_off(i)+1:at_off(i)+at_n(i),:,ee+3*(i-1),is) = &
                        Zt(at_off(i)+1:at_off(i)+at_n(i),:,ee+3*(i-1),is) &
                        - Zloc(at_off(i)+1:at_off(i)+at_n(i),:,ee,is)
                  enddo
               enddo
            enddo
         endif
         grad1 = grad1 + gloc
      endif
   enddo
end block
deallocate(Hloc, gloc, buf, dbuf, dbuf2, dbpp, dbvp, dbpp2, int_cache, env_g, Fmax, Tmax, Fu, Tu, Abuf, Rg, Og)
nullify(Asub)
deallocate(pok, pd2, shon, aoc, slist, Csl, Bg, bnm)
if (allocated(Zloc)) deallocate(Zloc)
deallocate(XWloc)
!$omp end parallel
deallocate(sqa, ch_ib, ch_st, ch_nb)
block
   real(8),allocatable :: XKr(:,:,:)
   integer :: y0, ny
   allocate(XKr(n,nop,8))
   do is = 1,ns
      do y0 = 1,N3,8
         ny = min(8, N3-y0+1)
         XKr(:,:,1:ny) = from_q(XKsum(:,:,y0:y0+ny-1,is))
         call dgemm('T','N',n,nop*ny,n,1.0d0,Pt,n,XKr,n,1.0d0,Kx(1,1,y0,is),n)
      enddo
   enddo
   deallocate(XKr)
end block
deallocate(XKsum)
do i = 1,size(xklk)
   call omp_destroy_lock(xklk(i))
enddo
deallocate(xklk)
if (allocated(zlock)) then
   do i = 1,natoms
      call omp_destroy_lock(zlock(i))
   enddo
   deallocate(zlock)
endif
if (verbose) then
   print '(A,5F10.2)', " COSX HESS main (thread-summed s): setup+tangents / ESP int / pair acc / dense A.R / post:", &
         tpre, tint, tloop-tint-tdns, tdns, tpost
   print '(A,2F10.2,A,2ES11.3)', " COSX HESS ESP (thread-s) val+ip+ip2 / ipip+ipvip+ipip2:", te1, te2, &
         "   calls:", cn1, cn2
   print '(A,F8.2,A,I0)', " COSX HESS A.R tangent columns per point:", ncol_ar/max(nptall,1.0d0), " of ", N3
   print '(A,6F9.2)', " COSX HESS Kx (thread-s): small dgemms / Q1C / KCa / RR build / 3 big dgemms / qT wGC+RS+scatter:", tkx
   print '(A,8F9.2)', " COSX HESS sub (thread-s): AOeval/tangent/fitA1C/Kxbound | Kx/T2eps/T3/T4+fitZ:", ts(1:8)
   call flush(6)
endif

end subroutine hess_main_pass

end submodule exchange_cosx_hess_main_impl
