! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic Hessian - S_num^y and P^y (full, by atom range, S^y.M).

submodule (mod_integrals:exchange_cosx_hess_impl) exchange_cosx_hess_setup_impl
implicit none
contains

module subroutine hess_sy_times(M, SM)
real(8),intent(in) :: M(n,nop)
real(8),intent(inout) :: SM(n,nop,N3)
integer,parameter :: AG = 16
real(8),allocatable :: xb(:,:), x1b(:,:,:), XM(:,:), WXM(:,:), X1M(:,:,:), Q(:,:,:), L(:,:,:), xc(:,:), x1c(:,:,:), &
                       xm1(:,:), rg(:,:), xg(:,:), xw1(:,:,:)
integer(omp_lock_kind) :: slock(natoms)
integer,allocatable :: aol(:), aoi(:)
integer :: ib, nb, bstart, nloc, loc(natoms), nl, kk, si, a, g, ee, A1, a_lo, a_hi, ny, yl, il, pp, np_, pg, i2, ia
integer :: idx(BATCH)
logical :: inl(natoms)
do A1 = 1,natoms
   call omp_init_lock(slock(A1))
enddo
!$omp parallel default(shared) private(xb,x1b,XM,WXM,X1M,Q,L,xc,x1c,xm1,rg,xg,xw1,aol,aoi,ib,nb,bstart,nloc,loc,nl,kk,si,a) &
!$omp&  private(g,ee,A1,a_lo,a_hi,ny,yl,il,pp,np_,pg,i2,idx,inl,ia)
allocate(aol(n), aoi(n))
!$omp do schedule(dynamic)
do ib = 0,nblocks-1
   bstart = cosx_blkstart_hi(ib) + 1
   nb = cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib)
   allocate(xb(n,nb), x1b(n,3,nb))
   call hess_block_ao(ib, nb, bstart, xb, x1b)
   call hess_block_atoms(ib, nloc, loc)
   inl = .false.
   inl(loc(1:nloc)) = .true.
   nl = 0
   aoi = 0
   do kk = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(kk)
      do a = 1,cosx_shell_dim(si)
         nl = nl + 1
         aol(nl) = cosx_ao_offset(si) + a
         aoi(aol(nl)) = nl
      enddo
   enddo
   allocate(XM(nb,nop), WXM(nb,nop), X1M(nb,nop,3), xc(nl,nb), x1c(nl,nb,3))
   call dgemm('T','N',nb,nop,n,1.0d0,xb,n,M,n,0.0d0,XM,nb)
   do ee = 1,3
      call dgemm('T','N',nb,nop,n,1.0d0,x1b(1,ee,1),3*n,M,n,0.0d0,X1M(1,1,ee),nb)
   enddo
   do g = 1,nb
      WXM(g,:) = cosx_w_hi(bstart+g-1)*XM(g,:)
      xc(:,g) = xb(aol(1:nl),g)
      do ee = 1,3
         x1c(:,g,ee) = x1b(aol(1:nl),ee,g)
      enddo
   enddo
   allocate(xw1(nl,nop,3))
   do ee = 1,3
      call dgemm('N','N',nl,nop,nb,1.0d0,x1c(1,1,ee),nl,WXM,nb,0.0d0,xw1(1,1,ee),nl)
   enddo
   do a_lo = 1,natoms,AG
      a_hi = min(natoms, a_lo+AG-1)
      ny = 3*(a_hi-a_lo+1)
      allocate(Q(nb,nop,ny), L(nl,nop,ny))
      do A1 = a_lo,a_hi
         do ee = 1,3
            yl = ee + 3*(A1-a_lo)
            do g = 1,nb
               Q(g,:,yl) = dwdr(ee,A1,bstart+g-1)*XM(g,:)
               if (cosx_atom_hi(bstart+g-1) .eq. A1) Q(g,:,yl) = Q(g,:,yl) + cosx_w_hi(bstart+g-1)*X1M(g,:,ee)
            enddo
            if (inl(A1)) then
               allocate(xm1(nb,nop))
               call dgemm('T','N',nb,nop,at_n(A1),1.0d0,x1b(at_off(A1)+1,ee,1),3*n,M(at_off(A1)+1,1),n, &
                          0.0d0,xm1,nb)
               do g = 1,nb
                  Q(g,:,yl) = Q(g,:,yl) - cosx_w_hi(bstart+g-1)*xm1(g,:)
               enddo
               deallocate(xm1)
            endif
         enddo
      enddo
      call dgemm('N','N',nl,nop*ny,nb,1.0d0,xc,nl,Q,nb,0.0d0,L,nl)
      do ee = 1,3
         do A1 = a_lo,a_hi
            if (.not. inl(A1)) cycle
            yl = ee + 3*(A1-a_lo)
            do i2 = at_off(A1)+1, at_off(A1)+at_n(A1)
               if (aoi(i2) .gt. 0) L(aoi(i2),:,yl) = L(aoi(i2),:,yl) - xw1(aoi(i2),:,ee)
            enddo
         enddo
      enddo
      allocate(rg(nb,nop), xm1(nl,nop), xg(nl,nb))
      do pp = a_lo,a_hi
         np_ = 0
         do g = 1,nb
            if (cosx_atom_hi(bstart+g-1) .eq. pp) then
               np_ = np_ + 1
               idx(np_) = g
            endif
         enddo
         if (np_ .eq. 0) cycle
         do pg = 1,np_
            rg(pg,:) = WXM(idx(pg),:)
         enddo
         do ee = 1,3
            do pg = 1,np_
               xg(:,pg) = x1c(:,idx(pg),ee)
            enddo
            call dgemm('N','N',nl,nop,np_,1.0d0,xg,nl,rg,nb,0.0d0,xm1,nl)
            yl = ee + 3*(pp-a_lo)
            L(:,:,yl) = L(:,:,yl) + xm1
         enddo
      enddo
      deallocate(rg, xm1, xg)
      do ia = 0, a_hi-a_lo
         A1 = a_lo + mod(ia + ib, a_hi-a_lo+1)
         call omp_set_lock(slock(A1))
         do yl = 3*(A1-a_lo)+1, 3*(A1-a_lo)+3
            do a = 1,nop
               do i2 = 1,nl
                  SM(aol(i2),a,3*(a_lo-1)+yl) = SM(aol(i2),a,3*(a_lo-1)+yl) + L(i2,a,yl)
               enddo
            enddo
         enddo
         call omp_unset_lock(slock(A1))
      enddo
      deallocate(Q, L)
   enddo
   deallocate(xb, x1b, XM, WXM, X1M, xc, x1c, xw1)
enddo
!$omp end do
deallocate(aol, aoi)
!$omp end parallel
do A1 = 1,natoms
   call omp_destroy_lock(slock(A1))
enddo
end subroutine hess_sy_times

module subroutine hess_snum_first(Sout)
use omp_lib, only: omp_get_max_threads
real(8),intent(out) :: Sout(n,n,N3)
real(8),allocatable :: Sall(:,:), xb(:,:), x1b(:,:,:), wx(:,:)
real(8),allocatable :: xc(:,:), xp(:,:), xm(:,:), cc(:,:), xm2(:), xg(:,:), rg(:,:)
integer,allocatable :: aol(:)
integer :: ib, nb, bstart, g, nloc, loc(natoms), yy, nl, kk, si, a, b, dd, A1, np_, nm_, tid, nth, t
integer(8) :: c, nel, nuse
real(8) :: dwv, st(6), sq0, tw0, tw(4), bmx, bmn, snmem
integer,allocatable :: bord(:)
real(8),allocatable :: bcst(:)
integer :: ibo, nab, a_lo, a_hi, nyr, y0, il, ee, pp, pg, idx(BATCH)
logical :: done(natoms)
character(len=32) :: envs3
nth = omp_get_max_threads()
st = 0.0d0 ; tw0 = omp_get_wtime() ; bmx = 0.0d0 ; bmn = 1.0d30 ; tw = 0.0d0
snmem = 1.0d9
call get_environment_variable("ENGINE_COSX_HESS_SNMEM", envs3)
if (len_trim(envs3) .gt. 0) read(envs3,*) snmem
nab = int(snmem/(24.0d0*dble(n)**2*dble(nth)))
nab = max(1, min(natoms, nab))
nab = (natoms + (natoms+nab-1)/nab - 1)/((natoms+nab-1)/nab)
nel = int(n,8)*int(n,8)*int(3*nab,8)
allocate(Sall(nel,0:nth-1))
allocate(bord(nblocks), bcst(nblocks))
do ib = 0,nblocks-1
   nl = 0
   do kk = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      nl = nl + cosx_shell_dim(cosx_blkidx_hi(kk))
   enddo
   bcst(ib+1) = dble(cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib))*dble(nl)**2
   bord(ib+1) = ib
enddo
do ib = 2,nblocks
   kk = bord(ib) ; dwv = bcst(ib) ; a = ib - 1
   do while (a .ge. 1)
      if (bcst(a) .ge. dwv) exit
      bord(a+1) = bord(a) ; bcst(a+1) = bcst(a)
      a = a - 1
   enddo
   bord(a+1) = kk ; bcst(a+1) = dwv
enddo
do a_lo = 1,natoms,nab
a_hi = min(natoms, a_lo+nab-1)
nyr = 3*(a_hi-a_lo+1)
y0 = 3*(a_lo-1)
tw0 = omp_get_wtime()
!$omp parallel default(shared) private(ib,ibo,nb,bstart,g,nloc,loc,yy,xb,x1b,wx,xc,xp,xm,cc,xm2,aol,nl,kk,si,a,b,dd,A1) &
!$omp&  private(np_,nm_,tid,dwv,sq0,xg,rg,il,ee,pp,pg,idx,done) reduction(+:st) reduction(max:bmx) reduction(min:bmn)
tid = omp_get_thread_num()
Sall(1:int(n,8)*int(n,8)*int(nyr,8),tid) = 0.0d0
!$omp barrier
!$omp master
tw(1) = tw(1) + (omp_get_wtime() - tw0) ; tw0 = omp_get_wtime()
!$omp end master
sq0 = omp_get_wtime()
allocate(aol(n))
!$omp do schedule(dynamic)
do ibo = 1,nblocks
   ib = bord(ibo)
   bstart = cosx_blkstart_hi(ib) + 1
   nb = cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib)
   sq0 = omp_get_wtime()
   allocate(xb(n,nb), x1b(n,3,nb), wx(nb,n))
   call hess_block_ao(ib, nb, bstart, xb, x1b)
   call hess_block_atoms(ib, nloc, loc)
   st(1) = st(1) + (omp_get_wtime()-sq0) ; sq0 = omp_get_wtime()
   do g = 1,nb
      wx(g,:) = cosx_w_hi(bstart+g-1)*xb(:,g)
   enddo
   nl = 0
   do kk = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(kk)
      do a = 1,cosx_shell_dim(si)
         nl = nl + 1
         aol(nl) = cosx_ao_offset(si) + a
      enddo
   enddo
   allocate(xc(nl,nb), xp(nl,nb), xm(nl,nb), cc(nl,nl), xm2(nb))
   do g = 1,nb
      xc(:,g) = xb(aol(1:nl),g)
      xm2(g) = maxval(abs(xc(:,g)))**2
   enddo
   do A1 = a_lo,a_hi
      do dd = 1,3
         yy = dd + 3*(A1-a_lo)
         np_ = 0 ; nm_ = 0
         do g = 1,nb
            dwv = dwdr(dd,A1,bstart+g-1)
            st(5) = st(5) + 1.0d0
            if (abs(dwv)*xm2(g) .lt. 1.0d-14) cycle
            st(6) = st(6) + 1.0d0
            if (dwv .gt. 0.0d0) then
               np_ = np_ + 1
               xp(:,np_) = xc(:,g)*sqrt(dwv)
            else if (dwv .lt. 0.0d0) then
               nm_ = nm_ + 1
               xm(:,nm_) = xc(:,g)*sqrt(-dwv)
            endif
         enddo
         if (np_ + nm_ .eq. 0) cycle
         cc = 0.0d0
         if (np_ .gt. 0) call dsyrk('U','N',nl,np_,1.0d0,xp,nl,0.0d0,cc,nl)
         if (nm_ .gt. 0) call dsyrk('U','N',nl,nm_,-1.0d0,xm,nl,1.0d0,cc,nl)
         do b = 1,nl
            do a = 1,b
               c = int(aol(a),8) + int(n,8)*(int(aol(b)-1,8) + int(n,8)*int(yy-1,8))
               Sall(c,tid) = Sall(c,tid) + merge(0.5d0*cc(a,b), cc(a,b), a .eq. b)
            enddo
         enddo
      enddo
   enddo
   deallocate(xc, xp, xm, cc, xm2)
   st(2) = st(2) + (omp_get_wtime()-sq0) ; sq0 = omp_get_wtime()
   do il = 1,nloc
      A1 = loc(il)
      if (A1 .lt. a_lo .or. A1 .gt. a_hi) cycle
      do ee = 1,3
         c = int(at_off(A1)+1,8) + int(n,8)*int(n,8)*int(ee+3*(A1-a_lo)-1,8)
         call dgemm('N','N',at_n(A1),n,nb,-1.0d0,x1b(at_off(A1)+1,ee,1),3*n,wx,nb,1.0d0,Sall(c,tid),n)
      enddo
   enddo
   allocate(xg(n,nb), rg(nb,n))
   done = .false.
   do g = 1,nb
      pp = cosx_atom_hi(bstart+g-1)
      if (done(pp)) cycle
      done(pp) = .true.
      if (pp .lt. a_lo .or. pp .gt. a_hi) cycle
      np_ = 0
      do pg = g,nb
         if (cosx_atom_hi(bstart+pg-1) .eq. pp) then
            np_ = np_ + 1
            idx(np_) = pg
         endif
      enddo
      do pg = 1,np_
         rg(pg,:) = wx(idx(pg),:)
      enddo
      do ee = 1,3
         do pg = 1,np_
            xg(:,pg) = x1b(:,ee,idx(pg))
         enddo
         c = 1_8 + int(n,8)*int(n,8)*int(ee+3*(pp-a_lo)-1,8)
         call dgemm('N','N',n,n,np_,1.0d0,xg,n,rg,nb,1.0d0,Sall(c,tid),n)
      enddo
   enddo
   deallocate(xg, rg)
   st(3) = st(3) + (omp_get_wtime()-sq0)
   deallocate(xb, x1b, wx)
enddo
!$omp end do nowait
bmx = max(bmx, omp_get_wtime() - tw0) ; bmn = min(bmn, omp_get_wtime() - tw0)
deallocate(aol)
!$omp end parallel
tw(2) = tw(2) + (omp_get_wtime() - tw0) ; tw0 = omp_get_wtime()
nuse = int(n,8)*int(n,8)*int(nyr,8)
!$omp parallel do private(c,t) schedule(static)
do c = 1,nuse
   Sall(c,0) = Sall(c,0) + sum(Sall(c,1:nth-1))
enddo
!$omp end parallel do
!$omp parallel do private(yy,kk) schedule(static)
do yy = 1,nyr
   do kk = 1,n
      Sout(:,kk,y0+yy) = Sall(int(n,8)*(int(kk-1,8)+int(n,8)*int(yy-1,8))+1: &
                              int(n,8)*(int(kk,8)+int(n,8)*int(yy-1,8)), 0)
   enddo
   Sout(:,:,y0+yy) = Sout(:,:,y0+yy) + transpose(Sout(:,:,y0+yy))
enddo
!$omp end parallel do
tw(3) = tw(3) + (omp_get_wtime() - tw0)
enddo
deallocate(Sall, bord, bcst)
if (verbose) print '(A,3F9.2,A,F7.4)', " COSX HESS S^y (thread-s): AO / w^y dsyrk+scatter / x^y rows:", st(1:3), &
      "   frac (g,y) |dw| max|x|^2 >= 1e-14:", st(6)/max(st(5),1.0d0)
if (verbose) print '(A,3F8.3,A,2F8.3,A,I0,A)', " COSX HESS S^y wall: zero / block loop / merge+sym:", tw(1:3), &
      "   thread busy max/min", bmx, bmn, "   (", nab, " atoms per range)"
end subroutine hess_snum_first

module subroutine hess_fit_py(Sin, Pout)
real(8),intent(in) :: Sin(n,n,N3)
real(8),intent(out) :: Pout(n,n,N3)
real(8),allocatable :: dS1e(:,:,:), dH1e(:,:,:), Tall(:,:,:)
integer :: yy, A, dd, k, l
allocate(dS1e(n,n,3), dH1e(n,n,3), Tall(n,n,N3))
call integrals_compute_force_1e(n, dS1e, dH1e)
!$omp parallel do private(yy,A,dd,k,l)
do yy = 1,N3
   A = (yy-1)/3 + 1
   dd = yy - 3*(A-1)
   Tall(:,:,yy) = 0.0d0
   do l = 1,n
      do k = 1,n
         if (ao_atom(k) .eq. A) Tall(k,l,yy) = Tall(k,l,yy) - dS1e(k,l,dd)
         if (ao_atom(l) .eq. A) Tall(k,l,yy) = Tall(k,l,yy) - dS1e(l,k,dd)
      enddo
   enddo
enddo
!$omp end parallel do
call dgemm('N','N',n,n*N3,n,-1.0d0,Pf,n,Sin,n,1.0d0,Tall,n)
!$omp parallel do private(yy)
do yy = 1,N3
   Pout(:,:,yy) = transpose(Tall(:,:,yy))
enddo
!$omp end parallel do
call dgetrs('N', n, n*N3, cosx_Snum_lu_hi, n, cosx_ipiv_hi, Pout, n, info)
!$omp parallel do private(yy)
do yy = 1,N3
   Tall(:,:,yy) = transpose(Pout(:,:,yy))
enddo
!$omp end parallel do
Pout = Tall
deallocate(dS1e, dH1e, Tall)
end subroutine hess_fit_py

module subroutine hess_snum_rng(a_lo, a_hi, Sout)
integer,intent(in) :: a_lo, a_hi
real(8),intent(out) :: Sout(n,n,3*(a_hi-a_lo+1))
real(8),allocatable :: xb(:,:), x1b(:,:,:), wx(:,:)
real(8),allocatable :: xc(:,:), xp(:,:), xm(:,:), cc(:,:,:), xm2(:)
integer,allocatable :: aol(:)
integer :: ib, nb, bstart, g, nloc, loc(natoms), yy, nl, kk, si, a, b, dd, A1, np_, nm_, nyb, ia
logical :: hasc(3)
real(8) :: dwv
integer,allocatable :: bord(:)
real(8),allocatable :: bcst(:)
integer :: ibo
integer(omp_lock_kind) :: slk(a_lo:a_hi)
nyb = 3*(a_hi-a_lo+1)
Sout = 0.0d0
do A1 = a_lo,a_hi
   call omp_init_lock(slk(A1))
enddo
allocate(bord(nblocks), bcst(nblocks))
do ib = 0,nblocks-1
   nl = 0
   do kk = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      nl = nl + cosx_shell_dim(cosx_blkidx_hi(kk))
   enddo
   bcst(ib+1) = dble(cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib))*dble(nl)**2
   bord(ib+1) = ib
enddo
do ib = 2,nblocks
   kk = bord(ib) ; dwv = bcst(ib) ; a = ib - 1
   do while (a .ge. 1)
      if (bcst(a) .ge. dwv) exit
      bord(a+1) = bord(a) ; bcst(a+1) = bcst(a)
      a = a - 1
   enddo
   bord(a+1) = kk ; bcst(a+1) = dwv
enddo
!$omp parallel default(shared) private(ib,ibo,nb,bstart,g,nloc,loc,yy,xb,x1b,wx,xc,xp,xm,cc,xm2,aol,nl,kk,si,a,b,dd,A1) &
!$omp&  private(np_,nm_,dwv,hasc,ia)
allocate(aol(n))
!$omp do schedule(dynamic)
do ibo = 1,nblocks
   ib = bord(ibo)
   bstart = cosx_blkstart_hi(ib) + 1
   nb = cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib)
   allocate(xb(n,nb), x1b(n,3,nb), wx(nb,n))
   call hess_block_ao(ib, nb, bstart, xb, x1b)
   call hess_block_atoms(ib, nloc, loc)
   do g = 1,nb
      wx(g,:) = cosx_w_hi(bstart+g-1)*xb(:,g)
   enddo
   nl = 0
   do kk = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
      si = cosx_blkidx_hi(kk)
      do a = 1,cosx_shell_dim(si)
         nl = nl + 1
         aol(nl) = cosx_ao_offset(si) + a
      enddo
   enddo
   allocate(xc(nl,nb), xp(nl,nb), xm(nl,nb), cc(nl,nl,3), xm2(nb))
   do g = 1,nb
      xc(:,g) = xb(aol(1:nl),g)
      xm2(g) = maxval(abs(xc(:,g)))**2
   enddo
   do ia = 0, a_hi-a_lo
      A1 = a_lo + mod(ia + ib, a_hi-a_lo+1)
      hasc = .false.
      do dd = 1,3
         np_ = 0 ; nm_ = 0
         do g = 1,nb
            dwv = dwdr(dd,A1,bstart+g-1)
            if (abs(dwv)*xm2(g) .lt. 1.0d-14) cycle
            if (dwv .gt. 0.0d0) then
               np_ = np_ + 1
               xp(:,np_) = xc(:,g)*sqrt(dwv)
            else if (dwv .lt. 0.0d0) then
               nm_ = nm_ + 1
               xm(:,nm_) = xc(:,g)*sqrt(-dwv)
            endif
         enddo
         if (np_ + nm_ .eq. 0) cycle
         hasc(dd) = .true.
         cc(:,:,dd) = 0.0d0
         if (np_ .gt. 0) call dsyrk('U','N',nl,np_,1.0d0,xp,nl,0.0d0,cc(1,1,dd),nl)
         if (nm_ .gt. 0) call dsyrk('U','N',nl,nm_,-1.0d0,xm,nl,1.0d0,cc(1,1,dd),nl)
      enddo
      if (.not. any(hasc)) cycle
      call omp_set_lock(slk(A1))
      do dd = 1,3
         if (.not. hasc(dd)) cycle
         yy = dd + 3*(A1-a_lo)
         do b = 1,nl
            do a = 1,b
               Sout(aol(a),aol(b),yy) = Sout(aol(a),aol(b),yy) + merge(0.5d0*cc(a,b,dd), cc(a,b,dd), a .eq. b)
            enddo
         enddo
      enddo
      call omp_unset_lock(slk(A1))
   enddo
   deallocate(xc, xp, xm, cc, xm2)
   call hess_rows_rng(nb, bstart, nloc, loc, x1b, wx, Sout, a_lo, a_hi, slk)
   deallocate(xb, x1b, wx)
enddo
!$omp end do
deallocate(aol)
!$omp end parallel
!$omp parallel do private(yy) schedule(static)
do yy = 1,nyb
   Sout(:,:,yy) = Sout(:,:,yy) + transpose(Sout(:,:,yy))
enddo
!$omp end parallel do
do A1 = a_lo,a_hi
   call omp_destroy_lock(slk(A1))
enddo
deallocate(bord, bcst)
end subroutine hess_snum_rng

module subroutine hess_rows_rng(nb, bstart, nloc, loc, x1b, Rm, Zout, a_lo, a_hi, slk)
integer,intent(in) :: nb, bstart, nloc, loc(natoms), a_lo, a_hi
real(8),intent(in) :: x1b(n,3,nb), Rm(nb,n)
real(8),intent(inout) :: Zout(n,n,3*(a_hi-a_lo+1))
integer(omp_lock_kind),intent(inout) :: slk(a_lo:a_hi)
real(8),allocatable :: xg(:,:), rg(:,:), tmp(:,:,:)
integer :: ee, il, Y, g, pg, np_, idx(nb), pp
logical :: done(natoms)
allocate(tmp(n,n,3))
do il = 1,nloc
   Y = loc(il)
   if (Y .lt. a_lo .or. Y .gt. a_hi) cycle
   do ee = 1,3
      call dgemm('N','N',at_n(Y),n,nb,-1.0d0,x1b(at_off(Y)+1,ee,1),3*n,Rm,nb,0.0d0,tmp(1,1,ee),n)
   enddo
   call omp_set_lock(slk(Y))
   do ee = 1,3
      Zout(at_off(Y)+1:at_off(Y)+at_n(Y),:,ee+3*(Y-a_lo)) = Zout(at_off(Y)+1:at_off(Y)+at_n(Y),:,ee+3*(Y-a_lo)) &
                                                             + tmp(1:at_n(Y),:,ee)
   enddo
   call omp_unset_lock(slk(Y))
enddo
allocate(xg(n,nb), rg(nb,n))
done = .false.
do g = 1,nb
   pp = cosx_atom_hi(bstart+g-1)
   if (done(pp)) cycle
   done(pp) = .true.
   if (pp .lt. a_lo .or. pp .gt. a_hi) cycle
   np_ = 0
   do pg = g,nb
      if (cosx_atom_hi(bstart+pg-1) .eq. pp) then
         np_ = np_ + 1
         idx(np_) = pg
      endif
   enddo
   do pg = 1,np_
      rg(pg,:) = Rm(idx(pg),:)
   enddo
   do ee = 1,3
      do pg = 1,np_
         xg(:,pg) = x1b(:,ee,idx(pg))
      enddo
      call dgemm('N','N',n,n,np_,1.0d0,xg,n,rg,nb,0.0d0,tmp(1,1,ee),n)
   enddo
   call omp_set_lock(slk(pp))
   Zout(:,:,3*(pp-a_lo)+1:3*(pp-a_lo)+3) = Zout(:,:,3*(pp-a_lo)+1:3*(pp-a_lo)+3) + tmp
   call omp_unset_lock(slk(pp))
enddo
deallocate(xg, rg, tmp)
end subroutine hess_rows_rng

module subroutine hess_fit_py_rng(Sin, Pout, a_lo, a_hi)
integer,intent(in) :: a_lo, a_hi
real(8),intent(in) :: Sin(n,n,3*(a_hi-a_lo+1))
real(8),intent(out) :: Pout(n,n,3*(a_hi-a_lo+1))
real(8),allocatable :: Tall(:,:,:)
integer :: yl, A, dd, k, l, nyb
nyb = 3*(a_hi-a_lo+1)
allocate(Tall(n,n,nyb))
!$omp parallel do private(yl,A,dd,k,l)
do yl = 1,nyb
   A = a_lo + (yl-1)/3
   dd = yl - 3*((yl-1)/3)
   Tall(:,:,yl) = 0.0d0
   do l = 1,n
      do k = 1,n
         if (ao_atom(k) .eq. A) Tall(k,l,yl) = Tall(k,l,yl) - dS1c(k,l,dd)
         if (ao_atom(l) .eq. A) Tall(k,l,yl) = Tall(k,l,yl) - dS1c(l,k,dd)
      enddo
   enddo
enddo
!$omp end parallel do
call dgemm('N','N',n,n*nyb,n,-1.0d0,Pf,n,Sin,n,1.0d0,Tall,n)
!$omp parallel do private(yl)
do yl = 1,nyb
   Pout(:,:,yl) = transpose(Tall(:,:,yl))
enddo
!$omp end parallel do
call dgetrs('N', n, n*nyb, cosx_Snum_lu_hi, n, cosx_ipiv_hi, Pout, n, info)
!$omp parallel do private(yl)
do yl = 1,nyb
   Tall(:,:,yl) = transpose(Pout(:,:,yl))
enddo
!$omp end parallel do
Pout = Tall
deallocate(Tall)
end subroutine hess_fit_py_rng

end submodule exchange_cosx_hess_setup_impl
