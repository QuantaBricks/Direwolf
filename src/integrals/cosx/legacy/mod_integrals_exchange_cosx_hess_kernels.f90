! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic Hessian - per-block kernels: AO values, tangents, row accumulations, dK/dR legs, shell maxima, sub-batch post.

submodule (mod_integrals:exchange_cosx_hess_impl) exchange_cosx_hess_kernels_impl
implicit none
contains

module subroutine hmem(tag)
use mod_meminfo, only: get_process_memory_bytes
character(len=*),intent(in) :: tag
integer(8) :: rss, hwm
integer :: u, ios
if (.not. verbose) return
call get_process_memory_bytes(rss, hwm)
print '(A,A,2F8.3)', " COSX HESS mem (GB) now / peak since last: ", tag, dble(rss)/1.0d9, dble(hwm)/1.0d9
open(newunit=u, file="/proc/self/clear_refs", status="old", action="write", iostat=ios)
if (ios .eq. 0) then
   write(u,'(A)',iostat=ios) "5"
   close(u)
endif
end subroutine hmem

module subroutine hess_block_ao(ib, nb, bstart, xb, x1b, x2b)
integer,intent(in) :: ib, nb, bstart
real(8),intent(out) :: xb(n,nb), x1b(n,3,nb)
real(8),intent(out),optional :: x2b(n,3,3,nb)
integer :: sidx, si
xb = 0.0d0
x1b = 0.0d0
if (present(x2b)) x2b = 0.0d0
do sidx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
   si = cosx_blkidx_hi(sidx)
   if (present(x2b)) then
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
           cosx_coor_hi(:,bstart:bstart+nb-1), nb, cosx_ao_offset(si), n, xb, x1b, x2b)
   else
      call GTOeval_shell_batch(cosx_shell_atom(si), cosx_shell_local(si), &
           cosx_coor_hi(:,bstart:bstart+nb-1), nb, cosx_ao_offset(si), n, xb, x1b)
   endif
enddo
end subroutine hess_block_ao

module subroutine hess_block_atoms(ib, nloc, loc)
integer,intent(in) :: ib
integer,intent(out) :: nloc, loc(natoms)
logical :: hit(natoms)
integer :: sidx
hit = .false.
do sidx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
   hit(cosx_shell_atom(cosx_blkidx_hi(sidx))) = .true.
enddo
nloc = 0
do sidx = 1,natoms
   if (hit(sidx)) then
      nloc = nloc + 1
      loc(nloc) = sidx
   endif
enddo
end subroutine hess_block_atoms

module subroutine hess_tangent_lr(nb, bstart, nloc, loc, x1b, M, L, R, T)
integer,intent(in) :: nb, bstart, nloc, loc(natoms)
real(8),intent(in) :: x1b(n,3,nb), M(n,n), L(n,nop), R(n,nop)
real(8),intent(out) :: T(nb,N3,n)
real(8),allocatable :: full(:,:), xl(:,:)
integer :: il, Y, ee, g, yy, pg
T = 0.0d0
allocate(full(nb,n), xl(nb,nop))
do ee = 1,3
   do il = 1,nloc
      Y = loc(il)
      yy = ee + 3*(Y-1)
      call dgemm('T','N',nb,n,at_n(Y),-1.0d0,x1b(at_off(Y)+1,ee,1),3*n, &
                 M(at_off(Y)+1,1),n,0.0d0,T(1,yy,1),nb*N3)
   enddo
   call dgemm('T','N',nb,nop,n,1.0d0,x1b(1,ee,1),3*n,L,n,0.0d0,xl,nb)
   call dgemm('N','T',nb,n,nop,1.0d0,xl,nb,R,n,0.0d0,full,nb)
   do g = 1,nb
      pg = cosx_atom_hi(bstart+g-1)
      yy = ee + 3*(pg-1)
      T(g,yy,:) = T(g,yy,:) + full(g,:)
   enddo
enddo
deallocate(full, xl)
end subroutine hess_tangent_lr

module subroutine hess_tangent_occ(nb, bstart, nloc, loc, x1b, L, T)
integer,intent(in) :: nb, bstart, nloc, loc(natoms)
real(8),intent(in) :: x1b(n,3,nb), L(n,nop)
real(8),intent(inout) :: T(nb,nop,N3)
real(8),allocatable :: xl(:,:)
integer :: il, Y, ee, g, yy
allocate(xl(nb,nop))
do ee = 1,3
   do il = 1,nloc
      Y = loc(il)
      yy = ee + 3*(Y-1)
      call dgemm('T','N',nb,nop,at_n(Y),-1.0d0,x1b(at_off(Y)+1,ee,1),3*n, &
                 L(at_off(Y)+1,1),n,1.0d0,T(1,1,yy),nb)
   enddo
   call dgemm('T','N',nb,nop,n,1.0d0,x1b(1,ee,1),3*n,L,n,0.0d0,xl,nb)
   do g = 1,nb
      yy = ee + 3*(cosx_atom_hi(bstart+g-1)-1)
      T(g,:,yy) = T(g,:,yy) + xl(g,:)
   enddo
enddo
deallocate(xl)
end subroutine hess_tangent_occ

module subroutine hess_tangent(nb, bstart, nloc, loc, x1b, M, T)
integer,intent(in) :: nb, bstart, nloc, loc(natoms)
real(8),intent(in) :: x1b(n,3,nb), M(n,n)
real(8),intent(out) :: T(nb,N3,n)
real(8),allocatable :: full(:,:)
integer :: il, Y, ee, g, yy, pg
T = 0.0d0
allocate(full(nb,n))
do ee = 1,3
   do il = 1,nloc
      Y = loc(il)
      yy = ee + 3*(Y-1)
      call dgemm('T','N',nb,n,at_n(Y),-1.0d0,x1b(at_off(Y)+1,ee,1),3*n, &
                 M(at_off(Y)+1,1),n,0.0d0,T(1,yy,1),nb*N3)
   enddo
   call dgemm('T','N',nb,n,n,1.0d0,x1b(1,ee,1),3*n,M,n,0.0d0,full,nb)
   do g = 1,nb
      pg = cosx_atom_hi(bstart+g-1)
      yy = ee + 3*(pg-1)
      T(g,yy,:) = T(g,yy,:) + full(g,:)
   enddo
enddo
deallocate(full)
end subroutine hess_tangent

module subroutine hess_xxy_contract(nb, bstart, ib, x2b, wv, Rv, Hloc, Hp, Hpp)
integer,intent(in) :: nb, bstart, ib
real(8),intent(in) :: x2b(n,3,3,nb), wv(nb), Rv(nb,n)
real(8),intent(inout) :: Hloc(N3,N3), Hp(nb,N3,3), Hpp(nb,3,3)
integer :: sidx, si, Y, kk, kap, dd, ee2, g
real(8) :: t, s
do sidx = cosx_blkptr_hi(ib)+1, cosx_blkptr_hi(ib+1)
   si = cosx_blkidx_hi(sidx)
   Y = cosx_shell_atom(si)
   do kk = 1,cosx_shell_dim(si)
      kap = cosx_ao_offset(si) + kk
      do ee2 = 1,3
         do dd = 1,3
            s = 0.0d0
            do g = 1,nb
               t = wv(g)*x2b(kap,dd,ee2,g)*Rv(g,kap)
               s = s + t
               Hp(g,dd+3*(Y-1),ee2) = Hp(g,dd+3*(Y-1),ee2) - t
               Hpp(g,dd,ee2) = Hpp(g,dd,ee2) + t
            enddo
            Hloc(dd+3*(Y-1),ee2+3*(Y-1)) = Hloc(dd+3*(Y-1),ee2+3*(Y-1)) + s
         enddo
      enddo
   enddo
enddo
end subroutine hess_xxy_contract

module subroutine hess_scatter_par(nb, bstart, Hp, Hpp, Hloc)
integer,intent(in) :: nb, bstart
real(8),intent(in) :: Hp(nb,N3,3), Hpp(nb,3,3)
real(8),intent(inout) :: Hloc(N3,N3)
integer :: g, pg, yy, ee
do g = 1,nb
   pg = cosx_atom_hi(bstart+g-1)
   do ee = 1,3
      do yy = 1,N3
         Hloc(yy, ee+3*(pg-1)) = Hloc(yy, ee+3*(pg-1)) + Hp(g,yy,ee)
         Hloc(ee+3*(pg-1), yy) = Hloc(ee+3*(pg-1), yy) + Hp(g,yy,ee)
      enddo
   enddo
   Hloc(3*(pg-1)+1:3*pg, 3*(pg-1)+1:3*pg) = Hloc(3*(pg-1)+1:3*pg, 3*(pg-1)+1:3*pg) + Hpp(g,:,:)
enddo
end subroutine hess_scatter_par

module subroutine hess_rows_par(nb, bstart, nloc, loc, x1b, Rm, Zout)
integer,intent(in) :: nb, bstart, nloc, loc(natoms)
real(8),intent(in) :: x1b(n,3,nb), Rm(nb,n)
real(8),intent(inout) :: Zout(n,n,N3)
real(8),allocatable :: Z1(:,:), xg(:,:), rg(:,:)
integer :: ee, il, Y, g, pg, np_, idx(nb), pp
logical :: done(natoms)
allocate(Z1(n,n), xg(n,nb), rg(nb,n))
do ee = 1,3
   call dgemm('N','N',n,n,nb,1.0d0,x1b(1,ee,1),3*n,Rm,nb,0.0d0,Z1,n)
   do il = 1,nloc
      Y = loc(il)
      Zout(at_off(Y)+1:at_off(Y)+at_n(Y),:,ee+3*(Y-1)) = &
         Zout(at_off(Y)+1:at_off(Y)+at_n(Y),:,ee+3*(Y-1)) - Z1(at_off(Y)+1:at_off(Y)+at_n(Y),:)
   enddo
enddo
done = .false.
do g = 1,nb
   pp = cosx_atom_hi(bstart+g-1)
   if (done(pp)) cycle
   done(pp) = .true.
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
      call dgemm('N','N',n,n,np_,1.0d0,xg,n,rg,nb,1.0d0,Zout(1,1,ee+3*(pp-1)),n)
   enddo
enddo
deallocate(Z1, xg, rg)
end subroutine hess_rows_par

module subroutine hess_rows_t7(nb, bstart, nloc, loc, x1b, Rm, T7)
integer,intent(in) :: nb, bstart, nloc, loc(natoms)
real(8),intent(in) :: x1b(n,3,nb), Rm(nb,n)
real(8),intent(inout) :: T7(N3,N3)
real(8),allocatable :: ZT(:,:,:), xg(:,:), rg(:,:)
integer :: ee, il, Y, g, pg, np_, idx(nb), pp
logical :: done(natoms)
allocate(ZT(n,n,3), xg(n,nb), rg(nb,n))
do ee = 1,3
   call dgemm('T','T',n,n,nb,1.0d0,Rm,nb,x1b(1,ee,1),3*n,0.0d0,ZT(1,1,ee),n)
enddo
do il = 1,nloc
   Y = loc(il)
   call dgemm('T','N',N3,3,n*at_n(Y),-1.0d0,Py(1,at_off(Y)+1,1),n*n,ZT(1,at_off(Y)+1,1),n*n, &
              1.0d0,T7(1,3*(Y-1)+1),N3)
enddo
done = .false.
do g = 1,nb
   pp = cosx_atom_hi(bstart+g-1)
   if (done(pp)) cycle
   done(pp) = .true.
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
      call dgemm('T','T',n,n,np_,1.0d0,rg,nb,xg,n,0.0d0,ZT(1,1,ee),n)
   enddo
   call dgemm('T','N',N3,3,n*n,1.0d0,Py,n*n,ZT,n*n,1.0d0,T7(1,3*(pp-1)+1),N3)
enddo
deallocate(ZT, xg, rg)
end subroutine hess_rows_t7

module subroutine hess_rows_shared(nb, bstart, is, x1b, Rm, Z1)
integer,intent(in) :: nb, bstart, is
real(8),intent(in) :: x1b(n,3,nb), Rm(nb,n)
real(8),intent(inout) :: Z1(n,n,3)
real(8),allocatable :: Zp(:,:,:), xg(:,:), rg(:,:)
integer :: ee, g, pg, np_, idx(nb), pp
logical :: done(natoms)
allocate(Zp(n,n,3), xg(n,nb), rg(nb,n))
do ee = 1,3
   call dgemm('N','N',n,n,nb,1.0d0,x1b(1,ee,1),3*n,Rm,nb,1.0d0,Z1(1,1,ee),n)
enddo
done = .false.
do g = 1,nb
   pp = cosx_atom_hi(bstart+g-1)
   if (done(pp)) cycle
   done(pp) = .true.
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
      call dgemm('N','N',n,n,np_,1.0d0,xg,n,rg,nb,0.0d0,Zp(1,1,ee),n)
   enddo
   call omp_set_lock(zlock(pp))
   Zy(:,:,3*(pp-1)+1:3*pp,is) = Zy(:,:,3*(pp-1)+1:3*pp,is) + Zp
   call omp_unset_lock(zlock(pp))
enddo
deallocate(Zp, xg, rg)
end subroutine hess_rows_shared

module subroutine hess_rows_shared_c(nb, bstart, is, x1b, Rm, Z1)
integer,intent(in) :: nb, bstart, is
real(8),intent(in) :: x1b(n,3,nb), Rm(nb,nop)
real(8),intent(inout) :: Z1(n,nop,3)
real(8),allocatable :: Zp(:,:,:), xg(:,:), rg(:,:)
integer :: ee, g, pg, np_, idx(nb), pp
logical :: done(natoms)
allocate(Zp(n,nop,3), xg(n,nb), rg(nb,nop))
do ee = 1,3
   call dgemm('N','N',n,nop,nb,1.0d0,x1b(1,ee,1),3*n,Rm,nb,1.0d0,Z1(1,1,ee),n)
enddo
done = .false.
do g = 1,nb
   pp = cosx_atom_hi(bstart+g-1)
   if (done(pp)) cycle
   done(pp) = .true.
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
      call dgemm('N','N',n,nop,np_,1.0d0,xg,n,rg,nb,0.0d0,Zp(1,1,ee),n)
   enddo
   call omp_set_lock(zlock(pp))
   Zt(:,:,3*(pp-1)+1:3*pp,is) = Zt(:,:,3*(pp-1)+1:3*pp,is) + Zp
   call omp_unset_lock(zlock(pp))
enddo
deallocate(Zp, xg, rg)
end subroutine hess_rows_shared_c

module subroutine kx_free_leg_max(nb, bstart, xb, x1b, wv, qv, Lpt)
integer,intent(in) :: nb, bstart
real(8),intent(in) :: xb(n,nb), x1b(n,3,nb), wv(nb), qv(nb,n)
real(8),intent(out) :: Lpt(nb)
real(8) :: mq(nb), dwm, x1m
integer :: g, i, yy
mq = 0.0d0
do i = 1,n
   do g = 1,nb
      mq(g) = max(mq(g), abs(qv(g,i)))
   enddo
enddo
do g = 1,nb
   if (var_k) mq(g) = max(mq(g), maxval(abs(xb(:,g))))
   dwm = 0.0d0
   do yy = 1,N3
      dwm = max(dwm, abs(dwdr(mod(yy-1,3)+1,(yy-1)/3+1,bstart+g-1)))
   enddo
   x1m = maxval(abs(x1b(:,:,g)))
   Lpt(g) = max((abs(wv(g)) + dwm)*mq(g), abs(wv(g))*x1m*pcol)
   if (fit) Lpt(g) = max(Lpt(g), abs(wv(g))*maxval(abs(xb(:,g)))*pymax)
enddo
end subroutine kx_free_leg_max

module subroutine kx_proj_leg(nb, bstart, is, xb, wv, wq, qq, x1b, Mm, MC, nloc, loc, Gl, kk, fitleg, &
                       QVc, GC, WC, wGC, Q1C, KCa, RR, RS, XP, Kloc, XWloc, tk, usexk, XK, kcg)
integer,intent(in) :: nb, bstart, is, nloc, loc(natoms)
real(8),intent(in) :: xb(n,nb), wv(nb), wq(nb,n), x1b(n,3,nb), Mm(n,n), MC(n,nop)
real(8),intent(in) :: qq(nb,n), Gl(nb,n), kk(nb,N3,*)
logical,intent(in) :: fitleg
real(8),intent(inout) :: QVc(nb,nop), GC(nb,nop), WC(nb,nop), wGC(nb,nop), Q1C(nb,N3,nop), KCa(nb,N3,nop)
real(8),intent(inout) :: RR(nb,nop,N3,3), RS(N3*n,nop), XP(*), Kloc(n,nop,N3,ns), XWloc(n,nop,*)
integer :: o, y, g2, i2, ee, il, Y1, pa, np_
integer :: plist(natoms)
logical :: pseen(natoms)
real(8),allocatable :: FC(:,:), XWe(:,:), wGCp(:,:), dwl(:,:), x1p(:,:)
integer :: npp
logical :: ident
real(8),intent(inout) :: tk(6)
logical,intent(in) :: usexk
real(8),intent(inout) :: XK(n,nop,N3,*)
logical,intent(in) :: kcg
real(8) :: tt0
integer :: nsig, ks(n), c2, y0, ny, yk, nyk, yq
real(8),allocatable :: xs(:,:), tmp(:,:)
allocate(FC(nb,nop), XWe(n,nop), wGCp(nb,nop), dwl(nb,N3), x1p(n,nb))
ident = .not. fitleg
do y = 1,N3
   do g2 = 1,nb
      dwl(g2,y) = dwdr(mod(y-1,3)+1,(y-1)/3+1,bstart+g2-1)
   enddo
enddo
np_ = 0
pseen = .false.
do g2 = 1,nb
   pa = cosx_atom_hi(bstart+g2-1)
   if (.not. pseen(pa)) then
      pseen(pa) = .true.
      np_ = np_ + 1
      plist(np_) = pa
   endif
enddo
tt0 = omp_get_wtime()
if (kxmir) call dgemm('N','N',nb,nop,n,1.0d0,qq,nb,Cp(1,1,is),n,0.0d0,QVc,nb)
call dgemm('N','N',nb,nop,n,1.0d0,Gl,nb,Cp(1,1,is),n,0.0d0,GC,nb)
if (kxmir) call dgemm('N','N',nb,nop,n,1.0d0,wq,nb,Cp(1,1,is),n,0.0d0,WC,nb)
tk(1) = tk(1) + (omp_get_wtime()-tt0) ; tt0 = omp_get_wtime()
if (kxmir) then
Q1C = 0.0d0
do ee = 1,3
   do il = 1,nloc
      Y1 = loc(il)
      call dgemm('T','N',nb,nop,at_n(Y1),-1.0d0,x1b(at_off(Y1)+1,ee,1),3*n, &
                 MC(at_off(Y1)+1,1),n,0.0d0,Q1C(1,ee+3*(Y1-1),1),nb*N3)
   enddo
   call dgemm('T','N',nb,nop,n,1.0d0,x1b(1,ee,1),3*n,MC,n,0.0d0,FC,nb)
   do g2 = 1,nb
      y = ee + 3*(cosx_atom_hi(bstart+g2-1)-1)
      Q1C(g2,y,:) = Q1C(g2,y,:) + FC(g2,:)
   enddo
enddo
endif
tk(2) = tk(2) + (omp_get_wtime()-tt0) ; tt0 = omp_get_wtime()
if (.not. kcg) call dgemm('N','N',nb*N3,nop,n,1.0d0,kk,nb*N3,Cp(1,1,is),n,0.0d0,KCa,nb*N3)
tk(3) = tk(3) + (omp_get_wtime()-tt0) ; tt0 = omp_get_wtime()
if (kxmir .and. fitleg .and. fit) call dgemm('T','N',nb,nop*N3,n,1.0d0,xb,n,PtC(1,1,1,is),n,0.0d0,XP,nb)
do o = 1,nop
   do g2 = 1,nb
      wGC(g2,o) = wv(g2)*GC(g2,o)
   enddo
enddo
do y = 1,N3
   do o = 1,nop
      if (usexk) then
         do g2 = 1,nb
            RR(g2,o,y,1) = dwl(g2,y)*GC(g2,o) + wv(g2)*KCa(g2,y,o)
         enddo
      else
      do g2 = 1,nb
         RR(g2,o,y,1) = dwl(g2,y)*GC(g2,o)
         RR(g2,o,y,3) = KCa(g2,y,o)
      enddo
      endif
      if (.not. kxmir) cycle
      do g2 = 1,nb
         RR(g2,o,y,2) = dwl(g2,y)*QVc(g2,o) + wv(g2)*Q1C(g2,y,o)
      enddo
      if (fitleg .and. fit) then
         do g2 = 1,nb
            RR(g2,o,y,2) = RR(g2,o,y,2) + wv(g2)*XP(g2+nb*((o-1)+nop*(y-1)))
         enddo
      endif
   enddo
enddo
tk(4) = tk(4) + (omp_get_wtime()-tt0) ; tt0 = omp_get_wtime()
if (kxd) then
   block
      real(8) :: lm(nb,4), rm, loc4(12)
      integer :: t4, it
      lm = 0.0d0
      do i2 = 1,n
         do g2 = 1,nb
            lm(g2,1) = max(lm(g2,1), abs(qq(g2,i2)))
            lm(g2,2) = max(lm(g2,2), abs(Gl(g2,i2)))
            lm(g2,3) = max(lm(g2,3), abs(wq(g2,i2)))
         enddo
      enddo
      do g2 = 1,nb
         lm(g2,4) = maxval(abs(WC(g2,:)))
      enddo
      loc4 = 0.0d0
      do y = 1,N3
         do g2 = 1,nb
            do t4 = 1,4
               if (t4 .le. 3) then
                  rm = lm(g2,t4)*maxval(abs(RR(g2,:,y,t4)))
               else
                  rm = lm(g2,4)*maxval(abs(kk(g2,y,1:n)))
               endif
               loc4(t4) = loc4(t4) + 1.0d0
               do it = 1,2
                  if (rm .ge. 10.0d0**(-7-it)) loc4(4*it+t4) = loc4(4*it+t4) + 1.0d0
               enddo
            enddo
         enddo
      enddo
      !$omp critical (kxdiag)
      kxdg = kxdg + loc4
      !$omp end critical (kxdiag)
   end block
endif
if (usexk) then
   nsig = 0
   do i2 = 1,n
      if (maxval(abs(xb(i2,:))) .ge. 1.0d-12) then
         nsig = nsig + 1
         ks(nsig) = i2
      endif
   enddo
   allocate(xs(nsig,nb), tmp(nsig,nop*8))
   do g2 = 1,nb
      xs(:,g2) = xb(ks(1:nsig),g2)
   enddo
   nyk = (N3+7)/8
   do yk = 0,nyk-1
      y0 = 8*yk
      if (xkshr) y0 = 8*mod(yk + omp_get_thread_num(), nyk)
      ny = min(8, N3-y0)
      call dgemm('N','N',nsig,nop*ny,nb,1.0d0,xs,nsig,RR(1,1,y0+1,1),nb,0.0d0,tmp,nsig)
      if (xkshr) then
         call omp_set_lock(xklk(y0/8+1))
         do y = y0+1,y0+ny
            do o = 1,nop
               c2 = o + nop*(y-y0-1)
               do i2 = 1,nsig
                  XKsum(ks(i2),o,y,is) = XKsum(ks(i2),o,y,is) + tmp(i2,c2)
               enddo
            enddo
         enddo
         call omp_unset_lock(xklk(y0/8+1))
      else
      do y = y0+1,y0+ny
         do o = 1,nop
            c2 = o + nop*(y-y0-1)
            do i2 = 1,nsig
               XK(ks(i2),o,y,is) = XK(ks(i2),o,y,is) + tmp(i2,c2)
            enddo
         enddo
      enddo
      endif
   enddo
   deallocate(xs, tmp)
else
call dgemm('T','N',n,nop*N3,nb,1.0d0,qq,nb,RR(1,1,1,1),nb,1.0d0,Kloc(1,1,1,is),n)
if (kxmir) call dgemm('T','N',n,nop*N3,nb,1.0d0,Gl,nb,RR(1,1,1,2),nb,1.0d0,Kloc(1,1,1,is),n)
call dgemm('T','N',n,nop*N3,nb,1.0d0,wq,nb,RR(1,1,1,3),nb,1.0d0,Kloc(1,1,1,is),n)
endif
tk(5) = tk(5) + (omp_get_wtime()-tt0) ; tt0 = omp_get_wtime()
do ee = 1,3
   call dgemm('N','N',n,nop,nb,1.0d0,x1b(1,ee,1),3*n,wGC,nb,0.0d0,XWe,n)
   do il = 1,nloc
      Y1 = loc(il)
      if (usexk .and. xkshr) then
         yq = ee+3*(Y1-1)
         call omp_set_lock(xklk((yq-1)/8+1))
         XKsum(at_off(Y1)+1:at_off(Y1)+at_n(Y1),:,yq,is) = &
            XKsum(at_off(Y1)+1:at_off(Y1)+at_n(Y1),:,yq,is) - XWe(at_off(Y1)+1:at_off(Y1)+at_n(Y1),:)
         call omp_unset_lock(xklk((yq-1)/8+1))
      else if (usexk) then
         XK(at_off(Y1)+1:at_off(Y1)+at_n(Y1),:,ee+3*(Y1-1),is) = &
            XK(at_off(Y1)+1:at_off(Y1)+at_n(Y1),:,ee+3*(Y1-1),is) - XWe(at_off(Y1)+1:at_off(Y1)+at_n(Y1),:)
      else
      call dgemm('T','N',n,nop,at_n(Y1),-1.0d0,Mm(at_off(Y1)+1,1),n,XWe(at_off(Y1)+1,1),n, &
                 1.0d0,Kloc(1,1,ee+3*(Y1-1),is),n)
      endif
   enddo
   do il = 1,np_
      pa = plist(il)
      npp = 0
      do g2 = 1,nb
         if (cosx_atom_hi(bstart+g2-1) .eq. pa) then
            npp = npp + 1
            x1p(:,npp) = x1b(:,ee,g2)
            wGCp(npp,:) = wGC(g2,:)
         endif
      enddo
      call dgemm('N','N',n,nop,npp,1.0d0,x1p,n,wGCp,nb,0.0d0,XWe,n)
      if (usexk .and. xkshr) then
         yq = ee+3*(pa-1)
         call omp_set_lock(xklk((yq-1)/8+1))
         XKsum(:,:,yq,is) = XKsum(:,:,yq,is) + XWe
         call omp_unset_lock(xklk((yq-1)/8+1))
      else if (usexk) then
         XK(:,:,ee+3*(pa-1),is) = XK(:,:,ee+3*(pa-1),is) + XWe
      else if (ident) then
         Kloc(:,:,ee+3*(pa-1),is) = Kloc(:,:,ee+3*(pa-1),is) + XWe
      else
         call dgemm('T','N',n,nop,n,1.0d0,Mm,n,XWe,n,1.0d0,Kloc(1,1,ee+3*(pa-1),is),n)
      endif
   enddo
enddo
if (kxmir) then
   do y = 1,N3
      call dgemm('T','N',n,nop,nb,1.0d0,kk(1,y,1),nb*N3,WC,nb,1.0d0,Kloc(1,1,y,is),n)
   enddo
endif
if (fitleg .and. fit) call dgemm('N','N',n,nop,nb,1.0d0,xb,n,wGC,nb,1.0d0,XWloc(1,1,is),n)
tk(6) = tk(6) + (omp_get_wtime()-tt0)
deallocate(FC, XWe, wGCp, dwl, x1p)
end subroutine kx_proj_leg

module subroutine hess_post_sub(ng, bst, ib, is, nloc, loc, xs, x1s, x2s, ws, avs, Gvs, Gtvs, qvs, KCs_, MCs_, &
                         XAs, A1s, B1s, Hl, gl, Zl, XWl, Kl, XKl, tk, tsx)
integer,intent(in) :: ng, bst, ib, is, nloc, loc(natoms)
real(8),intent(in) :: xs(n,ng), x1s(n,3,ng), x2s(n,3,3,ng), ws(ng), avs(ng,n), Gvs(ng,n), Gtvs(ng,n)
real(8),intent(in) :: qvs(ng,n), MCs_(ng,N3,nop), XAs(ng,nop), A1s(ng,nop,N3), B1s(ng,nop,N3)
real(8),intent(inout) :: KCs_(ng,N3,nop)
real(8),intent(inout) :: Hl(N3,N3,NT), gl(N3), Zl(n,nop,3), XWl(n,nop,*), Kl(n,nop,N3,*), XKl(n,nop,N3,*)
real(8),intent(inout) :: tk(6), tsx(10)
real(8),allocatable :: QVc(:,:), GC(:,:), WC(:,:), wGC(:,:), Q1C(:,:,:), KCa(:,:,:), RR(:,:,:,:), RS(:,:), XP(:)
real(8),allocatable :: wq(:,:), epsy(:,:), GCt(:,:), KT(:,:,:), wA(:,:,:), hv(:,:), htv(:,:), Psi(:,:), HpT(:,:,:), HppT(:,:,:)
integer :: g, i, yy, oo
real(8) :: t0
t0 = omp_get_wtime()
allocate(QVc(ng,nop), GC(ng,nop), WC(ng,nop), wGC(ng,nop), Q1C(ng,N3,nop), KCa(ng,N3,nop))
allocate(RR(ng,nop,N3,3), RS(N3*n,nop), XP(1), wq(ng,n))
do i = 1,n
   wq(:,i) = ws(:)*qvs(:,i)
enddo
call kx_proj_leg(ng, bst, is, xs, ws, wq, qvs, x1s, Pt, PtCp(:,:,is), nloc, loc, &
                 Gvs, KCs_, .true., QVc, GC, WC, wGC, Q1C, KCs_, RR, RS, XP, Kl, XWl, tk, .true., XKl, .true.)
deallocate(QVc, GC, WC, wGC, Q1C, KCa, RR, RS, XP, wq)
tsx(5) = tsx(5) + (omp_get_wtime() - t0) ; t0 = omp_get_wtime()
allocate(epsy(ng,N3), GCt(ng,nop))
do g = 1,ng
   epsall(bst+g-1,is) = dot_product(avs(g,:), Gvs(g,:))
enddo
call dgemm('N','N',ng,nop,n,1.0d0,Gvs,ng,Cp(1,1,is),n,0.0d0,GCt,ng)
epsy = 0.0d0
do yy = 1,N3
   do oo = 1,nop
      do g = 1,ng
         epsy(g,yy) = epsy(g,yy) + A1s(g,oo,yy)*GCt(g,oo) + XAs(g,oo)*KCs_(g,yy,oo)
      enddo
   enddo
enddo
do yy = 1,N3
   do g = 1,ng
      gl(yy) = gl(yy) + dwdr(mod(yy-1,3)+1,(yy-1)/3+1,bst+g-1)*epsall(bst+g-1,is) + ws(g)*epsy(g,yy)
   enddo
enddo
call dgemm('N','N',N3,N3,ng,1.0d0,dwdr(1,1,bst),N3,epsy,ng,1.0d0,Hl(1,1,2),N3)
deallocate(epsy)
tsx(6) = tsx(6) + (omp_get_wtime() - t0) ; t0 = omp_get_wtime()
allocate(KT(ng,nop,N3), wA(ng,nop,N3))
do yy = 1,N3
   do oo = 1,nop
      do g = 1,ng
         KT(g,oo,yy) = KCs_(g,yy,oo)
         wA(g,oo,yy) = ws(g)*A1s(g,oo,yy)
      enddo
   enddo
enddo
call dgemm('T','N',N3,N3,ng*nop,1.0d0,wA,ng*nop,KT,ng*nop,1.0d0,Hl(1,1,3),N3)
do yy = 1,N3
   do oo = 1,nop
      do g = 1,ng
         KT(g,oo,yy) = ws(g)*MCs_(g,yy,oo)
      enddo
   enddo
enddo
call dgemm('T','N',N3,N3,ng*nop,1.0d0,KT,ng*nop,B1s,ng*nop,1.0d0,Hl(1,1,3),N3)
deallocate(KT, wA)
tsx(7) = tsx(7) + (omp_get_wtime() - t0) ; t0 = omp_get_wtime()
allocate(hv(ng,n), htv(ng,n), Psi(ng,n), HpT(ng,N3,3), HppT(ng,3,3))
call dgemm('N','N',ng,n,n,1.0d0,Gvs,ng,Dsp(:,:,is),n,0.0d0,hv,ng)
call dgemm('N','N',ng,n,n,1.0d0,Gtvs,ng,Dsp(:,:,is),n,0.0d0,htv,ng)
if (fit) then
   Psi = htv
   call dgemm('N','N',ng,n,n,1.0d0,hv,ng,Pf,n,1.0d0,Psi,ng)
else
   Psi = hv + htv
endif
HpT = 0.0d0 ; HppT = 0.0d0
call hess_xxy_contract(ng, bst, ib, x2s, ws, Psi, Hl(:,:,4), HpT, HppT)
call hess_scatter_par(ng, bst, HpT, HppT, Hl(:,:,4))
deallocate(hv, htv, Psi, HpT, HppT)
if (fit) then
   do oo = 1,nop
      GCt(:,oo) = ws(:)*GCt(:,oo)
   enddo
   call hess_rows_shared_c(ng, bst, is, x1s, GCt, Zl)
endif
deallocate(GCt)
tsx(8) = tsx(8) + (omp_get_wtime() - t0)
end subroutine hess_post_sub

module subroutine a1_shell_max(nb, is, sw, A1, Tu, Tmax)
integer,intent(in) :: nb, is
real(8),intent(in) :: sw(nb), A1(nb,nop,N3)
real(8),intent(inout) :: Tu(0:nBases-1), Tmax(0:nBases-1)
integer,parameter :: NBK = 60, NRC = 256
real(8),allocatable :: rn(:), Rc(:,:), Cg(:,:), Oc(:,:)
integer,allocatable :: ord(:), bk(:), act(:)
integer :: cnt(0:NBK+1), pos(0:NBK+1), nr, r, g, yy, i, c, k, m, p, bn, emax, nact, nkeep, s
real(8) :: bw(0:NBK+1), bu(0:NBK+1), v
nr = nb*N3
allocate(rn(nr), bk(nr), ord(nr))
emax = -huge(1)
do yy = 1,N3
   do g = 1,nb
      r = g + nb*(yy-1)
      rn(r) = sqrt(sum(A1(g,:,yy)**2))
      if (rn(r) .gt. 0.0d0) emax = max(emax, exponent(rn(r)))
   enddo
enddo
if (emax .eq. -huge(1)) return
cnt = 0 ; bw = 0.0d0
do r = 1,nr
   if (rn(r) .gt. 0.0d0) then
      bk(r) = min(NBK, emax - exponent(rn(r)))
   else
      bk(r) = NBK + 1
   endif
   cnt(bk(r)) = cnt(bk(r)) + 1
   bw(bk(r)) = max(bw(bk(r)), sw(mod(r-1,nb)+1)*rn(r))
enddo
pos(0) = 0
do k = 1,NBK+1
   pos(k) = pos(k-1) + cnt(k-1)
enddo
do r = 1,nr
   pos(bk(r)) = pos(bk(r)) + 1
   ord(pos(bk(r))) = r
enddo
do k = 0,NBK
   bu(k) = 2.0d0**(emax-k)
enddo
bu(NBK+1) = 0.0d0 ; bw(NBK+1) = 0.0d0
do k = NBK,0,-1
   bw(k) = max(bw(k), bw(k+1))
enddo
allocate(act(n))
nact = 0
do i = 1,n
   if (cnrm(i,is) .gt. 0.0d0) then
      nact = nact + 1
      act(nact) = i
   endif
enddo
allocate(Rc(NRC+nr/4,nop), Cg(n,nop), Oc(NRC+nr/4,n))
p = 0
bn = 0
do while (bn .le. NBK)
   nkeep = 0
   do c = 1,nact
      i = act(c)
      s = ao_sh(i)
      if (bu(bn)*cnrm(i,is) .gt. Tu(s) .or. bw(bn)*cnrm(i,is) .gt. Tmax(s)) then
         nkeep = nkeep + 1
         act(nkeep) = i
      endif
   enddo
   nact = nkeep
   if (nact .eq. 0) exit
   m = 0
   do while (bn .le. NBK .and. m .lt. NRC)
      m = m + cnt(bn)
      bn = bn + 1
   enddo
   if (m .eq. 0) cycle
   if (m .gt. size(Rc,1)) then
      deallocate(Rc, Oc)
      allocate(Rc(m,nop), Oc(m,n))
   endif
   do k = 1,m
      r = ord(p+k)
      Rc(k,:) = A1(mod(r-1,nb)+1,:,(r-1)/nb+1)
   enddo
   do c = 1,nact
      Cg(c,:) = Cp(act(c),1:nop,is)
   enddo
   call dgemm('N','T',m,nact,nop,1.0d0,Rc,size(Rc,1),Cg,n,0.0d0,Oc,size(Oc,1))
   do c = 1,nact
      s = ao_sh(act(c))
      do k = 1,m
         v = abs(Oc(k,c))
         Tu(s) = max(Tu(s), v)
         Tmax(s) = max(Tmax(s), sw(mod(ord(p+k)-1,nb)+1)*v)
      enddo
   enddo
   p = p + m
enddo
deallocate(rn, bk, ord, act, Rc, Cg, Oc)
end subroutine a1_shell_max

end submodule exchange_cosx_hess_kernels_impl
