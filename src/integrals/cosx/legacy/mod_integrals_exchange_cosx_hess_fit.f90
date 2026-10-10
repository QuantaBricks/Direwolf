! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic Hessian - the overlap-fit second-order terms and the u pass.

submodule (mod_integrals:exchange_cosx_hess_impl) exchange_cosx_hess_fit_impl
implicit none
contains

module subroutine hess_fit_second()
real(8),allocatable :: s1aa(:,:,:), s1ab(:,:,:), h1aa(:,:,:), h1ab(:,:,:)
real(8),allocatable :: CT(:,:,:), Tm(:,:), Ysym(:,:)
real(8),allocatable :: hwu(:,:,:,:), Huu(:,:,:)
integer :: yy, k, l, A, B, dd, ee
real(8) :: blk(3,3), tf0
tf0 = omp_get_wtime()
allocate(s1aa(n,n,9), s1ab(n,n,9), h1aa(n,n,9), h1ab(n,n,9), Ysym(n,n))
call integrals_hessian_1e_shell(n, s1aa, s1ab, h1aa, h1ab)
deallocate(h1aa, h1ab)
do is = 1,ns
   Ysym = 0.5d0*(Ym(:,:,is) + transpose(Ym(:,:,is)))
   do k = 1,n
      A = ao_atom(k)
      do l = 1,n
         B = ao_atom(l)
         blk = transpose(reshape(s1aa(k,l,:)*Ysym(k,l), (/3,3/)))
         Hterm(3*(A-1)+1:3*A, 3*(A-1)+1:3*A, 6, is) = Hterm(3*(A-1)+1:3*A, 3*(A-1)+1:3*A, 6, is) + blk
         blk = transpose(reshape(s1aa(l,k,:)*Ysym(k,l), (/3,3/)))
         Hterm(3*(B-1)+1:3*B, 3*(B-1)+1:3*B, 6, is) = Hterm(3*(B-1)+1:3*B, 3*(B-1)+1:3*B, 6, is) + blk
         blk = transpose(reshape(s1ab(k,l,:)*Ysym(k,l), (/3,3/)))
         Hterm(3*(A-1)+1:3*A, 3*(B-1)+1:3*B, 6, is) = Hterm(3*(A-1)+1:3*A, 3*(B-1)+1:3*B, 6, is) + blk
         Hterm(3*(B-1)+1:3*B, 3*(A-1)+1:3*A, 6, is) = Hterm(3*(B-1)+1:3*B, 3*(A-1)+1:3*A, 6, is) &
                                                      + transpose(blk)
      enddo
   enddo
enddo
deallocate(s1aa, s1ab, Ysym)
if (verbose) print '(A,F8.3)', " COSX HESS fit: T6a (1e 2nd derivs)  ", omp_get_wtime()-tf0

if (lrfit) then
   block
      real(8),allocatable :: SM(:,:,:), SR(:,:), Huu2(:,:)
      integer :: A2, d2, y2
      allocate(SM(n,nop,N3), SR(n,nop), Huu2(N3,N3))
      do is = 1,ns
         SM = 0.0d0
         call hess_sy_times(XWg(:,:,is), SM)
         Zt(:,:,:,is) = Zt(:,:,:,is) - SM
         call dgemm('T','N',N3,N3,n*nop,1.0d0,PtC(1,1,1,is),n*nop,Zt(1,1,1,is),n*nop,0.0d0,Huu2,N3)
         Hterm(:,:,7,is) = Hterm(:,:,7,is) + Huu2
         if (kproj) then
            call dgemm('N','N',n,nop*N3,n,-1.0d0,Pf,n,SM,n,1.0d0,Kx(1,1,1,is),n)
            do A2 = 1,natoms
               do d2 = 1,3
                  y2 = d2 + 3*(A2-1)
                  call dgemm('N','N',at_n(A2),nop,n,-1.0d0,dS1c(at_off(A2)+1,1,d2),n,XWg(1,1,is),n, &
                             1.0d0,Kx(at_off(A2)+1,1,y2,is),n)
                  call dgemm('T','N',n,nop,at_n(A2),-1.0d0,dS1c(at_off(A2)+1,1,d2),n,XWg(at_off(A2)+1,1,is),n, &
                             1.0d0,Kx(1,1,y2,is),n)
               enddo
            enddo
         endif
      enddo
      deallocate(SM, SR, Huu2)
   end block
   if (verbose) print '(A,F8.3)', " COSX HESS fit: + S^y R, P''/Z terms  ", omp_get_wtime()-tf0
   call hess_fit_u_pass()
   if (verbose) print '(A,F8.3)', " COSX HESS fit: + u pass            ", omp_get_wtime()-tf0
   return
endif
if (pstream) then
   block
      real(8),allocatable :: Sb(:,:,:), Pb(:,:,:), CTb(:,:,:), Hb(:,:), Wt(:,:)
      integer :: a_lo, a_hi, nyb, y0, yl
      allocate(Wt(n,n))
      do a_lo = 1,natoms,nab
         a_hi = min(natoms, a_lo+nab-1)
         nyb = 3*(a_hi-a_lo+1)
         y0 = 3*(a_lo-1)
         allocate(Sb(n,n,nyb), CTb(n,n,nyb))
         call hess_snum_rng(a_lo, a_hi, Sb)
         do is = 1,ns
            call dgemm('T','N',n,n*nyb,n,1.0d0,Ym(1,1,is),n,Sb,n,0.0d0,CTb,n)
            do yl = 1,nyb
               Wt = transpose(Zy(:,:,y0+yl,is))
               Zy(:,:,y0+yl,is) = Wt - CTb(:,:,yl)
            enddo
         enddo
         deallocate(Sb, CTb)
      enddo
      deallocate(Wt)
      if (verbose) print '(A,F8.3)', " COSX HESS fit: + W = Z^T - Y^T S^y ", omp_get_wtime()-tf0
      do a_lo = 1,natoms,nab
         a_hi = min(natoms, a_lo+nab-1)
         nyb = 3*(a_hi-a_lo+1)
         y0 = 3*(a_lo-1)
         allocate(Sb(n,n,nyb), Pb(n,n,nyb), Hb(nyb,N3))
         call hess_snum_rng(a_lo, a_hi, Sb)
         call hess_fit_py_rng(Sb, Pb, a_lo, a_hi)
         deallocate(Sb)
         do is = 1,ns
            call dgemm('T','N',nyb,N3,n*n,1.0d0,Pb,n*n,Zy(1,1,1,is),n*n,0.0d0,Hb,nyb)
            Hterm(y0+1:y0+nyb,:,7,is) = Hterm(y0+1:y0+nyb,:,7,is) + Hb
            if (kproj) then
               do yl = 1,nyb
                  call dgemm('N','N',n,nop,n,1.0d0,Pb(:,:,yl),n,XWg(1,1,is),n,1.0d0,Kx(1,1,y0+yl,is),n)
               enddo
            endif
         enddo
         deallocate(Pb, Hb)
      enddo
   end block
   if (verbose) print '(A,F8.3)', " COSX HESS fit: + P''/Z terms        ", omp_get_wtime()-tf0
   call hess_fit_u_pass()
   if (verbose) print '(A,F8.3)', " COSX HESS fit: + u pass            ", omp_get_wtime()-tf0
   return
endif
allocate(CT(n,n,N3), Tm(n,n))
do is = 1,ns
   call dgemm('T','N',n,n*N3,n,1.0d0,Ym(1,1,is),n,Sy,n,0.0d0,CT,n)
   allocate(Huu(N3,N3,1))
   call dgemm('T','N',N3,N3,n*n,1.0d0,Py,n*n,CT,n*n,0.0d0,Huu,N3)
   Hterm(:,:,6,is) = Hterm(:,:,6,is) - Huu(:,:,1) - transpose(Huu(:,:,1))
   if (lrz) then
      call dgemm('T','N',N3,N3,n*nop,1.0d0,PtC(1,1,1,is),n*nop,Zt(1,1,1,is),n*nop,0.0d0,Huu,N3)
      Hterm(:,:,7,is) = Hterm(:,:,7,is) + Huu(:,:,1)
   else if (.not. zt7) then
   do yy = 1,N3
      CT(:,:,yy) = transpose(Zy(:,:,yy,is))
   enddo
   call dgemm('T','N',N3,N3,n*n,1.0d0,Py,n*n,CT,n*n,0.0d0,Huu,N3)
   Hterm(:,:,7,is) = Hterm(:,:,7,is) + Huu(:,:,1)
   endif
   deallocate(Huu)
enddo
deallocate(CT, Tm)
if (verbose) print '(A,F8.3)', " COSX HESS fit: + P''/Z terms        ", omp_get_wtime()-tf0

call hess_fit_u_pass()
if (verbose) print '(A,F8.3)', " COSX HESS fit: + u pass            ", omp_get_wtime()-tf0
end subroutine hess_fit_second

module subroutine hess_fit_u_pass()
real(8),allocatable :: xb(:,:), x1b(:,:,:), x2b(:,:,:,:), wv(:), vv(:,:), vvT(:,:), WxV(:,:,:), Vw(:,:,:)
real(8),allocatable :: uy(:,:), Hloc(:,:,:), Hp(:,:,:), Hpp(:,:,:), Hdw(:,:,:), Cl(:,:), Cs(:,:)
real(8) :: sA(natoms,3,3), tot(3,3), su, w2
integer :: ib, nb, bstart, g, yy, i2, ee, ff, A, pg, B
allocate(Cs(3*n,3*n))
Cs = 0.0d0
!$omp parallel default(shared) &
!$omp&  private(xb,x1b,x2b,wv,vv,vvT,WxV,Vw,uy,Hloc,Hp,Hpp,Hdw,Cl,ib,nb,bstart,g,yy,i2,is,ee,ff,A,pg,sA,tot,su,w2)
allocate(Hloc(N3,N3,ns), Hdw(N3,N3,ns), Cl(3*n,3*n))
Hloc = 0.0d0
Hdw = 0.0d0
Cl = 0.0d0
!$omp do schedule(dynamic)
do ib = 0,nblocks-1
   bstart = cosx_blkstart_hi(ib) + 1
   nb = cosx_blkstart_hi(ib+1) - cosx_blkstart_hi(ib)
   allocate(xb(n,nb), x1b(n,3,nb), x2b(n,3,3,nb), wv(nb), vv(nb,n), vvT(n,nb), WxV(n,3,nb), Vw(n,3,nb))
   allocate(uy(nb,N3), Hp(nb,N3,3), Hpp(nb,3,3))
   call hess_block_ao(ib, nb, bstart, xb, x1b, x2b)
   wv = cosx_w_hi(bstart:bstart+nb-1)
   do g = 1,nb
      Vw(:,:,g) = 2.0d0*wv(g)*x1b(:,:,g)
   enddo
   call dgemm('N','T',3*n,3*n,nb,1.0d0,x1b,3*n,Vw,3*n,1.0d0,Cl,3*n)
   do is = 1,ns
      call dgemm('N','N',n,nb,n,1.0d0,Wbar(:,:,is),n,xb,n,0.0d0,vvT,n)
      do g = 1,nb
         uall(bstart+g-1,is) = dot_product(xb(:,g), vvT(:,g))
      enddo
      uy = 0.0d0
      do g = 1,nb
         pg = cosx_atom_hi(bstart+g-1)
         do ee = 1,3
            su = 0.0d0
            do A = 1,natoms
               tot(1,1) = dot_product(x1b(at_off(A)+1:at_off(A)+at_n(A),ee,g), vvT(at_off(A)+1:at_off(A)+at_n(A),g))
               uy(g,ee+3*(A-1)) = uy(g,ee+3*(A-1)) - 2.0d0*tot(1,1)
               su = su + tot(1,1)
            enddo
            uy(g,ee+3*(pg-1)) = uy(g,ee+3*(pg-1)) + 2.0d0*su
         enddo
      enddo
      call dgemm('N','N',N3,N3,nb,1.0d0,dwdr(1,1,bstart),N3,uy,nb,1.0d0,Hdw(1,1,is),N3)
      Hp = 0.0d0 ; Hpp = 0.0d0
      vv = 2.0d0*transpose(vvT)
      call hess_xxy_contract(nb, bstart, ib, x2b, wv, vv, Hloc(:,:,is), Hp, Hpp)
      call hess_scatter_par(nb, bstart, Hp, Hpp, Hloc(:,:,is))
      call dgemm('N','N',n,3*nb,n,1.0d0,Wbar(:,:,is),n,x1b,n,0.0d0,WxV,n)
      do g = 1,nb
         pg = cosx_atom_hi(bstart+g-1)
         w2 = 2.0d0*wv(g)
         tot = 0.0d0
         do A = 1,natoms
            do ff = 1,3
               do ee = 1,3
                  sA(A,ee,ff) = dot_product(x1b(at_off(A)+1:at_off(A)+at_n(A),ee,g), &
                                            WxV(at_off(A)+1:at_off(A)+at_n(A),ff,g))
                  tot(ee,ff) = tot(ee,ff) + sA(A,ee,ff)
               enddo
            enddo
         enddo
         do A = 1,natoms
            do ff = 1,3
               do ee = 1,3
                  Hloc(ee+3*(A-1), ff+3*(pg-1), is) = Hloc(ee+3*(A-1), ff+3*(pg-1), is) - w2*sA(A,ee,ff)
                  Hloc(ee+3*(pg-1), ff+3*(A-1), is) = Hloc(ee+3*(pg-1), ff+3*(A-1), is) - w2*sA(A,ff,ee)
               enddo
            enddo
         enddo
         do ff = 1,3
            do ee = 1,3
               Hloc(ee+3*(pg-1), ff+3*(pg-1), is) = Hloc(ee+3*(pg-1), ff+3*(pg-1), is) + w2*tot(ee,ff)
            enddo
         enddo
      enddo
   enddo
   deallocate(xb, x1b, x2b, wv, vv, vvT, WxV, Vw, uy, Hp, Hpp)
enddo
!$omp end do
block
   integer :: tid
   do tid = 0, omp_get_num_threads()-1
      !$omp barrier
      if (omp_get_thread_num() .eq. tid) then
         do is = 1,ns
            Hterm(:,:,6,is) = Hterm(:,:,6,is) - Hloc(:,:,is) - Hdw(:,:,is) - transpose(Hdw(:,:,is))
         enddo
         Cs = Cs + Cl
      endif
   enddo
end block
deallocate(Hloc, Hdw, Cl)
!$omp end parallel
do is = 1,ns
   do B = 1,natoms
      do A = 1,natoms
         do ff = 1,3
            do ee = 1,3
               su = 0.0d0
               do i2 = 1,at_n(B)
                  su = su + dot_product(Wbar(at_off(A)+1:at_off(A)+at_n(A), at_off(B)+i2, is), &
                                        Cs(n*(ee-1)+at_off(A)+1:n*(ee-1)+at_off(A)+at_n(A), n*(ff-1)+at_off(B)+i2))
               enddo
               Hterm(ee+3*(A-1), ff+3*(B-1), 6, is) = Hterm(ee+3*(A-1), ff+3*(B-1), 6, is) - su
            enddo
         enddo
      enddo
   enddo
enddo
deallocate(Cs)
end subroutine hess_fit_u_pass

end submodule exchange_cosx_hess_fit_impl
