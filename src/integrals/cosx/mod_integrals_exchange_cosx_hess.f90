! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! COSX analytic exchange Hessian at FIXED density - descendant submodule of

submodule (mod_integrals:exchange_cosx_impl) exchange_cosx_hess_impl
use mod_fixq, only: to_q, from_q
use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_wtime, omp_get_max_threads, &
                   omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, omp_destroy_lock
use MOL_info, only: natoms, engine_verbose
use GRID_info, only: cosx_becke_treutler
interface
   subroutine vv10_nlc_weight_gradient(np, coor, atom_of, nper, rho, eps, wquad, beta_v, gv_wt, dwdr_out)
   use MOL_info, only: natoms
   implicit none
   integer,intent(in) :: np, nper
   real(8),intent(in) :: coor(3,np)
   integer,intent(in) :: atom_of(np)
   real(8),intent(in) :: rho(np), eps(np), wquad(np), beta_v
   real(8),intent(out) :: gv_wt(natoms,3)
   real(8),optional,intent(out) :: dwdr_out(3,natoms,np)
   end subroutine vv10_nlc_weight_gradient
   subroutine vv10_nlc_weight_hessian_multi(np, coor, atom_of, wquad, nf, fct, hess_out)
     use MOL_info, only: natoms
     implicit none
     integer,intent(in) :: np, nf
     real(8),intent(in) :: coor(3,np)
     integer,intent(in) :: atom_of(np)
     real(8),intent(in) :: wquad(np), fct(np,nf)
     real(8),intent(out) :: hess_out(natoms,3,natoms,3,nf)
   end subroutine vv10_nlc_weight_hessian_multi
   subroutine vv10_nlc_weight_hessian(np, coor, atom_of, wquad, fct, hess_out)
   use MOL_info, only: natoms
   implicit none
   integer,intent(in) :: np
   real(8),intent(in) :: coor(3,np)
   integer,intent(in) :: atom_of(np)
   real(8),intent(in) :: wquad(np), fct(np)
   real(8),intent(out) :: hess_out(natoms,3,natoms,3)
   end subroutine vv10_nlc_weight_hessian
end interface

integer,parameter :: NT = 7
integer,parameter :: INT_CACHE_SIZE = 4000000
integer :: N3, n, npts, nblocks, ns, is, maxd, info
logical :: fit, same_spin, verbose
real(8),pointer,contiguous :: Pf(:,:)
real(8),allocatable :: Dsp(:,:,:), DPt(:,:,:)
real(8),allocatable :: dwdr(:,:,:), epsall(:,:), uall(:,:)
real(8),allocatable :: Hterm(:,:,:,:), grad1(:,:)
real(8),allocatable :: Sy(:,:,:), Py(:,:,:)
real(8),allocatable :: Ym(:,:,:), Wbar(:,:,:)
integer,allocatable :: ao_atom(:), at_off(:), at_n(:)
real(8) :: kscreen_h, kscreen_h2, kscreen_hk
character(len=32) :: envs
integer :: y1, y2, iat, d, e, mu
real(8),allocatable :: zero_np(:), gvdum(:,:), hw(:,:,:,:)
real(8) :: t0, t_setup, t_main, t_fit
logical :: want_kx, kproj
real(8),allocatable :: Kx(:,:,:,:), Pt(:,:), Cp(:,:,:)
integer :: nop, noc(2)
real(8),allocatable :: PtC(:,:,:,:), XWg(:,:,:), PtCp(:,:,:)
real(8) :: pymax, pcol
real(8),allocatable :: cmx(:)
logical :: lowrank
integer(omp_lock_kind),allocatable :: zlock(:)
integer(16),allocatable :: XKsum(:,:,:,:)
integer(omp_lock_kind),allocatable :: xklk(:)
real(8),allocatable :: cnrm(:,:)
integer,allocatable :: ao_sh(:)
integer :: nab
real(8),allocatable :: dS1c(:,:,:)
logical :: lrfit
real(8),allocatable :: Zt(:,:,:,:)
integer(16),allocatable :: Ztq(:,:,:,:)
real(8) :: hs_rsomega

interface
module subroutine hmem(tag)
implicit none
character(len=*),intent(in) :: tag
end subroutine hmem

module subroutine hess_block_ao(ib, nb, bstart, xb, x1b, x2b)
implicit none
integer,intent(in) :: ib, nb, bstart
real(8),intent(out) :: xb(n,nb), x1b(n,3,nb)
real(8),intent(out),optional :: x2b(n,3,3,nb)
end subroutine hess_block_ao

module subroutine hess_block_atoms(ib, nloc, loc)
implicit none
integer,intent(in) :: ib
integer,intent(out) :: nloc, loc(natoms)
end subroutine hess_block_atoms

module subroutine hess_tangent_occ(nb, bstart, nloc, loc, x1b, L, T)
implicit none
integer,intent(in) :: nb, bstart, nloc, loc(natoms)
real(8),intent(in) :: x1b(n,3,nb), L(n,nop)
real(8),intent(inout) :: T(nb,nop,N3)
end subroutine hess_tangent_occ

module subroutine hess_xxy_contract(nb, bstart, ib, x2b, wv, Rv, Hloc, Hp, Hpp)
implicit none
integer,intent(in) :: nb, bstart, ib
real(8),intent(in) :: x2b(n,3,3,nb), wv(nb), Rv(nb,n)
real(8),intent(inout) :: Hloc(N3,N3), Hp(nb,N3,3), Hpp(nb,3,3)
end subroutine hess_xxy_contract

module subroutine hess_scatter_par(nb, bstart, Hp, Hpp, Hloc)
implicit none
integer,intent(in) :: nb, bstart
real(8),intent(in) :: Hp(nb,N3,3), Hpp(nb,3,3)
real(8),intent(inout) :: Hloc(N3,N3)
end subroutine hess_scatter_par

module subroutine hess_rows_shared_c(nb, bstart, is, x1b, Rm, Z1)
implicit none
integer,intent(in) :: nb, bstart, is
real(8),intent(in) :: x1b(n,3,nb), Rm(nb,nop)
real(8),intent(inout) :: Z1(n,nop,3)
end subroutine hess_rows_shared_c

module subroutine kx_free_leg_max(nb, bstart, xb, x1b, wv, qv, Lpt)
implicit none
integer,intent(in) :: nb, bstart
real(8),intent(in) :: xb(n,nb), x1b(n,3,nb), wv(nb), qv(nb,n)
real(8),intent(out) :: Lpt(nb)
end subroutine kx_free_leg_max

module subroutine kx_proj_leg(nb, bstart, is, xb, wv, x1b, nloc, loc, Gl, KCa, XWl, tk)
implicit none
integer,intent(in) :: nb, bstart, is, nloc, loc(natoms)
real(8),intent(in) :: xb(n,nb), wv(nb), x1b(n,3,nb), Gl(nb,n), KCa(nb,N3,nop)
real(8),intent(inout) :: XWl(n,nop,*), tk(6)
end subroutine kx_proj_leg

module subroutine hess_post_sub(ng, bst, ib, is, nloc, loc, xs, x1s, x2s, ws, avs, Gvs, Gtvs, KCs_, MCs_, &
                         XAs, A1s, B1s, Hl, gl, Zl, XWl, tk, tsx)
implicit none
integer,intent(in) :: ng, bst, ib, is, nloc, loc(natoms)
real(8),intent(in) :: xs(n,ng), x1s(n,3,ng), x2s(n,3,3,ng), ws(ng), avs(ng,n), Gvs(ng,n), Gtvs(ng,n)
real(8),intent(in) :: KCs_(ng,N3,nop), MCs_(ng,N3,nop), XAs(ng,nop), A1s(ng,nop,N3), B1s(ng,nop,N3)
real(8),intent(inout) :: Hl(N3,N3,NT), gl(N3), Zl(n,nop,3), XWl(n,nop,*)
real(8),intent(inout) :: tk(6), tsx(10)
end subroutine hess_post_sub

module subroutine a1_shell_max(nb, is, sw, A1, Tu, Tmax)
implicit none
integer,intent(in) :: nb, is
real(8),intent(in) :: sw(nb), A1(nb,nop,N3)
real(8),intent(inout) :: Tu(0:nBases-1), Tmax(0:nBases-1)
end subroutine a1_shell_max

module subroutine hess_sy_times(M, SM)
implicit none
real(8),intent(in) :: M(n,nop)
real(8),intent(inout) :: SM(n,nop,N3)
end subroutine hess_sy_times

module subroutine hess_snum_first(Sout)
implicit none
real(8),intent(out) :: Sout(n,n,N3)
end subroutine hess_snum_first

module subroutine hess_fit_py(Sin, Pout)
implicit none
real(8),intent(in) :: Sin(n,n,N3)
real(8),intent(out) :: Pout(n,n,N3)
end subroutine hess_fit_py

module subroutine hess_snum_rng(a_lo, a_hi, Sout)
implicit none
integer,intent(in) :: a_lo, a_hi
real(8),intent(out) :: Sout(n,n,3*(a_hi-a_lo+1))
end subroutine hess_snum_rng

module subroutine hess_rows_rng(nb, bstart, nloc, loc, x1b, Rm, Zout, a_lo, a_hi, slk)
implicit none
integer,intent(in) :: nb, bstart, nloc, loc(natoms), a_lo, a_hi
real(8),intent(in) :: x1b(n,3,nb), Rm(nb,n)
integer(16),intent(inout) :: Zout(n,n,3*(a_hi-a_lo+1))
integer(omp_lock_kind),intent(inout) :: slk(a_lo:a_hi)
end subroutine hess_rows_rng

module subroutine hess_fit_py_rng(Sin, Pout, a_lo, a_hi)
implicit none
integer,intent(in) :: a_lo, a_hi
real(8),intent(in) :: Sin(n,n,3*(a_hi-a_lo+1))
real(8),intent(out) :: Pout(n,n,3*(a_hi-a_lo+1))
end subroutine hess_fit_py_rng

module subroutine hess_main_pass()
implicit none
end subroutine hess_main_pass

module subroutine hess_fit_second()
implicit none
end subroutine hess_fit_second

module subroutine hess_fit_u_pass()
implicit none
end subroutine hess_fit_u_pass

end interface

contains

module subroutine cosx_hess_exchange(nConts, natoms, Da, Db, RSomega, d2Ea, d2Eb, Kxa, Kxb, Cocca, Coccb)
use omp_lib, only: omp_get_thread_num, omp_get_num_threads, omp_get_wtime, omp_get_max_threads, &
                   omp_lock_kind, omp_init_lock, omp_set_lock, omp_unset_lock, omp_destroy_lock
use MOL_info, only: engine_verbose
use GRID_info, only: cosx_becke_treutler
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
real(8),intent(in) :: RSomega
real(8),intent(out) :: d2Ea(3,natoms,3,natoms), d2Eb(3,natoms,3,natoms)
real(8),intent(out),optional :: Kxa(nConts,nConts,3,natoms), Kxb(nConts,nConts,3,natoms)
real(8),intent(in),optional :: Cocca(:,:), Coccb(:,:)
hs_rsomega = RSomega
if (allocated(Dsp)) deallocate(Dsp)
if (allocated(DPt)) deallocate(DPt)
if (allocated(dwdr)) deallocate(dwdr)
if (allocated(epsall)) deallocate(epsall)
if (allocated(uall)) deallocate(uall)
if (allocated(Hterm)) deallocate(Hterm)
if (allocated(grad1)) deallocate(grad1)
if (allocated(Sy)) deallocate(Sy)
if (allocated(Py)) deallocate(Py)
if (allocated(Ym)) deallocate(Ym)
if (allocated(Wbar)) deallocate(Wbar)
if (allocated(ao_atom)) deallocate(ao_atom)
if (allocated(at_off)) deallocate(at_off)
if (allocated(at_n)) deallocate(at_n)
if (allocated(zero_np)) deallocate(zero_np)
if (allocated(gvdum)) deallocate(gvdum)
if (allocated(hw)) deallocate(hw)
if (allocated(Kx)) deallocate(Kx)
if (allocated(Pt)) deallocate(Pt)
if (allocated(Cp)) deallocate(Cp)
if (allocated(PtC)) deallocate(PtC)
if (allocated(XWg)) deallocate(XWg)
if (allocated(PtCp)) deallocate(PtCp)
if (allocated(cmx)) deallocate(cmx)
if (allocated(zlock)) deallocate(zlock)
if (allocated(XKsum)) deallocate(XKsum)
if (allocated(xklk)) deallocate(xklk)
if (allocated(cnrm)) deallocate(cnrm)
if (allocated(ao_sh)) deallocate(ao_sh)
if (allocated(dS1c)) deallocate(dS1c)
if (allocated(Zt)) deallocate(Zt)
if (allocated(Ztq)) deallocate(Ztq)

cosx_becke_treutler = cosx_hi_treutler
N3 = 3*natoms
n = nConts
npts = cosx_npts_hi
nblocks = cosx_nblk_hi
same_spin = all(Da .eq. Db)
ns = merge(1, 2, same_spin)
call get_environment_variable("ENGINE_COSX_HESS_VERBOSE", envs)
verbose = (trim(envs) .eq. "1")
kscreen_h = 1.0d-10
call get_environment_variable("ENGINE_COSX_HESS_KSCREEN", envs)
if (len_trim(envs) .gt. 0) read(envs,*) kscreen_h
kscreen_h2 = 10.0d0*kscreen_h
call get_environment_variable("ENGINE_COSX_HESS_KSCREEN2", envs)
if (len_trim(envs) .gt. 0) read(envs,*) kscreen_h2
kscreen_hk = 1.0d-8
call get_environment_variable("ENGINE_COSX_HESS_KSCREENK", envs)
if (len_trim(envs) .gt. 0) read(envs,*) kscreen_hk

fit = cosx_overlap_fit .and. allocated(cosx_QfullT_hi) .and. cosx_Snum_lu_ready
if (fit) fit = (size(cosx_QfullT_hi,1) .eq. n)
allocate(Dsp(n,n,ns), DPt(n,n,ns))
Dsp(:,:,1) = Da
if (ns .eq. 2) Dsp(:,:,2) = Db
if (fit) then
   Pf => cosx_QfullT_hi
   do is = 1,ns
      call dgemm('T','N',n,n,n,1.0d0,Pf,n,Dsp(:,:,is),n,0.0d0,DPt(:,:,is),n)
   enddo
else
   DPt = Dsp
endif

maxd = 0
do iat = 0,nBases-1
   maxd = max(maxd, cosx_shell_dim(iat))
enddo
allocate(ao_atom(n), at_off(natoms), at_n(natoms))
at_n = 0
at_off = n
do iat = 0,nBases-1
   do mu = 1,cosx_shell_dim(iat)
      ao_atom(cosx_ao_offset(iat)+mu) = cosx_shell_atom(iat)
   enddo
enddo
do mu = 1,n
   at_n(ao_atom(mu)) = at_n(ao_atom(mu)) + 1
   at_off(ao_atom(mu)) = min(at_off(ao_atom(mu)), mu-1)
enddo
do mu = 1,n
   if (mu-1 .lt. at_off(ao_atom(mu)) .or. mu-1 .ge. at_off(ao_atom(mu))+at_n(ao_atom(mu))) then
      print *, "cosx_hess_exchange: AO order is not atom-major"
      stop
   endif
enddo

t0 = omp_get_wtime()
allocate(dwdr(3,natoms,npts), zero_np(npts), gvdum(natoms,3))
zero_np = 0.0d0
call vv10_nlc_weight_gradient(npts, cosx_coor_hi, cosx_atom_hi, max(npts/max(natoms,1),1), &
     zero_np, zero_np, cosx_w_hi, 0.0d0, gvdum, dwdr)
if (verbose) print '(A,F8.3)', " COSX HESS setup: weight gradient ", omp_get_wtime()-t0
call hmem("start  ")
if (verbose) then
   block
      integer :: ib2, yq, cnt(4), tot2
      real(8) :: bm, th4(4)
      th4 = (/1.0d-8, 1.0d-10, 1.0d-12, 1.0d-14/)
      cnt = 0 ; tot2 = 0
      do ib2 = 0, cosx_nblk_hi-1
         do yq = 1,N3
            bm = maxval(abs(dwdr(mod(yq-1,3)+1,(yq-1)/3+1,cosx_blkstart_hi(ib2)+1:cosx_blkstart_hi(ib2+1))))
            tot2 = tot2 + 1
            do d = 1,4
               if (bm .ge. th4(d)) cnt(d) = cnt(d) + 1
            enddo
         enddo
      enddo
      print '(A,4F7.3)', " COSX HESS (block,y) with max|dw| >= 1e-8/-10/-12/-14 (fraction):", dble(cnt)/dble(tot2)
   end block
endif

allocate(Hterm(N3,N3,NT,ns), grad1(N3,ns), epsall(npts,ns))
Hterm = 0.0d0
grad1 = 0.0d0
epsall = 0.0d0

t_setup = omp_get_wtime() - t0

want_kx = present(Kxa)
allocate(Pt(n,n))
if (fit) then
   Pt = transpose(Pf)
else
   Pt = 0.0d0
   do mu = 1,n
      Pt(mu,mu) = 1.0d0
   enddo
endif
kproj = want_kx .and. present(Cocca)
if (kproj) then
   noc(1) = size(Cocca,2)
   noc(2) = noc(1)
   if (present(Coccb)) noc(2) = size(Coccb,2)
   nop = max(noc(1), noc(2))
   allocate(Cp(n,nop,ns))
   Cp = 0.0d0
   Cp(:,1:noc(1),1) = Cocca
   if (ns .eq. 2) Cp(:,1:noc(2),2) = Coccb
else
   nop = n
endif
if (want_kx) then
   allocate(Kx(n,nop,N3,ns))
   Kx = 0.0d0
endif
allocate(cmx(0:nBases-1))
cmx = 1.0d0
if (kproj) then
   do iat = 0,nBases-1
      cmx(iat) = maxval(abs(Cp(cosx_ao_offset(iat)+1:cosx_ao_offset(iat)+cosx_shell_dim(iat),:,:)))
   enddo
endif
pymax = 0.0d0
pcol = 1.0d0
do mu = 1,n
   pcol = max(pcol, sum(abs(Pt(:,mu))))
enddo
if (kproj) then
   allocate(PtCp(n,nop,ns))
   do is = 1,ns
      call dgemm('N','N',n,nop,n,1.0d0,Pt,n,Cp(1,1,is),n,0.0d0,PtCp(1,1,is),n)
   enddo
endif
lowrank = .false.
if (kproj) then
   block
      real(8),allocatable :: dcc(:,:)
      real(8) :: derr
      allocate(dcc(n,n))
      derr = 0.0d0
      do is = 1,ns
         dcc = Dsp(:,:,is)
         call dgemm('N','T',n,n,nop,-1.0d0,Cp(1,1,is),n,Cp(1,1,is),n,1.0d0,dcc,n)
         derr = max(derr, maxval(abs(dcc)))
      enddo
      lowrank = (derr .lt. 1.0d-10)
      if (verbose) print '(A,ES10.3,A,L2)', " COSX HESS max|D - C C^T| =", derr, "  low-rank fit tangents:", lowrank
      deallocate(dcc)
   end block
endif
t0 = omp_get_wtime()
lrfit = .false.
if (fit) then
   call get_environment_variable("ENGINE_COSX_HESS_LRFIT", envs)
   lrfit = (trim(envs) .eq. "1")
endif
if (fit .and. lrfit) then
   block
      real(8),allocatable :: Sb(:,:,:), Pb(:,:,:), dH1e(:,:,:)
      real(8) :: symem
      integer :: a_lo, a_hi, nyb, y0, yl
      symem = 512.0d6
      call get_environment_variable("ENGINE_COSX_HESS_SYMEM", envs)
      if (len_trim(envs) .gt. 0) read(envs,*) symem
      nab = int(symem/(120.0d0*dble(n)**2))
      nab = max(1, min(natoms, nab))
      nab = (natoms + (natoms+nab-1)/nab - 1)/((natoms+nab-1)/nab)
      allocate(dS1c(n,n,3), dH1e(n,n,3))
      call integrals_compute_force_1e(n, dS1c, dH1e)
      deallocate(dH1e)
      allocate(PtC(n,nop,N3,ns), XWg(n,nop,ns))
      XWg = 0.0d0
      do a_lo = 1,natoms,nab
         a_hi = min(natoms, a_lo+nab-1)
         nyb = 3*(a_hi-a_lo+1)
         y0 = 3*(a_lo-1)
         allocate(Sb(n,n,nyb), Pb(n,n,nyb))
         call hess_snum_rng(a_lo, a_hi, Sb)
         call hess_fit_py_rng(Sb, Pb, a_lo, a_hi)
         do is = 1,ns
            do yl = 1,nyb
               call dgemm('T','N',n,nop,n,1.0d0,Pb(:,:,yl),n,Cp(1,1,is),n,0.0d0,PtC(1,1,y0+yl,is),n)
            enddo
         enddo
         do yl = 1,nyb
            do mu = 1,n
               pymax = max(pymax, sum(abs(Pb(mu,:,yl))))
            enddo
         enddo
         deallocate(Sb, Pb)
      enddo
      if (verbose) print '(A,F8.3,A,I0,A)', " COSX HESS setup: S^y, P^y -> PtC   ", omp_get_wtime()-t0, &
            "   (", nab, " atoms per range)"
   end block
else if (fit) then
   allocate(Sy(n,n,N3), Py(n,n,N3))
   call hess_snum_first(Sy)
   if (verbose) print '(A,F8.3)', " COSX HESS setup: + S^y numeric    ", omp_get_wtime()-t0
   call hmem("S^y    ")
   call hess_fit_py(Sy, Py)
   if (verbose) print '(A,F8.3)', " COSX HESS setup: + P^y            ", omp_get_wtime()-t0
   call hmem("P^y    ")
   allocate(PtC(n,nop,N3,ns), XWg(n,nop,ns))
   XWg = 0.0d0
   do is = 1,ns
      do y1 = 1,N3
         call dgemm('T','N',n,nop,n,1.0d0,Py(:,:,y1),n,Cp(1,1,is),n,0.0d0,PtC(1,1,y1,is),n)
      enddo
   enddo
   do y1 = 1,N3
      do mu = 1,n
         pymax = max(pymax, sum(abs(Py(mu,:,y1))))
      enddo
   enddo
endif
if (fit) then
   allocate(Zt(n,nop,N3,ns), Ztq(n,nop,N3,ns))
   Zt = 0.0d0
   Ztq = 0_16
endif
t_setup = t_setup + (omp_get_wtime() - t0)
t0 = omp_get_wtime()
call hess_main_pass()
if (allocated(Ztq)) then
   Zt = Zt + from_q(Ztq)
   deallocate(Ztq)
endif
call hmem("main   ")
t_main = omp_get_wtime() - t0
t0 = omp_get_wtime()
if (fit) then
   allocate(Ym(n,n,ns), Wbar(n,n,ns), uall(npts,ns))
   if (lrfit) then
      do is = 1,ns
         call dgetrs('N', n, nop, cosx_Snum_lu_hi, n, cosx_ipiv_hi, XWg(1,1,is), n, info)
      enddo
   endif
   do is = 1,ns
      if (lrfit) then
         call dgemm('N','T',n,n,nop,1.0d0,XWg(1,1,is),n,Cp(1,1,is),n,0.0d0,Ym(1,1,is),n)
      else
         block
            real(8),allocatable :: Rt(:,:)
            allocate(Rt(n,nop))
            Rt = XWg(:,:,is)
            call dgetrs('N', n, nop, cosx_Snum_lu_hi, n, cosx_ipiv_hi, Rt, n, info)
            call dgemm('N','T',n,n,nop,1.0d0,Rt,n,Cp(1,1,is),n,0.0d0,Ym(1,1,is),n)
            deallocate(Rt)
         end block
      endif
      call dgemm('N','N',n,n,n,1.0d0,Ym(:,:,is),n,Pf,n,0.0d0,Wbar(:,:,is),n)
      Wbar(:,:,is) = 0.5d0*(Wbar(:,:,is) + transpose(Wbar(:,:,is)))
   enddo
   call hess_fit_second()
   deallocate(Ym, Wbar)
endif
t_fit = omp_get_wtime() - t0
t0 = omp_get_wtime()
if (kproj .and. fit) then
   if (lrfit) then
      deallocate(PtC)
   else
      do is = 1,ns
         do y1 = 1,N3
            call dgemm('N','N',n,nop,n,1.0d0,Py(:,:,y1),n,XWg(1,1,is),n,1.0d0,Kx(1,1,y1,is),n)
         enddo
      enddo
      deallocate(PtC)
   endif
   deallocate(XWg)
endif
if (kproj) then
   do iat = 1,natoms
      do d = 1,3
         Kxa(:,1:noc(1),d,iat) = Kx(:,1:noc(1),d+3*(iat-1),1)
         if (present(Kxb)) Kxb(:,1:noc(2),d,iat) = Kx(:,1:noc(2),d+3*(iat-1),ns)
      enddo
   enddo
   deallocate(Kx, Cp, PtCp)
endif
deallocate(Pt)
t_main = t_main + (omp_get_wtime() - t0)

allocate(hw(natoms,3,natoms,3))
block
   real(8),allocatable :: fm(:,:), hw2(:,:,:,:,:)
   integer :: nfm
   nfm = merge(2, 1, fit)
   allocate(fm(npts,nfm), hw2(natoms,3,natoms,3,nfm))
   do is = 1,ns
      fm(:,1) = epsall(:,is)
      if (fit) fm(:,2) = uall(:,is)
      call vv10_nlc_weight_hessian_multi(npts, cosx_coor_hi, cosx_atom_hi, cosx_w_hi, nfm, fm, hw2)
      do iat = 1,natoms
         do d = 1,3
            do e = 1,3
               do mu = 1,natoms
                  Hterm(d+3*(iat-1), e+3*(mu-1), 1, is) = hw2(iat,d,mu,e,1)
                  if (fit) Hterm(d+3*(iat-1), e+3*(mu-1), 6, is) = Hterm(d+3*(iat-1), e+3*(mu-1), 6, is) &
                                                                   - hw2(iat,d,mu,e,2)
               enddo
            enddo
         enddo
      enddo
   enddo
   deallocate(fm, hw2)
end block
if (fit) deallocate(uall)
t_fit = t_fit + (omp_get_wtime() - t0)
if (verbose) print '(A,F8.3)', " COSX HESS fit + weight hessians (T1, T6b):", t_fit
call hmem("fit    ")

do is = 1,ns
   Hterm(:,:,2,is) = Hterm(:,:,2,is) + transpose(Hterm(:,:,2,is))
   Hterm(:,:,3,is) = Hterm(:,:,3,is) + transpose(Hterm(:,:,3,is))
   Hterm(:,:,7,is) = Hterm(:,:,7,is) + transpose(Hterm(:,:,7,is))
enddo

block
   real(8) :: Htot(N3,N3)
   do is = 1,ns
      Htot = sum(Hterm(:,:,:,is), dim=3)
      if (is .eq. 1) then
         d2Ea = reshape(-0.25d0*Htot, (/3,natoms,3,natoms/))
      else
         d2Eb = reshape(-0.25d0*Htot, (/3,natoms,3,natoms/))
      endif
   enddo
end block
if (ns .eq. 1) d2Eb = d2Ea

if (verbose) then
   print '(A,3F9.3)', " COSX HESS timing setup/main/fit (s):", t_setup, t_main, t_fit
   do is = 1,ns
      print '(A,I2,A,ES24.16)', " COSX HESS spin",is,"  E = sum w eps =", sum(cosx_w_hi(1:npts)*epsall(:,is))
      print '(A,I2,A)', " COSX HESS spin",is,"  first derivative -1/4 dE/dR (from tangents):"
      do iat = 1,natoms
         print '(3ES24.16)', -0.25d0*grad1(3*(iat-1)+1:3*iat, is)
      enddo
      do d = 1,NT
         print '(A,I2,A,ES12.4,A,ES12.4)', "   term T",d,"  max|.| =", maxval(abs(Hterm(:,:,d,is))), &
              "  transl =", maxval(abs(sum(reshape(Hterm(:,:,d,is),(/3,natoms,N3/)),dim=2)))
      enddo
   enddo
   call flush(6)
endif
call get_environment_variable("ENGINE_COSX_HESS_DUMP", envs)
if (trim(envs) .eq. "1") then
   block
      integer :: u
      open(newunit=u, file='cosx_hess_terms.txt', status='replace', action='write')
      write(u,'(3I8)') natoms, NT, ns
      do is = 1,ns
         do d = 1,NT
            do y2 = 1,N3
               write(u,'(*(ES25.17))') Hterm(:,y2,d,is)
            enddo
         enddo
      enddo
      close(u)
   end block
endif

deallocate(Dsp, DPt, dwdr, zero_np, gvdum, Hterm, grad1, epsall, hw, ao_atom, at_off, at_n, cmx)
if (allocated(Sy)) deallocate(Sy, Py)
if (allocated(dS1c)) deallocate(dS1c)
if (allocated(Zt)) deallocate(Zt)
if (allocated(Ztq)) deallocate(Ztq)
cosx_becke_treutler = .false.
call cosx_malloc_trim()

if (allocated(Dsp)) deallocate(Dsp)
if (allocated(DPt)) deallocate(DPt)
if (allocated(dwdr)) deallocate(dwdr)
if (allocated(epsall)) deallocate(epsall)
if (allocated(uall)) deallocate(uall)
if (allocated(Hterm)) deallocate(Hterm)
if (allocated(grad1)) deallocate(grad1)
if (allocated(Sy)) deallocate(Sy)
if (allocated(Py)) deallocate(Py)
if (allocated(Ym)) deallocate(Ym)
if (allocated(Wbar)) deallocate(Wbar)
if (allocated(ao_atom)) deallocate(ao_atom)
if (allocated(at_off)) deallocate(at_off)
if (allocated(at_n)) deallocate(at_n)
if (allocated(zero_np)) deallocate(zero_np)
if (allocated(gvdum)) deallocate(gvdum)
if (allocated(hw)) deallocate(hw)
if (allocated(Kx)) deallocate(Kx)
if (allocated(Pt)) deallocate(Pt)
if (allocated(Cp)) deallocate(Cp)
if (allocated(PtC)) deallocate(PtC)
if (allocated(XWg)) deallocate(XWg)
if (allocated(PtCp)) deallocate(PtCp)
if (allocated(cmx)) deallocate(cmx)
if (allocated(zlock)) deallocate(zlock)
if (allocated(XKsum)) deallocate(XKsum)
if (allocated(xklk)) deallocate(xklk)
if (allocated(cnrm)) deallocate(cnrm)
if (allocated(ao_sh)) deallocate(ao_sh)
if (allocated(dS1c)) deallocate(dS1c)
if (allocated(Zt)) deallocate(Zt)
if (allocated(Ztq)) deallocate(Ztq)
end subroutine cosx_hess_exchange

end submodule exchange_cosx_hess_impl
