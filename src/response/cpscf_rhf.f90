! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: the closed-shell (RHF/RKS) response solvers - the PCG and

submodule (mod_cpscf) cpscf_rhf_impl
implicit none
contains

subroutine cphf_solve_vo_pcg(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                             RHSfix, Uvo_all, niter_total, conv_thr)
use mod_xc_hess, only: xc_fxc_project_vo
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc, nvir, npert
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac, conv_thr
real(8),intent(in) :: RHSfix(nvir,nocc,npert)
real(8),intent(out) :: Uvo_all(nvir,nocc,npert)
integer,intent(out) :: niter_total
real(8),allocatable :: dlt(:,:), R(:,:,:), Z(:,:,:), Pd(:,:,:), Q(:,:,:), rz(:), Qoo(:,:,:)
logical,allocatable :: conv(:)
integer,allocatable :: act(:)
integer :: p, it, nact, k, a, i
real(8) :: alpha, beta, rznew, pq

allocate(dlt(nvir,nocc), R(nvir,nocc,npert), Z(nvir,nocc,npert), Pd(nvir,nocc,npert))
allocate(Q(nvir,nocc,npert), rz(npert), conv(npert), act(npert))
if (cphf_voo_track) allocate(Qoo(nocc,nocc,npert))
do i = 1,nocc
   do a = 1,nvir
      dlt(a,i) = eLev_a(nocc+a) - eLev_a(i)
      if (cphf_omega2 /= 0.0d0) dlt(a,i) = dlt(a,i) - cphf_omega2/dlt(a,i)
   enddo
enddo
block
   character(len=16) :: solv
   call get_environment_variable("ENGINE_CPHF_SOLVER", solv)
   if (trim(solv) .ne. "pcg") then
      call bcg_rhf()
      deallocate(dlt, R, Z, Pd, Q, rz, conv, act)
      if (cphf_voo_track) then
         deallocate(Qoo)
         cphf_voo_ready = .true.
      endif
      return
   endif
end block
do p = 1,npert
   Uvo_all(:,:,p) = -RHSfix(:,:,p)/dlt
   act(p) = p
enddo
call apply_A(npert, act, Uvo_all, Q)
if (cphf_voo_track) cphf_voo_u = Qoo
niter_total = npert
do p = 1,npert
   R(:,:,p) = -RHSfix(:,:,p) - Q(:,:,p)
   Z(:,:,p) = R(:,:,p)/dlt
   Pd(:,:,p) = Z(:,:,p)
   rz(p) = sum(R(:,:,p)*Z(:,:,p))
   conv(p) = maxval(abs(Z(:,:,p))) < conv_thr
enddo
do it = 1,200
   nact = 0
   do p = 1,npert
      if (.not. conv(p)) then
         nact = nact + 1
         act(nact) = p
      endif
   enddo
   if (nact .eq. 0) exit
   niter_total = niter_total + nact
   call apply_A(nact, act, Pd, Q)
   do k = 1,nact
      p = act(k)
      pq = sum(Pd(:,:,p)*Q(:,:,p))
      if (cphf_omega2 /= 0.0d0 .and. pq <= 0.0d0) then
         print '(A,F10.6,A)', 'ERROR: dynamic polarizability: omega = ', sqrt(cphf_omega2), &
               ' a.u. is not below the lowest excitation (negative curvature)'
         error stop 'dynamic polarizability: omega beyond the first excitation (resonance region)'
      endif
      alpha = rz(p)/pq
      Uvo_all(:,:,p) = Uvo_all(:,:,p) + alpha*Pd(:,:,p)
      if (cphf_voo_track) cphf_voo_u(:,:,p) = cphf_voo_u(:,:,p) + alpha*Qoo(:,:,p)
      R(:,:,p) = R(:,:,p) - alpha*Q(:,:,p)
      Z(:,:,p) = R(:,:,p)/dlt
      rznew = sum(R(:,:,p)*Z(:,:,p))
      beta = rznew/rz(p)
      rz(p) = rznew
      Pd(:,:,p) = Z(:,:,p) + beta*Pd(:,:,p)
      conv(p) = maxval(abs(Z(:,:,p))) < conv_thr
   enddo
enddo
deallocate(dlt, R, Z, Pd, Q, rz, conv, act)
if (cphf_voo_track) then
   deallocate(Qoo)
   cphf_voo_ready = .true.
endif

contains

subroutine bcg_rhf()
integer :: nt_, it_, k_, p_, q_, nk_, info_, lwork_
real(8),allocatable :: G_(:,:), ev_(:), work_(:), T_(:,:), W_(:,:), Wo(:,:)
logical :: verb_
character(len=8) :: e_
verb_ = (engine_verbose .ge. 2)
nt_ = nvir*nocc
do p_ = 1,npert
   Uvo_all(:,:,p_) = -RHSfix(:,:,p_)/dlt
   act(p_) = p_
enddo
call apply_A(npert, act, Uvo_all, Q)
if (cphf_voo_track) cphf_voo_u = Qoo
niter_total = npert
R = -RHSfix - Q
nk_ = npert
do it_ = 1,200
   do p_ = 1,npert
      Z(:,:,p_) = R(:,:,p_)/dlt
      conv(p_) = maxval(abs(Z(:,:,p_))) < conv_thr
      if (conv(p_)) Z(:,:,p_) = 0.0d0
   enddo
   if (verb_) print '(A,I4,A,I5)', '  RHF BCG it', it_, '  unconverged', count(.not. conv)
   if (all(conv)) exit
   if (it_ .eq. 1) then
      Pd = Z
   else
      allocate(T_(nk_,npert), W_(nt_,npert))
      call dgemm('T','N',nk_,npert,nt_,1.0d0,Q,nt_,Z,nt_,0.0d0,T_,nk_)
      W_ = reshape(Z,(/nt_,npert/))
      call dgemm('N','N',nt_,npert,nk_,-1.0d0,Pd,nt_,T_,nk_,1.0d0,W_,nt_)
      Pd = reshape(W_,(/nvir,nocc,npert/))
      deallocate(T_, W_)
   endif
   q_ = 0
   do p_ = 1,npert
      if (.not. conv(p_)) then
         q_ = q_ + 1
         if (q_ .ne. p_) Pd(:,:,q_) = Pd(:,:,p_)
      endif
   enddo
   nk_ = q_
   do k_ = 1,nk_
      act(k_) = k_
   enddo
   call apply_A(nk_, act, Pd, Q)
   niter_total = niter_total + nk_
   allocate(G_(nk_,nk_), ev_(nk_), work_(1))
   call dgemm('T','N',nk_,nk_,nt_,1.0d0,Pd,nt_,Q,nt_,0.0d0,G_,nk_)
   G_ = 0.5d0*(G_ + transpose(G_))
   call dsyev('V','U',nk_,G_,nk_,ev_,work_,-1,info_)
   lwork_ = int(work_(1)) ; deallocate(work_) ; allocate(work_(lwork_))
   call dsyev('V','U',nk_,G_,nk_,ev_,work_,lwork_,info_)
   q_ = 0
   do k_ = 1,nk_
      if (ev_(k_) .gt. 1.0d-12*ev_(nk_)) then
         q_ = q_ + 1
         G_(:,q_) = G_(:,k_)/sqrt(ev_(k_))
      endif
   enddo
   allocate(W_(nt_,q_))
   call dgemm('N','N',nt_,q_,nk_,1.0d0,Pd,nt_,G_,nk_,0.0d0,W_,nt_)
   Pd(:,:,1:q_) = reshape(W_,(/nvir,nocc,q_/))
   call dgemm('N','N',nt_,q_,nk_,1.0d0,Q,nt_,G_,nk_,0.0d0,W_,nt_)
   Q(:,:,1:q_) = reshape(W_,(/nvir,nocc,q_/))
   deallocate(W_)
   if (cphf_voo_track) then
      allocate(Wo(nocc*nocc,q_))
      call dgemm('N','N',nocc*nocc,q_,nk_,1.0d0,Qoo,nocc*nocc,G_,nk_,0.0d0,Wo,nocc*nocc)
      Qoo(:,:,1:q_) = reshape(Wo,(/nocc,nocc,q_/))
      deallocate(Wo)
   endif
   deallocate(G_, ev_, work_)
   nk_ = q_
   allocate(T_(nk_,npert))
   call dgemm('T','N',nk_,npert,nt_,1.0d0,Pd,nt_,R,nt_,0.0d0,T_,nk_)
   call dgemm('N','N',nt_,npert,nk_,1.0d0,Pd,nt_,T_,nk_,1.0d0,Uvo_all,nt_)
   call dgemm('N','N',nt_,npert,nk_,-1.0d0,Q,nt_,T_,nk_,1.0d0,R,nt_)
   if (cphf_voo_track) call dgemm('N','N',nocc*nocc,npert,nk_,1.0d0,Qoo,nocc*nocc,T_,nk_,1.0d0,cphf_voo_u,nocc*nocc)
   deallocate(T_)
enddo
end subroutine bcg_rhf

subroutine apply_A(n_act, alist, X, Qout)
integer,intent(in) :: n_act, alist(n_act)
real(8),intent(in) :: X(nvir,nocc,npert)
real(8),intent(inout) :: Qout(nvir,nocc,npert)
real(8),allocatable :: Dtr(:,:,:), Jtr(:,:,:), Katr(:,:,:), Xp(:,:,:), Ua(:,:,:), Tmp(:,:), T2(:,:), CU(:,:)
real(8),allocatable :: VC(:,:,:), LU(:,:,:)
logical :: vcm
integer :: c, pp
real(8) :: tv0
vcm = cpscf_use_vind_c(HF_exchange_frac)
if (vcm) then
   allocate(Dtr(1,1,1), Jtr(1,1,1), Katr(1,1,1), LU(nConts,nocc,n_act))
else
   allocate(Dtr(nConts,nConts,n_act), Jtr(nConts,nConts,n_act), Katr(nConts,nConts,n_act))
endif
allocate(Tmp(nConts,nocc), T2(nConts,nConts), CU(nConts,nocc))
do c = 1,n_act
   pp = alist(c)
   call dgemm('N','N',nConts,nocc,nvir,1.0d0,C_a(1,nocc+1),nConts,X(1,1,pp),nvir,0.0d0,CU,nConts)
   if (vcm) then
      LU(:,:,c) = 2.0d0*CU
   else
      call dgemm('N','T',nConts,nConts,nocc,1.0d0,CU,nConts,C_a,nConts,0.0d0,T2,nConts)
      Dtr(:,:,c) = 2.0d0*(T2 + transpose(T2))
   endif
enddo
if (vcm) then
   allocate(VC(nConts,nocc,n_act))
   call cpscf_vind_c(nConts, n_act, HF_exchange_frac, nocc, C_a, LU, VC)
else
   call cpscf_vind(nConts, n_act, HF_exchange_frac, Dtr, Jtr, Katr, skip_xc=.true.)
endif
if (cpscf_xc_active) then
   allocate(Xp(nvir,nocc,n_act), Ua(nvir,nocc,n_act))
   do c = 1,n_act
      Ua(:,:,c) = X(:,:,alist(c))
   enddo
   tv0 = omp_get_wtime()
   call xc_fxc_project_vo(nConts, nocc, nvir, n_act, C_a, Ua, cpscf_xc_Ptot, Xp)
   cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
endif
do c = 1,n_act
   pp = alist(c)
   if (vcm) then
      Tmp = VC(:,:,c)
   else
      call dgemm('N','N',nConts,nocc,nConts,1.0d0,Jtr(1,1,c),nConts,C_a,nConts,0.0d0,Tmp,nConts)
   endif
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,Tmp,nConts,0.0d0,Qout(1,1,pp),nvir)
   if (cphf_voo_track) call dgemm('T','N',nocc,nocc,nConts,1.0d0,C_a,nConts,Tmp,nConts,0.0d0,Qoo(1,1,pp),nocc)
   if (cpscf_xc_active) Qout(:,:,pp) = Qout(:,:,pp) + Xp(:,:,c)
   Qout(:,:,pp) = dlt*X(:,:,pp) + Qout(:,:,pp)
enddo
deallocate(Dtr, Jtr, Katr, Tmp, T2, CU)
if (allocated(Xp)) deallocate(Xp, Ua)
if (allocated(VC)) deallocate(VC, LU)
end subroutine apply_A
end subroutine cphf_solve_vo_pcg

module subroutine cphf_solve_vo_direct(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                                 RHSfix, Uvo_all, niter_total)
use mod_integrals, only: direct_mode
use mod_xc_hess, only: xc_fxc_project_vo
use omp_lib, only: omp_get_wtime
implicit none
integer,intent(in) :: nConts, nocc, nvir, npert
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: RHSfix(nvir,nocc,npert)
real(8),intent(out) :: Uvo_all(nvir,nocc,npert)
integer,intent(out) :: niter_total

integer(8), parameter :: cache_budget_bytes = 2_8*1024_8*1024_8
integer(8), parameter :: diis_budget_bytes  = 256_8*1024_8*1024_8
integer :: bchunk, ndiis
integer :: iter, p, i, a, cstart, cend, nb, cidx, nact
integer :: info, nuse, ihist, r, s
real(8) :: dmax_p, tv0
real(8) :: conv_thr
integer, parameter :: NTHR = 5
real(8), parameter :: thr_list(NTHR) = (/1.0d-4,1.0d-6,1.0d-8,1.0d-10,1.0d-12/)
integer :: thr_passes(NTHR), it
logical,allocatable :: thr_hit(:,:)
logical,allocatable :: converged(:)
integer,allocatable :: actlist(:), nhist(:), ipiv(:)
real(8),allocatable :: Dtr(:,:,:), Jtr(:,:,:), Katr(:,:,:)
real(8),allocatable :: Tmp(:,:), Unew(:,:,:)
real(8),allocatable :: dsU(:,:,:,:), dsE(:,:,:,:)
real(8),allocatable :: Bmat(:,:), rhs(:)
real(8),allocatable :: Xproj(:,:,:), Uact(:,:,:)

if (direct_mode) then
   bchunk = npert
else
   bchunk = int(max(1_8, cache_budget_bytes / (2_8*int(nConts,8)*int(nConts,8)*8_8)))
   bchunk = max(1, min(npert, bchunk))
endif

ndiis = int(max(1_8, diis_budget_bytes / &
         (2_8*int(nvir,8)*int(nocc,8)*int(npert,8)*8_8)))
ndiis = max(2, min(10, ndiis))

allocate(Dtr(nConts,nConts,bchunk), Jtr(nConts,nConts,bchunk))
allocate(Katr(nConts,nConts,bchunk))
allocate(Tmp(nConts,nocc))
allocate(Unew(nvir,nocc,npert), converged(npert), actlist(npert), nhist(npert))
allocate(Xproj(nvir,nocc,bchunk), Uact(nvir,nocc,bchunk))
allocate(dsU(nvir,nocc,ndiis,npert), dsE(nvir,nocc,ndiis,npert))

Uvo_all = 0.0d0
converged = .false.
nhist = 0
niter_total = 0
allocate(thr_hit(NTHR,npert))
thr_hit = .false.
thr_passes = 0

conv_thr = 5.0d-5
block
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_CPHF_CONV", envval, status=envstat)
   if (envstat .eq. 0) read(envval,*) conv_thr
   if (cphf_conv_override > 0.0d0) conv_thr = cphf_conv_override
end block
block
   character(16) :: solv
   call get_environment_variable("ENGINE_CPHF_SOLVER", solv)
   if (trim(solv) .ne. "diis") then
      call cphf_solve_vo_pcg(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                             RHSfix, Uvo_all, niter_total, conv_thr)
      return
   endif
end block

do iter = 1,200
   nact = 0
   do p = 1,npert
      if (.not. converged(p)) then
         nact = nact + 1
         actlist(nact) = p
      endif
   enddo
   if (nact .eq. 0) exit
   niter_total = niter_total + nact

   cstart = 1
   do while (cstart .le. nact)
      cend = min(nact, cstart+bchunk-1)
      nb = cend - cstart + 1
      do cidx = 1,nb
         p = actlist(cstart+cidx-1)
         Dtr(:,:,cidx) = 2.0d0*( matmul(C_a(:,nocc+1:nConts), matmul(Uvo_all(:,:,p), transpose(C_a(:,1:nocc)))) &
                       + matmul(C_a(:,1:nocc), matmul(transpose(Uvo_all(:,:,p)), transpose(C_a(:,nocc+1:nConts)))) )
      enddo
      call cpscf_vind(nConts, nb, HF_exchange_frac, Dtr(:,:,1:nb), Jtr(:,:,1:nb), Katr(:,:,1:nb), skip_xc=.true.)
      if (cpscf_xc_active) then
         do cidx = 1,nb
            p = actlist(cstart+cidx-1)
            Uact(:,:,cidx) = Uvo_all(:,:,p)
         enddo
         tv0 = omp_get_wtime()
         call xc_fxc_project_vo(nConts, nocc, nvir, nb, C_a, Uact(:,:,1:nb), cpscf_xc_Ptot, Xproj(:,:,1:nb))
         cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
      endif
      do cidx = 1,nb
         p = actlist(cstart+cidx-1)
         Tmp = matmul(Jtr(:,:,cidx), C_a(:,1:nocc))
         do i = 1,nocc
            do a = 1,nvir
               if (cpscf_xc_active) then
                  Unew(a,i,p) = ( RHSfix(a,i,p) + dot_product(C_a(:,nocc+a), Tmp(:,i)) &
                                  + Xproj(a,i,cidx) ) / (eLev_a(i)-eLev_a(nocc+a))
               else
                  Unew(a,i,p) = ( RHSfix(a,i,p) + dot_product(C_a(:,nocc+a), Tmp(:,i)) ) &
                                / (eLev_a(i)-eLev_a(nocc+a))
               endif
            enddo
         enddo
      enddo
      cstart = cend + 1
   enddo

   do cidx = 1,nact
      p = actlist(cidx)
      dmax_p = maxval(abs(Unew(:,:,p)-Uvo_all(:,:,p)))

      ihist = mod(nhist(p), ndiis) + 1
      dsU(:,:,ihist,p) = Unew(:,:,p)
      dsE(:,:,ihist,p) = Unew(:,:,p) - Uvo_all(:,:,p)
      nhist(p) = nhist(p) + 1
      Uvo_all(:,:,p) = Unew(:,:,p)

      nuse = min(nhist(p), ndiis)
      if (nuse .ge. 2) then
         allocate(Bmat(nuse+1,nuse+1), rhs(nuse+1), ipiv(nuse+1))
         do r = 1,nuse
            do s = 1,nuse
               Bmat(r,s) = sum(dsE(:,:,r,p)*dsE(:,:,s,p))
            enddo
         enddo
         Bmat(1:nuse,nuse+1) = -1.0d0
         Bmat(nuse+1,1:nuse) = -1.0d0
         Bmat(nuse+1,nuse+1) = 0.0d0
         rhs = 0.0d0
         rhs(nuse+1) = -1.0d0
         call dgesv(nuse+1, 1, Bmat, nuse+1, ipiv, rhs, nuse+1, info)
         if (info .eq. 0) then
            Uvo_all(:,:,p) = 0.0d0
            do r = 1,nuse
               Uvo_all(:,:,p) = Uvo_all(:,:,p) + rhs(r)*dsU(:,:,r,p)
            enddo
         endif
         deallocate(Bmat, rhs, ipiv)
      endif

      do it = 1,NTHR
         if (.not. thr_hit(it,p)) then
            thr_passes(it) = thr_passes(it) + 1
            if (dmax_p < thr_list(it)) thr_hit(it,p) = .true.
         endif
      enddo

      if (dmax_p < conv_thr) converged(p) = .true.
   enddo
enddo

block
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_CPHF_THRESHOLD_SCAN", envval, status=envstat)
   if (envstat .eq. 0) then
      print *, ' CPHF passes vs convergence threshold (total over', npert, 'perturbations):'
      do it = 1,NTHR
         print *, '   thr=', thr_list(it), ' passes=', thr_passes(it)
      enddo
   endif
end block
deallocate(thr_hit)

deallocate(Dtr,Jtr,Katr,Tmp,Unew,converged,actlist,nhist,dsU,dsE)
deallocate(Xproj,Uact)
end subroutine cphf_solve_vo_direct

module subroutine cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, &
                          U_full, h1_MO_all, s1_MO_all, mo_e1_all)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force, integrals_compute_force_1e, integrals_compute_force_exchange_lr
use mod_exchange, only: exchange_build, RS_omega, cosx_enabled
use mod_integrals, only: direct_mode
use mod_profile, only: prof_start, prof_stop, prof_add, prof_add_cpu
use mod_xc_hess, only: xc_fxc_project, xc_vxc_deriv_occ, xc_fxc_c, xc_fxc_set_loc, xc_fxc_clear_loc
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: U_full(nConts,nocc,3,natoms)
real(8),intent(out) :: h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms)
real(8),intent(out) :: mo_e1_all(nocc,nocc,3,natoms)

integer :: nvir
integer :: npert, ipert
integer :: pchunk, pstart, pend, pn, pidx
integer(8), parameter :: pchunk_cache_bytes = 2_8*1024_8*1024_8
real(8),allocatable :: s1oo_all(:,:,:)
real(8),allocatable :: Dbat(:,:,:), Jbat(:,:,:), Kabat(:,:,:), VCb(:,:,:), Lb(:,:,:)
logical :: vcm1
real(8) :: tstage
real(8) :: omp_get_wtime
external :: omp_get_wtime
real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: RHSfix(:,:,:), Uvo_all(:,:,:)
real(8),allocatable :: dF_full(:,:), dS_full(:,:), h1_MO(:,:), s1_MO(:,:)
real(8),allocatable :: s1oo(:,:), Pconn(:,:), Kbc(:,:), Tmp(:,:)
real(8),allocatable :: Uvo(:,:), Pvo(:,:), Ptot2(:,:), v1oo(:,:)
real(8),allocatable :: Din1(:,:,:), Jout1(:,:,:)
real(8),allocatable :: Xd_oo(:,:,:), Xd_vo(:,:,:), Xp_vo(:,:,:), Xo1(:,:,:), Xo2(:,:,:)
logical :: xc_proj_mo
integer :: iat,ic,i,j,a,b,mu

nvir = nConts - nocc

cpscf_xc_active = (HF_exchange_frac < 1.0d0)
if (cpscf_xc_active) then
   if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
   allocate(cpscf_xc_Ptot(nConts,nConts))
   cpscf_xc_Ptot = Pa + Pb
endif
if (cpscf_xc_active) then
   block
      use mod_xc_hess, only: xc_resp_grid_build
      use GRID_info, only: XCGRID_TIER_LO
      character(8) :: e
      call get_environment_variable("ENGINE_CPSCF_XC_GRID", e)
      if (trim(e) .ne. "full") call xc_resp_grid_build(XCGRID_TIER_LO)
   end block
endif

cpscf_ecp_atom = -1
allocate(dS(nConts,nConts,3), dHcore(nConts,nConts,3))
allocate(dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3))
tstage = omp_get_wtime()
call cpu_time(tstage_c)
call prof_start("hess_force_int")
if (allocated(cpscf_dfj_jx) .and. (allocated(cpscf_cosx_kxa) .or. HF_exchange_frac == 0.0d0)) then
   call integrals_compute_force_1e(nConts, dS, dHcore)
   dJi = 0.0d0 ; dKa = 0.0d0 ; dKb = 0.0d0
else
   call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
endif
call prof_stop("hess_force_int")
if (RS_omega > 0.0d0 .and. .not. cosx_enabled) then
   if (allocated(cpscf_lr_dKa)) deallocate(cpscf_lr_dKa, cpscf_lr_dKb)
   allocate(cpscf_lr_dKa(nConts,nConts,3), cpscf_lr_dKb(nConts,nConts,3))
   call integrals_compute_force_exchange_lr(nConts, Pa, Pb, cpscf_lr_dKa, cpscf_lr_dKb)
endif
cphf_t_forceint = omp_get_wtime() - tstage

allocate(ao_atom(nConts))
mu = 0
do iat = 1,natoms
   do i = 1,atoms(iat)%nconts
      mu = mu+1
      ao_atom(mu) = iat
   enddo
enddo

cphf_t_amat = 0.0d0

allocate(dF_full(nConts,nConts), dS_full(nConts,nConts))
allocate(h1_MO(nConts,nConts), s1_MO(nConts,nConts))
allocate(s1oo(nocc,nocc), Pconn(nConts,nConts))
allocate(Kbc(nConts,nConts), Tmp(nConts,nocc))
allocate(Uvo(nvir,nocc), Pvo(nConts,nConts), Ptot2(nConts,nConts), v1oo(nocc,nocc))
allocate(Din1(nConts,nConts,1), Jout1(nConts,nConts,1))

U_full = 0.0d0
cphf_t_fockderiv = 0.0d0 ; cphf_c_fockderiv = 0.0d0
cphf_t_iter = 0.0d0
cphf_t_pass1 = 0.0d0
cphf_t_pass3 = 0.0d0
cphf_t_mo = 0.0d0
cphf_n_iter_total = 0

npert = 3*natoms
allocate(RHSfix(nvir,nocc,npert), Uvo_all(nvir,nocc,npert))
allocate(s1oo_all(nocc,nocc,npert))
block
   character(16) :: solv0
   call get_environment_variable("ENGINE_CPHF_SOLVER", solv0)
   cphf_voo_track = (trim(solv0) .ne. "diis") .and. .not. cpscf_xc_active
end block
vcm1 = cpscf_use_vind_c(HF_exchange_frac)
cphf_voo_ready = .false.
if (allocated(cphf_voo_conn)) deallocate(cphf_voo_conn)
if (allocated(cphf_voo_u)) deallocate(cphf_voo_u)
if (cphf_voo_track) allocate(cphf_voo_conn(nocc,nocc,npert), cphf_voo_u(nocc,nocc,npert))
if (direct_mode) then
   pchunk = npert
else
   pchunk = int(max(1_8, pchunk_cache_bytes / (2_8*int(nConts,8)*int(nConts,8)*8_8)))
   pchunk = max(1, min(npert, pchunk))
endif
if (vcm1) then
   allocate(Dbat(1,1,1), Jbat(1,1,1), Kabat(1,1,1))
else
   allocate(Dbat(nConts,nConts,pchunk), Jbat(nConts,nConts,pchunk))
   allocate(Kabat(nConts,nConts,pchunk))
endif
allocate(Xd_oo(nocc,nocc,pchunk), Xd_vo(nvir,nocc,pchunk), Xp_vo(nvir,nocc,pchunk))
allocate(Xo1(nocc,nocc,pchunk), Xo2(nocc,nocc,pchunk))

block
   character(len=8) :: e1, e2
   call get_environment_variable("ENGINE_XC_PROJ_GRP", e1)
   call get_environment_variable("ENGINE_XC_PROJ_CHK", e2)
   xc_proj_mo = (trim(e1) .eq. "1" .or. trim(e2) .eq. "1")
   call get_environment_variable("ENGINE_XC_PROJ_LOC", e1)
   if (HF_exchange_frac < 1.0d0 .and. trim(e1) .eq. "1") call xc_fxc_set_loc(nConts, nocc, C_a)
end block
ipert = 0
if (HF_exchange_frac < 1.0d0) then
   block
      character(len=8) :: e
      call get_environment_variable("ENGINE_XC_VXCDERIV_FULL", e)
      if (trim(e) .ne. "1") then
         allocate(cpscf_xc_vxco(nConts,nocc,3,natoms))
         tstage = omp_get_wtime()
         call cpu_time(tstage_c)
         call prof_start("fd_xc_vxc_deriv")
         call xc_vxc_deriv_occ(nConts, natoms, nocc, Pa+Pb, C_a(:,1:nocc), cpscf_xc_vxco)
         call prof_stop("fd_xc_vxc_deriv")
         cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
         call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
      endif
   end block
endif
do iat = 1,natoms
   do ic = 1,3
      ipert = ipert + 1
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, dF_full, dS_full)
      cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
      call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      block
         real(8) :: FC(nConts,nocc)
         call dgemm('N','N',nConts,nocc,nConts,1.0d0,dF_full,nConts,C_a,nConts,0.0d0,FC,nConts)
         if (cpscf_cosx_kx_proj) FC = FC - HF_exchange_frac*cpscf_cosx_kxa(:,1:nocc,ic,iat)
         if (allocated(cpscf_xc_vxco)) FC = FC + cpscf_xc_vxco(:,:,ic,iat)
         call dgemm('T','N',nConts,nocc,nConts,1.0d0,C_a,nConts,FC,nConts,0.0d0,h1_MO_all(1,1,ic,iat),nConts)
         call dgemm('N','N',nConts,nocc,nConts,1.0d0,dS_full,nConts,C_a,nConts,0.0d0,FC,nConts)
         call dgemm('T','N',nConts,nocc,nConts,1.0d0,C_a,nConts,FC,nConts,0.0d0,s1_MO_all(1,1,ic,iat),nConts)
      end block
      s1_MO(:,1:nocc) = s1_MO_all(:,:,ic,iat)
      cphf_t_mo = cphf_t_mo + (omp_get_wtime() - tstage)

      s1oo_all(:,:,ipert) = s1_MO(1:nocc,1:nocc)
   enddo
enddo

tstage = omp_get_wtime()
call cpu_time(tstage_c)
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      if (.not. vcm1) Dbat(:,:,pidx) = -2.0d0*matmul(C_a(:,1:nocc), &
                        matmul(s1oo_all(:,:,ipert), transpose(C_a(:,1:nocc))))
      Xd_oo(:,:,pidx) = -s1oo_all(:,:,ipert)
   enddo
   if (vcm1) then
      allocate(VCb(nConts,nocc,pn), Lb(nConts,nocc,pn))
      do pidx = 1,pn
         Lb(:,:,pidx) = -matmul(C_a(:,1:nocc), s1oo_all(:,:,pstart+pidx-1))
      enddo
      call cpscf_vind_c(nConts, pn, HF_exchange_frac, nocc, C_a(:,1:nocc), Lb, VCb)
      deallocate(Lb)
   else
      call cpscf_vind(nConts, pn, HF_exchange_frac, Dbat(:,:,1:pn), Jbat(:,:,1:pn), &
                      Kabat(:,:,1:pn), skip_xc=.true.)
   endif
   if (cpscf_xc_active) then
      call xc_fxc_project(nConts, nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                          nvir, nocc, C_a(:,nocc+1:nConts), C_a(:,1:nocc), &
                          pn, Xd_oo(:,:,1:pn), cpscf_xc_Ptot, Xp_vo(:,:,1:pn))
   endif
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO(:,1:nocc) = h1_MO_all(:,:,ic,iat)
      s1_MO(:,1:nocc) = s1_MO_all(:,:,ic,iat)
      if (vcm1) then
         Tmp = VCb(:,:,pidx)
      else
         Tmp = matmul(Jbat(:,:,pidx), C_a(:,1:nocc))
      endif
      if (cphf_voo_track) cphf_voo_conn(:,:,ipert) = matmul(transpose(C_a(:,1:nocc)), Tmp)

      do i = 1,nocc
         do a = 1,nvir
            RHSfix(a,i,ipert) = h1_MO(nocc+a,i) - s1_MO(nocc+a,i)*eLev_a(i) &
                                + dot_product(C_a(:,nocc+a), Tmp(:,i))
            if (cpscf_xc_active) RHSfix(a,i,ipert) = RHSfix(a,i,ipert) + Xp_vo(a,i,pidx)
         enddo
      enddo
   enddo
   if (allocated(VCb)) deallocate(VCb)
   pstart = pend + 1
enddo
cphf_t_pass1 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass1 = ctmp_c - tstage_c
call prof_add("hess_cphf_rhs", cphf_t_pass1)

call prof_add_cpu("hess_cphf_rhs", cphf_c_pass1)
tstage = omp_get_wtime()
call cpu_time(tstage_c)
call cphf_solve_vo_direct(nConts, nocc, nvir, npert, C_a, eLev_a, HF_exchange_frac, &
                          RHSfix, Uvo_all, cphf_n_iter_total)
cphf_t_iter = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_iter = ctmp_c - tstage_c
call prof_add("hess_cphf_solve", cphf_t_iter)

call prof_add_cpu("hess_cphf_solve", cphf_c_iter)
tstage = omp_get_wtime()
call cpu_time(tstage_c)
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      Uvo = Uvo_all(:,:,ipert)
      s1oo = s1oo_all(:,:,ipert)
      Pconn = -2.0d0*matmul(C_a(:,1:nocc), matmul(s1oo, transpose(C_a(:,1:nocc))))
      Pvo = 2.0d0*( matmul(C_a(:,nocc+1:nConts), matmul(Uvo, transpose(C_a(:,1:nocc)))) &
                  + matmul(C_a(:,1:nocc), matmul(transpose(Uvo), transpose(C_a(:,nocc+1:nConts)))) )
      if (.not. vcm1) Dbat(:,:,pidx) = Pconn + Pvo
      Xd_oo(:,:,pidx) = -s1oo
      Xd_vo(:,:,pidx) = 2.0d0*Uvo

      do i = 1,nocc
         do j = 1,nocc
            U_full(j,i,ic,iat) = -0.5d0*s1oo(j,i)
         enddo
         do a = 1,nvir
            U_full(nocc+a,i,ic,iat) = Uvo(a,i)
         enddo
      enddo
   enddo
   if (.not. cphf_voo_ready) then
      if (vcm1) then
         allocate(VCb(nConts,nocc,pn), Lb(nConts,nocc,pn))
         do pidx = 1,pn
            ipert = pstart + pidx - 1
            Lb(:,:,pidx) = -matmul(C_a(:,1:nocc), s1oo_all(:,:,ipert)) &
                           + 2.0d0*matmul(C_a(:,nocc+1:nConts), Uvo_all(:,:,ipert))
         enddo
         call cpscf_vind_c(nConts, pn, HF_exchange_frac, nocc, C_a(:,1:nocc), Lb, VCb)
         deallocate(Lb)
      else
         call cpscf_vind(nConts, pn, HF_exchange_frac, Dbat(:,:,1:pn), Jbat(:,:,1:pn), &
                         Kabat(:,:,1:pn), skip_xc=.true.)
      endif
   endif
   if (cpscf_xc_active) then
      if (xc_proj_mo) then
         call xc_fxc_project(nConts, nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                             nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                             pn, Xd_oo(:,:,1:pn), cpscf_xc_Ptot, Xo1(:,:,1:pn))
         call xc_fxc_project(nConts, nvir, nocc, C_a(:,nocc+1:nConts), C_a(:,1:nocc), &
                             nocc, nocc, C_a(:,1:nocc), C_a(:,1:nocc), &
                             pn, Xd_vo(:,:,1:pn), cpscf_xc_Ptot, Xo2(:,:,1:pn))
      else
         block
            real(8),allocatable :: Lx(:,:,:), Xt(:,:,:)
            allocate(Lx(nConts,nocc,pn), Xt(nConts,nocc,pn))
            do pidx = 1,pn
               call dgemm('N','N',nConts,nocc,nocc,1.0d0,C_a,nConts,Xd_oo(1,1,pidx),nocc,0.0d0,Lx(1,1,pidx),nConts)
               call dgemm('N','N',nConts,nocc,nvir,1.0d0,C_a(1,nocc+1),nConts,Xd_vo(1,1,pidx),nvir,1.0d0,Lx(1,1,pidx),nConts)
            enddo
            call xc_fxc_c(nConts, nocc, C_a(:,1:nocc), pn, Lx, cpscf_xc_Ptot, Xt)
            do pidx = 1,pn
               call dgemm('T','N',nocc,nocc,nConts,1.0d0,C_a,nConts,Xt(1,1,pidx),nConts,0.0d0,Xo1(1,1,pidx),nocc)
            enddo
            Xo2(:,:,1:pn) = 0.0d0
            deallocate(Lx, Xt)
         end block
      endif
   endif
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO(:,1:nocc) = h1_MO_all(:,:,ic,iat)
      s1_MO(:,1:nocc) = s1_MO_all(:,:,ic,iat)
      if (cphf_voo_ready) then
         v1oo = cphf_voo_conn(:,:,ipert) + cphf_voo_u(:,:,ipert)
      else if (vcm1) then
         v1oo = matmul(transpose(C_a(:,1:nocc)), VCb(:,:,pidx))
      else
         v1oo = matmul(transpose(C_a(:,1:nocc)), matmul(Jbat(:,:,pidx), C_a(:,1:nocc)))
      endif
      if (cpscf_xc_active) v1oo = v1oo + Xo1(:,:,pidx) + Xo2(:,:,pidx)
      do i = 1,nocc
         do j = 1,nocc
            mo_e1_all(i,j,ic,iat) = h1_MO(i,j) - 0.5d0*(eLev_a(i)+eLev_a(j))*s1_MO(i,j) + v1oo(i,j)
         enddo
      enddo
   enddo
   if (allocated(VCb)) deallocate(VCb)
   pstart = pend + 1
enddo
cphf_t_pass3 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass3 = ctmp_c - tstage_c
call prof_add("hess_cphf_obs", cphf_t_pass3)
call prof_add_cpu("hess_cphf_obs", cphf_c_pass3)
block
   use mod_xc_hess, only: xc_resp_grid_clear
   call xc_resp_grid_clear()
end block

call prof_add("hess_fock_deriv", cphf_t_fockderiv)

call prof_add_cpu("hess_fock_deriv", cphf_c_fockderiv)
deallocate(dS,dHcore,dJi,dKa,dKb,ao_atom,RHSfix,Uvo_all,s1oo_all)
if (allocated(cpscf_lr_dKa)) deallocate(cpscf_lr_dKa, cpscf_lr_dKb)
deallocate(Dbat,Jbat,Kabat)
deallocate(Xd_oo,Xd_vo,Xp_vo,Xo1,Xo2)
deallocate(dF_full,dS_full,h1_MO,s1_MO,s1oo,Pconn,Kbc,Tmp)
if (allocated(cpscf_xc_vxco)) deallocate(cpscf_xc_vxco)
call xc_fxc_clear_loc()
deallocate(Uvo,Pvo,Ptot2,v1oo)
deallocate(Din1,Jout1)
if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
cpscf_xc_active = .false.
end subroutine cphf_solve_rhf

end submodule cpscf_rhf_impl
