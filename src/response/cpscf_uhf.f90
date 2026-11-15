! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: the unrestricted (UHF/UKS) response - the spin-resolved

submodule (mod_cpscf) cpscf_uhf_impl
implicit none
contains

subroutine cphf_fock_deriv_atom_spin(nConts, ao_atom, iatom, ic, &
                                      dS, dHcore, dJi, dKa, dKb, HF_exchange_frac, Pa, Pb, &
                                      dF_a, dF_b, dS_full)
use MOL_info, only: atoms
use mod_integrals, only: integrals_nuc_attraction_deriv, integrals_cphf_2e_density_role
use mod_profile, only: prof_start, prof_stop
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_Fx
implicit none
integer,intent(in) :: nConts, ao_atom(nConts), iatom, ic
real(8),intent(in) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)
real(8),intent(in) :: dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: dF_a(nConts,nConts), dF_b(nConts,nConts), dS_full(nConts,nConts)
INCLUDE 'parameter.h'

integer :: mu,nu
real(8),allocatable :: Dr(:,:,:)
real(8) :: Zeff, kfrac, jfac

kfrac = HF_exchange_frac
if (allocated(cpscf_cosx_kxa)) kfrac = 0.0d0
jfac = 1.0d0
if (allocated(cpscf_dfj_jx)) jfac = 0.0d0

dF_a = 0.0d0
dF_b = 0.0d0
dS_full = 0.0d0
do nu = 1,nConts
   do mu = 1,nConts
      if (ao_atom(mu) == iatom) then
         dF_a(mu,nu) = dF_a(mu,nu) - (dHcore(mu,nu,ic) + jfac*dJi(mu,nu,ic) - kfrac*dKa(mu,nu,ic))
         dF_b(mu,nu) = dF_b(mu,nu) - (dHcore(mu,nu,ic) + jfac*dJi(mu,nu,ic) - kfrac*dKb(mu,nu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(mu,nu,ic)
      endif
      if (ao_atom(nu) == iatom) then
         dF_a(mu,nu) = dF_a(mu,nu) - (dHcore(nu,mu,ic) + jfac*dJi(nu,mu,ic) - kfrac*dKa(nu,mu,ic))
         dF_b(mu,nu) = dF_b(mu,nu) - (dHcore(nu,mu,ic) + jfac*dJi(nu,mu,ic) - kfrac*dKb(nu,mu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(nu,mu,ic)
      endif
   enddo
enddo

allocate(Dr(nConts,nConts,3))
call prof_start("fd_nuc_attr_deriv")
call integrals_nuc_attraction_deriv(atoms(iatom)%coor*ans2bohr, nConts, Dr)
call prof_stop("fd_nuc_attr_deriv")
Zeff = atoms(iatom)%charge - atoms(iatom)%ecpCoreElec
dF_a = dF_a - Zeff*(Dr(:,:,ic) + transpose(Dr(:,:,ic)))
dF_b = dF_b - Zeff*(Dr(:,:,ic) + transpose(Dr(:,:,ic)))
call cpscf_ecp_centre(nConts, iatom, Dr)
if (cpscf_ecp_has) then
   dF_a = dF_a + Dr(:,:,ic) + transpose(Dr(:,:,ic))
   dF_b = dF_b + Dr(:,:,ic) + transpose(Dr(:,:,ic))
endif
deallocate(Dr)

if (.not. (jfac == 0.0d0 .and. kfrac == 0.0d0)) then
block
   integer, save :: cache_iatom = -1, cache_nConts = -1
   real(8), allocatable, save :: dJ2_cache(:,:,:), dK2a_cache(:,:,:), dK2b_cache(:,:,:)
   real(8), allocatable, save :: Pa_cache(:,:), Pb_cache(:,:)
   logical :: cache_valid

   cache_valid = allocated(dJ2_cache) .and. (cache_iatom == iatom) .and. (cache_nConts == nConts)
   if (cache_valid) cache_valid = all(Pa_cache == Pa) .and. all(Pb_cache == Pb)

   if (.not. cache_valid) then
      if (allocated(dJ2_cache)) deallocate(dJ2_cache, dK2a_cache, dK2b_cache, Pa_cache, Pb_cache)
      allocate(dJ2_cache(nConts,nConts,3), dK2a_cache(nConts,nConts,3), dK2b_cache(nConts,nConts,3))
      allocate(Pa_cache(nConts,nConts), Pb_cache(nConts,nConts))
      call integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2_cache, dK2a_cache, dK2b_cache)
      cache_iatom = iatom
      cache_nConts = nConts
      Pa_cache = Pa
      Pb_cache = Pb
   endif

   do nu = 1,nConts
      do mu = 1,nConts
         dF_a(mu,nu) = dF_a(mu,nu) &
            - (jfac*dJ2_cache(mu,nu,ic) - kfrac*dK2a_cache(mu,nu,ic)) &
            - (jfac*dJ2_cache(nu,mu,ic) - kfrac*dK2a_cache(nu,mu,ic))
         dF_b(mu,nu) = dF_b(mu,nu) &
            - (jfac*dJ2_cache(mu,nu,ic) - kfrac*dK2b_cache(mu,nu,ic)) &
            - (jfac*dJ2_cache(nu,mu,ic) - kfrac*dK2b_cache(nu,mu,ic))
      enddo
   enddo
end block
endif
if (allocated(cpscf_lr_dKa)) then
   block
      use mod_exchange, only: RS_beta
      real(8),allocatable :: dK2a(:,:,:), dK2b(:,:,:)
      allocate(dK2a(nConts,nConts,3), dK2b(nConts,nConts,3))
      call cpscf_lr_density_role(nConts, iatom, Pa, Pb, dK2a, dK2b)
      do nu = 1,nConts
         do mu = 1,nConts
            if (ao_atom(mu) == iatom) then
               dF_a(mu,nu) = dF_a(mu,nu) + RS_beta*cpscf_lr_dKa(mu,nu,ic)
               dF_b(mu,nu) = dF_b(mu,nu) + RS_beta*cpscf_lr_dKb(mu,nu,ic)
            endif
            if (ao_atom(nu) == iatom) then
               dF_a(mu,nu) = dF_a(mu,nu) + RS_beta*cpscf_lr_dKa(nu,mu,ic)
               dF_b(mu,nu) = dF_b(mu,nu) + RS_beta*cpscf_lr_dKb(nu,mu,ic)
            endif
            dF_a(mu,nu) = dF_a(mu,nu) + RS_beta*(dK2a(mu,nu,ic) + dK2a(nu,mu,ic))
            dF_b(mu,nu) = dF_b(mu,nu) + RS_beta*(dK2b(mu,nu,ic) + dK2b(nu,mu,ic))
         enddo
      enddo
      deallocate(dK2a, dK2b)
   end block
endif
if (allocated(cpscf_dfj_jx)) then
   dF_a = dF_a + cpscf_dfj_jx(:,:,ic,iatom)
   dF_b = dF_b + cpscf_dfj_jx(:,:,ic,iatom)
endif
if (allocated(cpscf_cosx_kxa) .and. .not. cpscf_cosx_kx_proj) then
   dF_a = dF_a - HF_exchange_frac*cpscf_cosx_kxa(:,:,ic,iatom)
   dF_b = dF_b - HF_exchange_frac*cpscf_cosx_kxb(:,:,ic,iatom)
endif
if (cosmo_hess_active) then
   dF_a = dF_a + cosmo_hess_Fx(:,:,ic,iatom)
   dF_b = dF_b + cosmo_hess_Fx(:,:,ic,iatom)
endif
end subroutine cphf_fock_deriv_atom_spin

module subroutine cpscf_vind_spin(nConts, nbatch, nspin, HF_exchange_frac, Dtrial, Dtot, Vout, Kscr)
use mod_integrals, only: integrals_build_exchange_lr_batch, integrals_build_jk_batch_uhf, integrals_build_coulomb_df_batch, &
                         integrals_build_coulomb_batch, integrals_build_exchange_batch, &
                         integrals_build_exchange_cosx, cosx_set_response_mode, engine_use_df_j
use mod_exchange, only: cosx_enabled
use mod_xc_hess, only: xc_fxc_apply_u
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_add_spin
implicit none
integer,intent(in) :: nConts, nbatch, nspin
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Dtrial(nConts,nConts,nbatch,nspin)
real(8),intent(in) :: Dtot(nConts,nConts,nbatch)
real(8),intent(inout) :: Vout(nConts,nConts,nbatch,nspin)
real(8),intent(inout) :: Kscr(nConts,nConts,nbatch,nspin)
integer :: b

if (nspin .ne. 2) then
   print *, 'ERROR: cpscf_vind_spin only implements nspin=2 (UHF); got', nspin
   stop 1
endif

if ((.not. cosx_enabled .or. HF_exchange_frac == 0.0d0) .and. &
    (.not. cpscf_dfj_response())) then
   call integrals_build_jk_batch_uhf(nConts, nbatch, Dtot, &
             Dtrial(:,:,:,1), Dtrial(:,:,:,2), &
             Vout(:,:,:,1), Kscr(:,:,:,1), Kscr(:,:,:,2))
   Vout(:,:,:,2) = Vout(:,:,:,1)
   Vout(:,:,:,1) = Vout(:,:,:,1) - HF_exchange_frac*Kscr(:,:,:,1)
   Vout(:,:,:,2) = Vout(:,:,:,2) - HF_exchange_frac*Kscr(:,:,:,2)
   if (allocated(cpscf_lr_dKa)) then
      block
         use mod_exchange, only: RS_beta
         real(8),allocatable :: Kla(:,:,:), Klb(:,:,:)
         allocate(Kla(nConts,nConts,nbatch), Klb(nConts,nConts,nbatch))
         call integrals_build_exchange_lr_batch(nConts, nbatch, Dtrial(:,:,1:nbatch,1), Dtrial(:,:,1:nbatch,2), Kla, Klb)
         Vout(:,:,1:nbatch,1) = Vout(:,:,1:nbatch,1) - RS_beta*Kla
         Vout(:,:,1:nbatch,2) = Vout(:,:,1:nbatch,2) - RS_beta*Klb
         deallocate(Kla, Klb)
      end block
   endif
   if (cpscf_xc_uks) call xc_fxc_apply_u(nConts, nbatch, cpscf_xc_Pa, cpscf_xc_Pb, &
                                         Dtrial(:,:,:,1), Dtrial(:,:,:,2), Vout(:,:,:,1), Vout(:,:,:,2))
if (cosmo_hess_active) call cosmo_hess_add_spin(nConts, nbatch, Dtot, Vout)
   return
endif

if (cpscf_dfj_response()) then
   call integrals_build_coulomb_df_batch(nConts, nbatch, Dtot, Vout(:,:,:,1))
else
   call integrals_build_coulomb_batch(nConts, nbatch, Dtot, Vout(:,:,:,1))
endif

if (cosx_enabled) then
   if (HF_exchange_frac /= 0.0d0) then
      call cpscf_cosx_k(nConts, 2*nbatch, Dtrial, Kscr)
   else
      Kscr = 0.0d0
   endif
else
   if (HF_exchange_frac /= 0.0d0) then
      call integrals_build_exchange_batch(nConts, nbatch, Dtrial(:,:,:,1), Dtrial(:,:,:,2), &
                                           Kscr(:,:,:,1), Kscr(:,:,:,2))
   else
      Kscr = 0.0d0
   endif
endif
Vout(:,:,:,2) = Vout(:,:,:,1)
Vout(:,:,:,1) = Vout(:,:,:,1) - HF_exchange_frac*Kscr(:,:,:,1)
Vout(:,:,:,2) = Vout(:,:,:,2) - HF_exchange_frac*Kscr(:,:,:,2)
if (cpscf_xc_uks) call xc_fxc_apply_u(nConts, nbatch, cpscf_xc_Pa, cpscf_xc_Pb, &
                                      Dtrial(:,:,:,1), Dtrial(:,:,:,2), Vout(:,:,:,1), Vout(:,:,:,2))
if (cosmo_hess_active) call cosmo_hess_add_spin(nConts, nbatch, Dtot, Vout)
end subroutine cpscf_vind_spin

subroutine cpscf_solve_vo_direct_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, &
                                       C_a, eLev_a, C_b, eLev_b, HF_exchange_frac, &
                                       RHSfix_a, RHSfix_b, Uvo_a, Uvo_b, niter_total)
use mod_integrals, only: direct_mode
implicit none
integer,intent(in) :: nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: C_b(nConts,nConts), eLev_b(nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: RHSfix_a(nvir_a,nocc_a,npert), RHSfix_b(nvir_b,nocc_b,npert)
real(8),intent(out) :: Uvo_a(nvir_a,nocc_a,npert), Uvo_b(nvir_b,nocc_b,npert)
integer,intent(out) :: niter_total

integer(8), parameter :: cache_budget_bytes = 2_8*1024_8*1024_8
integer(8), parameter :: diis_budget_bytes  = 256_8*1024_8*1024_8
integer :: bchunk, ndiis, nuvo
integer :: iter, p, i, a, cstart, cend, nb, cidx, nact
integer :: info, nuse, ihist, r, s, off
real(8) :: dmax_p
real(8) :: conv_thr
integer, parameter :: NTHR = 5
real(8), parameter :: thr_list(NTHR) = (/1.0d-4,1.0d-6,1.0d-8,1.0d-10,1.0d-12/)
integer :: thr_passes(NTHR), it
logical,allocatable :: thr_hit(:,:)
logical,allocatable :: converged(:)
integer,allocatable :: actlist(:), nhist(:), ipiv(:)
real(8),allocatable :: Dtr(:,:,:,:), Dtot(:,:,:), Vtr(:,:,:,:), Kscr(:,:,:,:)
real(8),allocatable :: Tmpa(:,:), Tmpb(:,:), Unew_a(:,:,:), Unew_b(:,:,:)
real(8),allocatable :: dsU(:,:,:), dsE(:,:,:), flat_u(:), flat_unew(:), flat_err(:)
real(8),allocatable :: Bmat(:,:), rhs(:)

nuvo = nvir_a*nocc_a + nvir_b*nocc_b

if (direct_mode) then
   bchunk = npert
else
   bchunk = int(max(1_8, cache_budget_bytes / (7_8*int(nConts,8)*int(nConts,8)*8_8)))
   bchunk = max(1, min(npert, bchunk))
endif

ndiis = int(max(1_8, diis_budget_bytes / (2_8*int(nuvo,8)*int(npert,8)*8_8)))
ndiis = max(2, min(10, ndiis))

allocate(Dtr(nConts,nConts,bchunk,2), Dtot(nConts,nConts,bchunk))
allocate(Vtr(nConts,nConts,bchunk,2), Kscr(nConts,nConts,bchunk,2))
allocate(Tmpa(nConts,nocc_a), Tmpb(nConts,nocc_b))
allocate(Unew_a(nvir_a,nocc_a,npert), Unew_b(nvir_b,nocc_b,npert))
allocate(converged(npert), actlist(npert), nhist(npert))
allocate(dsU(nuvo,ndiis,npert), dsE(nuvo,ndiis,npert))
allocate(flat_u(nuvo), flat_unew(nuvo), flat_err(nuvo))

Uvo_a = 0.0d0
Uvo_b = 0.0d0
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
         Dtr(:,:,cidx,1) = matmul(C_a(:,nocc_a+1:nConts), matmul(Uvo_a(:,:,p), transpose(C_a(:,1:nocc_a)))) &
                         + matmul(C_a(:,1:nocc_a), matmul(transpose(Uvo_a(:,:,p)), transpose(C_a(:,nocc_a+1:nConts))))
         Dtr(:,:,cidx,2) = matmul(C_b(:,nocc_b+1:nConts), matmul(Uvo_b(:,:,p), transpose(C_b(:,1:nocc_b)))) &
                         + matmul(C_b(:,1:nocc_b), matmul(transpose(Uvo_b(:,:,p)), transpose(C_b(:,nocc_b+1:nConts))))
         Dtot(:,:,cidx) = Dtr(:,:,cidx,1) + Dtr(:,:,cidx,2)
      enddo
      call cpscf_vind_spin(nConts, nb, 2, HF_exchange_frac, Dtr(:,:,1:nb,:), &
                           Dtot(:,:,1:nb), Vtr(:,:,1:nb,:), Kscr(:,:,1:nb,:))
      do cidx = 1,nb
         p = actlist(cstart+cidx-1)
         Tmpa = matmul(Vtr(:,:,cidx,1), C_a(:,1:nocc_a))
         Tmpb = matmul(Vtr(:,:,cidx,2), C_b(:,1:nocc_b))
         do i = 1,nocc_a
            do a = 1,nvir_a
               Unew_a(a,i,p) = ( RHSfix_a(a,i,p) + dot_product(C_a(:,nocc_a+a), Tmpa(:,i)) ) &
                               / (eLev_a(i)-eLev_a(nocc_a+a))
            enddo
         enddo
         do i = 1,nocc_b
            do a = 1,nvir_b
               Unew_b(a,i,p) = ( RHSfix_b(a,i,p) + dot_product(C_b(:,nocc_b+a), Tmpb(:,i)) ) &
                               / (eLev_b(i)-eLev_b(nocc_b+a))
            enddo
         enddo
      enddo
      cstart = cend + 1
   enddo

   do cidx = 1,nact
      p = actlist(cidx)
      dmax_p = max( maxval(abs(Unew_a(:,:,p)-Uvo_a(:,:,p))), &
                    maxval(abs(Unew_b(:,:,p)-Uvo_b(:,:,p))) )

      flat_unew(1:nvir_a*nocc_a) = reshape(Unew_a(:,:,p), (/nvir_a*nocc_a/))
      flat_u(1:nvir_a*nocc_a)    = reshape(Uvo_a(:,:,p),  (/nvir_a*nocc_a/))
      off = nvir_a*nocc_a
      flat_unew(off+1:nuvo) = reshape(Unew_b(:,:,p), (/nvir_b*nocc_b/))
      flat_u(off+1:nuvo)    = reshape(Uvo_b(:,:,p),  (/nvir_b*nocc_b/))
      flat_err = flat_unew - flat_u

      ihist = mod(nhist(p), ndiis) + 1
      dsU(:,ihist,p) = flat_unew
      dsE(:,ihist,p) = flat_err
      nhist(p) = nhist(p) + 1
      Uvo_a(:,:,p) = Unew_a(:,:,p)
      Uvo_b(:,:,p) = Unew_b(:,:,p)

      nuse = min(nhist(p), ndiis)
      if (nuse .ge. 2) then
         allocate(Bmat(nuse+1,nuse+1), rhs(nuse+1), ipiv(nuse+1))
         do r = 1,nuse
            do s = 1,nuse
               Bmat(r,s) = dot_product(dsE(:,r,p), dsE(:,s,p))
            enddo
         enddo
         Bmat(1:nuse,nuse+1) = -1.0d0
         Bmat(nuse+1,1:nuse) = -1.0d0
         Bmat(nuse+1,nuse+1) = 0.0d0
         rhs = 0.0d0
         rhs(nuse+1) = -1.0d0
         call dgesv(nuse+1, 1, Bmat, nuse+1, ipiv, rhs, nuse+1, info)
         if (info .eq. 0) then
            flat_u = 0.0d0
            do r = 1,nuse
               flat_u = flat_u + rhs(r)*dsU(:,r,p)
            enddo
            Uvo_a(:,:,p) = reshape(flat_u(1:nvir_a*nocc_a), (/nvir_a,nocc_a/))
            Uvo_b(:,:,p) = reshape(flat_u(off+1:nuvo), (/nvir_b,nocc_b/))
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

deallocate(Dtr,Dtot,Vtr,Kscr,Tmpa,Tmpb,Unew_a,Unew_b,converged,actlist,nhist,dsU,dsE)
deallocate(flat_u,flat_unew,flat_err)
end subroutine cpscf_solve_vo_direct_spin

module subroutine cpscf_vind_spin_c(nConts, nb, HF_exchange_frac, noa, nob, Ca, Cb, La, Lb, VCa, VCb)
use mod_exchange, only: RS_beta
use mod_integrals, only: integrals_build_exchange_lr_batch, &
   integrals_build_coulomb_df_batch_c, integrals_build_exchange_cosx_batch_c, &
                         integrals_build_jk_batch_uhf, integrals_build_coulomb_batch, integrals_build_exchange_batch, &
                         integrals_build_coulomb_df_batch
use mod_exchange, only: cosx_enabled
use mod_xc_hess, only: xc_fxc_c_u
use omp_lib, only: omp_get_wtime
use mod_cosmo_hessian, only: cosmo_hess_active, cosmo_hess_vind_spin_c
implicit none
integer,intent(in) :: nConts, nb, noa, nob
real(8),intent(in) :: HF_exchange_frac, Ca(nConts,noa), Cb(nConts,nob)
real(8),intent(in) :: La(nConts,noa,nb), Lb(nConts,nob,nb)
real(8),intent(out) :: VCa(nConts,noa,nb), VCb(nConts,nob,nb)
real(8),allocatable :: Lc(:,:,:), Cc(:,:), JC(:,:,:), KC(:,:,:), Xa(:,:,:), Xb(:,:,:)
real(8),allocatable :: Da(:,:,:), Db(:,:,:), Dt(:,:,:), Jt(:,:,:), Ka(:,:,:), Kb(:,:,:), Kla(:,:,:), Klb(:,:,:)
integer :: b, c0, nch, ib2
logical :: rij, kcosx
real(8) :: tv0
rij = cpscf_dfj_response()
kcosx = cosx_enabled .and. HF_exchange_frac /= 0.0d0
if (rij .and. (kcosx .or. HF_exchange_frac == 0.0d0)) then
tv0 = omp_get_wtime()
allocate(Lc(nConts,noa+nob,nb), Cc(nConts,noa+nob), JC(nConts,noa+nob,nb))
Cc(:,1:noa) = Ca ; Cc(:,noa+1:noa+nob) = Cb
do b = 1,nb
   Lc(:,1:noa,b) = La(:,:,b) ; Lc(:,noa+1:noa+nob,b) = Lb(:,:,b)
enddo
call integrals_build_coulomb_df_batch_c(nConts, nb, Lc, noa+nob, Cc, JC)
do b = 1,nb
   VCa(:,:,b) = JC(:,1:noa,b) ; VCb(:,:,b) = JC(:,noa+1:noa+nob,b)
enddo
deallocate(Lc, Cc, JC)
cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
if (HF_exchange_frac /= 0.0d0) then
   tv0 = omp_get_wtime()
   allocate(KC(nConts,max(noa,nob),nb))
   call integrals_build_exchange_cosx_batch_c(nConts, nb, La, noa, Ca, KC(:,1:noa,:))
   VCa = VCa - HF_exchange_frac*KC(:,1:noa,:)
   call integrals_build_exchange_cosx_batch_c(nConts, nb, Lb, nob, Cb, KC(:,1:nob,:))
   VCb = VCb - HF_exchange_frac*KC(:,1:nob,:)
   deallocate(KC)
   cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   cphf_n_vk = cphf_n_vk + 2
   cphf_nd_vk = cphf_nd_vk + 2*nb
endif
else
   tv0 = omp_get_wtime()
   allocate(Da(nConts,nConts,16), Db(nConts,nConts,16), Dt(nConts,nConts,16), Jt(nConts,nConts,16))
   allocate(Ka(nConts,nConts,16), Kb(nConts,nConts,16), Kla(nConts,nConts,16), Klb(nConts,nConts,16))
   do c0 = 0, nb-1, 16
      nch = min(16, nb-c0)
      do ib2 = 1,nch
         call dgemm('N','T',nConts,nConts,noa,1.0d0,La(1,1,c0+ib2),nConts,Ca,nConts,0.0d0,Da(1,1,ib2),nConts)
         Da(:,:,ib2) = Da(:,:,ib2) + transpose(Da(:,:,ib2))
         call dgemm('N','T',nConts,nConts,nob,1.0d0,Lb(1,1,c0+ib2),nConts,Cb,nConts,0.0d0,Db(1,1,ib2),nConts)
         Db(:,:,ib2) = Db(:,:,ib2) + transpose(Db(:,:,ib2))
         Dt(:,:,ib2) = Da(:,:,ib2) + Db(:,:,ib2)
      enddo
      if (.not. rij .and. .not. kcosx) then
         call integrals_build_jk_batch_uhf(nConts, nch, Dt(:,:,1:nch), Da(:,:,1:nch), Db(:,:,1:nch), &
                                           Jt(:,:,1:nch), Ka(:,:,1:nch), Kb(:,:,1:nch))
      else
         if (rij) then
            call integrals_build_coulomb_df_batch(nConts, nch, Dt(:,:,1:nch), Jt(:,:,1:nch))
         else
            call integrals_build_coulomb_batch(nConts, nch, Dt(:,:,1:nch), Jt(:,:,1:nch))
         endif
         Ka(:,:,1:nch) = 0.0d0 ; Kb(:,:,1:nch) = 0.0d0
         if (HF_exchange_frac /= 0.0d0 .and. .not. kcosx) &
            call integrals_build_exchange_batch(nConts, nch, Da(:,:,1:nch), Db(:,:,1:nch), Ka(:,:,1:nch), Kb(:,:,1:nch))
      endif
      if (allocated(cpscf_lr_dKa)) &
         call integrals_build_exchange_lr_batch(nConts, nch, Da(:,:,1:nch), Db(:,:,1:nch), Kla(:,:,1:nch), Klb(:,:,1:nch))
      do ib2 = 1,nch
         Jt(:,:,ib2) = Jt(:,:,ib2) - HF_exchange_frac*Ka(:,:,ib2)
         if (allocated(cpscf_lr_dKa)) Jt(:,:,ib2) = Jt(:,:,ib2) - RS_beta*Kla(:,:,ib2)
         call dgemm('N','N',nConts,noa,nConts,1.0d0,Jt(1,1,ib2),nConts,Ca,nConts,0.0d0,VCa(1,1,c0+ib2),nConts)
         Jt(:,:,ib2) = Jt(:,:,ib2) + HF_exchange_frac*(Ka(:,:,ib2) - Kb(:,:,ib2))
         if (allocated(cpscf_lr_dKa)) Jt(:,:,ib2) = Jt(:,:,ib2) + RS_beta*(Kla(:,:,ib2) - Klb(:,:,ib2))
         call dgemm('N','N',nConts,nob,nConts,1.0d0,Jt(1,1,ib2),nConts,Cb,nConts,0.0d0,VCb(1,1,c0+ib2),nConts)
      enddo
   enddo
   deallocate(Da, Db, Dt, Jt, Ka, Kb, Kla, Klb)
   cphf_t_vj = cphf_t_vj + (omp_get_wtime() - tv0)
   if (kcosx) then
      tv0 = omp_get_wtime()
      allocate(KC(nConts,max(noa,nob),nb))
      call integrals_build_exchange_cosx_batch_c(nConts, nb, La, noa, Ca, KC(:,1:noa,:))
      VCa = VCa - HF_exchange_frac*KC(:,1:noa,:)
      call integrals_build_exchange_cosx_batch_c(nConts, nb, Lb, nob, Cb, KC(:,1:nob,:))
      VCb = VCb - HF_exchange_frac*KC(:,1:nob,:)
      deallocate(KC)
      cphf_t_vk = cphf_t_vk + (omp_get_wtime() - tv0)
   endif
endif
if (cpscf_xc_uks) then
   tv0 = omp_get_wtime()
   allocate(Xa(nConts,noa,nb), Xb(nConts,nob,nb))
   call xc_fxc_c_u(nConts, noa, nob, Ca, Cb, nb, La, Lb, cpscf_xc_Pa, cpscf_xc_Pb, Xa, Xb)
   VCa = VCa + Xa ; VCb = VCb + Xb
   deallocate(Xa, Xb)
   cphf_t_vxc = cphf_t_vxc + (omp_get_wtime() - tv0)
endif
if (cosmo_hess_active) call cosmo_hess_vind_spin_c(nConts, nb, noa, nob, Ca, Cb, La, Lb, VCa, VCb)
end subroutine cpscf_vind_spin_c

subroutine cpscf_solve_vo_pcg_spin(nConts, noa, nva, nob, nvb, npert, Ca, ea, Cb, eb, &
                                   HF_exchange_frac, RHSa, RHSb, Ua, Ub, niter_total, conv_thr)
implicit none
integer,intent(in) :: nConts, noa, nva, nob, nvb, npert
real(8),intent(in) :: Ca(nConts,nConts), ea(nConts), Cb(nConts,nConts), eb(nConts)
real(8),intent(in) :: HF_exchange_frac, conv_thr
real(8),intent(in) :: RHSa(nva,noa,npert), RHSb(nvb,nob,npert)
real(8),intent(out) :: Ua(nva,noa,npert), Ub(nvb,nob,npert)
integer,intent(out) :: niter_total
real(8),allocatable :: ma(:,:), mb(:,:), da(:,:), db(:,:), Ra(:,:,:), Rb(:,:,:), Za(:,:,:), Zb(:,:,:)
real(8),allocatable :: Pa_(:,:,:), Pb_(:,:,:), Qa(:,:,:), Qb(:,:,:), rz(:)
logical,allocatable :: conv(:)
integer,allocatable :: act(:)
integer :: p, it, nact, k, a, i
real(8) :: alpha, beta, rznew, pq
logical :: verb
character(len=8) :: ev
verb = (engine_verbose .ge. 2)
allocate(da(nva,noa), db(nvb,nob), rz(npert), conv(npert), act(npert))
allocate(Ra(nva,noa,npert), Za(nva,noa,npert), Pa_(nva,noa,npert), Qa(nva,noa,npert))
allocate(Rb(nvb,nob,npert), Zb(nvb,nob,npert), Pb_(nvb,nob,npert), Qb(nvb,nob,npert))
do i = 1,noa
   do a = 1,nva
      da(a,i) = ea(noa+a) - ea(i)
   enddo
enddo
do i = 1,nob
   do a = 1,nvb
      db(a,i) = eb(nob+a) - eb(i)
   enddo
enddo
if (verb) print '(A,2ES11.3)', '  UHF PCG min orbital-energy gaps a / b:', minval(da), minval(db)
block
   real(8) :: pfl
   character(len=16) :: e
   pfl = 0.0d0
   call get_environment_variable("ENGINE_CPHF_PFLOOR", e)
   if (len_trim(e) .gt. 0) read(e,*) pfl
   allocate(ma(nva,noa), mb(nvb,nob))
   ma = max(da, pfl) ; mb = max(db, pfl)
end block
block
   real(8) :: gthr
   character(len=16) :: e
   integer,parameter :: NMAXD = 64
   integer :: nd, ks, ka_(NMAXD), ki_(NMAXD), kspin(NMAXD), q
   real(8),allocatable :: La1(:,:,:), Lb1(:,:,:), V1a(:,:,:), V1b(:,:,:)
   gthr = 0.05d0
   call get_environment_variable("ENGINE_CPHF_DIAGGAP", e)
   if (len_trim(e) .gt. 0) read(e,*) gthr
   nd = 0
   do i = 1,noa
      do a = 1,nva
         if (da(a,i) .lt. gthr .and. nd .lt. NMAXD) then
            nd = nd + 1 ; ka_(nd) = a ; ki_(nd) = i ; kspin(nd) = 1
         endif
      enddo
   enddo
   do i = 1,nob
      do a = 1,nvb
         if (db(a,i) .lt. gthr .and. nd .lt. NMAXD) then
            nd = nd + 1 ; ka_(nd) = a ; ki_(nd) = i ; kspin(nd) = 2
         endif
      enddo
   enddo
   if (nd .gt. 0) then
      allocate(La1(nConts,noa,nd), Lb1(nConts,nob,nd), V1a(nConts,noa,nd), V1b(nConts,nob,nd))
      La1 = 0.0d0 ; Lb1 = 0.0d0
      do q = 1,nd
         if (kspin(q) .eq. 1) then
            La1(:,ki_(q),q) = Ca(:,noa+ka_(q))
         else
            Lb1(:,ki_(q),q) = Cb(:,nob+ka_(q))
         endif
      enddo
      call cpscf_vind_spin_c(nConts, nd, HF_exchange_frac, noa, nob, Ca(:,1:noa), Cb(:,1:nob), La1, Lb1, V1a, V1b)
      do q = 1,nd
         if (kspin(q) .eq. 1) then
            ks = ka_(q)
            ma(ks,ki_(q)) = max(ma(ks,ki_(q)), da(ks,ki_(q)) + sum(Ca(:,noa+ks)*V1a(:,ki_(q),q)))
         else
            ks = ka_(q)
            mb(ks,ki_(q)) = max(mb(ks,ki_(q)), db(ks,ki_(q)) + sum(Cb(:,nob+ks)*V1b(:,ki_(q),q)))
         endif
      enddo
      if (verb) print '(A,I0,A)', '  UHF PCG exact diagonal for ', nd, ' small-gap pairs'
      deallocate(La1, Lb1, V1a, V1b)
   endif
end block
do p = 1,npert
   Ua(:,:,p) = -RHSa(:,:,p)/ma
   Ub(:,:,p) = -RHSb(:,:,p)/mb
   act(p) = p
enddo
call apply_A_spin(npert, act, Ua, Ub, Qa, Qb)
niter_total = npert
do p = 1,npert
   Ra(:,:,p) = -RHSa(:,:,p) - Qa(:,:,p) ; Rb(:,:,p) = -RHSb(:,:,p) - Qb(:,:,p)
   Za(:,:,p) = Ra(:,:,p)/ma ; Zb(:,:,p) = Rb(:,:,p)/mb
   Pa_(:,:,p) = Za(:,:,p) ; Pb_(:,:,p) = Zb(:,:,p)
   rz(p) = sum(Ra(:,:,p)*Za(:,:,p)) + sum(Rb(:,:,p)*Zb(:,:,p))
   conv(p) = max(maxval(abs(Za(:,:,p))), maxval(abs(Zb(:,:,p)))) < conv_thr
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
   if (verb) print '(A,I4,A,I5,A,ES10.3)', '  UHF PCG it', it, '  active', nact, '  max|z|', &
      maxval((/ (max(maxval(abs(Za(:,:,act(k)))), maxval(abs(Zb(:,:,act(k))))), k=1,nact) /))
   niter_total = niter_total + nact
   call apply_A_spin(nact, act, Pa_, Pb_, Qa, Qb)
   do k = 1,nact
      p = act(k)
      pq = sum(Pa_(:,:,p)*Qa(:,:,p)) + sum(Pb_(:,:,p)*Qb(:,:,p))
      alpha = rz(p)/pq
      Ua(:,:,p) = Ua(:,:,p) + alpha*Pa_(:,:,p) ; Ub(:,:,p) = Ub(:,:,p) + alpha*Pb_(:,:,p)
      Ra(:,:,p) = Ra(:,:,p) - alpha*Qa(:,:,p) ; Rb(:,:,p) = Rb(:,:,p) - alpha*Qb(:,:,p)
      Za(:,:,p) = Ra(:,:,p)/ma ; Zb(:,:,p) = Rb(:,:,p)/mb
      rznew = sum(Ra(:,:,p)*Za(:,:,p)) + sum(Rb(:,:,p)*Zb(:,:,p))
      beta = rznew/rz(p)
      rz(p) = rznew
      Pa_(:,:,p) = Za(:,:,p) + beta*Pa_(:,:,p) ; Pb_(:,:,p) = Zb(:,:,p) + beta*Pb_(:,:,p)
      conv(p) = max(maxval(abs(Za(:,:,p))), maxval(abs(Zb(:,:,p)))) < conv_thr
   enddo
enddo
deallocate(ma, mb, da, db, Ra, Rb, Za, Zb, Pa_, Pb_, Qa, Qb, rz, conv, act)

contains

subroutine apply_A_spin(n_act, alist, Xa, Xb, Qa_, Qb_)
integer,intent(in) :: n_act, alist(n_act)
real(8),intent(in) :: Xa(nva,noa,npert), Xb(nvb,nob,npert)
real(8),intent(inout) :: Qa_(nva,noa,npert), Qb_(nvb,nob,npert)
real(8),allocatable :: La(:,:,:), Lb(:,:,:), VCa(:,:,:), VCb(:,:,:)
integer :: c, pp
allocate(La(nConts,noa,n_act), Lb(nConts,nob,n_act), VCa(nConts,noa,n_act), VCb(nConts,nob,n_act))
!$omp parallel do private(c,pp)
do c = 1,n_act
   pp = alist(c)
   call dgemm('N','N',nConts,noa,nva,1.0d0,Ca(1,noa+1),nConts,Xa(1,1,pp),nva,0.0d0,La(1,1,c),nConts)
   call dgemm('N','N',nConts,nob,nvb,1.0d0,Cb(1,nob+1),nConts,Xb(1,1,pp),nvb,0.0d0,Lb(1,1,c),nConts)
enddo
!$omp end parallel do
call cpscf_vind_spin_c(nConts, n_act, HF_exchange_frac, noa, nob, Ca(:,1:noa), Cb(:,1:nob), La, Lb, VCa, VCb)
!$omp parallel do private(c,pp)
do c = 1,n_act
   pp = alist(c)
   call dgemm('T','N',nva,noa,nConts,1.0d0,Ca(1,noa+1),nConts,VCa(1,1,c),nConts,0.0d0,Qa_(1,1,pp),nva)
   call dgemm('T','N',nvb,nob,nConts,1.0d0,Cb(1,nob+1),nConts,VCb(1,1,c),nConts,0.0d0,Qb_(1,1,pp),nvb)
   Qa_(:,:,pp) = Qa_(:,:,pp) + da*Xa(:,:,pp)
   Qb_(:,:,pp) = Qb_(:,:,pp) + db*Xb(:,:,pp)
enddo
!$omp end parallel do
deallocate(La, Lb, VCa, VCb)
end subroutine apply_A_spin
end subroutine cpscf_solve_vo_pcg_spin

subroutine cpscf_solve_vo_bcg_spin(nConts, noa, nva, nob, nvb, npert, Ca, ea, Cb, eb, &
                                   HF_exchange_frac, RHSa, RHSb, Ua, Ub, niter_total, conv_thr)
implicit none
integer,intent(in) :: nConts, noa, nva, nob, nvb, npert
real(8),intent(in) :: Ca(nConts,nConts), ea(nConts), Cb(nConts,nConts), eb(nConts)
real(8),intent(in) :: HF_exchange_frac, conv_thr
real(8),intent(in) :: RHSa(nva,noa,npert), RHSb(nvb,nob,npert)
real(8),intent(out) :: Ua(nva,noa,npert), Ub(nvb,nob,npert)
integer,intent(out) :: niter_total
integer :: na, nt, it, k, p, i, a, nk, info, lwork, q
real(8),allocatable :: dm(:), X(:,:), R(:,:), Z(:,:), P_(:,:), Q_(:,:), G(:,:), ev(:), work(:), T(:,:), W(:,:)
logical,allocatable :: conv(:)
logical :: verb
character(len=8) :: e
external :: dgemm, dsyev
verb = (engine_verbose .ge. 2)
na = nva*noa ; nt = na + nvb*nob
allocate(dm(nt), X(nt,npert), R(nt,npert), Z(nt,npert), P_(nt,npert), Q_(nt,npert), conv(npert))
do i = 1,noa
   do a = 1,nva
      dm(a + (i-1)*nva) = ea(noa+a) - ea(i)
   enddo
enddo
do i = 1,nob
   do a = 1,nvb
      dm(na + a + (i-1)*nvb) = eb(nob+a) - eb(i)
   enddo
enddo
do p = 1,npert
   R(1:na,p) = -reshape(RHSa(:,:,p),(/na/)) ; R(na+1:nt,p) = -reshape(RHSb(:,:,p),(/nt-na/))
   X(:,p) = R(:,p)/dm
enddo
call applyA(npert, X, Q_)
niter_total = npert
R = R - Q_
nk = npert
do it = 1,200
   do p = 1,npert
      Z(:,p) = R(:,p)/dm
      conv(p) = maxval(abs(Z(:,p))) < conv_thr
      if (conv(p)) Z(:,p) = 0.0d0
   enddo
   if (verb) print '(A,I4,A,I5,A,ES10.3)', '  UHF BCG it', it, '  unconverged', count(.not. conv), &
      '  max|z|', maxval(abs(R/spread(dm,2,npert)))
   if (all(conv)) exit
   if (it .eq. 1) then
      P_(:,1:npert) = Z ; nk = npert
   else
      allocate(T(nk,npert))
      call dgemm('T','N',nk,npert,nt,1.0d0,Q_,nt,Z,nt,0.0d0,T,nk)
      allocate(W(nt,npert))
      W = Z
      call dgemm('N','N',nt,npert,nk,-1.0d0,P_,nt,T,nk,1.0d0,W,nt)
      P_(:,1:npert) = W ; nk = npert
      deallocate(T, W)
   endif
   q = 0
   do p = 1,npert
      if (.not. conv(p)) then
         q = q + 1
         if (q .ne. p) P_(:,q) = P_(:,p)
      endif
   enddo
   nk = q
   call applyA(nk, P_, Q_)
   niter_total = niter_total + nk
   allocate(G(nk,nk), ev(nk), work(1))
   call dgemm('T','N',nk,nk,nt,1.0d0,P_,nt,Q_,nt,0.0d0,G,nk)
   G = 0.5d0*(G + transpose(G))
   call dsyev('V','U',nk,G,nk,ev,work,-1,info)
   lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
   call dsyev('V','U',nk,G,nk,ev,work,lwork,info)
   q = 0
   do k = 1,nk
      if (ev(k) .gt. 1.0d-12*ev(nk)) then
         q = q + 1
         G(:,q) = G(:,k)/sqrt(ev(k))
      endif
   enddo
   allocate(W(nt,q))
   call dgemm('N','N',nt,q,nk,1.0d0,P_,nt,G,nk,0.0d0,W,nt)
   P_(:,1:q) = W
   call dgemm('N','N',nt,q,nk,1.0d0,Q_,nt,G,nk,0.0d0,W,nt)
   Q_(:,1:q) = W
   deallocate(W, G, ev, work)
   nk = q
   allocate(T(nk,npert))
   call dgemm('T','N',nk,npert,nt,1.0d0,P_,nt,R,nt,0.0d0,T,nk)
   call dgemm('N','N',nt,npert,nk,1.0d0,P_,nt,T,nk,1.0d0,X,nt)
   call dgemm('N','N',nt,npert,nk,-1.0d0,Q_,nt,T,nk,1.0d0,R,nt)
   deallocate(T)
enddo
do p = 1,npert
   Ua(:,:,p) = reshape(X(1:na,p),(/nva,noa/)) ; Ub(:,:,p) = reshape(X(na+1:nt,p),(/nvb,nob/))
enddo
deallocate(dm, X, R, Z, P_, Q_, conv)

contains

subroutine applyA(n, Pin, Qo)
integer,intent(in) :: n
real(8),intent(in) :: Pin(nt,*)
real(8),intent(inout) :: Qo(nt,*)
real(8),allocatable :: La(:,:,:), Lb(:,:,:), VCa(:,:,:), VCb(:,:,:), Ta(:,:), Tb(:,:)
integer :: c
allocate(La(nConts,noa,n), Lb(nConts,nob,n), VCa(nConts,noa,n), VCb(nConts,nob,n))
!$omp parallel do private(c)
do c = 1,n
   call dgemm('N','N',nConts,noa,nva,1.0d0,Ca(1,noa+1),nConts,Pin(1,c),nva,0.0d0,La(1,1,c),nConts)
   call dgemm('N','N',nConts,nob,nvb,1.0d0,Cb(1,nob+1),nConts,Pin(na+1,c),nvb,0.0d0,Lb(1,1,c),nConts)
enddo
!$omp end parallel do
call cpscf_vind_spin_c(nConts, n, HF_exchange_frac, noa, nob, Ca(:,1:noa), Cb(:,1:nob), La, Lb, VCa, VCb)
!$omp parallel do private(c)
do c = 1,n
   call dgemm('T','N',nva,noa,nConts,1.0d0,Ca(1,noa+1),nConts,VCa(1,1,c),nConts,0.0d0,Qo(1,c),nva)
   call dgemm('T','N',nvb,nob,nConts,1.0d0,Cb(1,nob+1),nConts,VCb(1,1,c),nConts,0.0d0,Qo(na+1,c),nvb)
   Qo(1:nt,c) = Qo(1:nt,c) + dm*Pin(1:nt,c)
enddo
!$omp end parallel do
deallocate(La, Lb, VCa, VCb)
end subroutine applyA
end subroutine cpscf_solve_vo_bcg_spin

module subroutine cpscf_solve_uhf(nConts, natoms, nocc_a, nocc_b, C_a, eLev_a, C_b, eLev_b, &
                           Pa, Pb, HF_exchange_frac, &
                           U_full_a, U_full_b, h1_MO_all_a, s1_MO_all_a, mo_e1_all_a, &
                           h1_MO_all_b, s1_MO_all_b, mo_e1_all_b)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force, integrals_compute_force_1e, integrals_compute_force_exchange_lr
use mod_integrals, only: direct_mode
use mod_profile, only: prof_start, prof_stop, prof_add, prof_add_cpu
use mod_xc_hess, only: xc_vxc_deriv_occ_u
use mod_exchange, only: cosx_enabled, RS_omega
implicit none
integer,intent(in) :: nConts, natoms, nocc_a, nocc_b
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: C_b(nConts,nConts), eLev_b(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: U_full_a(nConts,nocc_a,3,natoms), U_full_b(nConts,nocc_b,3,natoms)
real(8),intent(out) :: h1_MO_all_a(nConts,nocc_a,3,natoms), s1_MO_all_a(nConts,nocc_a,3,natoms)
real(8),intent(out) :: mo_e1_all_a(nocc_a,nocc_a,3,natoms)
real(8),intent(out) :: h1_MO_all_b(nConts,nocc_b,3,natoms), s1_MO_all_b(nConts,nocc_b,3,natoms)
real(8),intent(out) :: mo_e1_all_b(nocc_b,nocc_b,3,natoms)

integer :: nvir_a, nvir_b
integer :: npert, ipert
integer :: pchunk, pstart, pend, pn, pidx
integer(8), parameter :: pchunk_cache_bytes = 2_8*1024_8*1024_8
real(8),allocatable :: s1oo_all_a(:,:,:), s1oo_all_b(:,:,:)
real(8) :: tstage
real(8) :: omp_get_wtime
external :: omp_get_wtime
real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: RHSfix_a(:,:,:), RHSfix_b(:,:,:), Uvo_all_a(:,:,:), Uvo_all_b(:,:,:)
real(8),allocatable :: dF_a(:,:), dF_b(:,:), dS_full(:,:)
real(8),allocatable :: h1_MO_a(:,:), s1_MO_a(:,:), h1_MO_b(:,:), s1_MO_b(:,:)
real(8),allocatable :: s1oo_a(:,:), s1oo_b(:,:), Tmpa(:,:), Tmpb(:,:)
real(8),allocatable :: Uvo_a(:,:), Uvo_b(:,:), Pvo_a(:,:), Pvo_b(:,:), v1oo_a(:,:), v1oo_b(:,:)
real(8),allocatable :: Dtr(:,:,:,:), Dtot(:,:,:), Vtr(:,:,:,:), Kscr(:,:,:,:)
integer :: iat,ic,i,j,a,mu
logical :: fm
real(8),allocatable :: Lfa(:,:,:), Lfb(:,:,:), VCa(:,:,:), VCb(:,:,:)
real(8) :: conv_thr

nvir_a = nConts - nocc_a
nvir_b = nConts - nocc_b

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
allocate(dF_a(nConts,nConts), dF_b(nConts,nConts), dS_full(nConts,nConts))
allocate(h1_MO_a(nConts,nConts), s1_MO_a(nConts,nConts))
allocate(h1_MO_b(nConts,nConts), s1_MO_b(nConts,nConts))
allocate(s1oo_a(nocc_a,nocc_a), s1oo_b(nocc_b,nocc_b))
allocate(Tmpa(nConts,nocc_a), Tmpb(nConts,nocc_b))

U_full_a = 0.0d0
U_full_b = 0.0d0
cphf_t_fockderiv = 0.0d0 ; cphf_c_fockderiv = 0.0d0
cphf_t_iter = 0.0d0
if (HF_exchange_frac < 1.0d0) then
   cpscf_xc_uks = .true.
   if (allocated(cpscf_xc_Pa)) deallocate(cpscf_xc_Pa, cpscf_xc_Pb)
   allocate(cpscf_xc_Pa(nConts,nConts), cpscf_xc_Pb(nConts,nConts))
   cpscf_xc_Pa = Pa ; cpscf_xc_Pb = Pb
   if (allocated(cpscf_xc_vxca)) deallocate(cpscf_xc_vxca, cpscf_xc_vxcb)
   allocate(cpscf_xc_vxca(nConts,nocc_a,3,natoms), cpscf_xc_vxcb(nConts,nocc_b,3,natoms))
   tstage = omp_get_wtime()
   call cpu_time(tstage_c)
   call prof_start("fd_xc_vxc_deriv")
   call xc_vxc_deriv_occ_u(nConts, natoms, nocc_a, nocc_b, Pa, Pb, C_a(:,1:nocc_a), C_b(:,1:nocc_b), &
                           cpscf_xc_vxca, cpscf_xc_vxcb)
   call prof_stop("fd_xc_vxc_deriv")
   cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
   call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
endif
if (cpscf_xc_uks) then
   block
      use mod_xc_hess, only: xc_resp_grid_build
      use GRID_info, only: XCGRID_TIER_LO
      character(8) :: e
      call get_environment_variable("ENGINE_CPSCF_XC_GRID", e)
      if (trim(e) .ne. "full") call xc_resp_grid_build(XCGRID_TIER_LO)
   end block
endif
cphf_t_pass1 = 0.0d0
cphf_t_pass3 = 0.0d0
cphf_t_mo = 0.0d0
cphf_n_iter_total = 0

npert = 3*natoms
allocate(RHSfix_a(nvir_a,nocc_a,npert), RHSfix_b(nvir_b,nocc_b,npert))
allocate(Uvo_all_a(nvir_a,nocc_a,npert), Uvo_all_b(nvir_b,nocc_b,npert))
allocate(s1oo_all_a(nocc_a,nocc_a,npert), s1oo_all_b(nocc_b,nocc_b,npert))

if (direct_mode) then
   pchunk = npert
else
   pchunk = int(max(1_8, pchunk_cache_bytes / (7_8*int(nConts,8)*int(nConts,8)*8_8)))
   pchunk = max(1, min(npert, pchunk))
endif
block
   character(len=8) :: e
   call get_environment_variable("ENGINE_CPSCF_UHF_DENSE", e)
   fm = (trim(e) .ne. "1")
end block
if (fm) then
   allocate(Dtr(1,1,1,2), Dtot(1,1,1), Vtr(1,1,1,2), Kscr(1,1,1,2))
else
   allocate(Dtr(nConts,nConts,pchunk,2), Dtot(nConts,nConts,pchunk))
   allocate(Vtr(nConts,nConts,pchunk,2), Kscr(nConts,nConts,pchunk,2))
endif

ipert = 0
do iat = 1,natoms
   do ic = 1,3
      ipert = ipert + 1
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      call cphf_fock_deriv_atom_spin(nConts, ao_atom, iat, ic, dS, dHcore, dJi, dKa, dKb, &
                                     HF_exchange_frac, Pa, Pb, dF_a, dF_b, dS_full)
      cphf_t_fockderiv = cphf_t_fockderiv + (omp_get_wtime() - tstage)
      call cpu_time(ctmp_c); cphf_c_fockderiv = cphf_c_fockderiv + (ctmp_c - tstage_c)
      tstage = omp_get_wtime()
      call cpu_time(tstage_c)
      h1_MO_a = matmul(transpose(C_a), matmul(dF_a, C_a))
      s1_MO_a = matmul(transpose(C_a), matmul(dS_full, C_a))
      h1_MO_b = matmul(transpose(C_b), matmul(dF_b, C_b))
      if (cpscf_xc_uks) then
         h1_MO_a(:,1:nocc_a) = h1_MO_a(:,1:nocc_a) + matmul(transpose(C_a), cpscf_xc_vxca(:,:,ic,iat))
         h1_MO_b(:,1:nocc_b) = h1_MO_b(:,1:nocc_b) + matmul(transpose(C_b), cpscf_xc_vxcb(:,:,ic,iat))
      endif
      if (cpscf_cosx_kx_proj) then
         h1_MO_a(:,1:nocc_a) = h1_MO_a(:,1:nocc_a) &
            - HF_exchange_frac*matmul(transpose(C_a), cpscf_cosx_kxa(:,1:nocc_a,ic,iat))
         h1_MO_b(:,1:nocc_b) = h1_MO_b(:,1:nocc_b) &
            - HF_exchange_frac*matmul(transpose(C_b), cpscf_cosx_kxb(:,1:nocc_b,ic,iat))
      endif
      s1_MO_b = matmul(transpose(C_b), matmul(dS_full, C_b))
      h1_MO_all_a(:,:,ic,iat) = h1_MO_a(:,1:nocc_a)
      s1_MO_all_a(:,:,ic,iat) = s1_MO_a(:,1:nocc_a)
      h1_MO_all_b(:,:,ic,iat) = h1_MO_b(:,1:nocc_b)
      s1_MO_all_b(:,:,ic,iat) = s1_MO_b(:,1:nocc_b)
      cphf_t_mo = cphf_t_mo + (omp_get_wtime() - tstage)

      s1oo_all_a(:,:,ipert) = s1_MO_a(1:nocc_a,1:nocc_a)
      s1oo_all_b(:,:,ipert) = s1_MO_b(1:nocc_b,1:nocc_b)
   enddo
enddo

tstage = omp_get_wtime()

call cpu_time(tstage_c)
if (fm) then
   allocate(Lfa(nConts,nocc_a,npert), Lfb(nConts,nocc_b,npert), VCa(nConts,nocc_a,npert), VCb(nConts,nocc_b,npert))
   do ipert = 1,npert
      call dgemm('N','N',nConts,nocc_a,nocc_a,-0.5d0,C_a,nConts,s1oo_all_a(1,1,ipert),nocc_a,0.0d0,Lfa(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_b,nocc_b,-0.5d0,C_b,nConts,s1oo_all_b(1,1,ipert),nocc_b,0.0d0,Lfb(1,1,ipert),nConts)
   enddo
   call cpscf_vind_spin_c(nConts, npert, HF_exchange_frac, nocc_a, nocc_b, C_a(:,1:nocc_a), C_b(:,1:nocc_b), &
                          Lfa, Lfb, VCa, VCb)
   do ipert = 1,npert
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      call dgemm('T','N',nvir_a,nocc_a,nConts,1.0d0,C_a(1,nocc_a+1),nConts,VCa(1,1,ipert),nConts,0.0d0,RHSfix_a(1,1,ipert),nvir_a)
      call dgemm('T','N',nvir_b,nocc_b,nConts,1.0d0,C_b(1,nocc_b+1),nConts,VCb(1,1,ipert),nConts,0.0d0,RHSfix_b(1,1,ipert),nvir_b)
      do i = 1,nocc_a
         do a = 1,nvir_a
            RHSfix_a(a,i,ipert) = RHSfix_a(a,i,ipert) + h1_MO_all_a(nocc_a+a,i,ic,iat) &
                                  - s1_MO_all_a(nocc_a+a,i,ic,iat)*eLev_a(i)
         enddo
      enddo
      do i = 1,nocc_b
         do a = 1,nvir_b
            RHSfix_b(a,i,ipert) = RHSfix_b(a,i,ipert) + h1_MO_all_b(nocc_b+a,i,ic,iat) &
                                  - s1_MO_all_b(nocc_b+a,i,ic,iat)*eLev_b(i)
         enddo
      enddo
   enddo
else
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      Dtr(:,:,pidx,1) = -matmul(C_a(:,1:nocc_a), &
                        matmul(s1oo_all_a(:,:,ipert), transpose(C_a(:,1:nocc_a))))
      Dtr(:,:,pidx,2) = -matmul(C_b(:,1:nocc_b), &
                        matmul(s1oo_all_b(:,:,ipert), transpose(C_b(:,1:nocc_b))))
      Dtot(:,:,pidx) = Dtr(:,:,pidx,1) + Dtr(:,:,pidx,2)
   enddo
   call cpscf_vind_spin(nConts, pn, 2, HF_exchange_frac, Dtr(:,:,1:pn,:), &
                        Dtot(:,:,1:pn), Vtr(:,:,1:pn,:), Kscr(:,:,1:pn,:))
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO_a(:,1:nocc_a) = h1_MO_all_a(:,:,ic,iat)
      s1_MO_a(:,1:nocc_a) = s1_MO_all_a(:,:,ic,iat)
      h1_MO_b(:,1:nocc_b) = h1_MO_all_b(:,:,ic,iat)
      s1_MO_b(:,1:nocc_b) = s1_MO_all_b(:,:,ic,iat)
      Tmpa = matmul(Vtr(:,:,pidx,1), C_a(:,1:nocc_a))
      Tmpb = matmul(Vtr(:,:,pidx,2), C_b(:,1:nocc_b))
      do i = 1,nocc_a
         do a = 1,nvir_a
            RHSfix_a(a,i,ipert) = h1_MO_a(nocc_a+a,i) - s1_MO_a(nocc_a+a,i)*eLev_a(i) &
                                  + dot_product(C_a(:,nocc_a+a), Tmpa(:,i))
         enddo
      enddo
      do i = 1,nocc_b
         do a = 1,nvir_b
            RHSfix_b(a,i,ipert) = h1_MO_b(nocc_b+a,i) - s1_MO_b(nocc_b+a,i)*eLev_b(i) &
                                  + dot_product(C_b(:,nocc_b+a), Tmpb(:,i))
         enddo
      enddo
   enddo
   pstart = pend + 1
enddo
endif
cphf_t_pass1 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass1 = ctmp_c - tstage_c
call prof_add("hess_cphf_rhs", cphf_t_pass1)

call prof_add_cpu("hess_cphf_rhs", cphf_c_pass1)
tstage = omp_get_wtime()

call cpu_time(tstage_c)
if (fm) then
   conv_thr = 5.0d-5
   block
      character(32) :: envval
      integer :: envstat
      call get_environment_variable("ENGINE_CPHF_CONV", envval, status=envstat)
      if (envstat .eq. 0) read(envval,*) conv_thr
   end block
   block
      character(len=16) :: solv
      call get_environment_variable("ENGINE_CPHF_SOLVER", solv)
      if (trim(solv) .eq. "pcg") then
         call cpscf_solve_vo_pcg_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, C_a, eLev_a, C_b, eLev_b, &
                                      HF_exchange_frac, RHSfix_a, RHSfix_b, Uvo_all_a, Uvo_all_b, &
                                      cphf_n_iter_total, conv_thr)
      else
         call cpscf_solve_vo_bcg_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, C_a, eLev_a, C_b, eLev_b, &
                                      HF_exchange_frac, RHSfix_a, RHSfix_b, Uvo_all_a, Uvo_all_b, &
                                      cphf_n_iter_total, conv_thr)
      endif
   end block
else
   call cpscf_solve_vo_direct_spin(nConts, nocc_a, nvir_a, nocc_b, nvir_b, npert, &
                                   C_a, eLev_a, C_b, eLev_b, HF_exchange_frac, &
                                   RHSfix_a, RHSfix_b, Uvo_all_a, Uvo_all_b, cphf_n_iter_total)
endif
cphf_t_iter = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_iter = ctmp_c - tstage_c
call prof_add("hess_cphf_solve", cphf_t_iter)

call prof_add_cpu("hess_cphf_solve", cphf_c_iter)
allocate(Uvo_a(nvir_a,nocc_a), Uvo_b(nvir_b,nocc_b))
allocate(Pvo_a(nConts,nConts), Pvo_b(nConts,nConts))
allocate(v1oo_a(nocc_a,nocc_a), v1oo_b(nocc_b,nocc_b))

tstage = omp_get_wtime()

call cpu_time(tstage_c)
if (fm) then
   do ipert = 1,npert
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      call dgemm('N','N',nConts,nocc_a,nocc_a,-0.5d0,C_a,nConts,s1oo_all_a(1,1,ipert),nocc_a,0.0d0,Lfa(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_a,nvir_a,1.0d0,C_a(1,nocc_a+1),nConts,Uvo_all_a(1,1,ipert),nvir_a,1.0d0,Lfa(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_b,nocc_b,-0.5d0,C_b,nConts,s1oo_all_b(1,1,ipert),nocc_b,0.0d0,Lfb(1,1,ipert),nConts)
      call dgemm('N','N',nConts,nocc_b,nvir_b,1.0d0,C_b(1,nocc_b+1),nConts,Uvo_all_b(1,1,ipert),nvir_b,1.0d0,Lfb(1,1,ipert),nConts)
      do i = 1,nocc_a
         do j = 1,nocc_a
            U_full_a(j,i,ic,iat) = -0.5d0*s1oo_all_a(j,i,ipert)
         enddo
         do a = 1,nvir_a
            U_full_a(nocc_a+a,i,ic,iat) = Uvo_all_a(a,i,ipert)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            U_full_b(j,i,ic,iat) = -0.5d0*s1oo_all_b(j,i,ipert)
         enddo
         do a = 1,nvir_b
            U_full_b(nocc_b+a,i,ic,iat) = Uvo_all_b(a,i,ipert)
         enddo
      enddo
   enddo
   call cpscf_vind_spin_c(nConts, npert, HF_exchange_frac, nocc_a, nocc_b, C_a(:,1:nocc_a), C_b(:,1:nocc_b), &
                          Lfa, Lfb, VCa, VCb)
   do ipert = 1,npert
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      call dgemm('T','N',nocc_a,nocc_a,nConts,1.0d0,C_a,nConts,VCa(1,1,ipert),nConts,0.0d0,v1oo_a,nocc_a)
      call dgemm('T','N',nocc_b,nocc_b,nConts,1.0d0,C_b,nConts,VCb(1,1,ipert),nConts,0.0d0,v1oo_b,nocc_b)
      do i = 1,nocc_a
         do j = 1,nocc_a
            mo_e1_all_a(i,j,ic,iat) = h1_MO_all_a(i,j,ic,iat) - 0.5d0*(eLev_a(i)+eLev_a(j))*s1_MO_all_a(i,j,ic,iat) &
                                      + v1oo_a(i,j)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            mo_e1_all_b(i,j,ic,iat) = h1_MO_all_b(i,j,ic,iat) - 0.5d0*(eLev_b(i)+eLev_b(j))*s1_MO_all_b(i,j,ic,iat) &
                                      + v1oo_b(i,j)
         enddo
      enddo
   enddo
   deallocate(Lfa, Lfb, VCa, VCb)
else
pstart = 1
do while (pstart .le. npert)
   pend = min(npert, pstart+pchunk-1)
   pn = pend - pstart + 1
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      Uvo_a = Uvo_all_a(:,:,ipert)
      Uvo_b = Uvo_all_b(:,:,ipert)
      Pvo_a = matmul(C_a(:,nocc_a+1:nConts), matmul(Uvo_a, transpose(C_a(:,1:nocc_a)))) &
            + matmul(C_a(:,1:nocc_a), matmul(transpose(Uvo_a), transpose(C_a(:,nocc_a+1:nConts))))
      Pvo_b = matmul(C_b(:,nocc_b+1:nConts), matmul(Uvo_b, transpose(C_b(:,1:nocc_b)))) &
            + matmul(C_b(:,1:nocc_b), matmul(transpose(Uvo_b), transpose(C_b(:,nocc_b+1:nConts))))
      Dtr(:,:,pidx,1) = -matmul(C_a(:,1:nocc_a), &
                        matmul(s1oo_all_a(:,:,ipert), transpose(C_a(:,1:nocc_a)))) + Pvo_a
      Dtr(:,:,pidx,2) = -matmul(C_b(:,1:nocc_b), &
                        matmul(s1oo_all_b(:,:,ipert), transpose(C_b(:,1:nocc_b)))) + Pvo_b
      Dtot(:,:,pidx) = Dtr(:,:,pidx,1) + Dtr(:,:,pidx,2)

      do i = 1,nocc_a
         do j = 1,nocc_a
            U_full_a(j,i,ic,iat) = -0.5d0*s1oo_all_a(j,i,ipert)
         enddo
         do a = 1,nvir_a
            U_full_a(nocc_a+a,i,ic,iat) = Uvo_a(a,i)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            U_full_b(j,i,ic,iat) = -0.5d0*s1oo_all_b(j,i,ipert)
         enddo
         do a = 1,nvir_b
            U_full_b(nocc_b+a,i,ic,iat) = Uvo_b(a,i)
         enddo
      enddo
   enddo
   call cpscf_vind_spin(nConts, pn, 2, HF_exchange_frac, Dtr(:,:,1:pn,:), &
                        Dtot(:,:,1:pn), Vtr(:,:,1:pn,:), Kscr(:,:,1:pn,:))
   do pidx = 1,pn
      ipert = pstart + pidx - 1
      iat = (ipert-1)/3 + 1
      ic  = ipert - 3*(iat-1)
      h1_MO_a(:,1:nocc_a) = h1_MO_all_a(:,:,ic,iat)
      s1_MO_a(:,1:nocc_a) = s1_MO_all_a(:,:,ic,iat)
      h1_MO_b(:,1:nocc_b) = h1_MO_all_b(:,:,ic,iat)
      s1_MO_b(:,1:nocc_b) = s1_MO_all_b(:,:,ic,iat)
      v1oo_a = matmul(transpose(C_a(:,1:nocc_a)), matmul(Vtr(:,:,pidx,1), C_a(:,1:nocc_a)))
      v1oo_b = matmul(transpose(C_b(:,1:nocc_b)), matmul(Vtr(:,:,pidx,2), C_b(:,1:nocc_b)))
      do i = 1,nocc_a
         do j = 1,nocc_a
            mo_e1_all_a(i,j,ic,iat) = h1_MO_a(i,j) - 0.5d0*(eLev_a(i)+eLev_a(j))*s1_MO_a(i,j) + v1oo_a(i,j)
         enddo
      enddo
      do i = 1,nocc_b
         do j = 1,nocc_b
            mo_e1_all_b(i,j,ic,iat) = h1_MO_b(i,j) - 0.5d0*(eLev_b(i)+eLev_b(j))*s1_MO_b(i,j) + v1oo_b(i,j)
         enddo
      enddo
   enddo
   pstart = pend + 1
enddo
endif
cphf_t_pass3 = omp_get_wtime() - tstage
call cpu_time(ctmp_c); cphf_c_pass3 = ctmp_c - tstage_c
call prof_add("hess_cphf_obs", cphf_t_pass3)

call prof_add_cpu("hess_cphf_obs", cphf_c_pass3)
call prof_add("hess_fock_deriv", cphf_t_fockderiv)
call prof_add_cpu("hess_fock_deriv", cphf_c_fockderiv)
deallocate(dS,dHcore,dJi,dKa,dKb,ao_atom)
if (allocated(cpscf_lr_dKa)) deallocate(cpscf_lr_dKa, cpscf_lr_dKb)
deallocate(RHSfix_a,RHSfix_b,Uvo_all_a,Uvo_all_b,s1oo_all_a,s1oo_all_b)
deallocate(dF_a,dF_b,dS_full,h1_MO_a,s1_MO_a,h1_MO_b,s1_MO_b,s1oo_a,s1oo_b,Tmpa,Tmpb)
deallocate(Uvo_a,Uvo_b,Pvo_a,Pvo_b,v1oo_a,v1oo_b)
deallocate(Dtr,Dtot,Vtr,Kscr)
if (allocated(cpscf_xc_vxca)) deallocate(cpscf_xc_vxca, cpscf_xc_vxcb)
if (allocated(cpscf_xc_Pa)) deallocate(cpscf_xc_Pa, cpscf_xc_Pb)
if (cpscf_xc_uks) then
   block
      use mod_xc_hess, only: xc_resp_grid_clear
      call xc_resp_grid_clear()
   end block
endif
cpscf_xc_uks = .false.
end subroutine cpscf_solve_uhf

end submodule cpscf_uhf_impl
