! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: debug and validation hooks (explicit Amat, residual, dP/dR

submodule (mod_cpscf) cpscf_debug_impl
implicit none
contains

module subroutine cphf_build_amat(nConts, nocc, C_a, HF_exchange_frac, Amat)
use mod_integrals, only: integrals_build_jk_batch
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: Amat(nConts-nocc,nocc,nConts-nocc,nocc)

integer :: nvir
integer :: i,j,a,b,mu,nu,bstart,bend,bn,bidx,bchunk
integer(8), parameter :: cache_budget_bytes = 2_8*1024_8*1024_8
real(8),allocatable :: Dtrial(:,:,:), Jtr(:,:,:), Katr(:,:,:), Rtr(:,:), Tmp(:,:)

nvir = nConts - nocc

bchunk = int(max(1_8, cache_budget_bytes / (2_8*int(nConts,8)*int(nConts,8)*8_8)))
bchunk = max(1, min(nvir, bchunk))

allocate(Dtrial(nConts,nConts,bchunk), Jtr(nConts,nConts,bchunk))
allocate(Katr(nConts,nConts,bchunk))
allocate(Rtr(nConts,nConts), Tmp(nConts,nocc))

do j = 1,nocc
   bstart = 1
   do while (bstart .le. nvir)
      bend = min(nvir, bstart+bchunk-1)
      bn = bend - bstart + 1
      do bidx = 1,bn
         b = bstart + bidx - 1
         do nu = 1,nConts
            do mu = 1,nConts
               Dtrial(mu,nu,bidx) = C_a(mu,nocc+b)*C_a(nu,j) + C_a(nu,nocc+b)*C_a(mu,j)
            enddo
         enddo
      enddo
      call integrals_build_jk_batch(nConts, bn, Dtrial(:,:,1:bn), Jtr(:,:,1:bn), Katr(:,:,1:bn))
      do bidx = 1,bn
         b = bstart + bidx - 1
         Rtr = 2.0d0*Jtr(:,:,bidx) - HF_exchange_frac*Katr(:,:,bidx)
         Tmp = matmul(Rtr, C_a(:,1:nocc))
         do i = 1,nocc
            do a = 1,nvir
               Amat(a,i,b,j) = dot_product(C_a(:,nocc+a), Tmp(:,i))
            enddo
         enddo
      enddo
      bstart = bend + 1
   enddo
enddo
deallocate(Dtrial,Jtr,Katr,Rtr,Tmp)
end subroutine cphf_build_amat

module subroutine cpscf_debug_amat_cosx_symmetry(nConts, nocc, C_a, HF_exchange_frac)
use mod_exchange, only: exchange_build, exchange_build_cosx, cosx_enabled
use mod_integrals, only: integrals_build_coulomb_batch, cosx_set_response_mode
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

integer :: nvir, nov, i, j, a, b, mu, nu
integer :: amax, imax, bmax, jmax
real(8) :: nA_ex, nA_cx, nK_ex, nK_cx
real(8) :: as_ex, as_cx, as_kex, as_kcx
real(8) :: rel_raw, rel_sym, dev, devmax
real(8),allocatable :: A_ex(:,:,:,:), A_cx(:,:,:,:)
real(8),allocatable :: K_ex(:,:,:,:), K_cx(:,:,:,:)
real(8),allocatable :: D3(:,:,:), J3(:,:,:)
real(8),allocatable :: Kex(:,:), Kcx(:,:), Kdisc(:,:)
real(8),allocatable :: Rmat(:,:), Tmp(:,:)

nvir = nConts - nocc
nov  = nvir*nocc

print *, ''
print *, ' ===== COSX CPSCF response-operator symmetry diagnostic ====='
print *, '   nConts =', nConts, '  nocc =', nocc, '  nvir =', nvir
print *, '   explicit A is (nvir*nocc)^2 =', nov, '^2'
if (.not. cosx_enabled) then
   print *, '   ERROR: COSX is not enabled for this run - nothing to compare.'
   print *, '          Re-run with the COSX K-build active (K=''cosx'').'
   return
endif
if (nov > 1500) then
   print *, '   ERROR: nvir*nocc =', nov, ' is too large for an explicit A here.'
   print *, '          Use a small molecule (water/6-31G is the intended size).'
   return
endif

allocate(A_ex(nvir,nocc,nvir,nocc), A_cx(nvir,nocc,nvir,nocc))
allocate(K_ex(nvir,nocc,nvir,nocc), K_cx(nvir,nocc,nvir,nocc))
allocate(D3(nConts,nConts,1), J3(nConts,nConts,1))
allocate(Kex(nConts,nConts), Kcx(nConts,nConts), Kdisc(nConts,nConts))
allocate(Rmat(nConts,nConts), Tmp(nConts,nocc))

do j = 1,nocc
   do b = 1,nvir
      do nu = 1,nConts
         do mu = 1,nConts
            D3(mu,nu,1) = C_a(mu,nocc+b)*C_a(nu,j) + C_a(nu,nocc+b)*C_a(mu,j)
         enddo
      enddo
      call integrals_build_coulomb_batch(nConts, 1, D3, J3)
      call exchange_build(nConts, D3(:,:,1), D3(:,:,1), Kex, Kdisc)
      call cosx_set_response_mode(.true.)
      call exchange_build_cosx(nConts, D3(:,:,1), D3(:,:,1), Kcx, Kdisc, .false.)
      call cosx_set_response_mode(.false.)

      Rmat = 2.0d0*J3(:,:,1) - HF_exchange_frac*Kex
      Tmp = matmul(Rmat, C_a(:,1:nocc))
      do i = 1,nocc
         do a = 1,nvir
            A_ex(a,i,b,j) = dot_product(C_a(:,nocc+a), Tmp(:,i))
         enddo
      enddo
      Rmat = 2.0d0*J3(:,:,1) - HF_exchange_frac*Kcx
      Tmp = matmul(Rmat, C_a(:,1:nocc))
      do i = 1,nocc
         do a = 1,nvir
            A_cx(a,i,b,j) = dot_product(C_a(:,nocc+a), Tmp(:,i))
         enddo
      enddo
      Tmp = matmul(Kex, C_a(:,1:nocc))
      do i = 1,nocc
         do a = 1,nvir
            K_ex(a,i,b,j) = dot_product(C_a(:,nocc+a), Tmp(:,i))
         enddo
      enddo
      Tmp = matmul(Kcx, C_a(:,1:nocc))
      do i = 1,nocc
         do a = 1,nvir
            K_cx(a,i,b,j) = dot_product(C_a(:,nocc+a), Tmp(:,i))
         enddo
      enddo
   enddo
enddo

nA_ex = sqrt(sum(A_ex*A_ex))
nA_cx = sqrt(sum(A_cx*A_cx))
nK_ex = sqrt(sum(K_ex*K_ex))
nK_cx = sqrt(sum(K_cx*K_cx))
as_ex  = amat_asym(A_ex) / max(nA_ex, 1.0d-300)
as_cx  = amat_asym(A_cx) / max(nA_cx, 1.0d-300)
as_kex = amat_asym(K_ex) / max(nK_ex, 1.0d-300)
as_kcx = amat_asym(K_cx) / max(nK_cx, 1.0d-300)

rel_raw = 0.0d0
rel_sym = 0.0d0
devmax = -1.0d0
amax = 0; imax = 0; bmax = 0; jmax = 0
do j = 1,nocc
   do b = 1,nvir
      do i = 1,nocc
         do a = 1,nvir
            dev = A_cx(a,i,b,j) - A_ex(a,i,b,j)
            rel_raw = rel_raw + dev*dev
            dev = 0.5d0*(A_cx(a,i,b,j) + A_cx(b,j,a,i)) - A_ex(a,i,b,j)
            rel_sym = rel_sym + dev*dev
            dev = abs(A_cx(a,i,b,j) - A_cx(b,j,a,i))
            if (dev > devmax) then
               devmax = dev
               amax = a; imax = i; bmax = b; jmax = j
            endif
         enddo
      enddo
   enddo
enddo
rel_raw = sqrt(rel_raw) / max(nA_ex, 1.0d-300)
rel_sym = sqrt(rel_sym) / max(nA_ex, 1.0d-300)

print *, ''
print *, '   ||A||_F            exact =', nA_ex, '   cosx =', nA_cx
print *, '   asym(A)   = ||A-A^T||/||A||'
print *, '                      exact =', as_ex, '   cosx =', as_cx
print *, '   asym(K-only)       exact =', as_kex, '   cosx =', as_kcx
print *, ''
print *, '   rel(A_cosx)        = ||A_cosx - A_exact|| / ||A_exact||  =', rel_raw
print *, '   rel(sym(A_cosx))   = ||sym(A_cosx) - A_exact|| / ||A_exact|| =', rel_sym
print *, '        ^^^ THE DECIDING NUMBER (see this routine''s header)'
print *, ''
print *, '   max |A(a,i,b,j)-A(b,j,a,i)| =', devmax
print *, '        at (a,i,b,j) =', amax, imax, bmax, jmax
print *, ''
print *, '   READING IT:'
print *, '     asym(A) exact should be ~1e-15 (control). If it is not,'
print *, '       this diagnostic itself is broken - fix that first.'
print *, '     rel(sym) << rel(raw), and rel(sym) at COSX energy accuracy'
print *, '       (~1e-4): symmetrization IS the fix - write the adjoint build.'
print *, '     rel(sym) ~ rel(raw): asymmetry is a symptom, not the cause.'
print *, '       Do NOT write the adjoint build; keep diagnosing.'
print *, '   Re-run with ENGINE_COSX_KSCREEN=1e-14 and compare asym(A):'
print *, '     barely moves => screening is not what drives the downstream'
print *, '     Hessian antisymmetry change, and that change needs its own'
print *, '     explanation (possible second, independent bug).'
print *, ' ==========================================================='

deallocate(A_ex,A_cx,K_ex,K_cx,D3,J3,Kex,Kcx,Kdisc,Rmat,Tmp)

contains

function amat_asym(M) result(nrm)
implicit none
real(8),intent(in) :: M(nvir,nocc,nvir,nocc)
real(8) :: nrm, d
integer :: aa,ii,bb,jj
nrm = 0.0d0
do jj = 1,nocc
   do bb = 1,nvir
      do ii = 1,nocc
         do aa = 1,nvir
            d = M(aa,ii,bb,jj) - M(bb,jj,aa,ii)
            nrm = nrm + d*d
         enddo
      enddo
   enddo
enddo
nrm = sqrt(nrm)
end function amat_asym

end subroutine cpscf_debug_amat_cosx_symmetry

module subroutine cphf_debug_amat_batch_check(nConts, nocc, C_a, HF_exchange_frac)
use mod_exchange, only: exchange_build
use mod_integrals, only: integrals_build_coulomb
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
integer :: nvir,i,j,a,b,mu,nu
real(8),allocatable :: Dtrial(:,:), Jtr(:,:), Katr(:,:), Kbtr(:,:), Rtr(:,:), Tmp(:,:)
real(8),allocatable :: Amat_old(:,:,:,:), Amat_new(:,:,:,:)
real(8) :: maxdiff

nvir = nConts - nocc
allocate(Dtrial(nConts,nConts), Jtr(nConts,nConts), Katr(nConts,nConts), Kbtr(nConts,nConts))
allocate(Rtr(nConts,nConts), Tmp(nConts,nocc))
allocate(Amat_old(nvir,nocc,nvir,nocc), Amat_new(nvir,nocc,nvir,nocc))

do j = 1,nocc
   do b = 1,nvir
      do nu = 1,nConts
         do mu = 1,nConts
            Dtrial(mu,nu) = C_a(mu,nocc+b)*C_a(nu,j) + C_a(nu,nocc+b)*C_a(mu,j)
         enddo
      enddo
      call integrals_build_coulomb(nConts, Dtrial, Jtr)
      call exchange_build(nConts, Dtrial, Dtrial, Katr, Kbtr)
      Rtr = 2.0d0*Jtr - HF_exchange_frac*Katr
      Tmp = matmul(Rtr, C_a(:,1:nocc))
      do i = 1,nocc
         do a = 1,nvir
            Amat_old(a,i,b,j) = dot_product(C_a(:,nocc+a), Tmp(:,i))
         enddo
      enddo
   enddo
enddo
deallocate(Dtrial,Jtr,Katr,Kbtr,Rtr,Tmp)

call cphf_build_amat(nConts, nocc, C_a, HF_exchange_frac, Amat_new)

maxdiff = maxval(abs(Amat_old-Amat_new))
print *, "CPHF DEBUG: Amat old-loop vs new-batch max abs diff =", maxdiff
deallocate(Amat_old, Amat_new)
end subroutine cphf_debug_amat_batch_check

module subroutine cphf_debug_batch_scaling(nConts, nocc, C_a)
use mod_integrals, only: integrals_build_exchange_batch
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
integer :: nb, m, mu, nu, bidx
integer, parameter :: nbatch_list(7) = (/1,2,4,8,16,32,64/)
real(8),allocatable :: Da(:,:,:), Db(:,:,:), Ka(:,:,:), Kb(:,:,:)
integer(8) :: wc1,wc2,wr
real(8) :: dt

do nb = 1,7
   m = nbatch_list(nb)
   if (m .gt. nConts-nocc) cycle
   allocate(Da(nConts,nConts,m), Db(nConts,nConts,m), Ka(nConts,nConts,m), Kb(nConts,nConts,m))
   do bidx = 1,m
      do nu = 1,nConts
         do mu = 1,nConts
            Da(mu,nu,bidx) = C_a(mu,nocc+bidx)*C_a(nu,1) + C_a(nu,nocc+bidx)*C_a(mu,1)
         enddo
      enddo
   enddo
   Db = Da
   call system_clock(wc1)
   call integrals_build_exchange_batch(nConts, m, Da, Db, Ka, Kb)
   call system_clock(wc2,wr)
   dt = dble(wc2-wc1)/dble(wr)
   print *, "CPHF DEBUG BATCH SCALING: nbatch=", m, " total(s)=", dt, " per-density(s)=", dt/dble(m)
   deallocate(Da,Db,Ka,Kb)
enddo
end subroutine cphf_debug_batch_scaling

module subroutine cphf_debug_check_residual(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force, integrals_build_coulomb
use mod_exchange, only: exchange_build
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

real(8),allocatable :: U_full(:,:,:,:), h1_MO_all(:,:,:,:), s1_MO_all(:,:,:,:), mo_e1_all(:,:,:,:)
real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: dF_full(:,:), dS_full(:,:), h1_MO(:,:), s1_MO(:,:)
real(8),allocatable :: s1oo(:,:), Pconn(:,:), Pvo(:,:), Jc(:,:), Kac(:,:), Kbc(:,:), Rc(:,:), Tmp(:,:)
integer :: nvir,iat,ic,i,a,mu,nu
real(8) :: maxres, res, Bai

nvir = nConts-nocc
allocate(U_full(nConts,nocc,3,natoms))
allocate(h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms), mo_e1_all(nocc,nocc,3,natoms))
call cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, U_full, h1_MO_all, s1_MO_all, mo_e1_all)
deallocate(h1_MO_all, s1_MO_all, mo_e1_all)

allocate(dS(nConts,nConts,3), dHcore(nConts,nConts,3))
allocate(dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3))
call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
allocate(ao_atom(nConts))
mu = 0
do iat = 1,natoms
   do i = 1,atoms(iat)%nconts
      mu = mu+1
      ao_atom(mu) = iat
   enddo
enddo
allocate(dF_full(nConts,nConts), dS_full(nConts,nConts))
allocate(h1_MO(nConts,nConts), s1_MO(nConts,nConts))
allocate(s1oo(nocc,nocc), Pconn(nConts,nConts), Pvo(nConts,nConts))
allocate(Jc(nConts,nConts), Kac(nConts,nConts), Kbc(nConts,nConts), Rc(nConts,nConts), Tmp(nConts,nocc))

maxres = 0.0d0
do iat = 1,natoms
   do ic = 1,3
      call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, dF_full, dS_full)
      h1_MO = matmul(transpose(C_a), matmul(dF_full, C_a))
      s1_MO = matmul(transpose(C_a), matmul(dS_full, C_a))
      s1oo = s1_MO(1:nocc,1:nocc)
      Pconn = -2.0d0*matmul(C_a(:,1:nocc), matmul(s1oo, transpose(C_a(:,1:nocc))))

      Pvo = 0.0d0
      do i = 1,nocc
         do a = 1,nvir
            do nu=1,nConts
               do mu=1,nConts
                  Pvo(mu,nu) = Pvo(mu,nu) + 2.0d0*U_full(nocc+a,i,ic,iat)* &
                     (C_a(mu,nocc+a)*C_a(nu,i)+C_a(nu,nocc+a)*C_a(mu,i))
               enddo
            enddo
         enddo
      enddo

      call integrals_build_coulomb(nConts, Pconn+Pvo, Jc)
      call exchange_build(nConts, Pconn+Pvo, Pconn+Pvo, Kac, Kbc)
      Rc = Jc - 0.5d0*HF_exchange_frac*Kac
      Tmp = matmul(Rc, C_a(:,1:nocc))

      do i = 1,nocc
         do a = 1,nvir
            Bai = h1_MO(nocc+a,i) - s1_MO(nocc+a,i)*eLev_a(i)
            res = (eLev_a(nocc+a)-eLev_a(i))*U_full(nocc+a,i,ic,iat) + dot_product(C_a(:,nocc+a),Tmp(:,i)) + Bai
            maxres = max(maxres, abs(res))
         enddo
      enddo
   enddo
enddo
print *, 'CPHF DEBUG: independent max residual over all (atom,comp,a,i) = ', maxres

deallocate(U_full,dS,dHcore,dJi,dKa,dKb,ao_atom,dF_full,dS_full,h1_MO,s1_MO)
deallocate(s1oo,Pconn,Pvo,Jc,Kac,Kbc,Rc,Tmp)
end subroutine cphf_debug_check_residual

module subroutine cphf_debug_dump_dP(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
use MOL_info, only: atoms, S
use mod_integrals, only: integrals_compute_force
implicit none
integer,intent(in) :: nConts, natoms, nocc, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

real(8),allocatable :: U_full(:,:,:,:), Cocc_resp(:,:), dP(:,:)
real(8),allocatable :: h1_MO_all(:,:,:,:), s1_MO_all(:,:,:,:), mo_e1_all(:,:,:,:)
real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: dF_full(:,:), dS_full(:,:)
integer :: mu,nu,i

allocate(U_full(nConts,nocc,3,natoms))
allocate(h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms), mo_e1_all(nocc,nocc,3,natoms))
call cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, U_full, h1_MO_all, s1_MO_all, mo_e1_all)
print *, 'CPHF DEBUG h1_MO_all(:,:,ic,iatom):'
do mu = 1,nConts
   print *, h1_MO_all(mu,:,ic,iatom)
enddo
print *, 'CPHF DEBUG s1_MO_all(:,:,ic,iatom):'
do mu = 1,nConts
   print *, s1_MO_all(mu,:,ic,iatom)
enddo
deallocate(h1_MO_all, s1_MO_all, mo_e1_all)

allocate(Cocc_resp(nConts,nocc), dP(nConts,nConts))
Cocc_resp = matmul(C_a, U_full(:,:,ic,iatom))
dP = 2.0d0*(matmul(Cocc_resp, transpose(C_a(:,1:nocc))) + matmul(C_a(:,1:nocc), transpose(Cocc_resp)))

print *, 'CPHF DEBUG dPtot/dR analytic, atom=',iatom,' comp=',ic
do mu = 1,nConts
   print *, dP(mu,:)
enddo

block
   real(8) :: h1_11, s1_11
   real(8),allocatable :: dS2(:,:,:), dHc2(:,:,:), dJ2(:,:,:), dKa2(:,:,:), dKb2(:,:,:)
   integer,allocatable :: ao_atom2(:)
   real(8),allocatable :: dF2(:,:), dS2f(:,:)
   integer :: mu2,i2,nu2
   allocate(dS2(nConts,nConts,3), dHc2(nConts,nConts,3))
   allocate(dJ2(nConts,nConts,3), dKa2(nConts,nConts,3), dKb2(nConts,nConts,3))
   call integrals_compute_force(nConts, Pa, Pb, dS2, dHc2, dJ2, dKa2, dKb2)
   allocate(ao_atom2(nConts))
   mu2 = 0
   do i2 = 1,natoms
      do nu2 = 1,atoms(i2)%nconts
         mu2 = mu2+1
         ao_atom2(mu2) = i2
      enddo
   enddo
   allocate(dF2(nConts,nConts), dS2f(nConts,nConts))
   call cphf_fock_deriv_atom(nConts, ao_atom2, iatom, ic, dS2, dHc2, dJ2, dKa2, HF_exchange_frac, Pa, Pb, dF2, dS2f)
   h1_11 = dot_product(C_a(:,1), matmul(dF2, C_a(:,1)))
   s1_11 = dot_product(C_a(:,1), matmul(dS2f, C_a(:,1)))
   print *, 'CPHF DEBUG deps_1/dR = h1_11 - s1_11*eps_1 =', h1_11 - s1_11*eLev_a(1)
   deallocate(dS2,dHc2,dJ2,dKa2,dKb2,ao_atom2,dF2,dS2f)
end block

allocate(dS(nConts,nConts,3), dHcore(nConts,nConts,3))
allocate(dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3))
call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
allocate(ao_atom(nConts))
mu = 0
do i = 1,natoms
   do nu = 1,atoms(i)%nconts
      mu = mu+1
      ao_atom(mu) = i
   enddo
enddo
allocate(dF_full(nConts,nConts), dS_full(nConts,nConts))
call cphf_fock_deriv_atom(nConts, ao_atom, iatom, ic, dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, dF_full, dS_full)
print *, 'CPHF DEBUG trace check: tr(dP*S)=', sum(dP*S), '  -tr(P*dS_full)=', -sum((Pa+Pb)*dS_full)

deallocate(U_full,Cocc_resp,dP,dS,dHcore,dJi,dKa,dKb,ao_atom,dF_full,dS_full)
end subroutine cphf_debug_dump_dP

module subroutine cphf_debug_grad_check(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: dF_full(:,:), dS_full(:,:), Wmat(:,:), Ptot(:,:)
integer :: iat,ic,i,j,k,mu
real(8) :: g

allocate(dS(nConts,nConts,3), dHcore(nConts,nConts,3))
allocate(dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3))
call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
allocate(ao_atom(nConts))
mu = 0
do iat = 1,natoms
   do i = 1,atoms(iat)%nconts
      mu = mu+1
      ao_atom(mu) = iat
   enddo
enddo

allocate(Wmat(nConts,nConts), Ptot(nConts,nConts))
Ptot = Pa+Pb
Wmat = 0.0d0
do i = 1,nConts
   do j = 1,nConts
      do k = 1,nocc
         Wmat(i,j) = Wmat(i,j) + 2.0d0*eLev_a(k)*C_a(i,k)*C_a(j,k)
      enddo
   enddo
enddo

allocate(dF_full(nConts,nConts), dS_full(nConts,nConts))
block
   real(8),allocatable :: dHc_full(:,:), dJ_full(:,:), dKa_full(:,:), dZero(:,:,:)
   allocate(dHc_full(nConts,nConts), dJ_full(nConts,nConts), dKa_full(nConts,nConts))
   allocate(dZero(nConts,nConts,3)); dZero = 0.0d0
   do iat = 1,natoms
      do ic = 1,3
         call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dHcore, dZero, dZero, 0.0d0, Pa, Pb, dHc_full, dS_full)
         call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dZero, dJi, dZero, 0.0d0, Pa, Pb, dJ_full, dS_full)
         call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dZero, dZero, dKa, -1.0d0, Pa, Pb, dKa_full, dS_full)
         if (iat==1) then
            print *, 'CPHF DEBUG PART tr(P*dHcore_full)=', sum(Ptot*dHc_full), ' comp=',ic
            print *, 'CPHF DEBUG PART tr(W*dS_full)     =', sum(Wmat*dS_full), ' comp=',ic
            print *, 'CPHF DEBUG PART tr(P*dJi_full)    =', sum(Ptot*dJ_full), ' comp=',ic
            print *, 'CPHF DEBUG PART tr(P*dKa_full)    =', sum(Ptot*dKa_full), ' comp=',ic
         endif
      enddo
   enddo
   deallocate(dHc_full,dJ_full,dKa_full,dZero)
end block
do iat = 1,natoms
   do ic = 1,3
      call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, dF_full, dS_full)
      g = sum(Ptot*dF_full) - sum(Wmat*dS_full)
      print *, 'CPHF DEBUG GRAD (electronic, fixed-P, no Enn) atom=',iat,' comp=',ic,' dE/dR=',g
   enddo
enddo

deallocate(dS,dHcore,dJi,dKa,dKb,ao_atom,Wmat,Ptot,dF_full,dS_full)
end subroutine cphf_debug_grad_check

module subroutine cphf_debug_e1_fixed(nConts, Pa, Pb)
use MOL_info, only: Hcore
use mod_integrals, only: integrals_hessian_e1e_fixed_density
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8) :: Eoracle
real(8),allocatable :: Wzero(:,:)
print *, 'CPHF DEBUG E1_FIXED (module Hcore) = ', sum((Pa+Pb)*Hcore)
allocate(Wzero(nConts,nConts)); Wzero = 0.0d0
call integrals_hessian_e1e_fixed_density(nConts, Pa+Pb, Wzero, Eoracle)
print *, 'CPHF DEBUG E1_FIXED (oracle, W=0)   = ', Eoracle
deallocate(Wzero)
end subroutine cphf_debug_e1_fixed

module subroutine cphf_debug_ccbk_raw(nConts, natoms, iatom_c, mu, nu)
use MOL_info, only: atoms
use mod_integrals, only: integrals_hessian_nuc_atom, integrals_hessian_1e_shell, integrals_hessian_nuc_atom_raw
implicit none
integer,intent(in) :: nConts, natoms, iatom_c, mu, nu
INCLUDE 'parameter.h'
real(8),allocatable :: CC(:,:,:), CB(:,:,:), CK(:,:,:)
real(8),allocatable :: s1aa(:,:,:), s1ab(:,:,:), h1aa(:,:,:), h1ab(:,:,:)
real(8),allocatable :: AA(:,:,:), AB(:,:,:)

allocate(CC(nConts,nConts,9), CB(nConts,nConts,9), CK(nConts,nConts,9))
call integrals_hessian_nuc_atom(atoms(iatom_c)%coor*ans2bohr, nConts, CC, CB, CK)
print *, 'CCBK_RAW CC(mu,nu,:)=', CC(mu,nu,:)
print *, 'CCBK_RAW CB(mu,nu,:)=', CB(mu,nu,:)
print *, 'CCBK_RAW CK(mu,nu,:)=', CK(mu,nu,:)
print *, 'CCBK_RAW CC(nu,nu,:)=', CC(nu,nu,:)
deallocate(CC,CB,CK)

allocate(AA(nConts,nConts,9), AB(nConts,nConts,9))
call integrals_hessian_nuc_atom_raw(atoms(iatom_c)%coor*ans2bohr, nConts, AA, AB)
print *, 'CCBK_RAW AA(mu,nu,:)=', AA(mu,nu,:)
print *, 'CCBK_RAW AB(mu,nu,:)=', AB(mu,nu,:)
print *, 'CCBK_RAW AA(nu,mu,:)=', AA(nu,mu,:)
deallocate(AA,AB)

allocate(s1aa(nConts,nConts,9), s1ab(nConts,nConts,9), h1aa(nConts,nConts,9), h1ab(nConts,nConts,9))
call integrals_hessian_1e_shell(nConts, s1aa, s1ab, h1aa, h1ab)
print *, 'CCBK_RAW h1ab(mu,nu,:)=', h1ab(mu,nu,:)
print *, 'CCBK_RAW h1ab(nu,mu,:)=', h1ab(nu,mu,:)
print *, 'CCBK_RAW h1aa(mu,nu,:)=', h1aa(mu,nu,:)
print *, 'CCBK_RAW h1aa(nu,mu,:)=', h1aa(nu,mu,:)
deallocate(s1aa,s1ab,h1aa,h1ab)
end subroutine cphf_debug_ccbk_raw

module subroutine cphf_debug_dF_raw(nConts, natoms, iatom, ic, Pa, Pb, HF_exchange_frac)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force
implicit none
integer,intent(in) :: nConts, natoms, iatom, ic
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: dF_full(:,:), dS_full(:,:)
integer :: mu,i,nu

allocate(dS(nConts,nConts,3), dHcore(nConts,nConts,3))
allocate(dJi(nConts,nConts,3), dKa(nConts,nConts,3), dKb(nConts,nConts,3))
call integrals_compute_force(nConts, Pa, Pb, dS, dHcore, dJi, dKa, dKb)
allocate(ao_atom(nConts))
mu = 0
do i = 1,natoms
   do nu = 1,atoms(i)%nconts
      mu = mu+1
      ao_atom(mu) = i
   enddo
enddo
allocate(dF_full(nConts,nConts), dS_full(nConts,nConts))
call cphf_fock_deriv_atom(nConts, ao_atom, iatom, ic, dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, dF_full, dS_full)
print *, 'CPHF DEBUG dF_full raw row1=', dF_full(1,:)
print *, 'CPHF DEBUG dF_full raw row2=', dF_full(2,:)
print *, 'CPHF DEBUG RAW dS (1,6,ic)=', dS(1,6,ic), ' dS(6,1,ic)=', dS(6,1,ic)
print *, 'CPHF DEBUG RAW pieces (1,1,ic)=', dHcore(1,1,ic), dJi(1,1,ic), dKa(1,1,ic)
print *, 'CPHF DEBUG RAW pieces (1,6,ic)=', dHcore(1,6,ic), dJi(1,6,ic), dKa(1,6,ic)
print *, 'CPHF DEBUG RAW pieces (6,1,ic)=', dHcore(6,1,ic), dJi(6,1,ic), dKa(6,1,ic)
print *, 'CPHF DEBUG RAW pieces (1,9,ic)=', dHcore(1,9,ic), dJi(1,9,ic), dKa(1,9,ic)
print *, 'CPHF DEBUG RAW pieces (9,1,ic)=', dHcore(9,1,ic), dJi(9,1,ic), dKa(9,1,ic)
deallocate(dS,dHcore,dJi,dKa,dKb,ao_atom,dF_full,dS_full)
end subroutine cphf_debug_dF_raw

end submodule cpscf_debug_impl
