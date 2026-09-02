! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cphf: closed-shell (RHF) coupled-perturbed Hartree-Fock solver -

module mod_cphf
implicit none
private
public :: cphf_solve_rhf, cphf_debug_check_residual, cphf_debug_dump_dP, cphf_debug_grad_check
public :: cphf_hess_response, cphf_hess_nuc, cphf_debug_hessian_row, cphf_debug_e1_fixed
public :: cphf_debug_ccbk_raw, cphf_debug_dF_raw

contains

subroutine cphf_fock_deriv_atom(nConts, ao_atom, iatom, ic, &
                                 dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, &
                                 dF_full, dS_full)
use MOL_info, only: atoms
use mod_integrals, only: integrals_nuc_attraction_deriv, integrals_cphf_2e_density_role
implicit none
integer,intent(in) :: nConts, ao_atom(nConts), iatom, ic
real(8),intent(in) :: dS(nConts,nConts,3), dHcore(nConts,nConts,3)
real(8),intent(in) :: dJi(nConts,nConts,3), dKa(nConts,nConts,3)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(out) :: dF_full(nConts,nConts), dS_full(nConts,nConts)
INCLUDE 'parameter.h'

integer :: mu,nu
real(8),allocatable :: Dr(:,:,:)
real(8),allocatable :: dJ2(:,:,:), dK2(:,:,:)
real(8) :: Zeff

dF_full = 0.0d0
dS_full = 0.0d0
do nu = 1,nConts
   do mu = 1,nConts
      if (ao_atom(mu) == iatom) then
         dF_full(mu,nu) = dF_full(mu,nu) - (dHcore(mu,nu,ic) + dJi(mu,nu,ic) - HF_exchange_frac*dKa(mu,nu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(mu,nu,ic)
      endif
      if (ao_atom(nu) == iatom) then
         dF_full(mu,nu) = dF_full(mu,nu) - (dHcore(nu,mu,ic) + dJi(nu,mu,ic) - HF_exchange_frac*dKa(nu,mu,ic))
         dS_full(mu,nu) = dS_full(mu,nu) - dS(nu,mu,ic)
      endif
   enddo
enddo

allocate(Dr(nConts,nConts,3))
call integrals_nuc_attraction_deriv(atoms(iatom)%coor*ans2bohr, nConts, Dr)
Zeff = atoms(iatom)%charge - atoms(iatom)%ecpCoreElec
dF_full = dF_full - Zeff*(Dr(:,:,ic) + transpose(Dr(:,:,ic)))
deallocate(Dr)

allocate(dJ2(nConts,nConts,3), dK2(nConts,nConts,3))
call integrals_cphf_2e_density_role(nConts, iatom, Pa, Pb, dJ2, dK2)
do nu = 1,nConts
   do mu = 1,nConts
      dF_full(mu,nu) = dF_full(mu,nu) &
         - (dJ2(mu,nu,ic) - HF_exchange_frac*dK2(mu,nu,ic)) &
         - (dJ2(nu,mu,ic) - HF_exchange_frac*dK2(nu,mu,ic))
   enddo
enddo
deallocate(dJ2, dK2)
end subroutine cphf_fock_deriv_atom

subroutine cphf_build_amat(nConts, nocc, C_a, HF_exchange_frac, Amat)
use mod_exchange, only: exchange_build
use mod_integrals, only: integrals_build_coulomb
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: Amat(nConts-nocc,nocc,nConts-nocc,nocc)

integer :: nvir
integer :: i,j,a,b,mu,nu
real(8),allocatable :: Dtrial(:,:), Jtr(:,:), Katr(:,:), Kbtr(:,:), Rtr(:,:), Tmp(:,:)

nvir = nConts - nocc
allocate(Dtrial(nConts,nConts), Jtr(nConts,nConts), Katr(nConts,nConts), Kbtr(nConts,nConts))
allocate(Rtr(nConts,nConts), Tmp(nConts,nocc))

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
            Amat(a,i,b,j) = dot_product(C_a(:,nocc+a), Tmp(:,i))
         enddo
      enddo
   enddo
enddo
deallocate(Dtrial,Jtr,Katr,Kbtr,Rtr,Tmp)
end subroutine cphf_build_amat

subroutine cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, &
                          U_full, h1_MO_all, s1_MO_all, mo_e1_all)
use MOL_info, only: atoms
use mod_integrals, only: integrals_compute_force
use mod_exchange, only: exchange_build
use mod_integrals, only: integrals_build_coulomb
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: U_full(nConts,nocc,3,natoms)
real(8),intent(out) :: h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms)
real(8),intent(out) :: mo_e1_all(nocc,nocc,3,natoms)

integer :: nvir
real(8),allocatable :: dS(:,:,:), dHcore(:,:,:), dJi(:,:,:), dKa(:,:,:), dKb(:,:,:)
integer,allocatable :: ao_atom(:)
real(8),allocatable :: Amat(:,:,:,:)
real(8),allocatable :: dF_full(:,:), dS_full(:,:), h1_MO(:,:), s1_MO(:,:)
real(8),allocatable :: s1oo(:,:), Pconn(:,:), Jc(:,:), Kac(:,:), Kbc(:,:), Rc(:,:), Tmp(:,:)
real(8),allocatable :: Uvo(:,:), Pvo(:,:), Ptot2(:,:), v1oo(:,:)
integer :: iat,ic,i,j,a,b,mu

nvir = nConts - nocc

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

allocate(Amat(nvir,nocc,nvir,nocc))
call cphf_build_amat(nConts, nocc, C_a, HF_exchange_frac, Amat)

allocate(dF_full(nConts,nConts), dS_full(nConts,nConts))
allocate(h1_MO(nConts,nConts), s1_MO(nConts,nConts))
allocate(s1oo(nocc,nocc), Pconn(nConts,nConts))
allocate(Jc(nConts,nConts), Kac(nConts,nConts), Kbc(nConts,nConts), Rc(nConts,nConts), Tmp(nConts,nocc))
allocate(Uvo(nvir,nocc), Pvo(nConts,nConts), Ptot2(nConts,nConts), v1oo(nocc,nocc))

U_full = 0.0d0

do iat = 1,natoms
   do ic = 1,3
      call cphf_fock_deriv_atom(nConts, ao_atom, iat, ic, dS, dHcore, dJi, dKa, HF_exchange_frac, Pa, Pb, dF_full, dS_full)
      h1_MO = matmul(transpose(C_a), matmul(dF_full, C_a))
      s1_MO = matmul(transpose(C_a), matmul(dS_full, C_a))
      h1_MO_all(:,:,ic,iat) = h1_MO(:,1:nocc)
      s1_MO_all(:,:,ic,iat) = s1_MO(:,1:nocc)

      s1oo = s1_MO(1:nocc,1:nocc)
      Pconn = -2.0d0*matmul(C_a(:,1:nocc), matmul(s1oo, transpose(C_a(:,1:nocc))))
      call integrals_build_coulomb(nConts, Pconn, Jc)
      call exchange_build(nConts, Pconn, Pconn, Kac, Kbc)
      Rc = Jc - 0.5d0*HF_exchange_frac*Kac
      Tmp = matmul(Rc, C_a(:,1:nocc))

      Uvo = 0.0d0
      block
         real(8) :: Uvo_new(nvir,nocc), Amat_dot, diffmax
         integer :: iterCPHF
         do iterCPHF = 1,200
            do i = 1,nocc
               do a = 1,nvir
                  Amat_dot = 0.0d0
                  do j = 1,nocc
                     do b = 1,nvir
                        Amat_dot = Amat_dot + Amat(a,i,b,j)*Uvo(b,j)
                     enddo
                  enddo
                  Uvo_new(a,i) = ( h1_MO(nocc+a,i) - s1_MO(nocc+a,i)*eLev_a(i) &
                                  + dot_product(C_a(:,nocc+a), Tmp(:,i)) + Amat_dot ) &
                                 / (eLev_a(i)-eLev_a(nocc+a))
               enddo
            enddo
            diffmax = maxval(abs(Uvo_new-Uvo))
            Uvo = Uvo_new
            if (diffmax < 1.0d-12) exit
         enddo
      end block

      do i = 1,nocc
         do j = 1,nocc
            U_full(j,i,ic,iat) = -0.5d0*s1oo(j,i)
         enddo
         do a = 1,nvir
            U_full(nocc+a,i,ic,iat) = Uvo(a,i)
         enddo
      enddo

      Pvo = 2.0d0*( matmul(C_a(:,nocc+1:nConts), matmul(Uvo, transpose(C_a(:,1:nocc)))) &
                  + matmul(C_a(:,1:nocc), matmul(transpose(Uvo), transpose(C_a(:,nocc+1:nConts)))) )
      Ptot2 = Pconn + Pvo
      call integrals_build_coulomb(nConts, Ptot2, Jc)
      call exchange_build(nConts, Ptot2, Ptot2, Kac, Kbc)
      Rc = Jc - 0.5d0*HF_exchange_frac*Kac
      v1oo = matmul(transpose(C_a(:,1:nocc)), matmul(Rc, C_a(:,1:nocc)))
      do i = 1,nocc
         do j = 1,nocc
            mo_e1_all(i,j,ic,iat) = h1_MO(i,j) - 0.5d0*(eLev_a(i)+eLev_a(j))*s1_MO(i,j) + v1oo(i,j)
         enddo
      enddo
   enddo
enddo

deallocate(dS,dHcore,dJi,dKa,dKb,ao_atom,Amat)
deallocate(dF_full,dS_full,h1_MO,s1_MO,s1oo,Pconn,Jc,Kac,Kbc,Rc,Tmp)
deallocate(Uvo,Pvo,Ptot2,v1oo)
end subroutine cphf_solve_rhf

subroutine cphf_debug_check_residual(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
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

subroutine cphf_debug_dump_dP(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
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

subroutine cphf_debug_grad_check(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
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

subroutine cphf_hess_response(nConts, natoms, nocc, U_full, h1_MO_all, s1_MO_all, mo_e1_all, eLev_a, d2E_response)
implicit none
integer,intent(in) :: nConts, natoms, nocc
real(8),intent(in) :: U_full(nConts,nocc,3,natoms)
real(8),intent(in) :: h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms)
real(8),intent(in) :: mo_e1_all(nocc,nocc,3,natoms)
real(8),intent(in) :: eLev_a(nConts)
real(8),intent(out) :: d2E_response(3,natoms,3,natoms)

integer :: ia,ja,x,y,i
real(8),allocatable :: hs(:,:)
real(8) :: t1,t3

allocate(hs(nConts,nocc))
do ia = 1,natoms
   do x = 1,3
      do i = 1,nocc
         hs(:,i) = h1_MO_all(:,i,x,ia) - eLev_a(i)*s1_MO_all(:,i,x,ia)
      enddo
      do ja = 1,natoms
         do y = 1,3
            t1 = sum(hs*U_full(:,:,y,ja))
            t3 = sum(s1_MO_all(1:nocc,:,x,ia)*mo_e1_all(:,:,y,ja))
            d2E_response(x,ia,y,ja) = 4.0d0*t1 - 2.0d0*t3
         enddo
      enddo
   enddo
enddo
deallocate(hs)
end subroutine cphf_hess_response

subroutine cphf_hess_nuc(natoms, znum, ecpce, coor_bohr, d2Enn)
implicit none
integer,intent(in) :: natoms, znum(natoms), ecpce(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(out) :: d2Enn(3,natoms,3,natoms)

integer :: i,k,x,y
real(8) :: r(3), r1, qi, qk, blk(3,3)

d2Enn = 0.0d0
do i = 1,natoms
   qi = real(znum(i)-ecpce(i),8)
   do k = 1,natoms
      if (k == i) cycle
      qk = real(znum(k)-ecpce(k),8)
      r = coor_bohr(:,i)-coor_bohr(:,k)
      r1 = sqrt(sum(r**2))
      do y = 1,3
         do x = 1,3
            blk(x,y) = qi*qk*3.0d0*r(x)*r(y)/r1**5
         enddo
         blk(y,y) = blk(y,y) - qi*qk/r1**3
      enddo
      d2Enn(:,i,:,i) = d2Enn(:,i,:,i) + blk
      d2Enn(:,i,:,k) = d2Enn(:,i,:,k) - blk
   enddo
enddo
end subroutine cphf_hess_nuc

subroutine cphf_debug_hessian_row(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
use MOL_info, only: atoms
use mod_integrals, only: integrals_hessian_1e_shell, integrals_hessian_2e_direct
implicit none
integer,intent(in) :: nConts, natoms, nocc, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

real(8),allocatable :: U_full(:,:,:,:), h1_MO_all(:,:,:,:), s1_MO_all(:,:,:,:), mo_e1_all(:,:,:,:)
real(8),allocatable :: d2E_1e(:,:,:,:), d2E_2e(:,:,:,:), d2E_resp(:,:,:,:), d2Enn(:,:,:,:), d2E_tot(:,:,:,:)
real(8),allocatable :: Ptot(:,:), Wtot(:,:)
integer,allocatable :: ao_atom(:), znum(:), ecpce(:)
real(8),allocatable :: coor_bohr(:,:)
INCLUDE 'parameter.h'
integer :: mu,i,ja,y

allocate(U_full(nConts,nocc,3,natoms))
allocate(h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms), mo_e1_all(nocc,nocc,3,natoms))
call cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, U_full, h1_MO_all, s1_MO_all, mo_e1_all)

allocate(d2E_resp(3,natoms,3,natoms))
call cphf_hess_response(nConts, natoms, nocc, U_full, h1_MO_all, s1_MO_all, mo_e1_all, eLev_a, d2E_resp)

allocate(ao_atom(nConts))
mu = 0
do i = 1,natoms
   do y = 1,atoms(i)%nconts
      mu = mu+1
      ao_atom(mu) = i
   enddo
enddo

allocate(Ptot(nConts,nConts), Wtot(nConts,nConts))
Ptot = Pa+Pb
Wtot = 0.0d0
do mu = 1,nConts
   do y = 1,nConts
      do i = 1,nocc
         Wtot(mu,y) = Wtot(mu,y) + 2.0d0*eLev_a(i)*C_a(mu,i)*C_a(y,i)
      enddo
   enddo
enddo

allocate(znum(natoms), ecpce(natoms), coor_bohr(3,natoms))
do i = 1,natoms
   znum(i) = atoms(i)%charge
   ecpce(i) = atoms(i)%ecpCoreElec
   coor_bohr(:,i) = atoms(i)%coor*ans2bohr
enddo

allocate(d2E_1e(3,natoms,3,natoms))
call calc_hessian_1e(nConts, natoms, ao_atom, Ptot, Wtot, znum, ecpce, coor_bohr, d2E_1e)

allocate(d2E_2e(3,natoms,3,natoms))
call integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e)

allocate(d2Enn(3,natoms,3,natoms))
call cphf_hess_nuc(natoms, znum, ecpce, coor_bohr, d2Enn)

allocate(d2E_tot(3,natoms,3,natoms))
d2E_tot = d2E_1e + d2E_2e + d2E_resp + d2Enn

print *, 'CPHF DEBUG HESSIAN ROW: perturbation atom=',iatom,' comp=',ic
do ja = 1,natoms
   do y = 1,3
      print *, 'target atom=',ja,' comp=',y,' d2E=', d2E_tot(ic,iatom,y,ja), &
               ' (1e=',d2E_1e(ic,iatom,y,ja),' 2e=',d2E_2e(ic,iatom,y,ja), &
               ' resp=',d2E_resp(ic,iatom,y,ja),' nn=',d2Enn(ic,iatom,y,ja),')'
   enddo
enddo
print *, 'CPHF DEBUG HESSIAN symmetry check: max|d2E(x,a,y,b)-d2E(y,b,x,a)| =', &
   maxval(abs(d2E_tot - reshape(d2E_tot,(/3,natoms,3,natoms/),order=(/3,4,1,2/))))

deallocate(U_full,h1_MO_all,s1_MO_all,mo_e1_all,d2E_resp,ao_atom,Ptot,Wtot,znum,ecpce,coor_bohr)
deallocate(d2E_1e,d2E_2e,d2Enn,d2E_tot)
end subroutine cphf_debug_hessian_row

subroutine cphf_debug_e1_fixed(nConts, Pa, Pb)
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

subroutine cphf_debug_ccbk_raw(nConts, natoms, iatom_c, mu, nu)
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

subroutine cphf_debug_dF_raw(nConts, natoms, iatom, ic, Pa, Pb, HF_exchange_frac)
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

end module mod_cphf
