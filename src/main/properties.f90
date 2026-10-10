! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Post-SCF property evaluation (population analysis, dipole, orbital printing).

submodule (mod_integrals) properties_impl
implicit none
contains

module subroutine calc_mayer_bond_orders(Natoms_in)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: Natoms_in
integer,allocatable :: atom_ao_off(:)
real(8),allocatable :: PaS(:,:), PbS(:,:)
real(8) :: bsum
integer :: ia, ib, i0, j0

allocate(atom_ao_off(0:Natoms_in))
atom_ao_off(0) = 0
do ia = 1,Natoms_in
   atom_ao_off(ia) = atom_ao_off(ia-1) + atoms(ia)%nconts
enddo
allocate(PaS(nConts,nConts), PbS(nConts,nConts))
PaS = matmul(Pa,S)
PbS = matmul(Pb,S)
if (allocated(mayer_bo)) deallocate(mayer_bo)
allocate(mayer_bo(Natoms_in,Natoms_in))
mayer_bo = 0.0d0
do ia = 1,Natoms_in
   do ib = ia+1,Natoms_in
      bsum = 0.0d0
      do j0 = atom_ao_off(ib-1)+1,atom_ao_off(ib)
         do i0 = atom_ao_off(ia-1)+1,atom_ao_off(ia)
            bsum = bsum + PaS(i0,j0)*PaS(j0,i0) + PbS(i0,j0)*PbS(j0,i0)
         enddo
      enddo
      mayer_bo(ia,ib) = 2.0d0*bsum
      mayer_bo(ib,ia) = 2.0d0*bsum
   enddo
enddo
deallocate(PaS, PbS, atom_ao_off)
end subroutine calc_mayer_bond_orders

module subroutine calc_properties(Natoms_in, MLcharge_out, RESPcharge_in)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: Natoms_in
real(8),intent(out) :: MLcharge_out(Natoms_in)
real(8),intent(in),optional :: RESPcharge_in(Natoms_in)

integer,allocatable :: ao_offset(:), atom_ao_off(:)
real(8),allocatable :: Rx(:,:), Ry(:,:), Rz(:,:)
real(8),allocatable :: Ptot(:,:), PS(:,:), Shalf(:,:), PSh(:,:)
real(8),allocatable :: buf(:,:,:)
real(8),allocatable :: q_lowdin(:)
real(8) :: dip_nuc(3), dip_elec(3), dip_au(3)
real(8) :: dip_au_norm, dip_debye_norm
real(8),parameter :: AU2DEBYE = 2.541746d0
integer :: ia, i0, j0, di, dj, n, num, e1, e2, mo
integer :: shls(4)
real(8) :: Zeff

real(8),allocatable :: S_e(:), S_tmp(:,:), S_diag(:,:), WORKs(:)
integer :: LWORKs, INFOs

allocate(ao_offset(0:nBases-1))
num = 0
do i0 = 0,nBases-1
   ao_offset(i0) = num
   num = num + cgto_engine(i0,bas)
enddo

allocate(atom_ao_off(0:Natoms_in))
atom_ao_off(0) = 0
do ia = 1,Natoms_in
   atom_ao_off(ia) = atom_ao_off(ia-1) + atoms(ia)%nconts
enddo

n = 0
do i0 = 0,nBases-1
   n = max(n, cgto_engine(i0,bas))
enddo
allocate(buf(n,n,3))
allocate(Rx(nConts,nConts),Ry(nConts,nConts),Rz(nConts,nConts))
Rx = 0.0d0; Ry = 0.0d0; Rz = 0.0d0
do i0 = 0,nBases-1
   do j0 = i0,nBases-1
      shls(1) = i0; shls(2) = j0
      di = cgto_engine(i0,bas); dj = cgto_engine(j0,bas)
      call r1e_engine(buf(1:di,1:dj,1:3), shls, atm, size(atm,2), bas, nBases, env, 0_8)
      do e1 = 1,di
         do e2 = 1,dj
            Rx(ao_offset(i0)+e1,ao_offset(j0)+e2) = buf(e1,e2,1)*NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
            Ry(ao_offset(i0)+e1,ao_offset(j0)+e2) = buf(e1,e2,2)*NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
            Rz(ao_offset(i0)+e1,ao_offset(j0)+e2) = buf(e1,e2,3)*NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
            Rx(ao_offset(j0)+e2,ao_offset(i0)+e1) = Rx(ao_offset(i0)+e1,ao_offset(j0)+e2)
            Ry(ao_offset(j0)+e2,ao_offset(i0)+e1) = Ry(ao_offset(i0)+e1,ao_offset(j0)+e2)
            Rz(ao_offset(j0)+e2,ao_offset(i0)+e1) = Rz(ao_offset(i0)+e1,ao_offset(j0)+e2)
         enddo
      enddo
   enddo
enddo
deallocate(buf)

allocate(Ptot(nConts,nConts))
Ptot = Pa + Pb

dip_nuc = 0.0d0
do ia = 1,Natoms_in
   Zeff = real(atoms(ia)%charge - atoms(ia)%ecpCoreElec,8)
   dip_nuc = dip_nuc + Zeff*atoms(ia)%coor*ans2bohr
enddo
dip_elec(1) = -sum(Ptot*Rx)
dip_elec(2) = -sum(Ptot*Ry)
dip_elec(3) = -sum(Ptot*Rz)
dip_au = dip_nuc + dip_elec
dip_au_norm = sqrt(sum(dip_au**2))
dip_debye_norm = dip_au_norm*AU2DEBYE
deallocate(Rx,Ry,Rz)

allocate(PS(nConts,nConts))
PS = matmul(Ptot,S)
do ia = 1,Natoms_in
   Zeff = real(atoms(ia)%charge - atoms(ia)%ecpCoreElec,8)
   MLcharge_out(ia) = Zeff - sum((/(PS(i0,i0),i0=atom_ao_off(ia-1)+1,atom_ao_off(ia))/))
enddo
deallocate(PS)

allocate(S_e(nConts),S_tmp(nConts,nConts),S_diag(nConts,nConts))
S_tmp = S
LWORKs = 1 + 6*nConts + nConts**2
allocate(WORKs(LWORKs))
call DSYEV('V','U',nConts,S_tmp,nConts,S_e,WORKs,LWORKs,INFOs)
S_diag = 0.0d0
do i0 = 1,nConts
   S_diag(i0,i0) = sqrt(max(S_e(i0),0.0d0))
enddo
allocate(Shalf(nConts,nConts))
Shalf = matmul(S_tmp,matmul(S_diag,transpose(S_tmp)))
deallocate(S_e,S_tmp,S_diag,WORKs)

allocate(PSh(nConts,nConts),q_lowdin(Natoms_in))
PSh = matmul(Shalf,matmul(Ptot,Shalf))
do ia = 1,Natoms_in
   Zeff = real(atoms(ia)%charge - atoms(ia)%ecpCoreElec,8)
   q_lowdin(ia) = Zeff - sum((/(PSh(i0,i0),i0=atom_ao_off(ia-1)+1,atom_ao_off(ia))/))
enddo
deallocate(Shalf,PSh,Ptot)

print *, '[RESULTS]'
print *,"===== Molecular Properties ====="
print '(" Dipole moment (a.u.):   X=",F10.5,"  Y=",F10.5,"  Z=",F10.5,"  |mu|=",F10.5)', &
      dip_au(1),dip_au(2),dip_au(3),dip_au_norm
print '(" Dipole moment (Debye):  X=",F10.5,"  Y=",F10.5,"  Z=",F10.5,"  |mu|=",F10.5)', &
      dip_au(1)*AU2DEBYE,dip_au(2)*AU2DEBYE,dip_au(3)*AU2DEBYE,dip_debye_norm
print *
if (present(RESPcharge_in)) then
   print '(A6,A5,3A12)', 'Atom','Z','Mulliken_q','Lowdin_q','RESP_q'
   do ia = 1,Natoms_in
      print '(I6,I5,3F12.6)', ia, atoms(ia)%charge, MLcharge_out(ia), q_lowdin(ia), RESPcharge_in(ia)
   enddo
else
   print '(A6,A5,2A12)', 'Atom','Z','Mulliken_q','Lowdin_q'
   do ia = 1,Natoms_in
      print '(I6,I5,2F12.6)', ia, atoms(ia)%charge, MLcharge_out(ia), q_lowdin(ia)
   enddo
endif
deallocate(q_lowdin)
print *
print *,"  Orbital energies (Hartree):"
if (multi .eq. 1) then
   call print_orbital_block("  Occupied (alpha=beta)", eLev_a, 1, n_alpha)
   call print_orbital_block("  Virtual  (alpha=beta)", eLev_a, n_alpha+1, nConts)
else
   call print_orbital_block("  Occupied alpha", eLev_a, 1, n_alpha)
   call print_orbital_block("  Virtual  alpha", eLev_a, n_alpha+1, nConts)
   call print_orbital_block("  Occupied beta ", eLev_b, 1, n_beta)
   call print_orbital_block("  Virtual  beta ", eLev_b, n_beta+1, nConts)
endif
print *,"================================="
print *

deallocate(ao_offset,atom_ao_off)

end subroutine calc_properties

subroutine print_orbital_block(label, eLev, lo, hi)
implicit none
character(len=*),intent(in) :: label
real(8),intent(in) :: eLev(:)
integer,intent(in) :: lo, hi
integer,parameter :: PER_ROW = 5
integer :: mo, row_end

if (hi .lt. lo) return
print '(A,":")', trim(label)
do mo = lo,hi,PER_ROW
   row_end = min(mo+PER_ROW-1,hi)
   print '(4X,5F13.6)', eLev(mo:row_end)
enddo
end subroutine print_orbital_block

module subroutine integrals_dipole_ao(nbf, R)
implicit none
integer,intent(in) :: nbf
real(8),intent(out) :: R(nbf,nbf,3)
integer,allocatable :: ao_offset(:)
real(8),allocatable :: buf(:,:,:)
integer :: i0, j0, di, dj, n, num, e1, e2, c, shls(4)
allocate(ao_offset(0:nBases-1))
num = 0
n = 0
do i0 = 0,nBases-1
   ao_offset(i0) = num
   num = num + cgto_engine(i0,bas)
   n = max(n, cgto_engine(i0,bas))
enddo
allocate(buf(n,n,3))
R = 0.0d0
do i0 = 0,nBases-1
   do j0 = 0,nBases-1
      shls(1) = i0; shls(2) = j0
      di = cgto_engine(i0,bas); dj = cgto_engine(j0,bas)
      call r1e_engine(buf(1:di,1:dj,1:3), shls, atm, size(atm,2), bas, nBases, env, 0_8)
      do c = 1,3
         do e2 = 1,dj
            do e1 = 1,di
               R(ao_offset(i0)+e1,ao_offset(j0)+e2,c) = buf(e1,e2,c)*NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
            enddo
         enddo
      enddo
   enddo
enddo
end subroutine integrals_dipole_ao

module subroutine integrals_dipole_pulay(nbf, natoms, ao_atom, Ptot, dpul)
implicit none
integer,intent(in) :: nbf, natoms, ao_atom(nbf)
real(8),intent(in) :: Ptot(nbf,nbf)
real(8),intent(out) :: dpul(3,3,natoms)
integer,allocatable :: ao_offset(:)
real(8),allocatable :: buf(:,:,:)
real(8) :: acc(3,3)
integer :: i0, j0, di, dj, n, num, e1, e2, k, c, ix, a, shls(4)
allocate(ao_offset(0:nBases-1))
num = 0
n = 0
do i0 = 0,nBases-1
   ao_offset(i0) = num
   num = num + cgto_engine(i0,bas)
   n = max(n, cgto_engine(i0,bas))
enddo
dpul = 0.0d0
!$omp parallel private(i0,j0,di,dj,buf,shls,e1,e2,c,ix,a,acc,k) reduction(+:dpul)
allocate(buf(n,n,9))
!$omp do schedule(dynamic)
do j0 = 0,nBases-1
   do i0 = 0,nBases-1
      shls(1) = i0; shls(2) = j0
      di = cgto_engine(i0,bas); dj = cgto_engine(j0,bas)
      call r1e_irp_engine(buf(1:di,1:dj,1:9), shls, atm, size(atm,2), bas, nBases, env, 0_8)
      a = ao_atom(ao_offset(j0)+1)
      acc = 0.0d0
      do c = 1,3
         do ix = 1,3
            k = 3*(c-1)+ix
            do e2 = 1,dj
               do e1 = 1,di
                  acc(c,ix) = acc(c,ix) + Ptot(ao_offset(i0)+e1,ao_offset(j0)+e2) * buf(e1,e2,k) &
                             * NorVEC(ao_offset(i0)+e1)*NorVEC(ao_offset(j0)+e2)
               enddo
            enddo
         enddo
      enddo
      dpul(:,:,a) = dpul(:,:,a) + 2.0d0*acc
   enddo
enddo
!$omp end do
deallocate(buf)
!$omp end parallel
end subroutine integrals_dipole_pulay

end submodule properties_impl
