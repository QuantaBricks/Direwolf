! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! One-electron piece of the analytic Hessian (S + kinetic + nuclear

subroutine calc_hessian_1e(nConts, natoms, ao_atom, Ptot, Wtot, znum, ecpce, coor_bohr, d2E_1e)
use mod_integrals, only: integrals_hessian_1e_shell, integrals_hessian_nuc_atom, integrals_hessian_nuc_atom_raw
implicit none
integer,intent(in) :: nConts, natoms
integer,intent(in) :: ao_atom(nConts)
real(8),intent(in) :: Ptot(nConts,nConts), Wtot(nConts,nConts)
integer,intent(in) :: znum(natoms), ecpce(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
real(8),intent(out) :: d2E_1e(3,natoms,3,natoms)

real(8),allocatable :: s1aa(:,:,:), s1ab(:,:,:), h1aa(:,:,:), h1ab(:,:,:)
real(8),allocatable :: CC(:,:,:), CB(:,:,:), CK(:,:,:)
integer :: mu,nu,ia,ja,ic
real(8) :: blk(3,3), Zeff

allocate(s1aa(nConts,nConts,9), s1ab(nConts,nConts,9))
allocate(h1aa(nConts,nConts,9), h1ab(nConts,nConts,9))
call integrals_hessian_1e_shell(nConts, s1aa, s1ab, h1aa, h1ab)

d2E_1e = 0.0d0

do mu = 1,nConts
   ia = ao_atom(mu)
   do nu = 1,nConts
      ja = ao_atom(nu)
      blk = transpose(reshape(h1aa(mu,nu,:)*Ptot(mu,nu) - s1aa(mu,nu,:)*Wtot(mu,nu), (/3,3/)))
      d2E_1e(:,ia,:,ia) = d2E_1e(:,ia,:,ia) + blk
      blk = transpose(reshape(h1ab(mu,nu,:)*Ptot(mu,nu) - s1ab(mu,nu,:)*Wtot(mu,nu), (/3,3/)))
      d2E_1e(:,ia,:,ja) = d2E_1e(:,ia,:,ja) + blk
      d2E_1e(:,ja,:,ia) = d2E_1e(:,ja,:,ia) + transpose(blk)
   enddo
enddo
deallocate(s1aa, s1ab, h1aa, h1ab)

allocate(CC(nConts,nConts,9), CB(nConts,nConts,9), CK(nConts,nConts,9))
do ic = 1,natoms
   Zeff = real(znum(ic)-ecpce(ic),8)
   call integrals_hessian_nuc_atom(coor_bohr(:,ic), nConts, CC, CB, CK)
   do mu = 1,nConts
      ia = ao_atom(mu)
      do nu = 1,nConts
         ja = ao_atom(nu)
         d2E_1e(:,ic,:,ic) = d2E_1e(:,ic,:,ic) &
              - Zeff*transpose(reshape(CC(mu,nu,:)*Ptot(mu,nu),(/3,3/)))
         blk = -Zeff*reshape(CB(mu,nu,:)*Ptot(mu,nu),(/3,3/))
         d2E_1e(:,ic,:,ia) = d2E_1e(:,ic,:,ia) + blk
         d2E_1e(:,ia,:,ic) = d2E_1e(:,ia,:,ic) + transpose(blk)
         blk = -Zeff*transpose(reshape(CK(mu,nu,:)*Ptot(mu,nu),(/3,3/)))
         d2E_1e(:,ic,:,ja) = d2E_1e(:,ic,:,ja) + blk
         d2E_1e(:,ja,:,ic) = d2E_1e(:,ja,:,ic) + transpose(blk)
      enddo
   enddo
enddo
deallocate(CC, CB, CK)

block
   real(8),allocatable :: h1aa2(:,:,:), h1ab2(:,:,:), s1aa2(:,:,:), s1ab2(:,:,:)
   real(8),allocatable :: AA(:,:,:), AB(:,:,:)
   real(8) :: d2diag(3,3), hblk(9), sblk(9)
   integer :: p,q
   allocate(h1aa2(nConts,nConts,9), h1ab2(nConts,nConts,9))
   allocate(s1aa2(nConts,nConts,9), s1ab2(nConts,nConts,9))
   call integrals_hessian_1e_shell(nConts, s1aa2, s1ab2, h1aa2, h1ab2)
   allocate(AA(nConts,nConts,9), AB(nConts,nConts,9))
   do ic = 1,natoms
      Zeff = real(znum(ic)-ecpce(ic),8)
      call integrals_hessian_nuc_atom_raw(coor_bohr(:,ic), nConts, AA, AB)
      d2diag = 0.0d0
      block
         real(8) :: d2diag_h(3,3), d2diag_s(3,3)
         d2diag_h = 0.0d0; d2diag_s = 0.0d0
         do p = 1,nConts
            do q = 1,nConts
               hblk = -Zeff*AA(p,q,:) - Zeff*AB(p,q,:)
               sblk = 0.0d0
               if (ao_atom(p) == ic) then
                  hblk = hblk + h1aa2(p,q,:) + Zeff*AA(p,q,:) + Zeff*AB(p,q,:)
                  sblk = sblk + s1aa2(p,q,:)
               endif
               if (ao_atom(q) == ic) then
                  hblk = hblk + Zeff*AA(q,p,:) + Zeff*AB(p,q,:)
               endif
               if (ao_atom(p) == ic .and. ao_atom(q) == ic) then
                  hblk = hblk + h1ab2(p,q,:)
                  sblk = sblk + s1ab2(p,q,:)
               endif
               d2diag_h = d2diag_h + transpose(reshape(hblk*Ptot(p,q),(/3,3/)))
               d2diag_s = d2diag_s + transpose(reshape(sblk*Wtot(p,q),(/3,3/)))
            enddo
         enddo
         d2diag = 2.0d0*d2diag_h - 2.0d0*d2diag_s
      end block
      d2E_1e(:,ic,:,ic) = d2diag
   enddo
   deallocate(h1aa2,h1ab2,s1aa2,s1ab2,AA,AB)
end block

end subroutine calc_hessian_1e
