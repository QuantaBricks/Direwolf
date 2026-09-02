! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_localization: orbital localization, decoupled from the integral

module mod_localization
implicit none
private
public :: pipek_mezey_localize

integer,parameter :: PM_MAX_SWEEPS = 100
real(8),parameter :: PM_CONV_TOL = 1.0d-4

contains

subroutine pipek_mezey_localize(nConts, S, C, nOcc, C_loc, n_sweeps_used)
use MOL_info, only: atoms, nAtoms
implicit none
integer,intent(in) :: nConts, nOcc
real(8),intent(in) :: S(nConts,nConts), C(nConts,nOcc)
real(8),intent(out) :: C_loc(nConts,nOcc)
integer,intent(out),optional :: n_sweeps_used
real(8),allocatable :: SC(:,:), tmp1(:), tmp2(:)
integer,allocatable :: atom_start(:), atom_count(:)
integer :: iat, i, j, sweep, ipos
real(8) :: Aij, Bij, gamma, cg, sg, Qi, Qj, Qij
real(8) :: func_old, func_new

allocate(atom_start(nAtoms), atom_count(nAtoms))
ipos = 0
do iat = 1,nAtoms
   atom_start(iat) = ipos+1
   atom_count(iat) = atoms(iat)%nconts
   ipos = ipos + atoms(iat)%nconts
enddo

C_loc = C
allocate(SC(nConts,nOcc))
call dgemm('N','N', nConts, nOcc, nConts, 1.0d0, S, nConts, C_loc, nConts, 0.0d0, SC, nConts)

allocate(tmp1(nConts), tmp2(nConts))
func_old = pm_functional(nConts, nOcc, nAtoms, atom_start, atom_count, C_loc, SC)

do sweep = 1,PM_MAX_SWEEPS
   do i = 1,nOcc-1
      do j = i+1,nOcc
         Aij = 0.0d0
         Bij = 0.0d0
         do iat = 1,nAtoms
            if (atom_count(iat) .eq. 0) cycle
            associate(a0 => atom_start(iat), a1 => atom_start(iat)+atom_count(iat)-1)
               Qi  = dot_product(C_loc(a0:a1,i), SC(a0:a1,i))
               Qj  = dot_product(C_loc(a0:a1,j), SC(a0:a1,j))
               Qij = 0.5d0*(dot_product(C_loc(a0:a1,i), SC(a0:a1,j)) &
                            +dot_product(C_loc(a0:a1,j), SC(a0:a1,i)))
            end associate
            Aij = Aij + Qij*Qij - 0.25d0*(Qi-Qj)**2
            Bij = Bij + Qij*(Qi-Qj)
         enddo
         if (abs(Aij) .lt. 1.0d-14 .and. abs(Bij) .lt. 1.0d-14) cycle
         gamma = 0.25d0*atan2(Bij, -Aij)
         if (abs(gamma) .lt. 1.0d-12) cycle
         cg = cos(gamma)
         sg = sin(gamma)
         tmp1 = cg*C_loc(:,i) + sg*C_loc(:,j)
         tmp2 = -sg*C_loc(:,i) + cg*C_loc(:,j)
         C_loc(:,i) = tmp1
         C_loc(:,j) = tmp2
         tmp1 = cg*SC(:,i) + sg*SC(:,j)
         tmp2 = -sg*SC(:,i) + cg*SC(:,j)
         SC(:,i) = tmp1
         SC(:,j) = tmp2
      enddo
   enddo
   func_new = pm_functional(nConts, nOcc, nAtoms, atom_start, atom_count, C_loc, SC)
   if (abs(func_new-func_old) .lt. PM_CONV_TOL*max(abs(func_new),1.0d0)) then
      if (present(n_sweeps_used)) n_sweeps_used = sweep
      deallocate(SC, tmp1, tmp2, atom_start, atom_count)
      return
   endif
   func_old = func_new
enddo
if (present(n_sweeps_used)) n_sweeps_used = PM_MAX_SWEEPS
deallocate(SC, tmp1, tmp2, atom_start, atom_count)
end subroutine pipek_mezey_localize

function pm_functional(nConts, nOcc, nAtoms, atom_start, atom_count, C_loc, SC) result(val)
implicit none
integer,intent(in) :: nConts, nOcc, nAtoms
integer,intent(in) :: atom_start(nAtoms), atom_count(nAtoms)
real(8),intent(in) :: C_loc(nConts,nOcc), SC(nConts,nOcc)
real(8) :: val
integer :: iat, i
real(8) :: Qi

val = 0.0d0
do i = 1,nOcc
   do iat = 1,nAtoms
      if (atom_count(iat) .eq. 0) cycle
      associate(a0 => atom_start(iat), a1 => atom_start(iat)+atom_count(iat)-1)
         Qi = dot_product(C_loc(a0:a1,i), SC(a0:a1,i))
      end associate
      val = val + Qi*Qi
   enddo
enddo
end function pm_functional

end module mod_localization
