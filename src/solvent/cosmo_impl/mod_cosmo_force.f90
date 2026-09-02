! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_force: analytic COSMO force (gepol cavity only, see

module mod_cosmo_force
use mod_cosmo_cavity, only: cosmo_build_cavity_deriv, bohr2ang, PI_cosmo
use mod_cosmo_state
implicit none

contains

subroutine cosmo_force_step(natoms, coor_ang, ao_atom, nConts, Ptot, force_cosmo)
use mod_cosmo_cavity, only: cosmo_choose_nf
use mod_integrals, only: cosmo_build_one_tess_matrix_shellderiv, integrals_nuc_attraction_deriv
implicit none
integer,intent(in) :: natoms
real(8),intent(in) :: coor_ang(3,natoms)
integer,intent(in) :: ao_atom(nConts)
integer,intent(in) :: nConts
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: force_cosmo(3,natoms)

real(8),allocatable :: coor_bohr(:,:)
real(8),allocatable :: Rvdw_smd(:)
real(8),allocatable :: d_tess_coor(:,:), d_tess_area(:)
integer,allocatable :: d_tess_atom(:)
integer :: d_n_tess
integer,allocatable :: d_tess_nslot(:), d_tess_atom_of_slot(:,:)
real(8),allocatable :: d_tess_dArea(:,:,:), d_tess_dCentroid(:,:,:,:)
real(8),allocatable :: dBk_shell(:,:,:), Dr_pos(:,:,:)
real(8) :: Vk_grad_s(3), rvec(3), dist, Zeff, dArea_pref
integer :: k, s, ia, ja, mu, nu
real(8) :: qk, contrib(3), shell_grad(3)

force_cosmo = 0.0d0
if (.not. cosmo_enabled) return

if (trim(saved_cavity_type) /= 'gepol') then
   print '("COSMO force error: analytic gradient only implemented for cavity_type=''gepol'' (got ''",A,"'')")', &
         trim(saved_cavity_type)
   stop 1
endif

allocate(coor_bohr(3,natoms))
coor_bohr = coor_ang * (1.0d0/bohr2ang)

if (use_smd_radii) then
   block
      use mod_cosmo_cavity, only: cosmo_smd_radius_bohr
      allocate(Rvdw_smd(natoms))
      do ia = 1,natoms
         Rvdw_smd(ia) = cosmo_smd_radius_bohr(saved_Z(ia), saved_smd_alpha)
      enddo
   end block
   call cosmo_build_cavity_deriv(natoms, saved_Z, coor_bohr, saved_radii_scale, saved_avg_area_ang2, &
                                  d_tess_coor, d_tess_area, d_tess_atom, d_n_tess, &
                                  d_tess_nslot, d_tess_atom_of_slot, d_tess_dArea, d_tess_dCentroid, &
                                  Rvdw_smd)
   deallocate(Rvdw_smd)
else
   call cosmo_build_cavity_deriv(natoms, saved_Z, coor_bohr, saved_radii_scale, saved_avg_area_ang2, &
                                  d_tess_coor, d_tess_area, d_tess_atom, d_n_tess, &
                                  d_tess_nslot, d_tess_atom_of_slot, d_tess_dArea, d_tess_dCentroid)
endif

if (d_n_tess /= n_tess) then
   print '("COSMO force error: tessera count changed between energy (",I0,") and force (",I0,&
         &") passes - cavity not reproducible at this geometry")', n_tess, d_n_tess
   stop 1
endif

allocate(Dr_pos(nConts,nConts,3), dBk_shell(nConts,nConts,3))

do k = 1,n_tess
   qk = q_last(k)

   Vk_grad_s = 0.0d0
   do ia = 1,natoms
      Zeff = real(saved_Z(ia)-saved_ecpCoreElec(ia),8)
      rvec = d_tess_coor(:,k) - coor_bohr(:,ia)
      dist = sqrt(sum(rvec*rvec))
      force_cosmo(:,ia) = force_cosmo(:,ia) + qk*Zeff*rvec/dist**3
      Vk_grad_s = Vk_grad_s - Zeff*rvec/dist**3
   enddo
   do s = 1,d_tess_nslot(k)
      ia = d_tess_atom_of_slot(s,k)
      contrib = matmul(Vk_grad_s, d_tess_dCentroid(:,:,s,k))
      force_cosmo(:,ia) = force_cosmo(:,ia) + qk*contrib
   enddo

   call cosmo_build_one_tess_matrix_shellderiv(d_tess_coor(:,k), dBk_shell)
   do mu = 1,nConts
      ia = ao_atom(mu)
      shell_grad = 0.0d0
      do nu = 1,nConts
         shell_grad = shell_grad + Ptot(mu,nu)*dBk_shell(mu,nu,:)
      enddo
      force_cosmo(:,ia) = force_cosmo(:,ia) - qk*2.0d0*shell_grad
   enddo

   call integrals_nuc_attraction_deriv(d_tess_coor(:,k), nConts, Dr_pos)
   contrib = 0.0d0
   do mu = 1,nConts
      do nu = 1,nConts
         contrib = contrib + Ptot(mu,nu)*Dr_pos(mu,nu,:)
      enddo
   enddo
   contrib = -2.0d0*contrib
   do s = 1,d_tess_nslot(k)
      ia = d_tess_atom_of_slot(s,k)
      force_cosmo(:,ia) = force_cosmo(:,ia) + qk*matmul(contrib, d_tess_dCentroid(:,:,s,k))
   enddo

   dArea_pref = -0.5d0*1.07d0*sqrt(4.0d0*PI_cosmo)*d_tess_area(k)**(-1.5d0)
   do s = 1,d_tess_nslot(k)
      ia = d_tess_atom_of_slot(s,k)
      force_cosmo(:,ia) = force_cosmo(:,ia) + (qk*qk/(2.0d0*f_eps))*dArea_pref*d_tess_dArea(:,s,k)
   enddo
   do ja = 1,n_tess
      if (ja == k) cycle
      rvec = d_tess_coor(:,k) - d_tess_coor(:,ja)
      dist = sqrt(sum(rvec*rvec))
      do s = 1,d_tess_nslot(k)
         ia = d_tess_atom_of_slot(s,k)
         contrib = matmul(-rvec/dist**3, d_tess_dCentroid(:,:,s,k))
         force_cosmo(:,ia) = force_cosmo(:,ia) + (qk*q_last(ja)/(2.0d0*f_eps))*contrib
      enddo
      do s = 1,d_tess_nslot(ja)
         ia = d_tess_atom_of_slot(s,ja)
         contrib = matmul(rvec/dist**3, d_tess_dCentroid(:,:,s,ja))
         force_cosmo(:,ia) = force_cosmo(:,ia) + (qk*q_last(ja)/(2.0d0*f_eps))*contrib
      enddo
   enddo
enddo

deallocate(coor_bohr, Dr_pos, dBk_shell)
deallocate(d_tess_coor, d_tess_area, d_tess_atom, d_tess_nslot, d_tess_atom_of_slot, &
           d_tess_dArea, d_tess_dCentroid)
end subroutine cosmo_force_step

end module mod_cosmo_force
