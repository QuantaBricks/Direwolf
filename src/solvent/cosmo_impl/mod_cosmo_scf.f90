! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cosmo_scf: the per-SCF-iteration COSMO step (solve for tessera

module mod_cosmo_scf
use mod_cosmo_cavity, only: bohr2ang
use mod_cosmo_state
implicit none

contains

subroutine cosmo_set_sigma_profile_debug(r_av_ang, output_file)
implicit none
real(8),intent(in) :: r_av_ang
character(*),intent(in) :: output_file
debug_sigma_rav = r_av_ang
debug_sigma_profile_file = trim(output_file)
end subroutine cosmo_set_sigma_profile_debug

subroutine cosmo_set_cosmors_request(on, cosmo_file, compound_name)
implicit none
logical,intent(in) :: on
character(*),intent(in) :: cosmo_file, compound_name
cosmors_request = on
cosmors_cosmo_file_req = trim(cosmo_file)
cosmors_compound_name_req = trim(compound_name)
end subroutine cosmo_set_cosmors_request

subroutine cosmo_scf_step(nConts, Ptot, Fock_cosmo, E_cosmo_nuc, iter)
use mod_integrals, only: cosmo_build_one_tess_matrix
use omp_lib, only: omp_get_thread_num, omp_get_num_threads
implicit none
integer,intent(in) :: nConts
real(8),intent(in) :: Ptot(nConts,nConts)
real(8),intent(out) :: Fock_cosmo(nConts,nConts)
real(8),intent(out) :: E_cosmo_nuc
integer,intent(in) :: iter

real(8),allocatable :: V_elec(:), Vtot(:), q0(:)
real(8),allocatable :: Bk(:,:), Fock_local(:,:)
integer :: k, info
real(8) :: lambda

allocate(V_elec(n_tess), Vtot(n_tess), q0(n_tess))

if (store_mode) then
   call dgemv('T', nConts*nConts, n_tess, 1.0d0, B_store, nConts*nConts, &
              Ptot, 1, 0.0d0, V_elec, 1)
else
   !$omp parallel private(Bk)
   allocate(Bk(nConts,nConts))
   !$omp do schedule(dynamic)
   do k = 1,n_tess
      call cosmo_build_one_tess_matrix(tess_coor(:,k), Bk)
      V_elec(k) = sum(Ptot*Bk)
   enddo
   !$omp end do
   deallocate(Bk)
   !$omp end parallel
endif

Vtot = V_nuc + V_elec

q0 = -Vtot
call dsytrs('U', n_tess, 1, A_fac, n_tess, A_ipiv, q0, n_tess, info)

lambda = (saved_Qm + sum(q0))/saved_sum_Ainv1
q_last = f_eps*(q0 - lambda*saved_Ainv1)

if (store_mode) then
   call dgemv('N', nConts*nConts, n_tess, 1.0d0, B_store, nConts*nConts, &
              q_last, 1, 0.0d0, Fock_cosmo, 1)
else
   Fock_cosmo = 0.0d0
   !$omp parallel private(Bk,Fock_local)
   allocate(Bk(nConts,nConts), Fock_local(nConts,nConts))
   Fock_local = 0.0d0
   !$omp do schedule(dynamic)
   do k = 1,n_tess
      call cosmo_build_one_tess_matrix(tess_coor(:,k), Bk)
      Fock_local = Fock_local + q_last(k)*Bk
   enddo
   !$omp end do
   block
      integer :: merge_tid, merge_nthreads
      merge_nthreads = omp_get_num_threads()
      do merge_tid = 0, merge_nthreads-1
         !$omp barrier
         if (omp_get_thread_num() .eq. merge_tid) then
            Fock_cosmo = Fock_cosmo + Fock_local
         endif
      enddo
   end block
   deallocate(Bk, Fock_local)
   !$omp end parallel
endif

E_cosmo_nuc = 0.5d0*sum(q_last*V_nuc)
if (allocated(Vtot_last)) deallocate(Vtot_last)
call move_alloc(Vtot, Vtot_last)

deallocate(V_elec, q0)
end subroutine cosmo_scf_step

subroutine cosmo_report_sigma_profile(r_av_ang, output_file)
implicit none
real(8),intent(in) :: r_av_ang
character(*),intent(in) :: output_file

real(8),allocatable :: sigma_bins(:), hist_area(:)
integer :: ib, iu
real(8),parameter :: sigma_lo = -0.025d0, sigma_hi = 0.025d0, sigma_step = 0.0006d0

if (len_trim(output_file) == 0) return

end subroutine cosmo_report_sigma_profile

subroutine cosmo_write_dot_cosmo_file(natoms, coor_ang, method_functional, method_basis, E_tot)
implicit none
integer,intent(in) :: natoms
real(8),intent(in) :: coor_ang(3,natoms)
character(*),intent(in) :: method_functional, method_basis
real(8),intent(in) :: E_tot

end subroutine cosmo_write_dot_cosmo_file

subroutine cosmo_finalize()
implicit none
cosmo_enabled = .false.
n_tess = 0
if (allocated(tess_coor)) deallocate(tess_coor)
if (allocated(tess_area)) deallocate(tess_area)
if (allocated(tess_atom)) deallocate(tess_atom)
if (allocated(V_nuc)) deallocate(V_nuc)
if (allocated(A_fac)) deallocate(A_fac)
if (allocated(A_ipiv)) deallocate(A_ipiv)
if (allocated(saved_Ainv1)) deallocate(saved_Ainv1)
if (allocated(B_store)) deallocate(B_store)
if (allocated(q_last)) deallocate(q_last)
if (allocated(Vtot_last)) deallocate(Vtot_last)
if (allocated(saved_Z)) deallocate(saved_Z)
if (allocated(saved_ecpCoreElec)) deallocate(saved_ecpCoreElec)
if (allocated(debug_prev_hist)) deallocate(debug_prev_hist)
debug_sigma_profile_file = ''
saved_natoms = 0
saved_cavity_type = ''
use_smd_radii = .false.
saved_smd_alpha = 0.0d0
store_mode = .false.
end subroutine cosmo_finalize

end module mod_cosmo_scf
