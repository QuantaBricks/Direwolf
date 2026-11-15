! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: the dipole polarizability (static and, for pure

submodule (mod_cpscf) cpscf_polar_impl
implicit none
contains

module subroutine cphf_solve_efield(nConts, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, alpha, niter, omega)
use mod_integrals, only: integrals_dipole_ao, engine_use_df_j
use mod_exchange, only: RS_omega, cosx_enabled
implicit none
integer,intent(in) :: nConts, nocc
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac
real(8),intent(out) :: alpha(3,3)
integer,intent(out) :: niter
real(8),intent(in),optional :: omega
integer :: nvir, x, y
real(8) :: w
character(len=16) :: solv
real(8),allocatable :: R(:,:,:), T(:,:), RHS(:,:,:), U(:,:,:)

if (RS_omega > 0.0d0 .and. .not. cosx_enabled) &
   error stop 'polarizability: range-separated hybrids with exact exchange are not implemented'
nvir = nConts - nocc
w = 0.0d0
if (present(omega)) w = omega
if (w /= 0.0d0) then
   if (HF_exchange_frac /= 0.0d0 .or. RS_omega > 0.0d0) &
      error stop 'dynamic polarizability: functionals with exact exchange are not implemented (pure functionals only)'
   call get_environment_variable("ENGINE_CPHF_SOLVER", solv)
   if (trim(solv) == 'diis') error stop 'dynamic polarizability: ENGINE_CPHF_SOLVER=diis has no frequency shift'
   if (w >= eLev_a(nocc+1) - eLev_a(nocc)) &
      error stop 'dynamic polarizability: omega is not below the HOMO-LUMO gap (resonance region)'
endif
cphf_omega2 = w*w
cpscf_xc_active = (HF_exchange_frac < 1.0d0)
if (cpscf_xc_active) then
   if (allocated(cpscf_xc_Ptot)) deallocate(cpscf_xc_Ptot)
   allocate(cpscf_xc_Ptot(nConts,nConts))
   cpscf_xc_Ptot = Pa + Pb
endif
cpscf_force_dfj = engine_use_df_j
cphf_conv_override = 1.0d-8
cphf_voo_track = .false.
cphf_t_vj = 0.0d0 ; cphf_t_vk = 0.0d0 ; cphf_t_vxc = 0.0d0
cphf_n_vk = 0 ; cphf_nd_vk = 0

allocate(R(nConts,nConts,3), T(nConts,nocc), RHS(nvir,nocc,3), U(nvir,nocc,3))
call integrals_dipole_ao(nConts, R)
do x = 1,3
   call dgemm('N','N',nConts,nocc,nConts,1.0d0,R(1,1,x),nConts,C_a,nConts,0.0d0,T,nConts)
   call dgemm('T','N',nvir,nocc,nConts,1.0d0,C_a(1,nocc+1),nConts,T,nConts,0.0d0,RHS(1,1,x),nvir)
enddo
call cphf_solve_vo_direct(nConts, nocc, nvir, 3, C_a, eLev_a, HF_exchange_frac, RHS, U, niter)
do x = 1,3
   do y = 1,3
      alpha(x,y) = -4.0d0*sum(RHS(:,:,x)*U(:,:,y))
   enddo
enddo
alpha = 0.5d0*(alpha + transpose(alpha))
cpscf_force_dfj = .false.
cphf_conv_override = -1.0d0
cphf_omega2 = 0.0d0
cpscf_xc_active = .false.
deallocate(R, T, RHS, U)
end subroutine cphf_solve_efield

end submodule cpscf_polar_impl
