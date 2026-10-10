! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_cpscf: the analytic Hessian drivers (RHF/UHF) and the fixed-density

submodule (mod_cpscf) cpscf_hessian_impl
implicit none
contains

module subroutine cpscf_hess_2e(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e, Cocc_a, Cocc_b)
use mod_integrals, only: integrals_hessian_2e_direct, integrals_hessian_2e_df_j, &
                         engine_use_df_j, cosx_hess_exchange, integrals_df_j_fock_deriv
use mod_exchange, only: cosx_enabled, RS_omega, RS_beta
implicit none
integer,intent(in) :: nConts, natoms
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts), HF_exchange_frac
real(8),intent(out) :: d2E_2e(3,natoms,3,natoms)
real(8),intent(in),optional :: Cocc_a(:,:), Cocc_b(:,:)
real(8),allocatable :: d2a(:,:,:,:), d2b(:,:,:,:)
cpscf_cosx_kx_proj = .false.
if (cosx_enabled .and. HF_exchange_frac > 0.0d0) then
   if (engine_use_df_j .and. (cpscf_match_scf_requested() .or. cpscf_rij_mode(HF_exchange_frac))) then
      call integrals_hessian_2e_df_j(nConts, natoms, Pa+Pb, d2E_2e)
   else
      call integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, 0.0d0, d2E_2e)
   endif
   allocate(d2a(3,natoms,3,natoms), d2b(3,natoms,3,natoms))
   if (allocated(cpscf_cosx_kxa)) deallocate(cpscf_cosx_kxa, cpscf_cosx_kxb)
   allocate(cpscf_cosx_kxa(nConts,nConts,3,natoms), cpscf_cosx_kxb(nConts,nConts,3,natoms))
   cpscf_cosx_kx_proj = present(Cocc_a)
   if (cpscf_cosx_kx_proj) then
      cpscf_cosx_kxa = 0.0d0 ; cpscf_cosx_kxb = 0.0d0
      if (present(Cocc_b)) then
         call cosx_hess_exchange(nConts, natoms, Pa, Pb, 0.0d0, d2a, d2b, cpscf_cosx_kxa, cpscf_cosx_kxb, &
                                 Cocca=Cocc_a, Coccb=Cocc_b)
      else
         call cosx_hess_exchange(nConts, natoms, Pa, Pb, 0.0d0, d2a, d2b, cpscf_cosx_kxa, cpscf_cosx_kxb, &
                                 Cocca=Cocc_a, Coccb=Cocc_a)
      endif
   else
      call cosx_hess_exchange(nConts, natoms, Pa, Pb, 0.0d0, d2a, d2b, cpscf_cosx_kxa, cpscf_cosx_kxb)
   endif
   d2E_2e = d2E_2e + 2.0d0*HF_exchange_frac*(d2a + d2b)
   if (RS_omega > 0.0d0) then
      block
         real(8),allocatable :: d2al(:,:,:,:), d2bl(:,:,:,:), kal(:,:,:,:), kbl(:,:,:,:)
         allocate(d2al(3,natoms,3,natoms), d2bl(3,natoms,3,natoms))
         allocate(kal(nConts,nConts,3,natoms), kbl(nConts,nConts,3,natoms))
         kal = 0.0d0 ; kbl = 0.0d0
         if (cpscf_cosx_kx_proj) then
            if (present(Cocc_b)) then
               call cosx_hess_exchange(nConts, natoms, Pa, Pb, RS_omega, d2al, d2bl, kal, kbl, Cocca=Cocc_a, Coccb=Cocc_b)
            else
               call cosx_hess_exchange(nConts, natoms, Pa, Pb, RS_omega, d2al, d2bl, kal, kbl, Cocca=Cocc_a, Coccb=Cocc_a)
            endif
         else
            call cosx_hess_exchange(nConts, natoms, Pa, Pb, RS_omega, d2al, d2bl, kal, kbl)
         endif
         d2E_2e = d2E_2e + 2.0d0*RS_beta*(d2al + d2bl)
         cpscf_cosx_kxa = cpscf_cosx_kxa + (RS_beta/HF_exchange_frac)*kal
         cpscf_cosx_kxb = cpscf_cosx_kxb + (RS_beta/HF_exchange_frac)*kbl
      end block
   endif
   deallocate(d2a, d2b)
else if (engine_use_df_j .and. HF_exchange_frac == 0.0d0) then
   call integrals_hessian_2e_df_j(nConts, natoms, Pa+Pb, d2E_2e)
else
   call integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e)
   if (RS_omega > 0.0d0) then
      block
         real(8),allocatable :: d2lr(:,:,:,:)
         allocate(d2lr(3,natoms,3,natoms))
         call integrals_hessian_2e_direct(nConts, natoms, Pa, Pb, RS_beta, d2lr, omega=RS_omega)
         d2E_2e = d2E_2e + d2lr
         deallocate(d2lr)
      end block
   endif
endif
if (allocated(cpscf_dfj_jx)) deallocate(cpscf_dfj_jx)
if (cpscf_rij_mode(HF_exchange_frac)) then
   allocate(cpscf_dfj_jx(nConts,nConts,3,natoms))
   call integrals_df_j_fock_deriv(nConts, natoms, Pa+Pb, cpscf_dfj_jx)
endif
end subroutine cpscf_hess_2e

module subroutine cpscf_hess_2e_release()
if (allocated(cpscf_cosx_kxa)) deallocate(cpscf_cosx_kxa, cpscf_cosx_kxb)
cpscf_cosx_kx_proj = .false.
if (allocated(cpscf_dfj_jx)) deallocate(cpscf_dfj_jx)
end subroutine cpscf_hess_2e_release

module subroutine cphf_debug_hessian_row(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
use MOL_info, only: atoms
use mod_integrals, only: integrals_hessian_1e_shell, integrals_hessian_2e_direct, &
                         integrals_hessian_2e_df_j, engine_use_df_j, df_store_release, df_store_restore
use mod_xc_hess, only: xc_hess_direct
use mod_cosmo_hessian, only: cosmo_hess_release
use mod_xc, only: xc_is_rsh => is_rsh
use mod_exchange, only: cosx_enabled_g => cosx_enabled
use mod_integrals, only: engine_use_df_j_g => engine_use_df_j
use mod_profile, only: prof_start, prof_stop, prof_report, profile_enabled
use mod_meminfo, only: release_free_memory
implicit none
integer,intent(in) :: nConts, natoms, nocc, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

real(8),allocatable :: U_full(:,:,:,:), h1_MO_all(:,:,:,:), s1_MO_all(:,:,:,:), mo_e1_all(:,:,:,:)
real(8),allocatable :: d2E_1e(:,:,:,:), d2E_2e(:,:,:,:), d2E_resp(:,:,:,:), d2Enn(:,:,:,:), d2E_tot(:,:,:,:)
real(8),allocatable :: d2E_xc(:,:,:,:), d2E_cosmo(:,:,:,:)
real(8),allocatable :: Ptot(:,:), Wtot(:,:)
integer,allocatable :: ao_atom(:), znum(:), ecpce(:)
real(8),allocatable :: coor_bohr(:,:)
INCLUDE 'parameter.h'
integer :: mu,i,ja,y
real(8) :: t0, t_cphf, t_resp, t_1e, t_2e, t_nn
real(8) :: omp_get_wtime
external :: omp_get_wtime

t0 = omp_get_wtime()
if (xc_is_rsh .and. .not. cosx_enabled_g .and. engine_use_df_j_g) error stop &
   'analytic Hessian: range-separated hybrids with exact K need exact J (RI-J not implemented)'
allocate(d2E_2e(3,natoms,3,natoms))
call prof_start("hess_2e_direct")
call df_store_release()
call release_free_memory()
call cpscf_hess_2e(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e, Cocc_a=C_a(:,1:nocc))
call prof_stop("hess_2e_direct")
t_2e = omp_get_wtime() - t0

t0 = omp_get_wtime()
block
   use mod_cosmo_state, only: cosmo_enabled, saved_cavity_type
   use mod_cosmo_hessian, only: cosmo_hess_prepare
   integer,allocatable :: ao_c(:)
   real(8),allocatable :: co_c(:,:)
   integer :: ic0, jc0, mc0
   if (cosmo_enabled) then
      if (trim(saved_cavity_type) /= 'gepol') error stop 'analytic COSMO Hessian: only cavity_type=gepol is implemented'
      allocate(ao_c(nConts), co_c(3,natoms), d2E_cosmo(3,natoms,3,natoms))
      mc0 = 0
      do ic0 = 1,natoms
         co_c(:,ic0) = atoms(ic0)%coor*ans2bohr
         do jc0 = 1,atoms(ic0)%nconts
            mc0 = mc0+1
            ao_c(mc0) = ic0
         enddo
      enddo
      call cosmo_hess_prepare(nConts, natoms, ao_c, co_c, Pa+Pb, d2E_cosmo)
   endif
end block
allocate(U_full(nConts,nocc,3,natoms))
allocate(h1_MO_all(nConts,nocc,3,natoms), s1_MO_all(nConts,nocc,3,natoms), mo_e1_all(nocc,nocc,3,natoms))
block
   character(len=8) :: e
   call get_environment_variable("ENGINE_CPSCF_DF_STORE", e)
   if (trim(e) .eq. "1") call df_store_restore()
end block
call cphf_solve_rhf(nConts, natoms, nocc, C_a, eLev_a, Pa, Pb, HF_exchange_frac, U_full, h1_MO_all, s1_MO_all, mo_e1_all)
call df_store_release()
call release_free_memory()
t_cphf = omp_get_wtime() - t0

t0 = omp_get_wtime()
allocate(d2E_resp(3,natoms,3,natoms))
call cphf_hess_response(nConts, natoms, nocc, U_full, h1_MO_all, s1_MO_all, mo_e1_all, eLev_a, d2E_resp)
t_resp = omp_get_wtime() - t0

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

t0 = omp_get_wtime()
allocate(d2E_1e(3,natoms,3,natoms))
call prof_start("hess_1e_direct")
call calc_hessian_1e(nConts, natoms, ao_atom, Ptot, Wtot, znum, ecpce, coor_bohr, d2E_1e)
call prof_stop("hess_1e_direct")
t_1e = omp_get_wtime() - t0

call cpscf_hess_2e_release()

t0 = omp_get_wtime()
allocate(d2Enn(3,natoms,3,natoms))
call cphf_hess_nuc(natoms, znum, ecpce, coor_bohr, d2Enn)
t_nn = omp_get_wtime() - t0

allocate(d2E_xc(3,natoms,3,natoms))
d2E_xc = 0.0d0
if (HF_exchange_frac < 1.0d0) then
   t0 = omp_get_wtime()
   call prof_start("hess_xc_direct")
   call xc_hess_direct(nConts, natoms, Ptot, d2E_xc)
   call prof_stop("hess_xc_direct")
   t_2e = t_2e
endif

allocate(d2E_tot(3,natoms,3,natoms))
d2E_tot = d2E_1e + d2E_2e + d2E_resp + d2Enn + d2E_xc
if (allocated(d2E_cosmo)) then
   d2E_tot = d2E_tot + d2E_cosmo
   deallocate(d2E_cosmo)
   call cosmo_hess_release()
endif

call cpscf_warn_mixed_direct()
call prof_report()
if (profile_enabled()) then
print *, 'HESSIAN STAGE TIMING (s):'
print *, '   cphf_solve_rhf            =', t_cphf
print *, '      build_amat             =', cphf_t_amat
print *, '      force integrals (1x)   =', cphf_t_forceint
print *, '      fock_deriv (36 pert)   =', cphf_t_fockderiv
print *, '      MO transforms (36 pert)=', cphf_t_mo
print *, '      pass1b occ-occ response=', cphf_t_pass1
print *, '      fixed-point solves     =', cphf_t_iter, '  iters total =', cphf_n_iter_total
print *, '        response J / COSX K / XC =', cphf_t_vj, cphf_t_vk, cphf_t_vxc
print *, '        COSX K calls / densities =', cphf_n_vk, cphf_nd_vk
print *, '      pass3 observables      =', cphf_t_pass3
print *, '   cphf_hess_response        =', t_resp
print *, '   calc_hessian_1e           =', t_1e
print *, '   integrals_hessian_2e_direct=', t_2e
print *, '   cphf_hess_nuc             =', t_nn
print *, '   stage total               =', t_cphf+t_resp+t_1e+t_2e+t_nn
endif
if (profile_enabled()) then
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
endif

block
   use mod_hessian_freq, only: hessian_write_full, hessian_frequencies, hess_classical, freq_unit, freq_base
   real(8),allocatable :: freq_cm(:), disp_m(:,:)
   integer :: nfreq, uh, kk, nzero_m
   character(264) :: hessfile
   if (allocated(hess_classical)) then
      if (all(shape(hess_classical) == shape(d2E_tot))) d2E_tot = d2E_tot + hess_classical
   endif
   d2E_tot = 0.5d0*(d2E_tot + reshape(d2E_tot,(/3,natoms,3,natoms/),order=(/3,4,1,2/)))

   hessfile = trim(freq_base)//'.hess'
   open(newunit=uh, file=trim(hessfile), status='replace', action='write')
   write(uh,'(A)') ''
   write(uh,'(A)') '$direwolf_hessian_file'
   call hessian_write_full(uh, natoms, d2E_tot)
   close(uh)
   print *, 'Full Hessian written to ', trim(hessfile)
   write(freq_unit,'(A)') ''
   write(freq_unit,'(A)') ' Hessian (Cartesian, Eh/bohr^2, 3N x 3N) is in '//trim(hessfile)

   allocate(freq_cm(3*natoms))
   allocate(disp_m(3*natoms,3*natoms))
   call hessian_frequencies(natoms, znum, coor_bohr, d2E_tot, freq_cm, nfreq, disp_m, nzero_m)
   write(freq_unit,'(A)') ' [FREQUENCIES]'
   write(freq_unit,'(A)') ' HARMONIC FREQUENCIES (cm-1), translations/rotations first, then vibrations:'
   do kk = 1,nfreq
      write(freq_unit,'(A,I5,F14.4)') '   mode', kk, freq_cm(kk)
   enddo
   write(freq_unit,'(A)') ' [FREQUENCIESEND]'
   call hess_ir_report(nConts, natoms, znum, ecpce, coor_bohr, ao_atom, Ptot, 4.0d0, &
                       C_a, nocc, U_full, C_a, 0, U_full, freq_cm, disp_m, nzero_m)
   deallocate(freq_cm, disp_m)
end block

deallocate(U_full,h1_MO_all,s1_MO_all,mo_e1_all,d2E_resp,ao_atom,Ptot,Wtot,znum,ecpce,coor_bohr)
deallocate(d2E_1e,d2E_2e,d2Enn,d2E_tot)
end subroutine cphf_debug_hessian_row

module subroutine cphf_debug_hessian_row_uhf(nConts, natoms, nocc_a, nocc_b, C_a, eLev_a, C_b, eLev_b, &
                                       Pa, Pb, HF_exchange_frac, iatom, ic)
use MOL_info, only: atoms
use mod_cosmo_hessian, only: cosmo_hess_release
use mod_integrals, only: integrals_hessian_1e_shell, integrals_hessian_2e_direct, &
                         integrals_hessian_2e_df_j, engine_use_df_j, df_store_release, df_store_restore
use mod_profile, only: prof_start, prof_stop, prof_report, profile_enabled
use mod_meminfo, only: release_free_memory
use mod_xc_hess, only: xc_hess_direct_u
use mod_xc, only: xc_is_rsh => is_rsh
use mod_exchange, only: cosx_enabled_g => cosx_enabled
use mod_integrals, only: engine_use_df_j_g => engine_use_df_j
implicit none
integer,intent(in) :: nConts, natoms, nocc_a, nocc_b, iatom, ic
real(8),intent(in) :: C_a(nConts,nConts), eLev_a(nConts)
real(8),intent(in) :: C_b(nConts,nConts), eLev_b(nConts)
real(8),intent(in) :: Pa(nConts,nConts), Pb(nConts,nConts)
real(8),intent(in) :: HF_exchange_frac

real(8),allocatable :: U_full_a(:,:,:,:), U_full_b(:,:,:,:)
real(8),allocatable :: h1_MO_all_a(:,:,:,:), s1_MO_all_a(:,:,:,:), mo_e1_all_a(:,:,:,:)
real(8),allocatable :: h1_MO_all_b(:,:,:,:), s1_MO_all_b(:,:,:,:), mo_e1_all_b(:,:,:,:)
real(8),allocatable :: d2E_1e(:,:,:,:), d2E_2e(:,:,:,:), d2E_resp(:,:,:,:), d2Enn(:,:,:,:), d2E_tot(:,:,:,:)
real(8),allocatable :: d2E_xc(:,:,:,:), d2E_cosmo(:,:,:,:)
real(8),allocatable :: Ptot(:,:), Wtot(:,:)
integer,allocatable :: ao_atom(:), znum(:), ecpce(:)
real(8),allocatable :: coor_bohr(:,:)
INCLUDE 'parameter.h'
integer :: mu,i,ja,y
real(8) :: t0, t_cphf, t_resp, t_1e, t_2e, t_nn
real(8) :: omp_get_wtime
external :: omp_get_wtime

t0 = omp_get_wtime()
if (xc_is_rsh .and. .not. cosx_enabled_g .and. engine_use_df_j_g) error stop &
   'analytic Hessian: range-separated hybrids with exact K need exact J (RI-J not implemented)'
allocate(d2E_2e(3,natoms,3,natoms))
call prof_start("hess_2e_direct")
call df_store_release()
call release_free_memory()
call cpscf_hess_2e(nConts, natoms, Pa, Pb, HF_exchange_frac, d2E_2e, Cocc_a=C_a(:,1:nocc_a), Cocc_b=C_b(:,1:nocc_b))
call prof_stop("hess_2e_direct")
t_2e = omp_get_wtime() - t0

t0 = omp_get_wtime()
block
   use mod_cosmo_state, only: cosmo_enabled, saved_cavity_type
   use mod_cosmo_hessian, only: cosmo_hess_prepare
   integer,allocatable :: ao_c(:)
   real(8),allocatable :: co_c(:,:)
   integer :: ic0, jc0, mc0
   if (cosmo_enabled) then
      if (trim(saved_cavity_type) /= 'gepol') error stop 'analytic COSMO Hessian: only cavity_type=gepol is implemented'
      allocate(ao_c(nConts), co_c(3,natoms), d2E_cosmo(3,natoms,3,natoms))
      mc0 = 0
      do ic0 = 1,natoms
         co_c(:,ic0) = atoms(ic0)%coor*ans2bohr
         do jc0 = 1,atoms(ic0)%nconts
            mc0 = mc0+1
            ao_c(mc0) = ic0
         enddo
      enddo
      call cosmo_hess_prepare(nConts, natoms, ao_c, co_c, Pa+Pb, d2E_cosmo)
   endif
end block
allocate(U_full_a(nConts,nocc_a,3,natoms), U_full_b(nConts,nocc_b,3,natoms))
allocate(h1_MO_all_a(nConts,nocc_a,3,natoms), s1_MO_all_a(nConts,nocc_a,3,natoms), mo_e1_all_a(nocc_a,nocc_a,3,natoms))
allocate(h1_MO_all_b(nConts,nocc_b,3,natoms), s1_MO_all_b(nConts,nocc_b,3,natoms), mo_e1_all_b(nocc_b,nocc_b,3,natoms))
block
   character(len=8) :: e
   call get_environment_variable("ENGINE_CPSCF_DF_STORE", e)
   if (trim(e) .eq. "1") call df_store_restore()
end block
call cpscf_solve_uhf(nConts, natoms, nocc_a, nocc_b, C_a, eLev_a, C_b, eLev_b, Pa, Pb, &
                     HF_exchange_frac, U_full_a, U_full_b, &
                     h1_MO_all_a, s1_MO_all_a, mo_e1_all_a, h1_MO_all_b, s1_MO_all_b, mo_e1_all_b)
call df_store_release()
call release_free_memory()
t_cphf = omp_get_wtime() - t0

t0 = omp_get_wtime()
allocate(d2E_resp(3,natoms,3,natoms))
call cphf_hess_response_spin(nConts, natoms, nocc_a, nocc_b, &
                             U_full_a, U_full_b, &
                             h1_MO_all_a, s1_MO_all_a, mo_e1_all_a, eLev_a, &
                             h1_MO_all_b, s1_MO_all_b, mo_e1_all_b, eLev_b, d2E_resp)
t_resp = omp_get_wtime() - t0

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
      do i = 1,nocc_a
         Wtot(mu,y) = Wtot(mu,y) + eLev_a(i)*C_a(mu,i)*C_a(y,i)
      enddo
      do i = 1,nocc_b
         Wtot(mu,y) = Wtot(mu,y) + eLev_b(i)*C_b(mu,i)*C_b(y,i)
      enddo
   enddo
enddo

allocate(znum(natoms), ecpce(natoms), coor_bohr(3,natoms))
do i = 1,natoms
   znum(i) = atoms(i)%charge
   ecpce(i) = atoms(i)%ecpCoreElec
   coor_bohr(:,i) = atoms(i)%coor*ans2bohr
enddo

t0 = omp_get_wtime()
allocate(d2E_1e(3,natoms,3,natoms))
call prof_start("hess_1e_direct")
call calc_hessian_1e(nConts, natoms, ao_atom, Ptot, Wtot, znum, ecpce, coor_bohr, d2E_1e)
call prof_stop("hess_1e_direct")
t_1e = omp_get_wtime() - t0

call cpscf_hess_2e_release()

t0 = omp_get_wtime()
allocate(d2Enn(3,natoms,3,natoms))
call cphf_hess_nuc(natoms, znum, ecpce, coor_bohr, d2Enn)
t_nn = omp_get_wtime() - t0

allocate(d2E_xc(3,natoms,3,natoms))
d2E_xc = 0.0d0
if (HF_exchange_frac < 1.0d0) then
   call prof_start("hess_xc_direct")
   call xc_hess_direct_u(nConts, natoms, Pa, Pb, d2E_xc)
   call prof_stop("hess_xc_direct")
endif
allocate(d2E_tot(3,natoms,3,natoms))
d2E_tot = d2E_1e + d2E_2e + d2E_resp + d2Enn + d2E_xc
if (allocated(d2E_cosmo)) then
   d2E_tot = d2E_tot + d2E_cosmo
   deallocate(d2E_cosmo)
   call cosmo_hess_release()
endif

call cpscf_warn_mixed_direct()
call prof_report()
if (profile_enabled()) then
print *, 'HESSIAN STAGE TIMING (s):'
print *, '   cpscf_solve_uhf           =', t_cphf
print *, '      force integrals (1x)   =', cphf_t_forceint
print *, '      fock_deriv (36 pert)   =', cphf_t_fockderiv
print *, '      MO transforms (36 pert)=', cphf_t_mo
print *, '      pass1b occ-occ response=', cphf_t_pass1
print *, '      fixed-point solves     =', cphf_t_iter, '  iters total =', cphf_n_iter_total
print *, '        response J / COSX K / XC =', cphf_t_vj, cphf_t_vk, cphf_t_vxc
print *, '        COSX K calls / densities =', cphf_n_vk, cphf_nd_vk
print *, '      pass3 observables      =', cphf_t_pass3
print *, '   cphf_hess_response_spin   =', t_resp
print *, '   calc_hessian_1e           =', t_1e
print *, '   integrals_hessian_2e_direct=', t_2e
print *, '   cphf_hess_nuc             =', t_nn
print *, '   stage total               =', t_cphf+t_resp+t_1e+t_2e+t_nn
endif
if (profile_enabled()) then
print *, 'CPHF DEBUG HESSIAN ROW (UHF): perturbation atom=',iatom,' comp=',ic
do ja = 1,natoms
   do y = 1,3
      print *, 'target atom=',ja,' comp=',y,' d2E=', d2E_tot(ic,iatom,y,ja), &
               ' (1e=',d2E_1e(ic,iatom,y,ja),' 2e=',d2E_2e(ic,iatom,y,ja), &
               ' resp=',d2E_resp(ic,iatom,y,ja),' nn=',d2Enn(ic,iatom,y,ja),')'
   enddo
enddo
print *, 'CPHF DEBUG HESSIAN symmetry check: max|d2E(x,a,y,b)-d2E(y,b,x,a)| =', &
   maxval(abs(d2E_tot - reshape(d2E_tot,(/3,natoms,3,natoms/),order=(/3,4,1,2/))))
endif

block
   use mod_hessian_freq, only: hessian_write_full, hessian_frequencies, hess_classical, freq_unit, freq_base
   real(8),allocatable :: freq_cm(:), disp_m(:,:)
   integer :: nfreq, uh, kk, nzero_m
   character(264) :: hessfile
   if (allocated(hess_classical)) then
      if (all(shape(hess_classical) == shape(d2E_tot))) d2E_tot = d2E_tot + hess_classical
   endif
   d2E_tot = 0.5d0*(d2E_tot + reshape(d2E_tot,(/3,natoms,3,natoms/),order=(/3,4,1,2/)))

   hessfile = trim(freq_base)//'.hess'
   open(newunit=uh, file=trim(hessfile), status='replace', action='write')
   write(uh,'(A)') ''
   write(uh,'(A)') '$direwolf_hessian_file'
   call hessian_write_full(uh, natoms, d2E_tot)
   close(uh)
   print *, 'Full Hessian written to ', trim(hessfile)
   write(freq_unit,'(A)') ''
   write(freq_unit,'(A)') ' Hessian (Cartesian, Eh/bohr^2, 3N x 3N) is in '//trim(hessfile)

   allocate(freq_cm(3*natoms))
   allocate(disp_m(3*natoms,3*natoms))
   call hessian_frequencies(natoms, znum, coor_bohr, d2E_tot, freq_cm, nfreq, disp_m, nzero_m)
   write(freq_unit,'(A)') ' [FREQUENCIES]'
   write(freq_unit,'(A)') ' HARMONIC FREQUENCIES (cm-1), translations/rotations first, then vibrations:'
   do kk = 1,nfreq
      write(freq_unit,'(A,I5,F14.4)') '   mode', kk, freq_cm(kk)
   enddo
   write(freq_unit,'(A)') ' [FREQUENCIESEND]'
   call hess_ir_report(nConts, natoms, znum, ecpce, coor_bohr, ao_atom, Ptot, 2.0d0, &
                       C_a, nocc_a, U_full_a, C_b, nocc_b, U_full_b, freq_cm, disp_m, nzero_m)
   deallocate(freq_cm, disp_m)
end block

deallocate(U_full_a,U_full_b,h1_MO_all_a,s1_MO_all_a,mo_e1_all_a)
deallocate(h1_MO_all_b,s1_MO_all_b,mo_e1_all_b,d2E_resp,ao_atom,Ptot,Wtot,znum,ecpce,coor_bohr)
deallocate(d2E_1e,d2E_2e,d2Enn,d2E_tot)
end subroutine cphf_debug_hessian_row_uhf

subroutine hess_ir_report(nConts, natoms, znum, ecpce, coor_bohr, ao_atom, Ptot, fac, &
                          C1, nocc1, U1, C2, nocc2, U2, freq_cm, disp, nzero)
use mod_integrals, only: integrals_dipole_ao, integrals_dipole_pulay
use mod_hessian_freq, only: print_normal_modes, ir_intensities, print_ir, write_molden_freq, write_modes_csv, freq_unit, &
                           print_thermo, freq_mult, freq_base, freq_thermo
implicit none
integer,intent(in) :: nConts, natoms, znum(natoms), ecpce(natoms), ao_atom(nConts), nocc1, nocc2, nzero
real(8),intent(in) :: coor_bohr(3,natoms), Ptot(nConts,nConts), fac
real(8),intent(in) :: C1(nConts,nConts), C2(nConts,nConts)
real(8),intent(in) :: U1(nConts,max(nocc1,1),3,natoms), U2(nConts,max(nocc2,1),3,natoms)
real(8),intent(in) :: freq_cm(3*natoms), disp(3*natoms,3*natoms)
real(8),allocatable :: R(:,:,:), M(:,:), RC(:,:), dmu(:,:), inten(:), redmass(:)
real(8) :: dpul(3,3,natoms)
integer :: c, x, a, i, ispin, nocc, i_at
n3: block
   integer :: n3
   n3 = 3*natoms
   allocate(R(nConts,nConts,3), dmu(3,n3), inten(n3), redmass(n3))
   call integrals_dipole_ao(nConts, R)
   call integrals_dipole_pulay(nConts, natoms, ao_atom, Ptot, dpul)
   dmu = 0.0d0
   do a = 1,natoms
      do x = 1,3
         i = 3*(a-1)+x
         do c = 1,3
            dmu(c,i) = dpul(c,x,a)
         enddo
         dmu(x,i) = dmu(x,i) + real(znum(a)-ecpce(a),8)
      enddo
   enddo
   do ispin = 1,2
      nocc = nocc1
      if (ispin == 2) nocc = nocc2
      if (nocc <= 0) cycle
      allocate(M(nConts,nocc), RC(nConts,nocc))
      do c = 1,3
         if (ispin == 1) then
            call dgemm('N','N',nConts,nocc,nConts,1.0d0,R(:,:,c),nConts,C1,nConts,0.0d0,RC,nConts)
            call dgemm('T','N',nConts,nocc,nConts,1.0d0,C1,nConts,RC,nConts,0.0d0,M,nConts)
         else
            call dgemm('N','N',nConts,nocc,nConts,1.0d0,R(:,:,c),nConts,C2,nConts,0.0d0,RC,nConts)
            call dgemm('T','N',nConts,nocc,nConts,1.0d0,C2,nConts,RC,nConts,0.0d0,M,nConts)
         endif
         do a = 1,natoms
            do x = 1,3
               i = 3*(a-1)+x
               if (ispin == 1) then
                  dmu(c,i) = dmu(c,i) - fac*sum(U1(:,1:nocc,x,a)*M)
               else
                  dmu(c,i) = dmu(c,i) - fac*sum(U2(:,1:nocc,x,a)*M)
               endif
            enddo
         enddo
      enddo
      deallocate(M, RC)
   enddo
   call ir_intensities(natoms, dmu, disp, inten)
   call print_normal_modes(natoms, znum, freq_cm, disp, nzero, redmass, inten)
   call print_ir(n3, nzero, freq_cm, inten)
   if (freq_thermo) call print_thermo(natoms, znum, coor_bohr, freq_cm, nzero, freq_mult)
   write(freq_unit,'(A)') ' [DIPOLEDERIV] Dipole derivatives d(mu_c)/d(R_atom,k), unit e'
   write(freq_unit,'(A)') '  Atom  mu_c           d/dx        d/dy        d/dz'
   do i_at = 1,natoms
      do c = 1,3
         write(freq_unit,'(I6,2X,A,3F12.4)') i_at, ('mu_'//achar(119+c)//'  '), dmu(c,3*(i_at-1)+1:3*(i_at-1)+3)
      enddo
   enddo
   write(freq_unit,'(A)') ' [DIPOLEDERIVEND]'
   call write_molden_freq(trim(freq_base)//'.freq.molden', natoms, znum, coor_bohr, freq_cm, disp, &
                          redmass, nzero, inten)
   call write_modes_csv(trim(freq_base)//'.modes.csv', natoms, znum, freq_cm, disp, redmass, nzero, inten)
   print *, 'Vibrational modes written to '//trim(freq_base)//'.freq.molden and '// &
            trim(freq_base)//'.modes.csv'
end block n3
end subroutine hess_ir_report

end submodule cpscf_hessian_impl
