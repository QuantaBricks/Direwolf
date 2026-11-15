! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_engine_report: the structured results EngineUp produces besides the energy/force arguments

module mod_engine_report
implicit none
private
public :: engine_report_t, tddft_set_t, polar_rec_t
public :: rpt, report_reset, report_hook_i
public :: report_energy_only, report_properties, report_polar, report_polar_dynamic
public :: report_tddft_state, report_freq, report_thermo, report_dipole_deriv
public :: RPT_TD_SINGLET, RPT_TD_TRIPLET, RPT_TD_UKS

integer, parameter :: RPT_TD_SINGLET = 1, RPT_TD_TRIPLET = 2, RPT_TD_UKS = 3
integer, parameter :: MAX_DYN = 16

type :: tddft_set_t
   logical :: used = .false.
   integer :: napply = 0
   character(len=256) :: method = ''
   integer :: nstate = 0
   real(8), allocatable :: w(:), f(:), mu(:,:)
   integer, allocatable :: ntr(:), off(:)
   integer :: ntot = 0
   integer, allocatable :: ifrom(:), ito(:)
   real(8), allocatable :: wt(:)
   character(len=1), allocatable :: spin(:)
   logical :: has_ct = .false.
   real(8), allocatable :: ct(:,:)
end type tddft_set_t

type :: polar_rec_t
   real(8) :: omega = 0.0d0, alpha(3,3) = 0.0d0, aiso = 0.0d0, aani = 0.0d0
   integer :: niter = 0
end type polar_rec_t

type :: engine_report_t
   logical :: energy_only = .false.
   real(8) :: energy = 0.0d0
   logical :: has_props = .false.
   real(8) :: dip_au(3) = 0.0d0
   real(8), allocatable :: q_mulliken(:), q_lowdin(:), q_resp(:), e_alpha(:), e_beta(:)
   logical :: have_resp = .false., restricted = .true.
   integer :: n_alpha = 0, n_beta = 0
   logical :: has_polar = .false.
   type(polar_rec_t) :: stat
   integer :: n_dyn = 0
   type(polar_rec_t) :: dyn(MAX_DYN)
   type(tddft_set_t) :: tddft(3)
   logical :: has_freq = .false.
   real(8), allocatable :: freq(:), inten(:), redmass(:), fc(:), disp(:,:,:)
   logical :: has_thermo = .false.
   integer :: th_mult = 1, th_sigma = 1
   real(8) :: th_T = 0.0d0, th_P = 0.0d0, th_mass = 0.0d0, th_zpe = 0.0d0, th_e = 0.0d0, th_h = 0.0d0, th_g = 0.0d0
   real(8) :: th_s = 0.0d0, th_cv = 0.0d0, th_s_el = 0.0d0, th_s_trans = 0.0d0, th_s_rot = 0.0d0, th_s_vib = 0.0d0
   real(8) :: th_e0 = 0.0d0
   logical :: has_dipder = .false.
   integer :: dipder_natoms = 0
   real(8), allocatable :: dmu(:,:)
end type engine_report_t

abstract interface
   subroutine report_hook_i(rep)
   import :: engine_report_t
   type(engine_report_t), intent(in) :: rep
   end subroutine report_hook_i
end interface
procedure(report_hook_i), pointer, save :: hook => null()

type(engine_report_t), save, target :: rpt

contains

subroutine report_reset(hook_in)
procedure(report_hook_i), pointer, intent(in) :: hook_in
type(engine_report_t) :: empty
rpt = empty
hook => hook_in
end subroutine report_reset

subroutine notify()
if (associated(hook)) call hook(rpt)
end subroutine notify

subroutine report_energy_only(energy)
real(8), intent(in) :: energy
rpt%energy_only = .true.
rpt%energy = energy
call notify()
end subroutine report_energy_only

subroutine report_properties(dip_au, q_mulliken, q_lowdin, have_resp, q_resp, restricted, e_alpha, n_alpha, e_beta, n_beta)
real(8), intent(in) :: dip_au(3), q_mulliken(:), q_lowdin(:), q_resp(:), e_alpha(:), e_beta(:)
logical, intent(in) :: have_resp, restricted
integer, intent(in) :: n_alpha, n_beta
rpt%has_props = .true.
rpt%dip_au = dip_au
rpt%q_mulliken = q_mulliken
rpt%q_lowdin = q_lowdin
rpt%q_resp = q_resp
rpt%have_resp = have_resp
rpt%restricted = restricted
rpt%e_alpha = e_alpha
rpt%e_beta = e_beta
rpt%n_alpha = n_alpha
rpt%n_beta = n_beta
call notify()
end subroutine report_properties

subroutine report_polar(alpha, aiso, aani, niter)
real(8), intent(in) :: alpha(3,3), aiso, aani
integer, intent(in) :: niter
rpt%has_polar = .true.
rpt%stat%alpha = alpha
rpt%stat%aiso = aiso
rpt%stat%aani = aani
rpt%stat%niter = niter
call notify()
end subroutine report_polar

subroutine report_polar_dynamic(omega, alpha, aiso, aani, niter)
real(8), intent(in) :: omega, alpha(3,3), aiso, aani
integer, intent(in) :: niter
if (rpt%n_dyn >= MAX_DYN) return
rpt%n_dyn = rpt%n_dyn + 1
rpt%dyn(rpt%n_dyn)%omega = omega
rpt%dyn(rpt%n_dyn)%alpha = alpha
rpt%dyn(rpt%n_dyn)%aiso = aiso
rpt%dyn(rpt%n_dyn)%aani = aani
rpt%dyn(rpt%n_dyn)%niter = niter
call notify()
end subroutine report_polar_dynamic

subroutine report_tddft_state(kind, istate, w, f, mu, n, ifrom, ito, wt, napply, method, spins, ct &
   )
integer, intent(in) :: kind, istate, n, ifrom(*), ito(*), napply
real(8), intent(in) :: w, f, mu(3), wt(*)
character(len=*), intent(in) :: method
character(len=1), intent(in), optional :: spins(*)
real(8), intent(in), optional :: ct(5)
type(tddft_set_t), pointer :: s
s => rpt%tddft(kind)
if (istate == 1 .or. .not. s%used) then
   s%used = .true.
   s%nstate = 0
   s%ntot = 0
   if (allocated(s%w)) deallocate(s%w, s%f, s%mu, s%ntr, s%off, s%ifrom, s%ito, s%wt)
   if (allocated(s%spin)) deallocate(s%spin)
   if (allocated(s%ct)) deallocate(s%ct)
   s%has_ct = .false.
   allocate(s%ct(5,16))
   allocate(s%w(16), s%f(16), s%mu(3,16), s%ntr(16), s%off(16), s%ifrom(64), s%ito(64), s%wt(64))
   if (present(spins)) allocate(s%spin(64))
endif
s%napply = napply
s%method = method
if (s%nstate + 1 > size(s%w)) call grow_states(s)
if (s%ntot + n > size(s%ifrom)) call grow_trans(s, s%ntot + n)
s%nstate = s%nstate + 1
s%w(s%nstate) = w ; s%f(s%nstate) = f ; s%mu(:,s%nstate) = mu
s%ct(:,s%nstate) = 0.0d0
if (present(ct)) then
   s%has_ct = .true. ; s%ct(:,s%nstate) = ct
endif
s%ntr(s%nstate) = n ; s%off(s%nstate) = s%ntot
s%ifrom(s%ntot+1:s%ntot+n) = ifrom(1:n)
s%ito(s%ntot+1:s%ntot+n) = ito(1:n)
s%wt(s%ntot+1:s%ntot+n) = wt(1:n)
if (present(spins)) s%spin(s%ntot+1:s%ntot+n) = spins(1:n)
s%ntot = s%ntot + n
call notify()
end subroutine report_tddft_state

subroutine grow_states(s)
type(tddft_set_t), intent(inout) :: s
integer :: n0
real(8), allocatable :: r1(:), r2(:), r3(:,:)
integer, allocatable :: i1(:), i2(:)
n0 = size(s%w)
allocate(r1(2*n0), r2(2*n0), r3(3,2*n0), i1(2*n0), i2(2*n0))
r1(1:n0) = s%w ; r2(1:n0) = s%f ; r3(:,1:n0) = s%mu ; i1(1:n0) = s%ntr ; i2(1:n0) = s%off
call move_alloc(r1, s%w) ; call move_alloc(r2, s%f) ; call move_alloc(r3, s%mu)
block
   real(8), allocatable :: r6(:,:)
   allocate(r6(5,2*n0)) ; r6(:,1:n0) = s%ct
   call move_alloc(r6, s%ct)
end block
call move_alloc(i1, s%ntr) ; call move_alloc(i2, s%off)
end subroutine grow_states

subroutine grow_trans(s, need)
type(tddft_set_t), intent(inout) :: s
integer, intent(in) :: need
integer :: n0, n1
integer, allocatable :: i1(:), i2(:)
real(8), allocatable :: r1(:)
character(len=1), allocatable :: c1(:)
n0 = size(s%ifrom)
n1 = max(2*n0, need)
allocate(i1(n1), i2(n1), r1(n1))
i1(1:n0) = s%ifrom ; i2(1:n0) = s%ito ; r1(1:n0) = s%wt
call move_alloc(i1, s%ifrom) ; call move_alloc(i2, s%ito) ; call move_alloc(r1, s%wt)
if (allocated(s%spin)) then
   allocate(c1(n1))
   c1(1:n0) = s%spin
   call move_alloc(c1, s%spin)
endif
end subroutine grow_trans

subroutine report_freq(freq, inten, redmass, fc, disp)
real(8), intent(in) :: freq(:), inten(:), redmass(:), fc(:), disp(:,:,:)
rpt%has_freq = .true.
rpt%freq = freq
rpt%inten = inten
rpt%redmass = redmass
rpt%fc = fc
rpt%disp = disp
call notify()
end subroutine report_freq

subroutine report_thermo(T, P, mult, sigma, mass, zpe, e_corr, h_corr, g_corr, s_tot, cv_tot, s_el, s_trans, s_rot, s_vib, e0)
real(8), intent(in) :: T, P, mass, zpe, e_corr, h_corr, g_corr, s_tot, cv_tot, s_el, s_trans, s_rot, s_vib, e0
integer, intent(in) :: mult, sigma
rpt%has_thermo = .true.
rpt%th_T = T ; rpt%th_P = P ; rpt%th_mult = mult ; rpt%th_sigma = sigma ; rpt%th_mass = mass
rpt%th_zpe = zpe ; rpt%th_e = e_corr ; rpt%th_h = h_corr ; rpt%th_g = g_corr
rpt%th_s = s_tot ; rpt%th_cv = cv_tot ; rpt%th_s_el = s_el ; rpt%th_s_trans = s_trans
rpt%th_s_rot = s_rot ; rpt%th_s_vib = s_vib ; rpt%th_e0 = e0
call notify()
end subroutine report_thermo

subroutine report_dipole_deriv(dmu, natoms)
integer, intent(in) :: natoms
real(8), intent(in) :: dmu(3, 3*natoms)
rpt%has_dipder = .true.
rpt%dipder_natoms = natoms
rpt%dmu = dmu
call notify()
end subroutine report_dipole_deriv

end module mod_engine_report
