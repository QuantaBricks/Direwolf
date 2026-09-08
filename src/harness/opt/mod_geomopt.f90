! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Geometry optimizer for the Direwolf harness.

module mod_geomopt
use mod_engine_input_types, only: engine_input_t
use mod_engineup_interface, only: EngineUp
use mod_lindh_hessian, only: lindh_cart_hessian
use mod_engine_results, only: write_results_block
use mod_engine_input_elements, only: elem_to_name
use mod_internal_coords, only: ic_set_t, ic_build, ic_rebuild_if_changed, ic_nprim, &
                                ic_values, ic_bmatrix, ic_ginv, ic_delta_q, &
                                ic_grad_to_internal, ic_back_transform, ic_diag_hessian
implicit none
private
public :: geomopt_run

real(8), parameter :: ANG2BOHR = 1.88972612456506d0
real(8), parameter :: BOHR2ANG = 1.0d0 / ANG2BOHR

real(8), parameter :: PROJ_SHIFT = 1.0d3

integer, parameter :: MAXRETRY = 3
real(8), parameter :: E_REJECT = 1.0d-4

contains

subroutine geomopt_run(spec, outfile, out_unit)
type(engine_input_t), intent(inout) :: spec
character(len=*), intent(in)        :: outfile
integer, intent(in)                 :: out_unit

integer :: nat, n, i, cyc, n_retry, ni
logical :: ric, ric_have_prev, topo_changed
real(8) :: t_gmax, t_grms, t_smax, t_srms, t_de
real(8) :: trust, trust_min, trust_max
real(8) :: E, E_old, dE, gmax, grms, smax, srms, sn, pred, pred_prev, ared, ratio
logical :: converged, have_prev, did_restart

real(8), allocatable :: g3(:,:), g3_prev(:,:)
real(8), allocatable :: grad(:), grad_old(:), delta(:), step_prev(:)
real(8), allocatable :: Hx(:,:)
real(8), allocatable :: Hq(:,:)
real(8), allocatable :: q_prev(:), gq_prev(:), q_cur(:), sq(:), gq_cur(:)
type(ic_set_t)       :: ic
character(len=256)   :: trajfile

integer, allocatable :: atomchg_a(:)
real(8), allocatable :: coord_a(:,:), force_out(:,:), MLcharge_out(:)
real(8), allocatable :: pc_charge_a(:), pc_coord_a(:,:)
real(8), allocatable :: dens_in_a(:,:), dens_in_b(:,:), dens_out_a(:,:), dens_out_b(:,:)
integer :: iconv
real(8) :: energy_out, econv, mem_grid_gb, mem_2e_gb

nat = spec%ncenters
n   = 3 * nat
ric = (trim(spec%opt_coord) == 'ric') .or. (trim(spec%opt_coord) == 'internal') &
      .or. (trim(spec%opt_coord) == 'redundant')

allocate(g3(3,nat), g3_prev(3,nat))
allocate(grad(n), grad_old(n), delta(n), step_prev(n))
allocate(atomchg_a(nat), coord_a(nat,3), force_out(nat,3), MLcharge_out(nat))
allocate(pc_charge_a(max(spec%npc,1)), pc_coord_a(max(spec%npc,1),3))
allocate(dens_in_a(0,0), dens_in_b(0,0))

atomchg_a = spec%atomchg
do i = 1, spec%npc
   pc_charge_a(i)  = spec%pc_q(i)
   pc_coord_a(i,1) = spec%pc_x(i)
   pc_coord_a(i,2) = spec%pc_y(i)
   pc_coord_a(i,3) = spec%pc_z(i)
enddo

if (spec%opt_conv_level >= 1) then
   t_gmax = 1.5d-5 ; t_grms = 1.0d-5 ; t_smax = 6.0d-5 ; t_srms = 4.0d-5 ; t_de = 1.0d-8
else
   t_gmax = 4.5d-4 ; t_grms = 3.0d-4 ; t_smax = 1.8d-3 ; t_srms = 1.2d-3 ; t_de = 1.0d-6
endif

trust     = spec%opt_trust
trust_min = 1.0d-2
trust_max = max(spec%opt_trust, 3.0d-1) * 2.0d0

trajfile = trajname(outfile)
did_restart = .false.
if (spec%opt_restart) call traj_read_last(trajfile, nat, g3, did_restart)
if (.not. did_restart) then
   do i = 1, nat
      g3(1,i) = spec%x(i) ; g3(2,i) = spec%y(i) ; g3(3,i) = spec%z(i)
   enddo
endif
g3 = g3 * ANG2BOHR

if (did_restart) then
   write(out_unit,'(A)') '[OPTRESTART]  resumed from '//trim(trajfile)//' (model Hessian rebuilt)'
else
   call traj_truncate(trajfile)
endif

if (ric) then
   call ic_build(ic, nat, spec%atomchg, g3)
   ni = ic_nprim(ic)
   if (ni < 1) then
      write(out_unit,'(A)') '[OPT] no internal coordinates generated - falling back to cart'
      ric = .false.
   endif
endif

if (ric) then
   allocate(Hq(ni,ni), q_prev(ni), gq_prev(ni), q_cur(ni), sq(ni), gq_cur(ni))
   call ic_diag_hessian(ic, Hq)
   write(out_unit,'(A,I0,A,I0,A,I0,A,F5.3)') '[OPT]  ric (', ni, &
        ' primitives) / diagonal model Hessian  maxcyc=', spec%opt_maxcyc, &
        '  conv=', spec%opt_conv_level, '  trust0=', trust
else
   allocate(Hx(n,n))
   call lindh_cart_hessian(nat, spec%atomchg, g3, Hx)
   write(out_unit,'(A,I0,A,I0,A,F5.3)') '[OPT]  cart / Lindh model Hessian  maxcyc=', &
        spec%opt_maxcyc, '  conv=', spec%opt_conv_level, '  trust0=', trust
endif

have_prev     = .false.
ric_have_prev = .false.
converged     = .false.
n_retry       = 0
E_old = 0.0d0 ; pred_prev = 0.0d0
step_prev = 0.0d0 ; g3_prev = 0.0d0 ; grad_old = 0.0d0

do cyc = 1, spec%opt_maxcyc

   do i = 1, nat
      coord_a(i,1) = g3(1,i) * BOHR2ANG
      coord_a(i,2) = g3(2,i) * BOHR2ANG
      coord_a(i,3) = g3(3,i) * BOHR2ANG
   enddo

   call EngineUp(spec%ncenters, spec%imult, spec%icharge, spec%functional, &
                 coord_a, atomchg_a, spec%baselabel, spec%ecplabel, &
                 spec%atom_basis, spec%atom_ecp, &
                 trim(spec%J), trim(spec%K), trim(spec%ri_aux_basis), spec%spherical, spec%harris_guess, &
                 .true., spec%vv10_nonself, spec%mem_cap_gb, &
                 .false., spec%n_threads, &
                 spec%npc, pc_charge_a, pc_coord_a, &
                 spec%molden_write, trim(spec%molden_file), spec%molden_read, trim(spec%molden_read_file), &
                 spec%cosmo_on, spec%cosmo_epsilon, spec%cosmo_radii_scale, spec%cosmo_avg_area, spec%cosmo_sigma_rav, &
                 trim(spec%cosmo_cavity_type), spec%cosmo_rsolv, spec%cosmo_ks_nseg, spec%cosmo_ks_nface, &
                 trim(spec%cosmo_sigma_profile_file), &
                 spec%cosmo_smd, trim(spec%cosmo_solvent), &
                 force_out, energy_out, MLcharge_out, iconv, econv, &
                 mem_grid_gb, mem_2e_gb, max(spec%scf_conv_level,1), trim(spec%basedir), &
                 dens_in_a, dens_in_b, dens_out_a, dens_out_b)

   if (iconv /= 1) then
      write(out_unit,'(A,I0,A)') '[OPTFAILED]  SCF did not converge at cycle ', cyc, &
           ' - aborting optimization'
      call finish(spec, g3, nat, out_unit, .false., cyc, energy_out, force_out, iconv)
      return
   endif

   E    = energy_out
   grad = reshape(transpose(force_out), [n])
   gmax = maxval(abs(grad))
   grms = sqrt(sum(grad*grad) / n)

   if (have_prev) then
      ared = E - E_old
      if (pred_prev < -1.0d-14) then
         ratio = ared / pred_prev
      else
         ratio = -1.0d0
      endif
      if (ratio < 0.25d0) then
         trust = max(0.25d0 * sn, trust_min)
      else if (ratio > 0.75d0 .and. sn > 0.99d0*trust) then
         trust = min(2.0d0 * trust, trust_max)
      endif
   endif

   call write_cycle(out_unit, cyc, E, have_prev, E_old, gmax, grms, have_prev, smax, srms, trust, iconv)
   call traj_append(trajfile, nat, spec%atomchg, g3, cyc, E, gmax)
   if (spec%chk_write) call chk_write_safe(spec, dens_out_a, dens_out_b)

   if (have_prev) then
      dE = E - E_old
      converged = (gmax < t_gmax) .and. (grms < t_grms) .and. &
                  ( ((smax < t_smax) .and. (srms < t_srms)) .or. (abs(dE) < t_de) )
   else
      converged = (gmax < t_gmax) .and. (grms < t_grms)
   endif
   if (converged) then
      write(out_unit,'(A)')       '[OPTCONVERGED]'
      write(out_unit,'(A,I0)')    '  cycles = ', cyc
      write(out_unit,'(A,F20.9)') '  final E (Hartree) = ', E
      call finish(spec, g3, nat, out_unit, .true., cyc, E, force_out, iconv)
      call profile_tail(out_unit)
      return
   endif

   if (have_prev .and. (E - E_old) > E_REJECT .and. n_retry < MAXRETRY) then
      n_retry = n_retry + 1
      write(out_unit,'(A,I0,A,ES9.2,A)') '  [reject #', n_retry, ']  dE=+', E-E_old, &
           ' - backtracking with smaller trust'
      g3    = g3_prev
      grad  = grad_old
      E     = E_old
      trust = max(0.25d0 * sn, trust_min)
      if (ric) then
         call ric_step(ic, nat, ni, g3, grad, Hq, trust, delta, pred, gq_cur, out_unit)
      else
         call rfo_step(n, nat, Hx, grad, g3, trust, delta, pred, out_unit)
      endif
      sn = sqrt(sum(delta*delta))
      step_prev = delta
      pred_prev = pred
      g3 = g3_prev + reshape(delta, [3,nat])
      cycle
   endif
   n_retry = 0

   if (ric) then
      call ic_rebuild_if_changed(ic, nat, spec%atomchg, g3, topo_changed)
      if (topo_changed) then
         ni = ic_nprim(ic)
         deallocate(Hq, q_prev, gq_prev, q_cur, sq, gq_cur)
         allocate(Hq(ni,ni), q_prev(ni), gq_prev(ni), q_cur(ni), sq(ni), gq_cur(ni))
         call ic_diag_hessian(ic, Hq)
         ric_have_prev = .false.
         write(out_unit,'(A,I0,A)') '  [ric] connectivity changed - rebuilt ', ni, &
              ' primitives, Hessian reset'
      endif
      call ic_grad_to_internal(ic, g3, grad, gq_cur)
      if (ric_have_prev) then
         call ic_values(ic, g3, q_cur)
         call ic_delta_q(ic, q_cur, q_prev, sq)
         call bfgs_update(ni, Hq, sq, gq_cur - gq_prev)
      endif
      call ric_step(ic, nat, ni, g3, grad, Hq, trust, delta, pred, gq_cur, out_unit)
      call ic_values(ic, g3, q_prev)
      gq_prev = gq_cur
      ric_have_prev = .true.
   else
      if (have_prev) call bfgs_update(n, Hx, step_prev, grad - grad_old)
      call rfo_step(n, nat, Hx, grad, g3, trust, delta, pred, out_unit)
   endif
   sn = sqrt(sum(delta*delta))

   g3_prev   = g3
   E_old     = E
   grad_old  = grad
   step_prev = delta
   pred_prev = pred
   smax = maxval(abs(delta))
   srms = sqrt(sum(delta*delta) / n)

   g3 = g3 + reshape(delta, [3,nat])
   have_prev = .true.

   call move_alloc(dens_out_a, dens_in_a)
   call move_alloc(dens_out_b, dens_in_b)
enddo

write(out_unit,'(A,I0,A)') '[OPTFAILED]  not converged in ', spec%opt_maxcyc, ' cycles'
call finish(spec, g3, nat, out_unit, .false., spec%opt_maxcyc, E, force_out, iconv)

end subroutine geomopt_run

subroutine rfo_step(n, nat, Hx, grad, coord3, trust, delta, pred, out_unit)
integer, intent(in)  :: n, nat, out_unit
real(8), intent(in)  :: Hx(n,n), grad(n), coord3(3,nat), trust
real(8), intent(out) :: delta(n), pred
real(8), allocatable :: P(:,:), T1(:,:), Hp(:,:), gp(:), tr(:,:), Hd(:)
integer :: ntr, i, j, k
real(8) :: sn

allocate(P(n,n), T1(n,n), Hp(n,n), gp(n), tr(n,6), Hd(n))
call build_tr_basis(nat, coord3, tr, ntr)

P = 0.0d0
do i = 1, n
   P(i,i) = 1.0d0
enddo
do k = 1, ntr
   do j = 1, n
      do i = 1, n
         P(i,j) = P(i,j) - tr(i,k)*tr(j,k)
      enddo
   enddo
enddo

call dgemm('N','N', n, n, n, 1.0d0, Hx, n, P, n, 0.0d0, T1, n)
call dgemm('N','N', n, n, n, 1.0d0, P, n, T1, n, 0.0d0, Hp, n)
do k = 1, ntr
   do j = 1, n
      do i = 1, n
         Hp(i,j) = Hp(i,j) + PROJ_SHIFT * tr(i,k) * tr(j,k)
      enddo
   enddo
enddo
call dgemv('N', n, n, 1.0d0, P, n, grad, 1, 0.0d0, gp, 1)

call solve_rfo(n, Hp, gp, delta, out_unit)

sn = sqrt(sum(delta*delta))
if (sn > trust .and. sn > 0.0d0) delta = delta * (trust / sn)

call dgemv('N', n, n, 1.0d0, Hp, n, delta, 1, 0.0d0, Hd, 1)
pred = dot_product(gp, delta) + 0.5d0 * dot_product(delta, Hd)
end subroutine rfo_step

subroutine ric_step(ic, nat, ni, coord3, gx, Hq, trust, dx, pred, gq_out, out_unit)
type(ic_set_t), intent(in) :: ic
integer, intent(in)  :: nat, ni, out_unit
real(8), intent(in)  :: coord3(3,nat), gx(3*nat), Hq(ni,ni), trust
real(8), intent(out) :: dx(3*nat), pred, gq_out(ni)
real(8), allocatable :: B(:,:), Ginv(:,:), G(:,:), P(:,:), T1(:,:), Hp(:,:)
real(8), allocatable :: gq(:), gp(:), dq(:), dqp(:), Hd(:), bg(:)
integer :: ndof, i, j
real(8) :: sn, dxn
logical :: ok

ndof = 3*nat
allocate(B(ni,ndof), Ginv(ni,ni), G(ni,ni), P(ni,ni), T1(ni,ni), Hp(ni,ni))
allocate(gq(ni), gp(ni), dq(ni), dqp(ni), Hd(ni), bg(ni))

call ic_bmatrix(ic, coord3, B)
call ic_ginv(B, ni, ndof, Ginv)

call dgemv('N', ni, ndof, 1.0d0, B, ni, gx, 1, 0.0d0, bg, 1)
call dgemv('N', ni, ni, 1.0d0, Ginv, ni, bg, 1, 0.0d0, gq, 1)
gq_out = gq

call dgemm('N','T', ni, ni, ndof, 1.0d0, B, ni, B, ni, 0.0d0, G, ni)
call dgemm('N','N', ni, ni, ni, 1.0d0, Ginv, ni, G, ni, 0.0d0, P, ni)

call dgemm('N','N', ni, ni, ni, 1.0d0, Hq, ni, P, ni, 0.0d0, T1, ni)
call dgemm('N','N', ni, ni, ni, 1.0d0, P, ni, T1, ni, 0.0d0, Hp, ni)
do j = 1, ni
   do i = 1, ni
      Hp(i,j) = Hp(i,j) + PROJ_SHIFT * (merge(1.0d0, 0.0d0, i == j) - P(i,j))
   enddo
enddo
call dgemv('N', ni, ni, 1.0d0, P, ni, gq, 1, 0.0d0, gp, 1)

call solve_rfo(ni, Hp, gp, dq, out_unit)
call dgemv('N', ni, ni, 1.0d0, P, ni, dq, 1, 0.0d0, dqp, 1)
dq = dqp

sn = sqrt(sum(dq*dq))
if (sn > trust .and. sn > 0.0d0) dq = dq * (trust / sn)

call dgemv('N', ni, ni, 1.0d0, Hp, ni, dq, 1, 0.0d0, Hd, 1)
pred = dot_product(gp, dq) + 0.5d0 * dot_product(dq, Hd)

call ic_back_transform(ic, coord3, dq, dx, ok)
if (.not. ok) write(out_unit,'(A)') &
     '  [ric] back-transform not fully converged - using linear estimate'

dxn = sqrt(sum(dx*dx))
if (dxn > 2.0d0*trust .and. dxn > 0.0d0) dx = dx * (2.0d0*trust / dxn)
end subroutine ric_step

subroutine solve_rfo(n, H, g, step, out_unit)
integer, intent(in)  :: n, out_unit
real(8), intent(in)  :: H(n,n), g(n)
real(8), intent(out) :: step(n)
real(8), allocatable :: aug(:,:), w(:), work(:)
integer :: lwork, info
real(8) :: vlast
allocate(aug(n+1,n+1), w(n+1))
aug = 0.0d0
aug(1:n,1:n) = H
aug(1:n,n+1) = g
aug(n+1,1:n) = g
lwork = max(1, 3*(n+1) + (n+1)*(n+1))
allocate(work(lwork))
call dsyev('V', 'U', n+1, aug, n+1, w, work, lwork, info)
if (info /= 0) then
   write(out_unit,'(A,I0)') '  [rfo] dsyev info=', info
   step = -g
   return
endif
vlast = aug(n+1,1)
if (abs(vlast) > 1.0d-8) then
   step = aug(1:n,1) / vlast
else
   step = -g
endif
end subroutine solve_rfo

subroutine build_tr_basis(nat, coord3, tr, ntr)
integer, intent(in)  :: nat
real(8), intent(in)  :: coord3(3,nat)
real(8), intent(out) :: tr(3*nat,6)
integer, intent(out) :: ntr
real(8) :: com(3), raw(3*nat,6), d(3), nrm, dotv
integer :: i, ax, k, m

com = 0.0d0
do i = 1, nat
   com = com + coord3(:,i)
enddo
com = com / real(nat,8)

raw = 0.0d0
do ax = 1, 3
   do i = 1, nat
      raw(3*(i-1)+ax, ax) = 1.0d0
   enddo
enddo
do ax = 1, 3
   do i = 1, nat
      d = coord3(:,i) - com
      select case (ax)
      case (1) ; raw(3*(i-1)+2, 3+ax) = -d(3) ; raw(3*(i-1)+3, 3+ax) =  d(2)
      case (2) ; raw(3*(i-1)+1, 3+ax) =  d(3) ; raw(3*(i-1)+3, 3+ax) = -d(1)
      case (3) ; raw(3*(i-1)+1, 3+ax) = -d(2) ; raw(3*(i-1)+2, 3+ax) =  d(1)
      end select
   enddo
enddo

ntr = 0
do k = 1, 6
   do m = 1, ntr
      dotv = dot_product(raw(:,k), tr(:,m))
      raw(:,k) = raw(:,k) - dotv * tr(:,m)
   enddo
   nrm = sqrt(sum(raw(:,k)*raw(:,k)))
   if (nrm > 1.0d-8) then
      ntr = ntr + 1
      tr(:,ntr) = raw(:,k) / nrm
   endif
enddo
end subroutine build_tr_basis

subroutine bfgs_update(n, H, s, y)
integer, intent(in)    :: n
real(8), intent(inout) :: H(n,n)
real(8), intent(in)    :: s(n), y(n)
real(8), allocatable   :: Hs(:)
real(8) :: sy, sHs
integer :: i, j
sy = dot_product(s, y)
if (sy <= 1.0d-10) return
allocate(Hs(n))
call dgemv('N', n, n, 1.0d0, H, n, s, 1, 0.0d0, Hs, 1)
sHs = dot_product(s, Hs)
if (sHs <= 1.0d-12) return
do j = 1, n
   do i = 1, n
      H(i,j) = H(i,j) - Hs(i)*Hs(j)/sHs + y(i)*y(j)/sy
   enddo
enddo
end subroutine bfgs_update

subroutine write_cycle(u, cyc, E, have_dE, E_old, gmax, grms, have_step, smax, srms, trust, iconv)
integer, intent(in) :: u, cyc, iconv
real(8), intent(in) :: E, E_old, gmax, grms, smax, srms, trust
logical, intent(in) :: have_dE, have_step
write(u,'(A,I0,A)') '[OPTCYCLE ', cyc, ']'
if (have_dE) then
   write(u,'(A,F18.9,A,ES11.3)') '  E = ', E, '   dE = ', E - E_old
else
   write(u,'(A,F18.9,A)')        '  E = ', E, '   dE =        --'
endif
write(u,'(A,ES11.3,A,ES11.3)') '  gmax = ', gmax, '   grms = ', grms
if (have_step) write(u,'(A,ES11.3,A,ES11.3)') '  smax = ', smax, '   srms = ', srms
write(u,'(A,F6.3,A,I0,A)') '  trust = ', trust, '   scf = ', iconv, ' (1=converged)'
end subroutine write_cycle

subroutine finish(spec, g3, nat, u, ok, ncyc, E, force_out, iconv)
type(engine_input_t), intent(inout) :: spec
integer, intent(in) :: nat, u, ncyc, iconv
real(8), intent(in) :: g3(3,nat), E, force_out(nat,3)
logical, intent(in) :: ok
integer :: i
do i = 1, nat
   spec%x(i) = g3(1,i) * BOHR2ANG
   spec%y(i) = g3(2,i) * BOHR2ANG
   spec%z(i) = g3(3,i) * BOHR2ANG
enddo
write(u,'(A)') ''
call write_results_block(u, nat, E, force_out, iconv, spec%calc_force)
write(u,'(A)') ''
write(u,'(A)') '[FINAL GEOMETRY] (angstrom)'
do i = 1, nat
   write(u,'(2X,A2,3F16.9)') elem_to_name(spec%atomchg(i)), &
        g3(1,i)*BOHR2ANG, g3(2,i)*BOHR2ANG, g3(3,i)*BOHR2ANG
enddo
if (.not. ok) write(u,'(A)') '[OPT] geometry above is the last step, NOT converged'
end subroutine finish

subroutine profile_tail(u)
use mod_profile, only: prof_report
integer, intent(in) :: u
if (u == 6) then
   write(u,'(A)') '[PROFILE] (final optimization cycle)'
   call prof_report(u)
   write(u,'(A)') '[RESULTSEND]'
endif
end subroutine profile_tail

function trajname(outfile) result(tf)
character(len=*), intent(in) :: outfile
character(len=256) :: tf
integer :: l
l = len_trim(outfile)
if (l > 4 .and. outfile(max(1,l-3):l) == '.out') then
   tf = outfile(1:l-4) // '.opt.xyz'
else
   tf = trim(outfile) // '.opt.xyz'
endif
end function trajname

subroutine traj_truncate(tf)
character(len=*), intent(in) :: tf
integer :: u, ios
open(newunit=u, file=trim(tf), status='replace', action='write', iostat=ios)
if (ios == 0) close(u)
end subroutine traj_truncate

subroutine traj_append(tf, nat, Z, g3, cyc, E, gmax)
character(len=*), intent(in) :: tf
integer, intent(in) :: nat, Z(nat), cyc
real(8), intent(in) :: g3(3,nat), E, gmax
integer :: u, i, ios
open(newunit=u, file=trim(tf), status='unknown', position='append', action='write', iostat=ios)
if (ios /= 0) return
write(u,'(I0)') nat
write(u,'(A,I0,A,F18.9,A,ES11.3)') 'cycle ', cyc, '  E = ', E, '  gmax = ', gmax
do i = 1, nat
   write(u,'(A2,3F16.9)') elem_to_name(Z(i)), g3(1,i)*BOHR2ANG, g3(2,i)*BOHR2ANG, g3(3,i)*BOHR2ANG
enddo
close(u)
end subroutine traj_append

subroutine traj_read_last(tf, nat, g3, ok)
character(len=*), intent(in)  :: tf
integer, intent(in)           :: nat
real(8), intent(out)          :: g3(3,nat)
logical, intent(out)          :: ok
integer :: u, ios, fnat, i
character(len=2)   :: sym
character(len=512) :: line
real(8) :: last(3,nat), x, y, z
ok = .false.
open(newunit=u, file=trim(tf), status='old', action='read', iostat=ios)
if (ios /= 0) return
do
   read(u,'(A)', iostat=ios) line
   if (ios /= 0) exit
   if (len_trim(line) == 0) cycle
   read(line,*,iostat=ios) fnat
   if (ios /= 0 .or. fnat /= nat) then
      close(u) ; return
   endif
   read(u,'(A)', iostat=ios) line
   if (ios /= 0) exit
   do i = 1, nat
      read(u,*,iostat=ios) sym, x, y, z
      if (ios /= 0) then
         close(u) ; return
      endif
      last(1,i) = x ; last(2,i) = y ; last(3,i) = z
   enddo
   g3 = last
   ok = .true.
enddo
close(u)
end subroutine traj_read_last

subroutine chk_write_safe(spec, da, db)
use mod_checkpoint, only: checkpoint_write
type(engine_input_t), intent(in) :: spec
real(8), allocatable, intent(in) :: da(:,:), db(:,:)
if (allocated(da) .and. allocated(db)) call checkpoint_write(trim(spec%chk_file), da, db)
end subroutine chk_write_safe

end module mod_geomopt
