! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Geometry optimizer for the Direwolf harness.

module mod_geomopt
use mod_engine_input_types, only: engine_input_t
use mod_engineup_interface, only: EngineUp
use mod_hessian_freq, only: freq_unit, freq_mult, freq_base, freq_thermo
use mod_integrals, only: freq_request, hessian_return, cosx_force_loose
use GRID_info, only: force_dense_allow
use mod_lindh_hessian, only: lindh_cart_hessian
use mod_engine_results, only: write_results_block
use mod_engine_input_elements, only: elem_to_name
use mod_scf_history, only: scf_hist_niter, scf_hist_last_prms, scf_hist_last_de
use mod_internal_coords, only: ic_set_t, ic_build, ic_rebuild_if_changed, ic_nprim, &
                                ic_values, ic_bmatrix, ic_ginv, ic_delta_q, &
                                ic_grad_to_internal, ic_back_transform, ic_lindh_hessian, &
                                ic_write, &
                                ic_rank
use mod_internal_coords, only: ic_n_dropped
use MOL_info, only: engine_quiet
use mod_opt_ts, only: ts_on, ts_nneg, ts_vprev, TS_MAX_RECALC, ts_start, ts_forget_mode, &
                      ts_endgame_refresh_due, ts_lost_mode_recalc_due, solve_prfo, bofill_update
use mod_opt_hessian, only: ANG2BOHR, BOHR2ANG, read_ext_hessian, ric_model_hessian, &
                           cart_hessian_to_internal, ric_grad_curv, hessian_from_xtb, hessian_analytic, &
                           ric_hessian_to_cart, write_hessian_file, remove_file
use mod_opt_diis, only: gediis_solve
use mod_opt_constraints, only: cons_init, cons_active, cons_has_frozen, cons_basis, cons_project_grad, cons_correct, &
                                orth_cols, cons_max_dev
use mod_opt_min, only: regularize_hessian, ric_floor_hessian, bfgs_update, solve_rfo
implicit none
private
public :: geomopt_run

real(8), parameter :: PROJ_SHIFT = 1.0d3

real(8), parameter :: CONS_DEV_TOL = 2.0d-3
integer, parameter :: MAXRETRY = 3
real(8), parameter :: E_REJECT = 1.0d-4
real(8), parameter :: TRUST_SHRINK = 0.5d0
real(8), parameter :: E_REJECT_NOISE = 1.0d-5
real(8), parameter :: E_REJECT_RATIO = -1.0d0
logical, parameter :: E_REJECT_BY_RATIO = .false.

contains

subroutine geomopt_run(spec, outfile, out_unit)
type(engine_input_t), intent(inout) :: spec
character(len=*), intent(in)        :: outfile
integer, intent(in)                 :: out_unit

integer :: nat, n, i, cyc, n_retry, ni, nrank, nexp
logical :: ext_hess_ok, hess_diag_on
character(32) :: hd_env
character(len=260) :: hessfile
logical :: hess_restored, hess_exist
integer :: hd_stat
real(8), allocatable :: Hext(:,:)
logical :: ric, ric_have_prev, topo_changed
real(8) :: t_gmax, t_grms, t_smax, t_srms, t_de
logical :: scf_adapt
logical :: hess_fresh, want_hess_now, hok
real(8), allocatable :: Hnew(:,:)
real(8), allocatable :: Hcarry(:,:), Bold(:,:), Tcarry(:,:)
character(len=8) :: hsrc0, hsrc_re
integer :: hess_every, n_hess_calls
integer :: ts_nrecalc
logical :: ts_multi_done
logical :: hq_pending
integer :: scf_lvl, scf_lvl_prev
logical :: lvl_switch
integer :: ged_m, ged_nh, ged_i, ged_ios
character(len=8) :: ged_env
real(8), allocatable :: ged_q(:,:), ged_g(:,:), ged_E(:), ged_dQ(:,:), ged_c(:), ged_qoff(:), ged_gstar(:)
real(8) :: ged_Estar
logical :: ged_ok, ged_used
real(8), parameter :: SCF_ADAPT_GMAX = 1.0d-3
logical, parameter :: COSX_FORCE_LOOSE_EARLY = .false.
real(8) :: trust, trust_min, trust_max
real(8) :: E, E_old, dE, gmax, grms, smax, srms, sn, pred, pred_prev, ared, ratio
logical :: converged, have_prev, did_restart

real(8), allocatable :: g3(:,:), g3_prev(:,:)
real(8), allocatable :: grad(:), grad_old(:), delta(:), step_prev(:), gproj(:)
real(8), allocatable :: Hx(:,:)
real(8), allocatable :: Hq(:,:)
real(8), allocatable :: q_prev(:), gq_prev(:), q_cur(:), sq(:), gq_cur(:)
real(8), allocatable :: Bmat(:,:), Ginvm(:,:), bgv(:)
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
call ts_start(spec%opt_ts)
ric = (trim(spec%opt_coord) == 'ric') .or. (trim(spec%opt_coord) == 'internal') &
      .or. (trim(spec%opt_coord) == 'redundant')

engine_quiet = .true.

force_dense_allow = .false.

allocate(g3(3,nat), g3_prev(3,nat))
allocate(grad(n), grad_old(n), delta(n), step_prev(n), gproj(n))
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
else if (spec%opt_conv_level < 0) then
   t_gmax = 2.0d-3 ; t_grms = 5.0d-4 ; t_smax = 1.0d-2 ; t_srms = 7.0d-3 ; t_de = 3.0d-5
else
   t_gmax = 4.5d-4 ; t_grms = 3.0d-4 ; t_smax = 1.8d-3 ; t_srms = 1.2d-3 ; t_de = 1.0d-6
endif
if (spec%opt_gmax > 0.0d0) t_gmax = spec%opt_gmax
if (spec%opt_grms > 0.0d0) t_grms = spec%opt_grms
if (spec%opt_smax > 0.0d0) t_smax = spec%opt_smax
if (spec%opt_srms > 0.0d0) t_srms = spec%opt_srms
if (spec%opt_de   > 0.0d0) t_de   = spec%opt_de
write(out_unit,'(A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2)') '[OPT]  convergence: gmax<', t_gmax, &
     ' grms<', t_grms, ' and (smax<', t_smax, ' srms<', t_srms, ' or |dE|<', t_de

scf_adapt = .false.
if (spec%scf_conv_level .lt. 1) then
   spec%scf_conv_level = 1
   write(out_unit,'(A)') &
      '[OPT]  SCF convergence raised to fine for the optimization (gradient quality)'
   scf_adapt = (spec%opt_conv_level .lt. 1)
   if (scf_adapt) write(out_unit,'(A,ES8.1)') &
      '[OPT]  early steps at regular SCF convergence while previous gmax >', SCF_ADAPT_GMAX
endif
if (spec%opt_conv_level .ge. 1 .and. spec%scf_conv_level .lt. 2) then
   spec%scf_conv_level = 2
   write(out_unit,'(A)') &
      '[OPT]  opt_conv=tight -> SCF convergence raised to tight'
endif

trust     = spec%opt_trust
trust_min = 5.0d-2
trust_max = max(spec%opt_trust, 3.0d-1) * 2.0d0
ts_nrecalc = 0
hq_pending = .false.
ts_multi_done = .false.

trajfile = trajname(outfile)
did_restart = .false.
if (spec%opt_restart) call traj_read_last(trajfile, nat, g3, did_restart)
if (.not. did_restart) then
   do i = 1, nat
      g3(1,i) = spec%x(i) ; g3(2,i) = spec%y(i) ; g3(3,i) = spec%z(i)
   enddo
endif
g3 = g3 * ANG2BOHR
block
   logical :: cons_ok
   call cons_init(spec%opt_freeze, spec%opt_constraints, nat, g3, out_unit, cons_ok)
   if (.not. cons_ok) then
      call finish(spec, g3, nat, out_unit, .false., 0, 0.0d0, force_out, 0)
      return
   endif
end block

if (did_restart) then
   write(out_unit,'(A)') '[OPTRESTART]  resumed from '//trim(trajfile)
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
   nexp = 3*nat - merge(5, 6, is_linear(nat, g3))
   if (nat < 3) nexp = max(3*nat - 5, 1)
   nrank = ic_rank(ic, g3)
   if (nrank < nexp) then
      write(out_unit,'(A,I0,A,I0,A)') &
         '[OPT] internal coordinates are incomplete (rank ', nrank, ' of ', nexp, &
         ' internal degrees of freedom) - falling back to cart'
      ric = .false.
   endif
   if (spec%opt_write_ric) call dump_ric(ic, g3, outfile, 'initial')
   if (ic_n_dropped .gt. 0) write(out_unit,'(A,I0,A)') &
      '[OPT] WARNING: ', ic_n_dropped, ' primitive(s)/bond(s) discarded for lack of room'// &
      ' - the coordinate set may be incomplete'
endif

call get_environment_variable("DIREWOLF_OPT_HESS_DIAG", hd_env, status=hd_stat)
hess_diag_on = (hd_stat == 0)
ged_m = 0
call get_environment_variable("DIREWOLF_OPT_GEDIIS", ged_env)
if (len_trim(ged_env) > 0) read(ged_env, *, iostat=ged_ios) ged_m
ged_m = max(0, min(ged_m, 12))
ged_nh = 0
hsrc0 = lower8(spec%opt_hess)
if (len_trim(hsrc0) == 0) hsrc0 = 'model'
if (hsrc0 == 'model' .and. len_trim(spec%opt_hessian_file) > 0) hsrc0 = 'read'
if (hsrc0 == 'model' .and. ts_on) hsrc0 = 'calc'
hessfile = trajfile(1:len_trim(trajfile)-4)//'.hess'
hess_restored = .false.
if (did_restart) then
   inquire(file=trim(hessfile), exist=hess_exist)
   if (hess_exist .and. hsrc0 == 'model' .and. .not. ts_on) hess_restored = .true.
else
   call remove_file(hessfile)
endif
hess_every = spec%opt_hess_every
if (hess_every <= 0 .and. ts_on) hess_every = spec%opt_ts_hess_every
hsrc_re = lower8(spec%opt_hess_recalc)
if (len_trim(hsrc_re) == 0) hsrc_re = merge('xtb ', 'calc', hsrc0 == 'xtb')
if (hsrc0 /= 'model' .and. hsrc0 /= 'calc' .and. hsrc0 /= 'read' .and. hsrc0 /= 'xtb') then
   write(out_unit,'(A)') "[OPT]  ERROR: opt_hess = '"//trim(spec%opt_hess)// &
        "' - expected 'model', 'calc', 'read' or 'xtb'"
   call finish(spec, g3, nat, out_unit, .false., 0, 0.0d0, force_out, 0)
   return
endif
if (hsrc_re /= 'calc' .and. hsrc_re /= 'xtb') then
   write(out_unit,'(A)') "[OPT]  ERROR: opt_hess_recalc = '"//trim(spec%opt_hess_recalc)// &
        "' - expected 'calc' or 'xtb'"
   call finish(spec, g3, nat, out_unit, .false., 0, 0.0d0, force_out, 0)
   return
endif
if (hsrc0 == 'read' .and. len_trim(spec%opt_hessian_file) == 0) then
   write(out_unit,'(A)') "[OPT]  ERROR: opt_hess = 'read' needs opt_hessian_file"
   call finish(spec, g3, nat, out_unit, .false., 0, 0.0d0, force_out, 0)
   return
endif
n_hess_calls = 0
if (hess_every > 0) write(out_unit,'(A,I0,A)') '[OPT]  Hessian recomputed every ', hess_every, &
     ' cycle(s) ('//trim(hsrc_re)//')'
ext_hess_ok = .false.
if (hess_restored) then
   allocate(Hext(n,n))
   call read_ext_hessian(trim(hessfile), n, Hext, ext_hess_ok, out_unit, raw=(ts_on .or. ric))
   if (ext_hess_ok) then
      write(out_unit,'(A)') '[OPT]  initial Hessian restored from '//trim(hessfile)
   else
      deallocate(Hext)
   endif
endif
select case (hsrc0)
case ('read')
   allocate(Hext(n,n))
   call read_ext_hessian(trim(spec%opt_hessian_file), n, Hext, ext_hess_ok, out_unit, raw=(ts_on .or. ric))
   if (ext_hess_ok) then
      write(out_unit,'(A)') '[OPT]  initial Hessian read from '// &
           trim(spec%opt_hessian_file)
   else
      deallocate(Hext)
   endif
case ('calc', 'xtb')
   write(out_unit,'(A)') '[OPT]  initial Hessian: '//trim(merge('analytic', 'xtb     ', hsrc0 == 'calc'))// &
        ' at the starting geometry'
   allocate(Hext(n,n))
   call hess_get(hsrc0, Hext, ext_hess_ok)
   if (.not. ext_hess_ok) deallocate(Hext)
end select
if (ts_on .and. .not. ext_hess_ok) then
   write(out_unit,'(A)') '[OPT]  the Hessian could not be obtained - transition-state search aborted'
   call finish(spec, g3, nat, out_unit, .false., 0, 0.0d0, force_out, 0)
   return
endif
if (ts_on) write(out_unit,'(A)') '[OPT]  TS search (P-RFO), start Hessian: '//trim(hsrc0)
if (len_trim(hsrc0) > 0 .and. hsrc0 /= 'model' .and. .not. ext_hess_ok) &
   write(out_unit,'(A)') '[OPT]  WARNING: start Hessian ('//trim(hsrc0)//') not available - Lindh model used'

if (ric) then
   allocate(Hq(ni,ni), q_prev(ni), gq_prev(ni), q_cur(ni), sq(ni), gq_cur(ni))
   allocate(bgv(ni))
   call ric_model_hessian(ic, nat, spec%atomchg, g3, ni, Hq, out_unit)
   hq_pending = ext_hess_ok
   write(out_unit,'(A,I0,A,A,A,I0,A,I0,A,F5.3)') '[OPT]  ric (', ni, &
        ' primitives) / ', &
        trim(merge('external Hessian   ', 'Lindh diagonal Hess', ext_hess_ok)), &
        '  maxcyc=', spec%opt_maxcyc, &
        '  conv=', spec%opt_conv_level, '  trust0=', trust
else
   allocate(Hx(n,n))
   if (ext_hess_ok) then
      Hx = Hext
   else
      call lindh_cart_hessian(nat, spec%atomchg, g3, Hx)
   endif
   write(out_unit,'(A,A,A,I0,A,I0,A,F5.3)') '[OPT]  cart / ', &
        trim(merge('external Hessian   ', 'Lindh model Hessian', ext_hess_ok)), &
        '  maxcyc=', spec%opt_maxcyc, '  conv=', spec%opt_conv_level, '  trust0=', trust
endif

have_prev     = .false.
ric_have_prev = .false.
hess_fresh    = .false.
scf_lvl_prev  = -1
lvl_switch    = .false.
converged     = .false.
n_retry       = 0
E_old = 0.0d0 ; pred_prev = 0.0d0 ; ratio = 0.0d0 ; sn = 0.0d0
step_prev = 0.0d0 ; g3_prev = 0.0d0 ; grad_old = 0.0d0

do cyc = 1, spec%opt_maxcyc

   do i = 1, nat
      coord_a(i,1) = g3(1,i) * BOHR2ANG
      coord_a(i,2) = g3(2,i) * BOHR2ANG
      coord_a(i,3) = g3(3,i) * BOHR2ANG
   enddo

   scf_lvl = spec%scf_conv_level
   if (scf_adapt .and. (cyc .eq. 1 .or. gmax .gt. SCF_ADAPT_GMAX)) scf_lvl = 0
   lvl_switch = (cyc .gt. 1 .and. scf_lvl .ne. scf_lvl_prev)
   scf_lvl_prev = scf_lvl
   cosx_force_loose = COSX_FORCE_LOOSE_EARLY .and. (scf_lvl .eq. 0) .and. scf_adapt
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
                 mem_grid_gb, mem_2e_gb, scf_lvl, trim(spec%basedir), &
                 dens_in_a, dens_in_b, dens_out_a, dens_out_b)

   if (iconv /= 1) then
      write(out_unit,'(A,I0,A)') '[OPTFAILED]  SCF did not converge at cycle ', cyc, &
           ' - aborting optimization'
      call finish(spec, g3, nat, out_unit, .false., cyc, energy_out, force_out, iconv)
      return
   endif

   E    = energy_out
   grad = reshape(transpose(force_out), [n])
   if (spec%npc == 0) then
      block
         integer :: ia, ix, k, j, nv
         real(8) :: cen(3), v(n,6), w(n), dot
         cen = sum(g3, dim=2) / dble(nat)
         v = 0.0d0
         do ia = 1,nat
            do ix = 1,3
               v(3*(ia-1)+ix, ix) = 1.0d0
            enddo
            v(3*(ia-1)+2, 4) = -(g3(3,ia)-cen(3)) ; v(3*(ia-1)+3, 4) =  (g3(2,ia)-cen(2))
            v(3*(ia-1)+3, 5) = -(g3(1,ia)-cen(1)) ; v(3*(ia-1)+1, 5) =  (g3(3,ia)-cen(3))
            v(3*(ia-1)+1, 6) = -(g3(2,ia)-cen(2)) ; v(3*(ia-1)+2, 6) =  (g3(1,ia)-cen(1))
         enddo
         nv = 0
         do k = 1,6
            w = v(:,k)
            do j = 1,nv
               w = w - dot_product(v(:,j), w) * v(:,j)
            enddo
            dot = sqrt(sum(w*w))
            if (dot > 1.0d-6 * max(1.0d0, sqrt(sum(v(:,k)**2)))) then
               nv = nv + 1
               v(:,nv) = w / dot
            endif
         enddo
         do j = 1,nv
            grad = grad - dot_product(v(:,j), grad) * v(:,j)
         enddo
      end block
   endif
   if (ric) then
      gproj = grad
      call cons_project_grad(nat, g3, gproj)
      gmax = maxval(abs(gproj))
      grms = sqrt(sum(gproj*gproj) / n)
   else
      call cons_project_grad(nat, g3, grad)
      gmax = maxval(abs(grad))
      grms = sqrt(sum(grad*grad) / n)
   endif

   if (have_prev) then
      ared = E - E_old
      if ((.not. ts_on .and. pred_prev < -1.0d-14) .or. (ts_on .and. abs(pred_prev) > 1.0d-10)) then
         ratio = ared / pred_prev
      else
         ratio = -1.0d0
      endif
      if (ts_on) then
         if (ratio < 0.25d0 .or. ratio > 1.75d0) then
            trust = max(0.5d0 * sn, trust_min)
         else if (ratio > 0.75d0 .and. ratio < 1.25d0 .and. sn > 0.99d0*trust) then
            trust = min(2.0d0 * trust, trust_max)
         endif
      else if (ratio < 0.25d0) then
         trust = max(TRUST_SHRINK * sn, trust_min)
      else if (ratio > 0.75d0 .and. sn > 0.99d0*trust) then
         trust = min(2.0d0 * trust, trust_max)
      endif
   endif

   call write_cycle(out_unit, cyc, E, have_prev, E_old, gmax, grms, have_prev, smax, srms, trust, iconv, &
                    have_prev, ratio, sn)
   call write_cycle_diag(out_unit, nat, spec%atomchg, grad)
   call traj_append(trajfile, nat, spec%atomchg, g3, cyc, E, gmax)
   if (spec%chk_write) call chk_write_safe(spec, dens_out_a, dens_out_b)

   if (have_prev) then
      dE = E - E_old
      converged = (gmax < t_gmax) .and. (grms < t_grms) .and. &
                  ( ((smax < t_smax) .and. (srms < t_srms)) .or. (abs(dE) < t_de) )
   else
      converged = (gmax < t_gmax) .and. (grms < t_grms)
   endif
   if (converged .and. cons_active()) converged = (cons_max_dev(nat, g3) < CONS_DEV_TOL)
   if (converged) then
      write(out_unit,'(A)')       '[OPTCONVERGED]'
      write(out_unit,'(A,I0)')    '  cycles = ', cyc
      write(out_unit,'(A,F20.9)') '  final E (Hartree) = ', E
      if (ts_on .and. spec%npc == 0) then
         block
            real(8), allocatable :: Hchk(:,:), wchk(:), wkchk(:)
            integer :: ichk, nneg_chk
            logical :: okchk
            allocate(Hchk(n,n), wchk(n), wkchk(66*n))
            call hess_get('calc', Hchk, okchk)
            engine_quiet = .true.
            if (okchk) then
               Hchk = 0.5d0 * (Hchk + transpose(Hchk))
               if (.not. cons_has_frozen()) then
                  block
                     real(8) :: trchk(n,6)
                     integer :: ntrchk, mchk
                     call build_tr_basis(nat, g3, trchk, ntrchk)
                     do mchk = 1, ntrchk
                        call dger(n, n, PROJ_SHIFT, trchk(:,mchk), 1, trchk(:,mchk), 1, Hchk, n)
                     enddo
                  end block
               endif
               call dsyev('N', 'U', n, Hchk, n, wchk, wkchk, 66*n, ichk)
               if (ichk == 0) then
                  nneg_chk = count(wchk < -2.0d-5)
                  if (nneg_chk == 1) then
                     write(out_unit,'(A)') '[OPT] Hessian check: one imaginary mode - a transition state'
                  else
                     write(out_unit,'(A,I0,A)') '[OPT] WARNING: the analytic Hessian has ', nneg_chk, &
                          ' imaginary mode(s), a transition state has exactly one - this is not a transition state'
                  endif
               endif
            endif
         end block
      endif
      call finish(spec, g3, nat, out_unit, .true., cyc, E, force_out, iconv)
      call profile_tail(out_unit)
      if (spec%freq) then
         if (spec%npc > 0) error stop 'freq=.true. with external point charges is not supported'
         write(out_unit,'(A)') ''
         write(out_unit,'(A)') '[FREQ] analytic Hessian at the optimized geometry'
         freq_unit = out_unit
         freq_mult = spec%imult
         freq_request = .true.
         call move_alloc(dens_out_a, dens_in_a)
         call move_alloc(dens_out_b, dens_in_b)
         call EngineUp(spec%ncenters, spec%imult, spec%icharge, spec%functional, &
                       coord_a, atomchg_a, spec%baselabel, spec%ecplabel, &
                       spec%atom_basis, spec%atom_ecp, &
                       trim(spec%J), trim(spec%K), trim(spec%ri_aux_basis), spec%spherical, spec%harris_guess, &
                       .false., spec%vv10_nonself, spec%mem_cap_gb, &
                       .false., spec%n_threads, &
                       spec%npc, pc_charge_a, pc_coord_a, &
                       spec%molden_write, trim(spec%molden_file), spec%molden_read, trim(spec%molden_read_file), &
                       spec%cosmo_on, spec%cosmo_epsilon, spec%cosmo_radii_scale, spec%cosmo_avg_area, spec%cosmo_sigma_rav, &
                       trim(spec%cosmo_cavity_type), spec%cosmo_rsolv, spec%cosmo_ks_nseg, spec%cosmo_ks_nface, &
                       trim(spec%cosmo_sigma_profile_file), &
                       spec%cosmo_smd, trim(spec%cosmo_solvent), &
                       force_out, energy_out, MLcharge_out, iconv, econv, &
                       mem_grid_gb, mem_2e_gb, spec%scf_conv_level, trim(spec%basedir), &
                       dens_in_a, dens_in_b, dens_out_a, dens_out_b)
      endif
      return
   endif

   if (have_prev .and. .not. ts_on .and. n_retry < MAXRETRY &
    .and. ((E - E_old) > E_REJECT .or. &
         (E_REJECT_BY_RATIO .and. (E - E_old) > E_REJECT_NOISE .and. ratio < E_REJECT_RATIO)) &
    .and. max(0.25d0 * sn, trust_min) < 0.999d0 * sn) then
      n_retry = n_retry + 1
      write(out_unit,'(A,I0,A,ES9.2,A)') '  [reject #', n_retry, ']  dE=+', E-E_old, &
           ' - backtracking with smaller trust'
      g3    = g3_prev
      grad  = grad_old
      E     = E_old
      trust = max(0.25d0 * sn, trust_min)
      if (ric) then
         call ic_bmatrix(ic, g3, Bmat)
         call ic_ginv(Bmat, ni, n, Ginvm)
         call dgemv('N', ni, n, 1.0d0, Bmat, ni, grad, 1, 0.0d0, bgv, 1)
         call dgemv('N', ni, ni, 1.0d0, Ginvm, ni, bgv, 1, 0.0d0, gq_cur, 1)
         call ric_step(ic, nat, ni, g3, Bmat, Ginvm, gq_cur, Hq, trust, delta, pred, sn, out_unit)
      else
         call rfo_step(n, nat, Hx, grad, g3, trust, delta, pred, out_unit)
         sn = sqrt(sum(delta*delta))
      endif
      pred_prev = pred
      g3 = g3_prev + reshape(delta, [3,nat])
      if (cons_active()) then
         call cons_correct(nat, g3, out_unit)
         delta = reshape(g3 - g3_prev, [n])
      endif
      step_prev = delta
      cycle
   endif
   n_retry = 0

   want_hess_now = hess_every > 0 .and. cyc > 1 .and. &
                   cyc < spec%opt_maxcyc .and. mod(cyc-1, hess_every) == 0
   if (ts_endgame_refresh_due(cyc, gmax, ts_nrecalc, ts_multi_done)) then
      want_hess_now = .true.
      ts_multi_done = .true.
      ts_nrecalc = ts_nrecalc + 1
   endif
   hess_fresh = .false.

   if (ric) then
      if (ts_on) then
         allocate(Hcarry(n,n), Bold(ni,n), Tcarry(ni,n))
         call ic_bmatrix(ic, g3, Bold)
         call dgemm('N','N', ni, n, ni, 1.0d0, Hq, ni, Bold, ni, 0.0d0, Tcarry, ni)
         call dgemm('T','N', n, n, ni, 1.0d0, Bold, ni, Tcarry, ni, 0.0d0, Hcarry, n)
         deallocate(Bold, Tcarry)
         block
            real(8), allocatable :: Kc(:,:)
            allocate(Kc(n,n))
            call ric_grad_curv(ic, g3, nat, ni, grad, Kc)
            Hcarry = Hcarry + Kc
         end block
      endif
      call ic_rebuild_if_changed(ic, nat, spec%atomchg, g3, topo_changed)
      if (topo_changed) then
         ni = ic_nprim(ic)
         deallocate(Hq, q_prev, gq_prev, q_cur, sq, gq_cur, bgv)
         allocate(Hq(ni,ni), q_prev(ni), gq_prev(ni), q_cur(ni), sq(ni), gq_cur(ni), bgv(ni))
         if (ts_on) then
            call cart_hessian_to_internal(ic, g3, nat, ni, Hcarry, Hq, grad)
         else
            call ric_model_hessian(ic, nat, spec%atomchg, g3, ni, Hq, out_unit)
         endif
         ric_have_prev = .false.
         ged_nh = 0
         call ts_forget_mode()
         write(out_unit,'(A,I0,A,A)') '  [ric] coordinate set rebuilt - ', ni, &
              ' primitives, Hessian ', trim(merge('carried over', 'reset       ', ts_on))
         if (spec%opt_write_ric) call dump_ric(ic, g3, outfile, 'rebuilt')
      endif
      if (allocated(Bmat)) then
         if (size(Bmat,1) /= ni) deallocate(Bmat, Ginvm)
      endif
      if (.not. allocated(Bmat)) allocate(Bmat(ni,n), Ginvm(ni,ni))
      call ic_bmatrix(ic, g3, Bmat)
      call ic_ginv(Bmat, ni, n, Ginvm, nrank)
      call dgemv('N', ni, n, 1.0d0, Bmat, ni, grad, 1, 0.0d0, bgv, 1)
      call dgemv('N', ni, ni, 1.0d0, Ginvm, ni, bgv, 1, 0.0d0, gq_cur, 1)
      if (nrank < nexp) then
         write(out_unit,'(A,I0,A,I0,A)') &
            '  [ric] coordinate set became incomplete (rank ', nrank, ' of ', nexp, &
            ') - switching to cart for the rest of the run'
         ric = .false.
         allocate(Hx(n,n))
         if (hq_pending) then
            Hx = Hext
            if (.not. ts_on) call regularize_hessian(Hx, n, out_unit)
            hq_pending = .false.
         else if (ts_on) then
            Hx = Hcarry
         else
            call lindh_cart_hessian(nat, spec%atomchg, g3, Hx)
         endif
         have_prev = .false.
         call ts_forget_mode()
      endif
      if (allocated(Hcarry)) deallocate(Hcarry)
   endif
   if (ric) then
      if (hq_pending) then
         call cart_hessian_to_internal(ic, g3, nat, ni, Hext, Hq, grad)
         if (.not. ts_on) call ric_floor_hessian(ic, g3, nat, ni, Hq, out_unit)
         hq_pending = .false.
      endif
      if (want_hess_now) then
         allocate(Hnew(n,n))
         call hess_get(hsrc_re, Hnew, hok)
         if (hok) then
            call cart_hessian_to_internal(ic, g3, nat, ni, Hnew, Hq, grad)
            if (.not. ts_on) call ric_floor_hessian(ic, g3, nat, ni, Hq, out_unit)
            hess_fresh = .true.
            write(out_unit,'(A,I0)') '  [opt] Hessian ('//trim(hsrc_re)//') recomputed at cycle ', cyc
         endif
         deallocate(Hnew)
         if (allocated(dens_out_a)) then
            call move_alloc(dens_out_a, dens_in_a)
            call move_alloc(dens_out_b, dens_in_b)
         endif
      endif
      if (ric_have_prev .and. .not. hess_fresh .and. .not. lvl_switch) then
         call ic_values(ic, g3, q_cur)
         call ic_delta_q(ic, q_cur, q_prev, sq)
         call hess_update(ni, Hq, sq, gq_cur - gq_prev)
      endif
      block
         real(8), allocatable :: Hs(:,:)
         allocate(Hs(n,n))
         call ric_hessian_to_cart(ic, g3, nat, ni, Bmat, Hq, grad, Hs)
         call write_hessian_file(hessfile, n, Hs)
      end block
      if (hess_diag_on) call report_hq_spectrum(ni, Bmat, Ginvm, n, Hq, out_unit)
      ged_used = .false.
      if (ged_m >= 2 .and. .not. ts_on .and. .not. cons_active()) then
         if (hess_fresh .or. lvl_switch) ged_nh = 0
         if (.not. allocated(ged_q)) allocate(ged_q(ni,ged_m+1), ged_g(ni,ged_m+1), ged_E(ged_m+1))
         if (size(ged_q,1) /= ni) then
            deallocate(ged_q, ged_g, ged_E)
            allocate(ged_q(ni,ged_m+1), ged_g(ni,ged_m+1), ged_E(ged_m+1))
            ged_nh = 0
         endif
         if (ged_nh == ged_m) then
            ged_q(:,1:ged_m-1) = ged_q(:,2:ged_m)
            ged_g(:,1:ged_m-1) = ged_g(:,2:ged_m)
            ged_E(1:ged_m-1)   = ged_E(2:ged_m)
            ged_nh = ged_m - 1
         endif
         ged_nh = ged_nh + 1
         call ic_values(ic, g3, ged_q(:,ged_nh))
         ged_g(:,ged_nh) = gq_cur
         ged_E(ged_nh)   = E
         if (ged_nh >= 2) then
            allocate(ged_dQ(ni,ged_nh), ged_c(ged_nh), ged_qoff(ni), ged_gstar(ni))
            do ged_i = 1, ged_nh
               call ic_delta_q(ic, ged_q(:,ged_i), ged_q(:,ged_nh), ged_dQ(:,ged_i))
            enddo
            call gediis_solve(ged_nh, ni, ged_E(1:ged_nh), ged_dQ, ged_g(:,1:ged_nh), 3.0d0, &
                              ged_c, ged_qoff, ged_gstar, ged_Estar, ged_ok)
            if (ged_ok) ged_ok = (ged_Estar < E - 1.0d-10 .and. sqrt(sum(ged_qoff*ged_qoff)) <= trust)
            if (ged_ok) then
               call ric_step(ic, nat, ni, g3, Bmat, Ginvm, ged_gstar, Hq, trust, delta, pred, sn, out_unit, &
                             dq_pre=ged_qoff, e_off=ged_Estar - E)
               ged_used = .true.
               write(out_unit,'(A,I0,A,ES10.3,A,F7.3)') '  [gediis] ', ged_nh, ' points, predicted E(q*)-E = ', &
                    ged_Estar - E, ', max |c| = ', maxval(abs(ged_c))
            endif
            deallocate(ged_dQ, ged_c, ged_qoff, ged_gstar)
         endif
      endif
      if (.not. ged_used) call ric_step(ic, nat, ni, g3, Bmat, Ginvm, gq_cur, Hq, trust, delta, pred, sn, out_unit)
      if (ts_lost_mode_recalc_due(hess_fresh, ts_nrecalc)) then
         allocate(Hnew(n,n))
         call hess_get(hsrc_re, Hnew, hok)
         if (hok) then
            ts_nrecalc = ts_nrecalc + 1
            call cart_hessian_to_internal(ic, g3, nat, ni, Hnew, Hq, grad)
            if (.not. ts_on) call ric_floor_hessian(ic, g3, nat, ni, Hq, out_unit)
            hess_fresh = .true.
            call ts_forget_mode()
            write(out_unit,'(A,I0)') '  [ts] no negative curvature left - Hessian ('//trim(hsrc_re)//') recomputed at cycle ', cyc
            call ric_step(ic, nat, ni, g3, Bmat, Ginvm, gq_cur, Hq, trust, delta, pred, sn, out_unit)
         endif
         deallocate(Hnew)
         if (allocated(dens_out_a)) then
            call move_alloc(dens_out_a, dens_in_a)
            call move_alloc(dens_out_b, dens_in_b)
         endif
      endif
      call ic_values(ic, g3, q_prev)
      gq_prev = gq_cur
      ric_have_prev = .true.
   else
      if (want_hess_now) then
         allocate(Hnew(n,n))
         call hess_get(hsrc_re, Hnew, hok)
         if (hok) then
            Hx = Hnew
            hess_fresh = .true.
            write(out_unit,'(A,I0)') '  [opt] Hessian ('//trim(hsrc_re)//') recomputed at cycle ', cyc
         endif
         deallocate(Hnew)
         if (allocated(dens_out_a)) then
            call move_alloc(dens_out_a, dens_in_a)
            call move_alloc(dens_out_b, dens_in_b)
         endif
      endif
      if (have_prev .and. .not. hess_fresh .and. .not. lvl_switch) then
         call hess_update(n, Hx, step_prev, grad - grad_old)
      endif
      call write_hessian_file(hessfile, n, Hx)
      call rfo_step(n, nat, Hx, grad, g3, trust, delta, pred, out_unit)
      if (ts_lost_mode_recalc_due(hess_fresh, ts_nrecalc)) then
         allocate(Hnew(n,n))
         call hess_get(hsrc_re, Hnew, hok)
         if (hok) then
            ts_nrecalc = ts_nrecalc + 1
            Hx = Hnew
            hess_fresh = .true.
            call ts_forget_mode()
            write(out_unit,'(A,I0)') '  [ts] no negative curvature left - Hessian ('//trim(hsrc_re)//') recomputed at cycle ', cyc
            call rfo_step(n, nat, Hx, grad, g3, trust, delta, pred, out_unit)
         endif
         deallocate(Hnew)
         if (allocated(dens_out_a)) then
            call move_alloc(dens_out_a, dens_in_a)
            call move_alloc(dens_out_b, dens_in_b)
         endif
      endif
      sn = sqrt(sum(delta*delta))
   endif

   g3_prev   = g3
   E_old     = E
   grad_old  = grad
   pred_prev = pred
   g3 = g3 + reshape(delta, [3,nat])
   if (cons_active()) then
      call cons_correct(nat, g3, out_unit)
      delta = reshape(g3 - g3_prev, [n])
   endif
   step_prev = delta
   smax = maxval(abs(delta))
   srms = sqrt(sum(delta*delta) / n)
   have_prev = .true.

   if (allocated(dens_out_a)) then
      call move_alloc(dens_out_a, dens_in_a)
      call move_alloc(dens_out_b, dens_in_b)
   endif
enddo

write(out_unit,'(A,I0,A)') '[OPTFAILED]  not converged in ', spec%opt_maxcyc, ' cycles'
call finish(spec, g3, nat, out_unit, .false., spec%opt_maxcyc, E, force_out, iconv)

contains

pure function lower8(x) result(y)
character(len=*), intent(in) :: x
character(len=8) :: y
integer :: i, c
y = adjustl(x)
do i = 1, len(y)
   c = iachar(y(i:i))
   if (c >= 65 .and. c <= 90) y(i:i) = achar(c+32)
enddo
end function lower8

subroutine hess_get(src, Hout, ok)
character(len=*), intent(in) :: src
real(8), intent(out) :: Hout(n,n)
logical, intent(out) :: ok
n_hess_calls = n_hess_calls + 1
if (trim(src) == 'xtb') then
   call hessian_from_xtb(spec, g3, nat, outfile, out_unit, (ts_on .or. ric), Hout, ok)
else
   call hessian_analytic(spec, g3, nat, out_unit, max(spec%scf_conv_level, 1), &
                         dens_in_a, dens_in_b, dens_out_a, dens_out_b, Hout, ok)
endif
end subroutine hess_get

end subroutine geomopt_run

subroutine write_cycle_diag(u, nat, Z, grad)
integer, intent(in) :: u, nat, Z(nat)
real(8), intent(in) :: grad(3*nat)
integer :: i, k, ord(nat), tmp
real(8) :: gat(nat)
character(len=2) :: nm
integer, parameter :: NTOP = 4
do i = 1, nat
   gat(i) = sqrt(sum(grad(3*i-2:3*i)**2))
   ord(i) = i
enddo
do k = 1, min(NTOP, nat)
   do i = k+1, nat
      if (gat(ord(i)) > gat(ord(k))) then
         tmp = ord(k) ; ord(k) = ord(i) ; ord(i) = tmp
      endif
   enddo
enddo
write(u,'(A,I0,A,ES10.3,A,ES10.3)') '  scf iters = ', scf_hist_niter(), &
     '   scf dP = ', scf_hist_last_prms(), '   scf dE = ', scf_hist_last_de()
write(u,'(A)',advance='no') '  largest |g| on: '
do k = 1, min(NTOP, nat)
   nm = elem_to_name(Z(ord(k)))
   write(u,'(A,A,I0,A,ES9.2)',advance='no') ' ', trim(nm), ord(k), '=', gat(ord(k))
enddo
write(u,'(A)') ''
end subroutine write_cycle_diag

subroutine report_hq_spectrum(ni, B, Ginv, ndof, Hq, out_unit)
integer, intent(in) :: ni, ndof, out_unit
real(8), intent(in) :: B(ni,ndof), Ginv(ni,ni), Hq(ni,ni)
real(8), allocatable :: G(:,:), P(:,:), T(:,:), Hp(:,:), w(:), work(:)
integer :: lwork, info, nnz
allocate(G(ni,ni), P(ni,ni), T(ni,ni), Hp(ni,ni), w(ni))
call dgemm('N','T', ni, ni, ndof, 1.0d0, B, ni, B, ni, 0.0d0, G, ni)
call dgemm('N','N', ni, ni, ni, 1.0d0, Ginv, ni, G, ni, 0.0d0, P, ni)
call dgemm('N','N', ni, ni, ni, 1.0d0, Hq, ni, P, ni, 0.0d0, T, ni)
call dgemm('N','N', ni, ni, ni, 1.0d0, P, ni, T, ni, 0.0d0, Hp, ni)
lwork = max(1, 3*ni + 2*ni*ni)
allocate(work(lwork))
call dsyev('N', 'U', ni, Hp, ni, w, work, lwork, info)
if (info /= 0) return
nnz = count(abs(w) > 1.0d-8)
write(out_unit,'(A,I0,A,ES10.3,A,ES10.3,A,I0)') '  [hess] P.Hq.P nonzero = ', nnz, &
     '   min = ', minval(w, mask = abs(w) > 1.0d-8), &
     '   max = ', maxval(w), '   negative = ', count(w < -1.0d-8)
end subroutine report_hq_spectrum

subroutine dump_ric(ic, coord3, outfile, tag)
type(ic_set_t), intent(in)   :: ic
real(8), intent(in)          :: coord3(3,ic%nat)
character(len=*), intent(in) :: outfile, tag
character(len=256) :: rf
integer :: u, l
l = len_trim(outfile)
if (l > 4 .and. outfile(max(l-3,1):l) == '.out') then
   rf = outfile(1:l-4) // '.ric'
else
   rf = trim(outfile) // '.ric'
endif
open(newunit=u, file=trim(rf), status='unknown', position='append', action='write')
call ic_write(ic, coord3, u, tag)
write(u,'(A)') ''
close(u)
end subroutine dump_ric

logical function is_linear(nat, coord)
integer, intent(in) :: nat
real(8), intent(in) :: coord(3,nat)
integer :: i
real(8) :: c(3), v(3), ax(3), nv, cross(3)
is_linear = .false.
if (nat < 3) then
   is_linear = .true. ; return
endif
c = 0.0d0
do i = 1, nat
   c = c + coord(:,i)
enddo
c = c / real(nat,8)
ax = 0.0d0
do i = 1, nat
   v = coord(:,i) - c
   nv = sqrt(sum(v*v))
   if (nv > 1.0d-4) then
      ax = v / nv ; exit
   endif
enddo
if (sum(ax*ax) < 0.5d0) return
do i = 1, nat
   v = coord(:,i) - c
   cross(1) = ax(2)*v(3) - ax(3)*v(2)
   cross(2) = ax(3)*v(1) - ax(1)*v(3)
   cross(3) = ax(1)*v(2) - ax(2)*v(1)
   if (sqrt(sum(cross*cross)) > 1.0d-4) return
enddo
is_linear = .true.
end function is_linear

subroutine rfo_step(n, nat, Hx, grad, coord3, trust, delta, pred, out_unit)
integer, intent(in)  :: n, nat, out_unit
real(8), intent(in)  :: Hx(n,n), grad(n), coord3(3,nat), trust
real(8), intent(out) :: delta(n), pred
real(8), allocatable :: P(:,:), T1(:,:), Hp(:,:), gp(:), tr(:,:), Hd(:), Vc(:,:), Wc(:,:)
integer :: ntr, i, j, k, nvc, nw
real(8) :: sn

allocate(P(n,n), T1(n,n), Hp(n,n), gp(n), Hd(n))
call cons_basis(nat, coord3, .false., Vc, nvc)
allocate(tr(n,6), Wc(n,6+max(nvc,1)))
call build_tr_basis(nat, coord3, tr, ntr)
if (cons_active()) then
   if (cons_has_frozen()) ntr = 0
   Wc(:,1:ntr) = tr(:,1:ntr)
   if (nvc > 0) Wc(:,ntr+1:ntr+nvc) = Vc(:,1:nvc)
   deallocate(tr)
   allocate(tr(n,ntr+nvc))
   call orth_cols(Wc(:,1:ntr+nvc), n, ntr+nvc, tr, nw)
   ntr = nw
endif

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

call solve_step(n, Hp, gp, delta, out_unit)

sn = sqrt(sum(delta*delta))
if (sn > trust .and. sn > 0.0d0) delta = delta * (trust / sn)

call dgemv('N', n, n, 1.0d0, Hp, n, delta, 1, 0.0d0, Hd, 1)
pred = dot_product(gp, delta) + 0.5d0 * dot_product(delta, Hd)
end subroutine rfo_step

subroutine ric_step(ic, nat, ni, coord3, B, Ginv, gq, Hq, trust, dx, pred, qnorm, out_unit, dq_pre, e_off)
type(ic_set_t), intent(in) :: ic
integer, intent(in)  :: nat, ni, out_unit
real(8), intent(in)  :: coord3(3,nat), B(ni,3*nat), Ginv(ni,ni), gq(ni)
real(8), intent(in)  :: Hq(ni,ni), trust
real(8), intent(out) :: dx(3*nat), pred
real(8), intent(out) :: qnorm
real(8), intent(in), optional :: dq_pre(ni), e_off
real(8), allocatable :: G(:,:), P(:,:), T1(:,:), Hp(:,:)
real(8), allocatable :: gp(:), dq(:), dqp(:), Hd(:)
real(8), allocatable :: Vc(:,:), BV(:,:), U(:,:), Uo(:,:), Qm(:,:), dqpre_p(:)
integer :: ndof, i, j, nvc, nu
real(8) :: sn, dxn
logical :: ok

ndof = 3*nat
allocate(G(ni,ni), P(ni,ni), T1(ni,ni), Hp(ni,ni))
allocate(gp(ni), dq(ni), dqp(ni), Hd(ni))

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

call cons_basis(nat, coord3, .true., Vc, nvc)
if (nvc > 0) then
   allocate(BV(ni,nvc), U(ni,nvc), Uo(ni,nvc))
   call dgemm('N','N', ni, nvc, ndof, 1.0d0, B, ni, Vc, ndof, 0.0d0, BV, ni)
   call dgemm('N','N', ni, nvc, ni, 1.0d0, Ginv, ni, BV, ni, 0.0d0, U, ni)
   call orth_cols(U, ni, nvc, Uo, nu)
   if (nu > 0) then
      allocate(Qm(ni,ni))
      Qm = 0.0d0
      do i = 1, ni
         Qm(i,i) = 1.0d0
      enddo
      call dgemm('N','T', ni, ni, nu, -1.0d0, Uo, ni, Uo, ni, 1.0d0, Qm, ni)
      call dgemm('N','N', ni, ni, ni, 1.0d0, Hp, ni, Qm, ni, 0.0d0, T1, ni)
      call dgemm('N','N', ni, ni, ni, 1.0d0, Qm, ni, T1, ni, 0.0d0, Hp, ni)
      call dgemm('N','T', ni, ni, nu, PROJ_SHIFT, Uo, ni, Uo, ni, 1.0d0, Hp, ni)
      call dgemv('N', ni, ni, 1.0d0, Qm, ni, gp, 1, 0.0d0, T1(:,1), 1)
      gp = T1(:,1)
   endif
endif

call solve_step(ni, Hp, gp, dq, out_unit)
call dgemv('N', ni, ni, 1.0d0, P, ni, dq, 1, 0.0d0, dqp, 1)
dq = dqp
if (allocated(Qm)) then
   call dgemv('N', ni, ni, 1.0d0, Qm, ni, dq, 1, 0.0d0, dqp, 1)
   dq = dqp
endif

sn = sqrt(sum(dq*dq))
if (sn > trust .and. sn > 0.0d0) dq = dq * (trust / sn)

if (present(dq_pre)) then
   allocate(dqpre_p(ni))
   call dgemv('N', ni, ni, 1.0d0, P, ni, dq_pre, 1, 0.0d0, dqpre_p, 1)
   dq = dq + dqpre_p
endif
qnorm = sqrt(sum(dq*dq))

call ic_back_transform(ic, coord3, dq, dx, ok)
if (.not. ok) write(out_unit,'(A)') &
     '  [ric] back-transform not fully converged - using linear estimate'

dxn = 0.0d0
do i = 1, nat
   dxn = max(dxn, sqrt(dx(3*i-2)**2 + dx(3*i-1)**2 + dx(3*i)**2))
enddo
if (dxn > 2.0d0*trust .and. dxn > 0.0d0) then
   dq = dq * (2.0d0*trust / dxn)
   dx = dx * (2.0d0*trust / dxn)
   if (allocated(dqpre_p)) dqpre_p = dqpre_p * (2.0d0*trust / dxn)
endif

if (present(dq_pre)) then
   dq = dq - dqpre_p
   call dgemv('N', ni, ni, 1.0d0, Hp, ni, dq, 1, 0.0d0, Hd, 1)
   pred = e_off + dot_product(gp, dq) + 0.5d0 * dot_product(dq, Hd)
else
   call dgemv('N', ni, ni, 1.0d0, Hp, ni, dq, 1, 0.0d0, Hd, 1)
   pred = dot_product(gp, dq) + 0.5d0 * dot_product(dq, Hd)
endif
end subroutine ric_step

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

subroutine solve_step(n, H, g, step, out_unit)
integer, intent(in)  :: n, out_unit
real(8), intent(in)  :: H(n,n), g(n)
real(8), intent(out) :: step(n)
if (ts_on) then
   call solve_prfo(n, H, g, step, out_unit)
else
   call solve_rfo(n, H, g, step, out_unit)
endif
end subroutine solve_step

subroutine hess_update(n, H, s, y)
integer, intent(in)    :: n
real(8), intent(inout) :: H(n,n)
real(8), intent(in)    :: s(n), y(n)
if (ts_on) then
   call bofill_update(n, H, s, y)
else
   call bfgs_update(n, H, s, y)
endif
end subroutine hess_update

subroutine write_cycle(u, cyc, E, have_dE, E_old, gmax, grms, have_step, smax, srms, trust, iconv, &
                       have_ratio, ratio, sn)
integer, intent(in) :: u, cyc, iconv
real(8), intent(in) :: E, E_old, gmax, grms, smax, srms, trust
logical, intent(in) :: have_dE, have_step
logical, intent(in) :: have_ratio
real(8), intent(in) :: ratio, sn
write(u,'(A,I0,A)') '[OPTCYCLE ', cyc, ']'
if (have_dE) then
   write(u,'(A,F18.9,A,ES11.3)') '  E = ', E, '   dE = ', E - E_old
else
   write(u,'(A,F18.9,A)')        '  E = ', E, '   dE =        --'
endif
write(u,'(A,ES11.3,A,ES11.3)') '  gmax = ', gmax, '   grms = ', grms
if (have_step) write(u,'(A,ES11.3,A,ES11.3)') '  smax = ', smax, '   srms = ', srms
if (have_ratio) then
   write(u,'(A,F6.3,A,ES10.3,A,F7.3,A,I0,A)') '  trust = ', trust, &
        '   |step| = ', sn, '   ratio = ', ratio, '   scf = ', iconv, ' (1=converged)'
else
   write(u,'(A,F6.3,A,I0,A)') '  trust = ', trust, '   scf = ', iconv, ' (1=converged)'
endif
end subroutine write_cycle

subroutine finish(spec, g3, nat, u, ok, ncyc, E, force_out, iconv)
type(engine_input_t), intent(inout) :: spec
integer, intent(in) :: nat, u, ncyc, iconv
real(8), intent(in) :: g3(3,nat), E, force_out(nat,3)
logical, intent(in) :: ok
integer :: i
engine_quiet = .false.
cosx_force_loose = .false.
ts_on = .false.
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
use mod_profile, only: prof_report, profile_enabled
integer, intent(in) :: u
if (u == 6) then
   if (profile_enabled()) then
      write(u,'(A)') '[PROFILE] (final optimization cycle)'
      call prof_report(u)
   endif
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
