! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Geometry optimizer for the Direwolf harness.

module mod_geomopt
use mod_engine_input_types, only: engine_input_t
use mod_engineup_interface, only: EngineUp
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

integer :: nat, n, i, cyc, n_retry, ni, nrank, nexp
logical :: ext_hess_ok, hess_diag_on
character(32) :: hd_env
integer :: hd_stat
real(8), allocatable :: Hext(:,:)
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
ric = (trim(spec%opt_coord) == 'ric') .or. (trim(spec%opt_coord) == 'internal') &
      .or. (trim(spec%opt_coord) == 'redundant')

engine_quiet = .true.

force_dense_allow = .false.

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

if (spec%scf_conv_level .lt. 1) then
   spec%scf_conv_level = 1
   write(out_unit,'(A)') &
      '[OPT]  SCF convergence raised to fine for the optimization (gradient quality)'
endif
if (spec%opt_conv_level .ge. 1 .and. spec%scf_conv_level .lt. 2) then
   spec%scf_conv_level = 2
   write(out_unit,'(A)') &
      '[OPT]  opt_conv=tight -> SCF convergence raised to tight'
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
ext_hess_ok = .false.
if (len_trim(spec%opt_hessian_file) > 0) then
   allocate(Hext(n,n))
   call read_ext_hessian(trim(spec%opt_hessian_file), n, Hext, ext_hess_ok, out_unit)
   if (ext_hess_ok) then
      write(out_unit,'(A)') '[OPT]  initial Hessian read from '// &
           trim(spec%opt_hessian_file)
   else
      deallocate(Hext)
   endif
endif

if (ric) then
   allocate(Hq(ni,ni), q_prev(ni), gq_prev(ni), q_cur(ni), sq(ni), gq_cur(ni))
   allocate(bgv(ni))
   if (ext_hess_ok) then
      call cart_hessian_to_internal(ic, g3, nat, ni, Hext, Hq)
   else
      call ric_model_hessian(ic, nat, spec%atomchg, g3, ni, Hq, out_unit)
   endif
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
                 mem_grid_gb, mem_2e_gb, spec%scf_conv_level, trim(spec%basedir), &
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
         call ic_bmatrix(ic, g3, Bmat)
         call ic_ginv(Bmat, ni, n, Ginvm)
         call dgemv('N', ni, n, 1.0d0, Bmat, ni, grad, 1, 0.0d0, bgv, 1)
         call dgemv('N', ni, ni, 1.0d0, Ginvm, ni, bgv, 1, 0.0d0, gq_cur, 1)
         call ric_step(ic, nat, ni, g3, Bmat, Ginvm, gq_cur, Hq, trust, delta, pred, sn, out_unit)
      else
         call rfo_step(n, nat, Hx, grad, g3, trust, delta, pred, out_unit)
         sn = sqrt(sum(delta*delta))
      endif
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
         deallocate(Hq, q_prev, gq_prev, q_cur, sq, gq_cur, bgv)
         allocate(Hq(ni,ni), q_prev(ni), gq_prev(ni), q_cur(ni), sq(ni), gq_cur(ni), bgv(ni))
         call ric_model_hessian(ic, nat, spec%atomchg, g3, ni, Hq, out_unit)
         ric_have_prev = .false.
         write(out_unit,'(A,I0,A)') '  [ric] coordinate set rebuilt - ', ni, &
              ' primitives, Hessian reset'
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
         call lindh_cart_hessian(nat, spec%atomchg, g3, Hx)
         have_prev = .false.
      endif
      if (ric_have_prev) then
         call ic_values(ic, g3, q_cur)
         call ic_delta_q(ic, q_cur, q_prev, sq)
         call bfgs_update(ni, Hq, sq, gq_cur - gq_prev)
      endif
      if (hess_diag_on) call report_hq_spectrum(ni, Bmat, Ginvm, n, Hq, out_unit)
      call ric_step(ic, nat, ni, g3, Bmat, Ginvm, gq_cur, Hq, trust, delta, pred, sn, out_unit)
      call ic_values(ic, g3, q_prev)
      gq_prev = gq_cur
      ric_have_prev = .true.
   else
      if (have_prev) call bfgs_update(n, Hx, step_prev, grad - grad_old)
      call rfo_step(n, nat, Hx, grad, g3, trust, delta, pred, out_unit)
      sn = sqrt(sum(delta*delta))
   endif

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

subroutine read_ext_hessian(path, n, H, ok, out_unit)
character(len=*), intent(in) :: path
integer, intent(in)  :: n, out_unit
real(8), intent(out) :: H(n,n)
logical, intent(out) :: ok
integer :: u, ios, got, i, j
real(8) :: buf(n*n), asym
character(len=512) :: line
ok = .false.
open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
if (ios /= 0) then
   write(out_unit,'(A)') '[OPT] cannot open Hessian file '//trim(path)
   return
endif
got = 0
do
   read(u,'(A)',iostat=ios) line
   if (ios /= 0) exit
   line = adjustl(line)
   if (len_trim(line) == 0) cycle
   if (line(1:1) == '$' .or. line(1:1) == '#' .or. line(1:1) == '!') cycle
   read(line,*,iostat=ios) (buf(got+i), i = 1, count_reals(line))
   if (ios /= 0) cycle
   got = got + count_reals(line)
   if (got >= n*n) exit
enddo
close(u)
if (got < n*n) then
   write(out_unit,'(A,I0,A,I0,A)') '[OPT] Hessian file has ', got, &
        ' values, expected ', n*n, ' - ignoring it'
   return
endif
do j = 1, n
   do i = 1, n
      H(i,j) = buf((j-1)*n + i)
   enddo
enddo
asym = maxval(abs(H - transpose(H)))
if (asym > 1.0d-6) write(out_unit,'(A,ES10.3,A)') &
     '[OPT] warning: supplied Hessian is asymmetric by ', asym, ' - symmetrizing'
H = 0.5d0 * (H + transpose(H))
call regularize_hessian(H, n, out_unit)
ok = .true.
end subroutine read_ext_hessian

subroutine regularize_hessian(H, n, out_unit)
integer, intent(in)    :: n, out_unit
real(8), intent(inout) :: H(n,n)
real(8), allocatable :: V(:,:), w(:), work(:), T(:,:)
integer :: i, j, lwork, info, nneg, nlift
real(8), parameter :: HMIN = 5.0d-3
allocate(V(n,n), w(n))
V = H
lwork = max(1, 3*n + 2*n*n)
allocate(work(lwork))
call dsyev('V', 'U', n, V, n, w, work, lwork, info)
if (info /= 0) then
   write(out_unit,'(A)') '[OPT] could not diagonalize the supplied Hessian - using it as given'
   return
endif
nneg  = count(w < -1.0d-6)
nlift = count(w < HMIN)
do i = 1, n
   w(i) = max(w(i), HMIN)
enddo
allocate(T(n,n))
do j = 1, n
   do i = 1, n
      T(i,j) = V(i,j) * w(j)
   enddo
enddo
call dgemm('N','T', n, n, n, 1.0d0, T, n, V, n, 0.0d0, H, n)
if (nneg > 0) then
   write(out_unit,'(A,I0,A,I0,A,ES9.2)') &
        '[OPT]  supplied Hessian had ', nneg, ' negative eigenvalue(s); ', nlift, &
        ' eigenvalue(s) raised to ', HMIN
else
   write(out_unit,'(A,I0,A,ES9.2)') '[OPT]  supplied Hessian: ', nlift, &
        ' eigenvalue(s) raised to ', HMIN
endif
end subroutine regularize_hessian

integer function count_reals(line)
character(len=*), intent(in) :: line
integer :: i
logical :: inword
count_reals = 0
inword = .false.
do i = 1, len_trim(line)
   if (line(i:i) == ' ' .or. line(i:i) == char(9)) then
      inword = .false.
   else if (.not. inword) then
      inword = .true. ; count_reals = count_reals + 1
   endif
enddo
end function count_reals

subroutine cart_hessian_to_internal(ic, coord3, nat, ni, Hx, Hq)
type(ic_set_t), intent(in) :: ic
integer, intent(in)  :: nat, ni
real(8), intent(in)  :: coord3(3,nat), Hx(3*nat,3*nat)
real(8), intent(out) :: Hq(ni,ni)
real(8), allocatable :: B(:,:), Ginv(:,:), A(:,:), T(:,:)
integer :: ndof
ndof = 3*nat
allocate(B(ni,ndof), Ginv(ni,ni), A(ni,ndof), T(ni,ndof))
call ic_bmatrix(ic, coord3, B)
call ic_ginv(B, ni, ndof, Ginv)
call dgemm('N','N', ni, ndof, ni, 1.0d0, Ginv, ni, B, ni, 0.0d0, A, ni)
call dgemm('N','N', ni, ndof, ndof, 1.0d0, A, ni, Hx, ndof, 0.0d0, T, ni)
call dgemm('N','T', ni, ni, ndof, 1.0d0, T, ni, A, ni, 0.0d0, Hq, ni)
end subroutine cart_hessian_to_internal

subroutine ric_model_hessian(ic, nat, Z, coord3, ni, Hq, out_unit)
type(ic_set_t), intent(in) :: ic
integer, intent(in)  :: nat, ni, out_unit
integer, intent(in)  :: Z(nat)
real(8), intent(in)  :: coord3(3,nat)
real(8), intent(out) :: Hq(ni,ni)
real(8), allocatable :: Hx(:,:), V(:,:), w(:), work(:), T(:,:)
integer :: ndof, i, j, lwork, info, nlift
character(len=32) :: hessenv
real(8), parameter :: HQ_MIN = 1.0d-3

hessenv = ''
call get_environment_variable("ENGINE_OPT_HESS", hessenv)
if (trim(hessenv) /= 'full') then
   call ic_lindh_hessian(ic, coord3, Hq)
   return
endif
write(out_unit,'(A)') '  [opt] ENGINE_OPT_HESS=full - full Lindh model Hessian in the primitive basis'

ndof = 3*nat
allocate(Hx(ndof,ndof))
call lindh_cart_hessian(nat, Z, coord3, Hx)
call cart_hessian_to_internal(ic, coord3, nat, ni, Hx, Hq)
deallocate(Hx)

do j = 1, ni
   do i = j+1, ni
      Hq(i,j) = 0.5d0*(Hq(i,j) + Hq(j,i))
      Hq(j,i) = Hq(i,j)
   enddo
enddo

allocate(V(ni,ni), w(ni))
V = Hq
lwork = max(1, 3*ni + 2*ni*ni)
allocate(work(lwork))
call dsyev('V', 'U', ni, V, ni, w, work, lwork, info)
if (info /= 0) then
   write(out_unit,'(A)') '  [opt] model-Hessian diagonalization failed - using it unfloored'
   return
endif
nlift = count(w < HQ_MIN)
do i = 1, ni
   w(i) = max(w(i), HQ_MIN)
enddo
allocate(T(ni,ni))
do j = 1, ni
   do i = 1, ni
      T(i,j) = V(i,j) * w(j)
   enddo
enddo
call dgemm('N','T', ni, ni, ni, 1.0d0, T, ni, V, ni, 0.0d0, Hq, ni)
write(out_unit,'(A,I0,A,I0,A,ES9.2)') '  [opt] full Lindh model Hessian: ', &
     nlift, ' of ', ni, ' eigenvalue(s) raised to ', HQ_MIN
end subroutine ric_model_hessian

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

subroutine ric_step(ic, nat, ni, coord3, B, Ginv, gq, Hq, trust, dx, pred, qnorm, out_unit)
type(ic_set_t), intent(in) :: ic
integer, intent(in)  :: nat, ni, out_unit
real(8), intent(in)  :: coord3(3,nat), B(ni,3*nat), Ginv(ni,ni), gq(ni)
real(8), intent(in)  :: Hq(ni,ni), trust
real(8), intent(out) :: dx(3*nat), pred
real(8), intent(out) :: qnorm
real(8), allocatable :: G(:,:), P(:,:), T1(:,:), Hp(:,:)
real(8), allocatable :: gp(:), dq(:), dqp(:), Hd(:)
integer :: ndof, i, j
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

call solve_rfo(ni, Hp, gp, dq, out_unit)
call dgemv('N', ni, ni, 1.0d0, P, ni, dq, 1, 0.0d0, dqp, 1)
dq = dqp

sn = sqrt(sum(dq*dq))
if (sn > trust .and. sn > 0.0d0) dq = dq * (trust / sn)

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
endif

call dgemv('N', ni, ni, 1.0d0, Hp, ni, dq, 1, 0.0d0, Hd, 1)
pred = dot_product(gp, dq) + 0.5d0 * dot_product(dq, Hd)
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
