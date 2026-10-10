! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Intrinsic reaction coordinate (IRC) for the Direwolf harness.

module mod_irc
use mod_engine_input_types, only: engine_input_t
use mod_engineup_interface, only: EngineUp
use mod_hessian_freq, only: atomic_mass
use mod_integrals, only: mayer_request, mayer_bo
use mod_opt_hessian, only: ANG2BOHR, BOHR2ANG, hessian_analytic, remove_file
use mod_opt_ts, only: bofill_update
use mod_engine_input_elements, only: elem_to_name
use GRID_info, only: force_dense_allow
use MOL_info, only: engine_quiet
implicit none
private
public :: irc_run

real(8), parameter :: HARTREE_TO_KCAL = 627.5094740631d0
real(8), parameter :: AMU_TO_ME = 1822.888486209d0
real(8), parameter :: AU_TO_CM = 219474.6313632d0
real(8), parameter :: E_RISE = 1.0d-6
real(8), parameter :: NEG_EIG = -1.0d-6
integer, parameter :: MAXRETRY = 4
integer, parameter :: BO_WINDOW = 3
real(8), parameter :: BO_GMAX = 2.0d-3
real(8), parameter :: BO_E_MIN = 1.0d0 / 627.5094740631d0
real(8), parameter :: DS_MIN = 0.01d0
real(8), parameter :: BS_ERR_MAX = 1.0d-3
real(8), parameter :: BS_SHIFT_MAX = 0.25d0

contains

subroutine irc_run(spec, outfile, out_unit)
type(engine_input_t), intent(inout) :: spec
character(len=*), intent(in)        :: outfile
integer, intent(in)                 :: out_unit

integer :: nat, n, i, maxp, nf, nb
real(8), allocatable :: x_ts(:,:), g_ts(:), Hy_ts(:,:), v(:), sm(:)
real(8), allocatable :: Xf(:,:), Ef(:), Sf(:), Gf(:), Xb(:,:), Eb(:), Sb(:), Gb(:)
real(8) :: E_ts, ds0, lam, freq_cm
integer :: npair
real(8), allocatable :: bo_last(:), bo_ts(:)
integer :: nneg, imode, ntr
real(8), allocatable :: lamv(:), Uv(:,:)
logical :: ok, do_fwd, do_bwd
character(len=24) :: reason_f, reason_b
character(len=16) :: dirn

integer, allocatable :: atomchg_a(:)
real(8), allocatable :: coord_a(:,:), force_out(:,:), MLcharge_out(:)
real(8), allocatable :: pc_charge_a(:), pc_coord_a(:,:)
real(8), allocatable :: dens_in_a(:,:), dens_in_b(:,:), dens_out_a(:,:), dens_out_b(:,:)
integer :: iconv, n_hess_calls, n_eg_calls
real(8) :: energy_out, econv, mem_grid_gb, mem_2e_gb

nat = spec%ncenters
n   = 3 * nat
ds0 = spec%irc_step
maxp = max(spec%irc_maxpoints, 1)
n_hess_calls = 0
n_eg_calls = 0
npair = nat*(nat-1)/2
allocate(bo_last(npair), bo_ts(npair))
bo_last = 0.0d0
mayer_request = .true.

dirn = lower16(spec%irc_direction)
if (len_trim(dirn) == 0) dirn = 'both'
if (trim(dirn) /= 'both' .and. trim(dirn) /= 'forward' .and. trim(dirn) /= 'backward') then
   write(out_unit,'(A)') '[IRCFAILED]  irc_direction = '''//trim(dirn)//''': use both, forward or backward'
   return
endif
if (ds0 <= 0.0d0) then
   write(out_unit,'(A)') '[IRCFAILED]  irc_step must be positive'
   return
endif
if (spec%npc > 0) then
   write(out_unit,'(A)') '[IRCFAILED]  external point charges are not supported with irc'
   return
endif
do_fwd = (trim(dirn) == 'both' .or. trim(dirn) == 'forward')
do_bwd = (trim(dirn) == 'both' .or. trim(dirn) == 'backward')

engine_quiet = .true.
force_dense_allow = .false.

allocate(atomchg_a(nat), coord_a(nat,3), force_out(nat,3), MLcharge_out(nat))
allocate(pc_charge_a(1), pc_coord_a(1,3))
allocate(dens_in_a(0,0), dens_in_b(0,0))
atomchg_a = spec%atomchg
pc_charge_a = 0.0d0 ; pc_coord_a = 0.0d0

allocate(x_ts(3,nat), g_ts(n), Hy_ts(n,n), v(n), sm(n))
allocate(lamv(n), Uv(n,n))
allocate(Xf(n,maxp+1), Ef(maxp+1), Sf(maxp+1), Gf(maxp+1))
allocate(Xb(n,maxp+1), Eb(maxp+1), Sb(maxp+1), Gb(maxp+1))
nf = 0 ; nb = 0
reason_f = '' ; reason_b = ''

do i = 1, nat
   sm(3*(i-1)+1:3*(i-1)+3) = sqrt(atomic_mass(spec%atomchg(i)))
   x_ts(1,i) = spec%x(i) * ANG2BOHR
   x_ts(2,i) = spec%y(i) * ANG2BOHR
   x_ts(3,i) = spec%z(i) * ANG2BOHR
enddo

write(out_unit,'(A)') '[IRC]'
write(out_unit,'(A,A)')       '  direction      = ', trim(dirn)
write(out_unit,'(A,F8.4,A)')  '  step           = ', ds0, '  (amu^1/2 bohr)'
write(out_unit,'(A,I0)')      '  max points     = ', maxp
write(out_unit,'(A,I0)')      '  Hessian every  = ', spec%irc_hess_every
write(out_unit,'(A,ES10.2)')  '  gradient stop  = ', spec%irc_stop(1)
if (spec%irc_stop(2) > 0.0d0) write(out_unit,'(A,F7.4,A)') '  bond-order stop = ', spec%irc_stop(2), &
     '  (max Mayer change over 3 points)'
write(out_unit,'(A)') '[IRCEND]'
write(out_unit,'(A)') ''

call eg(x_ts, E_ts, g_ts, ok)
bo_ts = bo_last
if (.not. ok) then
   write(out_unit,'(A)') '[IRCFAILED]  SCF did not converge at the transition state'
   call restore_state()
   return
endif
block
   real(8), allocatable :: Hc(:,:)
   allocate(Hc(n,n))
   call hes(x_ts, Hc, ok)
   if (.not. ok) then
      write(out_unit,'(A)') '[IRCFAILED]  the analytic Hessian at the transition state could not be obtained'
      call restore_state()
      return
   endif
   call to_mass_weighted_hessian(Hc, Hy_ts)
   call project_hessian(x_ts, Hy_ts)
end block

Uv = Hy_ts
call eigen(n, Uv, lamv, ok)
if (.not. ok) then
   write(out_unit,'(A)') '[IRCFAILED]  Hessian diagonalization failed'
   call restore_state()
   return
endif
nneg = count(lamv < NEG_EIG)
imode = 1
lam = lamv(1)
freq_cm = sign(sqrt(abs(lam)/AMU_TO_ME) * AU_TO_CM, lam)
write(out_unit,'(A)') '[IRCSTART]'
write(out_unit,'(A,F18.9)')  '  E(TS)                = ', E_ts
write(out_unit,'(A,F10.2,A)') '  lowest mode          = ', freq_cm, ' cm-1  (negative = imaginary)'
write(out_unit,'(A,I0)')     '  imaginary modes      = ', nneg
write(out_unit,'(A,ES10.2)') '  gradient max (TS)    = ', maxval(abs(g_ts))
write(out_unit,'(A)') '[IRCSTARTEND]'
write(out_unit,'(A)') ''
if (nneg /= 1) then
   write(out_unit,'(A,I0,A)') '[IRCFAILED]  the structure has ', nneg, &
        ' imaginary modes, a transition state has exactly one - optimize it with opt_ts first'
   call restore_state()
   return
endif
v = Uv(:, imode)
i = maxloc(abs(v), 1)
if (v(i) < 0.0d0) v = -v

if (do_fwd) call walk(+1.0d0, 'forward ', Xf, Ef, Sf, Gf, nf, reason_f)
if (do_bwd) call walk(-1.0d0, 'backward', Xb, Eb, Sb, Gb, nb, reason_b)

call write_summary()
call write_path_files()
call restore_state()

contains

   subroutine restore_state()
   engine_quiet = .false.
   mayer_request = .false.
   end subroutine restore_state

   subroutine eg(x3, E, g, okk)
   real(8), intent(in)  :: x3(3,nat)
   real(8), intent(out) :: E, g(n)
   logical, intent(out) :: okk
   integer :: k
   do k = 1, nat
      coord_a(k,1) = x3(1,k) * BOHR2ANG
      coord_a(k,2) = x3(2,k) * BOHR2ANG
      coord_a(k,3) = x3(3,k) * BOHR2ANG
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
                 mem_grid_gb, mem_2e_gb, max(spec%scf_conv_level, 1), trim(spec%basedir), &
                 dens_in_a, dens_in_b, dens_out_a, dens_out_b)
   n_eg_calls = n_eg_calls + 1
   okk = (iconv == 1)
   if (.not. okk) return
   E = energy_out
   g = reshape(transpose(force_out), [n])
   if (allocated(mayer_bo)) then
      block
         integer :: ia, ib, ip
         ip = 0
         do ia = 1, nat
            do ib = ia+1, nat
               ip = ip + 1
               bo_last(ip) = mayer_bo(ia,ib)
            enddo
         enddo
      end block
   endif
   if (allocated(dens_out_a) .and. allocated(dens_out_b)) then
      call move_alloc(dens_out_a, dens_in_a)
      call move_alloc(dens_out_b, dens_in_b)
   endif
   end subroutine eg

   subroutine hes(x3, H, okk)
   real(8), intent(in)  :: x3(3,nat)
   real(8), intent(out) :: H(n,n)
   logical, intent(out) :: okk
   n_hess_calls = n_hess_calls + 1
   call hessian_analytic(spec, x3, nat, out_unit, max(spec%scf_conv_level, 1), &
                         dens_in_a, dens_in_b, dens_out_a, dens_out_b, H, okk)
   engine_quiet = .true.
   end subroutine hes

   subroutine to_mass_weighted_hessian(Hc, Hy)
   real(8), intent(in)  :: Hc(n,n)
   real(8), intent(out) :: Hy(n,n)
   integer :: a, b
   do b = 1, n
      do a = 1, n
         Hy(a,b) = 0.5d0 * (Hc(a,b) + Hc(b,a)) / (sm(a) * sm(b))
      enddo
   enddo
   end subroutine to_mass_weighted_hessian

   subroutine tr_basis(x3, tr, ntr_)
   real(8), intent(in)  :: x3(3,nat)
   real(8), intent(out) :: tr(n,6)
   integer, intent(out) :: ntr_
   real(8) :: com(3), mtot, raw(n,6), d(3), w, nrm, dotv
   integer :: a, ax, k, m
   com = 0.0d0 ; mtot = 0.0d0
   do a = 1, nat
      w = sm(3*(a-1)+1)**2
      com = com + w * x3(:,a)
      mtot = mtot + w
   enddo
   com = com / mtot
   raw = 0.0d0
   do ax = 1, 3
      do a = 1, nat
         raw(3*(a-1)+ax, ax) = sm(3*(a-1)+ax)
      enddo
   enddo
   do a = 1, nat
      w = sm(3*(a-1)+1)
      d = x3(:,a) - com
      raw(3*(a-1)+2, 4) = -w*d(3) ; raw(3*(a-1)+3, 4) =  w*d(2)
      raw(3*(a-1)+1, 5) =  w*d(3) ; raw(3*(a-1)+3, 5) = -w*d(1)
      raw(3*(a-1)+1, 6) = -w*d(2) ; raw(3*(a-1)+2, 6) =  w*d(1)
   enddo
   ntr_ = 0
   do k = 1, 6
      do m = 1, ntr_
         dotv = dot_product(raw(:,k), tr(:,m))
         raw(:,k) = raw(:,k) - dotv * tr(:,m)
      enddo
      nrm = sqrt(sum(raw(:,k)*raw(:,k)))
      if (nrm > 1.0d-8) then
         ntr_ = ntr_ + 1
         tr(:,ntr_) = raw(:,k) / nrm
      endif
   enddo
   end subroutine tr_basis

   subroutine project_vec(x3, vec)
   real(8), intent(in)    :: x3(3,nat)
   real(8), intent(inout) :: vec(n)
   real(8) :: tr(n,6)
   integer :: m, nt
   call tr_basis(x3, tr, nt)
   do m = 1, nt
      vec = vec - dot_product(tr(:,m), vec) * tr(:,m)
   enddo
   end subroutine project_vec

   subroutine project_hessian(x3, Hy)
   real(8), intent(in)    :: x3(3,nat)
   real(8), intent(inout) :: Hy(n,n)
   real(8) :: tr(n,6), P(n,n), T(n,n)
   integer :: m, nt, a
   call tr_basis(x3, tr, nt)
   P = 0.0d0
   do a = 1, n
      P(a,a) = 1.0d0
   enddo
   do m = 1, nt
      call dger(n, n, -1.0d0, tr(:,m), 1, tr(:,m), 1, P, n)
   enddo
   call dgemm('N', 'N', n, n, n, 1.0d0, P, n, Hy, n, 0.0d0, T, n)
   call dgemm('N', 'N', n, n, n, 1.0d0, T, n, P, n, 0.0d0, Hy, n)
   ntr = nt
   end subroutine project_hessian

   subroutine lqa_step(H, g, ds, dy, reached_min)
   real(8), intent(in)  :: H(n,n), g(n), ds
   real(8), intent(out) :: dy(n)
   logical, intent(out) :: reached_min
   real(8) :: U(n,n), lamL(n), gp(n), gn, dt, t, s, w, dsl, frac, Ri
   integer, parameter :: NEULER = 5000, MAXDBL = 40
   integer :: it, k, ndbl
   logical :: okk
   dy = 0.0d0
   reached_min = .false.
   gn = sqrt(sum(g*g))
   if (gn < 1.0d-12) then
      reached_min = .true.
      return
   endif
   U = H
   call eigen(n, U, lamL, okk)
   if (.not. okk) then
      reached_min = .true.
      return
   endif
   call dgemv('T', n, n, 1.0d0, U, n, g, 1, 0.0d0, gp, 1)
   dt = ds / (NEULER * gn)
   t = 0.0d0 ; s = 0.0d0 ; ndbl = 0
   reached_min = .true.
   do it = 1, NEULER * (MAXDBL + 1)
      w = 0.0d0
      do k = 1, n
         w = w + gp(k)*gp(k) * exp(min(-2.0d0*lamL(k)*t, 600.0d0))
      enddo
      dsl = sqrt(w) * dt
      if (s + dsl >= ds) then
         frac = (ds - s) / dsl
         t = t + frac * dt
         s = ds
         reached_min = .false.
         exit
      endif
      s = s + dsl
      t = t + dt
      if (sqrt(w) < 1.0d-14) exit
      if (mod(it, NEULER) == 0) then
         ndbl = ndbl + 1
         if (ndbl > MAXDBL) exit
         dt = dt * 2.0d0
      endif
   enddo
   do k = 1, n
      if (abs(lamL(k)*t) < 1.0d-8) then
         Ri = -t * (1.0d0 - 0.5d0*lamL(k)*t)
      else
         Ri = (exp(min(-lamL(k)*t, 300.0d0)) - 1.0d0) / lamL(k)
      endif
      gp(k) = Ri * gp(k)
   enddo
   call dgemv('N', n, n, 1.0d0, U, n, gp, 1, 0.0d0, dy, 1)
   end subroutine lqa_step

   subroutine dwi_eval(y, y1, E1, g1, H1, y2, E2, g2, H2, E, g)
   real(8), intent(in)  :: y(n), y1(n), E1, g1(n), H1(n,n), y2(n), E2, g2(n), H2(n,n)
   real(8), intent(out) :: E, g(n)
   real(8) :: d1(n), d2(n), Hd1(n), Hd2(n), T1, T2, gT1(n), gT2(n), r1, r2, ssum, w1, w2, gw1(n)
   d1 = y - y1 ; d2 = y - y2
   call dgemv('N', n, n, 1.0d0, H1, n, d1, 1, 0.0d0, Hd1, 1)
   call dgemv('N', n, n, 1.0d0, H2, n, d2, 1, 0.0d0, Hd2, 1)
   T1 = E1 + dot_product(g1, d1) + 0.5d0*dot_product(d1, Hd1)
   T2 = E2 + dot_product(g2, d2) + 0.5d0*dot_product(d2, Hd2)
   gT1 = g1 + Hd1 ; gT2 = g2 + Hd2
   r1 = dot_product(d1, d1) ; r2 = dot_product(d2, d2) ; ssum = r1 + r2
   if (ssum < 1.0d-24) then
      E = 0.5d0*(T1 + T2) ; g = 0.5d0*(gT1 + gT2)
      return
   endif
   w1 = r2/ssum ; w2 = r1/ssum
   gw1 = (2.0d0*d2*r1 - 2.0d0*d1*r2) / (ssum*ssum)
   E = w1*T1 + w2*T2
   g = w1*gT1 + w2*gT2 + (T1 - T2)*gw1
   end subroutine dwi_eval

   subroutine dwi_euler(nstep, ds, y1, E1, g1, H1, y2, E2, g2, H2, yout)
   integer, intent(in)  :: nstep
   real(8), intent(in)  :: ds, y1(n), E1, g1(n), H1(n,n), y2(n), E2, g2(n), H2(n,n)
   real(8), intent(out) :: yout(n)
   real(8) :: y(n), Ex, gx(n), gnrm, h
   integer :: k
   y = y1 ; h = ds / nstep
   do k = 1, nstep
      call dwi_eval(y, y1, E1, g1, H1, y2, E2, g2, H2, Ex, gx)
      gnrm = sqrt(sum(gx*gx))
      if (gnrm < 1.0d-10) exit
      y = y - h * gx / gnrm
   enddo
   yout = y
   end subroutine dwi_euler

   subroutine dwi_corrector(ds, y1, E1, g1, H1, y2, E2, g2, H2, ynew, err)
   real(8), intent(in)  :: ds, y1(n), E1, g1(n), H1(n,n), y2(n), E2, g2(n), H2(n,n)
   real(8), intent(out) :: ynew(n), err
   integer, parameter :: KMAX = 8
   integer :: nseq(KMAX), k, j
   real(8), allocatable :: Tb(:,:,:)
   real(8) :: fac
   nseq = [2, 4, 6, 8, 10, 12, 14, 16]
   allocate(Tb(n, KMAX, KMAX))
   err = huge(1.0d0)
   do k = 1, KMAX
      call dwi_euler(nseq(k), ds, y1, E1, g1, H1, y2, E2, g2, H2, Tb(:,k,1))
      do j = 2, k
         fac = dble(nseq(k)) / dble(nseq(k-j+1)) - 1.0d0
         Tb(:,k,j) = Tb(:,k,j-1) + (Tb(:,k,j-1) - Tb(:,k-1,j-1)) / fac
      enddo
      ynew = Tb(:,k,k)
      if (k > 1) then
         err = sqrt(sum((Tb(:,k,k) - Tb(:,k-1,k-1))**2))
         if (err < 1.0d-6) exit
      endif
   enddo
   end subroutine dwi_corrector

   subroutine walk(sgn, label, X, Es, Ss, Gm, np, reason)
   real(8), intent(in)  :: sgn
   character(len=*), intent(in) :: label
   real(8), intent(out) :: X(n,maxp+1), Es(maxp+1), Ss(maxp+1), Gm(maxp+1)
   integer, intent(out) :: np
   character(len=*), intent(out) :: reason
   real(8) :: x3(3,nat), x3n(3,nat), gn(n), gy(n), gyn(n), Hy(n,n), Hold(n,n), Hc(n,n)
   real(8) :: yc(n), yp(n), ynew(n), gnew(n), Ecur, Enew
   real(8) :: dy(n), E, En, ds, s_tot, dE_pred, gmaxc, gyp(n), dbo, cerr, arc
   real(8) :: BOH(npair,maxp+1)
   integer :: ip, nretry, since_hess, ncorr_skip
   logical :: okk, atmin, hess_fresh, accepted, use_corr

   use_corr = spec%irc_stop(3) > 0.0d0
   ncorr_skip = 0
   np = 1
   BOH(:,1) = bo_ts
   X(:,1) = reshape(x_ts, [n]) ; Es(1) = E_ts ; Ss(1) = 0.0d0 ; Gm(1) = maxval(abs(g_ts))
   x3 = x_ts ; E = E_ts ; Ecur = E_ts
   yc = reshape(x3, [n]) * sm
   Hy = Hy_ts
   gy = g_ts / sm ; call project_vec(x3, gy)
   s_tot = 0.0d0
   since_hess = 0
   hess_fresh = .true.
   reason = 'max points'

   write(out_unit,'(A)') '[IRCWALK '//trim(adjustl(label))//']'
   write(out_unit,'(A)') '  point        s          E (Ha)       E-E(TS) (kcal/mol)   gmax (Ha/bohr)   step   dBO(3pt)  Hessian'

   do ip = 1, maxp
      ds = ds0
      nretry = 0
      accepted = .false.
      do while (.not. accepted)
         if (ip == 1) then
            dy = sgn * ds * v
            atmin = .false.
         else
            call lqa_step(Hy, gy, ds, dy, atmin)
         endif
         yp = yc + dy
         x3n = reshape(yp / sm, [3, nat])
         call eg(x3n, En, gn, okk)
         if (.not. okk) then
            reason = 'SCF failed'
            goto 900
         endif
         if (En > E + E_RISE) then
            nretry = nretry + 1
            if (nretry > MAXRETRY) then
               reason = 'energy rises'
               goto 900
            endif
            if (.not. hess_fresh) then
               call hes(x3, Hc, okk)
               if (.not. okk) then
                  reason = 'Hessian failed'
                  goto 900
               endif
               call to_mass_weighted_hessian(Hc, Hy)
               call project_hessian(x3, Hy)
               hess_fresh = .true. ; since_hess = 0
            else
               ds = max(0.5d0*ds, DS_MIN)
            endif
            cycle
         endif
         accepted = .true.
      enddo

      gyn = gn / sm
      call project_vec(x3n, gyn)
      dE_pred = dot_product(gy, dy) + 0.5d0 * dot_product(dy, matmul(Hy, dy))
      arc = sqrt(sum(dy*dy))
      s_tot = s_tot + arc
      gyp = gyn * sm
      gmaxc = maxval(abs(gyp))
      dbo = -1.0d0

      Hold = Hy
      call bofill_update(n, Hy, dy, gyn - gy)
      call project_hessian(x3n, Hy)
      hess_fresh = .false.
      since_hess = since_hess + 1

      cerr = 0.0d0
      if (use_corr .and. ip > 1) then
         call dwi_corrector(arc, yc, Ecur, gy, Hold, yp, En, gyn, Hy, ynew, cerr)
         if (cerr > BS_ERR_MAX .or. sqrt(sum((ynew - yp)**2)) > BS_SHIFT_MAX * arc) then
            ynew = yp ; Enew = En ; gnew = gyn
            ncorr_skip = ncorr_skip + 1
         else
            call dwi_eval(ynew, yc, Ecur, gy, Hold, yp, En, gyn, Hy, Enew, gnew)
            call project_vec(reshape(ynew / sm, [3, nat]), gnew)
         endif
      else
         ynew = yp ; Enew = En ; gnew = gyn
      endif

      np = np + 1
      X(:,np) = ynew / sm ; Es(np) = Enew ; Ss(np) = s_tot ; Gm(np) = gmaxc
      BOH(:,np) = bo_last
      if (np - BO_WINDOW >= 1) dbo = maxval(abs(bo_last - BOH(:,np-BO_WINDOW)))
      write(out_unit,'(2X,I5,F12.4,F18.9,F18.4,ES18.3,F9.4,F10.4,3X,A)') ip, s_tot, Enew, (Enew - E_ts)*HARTREE_TO_KCAL, &
           gmaxc, arc, dbo, merge('analytic', 'updated ', hess_fresh)
      if (engine_verbose_irc()) write(out_unit,'(A,ES11.3,A,ES11.3,A,ES10.2,A,ES9.2)') '    predicted dE = ', dE_pred, &
           '   actual dE = ', En - E, '   corrector shift = ', sqrt(sum((ynew - yp)**2)), '   BS error = ', cerr

      yc = ynew ; Ecur = Enew ; gy = gnew
      x3 = reshape(yc / sm, [3, nat]) ; E = En
      if (gmaxc < spec%irc_stop(1)) then
         reason = 'gradient converged'
         exit
      endif
      if (atmin) then
         reason = 'minimum reached'
         exit
      endif
      if (spec%irc_stop(2) > 0.0d0 .and. dbo >= 0.0d0 .and. E_ts - En >= BO_E_MIN .and. gmaxc < BO_GMAX &
          .and. dbo < spec%irc_stop(2)) then
         reason = 'bond orders stable'
         exit
      endif

      if (spec%irc_hess_every > 0 .and. since_hess >= spec%irc_hess_every) then
         call hes(x3n, Hc, okk)
         if (.not. okk) then
            reason = 'Hessian failed'
            exit
         endif
         call to_mass_weighted_hessian(Hc, Hy)
         call project_hessian(x3n, Hy)
         hess_fresh = .true. ; since_hess = 0
      endif
   enddo
900 continue
   if (use_corr .and. ncorr_skip > 0) write(out_unit,'(A,I0,A)') '  corrector skipped at ', ncorr_skip, ' point(s)'
   write(out_unit,'(A)') '  stop: '//trim(reason)
   write(out_unit,'(A)') '[IRCWALKEND]'
   write(out_unit,'(A)') ''
   end subroutine walk

   logical function engine_verbose_irc()
   engine_verbose_irc = spec%verbose >= 2
   end function engine_verbose_irc

   subroutine write_summary()
   real(8) :: Emin_f, Emin_b
   write(out_unit,'(A)') '[IRCRESULT]'
   write(out_unit,'(A,F18.9)') '  E(TS)             = ', E_ts
   if (do_fwd) then
      write(out_unit,'(A,I0,A,A)') '  forward points    = ', nf-1, '   stop: ', trim(reason_f)
      write(out_unit,'(A,F18.9,A,F10.4,A)') '  forward end E     = ', Ef(nf), &
           '   barrier ', (E_ts - Ef(nf))*HARTREE_TO_KCAL, ' kcal/mol'
   endif
   if (do_bwd) then
      write(out_unit,'(A,I0,A,A)') '  backward points   = ', nb-1, '   stop: ', trim(reason_b)
      write(out_unit,'(A,F18.9,A,F10.4,A)') '  backward end E    = ', Eb(nb), &
           '   barrier ', (E_ts - Eb(nb))*HARTREE_TO_KCAL, ' kcal/mol'
   endif
   write(out_unit,'(A,I0,A,I0)') '  energy+gradient calls = ', n_eg_calls, '   analytic Hessians = ', n_hess_calls
   write(out_unit,'(A)') '[IRCRESULTEND]'
   write(out_unit,'(A)') ''
   if (do_fwd .and. nf > 0) then
      if (trim(reason_f) /= 'gradient converged' .and. trim(reason_f) /= 'minimum reached' .and. &
          trim(reason_f) /= 'bond orders stable') &
         write(out_unit,'(A)') '[IRC] forward path did not reach a minimum: raise irc_maxpoints '// &
            'or irc_step, then optimize the end point'
   endif
   if (do_bwd .and. nb > 0) then
      if (trim(reason_b) /= 'gradient converged' .and. trim(reason_b) /= 'minimum reached' .and. &
          trim(reason_b) /= 'bond orders stable') &
         write(out_unit,'(A)') '[IRC] backward path did not reach a minimum: raise irc_maxpoints '// &
            'or irc_step, then optimize the end point'
   endif
   end subroutine write_summary

   subroutine write_path_files()
   character(len=256) :: base, fx, fd
   integer :: l, u, ud, k, ntot
   real(8) :: s_here
   l = len_trim(outfile)
   if (l > 4 .and. outfile(max(1,l-3):l) == '.out') then
      base = outfile(1:l-4)
   else
      base = trim(outfile)
   endif
   fx = trim(base)//'.irc.xyz'
   fd = trim(base)//'.irc.dat'
   call remove_file(fx) ; call remove_file(fd)
   open(newunit=u, file=trim(fx), status='replace', action='write')
   open(newunit=ud, file=trim(fd), status='replace', action='write')
   write(ud,'(A)') '# s (amu^1/2 bohr, negative = backward)   E (Ha)   E-E(TS) (kcal/mol)   gmax (Ha/bohr)'
   ntot = 0
   do k = nb, 2, -1
      s_here = -Sb(k)
      call frame(u, ud, Xb(:,k), Eb(k), s_here, Gb(k)) ; ntot = ntot + 1
   enddo
   call frame(u, ud, reshape(x_ts, [n]), E_ts, 0.0d0, maxval(abs(g_ts))) ; ntot = ntot + 1
   do k = 2, nf
      call frame(u, ud, Xf(:,k), Ef(k), Sf(k), Gf(k)) ; ntot = ntot + 1
   enddo
   close(u) ; close(ud)
   write(out_unit,'(A)') '[IRCFILES]'
   write(out_unit,'(A,A,A,I0,A)') '  path        = ', trim(fx), '  (', ntot, ' frames, backward end to forward end)'
   write(out_unit,'(A,A)')        '  profile     = ', trim(fd)
   write(out_unit,'(A)') '[IRCFILESEND]'
   write(out_unit,'(A)') ''
   end subroutine write_path_files

   subroutine frame(u, ud, xv, E, s, gm)
   integer, intent(in) :: u, ud
   real(8), intent(in) :: xv(n), E, s, gm
   integer :: a
   write(u,'(I0)') nat
   write(u,'(A,F10.4,A,F18.9)') 's = ', s, '  E = ', E
   do a = 1, nat
      write(u,'(A2,3F16.9)') elem_to_name(spec%atomchg(a)), xv(3*(a-1)+1)*BOHR2ANG, &
           xv(3*(a-1)+2)*BOHR2ANG, xv(3*(a-1)+3)*BOHR2ANG
   enddo
   write(ud,'(F12.4,F20.9,F14.4,ES14.3)') s, E, (E - E_ts)*HARTREE_TO_KCAL, gm
   end subroutine frame

end subroutine irc_run

subroutine eigen(n, A, w, ok)
integer, intent(in)    :: n
real(8), intent(inout) :: A(n,n)
real(8), intent(out)   :: w(n)
logical, intent(out)   :: ok
real(8), allocatable :: work(:)
integer :: lwork, info
lwork = max(1, 66*n)
allocate(work(lwork))
call dsyev('V', 'U', n, A, n, w, work, lwork, info)
ok = (info == 0)
end subroutine eigen

function lower16(s) result(r)
character(len=*), intent(in) :: s
character(len=16) :: r
integer :: i, c
r = ''
do i = 1, min(len_trim(s), 16)
   c = ichar(s(i:i))
   if (c >= ichar('A') .and. c <= ichar('Z')) c = c + 32
   r(i:i) = char(c)
enddo
r = adjustl(r)
end function lower16

end module mod_irc
