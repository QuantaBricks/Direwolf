! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Fragment embedding for a cluster or complex (opt embed = true, needs

module mod_opt_embed
use mod_engineup_interface, only: EngineUp
use mod_engine_options, only: engine_options_t
use mod_engine_input_types, only: engine_input_t
use mod_opt_fragment, only: frag_get_fid, frag_nfrag, lj_params, frag_classical_opt
implicit none
private
public :: embed_init, embed_evaluate, embed_reset, embed_prerelax, embed_prepare_micro

real(8), parameter :: BOHR2ANG = 0.529177210903d0
real(8), parameter :: KCAL2HA = 1.0d0/627.5094740631d0
real(8), parameter :: Q_TOL = 1.0d-3
real(8), parameter :: Q_TOL_CARRY = 1.0d-3
integer, parameter :: SWEEP_MAX = 5
logical, parameter :: USE_DH = .false.

type :: dens_t
   real(8), allocatable :: a(:,:), b(:,:)
end type dens_t

integer :: nf = 0, nat_e = 0
integer, allocatable :: fid(:), nfa(:), fchg(:), fmult(:)
integer, allocatable :: flist(:,:)
real(8), allocatable :: qchg(:)
logical :: have_q = .false.
logical :: q_frozen_on = .false.
real(8), allocatable :: dg(:,:), xref(:,:)
real(8), allocatable :: dH(:,:)
real(8), allocatable :: dg_prev_f(:), xref_prev_f(:)
logical :: have_macro_prev = .false.
type(dens_t), allocatable :: dens(:)

contains

subroutine embed_reset()
nf = 0 ; nat_e = 0 ; have_q = .false. ; q_frozen_on = .false.
if (allocated(dg)) deallocate(dg)
if (allocated(xref)) deallocate(xref)
if (allocated(dH)) deallocate(dH)
if (allocated(dg_prev_f)) deallocate(dg_prev_f)
if (allocated(xref_prev_f)) deallocate(xref_prev_f)
have_macro_prev = .false.
if (allocated(fid)) deallocate(fid)
if (allocated(nfa)) deallocate(nfa)
if (allocated(fchg)) deallocate(fchg)
if (allocated(fmult)) deallocate(fmult)
if (allocated(flist)) deallocate(flist)
if (allocated(qchg)) deallocate(qchg)
if (allocated(dens)) deallocate(dens)
end subroutine embed_reset

subroutine embed_init(spec, nat, out_unit, ok, show)
type(engine_input_t), intent(in) :: spec
integer, intent(in)  :: nat, out_unit
logical, intent(out) :: ok
logical, intent(in)  :: show
integer :: f, i, k, ne, ios
character(len=256) :: s
call embed_reset()
ok = .false.
nf = frag_nfrag()
if (nf < 2) then
   write(out_unit,'(A)') '[EMBED] ERROR: embed needs at least two covalent fragments'
   return
endif
nat_e = nat
allocate(fid(nat), nfa(nf), fchg(nf), fmult(nf), qchg(nat), dens(nf))
call frag_get_fid(nat, fid)
nfa = 0
do i = 1, nat
   nfa(fid(i)) = nfa(fid(i)) + 1
enddo
allocate(flist(maxval(nfa), nf))
flist = 0 ; nfa = 0
do i = 1, nat
   nfa(fid(i)) = nfa(fid(i)) + 1
   flist(nfa(fid(i)), fid(i)) = i
enddo
fchg = 0 ; fmult = 1
s = spec%opt_frag_charges
if (len_trim(s) > 0) then
   call tolist(s)
   read(s, *, iostat=ios) fchg
   if (ios /= 0) then
      write(out_unit,'(A,I0,A)') '[EMBED] ERROR: opt_frag_charges needs ', nf, ' integers'
      return
   endif
endif
s = spec%opt_frag_mults
if (len_trim(s) > 0) then
   call tolist(s)
   read(s, *, iostat=ios) fmult
   if (ios /= 0) then
      write(out_unit,'(A,I0,A)') '[EMBED] ERROR: opt_frag_mults needs ', nf, ' integers'
      return
   endif
endif
if (sum(fchg) /= spec%icharge) then
   write(out_unit,'(A,I0,A,I0)') '[EMBED] ERROR: fragment charges add up to ', sum(fchg), ', total charge is ', spec%icharge
   return
endif
do f = 1, nf
   ne = -fchg(f)
   do k = 1, nfa(f)
      ne = ne + spec%atomchg(flist(k,f))
   enddo
   if (mod(ne, 2) /= mod(fmult(f) - 1, 2)) then
      write(out_unit,'(A,I0,A,I0,A,I0,A)') '[EMBED] ERROR: fragment ', f, ' has ', ne, &
           ' electrons, multiplicity ', fmult(f), ' does not fit'
      return
   endif
   if (fmult(f) /= 1) then
      write(out_unit,'(A,I0,A)') '[EMBED] ERROR: fragment ', f, ' is open-shell; only closed-shell fragments are supported'
      return
   endif
enddo
qchg = 0.0d0 ; have_q = .false.
if (.not. show) then
   ok = .true.
   return
endif
write(out_unit,'(A)') '[EMBED]  every fragment is computed alone in the field of the others'' ESP charges (+ LJ)'
write(out_unit,'(A,I0,A,I0,A)') '  fragments = ', nf, '   atoms = ', nat, '   (no calculation of the whole system)'
do f = 1, nf
   write(out_unit,'(A,I0,A,I0,A,I0,A,I0)') '  fragment ', f, ':  atoms ', nfa(f), '   charge ', fchg(f), '   mult ', fmult(f)
enddo
write(out_unit,'(A)') '[EMBEDEND]'
ok = .true.
contains
   subroutine tolist(str)
   character(len=*), intent(inout) :: str
   integer :: c
   do c = 1, len_trim(str)
      if (str(c:c) == ',' .or. str(c:c) == ';') str(c:c) = ' '
   enddo
   end subroutine tolist
end subroutine embed_init

subroutine embed_evaluate(spec, opts, coord_a, atomchg_a, scf_lvl, cyc, out_unit, E_tot, force_out, iconv)
type(engine_input_t), intent(in)    :: spec
type(engine_options_t), intent(inout) :: opts
integer, intent(in)  :: atomchg_a(nat_e), scf_lvl, cyc, out_unit
real(8), intent(in)  :: coord_a(nat_e,3)
real(8), intent(out) :: E_tot, force_out(nat_e,3)
integer, intent(out) :: iconv

integer :: f, k, i, j, npcf, sweep, smax_, iconv_f, nn
real(8) :: Ef(nf), Eqm, E_es, E_lj, dq, econv, mem_grid_gb, mem_2e_gb, t_qm
integer(8) :: c0, c1, crate
real(8), allocatable :: qin(:), qnew(:), cf(:,:), fo(:,:), mlq(:), pcq(:), pcx(:,:), flj(:,:), fes(:,:), fown(:,:), fpc(:,:)
integer, allocatable :: zf(:)
character(len=30), allocatable :: bas(:), ecp(:)
real(8), allocatable :: da(:,:), db(:,:), oa(:,:), ob(:,:)
integer :: npc_ext, uq
logical :: qfrozen
integer, allocatable :: pcatom(:)

allocate(qin(nat_e), qnew(nat_e), flj(nat_e,3), fes(nat_e,3), fown(nat_e,3), fpc(nat_e,3))
qfrozen = q_frozen_on
if (spec%verbose >= 3 .and. .not. q_frozen_on) then
   inquire(file='embed_q.dat', exist=qfrozen)
   if (qfrozen) then
      open(newunit=uq, file='embed_q.dat', status='old', action='read')
      read(uq, *) qchg
      close(uq)
   endif
endif
E_tot = 0.0d0 ; force_out = 0.0d0
npc_ext = spec%npc
qin = qchg
smax_ = SWEEP_MAX
iconv = 1
t_qm = 0.0d0
opts%esp_charges = .not. qfrozen
opts%pc_grad_request = .true.
do sweep = 1, smax_
   Ef = 0.0d0 ; force_out = 0.0d0 ; qnew = 0.0d0 ; fown = 0.0d0 ; fpc = 0.0d0
   do f = 1, nf
      nn = nfa(f)
      allocate(cf(nn,3), fo(nn,3), mlq(nn), zf(nn), bas(nn), ecp(nn))
      do k = 1, nn
         i = flist(k,f)
         cf(k,:) = coord_a(i,:) ; zf(k) = atomchg_a(i)
         bas(k) = spec%atom_basis(i) ; ecp(k) = spec%atom_ecp(i)
      enddo
      npcf = (nat_e - nn) + npc_ext
      allocate(pcq(max(npcf,1)), pcx(max(npcf,1),3), pcatom(max(npcf,1)))
      pcatom = 0
      j = 0
      do i = 1, nat_e
         if (fid(i) == f) cycle
         j = j + 1 ; pcq(j) = qin(i) ; pcx(j,:) = coord_a(i,:) ; pcatom(j) = i
      enddo
      do i = 1, npc_ext
         j = j + 1 ; pcq(j) = spec%pc_q(i)
         pcx(j,1) = spec%pc_x(i) ; pcx(j,2) = spec%pc_y(i) ; pcx(j,3) = spec%pc_z(i)
      enddo
      if (allocated(dens(f)%a)) call move_alloc(dens(f)%a, da)
      if (allocated(dens(f)%b)) call move_alloc(dens(f)%b, db)
      call system_clock(c0, crate)
      call EngineUp(nn, fmult(f), fchg(f), spec%functional, &
                    cf, zf, spec%baselabel, spec%ecplabel, bas, ecp, &
                    trim(spec%J), trim(spec%K), trim(spec%ri_aux_basis), spec%spherical, spec%harris_guess, &
                    .true., spec%vv10_nonself, spec%mem_cap_gb, &
                    .false., spec%n_threads, &
                    npcf, pcq, pcx, &
                    .false., '', .false., '', &
                    .false., spec%cosmo_epsilon, spec%cosmo_radii_scale, spec%cosmo_avg_area, spec%cosmo_sigma_rav, &
                    trim(spec%cosmo_cavity_type), spec%cosmo_rsolv, spec%cosmo_ks_nseg, spec%cosmo_ks_nface, &
                    '', .false., '', &
                    fo, Ef(f), mlq, iconv_f, econv, &
                    mem_grid_gb, mem_2e_gb, scf_lvl, trim(spec%basedir), &
                    da, db, oa, ob, opts)
      call system_clock(c1) ; t_qm = t_qm + real(c1 - c0, 8) / real(crate, 8)
      if (allocated(oa)) call move_alloc(oa, dens(f)%a)
      if (allocated(ob)) call move_alloc(ob, dens(f)%b)
      if (allocated(da)) deallocate(da)
      if (allocated(db)) deallocate(db)
      if (iconv_f /= 1) iconv = 0
      do k = 1, nn
         i = flist(k,f)
         force_out(i,:) = force_out(i,:) + fo(k,:) ; fown(i,:) = fo(k,:)
         if (allocated(opts%esp_q)) qnew(i) = opts%esp_q(k)
      enddo
      if (allocated(opts%pc_grad)) then
         do j = 1, npcf
            if (pcatom(j) > 0) then
               force_out(pcatom(j),:) = force_out(pcatom(j),:) + opts%pc_grad(j,:)
               fpc(pcatom(j),:) = fpc(pcatom(j),:) + opts%pc_grad(j,:)
            endif
         enddo
      endif
      deallocate(cf, fo, mlq, zf, bas, ecp, pcq, pcx, pcatom)
   enddo
   if (iconv /= 1) return
   dq = maxval(abs(qnew - qin))
   E_es = es_energy(atomchg_a, coord_a, qin)
   qchg = qnew
   if (qfrozen) then
      qchg = qin
      exit
   endif
   if (dq < Q_TOL .or. (cyc > 1 .and. dq < Q_TOL_CARRY)) exit
   if (sweep < smax_) qin = qnew
enddo
Eqm = sum(Ef)
if (spec%verbose >= 3 .and. .not. qfrozen) then
   open(newunit=uq, file='embed_q.dat', status='replace', action='write')
   write(uq, '(ES24.16)') qchg
   close(uq)
endif
call lj_term(atomchg_a, coord_a, E_lj, flj)
call es_grad(coord_a, qin, fes)
force_out = force_out + flj - fes
E_tot = Eqm - E_es + E_lj
if (q_frozen_on .and. allocated(dg)) then
   block
      real(8) :: dxf(3*nat_e), hd(3*nat_e)
      integer :: ia
      do ia = 1, nat_e
         dxf(3*(ia-1)+1:3*(ia-1)+3) = (coord_a(ia,:) - xref(ia,:)) / BOHR2ANG
      enddo
      hd = 0.0d0
      if (allocated(dH)) call dgemv('N', 3*nat_e, 3*nat_e, 1.0d0, dH, 3*nat_e, dxf, 1, 0.0d0, hd, 1)
      do ia = 1, nat_e
         force_out(ia,:) = force_out(ia,:) + dg(ia,:) + hd(3*(ia-1)+1:3*(ia-1)+3)
         E_tot = E_tot + dot_product(dg(ia,:), dxf(3*(ia-1)+1:3*(ia-1)+3)) &
                       + 0.5d0 * dot_product(hd(3*(ia-1)+1:3*(ia-1)+3), dxf(3*(ia-1)+1:3*(ia-1)+3))
      enddo
   end block
endif
have_q = .true.
write(out_unit,'(A,I0,A)') '[EMBED cycle ', cyc, ']'
write(out_unit,'(A,I0,A,ES9.2,A,I0,A,F7.2,A)') '  sweeps = ', sweep, '   max |dq| = ', dq, &
     '   fragment calculations = ', sweep*nf, '   wall = ', t_qm, ' s'
write(out_unit,'(A,F18.9,A,F14.9,A,F14.9)') '  E_fragments = ', Eqm, '   E_charge-charge = ', E_es, '   E_LJ = ', E_lj
if (spec%verbose >= 2) then
   write(out_unit,'(A,10F16.9)') '  E per fragment: ', Ef
   write(out_unit,'(A)') '  atom: own-fragment force | point-charge force from the other fragments | classical es | LJ'
   do i = 1, nat_e
      write(out_unit,'(I5,12F14.8)') i, fown(i,:), fpc(i,:), fes(i,:), flj(i,:)
   enddo
   write(out_unit,'(A)') '  dE/dR (Ha/bohr) per atom: x y z'
   do i = 1, nat_e
      write(out_unit,'(I5,3F16.10)') i, force_out(i,:)
   enddo
endif
write(out_unit,'(A)') '[EMBEDEND]'
end subroutine embed_evaluate

real(8) function es_energy(Z, coord_a, q)
integer, intent(in) :: Z(nat_e)
real(8), intent(in) :: coord_a(nat_e,3), q(nat_e)
integer :: i, j
real(8) :: r
es_energy = 0.0d0
do i = 1, nat_e
   do j = i+1, nat_e
      if (fid(i) == fid(j)) cycle
      r = sqrt(sum((coord_a(i,:) - coord_a(j,:))**2)) / BOHR2ANG
      es_energy = es_energy + q(i)*q(j)/r
   enddo
enddo
end function es_energy

subroutine embed_charges(spec, opts, coord_a, atomchg_a, scf_lvl, qin, qout, iconv)
type(engine_input_t), intent(in)      :: spec
type(engine_options_t), intent(inout) :: opts
integer, intent(in)  :: atomchg_a(nat_e), scf_lvl
real(8), intent(in)  :: coord_a(nat_e,3), qin(nat_e)
real(8), intent(out) :: qout(nat_e)
integer, intent(out) :: iconv
integer :: f, k, i, j, npcf, iconv_f, nn
real(8) :: Ef, econv, mem_grid_gb, mem_2e_gb
real(8), allocatable :: cf(:,:), fo(:,:), mlq(:), pcq(:), pcx(:,:)
integer, allocatable :: zf(:)
character(len=30), allocatable :: bas(:), ecp(:)
real(8), allocatable :: da(:,:), db(:,:), oa(:,:), ob(:,:)
iconv = 1 ; qout = 0.0d0
opts%esp_charges = .true.
opts%pc_grad_request = .false.
do f = 1, nf
   nn = nfa(f)
   allocate(cf(nn,3), fo(nn,3), mlq(nn), zf(nn), bas(nn), ecp(nn))
   do k = 1, nn
      i = flist(k,f)
      cf(k,:) = coord_a(i,:) ; zf(k) = atomchg_a(i)
      bas(k) = spec%atom_basis(i) ; ecp(k) = spec%atom_ecp(i)
   enddo
   npcf = (nat_e - nn) + spec%npc
   allocate(pcq(max(npcf,1)), pcx(max(npcf,1),3))
   j = 0
   do i = 1, nat_e
      if (fid(i) == f) cycle
      j = j + 1 ; pcq(j) = qin(i) ; pcx(j,:) = coord_a(i,:)
   enddo
   do i = 1, spec%npc
      j = j + 1 ; pcq(j) = spec%pc_q(i)
      pcx(j,1) = spec%pc_x(i) ; pcx(j,2) = spec%pc_y(i) ; pcx(j,3) = spec%pc_z(i)
   enddo
   allocate(da(0,0), db(0,0))
   call EngineUp(nn, fmult(f), fchg(f), spec%functional, &
                 cf, zf, spec%baselabel, spec%ecplabel, bas, ecp, &
                 trim(spec%J), trim(spec%K), trim(spec%ri_aux_basis), spec%spherical, spec%harris_guess, &
                 .false., spec%vv10_nonself, spec%mem_cap_gb, &
                 .false., spec%n_threads, &
                 npcf, pcq, pcx, &
                 .false., '', .false., '', &
                 .false., spec%cosmo_epsilon, spec%cosmo_radii_scale, spec%cosmo_avg_area, spec%cosmo_sigma_rav, &
                 trim(spec%cosmo_cavity_type), spec%cosmo_rsolv, spec%cosmo_ks_nseg, spec%cosmo_ks_nface, &
                 '', .false., '', &
                 fo, Ef, mlq, iconv_f, econv, &
                 mem_grid_gb, mem_2e_gb, scf_lvl, trim(spec%basedir), &
                 da, db, oa, ob, opts)
   if (allocated(oa)) deallocate(oa)
   if (allocated(ob)) deallocate(ob)
   if (allocated(da)) deallocate(da)
   if (allocated(db)) deallocate(db)
   if (iconv_f /= 1) iconv = 0
   if (iconv_f == 1 .and. allocated(opts%esp_q)) then
      do k = 1, nn
         qout(flist(k,f)) = opts%esp_q(k)
      enddo
   endif
   deallocate(cf, fo, mlq, zf, bas, ecp, pcq, pcx)
enddo
end subroutine embed_charges

subroutine embed_prerelax(spec, opts, atomchg_a, scf_lvl, nat, g3, out_unit)
type(engine_input_t), intent(in)      :: spec
type(engine_options_t), intent(inout) :: opts
integer, intent(in)    :: atomchg_a(:), scf_lvl, nat, out_unit
real(8), intent(inout) :: g3(3,nat)
integer, parameter :: NROUND = 2
integer :: round, iconv, nit, i, f
real(8) :: q(nat), qin(nat), E0, E1, dmax, ca(nat,3), qf
qin = 0.0d0
write(out_unit,'(A)') '[PRERELAX]  rigid-body arrangement of the fragments in a classical model'
write(out_unit,'(A)') '  (ESP charges + LJ), fragments frozen, no DFT calculation of the whole system'
do round = 1, NROUND
   do i = 1, nat
      ca(i,:) = g3(:,i) * BOHR2ANG
   enddo
   call embed_charges(spec, opts, ca, atomchg_a, scf_lvl, qin, q, iconv)
   if (iconv /= 1) then
      write(out_unit,'(A)') '  an SCF did not converge - pre-relaxation stopped'
      exit
   endif
   call frag_classical_opt(nat, atomchg_a, q, g3, 500, E0, E1, nit, dmax)
   write(out_unit,'(A,I0,A,I0,A,F10.5,A,F10.5,A,F7.3,A)') '  round ', round, ':  steps = ', nit, &
        '   model E = ', E0, ' -> ', E1, ' Ha   largest atom shift = ', dmax*BOHR2ANG, ' A'
   qin = q
enddo
do f = 1, nf
   qf = 0.0d0
   do i = 1, nat
      if (fid(i) == f) qf = qf + qin(i)
   enddo
   write(out_unit,'(A,I0,A,F9.4)') '  fragment ', f, ':  sum of ESP charges = ', qf
enddo
write(out_unit,'(A)') '[PRERELAXEND]'
end subroutine embed_prerelax

subroutine embed_prepare_micro(spec, opts, coord_a, atomchg_a, scf_lvl, g_full, out_unit, iconv, g_model_out)
type(engine_input_t), intent(in)      :: spec
type(engine_options_t), intent(inout) :: opts
integer, intent(in)  :: atomchg_a(nat_e), scf_lvl, out_unit
real(8), intent(in)  :: coord_a(nat_e,3), g_full(nat_e,3)
integer, intent(out) :: iconv
real(8), intent(out) :: g_model_out(nat_e,3)
real(8) :: Em, fm(nat_e,3), dmax
q_frozen_on = .false.
if (allocated(dg)) dg = 0.0d0
call embed_evaluate(spec, opts, coord_a, atomchg_a, scf_lvl, 1, out_unit, Em, fm, iconv)
if (iconv /= 1) return
if (.not. allocated(dg)) allocate(dg(nat_e,3))
if (.not. allocated(xref)) allocate(xref(nat_e,3))
dg = g_full - fm
xref = coord_a
block
   real(8) :: sv(3*nat_e), yv(3*nat_e), rv(3*nat_e), den, dgf(3*nat_e), xf(3*nat_e)
   integer :: ia
   do ia = 1, nat_e
      dgf(3*(ia-1)+1:3*(ia-1)+3) = dg(ia,:)
      xf(3*(ia-1)+1:3*(ia-1)+3) = coord_a(ia,:) / BOHR2ANG
   enddo
   if (.not. allocated(dH)) then
      allocate(dH(3*nat_e,3*nat_e)) ; dH = 0.0d0
   endif
   if (have_macro_prev) then
      sv = xf - xref_prev_f ; yv = dgf - dg_prev_f
      call dgemv('N', 3*nat_e, 3*nat_e, 1.0d0, dH, 3*nat_e, sv, 1, 0.0d0, rv, 1)
      rv = yv - rv
      den = dot_product(rv, sv)
      if (USE_DH .and. sqrt(sum(sv*sv)) > 1.0d-3 .and. abs(den) > 1.0d-6 * sqrt(sum(rv*rv)) * sqrt(sum(sv*sv))) then
         call dger(3*nat_e, 3*nat_e, 1.0d0/den, rv, 1, rv, 1, dH, 3*nat_e)
      endif
   endif
   if (.not. allocated(dg_prev_f)) allocate(dg_prev_f(3*nat_e), xref_prev_f(3*nat_e))
   dg_prev_f = dgf ; xref_prev_f = xf ; have_macro_prev = .true.
end block
g_model_out = fm
q_frozen_on = .true.
dmax = maxval(abs(dg))
write(out_unit,'(A,ES10.3,A)') '  [micro] gradient correction: largest component ', dmax, ' Ha/bohr; charges frozen'
end subroutine embed_prepare_micro

subroutine es_grad(coord_a, q, g)
real(8), intent(in)  :: coord_a(nat_e,3), q(nat_e)
real(8), intent(out) :: g(nat_e,3)
integer :: i, j
real(8) :: d(3), r
g = 0.0d0
do i = 1, nat_e
   do j = i+1, nat_e
      if (fid(i) == fid(j)) cycle
      d = (coord_a(i,:) - coord_a(j,:)) / BOHR2ANG
      r = sqrt(sum(d*d))
      g(i,:) = g(i,:) - q(i)*q(j)*d/r**3
      g(j,:) = g(j,:) + q(i)*q(j)*d/r**3
   enddo
enddo
end subroutine es_grad

subroutine lj_term(Z, coord_a, E, grad)
integer, intent(in) :: Z(nat_e)
real(8), intent(in) :: coord_a(nat_e,3)
real(8), intent(out) :: E, grad(nat_e,3)
integer :: i, j
real(8) :: si, ei, sj, ej, sg, ep, d(3), r, s6, dEdr
E = 0.0d0 ; grad = 0.0d0
do i = 1, nat_e
   call lj_params(Z(i), si, ei)
   if (ei <= 0.0d0) cycle
   do j = i+1, nat_e
      if (fid(i) == fid(j)) cycle
      call lj_params(Z(j), sj, ej)
      if (ej <= 0.0d0) cycle
      sg = 0.5d0*(si + sj) ; ep = sqrt(ei*ej) * KCAL2HA
      d = coord_a(i,:) - coord_a(j,:)
      r = sqrt(sum(d*d))
      s6 = (sg/r)**6
      E = E + 4.0d0*ep*(s6*s6 - s6)
      dEdr = 4.0d0*ep*(-12.0d0*s6*s6 + 6.0d0*s6)/r
      grad(i,:) = grad(i,:) + dEdr * d / r
      grad(j,:) = grad(j,:) - dEdr * d / r
   enddo
enddo
grad = grad * BOHR2ANG
end subroutine lj_term

end module mod_opt_embed
