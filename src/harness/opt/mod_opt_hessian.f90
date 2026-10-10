! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Hessian sources for the geometry optimizer: the model Hessian (Lindh,

module mod_opt_hessian
use mod_engine_input_types, only: engine_input_t
use mod_engineup_interface, only: EngineUp
use mod_hessian_freq, only: freq_unit, freq_mult, freq_base, freq_thermo
use mod_integrals, only: freq_request, hessian_return
use mod_engine_input_elements, only: elem_to_name
use mod_internal_coords, only: ic_set_t, ic_nprim, ic_values, ic_bmatrix, ic_ginv, &
                                ic_lindh_hessian
use mod_lindh_hessian, only: lindh_cart_hessian
use mod_opt_min, only: regularize_hessian
implicit none
private
public :: ANG2BOHR, BOHR2ANG
public :: read_ext_hessian, ric_model_hessian, cart_hessian_to_internal, ric_grad_curv
public :: hessian_from_xtb, hessian_analytic
public :: ric_hessian_to_cart, write_hessian_file, remove_file

real(8), parameter :: ANG2BOHR = 1.88972612456506d0
real(8), parameter :: BOHR2ANG = 1.0d0 / ANG2BOHR

contains

subroutine ric_hessian_to_cart(ic, coord3, nat, ni, B, Hq, gx, Hx)
type(ic_set_t), intent(in) :: ic
integer, intent(in) :: nat, ni
real(8), intent(in) :: coord3(3,nat), B(ni,3*nat), Hq(ni,ni), gx(3*nat)
real(8), intent(out) :: Hx(3*nat,3*nat)
real(8), allocatable :: T(:,:), K(:,:)
integer :: n
n = 3*nat
allocate(T(ni,n), K(n,n))
call dgemm('N','N', ni, n, ni, 1.0d0, Hq, ni, B, ni, 0.0d0, T, ni)
call dgemm('T','N', n, n, ni, 1.0d0, B, ni, T, ni, 0.0d0, Hx, n)
call ric_grad_curv(ic, coord3, nat, ni, gx, K)
Hx = Hx + K
end subroutine ric_hessian_to_cart

subroutine write_hessian_file(path, n, H)
character(len=*), intent(in) :: path
integer, intent(in) :: n
real(8), intent(in) :: H(n,n)
integer :: u, ios, i, j0
open(newunit=u, file=trim(path), status='replace', action='write', iostat=ios)
if (ios /= 0) return
write(u,'(A)') '# Direwolf optimizer Hessian, Cartesian, Hartree/bohr^2 (restart)'
do i = 1, n
   do j0 = 1, n, 5
      write(u,'(*(1X,ES24.16E3))') H(i,j0:min(j0+4,n))
   enddo
enddo
close(u)
end subroutine write_hessian_file

subroutine remove_file(path)
character(len=*), intent(in) :: path
integer :: u, ios
open(newunit=u, file=trim(path), status='old', iostat=ios)
if (ios == 0) close(u, status='delete')
end subroutine remove_file

subroutine hessian_from_xtb(spec, g3, nat, outfile, out_unit, raw, Hout, ok)
type(engine_input_t), intent(in) :: spec
integer, intent(in) :: nat, out_unit
real(8), intent(in) :: g3(3,nat)
character(len=*), intent(in) :: outfile
logical, intent(in) :: raw
real(8), intent(out) :: Hout(3*nat,3*nat)
logical, intent(out) :: ok
character(len=1024) :: dir, cmd, meth
integer :: k, u, es, cs
ok = .false.
k = len_trim(outfile)
if (k > 4 .and. outfile(max(k-3,1):k) == '.out') then
   dir = outfile(1:k-4)//'.xtb'
else
   dir = trim(outfile)//'.xtb'
endif
call execute_command_line('mkdir -p "'//trim(dir)//'" && rm -f "'//trim(dir)//'/hessian"', &
                          exitstat=es, cmdstat=cs)
open(newunit=u, file=trim(dir)//'/geom.xyz', status='replace', action='write', iostat=es)
if (es /= 0) then
   write(out_unit,'(A)') '[OPT] xtb Hessian: cannot write '//trim(dir)//'/geom.xyz'
   return
endif
write(u,'(I0)') nat
write(u,'(A)') 'Direwolf geometry for the xtb Hessian'
do k = 1, nat
   write(u,'(A3,3F20.12)') elem_to_name(spec%atomchg(k)), g3(:,k)*BOHR2ANG
enddo
close(u)
select case (trim(spec%opt_xtb_method))
case ('gfnff')
   meth = '--gfnff'
case ('gfn1')
   meth = '--gfn 1'
case ('gfn0')
   meth = '--gfn 0'
case default
   meth = '--gfn 2'
end select
write(cmd,'(A,I0,A,I0,A)') 'cd "'//trim(dir)//'" && '//trim(spec%opt_xtb_cmd)//' geom.xyz --hess '// &
     trim(meth)//' --chrg ', spec%icharge, ' --uhf ', spec%imult-1, ' '//trim(spec%opt_xtb_args)// &
     ' > xtb.out 2>&1'
call execute_command_line(trim(cmd), exitstat=es, cmdstat=cs)
if (cs /= 0 .or. es /= 0) then
   write(out_unit,'(A,I0,A)') '[OPT] xtb Hessian: command failed (exit ', es, '): '//trim(cmd)
   return
endif
call read_ext_hessian(trim(dir)//'/hessian', 3*nat, Hout, ok, out_unit, raw=raw)
if (.not. ok) write(out_unit,'(A)') '[OPT] xtb Hessian: '//trim(dir)//'/hessian not usable'
end subroutine hessian_from_xtb

subroutine hessian_analytic(spec, g3, nat, out_unit, scf_lvl_in, dens_in_a, dens_in_b, &
                            dens_out_a, dens_out_b, Hout, ok)
type(engine_input_t), intent(in) :: spec
integer, intent(in) :: nat, out_unit, scf_lvl_in
real(8), intent(in) :: g3(3,nat)
real(8), allocatable, intent(inout) :: dens_in_a(:,:), dens_in_b(:,:), dens_out_a(:,:), dens_out_b(:,:)
real(8), intent(out) :: Hout(3*nat,3*nat)
logical, intent(out) :: ok
integer, allocatable :: atomchg_a(:)
real(8), allocatable :: coord_a(:,:), force_out(:,:), MLcharge_out(:)
real(8), allocatable :: pc_charge_a(:), pc_coord_a(:,:)
integer :: k, usave, unull, iconv
real(8) :: energy_out, econv, mem_grid_gb, mem_2e_gb
allocate(atomchg_a(nat), coord_a(nat,3), force_out(nat,3), MLcharge_out(nat))
allocate(pc_charge_a(max(spec%npc,1)), pc_coord_a(max(spec%npc,1),3))
atomchg_a = spec%atomchg
pc_charge_a = 0.0d0
pc_coord_a = 0.0d0
do k = 1, spec%npc
   pc_charge_a(k)  = spec%pc_q(k)
   pc_coord_a(k,1) = spec%pc_x(k)
   pc_coord_a(k,2) = spec%pc_y(k)
   pc_coord_a(k,3) = spec%pc_z(k)
enddo
do k = 1, nat
   coord_a(k,1) = g3(1,k) * BOHR2ANG
   coord_a(k,2) = g3(2,k) * BOHR2ANG
   coord_a(k,3) = g3(3,k) * BOHR2ANG
enddo
usave = freq_unit
open(newunit=unull, file='/dev/null', status='old', action='write')
freq_unit = unull
freq_mult = spec%imult
freq_thermo = .false.
freq_request = .true.
hessian_return = .true.
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
              mem_grid_gb, mem_2e_gb, scf_lvl_in, trim(spec%basedir), &
              dens_in_a, dens_in_b, dens_out_a, dens_out_b)
freq_request = .false.
hessian_return = .false.
freq_thermo = .true.
freq_unit = usave
close(unull)
ok = .false.
if (iconv /= 1) return
call read_ext_hessian(trim(freq_base)//'.hess', 3*nat, Hout, ok, out_unit, raw=.true.)
end subroutine hessian_analytic

subroutine read_ext_hessian(path, n, H, ok, out_unit, raw)
character(len=*), intent(in) :: path
integer, intent(in)  :: n, out_unit
real(8), intent(out) :: H(n,n)
logical, intent(out) :: ok
logical, intent(in), optional :: raw
integer :: u, ios, got, i, j
real(8) :: buf(n*n), asym
character(len=512) :: line
logical :: blk
ok = .false.
call read_block_hessian(path, n, H, blk)
if (.not. blk) then
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
endif
asym = maxval(abs(H - transpose(H)))
if (asym > 1.0d-6) write(out_unit,'(A,ES10.3,A)') &
     '[OPT] warning: supplied Hessian is asymmetric by ', asym, ' - symmetrizing'
H = 0.5d0 * (H + transpose(H))
if (.not. present(raw)) then
   call regularize_hessian(H, n, out_unit)
else if (.not. raw) then
   call regularize_hessian(H, n, out_unit)
endif
ok = .true.
end subroutine read_ext_hessian

subroutine read_block_hessian(path, n, H, ok)
character(len=*), intent(in) :: path
integer, intent(in)  :: n
real(8), intent(out) :: H(n,n)
logical, intent(out) :: ok
integer :: u, ios, m, c0, nc, r, k, idx(5)
character(len=512) :: line
real(8) :: v(5)
ok = .false.
H = 0.0d0
open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
if (ios /= 0) return
do
   read(u,'(A)',iostat=ios) line
   if (ios /= 0) then
      close(u) ; return
   endif
   if (index(adjustl(line), '$hessian') == 1) exit
enddo
read(u,*,iostat=ios) m
if (ios /= 0 .or. m /= n) then
   close(u) ; return
endif
c0 = 1
do while (c0 <= n)
   nc = min(5, n - c0 + 1)
   read(u,*,iostat=ios) (idx(k), k = 1, nc)
   if (ios /= 0 .or. idx(1) /= c0 - 1) then
      close(u) ; return
   endif
   do r = 1, n
      read(u,*,iostat=ios) k, (v(m), m = 1, nc)
      if (ios /= 0 .or. k /= r - 1) then
         close(u) ; return
      endif
      H(r, c0:c0+nc-1) = v(1:nc)
   enddo
   c0 = c0 + nc
enddo
close(u)
ok = .true.
end subroutine read_block_hessian

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

subroutine ric_grad_curv(ic, coord3, nat, ni, gx, K)
type(ic_set_t), intent(in) :: ic
integer, intent(in)  :: nat, ni
real(8), intent(in)  :: coord3(3,nat), gx(3*nat)
real(8), intent(out) :: K(3*nat,3*nat)
real(8), parameter :: h = 1.0d-4
real(8), allocatable :: B(:,:), Ginv(:,:), gq(:), bg(:), Bp(:,:), Bm(:,:), c2(:,:)
integer :: ndof, y, a, r
ndof = 3*nat
allocate(B(ni,ndof), Ginv(ni,ni), gq(ni), bg(ni), Bp(ni,ndof), Bm(ni,ndof), c2(3,nat))
call ic_bmatrix(ic, coord3, B)
call ic_ginv(B, ni, ndof, Ginv)
call dgemv('N', ni, ndof, 1.0d0, B, ni, gx, 1, 0.0d0, bg, 1)
call dgemv('N', ni, ni, 1.0d0, Ginv, ni, bg, 1, 0.0d0, gq, 1)
do y = 1, ndof
   a = (y-1)/3 + 1 ; r = mod(y-1,3) + 1
   c2 = coord3 ; c2(r,a) = coord3(r,a) + h
   call ic_bmatrix(ic, c2, Bp)
   c2(r,a) = coord3(r,a) - h
   call ic_bmatrix(ic, c2, Bm)
   call dgemv('T', ni, ndof, 1.0d0/(2.0d0*h), Bp - Bm, ni, gq, 1, 0.0d0, K(:,y), 1)
enddo
K = 0.5d0*(K + transpose(K))
end subroutine ric_grad_curv

subroutine cart_hessian_to_internal(ic, coord3, nat, ni, Hx, Hq, gx)
type(ic_set_t), intent(in) :: ic
integer, intent(in)  :: nat, ni
real(8), intent(in)  :: coord3(3,nat), Hx(3*nat,3*nat)
real(8), intent(out) :: Hq(ni,ni)
real(8), intent(in), optional :: gx(3*nat)
real(8), allocatable :: B(:,:), Ginv(:,:), A(:,:), T(:,:), Hw(:,:)
integer :: ndof
ndof = 3*nat
allocate(B(ni,ndof), Ginv(ni,ni), A(ni,ndof), T(ni,ndof), Hw(ndof,ndof))
Hw = Hx
if (present(gx)) then
   block
      real(8), allocatable :: K(:,:)
      allocate(K(ndof,ndof))
      call ric_grad_curv(ic, coord3, nat, ni, gx, K)
      Hw = Hw - K
   end block
endif
call ic_bmatrix(ic, coord3, B)
call ic_ginv(B, ni, ndof, Ginv)
call dgemm('N','N', ni, ndof, ni, 1.0d0, Ginv, ni, B, ni, 0.0d0, A, ni)
call dgemm('N','N', ni, ndof, ndof, 1.0d0, A, ni, Hw, ndof, 0.0d0, T, ni)
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

end module mod_opt_hessian
