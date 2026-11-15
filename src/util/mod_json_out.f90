! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! <base>.json: the machine-readable twin of the .out file. The .out text is

module mod_json_out
use mod_engine_report, only: engine_report_t
implicit none
private
public :: json_open, json_active, json_set_meta, json_set_geometry
public :: json_result, json_opt_cycle, json_opt_status, json_irc_start, json_irc_point, json_irc_end
public :: json_result_energy, json_freq, json_thermo, json_flush, json_properties, json_polar, json_polar_dynamic, json_dipole_deriv
public :: json_tddft, json_write_report, json_nmr, json_jcoupling, json_nics, json_scan

character(len=:), allocatable, save :: path, s_meta, s_result, s_opt_cycles, s_opt_head, s_irc_head, &
                                       s_irc_points, s_freq, s_thermo, s_props, s_polar, s_dipder, s_tddft, s_tddft_t, s_nmr, s_scan
character(len=:), allocatable, save :: geom_json
character(len=:), allocatable, save :: runtype_name
logical, save :: opened = .false.
integer, save :: n_opt_cycles = 0, n_irc_points = 0
character(len=3), allocatable, save :: sym(:)

logical :: jt_has_ct = .false.
real(8) :: jt_ct(5) = 0.0d0

contains

logical function json_active()
json_active = opened
end function json_active

subroutine json_open(filename)
character(len=*), intent(in) :: filename
path = trim(filename)
s_meta = ''; s_result = ''; s_opt_cycles = ''; s_opt_head = ''; s_irc_head = ''
s_scan = ''
s_irc_points = ''; s_freq = ''; s_thermo = ''; s_props = ''; s_polar = ''; s_dipder = ''; s_tddft = ''; s_tddft_t = ''; s_nmr = ''
geom_json = ''; runtype_name = ''
n_opt_cycles = 0; n_irc_points = 0
opened = .true.
end subroutine json_open

function jr(x) result(s)
real(8), intent(in) :: x
character(len=:), allocatable :: s
character(len=32) :: b
if (x /= x .or. abs(x) > huge(x)) then
   s = 'null'
else
   write(b, '(ES24.16E3)') x
   s = trim(adjustl(b))
endif
end function jr

function ji(i) result(s)
integer, intent(in) :: i
character(len=:), allocatable :: s
character(len=16) :: b
write(b, '(I0)') i
s = trim(b)
end function ji

function js(str) result(s)
character(len=*), intent(in) :: str
character(len=:), allocatable :: s
integer :: i
s = '"'
do i = 1, len_trim(str)
   select case (str(i:i))
   case ('"');  s = s // '\"'
   case ('\');  s = s // '\\'
   case default; s = s // str(i:i)
   end select
enddo
s = s // '"'
end function js

function jb(l) result(s)
logical, intent(in) :: l
character(len=:), allocatable :: s
s = merge('true ', 'false', l)
s = trim(s)
end function jb

function jvec(x) result(s)
real(8), intent(in) :: x(:)
character(len=:), allocatable :: s
integer :: i
s = '['
do i = 1, size(x)
   if (i > 1) s = s // ','
   s = s // jr(x(i))
enddo
s = s // ']'
end function jvec

function jxyz(x) result(s)
real(8), intent(in) :: x(:,:)
character(len=:), allocatable :: s
integer :: i
s = '['
do i = 1, size(x, 2)
   if (i > 1) s = s // ','
   s = s // jvec(x(:, i))
enddo
s = s // ']'
end function jxyz

function jgeom(xyz) result(s)
real(8), intent(in) :: xyz(:,:)
character(len=:), allocatable :: s
s = jxyz(xyz)
end function jgeom

subroutine json_set_meta(version, build, functional, basis, charge, mult, runtype, symbols, znum, xyz)
character(len=*), intent(in) :: version, build, functional, basis, runtype
integer, intent(in) :: charge, mult, znum(:)
character(len=*), intent(in) :: symbols(:)
real(8), intent(in) :: xyz(:,:)
integer :: i
if (.not. opened) return
runtype_name = trim(runtype)
s_meta = '{"program":"Direwolf","version":'//js(version)//',"build":'//js(build)// &
         ',"functional":'//js(functional)//',"basis":'//js(basis)// &
         ',"charge":'//ji(charge)//',"multiplicity":'//ji(mult)//',"runtype":'//js(runtype)// &
         ',"natoms":'//ji(size(znum))//',"units":{"energy":"hartree","gradient":"hartree/bohr",'// &
         '"geometry":"angstrom","frequency":"cm-1","ir_intensity":"km/mol"},"elements":['
do i = 1, size(znum)
   if (i > 1) s_meta = s_meta // ','
   s_meta = s_meta // js(symbols(i))
enddo
s_meta = s_meta // '],"atomic_numbers":['
do i = 1, size(znum)
   if (i > 1) s_meta = s_meta // ','
   s_meta = s_meta // ji(znum(i))
enddo
s_meta = s_meta // '],"input_geometry":' // jgeom(xyz) // '}'
end subroutine json_set_meta

subroutine json_set_geometry(xyz)
real(8), intent(in) :: xyz(:,:)
if (.not. opened) return
geom_json = jgeom(xyz)
end subroutine json_set_geometry

subroutine json_result(converged, energy, hcore, ecoul, exx, exc, enuc, vv10_used, vv10, want_forces, forces)
logical, intent(in) :: converged, vv10_used, want_forces
real(8), intent(in) :: energy, hcore, ecoul, exx, exc, enuc, vv10
real(8), intent(in) :: forces(:,:)
if (.not. opened) return
s_result = '{"converged":'//jb(converged)//',"energy":'//jr(energy)// &
           ',"energy_components":{"hcore":'//jr(hcore)//',"coulomb":'//jr(ecoul)// &
           ',"exact_exchange":'//jr(exx)//',"xc":'//jr(exc)//',"nuclear_repulsion":'//jr(enuc)
if (vv10_used) s_result = s_result // ',"vv10_nlc":'//jr(vv10)
s_result = s_result // '}'
if (want_forces) s_result = s_result // ',"forces":'//jxyz(transpose(forces))
if (len(geom_json) > 0) s_result = s_result // ',"geometry":'//geom_json
s_result = s_result // '}'
call json_flush()
end subroutine json_result

subroutine json_result_energy(energy)
real(8), intent(in) :: energy
if (.not. opened) return
if (len(s_result) > 0) return
s_result = '{"converged":true,"energy":'//jr(energy)
if (len(geom_json) > 0) s_result = s_result // ',"geometry":'//geom_json
s_result = s_result // '}'
call json_flush()
end subroutine json_result_energy

subroutine json_opt_cycle(cyc, energy, gmax, grms, have_step, smax, srms, trust, xyz)
integer, intent(in) :: cyc
real(8), intent(in) :: energy, gmax, grms, smax, srms, trust, xyz(:,:)
logical, intent(in) :: have_step
if (.not. opened) return
if (n_opt_cycles > 0) s_opt_cycles = s_opt_cycles // ','//new_line('a')
s_opt_cycles = s_opt_cycles // '{"cycle":'//ji(cyc)//',"energy":'//jr(energy)// &
               ',"gmax":'//jr(gmax)//',"grms":'//jr(grms)
if (have_step) s_opt_cycles = s_opt_cycles // ',"smax":'//jr(smax)//',"srms":'//jr(srms)
s_opt_cycles = s_opt_cycles // ',"trust":'//jr(trust)//',"geometry":'//jgeom(xyz)//'}'
n_opt_cycles = n_opt_cycles + 1
call json_flush()
end subroutine json_opt_cycle

subroutine json_opt_status(converged, ncycles, message)
logical, intent(in) :: converged
integer, intent(in) :: ncycles
character(len=*), intent(in) :: message
if (.not. opened) return
s_opt_head = '"converged":'//jb(converged)//',"ncycles":'//ji(ncycles)//',"message":'//js(message)
call json_flush()
end subroutine json_opt_status

subroutine json_scan(coords, nd, npt, target, value, energy, rel_kcal, converged, ncycles, xyz)
character(len=*), intent(in) :: coords(:)
integer, intent(in) :: nd, npt, ncycles(:)
real(8), intent(in) :: target(:,:), value(:,:), energy(:), rel_kcal(:), xyz(:,:,:)
logical, intent(in) :: converged(:)
integer :: i, j
if (.not. opened) return
if (nd == 1) then
   s_scan = '{"coordinate":'//js(coords(1))//',"points":['//new_line('a')
else
   s_scan = '{"coordinates":['
   do j = 1, nd
      if (j > 1) s_scan = s_scan // ','
      s_scan = s_scan // js(coords(j))
   enddo
   s_scan = s_scan // '],"points":['//new_line('a')
endif
do i = 1, npt
   if (i > 1) s_scan = s_scan // ','//new_line('a')
   if (nd == 1) then
      s_scan = s_scan // '{"target":'//jr(target(1,i))//',"value":'//jr(value(1,i))
   else
      s_scan = s_scan // '{"target":'//jvec(target(:,i))//',"value":'//jvec(value(:,i))
   endif
   s_scan = s_scan // ',"energy":'//jr(energy(i))// &
            ',"relative_energy_kcal_mol":'//jr(rel_kcal(i))//',"converged":'//jb(converged(i))// &
            ',"ncycles":'//ji(ncycles(i))//',"geometry":'//jgeom(xyz(:,:,i))//'}'
enddo
s_scan = s_scan // ']}'
s_opt_cycles = '' ; s_opt_head = '' ; n_opt_cycles = 0
call json_flush()
end subroutine json_scan

subroutine json_irc_start(ts_energy, lowest_freq, nimag)
real(8), intent(in) :: ts_energy, lowest_freq
integer, intent(in) :: nimag
if (.not. opened) return
s_irc_head = '"ts_energy":'//jr(ts_energy)//',"ts_lowest_frequency":'//jr(lowest_freq)// &
             ',"ts_imaginary_modes":'//ji(nimag)
end subroutine json_irc_start

subroutine json_irc_point(s, energy, rel_kcal, gmax, xyz)
real(8), intent(in) :: s, energy, rel_kcal, gmax, xyz(:,:)
if (.not. opened) return
if (n_irc_points > 0) s_irc_points = s_irc_points // ','//new_line('a')
s_irc_points = s_irc_points // '{"s":'//jr(s)//',"energy":'//jr(energy)// &
               ',"relative_energy_kcal_mol":'//jr(rel_kcal)//',"gmax":'//jr(gmax)//',"geometry":'//jgeom(xyz)//'}'
n_irc_points = n_irc_points + 1
end subroutine json_irc_point

subroutine json_irc_end(fwd_reason, bwd_reason, n_eg, n_hess)
character(len=*), intent(in) :: fwd_reason, bwd_reason
integer, intent(in) :: n_eg, n_hess
if (.not. opened) return
s_irc_head = s_irc_head // ',"forward_stop":'//js(fwd_reason)//',"backward_stop":'//js(bwd_reason)// &
             ',"energy_gradient_calls":'//ji(n_eg)//',"hessian_calls":'//ji(n_hess)
call json_flush()
end subroutine json_irc_end

subroutine json_freq(freq, inten, redmass, fc, disp)
real(8), intent(in) :: freq(:), inten(:), redmass(:), fc(:)
real(8), intent(in) :: disp(:,:,:)
integer :: k
if (.not. opened) return
s_freq = '['
do k = 1, size(freq)
   if (k > 1) s_freq = s_freq // ','//new_line('a')
   s_freq = s_freq // '{"mode":'//ji(k)//',"frequency":'//jr(freq(k))//',"ir_intensity":'//jr(inten(k))// &
            ',"reduced_mass":'//jr(redmass(k))//',"force_constant":'//jr(fc(k))// &
            ',"displacements":'//jxyz(disp(:, :, k))//'}'
enddo
s_freq = s_freq // ']'
call json_flush()
end subroutine json_freq

subroutine json_thermo(T, P, mult, sigma, mass, zpe, e_corr, h_corr, g_corr, s_tot, cv_tot, &
                       s_el, s_trans, s_rot, s_vib, e0)
real(8), intent(in) :: T, P, mass, zpe, e_corr, h_corr, g_corr, s_tot, cv_tot, s_el, s_trans, s_rot, s_vib, e0
integer, intent(in) :: mult, sigma
if (.not. opened) return
s_thermo = '{"temperature":'//jr(T)//',"pressure":'//jr(P)//',"multiplicity":'//ji(mult)// &
           ',"rotational_symmetry_number":'//ji(sigma)//',"mass_amu":'//jr(mass)// &
           ',"zero_point_correction":'//jr(zpe)//',"thermal_correction_E":'//jr(e_corr)// &
           ',"thermal_correction_H":'//jr(h_corr)//',"thermal_correction_G":'//jr(g_corr)// &
           ',"entropy_total":'//jr(s_tot)//',"entropy_electronic":'//jr(s_el)//',"entropy_translational":'//jr(s_trans)// &
           ',"entropy_rotational":'//jr(s_rot)//',"entropy_vibrational":'//jr(s_vib)//',"cv_total":'//jr(cv_tot)
if (e0 /= 0d0) s_thermo = s_thermo // ',"electronic_energy":'//jr(e0)// &
           ',"sum_electronic_zero_point":'//jr(e0 + zpe)//',"sum_electronic_thermal_E":'//jr(e0 + e_corr)// &
           ',"sum_electronic_thermal_H":'//jr(e0 + h_corr)//',"sum_electronic_thermal_G":'//jr(e0 + g_corr)
s_thermo = s_thermo // '}'
call json_flush()
end subroutine json_thermo

subroutine json_properties(dip_au, q_mulliken, q_lowdin, have_resp, q_resp, restricted, &
                           e_alpha, n_alpha, e_beta, n_beta)
real(8), intent(in) :: dip_au(3), q_mulliken(:), q_lowdin(:), q_resp(:), e_alpha(:), e_beta(:)
logical, intent(in) :: have_resp, restricted
integer, intent(in) :: n_alpha, n_beta
real(8), parameter :: AU2DEBYE = 2.541746473d0
if (.not. opened) return
s_props = '{"dipole_au":'//jvec(dip_au)//',"dipole_debye":'//jvec(dip_au*AU2DEBYE)// &
          ',"dipole_norm_debye":'//jr(sqrt(sum(dip_au**2))*AU2DEBYE)// &
          ',"charges":{"mulliken":'//jvec(q_mulliken)//',"lowdin":'//jvec(q_lowdin)
if (have_resp) s_props = s_props // ',"resp":'//jvec(q_resp)
s_props = s_props // '},"orbital_energies":{"restricted":'//jb(restricted)// &
          ',"n_occupied_alpha":'//ji(n_alpha)//',"alpha":'//jvec(e_alpha)
if (.not. restricted) s_props = s_props // ',"n_occupied_beta":'//ji(n_beta)//',"beta":'//jvec(e_beta)
if (n_alpha >= 1 .and. n_alpha < size(e_alpha)) s_props = s_props // &
          ',"homo_alpha":'//jr(e_alpha(n_alpha))//',"lumo_alpha":'//jr(e_alpha(n_alpha+1))
s_props = s_props // '}}'
call json_flush()
end subroutine json_properties

subroutine json_polar(alpha, aiso, aani, niter)
real(8), intent(in) :: alpha(3,3), aiso, aani
integer, intent(in) :: niter
real(8), parameter :: AU2ANG3 = 0.148184711d0
if (.not. opened) return
s_polar = '{"tensor_au":'//jxyz(alpha)//',"isotropic_au":'//jr(aiso)//',"anisotropy_au":'//jr(aani)// &
          ',"isotropic_angstrom3":'//jr(aiso*AU2ANG3)//',"anisotropy_angstrom3":'//jr(aani*AU2ANG3)// &
          ',"cphf_iterations":'//ji(niter)//'}'
call json_flush()
end subroutine json_polar

subroutine haeberlen(t, ani)
real(8), intent(in) :: t(3,3)
real(8), intent(out) :: ani
real(8) :: s(3,3), w(3), work(20), tmp
integer :: info, i, j
s = 0.5d0*(t + transpose(t))
call dsyev('N', 'U', 3, s, 3, w, work, 20, info)
do i = 1,2
   do j = i+1,3
      if (abs(w(j)) < abs(w(i))) then
         tmp = w(i) ; w(i) = w(j) ; w(j) = tmp
      endif
   enddo
enddo
ani = w(3) - 0.5d0*(w(1) + w(2))
end subroutine haeberlen

subroutine json_nmr(nat, sym, sdia, spara, niter)
integer, intent(in) :: nat, niter
character(len=*), intent(in) :: sym(nat)
real(8), intent(in) :: sdia(3,3,nat), spara(3,3,nat)
real(8) :: tot(3,3), iso, ani
integer :: ia
character(len=:), allocatable :: rec
if (.not. opened) return
rec = '{"unit":"ppm","cphf_iterations":'//ji(niter)//',"atoms":['
do ia = 1, nat
   tot = sdia(:,:,ia) + spara(:,:,ia)
   iso = (tot(1,1) + tot(2,2) + tot(3,3))/3.0d0
   call haeberlen(tot, ani)
   if (ia > 1) rec = rec // ','
   rec = rec // '{"index":'//ji(ia)//',"element":"'//trim(sym(ia))//'","isotropic":'//jr(iso)//',"anisotropy":'//jr(ani)// &
         ',"diamagnetic_isotropic":'//jr((sdia(1,1,ia)+sdia(2,2,ia)+sdia(3,3,ia))/3.0d0)// &
         ',"paramagnetic_isotropic":'//jr((spara(1,1,ia)+spara(2,2,ia)+spara(3,3,ia))/3.0d0)// &
         ',"tensor":'//jxyz(transpose(tot))//',"diamagnetic":'//jxyz(transpose(sdia(:,:,ia)))// &
         ',"paramagnetic":'//jxyz(transpose(spara(:,:,ia)))//'}'
enddo
s_nmr = rec // ']}'
call json_flush()
end subroutine json_nmr

subroutine json_nics(np, label, xyz_ang, niso, nzz)
integer, intent(in) :: np
character(len=*), intent(in) :: label(np)
real(8), intent(in) :: xyz_ang(3,np), niso(np), nzz(np)
integer :: ip, l
character(len=:), allocatable :: rec
if (.not. opened) return
if (len(s_nmr) == 0) return
rec = ',"nics":{"unit":"ppm","probes":['
do ip = 1, np
   if (ip > 1) rec = rec // ','
   rec = rec // '{"label":"'//trim(label(ip))//'","xyz_angstrom":['//jr(xyz_ang(1,ip))//','//jr(xyz_ang(2,ip))//','// &
         jr(xyz_ang(3,ip))//'],"nics":'//jr(niso(ip))//',"nics_zz":'//jr(nzz(ip))//'}'
enddo
rec = rec // ']}'
l = len(s_nmr)
s_nmr = s_nmr(1:l-1) // rec // '}'
call json_flush()
end subroutine json_nics

subroutine json_jcoupling(nat, J, Kred, Jfc, Jsd, Jpso, Jdso, terms)
integer, intent(in) :: nat
real(8), intent(in) :: J(nat,nat), Kred(nat,nat), Jfc(nat,nat), Jsd(nat,nat), Jpso(nat,nat), Jdso(nat,nat)
character(len=*), intent(in) :: terms
integer :: ia, ib, l
character(len=:), allocatable :: rec
logical :: first
if (.not. opened) return
if (len(s_nmr) == 0) return
rec = ',"jcoupling":{"terms":"'//trim(terms)//'","unit":"Hz","pairs":['
first = .true.
do ia = 1, nat
   do ib = ia+1, nat
      if (J(ia,ib) == 0.0d0) cycle
      if (.not. first) rec = rec // ','
      first = .false.
      rec = rec // '{"a":'//ji(ia)//',"b":'//ji(ib)//',"J":'//jr(J(ia,ib))//',"K_reduced_1e19_T2_per_J":'//jr(Kred(ia,ib)) &
            //',"FC":'//jr(Jfc(ia,ib))//',"SD":'//jr(Jsd(ia,ib))//',"PSO":'//jr(Jpso(ia,ib))//',"DSO":'//jr(Jdso(ia,ib))//'}'
   enddo
enddo
rec = rec // ']}'
l = len(s_nmr)
s_nmr = s_nmr(1:l-1) // rec // '}'
call json_flush()
end subroutine json_jcoupling

subroutine json_polar_dynamic(omega, alpha, aiso, aani, niter)
real(8), intent(in) :: omega, alpha(3,3), aiso, aani
integer, intent(in) :: niter
real(8), parameter :: AU2ANG3 = 0.148184711d0, HARTREE_EV = 27.211386245988d0, HARTREE_NM = 45.56335252760d0
character(len=:), allocatable :: rec
integer :: l
if (.not. opened) return
if (len(s_polar) == 0) return
rec = '{"omega_au":'//jr(omega)//',"omega_ev":'//jr(omega*HARTREE_EV)
if (omega > 0.0d0) rec = rec // ',"wavelength_nm":'//jr(HARTREE_NM/omega)
rec = rec // ',"tensor_au":'//jxyz(alpha)//',"isotropic_au":'//jr(aiso)//',"anisotropy_au":'//jr(aani)// &
      ',"isotropic_angstrom3":'//jr(aiso*AU2ANG3)//',"anisotropy_angstrom3":'//jr(aani*AU2ANG3)// &
      ',"cphf_iterations":'//ji(niter)//'}'
l = len(s_polar)
if (index(s_polar, '"dynamic"') == 0) then
   s_polar = s_polar(1:l-1) // ',"dynamic":[' // rec // ']}'
else
   s_polar = s_polar(1:l-2) // ',' // rec // ']}'
endif
call json_flush()
end subroutine json_polar_dynamic

subroutine json_tddft(istate, w, f, mu, n, ifrom, ito, wt, napply, triplet, spins, method &
   )
integer, intent(in) :: istate, n, ifrom(*), ito(*), napply
character(len=*), intent(in), optional :: method
logical, intent(in), optional :: triplet
character(len=1), intent(in), optional :: spins(*)
logical :: trip
real(8), intent(in) :: w, f, mu(3), wt(*)
real(8), parameter :: HARTREE_EV = 27.211386245988d0, HARTREE_NM = 45.56335252760d0
character(len=:), allocatable :: rec, meth
integer :: k
if (.not. opened) return
meth = ''
if (present(method)) meth = method
rec = '{"state":'//ji(istate)//',"energy_au":'//jr(w)//',"energy_ev":'//jr(w*HARTREE_EV)// &
      ',"wavelength_nm":'//jr(HARTREE_NM/w)//',"oscillator_strength":'//jr(f)// &
      ',"transition_dipole_au":['//jr(mu(1))//','//jr(mu(2))//','//jr(mu(3))//']'
if (jt_has_ct) rec = rec // ',"hole_electron":{"d_ct_angstrom":'//jr(jt_ct(1))//',"q_ct":'//jr(jt_ct(2))// &
      ',"sigma_hole_angstrom":'//jr(jt_ct(3))//',"sigma_electron_angstrom":'//jr(jt_ct(4))//',"lambda":'//jr(jt_ct(5))//'}'
jt_has_ct = .false.
rec = rec // ',"transitions":['
do k = 1, n
   if (k > 1) rec = rec // ','
   if (present(spins)) then
      rec = rec // '{"from":'//ji(ifrom(k))//',"to":'//ji(ito(k))//',"spin":"'//spins(k)//'","weight":'//jr(wt(k))//'}'
   else
      rec = rec // '{"from":'//ji(ifrom(k))//',"to":'//ji(ito(k))//',"weight":'//jr(wt(k))//'}'
   endif
enddo
rec = rec // ']}'
trip = .false.
if (present(triplet)) trip = triplet
if (trip) then
   if (istate == 1 .or. len(s_tddft_t) == 0) then
      s_tddft_t = '{"singlet":false,"method":"'//meth//'","response_vectors":'//ji(napply)//',"states":['//rec//']}'
   else
      s_tddft_t = s_tddft_t(1:len(s_tddft_t)-2) // ',' // rec // ']}'
   endif
else if (istate == 1 .or. len(s_tddft) == 0) then
   if (present(spins)) then
      s_tddft = '{"unrestricted":true,"method":"'//meth//'","response_vectors":'//ji(napply)//',"states":['//rec//']}'
   else
      s_tddft = '{"singlet":true,"method":"'//meth//'","response_vectors":'//ji(napply)//',"states":['//rec//']}'
   endif
else
   s_tddft = s_tddft(1:len(s_tddft)-2) // ',' // rec // ']}'
endif
call json_flush()
end subroutine json_tddft

subroutine json_dipole_deriv(dmu, natoms)
integer, intent(in) :: natoms
real(8), intent(in) :: dmu(3, 3*natoms)
integer :: a
if (.not. opened) return
s_dipder = '['
do a = 1, natoms
   if (a > 1) s_dipder = s_dipder // ','
   s_dipder = s_dipder // jxyz(dmu(:, 3*(a-1)+1:3*(a-1)+3))
enddo
s_dipder = s_dipder // ']'
call json_flush()
end subroutine json_dipole_deriv

subroutine write_state(k, s, o1, o2, spins)
use mod_engine_report, only: tddft_set_t
integer, intent(in) :: k, o1, o2
type(tddft_set_t), intent(in) :: s
logical, intent(in) :: spins
jt_has_ct = s%has_ct
if (s%has_ct) jt_ct = s%ct(:,k)
   if (spins) then
      call json_tddft(k, s%w(k), s%f(k), s%mu(:,k), s%ntr(k), s%ifrom(o1:o2), s%ito(o1:o2), s%wt(o1:o2), s%napply, &
                      spins=s%spin(o1:o2), method=trim(s%method))
   else
      call json_tddft(k, s%w(k), s%f(k), s%mu(:,k), s%ntr(k), s%ifrom(o1:o2), s%ito(o1:o2), s%wt(o1:o2), s%napply, &
                      method=trim(s%method))
   endif
end subroutine write_state

subroutine json_write_report(rep)
type(engine_report_t), intent(in) :: rep
integer :: k, ks, o1, o2
if (.not. opened) return
if (rep%energy_only) call json_result_energy(rep%energy)
if (rep%has_props) call json_properties(rep%dip_au, rep%q_mulliken, rep%q_lowdin, rep%have_resp, rep%q_resp, &
                                        rep%restricted, rep%e_alpha, rep%n_alpha, rep%e_beta, rep%n_beta)
if (rep%has_polar) then
   call json_polar(rep%stat%alpha, rep%stat%aiso, rep%stat%aani, rep%stat%niter)
   do k = 1, rep%n_dyn
      call json_polar_dynamic(rep%dyn(k)%omega, rep%dyn(k)%alpha, rep%dyn(k)%aiso, rep%dyn(k)%aani, rep%dyn(k)%niter)
   enddo
endif
do ks = 1, 3
   if (.not. rep%tddft(ks)%used) cycle
   do k = 1, rep%tddft(ks)%nstate
      o1 = rep%tddft(ks)%off(k) + 1
      o2 = rep%tddft(ks)%off(k) + rep%tddft(ks)%ntr(k)
      associate (s => rep%tddft(ks))
      select case (ks)
      case (2)
         call json_tddft(k, s%w(k), s%f(k), s%mu(:,k), s%ntr(k), s%ifrom(o1:o2), s%ito(o1:o2), s%wt(o1:o2), &
                         s%napply, triplet=.true., method=trim(s%method))
      case (3)
         call write_state(k, s, o1, o2, spins=.true.)
      case default
         call write_state(k, s, o1, o2, spins=.false.)
      end select
      end associate
   enddo
enddo
if (rep%has_freq) call json_freq(rep%freq, rep%inten, rep%redmass, rep%fc, rep%disp)
if (rep%has_thermo) call json_thermo(rep%th_T, rep%th_P, rep%th_mult, rep%th_sigma, rep%th_mass, rep%th_zpe, &
                                     rep%th_e, rep%th_h, rep%th_g, rep%th_s, rep%th_cv, rep%th_s_el, &
                                     rep%th_s_trans, rep%th_s_rot, rep%th_s_vib, rep%th_e0)
if (rep%has_dipder) call json_dipole_deriv(rep%dmu, rep%dipder_natoms)
end subroutine json_write_report

subroutine json_flush()
integer :: u, ios
logical :: first
if (.not. opened) return
open(newunit=u, file=path, status='replace', action='write', iostat=ios)
if (ios /= 0) return
write(u, '(A)') '{'
first = .true.
call item(u, first, '"meta"', s_meta)
call item(u, first, '"result"', s_result)
if (len(s_opt_cycles) > 0 .or. len(s_opt_head) > 0) &
   call item(u, first, '"opt"', '{'//s_opt_head//merge(',', ' ', len(s_opt_head) > 0)// &
             '"cycles":['//new_line('a')//s_opt_cycles//']}')
if (len(s_irc_head) > 0) &
   call item(u, first, '"irc"', '{'//s_irc_head//',"points":['//new_line('a')//s_irc_points//']}')
call item(u, first, '"scan"', s_scan)
call item(u, first, '"properties"', s_props)
call item(u, first, '"polarizability"', s_polar)
call item(u, first, '"nmr"', s_nmr)
call item(u, first, '"tddft"', s_tddft)
call item(u, first, '"tddft_triplet"', s_tddft_t)
call item(u, first, '"freq"', s_freq)
call item(u, first, '"dipole_derivatives"', s_dipder)
call item(u, first, '"thermo"', s_thermo)
write(u, '(A)') '}'
close(u)
end subroutine json_flush

subroutine item(u, first, key, val)
integer, intent(in) :: u
logical, intent(inout) :: first
character(len=*), intent(in) :: key, val
if (len(val) == 0) return
if (.not. first) write(u, '(A)') ','
first = .false.
write(u, '(A)', advance='no') key//': '//val
end subroutine item

end module mod_json_out
