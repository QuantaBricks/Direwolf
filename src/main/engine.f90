! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Top-level driver: EngineUp entry point wiring molecule/basis setup through SCF, force, and property evaluation.


subroutine EngineUp(ncenters,imult,icharge,functional_in,&
                    coord,atomchg,baselable,ecplabel,atombasis,atomecp,&
                    j_mode,k_mode,ri_aux_basis,puream,harris_guess,do_force,vv10_nonself,mem_cap_gb,&
                    estimate_only,n_threads,&
                    npc,pc_charge,pc_coord,&
                    molden_write,molden_file,molden_read,molden_read_file,&
                    cosmo_on,cosmo_epsilon,cosmo_radii_scale,cosmo_avg_area,cosmo_sigma_rav,&
                    cosmo_cavity_type,cosmo_rsolv,cosmo_ks_nseg,cosmo_ks_nface,cosmo_sigma_profile_file,&
                    cosmo_smd,cosmo_solvent,&
                    force_out,energy_out,MLcharge_out,iconv,econv,&
                    mem_grid_gb,mem_2e_gb,scf_conv_level,basedir,&
                    dens_in_a,dens_in_b,dens_out_a,dens_out_b,&
                    opts)

use MOL_info
use GRID_info, only: Grids
use GRID_info, only: xcgrid_level, &
                     xcgrid_theory_dynamic, xcgrid_theory_refined, &
                     force_dense_allow, &
                     XCGRID_RSCALE, XCGRID_SPH_IN, XCGRID_SPH_EDGE, &
                     xcgrid_switch_prms, XCGRID_COARSE, XCGRID_FINE, &
                     XCGRID_L4, XCGRID_L5, XCGRID_L6, XCGRID_L7, &
                     xc_direct_mode, XC_DIRECT_NCONTS, force_dense, force_dense_mgga, &
                     xcgrid_tune_from_env, xcgrid_use_theory_table, xcgrid_nmr_active, &
                     xcgrid_audit_only, xcgrid_int_rho, &
                     xcgrid_radpot_rowscale, xcgrid_prune_orders, &
                     xcgrid_prune_dmin_bounds, xcgrid_prune_dmin_orders, &
                     xcgrid_lebedev_supported, &
                     grid_cache_mode, max_batch_nsig, grid_cache_capacity, ngrids
use mod_integrals, only: resp_charges, mayer_bo, calc_mayer_bond_orders, esp_at_grid_batch, &
                          integrals_init, engine_use_df, engine_use_df_j, engine_use_df_k, &
                          engine_df_aux_basis, &
                          engine_estimate_only, integrals_estimate_aux_size, &
                          nRec, nContsAux, decide_df_force_mode, engine_puream, &
                          engine_harris_guess, calc_properties, engine_do_force, &
                          hessian_run
use mod_engine_options, only: engine_options_t
use mod_engine_report, only: rpt, report_reset
use mod_hessian_freq, only: freq_unit, freq_mult, freq_thermo, freq_base, freq_temp, freq_pres
use mod_meminfo, only: engine_mem_cap_bytes, engine_avail_at_start_bytes, get_memory_budget_bytes, get_process_memory_bytes, &
    release_free_memory
use mod_xc, only: xc_select_functional, xc_uses_tau, xc_uses_vv10
use mod_vv10, only: engine_vv10_nonself, vv10_dynamic, vv10_active_now, &
                     vv10_just_activated, vv10_switch_prms, vv10_flush_on_switch, &
                     VV10_DYNGRID_NCONTS
use mod_exchange, only: HF_exchange_frac, RS_omega, RS_beta, cosx_enabled
use mod_dispersion, only: dispersion_d2, dispersion_chg
use mod_dispersion_d3, only: dispersion_d3
use mod_dispersion_d4, only: dispersion_d4
use mod_gcp, only: gcp_correction
use mod_profile, only: prof_reset, prof_start, prof_stop, prof_report, fmt_gb
use mod_output_format, only: finalize_output_file, print_distance_matrix
use mod_molden, only: write_molden, read_molden, elem_symbol
use mod_cube, only: write_cube_files
use mod_cosmo, only: cosmo_enabled, cosmo_init, cosmo_report_sigma_profile, cosmo_finalize, &
                      cosmo_solvent_epsilon, cosmo_solvent_smd_params, cosmo_set_sigma_profile_debug, &
                      cosmo_write_dot_cosmo_file
use mod_smd_cds, only: smd_cds_energy_force
use omp_lib, only: omp_set_num_threads, omp_get_max_threads, omp_get_num_procs
implicit none
include "parameter.h"
interface
   subroutine gridgen_count_points(totGrid, grdRec)
   integer,intent(out) :: totGrid
   integer,allocatable,intent(out) :: grdRec(:,:)
   end subroutine gridgen_count_points
   subroutine openblas_set_num_threads(num_threads) bind(c, name="openblas_set_num_threads")
   use iso_c_binding, only: c_int
   integer(c_int),value :: num_threads
   end subroutine openblas_set_num_threads
end interface
integer,external :: ecp_core_electrons_prescan
character(len=30),external :: resolve_basis_name
character(len=30),external :: auto_ri_aux_basis
logical,external :: aux_basis_covers_all_atoms
integer       :: ncenters,imult,icharge,atomchg(ncenters)
character     :: baselable*30,functional_in*30
character  :: ecplabel*30
character  :: atombasis(ncenters)*30, atomecp(ncenters)*30
real(8)       :: coord(ncenters,3)
character(len=*),intent(in) :: j_mode, k_mode
logical,intent(in) :: puream
logical,intent(in) :: harris_guess
character(len=*),intent(in) :: ri_aux_basis
logical,intent(in) :: do_force
logical,intent(in) :: vv10_nonself
real(8),intent(in) :: mem_cap_gb
logical,intent(in) :: estimate_only
integer,intent(in) :: n_threads
integer,intent(in) :: npc
real(8),intent(in) :: pc_charge(npc),pc_coord(npc,3)
logical,intent(in) :: molden_write, molden_read
character(len=*),intent(in) :: molden_file, molden_read_file
logical,intent(in) :: cosmo_on
real(8),intent(in) :: cosmo_epsilon, cosmo_radii_scale, cosmo_sigma_rav
real(8),intent(in) :: cosmo_avg_area
character(*),intent(in) :: cosmo_sigma_profile_file
character(len=*),intent(in) :: cosmo_cavity_type
real(8),intent(in) :: cosmo_rsolv
integer,intent(in) :: cosmo_ks_nseg, cosmo_ks_nface
logical,intent(in) :: cosmo_smd
character(len=*),intent(in) :: cosmo_solvent
integer       :: iconv
real(8)       :: energy_out,econv
real(8)       :: force_out(ncenters,3),MLcharge_out(ncenters)
real(8),intent(out) :: mem_grid_gb
real(8),intent(out) :: mem_2e_gb
integer,intent(in) :: scf_conv_level
character(len=*),intent(in) :: basedir
real(8),intent(in) :: dens_in_a(:,:), dens_in_b(:,:)
real(8),intent(out),allocatable :: dens_out_a(:,:), dens_out_b(:,:)
type(engine_options_t),intent(inout) :: opts
real(8)    :: Emax,Erms,Pmax,Prms
integer    :: itmax

integer       :: Gtype
integer :: i, j
integer :: info
real(8),allocatable :: Pa_chk(:,:),Pb_chk(:,:)
real(8) :: disp_s6
logical :: use_chg_dispersion
character(8) :: disp_d3_damping
character(20) :: disp_d4_method
character(len=20) :: gcp_method
character(len=30) :: functional_grid
logical :: chk_ok
logical :: has_warm_start

real(8)   :: dist

real(8) :: smd_n_val, smd_alpha_val, smd_beta_val, smd_gamma_val, smd_phi_val, smd_psi_val
logical :: smd_is_water

if (n_threads .gt. 0) call omp_set_num_threads(n_threads)

call openblas_set_num_threads(1)

engine_verbose = opts%verbose ; engine_quiet = opts%quiet
force_dense_allow = opts%force_dense_allow
resp_charges_on = opts%resp_charges
block
   use mod_integrals, only: cosx_force_loose
   use mod_cosmo, only: cosmo_set_cosmors_request
   cosx_force_loose = opts%cosx_force_loose
   call cosmo_set_cosmors_request(opts%cosmors_on, trim(opts%cosmors_file), trim(opts%cosmors_compound))
end block
block
   use mod_cosmo_state, only: cosmo_cpcm
   use mod_cosmo_radii, only: cavity_h_radius_ang
   cosmo_cpcm = opts%cosmo_cpcm
   cavity_h_radius_ang = merge(1.10d0, 0.0d0, opts%cosmo_cpcm)
end block
if (.not. engine_quiet) print *,"System cores available:",omp_get_num_procs(), &
        " Direwolf OpenMP threads in use:",omp_get_max_threads()
call flush(6)

call reset_engine_state()
opts%uks_restart = .false.
opts%done = .false.
call report_reset(opts%progress)
opts%rep = rpt
freq_unit = opts%freq_unit ; freq_mult = opts%freq_mult ; freq_thermo = opts%freq_thermo
freq_temp = opts%freq_temp ; freq_pres = opts%freq_pres
freq_base = opts%freq_base
call engine_startup_mark("reset_engine_state")
call prof_reset()

Natoms = ncenters
Multi = imult
Charge = icharge
Functional = functional_in
block
   integer :: fic, fii
   do fii = 1, len_trim(Functional)
      fic = ichar(Functional(fii:fii))
      if (fic >= ichar('a') .and. fic <= ichar('z')) &
         Functional(fii:fii) = char(fic - 32)
   enddo
end block
functional_grid = Functional
disp_s6 = 0.0d0
block
   integer :: dpos
   dpos = index(Functional, '-D2')
   if (dpos .eq. 0) dpos = index(Functional, '_D2')
   if (dpos .gt. 0) then
      select case (trim(Functional(1:dpos-1)))
      case ("B3LYP","B3LYP_HYB","B3LYPHYB")
         disp_s6 = 1.05d0
      case ("PBE_PBE","PBE","PBEPBE")
         disp_s6 = 0.75d0
      case default
         disp_s6 = 1.0d0
      end select
      Functional = Functional(1:dpos-1)
   endif
end block
disp_d3_damping = ""
block
   integer :: dpos
   dpos = index(Functional, '-D3BJ')
   if (dpos .eq. 0) dpos = index(Functional, '_D3BJ')
   if (dpos .gt. 0) then
      disp_d3_damping = "BJ"
      Functional = Functional(1:dpos-1)
   else
      dpos = index(Functional, '-D3')
      if (dpos .eq. 0) dpos = index(Functional, '_D3')
      if (dpos .gt. 0) then
         disp_d3_damping = "ZERO"
         Functional = Functional(1:dpos-1)
      endif
   endif
end block
gcp_method = ""
disp_d4_method = ""
block
   character(len=30) :: req_3c_basis
   req_3c_basis = ""
   select case (trim(Functional))
   case ("B97-3C","B97_3C","B973C")
      disp_d3_damping = "BJ"
      gcp_method = "b973c"
      req_3c_basis = "mtzvp"
   case ("R2SCAN-3C","R2SCAN_3C","R2SCAN3C")
      disp_d4_method = "r2scan-3c"
      gcp_method = "r2scan3c"
      req_3c_basis = "mtzvpp"
   case ("WB97X-3C","WB97X_3C","WB97X3C")
      disp_d4_method = "wb97x-3c"
      req_3c_basis = "vdzp"
   end select
   if (len_trim(req_3c_basis) .gt. 0) then
      block
         character(len=30) :: cur_u
         integer :: ic2, ii2
         cur_u = adjustl(baselable)
         do ii2 = 1, len_trim(cur_u)
            ic2 = ichar(cur_u(ii2:ii2))
            if (ic2 >= ichar('A') .and. ic2 <= ichar('Z')) cur_u(ii2:ii2) = char(ic2 + 32)
         enddo
         if (trim(cur_u) .ne. trim(req_3c_basis)) then
            write(*,'(A,A,A,A,A,A)') ' [Notice] ', trim(Functional), &
               ' is a fixed-basis composite method - overriding baselabel "', &
               trim(adjustl(baselable)), '" with its required basis "', trim(req_3c_basis)//'"'
            baselable = req_3c_basis
         endif
      end block
   endif
end block
if (trim(Functional) .eq. "R2SCAN-3C" .or. trim(Functional) .eq. "R2SCAN_3C" .or. &
    trim(Functional) .eq. "R2SCAN3C") &
   Functional = "R2SCAN"
use_chg_dispersion = (trim(Functional) .eq. "WB97X-D" .or. trim(Functional) .eq. "WB97X_D")
call xc_select_functional(trim(Functional), HF_exchange_frac, RS_omega, RS_beta)
call engine_startup_mark("xc_select")

block
   character(len=16) :: j_norm, k_norm
   integer :: ic, ii
   j_norm = adjustl(j_mode)
   k_norm = adjustl(k_mode)
   do ii = 1,len(j_norm)
      ic = iachar(j_norm(ii:ii))
      if (ic .ge. iachar('a') .and. ic .le. iachar('z')) j_norm(ii:ii) = achar(ic - 32)
   enddo
   do ii = 1,len(k_norm)
      ic = iachar(k_norm(ii:ii))
      if (ic .ge. iachar('a') .and. ic .le. iachar('z')) k_norm(ii:ii) = achar(ic - 32)
   enddo
   select case (trim(j_norm))
   case ('EXACT')
      engine_use_df_j = .false.
   case ('RI')
      engine_use_df_j = .true.
   case default
      print *, 'Error: J must be "exact" or "RI", got "', trim(j_mode), '"'
      stop 1
   end select
   engine_do_force = do_force
   select case (trim(k_norm))
   case ('EXACT')
      engine_use_df_k = .false.
      cosx_enabled = .false.
   case ('RI')
      engine_use_df_k = .true.
      cosx_enabled = .false.
   case ('COSX')
      engine_use_df_k = .false.
      cosx_enabled = .true.
   case default
      print *, 'Error: K must be "exact", "RI", or "cosx", got "', trim(k_mode), '"'
      stop 1
   end select
end block
engine_use_df = engine_use_df_j .or. engine_use_df_k
engine_vv10_nonself = vv10_nonself
vv10_dynamic = .false.
vv10_active_now = .true.
vv10_just_activated = .false.
vv10_flush_on_switch = .true.
engine_puream = puream
engine_harris_guess = harris_guess
if (len_trim(ri_aux_basis) .gt. 0) then
   engine_df_aux_basis = resolve_basis_name(ri_aux_basis)
else if (engine_use_df) then
   engine_df_aux_basis = auto_ri_aux_basis(resolve_basis_name(baselable), engine_use_df_k)
   if (.not. aux_basis_covers_all_atoms(trim(engine_df_aux_basis), ncenters, atomchg)) then
      if (engine_use_df_k) then
         engine_df_aux_basis = 'def2universaljkfit'
      else
         engine_df_aux_basis = 'def2universaljfit'
      endif
      if (.not. aux_basis_covers_all_atoms(trim(engine_df_aux_basis), ncenters, atomchg)) then
         write(*,'(A)') ' [FATAL] no published RI auxiliary basis covers every element in this'// &
                        ' molecule (tried def2-universal-J/JKFIT). Set ri_aux_basis explicitly'// &
                        ' or run with J=''exact''.'
         stop 1
      endif
   endif
else
   engine_df_aux_basis = ""
endif
engine_estimate_only = estimate_only
engine_mem_cap_bytes = 0_8
if (mem_cap_gb .gt. 0.0d0) engine_mem_cap_bytes = int(mem_cap_gb*1024_8*1024_8*1024_8, 8)
engine_avail_at_start_bytes = get_memory_budget_bytes()
call engine_startup_mark("memory_budget")

block
use mod_basis_files, only: set_basis_dir_override
call set_basis_dir_override(basedir)
call engine_startup_mark("basis_dir")
end block

mem_grid_gb = 0.0d0
mem_2e_gb = 0.0d0
energy_out = 0.0d0
econv = 0.0d0
iconv = 0
force_out = 0.0d0
MLcharge_out = 0.0d0

E = 0

allocate(atoms(Natoms))

nPointCharges = npc
allocate(pointcharge_q(nPointCharges))
allocate(pointcharge_coor(3,nPointCharges))
do i = 1,nPointCharges
   pointcharge_q(i) = pc_charge(i)
   pointcharge_coor(:,i) = pc_coord(i,:)
enddo

coreChg = 0
do i = 1,Natoms
  atoms(i)%coor = coord(i,:)
  atoms(i)%charge = atomchg(i)
  if (len_trim(atombasis(i)) == 0) then
     write (atoms(i)%base,"(I0.2,A1,A27)")  atomchg(i),'-',resolve_basis_name(baselable)
  else
     write (atoms(i)%base,"(I0.2,A1,A27)")  atomchg(i),'-',resolve_basis_name(atombasis(i))
  endif
  if (len_trim(atomecp(i)) .gt. 0) then
     write (atoms(i)%ecpbase,"(I0.2,A1,A27)")  atomchg(i),'-',resolve_basis_name(atomecp(i))
  elseif (len_trim(ecplabel) .gt. 0) then
     write (atoms(i)%ecpbase,"(I0.2,A1,A27)")  atomchg(i),'-',resolve_basis_name(ecplabel)
  else
     atoms(i)%ecpbase = atoms(i)%base
     if (index(adjustl(atoms(i)%base(4:)), 'def2') .eq. 1) &
        atoms(i)%ecpbase = atoms(i)%base(1:3)//'def2tzvpp'
  endif
  if (index(adjustl(atoms(i)%ecpbase(4:)), 'def2') .eq. 1) &
     atoms(i)%ecpbase = atoms(i)%ecpbase(1:3)//'def2tzvpp'
  atoms(i)%ecpCoreElec = ecp_core_electrons_prescan(atoms(i)%ecpbase)
  if (atoms(i)%ecpCoreElec > 0 .and. engine_verbose .ge. 2) then
     print *, "DEBUG: atom", i, "Z=", atoms(i)%charge, "base=", trim(atoms(i)%ecpbase), &
              "n_core=", atoms(i)%ecpCoreElec, "Z_eff=", atoms(i)%charge - atoms(i)%ecpCoreElec
  endif
  coreChg = atoms(i)%charge - atoms(i)%ecpCoreElec + coreChg
enddo

if (Multi .lt. 1) then
   print '("Charge/multiplicity error: imult = ",I0," is not a valid spin multiplicity (must be >= 1)")', Multi
   stop 1
endif
if (coreChg-Charge .lt. 0) then
   print '("Charge/multiplicity error: icharge = ",I0," leaves a negative electron count (",I0,&
         &" active electrons available after ECP core removal)")', Charge, coreChg
   stop 1
endif
if (mod(coreChg-Charge-(Multi-1), 2) .ne. 0) then
   print '("Charge/multiplicity error: icharge = ",I0," / imult = ",I0," are inconsistent - ",&
         &I0," active electrons cannot supply ",I0," unpaired electron(s) (parity mismatch, try imult = ",I0,")")', &
         Charge, Multi, coreChg-Charge, Multi-1, mod(coreChg-Charge,2)+1
   stop 1
endif
if (Multi-1 .gt. coreChg-Charge) then
   print '("Charge/multiplicity error: imult = ",I0," requests more unpaired electrons than the ",&
         &I0," active electrons available (icharge = ",I0,")")', Multi, coreChg-Charge, Charge
   stop 1
endif

n_alpha = 0.5*(Multi-Charge+coreChg-1)
n_beta = 0.5*(-Multi-Charge+coreChg+1)
if (opts%uks .and. n_alpha == n_beta) Multi = 2

if (.not. engine_quiet) then
   call engine_startup_mark("before n_alpha")
   print '("n_alpha = ", I0, "   n_beta = ", I0, "   (# occupied alpha/beta electrons)")', n_alpha, n_beta
endif

if (.not. engine_quiet) then
   block
      character(len=2) :: dsym(ncenters)
      integer :: ids
      do ids = 1, ncenters
         dsym(ids) = elem_symbol(atomchg(ids))
      enddo
      call print_distance_matrix(6, ncenters, dsym, transpose(coord))
   end block
endif

allocate(linkMat(natoms,natoms))
linkMat = 0
do i = 1,natoms
   do j = i+1,natoms
       dist = sum((atoms(i)%coor - atoms(j)%coor)**2)
       if (dist**0.5 < 1.2*(covrad(atoms(i)%charge) + covrad(atoms(j)%charge))) then
         linkMat(i,j) = 1
         linkMat(j,i) = 1
       endif
   enddo
enddo

block
   character(len=64) :: hessenv
   integer :: hessstat
   hessenv = ""
   call get_environment_variable("ENGINE_HESSIAN_TOTAL_DEBUG", hessenv, status=hessstat)
   hessian_run = (hessstat == 0) .or. opts%freq .or. opts%polar .or. opts%nmr .or. opts%tddft_nroots > 0
end block

call engine_startup_mark("before integrals_init")
call prof_start("basis_and_integrals")
call integrals_init(info)
call prof_stop("basis_and_integrals")
call engine_startup_mark("integrals_init")

if (engine_use_df) then
   call integrals_estimate_aux_size()
   mem_2e_gb = real(nconts,8)*real(nconts,8)*real(nContsAux,8)*8.0d0/1024.0d0**3
else
   mem_2e_gb = real(nRec,8)*8.0d0/1024.0d0**3
endif

if (HF_exchange_frac .lt. 1.0d0) then
   if (.not. engine_quiet) print *, '[GRID]'
   block
      block
         character(len=32) :: dyngridenv, cacheenv
         call get_environment_variable("ENGINE_GRID_CACHE", cacheenv)
         if (len_trim(cacheenv) .gt. 0) then
            xc_direct_mode = (trim(cacheenv) .eq. "0")
         else
            xc_direct_mode = .true.
         endif
         xcgrid_nmr_active = opts%nmr
         call xcgrid_tune_from_env(functional_grid)
         has_warm_start = (size(dens_in_a,1) .gt. 0) .or. molden_read
         xcgrid_theory_refined = .false.
         xcgrid_level = XCGRID_FINE
         if (xc_uses_vv10()) then
            call get_environment_variable("ENGINE_VV10_DYNGRID", dyngridenv)
            if (len_trim(dyngridenv) .gt. 0) then
               vv10_dynamic = (trim(dyngridenv) .eq. "1") .and. .not. engine_vv10_nonself
               if (trim(dyngridenv) .eq. "1" .and. engine_vv10_nonself) &
                  print *,"ENGINE_VV10_DYNGRID ignored: not compatible with vv10_nonself"
            else
               vv10_dynamic = (nconts .gt. VV10_DYNGRID_NCONTS) .and. .not. engine_vv10_nonself
               if (vv10_dynamic) print *,"VV10 delayed self-consistent switch by default: nconts=",nconts, &
                    " >",VV10_DYNGRID_NCONTS," (override with ENGINE_VV10_DYNGRID=0)"
            endif
            vv10_active_now = (.not. vv10_dynamic) .or. has_warm_start
            vv10_just_activated = .false.
            vv10_flush_on_switch = .true.
            if (vv10_dynamic .and. .not. has_warm_start) then
               call get_environment_variable("ENGINE_VV10_DYNGRID_PRMS", dyngridenv)
               if (len_trim(dyngridenv) .gt. 0) read(dyngridenv,*) vv10_switch_prms
               call get_environment_variable("ENGINE_VV10_DYNGRID_NOFLUSH", dyngridenv)
               if (len_trim(dyngridenv) .gt. 0) vv10_flush_on_switch = (trim(dyngridenv) .ne. "1")
               print *,"VV10 dynamic: starts OFF, activates self-consistently at dP <",vv10_switch_prms, &
                       " flush_on_switch=",vv10_flush_on_switch
            endif
         else
            vv10_dynamic    = .false.
            vv10_active_now = .true.
         endif
      end block
      if (engine_verbose .ge. 2) then
         block
            integer(8) :: diag_rss, diag_hwm
            call get_process_memory_bytes(diag_rss, diag_hwm)
            print '(A)', "  DIAG before initial grid_gen: RSS="//trim(fmt_gb(real(diag_rss,8)/1024**3))//"GB"
         end block
      endif
      call prof_start("grid_gen")
      call load_or_generate_grid(info)
      call prof_stop("grid_gen")
      mem_grid_gb = real(size(Grids),8)*real(nconts,8)*4.0d0*4.0d0/1024.0d0**3
      call prof_start("gto_eval")
      call GTOeval(info)
      call prof_stop("gto_eval")
      if (engine_verbose .ge. 2) then
         block
            integer(8) :: diag_rss, diag_hwm
            call get_process_memory_bytes(diag_rss, diag_hwm)
            print '(A)', "  DIAG after initial grid_gen+GTOeval (coarse grid): RSS="//trim(fmt_gb(real(diag_rss,8)/1024**3))//"GB"
         end block
      endif
   end block
   if (.not. engine_quiet) then
      print *, '[GRIDEND]'
      print *
   endif
endif

if (estimate_only) then
   if (engine_use_df) then
      block
         use mod_mem_predict, only: predict_energy_peak_bytes, predict_force_extra_bytes, &
                                     GRID_RAGGED_FRACTION_EST
         integer(8) :: grid_cache_bytes_est, energy_need_bytes, force_need_bytes
         logical :: need_k_est
         need_k_est = (HF_exchange_frac .gt. 1.0d-12)
         grid_cache_bytes_est = 0_8
         if (grid_cache_mode) &
            grid_cache_bytes_est = int(int(max_batch_nsig,8)*int(grid_cache_capacity,8)*4_8*4_8 &
                                       *GRID_RAGGED_FRACTION_EST, 8) &
                                  + 8_8*int(ngrids,8)*8_8
         energy_need_bytes = predict_energy_peak_bytes(nConts, Natoms, nContsAux, need_k_est, &
                                                         .true., grid_cache_bytes_est, &
                                                         rs_omega_active=(RS_omega .gt. 0.0d0 .and. engine_use_df_k))
         print '(A)', "  [estimate_only] predicted total peak, energy phase STORE="// &
                 trim(fmt_gb(real(energy_need_bytes,8)/1024**3))// &
                 "GB  available="//trim(fmt_gb(real(engine_avail_at_start_bytes,8)/1024**3))// &
                 "GB  would pick="//trim(merge("DIRECT","STORE ",energy_need_bytes > engine_avail_at_start_bytes))
         if (do_force) then
            force_need_bytes = energy_need_bytes &
                              + predict_force_extra_bytes(nConts, nContsAux, need_k_est, .true.)
            print '(A)', "  [estimate_only] predicted total peak, force phase STORE="// &
                    trim(fmt_gb(real(force_need_bytes,8)/1024**3))// &
                    "GB  available="//trim(fmt_gb(real(engine_avail_at_start_bytes,8)/1024**3))// &
                    "GB  would pick="//trim(merge("DIRECT","STORE ",force_need_bytes > engine_avail_at_start_bytes))
         endif
      end block
   endif
   call prof_report()
   return
endif

if (cosmo_on .and. cosmo_smd) then
   if (trim(cosmo_cavity_type) /= 'gepol') then
      print '("COSMO error: cosmo_smd=.true. requires cosmo_cavity_type=''gepol'' (got ''",A,"'')")', &
            trim(cosmo_cavity_type)
      stop 1
   endif
   if (len_trim(cosmo_solvent) == 0) then
      print '("COSMO error: cosmo_smd=.true. requires a non-empty cosmo_solvent (named solvent)")'
      stop 1
   endif
endif
if (cosmo_on) then
   block
      integer :: znum(Natoms), ecpce(Natoms), ia
      real(8) :: coor_ang(3,Natoms)
      real(8) :: eps_use
      do ia = 1,Natoms
         znum(ia) = atoms(ia)%charge
         ecpce(ia) = atoms(ia)%ecpCoreElec
         coor_ang(:,ia) = atoms(ia)%coor
      enddo
      if (cosmo_smd) then
         call cosmo_solvent_smd_params(trim(cosmo_solvent), smd_n_val, smd_alpha_val, smd_beta_val, &
                                        smd_gamma_val, smd_phi_val, smd_psi_val)
         smd_is_water = (trim(cosmo_solvent) == 'water')
         eps_use = cosmo_solvent_epsilon(trim(cosmo_solvent))
         call cosmo_init(Natoms, znum, ecpce, coor_ang, eps_use, cosmo_radii_scale, &
                          cosmo_avg_area, trim(cosmo_cavity_type), cosmo_rsolv, &
                          cosmo_ks_nseg, cosmo_ks_nface, nConts, real(icharge,8), smd_alpha_val)
      else
         call cosmo_init(Natoms, znum, ecpce, coor_ang, cosmo_epsilon, cosmo_radii_scale, &
                          cosmo_avg_area, trim(cosmo_cavity_type), cosmo_rsolv, &
                          cosmo_ks_nseg, cosmo_ks_nface, nConts, real(icharge,8))
      endif
      call cosmo_set_sigma_profile_debug(cosmo_sigma_rav, cosmo_sigma_profile_file)
   end block
endif

if (.not. engine_quiet) print *, '[SCF]'
Gtype = 3
allocate(Pa_chk(nConts,nConts),Pb_chk(nConts,nConts))
chk_ok = .false.
if (size(dens_in_a,1) .eq. nConts .and. size(dens_in_a,2) .eq. nConts .and. &
    size(dens_in_b,1) .eq. nConts .and. size(dens_in_b,2) .eq. nConts) then
   Pa_chk = dens_in_a
   Pb_chk = dens_in_b
   chk_ok = .true.
else if (size(dens_in_a,1) .gt. 0) then
   print *,"dens_in_a/dens_in_b shape mismatch with current basis (need ",nConts,"x",nConts, &
           "), falling back to molden_read/SAD"
endif
if (.not. chk_ok .and. molden_read) call read_molden(molden_read_file, Pa_chk, Pb_chk, chk_ok)
if (chk_ok) Gtype = 4
call engine_startup_mark("before guess")
call guess(info,Gtype,Pa_chk,Pb_chk)
call engine_startup_mark("guess")
deallocate(Pa_chk,Pb_chk)

if (xcgrid_audit_only) then
   if (.not. chk_ok) then
      print '(A)', "[GRIDAUDIT] ERROR: no density supplied "// &
                   "(need &checkpoint chk_read / &molden molden_read) - nothing to integrate"
      iconv = 0
   else
      block
         real(8) :: exc_audit
         real(8), allocatable :: fdum_a(:,:), fdum_b(:,:)
         integer :: nel, u, ios2, jj, dmord2(4)
         real(8) :: rp, dmb2(3)
         character(len=512) :: scanfile, line
         logical :: have_scan
         nel = n_alpha + n_beta
         allocate(fdum_a(nConts,nConts), fdum_b(nConts,nConts))
         call get_environment_variable("ENGINE_XC_GRID_SCAN", scanfile)
         have_scan = .false.
         u = -1
         if (len_trim(scanfile) > 0) then
            open(newunit=u, file=trim(scanfile), status='old', &
                 action='read', iostat=ios2)
            if (ios2 == 0) then
               have_scan = .true.
            else
               print '(A)', "[GRIDAUDIT] ERROR: cannot open "// &
                            trim(scanfile)//" - auditing the current grid only"
            endif
         endif
         do
            if (have_scan) then
               read(u, '(A)', iostat=ios2) line
               if (ios2 /= 0) exit
               if (len_trim(line) == 0) cycle
               read(line, *, iostat=ios2) rp, dmb2, dmord2
               if (ios2 /= 0 .or. rp <= 0.0d0) cycle
               if (.not. all(dmb2(2:3) > dmb2(1:2))) cycle
               if (any(.not. xcgrid_lebedev_supported(dmord2))) cycle
               xcgrid_radpot_rowscale = rp
               xcgrid_prune_dmin_bounds = dmb2
               do jj = 1, 7
                  xcgrid_prune_dmin_orders(:, jj) = dmord2
               enddo
               call xcgrid_free_derived()
               call load_or_generate_grid(info)
               call GTOeval(info)
            endif
            call DFT_calc(exc_audit, fdum_a, fdum_b, info, .false.)
            print '(A)', "[GRIDAUDIT]"
            if (have_scan) &
               print '(A,F8.4,A,3(1X,F8.4),A,4(1X,I0))', &
                  "  cand          = ", rp, "  |", dmb2, "  |", dmord2
            print '(A,I0)',        "  grid points   = ", ngrids
            print '(A,F18.10)',    "  int(rho)      = ", xcgrid_int_rho
            print '(A,I0)',        "  N_elec        = ", nel
            print '(A,ES14.6)',    "  d_rho         = ", xcgrid_int_rho - real(nel,8)
            print '(A,ES14.6)',    "  d_rho_rel     = ", (xcgrid_int_rho - real(nel,8))/max(1,nel)
            print '(A,ES14.6)',    "  d_trPS        = ", &
               sum(Pa*S) + sum(Pb*S) - real(nel,8)
            print '(A,F18.10)',    "  Exc(this grid)= ", exc_audit
            print '(A)', "[GRIDAUDITEND]"
            if (.not. have_scan) exit
         enddo
         if (have_scan) close(u)
         deallocate(fdum_a, fdum_b)
         iconv = 1
      end block
   endif
   energy_out = 0.0d0
   force_out = 0.0d0
   MLcharge_out = 0.0d0
   econv = 0.0d0
   mem_grid_gb = 0.0d0
   mem_2e_gb = 0.0d0
   if (.not. allocated(dens_out_a)) allocate(dens_out_a(nConts,nConts))
   if (.not. allocated(dens_out_b)) allocate(dens_out_b(nConts,nConts))
   dens_out_a = Pa ; dens_out_b = Pb
   return
endif

select case (scf_conv_level)
case (1)
   Emax = 1.00E-08
   Pmax = 1.00E-06
case (2)
   Emax = 1.00E-09
   Pmax = 1.00E-08
case default
   Emax = 1.00E-06
   Pmax = 1.00E-05
end select
itmax =150
call engine_startup_mark("before SCFcycle")
call prof_start("scf_total")
block
   call SCFcycle(info,Emax,Pmax,itmax,iconv,econv,do_force,opts%nmr .or. opts%cosx_hi_reconverge)
end block
call prof_stop("scf_total")
block
   use mod_scf_history, only: scf_hist_niter, scf_hist_last_prms, scf_hist_last_de, scf_hist_last_exc, &
        scf_hist_last_ecoul, scf_hist_last_exact_exchange, scf_hist_last_hcore, scf_hist_last_erep
   use mod_vv10, only: vv10_report
   opts%scf_niter = scf_hist_niter() ; opts%scf_dP = scf_hist_last_prms() ; opts%scf_dE = scf_hist_last_de()
   opts%e_hcore = scf_hist_last_hcore() ; opts%e_coul = scf_hist_last_ecoul()
   opts%e_exx = scf_hist_last_exact_exchange() ; opts%e_xc = scf_hist_last_exc() ; opts%e_nuc = scf_hist_last_erep()
   call vv10_report(opts%e_vv10, opts%vv10_used)
end block
block
   use mod_tddft, only: tddft_stability
   use mod_cosmo_state, only: cosmo_enabled
   character(len=8) :: sme
   integer :: smode, smios, snit
   real(8) :: lam_t
   smode = 0
   call get_environment_variable("ENGINE_STABILITY", sme)
   if (len_trim(sme) > 0) read(sme, *, iostat=smios) smode
   if (smode >= 1 .and. Multi == 1 .and. iconv == 1 .and. .not. cosmo_enabled) then
      call tddft_stability(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, .true., lam_t, snit)
      print '(A,ES11.3,A,I2,A)', ' STABILITY (RKS -> UKS, triplet-type): lowest eigenvalue of A+B =', lam_t, &
            '  (', snit, ' Davidson steps)'
      if (lam_t < -1.0d-3) then
         print '(A)', ' STABILITY (RKS -> UKS): the restricted solution is UNSTABLE toward an unrestricted one'
         if (smode >= 2) then
            opts%uks_restart = .true.
            opts%uks_lambda = lam_t
            energy_out = E
            force_out = 0.0d0 ; MLcharge_out = 0.0d0
            mem_grid_gb = 0.0d0 ; mem_2e_gb = 0.0d0
            allocate(dens_out_a(nConts,nConts), dens_out_b(nConts,nConts))
            dens_out_a = Pa ; dens_out_b = Pb
            return
         endif
         print '(A)', '                        (stability = 2 redoes the run unrestricted)'
      else
         print '(A)', ' STABILITY (RKS -> UKS): stable'
      endif
   endif
end block
block
   character(len=8) :: esm
   character(len=16) :: spid
   call get_environment_variable("ENGINE_SMAPS_DUMP", esm)
   if (trim(esm) .eq. "1") then
      write(spid,'(I0)') getpid()
      call execute_command_line("cp /proc/"//trim(spid)//"/smaps smaps_scf_end.txt")
      call smaps_malloc_info("malloc_info_scf_end.xml")
      call execute_command_line("cp /proc/"//trim(spid)//"/smaps smaps_scf_end_trim.txt")
      call release_free_memory()
      call execute_command_line("cp /proc/"//trim(spid)//"/smaps smaps_scf_end_rel.txt")
   endif
end block
block
   use mod_integrals, only: integrals_build_exchange_cosx_batch_c
   use omp_lib, only: omp_get_wtime
   character(len=16) :: eb
   integer :: ndb, ib2, no2, kq, oq
   real(8),allocatable :: Lb(:,:,:), KCb(:,:,:)
   real(8) :: tb2
   call get_environment_variable("ENGINE_COSXK_BENCH", eb)
   if (len_trim(eb) .gt. 0) then
      read(eb,*) ndb
      no2 = n_alpha
      allocate(Lb(nConts,no2,ndb), KCb(nConts,no2,ndb))
      do ib2 = 1,ndb
         do oq = 1,no2
            do kq = 1,nConts
               Lb(kq,oq,ib2) = 1.0d-2*C_a(kq,no2+1+mod(oq+ib2,nConts-no2))*sin(0.1d0*(kq+oq+ib2))
            enddo
         enddo
      enddo
      do ib2 = 1,2
         tb2 = omp_get_wtime()
         call integrals_build_exchange_cosx_batch_c(nConts, ndb, Lb, no2, C_a(:,1:no2), KCb)
         print '(A,I0,A,F9.3,A,ES24.16)', " COSXK BENCH nd=", ndb, "  wall ", omp_get_wtime()-tb2, &
               "  sum|KC| ", sum(abs(KCb))
      enddo
      call flush(6)
      stop
   endif
end block
block
   use mod_integrals, only: integrals_build_coulomb_df_batch
   use omp_lib, only: omp_get_wtime
   character(len=16) :: eb
   integer :: ndb, ib2, tb2i, tb2j
   real(8),allocatable :: Db(:,:,:), Jb(:,:,:)
   real(8) :: tb2
   call get_environment_variable("ENGINE_DFJ_BENCH", eb)
   if (len_trim(eb) .gt. 0) then
      read(eb,*) ndb
      allocate(Db(nConts,nConts,ndb), Jb(nConts,nConts,ndb))
      do ib2 = 1,ndb
         do tb2i = 1,nConts
            Db(:,tb2i,ib2) = Pa(:,tb2i)*sin(0.1d0*ib2*[(dble(tb2i+tb2j), tb2j=1,nConts)])
         enddo
      enddo
      do ib2 = 1,2
         tb2 = omp_get_wtime()
         call integrals_build_coulomb_df_batch(nConts, ndb, Db, Jb)
         print '(A,I0,A,F9.3,A,ES24.16)', " DFJ BENCH nd=", ndb, "  wall ", omp_get_wtime()-tb2, &
               "  sum|J| ", sum(abs(Jb))
      enddo
      call flush(6)
      stop
   endif
end block

block
   use GRID_info, only: Grids, ngrids
   character(len=256) :: gf
   integer :: ug, ig2
   call get_environment_variable("ENGINE_XC_GRID_DUMP", gf)
   if (len_trim(gf) .gt. 0) then
      open(newunit=ug, file=trim(gf), status='replace', action='write')
      do ig2 = 1,ngrids
         write(ug,'(4ES25.16)') Grids(ig2)%coor, Grids(ig2)%weight
      enddo
      close(ug)
   endif
end block

block
   use GRID_info, only: Grids, ngrids, xc_dbg_pt
   character(len=256) :: gf
   integer :: ug, ig2
   call get_environment_variable("ENGINE_XC_POINT_DUMP", gf)
   if (len_trim(gf) .gt. 0 .and. allocated(xc_dbg_pt)) then
      open(newunit=ug, file=trim(gf), status='replace', action='write')
      do ig2 = 1,ngrids
         write(ug,'(9ES25.16E3)') Grids(ig2)%coor, Grids(ig2)%weight, xc_dbg_pt(:,ig2)
      enddo
      close(ug)
   endif
end block

block
   use mod_integrals, only: integrals_build_coulomb_df_batch_c
   use omp_lib, only: omp_get_wtime
   character(len=16) :: eb
   integer :: ndb, ib2, nop2
   real(8),allocatable :: Lb2(:,:,:), Jc2(:,:,:)
   real(8) :: tb2
   call get_environment_variable("ENGINE_DFJC_BENCH", eb)
   if (len_trim(eb) .gt. 0) then
      read(eb,*) ndb
      nop2 = nConts/5
      allocate(Lb2(nConts,nop2,ndb), Jc2(nConts,nop2,ndb))
      do ib2 = 1,ndb
         Lb2(:,:,ib2) = Pa(:,1:nop2)*sin(0.1d0*ib2)
      enddo
      do ib2 = 1,2
         tb2 = omp_get_wtime()
         call integrals_build_coulomb_df_batch_c(nConts, ndb, Lb2, nop2, Pa(:,nConts-nop2+1:nConts), Jc2)
         print '(A,I0,A,F9.3,A,ES24.16)', " DFJC BENCH nd=", ndb, "  wall ", omp_get_wtime()-tb2, &
               "  sum|JC| ", sum(abs(Jc2))
      enddo
      call flush(6)
      stop
   endif
end block

allocate(dens_out_a(nConts,nConts), dens_out_b(nConts,nConts))
dens_out_a = Pa
dens_out_b = Pb
if (molden_write .and. iconv .eq. 1) call write_molden(molden_file)
if (iconv .eq. 1 .and. (.not. engine_quiet .or. opts%cube_every) .and. (opts%cube_den .or. opts%cube_esp)) &
   call write_cube_files(trim(opts%freq_base), opts%cube_den, opts%cube_esp, opts%cube_esp_exact)

block
   use mod_integrals, only: integrals_hessian_debug_check2e_direct
   character(32) :: envval
   integer :: envstat
   call get_environment_variable("ENGINE_HESSIAN_DEBUG2E_DIRECT", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call integrals_hessian_debug_check2e_direct(Nconts, Natoms, Pa, Pb, HF_exchange_frac)
      stop
   endif
end block

block
   interface
      subroutine vv10_nlc_weight_hessian_selftest()
      end subroutine vv10_nlc_weight_hessian_selftest
   end interface
   character(32) :: envwh
   integer :: envstatwh
   call get_environment_variable("ENGINE_VV10_WHESS_TEST", envwh, status=envstatwh)
   if (envstatwh == 0) then
      call vv10_nlc_weight_hessian_selftest()
      stop
   endif
end block

block
   use mod_integrals, only: integrals_build_exchange_cosx, integrals_build_exchange_cosx_lr, &
                            cosx_set_response_mode
   use mod_exchange, only: cosx_enabled, RS_omega
   character(32) :: ksenv
   integer :: ksstat, ki, kj, nseed
   integer,allocatable :: seedv(:)
   real(8),allocatable :: D1(:,:), D2(:,:), K1a(:,:), K1b(:,:), K2a(:,:), K2b(:,:)
   real(8) :: t12, t21
   call get_environment_variable("ENGINE_COSX_KSYM_TEST", ksenv, status=ksstat)
   if (ksstat == 0 .and. cosx_enabled .and. iconv == 1) then
      call random_seed(size=nseed)
      allocate(seedv(nseed)) ; seedv = 12345
      call random_seed(put=seedv)
      allocate(D1(Nconts,Nconts), D2(Nconts,Nconts), K1a(Nconts,Nconts), K1b(Nconts,Nconts))
      allocate(K2a(Nconts,Nconts), K2b(Nconts,Nconts))
      call random_number(D1) ; call random_number(D2)
      D1 = 0.5d0*(D1 + transpose(D1)) - 0.5d0 ; D2 = 0.5d0*(D2 + transpose(D2)) - 0.5d0
      do ki = 1,Nconts
         do kj = 1,Nconts
            D1(ki,kj) = D1(ki,kj)*Pa(ki,kj) ; D2(ki,kj) = D2(ki,kj)*Pa(ki,kj)
         enddo
      enddo
      call cosx_set_response_mode(.true.)
      call integrals_build_exchange_cosx(Nconts, D1, D1, K1a, K1b, .false.)
      call integrals_build_exchange_cosx(Nconts, D2, D2, K2a, K2b, .false.)
      t12 = sum(D1*K2a) ; t21 = sum(D2*K1a)
      print '(A,2ES24.16,A,ES10.2)', ' COSX KSYM full: Tr(D1K2), Tr(D2K1) =', t12, t21, &
            '  rel asym =', abs(t12-t21)/max(abs(t12),1.0d-300)
      if (RS_omega > 0.0d0) then
         call integrals_build_exchange_cosx_lr(Nconts, D1, D1, K1a, K1b)
         call integrals_build_exchange_cosx_lr(Nconts, D2, D2, K2a, K2b)
         t12 = sum(D1*K2a) ; t21 = sum(D2*K1a)
         print '(A,2ES24.16,A,ES10.2)', ' COSX KSYM LR:   Tr(D1K2), Tr(D2K1) =', t12, t21, &
               '  rel asym =', abs(t12-t21)/max(abs(t12),1.0d-300)
      endif
      call cosx_set_response_mode(.false.)
      call flush(6)
      stop
   endif
end block

block
   use mod_integrals, only: integrals_force_exchange_cosx, integrals_build_exchange_cosx, &
                            cosx_set_response_mode, integrals_build_exchange_cosx_lr
   use mod_exchange, only: cosx_enabled, RS_omega
   character(256) :: fxfile
   real(8) :: fx_omega
   integer :: fxstat, ufx, ifx, jfx
   real(8),allocatable :: Dfa(:,:), Dfb(:,:), fxa(:,:), fxb(:,:), Kfa(:,:), Kfb(:,:)
   call get_environment_variable("ENGINE_COSX_FIXD_ORACLE", fxfile, status=fxstat)
   if (fxstat == 0 .and. cosx_enabled .and. iconv == 1) then
      allocate(Dfa(Nconts,Nconts), Dfb(Nconts,Nconts), fxa(Natoms,3), fxb(Natoms,3))
      allocate(Kfa(Nconts,Nconts), Kfb(Nconts,Nconts))
      call get_environment_variable("ENGINE_COSX_FIXD_READ", fxfile, status=fxstat)
      if (fxstat == 0) then
         open(newunit=ufx, file=trim(fxfile), status='old', action='read')
         do ifx = 1,Nconts
            read(ufx,*) (Dfa(ifx,jfx), jfx=1,Nconts)
         enddo
         do ifx = 1,Nconts
            read(ufx,*) (Dfb(ifx,jfx), jfx=1,Nconts)
         enddo
         close(ufx)
      else
         Dfa = Pa ; Dfb = Pb
         open(newunit=ufx, file='cosx_fixd_dens.txt', status='replace', action='write')
         do ifx = 1,Nconts
            write(ufx,'(*(ES25.17))') Dfa(ifx,:)
         enddo
         do ifx = 1,Nconts
            write(ufx,'(*(ES25.17))') Dfb(ifx,:)
         enddo
         close(ufx)
      endif
      fx_omega = 0.0d0
      call get_environment_variable("ENGINE_COSX_FIXD_OMEGA", fxfile, status=fxstat)
      if (fxstat == 0) fx_omega = RS_omega
      call integrals_force_exchange_cosx(Nconts, Natoms, Dfa, Dfb, fx_omega, fxa, fxb)
      call get_environment_variable("ENGINE_COSX_GRAD_XCHK", fxfile, status=fxstat)
      if (fxstat == 0) then
         block
            use mod_integrals, only: integrals_compute_force
            use mod_exchange, only: HF_exchange_frac
            real(8),allocatable :: dS_(:,:,:), dH_(:,:,:), dJ_(:,:,:), dKa_(:,:,:), dKb_(:,:,:), gex(:,:)
            integer :: ug, ia, ll, kk
            allocate(dS_(Nconts,Nconts,3), dH_(Nconts,Nconts,3), dJ_(Nconts,Nconts,3))
            allocate(dKa_(Nconts,Nconts,3), dKb_(Nconts,Nconts,3), gex(Natoms,3))
            call integrals_compute_force(Nconts, Dfa, Dfb, dS_, dH_, dJ_, dKa_, dKb_)
            gex = 0.0d0
            ll = 0
            do ia = 1,Natoms
               do ifx = 1,atoms(ia)%nconts
                  ll = ll + 1
                  do kk = 1,Nconts
                     gex(ia,:) = gex(ia,:) + 2.0d0*HF_exchange_frac* &
                        (Dfa(ll,kk)*dKa_(ll,kk,:) + Dfb(ll,kk)*dKb_(ll,kk,:))
                  enddo
               enddo
            enddo
            open(newunit=ug, file='cosx_gradxchk.txt', status='replace', action='write')
            do ia = 1,Natoms
               write(ug,'(6ES25.17)') 2.0d0*HF_exchange_frac*(fxa(ia,:)+fxb(ia,:)), gex(ia,:)
            enddo
            close(ug)
            print '(A,2ES12.4)', ' COSX GRAD XCHK fixed-D: RMS/max |cosx-exact| =', &
                 sqrt(sum((2.0d0*HF_exchange_frac*(fxa+fxb)-gex)**2)/dble(3*Natoms)), &
                 maxval(abs(2.0d0*HF_exchange_frac*(fxa+fxb)-gex))
            deallocate(dS_, dH_, dJ_, dKa_, dKb_, gex)
         end block
      endif
      call get_environment_variable("ENGINE_COSX_HESS_TEST", fxfile, status=fxstat)
      if (fxstat == 0) then
         block
            use mod_integrals, only: cosx_hess_exchange
            real(8),allocatable :: h2a(:,:,:,:), h2b(:,:,:,:), kxa_(:,:,:,:), kxb_(:,:,:,:)
            integer :: ukx
            allocate(h2a(3,Natoms,3,Natoms), h2b(3,Natoms,3,Natoms))
            allocate(kxa_(Nconts,Nconts,3,Natoms), kxb_(Nconts,Nconts,3,Natoms))
            call cosx_hess_exchange(Nconts, Natoms, Dfa, Dfb, 0.0d0, h2a, h2b, kxa_, kxb_)
            open(newunit=ukx, file='cosx_hess_kx.txt', status='replace', action='write')
            do jfx = 1,Natoms
               do ifx = 1,3
                  write(ukx,'(*(ES25.17))') kxa_(:,:,ifx,jfx)
                  write(ukx,'(*(ES25.17))') kxb_(:,:,ifx,jfx)
               enddo
            enddo
            close(ukx)
            deallocate(h2a, h2b, kxa_, kxb_)
         end block
      endif
      call cosx_set_response_mode(.true.)
      if (fx_omega > 0.0d0) then
         call integrals_build_exchange_cosx_lr(Nconts, Dfa, Dfb, Kfa, Kfb)
      else
         call integrals_build_exchange_cosx(Nconts, Dfa, Dfb, Kfa, Kfb, .false.)
      endif
      call cosx_set_response_mode(.false.)
      open(newunit=ufx, file='cosx_fixd_out.txt', status='replace', action='write')
      write(ufx,'(2I8)') Natoms, Nconts
      write(ufx,'(2ES25.17)') sum(Dfa*Kfa), sum(Dfb*Kfb)
      do ifx = 1,Natoms
         write(ufx,'(6ES25.17)') fxa(ifx,:), fxb(ifx,:)
      enddo
      do ifx = 1,Nconts
         write(ufx,'(*(ES25.17))') Kfa(ifx,:)
      enddo
      do ifx = 1,Nconts
         write(ufx,'(*(ES25.17))') Kfb(ifx,:)
      enddo
      close(ufx)
      deallocate(Dfa,Dfb,fxa,fxb,Kfa,Kfb)
      stop
   endif
end block

block
   use mod_integrals, only: integrals_force_exchange_cosx, cosx_hess_weight_term, &
                            cosx_hess_a_term
   use mod_exchange, only: cosx_enabled
   character(32) :: envch
   integer :: envstatch
   real(8),allocatable :: fka(:,:), fkb(:,:), d2w(:,:,:,:), d2a(:,:,:,:), d2b(:,:,:,:)
   call get_environment_variable("ENGINE_COSX_HESS_W_TEST", envch, status=envstatch)
   if (envstatch == 0 .and. cosx_enabled .and. iconv == 1) then
      allocate(fka(Natoms,3), fkb(Natoms,3), d2w(3,Natoms,3,Natoms))
      allocate(d2a(3,Natoms,3,Natoms), d2b(3,Natoms,3,Natoms))
      call integrals_force_exchange_cosx(Nconts, Natoms, Pa, Pb, 0.0d0, fka, fkb)
      call cosx_hess_weight_term(Natoms, d2w)
      call cosx_hess_a_term(Natoms, d2a, d2b)
      print '(A,ES12.4)', " COSX HESS W-TERM: max|d2E_w| = ", maxval(abs(d2w))
      print '(A,3ES12.4)', " COSX HESS W-TERM: transl sum (atom1) = ", &
            sum(d2w(1,1,:,:)), sum(d2w(2,1,:,:)), sum(d2w(3,1,:,:))
      print '(A,ES12.4)', " COSX HESS A-TERM: max|d2E_a+d2E_b| = ", maxval(abs(d2a+d2b))
      print '(A,3ES12.4)', " COSX HESS A-TERM: transl sum (atom1) = ", &
            sum(d2a(1,1,:,:)+d2b(1,1,:,:)), sum(d2a(2,1,:,:)+d2b(2,1,:,:)), &
            sum(d2a(3,1,:,:)+d2b(3,1,:,:))
      deallocate(fka,fkb,d2w,d2a,d2b)
      stop
   endif
end block

block
   use mod_integrals, only: integrals_hessian_2e_df_j
   use mod_hessian_freq, only: hessian_write_full
   character(32) :: envdfj
   integer :: envstatdfj, uhdfj
   real(8),allocatable :: Ptdfj(:,:), d2Edfj(:,:,:,:)
   call get_environment_variable("ENGINE_HESSIAN_DFJ_DUMP", envdfj, status=envstatdfj)
   if (envstatdfj == 0 .and. iconv == 1) then
      allocate(Ptdfj(Nconts,Nconts), d2Edfj(3,Natoms,3,Natoms))
      Ptdfj = Pa + Pb
      block
         use mod_integrals, only: hessian_dfj_read_density
         character(256) :: dfile
         integer :: dstat
         call get_environment_variable("ENGINE_DFJ_DENSITY_FILE", dfile, status=dstat)
         if (dstat == 0) then
            call hessian_dfj_read_density(Nconts, trim(dfile), Ptdfj)
            print *, 'DFJ: read trial density from ', trim(dfile)
         endif
      end block
      call integrals_hessian_2e_df_j(Nconts, Natoms, Ptdfj, d2Edfj)
      open(newunit=uhdfj, file='direwolf_dfj.hess', status='replace', action='write')
      write(uhdfj,'(A)') ''
      write(uhdfj,'(A)') '$direwolf_hessian_file'
      call hessian_write_full(uhdfj, Natoms, d2Edfj)
      close(uhdfj)
      print *, 'RI-J direct Hessian written to direwolf_dfj.hess'
      stop
   endif
end block

block
   use mod_xc_hess, only: xc_hess_direct
   use mod_hessian_freq, only: hessian_write_full
   use mod_integrals, only: hessian_dfj_read_density
   character(32) :: envxch, envxcdump
   character(256) :: dfilexch
   integer :: envstatxch, envstatxcd, dstatxch, uhxch, ixch, iych, iach, ibch
   real(8),allocatable :: Ptxch(:,:), d2Exch(:,:,:,:)
   real(8) :: colsumxch, xcmax
   call get_environment_variable("ENGINE_XC_HESS_DEBUG", envxch, status=envstatxch)
   if (envstatxch == 0 .and. iconv == 1) then
      allocate(Ptxch(Nconts,Nconts), d2Exch(3,Natoms,3,Natoms))
      Ptxch = Pa + Pb
      call get_environment_variable("ENGINE_DFJ_DENSITY_FILE", dfilexch, status=dstatxch)
      if (dstatxch == 0) then
         call hessian_dfj_read_density(Nconts, trim(dfilexch), Ptxch)
         print *, 'XC HESS: read trial density from ', trim(dfilexch)
      endif
      call xc_hess_direct(Nconts, Natoms, Ptxch, d2Exch)
      xcmax = 0.0d0
      do iach = 1,Natoms
         do ixch = 1,3
            do iych = 1,3
               colsumxch = 0.0d0
               do ibch = 1,Natoms
                  colsumxch = colsumxch + d2Exch(ixch,iach,iych,ibch)
               enddo
               xcmax = max(xcmax, abs(colsumxch))
            enddo
         enddo
      enddo
      print *, 'XC HESS: max|x-sum_B| = ', xcmax, ' max|d2E_xc| = ', maxval(abs(d2Exch))
      call get_environment_variable("ENGINE_XC_HESS_DUMP", envxcdump, status=envstatxcd)
      if (envstatxcd == 0) then
         open(newunit=uhxch, file='direwolf_xc.hess', status='replace', action='write')
         write(uhxch,'(A)') ''
         write(uhxch,'(A)') '$direwolf_hessian_file'
         call hessian_write_full(uhxch, Natoms, d2Exch)
         close(uhxch)
         print *, 'fixed-density XC Hessian written to direwolf_xc.hess'
      endif
      stop
   endif
end block

block
   use mod_xc_hess, only: xc_fxc_apply
   use mod_integrals, only: hessian_dfj_read_density
   use GRID_info, only: Grids, ngrids
   character(32) :: envxf
   character(256) :: dfxf
   integer :: envstatf, dstatf, uf, mf, nf, gfi
   real(8),allocatable :: Ptf(:,:), Dtf(:,:,:), Vtf(:,:,:)
   call get_environment_variable("ENGINE_XC_FXC_TEST", envxf, status=envstatf)
   if (envstatf == 0 .and. iconv == 1) then
      allocate(Ptf(Nconts,Nconts), Dtf(Nconts,Nconts,2), Vtf(Nconts,Nconts,2))
      Ptf = Pa + Pb
      call get_environment_variable("ENGINE_DFJ_DENSITY_FILE", dfxf, status=dstatf)
      if (dstatf == 0) call hessian_dfj_read_density(Nconts, trim(dfxf), Ptf)
      Dtf = 0.0d0
      do mf = 1,Nconts
         Dtf(mf,mf,1) = 1.0d0
      enddo
      do mf = 1,Nconts
         do nf = 1,mf
            Dtf(mf,nf,2) = 0.01d0*dble(mf-nf)
            Dtf(nf,mf,2) = Dtf(mf,nf,2)
         enddo
      enddo
      Vtf = 0.0d0
      call xc_fxc_apply(Nconts, 2, Ptf, Dtf, Vtf)
      open(newunit=uf, file='direwolf_xc_fxc.txt', status='replace', action='write')
      write(uf,*) ngrids, Nconts
      do gfi = 1,ngrids
         write(uf,'(4ES24.16)') Grids(gfi)%coor, Grids(gfi)%weight
      enddo
      do mf = 1,Nconts
         write(uf,'(*(ES24.16))') Dtf(mf,:,1)
      enddo
      do mf = 1,Nconts
         write(uf,'(*(ES24.16))') Dtf(mf,:,2)
      enddo
      do mf = 1,Nconts
         write(uf,'(*(ES24.16))') Vtf(mf,:,1)
      enddo
      do mf = 1,Nconts
         write(uf,'(*(ES24.16))') Vtf(mf,:,2)
      enddo
      close(uf)
      print *, 'XC FXC test written to direwolf_xc_fxc.txt'
      stop
   endif
end block

block
   use mod_xc_hess, only: xc_vxc_deriv_atom
   use mod_integrals, only: hessian_dfj_read_density
   use GRID_info, only: Grids, ngrids
   character(32) :: envxv
   character(256) :: dfxv
   integer :: envstatv, dstatv, uv, mv, nv, gvi, iav, icv
   real(8),allocatable :: Ptv(:,:), dVv(:,:)
   call get_environment_variable("ENGINE_XC_VXCDERIV_TEST", envxv, status=envstatv)
   if (envstatv == 0 .and. iconv == 1) then
      allocate(Ptv(Nconts,Nconts), dVv(Nconts,Nconts))
      Ptv = Pa + Pb
      call get_environment_variable("ENGINE_DFJ_DENSITY_FILE", dfxv, status=dstatv)
      if (dstatv == 0) call hessian_dfj_read_density(Nconts, trim(dfxv), Ptv)
      open(newunit=uv, file='direwolf_xc_vxcderiv.txt', status='replace', action='write')
      write(uv,*) ngrids, Nconts, Natoms
      do gvi = 1,ngrids
         write(uv,'(4ES24.16)') Grids(gvi)%coor, Grids(gvi)%weight
      enddo
      do mv = 1,Nconts
         write(uv,'(*(ES24.16))') Ptv(mv,:)
      enddo
      do iav = 1,Natoms
         do icv = 1,3
            dVv = 0.0d0
            call xc_vxc_deriv_atom(Nconts, Natoms, Ptv, iav, icv, dVv)
            do mv = 1,Nconts
               write(uv,'(*(ES24.16))') dVv(mv,:)
            enddo
         enddo
      enddo
      close(uv)
      print *, 'XC Vxc-deriv test written to direwolf_xc_vxcderiv.txt'
      stop
   endif
end block

block
   use mod_integrals, only: integrals_hessian_2e_df_j
   character(32) :: envtr
   integer :: envstattr, itr_i, itr_j, iatr, ibtr, icm, jcm
   real(8),allocatable :: Ptr(:,:), d2Etr(:,:,:,:)
   real(8) :: str, seedtr, hmaxtr, rowmaxtr
   call get_environment_variable("ENGINE_HESSIAN_DFJ_TRANSCHECK", envtr, status=envstattr)
   if (envstattr == 0 .and. iconv == 1) then
      allocate(Ptr(Nconts,Nconts), d2Etr(3,Natoms,3,Natoms))
      seedtr = 0.12345d0
      do itr_i = 1,Nconts
         do itr_j = 1,itr_i
            seedtr = mod(seedtr*16807.0d0, 2147483647.0d0)
            if (abs(seedtr) < 1.0d-30) seedtr = 12345.0d0
            str = seedtr/2147483647.0d0 - 0.5d0
            Ptr(itr_i,itr_j) = str
            Ptr(itr_j,itr_i) = str
         enddo
      enddo
      call integrals_hessian_2e_df_j(Nconts, Natoms, Ptr, d2Etr)
      hmaxtr = maxval(abs(d2Etr))
      rowmaxtr = 0.0d0
      do iatr = 1,Natoms
         do icm = 1,3
            do jcm = 1,3
               str = 0.0d0
               do ibtr = 1,Natoms
                  str = str + d2Etr(icm,iatr,jcm,ibtr)
               enddo
               if (abs(str) > rowmaxtr) then
                  rowmaxtr = abs(str)
                  print *, 'DFJ TRANSCHECK worst row: atom',iatr,' comp',icm,'/',jcm, &
                           ' sum=',str
               endif
            enddo
         enddo
      enddo
      print *, 'DFJ TRANSCHECK: max|sum_B H[(A,x),(B,y)]| =', rowmaxtr
      print *, 'DFJ TRANSCHECK: max|H| =', hmaxtr, '  ratio =', rowmaxtr/max(hmaxtr,1d-300)
      stop
   endif
end block

block
   use mod_integrals, only: integrals_hessian_debug_check1e
   character(32) :: envval1e
   integer :: envstat1e, i1e, j1e
   real(8),allocatable :: Ptot1e(:,:), Wtot1e(:,:)
   call get_environment_variable("ENGINE_HESSIAN_DEBUG1E", envval1e, status=envstat1e)
   if (envstat1e == 0 .and. iconv == 1) then
      allocate(Ptot1e(Nconts,Nconts), Wtot1e(Nconts,Nconts))
      Ptot1e = Pa+Pb
      Wtot1e = 0.0d0
      do i1e = 1,Nconts
         do j1e = 1,Nconts
            block
               integer :: k1e
               do k1e = 1,n_alpha
                  Wtot1e(i1e,j1e) = Wtot1e(i1e,j1e) + 2.0d0*eLev_a(k1e)*C_a(i1e,k1e)*C_a(j1e,k1e)
               enddo
            end block
         enddo
      enddo
      call integrals_hessian_debug_check1e(Nconts, Natoms, Ptot1e, Wtot1e)
      stop
   endif
end block

block
   use mod_cpscf, only: cphf_debug_check_residual, cphf_debug_dump_dP, cphf_debug_grad_check, &
                        cphf_debug_hessian_row, cphf_debug_hessian_row_uhf, cphf_debug_e1_fixed, &
                        cphf_debug_ccbk_raw, cphf_debug_dF_raw, &
                        cphf_debug_amat_batch_check, cphf_debug_batch_scaling, &
                        cpscf_debug_amat_cosx_symmetry, cphf_solve_efield
   character(32) :: envval
   character(len=32) :: cphf_envs(8)
   character(len=256) :: cphf_eval
   integer :: envstat, iatom, ic, i, j, ie
   logical :: any_cphf_env
   cphf_envs = [ character(len=32) :: &
      "ENGINE_CPHF_DEBUG", "ENGINE_CPHF_AMAT_BATCH_CHECK", &
      "ENGINE_CPHF_BATCH_SCALING", "ENGINE_CPHF_DUMP_DP", &
      "ENGINE_CPHF_GRAD_CHECK", "ENGINE_DF_RAW", "ENGINE_CCBK_RAW", &
      "ENGINE_CPSCF_AMAT_COSX_SYM" ]
   any_cphf_env = .false.
   do ie = 1, size(cphf_envs)
      call get_environment_variable(trim(cphf_envs(ie)), cphf_eval, status=envstat)
      if (envstat == 0) any_cphf_env = .true.
   enddo
   if (any_cphf_env .and. (n_alpha /= n_beta)) then
      print *, 'ERROR: this CPSCF debug hook is closed-shell (RHF) only.'
      print *, '       n_alpha =', n_alpha, ', n_beta =', n_beta, ', multiplicity =', Multi
      print *, '       Refusing to return an alpha-only response.'
      stop 1
   endif
   call get_environment_variable("ENGINE_CPHF_DEBUG", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call cphf_debug_check_residual(Nconts, Natoms, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
      stop
   endif
   call get_environment_variable("ENGINE_CPHF_AMAT_BATCH_CHECK", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call cphf_debug_amat_batch_check(Nconts, n_alpha, C_a, HF_exchange_frac)
      stop
   endif
   call get_environment_variable("ENGINE_CPSCF_AMAT_COSX_SYM", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call cpscf_debug_amat_cosx_symmetry(Nconts, n_alpha, C_a, HF_exchange_frac)
      stop
   endif
   call get_environment_variable("ENGINE_CPHF_BATCH_SCALING", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call cphf_debug_batch_scaling(Nconts, n_alpha, C_a)
      stop
   endif
   call get_environment_variable("ENGINE_CPHF_DUMP_DP", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      read(envval,*) iatom, ic
      call cphf_debug_dump_dP(Nconts, Natoms, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
      stop
   endif
   call get_environment_variable("ENGINE_CPHF_GRAD_CHECK", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call cphf_debug_grad_check(Nconts, Natoms, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
      stop
   endif
   if (opts%nmr .and. iconv == 1) then
      block
         use mod_cosmo_state, only: cosmo_enabled
         use mod_nmr, only: nmr_shielding, nmr_isotropy_anisotropy
         use mod_nics, only: nics_setup
         use mod_molden, only: elem_symbol
         use mod_json_out, only: json_nmr, json_jcoupling, json_nics
         use mod_engine_report, only: report_energy_only
         use omp_lib, only: omp_get_wtime
         real(8),allocatable :: sdia(:,:,:), spara(:,:,:), cn(:,:)
         character(len=2),allocatable :: sym(:)
         real(8) :: tot(3,3), iso, ani, tn0, pvn(3)
         integer :: ia_n, ib_n, niter_n, kk_n
         integer :: np_n, nring_n, zn_n(Natoms)
         real(8),allocatable :: pos_n(:,:), nrm_n(:,:), sdp_n(:,:,:), spp_n(:,:,:)
         character(len=40),allocatable :: lab_n(:)
         integer,allocatable :: ring_n(:), ratm_n(:,:), rsz_n(:)
         logical :: ok_n
         character(len=160) :: msg_n
         tn0 = omp_get_wtime()
         allocate(sdia(3,3,Natoms), spara(3,3,Natoms), cn(3,Natoms), sym(Natoms))
         do ia_n = 1,Natoms
            cn(:,ia_n) = atoms(ia_n)%coor*ans2bohr
            sym(ia_n) = elem_symbol(atoms(ia_n)%charge)
            zn_n(ia_n) = atoms(ia_n)%charge
         enddo
         np_n = 0
         if (len_trim(opts%nics) > 0) then
            call nics_setup(opts%nics, Natoms, zn_n, cn, np_n, pos_n, lab_n, nrm_n, ring_n, ratm_n, rsz_n, nring_n, ok_n, msg_n)
            if (.not. ok_n) then
               print '(A)', 'Input error: '//trim(msg_n)
               stop 1
            endif
         else
            allocate(pos_n(3,1))
         endif
         allocate(sdp_n(3,3,max(np_n,1)), spp_n(3,3,max(np_n,1)))
         if (n_alpha /= n_beta .or. Multi /= 1) then
            call nmr_shielding(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, Natoms, cn, sdia, spara, niter_n, &
                               HF_exchange_frac < 1.0d0, engine_use_df_j, nocc_b=n_beta, C_b=C_b, eLev_b=eLev_b, &
                               nprobe=np_n, probe_bohr=pos_n, sdia_p=sdp_n, spara_p=spp_n)
         else
            call nmr_shielding(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, Natoms, cn, sdia, spara, niter_n, &
                               HF_exchange_frac < 1.0d0, engine_use_df_j, &
                               nprobe=np_n, probe_bohr=pos_n, sdia_p=sdp_n, spara_p=spp_n)
         endif
         print '(A)', '[NMR]'
         print '(A)', 'Magnetic shielding tensors (ppm), GIAO; sigma_ab = d2E/dB_a dm_b'
         print '(A)', ''
         print '(A)', '  atom  elem   isotropic   anisotropy   diamagnetic   paramagnetic'
         do ia_n = 1,Natoms
            tot = sdia(:,:,ia_n) + spara(:,:,ia_n)
            call nmr_isotropy_anisotropy(tot, iso, ani, pvn)
            print '(I6,2X,A4,4F13.3)', ia_n, sym(ia_n), iso, ani, (sdia(1,1,ia_n)+sdia(2,2,ia_n)+sdia(3,3,ia_n))/3.0d0, &
                  (spara(1,1,ia_n)+spara(2,2,ia_n)+spara(3,3,ia_n))/3.0d0
         enddo
         if (engine_verbose >= 2) then
            do ia_n = 1,Natoms
               tot = sdia(:,:,ia_n) + spara(:,:,ia_n)
               print '(A,I0,A,A)', 'Atom ', ia_n, ' ', sym(ia_n)
               do kk_n = 1,3
                  print '(A,3F14.4)', '   ', (tot(kk_n,ib_n), ib_n = 1,3)
               enddo
            enddo
         endif
         print '(A,I0)', 'Coupled-perturbed iterations: ', niter_n
         print '(A,F10.1)', 'Wall time (s): ', omp_get_wtime() - tn0
         call report_energy_only(E)
         call json_nmr(Natoms, sym, sdia, spara, niter_n)
         print '(A)', '[NMREND]'
         call flush(6)
         if (np_n > 0) then
            block
               real(8) :: tp(3,3), niso(np_n), nzz(np_n)
               integer :: ip, ir
               do ip = 1,np_n
                  tp = sdp_n(:,:,ip) + spp_n(:,:,ip)
                  niso(ip) = -(tp(1,1) + tp(2,2) + tp(3,3))/3.0d0
                  nzz(ip) = -sum(matmul(nrm_n(:,ip), tp)*nrm_n(:,ip))
               enddo
               print '(A)', '[NICS]'
               print '(A)', 'Nucleus-independent chemical shifts, NICS = -sigma (ppm); zz along the ring normal (lab z for a point)'
               print '(A)', ''
               if (nring_n > 0) then
                  print '(A)', '  ring  size  atoms                      NICS(0)  NICS(0)zz   NICS(+1)   NICS(-1)  '// &
                        'NICS(1)avg  NICS(1)zz'
                  do ir = 1,nring_n
                     block
                        character(len=60) :: al
                        integer :: ja
                        al = ''
                        do ja = 1,rsz_n(ir)
                           write(al(len_trim(al)+merge(1,2,len_trim(al) == 0):), '(I0)') ratm_n(ja,ir)
                        enddo
                        print '(I6,I6,2X,A24,6F11.3)', ir, rsz_n(ir), al(1:24), niso(3*ir-2), nzz(3*ir-2), &
                              niso(3*ir-1), niso(3*ir), &
                              0.5d0*(niso(3*ir-1) + niso(3*ir)), 0.5d0*(nzz(3*ir-1) + nzz(3*ir))
                     end block
                  enddo
               else
                  print '(A)', '  point       x         y         z (Angstrom)       NICS      NICS_zz'
                  do ip = 1,np_n
                     print '(I7,3F10.4,2F12.3)', ip, pos_n(:,ip)/ans2bohr, niso(ip), nzz(ip)
                  enddo
               endif
               print '(A)', '[NICSEND]'
               call json_nics(np_n, lab_n, pos_n/ans2bohr, niso, nzz)
            end block
            call flush(6)
         endif
         if (opts%nmr_j .and. n_alpha == n_beta .and. Multi == 1) then
            block
               use mod_jcoupling, only: jcoupling_fc, jcoupling_pso_sd, jcoupling_dso, jcoupling_to_hz_all, gyromagnetic_ratio
               real(8),allocatable :: Sj(:,:), Spso(:,:), Ssd(:,:), Dso(:,:), Jf(:,:), Js(:,:), Jp(:,:), Jd(:,:)
               real(8),allocatable :: Jt(:,:), Kr(:,:)
               integer :: nj, nj2, zj(Natoms)
               logical :: w_fc, w_sd, w_pso, w_dso
               character(len=40) :: tl
               character(len=64) :: al, tok
               logical, allocatable :: insel(:)
               integer, allocatable :: sel(:)
               real(8), allocatable :: cs_(:,:), Ss_(:,:), Sp_(:,:), Sd_(:,:)
               integer :: ns, ic, ir, i1, i2, ios
               tl = ',' // trim(opts%nmr_j_terms) // ','
               w_fc = trim(opts%nmr_j_terms) == 'all' .or. index(tl, ',fc,') > 0
               w_sd = trim(opts%nmr_j_terms) == 'all' .or. index(tl, ',sd,') > 0
               w_pso = trim(opts%nmr_j_terms) == 'all' .or. index(tl, ',pso,') > 0
               w_dso = trim(opts%nmr_j_terms) == 'all' .or. index(tl, ',dso,') > 0
               allocate(Sj(Natoms,Natoms), Spso(Natoms,Natoms), Ssd(Natoms,Natoms), Dso(Natoms,Natoms))
               allocate(Jf(Natoms,Natoms), Js(Natoms,Natoms), Jp(Natoms,Natoms), Jd(Natoms,Natoms))
               allocate(Jt(Natoms,Natoms), Kr(Natoms,Natoms))
               Sj = 0.0d0 ; Spso = 0.0d0 ; Ssd = 0.0d0 ; Dso = 0.0d0 ; nj = 0 ; nj2 = 0
               do ia_n = 1,Natoms
                  zj(ia_n) = atoms(ia_n)%charge
               enddo
               allocate(insel(Natoms))
               insel = (trim(opts%nmr_j_atoms) == 'all')
               if (trim(opts%nmr_j_atoms) /= 'all') then
                  al = trim(opts%nmr_j_atoms)
                  do while (len_trim(al) > 0)
                     ic = index(al, ',')
                     if (ic == 0) ic = len_trim(al) + 1
                     tok = al(1:ic-1)
                     al = al(min(ic+1, len(al)):)
                     if (len_trim(tok) == 0) cycle
                     ir = index(tok, '-')
                     if (ir > 0) then
                        read(tok(1:ir-1), *, iostat=ios) i1
                        if (ios == 0) read(tok(ir+1:), *, iostat=ios) i2
                     else
                        read(tok, *, iostat=ios) i1
                        i2 = i1
                     endif
                     if (ios /= 0 .or. i1 < 1 .or. i2 > Natoms .or. i1 > i2) &
                        error stop "jcoupling_atoms: 'all' or atom numbers and ranges such as 1,3-5"
                     insel(i1:i2) = .true.
                  enddo
               endif
               do ia_n = 1,Natoms
                  if (gyromagnetic_ratio(zj(ia_n)) == 0.0d0) insel(ia_n) = .false.
               enddo
               ns = count(insel)
               allocate(sel(ns))
               ns = 0
               do ia_n = 1,Natoms
                  if (insel(ia_n)) then
                     ns = ns + 1
                     sel(ns) = ia_n
                  endif
               enddo
               if (ns > 0) then
                  allocate(cs_(3,ns), Ss_(ns,ns))
                  do ia_n = 1,ns
                     cs_(:,ia_n) = cn(:,sel(ia_n))
                  enddo
                  if (w_fc) then
                     call jcoupling_fc(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, ns, cs_, Ss_, nj)
                     Sj(sel,sel) = Ss_
                  endif
                  if (w_pso .or. w_sd) then
                     allocate(Sp_(ns,ns), Sd_(ns,ns))
                     call jcoupling_pso_sd(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, ns, cs_, &
                                           w_pso, w_sd, Sp_, Sd_, nj2)
                     Spso(sel,sel) = Sp_ ; Ssd(sel,sel) = Sd_
                  endif
                  if (w_dso) then
                     call jcoupling_dso(Nconts, Pa, Pb, ns, cs_, Ss_)
                     Dso(sel,sel) = Ss_
                  endif
               endif
               call jcoupling_to_hz_all(Natoms, zj, Sj, Ssd, Spso, Dso, Jf, Js, Jp, Jd, Jt, Kr)
               print '(A)', '[JCOUPLING]'
               print '(A)', 'Indirect spin-spin coupling, Ramsey terms: FC Fermi contact, SD spin-dipole,'
               print '(A)', 'PSO paramagnetic spin-orbit, DSO diamagnetic spin-orbit'
               print '(A)', 'Terms included: '//trim(opts%nmr_j_terms)
               print '(A,I0,A)', 'Nuclei taking part: ', ns, ' (jcoupling_atoms = '//trim(opts%nmr_j_atoms)//')'
               print '(A)', 'Isotopes: 1H 13C 15N 17O 19F 29Si 31P 33S 35Cl and 7Li; J in Hz, K reduced, 1e19 T^2/J'
               print '(A,I0)', 'Response iterations (triplet, singlet): ', max(nj, 0)
               print '(A)', ''
               print '(A)', '  atom  atom   elem  elem   J_total  K_reduced        FC        SD       PSO       DSO'
               do ia_n = 1,Natoms
                  do ib_n = ia_n+1,Natoms
                     if (Jt(ia_n,ib_n) /= 0.0d0 .or. Jf(ia_n,ib_n) /= 0.0d0) &
                        print '(2I6,2X,A4,2X,A4,2F11.3,4F10.3)', ia_n, ib_n, sym(ia_n), sym(ib_n), Jt(ia_n,ib_n), Kr(ia_n,ib_n), &
                              Jf(ia_n,ib_n), Js(ia_n,ib_n), Jp(ia_n,ib_n), Jd(ia_n,ib_n)
                  enddo
               enddo
               call json_jcoupling(Natoms, Jt, Kr, Jf, Js, Jp, Jd, opts%nmr_j_terms)
               print '(A)', '[JCOUPLINGEND]'
            end block
         endif
      end block
   endif
   if (opts%polar .and. iconv == 1) then
      block
         use mod_cosmo_state, only: cosmo_enabled
         use mod_engine_report, only: report_polar, report_polar_dynamic, report_energy_only
         use mod_engine_input_types, only: parse_polar_omega
         use mod_cosmo_state, only: cosmo_enabled, f_eps, cosmo_cpcm
         use mod_cosmo_hessian, only: cosmo_hess_active
         use omp_lib, only: omp_get_wtime
         real(8) :: alpha(3,3), aiso, aani
         real(8), parameter :: AU2ANG3 = 0.148184711d0
         integer :: niter_p, ip, jp, iom, nom
         real(8) :: om(16), ad(3,3), adiso, adani
         logical :: omok
         real(8) :: tp0, f_eps_static_p
         character(len=1), parameter :: ax(3) = (/ 'x', 'y', 'z' /)
         tp0 = omp_get_wtime()
         if (Multi /= 1) error stop 'polarizability: open-shell references are not implemented'
         f_eps_static_p = f_eps
         if (cosmo_enabled) cosmo_hess_active = .true.
         if (RS_omega > 0.0d0 .and. .not. cosx_enabled) then
            block
               use mod_tddft, only: tddft_polar_hybrid
               call tddft_polar_hybrid(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, 0.0d0, alpha, niter_p)
            end block
         else
            call cphf_solve_efield(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, alpha, niter_p)
         endif
         aiso = (alpha(1,1) + alpha(2,2) + alpha(3,3))/3.0d0
         aani = sqrt(0.5d0*((alpha(1,1)-alpha(2,2))**2 + (alpha(2,2)-alpha(3,3))**2 + (alpha(3,3)-alpha(1,1))**2 &
                    + 6.0d0*(alpha(1,2)**2 + alpha(1,3)**2 + alpha(2,3)**2)))
         print '(A)', '[POLARIZABILITY]'
         print '(A)', 'Static dipole polarizability (a.u.), alpha_ij = -d2E/dF_i dF_j'
         print '(A,3A14)', '       ', ('             '//ax(jp), jp = 1,3)
         do ip = 1,3
            print '(A,A,3F14.4)', '    ', ax(ip), (alpha(ip,jp), jp = 1,3)
         enddo
         print '(A,F12.4,A,F12.4,A)', 'Isotropic   (xx+yy+zz)/3 = ', aiso, ' a.u. = ', aiso*AU2ANG3, ' A^3'
         print '(A,F12.4,A,F12.4,A)', 'Anisotropy               = ', aani, ' a.u. = ', aani*AU2ANG3, ' A^3'
         if (cosmo_enabled) print '(A)', 'Solvent: equilibrium (static eps) for the static tensor'
         print '(A,I0)', 'Coupled-perturbed iterations: ', niter_p
         print '(A,F10.1)', 'Wall time (s): ', omp_get_wtime() - tp0
         call report_energy_only(E)
         call report_polar(alpha, aiso, aani, niter_p)
         if (cosmo_enabled) then
            if (cosmo_cpcm) then
               f_eps = (opts%tddft_eps_opt-1.0d0)/opts%tddft_eps_opt
            else
               f_eps = (opts%tddft_eps_opt-1.0d0)/(opts%tddft_eps_opt+0.5d0)
            endif
         endif
         nom = 0
         if (len_trim(opts%polar_omega) > 0) call parse_polar_omega(trim(opts%polar_omega), om, nom, omok)
         do iom = 1, nom
            tp0 = omp_get_wtime()
            if (HF_exchange_frac /= 0.0d0 .or. RS_omega > 0.0d0) then
               block
                  use mod_tddft, only: tddft_polar_hybrid
                  call tddft_polar_hybrid(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, om(iom), ad, niter_p)
               end block
            else
               call cphf_solve_efield(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, ad, niter_p, om(iom))
            endif
            adiso = (ad(1,1) + ad(2,2) + ad(3,3))/3.0d0
            adani = sqrt(0.5d0*((ad(1,1)-ad(2,2))**2 + (ad(2,2)-ad(3,3))**2 + (ad(3,3)-ad(1,1))**2 &
                       + 6.0d0*(ad(1,2)**2 + ad(1,3)**2 + ad(2,3)**2)))
            print '(A)', ''
            print '(A,F10.6,A,F9.4,A,F10.2,A)', 'Dynamic dipole polarizability (a.u.) at omega = ', om(iom), ' a.u. = ', &
                  om(iom)*27.211386245988d0, ' eV = ', merge(45.56335252760d0/om(iom), 0.0d0, om(iom) > 0.0d0), ' nm'
            print '(A,3A14)', '       ', ('             '//ax(jp), jp = 1,3)
            do ip = 1,3
               print '(A,A,3F14.4)', '    ', ax(ip), (ad(ip,jp), jp = 1,3)
            enddo
            print '(A,F12.4,A,F12.4,A)', 'Isotropic   (xx+yy+zz)/3 = ', adiso, ' a.u. = ', adiso*AU2ANG3, ' A^3'
            print '(A,F12.4,A,F12.4,A)', 'Anisotropy               = ', adani, ' a.u. = ', adani*AU2ANG3, ' A^3'
            print '(A,I0)', 'Coupled-perturbed iterations: ', niter_p
            print '(A,F10.1)', 'Wall time (s): ', omp_get_wtime() - tp0
            call report_polar_dynamic(om(iom), ad, adiso, adani, niter_p)
         enddo
         if (cosmo_enabled) then
            if (nom > 0) print '(A,F7.3)', 'Solvent: dynamic tensors with the optical eps ', opts%tddft_eps_opt
            f_eps = f_eps_static_p
            cosmo_hess_active = .false.
         endif
         if (xc_uses_vv10()) then
            print '(A)', 'Note: the VV10 nonlocal kernel is not in the response '// &
                         '(0.04% of alpha for water, not quantified elsewhere)'
         endif
         print '(A)', '[POLARIZABILITYEND]'
         call flush(6)
      end block
   endif
   if (opts%tddft_nroots > 0 .and. iconv == 1) then
      block
         use mod_cosmo_state, only: cosmo_enabled, f_eps, cosmo_cpcm
         use mod_cosmo_hessian, only: cosmo_hess_active
         use mod_tddft, only: tddft_rks, tddft_main_transitions
         use mod_tddft_density, only: tddft_hole_electron, tddft_hole_electron_uks, tddft_he_wanted
         use mod_tddft_nto, only: tddft_nto_rks, tddft_nto_uks
         use mod_engine_report, only: report_tddft_state, report_energy_only, RPT_TD_SINGLET, RPT_TD_TRIPLET, RPT_TD_UKS
         use mod_cpscf, only: cpscf_triplet
         use mod_tddft, only: tddft_uks, tddft_tda, tddft_note
         real(8) :: f_eps_static
         logical :: tdp
         real(8), parameter :: HARTREE_EV = 27.211386245988d0, HARTREE_NM = 45.56335252760d0
         real(8), allocatable :: wex(:), fosc(:), tdm(:,:), xpy(:,:,:), xmy(:,:,:), ctm(:,:), ctv(:)
         integer :: nr, nv, napp, ir, ii, nlist, ispin, nspin
         integer, allocatable :: tfrom(:), tto(:)
         real(8), allocatable :: twt(:)
         real(8) :: tdconv
         tdconv = 1.0d-4
         f_eps_static = f_eps
         if (cosmo_enabled) then
            if (cosmo_cpcm) then
               f_eps = (opts%tddft_eps_opt-1.0d0)/opts%tddft_eps_opt
            else
               f_eps = (opts%tddft_eps_opt-1.0d0)/(opts%tddft_eps_opt+0.5d0)
            endif
            cosmo_hess_active = .true.
         endif
         tddft_tda = opts%tddft_tda
         tdp = opts%tddft_print
         if (Multi /= 1) then
            block
               real(8), allocatable :: uw(:), uf(:), ut(:,:), upy(:,:), umy(:,:), uct(:,:), uctv(:)
               integer :: nva, nvb, nau, ntu, s, nj
               character(len=1) :: sl
               integer, allocatable :: jf(:), jt(:)
               real(8), allocatable :: jw(:)
               character(len=1), allocatable :: js(:)
               if (tdp .and. (opts%tddft_triplets)) &
                  print '(A)', 'Note: tddft_triplets applies to closed-shell references only; '// &
                  'the states of an open-shell reference are all in the one list'
               nva = Nconts - n_alpha ; nvb = Nconts - n_beta
               nau = nva*n_alpha ; ntu = nau + nvb*n_beta
               nr = min(opts%tddft_nroots, ntu)
               allocate(uw(nr), uf(nr), ut(3,nr), upy(ntu,nr), umy(ntu,nr))
               allocate(tfrom(ntu), tto(ntu), twt(ntu))
               call report_energy_only(E)
               call tddft_uks(Nconts, n_alpha, n_beta, C_a, eLev_a, C_b, eLev_b, Pa, Pb, HF_exchange_frac, nr, &
                              uw, uf, ut, upy, umy, napp, tdconv, opts%tddft_robust)
               if (tddft_he_wanted(opts%tddft_cube)) then
                  allocate(uct(5,nr))
                  call tddft_hole_electron_uks(Nconts, n_alpha, n_beta, C_a, C_b, nr, upy, umy, trim(opts%tddft_cube), &
                                               trim(opts%freq_base), uct)
               endif
               if (tdp) print '(A)', '[TDDFT]'
               if (tdp) print '(A)', 'Excited states (UKS linear response, adiabatic kernel)'
               if (tdp) print '(A)', 'Method: '//trim(tddft_note)
               if (tdp .and. (cosmo_enabled)) print '(A,F7.3)', 'Solvent: non-equilibrium, optical eps ', opts%tddft_eps_opt
               if (tdp) print '(A)', ''
               if (tdp) print '(A)', 'state      E (eV)    E (Ha)       nm        f        mu_x      mu_y      mu_z'
               do ir = 1,nr
                  if (tdp) print '(I5,F12.4,F10.4,F10.2,F9.4,3F10.4)', ir, uw(ir)*HARTREE_EV, uw(ir), HARTREE_NM/uw(ir), &
                        uf(ir), ut(:,ir)
               enddo
               if (tdp) print '(A)', ''
               if (tdp) &
                  print '(A,A,I0,A,I0,A)', 'Main transitions, weight X^2 - Y^2 >= 0.05 (orbitals numbered from 1 per spin, ', &
                     'HOMO = ', n_alpha, 'a / ', n_beta, 'b)'
               allocate(jf(2*ntu), jt(2*ntu), jw(2*ntu), js(2*ntu))
               do ir = 1,nr
                  if (tdp) print '(A,I0,A,F9.4,A)', 'state ', ir, '  (', uw(ir)*HARTREE_EV, ' eV)'
                  nj = 0
                  do s = 1,2
                     if (s == 1) then
                        call tddft_main_transitions(nva, n_alpha, upy(1:nau,ir), umy(1:nau,ir), 0.05d0, &
                                                    nlist, tfrom, tto, twt)
                        sl = 'a'
                     else
                        if (n_beta == 0) cycle
                        call tddft_main_transitions(nvb, n_beta, upy(nau+1:ntu,ir), umy(nau+1:ntu,ir), 0.05d0, &
                                                    nlist, tfrom, tto, twt)
                        sl = 'b'
                     endif
                     do ii = 1,nlist
                        if (tdp) print '(A,I6,A,A,I6,A,F10.4)', '  ', tfrom(ii), sl, ' ->', tto(ii), sl, twt(ii)
                        nj = nj + 1
                        jf(nj) = tfrom(ii) ; jt(nj) = tto(ii) ; jw(nj) = twt(ii) ; js(nj) = sl
                     enddo
                  enddo
                  if (allocated(uct)) then
                     uctv = uct(:,ir)
                  else if (allocated(uctv)) then
                     deallocate(uctv)
                  endif
                     call report_tddft_state(RPT_TD_UKS, ir, uw(ir), uf(ir), ut(:,ir), nj, jf, jt, jw, napp, &
                                             trim(tddft_note), spins=js, ct=uctv)
               enddo
               deallocate(jf, jt, jw, js)
               if (allocated(uct)) then
                  if (tdp) print '(A)', ''
                  if (tdp) print '(A)', 'Hole and electron (unrelaxed densities of the response vectors, both spins, on a grid)'
                  if (tdp) print '(A)', 'state   D_CT (A)   q_CT (e)  sigma_h (A)  sigma_e (A)   Lambda'
                  do ir = 1,nr
                     if (tdp) print '(I5,F11.3,F11.3,F13.3,F13.3,F9.3)', ir, uct(:,ir)
                  enddo
               endif
               if (tdp) print '(A,I0)', 'Response vectors applied (Davidson): ', napp
               if (xc_uses_vv10()) then
                  if (tdp) print '(A)', 'Note: the VV10 nonlocal kernel is not in the response (the ground state includes VV10)'
               endif
               if (tdp) print '(A)', '[TDDFTEND]'
               call flush(6)
               if (index(opts%tddft_cube, 'nto') > 0) &
                  call tddft_nto_uks(Nconts, n_alpha, n_beta, C_a, C_b, nr, upy, trim(opts%tddft_cube), trim(opts%freq_base))
               deallocate(uw, uf, ut, upy, umy, tfrom, tto, twt)
               if (allocated(uct)) deallocate(uct)
               if (allocated(uctv)) deallocate(uctv)
            end block
         else
         nr = opts%tddft_nroots
         nv = Nconts - n_alpha
         if (nr > nv*n_alpha) nr = nv*n_alpha
         allocate(wex(nr), fosc(nr), tdm(3,nr), xpy(nv,n_alpha,nr), xmy(nv,n_alpha,nr))
         allocate(tfrom(nv*n_alpha), tto(nv*n_alpha), twt(nv*n_alpha))
         call report_energy_only(E)
         nspin = merge(2, 1, opts%tddft_triplets)
         do ispin = 1,nspin
            cpscf_triplet = (ispin == 2)
            call tddft_rks(Nconts, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, nr, wex, fosc, tdm, &
                           xpy, xmy, napp, tdconv, opts%tddft_robust)
            cpscf_triplet = .false.
            if (ispin == 2) then
               fosc = 0.0d0 ; tdm = 0.0d0
            else
               if (allocated(ctm)) deallocate(ctm)
               if (tddft_he_wanted(opts%tddft_cube)) then
                  allocate(ctm(5,nr))
                  call tddft_hole_electron(Nconts, n_alpha, C_a, nr, xpy, xmy, trim(opts%tddft_cube), trim(opts%freq_base), ctm)
               endif
            endif
            if (tdp) print '(A)', '[TDDFT]'
            if (ispin == 1) then
               if (tdp) print '(A)', 'Singlet excited states (RKS linear response, adiabatic kernel)'
            else
               if (tdp) print '(A)', 'Triplet excited states (RKS linear response, adiabatic kernel)'
            endif
            if (tdp) print '(A)', 'Method: '//trim(tddft_note)
            if (tdp .and. (cosmo_enabled)) print '(A,F7.3)', 'Solvent: non-equilibrium, optical eps ', opts%tddft_eps_opt
            if (tdp) print '(A)', ''
            if (tdp) print '(A)', 'state      E (eV)    E (Ha)       nm        f        mu_x      mu_y      mu_z'
            do ir = 1,nr
               if (tdp) print '(I5,F12.4,F10.4,F10.2,F9.4,3F10.4)', ir, wex(ir)*HARTREE_EV, wex(ir), HARTREE_NM/wex(ir), &
                     fosc(ir), tdm(:,ir)
            enddo
            if (tdp) print '(A)', ''
            if (tdp) print '(A,I0,A)', 'Main transitions, weight X^2 - Y^2 >= 0.05 (orbitals numbered from 1, HOMO = ', n_alpha, ')'
            do ir = 1,nr
               call tddft_main_transitions(nv, n_alpha, xpy(:,:,ir), xmy(:,:,ir), 0.05d0, nlist, tfrom, tto, twt)
               if (tdp) print '(A,I0,A,F9.4,A)', 'state ', ir, '  (', wex(ir)*HARTREE_EV, ' eV)'
               do ii = 1,nlist
                  if (tdp) print '(A,I6,A,I6,F10.4)', '  ', tfrom(ii), ' ->', tto(ii), twt(ii)
               enddo
               if (ispin == 1 .and. allocated(ctm)) then
                  ctv = ctm(:,ir)
               else if (allocated(ctv)) then
                  deallocate(ctv)
               endif
                  call report_tddft_state(merge(RPT_TD_TRIPLET, RPT_TD_SINGLET, ispin == 2), ir, wex(ir), fosc(ir), tdm(:,ir), &
                                          nlist, tfrom, tto, twt, napp, trim(tddft_note), ct=ctv)
            enddo
            if (ispin == 1 .and. allocated(ctm)) then
               if (tdp) print '(A)', ''
               if (tdp) print '(A)', 'Hole and electron (unrelaxed densities of the response vectors, on a grid)'
               if (tdp) print '(A)', 'state   D_CT (A)   q_CT (e)  sigma_h (A)  sigma_e (A)   Lambda'
               do ir = 1,nr
                  if (tdp) print '(I5,F11.3,F11.3,F13.3,F13.3,F9.3)', ir, ctm(:,ir)
               enddo
            endif
            if (tdp) print '(A,I0)', 'Response vectors applied (Davidson): ', napp
            if (xc_uses_vv10()) then
               if (tdp) print '(A)', 'Note: the VV10 nonlocal kernel is not in the response (the ground state includes VV10)'
            endif
            if (tdp) print '(A)', '[TDDFTEND]'
            call flush(6)
            if (ispin == 1 .and. index(opts%tddft_cube, 'nto') > 0) &
               call tddft_nto_rks(Nconts, n_alpha, C_a, nr, xpy, trim(opts%tddft_cube), trim(opts%freq_base))
         enddo
         deallocate(wex, fosc, tdm, xpy, xmy, tfrom, tto, twt)
         if (allocated(ctm)) deallocate(ctm)
         if (allocated(ctv)) deallocate(ctv)
         endif
         f_eps = f_eps_static
         cosmo_hess_active = .false.
      end block
   endif
   call get_environment_variable("ENGINE_HESSIAN_TOTAL_DEBUG", envval, status=envstat)
   if ((envstat == 0 .or. opts%freq) .and. iconv == 1) then
      iatom = 1; ic = 1
      if (envstat == 0) read(envval,*) iatom, ic
      block
         use mod_hessian_freq, only: freq_unit, freq_energy
         use mod_engine_report, only: report_energy_only
         freq_energy = E
         call report_energy_only(E)
         write(freq_unit,'(A,F20.10)') ' FINAL SINGLE POINT ENERGY (Eh): ', E
      end block
      if (engine_use_df_k .and. (HF_exchange_frac > 1.0d-12 .or. RS_omega > 0.0d0)) then
         print *, "ERROR: analytic Hessian is not implemented for K='RI' (RI-K): the exchange second"
         print *, "       derivative would need third-derivative 3-index integrals and a derivative of the"
         print *, "       fitting metric inverse. Use K='cosx' (RIJCOSX) or K='exact'."
         error stop "analytic Hessian: K='RI' not supported"
      endif
      if (xc_uses_vv10()) then
         print *, "WARNING: analytic Hessian omits the VV10 nonlocal term (second derivative of E_VV10);"
         print *, "         frequencies are off by up to ~14 cm-1 (methylamine wB97M-V/def2-SVP)."
      endif
      block
         use mod_hessian_freq, only: hess_classical
         real(8) :: cb(3,Natoms), cp(3,Natoms), cm(3,Natoms), hstep
         integer :: zc(Natoms), ja, jc, ka
         logical :: smd_on
         smd_on = cosmo_on .and. cosmo_smd
         do ka = 1,Natoms
            cb(:,ka) = atoms(ka)%coor*ans2bohr
            zc(ka) = atoms(ka)%charge
         enddo
         if (allocated(hess_classical)) deallocate(hess_classical)
         allocate(hess_classical(3,Natoms,3,Natoms))
         hess_classical = 0.0d0
         if (disp_s6 > 0.0d0 .or. use_chg_dispersion .or. len_trim(disp_d3_damping) > 0 .or. &
             len_trim(disp_d4_method) > 0 .or. len_trim(gcp_method) > 0 .or. smd_on) then
            hstep = 1.0d-3
            do ja = 1,Natoms
               do jc = 1,3
                  cp = cb; cp(jc,ja) = cp(jc,ja) + hstep
                  cm = cb; cm(jc,ja) = cm(jc,ja) - hstep
                  block
                     real(8) :: gp(3,Natoms), gm(3,Natoms)
                     call classical_gradient_total(Natoms, zc, cp, icharge, disp_s6, use_chg_dispersion, &
                        trim(disp_d3_damping), trim(Functional), trim(disp_d4_method), trim(gcp_method), &
                        smd_on, smd_n_val, smd_alpha_val, smd_beta_val, smd_gamma_val, smd_phi_val, &
                        smd_psi_val, smd_is_water, gp)
                     call classical_gradient_total(Natoms, zc, cm, icharge, disp_s6, use_chg_dispersion, &
                        trim(disp_d3_damping), trim(Functional), trim(disp_d4_method), trim(gcp_method), &
                        smd_on, smd_n_val, smd_alpha_val, smd_beta_val, smd_gamma_val, smd_phi_val, &
                        smd_psi_val, smd_is_water, gm)
                     hess_classical(:,:,jc,ja) = (gp - gm)/(2.0d0*hstep)
                  end block
               enddo
            enddo
            hess_classical = 0.5d0*(hess_classical + reshape(hess_classical,(/3,Natoms,3,Natoms/),order=(/3,4,1,2/)))
            print *, 'classical (dispersion/gCP/CDS) Hessian added, max|H| =', maxval(abs(hess_classical))
            block
               use mod_hessian_freq, only: hessian_write_full
               integer :: uc
               open(newunit=uc, file='direwolf_classical.hess', status='replace', action='write')
               write(uc,'(A)') ''
               write(uc,'(A)') '$direwolf_hessian_file'
               call hessian_write_full(uc, Natoms, hess_classical)
               close(uc)
            end block
         endif
      end block
      Pa = matmul(C_a(:,1:n_alpha), transpose(C_a(:,1:n_alpha)))
      if (n_alpha /= n_beta .or. Multi /= 1) then
         Pb = matmul(C_b(:,1:n_beta), transpose(C_b(:,1:n_beta)))
      else
         Pb = Pa
      endif
      call get_environment_variable("ENGINE_HESSIAN_FORCE_UHF", envval, status=envstat)
      if (n_alpha /= n_beta .or. Multi /= 1) then
         call cphf_debug_hessian_row_uhf(Nconts, Natoms, n_alpha, n_beta, C_a, eLev_a, C_b, eLev_b, &
                                         Pa, Pb, HF_exchange_frac, iatom, ic &
                                         )
      else if (envstat == 0 .and. trim(envval) == "1") then
         call cphf_debug_hessian_row_uhf(Nconts, Natoms, n_alpha, n_alpha, C_a, eLev_a, C_a, eLev_a, &
                                         Pa, Pb, HF_exchange_frac, iatom, ic &
                                         )
      else
         call cphf_debug_hessian_row(Nconts, Natoms, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic &
                                     )
      endif
      if (.not. opts%hessian_return) then
         opts%done = .true.
         energy_out = E
         force_out = 0.0d0 ; MLcharge_out = 0.0d0
         mem_grid_gb = 0.0d0 ; mem_2e_gb = 0.0d0
         if (.not. allocated(dens_out_a)) allocate(dens_out_a(nConts,nConts))
         if (.not. allocated(dens_out_b)) allocate(dens_out_b(nConts,nConts))
         dens_out_a = Pa ; dens_out_b = Pb
         opts%rep = rpt
         return
      endif
   endif
   call get_environment_variable("ENGINE_DF_RAW", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      block
         integer :: iatc2, icc2
         read(envval,*) iatc2, icc2
         call cphf_debug_dF_raw(Nconts, Natoms, iatc2, icc2, Pa, Pb, HF_exchange_frac)
      end block
      stop
   endif
   call get_environment_variable("ENGINE_CCBK_RAW", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      block
         integer :: iatc, imu, inu
         read(envval,*) iatc, imu, inu
         call cphf_debug_ccbk_raw(Nconts, Natoms, iatc, imu, inu)
      end block
      stop
   endif
   call get_environment_variable("ENGINE_DUMP_S", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      block
         integer :: i1s,j1s
         print *, 'S DUMP:'
         do i1s = 1,Nconts
            print *, (S(i1s,j1s), j1s=1,Nconts)
         enddo
      end block
      stop
   endif
   call get_environment_variable("ENGINE_DUMP_E1FIXED", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call cphf_debug_e1_fixed(Nconts, Pa, Pb)
      stop
   endif
   call get_environment_variable("ENGINE_DUMP_PTOT", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      print *, 'PTOT DUMP:'
      do i = 1,Nconts
         print *, (Pa(i,j)+Pb(i,j), j=1,Nconts)
      enddo
      stop
   endif
end block

if (do_force) then
   if (.not. engine_quiet) print *, '[FORCE]'
   if (engine_use_df) call decide_df_force_mode()
   if (do_force) then
      force_dense = force_dense_allow
      block
         character(len=8) :: nfd
         call get_environment_variable("ENGINE_NO_FORCE_DENSE", nfd)
         if (len_trim(nfd) > 0) force_dense = .false.
      end block
      force_dense_mgga = force_dense .and. xc_uses_tau()
      if (force_dense) then
         call xcgrid_free_derived()
         call gridgen(info)
         call GTOeval(info)
      endif
      block
         real(8) :: excf
         real(8),allocatable :: fxca(:,:), fxcb(:,:)
         allocate(fxca(nconts,nconts), fxcb(nconts,nconts))
         call DFT_calc(excf, fxca, fxcb, info, .true.)
      end block
   endif
   call prof_start("force")
   call calc_force(info)
   call prof_stop("force")
   if (.not. engine_quiet) then
      print *, '[FORCEEND]'
      print *
   endif
endif

if (iconv .eq. 1 .and. .not. engine_quiet) then
   if (resp_charges_on) then
      block
      real(8) :: RESPcharge_out(Natoms)
      call resp_charges(Natoms, RESPcharge_out)
      call calc_properties(Natoms, MLcharge_out, RESPcharge_out)
      end block
   else
      call calc_properties(Natoms, MLcharge_out)
   endif
endif
if (iconv .eq. 1 .and. opts%esp_charges) then
   if (allocated(opts%esp_q)) deallocate(opts%esp_q)
   allocate(opts%esp_q(Natoms))
   call resp_charges(Natoms, opts%esp_q)
endif
if (iconv .eq. 1 .and. opts%pc_grad_request .and. nPointCharges > 0) then
   block
   real(8), parameter :: hpc = 2.0d-3
   real(8), allocatable :: pts(:,:), Vp(:), Ptot_pc(:,:)
   real(8) :: rp(3), dist_pc
   integer :: ic_pc, k_pc, ia_pc
   allocate(pts(3,6*nPointCharges), Vp(6*nPointCharges), Ptot_pc(nConts,nConts))
   Ptot_pc = Pa + Pb
   do ic_pc = 1, nPointCharges
      rp = pointcharge_coor(:,ic_pc) * ans2bohr
      do k_pc = 1, 3
         pts(:,6*(ic_pc-1)+2*k_pc-1) = rp ; pts(k_pc,6*(ic_pc-1)+2*k_pc-1) = rp(k_pc) + hpc
         pts(:,6*(ic_pc-1)+2*k_pc)   = rp ; pts(k_pc,6*(ic_pc-1)+2*k_pc)   = rp(k_pc) - hpc
      enddo
   enddo
   call esp_at_grid_batch(6*nPointCharges, pts, Ptot_pc, Vp)
   do ic_pc = 1, 6*nPointCharges
      do ia_pc = 1, Natoms
         dist_pc = sqrt(sum((pts(:,ic_pc) - atoms(ia_pc)%coor*ans2bohr)**2))
         Vp(ic_pc) = Vp(ic_pc) + real(atoms(ia_pc)%charge - atoms(ia_pc)%ecpCoreElec,8)/dist_pc
      enddo
   enddo
   if (allocated(opts%pc_grad)) deallocate(opts%pc_grad)
   allocate(opts%pc_grad(nPointCharges,3))
   do ic_pc = 1, nPointCharges
      do k_pc = 1, 3
         opts%pc_grad(ic_pc,k_pc) = pointcharge_q(ic_pc) * (Vp(6*(ic_pc-1)+2*k_pc-1) - Vp(6*(ic_pc-1)+2*k_pc)) / (2.0d0*hpc)
      enddo
   enddo
   end block
endif
if (iconv .eq. 1 .and. opts%mayer) then
   call calc_mayer_bond_orders(Natoms)
   if (allocated(opts%mayer_bo)) deallocate(opts%mayer_bo)
   allocate(opts%mayer_bo(Natoms,Natoms))
   opts%mayer_bo = mayer_bo
endif
if (iconv .eq. 1 .and. cosmo_enabled) call cosmo_report_sigma_profile(cosmo_sigma_rav, cosmo_sigma_profile_file)
if (iconv .eq. 1 .and. cosmo_enabled) then
   block
   real(8) :: coor_ang_cosmors(3,Natoms)
   integer :: ia_cosmors
   do ia_cosmors = 1,Natoms
      coor_ang_cosmors(:,ia_cosmors) = atoms(ia_cosmors)%coor
   enddo
   call cosmo_write_dot_cosmo_file(Natoms, coor_ang_cosmors, trim(Functional), trim(baselable), E)
   end block
endif

if (disp_s6 .gt. 0.0d0) then
   block
      real(8) :: Edisp, dEdisp(3,Natoms)
      real(8) :: coor_bohr(3,Natoms)
      integer :: znum(Natoms), ia
      do ia = 1,Natoms
         coor_bohr(:,ia) = atoms(ia)%coor*ans2bohr
         znum(ia) = atoms(ia)%charge
      enddo
      call dispersion_d2(Natoms, znum, coor_bohr, disp_s6, Edisp, dEdisp)
      E = E + Edisp
      do ia = 1,Natoms
         atoms(ia)%atmForce = atoms(ia)%atmForce + dEdisp(:,ia)
      enddo
   end block
endif

if (use_chg_dispersion) then
   block
      real(8) :: Edisp, dEdisp(3,Natoms)
      real(8) :: coor_bohr(3,Natoms)
      integer :: znum(Natoms), ia
      do ia = 1,Natoms
         coor_bohr(:,ia) = atoms(ia)%coor*ans2bohr
         znum(ia) = atoms(ia)%charge
      enddo
      call dispersion_chg(Natoms, znum, coor_bohr, Edisp, dEdisp)
      E = E + Edisp
      do ia = 1,Natoms
         atoms(ia)%atmForce = atoms(ia)%atmForce + dEdisp(:,ia)
      enddo
   end block
endif

if (len_trim(disp_d3_damping) .gt. 0) then
   block
      real(8) :: Edisp, dEdisp(3,Natoms)
      real(8) :: coor_bohr(3,Natoms)
      integer :: znum(Natoms), ia
      do ia = 1,Natoms
         coor_bohr(:,ia) = atoms(ia)%coor*ans2bohr
         znum(ia) = atoms(ia)%charge
      enddo
      call dispersion_d3(Natoms, znum, coor_bohr, trim(Functional), &
                          trim(disp_d3_damping), Edisp, dEdisp)
      E = E + Edisp
      do ia = 1,Natoms
         atoms(ia)%atmForce = atoms(ia)%atmForce + dEdisp(:,ia)
      enddo
   end block
endif

if (len_trim(disp_d4_method) .gt. 0) then
   block
      real(8) :: Edisp, dEdisp(3,Natoms)
      real(8) :: coor_bohr(3,Natoms)
      integer :: znum(Natoms), ia
      do ia = 1,Natoms
         coor_bohr(:,ia) = atoms(ia)%coor*ans2bohr
         znum(ia) = atoms(ia)%charge
      enddo
      call dispersion_d4(Natoms, znum, coor_bohr, icharge, trim(disp_d4_method), Edisp, dEdisp)
      E = E + Edisp
      do ia = 1,Natoms
         atoms(ia)%atmForce = atoms(ia)%atmForce + dEdisp(:,ia)
      enddo
   end block
endif

if (len_trim(gcp_method) .gt. 0) then
   block
      real(8) :: Egcp, dEgcp(3,Natoms)
      real(8) :: coor_bohr(3,Natoms)
      integer :: znum(Natoms), ia
      do ia = 1,Natoms
         coor_bohr(:,ia) = atoms(ia)%coor*ans2bohr
         znum(ia) = atoms(ia)%charge
      enddo
      call gcp_correction(Natoms, znum, coor_bohr, trim(gcp_method), Egcp, dEgcp)
      E = E + Egcp
      do ia = 1,Natoms
         atoms(ia)%atmForce = atoms(ia)%atmForce + dEgcp(:,ia)
      enddo
   end block
endif

if (cosmo_on .and. cosmo_smd) then
   block
      real(8) :: E_cds, force_cds(3,Natoms)
      real(8) :: coor_bohr(3,Natoms)
      integer :: znum(Natoms), ia
      do ia = 1,Natoms
         coor_bohr(:,ia) = atoms(ia)%coor*ans2bohr
         znum(ia) = atoms(ia)%charge
      enddo
      call smd_cds_energy_force(Natoms, znum, coor_bohr, &
           smd_n_val, smd_alpha_val, smd_beta_val, smd_gamma_val, smd_phi_val, smd_psi_val, smd_is_water, &
           E_cds, force_cds)
      E = E + E_cds
      do ia = 1,Natoms
         atoms(ia)%atmForce = atoms(ia)%atmForce + force_cds(:,ia)
      enddo
   end block
endif

energy_out = E
do i = 1,Natoms
   force_out(i,:) = atoms(i)%atmForce
enddo
if (iconv .ne. 1) MLcharge_out = 0.0d0
opts%rep = rpt

end subroutine EngineUp

subroutine classical_gradient_total(natoms, znum, coor_bohr, icharge, s6, use_chg, d3damp, functional, &
                                    d4method, gcpmethod, smd_on, smd_n, smd_a, smd_b, smd_g, smd_p, smd_s, &
                                    smd_water, grad)
use mod_dispersion, only: dispersion_d2, dispersion_chg
use mod_dispersion_d3, only: dispersion_d3
use mod_dispersion_d4, only: dispersion_d4
use mod_gcp, only: gcp_correction
use mod_smd_cds, only: smd_cds_energy_force
implicit none
integer,intent(in) :: natoms, znum(natoms), icharge
real(8),intent(in) :: coor_bohr(3,natoms), s6, smd_n, smd_a, smd_b, smd_g, smd_p, smd_s
logical,intent(in) :: use_chg, smd_on, smd_water
character(*),intent(in) :: d3damp, functional, d4method, gcpmethod
real(8),intent(out) :: grad(3,natoms)
real(8) :: e, g(3,natoms)
grad = 0.0d0
if (s6 > 0.0d0) then
   call dispersion_d2(natoms, znum, coor_bohr, s6, e, g); grad = grad + g
endif
if (use_chg) then
   call dispersion_chg(natoms, znum, coor_bohr, e, g); grad = grad + g
endif
if (len(d3damp) > 0) then
   call dispersion_d3(natoms, znum, coor_bohr, functional, d3damp, e, g); grad = grad + g
endif
if (len(d4method) > 0) then
   call dispersion_d4(natoms, znum, coor_bohr, icharge, d4method, e, g); grad = grad + g
endif
if (len(gcpmethod) > 0) then
   call gcp_correction(natoms, znum, coor_bohr, gcpmethod, e, g); grad = grad + g
endif
if (smd_on) then
   call smd_cds_energy_force(natoms, znum, coor_bohr, smd_n, smd_a, smd_b, smd_g, smd_p, smd_s, smd_water, e, g)
   grad = grad + g
endif
end subroutine classical_gradient_total

subroutine reset_engine_state()
use MOL_info
use GRID_info
use mod_integrals, only: integrals_finalize, cosx_needs_rebuild
use mod_cosmo, only: cosmo_finalize
use mod_vv10, only: vv10_nlc_grid_reset
implicit none
call cosmo_finalize()
call vv10_nlc_grid_reset()
if (allocated(atoms))   deallocate(atoms)
if (allocated(Fa))      deallocate(Fa)
if (allocated(Fb))      deallocate(Fb)
if (allocated(S))       deallocate(S)
if (allocated(Hcore))   deallocate(Hcore)
if (allocated(Pa))      deallocate(Pa)
if (allocated(Pb))      deallocate(Pb)
if (allocated(C_a))     deallocate(C_a)
if (allocated(C_b))     deallocate(C_b)
if (allocated(eLev_a))  deallocate(eLev_a)
if (allocated(eLev_b))  deallocate(eLev_b)
if (allocated(X))       deallocate(X)
if (allocated(linkMat)) deallocate(linkMat)
if (allocated(Grids))   deallocate(Grids)
call xcgrid_free_derived()
if (allocated(pointcharge_q)) deallocate(pointcharge_q)
if (allocated(pointcharge_coor)) deallocate(pointcharge_coor)
grid_cache_capacity = 0
Natoms = 0
Nconts = 0
ncontssph = 0
nPointCharges = 0
E = 0
E_rep = 0
call integrals_finalize()
cosx_needs_rebuild = .true.
force_dense = .false.
force_dense_mgga = .false.
end subroutine reset_engine_state

integer function ecp_core_electrons_prescan(base_label)
use mod_basis_files, only: basis_file_for_label
implicit none
character(len=*),intent(in) :: base_label
character(len=40) :: ecp_label, tok
character(len=200) :: bpath
integer :: ios
logical :: found

ecp_core_electrons_prescan = 0
ecp_label = trim(base_label)//".ecp"
found = .false.
bpath = basis_file_for_label(ecp_label)
if (len_trim(bpath) == 0) return
open(unit=203, file=trim(bpath))
do
   read(203,*,iostat=ios) tok
   if (ios /= 0) exit
   if (trim(tok) == trim(ecp_label)) then
      found = .true.
      exit
   endif
enddo
if (found) read(203,*) ecp_core_electrons_prescan
close(203)
end function ecp_core_electrons_prescan

character(len=30) function resolve_basis_name(name_in)
implicit none
character(len=*),intent(in) :: name_in
character(len=30) :: u
integer :: i, ic

u = adjustl(name_in)
do i = 1, len_trim(u)
   ic = ichar(u(i:i))
   if (ic >= ichar('a') .and. ic <= ichar('z')) u(i:i) = char(ic - 32)
enddo

select case (trim(u))
case ('3-21G');                      resolve_basis_name = '321g'
case ('6-31G');                      resolve_basis_name = '6-31g'
case ('6-31G*','6-31G(D)');          resolve_basis_name = '631gs'
case ('6-31G**','6-31G(D,P)');       resolve_basis_name = '631gdp'
case ('6-31+G*','6-31+G(D)');        resolve_basis_name = '631pgs'
case ('6-31+G**','6-31+G(D,P)');     resolve_basis_name = '631pgss'
case ('6-311G');                     resolve_basis_name = '6311g'
case ('6-311G*','6-311G(D)');        resolve_basis_name = '6311gs'
case ('6-311G**','6-311G(D,P)');     resolve_basis_name = '6311gdp'
case ('6-311+G**','6-311+G(D,P)');   resolve_basis_name = '6311pgss'
case ('6-311++G**','6-311++G(D,P)'); resolve_basis_name = '6311ppgss'
case ('6-311+G(2D,P)');              resolve_basis_name = '6311pg2dp'
case ('6-311++G(2D,P)');             resolve_basis_name = '6311ppg2dp'
case ('CC-PVDZ');                    resolve_basis_name = 'ccpvdz'
case ('CC-PVTZ');                    resolve_basis_name = 'ccpvtz'
case ('CC-PVQZ');                    resolve_basis_name = 'ccpvqz'
case ('AUG-CC-PVDZ');                resolve_basis_name = 'augccpvdz'
case ('AUG-CC-PVTZ');                resolve_basis_name = 'augccpvtz'
case ('DEF2-SVP');                   resolve_basis_name = 'def2svp'
case ('DEF2-TZVP');                  resolve_basis_name = 'def2tzvp'
case ('DEF2-TZVPP');                 resolve_basis_name = 'def2tzvpp'
case ('DEF2-QZVP');                  resolve_basis_name = 'def2qzvp'
case ('DEF2-TZVPD');                 resolve_basis_name = 'def2tzvpd'
case ('SVP');                        resolve_basis_name = 'def2svp'
case ('TZVP');                       resolve_basis_name = 'def2tzvp'
case ('TZVPP');                      resolve_basis_name = 'def2tzvpp'
case ('QZVP');                       resolve_basis_name = 'def2qzvp'
case ('LANL2DZ');                    resolve_basis_name = 'lanl2dz'
case ('LANL08');                     resolve_basis_name = 'lanl08'
case ('CC-PVDZ-JKFIT');              resolve_basis_name = 'ccpvdzjkfit'
case ('CC-PVTZ-JKFIT');              resolve_basis_name = 'ccpvtzjkfit'
case ('CC-PVQZ-JKFIT');              resolve_basis_name = 'ccpvqzjkfit'
case ('DEF2-UNIVERSAL-JKFIT');       resolve_basis_name = 'def2universaljkfit'
case default;                        resolve_basis_name = name_in
end select
end function resolve_basis_name

character(len=30) function auto_ri_aux_basis(resolved_base, need_k_fit)
implicit none
character(len=*),intent(in) :: resolved_base
logical,intent(in) :: need_k_fit
select case (trim(resolved_base))
case ('ccpvdz','augccpvdz'); auto_ri_aux_basis = 'ccpvdzjkfit'
case ('ccpvtz','augccpvtz'); auto_ri_aux_basis = 'ccpvtzjkfit'
case ('ccpvqz');             auto_ri_aux_basis = 'ccpvqzjkfit'
case ('mtzvpp')
   if (need_k_fit) then
      auto_ri_aux_basis = 'def2universaljkfit'
   else
      auto_ri_aux_basis = 'def2mtzvpprij'
   endif
case default
   if (need_k_fit) then
      auto_ri_aux_basis = 'def2universaljkfit'
   else
      auto_ri_aux_basis = 'def2universaljfit'
   endif
end select
end function auto_ri_aux_basis

logical function aux_basis_covers_all_atoms(auxname, natoms_in, z_in)
use mod_basis_files, only: basis_file_for_label
implicit none
character(len=*),intent(in) :: auxname
integer,intent(in) :: natoms_in, z_in(natoms_in)
character(len=40) :: lbl
integer :: i

aux_basis_covers_all_atoms = .true.
do i = 1, natoms_in
   write(lbl,"(I0.2,A1,A)") z_in(i), '-', trim(auxname)
   if (len_trim(basis_file_for_label(trim(lbl))) == 0) then
      aux_basis_covers_all_atoms = .false.
      return
   endif
enddo
end function aux_basis_covers_all_atoms

subroutine engine_startup_mark(label)
use omp_lib, only: omp_get_wtime
use MOL_info, only: engine_verbose
implicit none
character(len=*),intent(in) :: label
real(8),save :: t0 = -1.0d0, tl = 0.0d0
character(len=8) :: e
real(8) :: t
t = omp_get_wtime()
if (t0 .lt. 0.0d0) then
   t0 = t ; tl = t
endif
if (engine_verbose .ge. 2) then
   block
      use mod_meminfo, only: get_process_memory_bytes
      integer(8) :: rss, hwm
      call get_process_memory_bytes(rss, hwm)
      print '(A,A30,2F9.3,A,F7.2,A,F7.2,A)', " STARTUP ", label, t - t0, t - tl, "   RSS ", real(rss,8)/1024.0d0**3, &
            " GB, peak ", real(hwm,8)/1024.0d0**3, " GB"
   end block
   call flush(6)
else if (t - tl .gt. 60.0d0) then
   print '(A,A,F8.1,A)', " [startup] step '", trim(label), t - tl, " s'"
   call flush(6)
endif
tl = t
end subroutine engine_startup_mark

subroutine smaps_malloc_info(fname)
use iso_c_binding
implicit none
character(len=*),intent(in) :: fname
interface
   function c_fopen(path, mode) bind(C, name="fopen") result(fp)
   import :: c_ptr, c_char
   character(kind=c_char),intent(in) :: path(*), mode(*)
   type(c_ptr) :: fp
   end function c_fopen
   function c_fclose(fp) bind(C, name="fclose") result(r)
   import :: c_ptr, c_int
   type(c_ptr),value :: fp
   integer(c_int) :: r
   end function c_fclose
   function c_malloc_info(opt, fp) bind(C, name="malloc_info") result(r)
   import :: c_ptr, c_int
   integer(c_int),value :: opt
   type(c_ptr),value :: fp
   integer(c_int) :: r
   end function c_malloc_info
   function c_malloc_trim(pad) bind(C, name="malloc_trim") result(r)
   import :: c_int, c_size_t
   integer(c_size_t),value :: pad
   integer(c_int) :: r
   end function c_malloc_trim
end interface
type(c_ptr) :: fp
integer(c_int) :: r
fp = c_fopen(trim(fname)//c_null_char, "w"//c_null_char)
r = c_malloc_info(0_c_int, fp)
r = c_malloc_trim(0_c_size_t)
r = c_malloc_info(0_c_int, fp)
r = c_fclose(fp)
end subroutine smaps_malloc_info
