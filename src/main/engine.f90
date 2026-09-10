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
                    dens_in_a,dens_in_b,dens_out_a,dens_out_b)

use MOL_info
use GRID_info, only: Grids
use GRID_info, only: xcgrid_dynamic, xcgrid_refined, xcgrid_level, &
                     force_dense_allow, &
                     XCGRID_RSCALE, XCGRID_SPH_IN, XCGRID_SPH_EDGE, &
                     xcgrid_switch_prms, XCGRID_COARSE, XCGRID_FINE, &
                     XCGRID_L4, XCGRID_L5, XCGRID_L6, XCGRID_L7, &
                     xcgrid_active_coarse, xcgrid_active_fine, &
                     xc_direct_mode, XC_DIRECT_NCONTS, force_dense, force_dense_mgga, &
                     xcgrid_fine_level_for_functional, &
                     grid_cache_mode, max_batch_nsig, grid_cache_capacity, ngrids
use mod_integrals, only: resp_charges, &
                          integrals_init, engine_use_df, engine_use_df_j, engine_use_df_k, &
                          engine_df_aux_basis, &
                          engine_estimate_only, integrals_estimate_aux_size, &
                          nRec, nContsAux, decide_df_force_mode, engine_puream, &
                          engine_harris_guess, calc_properties, engine_do_force
use mod_meminfo, only: engine_mem_cap_bytes, engine_avail_at_start_bytes, get_memory_budget_bytes, get_process_memory_bytes
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
use mod_molden, only: write_molden, read_molden
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
logical :: chk_ok
logical :: has_warm_start

real(8)   :: dist

real(8) :: smd_n_val, smd_alpha_val, smd_beta_val, smd_gamma_val, smd_phi_val, smd_psi_val
logical :: smd_is_water

if (n_threads .gt. 0) call omp_set_num_threads(n_threads)

call openblas_set_num_threads(1)

print *,"System cores available:",omp_get_num_procs(), &
        " Direwolf OpenMP threads in use:",omp_get_max_threads()
call flush(6)

call reset_engine_state()
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

block
use mod_basis_files, only: set_basis_dir_override
call set_basis_dir_override(basedir)
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
  endif
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

print '("n_alpha = ", I0, "   n_beta = ", I0, "   (# occupied alpha/beta electrons)")', n_alpha, n_beta

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

call prof_start("basis_and_integrals")
call integrals_init(info)
call prof_stop("basis_and_integrals")

if (engine_use_df) then
   call integrals_estimate_aux_size()
   mem_2e_gb = real(nconts,8)*real(nconts,8)*real(nContsAux,8)*8.0d0/1024.0d0**3
else
   mem_2e_gb = real(nRec,8)*8.0d0/1024.0d0**3
endif

if (HF_exchange_frac .lt. 1.0d0) then
   print *, '[GRID]'
   block
      block
         character(len=32) :: dyngridenv, cacheenv
         call get_environment_variable("ENGINE_GRID_CACHE", cacheenv)
         if (len_trim(cacheenv) .gt. 0) then
            xc_direct_mode = (trim(cacheenv) .eq. "0")
         else
            xc_direct_mode = (nconts .gt. XC_DIRECT_NCONTS)
            if (xc_direct_mode) print *,"DIRECT XC by default: nconts=",nconts, &
                 " >",XC_DIRECT_NCONTS," (override with ENGINE_GRID_CACHE=1)"
         endif
         xcgrid_active_coarse = XCGRID_COARSE
         xcgrid_active_fine   = xcgrid_fine_level_for_functional(Functional)
         if (xcgrid_active_fine .eq. XCGRID_L5) xcgrid_active_coarse = XCGRID_L4
         call get_environment_variable("ENGINE_XC_GRID_LEVEL", dyngridenv)
         if (len_trim(dyngridenv) .gt. 0) then
            block
              integer :: reqlvl
              read(dyngridenv,*) reqlvl
              select case (reqlvl)
              case (2:7)
                 xcgrid_active_fine   = reqlvl
                 xcgrid_active_coarse = max(XCGRID_COARSE, reqlvl - 1)
              case default
                 print *,"ENGINE_XC_GRID_LEVEL: must be 2-7, ignoring"
              end select
            end block
            print *,"ENGINE_XC_GRID_LEVEL override: active pair = (", &
                    xcgrid_active_coarse,",",xcgrid_active_fine,")"
         endif
         block
           logical :: split_override
           split_override = .false.
           call get_environment_variable("ENGINE_XC_GRID_LEVEL_COARSE", dyngridenv)
           if (len_trim(dyngridenv) .gt. 0) then
              block
                integer :: reqlvl
                read(dyngridenv,*) reqlvl
                if (reqlvl .ge. 2 .and. reqlvl .le. 7) then
                   xcgrid_active_coarse = reqlvl
                   split_override = .true.
                else
                   print *,"ENGINE_XC_GRID_LEVEL_COARSE: must be 2-7, ignoring"
                endif
              end block
           endif
           call get_environment_variable("ENGINE_XC_GRID_LEVEL_FINE", dyngridenv)
           if (len_trim(dyngridenv) .gt. 0) then
              block
                integer :: reqlvl
                read(dyngridenv,*) reqlvl
                if (reqlvl .ge. 2 .and. reqlvl .le. 7) then
                   xcgrid_active_fine = reqlvl
                   split_override = .true.
                else
                   print *,"ENGINE_XC_GRID_LEVEL_FINE: must be 2-7, ignoring"
                endif
              end block
           endif
           if (split_override) &
              print *,"ENGINE_XC_GRID_LEVEL_COARSE/_FINE override: active pair = (", &
                      xcgrid_active_coarse,",",xcgrid_active_fine,")"
         end block
         call get_environment_variable("ENGINE_XC_DYNGRID", dyngridenv)
         if (len_trim(dyngridenv) .gt. 0) then
            xcgrid_dynamic = (trim(dyngridenv) .eq. "1")
         else
            xcgrid_dynamic = (nconts .gt. XC_DIRECT_NCONTS)
         endif
         has_warm_start = (size(dens_in_a,1) .gt. 0) .or. molden_read
         xcgrid_refined = .false.
         if (xcgrid_dynamic .and. .not. has_warm_start) then
            xcgrid_level = xcgrid_active_coarse
            call get_environment_variable("ENGINE_XC_DYNGRID_RSCALE", dyngridenv)
            if (len_trim(dyngridenv) .gt. 0) read(dyngridenv,*) XCGRID_RSCALE(xcgrid_active_coarse)
            call get_environment_variable("ENGINE_XC_DYNGRID_SPH", dyngridenv)
            if (len_trim(dyngridenv) .gt. 0) then
               read(dyngridenv,*) XCGRID_SPH_IN(xcgrid_active_coarse)
               select case (XCGRID_SPH_IN(xcgrid_active_coarse))
               case (974) ; XCGRID_SPH_EDGE(xcgrid_active_coarse) = 590
               case (434) ; XCGRID_SPH_EDGE(xcgrid_active_coarse) = 230
               case (302) ; XCGRID_SPH_EDGE(xcgrid_active_coarse) = 170
               case (230) ; XCGRID_SPH_EDGE(xcgrid_active_coarse) = 110
               case (170) ; XCGRID_SPH_EDGE(xcgrid_active_coarse) = 110
               case default
                  print *,"ENGINE_XC_DYNGRID_SPH: unsupported order, ignoring"
                  XCGRID_SPH_IN(xcgrid_active_coarse) = 302
               end select
            endif
            call get_environment_variable("ENGINE_XC_PROD_SPH", dyngridenv)
            if (len_trim(dyngridenv) .gt. 0) then
               read(dyngridenv,*) XCGRID_SPH_IN(xcgrid_active_fine)
               select case (XCGRID_SPH_IN(xcgrid_active_fine))
               case (974) ; XCGRID_SPH_EDGE(xcgrid_active_fine) = 590
               case (590) ; XCGRID_SPH_EDGE(xcgrid_active_fine) = 434
               case (434) ; XCGRID_SPH_EDGE(xcgrid_active_fine) = 230
               case (302) ; XCGRID_SPH_EDGE(xcgrid_active_fine) = 170
               case (230) ; XCGRID_SPH_EDGE(xcgrid_active_fine) = 110
               case (194) ; XCGRID_SPH_EDGE(xcgrid_active_fine) = 110
               case default
                  print *,"ENGINE_XC_PROD_SPH: unsupported order, ignoring"
                  XCGRID_SPH_IN(xcgrid_active_fine) = 434
               end select
            endif
            call get_environment_variable("ENGINE_XC_PROD_RSCALE", dyngridenv)
            if (len_trim(dyngridenv) .gt. 0) read(dyngridenv,*) XCGRID_RSCALE(xcgrid_active_fine)
            call get_environment_variable("ENGINE_XC_DYNGRID_PRMS", dyngridenv)
            if (len_trim(dyngridenv) .gt. 0) read(dyngridenv,*) xcgrid_switch_prms
            if (engine_verbose .ge. 2) &
               print '(A,F5.3,A,I0,A,ES9.2)',"  Two-stage XC grid: coarse rscale=", &
                       XCGRID_RSCALE(xcgrid_active_coarse)," sph=",XCGRID_SPH_IN(xcgrid_active_coarse), &
                       " promote at dP <",xcgrid_switch_prms
         else
            xcgrid_level = xcgrid_active_fine
            xcgrid_refined = .true.
         endif
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
      call gridgen(info)
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
   print *, '[GRIDEND]'
   print *
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

print *, '[SCF]'
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
call guess(info,Gtype,Pa_chk,Pb_chk)
deallocate(Pa_chk,Pb_chk)

select case (scf_conv_level)
case (1)
   Emax = 1.00E-08
   Pmax = 1.00E-07
case (2)
   Emax = 1.00E-09
   Pmax = 1.00E-08
case default
   Emax = 1.00E-06
   Pmax = 1.00E-05
end select
itmax =150
call prof_start("scf_total")
call SCFcycle(info,Emax,Pmax,itmax,iconv,econv,do_force)
call prof_stop("scf_total")

allocate(dens_out_a(nConts,nConts), dens_out_b(nConts,nConts))
dens_out_a = Pa
dens_out_b = Pb
if (molden_write .and. iconv .eq. 1) call write_molden(molden_file)

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
   use mod_cphf, only: cphf_debug_check_residual, cphf_debug_dump_dP, cphf_debug_grad_check, &
                        cphf_debug_hessian_row, cphf_debug_e1_fixed, cphf_debug_ccbk_raw, cphf_debug_dF_raw
   character(32) :: envval
   integer :: envstat, iatom, ic, i, j
   call get_environment_variable("ENGINE_CPHF_DEBUG", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      call cphf_debug_check_residual(Nconts, Natoms, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac)
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
   call get_environment_variable("ENGINE_HESSIAN_TOTAL_DEBUG", envval, status=envstat)
   if (envstat == 0 .and. iconv == 1) then
      read(envval,*) iatom, ic
      call cphf_debug_hessian_row(Nconts, Natoms, n_alpha, C_a, eLev_a, Pa, Pb, HF_exchange_frac, iatom, ic)
      stop
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
   print *, '[FORCE]'
   if (engine_use_df) call decide_df_force_mode()
   if (do_force) then
      force_dense = force_dense_allow
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
   print *, '[FORCEEND]'
   print *
endif

if (iconv .eq. 1) then
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

end subroutine EngineUp

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
