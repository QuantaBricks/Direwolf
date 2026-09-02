! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Dynamic two-stage XC grid: frees/refines the coarse-to-fine DFT quadrature grid mid-SCF.

subroutine xcgrid_free_derived()
use GRID_info
implicit none

if (allocated(val_blocks))     deallocate(val_blocks)
grid_cache_capacity = 0

if (allocated(TempD_all))      deallocate(TempD_all)
if (allocated(rcut2_shared))   deallocate(rcut2_shared)
if (allocated(batch_nsig))     deallocate(batch_nsig)
if (allocated(batch_sig_idx))  deallocate(batch_sig_idx)
n_cache_batches = 0
max_batch_nsig  = 0

if (allocated(dft_batch_shells))     deallocate(dft_batch_shells)
if (allocated(dft_batch_shells_lvl)) deallocate(dft_batch_shells_lvl)
n_dft_batches          = 0
max_dft_batch_nsig     = 0
max_dft_batch_nsig_loose = 0
dft_batch_shells_built = .false.
dft_loose_built        = .false.

if (allocated(vxc_sh_atom))   deallocate(vxc_sh_atom)
if (allocated(vxc_sh_local))  deallocate(vxc_sh_local)
if (allocated(vxc_sh_ao0))    deallocate(vxc_sh_ao0)
if (allocated(vxc_sh_ndim))   deallocate(vxc_sh_ndim)
if (allocated(batch_natom))   deallocate(batch_natom)
if (allocated(batch_atoms))   deallocate(batch_atoms)
if (allocated(batch_phimax))  deallocate(batch_phimax)
if (allocated(ao_atom))       deallocate(ao_atom)
if (allocated(sh_in_batch))   deallocate(sh_in_batch)
if (allocated(shell_batches)) deallocate(shell_batches)
vxc_nshell_tot     = 0
shell_tables_built = .false.

if (allocated(xcr_Pa))     deallocate(xcr_Pa)
if (allocated(xcr_rho))    deallocate(xcr_rho)
if (allocated(xcr_grad))   deallocate(xcr_grad)
if (allocated(xcr_Fxc))    deallocate(xcr_Fxc)
if (allocated(atpair_dp))  deallocate(atpair_dp)
if (allocated(atpair_ddp)) deallocate(atpair_ddp)
xcr_Exc         = 0.0d0
xc_incr_primed  = .false.

end subroutine xcgrid_free_derived

subroutine xcgrid_refine()
use GRID_info
use MOL_info, only: engine_verbose
use mod_meminfo, only: get_process_memory_bytes
use mod_profile, only: prof_start, prof_stop, fmt_gb
implicit none
integer :: info
integer(8) :: diag_rss, diag_hwm

if (.not. xcgrid_dynamic)  return
if (xcgrid_refined)        return
if (xcgrid_level .eq. xcgrid_active_fine) return

call xcgrid_free_derived()

xcgrid_level        = xcgrid_active_fine
xcgrid_refined      = .true.
xcgrid_just_refined = .true.

print *,"XC grid: promoting to production grid (coarse stage done)"
if (engine_verbose .ge. 2) then
   call get_process_memory_bytes(diag_rss, diag_hwm)
   print '(A)', "  DIAG before production grid_gen: RSS="//trim(fmt_gb(real(diag_rss,8)/1024**3))//"GB"
   call flush(6)
endif

call prof_start("grid_gen")
call gridgen(info)
call prof_stop("grid_gen")
call prof_start("gto_eval")
call GTOeval(info)
call prof_stop("gto_eval")
if (engine_verbose .ge. 2) then
   call get_process_memory_bytes(diag_rss, diag_hwm)
   print '(A)', "  DIAG after production grid_gen+GTOeval: RSS="//trim(fmt_gb(real(diag_rss,8)/1024**3))//"GB"
endif

end subroutine xcgrid_refine
