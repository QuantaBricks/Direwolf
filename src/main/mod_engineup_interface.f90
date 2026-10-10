! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Explicit interface for EngineUp (src/main/engine.f90), for use ONLY by

module mod_engineup_interface
implicit none

interface
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
   implicit none
   integer       :: ncenters,imult,icharge,atomchg(ncenters)
   character     :: baselable*30,functional_in*30
   character     :: ecplabel*30
   character     :: atombasis(ncenters)*30, atomecp(ncenters)*30
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
   end subroutine EngineUp
end interface

end module mod_engineup_interface
