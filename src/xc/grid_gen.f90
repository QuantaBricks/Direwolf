! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Builds the DFT quadrature grid (Becke partitioning, radial/angular point generation, including NLC/VV10 variants).

subroutine load_or_generate_grid(info)
use MOL_info
use GRID_info
implicit none
integer :: info
character(len=512) :: gpath
logical :: ok
call get_environment_variable("ENGINE_XC_GRID_FILE", gpath)
if (len_trim(gpath) > 0) then
   call gridgen_from_file(trim(gpath), ok)
   if (ok) return
   print '(A)', " [GRIDFILE] could not read "//trim(gpath)//" - generating instead"
endif
call gridgen(info)
end subroutine load_or_generate_grid

subroutine gridgen_from_file(path, ok)
use MOL_info
use GRID_info
implicit none
character(len=*), intent(in) :: path
logical, intent(out) :: ok
integer :: u, ios, n, k
real(8) :: gx, gy, gz, gw
ok = .false.
open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
if (ios /= 0) return
n = 0
do
   read(u, *, iostat=ios) gx, gy, gz, gw
   if (ios /= 0) exit
   n = n + 1
enddo
if (n <= 0) then
   close(u)
   return
endif
rewind(u)
if (allocated(Grids)) deallocate(Grids)
allocate(Grids(n))
do k = 1, n
   read(u, *, iostat=ios) gx, gy, gz, gw
   if (ios /= 0) exit
   Grids(k)%coor(1) = gx
   Grids(k)%coor(2) = gy
   Grids(k)%coor(3) = gz
   Grids(k)%weight  = gw
enddo
close(u)
call spatial_sort_grid(n)
ngrids = n
print '(A,I0,A)', " [GRIDFILE] loaded ", n, " points from "//trim(path)
ok = .true.
end subroutine gridgen_from_file

subroutine gridgen(info)
use MOL_info
use GRID_info
    implicit none
INCLUDE 'parameter.h'
    interface
       subroutine gridgen_count_points(totGrid, grdRec)
       integer,intent(out) :: totGrid
       integer,allocatable,intent(out) :: grdRec(:,:)
       end subroutine gridgen_count_points
    end interface
    integer    :: i,j,k,label,info,Ntemp
    integer    :: ii,jj,kk
    integer    :: iatm
    integer    :: totGrid,radpot,sphpot
    integer    :: nr_quarter,cursphpot,row
    integer    :: acap
    integer,allocatable    :: grdRec(:,:)
    real(8),allocatable    :: potx(:),poty(:),potz(:),potw(:)
    real(8)    :: radx,radr,radw,d_min,ratio
    real(8)    :: parm
    real(8)    :: radr_pre,radw_pre
    real(8)    :: rmiu,tmps,sij
    integer    :: n_pairs,ip
    integer,allocatable :: pair_i(:),pair_j(:)
    integer,allocatable :: far_j_list(:), far_j_start(:), far_j_count(:)
    integer :: n_far, far_ii, far_jj
    real(8) :: PAIR_CUTOFF = 30.0d0
    character(len=32) :: pcenv
    real(8) :: pcval
    integer :: pcios
    real(8),allocatable    :: Rij(:,:),aij_mat(:,:)
    real(8),allocatable    :: rdist(:)
    real(8)    :: chi,uij,aij

    real(8),allocatable    :: Pvec(:)
    integer,allocatable    :: atom_offset(:)
    integer :: atom_npts

call get_environment_variable("ENGINE_XC_PAIR_CUTOFF", pcenv)
if (len_trim(pcenv) > 0) then
   read(pcenv, *, iostat=pcios) pcval
   if (pcios == 0 .and. pcval > 0.0d0) PAIR_CUTOFF = pcval
endif

   call gridgen_count_points(totGrid, grdRec)
   if (.not. engine_quiet) print '(A,I0)',"  VXC grid points: ",totGrid
   call flush(6)
   nGrids = totGrid
   if (allocated(Grids)) deallocate(Grids)
   allocate(Grids(totGrid))
   allocate(Rij(natoms,natoms))
   allocate(aij_mat(natoms,natoms))
   Rij = 0
   aij_mat = 0
   do ii = 1,natoms
      do jj = 1,natoms
         if (ii .eq. jj) cycle
         Rij(ii,jj) = dsqrt(sum((atoms(ii)%coor*ans2bohr - atoms(jj)%coor*ans2bohr)**2))
         if (xcgrid_use_treutler_aij) then
            chi = dsqrt(covrad(atoms(ii)%charge))/dsqrt(covrad(atoms(jj)%charge))
         else
            chi = covrad(atoms(ii)%charge)/covrad(atoms(jj)%charge)
         endif
         uij = (chi-1)/(chi+1)
         aij = uij/(uij**2 -1)
         if (aij >  0.5)  aij = 0.5
         if (aij < -0.5)  aij = -0.5
         aij_mat(ii,jj) = aij
      enddo
   enddo
   n_pairs = 0
   do ii = 1,natoms-1
      do jj = ii+1,natoms
         if (Rij(ii,jj) .lt. PAIR_CUTOFF) n_pairs = n_pairs + 1
      enddo
   enddo
   allocate(pair_i(n_pairs),pair_j(n_pairs))
   n_pairs = 0
   do ii = 1,natoms-1
      do jj = ii+1,natoms
         if (Rij(ii,jj) .lt. PAIR_CUTOFF) then
            n_pairs = n_pairs + 1
            pair_i(n_pairs) = ii
            pair_j(n_pairs) = jj
         endif
      enddo
   enddo
   allocate(far_j_start(natoms), far_j_count(natoms))
   n_far = 0
   do far_ii = 1,natoms
      far_j_start(far_ii) = n_far+1
      do far_jj = 1,natoms
         if (far_jj .eq. far_ii) cycle
         if (Rij(far_ii,far_jj) .ge. PAIR_CUTOFF) n_far = n_far+1
      enddo
      far_j_count(far_ii) = n_far - far_j_start(far_ii) + 1
   enddo
   allocate(far_j_list(n_far))
   n_far = 0
   do far_ii = 1,natoms
      do far_jj = 1,natoms
         if (far_jj .eq. far_ii) cycle
         if (Rij(far_ii,far_jj) .ge. PAIR_CUTOFF) then
            n_far = n_far+1
            far_j_list(n_far) = far_jj
         endif
      enddo
   enddo
   allocate(atom_offset(natoms))
   atom_npts = 0
   do iatm = 1,natoms
      atom_offset(iatm) = atom_npts
      select case (int(atoms(iatm)%charge))
      case (:2);   row = 1
      case (3:10); row = 2
      case default; row = 3
      end select
      acap = xcgrid_ang_cap(atoms(iatm)%npShell, atoms(iatm)%ndShell, &
                            atoms(iatm)%nfShell, atoms(iatm)%ngShell)
      if (force_dense) acap = 100000
      do i = 1,grdRec(iatm,1)
         call xcgrid_m3_radial(int(atoms(iatm)%charge), i, grdRec(iatm,1), radr_pre, radw_pre)
         ratio = radr_pre / (covrad(atoms(iatm)%charge)*ans2bohr)
         d_min = 0.0d0
         if (xcgrid_use_dmin_table .and. natoms .gt. 1) &
            call xcgrid_dmin_shell(radr_pre, Rij(iatm,:), iatm, natoms, d_min)
         atom_npts = atom_npts + xcgrid_atom_sphpot(grdRec(iatm,2), xcgrid_theory_anggrid(row), &
                        ratio, row, int(atoms(iatm)%charge), acap, &
                        xcgrid_use_dmin_table .and. natoms .gt. 1, d_min)
      enddo
   enddo

   !$omp parallel &
   !$omp&   private(iatm,radpot,sphpot,nr_quarter,cursphpot,row,ratio,d_min) &
   !$omp&   private(parm,i,j,radx,radr,radw,Ntemp,potx,poty,potz,potw) &
   !$omp&   private(label,ii,jj,ip,rdist,Pvec,rmiu,tmps,sij,far_jj)
   allocate(potx(1202),poty(1202),potz(1202),potw(1202))
   allocate(Pvec(natoms))
   allocate(rdist(natoms))
   !$omp do schedule(dynamic)
   do iatm = 1,natoms
      radpot = grdRec(iatm,1)
      sphpot = grdRec(iatm,2)
      label = atom_offset(iatm) + 1
      select case (int(atoms(iatm)%charge))
      case (:2);   row = 1
      case (3:10); row = 2
      case default; row = 3
      end select
      acap = xcgrid_ang_cap(atoms(iatm)%npShell, atoms(iatm)%ndShell, &
                            atoms(iatm)%nfShell, atoms(iatm)%ngShell)
      if (force_dense) acap = 100000
      do i = 1,radpot
         call xcgrid_m3_radial(int(atoms(iatm)%charge), i, radpot, radr, radw)
         ratio = radr / (covrad(atoms(iatm)%charge)*ans2bohr)
         d_min = 0.0d0
         if (xcgrid_use_dmin_table .and. natoms .gt. 1) &
            call xcgrid_dmin_shell(radr, Rij(iatm,:), iatm, natoms, d_min)
         cursphpot = xcgrid_atom_sphpot(sphpot, xcgrid_theory_anggrid(row), &
                        ratio, row, int(atoms(iatm)%charge), acap, &
                        xcgrid_use_dmin_table .and. natoms .gt. 1, d_min)
         if (xcgrid_dump_shells) then
            !$omp critical
            print '(A,I4,A,I4,A,I5,A,F12.6,A,ES12.4,A,I5)', &
               " [SHELL] atom=", iatm, " Z=", int(atoms(iatm)%charge), &
               " i=", i, " r=", radr, " w=", radw, " ang=", cursphpot
            !$omp end critical
         endif
         select case (cursphpot)
         case (974); call LD0974(potx,poty,potz,potw,Ntemp)
         case (770); call LD0770(potx,poty,potz,potw,Ntemp)
         case (590); call LD0590(potx,poty,potz,potw,Ntemp)
         case (434); call LD0434(potx,poty,potz,potw,Ntemp)
         case (350); call LD0350(potx,poty,potz,potw,Ntemp)
         case (302); call LD0302(potx,poty,potz,potw,Ntemp)
         case (266); call LD0266(potx,poty,potz,potw,Ntemp)
         case (230); call LD0230(potx,poty,potz,potw,Ntemp)
         case (194); call LD0194(potx,poty,potz,potw,Ntemp)
         case (170); call LD0170(potx,poty,potz,potw,Ntemp)
         case (146); call LD0146(potx,poty,potz,potw,Ntemp)
         case (110); call LD0110(potx,poty,potz,potw,Ntemp)
         case (86);  call LD0086(potx,poty,potz,potw,Ntemp)
         case (74);  call LD0074(potx,poty,potz,potw,Ntemp)
         case (50);  call LD0050(potx,poty,potz,potw,Ntemp)
         case (38);  call LD0038(potx,poty,potz,potw,Ntemp)
         case (26);  call LD0026(potx,poty,potz,potw,Ntemp)
         case (14);  call LD0014(potx,poty,potz,potw,Ntemp)
         case (6);   call LD0006(potx,poty,potz,potw,Ntemp)
         case default
            write(*,'(A,I0,A)') " grid_gen: no Lebedev grid for order ", &
                 cursphpot, " - aborting (pruning emitted an unsupported order)"
            stop 1
         end select
         do j = 1,Ntemp
            Grids(label)%coor(1)=radr*potx(j) + atoms(iatm)%coor(1)*ans2bohr
            Grids(label)%coor(2)=radr*poty(j) + atoms(iatm)%coor(2)*ans2bohr
            Grids(label)%coor(3)=radr*potz(j) + atoms(iatm)%coor(3)*ans2bohr
            do ii = 1,natoms
               rdist(ii) = dsqrt(sum((Grids(label)%coor - atoms(ii)%coor*ans2bohr)**2))
            enddo
            Pvec = 1.0D0
            do ip = 1,n_pairs
               ii = pair_i(ip)
               jj = pair_j(ip)
               rmiu = (rdist(ii)-rdist(jj))/Rij(ii,jj)
               rmiu = rmiu + aij_mat(ii,jj)*(1-rmiu**2)

               tmps = 1.5D0 *(rmiu) - 0.5D0*(rmiu)**3
               tmps = 1.5D0 *(tmps) - 0.5D0*(tmps)**3
               tmps = 1.5D0 *(tmps) - 0.5D0*(tmps)**3
               sij = 0.5D0 * (1-tmps)
               Pvec(ii) = Pvec(ii) * sij
               Pvec(jj) = Pvec(jj) * (1.0D0 - sij)
            enddo
            do ip = far_j_start(iatm),far_j_start(iatm)+far_j_count(iatm)-1
               far_jj = far_j_list(ip)
               rmiu = (rdist(iatm)-rdist(far_jj))/Rij(iatm,far_jj)
               rmiu = rmiu + aij_mat(iatm,far_jj)*(1-rmiu**2)

               tmps = 1.5D0 *(rmiu) - 0.5D0*(rmiu)**3
               tmps = 1.5D0 *(tmps) - 0.5D0*(tmps)**3
               tmps = 1.5D0 *(tmps) - 0.5D0*(tmps)**3
               sij = 0.5D0 * (1-tmps)
               Pvec(iatm) = Pvec(iatm) * sij
               Pvec(far_jj) = Pvec(far_jj) * (1.0D0 - sij)
            enddo
            Grids(label)%weight =radw*potw(j)*Pvec(iatm)/sum(Pvec)
            label = 1 + label
         enddo
      enddo
   enddo
   !$omp end do
   deallocate(potx,poty,potz,potw)
   deallocate(Pvec)
   deallocate(rdist)
   !$omp end parallel
   deallocate(grdRec)
   deallocate(atom_offset)
   deallocate(Rij)
   deallocate(aij_mat)
   deallocate(pair_i,pair_j)
   deallocate(far_j_list,far_j_start,far_j_count)

   call spatial_sort_grid(totGrid)
end subroutine

subroutine gridgen_nlc(nrad, nsph, npts, coor_out, weight_out, coor_in, per_atom_period_scale, intacc_eps)
use MOL_info
use GRID_info, only: xcgrid_prune_sphpot, gridgen_nlc_atom_of, xcgrid_m3_radial
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: nrad, nsph
integer,intent(out) :: npts
real(8),allocatable,intent(out) :: coor_out(:,:), weight_out(:)
real(8),optional,intent(in) :: coor_in(3,natoms)
logical,optional,intent(in) :: per_atom_period_scale
real(8),optional,intent(in) :: intacc_eps
logical :: use_per_atom, use_intacc
integer :: row_i
integer :: nrad_atom(natoms), nsph_atom(natoms), nr_quarter_atom(natoms)
integer :: atom_pts_i(natoms), atom_off(natoms), period_i
real(8) :: Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8) :: acoor(3,natoms)
real(8) :: parm, radx, radr, radw, chi, uij, aij
real(8) :: potx(1202),poty(1202),potz(1202),potw(1202)
real(8) :: rdist(natoms), Pvec(natoms), rmiu, tmps, sij
integer :: iatm, i, j, Ntemp, ii, jj, ip, label
integer :: nr_quarter, cursphpot, atom_pts

if (present(coor_in)) then
   if (engine_verbose .ge. 2) then
      print *, "NLC gen present coor_in, natoms=",natoms
      call flush(6)
   endif
   acoor = coor_in
else
   do iatm = 1,natoms
      acoor(:,iatm) = atoms(iatm)%coor*ans2bohr
   enddo
endif

use_per_atom = present(per_atom_period_scale)
if (use_per_atom) use_per_atom = per_atom_period_scale
use_intacc = present(intacc_eps)
do iatm = 1,natoms
   nrad_atom(iatm) = nrad
   nsph_atom(iatm) = nsph
   if (use_intacc) then
      if (atoms(iatm)%charge .le. 2) then
         row_i = 1
      else if (atoms(iatm)%charge .le. 10) then
         row_i = 2
      else if (atoms(iatm)%charge .le. 18) then
         row_i = 3
      else if (atoms(iatm)%charge .le. 36) then
         row_i = 4
      else if (atoms(iatm)%charge .le. 54) then
         row_i = 5
      else if (atoms(iatm)%charge .le. 86) then
         row_i = 6
      else
         row_i = 7
      endif
      nrad_atom(iatm) = max(1, int(intacc_eps*15.0d0 - 40.0d0 + 5.0d0*row_i))
   endif
   if (use_per_atom) then
      if (atoms(iatm)%charge .le. 10) then
         period_i = 2
      else if (atoms(iatm)%charge .le. 18) then
         period_i = 3
      else if (atoms(iatm)%charge .le. 36) then
         period_i = 4
      else if (atoms(iatm)%charge .le. 54) then
         period_i = 5
      else
         period_i = 6
      endif
      if (period_i .ge. 3) nrad_atom(iatm) = nrad_atom(iatm)*3
      if (period_i .ge. 4) nsph_atom(iatm) = 974
   endif
   nr_quarter_atom(iatm) = max(1,nrad_atom(iatm)/4)
enddo
npts = 0
do iatm = 1,natoms
   atom_pts_i(iatm) = 0
   do i = 1,nrad_atom(iatm)
      atom_pts_i(iatm) = atom_pts_i(iatm) + xcgrid_prune_sphpot(nsph_atom(iatm), i, nr_quarter_atom(iatm))
   enddo
   atom_off(iatm) = npts
   npts = npts + atom_pts_i(iatm)
enddo
if (engine_verbose .ge. 2) then
   print *, "NLC gen nrad=",nrad," nsph=",nsph," npts=",npts," (pruned from ",natoms*nrad*nsph,")"
   if (use_per_atom) print *, "NLC gen per-atom period scaling active"
   call flush(6)
endif
allocate(coor_out(3,npts))
allocate(weight_out(npts))
if (allocated(gridgen_nlc_atom_of)) deallocate(gridgen_nlc_atom_of)
allocate(gridgen_nlc_atom_of(npts))

Rij = 0.0d0
aij_mat = 0.0d0
do ii = 1,natoms
   do jj = 1,natoms
      if (ii .eq. jj) cycle
      Rij(ii,jj) = dsqrt(sum((acoor(:,ii) - acoor(:,jj))**2))
      chi = covrad(atoms(ii)%charge)/covrad(atoms(jj)%charge)
      uij = (chi-1)/(chi+1)
      aij = uij/(uij**2 -1)
      if (aij >  0.5)  aij = 0.5
      if (aij < -0.5)  aij = -0.5
      aij_mat(ii,jj) = aij
   enddo
enddo

!$omp parallel do default(shared) private(iatm,parm,i,cursphpot,potx,poty,potz,potw,Ntemp) &
!$omp&   private(radx,radr,radw,j,label,ii,jj,rdist,Pvec,rmiu,tmps,sij) schedule(dynamic)
do iatm = 1,natoms
   label = atom_off(iatm)
   do i = 1,nrad_atom(iatm)
      cursphpot = xcgrid_prune_sphpot(nsph_atom(iatm), i, nr_quarter_atom(iatm))
      call ld_dispatch(cursphpot, potx,poty,potz,potw, Ntemp)
      call xcgrid_m3_radial(int(atoms(iatm)%charge), i, nrad_atom(iatm), radr, radw)
      do j = 1,Ntemp
         label = label + 1
         coor_out(1,label) = radr*potx(j) + acoor(1,iatm)
         coor_out(2,label) = radr*poty(j) + acoor(2,iatm)
         coor_out(3,label) = radr*potz(j) + acoor(3,iatm)
         gridgen_nlc_atom_of(label) = iatm
         do ii = 1,natoms
            rdist(ii) = dsqrt(sum((coor_out(:,label) - acoor(:,ii))**2))
         enddo
         Pvec = 1.0d0
         do ii = 1,natoms-1
            do jj = ii+1,natoms
               rmiu = (rdist(ii)-rdist(jj))/Rij(ii,jj)
               rmiu = rmiu + aij_mat(ii,jj)*(1-rmiu**2)
               tmps = 1.5d0*rmiu - 0.5d0*rmiu**3
               tmps = 1.5d0*tmps  - 0.5d0*tmps**3
               tmps = 1.5d0*tmps  - 0.5d0*tmps**3
               sij = 0.5d0*(1-tmps)
               Pvec(ii) = Pvec(ii)*sij
               Pvec(jj) = Pvec(jj)*(1.0d0-sij)
            enddo
         enddo
         weight_out(label) = radw*potw(j)*Pvec(iatm)/sum(Pvec)
      enddo
   enddo
enddo
!$omp end parallel do

contains

subroutine ld_dispatch(nsph_in, px,py,pz,pw, nt)
   implicit none
   integer,intent(in) :: nsph_in
   real(8),intent(out) :: px(974),py(974),pz(974),pw(974)
   integer,intent(out) :: nt
   select case (nsph_in)
   case (6);   call LD0006(px,py,pz,pw,nt)
   case (14);  call LD0014(px,py,pz,pw,nt)
   case (26);  call LD0026(px,py,pz,pw,nt)
   case (38);  call LD0038(px,py,pz,pw,nt)
   case (50);  call LD0050(px,py,pz,pw,nt)
   case (74);  call LD0074(px,py,pz,pw,nt)
   case (86);  call LD0086(px,py,pz,pw,nt)
   case (110); call LD0110(px,py,pz,pw,nt)
   case (146); call LD0146(px,py,pz,pw,nt)
   case (170); call LD0170(px,py,pz,pw,nt)
   case (194); call LD0194(px,py,pz,pw,nt)
   case (230); call LD0230(px,py,pz,pw,nt)
   case (266); call LD0266(px,py,pz,pw,nt)
   case (302); call LD0302(px,py,pz,pw,nt)
   case (350); call LD0350(px,py,pz,pw,nt)
   case (434); call LD0434(px,py,pz,pw,nt)
   case (590); call LD0590(px,py,pz,pw,nt)
   case (770); call LD0770(px,py,pz,pw,nt)
   case (974); call LD0974(px,py,pz,pw,nt)
   case default
      print *, "gridgen_nlc: unsupported spherical order", nsph_in, "- using 110"
      call LD0110(px,py,pz,pw,nt)
   end select
end subroutine ld_dispatch

end subroutine gridgen_nlc

subroutine gridgen_nlc_subset(nrad, nsph, atom_list, n_list, npts, coor_out, weight_out, atom_of_out, coor_in)
use MOL_info
use GRID_info, only: xcgrid_prune_sphpot, gridgen_nlc_atom_of, xcgrid_m3_radial
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: nrad, nsph, n_list
integer,intent(in) :: atom_list(n_list)
integer,intent(out) :: npts
real(8),allocatable,intent(out) :: coor_out(:,:), weight_out(:)
integer,allocatable,intent(out) :: atom_of_out(:)
real(8),optional,intent(in) :: coor_in(3,natoms)
real(8) :: Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8) :: acoor(3,natoms)
real(8) :: parm, radx, radr, radw, chi, uij, aij
real(8) :: potx(1202),poty(1202),potz(1202),potw(1202)
real(8) :: rdist(natoms), Pvec(natoms), rmiu, tmps, sij
integer :: iatm, i, j, Ntemp, ii, jj, ip, label, ial
integer :: nr_quarter, cursphpot, atom_pts

if (present(coor_in)) then
   acoor = coor_in
else
   do iatm = 1,natoms
      acoor(:,iatm) = atoms(iatm)%coor*ans2bohr
   enddo
endif

nr_quarter = max(1,nrad/4)
atom_pts = 0
do i = 1,nrad
   atom_pts = atom_pts + xcgrid_prune_sphpot(nsph, i, nr_quarter)
enddo
npts = n_list*atom_pts
allocate(coor_out(3,npts))
allocate(weight_out(npts))
allocate(atom_of_out(npts))

Rij = 0.0d0
aij_mat = 0.0d0
do ii = 1,natoms
   do jj = 1,natoms
      if (ii .eq. jj) cycle
      Rij(ii,jj) = dsqrt(sum((acoor(:,ii) - acoor(:,jj))**2))
      chi = covrad(atoms(ii)%charge)/covrad(atoms(jj)%charge)
      uij = (chi-1)/(chi+1)
      aij = uij/(uij**2 -1)
      if (aij >  0.5)  aij = 0.5
      if (aij < -0.5)  aij = -0.5
      aij_mat(ii,jj) = aij
   enddo
enddo

label = 0
do ial = 1,n_list
   iatm = atom_list(ial)
   do i = 1,nrad
      cursphpot = xcgrid_prune_sphpot(nsph, i, nr_quarter)
      call ld_dispatch_sub(cursphpot, potx,poty,potz,potw, Ntemp)
      call xcgrid_m3_radial(int(atoms(iatm)%charge), i, nrad, radr, radw)
      do j = 1,Ntemp
         label = label + 1
         coor_out(1,label) = radr*potx(j) + acoor(1,iatm)
         coor_out(2,label) = radr*poty(j) + acoor(2,iatm)
         coor_out(3,label) = radr*potz(j) + acoor(3,iatm)
         atom_of_out(label) = iatm
         do ii = 1,natoms
            rdist(ii) = dsqrt(sum((coor_out(:,label) - acoor(:,ii))**2))
         enddo
         Pvec = 1.0d0
         do ii = 1,natoms-1
            do jj = ii+1,natoms
               rmiu = (rdist(ii)-rdist(jj))/Rij(ii,jj)
               rmiu = rmiu + aij_mat(ii,jj)*(1-rmiu**2)
               tmps = 1.5d0*rmiu - 0.5d0*rmiu**3
               tmps = 1.5d0*tmps  - 0.5d0*tmps**3
               tmps = 1.5d0*tmps  - 0.5d0*tmps**3
               sij = 0.5d0*(1-tmps)
               Pvec(ii) = Pvec(ii)*sij
               Pvec(jj) = Pvec(jj)*(1.0d0-sij)
            enddo
         enddo
         weight_out(label) = radw*potw(j)*Pvec(iatm)/sum(Pvec)
      enddo
   enddo
enddo

contains

subroutine ld_dispatch_sub(nsph_in, px,py,pz,pw, nt)
   implicit none
   integer,intent(in) :: nsph_in
   real(8),intent(out) :: px(974),py(974),pz(974),pw(974)
   integer,intent(out) :: nt
   select case (nsph_in)
   case (6);   call LD0006(px,py,pz,pw,nt)
   case (14);  call LD0014(px,py,pz,pw,nt)
   case (26);  call LD0026(px,py,pz,pw,nt)
   case (38);  call LD0038(px,py,pz,pw,nt)
   case (50);  call LD0050(px,py,pz,pw,nt)
   case (74);  call LD0074(px,py,pz,pw,nt)
   case (86);  call LD0086(px,py,pz,pw,nt)
   case (110); call LD0110(px,py,pz,pw,nt)
   case (146); call LD0146(px,py,pz,pw,nt)
   case (170); call LD0170(px,py,pz,pw,nt)
   case (194); call LD0194(px,py,pz,pw,nt)
   case (230); call LD0230(px,py,pz,pw,nt)
   case (266); call LD0266(px,py,pz,pw,nt)
   case (302); call LD0302(px,py,pz,pw,nt)
   case (350); call LD0350(px,py,pz,pw,nt)
   case (434); call LD0434(px,py,pz,pw,nt)
   case (590); call LD0590(px,py,pz,pw,nt)
   case (770); call LD0770(px,py,pz,pw,nt)
   case (974); call LD0974(px,py,pz,pw,nt)
   case default
      print *, "gridgen_nlc_subset: unsupported spherical order", nsph_in, "- using 110"
      call LD0110(px,py,pz,pw,nt)
   end select
end subroutine ld_dispatch_sub

end subroutine gridgen_nlc_subset

subroutine vv10_nlc_weight_gradient(np, coor, atom_of, nper, rho, eps, wquad, beta_v, gv_wt, dwdr_out)
use MOL_info
use GRID_info, only: cosx_becke_treutler
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: np, nper
real(8),intent(in) :: coor(3,np)
integer,intent(in) :: atom_of(np)
real(8),intent(in) :: rho(np), eps(np), wquad(np), beta_v
real(8),intent(out) :: gv_wt(natoms,3)
real(8),optional,intent(out) :: dwdr_out(3,natoms,np)

real(8) :: Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8) :: acoor(3,natoms)
real(8) :: chi, uij_v, aij
integer :: i2, jj2
real(8),parameter :: SAFE_S_FLOOR = 1.0d-280
real(8),allocatable :: invR(:,:), ge(:,:,:)

do i2 = 1,natoms
   acoor(:,i2) = atoms(i2)%coor*ans2bohr
enddo
Rij = 0.0d0
aij_mat = 0.0d0
do i2 = 1,natoms
   do jj2 = 1,natoms
      if (i2 .eq. jj2) cycle
      Rij(i2,jj2) = dsqrt(sum((acoor(:,i2)-acoor(:,jj2))**2))
      chi = covrad(atoms(i2)%charge)/covrad(atoms(jj2)%charge)
      if (cosx_becke_treutler) chi = sqrt(chi)
      uij_v = (chi-1.0d0)/(chi+1.0d0)
      aij = uij_v/(uij_v**2-1.0d0)
      if (aij >  0.5d0) aij =  0.5d0
      if (aij < -0.5d0) aij = -0.5d0
      aij_mat(i2,jj2) = aij
   enddo
enddo

gv_wt = 0.0d0
if (present(dwdr_out)) dwdr_out = 0.0d0

allocate(invR(natoms,natoms), ge(natoms,natoms,3))
invR = 0.0d0
ge = 0.0d0
do jj2 = 1,natoms
   do i2 = 1,natoms
      if (i2 .eq. jj2) cycle
      invR(i2,jj2) = 1.0d0/Rij(i2,jj2)
      ge(i2,jj2,:) = (acoor(:,i2)-acoor(:,jj2))*invR(i2,jj2)**2
   enddo
enddo

!$omp parallel do schedule(dynamic) reduction(+:gv_wt)
do i2 = 1,np
   block
      integer :: O, A, B, C, dd
      real(8) :: rdist(natoms), uvec(3,natoms)
      real(8) :: mu(natoms,natoms), sAB(natoms,natoms), gfac(natoms,natoms)
      real(8) :: x1, x2, x3
      real(8) :: Pvec(natoms), sumP
      real(8) :: factor, t, tm, a1, a2, bs(3), cc(3), S1, S2, dsumP, dlogPO, val, acc(3)

      if (dabs(wquad(i2)) .lt. SAFE_S_FLOOR) cycle
      O = atom_of(i2)
      do A = 1,natoms
         rdist(A) = dsqrt(sum((coor(:,i2)-acoor(:,A))**2))
         uvec(:,A) = (coor(:,i2)-acoor(:,A))/rdist(A)
      enddo
      do B = 1,natoms
         do A = 1,natoms
            if (A .eq. B) then
               mu(A,B) = 0.0d0 ; sAB(A,B) = 1.0d0 ; gfac(A,B) = 0.0d0
               cycle
            endif
            mu(A,B) = (rdist(A)-rdist(B))/Rij(A,B)
            x1 = mu(A,B) + aij_mat(A,B)*(1.0d0-mu(A,B)**2)
            x2 = 1.5d0*x1 - 0.5d0*x1**3
            x3 = 1.5d0*x2 - 0.5d0*x2**3
            sAB(A,B) = 0.5d0*(1.0d0 - (1.5d0*x3-0.5d0*x3**3))
            if (sAB(A,B) .lt. SAFE_S_FLOOR) then
               gfac(A,B) = 0.0d0
            else
               gfac(A,B) = -1.6875d0*(1.0d0-x1**2)*(1.0d0-x2**2)*(1.0d0-x3**2) &
                           *(1.0d0 - 2.0d0*aij_mat(A,B)*mu(A,B))/sAB(A,B)
            endif
         enddo
      enddo
      Pvec = 1.0d0
      do B = 1,natoms
         do A = 1,natoms
            Pvec(A) = Pvec(A)*sAB(A,B)
         enddo
      enddo
      sumP = sum(Pvec)
      factor = rho(i2)*(2.0d0*eps(i2)-beta_v)

      acc = 0.0d0
      do C = 1,natoms
         if (C .eq. O) cycle
         a1 = 0.0d0 ; bs = 0.0d0
         do A = 1,natoms
            t = Pvec(A)*gfac(A,C)
            a1 = a1 + t*invR(A,C)
            tm = t*mu(A,C)
            bs(1) = bs(1) + tm*ge(A,C,1)
            bs(2) = bs(2) + tm*ge(A,C,2)
            bs(3) = bs(3) + tm*ge(A,C,3)
         enddo
         a2 = 0.0d0 ; cc = 0.0d0
         do B = 1,natoms
            t = gfac(C,B)
            a2 = a2 + t*invR(C,B)
            tm = t*mu(C,B)
            cc(1) = cc(1) + tm*ge(C,B,1)
            cc(2) = cc(2) + tm*ge(C,B,2)
            cc(3) = cc(3) + tm*ge(C,B,3)
         enddo
         do dd = 1,3
            S1 = uvec(dd,C)*a1 + bs(dd)
            S2 = -uvec(dd,C)*a2 - cc(dd)
            dsumP = S1 + Pvec(C)*S2
            dlogPO = gfac(O,C)*(uvec(dd,C)*invR(O,C) + mu(O,C)*ge(O,C,dd))
            if (sumP .gt. 0.0d0) then
               val = wquad(i2)*(dlogPO - dsumP/sumP)
            else
               val = wquad(i2)*dlogPO
            endif
            gv_wt(C,dd) = gv_wt(C,dd) + factor*val
            acc(dd) = acc(dd) + val
            if (present(dwdr_out)) dwdr_out(dd,C,i2) = val
         enddo
      enddo
      do dd = 1,3
         gv_wt(O,dd) = gv_wt(O,dd) - factor*acc(dd)
         if (present(dwdr_out)) dwdr_out(dd,O,i2) = -acc(dd)
      enddo
   end block
enddo
!$omp end parallel do
deallocate(invR, ge)

end subroutine vv10_nlc_weight_gradient

subroutine atom_pair_geom(acoor, Rij, aij_mat)
use MOL_info
use GRID_info, only: cosx_becke_treutler
implicit none
INCLUDE 'parameter.h'
real(8),intent(in) :: acoor(3,natoms)
real(8),intent(out) :: Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8) :: chi, uij_v, aij
integer :: i, p
Rij = 0.0d0 ; aij_mat = 0.0d0
do i = 1,natoms
   do p = 1,natoms
      if (i .eq. p) cycle
      Rij(i,p) = dsqrt(sum((acoor(:,i)-acoor(:,p))**2))
      chi = covrad(atoms(i)%charge)/covrad(atoms(p)%charge)
      if (cosx_becke_treutler) chi = sqrt(chi)
      uij_v = (chi-1.0d0)/(chi+1.0d0)
      aij = uij_v/(uij_v**2-1.0d0)
      if (aij >  0.5d0) aij =  0.5d0
      if (aij < -0.5d0) aij = -0.5d0
      aij_mat(i,p) = aij
   enddo
enddo
end subroutine atom_pair_geom

subroutine becke_fractions(np, coor, atom_of, acoor, Rij, aij_mat, frac)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: np
real(8),intent(in) :: coor(3,*)
integer,intent(in) :: atom_of(*)
real(8),intent(in) :: acoor(3,natoms), Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8),intent(out) :: frac(*)
integer :: p, O, A, B
real(8) :: rdist(natoms), Pvec(natoms), mu, x1, x2, x3, sv
do p = 1,np
   O = atom_of(p)
   do A = 1,natoms
      rdist(A) = dsqrt(sum((coor(:,p)-acoor(:,A))**2))
   enddo
   Pvec = 1.0d0
   do A = 1,natoms
      do B = 1,natoms
         if (A .eq. B) cycle
         mu = (rdist(A)-rdist(B))/Rij(A,B)
         x1 = mu + aij_mat(A,B)*(1.0d0-mu**2)
         x2 = 1.5d0*x1 - 0.5d0*x1**3
         x3 = 1.5d0*x2 - 0.5d0*x2**3
         sv = 0.5d0*(1.0d0-(1.5d0*x3-0.5d0*x3**3))
         Pvec(A) = Pvec(A)*sv
      enddo
   enddo
   frac(p) = Pvec(O)/sum(Pvec)
enddo
end subroutine becke_fractions

subroutine vv10_nlc_weight_hessian(np, coor, atom_of, wquad, fct, hess_out)
use MOL_info
implicit none
integer,intent(in) :: np
real(8),intent(in) :: coor(3,np)
integer,intent(in) :: atom_of(np)
real(8),intent(in) :: wquad(np), fct(np)
real(8),intent(out) :: hess_out(natoms,3,natoms,3)
call vv10_nlc_weight_hessian_multi(np, coor, atom_of, wquad, 1, fct, hess_out)
end subroutine vv10_nlc_weight_hessian

subroutine vv10_nlc_weight_hessian_multi(np, coor, atom_of, wquad, nf, fct, hess_out)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: np, nf
real(8),intent(in) :: coor(3,np)
integer,intent(in) :: atom_of(np)
real(8),intent(in) :: wquad(np), fct(np,nf)
real(8),intent(out) :: hess_out(natoms,3,natoms,3,nf)
real(8) :: acoor(3,natoms), Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8),allocatable :: hloc(:,:,:,:,:)
integer :: i2, ii, tid
integer :: omp_get_thread_num, omp_get_num_threads
external :: omp_get_thread_num, omp_get_num_threads
hess_out = 0.0d0
do ii = 1,natoms ; acoor(:,ii) = atoms(ii)%coor*ans2bohr ; enddo
call atom_pair_geom(acoor, Rij, aij_mat)
!$omp parallel default(shared) private(i2, hloc, tid)
allocate(hloc(natoms,3,natoms,3,nf))
hloc = 0.0d0
!$omp do schedule(static)
do i2 = 1,np
   if (all(fct(i2,:) .eq. 0.0d0)) cycle
   call weight_hess_point(i2, np, nf, coor, atom_of, wquad, fct, acoor, Rij, aij_mat, hloc)
enddo
!$omp end do
do tid = 0, omp_get_num_threads()-1
   !$omp barrier
   if (omp_get_thread_num() .eq. tid) hess_out = hess_out + hloc
enddo
deallocate(hloc)
!$omp end parallel
end subroutine vv10_nlc_weight_hessian_multi

subroutine weight_hess_point(i2, np, nf, coor, atom_of, wquad, fct, acoor, Rij, aij_mat, hess_out)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: i2, np, nf
real(8),intent(in) :: coor(3,*)
integer,intent(in) :: atom_of(*)
real(8),intent(in) :: wquad(*), fct(np,nf)
real(8),intent(in) :: acoor(3,natoms), Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8),intent(inout) :: hess_out(natoms,3,natoms,3,nf)
real(8),parameter :: SAFE_S_FLOOR = 1.0d-280
real(8),allocatable,save :: rdist(:), uvec(:,:), wAB(:,:,:), mu(:,:)
real(8),allocatable,save :: gAB(:,:), gpAB(:,:), sAB(:,:)
real(8),allocatable,save :: Pvec(:), DlogP(:,:,:), DsumP(:,:)
real(8),allocatable,save :: D2logPO(:,:,:,:), D2sumP(:,:,:,:), PD(:,:,:), DLc(:,:,:), gW(:,:)
!$omp threadprivate(rdist, uvec, wAB, mu, gAB, gpAB, sAB, Pvec, DlogP, DsumP, D2logPO, D2sumP, PD, DLc, gW)
real(8) :: ds(3,3), dm(3,3), d2s(3,3,3,3), prA(3,3), prB(3,3), prR(3,3), vDel(3,3), vR(3,3), rinv
real(8) :: caA(3), caB(3), cABv(3)
integer :: idx(3)
real(8) :: rA, rB, muAB, x1, x2, x3, f3, ss, q, p1, p2, p3, x1p, x2p, x3p, dp
real(8) :: gv, gpv, d2rA, d2rB, d2Del, d2R, dDel, dDelp, dR, dRp, d2mu, dde
real(8) :: sumP, fac, d2logW, dlogW_C, dlogW_Cp, term
integer :: Oi, A, B, C, Cp, dd, ee, u, v, m
logical :: pact(natoms)
if (allocated(rdist)) then
   if (size(rdist) .ne. natoms) deallocate(rdist, uvec, wAB, mu, gAB, gpAB, sAB, Pvec, DlogP, DsumP, &
                                           D2logPO, D2sumP, PD, DLc, gW)
endif
if (.not. allocated(rdist)) then
   allocate(rdist(natoms), uvec(3,natoms), wAB(3,natoms,natoms), mu(natoms,natoms))
   allocate(gAB(natoms,natoms), gpAB(natoms,natoms), sAB(natoms,natoms))
   allocate(Pvec(natoms), DlogP(natoms,natoms,3), DsumP(natoms,3))
   allocate(D2logPO(natoms,3,natoms,3), D2sumP(natoms,3,natoms,3), PD(natoms,natoms,3), DLc(natoms,natoms,3), gW(natoms,3))
endif
Oi = atom_of(i2)
do A = 1,natoms
   rdist(A) = dsqrt(sum((coor(:,i2)-acoor(:,A))**2))
   uvec(:,A) = (coor(:,i2)-acoor(:,A))/rdist(A)
enddo
do A = 1,natoms
   do B = 1,natoms
      if (A .eq. B) cycle
      wAB(:,A,B) = (acoor(:,A)-acoor(:,B))/Rij(A,B)
      muAB = (rdist(A)-rdist(B))/Rij(A,B)
      mu(A,B) = muAB
      x1 = muAB + aij_mat(A,B)*(1.0d0-muAB**2)
      x2 = 1.5d0*x1 - 0.5d0*x1**3
      x3 = 1.5d0*x2 - 0.5d0*x2**3
      f3 = 1.5d0*x3 - 0.5d0*x3**3
      sAB(A,B) = 0.5d0*(1.0d0-f3)
      q = 1.0d0 - 2.0d0*aij_mat(A,B)*muAB
      p1 = 1.0d0-x1*x1 ; p2 = 1.0d0-x2*x2 ; p3 = 1.0d0-x3*x3
      gAB(A,B) = -1.6875d0*p1*p2*p3*q
      x1p = q
      x2p = 1.5d0*p1*q
      x3p = 1.5d0*p2*x2p
      dp = -2.0d0*( x1*x1p*p2*p3 + p1*x2*x2p*p3 + p1*p2*x3*x3p )
      gpAB(A,B) = -1.6875d0*(dp*q - 2.0d0*aij_mat(A,B)*p1*p2*p3)
   enddo
enddo
Pvec = 1.0d0
do A = 1,natoms
   do B = 1,natoms
      if (A .eq. B) cycle
      Pvec(A) = Pvec(A)*sAB(A,B)
   enddo
enddo
sumP = sum(Pvec)
if (dabs(wquad(i2)) .lt. SAFE_S_FLOOR) then
   return
endif
DlogP = 0.0d0 ; D2logPO = 0.0d0 ; D2sumP = 0.0d0
do A = 1,natoms
   pact(A) = (A .eq. Oi) .or. (Pvec(A) .ge. 1.0d-14*sumP)
enddo
do A = 1,natoms
   do B = A+1,natoms
      if (.not. (pact(A) .or. pact(B))) cycle
      idx(1) = A ; idx(2) = B ; m = 2
      if (Oi .ne. A .and. Oi .ne. B) then
         m = 3 ; idx(3) = Oi
      endif
      do u = 1,m
         C = idx(u)
         caA(u) = merge(1.0d0,0.0d0,C .eq. Oi) - merge(1.0d0,0.0d0,C .eq. A)
         caB(u) = merge(1.0d0,0.0d0,C .eq. Oi) - merge(1.0d0,0.0d0,C .eq. B)
         cABv(u)= merge(1.0d0,0.0d0,C .eq. A) - merge(1.0d0,0.0d0,C .eq. B)
      enddo
      rA = rdist(A) ; rB = rdist(B) ; muAB = mu(A,B)
      gv = gAB(A,B) ; gpv = gpAB(A,B)
      rinv = 1.0d0/Rij(A,B)
      do dd = 1,3
         do ee = 1,3
            dde = merge(1.0d0,0.0d0,dd .eq. ee)
            prA(dd,ee) = (dde - uvec(dd,A)*uvec(ee,A))/rA
            prB(dd,ee) = (dde - uvec(dd,B)*uvec(ee,B))/rB
            prR(dd,ee) = (dde - wAB(dd,A,B)*wAB(ee,A,B))*rinv
         enddo
      enddo
      do u = 1,m
         do dd = 1,3
            vDel(u,dd) = caA(u)*uvec(dd,A) - caB(u)*uvec(dd,B)
            vR(u,dd) = cABv(u)*wAB(dd,A,B)
            dm(u,dd) = (vDel(u,dd) - muAB*vR(u,dd))*rinv
            ds(u,dd) = gv*dm(u,dd)
         enddo
      enddo
      d2s = 0.0d0
      do v = 1,m ; do u = 1,m
         do ee = 1,3 ; do dd = 1,3
            d2mu = (caA(u)*caA(v)*prA(dd,ee) - caB(u)*caB(v)*prB(dd,ee))*rinv &
                   - (vDel(v,ee)*vR(u,dd) + vDel(u,dd)*vR(v,ee))*rinv*rinv &
                   - muAB*cABv(u)*cABv(v)*prR(dd,ee)*rinv &
                   + 2.0d0*muAB*vR(u,dd)*vR(v,ee)*rinv*rinv
            d2s(u,dd,v,ee) = gpv*dm(u,dd)*dm(v,ee) + gv*d2mu
         enddo ; enddo
      enddo ; enddo
      ss = sAB(A,B)
      if (ss .ge. SAFE_S_FLOOR) then
         do u = 1,m
            C = idx(u)
            do dd = 1,3
               DlogP(A,C,dd) = DlogP(A,C,dd) + ds(u,dd)/ss
            enddo
         enddo
         do u = 1,m ; do v = 1,m ; do dd = 1,3 ; do ee = 1,3
            term = d2s(u,dd,v,ee)/ss - ds(u,dd)*ds(v,ee)/ss**2
            D2sumP(idx(u),dd,idx(v),ee) = D2sumP(idx(u),dd,idx(v),ee) + Pvec(A)*term
            if (A .eq. Oi) D2logPO(idx(u),dd,idx(v),ee) = D2logPO(idx(u),dd,idx(v),ee) + term
         enddo ; enddo ; enddo ; enddo
      endif
      ss = 1.0d0 - sAB(A,B)
      if (ss .ge. SAFE_S_FLOOR) then
         do u = 1,m
            C = idx(u)
            do dd = 1,3
               DlogP(B,C,dd) = DlogP(B,C,dd) - ds(u,dd)/ss
            enddo
         enddo
         do u = 1,m ; do v = 1,m ; do dd = 1,3 ; do ee = 1,3
            term = -d2s(u,dd,v,ee)/ss - ds(u,dd)*ds(v,ee)/ss**2
            D2sumP(idx(u),dd,idx(v),ee) = D2sumP(idx(u),dd,idx(v),ee) + Pvec(B)*term
            if (B .eq. Oi) D2logPO(idx(u),dd,idx(v),ee) = D2logPO(idx(u),dd,idx(v),ee) + term
         enddo ; enddo ; enddo ; enddo
      endif
   enddo
enddo
DsumP = 0.0d0
do A = 1,natoms
   do C = 1,natoms
      do dd = 1,3
         DsumP(C,dd) = DsumP(C,dd) + Pvec(A)*DlogP(A,C,dd)
      enddo
   enddo
enddo
block
   integer :: na
   na = 0
   do A = 1,natoms
      if (.not. pact(A)) cycle
      na = na + 1
      PD(na,:,:) = Pvec(A)*DlogP(A,:,:)
      DLc(na,:,:) = DlogP(A,:,:)
   enddo
   call dgemm('T','N',3*natoms,3*natoms,na,1.0d0,PD,natoms,DLc,natoms,1.0d0,D2sumP,3*natoms)
end block
fac = 1.0d0/sumP
do ee = 1,3 ; do C = 1,natoms
   gW(C,ee) = DlogP(Oi,C,ee) - DsumP(C,ee)*fac
enddo ; enddo
do ee = 1,3 ; do Cp = 1,natoms ; do dd = 1,3 ; do C = 1,natoms
   D2sumP(C,dd,Cp,ee) = D2logPO(C,dd,Cp,ee) &
          - (D2sumP(C,dd,Cp,ee)*fac - DsumP(C,dd)*DsumP(Cp,ee)*fac*fac) + gW(C,dd)*gW(Cp,ee)
enddo ; enddo ; enddo ; enddo
do m = 1,nf
   dlogW_C = fct(i2,m)*wquad(i2)
   hess_out(:,:,:,:,m) = hess_out(:,:,:,:,m) + dlogW_C*D2sumP
enddo
end subroutine weight_hess_point

subroutine weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, Fsum)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: np
real(8),intent(in) :: coor0(3,np), baseat(3,natoms), shift(3,natoms), fct(np)
integer,intent(in) :: atom_of(np)
real(8),intent(inout) :: acoor(3,natoms), Rij(natoms,natoms), aij(natoms,natoms), frac(np)
real(8),intent(out) :: Fsum
real(8) :: cc(3,np)
integer :: a, pp
do a = 1,natoms
   acoor(:,a) = (baseat(:,a) + shift(:,a)/ans2bohr)*ans2bohr
enddo
do pp = 1,np
   cc(:,pp) = coor0(:,pp) + shift(:,atom_of(pp))
enddo
call atom_pair_geom(acoor, Rij, aij)
call becke_fractions(np, cc, atom_of, acoor, Rij, aij, frac)
Fsum = sum(fct*frac)
end subroutine weight_fraction_F

subroutine vv10_nlc_weight_hessian_selftest()
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,parameter :: NPPA = 2
real(8),parameter :: HFDTEST = 1.0d-4
real(8),allocatable :: coor0(:,:), wquad(:), frac(:), fct(:), baseat(:,:)
integer,allocatable :: atom_of(:)
real(8),allocatable :: acoor(:,:), Rij(:,:), aij(:,:), hess(:,:,:,:), hfd(:,:,:,:)
real(8) :: shift(3,natoms), Fp, Fm, F0, Fpp, Fpm, Fmp, Fmm, gfd, xmax, xerr, uni
integer :: np, p, i, sgn, C, Cp, dd, ee, a1, a2, a3, a4
real(8) :: asymfd, asyman, trfd, tran, ae, af
np = natoms*NPPA
allocate(coor0(3,np), atom_of(np), wquad(np), fct(np), frac(np), baseat(3,natoms))
p = 0
do i = 1,natoms
   do sgn = 1,NPPA
      p = p + 1
      coor0(1,p) = atoms(i)%coor(1)*ans2bohr + 0.25d0*sgn
      coor0(2,p) = atoms(i)%coor(2)*ans2bohr - 0.15d0*sgn
      coor0(3,p) = atoms(i)%coor(3)*ans2bohr + 0.10d0*(3-sgn)
      atom_of(p) = i
   enddo
enddo
do i = 1,natoms ; baseat(:,i) = atoms(i)%coor ; enddo
do p = 1,np
   i = atom_of(p)
   fct(p) = exp(-dsqrt(sum((coor0(:,p)-baseat(:,i)*ans2bohr)**2)))
enddo
allocate(acoor(3,natoms), Rij(natoms,natoms), aij(natoms,natoms))
shift = 0.0d0
call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, F0)
wquad = frac
allocate(hess(natoms,3,natoms,3), hfd(natoms,3,natoms,3))
call vv10_nlc_weight_hessian(np, coor0, atom_of, wquad, fct, hess)

shift = 0.0d0
do i = 1,natoms ; shift(1,i) = 0.70d0 ; shift(2,i) = -0.30d0 ; shift(3,i) = 0.20d0 ; enddo
call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, uni)
print '(A,ES12.4)', " VV10 WHESS SELFTEST: |F(uniform)-F(0)| = ", dabs(uni-F0)

hfd = 0.0d0 ; xmax = 0.0d0 ; xerr = 0.0d0
do C = 1,natoms
   do ee = 1,3
      do i = 1,natoms
         do dd = 1,3
            shift = 0.0d0
            if (i .eq. C .and. dd .eq. ee) then
               shift(ee,C) = HFDTEST
               call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, Fp)
               shift = 0.0d0 ; shift(ee,C) = -HFDTEST
               call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, Fm)
               gfd = (Fp - 2.0d0*F0 + Fm)/HFDTEST**2
            else
               shift = 0.0d0 ; shift(dd,i) = HFDTEST ; shift(ee,C) = HFDTEST
               call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, Fpp)
               shift = 0.0d0 ; shift(dd,i) = HFDTEST ; shift(ee,C) = -HFDTEST
               call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, Fpm)
               shift = 0.0d0 ; shift(dd,i) = -HFDTEST ; shift(ee,C) = HFDTEST
               call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, Fmp)
               shift = 0.0d0 ; shift(dd,i) = -HFDTEST ; shift(ee,C) = -HFDTEST
               call weight_fraction_F(np, coor0, atom_of, baseat, shift, fct, acoor, Rij, aij, frac, Fmm)
               gfd = (Fpp - Fpm - Fmp + Fmm)/(4.0d0*HFDTEST**2)
            endif
            hfd(i,dd,C,ee) = gfd
            xerr = max(xerr, dabs(hess(i,dd,C,ee)-gfd))
            xmax = max(xmax, dabs(gfd))
         enddo
      enddo
   enddo
enddo
print '(A,ES12.4)', " VV10 WHESS SELFTEST: max|analytic - FD| = ", xerr
print '(A,ES12.4)', " VV10 WHESS SELFTEST: relative            = ", xerr/max(xmax,1.0d-300)
asymfd = 0.0d0 ; asyman = 0.0d0 ; trfd = 0.0d0 ; tran = 0.0d0
ae = 0.0d0 ; af = 0.0d0 ; a1=0;a2=0;a3=0;a4=0
do i = 1,natoms ; do dd = 1,3 ; do C = 1,natoms ; do ee = 1,3
   asymfd = max(asymfd, dabs(hfd(i,dd,C,ee)-hfd(C,ee,i,dd)))
   asyman = max(asyman, dabs(hess(i,dd,C,ee)-hess(C,ee,i,dd)))
   trfd = max(trfd, dabs(sum(hfd(i,dd,:,:))))
   tran = max(tran, dabs(sum(hess(i,dd,:,:))))
   if (dabs(hess(i,dd,C,ee)-hfd(i,dd,C,ee)) .gt. ae) then
      ae = dabs(hess(i,dd,C,ee)-hfd(i,dd,C,ee))
      af = hfd(i,dd,C,ee)
      a1=i;a2=dd;a3=C;a4=ee
   endif
enddo ; enddo ; enddo ; enddo
print '(A,2ES12.4)', " VV10 WHESS SELFTEST: asym FD / analytic = ", asymfd, asyman
print '(A,2ES12.4)', " VV10 WHESS SELFTEST: transl FD / analytic = ", trfd, tran
print '(A,4I4,2ES16.6)', " VV10 WHESS SELFTEST: worst (i,dd,C,ee) analytic/FD =", a1, a2, a3, a4, &
      hess(a1,a2,a3,a4), af
deallocate(coor0,atom_of,wquad,fct,frac,baseat,acoor,Rij,aij,hess,hfd)
end subroutine vv10_nlc_weight_hessian_selftest

subroutine weight_fraction_dummy(np, coor, atom_of, acoor, Rij, aij, fct, wq)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: np
real(8),intent(in) :: coor(3,*)
integer,intent(in) :: atom_of(*)
real(8),intent(in) :: acoor(3,natoms), Rij(natoms,natoms), aij(natoms,natoms), fct(*)
real(8),intent(out) :: wq(*)
call becke_fractions(np, coor, atom_of, acoor, Rij, aij, wq)
end subroutine weight_fraction_dummy

subroutine spatial_sort_grid(n)
use GRID_info, only: Grid, Grids
implicit none
interface
   recursive subroutine mergesort_morton(key, idx, key_tmp, idx_tmp, lo, hi)
   integer(8),intent(inout) :: key(:)
   integer,intent(inout) :: idx(:)
   integer(8),intent(inout) :: key_tmp(:)
   integer,intent(inout) :: idx_tmp(:)
   integer,intent(in) :: lo,hi
   end subroutine mergesort_morton
end interface
integer,intent(in) :: n
integer(8),allocatable :: morton(:)
integer,allocatable :: idx(:),idx_tmp(:)
integer(8),allocatable :: morton_tmp(:)
type(Grid),allocatable :: Grids_sorted(:)
real(8) :: xmin,xmax,ymin,ymax,zmin,zmax
integer(8) :: qx,qy,qz
integer,parameter :: QBITS = 20
integer(8),parameter :: QMAX = 2_8**QBITS - 1
integer :: i,b

allocate(morton(n),idx(n))
xmin = minval(Grids(1:n)%coor(1)); xmax = maxval(Grids(1:n)%coor(1))
ymin = minval(Grids(1:n)%coor(2)); ymax = maxval(Grids(1:n)%coor(2))
zmin = minval(Grids(1:n)%coor(3)); zmax = maxval(Grids(1:n)%coor(3))
!$omp parallel do private(i,qx,qy,qz,b)
do i = 1,n
   if (xmax .gt. xmin) then
      qx = int((Grids(i)%coor(1)-xmin)/(xmax-xmin) * QMAX, 8)
   else
      qx = 0_8
   endif
   if (ymax .gt. ymin) then
      qy = int((Grids(i)%coor(2)-ymin)/(ymax-ymin) * QMAX, 8)
   else
      qy = 0_8
   endif
   if (zmax .gt. zmin) then
      qz = int((Grids(i)%coor(3)-zmin)/(zmax-zmin) * QMAX, 8)
   else
      qz = 0_8
   endif
   morton(i) = 0_8
   do b = 0,QBITS-1
      if (btest(qx,b)) morton(i) = ibset(morton(i), 3*b)
      if (btest(qy,b)) morton(i) = ibset(morton(i), 3*b+1)
      if (btest(qz,b)) morton(i) = ibset(morton(i), 3*b+2)
   enddo
   idx(i) = i
enddo
!$omp end parallel do

allocate(morton_tmp(n),idx_tmp(n))
call mergesort_morton(morton, idx, morton_tmp, idx_tmp, 1, n)
deallocate(morton_tmp,idx_tmp,morton)

allocate(Grids_sorted(n))
do i = 1,n
   Grids_sorted(i) = Grids(idx(i))
enddo
Grids(1:n) = Grids_sorted(1:n)
deallocate(Grids_sorted,idx)
end subroutine spatial_sort_grid

recursive subroutine mergesort_morton(key, idx, key_tmp, idx_tmp, lo, hi)
implicit none
integer(8),intent(inout) :: key(:)
integer,intent(inout) :: idx(:)
integer(8),intent(inout) :: key_tmp(:)
integer,intent(inout) :: idx_tmp(:)
integer,intent(in) :: lo,hi
integer :: mid,i,j,k

if (hi .le. lo) return
mid = (lo+hi)/2
call mergesort_morton(key, idx, key_tmp, idx_tmp, lo, mid)
call mergesort_morton(key, idx, key_tmp, idx_tmp, mid+1, hi)

i = lo; j = mid+1; k = lo
do while (i .le. mid .and. j .le. hi)
   if (key(i) .le. key(j)) then
      key_tmp(k) = key(i); idx_tmp(k) = idx(i); i = i+1
   else
      key_tmp(k) = key(j); idx_tmp(k) = idx(j); j = j+1
   endif
   k = k+1
enddo
do while (i .le. mid)
   key_tmp(k) = key(i); idx_tmp(k) = idx(i); i = i+1; k = k+1
enddo
do while (j .le. hi)
   key_tmp(k) = key(j); idx_tmp(k) = idx(j); j = j+1; k = k+1
enddo
key(lo:hi) = key_tmp(lo:hi)
idx(lo:hi) = idx_tmp(lo:hi)
end subroutine mergesort_morton

subroutine gridgen_count_points(totGrid, grdRec)
use MOL_info
use GRID_info, only: xcgrid_level, xcgrid_m3_radial, xcgrid_prune_angular_dmin, &
                     xcgrid_dmin_shell, xcgrid_prune_angular, xcgrid_ang_cap, &
                     xcgrid_prune_sphpot_zone, xcgrid_use_zone_table, &
                     xcgrid_zone_calibrated_for, xcgrid_atom_sphpot, &
                     xcgrid_radpot_ladder, xcgrid_radpot_theory, &
                     xcgrid_use_dmin_table, &
                     xcgrid_lebedev_floor_to_supported, &
                     XCGRID_M3_ANGGRID, XCGRID_M3_RSCALE, XCGRID_M3_NRAD, &
                     xcgrid_radpot_scale, xcgrid_radpot_abs, &
                     xcgrid_radpot_rowscale, xcgrid_radpot_rowscale_h, &
                     xcgrid_radpot_rowscale_z, xcgrid_theory_rowscale, &
                     xcgrid_theory_anggrid, xcgrid_dump_shells, &
                     force_dense, force_dense_mgga, XCGRID_FINE, XCGRID_COARSE
implicit none
INCLUDE 'parameter.h'
integer,intent(out) :: totGrid
integer,allocatable,intent(out) :: grdRec(:,:)
integer :: iatm, jatm, radpot, sphpot, i, row, nr_quarter, acap
real(8) :: radr, radw, d_min, ratio
real(8),allocatable :: Rij(:,:)

allocate(Rij(natoms,natoms))
Rij = 0.0d0
do iatm = 1,natoms
   do jatm = 1,natoms
      if (iatm .eq. jatm) cycle
      Rij(iatm,jatm) = dsqrt(sum((atoms(iatm)%coor*ans2bohr - atoms(jatm)%coor*ans2bohr)**2))
   enddo
enddo

totGrid = 0
allocate(grdRec(natoms,2))
do iatm = 1,natoms
   select case (int(atoms(iatm)%charge))
   case (:2);    row = 1
   case (3:10);  row = 2
   case (11:18); row = 3
   case (19:36); row = 4
   case (37:54); row = 5
   case (55:86); row = 6
   case default; row = 7
   end select
   if (xcgrid_radpot_rowscale_z(max(1,min(103,int(atoms(iatm)%charge)))) > 0.0d0) then
      radpot = max(10, nint(XCGRID_M3_NRAD(row) * &
               xcgrid_radpot_rowscale_z(max(1,min(103,int(atoms(iatm)%charge))))))
   else if (xcgrid_theory_rowscale(row) > 0.0d0) then
      radpot = xcgrid_radpot_theory(row)
   else if (row .eq. 1 .and. xcgrid_radpot_rowscale_h > 0.0d0) then
      radpot = max(10, nint(XCGRID_M3_NRAD(row) * xcgrid_radpot_rowscale_h))
   else if (xcgrid_radpot_rowscale > 0.0d0) then
      radpot = max(10, nint(XCGRID_M3_NRAD(row) * xcgrid_radpot_rowscale))
   else if (xcgrid_radpot_abs > 0) then
      radpot = max(10, xcgrid_radpot_abs)
   else
      radpot = xcgrid_radpot_ladder(row)
   endif
   sphpot  = XCGRID_M3_ANGGRID(xcgrid_level)
   if (force_dense) then
      radpot = max(radpot, 99)
      sphpot = max(sphpot, 7)
   endif
   grdRec(iatm,1) = radpot
   grdRec(iatm,2) = sphpot
   nr_quarter = min(row, 3)
   acap = xcgrid_ang_cap(atoms(iatm)%npShell, atoms(iatm)%ndShell, &
                         atoms(iatm)%nfShell, atoms(iatm)%ngShell)
   if (force_dense) acap = 100000
   do i = 1,radpot
      call xcgrid_m3_radial(int(atoms(iatm)%charge), i, radpot, radr, radw)
      ratio = radr / (covrad(atoms(iatm)%charge)*ans2bohr)
      d_min = 0.0d0
      if (xcgrid_use_dmin_table .and. natoms .gt. 1) &
         call xcgrid_dmin_shell(radr, Rij(iatm,:), iatm, natoms, d_min)
      totGrid = totGrid + xcgrid_atom_sphpot(sphpot, xcgrid_theory_anggrid(row), &
                     ratio, nr_quarter, int(atoms(iatm)%charge), acap, &
                     xcgrid_use_dmin_table .and. natoms .gt. 1, d_min)
   enddo
enddo
deallocate(Rij)
end subroutine gridgen_count_points
