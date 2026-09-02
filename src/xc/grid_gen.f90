! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Builds the DFT quadrature grid (Becke partitioning, radial/angular point generation, including NLC/VV10 variants).

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
    integer,allocatable    :: grdRec(:,:)
    real(8),allocatable    :: potx(:),poty(:),potz(:),potw(:)
    real(8)    :: radx,radr,radw,ratio
    real(8)    :: parm
    real(8)    :: radr_pre,radw_pre
    real(8)    :: rmiu,tmps,sij
    integer    :: n_pairs,ip
    integer,allocatable :: pair_i(:),pair_j(:)
    integer,allocatable :: far_j_list(:), far_j_start(:), far_j_count(:)
    integer :: n_far, far_ii, far_jj
    real(8),parameter :: PAIR_CUTOFF = 30.0d0
    real(8),allocatable    :: Rij(:,:),aij_mat(:,:)
    real(8),allocatable    :: rdist(:)
    real(8)    :: chi,uij,aij

    real(8),allocatable    :: Pvec(:)
    integer,allocatable    :: atom_offset(:)
    integer :: atom_npts

   call gridgen_count_points(totGrid, grdRec)
   print '(A,I0)',"  VXC grid points: ",totGrid
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
         chi = covrad(atoms(ii)%charge)/covrad(atoms(jj)%charge)
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
      do i = 1,grdRec(iatm,1)
         call xcgrid_m3_radial(int(atoms(iatm)%charge), i, grdRec(iatm,1), radr_pre, radw_pre)
         ratio = radr_pre / (covrad(atoms(iatm)%charge)*ans2bohr)
         atom_npts = atom_npts + xcgrid_prune_orca(grdRec(iatm,2), ratio, row)
      enddo
   enddo

   !$omp parallel &
   !$omp&   private(iatm,radpot,sphpot,nr_quarter,cursphpot,row,ratio) &
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
      do i = 1,radpot
         call xcgrid_m3_radial(int(atoms(iatm)%charge), i, radpot, radr, radw)
         ratio = radr / (covrad(atoms(iatm)%charge)*ans2bohr)
         cursphpot = xcgrid_prune_orca(sphpot, ratio, row)
         if (cursphpot .eq. 974) call LD0974(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 770) call LD0770(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 590) call LD0590(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 434) call LD0434(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 350) call LD0350(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 302) call LD0302(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 266) call LD0266(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 230) call LD0230(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 194) call LD0194(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 170) call LD0170(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 110) call LD0110(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 86)  call LD0086(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 74)  call LD0074(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 50)  call LD0050(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 26)  call LD0026(potx,poty,potz,potw,Ntemp)
         if (cursphpot .eq. 14)  call LD0014(potx,poty,potz,potw,Ntemp)
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
use GRID_info, only: xcgrid_prune_sphpot
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
   parm = covrad(atoms(iatm)%charge)/2*ans2bohr
   label = atom_off(iatm)
   do i = 1,nrad_atom(iatm)
      cursphpot = xcgrid_prune_sphpot(nsph_atom(iatm), i, nr_quarter_atom(iatm))
      call ld_dispatch(cursphpot, potx,poty,potz,potw, Ntemp)
      radx = cos(i*PI/(nrad_atom(iatm)+1))
      radr = (1+radx)/(1-radx)*parm
      radw = 2*PI/(nrad_atom(iatm)+1)*parm**3*(1+radx)**2.5D0/(1-radx)**3.5D0*4*PI
      do j = 1,Ntemp
         label = label + 1
         coor_out(1,label) = radr*potx(j) + acoor(1,iatm)
         coor_out(2,label) = radr*poty(j) + acoor(2,iatm)
         coor_out(3,label) = radr*potz(j) + acoor(3,iatm)
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
use GRID_info, only: xcgrid_prune_sphpot
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
   parm = covrad(atoms(iatm)%charge)/2*ans2bohr
   do i = 1,nrad
      cursphpot = xcgrid_prune_sphpot(nsph, i, nr_quarter)
      call ld_dispatch_sub(cursphpot, potx,poty,potz,potw, Ntemp)
      radx = cos(i*PI/(nrad+1))
      radr = (1+radx)/(1-radx)*parm
      radw = 2*PI/(nrad+1)*parm**3*(1+radx)**2.5D0/(1-radx)**3.5D0*4*PI
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

subroutine vv10_nlc_weight_gradient(np, coor, atom_of, nper, rho, eps, wquad, beta_v, gv_wt)
use MOL_info
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: np, nper
real(8),intent(in) :: coor(3,np)
integer,intent(in) :: atom_of(np)
real(8),intent(in) :: rho(np), eps(np), wquad(np), beta_v
real(8),intent(out) :: gv_wt(natoms,3)

real(8) :: Rij(natoms,natoms), aij_mat(natoms,natoms)
real(8) :: acoor(3,natoms)
real(8) :: chi, uij_v, aij
integer :: i2, jj2
real(8),parameter :: SAFE_S_FLOOR = 1.0d-280

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
      uij_v = (chi-1.0d0)/(chi+1.0d0)
      aij = uij_v/(uij_v**2-1.0d0)
      if (aij >  0.5d0) aij =  0.5d0
      if (aij < -0.5d0) aij = -0.5d0
      aij_mat(i2,jj2) = aij
   enddo
enddo

gv_wt = 0.0d0

!$omp parallel do schedule(dynamic) reduction(+:gv_wt)
do i2 = 1,np
   block
      integer :: O, A, B, C, dd
      real(8) :: rdist(natoms), uvec(3,natoms)
      real(8) :: mu(natoms,natoms), x1(natoms,natoms), x2(natoms,natoms), x3(natoms,natoms), sAB(natoms,natoms)
      real(8) :: Pvec(natoms), sumP
      real(8) :: dPA_sum(natoms), dsumP, dPO
      real(8) :: drA, drB, dRAB, dmu, dmuhat, dsAB_val
      real(8) :: factor

      O = atom_of(i2)
      do A = 1,natoms
         rdist(A) = dsqrt(sum((coor(:,i2)-acoor(:,A))**2))
         uvec(:,A) = (coor(:,i2)-acoor(:,A))/rdist(A)
      enddo
      do A = 1,natoms
         do B = 1,natoms
            if (A .eq. B) cycle
            mu(A,B) = (rdist(A)-rdist(B))/Rij(A,B)
            x1(A,B) = mu(A,B) + aij_mat(A,B)*(1.0d0-mu(A,B)**2)
            x2(A,B) = 1.5d0*x1(A,B) - 0.5d0*x1(A,B)**3
            x3(A,B) = 1.5d0*x2(A,B) - 0.5d0*x2(A,B)**3
            sAB(A,B) = 0.5d0*(1.0d0 - (1.5d0*x3(A,B)-0.5d0*x3(A,B)**3))
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
      factor = rho(i2)*(2.0d0*eps(i2)-beta_v)
      if (dabs(wquad(i2)) .lt. SAFE_S_FLOOR) cycle

      do C = 1,natoms
         do dd = 1,3
            dPA_sum = 0.0d0
            if (C .eq. O) then
               do A = 1,natoms
                  do B = 1,natoms
                     if (A .eq. B) cycle
                     if (sAB(A,B) .lt. SAFE_S_FLOOR) cycle
                     if (A .ne. C) then; drA = uvec(dd,A); else; drA = 0.0d0; endif
                     if (B .ne. C) then; drB = uvec(dd,B); else; drB = 0.0d0; endif
                     dRAB = 0.0d0
                     if (A .eq. C) dRAB = dRAB + (acoor(dd,A)-acoor(dd,B))/Rij(A,B)
                     if (B .eq. C) dRAB = dRAB - (acoor(dd,A)-acoor(dd,B))/Rij(A,B)
                     dmu = (drA-drB)/Rij(A,B) - mu(A,B)/Rij(A,B)*dRAB
                     dmuhat = dmu*(1.0d0 - 2.0d0*aij_mat(A,B)*mu(A,B))
                     dsAB_val = -1.6875d0*(1.0d0-x1(A,B)**2)*(1.0d0-x2(A,B)**2)*(1.0d0-x3(A,B)**2)*dmuhat
                     dPA_sum(A) = dPA_sum(A) + dsAB_val/sAB(A,B)
                  enddo
               enddo
               dPA_sum = dPA_sum*Pvec
            else
               do A = 1,natoms
                  if (A .eq. C) cycle
                  if (sAB(A,C) .lt. SAFE_S_FLOOR) cycle
                  drB = -uvec(dd,C)
                  dRAB = -(acoor(dd,A)-acoor(dd,C))/Rij(A,C)
                  dmu = (0.0d0-drB)/Rij(A,C) - mu(A,C)/Rij(A,C)*dRAB
                  dmuhat = dmu*(1.0d0 - 2.0d0*aij_mat(A,C)*mu(A,C))
                  dsAB_val = -1.6875d0*(1.0d0-x1(A,C)**2)*(1.0d0-x2(A,C)**2)*(1.0d0-x3(A,C)**2)*dmuhat
                  dPA_sum(A) = Pvec(A)*dsAB_val/sAB(A,C)
               enddo
               block
                  real(8) :: dsum_c, drC2, dRCB, dmu3, dmuhat3, dsCB
                  dsum_c = 0.0d0
                  drC2 = -uvec(dd,C)
                  do B = 1,natoms
                     if (B .eq. C) cycle
                     if (sAB(C,B) .lt. SAFE_S_FLOOR) cycle
                     dRCB = (acoor(dd,C)-acoor(dd,B))/Rij(C,B)
                     dmu3 = (drC2-0.0d0)/Rij(C,B) - mu(C,B)/Rij(C,B)*dRCB
                     dmuhat3 = dmu3*(1.0d0 - 2.0d0*aij_mat(C,B)*mu(C,B))
                     dsCB = -1.6875d0*(1.0d0-x1(C,B)**2)*(1.0d0-x2(C,B)**2)*(1.0d0-x3(C,B)**2)*dmuhat3
                     dsum_c = dsum_c + dsCB/sAB(C,B)
                  enddo
                  dPA_sum(C) = Pvec(C)*dsum_c
               end block
            endif
            dsumP = sum(dPA_sum)
            dPO = dPA_sum(O)
            gv_wt(C,dd) = gv_wt(C,dd) + factor*wquad(i2)*(dPO/Pvec(O) - dsumP/sumP)
         enddo
      enddo
   end block
enddo
!$omp end parallel do

end subroutine vv10_nlc_weight_gradient

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
use GRID_info, only: xcgrid_level, xcgrid_m3_radial, xcgrid_prune_orca, &
                     XCGRID_M3_ANGGRID, XCGRID_M3_RSCALE, XCGRID_M3_NRAD, &
                     force_dense, force_dense_mgga, XCGRID_FINE, XCGRID_COARSE
implicit none
INCLUDE 'parameter.h'
integer,intent(out) :: totGrid
integer,allocatable,intent(out) :: grdRec(:,:)
integer :: iatm, radpot, sphpot, i, row, nr_quarter
real(8) :: ratio, radr, radw

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
   radpot  = max(10, nint(XCGRID_M3_NRAD(row) * XCGRID_M3_RSCALE(xcgrid_level)))
   sphpot  = XCGRID_M3_ANGGRID(xcgrid_level)
   if (force_dense) then
      radpot = max(radpot, 99)
      sphpot = max(sphpot, 7)
   endif
   grdRec(iatm,1) = radpot
   grdRec(iatm,2) = sphpot
   nr_quarter = min(row, 3)
   do i = 1,radpot
      call xcgrid_m3_radial(int(atoms(iatm)%charge), i, radpot, radr, radw)
      ratio = radr / (covrad(atoms(iatm)%charge)*ans2bohr)
      totGrid = totGrid + xcgrid_prune_orca(sphpot, ratio, nr_quarter)
   enddo
enddo
end subroutine gridgen_count_points
