! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! RESP (Restrained ElectroStatic Potential) atomic charges, implemented to

submodule (mod_integrals) resp_impl
implicit none

integer,parameter :: RESP_MAXDEG = 8
real(8),parameter :: RESP_A1 = 0.0005d0
real(8),parameter :: RESP_A2 = 0.001d0
real(8),parameter :: RESP_B  = 0.1d0
real(8),parameter :: RESP_TOLER = 1.0d-5
integer,parameter :: RESP_MAXIT = 25
real(8),parameter :: RESP_DENSITY = 1.0d0
real(8),parameter :: RESP_SCALES(4) = (/1.4d0,1.6d0,1.8d0,2.0d0/)

contains

module subroutine resp_charges(Natoms_in, RESPcharge_out)
use MOL_info
use mod_profile, only: prof_start, prof_stop
implicit none
INCLUDE 'parameter.h'
integer,intent(in) :: Natoms_in
real(8),intent(out) :: RESPcharge_out(Natoms_in)

real(8),allocatable :: grid(:,:), Vgrid(:), Dmat(:,:), sphere(:,:)
real(8),allocatable :: Ptot(:,:), Bk(:,:)
real(8),allocatable :: con_coef(:,:), con_rhs(:)
real(8) :: q_stage1(Natoms_in), q_final(Natoms_in)
real(8) :: rad_scaled(Natoms_in), p_ang(3), dist, r_ang
integer :: deg(Natoms_in), nbr(RESP_MAXDEG,Natoms_in)
integer :: stage2_group(Natoms_in)
integer :: ia, ib, ilayer, ip, ntot_cap, npts, nsph, ncon, nsph_want
character(64) :: espmode

ntot_cap = 0
do ia = 1,Natoms_in
   do ilayer = 1,4
      r_ang = resp_vdw_radius_ang(atoms(ia)%charge)*RESP_SCALES(ilayer)
      ntot_cap = ntot_cap + int(RESP_DENSITY*4.0d0*acos(-1.0d0)*r_ang*r_ang) + 4
   enddo
enddo
allocate(grid(3,ntot_cap), sphere(3,ntot_cap))

call prof_start("resp_grid_and_esp")
npts = 0
do ilayer = 1,4
   do ia = 1,Natoms_in
      rad_scaled(ia) = resp_vdw_radius_ang(atoms(ia)%charge)*RESP_SCALES(ilayer)
   enddo
   do ia = 1,Natoms_in
      r_ang = rad_scaled(ia)
      nsph_want = int(RESP_DENSITY*4.0d0*acos(-1.0d0)*r_ang*r_ang)
      call resp_gamess_surface(nsph_want, sphere, nsph)
      do ip = 1,nsph
         p_ang = atoms(ia)%coor + r_ang*sphere(:,ip)
         do ib = 1,Natoms_in
            if (ib .eq. ia) cycle
            dist = sqrt(sum((p_ang - atoms(ib)%coor)**2))
            if (dist .lt. rad_scaled(ib)) goto 100
         enddo
         npts = npts + 1
         grid(:,npts) = p_ang*ans2bohr
100      continue
      enddo
   enddo
enddo
deallocate(sphere)

allocate(Ptot(nConts,nConts))
Ptot = Pa + Pb
allocate(Vgrid(npts))

call get_environment_variable("ENGINE_RESP_ESP_MODE", espmode)
if (trim(espmode) .eq. 'legacy') then
   !$omp parallel private(Bk)
   allocate(Bk(nConts,nConts))
   !$omp do schedule(dynamic)
   do ip = 1,npts
      call cosmo_build_one_tess_matrix(grid(:,ip), Bk)
      Vgrid(ip) = sum(Ptot*Bk)
   enddo
   !$omp end do
   deallocate(Bk)
   !$omp end parallel
else
   call esp_at_grid_batch(npts, grid(:,1:npts), Ptot, Vgrid)
endif
deallocate(Ptot)

!$omp parallel do private(ia,dist) schedule(static)
do ip = 1,npts
   do ia = 1,Natoms_in
      dist = sqrt(sum((grid(:,ip) - atoms(ia)%coor*ans2bohr)**2))
      Vgrid(ip) = Vgrid(ip) + real(atoms(ia)%charge - atoms(ia)%ecpCoreElec,8)/dist
   enddo
enddo
!$omp end parallel do

block
character(256) :: dumpenv
integer :: iu
call get_environment_variable("ENGINE_RESP_DUMP", dumpenv)
if (len_trim(dumpenv) > 0) then
   open(newunit=iu, file=trim(dumpenv)//'.grid', status='replace', action='write')
   do ip = 1,npts
      write(iu,'(3F20.12)') grid(:,ip)/ans2bohr
   enddo
   close(iu)
   open(newunit=iu, file=trim(dumpenv)//'.esp', status='replace', action='write')
   do ip = 1,npts
      write(iu,'(F20.12)') Vgrid(ip)
   enddo
   close(iu)
endif
end block

allocate(Dmat(npts,Natoms_in))
do ip = 1,npts
   do ia = 1,Natoms_in
      Dmat(ip,ia) = 1.0d0/sqrt(sum((grid(:,ip) - atoms(ia)%coor*ans2bohr)**2))
   enddo
enddo
call prof_stop("resp_grid_and_esp")

call prof_start("resp_two_stage_fit")

allocate(con_coef(Natoms_in,2*Natoms_in), con_rhs(2*Natoms_in))

ncon = 0
call resp_fit(Natoms_in, npts, Dmat, Vgrid, real(Charge,8), &
              RESP_A1, RESP_B, ncon, con_coef, con_rhs, q_stage1)

call resp_build_bonds(Natoms_in, deg, nbr)
call resp_stage2_constraints(Natoms_in, deg, nbr, q_stage1, ncon, con_coef, con_rhs, stage2_group)
if (ncon .gt. 0) then
   call resp_fit(Natoms_in, npts, Dmat, Vgrid, real(Charge,8), &
                 RESP_A2, RESP_B, ncon, con_coef, con_rhs, q_final)
else
   q_final = q_stage1
endif

call prof_stop("resp_two_stage_fit")

RESPcharge_out = q_final

print *
print *,"===== RESP fit (Psi4 cdsgroup/resp protocol, Bayly:93) ====="
print '(" Grid points kept: ",I0," (of ",I0," candidates); density=",F4.1,"/Ang^2")', &
      npts, ntot_cap, RESP_DENSITY
print '(" VDW scale factors: ",4F6.2,"   radii: GAMESS")', RESP_SCALES
print '(" resp_a: stage1=",F7.4,"  stage2=",F7.4,"   resp_b=",F5.2,"   IHFREE=T")', &
      RESP_A1, RESP_A2, RESP_B
print '(" Equivalence groups: ",I0,"   Total charge: ",F12.8," (target ",F12.8,")")', &
      maxval(stage2_group), sum(q_final), real(Charge,8)
print *,"=============================================================="

deallocate(grid, Vgrid, Dmat, con_coef, con_rhs)

end subroutine resp_charges

function resp_vdw_radius_ang(Z) result(r)
implicit none
integer,intent(in) :: Z
real(8) :: r
select case (Z)
case (1);  r = 1.20d0
case (2);  r = 1.20d0
case (3);  r = 1.37d0
case (4);  r = 1.45d0
case (5);  r = 1.45d0
case (6);  r = 1.50d0
case (7);  r = 1.50d0
case (8);  r = 1.40d0
case (9);  r = 1.35d0
case (10); r = 1.30d0
case (11); r = 1.57d0
case (12); r = 1.36d0
case (13); r = 1.24d0
case (14); r = 1.17d0
case (15); r = 1.80d0
case (16); r = 1.75d0
case (17); r = 1.70d0
case default
   r = 1.80d0
end select
end function resp_vdw_radius_ang

subroutine resp_gamess_surface(n_want, u, nsph_out)
implicit none
integer,intent(in) :: n_want
real(8),intent(out) :: u(:,:)
integer,intent(out) :: nsph_out
real(8),parameter :: EPS = 1.0d-10
real(8) :: pi, fi, fj, z, xy
integer :: nequat, nvert, nhor, i, j

pi = acos(-1.0d0)
nsph_out = 0
if (n_want .le. 0) return
nequat = int(sqrt(pi*real(n_want,8)))
nvert = nequat/2
if (nvert .lt. 1) nvert = 1
do i = 0,nvert
   fi = pi*real(i,8)/real(nvert,8)
   z = cos(fi)
   xy = sin(fi)
   nhor = int(real(nequat,8)*xy + EPS)
   if (nhor .lt. 1) nhor = 1
   do j = 0,nhor-1
      fj = 2.0d0*pi*real(j,8)/real(nhor,8)
      if (nsph_out .ge. n_want) return
      nsph_out = nsph_out + 1
      u(1,nsph_out) = cos(fj)*xy
      u(2,nsph_out) = sin(fj)*xy
      u(3,nsph_out) = z
   enddo
enddo
end subroutine resp_gamess_surface

subroutine resp_fit(natoms, npts, invr, V, mol_charge, resp_a, resp_b, &
                    ncon, con_coef, con_rhs, q_out)
use MOL_info, only: atoms
implicit none
integer,intent(in) :: natoms, npts, ncon
real(8),intent(in) :: invr(npts,natoms), V(npts), mol_charge, resp_a, resp_b
real(8),intent(in) :: con_coef(natoms,*), con_rhs(*)
real(8),intent(out) :: q_out(natoms)

real(8),allocatable :: A0(:,:), A(:,:), B(:), rhs(:)
integer,allocatable :: ipiv(:)
integer :: ndim, i, k, iter, info
real(8) :: q(natoms), qlast(natoms), dif

ndim = natoms + 1 + ncon
allocate(A0(ndim,ndim), A(ndim,ndim), B(ndim), rhs(ndim), ipiv(ndim))

A0 = 0.0d0
B = 0.0d0
A0(1:natoms,1:natoms) = matmul(transpose(invr),invr)
B(1:natoms) = matmul(transpose(invr),V)

A0(1:natoms,natoms+1) = 1.0d0
A0(natoms+1,1:natoms) = 1.0d0
B(natoms+1) = mol_charge

do k = 1,ncon
   B(natoms+1+k) = con_rhs(k)
   do i = 1,natoms
      A0(natoms+1+k,i) = con_coef(i,k)
      A0(i,natoms+1+k) = con_coef(i,k)
   enddo
enddo

A = A0
rhs = B
call dgesv(ndim, 1, A, ndim, ipiv, rhs, ndim, info)
q = rhs(1:natoms)

if (resp_a .gt. 0.0d0) then
   qlast = q
   do iter = 1,RESP_MAXIT
      A = A0
      do i = 1,natoms
         if (atoms(i)%charge .ne. 1) &
            A(i,i) = A0(i,i) + resp_a/sqrt(q(i)*q(i) + resp_b*resp_b)
      enddo
      rhs = B
      call dgesv(ndim, 1, A, ndim, ipiv, rhs, ndim, info)
      q = rhs(1:natoms)
      dif = maxval(abs(q - qlast))
      qlast = q
      if (dif .le. RESP_TOLER) exit
   enddo
endif

q_out = q
deallocate(A0, A, B, rhs, ipiv)
end subroutine resp_fit

subroutine resp_build_bonds(Natoms_in, deg, nbr)
use MOL_info, only: atoms
implicit none
integer,intent(in) :: Natoms_in
integer,intent(out) :: deg(Natoms_in), nbr(RESP_MAXDEG,Natoms_in)
real(8),parameter :: TOL = 1.3d0
integer :: ia, ib
real(8) :: rcut, dist
deg = 0
do ia = 1,Natoms_in
   do ib = ia+1,Natoms_in
      rcut = TOL*(resp_cov_radius_ang(atoms(ia)%charge) + resp_cov_radius_ang(atoms(ib)%charge))
      dist = sqrt(sum((atoms(ia)%coor - atoms(ib)%coor)**2))
      if (dist .lt. rcut) then
         if (deg(ia) .lt. RESP_MAXDEG) then
            deg(ia) = deg(ia) + 1
            nbr(deg(ia),ia) = ib
         endif
         if (deg(ib) .lt. RESP_MAXDEG) then
            deg(ib) = deg(ib) + 1
            nbr(deg(ib),ib) = ia
         endif
      endif
   enddo
enddo
end subroutine resp_build_bonds

function resp_cov_radius_ang(Z) result(r)
implicit none
integer,intent(in) :: Z
real(8) :: r
real(8),parameter :: tbl(36) = (/ &
   0.31d0,0.28d0, 1.28d0,0.96d0,0.84d0,0.76d0,0.71d0,0.66d0,0.57d0,0.58d0, &
   1.66d0,1.41d0,1.21d0,1.11d0,1.07d0,1.05d0,1.02d0,1.06d0, &
   2.03d0,1.76d0,1.70d0,1.60d0,1.53d0,1.39d0,1.39d0,1.32d0,1.26d0,1.24d0, &
   1.32d0,1.22d0,1.22d0,1.20d0,1.19d0,1.20d0,1.20d0,1.16d0 /)
if (Z .ge. 1 .and. Z .le. 36) then
   r = tbl(Z)
else if (Z .eq. 53) then
   r = 1.39d0
else
   r = 1.50d0
endif
end function resp_cov_radius_ang

subroutine resp_stage2_constraints(Natoms_in, deg, nbr, q_stage1, ncon, con_coef, con_rhs, stage2_group)
use MOL_info, only: atoms
implicit none
integer,intent(in) :: Natoms_in, deg(Natoms_in), nbr(RESP_MAXDEG,Natoms_in)
real(8),intent(in) :: q_stage1(Natoms_in)
integer,intent(out) :: ncon
real(8),intent(out) :: con_coef(Natoms_in,*), con_rhs(*)
integer,intent(out) :: stage2_group(Natoms_in)

logical :: refit(Natoms_in)
integer :: hlist(RESP_MAXDEG), nh
integer :: ia, k, j, ig

refit = .false.
stage2_group = 0
ncon = 0
ig = 0

do ia = 1,Natoms_in
   if (atoms(ia)%charge .ne. 6) cycle
   if (deg(ia) .ne. 4) cycle
   nh = 0
   do k = 1,deg(ia)
      if (atoms(nbr(k,ia))%charge .eq. 1) then
         nh = nh + 1
         hlist(nh) = nbr(k,ia)
      endif
   enddo
   if (nh .eq. 0) cycle
   ig = ig + 1
   refit(ia) = .true.
   stage2_group(ia) = -1
   do k = 1,nh
      refit(hlist(k)) = .true.
      stage2_group(hlist(k)) = ig
   enddo
   do k = 2,nh
      ncon = ncon + 1
      con_rhs(ncon) = 0.0d0
      do j = 1,Natoms_in
         con_coef(j,ncon) = 0.0d0
      enddo
      con_coef(hlist(k-1),ncon) = -1.0d0
      con_coef(hlist(k),ncon)   =  1.0d0
   enddo
enddo

if (ig .eq. 0) then
   ncon = 0
   return
endif

do ia = 1,Natoms_in
   if (refit(ia)) cycle
   ncon = ncon + 1
   con_rhs(ncon) = q_stage1(ia)
   do j = 1,Natoms_in
      con_coef(j,ncon) = 0.0d0
   enddo
   con_coef(ia,ncon) = 1.0d0
enddo
end subroutine resp_stage2_constraints

end submodule resp_impl
