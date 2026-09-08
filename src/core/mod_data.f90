! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! MOL_info (molecule/basis global state) and GRID_info (DFT quadrature grid global state) - Engine's two core shared-state modules.

module MOL_info
implicit none
type shells
    integer    :: angMoment
    integer    :: nGauss
    real(8),allocatable    :: exponents(:)
    real(8),allocatable    :: contrCoeff(:)
    real(8),allocatable    :: cnVal(:,:)
end type shells

type atom
    integer    :: charge
    real(8)    :: coor(3)
    character  :: base*30
    integer    :: nsShell
    integer    :: npShell
    integer    :: ndShell
    integer    :: nfShell
    integer    :: ngShell
    integer    :: nShell
    integer    :: nconts
    real(8)    :: atmForce(3)
    integer    :: ecpCoreElec = 0
    character  :: ecpbase*30
    type(shells),allocatable :: shell(:)
end type atom
type(atom),allocatable :: atoms(:)

integer Natoms,              Nconts ,            ncontssph

integer :: nPointCharges = 0
real(8),allocatable :: pointcharge_q(:)
real(8),allocatable :: pointcharge_coor(:,:)

real(8),allocatable  :: Fa(:,:),  Fb(:,:), S(:,:),   Hcore(:,:)

real(8),allocatable  ::  Pa(:,:),    Pb(:,:)

real(8),allocatable  :: C_a(:,:),   C_b(:,:)

real(8),allocatable  :: eLev_a(:),eLev_b(:)

real(8),allocatable  :: X(:,:)

integer              :: Charge
integer              :: Multi
logical               :: resp_charges_on = .false.
real(8)              :: E,E_rep

integer              :: n_alpha
integer              :: n_beta
integer              :: coreChg

integer,allocatable  :: linkMat(:,:)

character Functional*30

integer :: engine_verbose = 1

end module MOL_info

module GRID_info
implicit none
type Grid
    real(8)    :: coor(3)
    real(8)    :: weight
    real(8)    :: rho_a,rho_b
    real(8)    :: sigma_aa,sigma_ab,sigma_bb
    real(8),pointer :: TempD(:) => null()
end type Grid

integer :: ngrids

type(Grid),allocatable :: Grids(:)

real(8),allocatable,target :: TempD_all(:,:)

type GridValBlock
   real(4),allocatable :: val0(:,:)
   real(4),allocatable :: val1(:,:,:)
end type GridValBlock
type(GridValBlock),allocatable :: val_blocks(:)

logical :: grid_cache_mode = .true.
integer :: grid_cache_capacity = 0
real(8),allocatable :: rcut2_shared(:,:)
integer :: maxShell_shared
integer :: maxGauss_shared

integer,parameter :: CACHE_BATCH = 128
integer :: n_cache_batches = 0
integer :: max_batch_nsig = 0
integer,allocatable :: batch_nsig(:)
integer,allocatable :: batch_sig_idx(:,:)

type BatchShellList
   integer :: n_shells = 0
   integer :: n_sig = 0
   integer,allocatable :: sh_atom(:)
   integer,allocatable :: sh_idx(:)
   integer,allocatable :: sh_off(:)
   integer,allocatable :: sh_ndim(:)
   integer,allocatable :: sig_idx(:)
end type BatchShellList
integer,parameter :: DFT_BATCH = 64
real(8),parameter :: XC_SIG_TIGHT = 1.0d-7
real(8),parameter :: XC_SIG_LOOSE = 1.0d-5
real(8) :: xc_sig_now = XC_SIG_TIGHT

integer :: n_dft_batches = 0
integer :: max_dft_batch_nsig = 0
type(BatchShellList),allocatable,target :: dft_batch_shells(:)
logical :: dft_batch_shells_built = .false.
integer,parameter :: XC_NLEVEL = 3
real(8),parameter :: XC_LEVEL_CUT(XC_NLEVEL) = (/1.0d-3, 1.0d-5, 1.0d-7/)
real(8),parameter :: XC_LEVEL_PRMS(XC_NLEVEL) = (/1.0d-3, 1.0d-4, 0.0d0/)
type(BatchShellList),allocatable,target :: dft_batch_shells_lvl(:,:)
integer :: max_dft_batch_nsig_loose = 0
integer :: xc_geo_level = XC_NLEVEL
integer,allocatable :: gridgen_nlc_atom_of(:)
logical :: dft_loose_built = .false.

integer,parameter :: XCGRID_COARSE = 2
integer,parameter :: XCGRID_FINE   = 3
integer,parameter :: XCGRID_L4     = 4
integer,parameter :: XCGRID_L5     = 5
integer,parameter :: XCGRID_L6     = 6
integer,parameter :: XCGRID_L7     = 7
real(8) :: XCGRID_RSCALE(7)  = (/0.0d0, 0.60d0, 1.848d0, 0.845d0, 2.6d0, 1.30d0, 4.0d0/)
integer :: XCGRID_M3_ANGGRID(7) = (/0, 3, 4, 5, 6, 7, 7/)
real(8) :: XCGRID_M3_RSCALE(7)  = (/0.0d0, 0.65d0, 0.95d0, 1.55d0, 2.10d0, 2.70d0, 3.40d0/)
integer :: XCGRID_M3_NRAD(7)    = (/45, 55, 60, 65, 70, 75, 80/)
integer :: XCGRID_SPH_IN(7)  = (/0, 302, 302, 974, 974, 974, 974/)
integer :: XCGRID_SPH_EDGE(7) = (/0, 302, 302, 590, 590, 590, 974/)
integer :: xcgrid_level = XCGRID_FINE
integer :: xcgrid_active_coarse = XCGRID_COARSE
integer :: xcgrid_active_fine   = XCGRID_FINE
logical :: xcgrid_dynamic = .false.
logical :: xcgrid_refined = .false.
logical :: force_dense = .false.
logical :: force_dense_mgga = .false.
integer,parameter :: XC_DIRECT_NCONTS = 500
logical :: xc_direct_mode = .false.
logical :: xcgrid_just_refined = .false.
real(8) :: xcgrid_switch_prms = 3.0d-4

type ShellBatchList
   integer :: n = 0
   integer,allocatable :: batch(:)
   integer,allocatable :: off(:)
end type ShellBatchList
type(ShellBatchList),allocatable :: shell_batches(:)
integer,allocatable :: sh_in_batch(:,:)
logical :: shell_tables_built = .false.
integer :: vxc_nshell_tot = 0
type PhimaxVec
   real(8),allocatable :: v(:)
end type PhimaxVec
type(PhimaxVec),allocatable :: batch_phimax(:)
real(8),allocatable :: atpair_dp(:,:)
integer,allocatable :: ao_atom(:)

logical :: xc_incr_on = .false.
logical :: xc_incr_primed = .false.
real(8),allocatable :: xcr_Pa(:,:),xcr_rho(:),xcr_grad(:,:)
real(8),allocatable :: xcr_Fxc(:,:)
real(8) :: xcr_Exc = 0.0d0
real(8),allocatable :: atpair_ddp(:,:)
integer,allocatable :: batch_atoms(:,:),batch_natom(:)
integer,allocatable :: vxc_sh_atom(:),vxc_sh_local(:),vxc_sh_ao0(:),vxc_sh_ndim(:)

contains

integer function xcgrid_m3_row(Zin) result(row)
implicit none
integer, intent(in) :: Zin
select case (max(1, Zin))
case (:2);    row = 1
case (3:10);  row = 2
case (11:18); row = 3
case (19:36); row = 4
case (37:54); row = 5
case (55:86); row = 6
case default; row = 7
end select
end function xcgrid_m3_row

real(8) function xcgrid_m3_xi(Zin) result(xi)
implicit none
integer, intent(in) :: Zin
integer :: z
real(8), parameter :: BR(103) = (/ &
   0.35d0, 1.40d0, &
   1.45d0, 1.05d0, 0.85d0, 0.70d0, 0.65d0, 0.60d0, 0.50d0, 1.50d0, &
   1.80d0, 1.50d0, 1.25d0, 1.10d0, 1.00d0, 1.00d0, 1.00d0, 1.80d0, &
   2.20d0, 1.80d0, &
   1.60d0, 1.40d0, 1.35d0, 1.40d0, 1.40d0, 1.40d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0, &
   1.30d0, 1.25d0, 1.15d0, 1.15d0, 1.15d0, 1.90d0, &
   2.35d0, 2.00d0, &
   1.80d0, 1.55d0, 1.45d0, 1.45d0, 1.35d0, 1.30d0, 1.35d0, 1.40d0, 1.60d0, 1.55d0, &
   1.55d0, 1.45d0, 1.45d0, 1.40d0, 1.40d0, 2.10d0, &
   2.60d0, 2.15d0, &
   1.95d0, 1.85d0, 1.85d0, 1.85d0, 1.85d0, 1.85d0, 1.85d0, &
   1.80d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, &
   1.55d0, 1.45d0, 1.35d0, 1.35d0, 1.30d0, 1.35d0, 1.35d0, 1.35d0, 1.50d0, &
   1.90d0, 1.80d0, 1.60d0, 1.90d0, 1.45d0, 2.10d0, &
   1.80d0, 2.15d0, &
   1.95d0, 1.80d0, 1.80d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, &
   1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0 /)
real(8), parameter :: XCGRID_M3_KROW(7) = &
   (/ 1.80d0, 1.40d0, 0.90d0, 0.75d0, 0.70d0, 0.65d0, 0.65d0 /)
z = max(1, min(103, Zin))
xi = XCGRID_M3_KROW(xcgrid_m3_row(z)) * BR(z) * 1.88972612456506d0
end function xcgrid_m3_xi

subroutine xcgrid_m3_radial(Zin, i, nrad, radr, radw)
implicit none
integer, intent(in) :: Zin, i, nrad
real(8), intent(out) :: radr, radw
real(8) :: pig, ln2, theta, x, s, xi, drdx
pig  = acos(-1.0d0)
ln2  = log(2.0d0)
xi   = xcgrid_m3_xi(Zin)
theta = i * pig / (nrad + 1)
x = cos(theta)
s = sin(theta)
radr = (xi/ln2) * log(2.0d0/(1.0d0 - x))
drdx = (xi/ln2) / (1.0d0 - x)
radw = 4.0d0*pig * (pig/(nrad+1)) * s * radr*radr * drdx
end subroutine xcgrid_m3_radial

integer function xcgrid_prune_sphpot(sphpot, i, nr_quarter) result(cursphpot)
implicit none
integer,intent(in) :: sphpot, i, nr_quarter
if (i .le. 2*nr_quarter) then
   cursphpot = sphpot
else if (i .le. 3*nr_quarter) then
   select case (sphpot)
   case (434) ; cursphpot = 302
   case (302) ; cursphpot = 170
   case (230) ; cursphpot = 110
   case (170) ; cursphpot = 74
   case default ; cursphpot = sphpot
   end select
else
   select case (sphpot)
   case (434) ; cursphpot = 110
   case (302) ; cursphpot = 74
   case (230) ; cursphpot = 26
   case (170) ; cursphpot = 26
   case default ; cursphpot = 26
   end select
endif
end function xcgrid_prune_sphpot

integer function xcgrid_prune_sphpot_pyscf(n_ang, ratio, row) result(cursphpot)
implicit none
integer,intent(in) :: n_ang, row
real(8),intent(in) :: ratio
real(8) :: alphas(4)
integer :: place
integer :: pat302(0:4), pat434(0:4)
pat302 = (/50,86,266,302,266/)
pat434 = (/50,86,350,434,350/)
select case (row)
case (1); alphas = (/0.25d0, 0.5d0, 1.0d0, 4.5d0/)
case (2); alphas = (/0.1667d0, 0.5d0, 0.9d0, 3.5d0/)
case default; alphas = (/0.1d0, 0.4d0, 0.8d0, 2.5d0/)
end select
place = count(ratio .gt. alphas)
if (n_ang .eq. 302) then
   cursphpot = pat302(place)
else
   cursphpot = pat434(place)
endif
end function xcgrid_prune_sphpot_pyscf

integer function xcgrid_prune_orca(anggrid, ratio, row) result(cursphpot)
implicit none
integer,intent(in) :: anggrid, row
real(8),intent(in) :: ratio
integer :: place, ag
real(8) :: alphas(4)
integer,parameter :: ORCA_ANG(5,7) = reshape((/ &
    14,  26,  50,  50,  26,   &
    14,  26,  50, 110,  50,   &
    26,  50, 110, 194, 110,   &
    26, 110, 194, 302, 194,   &
    26, 194, 302, 434, 302,   &
    50, 302, 434, 590, 434,   &
   110, 434, 590, 770, 590 /), (/5,7/))
select case (row)
case (1); alphas = (/0.25d0,   0.5d0, 1.0d0, 4.5d0/)
case (2); alphas = (/0.1667d0, 0.5d0, 0.9d0, 3.5d0/)
case default; alphas = (/0.1d0, 0.4d0, 0.8d0, 2.5d0/)
end select
ag = max(1, min(7, anggrid))
place = count(ratio .gt. alphas)
cursphpot = ORCA_ANG(place+1, ag)
end function xcgrid_prune_orca

integer function xcgrid_fine_level_for_functional(functional) result(lvl)
implicit none
character(len=*),intent(in) :: functional
lvl = XCGRID_FINE
select case (trim(functional))
case ("R2SCAN","R2SCAN0","R2SCAN-3C","R2SCAN_3C", &
      "M06L","M06-L","M06_L","M06","M06_HYB","M06-2X","M06-2X_HYB","M05-2X", &
      "MN15","MN15L","MN15-L")
   lvl = XCGRID_L5
end select
end function xcgrid_fine_level_for_functional

end module GRID_info

