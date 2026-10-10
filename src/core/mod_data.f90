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

logical :: engine_quiet = .false.

end module MOL_info

module GRID_info
use xcgrid_defaults
implicit none
real(8), allocatable, save :: xc_dbg_pt(:,:)
logical, save :: xc_dbg_checked = .false., xc_dbg_on = .false.
INCLUDE 'xcgrid_theory_levels.h'
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

real(8) :: XCGRID_RSCALE(7)  = (/0.0d0, 0.60d0, 1.848d0, 0.845d0, 2.6d0, 1.30d0, 4.0d0/)
integer :: XCGRID_M3_ANGGRID(7) = (/0, 3, 4, 5, 6, 7, 7/)
real(8) :: XCGRID_M3_RSCALE(7)  = (/0.0d0, 0.65d0, 0.95d0, 1.55d0, 2.10d0, 2.70d0, 3.40d0/)
integer :: XCGRID_M3_NRAD(7)    = (/45, 55, 60, 65, 70, 75, 80/)
integer :: XCGRID_SPH_IN(7)  = (/0, 302, 302, 974, 974, 974, 974/)
integer :: XCGRID_SPH_EDGE(7) = (/0, 302, 302, 590, 590, 590, 974/)

real(8) :: xcgrid_prune_alphas(4,3) = reshape((/ &
   0.25d0,   0.5d0, 1.0d0, 4.5d0,   &
   0.1667d0, 0.5d0, 0.9d0, 3.5d0,   &
   0.1d0,    0.4d0, 0.8d0, 2.5d0 /), (/4,3/))

integer :: xcgrid_prune_orders(5,7) = reshape((/ &
    14,  26,  50,  50,  26,   &
    14,  26,  50, 110,  50,   &
    26,  50, 110, 194, 110,   &
    26, 110, 194, 302, 194,   &
    26, 194, 302, 434, 302,   &
    50, 302, 434, 590, 434,   &
   110, 434, 590, 770, 590 /), (/5,7/))

real(8) :: xcgrid_prune_dmin_bounds(3) = (/ 1.0d0, 2.5d0, 5.0d0 /)

integer :: xcgrid_prune_dmin_orders(4,7) = reshape((/ &
    50,  50,  26,  14,   &
   110,  50,  26,  14,   &
   194, 110,  50,  26,   &
   302, 194, 110,  26,   &
   434, 302, 194,  26,   &
   590, 434, 302,  50,   &
   770, 590, 434, 110 /), (/4,7/))

real(8) :: xcgrid_radpot_scale = 1.0d0

integer :: xcgrid_radpot_abs = 0

real(8) :: xcgrid_radpot_rowscale = 0.0d0

real(8) :: xcgrid_radpot_rowscale_h = 0.0d0

real(8) :: xcgrid_radpot_rowscale_z(103) = 0.0d0

real(8) :: xcgrid_theory_rowscale(7) = 0.0d0

integer :: xcgrid_theory_anggrid(7) = 0

logical :: xcgrid_dump_shells = .false.

logical :: xcgrid_use_zone_table = .false.

logical :: xcgrid_use_dmin_table = .false.

logical :: xcgrid_use_treutler_aij = .false.
logical :: cosx_becke_treutler = .false.

logical :: xcgrid_use_theory_table = .false.
integer, parameter :: XCGRID_LVL_SMOOTHDFT = 1, XCGRID_LVL_METAGGA = 2, &
                       XCGRID_LVL_R2SCAN   = 3, XCGRID_LVL_M06     = 4, &
                       XCGRID_LVL_WB97XD   = 5, XCGRID_LVL_WB97MV  = 6, &
                       XCGRID_LVL_R2SCAN3C = 7
integer, parameter :: XCGRID_TIER_LO = 1, XCGRID_TIER_REGULAR = 2, XCGRID_TIER_HIGH = 3
integer :: xcgrid_theory_level = XCGRID_LVL_SMOOTHDFT
integer :: xcgrid_theory_tier = XCGRID_TIER_REGULAR

integer :: xcgrid_zone1 = 26
integer :: xcgrid_zone2 = 50
integer :: xcgrid_zone3 = 110
integer :: xcgrid_zone5 = 146
integer :: xcgrid_zone4 = -1
integer :: xcgrid_zone4_h = -1
integer :: xcgrid_zone4_z(103) = -1
integer :: xcgrid_zone1_z(103) = -1
integer :: xcgrid_zone2_z(103) = -1
integer :: xcgrid_zone3_z(103) = -1
real(8) :: xcgrid_zone_alpha3 = -1.0d0

logical :: xcgrid_rho_dump = .false.

integer :: xcgrid_ang_margin = 0

logical :: xcgrid_audit_only = .false.
real(8) :: xcgrid_int_rho = 0.0d0
integer :: xcgrid_level = XCGRID_FINE
integer :: xcgrid_active_coarse = XCGRID_COARSE
integer :: xcgrid_active_fine   = XCGRID_FINE
logical :: xcgrid_dynamic = .false.
logical :: xcgrid_refined = .false.
logical :: xcgrid_theory_dynamic = .false.
logical :: xcgrid_theory_refined = .false.
integer :: xcgrid_theory_dynamic_level = XCGRID_LVL_SMOOTHDFT
integer :: xcgrid_theory_dynamic_tier  = XCGRID_TIER_REGULAR
logical :: force_dense = .false.
logical :: force_dense_mgga = .false.
logical :: force_dense_allow = .true.
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
real(8), parameter :: TA_XI(103) = (/ &
   0.80d0, 0.90d0, 1.80d0, 1.40d0, 1.30d0, 1.10d0, 0.90d0, 0.90d0, &
   0.90d0, 0.90d0, 1.40d0, 1.30d0, 1.30d0, 1.20d0, 1.10d0, 1.00d0, &
   1.00d0, 1.00d0, 1.50d0, 1.40d0, 1.30d0, 1.20d0, 1.20d0, 1.20d0, &
   1.20d0, 1.20d0, 1.20d0, 1.10d0, 1.10d0, 1.10d0, 1.10d0, 1.00d0, &
   0.90d0, 0.90d0, 0.90d0, 0.90d0, 2.00d0, 1.70d0, 1.50d0, 1.50d0, &
   1.35d0, 1.35d0, 1.25d0, 1.20d0, 1.25d0, 1.30d0, 1.50d0, 1.50d0, &
   1.30d0, 1.20d0, 1.20d0, 1.15d0, 1.15d0, 1.15d0, 2.50d0, 2.20d0, &
   2.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, &
   1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, &
   1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, &
   1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 2.50d0, 2.10d0, &
   3.69d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, &
   1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0, 1.50d0 /)
z = max(1, min(103, Zin))
xi = TA_XI(z)
end function xcgrid_m3_xi

subroutine xcgrid_m3_radial(Zin, i, nrad, radr, radw)
implicit none
integer, intent(in) :: Zin, i, nrad
real(8), intent(out) :: radr, radw
real(8) :: pig, ln2, theta, x, s, xi, drdx, pa, plog
real(8), parameter :: ALPHA = 0.6d0
pig  = acos(-1.0d0)
ln2  = log(2.0d0)
xi   = xcgrid_m3_xi(Zin)
theta = i * pig / (nrad + 1)
x = cos(theta)
s = sin(theta)
pa   = (1.0d0 + x)**ALPHA
plog = log(2.0d0/(1.0d0 - x))
radr = (xi/ln2) * pa * plog
drdx = (xi/ln2) * pa * (ALPHA/(1.0d0 + x)*plog + 1.0d0/(1.0d0 - x))
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

integer function xcgrid_prune_sphpot_zone(n_ang, ratio, row, z) result(cursphpot)
implicit none
integer,intent(in) :: n_ang, row, z
real(8),intent(in) :: ratio
real(8) :: alphas(4)
integer :: place
integer :: pat302(0:4), pat434(0:4)
pat302 = (/50,86,266,302,266/)
pat434 = (/50,86,350,434,350/)
if (xcgrid_zone1 > 0) pat434(0) = xcgrid_zone1
if (xcgrid_zone2 > 0) pat434(1) = xcgrid_zone2
if (xcgrid_zone3 > 0) pat434(2) = xcgrid_zone3
if (xcgrid_zone5 > 0) pat434(4) = xcgrid_zone5
if (xcgrid_zone4 > 0) pat434(3) = xcgrid_zone4
if (xcgrid_zone1 > 0) pat302(0) = xcgrid_zone1
if (xcgrid_zone2 > 0) pat302(1) = xcgrid_zone2
if (xcgrid_zone3 > 0) pat302(2) = xcgrid_zone3
if (xcgrid_zone5 > 0) pat302(4) = xcgrid_zone5
if (xcgrid_zone4_h > 0) pat302(3) = xcgrid_zone4_h
select case (row)
case (1); alphas = (/0.25d0, 0.5d0, 1.0d0, 4.5d0/)
case (2); alphas = (/0.1667d0, 0.5d0, 0.9d0, 3.5d0/)
case default; alphas = (/0.1d0, 0.4d0, 0.8d0, 2.5d0/)
end select
if (row .eq. 1 .and. xcgrid_zone_alpha3 > 0.0d0) alphas(4) = xcgrid_zone_alpha3
place = count(ratio .gt. alphas)
if (n_ang .eq. 302) then
   cursphpot = pat302(place)
else
   cursphpot = pat434(place)
endif
if (z .ge. 1 .and. z .le. 103) then
   if (place .eq. 0 .and. xcgrid_zone1_z(z) > 0) cursphpot = xcgrid_zone1_z(z)
   if (place .eq. 1 .and. xcgrid_zone2_z(z) > 0) cursphpot = xcgrid_zone2_z(z)
   if (place .eq. 2 .and. xcgrid_zone3_z(z) > 0) cursphpot = xcgrid_zone3_z(z)
   if (place .eq. 3 .and. xcgrid_zone4_z(z) > 0) cursphpot = xcgrid_zone4_z(z)
endif
end function xcgrid_prune_sphpot_zone

logical function xcgrid_zone_calibrated_for(z) result(ok)
integer, intent(in) :: z
ok = .false.
if (z .ge. 1 .and. z .le. 10) then
   ok = .true.
else if (z .ge. 1 .and. z .le. 103) then
   ok = (xcgrid_zone1_z(z) > 0) .or. (xcgrid_zone2_z(z) > 0) .or. &
        (xcgrid_zone3_z(z) > 0) .or. (xcgrid_radpot_rowscale_z(z) > 0.0d0)
endif
end function xcgrid_zone_calibrated_for

integer function xcgrid_sphpot_ladder(sphpot, anggrid_legacy, ratio, row, acap, &
                                       use_dmin, d_min) result(cursphpot)
implicit none
integer, intent(in) :: sphpot, anggrid_legacy, row, acap
real(8), intent(in) :: ratio, d_min
logical, intent(in) :: use_dmin
if (anggrid_legacy > 0) then
   cursphpot = xcgrid_prune_angular(anggrid_legacy, ratio, row, acap)
else if (use_dmin) then
   cursphpot = xcgrid_prune_angular_dmin(sphpot, d_min, acap)
else
   cursphpot = xcgrid_prune_angular(sphpot, ratio, row, acap)
endif
end function xcgrid_sphpot_ladder

integer function xcgrid_sphpot_theory(ratio, row, z, acap) result(cursphpot)
implicit none
integer, intent(in) :: row, z, acap
real(8), intent(in) :: ratio
cursphpot = min(xcgrid_prune_sphpot_zone(merge(302, 434, row .eq. 1), ratio, row, z), &
                 xcgrid_lebedev_floor_to_supported(acap))
end function xcgrid_sphpot_theory

integer function xcgrid_atom_sphpot(sphpot, anggrid_legacy, ratio, row, z, acap, &
                                     use_dmin, d_min) result(cursphpot)
implicit none
integer, intent(in) :: sphpot, anggrid_legacy, row, z, acap
real(8), intent(in) :: ratio, d_min
logical, intent(in) :: use_dmin
if (anggrid_legacy .eq. 0 .and. xcgrid_use_zone_table .and. xcgrid_zone_calibrated_for(z)) then
   cursphpot = xcgrid_sphpot_theory(ratio, row, z, acap)
else
   cursphpot = xcgrid_sphpot_ladder(sphpot, anggrid_legacy, ratio, row, acap, use_dmin, d_min)
endif
end function xcgrid_atom_sphpot

integer function xcgrid_radpot_ladder(row) result(radpot)
implicit none
integer, intent(in) :: row
radpot = max(10, nint(XCGRID_M3_NRAD(row) * XCGRID_M3_RSCALE(xcgrid_level) * xcgrid_radpot_scale))
end function xcgrid_radpot_ladder

integer function xcgrid_radpot_theory(row) result(radpot)
implicit none
integer, intent(in) :: row
radpot = max(10, nint(XCGRID_M3_NRAD(row) * xcgrid_theory_rowscale(row)))
end function xcgrid_radpot_theory

integer function xcgrid_prune_sphpot_cell(n_ang, ratio, row, z, zone434, zone302) result(cursphpot)
implicit none
integer,intent(in) :: n_ang, row, z, zone434(5), zone302(5)
real(8),intent(in) :: ratio
real(8) :: alphas(4)
integer :: place, pat302(0:4), pat434(0:4)
pat434 = (/zone434(1), zone434(2), zone434(3), zone434(4), zone434(5)/)
pat302 = (/zone302(1), zone302(2), zone302(3), zone302(4), zone302(5)/)
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
if (z .ge. 1 .and. z .le. 103) then
   if (place .eq. 0 .and. xcgrid_zone1_z(z) > 0) cursphpot = xcgrid_zone1_z(z)
   if (place .eq. 1 .and. xcgrid_zone2_z(z) > 0) cursphpot = xcgrid_zone2_z(z)
   if (place .eq. 2 .and. xcgrid_zone3_z(z) > 0) cursphpot = xcgrid_zone3_z(z)
   if (place .eq. 3 .and. xcgrid_zone4_z(z) > 0) cursphpot = xcgrid_zone4_z(z)
endif
end function xcgrid_prune_sphpot_cell

integer function xcgrid_lvl_from_str(s) result(lvl)
implicit none
character(len=*), intent(in) :: s
character(len=32) :: u
integer :: i
u = adjustl(s)
do i = 1, len_trim(u)
   if (u(i:i) >= 'a' .and. u(i:i) <= 'z') u(i:i) = achar(iachar(u(i:i)) - 32)
enddo
select case (trim(u))
case ("METAGGA", "MGGA");            lvl = XCGRID_LVL_METAGGA
case ("R2SCAN");                     lvl = XCGRID_LVL_R2SCAN
case ("R2SCAN3C", "R2SCAN-3C", "R2SCAN_3C");  lvl = XCGRID_LVL_R2SCAN3C
case ("M06");                        lvl = XCGRID_LVL_M06
case ("WB97XD", "WB97X-D");          lvl = XCGRID_LVL_WB97XD
case ("WB97MV", "WB97M-V");          lvl = XCGRID_LVL_WB97MV
case default;                        lvl = XCGRID_LVL_SMOOTHDFT
end select
end function xcgrid_lvl_from_str

integer function xcgrid_theory_level_for_functional(functional) result(lvl)
implicit none
character(len=*), intent(in) :: functional
character(len=32) :: u
integer :: i
u = adjustl(functional)
do i = 1, len_trim(u)
   if (u(i:i) >= 'a' .and. u(i:i) <= 'z') u(i:i) = achar(iachar(u(i:i)) - 32)
enddo
select case (trim(u))
case ("R2SCAN-3C", "R2SCAN_3C", "R2SCAN3C");  lvl = XCGRID_LVL_R2SCAN3C
case ("R2SCAN", "R2SCAN0");  lvl = XCGRID_LVL_R2SCAN
case ("M06L", "M06-L", "M06_L", "M06", "M06_HYB", "M06-2X", "M06-2X_HYB", &
      "M062X", "M06HYB", "M062XHYB", &
      "M05-2X", "M052X", "MN15", "MN15L", "MN15-L");
   lvl = XCGRID_LVL_M06
case ("WB97XD", "WB97X-D", "WB97X");          lvl = XCGRID_LVL_WB97XD
case ("WB97MV", "WB97M-V", "B97MV", "B97M-V");  lvl = XCGRID_LVL_WB97MV
case default;                        lvl = XCGRID_LVL_SMOOTHDFT
end select
end function xcgrid_theory_level_for_functional

logical function xcgrid_theory_default_on(functional) result(ok)
implicit none
character(len=*), intent(in) :: functional
character(len=32) :: u
integer :: i
u = adjustl(functional)
do i = 1, len_trim(u)
   if (u(i:i) >= 'a' .and. u(i:i) <= 'z') u(i:i) = achar(iachar(u(i:i)) - 32)
enddo
select case (trim(u))
case default
   ok = .true.
end select
end function xcgrid_theory_default_on

integer function xcgrid_tier_from_str(s) result(tier)
implicit none
character(len=*), intent(in) :: s
character(len=32) :: u
integer :: i
u = adjustl(s)
do i = 1, len_trim(u)
   if (u(i:i) >= 'a' .and. u(i:i) <= 'z') u(i:i) = achar(iachar(u(i:i)) - 32)
enddo
select case (trim(u))
case ("LO", "LOW");                  tier = XCGRID_TIER_LO
case ("HIGH", "HI");                 tier = XCGRID_TIER_HIGH
case default;                        tier = XCGRID_TIER_REGULAR
end select
end function xcgrid_tier_from_str

function xcgrid_lvl_name(lvl) result(nm)
implicit none
integer, intent(in) :: lvl
character(len=9) :: nm
select case (lvl)
case (XCGRID_LVL_METAGGA); nm = "METAGGA  "
case (XCGRID_LVL_R2SCAN);  nm = "R2SCAN   "
case (XCGRID_LVL_R2SCAN3C); nm = "R2SCAN3C "
case (XCGRID_LVL_M06);     nm = "M06      "
case (XCGRID_LVL_WB97XD);  nm = "WB97XD   "
case (XCGRID_LVL_WB97MV);  nm = "WB97MV   "
case default;              nm = "SMOOTHDFT"
end select
end function xcgrid_lvl_name

function xcgrid_tier_name(tier) result(nm)
implicit none
integer, intent(in) :: tier
character(len=7) :: nm
select case (tier)
case (XCGRID_TIER_LO);   nm = "LO     "
case (XCGRID_TIER_HIGH); nm = "HIGH   "
case default;            nm = "REGULAR"
end select
end function xcgrid_tier_name

subroutine xcgrid_theory_cell(level, tier, rs, zone434, zone302, aij, ovn, ovz, ovs)
implicit none
integer, intent(in)  :: level, tier
real(8), intent(out) :: rs(7)
integer, intent(out) :: zone434(5), zone302(5), aij
integer, intent(out) :: ovn, ovz(16)
real(8), intent(out) :: ovs(16)
ovn = 0
select case (level)
case (XCGRID_LVL_METAGGA)
   select case (tier)
   case (XCGRID_TIER_LO)
      rs = XCGRID_ROWSCALE_METAGGA_LO
      zone434 = XCGRID_ZONE434_METAGGA_LO; zone302 = XCGRID_ZONE302_METAGGA_LO
      aij = XCGRID_AIJ_SCHEME_METAGGA_LO
      ovn = XCGRID_ELEM_OVERRIDE_METAGGA_LO_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_METAGGA_LO_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_METAGGA_LO_SCALE(1:ovn)
      endif
   case (XCGRID_TIER_HIGH)
      rs = XCGRID_ROWSCALE_METAGGA_HIGH
      zone434 = XCGRID_ZONE434_METAGGA_HIGH; zone302 = XCGRID_ZONE302_METAGGA_HIGH
      aij = XCGRID_AIJ_SCHEME_METAGGA_HIGH
      ovn = XCGRID_ELEM_OVERRIDE_METAGGA_HIGH_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_METAGGA_HIGH_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_METAGGA_HIGH_SCALE(1:ovn)
      endif
   case default
      rs = XCGRID_ROWSCALE_METAGGA_REGULAR
      zone434 = XCGRID_ZONE434_METAGGA_REGULAR; zone302 = XCGRID_ZONE302_METAGGA_REGULAR
      aij = XCGRID_AIJ_SCHEME_METAGGA_REGULAR
      ovn = XCGRID_ELEM_OVERRIDE_METAGGA_REGULAR_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_METAGGA_REGULAR_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_METAGGA_REGULAR_SCALE(1:ovn)
      endif
   end select
case (XCGRID_LVL_R2SCAN)
   select case (tier)
   case (XCGRID_TIER_LO)
      rs = XCGRID_ROWSCALE_R2SCAN_LO
      zone434 = XCGRID_ZONE434_R2SCAN_LO; zone302 = XCGRID_ZONE302_R2SCAN_LO
      aij = XCGRID_AIJ_SCHEME_R2SCAN_LO
      ovn = XCGRID_ELEM_OVERRIDE_R2SCAN_LO_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN_LO_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN_LO_SCALE(1:ovn)
      endif
   case (XCGRID_TIER_HIGH)
      rs = XCGRID_ROWSCALE_R2SCAN_HIGH
      zone434 = XCGRID_ZONE434_R2SCAN_HIGH; zone302 = XCGRID_ZONE302_R2SCAN_HIGH
      aij = XCGRID_AIJ_SCHEME_R2SCAN_HIGH
      ovn = XCGRID_ELEM_OVERRIDE_R2SCAN_HIGH_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN_HIGH_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN_HIGH_SCALE(1:ovn)
      endif
   case default
      rs = XCGRID_ROWSCALE_R2SCAN_REGULAR
      zone434 = XCGRID_ZONE434_R2SCAN_REGULAR; zone302 = XCGRID_ZONE302_R2SCAN_REGULAR
      aij = XCGRID_AIJ_SCHEME_R2SCAN_REGULAR
      ovn = XCGRID_ELEM_OVERRIDE_R2SCAN_REGULAR_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN_REGULAR_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN_REGULAR_SCALE(1:ovn)
      endif
   end select
case (XCGRID_LVL_R2SCAN3C)
   if (tier .eq. XCGRID_TIER_REGULAR) then
      rs = XCGRID_ROWSCALE_R2SCAN3C_REGULAR
      zone434 = XCGRID_ZONE434_R2SCAN3C_REGULAR; zone302 = XCGRID_ZONE302_R2SCAN3C_REGULAR
      aij = XCGRID_AIJ_SCHEME_R2SCAN_REGULAR
      ovn = XCGRID_ELEM_OVERRIDE_R2SCAN3C_REGULAR_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN3C_REGULAR_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_R2SCAN3C_REGULAR_SCALE(1:ovn)
      endif
   else
      call xcgrid_theory_cell(XCGRID_LVL_R2SCAN, tier, rs, zone434, zone302, aij, ovn, ovz, ovs)
   endif
case (XCGRID_LVL_M06)
   select case (tier)
   case (XCGRID_TIER_LO)
      rs = XCGRID_ROWSCALE_M06_LO
      zone434 = XCGRID_ZONE434_M06_LO; zone302 = XCGRID_ZONE302_M06_LO
      aij = XCGRID_AIJ_SCHEME_M06_LO
      ovn = XCGRID_ELEM_OVERRIDE_M06_LO_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_M06_LO_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_M06_LO_SCALE(1:ovn)
      endif
   case (XCGRID_TIER_HIGH)
      rs = XCGRID_ROWSCALE_M06_HIGH
      zone434 = XCGRID_ZONE434_M06_HIGH; zone302 = XCGRID_ZONE302_M06_HIGH
      aij = XCGRID_AIJ_SCHEME_M06_HIGH
      ovn = XCGRID_ELEM_OVERRIDE_M06_HIGH_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_M06_HIGH_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_M06_HIGH_SCALE(1:ovn)
      endif
   case default
      rs = XCGRID_ROWSCALE_M06_REGULAR
      zone434 = XCGRID_ZONE434_M06_REGULAR; zone302 = XCGRID_ZONE302_M06_REGULAR
      aij = XCGRID_AIJ_SCHEME_M06_REGULAR
      ovn = XCGRID_ELEM_OVERRIDE_M06_REGULAR_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_M06_REGULAR_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_M06_REGULAR_SCALE(1:ovn)
      endif
   end select
case (XCGRID_LVL_WB97XD)
   select case (tier)
   case (XCGRID_TIER_LO)
      rs = XCGRID_ROWSCALE_WB97XD_LO
      zone434 = XCGRID_ZONE434_WB97XD_LO; zone302 = XCGRID_ZONE302_WB97XD_LO
      aij = XCGRID_AIJ_SCHEME_WB97XD_LO
      ovn = XCGRID_ELEM_OVERRIDE_WB97XD_LO_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97XD_LO_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97XD_LO_SCALE(1:ovn)
      endif
   case (XCGRID_TIER_HIGH)
      rs = XCGRID_ROWSCALE_WB97XD_HIGH
      zone434 = XCGRID_ZONE434_WB97XD_HIGH; zone302 = XCGRID_ZONE302_WB97XD_HIGH
      aij = XCGRID_AIJ_SCHEME_WB97XD_HIGH
      ovn = XCGRID_ELEM_OVERRIDE_WB97XD_HIGH_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97XD_HIGH_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97XD_HIGH_SCALE(1:ovn)
      endif
   case default
      rs = XCGRID_ROWSCALE_WB97XD_REGULAR
      zone434 = XCGRID_ZONE434_WB97XD_REGULAR; zone302 = XCGRID_ZONE302_WB97XD_REGULAR
      aij = XCGRID_AIJ_SCHEME_WB97XD_REGULAR
      ovn = XCGRID_ELEM_OVERRIDE_WB97XD_REGULAR_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97XD_REGULAR_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97XD_REGULAR_SCALE(1:ovn)
      endif
   end select
case (XCGRID_LVL_WB97MV)
   select case (tier)
   case (XCGRID_TIER_LO)
      rs = XCGRID_ROWSCALE_WB97MV_LO
      zone434 = XCGRID_ZONE434_WB97MV_LO; zone302 = XCGRID_ZONE302_WB97MV_LO
      aij = XCGRID_AIJ_SCHEME_WB97MV_LO
      ovn = XCGRID_ELEM_OVERRIDE_WB97MV_LO_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97MV_LO_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97MV_LO_SCALE(1:ovn)
      endif
   case (XCGRID_TIER_HIGH)
      rs = XCGRID_ROWSCALE_WB97MV_HIGH
      zone434 = XCGRID_ZONE434_WB97MV_HIGH; zone302 = XCGRID_ZONE302_WB97MV_HIGH
      aij = XCGRID_AIJ_SCHEME_WB97MV_HIGH
      ovn = XCGRID_ELEM_OVERRIDE_WB97MV_HIGH_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97MV_HIGH_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97MV_HIGH_SCALE(1:ovn)
      endif
   case default
      rs = XCGRID_ROWSCALE_WB97MV_REGULAR
      zone434 = XCGRID_ZONE434_WB97MV_REGULAR; zone302 = XCGRID_ZONE302_WB97MV_REGULAR
      aij = XCGRID_AIJ_SCHEME_WB97MV_REGULAR
      ovn = XCGRID_ELEM_OVERRIDE_WB97MV_REGULAR_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97MV_REGULAR_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_WB97MV_REGULAR_SCALE(1:ovn)
      endif
   end select
case default
   select case (tier)
   case (XCGRID_TIER_LO)
      rs = XCGRID_ROWSCALE_SMOOTHDFT_LO
      zone434 = XCGRID_ZONE434_SMOOTHDFT_LO; zone302 = XCGRID_ZONE302_SMOOTHDFT_LO
      aij = XCGRID_AIJ_SCHEME_SMOOTHDFT_LO
      ovn = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_LO_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_LO_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_LO_SCALE(1:ovn)
      endif
   case (XCGRID_TIER_HIGH)
      rs = XCGRID_ROWSCALE_SMOOTHDFT_HIGH
      zone434 = XCGRID_ZONE434_SMOOTHDFT_HIGH; zone302 = XCGRID_ZONE302_SMOOTHDFT_HIGH
      aij = XCGRID_AIJ_SCHEME_SMOOTHDFT_HIGH
      ovn = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_HIGH_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_HIGH_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_HIGH_SCALE(1:ovn)
      endif
   case default
      rs = XCGRID_ROWSCALE_SMOOTHDFT_REGULAR
      zone434 = XCGRID_ZONE434_SMOOTHDFT_REGULAR; zone302 = XCGRID_ZONE302_SMOOTHDFT_REGULAR
      aij = XCGRID_AIJ_SCHEME_SMOOTHDFT_REGULAR
      ovn = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_REGULAR_N
      if (ovn > 0) then
         ovz(1:ovn) = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_REGULAR_Z(1:ovn)
         ovs(1:ovn) = XCGRID_ELEM_OVERRIDE_SMOOTHDFT_REGULAR_SCALE(1:ovn)
      endif
   end select
end select
end subroutine xcgrid_theory_cell

integer function xcgrid_theory_legacy_anggrid(level, tier) result(ag)
implicit none
integer, intent(in) :: level, tier
if (tier .eq. XCGRID_TIER_LO) then
   ag = XCGRID_M3_ANGGRID(XCGRID_COARSE)
else if (level .eq. XCGRID_LVL_SMOOTHDFT) then
   ag = XCGRID_M3_ANGGRID(XCGRID_FINE)
else
   ag = XCGRID_M3_ANGGRID(XCGRID_L5)
endif
end function xcgrid_theory_legacy_anggrid

subroutine xcgrid_theory_apply()
implicit none
real(8) :: rs(7)
integer :: zone434(5), zone302(5), aij, ovn, ovz(16), k
real(8) :: ovs(16)
integer :: eff_level, eff_tier, ag

eff_level = xcgrid_theory_level
eff_tier  = xcgrid_theory_tier
call xcgrid_theory_cell(xcgrid_theory_level, xcgrid_theory_tier, &
                        rs, zone434, zone302, aij, ovn, ovz, ovs)

if (rs(1) < 0.0d0 .or. zone434(1) < 0 .or. zone302(1) < 0 .or. aij < 0) then
   eff_level = XCGRID_LVL_SMOOTHDFT
   eff_tier  = XCGRID_TIER_REGULAR
   call xcgrid_theory_cell(eff_level, eff_tier, &
                           rs, zone434, zone302, aij, ovn, ovz, ovs)
endif

block
   character(len=256) :: tn
   integer :: tios
   call get_environment_variable("ENGINE_XC_TUNE_RS", tn)
   if (len_trim(tn) > 0) read(tn, *, iostat=tios) rs
   tn = ''
   call get_environment_variable("ENGINE_XC_TUNE_ELEM", tn)
   if (len_trim(tn) > 0 .and. ovn > 0) read(tn, *, iostat=tios) ovs(1:ovn)
   tn = ''
   call get_environment_variable("ENGINE_XC_TUNE_Z434", tn)
   if (len_trim(tn) > 0) read(tn, *, iostat=tios) zone434
   tn = ''
   call get_environment_variable("ENGINE_XC_TUNE_Z302", tn)
   if (len_trim(tn) > 0) read(tn, *, iostat=tios) zone302
end block

xcgrid_theory_anggrid = 0

xcgrid_theory_rowscale   = rs
xcgrid_use_zone_table   = .true.
xcgrid_use_treutler_aij  = (aij .eq. 2)
xcgrid_zone1 = zone434(1)
xcgrid_zone2 = zone434(2)
xcgrid_zone3 = zone434(3)
xcgrid_zone4 = zone434(4)
xcgrid_zone5 = zone434(5)
xcgrid_zone4_h = zone302(4)

xcgrid_radpot_rowscale_z = 0.0d0
do k = 1, ovn
   xcgrid_radpot_rowscale_z(ovz(k)) = ovs(k)
enddo

xcgrid_zone4_z = -1
xcgrid_zone4_z(11:18)  = 434
xcgrid_zone4_z(19:36)  = 590
xcgrid_zone4_z(37:54)  = 770
xcgrid_zone4_z(55:103) = 974
do k = 1, ovn
   xcgrid_zone4_z(ovz(k)) = -1
enddo

if (eff_level .eq. XCGRID_LVL_SMOOTHDFT .and. eff_tier .eq. XCGRID_TIER_HIGH) then
   xcgrid_zone4_z(26) = 770
endif

xcgrid_zone1_z = -1
xcgrid_zone2_z = -1
xcgrid_zone3_z = -1
if ((eff_level .eq. XCGRID_LVL_R2SCAN .or. eff_level .eq. XCGRID_LVL_R2SCAN3C) &
    .and. eff_tier .eq. XCGRID_TIER_LO) then
   xcgrid_zone1_z(26) = 86
   xcgrid_zone2_z(26) = 266
   xcgrid_zone3_z(26) = 230
   xcgrid_radpot_rowscale_z(26) = 1.4d0
else if (eff_level .eq. XCGRID_LVL_R2SCAN .or. eff_level .eq. XCGRID_LVL_R2SCAN3C) then
   xcgrid_zone1_z(26) = 110
   xcgrid_zone2_z(26) = 302
   xcgrid_zone3_z(26) = 266
   xcgrid_radpot_rowscale_z(26) = 3.0d0
endif
if (eff_level .eq. XCGRID_LVL_WB97XD) then
   xcgrid_zone1_z(26) = 110
   xcgrid_zone2_z(26) = 302
   xcgrid_zone3_z(26) = 266
   xcgrid_radpot_rowscale_z(26) = 3.0d0
endif
if (eff_level .eq. XCGRID_LVL_WB97MV) then
   xcgrid_zone1_z(26) = 110
   xcgrid_zone2_z(26) = 302
   xcgrid_zone3_z(26) = 266
   xcgrid_radpot_rowscale_z(26) = 3.0d0
endif
if (eff_level .eq. XCGRID_LVL_SMOOTHDFT) then
   xcgrid_zone1_z(26) = 86
   xcgrid_zone2_z(26) = 266
   xcgrid_zone3_z(26) = 230
endif

if (eff_level .eq. XCGRID_LVL_WB97XD) then
   xcgrid_zone1_z(44) = 110
   xcgrid_zone2_z(44) = 302
   xcgrid_zone3_z(44) = 266
   xcgrid_radpot_rowscale_z(44) = 3.0d0
endif

if (eff_level .eq. XCGRID_LVL_WB97MV) then
   xcgrid_radpot_rowscale_z(44) = 1.4d0
endif

if (eff_level .eq. XCGRID_LVL_SMOOTHDFT) then
   xcgrid_zone2_z(19:36) = 86
   xcgrid_zone2_z(26) = 266
endif

if (eff_level .eq. XCGRID_LVL_SMOOTHDFT) then
   xcgrid_radpot_rowscale_z(44) = 1.4d0
endif

if (eff_level .eq. XCGRID_LVL_SMOOTHDFT .and. eff_tier .eq. XCGRID_TIER_HIGH) then
   xcgrid_radpot_rowscale_z(26) = 1.4d0
endif

if (eff_level .eq. XCGRID_LVL_SMOOTHDFT) then
   xcgrid_radpot_rowscale_z(11:18) = 1.5d0
endif

end subroutine xcgrid_theory_apply

integer function xcgrid_prune_angular(anggrid, ratio, row, ang_cap) result(cursphpot)
implicit none
integer,intent(in) :: anggrid, row
real(8),intent(in) :: ratio
integer,intent(in) :: ang_cap
integer :: place, ag
ag = max(1, min(7, anggrid))
place = count(ratio .gt. xcgrid_prune_alphas(:, min(3, max(1, row))))
cursphpot = xcgrid_prune_orders(place+1, ag)
cursphpot = min(cursphpot, xcgrid_lebedev_floor_to_supported(ang_cap))
end function xcgrid_prune_angular

subroutine xcgrid_dmin_shell(radr, Rij_row, iatm, natoms, d_min)
implicit none
integer, intent(in) :: iatm, natoms
real(8), intent(in) :: radr, Rij_row(natoms)
real(8), intent(out) :: d_min
integer :: jj
d_min = huge(1.0d0)
do jj = 1, natoms
   if (jj .eq. iatm) cycle
   d_min = min(d_min, abs(radr - 0.5d0*Rij_row(jj)))
enddo
end subroutine xcgrid_dmin_shell

integer function xcgrid_prune_angular_dmin(anggrid, d_min, ang_cap) result(cursphpot)
implicit none
integer,intent(in) :: anggrid
real(8),intent(in) :: d_min
integer,intent(in) :: ang_cap
integer :: place, ag
ag = max(1, min(7, anggrid))
place = count(d_min .gt. xcgrid_prune_dmin_bounds)
cursphpot = xcgrid_prune_dmin_orders(place+1, ag)
cursphpot = min(cursphpot, xcgrid_lebedev_floor_to_supported(ang_cap))
end function xcgrid_prune_angular_dmin

integer function xcgrid_lebedev_for_degree(degree) result(order)
implicit none
integer, intent(in) :: degree
integer, parameter :: ORD(19) = (/ 6, 14, 26, 38, 50, 74, 86, 110, 146, &
     170, 194, 230, 266, 302, 350, 434, 590, 770, 974 /)
integer, parameter :: LDEG(19) = (/ 3,  5,  7,  9, 11, 13, 15,  17,  19, &
       21,  23,  25,  27,  29,  31,  35,  41,  47,  53 /)
integer :: k
order = ORD(19)
do k = 1, 19
   if (LDEG(k) >= degree) then
      order = ORD(k)
      return
   endif
enddo
end function xcgrid_lebedev_for_degree

integer function xcgrid_lebedev_floor_to_supported(n) result(order)
implicit none
integer, intent(in) :: n
integer, parameter :: ORD(19) = (/ 6, 14, 26, 38, 50, 74, 86, 110, 146, &
     170, 194, 230, 266, 302, 350, 434, 590, 770, 974 /)
integer :: k
order = ORD(1)
do k = 19, 1, -1
   if (ORD(k) <= n) then
      order = ORD(k)
      return
   endif
enddo
end function xcgrid_lebedev_floor_to_supported

integer function xcgrid_ang_cap(np_shell, nd_shell, nf_shell, ng_shell) result(cap)
implicit none
integer, intent(in) :: np_shell, nd_shell, nf_shell, ng_shell
integer :: lmax
if (xcgrid_ang_margin <= 0) then
   cap = 100000
   return
endif
lmax = 0
if (np_shell > 0) lmax = 1
if (nd_shell > 0) lmax = 2
if (nf_shell > 0) lmax = 3
if (ng_shell > 0) lmax = 4
cap = xcgrid_lebedev_for_degree(2*lmax + xcgrid_ang_margin)
end function xcgrid_ang_cap

subroutine xcgrid_tune_from_env(functional)
implicit none
character(len=*), intent(in) :: functional
character(len=512) :: buf
integer :: ios, i, j, ord(5), dmord(4)
real(8) :: a(12), sc, dmb(3)
logical :: level_auto_detected

call get_environment_variable("ENGINE_XC_GRID", buf)
if (len_trim(buf) > 0) then
   xcgrid_use_theory_table = .true.
   block
      character(len=512) :: seg1, seg2
      integer :: slashpos, dashpos
      buf = adjustl(buf)
      slashpos = index(trim(buf), "/")
      if (slashpos > 0) then
         seg1 = buf(1:slashpos-1); seg2 = buf(slashpos+1:)
      else
         seg1 = buf; seg2 = ""
      endif
      dashpos = index(trim(seg1), "-")
      if (dashpos > 0) then
         xcgrid_theory_level = xcgrid_lvl_from_str(seg1(1:dashpos-1))
         xcgrid_theory_tier  = xcgrid_tier_from_str(seg1(dashpos+1:))
      else
         xcgrid_theory_level = xcgrid_lvl_from_str(seg1)
         xcgrid_theory_tier  = XCGRID_TIER_REGULAR
      endif
      if (len_trim(seg2) > 0) then
         dashpos = index(trim(seg2), "-")
         if (dashpos > 0) then
            xcgrid_theory_dynamic_level = xcgrid_lvl_from_str(seg2(1:dashpos-1))
            xcgrid_theory_dynamic_tier  = xcgrid_tier_from_str(seg2(dashpos+1:))
         else
            xcgrid_theory_dynamic_level = xcgrid_lvl_from_str(seg2)
            xcgrid_theory_dynamic_tier  = XCGRID_TIER_REGULAR
         endif
         xcgrid_theory_dynamic = .true.
      endif
   end block
   if (xcgrid_theory_dynamic) then
      print '(A,ES9.2,A)', " ENGINE_XC_GRID active: start="// &
         trim(xcgrid_lvl_name(xcgrid_theory_level))//"-"// &
         trim(xcgrid_tier_name(xcgrid_theory_tier))//"  target="// &
         trim(xcgrid_lvl_name(xcgrid_theory_dynamic_level))//"-"// &
         trim(xcgrid_tier_name(xcgrid_theory_dynamic_tier))// &
         " (promote at dP <", xcgrid_switch_prms, ")"
   else
      print '(A)', " ENGINE_XC_GRID active: "//trim(xcgrid_lvl_name(xcgrid_theory_level))// &
         "-"//trim(xcgrid_tier_name(xcgrid_theory_tier))
   endif
else
   call get_environment_variable("ENGINE_XC_USE_THEORY_TABLE", buf)
   if (len_trim(buf) > 0) then
      xcgrid_use_theory_table = (trim(adjustl(buf)) == "1")
   else
      xcgrid_use_theory_table = xcgrid_theory_default_on(functional)
   endif
   if (xcgrid_use_theory_table) then
      call get_environment_variable("ENGINE_XC_THEORY_LEVEL", buf)
      level_auto_detected = (len_trim(buf) == 0)
      if (level_auto_detected) then
         xcgrid_theory_level = xcgrid_theory_level_for_functional(functional)
      else
         xcgrid_theory_level = xcgrid_lvl_from_str(buf)
      endif
      call get_environment_variable("ENGINE_XC_THEORY_TIER", buf)
      xcgrid_theory_tier = xcgrid_tier_from_str(buf)
      print '(A)', " XC grid: "//trim(xcgrid_lvl_name(xcgrid_theory_level))//" table, "// &
         trim(xcgrid_tier_name(xcgrid_theory_tier))//" tier"
      call get_environment_variable("ENGINE_XC_THEORY_DYNAMIC", buf)
      xcgrid_theory_dynamic = (len_trim(buf) > 0 .and. trim(adjustl(buf)) == "1")
      if (xcgrid_theory_dynamic) then
         xcgrid_theory_dynamic_level = xcgrid_theory_level
         xcgrid_theory_dynamic_tier  = XCGRID_TIER_REGULAR
         print '(A,ES9.2,A)', " ENGINE_XC_THEORY_DYNAMIC active (promote at dP <", &
            xcgrid_switch_prms, ")"
      endif
   endif
endif
if (xcgrid_use_theory_table) call xcgrid_theory_apply()

call get_environment_variable("ENGINE_XC_PRUNE_ALPHAS", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) (a(i), i = 1, 12)
   if (ios == 0) then
      do j = 1, 3
         do i = 1, 4
            xcgrid_prune_alphas(i, j) = a((j-1)*4 + i)
         enddo
      enddo
      print '(A)', " ENGINE_XC_PRUNE_ALPHAS override active"
   else
      print '(A)', " ENGINE_XC_PRUNE_ALPHAS: need 12 reals, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_PRUNE_ORDERS", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) (ord(i), i = 1, 5)
   if (ios == 0) then
      if (any(.not. xcgrid_lebedev_supported(ord))) then
         print '(A)', " ENGINE_XC_PRUNE_ORDERS: unsupported Lebedev order, ignoring"
      else
         do j = 1, 7
            xcgrid_prune_orders(:, j) = ord
         enddo
         print '(A,5I5)', " ENGINE_XC_PRUNE_ORDERS override active:", ord
      endif
   else
      print '(A)', " ENGINE_XC_PRUNE_ORDERS: need 5 integers, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_PRUNE_DMIN_BOUNDS", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) (dmb(i), i = 1, 3)
   if (ios == 0 .and. all(dmb(2:3) > dmb(1:2))) then
      xcgrid_prune_dmin_bounds = dmb
      print '(A,3F8.4)', " ENGINE_XC_PRUNE_DMIN_BOUNDS override active:", dmb
   else
      print '(A)', " ENGINE_XC_PRUNE_DMIN_BOUNDS: need 3 increasing reals, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_PRUNE_DMIN_ORDERS", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) (dmord(i), i = 1, 4)
   if (ios == 0) then
      if (any(.not. xcgrid_lebedev_supported(dmord))) then
         print '(A)', " ENGINE_XC_PRUNE_DMIN_ORDERS: unsupported Lebedev order, ignoring"
      else
         do j = 1, 7
            xcgrid_prune_dmin_orders(:, j) = dmord
         enddo
         print '(A,4I5)', " ENGINE_XC_PRUNE_DMIN_ORDERS override active:", dmord
      endif
   else
      print '(A)', " ENGINE_XC_PRUNE_DMIN_ORDERS: need 4 integers, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_RADPOT_SCALE", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) sc
   if (ios == 0 .and. sc > 0.0d0) then
      xcgrid_radpot_scale = sc
      print '(A,F8.4)', " ENGINE_XC_RADPOT_SCALE override active:", sc
   else
      print '(A)', " ENGINE_XC_RADPOT_SCALE: need one positive real, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_RADPOT", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) i
   if (ios == 0 .and. i > 0) then
      xcgrid_radpot_abs = i
      print '(A,I0)', " ENGINE_XC_RADPOT override active (absolute, tier-free): ", i
   else
      print '(A)', " ENGINE_XC_RADPOT: need a positive integer, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_RADPOT_ROWSCALE", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) sc
   if (ios == 0 .and. sc > 0.0d0) then
      xcgrid_radpot_rowscale = sc
      print '(A,F8.4)', " ENGINE_XC_RADPOT_ROWSCALE override active " // &
                        "(tier-free, x NRAD(row)):", sc
   else
      print '(A)', " ENGINE_XC_RADPOT_ROWSCALE: need a positive real, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_RADPOT_ROWSCALE_H", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) sc
   if (ios == 0 .and. sc > 0.0d0) then
      xcgrid_radpot_rowscale_h = sc
      print '(A,F8.4)', " ENGINE_XC_RADPOT_ROWSCALE_H override active " // &
                        "(row=1 only):", sc
   else
      print '(A)', " ENGINE_XC_RADPOT_ROWSCALE_H: need a positive real, ignoring"
   endif
endif

block
   integer :: zpairs(10), kk, nz
   real(8) :: zscales(10)
   call get_environment_variable("ENGINE_XC_RADPOT_ROWSCALE_Z", buf)
   if (len_trim(buf) > 0) then
      nz = 0
      zpairs = 0
      zscales = 0.0d0
      read(buf, *, iostat=ios) (zpairs(kk), zscales(kk), kk = 1, 10)
      do kk = 1, 10
         if (zpairs(kk) >= 1 .and. zpairs(kk) <= 103 .and. zscales(kk) > 0.0d0) then
            xcgrid_radpot_rowscale_z(zpairs(kk)) = zscales(kk)
            nz = nz + 1
         endif
      enddo
      if (nz > 0) print '(A,I0,A)', " ENGINE_XC_RADPOT_ROWSCALE_Z override active (", nz, " elements)"
   endif
end block

call get_environment_variable("ENGINE_XC_DUMP_SHELLS", buf)
xcgrid_dump_shells = (len_trim(buf) > 0 .and. trim(adjustl(buf)) == "1")

call get_environment_variable("ENGINE_XC_USE_ZONE_TABLE", buf)
if (len_trim(buf) > 0) xcgrid_use_zone_table = (trim(adjustl(buf)) == "1")

call get_environment_variable("ENGINE_XC_USE_DMIN_TABLE", buf)
if (len_trim(buf) > 0) xcgrid_use_dmin_table = (trim(adjustl(buf)) == "1")

call get_environment_variable("ENGINE_XC_USE_TREUTLER_AIJ", buf)
if (len_trim(buf) > 0) then
   xcgrid_use_treutler_aij = (trim(adjustl(buf)) == "1")
else
   xcgrid_use_treutler_aij = xcgrid_use_zone_table
endif

call get_environment_variable("ENGINE_XC_ZONE_NEARCORE", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) i, j
   if (ios == 0 .and. i > 0 .and. j > 0) then
      xcgrid_zone1 = i
      xcgrid_zone2 = j
      print '(A,2I5)', " ENGINE_XC_ZONE_NEARCORE override active:", i, j
   else
      print '(A)', " ENGINE_XC_ZONE_NEARCORE: need 2 positive integers, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_ZONE_FLANKS", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) i, j
   if (ios == 0 .and. i > 0 .and. j > 0) then
      xcgrid_zone3 = i
      xcgrid_zone5 = j
      print '(A,2I5)', " ENGINE_XC_ZONE_FLANKS override active:", i, j
   else
      print '(A)', " ENGINE_XC_ZONE_FLANKS: need 2 positive integers, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_ZONE_PEAK", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) i
   if (ios == 0 .and. i > 0) then
      xcgrid_zone4 = i
      print '(A,I5)', " ENGINE_XC_ZONE_PEAK override active:", i
   else
      print '(A)', " ENGINE_XC_ZONE_PEAK: need 1 positive integer, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_ZONE_PEAK_H", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) i
   if (ios == 0 .and. i > 0) then
      xcgrid_zone4_h = i
      print '(A,I5)', " ENGINE_XC_ZONE_PEAK_H override active:", i
   else
      print '(A)', " ENGINE_XC_ZONE_PEAK_H: need 1 positive integer, ignoring"
   endif
endif

block
   integer :: zpairs(10), zord(10), kk, nz
   call get_environment_variable("ENGINE_XC_ZONE_PEAK_Z", buf)
   if (len_trim(buf) > 0) then
      nz = 0
      zpairs = 0
      zord = 0
      read(buf, *, iostat=ios) (zpairs(kk), zord(kk), kk = 1, 10)
      do kk = 1, 10
         if (zpairs(kk) >= 1 .and. zpairs(kk) <= 103 .and. zord(kk) > 0) then
            xcgrid_zone4_z(zpairs(kk)) = zord(kk)
            nz = nz + 1
         endif
      enddo
      if (nz > 0) print '(A,I0,A)', " ENGINE_XC_ZONE_PEAK_Z override active (", nz, " elements)"
   endif
end block

call get_environment_variable("ENGINE_XC_RHO_DUMP", buf)
xcgrid_rho_dump = (len_trim(buf) > 0 .and. trim(adjustl(buf)) == "1")

call get_environment_variable("ENGINE_XC_ZONE_ALPHA3", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) sc
   if (ios == 0 .and. sc > 0.0d0) then
      xcgrid_zone_alpha3 = sc
      print '(A,F8.4)', " ENGINE_XC_ZONE_ALPHA3 override active:", sc
   else
      print '(A)', " ENGINE_XC_ZONE_ALPHA3: need a positive real, ignoring"
   endif
endif

call get_environment_variable("ENGINE_XC_GRID_AUDIT_ONLY", buf)
xcgrid_audit_only = (len_trim(buf) > 0 .and. trim(adjustl(buf)) == "1")

call get_environment_variable("ENGINE_XC_ANG_MARGIN", buf)
if (len_trim(buf) > 0) then
   read(buf, *, iostat=ios) i
   if (ios == 0 .and. i >= 0) then
      xcgrid_ang_margin = i
      print '(A,I0)', " ENGINE_XC_ANG_MARGIN override active: ", i
   else
      print '(A)', " ENGINE_XC_ANG_MARGIN: need a non-negative integer, ignoring"
   endif
endif
end subroutine xcgrid_tune_from_env

function xcgrid_lebedev_supported(ord) result(ok)
implicit none
integer, intent(in) :: ord(:)
logical :: ok(size(ord))
integer, parameter :: LEGAL(19) = (/ 6, 14, 26, 38, 50, 74, 86, 110, 146, &
                                     170, 194, 230, 266, 302, 350, 434, &
                                     590, 770, 974 /)
integer :: i
do i = 1, size(ord)
   ok(i) = any(LEGAL == ord(i))
enddo
end function xcgrid_lebedev_supported

end module GRID_info

