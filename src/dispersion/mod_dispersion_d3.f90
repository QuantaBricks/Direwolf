! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_dispersion_d3: Grimme's -D3 dispersion correction (S. Grimme, J.

module mod_dispersion_d3
use dftd3, only: d3_model, new_d3_model, d3_param, get_rational_damping, &
                  get_zero_damping, get_dispersion, realspace_cutoff, &
                  damping_param, rational_damping_param, new_rational_damping, &
                  zero_damping_param, new_zero_damping
use mctc_io, only: structure_type, new
use mctc_env, only: error_type
implicit none
private
public :: dispersion_d3

contains

function normalize_functional_name(functional) result(fname)
implicit none
character(*),intent(in) :: functional
character(len=:),allocatable :: fname
select case (trim(functional))
case ("PBE_PBE");    fname = "PBE"
case ("B3LYP_HYB");  fname = "B3LYP"
case ("CAM_B3LYP");  fname = "CAM-B3LYP"
case default;        fname = trim(functional)
end select
end function normalize_functional_name

subroutine dispersion_d3(natoms, atomic_number, coor_bohr, functional, damping, Edisp, dEdisp)
implicit none
integer,intent(in) :: natoms
integer,intent(in) :: atomic_number(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
character(*),intent(in) :: functional, damping
real(8),intent(out) :: Edisp
real(8),intent(out) :: dEdisp(3,natoms)

type(structure_type) :: mol
type(d3_model) :: disp
type(d3_param) :: param
class(damping_param),allocatable :: dmp
type(error_type),allocatable :: error
real(8) :: sigma_unused(3,3)

call new(mol, atomic_number, coor_bohr)
call new_d3_model(disp, mol)

if (trim(damping) .eq. 'ZERO') then
   call get_zero_damping(param, normalize_functional_name(functional), error)
else
   call get_rational_damping(param, normalize_functional_name(functional), error)
endif

if (allocated(error)) then
   print '(A,A,A,A,A)', " Warning: simple-dftd3 has no published DFT-D3(", &
         trim(damping), ") fit for functional '", trim(functional), &
         "' - dispersion correction skipped for this run."
   Edisp = 0.0d0
   dEdisp = 0.0d0
   return
endif

if (trim(damping) .eq. 'ZERO') then
   block
   type(zero_damping_param) :: zd
   call new_zero_damping(zd, param)
   allocate(dmp, source=zd)
   end block
else
   block
   type(rational_damping_param) :: rd
   call new_rational_damping(rd, param)
   allocate(dmp, source=rd)
   end block
endif

call get_dispersion(mol, disp, dmp, realspace_cutoff(), Edisp, dEdisp, sigma_unused)

end subroutine dispersion_d3

end module mod_dispersion_d3
