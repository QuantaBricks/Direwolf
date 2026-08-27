! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_dispersion_d4: Grimme's -D4 dispersion correction (E. Caldeweyher,

module mod_dispersion_d4
use dftd4, only: d4_model, new_d4_model, get_rational_damping, get_dispersion, &
                  realspace_cutoff, damping_param, rational_damping_param
use mctc_io, only: structure_type, new
use mctc_env, only: error_type
implicit none
private
public :: dispersion_d4

contains

subroutine dispersion_d4(natoms, atomic_number, coor_bohr, icharge, functional, Edisp, dEdisp)
implicit none
integer,intent(in) :: natoms
integer,intent(in) :: atomic_number(natoms)
real(8),intent(in) :: coor_bohr(3,natoms)
integer,intent(in) :: icharge
character(*),intent(in) :: functional
real(8),intent(out) :: Edisp
real(8),intent(out) :: dEdisp(3,natoms)

type(structure_type) :: mol
type(d4_model) :: disp
class(damping_param),allocatable :: param
type(error_type),allocatable :: error
real(8) :: sigma_unused(3,3)

call new(mol, atomic_number, coor_bohr, charge=real(icharge,8))
call new_d4_model(error, disp, mol)
if (allocated(error)) then
   print '(A,A)', " Warning: dftd4's D4 model setup failed for this geometry - ", &
         "dispersion correction skipped for this run."
   Edisp = 0.0d0
   dEdisp = 0.0d0
   return
endif

call get_rational_damping(trim(functional), param)
if (.not. allocated(param)) then
   print '(A,A,A)', " Warning: dftd4 has no published fit for functional '", &
         trim(functional), "' - dispersion correction skipped for this run."
   Edisp = 0.0d0
   dEdisp = 0.0d0
   return
endif

call get_dispersion(mol, disp, param, realspace_cutoff(), Edisp, dEdisp, sigma_unused)

end subroutine dispersion_d4

end module mod_dispersion_d4
