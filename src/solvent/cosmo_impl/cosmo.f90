! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: Apache-2.0

! mod_cosmo: top-level COSMO driver - ties mod_cosmo_cavity (pure

module mod_cosmo
use mod_cosmo_state, only: cosmo_enabled
use mod_cosmo_init, only: cosmo_init
use mod_cosmo_scf, only: cosmo_scf_step, cosmo_report_sigma_profile, cosmo_finalize, &
                         cosmo_set_sigma_profile_debug
use mod_cosmo_force, only: cosmo_force_step
use mod_cosmo_solvents, only: cosmo_solvent_epsilon, cosmo_solvent_smd_params
implicit none
private
public :: cosmo_enabled
public :: cosmo_init
public :: cosmo_scf_step
public :: cosmo_report_sigma_profile
public :: cosmo_set_sigma_profile_debug
public :: cosmo_force_step
public :: cosmo_finalize
public :: cosmo_solvent_epsilon
public :: cosmo_solvent_smd_params

end module mod_cosmo
