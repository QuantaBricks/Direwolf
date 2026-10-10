! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_exchange: public-facing entry point for HF exact-exchange (K

module mod_exchange
use mod_integrals, only: integrals_build_exchange, integrals_build_exchange_lr, &
                          integrals_build_exchange_df, integrals_build_exchange_df_lr, &
                          integrals_build_exchange_cosx, integrals_build_exchange_cosx_lr, &
                          integrals_build_exchange_cosx_sr, &
                          integrals_force_exchange_cosx
implicit none
private
public :: HF_exchange_frac
public :: RS_omega
public :: RS_beta
public :: cosx_enabled
public :: exchange_build
public :: exchange_build_cosx
public :: exchange_build_cosx_lr
public :: exchange_build_cosx_sr
public :: exchange_build_lr
public :: exchange_build_df
public :: exchange_build_df_lr
public :: exchange_force_cosx

logical :: cosx_enabled = .false.

real(8) :: HF_exchange_frac = 0.0d0
real(8) :: RS_omega = 0.0d0
real(8) :: RS_beta  = 0.0d0

contains

subroutine exchange_build(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    call integrals_build_exchange(nConts, Da, Db, Ka, Kb)
end subroutine exchange_build

subroutine exchange_build_cosx(nConts, Da, Db, Ka, Kb, need_force)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    logical,intent(in) :: need_force
    call integrals_build_exchange_cosx(nConts, Da, Db, Ka, Kb, need_force)
end subroutine exchange_build_cosx

subroutine exchange_build_cosx_lr(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    call integrals_build_exchange_cosx_lr(nConts, Da, Db, Ka, Kb)
end subroutine exchange_build_cosx_lr

subroutine exchange_build_cosx_sr(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    call integrals_build_exchange_cosx_sr(nConts, Da, Db, Ka, Kb)
end subroutine exchange_build_cosx_sr

subroutine exchange_force_cosx(nConts, natoms, Da, Db, RSomega, forceKa, forceKb)
    implicit none
    integer,intent(in) :: nConts, natoms
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(in) :: RSomega
    real(8),intent(out) :: forceKa(natoms,3), forceKb(natoms,3)
    call integrals_force_exchange_cosx(nConts, natoms, Da, Db, RSomega, forceKa, forceKb)
end subroutine exchange_force_cosx

subroutine exchange_build_lr(nConts, Da, Db, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Da(nConts,nConts), Db(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    call integrals_build_exchange_lr(nConts, Da, Db, Ka, Kb)
end subroutine exchange_build_lr

subroutine exchange_build_df(nConts, Ca, nOccA, Cb, nOccB, Ka, Kb)
    implicit none
    integer,intent(in) :: nConts, nOccA, nOccB
    real(8),intent(in) :: Ca(nConts,nConts), Cb(nConts,nConts)
    real(8),intent(out) :: Ka(nConts,nConts), Kb(nConts,nConts)
    call integrals_build_exchange_df(nConts, Ca, nOccA, Cb, nOccB, Ka, Kb)
end subroutine exchange_build_df

subroutine exchange_build_df_lr(nConts, Ca, nOccA, Cb, nOccB, Ka_lr, Kb_lr)
    implicit none
    integer,intent(in) :: nConts, nOccA, nOccB
    real(8),intent(in) :: Ca(nConts,nConts), Cb(nConts,nConts)
    real(8),intent(out) :: Ka_lr(nConts,nConts), Kb_lr(nConts,nConts)
    call integrals_build_exchange_df_lr(nConts, Ca, nOccA, Cb, nOccB, Ka_lr, Kb_lr)
end subroutine exchange_build_df_lr

end module mod_exchange
