! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! mod_density_fitting: thin public wrapper around mod_integrals' RI-J

module mod_density_fitting
use mod_integrals, only: integrals_build_coulomb_df
implicit none
private
public :: df_build_coulomb

contains

subroutine df_build_coulomb(nConts, Ptot, J)
    implicit none
    integer,intent(in) :: nConts
    real(8),intent(in) :: Ptot(nConts,nConts)
    real(8),intent(out) :: J(nConts,nConts)
    call integrals_build_coulomb_df(nConts, Ptot, J)
end subroutine df_build_coulomb

end module mod_density_fitting
