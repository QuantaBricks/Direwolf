! Copyright (c) 2026 QuantaBricks
! SPDX-License-Identifier: AGPL-3.0-or-later

! Default XC integration-grid tier, per functional.

module xcgrid_defaults
implicit none

integer, parameter :: XCGRID_COARSE = 2
integer, parameter :: XCGRID_FINE   = 3
integer, parameter :: XCGRID_L4     = 4
integer, parameter :: XCGRID_L5     = 5
integer, parameter :: XCGRID_L6     = 6
integer, parameter :: XCGRID_L7     = 7

contains

integer function xcgrid_default_level(functional) result(lvl)
character(len=*), intent(in) :: functional
lvl = XCGRID_FINE
select case (trim(functional))
case ("R2SCAN-3C", "R2SCAN_3C", "R2SCAN3C")
   lvl = XCGRID_L4
case ("R2SCAN", "R2SCAN0", &
      "M06L", "M06-L", "M06_L", "M06", "M06_HYB", "M06-2X", "M06-2X_HYB", "M05-2X", &
      "MN15", "MN15L", "MN15-L")
   lvl = XCGRID_L5
end select
end function xcgrid_default_level

end module xcgrid_defaults
