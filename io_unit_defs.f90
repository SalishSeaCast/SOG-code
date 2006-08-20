! $Id$
! $Source$

module io_unit_defs
  ! Input/output unit definitions

  implicit none

  ! Unix standards (set to g95 default values)
  integer, parameter :: stdin = 5
  integer, parameter :: stdout = 6
  integer, parameter :: stderr = 0

  ! Input data file units
  integer, parameter :: met_data = 12
  integer, parameter :: river_data = 12

end module io_unit_defs
