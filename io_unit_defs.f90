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
  integer, parameter :: stripped_infile = 10
  integer, parameter :: met_data = 12
  integer, parameter :: river_data = 12

  ! Results output file units
  integer, parameter :: profiles = 200
  integer, parameter :: haloclines = 201
  integer, parameter :: Hoffmueller = 202
  integer, parameter :: std_phys_timeseries = 300
  integer, parameter :: user_phys_timeseries = 301
  integer, parameter :: std_bio_timeseries = 400
  integer, parameter :: user_bio_timeseries = 401

end module io_unit_defs
