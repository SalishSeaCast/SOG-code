! $Id$
! $Source$

module fundamental_constants
  ! Parameter value declarations of fundamental constants widely used
  ! in other modules of the SOG code.
  !
  ! Public Parameters:
  !
  !   g -- Acceleration due to gravity [m/s^2]
  !
  !   latitude -- Station S3 latitude [deg]
  !
  !   pi -- Ratio of circumference to diameter of a circle [-]

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       g,        &  ! Acceleration due to gravity [m/s^2]
       latitude, &  ! Station S3 latitude [deg]
       pi           ! Ratio of circumference to diameter of a circle [-]

  ! Parameter Value Declarations:
  !
  ! Public:
  real(kind=dp), parameter :: &
       g = 9.80665, &                  ! Acceleration due to gravity [m/s^2]
       latitude = 49. + 7.517 / 60., & ! Station S3 latitude [deg]
       pi = 3.141592653589793

end module fundamental_constants
