module fundamental_constants
  ! Parameter value declarations of fundamental constants widely used
  ! in other modules of the SOG code.
  !
  ! Public Parameters:
  !
  !   f -- Coriolis factor
  !
  !   g -- Acceleration due to gravity [m/s^2]
  !
  !   latitude -- Latitude of location being modeled [deg]
  !
  !   pi -- Ratio of circumference to diameter of a circle [-]

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       f,        &  !  Coriolis factor
       g,        &  ! Acceleration due to gravity [m/s^2]
       latitude, &  ! Latitude of location being modeled [deg]
       pi,       &  ! Ratio of circumference to diameter of a circle [-]
       ! Subroutine:
       init_constants

  ! Public parameter declarations:
  real(kind=dp) :: &
       f, &      ! Coriolis factor
       latitude  ! latitude [deg]
  real(kind=dp), parameter :: &
       g = 9.80665d0, &                  ! Acceleration due to gravity [m/s^2]
       pi = 3.141592653589793d0

contains

  subroutine init_constants()
    use input_processor, only: getpard
    implicit none

    ! Latitude of location being modeled from infile
    latitude = getpard('latitude')
    ! Coriolis factor
    f = 2.0d0 * (2.0d0 * pi / 86400.0d0) * sin(pi * latitude / 180.0d0)
  end subroutine init_constants
  
end module fundamental_constants
