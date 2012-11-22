module fundamental_constants
  ! Parameter value declarations of fundamental constants widely used
  ! in other modules of the SOG code.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       f,           &  !  Coriolis factor
       g,           &  ! Acceleration due to gravity [m/s^2]
       latitude,    &  ! Latitude of location being modeled [deg]
       pCO2_atm,    &  ! Partial pressure of atmospheric CO2 [atm]
       pi,          &  ! Ratio of circumference to diameter of a circle [-]
       R_gas,       &  ! Gas constant
       Redfield_C,  &  ! Biological uptake ratio Carbon to Nitrogen
       Redfield_O,  &  ! Biological uptake ratio O2 to Nitrogen
       Redfield_NP, &  ! Biological uptake ratio Nitrogen to Phosphorus
       ! Subroutine:
       init_constants

  ! Public parameter declarations:
  real(kind=dp) :: &
       f,          &  ! Coriolis factor
       latitude,   &  ! latitude [deg]
       pCO2_atm       ! Partial pressure of atmospheric CO2 [atm]
  real(kind=dp), parameter :: &
       g = 9.80665d0,            &    ! Acceleration due to gravity [m/s^2]
       pi = 3.141592653589793d0, &
       R_gas = 83.1451d0,        &    ! Gas const [ml bar-1 K-1 mol-1] DOEv2
       Redfield_C = 6.625d0,     &    ! 106:16:1:-150 (C:N:P:O2)
       Redfield_O = 9.375d0,     &    ! 106:16:1:-150 (C:N:P:O2)
       Redfield_NP = 16.0d0           ! 106:16:1:-150 (C:N:P:O2)

contains

  subroutine init_constants()
    use input_processor, only: getpard
    implicit none

    latitude = getpard('latitude')
    f = 2.0d0 * (2.0d0 * pi / 86400.0d0) * sin(pi * latitude / 180.0d0)
    pCO2_atm = getpard('pCO2_atm')
  end subroutine init_constants

end module fundamental_constants
