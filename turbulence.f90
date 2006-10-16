! $Id$
! $Source$

module turbulence
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the KPP turbulence parameterization in the
  ! SOG code.
  !
  ! Public Parameters:
  !
  !   Constants for the nondimensional flux profiles.  See Large, et al
  !   (1994) pg 392, eqn B2
  !      xsi_m -- Maximum xsi value of the -1/3 power law regime of phi%m
  !      xsi_s -- Maximum xsi value of the -1/3 power law regime of phi%s
  !      a_m --  Coefficient of phi%m in 1/3 power law regime
  !      a_s --  Coefficient of phi%s in 1/3 power law regime
  !      c_m  -- Coefficient of phi%m in 1/3 power law regime
  !      c_s --  Coefficient of phi%s in 1/3 power law regime
  !
  ! Public Variables:
  !
  ! Public Subroutines:
  !

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       xsi_m, &  ! Max xsi value of the -1/3 power law regime of phi%m
       xsi_s, &  ! Max xsi value of the -1/3 power law regime of phi%s
       a_m,   &  ! Coefficient of phi%m in 1/3 power law regime
       a_s,   &  ! Coefficient of phi%s in 1/3 power law regime
       c_m,   &  ! Coefficient of phi%m in 1/3 power law regime
       c_s       ! Coefficient of phi%s in 1/3 power law regime
       ! Variables:
       ! Subroutines:

  ! Type Definitions:
  !
  ! Public:
  !
  ! Private to module:

  ! Parameter Value Declarations:
  !
  ! Public:
  real(kind=dp), parameter :: &
       xsi_m = -0.20, &  ! Max xsi value of the -1/3 power law regime of phi%m
       xsi_s = -1.0,  &  ! Max xsi value of the -1/3 power law regime of phi%s
       a_m = 1.26,    &  ! Coefficient of phi%m in 1/3 power law regime
       a_s = -28.86,  &  ! Coefficient of phi%s in 1/3 power law regime
       c_m = 8.38,    &  ! Coefficient of phi%m in 1/3 power law regime
       c_s = 98.96       ! Coefficient of phi%s in 1/3 power law regime
  !
  ! Private to module:

!!$contains

end module turbulence
