! $Id$
! $Source$

module unit_conversions
  ! A collection of functions that do unit conversions:
  !  CtoK -- degrees Celcius to Kelvins
  !  KtoC -- Kelvins to degrees Celcius

  use precision_defs

  implicit none

  interface CtoK
     module procedure CtoK_sp, CtoK_dp, CtoK_array
  end interface
  interface KtoC
     module procedure KtoC_sp, KtoC_dp, KtoC_array
  end interface

  private
  public :: CtoK, KtoC

contains

  real(kind=sp) function CtoK_sp(degC) result(K)
    ! Convert temperature from Celcius to Kelvin
    implicit none
    real(kind=sp), intent(in)  :: degC

    K = degC + 273.15
  end function CtoK_sp


  real(kind=dp) function CtoK_dp(degC) result(K)
    ! Convert temperature from Celcius to Kelvin
    implicit none
    real(kind=dp), intent(in)  :: degC

    K = degC + 273.15
  end function CtoK_dp


  function CtoK_array(degC) result(K)
    ! Convert temperature from Celcius to Kelvin
    ! kind=dp array version
    implicit none
    real(kind=dp), dimension(:), intent(in)  :: degC
    real(kind=dp), dimension(size(degC)) :: K

    K = degC + 273.15
  end function CtoK_array


  real(kind=sp) function KtoC_sp(K) result(degC)
    ! Convert temperature from Kelvin to Celcius
    ! kind=sp version
    implicit none
    real(kind=sp), intent(in)  :: K

    degC = K - 273.15
  end function KtoC_sp


  real(kind=dp) function KtoC_dp(K) result(degC)
    ! Convert temperature from Kelvin to Celcius
    ! kind=dp version
    implicit none
    real(kind=dp), intent(in)  :: K

    degC = K - 273.15
  end function KtoC_dp


  function KtoC_array(K) result(degC)
    ! Convert temperature from Kelvin to Celcius
    ! kind=dp array version
    implicit none
    real(kind=dp), dimension(:), intent(in)  :: K
    real(kind=dp), dimension(size(K)) :: degC

    degC = K - 273.15
  end function KtoC_array

end module unit_conversions
