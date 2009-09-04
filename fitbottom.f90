! $Id$
! $Source$

module fitbottom
  ! Constants, functions, and subroutines for calculation of boundary
  ! condition values at bottom of grid.
  !
  ! Public Subroutines:
  !
  ! bot__bound_time (day, day_time, Tbot, Sbot, Nobot, Sibot, Pmbot, Pnbot)
  !  -- Calculate the values at the bottom of the grid for those
  !     quantities that we have data for from an annual fit.
  !
  ! bot_bound_uniform (M, Z, NH, Detritus)
  !  -- Set the values at the bottom of the grid for those quantities
  ! that we don't have time series data for.

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Subroutines:
       init_fitbottom, bot_bound_time, bot_bound_uniform

  ! Private module variable declarations:
  ! 
  ! Number of different quantities
  integer, parameter :: NQ = 7
  ! Quantity names
  character(len=3), dimension(NQ), parameter :: quantity &
       = ['sal', 'tmp', 'chl', 'nit', 'sil', 'NH4', 'prt']

    ! Unless otherwise specified: tit coefficients (see fitbottom.py and 
    ! fitted.m in
    ! sog-forcing/bottom for details of how these coefficient values
    ! were derived. - Need to read in values from infile
    real(kind=dp), dimension(7, NQ) :: c 

    ! Flag for constant bottom temperature
    logical :: temp_constant  !is the temperature coming up through the 
                              !bottom constant?
     
contains

  subroutine init_fitbottom()
    ! Read the bottom temperature constraints
    use input_processor, only: getparl, getpardv
    implicit none

    ! Is the temperature coming up through the bottom constant?
    temp_constant = getparl('temp_constant')
    ! Read in values for sal, temp, chl, nit, sil, NH4, prt and set up
    ! c matrix
    call getpardv("salinity", 7, c(:,1))
    call getpardv("temperature", 7, c(:,2))
    call getpardv("Phytoplankton", 7, c(:,3))
    call getpardv("Nitrate", 7, c(:,4))
    call getpardv("Silicon", 7, c(:,5))
    call getpardv("Ammonium", 7, c(:,6))
    call getpardv("Ratio", 7, c(:,7))
  end subroutine init_fitbottom

  function bottom_value (arg, qty) result(value)
    ! Calculate the value of the specified quantitity at the bottom
    ! grid boundary for the given value of arg.
    use precision_defs, only: dp
    use io_unit_defs, only: stdout
    use fundamental_constants, only: pi
    
    implicit none
    ! Arguments:
    real(kind=dp), intent(in):: arg
    character(len=*), intent(in):: qty
    ! Result:
    real(kind=dp) :: value
    ! Local variable:
    integer:: index
    
    ! Set the coefficient matrix 2nd index for the quantity that we've
    ! calculating
    if (qty == quantity(1)) then
       index = 1
    elseif (qty == quantity(2)) then
       index = 2
       ! Apply constant temperature bottom boundary condition, if requested
       if(temp_constant) then
          c(2,index) = 0.d0
          c(3,index) = 0.d0
          c(4,index) = 0.d0
          c(5,index) = 0.d0
          c(6,index) = 0.d0
          c(7,index) = 0.d0
       endif
    elseif (qty == quantity(3)) then
       index = 3
    elseif (qty == quantity(4)) then
       index = 4
    elseif (qty == quantity(5)) then
       index = 5
    elseif (qty == quantity(6)) then
       index = 6
    elseif (qty == quantity(7)) then
       index = 7
    else
       write (stdout,*) 'bottom_value in fitbottom.f90:', &
            'Unexpected quantity: ', qty
       call exit(1)
    endif
    ! Calculate the bottom boundary condition value
    value = c(1,index) &
         + c(2,index) * cos(arg) + c(3,index) * sin(arg + c(4,index)) &
         + c(5,index) * cos(3.0d0 * arg) &
         + c(6,index) * sin(3.0d0 * arg) &
         + c(7,index)*((arg*365.25) /( 2 * pi))
  end function bottom_value


  subroutine bot_bound_time (year, day, day_time, &
       Tbot, Sbot, Nobot, Sibot, Nhbot, Pmbot, Pnbot, Ppbot, uZbot)
    ! Calculate the values at the bottom of the grid for those
    ! quantities that we have data for from an annual fit.
    use precision_defs, only: dp
    use unit_conversions, only: CtoK
    use fundamental_constants, only: pi
    use forcing, only: vary_forcing, vary
    implicit none
    ! Arguments:
    integer, intent(in) :: year, day
    real(kind=dp), intent(in) :: day_time
    real(kind=dp), intent(out) :: Tbot, Sbot, Nobot, Sibot, Nhbot, Pmbot, &
         Pnbot, Ppbot, uZbot
    ! Local variables:
    real(kind=dp) :: arg, chl, ratio, yeartime

    yeartime = (year - 2002) + day / 365.25d0 + day_time / 86400.0d0 / 365.25d0
    arg = 2.0d0 * pi * yeartime
    Tbot = bottom_value(arg, 'tmp')          ! in Celcius
    if (vary%temperature%enabled .and. .not. vary%temperature%fixed) then
       Tbot = CtoK(Tbot) + vary%temperature%addition
    else
       Tbot = CtoK(Tbot)
    endif
    Sbot = bottom_value(arg, 'sal')
    Nobot = bottom_value(arg, 'nit')
    Sibot = bottom_value(arg, 'sil')
    Nhbot = bottom_value(arg, 'NH4')
    chl = bottom_value(arg, 'chl')
    ratio = bottom_value(arg, 'prt')
    Pmbot = chl / (ratio + 1)
    Pnbot = ratio * Pmbot
    Ppbot = Pnbot
    uZbot = Pnbot
  end subroutine bot_bound_time


  subroutine bot_bound_uniform (M, D_DON, D_PON, D_refr, D_bSi)
    ! Set the values at the bottom of the grid for those
    ! quantities that we don't have time series data for.
    use precision_defs, only: dp
    implicit none
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    real(kind=dp), intent(inout), dimension(0:) :: &
         D_DON,  &  ! Dissolved organic nitrogen concentration profile array
         D_PON,  &  ! Particulate organic nitrogen concentration profile array
         D_refr, &  ! Refractory nitrogen concentration profile array
         D_bSi      ! Biogenic silicon concentration profile array

    D_DON(M+1) = D_DON(M)
    D_PON(M+1) = D_PON(M)
    D_refr(M+1) = D_refr(M)
    D_bSi(M+1) = D_bSi(M)
  end subroutine bot_bound_uniform
  
end module fitbottom
