module irradiance
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to light and heat influences in SOG code.
  !
  ! Public Subroutines:
  !
  !    init_irradiance -- Allocate memory for parameters for Kpar fit
  !                       and read parameter values from infile.
  !
  !    calc_irradiance -- Contains cloud model which calculates the
  !                       light filtered from sun, and the
  !                       parameterizations for the deposition of
  !                       photosynthetic available radiation (PAR) and
  !                       total light available for water heating.
   
  use precision_defs, only: dp, sp
  use grid_mod, only: interp_g
  use fundamental_constants, only: pi, latitude
  use core_variables, only: &
       N2chl ! ratio of chl mg/m3 to N uMol for phytoplankton
  use input_processor, only: getpard

  implicit none

  private
  public :: &
       ! Variables:
       I_par,   & ! Photosynthetic available radiation (PAR) at grid points
       I_total, & ! Total light available for water heating at grid interfaces
       Q_n,     & ! Non-turbulent heat flux profile
       ! Subroutines:
       init_irradiance, calc_irradiance, dalloc_irradiance_variables
  ! Parameter Value Declarations:
  !
  ! Private to module:
  real(kind=dp), parameter :: &
       Q_o = 1368.0, & !1367.0? W/m^2  Solar constant
       albedo = 0.18   !KC 17% OCT.22 2004
  !
  ! Variable Declarations:
  !
  ! Public:
  real(kind=dp), dimension(:), allocatable :: &
       I_par, &    ! Photosynthetic available radiation (PAR) at grid points
       I_total, &  ! Total light available for heating at grid interfaces
       Q_n         ! Non-turbulent heat flux profile
  !
  ! Private to module:
  real(kind=dp) :: &     ! Values for the Kpar fit
       ialpha, ibeta, igamma, isigma, itheta, idl
  real(kind=dp), dimension(:), allocatable :: &
       I_par_i  ! Photosynthetic available radiation (PAR) at grid interfaces

contains

  subroutine init_irradiance(M)
    ! Allocate memory for Kpar variables, and read parameter values
    ! from infile
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Allocate memory for upwelling quantity arrays
    call alloc_irradiance_variables(M)
    ! Read Kpar parameter values from infile
    call read_irradiance_params()
  end subroutine init_irradiance


  subroutine read_irradiance_params
    ! Read the Kpar parameter values from the infile
    use input_processor, only: getpard
    implicit none 

    ! Read values for Kpar fit- Collins et al 2008 equation 8.
    ialpha = getpard('ialpha')
    ibeta = getpard('ibeta')
    igamma = getpard('igamma')
    isigma = getpard('isigma')
    itheta = getpard('itheta')
    idl = getpard('idl')
  end subroutine read_irradiance_params


  subroutine calc_irradiance(cf, day_time, day, Qriver, &
       Pmicro, Pnano, Ppico)
    ! Calculates the amount of light (total and PAR) as a function of depth
    use grid_mod, only: grid
    implicit none
    ! Arguments:
    real(kind=sp), intent(in) :: cf        ! cloud fraction
    real(kind=dp), intent(in) :: day_time  ! day-second
    integer, intent(in) :: day             ! year-day 
    real(kind=dp), intent(in) :: Qriver    ! major river flow
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano, &   ! Nano phytoplankton
         Ppico      ! Pico phytoplankton
    
    ! Local:
    ! Type declarations:
    type :: cloudy
       real(kind=dp) :: A, B    !Regression coefficients
    end type cloudy
    type :: tenths
       type(cloudy), dimension(0:10) :: type
    end type tenths
    ! Variable declarations:
    type(tenths) :: cloud
    real(kind=dp) :: lat  ! Latitude of centre of model domain [rad]
    integer :: k, fcf, ccf
    real(kind=dp) :: declination, hour, cos_Z, day_length, hour_angle, &
         sunrise, sunset, Qso, a, b, KK
    real(kind=dp):: I_incident
    real(kind=dp) :: Ired, Iblue ! light in red and blue parts of spectrum
    
    ! Cloud model based on Dobson and Smith, table 5
    !
    !!! SEA -- May 2010 : redid the cloud parameterization based on UBC
    !!! Solar data (/ocean/shared/SoG/met/solar/) fitting Q to cos_Z
    !!! (not Q/cos_Z as Kate did).  Page 61 in Lab Book.  (0) no
    !!! clouds, (1) 1/10 cloud fraction (10) 100% clouds.  Four sig
    !!! figs are what comes out of matlab but standard deviations are
    !!! 40W/m2 for low cloud fraction to 120 W/m2 for 6-9 cloud
    !!! fraction to 85 W/m2 for completely cloudy.
    cloud%type(0) = cloudy(0.6337, 0.1959) 
    cloud%type(1) = cloudy(0.6149, 0.2119)  
    cloud%type(2) = cloudy(0.5861, 0.2400)
    cloud%type(3) = cloudy(0.5512, 0.2859)
    cloud%type(4) = cloudy(0.5002, 0.3192)
    cloud%type(5) = cloudy(0.4649, 0.3356)
    cloud%type(6) = cloudy(0.4225, 0.3339)
    cloud%type(7) = cloudy(0.3669, 0.3490)
    cloud%type(8) = cloudy(0.2468, 0.4427)
    cloud%type(9) = cloudy(0.1981, 0.3116)
    cloud%type(10) = cloudy(0.0841, 0.2283)


    hour = (day_time / 3600.0d0 - 12.0d0) * 15.0d0  !degrees
    declination = 23.45*PI/180.0*SIN((284.0+DBLE(day))/365.25*2.0*PI)  !radians

    ! Convert latitude of centre of model domain from degrees to radians
    lat = pi * latitude / 180.0d0
    
    a = sin(declination) * sin(lat) 
    b = cos(declination) * cos(lat)
    cos_Z = a + b * cos(pi / 180.0 * hour)      !solar elevation
    hour_angle = tan(lat)*tan(declination)  ! cos of -hour_angle in radians
    if (hour_angle > 1) then   ! so far North in summer that there is no night
       day_length = 24.+0.0001
    elseif (hour_angle < -1) then ! so far North in winter that there is no day
       day_length = -0.0001
    else
       day_length = acos(-hour_angle) / 15.0 * 2.0 * 180.0 / pi ! hours: 15 = 360/24
    endif
    sunrise = 12.0 - 0.5 * day_length  !hours
    sunset = 12.0 + 0.5 * day_length   !hours

    Qso = Q_o*(1.0+0.033*COS(DBLE(day)/365.25*2.0*PI))*(1.0-albedo) !*(1.0-insol)

    fcf = floor(cf)   ! integer below cf value
    ccf = ceiling(cf) ! integer above cf value
    if (fcf.eq.ccf) then
       if (fcf.eq.10) then
          fcf = 9
       else
          ccf = fcf+1
       endif
    endif

    if (day_time / 3600.0 > sunrise .and. day_time / 3600.0 < sunset) then
       I_incident = Qso * (  &
            cloud%type(fcf)%A * (ccf-cf) + cloud%type(ccf)%A * (cf-fcf) &
            + (cloud%type(fcf)%B * (ccf-cf) + cloud%type(ccf)%B * (cf-fcf)) &
            * cos_Z) * cos_Z   
    else
       I_incident = 0.
    end if

    ! Light is defined on interfaces
    ! Parameterized for a particular basin
    !
    ! Parameterization of Dec 2006 by S.Allen for the Strait of
    ! Georgia based on SOG data:
    !  * KK is evaluated on the grid points         
    !  * N2chl is correction uM to mg/m3 chl
    !  * limit of 2.5/m is set from Cruise 02-03 and 5 W/m2 seen at 2 m
    !  * re-fit Jun 13, 2007 pg 113-114 lab-book and in fitlight
    !  * 1.5d0 is correction for CTD fluoresence which seems high
    !
    ! Parameterization of Mar 2009 by M. Wolfe for Rivers Inlet based
    ! on RI 2007 hydrolab data:
    !  * limit of 2.5/m is set from Cruise 02-03 and 5 W/m2 seen at 2 m
    I_total(0) =  I_incident
    ! 44% of total light is PAR at surface (Jerlov)
    I_par_i(0) = I_incident * 0.44
    Iblue = I_total(0) * 0.7
    Ired = I_total(0) * 0.3
    do k = 1, grid%M    
       KK = ialpha +ibeta * (N2chl &
            * (Pmicro(k) + Pnano(k) + Ppico(k))) ** 0.665 &
            + (igamma * Qriver ** isigma + itheta) * exp(-grid%d_g(k) / idl)
       KK = min(2.5d0, KK)
       I_par_i(k) = I_par_i(k-1) * exp(-grid%i_space(k) * KK)
       ! Total light for heat budget
       Iblue = Iblue * exp(-grid%i_space(k) * (0.8102d0 * KK + 1.1854d0))
       Ired = Ired * exp(-grid%i_space(k) * (0.8226d0 * KK - 0.0879d0))
       I_total(k) = Iblue + Ired
    end do
    I_par(0) = I_par_i(0)
    I_par(1:) = interp_g(I_par_i)
   end subroutine calc_irradiance


   subroutine alloc_irradiance_variables(M)
     ! Allocate memory for irradiance quantity arrays
     use malloc, only: alloc_check
     implicit none
     ! Argument:
     integer, intent(in) :: M  ! Number of grid points
     ! Local variables:
     integer           :: allocstat  ! Allocation return status
     character(len=80) :: msg        ! Allocation failure message prefix

     msg = "Vertical profiles of uwelling velocity arrays"
     allocate(I_par(0:M), I_par_i(0:M), I_total(0:M), Q_n(0:M), &
          stat=allocstat) 
     call alloc_check(allocstat, msg)
   end subroutine alloc_irradiance_variables


   subroutine dalloc_irradiance_variables
     ! Deallocate memory for irradiance quantity arrays
     use malloc, only: dalloc_check
     implicit none
     ! Local variables:
     integer           :: dallocstat  ! Deallocation return status
     character(len=80) :: msg         ! Deallocation failure message prefix

     msg = "Vertical profiles of upwelling velocity arrays"
     deallocate(I_par, I_par_i, I_total, Q_n, &
          stat=dallocstat)
     call dalloc_check(dallocstat, msg)
   end subroutine dalloc_irradiance_variables

end module irradiance
