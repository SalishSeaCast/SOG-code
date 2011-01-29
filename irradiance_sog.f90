module irradiance
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to light and heat influences in SOG code.
  !
  ! Public Subroutines:
  !
  !    init_irradiance -- Allocate memory for parameters for Kpar fit
  !                       and read parameter values from infile.
  !
  !    irradiance_sog -- Contains cloud model which calculates the
  !                      light filtered from sun, and the
  !                      parameterizations for the deposition of PAR
  !                      and total light available for heat.
   
  use precision_defs, only: dp, sp
  use grid_mod, only: grid_, interp_g
  use fundamental_constants, only: pi, latitude
  use core_variables, only: N2chl ! ratio of chl mg/m3 to N uMol for
                                  ! phytoplankton
  use input_processor, only: getpard

  implicit none

  private
  public :: &
       ! Subroutines:
       init_irradiance, irradiance_sog

  ! Private module variable declarations:
  real(kind=dp) :: &     ! Values for the Kpar fit
       ialpha, ibeta, igamma, isigma, itheta, idl
  real(kind=dp), parameter :: &
       Q_o = 1368.0, & !1367.0? W/m^2  Solar constant
       albedo = 0.18   !KC 17% OCT.22 2004

contains

  subroutine init_irradiance()
    ! Initialize irradiance model
    implicit none

    ! Read Kpar parameter values from infile
    call read_irradiance_params()
  end subroutine init_irradiance


  subroutine read_irradiance_params()
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


  subroutine irradiance_sog(cf, day_time, day, In, Ipar, d, Qriver, &
       Pmicro, Pnano, Ppico, I_k)

    implicit none

    ! Arguments:
    real(kind=sp), intent(in) :: cf        ! cloud fraction
    real(kind=dp), intent(in) :: day_time  ! day-second
    integer, intent(in) :: day             ! year-day 
    type(grid_), intent(in) :: d           ! grid
    real(kind=dp), intent(in) :: Qriver    ! dominant river flow
    real(kind=dp), dimension(0:d%M), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano, &   ! Nano phytoplankton
         Ppico      ! Pico phytoplankton
    real(kind=dp), dimension(0:d%M), intent(out) :: In, Ipar  
    integer, intent(out) :: I_k
    
    ! Local variables:
    real(kind=dp) :: lat  ! Latitude of centre of model domain [rad]
    integer :: k, check, of                         
    real(kind=dp) :: declination, hour, cos_Z, day_length, hour_angle, &
         sunrise, sunset, Qso, a, b,KK
    real(kind=dp):: II    
    real(kind=dp), dimension(0:d%M) :: Ipar_i ! Ipar values on interfaces
    real(kind=dp) :: Ired, Iblue ! light in red and blue part of spectrum
    
 

!!!Define Okta Cloud Model!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! based on Dobson and Smith, table 5

    type :: cloudy
       real(kind=dp) :: A, B    !Regression coefficients
    end type cloudy

    type :: tenths
       type(cloudy), dimension(0:10) :: type
    end type tenths


    !!! SEA -- May 2010 : redid the cloud parameterization based on UBC
    !!! Solar data (/ocean/shared/SoG/met/solar/) fitting Q to cos_Z
    !!! (not Q/cos_Z as Kate did).  Page 61 in Lab Book.  (0) no
    !!! clouds, (1) 1/10 cloud fraction (10) 100% clouds.  Four sig
    !!! figs are what comes out of matlab but standard deviations are
    !!! 40W/m2 for low cloud fraction to 120 W/m2 for 6-9 cloud
    !!! fraction to 85 W/m2 for completely cloudy.
    type(tenths)::cloud
    cloud%type(0) = cloudy(0.6337,0.1959) 
    cloud%type(1) = cloudy(0.6149,0.2119)  
    cloud%type(2) = cloudy(0.5861,0.2400)
    cloud%type(3) = cloudy(0.5512,0.2859)
    cloud%type(4) = cloudy(0.5002,0.3192)
    cloud%type(5) = cloudy(0.4649,0.3356)
    cloud%type(6) = cloudy(0.4225,0.3339)
    cloud%type(7) = cloudy(0.3669,0.3490)
    cloud%type(8) = cloudy(0.2468,0.4427)
    cloud%type(9) = cloudy(0.1981,0.3116)
    cloud%type(10) = cloudy(0.0841, 0.2283)



    check = 0
    hour = (day_time/3600.0-12.0)*15.  !degrees
    declination = 23.45*PI/180.0*SIN((284.0+DBLE(day))/365.25*2.0*PI)  !radians

    ! Convert latitude of centre of model domain from degrees to radians
    lat = pi * latitude / 180.
    
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

    of = floor(cf+0.2)  ! "of" is integer type

    IF (day_time/3600.0 > sunrise .AND. day_time/3600.0 < sunset) THEN    
       II = Qso*(cloud%type(of)%A + cloud%type(of)%B*cos_Z)*cos_Z   
    ELSE
       II = 0.
    END IF
    ! so II is the incoming 

    Ipar_i = 0.        !PAR on interfaces
    Ipar = 0.        !PAR on grid points
    In = 0.          !total light on interfaces

    In(0) =  II        
    Ipar_i(0) = II*0.44  !44% of total light is PAR at surface (Jerlov)
    Iblue = 0.70*In(0)
    Ired = 0.30*In(0)


    ! Light is defined on interfaces
    ! Parameterized for a particular basin


    ! parameterization of Dec 2006 by S.Allen for the Strait of Georgia based on SOG data
         ! KK is evaluated on the grid points         
         ! N2chl is correction uM to mg/m3 chl
         ! limit of 2.5/m is set from Cruise 02-03 and 5 W/m2 seen at 2 m
         ! re-fit Jun 13, 2007 pg 113-114 lab-book and in fitlight
         ! 1.5d0 is correction for CTD fluoresence which seems high
 
    ! parameterization of Mar 2009  by M. Wolfe for Rivers Inlet based on RI 2007 hydrolab data
         ! limit of 2.5/m is set from Cruise 02-03 and 5 W/m2 seen at 2 m



    do k = 1, d%M    
         
       KK = ialpha +ibeta * (N2chl &
            * (Pmicro(k) + Pnano(k) + Ppico(k))) ** 0.665 &
            + (igamma * Qriver ** isigma + itheta) * exp(-d%d_g(k) / idl)
       KK = min(2.5d0, KK)
       Ipar_i(k) = Ipar_i(k-1) * exp(-d%i_space(k) * KK)
       
       if (Ipar_i(k) < 0.01d0 * Ipar_i(0) .and. check == 0) then            
          I_k = k                                          
          check = 1                 
       end if
       ! Total light for heat budget
       Iblue = Iblue * exp(-d%i_space(k) &
            * (0.8102d0 * KK + 1.1854d0))
       Ired = Ired * exp(-d%i_space(k) &
            * (0.8226d0 * KK - 0.0879d0))
       In(k) = Iblue + Ired

    end do

    Ipar(0) = Ipar_i(0)
    Ipar(1:d%M) = interp_g(Ipar_i)

   END SUBROUTINE irradiance_sog

end module irradiance
