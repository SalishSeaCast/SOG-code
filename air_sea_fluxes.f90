module air_sea_fluxes
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to air-sea flux calculations in the SOG code.
  !
  ! Public Subroutines:
  !
  !   wind_stress -- Calculate the wind-stress.

  use precision_defs, only: dp, sp
  implicit none

  private
  public :: &
       ! Variables:
       ! Subroutines:
       wind_stress, longwave_latent_sensible_heat, gas_flux, solubility

  ! Parameter Value Declarations:
  !
  ! Private to module:
  real(kind=dp), parameter :: &
       rho_atm = 1.25d0,      &  ! Density of air [kg/m^3]
       pCO2_atm = 3.8d-4,     &  ! Partial pressure of CO2 [atm]
       pO2_atm = 0.20946d0       ! Partial pressure of O2 [atm]

  ! Variable Declarations:
  !
  ! Private to module:
  real(kind=dp) :: &
       UU  ! Surface wind-speed (calculated in wind-stress; used in heat flux)

contains

  subroutine wind_stress(unow, vnow, rho)
    ! Calculate the wind-stress.
    !
    ! This calculates the values of the wbar%u(0), and wbar%v(0)
    ! variables that are declared in the turbulence module.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variable Declarations:
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         unow, &  ! 35 degree wind component (cross-strait)
         vnow, &  ! 325 degree wind component (along-strait)
         rho      ! surface water density (doesn't have to be current
                  ! value)
    ! Local variables:
    real(kind=dp) :: &
         C_D,      &  ! drag coefficient (Large and Pond)
         stress_u, &  ! 35 degree wind-stress component (cross-strait)
         stress_v     ! 325 degree wind-stress component (along-strait)

    ! Calculate the surface wind-speed
    UU = sqrt(unow ** 2 + vnow ** 2)
    ! We need to handle dead, flat, calm as a special case to avoid a
    ! divide-by-zero in the drag coefficient formula.  Note that
    ! abs(x) < epsilon(x) is a real-number-robust test for x == 0.
    if (abs(UU) < epsilon(0.0d0)) then
       stress_u = 0.
       stress_v = 0.
    else
       C_D = 1.0d-3 * (2.70 / UU + 0.142 + 0.0764 * UU)
       stress_u = unow / UU * C_D * rho_atm * UU ** 2
       stress_v = vnow / UU * C_D * rho_atm * UU ** 2
    endif
    ! Calculate the surface turbulent kinematic flux momentum
    ! components (Large, et al, (1994), eq'ns A2a & A2b).
    wbar%u(0) = -stress_u / rho
    wbar%v(0) = -stress_v / rho
  end subroutine wind_stress


  subroutine longwave_latent_sensible_heat(cf, atemp, humid, T_o, rho_o, Cp_o)
    ! Calculate the surface turbulent heat flux.
    !
    ! This calculates the values of the wbar%t(0) variable that is
    ! declared in the turbulence module.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variable Declarations:
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays
    ! Functions:
    use unit_conversions, only: KtoC

    implicit none

    ! Arguments
    real(kind=sp),intent(in) :: &
         cf,    & ! Cloud fraction [0 to 1]
         atemp, & ! Air temperature [K]
         humid    ! humidity [%]
    real(kind=dp),intent(in) :: &
         T_o,   & ! Sea surface water temperature [K]
         rho_o, & ! Sea surface water density  [kg/m^3]
         Cp_o     ! Sea surface water specific heat [kJ/kg-K]
    ! Local parameter values and variables:
    !
    ! Longwave radiation
    real(kind=dp), parameter :: &
         r = 0.03,         &
         Ce = 9.37e-6,      &
         sigma = 5.6697e-8, &  ! Stefan Boltzmann constant
         epsilon_w = 0.96     ! surface emissivity in IR portion of
                                ! the spectrum
    real(kind=dp) :: &
         lw_in,  &  ! Downward long wave radiation
         lw_out, &  ! Upward emission of long wave radiation
         lw_net     ! Net longwave radiation
    !
    ! Sensible heat flux
    real(kind=dp), parameter :: &
         Cs = 1.3e-3, &  ! Sensible heat transfer coefficient
         Cp = 1003   ! Specific heat of air
    real(kind=dp) :: &
         h_sens  ! Sensible heat flux
    !
    ! Latent heat flux
    real(kind=dp), parameter :: &
         a = 7.5,      &  ! Vapour pressure consts
         b = 237.3,    &
         c = 0.7858,   &
         CL = 1.3e-3,    &  ! Latent heat consts
         LE = 2.453e6
    real(kind=dp) :: &
         ea,       &  ! Vapour pressure [mb]
         es,       &  ! Saturated vapous pressure [mb]
         h_latent, &  ! Latent heat
         Q_t          ! Surface heat flux [W/m^2]

    ! Downward radiation from atmosphere
    ! *** There are several different ways to calculate this
    lw_in = (1 - r) * (1 + 0.170 * cf**2) * Ce * atemp**2 &
         * sigma * atemp**4
    ! Upward emission of radiation from earth's surface, stull page 48
    lw_out = -epsilon_w * sigma * T_o**4
    lw_net = lw_in + lw_out
    ! Sensible heat flux
    h_sens = Cs * rho_atm * Cp * UU * (atemp - T_o)
    ! Vapour pressure in mb
    ea = (humid / 100) * &
         exp(2.303 * ((a * KtoC(atemp) / ( KtoC(atemp) + b)) + c))
    !Saturated vapour pressure
    ! *** Same const as above??
    es = exp(2.3026 * ((a * KtoC(T_o) / (KtoC(T_o) + b)) + c))
    ! Latent heat
    h_latent = (0.622 / 1013) * CL * rho_atm * LE * UU * (ea - es)
    h_latent = max(h_latent, 0.0d0)
    ! Surface heat flux in W/m^2
    Q_t = lw_net + h_sens + h_latent
    ! Surface turbulent heat flux (Large, et al (1994), eq'n A2c)
    wbar%t(0) = -Q_t / (rho_o * Cp_o)
  end subroutine longwave_latent_sensible_heat


  subroutine gas_flux(type, T, S, C_water, unow, vnow, gasflux, p_gas)
    ! Calculate the gas flux.
    !
    ! This calculates the value of the gasflux variable determined
    ! somewhere.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use io_unit_defs, only: stdout
    use unit_conversions, only: KtoC

    ! Variable Declarations:

    implicit none

    ! Arguments:
    character(len=*), intent(in):: &
         type         ! Either 'CO2' or 'Oxy'
    real(kind=dp), intent(in) :: &
         T,        &  ! Sea surface water temperature [K]
         S,        &  ! Sea surface practical salinity [PSU]
         C_water,  &  ! Sea surface gas concentration [uM]
         unow,     &  ! 35 degree wind component (cross-strait)
         vnow         ! 325 degree wind component (along-strait)
    real(kind=dp), intent(out) :: &
         gasflux,  &  ! Gas flux [umol m-2 s-1]
         p_gas        ! Partial pressure of gas [ppm]

    ! Local variables:
    real(kind=dp) :: &
         T_c,      &  ! Temperature in celcius
         sol,      &  ! Solubility [mol m-3 atm-1]
         sc,       &  ! Schmidt number (dimensionless)
         kps,      &  ! Transfer velocity [m/s]
         p_air,    &  ! Partial pressure in air [atm]
         ws10         ! Surface windspeed

    ! Temperature in celcius
    T_c = KtoC(T)

    ! Calculate the surface wind-speed
    ws10 = sqrt(unow ** 2 + vnow ** 2)

    ! Solubility [mol m-3 atm-1]
    call solubility(type, T, S, sol)

    ! Solubility and Schmidt Number
    if (type == 'Oxy') then
       ! from Wanninkhof 1992
       ! (T in deg C)
       sc = 1953.4d0 - 128.00d0 * T_c + 3.9918d0 * T_c**2 - 0.050091d0 * T_c**3

       ! Define partial pressure in air
       p_air = pO2_atm

    elseif (type == 'CO2') then
       ! from Wanninkhof et al. 1992
       ! (T in deg C)
       sc = 2073.1d0 - 125.62d0 * T_c + 3.6276d0 * T_c**2 - 0.043219d0 * T_c**3

       ! Define partial pressure in air
       p_air = pCO2_atm

    else
       write (stdout,*) 'gas_flux in air_sea_fluxes.f90: ', &
            'Unexpected gas type: ', type
       call exit(1)
    endif

    ! Transfer velocity, Nightingale et al. 2000
    kps = (0.22d0 * ws10**2 + 0.33d0 * ws10) * (sc/600.0d0)**(-0.5d0)
    kps = kps * 2.778d-6   ! m/s

    ! Partial Pressure
    p_gas = 1.0d3 * (C_water/sol)

    ! Gas flux
    ! (m/s * mol/m^3atm * (umol/mol * m^3/L) * atm = umol m-2 s-1)
    gasflux = kps * sol * 1.0d3 * (1.0d-6 * p_gas - p_air)

  end subroutine gas_flux


  subroutine solubility(type, T, S, sol)
    ! Calculate the solubility of a gas type, 'type'
    !
    ! This calculates the value of the 'sol' variable needed in
    ! the gas_flux subroutine and the solve_gas_flux subroutine in
    ! chemistry_fluxes

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use io_unit_defs, only: stdout

    ! Variable Declarations:

    implicit none

    ! Arguments:
    character(len=*), intent(in):: &
         type         ! Either 'CO2' or 'Oxy'
    real(kind=dp), intent(in) :: &
         T,        &  ! Sea surface water temperature [K]
         S            ! Sea surface practical salinity [PSU]
    real(kind=dp), intent(out) :: &
         sol          ! Solubility [mol m-3 atm-1]

    if (type == 'Oxy') then
       ! from Wanninkhof 1992
       ! Divide by molar volume of O2 (22.3914 from SeaBird Electronics)
       ! (multiply by rho/rho = 1 for mol m-3 atm-1)
       sol = exp(85.8079d0/(T/100.0d0) - 58.3877d0 + 23.8439d0 * &
            log(T/100.0d0) + S * (-0.034892d0 + 0.015568d0 * (T/100.0d0) - &
            0.0019387d0 * (T/100.0d0)**2))/22.3914d-3
    elseif (type == 'CO2') then
       ! from Weiss 1974 (M/atm)
       ! (multiply by 1.0e3 for mol m-3 atm-1)
       sol = exp(90.5069d0/(T/100.0d0) - 58.0931d0 + 22.2940d0 * &
            log(T/100.0d0) + S * (0.027766d0 - 0.025888d0 * (T/100.0d0) + &
            0.0050578d0 * (T/100.0d0)**2)) * 1.0d3
    else
       write (stdout,*) 'solubility in air_sea_fluxes.f90: ', &
            'Unexpected gas type: ', type
       call exit(1)
    endif

  end subroutine solubility

end module air_sea_fluxes
