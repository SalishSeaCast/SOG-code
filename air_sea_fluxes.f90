module air_sea_fluxes

  use precision_defs, only: dp, sp
  implicit none
  private
  public :: wind_stress
  
  real(kind=dp), parameter :: rho_atm = 1.25 ! kg/m3
  real(kind=dp) :: UU ! surface wind-speed (calc in wind-stress used in heat flux)
contains
  subroutine wind_stress (unow, vnow, rho, &
       w_u, w_v)
    ! subroutine to calculate the wind-stress

    implicit none
    ! Arguments :
    ! 35 degree and 325 degree wind components
    real(kind=dp), intent(in) :: unow, vnow, &
         ! surface water density (doesn't have to be current value)
         rho
    real(kind=dp), intent(out) :: w_u, w_v ! surface momentum flux components

    ! local variables
    real(kind=dp) :: C_D, & ! drag coefficient (Large and Pond)
         stress_u, stress_v  ! wind stress
    UU = SQRT(unow**2 + vnow**2)   

    if (UU /= 0.) then
       C_D = 1.0D-03 * (2.70 / UU + 0.142 + 0.0764 * UU)
       stress_u = unow / UU * C_D * rho_atm * UU**2
       stress_v = vnow / UU * C_D * rho_atm * UU**2
    else
       stress_u = 0.
       stress_v = 0.
    endif

    ! Momentum (eq'ns A2a & A2b)
    w_u = -stress_u / rho  ! w%u(0)
    w_v = -stress_v / rho  ! w%v(0)

  end subroutine wind_stress

  subroutine longwave_latent_sensible_heat(cf, atemp, humid, T_o, rho_o, Cp_o,&
       w_t)

    use unit_conversions, only: KtoC
    implicit none

    ! Arguments
    real(kind=sp),intent(in) :: cf, & ! cloud fraction in 0-1
         atemp, & ! air temperature
         humid ! humidity (%)
    real(kind=dp), intent(in) :: T_o, & ! sea surface water temperature
         rho_o, & ! surface water density
         Cp_o ! surface water specific heat
    real(kind=dp), intent(out) :: w_t ! heat flux at ocean surface
    ! local variables
    ! Longwave radiation
    real(kind=dp), parameter :: r = 0.03, &
         Ce = 9.37d-06, &
         sigma = 5.6697d-08, & ! Stefan Boltzmann constant
         epsilon_w = 0.96 ! surface emissivity in IR portion of the spectrum
    real(kind=dp) :: lw_in, & ! downward long wave radiation
         lw_out, & ! upward emission of long wave radiation 
         lw_net ! net longwave radiation

    ! Sensible heat flux
    real(kind=dp), parameter :: Cs = 1.3e-3, & ! sensible heat transfer coeff.
         Cp = 1003 ! specific heat of air
    real(kind=dp), h_sens ! sensible heat flux

    ! Latent heat flux 
    ! vapour pressure consts
    real(kind=dp), parameter :: a = 7.5, &
         b = 237.3, &
         c = 0.7858
    ! latent heat consts
    real (kind=dp), parameter :: CL = 1.3e-3, &
         LE = 2.453d06

    real(kind=dp) :: ea, & ! vapour pressure
         es, & ! saturated vapous pressure
         h_latent ! latent heat

    real(kind=dp) :: Q_t ! surface heat flux

    !downward radiation from atmosphere
    !there are several different ways to calculate this
    lw_in = (1 - r) * (1 +0.17 * cf**2) * Ce * atemp**2 &
         * sigma * atemp**4  

    !upward emission of radiation from earth's surface, stull page 48
    lw_out = -epsilon_w * sigma * T_o**4                      

    lw_net = lw_in + lw_out

    ! Sensible heat flux
    h_sens = Cs * rho_atm * Cp * UU * (atemp - T_o) 

    ! vapour pressure
    ea = (humid / 100.) * &
         exp(2.303 * ((a * KtoC(atemp) / ( KtoC(atemp) + b)) + c))  ! mb
    !saturated vapour pressure !*** same const as above??
    es = exp(2.3026 * ((a * KtoC(T_o) / (KtoC(T_o) + b)) + c))

    ! latent heat
    h_latent = (0.622 / 1013.) * CL * rho_atm * LE * UU * (ea-es)
    h_latent = max(h_latent,0.d0)

    Q_t = lw_net + h_sens + h_latent    ! W/m2

    ! Temperature (eq'n A2c)
    w_t = -Q_t / (rho_o * Cp_o) ! w%t(0)

  end subroutine longwave_latent_sensible_heat

end module air_sea_fluxes

