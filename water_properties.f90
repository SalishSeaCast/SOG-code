module water_properties
  ! Type definition, variable declarations, and subroutine to
  ! calculate water properties: density (rho), thermal expansion
  ! coefficient (alpha), salinity expansion coefficient (beta), and
  ! specific heat capacity (Cp).
  !
  ! Public type:
  !
  ! water_property -- Density, thermal expansion and salinity contraction
  !                   coefficients, and specific heat capacity profiles.
  !
  ! Public Variables:
  ! 
  ! Cp -- Specific heat capacity [J/kg.K]
  !
  ! rho -- Density [kg/m^3]
  ! 
  ! alpha -- Thermal expansion coefficient [K^-1]
  !
  ! beta -- Saline contraction coefficient [K^-1]
  ! 
  ! Public Subroutines:
  !
  ! calc_rho_alpha_beta_Cp_profiles -- Calculate profiles of rho, alpha,
  !                                    beta, and Cp at grid point depths.
  !
  ! alloc_water_props -- Allocate memory for water property profiles.
  !
  ! dalloc_water_props -- Deallocate memory from water property
  !                       profiles.

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Type:
       water_property, &
       ! Variables:
       rho,   &  ! Density [kg/m^3]
       alpha, &  ! Thermal expansion coefficient [K^-1]
       beta,  &  ! Saline contraction coefficient [K^-1]
       Cp,    &  ! Specific heat capacity [J/kg.K]
       ! Subroutines:
       calc_rho_alpha_beta_Cp_profiles, &
       alloc_water_props,  &
       dalloc_water_props, &
       oxygen_saturation

  ! Public type definition:
  !
  ! Density, thermal expansion and saline contraction coefficients,
  ! and specific heat capacity profiles
  type :: water_property
     real(kind=dp), dimension(:), pointer :: &
          g,      &  ! Value at grid point depth
          grad_g, &  ! Rate of change (d/dz) at grid point depth
          i,      &  ! Value at grid cell interface depth
          grad_i     ! Rate of change (d/dz) at grid cell interface depth
  end type water_property

  ! Public variable declarations:
  !
  ! Density, thermal expansion and saline contraction coefficients,
  ! and specific heat capacity profiles
  type(water_property) :: &
       rho,   &  ! Density [kg/m^3]
       alpha, &  ! Thermal expansion coefficient [K^-1]
       beta,  &  ! Saline contraction coefficient [K^-1]
       Cp        ! Specific heat capacity [J/kg.K]

contains

  subroutine alloc_water_props(M)
    ! Allocate memory for water property profiles.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in)               :: M   ! Number of grid points
    ! Local variables:
    integer              :: allocstat  ! Allocation return status
    character(len=80)    :: msg        ! Allocation failure message prefix

    msg = "Water density profile arrays"
    allocate(rho%g(0:M+1), rho%i(0:M), rho%grad_g(1:M), rho%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Water thermal expansion coefficient profile arrays"
    allocate(alpha%g(0:M+1), alpha%i(0:M), alpha%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Water saline contraction coefficient profile arrays"
    allocate(beta%g(0:M+1), beta%i(0:M), beta%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Water heat capacity profile (Cp) arrays"
    allocate(Cp%g(0:M+1), Cp%i(0:M), stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_water_props


  subroutine dalloc_water_props
    ! Deallocate memory for water property profiles.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer              :: dallocstat  ! Allocation return status
    character(len=80)    :: msg         ! Allocation failure message prefix

    msg = "Water density profile arrays"
    deallocate(rho%g, rho%i, rho%grad_g, rho%grad_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Water thermal expansion coefficient profile arrays"
    deallocate(alpha%g, alpha%i, alpha%grad_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Water saline contraction coefficient profile arrays"
    deallocate(beta%g, beta%i, beta%grad_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Water specific heat capacity profile (Cp) arrays"
    deallocate(Cp%g, Cp%i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_water_props


  subroutine calc_rho_alpha_beta_Cp_profiles(T, S)
    ! Calculate the profiles of density, thermal expansion
    ! coefficient, salinity expansion coefficient, and specific heat
    ! capacity through the water column at the grid point depths based
    ! on the supplied water temperature and salinity profiles.
    !
    ! This calculates the values of the grid point depth arrays (*%g) of
    ! the 4 water properties.
    !
    ! Based on Rich Pawlowicz's Matlab functions in his OCEANS
    ! toolbox (http://www.eos.ubc.ca/~rich/#OCEANS) with the pressure
    ! dependency removed.  rho, alpha, and beta calculations are based
    ! on swstate.m.  Cp calculation is based on swcp.m.  In all cases the
    ! pressure dependency has been removed.  Variable names are, as
    ! far as possible, the same as those in the .m files.
    !
    ! The calculations have been combined into 1 subroutine because
    ! they all require the evaluation of the square root of the
    ! salinity profile, and this way we can minimize the number of
    ! times that CPU-cycle intensive operation is done.
    !
    ! rho and Cp results have been verified driectly against
    ! http://fermi.jhuapl.edu/denscalc.html, and alpha and beta
    ! results have been verified against the same calculator by
    ! calculating slopes over 0.1 K and 0.1 PSU intervals.
    !
    ! Check values: 
    !   rho(283.15 K, 24 PSU) = 1018.379289 kg/m^3
    ! alpha(283.15 K, 24 PSU) = 0.000144174 K^-1
    !  beta(283.15 K, 24 PSU) = 0.000763851 K^-1
    !    Cp(283.15 K, 24 PSU) = 4047.5 J/kg-K
    use unit_conversions, only: KtoC
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(in) :: &
         T, &    ! Temperature [K]
         S       ! Salinity [-]
    ! Local variables:
    real(kind=dp), dimension(0:size(T) - 1) :: &
         TC,         &  ! Temperature profile in degrees Celcius
         sqrtS,      &  ! Square root of salinity profile
         A, B, C,    &  ! Coefficients for specific heat calculation
         R1, R2, R3, R4S,              &  ! Coefficients for density calc
         SIG, SIGMA, V350P, SVA, DVAN, &  ! Temporary vars for density calc
         RHO1, RHOS, V0S, drho_dS,     &  ! Temporary variables for beta calc
         RHOT, V0T, V, drho_dT            ! Temporary variables for alpha calc
    ! Parameters for density calculation
    real(kind=dp), parameter :: &
         R3500 = 1028.1063, &
         R4 = 4.8314d-4,    &
         DR350 = 28.106331

    ! Convert the temperature to Celsius
    TC = KtoC(T)
    ! Calculate the square root of the salinities
    sqrtS = sqrt(S)

    ! Calculate the density profile at the grid point depths
    ! Pure water density at atmospheric pressure
    ! (Bigg P.H., (1967) Br. J. Applied Physics 8 pp 521-537)
    R1 = ((((6.536332d-9 * TC - 1.120083d-6) * TC + 1.001685d-4) &
            * TC - 9.095290d-3) * TC + 6.793952d-2) * TC - 28.263737
    ! Seawater density at atmospheric pressure
    !   Coefficients of salinity
    R2 = (((5.3875d-9 * TC - 8.2467d-7) * TC + 7.6438d-5) &
           * TC - 4.0899d-3) * TC + 8.24493d-1
    R3 = (-1.6546d-6 * TC + 1.0227d-4) * TC - 5.72466d-3
    ! International one-atmosphere equation of state of seawater
    SIG = (R4 * S + R3 * sqrtS + R2) * S + R1
    ! Specific volume at atmospheric pressure
    V350P = 1.0 / R3500
    SVA = -SIG * V350P / (R3500 + SIG)
    ! Density at atmospheric pressure
    SIGMA = SIG + DR350
    ! Derivative drho/dS at atmospheric pressure
    R4S = 9.6628E-4
    RHO1 = 1000.0 + SIGMA
    RHOS = R4S * S + 1.5 * R3 * sqrtS + R2
    V0S = -RHOS / (RHO1 ** 2)
    ! Derivative drho/dT at atmospheric pressure
    R1 =(((3.268166d-8 * TC - 4.480332d-6) * TC + 3.005055d-4) &
         * TC - 1.819058d-2) * TC + 6.793952d-2
    R2 = ((2.155d-8 * TC - 2.47401d-6) * TC + 1.52876d-4) * TC - 4.0899d-3
    R3 = -3.3092d-6 * TC + 1.0227d-4
    RHO1 = 1000.0 + SIGMA
    RHOT = (R3 * sqrtS + R2) * S + R1
    V0T = -RHOT / (RHO1 ** 2)
    ! Density anomoly at atmospheric pressure
    DVAN = SVA / (V350P * (V350P + SVA))
    SIGMA = DR350 - DVAN
    ! Density at atmospheric pressure
    rho%g = 1000. + SIGMA
    ! Specific volume at atmospheric pressure
    V = 1. / rho%g
    ! Thermal expansion coefficient at atmostpheric pressure
    drho_dT = -V0T / (V ** 2)
    alpha%g = -(1. / rho%g) * drho_dT
    ! Saline contraction coefficient at atmospheric pressure
    drho_dS = -V0S / (V ** 2)
    beta%g = (1. /rho%g) * drho_dS

    ! Calculate the heat capacity profile at the grid point depths
    ! Specific heat Cp0 for P=0 (Millero, et al; UNESCO 1981)
    A = (-1.38385d-3 * TC + 0.1072763) * TC - 7.643575
    B = (5.148d-5 * TC - 4.07718d-3) * TC + 0.1770383
    C = (((2.093236d-5 * TC - 2.654387d-3) * TC + 0.1412855) &
         * TC - 3.720283) * TC + 4217.4
    Cp%g = (A + B * sqrtS) * S + C
  end subroutine calc_rho_alpha_beta_Cp_profiles

  subroutine oxygen_saturation(sal,temp,Osat)
  ! subroutine to calculate oxygen saturation values in [umol/L]
  ! sal is salinity, temp is temperature in K
  ! Conversion between mg and ml:  OST[ml/l] = OST[mg/l] * 0.7 [ml/mg]
  ! This is a modified version of the function oxsat by Christian Mertens, IfM Kiel
  real(kind=dp), intent(in) :: sal, temp
  real(kind=dp), intent(out) :: Osat
  ! local variable
  real(kind=dp) :: T100
  T100 = 0.01d0*temp
  Osat = -173.4292d0 + 249.6339d0/T100 + 143.3483d0*log(T100) - 21.8492d0*T100 &
        + sal*(-0.033096d0 + (0.014259d0 - 0.0017d0*T100)*T100)
  Osat = exp(Osat)            ! oxygen saturation in [ml/L]
  Osat = Osat*1000.0d0/22.4d0 ! [umol/L] or [uM]
  end subroutine oxygen_saturation

end module water_properties
