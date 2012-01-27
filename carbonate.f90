module carbonate
  ! Variable & parameter value declarations, and subroutines related
  ! to carbonate system calculations in the SOG code
  !
  ! Public Variables:
  !
  !   DIC      -- Dissolved inorganic carbon in surface
  !
  !   CO2_star -- Total CO2 (CO2*) in surface water
  !
  !   Alk      -- Total Alkalinity based on salinity fit
  !
  ! Public Subroutines:
  !
  !   calculate_co2 -- Calculate between CO2* and DIC at T, S and 1 atm.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Carbonate system constants
       K0,          &  ! CO2 solubility
       K1,          &  ! First equilibrium constant [H][HCO3]/[H2CO3]
       K2,          &  ! Second equilibrium constant [H][CO3]/[HCO3]
       Kb,          &  ! Borate equilibrium constant [H][BO2]/[HBO2]
       Kw,          &  ! Dissociation of water [H][OH]
       BO3_total,   &  ! Total borate
       ! Carbonate variables
       Alk,         &  ! Total alkalinity [uM]
       ! Subroutines
       calculate_co2

  ! Variable declarations:
  real(kind=dp) :: &
       ! Carbonate system constants
       K0,          &  ! CO2 solubility
       K1,          &  ! First equilibrium constant [H][HCO3]/[H2CO3]
       K2,          &  ! Second equilibrium constant [H][CO3]/[HCO3]
       Kb,          &  ! Borate equilibrium constant [H][BO2]/[HBO2]
       Kw,          &  ! Dissociation of water [H][OH]
       BO3_total,   &  ! Total borate
       ! Carbonate variables
       Alk             ! Total alkalinity [uM]

contains

  subroutine set_constants(T, S)
    ! Calculate the carbonate system constants K0, K1, K2, Kb, Kw,
    ! and total borate from T and S
    !
    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use io_unit_defs, only: stdout

    ! Variable Declarations:
    implicit none

    ! Arguments
    real(kind=dp), intent(in) :: &
         T,       &  ! Temperature [K]
         S           ! Practical salinity [PSU]

    ! Calculate solubility and equilibrium constants

    ! From Weiss 1974
    K0 = exp(93.4517d0/(T/100.0d0) - 60.2409d0 + 23.3585d0 * log(T/100.0d0) + &
         S * (0.023517d0 - 0.023656d0 * (T/100.0d0) + 0.0047036d0 * &
         (T/100.0d0)**2))

    ! Millero p.664 (1995) using Mehrbach et al. data on seawater scale 
    K1 = 10.0d0**(-1.0d0*((3670.7d0/T) - 62.008d0 + 9.7944d0 * log(T) - & 
         0.0118d0 * S + 0.000116d0 * S**2))

    K2 = 10.0d0**(-1.0d0*((1394.7d0/T) + 4.777d0 - 0.0184d0 * S + &
         0.000118d0 * S**2))

    ! Millero p.669 (1995) using data from Dickson (1990)
    Kb = exp((-8966.90d0 - 2890.53d0 * sqrt(S) - 77.942d0 * S + &
         1.728d0 * S**1.5d0 - 0.0996d0 * S**2)/T + &
         148.0248d0 + 137.1942d0 * sqrt(S) + 1.62142d0 * S + &
         (-24.4344d0 - 25.085d0 * sqrt(S) - 0.2474d0 * S) * log(T) + &
         0.053105d0 * sqrt(S) * T) 

    ! Millero p.670 (1995) using composite data
    Kw = exp(-13847.26d0/T + 148.9652d0 - 23.6521d0 * log(T) + &
         (118.67d0/T - 5.977d0 + 1.0495d0 * log(T)) * sqrt(S) - &
         0.01615d0 * S)

    ! Uppstrom (1974)
    BO3_total = 0.000232d0 * (S/19.53061205d0)

  end subroutine set_constants


  subroutine calculate_co2(T, S, rho, DIC, CO2_star)
    ! Calculate CO2* from DIC at T, S and 1 atm. Total alkalinity
    ! is required for this calculation, but for our purposes will
    ! obtained using a linear fit to salinity.
    ! (2010 SoG cruise data)
    !
    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp

    ! Variable Declarations:
    implicit none
    
    ! Arguments
    real(kind=dp), intent(in) :: &
         T,          &    ! Temperature [K]
         S,          &    ! Practical Salinity [PSU]
         rho,        &    ! Surface seawater density [kg/m^3]
         DIC              ! Dissolved inorganic carbon [uM]
    real(kind=dp), intent(out) :: &
         CO2_star         ! CO2* [uM]

    ! Local variables
    ! Iteration parameters
    real(kind=dp) :: &
         x1,         &    ! Lower [H+] iteration boundary
         x2,         &    ! Upper [H+] iteration boundary
         x_guess,    &    ! Iteration [H+] initialization value 
         x_acc            ! Iteration [H+] difference limit
    ! Carbonate system properties
    real(kind=dp) :: &
         H_total,    &    ! Total [H+]
         DIC_molal        ! DIC converted to mol/kg

    ! Total alkalinity in terms of salinity (mol/kg)
    ! Fit from 2010 SOG cruises
    Alk = (47.3210d0 * S + 653.8105d0) * 1.0d-6

    ! DIC from uM to mol/kg
    DIC_molal = DIC * (1.0d-3/rho)

    ! Set equilibrium constants
    call set_constants(T, S)

    ! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
    !
    ! The solution converges to err of x_acc. The solution must be within
    ! the range x1 to x2.
    !
    ! If DIC and TA are known then either a root finding or iterative method
    ! must be used to calculate H_total. In this case we use the Newton-Raphson
    ! "safe" method taken from "Numerical Recipes" (function "rtsafe.f" with
    ! error trapping removed).
    !
    ! As currently set, this procedure iterates about 12 times. The x1 and x2
    ! values set below will accomodate ANY oceanographic values. If an initial
    ! guess of the pH is known, then the number of iterations can be reduced to
    ! about 5 by narrowing the gap between x1 and x2. It is recommended that
    ! the first few time steps be run with x1 and x2 set as below. After that,
    ! set x1 and x2 to the previous value of the pH +/- ~0.5. The current
    ! setting of x_acc will result in co2_star accurate to 3 significant figures
    ! (xx.y). Making x_acc bigger will result in faster convergence also, but
    ! this is not recommended (x_acc of 1e-9 drops precision to 2 significant
    ! figures)

    x1 = 1.0d-7
    x2 = 1.0d-9
    x_guess = 1.0d-8
    x_acc = 1.0d-10

    call drt_safe(DIC_molal, x1, x2, x_guess, x_acc, H_total)

    ! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
    ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    CO2_star = DIC_molal * H_total**2/(H_total**2 + K1 * H_total + K1 * K2)

    ! Convert units of output arguments
    ! Note: CO2_star and dCO2_star are calculated in mol/kg within
    ! this routine, thus convert now from mol/kg -> uM
    CO2_star = CO2_star * rho * 1.0d3

  end subroutine calculate_co2


  subroutine ta_iter(DIC, H, F, dFdH)
    ! Expresses TA as a function of DIC, H_total and constants,
    ! and calculates the derivative dTA([H+])/d[H+], which is used
    ! in the iterative solution for H_total in the drt_safe subroutine.
    !
    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp

    ! Variable Declarations:
    ! Arguments
    real(kind=dp), intent(in) :: &
         DIC,        &  ! Surface dissolved inorganic carbon [uM]
         H              ! Total hydrogen ion concentration, [H+]
    real(kind=dp), intent(out) :: &
         F,          &  ! Alkalinity difference function
         dFdH           ! Derivative of F with respect to [H+]

    ! Local variables
    real(kind=dp) :: &
         Beta           ! This is the factor [DIC]*K1*K2/[CO3]

    ! Compute Beta = [DIC]*K1*K2/[CO3]
    Beta = H**2 + K1 * H + K1 * K2

    ! Compute the difference, F, between the alkalinity calculated
    ! from a given [DIC], [H+], and [BT], and the given alkalinity
    F = K1 * H * DIC/Beta + &
         2.0d0 * DIC * K1 * K2/Beta + & 
         BO3_total/(1.0d0 + H/Kb) + Kw/H - H - Alk
       
    ! Compute dTAdH = dTA([H+])/d[H+]
    dFdH = ((K1 * DIC * Beta) - &
         K1 * H * DIC * (2.0d0 * H + K1))/Beta**2 - &
         2.0d0 * DIC * K1 * K2 * (2.0d0 * H + K1)/Beta**2 - &
         BO3_total/(Kb * (1.0d0 + H/Kb)**2) - Kw/H**2 - 1.0d0

  end subroutine ta_iter


  subroutine drt_safe(DIC, x1, x2, x_guess, x_acc, x)
    ! Safe iteration to find the zero of a function
    !
    ! Returns the correct zero of the TA function using the ta_iter
    ! subroutine by minimizing the difference between x1, x2, and x_guess
    ! specified in the calculate_pco2 subroutine
    !
    ! File taken from Numerical Recipes. Modified  R.M.Key 4/94
    ! Transformed into FORTRAN90 (only minimal changes) by C.Voelker, 4/05
    ! Reformatted for SOG deep estuary model by B.Moore-Maley, 8/2011
    !
    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp

    ! Variable Declarations:
    implicit none

    ! Arguments
    real(kind=dp), intent(in) :: &
         DIC,        &    ! Dissolved inorganic carbon [uM]
         x1,         &    ! First zero boundary
         x2,         &    ! Second zero boundary
         x_guess,    &    ! First zero guess
         x_acc            ! Zero difference limit
    real(kind=dp), intent(out) :: &
         x                ! Function zero
    ! Local Variables
    real(kind=dp) :: &
         F_low,      &    ! F(x_low)
         F_high,     &    ! F(x_high)
         F,          &    ! F(x)
         dFdx,       &    ! dF(x)/dx
         x_low,      &    ! Lower zero boundary
         x_high,     &    ! Upper zero boundary
         dx,         &    ! Zero difference (i)
         dx_old,     &    ! Zero difference (i-1)
         swap,       &    ! Temporary F_low, F_high exchange variable
         temp             ! Temporary zero comparison variable
    integer :: &
         n_iter,     &    ! Interation number
         max_iter         ! Maximum number of iterations

    ! Set maximum number of iterations to 100
    max_iter = 100

    ! Generate upper and lower boundaries
    call ta_iter(DIC, x1, F_low, dFdx)
    call ta_iter(DIC, x2, F_high, dFdx)
    if(F_low .lt. 0.0d0) then
       x_low = x1
       x_high = x2
    else
       x_high = x1
       x_low = x2
       swap = F_low
       F_low = F_high 
       F_high = swap
    end if

    ! Initialize zeros for iteration loop
    x = x_guess
    dx_old = abs(x2 - x1)
    dx = dx_old
    call ta_iter(DIC, x, F, dFdx)
    ! Begin iteration loop
    do n_iter = 1, max_iter
       if(((x - x_high) * dFdx - F) * ((x - x_low) * &
            dFdx - F) .ge. 0.0d0 .or. abs(2.0d0 * F) .gt. &
            abs(dx_old * dFdx)) then
          dx_old = dx
          dx = 0.5d0 * (x_high - x_low)
          x = x_low + dx
          if(x_low .eq. x) return
       else
          dx_old = dx
          dx = F/dFdx
          temp = x
          x = x - dx
          if(temp .eq. x) return
       end if
       if(abs(dx) .lt. x_acc) return
       call ta_iter(DIC, x, F, dFdx)
       if(F .lt. 0.0d0) then
          x_low = x
          F_low = F
       else
          x_high = x
          F_high = F
       end if
    end do

  end subroutine drt_safe

end module carbonate

