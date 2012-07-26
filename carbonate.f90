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
       ! Carbonate system constant arrays
       K0,          &  ! CO2 solubility
       K1,          &  ! First equilibrium constant [H][HCO3]/[H2CO3]
       K2,          &  ! Second equilibrium constant [H][CO3]/[HCO3]
       Kb,          &  ! Borate equilibrium constant [H][BO2]/[HBO2]
       Kw,          &  ! Dissociation of water [H][OH]
       BO3_total,   &  ! Total borate
       ! Subroutines
       DIC_to_carbonate, pH_to_carbonate

  ! Variable declarations:
  real(kind=dp) :: &
       ! Carbonate system constants
       K0,          &  ! CO2 solubility
       K1,          &  ! First equilibrium constant [H][HCO3]/[H2CO3]
       K2,          &  ! Second equilibrium constant [H][CO3]/[HCO3]
       Kb,          &  ! Borate equilibrium constant [H][BO2]/[HBO2]
       Kw,          &  ! Dissociation of water [H][OH]
       BO3_total       ! Total borate

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
    real(kind=dp) :: &
         T,       &  ! Temperature [K]
         S           ! Practical salinity [PSU]

    ! Calculate solubility and equilibrium constants

    ! From Weiss 1974
    K0 = exp(93.4517d0/(T/100.0d0) - 60.2409d0 + 23.3585d0 * log(T/100.0d0) &
         + S * (0.023517d0 - 0.023656d0 * (T/100.0d0) + 0.0047036d0 * &
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


  subroutine DIC_to_carbonate(T, S, rho, Alk, DIC, CO2_star)
    ! Calculate pCO2 from DIC at T, S and 1 atm. Total alkalinity
    ! is required for this calculation, but for our purposes will
    ! obtained using a linear fit to salinity.
    ! (2010 SoG cruise data)
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
         T,          &    ! Temperature [K]
         S,          &    ! Practical Salinity [PSU]
         rho,        &    ! Seawater density [kg/m^3]
         Alk,        &    ! Total alkalinity [ueq/L]
         DIC              ! DIC [uM]
    real(kind=dp), intent(out) :: &
         CO2_star         ! CO2_star [uM]
    ! Local variables
    ! Carbonate system properties
    real(kind=dp) :: &
         Alk_molal,  &    ! Total alkalinity [ueq/kg]
         DIC_molal,  &    ! DIC [mol/kg]
         H_total          ! Total [H+]
    ! Iteration parameters
    real(kind=dp) :: &
         x1,         &    ! Lower [H+] iteration boundary
         x2,         &    ! Upper [H+] iteration boundary
         x_guess,    &    ! Iteration [H+] initialization value 
         x_acc            ! Iteration [H+] difference limit

    ! Set equilibrium constants
    call set_constants(T, S)

    ! Convert from uM to mol/kg
    Alk_molal = Alk * (1.0d-3/rho)
    DIC_molal = DIC * (1.0d-3/rho) ! DIC [mol/kg]

    ! Calculate [H+] total when DIC and TA are known at T, S and 1 atm.
    !
    ! The solution converges to err of x_acc. The solution must be within
    ! the range x1 to x2.
    !
    ! If DIC and TA are known then either a root finding or iterative
    ! method must be used to calculate H_total. In this case we use the
    ! Newton-Raphson "safe" method taken from "Numerical Recipes"
    ! (function "rtsafe.f" with error trapping removed).
    !
    ! As currently set, this procedure iterates about 12 times. The x1
    ! and x2 values set below will accomodate ANY oceanographic values.
    ! If an initial guess of the pH is known, then the number of iterations
    ! can be reduced to about 5 by narrowing the gap between x1 and x2.
    ! It is recommended that the first few time steps be run with x1 and x2
    ! set as below. After that, set x1 and x2 to the previous value of the
    ! pH +/- ~0.5. The current setting of x_acc will result in co2_star
    ! accurate to 3 significant figures (xx.y). Making x_acc bigger will
    ! result in faster convergence also, but this is not recommended
    ! (x_acc of 1e-9 drops precision to 2 significant figures)

    x1 = 1.0d-7
    x2 = 1.0d-9
    x_guess = 1.0d-8
    x_acc = 1.0d-10

    call drt_safe('DIC', Alk_molal, DIC_molal, x1, x2, x_guess, x_acc, H_total)

    ! Calculate [CO2*] as defined in DOE Methods Handbook 1994 Ver.2, 
    ! ORNL/CDIAC-74, Dickson and Goyet, eds. (Ch 2 p 10, Eq A.49)
    CO2_star = DIC_molal * H_total**2/(H_total**2 + K1 * H_total + K1 * K2)

    ! Convert from mol/kg to uM
    CO2_star = CO2_star * rho * 1.0d3

  end subroutine DIC_to_carbonate


  subroutine pH_to_carbonate(T, S, rho, Alk, pH, DIC)
    ! Calculate DIC from total alkalinity and pH at T, S and 1 atm.
    !
    ! Elements from other modules:
    ! Type Definitions:
    use precision_defs, only: dp
    use io_unit_defs, only: stdout

    ! Variable Declarations:
    implicit none
    
    ! Arguments
    real(kind=dp), intent(in) :: &
         T,          &    ! Temperature [K]
         S,          &    ! Practical Salinity [PSU]
         rho,        &    ! Seawater density [kg/m^3]
         Alk,        &    ! Total alkalinity [ueq/L]
         pH               ! pH
    real(kind=dp), intent(out) :: &
         DIC              ! DIC [uM]
    ! Local variables
    ! Carbonate system properties
    real(kind=dp) :: &
         Alk_molal,  &    ! Total alkalinity [ueq/kg]
         H_total          ! Total [H+] [mol/kg]

    ! Set equilibrium constants
    call set_constants(T, S)

    ! Convert from uM to mol/kg
    Alk_molal = Alk * (1.0d-3/rho)

    ! Convert pH to [H+]
    H_total = 10**(-pH)

    ! Calculate DIC as defined in CDIAC Best Practices 2007, PICES Special
    ! Publication 3, Dickson et al. eds. (Ch 2 pp 7-9, Eq 17, 53-59)
    ! See BMM Lab book pg 49
    DIC = (1.0d0 / (H_total + 2.0d0 * K2)) &
        * (Alk_molal - BO3_total/(1.0d0 + H_total/Kb) - Kw/H_total + H_total) &
        * (H_total**2 / K1 + H_total + K2)

    ! Convert units of output arguments
    ! Note: CO2_star, dCO2_star, and DIC are calculated in mol/kg within
    ! this routine, thus convert now from mol/kg -> uM
    DIC = DIC * rho * 1.0d3

  end subroutine pH_to_carbonate


  subroutine ta_iter(qty, Alk, par, H, F, dFdH)
    ! Expresses TA as a function of parameter, PAR, indicated by PARTYPE,
    ! H_total and constants, and calculates the derivative dTA([H+])/d[H+],
    ! which is used in the iterative solution for H_total in the
    ! drt_safe subroutine.
    !
    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use io_unit_defs, only: stdout

    ! Variable Declarations:
    ! Arguments
    character(len=*), intent(in) :: &
         qty            ! 'DIC' or 'pCO2'
    real(kind=dp), intent(in) :: &
         Alk,        &  ! Total alkalinity [ueq/L]
         par,        &  ! Surface DIC [uM] or pCO2 (uatm)
         H              ! Total hydrogen ion concentration, [H+]
    real(kind=dp), intent(out) :: &
         F,          &  ! Alkalinity difference function
         dFdH           ! Derivative of F with respect to [H+]

    ! Local variables
    real(kind=dp) :: &
         H2,         &  ! [H+] squared
         Beta,       &  ! This is the factor [DIC]*K1*K2/[CO3]
         Beta2          ! Beta-squared

    ! [H+] squared
    H2 = H**2

    if (qty == 'DIC') then

       ! Compute Beta = [DIC]*K1*K2/[CO3]
       Beta = H2 + K1 * H + K1 * K2
       Beta2 = Beta**2

       ! The root equation, F, for Alkalinity in terms of [H+], DIC, and Alk
       F = K1 * H * par/Beta + 2.0d0 * par * K1 * K2/Beta + &
            BO3_total/(1.0d0 + H/Kb) + Kw/H - H - Alk

       ! Compute dTAdH = dTA([H+])/d[H+]
       dFdH = ((K1 * par * Beta) - K1 * H * par * (2.0d0 * H + K1))/Beta2 - &
            2.0d0 * par * K1 * K2 * (2.0d0 * H + K1)/Beta2 - &
            BO3_total/(Kb * (1.0d0 + H/Kb)**2) - Kw/H2 - 1.0d0

    elseif (qty == 'pCO2') then

       ! See BMM Lab book pg 49
       ! The root equation, F, for Alkalinity in terms of [H+], pCO2, and Alk
       F = (K1 * par) / H + (2.0d0 * K1 * K2 * par) / H2  &
            + BO3_total/(1.0d0 + H/Kb) + Kw/H - H - Alk

       ! Compute dTAdH = dTA([H+])/d[H+]
       dFdH = -1.0d0 * (K1 * par) / H2 - (4.0d0 * K1 * K2 * par) / (H2 * H) - &
            BO3_total/(Kb * (1.0d0 + H/Kb)**2) - Kw/H2 - 1.0d0

    else
       write (stdout, *) 'parameter type in carbonate.f90:', &
            'Unexpected value: ', qty
       call exit(1)
    endif

  end subroutine ta_iter


  subroutine drt_safe(qty, Alk, par, x1, x2, x_guess, x_acc, x)
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
    character(len=*), intent(in) :: &
         qty              ! 'DIC' or 'pCO2'
    real(kind=dp), intent(in) :: &
         Alk,        &    ! Total alkalinity [ueq/L]
         par,        &    ! DIC [uM] or pCO2 [uatm]
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
    call ta_iter(qty, Alk, par, x1, F_low, dFdx)
    call ta_iter(qty, Alk, par, x2, F_high, dFdx)
    if(F_low .lt. 0.0d0) then
       x_low = x1
       x_high = x2
    else
       x_high = x1
       x_low = x2
       swap = F_low
       F_low = F_high 
       F_high = swap
    endif

    ! Initialize zeros for iteration loop
    x = x_guess
    dx_old = abs(x2 - x1)
    dx = dx_old
    call ta_iter(qty, Alk, par, x, F, dFdx)
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
       endif
       if(abs(dx) .lt. x_acc) return
       call ta_iter(qty, Alk, par, x, F, dFdx)
       if(F .lt. 0.0d0) then
          x_low = x
          F_low = F
       else
          x_high = x
          F_high = F
       endif
    end do

  end subroutine drt_safe

  
  subroutine ca_solubility(S, T_K, P, rho, DIC, pH, Omega_ca, Omega_ar)
    ! Taken from CO2SYS subfunction CaSolubility
    ! ***********************************************************************
    ! SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
    ! Inputs: Sal, TempCi, Pdbari, TCi, pHi, K1, K2
    ! Outputs: OmegaCa, OmegaAr
    ! This calculates omega, the solubility ratio, for calcite and aragonite.
    ! This is defined by: Omega = [CO3--]*[Ca++]./Ksp,
    !       where Ksp is the solubility product (either KCa or KAr).
    ! ***********************************************************************
    ! These are from:
    ! Mucci, Alphonso, The solubility of calcite and aragonite in seawater
    !       at various salinities, temperatures, and one atmosphere total
    !       pressure, American Journal of Science 283:781-799, 1983.
    ! Ingle, S. E., Solubility of calcite in the ocean,
    !       Marine Chemistry 3:301-319, 1975,
    ! Millero, Frank, The thermodynamics of the carbonate system in seawater,
    !       Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
    ! Ingle et al, The solubility of calcite in seawater at atmospheric
    !       pressure and 35!o salinity, Marine Chemistry 1:295-307, 1973.
    ! Berner, R. A., The solubility of calcite and aragonite in seawater in
    !       atmospheric pressure and 34.5!o salinity, American Journal of
    !       Science 276:713-730, 1976.
    ! Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
    ! Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic
    !       acid, boric acid, and the pHi of seawater, Limnology and
    !       Oceanography, 13:403-417, 1968.
    ! ***********************************************************************
    ! Elements from other modules:
    !
    ! Constants:
    use fundamental_constants, only: R_gas
    ! Functions:
    use unit_conversions, only: KtoC

    ! Variable Declarations:
    implicit none

    ! Arguments
    real(kind=dp), intent(in) :: &
         S,          &  ! Salinity profile array
         T_K,        &  ! Temperature profile array [K]
         P,          &  ! Pressure profile array [dbar]
         rho,        &  ! Seawater density [kg/m^3]
         DIC,        &  ! Dissolved inorganic carbon profile array [uM]
         pH             ! pH profile array
    real(kind=dp), intent(out) :: &
         Omega_ca,   &  ! Calcite saturation state
         Omega_ar       ! Aragonite saturation state
    ! Local Variables
    real(kind=dp) :: &
         T_C,        &  ! Temperature profile array [C]
         logT_K,     &  ! Natural logarithm of T_K
         sqrtS,      &  ! Square-root of S
         CO3,        &  ! [CO3-] [mol/kg]
         Ca,         &  ! [Ca^2+] [mol/kg]
         H,          &  ! Free [H+] [mol/kg]
         KCa,        &  ! Calcite solubility [(mol/kg)^2]
         KAr,        &  ! Aragonite solubility [(mol/kg)^2]
         deltaV_KCa, &  ! delta factor for KCa pressure correction
         deltaV_KAr, &  ! delta factor for KAr pressure correction
         Kappa_KCa,  &  ! Kappa factor for KCa pressure correction
         Kappa_KAr      ! Kappa factor for KAr pressure correction

    ! Set equilibrium constants
    call set_constants(T_K, S)

    ! Precalculate quantities
    T_C = KtoC(T_K)
    logT_K = log(T_K)
    sqrtS = sqrt(S)

    ! Calculate Ca^2+:
    ! Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
    ! this is .010285 * Si / 35
    Ca = 0.02128d0 / 40.087d0 * (S / 1.80655d0)
 
    ! Calcite solubility:
    ! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    KCa = 10.0d0**(-171.9065d0 - 0.077993d0 * T_K + 2839.319d0 / T_K     &
         + 71.595d0 * logT_K / log(10.0d0)                           &
         + (-0.77712d0 + 0.0028426d0 * T_K + 178.34d0 / T_K) * sqrtS &
         - 0.07711d0 * S + 0.0041249d0 * sqrtS * S)

    ! Aragonite solubility:
    ! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    KAr = 10.0d0**(-171.945d0 - 0.077993d0 * T_K + 2903.293d0 / T_K &
         + 71.595d0 * logT_K / log(10.0d0) &
         + (-0.068393d0 + 0.0017276d0 * T_K + 88.135d0 / T_K) * sqrtS &
         - 0.10018d0 * S + 0.0059415d0 * sqrtS * S)

    ! Pressure correction for calcite:
    ! Ingle, Marine Chemistry 3:301-319, 1975
    ! same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
    ! has typos (-.5304, -.3692, and 10^3 for Kappa factor)
    deltaV_KCa = -48.76d0 + 0.5304d0 * T_C
    Kappa_KCa  = (-11.76d0 + 0.3692d0 * T_C) / 1000.0d0
    KCa = KCa * exp((-deltaV_KCa + 0.5d0 * Kappa_KCa * P) * P / (R_gas * T_K))

    ! Pressure correction for aragonite:
    ! Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
    ! same as Millero, GCA 1995 except for typos (-.5304, -.3692,
    ! and 10^3 for Kappa factor)
    deltaV_KAr = deltaV_KCa + 2.8d0
    Kappa_KAr  = Kappa_KCa
    KAr = KAr * exp((-deltaV_KAr + 0.5d0 * Kappa_KAr * P) * P / (R_gas * T_K))

    ! Calculate Omegas:
    H = 10.0d0**(-pH)
    CO3 = (1.0d-3/rho) * DIC * K1 * K2 / (K1 * H + H * H + K1 * K2)
    Omega_ca = CO3 * Ca / KCa
    Omega_ar = CO3 * Ca / KAr

  end subroutine ca_solubility

end module carbonate

