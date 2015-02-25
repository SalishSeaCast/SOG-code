module carbonate
  ! Variable & parameter value declarations, and subroutines related
  ! to carbonate system calculations in the SOG code
  !
  ! Constants and alkalinity parametrizations taken from CO2SYS v1.1
  !  Lewis, E., and D. W. R. Wallace. 1998. Program Developed for
  !  CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information
  !  Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy,
  !  Oak Ridge, Tennessee. 
  !  http://cdiac.ornl.gov/oceans/co2rprt.html
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
       K1,          &  ! First carbonate eq constant [H][HCO3]/[H2CO3]
       K2,          &  ! Second carbonate eq constant [H][CO3]/[HCO3]
       KB,          &  ! Borate equilibrium constant [H][BO2]/[HBO2]
       KW,          &  ! Dissociation of water [H][OH]
       KP1,         &  ! First phosphate eq constant [H][H2PO4]/[H3PO4]
       KP2,         &  ! Second phosphate eq constant [H][HPO4]/[H2PO4]
       KP3,         &  ! Third phosphate eq constant [H][PO4]/[HPO4]
       KSi,         &  ! Silicate equilibrium constant [H][SiO(OH)3]/[Si(OH)4]
       KS,          &  ! Sulfate equilibrium constant [H][SO4]/[HSO4]
       KF,          &  ! Fluoride equilibrium constant [H][F]/[HF]
       TB,          &  ! Total borate
       TS,          &  ! Total sulfate
       TF,          &  ! Total fluoride
       FugFac,      &  !
       ! Subroutines
       calc_carbonate
       ! Private to module
       ! set_contants, pressure_corrections

  ! Variable declarations:
  real(kind=dp) :: &
       ! Carbonate system constants
       K0,          &  ! CO2 solubility
       K1,          &  ! First carbonate eq constant [H][HCO3]/[H2CO3]
       K2,          &  ! Second carbonate eq constant [H][CO3]/[HCO3]
       KB,          &  ! Borate equilibrium constant [H][BO2]/[HBO2]
       KW,          &  ! Dissociation of water [H][OH]
       KP1,         &  ! First phosphate eq constant [H][H2PO4]/[H3PO4]
       KP2,         &  ! Second phosphate eq constant [H][HPO4]/[H2PO4]
       KP3,         &  ! Third phosphate eq constant [H][PO4]/[HPO4]
       KSi,         &  ! Silicate equilibrium constant [H][SiO(OH)3]/[Si(OH)4]
       KS,          &  ! Sulfate equilibrium constant [H][SO4]/[HSO4]
       KF,          &  ! Fluoride equilibrium constant [H][F]/[HF]
       TB,          &  ! Total borate
       TS,          &  ! Total sulfate
       TF,          &  ! Total fluoride
       FugFac

contains

  subroutine set_constants(watertype, Sal, TempK, Pdbar)
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
    character(len=*), intent(in) :: &
         watertype      ! 'sea' or 'fresh'?
    real(kind=dp), intent(in) :: &
         TempK,      &  ! Temperature [K]
         Sal,        &  ! Practical salinity [PSU]
         Pdbar          ! Pressure [dbar]
    ! Local Variables
    character(len=4) :: &
         seaconstants   ! 'L-00' or 'M-10'?
    real(kind=dp) :: &
         Pbar,     &  ! Pressure [bar]
         TempK100, &  ! Temperature [K] / 100
         logTempK, &  ! Natural log of temperature [K]
         sqrSal,   &  ! Square root salinity
         IonS,     &  ! Ionic strength
         sqrIonS,  &  ! Square root ionic strength
         SWStoTOT, &  ! Seawater scale to total scale conversion factor
         lnK0,     &  ! Natural log K0
         pK1,      &  ! -log10 K1
         pK2,      &  ! -log10 K2
         lnK1,     &  ! Natural log K1
         lnK2,     &  ! Natural log K2
         lnKB,     &  ! Natural log KB
         lnKW,     &  ! Natural log KW
         lnKP1,    &  ! Natural log KP1
         lnKP2,    &  ! Natural log KP2
         lnKP3,    &  ! Natural log KP3
         lnKSi,    &  ! Natural log KSi
         lnKS,     &  ! Natural log KS
         lnKF         ! Natural log KF
    ! Millero 2010 Variables
    real(kind=dp) :: &
         pK10,     &  !
         A1,       &  !
         B1,       &  !
         C1,       &  !
         pK20,     &  !
         A2,       &  !
         B2,       &  !
         C2

    ! Preallocate common operations
    Pbar = Pdbar / 10
    TempK100 = TempK / 100.0d0
    logTempK = log(TempK)
    sqrSal   = sqrt(Sal)

    ! Calculate IonS:
    ! This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
    IonS = 19.924d0 * Sal / (1000.0d0 - 1.005d0 * Sal)
    sqrIonS = sqrt(IonS)

    ! CALCULATE SEAWATER CONSTITUENTS USING EMPIRCAL FITS
    ! Calculate total borate:
    ! Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
    ! this is 0.000416 * Sali / 35 = 0.0000119 * Sali
    ! TB = (0.000232d0 / 10.811d0) * (Sal / 1.80655d0) ! in mol/kg-SW
    TB = 0.0004157d0 * Sal / 35.0d0    ! in mol/kg-SW

    ! Calculate total sulfate:
    ! Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
    ! this is .02824 * Sali / 35 = .0008067 * Sali
    TS = (0.14d0 / 96.062d0) * (Sal / 1.80655d0)     ! in mol/kg-SW

    ! Calculate total fluoride:
    ! Riley, J. P., Deep-Sea Research 12:219-220, 1965:
    ! this is .000068 * Sali / 35 = .00000195 * Sali
    ! Approximate [F-] of Fraser River is 3 umol/kg
    TF = max((0.000067d0 / 18.998d0) * (Sal / 1.80655d0), 3.0d-6) ! in mol/kg-SW
    
    ! CALCULATE EQUILIBRIUM CONSTANTS (SW scale)
    ! Calculate KS:
    ! Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
    ! The goodness of fit is .021.
    ! It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    ! TYPO on p. 121: the constant e9 should be e8.
    ! This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
    lnKS = -4276.1d0 / TempK + 141.328d0 - 23.093d0 * logTempK +              &
         (-13856.0d0 / TempK + 324.57d0 - 47.986d0 * logTempK) * sqrIonS + &
         (35474.0d0 / TempK - 771.54d0 + 114.723d0 * logTempK) * IonS +       &
         (-2698.0d0 / TempK) * sqrIonS * IonS + (1776.0d0 / TempK) * IonS**2 
    KS = exp(lnKS)     &     ! this is on the free pH scale in mol/kg-H2O
         * (1.0d0 - 0.001005d0 * Sal)     ! convert to mol/kg-SW

    ! Calculate KF:
    ! Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
    lnKF = 1590.2d0 / TempK - 12.641d0 + 1.525d0 * sqrIonS
    KF   = exp(lnKF)   &     ! this is on the free pH scale in mol/kg-H2O
         * (1.0d0 - 0.001005d0 * Sal)     ! convert to mol/kg-SW

    ! Calculate pH scale conversion factors ( NOT pressure-corrected)
    SWStoTOT  = (1 + TS / KS) / (1 + TS / KS + TF / KF)

    ! Calculate K0:
    ! Weiss, R. F., Marine Chemistry 2:203-215, 1974.
    lnK0 = -60.2409d0 + 93.4517d0 / TempK100 + 23.3585d0 * log(TempK100) + &
         Sal * (0.023517d0 - 0.023656d0 * TempK100 + 0.0047036d0 * TempK100**2)
    K0 = exp(lnK0)               ! this is in mol/kg-SW/atm

    ! Calculate K1 & K2:
    if (watertype .eq. 'sea') then

       ! Which constants?
       seaconstants = 'M-10'

       if (seaconstants .eq. 'L-00') then

          ! From Lueker, Dickson, Keeling, 2000
          ! This is Mehrbach's data refit after conversion to the total scale,
          ! for comparison with their equilibrator work. 
          ! Mar. Chem. 70 (2000) 105-119
          ! Total scale and kg-sw
          pK1 = 3633.86d0 / TempK - 61.2172d0 + 9.6777d0 * logTempK - &
               0.011555d0 * Sal + 0.0001152d0 * Sal**2
          K1 = 10**(-1.0d0 * pK1)  & ! this is on the total pH scale in mol/kg-SW
               / SWStoTOT            ! convert to SWS pH scale

          pK2 = 471.78d0 / TempK + 25.929d0 - 3.16967d0 * logTempK - &
               0.01781d0 * Sal + 0.0001122d0 * Sal**2
          K2 = 10**(-1.0d0 * pK2)  & ! this is on the total pH scale in mol/kg-SW
               / SWStoTOT            ! convert to SWS pH scale

       else if (seaconstants .eq. 'M-10') then

          ! From Millero, 2010, also for estuarine use.
          ! Marine and Freshwater Research, v. 61, p. 139–142.
          ! Fits through compilation of real seawater titration results:
          ! Mehrbach et al. (1973), Mojica-Prieto & Millero (2002),
          ! Millero et al. (2006)
          ! Constants for K's on the SWS;
          ! This is from page 141
          pK10 = -126.34048d0 + 6320.813d0 / TempK + 19.568224d0 * log(TempK)
          ! This is from their table 2, page 140.
          A1   = 13.4038d0 * Sal**0.5 + 0.03206d0 * Sal - 5.242d-5 * Sal**2
          B1   = -530.659d0 * Sal**0.5 - 5.8210d0 * Sal
          C1   = -2.0664d0 * Sal**0.5
          pK1  = pK10 + A1 + B1 / TempK + C1 * log(TempK)
          K1   = 10**(-pK1)
          ! This is from page 141
          pK20 =  -90.18333d0 + 5143.692d0 / TempK + 14.613358d0 * log(TempK)
          ! This is from their table 3, page 140.
          A2   = 21.3728d0 * Sal**0.5 + 0.1218d0 * Sal - 3.688d-4 * Sal**2
          B2   = -788.289d0 * Sal**0.5 - 19.189d0 * Sal
          C2   = -3.374d0 * Sal**0.5
          pK2  = pK20 + A2 + B2 / TempK + C2 * log(TempK)
          K2   = 10**(-pK2)

       else
          write (stdout, *) 'parameter type in carbonate.f90:', &
               'Unexpected value: ', seaconstants
          call exit(1)
       end if
       
    elseif (watertype .eq. 'fresh') then
       ! PURE WATER CASE
       ! Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
       ! K1 from refit data from Harned and Davis,
       ! J American Chemical Society, 65:2030-2037, 1943.
       ! K2 from refit data from Harned and Scholes,
       ! J American Chemical Society, 43:1706-1709, 1941.
       ! This is only to be used for Sal=0 water
       ! (note the absence of S in the below formulations)
       ! These are the thermodynamic Constants:
       lnK1 = 290.9097d0 - 14554.21d0 / TempK - 45.0575d0 * logTempK
       K1 = exp(lnK1)
       lnK2 = 207.6548d0 - 11843.79d0 / TempK - 33.6485d0 * logTempK
       K2 = exp(lnK2)

       ! Assign Millero 2010 vars to zero so compiler doesn't get angry
       pK10 = 0.0d0
       A1   = 0.0d0
       B1   = 0.0d0
       C1   = 0.0d0
       pK20 = 0.0d0
       A2   = 0.0d0
       B2   = 0.0d0
       C2   = 0.0d0

    else
       write (stdout, *) 'parameter type in carbonate.f90:', &
            'Unexpected value: ', watertype
       call exit(1)
    end if
    
    ! Calculate KW:
    ! Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
    ! his check value of 1.6 umol/kg-SW should be 6.2
    lnKW = 148.9802d0 - 13847.26d0 / TempK - 23.6521d0 * logTempK +      &
         (-5.977d0 + 118.67d0 / TempK + 1.0495d0 * logTempK) * sqrSal -  &
         0.01615d0 * Sal
    KW = exp(lnKW) ! this is on the SWS pH scale in (mol/kg-SW)^2
    
    ! KB, KP1, KP2, KP3, & KSi
    if (watertype .eq. 'sea') then
       ! Calculate KB:
       ! Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
       lnKB = (-8966.9d0 - 2890.53d0 * sqrSal - 77.942d0 * Sal +            &
            1.728d0 * sqrSal * Sal - 0.0996d0 * Sal**2) / TempK +           &
            148.0248d0 + 137.1942d0 * sqrSal + 1.62142d0 * Sal +            &
            (-24.4344d0 - 25.085d0 * sqrSal - 0.2474d0 * Sal) * logTempK +  &
            0.053105d0 * sqrSal * TempK
       KB = exp(lnKB)         &   ! this is on the total pH scale in mol/kg-SW
            / SWStoTOT            ! convert to SWS pH scale

       ! Calculate KP1, KP2, KP3, and KSi:
       ! Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
       ! KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
       ! KSi was given on the SWS pH scale in molal units.
       lnKP1 = -4576.752d0 / TempK + 115.54d0 - 18.453d0 * logTempK +    &
            (-106.736d0 / TempK + 0.69171d0) * sqrSal +                  &
            (-0.65643d0 / TempK - 0.01844d0) * Sal
       KP1 = exp(lnKP1)

       lnKP2 = -8814.715d0 / TempK + 172.1033d0 - 27.927d0 * logTempK +  &
            (-160.34d0 / TempK + 1.3566d0) * sqrSal +                    &
            (0.37335d0 / TempK - 0.05778d0) * Sal
       KP2 = exp(lnKP2)

       lnKP3 = -3070.75d0 / TempK - 18.126d0 + &
            (17.27039d0 / TempK + 2.81197d0) * sqrSal + &
            (-44.99486d0 / TempK - 0.09984d0) * Sal
       KP3 = exp(lnKP3)

       lnKSi = -8904.2d0 / TempK + 117.4d0 - 19.334d0 * logTempK +       &
            (-458.79d0 / TempK + 3.5913d0) * sqrIonS +                   &
            (188.74d0 / TempK - 1.5998d0) * IonS +                       &
            (-12.1652d0 / TempK + 0.07871d0) * IonS**2
       KSi = exp(lnKSi)   &      ! this is on the SWS pH scale in mol/kg-H2O
            * (1.0d0 - 0.001005d0 * Sal)      ! convert to mol/kg-SW
    elseif (watertype .eq. 'fresh') then
       KB  = 0.0d0
       KP1 = 0.0d0
       KP2 = 0.0d0
       KP3 = 0.0d0
       KSi = 0.0d0
    else
       write (stdout, *) 'parameter type in carbonate.f90:', &
            'Unexpected value: ', watertype
       call exit(1)
    end if

    call pressure_corrections(watertype, TempK, Pbar)

    ! Re-calculate pH scale conversion factors (now pressure-corrected)
    SWStoTOT  = (1 + TS / KS) / (1 + TS / KS + TF / KF)

    K1  = K1  * SWStoTOT
    K2  = K2  * SWStoTOT
    KW  = KW  * SWStoTOT
    KB  = KB  * SWStoTOT
    KP1 = KP1 * SWStoTOT
    KP2 = KP2 * SWStoTOT
    KP3 = KP3 * SWStoTOT
    KSi = KSi * SWStoTOT

  end subroutine set_constants


  subroutine pressure_corrections(watertype, TempK, Pbar)
    ! Calculate pressure corrections for constants defined in set_constants
    !
    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use io_unit_defs, only: stdout
    ! Parameters:
    use fundamental_constants, only: R_gas
    ! Functions:
    use unit_conversions, only: KtoC

    ! Variable Declarations:
    implicit none

    ! Arguments
    character(len=*), intent(in) :: &
         watertype      ! 'sea' or 'fresh'?
    real(kind=dp), intent(in) :: &
         TempK,      &  ! Temperature [K]
         Pbar           ! Pressure [dbar]
    ! Local Variables
    real(kind=dp) :: &
         RT,         &  ! Temperature * gas constant (ideal gas law)
         TempC,      &  ! Temperature [deg C]
         deltaV,     &  ! Molal volume
         Kappa,      &  ! Compressibility
         lnK1fac,    &  ! K1  pressure correction factor
         lnK2fac,    &  ! K2  pressure correction factor
         lnKBfac,    &  ! KB  pressure correction factor
         lnKWfac,    &  ! KW  pressure correction factor
         lnKFfac,    &  ! KF  pressure correction factor
         lnKSfac,    &  ! KS  pressure correction factor
         lnKP1fac,   &  ! KP1 pressure correction factor
         lnKP2fac,   &  ! KP2 pressure correction factor
         lnKP3fac,   &  ! KP3 pressure correction factor
         lnKSifac,   &  ! KSi pressure correction factor
         Delta,      &  !
         b              !

    ! Temperature and gas constant
    RT    = R_gas * TempK
    TempC = KtoC(TempK)

    ! Fugacity Factor
    Delta = 57.7d0 - 0.118d0 * TempK
    b = -1636.75d0 + 12.0408d0 * TempK - 0.0327957d0 * TempK**2 + &
         3.16528d0 * 1.0d-5 * TempK**3
    FugFac = exp((b + 2.0d0 * Delta) * 1.01325d0 / RT)

    if (watertype .eq. 'sea') then
       ! Pressure effects on K1 & K2:
       ! These are from Millero, 1995.
       ! They are the same as Millero, 1979 and Millero, 1992.
       ! They are from data of Culberson and Pytkowicz, 1968.
       deltaV = -25.5d0 + 0.1271d0 * TempC
       ! deltaV = deltaV - .151 * (Sali - 34.8)   ! Millero, 1979
       Kappa = (-3.08d0 + 0.0877d0 * TempC) / 1000.0d0
       ! Kappa = Kappa - .578 * (Sali - 34.8) / 1000  ! Millero, 1979
       lnK1fac = (-deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
       ! The fits given in Millero, 1983 are somewhat different.
       deltaV = -15.82d0 - 0.0219d0 * TempC
       ! deltaV = deltaV + .321 * (Sali - 34.8)  ! Millero, 1979
       Kappa = (1.13d0 - 0.1475d0 * TempC) / 1000.0d0
       ! Kappa = Kappa - .314 * (Sali - 34.8) / 1000  ! Millero, 1979
       lnK2fac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
       ! The fit given in Millero, 1983 is different.
       ! Not by a lot for deltaV, but by much for Kappa.

       ! Pressure effects on KW:
       ! This is from Millero, 1983 and his programs CO2ROY(T).BAS.
       deltaV = -20.02d0 + 0.1119d0 * TempC - 0.001409d0 * TempC**2
       ! Millero, 1992 and Millero, 1995 have:
       Kappa = (-5.13d0 + 0.0794d0 * TempC) / 1000.0d0  ! Millero, 1983
       ! Millero, 1995 has this too, but Millero, 1992 is different.
       lnKWfac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
       ! Millero, 1979 does not list values for these.
    elseif (watertype .eq. 'fresh') then
       ! K1, K2 and KW pressure corrections from Millero, 1983
       deltaV  = -30.54d0 + 0.1849d0 * TempC - 0.0023366d0 * TempC**2
       Kappa   = (-6.22d0 + 0.1368d0 * TempC - 0.001233d0 * TempC**2) / 1.0d3
       lnK1fac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
       deltaV  = -29.81d0 + 0.115d0 * TempC - 0.001816d0 * TempC**2
       Kappa   = (-5.74d0 + 0.093d0 * TempC - 0.001896d0 * TempC**2) / 1.0d3
       lnK2fac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
       deltaV  =  -25.6d0 + 0.2324d0 * TempC - 0.0036246d0 * TempC**2
       Kappa   = (-7.33d0 + 0.1368d0 * TempC - 0.001233d0 * TempC**2) / 1.0d3
       lnKWfac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
       ! NOTE the temperature dependence of KappaK1 and KappaKW
       ! for fresh water in Millero, 1983 are the same.
    else
       write (stdout, *) 'parameter type in carbonate.f90:', &
            'Unexpected value: ', watertype
       call exit(1)
    end if
    
    ! Pressure effects on KB:
    ! This is from Millero, 1979.
    ! It is from data of Culberson and Pytkowicz, 1968.
    deltaV = -29.48d0 + 0.1622d0 * TempC - 0.002608d0 * TempC**2
    ! deltaV = -28.56 + .1211*TempCi - .000321*TempCi*TempCi  ! Millero, 1983
    ! deltaV = -29.48 + .1622 * TempCi + .295 * (Sali - 34.8) ! Millero, 1992
    ! deltaV = -29.48 - .1622*TempCi - .002608*TempCi*TempCi  ! Millero, 1995
    ! deltaV = deltaV + .295 * (Sali - 34.8)                  ! Millero, 1979
    Kappa = -2.84d0 / 1000.0d0    ! Millero, 1979
    ! Millero, 1992 and Millero, 1995 also have this.
    ! Kappa = Kappa + .354 * (Sali - 34.8) / 1000  ! Millero, 1979
    ! Kappa = (-3.0d0 + .0427d0 * TempCi) / 1000   ! Millero, 1983
    lnKBfac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT

    ! Pressure effects on KF & KS:
    ! These are from Millero, 1995, which is the same as Millero, 1983.
    ! It is assumed that KF and KS are on the free pH scale.
    deltaV  = -9.78d0 - 0.009d0 * TempC - 0.000942d0 * TempC**2
    Kappa   = (-3.91d0 + 0.054d0 * TempC) / 1000.0d0
    lnKFfac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
    deltaV  = -18.03d0 + 0.0466d0 * TempC + 0.000316d0 * TempC**2
    Kappa   = (-4.53d0 + 0.09d0 * TempC) / 1000.0d0
    lnKSfac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT

    ! Correct KP1, KP2, & KP3 for pressure:
    ! The corrections for KP1, KP2, and KP3 are from Millero, 1995,
    ! which are the same as Millero, 1983.
    deltaV   = -14.51d0 + 0.1211d0 * TempC - 0.000321d0 * TempC**2
    Kappa    = (-2.67d0 + 0.0427d0 * TempC) / 1000.0d0
    lnKP1fac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
    deltaV   = -23.12d0 + 0.1758d0 * TempC - 0.002647d0 * TempC**2
    Kappa    = (-5.15d0 + 0.09d0 * TempC) / 1000.0d0
    lnKP2fac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT
    deltaV   = -26.57d0 + 0.202d0 * TempC - 0.003042d0 * TempC**2
    Kappa    = (-4.08d0 + 0.0714d0 * TempC) / 1000.0d0
    lnKP3fac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT

    ! Pressure effects on KSi:
    ! The only mention of this is Millero, 1995 where it is stated that the
    ! values have been estimated from the values of boric acid. HOWEVER,
    ! there is no listing of the values in the table.
    ! I used the values for boric acid from above.
    deltaV   = -29.48d0 + 0.1622d0 * TempC - 0.002608d0 * TempC**2
    Kappa    = -2.84 / 1000.0d0
    lnKSifac = (-1.0d0 * deltaV + 0.5d0 * Kappa * Pbar) * Pbar / RT

    ! Correct K's for pressure here:
    K1  = K1  * exp(lnK1fac)
    K2  = K2  * exp(lnK2fac)
    KW  = KW  * exp(lnKWfac)
    KB  = KB  * exp(lnKBfac)
    KF  = KF  * exp(lnKFfac)
    KS  = KS  * exp(lnKSfac)
    KP1 = KP1 * exp(lnKP1fac)
    KP2 = KP2 * exp(lnKP2fac)
    KP3 = KP3 * exp(lnKP3fac)
    KSi = KSi * exp(lnKSifac)

  end subroutine pressure_corrections


  subroutine calc_carbonate(watertype, typein, typeout, par1, par2, &
       rho, S, T, P, PO4, Si, parout)
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
    character(len=*), intent(in) :: &
         watertype,  &    ! 'sea' or 'fresh'?
         typein,     &    ! 'DIC' or 'pH'?
         typeout          ! 'DIC', 'pH', 'pCO2', or 'Omega'?
    real(kind=dp), intent(in) :: &
         par1,       &    ! Total alkalinity [ueq/L]
         par2,       &    ! Either DIC [uM C] or pH [total scale]
         rho,        &    ! Seawater/freshwater density [kg/m^3]
         S,          &    ! Practical Salinity [PSU]
         T,          &    ! Temperature [K]
         P,          &    ! Pressure [dbar]
         PO4,        &    ! Phosphate [uM]
         Si               ! Silicate [uM]
    real(kind=dp), intent(out) :: &
         parout           ! 'DIC', 'pH', 'pCO2', or 'Omega'?
    ! Local variables
    ! Carbonate system properties
    real(kind=dp) :: &
         Alk,        &    ! Total alkalinity [eq/kg]
         DIC,        &    ! DIC [mol/kg]
         TP,         &    ! Total phosphate [mol/kg]
         TSi,        &    ! Total silicate [mol/kg]
         pH,         &    ! pH [total scale]
         H_total          ! Total [H+]

    ! Set equilibrium constants
    call set_constants(watertype, S, T, P)

    ! Convert from uM to mol/kg
    Alk = par1 * 1.0d-3 / rho
    TP  = PO4 * 1.0d-3 / rho
    TSi = Si * 1.0d-3 / rho

    if (typein .eq. 'DIC') then

       if (typeout .eq. 'DIC') then
          write (stdout, *) 'output parameter type in carbonate.f90:', &
            typeout, 'is not appropriate for input type: ', typein
          call exit(1)
       else

          ! Convert from uM to mol/kg
          DIC = par2 * 1.0d-3 / rho

          ! pH from DIC and Alkalinity
          call CalculatepHfromTATC(Alk, DIC, TP, TSi, pH)
          
          if (typeout .eq. 'pH') then

             ! Parout equals pH
             parout = pH

          elseif (typeout .eq. 'pCO2') then

             ! [H+] from pH
             H_total = 10.0d0**(-1.0d0 * pH)

             ! Calculate pCO2 as defined in CDIAC Best Practices 2007,
             ! PICES Special Publication 3, Dickson et al. eds.
             ! (Ch 2 pp 9-10, EQs 60, 68, & 70)
             parout = DIC * H_total**2 / &
                  (H_total**2 + K1 * H_total + K1 * K2) / K0

             ! Parout equals pCO2
             parout = parout * 1e6 / FugFac

          elseif (typeout .eq. 'Omega_A') then

             ! Parout equals Omega_A
             call ca_solubility(S, T, P, DIC, pH, parout)
             
          else
             write (stdout, *) 'output parameter type in carbonate.f90:', &
                  'Unexpected value: ', typeout
             call exit(1)
          endif

       end if

    elseif (typein .eq. 'pH') then

       if (typeout .eq. 'pH') then
          write (stdout, *) 'output parameter type in carbonate.f90:', &
               typeout, 'is not appropriate for input type: ', typein
          call exit(1)
       else

          ! DIC from Alkalinity and pH
          call CalculateTCfromTApH(Alk, par2, TP, TSi, DIC)

          if (typeout .eq. 'DIC') then

             ! Parout equals DIC
             parout = DIC * rho * 1.0d3
          
          else if (typeout .eq. 'pCO2') then

             ! [H+] from pH
             H_total = 10.0d0**(-1.0d0 * pH)

             ! Calculate CO2* as defined in CDIAC Best Practices 2007,
             ! PICES Special Publication 3, Dickson et al. eds.
             ! (Ch 2 pp 9-10, EQs 60, 68, & 70)
             parout = DIC * H_total**2 / &
                  (H_total**2 + K1 * H_total + K1 * K2) / K0

             ! Parout equals pCO2
             parout = parout * 1e6 / FugFac

          elseif (typeout .eq. 'Omega_A') then

             ! Parout equals Omega_A
             call ca_solubility(S, T, P, DIC, pH, parout)

          else
             write (stdout, *) 'output parameter type in carbonate.f90:', &
                  'Unexpected value: ', typeout
             call exit(1)

          end if

       end if

    else
       write (stdout, *) 'parameter type in carbonate.f90:', &
            'Unexpected value: ', typein
       call exit(1)
    endif

  end subroutine calc_carbonate


  subroutine CalculateTCfromTApH(TA, pH, TP, TSi, TC)
    ! SUB CalculateTCfromTApH, version 02.03, 10-10-97, written by Ernie Lewis.
    ! Inputs: TA, pH, K(), T()
    ! Output: TC
    ! This calculates TC from TA and pH.
    ! Though it is coded for H on the total pH scale, for the pH values
    ! occuring in seawater (pH > 6) it will be equally valid on any pH scale
    ! (H terms negligible) as long as the K Constants are on that scale.
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
         TA,       &    ! Total alkalinity [eq kg-1]
         pH,       &    ! pH
         TP,       &    ! Phosphate [mol kg-1]
         TSi            ! Silicate [mol kg-1]
    real(kind=dp), intent(out) :: &
         TC             ! Either DIC [mol kg-1] or pH [total scale]
    ! Local variables
    ! Carbonate system properties
    real(kind=dp) :: &
         H, CAlk, BAlk, OH, PhosTop, PhosBot, PAlk, SiAlk, &
         FREEtoTOT, Hfree, HSO4, HF

    ! Calculate alkalinities
    H         = 10.0d0**(-1.0d0 * pH)
    BAlk      = TB * KB / (KB + H)
    OH        = KW / H
    PhosTop   = KP1 * KP2 * H + 2.0d0 * KP1 * KP2 * KP3 - H * H * H
    PhosBot   = H * H * H + KP1 * H * H + KP1 * KP2 * H + KP1 * KP2 * KP3
    PAlk      = TP * PhosTop / PhosBot
    SiAlk     = TSi * KSi / (KSi + H)
    FREEtoTOT = (1 + TS / KS) ! pH scale conversion factor
    Hfree     = H / FREEtoTOT ! for H on the total scale
    HSO4      = TS / (1 + KS / Hfree) ! since KS is on the free scale
    HF        = TF /(1 + KF / Hfree) ! since KF is on the free scale
    CAlk      = TA - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF
    TC        = CAlk * (H * H + K1 * H + K1 * K2) / (K1 * (H + 2.0d0 * K2))

  end subroutine CalculateTCfromTApH


  subroutine CalculatepHfromTATC(TA, TC, TP, TSi, pH)
    ! SUB CalculatepHfromTATC, version 04.01, 10-13-96, written by Ernie Lewis.
    ! Inputs: TA, TC, TP, TSi
    ! Output: pH
    ! This calculates pH from TA and TC using K1 and K2 by Newton's method.
    ! It tries to solve for the pH at which Residual = 0.
    ! The starting guess is pH = 8.
    ! Though it is coded for H on the total pH scale, for the pH values
    ! occuring in seawater (pH > 6) it will be equally valid on any pH scale
    ! (H terms negligible) as long as the K Constants are on that scale.
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
         TA,       &    ! Total alkalinity [ueq/L]
         TC,       &    ! Either DIC [uM C] or pH [total scale]
         TP,       &    ! Phosphate [uM]
         TSi            ! Silicate [uM]
    real(kind=dp), intent(out) :: &
         pH             ! pH
    ! Local variables
    ! Carbonate system properties
    real(kind=dp) :: &
         pHGuess, pHTol, ln10, deltapH, H, Denom, CAlk, BAlk, OH, &
         PhosTop, PhosBot, PAlk, SiAlk, FREEtoTOT, Hfree, HSO4, HF, &
         Residual, Slope

    ! Set iteration parameters
    pHGuess     = 8.0d0       ! this is the first guess
    pHTol       = 1.0d-4      ! tolerance for iterations end
    ln10        = log(10.0d0) !
    pH = pHGuess ! creates a vector holding the first guess for all samples
    deltapH     = pHTol + 1.0d0

    ! Begin iteration to find pH
    do while (abs(deltapH) > pHTol)
       H         = 10.d0**(-1.0d0 * pH)
       Denom     = (H * H + K1 * H + K1 * K2)
       CAlk      = TC * K1 * (H + 2.0d0 * K2) / Denom
       BAlk      = TB * KB / (KB + H)
       OH        = KW / H
       PhosTop   = KP1 * KP2 * H + 2 * KP1 * KP2 * KP3 - H * H * H
       PhosBot   = H * H * H + KP1 * H * H + KP1 * KP2 * H + KP1 * KP2 * KP3
       PAlk      = TP * PhosTop / PhosBot
       SiAlk     = TSi * KSi / (KSi + H)
       FREEtoTOT = (1 + TS / KS)         ! pH scale conversion factor
       Hfree     = H / FREEtoTOT         ! for H on the total scale
       HSO4      = TS / (1 + KS / Hfree) ! since KS is on the free scale
       HF        = TF / (1 + KF / Hfree) ! since KF is on the free scale
       Residual  = TA - CAlk - BAlk - OH - PAlk - SiAlk + Hfree + HSO4 + HF
       ! find Slope dTA/dpH (not exact, but keeps all important terms)
       Slope     = ln10 * (TC * K1 * H * (H * H + K1 * K2 + 4.0d0 * H * K2) &
            / Denom / Denom + BAlk * H / (KB + H) + OH + H)
       deltapH   = Residual / Slope ! this is Newton's method
       ! to keep the jump from being too big
       do while (abs(deltapH) > 1)
          deltapH = deltapH / 2.0d0
       end do
       pH = pH + deltapH ! Is on the same scale as K1 and K2 were calculated
    end do

  end subroutine CalculatepHfromTATC


  subroutine ca_solubility(S, TempK, P, DIC, pH, Omega_ar)
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
         TempK,      &  ! Temperature profile array [K]
         P,          &  ! Pressure profile array [dbar]
         DIC,        &  ! Dissolved inorganic carbon profile array [mol kg-1]
         pH             ! pH profile array
    real(kind=dp), intent(out) :: &
         Omega_ar       ! Aragonite saturation state
    ! Local Variables
    real(kind=dp) :: &
         TempC,      &  ! Temperature profile array [C]
         logTempK,   &  ! Natural logarithm of T_K
         sqrtS,      &  ! Square-root of S
         CO3,        &  ! [CO3-] [mol/kg]
         Ca,         &  ! [Ca^2+] [mol/kg]
         H,          &  ! Free [H+] [mol/kg]
         KCa,        &  ! Calcite solubility [(mol/kg)^2]
         KAr,        &  ! Aragonite solubility [(mol/kg)^2]
         deltaV_KCa, &  ! Molal volume for KCa pressure correction
         deltaV_KAr, &  ! Molal volume for KAr pressure correction
         Kappa_KCa,  &  ! Compressibility for KCa pressure correction
         Kappa_KAr,  &  ! Compressibility for KAr pressure correction
         Omega_ca       ! Calcite saturation state

    ! Precalculate quantities
    TempC = KtoC(TempK)
    logTempK = log(TempK)
    sqrtS = sqrt(S)

    ! Calculate Ca^2+:
    ! Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
    ! this is .010285 * Si / 35
    Ca = 0.02128d0 / 40.087d0 * (S / 1.80655d0)
 
    ! Calcite solubility:
    ! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    KCa = 10.0d0**(-171.9065d0 - 0.077993d0 * TempK + 2839.319d0 / TempK     &
         + 71.595d0 * logTempK / log(10.0d0)                           &
         + (-0.77712d0 + 0.0028426d0 * TempK + 178.34d0 / TempK) * sqrtS &
         - 0.07711d0 * S + 0.0041249d0 * sqrtS * S)

    ! Aragonite solubility:
    ! Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    KAr = 10.0d0**(-171.945d0 - 0.077993d0 * TempK + 2903.293d0 / TempK &
         + 71.595d0 * logTempK / log(10.0d0) &
         + (-0.068393d0 + 0.0017276d0 * TempK + 88.135d0 / TempK) * sqrtS &
         - 0.10018d0 * S + 0.0059415d0 * sqrtS * S)

    ! Pressure correction for calcite:
    ! Ingle, Marine Chemistry 3:301-319, 1975
    ! same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
    ! has typos (-.5304, -.3692, and 10^3 for Kappa factor)
    deltaV_KCa = -48.76d0 + 0.5304d0 * TempC
    Kappa_KCa  = (-11.76d0 + 0.3692d0 * TempC) / 1000.0d0
    KCa = KCa * exp((-deltaV_KCa + 0.5d0 * Kappa_KCa * P) * P / (R_gas * TempK))

    ! Pressure correction for aragonite:
    ! Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
    ! same as Millero, GCA 1995 except for typos (-.5304, -.3692,
    ! and 10^3 for Kappa factor)
    deltaV_KAr = deltaV_KCa + 2.8d0
    Kappa_KAr  = Kappa_KCa
    KAr = KAr * exp((-deltaV_KAr + 0.5d0 * Kappa_KAr * P) * P / (R_gas * TempK))

    ! Calculate Omegas:
    H = 10.0d0**(-pH)
    CO3 = DIC * K1 * K2 / (K1 * H + H * H + K1 * K2)
    Omega_ca = CO3 * Ca / KCa
    Omega_ar = CO3 * Ca / KAr

  end subroutine ca_solubility

end module carbonate

