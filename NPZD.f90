module NPZD
  ! Type definitions, variables, and subroutines for the biological
  ! model.

  use precision_defs, only: dp

  implicit none
  private

  public :: &
       ! Type definitions:
       bins,              &
       ! Variables:
       PZ,                &
       PZ_bins,           &
       flagellates,       &  ! Can flagellates can influence other biology?
       remineralization,  &  ! Is there a remineralization loop?
       microzooplankton,  &  ! Active or not
       strong_limitation, &  ! single species light limitation
       micro,             &  ! micro-plankton growth profile arrays
       ! Diagnostics:
       uptake,            &  ! NO, NH, and PC uptake arrays
       Mesozoo,           &  ! Mesozooplankton profile array
       f_ratio,           &  ! Ratio of new to total production profile
       ! Subroutines:
       init_NPZD,            &  ! Allocate memory for biology model
                                ! variables, and read parameter values
                                ! from infile.
       derivs,               &  ! Calculate the derivatives of the
                                ! biology quantities for odeint() to
                                ! use to calculate their values at the
                                ! next time step.
       dalloc_NPZD_variables    ! Deallocate memory used by NPZD variables.

  ! Private type definitions:
  !
  ! Indices for quantities in PZ vector
  type :: bins
     integer :: &
          Quant, micro, nano, pico, zoo, NO, NH, Si, det, &
          D_DON, D_PON, D_refr, D_bSi, &
          DIC, Oxy, Alk, D_DOC, D_POC  ! Chemistry bins
  end type bins
  !
  ! Rate parameters for phytoplankton
  type :: rate_para_phyto
     real(kind=dp) :: R, &  ! max growth rate for light limitation
          Iopt, & ! optimum light (Steeves scheme)
          maxtemp, & ! maximum temp for growth
          temprange, & ! range below tempmax that temp is limiting
          Q10exp, & ! exponent for Q10 tempeffect
          gamma, & ! loss parameter under light limitation
          k, & ! half saturation constant
          ! preference for NO3 over NH (always less than 1 as NH is preferred)
          kapa, &
          gamma_o, & ! exp strength of NH inhit of NO3 uptake
          N_o, & ! overall half saturation constant????
          N_x, & ! exponent in inhibition equation
          Si_ratio, & ! silicon to nitrogen ratio in phyto
          K_Si, & ! half-saturation for Si
          Rm      ! natural mortality rate
  end type rate_para_phyto
  !
  ! Rate parameters for Zooplankton
  type :: rate_para_zoo
     real(kind=dp) :: &
          winterconc, &      ! background, winter conc (mesozoo only)
          summerconc      ! relative summer level (mesozoo only)
     real(kind=dp), dimension(1:3) :: &  ! three summer Gaussians
          sumpeakval, &   ! peak values
          sumpeakpos, &   ! day of peak
          sumpeakwid      ! width of peak
     real(kind=dp) :: &
          R, &               ! max ingestion rate
          eff, &             ! assimilation efficiency
          Rm, &              ! natural mortality
          excr, &            ! excretion
          PredSlope, &       ! overall limit from predation
          HalfSat, &         ! overall half-saturation
          MicroPref, &       ! fractional preference for diatoms
          MicroPredSlope, &  ! limit from predation
          MicroHalfSat, &    ! half-saturation for getting eaten
          NanoPref, &       ! fractional preference for diatoms
          NanoPredSlope, &  ! limit from predation
          NanoHalfSat, &    ! half-saturation for getting eaten
          PicoPref, &       ! fractional preference for picoplankton
          PicoPredSlope, &  ! limit from predation
          PicoHalfSat, &    ! half-saturation for getting eaten
          PON_Pref, &        ! fraction preference for PON
          PON_PredSlope, &   ! limit from predation
          PON_HalfSat, &        ! half-saturation for PON
          Z_Pref, &        ! fraction preference for Z
          Z_PredSlope, &   ! limit from predation
          Z_HalfSat        ! half-saturation for Z
  end type rate_para_zoo
  !
  ! New parameters for waste
  type :: nloss_param
     real(kind=dp) :: NH, DON, PON, Ref, Bsi
  end type nloss_param
  !
  ! Remineralization rates of ammonium and detritus
  type :: remin_rates
     real(kind=dp) :: &
          NH,    &  ! Ammonium remineralization rate [1/s]
          D_DON, &  ! Dissolved organic nitrogen detritus remin rate [1/s]
          D_PON, &  ! Particulate organic nitrogen detritus remin rate [1/s]
          D_bSi     ! Biogenic silicon detritus remineralization rate [1/s]
  end type remin_rates
  !
  ! Nitrogen compound uptake diagnostics
  type :: uptake_
     real(kind=dp), dimension(:), pointer :: &
          NO, &  ! Nitrate uptake profile
          NH, &  ! Ammonium uptake profile
          PC     ! Differential carbon uptake profile
  end type uptake_
  !
  ! Plankton growth
  type :: grow
     real(kind=dp), dimension(:), pointer :: light, new
  end type grow
  type :: plankton_growth
     type(grow) :: growth
     real(kind=dp), dimension(:), pointer :: Nlimit
  end type plankton_growth

  ! Public variable declarations:
  !
  ! Biological model logicals (turned parts of the model on or off)
  logical :: &
       flagellates,      &  ! Can flagellates can influence other biology?
       remineralization, &  ! Is there a remineralization loop?
       microzooplankton, &  ! use a microzooplantkon pool?
       strong_limitation    ! impose single species strong light limitation
  real(kind=dp), dimension(:), allocatable :: Mesozoo
  type(plankton_growth) :: micro  ! Micro-plankton growth profile arrays

  ! Private variable declarations:
  !
  ! Remineralization rates:
  type(remin_rates) :: remin
  !
  ! Indices for quantities (e.g. phyto, nitrate, etc.) in PZ vector
  type(bins):: PZ_bins
  data PZ_bins%Quant /10/    ! Size of the biology (Quantities vs Detritus)
  data PZ_bins%micro /1/     ! Position of Diatoms (micro plankton)
  data PZ_bins%nano  /2/     ! Position of Meso rub (nano plankton)
  data PZ_bins%pico  /3/     ! Position of Flagellates (pico plankton)
  data PZ_bins%zoo   /4/     ! Position of microzooplankton
  data PZ_bins%NO    /5/     ! Position of Nitrate
  data PZ_bins%NH    /6/     ! Position of Ammonium
  data PZ_bins%Si    /7/     ! Position of Silicon
  data PZ_bins%DIC   /8/     ! Position of Dissolved Inorganic Carbon
  data PZ_bins%Oxy   /9/     ! Position of Dissolved Oxygen
  data PZ_bins%Alk   /10/    ! Position of Alkalinity
  data PZ_bins%det   /11/    ! Start of detritus
  data PZ_bins%D_DOC   /11/  ! Position of dissolved organic carbon detritus
  data PZ_bins%D_POC   /12/  ! Position of particulate organic carbon detritus
  data PZ_bins%D_DON   /13/  ! Position of dissolved organic nitrogen detritus
  data PZ_bins%D_PON   /14/  ! Position of particulate organic nitrogen detritus
  data PZ_bins%D_bSi   /15/  ! Position of dissolved biogenic silcon detritus
  data PZ_bins%D_refr  /16/  ! Position of refractory nitrogen detritus
  !
  ! PZ vector
  real(kind=dp), dimension(:), allocatable :: PZ
  !
  ! Biological rate parameters
  type(rate_para_phyto) :: rate_micro, rate_nano, rate_pico
  type(rate_para_zoo) :: rate_mesozoo, rate_mesorub, rate_uzoo
  type(nloss_param) :: &
       frac_waste_DNM, & ! Diatom natural mortality
       frac_waste_NNM, & ! Nano natural mortality
       frac_waste_FNM, & ! Flagellates (pico) natural mortality
       frac_waste_MNM, & ! Mesozoo natural mortality
       frac_waste_MEX, & ! Mesozoo excretion
       frac_waste_ZNM, & ! uZoo natural mortality
       frac_waste_ZEX, & ! uZoo excretion
       frac_waste_DEM, & ! Diatoms eaten by mesozoo
       frac_waste_NEM, & ! Nano eaten by mesozoo
       frac_waste_FEM, & ! Flagellates (pico) eaten by mesozoo
       frac_waste_PEM, & ! PON eaten by Mesozoo
       frac_waste_ZEM, & ! uZ eaten by Mesozoo
       frac_waste_PEZ, & ! PON eaten by uZoo
       frac_waste_DEZ, & ! Diatoms eaten by uZoo
       frac_waste_NEZ, & ! Nano eaten by uZoo
       frac_waste_FEZ, & ! Flagelattes (pico) eaten by uZoo
       frac_waste_ZEZ, & ! uZoo eaten by uZoo
       frac_waste_FEN    ! Flagelattes (pico) eaten by Mesorub (nano)
  !
  ! Nitrogen compound uptake diagnotics
  type(uptake_) :: uptake
  !
  ! Nitrogen remineralization diagnostics
  real(kind=dp), dimension(:), allocatable :: &
       remin_NH, &  ! Total remineralization to ammonium
       NH_oxid,  &  ! Bacterial oxidation of NH4 to NO3
       f_ratio  ! Ratio of new to total production profile
  !
  ! Nano & pico plankton growth profile arrays
  type(plankton_growth) :: nano, pico

  ! Meso and microzoo grazing terms
  real(kind=dp), dimension(:), allocatable ::  &
         Meso_mort_micro,   &  ! Mesozooplankton grazing microphyto
         Meso_mort_nano,    &  ! Mesozooplankton grazing nanophyto
         Meso_mort_pico,    &  ! Mesozooplankton grazing picophyto
         Meso_graz_PON,     &  ! Mesozooplankton grazing PON
         Meso_mort_Z,       &  ! Mesozooplankton grazing microzoo
         Mesorub_mort_pico, &  ! Mesodinium rubrum grazing picophyto
         uZoo_mort_micro,   &  ! Micozooplankton grazing microphyto
         uZoo_mort_nano,    &  ! Micozooplankton grazing nanophyto
         uZoo_mort_pico,    &  ! Micozooplankton grazing picophyto
         uZoo_graz_PON,     &  ! Micozooplankton grazing PON
         uZoo_graz_Z           ! Micozooplankton grazing microzoo

contains

  subroutine read_biology_params
    ! Read the biology model parameters from the input file
    use input_processor, only: getpard, getpardv, getpari, getparl

    implicit none

    ! Logicals that control biology model options
    flagellates = getparl('flagellates_on')
    remineralization = getparl('remineralization')
    microzooplankton = getparl('use microzooplankton')
    strong_limitation = getparl('single species light')

    ! Biological rate parameters
    ! zooplankton rates
    ! winter concentration/ summer concentration
    rate_mesozoo%winterconc = getpard('Mesozoo, winter conc')
    rate_mesozoo%summerconc = getpard('Mesozoo, summer conc')
    ! parameters governing the summer conc
    call getpardv('Mesozoo, summer peak mag',3,rate_mesozoo%sumpeakval)
    call getpardv('Mesozoo, summer peak pos', 3,rate_mesozoo%sumpeakpos)
    call getpardv('Mesozoo, summer peak wid',3, rate_mesozoo%sumpeakwid)
    ! max igestion rate
    rate_mesozoo%R = getpard('Mesozoo, max ingestion')
    ! assimilation efficiency
    rate_mesozoo%eff = getpard("Mesozoo, assimil. eff")
    ! natural mortality
    rate_mesozoo%Rm = getpard("Mesozoo, nat mort")
    ! excretion
    rate_mesozoo%excr = getpard("Mesozoo, excretion")
    ! limit from predation
    rate_mesozoo%PredSlope = getpard('Mesozoo, pred slope')
    ! half saturation
    rate_mesozoo%HalfSat = getpard('Mesozoo, half-sat')
    ! preference for diatoms
    rate_mesozoo%MicroPref = getpard('Mesozoo, pref for diatoms')
    ! limit from predation
    rate_mesozoo%MicroPredSlope = getpard('Mesozoo, micro pred slope')
    ! half saturation
    rate_mesozoo%MicroHalfSat = getpard('Mesozoo, micro half-sat')
    ! preference for nano
    rate_mesozoo%NanoPref = getpard('Mesozoo, pref for nano')
    ! limit from predation
    rate_mesozoo%NanoPredSlope = getpard('Mesozoo, nano pred slope')
    ! half saturation
    rate_mesozoo%NanoHalfSat = getpard('Mesozoo, nano half-sat')
    ! preference for pico
    rate_mesozoo%PicoPref = getpard('Mesozoo, pref for pico')
    ! limit from predation
    rate_mesozoo%PicoPredSlope = getpard('Mesozoo, pico pred slope')
    ! half saturation
    rate_mesozoo%PicoHalfSat = getpard('Mesozoo, pico half-sat')
    ! preference for PON
    rate_mesozoo%PON_Pref = getpard('Mesozoo, pref for PON')
    ! limit from predation
    rate_mesozoo%PON_PredSlope = getpard('Mesozoo, PON pred slope')
    ! half saturation
    rate_mesozoo%PON_HalfSat = getpard('Mesozoo, PON half-sat')
    ! uzoo preference
    rate_mesozoo%Z_Pref = getpard('Mesozoo, pref for uZoo')
    ! limit from predation
    rate_mesozoo%Z_PredSlope = getpard('Mesozoo, uZoo pred slope')
    ! half saturation
    rate_mesozoo%Z_HalfSat = getpard('Mesozoo, uZoo half-sat')
    ! Mesodinium rates
    rate_mesorub%R = getpard("Mesorub, max ingestion")
    rate_mesorub%eff = getpard("Mesorub, assimilation eff")
    rate_mesorub%PicoHalfSat = getpard('Mesorub, nano half-sat')
    rate_mesorub%PicoPredSlope = getpard('Mesorub, nano predslope')
    ! Microzoo rates
    rate_uzoo%R = getpard("Microzoo, max ingestion")
    rate_uzoo%eff = getpard("Microzoo, assimil. eff")
    rate_uzoo%Rm = getpard("Microzoo, nat mort")
    rate_uzoo%excr = getpard("Microzoo, excretion")
    rate_uzoo%PredSlope = getpard('Microzoo, pred slope')
    rate_uzoo%HalfSat = getpard('Microzoo, half-sat')
    rate_uzoo%PicoPref = getpard('Microzoo, pref for Pico')
    rate_uzoo%PicoPredSlope = getpard('uzoo, Pico pred slope')
    rate_uzoo%PicoHalfSat = getpard('uzoo, Pico half-sat')
    rate_uzoo%MicroPref = getpard('Microzoo, pref for Micro')
    rate_uzoo%MicroPredSlope = getpard('uzoo, Micro pred slope')
    rate_uzoo%MicroHalfSat = getpard('Microzoo, Micro half-sat')
    rate_uzoo%NanoPref = getpard('Microzoo, pref for nano')
    rate_uzoo%NanoPredSlope = getpard('Microzoo, nano pred slope')
    rate_uzoo%NanoHalfSat = getpard('Microzoo, nano half-sat')
    rate_uzoo%PON_Pref = getpard('Microzoo, pref for PON')
    rate_uzoo%PON_PredSlope = getpard('Microzoo, PON pred slope')
    rate_uzoo%PON_HalfSat = getpard('Microzoo, PON half-sat')
    rate_uzoo%Z_Pref = getpard('Microzoo, pref for uZoo')
    rate_uzoo%Z_PredSlope = getpard('Microzoo, uZoo pred slope')
    rate_uzoo%Z_HalfSat = getpard('Microzoo, uZoo half-sat')

    ! max growth rate for light limitation
    rate_micro%R = getpard('Micro, max growth')
    rate_nano%R = getpard('Nano, max growth')
    rate_pico%R = getpard('Pico, max growth')
    ! optimum light level
    rate_micro%Iopt = getpard('Micro, I_opt')
    rate_nano%Iopt = getpard('Nano, I_opt')
    rate_pico%Iopt = getpard('Pico, I_opt')
    ! maximum temp for growth
    rate_micro%maxtemp = getpard('Micro, max temp')
    rate_nano%maxtemp = getpard('Nano, max temp')
    rate_pico%maxtemp = getpard('Pico, max temp')
    ! temp range of limitation below max
    rate_micro%temprange = getpard('Micro, temp range')
    rate_nano%temprange = getpard('Nano, temp range')
    rate_pico%temprange = getpard('Pico, temp range')
    ! exponent for Q_10
    rate_micro%Q10exp = getpard('Micro, Q10 exp')
    rate_nano%Q10exp = getpard('Nano, Q10 exp')
    rate_pico%Q10exp = getpard('Pico, Q10 exp')
    ! loss parameter under light limitation
    rate_micro%gamma = getpard('Micro, gamma loss')
    rate_nano%gamma = getpard('Nano, gamma loss')
    rate_pico%gamma = getpard('Pico, gamma loss')
    ! NO3 half saturation constant
    rate_micro%k = getpard('Micro, NO3 k')
    rate_nano%k = getpard('Nano, NO3 k')
    rate_pico%k = getpard('Pico, NO3 k')
    ! preference for NO3 over NH
    rate_micro%kapa = getpard('Micro, kapa')
    rate_nano%kapa = getpard('Nano, kapa')
    rate_pico%kapa = getpard('Pico, kapa')
    ! exp strength of NH inhit of NO3 uptake
    rate_micro%gamma_o = getpard('Micro, NH inhib')
    rate_nano%gamma_o = getpard('Nano, NH inhib')
    rate_pico%gamma_o = getpard('Pico, NH inhib')
    ! overall half saturation constant????
    rate_micro%N_o = getpard('Micro, N_o')
    rate_nano%N_o = getpard('Nano, N_o')
    rate_pico%N_o = getpard('Pico, N_o')
    ! exponent in inhibition equation
    rate_micro%N_x = getpard('Micro, N_x')
    rate_nano%N_x = getpard('Nano, N_x')
    rate_pico%N_x = getpard('Pico, N_x')
    ! silicon to nitrogen ratio in phyto
    rate_micro%Si_ratio = getpard('Micro, Si ratio')
    rate_nano%Si_ratio = getpard('Nano, Si ratio')
    rate_pico%Si_ratio = getpard('Pico, Si ratio')
    ! half-saturation for Silicon
    rate_micro%K_Si = getpard('Micro, K Si')
    rate_nano%K_Si = getpard('Nano, K Si')
    rate_pico%K_Si = getpard('Pico, K Si')
    ! natural mortality rate
    rate_micro%Rm = getpard('Micro, nat mort')
    rate_nano%Rm = getpard('Nano, nat mort')
    rate_pico%Rm = getpard('Pico, nat mort')
    ! Remineralization rates:
    ! Ammonium
    remin%NH = getpard('NH remin rate')
    ! DON detritus
    remin%D_DON = getpard('DON remin rate')
    ! PON detritus
    remin%D_PON = getpard('PON remin rate')
    ! Biogenic silicon detritus
    remin%D_bSi = getpard('bSi remin rate')
    ! Microphyto(diatom) natural mortality
    frac_waste_DNM%NH = getpard('Waste, dnm, NH')
    frac_waste_DNM%DON = getpard('Waste, dnm, DON')
    frac_waste_DNM%PON = getpard('Waste, dnm, PON')
    frac_waste_DNM%Ref = getpard('Waste, dnm, Ref')
    frac_waste_DNM%Bsi = getpard('Waste, dnm, Bsi')
    ! Nanophyto natural mortality
    frac_waste_NNM%NH = getpard('Waste, nnm, NH')
    frac_waste_NNM%DON = getpard('Waste, nnm, DON')
    frac_waste_NNM%PON = getpard('Waste, nnm, PON')
    frac_waste_NNM%Ref = getpard('Waste, nnm, Ref')
    frac_waste_NNM%Bsi = getpard('Waste, nnm, Bsi')
    ! Picophyto (flagellates) natural mortality
    frac_waste_FNM%NH = getpard('Waste, fnm, NH')
    frac_waste_FNM%DON = getpard('Waste, fnm, DON')
    frac_waste_FNM%PON = getpard('Waste, fnm, PON')
    frac_waste_FNM%Ref = getpard('Waste, fnm, Ref')
    frac_waste_FNM%Bsi = getpard('Waste, fnm, Bsi')
    ! Mesozoo natural mortality
    frac_waste_MNM%NH = getpard('Waste, mnm, NH')
    frac_waste_MNM%DON = getpard('Waste, mnm, DON')
    frac_waste_MNM%PON = getpard('Waste, mnm, PON')
    frac_waste_MNM%Ref = getpard('Waste, mnm, Ref')
    frac_waste_MNM%Bsi = getpard('Waste, mnm, Bsi')
    ! Mesozoo excretion
    frac_waste_MEX%NH = getpard('Waste, mex, NH')
    frac_waste_MEX%DON = getpard('Waste, mex, DON')
    frac_waste_MEX%PON = getpard('Waste, mex, PON')
    frac_waste_MEX%Ref = getpard('Waste, mex, Ref')
    frac_waste_MEX%Bsi = getpard('Waste, mex, Bsi')
    ! uZoo natural mortality
    frac_waste_ZNM%NH = getpard('Waste, znm, NH')
    frac_waste_ZNM%DON = getpard('Waste, znm, DON')
    frac_waste_ZNM%PON = getpard('Waste, znm, PON')
    frac_waste_ZNM%Ref = getpard('Waste, znm, Ref')
    frac_waste_ZNM%Bsi = getpard('Waste, znm, Bsi')
    ! uZoo excretion
    frac_waste_ZEX%NH = getpard('Waste, zex, NH')
    frac_waste_ZEX%DON = getpard('Waste, zex, DON')
    frac_waste_ZEX%PON = getpard('Waste, zex, PON')
    frac_waste_ZEX%Ref = getpard('Waste, zex, Ref')
    frac_waste_ZEX%Bsi = getpard('Waste, zex, Bsi')
    ! Microphyto(diatom) eaten by Mesozoo
    frac_waste_DEM%NH = getpard('Waste, dem, NH')
    frac_waste_DEM%DON = getpard('Waste, dem, DON')
    frac_waste_DEM%PON = getpard('Waste, dem, PON')
    frac_waste_DEM%Ref = getpard('Waste, dem, Ref')
    frac_waste_DEM%Bsi = getpard('Waste, dem, Bsi')
    ! Nanophyto eaten by Mesozoo
    frac_waste_NEM%NH = getpard('Waste, nem, NH')
    frac_waste_NEM%DON = getpard('Waste, nem, DON')
    frac_waste_NEM%PON = getpard('Waste, nem, PON')
    frac_waste_NEM%Ref = getpard('Waste, nem, Ref')
    frac_waste_NEM%Bsi = getpard('Waste, nem, Bsi')
    ! Picophyto eaten by Mesozoo
    frac_waste_FEM%NH = getpard('Waste, fem, NH')
    frac_waste_FEM%DON = getpard('Waste, fem, DON')
    frac_waste_FEM%PON = getpard('Waste, fem, PON')
    frac_waste_FEM%Ref = getpard('Waste, fem, Ref')
    frac_waste_FEM%Bsi = getpard('Waste, fem, Bsi')
    ! PON eaten by Mesozoo
    frac_waste_PEM%NH = getpard('Waste, pem, NH')
    frac_waste_PEM%DON = getpard('Waste, pem, DON')
    frac_waste_PEM%PON = getpard('Waste, pem, PON')
    frac_waste_PEM%Ref = getpard('Waste, pem, Ref')
    frac_waste_PEM%Bsi = getpard('Waste, pem, Bsi')
    ! uZoo eaten by Mesozoo
    frac_waste_ZEM%NH = getpard('Waste, zem, NH')
    frac_waste_ZEM%DON = getpard('Waste, zem, DON')
    frac_waste_ZEM%PON = getpard('Waste, zem, PON')
    frac_waste_ZEM%Ref = getpard('Waste, zem, Ref')
    frac_waste_ZEM%Bsi = getpard('Waste, zem, Bsi')
    ! Diatoms eaten by uZoo
    frac_waste_DEZ%NH = getpard('Waste, dez, NH')
    frac_waste_DEZ%DON = getpard('Waste, dez, DON')
    frac_waste_DEZ%PON = getpard('Waste, dez, PON')
    frac_waste_DEZ%Ref = getpard('Waste, dez, Ref')
    frac_waste_DEZ%Bsi = getpard('Waste, dez, Bsi')
    ! Nano eaten by uZoo
    frac_waste_NEZ%NH = getpard('Waste, nez, NH')
    frac_waste_NEZ%DON = getpard('Waste, nez, DON')
    frac_waste_NEZ%PON = getpard('Waste, nez, PON')
    frac_waste_NEZ%Ref = getpard('Waste, nez, Ref')
    frac_waste_NEZ%Bsi = getpard('Waste, nez, Bsi')
    ! Pico (flagellates) eaten by uZoo
    frac_waste_FEZ%NH = getpard('Waste, fez, NH')
    frac_waste_FEZ%DON = getpard('Waste, fez, DON')
    frac_waste_FEZ%PON = getpard('Waste, fez, PON')
    frac_waste_FEZ%Ref = getpard('Waste, fez, Ref')
    frac_waste_FEZ%Bsi = getpard('Waste, fez, Bsi')
    ! PON eaten by uZoo
    frac_waste_PEZ%NH = getpard('Waste, pez, NH')
    frac_waste_PEZ%DON = getpard('Waste, pez, DON')
    frac_waste_PEZ%PON = getpard('Waste, pez, PON')
    frac_waste_PEZ%Ref = getpard('Waste, pez, Ref')
    frac_waste_PEZ%Bsi = getpard('Waste, pez, Bsi')
    ! uZoo eaten by uZoo
    frac_waste_ZEZ%NH = getpard('Waste, zez, NH')
    frac_waste_ZEZ%DON = getpard('Waste, zez, DON')
    frac_waste_ZEZ%PON = getpard('Waste, zez, PON')
    frac_waste_ZEZ%Ref = getpard('Waste, zez, Ref')
    frac_waste_ZEZ%Bsi = getpard('Waste, zez, Bsi')
    ! Pico (flagellates) eaten by Mesorub
    frac_waste_FEN%NH = getpard('Waste, fen, NH')
    frac_waste_FEN%DON = getpard('Waste, fen, DON')
    frac_waste_FEN%PON = getpard('Waste, fen, PON')
    frac_waste_FEN%Ref = getpard('Waste, fen, Ref')
    frac_waste_FEN%Bsi = getpard('Waste, fen, Bsi')
  end subroutine read_biology_params


  subroutine init_NPZD(M, PZ_length)
    ! Initialize biology model.

    implicit none

    ! Argument:
    integer, intent(in) :: &
         M, &       ! Number of grid points
         PZ_length  ! Length of consolidated array of biology quantity values

    call alloc_NPZD_variables(M, PZ_length)
    call read_biology_params
    ! If phytoplankton Si to nitrate ratio is zero set the half
    ! saturation limit to zero
    if (rate_micro%Si_ratio == 0.0d0) rate_micro%K_Si = 0.0d0
    if (rate_nano%Si_ratio == 0.0d0) rate_nano%K_Si = 0.0d0
    if (rate_pico%Si_ratio == 0.0d0) rate_pico%K_Si = 0.0d0
    ! If remineralization is turned off, set the rates to zero
    if (.not. remineralization) then
       remin%NH = 0.0d0
       remin%D_DON = 0.0d0
       remin%D_PON = 0.0d0
       remin%D_bSi = 0.0d0
    endif
  end subroutine init_NPZD


  subroutine alloc_NPZD_variables(M, PZ_length)
    ! Allocate memory for biological model arrays.
    use malloc, only: alloc_check
    implicit none

    ! Argument:
    integer :: &
         M, &  ! Number of grid points
         PZ_length
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Consolidated array of biology quantity values"
    allocate(PZ(1:PZ_length), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrogen compounds uptake diagnostic arrays"
    allocate(uptake%NO(1:M), uptake%NH(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Carbon compounds uptake diagnostic arrays"
    allocate(uptake%PC(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrogen remineralization diagnostic arrays"
    allocate(remin_NH(1:M), NH_oxid(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Mesozooplankton diagnostic array"
    allocate(Mesozoo(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ratio of new to total production profile array"
    allocate(f_ratio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Plankton growth arrays"
    allocate(micro%growth%light(M), micro%growth%new(M), &
         nano%growth%light(M), nano%growth%new(M), &
         pico%growth%light(M), pico%growth%new(M), &
         micro%Nlimit(M), nano%Nlimit(M), pico%Nlimit(M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Zooplankton grazing arrays"
    allocate(Meso_mort_micro(1:M), Meso_mort_nano(1:M),                  &
         Meso_mort_pico(1:M), Meso_graz_PON(1:M), Meso_mort_Z(1:M),      &
         Mesorub_mort_pico(1:M), uZoo_mort_micro(1:M), uZoo_mort_nano(1:M), &
         uZoo_mort_pico(1:M), uZoo_graz_PON(1:M), uZoo_graz_Z(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_NPZD_variables


  subroutine dalloc_NPZD_variables
    ! Deallocate memory for NPZD model arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Consolidated array of biology quantity values"
    deallocate(PZ, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrogen compounds uptake diagnostic arrays"
    deallocate(uptake%NO, uptake%NH, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Carbon compounds uptake diagnostic arrays"
    deallocate(uptake%PC, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrogen remineralization diagnostic arrays"
    deallocate(remin_NH, NH_oxid, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Mesozooplankton diagnostic array"
    deallocate(Mesozoo, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ratio of new to total production profile array"
    deallocate(f_ratio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Plankton growth arrays"
    deallocate(micro%growth%light, micro%growth%new, &
         nano%growth%light, nano%growth%new, &
         pico%growth%light, pico%growth%new, &
         micro%Nlimit, nano%Nlimit, pico%Nlimit, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Zooplankton grazing arrays"
    deallocate(Meso_mort_micro, Meso_mort_nano,           &
         Meso_mort_pico, Meso_graz_PON, Meso_mort_Z, Mesorub_mort_pico,     &
         uZoo_mort_micro, uZoo_mort_nano, uZoo_mort_pico, &
         uZoo_graz_PON, uZoo_graz_Z, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_NPZD_variables


  subroutine p_growth(M, NO, NH, Si, P, I_par, temp_Q10, Temp, rate, plank)
    ! Calculate the growth (light limited
    ! or nutrient limited)
    ! of either phytoplankton class which
    ! are functions of temperature
    use precision_defs, only: dp
    use unit_conversions, only: KtoC
    implicit none
    ! Arguments:
    integer, intent(in) :: M
    ! Nitrate, ammonium & silicon concentrations
    real(kind=dp), dimension(1:), intent(in) :: NO, NH, Si
    ! Plankton concentraton (either Pmicro or Pnano)
    real(kind=dp), dimension(1:), intent(in) :: P
    real(kind=dp), dimension(0:), intent(in) :: I_par  ! light
    real(kind=dp), dimension(1:), intent(in) :: temp_Q10  ! Q10 temp effect
    real(kind=dp), dimension(0:), intent(in) :: temp ! temperature
    ! Parameters of the growth equations
    type(rate_para_phyto), intent(in) :: rate
    ! Output is the growth values
    type(plankton_growth), intent(out) :: plank ! either micro or nano

    ! Local variables:
    integer :: j ! counter through depth
    real(kind=dp), dimension(1:M) :: Uc, & ! growth based on light
         Sc ! growth based on Si
    real(kind=dp), dimension(1:M) :: Oup_cell ! NO uptake assuming full light
    real(kind=dp), dimension(1:M) :: Hup_cell ! NH uptake assuming full light
    real(kind=dp), dimension(1:M) :: Rmax ! maximum growth rate for given Temp
    real(kind=dp) :: NH_effect ! ammonium effect on nutrient uptake

!!!!!!!!!!!Define growth Due to I_par!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    plank%growth%light = 0.

    Oup_cell = 0.
    Hup_cell = 0.

    DO j = 1,M

       ! maximum growth rate (before light/nutrient limitation)

       ! biological process impactedby Q10 temperature effect
       Rmax(j)= rate%R * temp_Q10(j) &
            ! and organisms can have a maximum growth temp
            * min(max(rate%maxtemp - KtoC(temp(j)), 0.0d0), &
                  rate%temprange) / &
            (rate%temprange + epsilon(rate%temprange))

       ! don't include Rmax effect until end

       if (strong_limitation) then

       plank%growth%light(j) = &
            ! Steeles scheme like but with extended high light range
            ! (Steeles scheme has built in light inhibition)
            ! 0.67, 2.7 and 1.8 are constants for making the fit
            ! Steeles like at small light and making it fit Durbin for
            ! Thalassosira nordelenski at higher light
            (1.0d0 - exp(-I_par(j) / (0.67d0 * rate%Iopt)) ) * &
            exp(-I_par(j) / (2.7d0 * rate%Iopt)) * 1.8d0

       else

       plank%growth%light(j) = &
            ! Steeles scheme like but with extended high light range
            ! (Steeles scheme has built in light inhibition)
            ! much broader pattern for a bucket of organisms
            (1.0d0 - exp(-I_par(j) / (0.33d0 * rate%Iopt)) ) * &
            (exp(-I_par(j) / (30.d0 * rate%Iopt))) * 1.06d0

    end if

       Uc(j) = (1.0d0 - rate%gamma) * plank%growth%light(j)
    END DO

!!!!!!!!!!!!!!!Define growth due to nutrients!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Silicon (vectorized)

    Sc = Si / (rate%K_Si + Si)

    ! Nitrate and Ammonium

    DO j = 1,M

       NH_effect = exp(-rate%gamma_o * NH(j))
       IF (NO(j) > epsilon(NO(j))) THEN
          Oup_cell(j) = NO(j) * rate%kapa * NH_effect / &
               (rate%k + NO(j) * rate%kapa * NH_effect + NH(j)) * &
               (NO(j) + NH(j))**rate%N_x / &
               (rate%N_o + (NO(j) + NH(j))**rate%N_x + epsilon(rate%N_o))
       ELSE
          Oup_cell(j) = 0.
       END IF
       IF (NH(j) > epsilon(NH(j))) THEN
          Hup_cell(j) = NH(j) / &
               (rate%k + NO(j) * rate%kapa * NH_effect + NH(j))* &
               (NO(j) + NH(j))**rate%N_x / &
               (rate%N_o + (NO(j) + NH(j))**rate%N_x + epsilon(rate%N_o))
       ELSE
          Hup_cell(j) = 0.
       END IF

       IF (Oup_cell(j) < 0.) THEN
          PRINT "(A)","Oup_cell(j) < 0. in NPZD.f90"
          PRINT *,Oup_cell(j)
          CALL EXIT(1)
       END IF

       IF (Hup_cell(j) < 0.) THEN
          PRINT "(A)","Hup_cell(j) < 0. in NPZD.f90"
          PRINT *,Hup_cell(j)
          CALL EXIT(1)
       END IF

       ! exponent of 1/5 follows Alain
       plank%Nlimit(j) = min((Oup_cell(j) + Hup_cell(j)), Sc(j))**0.2

    END DO

    ! Choose light limitation or nutrient limitation

    DO j = 1,M

       If (Uc(j) < 0. ) then
          plank%growth%new(j) = 0.
       else

          IF (min(Uc(j),Sc(j)) >= Oup_cell(j) + Hup_cell(j)) THEN

             !N LIMITING
             plank%growth%new(j) = Rmax(j) * (Oup_cell(j) + Hup_cell(j))

             IF (plank%growth%new(j) < 0.) THEN
                plank%growth%new(j) = 0.
             ENDIF

             uptake%NO(j) = Rmax(j) * Oup_cell(j) * P(j) + uptake%NO(j)
             uptake%NH(j) = Rmax(j) * Hup_cell(j) * P(j) + uptake%NH(j)

          else

             if (Uc(j) < Sc(j)) then
                !LIGHT LIMITING
                plank%growth%new(j) = Rmax(j) * Uc(j)
             else
                ! Si limitation
                plank%growth%new(j) = Rmax(j) * Sc(j)
             endif

             ! split the nitrogen uptake between NH and NO

             IF (plank%growth%new(j) <= Rmax(j) * Hup_cell(j)) THEN
                ! add to nutrient uptake so we combined the effects of
                ! different phyto classes
                uptake%NH(j) = plank%growth%new(j) * P(j) + uptake%NH(j)
                uptake%NO(j) = uptake%NO(j)
             ELSE
                uptake%NH(j) = Rmax(j) * Hup_cell(j) * P(j) + uptake%NH(j)
                uptake%NO(j) = (plank%growth%new(j) - Rmax(j) * Hup_cell(j)) &
                     * P(j) + uptake%NO(j)
             END IF

          END IF

          ! PC Differential Carbon Uptake
          ! Chlr Redux Factor = 0.2 (Ianson and Allen 2002)
          uptake%PC(j) = (Rmax(j) * Uc(j) - plank%growth%new(j)) &
                     * 0.2d0 * P(j) + uptake%PC(j)

       endif
    end do
  end subroutine p_growth


  subroutine derivs(M, PZ, Temp, I_par, day, dPZdt)
    ! Calculate the derivatives of the biology quantities for odeint()
    ! to use to calculate their values at the next time step.
    use precision_defs, only: dp
    use fundamental_constants, only: &
         pi, Redfield_C, Redfield_O, Redfield_NP
    use unit_conversions, only: KtoC
    use grid_mod, only: full_depth_average
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(1:), intent(in) :: &
         PZ  ! Consolidated array of biology quantity values
    real(kind=dp), dimension(0:), intent(in) :: &
         ! *** Should be able to eliminate upper bound here
         Temp,  &  ! Temperature profile [K]
         I_par     ! Photosynthetic available radiation profile [W/m^2]
    integer, intent(in) :: &
         day  ! Year-day of current time step
    real(kind=dp), dimension(1:), intent(out) :: &
         dPZdt  ! Consolidated array of time derivatives of bio quantities
    ! Local variables
    real(kind=dp), dimension(1:M) :: &
          Pmicro,         &  ! Micro phytoplankton (diatoms) profile
          Pnano,          &  ! Nano phytoplankton (mesorub) profile
          Ppico,          &  ! Pico phytoplankton (flagellates) profile
          Z,              &  ! Micro zooplankton profile
          NO,             &  ! Nitrate concentrationy
          NH,             &  ! Ammonium concentration profile
          Si,             &  ! Silicon concentration profile
          DIC,            &  ! Dissolved inorganic carbon profile
          Oxy,            &  ! Dissolved oxygen profile
          Alk,            &  ! Alkalinity profile
          D_DOC,          &  ! Dissolved organic carbon detritus profile
          D_POC,          &  ! Particulate organic carbon detritus profile
          D_DON,          &  ! Dissolved organic nitrogen detritus profile
          D_PON,          &  ! Particulate organic nitrogen detritus profile
          D_refr,         &  ! Refractory nitrogen detritus profile
          D_bSi,          &  ! Biogenic silicon detritus profile
          NatMort_micro,  &  ! Micro phytoplankton natural mortality profile
          GrazMort_micro, &  ! Micro phytoplankton grazing mortality profile
          NatMort_nano,   &  ! Nano phytoplankton natural mortality profile
          GrazMort_nano,  &  ! Nano phytoplankton grazing mortality profile
          GrazMort_pico,  &  ! Pico phytoplankton grazaing mortality profile
          NatMort_pico,   &  ! Pico phytoplankton natural mortality profile
          NatMort_uzoo,   &  ! uZoo natural mortality profile
          NatMort_mesozoo,&  ! mesozoo natural mortality profile - added by Tara
          Excr_uzoo,      &  ! uZoo excretion
          Excr_mesozoo,   &  ! mesozoo excretion - added by Tara
          Graz_PON,       &  ! PON grazing mortality profile
          GrazMort_Z,     &  ! Z grazing mortality profile
          Mesorub_eat,    &  ! Mesorub eat pico
          Microzoo_eat,   &  ! Micro zoo eat

! new waste variables
          was_NH, was_DON, was_PON, was_Ref, was_Bsi, &
          Si_remin           ! Profile of dissolution of biogenic Si detritus
    real(kind=dp), dimension(1:M) :: temp_Q10
    ! temporary variable to calculate mortality scheme
    real(kind=dp) :: food_limitation, denominator, average_prey

    integer :: &
         bPZ, &  ! Beginning index for a quantity in the PZ array
         ePZ     ! Ending index for a quantity in the PZ array

    integer :: jj ! counter through PZ


    ! Unload PZ array into local biology quantity arrays to make the
    ! formulation of the derivatives clearer.  Set any negative values
    ! to zero.
    !
    ! Micro phytoplankton
    bPZ = (PZ_bins%micro - 1) * M + 1
    ePZ = PZ_bins%micro * M
    Pmicro = PZ(bPZ:ePZ)
    where (Pmicro < 0.) Pmicro = 0.
    ! Nano phytoplankton
    if (flagellates) then
       bPZ = (PZ_bins%nano - 1) * M + 1
       ePZ = PZ_bins%nano * M
       Pnano = PZ(bPZ:ePZ)
       where (Pnano < 0.) Pnano = 0.
    else
       Pnano = 0.
    endif
    ! Pico phytoplankton
    if (flagellates) then
       bPZ = (PZ_bins%pico - 1) * M + 1
       ePZ = PZ_bins%pico * M
       Ppico = PZ(bPZ:ePZ)
       where (Ppico < 0.) Ppico = 0.
    else
       Ppico = 0.
    endif
    ! Microzooplankton
    if (microzooplankton) then
       bPZ = (PZ_bins%zoo - 1) * M + 1
       ePZ = PZ_bins%zoo * M
       Z = PZ(bPZ:ePZ)
       where (Z < 0.) Z = 0.
    else
       Z = 0.
    endif
    ! Nitrate
    bPZ = (PZ_bins%NO - 1) * M + 1
    ePZ = PZ_bins%NO * M
    NO = PZ(bPZ:ePZ)
    where (NO < 0.) NO = 0.
    ! Ammonium
    if (remineralization) then
       bPZ = (PZ_bins%NH - 1) * M + 1
       ePZ = PZ_bins%NH * M
       NH = PZ(bPZ:ePZ)
       where (NH < 0.) NH = 0.
    else
       NH = 0.
    endif
    ! Silicon
    bPZ = (PZ_bins%Si - 1) * M + 1
    ePZ = PZ_bins%Si * M
    Si = PZ(bPZ:ePZ)
    where (Si < 0.) Si = 0.
    ! Dissolved inorganic carbon
    bPZ = (PZ_bins%DIC - 1) * M + 1
    ePZ = PZ_bins%DIC * M
    DIC = PZ(bPZ:ePZ)
    where (DIC < 0.) DIC = 0.
   ! Dissolved oxygen
    bPZ = (PZ_bins%Oxy - 1) * M + 1
    ePZ = PZ_bins%Oxy * M
    Oxy = PZ(bPZ:ePZ)
    where (Oxy < 0.) Oxy = 0.
   ! Alkalinity
    bPZ = (PZ_bins%Alk - 1) * M + 1
    ePZ = PZ_bins%Alk * M
    Alk = PZ(bPZ:ePZ)
    where (Alk < 0.) Alk = 0.
    ! Dissolved organic carbon detritus
    if (remineralization) then
       bPZ = (PZ_bins%D_DOC - 1) * M + 1
       ePZ = PZ_bins%D_DOC * M
       D_DOC = PZ(bPZ:ePZ)
       where (D_DOC < 0.) D_DOC = 0.
    else
       D_DOC = 0.
    endif
    ! Particulate organic carbon detritus
    if (remineralization) then
       bPZ = (PZ_bins%D_POC - 1) * M + 1
       ePZ = PZ_bins%D_POC * M
       D_POC = PZ(bPZ:ePZ)
       where (D_POC < 0.) D_POC = 0.
    else
       D_POC = 0.
    endif
    ! Dissolved organic nitrogen detritus
    if (remineralization) then
       bPZ = (PZ_bins%D_DON - 1) * M + 1
       ePZ = PZ_bins%D_DON * M
       D_DON = PZ(bPZ:ePZ)
       where (D_DON < 0.) D_DON = 0.
    else
       D_DON = 0.
    endif
    ! Particulate organic nitrogen detritus
    if (remineralization) then
       bPZ = (PZ_bins%D_PON - 1) * M + 1
       ePZ = PZ_bins%D_PON * M
       D_PON = PZ(bPZ:ePZ)
       where (D_PON < 0.) D_PON = 0.
    else
       D_PON = 0.
    endif
    ! Refractory nitrogen detritus
    if (remineralization) then
       bPZ = (PZ_bins%D_refr - 1) * M + 1
       ePZ = PZ_bins%D_refr * M
       D_refr = PZ(bPZ:ePZ)
       where (D_refr < 0.) D_refr = 0.
    else
       D_refr = 0.
    endif
    ! Biogenic silicon detritus
    if (remineralization) then
       bPZ = (PZ_bins%D_bSi - 1) * M + 1
       ePZ = PZ_bins%D_bSi * M
       D_bSi = PZ(bPZ:ePZ)
       where (D_bSi < 0.) D_bSi = 0.
    else
       D_bSi = 0.
    endif

    ! Initialize transfer rates between the pools
    uptake%NO = 0.
    uptake%NH = 0.
    uptake%PC = 0.
    ! Initialize grazing terms
    GrazMort_micro = 0.
    GrazMort_nano  = 0.
    GrazMort_pico  = 0.
    Graz_PON       = 0.
    GrazMort_Z     = 0.
    ! new waste variables
    was_NH = 0.
    was_DON = 0.
    was_PON = 0.
    was_Ref = 0.
    was_Bsi = 0.

    remin_NH = 0.

    ! all biological processes are impacted by Q10 temp effect
    ! parameters are set for 20 degrees C.
    ! calculate it
    temp_Q10 = dexp (0.07 * (KtoC(Temp(1:M)) - 20.d0))

    ! phytoplankton growth:  growth is per cell, uptake is in uMol

    call p_growth(M, NO, NH, Si, Pmicro, I_par, temp_Q10**rate_micro%Q10exp, &
         Temp, & ! in
         rate_micro, micro)         ! in and out, in, out

    call p_growth(M, NO, NH, Si, Pnano, I_par, temp_Q10**rate_nano%Q10exp, &
         Temp, & ! in
         rate_nano, nano)              ! in and out, in, out

    call p_growth(M, NO, NH, Si, Ppico, I_par, temp_Q10**rate_pico%Q10exp, &
         Temp, & ! in
         rate_pico, pico)              ! in and out, in, out

    ! change growth to be total
    micro%growth%new = micro%growth%new * Pmicro
    nano%growth%new = nano%growth%new * Pnano
    pico%growth%new = pico%growth%new * Ppico

    ! phytoplankton natural mortality:

    NatMort_micro = (rate_micro%Rm) * temp_Q10 * Pmicro
    NatMort_nano  = (rate_nano%Rm)  * temp_Q10 * Pnano
    NatMort_pico  = (rate_pico%Rm)  * temp_Q10 * Ppico

    was_NH = was_NH + (frac_waste_DNM%NH * NatMort_micro &
         + frac_waste_NNM%NH * NatMort_nano &
         + frac_waste_FNM%NH * NatMort_pico) * (1-rate_mesozoo%eff)
    was_DON = was_DON + (frac_waste_DNM%DON * NatMort_micro &
         + frac_waste_NNM%DON * NatMort_nano &
         + frac_waste_FNM%DON * NatMort_pico) * (1-rate_mesozoo%eff)
    was_PON = was_PON + (frac_waste_DNM%PON * NatMort_micro &
         + frac_waste_NNM%PON * NatMort_nano &
         + frac_waste_FNM%PON * NatMort_pico) * (1-rate_mesozoo%eff)
    was_Ref = was_Ref + (frac_waste_DNM%Ref * NatMort_micro &
         + frac_waste_NNM%Ref * NatMort_nano &
         + frac_waste_FNM%Ref * NatMort_pico) * (1-rate_mesozoo%eff)
    was_BSi = was_BSi &
         + (frac_waste_DNM%BSi * NatMort_micro * rate_micro%Si_ratio &
         + frac_waste_NNM%Bsi * NatMort_nano * rate_nano%Si_ratio &
         + frac_waste_FNM%Bsi * NatMort_pico * rate_pico%Si_ratio) &
         * (1-rate_mesozoo%eff)

    ! Grazing processes: MesoZooplankton Fit

    Mesozoo(1:M) = rate_mesozoo%winterconc + &
         rate_mesozoo%summerconc * &
         ( sum ( rate_mesozoo%sumpeakval * &
            exp( -(day-rate_mesozoo%sumpeakpos)**2 / &
                                rate_mesozoo%sumpeakwid**2 ) ) &
         + sum ( rate_mesozoo%sumpeakval * &
            exp( -(day-rate_mesozoo%sumpeakpos-365.25)**2 / &
                                rate_mesozoo%sumpeakwid**2 ) ) &
         + sum ( rate_mesozoo%sumpeakval * &
            exp( -(day-rate_mesozoo%sumpeakpos+365.25)**2 / &
                                rate_mesozoo%sumpeakwid**2 ) ) )


    average_prey = full_depth_average(Pmicro + D_PON + Pnano + Ppico + Z)

    Mesozoo(1:M) = Mesozoo(1:M) &
         * (Pmicro(1:M) + D_PON(1:M) + Pnano(1:M) + Ppico(1:M) + Z(1:M)) &
         / ( average_prey + epsilon(average_prey))

    ! microzoo natural mortality:
    NatMort_uzoo = rate_uzoo%Rm * temp_Q10 * Z
    Excr_uzoo = rate_uzoo%excr * temp_Q10 * Z

    was_NH = was_NH + frac_waste_ZNM%NH * NatMort_uzoo &
         + frac_waste_ZEX%NH * Excr_uzoo
    was_DON = was_DON + frac_waste_ZNM%DON * NatMort_uzoo &
         + frac_waste_ZEX%DON * Excr_uzoo
    was_PON = was_PON + frac_waste_ZNM%PON * NatMort_uzoo &
         + frac_waste_ZEX%PON * Excr_uzoo
    was_Ref = was_Ref + frac_waste_ZNM%Ref * NatMort_uzoo &
         + frac_waste_ZEX%Ref * Excr_uzoo
    was_BSi = was_BSi &
         + frac_waste_ZNM%BSi * NatMort_uzoo * 0.d0 &
         + frac_waste_ZEX%Bsi * Excr_uzoo * 0.d0

    ! mesozoo natural mortality:
    NatMort_mesozoo = rate_mesozoo%Rm * temp_Q10 * Mesozoo
    Excr_mesozoo = rate_mesozoo%excr * temp_Q10 * Mesozoo

    was_NH = was_NH + frac_waste_MNM%NH * NatMort_mesozoo &
         + frac_waste_MEX%NH * Excr_mesozoo
    was_DON = was_DON + frac_waste_MNM%DON * NatMort_mesozoo &
         + frac_waste_MEX%DON * Excr_mesozoo
    was_PON = was_PON + frac_waste_MNM%PON * NatMort_mesozoo &
         + frac_waste_MEX%PON * Excr_mesozoo
    was_Ref = was_Ref + frac_waste_MNM%Ref * NatMort_mesozoo &
         + frac_waste_MEX%Ref * Excr_mesozoo
    was_BSi = was_BSi &
         + frac_waste_MNM%BSi * NatMort_mesozoo * 0.d0 &
         + frac_waste_MEX%Bsi * Excr_mesozoo * 0.d0

    do jj = 1,M
       ! global food limitation
       food_limitation = (Pmicro(jj) + D_PON(jj) + Pnano(jj) + Ppico(jj) &
            + Z(jj) - rate_mesozoo%PredSlope) / &
            (rate_mesozoo%HalfSat + Pmicro(jj) + D_PON(jj) + Pnano(jj) &
            + Ppico(jj) + Z(jj) - rate_mesozoo%PredSlope &
            + epsilon(rate_mesozoo%HalfSat))

       denominator = (rate_mesozoo%MicroPref * Pmicro(jj) + &
            rate_mesozoo%NanoPref * Pnano(jj) + &
            rate_mesozoo%PicoPref * Ppico(jj) + &
            rate_mesozoo%PON_Pref * D_PON(jj) + &
            rate_mesozoo%Z_Pref * Z(jj) + epsilon(Pmicro(jj)) )

       ! limitation based on microplankton
       Meso_mort_micro(jj) = min(rate_mesozoo%MicroPref * food_limitation &
            * Pmicro(jj) / denominator, &
            (Pmicro(jj) - rate_mesozoo%MicroPredslope) / &
            (rate_mesozoo%MicroHalfSat + Pmicro(jj) &
            - rate_mesozoo%MicroPredSlope + &
            epsilon(rate_mesozoo%MicroHalfSat)))

       ! limitation based on nanoplankton
       Meso_mort_nano(jj) = min(rate_mesozoo%NanoPref * food_limitation &
            * Pnano(jj) / denominator, &
            (Pnano(jj) - rate_mesozoo%NanoPredslope) / &
            (rate_mesozoo%NanoHalfSat + Pnano(jj) &
            - rate_mesozoo%NanoPredSlope + &
            epsilon(rate_mesozoo%NanoHalfSat)))

       ! limitation based on picoplankton
       Meso_mort_pico(jj) = min(rate_mesozoo%PicoPref * food_limitation &
            * Ppico(jj) / denominator, &
            (Ppico(jj) - rate_mesozoo%PicoPredslope) / &
            (rate_mesozoo%PicoHalfSat + Ppico(jj) &
            - rate_mesozoo%PicoPredSlope + &
            epsilon(rate_mesozoo%PicoHalfSat)))

       ! limitation based on PON
       Meso_graz_PON(jj) = min(rate_mesozoo%PON_Pref * food_limitation &
            * D_PON(jj) / denominator, &
            (D_PON(jj) - rate_mesozoo%PON_Predslope) / &
            (rate_mesozoo%PON_HalfSat + D_PON(jj) - rate_mesozoo%PON_PredSlope &
            + epsilon(rate_mesozoo%PON_HalfSat)))

       ! limitation based on Z
       Meso_mort_Z(jj) = min(rate_mesozoo%Z_Pref * food_limitation &
            * Z(jj) / denominator, &
            (Z(jj) - rate_mesozoo%Z_Predslope) / &
            (rate_mesozoo%Z_HalfSat + Z(jj) - rate_mesozoo%Z_PredSlope &
            + epsilon(rate_mesozoo%Z_HalfSat)))

       ! global corrected by individual
       food_limitation = Meso_mort_micro(jj) + Meso_mort_nano(jj) &
            + Meso_mort_pico(jj) + Meso_graz_PON(jj) + Meso_mort_Z(jj)

       Meso_mort_micro(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * &
            max(0.d0, Meso_mort_micro(jj))

       Meso_mort_nano(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * &
            max(0.d0, Meso_mort_nano(jj))

       Meso_mort_pico(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * &
            max(0.d0, Meso_mort_pico(jj))

       Meso_graz_PON(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * &
            max(0.d0, Meso_graz_PON(jj))

       Meso_mort_Z(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * &
            max(0.d0, Meso_mort_Z(jj))

    enddo

    GrazMort_micro = GrazMort_micro + Meso_mort_micro
    GrazMort_nano  = GrazMort_nano  + Meso_mort_nano
    GrazMort_pico  = GrazMort_pico  + Meso_mort_pico
    Graz_PON       = Graz_PON       + Meso_graz_PON
    GrazMort_Z     = GrazMort_Z     + Meso_mort_Z

    was_NH = was_NH + frac_waste_PEM%NH * Meso_graz_PON +     &
         frac_waste_DEM%NH  * Meso_mort_micro +               &
         frac_waste_NEM%NH  * Meso_mort_nano +                &
         frac_waste_FEM%NH  * Meso_mort_pico +                &
         frac_waste_ZEM%NH  * Meso_mort_Z
    was_DON = was_DON + frac_waste_PEM%DON * Meso_graz_PON +  &
         frac_waste_DEM%DON * Meso_mort_micro +               &
         frac_waste_NEM%DON * Meso_mort_nano +                &
         frac_waste_FEM%DON * Meso_mort_pico +                &
         frac_waste_ZEM%DON * Meso_mort_Z
    was_PON = was_PON + frac_waste_PEM%PON * Meso_graz_PON +  &
         frac_waste_DEM%PON * Meso_mort_micro +               &
         frac_waste_NEM%PON * Meso_mort_nano +                &
         frac_waste_FEM%PON * Meso_mort_pico +                &
         frac_waste_ZEM%PON * Meso_mort_Z
    was_Ref = was_Ref + frac_waste_PEM%Ref * Meso_graz_PON +  &
         frac_waste_DEM%Ref * Meso_mort_micro +               &
         frac_waste_NEM%Ref * Meso_mort_nano +                &
         frac_waste_FEM%Ref * Meso_mort_pico +                &
         frac_waste_ZEM%Ref * Meso_mort_Z
    was_BSi = was_BSi + frac_waste_PEM%BSi * Meso_graz_PON * 0.d0 +   &
         frac_waste_DEM%BSi * Meso_mort_micro * rate_micro%Si_ratio + &
         frac_waste_NEM%BSi * Meso_mort_nano * rate_nano%Si_ratio +   &
         frac_waste_FEM%BSi * Meso_mort_pico * rate_pico%Si_ratio +   &
         frac_waste_ZEM%BSi * Meso_mort_Z * 0.d0

    ! Mesodinium rubrum

    Mesorub_mort_pico = rate_mesorub%R * (Ppico - rate_mesorub%PicoPredSlope) &
         / (rate_mesorub%PicoHalfSat + Ppico - rate_mesorub%PicoPredSlope &
         + epsilon(rate_mesorub%PicoHalfSat)) &
         * Pnano * temp_Q10

    do jj=1,M
       Mesorub_mort_pico(jj) = max(Mesorub_mort_pico(jj), 0.d0)
       GrazMort_pico(jj) = GrazMort_pico(jj) + Mesorub_mort_pico(jj)
    enddo

    was_NH = was_NH + (frac_waste_FEN%NH) * Mesorub_mort_pico * &
         (1-rate_mesorub%eff)
    was_DON = was_DON + (frac_waste_FEN%DON) * Mesorub_mort_pico * &
         (1-rate_mesorub%eff)
    was_PON = was_PON + (frac_waste_FEN%PON) * Mesorub_mort_pico * &
         (1-rate_mesorub%eff)
    was_Ref = was_Ref + (frac_waste_FEN%Ref) * Mesorub_mort_pico * &
         (1-rate_mesorub%eff)
    was_BSi = was_BSi + (frac_waste_FEN%BSi) * Mesorub_mort_pico * &
         (1-rate_mesorub%eff)

    Mesorub_eat = Mesorub_mort_pico * rate_mesorub%eff

    ! Microzooplankton grazing
    do jj= 1, M
       food_limitation = (Pmicro(jj) + Pnano(jj) + Ppico(jj) + &
            D_PON(jj) + Z(jj) - rate_uzoo%PredSlope) / &
            (rate_uzoo%HalfSat + Pmicro(jj) + Pnano(jj) + Ppico(jj) &
            + D_PON(jj) + Z(jj) - rate_uzoo%PredSlope &
            + epsilon(rate_uzoo%HalfSat))

       denominator = (rate_uzoo%MicroPref * Pmicro(jj) + &
            rate_uzoo%NanoPref * Pnano(jj) + &
            rate_uzoo%PicoPref * Ppico(jj) + &
            rate_uzoo%PON_Pref * D_PON(jj) + &
            rate_uzoo%Z_Pref * Z(jj) + &
            epsilon(Pmicro(jj)) )

       uZoo_mort_micro(jj) = min(rate_uzoo%MicroPref * food_limitation &
            * Pmicro(jj) / denominator, &
            (Pmicro(jj) - rate_uzoo%MicroPredslope) / &
            (rate_uzoo%MicroHalfSat + Pmicro(jj) &
            - rate_uzoo%MicroPredSlope + &
            epsilon(rate_uzoo%MicroHalfSat)) )

       uZoo_mort_nano(jj) = min(rate_uzoo%NanoPref * food_limitation &
            * Pnano(jj) / denominator, &
            (Pnano(jj) - rate_uzoo%NanoPredslope) / &
            (rate_uzoo%NanoHalfSat + Pnano(jj) &
            - rate_uzoo%NanoPredSlope + &
            epsilon(rate_uzoo%NanoHalfSat)) )

       uZoo_mort_pico(jj) = min(rate_uzoo%PicoPref * food_limitation &
            * Ppico(jj) / denominator, &
            (Ppico(jj) - rate_uzoo%PicoPredslope) / &
            (rate_uzoo%PicoHalfSat + Ppico(jj) &
            - rate_uzoo%PicoPredSlope + &
            epsilon(rate_uzoo%PicoHalfSat)) )

       uZoo_graz_PON(jj) = min(rate_uzoo%PON_Pref * food_limitation &
            * D_PON(jj) / denominator, &
            (D_PON(jj) - rate_uzoo%PON_Predslope) / &
            (rate_uzoo%PON_HalfSat + D_PON(jj) &
            - rate_uzoo%PON_PredSlope + &
            epsilon(rate_uzoo%PON_HalfSat)) )

       uZoo_graz_Z(jj) = min(rate_uzoo%Z_Pref * food_limitation &
            * Z(jj) / denominator, &
            (Z(jj) - rate_uzoo%Z_Predslope) / &
            (rate_uzoo%Z_HalfSat + Z(jj) &
            - rate_uzoo%Z_PredSlope + &
            epsilon(rate_uzoo%Z_HalfSat)) )

       uZoo_mort_micro(jj) = rate_uzoo%R * temp_Q10(jj) &
            * Z(jj) * max(0.d0, uZoo_mort_micro(jj))
       uZoo_mort_nano(jj) = rate_uzoo%R * temp_Q10(jj) &
            * Z(jj) * max(0.d0, uZoo_mort_nano(jj))
       uZoo_mort_pico(jj) = rate_uzoo%R * temp_Q10(jj) &
            * Z(jj) * max(0.d0, uZoo_mort_pico(jj))
       uZoo_graz_PON(jj) = rate_uzoo%R * temp_Q10(jj) &
            * Z(jj) * max(0.d0, uZoo_graz_PON(jj))
       uZoo_graz_Z(jj) = rate_uzoo%R * temp_Q10(jj) &
            * Z(jj) * max(0.d0, uZoo_graz_Z(jj))

    enddo
    
    GrazMort_micro = GrazMort_micro + uZoo_mort_micro
    GrazMort_nano  = GrazMort_nano  + uZoo_mort_nano
    GrazMort_pico  = GrazMort_pico  + uZoo_mort_pico
    Graz_PON       = Graz_PON       + uZoo_graz_PON
    GrazMort_Z     = GrazMort_Z     + uZoo_graz_Z
    MicroZoo_eat = (uZoo_mort_micro + uZoo_mort_nano + uzoo_mort_pico &
         + uZoo_graz_PON + uZoo_graz_Z) * rate_uZoo%eff

    was_NH = was_NH + (frac_waste_PEZ%NH * UZoo_graz_PON + &
         frac_waste_DEZ%NH * uZoo_mort_micro + &
         frac_waste_NEZ%NH * uZoo_mort_nano + &
         frac_waste_FEZ%NH * uZoo_mort_pico + &
         frac_waste_ZEZ%NH * uZoo_graz_Z) * (1-rate_uZoo%eff)
    was_DON = was_DON + (frac_waste_PEZ%DON * UZoo_graz_PON + &
         frac_waste_DEZ%DON * uZoo_mort_micro + &
         frac_waste_NEZ%DON * uZoo_mort_nano + &
         frac_waste_FEZ%DON * uZoo_mort_pico + &
         frac_waste_ZEZ%DON * uZoo_graz_Z) * (1-rate_uZoo%eff)
    was_PON = was_PON + (frac_waste_PEZ%PON * UZoo_graz_PON + &
         frac_waste_DEZ%PON * uZoo_mort_micro + &
         frac_waste_NEZ%PON * uZoo_mort_nano + &
         frac_waste_FEZ%PON * uZoo_mort_pico + &
         frac_waste_ZEZ%PON * uZoo_graz_Z)  * (1-rate_uZoo%eff)
    was_Ref = was_Ref + (frac_waste_PEZ%Ref * UZoo_graz_PON + &
         frac_waste_DEZ%Ref * uZoo_mort_micro + &
         frac_waste_NEZ%Ref * uZoo_mort_nano + &
         frac_waste_FEZ%Ref * uZoo_mort_pico + &
         frac_waste_ZEZ%Ref * uZoo_graz_Z)  * (1-rate_uZoo%eff)
    was_BSi = was_BSi + frac_waste_PEZ%BSi * UZoo_graz_PON * 0.d0 + &
         frac_waste_DEZ%BSi * uZoo_mort_micro * rate_micro%Si_ratio + &
         frac_waste_NEZ%BSi * uZoo_mort_nano * rate_nano%Si_ratio + &
         frac_waste_FEZ%BSi * uZoo_mort_pico * rate_pico%Si_ratio + &
         frac_waste_ZEZ%BSi * uZoo_graz_Z * 0.d0

    ! Remineralization:
    !
    ! Bacterial oxidation of NH3 to NO pool proportional to NH^2
    NH_oxid = remin%NH * NH**2 * temp_Q10
    ! Dissolved organic nitrogen
    remin_NH = remin_NH + remin%D_DON * D_DON * temp_Q10
    ! Particulate organic nitrogen
    remin_NH = remin_NH + remin%D_PON * D_PON * temp_Q10
    ! Biogenic silcon
    Si_remin = remin%D_bSi * D_bSi * temp_Q10

    IF (MINVAL(remin_NH) < 0.) THEN
       PRINT "(A)","remin_NH < 0. in derivs.f90"
       PRINT *,remin_NH
       CALL EXIT(1)
    END IF

    ! calculate the f-ratio
    do jj = 1,M
       IF (uptake%NO(jj) + uptake%NH(jj) > 0.) THEN
          f_ratio(jj) = uptake%NO(jj) /  &
               (uptake%NO(jj) + uptake%NH(jj))
       ELSE
          f_ratio(jj) = 0.
       END IF
    END DO

    ! Build the derivatives array
    !
    ! Initialize the array
    dPZdt = 0.
    ! Micro phytoplankton
    bPZ = (PZ_bins%micro - 1) * M + 1
    ePZ = PZ_bins%micro * M
    where (Pmicro > 0.)
       dPZdt(bPZ:ePZ) = micro%growth%new - NatMort_micro &
               - GrazMort_micro
    endwhere
    ! Nano phytoplankton
    bPZ = (PZ_bins%nano - 1) * M + 1
    ePZ = PZ_bins%nano * M
    where (Pnano > 0.)
       dPZdt(bPZ:ePZ) = nano%growth%new - NatMort_nano &
               - GrazMort_nano + Mesorub_eat
    endwhere
    ! Pico phytoplankton
    bPZ = (PZ_bins%pico - 1) * M + 1
    ePZ = PZ_bins%pico * M
    where (Ppico > 0.)
       dPZdt(bPZ:ePZ) = pico%growth%new - NatMort_pico &
               - GrazMort_pico
    endwhere

    ! Microzooplantkon
    bPZ = (PZ_bins%zoo - 1) * M + 1
    ePZ = PZ_bins%zoo *M
    dPZdt(bPZ:ePZ) = Microzoo_eat - GrazMort_Z - NatMort_uzoo - Excr_uzoo

    ! Nitrate
    bPZ = (PZ_bins%NO - 1) * M + 1
    ePZ = PZ_bins%NO * M
    where (NO > 0.)
       dPZdt(bPZ:ePZ) = -uptake%NO + NH_oxid
    endwhere

    ! Ammonium
    bPZ = (PZ_bins%NH - 1) * M + 1
    ePZ = PZ_bins%NH * M
    where (NH > 0.)
       dPZdt(bPZ:ePZ) = -uptake%NH              &
               + was_NH &
               + remin_NH - NH_oxid
    endwhere

    ! Silicon
    bPZ = (PZ_bins%Si - 1) * M + 1
    ePZ = PZ_bins%Si * M
    where (Si > 0.)
       dPZdt(bPZ:ePZ) = -micro%growth%new * rate_micro%Si_ratio &
            - nano%growth%new * rate_nano%Si_ratio &
            - pico%growth%new * rate_pico%Si_ratio &
            + Si_remin
    endwhere

    ! Dissolved inorganic carbon
    bPZ = (PZ_bins%DIC - 1) * M + 1
    ePZ = PZ_bins%DIC * M
    where (DIC > 0.)
       dPZdt(bPZ:ePZ) = (remin_NH - uptake%NO - uptake%NH - uptake%PC) &
                  * Redfield_C
    endwhere

    ! Dissolved oxygen
    bPZ = (PZ_bins%Oxy - 1) * M + 1
    ePZ = PZ_bins%Oxy * M
    where (Oxy > 0.)
       dPZdt(bPZ:ePZ) = (uptake%NO + uptake%NH - remin_NH - NH_oxid) &
                  * Redfield_O
    endwhere

    ! Alkalinity
    bPZ = (PZ_bins%Alk - 1) * M + 1
    ePZ = PZ_bins%Alk * M
    where (Alk > 0.)
       dPZdt(bPZ:ePZ) = ((1 + Redfield_NP) * uptake%NO + (1 - Redfield_NP) &
            * (uptake%NH - remin_NH)) * (1 / Redfield_NP) - 2 * NH_oxid
    endwhere

    ! Dissolved organic carbon detritus
    bPZ = (PZ_bins%D_DOC - 1) * M + 1
    ePZ = PZ_bins%D_DOC * M
    where (D_DOC > 0.)
       dPZdt(bPZ:ePZ) = (was_DON &
                  - remin%D_DON * D_DON * temp_Q10) * Redfield_C
    endwhere
    ! Particulate organic carbon detritus
    bPZ = (PZ_bins%D_POC - 1) * M + 1
    ePZ = PZ_bins%D_POC * M
    where (D_POC > 0.)
       dPZdt(bPZ:ePZ) = (was_PON - Graz_PON &
                  - remin%D_PON * D_PON * temp_Q10) * Redfield_C
    endwhere

    ! Dissolved organic nitrogen detritus
    bPZ = (PZ_bins%D_DON - 1) * M + 1
    ePZ = PZ_bins%D_DON * M
    where (D_DON > 0.)
       dPZdt(bPZ:ePZ) = was_DON &
                  - remin%D_DON * D_DON * temp_Q10
    endwhere
    ! Particulate organic nitrogen detritus
    bPZ = (PZ_bins%D_PON - 1) * M + 1
    ePZ = PZ_bins%D_PON * M
    where (D_PON > 0.)
       dPZdt(bPZ:ePZ) = was_PON - Graz_PON &
                  - remin%D_PON * D_PON * temp_Q10
    endwhere
    ! Biogenic silicon detritus
    bPZ = (PZ_bins%D_bSi - 1) * M + 1
    ePZ = PZ_bins%D_bSi * M
    where (D_bSi > 0.)
       dPZdt(bPZ:ePZ) = was_Bsi &
               - Si_remin
    endwhere
  end subroutine derivs

end module NPZD
