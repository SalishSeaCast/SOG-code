! $Id$
! $Source$

module NPZD
  ! Type definitions, and subroutines for the biological model.
  !
  ! Public variables:
  !
  ! Public Variables:
  !
  !   flagellates -- Can flagellates can influence other biology?
  !
  !   remineralization -- Is there a remineralization loop?
  !
  !   microzooplankton -- are they eating things?
  !
  ! Public subroutines:
  !
  !   init_biology -- Initialize biology model.
  !
  !   derivs --
  !
  !   dalloc_biology_variables -- Deallocate memory from biology model
  !                               variables.

  use precision_defs, only: dp
  use declarations, only: M2 ! hopefully can get these out of everything else

  implicit none

  private
  public :: &
       ! Type definitions:
       bins, &
       ! Variables:
       PZ, &
       PZ_bins, &
       flagellates, &       ! Can flagellates can influence other biology?
       remineralization, &  ! Is there a remineralization loop?
       microzooplankton, &  ! Active or not
       ! diagnostics
       Mesozoo, &
       Ppico, &
       f_ratio, &  ! Ratio of new to total production profile
       ! Subroutines:
       init_NPZD, derivs, &
       dalloc_biology_variables

  ! Private type definitions:
  !
  ! Indices for quantities in PZ vector
  TYPE :: bins
     INTEGER :: &
          Quant, micro, nano, zoo, NO, NH, Si, det, D_DON, D_PON, D_refr, D_bSi
  END TYPE bins
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
          Rm, &! natural mortality rate
          winterconc, & ! winter concentration (pico only)
          summerconc  ! summer concentration (pico only)
  end type rate_para_phyto
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
          PredSlope, &       ! overall limit from predation
          HalfSat, &         ! overall half-saturation
          MicroPref, &       ! fractional preference for diatoms
          MicroPredSlope, &  ! limit from predation
          MicroHalfSat, &    ! half-saturation for getting eaten
          NanoPref, &       ! fractional preference for diatoms
          NanoPredSlope, &  ! limit from predation
          NanoHalfSat, &    ! half-saturation for getting eaten
          PicoHalfSat, &    ! half-saturation for getting eaten
          PON_Pref, &        ! fraction preference for PON
          PON_PredSlope, &   ! limit from predation
          PON_HalfSat        ! half-saturation for PON
  end type rate_para_zoo  

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
          NH     ! Ammonium uptake profile
  end type uptake_


  ! Public variable declarations:
  !
  ! Biological model logicals (turned parts of the model on or off)
  logical :: &
       flagellates, &    ! Can flagellates can influence other biology?
       remineralization, &  ! Is there a remineralization loop?
       microzooplankton ! use a microzooplantkon pool?
  real(kind=dp), dimension(:), allocatable :: Mesozoo, Ppico

  ! Private variable declarations:
  !
  ! Remineralization rates:
  type(remin_rates) :: remin
  !
  ! Indices for quantities (e.g. phyto, nitrate, etc.) in PZ vector
  type(bins):: PZ_bins
  data PZ_bins%Quant /6/    ! Size of the biology (Quantities and Detritus)
  data PZ_bins%micro /1/    ! Position of Diatoms (micro plankton)
  data PZ_bins%nano  /2/    ! Position of Flagellates (nano plankton)
  data PZ_bins%zoo   /3/    ! Position of microzooplankton
  data PZ_bins%NO    /4/    ! Position of Nitrate
  data PZ_bins%NH    /5/    ! Position of Ammonium
  data PZ_bins%Si    /6/    ! Position of Silicon
  data PZ_bins%det   /7/    ! Start of detritus
  data PZ_bins%D_DON   /7/  ! Position of dissolved organic nitrogen detritus
  data PZ_bins%D_PON   /8/  ! Position of particulate organic nitrogen detritus
  data PZ_bins%D_refr  /10/  ! Position of refractory nitrogen detritus
  data PZ_bins%D_bSi   /9/  ! Position of dissolved biogenic silcon detritus
  !
  ! PZ vector
  real(kind=dp), dimension(:), allocatable :: PZ
  !
  ! Biological rate parameters
  type(rate_para_phyto) :: rate_micro, rate_nano, rate_pico
  type(rate_para_zoo) :: rate_mesozoo, rate_mesorub
  type(nloss_param) :: frac_waste_DNM, & ! Diatom natural mortality
       frac_waste_NNM, & ! Nano natural mortality
       frac_waste_DEM, & ! Diatoms eaten by mesozoo
       frac_waste_NEM, & ! Nano eaten by mesozoo
       frac_waste_PEM   ! PON eaten by Mesozoo
  !
  ! Nitrogen compound uptake diagnotics
  type(uptake_) :: uptake
  !
  ! Nitrogen remineralization diagnostics
  real(kind=dp), dimension(:), pointer :: &
       remin_NH, &  ! Total remineralization to ammonium
       NH_oxid      ! Bacterial oxidation of NH4 to NO3
  real(kind=dp), dimension(:), allocatable :: &
       f_ratio  ! Ratio of new to total production profile
  integer :: D_bins

contains

  subroutine init_NPZD(M)
    ! Initialize biological model.
    ! *** Incomplete...
    use input_processor, only: getpard, getparl, getpari, getpardv
    use biology_eqn_builder, only: alloc_bio_RHS_variables

    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables
    integer :: xx ! counter through detritus bins
    integer :: binno ! detritus bin as read

    ! Number of detritus bins, dissolved, slow sink and fast sink
    D_bins = 4

    ! Size of the PZ vector for biological model
    M2 = (PZ_bins%Quant + D_bins) * M   !size of PZ in biology: 
    ! Allocate memory for biology model variables
    call alloc_biology_variables(M)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    call alloc_bio_RHS_variables(M)

    ! Read the values of the flagellates and remineralization loop
    ! selector flags
    flagellates = getparl('flagellates_on')
    remineralization = getparl('remineralization')
    microzooplankton = getparl('use microzooplankton')


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
    ! preference for PON
    rate_mesozoo%PON_Pref = getpard('Mesozoo, pref for PON')
    ! limit from predation
    rate_mesozoo%PON_PredSlope = getpard('Mesozoo, PON pred slope')
    ! half saturation
    rate_mesozoo%PON_HalfSat = getpard('Mesozoo, PON half-sat')
    ! Mesodinium rates
    rate_mesorub%R = getpard("Mesorub, max ingestion")
    rate_mesorub%eff = getpard("Mesorub, assimilation eff")
    rate_mesorub%PicoHalfSat = getpard('Mesorub, nano half-sat')

    ! phytoplankton rates
    rate_pico%winterconc = getpard('Pico, winter conc')
    rate_pico%summerconc = getpard('Pico, summer conc')
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
    ! temp range of limitation below max
    rate_micro%temprange = getpard('Micro, temp range')
    rate_nano%temprange = getpard('Nano, temp range')
    ! exponent for Q_10
    rate_micro%Q10exp = getpard('Micro, Q10 exp')
    rate_nano%Q10exp = getpard('Nano, Q10 exp')
    ! loss parameter under light limitation
    rate_micro%gamma = getpard('Micro, gamma loss')
    rate_nano%gamma = getpard('Nano, gamma loss')
    ! NO3 half saturation constant
    rate_micro%k = getpard('Micro, NO3 k')
    rate_nano%k = getpard('Nano, NO3 k')
    ! preference for NO3 over NH 
    rate_micro%kapa = getpard('Micro, kapa')
    rate_nano%kapa = getpard('Nano, kapa')
    ! exp strength of NH inhit of NO3 uptake
    rate_micro%gamma_o = getpard('Micro, NH inhib')
    rate_nano%gamma_o = getpard('Nano, NH inhib')
    ! overall half saturation constant????
    rate_micro%N_o = getpard('Micro, N_o')
    rate_nano%N_o = getpard('Nano, N_o')
    ! exponent in inhibition equation
    rate_micro%N_x = getpard('Micro, N_x')
    rate_nano%N_x = getpard('Nano, N_x')
    ! silicon to nitrogen ratio in phyto
    rate_micro%Si_ratio = getpard('Micro, Si ratio')
    rate_nano%Si_ratio = getpard('Nano, Si ratio')
    ! half-saturation for Silicon
    rate_micro%K_Si = getpard('Micro, K Si')
    if (rate_micro%Si_ratio.eq.0) rate_micro%K_Si = 0. ! to eliminate Si limitn
    rate_nano%K_Si = getpard('Nano, K Si')
    if (rate_nano%Si_ratio.eq.0) rate_nano%K_Si = 0. ! to eliminate Si limitn
    ! natural mortality rate
    rate_micro%Rm = getpard('Micro, nat mort')
    rate_nano%Rm = getpard('Nano, nat mort')
    ! Remineralization rates:
    ! Ammonium
    remin%NH = getpard('NH remin rate')
    ! DON detritus
    remin%D_DON = getpard('DON remin rate')
    ! PON detritus
    remin%D_PON = getpard('PON remin rate')
    ! Biogenic silicon detritus
    remin%D_bSi = getpard('bSi remin rate')
    ! If remineralization is turned off, set the rates to zero
    if (.not. remineralization) then
       remin%NH = 0.
       remin%D_DON = 0.
       remin%D_PON = 0.
       remin%D_bSi = 0.
    endif
    
    ! New waste code
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


    ! PON eaten by Mesozoo
    frac_waste_PEM%NH = getpard('Waste, pem, NH')
    frac_waste_PEM%DON = getpard('Waste, pem, DON')
    frac_waste_PEM%PON = getpard('Waste, pem, PON')
    frac_waste_PEM%Ref = getpard('Waste, pem, Ref')
    frac_waste_PEM%Bsi = getpard('Waste, pem, Bsi')

  end subroutine init_NPZD


  subroutine alloc_biology_variables(M)
    ! Allocate memory for biological model arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "PZ vector (biology model timestep initial conditions) array"
    allocate(PZ(1:(PZ_bins%Quant + D_bins) * M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrogen compounds uptake diagnostic arrays"
    allocate(uptake%NO(1:M), uptake%NH(1:M), &
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
    msg = "Pico phytoplankton diagnostic array"
    allocate(Ppico(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ratio of new to total production profile array"
    allocate(f_ratio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_biology_variables


  subroutine dalloc_biology_variables
    ! Deallocate memory for biological model arrays.
    use malloc, only: dalloc_check
    use biology_eqn_builder, only: dalloc_bio_RHS_variables
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "PZ vector (biology model timestep initial conditions) array"
    deallocate(PZ, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrogen compounds uptake diagnostic arrays"
    deallocate(uptake%NO, uptake%NH, &
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
    msg = "Pico phytoplankton diagnostic array"
    deallocate(Ppico, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ratio of new to total production profile array"
    deallocate(f_ratio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    call dalloc_bio_RHS_variables()
  end subroutine dalloc_biology_variables


  subroutine p_growth(M, NO, NH, Si, P, I_par, temp_Q10, Temp, rate, plank) 
    ! Calculate the growth (light limited
    ! or nutrient limited) 
    ! of either phytoplankton class which 
    ! are functions of temperature
    use precision_defs, only: dp
    use unit_conversions, only: KtoC
    use mean_param, only: plankton2
    use surface_forcing, only: small
    implicit none
    ! Arguments:
    integer, intent(in) :: M
    ! Nitrate, ammonium & silicon concentrations
    real(kind=dp), dimension(1:M), intent(in) :: NO, NH, Si 
    ! plankton concentraton (either Pmicro or Pnano)
    real(kind=dp), dimension(1:M), intent(in) :: P
    real(kind=dp), dimension(0:M), intent(in) :: I_par  ! light
    real(kind=dp), dimension(1:M), intent(in) :: temp_Q10  ! Q10 temp effect
    real(kind=dp), dimension(0:), intent(in) :: temp ! temperature
    ! parameters of the growth equations
    type(rate_para_phyto), intent(in) :: rate
    ! out are the growth values
    type(plankton2), intent(out) :: plank ! either micro or nano 

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
            rate%temprange

       ! don't include Rmax effect until end

       plank%growth%light(j) = &
            ! Steeles scheme like but with extended high light range
            ! (Steeles scheme has built in light inhibition)
            ! 0.67, 2.7 and 1.8 are constants for making the fit
            ! Steeles like at small light and making it fit Durbin for
            ! Thalassosira nordelenski at higher light
            (1.0d0 - exp(-I_par(j) / (0.67d0 * rate%Iopt)) ) * &
            exp(-I_par(j) / (2.7d0 * rate%Iopt)) * 1.8d0

       Uc(j) = (1.0d0 - rate%gamma) * plank%growth%light(j)
    END DO

!!!!!!!!!!!!!!!Define growth due to nutrients!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Silicon (vectorized)

    Sc = Si / (rate%K_Si + Si)

    ! Nitrate and Ammonium

    DO j = 1,M

       NH_effect = exp(-rate%gamma_o * NH(j))
       IF (NO(j) > small) THEN
          Oup_cell(j) = NO(j) * rate%kapa * NH_effect / &
               (rate%k + NO(j) * rate%kapa * NH_effect + NH(j)) * &
               (NO(j) + NH(j))**rate%N_x / &
               (rate%N_o + (NO(j) + NH(j))**rate%N_x)
       ELSE
          Oup_cell(j) = 0.
       END IF
       IF (NH(j) > small) THEN
          Hup_cell(j) = NH(j) / &
               (rate%k + NO(j) * rate%kapa * NH_effect + NH(j))* &
               (NO(j) + NH(j))**rate%N_x / &
               (rate%N_o + (NO(j) + NH(j))**rate%N_x)
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

       endif
    END DO
 
  END SUBROUTINE p_growth

  subroutine derivs(M, PZ, Temp, I_par, day, dPZdt)
    ! Calculate the derivatives of the biology quantities for odeint()
    ! to use to calculate their values at the next time step.
    use precision_defs, only: dp
    use fundamental_constants, only: pi
    use declarations, only: micro, nano
    use unit_conversions, only: KtoC
    use grid_mod, only: depth_average
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(1:), intent(in) :: &
         PZ  ! Consolidated array of biology quantities values
    real(kind=dp), dimension(0:M), intent(in) :: &
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
          Pnano,          &  ! Nano phytoplankton (flagellates) profile
          Z,              &  ! Micro zooplankton profile
          NO,             &  ! Nitrate concentrationy
          NH,             &  ! Ammonium concentration profile
          Si,             &  ! Silicon concentration profile
          D_DON,          &  ! Dissolved organic nitrogen detritus profile
          D_PON,          &  ! Particulate organic nitrogen detritus profile
          D_refr,         &  ! Refractory nitrogen detritus profile
          D_bSi,          &  ! Biogenic silicon detritus profile
          NatMort_micro,  &  ! Micro phytoplankton natural mortality profile
          GrazMort_micro, &  ! Micro phytoplankton grazing mortality profile
          NatMort_nano,   &  ! Nano phytoplankton natural mortality profile
          GrazMort_nano,  &  ! Nano phytoplankton grazing mortality profile
          Graz_PON,       &  ! PON grazing mortality profile
          Pico_growth,    &  ! maximum growth of pico
          Mesorub_eat,    &  ! Mesorub eat pico
! new waste variables
          was_NH, was_DON, was_PON, was_Ref, was_Bsi, &
          Si_remin           ! Profile of dissolution of biogenic Si detritus
    real(kind=dp), dimension(1:M) :: temp_Q10
    ! temporary variable to calculate mortality scheme
    real(kind=dp) :: food_limitation, denominator, &
         Meso_mort_micro, Meso_graz_PON , Meso_mort_nano, average_prey
   
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

    ! change growth to be total
    micro%growth%new = micro%growth%new * Pmicro
    nano%growth%new = nano%growth%new * Pnano

    ! phytoplankton natural mortality:

    NatMort_micro= (rate_micro%Rm) * temp_Q10 * Pmicro
    NatMort_nano = (rate_nano%Rm) * temp_Q10 * Pnano

    was_NH = was_NH + frac_waste_DNM%NH * NatMort_micro &
         + frac_waste_NNM%NH * NatMort_nano
    was_DON = was_DON + frac_waste_DNM%DON * NatMort_micro &
         + frac_waste_NNM%DON * NatMort_nano
    was_PON = was_PON + frac_waste_DNM%PON * NatMort_micro &
         + frac_waste_NNM%PON * NatMort_nano
    was_Ref = was_Ref + frac_waste_DNM%Ref * NatMort_micro &
         + frac_waste_NNM%Ref * NatMort_nano
    was_BSi = was_BSi &
         + frac_waste_DNM%BSi * NatMort_micro * rate_micro%Si_ratio &
         + frac_waste_NNM%Bsi * NatMort_nano * rate_nano%Si_ratio

    ! Grazing processes

    ! Mesozooplankton
    ! amount of Mesozoo : const. winter conc + 3 summer Gaussians
    Mesozoo(1:M) = rate_mesozoo%winterconc + & 
         rate_mesozoo%summerconc * &
         sum ( rate_mesozoo%sumpeakval * &
           exp( -(day-rate_mesozoo%sumpeakpos)**2 / &
                                rate_mesozoo%sumpeakwid**2 ) )

    average_prey = depth_average(Pmicro+D_PON+Pnano,0.25d0,37.75d0)

    Mesozoo(1:M) = Mesozoo(1:M) * (Pmicro(1:M) + D_PON(1:M) + Pnano(1:M)) &
         / average_prey

    do jj = 1,M
       ! global food limitation
       food_limitation = (Pmicro(jj) + D_PON(jj) + Pnano(jj) &
            - rate_mesozoo%PredSlope) / &
            (rate_mesozoo%HalfSat + Pmicro(jj) + D_PON(jj) + Pnano(jj) &
            - rate_mesozoo%PredSlope + epsilon(rate_mesozoo%HalfSat))
       
       denominator = (rate_mesozoo%MicroPref * Pmicro(jj) + &
            rate_mesozoo%NanoPref * Pnano(jj) + &
            rate_mesozoo%PON_Pref * D_PON(jj) + epsilon(Pmicro(jj)) )

       ! limitation based on microplankton
       Meso_mort_micro = min(rate_mesozoo%MicroPref * food_limitation &
            * Pmicro(jj) / denominator, &
            (Pmicro(jj) - rate_mesozoo%MicroPredslope) / &
            (rate_mesozoo%MicroHalfSat + Pmicro(jj) &
            - rate_mesozoo%MicroPredSlope + &
            epsilon(rate_mesozoo%MicroHalfSat)))

       ! limitation based on nanoplankton
       Meso_mort_nano = min(rate_mesozoo%NanoPref * food_limitation &
            * Pnano(jj) / denominator, &
            (Pnano(jj) - rate_mesozoo%NanoPredslope) / &
            (rate_mesozoo%NanoHalfSat + Pnano(jj) &
            - rate_mesozoo%NanoPredSlope + &
            epsilon(rate_mesozoo%NanoHalfSat)))

       ! limitation based on PON
       Meso_graz_PON = min(rate_mesozoo%PON_Pref * food_limitation &
            * D_PON(jj) / denominator, &
            (D_PON(jj) - rate_mesozoo%PON_Predslope) / &
            (rate_mesozoo%PON_HalfSat + D_PON(jj) - rate_mesozoo%PON_PredSlope &
            + epsilon(rate_mesozoo%PON_HalfSat)))

       ! global corrected by individual
       food_limitation = Meso_mort_micro + Meso_mort_nano + Meso_graz_PON

       GrazMort_micro(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * & 
            max(0., Meso_mort_micro)

       GrazMort_nano(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * &
            max(0., Meso_mort_nano)

       Graz_PON(jj) = rate_mesozoo%R * temp_Q10(jj) * Mesozoo(jj) * &
            max(0., Meso_graz_PON)
       
    enddo

    was_NH = was_NH + frac_waste_PEM%NH * Graz_PON + &
         frac_waste_DEM%NH * GrazMort_micro + &
         frac_waste_NEM%NH * GrazMort_nano
    was_DON = was_DON + frac_waste_PEM%DON * Graz_PON + &
         frac_waste_DEM%DON * GrazMort_micro + &
         frac_waste_NEM%DON * GrazMort_nano
    was_PON = was_PON + frac_waste_PEM%PON * Graz_PON + &
         frac_waste_DEM%PON * GrazMort_micro + &
         frac_waste_NEM%PON * GrazMort_nano
    was_Ref = was_Ref + frac_waste_PEM%Ref * Graz_PON + &
         frac_waste_DEM%Ref * GrazMort_micro + &
         frac_waste_NEM%Ref * GrazMort_nano
    was_BSi = was_BSi + frac_waste_PEM%BSi * Graz_PON * 0.d0 + &
         frac_waste_DEM%BSi * GrazMort_micro * rate_micro%Si_ratio + &
         frac_waste_NEM%BSi * GrazMort_nano * rate_nano%Si_ratio

    ! pico phytoplankton

    Ppico = rate_pico%winterconc + & 
         (rate_pico%summerconc - rate_pico%winterconc) * &
         (1.d0 - cos(2.d0*pi*day/365.25))/2.d0
 
    Pico_growth = &
            ! Steeles scheme like but with extended high light range
            ! (Steeles scheme has built in light inhibition)
            ! 0.67, 2.7 and 1.8 are constants for making the fit
            ! Steeles like at small light and making it fit Durbin for
            ! Thalassosira nordelenski at higher light
         (1.0d0 - exp(-I_par / (0.67d0 * rate_pico%Iopt))) * &
         exp(-I_par / (2.7d0 * rate_pico%Iopt)) * 1.8d0

    Pico_growth = rate_pico%R * Ppico * Pico_growth * temp_Q10

!   currently (because of p_growth limit above) only allowing them to eat 
! during the day

    Mesorub_eat = 2.d0 & ! to make up for eating day and night
         * rate_mesorub%R * Ppico / (rate_mesorub%PicoHalfSat + Ppico) &
         * Pnano * temp_Q10

    ! Mesorub get either max they want to eat or max prod of pico
    do jj=1,M
       Mesorub_eat(jj) = min(Pico_growth(jj), Mesorub_eat(jj)) * &
            rate_mesorub%eff
    enddo

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
    ! Microzooplantkon
    bPZ = (PZ_bins%zoo - 1) * M + 1
    ePZ = PZ_bins%zoo *M
!    dPZdt(bPZ:ePZ) = GrazMort_nano * Pnano - 0.1 * Z
    dPZdt(bPZ:ePZ) = 0.
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
            + Si_remin
    endwhere
    ! Dissolved organic ditrogen detritus
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
