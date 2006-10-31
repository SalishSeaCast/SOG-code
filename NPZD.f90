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
          sigma, & ! parameter in light limitation
          gamma, & ! loss parameter under light limitation
          inhib, & ! light inhibition
          k, & ! half saturation constant
          ! preference for NO3 over NH (always less than 1 as NH is preferred)
          kapa, &
          gamma_o, & ! exp strength of NH inhit of NO3 uptake
          N_o, & ! overall half saturation constant????
          N_x, & ! exponent in inhibition equation
          Si_ratio, & ! silicon to nitrogen ratio in phyto     
          Rm, & ! natural mortality rate
          M_z ! mortality rate
  end type rate_para_phyto
  !
  ! Parameters for mortality (where it goes)
  type :: loss_param
     real(kind=dp), dimension(:), pointer :: s,m  ! nano, micro
  end type loss_param
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
  type(rate_para_phyto) :: rate_micro, rate_nano
  type(loss_param) :: wastedestiny
  !
  ! Nitrogen compound uptake diagnotics
  type(uptake_) :: uptake
  !
  ! Nitrogen remineralization diagnostics
  real(kind=dp), dimension(:), pointer :: &
       remin_NH, &  ! Total remineralization to ammonium
       NH_oxid      ! Bacterial oxidation of NH4 to NO3

  integer :: D_bins

contains

  subroutine init_NPZD(M)
    ! Initialize biological model.
    ! *** Incomplete...
    use input_processor, only: getpard, getparl, getpari
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
    ! max growth rate for light limitation
    rate_micro%R = getpard('Micro, max growth')
    rate_nano%R = getpard('Nano, max growth')
    ! parameter in light limitation 
    rate_micro%sigma = getpard('Micro, sigma')
    rate_nano%sigma = getpard('Nano, sigma')
    ! loss parameter under light limitation
    rate_micro%gamma = getpard('Micro, gamma loss')
    rate_nano%gamma = getpard('Nano, gamma loss')
    ! light inhibition
    rate_micro%inhib = getpard('Micro, light inhib') 
    rate_nano%inhib = getpard('Nano, light inhib') 
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
    ! natural mortality rate
    rate_micro%Rm = getpard('Micro, nat mort')
    rate_nano%Rm = getpard('Nano, nat mort')
    ! mortality rate
    rate_micro%M_z = getpard('Micro, graze mort')
    rate_nano%M_z = getpard('Nano, graze mort')
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
    
    ! Detritus
    do xx= 1, D_bins
       binno = getpari("Bin No")
       if (binno .ne. xx) then
          write (*,*) "Bins must be in order, expecting"
          write (*,*) "bin number ",xx
          write (*,*) "but got bin number ",binno,"."
          stop
       endif
       wastedestiny%m(xx) = getpard('frac waste m')
       wastedestiny%s(xx) = getpard('frac waste s')
    enddo
    wastedestiny%m(0) = 1. - wastedestiny%m(1) - wastedestiny%m(2) - &
         wastedestiny%m(D_bins)
    wastedestiny%s(0) = 1. - wastedestiny%s(1) - wastedestiny%s(2) - &
         wastedestiny%s(D_bins)
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
    msg = "Waste destiny (biology model timestep initial conditions) array"
    allocate(wastedestiny%m(0:D_bins), wastedestiny%s(0:D_bins), &
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
    msg = "Waste destiny (biology model timestep initial conditions) array"
    deallocate(wastedestiny%m, wastedestiny%s, &
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
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    call dalloc_bio_RHS_variables()
  end subroutine dalloc_biology_variables


  subroutine p_growth(M, NO, NH, Si, P, I_par, Temp, rate, plank, & 
       NatMort, GrazMort) 
    ! Calculate the growth (light limited
    ! or nutrient limited) natural mortality and grazing mortality
    ! of either phytoplankton class.  All three
    ! are functions of temperature
    use precision_defs, only: dp
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
    real(kind=dp), dimension(0:M), intent(in) :: Temp  ! temperature
    ! parameters of the growth equations
    type(rate_para_phyto), intent(in) :: rate
    ! out are the growth values
    type(plankton2), intent(out) :: plank ! either micro or nano 
    real(kind=dp), dimension(1:M), intent(out) :: NatMort, GrazMort

    ! Local variables:
    integer :: j ! counter through depth
    real(kind=dp), dimension(1:M) :: Uc ! growth based on light
    real(kind=dp), dimension(1:M) :: Oup_cell ! NO uptake assuming full light
    real(kind=dp), dimension(1:M) :: Hup_cell ! NH uptake assuming full light
    real(kind=dp), dimension(1:M) :: Rmax ! maximum growth rate for given Temp
    real(kind=dp) :: temp_effect ! temperature limitation
    real(kind=dp) :: NH_effect ! ammonium effect on nutrient uptake

!!!!!!!!!!!Define growth Due to I_par!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    plank%growth%light = 0.

    Oup_cell = 0.
    Hup_cell = 0.

    ! flagella_jun06/ Next loop is different from V.16.1 => Rmax(replaces plank%R in KPP) 
    ! is temp dependent in SOG 

    DO j = 1,M

       ! biological process impactedby Q10 temperature effect
       temp_effect = 1.88**(0.1*(Temp(j)-273.15-20.)) 

       ! maximum growth rate (before light/nutrient limitation)
       Rmax(j)=rate%R*temp_effect 

       NatMort(j)=(rate%Rm)*temp_effect 
       GrazMort(j)=(rate%M_z)*temp_effect

       plank%growth%light(j) = Rmax(j)*(1.0-EXP(-rate%sigma*I_par(j)/Rmax(j)))

       Uc(j) = (1.0-rate%gamma)*plank%growth%light(j)* &
            (1/(1+rate%inhib*I_par(j)))
    END DO

!!!!!!!!!!!!!!!Define growth due to nutrients!!!!!!!!!!!!!!!!!!!!!!!!!

    DO j = 1,M

       NH_effect = exp(-rate%gamma_o * NH(j))
       IF (NO(j) > small) THEN
          Oup_cell(j) = Rmax(j) * NO(j) * rate%kapa * NH_effect / &
               (rate%k + NO(j) * rate%kapa * NH_effect + NH(j)) * &
               (NO(j) + NH(j))**rate%N_x / &
               (rate%N_o + (NO(j) + NH(j))**rate%N_x)
       ELSE
          Oup_cell(j) = 0.
       END IF
       IF (NH(j) > small) THEN
          Hup_cell(j) = Rmax(j) * NH(j) / &
               (rate%k + NO(j) * rate%kapa * NH_effect + NH(j))* &
               (NO(j) + NH(j))**rate%N_x / &
               (rate%N_o + (NO(j) + NH(j))**rate%N_x)
       ELSE
          Hup_cell(j) = 0.
       END IF

       plank%Nlimit(j) = (Oup_cell(j) + Hup_cell(j)) / Rmax(j)

       IF (Oup_cell(j) < 0.) THEN
          PRINT "(A)","Oup_cell(j) < 0. in reaction.f90"
          PRINT *,Oup_cell(j)
          STOP
       END IF

       IF (Hup_cell(j) < 0.) THEN
          PRINT "(A)","Hup_cell(j) < 0. in reaction.f90"
          PRINT *,Hup_cell(j)
          STOP
       END IF

    END DO

    ! Choose light limitation or nutrient limitation

    DO j = 1,M

       IF (Uc(j) >= 0. .AND. Uc(j) < Oup_cell(j) + Hup_cell(j)) THEN 

          !LIGHT LIMITING
          plank%growth%new(j) = Uc(j) 

          ! split the nutrient uptake between NH and NO

          IF (Uc(j) <= Hup_cell(j)) THEN
             ! add to nutrient uptake so we combined the effects of
             ! different phyto classes
             uptake%NH(j) = Uc(j) * P(j) + uptake%NH(j)
             uptake%NO(j) = uptake%NO(j)
          ELSE
             uptake%NH(j) = Hup_cell(j) * P(j) + uptake%NH(j)
             uptake%NO(j) = (Uc(j) - Hup_cell(j)) * P(j) + &
                  uptake%NO(j)
          END IF

       ELSE IF (Uc(j) >= 0. .AND. Uc(j) >= Oup_cell(j) + Hup_cell(j)) THEN

          !NUTRIENT LIMITING
          plank%growth%new(j) = Oup_cell(j) + Hup_cell(j)  

          IF (plank%growth%new(j) < 0.) THEN
             plank%growth%new(j) = 0.
          ENDIF

          uptake%NO(j) = Oup_cell(j) * P(j) + uptake%NO(j)
          uptake%NH(j) = Hup_cell(j) * P(j) + uptake%NH(j)

       ELSE  !No nutrient uptake, no growth
          plank%growth%new(j) =  0
       END IF

    END DO

  END SUBROUTINE p_growth

  subroutine derivs(M, PZ, Temp, I_par, day, dPZdt)
    ! Calculate the derivatives of the biology quantities for odeint()
    ! to use to calculate their values at the next time step.
    use precision_defs, only: dp
    use declarations, only: micro, nano, f_ratio
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
          WasteMicro,     &  ! Profile of micro phytos converted to waste
          WasteNano,      &  ! Profile of nano phytos converted to waste
          Si_remin           ! Profile of dissolution of biogenic Si detritus
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
    bPZ = (PZ_bins%nano - 1) * M + 1
    ePZ = PZ_bins%nano * M
    Pnano = PZ(bPZ:ePZ)
    where (Pnano < 0.) Pnano = 0.
    ! Nitrate
    bPZ = (PZ_bins%NO - 1) * M + 1
    ePZ = PZ_bins%NO * M
    NO = PZ(bPZ:ePZ)
    where (NO < 0.) NO = 0.
    ! Ammonium
    bPZ = (PZ_bins%NH - 1) * M + 1
    ePZ = PZ_bins%NH * M
    NH = PZ(bPZ:ePZ)
    where (NH < 0.) NH = 0.
    ! Silicon
    bPZ = (PZ_bins%Si - 1) * M + 1
    ePZ = PZ_bins%Si * M
    Si = PZ(bPZ:ePZ)
    where (Si < 0.) Si = 0.
    ! Dissolved organic nitrogen detritus
    bPZ = (PZ_bins%D_DON - 1) * M + 1
    ePZ = PZ_bins%D_DON * M
    D_DON = PZ(bPZ:ePZ)
    where (D_DON < 0.) D_DON = 0.
    ! Particulate organic nitrogen detritus
    bPZ = (PZ_bins%D_PON - 1) * M + 1
    ePZ = PZ_bins%D_PON * M
    D_PON = PZ(bPZ:ePZ)
    where (D_PON < 0.) D_PON = 0.
    ! Refractory nitrogen detritus
    bPZ = (PZ_bins%D_refr - 1) * M + 1
    ePZ = PZ_bins%D_refr * M
    D_refr = PZ(bPZ:ePZ)
    where (D_refr < 0.) D_refr = 0.
    ! Biogenic silicon detritus
    bPZ = (PZ_bins%D_bSi - 1) * M + 1
    ePZ = PZ_bins%D_bSi * M
    D_bSi = PZ(bPZ:ePZ)
    where (D_bSi < 0.) D_bSi = 0.

    ! Initialize transfer rates between the pools
    uptake%NO = 0.
    uptake%NH = 0.
    WasteMicro = 0.
    WasteNano = 0.
    remin_NH = 0.

    ! phytoplankton growth: Nitrate and Ammonimum, conc. of micro plankton
    ! I_par is light, Temp is temperature 
    ! N ammonium and nitrate uptake rates (IN and OUT) incremented
    ! micro is the growth parameters for micro plankton (IN) and the growth rates 
    ! (OUT)
    ! NatMort is physiological death, GrazMort is grazing mortality

    call p_growth(M, NO, NH, Si, Pmicro, I_par, Temp, & ! in
         rate_micro, micro, &         ! in and out, in, out
         NatMort_micro, GrazMort_micro)                    ! out

    ! put microplankton mortality into the medium detritus flux
!SEA    GrazMort_micro = (1.0 + 5.0 * exp(-(day-150.)**2/60.**2) + &
!SEA         2.0 * exp(-(day-230.)**2/60.**2) ) * GrazMort_micro
    ! copepods can't crop a full bloom (half-saturation 0.7 Alain)
!SEA    GrazMort_micro = GrazMort_micro / (2.0 + Pmicro)
    WasteMicro = WasteMicro + (GrazMort_micro+NatMort_micro)*Pmicro

    ! phytoplankton growth: Nitrate and Ammonimum, conc. of nano plankton
    ! I_par is light, Temp is temperature 
    ! N ammonium and nitrate uptake rates (IN and OUT) incremented
    ! nano is the growth parameters for micro plankton (IN) and the growth rates 
    ! (OUT)
    ! NatMort_nano is physiological death, GrazMort_nano is grazing mortality

    call p_growth(M, NO, NH, Si, Pnano, I_par, Temp, & ! in
         rate_nano, nano, &             ! in and out, in, out
         NatMort_nano, GrazMort_nano)                     ! out

!SEA    GrazMort_nano = GrazMort_nano*Pnano ! ie propto the square
    WasteNano = WasteNano + (GrazMort_nano+NatMort_nano)*Pnano

    ! Remineralization:
    !
    ! Bacterial oxidation of NH3 to NO pool proportional to NH^2
    NH_oxid = remin%NH * NH**2
    ! Dissolved organic nitrogen
    remin_NH = remin_NH + remin%D_DON * D_DON
    ! Particulate organic nitrogen
    remin_NH = remin_NH + remin%D_PON * D_PON
    ! Biogenic silcon
    Si_remin = remin%D_bSi * D_bSi

    IF (MINVAL(remin_NH) < 0.) THEN
       PRINT "(A)","remin_NH < 0. in derivs.f90"
       PRINT *,remin_NH
       STOP
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
       dPZdt(bPZ:ePZ) = (micro%growth%new - NatMort_micro &
               - GrazMort_micro) * Pmicro
    endwhere
    ! Nano phytoplankton
    bPZ = (PZ_bins%nano - 1) * M + 1
    ePZ = PZ_bins%nano * M
    where (Pnano > 0.) 
       dPZdt(bPZ:ePZ) = (nano%growth%new - NatMort_nano &
               - GrazMort_nano) * Pnano
    endwhere
    ! Microzooplantkon
    bPZ = (PZ_bins%zoo - 1) * M + 1
    ePZ = PZ_bins%zoo *M
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
               + WasteMicro * wastedestiny%m(0) &
               + WasteNano * wastedestiny%s(0)  &
               + remin_NH - NH_oxid
    endwhere
    ! Silicon
    bPZ = (PZ_bins%Si - 1) * M + 1
    ePZ = PZ_bins%Si * M
    where (Si > 0.) 
       dPZdt(bPZ:ePZ) = -micro%growth%new * Pmicro * rate_micro%Si_ratio &
            - nano%growth%new * Pnano * rate_nano%Si_ratio &
            + Si_remin
    endwhere
    ! Dissolved organic ditrogen detritus
    bPZ = (PZ_bins%D_DON - 1) * M + 1
    ePZ = PZ_bins%D_DON * M
    where (D_DON > 0.) 
       dPZdt(bPZ:ePZ) = WasteMicro * wastedestiny%m(1)  &
                  + WasteNano * wastedestiny%s(1) &
                  - remin%D_DON * D_DON
    endwhere
    ! Particulate organic nitrogen detritus
    bPZ = (PZ_bins%D_PON - 1) * M + 1
    ePZ = PZ_bins%D_PON * M
    where (D_PON > 0.) 
       dPZdt(bPZ:ePZ) = WasteMicro * wastedestiny%m(2)  &
                  + WasteNano * wastedestiny%s(2) &
                  - remin%D_PON * D_PON
    endwhere
    ! Biogenic silicon detritus
    bPZ = (PZ_bins%D_bSi - 1) * M + 1
    ePZ = PZ_bins%D_bSi * M
    where (D_bSi > 0.) 
       dPZdt(bPZ:ePZ) = (NatMort_micro + GrazMort_micro) &
               * Pmicro * rate_micro%Si_ratio &
               - Si_remin
    endwhere
  end subroutine derivs

end module NPZD
