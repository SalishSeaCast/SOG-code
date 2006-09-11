! $Id$
! $Source$

module biological_mod
  ! Type definitions, and subroutines for the biological model.

  use precision_defs, only: dp
  use declarations, only: D_bins, M2 ! hopefully can get these out of everything else

  implicit none

  private

  public :: derivs_sog, &
       define_PZ, & ! temporary
       reaction_p_sog, &
       init_biology

  ! Private type definitions:
  !
  ! Indices for quantities in PZ vector
  TYPE :: bins
     INTEGER :: micro, nano, NO, NH, Sil, det, Quant
  END TYPE bins

  !
  ! Rate parameters for phytoplankton
  type :: rate_para_phyto
     real(kind=dp) :: R, &  ! max growth rate for light limitation
          sigma, & ! parameter in light limitation
          gamma, & ! loss parameter under light limitation
          inhib, & ! light inhibition
          k, & ! NO3 half saturation constant
          ! preference for NO3 over NH (always less than 1 as NH is preferred)
          kapa, & 
          gamma_o, & ! exp strength of NH inhit of NO3 uptake
          N_o, & ! overall half saturation constant????
          N_x, & ! exponent in inhibition equation
          Sil_rat, & ! silicon to nitrogen ratio in phyto     
          Rm, & ! respiration rate
          M_z ! mortality rate
  end type rate_para_phyto

  ! rate parameters for nutrients !*** temporary for old N%r which got orphaned
  type :: rate_para_nutrient
     real(kind=dp) r ! reminerialization rate NH to NO3
  end type rate_para_nutrient
  
  ! Private variable declarations:
  !
  ! Indices for quantities (e.g. phyto, nitrate, etc.) in PZ vector
  type(bins):: PZ_bins
  data PZ_bins%Quant /5/  ! Size of the biology (Quantities and Detritus)
  data PZ_bins%micro /1/  ! Position of Diatoms (micro plankton)
  data PZ_bins%nano  /2/  ! Position of Flagellates (nano plankton)
  data PZ_bins%NO    /3/  ! Position of Nitrate
  data PZ_bins%NH    /4/  ! Position of Ammonium
  data PZ_bins%Sil   /5/  ! Position of Silicon
  data PZ_bins%det   /6/  ! Start of detritus
  !
  ! PZ vector
  real(kind=dp), dimension(:), allocatable :: PZ
  ! Biological model logicals (turned parts of the model on or off)
  logical :: getparl ! until getpar is a module
  real(kind=dp) :: getpard ! until getpar is a module
  external getparl, getpard
  logical :: flagellates  ! whether flagellates can influence other biology
  logical :: remineralization !whether there is a remineralization loop
  ! Biological rate parameters
  type(rate_para_phyto) :: rate_micro, rate_nano
  type(rate_para_nutrient) :: rate_N

contains

  subroutine alloc_biology
    ! Allocate memory for biological model arrays.
    use malloc, only: alloc_check
    implicit none
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "PZ vector (biology model timestep initial conditions) array"
    allocate(PZ(M2), stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_biology


  subroutine init_biology (M)
    ! Initialize biological model.
    ! *** Incomplete...
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points

    ! Number of detritus bins, dissolved, slow sink and fast sink
    D_bins = 3

    ! Size of the PZ vector for biological model
    M2 = (PZ_bins%Quant + D_bins) * M   !size of PZ in biology: 
    ! Allocate memory for PZ
    call alloc_biology

    flagellates = getparl('flagellates_on',1)
    remineralization = getparl('remineralization',1)

    ! Biological rate parameters
    ! max growth rate for light limitation
    rate_micro%R = getpard ('Micro, max growth',1)
    rate_nano%R = getpard ('Nano, max growth',1)
    ! parameter in light limitation 
    rate_micro%sigma = getpard ('Micro, sigma', 1)
    rate_nano%sigma = getpard ('Nano, sigma', 1)
    ! loss parameter under light limitation
    rate_micro%gamma = getpard ('Micro, gamma loss',1)
    rate_nano%gamma = getpard ('Nano, gamma loss',1)
    ! light inhibition
    rate_micro%inhib = getpard('Micro, light inhib',1) 
    rate_nano%inhib = getpard('Nano, light inhib',1) 
    ! NO3 half saturation constant
    rate_micro%k = getpard('Micro, NO3 k',1)
    rate_nano%k = getpard('Nano, NO3 k',1)
    ! preference for NO3 over NH 
    rate_micro%kapa = getpard('Micro, kapa',1)
    rate_nano%kapa = getpard('Nano, kapa',1)
    ! exp strength of NH inhit of NO3 uptake
    rate_micro%gamma_o = getpard('Micro, NH inhib',1)
    rate_nano%gamma_o = getpard('Nano, NH inhib',1)
    ! overall half saturation constant????
    rate_micro%N_o = getpard('Micro, N_o',1)
    rate_nano%N_o = getpard('Nano, N_o',1)
    ! exponent in inhibition equation
    rate_micro%N_x = getpard('Micro, N_x',1)
    rate_nano%N_x = getpard('Nano, N_x',1)
    ! silicon to nitrogen ratio in phyto
    rate_micro%Sil_rat = getpard('Micro, Sil ratio',1)
    rate_nano%Sil_rat = getpard('Nano, Sil ratio',1)
    ! respiration rate
    rate_micro%Rm = getpard('Micro, resp',1)
    rate_nano%Rm = getpard('Nano, resp',1)
    ! mortality rate
    rate_micro%M_z = getpard('Micro, mort',1)
    rate_nano%M_z = getpard('Nano, mort',1)
    ! nitrate remineralization rate
    rate_N%r = getpard('Nitrate, remineral',1)

  end subroutine init_biology


  SUBROUTINE define_PZ(M, Pmicro, Pnano, NO, NH, Sil, Detritus, &
       PZ)

    ! This subroutine takes all the separate variables (microplankton,
    ! nanoplankton, nitrate, ammonium and detritus and loads them
    ! sequentially into the PZ vector for the ODE solver to use.

    USE mean_param, only : snow 


    IMPLICIT NONE

    integer, intent(in):: M
    DOUBLE PRECISION, DIMENSION(0:), INTENT(IN)::Pmicro, Pnano, NO, NH, Sil
    TYPE(snow), DIMENSION(D_bins), INTENT(IN)::Detritus
    DOUBLE PRECISION, DIMENSION (M2), INTENT (OUT) :: PZ

    ! Local Variables

    INTEGER :: bPZ, ePZ ! start position and end position in PZ array
    INTEGER :: j

    PZ = 0.

    ! Microplankton

    bPz = (PZ_bins%micro-1) * M + 1
    ePZ = PZ_bins%micro * M
    PZ(bPZ:ePZ) = Pmicro(1:M)

    ! Nanoplankton

    bPz = (PZ_bins%nano-1) * M + 1
    ePZ = PZ_bins%nano * M
    if (flagellates) then
       PZ(bPZ:ePZ) = Pnano(1:M)
    else
       PZ(bPZ:ePZ) = 0
    endif

    ! Nitrate

    bPz = (PZ_bins%NO-1) * M + 1
    ePZ = PZ_bins%NO * M
    PZ(bPZ:ePZ) = NO(1:M)

    ! Ammonium

    bPz = (PZ_bins%NH-1) * M + 1
    ePZ = PZ_bins%NH * M
    if (remineralization) then
       PZ(bPZ:ePZ) = NH(1:M)
    else
       PZ(bPZ:ePZ) = 0
    endif

    ! Silicon

    bPz = (PZ_bins%Sil-1) * M + 1
    ePZ = PZ_bins%Sil * M
    PZ(bPZ:ePZ) = Sil(1:M)

    ! Detritus
    do j=1,D_bins
       bPz = (PZ_bins%Quant + (j-1) ) * M + 1
       ePz = (PZ_bins%Quant + j) * M
       if (remineralization) then
          PZ(bPz:ePz) = Detritus(j)%D%new(1:M)
       else
          PZ(bPz:ePz) = 0.
       endif
    enddo

  END SUBROUTINE define_PZ

  SUBROUTINE reaction_p_sog(M, PZ, Pmicro, Pnano, NO, &
       NH, Sil, Detritus, Gvector)

    USE mean_param, only: snow, UVST

    IMPLICIT NONE

    integer, intent(in):: M
    DOUBLE PRECISION, DIMENSION(M2), INTENT(IN)::PZ
    DOUBLE PRECISION, DIMENSION(0:M+1), INTENT(IN)::Pmicro,Pnano, NO, NH, Sil
    TYPE(snow), DIMENSION(D_bins), INTENT(IN)::Detritus
    TYPE(UVST), INTENT(IN OUT)::Gvector  ! IN only 'cause not setting all

    ! Local variables
    INTEGER :: bPZ, ePZ ! start position and end position in PZ array
    INTEGER::j

    ! Micro plankton
    bPZ = (PZ_bins%micro-1) * M + 1
    ePZ = PZ_bins%micro * M
    Gvector%p%micro = PZ(bPZ:ePZ) - Pmicro(1:M)

    ! Nano plankton
    bPZ = (PZ_bins%nano-1) * M + 1
    ePZ = PZ_bins%nano * M
    Gvector%p%nano = PZ(bPZ:ePZ) - Pnano(1:M)

    ! Nitrate
    bPZ = (PZ_bins%NO-1) * M + 1
    ePZ = PZ_bins%NO * M
    Gvector%N%O = PZ(bPZ:ePZ) - NO(1:M)

    ! Ammonimum
    bPZ = (PZ_bins%NH-1) * M + 1
    ePZ = PZ_bins%NH * M
    Gvector%N%H = PZ(bPZ:ePZ) - NH(1:M)

    ! Silicon
    bPZ = (PZ_bins%Sil-1) * M + 1
    ePZ = PZ_bins%Sil * M
    Gvector%Sil = PZ(bPZ:ePZ) - Sil(1:M)

    ! Detritus
    DO j = 1,D_bins
       bPz = (PZ_bins%Quant + (j-1) ) * M + 1
       ePz = (PZ_bins%Quant + j) * M
       Gvector%d(j)%bin = PZ(bPz:ePz) - Detritus(j)%D%old(1:M)
    END DO

  END SUBROUTINE reaction_p_sog

  SUBROUTINE p_growth(M, NO, NH, Sil, P, I_par, Temp, N, rate, plank, & 
       Resp, Mort) 

    ! this subroutine calculates the growth, respiration and mortality
    ! (light limited or nutrient limited) of either phytoplankton class.
    ! All three are functions of temperature

    use mean_param, only: plankton2, nutrient
    use surface_forcing, only: small

    implicit none

    integer, intent(in) :: M
    ! Nitrate and ammonium concentrations
    double precision, dimension(M), intent(in) :: NO, NH, Sil 
    ! plankton concentraton (either Pmicro or Pnano)
    DOUBLE PRECISION, DIMENSION(M), INTENT(IN)::P !V.flagella.01 note: either Pmicro or Pnano

    DOUBLE PRECISION, DIMENSION(0:M), INTENT(IN):: I_par  ! light
    DOUBLE PRECISION, DIMENSION(0:M), INTENT(IN):: Temp  ! temperature


    ! parameters of the growth equations
    type(rate_para_phyto), intent(in) :: rate
    ! out are the growth values
    TYPE(plankton2), INTENT(OUT):: plank ! either micro or nano 

    ! nutrient uptake values (incremented)
    TYPE(nutrient), INTENT(IN OUT):: N

    DOUBLE PRECISION, DIMENSION(M), intent(out) :: Resp, Mort


    ! internal variables

    INTEGER::j ! counter through depth

    DOUBLE PRECISION, DIMENSION(M) :: Uc ! growth based on light
    double precision, dimension(M) :: Oup_cell ! NO uptake assuming full light
    double precision, dimension(M) :: Hup_cell ! NH uptake assuming full light
    DOUBLE PRECISION, DIMENSION(M) :: Rmax ! maximum growth rate for given Temp
    double precision :: temp_effect ! temperature limitation
    double precision :: NH_effect ! ammonium effect on nutrient uptake

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

       Resp(j)=(rate%Rm)*temp_effect 
       Mort(j)=(rate%M_z)*temp_effect

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
             N%H_uptake%new(j) = Uc(j) * P(j) + N%H_uptake%new(j)
             N%O_uptake%new(j) = 0.
          ELSE
             N%H_uptake%new(j) = Hup_cell(j) * P(j) + N%H_uptake%new(j)
             N%O_uptake%new(j) = (Uc(j) - Hup_cell(j)) * P(j) + &
                  N%O_uptake%new(j)
          END IF

       ELSE IF (Uc(j) >= 0. .AND. Uc(j) >= Oup_cell(j) + Hup_cell(j)) THEN

          !NUTRIENT LIMITING
          plank%growth%new(j) = Oup_cell(j) + Hup_cell(j)  

          IF (plank%growth%new(j) < 0.) THEN
             plank%growth%new(j) = 0.
          ENDIF

          N%O_uptake%new(j) = Oup_cell(j) * P(j) + N%O_uptake%new(j)
          N%H_uptake%new(j) = Hup_cell(j) * P(j) + N%H_uptake%new(j)

       ELSE  !No nutrient uptake, no growth
          plank%growth%new(j) =  0
       END IF

    END DO

  END SUBROUTINE p_growth

  subroutine derivs_sog(M, PZ, dPZdt, Temp, I_par)
    ! Calculate the derivatives of the biological model for odeint to
    ! use to advance the biology to the next time step.

    use precision_defs, only: dp
    use mean_param, only: losses, plankton2, nutrient, snow
    use declarations, only: D_bins, N, micro, nano, &
         f_ratio, waste, Detritus

    implicit none

    ! Arguments:
    integer, intent(in) :: M ! grid size
    real(kind=dp), DIMENSION(M2), INTENT(IN):: PZ  ! values
    real(kind=dp), DIMENSION(M2), INTENT(OUT)::dPZdt ! derivatives
    real(kind=dp), dimension(0:M), INTENT(IN):: Temp ! temperature
    real(kind=dp), dimension(0:M), INTENT(IN):: I_par ! light

    ! Local variables

    integer :: ii ! counter through grid ie 1 to M
    integer :: jj ! counter through PZ
    integer :: kk ! counter through detritus

    real(kind=dp), DIMENSION(M):: Resp_micro, Mort_micro! respiration & mortality
    real(kind=dp), DIMENSION(M):: Resp_nano, Mort_nano  ! respiration & mortality

    real(kind=dp), dimension (M) :: Pmicro, Pnano ! micro/nano plankton conc.
    real(kind=dp), dimension (M) :: NO, NH,Sil ! nitrate and ammonium conc.
    real(kind=dp), dimension (D_bins, M) :: detr ! detritus

    ! Put PZ micro values into Pmicro variable, removing any negative values
    do ii = 1,M                        ! counter through grid
       jj = (PZ_bins%micro-1) * M + ii ! counter into PZ

       if (PZ(jj)>0) then
          Pmicro(ii) = PZ(jj)
       else
          Pmicro(ii) = 0
       endif
    enddo

    ! put PZ nano values into Pnano variable, removing any negative values

    do ii = 1,M                        ! counter through grid
       jj = (PZ_bins%nano-1) * M + ii ! counter into PZ

       if (PZ(jj)>0 .and. flagellates) then
          Pnano(ii) = PZ(jj)
       else
          Pnano(ii) = 0
       endif
    enddo


    ! put PZ nitrate values into NO variable, removing any negative values

    do ii = 1,M                        ! counter through grid
       jj = (PZ_bins%NO-1) * M + ii    ! counter into PZ

       if (PZ(jj)>0) then
          NO(ii) = PZ(jj)
       else
          NO(ii) = 0
       endif
    enddo

    ! put PZ ammonium values into NH variable, removing any negative values

    do ii = 1,M                        ! counter through grid
       jj = (PZ_bins%NH-1) * M + ii    ! counter into PZ

       if (PZ(jj)>0 .and. remineralization) then
          NH(ii) = PZ(jj)
       else
          NH(ii) = 0
       endif
    enddo

    ! put PZ silicon values into Sil variable, removing any negative values

    do ii = 1,M                        ! counter through grid
       jj = (PZ_bins%Sil-1) * M + ii    ! counter into PZ

       if (PZ(jj) > 0) then
          Sil(ii) = PZ(jj)
       else
          Sil(ii) = 0
       endif
    enddo


    ! put PZ detritus values into detr  variable, removing any negative values

    do kk = 1,D_bins
       do ii = 1,M                                 ! counter through grid
          jj = (PZ_bins%Quant + (kk-1)) * M + ii   ! counter through PZ

          if (PZ(jj)>0 .and. remineralization) then
             detr(kk,ii) = PZ(jj)
          else
             detr(kk,ii) = 0
          endif
       enddo
    enddo

    ! initialize transfer rates between the pools
    N%O_uptake%new = 0.
    N%H_uptake%new = 0.
    waste%medium = 0.
    waste%small = 0.
    N%remin = 0.

    ! phytoplankton growth: Nitrate and Ammonimum, conc. of micro plankton
    ! I_par is light, Temp is temperature 
    ! N ammonium and nitrate uptake rates (IN and OUT) incremented
    ! micro is the growth parameters for micro plankton (IN) and the growth rates 
    ! (OUT)
    ! Resp is respiration, Mort is mortality

    call p_growth(M, NO, NH, Sil, Pmicro, I_par, Temp, & ! in
         N, rate_micro, micro, &         ! in and out, in, out
         Resp_micro, Mort_micro)                    ! out

    ! put microplankton mortality into the medium detritus flux
    waste%medium = waste%medium + Mort_micro*Pmicro

    ! phytoplankton growth: Nitrate and Ammonimum, conc. of nano plankton
    ! I_par is light, Temp is temperature 
    ! N ammonium and nitrate uptake rates (IN and OUT) incremented
    ! nano is the growth parameters for micro plankton (IN) and the growth rates 
    ! (OUT)
    ! Resp_nano is respiration, Mort_nano is mortality

    call p_growth(M, NO, NH, Sil, Pnano, I_par, Temp, & ! in
         N, rate_nano, nano, &             ! in and out, in, out
         Resp_nano, Mort_nano)                     ! out

    waste%small = waste%small + Mort_nano*Pnano

!!!New quantity, bacterial 0xidation of NH to NO pool ==> NH^2
    N%bacteria(1:M) = rate_N%r * NH**2

    ! remineralization of detritus groups 1 and 2 (not last one)
    do kk = 1,D_bins-1
       N%remin(:) = N%remin(:) + Detritus(kk)%r * detr(kk,:)
    enddo

    IF (MINVAL(N%remin) < 0.) THEN
       PRINT "(A)","N%remin < 0. in derivs.f90"
       PRINT *,N%remin
       STOP
    END IF

    ! calculate the f-ratio
    do jj = 1,M 
       IF (N%O_uptake%new(jj) + N%H_uptake%new(jj) > 0.) THEN
          f_ratio(jj) = N%O_uptake%new(jj) /  &
               (N%O_uptake%new(jj) + N%H_uptake%new(jj))
       ELSE
          f_ratio(jj) = 0.
       END IF
    END DO

    ! now put it all together into the derivatives

    ! initialize derivatives
    dPZdt(1:M2) = 0.      

    ! microplankton

    do ii = 1,M ! index for Pmicro
       jj = (PZ_bins%micro-1) * M + ii ! index for PZ, derivatives etc

       if (Pmicro(ii) > 0.) then 
          dPZdt(jj) = (micro%growth%new(ii) - Resp_micro(ii) &
               - Mort_micro(ii)) * Pmicro(ii)
       endif
    enddo

    ! nanoplankton

    do ii = 1,M ! index for Pnano
       jj = (PZ_bins%nano-1) * M + ii ! index for PZ, derivatives etc

       if (Pnano(ii) > 0.) then 
          dPZdt(jj) = (nano%growth%new(ii) - Resp_nano(ii) &
               - Mort_nano(ii)) * Pnano(ii)
       endif
    enddo


    ! nitrate

    do ii = 1,M ! index for NO
       jj = (PZ_bins%NO-1) * M + ii ! index for PZ, derivatives etc

       if (NO(ii) > 0.) then  
          dPZdt(jj) = -N%O_uptake%new(ii) + N%bacteria(ii)
       endif
    enddo

    ! ammonium

    do ii = 1,M ! index for NH
       jj = (PZ_bins%NH-1) * M + ii ! index for PZ, derivatives etc

       if (NH(ii) > 0.) then  
          dPZdt(jj) = -N%H_uptake%new(ii) &
               + Resp_micro(ii)*Pmicro(ii) + Resp_nano(ii)*Pnano(ii)       & 
               + waste%medium(ii) * waste%m%destiny(0)                     &
               + waste%small(ii) * waste%s%destiny(0)                      &
               + N%remin(ii) - N%bacteria(ii)
       endif
    enddo

    ! silicon

    do ii = 1, M
       jj = (PZ_bins%Sil-1) * M + ii ! index for PZ, derivatives etc

       if (Sil(ii) > 0.) then
          dPZdt(jj) = - (micro%growth%new(ii) * Pmicro(ii) * &
                             rate_micro%Sil_rat &
                         + nano%growth%new(ii) * Pnano(ii) * &
                             rate_nano%Sil_rat)
       endif
    enddo

    ! detritus (not last bin)

    do ii = 1, M
       do kk = 1, D_bins-1
          jj = (PZ_bins%Quant + (kk-1)) * M + ii

          if (detr(kk,ii) > 0) then
             dPZdt(jj) = waste%medium(ii) * waste%m%destiny(kk)  &
                  + waste%small(ii) * waste%s%destiny(kk) &
                  - Detritus(kk)%r * detr(kk,ii) 
          end if
       end do
    end do

  end subroutine derivs_sog

end module biological_mod
