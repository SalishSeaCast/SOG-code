! $Id$
! $Source$

module reaction

  implicit none

  private
  public :: p_growth

contains

  SUBROUTINE p_growth(M, NO, NH, P, I_par, Temp, N, plank, Resp, Mort) 

! this subroutine calculates the growth, respiration and mortality
! (light limited or nutrient limited) of either phytoplankton class.
! All three are functions of temperature

  use mean_param, only: plankton2, nutrient
  use surface_forcing, only: small

  implicit none

    integer, INTENT(IN)::M  !grid size

    ! Nitrate and ammonium concentrations
    double precision, dimension(M), intent(in) :: NO, NH 
    ! plankton concentraton (either Pmicro or Pnano)
    DOUBLE PRECISION, DIMENSION(M), INTENT(IN)::P !V.flagella.01 note: either Pmicro or Pnano

    DOUBLE PRECISION, DIMENSION(0:M), INTENT(IN):: I_par  ! light
    DOUBLE PRECISION, DIMENSION(0:M), INTENT(IN):: Temp  ! temperature


    ! in are the parameters of the growth equations, out are the growth values
    TYPE(plankton2), INTENT(IN OUT):: plank ! either micro or nano 

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
       Rmax(j)=plank%R*temp_effect 

       Resp(j)=(plank%Rm)*temp_effect 
       Mort(j)=(plank%M_z)*temp_effect

       plank%growth%light(j) = Rmax(j)*(1.0-EXP(-plank%sigma*I_par(j)/Rmax(j)))

       Uc(j) = (1.0-plank%gamma)*plank%growth%light(j)* &
            (1/(1+plank%inhib*I_par(j)))
    END DO 

!!!!!!!!!!!!!!!Define growth due to nutrients!!!!!!!!!!!!!!!!!!!!!!!!!

    DO j = 1,M

       NH_effect = exp(-plank%gamma_o * NH(j))
       IF (NO(j) > small) THEN
          Oup_cell(j) = Rmax(j) * NO(j) * plank%kapa * NH_effect / &
               (plank%k + NO(j) * plank%kapa * NH_effect + NH(j)) * &
               (NO(j) + NH(j))**plank%N_x / &
               (plank%N_o + (NO(j) + NH(j))**plank%N_x)
       ELSE
          Oup_cell(j) = 0.
       END IF
       IF (NH(j) > small) THEN
          Hup_cell(j) = Rmax(j) * NH(j) / &
               (plank%k + NO(j) * plank%kapa * NH_effect + NH(j))* &
               (NO(j) + NH(j))**plank%N_x / &
               (plank%N_o + (NO(j) + NH(j))**plank%N_x)
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

end module reaction
