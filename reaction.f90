! $Id$
! $Source$

module reaction

  use mean_param
  use surface_forcing
  use pdf

  implicit none

contains

  SUBROUTINE p_growth(NO,NH,P,M2,I_par,mm,N,micro,TT,Resp) 

    double precision, dimension(mm%M) :: NO, NH ! nitrate and ammonium conc.

    TYPE(plankton2), INTENT(IN OUT)::micro  
    TYPE(gr_d), INTENT(IN)::mm  !grid
    DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::P !V.flagella.01 note: either Pmicro or Pnano
    INTEGER, INTENT(IN)::M2  
    TYPE(nutrient), INTENT(IN OUT)::N
    DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::I_par  
    DOUBLE PRECISION, DIMENSION(mm%M+1), INTENT(IN)::TT 

    INTEGER::j
    DOUBLE PRECISION, DIMENSION(mm%M)::Uc, Oup_cell, Hup_cell, ratio  ! Carbon uptake
    DOUBLE PRECISION, DIMENSION(mm%M)::Rmax
    DOUBLE PRECISION, DIMENSION(2*mm%M)::Resp
!!!!!!!!!!!Define growth Due to I_par!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    micro%growth%light = 0.
    ratio = 0.
    Oup_cell = 0.
    Hup_cell = 0.
    
! flagella_jun06/ Next loop is different from V.16.1 => Rmax(replaces micro%R in KPP) 
! is temp dependent in SOG 

    DO j = 1,mm%M
       Rmax(j)=micro%R*1.88**(0.1*(TT(j)-273.15-20))
       !       Resp(j)=(micro%Rm)*1.88**(0.1*(TT(j)-273.15-20))
       Resp(j)=(micro%Rm)*1.88**(0.1*(TT(j)-273.15-20))
       Resp(j+mm%M)=(micro%M_z)*1.88**(0.1*(TT(j)-273.15-20))
       micro%growth%light(j) = Rmax(j)*(1.0-EXP(-micro%sigma*I_par(j)/Rmax(j))) 
       Uc(j) = (1.0-micro%gamma)*micro%growth%light(j)*(1/(1+micro%inhib*I_par(j))) !- micro%Rm
    END DO 
!!!!!!!!!!!!!!!Define growth due to nutrients!!!!!!!!!!!!!!!!!!!!!!!!!

    DO j = 1,mm%M
       IF (NO(j) > small) THEN
          Oup_cell(j) = Rmax(j)*NO(j)*micro%kapa*EXP(-micro%gamma_o*NH(j))/&
               (micro%k + NO(j)*micro%kapa*EXP(-micro%gamma_o*NH(j))+NH(j))*&
               (NO(j)+NH(j))**micro%N_x/(micro%N_o+(NO(j)+NH(j))**micro%N_x)
       ELSE
          Oup_cell(j) = 0.
       END IF
       IF (NH(j) > small) THEN
          Hup_cell(j) = Rmax(j)*NH(j)/&
               (micro%k +NO(j)*micro%kapa*EXP(-micro%gamma_o*NH(j))+NH(j))*&
               (NO(j)+NH(j))**micro%N_x/(micro%N_o+(NO(j)+NH(j))**micro%N_x)
       ELSE
          Hup_cell(j) = 0.
       END IF
       IF (Oup_cell(j) > small) THEN
          ratio(j) = Hup_cell(j)/Oup_cell(j)
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

    DO j = 1,mm%M
       IF (Uc(j) >= 0. .AND. Uc(j) < Oup_cell(j) + Hup_cell(j)) THEN 
          !Carbon uptake (light) controls growth and nutrient uptake 
          !LIGHT LIMITING
          micro%growth%new(j) = Uc(j) !- micro%Rm
          IF (Uc(j) <= Hup_cell(j)) THEN
             N%H_uptake%new(j) = Uc(j)*P(j)+N%H_uptake%new(j)
             N%O_uptake%new(j) = 0.
          ELSE
             N%H_uptake%new(j) = Hup_cell(j)*P(j)+N%H_uptake%new(j)
             N%O_uptake%new(j) = (Uc(j) - Hup_cell(j))*P(j)+N%O_uptake%new(j)
          END IF
       ELSE IF (Uc(j) >= 0. .AND. Uc(j) >= Oup_cell(j) + Hup_cell(j)) THEN
          !nutrient uptake controls growth
          !NUTRIENT LIMITING
          micro%growth%new(j) = Oup_cell(j) + Hup_cell(j)  !- micro%Rm

          IF (micro%growth%new(j) < 0.) THEN
             micro%growth%new(j) = 0.
          ENDIF

          N%O_uptake%new(j) = Oup_cell(j)*P(j)+N%O_uptake%new(j)
          N%H_uptake%new(j) = Hup_cell(j)*P(j)+N%H_uptake%new(j)
       ELSE  !No nutrient uptake, no growth
          micro%growth%new(j) =  0!- micro%Rm
       END IF
    END DO


    !---mortality-------------------------
    !    micro%mort%new = Zoo%Gmax*Zoo%ks*Zoo%molt_wt(0)**Zoo%b2_Ex*&
    !           P**Zoo%nn/(Zoo%Gmax+Zoo%ks*P**Zoo%nn)*120000./80.*micro%M_z   !cevent%nauplii == 120000



  END SUBROUTINE p_growth


end module reaction
