! $Id$
! $Source$

module reaction

  use mean_param
  use surface_forcing
  use pdf

  implicit none

contains

  subroutine Copepod_growth(PZ, size, Zoo, species, cevent, mm, Csources, &
       C_types, waste, urea1, bin, bin_tot)
    ! *** What's it do?

    ! Arguments:
    INTEGER, INTENT(IN)::size,Csources,C_types,bin,bin_tot !size of PZ
    TYPE(gr_d), INTENT(IN)::mm !grid
    DOUBLE PRECISION, DIMENSION(size), INTENT(IN)::PZ
    TYPE(Cdata), DIMENSION(C_types), INTENT(IN)::Zoo
    TYPE(copepod), DIMENSION(Csources), INTENT(IN OUT)::species
    TYPE(event), DIMENSION(Csources), INTENT(IN)::cevent
    TYPE(losses),INTENT(IN OUT)::waste
    DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN OUT)::urea1  !N%urea%new

    INTEGER::j, k, x,total_prey,bin_tot2

    DOUBLE PRECISION, DIMENSION(prey+d_prey,mm%M,bin,Csources)::graze
    DOUBLE PRECISION, DIMENSION(mm%M,bin,Csources)::ingest
    DOUBLE PRECISION, DIMENSION(bin)::excrete
    DOUBLE PRECISION, DIMENSION(mm%M)::n_dz, n_z
    DOUBLE PRECISION, DIMENSION(prey+d_prey)::p_bin_size
    DOUBLE PRECISION::avg_Ex, Ntot
    DOUBLE PRECISION, DIMENSION(mm%M)::Pmicro,Zmicro,Pnano,SPN,totalprey
    DOUBLE PRECISION, DIMENSION(Csources,bin)::weight  !species(yy)%mature_pdf(1:Cevent(yy)%length)%wt_avg
    DOUBLE PRECISION, DIMENSION(Csources,mm%M)::Zcopepod !species(yy)%n(1:M)*species(yy)%Ntot = 
    !species(yy)%Z%new(1:M)

    total_prey = prey+d_prey   !Note: prey = 2, d_prey = 1, D_bins = 3, zprey = 3:   
    p_bin_size = 1.0  
    graze = 0.
    ingest = 0.
    bin_tot2 = 0
    weight = 0.

    Pmicro(1:mm%M) = PZ(1:mm%M)
    Zmicro(1:mm%M) = PZ(mm%M+1:2*mm%M)
    Pnano(1:mm%M) = PZ(2*mm%M+1:3*mm%M)
    !V.flagella.01=> should have been but since this sub is not called in this version does not matter=> Pnano(1:mm%M) = PZ(mm%M+1:2*mm%M)
    !sinking particulates
    SPN(1:mm%M) = PZ(5*mm%M+bin_tot+Csources*mm%M+mm%M+1:5*mm%M+bin_tot+Csources*mm%M+2*mm%M)
    DO k = 1, Csources  
       totalprey(1:mm%M) = Pmicro(1:mm%M)+Zoo(Cevent(k)%type)%q(1)**0.5*Zmicro(1:mm%M)+&
            Zoo(Cevent(k)%type)%q(2)**0.5*SPN(1:mm%M)
       bin_tot2 = bin_tot2 + Cevent(k)%length
       weight(k,1:Cevent(k)%length) = PZ(5*mm%M+(k-1)*mm%M + &
            bin_tot2-Cevent(k)%length+1:5*mm%M+bin_tot2+(k-1)*mm%M)
       Zcopepod(k,1:mm%M) = PZ(5*mm%M+bin_tot2+(k-1)*mm%M+1:5*mm%M+bin_tot2+k*mm%M)

       IF (species(k)%Ntot > small .AND. Cevent(k)%on /= 0) THEN

!!define ingestion averaged over z given n(z) (eg. of [P] for wt_j: 
       !!!= < (delta-eta)*(Gmax*ks*[P(z)]^2/(Gmax + ks*([P(z)]^2+q(1)*[Z]^2+q(2)*[SPN]^2))*n(z)*wt_j^(b_Ex-1)>_z
       !!!delta = assimilation efficiency ==> (1-delta)*grazing = egestion
       !!!eta = Excretion fraction ==> eta*grazing + Excretion(wt) = excretion

          !Find n(j):
          CALL pdf_avg(Zcopepod(k,1:mm%M),mm%i_space,mm%M,Ntot)
          IF (Ntot > 0.) THEN
             n_z(1:mm%M) = Zcopepod(k,1:mm%M)/Ntot
          ELSE
             PRINT "(A)","Ntot <= 0. in reaction.f90/Copepod_growth"
             PRINT *,Ntot
             n_z = 0.
          END IF
          DO j = 1,Cevent(k)%length
             IF (weight(k,j) > 0.) THEN
                ingest(:,j,k) = (Zoo(Cevent(k)%type)%delta-Zoo(cevent(k)%type)%eta)*&
                     weight(k,j)**(Zoo(Cevent(k)%type)%b2_Ex-1.)*Zoo(cevent(k)%type)%Gmax*&
                     Zoo(Cevent(k)%type)%ks*(Pmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(1)*Zmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(2)*SPN(:)**Zoo(Cevent(k)%type)%nn)/(Zoo(Cevent(k)%type)%Gmax+&
                     Zoo(Cevent(k)%type)%ks*(Pmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(1)*Zmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(2)*SPN(:)**Zoo(Cevent(k)%type)%nn)) 
                ! ingest(:,j,k) = (Zoo(Cevent(k)%type)%delta-Zoo(cevent(k)%type)%eta)*&
                !      weight(k,j)**(Zoo(Cevent(k)%type)%b_Ex-1.)*Zoo(cevent(k)%type)%Gmax*&
                !      Zoo(Cevent(k)%type)%ks*totalprey(:)**Zoo(Cevent(k)%type)%nn/&
                !      (Zoo(cevent(k)%type)%Gmax+Zoo(Cevent(k)%type)%ks*totalprey(:)**Zoo(Cevent(k)%type)%nn)
                n_dz = n_z(:)*mm%i_space(:)
                CALL pdf_avg(ingest(:,j,k),n_dz,mm%M,species(k)%ingest(j))
             END IF
          END DO
          DO j = 1,Cevent(k)%length
             IF (weight(k,j)*species(k)%mature_pdf(j)%f > 0.) THEN
                graze(1,:,j,k) = Zcopepod(k,:)*Zoo(cevent(k)%type)%Gmax*Zoo(Cevent(k)%type)%ks*&
                     weight(k,j)**(Zoo(Cevent(k)%type)%b2_Ex)*Pmicro(:)**Zoo(Cevent(k)%type)%nn/&
                     (Zoo(Cevent(k)%type)%Gmax+Zoo(Cevent(k)%type)%ks*(Pmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(1)*Zmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(2)*SPN(:)**Zoo(Cevent(k)%type)%nn))  
                graze(2,:,j,k) = Zcopepod(k,:)*Zoo(cevent(k)%type)%Gmax*Zoo(Cevent(k)%type)%ks*&
                     weight(k,j)**(Zoo(Cevent(k)%type)%b2_Ex)*Zoo(Cevent(k)%type)%q(1)*&
                     Zmicro(:)**Zoo(Cevent(k)%type)%nn/&
                     (Zoo(Cevent(k)%type)%Gmax+Zoo(Cevent(k)%type)%ks*(Pmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(1)*Zmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(2)*SPN(:)**Zoo(Cevent(k)%type)%nn))
                graze(3,:,j,k) = Zcopepod(k,:)*Zoo(cevent(k)%type)%Gmax*Zoo(Cevent(k)%type)%ks*&
                     weight(k,j)**(Zoo(Cevent(k)%type)%b2_Ex)*Zoo(Cevent(k)%type)%q(2)*&
                     SPN(:)**Zoo(Cevent(k)%type)%nn/&
                     (Zoo(Cevent(k)%type)%Gmax+Zoo(Cevent(k)%type)%ks*(Pmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(1)*Zmicro(:)**Zoo(Cevent(k)%type)%nn+&
                     Zoo(Cevent(k)%type)%q(2)*SPN(:)**Zoo(Cevent(k)%type)%nn))
                ! graze(1,:,j,k) = Zcopepod(k,:)*Zoo(cevent(k)%type)%Gmax*Zoo(Cevent(k)%type)%ks*&
                !      weight(k,j)**(Zoo(Cevent(k)%type)%b_Ex)*Pmicro(:)*totalprey(:)/&
                !      (Zoo(cevent(k)%type)%Gmax+&
                !      Zoo(Cevent(k)%type)%ks*totalprey(:)**Zoo(Cevent(k)%type)%nn)
                ! graze(2,:,j,k) = Zcopepod(k,:)*Zoo(cevent(k)%type)%Gmax*Zoo(Cevent(k)%type)%ks*&
                !      weight(k,j)**(Zoo(Cevent(k)%type)%b_Ex)*Zoo(Cevent(k)%type)%q(1)*&
                !      Zmicro(:)*totalprey(:)/(Zoo(cevent(k)%type)%Gmax+&
                !      Zoo(Cevent(k)%type)%ks*totalprey(:)**Zoo(Cevent(k)%type)%nn)
                ! graze(3,:,j,k) = Zcopepod(k,:)*Zoo(cevent(k)%type)%Gmax*Zoo(Cevent(k)%type)%ks*&
                !      weight(k,j)**(Zoo(Cevent(k)%type)%b_Ex)*Zoo(Cevent(k)%type)%q(2)*&
                !      SPN(:)*totalprey(:)/(Zoo(cevent(k)%type)%Gmax+&
                !      Zoo(Cevent(k)%type)%ks*totalprey(:)**Zoo(Cevent(k)%type)%nn)
             ELSE
                graze(:,:,j,k) = 0.
             END IF
          END DO
          DO x = 1,total_prey
             DO j = 1,mm%M  !Average over weights : species(k)%graze(x,j) loss of Prey l (gN/m^3) at grid j
                CALL pdf_avg(graze(x,j,1:Cevent(k)%length,k),species(k)%mature_pdf(:)%f,Cevent(k)%length,&
                     species(k)%graze(x,j))
             END DO
          END DO

!!!!!!Egestion of Copepods : summed over Csources (k) and prey (l) (source of detritus)!!!!!!

          DO x = 1,total_prey
             waste%medium(:) = waste%medium(:)+(1.0-Zoo(cevent(k)%type)%delta)*species(k)%graze(x,:)
          END DO

!!!!!!!Define excretion for individuals with weight wt_j:  Ex_j = a*(wt_j)^b

          DO j = 1,Cevent(k)%length
             IF (species(k)%mature_pdf(j)%f > 0.) THEN
                species(k)%Ex(j) = Zoo(cevent(k)%type)%a_Ex*weight(k,j)**Zoo(cevent(k)%type)%b_Ex
                excrete(j) = species(k)%Ex(j)
             ELSE
                species(k)%Ex(j) = 0.
                excrete(j) = 0.
             END IF
          END DO

          CALL pdf_avg(excrete(1:Cevent(k)%length),species(k)%mature_pdf(:)%f,Cevent(k)%length,avg_Ex)
          urea1(1:mm%M) = urea1(1:mm%M) + avg_Ex*Zcopepod(k,1:mm%M)
          !PRINT "(A)","avg_Ex"
          !PRINT *,avg_Ex
          DO x = 1,total_prey
             urea1(1:mm%M) = urea1(1:mm%M) + Zoo(cevent(k)%type)%eta*species(k)%graze(x,:)
          END DO
          !PRINT "(A)","urea1(1)"
          !PRINT *,urea1(1)
       END IF
    END DO

  END SUBROUTINE Copepod_growth

  SUBROUTINE Copepod_mortality(PZ,size,zoo,species,cevent,mm,Csources,C_types,waste,bin)

    INTEGER, INTENT(IN)::size,Csources,C_types,bin
    TYPE(gr_d),INTENT(IN)::mm
    DOUBLE PRECISION, DIMENSION(size), INTENT(IN)::PZ
    TYPE(Cdata), DIMENSION(C_types),INTENT(IN)::zoo
    TYPE(event), DIMENSION(Csources),INTENT(IN)::cevent
    TYPE(copepod), DIMENSION(Csources),INTENT(IN OUT)::species
    TYPE(losses),INTENT(IN OUT)::waste

    INTEGER::j,k,bin_tot2
    DOUBLE PRECISION::avg_wt, n2_avg_z
    DOUBLE PRECISION, DIMENSION(mm%M)::n_dz
    DOUBLE PRECISION, DIMENSION(Csources,bin)::weight  !species(yy)%mature_pdf(1:Cevent(yy)%length)%wt_avg
    DOUBLE PRECISION, DIMENSION(Csources,mm%M)::Zcopepod !species(yy)%n(1:M)*species(yy)%Ntot = 
    !species(yy)%Z%new(1:M)

    bin_tot2 = 0
    weight = 0.

    !d(n*Ntot)/dt = -species(k)%mort = -M*(n*Ntot)^2
    !d(D(3))/dt = M*(n*Ntot)^2*<wt> = waste%large
    !or for weight dependent mortality  d(n*Ntot)/dt = -M'*(n*Ntot)*sum([f(i)wavg(i)]**p)
    ! where p = 1 or -1
    !and for each cohort  f'(i) = Ntot/Ntot'(f'(i)-M'dt*(f(i)wavg(i)**p) where ' on f and Ntot indicates
    !next timestep value.  
    DO k = 1,Csources 
       bin_tot2 = bin_tot2 + cevent(k)%length

       weight(k,1:cevent(k)%length) = PZ(5*mm%M+(k-1)*mm%M + &
            bin_tot2-cevent(k)%length+1:5*mm%M+bin_tot2+(k-1)*mm%M)
       Zcopepod(k,1:mm%M) = PZ(5*mm%M+bin_tot2+(k-1)*mm%M+1:5*mm%M+bin_tot2+k*mm%M)

       CALL pdf_avg(weight(k,1:cevent(k)%length),species(k)%mature_pdf(:)%f,cevent(k)%length,avg_wt)

       IF (species(k)%Ntot > 0.) THEN
          DO j = 1,mm%M    
             species(k)%mort(j) = zoo(cevent(k)%type)%M*Zcopepod(k,j)  !**2.0
             waste%large(j) = waste%large(j) + species(k)%mort(j)*avg_wt
             IF (waste%large(j) < 0.) THEN
                PRINT "(A)","waste%large(j) < 0. in reaction.f90"
                PRINT *,waste%large,j
                STOP
             END IF
          END DO
       END IF

    END DO
  END SUBROUTINE Copepod_mortality

  SUBROUTINE p_growth(PZ,P,M2,I_par,mm,N,micro,TT,Resp) 

    TYPE(plankton2), INTENT(IN OUT)::micro  
    TYPE(gr_d), INTENT(IN)::mm  !grid
    DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::P !V.flagella.01 note: either Pmicro or Pnano
    INTEGER, INTENT(IN)::M2  
    DOUBLE PRECISION, DIMENSION(M2), INTENT(IN)::PZ 
    TYPE(nutrient), INTENT(IN OUT)::N
    DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::I_par  
    DOUBLE PRECISION, DIMENSION(mm%M+1), INTENT(IN)::TT 

    INTEGER::j
    DOUBLE PRECISION, DIMENSION(mm%M)::Uc, Oup_cell, Hup_cell, ratio  ! Carbon uptake
    DOUBLE PRECISION, DIMENSION(mm%M)::Nitrate,Ammonium,Rmax
    DOUBLE PRECISION, DIMENSION(2*mm%M)::Resp
!!!!!!!!!!!Define growth Due to I_par!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    Nitrate(1:mm%M) = PZ(mm%M+1:2*mm%M)
    Ammonium(1:mm%M) = PZ(2*mm%M+1:3*mm%M)
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
       IF (Nitrate(j) > small) THEN
          Oup_cell(j) = Rmax(j)*Nitrate(j)*micro%kapa*EXP(-micro%gamma_o*Ammonium(j))/&
               (micro%k + Nitrate(j)*micro%kapa*EXP(-micro%gamma_o*Ammonium(j))+Ammonium(j))*&
               (Nitrate(j)+Ammonium(j))**micro%N_x/(micro%N_o+(Nitrate(j)+Ammonium(j))**micro%N_x)
       ELSE
          Oup_cell(j) = 0.
       END IF
       IF (Ammonium(j) > small) THEN
          Hup_cell(j) = Rmax(j)*Ammonium(j)/&
               (micro%k +Nitrate(j)*micro%kapa*EXP(-micro%gamma_o*Ammonium(j))+Ammonium(j))*&
               (Nitrate(j)+Ammonium(j))**micro%N_x/(micro%N_o+(Nitrate(j)+Ammonium(j))**micro%N_x)
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

  subroutine p_graze(PZ, size, mm, zmicro, micro, N, waste, Csources, &
       bin_tot, Zoo, cevent)
    ! *** What's it do?

    ! Arguments:
    INTEGER, INTENT(IN)::size, Csources,bin_tot
    DOUBLE PRECISION, DIMENSION(size), INTENT(IN)::PZ   !PZ
    TYPE(gr_d), INTENT(IN)::mm
    TYPE(nutrient), INTENT(IN OUT)::N     !N
    TYPE(plankton2), INTENT(IN OUT)::zmicro, micro  !zmicro and micro
    TYPE(losses), INTENT(IN OUT)::waste
    TYPE(Cdata), INTENT(IN)::Zoo  !Zoo(1)
    TYPE(event),INTENT(IN)::cevent !Cevent(1)

    DOUBLE PRECISION, DIMENSION(zprey,mm%M)::graze
    DOUBLE PRECISION, DIMENSION(zprey)::p_bin_size
    INTEGER::kk, xx, j,total_prey    
    DOUBLE PRECISION, DIMENSION(mm%M)::Z_micro,Pnano,SUS,SPN  !suspended and sinking particulate nitrogen   

    Z_micro(1:mm%M) = PZ(mm%M+1:2*mm%M)
    Pnano(1:mm%M) = PZ(2*mm%M+1:3*mm%M)
    ! V.flagella.01 the same explanation as above=> Pnano(1:mm%M) = PZ(mm%M+1:2*mm%M)
    !suspended and sinking particulates
    SUS(1:mm%M) = PZ(5*mm%M+bin_tot+Csources*mm%M+1:5*mm%M+bin_tot+Csources*mm%M+mm%M)
    SPN(1:mm%M) = PZ(5*mm%M+bin_tot+Csources*mm%M+mm%M+1:5*mm%M+bin_tot+Csources*mm%M+2*mm%M)  

    p_bin_size = 1.0  

    DO kk = 1,mm%M
       zmicro%graze(1,kk) = zmicro%G*zmicro%ks*Pnano(kk)**zmicro%nn/&
            (zmicro%G+zmicro%ks*(Pnano(kk)**zmicro%nn+zmicro%q(1)*SUS(kk)**zmicro%nn+&
            zmicro%q(2)*SPN(kk)**zmicro%nn))
       zmicro%graze(2,kk) = zmicro%G*zmicro%ks*zmicro%q(1)*SUS(kk)**zmicro%nn/&
            (zmicro%G+zmicro%ks*(Pnano(kk)**zmicro%nn+zmicro%q(1)*SUS(kk)**zmicro%nn+&
            zmicro%q(2)*SPN(kk)**zmicro%nn))
       zmicro%graze(3,kk) = zmicro%G*zmicro%ks*zmicro%q(2)*SPN(kk)**zmicro%nn/&
            (zmicro%G+zmicro%ks*(Pnano(kk)**zmicro%nn+zmicro%q(1)*SUS(kk)**zmicro%nn+&
            zmicro%q(2)*SPN(kk)**zmicro%nn))
    END DO

    DO kk = 1,mm%M
       DO j = 1,zprey  !total growth per zmicro
          zmicro%growth%new(kk) = zmicro%growth%new(kk) + zmicro%graze(j,kk)*(zmicro%delta-zmicro%eta)
          !egested material
          waste%small(kk) = waste%small(kk) + (1.0-zmicro%delta)*zmicro%graze(j,kk)*Z_micro(kk)
          !excretion
          N%urea%new(kk) = N%urea%new(kk) + zmicro%eta*zmicro%graze(j,kk)*Z_micro(kk)
       END DO
       zmicro%mort%new(kk) = Zoo%Gmax*Zoo%ks*Zoo%molt_wt(0)**Zoo%b2_Ex*Zoo%q(1)*&
            Z_micro(kk)**Zoo%nn/(Zoo%Gmax+Zoo%ks*Zoo%q(1)*Z_micro(kk)**Zoo%nn)*120000./80.*&
            !cevent%nauplii = 120000
       zmicro%M_z 
       !=Z_micro(kk)%M_z is now a fraction of total arriving copepodites
    END DO                !total dead zmicro
    waste%medium(1:mm%M) = waste%medium(1:mm%M) + (1.-Zoo%delta)*zmicro%mort%new(1:mm%M) !*Z_micro(1:mm%M)
    ! N%urea%new(1:mm%M) =  N%urea%new(1:mm%M) + Zoo%eta*zmicro%mort%new(1:mm%M)

  END SUBROUTINE p_graze

end module reaction
