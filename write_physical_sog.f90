! $Id$
! $Source$

subroutine write_physical_sog(unow, vnow, euphotic)
! Output physical results.

  use mean_param
  use declarations
  use surface_forcing

  implicit none

  ! Arguments:
  double precision :: unow, vnow      
  type(entrain), intent(in) :: euphotic !euphotic%depth, euphotic%i 

  ! Local variables:
  integer :: jj, years, sizerun, jj_day, day_count
  double precision :: dh_dt_hour, Ks_ml, I_ml_hour
  integer :: isusan, nsal, isal(6), tsal
  double precision :: sal(6), dsal(6)

100 format(1X, E15.8)

  nsal = 6
  sal(1) = 28.75
  sal(2) = 29
  sal(3) = 29.25
  sal(4) = 28.75
  sal(5) = 29.
  sal(6) = 29.25
  do jj=1,nsal
     isal(jj) = 0
  enddo

  do isusan=1,grid%M
     do jj=1,nsal
        if (S%new(isusan).lt.sal(jj)) isal(jj) = isusan
     enddo
  enddo
  do jj=1,nsal
     tsal = isal(jj)
     if (tsal.lt.grid%M.and.S%new(tsal+1).ne.S%new(tsal)) then
        dsal(jj) = isal(jj) + (sal(jj)-S%new(tsal))/(S%new(tsal+1)-S%new(tsal))
     else
        dsal(jj) = isal(jj)
     endif
  enddo

111 format (17(1x,f10.4))   

  write(293, 111) time/3600., h%new, stress%u%new, stress%v%new,t%new(0), &
       s%new(0), u%new(1), Q_t(0), K%t%all(1), I_par(0), density%new(0), &
       Q_n(1)*100000.

  write(295, 111) time/3600., P%micro%new(0), P%nano%new(0), N%O%new(0), N%H%new(0), &
       Detritus(1)%D%new(0), Detritus(2)%D%new(0), Detritus(3)%D%new(0), &
       f_ratio(1)


  open(294,file="output/PON_flux.dat")
  write(294,*) NHflux_200,NOflux_200,PONflux_200

  jj_day = 0

  DO jj = 1,27
     IF (jj == 1 .AND. day < j_day(2)-7) THEN
        jj_day = 1
        EXIT
     ELSE IF (jj == 27 .AND. day >= j_day(26)+7) THEN
        jj_day = 27
        EXIT
     ELSE IF (day >= j_day(jj)-7 .AND. day < j_day(jj)+7) THEN
        jj_day = jj
        EXIT
     END IF
  END DO
  years = 1

  sizerun = 27*years  !number of data points interested in saving
  day_count = 0

  DO jj = 1,years
     IF (jj == year-1966) THEN  !year-1966
        day_count = (jj-1)*27+jj_day
        EXIT
     END IF
  END DO
  !   OPEN(UNIT = 160, FILE = "output/UV_20res.dat",STATUS = "REPLACE", &
  !       ACTION = "WRITE")

  IF (year == 1972 .AND. day_time <= 3600. .AND. day_time > 0.) THEN
     WRITE(160,100)REAL(U_ten),REAL(V_ten)
  END IF

  !Data to save at each timestep
  !  OPEN(UNIT = 152, FILE = "output/dhdt_ml_hourly_20res.dat",STATUS = "REPLACE", &
  !       ACTION = "WRITE")
  !  OPEN(UNIT = 153, FILE = "output/Ks_ml_hourly_20res.dat",STATUS = "REPLACE", &
  !       ACTION = "WRITE")
!!$  OPEN(UNIT = 154, FILE = "output/Ipar_ml_hourly_20res.dat",STATUS = "REPLACE", &
!!$       ACTION = "WRITE")
  !  OPEN(UNIT = 155, FILE = "output/time_hourly_20res.dat",STATUS = "REPLACE", &
  !       ACTION = "WRITE")

  IF (year == 1972 .AND. (day == 100 .OR. day == 101)) THEN
     dh_dt_hour = (h_m%new-h_m%old)/dt
     CALL average3(grid,K%s%all(1:grid%M),h_m%new,Ks_ml)
     CALL average3(grid,I_par(1:grid%M),h_m%new,I_ml_hour)
     WRITE(152,100)REAL(dh_dt_hour),REAL(h_m%new)
     WRITE(153,100)REAL(Ks_ml)
     WRITE(154,100)REAL(I_ml_hour),I_par(1)
     WRITE(155,100)REAL(time)
  END IF

  IF (year == 1972 .AND.  time_step  /= steps ) THEN  !.OR. year == 1962 
     IF (month == 3) THEN  !March  
        cntp_mar = cntp_mar + 1
        T_mar(1:M) = T_mar(1:M) + T%new(1:M)
        S_mar(1:M) = S_mar(1:M) + S%new(1:M)
        U_mar(1:M) = U_mar(1:M) + U%new(1:M)
        V_mar(1:M) = V_mar(1:M) + V%new(1:M)
        Ku_mar(1:M) = Ku_mar(1:M) + K%u%all(1:M)
        Ks_mar(1:M) = Ks_mar(1:M) + K%s%all(1:M)
        Kt_mar(1:M) = Kt_mar(1:M) + K%t%all(1:M)
     ELSE IF (month == 6) THEN !June
        cntp_jun = cntp_jun + 1
        T_jun(1:M) = T_jun(1:M) + T%new(1:M)
        S_jun(1:M) = S_jun(1:M) + S%new(1:M)
        U_jun(1:M) = U_jun(1:M) + U%new(1:M)
        V_jun(1:M) = V_jun(1:M) + V%new(1:M)
        Ku_jun(1:M) = Ku_jun(1:M) + K%u%all(1:M)
        Ks_jun(1:M) = Ks_jun(1:M) + K%s%all(1:M)
        Kt_jun(1:M) = Kt_jun(1:M) + K%t%all(1:M)
     ELSE IF (month == 9) THEN !Sept
        cntp_sep = cntp_sep + 1
        T_sep(1:M) = T_sep(1:M) + T%new(1:M)
        S_sep(1:M) = S_sep(1:M) + S%new(1:M)
        U_sep(1:M) = U_sep(1:M) + U%new(1:M)
        V_sep(1:M) = V_sep(1:M) + V%new(1:M)
        Ku_sep(1:M) = Ku_sep(1:M) + K%u%all(1:M)
        Ks_sep(1:M) = Ks_sep(1:M) + K%s%all(1:M)
        Kt_sep(1:M) = Kt_sep(1:M) + K%t%all(1:M)
     ELSE IF (month == 12) THEN !Dec
        cntp_dec = cntp_dec + 1
        T_dec(1:M) = T_dec(1:M) + T%new(1:M)
        S_dec(1:M) = S_dec(1:M) + S%new(1:M)
        U_dec(1:M) = U_dec(1:M) + U%new(1:M)
        V_dec(1:M) = V_dec(1:M) + V%new(1:M)
        Ku_dec(1:M) = Ku_dec(1:M) + K%u%all(1:M)
        Ks_dec(1:M) = Ks_dec(1:M) + K%s%all(1:M)
        Kt_dec(1:M) = Kt_dec(1:M) + K%t%all(1:M)
     END IF
  ELSE IF ( time_step == steps) THEN
     T_mar = T_mar/DBLE(cntp_mar) 
     S_mar = S_mar/DBLE(cntp_mar) 
     U_mar = U_mar/DBLE(cntp_mar) 
     V_mar = V_mar/DBLE(cntp_mar) 
     Ku_mar = Ku_mar/DBLE(cntp_mar) 
     Ks_mar = Ks_mar/DBLE(cntp_mar) 
     Kt_mar = Kt_mar/DBLE(cntp_mar) 
     T_jun = T_jun/DBLE(cntp_jun) 
     S_jun = S_jun/DBLE(cntp_jun) 
     U_jun = U_jun/DBLE(cntp_jun) 
     V_jun = V_jun/DBLE(cntp_jun) 
     Ku_jun = Ku_jun/DBLE(cntp_jun) 
     Ks_jun = Ks_jun/DBLE(cntp_jun) 
     Kt_jun = Kt_jun/DBLE(cntp_jun) 
     T_sep = T_sep/DBLE(cntp_sep) 
     S_sep = S_sep/DBLE(cntp_sep) 
     U_sep = U_sep/DBLE(cntp_sep) 
     V_sep = V_sep/DBLE(cntp_sep) 
     Ku_sep = Ku_sep/DBLE(cntp_sep) 
     Ks_sep = Ks_sep/DBLE(cntp_sep) 
     Kt_sep = Kt_sep/DBLE(cntp_sep) 
     T_dec = T_dec/DBLE(cntp_dec) 
     S_dec = S_dec/DBLE(cntp_dec) 
     U_dec = U_dec/DBLE(cntp_dec) 
     V_dec = V_dec/DBLE(cntp_dec) 
     Ku_dec = Ku_dec/DBLE(cntp_dec) 
     Ks_dec = Ks_dec/DBLE(cntp_dec) 
     Kt_dec = Kt_dec/DBLE(cntp_dec) 
     !  OPEN(UNIT = 156, FILE = "output/phys_mar_20res.dat",STATUS = "REPLACE", &
     !       ACTION = "WRITE")
     !  OPEN(UNIT = 157, FILE = "output/phys_jun_20res.dat",STATUS = "REPLACE", &
     !       ACTION = "WRITE") 
     !  OPEN(UNIT = 158, FILE = "output/phys_sep_20res.dat",STATUS = "REPLACE", &
     !       ACTION = "WRITE")
     !  OPEN(UNIT = 159, FILE = "output/phys_dec_20res.dat",STATUS = "REPLACE", &
     !       ACTION = "WRITE")
     WRITE(156,100)REAL(T_mar),REAL(S_mar),REAL(U_mar),REAL(V_mar),REAL(Ku_mar),REAL(Ks_mar),REAL(Kt_mar)
     WRITE(157,100)REAL(T_jun),REAL(S_jun),REAL(U_jun),REAL(V_jun),REAL(Ku_jun),REAL(Ks_jun),REAL(Kt_jun)
     WRITE(158, 100)REAL(T_sep),REAL(S_sep),REAL(U_sep),REAL(V_sep),REAL(Ku_sep),REAL(Ks_sep),REAL(Kt_sep)
     WRITE(159,100)REAL(T_dec),REAL(S_dec),REAL(U_dec),REAL(V_dec),REAL(Ku_dec),REAL(Ks_dec),REAL(Kt_dec)
  END IF

  IF (year >= 1967   .AND. year < 1968 .AND.   time_step /= steps  .AND. jj_day /= 0) THEN
     p_cnt(day_count) = p_cnt(day_count)+1
     SST(day_count) = SST(day_count) + T%new(0)
     SSS(day_count) = SSS(day_count) + S%new(0)
     hm_avg(day_count) = hm_avg(day_count) + h_m%new
     Ipar_o(day_count) = Ipar_o(day_count) + I_par(0)
!!!!!!!!! change to surface wind stress N/m^2 !!!!
     Uten_o(day_count) = Uten_o(day_count)+ stress%u%new !Uten_o(day_count) + U_ten 
     Vten_o(day_count) = Vten_o(day_count)+stress%v%new !Vten_o(day_count) + V_ten  !
     UVten_o(day_count) = UVten_o(day_count) + &
          SQRT(stress%u%new ** 2.0 + stress%v%new ** 2.0)

     Qflux(day_count) = Qflux(day_count) + Q_t(0)  !W/m^2
     Fflux(day_count) = Fflux(day_count) + w%s(0)  !Ft*S_o/rho_o  (PSU*m/s)
     CALL average3(grid,K%s%all(1:grid%M),h_m%new,Ks_ml)
     KsML(day_count) = KsML(day_count) + Ks_ml
     CALL average3(grid,I_par(1:grid%M),h_m%new,I_ml_hour)
     IparML(day_count) = IparML(day_count) + I_ml_hour 

  ELSE IF (time_step == steps) THEN  !time_step == steps
     DO xx = 1,sizerun        
        IF (DBLE(p_cnt(xx)) <= 0.) THEN
           SST(xx) = 0.
           SSS(xx) =  0.
           hm_avg(xx) = 0.
           Ipar_o(xx) =  0.
           Uten_o(xx) = 0.
           Vten_o(xx) =  0.
           UVten_o(xx) =  0.
           Qflux(xx) = 0.
           Fflux(xx) = 0.
           KsML(xx) = 0.
           IparML(xx) = 0.
        ELSE
           SST(xx) =  SST(xx)/DBLE(p_cnt(xx)) 
           SSS(xx) =  SSS(xx)/DBLE(p_cnt(xx))
           hm_avg(xx) = hm_avg(xx)/DBLE(p_cnt(xx)) 
           Ipar_o(xx) =Ipar_o(xx)/DBLE(p_cnt(xx)) 
           Uten_o(xx) = Uten_o(xx)/DBLE(p_cnt(xx))  
           Vten_o(xx) = Vten_o(xx)/DBLE(p_cnt(xx))  
           UVten_o(xx) = UVten_o(xx)/DBLE(p_cnt(xx))    
           Qflux(xx) = Qflux(xx)/DBLE(p_cnt(xx))
           Fflux(xx) = Fflux(xx)/DBLE(p_cnt(xx))  
           KsML(xx) = KsML(xx)/DBLE(p_cnt(xx))
           IparML(xx) = IparML(xx)/DBLE(p_cnt(xx))
        END IF
     END DO
     OPEN(UNIT = 150, FILE = "output/physical_o_20res.dat",STATUS = "REPLACE", &
          ACTION = "WRITE")
     WRITE(150,100)REAL(SST(1:sizerun)),REAL(SSS(1:sizerun)),REAL(hm_avg(1:sizerun))
     WRITE(150,100)REAL(Ipar_o(1:sizerun)),REAL(Uten_o(1:sizerun)),REAL(Vten_o(1:sizerun)),REAL(UVten_o(1:sizerun))
     WRITE(150,100)REAL(Qflux(1:sizerun)),REAL(Fflux(1:sizerun)),REAL(KsML(1:sizerun)),REAL(IparML(1:sizerun))
     p_cnt(1:sizerun) = 0
     SST(1:sizerun) = 0.
     SSS(1:sizerun) =  0.
     hm_avg(1:sizerun) = 0.
     Ipar_o(1:sizerun) =  0.
     Uten_o(1:sizerun) = 0.
     Vten_o(1:sizerun) =  0.
     UVten_o(1:sizerun) =  0.
     Qflux(1:sizerun) = 0.
     Fflux(1:sizerun) = 0.
     KsML(1:sizerun) = 0.
     IparML(1:sizerun) = 0.
  end if

end subroutine write_physical_sog
