! $Id$
! $Source$

SUBROUTINE irradiance_sog(cf, day_time, day, In, Ipar, d, &
     I_k, Qs, euphotic, Qriver, Pmicro_new, Pnano_new)

  use precision_defs, only: dp, sp
  use grid_mod, only: grid_
      USE mean_param, only: entrain
      USE surface_forcing

      IMPLICIT NONE
      ! Arguments:
      real(kind=sp), intent(in) :: cf        ! cloud fraction
      real(kind=dp), intent(in) :: day_time  ! day-second
      integer, intent(in) :: day             ! year-day 
      type(grid_), intent(in) :: d           ! grid
      real(kind=dp), dimension(0:d%M), intent(out) :: In, Ipar  
      integer, intent(out) :: I_k
      ! Integrated daily solar contribution to the heat flux
      real(kind=dp), intent(out) :: Qs   
      type(entrain), intent(out) :: euphotic !euphotic%depth, euphotic%i
      real, intent(in) :: Qriver
      real(kind=dp), dimension(0:d%M), intent(in) :: &
           Pmicro_new, &  ! Micro phytoplankton
           Pnano_new      ! Nano phytoplankton

      ! Local variables:
      INTEGER::k, check, of                          
      REAL(KIND=DP)::declination, hour, cos_Z, day_length, hour_angle, &
                        sunrise, sunset, Qso, a, b,KK
      REAL(KIND=DP):: II, cos_Z_max, IImax      
      REAL(KIND=DP), DIMENSION(0:d%M)::Iparmax

!!!Define Okta Cloud Model!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! based on Dobson and Smith, table 5

  type :: cloudy
     real(kind=dp) :: A, B, fraction   !Regression coefficients
  end type cloudy

  type :: okta
     type(cloudy), dimension(0:9) :: type
  end type okta


!!!KC-- new coefficients addded august, 2004. Check on the standard deviation values, but I don't think these are used anywhere, anyway
  type(okta)::cloud
      cloud%type(0) = cloudy(0.4641,0.3304,0.0945) 
      cloud%type(1) = cloudy(0.4062,0.3799,0.0186)  
      cloud%type(2) = cloudy(0.4129,0.3420,0.0501)
      cloud%type(3) = cloudy(0.4263,0.3212,0.0743)
      cloud%type(4) = cloudy(0.4083,0.3060,0.0723)
      cloud%type(5) = cloudy(0.3360,0.3775,0.0294)
      cloud%type(6) = cloudy(0.3448,0.3128,0.0226)
      cloud%type(7) = cloudy(0.3232,0.3259,0.0019)
      cloud%type(8) = cloudy(0.2835,0.2949,0.0081)
      cloud%type(9) = cloudy(0.1482,0.3384,0.1345)



      check = 0
      hour = (day_time/3600.0-12.0)*15.  !degrees
      declination = 23.45*PI/180.0*SIN((284.0+DBLE(day))/365.25*2.0*PI)  !radians

      a = SIN(declination)*SIN(Lat) 
      b = COS(declination)*COS(Lat)
      cos_Z = a+b*COS(PI/180.0*hour)      !solar elevation
      hour_angle = ACOS(-(TAN(Lat)*TAN(declination)))  !radians
      day_length = hour_angle/15.0*2.0*180.0/PI !hours
      sunrise = 12.0 - 0.5*day_length  !hours
      sunset = 12.0 + 0.5*day_length   !hours
      cos_Z_max = COS(declination-Lat)  !zenith angle


      Qso = Q_o*(1.0+0.033*COS(DBLE(day)/365.25*2.0*PI))*(1.0-albedo) !*(1.0-insol)
      if (cf .ge. 9) then
         of = int(cf) - 2
      elseif (cf .ge. 3) then
         of = int(cf) - 1
      else 
         of = int(cf)
      endif

      IImax = Qso*(cloud%type(of)%A + cloud%type(of)%B*cos_Z_max)*cos_Z_max  
      !!Jeffrey thesis page 124

      IF (day_time/3600.0 > sunrise .AND. day_time/3600.0 < sunset) THEN    
          II = Qso*(cloud%type(of)%A + cloud%type(of)%B*cos_Z)*cos_Z   
      ELSE
         II = 0.
      END IF

! so II is the incoming radiation

      Qs = Qso*(((cloud%type(of)%A+cloud%type(of)%B*a)*a+cloud%type(of)%B*b**2.0/2.0)*&
           day_length+&
           (cloud%type(of)%A*b+2.0*cloud%type(of)%B*a*b)*&
           180.0/PI/15.0*(SIN(PI/180.0*(sunset-12.0)*15.0)-SIN(PI/180.0*(sunrise-12.0)*15.0)) + &
           cloud%type(of)%B*b**2.0*&
           180.0/PI/60.0*(SIN(PI/180.0*(sunset-12.0)*30.0)-SIN(PI/180.0*(sunrise-12.0)*30.0)))/24.0 

! Qs is the daily integrated value

!     PRINT "(A)","II"
!      PRINT *,II
!      PRINT "(A)","Qs"
!      PRINT *,Qs
!pause
!----------------------------------------------------------
!KC Oct.12 2004

      Ipar = 0.        !PAR
      In = 0.          !total light!

      In(0) =  II        !total light begins attenuation at 0.5m
      Ipar(0) = II*0.44  !44% of total light is PAR at surface (Jerlov)
      Iparmax(0) = IImax*0.44

     ! plank=Pmicro_new(1)+Pmicro_new(2)+Pmicro_new(3)+Pmicro_new(4)+Pmicro_new(5)
     ! plank=plank/5;

      !added V.flagella.01
      KK = 0.0146 * (Pmicro_new(1) + Pnano_new(1)) &
           + Qriver * 3.6597e-05 - 0.0110 * 3 + 0.2697
      !*** set mixed layer depth in this parameterization to a constant 3 m -- 
      ! otherwise a shallow mixed layer in winter is actually a disadvantage.
      ! this parameterization needs to be revisited.
      do k = 1, d%M    
         Ipar(k) = Ipar(0) * exp(-d%d_i(k) * KK)
         Iparmax(k) = Iparmax(0) * exp(-d%d_i(k) * KK)
         if (Ipar(k) < 0.01 * Ipar(0) .and. check == 0) then            
            I_k = k                                          
            check = 1                 
         end if               
         ! Total light for heat budget
         In(k) = 0.70 * II * exp(-d%d_i(k) * (0.8102 * KK + 1.1854)) &
              + 0.30 * II * exp(-d%d_i(k) * (0.8226 * KK - 0.0879))
      end do

open(556,file="output/light.dat")
write(556,*)KK


!----------------------------------------------------------
 

     ! splitting up the incoming depending on the penetration

     euphotic%i = 0
     euphotic%g = 0
     euphotic%depth = 0.
     DO k = 1,d%M
!        IF (Iparmax(k) <= Iparmax(0)*0.01) THEN
        IF (Ipar(k) <= 1.4) THEN  ! saturation light intensity for Thalasosira
           euphotic%i = k
           euphotic%g = k
           euphotic%depth = d%d_g(k)
           if (Ipar(k) == 0) euphotic%depth = 0.
           EXIT
        END IF
     END DO


END SUBROUTINE irradiance_sog
