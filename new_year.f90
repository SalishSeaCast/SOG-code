subroutine new_year(day_ti, o_day, o_year, ti, ddt, day_ch, day_ch2, &
     month)

  USE mean_param
  USE surface_forcing      

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(IN OUT)::day_ti, ti   !day_time, time
  DOUBLE PRECISION, INTENT(IN)::ddt             !dt
  INTEGER, INTENT(IN OUT)::o_day,o_year,day_ch,day_ch2      !day, year, day_check,day_check2
  INTEGER, INTENT(IN OUT)::month

  INTEGER::zz

  day_ti = day_ti + ddt
  ti = ti + ddt
      
  IF (day_ti >= 86400.0) THEN
     day_ti = day_ti - 86400.0
     o_day = o_day + 1
     day_ch = 0
     day_ch2 = 0
     IF (o_day == 366) THEN 
        PRINT "(A)","time,day,year"
        PRINT *,ti,o_day,o_year
        IF (is_leap_year == 0) THEN
           was_leap_year = 0
           o_year = o_year + 1        
           o_day = 1
        PRINT *,ti,o_day,o_year
           DO zz = 1,SIZE(leap_year,1)
              IF (o_year == leap_year(zz)) THEN
                 is_leap_year = 1
                 EXIT
              ELSE
                 is_leap_year = 0
              END IF
           END DO
        END IF
     ELSE IF (o_day == 367 .AND. is_leap_year == 1) THEN
        was_leap_year = 1
        o_year = o_year + 1        
        o_day = 1 
        is_leap_year = 0
     ELSE IF (o_day == 367 .AND. is_leap_year == 0) THEN
        PRINT "(A)", "Error in new_year.f90. is_leap_year,o_day:"
        PRINT *,is_leap_year,o_day
        STOP
     ELSE IF (o_day > 367) THEN
        PRINT "(A)", "Day > 367, see new_year.f90"
        PRINT *,o_day
        STOP
     END IF
  END IF


  SELECT CASE (month)
  CASE(1)
     IF (o_day > 31) THEN
        month = 2
     END IF
  CASE(2)
     IF (o_day > 59+is_leap_year) THEN
        month = 3
     END IF
  CASE(3)
     IF (o_day > 90+is_leap_year) THEN
        month = 4
     END IF
  CASE(4)
     IF (o_day > 120+is_leap_year) THEN
        month = 5
     END IF
  CASE(5)
     IF (o_day > 151+is_leap_year) THEN
        month = 6
     END IF
  CASE(6)
     IF (o_day > 181+is_leap_year) THEN
        month = 7
     END IF
  CASE(7)
     IF (o_day > 212+is_leap_year) THEN
        month = 8
     END IF
  CASE(8)
     IF (o_day > 243+is_leap_year) THEN
        month = 9
     END IF
  CASE(9)
     IF (o_day > 273+is_leap_year) THEN
        month = 10
     END IF
  CASE(10)
     IF (o_day > 305+is_leap_year) THEN
        month = 11
     END IF
  CASE(11)
     IF (o_day > 334+is_leap_year) THEN
        month = 12
     END IF
  CASE(12)
     IF (o_day == 1) THEN
        month = 1
     END IF
  END SELECT


END SUBROUTINE new_year




