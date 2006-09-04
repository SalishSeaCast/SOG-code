! $Id$
! $Source$

SUBROUTINE Julian_day(leap,month,day,Jday,size)

  USE mean_param

  IMPLICIT NONE

  INTEGER, INTENT(IN)::size
  INTEGER, DIMENSION(size), INTENT(IN)::leap, month, day
  INTEGER, DIMENSION(size), INTENT(OUT)::Jday

  INTEGER::k
  DO k = 1, size
     IF (leap(k)== 1) THEN               ! leap year (2000, 1996, 1992, 1988, 1984...)
        SELECT CASE (month(k))        
        CASE (1)
           Jday(k) = day(k)
        CASE (2)
           Jday(k) = day(k) + 31
        CASE (3)                                
           Jday(k) = day(k) + 60
        CASE(4)
           Jday(k) = day(k) + 91
        CASE(5)
           Jday(k) = day(k) + 121
        CASE(6)          
           Jday(k) = day(k) + 152         
        CASE(7)          
           Jday(k) = day(k) + 182  
        CASE(8)          
           Jday(k) = day(k) + 213
        CASE(9)          
           Jday(k) = day(k) + 244
        CASE(10)          
           Jday(k) = day(k) + 274
        CASE(11)          
           Jday(k) = day(k) + 305
        CASE(12)          
           Jday(k) = day(k) + 335
        END SELECT
     ELSE 
        SELECT CASE (month(k))
        CASE (1)
           Jday(k) = day(k)
        CASE (2)
           Jday(k) = day(k) + 31
        CASE (3)                                
           Jday(k) = day(k) + 59
        CASE(4)
           Jday(k) = day(k) + 90
        CASE(5)
           Jday(k) = day(k) + 120
        CASE(6)          
           Jday(k) = day(k) + 151         
        CASE(7)          
           Jday(k) = day(k) + 181  
        CASE(8)          
           Jday(k) = day(k) + 212
        CASE(9)          
           Jday(k) = day(k) + 243
        CASE(10)          
           Jday(k) = day(k) + 273
        CASE(11)          
           Jday(k) = day(k) + 304
        CASE(12)          
           Jday(k) = day(k) + 334
        END SELECT
     END IF
  END DO
END SUBROUTINE Julian_day
