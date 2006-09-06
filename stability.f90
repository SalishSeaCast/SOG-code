! $Id$
! $Source$

SUBROUTINE stability

      USE mean_param
      USE declarations

      IMPLICIT NONE

            IF (Bf%b(0) < 0.) THEN                              !!!unstable
               stable = 0
            ELSE IF (Bf%b(0) > 0. .OR.(Bf%b(0) == 0. .AND. u_star /= 0.)) THEN
               stable = 1                                  !!!stable
            ELSE IF (Bf%b(0) == 0. .AND. u_star == 0.) THEN      !!!no forcing
               stable = 2
            ELSE
               PRINT "(A)","stability case not considered.  See stable.f, Bf%b(0), u_star:"
               PRINT *,Bf%b(0),u_star
               STOP
            END IF

END SUBROUTINE stability     



