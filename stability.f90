! $Id$
! $Source$

subroutine stability

  use turbulence, only: u_star
  use declarations, only: Bf, stable

      IMPLICIT NONE

            IF (Bf < 0.) THEN                              !!!unstable
               stable = 0
            ELSE IF (Bf > 0. .OR.(Bf == 0. .AND. u_star /= 0.)) THEN
               stable = 1                                  !!!stable
            ELSE IF (Bf == 0. .AND. u_star == 0.) THEN      !!!no forcing
               stable = 2
            ELSE
               PRINT "(A)","stability case not considered.  See stable.f, Bf, u_star:"
               PRINT *,Bf,u_star
               STOP
            END IF

END SUBROUTINE stability     



