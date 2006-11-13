! $Id$
! $Source$

SUBROUTINE def_gamma( L,  mm, ww, wtr, hh, gamm, Bf, omeg)

  use precision_defs, only: dp
  use grid_mod, only: grid_
  use turbulence, only: &
       c_s, &  ! Coefficient of phi%s in 1/3 power law regime
       kapa    ! von Karman constant
      USE mean_param, only: flux, height, MST, MS
      USE surface_forcing, only: ep, c_star

      IMPLICIT NONE

      REAL(KIND=DP), INTENT(IN)::L,  wtr, Bf !L_star, wt_r, Bf
      TYPE(flux), INTENT(IN)::ww            !w
      TYPE(grid_), INTENT(IN)::mm            !grid
      TYPE(height), INTENT(IN)::hh          !h 
      TYPE(MST), INTENT(OUT)::gamm          !gamma
      TYPE(MS), INTENT(IN)::omeg            !omega
      
      INTEGER::xx

      gamm%m = 0.0
      gamm%s = 0.0
      gamm%t = 0.0 
 

      DO xx = 1, hh%i  
         IF ((L /= 0.0 .AND. mm%d_i(xx)/L < 0) .OR. (L == 0.0 .AND. Bf < 0)) THEN
           IF (mm%d_i(xx) > hh%new) THEN   
              gamm%t(xx) = 0.
              gamm%s(xx) = 0.
           ELSE
               gamm%t(xx) = C_star*kapa*(c_s*kapa*ep)**(1.0/3.0)*(ww%t(0)+wtr)/&
                            omeg%s%value(xx)/hh%new
               gamm%s(xx) = C_star*kapa*(c_s*kapa*ep)**(1.0/3.0)*ww%s(0)/omeg%s%value(xx)/&
                            hh%new
           END IF
         END IF
      END DO
END SUBROUTINE def_gamma
