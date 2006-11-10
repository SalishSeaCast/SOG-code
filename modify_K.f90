! $Id$
! $Source$

SUBROUTINE modify_K(mm, hh, X)

  use precision_defs, only: dp
  use grid_mod, only: grid_
      USE mean_param, only: diff, height

      IMPLICIT NONE

      TYPE(grid_), INTENT(IN)::mm
      TYPE(diff), INTENT(IN OUT)::X
      TYPE(height), INTENT(IN)::hh
      REAL(KIND=DP)::delta, K_star, K_grid

      !!Find modified X%ML(hh%g-1) for case(2) and X%total(hh%g-1) for case(1)!!
      !!or X%ML(hh%i-1) for case(2) and X%total(hh%i) for case(1)!!


      delta = (hh%new-mm%d_g(hh%g-1))/mm%g_space(hh%g-1) 

      IF (hh%g == 1) THEN  !! Actually minimum hh%g == 2
         
         K_grid = 0.0
         K_star = 0.0

      ELSE IF (hh%g > 1) THEN
          IF (hh%i == hh%g) THEN              !!!Condition (2) change X%ML

             IF (hh%g == 2) THEN
                K_grid = X%ML(1) + (mm%d_g(1)-mm%d_i(1))*X%ML(1)/mm%i_space(1)
                K_star = (1.0 - delta)**2.0*K_grid + delta**2.0*X%ML(1)
                X%ML(1) = (1.0 - delta)*X%total(1) + delta*K_star
                IF (X%ML(1) < 0.) THEN
                   PRINT "(A)","X%ML(1) < 0. : See modify_K"
                   PRINT *,X%ML(1)
                   STOP
                END IF
             ELSE
                K_grid = X%ML(hh%g-1) + (mm%d_g(hh%g-1)-mm%d_i(hh%g-1))*&
                (X%ML(hh%g-1) - X%ML(hh%g-2))/mm%i_space(hh%g-1)                
                K_star = (1.0-delta)**2.0*K_grid + delta**2.0*X%ML(hh%g-1)
                X%ML(hh%g-1)= (1.0-delta)*X%total(hh%g-1)+delta*K_star
                IF (X%ML(hh%g-1) < 0.) THEN
                   PRINT "(A)","X%ML(hh%g-1) < 0. : See modify_K"
                   PRINT *,X%ML(hh%g-1)
                   STOP
                END IF
             END IF

          ELSE IF (hh%i == hh%g - 1) THEN !!!Condition (1) change X%total
 
             IF (hh%g == 2) THEN
                K_grid = X%total(1) + (mm%d_g(1) - mm%d_i(1))*X%total(1)/mm%d_i(1)
                K_star = (1.0 - delta)**2.0*K_grid + delta**2.0*X%total(1)
                X%total(1) = (1.0 - delta)*X%total(1) + delta*K_star
                IF (X%total(1) < 0.) THEN
                   PRINT "(A)","X%total(1) < 0. Check Modify_K.f"
                   PRINT *,X%total(1)
                   STOP
                END IF
             ELSE
                K_grid = X%ML(hh%g-2) + (mm%d_g(hh%g-1)-mm%d_i(hh%g-2))*&
                (X%total(hh%g-1)-X%ML(hh%g-2))/mm%i_space(hh%g-1)
                K_star = (1.0-delta)**2.0*K_grid + delta**2.0*X%total(hh%g-1)
                X%total(hh%g-1) = (1.0-delta)*X%total(hh%g-1)+delta*K_star
                IF (X%total(hh%g-1) < 0.) THEN
                   PRINT "(A)","X%total(hh%g-1) < 0. Check Modify_K.f"
                   PRINT *,X%total(hh%g-1)
                   STOP
                END IF
             END IF
          END IF
       ELSE
          PRINT "(A)","Error, modify_K.f"
          write (*,*) 'modify-K'
          STOP
       END IF

 
END SUBROUTINE modify_K

