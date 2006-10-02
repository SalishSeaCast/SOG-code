! $Id$
! $Source$

SUBROUTINE def_v_t_sog(mm, hh, N_2, omeg_s, vt_sq, betat,Lstar)

  use precision_defs, only: dp
  use grid_mod, only: grid_
      USE mean_param, only: height
      USE surface_forcing

      IMPLICIT NONE

      TYPE(grid_), INTENT(IN)::mm   !grid
      TYPE(height), INTENT(IN)::hh   !h
      REAL(KIND=DP), DIMENSION(mm%M), INTENT(IN)::N_2  !N^2
      REAL(KIND=DP), DIMENSION(0:mm%M), INTENT(IN)::omeg_s  !omega%s%value
      REAL(KIND=DP), DIMENSION(mm%M), INTENT(OUT)::vt_sq   !V_t_square
      REAL(KIND=DP), INTENT(IN)::betat, Lstar   !beta_t (-0.2) and L_star

      REAL(KIND=DP), DIMENSION(mm%M)::omeg_g
      INTEGER::xx

      vt_sq = 0.

      IF (Lstar <= 0.0) THEN   !Doesn't make much difference!!!!
         DO xx = 1, mm%M   ! (13) surface layer versus mixed layer
            IF (mm%d_g(xx) >= ep*hh%new) THEN
               omeg_g(xx) = omeg_s(xx)
            ELSE IF (mm%d_g(xx) < ep*hh%new .AND. ep*hh%new <= mm%d_i(xx)) THEN
               omeg_g(xx) = omeg_s(xx-1) + (omeg_s(xx)-omeg_s(xx-1))*&
               (mm%d_g(xx)-mm%d_i(xx-1))/(ep*hh%new-mm%d_i(xx-1))
            ELSE
               omeg_g(xx) = omeg_s(xx)-(omeg_s(xx-1)-omeg_s(xx))*&
               (mm%d_g(xx)-mm%d_i(xx))/mm%i_space(xx)
            END IF
          END DO
      ELSE
          DO xx = 1, mm%M
             omeg_g(xx) = omeg_s(xx)-(omeg_s(xx-1)-omeg_s(xx))*&
             (mm%d_g(xx)-mm%d_i(xx))/mm%i_space(xx)
           END DO
      END IF

      DO xx = 1, mm%M 
         IF (N_2(xx) >= 0) THEN
            vt_sq(xx) = Cv*SQRT(-betat*N_2(xx)/c_s/ep)/&
                 (Ri_c*kapa**2.0)*mm%d_g(xx)*omeg_g(xx) ! (23)
         ELSE IF (N_2(xx) < 0) THEN !Unstable  Ri_b < Ri_c
            vt_sq(xx) =  0.0  ! 1.0D+30   ! 0.0 
         END IF
      END DO



END SUBROUTINE def_v_t_sog
      
