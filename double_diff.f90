! $Id$
! $Source$

SUBROUTINE double_diff(mm,temp,sal,nus,alp_i, betai)

  use precision_defs, only: dp
  use grid_mod, only: grid_
      USE mean_param, only: prop, Knu, div_grid, div_interface
      USE surface_forcing

      IMPLICIT NONE

      TYPE(grid_),INTENT(IN)::mm  
      TYPE(prop), INTENT(IN OUT)::temp,sal      !T,S
      TYPE(Knu), INTENT(IN OUT)::nus
      REAL(KIND=DP), DIMENSION(0:mm%M), INTENT(IN)::alp_i, betai
      REAL(KIND=DP), DIMENSION(0:mm%M)::Ri_rho
      INTEGER::k

      CALL div_interface(mm, temp)
      CALL div_interface(mm, sal)
      CALL div_grid(mm, temp)
      CALL div_grid(mm, sal)

      Ri_rho = 0.0
      nus%s%dd = 0.0
      nus%t%dd = 0.0


! calculate double diffusion denisty ratio (30)
      DO k = 1, mm%M
         Ri_rho(k) = alp_i(k)*temp%div_i(k)/(betai(k)*sal%div_i(k))
      END DO
      
      DO k = 1, mm%M
         IF (1.0 < Ri_rho(k) .AND. Ri_rho(k) < 2.0 .AND. &
              alp_i(k)*temp%div_i(k) > 0 .AND. betai(k)*sal%div_i(k) > 0.) THEN  
                                !salt fingering!
              IF (1.0 < Ri_rho(k) .AND. Ri_rho(k) < Ri_rho_o) THEN

                 nus%s%dd(k) = nu_f*(1.0 - ((Ri_rho(k) - &
                 1.0)/(Ri_rho_o-1.0))**2.0)**p_2   ! (31a)

                 nus%t%dd(k) = 0.7*nus%s%dd(k) ! (31c)            

              END IF ! (31b) not needed, already set to 0

         ELSE IF (0. < Ri_rho(k) .AND. Ri_rho(k) < 1.0 .AND. &
            alp_i(k)*temp%div_i(k) < 0. .AND. betai(k)*sal%div_i(k) < 0.) THEN   
                                !Diffusive instability!

            nus%t%dd(k) =nu*0.909*EXP(4.6*EXP(-0.54*(1.0/Ri_rho(k)-1.0))) !(32)

            IF (Ri_rho(k) >= 0.5 .AND. Ri_rho(k) < 1.0) THEN
               nus%s%dd(k) = nus%t%dd(k)*(1.85 - 0.85/Ri_rho(k))*Ri_rho(k) !(34)
            ELSE IF (Ri_rho(k) < 0.5) THEN
               nus%s%dd(k) = nus%t%dd(k)*0.15*Ri_rho(k) ! (34)
            END IF
         ELSE
            nus%t%dd(k) = 0.0
            nus%s%dd(k) = 0.0
         END IF
      END DO


END SUBROUTINE double_diff
