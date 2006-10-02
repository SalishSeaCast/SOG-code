! $Id$
! $Source$

SUBROUTINE shear_diff(mm, u_vel, v_vel, rho, nu_s)

  use precision_defs, only: dp
  use grid_mod, only: grid_
  use water_properties, only: water_property
      USE mean_param, only: prop, div_interface, div_grid
      USE surface_forcing

      IMPLICIT NONE

      TYPE(grid_), INTENT(IN)::mm  !
      TYPE(prop), INTENT(IN OUT)::u_vel, v_vel  !U, V
      type(water_property), intent(in) :: rho 
      REAL(KIND=DP), DIMENSION(0:mm%M),INTENT(OUT)::nu_s

      REAL(KIND=DP), DIMENSION(0:mm%M)::Rig
      REAL(KIND=DP), DIMENSION(mm%M)::N_2  !buoyancy freq squared 
      REAL(KIND=DP), DIMENSION(mm%M)::V_grad_sq
      INTEGER::k
 
      CALL div_interface(mm, u_vel) ! calculate du/dz in u%div_i
      CALL div_interface(mm, v_vel) ! calculate dv/dz in v%div_i
      CALL div_grid(mm, u_vel) ! cclculate 1/2 grid point off
      CALL div_grid(mm, v_vel) ! calculate 1/2 grid point off

      Rig = 0.0
      nu_s = 0.0
      ! *** Vectorizable
      DO k = 1, mm%M
         N_2(k) = -(g / rho%i(k)) * rho%grad_i(k) ! calculate N2

         V_grad_sq(k) = u_vel%div_i(k)**2.0+v_vel%div_i(k)**2.0
         Rig(k) = N_2(k)/(V_grad_sq(k) + 1.0D-30) ! grad Ri number for interior (27)
      END DO
      
!PRINT*,'Rig',Rig
!pause
      DO k = 1, mm%M  !not 0
         IF (Rig(k) <= 0.) THEN 
            nu_s(k) = nu_o  ! nu_o is max diffusivity = 0.0050
!PRINT*,'1yes'
         ELSE IF (0. < Rig(k) .AND. Rig(k) < Ri_o) THEN
!PRINT*,'2yes'
            nu_s(k) = nu_o*(1.0 - (Rig(k)/Ri_o)**2.0)**p_1 ! (28b)
         ELSE 
            nu_s(k) = 0. ! stable
!PRINT*,'3yes'
         END IF
      END DO

!!! Rig = +inf, puts nu_s = 0!!!

END SUBROUTINE shear_diff


      
