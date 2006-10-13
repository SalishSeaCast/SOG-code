! $Id$
! $Source$

SUBROUTINE shear_diff(mm, U_grad_i, V_grad_i, rho, nu_s)

  use precision_defs, only: dp
  use grid_mod, only: grid_
  use water_properties, only: water_property
  use physics_model, only: g
      USE surface_forcing, only: nu_o, Ri_o, p_1

      IMPLICIT NONE

      TYPE(grid_), INTENT(IN)::mm  !
      real(kind=dp), dimension(1:mm%M), intent(in) :: &
           U_grad_i, &  ! Cross-strait vel component gradient at interfaces
           V_grad_i     ! Along-strait vel component gradient at interfaces
      type(water_property), intent(in) :: rho 
      REAL(KIND=DP), DIMENSION(0:mm%M),INTENT(OUT)::nu_s

      REAL(KIND=DP), DIMENSION(0:mm%M)::Rig
      REAL(KIND=DP), DIMENSION(mm%M)::N_2  !buoyancy freq squared 
      REAL(KIND=DP), DIMENSION(mm%M)::V_grad_sq
      INTEGER::k
 
      Rig = 0.0
      nu_s = 0.0
      ! *** Vectorizable
      DO k = 1, mm%M
         N_2(k) = -(g / rho%i(k)) * rho%grad_i(k) ! calculate N2
         V_grad_sq(k) = U_grad_i(k) ** 2 + V_grad_i(k) ** 2
         Rig(k) = N_2(k) / (V_grad_sq(k) + 1.0D-30) ! grad Ri number for interior (27)
      END DO
      
!PRINT*,'Rig',Rig
!pause
      DO k = 1, mm%M  !not 0
         IF (Rig(k) <= 0.) THEN 
            nu_s(k) = nu_o  ! nu_o is max diffusivity = 0.0050
!PRINT*,'1yes'
         ELSE IF (0. < Rig(k) .AND. Rig(k) < Ri_o) THEN
!PRINT*,'2yes'
            nu_s(k) = nu_o * (1.0 - (Rig(k) / Ri_o) ** 2) ** p_1 ! (28b)
         ELSE 
            nu_s(k) = 0. ! stable
!PRINT*,'3yes'
         END IF
      END DO

!!! Rig = +inf, puts nu_s = 0!!!

END SUBROUTINE shear_diff


      
