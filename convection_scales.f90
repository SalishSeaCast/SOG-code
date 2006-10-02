! $Id$
! $Source$

SUBROUTINE convection_scales(mm,omeg,hh,w_st)

  use precision_defs, only: dp
  use grid_mod, only: grid_
      USE mean_param, only: height, MS
      USE surface_forcing

      IMPLICIT NONE

      TYPE(grid_), INTENT(IN)::mm
      TYPE(height), INTENT(IN)::hh
      TYPE(MS), INTENT(OUT)::omeg
      REAL(KIND=DP), INTENT(IN)::w_st !convective velocity
      TYPE(height)::surface_h 
      INTEGER::k
 
      surface_h%new = ep*hh%new

      omeg%s%value = 0
      omeg%m%value = 0

      DO k = 0, mm%M !hh%i
         IF (mm%d_i(k) < surface_h%new) THEN
            omeg%s%value(k) = kapa*(c_s*kapa*mm%d_i(k)/hh%new)**(1.0/3.0)*w_st
            omeg%m%value(k) = kapa*(c_m*kapa*mm%d_i(k)/hh%new)**(1.0/3.0)*w_st
         ELSE IF (surface_h%new <= mm%d_i(k)) THEN
            omeg%s%value(k) = kapa*(c_s*kapa*ep)**(1.0/3.0)*w_st
            omeg%m%value(k) = kapa*(c_m*kapa*ep)**(1.0/3.0)*w_st
         ELSE
            PRINT "(A)","Error, convection_scales.f"
         END IF
      END DO


END SUBROUTINE convection_scales
