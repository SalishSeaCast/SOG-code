! $Id$
! $Source$

SUBROUTINE interior_match2(omeg, L, u_st, hh, mm)

  use precision_defs, only: dp
  use grid_mod, only: grid_
      USE mean_param, only: MS, height
      USE surface_forcing

      IMPLICIT NONE

      TYPE(MS), INTENT(IN OUT)::omeg       !omega
      REAL(KIND=DP), INTENT(IN)::L, u_st     !L_star, u_star
      TYPE(height), INTENT(IN)::hh         !h
      TYPE(grid_), INTENT(IN)::mm

      IF (L /= 0) THEN          
         IF (hh%new/L <= 0) THEN !!!Unstable and Neutral Conditions
            omeg%s%div = 0
            omeg%m%div = 0
            omeg%s%h = omeg%s%value(hh%i) 
            omeg%m%h = omeg%m%value(hh%i) 
         ELSE IF (hh%new/L > 0) THEN !!!Stable
            omeg%s%div =  &
            -kapa*u_st/((1.0+5.0*hh%new/L)**2.0)*5.0*hh%new/L ! d/dz of (B1a) in (13)
            omeg%s%h = kapa*u_st/(1.0+5.0*hh%new/L) ! (B1a) in (13)
            omeg%m%div = omeg%s%div 
            omeg%m%h = omeg%s%h
         END IF
      ELSE IF (L == 0.0) THEN
         omeg%s%div = 0
         omeg%m%div = 0
         omeg%s%h = omeg%s%value(hh%i) 
         omeg%m%h = omeg%m%value(hh%i)
      END IF

      
END SUBROUTINE interior_match2
