! $Id$
! $Source$

SUBROUTINE fun_constants(u_st,w_st,L_mo,ww,Bf,hh)

  use precision_defs, only: dp
      USE mean_param, only: flux
      USE turbulence, only: kapa

      IMPLICIT NONE

      REAL(KIND=DP),INTENT(OUT)::u_st, w_st, L_mo !friction velocity
                      !convective vel. scale, monin-obukov length
      TYPE(flux), INTENT(IN)::ww
      REAL(KIND=DP),INTENT(IN)::Bf !current forcing
      REAL(KIND=DP),INTENT(IN)::hh !current mld

      u_st = SQRT(SQRT(ww%u(0)**2.0+ww%v(0)**2.0))
      w_st = (-Bf*hh)**(1.0/3.0)
      L_mo = u_st**3.0/(kapa*Bf)



END SUBROUTINE
