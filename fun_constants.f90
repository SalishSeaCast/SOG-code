! $Id$
! $Source$

SUBROUTINE fun_constants(u_st,w_st,L_mo,ww,Bf,hh)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      DOUBLE PRECISION,INTENT(OUT)::u_st, w_st, L_mo !friction velocity
                      !convective vel. scale, monin-obukov length
      TYPE(flux), INTENT(IN)::ww
      DOUBLE PRECISION,INTENT(IN)::Bf !current forcing
      DOUBLE PRECISION,INTENT(IN)::hh !current mld

      u_st = SQRT(SQRT(ww%u(0)**2.0+ww%v(0)**2.0))
      w_st = (-Bf*hh)**(1.0/3.0)
      L_mo = u_st**3.0/(kapa*Bf)



END SUBROUTINE
