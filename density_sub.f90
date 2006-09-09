SUBROUTINE density_sub(Ti,Si,ro,mm)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      INTEGER, INTENT(IN)::mm  !M
      TYPE(prop), INTENT(IN)::Ti, Si
      DOUBLE PRECISION, DIMENSION(0:mm+1), INTENT(OUT)::ro  !density of salt water
      DOUBLE PRECISION::ro_w
      INTEGER::x

      ro = 0.0

      DO x = 1, mm+1
         ro_w = 999.842594 + (Ti%new(x)-273.15)*(6.793952E-02 + (Ti%new(x)- &
         273.15)*(-9.095290E-03 + &
                (Ti%new(x)-273.15)*(1.001685E-04 + (Ti%new(x)-273.15)*(-1.120083E-06 + &
                (Ti%new(x)-273.15)*6.536332E-09))))
         ro(x) = ro_w + Si%new(x)*(0.824493 + (Ti%new(x)-273.15)*(-4.0899E-03 + &
                 (Ti%new(x)-273.15)*(7.6438E-05 + &
                 (Ti%new(x)-273.15)*(-8.2467E-07 + (Ti%new(x)-273.15)*5.3875E-09)))) + &
                 Si%new(x)**(3.0/2.0)*(-5.72466E-03 + (Ti%new(x)-273.15)*(1.0227E-04 - &
                 (Ti%new(x)-273.15)*1.6546E-06)) &
                 + 4.8314E-04*Si%new(x)**(2.0)
      END DO
      ro(0) = ro(1)
!PRINT*,'density sub: ro(0),S(0),T(0)',ro(0),Si%new(1),Ti%new(1)

END SUBROUTINE density_sub      
