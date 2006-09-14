! $Id$
! $Source$

SUBROUTINE shape_parameters(X, omeg_x, hh, a2, a3)

      USE mean_param

      IMPLICIT NONE

      TYPE(diff), INTENT(IN)::X
      TYPE(boundary), INTENT(IN)::omeg_x
      DOUBLE PRECISION, INTENT(IN)::hh  !Mixed layer depth  h
      DOUBLE PRECISION, INTENT(OUT)::a2, a3  !parameters for shape function G

      DOUBLE PRECISION::G_1, div_G_1

      G_1 = X%h/(hh*omeg_x%h) ! (18) match interior
      div_G_1 = -X%div/omeg_x%h - X%h*omeg_x%div/(omeg_x%h**2*hh) ! (18)

      IF (G_1 < 0.) THEN  ! diffusivities should not be less than 0
         G_1 = 0.
         div_G_1 = 0.
      ELSE IF (div_G_1 > 0.) THEN  !Interior can only increase BL mixing
         div_G_1 = 0.
      END IF

      a2 = -2.0 + 3.0*G_1-div_G_1
      a3 = 1.0 - 2.0*G_1 + div_G_1


END SUBROUTINE shape_parameters
