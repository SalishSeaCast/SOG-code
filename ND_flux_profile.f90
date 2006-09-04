! $Id$
! $Source$

SUBROUTINE ND_flux_profile(mm,L,phi1)
  ! Calculate the nondimensional flux profiles per Large, etal (1994),
  ! App. B
  ! *** Cleanup needed.  Constants should be declared here instead of
  ! *** in surface forcing, etc.
     
      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::mm
      TYPE(MS), INTENT(OUT)::phi1
      DOUBLE PRECISION,INTENT(IN)::L  !L_star

      INTEGER::y

     !Define phi1%m%value
      
      phi1%m%value = 0.
      phi1%s%value = 0.

      DO y = 0,mm%M
         IF (mm%d_i(y)/L >= 0.)THEN ! stable (defined as functions of d, the grid)
            phi1%m%value(y) = 1.0+5.0*mm%d_i(y)/L ! (B1a)
         ELSE IF (mm%d_i(y)/L < 0. .AND. xsi_m <= mm%d_i(y)/L)THEN
            phi1%m%value(y) = (1.0-16.0*mm%d_i(y)/L)**(-1.0/4.0) ! (B1b)
         ELSE
            phi1%m%value(y) = (a_m - c_m*mm%d_i(y)/L)**(-1.0/3.0) ! (B1c)
         END IF

      !Define phi1%s%value

         IF (mm%d_i(y)/L >= 0.)THEN
           phi1%s%value(y) = 1.0+5.0*mm%d_i(y)/L ! (B1a)
         ELSE IF (mm%d_i(y)/L < 0. .AND. xsi_s <= mm%d_i(y)/L)THEN
           phi1%s%value(y) = (1.0-16.0*mm%d_i(y)/L)**(-1.0/2.0) ! (B1d)
         ELSE
            phi1%s%value(y) = (a_s - c_s*mm%d_i(y)/L)**(-1.0/3.0) ! (B1e)
         END IF
      END DO

END SUBROUTINE ND_flux_profile
 
