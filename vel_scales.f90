! $Id$
! $Source$

SUBROUTINE vel_scales(mm,omeg,phig,u_st,L,hh)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::mm
      TYPE(height), INTENT(IN)::hh
      TYPE(MS), INTENT(IN)::phig
      TYPE(MS), INTENT(OUT)::omeg
      DOUBLE PRECISION, INTENT(IN)::u_st, L !friction velocity scale and Monin- 
                                      !Obukov Length
      INTEGER::k
      DOUBLE PRECISION::phime,phise
      TYPE(height)::surface_h

      surface_h%new = ep*hh%new  ! height of the surface layer is 0.1 hh
  
      CALL find_jmax_i(surface_h,mm) ! find grid point

      phime = 0.
      phise = 0.

      phime = (phig%m%value(surface_h%i-1) - phig%m%value(surface_h%i))*&
      (mm%d_i(surface_h%i)-&
      surface_h%new)/mm%i_space(surface_h%i)+  phig%m%value(surface_h%i)
! interpolate phi(surface layer) for equation (13)
      phise = (phig%s%value(surface_h%i-1) - phig%s%value(surface_h%i))*&
      (mm%d_i(surface_h%i)-&
      surface_h%new)/mm%i_space(surface_h%i)+ phig%s%value(surface_h%i)

      omeg%m%value = 0.
      omeg%m%value = 0.
      DO k = 0, mm%M !hh%i
         IF (mm%d_i(k)/L < 0. .AND. k >= surface_h%i) THEN  
            omeg%s%value(k) = kapa*u_st/phise
            omeg%m%value(k) = kapa*u_st/phime
         ELSE
            omeg%s%value(k) = kapa*u_st/phig%s%value(k)
            omeg%m%value(k) = kapa*u_st/phig%m%value(k)
         END IF
      END DO


!!! In unstable conditions, the d_i < h condition is ignored so that omegax
!!!  changes smoothly across the boundary layer depth and the derivative at
!!!  the BL_depth can be calculated.  Maybe, omegax should be determined only
!!! up to hh%i - 1 (minimum omegax(0) = omegax)!!!!

END SUBROUTINE vel_scales 






