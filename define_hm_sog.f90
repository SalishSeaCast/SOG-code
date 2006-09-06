! $Id$
! $Source$

SUBROUTINE define_hm_sog(grid,S,h_m)

! define based on salinity rather than temperature
! use 0.02 rather than 0.1 C as density is more sensitive to salinity

  USE mean_param

  IMPLICIT NONE

  TYPE(gr_d), INTENT(IN)::grid
  DOUBLE PRECISION, DIMENSION(0:grid%M+1), INTENT(IN)::S  !S%new
  TYPE(height), INTENT(OUT)::h_m  !mixed

  INTEGER::kk
  DOUBLE PRECISION::d_S

  DO kk = 1,grid%M+1
     d_S = S(kk) - S(1)
!     write (*,*) S(1), S(kk),d_S
     IF (d_S >= 0.02) THEN
        h_m%g = kk
        h_m%new = grid%d_g(kk)-(grid%g_space(kk-1))/(S(kk-1)-S(kk))*(0.02-d_S)
        EXIT
     ELSE IF (kk == grid%M+1) THEN
        PRINT "(A)","Mixed layer bottomed out.  See define_hm.f90"
        h_m%g = kk
        h_m%new = grid%d_g(kk)
     END IF
  END DO

  IF (h_m%new < 0.) THEN
     PRINT "(A)","Negative h_m, see define_hm.f90"
     PRINT "(A)","h_m%new,kk,grid%M+1,grid%d_g(kk),T(kk),T(kk-1),d_T,grid%g_space(kk-1)"
     PRINT *,h_m%new,kk,grid%M+1,grid%d_g(kk),S(kk),S(kk-1),d_S,grid%g_space(kk-1)
     STOP
  END IF

END SUBROUTINE define_hm_sog

