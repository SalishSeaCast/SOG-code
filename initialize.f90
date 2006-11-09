! $Id$
! $Source$

SUBROUTINE initialize

      USE declarations

      IMPLICIT NONE

      time = dble(t_o)
      day = day_o   !initialize julian day
       year = year_o

      day_time = time

      j_gamma = 0  !h_gamma = 0 or wt_r contribution to gamma%t vanishes

      f_ratio = 0.
end subroutine initialize
