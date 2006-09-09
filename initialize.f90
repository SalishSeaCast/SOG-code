! $Id$
! $Source$

SUBROUTINE initialize

      USE mean_param
      USE surface_forcing
      USE declarations

      IMPLICIT NONE

      time = t_o
      dummy_time = 0.
      day = day_o   !initialize julian day
     ! month = month_o   !KC june 28, 2004
       year = year_o

      day_time = time
      time_step = 1

      j_gamma = 0  !h_gamma = 0 or wt_r contribution to gamma%t vanishes

      neg_count = 0  ! counts number of times scheme biology becomes negative

!!!!Initial CN ratios for P%micro and P%nano

      micro%Q_cn = Q_min !6.67  
      nano%Q_cn = Q_min !6.67 

      f_ratio = 0.
      null_vector = 0.
end subroutine initialize
