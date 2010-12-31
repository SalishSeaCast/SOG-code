module surface_forcing
  ! *** This is a very poorly named module!

  use precision_defs, only: dp

  implicit none


  DOUBLE PRECISION, PARAMETER:: &
       small = 1.D-15, &
       small2 = 0.001, &  !smallest number of copepods out
       min_out = 2.5D-04, &   !1.D-05, &
       cutoff = 10.0, &
       var_c = 5.0, &   !***  !11.0
       zero =  0.
  DOUBLE PRECISION, PARAMETER:: &
       Q_o = 1368.0, & !1367.0? W/m^2  Solar constant
       !           albedo = 0.061, &  !6% Large 1996
       albedo = 0.18  !KC 17% OCT.22 2004

  !***
  DOUBLE PRECISION, DIMENSION(13), PARAMETER::leap_year = (/ 1956, 1960, 1964, 1968, 1972,1976, 1980, &
       1984, 1988, 1992, 1996, 2000, 2004 /)
  INTEGER, PARAMETER:: runsize = 15 ! runsize = 15 !, &  !used in smoothdata.f90
  ! bin = 50 !50, &  !number of copepod weight bins in wt_pdf
  DOUBLE PRECISION, PARAMETER::min_bin = 0.005, & !minimum bin size for copepod wt_pdf is 0.0001*Zoo%wt
       sigma_max = 10.  !12 is too big ==> stage 1's increase slightly. see bye_zoo
  DOUBLE PRECISION, PARAMETER::desired_var = 6.267101794384283, &  !12.6, & !6.267101794384283, & !used for ECMWF smoothed winds
       desired_avg =   10.06734127539335 
  DOUBLE PRECISION :: &
       precision = 1.0D-4, &
       step_min = 3., &
       step_guess = 100.0
end module surface_forcing




