! $Id$
! $Source$

module surface_forcing
  ! *** This is a very poorly named module!

  use precision_defs, only: dp
!!$  use mean_param

  implicit none

  real(kind=dp), parameter :: PI = 3.141592653589793
  ! Coriolis factor
  ! *** This should be reworked so that the latitude is read as a run
  ! *** parameter, and f calculated for the run.
  real(kind=dp), parameter :: latitude = 49. + 7.517 / 60. ! station S3
!*** Correction pgf90 will not except an intrinsic in parameter statement  
!*** real(kind=dp), parameter :: f = 2. * (2. * PI / 86400.) * &
!***    sin(PI * latitude / 180.)
  real(kind=dp), parameter :: f = 1.1d-4
  ! Acceleration due to gravity
  real(kind=dp), parameter :: g = 9.81  ! m/s^2

  DOUBLE PRECISION, PARAMETER:: &
       small = 1.D-15, &
       small2 = 0.001, &  !smallest number of copepods out
       min_out = 2.5D-04, &   !1.D-05, &
       cutoff = 10.0, &
       var_c = 5.0, &   !***  !11.0
       zero =  0.
  DOUBLE PRECISION, PARAMETER:: Cv = 1.5  !1.99 !set to keep beta_t = -0.2
       !beta_t = -0.2 only for convection
  DOUBLE PRECISION, PARAMETER::S_ice = 0.0, & !ice salinity. unknown
       rho_ice = 1000.0, & !density ice, unknown
       kapa = 0.4 !Von Karman constant
  DOUBLE PRECISION, PARAMETER::a_s = -28.86, a_m = 1.26, &
       c_s =  98.96, c_m = 8.38, &
       xsi_m = -0.20, xsi_s = -1.0,&
       ep = 0.1  !epsilon
  DOUBLE PRECISION, PARAMETER::Ri_o = 0.7, & !critical grad Richardson #
       nu_o = 0.001, & !0.0050 !m^2/s interior diff constant
       p_1 = 3.0 ! shear diff power constant
  DOUBLE PRECISION, PARAMETER::C_star = 9.9
  DOUBLE PRECISION, PARAMETER::del_o = 0.10, &!tolerance
       A_stress = -0.6, & !N/m^2 Wind stress constant
       T_stress = 57600.0, &  !seconds or 16 hrs
       Q_o = 1368.0, & !1367.0? W/m^2  Solar constant
       Lat = 0.85521, &  ! 49 degrees  latitude
       Lon = 145., &  !(145oW)
       !           albedo = 0.061, &  !6% Large 1996
  albedo = 0.18, &  !KC 17% OCT.22 2004
       emiss = -1.0, & !emissivity Large 1996
       Stef_Boltz = 5.6697D-08, & !W/m^2/K^4
       R_v = 287.04, & !J/kg/K  gas constant for dry air
       mu_3 = 0.2, &   !(Denman)^(-1)  gN/m^2  attenuation coeff  
       handle = 1.0D-02  !gN(P)/gN(Z)/(handle time 's')   swim diffusivity
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




