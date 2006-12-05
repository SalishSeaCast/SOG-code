! $Id$
! $Source$

module declarations
! *** What's in here?

  use precision_defs, only: dp

  use mean_param

  implicit NONE

  ! Run parameters read from file specified as 'inputfile' in infile
  integer :: year_o         ! year in which run starts [yr]
  ! *** Why is month_o not declared here; rather it's local in SOG.f90
  integer :: day_o          ! year-day on which run starts [1-Jan = 1]
  integer :: t_o            ! start of run time on day_o [s]
                            ! *** t_o should match CTD profile time of day
  real(kind=dp) :: t_f      ! end of run time [s]
  real(kind=dp) :: dt       ! time step [s]
  character*4 :: cruise_id  ! four number code that labels the start cruise

  TYPE(height)::oh, surface_h
  TYPE(MS)::omega
  type(UVST):: Gvector
  TYPE(plankton2)::micro, nano
  TYPE(entrain)::euph

  DOUBLE PRECISION :: time    !*** core variable?

  ! Surface buoyancy forcing
  real(kind=dp) :: &
       Bf, &   ! Current iteration step value
       Bf_old  ! Previous iteration step value

  ! Fresh water flux quantities:
  ! *** Destined for a module (probably new surface_forcing) eventually
  logical :: Fw_surface  ! Add all of the fresh water on the surface?
  real(kind=dp) :: &
       Ft, &        ! Total fresh water flux
       Fw_scale, &  ! Fresh water scale factor for river flows
       Fw_depth     ! Depth to distribute fresh water flux over
  real(kind=dp), dimension(:), allocatable :: &
       Fw, &  ! Fresh water flux profile
       F_n    ! Fresh water contribution to salinity flux

  ! Heat fluxes
  real(kind=dp) :: Q_t  ! Turbulent surface heat flux
  real(kind=dp), dimension(:), allocatable :: Q_n  ! Nonturb heat flux profile

  double precision :: wt_r

  INTEGER::j_gamma !interface point corresponding to depth at which Radiation contributes to turbulent
  !surface heat flux
  DOUBLE PRECISION::Q_sol  !See irradiance.f90
  DOUBLE PRECISION::Br
  CHARACTER::ignored_input
  INTEGER, DIMENSION(20)::alloc_stat  !***
  INTEGER :: steps, time_step, index,xx, xx2, year, jmax_i,    &
       count, yy, yy2, gg, jmaxg, cloud_type, day, &
       water_type
  REAL::D_test
  DOUBLE PRECISION :: &
       day_time, Io, h_Ekman  !depth
  !friction vel, convective vel scale
  !monin_obukov length and julian day
  DOUBLE PRECISION::beta_t ! entrainment coefficient under convection
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::V_t_square, &
       Ri_b, N_2_g
  !I = intensity
  DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE:: I, I_par

  DOUBLE PRECISION, DIMENSION(1)::microQ1_p,nanoQ1_p
  ! Surface biological flux == 0
  DOUBLE PRECISION::pflux_o

  !Variables for printing test functions !

  DOUBLE PRECISION :: avg_T, read_var 

  !Variable for odeint.f,  rkqs.f and derivs.f

  DOUBLE PRECISION::next_time  !see Read_data.f90 and input file 
  ! input/forcing.dat
  INTEGER::M2, M3, N_ok, N_bad      

  DOUBLE PRECISION::vapour_pressure, prain ! (atmosphere at 17m) and current value for &

  ! Vertical profile of the entrainment velocity
  real(kind=dp), dimension(:), allocatable :: wupwell

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::f_ratio

end module declarations
