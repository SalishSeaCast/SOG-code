! $Id$
! $Source$

module declarations
! *** What's in here?

  use precision_defs, only: dp

  use mean_param

  implicit NONE

  ! Run parameters read from file specified as 'inputfile' in infile
  character*4 :: cruise_id  ! four number code that labels the start cruise

  type(UVST):: Gvector
  TYPE(plankton2)::micro, nano, pico
  TYPE(entrain)::euph

  ! Heat fluxes
  real(kind=dp) :: Q_t  ! Turbulent surface heat flux
  real(kind=dp), dimension(:), allocatable :: Q_n  ! Nonturb heat flux profile

  !surface heat flux
  DOUBLE PRECISION::Q_sol  !See irradiance.f90
  DOUBLE PRECISION::Br
  CHARACTER::ignored_input
  INTEGER, DIMENSION(20)::alloc_stat  !***
  INTEGER :: &
       time_step, index,xx, xx2, &
       jmax_i,    &
       count, yy, yy2, gg, jmaxg, cloud_type, &
       water_type
  REAL::D_test
  DOUBLE PRECISION :: &
       Io
  !I = intensity
  DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE:: I, I_par

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

end module declarations
