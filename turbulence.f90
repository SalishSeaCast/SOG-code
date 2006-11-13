! $Id$
! $Source$

module turbulence
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the KPP turbulence parameterization in the
  ! SOG code.
  !
  ! Public Parameters:
  !
  !   Constants for the nondimensional flux profiles.  See Large, et al
  !   (1994) pg 392, eqn B2
  !      xsi_m -- Maximum xsi value of the -1/3 power law regime of phi%m
  !      xsi_s -- Maximum xsi value of the -1/3 power law regime of phi%s
  !      a_m --  Coefficient of phi%m in 1/3 power law regime
  !      a_s --  Coefficient of phi%s in 1/3 power law regime
  !      c_m  -- Coefficient of phi%m in 1/3 power law regime
  !      c_s --  Coefficient of phi%s in 1/3 power law regime
  !      kapa -- von Karman constant
  !
  ! Public Variables:
  !
  !   K -- Overall diffusivity profile; a continuous profile of K_ML%*
  !        in the mixing layer, and K%*%total below it
  !
  !   wbar -- Turbulent kinematic flux profiles.  Note that only the
  !           surface values (index = 0) are used in the model, but that
  !           code to calculate the full profiles is available in ***.
  !           It is commented out, but can be included in a build if it
  !           is required for diagnostic purposes.
  !
  ! Public Subroutines:
  !
  !   init_turbulence --
  !
  !   dalloc_turbulence_variables -- 

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       xsi_m, &  ! Max xsi value of the -1/3 power law regime of phi%m
       xsi_s, &  ! Max xsi value of the -1/3 power law regime of phi%s
       a_m,   &  ! Coefficient of phi%m in 1/3 power law regime
       a_s,   &  ! Coefficient of phi%s in 1/3 power law regime
       c_m,   &  ! Coefficient of phi%m in 1/3 power law regime
       c_s,   &  ! Coefficient of phi%s in 1/3 power law regime
       kapa,  &  ! von Karman constant
       ! Variables:
       K, &       ! Overall diffusivity profile; a continuous profile of
                  ! K_ML%* in the mixing layer, and K%*%total below it
       wbar, & ! Turbulent kinematic flux profiles
       ! Subroutines:
       alloc_turbulence_variables, dalloc_turbulence_variables

  ! Type Definitions:
  !
  ! Public:
  !
  ! Profile arrays for momentum, temperature, and salinity components
  ! of diffusivities
  type :: mTS_arrays
     real(kind=dp), dimension(:), allocatable :: &
          m, &  ! Momentum diffusivity profile array
          T, &  ! Temperature diffusivity profile array
          S     ! Salinity diffusivity profile array
  end type mTS_arrays
  !
  ! Profile arrays for turbulent kinematic fluxes
  type :: turbulent_fluxes
     real(kind=dp), dimension(:), allocatable :: &
          u,     &  ! u velocity component (cross-strait)
          v,     &  ! v velocity component (along-strait)
          t,     &  ! temperature
          s,     &  ! salinity
          b,     &  ! buoyancy
          b_err     ! buoyancy error due to ???
  end type turbulent_fluxes
  !
  ! Private to module:
  !
  ! Interior diffusivity components
  type :: interior_diffusivity
     real(kind=dp), dimension(:), allocatable :: &
          shear, &  ! Diffusivity due to vertical shear
          dd        ! Diffusivity due to double diffusion
     real(kind=dp) :: &
          int_wave  ! Diffusivity due to internal wave breaking
     real(kind=dp), dimension(:), allocatable :: &
          total     ! Total interior diffusivity (sum of above 3 components)
     real(kind=dp) :: &
          tot_h, &    ! Total interior diffusivity at mixing layer depth
          tot_grad_h  ! Vertical gradient of total at mixing layer depth
  end type interior_diffusivity
  !
  ! Components for momentum, temperature, and salinity interior diffusivities
  type :: mTS_components
     type(interior_diffusivity) :: &
          m, &  ! Elements of interior diffusivity for momentum
          T, &  ! Elements of interior diffusivity for temperature
          S     ! Elements of interior diffusivity for salinity
  end type mTS_components

  ! Parameter Value Declarations:
  !
  ! Public:
  real(kind=dp), parameter :: &
       xsi_m = -0.20, &  ! Max xsi value of the -1/3 power law regime of phi%m
       xsi_s = -1.0,  &  ! Max xsi value of the -1/3 power law regime of phi%s
       a_m = 1.26,    &  ! Coefficient of phi%m in 1/3 power law regime
       a_s = -28.86,  &  ! Coefficient of phi%s in 1/3 power law regime
       c_m = 8.38,    &  ! Coefficient of phi%m in 1/3 power law regime
       c_s = 98.96,   &  ! Coefficient of phi%s in 1/3 power law regime
       kapa = 0.4        ! von Karman constant
  !
  ! Private to module:

  ! Variable Declarations:
  !
  ! Public:
  type(mTS_arrays) :: &
       K  ! Overall diffusivity of the water column; a continuous
          ! profile of K_ML%* in the mixing layer, and K%*%total below it
  type(turbulent_fluxes) :: &
       wbar  ! Turbulent kinematic flux profiles. 
             ! *** Note that only the surface values (index = 0) are used
             ! *** in the model, but that code to calculate the full
             ! *** profiles is available in ***.  It is commented out, but
             ! *** can be included in a build if it is required for diagnostic
             ! *** purposes.
  !
  ! Private to module:
  type(mTS_components) :: &
       nu  ! Interior diffusivity
           !   Components:
           !     %m  -- momentum components
           !     %T  -- temperature
           !     %S  -- salinity
           !       %shear      -- vertical shear
           !       %int_wave   -- internal wave breaking
           !       %dd         -- double diffusion
           !       %total      -- sum of above components
           !       %tot_h      -- value at mixing layer depth
           !       %tot_grad_h -- gradient at mixing layer depth
  type(mTS_arrays) :: &
       K_ML  ! Mixing layer diffusivity

contains

  subroutine init_turbulence()
    !
  end subroutine init_turbulence


  subroutine calc_KPP_diffusivity()
    ! Calculate the diffusivity profile using the K profile
    ! parameterization algorithm of Large, et al (1994).
    
  end subroutine calc_KPP_diffusivity


  subroutine alloc_turbulence_variables(M)
    ! Allocate memory for turbulence quantities.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Interior diffusivity profile arrays"
    allocate(nu%m%shear(1:M), nu%T%dd(1:M), nu%S%dd(1:M),   &
         nu%m%total(1:M), nu%T%total(1:M), nu%S%total(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Mixing layer diffusivity profile arrays"
    allocate(K_ML%m(1:M), K_ML%T(1:M), K_ML%S(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Overall diffusivity profile arrays"
    allocate(K%m(1:M), K%T(1:M), K%S(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Turbulent kinematic flux profile arrays"
    allocate(wbar%u(0:M), wbar%v(0:M), wbar%t(0:M), wbar%s(0:M), &
         wbar%b(0:M), wbar%b_err(0:M),                           &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_turbulence_variables


  subroutine dalloc_turbulence_variables()
    ! Deallocate memory for turbulence quantities.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Interior diffusivity profile arrays"
    deallocate(nu%m%shear, nu%T%dd, nu%S%dd, &
         nu%m%total, nu%T%total, nu%S%total, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Mixing layer diffusivity profile arrays"
    deallocate(K_ML%m, K_ML%T, K_ML%S, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Overall diffusivity profile arrays"
    deallocate(K%m, K%T, K%S, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Turbulent kinematic flux profile arrays"
    deallocate(wbar%u, wbar%v, wbar%t, wbar%s, &
         wbar%b, wbar%b_err,                   &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_turbulence_variables

end module turbulence
