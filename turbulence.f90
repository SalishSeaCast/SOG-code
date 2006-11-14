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
  !   nK -- Overall diffusivity profile; a continuous profile of K_ML%*
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
       nK, &      ! Overall diffusivity profile; a continuous profile of
                 ! K_ML%* in the mixing layer, and K%*%total below it
       wbar, &   ! Turbulent kinematic flux profiles
       ! *** Temporary until turbulence refactoring is completed
       nu, &  ! Interior diffusivity profile arrays
       u_star, &  ! Turbulent friction velocity
       w_star, &  ! Convective velocity scale
       L_mo, &    ! Monin-Obukhov length scale
       ! Subroutines:
       init_turbulence, calc_KPP_diffusivity, dalloc_turbulence_variables

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
       nK  ! Overall diffusivity of the water column; a continuous
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
  real(kind=dp) :: &
       u_star, &  ! Turbulent friction velocity
       w_star, &  ! Convective velocity scale
       L_mo       ! Monin-Obukhov length scale

contains

  subroutine init_turbulence(M)
    ! Allocate memory for turbulence model variables, and read
    ! parameter values from infile.

    ! Functions from other modules:
    use input_processor, only: getpard

    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Allocate memory for turbulence quantities.
    call alloc_turbulence_variables(M)
    ! Read parameter values from infile.
    !
    ! Internal wave breaking diffusivity for momentum and scalar quantities
    nu%m%int_wave = getpard('nu_w_m')
    nu%T%int_wave = getpard('nu_w_s')
    nu%S%int_wave = nu%T%int_wave
  end subroutine init_turbulence


  subroutine calc_KPP_diffusivity(Bf, h)
    ! Calculate the diffusivity profile using the K profile
    ! parameterization algorithm of Large, et al (1994).

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         Bf, &  ! Surface buoyancy forcing
         h      ! Mixing layer depth

    ! Calculate the turbulent momentum diffusivity due to shear
    ! instabilities (nu%m%shear) using the Large, et al (1994)
    ! parameterization based on the local gradient Richardson number.
    call shear_diffusivity()

    ! Calculate turbulence scales that characterize the mixing layer:
    ! turbulent friction velocity, convective velocity scale, and
    ! Monin-Obukhov length scale
    u_star = (wbar%u(0) ** 2 + wbar%v(0) ** 2) ** (1./4.)
    w_star = (-Bf * h) ** (1./3.)
    L_mo = u_star ** 3 / (kapa * Bf)
  end subroutine calc_KPP_diffusivity


  subroutine shear_diffusivity()
    ! Calculate the turbulent momentum diffusivity due to shear
    ! instabilities (nu%m%shear) using the Large, et al (1994)
    ! parameterization based on the local gradient Richardson number.

    ! Elements from other modules:
    ! Type Definitions:
    use precision_defs, only: dp
    ! Parameter Value Declarations:
    use fundamental_constants, only: &
         g  ! Acceleration due to gravity [m/s^2]
    use grid_mod, only: &
         grid  ! Grid parameters, and arrays
    use core_variables, only: &
         U, &  ! Cross-strait velocity component profile arrays; we
                                ! need the gradients at the grid layer interface
                                ! depths: U%grad_i
         V     ! Along-strait velocity component profile arrays; we
    ! need the gradients at the grid layer interface
    ! depths: V%grad_i
    use water_properties, only: &
         rho  ! Density profile arrays; we need the density, and its
    ! gradient at the grid layer interface depths: rho%i &
    ! rho%grad_i

    implicit none

    ! Local parameter value declarations:
    real(kind=dp), parameter :: &
         Ri_o = 0.7d0, &   ! Critical gradient Richardson number
         nu_o = 1.0d-3, &  ! Maximum shear diffusivity; Large, et al
                           ! (1994) recommends 50e-4
         p_1 = 3.0d0       ! Power constant for shear diffusivity
                           ! parameterization
    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: &
         Ri_g, &  ! Profile of Richardson number gradient at grid
                  ! layer interface depths
         N2,   &  ! Profile of buoyancy frequency squared at grid
                  ! layer interface depths
         V2       ! Profile of the square of the magnitude of the
                  ! velocity gradient at the grid layer interface
                  ! depths
    integer :: &
         i  ! Index over profile depth

    ! Calculate the profiles of the buoyancy frequency squared,
    ! velocity gradient squared, and Richardson number gradient at the
    ! grid layer interface depths.  (Large, et al (1994), eq'n (27))
    N2 = (-g / rho%i(1:grid%M)) * rho%grad_i(1:grid%M)
    V2 = U%grad_i(1:grid%M) ** 2 + v%grad_i(1:grid%M) ** 2
    Ri_g = N2 / (V2 + epsilon(V2))
    ! Apply the shear diffusivity parameterization (Large, et al
    ! (1994), eq'n(28))
    do i = 1, grid%M
       if (Ri_g(i) <= 0.) then
          ! Eq'n (28a)
          nu%m%shear(i) = nu_o
       elseif (0. < Ri_g(i) .and. Ri_g(i) < Ri_o) then
          ! Eq'n (28b)
          nu%m%shear(i) = nu_o * (1. - (Ri_g(i) / Ri_o) ** 2) ** p_1
       else
          ! Eq'n (28c)
          nu%m%shear(i) = 0.
       endif
    enddo
  end subroutine shear_diffusivity


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
    allocate(nK%m(1:M), nK%T(1:M), nK%S(1:M), &
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
    deallocate(nK%m, nK%T, nK%S, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Turbulent kinematic flux profile arrays"
    deallocate(wbar%u, wbar%v, wbar%t, wbar%s, &
         wbar%b, wbar%b_err,                   &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_turbulence_variables

end module turbulence
