module turbulence
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the KPP turbulence parameterization in the
  ! SOG code.
  !
  ! Public Parameters:
  !
  !   Constants for the nondimensional flux profiles.  See Large, et al
  !   (1994) pg 392, eqn B2
  !      zeta_m -- Maximum zeta value of the -1/3 power law regime of
  !                non-dimensional turbulent momentum flux profile
  !      zeta_s -- Maximum zeta value of the -1/3 power law regime of
  !                non-dimensional turbulent scalar flux profile
  !      a_m -- Coefficient of  non-dimensional turbulent momentum flux
  !             profile in 1/3 power law regime
  !      a_s -- Coefficient of  non-dimensional turbulent scalar flux
  !             profile in 1/3 power law regime
  !      c_m -- Coefficient of  non-dimensional turbulent momentum flux
  !             profile in 1/3 power law regime
  !      c_s -- Coefficient of  non-dimensional turbulent scalar flux
  !             profile in 1/3 power law regime
  !
  !   kapa -- von Karman constant
  !
  !   epsiln -- Non-dimensional extent of the surface layer
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
  !   u_star -- Turbulent friction velocity
  !
  !   L_mo -- Monin-Obukhov length scale
  !
  ! Public Subroutines:
  !
  !   init_turbulence -- Allocate memory for turbulence model
  !                      variables, and read parameter values from
  !                      infile.
  !
  !   calc_KPP_diffusivity -- Calculate the turbulent diffusivity
  !                           profile using the K profile
  !                           parameterization (KPP) algorithm of
  !                           Large, et al (1994).
  !
  !   calc_turbulent_fluxes -- Calculate the turbulent kinematic flux
  !                            profiles.  Note that only the surface
  !                            values (index = 0) are used in the
  !                            model.  This subroutine calculates the
  !                            full profiles.  It is not called, but
  !                            can be included in a build if the
  !                            profiles are required for diagnostic
  !                            purposes.
  !
  !   dalloc_turbulence_variables -- Deallocate memory for turbulence
  !                                  quantities.
  !
  ! Public Functions:
  !
  !   w_scale -- Return the vertical turbulent velocity scale value
  !              at the specified depth.
  !
  !   nondim_scalar_flux -- Return the value of the non-dimensional
  !                         scalar flux profile at the specified depth
  !                         (Large, etal (1994), App. B).
  !
  !   nonlocal_scalar_transport -- Return the value of the non-local
  !                                transport term for the specified
  !                                flux, and depth.  (Large, et al
  !                                (1994), eqn 20).

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       zeta_m, &  ! Max zeta value of the -1/3 power law regime of
                  ! non-dimensional turbulent momentum flux profile
       zeta_s, &  ! Max zeta value of the -1/3 power law regime of
                  ! non-dimensional turbulent scalar flux profile
       a_m,    &  ! Coefficient of non-dimensional turbulent momentum
                  ! flux profile in 1/3 power law regime
       a_s,    &  ! Coefficient of non-dimensional turbulent scalar
                  ! flux profile in 1/3 power law regime
       c_m,    &  ! Coefficient of non-dimensional turbulent momentum
                  ! flux profile in 1/3 power law regime
       c_s,    &  ! Coefficient of non-dimensional turbulent scalar
                  ! flux profile in 1/3 power law regime
       kapa,   &  ! von Karman constant
       epsiln, &  ! Non-dimensional extent of the surface layer
       ! Variables:
       K,      &  ! Overall diffusivity profile; a continuous profile
                  ! of K_ML%* in the mixing layer, and nu%*%total
                  ! below it
       wbar,   &  ! Turbulent kinematic flux profiles
       u_star, &  ! Turbulent friction velocity
       L_mo,   &  ! Monin-Obukhov length scale
       ! *** Temporary until turbulence refactoring is completed
       w, &  ! Turbulent velocity scale profile arrays
             ! values at the mixing layer depth
       ! Types (as required by new pg compiler)
       turbulent_fluxes, & ! type for wbar
       mTS_arrays, & ! type for k
       ! Subroutines:
       init_turbulence, calc_KPP_diffusivity, calc_turbulent_fluxes, &
       dalloc_turbulence_variables, &
       ! Functions:
       w_scale, nondim_scalar_flux, nonlocal_scalar_transport

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
  ! Interior diffusivity elements
  type :: interior_diffusivity_elements
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
  end type interior_diffusivity_elements
  !
  ! Components for momentum, temperature, and salinity interior diffusivities
  type :: interior_diffusivity
     type(interior_diffusivity_elements) :: &
          m, &  ! Elements of interior diffusivity for momentum
          T, &  ! Elements of interior diffusivity for temperature
          S     ! Elements of interior diffusivity for salinity
  end type interior_diffusivity
  !
  ! Elements of turbulent velocity scale profiles
  type :: turbulent_vel_scale_elements
     real(kind=dp) :: &
          h,      &  ! Value at the mixing layer depth
          grad_h     ! Vertical gradient value at the mixing layer
                     ! depth
  end type turbulent_vel_scale_elements
  !
  ! Components for turbulent velocity scale profiles
  type :: turbulent_vel_scales
     type(turbulent_vel_scale_elements) :: &
          m, &  ! Momentum profile
          s     ! Scalar (temperature, salinity) profile
  end type turbulent_vel_scales
  !
  ! Elements of non-dimensional vertical shape functions (G(sigma))
  type :: G_shape_elements
     real(kind=dp) :: &
          h,      &  ! Value at the mixing layer depth
          grad_h, &  ! Vertical gradient value at the mixing layer
                     ! depth
          a2,     &  ! a2 coefficient (Large, et al (1994), eqns 11 & 17)
          a3         ! a3 coefficient (Large, et al (1994), eqns 11 & 17)
  end type G_shape_elements
  !
  ! Components for momentum, temperature, and salinity non
  ! -dimensional vertical shape functions (G(sigma))
  type :: G_shape_functions
     type(G_shape_elements) :: &
          m, &  ! Elements of G(sigma) for momentum
          T, &  ! Elements of G(sigma) for temperature
          S     ! Elements of G(sigma) for salinity
  end type G_shape_functions

  ! Parameter Value Declarations:
  !
  ! Public:
  real(kind=dp), parameter :: &
       zeta_m = -0.20d0, &  ! Max zeta value of the -1/3 power law
                            ! regime of non-dimensional turbulent
                            ! momentum flux profile
       zeta_s = -1.0d0,  &  ! Max zeta value of the -1/3 power law
                            ! regime of non-dimensional turbulent
                            ! scalar flux profile
       a_m = 1.26d0,     &  ! Coefficient of non-dimensional
                            ! turbulent momentum flux profile in 1/3
                            ! power law regime
       a_s = -28.86d0,   &  ! Coefficient of non-dimensional
                            ! turbulent scalar flux profile in 1/3
                            ! power law regime
       c_m = 8.38d0,     &  ! Coefficient of non-dimensional
                            ! turbulent momentum flux profile in 1/3
                            ! power law regime
       c_s = 98.96d0,    &  ! Coefficient of non-dimensional
                            ! turbulent scalar flux profile in 1/3
                            ! power law regime
       kapa = 0.4d0,     &  ! von Karman constant
       epsiln = 0.1d0       ! Non-dimensional extent of the surface layer
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
  type(interior_diffusivity) :: &
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
  real(kind=dp), dimension(2) :: &
       shear_diff_smooth  ! Shear diffusivity smoothing parameters,
                          ! central value and value offset 1
                          ! value offset 2 is calculated
  real(kind=dp) :: &
       u_star, &  ! Turbulent friction velocity
       w_star, &  ! Convective velocity scale
       L_mo       ! Monin-Obukhov length scale
  type(turbulent_vel_scales) :: &
       w  ! Velocity scale profile arrays
  type(G_shape_functions) :: &
       G  ! Non-dimensional vertical shape function coefficients, and
          ! values at the mixing layer depth
          !   Components:
          !     %m  -- momentum components
          !     %T  -- temperature
          !     %S  -- salinity
          !       %h       -- value at mixing layer depth
          !       %grad_h  -- vertical gradient value at mixing layer
          !                   depth
          !       %a2      -- a2 coefficient value
          !       %a3      -- a3 coefficient value
  type(mTS_arrays) :: &
       K_ML  ! Mixing layer diffusivity

contains

  subroutine init_turbulence(M)
    ! Allocate memory for turbulence model variables, and read
    ! parameter values from infile.

    ! Functions from other modules:
    use input_processor, only: getpard, getpardv

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
    ! Shear diffusivity smoothing parameter
    call getpardv('shear smooth', 2, shear_diff_smooth)
  end subroutine init_turbulence


  subroutine calc_KPP_diffusivity(Bf, h, h_i, h_g)
    ! Calculate the diffusivity profile using the K profile
    ! parameterization algorithm of Large, et al (1994).

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters & arrays

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         Bf, &  ! Surface buoyancy forcing
         h      ! Mixing layer depth
    integer, intent(in) :: &
         h_i, &  ! Index of grid layer interface immediate below
                 ! mizing layer depth
         h_g     ! Index of grid point immediate below mizing layer
                 ! depth

    ! Local variables:
    real(kind=dp) :: &
         sigma  ! Non-dimensional depth coordinate within mixing layer
    integer :: &
         j  ! Loop index over profile depth

    ! Step 1: Calculate total turbulent momentum, thermal &
    ! salinity diffusivities for the whole water column based on
    ! shear, double diffusion, and internal wave breaking
    ! instabilities.

    ! Calculate the turbulent momentum diffusivity due to shear
    ! instabilities (nu%m%shear) using the Large, et al (1994)
    ! parameterization based on the local gradient Richardson number.
    call shear_diffusivity()

    ! Calculate turbulent temperature and salinity diffusivities due
    ! to double diffusive mixing using the Large, et al (1994)
    ! formulation that considers both salt fingering, and diffusive
    ! convection.
    call double_diffusion()

    ! Calculate total interior diffusivities: sum of vertical shear,
    ! internal wave breaking, and double diffusion diffusivities
    ! (Large, et al (1994), eq'n (25)).  The zero value at the surfae
    ! (index = 0) imposed the boundary condition of no diffusion
    ! across the air/water interface.
    !
    ! Momentum
    nu%m%total(0) = 0.0d0
    nu%m%total(1:) = nu%m%shear + nu%m%int_wave + nu%S%dd
    ! Temperature
    nu%T%total(0) = 0.0d0
    nu%T%total(1:) = nu%m%shear + nu%T%int_wave + nu%T%dd
    ! Salinity
    nu%S%total(0) = 0.0d0
    nu%S%total(1:) = nu%m%shear + nu%S%int_wave + nu%S%dd

    ! Step 2: Calculate the turbulent momentum, thermal &
    ! salinity diffusivities in the mixing layer

    ! Calculate turbulence scales that characterize the mixing layer:
    ! turbulent friction velocity, convective velocity scale, and
    ! Monin-Obukhov length scale
    u_star = (wbar%u(0) ** 2 + wbar%v(0) ** 2) ** (1.0d0/4.0d0)
    ! the cube root of a negative number is defined but F90 has trouble and
    ! gets a NaN: so do it explicitly
    if (-Bf*h < 0) then
       w_star = - (Bf * h) ** (1.0d0/3.0d0)
    else
       w_star = (-Bf * h) ** (1.0d0/3.0d0)
    endif
    L_mo = u_star ** 3 / (kapa * Bf + epsilon(Bf))

    ! there is only a mixing layer if there is wind or negative buoyancy
    if (abs(u_star) > epsilon(u_star) .or. Bf < epsilon(Bf)) then
       ! Calculate the coefficients of the non-dimension vertical shape
       ! functions (Large, et al (1994), eqn 11).
       call G_shape_parameters(h, h_i, Bf)

       ! Calculate the profiles of turbulent momentum, thermal & salinity
       ! diffusivities in the mixing layer (Large, et al (1994), eqn 10).
       K_ML%m = 0.0d0
       K_ML%T = 0.0d0
       K_ML%S = 0.0d0
       do j = 1, h_i
          if (grid%d_i(j) > h) then
             K_ML%m(j) = 0.0d0
             K_ML%T(j) = 0.0d0
             K_ML%S(j) = 0.0d0
          else
             sigma = grid%d_i(j) / h
             K_ML%m(j) = h &
                  * w_scale(h, grid%d_i(j), Bf, nondim_momentum_flux, c_m) &
                  * G_shape(G%m, sigma)
             K_ML%T(j) = h &
                  * w_scale(h, grid%d_i(j), Bf, nondim_scalar_flux, c_s) &
                  * G_shape(G%T, sigma)
             K_ML%S(j) = h &
                  * w_scale(h, grid%d_i(j), Bf, nondim_scalar_flux, c_s) &
                  * G_shape(G%S, sigma)
          endif
       enddo

       ! Modify the values of the diffusivities at the grid layer
       ! interface just above the mixing layer depth (Large, et al
       ! (1994), eqn D6).  This reduces the mixing layer depth bias
       ! discussed in Large, et al (1994), App. C.  See fig D1 to
       ! understand why we use (h_g - 1) as the index.
       !
       ! Momentum
       K_ML%m(h_g-1) = modify_K(h, h_i, h_g, Bf, nondim_momentum_flux, &
            c_m, G%m, nu%m%total(h_g-1), K_ML%m(h_g-1))
       ! Temperature
       K_ML%T(h_g-1) = modify_K(h, h_i, h_g, Bf, nondim_scalar_flux, &
            c_s, G%T, nu%T%total(h_g-1), K_ML%T(h_g-1))
       ! Salinity
       K_ML%S(h_g-1) = modify_K(h, h_i, h_g, Bf, nondim_scalar_flux, &
            c_s, G%S, nu%S%total(h_g-1), K_ML%S(h_g-1))
    else ! no mixing layer
       K_ML%m = 0.0d0
       K_ML%T = 0.0d0
       K_ML%S = 0.0d0
    endif

    ! Calculate the overall profiles of turbulent momentum, thermal & salinity
    ! diffusivities in the water column
    do j = 1, grid%M
       if (j <= h_i .and. grid%d_i(j) <= h) then
          ! Use the larger of the mixing layer diffusivity and the
          ! shear diffusivity.  This handles the situation of, for
          ! instance, the wind dying so the surface forcing produces
          ! minimal mixing layer diffusivity, but significant shear
          ! remains in the mixing layer.
          K%m(j) = max(K_ML%m(j), nu%m%total(j))
          K%T(j) = max(K_ML%T(j), nu%T%total(j))
          K%S(j) = max(K_ML%S(j), nu%S%total(j))
       else
          K%m(j) = nu%m%total(j)
          K%T(j) = nu%T%total(j)
          K%S(j) = nu%S%total(j)
       endif
    enddo
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
    ! Variable Declarations:
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
         Ri_o = 0.7d0,  &  ! Critical gradient Richardson number
         nu_o = 5.0d-3, &  ! Maximum shear diffusivity
         p_1 = 3.0d0       ! Power constant for shear diffusivity
                           ! parameterization
    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: &
         Ri_g, &   ! Profile of Richardson number gradient at grid
                   ! layer interface depths
         N2,   &   ! Profile of buoyancy frequency squared at grid
                   ! layer interface depths
         V2,   &   ! Profile of the square of the magnitude of the
                   ! velocity gradient at the grid layer interface
                   ! depths
         smoothed  ! Smoothed values of shear diffusivity
    real(kind=dp) :: &
         a1, a2, a3, a4, a5  ! Smoothing algorithm parameters
    integer :: &
         i  ! Index over profile depth

    ! Calculate the profiles of the buoyancy frequency squared,
    ! velocity gradient squared, and Richardson number gradient at the
    ! grid layer interface depths.  (Large, et al (1994), eq'n (27))
    N2 = (-g / rho%i(1:grid%M)) * rho%grad_i(1:grid%M)
    V2 = U%grad_i(1:grid%M) ** 2 + V%grad_i(1:grid%M) ** 2
    Ri_g = N2 / (V2 + epsilon(V2))
    ! Apply the shear diffusivity parameterization (Large, et al
    ! (1994), eq'n(28))
    do i = 1, grid%M
       if (Ri_g(i) <= 0.0d0) then
          ! Eq'n (28a)
          nu%m%shear(i) = nu_o
       elseif (0.0d0 < Ri_g(i) .and. Ri_g(i) < Ri_o) then
          ! Eq'n (28b)
          nu%m%shear(i) = nu_o * (1.0d0 - (Ri_g(i) / Ri_o) ** 2) ** p_1
       else
          ! Eq'n (28c)
          nu%m%shear(i) = 0.0d0
       endif
    enddo
    ! Smooth the shear diffusivity over adjacent grid layer interfaces
    ! because the estimation of shear diffusivity is noisy as it is
    ! the ratio of a pair of differences.
    !
    ! Calculate the coefficients:
    ! central points
    a1 = shear_diff_smooth(1)
    a2 = shear_diff_smooth(2) / 2.0d0
    a3 = (1.0d0 - a1 - a2 * 2) / 2.0d0
    ! points one from boundary
    a4 = a1 + 2 * a2 + a3
    ! boundary points
    a5 = a1 + a2 + a3
    !
    ! Calculate the smoothed sher diffusivity profile
    smoothed(1) = (a1 * nu%m%shear(1) + a2 * nu%m%shear(2) &
                  + a3 * nu%m%shear(3)) / a5
    smoothed(2) = (a2 * nu%m%shear(1) + a1 * nu%m%shear(2)       &
                  + a2 * nu%m%shear(3) + a3 * nu%m%shear(4)) / a4
    smoothed(3:grid%M-2) = a3 * nu%m%shear(1 : grid%M - 4)  &
         + a2 * nu%m%shear(2 : grid%M - 3)                  &
         + a1 * nu%m%shear(3 : grid%M - 2)                  &
         + a2 * nu%m%shear(4 : grid%M - 1)                  &
         + a3 * nu%m%shear(5 : grid%M)
    smoothed(grid%M-1) = (a3 * nu%m%shear(grid%M - 3)   &
                         + a2 * nu%m%shear(grid%M - 2)  &
                         + a1 * nu%m%shear(grid%M - 1)  &
                         + a2 * nu%m%shear(grid%M)) / a4
    smoothed(grid%M) = (a3 * nu%m%shear(grid%M - 2)    &
                        + a2 * nu%m%shear(grid%M - 1)  &
                        + a1 * nu%m%shear(grid%M)) / a5
    nu%m%shear = smoothed
  end subroutine shear_diffusivity


  subroutine double_diffusion()
    ! Calculate turbulent temperature and salinity diffusivities
    ! (nu%*%dd) due to double diffusive mixing using the Large, et al
    ! (1994) formulation that considers both salt fingering, and
    ! diffusive convection.  (Large, etal (1994), pp373-374, equations
    ! 30 to 34).

    ! Elements from other modules:
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variable Declarations:
    use grid_mod, only: &
         grid  ! Grid parameters, and arrays
    use core_variables, only: &
         T, &  ! Temperature profile arrays; we need the gradients at
               ! the grid layer interface depths: T%grad_i
         S     ! Salinity profile arrays; we need the gradients at the
               ! grid layer interface depths: S%grad_i
    use water_properties, only: &
         alpha, &  ! Thermal expansion coefficient profile arrays; we
                   ! need the gradient at the grid layer interface
                   ! depths: alpha%grad_i
         beta      ! Salinity contraction coefficient profile arrays; we
                   ! need the gradient at the grid layer interface
                   ! depths: beta%grad_i

    implicit none

    ! Local parameter value declarations:
    real(kind=dp), parameter :: &
         nu_f = 10.0d-4,  &  ! Maximum salt fingering diffusivity [m^2/s]
         R_rho_o = 1.9d0, &  ! Salt fingering limit (nu_s = 0)
         p_2 = 3,         &  ! Salt finger diff power constant
         nu_m = 1.5d-06      ! Molecular diffusivity (viscosity) [m^2/s]
    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: &
         R_rho  ! Double diffusion density ratio
    integer :: &
         j  ! Loop index over depth

    ! Calculate double diffusion density ratio (Large, et al, (1994) eq'n 30)
    R_rho = alpha%i(1:) * T%grad_i / (beta%i(1:) * S%grad_i + epsilon(S%grad_i(1)))

    do j = 1, grid%M
       ! Determine if there is double diffision, and if so, what type,
       ! salt fingering, or diffusive convection
       if (1.0d0 < R_rho(j) &
            .and. R_rho(j) < 2.0d0 &
            .and. alpha%i(j) * T%grad_i(j) > 0.0d0 &
            .and. beta%i(j) * S%grad_i(j) > 0.0d0) then
          ! Salt fingering
          if (1.0d0 < R_rho(j) .and. R_rho(j) < R_rho_o) then
             ! Large, et al, (1994) eq'n 31a
             nu%S%dd(j) = nu_f &
                  * (1.0d0 - ((R_rho(j) - 1.0d0) &
                             / (R_rho_o - 1.0d0)) ** 2) ** p_2
          else  ! R_rho >= R_rho_o
             ! Large, et al, (1994) eq'n 31b
             nu%S%dd = 0.
          endif
          ! Large, et al, (1994) eq'n 31c
          nu%T%dd(j) = 0.7 * nu%S%dd(j)
       else if (0. < R_rho(j) &
            .and. R_rho(j) < 1.0d0 &
            .and. alpha%i(j) * T%grad_i(j) < 0.0d0 &
            .and. beta%i(j) * S%grad_i(j) < 0.0d0) then
          ! Diffusive convection
          ! Large, et al, (1994) eq'n 32
          nu%T%dd(j) = nu_m * 0.909d0 &
               * exp(4.6d0 * exp(-0.54d0 * (1.0d0 / R_rho(j) - 1.0d0)))
          ! Large, et al, (1994) eq'n 34
          if (R_rho(j) >= 0.5d0 .and. R_rho(j) < 1.0d0) then
             nu%S%dd(j) = nu%T%dd(j) * (1.85d0 - 0.85d0 / R_rho(j)) * R_rho(j)
          else if (R_rho(j) < 0.5d0) then
             nu%S%dd(j) = nu%T%dd(j) * 0.15d0 * R_rho(j)
          endif
       else
          ! No double diffusion (i.e. temperature and salinity are
          ! both stabilizing, or water column is convectively
          ! unstable)
          nu%T%dd(j) = 0.0d0
          nu%S%dd(j) = 0.0d0
       endif
    enddo
  end subroutine double_diffusion


  function w_scale(h, d, Bf, nondim_flux, c) result(w_value)
    ! Calculate the vertical turbulent velocity scale value at the
    ! specified depth.


    ! Elements from other modules:
    ! Type Definitions:
    use precision_defs, only: dp

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         h,  &  ! Mixing layer extent [m]
         d,  &  ! Depth [m]
         Bf, &  ! Surface buoyancy forcing
         c      ! Coefficient of non-dimensional turbulent flux profile
                ! in 1/3 power law regime
    real(kind=dp) :: &
         nondim_flux  ! Non-dimensional flux function appropriate to
                      ! the scale (momentum or scalar) that we are
                      ! calculating.

    ! Result:
    real(kind=dp) :: &
         w_value  ! Value of the turbulent vertical velocity scale at
                  ! the specified depth.

    ! Local variables:
    real(kind=dp) :: &
         d_surf, &  ! Surface layer extent [m]
         zeta,   &  ! Stability parameter; ratio of depth to
                    ! Monin-Obukhov length scale
         sigma      ! Non-dimensional vertical coordinate in the
                    ! mixing layer

    ! Calculate extent of surface layer
    d_surf = epsiln * h

    ! Note that abs(x) > epsilon(x) is a real-number-robust test for x
    ! /= 0, and abs(x) < epsilon(x) is similarly for x == 0.
    if (abs(u_star) > epsilon(u_star)) then
       ! Wind forced mixing layer
       !
       ! Calculate the stability parameter value at the jth grid
       ! layer interface depth
       zeta = d / L_mo
       ! Calculate the turbulent velocity scale value using Large, et
       ! al (1994), eqn (13) that in turn uses the non-dimensional
       ! flux profiles (Large, et al (1994), App. B).
       if (d_surf < d .and. d < h .and. zeta < 0.0d0) then
          ! Special case of depths between the surface layer depth
          ! and the mixing layer depth, in a stable mixing layer
          w_value = kapa * u_star / nondim_flux(d_surf)
       else
          ! General value of the turbulent velocity scale
          w_value = kapa * u_star / nondim_flux(d)
       endif
    elseif (abs(u_star) < epsilon(u_star) .and. Bf < 0.0d0) then
       ! Special case of the limit of convectively unstable mixing
       ! layer in the absence of wind forcing (Large, et al (1994),
       ! eqn 15)
       if (d < d_surf) then
          ! In the surface layer
          sigma = d / h
       elseif (d_surf <= d ) then
          ! Between the surface layer and mixing layer depths
          ! and below the mixing layer as scales in unstable circumstances
          ! should continue (below Eqn 13)
          sigma = epsiln
       endif
       ! Calculate the turbulent velocity scale value
       w_value = kapa * (c * kapa * sigma) ** (1.0d0 / 3.0d0) * w_star
    else
       ! No wind forcing, and no convective forcing, so no turbulence
       ! in the mixing layer
       w_value = 0.0d0
    endif
  end function w_scale

  function w_scale_special(h, d, c) result(w_value)
    ! Calculate the vertical turbulent velocity scale value at the
    ! specified depth for use in the bottom of (20) when no wind
    ! note that abs(u_star) < epsilon(u_star) .and. Bf < 0.0d0

    ! Elements from other modules:
    ! Type Definitions:
    use precision_defs, only: dp
    use water_properties, only: alpha

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         h,  &  ! Mixing layer extent [m]
         d,  &  ! Depth [m]
         c      ! Coefficient of non-dimensional turbulent flux profile
                ! in 1/3 power law regime

    ! Result:
    real(kind=dp) :: &
         w_value  ! Value of the turbulent vertical velocity scale at
                  ! the specified depth.

    ! Local variables:
    real(kind=dp) :: &
         d_surf, &  ! Surface layer extent [m]
         wspecial, &! a special wstar based on heat flux only
         sigma      ! Non-dimensional vertical coordinate in the
                    ! mixing layer

    ! Calculate extent of surface layer
    d_surf = epsiln * h

    ! Special case of the limit of convectively unstable mixing
    ! layer in the absence of wind forcing (Large, et al (1994),
    ! eqn 15)
    if (d < d_surf) then
       ! In the surface layer
       sigma = d / h
    elseif (d_surf <= d ) then
       ! Between the surface layer and mixing layer depths
       ! and below the mixing layer as scales in unstable circumstances
       ! should continue (below Eqn 13)
       sigma = epsiln
    endif
    ! Calculate a special "wstar" based on wbar%t(0) not Bf
    wspecial = (alpha%g(0) * wbar%t(0) * h) ** (1.0d0 / 3.0d0)
    ! Calculate the turbulent velocity scale value
    w_value = kapa * (c * kapa * sigma) ** (1.0d0 / 3.0d0) * wspecial

  end function w_scale_special


  function nondim_momentum_flux(d) result(phi_m)
    ! Return the value of the non-dimensional momentum flux profile at
    ! the specified depth (Large, etal (1994), App. B)

    ! Type Definitions from Other Modules:
    use precision_defs, only: dp
    ! Variable Declarations:

    implicit none

    ! Argument:
    real(kind=dp), intent(in) :: &
         d  ! Depth [m]
    ! Result:
    real(kind=dp) :: &
         phi_m  ! Value of the non-dimensional momentum flux profile
                ! at the specified depth

    ! Local variables:
    real(kind=dp) :: &
         zeta  ! Stability parameter; ratio of depth to Monin-Obukhov
               ! length scale

    ! Calculate the stability parameter value at the specified depth
    zeta = d / (L_mo + epsilon(L_mo))
    if (0.0d0 <= zeta) then
       ! Stable (eqn B1a)
       phi_m = 1.0d0 + 5.0d0 * zeta
    endif
    ! Unstable
    if (zeta_m <= zeta .and. zeta < 0.) then
       ! Eqn B1b
       phi_m = (1.0d0 - 16.0d0 * zeta) ** (-1.0d0 / 4.0d0)
    elseif (zeta < zeta_m) then
       ! Eqn B1c
       phi_m = (a_m - c_m * zeta) ** (-1.0d0 / 3.0d0)
    endif
  end function nondim_momentum_flux


  function nondim_scalar_flux(d) result(phi_s)
    ! Return the value of the non-dimensional scalar flux profile at
    ! the specified depth (Large, etal (1994), App. B).

    ! Type Definitions from Other Modules:
    use precision_defs, only: dp

    implicit none

    ! Argument:
    real(kind=dp), intent(in) :: &
         d  ! Depth [m]
    ! Result:
    real(kind=dp) :: &
         phi_s  ! Value of the non-dimensional scalar flux profile
                ! at the specified depth

    ! Local variables:
    real(kind=dp) :: &
         zeta  ! Stability parameter; ratio of depth to Monin-Obukhov
               ! length scale

    ! Calculate the stability parameter value at the specified depth
    zeta = d / (L_mo + epsilon(L_mo))
    if (0.0d0 <= zeta) then
       ! Stable (eqn B1a)
       phi_s = 1.0d0 + 5.0d0 * zeta
    endif
    ! Unstable
    if (zeta_s <= zeta .and. zeta < 0.0d0) then
       ! Eqn B1b
       phi_s = (1.0d0 - 16.0d0 * zeta) ** (-1.0d0 / 2.0d0)
    elseif (zeta < zeta_s) then
       ! Eqn B1c
       phi_s = (a_s - c_s * zeta) ** (-1.0d0 / 3.0d0)
    endif
  end function nondim_scalar_flux


  subroutine G_shape_parameters(h, h_i, Bf)
    ! Calculate the coefficients of the non-dimensional vertical shape
    ! functions (Large, et al (1994), eqn 11).

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         h, &  ! Mixing layer depth
         Bf    ! Surface buoyancy forcing
    integer, intent(in) :: &
         h_i  ! Index of grid layer interface immediate below mizing
              ! layer depth

    ! Local Variables:
    real(kind=dp) :: &
         zeta_h  ! Value of stability parameter at mixing layer depth

    ! Calculate the values of the total interior diffusivities, and
    ! their vertical gradients at the mixing layer depth.  These
    ! values facilitate a continuous, smooth match between the
    ! interior diffusivity and that of the the mixing layer.  They are
    ! used to calculate the values of the mixing layer diffusivity
    ! shape functions at the mixing layer depth, and the coefficients
    ! of those shape function.  (Large, et al (1994), eqns 17 & 18 &
    ! App. D)
    !
    ! Momentum
    call nu_h(h, h_i, nu%m)
    ! Temperature
    call nu_h(h, h_i, nu%T)
    ! Salinity
    call nu_h(h, h_i, nu%S)

    ! Set the values of the turbulent velocity scales, and their
    ! vertical gradients at the mixing layer depth.
    !
    ! Momentum
    w%m%h = w_scale(h, h, Bf, nondim_momentum_flux, c_m)
    ! Calculate the stability parameter value at the mixing layer depth
    if (abs(L_mo) > epsilon(L_mo)) then
       zeta_h = h / L_mo
    else ! this gets the sign right, the value is infinite
       zeta_h = Bf
    endif
    if (zeta_h <= 0.d0) then
       ! Unstable and neutral forcing; gradients are zero (Large, et
       ! al (1994), below eqn 18)
       w%s%h = w_scale(h, h, Bf, nondim_scalar_flux, c_m)
       w%m%grad_h = 0.d0
       w%s%grad_h = 0.d0
    else
       ! Stable forcing; momentum and scalar value are equal (Large,
       ! et al (1994), eqn B1a)
       w%s%h = w%m%h
       ! d/d{sigma} of Eqn B1a substitutes into stable forcing form of
       ! eqn 13
       w%m%grad_h = -5.0d0 * kapa * u_star * zeta_h &
            / (1.0d0 + 5.0d0 * zeta_h) ** 2
       w%s%grad_h = w%m%grad_h
    endif

    ! Calculate the coefficients of the non-dimensional vertical shape
    ! functions (G(sigma)) for the mixing layer diffusivity (Large, et
    ! al (1994), eqn 17)
    !
    ! Momentum
    call G_coefficients(h, nu%m, w%m, G%m)
    ! Temperature
    call G_coefficients(h, nu%T, w%s, G%T)
    ! Salinity
    call G_coefficients(h, nu%S, w%s, G%S)
  end subroutine G_shape_parameters


  subroutine nu_h(h, h_i, nu)
    ! Calculate the values of the total interior diffusivities, and
    ! their vertical gradients at the mixing layer depth.  These
    ! values facilitate a continuous, smooth match between the
    ! interior diffusivity and that of the the mixing layer.  They are
    ! used to calculate the values of the mixing layer diffusivity
    ! shape functions at the mixing layer depth, and the coefficients
    ! of those shape function.  (Large, et al (1994), eqns 17 & 18 &
    ! App. D)
    !
    ! Note that the 2 situations: (1) mixing layer depth above grid
    ! layer interface (d_{k-1} < h < d_{k-1/2}), and (2) mixing layer
    ! depth below grid layer interface (d_{k-1/2} < h d_{k}) are both
    ! handled by the same equations because h_i is always the index of
    ! the grid layer interface immediately below the mixing layer
    ! depth.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters & arrays

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         h      ! Mixing layer depth
    integer, intent(in) :: &
         h_i  ! Index of grid layer interface immediate below mizing
              ! layer depth
    type(interior_diffusivity_elements), intent(inout) :: &
         nu  ! Interior diffusivity

    ! Local Variables:
    real(kind=dp) :: &
         R, Rdel_n, Rdel_np1  ! Weighting factor for calculation of
                              ! values of the total interior
                              ! diffusivities, and their vertical
                              ! gradients at the mixing layer depth.

    ! Calculate the weighting factor (Large, et al (1994), below eqn
    ! D5)
    R = (h - grid%d_i(h_i-1)) / grid%i_space(h_i)
    ! Common factors for gradient
    Rdel_n = (1.0d0 - R) / grid%i_space(h_i)
    Rdel_np1 = R / grid%i_space(h_i+1)
    ! Diffusivity gradients and interpolated values at the mixing
    ! layer depth (Large, et al (1994), eqn D5)
    nu%tot_grad_h = Rdel_n * (nu%total(h_i-1) - nu%total(h_i)) &
         + Rdel_np1 * (nu%total(h_i) - nu%total(h_i+1))
    nu%tot_h = nu%total(h_i) + nu%tot_grad_h * (grid%d_i(h_i) - h)
  end subroutine nu_h


  subroutine G_coefficients(h, nu, w, G)
    ! Calculate the values of the a2 and a3 coefficients of the
    ! non-dimensional vertical shape function (G(sigma)) (Large, et al
    ! (1994), eqn 17).

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         h      ! Mixing layer depth
    type(interior_diffusivity_elements), intent(in) :: &
         nu  ! Interior diffusivity
    type(turbulent_vel_scale_elements), intent(in) :: &
         w  ! Turbulent velocity scale profile
    type(G_shape_elements), intent(out) :: &
         G  ! Shape function coefficients

    ! Calculate mixing layer shape function value, and its gradient at
    ! mixing layer depth (eqn 18)
    G%h = nu%tot_h / (h * w%h + epsilon(h))
    G%grad_h = -nu%tot_grad_h / w%h - nu%tot_h * w%grad_h &
         / (h * w%h ** 2 + epsilon(h))
    if (G%h < 0.0d0) then
       ! Diffusivity must be positive
       G%h = 0.0d0
       G%grad_h = 0.0d0
    elseif (G%grad_h > 0.0d0) then
       ! Diffusivity gradient must be negative; i.e. interior can only
       ! contribute to increasing the diffusivity in the mixing layer
       G%grad_h = 0.0d0
    endif
    ! Calculate the shape function coefficients (eqn 17)
    G%a2 = -2.0d0 + 3.0d0 * G%h - G%grad_h
    G%a3 =  1.0d0 - 2.0d0 * G%h + G%grad_h
  end subroutine G_coefficients


  function G_shape(G, sigma) result(G_sigma)
    ! Return the value of the non-dimensional vertical shape function
    ! at the specified value of sigma, the non-dimensional coordinate
    ! in the mixing layer.

    ! Type definitions from other modules:
    use precision_defs, only: dp

    implicit none

    ! Arguments:
    type(G_shape_elements), intent(in) :: &
         G  ! Shape function coefficients
    real(kind=dp), intent(in) :: &
         sigma  ! Value of non-dimensional coordinate in the mixing
                ! layer to evaluate G at

    ! Results:
    real(kind=dp) :: &
         G_sigma  !  Value of value of the non-dimensional vertical
                  !  shape function at the specified value of sigma

    G_sigma = sigma + G%a2 * sigma ** 2 + G%a3 * sigma ** 3
  end function G_shape


  function modify_K(h, h_i, h_g, Bf, nondim_flux, c, G, nu_h_gm1, K_ML_h_gm1) &
       result(Lambda)
    ! Modify the value of the diffusivity at the grid layer interface
    ! just above the mixing layer depth (Large, et al (1994), eqn D6).
    ! This reduces the mixing layer depth bias discussed in Large, et
    ! al (1994), App. C.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters & arrays

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         h   ! Mixing layer depth
    integer, intent(in) :: &
         h_i, &  ! Index of grid layer interface immediately below
                 ! mixing layer depth.
         h_g     ! Index of grid point immediately below mixing layer
                 ! depth.
    real(kind=dp), intent(in) :: &
         Bf     ! Surface buoyancy flux
    real(kind=dp), external :: &
         nondim_flux  ! Non-dimensional flux function appropriate to
                      ! the quantity (momentum or scalar) that we are
                      ! modifying the diffusivity of.
    real(kind=dp), intent(in) :: &
         c      ! Coefficient of non-dimensional turbulent flux profile
                ! in 1/3 power law regime
    type(G_shape_elements), intent(in) :: &
         G  ! Shape function coefficients
    real(kind=dp), intent(in) :: &
         nu_h_gm1, &  ! Total interior diffusivity at the grid layer
                      ! interface depth immediately above h%g.
         K_ML_h_gm1   ! Mixing layer diffusivity at the grid layer
                      ! interface depth immediately above h%g.

    ! Result:
    real(kind=dp) :: &
         Lambda  ! Modified diffusivity value at the grid layer
                 ! interface just above the mixing layer depth.

    ! Local variables:
    real(kind=dp) :: &
         del, &       ! Relative location of the mixing layer depth
                      ! withing the grid layer
         K_ML_gm1, &  ! Value of the mixing layer diffusivity at the
                      ! grid point depth immediately above the mixing
                      ! layer depth
         K_star       ! Enhanced diffusivity

    ! Calculate the relative location within the grid layer of the
    ! mixing layer depth (Large, et al (1994), eqn D2).
    del = (h - grid%d_g(h_g-1)) / grid%g_space(h_g-1)
    ! Calculate the value of the mixing layer diffusivity at the grid
    ! point immediately above the mixing layer depth (Large, et al
    ! (1994), eqn D6).
    K_ML_gm1 = h &
         * w_scale(h, grid%d_g(h_g-1), Bf, nondim_flux, c) &
         * G_shape(G, 1.0d0)
    ! Calculate the enhanced diffusivity depending on the location of
    ! the mixing layer depth within the grid layer; 2 situations:
    if (grid%d_g(h_g-1) < h .and. h <= grid%d_i(h_i-1)) then
       ! Between the grid point above, and the grid layer interface
       K_star = K_ML_gm1 * (1.0d0 - del) ** 2 + nu_h_gm1 * del ** 2
    elseif (grid%d_i(h_i-1) < h .and. h <= grid%d_g(h_g)) then
       ! Between the grid layer interface above, and the grid point
       K_star = K_ML_gm1 * (1.0d0 - del) ** 2 + K_ML_h_gm1 * del ** 2
    endif
    ! Calculate the modified diffusivity (Large, et al (1994), eqn
    ! D6).
    Lambda = (1.0d0 - del) * nu_h_gm1 + del * K_star
  end function modify_K


  function nonlocal_scalar_transport(flux, h, d_i, Bf) result(gamma_value)
    ! Return the value of the non-local transport term for the
    ! specified flux, and depth.  (Large, et al (1994), eqn 20).

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Functions:

    ! Arguments:
    real(kind=dp), intent(in) :: &
         flux, &  ! Flux to calculate the non-local transport term
                  ! value for.
         h,    &  ! Mixing layer depth
         d_i,  &  ! Depth of grid layer interface to calculate the
                  ! non-local transport term value at.
         Bf       ! Surface bouyancy flux

    ! Result:
    real(kind=dp) :: &
         gamma_value  ! Value of the non-local transport term for the
                      ! specified flux, and depth.

    ! Local Parameters:
    real(kind=dp), parameter :: &
         C_star = 9.9d0  ! Proportionality coefficient for
                         ! parameterizatio of non-local transport
                         ! terms.  *** Large, et al (1994) recommends
                         ! a value of 10.

    !Local Variables:
    real(kind=dp) :: &
         w_S      ! Value of the scalar turbulent velocity scale at
                  ! the depth for which the non-local transport term
                  ! is being calculated.

    ! Note that abs(x) > epsilon(x) is a real-number-robust test for
    ! x /= 0.

!    if ((abs(u_star) > epsilon(u_star) .or. abs(Bf) > epsilon(Bf)) &
!         .and. (Bf < 0.0d0) &
!         .and. d_i < h) then
    if ( Bf < 0.0d0  .and. d_i < h) then
       ! it is unstable and we're in the mixed layer
       if (abs(u_star) > epsilon(u_star)) then
          ! there is wind forcing
          ! Calculate Value of the scalar turbulent velocity scale at the
          ! depth for which the non-local transport term is being
          ! calculated.
          w_S = w_scale(h, d_i, Bf, nondim_scalar_flux, c_s)
          ! Calculate value of the non-local transport term for the
          ! specified flux, and depth,
          gamma_value = C_star * kapa                     &
               * (c_s * kapa * epsiln) ** (1.0d0 / 3.0d0) &
               * flux / (w_S * h + epsilon(h))
       elseif (abs(Bf) > epsilon(Bf)) then
          ! there is no wind forcing
          ! because, typically the salinity is stabilizing and the temp
          ! is destabilizing, if we use (20) from Large we can get
          ! enormous forcing if Bf is close to zero because w* is in the
          ! denominator.  Instead make gamma propto wbar%t(0) to the 2/3
          ! K remains propto Bf**1/3
          w_S = w_scale_special(h, d_i, c_s)
          ! Calculate value of the non-local transport term for the
          ! specified flux, and depth,
          gamma_value = C_star * kapa                     &
               * (c_s * kapa * epsiln) ** (1.0d0 / 3.0d0) &
               * flux / (w_S * h + epsilon(h))
       else
          ! No forcing
          gamma_value = 0.0d0
       endif
    else
       ! it's stable, or we're below the mixing layer.
       gamma_value = 0.0d0
    endif
  end function nonlocal_scalar_transport


  subroutine calc_turbulent_fluxes(h, h_i, Bf)
    ! Calculate the turbulent kinematic flux profiles.
    !
    ! Note that only the surface values (index = 0) are used in the
    ! model.  This subroutine calculates the full profiles.  It is not
    ! called, but can be included in a build if the profiles are
    ! required for diagnostic purposes.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Parameters:
    use fundamental_constants, only: g
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and arrays
    use core_variables, only: &
         U,  &  ! Cross-strait (35 deg) velocity component arrays
         V,  &  ! Along-strait (305 deg) velocity component arrays
         T,  &  ! Temperature profile arrays
         S      ! Salinity profile arrays
    use water_properties, only: &
         alpha, &  ! Thermal expansion coefficient profile arrays
         beta      ! Salinity contraction coefficient profile arrays

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         h, &   ! Mixing layer depth
         Bf     ! Surface bouyancy flux
    integer, intent(in) :: &
         h_i  ! Index of grid layer interface immediately below mixing
              ! layer depth.

    ! Local variable:
    integer :: &
         j  ! Index over grid depth

    do j = 1, grid%M
       if (j <= h_i .and. grid%d_i(j) <= h) then
          ! Large, et al (1994), eqn 9
          !
          ! Momentum (non-local transport term is zero by definition;
          ! see Large, et al (1994), eqn 20)
          wbar%u(j) = -K_ML%m(j) * (U%grad_i(j) - 0.0d0)
          wbar%v(j) = -K_ML%m(j) * (V%grad_i(j) - 0.0d0)
          ! Temperature.  Note that the \bar{w\theta}_R term in Large,
          ! et al (1994), eqn 20 is not included here because it is
          ! zero.  See Large, et al (1994), pg. 379, and App. A for
          ! the explanation.
          wbar%t(j) = -K_ML%T(j) * (T%grad_i(j) &
               - nonlocal_scalar_transport(wbar%t(0), h, grid%d_i(j), Bf))
          wbar%s(j) = -K_ML%S(j) * (S%grad_i(j) &
               - nonlocal_scalar_transport(wbar%s(0), h, grid%d_i(j), Bf))
          ! Large, et al (1994), eqn A3b
          wbar%b(j) = g * (alpha%i(j) * wbar%t(j) - beta%i(j) * wbar%s(j))
          ! *** Buoyancy flux variation due to error in z ???
          wbar%b_err(j) = g * (alpha%grad_i(j) * wbar%t(j) &
               - beta%grad_i(j) * wbar%s(j))
       else
          wbar%u(j) = 0.0d0
          wbar%v(j) = 0.0d0
          wbar%t(j) = 0.0d0
          wbar%s(j) = 0.0d0
          wbar%b(j) = 0.0d0
          wbar%b_err(j) = 0.0d0
       endif
    enddo
  end subroutine calc_turbulent_fluxes


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
    ! nu%*%total is 0-based because its surface value (zero by
    ! definition) is required in subroutine nu_h() when the mixing
    ! layer depth is very shallow
    allocate(nu%m%shear(1:M), nu%T%dd(1:M), nu%S%dd(1:M),   &
         nu%m%total(0:M), nu%T%total(0:M), nu%S%total(0:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Mixing layer diffusivity profile arrays"
    ! K_ML%* is 0-based because its surface value (zero by definition)
    ! is required in subroutine modify_K() when the mixing layer depth
    ! is very shallow
    allocate(K_ML%m(0:M), K_ML%T(0:M), K_ML%S(0:M), &
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
