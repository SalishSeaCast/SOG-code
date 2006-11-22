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
       nK, &      ! Overall diffusivity profile; a continuous profile of
                 ! K_ML%* in the mixing layer, and K%*%total below it
       wbar, &   ! Turbulent kinematic flux profiles
       ! *** Temporary until turbulence refactoring is completed
       nu, &  ! Interior diffusivity profile arrays
       u_star, &  ! Turbulent friction velocity
       L_mo, &    ! Monin-Obukhov length scale
       w, &  ! Turbulent velocity scale profile arrays
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
  !
  ! Components for turbulent velocity scale profiles, and the related
  ! non-dimensional flux profiles
  type :: turbulent_vel_scales
     real(kind=dp), dimension(:), allocatable :: &
          m, &  ! Momentum profile
          s     ! Scalar (temperature, salinity) profile
  end type turbulent_vel_scales

  ! Parameter Value Declarations:
  !
  ! Public:
  real(kind=dp), parameter :: &
       zeta_m = -0.20, &  ! Max zeta value of the -1/3 power law
                          ! regime of non-dimensional turbulent
                          ! momentum flux profile
       zeta_s = -1.0,  &  ! Max zeta value of the -1/3 power law
                          ! regime of non-dimensional turbulent scalar
                          ! flux profile
       a_m = 1.26,    &   ! Coefficient of non-dimensional turbulent
                          ! momentum flux profile in 1/3 power law
                          ! regime
       a_s = -28.86,  &   ! Coefficient of non-dimensional turbulent
                          ! scalar flux profile in 1/3 power law
                          ! regime
       c_m = 8.38,    &   ! Coefficient of non-dimensional turbulent
                          ! momentum flux profile in 1/3 power law
                          ! regime
       c_s = 98.96,   &   ! Coefficient of non-dimensional turbulent
                          ! scalar flux profile in 1/3 power law
                          ! regime
       kapa = 0.4,    &   ! von Karman constant
       epsiln = 0.1       ! Non-dimensional extent of the surface layer
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
  real(kind=dp) :: &
       u_star, &  ! Turbulent friction velocity
       w_star, &  ! Convective velocity scale
       L_mo       ! Monin-Obukhov length scale
  type(turbulent_vel_scales) :: &
       w  ! velocity scale profile arrays
  type(mTS_arrays) :: &
       K_ML  ! Mixing layer diffusivity

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

    use precision_defs, only: dp

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         Bf, &  ! Surface buoyancy forcing
         h      ! Mixing layer depth

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
    ! (Large, et al (1994), eq'n (25))
    !
    ! Momentum
    nu%m%total(1:) = nu%m%shear(1:) + nu%m%int_wave + nu%S%dd(1:)
    ! Temperature
    nu%T%total(1:) = nu%m%shear(1:) + nu%T%int_wave + nu%T%dd(1:)
    ! Salinity
    nu%S%total(1:) = nu%m%shear(1:) + nu%S%int_wave + nu%S%dd(1:)

    ! Calculate the value of the total interior diffusivities, and
    ! their vertical gradients at the mixing layer depth
    call nu_h()

    ! Step 2: Calculate the turbulent momentum, thermal &
    ! salinity diffusivities in the mixing layer

    ! Calculate turbulence scales that characterize the mixing layer:
    ! turbulent friction velocity, convective velocity scale, and
    ! Monin-Obukhov length scale
    u_star = (wbar%u(0) ** 2 + wbar%v(0) ** 2) ** (1./4.)
    w_star = (-Bf * h) ** (1./3.)
    L_mo = u_star ** 3 / (kapa * Bf)

    ! Calculate the vertical turbulent velocity scale profiles (w%*)
    if (u_star /= 0.) then
       ! The mixing layer turbulence is driven by wind forcing, so
       ! calculate the velocity scales using Large, et al (1994), eqn
       ! (13) that in turn uses the non-dimensional flux profiles
       ! (Large, et al (1994), App. B),
       call wind_driven_turbulence(h)
    elseif (u_star == 0. .and. Bf < 0.) then
       ! The mixing layer turbulence is driven by convection, so
       ! calculate the velocity scales in their convective limits
       ! (Large, et al (1994), eqn 15)
       call convection_driven_turbulence(h)
    else
       ! No wind forcing, and no convective forcing, so no turbulence
       ! in the mixing layer
       w%m = 0.
       w%s = 0.
    endif
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
         nu_f = 10.0d-4, &  ! Maximum salt fingering diffusivity [m^2/s]
         R_rho_o = 1.9,  &  ! Salt fingering limit (nu_s = 0)
         p_2 = 3,        &  ! Salt finger diff power constant
         nu_m = 1.5d-06     ! Molecular diffusivity (viscosity) [m^2/s]
    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: &
         R_rho  ! Double diffusion density ratio
    integer :: &
         k  ! Loop index over depth

    ! Calculate double diffusion density ratio (Large, et al, (1994) eq'n 30)
    R_rho = alpha%i(1:) * T%grad_i / (beta%i(1:) * S%grad_i)
      
    do k = 1, grid%M
       ! Determine if there is double diffision, and if so, what type,
       ! salt fingering, or diffusive convection
       if (1.0 < R_rho(k) &
            .and. R_rho(k) < 2.0 &
            .and. alpha%i(k) * T%grad_i(k) > 0 &
            .and. beta%i(k) * S%grad_i(k) > 0.) then  
          ! Salt fingering
          if (1.0 < R_rho(k) .and. R_rho(k) < R_rho_o) then
             ! Large, et al, (1994) eq'n 31a
             nu%S%dd(k) = nu_f &
                  * (1.0 - ((R_rho(k) - 1.0) &
                             / (R_rho_o - 1.0)) ** 2) ** p_2
          else  ! R_rho >= R_rho_o
             ! Large, et al, (1994) eq'n 31b
             nu%S%dd = 0.
          endif
          ! Large, et al, (1994) eq'n 31c
          nu%T%dd(k) = 0.7 * nu%S%dd(k)
       else if (0. < R_rho(k) &
            .and. R_rho(k) < 1.0 &
            .and. alpha%i(k) * T%grad_i(k) < 0. &
            .and. beta%i(k) * S%grad_i(k) < 0.) then   
          ! Diffusive convection
          ! Large, et al, (1994) eq'n 32
          nu%T%dd(k) = nu_m * 0.909 &
               * exp(4.6 * exp(-0.54 * (1.0 / R_rho(k) - 1.0)))
          ! Large, et al, (1994) eq'n 34
          if (R_rho(k) >= 0.5 .and. R_rho(k) < 1.0) then
             nu%S%dd(k) = nu%T%dd(k) * (1.85 - 0.85 / R_rho(k)) * R_rho(k)
          else if (R_rho(k) < 0.5) then
             nu%S%dd(k) = nu%T%dd(k) * 0.15 * R_rho(k)
          endif
       else
          ! No double diffusion (i.e. temperature and salinity are
          ! both stabilizing, or water column is convectively
          ! unstable)
          nu%T%dd(k) = 0.0
          nu%S%dd(k) = 0.0
       endif
    enddo
  end subroutine double_diffusion


  subroutine nu_h()
    ! Calculate the value of the total interior diffusivities, and
    ! their vertical gradients at the mixing layer depth

  end subroutine nu_h


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
    zeta = d / L_mo
    if (0. <= zeta) then
       ! Stable (eqn B1a)
       phi_m = 1. + 5. * zeta
    endif
    ! Unstable
    if (zeta_m <= zeta .and. zeta < 0.) then
       ! Eqn B1b
       phi_m = (1. - 16. * zeta) ** (-1./4.)
    elseif (zeta < zeta_m) then
       ! Eqn B1c
       phi_m = (a_m - c_m * zeta) ** (-1./3.)
    endif
  end function nondim_momentum_flux


  function nondim_scalar_flux(d) result(phi_s)
    ! Return the value of the non-dimensional scalar flux profile at
    ! the specified depth (Large, etal (1994), App. B)

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
    zeta = d / L_mo
    if (0. <= zeta) then
       ! Stable (eqn B1a)
       phi_s = 1. + 5. * zeta
    endif
    ! Unstable
    if (zeta_s <= zeta .and. zeta < 0.) then
       ! Eqn B1b
       phi_s = (1. - 16. * zeta) ** (-1./2.)
    elseif (zeta < zeta_s) then
       ! Eqn B1c
       phi_s = (a_s - c_s * zeta) ** (-1./3.)
    endif
  end function nondim_scalar_flux


  subroutine wind_driven_turbulence(h)
    ! Calculate the turbulent velocity scale profiles under conditions
    ! of wind driven turbulence (u_star /= 0) (Large, etal (1994), eqn
    ! 13)

    ! Elements from other modules:
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variable Declarations:
    use grid_mod, only: &
         grid  ! Grid parameters, and arrays

    implicit none

    ! Argument:
    real(kind=dp), intent(in) :: &
         h  ! Mixing layer extent [m]

    ! Local variables:
    real(kind=dp) :: &
         d_surf, &  ! Surface layer extent [m]
         zeta       ! Stability parameter; ratio of depth to
                    ! Monin-Obukhov length scale
    integer :: &
         j  ! Loop index over grid depth

      ! Calculate extent of surface layer, and interpolate the values
      ! of the non-dimensional flux profiles at that depth
      d_surf = epsiln * h
      do j = 0, grid%M
         ! Calculate the stability parameter value at the jth grid
         ! layer interface depth
         zeta = grid%d_i(j) / L_mo
         if (d_surf < grid%d_i(j) .and. grid%d_i(j) < h &
              .and. zeta < 0.) then
            ! Special case of depths between the surface layer depth
            ! and the mixing layer depth, in a stable mixing layer
            w%m(j) = kapa * u_star / nondim_momentum_flux(d_surf)
            w%s(j) = kapa * u_star / nondim_scalar_flux(d_surf)
         else
            ! General values of the turbulent velocity scales
            w%m(j) = kapa * u_star / nondim_momentum_flux(grid%d_i(j))
            w%s(j) = kapa * u_star / nondim_scalar_flux(grid%d_i(j))
         endif
      enddo
  end subroutine wind_driven_turbulence


  subroutine convection_driven_turbulence(h)
    ! Calculate the turbulent velocity scale profiles under conditions
    ! of convection driven turbulence (u_star = 0 and Bf < 0) (Large,
    ! etal (1994), eqn 15)

    ! Elements from other modules:
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variable Declarations:
    use grid_mod, only: &
         grid  ! Grid parameters, and arrays

    implicit none

    ! Argument:
    real(kind=dp), intent(in) :: &
         h  ! Mixing layer extent [m]
    ! Local variables:
    real(kind=dp) :: &
         d_surf, &  ! Surface layer extent [m]
         sigma      ! Non-dimensional vertical coordinate in the
                    ! mixing layer
    integer :: &
         j  ! Loop index over grid depth

    ! Calculate extent of surface layer
    d_surf = epsiln * h

    do j = 0, grid%M
       if (grid%d_i(j) < d_surf) then
          ! In the surface layer
          sigma = grid%d_i(j) / h
       elseif (d_surf <= grid%d_i(j) .and. grid%d_i(j) < h) then
          ! Between the surface layer and mixing layer depths
          sigma = epsiln
       endif
       ! Calculate the turbulent velocity scales profiles
       w%m(j) = kapa * (c_m * kapa * sigma) ** (1./3.) * w_star
       w%s(j) = kapa * (c_s * kapa * sigma) ** (1./3.) * w_star
    enddo
  end subroutine convection_driven_turbulence


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
    msg = "Turbulent velocity scale profile arrays"
    allocate(w%m(0:M), w%s(0:M), &
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
    msg = "Turbulent velocity scale profile arrays"
    deallocate(w%m, w%s, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_turbulence_variables

end module turbulence
