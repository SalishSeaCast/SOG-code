! $Id$
! $Source$

module mixing_layer
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the mixing layer depth calculation in the
  ! SOG code.
  !
  ! Public Variables:
  !
  !   h -- Mixing layer depth values & indices
  !
  ! Public Subroutines:
  !
  !   init_mixing_layer --
  !
  !   find_mixing_layer_depth --
  !
  !   find_mixing_layer_indices -- Set the value of the indices of the
  !                                grid point & grid layer interface
  !                                immediately below the mixing layer
  !                                depth.
  !
  !   dalloc_mixing_layer_variables -- Deallocate memory from mixing
  !                                    layer depth variables.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       h, &  ! Mixing layer depth values & indices
       ! Subroutines:
       init_mixing_layer, find_mixing_layer_depth, &
       find_mixing_layer_indices, dalloc_mixing_layer_variables

  ! Type Definitions:
  !
  ! Mixing layer
  type mixing_layer_depth
     real(kind=dp) :: &
          new  ! Depth at current time step
     integer :: &
          i, &  ! Index of grid layer interface immediately below
                ! mixing layer depth
          g     ! Index of grid point immediately below mixing layer
                ! depth
  end type mixing_layer_depth

  ! Public Parameter Declarations:
  !
  real(kind=dp), parameter :: &
       Ri_c = 0.3d0  ! Critical value of Richardson number for mixed
                     ! layer depth determination

  ! Variable Declarations:
  !
  ! Public
  type(mixing_layer_depth) :: &
       h  ! Mixing layer depth values & indices
  
  ! Private:
  real(kind=dp), dimension(:), allocatable :: &
       N, &      ! Local buoyancy frequency profile
       V_t_2, &  ! Profile of the square fo the turbulent velocity shear
       Ri_b      ! Bulk Richardson number profile

contains

  subroutine init_mixing_layer(M)
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Allocate memory for mixing layer depth related arrays
    call alloc_mixing_layer_variables(M)
    ! Initialize the mixing layer depth, and the indices of the grid
    ! point & grid layer interface immediately below it
    h%new = 2.0d0
    call find_mixing_layer_indices()
  end subroutine init_mixing_layer


  subroutine find_mixing_layer_depth(year, day, day_time, count)
    ! Find the mixing layer depth.  See Large, etal (1994) pp 371-372.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Parameter Value Declarations:
    use io_unit_defs, only: stdout
    use fundamental_constants, only: f
    ! Variables
    use grid_mod, only: &
         grid  ! Grid parameters and arrays
    use turbulence, only: &
         u_star, &  ! Turbulent friction velocity
         L_mo       ! Monin_Obukhov length scale
    use declarations, only: Bf  ! *** Should come from somewhere else
    ! Subroutines:
    use grid_mod, only: interp_value

    implicit none

    ! Arguments:
    integer, intent(in) :: &
         year, &  ! Year for flagging mixing too deep events
         day,  &  ! Year-day  for flagging mixing too deep events
         count    ! Iteration count for flagging mixing too deep events
    real(kind=dp), intent(in) :: &
         day_time  ! Day-sec for flagging mixing too deep events

    ! Local variable:
    real(kind=dp) :: &
         d_Ekman  ! Ekman depth [m]

    ! Calculate the profile of the bulk Richardson number.
    call calc_Ri_b()
    ! Find the depth at which the bulk Richardson number exceeds the
    ! critical value
    call interp_value(Ri_c, 0, Ri_b, grid%d_g(0:grid%M), h%new, h%g)

    ! Apply the Ekman and Monin-Obukhov depth criteria to the mixing
    ! layer depth when stable forcing exists.  Note that abs(x) <
    ! epsilon(x) is a real-number-robust test for x == 0, and abs(x) >
    ! epsilon(x) is similarly for x /= 0.
    if (Bf > 0. &
         .or. (abs(Bf) < epsilon(Bf) &
               .and. abs(u_star) > epsilon(u_star))) then
       ! Calculate the Ekmann depth (Large, et al (1994), eqn 24)
       d_Ekman = 0.7d0 * u_star / f
       ! Under stable forcing, mixing layer depth is the minimum of
       ! the values from Richardson number, Ekman depth, and Monim
       ! -Obukhov length scale criteria
       h%new = min(h%new, d_Ekman, L_mo)
       ! But it also can't be shallower than the depth of the 1st grid
       ! point
       h%new = max(h%new, grid%d_g(1))
    endif

    ! Handle mixing layer extending nearly to the bottom of the grid
    if (h%new > grid%d_g(grid%M - 3)) then
       h%new = grid%d_g(grid%M - 3)
       write(stdout, *) "find_mixing_layer_depth: Mixing too deep. ", &
            "Set h%new = ", h%new
       write(stdout, *) "Iteration count = ", count, " Time: yr = ", &
            year, " day = ", day, " day_time = ", day_time
    endif

    ! Set the value of the indices of the grid point & grid layer
    ! interface immediately below the mixing layer depth.
    call find_mixing_layer_indices()
  end subroutine find_mixing_layer_depth


  subroutine calc_Ri_b()
    ! Calculate the profile of the bulk Richardson number.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Functions:
    use grid_mod, only: depth_average
    ! Parameter Values:
    use fundamental_constants, only: g
    use turbulence, only: epsiln
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and arrays
    use water_properties, only: &
         rho  ! Density profile
    use core_variables, only: &
         U, &  ! Cross-strait (35 deg) velocity component arrays
         V     ! Along-strait (305 deg) velocity component arrays
         

    ! Local Variables:
    real(kind=dp) :: &
         d_surf,   & ! Surface layer extent [m]
         U_r, V_r, & ! Average value of velocity components in the
                     ! surface layer.
         rho_r       ! Average value of the density in the surface
                     ! layer.
    integer :: &
         M  ! Number of grid points (to improve readability of
            ! vectorized code).


    ! Note that we are going to make use of the relation B = (alpha *
    ! T - beta * S) = -g * rho / rho_o here to express the buoyancy
    ! terms in Large, et al (1994), eqn 21 in terms of rho.

    ! Calculate the average values of the density, and the velocity
    ! components in the surface layer.
    d_surf = epsiln * h%new
    rho_r = depth_average(rho%g, 0.0d0, d_surf)
    U_r = depth_average(U%new, 0.0d0, d_surf)
    V_r = depth_average(V%new, 0.0d0, d_surf)

    ! Calculate the profile of turbulent velocity shear
    call calc_turbulent_vel_shear()

    ! Calculate the profile of the bulk Richardson number
    !
    ! Note that we are going to make use of the relation B = (alpha *
    ! T - beta * S) = -g * rho / rho_o here to express the buoyancy
    ! terms in Large, et al (1994), eqn 21 in terms of rho.
    M = grid%M
    Ri_b(0) = 0.0d0
    Ri_b(1:M) = (-g * (rho_r - rho%g(1:M)) / rho%g(0)) * grid%d_g(1:M)  &
         / ((U_r - U%new(1:M)) ** 2 + (V_r - V%new(1:M)) ** 2 &
            + V_t_2(1:M) + epsilon(U_r))
  end subroutine calc_Ri_b


  subroutine calc_turbulent_vel_shear()
    ! Calculate the profile of turbulent velocity shear for use in the
    ! bulk Richardson number profile calculation.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Functions:
    use turbulence, only: w_scale, nondim_scalar_flux
    ! Parameter Values:
    use fundamental_constants, only: &
         g  ! Acceleration due to gravity
    use turbulence, only: &
       c_s,    &  ! Coefficient of non-dimensional turbulent scalar
                  ! flux profile in 1/3 power law regime
       kapa,   &  ! von Karman constant
       epsiln     ! Non-dimensional extent of the surface layer
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and arrays
    use water_properties, only: &
         rho  ! Density profile
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays
    use declarations, only: Bf  ! *** Should come from somewhere else

    ! Local Parameters:
    real(kind=dp), parameter :: &
         Cv = 1.5d0  ! Ratio of interior local buoyancy frequency to
                     ! same at entrainment depth.  *** Large, et al
                     ! (1994), table 1 recommends a value of 1.5 to
                     ! 1.6

    ! Local Variables:
    real(kind=dp) :: &
         beta_t, &  ! Ratio of entrainment flux to surface buoyancy flux
         w_s        ! Value of scalar turbulent velocity scale
    integer :: &
         j  ! Index over grid depth
    
    ! Calculate the profile of the local buoyancy frequency at the
    ! grid point depths.
    N = sqrt(-(g / rho%g(1:grid%M)) * rho%grad_g(1:grid%M))
    ! Set the value of the ratio of entrainment flux to surface
    ! buoyancy flux, and calculate the turbulent velocity shear
    ! profile (Large, et al (1994), eqn 23).
    if (wbar%b(0) > 0.0d0) then
       ! Constant when we have convection (see Large, et al (1994),
       ! table 1).
       beta_t = -0.2d0
       do j = 1, grid%M
          w_s = w_scale(h%new, grid%d_g(j), Bf, nondim_scalar_flux, c_s)
          ! Large, et al (1994), eqn 23.
          V_t_2(j) = Cv * sqrt(-beta_t / (c_s * epsiln)) &
               * grid%d_g(j) * N(j) * w_s / (Ri_c * kapa ** 2)
       enddo
    else
       V_t_2 = 0.0d0
    endif
  end subroutine calc_turbulent_vel_shear


  subroutine find_mixing_layer_indices()
    ! Set the value of the indices of the grid point & grid layer
    ! interface immediately below the mixing layer depth.

    ! Elements from other modules:
    !
    ! Variable:
    use grid_mod, only: grid
    ! Subroutine:
    use grid_mod, only: interp_value

    implicit none

    ! Using interp_value() here is a convenient way of getting h%g
    ! without duplicating code that searches through the grid
    call interp_value(h%new, 0, grid%d_g, grid%d_g, h%new, h%g)
    if (grid%d_i(h%g - 1) > h%new) then
       ! Mixing layer depth is in the grid layer above the grid point
       h%i = h%g - 1
    else
       ! Mixing layer depth is in the same grid layer as the grid
       ! point
       h%i = h%g
    endif
  end subroutine find_mixing_layer_indices


  subroutine alloc_mixing_layer_variables(M)
    ! Allocate memory for core variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Local buoyancy frequency profile array"
    allocate(N(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Square of turbulent velocity shear profile array"
    allocate(V_t_2(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Bulk Richardson number profile array"
    allocate(Ri_b(0:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_mixing_layer_variables


  subroutine dalloc_mixing_layer_variables()
    ! Deallocate memory for turbulence quantities.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Local buoyancy frequency profile array"
    deallocate(N, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Square of turbulent velocity shear profile array"
    deallocate(V_t_2, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Bulk Richardson number profile array"
    deallocate(Ri_b, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_mixing_layer_variables

end module mixing_layer
