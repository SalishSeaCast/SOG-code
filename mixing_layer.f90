! $Id$
! $Source$

module mixing_layer
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the mixing layer depth calculation in the
  ! SOG code.
  !
  ! Public Parameters:
  !
  !   Ri_c -- Critical Richardson number
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
  !   new_to_old_mixing_layer -- Copy %new component of mixing layer
  !                              depth variables to %old component for
  !                              use at next time step.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       Ri_c, &    ! Critical Richardson number
       ! Variables:
       h, &  ! Mixing layer depth values & indices
       ! Subroutines:
       init_mixing_layer, find_mixing_layer_depth, &
       find_mixing_layer_indices, new_to_old_mixing_layer

  ! Type Definitions:
  !
  ! Mixing layer
  type mixing_layer_depth
     real(kind=dp) :: &
          new, &  ! Depth at current time step
          old     ! Depth at previous time step
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

contains

  subroutine init_mixing_layer(M)
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Initialize the mixing layer depth, and the indices of the grid
    ! point & grid layer interface immediately below it
    h%new = 2.0d0
    call find_mixing_layer_indices()
  end subroutine init_mixing_layer


  subroutine new_to_old_mixing_layer()
    ! Copy %new component of mixing layer depth variables to %old
    ! component for use at next time step.

    implicit none

    h%old = h%new
  end subroutine new_to_old_mixing_layer


  subroutine find_mixing_layer_depth(grid, Ri_b, Bf, year, day, day_time, &
       count)
    ! Find the mixing layer depth.  See Large, etal (1994) pp 371-372.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use grid_mod, only: grid_
    ! Parameter Value Declarations:
    use io_unit_defs, only: stdout
    use fundamental_constants, only: f
    ! Variables
    use turbulence, only: u_star, L_mo
    ! Subroutines:
    use grid_mod, only: interp_value

    implicit none

    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid parameters and arrays
    real(kind=dp), dimension(0:), intent(in) :: &
         Ri_b  ! Bulk Richardson no.
    real(kind=dp), intent(in) :: &
         Bf  ! Surface buoyancy flux
    integer, intent(in) :: &
         year, &  ! Year for flagging mixing too deep events
         day,  &  ! Year-day  for flagging mixing too deep events
         count    ! Iteration count for flagging mixing too deep events
    real(kind=dp), intent(in) :: &
         day_time  ! Day-sec for flagging mixing too deep events

    ! Local variable:
    real(kind=dp) :: &
         d_Ekman  ! Ekman depth [m]

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


  subroutine alloc_mixing_layer(M)
    ! Allocate memory for core variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

  end subroutine alloc_mixing_layer


  subroutine dalloc_mixing_layer
  end subroutine dalloc_mixing_layer

end module mixing_layer
