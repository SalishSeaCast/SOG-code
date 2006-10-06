! $Id$
! $Source$

module mixing_layer
  ! Type declaration, variables, and subroutines related to the mixing
  ! layer depth calculation.

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Parameter:
       Ri_c, &  ! Critical Richardson number
       ! Subroutine:
       init_mixing_layer, find_mixing_layer_depth

  ! Public parameter declaration:
  !
  ! Critical value of Richardson number for mixed layer depth
  ! determination
  ! *** Susan was surprised that this value was not 0.25
  real(kind=dp), parameter :: Ri_c = 0.3


contains

  subroutine init_mixing_layer(M)
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
  end subroutine init_mixing_layer


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


  subroutine find_mixing_layer_depth(grid, Ri_b, year, day, day_time, count, &
       h_new)
    ! Find the mixing layer depth.  See Large, etal (1994) pp 371-372.
    use precision_defs, only: dp
    use io_unit_defs, only: stdout
    use grid_mod, only: grid_, interp_value
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: grid 
    real(kind=dp), dimension(0:), intent(in) :: Ri_b  ! Bulk Richardson no.
    ! Time and iteration count values for flagging mixing too deep events
    integer, intent(in) :: year, day, count
    real(kind=dp), intent(in) :: day_time
    real(kind=dp), intent(out) :: h_new  ! Mixing layer depth [m]

    ! Local variable:
    integer :: j_below  ! Index of grid point immediately below mixing
                        ! layer depth

    ! Find the depth at which the bulk Richardson number exceeds the
    ! critical value
    call interp_value(Ri_c, Ri_b, grid%d_g(0:grid%M), h_new, j_below)
    ! Handle mixing layer extending nearly to the bottom of the grid
    if (h_new > grid%d_g(grid%M - 3)) then
       h_new = grid%d_g(grid%M - 3)
       write(stdout, *) "ML_height_sog: Mixing too deep. ", &
            "Set h%new = ", h_new
       write(stdout, *) "Iteration count = ", count, " Time: yr = ", &
            year, " day = ", day, " day_time = ", day_time
    endif

  end subroutine find_mixing_layer_depth

end module mixing_layer
