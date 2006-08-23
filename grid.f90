! $Id$
! $Source$

module grid_mod
  ! Type definitions, functions, and subroutines related to the grid.
  ! *** For the moment, just grid interpolation factors and function,
  ! *** but eventually all of the grid definition code will reside here.
  !
  ! Public Function:
  !
  ! interp_d(qty, d) -- Return the interpolated value of a quantity
  !                     stored at the grid points at the specified depth.
  !
  ! interp_i(qty_g) -- Return the interpolated values of a quantity at
  !                    the grid interface depths from its values at the 
  !                    grid point depths.
  ! 
  ! Public Subroutines:
  !
  ! init_grid(M, grid) -- Allocate memory for grid arrays and
  !                       initialize grid interpolation factor arrays.
  !
  ! dalloc_grid() -- Deallocate memory from grid arrays.

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Functions:
       interp_d, interp_i, &
       ! Subroutines:
       init_grid, dalloc_grid

  ! Private type definition:
  !
  ! Interpolation factors
  type :: interp_factor
     real(kind=dp), dimension(:), pointer :: factor
  end type interp_factor

  ! Private module variable declarations:
  !
  ! Fraction of grid point spacing that interface j is above grid point j+1
  type(interp_factor) :: above_g
  ! Fraction of grid point spacing that interface j is below grid point j
  type(interp_factor) :: below_g

contains

  subroutine init_grid(M, grid)
    ! Initialize grid interpolation factor values
    ! *** Eventually this will absorb define_grid() and do all grid
    ! *** initialization.
    use malloc, only: alloc_check
    ! *** This is temporary, but we need the grid type-def
    use mean_param, only: gr_d
    implicit none
    ! Argument:
    integer, intent(in)    :: M     ! Number of grid points
    type(gr_d), intent(in) :: grid  ! Grid depths & spacings
    ! Local variables:
    integer              :: allocstat  ! Allocation return status
    character(len=80)    :: msg        ! Allocation failure message prefix

    ! Allocate memory for grid interface interpolation factor arrays
    msg = "Grid interface interpolation factor arrays"
    allocate(above_g%factor(1:M+1), below_g%factor(0:M), stat=allocstat)
    call alloc_check(allocstat, msg)
    ! Calculate the grid interface interpolation factors
    above_g%factor(1:M) = abs(grid%d_i(0:M-1) - grid%d_g(1:M)) &
         / grid%g_space(0:M-1)
    above_g%factor(M+1) = 0.
    below_g%factor(0) = 0.
    below_g%factor(1:M) = abs(grid%d_g(1:M) - grid%d_i(1:M)) &
         / grid%g_space(1:M)
  end subroutine init_grid


  subroutine dalloc_grid
    ! Deallocate memory for grid interface interpolation factors
    ! *** and eventually the rest of the grid
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer              :: dallocstat  ! Allocation return status
    character(len=80)    :: msg         ! Allocation failure message prefix

    msg = "Grid interface interpolation factor arrays"
    deallocate(above_g%factor, below_g%factor, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_grid


  function interp_d(qty_g, d, grid) result(d_value)
    ! Return the interpolated value of a quantity stored at the grid
    ! points at the specified depth.
    ! *** use of mean_param can be removed when grid comes into this
    ! *** module
    use mean_param, only: gr_d
    use io_unit_defs, only: stderr
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    real(kind=dp), intent(in) :: d
    type(gr_d), intent(in) :: grid
    ! Result:
    real(kind=dp) :: d_value
    ! Local variables:
    integer :: j, j_above

    ! Make sure the requested depth is within the grid
    if (d < grid%d_g(0) .or. d > grid%d_g(grid%M+1)) then
       write(stderr, *) "Warning: d = ", d, " out of range in interp_d"
       d_value = 9999999999.
    endif
    ! Find the index of the grid point above the specified depth
    do j = 0, grid%M
       if (d > grid%d_g(j)) then
          j_above = j
       else
          exit
       endif
    enddo
    ! Interpolate the quantity value at the specified depth from its
    ! value at the grid points above and below
    d_value = qty_g(j_above) + (qty_g(j_above+1) - qty_g(j_above)) &
         * (d - grid%d_g(j_above)) / grid%g_space(j_above)
  end function interp_d


  function interp_i(qty_g) result(qty_i)
    ! Return the interpolated values of a quantity at the grid interface
    ! depths from its values at the grid point depths.
    implicit none
    ! Arguments:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    ! Result:
    ! Quantity values at grid interface depths
    real(kind=dp), dimension(0:size(qty_g) - 2) :: qty_i
    ! *** size(qty_g) - 2) can be replaced by grid%M when grid comes
    ! *** into this module
    qty_i = qty_g(0:size(qty_g) - 2) * above_g%factor &
         + qty_g(1:) * below_g%factor
  end function interp_i

end module grid_mod
