! $Id$
! $Source$

module grid_mod
  ! Type definitions, functions, and subroutines related to the grid.
  ! *** For the moment, just grid interpolation factors and function,
  ! *** but eventually all of the grid definition code will reside here.
  !
  ! Public Function:
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

  use precision_defs

  implicit none

  private
  public :: &
       ! Functions:
       interp_i, &
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
    use malloc
    ! *** This is temporary, but we need the grid type-def
    use mean_param
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
    use malloc
    implicit none
    ! Local variables:
    integer              :: dallocstat  ! Allocation return status
    character(len=80)    :: msg         ! Allocation failure message prefix

    msg = "Grid interface interpolation factor arrays"
    deallocate(above_g%factor, below_g%factor, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_grid


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
