! $Id$
! $Source$

module grid_mod
  ! Type definitions, variables, functions, and subroutines related to
  ! the grid.
  !
  ! Public Type:
  !
  !   grid_ -- Grid parameters:
  !              D -- depth of bottom of grid [m]
  !              M -- number of grid points
  !              d_g -- array of depths of grid points [m]
  !              d_i -- array of depths of grid cell interfaces [m]
  !              g_space -- array of depths of grid point spacings [m]
  !              i_space -- array of depths of grid cell interface spacings [m]
  !
  ! Public Variables:
  !
  !   grid   -- Grid point and interface depth and spacing arrays.
  !
  ! Public Functions:
  !
  !   interp_value -- Return the value of the quantity in the
  !                   interp_array that corresponds to the the specified
  !                   value in the search_array.  This can be used to
  !                   find the value of a quantity at a specified depth:
  !                      qty_g = find(d, d_g, qty_g) or
  !                      qty_i = find(d, d_i, qty_i)
  !                   or to find the shallowest depth at which a specified
  !                   quantity value occurs:
  !                      d = find(q, qty_g, d_g) or
  !                      d = find(q, qty_i, d_i)
  !
  !   interp_i(qty_g) -- Return the interpolated values of a quantity at
  !                      the grid interface depths from its values at the 
  !                      grid point depths.
  !
  !   gradient_i(qty_g) -- Return the values of the gradients at the
  !                        grid layer interface depths of a quantity stored
  !                        at the grid point depths.
  !
  !   gradient_g(qty_g) -- Return the values of the gradients at the
  !                        grid point depths of a quantity stored
  !                        at the grid point depths.
  ! 
  ! Public Subroutines:
  !
  !   init_grid(M, grid) -- Allocate memory for grid arrays and
  !                         initialize grid interpolation factor arrays.
  !
  !   dalloc_grid() -- Deallocate memory from grid arrays.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Types:
       grid_, &
       ! Variables:
       grid, &
       ! Functions:
       interp_value, interp_i, gradient_i, gradient_g, &
       ! Subroutines:
       init_grid, dalloc_grid

  ! Public type definition:
  !
  ! Grid parameters
  type :: grid_
     integer :: M               ! Number of grid points
     real(kind=dp) :: D         ! Depth of bottom of grid [m]
     real(kind=dp), dimension(:), pointer :: d_g  ! Grid point depths [m]
     real(kind=dp), dimension(:), pointer :: d_i  ! Grid interface depths [m]
     real(kind=dp), dimension(:), pointer :: g_space  ! Grid point spacing [m]
     real(kind=dp), dimension(:), pointer :: i_space  ! Grid i'face spacing [m]
  end type grid_

  ! Public variable declarations:
  !
  ! Grid point and interface depth and spacing arrays
  type(grid_) :: grid

  ! Private type definition:
  !
  ! Interpolation factors
  type :: interp_factor
     real(kind=dp), dimension(:), pointer :: factor
  end type interp_factor

  ! Private module variable declarations:
  !
  ! Grid spacing parameter (<0 concentrates resolution near
  ! surface, =0 produced uniform grid, >0 concentrates resolution 
  ! near bottom of grid).  See Large, et al (1994), App. D.
  real(kind=dp) :: lambda
  real(kind=dp), dimension(:), pointer :: xsi_i  ! Interface indices
  real(kind=dp), dimension(:), pointer :: xsi_g  ! Grid point indices
  ! Fraction of grid point spacing that interface j is above grid point j+1
  type(interp_factor) :: above_g
  ! Fraction of grid point spacing that interface j is below grid point j
  type(interp_factor) :: below_g
  ! Array of index numbers used by the inter_value function
  integer, dimension(:), pointer :: indices

contains

  subroutine init_grid
    ! Initialize grid.
    use precision_defs, only: dp, sp
    use io_unit_defs, only: stderr
    implicit none
    ! Local variables:
    integer :: j  ! Loop index over depth of grid

    ! Read the grid parameters
    call read_grid
    ! Allocate memory for grid arrays
    call alloc_grid

    ! Calculate the grid point and interface indices per Large, et al
    ! (1994) Appendix D (expressions are defined in the paragraph
    ! following the one containing eq'n D1)
    do j = 0, grid%M
       xsi_i(j) = dble(j) / dble(grid%M)
    enddo
    xsi_g(0) = 0.
    do j = 1, grid%M
       xsi_g(j) = (dble(j) - 0.5) / dble(grid%M)
    enddo
    xsi_g(grid%M + 1) = xsi_i(grid%M)

    ! Calculate the grid point and interface depths (Large, etal 1994 eq'n D1)
    if (abs(lambda) < epsilon(lambda)) then
       ! Uniform grid (lambda == 0)
       grid%d_g = grid%D * xsi_g
       grid%d_i = grid%D * xsi_i
    else
       ! Non-uniform grid (lambda < 0 concentrates resolution near
       ! surface, lambda > 0 concentrates resolution near D)
       ! *** Note that implementation of non-uniform grid is not
       ! *** consistent throughout the code
       write(stderr, *) "init_grid: Non-uniform grid is not ", &
            "fully implemented, lambda = ", lambda
       stop
       grid%d_g = (grid%D / lambda) * log(1. * xsi_g * (1. - exp(lambda)))
       grid%d_i = (grid%D / lambda) * log(1. * xsi_i * (1. - exp(lambda)))
    endif

    ! Calculate the grid point and interface spacings
    grid%g_space = grid%d_g(1:) - grid%d_g(0:grid%M)
    grid%i_space = grid%d_i(1:) - grid%d_i(0:grid%M-1)

    ! Calculate the grid interface interpolation factors
    above_g%factor(1:grid%M) = abs(grid%d_i(0:grid%M-1) - grid%d_g(1:grid%M)) &
         / grid%g_space(0:grid%M-1)
    above_g%factor(grid%M+1) = 0.
    below_g%factor(0) = 0.
    below_g%factor(1:grid%M) = abs(grid%d_g(1:grid%M) - grid%d_i(1:grid%M)) &
         / grid%g_space(1:grid%M)

    ! Set index numbers for use by interp_value() function
    forall (j = 1:size(indices)) indices(j) = j
  end subroutine init_grid


  subroutine read_grid
    ! Read the grid parameters from the input file
    use input_processor, only: getpari, getpard
    implicit none

    ! Read the model depths, the number of grid points, and the grid
    ! spacing factor values.
    grid%D = getpard("maxdepth")
    grid%M = getpari("gridsize")
    lambda = getpard("lambda")
  end subroutine read_grid


  subroutine alloc_grid
    ! Allocate memory for grid arrays.
    use malloc, only: alloc_check
    implicit none
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Grid point and interface index arrays"
    allocate(xsi_g(0:grid%M+1), xsi_i(0:grid%M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Grid point and interface depth and spacing arrays"
    allocate(grid%d_g(0:grid%M+1), grid%d_i(0:grid%M), &
         grid%g_space(0:grid%M), grid%i_space(1:grid%M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Grid interface interpolation factor arrays"
    allocate(above_g%factor(1:grid%M+1), below_g%factor(0:grid%M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Grid index number array"
    allocate(indices(1:grid%M+2), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_grid


  subroutine dalloc_grid
    ! Deallocate memory for grid arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Grid point and interface index arrays"
    deallocate(xsi_g, xsi_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Grid point and interface depth and spacing arrays"
    deallocate(grid%d_g, grid%d_i, grid%g_space, grid%i_space, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Grid interface interpolation factor arrays"
    deallocate(above_g%factor, below_g%factor, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Grid index number array"
    deallocate(indices, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_grid


  function interp_value(value, search_array, interp_array) result(interp_val)
    ! Return the value of the quantity in the interp_array that
    ! corresponds to the the specified value in the search_array.
    ! This can be used to find the value of a quantity at a specified
    ! depth:
    !        qty_g = find(d, d_g, qty_g) or
    !        qty_i = find(d, d_i, qty_i)
    ! or to find the shallowest depth at which a specified quantity
    ! value occurs:
    !        d = find(q, qty_g, d_g) or
    !        d = find(q, qty_i, d_i)
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: value  ! Depth or quantity value to find
    real(kind=dp), dimension(1:), intent(in) :: &
         search_array, &  ! Array to find value in
         interp_array     ! Array to interpolate over
    ! Result:
    real(kind=dp) :: interp_val
    ! Local variables:
    logical, dimension(1:size(search_array)) :: mask
    integer :: j_below

    ! Mask the search array for values less than the specified
    ! value, and confirm that the specified value is within the range
    ! of the search array
    mask = search_array > value
    if (any(mask) .and. .not. all(mask)) then
       ! Find the index of the search array element below the
       ! specified value
       j_below = minval(indices, mask)
       ! Interpolate the value of the quantity in interp_array that
       ! corresponds to the specified value in the search array
       interp_val = interp_array(j_below-1) &
            + (interp_array(j_below) - interp_array(j_below-1)) &
            * (value - search_array(j_below-1)) &
            / (search_array(j_below) - search_array(j_below-1))
    else
       write(stdout, *) "Warning: value = ", value, &
            " out of range in interp_value()"
       interp_val = 9999999999.
    endif
  end function interp_value


  function interp_i(qty_g) result(qty_i)
    ! Return the interpolated values of a quantity at the grid interface
    ! depths from its values at the grid point depths.
    implicit none
    ! Arguments:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    ! Result:
    ! Quantity values at grid interface depths
    real(kind=dp), dimension(0:grid%M) :: qty_i

    qty_i = qty_g(0:size(qty_g) - 2) * above_g%factor &
         + qty_g(1:) * below_g%factor
  end function interp_i


  function gradient_i(qty_g) result(grad_i)
    ! Return the values of the gradients at the grid layer interface
    ! depths of a quantity stored at the grid point depths.
    ! *** This function may not work for an unevenly spaced grid
    implicit none
    ! Argument:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    ! Result:
    ! Values of the gradients at grid interface depths
    real(kind=dp), dimension(1:grid%M) :: grad_i

    grad_i = (qty_g(1:grid%M) - qty_g(2:grid%M+1)) / grid%g_space(1:grid%M)
  end function gradient_i


  function gradient_g(qty_g) result(grad_g)
    ! Return the values of the gradients at the grid point
    ! depths of a quantity stored at the grid point depths.
    ! *** This function does not work for an unevenly spaced grid
    implicit none
    ! Argument:
    ! Quantity values at grid point depths
    real(kind=dp), dimension(0:), intent(in) :: qty_g
    ! Result:
    ! Values of the gradients at grid point depths
    real(kind=dp), dimension(1:grid%M) :: grad_g

    grad_g = (qty_g(0:grid%M-1) - qty_g(2:grid%M+1)) &
         / (grid%g_space(0:grid%M-1) + grid%g_space(1:))
  end function gradient_g

end module grid_mod
