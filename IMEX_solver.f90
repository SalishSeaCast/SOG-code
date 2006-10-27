! $Id$
! $Source$

module IMEX_solver
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the IMEX semi-implicit PDE solver in the
  ! SOG code.
  !
  ! Public Variables:
  !
  !   eggs -- Description
  !
  ! Public Subroutines:
  !
  !   ham -- Description

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       ! *** Temporary for refactoring
       a_IMEX1
!!$       ! Variables:
!!$       eggs, &
!!$       ! Subroutines:
!!$       ham

  ! Type Definitions:
  !
  ! Public:
  !
  ! Private to module:

  ! Parameter Value Declarations:
  !
  ! Public:
  !
  ! Private to module:
  !
  ! Implicit-explicit first order scheme:
  !   a_IMEX1 = 1.0 --> backward Euler in diffusion
  !           = 0.5 --> Crank-Nicolson in diffusion
  !           = 0.0 --> forward Euler in diffusion
  ! Both have forward Euler in explicit terms (i.e. Coriolis and
  ! boundary conditions)
  real(kind=dp), parameter :: a_IMEX1 = 1.0

contains

  subroutine alloc_IMEX_variables(M)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
!!$    allocate(diff_coeffs_bio%sub_diag(1:M), diff_coeffs_bio%diag(1:M), &
!!$         diff_coeffs_bio%super_diag(1:M), &
!!$         stat=allocstat)
!!$    call alloc_check(allocstat, msg)
  end subroutine alloc_IMEX_variables

end module IMEX_solver
