! $Id$
! $Source$

module IMEX_solver
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the IMEX semi-implicit PDE solver in the
  ! SOG code.
  !
  ! Public Subroutines:
  !
  !   alloc_IMEX_variables -- Allocate memory for arrays for tridiagonal
  !                           matrices and right-hand sides of
  !                           diffusion/advection PDEs.
  !
  !   dalloc_IMEX_variables -- Deallocate memory for arrays for tridiagonal
  !                            matrices and right-hand sides of
  !                            diffusion/advection PDEs.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       ! *** Temporary for refactoring
       a_IMEX1, &
       ! Subroutines:
       alloc_IMEX_variables, dalloc_IMEX_variables

  ! Type Definitions:
  !
  ! Public:
  !
  ! Private to module:
  !
  ! Tridiagnonal matrix vectors:
  type :: tridiag
     real(kind=dp), dimension(:), allocatable :: &
          sub_diag, &  ! Sub-diagonal vector of a tridiagonal matrix
          diag,     &  ! Diagonal vector of a tridiagonal matrix
          super_diag   ! Super-diagonal vector of a tridiagonal matrix
  end type tridiag
  !
  ! LHS matrix types:
  type :: matrix
     type(tridiag) :: &
          vel, &  ! Velocity LHS matrix
          T,   &  ! Temperature LHS matrix
          S,   &  ! Salinity LHS matrix
          bio     ! Biology quantity LHS matrix
  end type matrix
  !
  ! Semi-implicit diffusion/advection PDE right-hand side (RHS) arrays
  type :: RHS
     real(kind=dp), dimension(:), allocatable :: &
          U,      &  ! U (cross-strait) velocity RHS array
          V,      &  ! V (along-strait) velocity RHS array
          T,      &  ! Temperature RHS array
          S,      &  ! Salinity RHS array
          Pmicro, &  ! Micro phytoplankton (diatoms) RHS array
          Pnano,  &  ! Nano phytoplankton (flagellates) RHS array
          NO,     &  ! Nitrate concentration RHS array
          NH,     &  ! Ammonium concentration RHS array
          Si,     &  ! Silicon concentration RHS array
          D_DON,  &  ! Dissolved organic nitrogen detritus RHS array
          D_PON,  &  ! Particulate organic nitrogen detritus RHS array
          D_refr, &  ! Refractory nitrogen detritus RHS array
          D_bSi      ! Biogenic silicon detritus RHS array
  end type RHS

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

  ! Variable Declarations:
  !
  ! Public:
  !
  ! Private to module:
  !
  type(tridiag) :: &
       Amatrix  ! LHS matrices (A) for semi-implicit PDE matrix eqns
  type(RHS) :: &
       Hvector  ! RHS vectors (H) for semi-implicit PDE matrix eqns

contains

  subroutine alloc_IMEX_variables(M)
    ! Allocate memory for arrays for tridiagonal matrices and
    ! right-hand sides of diffusion/advection PDEs.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "LHS matrices (A) for semi-implicit PDE matrix eqns"
    allocate(Amatrix%vel(1:M), Amatrix%T(1:M), Amatrix%S(1:M), &
         Amatrix%bio(1:M),                                     &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    allocate(Hvector%U(1:M), Hvector%V(1:M), Hvector%T(1:M),          &
         Hvector%S(1:M), Hvector%Pmicro(1:M), Hvector%Pnano(1:M),     &
         Hvecctor%NO(1:M), Hvector%NH(1:M), Hvector%Si(1:M),          &
         Hvector%D_DON(1:M), Hvector%D_PON(1:M), Hvector%D_refr(1:M), &
         Hvector%D_bSi(1:M),                                          &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_IMEX_variables


  subroutine dalloc_IMEX_variables
    ! Deallocate memory for arrays for tridiagonal matrices and
    ! right-hand sides of diffusion/advection PDEs.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "LHS matrices (A) for semi-implicit PDE matrix eqns"
    deallocate(Amatrix%vel, Amatrix%T, Amatrix%S, &
         Amatrix%bio,                             &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    deallocate(Hvector%U, Hvector%V, Hvector%T,        &
         Hvector%S, Hvector%Pmicro, Hvector%Pnano,     &
         Hvecctor%NO, Hvector%NH, Hvector%Si,          &
         Hvector%D_DON, Hvector%D_PON, Hvector%D_refr, &
         Hvector%D_bSi,                                &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_IMEX_variables

end module IMEX_solver
