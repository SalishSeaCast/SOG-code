! $Id$
! $Source$

module IMEX_solver
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the IMEX semi-implicit PDE solver in the
  ! SOG code.
  !
  ! Public Subroutines:
  !
  !   init_IMEX_solver -- Initialize IMEX semi-implicit PDE solver variables.
  !
  !   dalloc_IMEX_variables -- Deallocate memory for arrays for tridiagonal
  !                            matrices and right-hand sides of
  !                            diffusion/advection PDEs.
  !
  !   solve_bio_eqns -- Solve the semi-implicit diffusion/advection
  !                     PDEs for the biology quantities.


  use precision_defs, only: dp
  use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
  implicit none

  private
  public :: &
       ! Parameter values:
       ! *** Temporary for refactoring
       a_IMEX1, &
       ! Variables:
       ! *** Temporary for refactoring
       Hvector, &
       ! Subroutines:
       init_IMEX_solver, solve_bio_eqns, dalloc_IMEX_variables

  ! Type Definitions:
  !
  ! Private to module:
  !
  ! LHS matrix types:
  type :: matrices
     type(tridiag) :: &
          vel, &  ! Velocity LHS matrix
          T,   &  ! Temperature LHS matrix
          S,   &  ! Salinity LHS matrix
          bio     ! Biology quantity LHS matrix
  end type matrices
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
          Z,      &  ! Micro zooplankton RHS array
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
  ! Private to module:
  !
  type(matrices) :: &
       nAmatrix  ! LHS matrices (A) for semi-implicit PDE matrix eqns
  type(RHS) :: &
       Hvector  ! RHS vectors (H) for semi-implicit PDE matrix eqns
  real(kind=dp), dimension(:), allocatable :: &
       null_vector  ! Vector of zeros; used for quantities that don't sink

contains

  subroutine init_IMEX_solver(M)
    ! Initialize IMEX semi-implicit PDE solver variables.
    implicit none
    ! Argument:
    integer :: M  ! Number of grid points

    ! Allocate memory for IMEX solver variables
    call alloc_IMEX_variables(M)
    ! Initialize the null vector of zeros used for non-sinking quantities
    null_vector = 0.
  end subroutine init_IMEX_solver


  subroutine solve_bio_eqns(M, Pmicro, Pnano, Z, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi, day, time)
    ! Solve the semi-implicit diffusion/advection PDEs for the
    ! biology quantities.
    use precision_defs, only: dp
    use biology_eqn_builder, only: nBmatrix
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Z,      &  ! Micro Zooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    integer, intent(in) :: &
         day   ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]

    
    ! Build the RHS vectors (h) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call build_bio_Hvectors(M, Pmicro, Pnano, Z, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi)

    ! Build the LHS matrix (A) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call build_Amatrix(nBmatrix%bio%new, nAmatrix%bio)

    ! Solve the discretized semi-implicit PDE matrix equations Aq = h
    call solve_bio_tridiags(M, Pmicro, Pnano, Z, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi)

    ! Check for negative values in results, and stop with a fatal
    ! error message if any are found
    call check_bio_negatives(Pmicro, Pnano, Z, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi, day, time, fatal=.true.)
  end subroutine solve_bio_eqns


  subroutine build_bio_Hvectors(M, Pmicro, Pnano, Z, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi)
    ! Build the RHS vectors (h) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    use precision_defs, only: dp
    use biology_eqn_builder, only: nBmatrix, Pmicro_RHS, Pnano_RHS, Z_RHS, &
       NO_RHS, NH_RHS, Si_RHS, D_DON_RHS, D_PON_RHS, D_refr_RHS, D_bSi_RHS
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Z,      &  ! Microzooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile

    ! Micro phytoplankton (diatoms)
    call bio_Hvector(M, Pmicro, Pmicro_RHS%diff_adv%new, &
         Pmicro_RHS%diff_adv%old, Pmicro_RHS%bio, Pmicro_RHS%sink, &
         nBmatrix%bio%old, &
         Hvector%Pmicro)
    ! Nano phytoplankton (flagellates); non-sinking
    call bio_Hvector(M, Pnano, Pnano_RHS%diff_adv%new, &
         Pnano_RHS%diff_adv%old, Pnano_RHS%bio, null_vector, &
         nBmatrix%bio%old, &
         Hvector%Pnano)
    ! Microzooplankton: non-sinking
    call bio_Hvector(M, Z, Z_RHS%diff_adv%new, &
         Z_RHS%diff_adv%old, Z_RHS%bio, null_vector, &
         nBmatrix%bio%old, &
         Hvector%Z)
    ! Nitrate; non-sinking
    call bio_Hvector(M, NO, NO_RHS%diff_adv%new, &
         NO_RHS%diff_adv%old, NO_RHS%bio, null_vector, &
         nBmatrix%bio%old, &
         Hvector%NO)
    ! Ammonium; non-sinking
    call bio_Hvector(M, NH, NH_RHS%diff_adv%new, &
         NH_RHS%diff_adv%old, NH_RHS%bio, null_vector, &
         nBmatrix%bio%old, &
         Hvector%NH)
    ! Silicon; non-sinking
    call bio_Hvector(M, Si, Si_RHS%diff_adv%new, &
         Si_RHS%diff_adv%old, Si_RHS%bio, null_vector, &
         nBmatrix%bio%old, &
         Hvector%Si)
    ! Dissolved organic nitrogen detritus; non-sinking
    call bio_Hvector(M, D_DON, D_DON_RHS%diff_adv%new, &
         D_DON_RHS%diff_adv%old, D_DON_RHS%bio, null_vector, &
         nBmatrix%bio%old, &
         Hvector%D_DON)
    ! Particulate organic nitrogen detritus
    call bio_Hvector(M, D_PON, D_PON_RHS%diff_adv%new, &
         D_PON_RHS%diff_adv%old, D_PON_RHS%bio, D_PON_RHS%sink, &
         nBmatrix%bio%old, &
         Hvector%D_PON)
    ! Refractory nitrogen detritus
    call bio_Hvector(M, D_refr, D_refr_RHS%diff_adv%new, &
         D_refr_RHS%diff_adv%old, D_refr_RHS%bio, D_refr_RHS%sink, &
         nBmatrix%bio%old, &
         Hvector%D_refr)
    ! Biogenic silicon detritus
    call bio_Hvector(M, D_bSi, D_bSi_RHS%diff_adv%new, &
         D_bSi_RHS%diff_adv%old, D_bSi_RHS%bio, D_bSi_RHS%sink, &
         nBmatrix%bio%old, &
         Hvector%D_bSi)
  end subroutine build_bio_Hvectors


    subroutine solve_bio_tridiags(M, Pmicro, Pnano, Z, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi)
    ! Solve the discretized semi-implicit PDE matrix equations Aq = h
    use precision_defs, only: dp
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Z,      &  ! Micro zooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    
    ! Micro phytoplankton (diatoms)
    call solve_tridiag(nAmatrix%bio, Hvector%Pmicro, Pmicro(1:M))
    ! Nano phytoplankton (flagellates)
    call solve_tridiag(nAmatrix%bio, Hvector%Pnano, Pnano(1:M))
    ! Microzooplantkon
    call solve_tridiag(nAmatrix%bio, Hvector%Z, Z(1:M))
    ! Nitrate
    call solve_tridiag(nAmatrix%bio, Hvector%NO, NO(1:M))
    ! Ammonium
    call solve_tridiag(nAmatrix%bio, Hvector%NH, NH(1:M))
    ! Silicon
    call solve_tridiag(nAmatrix%bio, Hvector%Si, Si(1:M))
    ! Dissolved organic nitrogen detritus
    call solve_tridiag(nAmatrix%bio, Hvector%D_DON, D_DON(1:M))
    ! Particulate organic nitrogen detritus
    call solve_tridiag(nAmatrix%bio, Hvector%D_PON, D_PON(1:M))
    ! Refractory nitrogen detritus
    call solve_tridiag(nAmatrix%bio, Hvector%D_refr, D_refr(1:M))
    ! Biogenic silicon detritus
    call solve_tridiag(nAmatrix%bio, Hvector%D_bSi, D_bSi(1:M))
  end subroutine solve_bio_tridiags


  subroutine check_bio_negatives(Pmicro, Pnano, Z, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi, day, time, fatal)
    ! Check for negative values in results.
    use precision_defs, only: dp
    use numerics, only: check_negative
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(inout) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Z,      &  ! Microzooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    integer, intent(in) :: &
         day   ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]
    logical, intent(in), optional :: &
         fatal

    ! Micro phytoplankton (diatoms)
    call check_negative(0, Pmicro, "P%micro after solve_bio_eqns()", &
         day, time, fatal)
    ! Nano phytoplankton (flagellates)
    call check_negative(0, Pnano, "P%nano after solve_bio_eqns()", &
         day, time, fatal)
    ! Micro zooplankton 
    call check_negative (0, Z, "Z after solve_bio_eqns()", &
         day, time, fatal)
    ! Nitrate
    call check_negative(0, NO, "N%O after solve_bio_eqns()", &
         day, time, fatal)
    ! Ammonium
    call check_negative(0, NH, "N%H after solve_bio_eqns()", &
         day, time, fatal)
    ! Silicon
    call check_negative(0, Si, "Si after solve_bio_eqns()", &
         day, time, fatal)
    ! Dissolved organic nitrogen detritus
    call check_negative(0, D_DON, "D%DON after solve_bio_eqns()", &
         day, time, fatal)
    ! Particulate organic nitrogen detritus
    call check_negative(0, D_PON, "D%PON after solve_bio_eqns()", &
         day, time, fatal)
    ! Refractory nitrogen detritus
    call check_negative(0, D_refr, "D%refr after solve_bio_eqns()", &
         day, time, fatal)
    ! Biogenic silicon detritus
    call check_negative(0, D_bSi, "D%bSi after solve_bio_eqns()", &
         day, time, fatal)
  end subroutine check_bio_negatives
  

  subroutine phys_Hvector(M, qty_old, flux, flux_old, C_pg, &
       C_pg_old, Bmatrix, Hvector)
    ! Calculate RHS (Hvector) of semi-implicit PDE for a physic quantity. 
    ! It includes the value of the quantity from the previous time step, and
    ! current time step non-local flux terms, bottom and surface 
    ! flux terms, and Coriolis and pressure gradient flux terms.  
    ! Depending on the IMEX scheme being used to solve the PDE values of
    ! the quanties above, and of the diffusion coefficient matrix from
    ! the previous time step may also be blended in.
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    implicit none
    ! Arguments:
    integer :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:) :: &
         qty_old  ! Profile of the quantity at the previous time step
    real(kind=dp), dimension(1:), intent(in):: &
         flux,     &  ! Non-local flux term at current time step
         flux_old, &  ! Non-local flux term at previous time step
         C_pg,     &  ! Coriolis and pressure gradient flux term (cur)
         C_pg_old     ! Coriolis and pressure gradient flux term (prev)
    type(tridiag), intent(in) :: &
         Bmatrix  ! Diffussion coefficients matrix at previous time step
    real(kind=dp), dimension(1:), intent(out) :: &
         Hvector  ! RHS of discretized semi-implicit PDE matrix equation

    Hvector = qty_old(1:M)                                      &
         + (1.0 - a_IMEX1) * (flux_old + C_pg_old               &
                            + Bmatrix%sub * qty_old(0:M-1)  &
                            + Bmatrix%diag * qty_old(1:M)   &
                            + Bmatrix%sup * qty_old(2:M+1)) &
         + a_IMEX1 * (flux + C_pg)
  end subroutine phys_Hvector

  
  subroutine bio_Hvector(M, qty_old, diff_adv, diff_adv_old, bio, sink, &
       Bmatrix, Hvector)
    ! Calculate RHS (Hvector) of semi-implicit PDE for a biology quantity. 
    ! It includes the value of the quantity from the previous time step, and
    ! current time step upwelling flux terms, bottom and surface 
    ! flux terms, growth - mortality terms, and sinking terms.  
    ! Depending on the IMEX scheme being used to solve the PDE values of
    ! the quanties above, and of the diffusion coefficient matrix from
    ! the previous time step may also be blended in.
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    implicit none
    ! Arguments:
    integer :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:) :: &
         qty_old  ! Profile of the quantity at the previous time step
    real(kind=dp), dimension(1:), intent(in):: &
         diff_adv,     &  ! Diffusion/advection term (current time step)
         diff_adv_old, &  ! Diffusion/advection term (previous time step)
         bio,          &  ! Growth - mortality term
         sink             ! Sinking term
    type(tridiag), intent(in) :: &
         Bmatrix  ! Diffussion coefficients matrix at previous time step
    real(kind=dp), dimension(1:), intent(out) :: &
         Hvector  ! RHS of discretized semi-implicit PDE matrix equation

    Hvector = qty_old(1:M) + bio + sink                           &
         + (1.0 - a_IMEX1) * (diff_adv_old                        &
                              + Bmatrix%sub * qty_old(0:M-1)  &
                              + Bmatrix%diag * qty_old(1:M)   &
                              + Bmatrix%sup * qty_old(2:M+1)) &
         + a_IMEX1 * diff_adv
    where (abs(Hvector) < epsilon(Hvector)) Hvector = 0.
  end subroutine bio_Hvector


  subroutine build_Amatrix(Bmatrix, Amatrix)
    ! Return the LHS matrix (A) for semi-implicit PDE matrix equation
    ! from the matrix of diffusion coefficients (B) for the quantity
    ! being solved.
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    implicit none
    ! Arguments:
    type(tridiag), intent(in) :: &
         Bmatrix  ! Diffusion coefficient matrix sub-diagonal (prev timestep)
    type(tridiag), intent(inout) :: &  ! inout 'cause not all components set
         Amatrix  ! LHS matrix for discretized semi-implicit PDE matrix eqn

    Amatrix%sub = -a_IMEX1 * Bmatrix%sub
    Amatrix%diag = 1.0 - a_IMEX1 * Bmatrix%diag
    Amatrix%sup = -a_IMEX1 * Bmatrix%sup
  end subroutine build_Amatrix


  subroutine solve_tridiag(A, h, q)
    ! Solve the matrix equation Aq = h, where A is a tridiagonal matrix,
    ! h is the RHS side vector, and q is the vector of the quantity to
    ! be solved for.
    use precision_defs, only: dp
    use io_unit_defs, only: stderr
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    implicit none
    ! Arguments:
    type(tridiag), intent(in) :: &
         A  ! Tridiagonal matrix in matrix eq'n Aq = h
    real(kind=dp), dimension(1:), intent(in) :: &
         h     ! Right-hand-side vector of matrix eq'm Aq = h
    real(kind=dp), dimension(1:), intent(out) :: &
         q  ! Vector of quantity to be solved for in matrix eq'n Aq = h
    ! Local variables:
    real(kind=dp), dimension(size(h)) :: gamma
    real(kind=dp) :: beta
    integer :: j

    ! Confirm that the matrix equation is properly expressed
    if (A%diag(1) == 0.) then
       write(stderr, *) "tridiag: Malformed matrix, Amatrix(1) = 0"
       stop
    endif
    ! Decomposition and forward substitution
    beta = A%diag(1)
    q(1) = h(1) / beta
    do j = 2, size(h)
       gamma(j) = A%sup(j-1) / beta
       beta = A%diag(j) - A%sub(J) * gamma(j)
       if (beta == 0.) then
          write(stderr, *) "tridiag: Solver failed, beta = 0 at j = ", j
          stop
       endif
       q(j) = (h(j) - A%sub(j) * q(j-1)) / beta
    enddo
    ! Back substitution
    do j = size(h) - 1, 1, -1
       q(j) = q(j) - gamma(j+1) * q(j+1)
    enddo
  end subroutine solve_tridiag


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
    allocate(nAmatrix%vel%sub(1:M), nAmatrix%vel%diag(1:M), &
         nAmatrix%vel%sup(1:M),                             &
         nAmatrix%T%sub(1:M), nAmatrix%T%diag(1:M),         &
         nAmatrix%T%sup(1:M),                               &
         nAmatrix%S%sub(1:M), nAmatrix%S%diag(1:M),         &
         nAmatrix%S%sup(1:M),                               &
         nAmatrix%bio%sub(1:M), nAmatrix%bio%diag(1:M),     &
         nAmatrix%bio%sup(1:M),                             &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    allocate(Hvector%U(1:M), Hvector%V(1:M), Hvector%T(1:M),          &
         Hvector%S(1:M), Hvector%Pmicro(1:M), Hvector%Pnano(1:M),     &
         Hvector%Z(1:M),                                              &
         Hvector%NO(1:M), Hvector%NH(1:M), Hvector%Si(1:M),           &
         Hvector%D_DON(1:M), Hvector%D_PON(1:M), Hvector%D_refr(1:M), &
         Hvector%D_bSi(1:M),                                          &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "IMEX solver null vector"
    allocate(null_vector(1:M), &
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
    deallocate(nAmatrix%vel%sub, nAmatrix%vel%diag, &
         nAmatrix%vel%sup,                          &
         nAmatrix%T%sub, nAmatrix%T%diag,           &
         nAmatrix%T%sup,                            &
         nAmatrix%S%sub, nAmatrix%S%diag,           &
         nAmatrix%S%sup,                            &
         nAmatrix%bio%sub, nAmatrix%bio%diag,       &
         nAmatrix%bio%sup,                          &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    deallocate(Hvector%U, Hvector%V, Hvector%T,        &
         Hvector%S, Hvector%Pmicro, Hvector%Pnano,     &
         Hvector%Z,                                    &
         Hvector%NO, Hvector%NH, Hvector%Si,           &
         Hvector%D_DON, Hvector%D_PON, Hvector%D_refr, &
         Hvector%D_bSi,                                &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "IMEX solver null vector"
    deallocate(null_vector, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_IMEX_variables

end module IMEX_solver
