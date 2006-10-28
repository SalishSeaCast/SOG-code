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
  !   solve_biology_PDEs -- 

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       ! *** Temporary for refactoring
       a_IMEX1, &
       Hvector, &
       nAmatrix, &
       ! Subroutines:
       init_IMEX_solver, solve_biology_PDEs, dalloc_IMEX_variables

  ! Type Definitions:
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
  type :: matrices
     type(tridiag) :: &
          vel, &  ! Velocity LHS matrix
          T,   &  ! Temperature LHS matrix
          S,   &  ! Salinity LHS matrix
          bio, &  ! Biology quantity LHS matrix
          null    ! Identity LHS matrix for non-diffusing quantities
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
    ! Initialize the identity Amatrix component vecctors used for
    ! non-diffusing quantities (e.g. refractory nitrogen detritus)
    nAmatrix%null%sub_diag = 0.
    nAmatrix%null%diag = 1.
    nAmatrix%null%super_diag = 0.
    ! Initialize the null vector of zeros used for non-sinking quantities
    null_vector = 0.
  end subroutine init_IMEX_solver


  subroutine solve_biology_PDEs(M, Pmicro, Pnano, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi)
    !
    use precision_defs, only: dp
    use biology_eqn_builder, only: nBmatrix, Pmicro_RHS, Pnano_RHS, &
       NO_RHS, NH_RHS, Si_RHS, D_DON_RHS, D_PON_RHS, D_refr_RHS, D_bSi_RHS
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    
    ! Build the RHS vectors (H) for the biology semi-implicit PDEs
    !
    ! Micro phytoplankton (diatoms)
    call bio_Hvector(M, Pmicro, Pmicro_RHS%diff_adv%new, &
         Pmicro_RHS%diff_adv%old, Pmicro_RHS%bio, Pmicro_RHS%sink, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%Pmicro)
    ! Nano phytoplankton (flagellates); non-sinking
    call bio_Hvector(M, Pnano, Pnano_RHS%diff_adv%new, &
         Pnano_RHS%diff_adv%old, Pnano_RHS%bio, null_vector, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%Pnano)
    ! Nitrate; non-sinking
    call bio_Hvector(M, NO, NO_RHS%diff_adv%new, &
         NO_RHS%diff_adv%old, NO_RHS%bio, null_vector, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%NO)
    ! Ammonium; non-sinking
    call bio_Hvector(M, NH, NH_RHS%diff_adv%new, &
         NH_RHS%diff_adv%old, NH_RHS%bio, null_vector, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%NH)
    ! Silicon; non-sinking
    call bio_Hvector(M, Si, Si_RHS%diff_adv%new, &
         Si_RHS%diff_adv%old, Si_RHS%bio, null_vector, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%Si)
    ! Dissolved organic nitrogen detritus; non-sinking
    call bio_Hvector(M, D_DON, D_DON_RHS%diff_adv%new, &
         D_DON_RHS%diff_adv%old, D_DON_RHS%bio, null_vector, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%D_DON)
    ! Particulate organic nitrogen detritus
    call bio_Hvector(M, D_PON, D_PON_RHS%diff_adv%new, &
         D_PON_RHS%diff_adv%old, D_PON_RHS%bio, D_PON_RHS%sink, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%D_PON)
    ! Refractory nitrogen detritus
    call bio_Hvector(M, D_refr, D_refr_RHS%diff_adv%new, &
         D_refr_RHS%diff_adv%old, D_refr_RHS%bio, D_refr_RHS%sink, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%D_refr)
    ! Biogenic silicon detritus
    call bio_Hvector(M, D_bSi, D_bSi_RHS%diff_adv%new, &
         D_bSi_RHS%diff_adv%old, D_bSi_RHS%bio, D_bSi_RHS%sink, &
         nBmatrix%bio%old%sub, nBmatrix%bio%old%diag, nBmatrix%bio%old%sup, &
         Hvector%D_bSi)
  end subroutine solve_biology_PDEs


  subroutine phys_Hvector(M, qty_old, flux, flux_old, C_pg, &
       C_pg_old, Bmatrix_old_sub, Bmatrix_old_diag, Bmatrix_old_sup, Hvector)
    ! Calculate RHS (Hvector) of semi-implicit PDE for a physic quantity. 
    ! It includes the value of the quantity from the previous time step, and
    ! current time step non-local flux terms, bottom and surface 
    ! flux terms, and Coriolis and pressure gradient flux terms.  
    ! Depending on the IMEX scheme being used to solve the PDE values of
    ! the quanties above, and of the diffusion coefficient matrix from
    ! the previous time step may also be blended in.
    use precision_defs, only: dp
    implicit none
    ! Arguments:
    integer :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:) :: &
         qty_old  ! Profile of the quantity at the previous time step
    real(kind=dp), dimension(1:), intent(in):: &
         flux,             &  ! Non-local flux term at current time step
         flux_old,         &  ! Non-local flux term at previous time step
         C_pg,             &  ! Coriolis and pressure gradient flux term (cur)
         C_pg_old,         &  ! Coriolis and pressure gradient flux term (prev)
         Bmatrix_old_sub,  &  ! Diff coeff matrix sub-diagonal (prev timestep)
         Bmatrix_old_diag, &  ! Diff coeff matrix sub-diagonal (prev timestep)
         Bmatrix_old_sup      ! Diff coeff matrix sub-diagonal (prev timestep)
    real(kind=dp), dimension(1:), intent(out) :: &
         Hvector  ! RHS of semi-implicit PDE matrix equation

    Hvector = qty_old(1:M)                                      &
         + (1.0 - a_IMEX1) * (flux_old + C_pg_old               &
                            + Bmatrix_old_sub * qty_old(0:M-1)  &
                            + Bmatrix_old_diag * qty_old(1:M)   &
                            + Bmatrix_old_sup * qty_old(2:M+1)) &
         + a_IMEX1 * (flux + C_pg)
  end subroutine phys_Hvector

  
  subroutine bio_Hvector(M, qty_old, diff_adv, diff_adv_old, bio, sink, &
       Bmatrix_old_sub, Bmatrix_old_diag, Bmatrix_old_sup, Hvector)
    ! Calculate RHS (Hvector) of semi-implicit PDE for a biology quantity. 
    ! It includes the value of the quantity from the previous time step, and
    ! current time step upwelling flux terms, bottom and surface 
    ! flux terms, growth - mortality terms, and sinking terms.  
    ! Depending on the IMEX scheme being used to solve the PDE values of
    ! the quanties above, and of the diffusion coefficient matrix from
    ! the previous time step may also be blended in.
    use precision_defs, only: dp
    implicit none
    ! Arguments:
    integer :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:) :: &
         qty_old  ! Profile of the quantity at the previous time step
    real(kind=dp), dimension(1:), intent(in):: &
         diff_adv,         &  ! Diffusion/advection term (current time step)
         diff_adv_old,     &  ! Diffusion/advection term (previous time step)
         bio,              &  ! Growth - mortality term
         sink,             &  ! Sinking term
         Bmatrix_old_sub,  &  ! Diff coeff matrix sub-diagonal (prev timestep)
         Bmatrix_old_diag, &  ! Diff coeff matrix sub-diagonal (prev timestep)
         Bmatrix_old_sup      ! Diff coeff matrix sub-diagonal (prev timestep)
    real(kind=dp), dimension(1:), intent(out) :: &
         Hvector  ! RHS of semi-implicit PDE matrix equation

    Hvector = qty_old(1:M) + bio + sink                           &
         + (1.0 - a_IMEX1) * (diff_adv_old                        &
                              + Bmatrix_old_sub * qty_old(0:M-1)  &
                              + Bmatrix_old_diag * qty_old(1:M)   &
                              + Bmatrix_old_sup * qty_old(2:M+1)) &
         + a_IMEX1 * diff_adv
    where (abs(Hvector) < epsilon(Hvector)) Hvector = 0.
  end subroutine bio_Hvector


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
    allocate(nAmatrix%vel%sub_diag(1:M), nAmatrix%vel%diag(1:M), &
         nAmatrix%vel%super_diag(1:M),                          &
         nAmatrix%T%sub_diag(1:M), nAmatrix%T%diag(1:M),         &
         nAmatrix%T%super_diag(1:M),                            &
         nAmatrix%S%sub_diag(1:M), nAmatrix%S%diag(1:M),         &
         nAmatrix%S%super_diag(1:M),                            &
         nAmatrix%bio%sub_diag(1:M), nAmatrix%bio%diag(1:M),     &
         nAmatrix%bio%super_diag(1:M),                          &
         nAmatrix%null%sub_diag(1:M), nAmatrix%null%diag(1:M),   &
         nAmatrix%null%super_diag(1:M),                         &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    allocate(Hvector%U(1:M), Hvector%V(1:M), Hvector%T(1:M),          &
         Hvector%S(1:M), Hvector%Pmicro(1:M), Hvector%Pnano(1:M),     &
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
    deallocate(nAmatrix%vel%sub_diag, nAmatrix%vel%diag, &
         nAmatrix%vel%super_diag,                       &
         nAmatrix%T%sub_diag, nAmatrix%T%diag,           &
         nAmatrix%T%super_diag,                         &
         nAmatrix%S%sub_diag, nAmatrix%S%diag,           &
         nAmatrix%S%super_diag,                         &
         nAmatrix%bio%sub_diag, nAmatrix%bio%diag,       &
         nAmatrix%bio%super_diag,                       &
         nAmatrix%null%sub_diag, nAmatrix%null%diag,     &
         nAmatrix%null%super_diag,                      &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "RHS vectors (H) for semi-implicit PDE matrix eqns"
    deallocate(Hvector%U, Hvector%V, Hvector%T,        &
         Hvector%S, Hvector%Pmicro, Hvector%Pnano,     &
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
