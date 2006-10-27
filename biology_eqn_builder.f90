! $Id$
! $Source$

module biology_eqn_builder
  ! Type definitions, variable declarations, and subroutines related
  ! to building the semi-implicit diffusion/advection equations for
  ! the biology quantities in the SOG code.
  !
  ! Public Variables:
  !
  !   diff_coeffs_bio -- Tridiagonal matrix of diffusion coefficient values
  !
  !   Pmicro_RHS -- Micro phytoplankton (diatoms) right-hand side arrays
  !
  !   Pnano_RHS -- Nano phytoplankton (flagellates) right-hand side arrays
  !
  !   NO_RHS -- Nitrate concentration right-hand side arrays
  !
  !   NH_RHS -- Ammonium concentration right-hand side arrays
  !
  !   Si_RHS -- Silicon concentration right-hand side arrays
  !
  !   D_DON_RHS -- Dissolved organic nitrogen detritus right-hand side arrays
  !
  !   D_PON_RHS -- Particulate organic nitro detritus right-hand side arrays
  !
  !   D_refr_RHS -- Refractory nitrogen detritus right-hand side arrays
  !
  !   D_bSi_RHS -- Biogenic silicon detritus right-hand side arrays
  !
  ! Public Subroutines:
  !
  !   build_biology_equations -- Build the right-hand side (RHS) arrays
  !                              for the diffusion/advection equations for
  !                              the biology quantities.
  !
  !   new_to_old_bio_RHS -- Copy %new component of the biology *_RHS%diff_adv
  !                         arrays to the %old component for use by the
  !                         IMEX semi-impllicit PDE solver.
  !
  !   alloc_bio_RHS_variables -- Allocate memory for biology RHS arrays
  !
  !   dalloc_bio_RHS_variables -- Deallocate memory for biology RHS arrays

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       diff_coeffs_bio, &  ! Tridiagonal matrix of diffusion coefficient values
       Pmicro_RHS,      &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,       &  ! Nano phytoplankton (flagellates) RHS arrays
       NO_RHS,          &  ! Nitrate concentration RHS arrays
       NH_RHS,          &  ! Ammonium concentration RHS arrays
       Si_RHS,          &  ! Silicon concentration RHS arrays
       D_DON_RHS,       &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,       &  ! Particulate organic nitrogen detritus RHS arrays
       D_refr_RHS,      &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS,       &  ! Biogenic silicon detritus RHS arrays
       ! Subroutines:
       build_biology_equations, new_to_old_bio_RHS, &
       alloc_bio_RHS_variables, dalloc_bio_RHS_variables

  ! Type Definitions:
  !
  ! Private to module:
  !
  ! Tridiagnonal matrix:
  type :: tridiag
     real(kind=dp), dimension(:), allocatable :: &
          sub_diag, &  ! Sub-diagonal vector of a tridiagonal matrix
          diag,     &  ! Diagonal vector of a tridiagonal matrix
          super_diag   ! Super-diagonal vector of a tridiagonal matrix
  end type tridiag
  !
  ! New/old array components:
  type :: new_old
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old
  !
  ! Semi-implicit diffusion/advection equation right-hand side (RHS) arrays
  type :: RHS
     type(new_old) :: &
          diff_adv  ! Diffusion/advection component of RHS
     real(kind=dp), dimension(:), allocatable :: &
          bio, &  ! Biology (growth - mortality) component of RHS
          sink    ! Sinking component of RHS
  end type RHS

  ! Variable Declarations:
  !
  ! Public:
  type(tridiag) :: &
       diff_coeffs_bio  ! Tridiagonal matrix of diffusion coefficient values
  type(RHS) :: &
       Pmicro_RHS, &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,  &  ! Nano phytoplankton (flagellates) RHS arrays
       NO_RHS,     &  ! Nitrate concentration RHS arrays
       NH_RHS,     &  ! Ammonium concentration RHS arrays
       Si_RHS,     &  ! Silicon concentration RHS arrays
       D_DON_RHS,  &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,  &  ! Particulate organic nitro detritus RHS arrays
       D_refr_RHS, &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS      ! Biogenic silicon detritus RHS arrays

contains

  subroutine build_biology_equations(grid, dt, Pmicro, Pnano, NO, NH,  & ! in
       Si, D_DON, D_PON, D_refr, D_bSi, Ft, K_all, wupwell)
    ! Build the right-hand side (RHS) arrays for the
    ! diffusion/advection equations for the biology quantities.
    !
    ! This calculates the values of the diffusion coefficients matrix
    ! (diff_coeff_bio%*), the RHS diffusion/advection term vectors
    ! (*_RHS%diff_adv%new), and the RHS sinking term vectors (*_RHS%sink).
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use diffusion, only: new_diffusion_coeff, diffusion_bot_surf_flux
    use find_upwell, only: vertical_advection
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
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
    real(kind=dp), intent(in) :: &
         Ft  ! Total fresh water flux
    real(kind=dp), dimension(0:), intent(in) :: &
         K_all  ! Total diffusion coefficient array
    real(kind=dp), dimension(1:), intent(in) :: &
         wupwell  ! Profile of vertical upwelling velocity [m/s]
    ! Local variables:
    real(kind=dp) :: &
         aflux  ! Surface nutrient flux

    ! Calculate the strength of the diffusion coefficients for biology
    ! model quantities.  They diffuse like salinity.
    call new_diffusion_coeff(grid, dt, K_all,            & ! in
         diff_coeffs_bio%sub_diag, diff_coeffs_bio%diag, & ! out
         diff_coeffs_bio%super_diag)                       ! out

    ! Initialize the RHS *%diff_adv%new arrays, and calculate the diffusive
    ! fluxes at the bottom and top of the grid
    aflux = -Ft * (0. - Pmicro(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         Pmicro(grid%M+1),                               &   ! in
         Pmicro_RHS%diff_adv%new)                            ! out
    aflux = -Ft * (0. - Pnano(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         Pnano(grid%M+1),                                &   ! in
         Pnano_RHS%diff_adv%new)                             ! out
    aflux = -Ft * (6.4 - NO(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         NO(grid%M+1),                                   &   ! in
         NO_RHS%diff_adv%new)                                ! out
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         NH(grid%M+1),                                   &   ! in
         NH_RHS%diff_adv%new)                                ! out
    aflux = -Ft * (60.0 - Si(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         Si(grid%M+1),                                   &   ! in
         Si_RHS%diff_adv%new)                                ! out
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         D_DON(grid%M+1),                                &   ! in
         D_DON_RHS%diff_adv%new)                             ! out
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         D_PON(grid%M+1),                                &   ! in
         D_PON_RHS%diff_adv%new)                             ! out
    D_refr_RHS%diff_adv%new = 0.
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         D_bSi(grid%M+1),                                &   ! in
         D_bSi_RHS%diff_adv%new)                             ! out

     ! Add vertical advection due to upwelling
     call vertical_advection(grid, dt, Pmicro, wupwell, &
          Pmicro_RHS%diff_adv%new)
     call vertical_advection(grid, dt, Pnano, wupwell, &
          Pnano_RHS%diff_adv%new)
     call vertical_advection(grid, dt, NO, wupwell, &
          NO_RHS%diff_adv%new)
     call vertical_advection(grid, dt, NH, wupwell, &
          NH_RHS%diff_adv%new)
     call vertical_advection(grid, dt, Si, wupwell, &
          Si_RHS%diff_adv%new)
     call vertical_advection(grid, dt, D_DON, wupwell, &
          D_DON_RHS%diff_adv%new)
     call vertical_advection(grid, dt, D_PON, wupwell, &
          D_PON_RHS%diff_adv%new)
     call vertical_advection(grid, dt, D_bSi, wupwell, &
          D_bSi_RHS%diff_adv%new)


  end subroutine build_biology_equations


  subroutine new_to_old_bio_RHS()
    ! Copy %new component of the biology *_RHS%diff_adv arrays to the
    ! %old component for use by the IMEX semi-impllicit PDE solver.
    implicit none

    Pmicro_RHS%diff_adv%old = Pmicro_RHS%diff_adv%new
    Pnano_RHS%diff_adv%old = Pnano_RHS%diff_adv%new
    NO_rhs%diff_adv%old = NO_RHS%diff_adv%new
    NH_rhs%diff_adv%old = NH_RHS%diff_adv%new
    Si_rhs%diff_adv%old = Si_RHS%diff_adv%new
    D_DON_rhs%diff_adv%old = D_DON_RHS%diff_adv%new
    D_PON_rhs%diff_adv%old = D_PON_RHS%diff_adv%new
    D_refr_rhs%diff_adv%old = D_refr_RHS%diff_adv%new
    D_bSi_rhs%diff_adv%old = D_bSi_RHS%diff_adv%new
  end subroutine new_to_old_bio_RHS


  subroutine alloc_bio_RHS_variables(M)
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
    allocate(diff_coeffs_bio%sub_diag(1:M), diff_coeffs_bio%diag(1:M), &
         diff_coeffs_bio%super_diag(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Micro phytoplankton RHS arrays"
    allocate(Pmicro_RHS%diff_adv%new(1:M), Pmicro_RHS%diff_adv%old(1:M), &
         Pmicro_RHS%bio(1:M), Pmicro_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    allocate(Pnano_RHS%diff_adv%new(1:M), Pnano_RHS%diff_adv%old(1:M), &
         Pnano_RHS%bio(1:M), Pnano_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    allocate(NO_RHS%diff_adv%new(1:M), NO_RHS%diff_adv%old(1:M), &
         NO_RHS%bio(1:M), NO_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    allocate(NH_RHS%diff_adv%new(1:M), NH_RHS%diff_adv%old(1:M), &
         NH_RHS%bio(1:M), NH_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Silicon concentration RHS arrays"
    allocate(Si_RHS%diff_adv%new(1:M), Si_RHS%diff_adv%old(1:M), &
         Si_RHS%bio(1:M), Si_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    allocate(D_DON_RHS%diff_adv%new(1:M), D_DON_RHS%diff_adv%old(1:M), &
         D_DON_RHS%bio(1:M), D_DON_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    allocate(D_PON_RHS%diff_adv%new(1:M), D_PON_RHS%diff_adv%old(1:M), &
         D_PON_RHS%bio(1:M), D_PON_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    allocate(D_refr_RHS%diff_adv%new(1:M), D_refr_RHS%diff_adv%old(1:M), &
         D_refr_RHS%bio(1:M), D_refr_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    allocate(D_bSi_RHS%diff_adv%new(1:M), D_bSi_RHS%diff_adv%old(1:M), &
         D_bSi_RHS%bio(1:M), D_bSi_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_bio_RHS_variables


  subroutine dalloc_bio_RHS_variables
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
    deallocate(diff_coeffs_bio%sub_diag, diff_coeffs_bio%diag, &
         diff_coeffs_bio%super_diag, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Micro phytoplankton RHS arrays"
    deallocate(Pmicro_RHS%diff_adv%new, Pmicro_RHS%diff_adv%old, &
         Pmicro_RHS%bio, Pmicro_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    deallocate(Pnano_RHS%diff_adv%new, Pnano_RHS%diff_adv%old, &
         Pnano_RHS%bio, Pnano_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    deallocate(NO_RHS%diff_adv%new, NO_RHS%diff_adv%old, &
         NO_RHS%bio, NO_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    deallocate(NH_RHS%diff_adv%new, NH_RHS%diff_adv%old, &
         NH_RHS%bio, NH_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Silicon concentration RHS arrays"
    deallocate(Si_RHS%diff_adv%new, Si_RHS%diff_adv%old, &
         Si_RHS%bio, Si_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    deallocate(D_DON_RHS%diff_adv%new, D_DON_RHS%diff_adv%old, &
         D_DON_RHS%bio, D_DON_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    deallocate(D_PON_RHS%diff_adv%new, D_PON_RHS%diff_adv%old, &
         D_PON_RHS%bio, D_PON_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    deallocate(D_refr_RHS%diff_adv%new, D_refr_RHS%diff_adv%old, &
         D_refr_RHS%bio, D_refr_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    deallocate(D_bSi_RHS%diff_adv%new, D_bSi_RHS%diff_adv%old, &
         D_bSi_RHS%bio, D_bSi_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_bio_RHS_variables

end module biology_eqn_builder
