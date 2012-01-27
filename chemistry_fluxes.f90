module chemistry_fluxes
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the chemistry air/sea gas fluxes in the SOG code.
  !
  ! Public Type Definitions:
 
  !
  ! Public Parameters:
  !
  !   spam -- Description
  !
  ! Public subroutine:
  !
  !   thing -- Description

  implicit none

  private
  public :: &
       ! Variables:
       ! Subroutines:
       solve_gas_flux

contains

  subroutine solve_gas_flux(grid, T, S, rho, unow, vnow, DIC, Oxy, &
       day, time, pCO2, pO2)
    ! Iteration for successful diffusion of CO2 and oxygen gas fluxes

    ! Variables:
    use numerics, only: chem_dt, chem_steps

    ! Type Definitions:
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use io_unit_defs, only: stdout

    ! Subroutines from other modules:
    use air_sea_fluxes, only: gas_flux
    use carbonate, only: calculate_co2
    use biology_eqn_builder, only: &
         new_to_old_chem_RHS, new_to_old_chem_Bmatrix

    implicit none

    ! Arguments:
    type(grid_), intent(in) :: &
         grid          ! Grid arrays
    integer, intent(in) :: &
         day           ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         T,         &  ! Sea surface water temperature [K]
         S,         &  ! Sea surface practical salinity [PSU]
         rho,       &  ! Sea surface water density  [kg/m^3]
         unow,      &  ! 35 degree wind component (cross-strait)
         vnow,      &  ! 325 degree wind component (along-strait)
         time          ! Time of current time step [s since start of run]
    real(kind=dp), dimension(0:), intent(inout) :: &
         DIC,       &  ! Surface dissolved inorganic carbon [uM]
         Oxy           ! Surface dissolved oxygen [uM]
    real(kind=dp), intent(out) :: &
         pCO2,      &  ! Partial pressure CO2 [ppm]
         pO2           ! Partial pressure O2 [ppm]

    ! Local variables:
    integer :: &
         count         ! Loop index
    real(kind=dp) :: &
         DIC_flux,  &  ! Dissolved inorganic carbon surface flux [umol m-2 s-1]
         Oxy_flux,  &  ! Oxygen surface flux [umol m-2 s-1]
         CO2           ! Surface carbon dioxide [uM]

    do count = 1, chem_steps !---- Begin Iteration Loop ----

       call new_to_old_chem_RHS()
       call new_to_old_chem_Bmatrix()

       ! Calculate surface CO2 from surface DIC
       call calculate_co2(T, S, rho, DIC(1), CO2)

       ! Calculate surface CO2 gas flux
       call gas_flux('CO2', T, S, CO2, unow, vnow, DIC_flux, pCO2)
       call gas_flux('Oxy', T, S, Oxy(1), unow, vnow, Oxy_flux, pO2)

       ! This calculates the values of the precursor diffusion
       ! coefficients matrix (Bmatrix%bio%*), the RHS diffusion/advection
       ! term vectors (*_RHS%diff_adv%new), and the RHS sinking term
       ! vectors (*_RHS%sink).
       call build_chem_equations(grid, chem_dt, DIC, DIC_flux, Oxy, Oxy_flux)

       ! Solve the semi-implicit diffusion/advection PDEs for the
       ! biology quantities.
       call solve_chem_equations(grid%M, DIC, Oxy, day, time)

    enddo !---- End Iteration Loop ----

  end subroutine solve_gas_flux


  subroutine build_chem_equations(grid, dt, DIC, DIC_flux, Oxy, Oxy_flux)
    ! Build the terms for the diffusion/advection equations for the
    ! gas flux dependent quantities DIC and Oxy.
    !
    ! This calculates the values of the precursor diffusion
    ! coefficients matrix (Bmatrix%bio%*), the RHS diffusion/advection
    ! term vectors (*_RHS%diff_adv%new), and the RHS sinking term
    ! vectors (*_RHS%sink).
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use biology_eqn_builder, only: &
         Bmatrix,    &  ! Precursor diffusion coefficient matrices
         DIC_RHS,    &  ! Dissolved inorganic carbon RHS arrays
         Oxy_RHS        ! Dissolved oxygen RHS arrays
    use turbulence, only: K
    use diffusion, only: & ! diffusion_coeff, diffusion_bot_surf_flux, &
         diffusion_nonlocal_fluxes, diffusion_coeff
    use upwelling, only: upwelling_advection
    use freshwater, only: freshwater_bio
    use buoyancy, only: Bf

    implicit none

    ! Arguments:
    type(grid_), intent(in) :: &
         grid          ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt,        &  ! Time step [s]
         DIC_flux,  &  ! Dissolved inorganic carbon surface flux [umol m-2 s-1]
         Oxy_flux      ! Oxygen surface flux [umol m-2 s-1]
    real(kind=dp), dimension(0:), intent(in) :: &
         DIC,       &  ! Dissolved inorganic carbon
         Oxy           ! Dissolved oxygen

    ! Local variables:
    real(kind=dp) :: &
         surf_flux  ! surface nutrient flux when all the river water on surface
    real(kind=dp), dimension(0:grid%M):: &
         distrib_flux  ! distributed nutrient flux

    ! Calculate the strength of the diffusion coefficients for biology
    ! model quantities.  They diffuse like salinity.
    call diffusion_coeff(dt, K%S, &      ! in
         Bmatrix%chem%new)               ! out

    ! Initialize the RHS *%diff_adv%new arrays, and calculate the diffusive
    ! fluxes at the bottom and top of the grid
    call freshwater_bio ('DIC', DIC(0:grid%M),               &
         surf_flux, distrib_flux)
    call diffusion_nonlocal_fluxes(dt, K%S, 0.0d0, Bf,       &  ! in
         surf_flux + DIC_flux, distrib_flux, DIC(grid%M+1),  &  ! in
         DIC_RHS%diff_adv%new)                                  ! out
    call freshwater_bio ('Oxy', Oxy(0:grid%M),               &
         surf_flux, distrib_flux)
    call diffusion_nonlocal_fluxes(dt, K%S, 0.0d0, Bf,       &  ! in
         surf_flux + Oxy_flux, distrib_flux, Oxy(grid%M+1),  &  ! in
         Oxy_RHS%diff_adv%new)                                  ! out

    ! Add vertical advection due to upwelling
    call upwelling_advection(dt, DIC, DIC_RHS%diff_adv%new)
    call upwelling_advection(dt, Oxy, Oxy_RHS%diff_adv%new)

  end subroutine build_chem_equations


  subroutine solve_chem_equations(M, DIC, Oxy, day, time)
    ! Solve the semi-implicit diffusion/advection PDEs for the
    ! gas flux dependent quantities DIC and Oxy.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use io_unit_defs, only: stdout
    use numerics, only: tridiag
    ! Variables:
    use IMEX_solver, only: Amatrix, Hvector, null_vector
    use biology_eqn_builder, only: &
         Bmatrix,    &  ! Precursor diffusion coefficient matrices
         DIC_RHS,    &  ! Dissolved inorganic carbon RHS arrays
         Oxy_RHS        ! Dissolved oxygen RHS arrays
    ! Subroutines:
    use IMEX_solver, only: bio_Hvector, build_Amatrix, solve_tridiag 
    use numerics, only: check_negative

    implicit none
    
    ! Arguments:
    integer, intent(in) :: &
         M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(inout) :: &
         DIC,    &  ! Dissolved inorganic carbon
         Oxy        ! Dissolved oxygen
    integer, intent(in) :: &
         day   ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]

    ! Build the RHS vectors (h) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call bio_Hvector(M, DIC, DIC_RHS%diff_adv%new, &
         DIC_RHS%diff_adv%old, DIC_RHS%bio, null_vector, &
         Bmatrix%chem%old, &
         Hvector%DIC)
    call bio_Hvector(M, Oxy, Oxy_RHS%diff_adv%new, &
         Oxy_RHS%diff_adv%old, Oxy_RHS%bio, null_vector, &
         Bmatrix%chem%old, &
         Hvector%Oxy)

    ! Build the LHS matrix (A) for the discretized semi-implicit PDE
    ! matrix equations Aq = h
    call build_Amatrix(Bmatrix%chem%new, Amatrix%bio)

    ! Solve the discretized semi-implicit PDE matrix equations Aq = h
    call solve_tridiag(Amatrix%bio, Hvector%DIC, DIC(1:M))
    call solve_tridiag(Amatrix%bio, Hvector%Oxy, Oxy(1:M))

    ! Check for negative values in results, and print with a warning
    ! message if any are found
    call check_negative(0, DIC, "DIC after solve_tridiag()", &
         day, time, fatal=.false.)
    call check_negative(0, Oxy, "Oxy after solve_tridiag()", &
         day, time, fatal=.false.)

  end subroutine solve_chem_equations

end module chemistry_fluxes
