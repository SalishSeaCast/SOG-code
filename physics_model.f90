module physics_model
  ! Wrapper subroutines for setting up and tearing down the physics
  ! model in the SOG code.
  !
  ! Public Subroutines:
  !
  !   init_physics -- Initialize physics model.
  !
  !   new_to_old_physics -- Copy %new components to %old for use in
  !                         next time step.
  !
  !   dalloc_physics_variables -- Deallocate memory used by physics 
  !                               variables.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Subroutines:
       init_physics, new_to_old_physics, dalloc_physics_variables
       ! Variables


contains

  subroutine init_physics(M)
    ! Initialize physics model.
    
    ! Subroutines from other modules:
    use water_properties, only: alloc_water_props
    use physics_eqn_builder, only: alloc_phys_RHS_variables
    use baroclinic_pressure, only: init_baroclinic_pressure
    use turbulence, only: init_turbulence
    use buoyancy, only: alloc_buoyancy_variables
    use mixing_layer, only: init_mixing_layer
    use freshwater, only: init_freshwater
    use irradiance, only: init_irradiance
    use find_upwell, only: init_upwelling
    implicit none
    
    ! Argument:
    integer :: M  ! Number of grid points

    ! Allocate memory for water property arrays
    call alloc_water_props(M)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection/Coriolis/baroclinic pressure gradient
    ! equations for the physics model.
    call alloc_phys_RHS_variables(M)
    ! Allocate memory for baroclinic pressure gradient calculation
    ! variables, and initialize them.
    call init_baroclinic_pressure(M)
    ! Allocate memory for turbulence model variables, and read
    ! parameter values from infile.
    call init_turbulence(M)
    ! Allocate memory for buoyancy variables.
    call alloc_buoyancy_variables(M)
    ! Allocate memory for mixing layer depth calculation variables,
    ! and initialize them.
    call init_mixing_layer(M)
    ! Allocate memory for fresh water variables, and read
    ! parameter values from infile.
    call init_freshwater(M)
    ! Allocate memory for Kpar variables, and read parameter values
    ! from infile
    call init_irradiance()
    ! Allocate memory for upwelling quantity arrays, and read
    ! parameter values from the infile
    call init_upwelling(M)

  end subroutine init_physics


  subroutine new_to_old_physics()
    ! Copy %new component of physics variables to %old component for
    ! use at next time step.

    ! Variables from other modules:
    use core_variables, only: &
         U,  &  ! Cross-strait (35 deg) velocity component arrays
         V,  &  ! Along-strait (305 deg) velocity component arrays
         T,  &  ! Temperature profile arrays
         S      ! Salinity profile arrays

    implicit none

    ! Physics core variables
    U%old = U%new
    V%old = V%new
    T%old = T%new
    S%old = S%new
  end subroutine new_to_old_physics


  subroutine dalloc_physics_variables
    ! Deallocate memory from physics model variables arrays.

    ! Subroutines from other modules:
    use malloc, only: dalloc_check
    use water_properties, only: dalloc_water_props
    use physics_eqn_builder, only: dalloc_phys_RHS_variables
    use baroclinic_pressure, only: dalloc_baro_press_variables
    use turbulence, only: dalloc_turbulence_variables
    use buoyancy, only: dalloc_buoyancy_variables
    use mixing_layer, only: dalloc_mixing_layer_variables
    use freshwater, only: dalloc_freshwater_variables
    
    implicit none
    
    ! Deallocate memory for water property arrays
    call dalloc_water_props
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection/Coriolis/baroclinic pressure gradient
    ! equations for the physics model.
    call dalloc_phys_RHS_variables()
    ! Deallocate memory for baroclinic pressure gradient calculation
    ! arrays.
    call dalloc_baro_press_variables()
    ! Deallocate memory for turbulence variables.
    call dalloc_turbulence_variables()
    ! Deallocate memory for buoyancy variables.
    call dalloc_buoyancy_variables()
    ! Deallocate memory for mixing layer variables.
    call dalloc_mixing_layer_variables()
    ! Deallocate memory for fresh water variables.
    call dalloc_freshwater_variables()
  end subroutine dalloc_physics_variables

end module physics_model
