! $Id$
! $Source$

module physics_model
  ! Type declarations, variables, parameters, and subroutines related
  ! to the physics model in the code.
  !
  ! Public Variables:
  !
  !   B -- Water column buoyancy [m/s^2]
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
       ! Variables:
       B,      &  ! Buoyancy profile array
       ! Subroutines:
       init_physics, new_to_old_physics, dalloc_physics_variables

  ! Public variable declarations:
  real(kind=dp), dimension(:), allocatable :: &
       B  ! Buoyancy profile array

contains

  subroutine init_physics(M)
    ! Initialize physics model.
    
    ! Elements from other modules:
    ! Parameter values:
    use fundamental_constants, only: latitude, pi
    ! Subroutines:
    use water_properties, only: alloc_water_props
    use physics_eqn_builder, only: alloc_phys_RHS_variables
    use baroclinic_pressure, only: init_baroclinic_pressure
    use turbulence, only: init_turbulence
    use mixing_layer, only: init_mixing_layer
    
    implicit none
    
    ! Argument:
    integer :: M  ! Number of grid points

    ! Allocate memory for physics model variables
    call alloc_physics_variables(M)
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
    ! Allocate memory for mixing layer depth calculation variables,
    ! and initialize them.
    call init_mixing_layer(M)
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


  subroutine alloc_physics_variables(M)
    ! Allocate memory for physics model variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Buoyancy profile array"
    allocate(B(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_physics_variables


  subroutine dalloc_physics_variables
    ! Deallocate memory from physics model variables arrays.

    ! Subroutines from other modules:
    use malloc, only: dalloc_check
    use water_properties, only: dalloc_water_props
    use physics_eqn_builder, only: dalloc_phys_RHS_variables
    use baroclinic_pressure, only: dalloc_baro_press_variables
    use turbulence, only: dalloc_turbulence_variables
    
    implicit none
    
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Buoyancy profile array"
    deallocate(B, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
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
  end subroutine dalloc_physics_variables

end module physics_model
