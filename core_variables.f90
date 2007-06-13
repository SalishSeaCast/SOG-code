! $Id$
! $Source$

module core_variables
  ! Type definitions, variable declarationss, and subroutines related
  ! to the core variables that the SOG code calculates.
  ! 
  ! Public Variables:
  !
  !   U -- Velocity component in the u (cross-strait, 35 deg) direction [m/s]
  !
  !   V -- Velocity component in the v (along-strait, 305 deg) direction [m/s]
  !
  !   T -- Water column temperature [K]
  !
  !   S -- Water column salinity [-]
  !
  !   P -- Phytoplankton biomasses:
  !          P%micro -- Micro phytoplankton (diatoms) [uM N]
  !          P%nano -- Nano phytoplankton (meso-rub) [uM N]
  !          P%pico -- Pico phytoplankton (flagellates) [uM N]
  !
  !   Z -- Microzooplankton biomass [uM N]
  !
  !   N -- Nitrogen compounds concentrations:
  !          N%O -- Nitrate concentration [uM N]
  !          N%H -- Ammonium concentration [uM N]
  !
  !   Si -- Silicon concentration [uM]
  !
  !   D -- Detritus concentrations:
  !          D%DON -- Dissolved organic nitrogen [uM N]
  !          D%PON -- Particulate organic nitrogen [uM N]
  !          D%refr -- Refractory nitrogen [uM N]
  !          D%bSi -- Biogenic silicon [uM Si]
  !
  ! Public Subroutines:
  !
  !   alloc_core_variables -- Allocate memory for core variables arrays.
  !
  !   dalloc_core_variables -- De-allocate memory for core variables arrays.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       U,  &  ! Cross-strait (35 deg) velocity component arrays
       V,  &  ! Along-strait (305 deg) velocity component arrays
       T,  &  ! Temperature profile arrays
       S,  &  ! Salinity profile arrays
       P,  &  ! Micro, nano & pico phytoplankton profile arrays
       Z,  &  ! Microzooplankton profile array
       N,  &  ! Nitrate & ammonium concentation profile arrays
       Si, &  ! Silicon concentration profile arrays
       D,  &  ! Detritus concentration profile arrays
       ! Types (as required by new pg compiler)
       profiles, & ! type for U, V, T, S
       nitrogen, & ! type for N
       plankton, & ! type for P
       detritus, & ! type for D
                   ! Z and Si are just dp real
       ! Subroutines:
       alloc_core_variables, dalloc_core_variables

  ! Private module type definitions:
  !
  ! Velocities, temperature, and salinity
  type :: profiles
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Profile of quantity at current time setp
          old, &  ! Profile of quantity at previous time step
          grad_i  ! Profile of gradient of quantity at grid layer interfaces
  end type profiles
  !
  ! Nitrogen compounds
  type :: nitrogen
     real(kind=dp), dimension(:), allocatable :: &
          O, &  ! N%O is nitrate (NO3) concentration profile
          H     ! H%H is ammonium (NH4) concentration profile
  end type nitrogen
  !
  ! Plankton
  type :: plankton
     real(kind=dp), dimension(:), allocatable :: &
          micro, &  ! P%micro is micro phytoplankton (diatoms) biomass profile
          nano, &   ! P%nano is nano phytoplankton (meso-rub) biomass profile
          pico      ! P%pico is pico phytoplankton (flagellate) biomass profile
  end type plankton
  !
  ! Detritus
  type :: detritus
     real(kind=dp), dimension(:), allocatable :: &
          DON,  &  ! Dissolved organic nitrogen [uM N]
          PON,  &  ! Particulate organic nitrogen [uM N]
          refr, &  ! Refractory nitrogen [uM N]
          bSi      ! Biogenic silicon [uM Si]
  end type detritus


  ! Public variable declarations:
  type(profiles) :: &
       U, &  ! Cross-strait (35 deg) velocity component arrays
       V, &  ! Along-strait (305 deg) velocity component arrays
       T, &  ! Temperature profile arrays
       S     ! Salinity profile arrays
  type(plankton) :: &
       P  ! Micro & nano & pico phytoplankton biomass profile arrays
  real(kind=dp), dimension(:), allocatable :: &
       Z  ! Microzooplankton concentration profile array
  type(nitrogen) :: &
       N  ! Nitrate & ammonium concentration profile arrays
  real(kind=dp), dimension(:), allocatable :: &
       Si ! Silicon concentration profile array
  type(detritus) :: &
       D  ! Detritus concentration profile arrays

contains

  subroutine alloc_core_variables(M)
    ! Allocate memory for core variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Cross-strait velocity component profile arrays"
    allocate(U%new(0:M+1), U%old(0:M+1), U%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Along-strait velocity component profile arrays"
    allocate(V%new(0:M+1), V%old(0:M+1), V%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Temperature profile arrays"
    allocate(T%new(0:M+1), T%old(0:M+1), T%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Salinity profile arrays"
    allocate(S%new(0:M+1), S%old(0:M+1), S%grad_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Micro phytoplankton biomass profile array"
    allocate(P%micro(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nano phytoplankton biomass profile array"
    allocate(P%nano(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Pico phytoplankton biomass profile array"
    allocate(P%pico(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Microzooplankton concentration profile array"
    allocate(Z(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrate concentration profile array"
    allocate(N%O(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ammonium concentration profile array"
    allocate(N%H(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Silicon concentration profile array"
    allocate(Si(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved organic nitrogen detritus concentration profile array"
    allocate(D%DON(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Particulate organic nitrogen detritus concentration profile array"
    allocate(D%PON(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Refractory nitrogen detritus concentration profile array"
    allocate(D%refr(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Biogenic silicon detritus concentration profile array"
    allocate(D%bSi(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_core_variables


  subroutine dalloc_core_variables
    ! Deallocate memory for core variables arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg        ! Deallocation failure message prefix

    msg = "Cross-strait velocity component profile arrays"
    deallocate(U%new, U%old, U%grad_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Along-strait velocity component profile arrays"
    deallocate(V%new, V%old, V%grad_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Temperature profile arrays"
    deallocate(T%new, T%old, T%grad_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Salinity profile arrays"
    deallocate(S%new, S%old, S%grad_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Micro phytoplankton biomass profile array"
    deallocate(P%micro, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nano phytoplankton biomass profile array"
    deallocate(P%nano, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Pico phytoplankton biomass profile array"
    deallocate(P%pico, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Microzooplankton concentration profile array"
    deallocate(Z, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrate concentration profile array"
    deallocate(N%O, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ammonium concentration profile array"
    deallocate(N%H, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Silicon concentration profile array"
    deallocate(Si, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved organic nitrogen detritus concentration profile array"
    deallocate(D%DON, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Particulate organic nitrogen detritus concentration profile array"
    deallocate(D%PON, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Refractory nitrogen detritus concentration profile array"
    deallocate(D%refr, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Biogenic silicon detritus concentration profile array"
    deallocate(D%bSi, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_core_variables

end module core_variables
