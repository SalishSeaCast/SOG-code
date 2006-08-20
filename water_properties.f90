! $Id$
! $Source$

module water_properties
  ! Type definition, constants, functions, and subroutines to
  ! calculate water properties: specific heat capacity (Cp), thermal
  ! expansion coefficient (alpha), salinity expansion coefficient
  ! (beta), and density.
  !
  ! Public Type:
  ! 
  ! water_property -- Thermal and salinity expansion coefficients,
  !                   heat capacity, and density profiles.
  ! 
  ! Public Functions:
  !
  ! Cp_profile(T, S) -- Return the heat capacity of the water column at
  !                     the grid point depths based on the supplied water 
  !                     temperature and salinity profiles.
  ! 
  ! Public Subroutines:
  !
  ! alloc_water_props() -- Allocate memory for water property profiles.
  !
  ! dalloc_water_props() -- Deallocate memory from water property
  !                         profiles.

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Types:
       water_property, &
       ! Functions:
       Cp_profile, &
       ! Subroutines:
       alloc_water_props, &
       dalloc_water_props

  ! Public type definition:
  !
  ! Thermal and salinity expansion coefficients, heat capacity, and
  ! density profiles
  type :: water_property
     real(kind=dp), dimension(:), pointer :: &
          g,    &  ! Value at grid point depth
          gdiv, &  ! Rate of change (d/dz) at grid point depth
          i,    &  ! Value at grid cell interface depth
          idiv     ! Rate of change (d/dz) at grid cell interface depth
     ! Rate of change values are used only for expansion coefficients
  end type water_property

contains

  subroutine alloc_water_props(M, Cp)
    ! Allocate memory for water property profiles.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in)               :: M   ! Number of grid points
    type(water_property), intent(out) :: Cp  ! Heat capacity profile [J/kg.K]
    ! Local variables:
    integer              :: allocstat  ! Allocation return status
    character(len=80)    :: msg        ! Allocation failure message prefix

    msg = "Water heat capacity profile (Cp) arrays"
    allocate(Cp%g(0:M+1), Cp%i(0:M), stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_water_props


  subroutine dalloc_water_props(Cp)
    ! Deallocate memory for water property profiles.
    use malloc, only: dalloc_check
    implicit none
    ! Argument:
    type(water_property), intent(out) :: Cp  ! Heat capacity profile [J/kg.K]
    ! Local variables:
    integer              :: dallocstat  ! Allocation return status
    character(len=80)    :: msg         ! Allocation failure message prefix

    msg = "Water heat capacity profile (Cp) arrays"
    deallocate(Cp%g, Cp%i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_water_props


  function Cp_profile(T, S) result(Cp_g)
    ! Return the specific heat capacity of the water column in J/kg-K
    ! at the grid point depths based on the supplied water temperature
    ! and salinity profiles.
    !
    ! Based on Rich Pawlowicz's swcp.m Matlab function in his OCEANS
    ! toolbox (http://www.eos.ubc.ca/~rich/#OCEANS) with the pressure
    ! dependency removed.
    !
    ! Check value: Cp(10., 24.) = 4047.5 J/kg-K
    use unit_conversions, only: KtoC
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(in) :: T    ! Temperature [K]
    real(kind=dp), dimension(0:), intent(in) :: S    ! Salinity [-]
    ! Result:
    real(kind=dp), dimension(0:size(T) - 1) :: Cp_g  ! Specific heat capacity
    ! Local variables:
    real(kind=dp), dimension(0:size(T) - 1) :: TC, A, B, C

    ! Convert the temperature to Celsius
    TC = KtoC(T)
    ! Calculate the heat capacity profile at the grid point depths
    ! Specific heat Cp0 for P=0 (Millero, et al; UNESCO 1981)
    A = (-1.38385E-3 * TC + 0.1072763) * TC - 7.643575
    B = (5.148E-5 * TC - 4.07718E-3) * TC + 0.1770383
    C = (((2.093236E-5 * TC - 2.654387E-3) * TC + 0.1412855) &
         * TC - 3.720283) * TC + 4217.4
    Cp_g = (A + B * sqrt(S)) * S + C
  end function Cp_profile

end module water_properties
