! $Id$
! $Source$

module water_properties
  ! Type definition, variable declarations, and subroutines to
  ! calculate water properties: specific heat capacity (Cp), thermal
  ! expansion coefficient (alpha), salinity expansion coefficient
  ! (beta), and density (rho).
  !
  ! Public type:
  !
  ! water_property -- Density, thermal and salinity expansion
  !                   coefficients, and specific heat capacity profiles.
  !
  ! Public Variables:
  ! 
  ! Cp -- Specific heat capacity [J/kg.K]
  !
  ! rho -- Density [kg/m^3]
  ! 
  ! alpha -- 
  !
  ! beta --
  ! 
  ! Public Subroutines:
  !
  ! calc_rho_alpha_beta_Cp_profiles --
  !
  ! alloc_water_props -- Allocate memory for water property profiles.
  !
  ! dalloc_water_props -- Deallocate memory from water property
  !                       profiles.

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Type:
       water_property, &
       ! Variables:
       Cp, rho, &
       ! Subroutines:
       calc_rho_alpha_beta_Cp_profiles, &
       alloc_water_props, &
       dalloc_water_props

  ! Public type definition:
  !
  ! Thermal and salinity expansion coefficients, heat capacity, and
  ! density profiles
  type :: water_property
     real(kind=dp), dimension(:), pointer :: &
          g,    &   ! Value at grid point depth
          div_g, &  ! Rate of change (d/dz) at grid point depth
          i,    &   ! Value at grid cell interface depth
          div_i     ! Rate of change (d/dz) at grid cell interface depth
     ! Rate of change values are used only for density & expansion coefficients
  end type water_property

  ! Public variable declarations:
  !
  ! Thermal and salinity expansion coefficients, heat capacity, and
  ! density profiles
  type(water_property) :: &
       Cp, &  ! Specific heat capacity [J/kg.K]
       rho    ! Density [kg/m^3]

contains

  subroutine alloc_water_props(M)
    ! Allocate memory for water property profiles.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in)               :: M   ! Number of grid points
    ! Local variables:
    integer              :: allocstat  ! Allocation return status
    character(len=80)    :: msg        ! Allocation failure message prefix

    msg = "Water heat capacity profile (Cp) arrays"
    allocate(Cp%g(0:M+1), Cp%i(0:M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Water density profile arrays"
    allocate(rho%g(0:M+1), rho%i(0:M), rho%div_g(1:M), rho%div_i(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_water_props


  subroutine dalloc_water_props
    ! Deallocate memory for water property profiles.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer              :: dallocstat  ! Allocation return status
    character(len=80)    :: msg         ! Allocation failure message prefix

    msg = "Water specific heat capacity profile (Cp) arrays"
    deallocate(Cp%g, Cp%i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Water density profile arrays"
    deallocate(rho%g, rho%i, rho%div_g, rho%div_i, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_water_props


  subroutine calc_rho_alpha_beta_Cp_profiles(T, S, Cp_g)
    ! Calculate the profiles of density, thermal expansion
    ! coefficient, salinity expansion coefficient, and specific heat
    ! capacity through the water column at the grid point depths based
    ! on the supplied water temperature and salinity profiles.
    !
    ! Based on Rich Pawlowicz's Matlab function in his OCEANS
    ! toolbox (http://www.eos.ubc.ca/~rich/#OCEANS) with the pressure
    ! dependency removed.  rho, alpha, and beta calculations are based
    ! on xxx.m.  Cp calculation is based on swcp.m.  In all cases the
    ! pressure dependency has been removed.
    !
    ! The calculations have been conbined into 1 subroutine because
    ! they all require the evaluation of the square root of the
    ! salinity profile, and this way we can minimize the number of
    ! times that CPU-cycle intensive operation is done.
    !
    ! Check values: 
    !   Cp(10., 24.) = 4047.5 J/kg-K
    use unit_conversions, only: KtoC
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(in) :: T    ! Temperature [K]
    real(kind=dp), dimension(0:), intent(in) :: S    ! Salinity [-]
    ! Result:
    real(kind=dp), dimension(0:size(T) - 1) :: Cp_g  ! Specific heat capacity
    ! Local variables:
    real(kind=dp), dimension(0:size(T) - 1) :: TC, sqrtS, A, B, C

    ! Convert the temperature to Celsius
    TC = KtoC(T)
    ! Calculate the square root of the salinities
    sqrtS = sqrt(S)
    ! Calculate the heat capacity profile at the grid point depths
    ! Specific heat Cp0 for P=0 (Millero, et al; UNESCO 1981)
    A = (-1.38385E-3 * TC + 0.1072763) * TC - 7.643575
    B = (5.148E-5 * TC - 4.07718E-3) * TC + 0.1770383
    C = (((2.093236E-5 * TC - 2.654387E-3) * TC + 0.1412855) &
         * TC - 3.720283) * TC + 4217.4
    Cp_g = (A + B * sqrtS) * S + C
  end subroutine calc_rho_alpha_beta_Cp_profiles

end module water_properties
