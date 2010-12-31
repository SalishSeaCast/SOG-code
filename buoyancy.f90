module buoyancy
  ! Variable declarations, and subroutines related to the buoyancy
  ! profile and surface buoyancy flux calculation in the SOG code.
  !
  ! Public Variables:
  !
  !   Bf -- Surface buoyancy flux
  !
  ! Diagnostic Variables (not used anywhere, but available for output):
  !
  !   B -- Buoyancy profile

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Variables:
       Bf, &  ! Surface buoyancy forcing
       ! Diagnostics:
       B, &  ! Buoyancy profile array
       ! Subroutines:
       alloc_buoyancy_variables, calc_buoyancy, dalloc_buoyancy_variables

  ! Variable Declarations:
  !
  ! Public:
  real(kind=dp) :: &
       Bf  ! Surface buoyancy forcing
  ! Diagnostic:
  real(kind=dp), dimension(:), allocatable ::&
       B  ! Buoyancy profile array

contains

subroutine calc_buoyancy(Tnew, Snew, hml, Itotal, rho, alpha, beta, Cp)
  ! Calculate the buoyancy profile, surface turbulent kinematic
  ! buoyancy flux, and the surface buoyancy forcing.

  ! Elements from other modules:
  ! Parameter values:
  use fundamental_constants, only: g
  ! Subroutines:
  use grid_mod, only: interp_value
  ! Variables:
  use grid_mod, only: &
       grid  ! Grid parameters and depth & spacing arrays
  use turbulence, only: &
       wbar  ! Turbulent kinematic flux profile arrays; we need the
             ! surface temperature and salinity values (wbar%t(0) &
             ! wbar%s(0)), and we set the surface buoyancy value
             ! wbar%b(0).
  use freshwater, only: &
       Fw_surface, &  ! Add all of the fresh water on the surface?
       F_n            ! Fresh water contribution to salinity flux

  implicit none

  ! Arguments:
  real(kind=dp), dimension(0:), intent(in) :: &
       Tnew, &   ! Temperature profile
       Snew      ! Salinity profile
  real(kind=dp), intent(in) :: &
       hml       ! Mixing layer depth
  real(kind=dp), dimension(0:), intent(in) :: &
       Itotal    ! Irradiance
  real(kind=dp), dimension(0:), intent(in) :: &
       rho,   &  ! Density profile
       alpha, &  ! Thermal expansion coefficient profile
       beta,  &  ! Saline expansion coefficient profile
       Cp        ! Specific heat capacity profile

  ! Local variables:
  ! Water properties, irradiance, and fresh water flux at mixing layer
  ! depth
  real(kind=dp) :: &
       rho_ml,   &  ! Density at mixing layer depth
       alpha_ml, &  ! Thermal expansion coefficient at mixing layer depth
       beta_ml,  &  ! Salinity contraction coefficient at mixing layer depth
       Cp_ml,    &  ! Specific heat capacity at mixing layer depth
       I_ml,     &  ! Total irradiance at mixing layer depth
       Fn_ml,    &  ! Fresh water contraction to salinity flux at
                    ! mixing layer depth
       Br,       &  ! Radiative contributions to surface buoyancy forcing
       Bfw          ! Fresh water flux contributions to surface buoyancy forcing
  integer :: j_junk ! Unused index returned by interp_value() subroutine

  ! Calculate buoyancy profile (Large, et al (1994), eqn A3a)
  B = g * (alpha * Tnew - beta * Snew)
  ! Interpolate to find water property and flux values at the mixing
  ! layer depth
  call interp_value(hml, 0, grid%d_g, rho, rho_ml, j_junk)
  call interp_value(hml, 0, grid%d_g, alpha, alpha_ml, j_junk)
  call interp_value(hml, 0, grid%d_g, beta, beta_ml, j_junk)
  call interp_value(hml, 0, grid%d_g, Cp, Cp_ml, j_junk)
  call interp_value(hml, 0, grid%d_i, Itotal, I_ml, j_junk)
  call interp_value(hml, 0, grid%d_i, F_n, Fn_ml, j_junk)

  ! Calculate contribution of heat and fresh water fluxes to surface
  ! buoyancy forcing (Surface values less the amount at (and thus
  ! below) the mixing layer depth)
  !
  ! Radiative contribution (see Large, etal (1994) eq'n A3c)
  Br = g * (alpha(0) * Itotal(0) / (rho(0) * Cp(0)) &
            - alpha_ml * I_ml / (rho_ml * Cp_ml))
  ! Fresh water salinity flux contribution
  if (Fw_surface) then
     Bfw = 0.
  else
     Bfw = g * (beta(0) * F_n(0) - beta_ml * Fn_ml)
  endif
  ! Calculate the surface turbulent buoyancy flux (Large, et al
  ! (1994), eqn A3b).
  wbar%b(0) = g * (alpha(0) * wbar%t(0) - beta(0) * wbar%s(0))
  ! Calculate surface buoyancy forcing (an extension of Large, et al
  ! (1994) eqn A3d)
  Bf = -wbar%b(0) + Br + Bfw
end subroutine calc_buoyancy


  subroutine alloc_buoyancy_variables(M)
    ! Allocate memory for buoyancy variables arrays.
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
  end subroutine alloc_buoyancy_variables


  subroutine dalloc_buoyancy_variables
    ! Deallocate memory from buoyancy variables arrays.

    ! Subroutines from other modules:
    use malloc, only: dalloc_check
    
    implicit none
    
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Buoyancy profile array"
    deallocate(B, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_buoyancy_variables

end module buoyancy
