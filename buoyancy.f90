! $Id$
! $Source$

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
!!$       B, &  ! Buoyancy profile
       ! Subroutines:
!!$       alloc_buoyancy_variables, &
       calc_buoyancy
!!$, &
!!$       dalloc_bouyancy variables

  ! Variable Declarations:
  !
  ! Public:
  real(kind=dp) :: &
       Bf  ! Surface buoyancy forcing

contains

subroutine calc_buoyancy(grid, Tnew, Snew, hml, Itotal, F_n, rho, &  ! in
     alpha, beta, Cp, Fw_surface,                            &  ! in
     B)                                                     ! out
  ! Calculate the buoyancy profile, and the surface buoyancy forcing.
  use precision_defs, only: dp
  use grid_mod, only: grid_, interp_value
  use turbulence, only: &
       wbar  ! Turbulent kinematic flux profile arrays; we need wbar%b(0)
  use fundamental_constants, only: g
  implicit none
  ! Arguments:
  type(grid_), intent(in) :: grid  ! Grid depths & spacings arrays
  real(kind=dp), dimension(0:grid%M+1), intent(in) :: &
       Tnew, &  ! Temperature profile
       Snew     ! Salinity profile
  real(kind=dp), intent(in) :: hml  ! Mixing layer depth
  real(kind=dp), dimension(0:grid%M), intent(in) :: &
       Itotal, &  ! Irradiance
       F_n        ! Fresh water contribution to salinity flux
  real(kind=dp), dimension(0:grid%M+1), intent(in) :: &
       rho,   &  ! Density
       alpha, &  ! Thermal expansion coefficient
       beta,  &  ! Salinity expansion coefficient
       Cp        ! Specific heat capacity
  logical, intent(in) :: Fw_surface
  ! Results:
  real(kind=dp), dimension(0:grid%M+1), intent(out) :: B  ! Buoyance profile

  ! Local variables:
  ! Water properties, irradiance, and fresh water flux at mixing layer
  ! depth
  real(kind=dp) :: rho_ml, alpha_ml, beta_ml, Cp_ml, I_ml, FN_ml, &
       ! Contributions to surface buoyancy forcing
       Br, &  ! Radiative
       Bfw    ! Fresh water flux
  integer :: j_below ! Index of grid point or interface immediately
                     ! below mixing layer depth

  ! Calculate buoyancy profile (Large, etal (1994), eq'n A3a)
  B = g * (alpha * Tnew - beta * Snew)

  ! Interpolate to find water property and flux values at the mixing
  ! layer depth
  call interp_value(hml, 0, grid%d_g, rho, rho_ml, j_below)
  call interp_value(hml, 0, grid%d_g, alpha, alpha_ml, j_below)
  call interp_value(hml, 0, grid%d_g, beta, beta_ml, j_below)
  call interp_value(hml, 0, grid%d_g, Cp, Cp_ml, j_below)
  call interp_value(hml, 0, grid%d_i, Itotal, I_ml, j_below)
  call interp_value(hml, 0, grid%d_i, F_n, Fn_ml, j_below)

  ! Contribution of heat and fresh water fluxes to surface buoyancy
  ! forcing (Surface values less the amount at (and thus below) the
  ! mixing layer depth)
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
  ! Calculate surface buoyancy forcing (an extension of Large, etal
  ! (1994) eq'n A3d
  Bf = -wbar%b(0) + Br + Bfw
end subroutine calc_buoyancy

end module buoyancy
