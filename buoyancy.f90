! $Id$
! $Source$

subroutine buoyancy(grid, Tnew, Snew, hml, Itotal, F_n, wb0, rho, &  ! in
     alpha, beta, Cp, Fw_surface,                                 &  ! in
     Bnew, Bf)                                                       ! out
  ! Calculate the buoyancy profile, and the surface buoyancy forcing.
  use precision_defs, only: dp
  use grid_mod, only: grid_, interp_g_d, interp_i_d
  use mean_param, only: height
  use surface_forcing, only: g
  implicit none
  ! Arguments:
  type(grid_), intent(in) :: grid  ! Grid depths & spacings arrays
  real(kind=dp), dimension(0:grid%M+1), intent(in) :: &
       Tnew, &  ! Temperature profile
       Snew     ! Salinity profile
  type(height), intent(in) :: hml  ! Mixing layer depth
  real(kind=dp), dimension(0:grid%M), intent(in) :: &
       Itotal, &  ! Irradiance
       F_n        ! Fresh water contribution to salinity flux
  real(kind=dp), intent(in) :: wb0  ! Surface buoyance flux
  real(kind=dp), dimension(0:grid%M+1), intent(in) :: &
       rho,   &  ! Density
       alpha, &  ! Thermal expansion coefficient
       beta,  &  ! Salinity expansion coefficient
       Cp        ! Specific heat capacity
  logical, intent(in) :: Fw_surface
  ! Results:
  real(kind=dp), dimension(0:grid%M+1), intent(out) :: Bnew  ! Buoyance profile
  real(kind=dp), intent(out) :: Bf  ! Surface buoyancy forcing

  ! Local variables:
  ! Water properties, irradiance, and fresh water flux at mixing layer
  ! depth
  real(kind=dp) :: rho_ml, alpha_ml, beta_ml, Cp_ml, I_ml, FN_ml, &
       ! Contributions to surface buoyancy forcing
       Br, &  ! Radiative
       Bfw    ! Fresh water flux

  ! Buoyancy profile (see Large, etal (1994), eq'n A3a)
  Bnew = g * (alpha * Tnew - beta * Snew)

  ! Interpolate to find water property and flux values at the mixing
  ! layer depth
  rho_ml = interp_g_d(rho, hml%new)
  alpha_ml = interp_g_d(alpha, hml%new)
  beta_ml = interp_g_d(beta, hml%new)
  Cp_ml = interp_g_d(Cp, hml%new)
  ! *** Irradiance and fresh water flux interpolations could be
  ! replaced with a function analogous to interp_d
  I_ml = Itotal(hml%i-1) + (Itotal(hml%i) - Itotal(hml%i-1)) &
       * (hml%new - grid%d_i(hml%i-1)) / grid%i_space(hml%i)
  Fn_ml = F_n(hml%i-1) + (F_n(hml%i) - F_n(hml%i-1)) &
       * (hml%new - grid%d_i(hml%i-1)) / grid%i_space(hml%i)

  ! Contribution of heat and fresh water fluxes to surface buoyancy
  ! forcing (Surface values less the amount at (and thus below) the
  ! mixing layer depth)
  !
  ! Radiative contribution (see Large, etal (1994) eq'n A3c)
  Br = g * (alpha(0) * Itotal(0) / (Cp(0) * rho(0)) &
            - alpha_ml * I_ml / (rho_ml * Cp_ml))
  ! Fresh water salinity flux contribution
  if (Fw_surface) then
     Bfw = 0.
  else
     Bfw = g * (beta(0) * F_n(0) - beta_ml * Fn_ml)
  endif
  ! Calculate surface buoyancy forcing (an extension of Large, etal
  ! (1994) eq'n A3d
  Bf = -wb0 + Br + Bfw
end subroutine buoyancy
