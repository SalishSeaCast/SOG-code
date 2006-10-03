! $Id$
! $Source$

module physics_model
  ! Type declarations, variables, and subroutines related to the physics
  ! model in the code.
  ! 
  ! Public Variables:
  !
  !
  ! Public Subroutines:
  !

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Subroutines:
       double_diffusion

  ! Private parameter declarations:

contains

  subroutine double_diffusion(M, T_grad_i, S_grad_i, alpha_i, beta_i, &
       nu_t_dd, nu_s_dd)
    ! Calculate double diffusion mixing.  
    ! See Large, etal (1994), pp373-374, equations 30 to 34.
    use precision_defs, only: dp
    use io_unit_defs, only: stderr
    implicit none
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    real(kind=dp), dimension(1:), intent(in) :: &
         T_grad_i, &  ! Temperature gradient profile at grid interface depths
         S_grad_i     ! Salinity gradient profile at grid interface depths
    real(kind=dp), dimension(0:), intent(in) :: &
         alpha_i, &  ! Thermal expansion coeff profile at interface depths
         beta_i      ! Salinity contraction coeff profile at interface depths
    real(kind=dp), dimension(1:), intent(out) :: &
         nu_t_dd, &  ! Thermal diffusivity for double diffusion mixing
         nu_s_dd     ! Saline diffusivity for double diffusion mixing

    ! Local variables:
    real(kind=dp), dimension(1:M) :: R_rho  ! Double diffusion density ratio
    integer :: k  ! Loop index over depth
    ! Double diffusion model parameters (see Large, etal (1994) pg 374, fig. 4)
    real(kind=dp), parameter :: &
         nu_f = 10.0d-4, &  ! Maximum salt fingering diffusivity [m^2/s]
         R_rho_o = 1.9, &   ! Salt fingering limit (nu_s = 0)
         p_2 = 3, &         ! Salt finger diff power constant
         nu = 1.5d-06       ! Molecular diffusivity (viscosity) [m^2/s]

    ! Calculate double diffusion density ratio (Large, et al, (1994) eq'n 30)
    R_rho = alpha_i(1:) * T_grad_i / (beta_i(1:) * S_grad_i)
      
    do k = 1, M
       ! Determine if there is double diffision, and if so, what type,
       ! salt fingering, or diffusive convection
       if (1.0 < R_rho(k) &
            .and. R_rho(k) < 2.0 &
            .and. alpha_i(k) * T_grad_i(k) > 0 &
            .and. beta_i(k) * S_grad_i(k) > 0.) then  
          ! Salt fingering
          if (1.0 < R_rho(k) .and. R_rho(k) < R_rho_o) then
             ! Large, et al, (1994) eq'n 31a
             nu_s_dd(k) = nu_f &
                  * (1.0 - ((R_rho(k) - 1.0) &
                             / (R_rho_o - 1.0)) ** 2) ** p_2
          else  ! R_rho >= R_rho_o
             ! Large, et al, (1994) eq'n 31b
             nu_s_dd = 0.
          endif
          ! Large, et al, (1994) eq'n 31c
          nu_t_dd(k) = 0.7 * nu_s_dd(k)
       else if (0. < R_rho(k) &
            .and. R_rho(k) < 1.0 &
            .and. alpha_i(k) * T_grad_i(k) < 0. &
            .and. beta_i(k) * S_grad_i(k) < 0.) then   
          ! Diffusive convection
          ! Large, et al, (1994) eq'n 32
          nu_t_dd(k) = nu * 0.909 &
               * exp(4.6 * exp(-0.54 * (1.0 / R_rho(k) - 1.0)))
          ! Large, et al, (1994) eq'n 34
          if (R_rho(k) >= 0.5 .and. R_rho(k) < 1.0) then
             nu_s_dd(k) = nu_t_dd(k) &
                  * (1.85 - 0.85 / R_rho(k)) * R_rho(k) !(34)
          else if (R_rho(k) < 0.5) then
             nu_s_dd(k) = nu_t_dd(k) * 0.15 * R_rho(k) ! (34)
          endif
       else
          ! No double diffusion (i.e. temperature and salinity are
          ! both stabilizing, or water column is convectively
          ! unstable)
          nu_t_dd(k) = 0.0
          nu_s_dd(k) = 0.0
       endif
    enddo
  end subroutine double_diffusion

end module physics_model
