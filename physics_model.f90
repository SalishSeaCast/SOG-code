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

  subroutine double_diffusion(mm, T_grad_i, sal, alpha_i, beta_i, nus)
    ! *** What's it do?
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use mean_param, only: prop, Knu, div_interface
!!$    use surface_forcing
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: mm  ! *** Replace this with M
    real(kind=dp), dimension(1:mm%M), intent(in) :: &
         T_grad_i  ! Temperature gradient profile at grid interface depths
    type(prop), intent(inout) :: sal
    real(kind=dp), dimension(0:mm%M), intent(in) :: &
         alpha_i, &  ! Thermal expansion coeff profile at interface depths
         beta_i      ! Salinity contraction coeff profile at interface depths
    ! *** Changed this to K%s%dd, and K%t%dd as output args
    type(Knu), intent(inout) :: nus

    ! Local variables:
    real(kind=dp), dimension(1:mm%M) :: Ri_rho
    integer :: k  ! Loop index over depth
    real(kind=dp), parameter :: &
         nu_f = 0.001, &    ! Salt fingering constant [m^2/s]
         Ri_rho_o = 1.9, &  ! *** Something to do with salt fingering
         p_2 = 3, &         ! Salt finger diff power constant
         nu = 1.5d-06       ! Molecular viscosity [m^2/s]

    CALL div_interface(mm, sal)

!!$    nus%s%dd = 0.0
!!$    nus%t%dd = 0.0

      ! Calculate double diffusion denisty ratio (30)
      Ri_rho = alpha_i(1:) * T_grad_i / (beta_i(1:) * sal%div_i)
      
      ! *** May be able to use a forall construct here
      do k = 1, mm%M
         if (1.0 < Ri_rho(k) &
              .and. Ri_rho(k) < 2.0 &
              .and. alpha_i(k) * T_grad_i(k) > 0 &
              .and. beta_i(k) * sal%div_i(k) > 0.) then  
            ! Salt fingering
            if (1.0 < Ri_rho(k) .and. Ri_rho(k) < Ri_rho_o) then
               nus%s%dd(k) = nu_f * (1.0 - ((Ri_rho(k) - 1.0) &
                    / (Ri_rho_o - 1.0)) ** 2) ** p_2   ! (31a)
               nus%t%dd(k) = 0.7 * nus%s%dd(k) ! (31c)            
            endif ! (31b) not needed, already set to 0
         else if (0. < Ri_rho(k) &
              .and. Ri_rho(k) < 1.0 &
              .and. alpha_i(k) * T_grad_i(k) < 0. &
              .and. beta_i(k) * sal%div_i(k) < 0.) then   
            ! Diffusive instability
            nus%t%dd(k) = nu * 0.909 &
                 * exp(4.6 * exp(-0.54 * (1.0 / Ri_rho(k) - 1.0))) !(32)
            if (Ri_rho(k) >= 0.5 .and. Ri_rho(k) < 1.0) then
               nus%s%dd(k) = nus%t%dd(k) &
                    * (1.85 - 0.85 / Ri_rho(k)) * Ri_rho(k) !(34)
            else if (Ri_rho(k) < 0.5) then
               nus%s%dd(k) = nus%t%dd(k) * 0.15 * Ri_rho(k) ! (34)
            endif
         else
            nus%t%dd(k) = 0.0
            nus%s%dd(k) = 0.0
         endif
      enddo
    end subroutine double_diffusion

end module physics_model
