! $Id$
! $Source$

module diffusion

  implicit none

  private
  public :: diffusion_coeff, new_diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux

contains

  subroutine diffusion_coeff(grid, dt, K_all, Bmatrix)
    ! Calculate the strength of the diffusion coefficients
    ! and put them in diagonal vectors of a tridiagonal matrix.
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use mean_param, only: trivector
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         K_all  ! Total diffusion coefficient array
    type(trivector), intent(out) :: &
         Bmatrix
    ! Local variables:
    integer :: index             ! counter through depth
    real(kind=dp), dimension(grid%M) :: O_minus, O_plus
    
    ! Calculate Omega+ and Omega-
    ! (Large, et al (1994), pg. 398, eqn D9)
    do index = 1, grid%M
       O_minus(index) = dt / grid%i_space(index) * &
            K_all(index-1) / grid%g_space(index-1)
       O_plus(index) = dt / grid%i_space(index) * &
            K_all(index) / grid%g_space(index)
    enddo

    ! Initialize diagonal vectors
    Bmatrix%A = 0.
    Bmatrix%B = 0.
    Bmatrix%C = 0.
    ! Put in the combinations of Omega+ and Omega- terms (eqn D9)
    Bmatrix%B = -O_plus -O_minus
    Bmatrix%C = O_plus
    Bmatrix%A = O_minus
    Bmatrix%C(grid%M) = 0.
  end subroutine diffusion_coeff


  subroutine new_diffusion_coeff(grid, dt, K_all, sub, diag, sup)
    ! Calculate the strength of the diffusion coefficients
    ! and put them in diagonal vectors of a tridiagonal matrix.
    use precision_defs, only: dp
    use grid_mod, only: grid_
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         K_all  ! Total diffusion coefficient array
    real(kind=dp), dimension(1:), intent(out) :: &
         sub,  &  ! Sub-diagonal vector of diffusion coeff matrix
         diag, &  ! Diagonal vector of diffusion coeff matrix
         sup      ! Super-diagonal vector of diffusion coeff matrix
    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: O_minus, O_plus
    
    ! Calculate Omega+ and Omega-
    ! (Large, et al (1994), pg. 398, eqn D9)
    O_minus = (dt / grid%i_space) * &
         (K_all(0:grid%M-1) / grid%g_space(0:grid%M-1))
    O_plus = (dt / grid%i_space) * &
         (K_all(1:grid%M) / grid%g_space(1:grid%M))

    ! Initialize diagonal vectors
    sub = 0.
    diag = 0.
    sup = 0.
    ! Put in the combinations of Omega+ and Omega- terms (eqn D9)
    sub = O_minus
    diag = -O_plus - O_minus
    sup = O_plus
    sup(grid%M) = 0.
  end subroutine new_diffusion_coeff


  subroutine diffusion_nonlocal_fluxes (grid, dt, Kall, gamma, surface_flux, &
       distrib_flux, bottom_value, Gvector)
    ! The subroutine calculates the diffusive flux into the bottom of
    ! the domain, non-local transport fluxes through the domain and
    ! and other transport flux through the domain.
    use precision_defs, only: dp
    use grid_mod, only: grid_

    implicit none

    ! Arguments:
    type(grid_),intent(in):: grid
    real(kind=dp), intent(in):: dt
    real(kind=dp), dimension(0:), intent(in):: Kall ! total mixing K
    real(kind=dp), dimension(0:), intent(in):: gamma ! nonlocal effect
    real(kind=dp), intent(in) :: surface_flux
    real(kind=dp), dimension(0:), intent(in):: distrib_flux
    real(kind=dp), intent(in) :: bottom_value ! of scalar
    real(kind=dp), dimension(1:), intent(out) :: Gvector
    ! Local variables:
    integer :: index                   ! counter through depth

    ! initialize Gvector to 0
    Gvector = 0.

    ! calculate rhs (the H vector) (D10)
    Gvector(1) = dt / grid%i_space(1) * &
         (Kall(1) * gamma(1) - surface_flux + distrib_flux(1)- distrib_flux(0))

    do index = 2, grid%M - 1
       Gvector(index) = dt / grid%i_space(index) * &
            (Kall(index) * gamma(index) - Kall(index-1) * gamma(index-1) &
             + distrib_flux(index) - distrib_flux(index-1))
    end do
    Gvector(grid%M) = dt / grid%i_space(grid%M) * &
         (distrib_flux(grid%M) - distrib_flux(grid%M-1) &
          + bottom_value * Kall(grid%M) / grid%g_space(grid%M))      
  end subroutine diffusion_nonlocal_fluxes


    subroutine diffusion_bot_surf_flux(grid, dt, K_all, surface_flux, &
       bottom_value, RHS)
      ! For variables without distributed fluxes and non-local fluxes.
      ! Initialize the RHS vector, and calculate the diffusive fluxes
      ! into the bottom of the grid and at the surface.
      use precision_defs, only: dp
      use grid_mod, only: grid_
      implicit none

      ! Arguments:
      type(grid_), intent(in) :: &
           grid  ! Grid arrays
      real(kind=dp), intent(in) :: &
           dt  ! Time step [s]
      real(kind=dp), dimension(0:), intent(in) :: &
           K_all  ! Total diffusion coefficient array
      real(kind=dp), intent(in) :: &
           surface_flux, &  ! Surface flux of quantity
           bottom_value     ! Value of quantity at bottom of grid
      real(kind=dp), dimension(1:), intent(out) :: &
           RHS  ! RHS term vector for semi-implicit diffusion/advection PDE

      ! Initialize the diffusion/advection vector of the right-hand
      ! side of the semi-implicit equations
      RHS = 0.

      ! Calculate grid boundary diffusive fluxes
      ! (Large, et al (1994), pg 398, eqn D10c)
      ! Surface
      RHS(1) = - dt / grid%i_space(1) * surface_flux
      ! Bottom of grid
      RHS(grid%M) = dt / grid%i_space(grid%M) * &
           (bottom_value * K_all(grid%M) / grid%g_space(grid%M))      
    end subroutine diffusion_bot_surf_flux

end module diffusion
