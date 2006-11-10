! $Id$
! $Source$

module diffusion

  implicit none

  private
  public :: diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux

contains

  subroutine diffusion_coeff(grid, dt, K_all, &  ! in
       Bmatrix)                                  ! out
    ! Calculate the strength of the diffusion coefficients
    ! and put them in diagonal vectors of a tridiagonal matrix.
    !
    ! This is an interim step on the way to the calculation of the A
    ! matrices in Large, et al (1994) App. D.  The difference is that
    ! we are using the IMEX scheme to solve the semi-implicit
    ! equations, and some additional manipulation of the diffusion
    ! coefficients is required to facilitate that.
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(1:), intent(in) :: &
         K_all  ! Total diffusion coefficient array
    type(tridiag) :: &
         Bmatrix  ! Diffusion coefficient matrix
    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: Omega_minus, Omega_plus
    
    ! Calculate Omega+ and Omega-
    ! (Large, et al (1994), pg. 398, eqn D9)
    Omega_minus(2:) = (dt / grid%i_space(2:)) * &
         (K_all(1:grid%M-1) / grid%g_space(1:grid%M-1))
    Omega_plus = (dt / grid%i_space) * &
         (K_all(1:grid%M) / grid%g_space(1:grid%M))

    ! Initialize diagonal vectors
    Bmatrix%sub = 0.
    Bmatrix%diag = 0.
    Bmatrix%sup = 0.
    ! Put in the combinations of Omega+ and Omega- terms (eqn D9)
    Bmatrix%diag(1) = -Omega_plus(1)
    Bmatrix%sub(2:) = Omega_minus(2:)
    Bmatrix%diag(2:) = -Omega_plus(2:) - Omega_minus(2:)
    Bmatrix%sup = Omega_plus
    Bmatrix%sup(grid%M) = 0.
  end subroutine diffusion_coeff


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
    real(kind=dp), dimension(1:), intent(in):: Kall ! total mixing K
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
      real(kind=dp), dimension(1:), intent(in) :: &
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
