! $Id$
! $Source$

module diffusion
  ! Type definitions & variable declarations, and subroutines related
  ! to the diffusion terms in the semi-implicit PDEs in the SOG code.
  !
  ! Public Subroutines:
  !
  !   diffusion_coeff -- Calculate the strength of the diffusion
  !                      coefficients and put them in diagonal vectors
  !                      of a tridiagonal matrix.  

  implicit none

  private
  public :: &
       ! Subroutines:
       diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux

contains

  subroutine diffusion_coeff(dt, K, &  ! in
       Bmatrix)                        ! out
    ! Calculate the strength of the diffusion coefficients
    ! and put them in diagonal vectors of a tridiagonal matrix.
    !
    ! This is an interim step on the way to the calculation of the A
    ! matrices in Large, et al (1994) App. D.  The difference is that
    ! we are using the IMEX scheme to solve the semi-implicit
    ! equations, and some additional manipulation of the diffusion
    ! coefficients is required to facilitate that.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use numerics, only: tridiag  ! Tridiagonal matrix arrays type def'n
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and arrays
    
    implicit none
    
    ! Arguments:
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(1:), intent(in) :: &
         K  ! Overall diffusivity profile
    type(tridiag), intent(inout) :: &
         Bmatrix  ! Diffusion coefficient matrix

    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: &
         Omega_minus, Omega_plus  ! Components of diffusion matrix
                                  ! elements (Large, et al (1994), eqn
                                  ! D9).
    integer :: M  ! Number of grid points (to improve readability of
                  ! vectorized code).
    
    ! Calculate Omega+ and Omega-
    ! (Large, et al (1994), pg. 398, eqn D9)
    M = grid%M
    Omega_minus(2:) = (dt / grid%i_space(2:)) * &
         (K(1:M-1) / grid%g_space(1:M-1))
    Omega_plus = (dt / grid%i_space) * &
         (K(1:M) / grid%g_space(1:M))

    ! Initialize diagonal vectors
    Bmatrix%sub = 0.
    Bmatrix%diag = 0.
    Bmatrix%sup = 0.
    ! Put in the combinations of Omega+ and Omega- terms (eqn D9)
    Bmatrix%diag(1) = -Omega_plus(1)
    Bmatrix%sub(2:) = Omega_minus(2:)
    Bmatrix%diag(2:) = -Omega_plus(2:) - Omega_minus(2:)
    Bmatrix%sup = Omega_plus
    Bmatrix%sup(M) = 0.
  end subroutine diffusion_coeff


  subroutine diffusion_nonlocal_fluxes(dt, K, flux, Bf, surface_flux, &
       distrib_flux, bottom_value, RHS)
    ! Calculate the diffusive flux into the bottom of the domain,
    ! non-local transport fluxes through the domain, and other
    ! transport flux through the domain.
    !
    ! This is an interim step on the way to the calculation of the H
    ! vectors in Large, et al (1994) App. D.  The difference is that
    ! we are using the IMEX scheme to solve the semi-implicit
    ! equations, and some additional manipulation of the RHS vectors
    ! is required to facilitate that.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Functions:
    use turbulence, only: nonlocal_scalar_transport
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and arrays
    use mixing_layer, only: &
         h  ! Mixing layer depth and indices

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(1:), intent(in) :: &
         K  ! Overall diffusivity profile
    real(kind=dp), intent(in) :: &
         flux, &  ! Flux to calculate the non-local transport term
                  ! value for.
         Bf       ! Surface bouyancy flux
    real(kind=dp), intent(in) :: &
         surface_flux  ! Surface flux of quantity
    real(kind=dp), dimension(0:), intent(in):: &
         distrib_flux  ! Distributed flux profile
    real(kind=dp), intent(in) :: &
         bottom_value  ! Value of quantity at bottom of grid
    real(kind=dp), dimension(1:), intent(out) :: &
         RHS  ! RHS term vector for semi-implicit diffusion/advection PDE

    ! Local variables:
    real(kind=dp) :: &
         gamma_j, gamma_jm1  ! Non-local transport term values at the
                             ! present grid layer interface depth, and
                             ! the one above.
    integer :: &
         j  ! Index over the grid depth

    ! Initialize RHS vector
    RHS = 0.0d0
    ! Calculate RHS (the H vector) (Large, et al (1994), eqn D10)
    gamma_j = nonlocal_scalar_transport(flux, h%new, grid%d_i(1), Bf)
    ! Surface
    RHS(1) = dt / grid%i_space(1) &
         * (K(1) * gamma_j - surface_flux + distrib_flux(1) - distrib_flux(0))
    ! Interior
    do j = 2, grid%M
       gamma_jm1 = gamma_j
       gamma_j = nonlocal_scalar_transport(flux, h%new, grid%d_i(j), Bf)
       RHS(j) = dt / grid%i_space(j) &
            * (K(j) * gamma_j - K(j-1) * gamma_jm1 &
            + distrib_flux(j) - distrib_flux(j-1))
    enddo
    ! Bottom of grid
    RHS(grid%M) = dt / grid%i_space(grid%M) * &
         (distrib_flux(grid%M) - distrib_flux(grid%M-1) &
          + bottom_value * K(grid%M) / grid%g_space(grid%M))      
  end subroutine diffusion_nonlocal_fluxes


  subroutine diffusion_bot_surf_flux(dt, K, surface_flux, &
       bottom_value, RHS)
    ! For variables without distributed fluxes and non-local fluxes.
    ! Initialize the RHS vector, and calculate the diffusive fluxes
    ! into the bottom of the grid and at the surface.
    !
    ! This is an interim step on the way to the calculation of the H
    ! vectors in Large, et al (1994) App. D.  The difference is that
    ! we are using the IMEX scheme to solve the semi-implicit
    ! equations, and some additional manipulation of the RHS vectors
    ! is required to facilitate that.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and arrays

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(1:), intent(in) :: &
         K  ! Total diffusion coefficient array
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
         (bottom_value * K(grid%M) / grid%g_space(grid%M))      
  end subroutine diffusion_bot_surf_flux

end module diffusion
