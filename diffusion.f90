module diffusion

  implicit none

  private
  public :: diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux

contains

  subroutine diffusion_coeff (grid, dt, Kall, Bx)

    ! Calculates the strength of the diffusion coefficients (K, dt and grid
    ! spacing) for diffusion and puts them in Bx

    use mean_param, only : gr_d, trivector

    implicit none

    ! arguments
    type(gr_d), intent(in) :: grid      
    double precision, intent(in) :: dt
    double precision, dimension(0:),intent(in)::Kall ! total mixing K

    type(trivector), intent(out)::Bx

    ! local variables
    integer :: index             ! counter through depth
    double precision, dimension(grid%M)::O_minus, O_plus
    
! Calculate Omega plus and minus below (D9)

    DO index = 1, grid%M
       O_minus(index) = dt / grid%i_space(index) * &
            Kall(index-1) / grid%g_space(index-1)
       O_plus(index) = dt / grid%i_space(index) * &
            Kall(index) / grid%g_space(index)
    END DO

    Bx%A = 0.
    Bx%B = 0.
    Bx%C = 0.
    
    ! put in three terms.(D9)

    Bx%B = -O_plus -O_minus
    Bx%C = O_plus
    Bx%A = O_minus
    Bx%C(grid%M) = 0.

  end subroutine diffusion_coeff

  subroutine diffusion_nonlocal_fluxes (grid, dt, Kall, gamma, surface_flux, &
       distrib_flux, bottom_value, Gvector)

! The subroutine calculates the diffusive flux into the bottom of the domain,
! non-local transport fluxes through the domain and and other transport
! flux through the domain.

    use mean_param, only: gr_d

    implicit none

    ! arguments
    type(gr_d),intent(in):: grid
    double precision, intent(in):: dt
    double precision, dimension(0:), intent(in):: Kall ! total mixing K
    double precision, dimension(0:), intent(in):: gamma ! nonlocal effect
    double precision, intent(in) :: surface_flux
    double precision, dimension(0:), intent(in):: distrib_flux
    double precision, intent(in) :: bottom_value ! of scalar

    double precision, dimension(1:), intent(out) :: Gvector

    ! local variables
    integer :: index                   ! counter through depth

    ! initialize Gvector to 0
    Gvector = 0.

    ! calculate rhs (the H vector) (D10)

    Gvector(1) = dt / grid%i_space(1) * &
         (Kall(1) * gamma(1) - surface_flux + distrib_flux(1)- distrib_flux(0))

      DO index = 2, grid%M-1
         Gvector(index) = dt / grid%i_space(index) * &
              (Kall(index) * gamma(index) - Kall(index-1) * gamma(index-1) &
               + distrib_flux(index) - distrib_flux(index-1))
      END DO

      Gvector(grid%M) = dt / grid%i_space(grid%M) * &
           (distrib_flux(grid%M) - distrib_flux(grid%M-1) &
            + bottom_value * Kall(grid%M) / grid%g_space(grid%M))      

    END SUBROUTINE diffusion_nonlocal_fluxes

    subroutine diffusion_bot_surf_flux(grid, dt, Kall, surface_flux, &
       bottom_value, Gvector)

! For variables without distributed fluxes and non-local fluxes
! The subroutine calculates the diffusive flux into the bottom of the domain
! and adds the surface flux

    use mean_param, only: gr_d

    implicit none

    ! arguments
    type(gr_d),intent(in):: grid
    double precision, intent(in):: dt
    double precision, dimension(0:), intent(in):: Kall ! total mixing K
    double precision, intent(in) :: surface_flux
    double precision, intent(in) :: bottom_value ! of scalar

    double precision, dimension(:), intent(out) :: Gvector

    ! initialize Gvector to 0
    Gvector = 0.

    ! calculate rhs (the H vector) (D10)

    Gvector(1) = - dt / grid%i_space(1) * surface_flux

    Gvector(grid%M) = dt / grid%i_space(grid%M) * &
           (bottom_value * Kall(grid%M) / grid%g_space(grid%M))      

  END SUBROUTINE diffusion_bot_surf_flux

end module diffusion
