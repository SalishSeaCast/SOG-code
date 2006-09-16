! $Id$
! $Source$

module define_flux_mod
  ! *** Temporary module wrapping to allow use of assumed-shape array for
  ! *** Cp_i argument

contains

  SUBROUTINE define_flux(rho_g, alpha, beta, Cp_i)
    ! *** rho_g should probably be rho_i
    ! *** Would it be more consistent to bring in all of rho & Cp too?

    use precision_defs, only: dp
    use water_properties, only: water_property

    USE mean_param
    USE declarations
    USE surface_forcing

    implicit none

    ! Arguments:
    type(water_property), intent(in) :: &
         alpha, &  ! Thermal expansion coefficient profile arrays
         beta      ! Saline contraction coefficient profile arrays
    real(kind=dp), dimension(0:), intent(in) :: &
         rho_g, &  ! Denisty profile at grid point depths
         Cp_i      ! Specific heat capacity profile at interface depths

    K%u%all = 0.0
    K%s%all = 0.0
    K%t%all = 0.0
    w%b_err(0) = 0.

    ! *** Vectorization would probably simplify this
    do xx = 1, h%i
       if (grid%d_i(xx) > h%new) then
          K%u%all(xx) = 0.0
          K%s%all(xx) = 0.0
          K%t%all(xx) = 0.0 
          w%t(xx) = 0.
          w%s(xx) = 0.
          w%v(xx) = 0.
          w%u(xx) = 0.
          w%b(xx) = 0.
          w%b_err(xx) = 0.
          Q_t(xx) = 0.
       else
          w%t(xx) = -K%t%ML(xx) * (T%div_i(xx) - gamma%t(xx)) ! (9)
          w%s(xx) = -K%s%ML(xx) * (S%div_i(xx) - gamma%s(xx))
          w%u(xx) = -K%u%ML(xx) * (U%div_i(xx) - gamma%m(xx))
          w%v(xx) = -K%u%ML(xx) * (V%div_i(xx) - gamma%m(xx))
          ! Buoyancy flux by definition given t and s flux
          w%b(xx) = g * (alpha%i(xx) * w%t(xx) - beta%i(xx) * w%s(xx))  
          ! Buoyancy flux variation due to error in z
          w%b_err(xx) = g * (alpha%grad_i(xx) * w%t(xx) - beta%grad_i(xx) * w%s(xx)) 
          K%u%all(xx) = K%u%ML(xx)
          K%s%all(xx) = K%s%ML(xx)
          K%t%all(xx) = K%t%ML(xx)
          Q_t(xx) = -w%t(xx) * rho_g(xx) * Cp_i(xx)
       endif
    enddo

    do xx = h%i, grid%M
       if (grid%d_i(xx) > h%new) then
          w%t(xx) = -K%t%total(xx) * T%div_i(xx)
          w%s(xx) = -K%s%total(xx) * S%div_i(xx)
          w%u(xx) = -K%u%total(xx) * U%div_i(xx)
          w%v(xx) = -K%u%total(xx) * V%div_i(xx)
          w%b(xx) = g * (alpha%i(xx) * w%t(xx) - beta%i(xx) * w%s(xx))
          w%b_err(xx) = g * (alpha%grad_i(xx) * w%t(xx) - beta%grad_i(xx) * w%s(xx))
          K%u%all(xx) = K%u%total(xx)
          K%s%all(xx) = K%s%total(xx)
          K%t%all(xx) = K%t%total(xx)
          Q_t(xx) = -w%t(xx) * rho_g(xx) * Cp_i(xx)
       endif
    enddo
  end subroutine define_flux

end module define_flux_mod
