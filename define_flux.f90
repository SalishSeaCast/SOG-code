! $Id$
! $Source$

module define_flux_mod
  ! *** Temporary module wrapping to allow use of assumed-shape array for
  ! *** Cpi argument

contains

  SUBROUTINE define_flux(Cpi)

    use precision_defs, only: dp

    USE mean_param
    USE declarations
    USE surface_forcing

    implicit none

    ! Argument:
    ! Heat capacity of water in J/kg-K at grid interface depths (i.e. Cp%i)
    real(kind=dp), dimension(0:), intent(in) :: Cpi

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
          w%b(xx) = g * (alph%i(xx) * w%t(xx) - beta%i(xx) * w%s(xx))  
          ! Buoyancy flux variation due to error in z
          w%b_err(xx) = g*(alph%idiv(xx)*w%t(xx) -beta%idiv(xx)*w%s(xx)) 
          K%u%all(xx) = K%u%ML(xx)
          K%s%all(xx) = K%s%ML(xx)
          K%t%all(xx) = K%t%ML(xx)
          Q_t(xx) = -w%t(xx) * density%new(xx) * Cpi(xx)
       endif
    enddo

    do xx = h%i, grid%M
       if (grid%d_i(xx) > h%new) then
          w%t(xx) = -K%t%total(xx) * T%div_i(xx)
          w%s(xx) = -K%s%total(xx) * S%div_i(xx)
          w%u(xx) = -K%u%total(xx) * U%div_i(xx)
          w%v(xx) = -K%u%total(xx) * V%div_i(xx)
          w%b(xx) = g * (alph%i(xx) * w%t(xx) - beta%i(xx) * w%s(xx))
          w%b_err(xx) = g * (alph%idiv(xx) * w%t(xx) - beta%idiv(xx) * w%s(xx))
          K%u%all(xx) = K%u%total(xx)
          K%s%all(xx) = K%s%total(xx)
          K%t%all(xx) = K%t%total(xx)
          Q_t(xx) = -w%t(xx) * density%new(xx) * Cpi(xx)
       endif
    enddo
  end subroutine define_flux

end module define_flux_mod
