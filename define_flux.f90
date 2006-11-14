! $Id$
! $Source$

module define_flux_mod
  ! *** Temporary module wrapping to allow use of assumed-shape array for
  ! *** arguments.

contains

  subroutine define_flux(grid, U_grad_i, V_grad_i, T_grad_i, S_grad_i, &
       alpha, beta)
    ! *** What's it do?

    use precision_defs, only: dp
    use grid_mod, only: grid_
    use water_properties, only: water_property
    use fundamental_constants, only: g
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays; we need
               ! wbar%u(0) & wbar%v(0)
    USE declarations, only: K, gamma, h

    implicit none

    ! Arguments:
    type(grid_) :: grid  ! Grid arrays
    real(kind=dp), dimension(1:) :: &
         U_grad_i, &  ! Cross-strait velocity gradient profile at interfaces
         V_grad_i, &  ! Along-strait velocity gradient profile at interfaces
         T_grad_i, &  ! Temperature gradient profile at grid interface depths
         S_grad_i     ! Salinity gradient profile at grid interface depths
    type(water_property), intent(in) :: &
         alpha, &  ! Thermal expansion coefficient profile arrays
         beta      ! Saline contraction coefficient profile arrays
    ! Local variable:
    integer :: xx  ! Loop index over depth

    K%u%all = 0.0
    K%s%all = 0.0
    K%t%all = 0.0
    wbar%b_err(0) = 0.

    ! *** Vectorization would probably simplify this
    do xx = 1, h%i
       if (grid%d_i(xx) > h%new) then
          K%u%all(xx) = 0.0
          K%s%all(xx) = 0.0
          K%t%all(xx) = 0.0 
!!$          wbar%t(xx) = 0.
!!$          wbar%s(xx) = 0.
!!$          wbar%v(xx) = 0.
!!$          wbar%u(xx) = 0.
!!$          wbar%b(xx) = 0.
!!$          wbar%b_err(xx) = 0.
       else
!!$          wbar%t(xx) = -K%t%ML(xx) * (T_grad_i(xx) - gamma%t(xx)) ! (9)
!!$          wbar%s(xx) = -K%s%ML(xx) * (S_grad_i(xx) - gamma%s(xx))
!!$          wbar%u(xx) = -K%u%ML(xx) * (U_grad_i(xx) - gamma%m(xx))
!!$          wbar%v(xx) = -K%u%ML(xx) * (V_grad_i(xx) - gamma%m(xx))
!!$          ! Buoyancy flux by definition given t and s flux
!!$          wbar%b(xx) = g * (alpha%i(xx) * wbar%t(xx) - beta%i(xx) * wbar%s(xx))  
!!$          ! Buoyancy flux variation due to error in z
!!$          wbar%b_err(xx) = g * (alpha%grad_i(xx) * wbar%t(xx) - beta%grad_i(xx) * wbar%s(xx)) 
          K%u%all(xx) = K%u%ML(xx)
          K%s%all(xx) = K%s%ML(xx)
          K%t%all(xx) = K%t%ML(xx)
       endif
    enddo

    do xx = h%i, grid%M
       if (grid%d_i(xx) > h%new) then
!!$          wbar%t(xx) = -K%t%total(xx) * T_grad_i(xx)
!!$          wbar%s(xx) = -K%s%total(xx) * S_grad_i(xx)
!!$          wbar%u(xx) = -K%u%total(xx) * U_grad_i(xx)
!!$          wbar%v(xx) = -K%u%total(xx) * V_grad_i(xx)
!!$          wbar%b(xx) = g * (alpha%i(xx) * wbar%t(xx) - beta%i(xx) * wbar%s(xx))
!!$          wbar%b_err(xx) = g * (alpha%grad_i(xx) * wbar%t(xx) - beta%grad_i(xx) * wbar%s(xx))
          K%u%all(xx) = K%u%total(xx)
          K%s%all(xx) = K%s%total(xx)
          K%t%all(xx) = K%t%total(xx)
       endif
    enddo
  end subroutine define_flux

end module define_flux_mod
