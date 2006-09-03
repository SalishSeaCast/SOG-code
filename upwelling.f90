! $Id$
! $Source$

module find_upwell
  ! Subroutines that calculate the vertical profile of entrainment
  ! velocity (upwelling), and the vertical advection of quantities (T,
  ! S, etc.).
  !
  ! Public Subroutines:
  !
  ! upwell_profile(grid, S, upwell, wupwell)
  !   -Calculate the vertical profile of the entrainment velocity based
  !    on the salinity profile, and the maximum upwelling velocity.
  !
  ! vertical_advection(grid, dt, property, wupwell, gvector)
  !   -Calculate the vertical advection of a quantity based on the
  !    vertical profile of the entrainment velocity.

  implicit none

  private
  public upwell_profile, vertical_advection

contains

  subroutine upwell_profile(grid, S, upwell, wupwell)
    ! Calculate the vertical profile of the entrainment velocity based
    ! on the salinity profile, and the maximum upwelling velocity.
    ! The latter is a function of the Fraser River flow, and is
    ! calculated in the surface_flux subroutine.  The details of the
    ! model to do this are in the the document entrainment.pdf written
    ! in late June 2006.

    use precision_defs, only: dp
    use mean_param, only: gr_d

    implicit none

    ! Arguments:
    type(gr_d), intent(in) :: grid                 ! Grid properties
    real(kind=dp), intent(in), dimension(0:) :: S  ! Salinity profile
    real(kind=dp), intent(in) :: upwell            ! Maximum upwelling velocity
    ! Vertical upwelling velocity profile
    real(kind=dp), intent(out), dimension(1:) :: wupwell 

    ! Local variables:
    real(kind=dp), dimension(1:grid%M) :: fwc  ! Fresh water content
    real(kind=dp) :: fwc68  ! 68% of total fresh water content
    real(kind=dp) :: d      ! depth in m, of 68% fwc
    real(kind=dp) :: d25    ! 2.5 d, depth of upwelling variation
    integer :: index        ! counter through depth

    ! Sum the freshwater content.  Note that fwc is defined on the interfaces.
    fwc_new = (30.0 - S(1)) * grid%i_space(1) ! 30 is the base salinity
    do index = 2, grid%M
       fwc(index)= (30.0 - S(index)) * grid%i_space(index) + fwc(index-1)
    enddo

    ! Depth of upwelling layer is defined as a multiple of d, the
    ! depth of 68% fwc (see entrainment.pdf)
    fwc68 = 0.68 * fwc(grid%M)
    ! *** This could be moved to a function in grid module that
    ! *** finds the depth at which the specified value of a quantity
    ! *** occurs, if that functionality is needed elsewhere.
    index = 1
    do while (fwc(index) < fwc68)
       index = index + 1
    end do
    ! Depth wanted is between index and index-1
    d = (fwc68 - fwc(index-1)) / (fwc(index) - fwc(index-1)) * &
         grid%i_space(index) + grid%d_i(index-1)

    ! Depth of upwelling variation is defined as 2.5*d (see entrainment.pdf)
    d25 = 2.5*d

    ! Set vertical velocity profile
    do index = 1, grid%M + 1
       if (grid%d_g(index) < d25) then
          wupwell(index)= upwell * (1. - ( (1.- grid%d_g(index) / d)**2) )
       else
          wupwell(index) = upwell
       endif
    enddo
  end subroutine upwell_profile


  subroutine vertical_advection (grid, dt, qty, wupwell, gvector)
    ! Calculate the vertical advection of a quantity based on the
    ! vertical profile of the entrainment velocity.  Note that the
    ! vertical velocity decreases as the surface is approaches and so
    ! water and property is squeezed out the sides of the water column.

    use precision_defs, only: dp
    use mean_param, only: gr_d

    implicit none

    ! Arguments:
    type(gr_d), intent(in) :: grid 
    real(kind=dp), intent(in) :: dt ! time step
    real(kind=dp), dimension (0:), intent(in) :: qty ! quantity to be advected
    real(kind=dp), intent(in), dimension(1:) :: wupwell
    real(kind=dp), intent(inout), dimension(1:) :: Gvector

    ! Local variables:
    integer :: index  ! counter through depth
    real(kind=dp) :: inbottom, outtop, squashing, valuelost

    do index = 1, grid%M
       inbottom = wupwell(index+1) * qty(index+1) ! upwind scheme
       outtop = wupwell(index) * qty(index)       ! upwind scheme
       squashing = wupwell(index+1) - wupwell(index)
       valuelost = qty(index) * (grid%i_space(index) &
            - wupwell(index) * dt) + qty(index+1) * (wupwell(index+1) * dt)
       Gvector(index) = Gvector(index) + dt / grid%i_space(index) &
            * (inbottom - outtop & 
            - squashing * valuelost / (grid%i_space(index) + squashing * dt))
    enddo
  end subroutine vertical_advection

end module find_upwell
