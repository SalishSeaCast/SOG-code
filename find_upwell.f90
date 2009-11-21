! $Id$
! $Source$

module find_upwell
  ! Subroutines that calculate the vertical profile of entrainment
  ! velocity (upwelling), and the vertical advection of quantities (T,
  ! S, etc.).
  !
  ! Public Subroutines:
  !
  ! init_find_upwell()
  !   -Read parameter values from the infile.
  !
  ! upwell_profile(grid, Qriver, upwell, wupwell)
  !   -Calculate the vertical profile of the entrainment velocity based
  !    on the strength of the River, and the maximum upwelling velocity.
  !
  ! vertical_advection(grid, dt, property, wupwell, gvector)
  !   -Calculate the vertical advection of a quantity based on the
  !    vertical profile of the entrainment velocity.
  use precision_defs, only: dp
  implicit none

  private
  public init_find_upwell, upwell_profile, vertical_advection
    
  ! Local variables:
  real(kind=dp) :: d      ! depth in m, of 68% fwc

contains
  subroutine init_find_upwell()
    ! Read parameter values from the infile.
    
    implicit none
  
    ! Read fresh water parameter values from the infile.
    call read_upwell_params()

     

  end subroutine init_find_upwell

  subroutine read_upwell_params()
    ! Read the fresh water depth parameter value from the infile.
    use input_processor, only: getpard
    implicit none
 
    ! Read values for d fit equation from infile
    ! fix d because fit is poor and not enough surface NO3 in winter
    d = getpard("d")
  end subroutine read_upwell_params

  subroutine upwell_profile(Qriver, wupwell)
    ! Calculate the vertical profile of the entrainment velocity based
    ! on the Fraser River flow, and the maximum upwelling velocity.
    ! The latter is a function of the Fraser River flow, and is
    ! calculated in the surface_flux subroutine.  The details of the
    ! model to do this are in the the document entrainment.pdf written
    ! in late June 2006.

    ! type definitions
    use precision_defs, only: dp
    use grid_mod, only: grid
    use freshwater, only: upwell
    use input_processor, only: getpard

    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: Qriver            ! River flow
    ! Vertical upwelling velocity profile
    real(kind=dp), intent(out), dimension(0:) :: wupwell 

    ! Local variables:
    real(kind=dp) :: d25    ! 2.5 d, depth of upwelling variation
    integer :: index        ! counter through depth

    ! Set the depth of 68% fresh water content to the fit found
    ! in entrainment.pdf. Using a value based on the current salinity
    ! profile was less stable. For 
    ! calculating d from the salinity, see code version 1.8
!    d = 11.7 * (Qriver/2720.)**(-0.23)

 
   ! Depth of upwelling variation is defined as 2.5*d (see entrainment.pdf)
    d25 = 2.5*d

    ! Set vertical velocity profile (on interfaces!)
    do index = 0, grid%M
       if (grid%d_g(index) < d25) then
          wupwell(index)= upwell * (1. - ( (1.- grid%d_i(index) / d25)**2) )
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
    use grid_mod, only: grid_

    implicit none

    ! Arguments:
    type(grid_), intent(in) :: grid 
    real(kind=dp), intent(in) :: dt ! time step
    real(kind=dp), dimension (0:), intent(in) :: qty ! quantity to be advected
    real(kind=dp), intent(in), dimension(0:) :: wupwell
    real(kind=dp), intent(inout), dimension(1:) :: Gvector

    ! Local variables:
    integer :: index  ! counter through depth
    real(kind=dp) :: inbottom, outtop, squashing, valuelost

    outtop = 0.
    do index = 1, grid%M
       inbottom = wupwell(index) * qty(index+1) ! upwind scheme
       squashing = wupwell(index) - wupwell(index-1)
       if (-squashing*dt > grid%i_space(index)) then
          write (*,*) "Problem in find_upwell, upwelling too strong for dt"
          stop
       endif
       valuelost = ( qty(index) * & 
            (grid%i_space(index) - wupwell(index-1) * dt) + & 
            qty(index+1) * (wupwell(index) * dt) ) / &
            (grid%i_space(index) + squashing * dt)
       Gvector(index) = Gvector(index) + dt / grid%i_space(index) &
            * (inbottom - outtop - squashing * valuelost)
       outtop = inbottom    ! of previous cell
    enddo
  end subroutine vertical_advection

end module find_upwell
