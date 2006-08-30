module find_upwell

  private
  public upwell_profile, vertical_advection

contains

  SUBROUTINE upwell_profile(grid,S,upwell,wupwell)

    use mean_param, only: gr_d

    IMPLICIT NONE
  
    type(gr_d),INTENT(IN):: grid      ! length of grid 
    double precision, intent(in), dimension(0:) :: S ! current salinity
    double precision, intent(in) :: upwell ! maximum upwelling
    double precision, intent(out), dimension(:) :: &
         wupwell ! vertical profile of upwelling

    ! internal variables

    DOUBLE PRECISION, DIMENSION(1:grid%M)::fwc ! fresh water content
    double precision :: fwc68 ! 68% of total fresh water content
    double precision :: d ! depth in m, of 68% fwc
    double precision :: d25 ! 2.5 d, depth of upwelling variation
    integer :: index ! counter through depth

    ! sum the freshwater content : note fwc is defined on the interfaces
    fwc(1) = (30.0 - S(1)) * grid%i_space(1) ! 30 is the base salinity
    DO index=2,grid%M
       fwc(index)= (30.0 - S(index)) * grid%i_space(index) + fwc(index-1)
    ENDDO

    ! depth of upwelling layer is defined as multiple of d, depth of 68% fwc
    ! see entrainment.pdf

    fwc68 = 0.68*fwc(grid%M)
    index = 1
    do while (fwc(index).lt.fwc68)
       index = index + 1
    end do
    ! depth wanted is between index and index-1
    d = (fwc68 - fwc(index-1)) / (fwc(index) - fwc(index-1)) * &
         grid%i_space(index) + grid%d_i(index-1)

    ! depth of upwelling variation is defined as 2.5*d see entrainment.pdf
    d25 = 2.5*d

    ! set vertical velocity profile
    do index = 1, grid%M+1
       if (grid%d_g(index) < d25) then
          wupwell(index)= upwell * (1. - ( (1.- grid%d_g(index) / d)**2) )
       else
          wupwell(index) = upwell
       endif
    enddo

  end SUBROUTINE upwell_profile

  subroutine vertical_advection (grid, dt, property, wupwell, gvector)

    use mean_param, only: gr_d

    implicit none

    type(gr_d), intent(in) :: grid 
    double precision, intent(in) :: dt ! time step
    double precision, dimension (0:), intent(in) &
         :: property ! quantity to be advected
    double precision, intent(in), dimension(:) :: wupwell
    double precision, intent(in out), dimension(:) :: Gvector

    ! internal variables
    integer :: index  ! counter through depth
    double precision :: inbottom, outtop, squashing, valuelost

    do index = 1, grid%M
       inbottom = wupwell(index+1) * property(index+1) ! upwind scheme
       outtop = wupwell(index) * property(index)       ! upwind scheme
       squashing = wupwell(index+1) - wupwell(index)
       valuelost = property(index) * (grid%i_space(index) - wupwell(index)*dt)&
            + property(index+1) * (wupwell(index+1)*dt)
       Gvector(index) = Gvector(index) + dt / grid%i_space(index) * &
            (inbottom - outtop & 
            - squashing * valuelost / (grid%i_space(index) + squashing * dt))
    enddo

  end subroutine vertical_advection

end module find_upwell
