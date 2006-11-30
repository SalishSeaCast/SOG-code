! $Id$
! $Source$

module baroclinic_pressure
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the baroclinic pressure gradient
  ! calculation in the SOG code.
  !
  ! Public Variables:
  !
  !   dPdx_b -- Baroclinic pressure gradient x (cross-strait) component 
  !             [kg/m^2 . s^2]
  !
  !   dPdy_b -- Baroclinic pressure gradient y (along-strait) component 
  !             [kg/m^2 . s^2]
  !
  !   ut -- Integral of u (cross-strait) velocity component over time
  !
  !   vt -- Integral of v (along-strait) velocity component over time
  !
  ! Public Subroutines:
  !
  !   init_baroclinic_pressure -- Initialize baroclinic pressure gradient
  !                               calculation quantities.
  !
  !   baroclinic_P_gradient -- Calculate baroclinic pressure gradient
  !                            components.
  !
  !   new_to_old_vel_integrals -- Copy %new components of the velocity
  !                               component integrals to %old for use in
  !                               next time step.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
       dPdy_b, &  ! Baroclinic pressure gradient y (along-strait) component
       ut,     &  ! Layer expansion (m) at eastern side
       vt,     &  ! layer expansion (m) at western side
       ! Subroutines:
       init_baroclinic_pressure, baroclinic_P_gradient, &
       new_to_old_vel_integrals, dalloc_baro_press_variables

  ! Type Definitions:
  !
  ! Private to module:
  !
  ! New/old array components:
  type :: new_old_arrays
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old_arrays

  ! Parameter Value Declarations:
  !
  ! Private parameter value declarations:
  real(kind=dp) :: &
       Lx = 30.0d3, &  ! Minor axis (cross-strait) of model domain [m]
       Ly = 120.0d3     ! Major axis (along-strait) of model domain [m]

  ! Variable Declarations:
  !
  ! Public:
  real(kind=dp), dimension(:), allocatable :: &
       dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
       dPdy_b     ! Baroclinic pressure gradient y (along-strait) component
  type(new_old_arrays) :: &
       ut, &  ! Layer expansion (m) on eastern side
       vt     ! Layer expansion (m) on northern side
  !
  ! Private variable declarations:
  real(kind=dp), dimension(:), allocatable :: &
       dze, &  ! Isopycnal depth profile array on eastern side
       dzw, &  ! Isopycnal depth profile array on western side
       dzn, &  ! Isopycnal depth profile array on northern side
       dzs     ! Isopycnal depth profile array on southern side

contains

  subroutine init_baroclinic_pressure(M)
    ! Initialize baroclinic pressure gradient calculation quantities.
    implicit none
    ! Argument:
    integer :: M  ! Number of grid points

    ! Allocate memory for baroclinic pressure gradient calculation
    ! quantities
    call alloc_baro_press_variables(M)
    ! Initialize velocity component integrals for baroclinic pressure
    ! gradient calculations
    ut%new = 0.
    vt%new = 0.
  end subroutine init_baroclinic_pressure


  subroutine baroclinic_P_gradient(grid, dt, U_new, V_new, rho_g)
    ! Calculate the baroclinic pressure gradient components.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use grid_mod, only: grid_
    ! Parameter Values:
    use fundamental_constants, only: g
    ! Variable Declarations:
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays; we need
               ! wbar%u(0) & wbar%v(0)
    ! Functions:
    use grid_mod, only: depth_average
    
    implicit none
    
    ! Agruments:
    type(grid_) :: &
         grid  ! Grid depth and spacing arrays
    real(kind=dp), intent(in) :: &
         dt    ! Time step [s]
    real(kind=dp), dimension(0:), intent(inout) :: &
         U_new, &  ! U component (cross-strait, 35 deg) velocity profile [m/s]
         V_new, &  ! V component (along-strait, 305 deg) velocity profile [m/s]
         rho_g     ! Density profile at grid point depths [kg/m^3]

    ! Local variables:
    real(kind=dp) :: &
         oLx, & ! 1/box size for ut
         oLy, & ! 1/box size for vt
         gorLx, & ! g/Lx
         gorLy, & ! g/Ly
         decayscale ! to keep the isopycnals less tilted than size of domain
    integer :: yy ! array variable

    oLx = 2./(Lx/2.) ! 2 infront because extra water
         ! goes into a triangular region, height ut, length 1/2 width of strait
    oLy = 2./(Ly/2.)

    gorLx = g / Lx
    gorLy = g / Ly

    ! Remove barotropic mode from velocity field
                
    U_new = U_new - depth_average(U_new, 0.0d0, grid%D)
    V_new = V_new - depth_average(V_new, 0.0d0, grid%D)

    ! Calculate the distortion of the isopycnals at the grid layer
    ! interfaces that the velocity field produces.
    !
    ! Calculate integrals of u and v over time.  This is from a
    ! geometrical argument that:
    !    u * grid%i_space * dt = ut * grid*i_space * Lx/ 2 / 2
    ! Thus, ut(j) is the amount in m that the u velocity
    ! deepens the isopycnal

    ! if the velocities get large enough the displacement of the
    ! isopycnals reaches the 40 m scale of the model.  To keep the
    ! isopycnals displacements reasonable, say 6 m, introduce a decay which
    ! particularly impacts high flow situations (by being proportional to
    ! the square) but only takes the 6 m displacements down slowly(90 days)

    decayscale = sum(abs(ut%old(1:grid%M)))**2/(90.*86400.)/(6.**2)
    
    ! calculate the added thickness of each layer at the "east" and
    ! "north" side of basin, respectively.  In units of m
    ut%new(1:grid%M) = ut%old(1:grid%M) * (1 - decayscale * dt) &
         + U_new(1:grid%M) * dt * oLx * grid%i_space(1:grid%M)
    vt%new(1:grid%M) = vt%old(1:grid%M) * (1 - decayscale * dt) &
         + V_new(1:grid%M) * dt * oLy * grid%i_space(1:grid%M)

    ! Calculate the depths of the distorted isopycnals
    ! in m
    call iso_distortion (grid%M, grid%i_space, ut%new, dze)
    call iso_distortion (grid%M, grid%i_space, -ut%new, dzw)
    call iso_distortion (grid%M, grid%i_space, vt%new, dzn)
    call iso_distortion (grid%M, grid%i_space, -vt%new, dzs)

    ! Calculate the baroclinic pressure gradient
    
    call delta_p (grid%M, dze, dzw, grid%d_g, rho_g, dPdx_b)
    call delta_p (grid%M, dzn, dzs, grid%d_g, rho_g, dPdy_b)
    
    ! at this point we have calculated delta p, to make it into the 
    ! pressure gradient term we need to multiply by g, divide by rho and
    ! the distance between the sides and the center

    ! also, to keep the flow purely baroclinic, take out the barotropic 
    ! effect of the wind forcing

    dPdx_b = dPdx_b * gorLx / rho_g(1:grid%M) - wbar%u(0) / grid%D
    dPdy_b = dPdy_b * gorLy / rho_g(1:grid%M) - wbar%v(0) / grid%D
    
  end subroutine baroclinic_P_gradient

  subroutine iso_distortion (M, grid_ispace, ut, dz)

    implicit none
    integer :: M
    real(kind=dp), intent(in), dimension(1:) :: grid_ispace, ut
    real(kind=dp), intent(out), dimension(1:) :: dz
    
    real(kind=dp) :: resid_u
    integer :: yy
    
    if (ut(1).gt.-grid_ispace(1)) then
       dz(1) = ut(1) + grid_ispace(1)
       resid_u = 0.
    else
       dz(1) = 0.
       resid_u = ut(1) + grid_ispace(1)
    endif
    do yy = 2,M
       if (ut(yy) + resid_u .gt. -grid_ispace(yy)) then
          dz(yy) = dz(yy-1) + ut(yy) + grid_ispace(yy) + resid_u
          resid_u = 0.
       else
          dz(yy) = dz(yy-1)
          resid_u = ut(yy) + grid_ispace(yy) + resid_u
       endif
    enddo
    
  end subroutine iso_distortion

  subroutine delta_p (M, dzpos, dzneg, dg, rho, dPres)

    implicit none
    integer, intent(in) :: M
    real(kind=dp), intent(in), dimension(1:) :: dzpos, dzneg
    real(kind=dp), intent(in), dimension(0:) :: dg, rho
    real(kind=dp), intent(out), dimension(1:) :: dPres

    real(kind=dp) :: sumpbx, cpos, cneg
    integer :: ii, jj, yy
    
    sumpbx = 0.
    cpos = 0.
    cneg = 0.
    ii = 1
    jj = 1

    ! calculated on grid points
    ! dPdx_b = P(eastern side) - P(western side)
    ! P(eastern side) = sum rho_g * eastern side thickness
    ! but note that this sum must go down to the depth grid%d_g(yy)
    ! and not below

    do yy = 1, M

       ! initialize
       if (yy == 1) then
          dPres(1) = 0.
       else
          dPres(yy) = dPres(yy-1)
       endif

       ! first part, add the pressure on the eastern side
       do while (dzpos(ii) - dg(yy) < epsilon(dzpos(ii) - dg(yy))&
                .and. ii < M)
          dPres(yy) = dPres(yy) + rho(ii) * (dzpos(ii)-cpos)
          cpos = dzpos(ii)
          ii = ii + 1
       enddo
       dPres(yy) = dPres(yy) + rho(ii) * (dg(yy)-cpos)
       cpos = dg(yy)
       
       ! second part, substract the pressure on the western side
       do while (dzneg(jj) - dg(yy) < epsilon(dzneg(jj) - dg(yy))&
            .and. jj < M)
          dPres(yy) = dPres(yy) - rho(jj)*(dzneg(jj)-cneg)
          cneg = dzneg(jj)
          jj = jj + 1
       enddo
       dPres(yy) = dPres(yy) - rho(jj)*(dg(yy)-cneg)
       cneg = dg(yy)
       
       sumpbx = sumpbx + dPres(yy)
    enddo
    
    ! at this point we have pressure on the grid points, assuming that the
    ! pressure at the surface is zero.
    
    ! sumpbx is the average pressure based on that calculation.  If we
    ! assume that the average pressure should be zero, then we need to
    ! subtract this value from the calculated pressure

    sumpbx = sumpbx/float(M)
    
    dPres = dPres - sumpbx
  end subroutine delta_p
  
  subroutine new_to_old_vel_integrals()
    ! Copy %new component of the velocity component integrals to %old
    ! component for use at next time step.

    ut%old = ut%new
    vt%old = vt%new
  end subroutine new_to_old_vel_integrals


  subroutine alloc_baro_press_variables(M)
    ! Allocate memory for baroclinic pressure gradient calculation
    ! quantities.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Baroclinic pressure gradient profile arrays"
    allocate(dPdx_b(1:M), dPdy_b(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Velocity component integral profile arrays"
    allocate(ut%new(0:M+1), ut%old(0:M+1), vt%new(0:M+1), vt%old(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Distorted isopicnal depth profile arrays"
    allocate(dze(1:M), dzw(1:M), dzn(1:M), dzs(1:M),  &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_baro_press_variables


  subroutine dalloc_baro_press_variables
    ! Deallocate memory for baroclinic pressure gradient calculation
    ! quantities.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Baroclinic pressure gradient profile arrays"
    deallocate(dPdx_b, dPdy_b, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Velocity component integral profile arrays"
    deallocate(ut%new, ut%old, vt%new, vt%old, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Distorted isopicnal depth profile arrays"
    deallocate(dze, dzw, dzn,dzs, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_baro_press_variables

end module baroclinic_pressure
