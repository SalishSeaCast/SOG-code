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
       ut,     &  ! Integral of u (cross-strait) velocity component over time
       vt,     &  ! Integral of v (along-strait) velocity component over time
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
       Lx = 20.0d3, &  ! Semi-minor axis (cross-strait) of model domain [m]
       Ly = 60.0d3     ! Semi-major axis (along-strait) of model domain [m]

  ! Variable Declarations:
  !
  ! Public:
  real(kind=dp), dimension(:), allocatable :: &
       dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
       dPdy_b     ! Baroclinic pressure gradient y (along-strait) component
  type(new_old_arrays) :: &
       ut, &  ! Integral of u (cross-strait) velocity component over time
       vt     ! Integral of v (along-strait) velocity component over time
  !
  ! Private variable declarations:
  real(kind=dp), dimension(:), allocatable :: &
       dzx, &  ! x component of distorted isopicnal depth profile array
       dzy     ! y component of distorted isopicnal depth profile array

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


  subroutine baroclinic_P_gradient(grid, dt, U_new, V_new, rho_g, &
       w_u, w_v)
    ! Calculate the baroclinic pressure gradient components.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: dp
    use grid_mod, only: grid_
    ! Parameter Values:
    use fundamental_constants, only: g
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
    real(kind=dp), intent(in) :: &
         w_u, &  ! Wind stress/rho or momentum flux
         w_v
    ! Local variables:
    real(kind=dp) :: &
         sumu, &
         sumv, &
         sumpbx, &
         sumpby, &
         tol, &
         oLx, &
         oLy, &
         gorLx, &
         gorLy, &
         cz, &
         decayscale
    integer :: &
         yy, &
         ii

    ! Remove barotropic mode from velocity field
        sumu = 0.
        sumv = 0.
        DO yy = 1, grid%M   !remove barotropic mode
           sumu = sumu+U_new(yy)
           sumv = sumv+V_new(yy)
        END DO
        sumu = sumu/grid%M
        sumv = sumv/grid%M

        DO yy = 1, grid%M   !remove barotropic mode
           U_new(yy) = U_new(yy)-sumu
           V_new(yy) = V_new(yy)-sumv
        END DO
!!$    U_new = U_new - depth_average(U_new, 0.0d0, grid%D)
!!$    V_new = V_new - depth_average(V_new, 0.0d0, grid%D)
    ! Calculate the distortion of the isopicnals at the grid layer
    ! interfaces that the velocity field produces.
    !
    ! Calculate integrals of u and v over time.  This is from a
    ! geometrical argument that:
    !    u * grid%i_space * dt = ut * grid*i_space * Lx / 2
    ! Thus, ut(j) is the fraction of the jth grid layer that the u velocity
    ! changes the isopicnal by.
    ! *** The decayscale relaxation factor is a hack to deal with the fact that
    ! *** ut & vt grow too large.
        oLx = 2./(20e3)
        oLy = 2./(60e3)
        gorLx = g / (1025.*20e3)
        gorLy = g / (1025.*60e3)

        decayscale = 1./(3*86400.)
        sumu = 0
        sumv = 0

        ! calculate the added thickness of each layer at the "east" and
        ! "north" side of basin, respectively.  In units of layer thickness
        ut%new = ut%old * (1 - decayscale * dt) + U_new * dt * oLx
        vt%new = vt%old * (1 - decayscale * dt) + V_new * dt * oLy 
!!$    ut%new = ut%old * 0.95 + U_new * dt / (Lx / 2.)
!!$    vt%new = vt%old * 0.95 + V_new * dt / (Ly / 2.)

        ! Calculate the depths of the distorted isopycnals
        ! in unit of layer thickness
        dzx(1) = ut%new(1)+1
        dzy(1) = vt%new(1)+1

        do yy=2,grid%M
           dzx(yy) = dzx(yy-1) + (ut%new(yy) + 1)
           dzy(yy) = dzy(yy-1) + (vt%new(yy) + 1)
        enddo

!!$    dpx(1) = (ut%new(1) + 1.) * grid%i_space(1)
!!$    dpy(1) = (vt%new(1) + 1.) * grid%i_space(1)
!!$    dpx(2:) = dpx(1:grid%M-1) + (ut%new(2:) + 1.) * grid%i_space(2:)
!!$    dpy(2:) = dpy(1:grid%M-1) + (vt%new(2:) + 1.) * grid%i_space(2:)
        ! Calculate the baroclinic pressure gradient
        ! *** Should this tolerance be read in as a run parameter?
        tol = 1e-6
        sumpbx = 0.
        cz = 0.
        ii = 1
        ! dPdx_b = P(eastern side) - P(center)
        ! P(center) = sum rho_g * layer depth
        ! P(eastern side) = sum rho_g * eastern side thickness
        ! but note that this second sum must go down to the depth yy and not
        ! below
        do yy = 1,grid%M
           if (yy == 1) then
              dPdx_b(yy) = -rho_g(yy)
           else
              dPdx_b(yy) = dPdx_b(yy-1) - rho_g(yy)
           endif
           do while ((dzx(ii) - yy) < -tol .and. ii < grid%M)
              dPdx_b(yy) = dPdx_b(yy) + rho_g(ii)*(dzx(ii)-cz)
              cz = dzx(ii)
              ii = ii + 1
           enddo
           dPdx_b(yy) = dPdx_b(yy) + rho_g(ii)*(yy-cz)
           sumpbx = sumpbx + dPdx_b(yy)
           cz = yy
        enddo

        sumpby = 0.
        cz = 0.
        ii=1
        do yy=1,grid%M
           if (yy == 1) then
              dPdy_b(yy) = -rho_g(yy)
           else
              dPdy_b(yy) = dPdy_b(yy-1)-rho_g(yy)
           endif
           do while ((dzy(ii)- yy) <-tol .and. ii < grid%M)
              dPdy_b(yy) = dPdy_b(yy) + rho_g(ii)*(dzy(ii)-cz)
              cz = dzy(ii)
              ii = ii + 1
           enddo
           dPdy_b(yy) = dPdy_b(yy) + rho_g(ii)*(yy-cz)
           sumpby = sumpby + dPdy_b(yy)
           cz = yy
        enddo

        sumpbx = sumpbx/grid%M
        sumpby = sumpby/grid%M

        do yy = 1, grid%M
           dPdx_b(yy) = (dPdx_b(yy) - sumpbx) * gorLx * grid%i_space(yy) &
                - w_u / (grid%M * grid%i_space(yy))
           dPdy_b(yy) = (dPdy_b(yy) - sumpby) * gorLy * grid%i_space(yy) &
                - w_v / (grid%M * grid%i_space(yy))
        enddo
    
  end subroutine baroclinic_P_gradient


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
    allocate(dzx(1:M), dzy(1:M), &
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
    deallocate(dzx, dzy, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_baro_press_variables

end module baroclinic_pressure
