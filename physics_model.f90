! $Id$
! $Source$

module physics_model
  ! Type declarations, variables, parameters, and subroutines related
  ! to the physics model in the code.
  !
  ! Public Parameters:
  !
  !   g -- Acceleration due to gravity [m/s^2]
  !
  !   f -- Coriolis factor
  !
  !   pi -- Ratio of circumference to diameter of a circle [-]
  ! 
  ! Public Variables:
  !
  !   B -- Water column buoyancy [m/s^2]
  !
  !   dPdx_b -- Baroclinic pressure gradient x (cross-strait) component 
  !             [kg/m^2 . s^2]
  !
  !   dPdy_b -- Baroclinic pressure gradient y (along-strait) component 
  !             [kg/m^2 . s^2]
  !
  ! Public Subroutines:
  !
  !   init_physics -- Initialize physics model.
  !
  !   baroclinic_P_gradient -- Calculate baroclinic pressure gradient
  !                            components.
  !
  !   double_diffusion -- Calculate double diffusion mixing.
  !
  !   new_to_old_physics -- Copy %new components to %old for use in
  !                         next time step.
  !
  !   dalloc_physics_variables -- Deallocate memory used by physics 
  !                               variables.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Parameter values:
       f,        &  ! Coriolis factor
       g,        &  ! Acceleration due to gravity [m/s^2]
       latitude, &  ! Station S3 latitude [deg]
       pi,       &
       ! Variables:
       B,      &  ! Buoyancy profile array
       dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
       dPdy_b, &  ! Baroclinic pressure gradient y (along-strait) component
       ut, vt, &
       ! Subroutines:
       init_physics, double_diffusion, baroclinic_P_gradient, &
       new_to_old_physics, dalloc_physics_variables

  ! Private module type definitions:
  !
  ! Profiles:
  type :: profiles
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Profile of quantity at current time setp
          old     ! Profile of quantity at previous time step
  end type profiles

  ! Public parameter declarations:
  real(kind=dp) :: &
       f  ! Coriolis factor (would be a parameter bit for a pgf90 bug)
  real(kind=dp), parameter :: &
       g = 9.80665, &  ! Acceleration due to gravity [m/s^2]
       latitude = 49. + 7.517 / 60., & ! Station S3 latitude [deg]
       pi = 3.141592653589793

  ! Public variable declarations:
  real(kind=dp), dimension(:), allocatable :: &
       B, &       ! Buoyancy profile array
       dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
       dPdy_b     ! Baroclinic pressure gradient y (along-strait) component

  ! Private parameter value declarations:
  real(kind=dp) :: &
       Lx = 20.0d3, &  ! Semi-minor axis (cross-strait) of model domain [m]
       Ly = 60.0d3     ! Semi-major axis (along-strait) of model domain [m]
  
  ! Private variable declarations:
  !
  ! Velocity component integrals for baroclinic pressure gradient calculation
  type(profiles) :: &
       ut, &  ! Integral of u (cross-strait) velocity component over time
       vt     ! Integral of v (along-strait) velocity component over time
  real(kind=dp), dimension(:), allocatable :: &
       dzx, &
       dzy

contains

  subroutine init_physics(M)
    ! Initialize physics model.
    use water_properties, only: alloc_water_props
    implicit none
    ! Argument:
    integer :: M  ! Number of grid points

    ! Allocate memory for physics model variables
    call alloc_physics_variables(M)
    ! Allocate memory for water property arrays
    call alloc_water_props(M)
    ! Initialize velocity component integrals for baroclinic pressure
    ! gradient calculations
    ut%new = 0.
    vt%new = 0.
    ! Coriolis factor
    ! *** This must be calculated because pgf90 will not accept an
    ! *** intrinsic in parameter statement
    f = 2. * (2. * pi / 86400.) * sin(pi * latitude / 180.)
!!$    f = 1.1d-4
  end subroutine init_physics


  subroutine new_to_old_physics()
    ! Copy %new component of physics variables to %old component for
    ! use at next time step.
    implicit none
    ! Velocity component integrals for baroclinic pressure gradient
    ! calculations
    ut%old = ut%new
    vt%old = vt%new
  end subroutine new_to_old_physics
  

  subroutine double_diffusion(M, T_grad_i, S_grad_i, alpha_i, beta_i, &
       nu_t_dd, nu_s_dd)
    ! Calculate double diffusion mixing.  
    ! See Large, etal (1994), pp373-374, equations 30 to 34.
    use precision_defs, only: dp
    use io_unit_defs, only: stderr
    implicit none
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    real(kind=dp), dimension(1:), intent(in) :: &
         T_grad_i, &  ! Temperature gradient profile at grid interface depths
         S_grad_i     ! Salinity gradient profile at grid interface depths
    real(kind=dp), dimension(0:), intent(in) :: &
         alpha_i, &  ! Thermal expansion coeff profile at interface depths
         beta_i      ! Salinity contraction coeff profile at interface depths
    real(kind=dp), dimension(1:), intent(out) :: &
         nu_t_dd, &  ! Thermal diffusivity for double diffusion mixing
         nu_s_dd     ! Saline diffusivity for double diffusion mixing

    ! Local variables:
    real(kind=dp), dimension(1:M) :: R_rho  ! Double diffusion density ratio
    integer :: k  ! Loop index over depth
    ! Double diffusion model parameters (see Large, etal (1994) pg 374, fig. 4)
    real(kind=dp), parameter :: &
         nu_f = 10.0d-4, &  ! Maximum salt fingering diffusivity [m^2/s]
         R_rho_o = 1.9, &   ! Salt fingering limit (nu_s = 0)
         p_2 = 3, &         ! Salt finger diff power constant
         nu = 1.5d-06       ! Molecular diffusivity (viscosity) [m^2/s]

    ! Calculate double diffusion density ratio (Large, et al, (1994) eq'n 30)
    R_rho = alpha_i(1:) * T_grad_i / (beta_i(1:) * S_grad_i)
      
    do k = 1, M
       ! Determine if there is double diffision, and if so, what type,
       ! salt fingering, or diffusive convection
       if (1.0 < R_rho(k) &
            .and. R_rho(k) < 2.0 &
            .and. alpha_i(k) * T_grad_i(k) > 0 &
            .and. beta_i(k) * S_grad_i(k) > 0.) then  
          ! Salt fingering
          if (1.0 < R_rho(k) .and. R_rho(k) < R_rho_o) then
             ! Large, et al, (1994) eq'n 31a
             nu_s_dd(k) = nu_f &
                  * (1.0 - ((R_rho(k) - 1.0) &
                             / (R_rho_o - 1.0)) ** 2) ** p_2
          else  ! R_rho >= R_rho_o
             ! Large, et al, (1994) eq'n 31b
             nu_s_dd = 0.
          endif
          ! Large, et al, (1994) eq'n 31c
          nu_t_dd(k) = 0.7 * nu_s_dd(k)
       else if (0. < R_rho(k) &
            .and. R_rho(k) < 1.0 &
            .and. alpha_i(k) * T_grad_i(k) < 0. &
            .and. beta_i(k) * S_grad_i(k) < 0.) then   
          ! Diffusive convection
          ! Large, et al, (1994) eq'n 32
          nu_t_dd(k) = nu * 0.909 &
               * exp(4.6 * exp(-0.54 * (1.0 / R_rho(k) - 1.0)))
          ! Large, et al, (1994) eq'n 34
          if (R_rho(k) >= 0.5 .and. R_rho(k) < 1.0) then
             nu_s_dd(k) = nu_t_dd(k) &
                  * (1.85 - 0.85 / R_rho(k)) * R_rho(k) !(34)
          else if (R_rho(k) < 0.5) then
             nu_s_dd(k) = nu_t_dd(k) * 0.15 * R_rho(k) ! (34)
          endif
       else
          ! No double diffusion (i.e. temperature and salinity are
          ! both stabilizing, or water column is convectively
          ! unstable)
          nu_t_dd(k) = 0.0
          nu_s_dd(k) = 0.0
       endif
    enddo
  end subroutine double_diffusion


  subroutine baroclinic_P_gradient(grid, dt, U_new, V_new, rho_g, &
       w_u, w_v, dPdx_b, dPdy_b)
    ! Calculate the baroclinic pressure gradient components.
    use precision_defs, only: dp
    use grid_mod, only: grid_, depth_average
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
         w_u, &  ! wind stress/rho or momentum flux
         w_v
    real(kind=dp), dimension(1:), intent(out) :: &
         dPdx_b, &  ! Baroclinic pressure gradient x (cross-strait) component
         dPdy_b     ! Baroclinic pressure gradient y (along-strait) component
    ! Local variables:
    real(kind=dp) :: sumu, sumv, sumpbx, sumpby, tol, oLx, oLy, gorLx, gorLy, &
         cz, decayscale
    integer :: yy, ii

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
    ! Calculate integrals of u and v over time.  This from a
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

  subroutine alloc_physics_variables(M)
    ! Allocate memory for physics model variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Buoyancy profile array"
    allocate(B(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
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
  end subroutine alloc_physics_variables


  subroutine dalloc_physics_variables
    ! Deallocate memory from physics model variables arrays.
    use malloc, only: dalloc_check
    use water_properties, only: dalloc_water_props
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Buoyancy profile array"
    deallocate(B, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
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
    call dalloc_water_props
  end subroutine dalloc_physics_variables

end module physics_model
