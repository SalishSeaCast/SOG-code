! $Id$
! $Source$

module biology_eqn_builder
  ! Type definitions, variable declarations, and subroutines related
  ! to building the semi-implicit diffusion/advection equations for
  ! the biology quantities in the SOG code.
  !
  ! Public Variables:
  !
  !   diff_coeffs_bio -- Tridiagonal matrix of diffusion coefficient values
  !
  !   Pmicro_RHS -- Micro phytoplankton (diatoms) right-hand side arrays
  !
  !   Pnano_RHS -- Nano phytoplankton (flagellates) right-hand side arrays
  !
  !   NO_RHS -- Nitrate concentration right-hand side arrays
  !
  !   NH_RHS -- Ammonium concentration right-hand side arrays
  !
  !   Si_RHS -- Silicon concentration right-hand side arrays
  !
  !   D_DON_RHS -- Dissolved organic nitrogen detritus right-hand side arrays
  !
  !   D_PON_RHS -- Particulate organic nitro detritus right-hand side arrays
  !
  !   D_refr_RHS -- Refractory nitrogen detritus right-hand side arrays
  !
  !   D_bSi_RHS -- Biogenic silicon detritus right-hand side arrays
  !
  ! Public Subroutines:
  !
  !   build_biology_equations -- Build the right-hand side (RHS) arrays
  !                              for the diffusion/advection equations for
  !                              the biology quantities.
  !
  !   new_to_old_bio_RHS -- Copy %new component of the biology *_RHS%diff_adv
  !                         arrays to the %old component for use by the
  !                         IMEX semi-impllicit PDE solver.
  !
  !   alloc_bio_RHS_variables -- Allocate memory for biology RHS arrays
  !
  !   dalloc_bio_RHS_variables -- Deallocate memory for biology RHS arrays

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       diff_coeffs_bio, &  ! Tridiagonal matrix of diffusion coefficient values
       Pmicro_RHS,      &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,       &  ! Nano phytoplankton (flagellates) RHS arrays
       NO_RHS,          &  ! Nitrate concentration RHS arrays
       NH_RHS,          &  ! Ammonium concentration RHS arrays
       Si_RHS,          &  ! Silicon concentration RHS arrays
       D_DON_RHS,       &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,       &  ! Particulate organic nitrogen detritus RHS arrays
       D_refr_RHS,      &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS,       &  ! Biogenic silicon detritus RHS arrays
       ! Subroutines:
       build_biology_equations, new_to_old_bio_RHS, &
       alloc_bio_RHS_variables, dalloc_bio_RHS_variables

  ! Type Definitions:
  !
  ! Private to module:
  !
  ! Tridiagnonal matrix:
  type :: tridiag
     real(kind=dp), dimension(:), allocatable :: &
          sub_diag, &  ! Sub-diagonal vector of a tridiagonal matrix
          diag,     &  ! Diagonal vector of a tridiagonal matrix
          super_diag   ! Super-diagonal vector of a tridiagonal matrix
  end type tridiag
  !
  ! New/old array components:
  type :: new_old
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old
  !
  ! Semi-implicit diffusion/advection equation right-hand side (RHS) arrays
  type :: RHS
     type(new_old) :: &
          diff_adv  ! Diffusion/advection component of RHS
     real(kind=dp), dimension(:), allocatable :: &
          bio, &  ! Biology (growth - mortality) component of RHS
          sink    ! Sinking component of RHS
  end type RHS
  !
  ! Sinking velocities of biology quantities
  type :: sink_vels
     real(kind=dp) :: &
          Pmicro_min, &  ! Minimum sinking velocity of micro phytos [m/s]
          Pmicro_max, &  ! Maximum sinking velocity of micro phytos [m/s]
          D_PON,      &  ! Sinking velocity of PON detritus [m/s]
          D_bSi          ! Sinking velocity of biogenic silicon detritus [m/s]
  end type sink_vels

  ! Variable Declarations:
  !
  ! Public:
  type(tridiag) :: &
       diff_coeffs_bio  ! Tridiagonal matrix of diffusion coefficient values
  type(RHS) :: &
       Pmicro_RHS, &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,  &  ! Nano phytoplankton (flagellates) RHS arrays
       NO_RHS,     &  ! Nitrate concentration RHS arrays
       NH_RHS,     &  ! Ammonium concentration RHS arrays
       Si_RHS,     &  ! Silicon concentration RHS arrays
       D_DON_RHS,  &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,  &  ! Particulate organic nitro detritus RHS arrays
       D_refr_RHS, &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS      ! Biogenic silicon detritus RHS arrays
  !
  ! Private to module:
  !
  ! Sinking velocities:
  type(sink_vels) :: &
       w_sink  ! Sinking velocities of biology quantities [m/s]

contains

  subroutine build_biology_equations(grid, dt, Pmicro, Pnano, NO, NH,  & ! in
       Si, D_DON, D_PON, D_refr, D_bSi, Ft, K_all, wupwell)
    ! Build the right-hand side (RHS) arrays for the
    ! diffusion/advection equations for the biology quantities.
    !
    ! This calculates the values of the diffusion coefficients matrix
    ! (diff_coeff_bio%*), the RHS diffusion/advection term vectors
    ! (*_RHS%diff_adv%new), and the RHS sinking term vectors (*_RHS%sink).
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use diffusion, only: new_diffusion_coeff, diffusion_bot_surf_flux
    use find_upwell, only: vertical_advection
    use declarations, only: micro  ! *** This definitely needs to go away
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent(in) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    real(kind=dp), intent(in) :: &
         Ft  ! Total fresh water flux
    real(kind=dp), dimension(0:), intent(in) :: &
         K_all  ! Total diffusion coefficient array
    real(kind=dp), dimension(1:), intent(in) :: &
         wupwell  ! Profile of vertical upwelling velocity [m/s]
    ! Local variables:
    real(kind=dp) :: &
         aflux, &  ! Surface nutrient flux
         Pmicro_w_sink  !

    ! Calculate the strength of the diffusion coefficients for biology
    ! model quantities.  They diffuse like salinity.
    call new_diffusion_coeff(grid, dt, K_all,            & ! in
         diff_coeffs_bio%sub_diag, diff_coeffs_bio%diag, & ! out
         diff_coeffs_bio%super_diag)                       ! out

    ! Initialize the RHS *%diff_adv%new arrays, and calculate the diffusive
    ! fluxes at the bottom and top of the grid
    aflux = -Ft * (0. - Pmicro(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         Pmicro(grid%M+1),                               &   ! in
         Pmicro_RHS%diff_adv%new)                            ! out
    aflux = -Ft * (0. - Pnano(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         Pnano(grid%M+1),                                &   ! in
         Pnano_RHS%diff_adv%new)                             ! out
    aflux = -Ft * (6.4 - NO(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         NO(grid%M+1),                                   &   ! in
         NO_RHS%diff_adv%new)                                ! out
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         NH(grid%M+1),                                   &   ! in
         NH_RHS%diff_adv%new)                                ! out
    aflux = -Ft * (60.0 - Si(1)) 
    call diffusion_bot_surf_flux(grid, dt, K_all, aflux, &   ! in
         Si(grid%M+1),                                   &   ! in
         Si_RHS%diff_adv%new)                                ! out
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         D_DON(grid%M+1),                                &   ! in
         D_DON_RHS%diff_adv%new)                             ! out
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         D_PON(grid%M+1),                                &   ! in
         D_PON_RHS%diff_adv%new)                             ! out
    D_refr_RHS%diff_adv%new = 0.
    call diffusion_bot_surf_flux(grid, dt, K_all, 0.d0,  &   ! in
         D_bSi(grid%M+1),                                &   ! in
         D_bSi_RHS%diff_adv%new)                             ! out

    ! Add vertical advection due to upwelling
    call vertical_advection(grid, dt, Pmicro, wupwell, &
         Pmicro_RHS%diff_adv%new)
    call vertical_advection(grid, dt, Pnano, wupwell, &
         Pnano_RHS%diff_adv%new)
    call vertical_advection(grid, dt, NO, wupwell, &
         NO_RHS%diff_adv%new)
    call vertical_advection(grid, dt, NH, wupwell, &
         NH_RHS%diff_adv%new)
    call vertical_advection(grid, dt, Si, wupwell, &
         Si_RHS%diff_adv%new)
    call vertical_advection(grid, dt, D_DON, wupwell, &
         D_DON_RHS%diff_adv%new)
    call vertical_advection(grid, dt, D_PON, wupwell, &
         D_PON_RHS%diff_adv%new)
    call vertical_advection(grid, dt, D_bSi, wupwell, &
         D_bSi_RHS%diff_adv%new)

    ! Calculate the sinking term for the quantities that sink
!!$     micro%sink_min = 0.3*1.1574D-05 ! 0.3/m per day
!!$     micro%sink_max = 1.2*1.1574D-05 ! 1.2/m per day
     w_sink%Pmicro_min = 0.3 / 86400.
     w_sink%Pmicro_max = 1.2 / 86400.
!!$     w_sink%D_PON = rate_det%sink(2)
!!$     w_sink%D_bSi = rate_det%sink(3)
     w_sink%D_PON = 1.157d-5
     w_sink%D_bSi = 2.3d-5
!!$     sflux = micro%sink_min * micro%Nlimit + micro%sink_max *(1 - micro%Nlimit)
     Pmicro_w_sink = w_sink%Pmicro_min * micro%Nlimit(1) &
          + w_sink%Pmicro_max * (1 - micro%Nlimit(1))
     call sinking_advection(grid, dt, Pmicro, Pmicro_w_sink, &
          Pmicro_RHS%sink)
     call sinking_advection(grid, dt, D_PON, w_sink%D_PON, &
          D_PON_RHS%sink)
     call sinking_advection(grid, dt, D_bSi, w_sink%D_bSi, &
          D_bSi_RHS%sink)
  end subroutine build_biology_equations


  subroutine sinking_advection(grid, dt, qty, w_sink, RHS)
    ! Calculate the sinking term of the semi-implicit PDEs for the biology
    ! quantities that sink (some classes of plankton, and some of detritus).
    use precision_defs, only: dp
    use grid_mod, only: grid_
    implicit none
    ! Arguments:
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp) :: &
         dt  ! Time step [s]
    real(kind=dp), dimension(0:), intent(in) :: &
         qty  ! Profile array of sinking quantity
    real(kind=dp), intent(in) :: &
         w_sink  ! Sinking velocity [m/s]
    real(kind=dp), dimension(1:), intent(out) :: &
         RHS  ! RHS term vector for semi-implicit diffusion/advection PDE
    ! Local variables:
    real(kind=dp), dimension(0:grid%M+1) :: &
         Ru,     &
         Ru_1,   &
         Ru_a,   &
         Ru_b,   &
         Au_1,   &
         Au_2,   &
         delta_1
    real(kind=dp), dimension(0:grid%M) :: &
         Fu
    real(kind=dp), dimension(-1:grid%M) :: &
         Ru_2,   &
         delta_2 
    real(kind=dp), dimension(0:grid%M+1, 0:3) :: &
         aa  ! for a third order advection scheme
    integer, dimension(0:grid%M+1) :: &
         L  ! order as a function of grid spacing
    INTEGER :: &
         index, &
         jk

    Ru_a = qty
    !Ru_a(0) = 0.

    !Use qty(0) = 0 for upper boundary condition in advection!!!!

    DO index = 0,grid%M+1
       IF (w_sink >= 0.) THEN
          Au_1(index) = w_sink   ! trivial for sink = constant
          Au_2(index) =  0.
       ELSE
          Au_1(index) = 0.
          Au_2(index) = -w_sink     
       END IF
    END DO

!!!coefficients for cubic polynomial (3C)***see Easter 1993, Monthly Weather Review 121: pp 297-304!!!

    aa = 0.

    DO index = 2,grid%M-1

       aa(index,0) = (-Ru_a(index+1) + 26.0*Ru_a(index) - Ru_a(index-1))/24.0
       aa(index,1) = (-5.0*Ru_a(index+2)+34.0*Ru_a(index+1)-34.0*Ru_a(index-1)+5.0*Ru_a(index-2))/48.0
       aa(index,2) = (Ru_a(index+1)-2.0*Ru_a(index)+Ru_a(index-1))/2.0
       aa(index,3) = (Ru_a(index+2)-2.0*Ru_a(index+1)+2.0*Ru_a(index-1)-Ru_a(index-2))/12.0
       L(index) = 3

    END DO



!!!!Boundary Conditions!!!!!!!!!!!!!!!!!!!!!

    L(0) = 0
    aa(0,0) = Ru_a(0)

    L(1) = 2
    aa(1,0) = (-Ru_a(2) + 26.0*Ru_a(1) - Ru_a(0))/24.0
    aa(1,1) = (Ru_a(2) - Ru_a(0))/2.0
    aa(1,2) = (Ru_a(2) - 2.0*Ru_a(1) + Ru_a(0))/2.0

    L(grid%M) = 2
    aa(grid%M,0) = (-Ru_a(grid%M+1) + 26.0*Ru_a(grid%M) - Ru_a(grid%M-1))/24.0
    aa(grid%M,1) = (Ru_a(grid%M+1) - Ru_a(grid%M-1))/2.0
    aa(grid%M,2) = (Ru_a(grid%M+1) - 2.0*Ru_a(grid%M) + Ru_a(grid%M-1))/2.0

    L(grid%M+1) = 0
    aa(grid%M+1,0) = Ru_a(grid%M+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    delta_1 = 0.  
    delta_2 = 0.  
    Ru_1 = 0.     
    Ru_2 = 0.
    Fu = 0.

    do index = 1, grid%M
       delta_1(index) = Au_1(index) * dt / grid%i_space(index)
       delta_2(index-1) = Au_2(index-1) * dt/ grid%i_space(index)
    end do

    do index = 0, grid%M + 1
       do jk = 0, L(index)
          if (jk == 0) then
             Ru_1(index) = Ru_1(index) + aa(index,0)
             Ru_2(index-1) = Ru_2(index-1) + aa(index,0) 
          else
             Ru_1(index) = Ru_1(index) + aa(index,jk) &
                  / (2.)**(jk+1) / (jk+1)             &
                  * (1. - (1. - 2. * (delta_1(index)  &
                  + epsilon(w_sink)))**(jk+1))          &
                  / (delta_1(index) + epsilon(w_sink))
             Ru_2(index-1) = Ru_2(index-1) + aa(index,jk) &
                  / (-2.)**(jk+1) / (jk+1)                &
                  * ((1. - 2. * (delta_2(index-1)         &
                  + epsilon(w_sink)))**(jk+1) - 1.)         &
                  / (delta_2(index-1) + epsilon(w_sink)) 
          end if
       end do
       if (Ru_1(index) < 0.) then
          Ru_1(index) =  0.
       end if
       ! *** Changing index to index-1 in the next 2 statements is a guess
       ! *** It's consistent with the use of Ru_2 above, and it avoids an
       ! *** array bound runtime error
       if (Ru_2(index-1) < 0.) then
          Ru_2(index-1) = 0.
       end if
    end do

    do index = 0, grid%M+1
       Ru_b(index) = delta_1(index) * Ru_1(index) &
            + delta_2(index-1) * Ru_2(index-1)
       Ru(index) = max(epsilon(w_sink), Ru_a(index), Ru_b(index))
    end do

    !Fu(0) = -Ru_2(0)*Au_2(0)*Ru_a(1)/Ru(1)
    !Fu(grid%M) = Ru_1(grid%M)*Au_1(grid%M)*Ru_a(grid%M)/Ru(grid%M)
    do index = 0, grid%M    !1, grid%M-1
       Fu(index) = Ru_1(index) * Au_1(index) * Ru_a(index) / Ru(index) &
            - Ru_2(index) * Au_2(index) * Ru_a(index+1) / Ru(index+1)
    end do




    !TRY!!!
    Fu(0) = 0.
    RHS = 0.


    DO index = 1, grid%M
       RHS(index) = -dt/grid%i_space(index)*(Fu(index)-Fu(index-1))
       !RHS(index) = -dt/grid%i_space(index)/2.0*w_sink*(qty(index+1)-qty(index-1))
       IF (ABS(RHS(index)) < EPSILON(RHS(index))) THEN
          RHS(index) = 0.
       end if
    end do
  end subroutine sinking_advection


  subroutine new_to_old_bio_RHS()
    ! Copy %new component of the biology *_RHS%diff_adv arrays to the
    ! %old component for use by the IMEX semi-impllicit PDE solver.
    implicit none

    Pmicro_RHS%diff_adv%old = Pmicro_RHS%diff_adv%new
    Pnano_RHS%diff_adv%old = Pnano_RHS%diff_adv%new
    NO_rhs%diff_adv%old = NO_RHS%diff_adv%new
    NH_rhs%diff_adv%old = NH_RHS%diff_adv%new
    Si_rhs%diff_adv%old = Si_RHS%diff_adv%new
    D_DON_rhs%diff_adv%old = D_DON_RHS%diff_adv%new
    D_PON_rhs%diff_adv%old = D_PON_RHS%diff_adv%new
    D_refr_rhs%diff_adv%old = D_refr_RHS%diff_adv%new
    D_bSi_rhs%diff_adv%old = D_bSi_RHS%diff_adv%new
  end subroutine new_to_old_bio_RHS


  subroutine alloc_bio_RHS_variables(M)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
    allocate(diff_coeffs_bio%sub_diag(1:M), diff_coeffs_bio%diag(1:M), &
         diff_coeffs_bio%super_diag(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Micro phytoplankton RHS arrays"
    allocate(Pmicro_RHS%diff_adv%new(1:M), Pmicro_RHS%diff_adv%old(1:M), &
         Pmicro_RHS%bio(1:M), Pmicro_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    allocate(Pnano_RHS%diff_adv%new(1:M), Pnano_RHS%diff_adv%old(1:M), &
         Pnano_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    allocate(NO_RHS%diff_adv%new(1:M), NO_RHS%diff_adv%old(1:M), &
         NO_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    allocate(NH_RHS%diff_adv%new(1:M), NH_RHS%diff_adv%old(1:M), &
         NH_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Silicon concentration RHS arrays"
    allocate(Si_RHS%diff_adv%new(1:M), Si_RHS%diff_adv%old(1:M), &
         Si_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    allocate(D_DON_RHS%diff_adv%new(1:M), D_DON_RHS%diff_adv%old(1:M), &
         D_DON_RHS%bio(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    allocate(D_PON_RHS%diff_adv%new(1:M), D_PON_RHS%diff_adv%old(1:M), &
         D_PON_RHS%bio(1:M), D_PON_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    allocate(D_refr_RHS%diff_adv%new(1:M), D_refr_RHS%diff_adv%old(1:M), &
         D_refr_RHS%bio(1:M), D_refr_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    allocate(D_bSi_RHS%diff_adv%new(1:M), D_bSi_RHS%diff_adv%old(1:M), &
         D_bSi_RHS%bio(1:M), D_bSi_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_bio_RHS_variables


  subroutine dalloc_bio_RHS_variables
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Diffusion coefficients tridiagonal matrix arrays"
    deallocate(diff_coeffs_bio%sub_diag, diff_coeffs_bio%diag, &
         diff_coeffs_bio%super_diag, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Micro phytoplankton RHS arrays"
    deallocate(Pmicro_RHS%diff_adv%new, Pmicro_RHS%diff_adv%old, &
         Pmicro_RHS%bio, Pmicro_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    deallocate(Pnano_RHS%diff_adv%new, Pnano_RHS%diff_adv%old, &
         Pnano_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    deallocate(NO_RHS%diff_adv%new, NO_RHS%diff_adv%old, &
         NO_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    deallocate(NH_RHS%diff_adv%new, NH_RHS%diff_adv%old, &
         NH_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Silicon concentration RHS arrays"
    deallocate(Si_RHS%diff_adv%new, Si_RHS%diff_adv%old, &
         Si_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    deallocate(D_DON_RHS%diff_adv%new, D_DON_RHS%diff_adv%old, &
         D_DON_RHS%bio, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    deallocate(D_PON_RHS%diff_adv%new, D_PON_RHS%diff_adv%old, &
         D_PON_RHS%bio, D_PON_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    deallocate(D_refr_RHS%diff_adv%new, D_refr_RHS%diff_adv%old, &
         D_refr_RHS%bio, D_refr_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    deallocate(D_bSi_RHS%diff_adv%new, D_bSi_RHS%diff_adv%old, &
         D_bSi_RHS%bio, D_bSi_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_bio_RHS_variables

end module biology_eqn_builder
