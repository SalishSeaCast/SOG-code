module northern_influence
  ! Subroutines that calculate the influence on scalars, including
  ! temperature T and nitrate N, of the underflow of water from the North
  !
  ! Public Subroutines
  !
  ! init_northern() 
  !   - set initial values for the sums
  !
  ! read_northern_params()
  !   - Read parameter values from the infile
  !
  ! integrate_northern(T, NO, NH, Si, DON, DIC, DOC, Oxy, dt)
  !   - Calculates the back average of the surface values of the parameters
  !   Does this by damping the average with a timescale, tauN, and adding the
  !   new value at an appropriate scale
  !
  ! northern_advection(dt, qty, advection)
  !   - Calculate the influence of the northern advection on a quantity

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! subroutines:
       init_northern, &
       integrate_northern, &
       northern_advection, &
       dalloc_northern_influence_variables, &
       ! variables:
       Northern_return, &  ! include return flow from North?
       Bnr_profile         ! Buoyancy profile of northern return flow

! Type and Variable Declarations:
  !
  ! Public:
  logical :: &
       Northern_return  ! include return flow from North?
    real(kind=dp), dimension(:), allocatable :: &
       Bnr_profile         ! Buoyancy profile of northern return flow
!
! Private type definitions:
  type :: quantities
     real(kind=dp) :: &
          T,       &  ! Temperature
          NO,      &  ! Nitrate
          NH,      &  ! Ammonium
          Si,      &  ! Silicon
          DON,     &  ! Dissolved organic nitrogen 
          DIC,     &  ! Dissolved inorganic carbon
          DOC,     &  ! Dissolved organic carbon 
          Oxy         ! Dissolved oxygen column number in CTD
  end type quantities

! Private variables:
  type(quantities) :: sum  ! current integrated value of the various quantities

  real(kind=dp) :: &
        tauN, &          ! integration timescale in s
        strength, &      ! strength of upwelling influence, in 1/m, estimated
                         ! May 18, 2012, page 173, Susan's lab book for initial
                         !estimate 1/43m
        central_depth, & ! depth of peak northern influence
        upper_width, &   ! depth range above central depth that is influenced
        lower_width      ! depth range below central depth that is influenced


  public :: quantities, sum  ! just for diagnostics
contains

  subroutine init_northern(M)
    ! initialize sum to initial values
    use core_variables, only: &
         T,  &  ! Temperature profile arrays
         N, &  ! current surface nitrate conc.
         Si, &  ! current surface silicon conc.
         D, & ! current surface Detritus %DON diss. org. N  %DOC diss. org. C
         DIC, & ! current surface DIC conc.
         Oxy    ! current surface oxygen conc.
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Allocate memory for fresh water quantity arrays
    call alloc_northern_influence_variables(M)
    ! Read fresh water parameter values from the infile.
    call read_northern_params()
    ! Initialize integrated values of various quantities
    sum%T = T%new(0)
    sum%NO = N%O(0)
    sum%NH = N%H(0)
    sum%Si = Si(0)
    sum%DON = D%DON(0)
    sum%DIC = DIC(0)
    sum%DOC = D%DOC(0)
    sum%Oxy = Oxy(0)
  end subroutine init_northern


  subroutine read_northern_params()
    use input_processor, only: getpard, getparl
    implicit none

    ! Include the return flow from the Northern Strait
    Northern_return = getparl('northern_return_flow_on')
    if (Northern_return) then
       ! Strength of Northern Influence (multiple of 1/43m)
       strength = getpard("strength_northern")/43.
       ! Integration time scale for Northern Influence (days)
       tauN = getpard("tau_northern")*86400.
       ! Central depth for Northern Influence (m)
       central_depth = getpard("depth_northern")
       ! Upper width for Northern Influence (m)
       upper_width = getpard("upper_northern")
       ! Lower width for Northern Influence (m)
       lower_width = getpard("lower_northern")
    endif
  end subroutine read_northern_params

  subroutine integrate_northern (T, NO, NH, Si, DON, DIC, DOC, Oxy, dt)
    ! updates the integrated value for each of the quantities
    implicit none
    ! Arguments:
    real(kind=dp) :: &
          T, &   ! current surface temperature
          NO, &  ! current surface nitrate conc.
          NH, &  ! current surface ammonium conc.
          Si, &  ! current surface silicon conc.
          DON, & ! current surface DON conc.
          DIC, & ! current surface DIC conc.
          DOC, & ! current surface DOC conc.
          Oxy, & ! current surface oxygen conc.
          dt     ! time step

    ! Local variables

    real(kind=dp) :: a, b ! parameters, calculated once

    a = dt/tauN
    b = (1-a)

    sum%T = b*sum%T + a*T
    sum%NO = b*sum%NO + a*NO
    sum%NH = b*sum%NH + a*NH
    sum%Si = b*sum%Si + a*Si
    sum%DON = b*sum%DON + a*DON
    sum%DIC = b*sum%DIC + a*DIC
    sum%DOC = b*sum%DOC + a*DOC
    sum%Oxy = b*sum%Oxy + a*Oxy
  end subroutine integrate_northern


  subroutine northern_advection (dt, qty, quantity, advection)
    !  calculate the size of the term for northern advection
    use precision_defs, only: dp
    use grid_mod, only: grid
    use freshwater, only: upwell  ! Maximum freshwater upwelling velocity
    use water_properties, only: alpha, rho, Cp ! expansion coefficient, density,
                                               ! specific heat capacity
    use fundamental_constants, only: Redfield_C
    implicit none
    ! Arguments
    real(kind=dp), intent(in) :: dt ! time step
    real(kind=dp), intent(in), dimension(0:) :: qty ! quantity to be advected
    character(len=3), intent(in) :: quantity ! name of quantity to be advected
    real(kind=dp), intent(inout), dimension(1:) :: &
         advection ! advection component of quantities RHS vector

    ! Local Variables
    integer :: index ! vertical grid step index
    real(kind=dp) :: requiredsum ! value of integrated sum for the quantity
    real(kind=dp) :: increment ! the increase in temperature at a given depth

    ! choose the sum we wish
    if (quantity .eq. 'T  ') then
       requiredsum = sum%T
    elseif (quantity .eq. 'NO ') then
       requiredsum = sum%NO
    elseif (quantity .eq. 'NH ') then
       requiredsum = sum%NH
    elseif (quantity .eq. 'DIC ') then
    !   requiredsum = sum%DIC
       requiredsum = sum%NO * Redfield_C + 1863.9274d0
    elseif (quantity .eq. 'Oxy ') then
       requiredsum = sum%Oxy
    elseif (quantity .eq. 'Si ') then
       requiredsum = sum%Si    
    elseif (quantity .eq. 'DON') then
       requiredsum = sum%DON
    elseif (quantity .eq. 'DOC') then
       requiredsum = sum%DOC
    else
       write (*,*) 'problems in northern_influence'
       stop
    endif

    do index = 1, grid%M
       if (grid%d_g(index) .lt. central_depth) then
       increment =  dt * strength * upwell * &
            (requiredsum - qty(index)) * &
            exp(-(grid%d_g(index)-central_depth)**2/(upper_width**2))
       else
          increment =  dt * strength * upwell * &
            (requiredsum - qty(index)) * &
            exp(-(grid%d_g(index)-central_depth)**2/(lower_width**2))
       endif
       advection(index) = advection(index) + increment
       if (quantity .eq. 'T  ') then
          Bnr_profile(index) = increment * &
               alpha%g(index) / (cp%g(index) * rho%g(index))
       endif
    enddo

  end subroutine northern_advection


  subroutine alloc_northern_influence_variables(M)
    ! Allocate memory for northern_influence variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Northern Advection Buoyancy Profile Array"
    allocate(Bnr_profile(0:M+1), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_northern_influence_variables


  subroutine dalloc_northern_influence_variables
    ! Deallocate memory from northern influence variables arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Northern Advection Buoyancy Profile Array"
    deallocate(Bnr_profile, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_northern_influence_variables

end module northern_influence

