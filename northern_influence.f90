module northern_influence
  ! Subroutines that calculate the influence on scalars, including
  ! temperature T and nitrate N, of the underflow of water from the North
  !
  ! Public Subroutines
  !
  ! init_northern() 
  !   - set initial values for the sums
  !
  ! read_northern_params (TO COME)
  !   - Read parameter values from the infile
  !
  ! integrate_northern()
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
       northern_advection

! Type and Variable Declarations:
!
! Private type definitions:
  type :: quantities
     real(kind=dp) :: &
          T,       &  ! Temperature 
          NO,      &  ! Nitrate 
          NH,      &  ! Ammonium
          Si,      &  ! Silicon 
          DIC,     &  ! Dissolved inorganic carbon 
          Oxy         ! Dissolved oxygen column number in CTD
  end type quantities

! Private variables:
  type(quantities) :: sum  ! current integrated value of the various quantities

  real(kind=dp) :: &
        tauN, &         ! integration timescale in s
        strength, &     ! strength of upwelling influence, in 1/m, estimated 
    ! May 18, 2012, page 173, Susan's lab book for initial estimate 1/43m
        central_depth, & ! depth of peak northern influence
        upper_width, &  ! depth range above central depth that is influenced
        lower_width     ! depth range below central depth that is influenced


  public :: quantities, sum  ! just for diagnostics
contains

  subroutine init_northern(T, NO, NH, Si, DIC, Oxy)
    ! initialize sum to initial values
    implicit none
    ! Arguments:
    real(kind=dp) :: &
          T, &   ! current surface temperature
          NO, &  ! current surface nitrate conc.
          NH, &  ! current surface ammonium conc.
          Si, &  ! current surface silicon conc.
          DIC, & ! current surface DIC conc.
          Oxy    ! current surface oxygen conc.

    sum%T = T
    sum%NO = NO
    sum%NH = NH
    sum%Si = Si
    sum%DIC = DIC
    sum%Oxy = Oxy

    call read_northern_params()

  end subroutine init_northern

  subroutine read_northern_params()
    use input_processor, only: getpard
    implicit none

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
    
  end subroutine read_northern_params

  subroutine integrate_northern (T, NO, NH, Si, DIC, Oxy, dt)
    ! updates the integrated value for each of the quantities
    implicit none
    ! Arguments:
    real(kind=dp) :: &
          T, &   ! current surface temperature
          NO, &  ! current surface nitrate conc.
          NH, &  ! current surface ammonium conc.
          Si, &  ! current surface silicon conc.
          DIC, & ! current surface DIC conc.
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
    sum%DIC = b*sum%DIC + a*DIC
    sum%Oxy = b*sum%Oxy + a*Oxy

  end subroutine integrate_northern

  subroutine northern_advection (dt, qty, quantity, advection)
  !  calculate the size of the term for northern advection

    use precision_defs, only: dp
    use grid_mod, only: grid
    use freshwater, only: upwell  ! Maximum freshwater upwelling velocity
    use buoyancy, only: Bnr_profile ! buoyancy profile, calculated here
    use water_properties, only: alpha, rho, Cp ! expansion coefficient, density, specific heat capacity

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

end module northern_influence

