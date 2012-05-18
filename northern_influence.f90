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

  end subroutine init_northern

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
    real(kind=dp), parameter :: &
        tauN = 10*86400.d0   ! integration timescale
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
    real(kind=dp), parameter :: & 
         f = 1/43.d0 ! strength of upwelling influence, in 1/m, estimated 
    ! May 18, 2012, page 173, Susan's lab book

    ! choose the sum we wish
    if (quantity .eq. 'T') then
       requiredsum = sum%T
    elseif (quantity .eq. 'NO') then
       requiredsum = sum%NO
    else
       write (*,*) 'problems in northern_influence'
       stop
    endif

    do index = 1, grid%M
       advection(index) = advection(index) + f * upwell * &
       (requiredsum - qty(index)) * &
       exp(-(grid%d_g(index)-15.d0)**2/(5.d0**2))
    enddo

  end subroutine northern_advection

end module northern_influence

