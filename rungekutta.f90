! $Id$
! $Source$

module rungekutta

! this module contains the subroutines that run the adaptive Runge-Kutta for
! the biological model.
! They are essentially the Numerical Recipes code with small variations to 
! allow temperature and light to passed through to the biological derivatives.

implicit none

private

public :: odeint

contains

  subroutine odeint (PZstart, M, M2, start_time, end_time, precision, &
       step_guess, &
       step_min, N_ok, N_bad, Temp, I_par, day)

    use precision_defs, only: dp
    use biological_mod, only: derivs_sog

    implicit none

    real(kind=dp), dimension(:), intent(in out) :: &
         PZstart ! starting and finishing value of PZ
    integer, intent(in) :: M, &    ! size of the grid
                           M2    ! total size of PZ
    real(kind=dp), intent(in) :: start_time, & ! start time of integration
                                 end_time, & ! end time of integration
                                 precision, & ! required precision
                                 step_guess, & ! guess for time step size
                                 step_min  ! minimum allowed step size
    integer, intent(out) :: N_ok, & ! number of successful steps
                           N_bad  ! number of bad steps
    real(kind=dp), dimension(0:), intent(in):: Temp  ! current temp field
    real(kind=dp), dimension(0:), intent(in):: I_par  ! current light field
    integer, intent(in):: day ! current day (for mort strength)

    ! internal parameters (alphabetical)
    integer :: maxstp ! maximum number of steps
    parameter (maxstp=1000)
    real(kind=dp) :: tiny
    parameter (tiny = 1.d-20)  ! a very small number

    ! internal variables (alphabetical)

    real(kind=dp), dimension(1:M2) :: dPZdt ! biological derivatives vector
    integer :: istp ! step counter
    real(kind=dp), dimension(1:M2) :: PZ ! biological conc. vector
    real(kind=dp), dimension(1:M2) :: PZscal ! scale for size of changes
    real(kind=dp) :: step  ! current step size
    real(kind=dp) :: step_done  ! completed step size
    real(kind=dp) :: step_next ! suggested next step size
    real(kind=dp) :: time ! time

    PZ = PZstart
    time = start_time
    ! set step to abs(step_guess) * sign of end_time-start_time
    step = sign(step_guess,end_time-start_time) 
    N_ok = 0
    N_bad = 0

    do istp = 1, maxstp

       call derivs_sog (M, PZ, dPZdt, Temp, I_par, day)

       PZscal(:) = abs(PZ(:)) + abs(step * dPZdt(:)) + tiny

       if ((time+step-end_time)*(time+step-start_time) > 0) &
            step = end_time - time

       call rkqs (PZ, dPZdt, M, M2, time, step, precision, PZscal, &
            step_done, step_next, Temp, I_par, day)

       if (step_done == step) then
          N_ok = N_ok + 1
       else
          N_bad = N_bad + 1
       endif

       if ( (time - end_time)*(end_time - start_time) .ge. 0) then
          PZstart = PZ
          return
       endif

       if ( abs(step_next) .lt. step_min) &
          write (*,*) 'Stepsize smaller than minimum in ODEINT'
       step = step_next

    enddo
    write (*,*) 'Too many steps in ODEINT'

    return
  end subroutine odeint

subroutine rkqs (PZ, dPZdt, M, M2, time, step_try, precision, PZscal, &
     step_done, step_next, Temp, I_par, day)

    use precision_defs, only: dp
    
    implicit none

    real(kind=dp), dimension(:), intent(in out) :: PZ ! biological conc. vector
    real(kind=dp), dimension(:), intent(in) :: &
                                          dPZdt ! biological derivatives vector
    integer, intent(in) :: M, &    ! size of the grid
                           M2    ! total size of PZ
    real(kind=dp), intent(in out) :: time ! time
    real(kind=dp), intent(in) :: step_try, & ! suggested step size
                                 precision ! required precision
    real(kind=dp), dimension(:), intent(in) :: PZscal !scale for size of change
    real(kind=dp), intent(out) :: step_done, & ! size of time step accomplished
                                  step_next  ! size for next time step
    real(kind=dp), dimension(0:), intent(in):: Temp  ! current temp field
    real(kind=dp), dimension(0:), intent(in):: I_par  ! current light field
    integer, intent(in):: day ! current day for mortality


    ! internal parameters
    real(kind=dp) :: errcon  ! definition for convergence error
    parameter (errcon=1.84e-4)
    ! exponents for how fast to grow or shrink the step
    real(kind=dp) :: Pgrow, Pshrnk 
    parameter (Pgrow=-0.2, Pshrnk=-0.25)
    real(kind=dp) :: safety ! drop down an extra 10% just to be sure
    parameter (safety=0.9)

    ! internal variables (alphabetical)

    real(kind=dp) :: errmax ! current maximum error PZ values
    integer:: i ! counter through PZ
    real(kind=dp), dimension(1:M2) :: PZerr ! error (delta) in PZ
    real(kind=dp), dimension(1:M2) :: PZtemp ! last value of PZ calculated  
    real(kind=dp) :: step  ! current step size
    real(kind=dp) :: time_new ! time after an unsuccessful step

    step = step_try

0001 call rkck (PZ, dPZdt, M, M2, step, PZtemp, PZerr, Temp, I_par, day)

    errmax = 0.

    do i = 1, M2
       errmax = max(errmax,abs(PZerr(i)/PZscal(i)))
    enddo

    errmax = errmax/precision

    if (errmax > 1) then
       step = safety * step * (errmax**Pshrnk)
       if (step < 0.1*step) then   !*** this needs to be fixed!
          step = 0.1*step
       endif
       time_new = time + step
       if (time_new == time) write (*,*) 'Stepsize underflow in rkqs'
       goto 1
    else
       if (errmax > errcon) then
          step_next = safety * step * (errmax**Pgrow)
       else
          step_next = 5. * step
       endif
       step_done = step
       time = time + step
       PZ = PZtemp
       return
    endif

  end subroutine rkqs

  subroutine rkck (PZ, dPZdt, M, M2, step, PZout, PZerr, Temp, I_par, day)

    use precision_defs, only: dp
    use biological_mod, only: derivs_sog
    
    implicit none

    real(kind=dp), dimension(:), intent(in) :: PZ ! biological conc. vector
    real(kind=dp), dimension(:), intent(in) :: &
                                   dPZdt ! biological derivatives vector
    integer, intent(in) :: M, &    ! size of the grid
                           M2    ! total size of PZ
    real(kind=dp), intent(in) :: step ! suggested step size
    real(kind=dp), dimension(:), intent(out) :: PZout ! calculated PZ value   
    real(kind=dp), dimension(:), intent(out) :: PZerr ! error (delta) in PZ
    real(kind=dp), dimension(0:), intent(in):: Temp  ! current temp field
    real(kind=dp), dimension(0:), intent(in):: I_par  ! current light field
    integer, intent(in) :: day ! for mortality fields
    
    ! the Runge-Kutta constants
    real(kind=dp) A2,A3,A4,A5,A6,B21,B31,B32,B41,B42,B43,B51, &
     B52,B53,B54,B61,B62,B63,B64,B65,C1,C3,C4,C6,DC1,DC3,DC4,DC5,DC6
    PARAMETER (A2=.2,A3=.3,A4=.6,A5=1.,A6=.875,B21=.2,B31=3./40., &
     B32=9./40.,B41=.3,B42=-.9,B43=1.2,B51=-11./54.,B52=2.5, &
     B53=-70./27.,B54=35./27.,B61=1631./55296.,B62=175./512., &
     B63=575./13824.,B64=44275./110592.,B65=253./4096.,C1=37./378., &
     C3=250./621.,C4=125./594.,C6=512./1771.,DC1=C1-2825./27648., &
     DC3=C3-18575./48384.,DC4=C4-13525./55296.,DC5=-277./14336., &
     DC6=C6-.25)
 
    ! other internal variables
    real(kind=dp), dimension(1:M2) :: AK2, AK3, AK4, AK5, &
                                      AK6  ! intermediate values of dPZdt
    real(kind=dp), dimension(1:M2) :: PZtemp ! intermediate values of PZ  

    PZtemp(:) = PZ(:) + B21 * step * dPZdt(:)

    call derivs_sog(M, PZtemp, AK2, Temp, I_par, day)

    PZtemp(:) = PZ(:) + step * (B31 * dPZdt(:) + B32 * AK2(:))

    call derivs_sog(M, PZtemp, AK3, Temp, I_par, day)

    PZtemp(:) = PZ(:) + step * (B41 * dPZdt(:) + B42 * AK2(:) + B43 * AK3(:))

    call derivs_sog(M, PZtemp, AK4, Temp, I_par, day)

    PZtemp(:) = PZ(:) + step * (B51 * dPZdt(:) + B52 * AK2(:) + B53 * AK3(:) &
         + B54 * AK4(:))

    call derivs_sog(M, PZtemp, AK5, Temp, I_par, day)

    PZtemp(:) = PZ(:) + step * (B61 * dPZdt(:) + B62 * AK2(:) + B63 * AK3(:) &
         + B64 * AK4(:) + B65*AK5(:))

    call derivs_sog(M, PZtemp, AK6, Temp, I_par, day)

    PZout(:) = PZ(:) + step * (C1 * dPZdt(:) + C3 * AK3(:) + C4 * AK4(:) &
         + C6 * AK6(:))

    PZerr(:) = step * (DC1 * dPZdt(:) + DC3 * AK3(:) + DC4 * AK4(:) &
         + DC5 * AK5(:) + DC6 * AK6(:))

    return

  end subroutine rkck

end module rungekutta
