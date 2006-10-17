! $Id$
! $Source$

module do_biology_mod
  ! Biological model execution module.  A wrapper around a bunch of
  ! subroutine calls that advance the biological quantities to the
  ! next time step.
  !
  ! Public subroutine:
  !
  ! do_biology -- A wrapper around a bunch of subroutine calls that
  !               advance the biological quantities to the next time step.

  implicit none

  private
  public :: do_biology

contains

  subroutine do_biology(time, day, dt, M, precision, step_guess, step_min,    &
       T_new, I_par, Pmicro, Pnano, NO, NH, Si, D_DON, D_PON, D_refr, D_bSi,  &
       Gvector)
    ! A wrapper around a bunch of subroutine calls that advance the
    ! biological quantities to the next time step.
    use precision_defs, only: dp
    use mean_param, only: UVST
    use declarations, only: M2   ! need to get rid of these
    use rungekutta, only: odeint
    use biological_mod, only: PZ, load_PZ, unload_PZ
    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: time  ! Current model time
    integer, intent(in) :: day         ! Current model year-day
    real(kind=dp), intent(in) :: dt    ! Time step
    integer, intent(in) :: M           ! Number of grid points
    ! Passed through to ODE solver:
    real(kind=dp), intent(in) :: precision, step_guess, step_min
    real(kind=dp), dimension(0:), intent(in) :: T_new   ! Temperature
    real(kind=dp), dimension(0:), intent(in) :: I_par  ! Photosynth avail rad
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano, &   ! Nano phytoplankton
         NO, &      ! Nitrate
         NH, &      ! Ammonium
         Si, &      ! Silicon
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    type(UVST), intent(in out) :: Gvector

    ! Local variables:
    integer :: N_ok, N_bad    ! counts bad and good steps in odeint
    real(kind=dp) :: next_time

    ! Load all of the biological quantities into the PZ vector for the
    ! ODE solver to operate on
    call load_PZ(M, Pmicro, Pnano, NO, NH, Si, &
         D_DON, D_PON, D_refr, D_bSi, &
         PZ)                                               ! out
    call check_negative(PZ, 'after define_PZ', time, day)

    ! Solve the biological model for values at the next time step
    next_time = time + dt
    call odeint(PZ, M, M2, time, next_time, precision, step_guess, &
         step_min, &
         N_ok, N_bad, T_new, I_par)
    call check_negative(PZ, 'after odeint', time, day)

    ! Unload the biological quantities from the PZ vector into the
    ! appropriate components of Gvector
    call unload_PZ(M, PZ, Pmicro, Pnano, NO, NH, Si, &
         D_DON, D_PON, D_refr, D_bSi, &
         Gvector)                                                ! out
  end subroutine do_biology


  subroutine check_negative(PZ, msg, time, day)
    ! Check for negative values in PZ vector, a fatal error.
    use precision_defs, only: dp
    use io_unit_defs, only: stderr
    implicit none
    ! Arguments:
    real(kind=dp), dimension(:), intent(in) :: PZ
    character(len=*) :: msg
    real(kind=dp), intent(in) :: time
    integer, intent(in) :: day

    integer :: i ! counter

    if (minval(PZ) < 0.) then
       write(stderr, *) "do_biology: Negative value in PZ ", &
            msg, " day = ", day, " time = ", time
       do i = 1, size(PZ)
          if (PZ(i) < 0.) then
             write (stderr, *) "Value at index, ",i
          endif
       enddo
       stop
    endif
  end subroutine check_negative

end module do_biology_mod
