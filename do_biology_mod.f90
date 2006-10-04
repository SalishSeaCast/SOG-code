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

  subroutine do_biology(time, day, dt, M, precision, step_guess, step_min, &
       T_new, I_par, P, NO_new, NH_new, Si_new, Detritus, Gvector)
    ! A wrapper around a bunch of subroutine calls that advance the
    ! biological quantities to the next time step.
    use precision_defs, only: dp
    use mean_param, only: plankton, snow, UVST
    use declarations, only: D_bins, M2   ! need to get rid of these
    use rungekutta, only: odeint
    use biological_mod, only: reaction_p_sog, define_PZ
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
    type(plankton), intent(in) :: P                   ! Plankton
    real(kind=dp), dimension(0:), intent(in) :: &
         NO_new, &  ! Nitrate
         NH_new, &  ! Ammonium
         Si_new     ! Silicon
    type(snow), dimension(D_bins), intent(in) :: Detritus
    type(UVST), intent(in out) :: Gvector

    ! Local variables:
    real(kind=dp), dimension(M2) :: PZ
    integer :: N_ok, N_bad    ! counts bad and good steps in odeint
    real(kind=dp) :: next_time

    ! Load all of the biological quantities into the PZ vector for the
    ! ODE solver to operate on
    call define_PZ(M, P%micro%new, P%nano%new, NO_new, NH_new, & ! in
         Si_new, Detritus,                                     & ! in
         PZ)                                                     ! out
    call check_negative(PZ, 'after define_PZ', time, day)

    ! Solve the biological model for values at the next time step
    next_time = time + dt
    call odeint(PZ, M, M2, time, next_time, precision, step_guess, &
         step_min, &
         N_ok, N_bad, T_new, I_par)
    call check_negative (PZ, 'after odeint', time, day)

    ! Unpack the biological quantities from the PZ vector into the
    ! appropriate components of Gvector
    ! *** This subroutine could have a more meaningful name...
    call reaction_p_sog (M, PZ, P%micro%old, P%nano%old, NO_new, & ! in
         NH_new, Si_new, Detritus,                               & ! in
         Gvector)                                                  ! out

  end subroutine do_biology


  subroutine check_negative (PZ, msg, time, day)
    ! Check for negative values in PZ vector, a fatal error.
    use precision_defs, only: dp
    use io_unit_defs, only: stderr
    implicit none
    ! Arguments:
    real(kind=dp), dimension(:), intent(in) :: PZ
    character(len=*) :: msg
    real(kind=dp), intent(in) :: time
    integer, intent(in) :: day

    if (minval(PZ) < 0.) then
       write(stderr, *) "do_biology: Negative value in PZ ", &
            msg, " day = ", day, " time = ", time
       stop
    endif
  end subroutine check_negative

end module do_biology_mod
