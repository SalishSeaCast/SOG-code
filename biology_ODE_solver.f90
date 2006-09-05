module do_biology_mod

implicit none

private

public :: do_biology

contains

  subroutine do_biology(time, day, dt, M, precision, step_guess, step_min, &
       Temp, I_par, P, N, Detritus, Gvector)

    use precision_defs, only: dp
    use mean_param, only: plankton, snow, nutrient, UVST
    use declarations, only: D_bins, M2   ! need to get rid of these
    use rungekutta, only: odeint
    use biological_mod, only: reaction_p_sog, define_PZ

    implicit none

    real(kind=dp), intent(in) :: time  ! current model time
    integer, intent(in) :: day
    real(kind=dp), intent(in) :: dt  ! time step
    integer, intent(in) :: M ! grid size
    ! for ode solver
    real(kind=dp), intent(in) :: precision, step_guess, step_min
    real(kind=dp), dimension(:), intent(in) :: Temp
    real(kind=dp), dimension(:), intent(in) :: I_par
    type(plankton), intent(in) :: P
    type(nutrient), intent(in) :: N
    type(snow), dimension(D_bins), intent(in) :: Detritus
    type(UVST), intent(in out) :: Gvector
    
    ! local variables
    real(kind=dp), dimension(M2) :: PZ
    integer :: N_ok, N_bad    ! counts bad and good steps in odeint
    real(kind=dp) :: next_time

    call define_PZ(M, &                    !in
          P%micro%new, P%nano%new, N%O%new, N%H%new, Detritus, & !in
          PZ)                                                    !out

    call check_negative (PZ,'PZ after define_PZ',time,day)

    next_time = time+dt

    ! not passing in all of T%new --- can be changed when derivs_sog modulized
    call odeint(PZ, M, M2, time, next_time, precision, &
         step_guess, step_min, &
         N_ok, N_bad, Temp, I_par)
    
    call check_negative (PZ,'PZ after odeint',time,day)

    call reaction_p_sog (M, PZ, & !in  (send in whole PZ, split in subroutine)
          P%micro%old, P%nano%old, N%O%old, &   !in
          N%H%old, Detritus, &                                              !in
          Gvector) 

     Gvector%Sil = 0 ! for now

   end subroutine do_biology

   subroutine check_negative (PZ,string,time,day)

     use precision_defs, only: dp

     implicit none
     
     real(kind=dp), dimension(:), intent(in) :: PZ
     character*(*) :: string
     real(kind=dp), intent(in) :: time
     integer, intent(in) :: day

     IF (MINVAL(PZ) < 0.) THEN
        PRINT "(A)","negative value in"
        PRINT "(A)","string"
        PRINT *,time,day
        STOP
     END IF

   end subroutine check_negative

     

 end module do_biology_mod

