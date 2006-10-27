! $Id$
! $Source$

module biology_ODE_solver
  ! Biological model execution module.  A wrapper around a bunch of
  ! subroutine calls that advance the biological quantities to the
  ! next time step.
  !
  ! Public subroutine:
  !
  !   solve_biology_ODEs -- A wrapper around a bunch of subroutine calls that
  !               advance the biological quantities to the next time step.

  implicit none

  private
  public :: solve_biology_ODEs

contains

  subroutine solve_biology_ODEs(time, day, dt, M, precision, step_guess, step_min,    &
       T_new, I_par, Pmicro, Pnano, NO, NH, Si, D_DON, D_PON, D_refr, D_bSi)
    ! Solve the biology model ODEs to advance the biology quantity values
    ! to the next time step, and calculate the growth - mortality terms
    ! (*_RHS%bio) of the semi-implicit diffusion/advection equations.
    use precision_defs, only: dp
    use declarations, only: M2   ! need to get rid of these
    use rungekutta, only: odeint
    use biological_mod, only: PZ, load_PZ
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
         N_ok, N_bad, T_new, I_par, day)
    call check_negative(PZ, 'after odeint', time, day)

    ! Unload the biological quantities from the PZ vector into the
    ! appropriate components of Gvector
    call unload_PZ(M, PZ, Pmicro, Pnano, NO, NH, Si, &
         D_DON, D_PON, D_refr, D_bSi)
  end subroutine solve_biology_ODEs


  subroutine unload_PZ(M, PZ, Pmicro, Pnano, NO, NH, Si, &
       D_DON, D_PON, D_refr, D_bSi)
    ! Unload the biological quantities from the PZ vector into the
    ! appropriate *_RHS%bio arrays.
    use precision_defs, only: dp
    use biological_mod, only: PZ_bins
    use biology_eqn_builder, only: Pmicro_RHS, Pnano_RHS, NO_RHS, NH_RHS, &
         Si_RHS, D_DON_RHS, D_PON_RHS, D_refr_RHS, D_bSi_RHS
    implicit none
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    real(kind=dp), dimension(1:), intent(in) :: &
         PZ  ! Array of biology variables for ODE solver
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton biomass profile array
         Pnano,  &  ! Nano phytoplankton biomass profile array
         NO,     &  ! Nitrate concentration profile array
         NH,     &  ! Ammonium concentration profile array
         Si,     &  ! Silicon concentration profile array
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    ! Local Variables
    integer :: &
         bPZ, &  ! Beginning index for a quantity in the PZ array
         ePZ     ! Ending index for a quantity in the PZ array

    ! Unload micro phytoplankton
    bPZ = (PZ_bins%micro - 1) * M + 1
    ePZ = PZ_bins%micro * M
    Pmicro_RHS%bio = PZ(bPZ:ePZ) - Pmicro(1:M)
    ! Unload nano phytoplankton
    bPZ = (PZ_bins%nano - 1) * M + 1
    ePZ = PZ_bins%nano * M
    Pnano_RHS%bio = PZ(bPZ:ePZ) - Pnano(1:M)
    ! Unload nitrate
    bPZ = (PZ_bins%NO - 1) * M + 1
    ePZ = PZ_bins%NO * M
    NO_RHS%bio = PZ(bPZ:ePZ) - NO(1:M)
    ! Unload ammonimum
    bPZ = (PZ_bins%NH - 1) * M + 1
    ePZ = PZ_bins%NH * M
    NH_RHS%bio = PZ(bPZ:ePZ) - NH(1:M)
    ! Unload silicon
    bPZ = (PZ_bins%Si - 1) * M + 1
    ePZ = PZ_bins%Si * M
    Si_RHS%bio = PZ(bPZ:ePZ) - Si(1:M)
    ! Unload dissolved organic nitrogen detritus
    bPz = (PZ_bins%DON - 1) * M + 1
    ePz = PZ_bins%DON * M
    D_DON_RHS%bio = PZ(bPz:ePz) - D_DON(1:M)
    ! Unload particulate organic nitrogen detritus
    bPz = (PZ_bins%PON - 1) * M + 1
    ePz = PZ_bins%PON * M
    D_PON_RHS%bio = PZ(bPz:ePz) - D_PON(1:M)
    ! Unload refractory nitrogen detritus
    bPz = (PZ_bins%refr - 1) * M + 1
    ePz = PZ_bins%refr * M
    D_refr_RHS%bio = PZ(bPz:ePz) - D_refr(1:M)
    ! Unload biogenic silicon detritus
    bPz = (PZ_bins%bSi - 1) * M + 1
    ePz = PZ_bins%bSi * M
    D_bSi_RHS%bio = PZ(bPz:ePz) - D_bSi(1:M)
  end subroutine unload_PZ


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

end module biology_ODE_solver
