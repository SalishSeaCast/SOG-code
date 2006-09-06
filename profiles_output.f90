! $Id$
! $Source$

module profiles_output
  ! Variables and subroutines to write profiles and halocline results
  ! to files specified in run parameters file.  The output files
  ! contain headers so that they can be used by the compareSOG package
  ! to create comparison graphs.
  !
  ! Public subroutines:
  !
  ! init_profiles_output -- Get the number of profiles output files
  !                         (different times) to be written, the base file 
  !                         name, and the times at which they are to be 
  !                         written from stdin.  Also get the name of the 
  !                         haloclines output file to be written, open it,
  !                         and write its header.
  !
  ! write_profiles -- Check to see if the time is right to write a
  !                   profiles output file.  If so, open the file,
  !                   write the profile results, and close it.  Also
  !                   write a line of data to the haloclines output
  !                   file.
  !
  ! profiles_output_close -- Close the haloclines output file.

  use precision_defs, only: dp
  use datetime, only: datetime_

  implicit none

  private
  public init_profiles_output, write_profiles, profiles_output_close

  ! Private variable declarations:
  !
  ! Maximum number of profiles (related to storage for date/times)
  integer, parameter :: maxprofiles = 20
  ! Dates/times to write profiles results at
  type(datetime_), dimension(maxprofiles) :: profileDatetime
  ! Temporary storage arrays for profile dates/times (year-days, and
  ! day-seconds)
  !  *** This should go away when we refactor to read
  ! calendar-dates, and clock-times
  integer, dimension(maxprofiles) :: profday
  real(kind=dp), dimension(maxprofiles) :: proftime
  ! Number of profiles output files to be written
  integer :: noprof  
  ! Profiles results files base file name (profile date/time gets appended)
  character(len=80) :: profilesBase_fn

  integer :: iprof ! counter, current profile

contains

  subroutine init_profiles_output(codeId, runDatetime)
    ! Get the number of profiles output files (different times) to be
    ! written, the base file name, and the times at which they are to
    ! be written from stdin.  Also get the name of the haloclines
    ! output file to be written, open it, and write its header.
    use io_unit_defs, only: haloclines, stderr
    implicit none
    ! Arguments:
    character(len=70), intent(in) :: codeId         ! Code identity string
    character(len=19), intent(in) :: runDatetime    ! Date/time of code run
    ! Local variables:
    character(len=80) :: haloclines_fn  ! halocline results file name
    ! Functions from getpar
    ! *** getpar should be a module so we can use these
    integer :: getpari
    character(len=80) :: getpars
    external getpari, getpars

    ! Read the number of profiles results files to be written, and
    ! validate it
    noprof =  getpari("noprof",1)
    if (noprof > maxprofiles) then
       write (stderr, *) "init_profiles_output: Limit on number of ", &
            "profiles is ", maxprofiles, ".  ", noprof, " requested."
       stop
    endif

    if (noprof.gt.0) then
       ! Read the arrays of profile dates/times (year-days, and
       ! day-seconds), and load them into the datetime array
       call getpariv ("profday", profday, noprof, 1)
       call getpardv ("proftime", proftime, noprof, 1)
       profileDatetime%yr_day = profday
       profileDatetime%day_sec = proftime

       ! Read the haloclines results file name, open it, and write its
       ! header
       haloclines_fn = getpars ("haloclinefile", 1)
       open(unit=haloclines, file=haloclines_fn, status="replace", &
            action="write")
       write(haloclines, 100) trim(codeId), runDatetime
100    format("! Halocline depths and magnitudes",/                    &
         "*FromCode: ", a/,                                            &
         "*RunDateTime: ", a/,                                         &
         "*FieldNames: year, month, day, year-day, hour, minute, ",    &
         "second, day-seconds, halocline depth, halocline magnitude"/, &
         "*FieldUnits: yr, mo, d, yrday, hr, min, s, d-s, m, m^-1"/,   &
         "*EndOfHeader")

       ! Read the profiles results file base-name
       profilesBase_fn = getpars("profile_base", 1)
    endif

    iprof = 1

  end subroutine init_profiles_output


  subroutine write_profiles(codeId, runDatetime, year, day, day_time, dt, &
       grid, T, S, rho, Pmicro, Pnano, NO, NH, Sil, dissolved_detritus,   &
       sinking_detritus, lost_Detritus, Ku, Kt, Ks, I_par, U, V)
    ! Check to see if the time is right to write a profiles output
    ! file.  If so, open the file, write the profile results, and
    ! close it.  Also write a line of data to the haloclines output
    ! file.
    use precision_defs, only: dp
    use io_unit_defs, only: profiles, haloclines
    use datetime, only: calendar_date, clock_time, datetime_str
    use unit_conversions, only: KtoC
    use grid_mod, only: grid_
    implicit none
    ! Arguments:
    character(len=70), intent(in) :: codeId         ! Code identity string
    character(len=19), intent(in) :: runDatetime    ! Date/time of code run
    integer, intent(in) :: year, day
    real(kind=dp), intent(in) :: day_time, dt ! can't expect exact time match
    type(grid_), intent(in) :: grid
    real(kind=dp), intent (in) :: &
         T(0:), &       ! Temperature [K]
         S(0:), &       ! Salinity [-]
         rho(0:), &     ! Density [kg/M^3]
         Pmicro(0:), &  ! Micro phytoplankton (diatoms) [uM N]
         Pnano(0:), &   ! Nano phytoplankton (flagellates) [uM N]
         NO(0:), &      ! Nitrates [uM N]
         NH(0:), &      ! Ammonium [uM N]
         Sil(0:), &     ! Silicon [uM]
         dissolved_detritus(0:), &  ! Dissolved detritus [uM N]
         sinking_detritus(0:),   &  ! Dissolved detritus [uM N]
         lost_detritus(0:),      &  ! Dissolved detritus [uM N]
         Ku(0:), &      ! Total momentum eddy diffusivity [m^2/s]
         Kt(0:), &      ! Total temperature eddy diffusivity [m^2/s]
         Ks(0:), &      ! Total salinity eddy diffusivity [m^2/s]
         I_par(0:), &   ! Photosynthetic available radiation [W/m^2]
         U(0:), &       ! U velocity component [m/s]
         V(0:)          ! U velocity component [m/s]
    ! Local variables:
    integer :: i  ! Loop index over grid depth
    real(kind=dp) :: &
         delS, &  ! delta S between two depth points
         dep,  &  ! depth half way between them
         derS     ! -dS/dz
  ! Temporary storage for formated datetime string.  Needed to work around
  ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
  ! calls in write statements
  character(len=19) :: str_proDatetime
  ! sigma-t quantity calculated from density for profile results output
  double precision :: sigma_t

    ! Check the day and the time
    if (iprof < noprof .and. day == profday(iprof)) then
       if (abs(day_time - proftime(iprof)) < 0.5 * dt) then
          ! Calculate the month number and month day for profile headers
          profileDatetime(iprof)%yr = year
          profileDatetime(iprof)%yr_day = day
          profileDatetime(iprof)%day_sec = day_time
          call calendar_date(profileDatetime(iprof))
          call clock_time(profileDatetime(iprof))
          ! Calculate the halocline depth and magnitude
          delS = 0.
          do i =1, grid%M
             if (S(i) - S(i-1) > delS) then ! find the maximum difference
                delS = S(i) - S(i-1)
                dep = 0.5 * (grid%d_g(i) + grid%d_g(i-1)) ! find the depth
                derS = delS / (grid%d_g(i) - grid%d_g(i-1)) ! find the gradient
             endif
          enddo
          ! Write the halocline results
          write(haloclines, 100) profileDatetime(iprof), dep, derS
100       format(i4, 2(1x, i2), 1x, i3, 3(1x, i2), 1x, i5, 2x, f6.2, 2x, f6.3)
          ! Open the profile results file
          open(unit=profiles, &
               file=trim(profilesBase_fn) // '-' &
               // datetime_str(profileDatetime(iprof), 'T'), &
               status='replace', action='write')
          ! Write the profile results file header
          ! Avoid a pgf90 idiocyncracy by getting datetimes formatted into
          ! string here rather than in the write statement
          str_proDatetime = datetime_str(profileDatetime(iprof))
          write(profiles, 200) trim(codeId), runDatetime, str_proDatetime
200       format("! Profiles of Temperature, Salinity, Density, ",           &
               "Phytoplankton (micro & nano),"/,                             &
               "! Nitrate, Ammonium, Silcion and Detritus (dissolved, ",     &
               "sinking, and mortality),"/,                                  &
               "! Total Eddy Diffusivities (momentum, temperature, & ",      &
               "salinity, Photosynthetic"/,                                  &
               "! Available Radiation, and Mean Velocity Components (u & v)"/,&
               "*FromCode: ", a/,                                            &
               "*RunDateTime: ", a/,                                         &
               "*FieldNames: depth, temperature, salinity, sigma-t, ",       &
               "micro phytoplankton, nano phytoplankton, nitrate, ",         &
               "ammonium, silicon, ",                                        &
               "dissolved detritus, ",                                       &
               "sinking detritus, lost detritus, ",                          &
               "total momentum eddy diffusivity, ",                          &
               "total temperature eddy diffusivity, ",                       &
               "total salinity eddy diffusivity, ",                          &
               "photosynthetic available radiation, ",                       &
               "u velocity, v velocity"/,                                    &
               "*FieldUnits: m, deg C, None, None, uM N, uM N, uM N, uM N,", &
               "uM, uM N, uM N, uM N, m^2/s, m^2/s, m^2/s, W/m^2, m/s, m/s"/,&
               "*ProfileDateTime: ", a/,                                     &
               "*EndOfHeader")
          ! Write the profile values at the surface, and at all grid points
          do i = 0, grid%M
             sigma_t = rho(i) - 1000.
             write(profiles, 201) grid%d_g(i), KtoC(T(i)), S(i), sigma_t, &
                  Pmicro(i), Pnano(i), NO(i), NH(i), Sil(i),              &
                  dissolved_detritus(i), sinking_detritus(i),             &
                  lost_detritus(i), Ku(i), Kt(i), Ks(i), I_par(i), U(i), V(i)
          enddo
          ! Write the values at the bottom grid boundary.  Some quantities are
          ! not defined there, so use their values at the Mth grid point.
          sigma_t = rho(grid%M+1) - 1000.
          write(profiles, 201) grid%d_g(grid%M+1), KtoC(T(grid%M+1)),    &
               S(grid%M+1), sigma_t, Pmicro(grid%M+1), Pnano(grid%M+1),  &
               NO(grid%M+1), NH(grid%M+1), Sil(grid%M+1),                &
               dissolved_detritus(grid%M+1), sinking_detritus(grid%M+1), &
               lost_detritus(grid%M+1), Ku(grid%M), Kt(grid%M),          &
               Ks(grid%M), I_par(grid%M), U(grid%M+1), V(grid%M+1)
201       format(f7.3, 80(2x, f8.4))
          close(profiles)          

          iprof = iprof+1

       endif
    endif
  end subroutine write_profiles


  subroutine profiles_output_close
    ! Close the haloclines output file
    use io_unit_defs, only: haloclines
    implicit none

    close(haloclines)
  end subroutine profiles_output_close

end module profiles_output
