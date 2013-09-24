module profiles_output
  ! Variables and subroutines to write profiles and halocline results
  ! to files specified in run parameters file.  The output files
  ! contain headers so that they can be used by the compareSOG package
  ! to create comparison graphs.

  use precision_defs, only: dp
  use datetime, only: datetime_

  implicit none

  private
  public :: &
       ! Variables:
       iprof, noprof, profday, proftime, profileDatetime,   &
       Hoff_startyr, Hoff_startday, Hoff_startsec,          &
       Hoff_endyr, Hoff_endday, Hoff_endsec, Hoff_interval, &
       Hoff_yr, Hoff_day, Hoff_sec,                         &
       ! Subroutines:
       init_profiles_output, &  ! Get the number of profiles output
                                ! files (different times) to be
                                ! written, the base file name, and the
                                ! times at which they are to be
                                ! written from stdin.  Get the name of
                                ! the haloclines output file to be
                                ! written, open it, and write its
                                ! header.  Get the start date/time,
                                ! time interval, and end date/time for
                                ! profile results output to the
                                ! Hoffmueller output file.
       write_std_profiles,   &  ! Check to see if the time is right to
                                ! write a profiles output file.  If
                                ! so, open the file, write the profile
                                ! results, and close it.  Also write a
                                ! line of data to the haloclines
                                ! output file.  Check to see if the
                                ! time is right to add results to the
                                ! Hoffmueller output file, if so,
                                ! write the results.
       profiles_output_close, & ! Close the haloclines output file.
       write_cmn_hdr_id

  ! Private parameter declarations:
  !
  ! Maximum number of profiles (related to storage for date/times)
  integer, parameter :: &
       maxprofiles = 50  ! Maximum number of profiles
  ! (related to storage for date/times)

  ! Private variable declarations:
  !
  integer :: &
       noprof  ! Number of profiles output files to be written
  character(len=80) :: &
       profilesBase_fn, &  ! Profiles results files base file name
       ! (profile date/time gets appended)
       Hoffmueller_fn      ! Hoffmueller results file name
  type(datetime_), dimension(maxprofiles) :: &
       profileDatetime  ! Dates/times to write profiles results at
  ! Temporary storage arrays for profile dates/times (year-days, and
  ! day-seconds)
  !  *** This should go away when we refactor to read
  !  *** calendar-dates, and clock-times
  integer, dimension(maxprofiles) :: &
       profday
  real(kind=dp), dimension(maxprofiles) :: &
       proftime
  integer :: &
       Hoff_startyr,  &  ! Year to start Hoffmueller results output in
       Hoff_startday, &  ! Year-day to start Hoffmueller results output at
       Hoff_startsec, &  ! Day-second to start Hoffmueller results output at
       Hoff_endyr,   &   ! Year to end Hoffmueller results output in
       Hoff_endday,   &  ! Year-day to end Hoffmueller results output at
       Hoff_endsec       ! Day-second to end Hoffmueller results output at
  real(kind=dp) :: &
       Hoff_interval  ! Time interval between Hoffmueller results [day]
  integer :: &
       iprof, &     ! Counter, current profile
       Hoff_yr, &   ! Counter: Next Hoffmueller output year
       Hoff_day, &  ! Counter: Next Hoffmueller output yr-day
       Hoff_sec     ! Counter: Next Hoffmueller output day-sec
  real(kind=dp) :: &
       sumHaloDepth, sumHaloStrength  ! To calculate average values

contains

  subroutine init_profiles_output(str_run_Datetime, CTD_Datetime)
    ! Get the number of profiles output files (different times) to be
    ! written, the base file name, and the times at which they are to
    ! be written from stdin.  Also get the name of the haloclines
    ! output file to be written, open it, and write its header.
    ! Get the start date/time, time interval, and end date/time for
    ! profile results output to the Hoffmueller output file.
    !
    ! The Hoffmueller results file is a series of profile results
    ! written to 1 file so that they can be plotted as a contour or
    ! colourmap plot.
    use io_unit_defs, only: haloclines, Hoffmueller, userHoff, stdout
    use datetime, only: datetime_, datetime_str
    use input_processor, only: getpars, getpari, getpariv, getpardv, getpard
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: &
         str_run_Datetime  ! Date/time of code run
    type(datetime_), intent(in)   :: &
         CTD_Datetime  ! Date/time of CTD profile that initialized the run
    ! Local variables:
    character(len=80) :: &
         haloclines_fn  ! halocline results file name
    ! Temporary storage for formated datetime string.  Needed to work around
    ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
    ! calls in write statements
    character(len=19) :: &
         str_CTD_Datetime  ! CTD profile date/time

    ! Initialize profile counter, and averaging variables
    iprof = 1
    sumHaloDepth = 0
    sumHaloStrength = 0
    ! Convert the initial CTD profile date/time to a string
    str_CTD_Datetime = datetime_str(CTD_Datetime)

    ! Read the number of profiles results files to be written, and
    ! validate it
    noprof =  getpari("noprof")
    if (noprof > maxprofiles) then
       write (stdout, *) "init_profiles_output: Limit on number of ", &
            "profiles is ", maxprofiles, ".  ", noprof, " requested."
       call exit(1)
    endif
    if (noprof > 0) then
       ! Read the arrays of profile dates/times (year-days, and
       ! day-seconds), and load them into the datetime array
       call getpariv("profday", noprof, profday)
       call getpardv("proftime", noprof, proftime)
       profileDatetime%yr_day = profday
       profileDatetime%day_sec = int(proftime)
       ! Read the haloclines results file name, open it, and write its
       ! header
       haloclines_fn = getpars("haloclinefile")
       call system('mkdir -p profiles')
       open(unit=haloclines, file=haloclines_fn, status="replace", &
            action="write")
       call write_haloclines_header(str_run_Datetime, str_CTD_Datetime)
       ! Read the profiles results file base-name
       profilesBase_fn = getpars("profile_base")
    endif

    ! Read the Hoffmueller results file name, its start/end
    ! years/days/times, and its output interval
    Hoffmueller_fn = getpars("Hoffmueller file")
    Hoff_startyr = getpari("Hoffmueller start yr")
    Hoff_startday = getpari("Hoffmueller start day")
    Hoff_startsec = getpari("Hoffmueller start sec")
    Hoff_endyr = getpari("Hoffmueller end yr")
    Hoff_endday = getpari("Hoffmueller end day")
    Hoff_endsec = getpari("Hoffmueller end sec")
    Hoff_interval = getpard("Hoffmueller interval")
    ! Initialize Hoffmueller output year, yr-day, and day-sec counters
    Hoff_yr = Hoff_startyr
    Hoff_day = Hoff_startday
    Hoff_sec = Hoff_startsec
    ! Open the Hoffmueller output results file, and write its header
    open(unit=Hoffmueller, file=Hoffmueller_fn, status="replace", &
         action="write")
    call write_Hoffmueller_header(str_run_Datetime, str_CTD_datetime)
  end subroutine init_profiles_output


  subroutine write_haloclines_header(str_run_Datetime, str_CTD_Datetime)
    ! Write the haloclines results file header
    use io_unit_defs, only: haloclines
    implicit none
    ! Arguments:
    character(len=*), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(haloclines, 100)
100 format("! Halocline depths and magnitudes")
    call write_cmn_hdr_id(haloclines, str_run_Datetime, str_CTD_Datetime)
    write(haloclines, 110)
110 format("*FieldNames: year, month, day, year-day, hour, minute, ",   &
         "second, day-seconds, halocline depth, halocline magnitude, ", &
         "1 m salinity",/,                                              &
         "*FieldUnits: yr, mo, d, yrday, hr, min, s, d-s, m, m^-1, none"/,   &
         "*EndOfHeader")
  end subroutine write_haloclines_header


  subroutine write_Hoffmueller_header(str_run_Datetime, str_CTD_datetime)
    ! Write the Hoffmueller results file header
    use io_unit_defs, only: Hoffmueller
    implicit none
    ! Arguments:
    character(len=*), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(Hoffmueller, 200)
200 format("! Profiles for Hoffmueller diagram")
    call write_cmn_hdr_id(Hoffmueller, str_run_Datetime, str_CTD_Datetime)
    call write_cmn_hdr_fields(Hoffmueller)
    write(Hoffmueller, 210) Hoff_startyr, Hoff_startday, Hoff_startsec, &
         Hoff_endyr, Hoff_endday, Hoff_endsec, Hoff_interval
210 format("*HoffmuellerStartYr: ", i4/,  &
         "*HoffmuellerStartDay: ", i3/,   &
         "*HoffmuellerStartSec: ", i5/,  &
         "*HoffmuellerEndYr: ", i4/,      &
         "*HoffmuellerEndDay: ", i3/,     &
         "*HoffmuellerEndSec: ", i5/,    &
         "*HoffmuellerInterval: ", f6.3/, &
         "*EndOfHeader")
  end subroutine write_Hoffmueller_header


  subroutine write_cmn_hdr_id(unit, str_run_Datetime, str_CTD_Datetime)
    ! Write the common code and run identification header lines.  This
    ! is broken out to reduce code duplication.
    implicit none
    ! Arguments:
    integer :: &
         unit  ! I/O unit to write to
    character(len=*), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(unit, 300) str_run_Datetime, str_CTD_Datetime
300 format("*RunDateTime: ", a/,                                         &
         "*InitialCTDDateTime: ", a)
  end subroutine write_cmn_hdr_id


  subroutine write_std_profiles(str_run_Datetime, str_CTD_Datetime,          &
       year, day, day_time, dt, grid, T, S, rho, Pmicro, Pnano, Ppico, Z,    &
       NO, NH, Si, DIC, Oxy, Alk, D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, &
       Ku, Kt, Ks, U, V)
    ! Check to see if the time is right to write a profiles output
    ! file.  If so, open the file, write the profile results, and
    ! close it.  Also write a line of data to the haloclines output
    ! file.
    use precision_defs, only: dp
    use io_unit_defs, only: profiles, haloclines, Hoffmueller, stdout
    use datetime, only: calendar_date, clock_time, datetime_str
    use unit_conversions, only: KtoC
    use grid_mod, only: grid_
    use irradiance, only: &
         I_par  ! Photosynthetic available radiation profile
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: str_run_Datetime  ! Date/time of code run
    character(len=19), intent(in) :: str_CTD_Datetime  ! Date/time of CTD init
    integer, intent(in) :: year, day
    real(kind=dp), intent(in) :: day_time, dt ! can't expect exact time match
    type(grid_), intent(in) :: grid
    real(kind=dp), intent (in) :: &
         T(0:),       &  ! Temperature [K]
         S(0:),       &  ! Salinity [-]
         rho(0:),     &  ! Density [kg/M^3]
         Pmicro(0:),  &  ! Micro phytoplankton (diatoms) [uM N]
         Pnano(0:),   &  ! Nano phytoplankton (meso-rub) [uM N]
         Ppico(0:),   &  ! Pico phytoplankton (flagellates) [uM N]
         Z(0:),       &  ! Micro zooplankton (uM N)
         NO(0:),      &  ! Nitrates [uM N]
         NH(0:),      &  ! Ammonium [uM N]
         Si(0:),      &  ! Silicon [uM]
         DIC(0:),     &  ! Dissolved inorganic carbon [uM C]
         Oxy(0:),     &  ! Dissolved oxygen [uM O2]
         Alk(0:),     &  ! Alkalinity [uM]
         D_DOC(0:),   &  ! Dissolved organic carbon detritus [uM C]
         D_POC(0:),   &  ! Particulate organic carbon detritus [uM C]
         D_DON(0:),   &  ! Dissolved organic nitrogen detritus [uM N]
         D_PON(0:),   &  ! Particulate organic nitrogen detritus [uM N]
         D_refr(0:),  &  ! Refractory nitrogen detritus [uM N]
         D_bSi(0:),   &  ! Biogenic silicon detritus [uM N]
         Ku(1:),      &  ! Total momentum eddy diffusivity [m^2/s]
         Kt(1:),      &  ! Total temperature eddy diffusivity [m^2/s]
         Ks(1:),      &  ! Total salinity eddy diffusivity [m^2/s]
         U(0:),       &  ! U velocity component [m/s]
         V(0:)           ! U velocity component [m/s]
    ! Local variables:
    integer :: &
         i  ! Loop index over grid depth
    real(kind=dp) :: &
         delS, &  ! delta S between two depth points
         dep,  &  ! depth half way between them
         derS     ! -dS/dz

    ! Check the day and the time
    if (iprof <= noprof .and. day == profday(iprof)) then
       if (abs(day_time - proftime(iprof)) < 0.5d0 * dt) then
          ! Calculate the month number and month day for profile headers
          profileDatetime(iprof)%yr = year
          profileDatetime(iprof)%yr_day = day
          profileDatetime(iprof)%day_sec = int(day_time)
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
          write(haloclines, 400) profileDatetime(iprof), dep, derS, &
               0.5*(S(2)+S(3))
400       format(i4, 2(1x, i2), 1x, i3, 3(1x, i2), 1x, i5, 2x, f6.2, 2x, &
               f6.3, 2x, f6.2)
          ! Add to the averaging variables
          sumHaloDepth = sumHaloDepth + dep
          sumHaloStrength = sumHaloStrength + derS
          ! Open the profile results file
          open(unit=profiles, &
               file=trim(profilesBase_fn) // '-' &
               // datetime_str(profileDatetime(iprof), datetime_sep='T', &
               time_sep='q'), &
               status='replace', action='write')
          ! Write the profile results file header
          call write_profiles_header(str_run_Datetime, str_CTD_Datetime, &
               profileDatetime(iprof), derS, dep)
          ! Write the profiles numbers, and close the profiles file
          call write_profiles_numbers(profiles, grid, T, S, rho,      &
               Pmicro, Pnano, Ppico, Z, NO, NH, Si, DIC, Oxy, Alk,    &
               D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, Ku, Kt, Ks, &
               I_par, U, V)
          close(profiles)
          ! Increment profile counter
          iprof = iprof + 1
       endif
    endif

    ! Check to see if we're in a Hoffmueller results output time step
    if (year == Hoff_yr .and. day == Hoff_day) then
       if ((abs(day_time - Hoff_sec) < 0.5d0 * dt) .or. &
           (abs(day_time - Hoff_sec) < dt .and. Hoff_sec < dt)) then
          ! Write the profiles numbers
          call write_profiles_numbers(Hoffmueller, grid, T, S, rho,   &
               Pmicro, Pnano, Ppico, Z, NO, NH, Si, DIC, Oxy, Alk,    &
               D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, Ku, Kt, Ks, &
               I_par, U, V)
          ! Add empty line as separator
          write(Hoffmueller, *)
          ! Increment the Hoffmueller output counters
          Hoff_sec = Hoff_sec + int(Hoff_interval * 86400.0d0)
          if (Hoff_sec .ge. 86400) then
             Hoff_sec = Hoff_sec - 86400
             Hoff_day = Hoff_day + 1
          endif
          if (Hoff_day > 365) then
             Hoff_yr = Hoff_yr + 1
             Hoff_day = Hoff_day - 365
          endif
          ! Check to see if we're beyond the end of the Hoffmueller
          ! output time interval
          if (Hoff_yr > Hoff_endyr &
               .or. (Hoff_yr == Hoff_endyr &
               .and. Hoff_day > Hoff_endday) &
               .or. (Hoff_yr == Hoff_endyr .and. Hoff_day == Hoff_endday &
               .and. Hoff_sec > Hoff_endsec)) then
             ! Set the Hoffmueller year counter to the year before we
             ! started so that we won't trigger anymore writing
             Hoff_yr = Hoff_startyr - 1
          endif
       endif
    endif
  end subroutine write_std_profiles


  subroutine write_profiles_header(str_run_Datetime, str_CTD_Datetime, &
       pro_Datetime, derS, dep)
    ! Write the profile results file header.
    use precision_defs, only: dp
    use io_unit_defs, only: profiles
    use datetime, only: datetime_, datetime_str
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: &
         str_run_Datetime, &  ! Date/time of code run
         str_CTD_Datetime     ! Date/time of CTD init
    type(datetime_) :: &
         pro_Datetime  ! Profile date/time
    real(kind=dp) :: &
         derS,  &  ! Halocline magnitude
         dep       ! halocline depth
    ! Local variable:
    ! Temporary storage for formated datetime string.  Needed to work around
    ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
    ! calls in write statements
    character(len=19) :: str_pro_Datetime

    str_pro_Datetime = datetime_str(pro_Datetime)
    write(profiles, 400)
400 format("! Profiles of Temperature, Salinity, Density, ",           &
         "Phytoplankton (micro & nano & pico),"/,                      &
         "! Micro zooplankton, Nitrate, Ammonium, Silicon, ",          &
         "Dissolved inorganic carbon, Dissolved Oxygen, Alkalinity,"/, &
         "! Detritus (DOC, POC, DON, PON, refractory N & ",            &
         "biogenic Si), Total Eddy Diffusivities (momentum,"/,         &
         "! temperature & salinity), Photosynthetically ",             &
         "Active Radiation, and"/,                                     &
         "! Mean Velocity Components (u & v)")
    call write_cmn_hdr_id(profiles, str_run_Datetime, str_CTD_Datetime)
    call write_cmn_hdr_fields(profiles)
    write(profiles, 410) str_pro_datetime, derS, dep
410 format("*ProfileDateTime: ", a/,              &
         "*HaloclineMagnitude: ", f6.3, " m^-1"/, &
         "*HaloclineDepth: ", f6.2, " m"/,        &
         "*EndOfHeader")
  end subroutine write_profiles_header


  subroutine write_cmn_hdr_fields(unit)
    ! Write the fieldNames, and FieldUnits part of the profiles
    ! header.  This is broken out to reduce code duplications.
    implicit none
    ! Argument:
    integer :: unit  ! I/O unit to write to

    write(unit, 500)
500 format("*FieldNames: depth, temperature, salinity, sigma-t, ",        &
         "micro phytoplankton, nano phytoplankton, pico phytoplankton, ", &
         "micro zooplankton, nitrate, ammonium, silicon, ",               &
         "dissolved inorganic carbon, dissolved oxygen, alkalinity, ",    &
         "DOC detritus, POC detritus, DON detritus, PON detritus, ",      &
         "refractory N detritus, biogenic Si detritus, ",                 &
         "total momentum eddy diffusivity, ",                             &
         "total temperature eddy diffusivity, ",                          &
         "total salinity eddy diffusivity, ",                             &
         "photosynthetic available radiation, ",                          &
         "u velocity, v velocity"/,                                       &
         "*FieldUnits: m, deg C, None, None, uM N, uM N, uM N, uM N, ",   &
         "uM N, uM N, uM, uM C, uM O2, uM, uM C, uM C, uM N, uM N, ",     &
         "uM N, uM, m^2/s, m^2/s, m^2/s, W/m^2, m/s, m/s")
  end subroutine write_cmn_hdr_fields


  subroutine write_profiles_numbers(unit, grid, T, S, rho, Pmicro,      &
       Pnano, Ppico, Z, NO, NH, Si, DIC, Oxy, Alk, D_DOC, D_POC, D_DON, &
       D_PON, D_refr, D_bSi, Ku, Kt, Ks, I_par, U, V)
    ! Write the profiles numbers.  This is broken out to reduce code
    ! duplications.
    use precision_defs, only: dp
    use io_unit_defs, only: profiles
    use unit_conversions, only: KtoC
    use grid_mod, only: grid_
    implicit none
    ! Arguments:
    integer :: &
         unit  ! I/O unit to write to
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent (in) :: &
         T(0:),       &  ! Temperature [K]
         S(0:),       &  ! Salinity [-]
         rho(0:),     &  ! Density [kg/M^3]
         Pmicro(0:),  &  ! Micro phytoplankton (diatoms) [uM N]
         Pnano(0:),   &  ! Nano phytoplankton (meso-rub) [uM N]
         Ppico(0:),   &  ! Pico phytoplankton (flagellates) [uM N]
         Z(0:),       &  ! Micro zooplankton (uM N)
         NO(0:),      &  ! Nitrates [uM N]
         NH(0:),      &  ! Ammonium [uM N]
         Si(0:),      &  ! Silicon [uM]
         DIC(0:),     &  ! Dissolved inorganic carbon [uM C]
         Oxy(0:),     &  ! Dissolved oxygen [uM O2]
         Alk(0:),     &  ! Alkalinity [uM]
         D_DOC(0:),   &  ! Dissolved organic carbon detritus [uM C]
         D_POC(0:),   &  ! Particulate organic carbon detritus [uM C]
         D_DON(0:),   &  ! Dissolved organic nitrogen detritus [uM N]
         D_PON(0:),   &  ! Particulate organic nitrogen detritus [uM N]
         D_refr(0:),  &  ! Refractory nitrogen detritus [uM N]
         D_bSi(0:),   &  ! Biogenic silicon detritus [uM N]
         Ku(1:),      &  ! Total momentum eddy diffusivity [m^2/s]
         Kt(1:),      &  ! Total temperature eddy diffusivity [m^2/s]
         Ks(1:),      &  ! Total salinity eddy diffusivity [m^2/s]
         I_par(0:),   &  ! Photosynthetic available radiation [W/m^2]
         U(0:),       &  ! U velocity component [m/s]
         V(0:)           ! V velocity component [m/s]
    ! Local variables:
    integer :: i  ! Loop index over grid depth
    real(kind=dp) :: &
         sigma_t  ! sigma-t quantity calculated from density

    ! Write the profile values at the surface.  Eddy diffusivity
    ! arrays don't have values there, so write zeros for them
    sigma_t = rho(0) - 1000.
    write(unit, 600) grid%d_g(0), KtoC(T(0)), S(0), sigma_t, &
         Pmicro(0), Pnano(0), Ppico(0), Z(0), NO(0), NH(0),  &
         Si(0), DIC(0), Oxy(0), Alk(0), D_DOC(0), D_POC(0),  &
         D_DON(0), D_PON(0), D_refr(0), D_bSi(0),            &
         0., 0., 0., I_par(0), U(0), V(0)
    ! Write the profile values at the interior grid points
    do i = 1, grid%M
       sigma_t = rho(i) - 1000.
       write(unit, 600) grid%d_g(i), KtoC(T(i)), S(i), sigma_t,    &
            Pmicro(i), Pnano(i), Ppico(i), Z(i), NO(i), NH(i),     &
            Si(i), DIC(i), Oxy(i), Alk(i), D_DOC(i), D_POC(i),     &
            D_DON(i), D_PON(i), D_refr(i), D_bSi(i), Ku(i), Kt(i), &
            Ks(i), I_par(i), U(i), V(i)
    enddo
    ! Write the values at the bottom grid boundary.  Some quantities are
    ! not defined there, so use their values at the Mth grid point.
    sigma_t = rho(grid%M+1) - 1000.
    write(unit, 600) grid%d_g(grid%M+1), KtoC(T(grid%M+1)),           &
         S(grid%M+1), sigma_t, Pmicro(grid%M+1), Pnano(grid%M+1),     &
         Ppico(grid%M+1), Z(grid%M+1), NO(grid%M+1), NH(grid%M+1),    &
         Si(grid%M+1), DIC(grid%M+1), Oxy(grid%M+1), Alk(grid%M + 1), &
         D_DOC(grid%M+1), D_POC(grid%M+1), D_DON(grid%M+1),           &
         D_PON(grid%M+1), D_refr(grid%M+1), D_bSi(grid%M+1),          &
         Ku(grid%M), Kt(grid%M), Ks(grid%M), I_par(grid%M),           &
         U(grid%M+1), V(grid%M+1)
600 format(f7.3, 10(2x, f8.4), 3(2x, f12.4), 80(2x, f8.4))
  end subroutine write_profiles_numbers


  subroutine profiles_output_close
    ! Close the haloclines and Hoffmueller output files.
    use io_unit_defs, only: haloclines, Hoffmueller, userHoff

    implicit none

    ! write out average values
    write (*,*) iprof
    write (*,*) "For w and mixing tuning"
    write (*,*) "Average halocline depth over profiles", SumHaloDepth/(iprof-1)
    write (*,*) "Average halocline strength over profiles", &
         SumHaloStrength/(iprof-1)

    close(haloclines)
    close(Hoffmueller)
    close(userHoff)
  end subroutine profiles_output_close

end module profiles_output
