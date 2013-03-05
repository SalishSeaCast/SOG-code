module user_output
  ! Subroutines for user-specified output in the SOG code.
  !
  ! Public Subroutines:
  !
  !   write_user_phys_timeseries_hdr --
  !
  !   write_user_bio_timeseries_hdr --
  !
  !   write_user_chem_timeseries_hdr --
  !
  !   write_user_timeseries --
  !
  !   init_user_profiles --
  !
  !   write_user_profiles --

  implicit none

  private
  public :: &
       ! Subroutines:
       write_user_phys_timeseries_hdr, write_user_bio_timeseries_hdr, &
       write_user_chem_timeseries_hdr, write_user_timeseries, &
       init_user_profiles, write_user_profiles

  ! Private variable declarations:
  !
  character(len=80) :: &
       userprofilesBase_fn, &  ! User profiles results files base file name
       ! (profile date/time gets appended)
       userHoffmueller_fn      ! User Hoffmueller results file name

contains

  subroutine init_user_profiles(str_run_Datetime, CTD_Datetime)
    ! Open the user-defined Hoffmueller output file for writing
    use io_unit_defs, only: userHoff
    use datetime, only: datetime_, datetime_str
    use profiles_output, only: noprof
    use input_processor, only: getpars
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: &
         str_run_Datetime  ! Date/time of code run
    type(datetime_), intent(in)   :: &
         CTD_Datetime  ! Date/time of CTD profile that initialized the run
    ! Local variables:
    ! Temporary storage for formated datetime string.  Needed to work around
    ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
    ! calls in write statements
    character(len=19) :: &
         str_CTD_Datetime  ! CTD profile date/time

    ! Convert the initial CTD profile date/time to a string
    str_CTD_Datetime = datetime_str(CTD_Datetime)

    if (noprof > 0) then
       ! Read the user profiles results file base-name
       userprofilesBase_fn = getpars("user_profile_base")
    endif

    ! Read the user Hoffmueller results file name
    userHoffmueller_fn = getpars("user Hoffmueller file")

    ! Open the user Hoffmueller output results file, and write its header
    open(unit=userHoff, file=userHoffmueller_fn, status="replace", &
         action="write")
    call write_user_Hoff_hdr(str_run_Datetime, str_CTD_datetime)
  end subroutine init_user_profiles


  subroutine write_user_phys_timeseries_hdr(str_run_Datetime, &
       str_CTD_Datetime, str_start_Datetime)
    ! Write user physics model time series results file header.
    !
    ! !!! This is the place to add exploratory, special, debugging, !!!
    ! !!! etc. output.  Please don't commit this file if you only   !!!
    ! !!! make personal changes here.                               !!!
    !
    ! !!! To include add a variable to this output, add a use statement !!!
    ! !!! for it at the top of the module.  (You may also have to add   !!!
    ! !!! it to the public variables list in the module where it is     !!!
    ! !!! declared.                                                     !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must  !!!
    ! !!! be kept in sync with the appropriate write statement in    !!!
    ! !!! write_user_timeseries(), or compareSOG plotting will fail. !!!
    use io_unit_defs, only: user_phys_timeseries
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_start_Datetime, &  ! Midnight of start day as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(user_phys_timeseries, 100) str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime
100 format("! User-defined time series output from physics model."/,      &
         "! Time series of eastern & northern components of layer "/,     &
         "! expansion at surface, and x & y components of baroclinic "/,  &
         "! pressure gradient at surface."/,                              &
         "*RunDateTime: ", a/,                                            &
         "*InitialCTDDateTime: ", a/,                                     &
         "*FieldNames: time, ut, vt, dPdx_b, dPdy_b"/,                    &
         "*FieldUnits: hr since ", a, " LST, m, m, uPa/m, uPa/m"/,        &
         "*EndOfHeader")
  end subroutine write_user_phys_timeseries_hdr


  subroutine write_user_bio_timeseries_hdr(str_run_Datetime, &
       str_CTD_Datetime, str_start_Datetime)
    ! User biology model time series results
    ! !!! This is the place to add exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    ! !!! Please don't commit this file if you only make personal !!!
    ! !!! changes here. !!!
    ! !!! changes here. !!!
    !
    ! !!! To include add a variable to this output, add a use statement !!!
    ! !!! for it at the top of the module.  (You may also have to add it !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    use io_unit_defs, only: user_bio_timeseries
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_start_Datetime, &  ! Midnight of start day as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(user_bio_timeseries, 103) str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime
103 format("! User-defined time series output from biology model"/,    &
         "! Time series of ..."/,                                      &
         "*RunDateTime: ", a/,                                         &
         "*InitialCTDDateTime: ", a/,                                  &
         "*FieldNames: time, integrated primary productivity"/,        &
         !SEA         "Avg (0-3m) nanoplankton biomass"/, &
         !SEA         "*FieldUnits: hr since ", a, " LST, uM N, uM N"/, &
         "*FieldUnits: hr since ", a, " LST, umol N m-2 hr-1"/, &
         "*EndOfHeader")
  end subroutine write_user_bio_timeseries_hdr


  subroutine write_user_chem_timeseries_hdr(str_run_Datetime, &
       str_CTD_Datetime, str_start_Datetime)
    ! User chemistry model time series results
    ! !!! This is the place to add exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    ! !!! Please don't commit this file if you only make personal !!!
    ! !!! changes here. !!!
    ! !!! changes here. !!!
    !
    ! !!! To include add a variable to this output, add a use statement !!!
    ! !!! for it at the top of the module.  (You may also have to add it !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    use io_unit_defs, only: user_chem_timeseries
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_start_Datetime, &  ! Midnight of start day as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(user_chem_timeseries, 103) str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime
103 format("! User-defined time series output from chemistry model"/,    &
         "! Time series of ..."/,                                      &
         "*RunDateTime: ", a/,                                         &
         "*InitialCTDDateTime: ", a/,                                  &
         !SEA         "*FieldNames: time, Avg (0-3m) microplankton biomass, "       &
    "*FieldNames: time"/,       &
         !SEA         "Avg (0-3m) nanoplankton biomass"/, &
         !SEA         "*FieldUnits: hr since ", a, " LST, uM N, uM N"/, &
         "*FieldUnits: hr since ", a, " LST"/, &
         "*EndOfHeader")
  end subroutine write_user_chem_timeseries_hdr


  subroutine write_user_timeseries(time)
    ! Write user-specified results of the current time step to the
    ! time series files.
    !
    ! !!! Please add use statements here to bring variables in from !!!
    ! !!! other module instead of adding them to the argument list. !!!
    ! !!! Doing the latter requires a change to SOG.f90 that should !!!
    ! !!! be taken out prior to commits on that file.               !!!
    !
    ! !!! You can freely change this subroutine without committing  !!!
    ! !!! changes.                                                  !!!
    use precision_defs, only: dp
    use io_unit_defs, only: user_phys_timeseries, user_bio_timeseries, &
         user_chem_timeseries
    use unit_conversions, only: KtoC
    use baroclinic_pressure, only: ut, vt, dPdx_b, dPdy_b
    use NPZD, only: uptake
    use grid_mod, only: trapz
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: &
         time  ! [hr after start midnight]
    ! Local variables:

    ! Write user-specified physics model time series results
    !
    ! !!! This is the place to add exploratory, special, debugging, !!!
    ! !!! etc. output.  Please don't commit this file if you only   !!!
    ! !!! make personal changes here.                               !!!
    !
    ! !!! To include add a variable to this output, add a use statement !!!
    ! !!! for it at the top of the module.  (You may also have to add   !!!
    ! !!! it to the public variables list in the module where it is     !!!
    ! !!! declared.                                                     !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in         !!!
    ! !!! above, or compareSOG plotting will fail.                        !!!
    write(user_phys_timeseries, 100) time, ut%new(1), vt%new(1), &
         dPdx_b(1)*1e6, dPdy_b(1)*1e6
100 format(f10.4, 80(2x, f8.4))

    ! Write user-specified biology model time series results
    !
    ! !!! This is the place to add exploratory, special, debugging, !!!
    ! !!! etc. output.  Please don't commit this file if you only   !!!
    ! !!! make personal changes here.                               !!!
    !
    ! !!! To include add a variable to this output, add a use statement !!!
    ! !!! for it at the top of the module.  (You may also have to add   !!!
    ! !!! it to the public variables list in the module where it is     !!!
    ! !!! declared.                                                     !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    write(user_bio_timeseries, 200) time, &
         trapz((uptake%NO + uptake%NH) * 3.6d6, 0.0d0, 39.0d0)
200 format(f10.4, 80(2x, f10.4))

    ! Write user-specified chemistry model time series results
    !
    ! !!! This is the place to add exploratory, special, debugging, !!!
    ! !!! etc. output.  Please don't commit this file if you only   !!!
    ! !!! make personal changes here.                               !!!
    !
    ! !!! To include add a variable to this output, add a use statement !!!
    ! !!! for it at the top of the module.  (You may also have to add   !!!
    ! !!! it to the public variables list in the module where it is     !!!
    ! !!! declared.                                                     !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    write(user_chem_timeseries, 300) time
300 format(f10.4)
  end subroutine write_user_timeseries


    subroutine write_user_profiles_hdr(str_run_Datetime, str_CTD_Datetime, &
       pro_Datetime)
    ! Write the user profile results file header.
    use precision_defs, only: dp
    use io_unit_defs, only: userprofiles
    use datetime, only: datetime_, datetime_str
    use profiles_output, only: write_cmn_hdr_id
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: &
         str_run_Datetime, &  ! Date/time of code run
         str_CTD_Datetime     ! Date/time of CTD init
    type(datetime_) :: &
         pro_Datetime  ! Profile date/time
    ! Local variable:
    ! Temporary storage for formated datetime string.  Needed to work around
    ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
    ! calls in write statements
    character(len=19) :: str_pro_Datetime

    ! Top of header
    str_pro_Datetime = datetime_str(pro_Datetime)
    write(userprofiles, 400)
400 format("! Profiles of Depth, Primary Productivity, ", &
         "and Nitrate Upwelling")
    call write_cmn_hdr_id(userprofiles, str_run_Datetime, str_CTD_Datetime)
    write(userprofiles, 500)
500 format("*FieldNames: Depth, Primary productivity, Nitrate Upwelling"/, &
         "*FieldUnits: m, uM C/hr, umol N m-2 s-1")
    write(userprofiles, 410) str_pro_datetime
410 format("*ProfileDateTime: ", a/, "*EndOfHeader")
  end subroutine write_user_profiles_hdr


  subroutine write_user_Hoff_hdr(str_run_Datetime, str_CTD_datetime)
    ! Write the user Hoffmueller results file header
    use io_unit_defs, only: userHoff
    use profiles_output, only: &
         Hoff_startyr, Hoff_startday, Hoff_startsec, &
         Hoff_endyr, Hoff_endday, Hoff_endsec, Hoff_interval
    implicit none
    ! Arguments:
    character(len=*), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    ! Top of header
    write(userHoff, 200)
200 format("! Profiles for user-defined Hoffmueller diagram")
    write(userHoff, 300) str_run_Datetime, str_CTD_Datetime
300 format("*RunDateTime: ", a/, "*InitialCTDDateTime: ", a)
    write(userHoff, 500)
500 format("*FieldNames: Depth, Primary productivity, Nitrate Upwelling"/, &
         "*FieldUnits: m, uM C/hr, umol N m-2 s-1")
    write(userHoff, 210) Hoff_startyr, Hoff_startday, Hoff_startsec, &
         Hoff_endyr, Hoff_endday, Hoff_endsec, Hoff_interval
210 format("*HoffmuellerStartYr: ", i4/,  &
         "*HoffmuellerStartDay: ", i3/,   &
         "*HoffmuellerStartSec: :", i5/,  &
         "*HoffmuellerEndYr: ", i4/,      &
         "*HoffmuellerEndDay: ", i3/,     &
         "*HoffmuellerEndSec: :", i5/,    &
         "*HoffmuellerInterval: ", f7.3/, &
         "*EndOfHeader")
  end subroutine write_user_Hoff_hdr


  subroutine write_user_profiles(str_run_Datetime, str_CTD_Datetime, &
       year, day, day_time, dt, grid)
    ! Write user-specified profiles results.
    !
    ! Recommended signature:
    !   subroutine write_user_profiles(str_run_Datetime, &
    !         str_CTD_Datetime, year, day, day_time, dt, grid)
    ! Arguments have been removed from the subroutine signature because
    ! they are not presently used.
    !
    ! !!! Please add use statements here to bring variables in from !!!
    ! !!! other module instead of adding them to the argument list. !!!
    ! !!! Doing the latter requires a change to SOG.f90 that should !!!
    ! !!! be taken out prior to commits on that file.               !!!
    !
    ! !!! You can freely change this subroutine without committing  !!!
    ! !!! changes.
    use precision_defs, only: dp
    use io_unit_defs, only: userprofiles, userHoff
    use datetime, only: calendar_date, clock_time, datetime_str
    use grid_mod, only: grid_
    use datetime, only: datetime_
    use profiles_output, only: &
         iprof, noprof, profday, proftime, profileDatetime, &
         Hoff_yr, Hoff_day, Hoff_sec
    use NPZD, only: uptake
    use fundamental_constants, only: Redfield_C
    use upwelling, only: w_upwell
    use core_variables, only: N
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: str_run_Datetime  ! Date/time of code run
    character(len=19), intent(in) :: str_CTD_Datetime  ! Date/time of CTD init
    integer, intent(in) :: year, day
    real(kind=dp), intent(in) :: day_time, dt ! can't expect exact time match
    type(grid_), intent(in) :: grid

    ! Check the day and the time
    if (iprof <= noprof .and. day == profday(iprof)) then
       if (abs(day_time - proftime(iprof)) < 0.5d0 * dt) then
          ! Calculate the month number and month day for profile headers
          profileDatetime(iprof)%yr = year
          profileDatetime(iprof)%yr_day = day
          profileDatetime(iprof)%day_sec = int(day_time)
          call calendar_date(profileDatetime(iprof))
          call clock_time(profileDatetime(iprof))

          ! Open the profile results file
          open(unit=userprofiles, &
               file=trim(userprofilesBase_fn) // '-' &
               // datetime_str(profileDatetime(iprof), datetime_sep='T', &
               time_sep='q'), &
               status='replace', action='write')
          ! Write the profile results file header
          call write_user_profiles_hdr(str_run_Datetime, str_CTD_Datetime, &
               profileDatetime(iprof))
          ! Write the profiles numbers, and close the profiles file
          call write_user_profiles_numbers(userprofiles, grid,          &
               (uptake%NO + uptake%NH + uptake%PC) * Redfield_C * 3600, &
               10.0d3 * (N%O(2:grid%M+1) + N%H(2:grid%M+1)) *           &
               w_upwell(1:grid%M))
          close(userprofiles)
       endif
    endif

    ! Check to see if we're in a Hoffmueller results output time step
    if (year == Hoff_yr .and. day == Hoff_day) then
       if (abs(day_time - Hoff_sec) < 0.5d0 * dt) then
          ! Write the profiles numbers
          call write_user_profiles_numbers(userHoff, grid,              &
               (uptake%NO + uptake%NH + uptake%PC) * Redfield_C * 3600, &
               10.0d3 * (N%O(2:grid%M+1) + N%H(2:grid%M+1)) *           &
               w_upwell(1:grid%M))
          ! Add empty line as separator
          write(userHoff, *)

       endif
    endif
  end subroutine write_user_profiles


  subroutine write_user_profiles_numbers(unit, grid, UptakeC, Nflux)
    ! Write the profiles numbers.  This is broken out to reduce code
    ! duplications.
    use precision_defs, only: dp
    use io_unit_defs, only: userprofiles
    use unit_conversions, only: KtoC
    use grid_mod, only: grid_
    implicit none
    ! Arguments:
    integer :: &
         unit  ! I/O unit to write to
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    real(kind=dp), intent (in) :: &
         UptakeC(1:), &  ! Carbon uptake [uM C/s]
         Nflux(1:)       ! Phytoplankton advection profile [umol m-2 s-1]
    ! Local variables:
    integer :: i  ! Loop index over grid depth
    
    ! Write the profile values at the surface.  Eddy diffusivity
    ! arrays don't have values there, so write zeros for them
    
    write(unit, 600) grid%d_g(0), UptakeC(1), Nflux(1)
    
    ! Write the profile values at the interior grid points
    
    do i = 1, grid%M
       write(unit, 600) grid%d_g(i), UptakeC(i), Nflux(i)
    enddo
    
    ! Write the values at the bottom grid boundary.  Some quantities are
    ! not defined there, so use their values at the Mth grid point.
    
    write(unit, 600) grid%d_g(grid%M+1), UptakeC(grid%M), Nflux(grid%M)
600 format(f7.3, 80(2x, f8.4))
  end subroutine write_user_profiles_numbers


end module user_output
