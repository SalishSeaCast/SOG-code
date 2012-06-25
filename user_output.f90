module user_output
  ! Subroutines for user-specified output in the SOG code.
  !
  ! Public Subroutines:
  !
  !   write_user_phys_timeseries_hdr --
  !
  !   write_user_bio_timeseries_hdr --
  !
  !   write_user_timeseries --
  !
  !   write_user_profiles --

  implicit none

  private
  public :: &
       ! Subroutines:
       write_user_phys_timeseries_hdr, write_user_bio_timeseries_hdr, &
       write_user_chem_timeseries_hdr, write_user_timeseries, &
       write_user_profiles

contains

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
         !SEA         "*FieldNames: time, Avg (0-3m) microplankton biomass, "       &
    "*FieldNames: time"/,       &
         !SEA         "Avg (0-3m) nanoplankton biomass"/, &
         !SEA         "*FieldUnits: hr since ", a, " LST, uM N, uM N"/, &
         "*FieldUnits: hr since ", a, " LST"/, &
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
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: &
         time  ! [hr aftr start midnight]
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
    write(user_bio_timeseries, 200) time
200 format(f10.4)

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


  subroutine write_user_profiles
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
    ! !!! changes.                                                  !!!
    use precision_defs, only: dp
    use io_unit_defs, only: profiles, haloclines
    use datetime, only: calendar_date, clock_time, datetime_str
    use unit_conversions, only: KtoC
    use grid_mod, only: grid_
    implicit none
    ! Recommended argument declarations:
    ! character(len=19), intent(in) :: str_run_Datetime  ! Date/time of code run
    ! character(len=19), intent(in) :: str_CTD_Datetime  ! Date/time of CTD init
    ! integer, intent(in) :: year, day
    ! real(kind=dp), intent(in) :: day_time, dt ! can't expect exact time match
    ! type(grid_), intent(in) :: grid

  end subroutine write_user_profiles

end module user_output
