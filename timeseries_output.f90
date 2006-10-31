! $Id$
! $Source$

module timeseries_output
  ! Subroutines to write time series results to files specified in run
  ! parameters file.  The output files contain headers so that they
  ! can be used by the compareSOG package to create comparison
  ! graphs.
  ! 
  ! Public subroutines:
  !
  ! init_timeseries_output(codeId, str_runDatetime, CTDdatetime)
  !   -- Get the names of the time series output files, open them, and
  !      write their headers.
  !
  ! write_timeseries(time, &
  !    ! Variables for standard physical model output
  !    h, T, S, iter_cnt, &
  !    ! User-defined physical model output variables
  !    &
  !    ! Variables for standard biological model output
  !    NO, NH, Pmicro, Pnano, Z &
  !    , D_bins, detritus, &
  !    ! User-defined biological model output variables
  !    )
  !   -- Write results of the current time step to the time series files.
  !
  ! subroutine timeseries_output_close
  !   -- Close the time series output files.

  use io_unit_defs, only: std_phys_timeseries, std_bio_timeseries, &
       user_phys_timeseries, user_bio_timeseries

  implicit none

  private
  public :: init_timeseries_output, write_timeseries, &
       timeseries_output_close

contains
  
  subroutine init_timeseries_output(codeId, str_runDatetime, CTDdatetime)
    ! Get the names of the time series output files from stdin using
    ! getpar(), open them, and write their headers.
    use input_processor, only: getpars
    use datetime, only: datetime_, datetime_str
    implicit none
    ! Arguments:
    character(len=70), intent(in) :: codeId             ! Code identity string
    character(len=19), intent(in) :: str_runDatetime    ! Date/time of code run
    ! Date/time of CTD profile that initialized the run
    type(datetime_), intent(in)   :: CTDdatetime      
    ! Local variable:
    character(len=80) :: fn           ! File name to open
    type(datetime_) :: startDatetime
    ! Temporary storage for formated datetime strings.  Needed to work around
    ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
    ! calls in write statements
    character(len=19) :: str_startDatetime  ! Midnight of start day
    character(len=19) :: str_CTDdatetime    ! CTD profile date/time

    ! Convert the initial CTD profile date/time to a string
    str_CTDdatetime = datetime_str(CTDdatetime)
    ! Midnight of the date of the initial CTD profile is t=0
    startDatetime = CTDdatetime
    startDatetime%hr = 0
    startDatetime%min = 0
    startDatetime%sec = 0
    str_startDatetime = datetime_str(startDatetime)

    ! Standard physics model time series results
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the user_phys_timeseries file below exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    fn = getpars("std_phys_ts_out")
    open(unit=std_phys_timeseries, file=fn, status="replace", action="write")
    write(std_phys_timeseries, 100) trim(codeId), str_runDatetime, &
         str_CTDdatetime, str_startDatetime
100 format("! Standard time series output from physics model"/,         &
         "! Time series of mixing layer depth; temperature, and, ",     &
         "salinity"/,                                                   &
         "! at surface, and averaged over top 5 m of water column;  ",  &
         " and interation "/,                                           &
         "! counts for each time step"/,                                & 
         "*FromCode: ", a/,                                             &
         "*RunDateTime: ", a/,                                          &
         "*InitialCTDDateTime: ", a/,                                   &
         "*FieldNames: time, iteration count, mixing layer depth, ",    &
         "surface temperature, surface salinity"/,                      &
         "*FieldUnits: hr since ", a, " LST, None, m, deg C, None"/,    &
         "*EndOfHeader")


    ! User physics model time series results
    ! !!! This is the place to add exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    ! !!! Please don't commit this file if you only make personal !!!
    ! !!! changes here. !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    fn = getpars("user_phys_ts_out")
    open(unit=user_phys_timeseries, file=fn, status="replace", action="write")
    write(user_phys_timeseries, 101) trim(codeId), str_runDatetime, &
         str_CTDdatetime, str_startDatetime
101 format("! User-defined time series output from physics model"/,    &
         "! Time series of ...", &
         "*FromCode: ", a/,                                            &
         "*RunDateTime: ", a/,                                         &
         "*InitialCTDDateTime: ", a/,                                  &
!SEA         "*FieldNames: time, ut, vt"/, &
         "*FieldNames: time"/, &
         "*FieldUnits: hr since ", a, " LST"/, &
         "*EndOfHeader")


    ! Standard biology model time series results
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the user_phys_timeseries file below exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    fn = getpars("std_bio_ts_out")
    open(unit=std_bio_timeseries, file=fn, status="replace", action="write")
    write(std_bio_timeseries, 102) trim(codeId), str_runDatetime, &
         str_CTDdatetime, str_startDatetime
102 format("! Standard time series output from biology model"/,         &
         "! Time series of nitrate, ammonium, and silicon ",            &
         "concentration, and "/,                                        &
         "! phytoplankton (micro & nano) biomass, at surface, and ",    &
         "averaged over ",                                              &
         "! top 5 m of water column; and detritus (dissolved, ",        &
         "sinking, and lost)"/,                                         &
         "! at 20 m depth."/,                                           &
         "*FromCode: ", a/,                                             &
         "*RunDateTime: ", a/,                                          &
         "*InitialCTDDateTime: ", a/,                                   &
         "*FieldNames: time, surface nitrate concentration, ",          &
         "surface ammonium concentration, ",                            &
         "surface silicon concentration, ",                             &
         "surface micro phytoplankton biomass, ",                       &
         "surface nano phytoplankton biomass, ",                        &
         "surface micro zooplankton biomass, ",                         &
         "dissolved detritus at 20 m, ",                                &
         "sinking detritus at 20 m, ",                                  &
         "lost detritus at 20 m"/,                                      &
         "*FieldUnits: hr since ", a, " LST, uM N, uM N, uM, uM N, ",   &
         "uM N, uM N, uM N, uM N, uM N"/,                                     &
         "*EndOfHeader")


    ! User biology model time series results
    ! !!! This is the place to add exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    ! !!! Please don't commit this file if you only make personal !!!
    ! !!! changes here. !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    fn = getpars("user_bio_ts_out")
    open(unit=user_bio_timeseries, file=fn, status="replace", action="write")
    write(user_bio_timeseries, 103) trim(codeId), str_runDatetime, &
         str_CTDdatetime, str_startDatetime
103 format("! User-defined time series output from biology model"/,    &
         "! Time series of ...", &
         "*FromCode: ", a/,                                            &
         "*RunDateTime: ", a/,                                         &
         "*InitialCTDDateTime: ", a/,                                  &
!SEA         "*FieldNames: time, Avg (0-3m) microplankton biomass, "       &
         "*FieldNames: time"/,       &
!SEA         "Avg (0-3m) nanoplankton biomass"/, &
!SEA         "*FieldUnits: hr since ", a, " LST, uM N, uM N"/, &
         "*FieldUnits: hr since ", a, " LST"/, &
         "*EndOfHeader")

  end subroutine init_timeseries_output


  subroutine write_timeseries(time, grid, &
       ! Variables for standard physical model output
       iter_cnt, h, T, S, &
       ! User-defined physical model output variables
!SEA       dPdx_b, dPdy_b, unow, vnow, u, v, &
       ! Variables for standard biological model output
       NO, NH, Si, Pmicro, Pnano, Z, remin_Detritus, sink_Detritus, &
       BSi_Detritus &
!!$       &
       ! User-defined biological model output variables
       )

    ! Write results of the current time step to the time series files.

    use precision_defs, only: dp
    use unit_conversions, only: KtoC
    use grid_mod, only: grid_, interp_value, depth_average
    use physics_model, only: ut, vt
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: time                ! [hr aftr start midnight]
    type(grid_), intent(in) :: grid                  ! Grid arrays
    integer, intent(in) :: iter_cnt                  ! Timestep iteration count
!SEA    real(kind=dp), intent(in) :: h, &                   ! Mixed layer depth [m]
    real(kind=dp), intent(in) :: h,                    ! Mixed layer depth [m]
!SEA         unow, &
!SEA         vnow
    real(kind=dp), dimension(0:), intent(in) :: &
         T, &               ! Temperature [K]
         S, &               ! Salinity [-]
!SEA         u, &
!SEA         v, &
         NO, &              ! Nitrate conc [uM N]
         NH, &              ! Ammonium conc [uM N]
         Si, &              ! Silicon conc [uM Si]
         Pmicro, &          ! Micro phytoplankton biomass [uM N]
         Pnano, &           ! Nano phytoplankton biomass [uM N]
         Z, &               ! Micro zooplankton biomass [uM N]
         remin_Detritus, &  ! Remineralized detritus biomass [uM N]
         sink_Detritus, &   ! Sinking detritus biomass [uM N]
         BSi_Detritus       ! Biogenic Si detritus biomass [uM Si]
!SEA    real(kind=dp), dimension(1:), intent(in) :: &
!SEA         dPdx_b, &
!SEA         dPdy_b


    ! Local variables:
    real(kind=dp) :: &
         remin_D_20m, &  ! Remineralized detritus at 20 m [uM N]
         sink_D_20m, &   ! Sinking detritus at 20 m [uM N]
         BSi_D_20m, &   ! Mortality detritus at 20 m [uM N]
!SEA         avg_Pmicro_0_3, & ! Average micro plankton biomass 0-3 m [uM N]
!SEA         avg_Pnano_0_3  ! Average nano plankton biomass 0-3 m [uM N]
    integer :: j_below   ! Index of result found by interp_value()
   

    ! Standard physics model time series results
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the user_phys_timeseries file below exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    !
    ! time, timestep iteration count, mixed layer depth, 
    ! surface temperature and salinity
    write(std_phys_timeseries, 100) time, iter_cnt, h, KtoC(T(0)), S(0)
100 format(f10.4, 2x, i3, 80(2x, f8.4))


    ! User physics model time series results
    ! !!! This is the place to add exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    ! !!! Please don't commit this file if you only make personal !!!
    ! !!! changes here. !!!
    !v
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    !
    ! time
!SEA    write(user_phys_timeseries, 102) time, & 
    write(user_phys_timeseries, 101) time,
!SEA    unow, vnow, u(1), v(1), ut%new(1), vt%new(1), dPdx_b(1)*1e5, dPdy_b(1)*1e5
101 format(f10.4)


    ! Standard biology model time series results
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the user_phys_timeseries file below exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    !
    ! time, surface nitrate, ammonium, and concentrations, 
    ! surface biomass of micro and nano phytoplankton,
    ! biomasses of detritus at 20 m depth
    call interp_value(20.0d0, grid%d_g, remin_Detritus, remin_D_20m, j_below)
    call interp_value(20.0d0, grid%d_g, sink_Detritus, sink_D_20m, j_below)
    call interp_value(20.0d0, grid%d_g, Bsi_Detritus, Bsi_D_20m, j_below)
    write(std_bio_timeseries, 102) time, NO(0), NH(0), Si(0),      &
         Pmicro(0), Pnano(0), Z(0), remin_D_20m, sink_D_20m, Bsi_D_20m
102 format(f10.4, 80(2x, f8.4))


    ! User biology model time series results
    ! !!! This is the place to add exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    ! !!! Please don't commit this file if you only make personal !!!
    ! !!! changes here. !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    !
    ! time, avg_Pmicro_0_3, avg_Pnano_0_3
!SEA    avg_Pmicro_0_3 = depth_average(Pmicro, 0.d0, 3.d0)
!SEA    avg_Pnano_0_3 = depth_average(Pnano, 0.d0, 3.d0)
!SEA    write(user_bio_timeseries, 102) time,avg_Pmicro_0_3, avg_Pnano_0_3
    write(user_bio_timeseries, 103) time
103 format(f10.4)

  end subroutine write_timeseries


  subroutine timeseries_output_close
    ! Close the time series output files.
    use io_unit_defs, only: std_phys_timeseries, user_phys_timeseries, &
         std_bio_timeseries, user_bio_timeseries
    implicit none
    close(std_phys_timeseries)
    close(user_phys_timeseries)
    close(std_bio_timeseries)
    close(user_bio_timeseries)
  end subroutine timeseries_output_close

end module timeseries_output
