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
  ! timeseries_output_open(codeId, runDatetime, startDatetime)
  !    -- Get the names of the time series output files from stdin 
  !       using getpars(), open them, and write their headers.
  !
  ! write_timeseries(time, &
  !    ! Variables for standard physical model output
  !    h, T, S, iter_cnt, &
  !    ! User-defined physical model output variables
  !    &
  !    ! Variables for standard biological model output
  !    NO, NH, Pmicro, Pnano &
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
  public :: timeseries_output_open, write_timeseries, &
       timeseries_output_close

contains
  
  subroutine timeseries_output_open(codeId, runDatetime, startDatetime)
    ! Get the names of the time series output files from stdin using
    ! getpars(), open them, and write their headers.
    ! *** Does reading the file names from stdin belong here, or in
    ! *** the input processor?

    ! *** getpars should come from the input processor module,
    ! *** eventually
    use initial_sog, only: getpars
    implicit none
    ! Arguments:
    character(len=70), intent(in) :: codeId         ! Code identity string
    character(len=19), intent(in) :: runDatetime    ! Date/time of code run
    character(len=19), intent(in) :: startDatetime  ! Midnight of start day
    ! Local variable:
    character(len=80) :: fn           ! File name to open

    ! Standard physics model time series results
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the user_phys_timeseries file below exploratory, !!!
    ! !!! special, debugging, etc. output !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    fn = getpars("std_phys_ts_out", 1)
    open(unit=std_phys_timeseries, file=fn, status="replace", action="write")
    write(std_phys_timeseries, 100) trim(codeId), runDatetime, &
         startDatetime
100 format("! Standard time series output from physics model"/,         &
         "! Time series of mixing layer depth; temperature, and, ",     &
         "salinity"/,                                                   &
         "! at surface, and averaged over top 5 m of water column;  ",  &
         " and interation "/,                                           &
         "! counts for each time step"/,                                & 
         "*FromCode: ", a/,                                             &
         "*RunDateTime: ", a/,                                          &
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
    fn = getpars("user_phys_ts_out", 1)
    open(unit=user_phys_timeseries, file=fn, status="replace", action="write")
    write(user_phys_timeseries, 101) trim(codeId), runDatetime, &
         startDatetime
101 format("! User-defined time series output from physics model"/,    &
         "! Time series of ...", &
         "*FromCode: ", a/,                                            &
         "*RunDateTime: ", a/,                                         &
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
    fn = getpars("std_bio_ts_out", 1)
    open(unit=std_bio_timeseries, file=fn, status="replace", action="write")
    write(std_bio_timeseries, 102) trim(codeId), runDatetime, &
         startDatetime
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
         "*FieldNames: time, surface nitrate concentration, ",          &
         "surface ammonium concentration, ",                            &
         "surface silicon concentration, ",                             &
         "surface micro phytoplankton biomass, ",                       &
         "surface nano phytoplankton biomass, ",                        &
         "dissolved detritus at 20 m, ",                                &
         "sinking detritus at 20 m, ",                                  &
         "lost detritus at 20 m"/,                                      &
         "*FieldUnits: hr since ", a, " LST, uM N, uM N, uM, uM N, ",   &
         "uM N, uM N, uM N, uM N"/,                                     &
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
    fn = getpars("user_bio_ts_out", 1)
    open(unit=user_bio_timeseries, file=fn, status="replace", action="write")
    write(user_bio_timeseries, 103) trim(codeId), runDatetime, &
         startDatetime
103 format("! User-defined time series output from biology model"/,    &
         "! Time series of ...", &
         "*FromCode: ", a/,                                            &
         "*RunDateTime: ", a/,                                         &
         "*FieldNames: time"/, &
         "*FieldUnits: hr since ", a, " LST"/, &
         "*EndOfHeader")

  end subroutine timeseries_output_open


  subroutine write_timeseries(time, &
       ! Variables for standard physical model output
       iter_cnt, h, T, S, &
       ! User-defined physical model output variables
!!$       &
       ! Variables for standard biological model output
       NO, NH, Sil, Pmicro, Pnano, remin_Detritus, sink_Detritus, &
       mort_Detritus &
!!$       &
       ! User-defined biological model output variables
       )

    ! Write results of the current time step to the time series files.

    use precision_defs, only: dp
    use unit_conversions, only: KtoC
    use grid_mod, only: interp_d
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: time                ! [hr aftr start midnight]
    integer, intent(in) :: iter_cnt                  ! Timestep iteration count
    real(kind=dp), intent(in) :: h                   ! Mixed layer depth [m]
    real(kind=dp), dimension(0:), intent(in) :: T    ! Temperature [K]
    real(kind=dp), dimension(0:), intent(in) :: S    ! Salinity [-]
    real(kind=dp), dimension(0:), intent(in) :: NO   ! Nitrate conc [uM N]
    real(kind=dp), dimension(0:), intent(in) :: NH   ! Ammonium conc [uM N]
    real(kind=dp), dimension(0:), intent(in) :: Sil  ! Silicon conc [uM Si]
    ! Micro phytoplankton biomass [uM N]
    real(kind=dp), dimension(0:), intent(in) :: Pmicro
    ! Nano phytoplankton biomass [uM N]
    real(kind=dp), dimension(0:), intent(in) :: Pnano
    ! Remineralized detritus biomass [uM N]
    real(kind=dp), dimension(0:), intent(in) :: remin_Detritus
    ! Sinking detritus biomass [uM N]
    real(kind=dp), dimension(0:), intent(in) :: sink_Detritus
    ! Mortality detritus biomass [uM N]
    real(kind=dp), dimension(0:), intent(in) :: mort_Detritus

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
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    !
    ! time
    write(user_phys_timeseries, 101) time
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
    write(std_bio_timeseries, 102) time, NO(0), NH(0), Sil(0), &
         Pmicro(0), Pnano(0), interp_d(remin_Detritus, 20.0d0),   &
         interp_d(sink_Detritus, 20.0d0), interp_d(mort_Detritus, 20.0d0)
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
    ! time
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
