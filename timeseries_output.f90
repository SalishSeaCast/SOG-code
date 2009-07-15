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
  !   init_timeseries_output -- Get the names of the time series output
  !                             files, open them, and write their headers.
  !
  !   write_std_timeseries -- Write results of the current time step to the
  !                           time series files.
  !
  !   timeseries_output_close -- Close the time series output files.

  implicit none

  private
  public :: &
       ! Subroutines:
       init_timeseries_output, write_std_timeseries, timeseries_output_close

contains
  
  subroutine init_timeseries_output(codeId, str_run_Datetime, CTD_Datetime)
    ! Get the names of the time series output files from stdin using
    ! getpar(), open them, and write their headers.
    use io_unit_defs, only: std_phys_timeseries, std_bio_timeseries, &
         user_phys_timeseries, user_bio_timeseries
    use input_processor, only: getpars
    use datetime, only: datetime_, datetime_str
    use user_output, only: write_user_phys_timeseries_hdr, &
         write_user_bio_timeseries_hdr
    implicit none
    ! Arguments:
    character(len=70), intent(in) :: &
         codeId  ! Code identity string
    character(len=19), intent(in) :: &
         str_run_Datetime  ! Date/time of code run
    type(datetime_), intent(in) :: &
         CTD_Datetime  ! Date/time of CTD profile that initialized the run
    ! Local variables:
    character(len=80) :: &
         fn  ! File name to open
    type(datetime_) :: &
         start_Datetime  ! Midnight of the run start day
    ! Temporary storage for formated datetime strings.  Needed to work around
    ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
    ! calls in write statements
    character(len=19) :: &
         str_start_Datetime, &  ! Midnight of start day as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    ! Convert the initial CTD profile date/time to a string
    str_CTD_Datetime = datetime_str(CTD_Datetime)
    ! Midnight of the date of the initial CTD profile is t=0
    start_Datetime = CTD_Datetime
    start_Datetime%hr = 0
    start_Datetime%min = 0
    start_Datetime%sec = 0
    str_start_Datetime = datetime_str(start_Datetime)

    ! Get time series output file names from infile, open them, and
    ! write their headers
    !
    ! Standard physics model time series results
    fn = getpars("std_phys_ts_out")
    open(unit=std_phys_timeseries, file=fn, status="replace", action="write")
    call write_std_phys_timeseries_hdr(codeId, str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime)
    ! User physics model time series results
    fn = getpars("user_phys_ts_out")
    open(unit=user_phys_timeseries, file=fn, status="replace", action="write")
    call write_user_phys_timeseries_hdr(codeId, str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime)
    ! Standard biology model time series results
    fn = getpars("std_bio_ts_out")
    open(unit=std_bio_timeseries, file=fn, status="replace", action="write")
    call write_std_bio_timeseries_hdr(codeId, str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime)
    ! User biology model time series results
    fn = getpars("user_bio_ts_out")
    open(unit=user_bio_timeseries, file=fn, status="replace", action="write")
    call write_user_bio_timeseries_hdr(codeId, str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime)
  end subroutine init_timeseries_output


  subroutine write_std_phys_timeseries_hdr(codeId, str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime)
    ! Write standard physics model time series results file header.
    !
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the write_user_phys_timeseries() subroutine in the    !!!
    ! !!! user_output module for exploratory, special, debugging,   !!!
    ! !!! etc. output !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in   !!!
    ! !!! write_std_timeseries(), or compareSOG plotting will fail. !!!
    use io_unit_defs, only: std_phys_timeseries
    implicit none
    ! Arguments:
    character(len=70), intent(in) :: &
         codeId  ! Code identity string
    character(len=19), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_start_Datetime, &  ! Midnight of start day as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(std_phys_timeseries, 100) trim(codeId), str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime
100 format("! SOG code standard time series output from physics model"/,     &
         "! Time series of iteration count for each time step; ",            &
         "mixing layer depth;"/,                                             &
         "! velocity components, temperature, & salinity; at surface, ",     &
         "and averaged over"/,                                               &
         "! top 3 m of water column"/,                                       &
         "and surface par"/,                                                 &
         "*FromCode: ", a/,                                                  &
         "*RunDateTime: ", a/,                                               &
         "*InitialCTDDateTime: ", a/,                                        &
         "*FieldNames: time, iteration count, mixing layer depth, ",         &
         "surface u velocity, 3 m avg u velocity, ",                         &
         "surface v velocity, 3 m avg v velocity, ",                         &
         "surface temperature, 3 m avg temperature, ",                       &
         "surface salinity, 3 m avg salinity, ",                             &
         "surface PAR"/,                                                     &
         "*FieldUnits: hr since ", a, " LST, None, m, m/s, m/s, m/s, m/s, ", &
         "deg C, deg C, None, None, W/m2"/,                                  &
         "*EndOfHeader")
  end subroutine write_std_phys_timeseries_hdr


  subroutine write_std_bio_timeseries_hdr(codeId, str_run_Datetime, &
       str_CTD_Datetime, str_start_Datetime)
    ! Write standard biology model time series results file header.
    !
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the write_user_bio_timeseries() subroutine in the     !!!
    ! !!! user_output module for exploratory, special, debugging,   !!!
    ! !!! etc. output !!!
    !
    ! !!! The *FieldNames, and *FieldUnits parts of the header must !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! write_timeseries(), or compareSOG plotting will fail. !!!
    use io_unit_defs, only: std_bio_timeseries
    implicit none
    ! Arguments:
    character(len=70), intent(in) :: &
         codeId  ! Code identity string
    character(len=19), intent(in) :: &
         str_run_Datetime,   &  ! Date/time of code run as a string
         str_start_Datetime, &  ! Midnight of start day as a string
         str_CTD_Datetime       ! CTD profile date/time as a string

    write(std_bio_timeseries, 200) trim(codeId), str_run_Datetime, &
         str_CTD_Datetime, str_start_Datetime
200 format("! SOG code standard time series output from biology model"/,     &
         "! Time series of nitrate, ammonium, & silicon ",                   &
         "concentration; and "/,                                             &
         "! phytoplankton (micro & nano & pico) & microzooplankton biomass, ",      &
         "at surface, and "/,                                                &
         "! averaged over top 3 m of water column; and detritus (DON, ",     &
         "PON, refractory N, and "/,                                         &
         "! biogenic Si) at 20 m depth."/,                                   &
         "*FromCode: ", a/,                                                  &
         "*RunDateTime: ", a/,                                               &
         "*InitialCTDDateTime: ", a/,                                        &
         "*FieldNames: time, ",                                              &
         "surface nitrate concentration, 3 m avg nitrate concentration, ",   &
         "surface ammonium concentration, 3 m avg ammonium concentration, ", &
         "surface silicon concentration, 3 m avg silicon concentration, ",   &
         "surface micro phytoplankton biomass, ",                            &
         "3 m avg micro phytoplankton biomass, ",                            &
         "surface nano phytoplankton biomass, ",                             &
         "3 m avg nano phytoplankton biomass, ",                             &
         "surface pico phytoplankton biomass, ",                             &
         "3 m avg pico phytoplankton biomass, ",                             &
         "surface micro zooplankton biomass, ",                              &
         "3 m avg micro zooplankton biomass, ",                              &
         "DON detritus at 20 m, PON detritus at 20 m, ",                     &
         "refractory N detritus at 20 m, biogenic Si detritus at 20 m"/,     &
         "*FieldUnits: hr since ", a, " LST, uM N, uM N, uM N, uM N, ",      &
         "uM, uM, uM N, uM N, uM N, uM N, uM N, uM N, uM N, uM N, uM N, uM N, uM N, ",   &
         "uM"/,                                                              &
         "*EndOfHeader")
  end subroutine write_std_bio_timeseries_hdr
  

  subroutine write_std_timeseries(time, grid, &
       ! Variables for standard physics model output
       iter_cnt, h, U, V, T, S, Ipar, &
       ! Variables for standard biology model output
       NO, NH, Si, Pmicro, Pnano, Ppico, Z, D_DON, D_PON, D_refr, D_bSi)
    ! Write results of the current time step to the time series files.
    use precision_defs, only: dp
    use io_unit_defs, only: std_phys_timeseries, std_bio_timeseries
    use unit_conversions, only: KtoC
    use grid_mod, only: grid_, depth_average, interp_value
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: &
         time  ! [hr aftr start midnight]
    type(grid_), intent(in) :: &
         grid  ! Grid arrays
    integer, intent(in) :: &
         iter_cnt  ! Timestep iteration count
    real(kind=dp), intent(in) :: &
         h, &  ! Mixing layer depth [m]
         Ipar  ! Surface PAR
    real(kind=dp), dimension(0:), intent(in) :: &
         U,      &  ! Cross-strait velocity component profile [m/s]
         V,      &  ! Along-strait velocity component profile [m/s]
         T,      &  ! Temperature profile [K]
         S,      &  ! Salinity profile [-]
         NO,     &  ! Nitrate conc profile [uM N]
         NH,     &  ! Ammonium conc profile [uM N]
         Si,     &  ! Silicon conc profile [uM]
         Pmicro, &  ! Micro phytoplankton biomass profile [uM N]
         Pnano,  &  ! Nano phytoplankton biomass profile [uM N]
         Ppico,  &  ! Pico phytoplankton biomass profile [uM N]
         Z,      &  ! Micro zooplankton biomass profile [uM N]
         D_DON,  &  ! Dissolved organic nitrogen detritus profile [uM N]
         D_PON,  &  ! Particulate organic nitrogen detritus profile [uM N]
         D_refr, &  ! Refractory nitrogen detritus profile [uM N]
         D_bSi      ! Biogenic silicondetritus profile [uM]
    ! Local variables:
    real(kind=dp) :: &
         U_avg_3m,       &   ! U velocity averaged over top 3 m [m/s]
         V_avg_3m,       &   ! Vvelocity averaged over top 3 m [m/s]
         T_avg_3m,       &   ! Temperature averaged over top 3 m [deg C]
         S_avg_3m,       &   ! Salinity averaged over top 3 m [-]
         NO_avg_3m,      &   ! Nitrate averaged over top 3 m [-]
         NH_avg_3m,      &   ! Ammonium averaged over top 3 m [-]
         Si_avg_3m,      &   ! Silicon averaged over top 3 m [-]
         Pmicro_avg_3m,  &   ! Micro phytoplankton averaged over top 3 m [-]
         Pnano_avg_3m,   &   ! Nano phytoplankton averaged over top 3 m [-]
         Ppico_avg_3m,   &   ! Pico phytoplankton averaged over top 3 m [-]
         Z_avg_3m,       &   ! Micro zooplankton averaged over top 3 m [-]
         D_DON_20m,      &   ! DON detritus at 20 m [uM N]
         D_PON_20m,      &   ! PON detritus at 20 m [uM N]
         D_refr_20m,     &   ! Refractory N detritus at 20 m [uM N]
         D_bSi_20m           ! Biogenic silicon detritus at 20 m [uM N]
    integer :: j_below   ! Index of result found by interp_value()

    ! Write standard physics model time series results.
    !
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the write_user_phys_timeseries() subroutine in the    !!!
    ! !!! user_output module for exploratory, special, debugging,   !!!
    ! !!! etc. output !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    U_avg_3m = depth_average(U, 0.0d0, 3.0d0)
    V_avg_3m = depth_average(V, 0.0d0, 3.0d0)
    T_avg_3m = depth_average(KtoC(T), 0.0d0, 3.0d0)
    S_avg_3m = depth_average(S, 0.0d0, 3.0d0)
    write(std_phys_timeseries, 100) time, iter_cnt, h, U(0), U_avg_3m, &
         V(0), V_avg_3m, KtoC(T(0)), T_avg_3m, S(0), S_avg_3m, Ipar
100 format(f10.4, 2x, i3, 80(2x, f8.4))

    ! Write standard biology model time series results file header.
    !
    ! !!! Please don't change this unless you have a good reason to !!!
    ! !!! Use the write_user_bio_timeseries() subroutine in the     !!!
    ! !!! user_output module for exploratory, special, debugging,   !!!
    ! !!! etc. output !!!
    !
    ! !!! This write statement must be kept in sync with the *FieldNames, !!!
    ! !!! and *FieldUnits parts of the header in timeseries_output_open() !!!
    ! !!! be kept in sync with the appropriate write statement in  !!!
    ! !!! above, or compareSOG plotting will fail. !!!
    NO_avg_3m = depth_average(NO, 0.0d0, 3.0d0)
    NH_avg_3m = depth_average(NH, 0.0d0, 3.0d0)
    Si_avg_3m = depth_average(Si, 0.0d0, 3.0d0)
    Pmicro_avg_3m = depth_average(Pmicro, 0.0d0, 3.0d0)
    Pnano_avg_3m = depth_average(Pnano, 0.0d0, 3.0d0)
    Ppico_avg_3m = depth_average(Ppico, 0.0d0, 3.0d0)
    Z_avg_3m = depth_average(Z, 0.0d0, 3.0d0)
    call interp_value(20.0d0, 0, grid%d_g, D_DON, D_DON_20m, j_below)
    call interp_value(20.0d0, 0, grid%d_g, D_PON, D_PON_20m, j_below)
    call interp_value(20.0d0, 0, grid%d_g, D_refr, D_refr_20m, j_below)
    call interp_value(20.0d0, 0, grid%d_g, D_bSi, D_bSi_20m, j_below)
    write(std_bio_timeseries, 200) time, NO(0), NO_avg_3m, NH(0), NH_avg_3m, &
         Si(0), Si_avg_3m, Pmicro(0), Pmicro_avg_3m, Pnano(0), Pnano_avg_3m, &
         Ppico(0), Ppico_avg_3m, &
         Z(0), Z_avg_3m, D_DON_20m, D_PON_20m, D_refr_20m, D_bSi_20m
200 format(f10.4, 80(2x, f8.4))
  end subroutine write_std_timeseries


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
