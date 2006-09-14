module forcing
  ! package to return forcing values and allow them to be shifted in
  ! time or multiply them up or set them constant

  ! we will assume that all data starts on Jan 1, 2001.  Wind data is in
  ! a single array with 1 hour per line, cf, atemp and humid are in array
  ! of day and time, and the river flows are one value per day.

  use precision_defs, only: dp, sp

  implicit none
  private
  public :: read_variation, read_forcing, get_forcing

  ! control variable
  TYPE :: vary_quantity
     logical :: enabled ! is the parameter (wind, cf, rivers) to be varied
     real(kind=dp) :: shift ! use wind/rivers shifted in time (in years)
     real(kind=dp) :: fraction ! multiply wind/rivers strength by this amount
     logical :: fixed ! use a fixed value 
     ! (if this is set to true, shift and fraction are not used and value is)
     real(kind=dp) :: value ! use this value if fixed is set true
  END TYPE vary_quantity

  TYPE :: vary_forcing
     TYPE(vary_quantity) :: wind, cf, rivers
  END TYPE vary_forcing

  TYPE (vary_forcing) :: vary

  ! Number of years of data (maximum for all the files)
  integer, parameter :: noyears=6

  ! Data variables
  real(kind=dp) :: wind_eastnorth(noyears*366*24), &
       wind_northwest(noyears*366*24)
  real(kind=sp) :: cf(noyears*366,24), atemp(noyears*366,24), &
       humid(noyears*366,24)
  real(kind=sp) :: Qriver(noyears*366), Eriver(noyears*366)

contains

  subroutine read_variation

    ! subroutine to read variation parameters from the inputfile.  

    IMPLICIT NONE

    ! read in vary parameters

    call read_vary ("vary%wind",vary%wind)

    call read_vary ("vary%cf",vary%cf)
    if (vary%cf%enabled .and. .not.vary%cf%fixed) then
       if (vary%cf%shift /= int(vary%cf%shift)) then
          write (*,*) "non-integer cf shifts not implemented"
          stop
       endif
    endif
    
    call read_vary ("vary%rivers",vary%rivers)
    if (vary%wind%enabled) then
       if  (vary%wind%fixed) then
          write (*,*) 'Fixed wind magnitude has not been implemented'
          stop
       elseif (vary%wind%shift /= dint(vary%wind%shift)) then
          write (*,*) 'Non-full year shifts have not been implemented ',&
               'for the wind'
          stop
       endif
    endif

  END SUBROUTINE read_variation

  subroutine read_forcing (wind_n, met_n, river_n)

    use unit_conversions, only: CtoK

    implicit none

    integer, intent(in) :: wind_n, met_n, river_n ! length of input files
    ! should switch this to being read in.

    ! Local variables:
    integer :: ic, j, para, stn, year, day, month
    real(kind=sp) :: hour

    
    ! WIND
    open(unit = 12, file = "../sog-forcing/wind/SHcompRotFmt.dat", &
         status = "OLD", action = "READ")

    do ic = 1,wind_n
       read(12, *) day, month, year, hour, &
            wind_eastnorth(ic), wind_northwest(ic)
    enddo
    ! checking day, month, year, time here would be a GOOD THING
    close (12)

    ! Cloud Fraction

    open(12, file="../sog-forcing/met/YVRCombCF.dat", &
         status = "OLD", action = "READ")
    do ic = 1, met_n
       read(12, *) stn, year, month, day, para, (cf(ic,j), j = 1, 24)
    enddo
    ! checking day, month, year, stn and para here would be a GOOD THING
    close (12)


    ! Air Temperature

    open(12, file="../sog-forcing/met/YVRCombATemp.dat", &
         status = "OLD", action = "READ")
    do ic = 1, met_n
       read(12, *) stn, year, month, day, para, (atemp(ic,j), j = 1, 24)
       do j= 1, 24
!          atemp(ic,j) = CtoK(atemp(ic,j)/10.)
          atemp(ic,j) = atemp(ic,j) / 10. + 273.15
       enddo
    enddo
    ! checking day, month, year, stn and para here would be a GOOD THING
    close(12)

    ! Humidity

    open(12, file="../sog-forcing/met/YVRCombHum.dat", &
         status = "OLD", action = "READ")
    do ic = 1, met_n
       read(12, *) stn, year, month, day, para, (humid(ic,j), j = 1, 24)
    enddo
    ! checking day, month, year, stn and para here would be a GOOD THING
    close(12)

    ! Fraser River
    open(12, file="../sog-forcing/rivers/Fraser_2001_2005.dat", &
         status = "OLD", action = "READ")
    do ic = 1, river_n
       read(12, *) year, month, day, Qriver(ic)
    enddo
    ! checking day, month, year, here would be a GOOD THING
    close(12)

    ! Englishman River
    open(12, file="../sog-forcing/rivers/eng200123456.dat", &
         status = "OLD", action = "READ")
    do ic = 1, river_n
       read(12, *) year, month, day, Eriver(ic)
    enddo
    ! checking day, month, year, here would be a GOOD THING
    close(12) 

  end subroutine read_forcing
       
  subroutine get_forcing (year, day, day_time, &
       Qinter, Einter, cf_value, atemp_value, humid_value, unow, vnow)

    ! given the year, day and day_time find the physical forcing values
    ! taking into account any shifts or other changes set by vary

    integer, intent(in) :: year ! current year
    integer, intent(in) :: day ! current year day
    real(kind=dp), intent(in) :: day_time ! current time of day in seconds
    
    real(kind=sp), intent(out) :: Qinter, Einter ! values of river flows
    real(kind=sp), intent(out) :: cf_value     ! cloud fraction
    real(kind=sp), intent(out) :: atemp_value  ! air temperature
    real(kind=sp), intent(out) :: humid_value  ! humidity
    real(kind=dp), intent(out) :: unow, vnow ! wind components

    ! we will assume that all data starts on Jan 1, 2001.  Wind data is in
    ! a single array with 1 hour per line, cf, atemp and humid are in array
    ! of day and time, and the river flows are one value per day.

    ! local variable

    integer:: accul_day ! accumulated day number with Jan 1, 2001 being 1
    integer:: shift_day, shift_year ! shifted day and year
    integer:: j ! index to hour, 1 am has j=2 
    integer:: index_day ! index into wind data for the end of the day before
    real(kind=dp) :: hr_frac ! fraction of the hour

    ! RIVERS

    if (vary%rivers%enabled .and. vary%rivers%fixed) then
       Qinter = SNGL(vary%rivers%value)
       Einter = 0.
    else
       if (vary%rivers%enabled) then
          call shift_time (day, year, vary%rivers%shift, &
               shift_day, shift_year)
          accul_day = accum_day(shift_year, shift_day)
       else
          accul_day = accum_day(year,day)
       endif 
       Qinter = (day_time * Qriver(accul_day) & 
            + (86400 - day_time) * Qriver(accul_day-1)) / 86400.
       Einter = (day_time * Eriver(accul_day) &
            + (86400 - day_time) * Eriver(accul_day-1)) / 86400.
       if (vary%rivers%enabled) then
          Qinter = Qinter * SNGL(vary%rivers%fraction)
          Einter = Einter * SNGL(vary%rivers%fraction)
       endif
    endif

    ! MET DATA (cloud fraction, air temperature and humidity)

    ! index into met data array for the appropriate hr
    j = floor(day_time/3600.0) + 1  
    
    if (vary%cf%enabled) then
       if (vary%cf%fixed) then
          cf_value = sngl(vary%cf%value)
       else
          accul_day = accum_day(year+int(vary%cf%shift), day)
          cf_value = sngl(cf(accul_day,j)*vary%cf%fraction)
       endif
    else
       accul_day = accum_day(year, day)
       cf_value = cf(accul_day,j)
    endif

    accul_day = accum_day(year, day)
    
    atemp_value = atemp(accul_day,j)
    humid_value = humid(accul_day,j)

    ! WIND

    if (vary%wind%enabled) then ! note no fixed wind values
       accul_day = accum_day(year+int(vary%wind%shift), day)
    else
       accul_day = accum_day(year, day)
    endif

    ! wind for the previous day, last value is at
    index_day = (accul_day-1) * 24

    ! frac hour is 
    hr_frac = day_time/3600. - (j-1)
    ! index the wind at the start of the hour is
    j = index_day + j

    vnow = wind_northwest(j) * (1 - hr_frac) &
         + wind_northwest(j+1) * hr_frac
    unow = wind_eastnorth(j) * (1 - hr_frac) &
         + wind_eastnorth(j+1) * hr_frac

    if (vary%wind%enabled) then
       vnow = vnow * vary%wind%fraction
       unow = unow * vary%wind%fraction
    endif

  end subroutine get_forcing

  SUBROUTINE read_vary (quantity_string,quantity)
    
    use input_processor, only: getparl, getpard

    IMPLICIT NONE
    
    ! the leading part of the variable quantity to be read in (eg vary%wind)
    CHARACTER, INTENT(IN) :: quantity_string*(*) 
    TYPE (vary_quantity), INTENT(OUT) :: quantity

    quantity%enabled = getparl(quantity_string//"%enabled")
    if (quantity%enabled) then
       quantity%fixed = getparl(quantity_string//"%fixed")
       if(quantity%fixed) then
          quantity%value = getpard(quantity_string//"%value")
       else
          quantity%shift = getpard(quantity_string//"%shift")
          quantity%fraction = getpard(quantity_string//"%fraction")
       endif
    endif

  END SUBROUTINE read_vary

  integer function accum_day (year, day) result(day_met)

    ! calculates the accumulated day based on Jan 1, 2001 being day 1.

    INTEGER,INTENT(IN) :: day, year ! year day and year
    
    IF (year==2001) then
       day_met=day
    else if (year==2002) then
       day_met=day + 365
    else if (year==2003) then
       day_met=day + 730
    else if (year==2004) then
       day_met=day + 1095
    else if (year==2005) then
       day_met=day + 1461
    else if (year==2006) then
       day_met=day + 1826
    else
       write (*,*) "unknown year in met data", year
       stop
    endif

  end function accum_day

  subroutine shift_time (day, year, shift, &
       shift_day, shift_year)

    ! subroutine to shift by part of a year

    integer, intent(in) :: day, year ! current day and year
    real(kind=dp), intent(in) :: shift ! shift by fraction of a year

    ! shifted day and year
    integer, intent(out) :: shift_day, shift_year


    shift_day = day + floor(365*shift)
    shift_year = year
    do while ((shift_day > 365 .and. .not.leapyear(shift_year)) &
         .or. shift_day .gt. 366)
       if (leapyear(year)) then
          shift_day = shift_day - 366
          shift_year = shift_year + 1
       else
          shift_day = shift_day - 365
          shift_year = shift_year + 1
       endif
    enddo
    
  end subroutine shift_time

  LOGICAL function leapyear(year)
    
    ! determines if year is a leap-year
    
    INTEGER, INTENT(IN) :: year
    
    if (year.gt.1953.and.year.lt.2007) then
       if (year.eq.1956.or.year.eq.1960.or.year.eq.1964.or. &
            year.eq.1968.or.year.eq.1972.or.year.eq.1976.or. &
            year.eq.1980.or.year.eq.1984.or.year.eq.1988.or. &
            year.eq.1992.or.year.eq.1996.or.year.eq.2000.or. &
            year.eq.2004) then
          leapyear = .TRUE.
       else
          leapyear = .FALSE.
       endif
    else
       write (*,*) 'Out of range for leap-year calculation. ', &
            'Add more years to function leapyear'
       stop
    endif
    
  end function leapyear

end module forcing
