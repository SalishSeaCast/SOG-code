module forcing
  ! package to return forcing values and allow them to be shifted in
  ! time or multiply them up or set them constant

  ! Data is coming from historic data files. We have ctd data for some 70,80,90's.
  ! Metar data goes back to 50's and river data goes back 10's.  
  ! Read in 2 years of data at a time. Wind data is in a single array 
  ! with 1 hour per line, cf, atemp and humid are in array
  ! of day and time, and the river flows are one value per day.

  use precision_defs, only: dp, sp

  implicit none
  private
  public :: read_variation, read_forcing, get_forcing, &
! and for varying bottom and initial temperatures
       vary_forcing, vary, &
       ! Types (as required by new pg compiler)
       vary_quantity ! sub-type of vary_forcing


  ! control variable
  TYPE :: vary_quantity
     logical :: enabled ! is the parameter (wind, cf, rivers, T) to be varied
     real(kind=dp) :: shift ! use wind/rivers shifted in time (in years)
     real(kind=dp) :: fraction ! multiply wind/rivers strength by this amount
     logical :: fixed ! use a fixed value 
     ! (if this is set to true, shift and fraction are not used and value is)
     real(kind=dp) :: value ! use this value if fixed is set true
     real(kind=dp) :: addition ! add this value (after fraction is multiplied)
  END TYPE vary_quantity

  TYPE :: vary_forcing
     TYPE(vary_quantity) :: wind, cf, rivers, temperature
  END TYPE vary_forcing

  TYPE (vary_forcing) :: vary


  ! Number of years of data  (set to 2 years)
  !
  integer, parameter :: wind_n=((366+366)*24), & ! number of days used for wind array
       met_n=(366+366), &  !# days used for met arrays (CF,AT,HUM)
       river_n=(366+366)   !# days used for river array
  ! Data variables

  real(kind=dp) :: wind_eastnorth(wind_n), &
       wind_northwest(wind_n)
  real(kind=sp) :: cf(met_n,24), atemp(met_n,24), &
       humid(732,24)
  real(kind=sp) :: Qriver(river_n), Eriver(river_n)

  !The starting year (startyear + 1 more year of data)

  integer :: startyear

  !File names for Met, wind and river data
  character* 80 :: &
       Wind,   &   ! file name for wind data
       AT,     &   ! file name for Air temp  data
       Cloud,  &   ! file name for Cloud data
       Hum,    &   ! file name for Humidity data
       MajorRiv, &     ! file name for Major River data
       MinorRiv         ! file name for Minor River data
 
  ! Is there data for a minor river?
  logical :: isMinorRiv    

contains

  subroutine read_variation

    ! subroutine to read variation parameters from the inputfile.  

    IMPLICIT NONE

    ! read in vary parameters

    call read_vary ("vary%wind",vary%wind)

    if (vary%wind%enabled) then
       if  (vary%wind%fixed) then
          write (*,*) 'Fixed wind magnitude has not been implemented'
          call exit(1)
       else
          if (vary%wind%shift /= dint(vary%wind%shift)) then
             write (*,*) 'Non-full year shifts have not been implemented ',&
                  'for the wind'
             call exit(1)
          endif
          if (vary%wind%addition /= 0) then
             write (*,*) "addtion not implemented for wind"
             call exit(1)
          endif
       endif
    endif

    call read_vary ("vary%cf",vary%cf)
    if (vary%cf%enabled .and. .not.vary%cf%fixed) then
       if (vary%cf%shift /= int(vary%cf%shift)) then
          write (*,*) "non-integer cf shifts not implemented"
          call exit(1)
       endif
       if (vary%cf%addition /= 0) then
          write (*,*) "addition not impleemented for cf"
          call exit(1)
       endif
    endif

    call read_vary ("vary%rivers",vary%rivers)
    if (vary%rivers%enabled .and. .not.vary%rivers%fixed) then
       if (vary%rivers%addition /= 0 ) then
          write (*,*) "addition not implemented for rivers"
          call exit(1)
       endif
    endif

    call read_vary ("vary%temperature", vary%temperature)
    if (vary%temperature%enabled) then
       if (vary%temperature%fixed) then
          write (*,*) "Fixed temperature not enabled"
          call exit (1)
       else
          if (vary%temperature%shift /= 0) then
             write (*,*) "Fixed temperature not enabled"
             call exit (1)
          endif
          if (vary%temperature%fraction /= 1) then
             write (*,*) "Multiplication of temperature not enabled, use addition"
             call exit(1)
          endif
       endif
    endif

  END SUBROUTINE read_variation

  subroutine read_forcing 

    use input_processor, only: getpari,getpars,getparl
    use unit_conversions, only: CtoK

    implicit none

    ! Local variables:
    integer :: ic, jc,j,  para, stn, year, day, month
    real(kind=sp) :: hour
    integer :: integration ! number of days to integrate the Englishman river
    integer, parameter :: Ieriver = 10000
    real(kind=sp) :: EriverI(Ieriver)
    logical ::  found_data
    real(kind=dp) :: wind_en, wind_nw, MajRiv, MinRiv ! temporary variables    
    ! startyears for the various variables that can be shifted
    integer :: wind_startyear, cf_startyear, rivers_startyear 
    integer :: rivers_startday ! startday for rivers allowing part year shifts

    if (river_n > Ieriver) then
       write (*,*) "in forcing.f90 need to increase size of EriverI array"
       stop
    endif
    
    ! read the number of days over which to integrate the Englishman river
    integration = getpari("Englishman integ days")
    
    ! read the start year in for initialization of runtime for model
    startyear = getpari("startyear")

    ! read the file names for input data (Wind, Met, Rivers)
    
    Wind = getpars("Wind")
    AT = getpars("Air temp")
    Cloud = getpars("Cloud")
    Hum = getpars("Humidity")
    MajorRiv = getpars("Major_River")
    isMinorRiv = getparl("isMinRiv")
    MinorRiv = getpars("Minor_River")
    

    ! WIND
    ! shift the start year if requested
    if (vary%wind%enabled) then
       wind_startyear = startyear + int(vary%wind%shift)
    else
       wind_startyear = startyear
    endif

    open(unit = 12, file = Wind, &
         status = "OLD", action = "READ")
    
    found_data = .false.
    do while (.not. found_data)
       read (12,*) day, month, year, hour, wind_en, wind_nw
       hour = hour  ! this is just to use "hour" to stop an unused warning
       month = month ! this is just to use "month" to stop an unused warning

       if (year == wind_startyear) then
          
          wind_eastnorth(1) = wind_en
          wind_northwest(1) = wind_nw

          do jc = 2, wind_n             
             read(12,*) day, month, year, hour, wind_eastnorth(jc), &
                  wind_northwest(jc)
          enddo          
          found_data = .true.
       endif
    enddo

    close (12)

    ! Cloud Fraction
    ! shift the start year if requested
    if (vary%cf%enabled .and. .not.(vary%cf%fixed) ) then
       cf_startyear = startyear + int(vary%cf%shift)
    else
       cf_startyear = startyear
    endif

    open(12, file= Cloud, &
         status = "OLD", action = "READ")
        
    found_data = .false.
    do while (.not. found_data)
       read (12,*) stn, year, month, day, para, (cf(1,j),j=1,24)
       para = para ! this is just to use "para" to stop an unused warning

       if (year == cf_startyear)then

          do jc = 2, met_n             
             read(12,*)stn,year,month,day,para,(cf(jc,j),j=1,24)
          enddo
          found_data = .true.
       endif
    enddo

    close (12)


    ! Air Temperature

    open (12, file= AT, &
         status = "OLD", action = "READ")
    
    found_data = .false.
    do while (.not. found_data)
       read (12,*) stn, year, month, day, para, (atemp(1,j),j=1,24)

       if (year == startyear)then

          do jc = 2, met_n
             read (12,*) stn, year, month, day, para,(atemp(jc,j),j=1,24)
             do j= 1, 24
                atemp(jc,j) = CtoK(atemp(jc,j) / 10.)
             enddo
          enddo
          found_data = .true.
       endif
    enddo

   

    close (5)

    close(12)

    ! Humidity

    open(12, file= Hum, &
         status = "OLD", action = "READ")
    
    found_data = .false.
    do while (.not. found_data)
       read (12,*) stn, year, month, day, para, (humid(1,j),j=1,24)

       if (year == startyear)then

          do jc = 2, met_n             
             read(12,*)stn,year,month,day,para,(humid(jc,j),j=1,24)
          enddo          
          found_data = .true.
       endif
    enddo

    close(12)

    ! Major River
    ! shift the start year if requested
    if (vary%rivers%enabled .and. .not. vary%rivers%fixed) then
       call shift_time (1, startyear, vary%rivers%shift, &
            rivers_startday, rivers_startyear)
    else
       rivers_startday = 1
       rivers_startyear = startyear
    endif



    open(12, file= MajorRiv, &
         status = "OLD", action = "READ")  

    found_data = .false.
    do while(.not. found_data)
       read(12,*) year, month, day, MajRiv

       if(year == rivers_startyear .and. day == rivers_startday) then
          
          Qriver(1) = MajRiv
          
          do jc = 2, river_n             
             read(12,*) year, month, day, Qriver(jc)
          enddo          
          found_data = .true.
          

       endif
    enddo

    close(12)

    if(isMinorRiv) then


       ! Minor  River
       ! uses same shift as Major
       open(12, file=MinorRiv, &
            status = "OLD", action = "READ")

       found_data = .false.
       do while (.not. found_data)
          read (12,*,END=175) year, month, day, MinRiv
       
          if (year == rivers_startyear .and. day == rivers_startday) then
          
             Eriver(1) = MinRiv

             do jc = 2, river_n             
                read (12,*) year, month, day, Eriver(jc)
             enddo
             found_data = .true.
          endif
       enddo


   !    FOR SOG, if there is not enough data in the englishman river, than use nanimo data.
175    if(.not. found_data) then
          open(12, file="../sog-forcing/rivers/NanimoNorm_historic.dat", &
               status = "OLD", action = "READ")
       
          do while (.not. found_data)
             read (12,*,END=175) year, month, day,MinRiv

             if (year == rivers_startyear .and. day == rivers_startday) then
          
                Eriver(1) = MinRiv

                do jc = 2, river_n             
                   read (12,*) year, month, day, Eriver(jc)
                enddo
                found_data = .true.
             endif
          enddo
       endif

       ! Want integrated Englishman River data over "integration" days
       do ic=1, integration
          EriverI(ic) = sum(Eriver(1:ic))/float(ic)
       enddo
       do ic=integration+1,river_n
          EriverI(ic) = sum(Eriver(ic-integration:ic))/integration
       enddo
       Eriver(1:river_n) = EriverI(1:river_n)
       close(12)  
    else 
       
       do jc = 1, river_n             
          Eriver(jc)=0
       enddo
  
    endif

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

    ! Wind data is in a single array with 1 hour per line, cf, atemp and
    ! hmid are in array of day and time, and the river flows are one value per day.

    ! local variable

    integer:: accul_day ! accumulated day number starting from startyear
    integer:: j ! index to hour, 1 am has j=2 
    integer:: index_day ! index into wind data for the end of the day before
    real(kind=dp) :: hr_frac ! fraction of the hour


    accul_day = accum_day(year,day)

    ! RIVERS

    if (vary%rivers%enabled .and. vary%rivers%fixed) then
       Qinter = SNGL(vary%rivers%value)
       Einter = 0.
    else
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
          cf_value = sngl(cf(accul_day,j)*vary%cf%fraction)
       endif
    else
       cf_value = cf(accul_day,j)
    endif

    ! TEMPERATURE

    if (vary%temperature%enabled .and. .not.vary%temperature%fixed) then
       ! multiply water temperature increment but 1.33 to get observed 
       ! Masson and Cummins 2004 air temperature increment
       atemp_value = atemp(accul_day,j) + sngl(vary%temperature%addition)*1.33
    else
       atemp_value = atemp(accul_day,j)
    endif

    humid_value = humid(accul_day,j)

    ! WIND

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
          quantity%addition = getpard(quantity_string//"%addition")
       endif
    endif

  END SUBROUTINE read_vary

  integer function accum_day (year, day) result(day_met)

    ! calculates the accumulated day from Jan 1,startyear.Based on the fact that 2 years of data is read in.

    INTEGER,INTENT(IN) :: day, year ! year day and year
    
    if (year == startyear) then
       day_met = day
    else
       if (leapyear(year-1)) then
          day_met = day + 366
       else
          day_met = day + 365
       endif
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
       call exit(1)
    endif
    
  end function leapyear
      
end module forcing
