! $Id$
! $Source$

module datetime
  ! A collection of type definitions, functions, and subroutines for
  ! manipulating dates and times.  Some are copied from the date_sub
  ! module (http://users.bigpond.net.au/amiller/datesub.f90 or
  ! http://ftp.cac.psu.edu/pub/ger/fortran/hdk/datesub.f90 or Google
  ! "date_sub module f90") compiled by H.D. Knoble in Jan-1972, and
  ! converted to f90 by Alan Miller on 1999-12-22.  Original
  ! references are cited in each of those routines.
  !
  ! Public Type Definitions:
  !
  !   datetime_ -- Date/time including year-day (1-Jan = 1), and
  !                day-second (midnight = 0).
  !
  !   timedelta -- Difference between 2 datetime_ structures.
  !
  ! Public Functions:
  !
  !   day_sec -- Return the day-sec value (seconds since midnight) for
  !              the specified datetime structure.
  !
  !   leapday -- Return and integer (0 or 1) for the number of leap
  !              days in the year of the specified datetime structure.
  !
  !   datetime_str -- Return the datetime value as a string.
  !
  !   year_day -- Given a datetime struct return the year-day.
  !
  !   datetime_incr -- Return the datetime incremented by the
  !                    specified timedelta.  Incremental timedelta
  !                    must be positive.
  !
  !   datetime_diff -- Return the difference between the datetime
  !                    structures as a timedelta strcuture.
  !
  ! Public Subroutines:
  !
  !   calendar_date -- Set the month, and day elements of the datetime
  !                    structure provided from its year and year-day
  !                    elements.
  !
  !   clock_time -- Set the hour, minute, and second elements of the
  !                 datetime structure provided from its day_sec
  !                 element.
  !
  !   os_datetime -- Read the operating system date and time and fill
  !                  in the yr, mo, day, hr, min & sec elements of the
  !                  datetime structure provided.

  implicit none

  private
  public :: &
       ! Type Definitions:
       datetime_, &  ! Date/time including year-day (1-Jan = 1), and
                     ! day-second (midnight = 0).
       timedelta, &  ! Difference between 2 datetime_ structures.
       ! Functions:
       day_sec,       &  ! Return the day-sec value (seconds since
                         ! midnight) for the specified datetime
                         ! structure.
       leapday,       &  ! Return and integer (0 or 1) for the number
                         ! of leap days in the year of the specified
                         ! datetime structure.
       datetime_str,  &  ! Return the datetime value as a string.
       year_day,      &  ! Given a datetime struct return the year-day.
       datetime_incr, &  ! Return the datetime incremented by the
                         ! specified timedelta.  Incremental timedelta
                         ! must be positive.
       datetime_diff, &  ! Return the difference between the datetime
                         ! structures as a timedelta strcuture.
       ! Subroutines:
       calendar_date, &  ! Set the month, and day elements of the
                         ! datetime structure provided from its year
                         ! and year-day elements.
       clock_time,    &  ! Set the hour, minute, and second elements
                         ! of the datetime structure provided from its
                         ! day_sec element.
       os_datetime       ! Read the operating system date and time and
                         ! fill in the yr, mo, day, hr, min & sec
                         ! elements of the datetime structure
                         ! provided.

  ! Public Type Definitions:
  !
  ! Date/time including year-day (1-Jan = 1), and day-second 
  ! (midnight = 0).
  type :: datetime_
     integer :: &
          yr,     &  ! Calendar year
          mo,     &  ! Calendar month (1 to 12)
          day,    &  ! Calendar day (1 to 31)
          yr_day, &  ! Year-day (1-Jan = 1; 1 to 365 or 366)
          hr,     &  ! Clock hour (0 to 23)
          min,    &  ! Clock minute (0 to 59)
          sec,    &  ! Clock second (0 to 59)
          day_sec    ! Day-second (midnight = 0; 0 to 86399)
  end type datetime_
  ! Difference between 2 datetime_ structures
  type :: timedelta
     integer :: &
          days, &  ! Number of days
          secs     ! Number of seconds (0 to 86399)
  end type timedelta

contains

  subroutine calendar_date(datetime)
    ! Set the month, and day elements of the datetime
    ! structure provided from its year and year-day elements.
    ! Based on calend subroutine from date_sub module.
    ! See ACM algorithm 398, Tableless Date Conversion, by Dick Stone,
    ! CACM 13(10):621.
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    type(datetime_), intent(inout) :: datetime

    ! Local variables:
    integer :: &
         yyyy, &
         ddd,  &
         mm,   &
         dd,   &
         t

    yyyy = datetime%yr
    ddd = datetime%yr_day
    t = leapday(yyyy)

    dd = ddd
    if (ddd > 59 + t) then
       dd = dd + 2 - t
    endif
    mm = ((dd + 91) * 100) / 3055
    dd = (dd + 91) - (mm * 3055) / 100
    mm = mm - 2

    datetime%mo = mm
    datetime%day = dd
    ! mm will be correct if, and only if, ddd is correct for yyyy
    if(mm < 1 .or. mm > 12) then
       write(stdout, *) "calendar_date(): Year-day out of range:", ddd
       call exit(1)
    endif
  end subroutine calendar_date


  subroutine clock_time(datetime)
    ! Set the hour, minute, and second elements of the datetime
    ! structure provided from its day_sec element.
    implicit none
    ! Argument:
    type(datetime_), intent(inout) :: datetime
    ! Local variable:
    integer :: rem

    ! Make use of integer division, and the mod() (remainder) function
    datetime%hr = datetime%day_sec / 3600
    rem = mod(datetime%day_sec, 3600)
    datetime%min = rem / 60
    datetime%sec = mod(rem, 60)
    return
  end subroutine clock_time


  function day_sec(datetime) result(day_s)
    ! Return the day-sec value (seconds since midnight) for the
    ! specified datetime structure.
    implicit none
    ! Argument:
    type(datetime_) :: datetime
    ! Result:
    integer :: day_s

    day_s = (datetime%hr * 3600) + (datetime%min * 60) + datetime%sec
  end function day_sec


  function leapday(year) result(leap_day)
    ! Return and integer (0 or 1) for the number of leap days in the
    ! year of the specified datetime structure.
    implicit none
    ! Argument:
    integer :: year
    ! Result:
    integer :: leap_day
    
    ! Most years aren't leapyears
    leap_day = 0
    ! Years that are evenly divisible by 4 are
    if(mod(year, 4) == 0) leap_day = 1
    ! But century years aren't unless they are evenly divisible by 400
    if(mod(year, 400) /= 0 &
         .and. mod(year, 100) == 0) leap_day = 0
  end function leapday


  function datetime_str(datetime, date_sep, datetime_sep, time_sep) result(str)
    ! Return the datetime value as a string.  Order of components is
    ! yr mo day, hr min s.  Date components are separated by the
    ! optional date_sep (defaults to '-'.  Date and time parts are
    ! separated by the optional datetime_sep (defaults to a space).
    ! Time components are separated by the optional time_sep
    ! (defaults to ':').  To get an ISO formatted date/time string use
    ! datetime_str(datetime, datetime_sep='T').  Use 'q' to get an
    ! empty separator, and no, there's no way to get q as a separator,
    ! because text processing in Fortran really sucks... :-)
    implicit none
    ! Arguments:
    type(datetime_), intent(in) :: datetime
    character(len=1), intent(in), optional :: date_sep, datetime_sep, time_sep
    ! Result:
    character(len=19) :: str
    ! Local variable:
    character(len=1)  :: dsep, dtsep, tsep
    character(len=10) :: date_str
    character(len=8)  :: time_str

    ! Establish what the separator characters are
    if (present(date_sep)) then
       dsep = date_sep
    else
       dsep = '-'
    endif
    if (present(datetime_sep)) then
       dtsep = datetime_sep
    else
       dtsep = ' '
    endif
    if (present(time_sep)) then
       tsep = time_sep
    else
       tsep = ':'
    endif
    ! Write the date/time to a string
    if (dsep /= 'q') then
       write(date_str, 100) datetime%yr, dsep, datetime%mo, dsep, datetime%day
100    format(i4, 2(a1,i2.2))
    else
       write(date_str, 101) datetime%yr, datetime%mo, datetime%day
101    format(i4, 2i2.2)
    endif
    if (tsep /= 'q') then
       write(time_str, 102) datetime%hr, tsep, datetime%min, tsep, datetime%sec
102    format(i2.2, 2(a1, i2.2))
    else
       write(time_str, 103) datetime%hr, datetime%min, datetime%sec
103    format(3i2.2)
    endif
    if (dtsep /= 'q') then
       str = trim(date_str) // dtsep // trim(time_str)
    else
       str = trim(date_str) // trim(time_str)
    endif
  end function datetime_str


  subroutine os_datetime(datetime)
    ! Read the operating system date and time and fill in the yr, mo,
    ! day, hr, min & sec elements of the datetime structure provided.
    implicit none
    ! Argument:
    type(datetime_), intent(out) :: datetime
    ! Local variable:
    integer, dimension(8) :: timeparts

    ! date_and_time is an intrinsic
    call date_and_time(values=timeparts)
    datetime%yr = timeparts(1)
    datetime%mo = timeparts(2)
    datetime%day = timeparts(3)
    ! timeparts(4) is "time difference in minutes" (whatever that means)
    datetime%hr = timeparts(5)
    datetime%min = timeparts(6)
    datetime%sec = timeparts(7)
    ! timeparts(8) is milliseconds

    ! Set the other 2 elements of datetime to their correct values
    datetime%yr_day = year_day(datetime)
    datetime%day_sec = day_sec(datetime)
  end subroutine os_datetime


  function year_day(datetime) result(yr_day)
    ! Given a datetime struct return the year-day.
    !
    ! Example: year_day(1984, 4, 22) = 113
    !
    ! Based on code from the date_sub module compiled by H.D. Knoble
    ! in Jan-1972, and converted to f90 by Alan Miller on 1999-12-22.
    implicit none
    ! Argument:
    type(datetime_), intent(in) :: datetime
    ! Result:
    integer :: yr_day
    ! Local variables:
    integer :: yyyy, mm, dd

    yyyy = datetime%yr
    mm = datetime%mo
    dd = datetime%day
    yr_day = 3055 * (mm + 2) / 100 - (mm + 10) / 13 * 2 - 91         &
         + (1 - (mod(yyyy, 4) + 3) / 4 + (mod(yyyy, 100) + 99) / 100 &
         - (mod(yyyy, 400) + 399) / 400) * (mm + 10) / 13 + dd
  end function year_day


  function datetime_incr(datetime, del) result(new_datetime)
    ! Return the datetime incremented by the specified timedelta.
    ! Incremental timedelta must be positive.
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    type(datetime_), intent(in) :: &
         datetime  ! Datetime to increment
    type(timedelta), intent(in) :: &
         del  ! Size of increment
    ! Result:
    type(datetime_) :: &
         new_datetime
    ! Local Variables:
    integer :: &
         days_in_yr  ! Number of days in the year

    ! Only positive increments allowed
    if (del%days < 0 .or. del%secs < 0) then
       write(stdout, *)  "datetime_incr(): Increment must be positive; ", &
            "days = ", del%days, " secs = ", del%secs
       call exit(1)
    end if
    ! Seconds increment must be less than a day
    if (del%secs > 86399) then
       write(stdout, *) "datetime_incr(): Increment seconds must be < 86400"
       call exit(1)
    end if
    ! Days increment must be less than a normal year
    !
    ! *** This restriction could be removed by figuring out how many
    ! *** leap days are included in an arbitrarily large days
    ! *** increment, but for present uses increments larger than 365
    ! *** days seem unlikely.
    if (del%days > 365) then
       write(stdout, *) "datetime_incr(): Increment days must be <= 365"
       call exit(1)
    end if
    ! Increment using day-seconds, year-day, and year, then calculate
    ! calendar date and clock time later.
    days_in_yr = 365 + leapday(datetime%yr)
    new_datetime%yr = datetime%yr
    new_datetime%yr_day = datetime%yr_day + del%days
    new_datetime%day_sec = datetime%day_sec + del%secs
    ! Handle roll-overs
    if (new_datetime%day_sec >= 86400) then
       ! Day-seconds
       new_datetime%day_sec = mod(new_datetime%day_sec, 86400)
       new_datetime%yr_day = datetime%yr_day + 1
    endif
    ! Year-days
    if (new_datetime%yr_day > days_in_yr) then
       new_datetime%yr_day = mod(new_datetime%yr_day, days_in_yr)
       ! Year
       new_datetime%yr = new_datetime%yr + 1
    endif
    ! Set the calendar date and clock time
    call calendar_date(new_datetime)
    call clock_time(new_datetime)
  end function datetime_incr


  function datetime_diff(dt_later, dt_earlier) result(del)
    ! Return the difference between the datetime structures as a
    ! timedelta strcuture.
    !
    ! *** For now we can only handle dt_later being after dt_earlier
    ! *** such that dt is positive.
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    type(datetime_), intent(in) :: &
         dt_later, &  ! Later in time datetime
         dt_earlier   ! Earlier in time datetime
    ! Result:
    type(timedelta) :: &
         del  ! Difference between specified datetimes
    ! Local Variables:
    integer :: &
         leap_days, &  ! Number of leap days between datetimes
         y             ! Index over years

    ! Years
    del%days = (dt_later%yr - dt_earlier%yr) * 365
    ! Leap years
    leap_days = 0
    ! Count the potential leap days in the interval
    do y = dt_earlier%yr, dt_later%yr
       leap_days = leap_days + leapday(y)
    end do
    ! Leap day in the last year of the interval are included in
    ! dt_later%yr_day.
    if (leapday(dt_later%yr) == 1) then
       leap_days = leap_days - 1
    endif
    del%days = del%days + leap_days
    ! Days
    del%days = del%days + (dt_later%yr_day - dt_earlier%yr_day)
    ! Seconds
    del%secs = dt_later%day_sec - dt_earlier%day_sec
    if (del%secs < 0) then
       del%secs = 86400 + del%secs
       del%days = del%days - 1
    endif
    ! This code has not been full tested for negative differences.  It
    ! certainly works for some cases, but until we're sure it's
    ! robust, let's not push our luck.
    if (del%days < 0) then
       write(stdout, *) "datetime_diff(): Only positive differences allowed"
       write(stdout, *) "  dt_later   = ", datetime_str(dt_later)
       write(stdout, *) "  dt_earlier = ", datetime_str(dt_earlier)
       call exit(1)
    end if
  end function datetime_diff
    
end module datetime
