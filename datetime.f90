! $Id$
! $Source$

module datetime
  ! A collection of functions and subroutines for manipulating dates
  ! and times.  Some are copied from the date_sub module (*** need a
  ! ref) compiled by H.D. Knoble in Jan-1972, and converted to f90 by
  ! Alan Miller on 1999-12-22.  Original references are cited in each
  ! of those routines.
  !
  ! Contents:
  !  calendar -- Subroutine.  
  !  year_day -- Function.  Given a calendar date (yyyy, mm, and dd) 
  !              return the day of the year (1-Jan = 1).  Oceanographers 
  !              may call this the Julian day (but it's not).  This 
  !              function is called iday in the date_sub module.

  implicit none

  ! *** We may need an interface block here to handle overloading of
  ! *** functions that need to be defined for more than 1 type.

  type :: datetime_
     integer :: yr
     integer :: mo
     integer :: day
     integer :: yr_day
     integer :: hr
     integer :: min
     integer :: sec
     integer :: day_sec
  end type datetime_

contains

  subroutine calendar_date(datetime)
    ! Set the month, and day elements of the datetime
    ! structure provided from its year and year-day elements.
    ! Based on calend subroutine from date_sub module.
    ! See ACM algorithm 398, Tableless Date Conversion, by Dick Stone,
    ! CACM 13(10):621.
    implicit none
    ! Arguments:
    type(datetime_), intent(inout) :: datetime

    ! Local variables:
    integer :: yyyy
    integer :: ddd
    integer :: mm
    integer :: dd
    integer :: t

    yyyy = datetime%yr
    ddd = datetime%yr_day
    t = 0
    if(mod(yyyy, 4) == 0) t = 1

    ! This statement is necessary if yyyy is < 1900 or > 2100
    if(mod(yyyy, 400) /= 0 .and. mod(yyyy, 100) == 0) t = 0

    dd = ddd
    if(ddd > 59 + t) dd = dd + 2 - t
    mm = ((dd + 91) * 100) / 3055
    dd = (dd + 91) - (mm * 3055) / 100
    mm = mm - 2
    datetime%mo = mm
    datetime%day = dd
    ! mm will be correct if, and only if, ddd is correct for yyyy
    if(mm >= 1 .and. mm <= 12) return
    ! *** Change this to write on stderr
    write(*,1) ddd
1   format('Error: calendar(): day of the year input =', i11, &
         ' is out of range.')
    stop
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


  integer function day_sec(datetime) result(day_s)
    ! 
    implicit none
    ! Argument:
    type(datetime_) :: datetime

    day_s = (datetime%hr * 3600) + (datetime%min * 60) + datetime%sec
  end function day_sec


  function datetime_str(datetime, separator) result(str)
    ! Return the datetime value as a string formatted as in ISO
    ! date/time with a space separating date and time
    implicit none
    ! Arguments:
    type(datetime_), intent(in) :: datetime
    character(len=*), intent(in), optional :: separator
    ! Result:
    character(len=18 + len(separator)) :: str
    ! Local variable:
    character(len=len(separator)) :: sep

    ! Establish what characters will separate date and time
    if (present(separator)) then
       sep = separator
    else
       sep = ' '
    endif
    ! Write the date/time to a string
    write(str, 100) datetime%yr, datetime%mo, datetime%day, &
         sep, datetime%hr, datetime%min, datetime%sec
100 format(i4, 2('-', i2.2), a, i2.2, 2(':', i2.2))
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

    ! Set the other 2 elements of datetime to zero
    datetime%yr_day = year_day(datetime)
    datetime%day_sec = day_sec(datetime)
  end subroutine os_datetime


  integer function year_day(datetime) result(yr_day)
    ! Given a datetime struct return the year-day.
    ! Example: year_day(1984, 4, 22) = 113
    ! Based on code from the date_sub module compiled by H.D. Knoble
    ! in Jan-1972, and converted to f90 by Alan Miller on 1999-12-22.
    implicit none
    ! Argument:
    type(datetime_), intent(in) :: datetime
    ! Local variables:
    integer :: yyyy, mm, dd
    yyyy = datetime%yr
    mm = datetime%mo
    dd = datetime%day

    yr_day = 3055 * (mm + 2) / 100 - (mm + 10) / 13 * 2 - 91         &
         + (1 - (mod(yyyy, 4) + 3) / 4 + (mod(yyyy, 100) + 99) / 100 &
         - (mod(yyyy, 400) + 399) / 400) * (mm + 10) / 13 + dd
  end function year_day

end module datetime
