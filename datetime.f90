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
    write(*, 1) ddd
1   format('Error: calendar(): day of the year input =', i11, &
         ' is out of range.')
    call exit(1)
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
