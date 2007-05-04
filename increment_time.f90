! $Id$
! $Source$

module increment_time

  use precision_defs, only: dp

  implicit none

  private
  public :: new_year

contains

  subroutine new_year(day_time, day, year, time, dt)

    implicit none

    ! Argument:
    real (kind=dp), INTENT(in out) :: day_time, time
    real (kind=dp), INTENT(in) :: dt            
    integer, INTENT(in out) :: day, year
    ! Local Parameters
    real (kind=dp), parameter :: nosecondperday = 86400.
    real (kind=dp), parameter :: nodayperyear = 365.
    
    day_time = day_time + dt
    time = time + dt
    
    ! have we crossed midnight?
    if (day_time >= nosecondperday) then
       day_time = day_time - nosecondperday
       day = day + 1
       
       ! have we crossed year end?
       if (day > nodayperyear) then 
          if ( .not. leapyear(year) .or. day > nodayperyear+1) then
             day = 1
             year = year + 1
          endif
       endif
    endif
    
  end subroutine new_year

  logical function leapyear (year)
    
    integer, intent (in) :: year
    integer, parameter :: leaps(15) = (/1952, 1956, 1960, 1964, 1966, 1972, 1976, 1980, 1984, 1988, 1992, 1996, 2000, 2004, 2008/)
  
    
    leapyear = .false.
    if (year < 1950 .or. year > 2010) then
       write (*,*) "Function leap-year out of bounds"
       stop
    endif
  if (any(leaps == year)) leapyear = .true.
  
end function leapyear
end module increment_time


