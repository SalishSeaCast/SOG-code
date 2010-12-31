module input_processor
  ! Subroutines and functions for processing the SOG infile (read from
  ! stdin).  This module provides the utility routines, but the actual
  ! reading of parameter values is done by the read_* subroutines of
  ! various modules, using these routines.
  !
  ! Public Subroutines:
  !
  ! init_input_processor(runDatetime, [showReport])
  !    -- Initialize the input processor.
  !
  ! Public Functions:
  !
  ! getpars(name, [showReport]) -- Read the named string value.
  !
  ! getpari(name, [showReport]) -- Read the named single integer value.
  !
  ! getpard(name, [showReport]) -- Read the named single real(kind=dp) value.
  !
  ! getparl(name, [showReport]) -- Read the named logical value.
  !
  ! getpariv(name, [showReport]) -- Read the named vector of integer values.
  !
  ! getpardv(name, [showReport]) -- Read the named vector of
  !                                 real(kind=dp) values.

  implicit none

  private
  public :: &
       ! Subroutines:
       init_input_processor, &
       ! Functions:
       getpars, getpari, getpard, getparl, getpariv, getpardv, &
       getpar_datetime

  ! Private variable declarations:
  logical :: show_report  ! getpar* routines report to stdout when TRUE
  integer, parameter :: &
       label_len = 25,  &  ! Max length of infile item label
       desc_len = 50,   &  ! Max length of infile item description
       str_len = 80        ! Max length of string value

contains

  subroutine init_input_processor(runDatetime, showReport)
    ! Initialize the input processor by:
    !   -Setting the state of the getpar* report to stdout
    !   -Starting the report, if requested
    !   -Read the infile from stdin, stripping comments from it, and
    !    writing the non-comment lines to a scratch file, ready for 
    !    the read_* subroutines of other modules the get values from.
    use io_unit_defs, only: stdout, stripped_infile
    implicit none
    ! Arguments:
    character(len=19), intent(in) :: runDatetime   ! Date/time of code run
    logical, intent(in), optional :: showReport    ! Report to stdout if TRUE
    ! Local variable:
    integer, parameter :: max_line_len = 240
    character(len=max_line_len) :: line  ! A line from the infile

    ! show_report defaults to TRUE
    if (present(showReport)) then
       show_report = showReport
    else
       show_report = .true.
    endif
    ! Start the report, if requested
    if (show_report) then
       write(stdout, '(a50, " = ", a)') "Date", runDatetime
    end if

    ! Open a scratch file to write the non-comment lines from the
    ! infile to
    open(stripped_infile, status='scratch', action='readwrite')
    ! Read the infile and write the non-comment lines to the scratch file
    do
       read(*, '(a240)', end=100) line
       if (line(1:1) /= '!') then
          write(stripped_infile, *) trim(line)
       endif
    enddo
    ! Rewind the scratch file to make it ready to be read by the
    ! read_* subroutines
100 rewind(stripped_infile)
  end subroutine init_input_processor


  function getpars(name, reportThis) result(str)
    ! Read the named string value, and report it to stdout, if
    ! requested.
    use io_unit_defs, only: stdout, stripped_infile
    implicit none
    ! Arguments:
    character(len=*), intent(in)  :: name
    logical, intent(in), optional :: reportThis
    ! Result:
    character(len=str_len) :: str
    ! Local variables:
    logical :: report_this
    character(len=label_len) :: label
    character(len=desc_len) :: desc

    ! Value reporting to stdout defaults to TRUE
    if (present(reportThis)) then
       report_this = reportThis
    else
       report_this = .true.
    endif
    ! Read and validate the value
    read(stripped_infile, *) label, str, desc
    if (name /= label) then
       write(stdout, *) "getpars: Expecting ", name, " but got ", &
            trim(label), " instead."
       call exit(1)
    end if
    ! Trim trailing blanks from the value
    str = trim(str)
    ! Report the value to stdout, if requested
    if (show_report .and. report_this) then
       write(stdout,'(a50, " = ", a)') trim(desc), trim(str)
    endif
  end function getpars


  function getpari(name, reportThis) result(int)
    ! Read the named single integer value, and report it to stdout, if
    ! requested.
    use io_unit_defs, only: stdout, stripped_infile
    implicit none
    ! Arguments:
    character(len=*), intent(in)  :: name
    logical, intent(in), optional :: reportThis
    ! Result:
    integer :: int
    ! Local variables:
    logical :: report_this             ! Report value to stdout?
    character(len=label_len) :: label  ! Item label in infile
    character(len=desc_len) :: desc    ! Item description in infile
    character(len=16) :: str           ! String representation of value

    ! Value reporting to stdout defaults to TRUE
    if (present(reportThis)) then
       report_this = reportThis
    else
       report_this = .true.
    endif
    ! Read and validate the value
    read(stripped_infile, *) label, int, desc
    if (name /= label) then
       write(stdout, *) "getpari: Expecting ", name, " but got ", &
            trim(label), " instead."
       call exit(1)
    end if
    ! Report the value to stdout, if requested
    if (show_report .and. report_this) then
       ! Create a left justified string representation of the value
       write(str, '(i16)'), int
       str = trim(adjustl(str))
       write(stdout,'(a50, " = ", a)') trim(desc), str
    endif
  end function getpari


  subroutine getpariv(name, len, vec, reportThis)
    ! Read the named vector of integer values, and report it to stdout, if
    ! requested.
    use io_unit_defs, only: stdout, stripped_infile
    implicit none
    ! Arguments:
    character(len=*), intent(in)  :: name
    integer, intent(in) :: len
    integer, dimension(:), intent(out) :: vec
    logical, intent(in), optional :: reportThis
    ! Local variables:
    logical :: report_this             ! Report value to stdout?
    character(len=label_len) :: label  ! Item label in infile
    character(len=desc_len) :: desc    ! Item description in infile
    character(len=16) :: str           ! String representation of value
    integer :: i                       ! Loop index

    ! Value reporting to stdout defaults to TRUE
    if (present(reportThis)) then
       report_this = reportThis
    else
       report_this = .true.
    endif
    ! Read and validate the value
    read(stripped_infile, *) label, (vec(i), i = 1, len), desc
    if (name /= label) then
       write(stdout, *) "getpariv: Expecting ", name, " but got ", &
            trim(label), " instead."
       call exit(1)
    end if
    ! Report the value to stdout, if requested
    if (show_report .and. report_this) then
       ! Create a left justified string representation of the value
       write(str, '(i16)'), vec(1)
       str = trim(adjustl(str))
       write(stdout,'(a46, " [1] = ", a)') trim(desc), str
       do i = 2, len
          write(str, '(i16)'), vec(i)
          str = trim(adjustl(str))
          if (i <= 9) then
             write(stdout, '(46x, " [", i1, "] = ", a)') i, str
          else
             write(stdout, '(45x, " [", i2, "] = ", a)') i, str
          endif
       enddo
    endif
  end subroutine getpariv


  function getpard(name, reportThis) result(num)
    ! Read the named single real(kind=dp) value, and report it to
    ! stdout, if requested.
    use io_unit_defs, only: stdout, stripped_infile
    use precision_defs, only: dp
    implicit none
    ! Arguments:
    character(len=*), intent(in)  :: name
    logical, intent(in), optional :: reportThis
    ! Result:
    real(kind=dp) :: num
    ! Local variables:
    logical :: report_this             ! Report value to stdout?
    character(len=label_len) :: label  ! Item label in infile
    character(len=desc_len) :: desc    ! Item description in infile
    character(len=16) :: str           ! String representation of value

    ! Value reporting to stdout defaults to TRUE
    if (present(reportThis)) then
       report_this = reportThis
    else
       report_this = .true.
    endif
    ! Read and validate the value
    read(stripped_infile, *) label, num, desc
    if (name /= label) then
       write(stdout, *) "getpard: Expecting ", name, " but got ", &
            trim(label), " instead."
       call exit(1)
    end if
    ! Report the value to stdout, if requested
    if (show_report .and. report_this) then
       ! Create a left justified string representation of the value
       write(str, '(1pg16.8)'), num
       str = trim(adjustl(str))
       write(stdout,'(a50, " = ", a)') trim(desc), str
    endif
  end function getpard


  subroutine getpardv(name, len, vec, reportThis)
    ! Read the named vector of real(kind=dp) values, and report it to
    ! stdout, if requested.
    use io_unit_defs, only: stdout, stripped_infile
    use precision_defs, only: dp
    implicit none
    ! Arguments:
    character(len=*), intent(in)  :: name
    integer, intent(in) :: len
    real(kind=dp), dimension(:), intent(out) :: vec
    logical, intent(in), optional :: reportThis
    ! Local variables:
    logical :: report_this             ! Report value to stdout?
    character(len=label_len) :: label  ! Item label in infile
    character(len=desc_len) :: desc    ! Item description in infile
    character(len=16) :: str           ! String representation of value
    integer :: i                       ! Loop index

    ! Value reporting to stdout defaults to TRUE
    if (present(reportThis)) then
       report_this = reportThis
    else
       report_this = .true.
    endif
    ! Read and validate the value
    read(stripped_infile, *) label, (vec(i), i = 1, len), desc
    if (name /= label) then
       write(stdout, *) "getpardv: Expecting ", name, " but got ", &
            trim(label), " instead."
       call exit(1)
    end if
    ! Report the value to stdout, if requested
    if (show_report .and. report_this) then
       ! Create a left justified string representation of the value
       write(str, '(1pg16.8)'), vec(1)
       str = trim(adjustl(str))
       write(stdout,'(a46, " [1] = ", a)') trim(desc), str
       do i = 2, len
          write(str, '(1pg16.8)'), vec(i)
          str = trim(adjustl(str))
          if (i <= 9) then
             write(stdout, '(46x, " [", i1, "] = ", a)') i, str
          else
             write(stdout, '(45x, " [", i2, "] = ", a)') i, str
          endif
       enddo
    endif
  end subroutine getpardv


  function getparl(name, reportThis) result(bool)
    ! Read the named logical value, and report it to
    ! stdout, if requested.
    use io_unit_defs, only: stdout, stripped_infile
    implicit none
    ! Arguments:
    character(len=*), intent(in)  :: name
    logical, intent(in), optional :: reportThis
    ! Result:
    logical :: bool
    ! Local variables:
    logical :: report_this             ! Report value to stdout?
    character(len=label_len) :: label  ! Item label in infile
    character(len=desc_len) :: desc    ! Item description in infile
    character(len=16) :: str           ! String representation of value

    ! Value reporting to stdout defaults to TRUE
    if (present(reportThis)) then
       report_this = reportThis
    else
       report_this = .true.
    endif
    ! Read and validate the value
    read(stripped_infile, *) label, bool, desc
    if (name /= label) then
       write(stdout, *) "getparl: Expecting ", name, " but got ", &
            trim(label), " instead."
       call exit(1)
    end if
    ! Report the value to stdout, if requested
    if (show_report .and. report_this) then
       ! Create a left justified string representation of the value
       write(str, '(l1)'), bool
       str = trim(adjustl(str))
       write(stdout,'(a50, " = ", a)') trim(desc), str
    endif
  end function getparl


  function getpar_datetime(name, reportThis) result(date_time)
    ! Read and validate the named datetime value, and report ti to
    ! stdout, if requested.

    ! Elements from other modules:
    ! Type Definitions:
    use datetime, only: datetime_
    ! Functions:
    use datetime, only: datetime_str, leapday, year_day, day_sec
    ! Parameter Values:
    use io_unit_defs, only: stdout, stripped_infile

    implicit none

    ! Arguments:
    character(len=*), intent(in)  :: name
    logical, intent(in), optional :: reportThis

    ! Result:
    type(datetime_) :: date_time

    ! Local variables:
    logical :: report_this             ! Report value to stdout?
    character(len=label_len) :: label  ! Item label in infile
    character(len=19) :: date_time_str ! Parameter value string in infile
    character(len=desc_len) :: desc    ! Item description in infile
    character(len=19) :: str           ! String representation of value
    character(len=1), dimension(5) :: sep  ! Datetime string separators

    ! Value reporting to stdout defaults to TRUE
    if (present(reportThis)) then
       report_this = reportThis
    else
       report_this = .true.
    endif
    ! Read and validate the value
    read(stripped_infile, *) label, date_time_str, desc
    ! Confirm we're reading the expected parameter
    if (name /= label) then
       write(stdout, *) "getpar_datetime: Expecting ", name, " but got ", &
            trim(label), " instead."
       call exit(1)
    end if
    ! Confirm that the string we read is the correct length
    if (len(trim(date_time_str)) /= 19) then
       write(stdout, *) "getpar_datetime: ", trim(label)
       write(stdout, *) "  Expecting yyyy-mm-dd hh:mm:ss but got ", &
            date_time_str
       call exit(1)
    endif
    ! Read the datetime elements from the parameter string value
    read(date_time_str, '(i4, 5(a1, i2))') &
         date_time%yr, sep(1), date_time%mo, sep(2), date_time%day, sep(3), &
         date_time%hr, sep(4), date_time%min, sep(5), date_time%sec
    ! Check the separators
    if (any(sep(1:2) /= '-') &
         .or. sep(3) /= ' ' &
         .or. any(sep(4:5) /= ':')) then
       write(stdout, *) "getpar_datetime: ", trim(label)
       write(stdout, *) "  Expecting yyyy-mm-dd hh:mm:ss but got ", &
            date_time_str
       call exit(1)
    endif
    ! Check the year (fairly lame, but better than nothing)
    if (date_time%yr < 1900 .or. date_time%yr > 2100) then
       write(stdout, *) "getpar_datetime:  ", trim(label)
       write(stdout, *) "  Year seems odd: ", date_time_str
       call exit(1)
    endif
    ! Check the month
    if (date_time%mo < 1 .or. date_time%mo > 12) then
       write(stdout, *) "getpar_datetime:  ", trim(label)
       write(stdout, *) "  Invalid month: ", date_time_str
       call exit(1)
    endif
    ! Check the day
    if (date_time%day < 1 &
         ! 31 day months
         .or. ((date_time%mo == 1 .or. date_time%mo == 3 &
               .or. date_time%mo == 5 .or. date_time%mo == 7 &
               .or. date_time%mo == 8 .or. date_time%mo == 10 &
               .or. date_time%mo == 12) .and. date_time%day > 31) &
         ! 30 day months
         .or. ((date_time%mo == 4 .or. date_time%mo == 6 &
               .or. date_time%mo == 9 .or. date_time%mo == 11) &
               .and. date_time%day > 30) &
         ! February
         .or. (date_time%mo == 2 &
              .and. ((leapday(date_time%yr) == 1 &
                      .and. date_time%day > 29) &
              .or. (leapday(date_time%yr) == 0 &
                    .and. date_time%day > 28)))) then
       write(stdout, *) "getpar_datetime:  ", trim(label)
       write(stdout, *) "  Invalid day: ", date_time_str
       call exit(1)
    endif
    ! Check the time
    if (date_time%hr < 0 .or. date_time%hr > 23 &
         .or. date_time%min < 0 .or. date_time%min > 59 &
         .or. date_time%sec < 0 .or. date_time%sec > 59) then
       write(stdout, *) "getpar_datetime:  ", trim(label)
       write(stdout, *) "  Invalid time: ", date_time_str
       call exit(1)
    endif
    ! Set the other 2 elements of datetime to their correct values
    date_time%yr_day = year_day(date_time)
    date_time%day_sec = day_sec(date_time)
    ! Report the value to stdout, if requested
    if (show_report .and. report_this) then
       ! Create a left justified string representation of the value
       str = datetime_str(date_time)
       write(stdout, 100) trim(desc), str, date_time%yr_day, date_time%day_sec
100    format(a50, " = ", a, " year-day = ", i3, " day-sec = ", i5) 
    endif
  end function getpar_datetime

end module input_processor
