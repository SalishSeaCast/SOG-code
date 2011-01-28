module test_datetime
  use fruit
contains
  subroutine test_increment_date_time_non_leapyear
    ! new_year increments end of non-leapyear correctly
    use precision_defs, only: dp
    use datetime, only: increment_date_time

    real (kind=dp) :: &
         day_time, &  ! Day-sec counter
         time         ! Time counter through run; seconds since midnight
                      ! of initDatetime
    integer :: &
         day, &    ! Year-day counter
         year, yr  ! Year counter
    real (kind=dp) :: &
         dt  ! Time step [s]

    day_time = 86350.0d0
    day = 365
    year = 1999
    time = 0.0d0
    call increment_date_time(day_time, day, year, time, 100.0d0)
    call assertEquals(2000, year)
    call assertEquals(1, day)
  end subroutine test_increment_date_time_non_leapyear


  subroutine test_increment_date_time_leapyear_2000
    ! new_year increments end of year 2000 correctly
    use precision_defs, only: dp
    use datetime, only: increment_date_time

    real (kind=dp) :: &
         day_time, &  ! Day-sec counter
         time         ! Time counter through run; seconds since midnight
                      ! of initDatetime
    integer :: &
         day, &    ! Year-day counter
         year, yr  ! Year counter
    real (kind=dp) :: &
         dt  ! Time step [s]

    day_time = 86350.0d0
    day = 365
    year = 2000
    time = 0.0d0
    call increment_date_time(day_time, day, year, time, 100.0d0)
    call assertEquals(2000, year)
    call assertEquals(366, day)
    day_time = 86350.0d0
    call increment_date_time(day_time, day, year, time, 100.0d0)
    call assertEquals(2001, year)
    call assertEquals(1, day)
  end subroutine  test_increment_date_time_leapyear_2000
end module test_datetime
