program fruit_driver
  use fruit
  use test_datetime

  call init_fruit
  call test_increment_date_time_non_leapyear
  call test_increment_date_time_leapyear_2000
  call fruit_summary
end program fruit_driver
