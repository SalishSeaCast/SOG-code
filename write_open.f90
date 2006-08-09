subroutine write_open 
! Get the names of 2 output files from stdin using getpars(), and open
! the files
  
  ! *** initial_sog supplies getpars(), but this should change
  use initial_sog

  implicit none

  ! File name to open
  character*80 :: str           

  ! *** File name tags could be more descriptive...
  ! Physics model time series results
  str = getpars("writefile1", 1)
  open (unit=293, file=str, status="replace", action="write")
  ! Biology model time sereis results
  str = getpars("writefile2", 1)
  open (unit=295, file=str, status="replace", action="write")

end subroutine write_open
