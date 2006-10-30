! $Id$
! $Source$

module numerics
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the numerical simulation aspects of the
  ! SOG code.
  !
  ! Public Type Definitions:
  !
  !   tridiag -- Tridiagonal matrix arrays.
  !
  ! Public Parameters:
  !
  !   spam -- Description
  !
  ! Public Variables:
  !
  !   eggs -- Description
  !
  ! Public Subroutines:
  !
  !   check_negative -- 

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Type definitions:
       tridiag, &
!!$       ! Parameter values:
!!$       spam, &
!!$       ! Variables:
!!$       eggs, &
       ! Subroutines:
       check_negative

  ! Type Definitions:
  !
  ! Public:
  !
  ! Tridiagnonal matrix vectors:
  type :: tridiag
     real(kind=dp), dimension(:), allocatable :: &
          sub,  &  ! Sub-diagonal vector of a tridiagonal matrix
          diag, &  ! Diagonal vector of a tridiagonal matrix
          sup      ! Super-diagonal vector of a tridiagonal matrix
  end type tridiag
  !
  ! Private to module:

  ! Parameter Value Declarations:
  !
  ! Public:
  !
  ! Private to module:

contains

  subroutine check_negative(lbound, vector, msg, day, time, fatal)
    ! Check for negative values in a vector.  Output a message if a
    ! negative value is found, and optionally stop execution; default
    ! is to stop.
    use precision_defs, only: dp
    use io_unit_defs, only: stderr, stdout
    implicit none
    ! Arguments:
    integer, intent(in) :: &
         lbound  ! Lower bound of vector
    real(kind=dp), dimension(lbound:), intent(inout) :: &
         vector  ! Vector of values to check; inout so -ves can be set to zero
    character(len=*) :: &
         msg  ! Vector specific part of message to print
    integer, intent(in) :: &
         day  ! Year-day of current time step
    real(kind=dp), intent(in) :: &
         time  ! Time of current time step [s since start of run]
    logical, intent(in), optional :: &
         fatal  ! Is a negative value a fatal error? [Default = .true.]
    ! Local variable
    logical :: &
         die    ! Stop execution if a -ve value is found?
    integer :: &
         stream, &  ! Output stream to send message to; stderr or stdout
         j          ! Indices of most negative values in vector

    ! Handle optional fatal argument
    if (present(fatal)) then
       die = fatal
    else
       die = .true.  ! Default
    endif
    ! Check for negative value(s)
    if (any(vector < 0.0d0)) then
       ! Build message depending on whether or not -ve values are fatal
       if (die) then
          stream = stderr
          write(stderr, *) "Error: Negative value(s) in " // msg, &
               " at day = ", day, " time = ", time
       else
          stream = stdout
          write(stdout, *) "Warning: Negative value(s) in " // msg &
               // " set to zero at day = ", day, " time = ", time
       endif
       ! Find the index of the most -ve value
       do j = lbound, size(vector) + lbound - 1
          if (vector(j) < 0.0d0) then
             ! Output the message
             write(stream, *) " j = ", j
          endif
       enddo
       ! Stop execution if -ve values are fatal
       if (die) then
          stop
       else
          ! If not a fatal errror, set the -ve value to zero
          where (vector < 0.0d0) vector = 0.
       endif
    endif
  end subroutine check_negative

end module numerics
       
