module malloc
  ! Memory allocation module.  Supplied memory allocation and
  ! deallocation status check subroutines.
  ! *** May also be the place for allocation subroutines for
  ! *** quantities that can't logically be allocated anywhere else.
  !
  ! Public Functions:
  !
  ! alloc_check(stat, msg) -- Check the stat result of an allocate
  !                           statement.  If it is non-zero, output 
  !                           the specified message on stdout and stop 
  !                           execution.
  !
  ! dalloc_check(stat, msg) -- Check the stat result of an deallocate
  !                            statement.  If it is non-zero, output 
  !                            the specified message on stdout and stop 
  !                            execution.

  implicit none

  private
  public :: alloc_check, dalloc_check

contains

  subroutine alloc_check(alloc_stat, message)
    ! Check the stat result of an allocate statement.
    ! If it is non-zero, output the specified message on stdout
    ! and stop execution.
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    integer, intent(in)           :: alloc_stat
    character(len=80), intent(in) :: message

    if (alloc_stat /= 0) then
       write(stdout, *) trim(message) // " memory allocation failed."
       call exit(1)
    end if
  end subroutine alloc_check


  subroutine dalloc_check(dalloc_stat, message)
    ! Check the stat result of an deallocate statement.
    ! If it is non-zero, output the specified message on stdout
    ! and stop execution.
    use io_unit_defs, only: stdout
    implicit none
    ! Arguments:
    integer, intent(in)           :: dalloc_stat
    character(len=80), intent(in) :: message

    if (dalloc_stat /= 0) then
       write(stdout, *) trim(message) // " memory deallocation failed."
       call exit(1)
    end if
  end subroutine dalloc_check

end module malloc
