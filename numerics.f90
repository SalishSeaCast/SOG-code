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
  !   ham -- Description

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Type definitions:
       tridiag
!!$       ! Parameter values:
!!$       spam, &
!!$       ! Variables:
!!$       eggs, &
!!$       ! Subroutines:
!!$       ham

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

!!$contains

end module numerics
       
