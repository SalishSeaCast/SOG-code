! $Id$
! $Source$

module tridiag_mod

contains

  subroutine tridiag(A, B, C, R, U) !!A,B,C,R,U are N dimensional vectors
    ! Solve the matrix equation MU = R, where M is a tridiagonal matrix
    ! with diagonal B, sub-diagonal A, and super-diagonal C.
    use precision_defs, only: dp
    use io_unit_defs, only: stderr
    implicit none
    ! Arguments:
    real(kind=dp), dimension(:), intent(in) :: &
         A, &  ! sub-diagonal vector of matrix M
         B, &  ! diagonal vector of matrix M
         C, &  ! super-diagonal vector of matrix M
         R     ! right-hand-side vector of matrix eq'm MU = R
    real(kind=dp), dimension(:), intent(out) :: U  ! solution of MU = R
    ! Local variables:
    real(kind=dp), dimension(size(R)) :: gamma
    real(kind=dp) :: beta
    integer :: j

    ! Confirm that the matrix equation is properly expressed
    if (B(1) == 0.) then
       write(stderr, *) "tridiag: Malformed matrix, B(1) = 0"
       stop
    endif
    ! Forward pass
    beta = B(1)
    U(1) = R(1) / beta
    do j = 2, size(R)
       gamma(j) = C(j-1) / beta
       beta = B(J) - A(J) * gamma(j)
       if (beta == 0.) then
          write(stderr, *) "tridiag: Solver failed, beta = 0 at j = ", j
          stop
       endif
       U(j) = (R(j) - A(j) * U(j-1)) / beta
    enddo
    ! Backward pass
    do j = size(R) - 1, 1, -1
       U(j) = U(j) - gamma(j+1) * U(j+1)
    enddo
  end subroutine tridiag

end module tridiag_mod
