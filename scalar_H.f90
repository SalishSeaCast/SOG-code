! $Id$
! $Source$

subroutine scalar_H(M, qty, Gvector, Gvector_o, Bmatrix_o, &
     Hvector)

  use precision_defs, only: dp
  use mean_param, only: trivector, prop
  use IMEX_constants, only: a_IMEX1

  implicit none

  integer, intent(in) :: M       ! Number of grid points
  type(prop), intent(in) :: qty  ! Model quantity, e.g. T, S
  real(kind=dp), dimension(M), intent(in) :: Gvector, Gvector_o
  type(trivector), intent(in) :: Bmatrix_o
  real(kind=dp), dimension(M), intent(inout) :: Hvector  ! Scalar vector

  integer :: index

  Hvector = 0.

  do index = 1, M
     Hvector(index) = qty%old(index) &
          + (1.0 - a_IMEX1) * (Gvector_o(index) &
          + Bmatrix_o%A(index) * qty%old(index-1) &
          + Bmatrix_o%B(index) * qty%old(index)                  &
          + Bmatrix_o%C(index) * qty%old(index+1)) &
          + a_IMEX1 * Gvector(index)
  end do
end subroutine scalar_H
