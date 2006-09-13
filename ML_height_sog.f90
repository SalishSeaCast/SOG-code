! $Id$
! $Source$

subroutine ML_height_sog(grid, Ri_b, year, day, day_time, count, hnew, jmaxg)
  ! Find the mixing layer depth.  See Large, etal (1994) pp 371-372.

  use precision_defs, only: dp
  use io_unit_defs, only: stdout
  use grid_mod, only: grid_
  use surface_forcing, only: Ri_c

  implicit none

  ! Arguments:
  type(grid_), intent(in) :: grid 
  real(kind=dp), dimension(grid%M), intent(in) :: Ri_b  ! Bulk Richardson no.
  ! *** Temporary arguments for debugging
  integer, intent(in) :: year, day, count
  real(kind=dp), intent(in) :: day_time
  real(kind=dp), intent(out) :: hnew  ! Mixing layer depth [m]
  integer, intent(out) :: jmaxg       ! Index of grid point below mixing depth
  ! Local variable:
  integer :: j  ! Index over grid depth

  ! Find the grid point at which the Richardson number exceeds the
  ! critical value
  jmaxg = grid%M             
  do j = 1, grid%M - 2
     if (Ri_b(j) > Ri_c) then
        jmaxg = j 
        exit
     end if
  end do

  if (jmaxg == 1) then
     ! Mixing layer is < 1 grid layer deep
     hnew = grid%d_g(1) * Ri_c / (Ri_b(1) + epsilon(Ri_b(1)))
  else if (jmaxg >= grid%M - 2) then
     ! Mixing layer extends nearly to bottom of grid
     jmaxg = grid%M - 3
     hnew = grid%d_g(jmaxg) - grid%g_space(jmaxg) / 4.0 !#
     write(stdout, *) "ML_height_sog: Mixing too deep. ", &
          "Set jmaxg = ", jmaxg, " and h%new = ", hnew
     write(stdout, *) "Iteration count = ", count, " Time: yr = ", &
          year, " day = ", day, " day_time = ", day_time
  else
     ! Interpolate to find mixing layer depth
     hnew = grid%d_g(jmaxg) &
          + (Ri_c - Ri_b(jmaxg)) / (Ri_b(jmaxg) - Ri_b(jmaxg-1) &
          + epsilon(Ri_b(jmaxg) - Ri_b(jmaxg-1))) &
          * grid%g_space(jmaxg-1)
  end if
end subroutine ML_height_sog
