! $Id$
! $Source$

subroutine ML_height_sog(grid, Ri_b, new_h, kmax)
  ! Find the mixed layer depth.

  use precision_defs, only: dp
  use io_unit_defs, only: stderr
  use grid_mod, only: grid_
  use surface_forcing, only: Ri_c

  implicit none

  ! Arguments:
  type(grid_), intent(in) :: grid 
  real(kind=dp), dimension(grid%M), intent(in) :: Ri_b
  real(kind=dp), intent(out) :: new_h
  integer, intent(out) :: kmax
  ! Local variable:
  integer :: k

  ! Find the grid point at which the Richardson number exceeds the
  ! critical value
  kmax = grid%M             
  new_h =grid%d_g(kmax)    
  do k = 1, grid%M-2
     if (Ri_b(k) > Ri_c) then
        kmax = k 
        exit
     end if
  end do

  if (kmax == 1) then
     new_h = grid%d_g(kmax) * Ri_c / (Ri_b(kmax) + epsilon(Ri_b(kmax)))
  else if (kmax >= grid%M - 2) then
     ! Mixing deeper than bottom of grid
     kmax = grid%M-3
     new_h = grid%d_g(kmax) - grid%g_space(kmax) / 4.0 !#
     write(stderr, *) "ML_height_sog: Mixing too deep. ", &
          "Set kmax = ", kmax, " and new_h = ", new_h
  else
     ! Interpolate to find mixed layer depth
     new_h = grid%d_g(kmax) &
          + (Ri_c - Ri_b(kmax)) / (Ri_b(kmax) - Ri_b(kmax-1) &
          + epsilon(Ri_b(kmax) - Ri_b(kmax-1))) &
          * grid%g_space(kmax-1)
  end if
end subroutine ML_height_sog
