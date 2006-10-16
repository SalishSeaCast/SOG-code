! $Id$
! $Source$

SUBROUTINE ND_flux_profile(grid, L_star, phi)
  ! Calculate the nondimensional flux profiles per Large, etal (1994),
  ! App. B
  ! *** Cleanup needed.
     
  use precision_defs, only: dp
  use grid_mod, only: grid_
  use turbulence, only: &
       xsi_m, &  ! Max xsi value of the -1/3 power law regime of phi%m
       xsi_s, &  ! Max xsi value of the -1/3 power law regime of phi%s
       a_m,   &  ! Coefficient of phi%m in 1/3 power law regime
       a_s,   &  ! Coefficient of phi%s in 1/3 power law regime
       c_m,   &  ! Coefficient of phi%m in 1/3 power law regime
       c_s       ! Coefficient of phi%s in 1/3 power law regime
  USE mean_param, only: MS
  implicit none
  ! Arguments:
  type(grid_), intent(in) :: grid
  real(kind=dp), intent(in) :: L_star
  type(MS), intent(out) :: phi
  ! Local variables:
  integer :: y

     !Define phi%m%value
      phi%m%value = 0.
      phi%s%value = 0.

      DO y = 0,grid%M
         IF (grid%d_i(y)/L_star >= 0.)THEN ! stable (defined as functions of d, the grid)
            phi%m%value(y) = 1.0+5.0*grid%d_i(y)/L_star ! (B1a)
         ELSE IF (grid%d_i(y)/L_star < 0. .AND. xsi_m <= grid%d_i(y)/L_star)THEN
            phi%m%value(y) = (1.0-16.0*grid%d_i(y)/L_star)**(-1.0/4.0) ! (B1b)
         ELSE
            phi%m%value(y) = (a_m - c_m*grid%d_i(y)/L_star)**(-1.0/3.0) ! (B1c)
         END IF

      !Define phi%s%value

         IF (grid%d_i(y)/L_star >= 0.)THEN
           phi%s%value(y) = 1.0+5.0*grid%d_i(y)/L_star ! (B1a)
         ELSE IF (grid%d_i(y)/L_star < 0. .AND. xsi_s <= grid%d_i(y)/L_star)THEN
           phi%s%value(y) = (1.0-16.0*grid%d_i(y)/L_star)**(-1.0/2.0) ! (B1d)
         ELSE
            phi%s%value(y) = (a_s - c_s*grid%d_i(y)/L_star)**(-1.0/3.0) ! (B1e)
         END IF
      END DO

END SUBROUTINE ND_flux_profile
 
