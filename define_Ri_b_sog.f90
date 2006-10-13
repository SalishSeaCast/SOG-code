! $Id$
! $Source$

subroutine define_Ri_b_sog(d, hh, surf_h, Uu, Vv, rho, Rib, &
     Vt_sq, N2)
  ! Calculate the profile of the bulk Richardson number, which is
  ! used to find the mixing layer depth.  See Large, etal (1994) pp
  ! 371-372, in particular, eq'n (21)
  ! *** Argument order needs to be standardized

  use precision_defs, only: dp
  use grid_mod, only: grid_, depth_average
  USE mean_param, only: height, prop
  USE surface_forcing, only: g, ep

  ! Arguments:
  type(grid_), intent(in) :: d
  type(height), intent(in) :: hh 
  type(height), intent(out) :: surf_h
  real(kind=dp), dimension(0:d%M+1), intent(in) :: rho
  type(prop), intent(in out) :: Uu, Vv
  real(kind=dp), dimension(1:d%M), intent(in) :: Vt_sq
  real(kind=dp), dimension(0:d%M), intent(out) :: Rib
  real(kind=dp), dimension(0:d%M+1), intent(in) :: N2
 
  ! Local variables:
  real(kind=dp), dimension(0:d%M+1) :: test_vector, test_vector2
  integer :: yy
  real(kind=dp) ::Ribmin, rho_avg, U_avg, V_avg

  surf_h%new = ep * hh%new
       
  U_avg = depth_average(Uu%new, 0.0d0, surf_h%new)
  V_avg = depth_average(Vv%new, 0.0d0, surf_h%new)
  rho_avg = depth_average(rho, 0.0d0, surf_h%new)
!!$print *, ''
!!$print *, "d%d_g(0:3) ", d%d_g(0:3)
!!$print *, "rho(0:3) ", rho(0:3)
!!$print *, ''
!!$rho_avg = depth_average(d%d_g, 0.0d0, 1.25d0)
!!$print *, "d%d_g_avg(0, 1.25) ", rho_avg
!!$print *, ''
!!$rho_avg = depth_average(rho, 0.0d0, 1.25d0)
!!$print *, "rho_avg(0, 1.25) ", rho_avg
!!$print *, ''
!!$rho_avg = depth_average(d%d_g, 0.75d0, 0.75d0)
!!$print *, "d%d_g_avg(0.75, 0.75) ", rho_avg
!!$print *, ''
!!$rho_avg = depth_average(rho, 0.75d0, 0.75d0)
!!$print *, "rho_avg(0.75, 0.75) ", rho_avg
!!$print *, ''
!!$rho_avg = depth_average(d%d_g, 0.0d0, 0.25d0)
!!$print *, "d%d_g_avg(0, 0.25) ", rho_avg
!!$print *, ''
!!$rho_avg = depth_average(rho, 0.0d0, 0.25d0)
!!$print *, "rho_avg(0, 0.25) ", rho_avg
!!$print *, ''
!!$rho_avg = depth_average(rho, 0.0d0, 0.1d0)
!!$print *, "rho_avg(0, 0.1) ", rho_avg
!!$stop

  test_vector = (U_avg - Uu%new) ** 2 + (V_avg - Vv%new) ** 2
  test_vector2 =  -g / rho(0) * (rho_avg - rho)
         
  Rib = 0.0
  Ribmin = 1000.
           
  do yy = 1, d%M !surf_h%g, d%M 
     if (N2(yy) >= 0.) then
        if (test_vector2(yy) > EPSILON(test_vector2(yy)) .OR. &
             test_vector2(yy) < -EPSILON(test_vector2(yy))) then
           Rib(yy) = test_vector2(yy) * d%d_g(yy) &
                / (test_vector(yy) + Vt_sq(yy) + 1.0D-30)
           if (Rib(yy) < Ribmin) then
               Ribmin = Rib(yy)
            endif
        end if
     end if
  end do
  
end subroutine define_Ri_b_sog
