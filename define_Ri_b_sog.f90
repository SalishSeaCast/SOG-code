! $Id$
! $Source$

subroutine define_Ri_b_sog(d, hh, surf_h, U_new, V_new, rho_g, Rib, &
     Vt_sq, N2)
  ! Calculate the profile of the bulk Richardson number, which is
  ! used to find the mixing layer depth.  See Large, etal (1994) pp
  ! 371-372, in particular, eq'n (21)
  ! *** Argument order needs to be standardized

  use precision_defs, only: dp
  use grid_mod, only: grid_, depth_average
  use fundamental_constants, only: g
  USE mean_param, only: height
  USE surface_forcing, only: ep

  ! Arguments:
  type(grid_), intent(in) :: d
  type(height), intent(in) :: hh 
  type(height), intent(out) :: surf_h
  real(kind=dp), dimension(0:d%M+1), intent(in) :: &
       U_new, &
       V_new, &
       rho_g
  real(kind=dp), dimension(1:d%M), intent(in) :: Vt_sq
  real(kind=dp), dimension(0:d%M), intent(out) :: Rib
  real(kind=dp), dimension(0:d%M+1), intent(in) :: N2
 
  ! Local variables:
  real(kind=dp), dimension(0:d%M+1) :: test_vector, test_vector2
  integer :: yy
  real(kind=dp) ::Ribmin, rho_avg, U_avg, V_avg

  surf_h%new = ep * hh%new
       
  U_avg = depth_average(U_new, 0.0d0, surf_h%new)
  V_avg = depth_average(V_new, 0.0d0, surf_h%new)
  rho_avg = depth_average(rho_g, 0.0d0, surf_h%new)

  test_vector = (U_avg - U_new) ** 2 + (V_avg - V_new) ** 2
  test_vector2 =  -g / rho_g(0) * (rho_avg - rho_g)
         
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
