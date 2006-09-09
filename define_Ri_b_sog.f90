! $Id$
! $Source$

subroutine define_Ri_b_sog(d, hh, surf_h, Bb, Uu, Vv, dens, Rib, &
     Vt_sq, N2)
  ! Calculate the profile of the bulk Richardson number, which is
  ! used to find the mixing layer depth.  See Large, etal (1994) pp
  ! 371-372, in particular, eq'n (21)
  ! *** Argument order needs to be standardized

  USE mean_param
  USE surface_forcing

  ! Arguments:
  type(gr_d), intent(in) :: d
  type(height), intent(in) :: hh 
  type(height), intent(out) :: surf_h
  type(prop), intent(in out) :: Bb, Uu, Vv, dens
  double precision, dimension(d%M), intent(in) :: Vt_sq
  double precision, dimension(d%M), intent(out) :: Rib
  double precision, dimension(0:d%M+1), intent(in) :: N2
 
  ! Local variables:
  double precision, dimension(0:d%M+1) :: test_vector, test_vector2
  integer :: yy, ymin
  double precision ::Ribmin

  surf_h%new = ep*hh%new
       
  call average(d,Uu,surf_h) !test conv  
  call average(d,Vv,surf_h) !test conv 
  call average(d,dens,surf_h)

  test_vector = (Uu%avg - Uu%new)**2.0 + (Vv%avg - Vv%new)**2.0 
  test_vector2 =  -g / dens%new(0) * (dens%avg - dens%new) !Bb%avg - Bb%new  
         
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
               ! *** This line doesn't appear to have any purpose,
               ! *** ymin is never used
               ymin = yy
            endif
        end if
     end if
  end do
  
end subroutine define_Ri_b_sog
