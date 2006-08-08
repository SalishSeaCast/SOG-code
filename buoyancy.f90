SUBROUTINE buoyancy(alp, Te, Sa,d, hh, Bu, II, Br, rho, &
                    Cpi, beta, Qn)     

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::d
      TYPE(constant), INTENT(IN)::alp, beta, Cpi
      TYPE(height), INTENT(IN)::hh
      DOUBLE PRECISION, DIMENSION(0:d%M+1), INTENT(IN)::Te, Sa !T, S
      DOUBLE PRECISION, DIMENSION(0:d%M+1), INTENT(OUT)::Bu !B
      DOUBLE PRECISION, DIMENSION(0:d%M), INTENT(IN)::II 
      DOUBLE PRECISION, DIMENSION(0:d%M), INTENT(OUT)::Qn 
      DOUBLE PRECISION, INTENT(OUT)::Br !used in forcing
      DOUBLE PRECISION, DIMENSION(0:d%M+1), INTENT(IN)::rho  !density
                                                                       

      INTEGER::k
      DOUBLE PRECISION:: rho_hb, Cp_hb, alp_hb, In_hb
      TYPE(height)::h_B

      h_B%new =   hh%new  !(0 if no solar heat contributes to surface buoyancy flux)
      h_B%g =  hh%g
      h_B%i = hh%i

      Bu = g*(alp%g*Te - beta%g*Sa)    !Large eq.A3a
 
!!! Intensity should be defined on the interface levels and not the grid levels!!!
 

     Qn(0) = II(0)/Cpi%i(0)/rho(0)        
     DO k = 1, d%M       
         Qn(k) = II(k)/Cpi%i(k)/(rho(k)+(d%d_i(k)-d%d_g(k))*(rho(k+1)-rho(k))/d%g_space(k))       
     END DO       

      IF (h_B%g == 1) THEN   !!!Shouldn't make a difference!!!
         rho_hb = rho(0)
         Cp_hb = Cpi%i(0)
         alp_hb = alp%i(0)
         In_hb = II(0)
      ELSE
         rho_hb = rho(h_B%g-1)+(rho(h_B%g)-rho(h_B%g-1))*(h_B%new-d%d_g(h_B%g-1))/&
                 d%g_space(h_B%g-1)
         Cp_hb = Cpi%g(h_B%g-1)+(Cpi%g(h_B%g)-Cpi%g(h_B%g-1))*(h_B%new-d%d_g(h_B%g-1))/&
                 d%g_space(h_B%g-1)
         alp_hb = alp%g(h_B%g-1)+(alp%g(h_B%g)-alp%g(h_B%g-1))*(h_B%new-d%d_g(h_B%g-1))/&
                 d%g_space(h_B%g-1)
         In_hb = II(h_B%i-1) + (h_B%new-d%d_i(h_B%i-1))*(II(h_B%i)-II(h_B%i-1))/&
                 d%i_space(h_B%i)
      END IF

      Br = g*(alp%i(0)*II(0)/(Cpi%i(0)*rho(0))- &         !!eq. A3c in Large
                        alp_hb*In_hb/(rho_hb*Cp_hb))

!      write (*,*) Br/g/alp%i(0)*Cpi%i(0)*rho(0)

!      Br = 0.
!      Qn = 0.



END SUBROUTINE buoyancy
     

