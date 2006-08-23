SUBROUTINE diffusion(mm,Bx,Gx,Kall,gamma,wso,J,step,scalar_end)

      USE mean_param
      USE surface_forcing

      IMPLICIT NONE

      TYPE(gr_d),INTENT(IN)::mm
      TYPE(trivector),INTENT(OUT)::Bx
      DOUBLE PRECISION,DIMENSION(mm%M),INTENT(IN OUT)::Gx
      DOUBLE PRECISION,DIMENSION(0:mm%M),INTENT(IN)::Kall
      DOUBLE PRECISION,DIMENSION(0:mm%M),INTENT(IN)::gamma
      DOUBLE PRECISION,INTENT(IN)::wso
      DOUBLE PRECISION,DIMENSION(0:mm%M),INTENT(IN)::J !interfaces  Q_n for temp and -F_n for salinity
      DOUBLE PRECISION, INTENT(IN)::step !dt
      DOUBLE PRECISION, INTENT(IN)::scalar_end  !eg. T(M+1)

      INTEGER::index
      DOUBLE PRECISION, DIMENSION(mm%M)::O_minus, O_plus

! Calculate Omega plus and minus below (D9)

      DO index = 1, mm%M
         O_minus(index) = step/mm%i_space(index)*Kall(index-1)/mm%g_space(index-1)
         O_plus(index) = step/mm%i_space(index)*Kall(index)/mm%g_space(index)
      END DO

      Bx%A = 0.
      Bx%B = 0.
      BX%C = 0.
 
! put in three terms.(D9)

      Bx%B = -O_plus -O_minus
      Bx%C = O_plus
      Bx%A = O_minus
      Bx%C(mm%M) = 0.
      
! calculate rhs (the H vector) (D10)
      Gx = 0.
!PRINT*,'in diff',Kall(1),gamma(1),wso,J(1),J(0)
      Gx(1) = step/mm%i_space(1)*(Kall(1)*gamma(1)-wso +J(1)-J(0))
!PRINT*,'Gx',Gx(1)
!pause
      DO index = 2, mm%M-1
         Gx(index) = step/mm%i_space(index)*(Kall(index)*gamma(index)-Kall(index-1)*&
         gamma(index-1)+ J(index)-J(index-1))
      END DO

      Gx(mm%M) = step/mm%i_space(mm%M)*(J(mm%M)-J(mm%M-1) + scalar_end*Kall(mm%M)/mm%g_space(mm%M))      

END SUBROUTINE diffusion



















