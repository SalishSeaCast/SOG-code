SUBROUTINE find_PON

  USE surface_forcing
  USE declarations

  IMPLICIT NONE

  DOUBLE PRECISION::oldavg,newavg,depth2,depth1
  DOUBLE PRECISION, DIMENSION(0:M+1)::newNO, newNH
  INTEGER:: isusan
  
  DO yy = 1, M
     IF (P1_p(yy) < zero) THEN
        P1_p(yy) = zero
     END IF
! Flagella.v.0
     IF (Pnano1_p(yy) < zero) THEN
        Pnano1_p(yy)= zero
     END IF

     IF (NO1_p(yy) < zero) THEN
        NO1_p(yy) = zero
     END IF
     IF (NH1_p(yy) < zero) THEN
        NH1_p(yy) = zero
     END IF
     DO xx = 1,D_bins
        IF (Detritus1_p(xx,yy) < zero) THEN
           Detritus1_p(xx,yy) = zero
        END IF
     END DO
     
  END DO
 
  DO xx = 1,M
     PON%new(xx) = P1_p(xx)+Pnano1_p(xx)+Detritus1_p(1,xx)+Detritus1_p(2,xx)
     !Flagella.v.0 PON%new(xx) = P1_p(xx)+Detritus1_p(1,xx)+Detritus1_p(2,xx)
     PON%old(xx) = P%nano%new(xx)+P%micro%new(xx)+Detritus(1)%D%new(xx)+Detritus(2)%D%new(xx)
     !Flagella.v.0 PON%old(xx) = P%micro%new(xx)+Detritus(1)%D%new(xx)+Detritus(2)%D%new(xx)
     newNO(xx) = NO1_p(xx)
     newNH(xx) = NH1_p(xx)
  END DO


  PON%new(0) = PON%new(1)
  newNO(0) = newNO(1)
  newNH(0) = newNH(1)
  PON%new(M+1) = 0.
  newNO(M+1) = N%O%new(M+1)
  newNH(M+1) = N%H%new(M+1)

  depth1 = 39.
  depth2 = 20.

  CALL average2(grid,PON%old(1:M),depth1,oldavg)  !PON%old is Pnano%new
  CALL average2(grid,PON%new(1:M),depth1,newavg)  !PON%new is Pnano1_p(xx)

  PONflux_200 = (oldavg-newavg)*depth1/dt*3600.*24.
  !IF (PONflux_200 < 0.) THEN
     !PRINT "(A)","bottom PON flux is negative??? PONflux_200,oldavg,newavg,depth1,dt"
     !PRINT *,PONflux_200,oldavg,newavg,depth1,dt
     !pause
  !END IF

  CALL average2(grid,N%O%new(1:M),depth1,oldavg)
  CALL average2(grid,newNO(1:M),depth1,newavg)

  NOflux_200 = (oldavg-newavg)*depth1/dt*3600.*24.

  CALL average2(grid,N%H%new(1:M),depth1,oldavg)
  CALL average2(grid,newNH(1:M),depth1,newavg)

  NHflux_200 = (oldavg-newavg)*depth1/dt*3600.*24.

  CALL average2(grid,PON%old(1:M),depth2,oldavg)
  CALL average2(grid,PON%new(1:M),depth2,newavg)

  PONflux_100 = (oldavg-newavg)*depth2/dt*3600.*24.
 
  CALL average2(grid,N%O%new(1:M),depth2,oldavg)
  CALL average2(grid,newNO(1:M),depth2,newavg)

  NOflux_100 = (oldavg-newavg)*depth2/dt*3600.*24.

  CALL average2(grid,N%H%new(1:M),depth2,oldavg)
  CALL average2(grid,newNH(1:M),depth2,newavg)

  NHflux_100 = (oldavg-newavg)*depth2/dt*3600.*24.

  CALL sum_g(grid,PON%old(1:M),h_m%old,oldavg)
  CALL sum_g(grid,PON%new(1:M),h_m%new,newavg)

  PONflux_ml = (oldavg-newavg)/dt*3600.*24.
 
  CALL sum_g(grid,N%O%new(1:M),h_m%old,oldavg)
  CALL sum_g(grid,newNO(1:M),h_m%new,newavg)

  NOflux_ml = (oldavg-newavg)/dt*3600.*24.
 
  CALL sum_g(grid,N%H%new(1:M),h_m%old,oldavg)
  CALL sum_g(grid,newNH(1:M),h_m%new,newavg)

  NHflux_ml = (oldavg-newavg)/dt*3600.*24.

  !PRINT "(A)","PONflux_100"
  !PRINT *,PONflux_100
  !PRINT "(A)","PONflux_200"
  !PRINT *,PONflux_200

END SUBROUTINE find_PON
