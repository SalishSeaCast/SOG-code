! $Id$
! $Source$

SUBROUTINE smoothdata(old, size, smooth)

  USE mean_param
  USE surface_forcing

  IMPLICIT NONE

  INTEGER, INTENT(IN)::size
  DOUBLE PRECISION, DIMENSION(size), INTENT(IN)::old
  DOUBLE PRECISION, DIMENSION(0:size+1), INTENT(OUT)::smooth

  DOUBLE PRECISION, DIMENSION(size + 2*runsize)::temp
  INTEGER::k

  temp(1+runsize:runsize+size) = old(1:size)
  temp(1:runsize) = old(size-runsize+1:size)
  temp(runsize+size+1:size+2*runsize) = old(1:runsize)

  DO k = 1,365
     smooth(k) = SUM(temp(k:2*runsize+k))/(2*runsize+1)
  END DO

  !PRINT "(A)","temp"
  !PRINT *,temp
  !DO k = 1, 365
  !   PRINT "(A)","old(k), smooth(k),k"
  !   PRINT *,old(k),smooth(k),k
  !END DO
  !STOP

END SUBROUTINE smoothdata
