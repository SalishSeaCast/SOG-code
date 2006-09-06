! $Id$
! $Source$

MODULE pdf

  IMPLICIT NONE

CONTAINS
  
  SUBROUTINE normalize(wt_pdf,bin_width,size)

    INTEGER, INTENT(IN)::size
    DOUBLE PRECISION, DIMENSION(size), INTENT(IN)::bin_width
    DOUBLE PRECISION, DIMENSION(size), INTENT(IN OUT)::wt_pdf
    
    INTEGER::j
    DOUBLE PRECISION::norm

    norm = 0.
    DO j = 1, size
       norm = norm + wt_pdf(j)*bin_width(j)
    END DO
    DO j = 1,size
       IF (norm /= 0.) THEN
          wt_pdf(j) = wt_pdf(j)/norm
       END IF
    END DO
    
  END SUBROUTINE normalize

  SUBROUTINE pdf_avg(wt,wt_pdf,size,wt_avg)

    INTEGER, INTENT(IN)::size
    DOUBLE PRECISION, DIMENSION(size), INTENT(IN)::wt,wt_pdf
    DOUBLE PRECISION, INTENT(IN OUT)::wt_avg

    INTEGER::j

    wt_avg = 0.
    DO j = 1,size
       wt_avg = wt_avg + wt(j)*wt_pdf(j)
    END DO
    
  END SUBROUTINE pdf_avg


END MODULE pdf

