! $Id$
! $Source$

subroutine find_new(M)
  ! Copies results from tridiagonal solver into real variable names,
  !removing any negative values.  Applies upper boundary conditions

  use surface_forcing, only: zero
  use declarations, only: P1_p, Pnano1_p, NO1_p, NH1_p, SIL1_p, Detritus1_p,&
       P, N, Sil, Detritus, D_bins, neg_count

  implicit none

  ! Argument:
  integer :: M  ! Number of grid points

  ! Local variables:
  INTEGER:: xx, yy

  DO yy = 1, M
     IF (P1_p(yy) < zero) THEN
        P1_p(yy) = zero
        neg_count = neg_count + 1
     END IF
     IF (Pnano1_p(yy) < zero) THEN
        Pnano1_p(yy)= zero
        neg_count = neg_count + 1
     END IF
     IF (NO1_p(yy) < zero) THEN
        NO1_p(yy) = zero
        neg_count = neg_count + 1
     END IF
     IF (NH1_p(yy) < zero) THEN
        NH1_p(yy) = zero
        neg_count = neg_count + 1
     END IF
     IF (SIL1_p(yy) < zero) THEN
        SIL1_p(yy) = zero
        neg_count = neg_count + 1
     END IF

     DO xx = 1, D_bins
        IF (Detritus1_p(xx,yy) < zero) THEN
           Detritus1_p(xx,yy) = zero
           neg_count = neg_count + 1
        END IF
        Detritus(xx)%D%new(yy) = Detritus1_p(xx,yy)
     END DO

     P%micro%new(yy) = P1_p(yy)
     P%nano%new(yy) = Pnano1_p(yy)
     N%O%new(yy) = NO1_p(yy)
     N%H%new(yy) = NH1_p(yy)
     Sil%new(yy) = SIL1_p(yy)

  END DO

  P%micro%new(0) = P%micro%new(1)
  P%nano%new(0) = P%nano%new(1)
  N%O%new(0) = N%O%new(1)
  N%H%new(0) = N%H%new(1)
  Sil%new(0) = Sil%new(1)

  do xx = 1, D_bins
     Detritus(xx)%D%new(0) = Detritus(xx)%D%new(1)
  end do

end subroutine find_new
