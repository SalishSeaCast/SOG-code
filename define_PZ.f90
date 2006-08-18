SUBROUTINE define_PZ(check)

!  USE mean_param
  USE declarations
!  USE surface_forcing

  IMPLICIT NONE

  INTEGER::check,isusan

  !Copepod variables
  !n(i), Ntot and Z_new(i) where i is grid point:
  !sum_i {Z_new(i)*i_space(i)}/ sum_i {n(i)*i_space(i)} = Ntot ==> average number per m^2
  !as well, n(i) is normalized such that sum_i {n(i)*i_space(i)} = 1


  PZ = 0.

  DO xx = 1, M2
     IF (xx <= M) THEN
        PZ(xx) = P%micro%new(xx)
     ELSE IF (xx > M .AND. xx <= 2*M) THEN
        PZ(xx) = P%nano%new(xx-M)
     ELSE IF (xx > 2*M .AND. xx <= 3*M) THEN
        PZ(xx) = N%O%new(xx-2*M)
     ELSE IF (xx > 3*M .AND. xx <= 4*M) THEN
        PZ(xx) = N%H%new(xx-3*M)
     ELSE IF (xx > 4*M .AND. xx <= 4*M+D_bins*M) THEN 
        DO yy = 1,D_bins    ! 1 is suspended particulate N (SUS), 2 is sinking PN (SPN), and 3 is fast sink
           IF (xx > 4*M+(yy-1)*M .AND. xx <= 4*M+yy*M) THEN
             PZ(xx) = Detritus(yy)%D%new(xx-(4*M+(yy-1)*M))
           END IF
        END DO
     END IF
  END DO

! Nanoplanktons are placed after microplanktons. 
! The whole scheme is different now. See above
!  DO xx = 1, M2
!     IF (xx <= M) THEN
!        PZ(xx) = P%micro%new(xx)
!     ELSE IF (xx > M .AND. xx <= 2*M) THEN
!        PZ(xx) = N%O%new(xx-M)
!     ELSE IF (xx > 2*M .AND. xx <= 3*M) THEN
!        PZ(xx) = N%H%new(xx-2*M)
!     ELSE IF (xx > 3*M .AND. xx <= 3*M+D_bins*M) THEN 
!        DO yy = 1,D_bins    ! 1 is suspended particulate N (SUS), 2 is sinking PN (SPN), and 3 is fast sink
!           IF (xx > 3*M+(yy-1)*M .AND. xx <= 3*M+yy*M) THEN
!             PZ(xx) = Detritus(yy)%D%new(xx-(3*M+(yy-1)*M))
!           END IF
!        END DO
!     END IF
!  END DO

!      open(558,file="output/PZ.dat")!
!      do isusan=1,M2
!      write(558,*)PZ(isusan)
!      enddo
!      close(558)      



END SUBROUTINE define_PZ

