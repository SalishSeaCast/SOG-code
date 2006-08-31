SUBROUTINE P_H (M, PropOld, Gvector, Gvector_o, Gvector_ro, Gvector_ao, &
               Bmatrix_o, Hvector) 
      
      USE mean_param, only: trivector
      USE IMEX_constants, only: a_IMEX1

      IMPLICIT NONE
      
      integer :: M ! length of the grid
      double precision, dimension(0:M+1), intent(in):: PropOld !eg P%micro%old
      DOUBLE PRECISION, DIMENSION(M), INTENT(IN):: &
           Gvector, Gvector_o, &  ! surface and deep fluxes
           Gvector_ro, & ! biological model
           Gvector_ao  ! sinking
      TYPE(trivector), INTENT(IN)::Bmatrix_o  ! old diffusion

      DOUBLE PRECISION, DIMENSION(M), INTENT(OUT)::Hvector ! eg Hvector%p%micro

      ! local variables
      INTEGER::index

      Hvector = 0.

      DO index = 1, M
         Hvector(index) = PropOld(index) &
              + Gvector_ao(index) + Gvector_ro(index) &
              + (1.0 - a_IMEX1) * ( Gvector_o(index) + &
              Bmatrix_o%A(index) * PropOld(index-1) + &
              Bmatrix_o%B(index) * PropOld(index) + &
              Bmatrix_o%C(index) * PropOld(index+1)) + &
              a_IMEX1 * Gvector(index)

         IF (ABS(Hvector(index)) < EPSILON(Hvector(index))) THEN
            Hvector(index) = 0.
         END IF

      END DO

    END SUBROUTINE P_H
