SUBROUTINE matrix_A(Amatrix, Bmatrix)

      USE mean_param, only: trivector
      USE IMEX_constants, only: a_IMEX1

      IMPLICIT NONE

      TYPE(trivector), INTENT(OUT)::Amatrix !integration matrix      
      TYPE(trivector), INTENT(IN):: Bmatrix !integration matrix

      Amatrix%A = -a_IMEX1 * Bmatrix%A
      Amatrix%B = 1.0 - a_IMEX1 * Bmatrix%B
      Amatrix%C = -a_IMEX1 * Bmatrix%C

END SUBROUTINE matrix_A





