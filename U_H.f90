SUBROUTINE U_H (M, velold, Gvector, Gvector_o, Gvector_c, Gvector_co, &
  Bmatrix_o, Hvector)
  ! subroutine puts vert adv, bottom and surface fluxes (Gvector),
  ! Coriolis and prressure gradient fluxes (Gvector_c) and appropriate old
  ! diffusion (Bmatrix_o) into Hvector

  USE mean_param, only: trivector
  USE IMEX_constants, only: a_IMEX1

  IMPLICIT NONE

  integer, intent(in) :: M
  double precision, dimension(0:M+1):: velold   ! U%old etc
  DOUBLE PRECISION, DIMENSION(M), INTENT(IN):: &
       Gvector, Gvector_o, &  ! vert adv & surf/bot fluxes
       Gvector_c, Gvector_co  ! Coriolis and Pressure gradient forces     
  TYPE(trivector),INTENT(IN)::Bmatrix_o ! old diffusion

  DOUBLE PRECISION, DIMENSION(M), INTENT(OUT)::Hvector  ! zonal vector
 
  ! local variables
  INTEGER::index ! count through vertical

  Hvector = 0.

  DO index = 1, M
     Hvector(index) = velold(index) + &
          (1.0 - a_IMEX1) * (Gvector_o(index) + Gvector_co(index) + &
           Bmatrix_o%A(index) * velold(index-1) + &
           Bmatrix_o%B(index) * velold(index) + &
           Bmatrix_o%C(index)*velold(index+1)) + &
           a_IMEX1 * (Gvector(index) + Gvector_c(index))
  END DO

END SUBROUTINE U_H




