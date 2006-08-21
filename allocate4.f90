SUBROUTINE allocate4

  USE mean_param
  USE surface_forcing
  USE declarations
  
  IMPLICIT NONE

  
  INTEGER::y_y
        
  DO y_y = 1,grid%M
     IF (grid%d_g(y_y)>=150) THEN
        g_150 = y_y
        EXIT
     END IF
  END DO
          
  DO y_y = 1,grid%M
     IF (grid%d_g(y_y)>=100) THEN
        g_100 = y_y
        EXIT
     END IF
  END DO
  
  DO y_y = 1,grid%M
     IF (grid%d_g(y_y) >= 80) THEN
        g_80 = y_y
        EXIT
     END IF
  END DO

  DO y_y = 1,Csources
     ALLOCATE(species(y_y)%ingest(Cevent(y_y)%length),species(y_y)%Ex(Cevent(y_y)%length),&
          Hvector%z%c(y_y)%wt(Cevent(y_y)%length),Gvector_ro%z%c(y_y)%wt(Cevent(y_y)%length),&
          Gvector_ro_o%z%c(y_y)%wt(Cevent(y_y)%length),Copepod_wt1_p(Csources,max_length),PZ(M2),&
          Amatrix%null2%A(max_length),Amatrix%null2%B(max_length))
  END DO


!  ALLOCATE(null_vector2(Cevent(1)%length))
  ALLOCATE(null_vector2(M2))
  null_vector2 = 0.

END SUBROUTINE allocate4

