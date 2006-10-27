! $Id$
! $Source$

SUBROUTINE define_sog(timestep)
  
  !!store previous time_steps in old (n)
  use core_variables, only: U, V, T, S
  USE declarations, only: D_bins, h, &
       Bmatrix, Bmatrix_o, & 
       Gvector, Gvector_o, Gvector_c, Gvector_co

      IMPLICIT NONE

      INTEGER, INTENT(IN)::timestep
      INTEGER::kk


      IF (timestep > 1) THEN
         Bmatrix_o%u%A = Bmatrix%u%A
         Bmatrix_o%u%B = Bmatrix%u%B
         Bmatrix_o%u%C = Bmatrix%u%C         
         Bmatrix_o%t%A = Bmatrix%t%A    !check
         Bmatrix_o%t%B = Bmatrix%t%B
         Bmatrix_o%t%C = Bmatrix%t%C         
         Bmatrix_o%s%A = Bmatrix%s%A
         Bmatrix_o%s%B = Bmatrix%s%B
         Bmatrix_o%s%C = Bmatrix%s%C

         Bmatrix_o%bio%A = Bmatrix%bio%A
         Bmatrix_o%bio%B = Bmatrix%bio%B
         Bmatrix_o%bio%C = Bmatrix%bio%C                 

         Gvector_o%u = Gvector%u
         Gvector_o%v = Gvector%v
         Gvector_o%s = Gvector%s
         Gvector_o%t = Gvector%t

         Gvector_o%p%micro = Gvector%p%micro
         Gvector_o%p%nano = Gvector%p%nano
         Gvector_o%n%o = Gvector%n%o
         Gvector_o%n%h = Gvector%n%h
         Gvector_o%si = Gvector%si
        
         Gvector_co%u = Gvector_c%u
         Gvector_co%v = Gvector_c%v

         DO kk = 1, D_bins
            Gvector_o%d(kk)%bin = Gvector%d(kk)%bin 
         END DO
         
      END IF

         U%old = U%new
         V%old = V%new
         S%old = S%new
         T%old = T%new
         h%old = h%new
end subroutine define_sog



