SUBROUTINE define_sog(timestep)

      USE mean_param
      USE declarations

      IMPLICIT NONE

      INTEGER, INTENT(IN)::timestep
      INTEGER::kk
                                !!store previous time_steps in old (n) and old_old (n-1)!!

      IF (timestep > 1) THEN
         U%old_old = U%old
         V%old_old = V%old
         S%old_old = S%old
         T%old_old = T%old
         h%old_old = h%old
         P%micro%old_old = P%micro%old
         P%nano%old_old = P%nano%old
         N%O%old_old = N%O%old
         N%H%old_old = N%H%old

         PZ_old = PZ       
         B%old = B%new   
         density%old = density%new
         K%t%old = K%t%all
         K%s%old = K%s%all
         K%u%old = K%u%all

         Bmatrix_o_o%u%A = Bmatrix_o%u%A
         Bmatrix_o_o%u%B = Bmatrix_o%u%B
         Bmatrix_o_o%u%C = Bmatrix_o%u%C         
         Bmatrix_o_o%t%A = Bmatrix_o%t%A
         Bmatrix_o_o%t%B = Bmatrix_o%t%B      !check
         Bmatrix_o_o%t%C = Bmatrix_o%t%C         
         Bmatrix_o_o%s%A = Bmatrix_o%s%A
         Bmatrix_o_o%s%B = Bmatrix_o%s%B
         Bmatrix_o_o%s%C = Bmatrix_o%s%C

                 
         Bmatrix_o_o%bio%A = Bmatrix_o%bio%A
         Bmatrix_o_o%bio%B = Bmatrix_o%bio%B
         Bmatrix_o_o%bio%C = Bmatrix_o%bio%C               
         Bmatrix_o_o%no%A = Bmatrix_o%no%A
         Bmatrix_o_o%no%B = Bmatrix_o%no%B
         Bmatrix_o_o%no%C = Bmatrix_o%no%C 

         Gvector_o_o%u = Gvector_o%u
         Gvector_o_o%v = Gvector_o%v
         Gvector_o_o%s = Gvector_o%s        !check
         Gvector_o_o%t = Gvector_o%t

         Gvector_o_o%p%micro = Gvector_o%p%micro
         Gvector_o_o%p%nano = Gvector_o%p%nano
         Gvector_o_o%n%o = Gvector_o%n%o
         Gvector_o_o%n%h = Gvector_o%n%h

         Gvector_co_o%u = Gvector_co%u
         Gvector_co_o%v = Gvector_co%v           !check

         Gvector_ao_o%p%micro = Gvector_ao%p%micro
!         Gvector_ao_o%p%nano = Gvector_ao%p%nano

         Gvector_ro_o%p%micro = Gvector_ro%p%micro
         Gvector_ro_o%p%nano = Gvector_ro%p%nano
         Gvector_ro_o%p%micro_Q = Gvector_ro%p%micro_Q
         Gvector_ro_o%p%nano_Q = Gvector_ro%p%nano_Q
         Gvector_ro_o%n%o = Gvector_ro%n%o
         Gvector_ro_o%n%h = Gvector_ro%n%h

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
         Bmatrix_o%no%A = Bmatrix%no%A
         Bmatrix_o%no%B = Bmatrix%no%B
         Bmatrix_o%no%C = Bmatrix%no%C 

         Gvector_o%u = Gvector%u
         Gvector_o%v = Gvector%v
         Gvector_o%s = Gvector%s
         Gvector_o%t = Gvector%t

         Gvector_o%p%micro = Gvector%p%micro
         Gvector_o%p%nano = Gvector%p%nano
         Gvector_o%n%o = Gvector%n%o
         Gvector_o%n%h = Gvector%n%h

        
         Gvector_co%u = Gvector_c%u
         Gvector_co%v = Gvector_c%v

         DO kk = 1, D_bins
            Detritus(kk)%D%old_old = Detritus(kk)%D%old
            Gvector_o_o%d(kk)%bin = Gvector_o%d(kk)%bin
            Gvector_ro_o%d(kk)%bin =Gvector_ro%d(kk)%bin
            Gvector_o%d(kk)%bin = Gvector%d(kk)%bin !crashes here
            Gvector_ao_o%d(kk)%bin = Gvector_ao%d(kk)%bin
         END DO
         nano%Q_old_old = nano%Q_old !V.flagella.01 not sure to keep it or not

         
      END IF

         U%old = U%new
         V%old = V%new
         S%old = S%new
         T%old = T%new
         h%old = h%new
         P%micro%old = P%micro%new
         P%nano%old = P%nano%new
         N%O%old = N%O%new
         N%H%old = N%H%new
         PON%old = PON%new

         ut%old = ut%new
         vt%old = vt%new

         DO kk = 1, D_bins
            Detritus(kk)%D%old = Detritus(kk)%D%new
         END DO

END SUBROUTINE define_sog



