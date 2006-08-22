SUBROUTINE P_H(mm,Hvector,Gvect,Gvect_o,Gvect_o_o,Gvect_ro,Gvect_ro_o,Gvect_ao,Gvect_ao_o,&
               Bmat_o,Bmat_o_o,PP)
      
      USE mean_param
      USE IMEX_constants

      IMPLICIT NONE
      
      TYPE(gr_d), INTENT(IN)::mm        !grid
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(OUT)::Hvector     !Hvector%p%micro
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::Gvect,Gvect_o,Gvect_o_o,&
                                                     Gvect_ro,Gvect_ro_o,Gvect_ao,Gvect_ao_o
      TYPE(trivector), INTENT(IN)::Bmat_o,Bmat_o_o
      TYPE(prop), INTENT(IN)::PP        !P%micro
 !     INTEGER, INTENT(IN)::tstep

      INTEGER::index

      Hvector = 0.

         DO index = 1, mm%M
            Hvector(index) = PP%old(index)+Gvect_ao(index)+Gvect_ro(index)+(1.0-a_IMEX1)*(Gvect_o(index) + &
            Bmat_o%A(index)*PP%old(index-1) + Bmat_o%B(index)*PP%old(index) + &
            Bmat_o%C(index)*PP%old(index+1)) + a_IMEX1*Gvect(index)

            IF (ABS(Hvector(index)) < EPSILON(Hvector(index))) THEN
               Hvector(index) = 0.
            END IF

         END DO

    !  ELSE                      !!!Second Order IMEX!!!

    !     DO index = 1, mm%M
    !        Hvector(index) = (2.0*a_IMEX2*PP%old(index) - (a_IMEX2-0.5)*PP%old_old(index) + &
    !        (a_IMEX2+1.0)*(Gvect_ao(index)+Gvect_ro(index))-a_IMEX2*(Gvect_ao_o(index)+Gvect_ro_o(index))+&
    !        (a_IMEX2+b_IMEX2/2.0)*Gvect(index) + &
    !        (1.0-a_IMEX2-b_IMEX2)*(Gvect_o(index) + Bmat_o%A(index)*PP%old(index-1) + &
    !        Bmat_o%B(index)*PP%old(index) + Bmat_o%C(index)*PP%old(index+1)) + &
    !        b_IMEX2/2.0*(Gvect_o_o(index) + Bmat_o_o%A(index)*PP%old_old(index-1) + &
    !        Bmat_o_o%B(index)*PP%old_old(index) + Bmat_o_o%C(index)*PP%old_old(index+1)))/(a_IMEX2+0.5)


            !Hvector(index) = PP%old(index) + Gvect_ao(index)  !BOTT's Scheme
            !Hvector(index) = PP%old(index) + Gvect_ro(index)  !First Order Scheme
            !Hvector(index) = PP%old(index) + 3.0/2.0*Gvect_ao(index) - 1.0/2.0*Gvect_ao_o(index)
            !Hvector(index) = PP%old_old(index) + 2.0*Gvect_ao(index)
            !Hvector(index) = 4.0/3.0*PP%old(index)-1.0/3.0*PP%old_old(index)+4.0/3.0*Gvect_ao(index)-&
            !     2.0/3.0*Gvect_ao_o(index)  
    !        IF (ABS(Hvector(index)) < EPSILON(Hvector(index))) THEN
    !           Hvector(index) = 0.
    !        END IF
    !     END DO

     ! END IF

END SUBROUTINE P_H
