SUBROUTINE N_H(mm,Hvector,Gvect,Gvect_o,Gvect_o_o,Gvect_ro,Gvect_ro_o,Bmat_o,Bmat_o_o,NN)

      USE mean_param
      USE IMEX_constants

      IMPLICIT NONE

      TYPE(gr_d), INTENT(IN)::mm        !grid
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(OUT)::Hvector     !Hvector%n%o
      DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::Gvect,Gvect_o,Gvect_o_o,&
                                                     Gvect_ro,Gvect_ro_o
      TYPE(trivector), INTENT(IN)::Bmat_o,Bmat_o_o
      TYPE(prop), INTENT(IN)::NN        !N%O
!      INTEGER, INTENT(IN)::tstep
      
      INTEGER::index

      Hvector = 0.

         DO index = 1, mm%M
            Hvector(index) = NN%old(index) + Gvect_ro(index)+ (1.0-a_IMEX1)*(Gvect_o(index) + &
            Bmat_o%A(index)*NN%old(index-1) + Bmat_o%B(index)*NN%old(index) + &
            Bmat_o%C(index)*NN%old(index+1)) + a_IMEX1*Gvect(index)            

            IF (ABS(Hvector(index)) < EPSILON(Hvector(index))) THEN
               Hvector(index) = 0.
            END IF
         END DO

    ! ELSE

    !     DO index = 1, mm%M
    !        Hvector(index) = (2.0*a_IMEX2*NN%old(index) - (a_IMEX2-0.5)*NN%old_old(index) + &
    !        (a_IMEX2+1.0)*Gvect_ro(index) - a_IMEX2*Gvect_ro_o(index)  + &
    !        (a_IMEX2+b_IMEX2/2.0)*Gvect(index) + &
    !        (1.0-a_IMEX2-b_IMEX2)*(Gvect_o(index) + Bmat_o%A(index)*NN%old(index-1) + &
    !        Bmat_o%B(index)*NN%old(index) + Bmat_o%C(index)*NN%old(index+1)) + &
    !        b_IMEX2/2.0*(Gvect_o_o(index) + Bmat_o_o%A(index)*NN%old_old(index-1) + &
    !        Bmat_o_o%B(index)*NN%old_old(index) + Bmat_o_o%C(index)*NN%old_old(index+1)))/(a_IMEX2+0.5)
 
           ! Hvector(index) = NN%old(index) + Gvect_ro(index)  !First order scheme
    !        IF (ABS(Hvector(index)) < EPSILON(Hvector(index))) THEN
    !           Hvector(index) = 0.
    !        END IF
    !     END DO

    !  END IF

END SUBROUTINE N_H







