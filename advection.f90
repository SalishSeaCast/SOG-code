! $Id$
! $Source$

subroutine advection(mm, sink, pp, step, Gvector)
  ! *** What's it do?

  use precision_defs, only: dp
  use grid_mod, only: grid_

  implicit none

  ! Arguments:
  type(grid_), intent(in) :: mm
  real(kind=dp), intent(in) :: sink                       ! sinking velocity
  real(kind=dp), dimension(0:mm%M+1), intent(in) :: pp    ! P%micro%old
  real(kind=dp), intent(in) :: step                       ! timestep
  real(kind=dp), dimension(mm%M), intent(out) :: Gvector  ! Gvector_ao%p%micro

  ! Local variables:
  REAL(KIND=DP), DIMENSION(0:mm%M+1)::Ru_1, Ru_a, &
       Au_1, Au_2, delta_1
  REAL(KIND=DP), DIMENSION(0:mm%M)::Fu
  REAL(KIND=DP), DIMENSION(0:mm%M+1)::Ru, Ru_b
  REAL(KIND=DP), DIMENSION(-1:mm%M)::Ru_2,delta_2 
  REAL(KIND=DP), DIMENSION(0:mm%M+1,0:3)::aa  !for a third order advection scheme
  INTEGER, DIMENSION(0:mm%M+1)::L    !order as a function of grid spacing
  INTEGER::index, jk


  Ru_a = pp
  !Ru_a(0) = 0.

  !Use pp(0) = 0 for upper boundary condition in advection!!!!

  DO index = 0,mm%M+1
     IF (sink >= 0.) THEN
        Au_1(index) = sink   ! trivial for sink = constant
        Au_2(index) =  0.
     ELSE
        Au_1(index) = 0.
        Au_2(index) = -sink     
     END IF
  END DO

!!!coefficients for cubic polynomial (3C)***see Easter 1993, Monthly Weather Review 121: pp 297-304!!!

  aa = 0.

  DO index = 2,mm%M-1

     aa(index,0) = (-Ru_a(index+1) + 26.0*Ru_a(index) - Ru_a(index-1))/24.0
     aa(index,1) = (-5.0*Ru_a(index+2)+34.0*Ru_a(index+1)-34.0*Ru_a(index-1)+5.0*Ru_a(index-2))/48.0
     aa(index,2) = (Ru_a(index+1)-2.0*Ru_a(index)+Ru_a(index-1))/2.0
     aa(index,3) = (Ru_a(index+2)-2.0*Ru_a(index+1)+2.0*Ru_a(index-1)-Ru_a(index-2))/12.0
     L(index) = 3

  END DO



!!!!Boundary Conditions!!!!!!!!!!!!!!!!!!!!!

  L(0) = 0
  aa(0,0) = Ru_a(0)

  L(1) = 2
  aa(1,0) = (-Ru_a(2) + 26.0*Ru_a(1) - Ru_a(0))/24.0
  aa(1,1) = (Ru_a(2) - Ru_a(0))/2.0
  aa(1,2) = (Ru_a(2) - 2.0*Ru_a(1) + Ru_a(0))/2.0

  L(mm%M) = 2
  aa(mm%M,0) = (-Ru_a(mm%M+1) + 26.0*Ru_a(mm%M) - Ru_a(mm%M-1))/24.0
  aa(mm%M,1) = (Ru_a(mm%M+1) - Ru_a(mm%M-1))/2.0
  aa(mm%M,2) = (Ru_a(mm%M+1) - 2.0*Ru_a(mm%M) + Ru_a(mm%M-1))/2.0

  L(mm%M+1) = 0
  aa(mm%M+1,0) = Ru_a(mm%M+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  delta_1 = 0.  
  delta_2 = 0.  
  Ru_1 = 0.     
  Ru_2 = 0.
  Fu = 0.

  do index = 1, mm%M
     delta_1(index) = Au_1(index) * step / mm%i_space(index)
     delta_2(index-1) = Au_2(index-1) * step/ mm%i_space(index)
  end do

  do index = 0, mm%M + 1
     do jk = 0, L(index)
        if (jk == 0) then
           Ru_1(index) = Ru_1(index) + aa(index,0)
           Ru_2(index-1) = Ru_2(index-1) + aa(index,0) 
        else
           Ru_1(index) = Ru_1(index) + aa(index,jk) &
                / (2.)**(jk+1) / (jk+1)             &
                * (1. - (1. - 2. * (delta_1(index)  &
                + epsilon(sink)))**(jk+1))          &
                / (delta_1(index) + epsilon(sink))
           Ru_2(index-1) = Ru_2(index-1) + aa(index,jk) &
                / (-2.)**(jk+1) / (jk+1)                &
                * ((1. - 2. * (delta_2(index-1)         &
                + epsilon(sink)))**(jk+1) - 1.)         &
                / (delta_2(index-1) + epsilon(sink)) 
        end if
     end do
     if (Ru_1(index) < 0.) then
        Ru_1(index) =  0.
     end if
     ! *** Changing index to index-1 in the next 2 statements is a guess
     ! *** It's consistent with the use of Ru_2 above, and it avoids an
     ! *** array bound runtime error
     if (Ru_2(index-1) < 0.) then
        Ru_2(index-1) = 0.
     end if
  end do

  do index = 0, mm%M+1
     Ru_b(index) = delta_1(index) * Ru_1(index) &
          + delta_2(index-1) * Ru_2(index-1)
     Ru(index) = max(epsilon(sink), Ru_a(index), Ru_b(index))
  end do

  !Fu(0) = -Ru_2(0)*Au_2(0)*Ru_a(1)/Ru(1)
  !Fu(mm%M) = Ru_1(mm%M)*Au_1(mm%M)*Ru_a(mm%M)/Ru(mm%M)
  do index = 0, mm%M    !1, mm%M-1
     Fu(index) = Ru_1(index) * Au_1(index) * Ru_a(index) / Ru(index) &
          - Ru_2(index) * Au_2(index) * Ru_a(index+1) / Ru(index+1)
  end do

  !        open(558,file="output/detritus.dat")
  !        do index=1,mm%M
  !        write(558,*)Fu(index)
  !        end do
  !       close(558)      



  !TRY!!!
  Fu(0) = 0.
  Gvector = 0.


  DO index = 1, mm%M
     Gvector(index) = -step/mm%i_space(index)*(Fu(index)-Fu(index-1))
     !Gvector(index) = -step/mm%i_space(index)/2.0*sink*(pp(index+1)-pp(index-1))
     IF (ABS(Gvector(index)) < EPSILON(Gvector(index))) THEN
        Gvector(index) = 0.
     end if
  end do

end subroutine advection
