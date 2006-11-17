! $Id$
! $Source$

subroutine allocate1(M, STAT)

      use mean_param
      use declarations

      implicit none
     
      ! Arguments:
      integer, intent(in) :: M  ! Number of grid points
      integer, dimension(20), intent(out) :: STAT  ! Memory allocation status

      ALLOCATE(&
           K%u%total(0:M),K%s%total(0:M),K%t%total(0:M),&
           K%t%all(1:M),K%u%all(1:M),K%s%all(1:M),&
           K%u%ML(1:M), K%s%ML(1:M), K%t%ML(1:M), &
           STAT = alloc_stat(5))

      ALLOCATE(Fw(0:M), F_n(0:M), &
           Q_n(0:M), &
           I(0:M), I_par(0:M), &
           omega%s%value(0:M), omega%m%value(0:M), &
           Ri_b(0:M), N_2_g(M), V_t_square(M), &
           STAT = alloc_stat(6)) 
                     
      ALLOCATE(gamma%m(0:M),gamma%s(0:M),gamma%t(0:M), &
           STAT = alloc_stat(7))

      ALLOCATE(Gvector%s(M),Gvector%t(M),Gvector%u(M),Gvector%v(M),&
               STAT = alloc_stat(8))
               
      ALLOCATE(micro%growth%light(M), micro%growth%new(M), &
           nano%growth%light(M), nano%growth%new(M), &
           micro%Nlimit(M), nano%Nlimit(M), &
           wupwell(0:M), &
           STAT = alloc_stat(10))

      ALLOCATE(f_ratio(M), &
           STAT = alloc_stat(11))

      ALLOCATE(G_shape%s(0:M+1), G_shape%m(0:M+1), G_shape%t(0:M+1), &
           STAT = alloc_stat(16))
      
END SUBROUTINE allocate1



