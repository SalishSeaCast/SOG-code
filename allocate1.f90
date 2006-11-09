! $Id$
! $Source$

subroutine allocate1(M, STAT)

      use mean_param
      use declarations

      implicit none
     
      ! Arguments:
      integer, intent(in) :: M  ! Number of grid points
      integer, dimension(20), intent(out) :: STAT  ! Memory allocation status

      ALLOCATE(w%b(0:M),w%b_err(0:M),&
           w%t(0:M),w%s(0:M),w%u(0:M),&
               w%v(0:M),w%p%micro(0:M),&
               STAT = alloc_stat(4))

      ALLOCATE(K%u%shear(0:M),&
           K%s%dd(1:M), K%t%dd(1:M),&
           K%u%total(0:M),K%s%total(0:M),K%t%total(0:M),&
           K%t%all(0:M),K%u%all(0:M),K%s%all(0:M),&
           STAT = alloc_stat(5))

      ALLOCATE(Fw(0:M), F_n(0:M), &
           Q_n(0:M), &
           I(0:M), I_par(0:M), &
           phi%m%value(0:M), phi%s%value(0:M), &
           Ri_b(0:M), N_2_g(M), V_t_square(M), &
           omega%s%value(0:M), omega%m%value(0:M), &
           STAT = alloc_stat(6)) 
                     
      ALLOCATE(gamma%m(0:M),gamma%s(0:M),gamma%t(0:M), &
           STAT = alloc_stat(7))

      ALLOCATE(&
               Gvector%s(M),Gvector%t(M),Gvector%u(M),Gvector%v(M),Gvector_o%s(M),Gvector_o%t(M),&
               Gvector_o%u(M),Gvector_o%v(M),&
               Gvector_c%u(M),Gvector_c%v(M),Gvector_co%u(M),&
               Gvector_co%v(M),&
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
      ALLOCATE(K%u%ML(0:M+1), K%s%ML(0:M+1), K%t%ML(0:M+1), &
           STAT = alloc_stat(17))
               
      
END SUBROUTINE allocate1



