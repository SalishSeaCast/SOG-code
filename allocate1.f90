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
           Q_n(0:M), &
           I(0:M), I_par(0:M), &
           STAT = alloc_stat(6)) 
                     
      ALLOCATE(Gvector%s(M),Gvector%t(M),Gvector%u(M),Gvector%v(M),&
               STAT = alloc_stat(8))
               
      ALLOCATE(micro%growth%light(M), micro%growth%new(M), &
           nano%growth%light(M), nano%growth%new(M), &
           pico%growth%light(M), pico%growth%new(M), &
           micro%Nlimit(M), nano%Nlimit(M), pico%Nlimit(M), &
           wupwell(0:M), &
           STAT = alloc_stat(10))

END SUBROUTINE allocate1



