! $Id$
! $Source$

subroutine allocate3(M)
  ! Allocate memory for biological model arrays (at least some of them...)

  use malloc, only: alloc_check
  use declarations, only: D_bins, Gvector, Gvector_o, &
       Gvector_ao, Hvector

  implicit none

  ! Argument:
  integer, intent(in)  :: M          ! Number of grid points
  ! Local variables:
  integer              :: allocstat  ! Allocation return status
  character(len=80)    :: msg        ! Allocation failure message prefix
  integer              :: y_y        ! Loop index

  do y_y = 1, D_bins
     msg = "Biology model Gvector bin components arrays"
     allocate(Gvector%d(y_y)%bin(M), Gvector_o%d(y_y)%bin(M), &
          Gvector_ao%d(y_y)%bin(M), &
          stat=allocstat)
     call alloc_check(allocstat, msg)
     msg = "Biology model Hvector bin component array"
     allocate(Hvector%d(y_y)%bin(M), stat=allocstat)
     call alloc_check(allocstat, msg)
  end do
      
end subroutine allocate3
