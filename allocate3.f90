! $Id$
! $Source$

subroutine allocate3
  ! Allocate memory for biological model arrays (at least some of them...)

  use malloc, only: alloc_check
  use declarations, only: M, M2, D_bins, PZ, Gvector, Gvector_o, &
       Gvector_o_o, Gvector_ao, Gvector_ao_o, Gvector_ro,        &
       Gvector_ro_o, Hvector, Detritus

  implicit none

  ! Local variables:
  integer              :: allocstat  ! Allocation return status
  character(len=80)    :: msg        ! Allocation failure message prefix
  integer              :: y_y        ! Loop index

  msg = "PZ vector (biology model timestep initial conditions) array"
  allocate(PZ(M2), stat=allocstat)
  call alloc_check(allocstat, msg)

  do y_y = 1, D_bins
     msg = "Biology model Gvector bin components arrays"
     allocate(Gvector%d(y_y)%bin(M), Gvector_o%d(y_y)%bin(M),     &
          Gvector_o_o%d(y_y)%bin(M), Gvector_ao%d(y_y)%bin(M),    &
          Gvector_ao_o%d(y_y)%bin(M), Gvector_ro%d(y_y)%bin(M),   &
          Gvector_ro_o%d(y_y)%bin(M), stat=allocstat)
     call alloc_check(allocstat, msg)
     msg = "Biology model Hvector bin component array"
     allocate(Hvector%d(y_y)%bin(M), stat=allocstat)
     call alloc_check(allocstat, msg)
     msg = "Detritus component arrays"
     allocate(Detritus(y_y)%D%new(0:M+1), Detritus(y_y)%D%old(0:M+1), &
          Detritus(y_y)%D%old_old(0:M+1), stat=allocstat)
     call alloc_check(allocstat, msg)
  end do
      
end subroutine allocate3
