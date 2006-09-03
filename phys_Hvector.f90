! $Id$
! $Source$

subroutine phys_Hvector(M, qty_old, Gvector, Gvector_o, Gvector_c, &
     Gvector_co, Bmatrix_o, Hvector)
  ! Put vert adv, bottom and surface fluxes (Gvector),
  ! Coriolis and pressure gradient fluxes (Gvector_c) and appropriate old
  ! diffusion (Bmatrix_o) into Hvector.

  use precision_defs, only: dp
  use mean_param, only: trivector
  use IMEX_constants, only: a_IMEX1

  implicit none

  ! Arguments:
  integer, intent(in) :: M  ! Number of grid points
  real(kind=dp), dimension(0:M+1):: qty_old   ! U%old, etc.
  real(kind=dp), dimension(M), intent(in):: &
       Gvector, Gvector_o, &  ! vert adv & surf/bot fluxes
       Gvector_c, Gvector_co  ! Coriolis and Pressure gradient forces     
  type(trivector),intent(in)::Bmatrix_o ! old diffusion
  ! Result:
  real(kind=dp), dimension(M), intent(out)::Hvector  ! zonal vector
 
  ! Local variable:
  integer :: j ! loop index over depth

  ! *** This can probably be vectorized, but no point until we're inside
  ! *** a module, and M goes away by virtual of assumed-shape arrays.
  do j = 1, M
     Hvector(j) = qty_old(j) + &
          (1.0 - a_IMEX1) * (Gvector_o(j) + Gvector_co(j) + &
           Bmatrix_o%A(j) * qty_old(j-1) + &
           Bmatrix_o%B(j) * qty_old(j) + &
           Bmatrix_o%C(j)*qty_old(j+1)) + &
          a_IMEX1 * (Gvector(j) + Gvector_c(j))
  end do
end subroutine phys_Hvector
