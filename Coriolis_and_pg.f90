! $Id$
! $Source$

module Coriolis_and_pg_mod

  use precision_defs

  implicit none

contains

  subroutine Coriolis_and_pg(dt, vel, pbgrad, Gvector_c)
    ! Calculate the Coriolis and baroclinic pressure gradient
    ! components of the G vector for the specified velocity component.
    use physics_model, only: f
    implicit none
    ! Arguments:
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(0:), intent(in)::vel
    real(kind=dp), dimension(1:), intent(in)::pbgrad
    real(kind=dp), dimension(1:), intent(out)::Gvector_c

    Gvector_c = (f * vel(1:) - pbgrad) * dt
  end subroutine Coriolis_and_pg

end module Coriolis_and_pg_mod
