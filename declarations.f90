module declarations
! *** What's in here?

  use precision_defs, only: dp

  implicit NONE

  ! Heat fluxes
  real(kind=dp) :: Q_t  ! Turbulent surface heat flux

  INTEGER :: &
       time_step, &
       count
  INTEGER::M2

end module declarations
