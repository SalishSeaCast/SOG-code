module declarations
! *** What's in here?

  use precision_defs, only: dp

  implicit NONE

  ! Plankton growth
  type :: grow               
     real(kind=dp), dimension(:), pointer :: light, new
  end type grow
  type :: plankton_growth
     type(grow) :: growth
     real(kind=dp), dimension(:), pointer :: Nlimit
  end type plankton_growth
  type(plankton_growth) :: micro, nano, pico

  ! Heat fluxes
  real(kind=dp) :: Q_t  ! Turbulent surface heat flux

  INTEGER, DIMENSION(20)::alloc_stat  !***
  INTEGER :: &
       time_step, &
       xx, &
       count
  INTEGER::M2

end module declarations
