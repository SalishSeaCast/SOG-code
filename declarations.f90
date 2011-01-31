module declarations
! *** What's in here?

  use precision_defs, only: dp

  implicit NONE

  ! Grow component of plankton2 type
  type :: grow               
     real(kind=dp), dimension(:), pointer :: light, new
  end type grow

  ! *** Plankton "behaviours" ??? type
  ! *** This will probably be refactored when Susan implements sinking
  type :: plankton2
     type(grow) :: growth
     real(kind=dp), dimension(:), pointer:: Nlimit
  end type plankton2

  TYPE(plankton2)::micro, nano, pico

  ! Heat fluxes
  real(kind=dp) :: Q_t  ! Turbulent surface heat flux

  INTEGER, DIMENSION(20)::alloc_stat  !***
  INTEGER :: &
       time_step, &
       xx, &
       count
  INTEGER::M2

end module declarations
