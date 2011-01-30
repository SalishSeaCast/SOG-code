module mean_param
  ! This should probably be renamed to something more descriptive like
  ! type_defs.
  ! It would also be good to re-organize somehow so that they types
  ! hierarchy is more evident.

  use precision_defs, only: dp

  implicit none

  ! Grow component of plankton2 type
  type :: grow               
     real(kind=dp), dimension(:), pointer :: light, new
  end type grow

  ! *** Plankton "behaviours" ??? type
  ! *** This will probably be refactored when Susan implements sinking
  type :: plankton2
     type(grow) :: growth
     real(kind=dp) :: sink_min, sink_max  ! Sinking velocity
     real(kind=dp), dimension(:), pointer:: Nlimit
  end type plankton2

  INTEGER::is_leap_year, was_leap_year

END MODULE mean_param
