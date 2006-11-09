! $Id$
! $Source$

module mean_param
  ! This should probably be renamed to something more descriptive like
  ! type_defs.
  ! It would also be good to re-organize somehow so that they types
  ! hierarchy is more evident.

  use precision_defs, only: dp

  implicit none

  ! Phytoplankton component of UVST type
  type :: phyto                     
     real(kind=dp), dimension(:), pointer :: micro, nano
  end type phyto

  ! Type for Gvectors and Hvector used in implicit solver
  type :: UVST              
     real(kind=dp), dimension(:), pointer :: u, v, s, t
  end type UVST

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



  TYPE :: diff
     DOUBLE PRECISION, DIMENSION(:), POINTER :: shear, total 
     DOUBLE PRECISION, DIMENSION(:), POINTER :: dd    !double diffusion
     DOUBLE PRECISION, DIMENSION(:), POINTER :: ML    !mixed layer
     DOUBLE PRECISION, DIMENSION(:), POINTER :: all
     DOUBLE PRECISION :: div      !only for interior
     DOUBLE PRECISION :: h        !diffusivities at h
  END TYPE diff

  TYPE :: Knu               !K and nu diffusivities
     TYPE(diff)::u, t, s       
  END TYPE Knu

  TYPE :: flux              !Reynolds fluxes
     DOUBLE PRECISION, DIMENSION(:), POINTER::u, v, s, t, b, b_err
     TYPE(phyto)::p                                         
  END TYPE flux

  TYPE :: entrain           !Entrainment or mixed layer depth
     DOUBLE PRECISION :: depth
     INTEGER :: i,g
  END TYPE entrain

  TYPE :: height            !boundary layer depth 
     DOUBLE PRECISION :: old, new
     INTEGER :: i, g        !interface and grid index
     TYPE(entrain)::e, ml   !entrainment depth or mixed layer depth
  END TYPE height

  TYPE :: MST               ! Momentum, Salinity and Temp component vectors
     DOUBLE PRECISION, DIMENSION(:), POINTER::m,s,t    
  END TYPE MST

  TYPE :: MSTscalar
     DOUBLE PRECISION::m,s,t
  END TYPE MSTscalar

  TYPE :: boundary          !boundary conditions at h
     DOUBLE PRECISION :: h, div
     DOUBLE PRECISION, DIMENSION(:), POINTER::value
  END TYPE boundary

  TYPE :: MS                !Momemtum and scalar
     TYPE(boundary)::m, s
  END TYPE MS

  TYPE :: old_new
     DOUBLE PRECISION::old,new
  END TYPE old_new

  TYPE :: windstress
     TYPE(old_new)::u,v
  END TYPE windstress

  INTEGER::is_leap_year, was_leap_year


CONTAINS

  SUBROUTINE find_jmax_g(hh, d)
    use grid_mod, only: grid_

    TYPE(height), INTENT(IN OUT)::hh !mixed layer depth
    TYPE(grid_), INTENT(IN)::d 

    INTEGER::i


    DO i = 1, d%M
       IF (d%d_g(i) > hh%new) THEN
          hh%g = i
          EXIT
       ELSE
          hh%g= i+1
       END IF
    END DO

    IF (hh%g >= d%M-2)THEN
       PRINT "(A, I5)","Mixing too deep. find_jmax_g. hh%g>=grid&M-2", hh%g
       write (*,*) 'mean_param'
       STOP
    END IF

  END SUBROUTINE find_jmax_g

  SUBROUTINE find_jmax_i(h1, d1)
    use grid_mod, only: grid_

    TYPE(height), INTENT(IN OUT)::h1 !depth
    TYPE(grid_), INTENT(IN)::d1

    INTEGER::i

    DO i = 1, d1%M
       IF (d1%d_i(i) >= h1%new) THEN
          h1%i = i
          EXIT
       ELSE
          h1%i = i+1
       END IF
    END DO

    IF (h1%i >= d1%M-2)THEN
       PRINT "(A, I5)","Mixing too deep. find_jmax_i:", h1%i
       STOP
    END IF

  END SUBROUTINE find_jmax_i

END MODULE mean_param
