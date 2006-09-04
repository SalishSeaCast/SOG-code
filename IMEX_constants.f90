! $Id$
! $Source$

module IMEX_constants
  ! Parameters that control the types of differencing used by the IMEX
  ! scheme.  See Ascher et al (1995), SIAM J. Numer Anal. 32(3),
  ! pp. 797-823

  use precision_defs, only: dp
  implicit none

  ! Implicit-explicit first order scheme:
  !   a_IMEX1 = 1.0 --> backward Euler in diffusion
  !           = 0.5 --> Crank-Nicolson in diffusion
  !           = 0.0 --> forward Euler in diffusion
  ! Both have forward Euler in explicit terms (i.e. Coriolis and
  ! boundary conditions)
  real(kind=dp), parameter :: a_IMEX1 = 1.0

  ! Implicit-explicit second order scheme: 
  !   (a_IMEX2, b_IMEX2) = (0.5, 0) --> CNAB (Crank-Nicolson, Adams-Bashforth)
  !			 = (0.5, 0.125) --> MCNAB (modified Crank-Nicolson, 
  !                                                Adams-Bashforth)
  !                      = (0, 1)--> CNLF (Crank-Nicolson, Leap frog)
  !                      = (1, 0)--> SBDF (second order backward differencing)
  real(kind=dp), parameter :: a_IMEX2 = 1.0, b_IMEX2 =0.0
end module IMEX_constants
