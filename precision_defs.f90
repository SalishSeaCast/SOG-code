module precision_defs
  ! Define numerical precision parameters for intrinsic types.
  ! Parameters defined here are for use as kind arguments in 
  ! declarations of intrinsic types.
  ! Note that x86, AMD, and PowerPC processors support only limited
  ! choices.

  implicit none

  ! precision = 6, range = 37
  ! Equivalent to the old single precision for reals
  integer, parameter :: sp = selected_real_kind(6, 37)
  ! precision = 15, range = 307
  ! Equivalent to the old double precision for reals
  integer, parameter :: dp = selected_real_kind(15, 307)

end module precision_defs
