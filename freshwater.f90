module freshwater

! subroutine to calculate the influences of the freshwater flow

  use precision_defs, only: dp

  implicit none

  private
  public :: freshwater_bio

  real(kind=dp), parameter :: river_nitrate = 6.4
  real(kind=dp), parameter :: river_Pmicro = 0.
  real(kind=dp), parameter :: river_Pnano = 0.

contains

  subroutine init_freshwater
! reads in required parameters
  end subroutine init_freshwater

  subroutine freshwater_phys
! calculates the freshwater flux Fw and the strength of upwelling/entrainment
  end subroutine freshwater_phys

  subroutine freshwater_bio (qty, current_value, aflux)
! calculates the freshwater biological fluxes

    use declarations, only: Fw
    
    implicit none
    
    character (len=*), intent(in) :: qty
    real (kind=dp), dimension(0:), intent(in) :: current_value
    real (kind=dp), dimension(0:), intent(out) :: aflux

    if (qty.eq."nitrate") then
       aflux = - Fw * (river_nitrate - current_value)
    elseif (qty.eq."Pmicro") then
       aflux = - Fw * (river_Pmicro - current_value)
    elseif (qty.eq."Pnano") then
       aflux = - Fw * (river_Pnano - current_value)
    else
       write (*,*) "problems in freshwater"
       stop
    endif

  end subroutine freshwater_bio

end module freshwater
