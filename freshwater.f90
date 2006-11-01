module freshwater

! subroutine to calculate the influences of the freshwater flow

  use precision_defs, only: dp

  implicit none

  private
  public :: freshwater_bio

  real(kind=dp), parameter :: river_nitrate = 6.4
  real(kind=dp), parameter :: river_Pmicro = 0.
  real(kind=dp), parameter :: river_Pnano = 0.
  real(kind=dp), parameter :: river_Zoo = 0.
  real(kind=dp), parameter :: river_silicon = 60.0

contains

  subroutine init_freshwater
! reads in required parameters
  end subroutine init_freshwater

  subroutine freshwater_phys
! calculates the freshwater flux Fw and the strength of upwelling/entrainment
  end subroutine freshwater_phys

  subroutine freshwater_bio (qty, current_value, surf_flux, distrib_flux)
! calculates the freshwater biological fluxes

    use declarations, only: Fw, Fw_surface 
    
    implicit none
    
    character (len=*), intent(in) :: qty
    real (kind=dp), dimension(0:), intent(in) :: current_value
    real (kind=dp), dimension(0:), intent(out) :: distrib_flux

    if (Fw_surface) then
       distrib_flux = 0.
       if (qty.eq."nitrate") then
          surf_flux = - Ft * (river_nitrate - current_value(1))
       elseif (qty.eq."Pmicro") then
          surf_flux = - Ft * (river_Pmicro - current_value(1))
       elseif (qty.eq."Pnano") then
          surf_flux = - Ft * (river_Pnano - current_value(1))
       elseif (qty.eq."Zoo") then
          surf_flux = - Ft * (river_Zoo - current_value(1)
       elseif (qty.eq."silicon") then
          surf_flux = - Ft * (river_silicon - current_value(1))
       else
          write (*,*) "problems in freshwater, river flux for ", qty, &
               " is not defined."
          stop
       endif
    else ! distributed fresh water and fluxes
       surf_flux = 0.
       if (qty.eq."nitrate") then
          distrib_flux = - Fw * (river_nitrate - current_value)
       elseif (qty.eq."Pmicro") then
          distrib_flux = - Fw * (river_Pmicro - current_value)
       elseif (qty.eq."Pnano") then
          distrib_flux = - Fw * (river_Pnano - current_value)
       elseif (qty.eq."Zoo") then
          distrib_flux = - Fw * (river_Zoo - current_value)
       elseif (qty.eq."silicon") then
          distrib_flux = - Fw * (river_silicon - current_value)
       else
          write (*,*) "problems in freshwater, river flux for ", qty, &
               " is not defined."
          stop
       endif
    endif
  end subroutine freshwater_bio

end module freshwater
