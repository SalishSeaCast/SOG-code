module freshwater

! subroutine to calculate the influences of the freshwater flow

  use precision_defs, only: dp, sp
!*** temporary will disappear as these variables become local here
  use declarations, only: Fw, Fw_surface, Ft, Fw_scale

  implicit none

  private
  ! public diagnostic variable
  public :: S_riv
  public :: freshwater_bio

  real(kind=dp), parameter :: river_nitrate = 6.4
  real(kind=dp), parameter :: river_Pmicro = 0.
  real(kind=dp), parameter :: river_Pnano = 0.
  real(kind=dp), parameter :: river_Zoo = 0.
  real(kind=dp), parameter :: river_silicon = 60.0

  real(kind=dp) :: S_riv ! surface salinity prediction from fit

contains

  subroutine init_freshwater
! reads in required parameters
  end subroutine init_freshwater

  subroutine freshwater_phys (Qriver, Eriver, upwell_const, S_o,& 
       upwell, w_s)
! calculates the freshwater flux Fw and the strength of upwelling/entrainment

    implicit none
    ! Arguments
    real(kind=sp), intent(in) :: Qriver, Eriver
    real(kind=dp), intent(in) :: upwell_const
    real(kind=dp), intent(in) :: S_o ! surface salinity
    real(kind=dp), intent(out) :: upwell
    real(kind=dp), intent(out) :: w_s ! surface salinity flux
    ! Local variables
    real(kind=dp), parameter :: Qmean=2720. ! mean fraser river flow from entrainment fit

    ! Parameterized fit of the surface salinity of the Strait of
    ! Georgia at station S3 based on the river flows.  (Derived by
    ! Kate Collins 16-Jun-2005)  This value is not directly used
    ! in the model but is used to make sure the tuned FT value
    ! is correct.
    S_riv = 29.1166 - Qriver * (0.0019) - Eriver * (0.0392)
    
    ! tuned fresh water flux value (to give, on average) the parameterized
    ! value above.  To keep the fit nice through river flow levels, the
    ! exponent here must match the upwell exponent.
    Ft = Fw_scale * (0.0019 * Qriver + 0.0392 * Eriver) * (Qriver/Qmean)**0.41
    
    ! The entrainment of deep water into the bottom of the
    ! grid is based on the parameterization derived by Susan Allen in
    ! Jun-2006 (See entrainment.pdf)
    upwell = upwell_const * (Qriver/Qmean)**0.41

    ! Salinity (eq'n A2d)
    ! Note that fresh water flux is added via Bf in buoyancy.f90
    ! *** Need to check the implications of w%s(0)=0 on def_gamma.f90
    if (Fw_surface) then
       w_s = Ft * S_o   ! w%s(0)
    else
       w_s = 0.  ! w%s(0)
    endif
    
  end subroutine freshwater_phys

  subroutine freshwater_bio (qty, current_value, surf_flux, distrib_flux)
! calculates the freshwater biological fluxes

    
    implicit none
    
    character (len=*), intent(in) :: qty
    real (kind=dp), dimension(0:), intent(in) :: current_value
    real (kind=dp), intent(out) :: surf_flux
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
          surf_flux = - Ft * (river_Zoo - current_value(1))
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
