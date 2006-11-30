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

  ! circulation strength of a scalar is equal to deep value minus river value
  real(kind=dp), parameter :: phys_circ_nitrate = 30.5d0 - 13.0d0
  real(kind=dp), parameter :: phys_circ_Pmicro = 0.0d0
  real(kind=dp), parameter :: phys_circ_Pnano = 0.0d0
  real(kind=dp), parameter :: phys_circ_Zoo = 0.0d0
  real(kind=dp), parameter :: phys_circ_silicon = 54.0d0 - 80.0d0

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
    real(kind=dp), parameter :: &
         Qmean=2720.0d0 ! mean fraser river flow from entrainment fit

    ! Parameterized fit of the surface salinity of the Strait of
    ! Georgia at station S3 based on the river flows.  (Derived by
    ! Kate Collins 16-Jun-2005)  This value is not directly used
    ! in the model but is used to make sure the tuned FT value
    ! is correct.
    S_riv = 29.1166d0 - Qriver * (0.0019d0) - Eriver * (0.0392d0)
    
    ! tuned fresh water flux value (to give, on average) the parameterized
    ! value above.  
    Ft = Fw_scale * (0.0019d0 * Qriver + 0.0392d0 * Eriver) 
    
    ! The entrainment of deep water into the bottom of the
    ! grid is based on the parameterization derived by Susan Allen in
    ! Jun-2006 (See entrainment.pdf)
    upwell = upwell_const * (Qriver/Qmean)**0.41d0

    ! Salinity (eq'n A2d)
    ! Note that fresh water flux is added via Bf in buoyancy.f90
    ! *** Need to check the implications of w%s(0)=0 on def_gamma.f90
    if (Fw_surface) then
       w_s = Ft * S_o   ! w%s(0)
    else
       w_s = 0.0d0  ! w%s(0)
    endif
    
  end subroutine freshwater_phys

  subroutine freshwater_bio (qty, current_value, surf_flux, distrib_flux)
! calculates the freshwater biological fluxes

    use grid_mod, only: grid
    implicit none
    
    character (len=*), intent(in) :: qty
    real (kind=dp), dimension(0:), intent(in) :: current_value
    real (kind=dp), intent(out) :: surf_flux
    real (kind=dp), dimension(0:), intent(out) :: distrib_flux

    integer :: i

    if (Fw_surface) then
       distrib_flux = 0.
       if (qty.eq."nitrate") then
          ! if phyto are using the nitrate the grad is zero'd. 
          ! 4.0 = 2 x (Michalis-Menton K)
          surf_flux = Ft * phys_circ_nitrate * &
               min(4.0, current_value(1))/ 4.0
       elseif (qty.eq."Pmicro") then
          surf_flux = Ft * (phys_circ_Pmicro)
       elseif (qty.eq."Pnano") then
          surf_flux = Ft * (phys_circ_Pnano)
       elseif (qty.eq."Zoo") then
          surf_flux = Ft * (phys_circ_Zoo)
       elseif (qty.eq."silicon") then
          surf_flux = Ft * phys_circ_silicon * &
               min(4.0, current_value(1)) / 4.0
       else
          write (*,*) "problems in freshwater, river flux for ", qty, &
               " is not defined."
          call exit(1)
       endif
    else ! distributed fresh water and fluxes
       surf_flux = 0.
       if (qty.eq."nitrate") then
          ! vector variation does not give the same results as the loop
          do i=0,grid%M
             ! if phyto are using the nitrate the grad is zero'd. 
             ! 4.0 = 2 x (Michalis-Menton K)
             distrib_flux(i) = Fw(i) * phys_circ_nitrate * &
                  min(4.0,current_value(i))/4.0
          enddo
       elseif (qty.eq."Pmicro") then
          distrib_flux = Fw * (phys_circ_Pmicro)
       elseif (qty.eq."Pnano") then
          distrib_flux = Fw * (phys_circ_Pnano)
       elseif (qty.eq."Zoo") then
          distrib_flux = Fw * (phys_circ_Zoo)
       elseif (qty.eq."silicon") then
          do i=0,grid%M
             distrib_flux = Fw(i) * phys_circ_silicon * &
                  min(4.0, current_value(i))/4.0
          enddo
       else
          write (*,*) "problems in freshwater, river flux for ", qty, &
               " is not defined."
          call exit(1)
       endif
    endif
  end subroutine freshwater_bio

end module freshwater
