module freshwater
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the fresh water flow influences SOG code.
  !
  ! Public Variables:
  !
  !   Fw_surface -- Add all of the fresh water on the surface?
  !
  !   Ft -- Total fresh water flux
  !
  !   F_n -- Fresh water contribution to salinity flux
  !
  ! Public Subroutines:
  !
  !   init_freshwater -- Allocate memory for fresh water quantity
  !                      arrays, and read parameter values from the
  !                      infile.
  !
  !   freshwater_phys -- Calculate the strength of
  !                      upwelling/entrainment, the freshwater flux,
  !                      the surface turbulent kinematic salinity
  !                      flux, and the profile of fresh water
  !                      contribution to the salinity flux.
  !
  !   freshwater_bio -- Calculate the freshwater biological fluxes.

  use precision_defs, only: dp

  implicit none

  private
  public :: &
       ! Variables:
       Fw_surface, &  ! Add all of the fresh water on the surface?
       Ft,         &  ! Total fresh water flux
       F_n,        &  ! Fresh water contribution to salinity flux
       upwell,     &  ! Upwelling velocity from river flows
                      ! parameterization.
       ! Diagnostics:
       S_riv, &  ! Surface salinity prediction from fit
       ! Subroutines:
       init_freshwater, freshwater_phys, freshwater_bio, &
       dalloc_freshwater_variables, &
       totalfresh

  ! Parameter Value Declarations:
  !
  ! Private:
  !
  ! Circulation strength of a scalar is equal to deep value minus
  ! river value
  real(kind=dp), parameter :: phys_circ_nitrate = 30.5d0 - 13.0d0
  real(kind=dp), parameter :: phys_circ_Pmicro = 0.0d0
  real(kind=dp), parameter :: phys_circ_Pnano = 0.0d0
  real(kind=dp), parameter :: phys_circ_Ppico = 0.0d0
  real(kind=dp), parameter :: phys_circ_Zoo = 0.0d0

  real(kind=dp), parameter :: phys_circ_silicon = 54.0d0 - 80.0d0
!--- BEGIN CHEMISTRY PARAMETERS
  ! From IOS cruise 2010-73 Oct 29th - Nov 2nd 2010
  real(kind=dp) :: phys_circ_DIC
  ! From Environment Canada Fraser River buoy October 2010
  real(kind=dp), parameter :: phys_circ_Oxy = 148.38d0 - 343.75d0
  ! **TODO**: Assign a sensible value for alkalinity
  real(kind=dp) :: phys_circ_Alk
!--- END CHEMISTRY PARAMETERS

  ! Variable Declarations:
  !
  ! Public:
  logical :: &
       Fw_surface   ! Add all of the fresh water on the surface?
  real(kind=dp) :: &
       Ft,  &  ! Total fresh water flux
       upwell  ! Upwelling velocity from river flows parameterization
  real(kind=dp), dimension(:), allocatable :: &
       F_n     ! Fresh water contribution to salinity flux
  real (kind=dp) :: totalfresh   !total freshwater -use to read in SOG

  ! Diagnostic:
  real(kind=dp) :: &
       S_riv  ! Surface salinity prediction from fit
  !
  ! Private:
  logical :: &
       use_Fw_nutrients  ! Include influence of Fw nutrients?
  integer :: &
       n_avg         ! Denominator for 30-day back-average
  integer, parameter :: &
       Ft_store_length = 30 * 86400 / 900
  real(kind=dp), dimension(:), allocatable :: &
       Fw    ! Fresh water flux profile
  real(kind=dp), dimension(Ft_store_length) :: &
       Ft_store    ! Ft storage vector
  real(kind=dp) :: &
       rho_riv,      &  ! Surface freshwater density [kg/m^3]
       Fresh_avg,    &  ! Running 30-day back-average
       Fw_scale,     &  ! Fresh water scale factor for river flows
       Fw_depth,     &  ! Depth to distribute fresh water flux over
       upwell_const, &  ! Maximum upwelling velocity (tuning parameter)  
       Qbar,         &  ! Mean total freshwater 
       F_SOG,        &  ! Exponential of SOG component of Ft     
       F_RI,         &  ! Exponential of RI component of Ft     
       cbottom,      &  ! Bottom salinity   
       ! Values for salinity fit
       calpha,       &        
       calpha2,      &
       cgamma,       &
       cbeta,        &
       ! Alkalinity fit parameters
       slope,        &  ! Slope of alkalinity versus discharge fit
       intercept,    &  ! Zero-discharge intercept of alk/discharge fit
       ! DIC fit parameters
       pCO2_riv         ! Freshwater pCO2 for DIC obtained from Alkalinity

contains

  subroutine init_freshwater(M)
    ! Allocate memory for fresh water quantity arrays, and read
    ! parameter values from the infile.
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points

    ! Initialize 30-day back average for alkalinity fit
    Ft_store = 0.0d0;
    n_avg = 1

    ! Allocate memory for fresh water quantity arrays
    call alloc_freshwater_variables(M)
    ! Read fresh water parameter values from the infile.
    call read_freshwater_params()

  end subroutine init_freshwater


  subroutine read_freshwater_params()
    ! Read the fresh water parameter values from the infile.
    use input_processor, only: getpard, getparl
    implicit none
    ! Maximum upwelling velocity (tuning parameter)
    upwell_const = getpard("upwell_const")
   ! Fresh water scale factor for river flows 
    Qbar = getpard("Qbar")
    F_SOG = getpard("F_SOG") 
    F_RI = getpard("F_RI")
    Fw_scale = getpard('Fw_scale')     
    ! Add all fresh water on surface?
    Fw_surface = getparl('Fw_surface') 
    ! Depth to distribute freshwater flux over
    if (.not. Fw_surface) then
       Fw_depth = getpard('Fw_depth')  
    endif
    ! Include effect of Fw nutrients?
    use_Fw_nutrients = getparl('use_Fw_nutrients')

    ! Values for salinity fit
    cbottom = getpard('cbottom')
    calpha = getpard('calpha')
    calpha2 = getpard('calpha2')
    cgamma = getpard('cgamma')
    cbeta = getpard('cbeta')

    ! Alkalinity fit parameters
    slope = getpard('slope_alk')
    intercept = getpard('intercept_alk')

    ! DIC fit parameters
    pCO2_riv = getpard('pCO2_river')
  end subroutine read_freshwater_params
  

  subroutine freshwater_phys(Qriver, Eriver, RiverTemp, S_old, Ts_old, Td_old, h)
    ! Calculate the strength of upwelling/entrainment, the freshwater
    ! flux, the surface turbulent kinematic salinity flux, and the
    ! profile of fresh water contribution to the salinity flux.

    ! Elements from other modules:
    !
    ! Type Definitions:
    use precision_defs, only: sp
    ! Variables:
    use grid_mod, only: &
         grid  ! Grid parameters and depth & spacing arrays
    use turbulence, only: &
         wbar  ! Turbulent kinematic flux profile arrays
    use irradiance, only: &
         Q_n   ! Non-turbulent heat flux profile array
    use forcing, only: &
         UseRiverTemp
    ! Functions and subroutines
    use unit_conversions, only: KtoC
    use carbonate, only: pCO2_to_carbonate
   
    implicit none

    ! Arguments
    real(kind=dp), intent(in) :: &
         Qriver, &  ! Fraser River flow
         Eriver, &  ! Englishman River flow
         RiverTemp  ! temperature of Major River
    real(kind=dp), intent(in) :: &
         S_old,        &  ! Surface salinity
         Ts_old,       &  ! Surface temperature
         Td_old,       &  ! Deep (bottom of grid) temperature        
         h                ! Mixing layer depth

    ! Local variables
    real(kind=dp) :: &
         RiverTC,    &    ! Major river temperature [deg C]
         river_alk,  &    ! River alkalinity [ueq L-1]
         river_DIC        ! River DIC [uM]

    ! Surface temperature to celsius
    RiverTC = KtoC(RiverTemp)

    ! fit to freshwater and entrainment pg 58-59, 29-Mar-2007
    totalfresh = Qriver + Eriver

    ! Parameterized fit of the surface salinity 

    ! For the  Strait of Georgia at station S3 based on the river flows.  
    !(Re-derived by S.Allen numerous times.  This time, late June 2007
    ! Labbook 125 to 131)

    ! For Rivers Inlet at Station DF02 based on river flows - derived by M.Wolfe, 
    ! March 2009 labbok March19th.
  
    ! This value is not 
    ! directly used in the model but is used to make sure the tuned Ft value
    ! is correct. tanh fn means no matter what the total fresh, a salinity
    ! below 0 is not predicted.  

    ! fit to freshwater and entrainment pg 58-59, 29-Mar-2007
    totalfresh = Qriver + 55.0*Eriver
    
    open(12,file="total_check")
    write(12,*)totalfresh, Qriver
    close(12)    

 
    S_riv = cbottom * ((exp(-totalfresh/calpha) + cbeta*exp(-totalfresh/calpha2)) / &
         ( cgamma + exp(-totalfresh/calpha) + cbeta*exp(-totalfresh/calpha2)))   

    open(12,file="S_riv_check")
    write(12,*)S_riv
    close(12)
 
    ! The entrainment of deep water into the bottom of the
    ! grid is based on the parameterization derived by Susan Allen in
    ! Jun-2006 (See entrainment.pdf) exponent with significant modifications
    ! in Mar-2007 (labbook pg 59)
    
    upwell = upwell_const * totalfresh/Qbar * exp(-totalfresh/Qbar) &
         /0.368
 
    ! Tuned fresh water flux value (to give, on average) the parameterized
    ! value above.
    Ft = Fw_scale * totalfresh * S_old/cbottom *  &
    ((1-0.35*exp(-(totalfresh-1.2*Qbar)**2/(4.*Qbar**2)))**F_RI * & !for SoG only
         (totalfresh/Qbar)**F_SOG) ! for both


!    Ft = Fw_scale * totalfresh * S_old/cbottom *  &
!    ((1-0.35*exp(-(totalfresh-1.2*Qbar)**2/(4.*Qbar**2)))*(totalfresh/Qbar)**0.2)**F_SOG * & ! for SoG
!    (1-exp(-totalfresh/Fm))**F_RI  ! for F_RI

    !-------------------------------------------------------------
    ! Calculate 30-day back average for alkalinity fit
    Ft_store(n_avg) = totalfresh
    Fresh_avg = sum(Ft_store) / n_avg
    if (n_avg .lt. Ft_store_length) then
       n_avg = n_avg + 1
    else
       Ft_store = cshift(Ft_store, 1)
    endif

    ! River density
    ! Density of pure water at p = 0 from Rich Pawlowicz limstate.m
    rho_riv = 999.842594d0 + RiverTC * ( 6.793952d-2      &
                           + RiverTC * (-9.095290d-3      &
                           + RiverTC * ( 1.001685d-4      &
                           + RiverTC * (-1.120083d-6      &
                           + RiverTC * ( 6.536332d-9)))))

    ! River alkalinity parametrization
    river_alk = intercept - slope * Fresh_avg

    ! River DIC parametrization
    call pCO2_to_carbonate(RiverTemp, 0.0d0, rho_riv, river_alk, &
         pCO2_riv, river_DIC)

    ! Circulation strength of a scalar is equal to deep value minus
    ! river value (see BMM Labbook pg 45 for deep values)
    phys_circ_Alk = 2092.98d0 - river_alk
    phys_circ_DIC = 2059.68d0 - river_DIC
    !-------------------------------------------------------------

    ! Calculate the surface turbulent kinematic salinity flux (Large,
    ! et al (1994), eqn A2d).  Note that fresh water flux is added via
    ! Bf in calc_buoyancy().
    if (Fw_surface) then
       wbar%s(0) = Ft * S_old
       if (UseRiverTemp) then
          wbar%t(0) = wbar%t(0) + Ft * (RiverTemp - Ts_old)
       endif
    else
       wbar%s(0) = 0.0d0
    endif

    ! Calculate the nonturbulent fresh water flux profile, and its
    ! contribution to the salinity profile
    if (Fw_surface) then
       F_n = 0.0d0
    else
       Fw = Ft * exp(-grid%d_i / (Fw_depth * h))
       F_n = cbottom * Fw  
       if (UseRiverTemp) then
          Q_n = Q_n + (RiverTemp - Td_old) * Fw  
       endif
    endif
  end subroutine freshwater_phys


  subroutine freshwater_bio(qty, current_value, surf_flux, distrib_flux)
    ! Calculate the freshwater biological fluxes.

    use grid_mod, only: grid
    implicit none
    
    character (len=*), intent(in) :: qty
    real (kind=dp), dimension(0:), intent(in) :: current_value
    real (kind=dp), intent(out) :: surf_flux
    real (kind=dp), dimension(0:), intent(out) :: distrib_flux

    integer :: i

    if (use_Fw_nutrients) then
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
          elseif (qty.eq."Ppico") then
             surf_flux = Ft * (phys_circ_Ppico)
          elseif (qty.eq."Zoo") then
             surf_flux = Ft * (phys_circ_Zoo)
          elseif (qty.eq."silicon") then
             surf_flux = Ft * phys_circ_silicon * &
                  min(4.0, current_value(1)) / 4.0
!--- BEGIN CHEMISTRY FRESHWATER FLUXES
          elseif (qty.eq."DIC") then
             surf_flux = Ft * (phys_circ_DIC)
          elseif (qty.eq."Oxy") then
             surf_flux = Ft * (phys_circ_Oxy)
          elseif (qty.eq."Alk") then
             surf_flux = Ft * (phys_circ_Alk)
!--- END CHEMISTRY FRESHWATER FLUXES
          else
             write (*,*) "problems in freshwater, river flux for ", qty, &
                  " is not defined."
             call exit(1)
          endif
       else ! distributed fresh water and fluxes
          surf_flux = 0.
          if (qty.eq."nitrate") then
             distrib_flux = Fw * (phys_circ_nitrate)
!             if (current_value(1).lt.10.) stop
          elseif (qty.eq."Pmicro") then
             distrib_flux = Fw * (phys_circ_Pmicro)
          elseif (qty.eq."Pnano") then
             distrib_flux = Fw * (phys_circ_Pnano)
          elseif (qty.eq."Ppico") then
             distrib_flux = Fw * (phys_circ_Ppico)
          elseif (qty.eq."Zoo") then
             distrib_flux = Fw * (phys_circ_Zoo)
          elseif (qty.eq."silicon") then
             do i=0,grid%M
                distrib_flux(i) = Fw(i) * phys_circ_silicon * &
                     min(4.0, current_value(i))/4.0
             enddo
!--- BEGIN CHEMISTRY FRESHWATER FLUXES
          elseif (qty.eq."DIC") then
             distrib_flux = Fw * (phys_circ_DIC)
          elseif (qty.eq."Oxy") then
             distrib_flux = Fw * (phys_circ_Oxy)
          elseif (qty.eq."Alk") then
             distrib_flux = Fw * (phys_circ_Alk)
!--- END CHEMISTRY FRESHWATER FLUXES
          else
             write (*,*) "problems in freshwater, river flux for ", qty, &
                  " is not defined."
             call exit(1)
          endif
       endif
    else
       distrib_flux = 0
       surf_flux = 0
    endif

  end subroutine freshwater_bio

  subroutine alloc_freshwater_variables(M)
    ! Allocate memory for fresh water variable arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Fresh water flux and salinity contribution profile arrays"
    allocate(Fw(0:M), F_n(0:M), &
         stat=allocstat) 
    call alloc_check(allocstat, msg)
  end subroutine alloc_freshwater_variables


  subroutine dalloc_freshwater_variables()
    ! Deallocate memory for fresh water variable arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg         ! Deallocation failure message prefix

    msg = "Fresh water flux and salinity contribution profile arrays"
    deallocate(Fw, F_n, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_freshwater_variables

end module freshwater
