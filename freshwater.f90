! $Id$
! $Source$

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

  ! Variable Declarations:
  !
  ! Public:
  logical :: &
       Fw_surface  ! Add all of the fresh water on the surface?
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
       use_Fw_nutrients, &  ! Include influence of Fw nutrients?
       Northern_return ! include return flow from North?
  real(kind=dp), dimension(:), allocatable :: &
       Fw, &  ! Fresh water flux profile
       FN ! northern return flow profile
  real(kind=dp) :: &
       Fw_scale, &   ! Fresh water scale factor for river flows
       Fw_depth, &   ! Depth to distribute fresh water flux over
       upwell_const, &  ! Maximum upwelling velocity (tuning parameter)  
       Qbar, &          ! mean total freshwater 
       F_SOG, &          ! exponential of SOG component of Ft     
       F_RI, &           ! exponential of RI component of Ft     
       Fm, &             ! scale factor for RI river flow
       Northern_frac, & ! fraction of outgoing freshwater returned from the North
       phyto_sum, &
       cbottom, &       ! bottom salinity   
       ! values for salinity fit
       calpha, &        
       calpha2, &
       cgamma, &
       cbeta
contains

  subroutine init_freshwater(M)
    ! Allocate memory for fresh water quantity arrays, and read
    ! parameter values from the infile.
    
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    

    ! Allocate memory for fresh water quantity arrays
    call alloc_freshwater_variables(M)
    ! Read fresh water parameter values from the infile.
    call read_freshwater_params()

    phyto_sum = 3.d0 

  end subroutine init_freshwater


  subroutine read_freshwater_params()
    ! Read the fresh water parameter values from the infile.
    use input_processor, only: getpard, getparl
    implicit none

    ! Maximum upwelling velocity (tuning parameter)
    upwell_const = getpard("upwell_const")
    ! Fresh water scale factor for river flows
   ! Fresh water scale factor for river flows 
    Qbar = getpard("Qbar")
    F_SOG = getpard("F_SOG") 
    F_RI = getpard("F_RI")
    Fm = getpard("Fm")
    Fw_scale = getpard('Fw_scale')     
    ! Add all fresh water on surface?
    Fw_surface = getparl('Fw_surface') 
    ! Depth to distribute freshwater flux over
    if (.not. Fw_surface) then
       Fw_depth = getpard('Fw_depth')  
    endif
    ! Include effect of Fw nutrients?
    use_Fw_nutrients = getparl('use_Fw_nutrients')
    ! Include the return flow from the Northern Strait
    Northern_return = getparl('northern_return_flow_on')
    if (Northern_return) then
       Northern_frac = getpard('northern_return_flow_frac')
    end if 

    ! Values for salinity fit
    cbottom = getpard('cbottom')
    calpha = getpard('calpha')
    calpha2 = getpard('calpha2')
    cgamma = getpard('cgamma')
    cbeta = getpard('cbeta')
  end subroutine read_freshwater_params
  

  subroutine freshwater_phys(Qriver, Eriver, S_old, h)
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
   
    implicit none

    ! Arguments
    real(kind=dp), intent(in) :: &
         Qriver, &  ! Fraser River flow
         Eriver     ! Englishman River flow
    real(kind=dp), intent(in) :: &
         S_old,        &  ! Surface salinity
         h                ! Mixing layer depth

    ! Local variables
    real(kind=dp), parameter :: &
         Qmean = 2720.0d0 ! Mean fraser river flow from entrainment fit
    ! totalfresh water into system (Fraser + rest multiplied up from Englishman
    !real (kind=dp) :: totalfresh 

 

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
    Ft = Fw_scale*totalfresh*S_old/cbottom * (1-exp(-totalfresh/Fm))**F_RI * &
         ((totalfresh)/Qbar)**F_SOG
    
    ! Calculate the surface turbulent kinematic salinity flux (Large,
    ! et al (1994), eqn A2d).  Note that fresh water flux is added via
    ! Bf in calc_buoyancy().
    if (Fw_surface) then
       wbar%s(0) = Ft * S_old
    else
       wbar%s(0) = 0.0d0
    endif

    ! Calculate the nonturbulent fresh water flux profile, and its
    ! contribution to the salinity profile
    if (Fw_surface) then
       F_n = 0.0d0
    else
       Fw = Ft * exp(-grid%d_i / (Fw_depth * h))
       if (Northern_return) then
          FN = Northern_frac * 0.5d0 * Ft * &
               (1.d0 - tanh((grid%d_i - 15.d0)/(7.5d0)))
       else
          FN = 0.
       endif
       F_n = cbottom * (Fw + FN*0.d0)  ! remove salinity effect of return flow
    endif

  end subroutine freshwater_phys


  subroutine freshwater_bio(qty, current_value, surf_flux, distrib_flux)
    ! Calculate the freshwater biological fluxes.

    use grid_mod, only: grid
    use core_variables, only: S, P
    implicit none
    
    character (len=*), intent(in) :: qty
    real (kind=dp), dimension(0:), intent(in) :: current_value
    real (kind=dp), intent(out) :: surf_flux
    real (kind=dp), dimension(0:), intent(out) :: distrib_flux

    integer :: i
    real (kind=dp) :: use_value

    phyto_sum = phyto_sum*(1.d0-900./(30.*86400.)) + 900./(30*86400.) &
         * sum(P%micro(1:10))/10.

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
          else
             write (*,*) "problems in freshwater, river flux for ", qty, &
                  " is not defined."
             call exit(1)
          endif
       else ! distributed fresh water and fluxes
          surf_flux = 0.
          if (qty.eq."nitrate") then
             ! vector variation does not give the same results as the loop
             if (current_value(1).gt.4.5) then
                use_value = phyto_sum
             elseif (current_value(1).lt.3.5) then
                use_value = (30.5-1.0-current_value(1))
             else
                use_value = phyto_sum*(4.5-current_value(1)) +  &
                (current_value(1)-3.5)*(30.5-1.0-current_value(1))
             endif
             do i=0,grid%M
                ! if phyto are using the nitrate the grad is zero'd. 
                ! 4.0 = 2 x (Michalis-Menton K)
                distrib_flux(i) = Fw(i) * phys_circ_nitrate * &
                     min(4.0,current_value(i))/4.0 + &
                     FN(i) * 2 * use_value *30.d0/  &
                     (30.d0-S%new(0)) * &
                     min(4.0,current_value(i))/4.0 
!                if (current_value(1).lt.10.) then
!                   write (*,134) i, 0.5*i, Fw(i) * phys_circ_nitrate * &
!                     min(4.0,current_value(i))/4.0, &
!                     FN(i) * 2 * (30.0 - current_value(1)) *30.d0/  &
!                     (30.d0-20.d0) * &
!                     min(4.0,current_value(i))/4.0 
!0134 format (1x,i3,1x,f5.2,1x,e13.5,1x,e13.5)
!                endif
             enddo
             do i=grid%M,1,-1
                if (distrib_flux(i).gt.distrib_flux(i-1)) &
                     distrib_flux(i-1) = distrib_flux(i)
             enddo
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
    allocate(Fw(0:M), FN(0:M), F_n(0:M), &
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
    deallocate(Fw, FN, F_n, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_freshwater_variables

end module freshwater
