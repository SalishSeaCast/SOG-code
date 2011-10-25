module biology_model
  ! Type definitions, variable & parameter value declarations, and
  ! subroutines related to the biolgy model in the SOG code.
  !
  ! Public Type Definitions:
  !
  !   cheese -- Description
  !
  ! Public Parameters:
  !
  !   spam -- Description
  !
  !
  ! Public subroutine:
  !
  !   calc_bio_rate -- A wrapper around a bunch of subroutine calls that
  !               advance the biological quantities to the next time step.

  implicit none

  private
  public :: &
       ! Type definitions:
!!$       cheese, &
       ! Parameter values:
!!$       spam, &
       ! Subroutines:
       init_biology, calc_bio_rate, dalloc_biology_variables

  ! Type Definitions:
  !
  ! Public:
  !
  ! Private to module:

  ! Parameter Value Declarations:
  !
  ! Public:
  !
  ! Private to module:
  integer, parameter :: &
       NPZD_bins = 9 + 6  ! Number of bins in biology model:
                          ! nutrients, phytoplankton, zooplankton, and
                          ! detritus (dissolved, slow sink, fast sink,
                          ! and silicon)
  !
  ! Variable Declarations:
  !
  ! Public:
  !
  ! Private to module:
  integer :: PZ_length

contains

  subroutine init_biology(M)
    ! Initialize the biology model.

    ! Subroutines from other modules:
    use NPZD, only: alloc_NPZD_variables, init_NPZD
    use biology_eqn_builder, only: alloc_bio_RHS_variables, read_sink_params
    implicit none

    ! Argument:
    integer, intent(in) :: &
         M  ! Number of grid points

    ! Length of PZ array
    PZ_length = NPZD_bins * M
    ! Allocate memory for NPZD model variables
    call alloc_NPZD_variables(M, PZ_length)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    call alloc_bio_RHS_variables(M)

!!$    call read_biology_params()
    call init_NPZD
    call read_sink_params
  end subroutine init_biology
  

!!$  subroutine read_biology_params
!!$    ! Read the biology model parameters from the input file
!!$    use input_processor, only: getpari, getpard
!!$    implicit none
!!$
!!$  end subroutine read_biology_params


  subroutine calc_bio_rate(time, day, dt, M, T_new, Pmicro, Pnano, Ppico, &
       Z, NO, NH, Si, DIC, Oxy, D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi)
    ! Solve the biology model ODEs to advance the biology quantity values
    ! to the next time step, and calculate the growth - mortality terms
    ! (*_RHS%bio) of the semi-implicit diffusion/advection equations.
    use precision_defs, only: dp
    use rungekutta, only: odeint
    use NPZD, only: PZ
    use numerics, only: check_negative
    use irradiance, only: &
         I_par  ! Photosynthetic available radiation profile
    implicit none

    ! Arguments:
    real(kind=dp), intent(in) :: time  ! Current model time
    integer, intent(in) :: day         ! Current model year-day
    real(kind=dp), intent(in) :: dt    ! Time step
    integer, intent(in) :: M           ! Number of grid points
    real(kind=dp), dimension(0:), intent(in) :: T_new   ! Temperature
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton
         Pnano,  &  ! Nano phytoplankton
         Ppico,  &  ! Pico phytoplankton
         Z,      &  ! Micro zooplankton
         NO,     &  ! Nitrate
         NH,     &  ! Ammonium
         Si,     &  ! Silicon
         DIC,    &  ! Dissolved inorganic carbon
         Oxy,    &  ! Dissolved oxygen
         D_DOC,  &  ! Dissolved organic carbon detritus profile
         D_POC,  &  ! Particulate organic carbon detritus profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    ! Local variables:
    integer :: N_ok, N_bad    ! counts bad and good steps in odeint
    real(kind=dp) :: next_time
    ! Biology model ODE solver control parameters
    real(kind=dp), parameter :: &
         precision = 1.0d-4,  &
         step_guess = 100.0d0,  &
         step_min = 3.0d0


    ! Load all of the biological quantities into the PZ vector for the
    ! ODE solver to operate on
    call load_PZ(M, Pmicro, Pnano, Ppico, Z, NO, NH, Si, DIC, Oxy, &
         D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, &
         PZ)                                               ! out
    call check_negative(1, PZ, "PZ after load_PZ()", day, time)

    ! Solve the biological model for values at the next time step
    next_time = time + dt
    call odeint(PZ, M, PZ_length, time, next_time, precision, step_guess, &
         step_min, &
         N_ok, N_bad, T_new, I_par, day)
    call check_negative(1, PZ, "PZ after odeint()", day, time)

    ! Unload the biological quantities from the PZ vector into the
    ! appropriate components of Gvector
    call unload_PZ(M, PZ, Pmicro, Pnano, Ppico, Z, NO, NH, Si, DIC, Oxy, &
         D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi)
  end subroutine calc_bio_rate


  subroutine load_PZ(M, Pmicro, Pnano, Ppico, Z, NO, NH, Si, DIC, Oxy, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi, &
       PZ)
    ! Load all of the separate biology variables (microplankton,
    ! nanoplankton, nitrate, ammonium, silicon, and detritus
    ! sequentially into the PZ vector for the ODE solver to use.
    use precision_defs, only: dp
    use NPZD, only: PZ_bins, flagellates, remineralization, microzooplankton
    implicit none
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton biomass profile array
         Pnano,  &  ! Nano phytoplankton biomass profile array
         Ppico,  &  ! Pico phytoplankton biomass profile array
         Z,      &  ! Micro zooplankton biomass profile array
         NO,     &  ! Nitrate concentration profile array
         NH,     &  ! Ammonium concentration profile array
         Si,     &  ! Silicon concentration profile array
         DIC,    &  ! Dissolved inorganic carbon profile array
         Oxy,    &  ! Dissolved oxygen profile array
         D_DOC,  &  ! Dissolved organic carbon detritus profile
         D_POC,  &  ! Particulate organic carbon detritus profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    real(kind=dp), dimension (1:), intent (out) :: &
         PZ  ! Array of biology variables for ODE solver
    ! Local Variables
    integer :: &
         bPZ, &  ! Beginning index for a quantity in the PZ array
         ePZ     ! Ending index for a quantity in the PZ array

    ! Initialize PZ array
    PZ = 0.
    ! Load micro phytoplankton
    bPz = (PZ_bins%micro - 1) * M + 1
    ePZ = PZ_bins%micro * M
    PZ(bPZ:ePZ) = Pmicro(1:M)
    ! Load nano phytoplankton
    bPz = (PZ_bins%nano - 1) * M + 1
    ePZ = PZ_bins%nano * M
    if (flagellates) then
       PZ(bPZ:ePZ) = Pnano(1:M)
    else
       PZ(bPZ:ePZ) = 0.
    endif
    ! Load pico phytoplankton
    bPz = (PZ_bins%pico - 1) * M + 1
    ePZ = PZ_bins%pico * M
    if (flagellates) then
       PZ(bPZ:ePZ) = Ppico(1:M)
    else
       PZ(bPZ:ePZ) = 0.
    endif
    ! Load micro zooplankton
    bPz = (PZ_bins%zoo - 1) * M + 1
    ePz = PZ_bins%zoo * M
    if (microzooplankton) then
       PZ(bPZ:ePz) = Z(1:M)
    else
       PZ(bPZ:ePz) = 0.
    endif
    ! Load nitrate
    bPz = (PZ_bins%NO - 1) * M + 1
    ePZ = PZ_bins%NO * M
    PZ(bPZ:ePZ) = NO(1:M)
    ! Load ammonium
    bPz = (PZ_bins%NH - 1) * M + 1
    ePZ = PZ_bins%NH * M
    if (remineralization) then
       PZ(bPZ:ePZ) = NH(1:M)
    else
       PZ(bPZ:ePZ) = 0
    endif
    ! Load silicon
    bPz = (PZ_bins%Si - 1) * M + 1
    ePZ = PZ_bins%Si * M
    PZ(bPZ:ePZ) = Si(1:M)
    ! Load dissolved inorganic carbon
    bPz = (PZ_bins%DIC - 1) * M + 1
    ePZ = PZ_bins%DIC * M
    PZ(bPZ:ePZ) = DIC(1:M)
   ! Load dissolved oxygen
    bPz = (PZ_bins%Oxy - 1) * M + 1
    ePZ = PZ_bins%Oxy * M
    PZ(bPZ:ePZ) = Oxy(1:M)
    ! Load dissolved organic carbon detritus
    bPz = (PZ_bins%D_DOC - 1) * M + 1
    ePz = PZ_bins%D_DOC * M
    if (remineralization) then
       PZ(bPz:ePz) = D_DOC(1:M)
    else
       PZ(bPz:ePz) = 0.
    endif
    ! Load particulate organic carbon detritus
    bPz = (PZ_bins%D_POC - 1) * M + 1
    ePz = PZ_bins%D_POC * M
    if (remineralization) then
       PZ(bPz:ePz) = D_POC(1:M)
    else
       PZ(bPz:ePz) = 0.
    endif
    ! Load dissolved organic nitrogen detritus
    bPz = (PZ_bins%D_DON - 1) * M + 1
    ePz = PZ_bins%D_DON * M
    if (remineralization) then
       PZ(bPz:ePz) = D_DON(1:M)
    else
       PZ(bPz:ePz) = 0.
    endif
    ! Load particulate organic nitrogen detritus
    bPz = (PZ_bins%D_PON - 1) * M + 1
    ePz = PZ_bins%D_PON * M
    if (remineralization) then
       PZ(bPz:ePz) = D_PON(1:M)
    else
       PZ(bPz:ePz) = 0.
    endif
    ! Load refractory nitrogen detritus
    bPz = (PZ_bins%D_refr - 1) * M + 1
    ePz = PZ_bins%D_refr * M
    if (remineralization) then
       PZ(bPz:ePz) = D_refr(1:M)
    else
       PZ(bPz:ePz) = 0.
    endif
    ! Load biogenic silicon detritus
    bPz = (PZ_bins%D_bSi - 1) * M + 1
    ePz = PZ_bins%D_bSi * M
    if (remineralization) then
       PZ(bPz:ePz) = D_bSi(1:M)
    else
       PZ(bPz:ePz) = 0.
    endif
  end subroutine load_PZ


  subroutine unload_PZ(M, PZ, Pmicro, Pnano, Ppico, Z, NO, NH, Si, DIC, Oxy, &
       D_DOC, D_POC, D_DON, D_PON, D_refr, D_bSi)
    ! Unload the biological quantities from the PZ vector into the
    ! appropriate *_RHS%bio arrays.
    use numerics, only: chem_steps
    use precision_defs, only: dp
    use NPZD, only: PZ_bins
    use biology_eqn_builder, only: Pmicro_RHS, Pnano_RHS, Ppico_RHS, Z_RHS, &
         NO_RHS, NH_RHS, Si_RHS, DIC_RHS, Oxy_RHS, &
         D_DOC_RHS, D_POC_RHS, D_DON_RHS, D_PON_RHS, D_refr_RHS, D_bSi_RHS
    implicit none
    ! Arguments:
    integer, intent(in) :: M  ! Number of grid points
    real(kind=dp), dimension(1:), intent(in) :: &
         PZ  ! Array of biology variables for ODE solver
    real(kind=dp), dimension(0:), intent(in) :: &
         Pmicro, &  ! Micro phytoplankton biomass profile array
         Pnano,  &  ! Nano phytoplankton biomass profile array
         Ppico,  &  ! Pico phytoplankton biomass profile array
         Z,      &  ! Micro zooplankton biomass profile array
         NO,     &  ! Nitrate concentration profile array
         NH,     &  ! Ammonium concentration profile array
         Si,     &  ! Silicon concentration profile array
         DIC,    &  ! Dissolved inorganic carbon profile array
         Oxy,    &  ! Dissolved oxygen profile array
         D_DOC,  &  ! Dissolved organic carbon detritus profile
         D_POC,  &  ! Particulate organic carbon detritus profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    ! Local Variables
    integer :: &
         bPZ, &  ! Beginning index for a quantity in the PZ array
         ePZ     ! Ending index for a quantity in the PZ array

    ! Unload micro phytoplankton
    bPZ = (PZ_bins%micro - 1) * M + 1
    ePZ = PZ_bins%micro * M
    Pmicro_RHS%bio = PZ(bPZ:ePZ) - Pmicro(1:M)
    ! Unload nano phytoplankton
    bPZ = (PZ_bins%nano - 1) * M + 1
    ePZ = PZ_bins%nano * M
    Pnano_RHS%bio = PZ(bPZ:ePZ) - Pnano(1:M)
    ! Unload pico phytoplankton
    bPZ = (PZ_bins%pico - 1) * M + 1
    ePZ = PZ_bins%pico * M
    Ppico_RHS%bio = PZ(bPZ:ePZ) - Ppico(1:M)
    ! Unlaod micro zooplankton
    bPZ = (PZ_bins%zoo - 1) * M + 1
    ePZ = PZ_bins%zoo * M
    Z_RHS%bio = PZ(bPZ:ePZ) - Z(1:M)
    ! Unload nitrate
    bPZ = (PZ_bins%NO - 1) * M + 1
    ePZ = PZ_bins%NO * M
    NO_RHS%bio = PZ(bPZ:ePZ) - NO(1:M)
    ! Unload ammonimum
    bPZ = (PZ_bins%NH - 1) * M + 1
    ePZ = PZ_bins%NH * M
    NH_RHS%bio = PZ(bPZ:ePZ) - NH(1:M)
    ! Unload silicon
    bPZ = (PZ_bins%Si - 1) * M + 1
    ePZ = PZ_bins%Si * M
    Si_RHS%bio = PZ(bPZ:ePZ) - Si(1:M)
    ! Unload dissolved inorganic carbon
    bPZ = (PZ_bins%DIC - 1) * M + 1
    ePZ = PZ_bins%DIC * M
    DIC_RHS%bio = (PZ(bPZ:ePZ) - DIC(1:M))/chem_steps
    ! Unload dissolved oxygen
    bPZ = (PZ_bins%Oxy - 1) * M + 1
    ePZ = PZ_bins%Oxy * M
    Oxy_RHS%bio = (PZ(bPZ:ePZ) - Oxy(1:M))/chem_steps
    ! Unload dissolved organic carbon detritus
    bPz = (PZ_bins%D_DOC - 1) * M + 1
    ePz = PZ_bins%D_DOC * M
    D_DOC_RHS%bio = PZ(bPz:ePz) - D_DOC(1:M)
    ! Unload particulate organic carbon detritus
    bPz = (PZ_bins%D_POC - 1) * M + 1
    ePz = PZ_bins%D_POC * M
    D_POC_RHS%bio = PZ(bPz:ePz) - D_POC(1:M)
    ! Unload dissolved organic nitrogen detritus
    bPz = (PZ_bins%D_DON - 1) * M + 1
    ePz = PZ_bins%D_DON * M
    D_DON_RHS%bio = PZ(bPz:ePz) - D_DON(1:M)
    ! Unload particulate organic nitrogen detritus
    bPz = (PZ_bins%D_PON - 1) * M + 1
    ePz = PZ_bins%D_PON * M
    D_PON_RHS%bio = PZ(bPz:ePz) - D_PON(1:M)
    ! Unload refractory nitrogen detritus
    bPz = (PZ_bins%D_refr - 1) * M + 1
    ePz = PZ_bins%D_refr * M
    D_refr_RHS%bio = PZ(bPz:ePz) - D_refr(1:M)
    ! Unload biogenic silicon detritus
    bPz = (PZ_bins%D_bSi - 1) * M + 1
    ePz = PZ_bins%D_bSi * M
    D_bSi_RHS%bio = PZ(bPz:ePz) - D_bSi(1:M)
  end subroutine unload_PZ


  subroutine dalloc_biology_variables
    ! Deallocate memory from physics model variables arrays.

    ! Subroutines from other modules:
    use malloc, only: dalloc_check
    use NPZD, only: dalloc_NPZD_variables
    use biology_eqn_builder, only: dalloc_bio_RHS_variables
    
    implicit none

    ! Deallocate memory from NPZD model arrays
    call dalloc_NPZD_variables
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    call dalloc_bio_RHS_variables
  end subroutine dalloc_biology_variables

end module biology_model
