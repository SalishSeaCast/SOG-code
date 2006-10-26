! $Id$
! $Source$

module biology_eqn_builder
  ! Type definitions, variable declarations, and subroutines related
  ! to building the semi-implicit diffusion/advection equations for
  ! the biology quantities in the SOG code.
  !
  ! Public Variables:
  !
  !   Pmicro_RHS -- Micro phytoplankton (diatoms) right-hand side arrays
  !
  !   Pnano_RHS -- Nano phytoplankton (flagellates) right-hand side arrays
  !
  !   NO_RHS -- Nitrate concentration right-hand side arrays
  !
  !   NH_RHS -- Ammonium concentration right-hand side arrays
  !
  !   Si_RHS -- Silicon concentration right-hand side arrays
  !
  !   D_DON_RHS -- Dissolved organic nitrogen detritus right-hand side arrays
  !
  !   D_PON_RHS -- Particulate organic nitro detritus right-hand side arrays
  !
  !   D_refr_RHS -- Refractory nitrogen detritus right-hand side arrays
  !
  !   D_bSi_RHS -- Biogenic silicon detritus right-hand side arrays
  !
  ! Public Subroutines:
  !
  !   alloc_bio_RHS_variables -- Allocate memory for biology RHS arrays
  !
  !   dalloc_bio_RHS_variables -- Deallocate memory for biology RHS arrays

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       Pmicro_RHS, &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,  &  ! Nano phytoplankton (flagellates) RHS arrays
       NO_RHS,     &  ! Nitrate concentration RHS arrays
       NH_RHS,     &  ! Ammonium concentration RHS arrays
       Si_RHS,     &  ! Silicon concentration RHS arrays
       D_DON_RHS,  &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,  &  ! Particulate organic nitrogen detritus RHS arrays
       D_refr_RHS, &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS,  &  ! Biogenic silicon detritus RHS arrays
       ! Subroutines:
       alloc_bio_RHS_variables, dalloc_bio_RHS_variables

  ! Type Definitions:
  !
  ! Private to module:
  !
  ! New/old array components:
  type :: new_old
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Current time step values
          old     ! Previous time step values
  end type new_old
  !
  ! Semi-implicit diffusion/advection equation right-hand side (RHS) arrays
  type :: RHS
     type(new_old) :: &
          diff_adv  ! Diffusion/advection component of RHS
     real(kind=dp), dimension(:), allocatable :: &
          bio, &  ! Biology (growth - mortality) component of RHS
          sink    ! Sinking component of RHS
  end type RHS

  ! Parameter Value Declarations:
  !
  ! Public:
  type(RHS) :: &
       Pmicro_RHS, &  ! Micro phytoplankton (diatoms) RHS arrays
       Pnano_RHS,  &  ! Nano phytoplankton (flagellates) RHS arrays
       NO_RHS,     &  ! Nitrate concentration RHS arrays
       NH_RHS,     &  ! Ammonium concentration RHS arrays
       Si_RHS,     &  ! Silicon concentration RHS arrays
       D_DON_RHS,  &  ! Dissolved organic nitrogen detritus RHS arrays
       D_PON_RHS,  &  ! Particulate organic nitro detritus RHS arrays
       D_refr_RHS, &  ! Refractory nitrogen detritus RHS arrays
       D_bSi_RHS      ! Biogenic silicon detritus RHS arrays
  !
  ! Private to module:

contains

  subroutine build_biology_equations()
    ! Build the right-hand side (RHS) arrays for the
    ! diffusion/advection equations for the biology quantities.

    ! Initialize the RHS *%bio%new arrays, and calculate the diffusive
    ! fluxes at the bottom and top of the grid
    

  end subroutine build_biology_equations


  subroutine alloc_bio_RHS_variables(M)
    ! Allocate memory for arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Micro phytoplankton RHS arrays"
    allocate(Pmicro_RHS%diff_adv%new(1:M), Pmicro_RHS%diff_adv%old(1:M), &
         Pmicro_RHS%bio(1:M), Pmicro_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    allocate(Pnano_RHS%diff_adv%new(1:M), Pnano_RHS%diff_adv%old(1:M), &
         Pnano_RHS%bio(1:M), Pnano_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    allocate(NO_RHS%diff_adv%new(1:M), NO_RHS%diff_adv%old(1:M), &
         NO_RHS%bio(1:M), NO_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    allocate(NH_RHS%diff_adv%new(1:M), NH_RHS%diff_adv%old(1:M), &
         NH_RHS%bio(1:M), NH_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Silicon concentration RHS arrays"
    allocate(Si_RHS%diff_adv%new(1:M), Si_RHS%diff_adv%old(1:M), &
         Si_RHS%bio(1:M), Si_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    allocate(D_DON_RHS%diff_adv%new(1:M), D_DON_RHS%diff_adv%old(1:M), &
         D_DON_RHS%bio(1:M), D_DON_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    allocate(D_PON_RHS%diff_adv%new(1:M), D_PON_RHS%diff_adv%old(1:M), &
         D_PON_RHS%bio(1:M), D_PON_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    allocate(D_refr_RHS%diff_adv%new(1:M), D_refr_RHS%diff_adv%old(1:M), &
         D_refr_RHS%bio(1:M), D_refr_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    allocate(D_bSi_RHS%diff_adv%new(1:M), D_bSi_RHS%diff_adv%old(1:M), &
         D_bSi_RHS%bio(1:M), D_bSi_RHS%sink(1:M), &
         stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_bio_RHS_variables


  subroutine dalloc_bio_RHS_variables
    ! Deallocate memory from arrays for right-hand sides of
    ! diffusion/advection equations for the biology model.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Allocation return status
    character(len=80) :: msg         ! Allocation failure message prefix

    msg = "Micro phytoplankton RHS arrays"
    deallocate(Pmicro_RHS%diff_adv%new, Pmicro_RHS%diff_adv%old, &
         Pmicro_RHS%bio, Pmicro_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nano phytoplankton RHS arrays"
    deallocate(Pnano_RHS%diff_adv%new, Pnano_RHS%diff_adv%old, &
         Pnano_RHS%bio, Pnano_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrate concentration RHS arrays"
    deallocate(NO_RHS%diff_adv%new, NO_RHS%diff_adv%old, &
         NO_RHS%bio, NO_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ammonium concentration RHS arrays"
    deallocate(NH_RHS%diff_adv%new, NH_RHS%diff_adv%old, &
         NH_RHS%bio, NH_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Silicon concentration RHS arrays"
    deallocate(Si_RHS%diff_adv%new, Si_RHS%diff_adv%old, &
         Si_RHS%bio, Si_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved organic nitrogen detritus RHS arrays"
    deallocate(D_DON_RHS%diff_adv%new, D_DON_RHS%diff_adv%old, &
         D_DON_RHS%bio, D_DON_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Particulate organic nitrogen detritus RHS arrays"
    deallocate(D_PON_RHS%diff_adv%new, D_PON_RHS%diff_adv%old, &
         D_PON_RHS%bio, D_PON_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Refractory nitrogen detritus RHS arrays"
    deallocate(D_refr_RHS%diff_adv%new, D_refr_RHS%diff_adv%old, &
         D_refr_RHS%bio, D_refr_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Biogenic silicon detritus RHS arrays"
    deallocate(D_bSi_RHS%diff_adv%new, D_bSi_RHS%diff_adv%old, &
         D_bSi_RHS%bio, D_bSi_RHS%sink, &
         stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_bio_RHS_variables

end module biology_eqn_builder
