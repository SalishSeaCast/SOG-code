module core_variables
  ! Type definitions, variable declarationss, and subroutines related
  ! to the core variables that the SOG code calculates.

  use precision_defs, only: dp
  implicit none

  private
  public :: &
       ! Variables:
       U,   &  ! Cross-estuary velocity component profile arrays [m/s]
       V,   &  ! Along-estuary velocity component profile arrays [m/s]
       T,   &  ! Water column temperature profile arrays [K]
       S,   &  ! Water column Salinity profile arrays [-]
       P,   &  ! Phytoplankton biomass profile arrays [uM N]
               !   P%micro: Diatoms
               !   P%nano: Mesodinium rubrum
               !   P%pico: Flagellates
       Z,   &  ! Microzooplankton biomass profile array [uM N]
       N,   &  ! Nitrogen compounds concentation profile arrays [uM N]
               !   N%O: Nitrate
               !   N%H: Ammonium
       Si,  &  ! Silicon concentration profile array [uM Si]
       DIC, &  ! Dissolved inorganic carbon concentration profile array [uM C]
       Oxy, &  ! Dissolved oxygen concentration profile array [uM O]
       Alk, &  ! Alkalinity profile array [uM]
       D,   &  ! Detritus concentration profile arrays
               !   D%DON -- Dissolved organic nitrogen [uM N]
               !   D%PON -- Particulate organic nitrogen [uM N]
               !   D%refr -- Refractory nitrogen [uM N]
               !   D%bSi -- Biogenic silicon [uM Si]
               !   D%DOC -- Dissolved organic carbon [uM C]
               !   D%POC -- Particulate organic carbon [uM C]
               !   D%reC -- Refractory carbon [uM C]

       ! Subroutines:
       init_core_variables,   &  ! Allocate memory for core variable
                                 ! arrays, and set their initial values
                                 ! to assumed values, and values
                                 ! interpolated from field data (CTD,
                                 ! STRATOGEM bottles, IOS bottles).
       dalloc_core_variables, &  ! De-allocate memory for core
                                 ! variables arrays.

       ! Parameter values:
       N2chl, &  ! ratio of chl mg/m3 to uMol N in phytoplankton
       ho,    &  ! Initial mixing layer depth [m]
       biology   ! No phytoplankton: test just the physics

  ! Public parameter declaration:
  real(kind=dp), parameter :: &
       ho = 2.0d0        ! Initial mixing layer depth [m]
  real(kind=dp) :: &
       N2chl    ! ratio of chlorophyl mg/m3 to uMol N in phytoplankton
  logical :: biology  ! No phytoplankton: test just the physics

  ! Public type definitions:
  !
  ! Velocities, temperature, and salinity
  type :: profiles
     real(kind=dp), dimension(:), allocatable :: &
          new, &  ! Profile of quantity at current time setp
          old, &  ! Profile of quantity at previous time step
          grad_i  ! Profile of gradient of quantity at grid layer interfaces
  end type profiles
  !
  ! Nitrogen compounds
  type :: nitrogen
     real(kind=dp), dimension(:), allocatable :: &
          O, &  ! N%O is nitrate (NO3) concentration profile
          H     ! N%H is ammonium (NH4) concentration profile
  end type nitrogen
  !
  ! Plankton
  type :: plankton
     real(kind=dp), dimension(:), allocatable :: &
          micro, &  ! P%micro is micro phytoplankton (diatoms) biomass profile
          nano, &   ! P%nano is nano phytoplankton (meso-rub) biomass profile
          pico      ! P%pico is pico phytoplankton (flagellate) biomass profile
  end type plankton
  !
  ! Detritus
  type :: detritus
     real(kind=dp), dimension(:), allocatable :: &
          DON,  &  ! Dissolved organic nitrogen [uM N]
          PON,  &  ! Particulate organic nitrogen [uM N]
          refr, &  ! Refractory nitrogen [uM N]
          bSi,  &  ! Biogenic silicon [uM Si]
          DOC,  &  ! Dissolved organic carbon [uM C]
          POC,  &  ! Particulate organic carbon [uM C]
          reC      ! Refractory carbon [uM C]
  end type detritus
  !
  ! Public variable declarations:
  type(profiles) :: &
       U, &  ! Cross-estuary velocity component arrays
       V, &  ! Along-estuary velocity component arrays
       T, &  ! Temperature profile arrays
       S     ! Salinity profile arrays
  type(plankton) :: &
       P  ! Micro & nano & pico phytoplankton biomass profile arrays
  real(kind=dp), dimension(:), allocatable :: &
       Z  ! Microzooplankton concentration profile array
  type(nitrogen) :: &
       N  ! Nitrate & ammonium concentration profile arrays
  real(kind=dp), dimension(:), allocatable :: &
       Si,  &  ! Silicon concentration profile array
       DIC, &  ! Dissolved inorganic carbon concentration profile array
       Oxy, &  ! Dissolved oxygen concentration profile array
       Alk     ! Alkalinity profile array
  type(detritus) :: &
       D  ! Detritus concentration profile arrays


  ! Private parameter declarations:
  real(kind=dp), parameter :: &
       Uo = 0.0d0,  &  ! Initial cross-estuary velocity component
       Vo = 0.0d0,  &  ! Initial along-estuary velocity component
       NHo = 1.0d0     ! Estimate of deep ammonium from
                       ! nitrate/salinity fits gives deep nitrate as
                       ! 31.5 uM but we measure 30.5 --- difference
                       ! must be remineralized NH
  !
  ! Private type definitions:
  type :: col_indices
     integer :: &
          depth,   &  ! Depth column number in CTD, Nuts or bottle data file
          T,       &  ! Temperature column number in CTD, Nuts or
                      ! bottle data file
          S,       &  ! Salinity column number in CTD, Nuts or
                      ! bottle data file
          Chloro,  &  ! Chlorophyl column number in CTD, Nuts or
                      ! bottle data file
          Fluores, &  ! Fluorescence column number in CTD, Nuts or
                      ! bottle data file
          NO,      &  ! Nitrate column number in CTD, Nuts or bottle
                      ! or Nuts data file
          Phyto,   &  ! Phytoplankton ?? column number in CTD, Nuts
                      ! or bottle data file
          Si,      &  ! Silicon column number in Nuts data file
          DIC,     &  ! Dissolved inorganic carbon column number in
                      ! bottle data file
          Alk,     &  ! Total alkalinity column number in
                      ! bottle data file
          Oxy         ! Dissolved oxygen column number in CTD
  end type col_indices


contains

  subroutine init_core_variables()
    ! Allocate memory for core variable arrays, and set their initial
    ! values to assumed values, and values interpolated from field data
    ! (CTD, STRATOGEM bottles, IOS bottles).
    use grid_mod, only: grid
    implicit none


    ! Allocate memory for core variable arrays
    call alloc_core_variables(grid%M)
    ! Initialize the values of the core variable profiles
    call init_state()
  end subroutine init_core_variables


  subroutine init_state()
    ! Initialize the values of the core variable profiles from assumed
    ! values and interpolated field data (CTD, STRATOGEM bottles, IOS
    ! bottles).
    use io_unit_defs, only: stdout
    use grid_mod, only: grid
    use input_processor, only: getpars, getpardv, getpard, getparl
    use unit_conversions, only: CtoK
    use forcing, only: vary
    use fundamental_constants, only: Redfield_C
    implicit none
    ! Local variables:
    character*80 :: fn  ! name of data file to read
    integer :: &
         i,            &  ! loop index
         ctd_records,  &  ! CTD data record counter
         nuts_records, &  ! STRATOGEM bottle data (Nuts*.txt) record counter
         botl_records, &  ! IOS bottle data record counter
         chem_records     ! IOS chem bottle data record counter
    logical :: got_Fluores, got_NO, got_Si, got_DIC, got_Alk, got_Oxy
    real(kind=dp), dimension(0:3*int(grid%M+1), 24) :: &

         data  ! Data records read
    real(kind=dp), dimension(3) :: &
         Psplit  ! Initial ratios of phytoplankton classes (micro, nano, pico)
    type(col_indices) :: col  ! Column numbers of quantities in data records

    ! Initializes N2chl ratio
    N2chl = getpard('N2chl')

    ! Biology 'on' or 'off'
    biology = getparl('biology')

    ! Initialize velocity profiles to assumed values
    U%new = Uo
    V%new = Vo


    ! Initialize the ammonium profile to the assumed deep water value
    ! except in the mixed layer where it is assumed to be zero.
    ! **Note: that assumption is not valid for runs starting in the
    ! winter.**
    do i = 0, grid%M + 1
       if(grid%d_g(i) <= ho) then
          N%H(i) = 0.0d0
       else
          N%H(i) = NHo
       endif
    enddo

    ! Read the CTD data file to get data to initialize the
    ! temperature, and salinity profiles, and also the
    ! microphytoplankton biomass and nitrate profiles if fluorescence
    ! and nitrate data are included in the file
    got_Fluores = .false.
    got_NO = .false.
    got_Oxy = .false.
    got_DIC = .false.
    got_Alk = .false.
    call read_init_data(getpars("ctd_in"), ctd_records, col, data)
    ! Interpolate CTD data to grid point depths
    if(col%T /= -1) then
       T%new = interp_to_grid(data(:ctd_records, col%depth), &
                              data(:ctd_records, col%T))
       ! Convert temperature to Kelvin, and apply variation, if
       ! enabled.
       ! *** TODO: Refactor initial temperature profile variation out
       ! *** of this module.
       if (vary%temperature%enabled .and. .not. vary%temperature%fixed) then
          T%new = CtoK(T%new + vary%temperature%addition)
       else
          T%new = CtoK(T%new)
       endif
    else


       ! Run can't proceed without temperature data
       write(stdout, *) "init_state: No temperature data found. Run aborted."
       call exit(1)
    endif
    ! Salinity
    if(col%S /= -1) then
       S%new = interp_to_grid(data(:ctd_records, col%depth), &
                              data(:ctd_records, col%S))
    else
       ! Run can't proceed without salinity data
       write(stdout, *) "init_state: No salinity data found. Run aborted."
       call exit(1)
    endif

    if(col%Fluores /= -1) then
       ! Microphytoplankton biomass comes from fluorescence if that
       ! data are in the CTD data file, otherwise it comes from the
       ! IOS bottle data file (typically only for historical runs)
       P%micro = interp_to_grid(data(:ctd_records, col%depth), &
                               data(:ctd_records, col%Fluores))
       got_Fluores = .true.
    endif

    if(col%NO /= -1) then
       ! Nitrate profile comes from CTD profile if the data are
       ! available, otherwise it comes from the STRATOGEM or IOS
       ! bottle data file
       N%O = interp_to_grid(data(:ctd_records, col%depth), &
                            data(:ctd_records, col%NO))
       got_NO = .true.
    endif

    if(col%Oxy /= -1) then
       ! Dissolved oxygen profile comes from CTD profile if the
       ! data are available, otherwise it comes from the IOS
       ! bottle data file
       ! 44.660 umol/ml O2 unit conversion from SeaBrid electronics
       Oxy = interp_to_grid(data(:ctd_records, col%depth), &
                            data(:ctd_records, col%Oxy) * 44.660d0)
       got_Oxy = .true.
    endif

    ! If a STRATOGEM bottle data file is specified, read it to get
    ! data to initialize the silicon profile, and the nitrate profile
    ! if those data weren't in the CTD data file
    got_Si = .false.
    fn = getpars("nuts_in")
    if(fn /= "N/A") then
       call read_init_data(fn, nuts_records, col, data)
       ! Nitrate
       if(col%NO /= -1) then
          N%O = interp_to_grid(data(:nuts_records, col%depth), &
                               data(:nuts_records, col%NO))
          got_NO = .true.
       endif
       ! Silicon
       if(col%Si /= -1) then
          Si = interp_to_grid(data(:nuts_records, col%depth), &
                              data(:nuts_records, col%Si))
          got_Si = .true.
       endif
    endif

    ! If an IOS bottle data file is specfied, read it to get data to
    ! initialize the nitrate, microphytoplankton, and silicon
    ! profiles, if those data weren't in the CTD data file.
    fn = getpars("botl_in")
    if(fn /= "N/A") then
       call read_init_data(fn, botl_records, col, data)
       ! Nitrate
       if(col%NO /= -1) then
          N%O = interp_to_grid(data(:botl_records, col%depth), &
                               data(:botl_records, col%NO))
          got_NO = .true.
       endif
       ! Silicon
       if(col%Si /= -1) then
          Si = interp_to_grid(data(:botl_records, col%depth), &
                              data(:botl_records, col%Si))
          got_Si = .true.
       endif
       ! Phytoplankton
       if(.not. got_Fluores) then
          ! First choice is fluorescence data
          if(col%Fluores /= -1) then
             P%micro = interp_to_grid(data(:botl_records, col%depth), &
                                      data(:botl_records, col%Fluores))
             got_Fluores = .true.
          elseif(col%Chloro /= -1) then
             ! Second choice is chlorophyl data
             P%micro = interp_to_grid(data(:botl_records, col%depth), &
                                      data(:botl_records, col%Chloro))
             got_Fluores = .true.
          elseif(col%Phyto /= -1) then
             ! Third choice is phytoplankton data
             P%micro = interp_to_grid(data(:botl_records, col%depth), &
                                      data(:botl_records, col%Phyto))
             got_Fluores = .true.
          endif
       endif
    endif

    ! If an IOS chemistry bottle data file is specfied, read it to
    ! get data to initialize the dissolved inorganic carbon and dissolved
    ! oxygen profiles, if those data weren't in the CTD data file.
    fn = getpars("chem_in")
    if(fn /= "N/A") then
       call read_init_data(fn, chem_records, col, data)
       ! Dissolved inorganic carbon
       if(col%DIC /= -1) then
          DIC = interp_to_grid(data(:chem_records, col%depth), &
                               data(:chem_records, col%DIC))
          got_DIC = .true.
       endif
       ! Total alkalinity
       if(col%Alk /= -1) then
          Alk = interp_to_grid(data(:chem_records, col%depth), &
                               data(:chem_records, col%Alk))
          got_Alk = .true.
       endif
       ! Dissolved oxygen
       if(col%Oxy /= -1) then
          Oxy = interp_to_grid(data(:chem_records, col%depth), &
                               data(:chem_records, col%Oxy))
          got_Oxy = .true.
       endif
    endif

    ! Run can't proceed without nitrate data
    if(.not. got_NO) then
       write(stdout, *) "init_state: No nitrate data found. Run aborted."
       call exit(1)
    endif

    ! No phytoplankton data is a degenerate case that someone might
    ! use to test just the physics
    if (.not. biology) then
       got_Fluores = .false.
    endif
    if(.not. got_Fluores) then
       write(stdout, *) "init_state: Phytoplankton initialized to zero. ", &
                        "Testing just the physics?"
       P%micro = 0.0d0
       endif

    ! No silicon data is okay, just initialize it to lots
    if(.not. got_Si) then
       Si = 50.0d0
    endif

    ! No dissolved inorganic carbon data? Initialize it to 2000 uM.
    if(.not. got_DIC) then
       DIC = 2.0d3
    endif

    ! No total alkalinity data? Initialize it to 2000 ueq/L.
    if(.not. got_Alk) then
       Alk = 2.0d3
    endif

    ! No dissolved oxygen data? Initialize it to 200 uM.
    if(.not. got_Oxy) then
       Oxy = 2.0d2
    endif

    ! Convert fluorescence to phytoplankton biomass expressed in uMol N
    P%micro = P%micro / N2chl
    ! Read the initial ratios of phytoplankton classes from infile,
    ! and apply them to get initial phytoplankton biomass profiles
    call getpardv("initial chl split", 3, Psplit)
    P%pico = P%micro * Psplit(3)
    P%nano = P%micro * Psplit(2)
    P%micro = P%micro * Psplit(1)
    ! Initial zooplankton biomass and detritus profiles
    Z = P%nano  ! *** hard value to estimate
    D%PON = P%micro / 5.0d0  ! estimate
    D%DON = D%PON * 10.0d0 ! estimate
    D%bSi = D%PON
    D%refr = 0.0d0
    D%DOC = D%DON * Redfield_C
    D%POC = D%PON * Redfield_C
  end subroutine init_state


  subroutine read_init_data(filename, n_records, col, data)
    ! Read records from the specified field data file to the model
    ! depth.  Return the number of records read, a column_indices data
    ! structure containing the column numbers of the various data
    ! quantities, and the array of data records themselves.
    !
    ! If the data file does not contain surface data, the data values
    ! at the surface are set to be the same as those at the shallowest
    ! depth for which data is read.  If the data file does not contain
    ! data at the model depth, the last record read will be that for
    ! the next deepest depth available.  If data doesn't reach model
    ! depth, create a record at model depth with the same data values
    ! as the deepest measurements.
    use io_unit_defs, only: stdout, field_data
    use grid_mod, only: grid
    implicit none
    ! Arguments:
    character*80, intent(in) :: filename  ! Name of data file to read
    integer, intent(out) :: n_records     ! Number of records read
    type(col_indices), intent(out):: col  ! Column number of various
                                          ! data quantities
    real(kind=dp), dimension(0:, :), intent(out) :: data  ! Data records read
    ! Local variables:
    integer :: &
         i,     &  ! Loop index
         n_cols    ! Number of data column to read per record
    logical :: data_to_model_depth

    open(field_data, file=filename)
    ! Read past header (11 lines)
    do i = 1, 11
       read(field_data, *)
    enddo
    ! Read data quantity column numbers, and set the number of column
    ! to read from each data record
    read(field_data, *) col%depth, col%T, col%S, col%Chloro, col%Fluores, &
                        col%NO, col%Phyto, col%Si, col%DIC, col%Alk, col%Oxy
    n_cols = max(col%depth, col%T, col%S, col%Chloro, col%Fluores, &
                 col%NO, col%Phyto, col%Si, col%DIC, col%Alk, col%Oxy)
    ! Read data to model depth, or next deeper record.  If data ends
    ! before model depth, read the whole file
    data_to_model_depth = .false.
    n_records = 1
    do while(.not. data_to_model_depth)
       read(field_data, *, end=100) (data(n_records, i), i = 1, n_cols)
       if(data(n_records, col%depth) >= grid%D) then
          data_to_model_depth = .true.
       endif
       if(n_records <= size(data, 1)) then
          n_records = n_records + 1
       else
          write(stdout, *) "read_init_data: "
          call exit(1)
       endif
    enddo
100 continue
    close(field_data)
    ! Surface value
    if(data(1, col%depth) == 0.0d0) then
       ! If the data contains a record at the surface, shift the array
       ! to put the surface values in record 0
       data = cshift(data, 1)
    else
       ! Otherwise, create a record at the surface that is a copy of
       ! the shallowest data record.  This will result in the profiles
       ! having constant values from the surface to the depth of the
       ! shallowest data.
       data(0, :) = data(1, :)
       data(0, col%depth) = 0.0d0
    endif
    ! If data doesn't reach model depth, create a record at model
    ! depth with the same data values as the deepest measurements.
    ! This will result in the profiles having constant values from the
    ! deepest measurement depth to the model depth.
    if(.not. data_to_model_depth) then
       data(n_records, :) = data(n_records - 1, :)
       data(n_records, col%depth) = grid%D
    else
       ! Undo loop-ending increment to get number of data records
       ! read
       n_records = n_records - 1
    endif
  end subroutine read_init_data


  function interp_to_grid(depth, data_qty) result(qty)
    ! Return an array of quantity values at the grid point depths
    ! calculated by linear interpolation from a data profile (CTD,
    ! STRATOGEM bottles, or IOS bottles data).
    use grid_mod, only: grid
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(in) :: &
         depth,   &  ! Depths from data file
         data_qty    ! Quantity values from data file
    ! Result:
    real(kind=dp), dimension(0:grid%M+1) :: &
         qty  ! Quantity values interpolated to grid point depths
    ! Local variables:
    logical, dimension(0:size(data_qty)-1) :: &
         mask  ! Mask array to filter out invalid data values
    integer :: &
         data_records, &  ! Number of valid data records
         i_g,          &  ! Index of grid point depths
         i_data           ! Index of data record
    real(kind=dp), dimension(0:size(data_qty)-1) :: &
         depth_clean, &  ! Depths at which there is valid data
         qty_clean,   &  ! Valid quantity values from data
         del_depth,   &  ! Depth differences from data
         del_qty         ! Quantity value differences from data


    ! Remove records with negative data values (typically -99.0 or
    ! -99999) because that indicates invalid data
    mask(0:size(data_qty >= 0.0d0)-1) = data_qty >= 0.0d0
    data_records = size(pack(data_qty, mask))
    qty_clean(0:data_records-1) = pack(data_qty, mask)
    depth_clean(0:data_records-1) = pack(depth, mask)
    ! Calculate depth and quantity differences from field data for use
    ! in interpolation
    del_depth(0:data_records-2) = depth_clean(1:data_records-1) &
                                  - depth_clean(0:data_records-2)
    del_qty(0:data_records-2) = qty_clean(1:data_records-1) &
                                - qty_clean(0:data_records-2)



    ! Interpolate quantity values at grid point depths
    i_data = 1
    do i_g = 0, grid%M + 1
       if(grid%d_g(i_g) > depth_clean(i_data)) then
          i_data = i_data + 1
       endif
       qty(i_g) = qty_clean(i_data-1) + del_qty(i_data-1) &
                  * ((grid%d_g(i_g) - depth_clean(i_data-1)) &
                  / del_depth(i_data-1))
    enddo
  end function interp_to_grid


  subroutine alloc_core_variables(M)
    ! Allocate memory for core variables arrays.
    use malloc, only: alloc_check
    implicit none
    ! Argument:
    integer, intent(in) :: M  ! Number of grid points
    ! Local variables:
    integer           :: allocstat  ! Allocation return status
    character(len=80) :: msg        ! Allocation failure message prefix

    msg = "Cross-estuary velocity component profile arrays"
    allocate(U%new(0:M+1), U%old(0:M+1), U%grad_i(1:M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Along-estuary velocity component profile arrays"
    allocate(V%new(0:M+1), V%old(0:M+1), V%grad_i(1:M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Temperature profile arrays"
    allocate(T%new(0:M+1), T%old(0:M+1), T%grad_i(1:M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Salinity profile arrays"
    allocate(S%new(0:M+1), S%old(0:M+1), S%grad_i(1:M), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Micro phytoplankton biomass profile array"
    allocate(P%micro(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nano phytoplankton biomass profile array"
    allocate(P%nano(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Pico phytoplankton biomass profile array"
    allocate(P%pico(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Microzooplankton concentration profile array"
    allocate(Z(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Nitrate concentration profile array"
    allocate(N%O(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Ammonium concentration profile array"
    allocate(N%H(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Silicon concentration profile array"
    allocate(Si(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved inorganic carbon concentration profile array"
    allocate(DIC(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved oxygen concentration profile array"
    allocate(Oxy(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Alkalinity profile array"
    allocate(Alk(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved organic nitrogen detritus concentration profile array"
    allocate(D%DON(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Particulate organic nitrogen detritus concentration profile array"
    allocate(D%PON(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Refractory nitrogen detritus concentration profile array"
    allocate(D%refr(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Biogenic silicon detritus concentration profile array"
    allocate(D%bSi(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Dissolved organic carbon detritus concentration profile array"
    allocate(D%DOC(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Particulate organic carbon detritus concentration profile array"
    allocate(D%POC(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
    msg = "Refractory carbon detritus concentration profile array"
    allocate(D%reC(0:M+1), stat=allocstat)
    call alloc_check(allocstat, msg)
  end subroutine alloc_core_variables


  subroutine dalloc_core_variables
    ! Deallocate memory for core variables arrays.
    use malloc, only: dalloc_check
    implicit none
    ! Local variables:
    integer           :: dallocstat  ! Deallocation return status
    character(len=80) :: msg        ! Deallocation failure message prefix

    msg = "Cross-estuary velocity component profile arrays"
    deallocate(U%new, U%old, U%grad_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Along-estuary velocity component profile arrays"
    deallocate(V%new, V%old, V%grad_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Temperature profile arrays"
    deallocate(T%new, T%old, T%grad_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Salinity profile arrays"
    deallocate(S%new, S%old, S%grad_i, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Micro phytoplankton biomass profile array"
    deallocate(P%micro, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nano phytoplankton biomass profile array"
    deallocate(P%nano, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Pico phytoplankton biomass profile array"
    deallocate(P%pico, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Microzooplankton concentration profile array"
    deallocate(Z, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Nitrate concentration profile array"
    deallocate(N%O, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Ammonium concentration profile array"
    deallocate(N%H, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Silicon concentration profile array"
    deallocate(Si, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved inorganic carbon concentration profile array"
    deallocate(DIC, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved oxygen concentration profile array"
    deallocate(Oxy, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Alkalinity profile array"
    deallocate(Alk, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved organic nitrogen detritus concentration profile array"
    deallocate(D%DON, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Particulate organic nitrogen detritus concentration profile array"
    deallocate(D%PON, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Refractory nitrogen detritus concentration profile array"
    deallocate(D%refr, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Biogenic silicon detritus concentration profile array"
    deallocate(D%bSi, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Dissolved organic carbon detritus concentration profile array"
    deallocate(D%DOC, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Particulate organic carobon detritus concentration profile array"
    deallocate(D%POC, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
    msg = "Refractory carbon detritus concentration profile array"
    deallocate(D%reC, stat=dallocstat)
    call dalloc_check(dallocstat, msg)
  end subroutine dalloc_core_variables

end module core_variables
