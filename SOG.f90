! $Id$
! $Source$

program SOG      
  ! Coupled physical and biological model of the Strait of Georgia

  ! Things that we will use from other modules
  !
  ! Type definitions:
  use datetime, only: datetime_
  !
  ! Parameter values:
  use fundamental_constants, only: g
  !
  ! Variables:
  use grid_mod, only: grid
  use core_variables, only: U, V, T, S, P, Z, N, Si, D
  use water_properties, only: rho, alpha, beta, Cp
  use physics_model, only: B
  use turbulence, only: K
  use mixing_layer, only: h
  ! *** Temporary until physics equations refactoring is completed
  use physics_eqn_builder, only: U_RHS, V_RHS, T_RHS, S_RHS
  ! *** Temporary until turbulence refactoring is completed
  use turbulence, only: u_star, L_mo, w, wbar
  !
  ! Subroutines and functions:
  use fundamental_constants, only: init_constants
  use core_variables, only: alloc_core_variables, dalloc_core_variables
  use grid_mod, only: init_grid, dalloc_grid, gradient_g, gradient_i, &
       interp_i
  use physics_model, only: init_physics, &
       new_to_old_physics, dalloc_physics_variables
  use water_properties, only: calc_rho_alpha_beta_Cp_profiles
  use physics_eqn_builder, only: build_physics_equations, &
       new_to_old_phys_RHS, new_to_old_phys_Bmatrix
  use turbulence, only: calc_KPP_diffusivity
  use baroclinic_pressure, only: new_to_old_vel_integrals, &
       baroclinic_P_gradient
  use air_sea_fluxes, only: wind_stress
  use biology_model, only: init_biology, calc_bio_rate
  use biology_eqn_builder, only: build_biology_equations, new_to_old_bio_RHS, &
       new_to_old_bio_Bmatrix
  use NPZD, only: dalloc_biology_variables
  use IMEX_solver, only: init_IMEX_solver, solve_phys_eqns, solve_bio_eqns, &
       dalloc_IMEX_variables
  use input_processor, only: init_input_processor, getpars, getpari, &
       getpard, getparl
  use timeseries_output, only: init_timeseries_output, write_std_timeseries, &
       timeseries_output_close
  use profiles_output, only: init_profiles_output, write_std_profiles, &
       profiles_output_close
  use user_output, only: write_user_timeseries, write_user_profiles
  use mixing_layer, only: find_mixing_layer_depth, &
       find_mixing_layer_indices, new_to_old_mixing_layer
  use find_upwell, only: upwell_profile, vertical_advection
  use diffusion, only: diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux
  use fitbottom, only: bot_bound_time, bot_bound_uniform
  use forcing, only: read_variation, read_forcing, get_forcing
  use precision_defs, only: dp, sp
  use io_unit_defs, only: stripped_infile, stdout
  use unit_conversions, only: KtoC
  use datetime, only: os_datetime, calendar_date, clock_time, datetime_str

  ! Inherited modules
  ! *** Goal is to make these go away
  use declarations
  use surface_forcing, only: precision, step_guess, step_min
  use initial_sog, only: initial_mean

  implicit none

  real(kind=dp) :: &
       upwell
  ! Upwelling constant (tuned parameter)
  !*** read in here, used by surface_flux_sog : eventually should be local
  ! to the surface_forcing module (not be be confused with current 
  ! surface_forcing module
  real(kind=dp) :: upwell_const, S_riv, sumS=0, sumSriv=0
  integer :: scount=0

  ! Interpolated river flows
  real(kind=sp) :: Qinter  ! Fraser River
  real(kind=sp) :: Einter  ! Englishman River
  ! Current time met data
  real(kind=sp) :: cf_value, atemp_value, humid_value
  ! Wind data
  real(kind=dp) unow, vnow

  integer :: wind_n, met_n, river_n ! length of various forcing files

  ! Initial month parameter read from run control input file
  ! and passed to new_year()
  integer :: month_o

  ! Iteration limit for inner loop that solves KPP turbulence model
  integer :: niter

  ! Scales required for convergence
  real(kind=dp) uscale, hscale

  ! Value of h%new, u%new and v%new from previous iteration for 
  ! blending with new value to stabilize convergence of the implicit 
  ! solver loop
  real(kind=dp) :: hprev
  real(kind=dp), dimension (:), allocatable :: uprev, vprev

  ! Variables for velocity component convergence metric
  double precision :: &
       & delh, delu, delv

  ! Code identification string (maintained by CVS), for output file headers
  character*70 :: &
       codeId = "$Source$"

  ! Date/time structures for output file headers
  type(datetime_) :: runDatetime     ! Date/time of code run
  type(datetime_) :: startDatetime   ! Date/time of initial conditions


  ! ---------- Beginning of Initialization Section ----------
  ! Get the current date/time from operating system to timestamp the
  ! run with
  call os_datetime(runDatetime)

  ! Initialize the input processor, including preparing the infile
  ! (read from stdin) to be read by the init_* subroutines called
  ! below.
  call init_input_processor(datetime_str(runDatetime))

  ! Initialize fundamental constants
  ! *** This is here because of the pgf90 bug that prevents f from
  ! *** being a parameter
  call init_constants()

  ! Read the run start date/time, duration, and time step
  ! *** Not sure where these should go?
  year_o = getpari("year_o")
  month_o = getpari("month_o")
  day_o = getpari("yr_day_o")
  t_o = getpari("t_o")
  t_f = getpard("run_dur")
  dt = getpard("dt")
  ! Calculate the number of time steps for the run (note that int()
  ! rounds down)
  steps = 1 + int((t_f - t_o) / dt)
  ! Calculate the month number and month day for output file headers
  startDatetime%yr = year_o
  startDatetime%yr_day = day_o
  startDatetime%day_sec = t_o
  call calendar_date(startDatetime)
  call clock_time(startDatetime)

  ! Initialize time series writing code
  call init_timeseries_output(codeId, datetime_str(runDatetime), &
       startDatetime)
  ! Initialize profiles writing code
  call init_profiles_output(codeId, datetime_str(runDatetime), &
       startDatetime)

  ! Initialize the grid
  call init_grid

  ! Allocate memory for the profiles arrays of the core variables that
  ! we are calculating, i.e. U, V, T, S, etc.
  call alloc_core_variables(grid%M)

  ! Initialize the physics model
  call init_physics(grid%M)

  ! Initialize the biology model
  call init_biology(grid%M)

  ! Initialize the IMEX semi-implicit PDE solver
  ! *** This should eventually go into init_numerics()
  call init_IMEX_solver(grid%M)

  ! Allocate memory
  ! *** It would be nice if everything in allocate[13] could end up in 
  ! *** alloc_* subroutines private to various modules, and called by 
  ! *** their init_* subroutines.
  CALL allocate1(grid%M, alloc_stat) 
  DO xx = 1,12
     IF (alloc_stat(xx) /= 0) THEN
        PRINT "(A)","ALLOCATION failed.  KPP.f  xx:"
        PRINT *,xx,alloc_stat(xx)
        STOP
     END IF
  END DO

  ! *** great refractory guru to fix this out of here
  allocate (uprev(0:grid%M+1), vprev(0:grid%M+1), &
       STAT = alloc_stat(1))
  if (alloc_stat(1) /= 0) then
     print "(a)", " Allocation failed for uprev, vprev"
     stop
  endif

  ! Length of forcing data files (move to be read in)
  wind_n = 46056 - 8 ! with wind shifted to LST we lose 8 records
  met_n = 1918
  river_n = 1826
  call read_forcing (wind_n, met_n, river_n)
  call read_variation

  ! Read the physic model parameter values
  ! *** These should go into the freshwater module.
  upwell_const = getpard("upwell_const")
  Fw_scale = getpard('Fw_scale')     ! Fresh water scale factor for river flows
  Fw_surface = getparl('Fw_surface') ! Add all fresh water on surface?
  if (.not. Fw_surface) then
     Fw_depth = getpard('Fw_depth')  ! Depth to distribute freshwater flux over
  endif

  CALL initialize ! initializes everything (biology too)

  ! Read the cruise id from stdin to use to build the file name for
  ! nutrient initial conditions data file
  cruise_id = getpars("cruise_id")
  CALL initial_mean(U%new, V%new, T%new, S%new, P%micro, P%nano, &
       Z, N%O, N%H, Si, D%DON, D%PON, D%refr, D%bSi, &
       h%new, grid, cruise_id)

  ! Read the iteration count limit for the physics model implicit
  ! solver
  ! *** This should go somewhere else - maybe init_physics() (yet to come)
  niter = getpari("niter")

  ! Initialize the profiles of the water column properties
  ! Density (rho), thermal expansion (alpha) and saline contraction (beta)
  ! coefficients, and specific heat capacity (Cp).
  ! This calculates the values of the grid point depth arrays (*%g) of
  ! the 4 water properties.
  call calc_rho_alpha_beta_Cp_profiles(T%new, S%new)
  ! Interpolate the values of density and specific heat capacity at
  ! the grid layer interface depths
  rho%i = interp_i(rho%g)
  alpha%i = interp_i(alpha%g)
  beta%i = interp_i(beta%g)
  Cp%i = interp_i(Cp%g)
  ! Calculate the gradients of density, and thermal expansion and saline
  ! contraction coefficients at the grid layer interface depths
  rho%grad_i = gradient_i(rho%g)
  alpha%grad_i = gradient_i(alpha%g)
  beta%grad_i = gradient_i(beta%g)

  ! Close the input parameters file
  close(stripped_infile)
  ! ---------- End of Initialization Section ----------

  do time_step = 1, steps  !---------- Beginning of the time loop ----------
     ! Store %new components of various variables from time step just
     ! completed in %old their components for use in the next time
     ! step

     h%old = h%new

     call new_to_old_physics()
     call new_to_old_vel_integrals()
     call new_to_old_mixing_layer()
     call new_to_old_phys_RHS()
     call new_to_old_phys_Bmatrix()
     call new_to_old_bio_RHS()
     call new_to_old_bio_Bmatrix()

     ! Get forcing data
     call get_forcing (year, day, day_time, &
          Qinter, Einter, cf_value, atemp_value, humid_value, &
          unow, vnow)

     CALL irradiance_sog(cf_value, day_time, day, &
          I, I_par, grid, jmax_i, Q_sol, euph, Qinter, P%micro, P%nano)

     DO count = 1, niter !------ Beginning of the implicit solver loop ------
        ! Calculate gradient pofiles of the velocity field and water column
        ! temperature, and salinity at the grid layer interface depths
        U%grad_i = gradient_i(U%new)
        V%grad_i = gradient_i(V%new)
        T%grad_i = gradient_i(T%new)
        S%grad_i = gradient_i(S%new)

        ! Calculate surface forcing components
        ! *** Confirm that all of these arguments are necessary
        CALL surface_flux_sog(grid%M, rho%g, wt_r, S%new(1),        &
             S%old(1), S_riv, T%new(0), j_gamma, I, Q_t,        &
             alpha%i(0), Cp%i(0), beta%i(0), unow, vnow, cf_value/10.,    &
             atemp_value, humid_value, Qinter, &
             day, dt/h%new, &
             upwell_const, upwell, Einter,       &
             u%new(1), dt, Fw_surface, Fw_scale, Ft, count) 

        ! Calculate the wind stress, and the turbulent kinenatic flux
        ! at the surface.
        !
        ! This calculates the values of the wbar%u(0), and wbar%v(0)
        ! variables that are declared in the turbulence module.
        call wind_stress (unow, vnow, rho%g(1))

        ! Calculate nonturbulent heat flux profile
        Q_n = I / (Cp%i * rho%i)

        ! Calculate the nonturbulent fresh water flux profile, and its
        ! contribution to the salinity profile
        ! *** Move to a subroutine
        if (Fw_surface) then
           F_n = 0.0d0
        else
           Fw = Ft * exp(-grid%d_i / (Fw_depth * h%new))
           F_n = 29.626 * Fw
        endif

        ! Calculate buoyancy profile, and surface buoyancy forcing
        CALL buoyancy(grid, T%new, S%new, h%new, I, F_n,   &  ! in
             rho%g, alpha%g, beta%g, Cp%g, Fw_surface,     &  ! in
             B, Bf)                                           ! out

        ! Calculate the turbulent diffusivity profile arrays using the
        ! K Profile Parameterization (KPP) of Large, et al (1994)
        !
        ! This calculates the values of the nu%*, K_ML%*, and K%*
        ! variables that are declared in the turbulence module so that
        ! they can be used by other modules.
        call calc_KPP_diffusivity(Bf, h%new, h%i, h%g)
!!$! *** Turbulence refactoring bridge code
omega%m%value = w%m%profile
omega%s%value = w%s%profile

! *** Mixing layer depth refactoring bridge code
oh%new = h%new
oh%g = h%g
oh%i = h%i
        
        ! Calculates non-local transport term (Large, et al (1994),
        ! eqn 20).  Only applicable to scalar quantities under
        ! convective forcing.
        IF (abs(u_star) > epsilon(u_star) .or. abs(Bf) > epsilon(Bf)) THEN
           CALL def_gamma(L_mo, grid, wt_r, oh, gamma, Bf, omega) 
        ELSE
           gamma%t = 0.0d0
           gamma%s = 0.0d0
           gamma%m = 0.0d0
        END IF

        ! Build the terms of the semi-implicit diffusion/advection
        ! PDEs with Coriolis and baroclinic pressure gradient terms for
        ! the physics quantities.
        !
        ! This calculates the values of the precursor diffusion
        ! coefficients matrices (Bmatrix%vel%*, Bmatrix%T%*, and
        ! Bmatrix%S%*), the RHS diffusion/advection term vectors
        ! (*_RHS%diff_adv%new), and the RHS Coriolis and barolcinic
        ! pressure gradient term vectors (*_RHS%C_pg).
        call build_physics_equations(grid, dt, U%new, V%new, T%new, S%new) ! in
     
! *** Physics equations refactoring bridge code
Gvector%u = U_RHS%diff_adv%new
Gvector%v = V_RHS%diff_adv%new
Gvector%t = T_RHS%diff_adv%new
Gvector%s = S_RHS%diff_adv%new

        ! Calculate profile of upwelling velocity
        call upwell_profile(grid, Qinter, upwell, wupwell)
        ! Upwell salinity, temperature, and u & v velocity components
        ! similarly to nitrates
        call vertical_advection (grid, dt, S%new, wupwell, &
             Gvector%s)
        call vertical_advection (grid, dt, T%new, wupwell, &
             Gvector%t)
        call vertical_advection (grid, dt, U%new, wupwell, &
             Gvector%u)
        call vertical_advection (grid, dt, V%new, wupwell, &
             Gvector%v)
! *** Physics equations refactoring bridge code
U_RHS%diff_adv%new = Gvector%u
V_RHS%diff_adv%new = Gvector%v
T_RHS%diff_adv%new = Gvector%t
S_RHS%diff_adv%new = Gvector%s

        ! Store %new components of RHS and Bmatrix variables in %old
        ! their components for use by the IMEX solver.  Necessary for the
        ! 1st time step because the values just calculated are a better
        ! estimate than zero.
        if (time_step == 1) then
           call new_to_old_phys_RHS()
           call new_to_old_phys_Bmatrix()
        endif

        uprev = u%new
        vprev = v%new

        ! Solve the semi-implicit diffusion/advection PDEs with
        ! Coriolis and baroclinic pressure gradient terms for the
        ! physics quantities.
        call solve_phys_eqns(grid%M, day, time, &  ! in
             U%old, V%old, T%old, S%old,        &  ! in
             U%new, V%new, T%new, S%new)           ! out

        ! Update boundary conditions at surface, and bottom of grid
        U%new(0) = U%new(1)
        V%new(0) = V%new(1)
        S%new(0) = S%new(1)   
        T%new(0) = T%new(1)
        U%new(grid%M+1) = U%new(grid%M) 
        V%new(grid%M+1) = V%new(grid%M) 
        T%new(grid%M+1) = T%old(grid%M+1)
        S%new(grid%M+1) = S%old(grid%M+1)

        ! Update gradient pofiles of the water column
        ! temperature, and salinity at the grid layer interface depths
        T%grad_i = gradient_i(T%new)
        S%grad_i = gradient_i(S%new)

        ! Update the profiles of the water column properties
        ! Density (rho), thermal expansion (alpha) and saline
        ! contraction (beta) coefficients, and specific heat capacity (Cp)
        ! This calculates the values of the grid point depth arrays (*%g) of
        ! the 4 water properties.
        call calc_rho_alpha_beta_Cp_profiles(T%new, S%new)
        ! Interpolate the values of density and specific heat capacity at
        ! the grid layer interface depths
        rho%i = interp_i(rho%g)
        alpha%i = interp_i(alpha%g)
        beta%i = interp_i(beta%g)
        Cp%i = interp_i(Cp%g)
        ! Calculate the gradients of thermal expansion and saline
        ! contraction coefficients at the grid layer interface depths
        rho%grad_i = gradient_i(rho%g)
        alpha%grad_i = gradient_i(alpha%g)
        beta%grad_i = gradient_i(beta%g)

        ! Update buoyancy profile, surface buoyancy forcing, and local
        ! buoyancy frequency
        ! *** Local buoyancy freqency may be part of mixing layer
        ! *** depth calculation
        CALL buoyancy(grid, T%new, S%new, h%new, I, F_n,   &  ! in
             rho%g, alpha%g, beta%g, Cp%g, Fw_surface,     &  ! in
             B, Bf)                                           ! out
        rho%grad_g = gradient_g(rho%g)
        DO xx = 1,grid%M
           N_2_g(xx) = -(g / rho%g(xx)) * rho%grad_g(xx)
        END DO

!!!Define turbulent velocity shear V_t_square for Ri_b  use last count values

!!!Calculate beta_t: need h%e and h%e%i and new w%b  w%i vertical turbulent flux of i


        wbar%b_err(0) = 0.

!!$        DO xx = 1, grid%M           !uses K%old and T%new
!!$           w%t(xx) = -K%t(xx)*(T%grad_i(xx) - gamma%t(xx))
!!$           w%s(xx) = -K%s(xx)*(S%grad_i(xx) - gamma%s(xx))
!!$           w%b(xx) = g * (alpha%i(xx) * w%t(xx) - beta%i(xx) * w%s(xx))   
!!$           w%b_err(xx) = g * (alpha%grad_i(xx) * w%t(xx) &
!!$                - beta%grad_i(xx) * w%s(xx))
!!$        END DO


        ! beta_t is the ratio of the entrainment flux to the surface buoyancy flux -- set equal to -0.2 when convection is occuring

        if (wbar%b(0)>0) then
           beta_t = -0.2
        else
           beta_t = 0.
        endif

        V_t_square = 0.


        IF (beta_t < 0.) THEN  !omega, vertical velocity scale
           CALL def_v_t_sog(grid, oh, N_2_g, omega%s%value, V_t_square, beta_t, L_mo) !test conv  
           ! V_t_square, the turbulent velocity shear, (23)
        END IF

        ! Calculate the profile of bulk Richardson number (Large, etal
        ! (1994) eqn 21)
        CALL define_Ri_b_sog(grid, oh, surface_h, U%new, V%new, rho%g, Ri_b, &
             V_t_square, N_2_g)

        hprev = h%new

        ! Find the mixing layer depth by comparing Ri_b to Ri_c
        ! This sets the values of the components of h%*.
        call find_mixing_layer_depth (grid, Ri_b, Bf, year, day, day_time, &
             count)
        ! Average the newly calculated mixing layer depth with that
        ! from the previous iteration to help stabilize the
        ! convergence of the implicit solver loop
        ! (and if things aren't converging, try something different)
        if (count.ge.10.and.count.le.13) then
           h%new = 0.75*h%new + 0.25*hprev
           u%new = 0.75*u%new + 0.25*uprev
           v%new = 0.75*v%new + 0.25*vprev
        elseif (count.ge.23.and.count.le.26) then
           h%new = 0.25*h%new + 0.75*hprev
           u%new = 0.25*u%new + 0.75*uprev
           v%new = 0.25*v%new + 0.75*vprev
        else
           h%new = (hprev + h%new) / 2.
           u%new = (uprev + u%new) / 2.
           v%new = (vprev + v%new) / 2.
        endif

        ! Set the value of the indices of the grid point & grid layer
        ! interface immediately below the mixing layer depth.
        call find_mixing_layer_indices()

        ! Calculate the mixing layer depth convergence metric for
        ! measuring the progress of the implicit solver
        ! Accuracy is 2% at large mixed layer 
        ! depths and about 4 cm at 1 m mixed layer depths
        hscale = 0.02 * hprev + 0.02
        delh = ABS(hprev - h%new)/hscale


        ! Calculate convergence metric for velocity field
        ! Accuracy is 2% at large 
        ! velocities and 0.7 cm/s at 10.0 cm/s (units are m/s)
        uscale = 0.02*sqrt(uprev(1)**2+vprev(1)**2) + 0.005
        delu = abs(u%new(1)-uprev(1))/uscale
        delv = abs(v%new(1)-vprev(1))/uscale

        ! Calculate baroclinic pressure gradient components
        !
        ! This calculates the values of the dPdx_b, and dPdy_b arrays.
        call baroclinic_P_gradient(grid, dt, U%new, V%new, rho%g)


        ! Test for convergence of implicit solver
        if (count >= 2 .and. max(delh,delu,delv) < 1) then
           exit
        endif
     enddo  !---------- End of the implicit solver loop ----------

     !---------- Biology Model ----------
     !
     ! Solve the biology model ODEs to advance the biology quantity values
     ! to the next time step, and calculate the growth - mortality terms
     ! (*_RHS%bio) of the semi-implicit diffusion/advection equations.
     call calc_bio_rate(time, day, dt, grid%M, precision, step_guess, step_min,  &
          T%new(0:grid%M), I_Par, P%micro, P%nano, Z, N%O, N%H, Si,              &
          D%DON, D%PON, D%refr, D%bSi)

     ! Build the rest of the terms of the semi-implicit diffusion/advection
     ! PDEs for the biology quantities.
     !
     ! This calculates the values of the precursor diffusion
     ! coefficients matrix (Bmatrix%bio%*), the RHS diffusion/advection
     ! term vectors (*_RHS%diff_adv%new), and the RHS sinking term
     ! vectors (*_RHS%sink).
     call build_biology_equations(grid, dt, P%micro, P%nano, Z, N%O, N%H, &! in
          Si, D%DON, D%PON, D%refr, D%bSi, Ft, wupwell)                    ! in

     ! Store %new components of RHS and Bmatrix variables in %old
     ! their components for use by the IMEX solver.  Necessary for the
     ! 1st time step because the values just calculated are a better
     ! estimate than zero.
     if (time_step == 1) then
        call new_to_old_bio_RHS()
        call new_to_old_bio_Bmatrix()
     endif

     ! Solve the semi-implicit diffusion/advection PDEs for the
     ! biology quantities.
     call solve_bio_eqns(grid%M, P%micro, P%nano, Z, N%O, N%H, Si, &
          D%DON, D%PON, D%refr, D%bSi, day, time)
     !
     !---------- End of Biology Model ----------

     ! Update boundary conditions at surface
     P%micro(0) = P%micro(1)
     P%nano(0) = P%nano(1)
     Z(0) = Z(1)
     N%O(0) = N%O(1)
     N%H(0) = N%H(1)
     Si(0) = Si(1)
     D%DON(0) = D%DON(1)
     D%PON(0) = D%PON(1)
     D%refr(0) = D%refr(1)
     D%bSi(0) = D%bSi(1)

     !--------bottom boundaries--------------------------
     ! Update boundary conditions at bottom of grid
     !
     ! For those variables that we have data, use the annual fit
     ! calculated from the data
     call bot_bound_time(day, day_time, &                             ! in
          T%new(grid%M+1), S%new(grid%M+1), N%O(grid%M+1), &          ! out
          Si(grid%M+1), P%micro(grid%M+1), P%nano(grid%M+1)) ! out
     ! For those variables that we have no data for, assume uniform at
     ! bottom of domain
     call bot_bound_uniform(grid%M, Z, N%H, D%DON, D%PON, D%refr, D%bSi)

     ! Write standard time series results
     ! !!! Please don't change this argument list without good reason. !!!
     ! !!! If it is changed, the change should be committed to CVS.    !!!
     ! !!! For exploratory, debugging, etc. output use                 !!!
     ! !!! write_user_timeseries() below.                              !!!
     call write_std_timeseries(time / 3600., grid,                    &
       ! Variables for standard physics model output
       count, h%new, U%new, V%new, T%new, S%new,                      &
       ! Variables for standard biology model output
       N%O , N%H, Si, P%micro, P%nano, Z, D%DON, D%PON, D%refr, D%bSi)

     ! Write user-specified time series results
     ! !!! Please don't add arguments to this call.           !!!
     ! !!! Instead put use statements in your local copy of   !!!
     ! !!! write_user_timeseries() in the user_output module. !!!
     call write_user_timeseries(time / 3600., grid)
     
     ! Write standard profiles results
     ! !!! Please don't change this argument list without good reason. !!!
     ! !!! If it is changed, the change should be committed to CVS.    !!!
     ! !!! For exploratory, debugging, etc. output use                 !!!
     ! !!! write_user_timeseries() below.                              !!!
     call write_std_profiles(codeId, datetime_str(runDatetime),       &
          datetime_str(startDatetime), year, day, day_time, dt, grid, &
          T%new, S%new, rho%g, P%micro, P%nano, Z, N%O, N%H, Si,      &
          D%DON, D%PON, D%refr, D%bSi, K%m, K%T, K%S,                 &
          I_par, U%new, V%new)

     ! Write user-specified profiles results
     ! !!! Please don't add arguments to this call.           !!!
     ! !!! Instead put use statements in your local copy of   !!!
     ! !!! write_user_profiles() in the user_output module.   !!!
     call write_user_profiles(codeId, datetime_str(runDatetime),      &
          datetime_str(startDatetime), year, day, day_time, dt, grid)

     ! Increment time, calendar date and clock time
     call new_year(day_time, day, year, time, dt, month_o)
     scount = scount + 1
     !*** Fix this to be grid independent
     sumS = sumS + 0.5*(S%new(2)+S%new(3))
     sumSriv = sumSriv + S_riv
  end do  !--------- End of time loop ----------

  write (stdout,*) "For Ft tuning"
  write (stdout,*) "Average SSS should be", sumSriv/float(scount)
  write (stdout,*) "Average SSS was", sumS/float(scount)

  ! Close output files
  call timeseries_output_close
  call profiles_output_close

  ! Deallocate memory
  call dalloc_grid
  call dalloc_core_variables
  call dalloc_physics_variables
  call dalloc_biology_variables
  call dalloc_IMEX_variables

end program SOG
