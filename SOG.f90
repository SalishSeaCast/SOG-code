! $Id$
! $Source$

program SOG      
  ! Coupled physical and biological model of the Strait of Georgia

  ! Things that we will use from other modules
  !
  ! Type definitions:
  use precision_defs, only: dp, sp
  use datetime, only: datetime_
  !
  ! Parameter values:
  use io_unit_defs, only: stripped_infile, stdout
  use fundamental_constants, only: g
  !
  ! Variables:
  use grid_mod, only: &
       grid, &  ! Grid parameters and depth and spacing arrays
       interp_value  ! interpolate salinity grid to find 1 m value
  use core_variables, only: &
       U,  &  ! Cross-strait (35 deg) velocity component arrays
       V,  &  ! Along-strait (305 deg) velocity component arrays
       T,  &  ! Temperature profile arrays
       S,  &  ! Salinity profile arrays
       P,  &  ! Micro & nano & pico phytoplankton profile arrays
       Z,  &  ! Microzooplankton profile array
       N,  &  ! Nitrate & ammonium concentration profile arrays
       Si, &  ! Silicon concentration profile arrays
       D      ! Detritus concentration profile arrays
  use water_properties, only: &
       rho,   &  ! Density [kg/m^3]
       alpha, &  ! Thermal expansion coefficient [K^-1]
       beta,  &  ! Saline contraction coefficient [K^-1]
       Cp        ! Specific heat capacity [J/kg.K]
  use buoyancy, only: &
       Bf  ! Surface buoyancy forcing
  use turbulence, only: &
       K  ! Overall diffusivity profile; a continuous profile of
          ! K_ML%* in the mixing layer, and nu%*%total below it
  use mixing_layer, only: &
       h  ! Mixing layer depth value & indices
  use numerics, only: &
       initDatetime, &   ! Date/time of initial conditions
       year, &   ! Year counter
       day,  &   ! Year-day counter
       day_time, &  ! Day-sec counter
       time, &  ! Time counter through run; seconds since midnight of 
                ! initDatetime.
       dt, &      ! Time step [s]
       steps, &   ! Number of time steps in the main time loop
       max_iter, &  ! Maximum number of iterations allowed for
                    ! implicit solver loop.
       hprev, &  ! Values of h%new, u%new and v%new from previous
       Uprev, &  ! iteration for blending with new value to stabilize
       Vprev, &  ! convergence of the implicit solver loop
       new_weight, prev_weight, &  ! Weighting factors for convergence
                                   ! stabilization blending.
       vel_scale, h_scale, &  ! Scale factors used to normalize the
                              ! velocity components, and mixing layer
                              ! depth values to calculate the
                              ! convergence metrics for the implicit
                              ! solver loop.
       del  ! Convergence metrics for implicit solver loop.

  ! *** Temporary until physics equations refactoring is completed
  use physics_eqn_builder, only: U_RHS, V_RHS, T_RHS, S_RHS

  !
  ! Subroutines and functions:
  use fundamental_constants, only: init_constants
  use input_processor, only: init_input_processor, getpars, &
       getpard, getparl
  use grid_mod, only: init_grid, dalloc_grid, gradient_g, gradient_i, &
       interp_i
  use numerics, only: init_numerics, dalloc_numerics_variables
  use core_variables, only: init_core_variables, dalloc_core_variables
  use physics_model, only: init_physics, &
       new_to_old_physics, dalloc_physics_variables
  use water_properties, only: calc_rho_alpha_beta_Cp_profiles
  use air_sea_fluxes, only: wind_stress
  use freshwater, only: freshwater_phys, S_riv
  use buoyancy, only: calc_buoyancy
  use baroclinic_pressure, only: new_to_old_vel_integrals, &
       baroclinic_P_gradient, w_wind
  use turbulence, only: calc_KPP_diffusivity
  use physics_eqn_builder, only: build_physics_equations, &
       new_to_old_phys_RHS, new_to_old_phys_Bmatrix
  use biology_model, only: init_biology, calc_bio_rate
  use biology_eqn_builder, only: build_biology_equations, new_to_old_bio_RHS, &
       new_to_old_bio_Bmatrix
  use NPZD, only: dalloc_biology_variables
  use IMEX_solver, only: init_IMEX_solver, solve_phys_eqns, solve_bio_eqns, &
       dalloc_IMEX_variables
  use timeseries_output, only: init_timeseries_output, write_std_timeseries, &
       timeseries_output_close
  use profiles_output, only: init_profiles_output, write_std_profiles, &
       profiles_output_close
  use user_output, only: write_user_timeseries, write_user_profiles
  use mixing_layer, only: find_mixing_layer_depth, &
       find_mixing_layer_indices
  use find_upwell, only: upwell_profile, vertical_advection
  use fitbottom, only: init_fitbottom, bot_bound_time, bot_bound_uniform
  use forcing, only: read_variation, read_forcing, get_forcing
  use unit_conversions, only: KtoC
  use datetime, only: os_datetime, calendar_date, clock_time, datetime_str
  use increment_time, only: new_year
  use irradiance, only: irradiance_sog

  ! Inherited modules
  ! *** Goal is to make these go away
  use declarations
  use surface_forcing, only: precision, step_guess, step_min

  implicit none

  ! Local variables:
  !
  ! Date/time structures for output file headers
  type(datetime_) :: runDatetime     ! Date/time of code run


  real(kind=dp) :: &
       sumS=0, sumSriv=0, Sone
  integer :: scount=0, junk

  ! Interpolated river flows
  real(kind=dp) :: Qinter  ! Fraser River
  real(kind=dp) :: Einter  ! Englishman River
  ! River Temperature (of major river)
  real(kind=dp) :: RiverTemp
  ! Current time met data
  real(kind=sp) :: cf_value, atemp_value, humid_value
  ! Wind data
  real(kind=dp) unow, vnow



  ! ---------- Beginning of Initialization Section ----------
  ! Get the current date/time from operating system to timestamp the
  ! run with
  call os_datetime(runDatetime)

  ! Initialize the input processor, including preparing the infile
  ! (read from stdin) to be read by the init_* subroutines called
  ! below.
  call init_input_processor(datetime_str(runDatetime))

  ! Initialize the salinity linearity file
  open (unit=129, file="salinity_check")

  ! Initialize fundamental constants
  call init_constants()

  ! Initialize the grid
  call init_grid()

  ! Initialize the numerics
  call init_numerics(grid%M)

  ! Read initial conditions and forcing variation parameters from
  ! infile
  ! *** TODO: Decouple initial conditions variation from core
  ! *** variables initialization
  call read_variation

  ! Initialize the profile arrays of the core variables that we are
  ! calculating, i.e. U, V, T, S, etc.
  call init_core_variables()

  ! Initialize time series writing code
  call init_timeseries_output(datetime_str(runDatetime), initDatetime)

  ! Initialize profiles writing code
  call init_profiles_output(datetime_str(runDatetime), initDatetime)

  ! Initialize fitbottom
  call init_fitbottom()

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
        CALL EXIT(1)
     END IF
  END DO

  ! Read forcing data files
  call read_forcing

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

  ! ---------- End of Initialization Section ---------

  do time_step = 1, steps  !---------- Beginning of the time loop ----------
     ! Store %new components of various variables from time step just
     ! completed in %old their components for use in the next time
     ! step.
     call new_to_old_physics()
     call new_to_old_vel_integrals()
     call new_to_old_phys_RHS()
     call new_to_old_phys_Bmatrix()
     call new_to_old_bio_RHS()
     call new_to_old_bio_Bmatrix()

    

     ! Get forcing data
     call get_forcing(year, day, day_time, &
          Qinter, Einter, RiverTemp, cf_value, atemp_value, humid_value, &
          unow, vnow)

     call irradiance_sog(cf_value, day_time, day, I, I_par, grid, &
          Qinter, P%micro, P%nano, P%pico, jmax_i)

     DO count = 1, max_iter !------ Beginning of the implicit solver loop ------
        ! Calculate gradient pofiles of the velocity field and water column
        ! temperature, and salinity at the grid layer interface depths
        U%grad_i = gradient_i(U%new)
        V%grad_i = gradient_i(V%new)
        T%grad_i = gradient_i(T%new)
        S%grad_i = gradient_i(S%new)

        ! Calculate surface forcing components
        ! *** Confirm that all of these arguments are necessary
        CALL surface_flux_sog(grid%M, rho%g, &
             T%new(0), I, Q_t,        &
             alpha%i(0), Cp%i(0), beta%i(0), unow, vnow, cf_value/10.,    &
             atemp_value, humid_value)

        ! Calculate the wind stress, and the turbulent kinenatic flux
        ! at the surface.
        !
        ! This calculates the values of the wbar%u(0), and wbar%v(0)
        ! variables that are declared in the turbulence module.
        call wind_stress (unow, vnow, rho%g(1))

        ! Calculate nonturbulent heat flux profile
        Q_n = I / (Cp%i * rho%i)

        ! Calculate the nonturbulent fresh water flux profile, and its
        ! contribution to the salinity profile.
        !
        ! This sets the values of the river salinity diagnostic
        ! (S_riv), the surface turbulent kinematic salinity flux
        ! (wbar%s(0)), and the profile of fresh water contribution to
        ! the salinity (F_n)
        call freshwater_phys(Qinter, Einter, RiverTemp, S%old(1), T%old(1), T%old(grid%M+1), h%new)

        ! Calculate the buoyancy profile, surface turbulent kinematic
        ! buoyancy flux, and the surface buoyancy forcing.
        !
        ! This sets the value of the diagnostic buoyancy profile (B),
        ! the surface turbulent kinematic buoyancy flux (wbar%b(0)),
        ! and the surface buoyancy flux (Bf).
        call calc_buoyancy(T%new, S%new, h%new, I, rho%g, alpha%g, &
             beta%g, Cp%g)

        ! Calculate the turbulent diffusivity profile arrays using the
        ! K Profile Parameterization (KPP) of Large, et al (1994)
        !
        ! This calculates the values of the nu%*, K_ML%*, and K%*
        ! variables that are declared in the turbulence module so that
        ! they can be used by other modules.
        call calc_KPP_diffusivity(Bf, h%new, h%i, h%g)

        ! Calculate baroclinic pressure gradient components
        !
        ! This calculates the values of the dPdx_b, and dPdy_b arrays.
        call baroclinic_P_gradient(dt, U%new, V%new, rho%g)

        ! Build the terms of the semi-implicit diffusion/advection
        ! PDEs with Coriolis and baroclinic pressure gradient terms for
        ! the physics quantities.
        !
        ! This calculates the values of the precursor diffusion
        ! coefficients matrices (Bmatrix%vel%*, Bmatrix%T%*, and
        ! Bmatrix%S%*), the RHS diffusion/advection term vectors
        ! (*_RHS%diff_adv%new), and the RHS Coriolis and barolcinic
        ! pressure gradient term vectors (*_RHS%C_pg).
        call build_physics_equations(dt, U%new, V%new, T%new, S%new)

! *** Physics equations refactoring bridge code
Gvector%u = U_RHS%diff_adv%new
Gvector%v = V_RHS%diff_adv%new
Gvector%t = T_RHS%diff_adv%new
Gvector%s = S_RHS%diff_adv%new

        ! Calculate profile of upwelling velocity
        call upwell_profile(Qinter, wupwell)
       
        ! Upwell salinity, temperature, and u & v velocity components
        ! similarly to nitrates
        call vertical_advection (grid, dt, S%new, wupwell+ w_wind, &
             Gvector%s)
        call vertical_advection (grid, dt, T%new, wupwell+ w_wind, &
             Gvector%t)
        call vertical_advection (grid, dt, U%new, wupwell+ w_wind, &
             Gvector%u)
        call vertical_advection (grid, dt, V%new, wupwell+ w_wind, &
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

        ! Preserve the profiles of the velocity components from the
        ! previous iteration to use in averaging below for convergence
        ! stabilization.
        Uprev = U%new
        Vprev = V%new

        ! Solve the semi-implicit diffusion/advection PDEs with
        ! Coriolis and baroclinic pressure gradient terms for the
        ! physics quantities.
        call solve_phys_eqns(grid%M, day, time, &  ! in
             U%old, V%old, T%old, S%old,        &  ! in
             U%new, V%new, T%new, S%new)           ! out
        if (minval(S%new) < 0) then
           write (*,*) 'Salinity less than 0'
           stop
        endif
    
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
        ! Calculate the gradients of density, thermal expansion, and
        ! saline contraction coefficients at the grid layer interface
        ! depths.
        rho%grad_i = gradient_i(rho%g)
        alpha%grad_i = gradient_i(alpha%g)
        beta%grad_i = gradient_i(beta%g)
        ! Calculate the gradient of denisty at the grid point depths.
        rho%grad_g = gradient_g(rho%g)

        ! Update the buoyancy profile, surface turbulent kinematic
        ! buoyancy flux, and the surface buoyancy forcing.
        !
        ! This sets the value of the diagnostic buoyancy profile (B),
        ! the surface turbulent kinematic buoyancy flux (wbar%b(0)),
        ! and the surface buoyancy flux (Bf).  .
        call calc_buoyancy(T%new, S%new, h%new, I, rho%g, alpha%g, &
             beta%g, Cp%g)

        ! Preserve the value of the mixing layer depth from the
        ! previous iteration to use in averaging below for convergence
        ! stabilization.
        hprev = h%new
        ! Find the mixing layer depth by comparing the profile of bulk
        ! Richardson number to a critical value.
        ! 
        ! This sets the values of the components of h%*.
        call find_mixing_layer_depth(year, day, day_time, count)

        ! Blend the newly calculated mixing layer depth, and velocity
        ! component profiles with those from the previous iteration to
        ! help stabilize the convergence of the implicit solver loop.
        ! Try "nudging things" if convergence isn't progressing well.
        !
        ! *** These count ranges for "nudging things" are base on a
        ! *** maximum iteration limit of max_iter = 30.  They should
        ! *** be adjusted if that infile parameter value is changed
        ! *** (and Susan has very definite ideas on *how* they should
        ! *** be changed).
        if (count >= 10 .and. count <= 13) then
           new_weight = 0.75d0
           prev_weight = 0.25d0
        elseif (count >= 23 .and. count <= 26) then
           new_weight = 0.25d0
           prev_weight = 0.75d0
        else
           new_weight = 0.5d0
           prev_weight = 0.5d0
        endif
        h%new = new_weight * h%new + prev_weight * hprev
        U%new = new_weight * U%new + prev_weight * Uprev
        V%new = new_weight * V%new + prev_weight * Vprev

        ! Set the value of the indices of the grid point & grid layer
        ! interface immediately below the mixing layer depth.
        call find_mixing_layer_indices()

        ! Calculate the mixing layer depth convergence metric for
        ! measuring the progress of the implicit solver
        ! Accuracy is 2% at large mixed layer 
        ! depths and about 4 cm at 1 m mixed layer depths
        h_scale = 0.02d0 * hprev + 0.02d0
        del%h = abs(h%new - hprev) / h_scale
        ! Calculate convergence metrics for the velocity field.
        ! Accuracy is 2% at large velocities and 0.7 cm/s at 10.0 cm/s
        ! (units are m/s)
        vel_scale = 0.02d0 * sqrt(Uprev(1) ** 2 + Vprev(1) ** 2) + 0.005d0
        del%U = abs(U%new(1) - Uprev(1)) / vel_scale
        del%V = abs(V%new(1) - Vprev(1)) / vel_scale

        ! Test for convergence of implicit solver
        if (count >= 2 .and. max(del%h, del%U, del%V) < 1.0d0) then
           exit
        endif
     enddo  !---------- End of the implicit solver loop ----------

     !---------- Biology Model ----------
     !
     ! Solve the biology model ODEs to advance the biology quantity values
     ! to the next time step, and calculate the growth - mortality terms
     ! (*_RHS%bio) of the semi-implicit diffusion/advection equations.
     call calc_bio_rate(time, day, dt, grid%M, precision, step_guess, step_min,  &
          T%new(0:grid%M), I_Par, P%micro, P%nano, P%pico, Z, N%O, N%H, Si,   &
          D%DON, D%PON, D%refr, D%bSi)
!write (142,*) time, S%new(158),  '% calc_bio_rate'
     ! Build the rest of the terms of the semi-implicit diffusion/advection
     ! PDEs for the biology quantities.
     !
     ! This calculates the values of the precursor diffusion
     ! coefficients matrix (Bmatrix%bio%*), the RHS diffusion/advection
     ! term vectors (*_RHS%diff_adv%new), and the RHS sinking term
     ! vectors (*_RHS%sink).
     call build_biology_equations(grid, dt, P%micro, P%nano, P%pico, Z, N%O, N%H, & ! in
          Si, D%DON, D%PON, D%refr, D%bSi, wupwell + w_wind)                        ! in

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
     call solve_bio_eqns(grid%M, P%micro, P%nano, P%pico, Z, N%O, N%H, Si, &
          D%DON, D%PON, D%refr, D%bSi, day, time)

   !
     !---------- End of Biology Model ----------

     ! Update boundary conditions at surface
     P%micro(0) = P%micro(1)
     P%nano(0) = P%nano(1)
     P%pico(0) = P%pico(1)
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
     ! For those variables that we have data for, use the annual fit
     ! calculated from the data
     call bot_bound_time(year, day, day_time, &                             ! in
          T%new(grid%M+1), S%new(grid%M+1), N%O(grid%M+1), &          ! out
          Si(grid%M+1), N%H(grid%M+1), P%micro(grid%M+1), &
          P%nano(grid%M+1), P%pico(grid%M+1), Z(grid%M+1)) ! out
     ! For those variables that we have no data for, assume uniform at
     ! bottom of domain



     call bot_bound_uniform(grid%M, D%DON, D%PON, D%refr, D%bSi)


     ! Write standard time series results
     ! !!! Please don't change this argument list without good reason. !!!
     ! !!! If it is changed, the change should be committed to CVS.    !!!
     ! !!! For exploratory, debugging, etc. output use                 !!!
     ! !!! write_user_timeseries() below.                              !!!
     call write_std_timeseries(time / 3600., grid,                    &
       ! Variables for standard physics model output
       count, h%new, U%new, V%new, T%new, S%new, I_par(0),             &
       ! Variables for standard biology model output
       N%O , N%H, Si, P%micro, P%nano, P%pico, Z, D%DON, D%PON, D%refr, D%bSi)


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
     call write_std_profiles(datetime_str(runDatetime),                  &
          datetime_str(initDatetime), year, day, day_time, dt, grid,     &
          T%new, S%new, rho%g, P%micro, P%nano, P%pico, Z, N%O, N%H, Si, &
          D%DON, D%PON, D%refr, D%bSi, K%m, K%T, K%S,                    &
          I_par, U%new, V%new)


     ! Write user-specified profiles results
     ! !!! Please don't add arguments to this call.           !!!
     ! !!! Instead put use statements in your local copy of   !!!
     ! !!! write_user_profiles() in the user_output module.   !!!
     call write_user_profiles(datetime_str(runDatetime), &
          datetime_str(initDatetime), year, day, day_time, dt, grid)



     ! Increment time, calendar date and clock time
     call new_year(day_time, day, year, time, dt)
     scount = scount + 1
     !*** Fix this to be grid independent
     sumS = sumS + 0.5 * (S%new(2) + S%new(3))
     call interp_value(1.0d0, 0, grid%d_g, S%new, Sone, junk)
     if (0.5 * (S%new(2) + S%new(3))- Sone .ne. 0) then
        write (*,*) 'Sdiff', 0.5 * (S%new(2) + S%new(3))- Sone
     endif
     sumSriv = sumSriv + S_riv
     ! Diagnostic, to check linearity of the freshwater forcing
     ! comment out for production runs
 
   
     write (129,*) S_riv,  S%new(1)

  end do  !--------- End of time loop ----------

  write (stdout,*) "For Ft tuning"
  write (stdout,*) "Average SSS should be", sumSriv / float(scount)
  write (stdout,*) "Average SSS was", sumS / float(scount)



  ! Close output files
  call timeseries_output_close
  call profiles_output_close

  close(129)

  ! Deallocate memory
  call dalloc_numerics_variables
  call dalloc_grid
  call dalloc_core_variables
  call dalloc_physics_variables
  call dalloc_biology_variables
  call dalloc_IMEX_variables

  ! Return a successful exit status to the operating system
  call exit(0)
  
end program SOG
