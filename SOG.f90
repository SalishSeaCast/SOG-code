! $Id$
! $Source$

program SOG      
  ! Coupled physical and biological model of the Strait of Georgia

  ! Utility modules:
  use io_unit_defs, only: stripped_infile
  use unit_conversions, only: KtoC
  use datetime, only: datetime_, os_datetime, calendar_date, &
       clock_time, datetime_str

  ! Inherited modules
  ! *** Goal is to make these go away
  use declarations
  use surface_forcing
  use initial_sog, only: initial_mean
  use pdf
  use IMEX_constants  

  ! Refactored modules
  use input_processor, only: init_input_processor, getpars, getpari, &
       getpard, getparl
  use timeseries_output, only: init_timeseries_output, write_timeseries, &
       timeseries_output_close
  use profiles_output, only: init_profiles_output, write_profiles, &
       profiles_output_close
  use water_properties, only: rho, alpha, beta, Cp, &
       calc_rho_alpha_beta_Cp_profiles, &
       alloc_water_props, dalloc_water_props
  use grid_mod, only: init_grid, dalloc_grid, interp_g_d, interp_i, &
       gradient_i, gradient_g
  use find_upwell, only: upwell_profile, vertical_advection
  use diffusion, only: diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux
  use fitbottom, only: bot_bound_time, bot_bound_uniform
  use do_biology_mod, only: do_biology
  use biological_mod, only: init_biology, rate_detritus, rate_det
  use forcing, only: read_variation, read_forcing, get_forcing

  ! Subroutine & function modules:
  ! (Wrapping subroutines and functions in modules provides compile-time
  !  checking of number and type of arguments - but not order!)
  ! *** These should eventually end up in refactored modules
  use Coriolis_and_pg_mod
  use define_flux_mod
  use tridiag_mod

  implicit none

  double precision:: cz, upwell
  ! Upwelling constant (tuned parameter)
  !*** read in here, used by surface_flux_sog : eventually should be local
  ! to the surface_forcing module (not be be confused with current 
  ! surface_forcing module
  real(kind=dp) :: upwell_const, S_riv, sumS=0, sumSriv=0
  integer :: scount=0

  ! Internal wave breaking eddy viscosity for momentum and scalars
  ! (tuned parameters)
  real(kind=dp) :: nu_w_m, nu_w_s

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

  ! Variables for baroclinic pressure gradient calculations
  double precision :: sumu, sumv, uprev, vprev, delu, delv
  double precision :: sumpbx, sumpby, tol
  integer :: ii       ! loop index

  ! Physical model domain parameters (should go elsewhere)
  double precision :: oLx, oLy, gorLx, gorLy

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

  ! Initialize the biology model
  call init_biology(grid%M)

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
  CALL allocate3(grid%M)
  ! Allocate memory for water property arrays
  call alloc_water_props(grid%M)

  ! Length of forcing data files (move to be read in)
  wind_n = 46056 - 8 ! with wind shifted to LST we lose 8 records
  met_n = 1918
  river_n = 1826
  call read_forcing (wind_n, met_n, river_n)
  call read_variation

  ! Read the physic model parameter values
  upwell_const = getpard("upwell_const")
  nu_w_m = getpard('nu_w_m')         ! Internal wave mixing momentum
  nu_w_s = getpard('nu_w_s')         ! Internal wave mixing scalar
  Fw_scale = getpard('Fw_scale')     ! Fresh water scale factor for river flows
  Fw_surface = getparl('Fw_surface') ! Add all fresh water on surface?
  if (.not. Fw_surface) then
     Fw_depth = getpard('Fw_depth')  ! Depth to distribute freshwater flux over
  endif

  CALL initialize ! initializes everything (biology too)

  ! Read the cruise id from stdin to use to build the file name for
  ! nutrient initial conditions data file
  cruise_id = getpars("cruise_id")
  CALL initial_mean(U, V, T, S, P, N%O%new, N%H%new, Sil%new, Detritus, &
       h%new, ut, vt, pbx, pby, &
       grid, D_bins, cruise_id)

  ! Read the iteration count limit for the physics model implicit
  ! solver
  ! *** This should go somewhere else - maybe init_physics() (yet to come)
  niter = getpari("niter")

  ! Initialize mixing, and mixed layer depths
  ! *** It would be cool to calculate these from the initial salinity
  ! *** profile read from the CTD file.
  IF(h%new < grid%d_g(1))THEN 
     h%new = grid%d_g(1)
  END IF
  h_m%new = 10.

  ! Initialize the profiles of the water column properties
  ! Density (rho), thermal expansion (alpha) and saline contraction (beta)
  ! coefficients, and specific heat capacity (Cp)
  call calc_rho_alpha_beta_Cp_profiles(T%new, S%new, rho%g, alpha%g, &
       beta%g, Cp%g)
  density%new = rho%g
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
     ! Store previous time_steps in %old components
     CALL define_sog(time_step) 

     ! get forcing data
     call get_forcing (year, day, day_time, &
          Qinter, Einter, cf_value, atemp_value, humid_value, &
          unow, vnow)

     CALL irradiance_sog(cf_value, day_time, day, &
          I, I_par, grid, jmax_i, Q_sol, euph, Qinter, P)

     DO count = 1, niter !------ Beginning of the implicit solver loop ------
        ! *** I think this is finding the depths of the grid point
        ! *** and grid interface that bound the mixed layer depth
        ! d_g(j_max_g) > h ie, hh%g is grid point just below h
        CALL find_jmax_g(h, grid)
        ! d_i(j_max_i) > h ie h%i is just above h
        CALL find_jmax_i(h, grid)

        ! *** Does this allocation really need to be done for every pass
        ! *** through the implicit loop?
        CALL allocate2(alloc_stat) ! allocate shape and K's
        DO xx = 16,17
           IF (alloc_stat(xx) /= 0) THEN
              PRINT "(A)","ALLOCATION failed.  KPP.f  xx:"
              PRINT *,xx
              STOP
           END IF
        END DO

        ! Calculate surface forcing components
        ! *** Confirm that all of these arguments are necessary
        CALL surface_flux_sog(grid%M, rho%g, w, wt_r, S%new(1),        &
             S%old(1), S_riv, T%new(0), j_gamma, I, Q_t,        &
             alpha%i(0), Cp%i(0), beta%i(0), unow, vnow, cf_value/10.,    &
             atemp_value, humid_value, Qinter,stress, &
             day, dt/h%new, h, upwell_const, upwell, Einter,       &
             u%new(1), dt, Fw_surface, Fw_scale, Ft, count) 

        ! Calculate nonturbulent heat flux profile
        ! *** Vectorize this and move it into a subroutine
        Q_n(0) = I(0) / (Cp%i(0) * rho%g(0))        
        do ii = 1, grid%M       
           Q_n(ii) = I(ii) / (Cp%i(ii) * rho%i(ii))
        enddo

        ! Calculate the nonturbulent fresh water flux profile, and its
        ! contribution to the salinity profile
        ! *** Move to a subroutine
        if (Fw_surface) then
           F_n = 0.
        else
!!$           Fw_depth = h%old
           Fw = Ft * exp(-grid%d_i / Fw_depth)
           F_n = S%new * Fw
        endif

!!$        ! Store the surface buoyancy forcing value from the previous
!!$        ! iteration so we can use it to blend with the new value to
!!$        ! help the implicit solver converge more quickly
!!$        Bf_old = Bf

        ! Calculate buoyancy profile, and surface buoyancy forcing
        CALL buoyancy(grid, T%new, S%new, h, I, F_n, w%b(0), &  ! in
             rho%g, alpha%g, beta%g, Cp%g, Fw_surface,  &  ! in
             B%new, Bf)                                         ! out

!!$        ! Blend the values of the surface buoyancy forcing from current
!!$        ! and previous iteration to help convergence
!!$        Bf = (count * Bf_old + (niter - count) * Bf) / niter

        CALL fun_constants(u_star, w_star, L_star,w, Bf, h%new)   !test conv

        CALL stability   !stable = 0 (unstable), stable = 1 (stable), stable = 2 (no forcing)  this is the stability of the water column.

        IF (u_star /= 0.)  THEN         !Wind stress /= 0.
           CALL ND_flux_profile(grid,L_star,phi)   ! define flux profiles aka (B1)
           CALL vel_scales(grid, omega, phi, u_star,L_star,h)
           ! calculates wx (13) as omega (not Omega's)
        ELSE IF (u_star == 0. .AND. Bf < 0.) THEN        !Convective unstable limit
           CALL convection_scales(grid,omega,h, w_star)   !test conv
           ! calculates wx (15) as omega (not Omega's)
        ELSE                !  No surface forcing or Bf > 0. 
           omega%m%value = 0.
           omega%s%value = 0.
           omega%m%div = 0.
           omega%s%div = 0.
           omega%m%h = 0.
           omega%m%h = 0.
        END IF

        CALL shear_diff(grid, U, V, rho, K%u%shear)  !test conv  !density instead of linear B  ! calculates ocean interior shear diffusion

        ! Calculates interior double diffusion    
        call double_diff(grid,T,S,K,alpha%i,beta%i)  !test conv

        ! Define interior diffusivity K%x%total, nu_w_m and nu_w_s
        ! constant internal wave mixing
        K%u%total = 0.0
        K%s%total = 0.0
        K%t%total = 0.0          
        do xx = 1, grid%M
           K%u%total(xx) = K%u%shear(xx) + nu_w_m + K%s%dd(xx) !test conv
           K%s%total(xx) = K%u%shear(xx) + nu_w_s + K%s%dd(xx) !test conv
           K%t%total(xx) = K%u%shear(xx) + nu_w_s + K%t%dd(xx) !test conv
        enddo

        CALL interior_match(grid, h, K%t, nu_w_s)  ! calculate nu (D5)
        CALL interior_match(grid, h, K%u, nu_w_m)
        CALL interior_match(grid, h, K%s, nu_w_s)   
        !test conv

        IF (u_star /= 0. .OR. Bf < 0.) THEN
           CALL interior_match2(omega, L_star, u_star, h, grid) !test conv
        END IF

        !Define shape functions G_shape%x

        IF (u_star /= 0. .OR. Bf < 0.) THEN
           CALL shape_parameters(K%u,omega%m, h%new, a2%m, a3%m)
           CALL shape_parameters(K%s,omega%s, h%new, a2%s, a3%s)
           CALL shape_parameters(K%t,omega%s, h%new, a2%t, a3%t) !test conv
        ELSE
           a2%m = -2.0
           a2%s = -2.0
           a2%t = -2.0
           a3%m = 1.0
           a3%s = 1.0
           a3%t = 1.0
        END IF

        G_shape%m = 0.   
        G_shape%s = 0.
        G_shape%t = 0.

        DO index = 0, h%i  !Just in Surface Layer  (h%i-1)?
           IF (grid%d_i(index) > h%new) THEN
              G_shape%m(index) = 0.   
              G_shape%s(index) = 0.
              G_shape%t(index) = 0.
           ELSE ! not sigma = d_i/h%new, a1 = 1 
              G_shape%m(index) = grid%d_i(index)/h%new*(1.0 + grid%d_i(index)/h%new*(a2%m + &
                   grid%d_i(index)/h%new*a3%m)) ! (11)
              G_shape%s(index) = grid%d_i(index)/h%new*(1.0 + grid%d_i(index)/h%new*(a2%s + &
                   grid%d_i(index)/h%new*a3%s)) 
              G_shape%t(index) = grid%d_i(index)/h%new*(1.0 + grid%d_i(index)/h%new*(a2%t + &
                   grid%d_i(index)/h%new*a3%t))
           END IF
        END DO   !test conv

        !Define K%x%ML profiles (10)

        K%u%ML = 0.0
        K%s%ML = 0.0
        K%t%ML = 0.0

        DO xx = 0, h%i  !use only up to h%i-1
           IF (grid%d_i(xx) > h%new) THEN
              K%u%ML(xx) = 0.0
              K%s%ML(xx) = 0.0
              K%t%ML(xx) = 0.0
           ELSE
              K%u%ML(xx) = h%new*omega%m%value(xx)*G_shape%m(xx)
              K%s%ML(xx) = h%new*omega%s%value(xx)*G_shape%s(xx) 
              K%t%ML(xx) = h%new*omega%s%value(xx)*G_shape%t(xx)   !test conv
           END IF
           IF (K%u%ML(xx) < 0. .OR. K%s%ML(xx) < 0. .OR. K%t%ML(xx) < 0.) THEN
              PRINT "(A)","Diffusivities < 0. See KPP.f90:  time, Jday, h%new"
              PRINT *,time,day,h%new
              PRINT "(A)","K%u%ML(xx),K%s%ML(xx),K%t%ML(xx)"
              PRINT *,K%u%ML(xx),K%s%ML(xx),K%t%ML(xx)
              STOP
           END IF
        END DO

        ! K star, enhance diffusion at grid point just below or above h (see (D6)
        CALL modify_K(grid, h, K%u) !test conv  !!must modify G_shape!!
        CALL modify_K(grid, h, K%s) !test conv  !!If G_shape is used again, use modified value!!
        CALL modify_K(grid, h, K%t) !test conv

        !!Define flux profiles, w,  and Q_t (vertical heat flux)!!
        !!Don-t need for running of the model!!  (not sure what this means)
        IF (u_star /= 0. .OR. Bf /= 0.) THEN
           CALL def_gamma(L_star, grid, w,  wt_r, h, gamma, Bf, omega) 
           ! calculates non-local transport (20)
        ELSE
           gamma%t = 0.
           gamma%s = 0.
           gamma%m = 0.
        END IF

        CALL div_interface(grid, T)
        CALL div_interface(grid, S)
        CALL div_interface(grid, U)
        CALL div_interface(grid, V)

        !defines w%x, K%x%all, K%x%old, Bf%b, and F_n 
        CALL define_flux(alpha, beta)

        ! Calculate matrix B (changed to Amatrix later) 
        call diffusion_coeff(grid, dt, K%t%all, &
             Bmatrix%t)
        call diffusion_coeff(grid, dt, K%s%all, &
             Bmatrix%s)
        call diffusion_coeff(grid, dt, K%u%all, &
             Bmatrix%u)

        ! Start calculation of RHS called Gvector (D9) (D10)
        call diffusion_nonlocal_fluxes(grid, dt, K%t%all, gamma%t, &   ! in
             w%t(0), -Q_n, T%new(grid%M+1), &                          ! in
             Gvector%t)                                                ! out
        call diffusion_nonlocal_fluxes(grid, dt, K%s%all, gamma%s, &   ! in
             w%s(0), F_n, S%new(grid%M+1), &                           ! in
             Gvector%s)                                                ! out
        call diffusion_bot_surf_flux (grid, dt, K%u%all, w%u(0), &    ! in
             U%new(grid%M+1), &                                       ! in
             Gvector%u)                                               ! out
        call diffusion_bot_surf_flux (grid, dt, K%u%all, w%v(0), &    ! in
             V%new(grid%M+1), &                                       ! in
             Gvector%v)                                               ! out

        ! Calculate profile of upwelling velocity
        call upwell_profile(grid, upwell, wupwell)
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

        ! Calculate the Coriolis and baroclinic pressure gradient
        ! components of the G vector for each velocity component
        call Coriolis_and_pg(f, dt, V%new, pbx, &
             Gvector_c%u)
        call Coriolis_and_pg(f, dt, -U%new, pby, &
             Gvector_c%v)      

        IF (time_step == 1 .AND. count  == 1) THEN

           Bmatrix_o%t%A = Bmatrix%t%A
           Bmatrix_o%t%B = Bmatrix%t%B
           Bmatrix_o%t%C = Bmatrix%t%C
           Bmatrix_o%s%A = Bmatrix%s%A
           Bmatrix_o%s%B = Bmatrix%s%B
           Bmatrix_o%s%C = Bmatrix%s%C         
           Bmatrix_o%u%A = Bmatrix%u%A
           Bmatrix_o%u%B = Bmatrix%u%B
           Bmatrix_o%u%C = Bmatrix%u%C

           Gvector_o%t = Gvector%t
           Gvector_o%s = Gvector%s
           Gvector_o%u = Gvector%u
           Gvector_o%v = Gvector%v

           Gvector_co%u = Gvector_c%u
           Gvector_co%v = Gvector_c%v
        END IF

!!!!! IMEX SCHEME !!!!   
        CALL matrix_A (Amatrix%u, Bmatrix%u) ! changed from Bmatrix to Amatrix
        CALL matrix_A (Amatrix%t, Bmatrix%t) ! now have (D9)
        CALL matrix_A (Amatrix%s, Bmatrix%s)

        ! Add Xt to H vector (D7)
        call phys_Hvector(grid%M, S%old, Gvector%s, Gvector_o%s, &  ! in
             null_vector, null_vector, Bmatrix_o%s,              &  ! in
             ! null_vector because no Coriolis or pressure gradients terms
             Hvector%s)                                             ! out
        call phys_Hvector(grid%M, T%old, Gvector%t, Gvector_o%t, &  ! in
             null_vector, null_vector, Bmatrix_o%t,              &  ! in
             ! null_vector because no Coriolis or pressure gradients terms
             Hvector%t)
        ! Add in Coriolis term (Gvector_c) and previous value to H vector (D7)
        call phys_Hvector(grid%M, U%old, Gvector%u, Gvector_o%u, &  ! in
             Gvector_c%u, Gvector_co%u, Bmatrix_o%u,             &  ! in
             Hvector%u)                                             ! out
        call phys_Hvector(grid%M, V%old, Gvector%v, Gvector_o%v, &  ! in
             Gvector_c%v, Gvector_co%v, Bmatrix_o%u,             &  ! in
             Hvector%v)                                             ! out

        ! Solves tridiagonal system
        call tridiag(Amatrix%u%A, Amatrix%u%B, Amatrix%u%C, Hvector%u, &
             U_p)
        call tridiag(Amatrix%u%A, Amatrix%u%B, Amatrix%u%C, Hvector%v, &
             V_p)
        call tridiag(Amatrix%s%A, Amatrix%s%B, Amatrix%s%C, Hvector%s, &
             S_p)
        call tridiag(Amatrix%t%A, Amatrix%t%B, Amatrix%t%C, Hvector%t, &
             T_p)

        DO yy = 1, grid%M   !remove diffusion!!!!!!!!!!  ? not sure
           U%new(yy) = U_p(yy)
           V%new(yy) = V_p(yy)
           S%new(yy) = S_p(yy)
           T%new(yy) = T_p(yy)
        END DO

        U%new(0) = U%new(1)
        V%new(0) = V%new(1)
        S%new(0) = S%new(1)   
        T%new(0) = T%new(1)

        U%new(grid%M+1) = U%new(grid%M) 
        V%new(grid%M+1) = V%new(grid%M) 
        T%new(grid%M+1) = T%old(grid%M+1)
        S%new(grid%M+1) = S%old(grid%M+1)

        ! Update the profiles of the water column properties
        ! Density (rho), thermal expansion (alpha) and saline
        ! contraction (beta) coefficients, and specific heat capacity (Cp)
        call calc_rho_alpha_beta_Cp_profiles(T%new, S%new, rho%g, alpha%g, &
             beta%g, Cp%g)
        density%new = rho%g
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

        ! Update buoyancy profile and local buoyancy frequency
        ! *** Figure out what N_2 designates (in contrast to N in Large, etal)
        ! *** (may be buoyancy frequency squared); can it go into buoyancy?
        ! *** Why not just call buoyancy again here?
        B%new = g * (alpha%g * T%new - beta%g * S%new)
        rho%grad_g = gradient_g(rho%g)
        DO xx = 1,grid%M
           N_2_g(xx) = -(g / rho%g(xx)) * rho%grad_g(xx)
        END DO

!!!Define turbulent velocity shear V_t_square for Ri_b  use last count values

!!!Calculate beta_t: need h%e and h%e%i and new w%b  w%i vertical turbulent flux of i

        CALL div_interface(grid, T)
        CALL div_interface(grid, S)

        w%b_err(0) = 0.

        DO xx = 1, grid%M           !uses K%old and T%new
           w%t(xx) = -K%t%all(xx)*(T%div_i(xx) - gamma%t(xx))
           w%s(xx) = -K%s%all(xx)*(S%div_i(xx) - gamma%s(xx))
           w%b(xx) = g * (alpha%i(xx) * w%t(xx) - beta%i(xx) * w%s(xx))   
           w%b_err(xx) = g * (alpha%grad_i(xx) * w%t(xx) &
                - beta%grad_i(xx) * w%s(xx))
        END DO


        ! beta_t is the ratio of the entrainment flux to the surface buoyancy flux -- set equal to -0.2 when convection is occuring

        if (w%b(0)>0) then
           beta_t = -0.2
        else
           beta_t = 0.
        endif

        V_t_square = 0.


        IF (beta_t < 0.) THEN  !omega, vertical velocity scale
           CALL def_v_t_sog(grid, h, N_2_g, omega%s%value, V_t_square, beta_t, L_star) !test conv  
           ! V_t_square, the turbulent velocity shear, (23)
        END IF

        ! Calculate the profile of bulk Richardson number (Large, etal
        ! (1994) eq'n(21))
        ! *** This needs to be refactored to change density to rho
        CALL define_Ri_b_sog(grid, h, surface_h, B, U, V, density, Ri_b, &
             V_t_square, N_2_g)
        ! Find the mixing layer depth by comparing Ri_b to Ri_c
        CALL ML_height_sog(grid, Ri_b, year, day, day_time, count, h_i, jmaxg)
        ! *** This block of code is never executed because
        ! *** ML_height_sog doesn't allow jamxg to exceed grid%M-3
        IF (jmaxg >= grid%M-2) THEN ! mixing down to bottom of grid
           PRINT "(A)","jmaxg >= grid%M-2. OR. h_i >= grid%d_g(grid%M-2) in KPP"
           PRINT "(A)","jmaxg,h_i"
           PRINT *,jmaxg,h_i
           PRINT "(A)","day,time"
           PRINT *,day,time
           PRINT "(A)","T%new"
           DO xx = 0,grid%M+1
              PRINT *,T%new(xx)
           END DO
           PRINT "(A)","S%new"
           DO xx = 0,grid%M+1
              PRINT *,S%new(xx)
           END DO
           PRINT "(A)","U%new"
           DO xx = 0,grid%M+1
              PRINT *,U%new(xx)
           END DO
           PRINT "(A)","V%new"
           DO xx = 0,grid%M+1
              PRINT *,V%new(xx)
           END DO
           PRINT "(A)","rho%g"

           DO xx = 0,grid%M+1
              PRINT *,rho%g(xx)
           END DO
           PRINT "(A)","Ri_b"
           PRINT *,Ri_b
           PRINT "(A)","K%u%all"
           PRINT *,K%u%all
           PRINT "(A)","h%old"
           PRINT *,h%old
           STOP
        END IF
        !  Apply the Ekman depth criterion to the mixing layer depth
        ! when stable forcing exists
        ! *** This code could probably go into ML_height_sog.
        h_i = (h_i + h%new)/2.0 !use the average value!!!!!!!!##
        IF (stable == 1) THEN          !Stability criteria 
           h_Ekman = 0.7*u_star/f         ! see (24) and surrounding
           IF (h_Ekman < L_star) THEN
              IF (h_i > h_Ekman) THEN
                 h_i = h_Ekman
                 write (*,*) 'Ekman', h_i
              END IF
           ELSE
              IF (h_i > L_star) THEN
                 h_i = L_star
              END IF
           END IF
        END IF

        IF (h_i < grid%d_g(1)) THEN  !minimum mixing
           h_i = grid%d_g(1)
        END IF

        CALL find_jmax_g(h,grid) !***! can't be less than the grid

        del = ABS(h_i - h%new)/grid%i_space(jmaxg)  !##
        h%new = h_i
        CALL find_jmax_g(h,grid)
        CALL define_hm_sog(grid,S%new,h_m) !***! definition of mixed layer depth

        !----pressure grads--------------------

        sumu = 0.
        sumv = 0.
        DO yy = 1, grid%M   !remove barotropic mode
           sumu = sumu+U%new(yy)
           sumv = sumv+V%new(yy)
        END DO
        sumu = sumu/grid%M
        sumv = sumv/grid%M

        DO yy = 1, grid%M   !remove barotropic mode
           U%new(yy) = U%new(yy)-sumu
           V%new(yy) = V%new(yy)-sumv
        END DO

        delu = abs(U%new(1)-uprev)/0.01
        delv = abs(V%new(1)-vprev)/0.01

        uprev = U%new(1)
        vprev = V%new(1)

        if (del.lt.delu) del = delu
        if (del.lt.delv) del = delv

        ! integrate ut and vt (note v is 305 degrees and u is 35 degrees)

        oLx = 2./(20e3)
        oLy = 2./(60e3)
        gorLx = 9.8/(1025.*20e3)
        gorLy = 9.8/(1025.*60e3)

        sumu = 0
        sumv = 0
        do yy=1,grid%M
           ut%new(yy) = ut%old(yy)*0.95+U%new(yy)*dt*oLx
           vt%new(yy) = vt%old(yy)*0.95+V%new(yy)*dt*oLy 
        enddo

        dzx(1) = ut%new(1)+1
        dzy(1) = vt%new(1)+1

        do yy=2,grid%M
           dzx(yy) = dzx(yy-1) + (ut%new(yy)+1)
           dzy(yy) = dzy(yy-1) + (vt%new(yy)+1)
        enddo

        ! Calculate the baroclinic pressure gradient
        ! *** This tolerance should be read in as a run parameter
        tol=1e-6
        sumpbx = 0.
        cz = 0.
        ii=1
        do yy=1,grid%M
           if (yy == 1) then
              pbx(yy) = -rho%g(yy)
           else
              pbx(yy) = pbx(yy-1)-rho%g(yy)
           endif
           do while ((dzx(ii)-yy) < -tol .and. ii < grid%M)
              pbx(yy) = pbx(yy) + rho%g(ii)*(dzx(ii)-cz)
              cz = dzx(ii)
              ii = ii + 1
           enddo
           pbx(yy) = pbx(yy) + rho%g(ii)*(yy-cz)
           sumpbx = sumpbx + pbx(yy)
           cz = yy
        enddo

        sumpby = 0.
        cz = 0.
        ii=1
        do yy=1,grid%M
           if (yy == 1) then
              pby(yy) = -rho%g(yy)
           else
              pby(yy) = pby(yy-1)-rho%g(yy)
           endif
           do while ((dzy(ii)- yy) <-tol .and. ii < grid%M)
              pby(yy) = pby(yy) + rho%g(ii)*(dzy(ii)-cz)
              cz = dzy(ii)
              ii = ii + 1
           enddo
           pby(yy) = pby(yy) + rho%g(ii)*(yy-cz)
           sumpby = sumpby + pby(yy)
           cz = yy
        enddo

        sumpbx = sumpbx/grid%M
        sumpby = sumpby/grid%M

        do yy = 1, grid%M
           pbx(yy) = (pbx(yy) - sumpbx) * gorLx * grid%i_space(yy) &
                + stress%u%new / (1025. * grid%M * grid%i_space(yy))
           pby(yy) = (pby(yy) - sumpby) * gorLy * grid%i_space(yy) &
                + stress%v%new / (1025. * grid%M * grid%i_space(yy))     
        enddo

        ! Test for convergence of implicit solver
        if (count >= 2 .and. del < del_o) then
           exit
        endif
     enddo  !---------- End of the implicit solver loop ----------

     ! Write time series results
     call write_timeseries(time / 3600., &
       ! Variables for standard physical model output
       count, h%new, T%new, S%new, &
       ! User-defined physical model output variables
!!$       &
       ! Variables for standard biological model output
       N%O%new , N%H%new, Sil%new, P%micro%new, P%nano%new, &
       Detritus(1)%D%new, Detritus(2)%D%new, Detritus(3)%D%new &
       ! User-defined biological model output variables
!!$       &
       )

     call write_profiles(codeId, datetime_str(runDatetime),            &
          datetime_str(startDatetime), year, day, day_time, dt, grid,  &
          T%new, S%new, rho%g, P%micro%new, P%nano%new, N%O%new, &
          N%H%new, Sil%new, Detritus(1)%D%new, Detritus(2)%D%new,      &
          Detritus(3)%D%new, K%u%all, K%t%all, K%s%all, I_par, U%new,  &
          V%new)

!------BIOLOGICAL MODEL--------------------------------------------

     call do_biology (time, day, dt, grid%M, precision, step_guess, step_min, &
          T%new(0:grid%M), I_Par, P, N, Sil, Detritus, &
          Gvector_ro)
!*** more of the below can be moved into the do_biology module

     ! Calculate diffusion matrix Bmatrix
     call diffusion_coeff (grid, dt, K%s%all, &
          Bmatrix%bio)
     Bmatrix%no = Bmatrix%bio          ! both diffuse like S

     ! Initialize Gvector and add the upward from the bottom (surface flux = 0)
     call diffusion_bot_surf_flux (grid, dt, K%s%all, 0.d0, &          ! in
          P%micro%new(grid%M+1), &                                     ! in
          Gvector%p%micro)                                             ! out
     call diffusion_bot_surf_flux (grid, dt, K%s%all, 0.d0, &          ! in
          P%nano%new(grid%M+1), &                                      ! in
          Gvector%p%nano)                                              ! out
     DO xx = 1, D_bins-1
        call diffusion_bot_surf_flux (grid, dt, K%s%all, 0.d0, &       ! in
             Detritus(xx)%d%new(grid%M+1), &                           ! in
             Gvector%d(xx)%bin)                                        ! out
     END DO
     Gvector%d(D_bins)%bin = 0.
     call diffusion_bot_surf_flux (grid, dt, K%s%all, 0.d0, &          ! in
          N%O%new(grid%M+1), &                                         ! in
          Gvector%n%o)                                                 ! out
     call diffusion_bot_surf_flux (grid, dt, K%s%all, 0.d0, &          ! in
          N%H%new(grid%M+1), &                                         ! in
          Gvector%n%h)                                                 ! out
     call diffusion_bot_surf_flux (grid, dt, K%s%all, 0.d0, &          ! in
          Sil%new(grid%M+1), &                                         ! in
          Gvector%sil)                                                 ! out

     ! vertical advection
     call vertical_advection (grid, dt, P%micro%new, wupwell, &
          Gvector%p%micro)
     call vertical_advection (grid, dt, P%nano%new, wupwell, &
          Gvector%p%nano)
     call vertical_advection (grid, dt, N%O%new, wupwell, &
          Gvector%n%o)
     call vertical_advection (grid, dt, N%H%new, wupwell, &
          Gvector%n%h)
     call vertical_advection (grid, dt, Sil%new, wupwell, &
          Gvector%sil)
     DO xx = 1,D_bins-1
        call vertical_advection (grid, dt, Detritus(xx)%D%new, wupwell, &
             Gvector%d(xx)%bin)
     END DO

     micro%sink = 0.*1.1574D-05 ! 0/m per day
!     CALL advection(grid,micro%sink,P%micro%old,dt,Gvector_ao%p%micro) !sinking phytoplankton
     Gvector_ao%p%micro = 0. ! sinking not implemented

     DO xx = 1,D_bins-1
        CALL advection(grid,rate_det%sink(xx),Detritus(xx)%D%old,dt,Gvector_ao%d(xx)%bin)
     END DO
     Gvector_ao%d(D_bins)%bin = 0.

     CALL matrix_A (Amatrix%bio, Bmatrix%bio)   !define Amatrix%A,%B,%C
     CALL matrix_A (Amatrix%no, Bmatrix%no)

     IF (time_step == 1) THEN ! initial estimate is better than 0.
        Bmatrix%null%A = 0. !no diffusion
        Bmatrix%null%B = 0.
        Bmatrix%null%C = 0.
        Amatrix%null%A = 0. !(M)
        Amatrix%null%B = 1.
        Bmatrix_o%bio%A = Bmatrix%bio%A
        Bmatrix_o%bio%B = Bmatrix%bio%B
        Bmatrix_o%bio%C = Bmatrix%bio%C
        Bmatrix_o%no%A = Bmatrix%no%A
        Bmatrix_o%no%B = Bmatrix%no%B
        Bmatrix_o%no%C = Bmatrix%no%C
        Gvector_o%p%micro = Gvector%p%micro 
        Gvector_o%p%nano = Gvector%p%nano !V.flagella.01
        Gvector_o%n%o = Gvector%n%o
        Gvector_o%n%h = Gvector%n%h
        Gvector_o%sil = Gvector%sil

        DO xx = 1,D_bins
           Gvector_o%d(xx)%bin = Gvector%d(xx)%bin
        END DO

     END IF ! time_step == 1

     ! Build the H vectors for the biological quantities
     CALL P_H (grid%M, P%micro%old, Gvector%p%micro, Gvector_o%p%micro, &
          Gvector_ro%p%micro, Gvector_ao%p%micro, Bmatrix_o%bio, &
          Hvector%p%micro)
     CALL P_H(grid%M, P%nano%old, Gvector%p%nano, Gvector_o%p%nano, &
          Gvector_ro%p%nano, null_vector, Bmatrix_o%bio, &
          Hvector%p%nano) ! null_vector 'cause no sinking
     DO xx = 1, D_bins - 1
        CALL P_H (grid%M, Detritus(xx)%D%old, Gvector%d(xx)%bin, &
             Gvector_o%d(xx)%bin, Gvector_ro%d(xx)%bin, &
             Gvector_ao%d(xx)%bin, Bmatrix_o%bio, &
             Hvector%d(xx)%bin)
     END DO
     CALL P_H(grid%M, Detritus(D_bins)%D%old, Gvector%d(D_bins)%bin, &
          Gvector_o%d(D_bins)%bin, Gvector_ro%d(D_bins)%bin, &
          null_vector, Bmatrix%null, Hvector%d(D_bins)%bin) 
     CALL P_H(grid%M, N%O%old, Gvector%n%o, Gvector_o%n%o, &
          Gvector_ro%n%o, null_vector, Bmatrix_o%no, &
          Hvector%n%o)  ! null_vector 'cause no sinking
     CALL P_H(grid%M,  N%H%old, Gvector%n%h, Gvector_o%n%h, &
          Gvector_ro%n%h, null_vector, Bmatrix_o%no, &
          Hvector%n%h)  ! null_vector 'cause no sinking
     CALL P_H(grid%M,  Sil%old, Gvector%sil, Gvector_o%sil, &
          Gvector_ro%sil, null_vector, Bmatrix_o%no, &
          Hvector%sil)  ! null_vector 'cause no sinking

     ! Solve the tridiagonal system for the biological quantities
     call tridiag(Amatrix%bio%A, Amatrix%bio%B, Amatrix%bio%C, &
          Hvector%p%micro, P1_p)
     call tridiag(Amatrix%bio%A, Amatrix%bio%B, Amatrix%bio%C, &
          Hvector%p%nano, Pnano1_p)
     call tridiag(Amatrix%no%A, Amatrix%no%B, Amatrix%no%C, Hvector%n%o, &
          NO1_p) 
     call tridiag(Amatrix%no%A, Amatrix%no%B, Amatrix%no%C, Hvector%n%h, &
          NH1_p) 
     call tridiag(Amatrix%no%A, Amatrix%no%B, Amatrix%no%C, Hvector%sil, &
          Sil1_p)
     do xx = 1,D_bins-1
        call tridiag(Amatrix%bio%A, Amatrix%bio%B, Amatrix%bio%C, &
             Hvector%d(xx)%bin, Detritus1_p(xx,:))
     enddo
     call tridiag(Amatrix%null%A, Amatrix%null%B, Amatrix%null%A, &
          Hvector%d(D_bins)%bin, Detritus1_p(D_bins,:))

     CALL find_new(grid%M)

     !-----END BIOLOGY------------------------------------------------

     !--------bottom boundaries--------------------------
     ! Update boundary conditions at bottom of grid
     !
     ! For those variables that we have data, use the annual fit
     ! calculated from the data
     call bot_bound_time (day, day_time, &                                ! in
          T%new(grid%M+1), S%new(grid%M+1), N%O%new(grid%M+1), &          ! out
          Sil%new(grid%M+1), P%micro%new(grid%M+1), P%nano%new(grid%M+1)) ! out
     ! For those variables that we have no data for, assume uniform at
     ! bottom of domain
     call bot_bound_uniform (grid%M, N%H%new, Detritus)

     ! Increment time, calendar date and clock time
     call new_year(day_time, day, year, time, dt, month_o)
     scount = scount + 1
     sumS = sumS + S%new(1)
     sumSriv = sumSriv + S_riv

  end do  !--------- End of time loop ----------

  write (*,*) "For Ft tuning"
  write (*,*) "Average SSS should be", sumSriv/float(scount)
  write (*,*) "Average SSS was", sumS/float(scount)

  ! Close output files
  call timeseries_output_close
  call profiles_output_close

  ! Deallocate memory
  call dalloc_water_props
  call dalloc_grid

end program SOG
