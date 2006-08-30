! $Id$
! $Source$

program SOG      
  ! Coupled physical and biological model of the Strait of Georgia

  ! Utility modules:
  use io_unit_defs, only: profiles
  use unit_conversions, only: KtoC
  use datetime, only: datetime_, os_datetime, calendar_date, &
       clock_time, datetime_str

  use declarations
  use surface_forcing
  use initial_sog
  use pdf
  use IMEX_constants  

  use water_properties, only: water_property, alloc_water_props, &
       dalloc_water_props, Cp_profile
  use grid_mod, only: init_grid, dalloc_grid, interp_i
  use find_upwell, only: upwell_profile, vertical_advection
  use diffusion, only: diffusion_coeff, diffusion_nonlocal_fluxes, &
       diffusion_bot_surf_flux
  ! Subroutine & function modules:
  ! (Wrapping subroutines and functions in modules provides compile-time
  !  checking of number and type of arguments - but not order!)
  use find_wind_mod
  use Coriolis_and_pg_mod
  use define_flux_mod

  implicit none

  external derivs_noflag, derivs_sog, rkqs

  ! Local variables:
  integer :: icheck
  integer :: ecmapp, day_met

  type(bins) :: PZ_bins  ! where quantities (eg. phyto, nitrate) are in PZ
  common /derivs/ PZ_bins
  integer :: bPZ, ePZ ! start position and end position in PZ array


  double precision:: cz, unow, vnow, upwell 

  ! Water column physical properties
  type(water_property) :: Cp    ! Heat capacity in J/kg-K

  ! Interpolated river flows
  real :: Qinter  ! Fraser River
  real :: Einter  ! Englishman River
  ! Filename returned by getpars()
  character*80 :: str	
  ! Initial month parameter read from run control input file
  ! and passed to new_year()
  integer :: month_o
  ! Iteration limit for inner loop that solves KPP turbulence model
  integer :: niter
  ! Hour of day index for cloud fraction data
  integer :: j
  ! Variables for baroclinic pressure gradient calculations
  double precision :: sumu, sumv, uprev, vprev, delu, delv
  double precision :: sumpbx, sumpby, tol
  integer :: ii       ! loop index
  ! Physical model domain parameter (should go elsewhere)
  double precision :: oLx, oLy, gorLx, gorLy
  ! Loop index for writing out profile results
  integer :: i_pro
  ! sigma-t quantity calculated from density for profile results output
  double precision :: sigma_t
  ! Code identification string (maintained by CVS), for output file headers
  character*70 :: &
       codeId = "$Source$"
  ! Date/time structures for output file headers
  type(datetime_) :: runDatetime     ! Date/time of code run
  type(datetime_) :: profileDatetime ! Date/time of profile
  ! Temporary storage for formated datetime strings.  Needed to work around
  ! an idiocyncracy in pgf90 that seems to disallow non-intrinsic function
  ! calls in write statements
  character*19 :: str_runDatetime, str_proDatetime
  ! *** Temporary flag to turn flagellates model on/off, and the
  ! *** function to read it
  logical :: flagellates
  logical :: getparl
  external getparl

  ! Get the current date/time from operating system to timestamp the
  ! run with
  call os_datetime(runDatetime)

  ! *** This needs to be expanded into a real input processor
  ! Initialize the parameter reader to output a report
  ! *** runDatetime and date/time that getpar_init() prints should be same
  call getpar_init(1)
  ! Read a flag to determine whether or not flagellates are included
  ! in the model
  flagellates = getparl('flagellates_on', 1)
  ! Open the time series output files
  CALL write_open
  ! Get the name of the main run parameters file from stdin, open it,
  ! and read them
  str = getpars("inputfile", 1)
  open(10, file=str, status="OLD", action="READ")
  read(10, *) M, D, lambda , t_o, t_f, dt, day_o, year_o, month_o
  cruise_id = getpars("cruise_id", 1)


  steps = 1 + int((t_f - t_o) / dt) !INT rounds down

  ! *** Why not read M & D directly into grid?
  ! *** D has an implicit type conversion problem too
  grid%M = M
  grid%D = D
  ! *** These constants should be set as parameter somewhere else, or read
  ! *** from the main run parameters file
  wind_n = 46056 - 8 ! with wind shifted to LST we lose 8 records
  stable = 1

  ! Size of the biology we are using (Quantities and Detritus)
  PZ_bins%Quant = 4
  ! Position of Diatoms (micro plankton)
  PZ_bins%micro = 1
  ! Position of Flagellates (nano plankton)
  PZ_bins%nano = 2
  ! Position of Nitrate
  PZ_bins%NO = 3
  ! Position of Ammonium
  PZ_bins%NH = 4
  ! Start of detritus
  PZ_bins%det = 5
  ! Number of detritus bins, dissolved, slow sink and fast sink
  D_bins = 3
  M2 = (PZ_bins%Quant+D_bins)*M   !size of PZ in biology: 

  IF (year_o==2001) then
     ecmapp = 1
  else if (year_o == 2002) then
     ecmapp = 8761
  else if (year_o == 2003) then
     ecmapp= 17521
  else if (year_o==2004) then
     ecmapp = 26281
  else if (year_o==2005) then
     ecmapp = 35065
  else if (year_o==2006) then
     ecmapp = 43825
  endif

  !print*,day,ecmapp,year_o,'sog 2'
  !pause

  icheck=346

  CALL allocate1(alloc_stat) 
  DO xx = 1,12
     IF (alloc_stat(xx) /= 0) THEN
        PRINT "(A)","ALLOCATION failed.  KPP.f  xx:"
        PRINT *,xx,alloc_stat(xx)
        STOP
     END IF
  END DO
  call alloc_water_props(grid%M, Cp)
  CALL read_sog
  CALL coefficients(alph, beta, dens, cloud,p_Knox)
  CALL allocate3

  CALL initialize ! initializes everything (biology too)
  CALL define_grid(grid, D, lambda) ! sets up the grid
  ! Set up grid interpolation factor arrays
  ! *** Eventually, define_grid will be refactored into init_grid
  call init_grid(M, grid)
  CALL initial_mean(U, V, T, S, P, N%O%new, N%H%new, Sil%new, Detritus, &
       h%new, ut, vt, pbx, pby, &
       grid, D_bins, cruise_id, flagellates)

  max_length = M2   !      max_length = MAXVAL(Cevent%length) Amatrix...

  !---End Biology from KPP--------------------------------------

  IF(h%new < grid%d_g(1))THEN 
     h%new = grid%d_g(1)
  END IF

  !    initialize h_m
  h_m%new = 10.

  ! Initialize the profiles of the water column properties
  ! Thermal and salinity expansion and their gradients
  CALL alpha_sub(T%new, S%new, alph, grid) 
  CALL alpha_sub(T%new, S%new, beta, grid)
  CALL div_i_param(grid, alph) ! alph%idiv = d alpha /dz
  CALL div_i_param(grid, beta) ! beta%idiv = d beta /dz
  ! Heat capacity
  Cp%g = Cp_profile(T%new, S%new)
  Cp%i = interp_i(Cp%g)
  ! Density with depth and density of fresh water
  CALL density_sub(T, S, density%new, M,rho_fresh_o) 
  
  do time_step = 1, steps  !---------- Beginning of the time loop ----------
     ! Store previous time_steps in old (n) and old_old (n-1)
     CALL define_sog(time_step) 

     ! Iteration limit for implicit solver that calculates physics
     ! *** This should be read from the run parameters data file
     niter = 30
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

        ! Initialization of the implicit loop
        IF (count == 1) THEN
           ! *** What do these comments mean?
           !only starts with zooplankton!! <- no idea what this means 
           !        
           ! fixed cloud_type

           ! *** Implicit type conversion problem
           j = day_time / 3600.0 + 1
           ! *** There has to be a better way...
           IF (year==2001) then
              day_met=day
           else if (year==2002) then
              day_met=day+365
           else if (year==2003) then
              day_met=day+730
           else if (year==2004) then
              day_met=day+1095
           else if (year==2005) then
              day_met=day+1461
           else if (year==2006) then
              day_met=day+1826
           endif

           ! Interpolate river flows for the second we're at
           ! *** Implicit type conversion problem here!!!
           Qinter = (day_time * Qriver(day_met) &
                + (86400. - day_time) * Qriver(day_met-1)) / 86400.
           Einter = (day_time * Eriver(day_met) &
                + (86400. - day_time) * Eriver(day_met-1)) / 86400.

           CALL irradiance_sog(cloud, cf(day_met, j), day_time, day, &
                I, I_par, grid, jmax_i, Q_sol, euph, Qinter, h, P)
        endif  !------ End of implicit loop initialization ------

        ! Br radiative contribution to the surface buoyancy forcing
        ! B%new buoyancy profile B(d)
        ! Q_n nonturbulent heat flux profile
        CALL buoyancy(alph, T%new, S%new, grid, h, B%new, I, Br, &
             density%new, Cp, beta, Q_n)

        CALL div_grid(grid, density)
        ! take density from grid points to interface or vice-versa

        ! Search and interpolate the wind data for the u and v components
        ! of the wind at the current time step
        call find_wind(year, day, time, ecmapp, wind_n, wind, unow, vnow)

        ! *** Implicit type conversion problem
        ! *** Identical statement to this above - why???
        j = day_time / 3600.0 + 1

        CALL surface_flux_sog(grid%M,density%new,w,wt_r,S%old(1),S%new(h%i),S%new(M),T%new(0),j_gamma, &
             I, Q_t(0), alph%i(0), Cp%i(0), beta%i(0),unow, vnow, cf(day_met,j)/10., atemp(day_met,j), humid(day_met,j), &
             Qinter,stress, rho_fresh_o,day,dt/grid%i_space(1),h,upwell,Einter,u%new(1), dt) 

        Bf%b(0) = -w%b(0)+Br   !surface buoyancy forcing *nonturbulent heat flux beta*F_n would also go here  Br is radiative contribution

        IF (time_step == 1 .AND. count == 1) THEN 
           !!Stores previous time_step for initial loop
           B%old = B%new   
           density%old = density%new

        END IF

        CALL fun_constants(u_star, w_star, L_star,w, Bf%b(0), h%new)   !test conv

        CALL stability   !stable = 0 (unstable), stable = 1 (stable), stable = 2 (no forcing)  this is the stability of the water column.

        IF (u_star /= 0.)  THEN         !Wind stress /= 0.
           CALL ND_flux_profile(grid,L_star,phi)   ! define flux profiles aka (B1)
           CALL vel_scales(grid, omega, phi, u_star,L_star,h)
           ! calculates wx (13) as omega (not Omega's)
        ELSE IF (u_star == 0. .AND. Bf%b(0) < 0.) THEN        !Convective unstable limit
           CALL convection_scales(grid,omega,h, w_star)   !test conv
           ! calculates wx (15) as omega (not Omega's)
        ELSE                !  No surface forcing or Bf%b(0) > 0. 
           omega%m%value = 0.
           omega%s%value = 0.
           omega%m%div = 0.
           omega%s%div = 0.
           omega%m%h = 0.
           omega%m%h = 0.
        END IF

        CALL shear_diff(grid,U,V,density,K%u%shear)  !test conv  !density instead of linear B  ! calculates ocean interior shear diffusion

        CALL double_diff(grid,T,S,K,alph%i,beta%i)  !test conv
        ! calculates interior double diffusion    

        !Define interior diffusivity K%x%total, nu_w_m and nu_w_s constant internal wave mixing

        K%u%total = 0.0
        K%s%total = 0.0
        K%t%total = 0.0          

        DO xx = 1,M

           K%u%total(xx) = K%u%shear(xx) + nu_w_m + K%s%dd(xx) !test conv
           K%s%total(xx) = K%u%shear(xx) + nu_w_s + K%s%dd(xx) !test conv
           K%t%total(xx) = K%u%shear(xx) + nu_w_s + K%t%dd(xx) !test conv          
        END DO

        CALL interior_match(grid, h, K%t, nu_w_s)  ! calculate nu (D5)
        CALL interior_match(grid, h, K%u, nu_w_m)
        CALL interior_match(grid, h, K%s, nu_w_s)   
        !test conv

        IF (u_star /= 0. .OR. Bf%b(0) < 0.) THEN
           CALL interior_match2(omega, L_star, u_star, h, grid) !test conv
        END IF

        !Define shape functions G_shape%x

        IF (u_star /= 0. .OR. Bf%b(0) < 0.) THEN
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
        IF (u_star /= 0. .OR. Bf%b(0) /= 0.) THEN
           CALL def_gamma(L_star, grid, w,  wt_r, h, gamma, Bf%b(0), omega) 
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
        CALL define_flux(Cp%i)

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
        call upwell_profile(grid, S%new, upwell, wupwell)
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

        ! add Xt to H vector (D7)

        CALL scalar_H(grid,Hvector%s,Gvector%s,Gvector_o%s,Gvector_o_o%s,Bmatrix_o%s,Bmatrix_o_o%s,&
             time_step,S)
        CALL scalar_H(grid,Hvector%t,Gvector%t,Gvector_o%t,Gvector_o_o%t,Bmatrix_o%t,Bmatrix_o_o%t,&
             time_step,T)

        ! add in Coriolis term (Gvector_c) and previous value to H vector (D7)
        CALL U_H(M, U%old, Gvector%u, Gvector_o%u, Gvector_c%u, &    ! in
             Gvector_co%u, Bmatrix_o%u, &                            ! in
             Hvector%u)                                              ! out
        CALL U_H(M, V%old, Gvector%v, Gvector_o%v, Gvector_c%v, &    ! in
             Gvector_co%v, Bmatrix_o%u, &                            ! in
             Hvector%v)                                              ! out
        ! solves tridiagonal system
        CALL TRIDAG(Amatrix%u%A,Amatrix%u%B,Amatrix%u%C,Hvector%u,U_p,M)
        CALL TRIDAG(Amatrix%u%A,Amatrix%u%B,Amatrix%u%C,Hvector%v,V_p,M)
        CALL TRIDAG(Amatrix%s%A,Amatrix%s%B,Amatrix%s%C,Hvector%s,S_p,M)
        CALL TRIDAG(Amatrix%t%A,Amatrix%t%B,Amatrix%t%C,Hvector%t,T_p,M)

        DO yy = 1, M   !remove diffusion!!!!!!!!!!  ? not sure
           U%new(yy) = U_p(yy)
           V%new(yy) = V_p(yy)
           S%new(yy) = S_p(yy)
           T%new(yy) = T_p(yy)
        END DO

        U%new(0) = U%new(1)
        V%new(0) = V%new(1)
        S%new(0) = S%new(1)   
        T%new(0) = T%new(1)

        U%new(M+1) = U%new(M) 
        V%new(M+1) = V%new(M) 
        T%new(M+1) = T%old(M+1)
        S%new(M+1) = S%old(M+1)

        ! Update the profiles of the water column properties
        ! Thermal and salinity expansion and their gradients
        CALL alpha_sub(T%new, S%new, alph, grid)
        CALL alpha_sub(T%new, S%new, beta, grid)
        ! Heat capacity
        Cp%g = Cp_profile(T%new, S%new)
        Cp%i = interp_i(Cp%g)
        
        ! Density with depth and density of fresh water
        CALL density_sub(T, S, dens_i, M,rho_fresh_o)
        density%new = dens_i

        CALL div_i_param(grid,alph)
        CALL div_i_param(grid,beta)

        B%new = g*(alph%g*T%new - beta%g*S%new) ! update buoyancy gradient

        CALL div_grid(grid,density)
        CALL div_grid(grid,B)

        DO xx = 1,grid%M
           N_2_g(xx) = -g/density%new(xx)*density%div_g(xx) !B%div_g(xx)
        END DO

!!!Define turbulent velocity shear V_t_square for Ri_b  use last count values

!!!Calculate beta_t: need h%e and h%e%i and new w%b  w%i vertical turbulent flux of i

        CALL div_interface(grid, T)
        CALL div_interface(grid, S)

        w%b_err(0) = 0.

        DO xx = 1, M           !uses K%old and T%new
           w%t(xx) = -K%t%all(xx)*(T%div_i(xx) - gamma%t(xx))
           w%s(xx) = -K%s%all(xx)*(S%div_i(xx) - gamma%s(xx))
           w%b(xx) = g*(alph%i(xx)*w%t(xx)-beta%i(xx)*w%s(xx))   
           w%b_err(xx) = g*(alph%idiv(xx)*w%t(xx)- beta%idiv(xx)*w%s(xx))
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

        CALL define_Ri_b_sog(grid, h, surface_h, B, U, V, density, Ri_b, V_t_square, N_2_g)
        ! (21)
        CALL ML_height_sog(grid, Ri_b, h_i, jmaxg) !test conv
        ! (21) tested for minimum value of d at which Ri_b = Ri_c

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
           PRINT "(A)","density%new"

           DO xx = 0,grid%M+1
              PRINT *,density%new(xx)
           END DO
           PRINT "(A)","Ri_b"
           PRINT *,Ri_b
           PRINT "(A)","K%u%all"
           PRINT *,K%u%all
           PRINT "(A)","h%old"
           PRINT *,h%old
           STOP
        END IF

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
        DO yy = 1, M   !remove barotropic mode
           sumu = sumu+U%new(yy)
           sumv = sumv+V%new(yy)
        END DO
        sumu = sumu/M
        sumv = sumv/M

        DO yy = 1, M   !remove barotropic mode
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
        do yy=1,M
           ut%new(yy) = ut%old(yy)+U%new(yy)*dt*oLx
           vt%new(yy) = vt%old(yy)+V%new(yy)*dt*oLy 
        enddo

        dzx(1) = ut%new(1)+1
        dzy(1) = vt%new(1)+1

        do yy=2,M
           dzx(yy) = dzx(yy-1) + (ut%new(yy)+1)
           dzy(yy) = dzy(yy-1) + (vt%new(yy)+1)
        enddo

        ! Calculate the baroclinic pressure gradient
        ! *** This tolerance should be read in as a run parameter
        tol=1e-6
        sumpbx = 0.
        cz = 0.
        ii=1
        do yy=1,M
           if (yy == 1) then
              pbx(yy) = -density%new(yy)
           else
              pbx(yy) = pbx(yy-1)-density%new(yy)
           endif
           do while ((dzx(ii)-yy) < -tol .and. ii < M)
              pbx(yy) = pbx(yy) + density%new(ii)*(dzx(ii)-cz)
              cz = dzx(ii)
              ii = ii + 1
           enddo
           pbx(yy) = pbx(yy) + density%new(ii)*(yy-cz)
           sumpbx = sumpbx + pbx(yy)
           cz = yy
        enddo

        sumpby = 0.
        cz = 0.
        ii=1
        do yy=1,M
           if (yy == 1) then
              pby(yy) = -density%new(yy)
           else
              pby(yy) = pby(yy-1)-density%new(yy)
           endif
           do while ((dzy(ii)- yy) <-tol .and. ii < M)
              pby(yy) = pby(yy) + density%new(ii)*(dzy(ii)-cz)
              cz = dzy(ii)
              ii = ii + 1
           enddo
           pby(yy) = pby(yy) + density%new(ii)*(yy-cz)
           sumpby = sumpby + pby(yy)
           cz = yy
        enddo

        sumpbx = sumpbx/M
        sumpby = sumpby/M

        do yy = 1, M
           pbx(yy) = (pbx(yy) - sumpbx) * gorLx * grid%i_space(yy) &
                + stress%u%new / (1025. * M * grid%i_space(yy))
           pby(yy) = (pby(yy) - sumpby) * gorLy * grid%i_space(yy) &
                + stress%v%new / (1025. * M * grid%i_space(yy))     
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


        IF (count >= 2 .AND. del < del_o) THEN
           IF (count == 2) THEN
              count_two = count_two + 1         
           END IF
           count_tot = count_tot + 1
           CALL write_physical_sog(unow,vnow,euph)
           EXIT
        ELSE IF (count == niter .AND. del >= del_o) THEN
           count_no =  count_no + 1 
           CALL write_physical_sog (unow,vnow,euph)
        END IF

     ENDDO

!------BIOLOGICAL MODEL--------------------------------------------


! load the PZ vector with all the biological quantities
     call define_PZ(M, PZ_bins, D_bins, M2, &                    !in
          P%micro%new, P%nano%new, N%O%new, N%H%new, Detritus, & !in
          PZ)                                                    !out

     IF (MINVAL(PZ) < 0.) THEN
        PRINT "(A)","PZ < 0. see SOG.f90"
        PRINT "(A)","time,day"
        PRINT *,time,day
        STOP
     END IF

     next_time = time+dt ! note, biology is calculated for the NEXT step

     ! *** Temporary code to allow flagellates model to be turned on/off
     !*** Size of T in odeint is hard-coded to 81
     if (flagellates) then
        call odeint(PZ, M2, time, next_time, precision, step_guess, step_min, &
             N_ok, N_bad, derivs_sog, rkqs, icheck, T%new(0:M))
     else
        call odeint(PZ, M2, time, next_time, precision, step_guess, step_min, &
             N_ok, N_bad, derivs_noflag, rkqs, icheck, T%new(0:M))
     endif

     ! check for negative NH values and then for negative Micro phyto values
     !*** add nanos and move into a subroutine in bio module 
     bPZ = (PZ_bins%NH - 1) * M + 1
     ePZ = PZ_bins%NH * M
     IF (MINVAL(PZ(bPZ:ePZ)) < 0.) THEN
        DO xx = bPZ,ePZ
           IF (PZ(xx) < 0.) THEN
              PRINT "(A)","N%H%new(xx-4*M) < 0."
              PRINT *,PZ(xx)
              PRINT "(A)","set equal to zero"
              PZ(xx) = 0.
           END IF
        END DO
     END IF
     bPZ = (PZ_bins%micro - 1) * M + 1
     ePZ = PZ_bins%micro * M
     IF (MINVAL(PZ(bPZ:ePZ)) < 0.) THEN
        PRINT "(A)","PZ < 0. After odeint.f see SOG.f90"
        PRINT "(A)","time,day"
        PRINT *,time,day,PZ(bPZ:ePZ)
        STOP
     END IF

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

     !      micro%sink=1.1574D-05
     !      CALL advection(grid,micro%sink,P%micro%old,dt,Gvector_ao%p%micro) !sinking phytoplankton

     !      Gvector_ao%d(D_bins)%bin = 0.

     DO xx = 1,D_bins-1
        CALL advection(grid,Detritus(2)%v,Detritus(2)%D%old,dt,Gvector_ao%d(2)%bin)
     END DO
     
     CALL reaction_p_sog (grid%M, PZ_bins, D_bins, PZ(1:PZ_bins%Quant*M), & !in
          PZ((PZ_bins%det-1)*M+1:M2), P%micro%old, P%nano%old, N%O%old, &   !in
          N%H%old, Detritus, &                                              !in
          Gvector_ro)                                       ! out
     Gvector_ro%Sil = 0 ! for now

     CALL matrix_A (Amatrix%bio, Bmatrix%bio)   !define Amatrix%A,%B,%C
     CALL matrix_A (Amatrix%no,Bmatrix%no)

     IF (time_step == 1) THEN
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

        DO xx = 1,D_bins
           Gvector_o%d(xx)%bin = Gvector%d(xx)%bin
        END DO

     END IF ! time_step == 1

     DO xx2 = 1,2
        IF (xx2 == 1) THEN
           CALL P_H(grid,Hvector%p%micro,Gvector%p%micro,Gvector_o%p%micro,Gvector_o_o%p%micro, &
                null_vector, null_vector,Gvector_ao%p%micro,Gvector_ao_o%p%micro, &
                Bmatrix_o%bio,Bmatrix_o_o%bio,P%micro)
           CALL N_H(grid,Hvector%p%nano,Gvector%p%nano,Gvector_o%p%nano,Gvector_o_o%p%nano, &
                null_vector, null_vector,Bmatrix_o%bio,Bmatrix_o_o%bio,P%nano) !V.flagella

           DO xx = 1,D_bins-1
              CALL P_H(grid,Hvector%d(xx)%bin,Gvector%d(xx)%bin,Gvector_o%d(xx)%bin, &
                   Gvector_o_o%d(xx)%bin, null_vector, null_vector, &
                   Gvector_ao%d(xx)%bin,Gvector_ao_o%d(xx)%bin,Bmatrix_o%bio,Bmatrix_o_o%bio,&
                   Detritus(xx)%D)
           END DO

           CALL N_H(grid,Hvector%d(D_bins)%bin,Gvector%d(D_bins)%bin,Gvector_o%d(D_bins)%bin, &
                Gvector_o_o%d(D_bins)%bin, null_vector, null_vector, &
                Bmatrix%null,Bmatrix%null,Detritus(D_bins)%D)

           CALL N_H(grid,Hvector%n%o,Gvector%n%o,Gvector_o%n%o,Gvector_o_o%n%o, &
                null_vector, null_vector,Bmatrix_o%no,Bmatrix_o_o%no,N%O)
           CALL N_H(grid,Hvector%n%h,Gvector%n%h,Gvector_o%n%h,Gvector_o_o%n%h, &
                null_vector, null_vector,Bmatrix_o%bio,Bmatrix_o_o%bio,N%H)

        ELSE IF (xx2 == 2) THEN
           CALL P_H(grid,Hvector%p%micro,Gvector%p%micro,Gvector_o%p%micro,Gvector_o_o%p%micro, &
                Gvector_ro%p%micro,Gvector_ro_o%p%micro,Gvector_ao%p%micro,Gvector_ao_o%p%micro, &
                Bmatrix_o%bio,Bmatrix_o_o%bio,P%micro)
           CALL N_H(grid,Hvector%p%nano,Gvector%p%nano,Gvector_o%p%nano,Gvector_o_o%p%nano, &
                Gvector_ro%p%nano,Gvector_ro_o%p%nano,Bmatrix_o%bio,Bmatrix_o_o%bio,P%nano) !V.flagella


           DO xx = 1,D_bins-1
              CALL P_H(grid,Hvector%d(xx)%bin,Gvector%d(xx)%bin,Gvector_o%d(xx)%bin, &
                   Gvector_o_o%d(xx)%bin,Gvector_ro%d(xx)%bin,Gvector_ro_o%d(xx)%bin, &
                   Gvector_ao%d(xx)%bin,Gvector_ao_o%d(xx)%bin,Bmatrix_o%bio,Bmatrix_o_o%bio,&
                   Detritus(xx)%D)
           END DO


           CALL N_H(grid,Hvector%d(D_bins)%bin,Gvector%d(D_bins)%bin,Gvector_o%d(D_bins)%bin, &
                Gvector_o_o%d(D_bins)%bin,Gvector_ro%d(D_bins)%bin,Gvector_ro_o%d(D_bins)%bin, &
                Bmatrix%null,Bmatrix%null,Detritus(D_bins)%D)

           CALL N_H(grid,Hvector%n%o,Gvector%n%o,Gvector_o%n%o,Gvector_o_o%n%o, &
                Gvector_ro%n%o,Gvector_ro_o%n%o,Bmatrix_o%no,Bmatrix_o_o%no,N%O)
           CALL N_H(grid,Hvector%n%h,Gvector%n%h,Gvector_o%n%h,Gvector_o_o%n%h, &
                Gvector_ro%n%h,Gvector_ro_o%n%h,Bmatrix_o%bio,Bmatrix_o_o%bio,N%H)


        END IF
        CALL TRIDAG(Amatrix%bio%A,Amatrix%bio%B,Amatrix%bio%C,Hvector%p%micro,P1_p,M)
        CALL TRIDAG(Amatrix%bio%A,Amatrix%bio%B,Amatrix%bio%C,Hvector%p%nano,Pnano1_p,M)
        CALL TRIDAG(Amatrix%no%A,Amatrix%no%B,Amatrix%no%C,Hvector%n%o,NO1_p,M) 
        CALL TRIDAG(Amatrix%bio%A,Amatrix%bio%B,Amatrix%bio%C,Hvector%n%h,NH1_p,M)



        DO xx = 1,D_bins-1
           CALL TRIDAG(Amatrix%bio%A,Amatrix%bio%B,Amatrix%bio%C,Hvector%d(xx)%bin,Detritus1_p(xx,:),M)
        END DO
        CALL TRIDAG(Amatrix%null%A,Amatrix%null%B,Amatrix%null%A,Hvector%d(D_bins)%bin,Detritus1_p(D_bins,:),M)

        IF (xx2 == 1) THEN
           CALL find_PON
        ELSE
           CALL find_new
        END IF
     END DO

     !-----END BIOLOGY------------------------------------------------

     !--------bottom boundaries--------------------------

     N%O%new(M+1) = ctd_bottom(day_met-281)%No
     P%micro%new(M+1) = ctd_bottom(day_met-281)%P
     P%nano%new(M+1) = ctd_bottom(day_met-281)%P !V.flagella ???
     S%new(M+1) = ctd_bottom(day_met-281)%sal
     T%new(M+1) = ctd_bottom(day_met-281)%temp+273.15
     N%H%new(M+1) = N%H%new(M)
     Detritus(1)%D%new(M+1)=Detritus(1)%D%new(M)
     Detritus(2)%D%new(M+1)=Detritus(2)%D%new(M)
     Detritus(3)%D%new(M+1)=Detritus(3)%D%new(M)


     dummy_time = dummy_time +dt

     ! Increment time, calendar date and clock time, unless this is
     ! the last time through the loop
     if(time_step < steps) then
        call new_year(day_time, day, year, time, dt, month_o)
     endif
  end do  !--------- End of time loop ----------

  ! Write profiles
  ! Calculate the month number and month day for profile headers
  profileDatetime%yr = year
  profileDatetime%yr_day = day
  profileDatetime%day_sec = day_time
  call calendar_date(profileDatetime)
  call clock_time(profileDatetime)
  ! Get the profile results file name, and open it
  ! *** Reading the file name should be done much earlier
  str = getpars("profile_out", 1)
  open(unit=profiles, file=str)
  ! Write the profile results file header
  ! Avoid a pgf90 idiocyncracy by getting datetimes formatted into
  ! string here rather than in the write statement
  str_runDatetime = datetime_str(runDatetime)
  str_proDatetime = datetime_str(profileDatetime)
  write(profiles, 200) trim(codeId), str_runDatetime, time, &
       str_proDatetime
200 format("! Profiles of Temperature, Salinity, Density, ",         &
       "Phytoplankton, Nitrate, Ammonium,"/,                         &
       "! and Detritus (remineralized, sinking, and mortality)"/,    &
       "*FromCode: ", a/,                                            &
       "*RunDateTime: ", a/,                                         &
       "*FieldNames: depth, temperature, salinity, sigma-t, ",       &
       "phytoplankton, nitrate, ammonium, remineralized detritus, ", &
       "sinking detritus, mortality detritus, ",                     &
       "total momentum eddy diffusivity, ",                          &
       "total temperature eddy diffusivity, ",                       &
       "total salinity eddy diffusivity, ",                          &
       "photosynthetic available radiation, ",                       &
       "u velocity, v velocity"/,                                    &
       "*FieldUnits: m, deg C, None, None, uM N, uM N, uM N, ",      &
       "uM N, uM N, uM N, m^2/s, m^2/s, m^2/s, W/m^2, m/s, m/s"/,    &
       "*ProfileTime: ", f9.0/                                       &
       "*ProfileDateTime: ", a/,                                     &
       "*EndOfHeader")
  ! Write the profile values at the surface, and at all grid points
  ! *** Change this index to i when we know common block effects have been
  ! *** cleaned up and i is safe to use as a local index
  do i_pro = 0, M
     sigma_t = density%new(i_pro) - 1000.
     write(profiles, 201) grid%d_g(i_pro), KtoC(T%new(i_pro)),          &
          S%new(i_pro), sigma_t, P%micro%new(i_pro), N%O%new(i_pro),    & 
          N%H%new(i_pro), Detritus(1)%D%new(i_pro),                     &
          Detritus(2)%D%new(i_pro), Detritus(3)%D%new(i_pro),           &
          K%u%all(i_pro), K%t%all(i_pro), K%s%all(i_pro), I_par(i_pro), &
          U%new(i_pro), V%new(i_pro)
  end do
  ! Write the values at the bottom grid boundary.  Some quantities are
  ! not defined there, so use their values at the Mth grid point.
  sigma_t = density%new(M+1) - 1000.
  write(profiles, 201) grid%d_g(M+1), KtoC(T%new(M+1)),        &
          S%new(M+1), sigma_t, P%micro%new(M+1), N%O%new(M+1), &
          N%H%new(M+1), Detritus(1)%D%new(M+1),                &
          Detritus(2)%D%new(M+1), Detritus(3)%D%new(M+1),      &
          K%u%all(M), K%t%all(M), K%s%all(M), I_par(M),        &
          U%new(M+1), V%new(M+1)
201 format(f7.3, 15(2x, f8.4))
  close(unit=profiles)

  ! Deallocate memory
  call dalloc_water_props(Cp)
  call dalloc_grid

end program SOG
