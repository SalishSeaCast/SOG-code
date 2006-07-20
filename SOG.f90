! $Id$
! $Source$

program SOG      
  ! Coupled physical and biological model of the Strait of Georgia

  use declarations
  use surface_forcing
  use initial_sog
  use pdf
  use IMEX_constants  

  implicit none

  external derivs_sog, rkqs

  ! Local variables:
  integer :: bin_tot2, smooth_i, icheck
  integer :: isusan, ecmapp, day_met
  double precision:: cz, unow, vnow, upwell, uplume, factor, Sa
  DOUBLE PRECISION, DIMENSION(0:80+1) :: dS
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
  ! Variables for pressure gradient calculations
  double precision :: sumu, sumv, uprev, vprev, delu, delv
  double precision :: sumpbx, sumpby, tol
  integer :: ii       ! loop index
  ! Physical model domain parameter (should go elsewhere)
  double precision :: oLx, oLy, gorLx, gorLy

  ! Open a bunch of output files
  CALL write_open


  ! infile contains the names of several other files to read
  ! It is opened on unit 5 to allow those names to be read by getpars()
  ! *** but that could be changed to use redirection of stdin
  open(5, file="infile")

  ! Get the name of the main run parameters file, open it, and read them
  str = getpars("inputfile", 1)
  open(10, file=str, status="OLD", action="READ")
  read(10, *) M, D, lambda , t_o, t_f, dt, day_o, year_o, month_o, &
       windyear, stormday, cruise_id

  steps = 1 + int((t_f - t_o) / dt) !INT rounds down

  ! These constants should be set as parameter somewhere else, or read
  ! from the main run parameters file
  Csources = 1
  C_types = 1
  D_bins = 3
  grid%M = M
  grid%D = D
  wind_n = 46056
  stable = 1
  Large_data_size = 63617

  !print*,day,ecmapp,year_o,'sog 1'

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
  CALL read_sog
  CALL coefficients(alph, beta, Cp, dens, cloud,p_Knox)
  CALL allocate3(alloc_stat)
  DO xx = 1,12
     IF (alloc_stat(xx) /= 0) THEN
        PRINT "(A)","ALLOCATION failed.  KPP.f  xx:"
        PRINT *,xx,alloc_stat(xx)
        STOP
     END IF
  END DO

  CALL initialize ! initializes everything (biology too)
  CALL define_grid(grid, D, lambda) ! sets up the grid
  CALL initial_mean(U,V,T,S,P,N,Detritus,h%new,ut,vt,pbx,pby,grid,D_bins,cruise_id)

  M2 = 6*M   !size of PZ in biology
  max_length = M2   !      max_length = MAXVAL(Cevent%length) Amatrix...

  CALL allocate4



  !---End Biology from KPP--------------------------------------

  IF(h%new < grid%d_g(1))THEN 
     h%new = grid%d_g(1)
  END IF

  !    initialize h_m
  h_m%new = 10.


  DO time_step = 1, steps    !main loop of program

     CALL define_sog(time_step) !!store previous time_steps in old (n) and old_old (n-1)!!     

     niter = 30
     DO count = 1, niter !20  !cutoff at 20 for iterations

        CALL find_jmax_g(h, grid) !d_g(j_max_g) > h ie, hh%g is grid point just below h
        CALL find_jmax_i(h, grid) !d_i(j_max_i) > h ie h%i is just above h

        CALL allocate2(alloc_stat) ! allocate shape and K's
        DO xx = 16,17
           IF (alloc_stat(xx) /= 0) THEN
              PRINT "(A)","ALLOCATION failed.  KPP.f  xx:"
              PRINT *,xx
              STOP
           END IF
        END DO


        IF (time_step == 1 .AND. count == 1) THEN
           CALL alpha_sub(T%new, S%new, alph, grid) ! find the actual alpha values over the grid 
           CALL alpha_sub(T%new, S%new, beta, grid)
           CALL Cp_sub(T%new, S%new, Cp, grid)   
           CALL density_sub(T, S, density%new, M,rho_fresh_o) ! density with depth and density of fresh water

           CALL div_i_param(grid, alph) ! alph%idiv = d alpha /dz
           CALL div_i_param(grid, beta) ! beta@idiv = d beta /dz
        END IF

        IF (count == 1) THEN   !only starts with zooplankton!! <- no idea what this means 
           !        

           ! fixed cloud_type
           ! *** Implicit type conversion problem
           j = day_time / 3600.0 + 1

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

           !PRINT*,cf(day_met,j),day_met,day
           !pause

           CALL irradiance_sog(cloud,cf(day_met,j),day_time,day,I,I_par,&
                grid, &
                water, jmax_i, Q_sol,euph,Qinter,h,P)

           !  CALL define_adv(grid,P_q,P_f,P_q_fraction(month_o),A_q,W_q,A_f,W_f,h) ! fudge advection

        ENDIF

        CALL buoyancy(alph, T%new, S%new, &
             grid, h, B%new, I, Br, density%new, Cp, beta, Q_n) 


        ! Br radiative contribution to the surface buoyancy forcing
        ! B%new buoyancy profile B(d)
        ! Q_n nonturbulent heat flux profile
        CALL div_grid(grid, density)
        ! take density from grid points to interface or vice-versa

        ! KPP uses smooth_x below but this is just set equal to 0 in input file 
        !print*,'here'
        call find_wind(year, day, time, ecmapp, wind_n, unow, vnow, wind)
        !print*,'here'
        !            write (*,*) year,day,time,ecmapp,wind_n,unow,vnow,' in SOG'

        ! *** Implicit type conversion problem
        ! *** Identical statement to this above - why???
        j = day_time / 3600.0 + 1
        !            write (*,*) day, j, atemp(day,j), 'air temperature'



        !PRINT*,day,day_met,cf(day_met,1)/10.,'day,day_met,cf(day_met,j)/10.'


        ! Interpolate river flows for the second we're at
        ! *** Implicit type conversion problem here!!!
        Qinter = (day_time * Qriver(day_met) &
             + (86400. - day_time) * Qriver(day_met-1)) / 86400.
        Einter = (day_time * Eriver(day_met) &
             + (86400. - day_time) * Eriver(day_met-1)) / 86400.

        !PRINT*,'humid',humid(day_met,j),j,day_met,cf(day_met,j),atemp(day_met,j)
        !pause

        CALL surface_flux_sog(grid%M,density%new,w,wt_r,S%new(0),S%new(h%i),S%new(M),T%new(0),j_gamma, &
             I, Q_t(0), alph%i(0), Cp%i(0), beta%i(0),unow, vnow, cf(day_met,j)/10., atemp(day_met,j), humid(day_met,j), &
             Qinter,stress, rho_fresh_o,day,dt/grid%i_space(1),h,upwell,Einter,u%new(1)) 

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

        CALL define_flux        !defines w%x, K%x%all, K%x%old, Bf%b, and F_n 



        ! calculate matrix B(changed to Amatrix later) and RHS called Gvector (D9) (D10)
        CALL diffusion(grid,Bmatrix%t,Gvector%t,K%t%all,gamma%t,w%t(0),-Q_n,dt,T%new(grid%M+1))
        CALL diffusion(grid,Bmatrix%s,Gvector%s,K%s%all,gamma%s,w%s(0),F_n,dt,S%new(grid%M+1))
        ! start with basic and all Coriolis later, extra scalar terms work out to 0
        CALL diffusion(grid,Bmatrix%u,Gvector%u,K%u%all,gamma%m,w%u(0),gamma%m,dt,U%new(grid%M+1))
        CALL diffusion(grid,Bmatrix%u,Gvector%v,K%u%all,gamma%m,w%v(0),gamma%m,dt,V%new(grid%M+1))           

        !      CALL horizontal_adv(grid,Gvector%t,P_q,density%new(1:M),Cp%g(1:M),dt)      
        !      CALL horizontal_adv(grid,Gvector%s,-P_f,density%new(1:M),1./S%new(1:M),dt)  !add salt to bottom 
        !!commented out by KC, july 20 because since P_q_fraction is 0, and define_adv is zero, P_q and P_f are zero so nothing happens in this function

        CALL find_upwell(grid,wupwell,upwell,S%new)
        !      CALL find_upwell(grid,S%new,Qriver(day),h_m,wupwell,upwell,factor,dS,U%new)

        CALL define_adv_bio(grid,S%new,Gvector%s,dt,P_sa,wupwell,grid%i_space(1))  !upwell salt similar to NO      
        CALL define_adv_bio(grid,T%new,Gvector%t,dt,P_ta,wupwell,grid%i_space(1))  !upwell temp too
        CALL define_adv_bio(grid,U%new,Gvector%u,dt,P_u,wupwell,grid%i_space(1))  !upwell u similar to NO
        CALL define_adv_bio(grid,V%new,Gvector%v,dt,P_v,wupwell,grid%i_space(1))  !upwell v similar to NO

        CALL Coriolis_and_pg(grid,V%new,pbx,Gvector_c%u,dt)
        CALL Coriolis_and_pg(grid,-U%new,pby,Gvector_c%v,dt)      

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
        CALL matrix_A(Amatrix%u,Bmatrix%u,time_step)   ! changed from Bmatrix to Amatrix
        CALL matrix_A(Amatrix%t,Bmatrix%t,time_step)   ! now have (D9)
        CALL matrix_A(Amatrix%s,Bmatrix%s,time_step)

        ! add Xt to H vector (D7)
        !PRINT*,'Hvec',Hvector%s(1)

        CALL scalar_H(grid,Hvector%s,Gvector%s,Gvector_o%s,Gvector_o_o%s,Bmatrix_o%s,Bmatrix_o_o%s,&
             time_step,S)
        !PRINT*,'Gvec',Gvector%s(1),Hvector%s(1)
        CALL scalar_H(grid,Hvector%t,Gvector%t,Gvector_o%t,Gvector_o_o%t,Bmatrix_o%t,Bmatrix_o_o%t,&
             time_step,T)

        ! add in Coriolis term (Gvector_c) and previous value to H vector (D7)
        CALL U_H(grid,Hvector%u,Gvector%u,Gvector_o%u,Gvector_o_o%u,Gvector_c%u,Gvector_co%u,&
             Gvector_co_o%u,Bmatrix_o%u,Bmatrix_o_o%u,time_step,U)
        CALL U_H(grid,Hvector%v,Gvector%v,Gvector_o%v,Gvector_o_o%v,Gvector_c%v,Gvector_co%v,&
             Gvector_co_o%v,Bmatrix_o%u,Bmatrix_o_o%u,time_step,V)
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

        !         U%new(M+1) = U%old(M+1) 
        !         V%new(M+1) = V%old(M+1) 
        U%new(M+1) = U%new(M) 
        V%new(M+1) = V%new(M) 
        T%new(M+1) = T%old(M+1)
        S%new(M+1) = S%old(M+1)

        CALL alpha_sub(T%new, S%new, alph, grid) ! find actual alpha values over grid
        CALL alpha_sub(T%new, S%new, beta, grid)
        CALL Cp_sub(T%new, S%new, Cp, grid)
        !PRINT*,'dens',dens_i(0),S%new(0),T%new(0)
        CALL density_sub(T, S, dens_i, M,rho_fresh_o) ! update density
        !PRINT*,'dens',dens_i(0)
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

        e_flux = MINVAL(w%b)   !entrainment flux

        ! beta_t is the ratio of the entrainment flux to the surface buoyancy flux

        beta_t = 0.
        V_t_square = 0.

        IF (w%b(0) /= 0.0 .AND. e_flux/w%b(0) < 0.) THEN
           beta_t = e_flux/w%b(0) !not -0.2 for beginning conv
        END IF

        IF (beta_t < 0.) THEN  !omega, vertical velocity scale
           CALL def_v_t_sog(grid, h, N_2_g, omega%s%value, V_t_square, beta_t, L_star) !test conv  
           ! V_t_square, the turbulent velocity shear, (23)
        END IF

        CALL define_Ri_b_sog(grid, h, surface_h, B, U, V, density, Ri_b, V_t_square, N_2_g)
        ! (21)
        CALL ML_height_sog(grid, Ri_b, h_i, jmaxg) !test conv
        ! (21) tested for minimum value of d at which Ri_b = Ri_c
        !PRINT*,'hi',h_i

        ! detrain for plume effect
        !         h_i = h_i - 0.1*S%new(jmaxg+2)*50/(S%new(jmaxg+2)-S%new(1))/ &
        !                  (h_i*100e3)*dt

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

        !PRINT*,'hi,hnew',h_i,h%new

        h_i = (h_i + h%new)/2.0 !use the average value!!!!!!!!##

        !PRINT*,'average, new h_i',h_i
        !pause
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
                 ! write (*,*) 'L*,u*, h_i', L_star, u_star, h_i
              END IF
           END IF
        END IF

        IF (h_i < grid%d_g(1)) THEN  !minimum mixing
           h_i = grid%d_g(1)
           write (*,*) 'minimum mixing', h_i
        END IF

        CALL find_jmax_g(h,grid) !***! can't be less than the grid

        del = ABS(h_i - h%new)/grid%i_space(jmaxg)  !##
        ! write (*,*), 'count,h_i,h%new,del',count,h_i,h%new,del 
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

        !         write (*,*) delu,delv,U%new(1),V%new(1)

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
           !            if (yy.lt.h%g) then
           !               vt%new(yy) = vt%new(yy) + dt*Qriver(day)*oLx*oLx/(2.*h%new)*S%new(M)/(S%new(M)-S%new(0))
           !            else
           !               vt%new(yy) = vt%new(yy) - dt*Qriver(day)*oLx*oLx/(2.*(80.-h%new))*S%new(M)/(S%new(M)-S%new(0))
           !            endif
           !            write (*,*) yy,vt%new(yy),h%g,h%new
        enddo

        dzx(1) = ut%new(1)+1
        dzy(1) = vt%new(1)+1
        !            write (*,*) dzx(1),dzy(1),' dzx'
        do yy=2,M
           dzx(yy) = dzx(yy-1) + (ut%new(yy)+1)
           dzy(yy) = dzy(yy-1) + (vt%new(yy)+1)
           !            write (*,*) yy,dzx(yy),dzy(yy),' dzx'
        enddo

        sumpbx = 0.
        tol=1e-6
        cz = 0.
        ii=1
        do yy=1,M
           if (yy == 1) then
              pbx(yy) = -density%new(yy)
              !               write (*,*) ii,yy,pbx(yy),cz,density%new(yy),'  **1'
           else
              pbx(yy) = pbx(yy-1)-density%new(yy)
              !               write (*,*) ii,yy,pbx(yy),cz,density%new(yy),'  **1'
           endif
           do while ((dzx(ii)-yy) < -tol)
              pbx(yy) = pbx(yy) + density%new(ii)*(dzx(ii)-cz)
              cz = dzx(ii)
              !               write (*,*) ii,yy,pbx(yy),cz,density%new(ii),'  **2'
              ii = ii + 1
           enddo
           pbx(yy) = pbx(yy) + density%new(ii)*(yy-cz)
           sumpbx = sumpbx + pbx(yy)
           !            write (*,*) ii,yy,pbx(yy),yy-cz,density%new(ii),'  **3'
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
           do while ((dzy(ii)- yy) <-tol)
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
                + stress%u%new / (1025 * M * grid%i_space(yy))
           pby(yy) = (pby(yy) - sumpby) * gorLy * grid%i_space(yy) &
                + stress%v%new / (1025 * M * grid%i_space(yy))     &
                + dS(yy) * uplume * factor * K%u%all(h_m%g) / h_m%new * 20.
           !            write (*,*) 'tadkfj', yy,h%i,uplume,pbx(yy)
           !            if (yy.le.h%i+1) then
           !               pbx(yy) = pbx(yy) + uplume/7200.
           !            else
           !               pbx(yy) = pbx(yy) - uplume*(h%i+1)/(M-h%i-1)/7200.
           !            endif
           !            write (*,*) 'tadkfj', yy,h%i,uplume,pbx(yy)
           !            write (*,173) v%new(yy),pby(yy),sumpby*gorLy*grid%i_space(yy),vt%new(yy),-stress%v%new/(1025*M*grid%i_space(yy))
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


        IF (count >= 2 .AND. del < del_o) THEN
           IF (count == 2) THEN
              count_two = count_two + 1         
           END IF
           count_tot = count_tot + 1
           !            write (*,*) 'count ',count,' Next step'
           CALL write_physical_sog(unow,vnow,euph)
           EXIT
        ELSE IF (count == niter .AND. del >= del_o) THEN
           !            write (*,*) 'count ',count,' Next step'
           count_no =  count_no + 1 
           CALL write_physical_sog (unow,vnow,euph)
        ELSE IF (count == 1) THEN
           !            CALL write_physical_sog
        END IF
        !         call write_physical_sog

     ENDDO

     !------BIOLOGICAL MODEL--------------------------------------------

     !      PRINT "(A)","KPP: start bio model ==> day,time,year"
     !      PRINT *,day,time,year

     CALL define_PZ(icheck)

     IF (MINVAL(PZ) < 0.) THEN
        PRINT "(A)","PZ < 0. see KPP.f90"
        PRINT "(A)","time,day"
        PRINT *,time,day
        STOP
     END IF

     next_time = time+dt

     CALL odeint(PZ,M2,time,next_time,precision,step_guess,step_min,N_ok,N_bad,derivs_sog,rkqs,icheck,T%new)

     IF (MINVAL(PZ(2*M+1:3*M)) < 0.) THEN
        DO xx = 2*M+1,3*M
           IF (PZ(xx) < 0.) THEN
              PRINT "(A)","N%H%new(xx-4*M) < 0."
              PRINT *,PZ(xx)
              PRINT "(A)","set equal to zero"
              PZ(xx) = 0.
           END IF
        END DO
     END IF

     IF (MINVAL(PZ(1:M)) < 0.) THEN
        PRINT "(A)","PZ < 0. After odeint.f see KPP.f90"
        PRINT "(A)","time,day"
        PRINT *,time,day,PZ(1:M)
        STOP
     END IF

     CALL diffusion(grid,Bmatrix%bio,Gvector%p%micro,K%s%all,gamma%m,pflux_o,gamma%m,dt,&
          P%micro%new(grid%M+1))

     Gvector%n%h = Gvector%p%micro
     Gvector%d(D_bins)%bin = 0.
     DO xx = 1, D_bins-1
        Gvector%d(xx)%bin = Gvector%p%micro
     END DO

     CALL diffusion(grid,Bmatrix%no,Gvector%n%o,K%s%all,gamma%m,pflux_o,gamma%m,dt,&
          N%O%new(grid%M+1))

     CALL define_adv_bio(grid,P%micro%new,Gvector%p%micro,dt,P_na,wupwell,grid%i_space(1))
     CALL define_adv_bio(grid,N%O%new,Gvector%n%o,dt,P_no,wupwell,grid%i_space(1))
     CALL define_adv_bio(grid,N%H%new,Gvector%n%h,dt,P_nh,wupwell,grid%i_space(1))


     DO xx = 1,D_bins-1
        IF (xx == 1) THEN
           CALL define_adv_bio(grid,Detritus(xx)%D%new,Gvector%d(xx)%bin,dt,P_d1,wupwell,grid%i_space(1))
        ELSE
           CALL define_adv_bio(grid,Detritus(xx)%D%new,Gvector%d(xx)%bin,dt,P_d2,wupwell,grid%i_space(1))
        END IF
     END DO

     !      micro%sink=1.1574D-05
     !      CALL advection(grid,micro%sink,P%micro%old,dt,Gvector_ao%p%micro) !sinking phytoplankton

     !      Gvector_ao%d(D_bins)%bin = 0.



     DO xx = 1,D_bins-1
        CALL advection(grid,Detritus(2)%v,Detritus(2)%D%old,dt,Gvector_ao%d(2)%bin)
     END DO

     CALL reaction_p_sog(grid,M2,D_bins,PZ(1:3*M),PZ(3*M+1:5*M),P%micro%old,N%O%old,N%H%old,Detritus,Gvector_ro)  !define Gevtor_ro
     bin_tot2 = 0.




     CALL matrix_A(Amatrix%bio,Bmatrix%bio,time_step)   !define Amatrix%A,%B,%C
     CALL matrix_A(Amatrix%no,Bmatrix%no,time_step)

     IF (time_step == 1) THEN
        Bmatrix%null%A = 0. !no diffusion
        Bmatrix%null%B = 0.
        Bmatrix%null%C = 0.
        Amatrix%null%A = 0. !(M)
        Amatrix%null%B = 1.
        Amatrix%null2%A = 0.!(bin)  
        Amatrix%null2%B = 1.
        Bmatrix_o%bio%A = Bmatrix%bio%A
        Bmatrix_o%bio%B = Bmatrix%bio%B
        Bmatrix_o%bio%C = Bmatrix%bio%C
        Bmatrix_o%no%A = Bmatrix%no%A
        Bmatrix_o%no%B = Bmatrix%no%B
        Bmatrix_o%no%C = Bmatrix%no%C
        Gvector_o%p%micro = Gvector%p%micro
        Gvector_o%n%o = Gvector%n%o
        Gvector_o%n%h = Gvector%n%h

        DO xx = 1,D_bins
           Gvector_o%d(xx)%bin = Gvector%d(xx)%bin
        END DO

     END IF


     !------------------- 


     DO xx2 = 1,2
        IF (xx2 == 1) THEN
           CALL P_H(grid,Hvector%p%micro,Gvector%p%micro,Gvector_o%p%micro,Gvector_o_o%p%micro, &
                null_vector, null_vector,Gvector_ao%p%micro,Gvector_ao_o%p%micro, &
                Bmatrix_o%bio,Bmatrix_o_o%bio,P%micro)

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


     !CALL write_biological_sog


     !-----END BIOLOGY------------------------------------------------

     !--------bottom boundaries--------------------------

     N%O%new(M+1) = ctd_bottom(day_met-281)%No
     P%micro%new(M+1) = ctd_bottom(day_met-281)%P
     S%new(M+1) = ctd_bottom(day_met-281)%sal
     T%new(M+1) = ctd_bottom(day_met-281)%temp+273.15
     N%H%new(M+1) = N%H%new(M)
     Detritus(1)%D%new(M+1)=Detritus(1)%D%new(M)
     Detritus(2)%D%new(M+1)=Detritus(2)%D%new(M)
     Detritus(3)%D%new(M+1)=Detritus(3)%D%new(M)


     dummy_time = dummy_time +dt
     CALL new_year(day_time,day,year,time,dt,day_check,day_check2,month_o,data_point_papmd)


  END DO  !main loop, ie time loop!

  ! Write profiles
131  format (13(1x,f10.4))
  ! Quantities without initial condition value value at bottom of the
  ! model domain (M values in the profile)
  ! *** K, I_par
  open(5, file="infile")
  str = getpars("profile_out1", 1)
  open(4, file=str)
  do isusan = 0, M
     write(4, 131) K%u%all(isusan), K%s%all(isusan), I_par(isusan), &
          K%t%all(isusan)
  end do
  close(4)

  ! Quantities that include an intial condition value at the bottom of the
  ! model domain (M + 1 values in the profile)
  ! Water temperature, *** salinity, density, micro-phytos, nitrate, ammonia,
  ! and detritus
  open (5, file="infile")
  str = getpars("profile_out2", 1)
  open(6, file=str)
  do isusan = 0, M + 1
     density%new(isusan) = density%new(isusan) - 1000
     write(6, 131) T%new(isusan), s%NEW(isusan), density%new(isusan), &
          P%micro%new(isusan), N%O%new(isusan), N%H%new(isusan),      &
          Detritus(1)%D%new(isusan), Detritus(2)%D%new(isusan),       &
          Detritus(3)%D%new(isusan)
     ! *** Removed f_ratio, ratio of nitrate uptake to total nitrogen uptake
     ! *** from above write, because it is dimensioned (1:M) and causes an
     ! *** array bound error in this loop
  end do
  close(6)

end program SOG
