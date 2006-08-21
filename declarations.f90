! $Id$
! $Source$

module declarations
! *** What's in here?

  use mean_param

  implicit NONE

  ! Run parameters read from file specified as 'inputfile' in infile
  double precision :: D         ! depth of model domain [m]
  integer :: M                  ! # of grid points in model domain
  double precision :: lambda    ! grid spacing parameter ???
  integer :: year_o             ! year in which run starts [yr]
  ! *** Why is month_o not declared here; rather it's local in SOG.f90
  integer :: day_o              ! year-day on which run starts [1-Jan = 1]
  double precision :: t_o       ! start of run time on day_o [s]
                                ! *** t_o should match CTD profile time of day
  double precision :: t_f       ! end of run time [s]
  double precision :: dt        ! time step [s]
  integer :: windyear           ! *** apparently not used
  integer :: stormday           ! *** apparently not used
  integer :: cruise_id          ! *** apparently not used

  TYPE(constant)::alph, beta, dens
  TYPE(gr_d)::grid
  TYPE(prop)::T, S, U, V, B, density, P_temp, PON, vt, ut
  TYPE(plankton):: P 
  TYPE(zplankton):: Z
  TYPE(nutrient):: N
  TYPE(Knu)::K 
  TYPE(flux)::w,Bf  
  TYPE(height)::h, surface_h, h_m
  TYPE(MST)::gamma,G_shape
  TYPE(MS)::phi,omega
  TYPE(MSTscalar)::a2,a3  !shape constants
  TYPE(UVST)::Hvector, Gvector, Gvector_o, Gvector_o_o, Gvector_c, Gvector_co, Gvector_co_o, &
       Gvector_ao, Gvector_ao_o, Gvector_ro , Gvector_ro_o          !**&
  TYPE(UVSTmatrix)::Amatrix, Bmatrix, Bmatrix_o, Bmatrix_o_o !**&
  TYPE(okta)::cloud
  TYPE(light)::water
  TYPE(windstress)::stress !wind stress (surface_flux.f90)
  TYPE(plankton2)::micro, nano, zmicro
  TYPE(losses)::waste
  TYPE(wind_ecmwf), DIMENSION(:), ALLOCATABLE::wind  
  TYPE(insol_daily), DIMENSION(:), ALLOCATABLE::insol   !2922 data points, about 3 years worth  some of 83, 84, 85, 86
  TYPE(Large1996_data), DIMENSION(14)::QN_1996, FN_1996
  TYPE(papmd), DIMENSION(:), ALLOCATABLE::Large_data  !Stores data from file PAPMD.60
  TYPE(entrain)::euph

  DOUBLE PRECISION, DIMENSION(365)::insol_avg, insol_actual_avg, insol_smooth, insol_actual_smooth
  !                                        avgdata_u, avgdata_v
  DOUBLE PRECISION, DIMENSION(0:366)::smooth_u, smooth_v
  REAL, DIMENSION(1919,24)::cf,atemp,humid
  REAL, DIMENSION(1919)::Qriver
  REAL, DIMENSION(1919)::Eriver
  DOUBLE PRECISION::begin_hour,end_hour  ! used in interpolate_dt
  DOUBLE PRECISION :: time, h_i, del, dummy_time, density_test, del_p 

!!!!For surface fluxes (surface_flux.f90): 
  DOUBLE PRECISION::U_ten,V_ten, & ! U and V velocities (m/s) at standard height (10m or 22m)
       T_atm,& !air temperature at standard height (17m) (K)
       Q_st,Q_atm, & !specific humidity of air in contact with salt water and at standard height (17m)
       Q_tot, F_tot, & !Total turbulent surface heat flux (W/m^2) and freshwater flux
       wt_r, &  !Radiative contribution to surface heat flux
       rho_fresh_o !surface density of pure water at SST 
  INTEGER::j_gamma !interface point corresponding to depth at which Radiation contributes to turbulent
  !surface heat flux
  INTEGER::max_length
  DOUBLE PRECISION::Q_sol  !See irradiance.f90
  DOUBLE PRECISION::Br
  CHARACTER::ignored_input
  INTEGER, DIMENSION(20)::alloc_stat  !***
  INTEGER :: steps, time_step, index,xx, xx2, year, month, jmax_i,    &
       count, yy, yy2, gg, jmaxg, cloud_type, day, stable,  &
       water_type, t_start, t_end, neg_count, current_day, day_check, &
       day_check2,                                         &
       migrate  ! migrate = 1 if timestep has copepod migration, else 0
  REAL::D_test
  DOUBLE PRECISION :: u_star, w_star,L_star, day_time, Io, h_Ekman  !depth
  !friction vel, convective vel scale
  !monin_obukov length and julian day
  DOUBLE PRECISION::beta_t ! entrainment coefficient under convection
  INTEGER, DIMENSION(:), ALLOCATABLE::Q_test
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::V_t_square, &
       Ri_b, Q_n, F_n, &
       N_2_i, N_2_g, N_2_dens_g
  double precision, dimension(:), allocatable:: pbx,pby,dzx,dzy
  DOUBLE PRECISION, DIMENSION(:),  ALLOCATABLE:: I, I_par, T_To, Q_t, dens_i       
  !I = intensity, Q_n is the
  !non_turbulent heat flux  (i.e. I)
  !F_n is freshwater flux (0 no ice)

  DOUBLE PRECISION, DIMENSION(1)::microQ1_p,nanoQ1_p
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::P1_p, Pnano1_p, Z1_p, &
       NO1_p, NH1_p,&
       U_p, V_p, S_p, T_p, density_i, &
       buoy_i, uv_square_i, TKE_rate       
  !Temporary vectors U_i , U_p...
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::Copepod1_p,Copepod_wt1_p
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::cgraze   !total grazing  of (prey+d_prey,M) by copepods
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE::Detritus1_p
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::m_TKE,b_TKE, TKE
  DOUBLE PRECISION::KE_o,P_m,P_b,PE_o,IKE,PE, KE_old,PE_old,Pm_old,Pb_old, Perr, &
       Perr_old, TKE_o
  !Surface Kinetic energy, Mechanical Production of TKE, Buoyant production of TKE,
  !Surface Potential energy
  DOUBLE PRECISION::pflux_o  !surface biological flux == 0

  !Variables for printing test functions !

  DOUBLE PRECISION :: avg_T, read_var 
  INTEGER :: counter,count_two,count_tot,count_no,count2_two,count2_tot,count2_no
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::ref_T, avg_12, tot_avg

  !Variable for odeint.f,  rkqs.f and derivs.f

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::PZ
  !     DOUBLE PRECISION::precision, step_min, step_guess, next_time  !see Read_data.f90 and input file 
  DOUBLE PRECISION::next_time  !see Read_data.f90 and input file 
  ! input/forcing.dat
  INTEGER::M2, M3, N_ok, N_bad      

  !Variable for wind and insol data
  DOUBLE PRECISION::last_wind_u, last_wind_v
  INTEGER::data_point_insol,data_point_ecmwf
  integer :: wind_n
  integer :: insol_n, smooth_x

  !current data number and Size of PAPMD.60 data array Large_data(:)

  INTEGER::data_point_papmd, Large_data_size  ! see initialize.f90 and read_data.f90 
  !     INTEGER::data_point_papmd  ! see initialize.f90 and read_data.f90 
  DOUBLE PRECISION::vapour_pressure, prain ! (atmosphere at 17m) and current value for &
  !precipitation from p_Knox (Large_param.f90) !
  DOUBLE PRECISION, DIMENSION(23)::p_Knox  !Piecewise linear function for precipitation
  !defined in coefficients.f90 !

  !Advection corrections to heat and salinity profiles

  DOUBLE PRECISION,DIMENSION(12)::P_q_fraction !defined in coefficients
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::P_q,P_f  !define in define_adv.f90
  !     DOUBLE PRECISION::A_q,W_q,A_f,W_f  !read in values in read_data2.f90

  !Upwelling corrections to biology and salinity

  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::P_no,P_nh,P_di,P_na,P_mi,P_d1,P_d2,P_sa,P_u,P_v,P_ta,wupwell

  !Copepod variables

  INTEGER::Csources, C_types,bin_tot
  TYPE(event), DIMENSION(:), ALLOCATABLE::Cevent    
  TYPE(copepod), DIMENSION(:), ALLOCATABLE::species 
  TYPE(Cdata), DIMENSION(:), ALLOCATABLE::Zoo  

  !Detritus variables

  INTEGER::D_bins  !number of detrital compartments  see input/biology.dat
  TYPE(snow), DIMENSION(:), ALLOCATABLE::Detritus  

  !!Variables for plotting:  see write_biological.f90
  INTEGER, DIMENSION(243)::o_cnt,hm_g,euph_g,species_cnt, species_cnt2
  INTEGER::nobins,bin_no
  DOUBLE PRECISION, DIMENSION(243)::NO_o,nano_ml,diatom_ml,zmicro_ml,copepod_ml,NO_ml,NH_ml,PN_ml,DN_ml,&
       copepod_150,NH_150,NO_150,PN_150,DN_150,PN_100,DN_100,nano_e,diatom_e,NO_e,NH_e,PN_e,DN_e, &
       NOup_e,NHup_e,fratio_e, NPP_e,NPP_ml,NPP_80,ngrow_o,ngrow_e,ngrow_ml,dgrow_o,dgrow_e,dgrow_ml,&
       species_b,species_mwt,species_mwt_o,zgrazen_ml,cgrazed_ml,cgrazez_ml,cgrazed2_ml,SPN_ml,&
       T2nano_ml,Thalfnano_ml,T2diatom_ml,Thalfdiatom_ml,T2zmicro_ml,&
       Thalfzmicro_ml,PONflux200,PONfluxml,PONflux100,feacalml,ureaf,ureac,migrateflux,NOflux200,&
       NOfluxml,&
       NOflux100,NOup_ml,NHup_ml,fratio_ml,hm_new,euph_new,species_avgwt,D3_ml,D3_200,NHflux200,&
       NHfluxml,NHflux100,PN_200,NH_200,NO_200
  DOUBLE PRECISION, DIMENSION(0:50)::bin_logwt,bin_cnt     
  DOUBLE PRECISION, DIMENSION(1:27,0:50)::bin_cnt_year
  INTEGER, DIMENSION(:),ALLOCATABLE::daybins
  TYPE(write_bio)::mixlay, depth150, ezone, depth80, total, depth100,doubleT,halfT,depth200
  INTEGER::g_150, g_80, g_100 !grid spacing near 150m defined in initialize
  INTEGER,DIMENSION(27)::j_day !1,15,29,...,365  every 14 days
  DOUBLE PRECISION,DIMENSION(365)::c_biomass
  DOUBLE PRECISION::c_biomass_150
  INTEGER,DIMENSION(365)::c_cnt
  !   INTEGER, DIMENSION(365)::o_cnt !12*5
  !   DOUBLE PRECISION, DIMENSION(365)::nano_o, diatom_o, zmicro_o,copepod_o,don_o,pon_o,&
  !     out_o,NH_o,stage1_f,stage2_f,stage3_f,stage4_f,stage5_f,stage1_n,stage2_n,stage3_n,&
  !     stage4_n,stage5_n,stage6_n,nano_avg,diatom_avg,zmicro_avg,copepod_avg,don_avg,pon_avg,&
  !     out_avg,NO_avg,NH_avg,Ntot_avg,n_loss,nano_gML,micro_gML,nano_g50,micro_g50,nano_go,micro_go,&
  !     nano_zML,micro_zML,nano_zo,micro_zo,NPPd,NPPn,NPPd_o,NPPn_o
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::nano_mar,nano_jun,nano_sep,nano_dec,diatom_mar,&
       diatom_jun,diatom_sep,diatom_dec,zmicro_mar,zmicro_sep,zmicro_jun,zmicro_dec,copepod_mar,&
       copepod_sep,copepod_jun,copepod_dec,NO_mar,NO_jun,NO_sep,NO_dec,NH_mar,NH_jun,NH_sep,NH_dec,&
       don_jun,don_mar,don_sep,don_dec,pon_jun,pon_mar,pon_sep,pon_dec,NPPn_mar,NPPn_jun,NPPn_dec,&
       NPPn_sep,NPPd_mar,NPPd_jun,NPPd_dec,NPPd_sep,ngrow_mar,ngrow_jun,ngrow_sep,ngrow_dec,dgrow_mar,&
       dgrow_jun,dgrow_sep,dgrow_dec,ngraz_mar,ngraz_jun,ngraz_sep,ngraz_dec,dgraz_mar,dgraz_jun,&
       dgraz_dec,dgraz_sep,NOup_mar,NOup_sep,NOup_jun,NOup_dec,NHup_mar,NHup_jun,NHup_sep,NHup_dec,&
       fratio_mar,fratio_jun,fratio_dec,fratio_sep
  DOUBLE PRECISION, DIMENSION(:,:),ALLOCATABLE::nano_pro,diatom_pro,zmicro_pro,&
       copepod_pro,NO_pro,NH_pro,fratio_pro,T_pro,S_pro,U_pro,V_pro,PON_pro  !60,M
  INTEGER, DIMENSION(84)::pro_cnt 
  INTEGER::cnt_mar,cnt_sep,cnt_jun,cnt_dec
  !  DOUBLE PRECISION, DIMENSION(365)::wt_stage1,wt_stage2,wt_stage3,wt_stage4,wt_stage5,molt_wt,avg_wt
  !  INTEGER, DIMENSION(365)::cnt_wt,cnt_avg_wt
  !  DOUBLE PRECISION, DIMENSION(365)::urea_cop,urea_fla,NO_flux,PO_flux, fratio
  DOUBLE PRECISION::urea_c, urea_f   !daily urea production rate 
  DOUBLE PRECISION::NO50_rate,PO50_rate,feacal50_rate,NO100_rate,PO100_rate,feacal100_rate,&
       new_NO,old_NO,new_PO,old_PO,feacalml_rate  
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::f_ratio
  !nitrate and Particulate organics integrated over upper 50m

  !!variables for plotting: see write_physical.f90
  INTEGER, DIMENSION(243)::p_cnt
  DOUBLE PRECISION, DIMENSION(243)::SST,SSS,hm_avg,Ipar_o,Uten_o,Vten_o,UVten_o,Qflux,Fflux,&
       stage1_no,stage2_no,stage3_no,stage4_no,stage5_no,out_no,KsML,IparML
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::T_mar,T_jun,T_sep,T_dec,S_mar,S_jun,S_sep,S_dec,&
       U_mar,U_dec,U_sep,U_jun,V_mar,V_jun,V_sep,V_dec,Ku_mar,Ku_jun,Ku_sep,Ku_dec,&
       Ks_mar,Ks_dec,Ks_sep,Ks_jun,Kt_mar,Kt_dec,Kt_sep,Kt_jun
  INTEGER::cntp_mar,cntp_sep,cntp_jun,cntp_dec

  !!variables for finding PON flux
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE::null_vector,null_vector2
  DOUBLE PRECISION::PONflux_200,PONflux_100,PONflux_ml,NOflux_200,NOflux_ml,NOflux_100,NHflux_200,&
       NHflux_ml,NHflux_100  !due to physical processes

!!!!!!!variables from KPP used in shift_vector.f90

  DOUBLE PRECISION::smooth_avg_u,smooth_avg_v,smooth_var_u,smooth_var_v,test_avg,test_var

  INTEGER, DIMENSION(366)::cloud_type_rand
  INTEGER :: cloud_day, ironday

  !TYPE(bottom_fit), DIMENSION(48):: ctd_bottom
  TYPE(bottom_fit), DIMENSION(1659):: ctd_bottom

end module declarations
