! $Id$
! $Source$

SUBROUTINE initialize

      USE mean_param
      USE surface_forcing
      USE declarations

      IMPLICIT NONE

      INTEGER::k_k

      j_day(1) = 1
      DO k_k = 2,27
         j_day(k_k) = j_day(k_k-1) + 14
      END DO      

      time = t_o
      dummy_time = 0.
      day = day_o   !initialize julian day
     ! month = month_o   !KC june 28, 2004
       year = year_o

      day_time = time
      time_step = 1

      j_gamma = 0  !h_gamma = 0 or wt_r contribution to gamma%t vanishes

      neg_count = 0  ! counts number of times scheme biology becomes negative

!!!!Initial CN ratios for P%micro and P%nano

      micro%Q_cn = Q_min !6.67  
      nano%Q_cn = Q_min !6.67 

!write_biology.f90
      o_cnt = 0
      NO_o = 0.
      nano_ml = 0.
      diatom_ml = 0.
      zmicro_ml = 0.
      copepod_ml = 0.
      NO_ml = 0.
      NH_ml = 0.
      PN_ml = 0.
      DN_ml = 0.
      copepod_150 = 0.
      NH_150 = 0.
      PN_150 = 0.
      DN_150 = 0.
      PN_100 = 0.
      DN_100 = 0.
      nano_e = 0.
      diatom_e = 0.
      NO_e = 0.
      NO_150 = 0.
      NH_e = 0.
      NOup_e = 0.
      NHup_e = 0.
      NOup_ml = 0.
      NHup_ml = 0.
      fratio_e = 0.
      fratio_ml = 0.
      NPP_e = 0.
      NPP_ml = 0.
      NPP_80 = 0.
      hm_g = 0
      euph_g = 0
      hm_new = 0.
      euph_new = 0.
      ngrow_ml = 0.
      ngrow_o = 0.
      ngrow_e = 0.
      dgrow_ml = 0.
      dgrow_o = 0.
      dgrow_e = 0.
      zgrazen_ml = 0.
      cgrazed_ml = 0.
      cgrazez_ml = 0.
      cgrazed2_ml = 0.
      SPN_ml = 0.
      T2nano_ml = 0.
      T2diatom_ml = 0.
      Thalfdiatom_ml = 0.
      T2zmicro_ml = 0.
      Thalfzmicro_ml = 0.
      feacalml = 0.
      D3_ml = 0.
      D3_200 = 0.
      PN_200 = 0.
      NO_200 = 0.
      NH_200 = 0.
      ureac = 0.
      ureaf = 0.
      migrateflux = 0.
      species_b = 0.
      species_mwt = 0.
      species_mwt_o = 0.
      species_avgwt = 0.
      DO k_k = 1,27
         bin_cnt_year(k_k,:) = 0.
      END DO
!      nano_avg = 0.
!      diatom_avg = 0.
!      zmicro_avg = 0.
!      copepod_avg = 0.
!      don_avg = 0.
!      pon_avg = 0.
 !     out_avg = 0.
!      NO_avg = 0.
!       NH_avg = 0.     
!       nano_o = 0.
!       diatom_o = 0.
!       zmicro_o = 0.
!       copepod_o = 0.
!       don_o = 0.
!       pon_o = 0.
!       out_o = 0.
!       NO_o = 0.
!       NH_o = 0.
!       stage1_f = 0.
!       stage2_f = 0.
!       stage3_f = 0.
!       stage4_f = 0.
!       stage5_f = 0.
!       stage1_n = 0.
!       stage2_n = 0.
!       stage3_n = 0.
!       stage4_n = 0.
!       stage5_n = 0.      
!       stage6_n = 0.
!       Ntot_avg = 0.
      pro_cnt = 0
      nano_pro = 0.
      diatom_pro = 0.
      copepod_pro = 0.
      zmicro_pro = 0.
      NO_pro = 0.
      NH_pro = 0.
      fratio_pro = 0.
      PON_pro = 0.
      T_pro = 0.
      S_pro = 0.
      U_pro = 0.
      V_pro = 0.
      nano_mar = 0.
      nano_jun = 0.
      nano_sep = 0.
      nano_dec = 0.
      diatom_mar = 0.
      diatom_jun = 0.
      diatom_sep = 0.
      diatom_dec = 0.
      zmicro_mar = 0.
      zmicro_sep = 0.
      zmicro_jun = 0.
      zmicro_dec = 0.
      copepod_mar = 0.
      copepod_sep = 0.
      copepod_jun = 0.
      copepod_dec = 0.
      NO_mar = 0.
      NO_jun = 0.
      NO_sep = 0.
      NO_dec = 0.
      NH_mar = 0.
      NH_jun = 0.
      NH_sep = 0.
      NH_dec = 0.
      don_jun = 0.
      don_mar = 0.
      don_sep = 0.
      don_dec = 0.
      pon_jun = 0.
      pon_mar = 0.
      pon_sep = 0.
      pon_dec = 0.     
      NPPn_jun = 0.
      NPPn_mar = 0.
      NPPn_sep = 0.
      NPPn_dec = 0.
      NPPd_jun = 0.
      NPPd_mar = 0.
      NPPd_sep = 0.
      NPPd_dec = 0. 
      ngrow_jun = 0.
      ngrow_mar = 0.
      ngrow_sep = 0.
      ngrow_dec = 0.  
      dgrow_jun = 0.
      dgrow_mar = 0.
      dgrow_sep = 0.
      dgrow_dec = 0.  
      ngraz_jun = 0.
      ngraz_mar = 0.
      ngraz_sep = 0.
      ngraz_dec = 0.   
      dgraz_jun = 0.
      dgraz_mar = 0.
      dgraz_sep = 0.
      dgraz_dec = 0. 
      NOup_jun = 0.
      NOup_mar = 0.
      NOup_sep = 0.
      NOup_dec = 0.  
      NHup_jun = 0.
      NHup_mar = 0.
      NHup_sep = 0.
      NHup_dec = 0. 
      fratio_jun = 0.
      fratio_mar = 0.
      fratio_sep = 0.
      fratio_dec = 0.     
      cnt_mar = 0
      cnt_sep = 0
      cnt_jun = 0
      cnt_dec = 0
      species_cnt = 0
      species_cnt2 = 0
      !       wt_stage1 = 0.
      !       wt_stage2 = 0.
!       wt_stage3 = 0.
!       wt_stage4 = 0.
 !      wt_stage5 = 0.
!       cnt_wt = 0
!       cnt_avg_wt = 0
!       molt_wt = 0.
!       avg_wt = 0.
!       urea_cop = 0.
!       urea_fla = 0.
!       n_loss = 0.
!       nano_go = 0.
!       nano_gML = 0.
! !       nano_g50 = 0.
!       micro_go = 0.
! !       micro_gML = 0.
!       micro_g50 = 0.
!       nano_zML = 0.
!       nano_zo = 0.
!       micro_zML = 0.
!       micro_zo = 0.
!       NPPn = 0.
!       NPPn_o = 0.
!       NPPd = 0.
!       NPPd_o = 0.
!       total%n_loss = 0.
     !  NO_rate = 0.
     !  PO_rate = 0.
       new_NO = 0.
       new_PO = 0.
       old_NO = 0.
       old_PO = 0.
!       PO_flux = 0.
!       NO_flux = 0.
!       fratio = 0.
       f_ratio = 0.
       stage1_no = 0.
       stage2_no = 0.
       stage3_no = 0.
       stage4_no = 0.
       stage5_no = 0.
       out_no = 0.
       c_cnt = 0
       c_biomass = 0.
!write_physical.f90
      p_cnt = 0
      SST = 0.
      SSS = 0.
      Ipar_o = 0.
      hm_avg = 0.
      Uten_o = 0.
      Vten_o = 0.
      UVten_o = 0.
      Qflux = 0.
      Fflux = 0.
      KsML = 0.
      IparML = 0.
      null_vector = 0.
      cntp_mar = 0
      T_mar = 0.
      S_mar = 0.
      U_mar =  0.
      V_mar =  0.
      Ku_mar =  0.
      Ks_mar =  0.
      Kt_mar =  0.
      cntp_jun = 0
      T_jun = 0.
      S_jun = 0.
      U_jun =  0.
      V_jun =  0.
      Ku_jun =  0.
      Ks_jun =  0.
      Kt_jun =  0.
      cntp_sep = 0
      T_sep = 0.
      S_sep = 0.
      U_sep =  0.
      V_sep =  0.
      Ku_sep =  0.
      Ks_sep =  0.
      Kt_sep =  0.
      cntp_dec = 0
      T_dec = 0.
      S_dec = 0.
      U_dec =  0.
      V_dec =  0.
      Ku_dec =  0.
      Ks_dec =  0.
      Kt_dec =  0.
      
END SUBROUTINE initialize




