! $Id$
! $Source$

SUBROUTINE initialize

      USE mean_param
      USE surface_forcing
      USE declarations

      IMPLICIT NONE

      INTEGER::k_k

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
!       PO_flux = 0.
!       NO_flux = 0.
!       fratio = 0.
       f_ratio = 0.
      null_vector = 0.
end subroutine initialize
