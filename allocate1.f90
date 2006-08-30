! $Id$
! $Source$

subroutine allocate1(STAT)

      use mean_param
      use declarations
      use surface_forcing

      implicit none
     
      ! Argument:
      integer, dimension(20), intent(OUT)::STAT

      ALLOCATE(U%new(0:M+1),U%old(0:M+1),U%old_old(0:M+1),U%last(0:M+1),V%new(0:M+1),V%old(0:M+1),&
               V%old_old(0:M+1),V%last(0:M+1),&
               T%new(0:M+1),T%old(0:M+1),T%old_old(0:M+1), S%new(0:M+1),&
               S%old(0:M+1),S%old_old(0:M+1),B%new(0:M+1),B%old(0:M+1),P%micro%old(0:M+1),&
               PON%new(0:M+1),PON%old(0:M+1),&
               P%micro%new(0:M+1),Z%micro%new(0:M+1),Z%micro%old(0:M+1),P%micro%old_old(0:M+1),&
               Z%micro%old_old(0:M+1), P%nano%new(0:M+1),P%nano%old(0:M+1),P%nano%old_old(0:M+1),&
               N%O%new(0:M+1),N%O%old(0:M+1),N%O%old_old(0:M+1),&
               N%H%new(0:M+1),N%H%old(0:M+1),N%H%old_old(0:M+1),&
               Sil%new(0:M+1), Sil%old(0:M+1), &            ! silicon
               P_temp%new(0:M+1),density%old(0:M+1),density%new(0:M+1),dens_i(0:M+1),&
               STAT = alloc_stat(1))

      ALLOCATE(Detritus(D_bins),Hvector%d(D_bins),Gvector%d(D_bins),Gvector_o%d(D_bins),&
              Gvector_o_o%d(D_bins),Gvector_ro%d(D_bins),Gvector_ro_o%d(D_bins),&
              Gvector_ao%d(D_bins),Gvector_ao_o%d(D_bins),Detritus1_p(D_bins,M),&
              STAT = alloc_stat(3))

      ALLOCATE(alph%g(0:M+1), alph%i(0:M), alph%idiv(M), &
           beta%idiv(M), beta%g(0:M+1), beta%i(0:M),     &
           w%b(0:M),w%b_err(0:M),Bf%b(0:M),&
               Bf%b_err(0:M),w%t(0:M),w%s(0:M),w%u(0:M),&
               w%v(0:M),w%p%micro(0:M),Amatrix%u%A(M),Amatrix%u%B(M),&
               Amatrix%u%C(M),Amatrix%t%A(M),Amatrix%t%B(M),Amatrix%t%C(M),Amatrix%s%A(M),Amatrix%s%B(M),&
               Amatrix%s%C(M),Bmatrix%u%A(M),Bmatrix%u%B(M),&
               Bmatrix%u%C(M),Bmatrix%t%A(M),Bmatrix%t%B(M),Bmatrix%t%C(M),Bmatrix%s%A(M),Bmatrix%s%B(M),&
               Bmatrix%s%C(M),Bmatrix_o%u%A(M),Bmatrix_o%u%B(M),&
               Bmatrix_o%u%C(M),Bmatrix_o%t%A(M),Bmatrix_o%t%B(M),Bmatrix_o%t%C(M),Bmatrix_o%s%A(M),&
               Bmatrix_o%s%B(M),Bmatrix_o%s%C(M),Bmatrix_o_o%u%A(M),Bmatrix_o_o%u%B(M),&
               Bmatrix_o_o%u%C(M),Bmatrix_o_o%t%A(M),Bmatrix_o_o%t%B(M),Bmatrix_o_o%t%C(M),&
               Bmatrix_o_o%s%A(M),Bmatrix_o_o%s%B(M),Bmatrix_o_o%s%C(M),&
               Amatrix%bio%A(M),Amatrix%bio%B(M),Amatrix%bio%C(M),&
               Bmatrix%bio%A(M),Bmatrix%bio%B(M),Bmatrix%bio%C(M),&
               Bmatrix%null%A(M),Bmatrix%null%B(M),Bmatrix%null%C(M),&
               Amatrix%null%A(M),Amatrix%null%B(M),&
               Bmatrix_o%bio%A(M),Bmatrix_o%bio%B(M),Bmatrix_o%bio%C(M),&
               Bmatrix_o_o%bio%A(M),Bmatrix_o_o%bio%B(M),Bmatrix_o_o%bio%C(M),&
               Amatrix%no%A(M),Amatrix%no%B(M),Amatrix%no%C(M),&
               Bmatrix%no%A(M),Bmatrix%no%B(M),Bmatrix%no%C(M),&
               Bmatrix_o%no%A(M),Bmatrix_o%no%B(M),Bmatrix_o%no%C(M),&
               Bmatrix_o_o%no%A(M),Bmatrix_o_o%no%B(M),Bmatrix_o_o%no%C(M),&
               STAT = alloc_stat(4)) !**&

      ALLOCATE(T%div_g(M),U%div_g(M),V%div_g(M),S%div_g(M),B%div_g(M),density%div_g(M),T%div_i(M),&
               U%div_i(M),V%div_i(M),S%div_i(M),density%div_i(M),B%div_i(M),P%micro%div_i(M), &
               K%u%shear(0:M),K%s%dd(0:M),K%t%dd(0:M),&
               K%u%total(0:M),K%s%total(0:M),K%t%total(0:M),K%t%all(0:M),K%u%all(0:M),K%s%all(0:M),&
               K%t%old(0:M),K%s%old(0:M),K%u%old(0:M),STAT = alloc_stat(5))

      ALLOCATE(ref_T(0:M+1), &
           avg_12(0:M+1), tot_avg(0:M+1), Q_t(0:M), T_To(0:M+1), &
           I(0:M), I_par(0:M), Q_n(0:M), F_n(0:M),  grid%d_g(0:M+1), &
           grid%d_i(0:M), grid%i_space(M), grid%g_space(0:M), &
           phi%m%value(0:M), phi%s%value(0:M), Ri_b(M), &
           N_2_i(M), N_2_g(M), N_2_dens_g(M), &
           Q_test(M), V_t_square(M), omega%s%value(0:M), omega%m%value(0:M), &
           STAT = alloc_stat(6)) 
                     
      ALLOCATE(gamma%m(0:M),gamma%s(0:M),gamma%t(0:M), T_mar(M),T_jun(M),T_sep(M),T_dec(M),&
           S_mar(M),S_jun(M),S_sep(M),S_dec(M),STAT = alloc_stat(7))

      ALLOCATE(Hvector%s(M),Hvector%t(M),Hvector%u(M),Hvector%v(M),Hvector%p%micro(M),Hvector%p%nano(M),&
               Hvector%z%micro(M),&
               Gvector%s(M),Gvector%t(M),Gvector%u(M),Gvector%v(M),Gvector_o%s(M),Gvector_o%t(M),&
               Gvector_o%u(M),Gvector_o%v(M),Gvector_o_o%s(M),Gvector_o_o%t(M),&
               Gvector_o_o%u(M),Gvector_o_o%v(M),Gvector_c%u(M),Gvector_c%v(M),Gvector_co%u(M),&
               Gvector_co%v(M),Gvector_co_o%u(M),Gvector_co_o%v(M),&
               Gvector%p%micro(M),Gvector_o%p%micro(M),Gvector_o_o%p%micro(M),&
               Gvector%p%nano(M),Gvector_o%p%nano(M),Gvector_o_o%p%nano(M),Gvector%z%micro(M),&
               Gvector_o%z%micro(M),Gvector_o_o%z%micro(M), Gvector_ao%p%micro(M), &
               Gvector_ao_o%p%micro(M), Gvector_ro%p%micro(M),Gvector_ro_o%p%micro(M),&
               Gvector_ro%p%nano(M),Gvector_ro_o%p%nano(M),&
               Gvector_ro%z%micro(M),Gvector_ro_o%z%micro(M), Gvector_ao%z%micro(M), &
               Gvector_ao_o%z%micro(M), &
               Hvector%n%o(M), Hvector%n%h(M), Gvector%n%o(M), Gvector%n%h(M), &
               Gvector%sil(M), &
               Gvector_o%n%o(M), Gvector_o%n%h(M), Gvector_o_o%n%o(M), Gvector_o_o%n%h(M),&
               Gvector_ro%n%o(M), Gvector_ro%n%h(M), Gvector_ro%sil(M), Gvector_ro_o%n%o(M), Gvector_ro_o%n%h(M),&
               null_vector(M),STAT = alloc_stat(8)) !**&
               
      ALLOCATE(U_p(M),V_p(M),S_p(M),T_p(M),P1_p(M),Pnano1_p(M),Z1_p(M),NO1_p(M),NH1_p(M), & 
           wind(wind_n), insol(insol_n), STAT = alloc_stat(9))

      ALLOCATE(micro%growth%light(M),micro%growth%new(M),nano%growth%light(M),nano%growth%new(M),&
           N%O_uptake%new(M),N%H_uptake%new(M),N%urea%new(M),N%remin(M),N%bacteria(M), &
           micro%mort%new(M),  nano%mort%new(M),&
           zmicro%growth%new(M), zmicro%mort%new(M),&
           zmicro%graze(zprey,M),zmicro%q(zprey-1), waste%small(M),waste%medium(M), &
           waste%large(M),waste%s%destiny(0:D_bins),waste%m%destiny(0:D_bins),waste%l%destiny(0:D_bins), &
           P_no(M),P_nh(M),P_di(M),P_na(M),P_mi(M),P_d1(M),P_d2(M),P_sa(M),P_u(M),P_v(M),P_ta(M),wupwell(M+1), &
           STAT = alloc_stat(10))

!Copepod allocations and plot allocations

      ALLOCATE(f_ratio(M), &
           nano_pro(84,M),diatom_pro(84,M), zmicro_pro(84,M), copepod_pro(84,M), &
           NO_pro(84,M), NH_pro(84,M), fratio_pro(84,M), T_pro(84,M), S_pro(84,M), &
           U_pro(84,M), V_pro(84,M), PON_pro(84,M), &
           nano_mar(M),nano_jun(M),nano_sep(M),nano_dec(M),diatom_mar(M),&
           diatom_jun(M),diatom_sep(M),diatom_dec(M),zmicro_mar(M),zmicro_sep(M),zmicro_jun(M),&
           zmicro_dec(M),copepod_mar(M),copepod_sep(M),copepod_jun(M),copepod_dec(M),NO_mar(M),&
           NO_jun(M),NO_sep(M),NO_dec(M),NH_mar(M),NH_jun(M),NH_sep(M),NH_dec(M),don_jun(M),don_mar(M),&
           don_sep(M),don_dec(M),pon_jun(M),pon_mar(M),pon_sep(M),pon_dec(M),&
           NPPn_mar(M),NPPn_jun(M),NPPn_dec(M),&
          NPPn_sep(M),NPPd_mar(M),NPPd_jun(M),NPPd_dec(M),NPPd_sep(M),ngrow_mar(M),ngrow_jun(M),&
          ngrow_sep(M),ngrow_dec(M),dgrow_mar(M),&
          dgrow_jun(M),dgrow_sep(M),dgrow_dec(M),ngraz_mar(M),ngraz_jun(M),ngraz_sep(M),ngraz_dec(M),&
          dgraz_mar(M),dgraz_jun(M),&
          dgraz_dec(M),dgraz_sep(M),NOup_mar(M),NOup_sep(M),NOup_jun(M),NOup_dec(M),NHup_mar(M),&
          NHup_jun(M),NHup_sep(M),NHup_dec(M),&
          fratio_mar(M),fratio_jun(M),fratio_dec(M),fratio_sep(M),&
          U_mar(M),U_dec(M),U_sep(M),U_jun(M),V_mar(M),V_jun(M),V_sep(M),V_dec(M),&
          Ku_mar(M),Ku_jun(M),Ku_sep(M),Ku_dec(M),&
          Ks_mar(M),Ks_dec(M),Ks_sep(M),Ks_jun(M),Kt_mar(M),Kt_dec(M),Kt_sep(M),&
          Kt_jun(M),STAT = alloc_stat(11))

! sea add
      ALLOCATE (ut%new(M),vt%new(M),pbx(M),pby(M),&
                ut%old(M),vt%old(M),dzx(M),dzy(M), STAT=alloc_stat(12))

               
      
END SUBROUTINE allocate1



