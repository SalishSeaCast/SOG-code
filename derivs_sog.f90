! $Id$
! $Source$

subroutine derivs_sog(time, M2, PZ, dPZdt, Temp)
  ! Calculate the derivatives of the biological model for odeint to
  ! use to advance the biology to the next time step.

  use precision_defs, only: dp
  use mean_param, only: bins, losses, plankton2, nutrient, snow
  use declarations, only: M, D_bins, N, micro, nano, f_ratio, waste, &
       Detritus, I_par
  use reaction, only: p_growth

  implicit none

  ! Arguments:
  real(kind=dp), INTENT(IN):: time ! not used
  INTEGER, INTENT(IN):: M2 ! size of PZ
  real(kind=dp), DIMENSION(M2), INTENT(IN):: PZ  ! values
  real(kind=dp), DIMENSION(M2), INTENT(OUT)::dPZdt ! derivatives
  real(kind=dp), dimension(0:M), INTENT(IN):: Temp ! temperature

  type (bins) :: PZ_bins
  common /derivs/  PZ_bins

  ! Local variables

  integer :: ii ! counter through grid ie 1 to M
  integer :: jj ! counter through PZ
  integer :: kk ! counter through detritus

  real(kind=dp), DIMENSION(M):: Resp_micro, Mort_micro! respiration & mortality
  real(kind=dp), DIMENSION(M):: Resp_nano, Mort_nano  ! respiration & mortality

  real(kind=dp), dimension (M) :: Pmicro, Pnano ! micro/nano plankton conc.
  real(kind=dp), dimension (M) :: NO, NH ! nitrate and ammonium conc.
  real(kind=dp), dimension (D_bins, M) :: detr ! detritus

! put PZ micro values into Pmicro variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%micro-1) * M + ii ! counter into PZ

     if (PZ(jj)>0) then
        Pmicro(ii) = PZ(jj)
     else
        Pmicro(ii) = 0
     endif
  enddo

! put PZ nano values into Pnano variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%nano-1) * M + ii ! counter into PZ

     if (PZ(jj)>0) then
        Pnano(ii) = PZ(jj)
     else
        Pnano(ii) = 0
     endif
  enddo


! put PZ nitrate values into NO variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%NO-1) * M + ii    ! counter into PZ
     
     if (PZ(jj)>0) then
        NO(ii) = PZ(jj)
     else
        NO(ii) = 0
     endif
  enddo

! put PZ ammonium values into NH variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%NH-1) * M + ii    ! counter into PZ
     
     if (PZ(jj)>0) then
        NH(ii) = PZ(jj)
     else
        NH(ii) = 0
     endif
  enddo

! put PZ detritus values into detr  variable, removing any negative values

  do kk = 1,D_bins
     do ii = 1,M                                 ! counter through grid
        jj = (PZ_bins%Quant + (kk-1)) * M + ii   ! counter through PZ

        if (PZ(jj)>0) then
           detr(kk,ii) = PZ(jj)
        else
           detr(kk,ii) = 0
        endif
     enddo
  enddo

! initialize transfer rates between the pools
  N%O_uptake%new = 0.
  N%H_uptake%new = 0.
  waste%medium = 0.
  waste%small = 0.
  N%remin = 0.

! phytoplankton growth: Nitrate and Ammonimum, conc. of micro plankton
! I_par is light, Temp is temperature 
! N ammonium and nitrate uptake rates (IN and OUT) incremented
! micro is the growth parameters for micro plankton (IN) and the growth rates 
! (OUT)
! Resp is respiration, Mort is mortality

  call p_growth(M, NO, NH, Pmicro, I_par, Temp, & ! in
       N, micro, &                                ! in and out
       Resp_micro, Mort_micro)                    ! out

! put microplankton mortality into the medium detritus flux
  waste%medium = waste%medium + Mort_micro*Pmicro

! phytoplankton growth: Nitrate and Ammonimum, conc. of nano plankton
! I_par is light, Temp is temperature 
! N ammonium and nitrate uptake rates (IN and OUT) incremented
! nano is the growth parameters for micro plankton (IN) and the growth rates 
! (OUT)
! Resp_nano is respiration, Mort_nano is mortality

  call p_growth(M, NO, NH, Pnano, I_par, Temp, & ! in
       N, nano, &                                ! in and out
       Resp_nano, Mort_nano)                     ! out

  waste%small = waste%small + Mort_nano*Pnano

!!!New quantity, bacterial 0xidation of NH to NO pool ==> NH^2
  N%bacteria(1:M) = N%r * NH**2

! remineralization of detritus groups 1 and 2 (not last one)
  do kk = 1,D_bins-1
     N%remin(:) = N%remin(:) + Detritus(kk)%r * detr(kk,:)
  enddo

  IF (MINVAL(N%remin) < 0.) THEN
     PRINT "(A)","N%remin < 0. in derivs.f90"
     PRINT *,N%remin
     STOP
  END IF

! calculate the f-ratio
  do jj = 1,M 
     IF (N%O_uptake%new(jj) + N%H_uptake%new(jj) > 0.) THEN
        f_ratio(jj) = N%O_uptake%new(jj) /  &
             (N%O_uptake%new(jj) + N%H_uptake%new(jj))
     ELSE
        f_ratio(jj) = 0.
     END IF
  END DO

! now put it all together into the derivatives

! initialize derivatives
  dPZdt(1:M2) = 0.      

! microplankton

  do ii = 1,M ! index for Pmicro
     jj = (PZ_bins%micro-1) * M + ii ! index for PZ, derivatives etc

     if (Pmicro(ii) > 0.) then 
        dPZdt(jj) = (micro%growth%new(ii) - Resp_micro(ii) &
             - Mort_micro(ii)) * Pmicro(ii)
! *** problem here, Resp and Mort are multiplied by Pmicro here but not in NH
! *** or detritus derivatives
! okay here, problem is below.
     endif
  enddo

! nanoplankton

  do ii = 1,M ! index for Pnano
     jj = (PZ_bins%nano-1) * M + ii ! index for PZ, derivatives etc

     if (Pnano(ii) > 0.) then 
        dPZdt(jj) = (nano%growth%new(ii) - Resp_nano(ii) &
             - Mort_nano(ii)) * Pnano(ii)
     endif
  enddo


! nitrate

  do ii = 1,M ! index for NO
     jj = (PZ_bins%NO-1) * M + ii ! index for PZ, derivatives etc

     if (NO(ii) > 0.) then  
        dPZdt(jj) = -N%O_uptake%new(ii) + N%bacteria(ii)
     endif
  enddo

! ammonium

  do ii = 1,M ! index for NH
     jj = (PZ_bins%NH-1) * M + ii ! index for PZ, derivatives etc
     
     if (NH(ii) > 0.) then  
        dPZdt(jj) = -N%H_uptake%new(ii) + Resp_micro(ii) + Resp_nano(ii) &
             + waste%medium(ii) * waste%m%destiny(0)                     &
             + waste%small(ii) * waste%s%destiny(0)                      &
             + N%remin(ii) - N%bacteria(ii)
     endif
  enddo

! detritus (not last bin)

  do ii = 1, M
     do kk = 1, D_bins-1
        jj = (PZ_bins%Quant + (kk-1)) * M + ii
        
        if (detr(kk,ii) > 0) then
           dPZdt(jj) = waste%medium(ii) * waste%m%destiny(kk)  &
                + waste%small(ii) * waste%s%destiny(kk) &
                - Detritus(kk)%r * detr(kk,ii) 
        end if
     end do
  end do

end subroutine derivs_sog
