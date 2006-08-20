! $Id$
! $Source$

subroutine derivs_noflag(time, M2, PZ, dPZdt, Temp)
  ! Calculate the derivatives of the biological model for odeint to
  ! use to advance the biology to the next time step.
  !
  ! *** This version does not contain flagellates code. ***
  ! *** It should not be needed once the flagellates model in
  ! *** derivs_sog is working.

  use mean_param, only: bins, losses, plankton2, nutrient, snow
  use declarations, only: M, D_bins, N, micro, f_ratio, waste, Detritus, &
       I_par, grid
  use reaction, only: p_growth

  implicit none

  ! Arguments:
  DOUBLE PRECISION, INTENT(IN):: time ! not used
  INTEGER, INTENT(IN):: M2 ! size of PZ
  DOUBLE PRECISION, DIMENSION(M2), INTENT(IN):: PZ  ! values
  DOUBLE PRECISION, DIMENSION(M2), INTENT(OUT)::dPZdt ! derivatives
  double precision, dimension(1:M), INTENT(IN):: Temp ! temperature

  type (bins) :: PZ_bins
  common /derivs/  PZ_bins

  ! Local variables

  integer :: ii ! counter through grid ie 1 to M
  integer :: jj ! counter through PZ
  integer :: kk ! counter through detritus

  DOUBLE PRECISION, DIMENSION(2*M):: Resp ! respiration and mortality

  double precision, dimension (M) :: Pmicro ! microplankton concentrations
  double precision, dimension (M) :: NO, NH ! nitrate and ammonium conc.
  double precision, dimension (D_bins, M) :: detr ! detritus

! put PZ micro values into Pmicro variable, removing any negative values

  do ii = 1,M                        ! counter through grid
     jj = (PZ_bins%micro-1) * M + ii ! counter into PZ

     if (PZ(jj)>0) then
        Pmicro(ii) = PZ(jj)
     else
        Pmicro(ii) = 0
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
  N%remin = 0.

! phytoplankton growth: Nitrate and Ammonimum
!*** I_par is light (IN) , grid is used for M only I think (should be changed)
! N ammonium and nitrate uptake rates (IN and OUT)
! micro is the growth numbers for micro plankton (IN and OUT)
! Temp is temperature (IN)
! Resp is respiration values (OUT)

  call p_growth(NO,NH,Pmicro,M2,I_par,grid,N,micro,Temp,Resp) !microplankton

! put microplankton mortality into the medium detritus flux
  waste%medium = waste%medium + Resp(M+1:2*M)*Pmicro

!!!New quantity, bacterial 0xidation of NH to NO pool ==> NH^2
  N%bacteria(1:M) = N%r*NH**2.0

! remineralization of detritus groups 1 and 2 (not last one)
  do kk = 1,D_bins-1
     N%remin(:) = N%remin(:) + Detritus(kk)%r*detr(kk,:)
  enddo

  IF (MINVAL(N%remin) < 0.) THEN
     PRINT "(A)","N%remin < 0. in derivs.f90"
     PRINT *,N%remin
     STOP
  END IF

! calculate the f-ratio
  do jj = 1,M 
     IF (N%O_uptake%new(jj)+N%H_uptake%new(jj) > 0.) THEN
        f_ratio(jj) = N%O_uptake%new(jj)/(N%O_uptake%new(jj)+N%H_uptake%new(jj))
     ELSE
        f_ratio(jj) = 0.
     END IF
  END DO

! now put it all together into the derivatives

! initialize derivatives
  dPZdt(1:M2) = 0.      

! microplankton

  do ii = 1,M ! index for Pmicro etc
     jj = (PZ_bins%micro-1) * M + ii ! index for PZ, derivatives etc

     if (Pmicro(ii) > 0.) then 
        dPZdt(jj) = (micro%growth%new(ii) - Resp(ii) &
             - Resp(ii+M)) * Pmicro(ii)
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
        dPZdt(jj) = -N%H_uptake%new(ii) + Resp(ii) &
             + waste%medium(ii) * waste%m%destiny(0)          &
             + N%remin(ii) - N%bacteria(ii)
     endif
  enddo

! detritus (not last bin)

  do ii = 1, M
     do kk = 1, D_bins-1
        jj = (PZ_bins%Quant + (kk-1)) * M + ii
        
        if (detr(kk,ii) > 0) then
           dPZdt(jj) = waste%medium(ii) * waste%m%destiny(kk)  &
                - Detritus(kk)%r * detr(kk,ii) 
        end if
     end do
  end do

end subroutine derivs_noflag
