! $Id$
! $Source$

subroutine derivs_noflag(x_var, nvar1, Y_var, DYDX_deriv, Temp)
  ! Calculate the derivatives of the biological model for odeint to
  ! use to advance the biology to the next time step.
  !
  ! *** This version does not contain flagellates code. ***
  ! *** It should not be needed once the flagellates model in
  ! *** derivs_sog is working.

  use surface_forcing
  use declarations
  use reaction
  use pdf

  implicit none

  ! Arguments:
  INTEGER, INTENT(IN)::nvar1 !5*M+bin_tot + Csources*M+D_bins*M
  DOUBLE PRECISION, DIMENSION(nvar1), INTENT(IN)::Y_var  !INTENT(IN OUT)
  DOUBLE PRECISION, DIMENSION(nvar1), INTENT(OUT)::DYDX_deriv
  DOUBLE PRECISION, INTENT(IN)::x_var  !Y_var known at x_var

  INTEGER::j_j, k_k,bin_tot2,bin,isusan
  DOUBLE PRECISION, DIMENSION(nvar1)::Yplus
  DOUBLE PRECISION, DIMENSION(1:grid%M)::PO_deriv,pellet_prod!,N_urea_save
  DOUBLE PRECISION, DIMENSION(1:grid%M)::Temp,fratioavg
  DOUBLE PRECISION::bottom,bottom2
  DOUBLE PRECISION, DIMENSION(2*grid%M)::Resp


  bin_tot2 = M2
  bin = M2

  DO j_j = 1,nvar1  
     IF (Y_var(j_j) >= 0.) THEN 
        Yplus(j_j) = Y_var(j_j)
     ELSE
        Yplus(j_j) = 0.
     END IF
  END DO

  N%O_uptake%new = 0.
  N%H_uptake%new = 0.
  waste%medium = 0.
  N%remin = 0.

  CALL p_growth(Yplus,Yplus(1:M),nvar1,I_par,grid,N,micro,Temp,Resp) !microplankton
  !      waste%medium = waste%medium + micro%M_z*Yplus(1:M)
  waste%medium = waste%medium + Resp(M+1:2*M)*Yplus(1:M)

  DYDX_deriv(1:nvar1) = 0.      

!!!New quantity, bacterial 0xidation of NH to NO pool ==> NH^2
  N%bacteria(1:M) = N%r*Yplus(2*M+1:3*M)**2.0


  DO k_k = 1,D_bins-1
     N%remin(:) = N%remin(:) + Detritus(k_k)%r*Yplus(3*M+(k_k-1)*M+1:3*M+k_k*M)
  END DO
  IF (MINVAL(N%remin) < 0.) THEN
     PRINT "(A)","N%remin < 0. in derivs.f90"
     PRINT *,N%remin
     STOP
  END IF

  DO j_j = 1,M 
     IF (N%O_uptake%new(j_j)+N%H_uptake%new(j_j) > 0.) THEN
        f_ratio(j_j) = N%O_uptake%new(j_j)/(N%O_uptake%new(j_j)+N%H_uptake%new(j_j))

        !fratioavg(j_j) = f_ratio(j_j)*(N%H_uptake%new(j_j)+N%O_uptake%new(j_j))
     ELSE
        f_ratio(j_j) = 0.
     END IF
  END DO


  do j_j = 1, nvar1
     if(j_j  <= M .AND. Yplus(j_j) > 0.) then 
        !Pmicro mortality
        DYDX_deriv(j_j) = (micro%growth%new(j_j) - Resp(j_j) &
             - Resp(j_j+M)) * Yplus(j_j)
        ! DYDX_deriv(j_j) = (micro%growth%new(j_j)-Resp(j_j)-micro%M_z)*Yplus(j_j)
     else if (j_j > M .AND. j_j <= 2*M .AND. Yplus(j_j) > 0.) then  
        ! NO
        DYDX_deriv(j_j) = -N%O_uptake%new(j_j-M) + N%bacteria(j_j-M)
     else if (j_j > 2*M .AND. j_j <= 3*M) then  
        ! NH
        DYDX_deriv(j_j) = -N%H_uptake%new(j_j-2*M) + Resp(j_j-2*M) &
             + waste%medium(j_j-2*M) * waste%m%destiny(0)          &
             + N%remin(j_j-2*M) - N%bacteria(j_j-2*M) 
     else if (j_j > 3*M .AND. j_j <= 3*M+D_bins*M) then
        ! Detritus
        do k_k = 1, D_bins-1
           if (j_j > 3 * M + (k_k-1) * M .AND. j_j <= 3 * M + k_k * M) then
              DYDX_deriv(j_j) = waste%medium(j_j-(3*M+(k_k-1)*M)) &
                   * waste%m%destiny(k_k) - Detritus(k_k)%r * Yplus(j_j) 
           end if
        end do
     end if
  end do


  PO_deriv(1:M) = DYDX_deriv(1:M)

  ! Integrate rate of urea production
  ! *** More hard-coded constants to get rid of
  ! *** These refer to the depth of the model domain and have caused
  ! *** trouble with index out of bounds runtime errors in sum_g()
  ! In the top 5 meters
  bottom = 5.  
  call sum_g(grid, DYDX_deriv(M+1:2*M), bottom, NO50_rate)
  call sum_g(grid, PO_deriv(1:M), bottom, PO50_rate)
  ! In the whole model domain
  bottom = 40.
  ! *** These calls produce array index out of bounds error in sum_g
  ! *** and they may have nothing to do with phytoplankton so 
  ! *** don't do them for now
!!$  call sum_g(grid, DYDX_deriv(M+1:2*M), bottom, NO100_rate)
!!$  call sum_g(grid, PO_deriv(1:M), bottom, PO100_rate)

end subroutine derivs_noflag
