! $Id$
! $Source$

module surface_forcing
  ! *** This is a very poorly named module!

  use precision_defs, only: dp
  use mean_param

  implicit none

  real(kind=dp), parameter :: PI = 3.141592653589793
  ! Coriolis factor
  ! *** This should be reworked so that the latitude is read as a run
  ! *** parameter, and f calculated for the run.
  real(kind=dp), parameter :: latitude = 49. + 7.517 / 60. ! station S3
!*** Correction pgf90 will not except an intrinsic in parameter statement  
!*** real(kind=dp), parameter :: f = 2. * (2. * PI / 86400.) * &
!***    sin(PI * latitude / 180.)
  real(kind=dp), parameter :: f = 1.1d-4
  ! Acceleration due to gravity
  real(kind=dp), parameter :: g = 9.81  ! m/s^2
  ! Critical value of Richardson number for mixed layer depth
  ! determination
  ! *** Susan was surprised that this value was not 0.25
  real(kind=dp), parameter :: Ri_c = 0.3

  DOUBLE PRECISION, PARAMETER:: &
       small = 1.D-15, &
       small2 = 0.001, &  !smallest number of copepods out
       min_out = 2.5D-04, &   !1.D-05, &
       cutoff = 10.0, &
       var_c = 5.0, &   !***  !11.0
       zero =  0.
  DOUBLE PRECISION, PARAMETER:: Cv = 1.5  !1.99 !set to keep beta_t = -0.2
       !beta_t = -0.2 only for convection
  DOUBLE PRECISION, PARAMETER::S_ice = 0.0, & !ice salinity. unknown
       rho_ice = 1000.0, & !density ice, unknown
       kapa = 0.4 !Von Karman constant
  DOUBLE PRECISION, PARAMETER::a_s = -28.86, a_m = 1.26, &
       c_s =  98.96, c_m = 8.38, &
       xsi_m = -0.20, xsi_s = -1.0,&
       ep = 0.1  !epsilon
  DOUBLE PRECISION, PARAMETER::Ri_o = 0.7, & !critical grad Richardson #
       nu_o = 0.001, & !0.0050 !m^2/s interior diff constant
       p_1 = 3.0, & ! shear diff power constant
       !           nu_w_m = 5.0D-04, & !winter m^2/s Wave breaking momentum
  nu_w_m = 1.0D-04, & !summer m^2/s Wave breaking momentum
       nu_w_s = 1.0D-05, & !summer 5.0D-05, &   !and scalar interior diff
       !           nu_w_s = 0.5D-05, & !winter 5.0D-05, &   !and scalar interior diff
  nu_f = 0.001, & ! m^2/s salt fingering constant 
       Ri_rho_o = 1.9, &
       p_2 = 3.0, & !salt finger diff power constant
       nu = 1.5D-06 !m^2/s molecular viscosity
  DOUBLE PRECISION, PARAMETER::C_star = 9.9
  DOUBLE PRECISION, PARAMETER::del_o = 0.10, &!tolerance
       A_stress = -0.6, & !N/m^2 Wind stress constant
       T_stress = 57600.0, &  !seconds or 16 hrs
       Q_o = 1368.0, & !1367.0? W/m^2  Solar constant
       Lat = 0.85521, &  ! 49 degrees  latitude
       Lon = 145., &  !(145oW)
       !           albedo = 0.061, &  !6% Large 1996
  albedo = 0.18, &  !KC 17% OCT.22 2004
       emiss = -1.0, & !emissivity Large 1996
       Stef_Boltz = 5.6697D-08, & !W/m^2/K^4
       R_v = 287.04, & !J/kg/K  gas constant for dry air
       mu_3 = 0.2, &   !(Denman)^(-1)  gN/m^2  attenuation coeff  
       handle = 1.0D-02  !gN(P)/gN(Z)/(handle time 's')   swim diffusivity
  !***
  DOUBLE PRECISION, DIMENSION(13), PARAMETER::leap_year = (/ 1956, 1960, 1964, 1968, 1972,1976, 1980, &
       1984, 1988, 1992, 1996, 2000, 2004 /)
  INTEGER, PARAMETER:: runsize = 15 ! runsize = 15 !, &  !used in smoothdata.f90
  ! bin = 50 !50, &  !number of copepod weight bins in wt_pdf
  DOUBLE PRECISION, PARAMETER::min_bin = 0.005, & !minimum bin size for copepod wt_pdf is 0.0001*Zoo%wt
       sigma_max = 10.  !12 is too big ==> stage 1's increase slightly. see bye_zoo
  DOUBLE PRECISION, PARAMETER::Q_max = 8.,&   !maximum and minimum allowed CN ratios
       Q_min = 5. !4.                             !

  INTEGER :: zprey  !number of Copepod animal prey (P%micro+Z%micro) and 
  !detrital prey (Detritus(1:D_bins-1)
  !See read_data.f90 and input/biology for values
  DOUBLE PRECISION, PARAMETER::desired_var = 6.267101794384283, &  !12.6, & !6.267101794384283, & !used for ECMWF smoothed winds
       desired_avg =   10.06734127539335 
  DOUBLE PRECISION :: &
       precision = 1.0D-3, &
       step_min = 30., &
       step_guess = 3600.0

CONTAINS

  SUBROUTINE average(mm,X,surf_h)

    TYPE(gr_d), INTENT(IN)::mm 
    TYPE(prop), INTENT(IN OUT)::X   !U, V, T ...
    TYPE(height), INTENT(IN OUT)::surf_h  !surface_height   

    DOUBLE PRECISION::X_sl,X_h 
    INTEGER::k

    CALL find_jmax_g(surf_h,mm)

    X%avg = 0.


    IF (surf_h%g <= 1) THEN
       X%avg = X%new(1)
    ELSE
       X_h = X%new(1)*mm%g_space(0)
       IF (surf_h%g > 2) THEN
          DO k = 2, surf_h%g - 1
             X_h = X_h + (X%new(k-1) + X%new(k))*mm%g_space(k-1)/2.0 
          END DO
       END IF
       X_sl = X%new(surf_h%g-1) + (surf_h%new-mm%d_g(surf_h%g-1))*&
            (X%new(surf_h%g)-X%new(surf_h%g-1))/mm%g_space(surf_h%g-1)
       X%avg = (X_h + (X%new(surf_h%g-1) + X_sl)*(surf_h%new-mm%d_g(surf_h%g-1))/2.0)/&
            surf_h%new
    END IF


  END SUBROUTINE average


  SUBROUTINE average3(mm,Xnew,bottom,avg_X)

    TYPE(gr_d), INTENT(IN)::mm 
    DOUBLE PRECISION, DIMENSION(mm%M), INTENT(IN)::Xnew   !K%new(mm%M)
    DOUBLE PRECISION, INTENT(IN)::bottom  !bottom depth
    DOUBLE PRECISION, INTENT(OUT)::avg_X 

    DOUBLE PRECISION::X_sl,X_h 
    INTEGER::k,bottom_g

    DO k = 1,mm%M
       IF (mm%d_i(k) > bottom) THEN
          bottom_g = k
          EXIT
       ELSE
          bottom_g = k+1
       END IF
    END DO

    avg_X = 0.

    IF (bottom_g <= 1) THEN
       avg_X = Xnew(1)
    ELSE
       X_h = Xnew(1)*mm%i_space(0)
       IF (bottom_g > 2) THEN
          DO k = 2, bottom_g - 1
             X_h = X_h + (Xnew(k-1) + Xnew(k))*mm%i_space(k-1)/2.0 
          END DO
       END IF
       X_sl = Xnew(bottom_g-1) + (bottom-mm%d_i(bottom_g-1))*&
            (Xnew(bottom_g)-Xnew(bottom_g-1))/mm%i_space(bottom_g-1)
       avg_X = (X_h + (Xnew(bottom_g-1) + X_sl)*(bottom-mm%d_i(bottom_g-1))/2.0)/&
            bottom
    end if
  end subroutine average3


  subroutine sum_g(mm, Xnew, bottom, sum_X)
    ! *** What's it do?

    ! Arguments:
    type(gr_d), intent(in) :: mm 
    double precision, dimension(mm%M), intent(in) :: Xnew   ! U%new(mm%M)
    double precision, intent(in) :: bottom                  ! bottom depth
    double precision, intent(out) :: sum_X 

    ! Local variables:
    double precision :: X_sl, X_h 
    integer :: k, bottom_g

    do k = 1, mm%M
       if (mm%d_g(k) > bottom) then
          bottom_g = k
          exit
       else
          bottom_g = k+1
       end if
    end do

    sum_X = 0.

    if (bottom_g <= 1) then
       sum_X = Xnew(1) * bottom
    else
       X_h = Xnew(1) * mm%g_space(0)
       if (bottom_g > 2) then
          do k = 2, bottom_g - 1
             X_h = X_h + (Xnew(k-1) + Xnew(k)) * mm%g_space(k-1) / 2.0 
          end do
       end if
       X_sl = Xnew(bottom_g - 1) + (bottom - mm%d_g(bottom_g - 1)) &
            * (Xnew(bottom_g) - Xnew(bottom_g - 1))                &
            / mm%g_space(bottom_g - 1)
       sum_X = (X_h + (Xnew(bottom_g - 1) + X_sl) &
            * (bottom - mm%d_g(bottom_g - 1)) / 2.0)
    end if
  end subroutine sum_g

end module surface_forcing




