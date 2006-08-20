! $Id$
! $Source$

module mean_param
  ! This should probably be renamed to something more descriptive like
  ! type_defs.
  ! It would also be good to re-organize somehow so that they types
  ! hierarchy is more evident.

  implicit none

  TYPE :: gr_d              !Grid and  interface points and spacing
     DOUBLE PRECISION, DIMENSION(:), POINTER :: d_i, d_g, g_space, i_space
     INTEGER :: M, D
  END TYPE gr_d

  TYPE :: alpha             !Interpolation data
     DOUBLE PRECISION :: value, TT, SS, dS
  END TYPE alpha

  TYPE :: wind_ecmwf
     INTEGER :: day, year, month, Jday, leap
     DOUBLE PRECISION :: time, numb, zonal, meridional
  END TYPE wind_ecmwf

  TYPE :: insol_daily
     INTEGER :: year, month, day, Jday, leap
     DOUBLE PRECISION :: actual, clear_sky, ratio
  END TYPE insol_daily

  TYPE :: Large1996_data
     DOUBLE PRECISION :: day, data
  END TYPE Large1996_data

  TYPE :: bins
     sequence
     INTEGER :: micro, nano, NO, NH, det, Quant
  END TYPE bins

  ! expansion coefficients, heat capacity, surface density
  TYPE :: constant          
     TYPE(alpha), DIMENSION(13) :: data
     DOUBLE PRECISION, DIMENSION(:), POINTER:: i, g, idiv, gdiv
  END TYPE constant

  TYPE :: grow               
     DOUBLE PRECISION, DIMENSION(:), POINTER :: light, new, net !**&    
  END TYPE grow

  TYPE :: grazing
     DOUBLE PRECISION, DIMENSION(:), POINTER :: new   !**&
  END TYPE grazing


  TYPE :: prop              !Properties (U, V, S...
     DOUBLE PRECISION,DIMENSION(:), POINTER :: new, old, old_old, div_i, div_g, last
     DOUBLE PRECISION :: avg      !surface layer average 
  END TYPE prop

  TYPE :: plankton                  
     TYPE(prop) :: micro, nano  !%new, old, old_old (0:M+1)
  END TYPE plankton

  TYPE :: zplankton
     TYPE(prop) :: micro, meso  !**&
  END TYPE zplankton

  !       TYPE :: Nconstants
  !        !  DOUBLE PRECISION::Vm, K, Ki
  !       END TYPE Nconstants

  TYPE :: plankton2                   
     ! TYPE(Nconstants)::NO,NH
     TYPE(grow)::growth
     TYPE(grazing):: mort   !**&
!!!phytoplankton constants!!!!!!!!!!!!!!!!!!
     DOUBLE PRECISION :: sink, R, sigma, &  !sink = sinking vel, R = max growth rate for light
          gamma, Rm, inhib, Rmax    ! sigma, attenuation coeff for light growth
     DOUBLE PRECISION :: k, kapa, gamma_o, N_o, N_x  !@@@
     DOUBLE PRECISION :: Q_cn, Q_old, Q_old_old, dlnQ_dt 
     !CN ratio (old copies) and 1/Q*dQ/dt
     !gamma = photorespiration param , Rm = maintenance respiration 
!!!zooplankton constants!!!!!!!!!!!!!!!!!
     DOUBLE PRECISION :: delta, eta, G, ks, nn !delta = assim efficiency, eta = excretion
     !G = max grazing rate, lambda = slope of linear grazing at small [prey]
     !ks and nn holling type nn+1 grazing function param.
     DOUBLE PRECISION :: M_z, beta     !M_z = specific mortality, and beta = mortality time constant
     DOUBLE PRECISION, DIMENSION(:), POINTER::q  !preference q(1+d_prey)
     DOUBLE PRECISION, DIMENSION(:,:),POINTER::graze !graze(1+d_prey,mm%M)
     !used in p_growth.dat
  END TYPE plankton2

  TYPE :: phyto                     
     DOUBLE PRECISION, DIMENSION(:), POINTER::micro, nano
     DOUBLE PRECISION, DIMENSION(1)::micro_q,nano_q
  END TYPE phyto

  TYPE :: macro_zoo
     DOUBLE PRECISION, DIMENSION(:), POINTER::macro, wt !macro(M), wt(bin)
  END TYPE macro_zoo

  TYPE :: zooplank
     DOUBLE PRECISION, DIMENSION(:), POINTER::micro
     TYPE(macro_zoo), DIMENSION(:), POINTER::c
  END TYPE zooplank

  TYPE :: nut
     DOUBLE PRECISION, DIMENSION(:), POINTER::o, h
  END TYPE nut

  TYPE :: det
     DOUBLE PRECISION, DIMENSION(:), POINTER::bin  !bin(M)
  END TYPE det

  TYPE :: nutrient
     TYPE(prop)::O,H
     TYPE(grazing)::O_uptake, H_uptake, urea
     DOUBLE PRECISION, DIMENSION(:), POINTER :: remin,bacteria  
     DOUBLE PRECISION::r 
  END TYPE nutrient

  TYPE :: diff
     DOUBLE PRECISION, DIMENSION(:), POINTER :: shear, total 
     DOUBLE PRECISION, DIMENSION(:), POINTER :: dd    !double diffusion
     DOUBLE PRECISION, DIMENSION(:), POINTER :: ML    !mixed layer
     DOUBLE PRECISION, DIMENSION(:), POINTER :: all, old   !total and ML combined
     DOUBLE PRECISION :: div      !only for interior
     DOUBLE PRECISION :: h        !diffusivities at h
  END TYPE diff

  TYPE :: Knu               !K and nu diffusivities
     TYPE(diff)::u, t, s       
  END TYPE Knu

  TYPE :: flux              !Reynolds fluxes
     DOUBLE PRECISION, DIMENSION(:), POINTER::u, v,s,t,b,b_old, &
          b_err,b_err_old            
     TYPE(phyto)::p                                         
  END TYPE flux

  TYPE :: entrain           !Entrainment or mixed layer depth
     DOUBLE PRECISION :: depth
     INTEGER :: i,g
  END TYPE entrain

  TYPE :: height            !boundary layer depth 
     DOUBLE PRECISION :: old, old_old, new
     INTEGER :: i, g        !interface and grid index
     TYPE(entrain)::e, ml   !entrainment depth or mixed layer depth
  END TYPE height

  TYPE :: MST               ! Momentum, Salinity and Temp component vectors
     DOUBLE PRECISION, DIMENSION(:), POINTER::m,s,t    
  END TYPE MST

  TYPE :: UVST              
     DOUBLE PRECISION, DIMENSION(:), POINTER::u,v,s,t
     TYPE(phyto)::p 
     TYPE(zooplank)::z
     TYPE(nut)::n
     TYPE(det), DIMENSION(:), POINTER::d !d(D_bins)
  END TYPE UVST

  TYPE :: trivector
     DOUBLE PRECISION, DIMENSION(:), POINTER::A,B,C
  END TYPE trivector

  TYPE :: UVSTmatrix
     TYPE(trivector)::u,s,t,bio,no,null,null2
     DOUBLE PRECISION, DIMENSION(1)::QA,QB
  END TYPE UVSTmatrix

  TYPE :: MSTscalar
     DOUBLE PRECISION::m,s,t
  END TYPE MSTscalar

  TYPE :: boundary          !boundary conditions at h
     DOUBLE PRECISION :: h, div
     DOUBLE PRECISION, DIMENSION(:), POINTER::value
  END TYPE boundary

  TYPE :: MS                !Momemtum and scalar
     TYPE(boundary)::m, s
  END TYPE MS

  TYPE :: cloudy
     DOUBLE PRECISION :: A, B, fraction   !Regression coefficients
  END TYPE cloudy

  TYPE :: okta
     TYPE(cloudy), DIMENSION(0:9) :: type
  END TYPE okta

  TYPE :: light
     !     DOUBLE PRECISION, DIMENSION (365) :: R, mu_1, mu_2
     DOUBLE PRECISION, DIMENSION (82) :: type1,type3,type5,type7,type9
  END TYPE light

  !      TYPE :: Jerlov
  !        TYPE(light), DIMENSION(82) :: type
  !      END TYPE Jerlov

  TYPE :: old_new
     DOUBLE PRECISION::old,new
  END TYPE old_new

  TYPE :: windstress
     TYPE(old_new)::u,v
  END TYPE windstress

  !stores data from Large's data file PAPMD.60
  TYPE :: papmd  ! Large_data
     INTEGER::ymdh, y, jday,leap !year-month-day-hour, year, jday, leap = 1 (yes) 
     !or 0 (no) if leap year 
     DOUBLE PRECISION::Uten, &! 10 m velocity (cm/s) ==> m/s
          theta, &  !direction (wind blowing from) measured east from north
          SST, & !sea surface temperature (oC) ==> K
          Ta, &!air temperature (oC) at 17 m ==> K
          Tw, & !wet bulb temperature (oC) ==> K
          Td, & !dew point temp (oC) ==> K
          cf, & !cloud fraction
          P, & !air pressure (mbar)
          t !time (s)
  END TYPE papmd


  !Copepod Types:

  TYPE :: event                      !Cevent(Csources)
     DOUBLE PRECISION::var, norm, var_67   !may not need norm!!!!  used in year 1967
     DOUBLE PRECISION, DIMENSION(:), POINTER::s, stage_time   !!both allocated
     DOUBLE PRECISION::nauplii,nauplii_67    !total number of arrival nauplii , used in year 1967           
     INTEGER::Cday, year_s, year_e, type, length, start, end, on, n, change_67,Cday_67
  END TYPE event

  !TYPE :: wt_dist
  !   DOUBLE PRECISION::wt, wt_old, wt_old_old, f  !weight (and old copies) and fraction
  !END TYPE wt_dist

  TYPE :: mat_dist
     DOUBLE PRECISION::M,sigma_M, f, o_f, pdf_out, o_pdf_out, & 
          !<M>, (<M^2>-<M>^2)^1/2, number of animals
     wt_node, wt_avg, o_wt_avg, o_wt_node, day_wt_avg, mwt_avg, M_old, &
          b, wt_m, norm, cnout  !b = condition number, wt_m = migration weight
     !in Mixed layer for a given event divided by current Ntot, 
     !pdf_out=frac in stage 6
     !o_pdf_out = previous days frac in stage 6 (for each initial
     !distribution === used to normalize pdf
     INTEGER::day  !arrival day
     DOUBLE PRECISION, DIMENSION(:), POINTER::stage_dur  !stage #-1 stage_dur(Zoo%s_number-1) allocated
  END TYPE mat_dist

  TYPE :: copepod                     !species(number of copepod events == Csources)
     !use Cevent to find copepod type. may be same for each species
     DOUBLE PRECISION, DIMENSION(:), POINTER::stage_f, o_stage_f,& !allocated  Total fraction in stage(yy)
          avg_stage_wt,o_avg_stage_wt !stage(zoo(Cevent(yy)%type)%s_number-1)
     DOUBLE PRECISION::n_out, wt_out !number of Copepods per m^2 in diapause and 
     DOUBLE PRECISION::new_wt_out, new_wt_in  !use to find daily migration flux (in or out) of ml 
     !total wt per m^2 in diapause(number*size)
     TYPE(prop)::Z   !species%n*species%Ntot = species%Z%new
     DOUBLE PRECISION, DIMENSION(:), POINTER::n !, PZ  
     !n(0:M+1) number density  !allocated  see define_PZ for precise definition
     !x(Zoo%s_number) ~ day-mature_pdf%day and y(Zoo%s_number) ~ Maturity
     ! interpolation points
     !for mature_pdf%M  !allocated
     !stores PZ(M) version of species%Z%new  !allocated
     DOUBLE PRECISION, DIMENSION(:), POINTER::x,y
     DOUBLE PRECISION::Ntot,gamma, avg_wt, node_wt, gamma_wt, a, b !Avg number of animals per m^2 
     !see define_PZ for precise definition
     !gamma*(t-t_o) = sigma for stage_pdf
     !gamma_wt*(wt(t+1)-wt(t)) = var for each bin of wt_pdf
     !a and b are best fit line coefficients for M = a*log(wt)+b
     INTEGER::current_arrival  !number from 1--cvent%length+1, If cevent%length + 1 ==> all have arrived
     !TYPE(wt_dist), DIMENSION(:), POINTER::wt_pdf !wt_pdf(bin)   
     TYPE(mat_dist), DIMENSION(:), POINTER::mature_pdf  !mature_pdf(Cevent%length)  !allocated
     DOUBLE PRECISION, DIMENSION(:), POINTER::ingest !ingest(Cevent%length)  !allocated
     DOUBLE PRECISION, DIMENSION(:,:), POINTER::graze !graze(prey,M) !allocated
     DOUBLE PRECISION, DIMENSION(:), POINTER::Ex  !excrete(Cevent%length) !allocated
     DOUBLE PRECISION, DIMENSION(:,:), POINTER::p !preferences(prey,M)  !allocated
     DOUBLE PRECISION, DIMENSION(:), POINTER::mort !mortality(M) m^3/#/s
  END TYPE copepod

  TYPE :: Cdata                  !Zoo(number of copepod types == C_types)
     INTEGER::s_number        !number of stages  
     DOUBLE PRECISION, DIMENSION(:),POINTER::stage_dur   !dimension: stage #-1  !allocated
     DOUBLE PRECISION, DIMENSION(:),POINTER::molt_wt !dimension: stage #  molt_wt(0:stage#-1)
     ! DOUBLE PRECISION::wt, min_wt  !weight of arrival nauplii/first stage and smallest wt
     DOUBLE PRECISION, DIMENSION(:),POINTER::q  !preference_coefficient(prey)  !allocated
     DOUBLE PRECISION::delta, eta, a_Ex, b_Ex, Gmax, ks, M, nn, b2_Ex ! delta = assimilation eff ;
     !eta, a_Ex and b_Ex = excretion coefficients
     !Gmax = max grazing rate; ks = in holling grazing
     !M = mortality coefficient
     !nn = holling type nn+1 exponent
     DOUBLE PRECISION::B,k  !M = LOG((1.+species%B)/(1.+species%B*EXP(-species%k*day)))
  END TYPE Cdata

  TYPE :: write_bio   !total
     DOUBLE PRECISION::nano,diatom,zmicro,copepod,out,NO,NH,size,PN,NPP,n_loss,ngrow,dgrow,DN,&
          zgrazen,cgrazed,cgrazez,cgrazed2
  END TYPE write_bio

  !Detritus Types

  TYPE :: snow  !Detritus(D_bins)
     DOUBLE PRECISION::r, & !regeneration rate
          v ! sinking velocity
     TYPE(prop)::D !%new, old, old_old (0:M+1)
  END TYPE snow

  TYPE :: loss_param
     DOUBLE PRECISION, DIMENSION(:), POINTER::destiny !waste%s%destiny(0:D_bins): 0==>NH, 
     !1:D_bins ==>Detritus 
  END TYPE loss_param

  TYPE :: losses  !waste
     DOUBLE PRECISION, DIMENSION(:),POINTER::small, medium, large  !new waste pools (M)
     TYPE(loss_param)::s,m,l
  END TYPE losses

  INTEGER::is_leap_year, was_leap_year

  TYPE :: bottom_fit
     DOUBLE PRECISION:: temp, sal, P, No, date
  END TYPE bottom_fit

CONTAINS


  SUBROUTINE div_grid (dm, X)

    TYPE(gr_d), INTENT(IN)::dm 
    TYPE(prop), INTENT(IN OUT)::X
    INTEGER::i

    X%div_g(1) = (X%new(1)-X%new(2))/(2*dm%g_space(1))

    DO i = 2, dm%M
       X%div_g(i) = (X%new(i-1)-X%new(i+1))/(dm%g_space(i-1)+ &
            dm%g_space(i))
    END DO

  END SUBROUTINE div_grid

  SUBROUTINE div_g_param (mm, alp)

    TYPE(gr_d), INTENT(IN)::mm
    TYPE(constant), INTENT(IN OUT)::alp  
    INTEGER::ii

    DO ii = 1,mm%M
       alp%gdiv(ii) = (alp%i(ii-1) - alp%i(ii))/mm%i_space(ii)
    END DO

  END SUBROUTINE div_g_param

  SUBROUTINE div_i_param (gr, bet)

    TYPE(gr_d), INTENT(IN)::gr
    TYPE(constant), INTENT(IN OUT)::bet  
    INTEGER::jj

    DO jj = 1,gr%M
       bet%idiv(jj) = (bet%g(jj) - bet%g(jj+1))/gr%g_space(jj)
    END DO

  END SUBROUTINE div_i_param


  SUBROUTINE div_interface(dm1, X1)

    TYPE(gr_d), INTENT(IN)::dm1         
    TYPE(prop), INTENT(IN OUT)::X1    !U, V...
    INTEGER::j

    DO j = 1, dm1%M
       X1%div_i(j) = (X1%new(j)-X1%new(j+1))/dm1%g_space(j)
    END DO

  END SUBROUTINE div_interface



  SUBROUTINE find_jmax_g(hh, d)

    TYPE(height), INTENT(IN OUT)::hh !mixed layer depth
    TYPE(gr_d), INTENT(IN)::d 

    INTEGER::i


    DO i = 1, d%M
       IF (d%d_g(i) > hh%new) THEN
          hh%g = i
          EXIT
       ELSE
          hh%g= i+1
       END IF
    END DO

    IF (hh%g >= d%M-2)THEN
       PRINT "(A, I5)","Mixing too deep. find_jmax_g. hh%g>=grid&M-2", hh%g
       write (*,*) 'mean_param'
       STOP
    END IF

  END SUBROUTINE find_jmax_g

  SUBROUTINE find_jmax_i(h1, d1)

    TYPE(height), INTENT(IN OUT)::h1 !depth
    TYPE(gr_d), INTENT(IN)::d1

    INTEGER::i

    DO i = 1, d1%M
       IF (d1%d_i(i) >= h1%new) THEN
          h1%i = i
          EXIT
       ELSE
          h1%i = i+1
       END IF
    END DO

    IF (h1%i >= d1%M-2)THEN
       PRINT "(A, I5)","Mixing too deep. find_jmax_i:", h1%i
       STOP
    END IF

  END SUBROUTINE find_jmax_i

END MODULE mean_param















