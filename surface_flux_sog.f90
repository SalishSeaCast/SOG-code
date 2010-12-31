SUBROUTINE surface_flux_sog(mm,ro, & 
                         temp_o, I, Q_t, alp, Cp_o, &
                         bet, U_ten, V_ten, cf, atemp, humid)
  
  use fundamental_constants, only: g
  use turbulence, only: wbar
      USE mean_param
  implicit none
  ! Arguments:
      INTEGER, INTENT(IN)::mm
      DOUBLE PRECISION,DIMENSION(0:mm+1), INTENT(IN)::ro ! density%new
      DOUBLE PRECISION, INTENT(IN):: alp, Cp_o, bet !alph%i(0), Cp%i(0), beta%i(0)
      DOUBLE PRECISION, INTENT(IN)::&
           temp_o, & 
           U_ten, V_ten
      REAL, INTENT(IN):: cf,atemp,humid
      DOUBLE PRECISION,DIMENSION(0:mm),INTENT(IN)::I
      DOUBLE PRECISION, INTENT(OUT)::Q_t
      
      ! Local Variables:
      DOUBLE PRECISION:: UU, rho_atm
      double precision:: r, Ce, sigma, lw_in, lw_out, lw_net
      double precision:: Cs, h_sens, h_latent, h_flux, Cp
      double precision :: epsilon_w,a,b,c,ea,es,cl,le

!           U_ten, V_ten used to find UU for drag coefficient calculation
!           U_ten, V_ten used to find wind stress (with drag coefficient)

!           fresh water flux is F_tot

!           wbar%t is Q_t/(ro*cp)
!           wbar%s is Ft*sal/rf
!           wbar%b is g,alp,bet* wbar%t and wbar%s 


         UU = SQRT(U_ten**2.0 + V_ten**2.0)   !note U_ten and V_ten at 22m height


!net longwave radiation.
!-----------------------------------------------
r=0.03
Ce=9.37e-6
sigma=5.6697e-8  !stefan boltzmann constant
epsilon_w=0.96   !surface emissivity in the IR portion of the spectrum
! cloud fraction cf 
! air temperature atemp 
! water temperature temp_o 

lw_in=(1-r)*(1+0.17*cf**2)*Ce*atemp**2*sigma*atemp**4  
!downward radiation from atmosphere
!there are several different ways to calculate this

lw_out=-epsilon_w*sigma*temp_o**4                      
!upward emission of radiation from earth's surface, stull page 48

lw_net=lw_in + lw_out


!sensible heat flux
!-----------------------------------------------
Cs=1.3e-3 ! sensible heat transfer coefficient
!fix rho_atm at 1.25 kg/m^3
rho_atm = 1.25
Cp=1003 ! specific heat of air
! use UU (was U=metday(ii,8)) ! is this wind speed at the surface?

h_sens = Cs*rho_atm*Cp*UU*(atemp-temp_o) ! power

!latent heat flux
!-----------------------------------------------
!vapour pressure
a = 7.5
b = 237.3
c = 0.7858
ea = (humid/100)*exp(2.303*((a*(atemp-273.15)/((atemp-273.15) + b))+c));  ! mb
!saturated vapour pressure
es = exp(2.3026*((7.5*(temp_o-273.15)/((temp_o-273.15)+237.3))+0.7858))
!latent heat flux
CL=1.3e-3
LE=2.453e6
h_latent=(0.622/1013)*CL*rho_atm*LE*UU*(ea-es) ! latent heat
if (h_latent.lt.0) then
   h_latent=0
endif

!----Total heat flux (Q_t)
!-------------------------
      
h_flux = lw_net+h_sens+h_latent    

      Q_t = h_flux ! W/m2

      
!-----Surface fluxes---------------------  equations A2a-d
      ! Calculate surface flux components (see Large, etal (1994),
      ! eq'ns A2 & A3b)
      !
      ! Temperature (eq'n A2c)
      wbar%t(0) = -Q_t / (ro(0) * Cp_o)
END SUBROUTINE surface_flux_sog
