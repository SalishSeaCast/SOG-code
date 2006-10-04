! $Id$
! $Source$

module initial_sog
  ! *** What's in here?
  implicit none

  ! *** These parameter values can probably be set in other, more
  ! *** appropriate modules
  DOUBLE PRECISION, PARAMETER:: &
       Uo = 0.0,     & ! m/s
       Vo = 0.0,     &   ! m/s
       hm =  2.0, & !28, & !75.0, & !Large1996
!!$       P_micro = 0.3D-3, &
! *** Parameter value setting of P_nano replaced by a variable version in 
! *** initial_mean below, so that initial value of flagellates biomass may
! *** be set to zero without recompiling
!       P_nano = 2.6D-3 * 0., & !V.flagella.01 add comm. 3.6D-3, &!2.6D-3 , & !7.5D-04 gN/m^3, winter estimate
       NHo = .5D-3 

contains

  subroutine initial_mean (Ui, Vi, Tnew, Snew, Pi, NO, NH, Si_new, Detritus, &
       hi, ut, vt, &
       pbx, pby, d, D_bins, cruise_id)       
    ! *** What's it do?
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use input_processor, only: getpars
    use mean_param, only: prop, plankton, snow
    implicit none
    ! Arguments:
    type(prop), intent(out) :: Ui, Vi
    real(kind=dp), dimension(0:), intent(out) :: &
         Tnew, &  ! Temperature profile
         Snew, &  ! Salinity profile
         Si_new   ! Silicon profile
    type(plankton), intent(out) :: Pi 
    real(kind=dp), dimension (0:), intent(out) :: NO  ! N%O%new,  nitrate
    real(kind=dp), dimension (0:), intent(out) :: NH  ! N%H%new,  ammonium
    integer, intent(in) :: D_bins
    type(snow), dimension(D_bins), intent(inout) :: Detritus
    real(kind=dp), intent(out) :: hi !h%new
    type(prop), intent(out) :: ut, vt
    type(grid_), intent(in) :: d
    real(kind=dp), dimension(d%M), intent(out) :: pbx, pby
    character*4  cruise_id           ! cruise_id

    ! Local variables:
    integer :: i, j      ! loop index
    ! File name to open
    character*80 :: fn           
    ! Place holders for reading CTD data file
    integer :: dum1
    real :: depth, dumc, dumt, dump, dumo
    ! *** P_nano variable replaces parameter version above,
    ! *** so that initial value of flagellates biomass may
    ! *** be set to zero without recompiling
    real(kind=dp) :: P_nano

       !V.flagella.01 add comm. 3.6D-3, &!2.6D-3 , & !7.5D-04 gN/m^3, winter estimate
       P_nano = 2.6D-3

    Ui%new(1) = Uo
    Vi%new(1) = Vo
!!$    Pi%micro%new(1) = P_micro
    Pi%nano%new(1) = P_nano !V.flagella.01
    NH(1) = NHo

    !-----detritus loop added march 2006---------------------------------

    open(unit=44, file="input/initial_Detritus.dat", &
         status="OLD", action="READ")
    do i = 1, d%M + 1
       read(44, *) Detritus(1)%D%new(i), Detritus(2)%D%new(i)
    end do
    close(44)

    ! *** The values read here are all overwritten, but removing this code
    ! *** causes unclean diffs (profiles are visually the same though)
    open(unit=49, file="input/NH4.dat", status="OLD", &
         action="READ")
    do i = 1, d%M + 1
       read(49, *) NH(i)
    end do
    close(49)


    !---biology jan 10 2005---------------------------------

    DO i = 2, d%M+1   
       IF (d%d_g(i) <= hm) THEN  !Large1996  March 1960 initial profile        

!!$          Pi%micro%new(i) = P_micro
          Pi%nano%new(i) = P_nano !V.flagella.01
          NH(i) = NHo
       ELSE
!!$          Pi%micro%new(i) = 0.
          Pi%nano%new(i) = 0. !V.flagella.01
          NH(i) = 0.
       END IF
       Ui%new(i) = Uo
       Vi%new(i) = Vo

    END DO

!!$    Pi%micro%new(d%M+1) = 0.
    Pi%nano%new(d%M+1) = 0.
    NH(d%M+1) = 0.
    Detritus(1)%D%new(d%M+1) = 0. ! (DON ==> Detritus(1), need some deep ocean value)
    Detritus(2)%D%new(d%M+1) =  0. !Detritus(2)%D%new(d%M) !PON needs a deep ocean value


    ! read in nutrients data
    write (fn,'("../sog-initial/Nuts_",a4,".txt")') cruise_id
    open(unit=66, file=fn, status="OLD", action="READ")

    do i = 0, d%M
       read(66, *) depth, NO(i), Si_new(i)
       if (depth.ne.i*0.5) then
          write (*,*) 'Expecting nutrients, NO3 and Silicon at 0.5 m intervals'
          stop
       endif
    enddo
    close (66)
    NO(d%M+1) = NO(d%M)
    Si_new(d%M+1) = Si_new(d%M)
    


    !---end biology------------------------------------------ 

    ! STRATOGEM CTD data from station S3 to initialize temperature,
    ! phytoplacnkton biomass, and salinity in water column
    fn = getpars("ctd_in")
    open(unit=46, file=fn, status="old")
    ! *** Maybe rework this so we can read orginal (not stripped) CTD
    ! *** data files?
    do i = 1, d%M + 1
       read(46, *) dum1, depth,Tnew(i), dumc, Pi%micro%new(i), &
            dumt, dump, dumo, Snew(i)  
       ! *** Maybe we should have a degC2degK function? 
       Tnew(i) = Tnew(i) + 273.15
       ! *** Does the next line actually do anything?
!!$       Pi%micro%new(i) = Pi%micro%new(i)
    enddo
    close(46)

    Tnew(0) = Tnew(1)  !Surface
    Snew(0) = Snew(1)  !Boundary
    Pi%micro%new(0) = Pi%micro%new(1)
    Pi%nano%new(0) = Pi%nano%new(1) !V.flagella.02
    NH(0) = NH(1)

    ! assuming dz = 0.5
    ! Interpolate CTD data values for grid points between measurements?
    do i = d%M + 1, 2, -1
       j = i / 2
       if (j * 2 == i) then
          Tnew(i) = Tnew(j)
          Snew(i) = Snew(j)
          Pi%micro%new(i) = Pi%micro%new(j) 
!!$          Pi%nano%new(i) = Pi%nano%new(j) !V.flagella.02
       else
          ! *** This looks like a job for a arith_mean function
          Tnew(i) = Tnew(j) * 0.5 + Tnew(j+1) * 0.5
          Snew(i) = Snew(j) * 0.5 + Snew(j+1) * 0.5
          Pi%micro%new(i) = Pi%micro%new(j) * 0.5 + Pi%micro%new(j+1) * 0.5
!!$          Pi%nano%new(i) = Pi%nano%new(j) * 0.5 + Pi%nano%new(j+1) * 0.5 !V.flagella.02
       endif
    enddo

    hi = hm  !Large1996 !ho 

    !If the bottom fluxes are fixed, use the following 
    !tp reevaluate M+1 values:
    Vi%new(0) = Vi%new(1)  !Conditions
    Ui%new(0) = Ui%new(1)

    !            initialize ut and vt
    do i=1,d%M
       ut%new(i) = 0.
       vt%new(i) = 0.
       pbx(i) = 0.
       pby(i) = 0.
    enddo
  end subroutine initial_mean

end module initial_sog
