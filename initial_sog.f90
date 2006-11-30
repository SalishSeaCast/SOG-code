! $Id$
! $Source$

module initial_sog
  ! *** What's in here?
  implicit none

  ! *** These parameter values can probably be set in other, more
  ! *** appropriate modules
  DOUBLE PRECISION, PARAMETER:: &
       Uo = 0.0d0,     &  ! m/s
       Vo = 0.0d0,     &  ! m/s
!!$       P_micro = 0.3D-3, &
! *** Parameter value setting of P_nano replaced by a variable version in 
! *** initial_mean below, so that initial value of flagellates biomass may
! *** be set to zero without recompiling
!       P_nano = 2.6D-3 * 0., & !V.flagella.01 add comm. 3.6D-3, &!2.6D-3 , & !7.5D-04 gN/m^3, winter estimate
       NHo = 0.5d-3 

contains

  subroutine initial_mean(U_new, V_new, T_new, S_new, Pmicro, Pnano, Z, &
       NO, NH, Si, D_DON, D_PON, D_refr, D_bSi, &
       hml, d, cruise_id)       
    ! *** What's it do?
    use precision_defs, only: dp
    use grid_mod, only: grid_
    use input_processor, only: getpars
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(out) :: &
         U_new,  &  ! Cross-strait velocity component profile
         V_new,  &  ! Along-strait velocity component profile
         T_new,  &  ! Temperature profile
         S_new,  &  ! Salinity profile
         Pmicro, &  ! Micro phytoplankton profile
         Pnano,  &  ! Nano phytoplankton profile
         Z,      &  ! Microzooplankton profile
         NO,     &  ! Nitrate profile
         NH,     &  ! Ammonium profile
         Si,     &  ! Silicon profile
         D_DON,  &  ! Dissolved organic nitrogen detritus profile
         D_PON,  &  ! Particulate organic nitrogen detritus profile
         D_refr, &  ! Refractory nitrogen detritus profile
         D_bSi      ! Biogenic silicon detritus profile
    real(kind=dp), intent(in) :: hml ! Mixing layer depth
    type(grid_), intent(in) :: d
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
       P_nano = 2.6d-3

    U_new(1) = Uo
    V_new(1) = Vo
!!$    Pmicro(1) = P_micro
!!$    Pnano(1) = P_nano
    NH(1) = NHo

    !-----detritus loop added march 2006---------------------------------

    open(unit=44, file="input/initial_Detritus.dat", &
         status="OLD", action="READ")
    do i = 1, d%M + 1
       read(44, *) D_DON(i), D_PON(i)
    end do
    close(44)
    D_bSi = D_PON
    D_refr = 0.
    ! Need deep ocean values for bottom boundary condition
    D_DON(d%M + 1) = 0.
    D_PON(d%M + 1) = 0.
    D_bSi(d%M + 1) = 0.

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
       IF (d%d_g(i) <= hml) THEN  !Large1996  March 1960 initial profile        

!!$          Pmicro(i) = P_micro
!!$          Pnano(i) = P_nano
          NH(i) = NHo
       ELSE
!!$          Pmicro(i) = 0.
!!$          Pnano(i) = 0.
          NH(i) = 0.
       END IF
       U_new(i) = Uo
       V_new(i) = Vo

    END DO

!!$    Pmicro(d%M+1) = 0.
!!$    Pnano(d%M+1) = 0.
    NH(d%M+1) = 0.


    ! read in nutrients data
    write (fn,'("../sog-initial/Nuts_",a4,".txt")') cruise_id
    open(unit=66, file=fn, status="OLD", action="READ")

    do i = 0, d%M
       read(66, *) depth, NO(i), Si(i)
       if (depth.ne.i*0.5) then
          write (*,*) 'Expecting nutrients, NO3 and Silicon at 0.5 m intervals'
          call exit(1)
       endif
    enddo
    close (66)
    NO(d%M+1) = NO(d%M)
    Si(d%M+1) = Si(d%M)
    


    !---end biology------------------------------------------ 

    ! STRATOGEM CTD data from station S3 to initialize temperature,
    ! phytoplacnkton biomass, and salinity in water column
    fn = getpars("ctd_in")
    open(unit=46, file=fn, status="old")
    ! *** Maybe rework this so we can read orginal (not stripped) CTD
    ! *** data files?
    do i = 1, d%M + 1
       read(46, *) dum1, depth, T_new(i), dumc, Pmicro(i), &
            dumt, dump, dumo, S_new(i)  
       ! *** Maybe we should have a degC2degK function? 
       T_new(i) = T_new(i) + 273.15
       Pmicro(i) = 5.0/6.0*Pmicro(i)
       Pnano(i) = 1.0/5.0*Pmicro(i)
       Z(i) = 0.5*Pmicro(i)
    enddo
    close(46)

    T_new(0) = T_new(1)  !Surface
    S_new(0) = S_new(1)  !Boundary
    Pmicro(0) = Pmicro(1)
    Pnano(0) = Pnano(1)
    Z(0) = Z(1)
    NH(0) = NH(1)

    ! assuming dz = 0.5
    ! Interpolate CTD data values for grid points between measurements?
    do i = d%M + 1, 2, -1
       j = i / 2
       if (j * 2 == i) then
          T_new(i) = T_new(j)
          S_new(i) = S_new(j)
          Pmicro(i) = Pmicro(j) 
          Pnano(i) = Pnano(j)
          Z(i) = Z(j)
       else
          ! *** This looks like a job for a arith_mean function
          T_new(i) = T_new(j) * 0.5 + T_new(j+1) * 0.5
          S_new(i) = S_new(j) * 0.5 + S_new(j+1) * 0.5
          Pmicro(i) = Pmicro(j) * 0.5 + Pmicro(j+1) * 0.5
          Pnano(i) = Pnano(j) * 0.5 + Pnano(j+1) * 0.5
          Z(i) = Z(j) * 0.5 + Z(j+1) * 0.5
       endif
    enddo

    !If the bottom fluxes are fixed, use the following 
    !tp reevaluate M+1 values:
    V_new(0) = V_new(1)  !Conditions
    U_new(0) = U_new(1)
  end subroutine initial_mean

end module initial_sog
