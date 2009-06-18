! $Id$
! $Source$

module initial_sog
  ! Initialization of core variables including the reading of the CTD start
  ! file
  !
  ! Public Parameters:     
  !
  ! N2chl -- ratio of mg/m3 of chl to uM of N in phytoplankton
  !

  use precision_defs, only: dp
  use unit_conversions, only: CtoK
  use forcing, only: vary_forcing, vary

  implicit none
  private
  public :: &
       ! Parameter values:
       N2chl, &  ! ratio of chl mg/m3 to uMol N in phytoplankton
       ! Subroutines:
       initial_mean

  ! Public parameter declarations :
  real(kind=dp), parameter :: &
       N2chl = 1.5d0  ! ratio of chl mg/m3 to uMol N in phytoplankton


  ! *** These parameter values can probably be set in other, more
  ! *** appropriate modules
  DOUBLE PRECISION, PARAMETER:: &
       Uo = 0.0d0,     &  ! m/s
       Vo = 0.0d0,     &  ! m/s

       ! estimate of deep NH from Nitrate/Salinity fits.  Gives deep
       ! nitrate as 31.5 uM but we measure 30.5 --- rest must be
       ! remineralized NH
       NHo = 1.0d0

contains

  subroutine initial_mean(U_new, V_new, T_new, S_new, Pmicro, Pnano, Ppico, &
       Z, NO, NH, Si, D_DON, D_PON, D_refr, D_bSi, &
       hml, d, cruise_id)       
    ! *** What's it do?
    use precision_defs, only: dp
    use grid_mod, only: grid_,interp_array
    use input_processor, only: getpars, getpardv, getparl
    implicit none
    ! Arguments:
    real(kind=dp), dimension(0:), intent(out) :: &
         U_new,       &  ! Cross-strait velocity component profile
         V_new,       &  ! Along-strait velocity component profile
         T_new,       &  ! Temperature profile
         S_new,       &  ! Salinity profile
         Pmicro,      &  ! Micro phytoplankton profile
         Pnano,       &  ! Nano phytoplankton profile
         Ppico,       &  ! Pico phytoplankton profile
         Z,           &  ! Microzooplankton profile
         NO,          &  ! Nitrate profile
         NH,          &  ! Ammonium profile
         Si,          &  ! Silicon profile
         D_DON,       &  ! Dissolved organic nitrogen detritus profile
         D_PON,       &  ! Particulate organic nitrogen detritus profile
         D_refr,      &  ! Refractory nitrogen detritus profile
         D_bSi          ! Biogenic silicon detritus profile

     real(kind=dp), dimension(0:81):: &    
         input,        &  ! Place holder for incoming data for different variables
         
         depth,         &  ! Depth profile from ctd data file 
         depth_index,   &  ! Depth profile from ctd data file
         T_interp,      &  ! Temperature interpolted at every .5m
         depth_interp,  &  ! Depth at every .5m
         S_interp,      &  ! Interpolated Salinity profile
         Pmicro_interp, &  ! Interpolated Micro phyto profile
         dummyNO,       &  ! Dummy Nitrate valu
         NO_interp,     &  ! Interpolated NO profile
         Pnano_interp,  &  ! Interpolated Nano phyto profile
         Ppico_interp,  &  ! Interpolated Pico phyto profile
         Z_interp,      &  ! Interpolated micro zooplankton profile
         D_PON_interp,  &  ! Interpolated PON profile
         D_DON_interp,  &  ! Interpolated DON profile
         D_BSi_interp,  &  ! Interpolated biogenic silicon profile
         D_Refr_interp   ! Interpolated refractory nitrogen profile



    real(kind=dp), intent(in) :: hml ! Mixing layer depth
    type(grid_), intent(in) :: d
    character*4  cruise_id           ! cruise_id

    ! Local variables:
    integer :: i, j      ! loop index
    integer :: position(7), position_bot(7) !variable placeholders
    integer :: index, index1        ! data index
    integer :: length, count, mcol_ctd, mcol_bot, header_length
    real :: data(50,20),databot(50,20)   !Place holders for ctd/bot data
    ! Used to find bottom of 40m and whether ctd file has fluores data 
    logical :: found_depth, noFluores, Nitrate
    ! File name to open
    character*80 :: fn           
    ! fraction to split Micro/Nano/Pico into
    real (kind=dp) :: split(3)
    real :: depth1


    U_new(1) = Uo
    V_new(1) = Vo

    DO i = 2, d%M+1   
       IF (d%d_g(i) <= hml) THEN               
          NH(i) = 0.   ! essentially zero in mixed layer except in winter
       ELSE
          NH(i) = NHo   ! deep value
       END IF
       U_new(i) = Uo
       V_new(i) = Vo
    END DO

    NH(d%M+1) = NHo


    ! read in nutrients data

    ! Do we have ctd nitrate data, i.e every .5 m?
    Nitrate  = getparl('Nitrate')
    
    write (fn,'("../sog-initial/Nuts_",a4,".txt")') cruise_id
    open(unit=66, file=fn, status="OLD", action="READ")
    
    do i = 0, d%M
       read(66, *) depth1, dummyNO(i), Si(i)
       if (depth1.ne.i*0.5) then
          write (*,*) 'Expecting nutrients, NO3 and Silicon at 0.5 m intervals'
          call exit(1)
       endif
    enddo
    close (66)
   
    Si(d%M+1) = Si(d%M)
   if(Nitrate)then 
      NO = dummyNO
      NO(d%M+1) = NO(d%M)
   endif

  

    ! STRATOGEM CTD data from station S3 to initialize temperature,
    ! phytoplacnkton biomass, and salinity in water column

    
   

   fn = getpars("ctd_in")
    
    open(unit=46, file=fn, status="old")

    header_length = 11
    
    ! Read past header
    do i =1, header_length
       read(46,*)
   
    enddo

    ! Read in column number of each variable (temp,sal etc.)
    read(46,*)(position(j), j = 1,7)
    mcol_ctd =max(position(1),position(2),position(3),position(4),position(5),position(6),position(7))

    
    !Read in data of upper 40m or if data ends before 40m,will read entire file
    found_depth =.false.
    index = 0
    do while (.not. found_depth)
       read(46,*,END=176)(data(index+1,i),i=1,mcol_ctd)
      
       if(data(index+1,position(1)) > 40)then
      
          found_depth =.true.

       endif
       index = index + 1

    enddo

!Check to see if there is fluores data in ctd file

    
176 noFluores = .false.

    ! if there is not enough data to reach below 40m, set
    ! a depth at 50m with same values as the previous depth  

    if(.not. found_depth) then         
       
       data(index+1, 1) = 50

       do i = 2,mcol_ctd
          data(index+1, i) = data(index, i)
       enddo

    endif

    if (position(5)==-1)then

       noFluores = .true.
  

    else
       do j=1,index
          Pmicro(j) = data(j,position(5))

       enddo

       length = index +1  ! for calculating C2K later on

    endif

    ! Importing Nitrate data
    if(.not. Nitrate) then
       
       if (position(6)==-1)then
          
          Nitrate = .false.
       else

          do j=1,index+1

             NO(j) = data(j,position(6))

          enddo
          Nitrate = .true.
  
       endif
     
    endif
    

    !Importing temperature,sal and depth 
    do i=1,3

       if(position(i) == -1)then
             
          input = 0
             
       else
          do j=1,index+1
             input(j) = data(j,position(i))
          enddo
     

       endif
   
    SELECT CASE (INT(i))   
         CASE (1)                          
            depth = input
         CASE (2)                       
            T_New = input
         CASE (3)
            S_New = input
     
     END SELECT
 
   enddo
   close(46)


   fn = getpars("bot_in")

   if(noFluores .OR. .not. Nitrate)then
      !Open up the bottle file corresponding to the ctd file

      
      open(unit=45, file=fn , status="old")

      ! Read past header
      do i =1, header_length
         read(45,*)
  
      enddo

      ! Read in column number of each variable (NO,Chlr etc.)
      read(45,*)(position_bot(j), j = 1,7)

      mcol_bot =max(position_bot(1),position_bot(2),position_bot(3),position_bot(4),position_bot(5),position_bot(6),position_bot(7))

      !Read in data of upper 40m
      found_depth =.false.
      index1 = 0
      do while(.not. found_depth)
         read(45,*,END=177)(databot(index1+1,i),i=1,mcol_bot)
   
         if(databot(index1+1,position_bot(1)) > 40)then
    
            found_depth =.true.

         endif
         index1 = index1 + 1
      
      enddo

177   input = 0  ! clear out input from ctd file
       
      if(.not. found_depth) then     ! if there is not enough data to reach below 40m, set
                                     ! a depth at 50m with same values as the previous depth
       
         
         databot(index1+1, 1) = 50

         do i = 2,mcol_bot
            databot(index1+1, i) = databot(index1, i)
         enddo
         
         

      endif
    ! Importing Nitrate data if not in CTD file
    if(.not. Nitrate) then
       
       do j=1,index1+1

          NO(j) = databot(j,position_bot(6))

       enddo
      
    endif
      ! Check to see if theres fluores data. If not, then check to see if theres chl data. If not, Check to see if theres phyto data. If not, Pmicro = 0

      
 if(noFluores) then

    if(position_bot(5) == -1)then
         
       if(position_bot(4) == -1)then
            
          if(position_bot(7) == -1)then

             input = 0
                
          else
                
             do j=1,index1+1
                      input(j) = databot(j,position_bot(7))
                      
             enddo

          endif
             
       else

          do j=1,index1+1
                input(j) = databot(j,position_bot(4))
    
          enddo

       endif

     else

          do j=1,index1+1
             input(j) = databot(j,position_bot(5))
    
          enddo

     endif
  endif
       
        
          Pmicro = input
            

       
 
      
   
   endif



    ! split between micro and nano and pico plankton
    call getpardv("initial chl split", 3, split)
  
     
    do i = 1, length 
       ! Change temperature to Kelvin
       if (vary%temperature%enabled .and. .not.vary%temperature%fixed) then
          T_new(i) = CtoK(T_new(i)) + vary%temperature%addition
       else
          T_new(i) = CtoK(T_new(i))
       endif
       ! convert fluorescence to uMol N
       Pmicro(i) = Pmicro(i)/N2chl

       ! split up photoplankton
       Ppico(i) = split(3)*Pmicro(i)
       Pnano(i) = split(2)*Pmicro(i)
       Pmicro(i) = split(1)*Pmicro(i)

       Z(i) = Pnano(i) !*** hard value to estimate
       D_PON(i) = Pmicro(i)/5. ! estimate
       D_DON(i) = D_PON(i)/10. ! estimate
       D_Bsi(i) = D_PON(i)
       D_Refr(i) = 0. 
    enddo
    close(45)
    

!Defining depth_interp array: assuming dz = 0.5m 
    
    count=0
    do i=0,81,1
       depth_interp(i)=count/2
       count = count +1
    enddo

!Define depth index to keep track of original depths

    depth_index = depth

    !Linear interpolation to obtain variables at every .5m depth


    call interp_array(depth,T_New,depth_interp,T_interp)
    depth = depth_index
    call interp_array(depth,S_New,depth_interp,S_interp)
    depth = depth_index
    call interp_array(depth,Pmicro,depth_interp,Pmicro_interp)
    depth = depth_index    
    if(.not.Nitrate) then
       call interp_array(depth,NO,depth_interp,NO_interp)
       depth = depth_index
       NO = NO_interp
    endif
    call interp_array(depth,Pnano,depth_interp,Pnano_interp)
    depth = depth_index
    call interp_array(depth,Ppico,depth_interp,Ppico_interp)
    depth = depth_index
    call interp_array(depth,Z,depth_interp,Z_interp)
    depth = depth_index
    call interp_array(depth,D_PON,depth_interp,D_PON_interp)
    depth = depth_index
    call interp_array(depth,D_DON,depth_interp,D_DON_interp)
    depth = depth_index
    call interp_array(depth,D_Bsi,depth_interp,D_Bsi_interp)
    depth = depth_index
    call interp_array(depth,D_Refr,depth_interp,D_Refr_interp)
    depth = depth_index

    T_new  = T_interp  !Surface
    S_new = S_interp  !Boundary
    Pmicro = Pmicro_interp
    Pnano = Pnano_interp
    Ppico = Ppico_interp
    Z = Z_interp
    D_PON = D_PON_interp
    D_DON = D_DON_interp
    D_Bsi = D_BSi_interp
    D_Refr = D_Refr_interp


    !If the bottom fluxes are fixed, use the following 
    !tp reevaluate M+1 values:
    V_new(0) = V_new(1)  !Conditions
    U_new(0) = U_new(1)
   
  end subroutine initial_mean

end module initial_sog
