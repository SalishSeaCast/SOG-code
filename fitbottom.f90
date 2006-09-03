module fitbottom

  private
  public :: init_fitbottom, bot_bound_time, bot_bound_uniform

  integer :: NQ ! number of different quantities
  parameter (NQ=6)
  character*12, dimension(NQ) :: quantity !*** how do I make this not 12?
  double precision, dimension(5,NQ) :: c ! fit coefficients

contains

  subroutine init_fitbottom

    implicit none

    ! data from fitting (see fitbottom.py and fitted.m
    
    ! Salinity (constant, seasonal and biseasonal components)
    quantity(1) = 'salinity'
    c(:,1) = [29.62605865, 0.10374454, -0.03562458, -0.14156091, -0.06348989]

    ! Temperature (constant, seasonal and biseasonal components)
    quantity(2) = 'temperature'
    c(:,2) = [ 9.34995751, -0.27245442, -0.90151197, 0.17933473, -0.05590939]
    
    ! Phytoplankton from fluor (constant and seasonal)
    quantity(3) = 'chl fluor'
    c(:,3) = [ 0.58316137, -0.11206845, 0.26241523, 0., 0.]
    
    ! Nitrate (constant)
    quantity(4) = 'nitrate'
    c(:,4) = [25.7263829787, 0., 0., 0., 0.]
    
    ! Silicon (constant)
    quantity(5) = 'silicon'
    c(:,5) = [46.5530851064, 0., 0., 0., 0.]
    
    ! Ratio of pico/micro plankton (constant and biseasonal)
    quantity(6) = 'plank ratio'
    c(:,6) = [ 1.25948868, 0., 0., 0.97697686, 0.46289294]
    
  end subroutine init_fitbottom

  double precision function bottom_value (arg, qty)

    implicit none

    double precision, intent(in):: arg
    character*(*), intent(in):: qty
    
    ! local variable
    integer:: index

    if (qty.eq.quantity(1)) then
       index = 1
    elseif (qty.eq.quantity(2)) then
       index = 2
    elseif (qty.eq.quantity(3)) then
       index = 3
    elseif (qty.eq.quantity(4)) then
       index = 4
    elseif (qty.eq.quantity(5)) then
       index = 5
    elseif (qty.eq.quantity(6)) then
       index = 6
    else
       write (*,*) 'Problems matching quantity for boundary conditions'
       write (*,*) 'See fitbottom.f90'
       write (*,*) '*',qty,'*',quantity(1),'*'
       stop
    endif

    bottom_value = c(1,index) + c(2,index)*cos(arg) + c(3,index)*sin(arg) &
         + c(4,index)*cos(2.*arg) + c(5,index)*sin(2.*arg)
    
  end function bottom_value

  subroutine bot_bound_time (day, day_time, &
       Tbot, Sbot, Nobot, Silbot, Pmbot, Pnbot)

! For those variables for which we have data, from the bottom value (M+1) from
! an annual fit.

    use surface_forcing, only: PI
    use unit_conversions, only: CtoK

    implicit none

    integer, intent(in):: day
    double precision, intent(in):: day_time
    double precision, intent(out):: Tbot, Sbot, Nobot, Silbot, Pmbot, Pnbot
    
    ! local variables
    double precision:: arg, chl, ratio
    
    arg = 2*PI*(day+day_time/86400.)/365.25
    
    Tbot = bottom_value(arg, 'temperature')          ! in Celcius
    Tbot = CtoK(Tbot)
    Sbot = bottom_value(arg, 'salinity')
    Nobot = bottom_value(arg, 'nitrate')
    Silbot = bottom_value(arg, 'silicon')
    
    chl = bottom_value(arg, 'chl fluor')
    ratio = bottom_value(arg, 'plank ratio')
    
    Pmbot = chl / (ratio + 1)
    Pnbot = ratio * Pmbot
    
  end subroutine bot_bound_time
  
  subroutine bot_bound_uniform (M, NH, Detritus)
    
    use mean_param, only: snow
    use declarations, only: D_bins
    
    implicit none
    
    integer:: M ! length of grid
    double precision, dimension(0:), intent(in out) :: NH
    type (snow), dimension(:), intent(in out) :: Detritus
    
    ! local variable
    integer xx
    
    NH(M+1) = NH(M)
    do xx = 1,D_bins
       Detritus(xx)%D%new(M+1) = Detritus(xx)%D%new(M)
    enddo
    
  end subroutine bot_bound_uniform
  
end module fitbottom

  
