! $Id$
! $Source$

! *** This is a collection of subroutines that should maybe be broken out
! *** into separate files ???

subroutine read_sog (upwell_const)
  ! Read in data from various files to initialize the run.
  ! Forcing data: 
  !               salinity at bottom of model domain
  ! *** And some other stuff that I need to figure out.
  ! *** Need to get rid of explicit file name references here.

  use precision_defs, only: dp
  use input_processor, only: getpard
  use declarations

  implicit none

  ! arguments
  real(kind=dp), intent(out) :: upwell_const 
                               ! tuned parameter for the strength of upwelling 

  ! Local variables:
  integer :: ndays, ic, j, para, stn, yr

  ! read the upwelling constant
  upwell_const = getpard("upwell_const")

  ! Preserve the value of year_o so the actual year does not change 
  ! to the last data year read when the data is read 
  ! ***huh?
  yr=year_o

  ! Read bottom salinity condition data
  ! *** Is this bottom of the model domain, or bottom of the Strait?
  open(unit=16, file="input/CTD/bottom_200123456.dat", &
       status = "OLD", action = "READ")
  ! *** Another hard-coded constant to det rid of **
  do xx = 1, 1659
     read(16, *) ctd_bottom(xx)%sal, ctd_bottom(xx)%temp, ctd_bottom(xx)%P, &
          ctd_bottom(xx)%No, ctd_bottom(xx)%date
  end do
  close(16)


  ! Read detritus model parameters
  ! *** Name of the detritus parameters file should be moved to the 
  ! *** run parameters file **
  open(unit=43, file="input/Detritus.dat", &
       status="OLD", action="READ")
  waste%m%destiny = 0.
  do xx = 1, D_bins - 1
     read(43, *) waste%m%destiny(xx), Detritus(xx)%r, Detritus(xx)%v
  end do
  waste%m%destiny(0) = 1. - waste%m%destiny(1) - waste%m%destiny(2)

end subroutine read_sog

