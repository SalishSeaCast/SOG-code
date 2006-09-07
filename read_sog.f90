! $Id$
! $Source$

! *** This is a collection of subroutines that should maybe be broken out
! *** into separate files ???

subroutine read_sog (upwell_const)
  ! Read in data from various files to initialize the run.
  ! Forcing data: wind, cloud fraction, air temperature, air humidity,
  !               Fraser River flow, Englishman River flow,
  !               salinity at bottom of model domain
  ! *** And some other stuff that I need to figure out.
  ! *** Need to get rid of explicit file name references here.

  use precision_defs, only: dp
  use declarations
  use surface_forcing

  implicit none

  ! arguments
  real(kind=dp), intent(out) :: upwell_const 
                               ! tuned parameter for the strength of upwelling 

  ! Local variables:
  integer :: ndays, ic, j, para, stn, yr

  ! *** Should be able to replace this external statment with uses of
  ! corresponding modules **
  external Julian_day  
  real(kind=dp) getpard

  ! read the upwelling constant
  upwell_const = getpard("upwell_const",1)

  ! Read the wind data
  ! Preserve the value of year_o so the actual year does not change 
  ! to the last data year read when the data is read 
  ! ***huh?
  yr=year_o
  ! *** Name of the wind data file should be moved to the run 
  ! *** parameters file **
  open(unit = 12, file = "../sog-forcing/wind/SHcompRotFmt.dat", &
       status = "OLD", action = "READ")
  !OPEN(UNIT = 12,FILE = "input/wind/wtwice.dat",STATUS = "OLD", &
  !    ACTION = "READ")
  do xx = 1, wind_n   !one year of wind data
     read(12, *) wind(xx)%day, wind(xx)%month, wind(xx)%year, wind(xx)%time, &
          wind(xx)%zonal, wind(xx)%meridional
     ! Adjust for leap years
     do yy = 1, size(leap_year, 1)
        if (wind(xx)%year == leap_year(yy)) then
           wind(xx)%leap = 1
           exit
        else if (yy == size(leap_year,1)) then
           wind(xx)%leap = 0
        end if
     end do
  end do
  close(12)
  ! Calculate the year-day (erroneously called Julian day) from the
  ! calendar dates read from the wind data file
  ! *** Change Julian day a less deceiving name (year_day?) and replace
  ! *** the Julian_day subroutine with a more robust algorithm **
  CALL Julian_day(wind%leap, wind%month, wind%day, wind%Jday, wind_n) 

  ! Read the cloud fraction data
  ! *** Name of the cloud faction data file should be moved to the 
  ! *** run parameters file **
  open(12, file="../sog-forcing/met/YVRCombCF.dat", &
       status = "OLD", action = "READ")
  ! *** This appears to be some number of year days ~(5.25 * 365)
  ! *** Should be defined elsewhere, or better, calculated from the data read
  ndays = 1918
  do ic = 1, ndays
     read(12, *) stn, year, month, day, para, (cf(ic,j), j = 1, 24)
  enddo

  close (12)

  ! Read the air temperature data
  ! *** Name of the air temperature data file should be moved to the 
  ! *** run parameters file **
  open(12, file="../sog-forcing/met/YVRCombATemp.dat", &
       status = "OLD", action = "READ")
  do ic = 1, ndays
     read(12, *) stn, year, month, day, para, (atemp(ic,j), j = 1, 24)
     !         write (*,*) day
  enddo
  close(12)
  ! *** Should confirm that dates are in sync with wind data **
  ! Take out the times 10 scaling of the raw temperature data and
  ! convert temperatures from Celius to Kelvin
  ! *** We need a C2K() utility function
  do ic = 1, ndays
     do j = 1, 24
        atemp(ic,j) = atemp(ic,j) / 10. + 273.15
     enddo
  enddo

  ! Read air humidity data
  ! *** Name of the humidity data file should be moved to the 
  ! *** run parameters file **
  open(12, file="../sog-forcing/met/YVRCombHum.dat", &
       status = "OLD", action = "READ")
  do ic = 1, ndays
     read(12, *) stn, year, month, day, para, (humid(ic,j), j = 1, 24)
  enddo
  close(12)
  ! *** Should confirm that dates are in sync with wind data **

  ! Read Fraser River flow data
  open(12, file="../sog-forcing/rivers/Fraser_2001_2005.dat", &
       status = "OLD", action = "READ")
  ! *** Another hard-coded constant to get rid of **
  do ic = 1, 1826
     read(12, *) year, month, day, Qriver(ic)
  enddo
  close(12)

  ! Read Englishman River flow data
  open(12, file="../sog-forcing/rivers/eng200123456.dat", &
       status = "OLD", action = "READ")
  ! *** Another hard-coded constant to get rid of **
  do ic = 1, 1919
     read(12, *) year, month, day, Eriver(ic)
  enddo
  close(12) 

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

