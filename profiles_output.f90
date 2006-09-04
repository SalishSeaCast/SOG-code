module profile_mod

implicit none

private
public init_profiles, profile, profile_close

! variables common to this module
integer :: maxprofiles
parameter (maxprofiles=20) ! don't worry about allocating for this, just set maximum to maxprofiles and check for it.
integer :: noprof, profday(maxprofiles) ! number of profiles, day of profiles
double precision :: proftime(maxprofiles) ! time of day for profiles
integer :: iprof ! counter, current profile

contains

subroutine init_profiles ()

! reads in from the control file, the days and times that profiles are required

character*80 haloclinefile

! function types
INTEGER getpari
CHARACTER*80 getpars
EXTERNAL getpari,getpars

 noprof =  getpari("noprof",1)

if (noprof.gt.20) then
   write (*,*) "Code written for a maximum of 20 profiles"
   stop
endif

if (noprof.gt.0) then

   call getpariv ("profday",profday,noprof,1)
   call getpardv ("proftime",proftime,noprof,1)

   haloclinefile = getpars ("haloclinefile",1)
   open(591,file=haloclinefile)

endif

iprof = 1

end subroutine init_profiles



subroutine profile (day, time, dt, depth, S, M)

! check and see if a profile is appropriate and do it

integer, intent(in):: day
double precision, intent(in):: time, dt ! can't expect exact time match
integer, intent(in):: M
double precision, intent (in) :: S(0:M), depth(0:M) ! salinity and depth of salinity

! internal variables

integer :: i ! depth counter
double precision :: delS, dep, derS ! delta S between two depth points, dep half way between them and -dS/dz

if (iprof.le.noprof.and.day.eq.profday(iprof)) then ! check the day
   if (abs(time-proftime(iprof)).le.0.5*dt) then ! check the time
      delS = 0.
      do i=1,M
         if (S(i)-S(i-1).gt.delS) then ! find the maximum difference
            delS = S(i)-S(i-1)
            dep = 0.5*(depth(i)+depth(i-1)) ! find the depth
            derS = delS/(depth(i)-depth(i-1)) ! find the gradient
         endif
      enddo
      write (591,1591) iprof,profday(iprof),proftime(iprof),derS,dep
1591  format (i2,1x,i3,1x,f13.5,1x,f13.5,1x,f7.2)
      iprof = iprof+1
   endif
endif

end subroutine profile

subroutine profile_close

close(591)

end subroutine profile_close


end module profile_mod
