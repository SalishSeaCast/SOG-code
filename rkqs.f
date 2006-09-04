      SUBROUTINE rkqs(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs_sog,
     *check,TT)
      INTEGER n,NMAX
      DOUBLE PRECISION eps,hdid,hnext,htry,x,dydx(n),y(n),yscal(n),
     *TT(1:81)
      PARAMETER(NMAX = 1000)
      EXTERNAL derivs_sog, rkck
      INTEGER check,ii
      DOUBLE PRECISION errmax,hh,xnew,yerr(NMAX),ytemp(NMAX),SAFETY,
     *PGROW,PSHRNK,ERRCON
      PARAMETER (SAFETY=0.9,PGROW=-.2,PSHRNK=-.25,ERRCON=1.84e-4)
      hh=htry

 1    call rkck(y,dydx,n,x,hh,ytemp,yerr,derivs_sog,TT)
      errmax=0.
      !open(555,file="output/errmax.dat")
      do 11 ii=1,n
        errmax=max(errmax,abs(yerr(ii)/yscal(ii)))
      !  write(555,*)errmax,hh
11    continue
      errmax=errmax/eps

      if(errmax.gt.1.)then
        hh=SAFETY*hh*(errmax**PSHRNK)
        if(hh.lt.0.1*hh)then
          hh=.1*hh
        endif
        xnew=x+hh
        if(xnew.eq.x)pause 'stepsize underflow in rkqs'
        goto 1
      else
        if(errmax.gt.ERRCON)then
          hnext=SAFETY*hh*(errmax**PGROW)
!          print*,'errmax.gt.ERRCON',errmax
        else
          hnext=5.*hh
        endif
        hdid=hh
        x=x+hh
        do 12 ii=1,n
          y(ii)=ytemp(ii)
12      continue
        return
      endif
      END



