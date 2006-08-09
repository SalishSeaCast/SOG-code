      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,
     *derivs_sog,rkqs,check,Ti)
      INTEGER nbad,nok,lup,nvar,KMAXX,MAXSTP,NMAX
      DOUBLE PRECISION eps,h1,hmin,x1,x2,ystart(nvar),TINY,Ti(1:81)
      ! *** Temporaily added derivs_noflag until flagellates code works
      EXTERNAL derivs_sog, derivs_noflag, rkqs
      PARAMETER (MAXSTP=1000,NMAX=1000,KMAXX=200,TINY=1.e-20) !KMAXX = 200,TINY=1.e-30
!      PARAMETER (MAXSTP=1000,NMAX=400,KMAXX=200,TINY=1.e-20) !KMAXX = 200,TINY=1.e-30
      INTEGER i,kmax,kount,nstp,check,isusan
      DOUBLE PRECISION dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),
     *y(NMAX),yyp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yyp
      kmax = 0
      dxsav = 10.0
      x=x1
      h=sign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.*dxsav 
      do 16 nstp=1,MAXSTP
        
      call derivs_sog(x,nvar,y,dydx,Ti)

        do 12 i=1,nvar
          yscal(i)=abs(y(i))+abs(h*dydx(i))+TINY
12      continue

       if(kmax.gt.0)then
          if(abs(x-xsav).gt.abs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yyp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
       
       call rkqs(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs_sog,
     *check,Ti)

      !open(554,file="output/dydx.dat")
      !do isusan=1,400
      !write(554,*)dydx(isusan)
      !enddo
      !close(554) 

        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yyp(i,kount)=y(i)
15          continue
          endif
          return
        endif
        if(abs(hnext).lt.hmin) pause 
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END





