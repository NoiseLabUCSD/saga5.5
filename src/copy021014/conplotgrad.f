c  Plotting contour plots for the GA optimization program
c  Peter Gerstoft
c

      subroutine conplot() 
c  Plots ambiguity surface between the two first parameters
      USE global
      INCLUDE 'comopt.h'
      integer maxtheta
      parameter(maxtheta=2)
c-----DESCRIPTION
c       Plotting of various parameters
c-----parameter vector
      real theta(maxtheta)
c---- help variables
      real v
c-----counters
c      integer i,j,k,n,m,i1
c---- object function
      integer mpoint
      parameter (mpoint=310)
      real fobj(mpoint,mpoint)
      integer n1,n2,ilog
      character*80 conopt      
      real XLEFT,XRIGHT,XSCALE,XINC,ydown,yup,ySCALE,yINC,
     &     dx,dy,ypointup,ypointdown,xpointleft,xpointright,
     &     zmin,zmax,zstep,zmax10,zdiv,xmax,ymax,xleft0,ydown0,
     &     xmin,ymin,zmin10
      real xdivdum,t1
      integer nxdifdum
      character*80 titlex,titley,TITLEH
      integer ixlay,iylay,ixpar,iypar,ix3,iy3,j,i
c-----calculate data and gradient
      zmin=10e20
      zmax=-10e20
      n1=min(mpoint,ndigit(1))
      n2=min(mpoint,ndigit(2))
      if (iopt(13).eq.1 .and. iopt(20).eq.0) then
c saclant conopt='UNI,COL,NCL,SCS                '
          conopt='COL,NCL,SCS                '
      else
c saclant conopt='UNI,COL,REV,NCL,SCS            '
          conopt='COL,REV,NCL,SCS            '
      endif
      ixlay= par2lay(1)
      ixpar= par2phy(1) 
      iylay= par2lay(2)
      iypar= par2phy(2)
      if (iopt(12).eq.1) then
         ix3=par3(1)
         iy3=par3(2)
      endif 
c     write(titlex,'(''Layer'',i2,'' parameter'',i2)')ixlay,ixpar
c     write(titley,'(''Layer'',i2,'' parameter'',i2)')iylay,iypar
      titlex=phystxt2(par2phy(1))
      titley=phystxt2(par2phy(2))
      xleft  = fmin(1)
      xleft0 = fmin(1)
      xright = fmax(1)
      
      ydown  = fmin(2)
      ydown0 = fmin(2)
      yup    = fmax(2)
      
      if (n1.gt.1) then
         dx=(xright-xleft)/(n1-1)
      else 
         dx=0
      endif
      if (n2.gt.1) then
         dy=(yup-ydown)/(n2-1)
      else 
         dy=0
      endif
      
      if (par2phy(1).eq.9) then
         xLEFT=xLEFT/1000
         xright=xright/1000
         titlex='Source range (km)'
      endif
      if (par2phy(2).eq.9) then
         ydown=ydown/1000
         yup=yup/1000
         titley='Source range (km)'
      endif
      CALL AUTOAX(Xleft,Xright,Xpointleft,Xpointright,
     1     XINC,XDIVdum,NXDIFdum)
      xinc=xinc/2
      yinc=(yup-ydown)/5
c     xinc=(xright-xleft)/5
      CALL AUTOAX(ydown,yup,Ypointdown,Ypointup,
     1     yINC,XDIVdum,NXDIFdum)
      yinc=yinc/2
      ypointup   = yup
      ypointdown  =ydown
      xpointleft =xleft
      xpointright=xright
      XSCALE= ABS(XRIGHT-XLEFT)/11.   ! x - axis are 11cm.
      YSCALE= ABS(yup-ydown)/11. ! Y- AXIS IS 11 CM.
c--   y-loop
      do j=1,n2
       if (j.eq.2) then
          CALL rdtime(t1)
          WRITE(*,320) t1*n2
 320      FORMAT(' >>> Estimated contour plot time:',F12.3,' secs.')
       endif
c---  x-loop 
         do i=1,n1
c     change of x axis          
            theta(1)=(i-1)*dx+xleft0        
c     change of y axis          
            theta(2)=(j-1)*dy+ydown0
            if (iopt(12).eq.0) then
               call setmodelx(ixlay,ixpar,theta(1))
               call setmodelx(iylay,iypar,theta(2))
            else
c      write(*,*)ixlay,ixpar,ix3,theta(1)
c      write(*,*)iylay,iypar,iy3,theta(2)
               call setmodelx(ixlay,ixpar,ix3,theta(1))
               call setmodelx(iylay,iypar,iy3,theta(2))
            endif
            call forw2grad
            call norm
            call cost(v)

c      write(*,*)'before:', i,j,v
c            if (iopt(13).eq.1 .and. iopt(20).eq.0) then
c               if (v.ge.0) then
c                  v=10.*log10(  abs(v)  )
c               else
c                  v=-300
c               endif
c            else
            if (iopt(8).eq.3) then
               if (v.ne.1) then
                  v=10.*log10(  abs(1.-v)  )
               else
                  v=-1000
               endif
            else
               if (v.ne.0) then
                  v=- 10.*log10(  abs(v)  )
               else
                  v=-1000
               endif
            endif
c      write(*,*)'after',i,j,v
c     if (iopt(30).eq.6) then
c     V=V*100
c     endif

c     v=10.*log10((1-v)*ndep)
c     if (v.lt.v0) 
c     zmin = min(v,zmin)
            if (v.gt.zmax) then
               xmax=theta(1)
               ymax=theta(2)
               zmax = v
            endif
            if (v.lt.zmin) then
               xmin=theta(1)
               ymin=theta(2)
               zmin = v
            endif
            fobj(i,n2+1-j)=v       
         enddo
      enddo

      if (zmax.eq.-1000) then
        write(*,*)' *** The data has been matched precisely. ***'
        write(*,*)' *** The dynamic scale must be found manually***'
      endif
      write(*,*)' The Unscaled  minimum and maximum are:'
      write(*,*)' zmin,ymin,xmin',zmin,ymin,xmin
      write(*,*)' zmax,ymax,xmax',zmax,ymax,xmax
      write(prtfil,*)' The Unscaled  minimum and maximum are:'
      write(prtfil,*)' zmin,ymin,xmin',zmin,ymin,xmin
      write(prtfil,*)' zmax,ymax,xmax',zmax,ymax,xmax
     
       if (iopt(8).eq.1) then
c pg feb98
       do j=1,n2
c---  x-loop 
         do i=1,n1
            fobj(i,n2+1-j)=  fobj(i,n2+1-j)-zmax      
        enddo
      enddo
      zmax=0
      zmin=zmin-zmax
      write(*,*)' The minimum and maximum are:'
      write(*,*)' zmin,ymin,xmin',zmin,ymin,xmin
      write(*,*)' zmax,ymax,xmax',zmax,ymax,xmax
      write(prtfil,*)' The minimum and maximum are:'
      write(prtfil,*)' zmin,ymin,xmin',zmin,ymin,xmin
      write(prtfil,*)' zmax,ymax,xmax',zmax,ymax,xmax
      endif
c     *** determine factor
      zmax=1.0*zmax
      ILOG=IFIX(ALOG10(abs(zMAX)))
      IF (zMAX.LT.1.0) ILOG=ILOG-1
      ilog  = ILOG-1
c---  for a plot in dB
      ilog=0
c      if (iopt(13).eq.1 .and. iopt(20).eq.0) then
c         zmax=zmin+5.5
c      else
         zmin=zmax-5.5
c      endif
c---  
      zDIV  =  10.**(-ilog)
      zmax10= zmax*zdiv
      zmin10= zmin*zdiv
      zstep = (zmax10-zmin10)/11 !/22
      zmin10=zmin10+zstep
      zmax10=zmax10-zstep
      write(*,*)' ilog',ilog
      write(titleH,'(a,i4,a)')'Object Function (10**', ilog ,')'
      write(*,*)'min and max for plotting are:',zmin10,zmax10 
C---  WRITEOUT CDR FILE
      CALL CONDRW1(TITLE,N1,n2,n1,n2,XLEFT,XRIGHT,XSCALE,XINC,
     1     yUP,yDOWN,YSCALE,yINC,ZMIN10,ZMAX10,
     2     ZSTEP,5.,1.,ypointup,ypointdown,xpointleft,xpointright,
     &     1,titlex,titley,conopt)
c--   y-loop
      do j=1,n2
         do i=1,n1
            fobj(i,j)=fobj(i,j)*zdiv ! *zmax10/zmax
         enddo
         CALL CONDRB(1,1,n1,fobj(1,j))
      enddo
      close(28)
      close(29)
      return
      end
c**********************************
      subroutine conppd2d(ppd2d) 
      USE global
      INCLUDE 'comopt.h'
      real xdivdum
      integer nxdifdum
      REAL*4 ppd2d(0:(mdig-1),0:(mdig-1))
c      integer ippd1,ippd2
c-----DESCRIPTION
c       Plotting of various parameters
c-----counters
c      integer i,j,k,n,m,i1,k1,l
c---- object function
      integer ilog,j,i
      real XLEFT,XRIGHT,XSCALE,XINC,ydown,yup,ySCALE,yINC,
     &     ypointup,ypointdown,xpointleft,xpointright
      real zmin,zmax,zstep,zmax10,zdiv
      character*80 titlex,titley
      character*80 conopt      
c      write(titlex,'(''Layer'',i2,'' parameter'',i2)')ixlay,ixpar
c      write(titley,'(''Layer'',i2,'' parameter'',i2)')iylay,iypar
       titlex=phystxt2(par2phy(ippd1))
       titley=phystxt2(par2phy(ippd2))
c saclant conopt='REV,COL,UNI,SCS,NCL       '
       conopt='REV,COL,SCS,NCL       '
       xleft  = fmin(ippd1)
       xright = fmax(ippd1)
       xinc=(xright-xleft)/5
       
       ydown  = fmin(ippd2)
       yup    = fmax(ippd2)
       yinc=(yup-ydown)/5

        CALL AUTOAX(Xleft,Xright,Xpointleft,Xpointright,
     1          XINC,XDIVdum,NXDIFdum)
        xinc=xinc/2
        CALL AUTOAX(ydown,yup,Ypointdown,Ypointup,
     1          yINC,XDIVdum,NXDIFdum)
      yinc=yinc/2

      ypointup   = yup
      ypointdown =ydown
      xpointleft =xleft
      xpointright=xright
      XSCALE= ABS(XRIGHT-XLEFT)/11.   ! x - axis is 11cm.
      YSCALE= ABS(yup-ydown)/11.      ! Y- AXIS IS 11 CM.
        zmin=10e10
        zmax=-10e10
c--   x loop
      do j=0,ndigit(ippd1)-1
c--      yloop
         do i=0,ndigit(ippd2)-1
            zmax = max(ppd2d(j,i),zmax)
         enddo
      enddo
c *** determine factor
         zmax=0.25*zmax
        ILOG=IFIX(ALOG10(zMAX))
        IF (zMAX.LT.1.0) ILOG=ILOG-1
        ilog  = ILOG-1
        ilog=0
        zmin=0
        zDIV  =  10.**(-ilog)
        zmax10= zmax*zdiv
        zstep = (zmax10-zmin)/10       !/22 or /11
        zmin=zmin+zstep
        zmax10=zmax10-zstep
C        write(*,*)' ilog',ilog
c        write(titleH,'(a,i4,a)')'2D marginal ppd (10**', ilog ,')'
C        write(*,*)'min and max are:',zmin,zmax 
C--- WRITEOUT CDR FILE
        CALL CONDRW1(TITLE,ndigit(ippd1),ndigit(ippd2),
     &     ndigit(ippd1),ndigit(ippd2),XLEFT,XRIGHT,XSCALE,XINC,
     1     yUP,yDOWN,YSCALE,yINC,ZMIN,ZMAX10,
     2     ZSTEP,5.,1.,ypointup,ypointdown,xpointleft,xpointright,
     &     1,titlex,titley,conopt)
c
c-- writeout BDR file
c
c--   y-loop
        do j=ndigit(ippd2)-1,0,-1
           do i=0,ndigit(ippd1)-1
              ppd2d(i,j)=ppd2d(i,j)*zmax10/zmax
           enddo
          CALL CONDRB(1,1,ndigit(ippd1),ppd2d(0,j))
        enddo
      close(28)
      close(29)
      return
      end
