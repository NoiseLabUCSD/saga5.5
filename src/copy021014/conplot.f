c  Plotting contour plots for the GA optimization PROGRAM
c  Peter Gerstoft
c

      SUBROUTINE conplot() 
c  Plots ambiguity surface between the two first parameters
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER maxtheta
      PARAMETER(maxtheta=2)
c-----DESCRIPTION
c       Plotting of various parameters
c-----PARAMETER vector
      REAL theta(maxtheta)
c---- help variables
      REAL v
c-----counters
c      INTEGER i,j,k,n,m,i1
c---- object FUNCTION
      INTEGER mpoint
c      PARAMETER (mpoint=310)
c      REAL fobj(mpoint,mpoint)
      REAL,    DIMENSION(:,:), ALLOCATABLE ::  fobj
      INTEGER n1,n2,ilog
      CHARACTER*80 conopt      
      REAL XLEFT,XRIGHT,XSCALE,XINC,ydown,yup,ySCALE,yINC,
     &     dx,dy,ypointup,ypointdown,xpointleft,xpointright,
     &     zmin,zmax,zstep,zmax10,zdiv,xmax,ymax,xleft0,ydown0,
     &     xmin,ymin,zmin10
      REAL xdivdum,t1
      INTEGER nxdifdum
      CHARACTER*80 titlex,titley,TITLEH
      INTEGER ixlay,iylay,ixpar,iypar,ix3,iy3,j,i
c-----calculate DATA and gradient
      zmin=10e20
      zmax=-10e20
      mpoint=max(ndigit(1),ndigit(2))
      allocate(fobj(mpoint,mpoint))
      n1=MIN(mpoint,ndigit(1))
      n2=MIN(mpoint,ndigit(2))
      IF (iopt(13).EQ.1 .AND. iopt(20).EQ.0) THEN
c saclant conopt='UNI,COL,NCL,SCS                '
         conopt='COL,NCL,SCS                '
      ELSE
c saclant conopt='UNI,COL,REV,NCL,SCS            '
         conopt='COL,REV,NCL,SCS            '
      ENDIF
      ixlay= par2lay(1)
      ixpar= par2phy(1) 
      iylay= par2lay(2)
      iypar= par2phy(2)
      IF (iopt(12).EQ.1) THEN
         ix3=par3(1)
         iy3=par3(2)
      ENDIF 
c     WRITE(titlex,'(''Layer'',i2,'' parameter'',i2)')ixlay,ixpar
c     WRITE(titley,'(''Layer'',i2,'' parameter'',i2)')iylay,iypar
      titlex=phystxt2(par2phy(1))
      titley=phystxt2(par2phy(2))
      xleft  = fmin(1)
      xleft0 = fmin(1)
      xright = fmax(1)
      
      ydown  = fmin(2)
      ydown0 = fmin(2)
      yup    = fmax(2)
      
      IF (n1.GT.1) THEN
         dx=(xright-xleft)/(n1-1)
      ELSE 
         dx=0
      ENDIF
      IF (n2.GT.1) THEN
         dy=(yup-ydown)/(n2-1)
      ELSE 
         dy=0
      ENDIF
      
      IF (par2phy(1).EQ.9) THEN
         xLEFT=xLEFT/1000
         xright=xright/1000
         titlex='Source range (km)'
      ENDIF
      IF (par2phy(2).EQ.9) THEN
         ydown=ydown/1000
         yup=yup/1000
         titley='Source range (km)'
      ENDIF
      CALL AUTOAX(Xleft,Xright,Xpointleft,Xpointright,
     1     XINC,XDIVdum,NXDIFdum)
      xinc=xinc/2
      yinc=(yup-ydown)/5
c     xinc=(xright-xleft)/5
      CALL AUTOAX(ydown,yup,Ypointdown,Ypointup,
     1     yINC,XDIVdum,NXDIFdum)
      yinc       = yinc/2
      ypointup   = yup
      ypointdown = ydown
      xpointleft = xleft
      xpointright= xright
      XSCALE= ABS(XRIGHT-XLEFT)/11. ! x - axis are 11cm.
      YSCALE= ABS(yup-ydown)/11.    ! Y- AXIS IS 11 CM.
c--   y-loop
      DO j=1,n2
         IF (j.EQ.2) THEN
            CALL rdtime(t1)
            WRITE(*,320) t1*n2
 320        FORMAT(' >>> Estimated contour plot time:',F12.3,' secs.')
         ENDIF
c---  x-loop 
         DO i=1,n1
c     change of x axis          
            theta(1)=(i-1)*dx+xleft0        
c     change of y axis          
            theta(2)=(j-1)*dy+ydown0
            IF (iopt(12).EQ.0) THEN
               CALL setmodelx(ixlay,ixpar,theta(1))
               CALL setmodelx(iylay,iypar,theta(2))
            ELSE
c      WRITE(*,*)ixlay,ixpar,ix3,theta(1)
c      WRITE(*,*)iylay,iypar,iy3,theta(2)
               CALL setmodelx(ixlay,ixpar,ix3,theta(1))
               CALL setmodelx(iylay,iypar,iy3,theta(2))
            ENDIF
            CALL forw2
            CALL norm
            CALL cost(v)
c            write(*,*) 'unc(1)',Aun(1)
c      WRITE(*,*)'before:', i,j,v
c            IF (iopt(13).EQ.1 .AND. iopt(20).EQ.0) THEN
c               IF (v.GE.0) THEN
c                  v=10.*LOG10(  ABS(v)  )
c               ELSE
c                  v=-300
c               ENDIF
c            ELSE
            IF (iopt(8).EQ.3) THEN
               IF (v.NE.1) THEN
                  v = 10.*LOG10(  ABS(1.-v)  )
               ELSE
                  v = -1000
               ENDIF
            ELSEIF (iopt(8).EQ.4) THEN
               IF (v.NE.(xcov_sum/nbart)) THEN
                  v = 10.*LOG10(  ABS(xcov_sum/nbart-v)  )
               ELSE
                  v = -1000
               ENDIF
            ELSE
               IF (v.NE.0) THEN
                  v = - 10.*LOG10(  ABS(v)  )
               ELSE
                  v = -1000
               ENDIF
            ENDIF
c     WRITE(*,*)'after',i,j,v
c     IF (iopt(30).EQ.6) THEN
c     V=V*100
c     ENDIF

c     v=10.*LOG10((1-v)*ndep)
c     IF (v.LT.v0) 
c     zmin = MIN(v,zmin)
            IF (v.GT.zmax) THEN
               xmax=theta(1)
               ymax=theta(2)
               zmax = v
            ENDIF
            IF (v.LT.zmin) THEN
               xmin=theta(1)
               ymin=theta(2)
               zmin = v
            ENDIF
            fobj(i,n2+1-j)=v       
         ENDDO
      ENDDO

      IF (zmax.EQ.-1000) THEN
        WRITE(*,*)' *** The data has been matched precisely. ***'
        WRITE(*,*)' *** The dynamic scale must be found manually***'
      ENDIF
      WRITE(*,*)' The Unscaled  minimum and maximum are:'
      WRITE(*,*)' zmin,ymin,xmin',zmin,ymin,xmin
      WRITE(*,*)' zmax,ymax,xmax',zmax,ymax,xmax
      WRITE(prtfil,*)' The Unscaled  minimum and maximum are:'
      WRITE(prtfil,*)' zmin,ymin,xmin',zmin,ymin,xmin
      WRITE(prtfil,*)' zmax,ymax,xmax',zmax,ymax,xmax
     
      IF (iopt(8).EQ.1) THEN
c pg feb98
         DO j=1,n2
c---  x-loop 
            DO i=1,n1
               fobj(i,n2+1-j)=  fobj(i,n2+1-j)-zmax      
            ENDDO
         ENDDO
         zmax=0
         zmin=zmin-zmax
         WRITE(*,*)' The minimum and maximum are:'
         WRITE(*,*)' zmin,ymin,xmin',zmin,ymin,xmin
         WRITE(*,*)' zmax,ymax,xmax',zmax,ymax,xmax
         WRITE(prtfil,*)' The minimum and maximum are:'
         WRITE(prtfil,*)' zmin,ymin,xmin',zmin,ymin,xmin
         WRITE(prtfil,*)' zmax,ymax,xmax',zmax,ymax,xmax
      ENDIF
c     *** determine factor
      zmax=1.0*zmax
      ILOG=IFIX(ALOG10(ABS(zMAX)))
      IF (zMAX.LT.1.0) ILOG=ILOG-1
      ilog  = ILOG-1
c---  for a plot in dB
      ilog=0
c     IF (iopt(13).EQ.1 .AND. iopt(20).EQ.0) THEN
c     zmax=zmin+5.5
c     ELSE
      zmin=zmax-5.5
c     ENDIF
c---  
      zDIV  =  10.**(-ilog)
      zmax10= zmax*zdiv
      zmin10= zmin*zdiv
      zstep = (zmax10-zmin10)/11 !/22
      zmin10=zmin10+zstep
      zmax10=zmax10-zstep
      WRITE(*,*)' ilog',ilog
      WRITE(titleH,'(a,i4,a)')'Object Function (10**', ilog ,')'
      WRITE(*,*)'min and max for plotting are:',zmin10,zmax10 
C---  WRITEOUT CDR FILE
      CALL CONDRW1(TITLE,N1,n2,n1,n2,XLEFT,XRIGHT,XSCALE,XINC,
     1     yUP,yDOWN,YSCALE,yINC,ZMIN10,ZMAX10,
     2     ZSTEP,5.,1.,ypointup,ypointdown,xpointleft,xpointright,
     &     1,titlex,titley,conopt)
c--   y-loop
      DO j=1,n2
         DO i=1,n1
            fobj(i,j)=fobj(i,j)*zdiv ! *zmax10/zmax
         ENDDO
         CALL CONDRB(1,1,n1,fobj(1,j))
      ENDDO
      CLOSE(28)
      CLOSE(29)
      RETURN
      END
c**********************************
      SUBROUTINE conppd2d(ppd2d) 
      USE global
      INCLUDE 'comopt.h'
      REAL xdivdum
      INTEGER nxdifdum
      REAL*4 ppd2d(0:(mdig-1),0:(mdig-1))
c      INTEGER ippd1,ippd2
c-----DESCRIPTION
c       Plotting of various parameters
c-----counters
c      INTEGER i,j,k,n,m,i1,k1,l
c---- object FUNCTION
      INTEGER ilog,j,i
      REAL XLEFT,XRIGHT,XSCALE,XINC,ydown,yup,ySCALE,yINC,
     &     ypointup,ypointdown,xpointleft,xpointright
      REAL zmin,zmax,zstep,zmax10,zdiv
      CHARACTER*80 titlex,titley
      CHARACTER*80 conopt      
c     WRITE(titlex,'(''Layer'',i2,'' parameter'',i2)')ixlay,ixpar
c     WRITE(titley,'(''Layer'',i2,'' parameter'',i2)')iylay,iypar
      titlex=phystxt2(par2phy(ippd1))
      titley=phystxt2(par2phy(ippd2))
c     saclant conopt='REV,COL,UNI,SCS,NCL       '
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
     1     yINC,XDIVdum,NXDIFdum)
      yinc=yinc/2
      
      ypointup   = yup
      ypointdown =ydown
      xpointleft =xleft
      xpointright=xright
      XSCALE= ABS(XRIGHT-XLEFT)/11. ! x - axis is 11cm.
      YSCALE= ABS(yup-ydown)/11.      ! Y- AXIS IS 11 CM.
      zmin=10e10
      zmax=-10e10
c--   x loop
      DO j=0,ndigit(ippd1)-1
c--   yloop
         DO i=0,ndigit(ippd2)-1
            zmax = MAX(ppd2d(j,i),zmax)
         ENDDO
      ENDDO
c *** determine factor
      zmax=0.25*zmax
      ILOG=IFIX(ALOG10(zMAX))
      IF (zMAX.LT.1.0) ILOG=ILOG-1
      ilog  = ILOG-1
      ilog=0
      zmin=0
      zDIV  =  10.**(-ilog)
      zmax10= zmax*zdiv
      zstep = (zmax10-zmin)/10  !/22 or /11
      zmin=zmin+zstep
      zmax10=zmax10-zstep
C        WRITE(*,*)' ilog',ilog
c        WRITE(titleH,'(a,i4,a)')'2D marginal ppd (10**', ilog ,')'
C        WRITE(*,*)'min and max are:',zmin,zmax 
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
      DO j=ndigit(ippd2)-1,0,-1
         DO i=0,ndigit(ippd1)-1
            ppd2d(i,j)=ppd2d(i,j)*zmax10/zmax
         ENDDO
         CALL CONDRB(1,1,ndigit(ippd1),ppd2d(0,j))
      ENDDO
      CLOSE(28)
      CLOSE(29)
      RETURN
      END
