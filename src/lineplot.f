c     Plotting contour plots for the GA optimization PROGRAM
c     Peter Gerstoft
c     

      SUBROUTINE lineplot() 
c     Plots ambiguity surface between the two first parameters
      USE global
      INCLUDE 'comopt.h'

c-----DESCRIPTION
c     Plotting of various parameters
c-----PARAMETER vector
      REAL theta
c---- help variables
      REAL v
c-----counters
c     INTEGER i,j,k,n,m,i1
c---- object FUNCTION
      INTEGER mpoint
      PARAMETER (mpoint=2000)
      REAL fobj(mpoint),thet(mpoint)
      INTEGER n1
c     CHARACTER*80 conopt      
      REAL XLEFT,XRIGHT,dx,zmin,zmax,xmax,xmin
      REAL xstart(100)
      CHARACTER*80 titlex
      INTEGER ixlay,ixpar,ix3,i,ipar
c-----calculate DATA and gradient

      CALL forw2
      CALL norm
      CALL cost(v)
      IF (iopt(32).EQ.3) THEN
         IF (v.NE.1) THEN
            v=10.*LOG10(  ABS(1.-v)  )
         ELSE
            v=-1000
         ENDIF
      ELSE
         IF (v.NE.0) THEN
            v=- 10.*LOG10(  ABS(v)  )
         ELSE
            v=-1000
         ENDIF
      ENDIF

      CALL getmodelreal(xstart)

      WRITE(13,*)'figure;'
      WRITE(13,*)'fobj=[];'
      WRITE(13,*)'theta=[];'
      WRITE(13,*)'xtitle=[]; xl=[]; xr=[]; ixpar=[];'
      WRITE(13,*)'v0=',v,';'
      WRITE(13,*)'xstart=['
      DO ipar=1,nparm
         IF (par2phy(ipar).EQ.9) THEN
             WRITE(13,*)xstart(ipar)/1000
          Else
             WRITE(13,*)xstart(ipar)
          ENDIF
      ENDDO
      WRITE(13,*)'];'
      DO 100 ipar=1,nparm
         WRITE(*,*)' Computing objective function for parameter',ipar
         zmin=10e20
         zmax=-10e20
         n1=MIN(mpoint,ndigit(ipar))
         ixlay= par2lay(ipar)
         ixpar= par2phy(ipar) 
         IF (iopt(12).EQ.1) THEN
            ix3=par3(ipar)
         ENDIF 
c     WRITE(titlex,'(''Layer'',i2,'' parameter'',i2)')ixlay,ixpar
c     WRITE(titley,'(''Layer'',i2,'' parameter'',i2)')iylay,iypar
         titlex=phystxt2(par2phy(ipar))
         xleft  = fmin(ipar)
         xright = fmax(ipar)
         WRITE(*,*)' ixlay, ixpar, xleft, xright,ipar,n1'
         WRITE(*,*)ixlay, ixpar, xleft, xright,ipar,n1
         IF (n1.GT.1) THEN
            dx=(xright-xleft)/(n1-1)
         ELSE 
            dx=0
         ENDIF
         
         IF (par2phy(ipar).EQ.9) THEN
            xLEFT=xLEFT/1000
            xright=xright/1000
            dx=dx/1000
            titlex='Source range (km)'
         ENDIF
c     CALL AUTOAX(Xleft,Xright,Xpointleft,Xpointright,
c     1        XINC,XDIVdum,NXDIFdum)
c     xinc=xinc/2
c     XSCALE= ABS(XRIGHT-XLEFT)/11. ! x - axis are 11cm.
c---  x-loop 
         WRITE(*,*)'xleft,xright,dx', xleft,xright,dx        
         DO i=1,n1
c     change of x axis          
            WRITE(*,*)' Computing objective function for parameter',
     1           ipar,i
            theta=(i-1)*dx+xleft        
            thet(i)=theta
            IF (par2phy(ipar).EQ.9) THEN
               theta=theta*1000
            ENDIF
c     WRITE(*,*)'ipar,i,theta,n1',ipar,i,theta,n1
            IF (iopt(12).EQ.0) THEN
               CALL setmodelx(ixlay,ixpar,theta)
            ELSE
               WRITE(*,*)ixlay,ixpar,ix3,theta
               CALL setmodelx(ixlay,ixpar,ix3,theta)
            ENDIF
            CALL forw2
            CALL norm
            CALL cost(v)
            IF (iopt(6).EQ.1) THEN
               WRITE(*,*) 'for parameter, theta:',i, theta
               WRITE(*,*)'value of obj fct',v
            ENDIF
            IF (iopt(32).EQ.3) THEN
               IF (v.NE.1) THEN
                  v=10.*LOG10(  ABS(1.-v)  )
               ELSE
                  v=-1000
               ENDIF
            ELSE
               IF (v.NE.0) THEN
                  v=- 10.*LOG10(  ABS(v)  )
               ELSE
                  v=-1000
               ENDIF
            ENDIF


            IF (v.GT.zmax) THEN
               xmax=theta
               zmax = v
            ENDIF
            IF (v.LT.zmin) THEN
               xmin=theta
               zmin = v
            ENDIF
            fobj(i)=v       
         ENDDO

         IF (zmax.EQ.-1000) THEN
            WRITE(*,*)' *** The data has been matched precisely. '
            WRITE(*,*)' *** The dynamic scale must be found manually'
         ENDIF
         WRITE(*,*)' The Unscaled  minimum and maximum are:'
         WRITE(*,*)' zmin,xmin',zmin,xmin
         WRITE(*,*)' zmax,xmax',zmax,xmax
         WRITE(prtfil,*)' The Unscaled  minimum and maximum are:'
         WRITE(prtfil,*)' zmin,xmin',zmin,xmin
         WRITE(prtfil,*)' zmax,xmax',zmax,xmax
         
c--   scaling relative to max
         IF (iopt(32).EQ.1) THEN
c     pg feb98
c---  x-loop 
            DO i=1,n1
               fobj(i)=  fobj(i)-zmax      
            ENDDO
            zmax=0
            zmin=zmin-zmax
            WRITE(*,*)' The minimum and maximum are:'
            WRITE(*,*)' zmin,xmin',zmin,xmin
            WRITE(*,*)' zmax,xmax',zmax,xmax
            WRITE(prtfil,*)' The minimum and maximum are:'
            WRITE(prtfil,*)' zmin,xmin',zmin,xmin
            WRITE(prtfil,*)' zmax,xmax',zmax,xmax
         ENDIF
c--   scaling relative to max
c     *** determine factor
c     zmax=1.0*zmax
c     ILOG=IFIX(ALOG10(ABS(zMAX)))
c     IF (zMAX.LT.1.0) ILOG=ILOG-1
c     ilog  = ILOG-1
c---  for a plot in dB
c     ilog=0
c     IF (iopt(13).EQ.1 .AND. iopt(20).EQ.0) THEN
c     zmax=zmin+5.5
c     ELSE
c     zmin=zmax-5.5
c     ENDIF
c---  
c     zDIV  =  10.**(-ilog)
c     zmax10= zmax*zdiv
c     zmin10= zmin*zdiv
c     zstep = (zmax10-zmin10)/11 !/22
c     zmin10=zmin10+zstep
c     zmax10=zmax10-zstep
c     WRITE(*,*)' ilog',ilog
c     WRITE(titleH,'(a,i4,a)')'Object Function (10**', ilog ,')'
c     WRITE(*,*)'min and max for plotting are:',zmin10,zmax10 
C---  WRITEOUT CDR FILE
         WRITE(13,*) ' fobj = [ fobj; ['
         DO i=1,n1
c     fobj(i)=fobj(i)*zdiv                ! *zmax10/zmax
            WRITE(13,*)  fobj(i)
c     WRITE(*,*)'fobj:', fobj(i)
         ENDDO
         WRITE(13,*) ']''];'
         WRITE(13,*) ' theta = [ theta; ['
         DO i=1,n1
            WRITE(13,*)  thet(i)
         ENDDO
         WRITE(13,*) ']''];'
         WRITE(13,*) 'xl=[xl',xleft,'];'
         WRITE(13,*) 'xr=[xr',xright,'];'
         WRITE(13,*) 'ixpar=[ixpar',ixpar,'];'

C     
C---- Reset PARAMETER
C     
         WRITE(*,*)'ipar,xstart(ipar)',ipar,xstart(ipar)
         IF (iopt(12).EQ.0) THEN
            CALL setmodelx(ixlay,ixpar,xstart(ipar))
         ELSE
c     WRITE(*,*)ixlay,ixpar,ix3,theta(1)
            CALL setmodelx(ixlay,ixpar,ix3,xstart(ipar))
         ENDIF
 100  CONTINUE
      WRITE(13,*)'yrange=[min(min(fobj)) max(max(fobj))];'

      WRITE(13,*)'for jj=1:',nparm
      WRITE(13,*) '   subplot(3,',INT(nparm/3.+0.99)
     1     ,',jj);'
      WRITE(13,*) 'plot(theta(jj,:),fobj(jj,:),''linewidth'',1.2)'
      IF (iopt(32).NE.1) THEN
         WRITE(13,*) 'hold on; plot([xl(jj), xr(jj)',
     1        '],[ v0 v0],''r--'' )'
      ENDIF

      WRITE(13,*)'set(gca,''xlim'',[xl(jj), xr(jj)])'
      WRITE(13,*) 'xlabel( xtitles(',iopt(30),
     1     ',ixpar(jj)));'
      WRITE(13,*)' plot([xstart(jj) xstart(jj)],yrange,''r--'');'
      WRITE(13,*)'bordure'
      WRITE(13,*) 'end'
      RETURN
      END
c**********************************







