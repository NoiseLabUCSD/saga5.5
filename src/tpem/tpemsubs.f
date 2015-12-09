c **************** THIS FILE CONTAINS TPEM MODEL SUBROUTINES *****************

c Author: Amalia E. Barrios
c         NCCOSC RDT&E DIV 543
c         49170 Propagation Path
c         San Diego, CA  92152-7385
c         e-mail: barrios@nosc.mil
c         phone: (619) 553-1429
c         fax: (619) 553-1417

c Summary: These routines model tropospheric radiowave propagation over 
c          variable 
c          terrain and calculates propagation loss vs. height and range.
c          Propagation loss is displayed in dB contours on a height vs. range  
c          plot. TPEM is based on the split-step Fourier PE method and was 
c          originally developed from an early PE model called PEPC, written 
c          by Fred Tappert.  Propagation loss over variable terrain is modeled 
c      by shifting the field an appropriate number of bin widths corresponding
c       to the height of the ground.  The field is determined using the smooth
c          earth PE method. 

c ****************************************************************************

c Variables in small letters in parameter lists are variables that are input
c or passed to called subroutines.  
c Variables in CAPS in parameter lists are returned
c from the called subroutines.
      
c **************************** SUBROUTINE PESTEP ******************************

c Purpose:  Propagates the field by one output range step DROUT.

c Called from MAIN DRIVER PROGRAM

c Routines called: DOSHIFT, SINFFT, REFINTER, PHASE2, CALCLOS

c Parameter list:

c HMINTER = height above mean sea level (meters) at which MLOSS is referenced.
c ROUT = output range in meters.
c MLOSS() = Array containing the propagation loss values in centibels, at each 
c           output range point ROUT.
c JSTART = index at which the valid propagation loss values begin.
c JEND = index at which the valid propagation loss values end.
 
      subroutine pestep( hminter, vnp, rf, tr, rp, sv,
     +              ROUT, MLOSS, closs,
     +              JSTART, JEND, ifieldtype ,itpem_opt,iperror)  
      
      integer itpem_opt(*),iperror
      include 'tpem.inc'
      
      common / htvar / ylast, ycur, ycurm, kbin,iylast,iycur
      common / arrays / X(0:maxpts), Y(0:maxpts), filt(0:maxn4),
     +                  SK(0:maxpts), CK(0:maxpts), CM(0:maxpts),
     +                  SM(0:maxpts), xlst(0:maxpts), ylst(0:maxpts)
      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1,
     +                   antref, rpe, hlim(mxrout), slp(mxter), fter
      TYPE(inputvar) vnp
      TYPE(refractivity) rf
      TYPE(terrain) tr
      TYPE(refparam) rp
      TYPE(systemvar) sv 

      integer ifieldtype             
      logical fter
      complex closs(*)
      integer mloss(*)
      real*8 r,rout,routhelp
      integer istep
      save r, kt,istep                                             

      if( rout .eq. 0. ) then
         r = 0.
         istep=0
         rout=vnp%rmin
      else
        rout = rout + drout
      endif
      routhelp=rout-dr*0.1      
c      write(*,*)'dr,routhelp,drout,fter',r,rout,dr,routhelp,drout,fter
      DO while( r .lt. routhelp )
         ylast = ycur 
         iylast = iycur 
         rlast = r
c      if (r.gt.110000) then
c      write(*,*)'dr,routhelp,rout/dr,r/dr',
c     &           dr,routhelp,rout/dr,r,r/dr,r-rout+dr,rlast/dr
c      endif   
c Store the field arrays of the previous range step for subsequent horizontal
c interpolation at range ROUT.
                                                     
         do i = 0, n
            xlst(i) = x(i)
            ylst(i) = y(i)
         end do 
                   
c         r = r + dr
         istep =istep+1
         r=dr*istep
c         write(*,*)'r,r1',r,r1

         rmid = r - dr2 
c         write(*,*)'r,rout',r,rout
      
         if( fter ) then
            if( r .eq. dr ) then
               slope = slp(1)
               kt = 1
            end if
            
c Check to see if current range is past a range point in terrain profile.
c If so, increment counter, determine terrain height at current range.
            
            do while((r .ge. tr%terx(kt+1)) .and. (kt .lt. tr%itp))
               kt = kt + 1
            end do        
            slope = slp(kt)
            ycur = tr%tery(kt) + slope * ( r - tr%terx(kt) )
            if(( kt .eq. tr%itp ) .and. (ycur .lt. 0.)) ycur = 0.
            ydif = ycur - ylast
            iycur = nint(ycur/delz)  
            kbin = abs(  iycur - iylast)
             
c Determine height at 1/2 range step - for interpolation on refractivity 
c profiles.
           
            kp = kt
            do while( rmid .lt. tr%terx(kp) ) 
               kp = kp-1
            end do
            ycurm = tr%tery(kp) + slp(kp) * (rmid - tr%terx(kp))
            if(( kp .eq. tr%itp ) .and. (ycurm .lt. 0.)) ycurm = 0.
         end if   
                                          
c Perform boundary shift for terrain case.
                                          
         if(( fter ) .and. ( slope .lt. 0. )) then
            if( kbin .gt. 0 ) call doshift( ydif, X, Y, n, nm1, kbin )
         end if

C  TRANSFORM TO FOURIER SPACE    
      
         call sinfft( ln, X )
         call sinfft( ln, Y )   
      
C  Multiply by free-space propagator.
                       
         DO I = 1, NM1
            T=X(I)*CK(I)+Y(I)*SK(I)
            Y(I)=CK(I)*Y(I)-SK(I)*X(I)
            X(I)=T
         end do
             
C  TRANSFORM BACK TO Z-SPACE  

         call sinfft( ln, X )
         call sinfft( ln, Y )

c If range-dependent and/or terrain case, then interpolate on profile.

         if(( rf%nprof .gt. 1 ) .or. ( fter ) 
     1     .or. (itpem_opt(6).eq.1)) then
            if (itpem_opt(6).eq.1)  then
c this reads in the markov chain base height
              call rdpara( rf,rp, sv, rmid,vnp%rmax,iperror)
            elseif (itpem_opt(5).eq.1)  then
c this interpolates in the trilinear profile.
               call paraminter( rf,rp,sv, rmid,iperror)
            else      
               call refinter( rf, rmid, hminter)
            end if
            CALL PHASE2
         end if
      if (iperror.lt.0) then
         return
      endif  

c Multiply by environment term.

         DO I = 1, nm1
            T=X(I)*CM(I)-Y(I)*SM(I)
            Y(I)=X(I)*SM(I)+Y(I)*CM(I)
            X(I)=T
         end do

c Perform boundary shift for terrain case.

         if(( fter ) .and. ( slope .ge. 0. )) then
            if( kbin .gt. 0 ) call doshift( ydif, X, Y, n, nm1, kbin )  
         end if
         
      end do  

c Calculate propagation loss at range ROUT.
c      write(*,*)'calculating loss for range',r,rlast,ifieldtype
      if (ifieldtype.eq.1)then
        call calclosc( r, real(rout), rlast, vnp, cLOSS, JSTART, JEND ) 
      else
        call calclos ( r, real(rout), rlast, vnp, MLOSS, JSTART, JEND ) 
      endif   
      end
           
c ************************** SUBROUTINE ANTPAT ******************************

c Purpose: Determines the antenna pattern factor based on angle passed.

c Parameter list:

c SANG = sine of angle
c PATFAC = antenna pattern factor

c Called from: XYINIT

c Routines called:  NONE
 
      subroutine antpat( iptrn, sang, PATFAC )
      
      common / pattern / pelev, afac, bw, elv, umax, sbw

c In the following pattern definitions, "u" refers to the angle for which 
c the antenna pattern is sought, and "u0" refers to the elevation angle.

c  IPTRN = 0 gives Omnidirectional antenna pattern factor : f(u) = 1 

      patfac = 1.

      if( iptrn .gt. 1 ) then
         u = asin( sang ) 
         udif = u - elv
      end if
         
c  IPTRN = 1 gives Gaussian antenna pattern based on
c  f(p-p0) = exp(-w**2 * ( p-p0 )**2 ) / 4, where p = sin(u) and 
c  p0 = sin(u0)

      if( iptrn .eq. 1 ) then
         pr = sang - pelev
         patfac = exp(-pr * pr * afac) 

c IPTRN = 2 gives sin(x)/x pattern based on 
c f(u-u0) = sin(x) / x where x = afac * sin(u-u0) for |u-u0| <= umax
c f(u-u0) = .03 for |u-u0| > umax
c IPTRN = 4 gives height-finder pattern which is a special case of sin(x)/x
        
      elseif(( iptrn .eq. 2 ) .or. ( iptrn .eq. 4 )) then
         if( iptrn .eq. 4 ) then
            dirang = abs( sang )
            if( dirang .gt. elv ) udif = u - dirang
         end if
         if( abs(udif) .le. 1.e-6 ) then
            patfac = 1. 
         elseif( abs( udif ) .gt. umax ) then
            patfac = .03 
         else
            arg = afac * sin( udif )                
            patfac = amin1( 1., amax1( .03, sin( arg ) / arg ) )
         end if
 
c IPTRN = 3 gives csc-sq pattern based on
c f(u) = 1 for u-u0 <= bw
c f(u) = sin(bw) / sin(u-u0) for u-u0 > bw
c f(u) = maximum of .03 or [1+(u-u0)/bw] for u-u0 < 0

      elseif( iptrn .eq. 3 ) then
         if( udif .gt. bw ) then
            patfac = sbw / sin( udif )
         elseif( udif .lt. 0 ) then
            patfac = amin1( 1., amax1( .03, (1. + udif/bw) ) )
         end if               
      end if

      end
                                                                             
c ******************** SUBROUTINE CALCLOS ************************************

c Purpose: Determines the propagation loss at each output range step ROUT.

c Called from subroutines: PESTEP
   
c Variables input:
c   R = PE range step in meters
c   ROUT = output range in meters
c   RLAST = last PE range in meters
c   FTER = logical flag - .TRUE.=terrain problem, .FALSE.=smooth surface

c Variables returned:
c   MLOSS() = integer array containing propagation loss values in centibels.
c   JSTART = index at which valid loss values in MLOSS begin.
c   JEND = index at which valid loss values in MLOSS ends.

      subroutine calclos( r, rout, rlast, vnp, MLOSS, JSTART, JEND )

      include 'tpem.inc'
                                  
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1,
     +                   antref, rpe, hlim(mxrout), slp(mxter), fter
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      common / arrays / X(0:maxpts), Y(0:maxpts), filt(0:maxn4),
     +                  SK(0:maxpts), CK(0:maxpts), CM(0:maxpts),
     +                  SM(0:maxpts), xlst(0:maxpts), ylst(0:maxpts)
      common / htvar / ylast, ycur, ycurm, kbin,iylast,iycur
      
      TYPE(inputvar) vnp
      
      integer mloss(*)
      real*8 r
      logical fter
      
      dimension rfac1(mxzout), rfac2(mxzout)
      
      save ic
                      
      data pfacmin / 300. /

c Define in-line function for linear interpolation.
                         
      plint(pl1, pl2, frac) = pl1 + frac * ( pl2 - pl1 )

c Initialize counter for HLIM array.
      
      if( rout .eq. drout ) ic = 1
      
c Get height of ground at output range ROUT and determine number of vertical
c output points that correspond to the ground height.  Fill the loss array MLOSS with
c zeros to represent ground for those vertical output points.
      
      xx = (rout - rlast) / dr
      zint = plint( ylast, ycur, xx )  
      izg = int( zint / dzout  )
      do i = 1, izg
         mloss(i) = 0  
      end do    
      jstart = izg + 1 
 
      if( rout .ge. rpe ) then
      
c If current output range is greater than RPE then begin calculation of 
c loss values and return them in MLOSS().
      
         rloglst = 0.                    
         if( rlast .gt. 0. ) rloglst = 10. * log10( rlast )
         rlog = 10. * log10( r )  
         fslrout = 20. * alog10(rout) + plcnst       !free space loss at ROUT
 
c Determine values of array elements corresponding to the ground and set these
c to the minimum propagation factor (-300) for later interpolation.

         if( fter ) then
            ip1 = int( ylast / dzout )
            ip2 = int( ycur / dzout ) 

            do i = 1, ip1
               rfac1(i) = pfacmin
            end do
            do i = 1, ip2
               rfac2(i) = pfacmin
            end do  
            ip1 = ip1 + 1
            ip2 = ip2 + 1 
         else
            ip1 = 1
            ip2 = 1
         end if

c Determine height/integer value at which to stop calculating loss.  
c NOTE: For terrain cases, ray tracing was performed for smooth (flat) surface
c reflection and sometimes HLIM(i) may be less than the local ground height.
c The GOTO statement is used just as a safety factor in this case.
                                     
c         zend1 = amax1( zint, hlim(ic) )
c         zend2 = amin1( real(vnp.hmax), zend1 ) 
c         jend = nint( zend2 / dzout )
c         write(*,*)'  zend1,zend2,zint, hlim(ic), dzout', 
c     1                zend1,zend2,zint, hlim(ic), dzout
         jend=vnp%nzout
c         if( jend .le. jstart ) goto 5

c Get propagation factor at valid heights from field at previous range step.
           
         if( rloglst .gt. 0. ) then
            do i = ip1, jend   
               zht = zout(i) - ylast
               rfac1(i) = getpfac( xlst, ylst, rloglst, zht )
            end do
         end if
         
c Get propagation factor at valid heights from field at current range step.
         
         do i = ip2, jend
            zht = zout(i) - ycur  
            rfac2(i) = getpfac( x, y, rlog, zht )
         end do

c Interpolate between the two PE range steps to get loss at range ROUT.
                                       
c        write(*,*)'fslrout,rloglst',fslrout,rloglst,xx,r,rout,dr
         do k = jstart, jend                        
            if( rloglst .gt. 0. ) then
               loss = 10.*( plint( rfac1(k), rfac2(k), xx ) + fslrout )
               mloss(k) = int( loss ) 
            else
               mloss(k) = int( 10. * ( rfac2(k) + fslrout ) )
            end if
         end do
      
 5       continue 
 
c Fill remainder of array with -1 indicating non-valid loss values.
       
         jn = jend + 1
         do i = jn, vnp%nzout
            mloss(i) = -1
         end do
      
      else

c If output range is less than RPE then there are no current valid loss values
c at any height - fill MLOSS with -1.
c JSTART and JEND will be equal and will have a value
c of 1 if smooth surface case, otherwise will have a value of the nearest
c integer multiple of DZOUT corresponding to the height of the local ground.

         jend = jstart
         do i = jstart, vnp%nzout
            mloss(i) = -1
         end do
         
      end if
      
      ic = ic + 1
      
      end 
                                                                             
c ***************************** SUBROUTINE CALCLOSc ************************

c Purpose: Determines the propagation loss at each output range step ROUT.

c Called from subroutines: PESTEP
   
c Variables input:
c   R = PE range step in meters
c   ROUT = output range in meters
c   RLAST = last PE range in meters
c   FTER = logical flag - .TRUE.=terrain problem, .FALSE.=smooth surface

c Variables returned:
c   MLOSS() = integer array containing propagation loss values in centibels.
c   JSTART = index at which valid loss values in MLOSS begin.
c   JEND = index at which valid loss values in MLOSS ends.

      subroutine calclosc( r, rout, rlast, vnp, CLOSS, JSTART, JEND )

      include 'tpem.inc'
                                  
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1,
     +                   antref, rpe, hlim(mxrout), slp(mxter), fter
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      real*8  dr, drout,dzout, dr2
      common / irhstps / izout(mxzout)
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      common / arrays / X(0:maxpts), Y(0:maxpts), filt(0:maxn4),
     +                  SK(0:maxpts), CK(0:maxpts), CM(0:maxpts),
     +                  SM(0:maxpts), xlst(0:maxpts), ylst(0:maxpts)
      common / htvar / ylast, ycur, ycurm, kbin,iylast,iycur
      
      TYPE(inputvar) vnp 
      
      complex closs(*),c_field
      real*8 r
      
      logical fter
      
      dimension rfac1(mxzout), rfac2(mxzout)
      
      save ic
                      
      data pfacmin / 300. /
      real rmult

c Get height of ground at output range ROUT and determine number of vertical
c output points that correspond to the ground height.  Fill the loss array MLOSS with
c zeros to represent ground for those vertical output points.
     
c      write(*,*) 'entering calclosc' 
      xx = (rout - rlast) / dr
       if( xx .lt. 0.995 ) then
          write(*,*)' xx,rout , rlast, zint, ylast,ycur'
          write(*,*)  xx,rout , rlast, zint, ylast,ycur
          write(*,*)' xx ',xx 
          stop ' only range-step size to grid size for complex field'
       end if

c      izg = int( (ycur) / dzout) 
      izg = 1+int( (ycur-vnp%hmin) / dzout  -.001)
c      write(*,*)'izg,ycur,vnp.hmin, dzout'
c      write(*,*)izg,ycur,vnp.hmin, dzout,
c     .        max1(ycur-vnp.hmin,delz) / dzout 
      if (izg.gt.0) then
        do i = 1, izg
           closs(i) = 0  
        end do    
        jstart = izg + 1 
      else
        jstart=1
      endif 
      
      rmult=1/(2*fko*sqrt(r))
 
c Determine height/integer value at which to stop calculating loss.  
c NOTE: For terrain cases, ray tracing was performed for smooth (flat) surface
c reflection and sometimes HLIM(i) may be less than the local ground height.
c The GOTO statement is used just as a safety factor in this case.
                                     
         jend=vnp%nzout
c Interpolate between the two PE range steps to get loss at range ROUT.
                                       
c         write(*,*)'xx,r,rout,dr',xx,r,rout,dr
c          write(*,*)'jstart,jend',jstart,jend
c         write(*,*)'ycur,zout(jstart)',ycur,zout(jstart)
         do i = jstart, jend                        
               zht = zout(i) - ycur  
               izht = izout(i) - iycur  
c                call  complexfield(closs(i), x, y, zht)
c                closs(i) =closs(i)*rmult
                closs(i) =cmplx(x(izht),y(izht))*rmult
c          write(*,*)' i,izht,zht,ycur,iycur,izout(i),closs(i)'
c          write(*,*) i,izht,closs(i)
         end do
      end 
      
c ************************* SUBROUTINE DOSHIFT *****************************

c Purpose: Shifts the field by the # of bins corresponding to height of 
c          the ground.

c Called from: PESTEP

c Routines called: NONE
 
c Parameter list:

c YDIF = height difference between YCUR and YLAST
c X() = Real part of field
c Y() = Imaginary part of field
c N = Number of field points
c INCR = Increment multiplier: -1=negative slope, array elements shifted up.
c                               1=positive slope, array elements shifted down.
c KBIN = Number of bins the field is to be shifted.
 
      subroutine doshift( ydif, X, Y, n, nm1, kbin )

      dimension x(0:*), y(0:*)

c If slope is positive then shift array elements down.
                         
      if( ydif .ge. 0. ) then
         incr = 1
         jst = 1
         jend = nm1 - kbin
      else            
      
c If slope is negative then shift array elements up.

         incr = -1
         jst = nm1
         jend = kbin + 1
      end if

      kinc = incr * kbin
      do j = jst, jend, incr
         jk = j + kinc
         x(j) = x(jk)
         y(j) = y(jk)
      end do

      if( incr .gt. 0 ) then
         nst = n - kbin
         do j = nst, nm1
            x(j) = 0.
            y(j) = 0.
         end do
      else
         do j = 1, kbin
            x(j) = 0.
            y(j) = 0.
         end do
      end if

      x(0) = 0.
      y(0) = 0.

      end

c ***************************** FUNCTION GETPFAC ****************************

      function GETPFAC( x, y, rlog, height )

c Purpose:  Performs linear interpolation on the power and then
c           calculates propagation factor in dB. 

c Called from: CALCLOS

c Routines called: NONE

c Parameter list:
 
c X() = Real part of field
c Y() = Imaginary part of field
c RLOG = 10. * log( range )
c HEIGHT = receiver height in meters 

      common / pevar / wl, fko, delz, n, ln, zmax, n34, con

      dimension x(0:*), y(0:*)
      
      data powmin/1.e-13/

      fb = height / delz
      nb=int(fb)
      fr=fb-float(nb)
      nbp1=nb+1
      
      x0=x(nb)
      y0=y(nb)
      x1=x(nbp1)
      y1=y(nbp1)

       
      pow0 = sqrt( x0*x0+y0*y0 )
      pow1 = sqrt( x1*x1+y1*y1 )
      
c      write(*,*)'fr,nb,nbp1',fr,nb,nbp1

      pow = pow0 + fr * (pow1 - pow0)

      rpow = amax1( pow, powmin )
      getpfac = -20.*alog10( rpow ) - rlog
      
      end  

c     ************************ FUNCTION COMPLEXFIELD *************************
c      complex function complexfield( x, y,  height )
      subroutine complexfield( c_field, x, y,  height )

c Purpose:  Performs linear interpolation on the complex field
c Called from: CALCLOS
c Routines called: NONE

c Parameter list:
c X() = Real part of field
c Y() = Imaginary part of field
c RLOG = 10. * log( range )
c HEIGHT = receiver height in meters 

      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      complex c_field
      dimension x(0:*), y(0:*)
      data powmin/1.e-13/

      fb = height / delz
      nb=int(fb)
      fr=fb-float(nb)
      nbp1=nb+1
      
      x0=x(nb)
      y0=y(nb)
      x1=x(nbp1)
      y1=y(nbp1)
c      writE(*,*)' Does it crash here ?'
c      write(*,*)' x0,fr,x1,nb,height,delz',  x0,fr,x1,nb,height,delz
c      write(*,*)' x0 +fr*(x0-x1)', x0 +fr*(x0-x1)
c      write(*,*)'y0 +fr*(y0-y1 ', y0 +fr*(y0-y1)
c      complexfield = cmplx( x0 +fr*(x0-x1),  y0 +fr*(y0-y1) )
      c_field = cmplx( x0 +fr*(x1-x0),  y0 +fr*(y1-y0) )
c      c_field = cmplx( x0 +fr*(x0-x1),  y0 +fr*(y0-y1) )
c      write(*,*)'fr,nb,nbp1',fr,nb,nbp1,x0,y0,x1,y1,c_field
      return
      end  
