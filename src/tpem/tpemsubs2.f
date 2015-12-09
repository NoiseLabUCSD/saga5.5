c     ********************** SUBROUTINE PEINIT ***********************************

c Purpose:  Initializes all variables used in TPEM subroutines for calculations

c Parameter list:

c HMINTER = Height of the minimum elevation of terrain profile.  This is used
c           to adjust entire terrain profile so subsequent loss values returned
c           will be referenced to this height.
c ROUT = Output range point (meters) - initialized in this routine
c IERROR = Integer value that is returned if any errors exist in input data: 
c            -6 : Last range in terrain profile is less than RMAX. (Will only return
c                 this error if error flag EF.LERR6 is set to .TRUE.).
c            -8 : HMAX is less than maximum height of terrain profile.
c           -10 : Transform size needed to do the given problem is greater than
c                 the maximum transform size allowed. (Return from subroutine 
c                 GETFFTSZ).
c           -12 : Range of last refrac profile  (for range dependent case)
c                 is less than RMAX. (returned from subroutine REFINIT).  
c                 Will only return this error if  EF.LERR12 is .TRUE. 
c           -14 : Last gradient in any refract profile entered is negative.
c                 (This is returned from REFINIT).
c           -17 : Range points of terrain profile is not increasing.
c           -18 : First range point is not 0.

c Called from MAIN DRIVER PROGRAM 

c Routines called: REFINIT, TRACEA, GETFFTSZ, XYINIT, SINFFT, TRACEH,
C                  PHASE1, PROFREF, INTPROF, PHASE2
                 
      subroutine peinit( ef, vnp, rf, sv, tr, rp,
     1             HMINTER, ROUT, IERROR ,itpem_opt,iperror)
      
      integer itpem_opt(*),iperror
      
      include 'tpem.inc'
                                     
c Common Blocks

c ARRAYS:
c   X() = Real part of field solution.
c   Y() = Imaginary part of field solution.
c   FILT() = Cosine-tapered (Tukey) filter array.
c   CK() = Real part (cosines) of free-space propagator exponential term.
c   SK() = Imaginary part (sines) of free-space propagator exponential term.
c     i.e., CK() + i*SK() = exp[-i * dr * (k - sqrt(k**2 - p**2)) ]
c   CM() = Real part (cosines) of refractivity-exponential term.
c   SM() = Imaginary part (sines) of refractivity-exponential term.
c     i.e. CM() + i*SM() = exp[i * dr * k * 1e-6 * M(z) ]

c MISCVAR:
c   PI2 = value of pi / 2.
c   FNORM = normalization factor used for DFT.
c   CNST = used in calculating CK() and SK() in routine PHASE1.
c          CNST = DELP/FKO.
c   DELP = mesh size in angle- (or p-) space.
c   THETAMAX = maximum propagation angle in PE calculations.
c   PLCNST = constant used in determining propagation loss 
c            PLCNST = 20log(2*FKO).
c   ANTREF = transmitting antenna height relative to the reference
c            height HMINTER.
c   RPE = range at which valid loss values will begin to be calculated.
c   HLIM() = array containing height at each output range at which the
c            last valid loss value exists.
c   SLP() = slope of each segment of terrain.
c   FTER = logical flag - .TRUE.=terrain case, .FALSE.=smooth surface case

c HTVAR:
c   YLAST = height of ground at the last range step
c   YCUR = height of ground at the current range step
c   YCURM = height of ground midway between last and current range step.
c           For use when shifting profiles to be relative to the local ground
c           height. 
c   KBIN = number of bins that the field is shifted

c PARINIT: 
c  RV2 = range of the next refractivity profile (for range-dependent cases)
c REFDUM() = array containing M-unit values for current (interpolated) profile.
c              When using BS method this represents refractivity profile with 
c              respect to the local ground height.
c  HTDUM() = array containing height values for current (interpolated) profile.
c            When using BS method this represents height values with respect to
c             the local ground height.
c   PROFINT() = M-unit profile interpolated to every DELZ in height 
c   HT() = height array of size N
c   IS = counter for current profile (for range-dependent cases)

c PATTERN:
c   PELEV = sine of elevation angle
c   AFAC = constant used in determining antenna pattern factors
c          AFAC = 1.39157 / sin( bw / 2 ) for SIN(X)/X and height-finder
c          AFAC = (.5*ln(2))/(sin(bw/2))**2 for GAUSSIAN 
c   BW = antenna pattern beamwidth in radians 
c   ELV = antenna pattern elevation angle in radians
c   UMAX = limiting angle used in 30 dB cut-off point for SIN(X)/X and 
c          generic height-finder antenna pattern factors
c   SBW = sine of the beamwidth

c PEVAR:
c   FKO = free-space wavenumber = WL / (2*pi)
c   WL = Wavelength in meters
c   DELZ = Bin width in z-space = WL / (2*sin(THETAMAX))
c   N = Transform size
c   LN = Power of 2 transform size, i.e. N = 2**LN 


c RHSTPS:
c   DR = PE range step in meters
c   DR2 = 1/2 PE range step in meters
c   RDOUT = Output range step in meters
c   DZOUT = output height increment in meters
c   ZOUT() = array containing all output height points


c TRVAR:
c   DMDH() = gradients of first profile in M-units/meters
c   ZLIM = height limit for ray trace
c   JLS = index of refractivity array at which antenna height is located
c   THETALAUNCH = angle in radians of launch angle for which, when traced,
c                 height of the ray at each output range step is stored.
c   RLIM = 90% of maximum range RMAX - used for ray tracing
      integer flagpu
      common /flagpu/flagpu


      common / arrays / X(0:maxpts), Y(0:maxpts), filt(0:maxn4),   
     +                  SK(0:maxpts), CK(0:maxpts), CM(0:maxpts),
     +                  SM(0:maxpts), xlst(0:maxpts), ylst(0:maxpts)
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      real*8  dr, drout,dzout, dr2
      real*8  rout
      common / irhstps / izout(mxzout)
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1,
     +                   antref, rpe, hlim(mxrout), slp(mxter), fter
      common / pattern / pelev, afac, bw, elv, umax, sbw
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is
      common / htvar / ylast, ycur, ycurm, kbin,iylast,iycur 
      common / trvar / dmdh(mxlvls), zlim, jls, thetalaunch, rlim 
                                      
      type(errorflag) ef
      type(inputvar) vnp
      type(refractivity) rf
      type(systemvar) sv
      type(terrain) tr
      type(refparam) rp
                                     
      dimension href(mxlvls), refref(mxlvls)

      logical fter, lopen, lang, lflag
      
      data radc / 1.74533e-2 /      !degree to radian conversion factor
      data iscr / 20 /                 ! Unit number for scratch file
      data hdeg / 8.726646e-3 /        ! 1/2 degree in radians
      data c0 / 299.79245 /         !speed of light x 1e.-6 m/s
      data sdeg20 / .3420201433 /       ! Sine of 20 degrees
      data sdeg15 / .258819045 /       ! Sine of 15 degrees
      


      iperror=0
      ierror = 0
      fter = .false.
      angu = 0.
      rangu = 0.
      thetamax = 0.
      hminter = 0.
      antref = sv%antht                  
c pg check dimensions:
      if (rf%lvlep .gt. mxlvls) then
         write(*,*)'rf%lvlep,mxlvls=',rf%lvlep,mxlvls
         stop 'dimension of mxlvls is too small'
      endif
      
c Put lower limit on HMAX and RMAX
      
c    vnp.rmax = amax1( vnp.rmax, 5000. ) !Set max. range to no less than 5 km.
c    vnp.hmax = amax1( vnp.hmax, 100.)  !Set max. height to no less than 100m.
      
c      dzout = vnp.hmax / float( vnp.nzout )
      dzout=vnp%dz
c      drout = vnp.rmax / float( vnp.nrout )
      drout = vnp%dr 
      
      WL = c0 / sv%freq
      FKo = 2. * pi / WL
      con = 1.e-6 * fko
      pi2=pi/2. 

c  Calculate constants used to determine antenna pattern factor
c SV.IPAT = 0 -> omni
c SV.IPAT = 1 -> gaussian
c SV.IPAT = 2 -> sinc x
c SV.IPAT = 3 -> csc**2 x
c SV.IPAT = 4 -> generic height-finder      
      
      sv%bwidth = amax1( sv%bwidth, .05 )   !Don't let vertical beamwidth go
      sv%bwidth = amin1( sv%bwidth, 45. )  !outside .5 to 45 degree limit.
      
      sv%elev = amax1( sv%elev, -10. ) !Don't let elevation angle go
      sv%elev = amin1( sv%elev, 10. )  !outside -10 to 10 degree limit.
      
      bw = sv%bwidth * radc
      elv = sv%elev * radc
      bw2 = .5 * bw
      if( sv%ipat .eq. 1 ) then
         afac = .34657359 / (sin( bw2 ))**2
         pelev = sin( elv )
      elseif( sv%ipat .eq. 3 ) then
         sbw = sin( bw )
      elseif( sv%ipat .ne. 0 ) then
         afac = 1.39157 / sin( bw2 )
         a = pi / afac
         umax = atan( a / sqrt(1. - a**2)  )
      end if

c Discard any unnecessary terrain points.  Test on the rate of change of slope,
c i.e., second derivative.  If that is below 1.e-3 then discard that point.
      
      if( tr%itp .gt. 0 ) fter = .true. 
      if( fter ) then             
         
c Check that all terrain range points are increasing.

         do i = 1, tr%itp-1
            ip1 = i + 1                      
            if( tr%terx(ip1) .lt. tr%terx(i) ) then
               ierror = -17
               return
            end if
         end do

c Test to see that first range value is 0.         

         if( tr%terx(1) .gt. 0. ) then
            ierror = -18
            return
         end if
         
c Test to see if the last range point in the profile meets or exceeds RMAX.  If
c not then return error code.
   
         if( tr%terx(tr%itp) .lt. vnp%rmax ) then
            if( ef%lerr6 ) then
               ierror = -6
               return
            else
               tr%itp = min(tr%itp + 1, mxter)
               tr%terx(tr%itp) = vnp%rmax*1.01
               tr%tery(tr%itp) = tr%tery(tr%itp - 1)
            end if
         end if
                
c Test to see if the unit number for the scratch file is already attached to 
c another file.  If so, search for a unit number that is unattached.

         inquire( iscr, OPENED = lopen )  
         do while (lopen)
            iscr = iscr + 1
            inquire( iscr, OPENED = lopen )
         end do

         open( iscr, status = 'SCRATCH')
         write(iscr,*) tr%terx(1), tr%tery(1) ! Keep 1st point of terrain prof
         do i = 2, tr%itp-1 
            im1 = i - 1
            ip1 = i + 1
            xm1 = tr%terx(im1)
            ym1 = tr%tery(im1)
            xi = tr%terx(i)
            yi = tr%tery(i)
            xp1 = tr%terx(ip1)
            yp1 = tr%tery(ip1)
            
            dx1 = amax1( 1.e-3, xi - xm1 )
            dx2 = amax1( 1.e-3, xp1 - xi )
            
            sl1 = (yi - ym1) / dx1
            sl2 = (yp1 - yi) / dx2
         
            scd = sl2 - sl1          ! dx is taken to be 1 m

c If second derivative is large enough then keep this point. 
            
            if( abs(scd) .GT. 1.e-3 ) write(iscr,*) xi, yi    
         end do
      
         write(iscr,*) tr%terx(tr%itp), tr%tery(tr%itp) !Keep last point in
         rewind( iscr )

c Now the scratch file contains all the necessary points for the terrain prof
c Go back and read them in the arrays TR.TERX, TR.TERY.
         
         bugfix = 0.
         lflag =.false.                              !eof(iscr)   
         tr%itp = 0
         do while( .not. lflag )                       
            tr%itp = tr%itp + 1
            read(iscr,*,end=999) tr%terx(tr%itp), tr%tery(tr%itp)
            if( tr%terx(tr%itp) .ge. vnp%rmax )  goto 999    !exit
            bugfix = 0.
c            lflag = eof(iscr)
         end do
999      continue
      
         close(iscr)
                    
c Determine minimum height of terrain profile.  Then adjust entire terrain prof
c by this minimum height HMINTER such that this is the new 0 reference.

c         hminter = vnp.hmax
c          do i = 1, tr.itp 
c            yi = tr.tery(i)
c            if( yi .lt. hminter ) hminter = yi
c         end do

             hminter = 0

         antref = sv%antht + tr%tery(1) 
         ipeak = 0
         ipeaklst = 0
         lang = .true.
         tamax = 0.
          
         do i = 1, tr%itp-1   
            
            ip1 = i + 1
            y1 = tr%tery(i)
            x1 = tr%terx(i)
            y2 = tr%tery(ip1)
            x2 = tr%terx(ip1) 
            
            xdif = x2 - x1
            ydif = y2 - y1
            xdif = amax1( xdif, 1.e-5 )
            slope = ydif / xdif                          

            slp(i) = slope
                      
c Calculate angle from each terrain point height to antenna height above
c reference(HMINTER).Determine maximum propagation angle so that direct ray
c angle will clear highest peak. 

            if( y1 .gt. antref ) then
               angle = atan( (y1-antref) / x1 )
               if( angle .gt. angu ) then
                  angu = angle
                  rangu = x1
               end if
            end if
            
         end do 
         
c  Add 1/2 degree to the angle that clears the highest peak.
         
         angu = angu + hdeg
         slp( tr%itp ) = 0.
         
      end if
       
      zlim = amax1( real(vnp%hmax), antref )

c  introduced by pg  24/1 97:
       if (itpem_opt(5).eq.1)  then
          do j=1, rf%nprof    
             zlim = amax1( zlim, rp%maxref(j))
          enddo
          is = 1
          rv2=rf%rngprof(is)
       else       
          do j=1, rf%nprof    
             zlim = amax1( zlim,      rf%hmsl(rf%lvlep,j))
          enddo
c Initialize refractivity arrays.
      
          call refinit( ef%lerr12, vnp%rmax, rf, IERROR ) 
          if( ierror .ne. 0 ) return

c Calculate gradients and other variables for use in ray tracing.
c pg raytracing is not used...       
cpg      if( rf.nprof .gt. 1 ) call remdup( REFDUM, HTDUM, RF.LVLEP)   
      endif

 
C*******************************
c      write(*,*)'dr=dr',dr
      if (dr.eq.0) then 

c---- pg changed experimentally
         thetamax=40*radc
         thetamax=7*radc
         call getfftsz( fter, thetamax, wl, mxnfft, zlim,  dzout,
     +        DELZ, ZMAX, LN,  N, IERROR )
c If error then exit routine.
         if( ierror .ne. 0 ) then
            write(*,*) 'Error in determining FFT size'
            stop
         endif
      
      if( fter ) then
c Use 74% of ZMAX instead of 75% to leave some slop and ensure the 
c FFT size is not surpassed.
      
         if( .74*zmax .gt. zlim ) then
            thetafrac = thetalaunch / thetamax
            zmax = zlim / .74      
            sthetamax = float(n) * wl * .5 / zmax 

            sthetamax = amin1( sthetamax, .9 ) 
            delz = wl * .5 / sthetamax           
            ifrac=(int(dzout/delz)+1)
            delz=dzout/ifrac
            do while (delz.lt.(wl/2)) 
               ifrac=ifrac-1
               delz=dzout/ifrac
c     write(*,*)'ifrac,delz',ifrac,delz
            enddo
            zmax = float(n) * delz 
         end if
      end if ! fter-loop
 
c Determine horizon range based on transmitter height and 0 receiver height
c by RHOR = 3572. * sqrt( 1.3333 * antref )                           
                     
c      rhor = 4124.5387 * sqrt( sv.antht )
      dr = 2. * fko * delz**2  
      rkm = vnp%rmax * 1.e-3
      
c Determine range step.
      
      if( fter ) then
          dr = amin1( real(dr), 700. ) 
          dr = amin1( real(dr),50*wl)
c          dr =  10*wl
      else   
          dr = amin1( real(dr), 700. ) 
          dr = amin1( 1000., 400*wl)
      end if   ! fter-loop 
      if (real(vnp%hmax).lt.100) then
c pg 22 sep 98: emirically a small antenna height requires smaller step !
c         dr=dr*vnp.hmax/100
         if (flagpu.lt.0) then
            write(*,*)' The antenna heigth is less than 100 m. '
            write(*,*)' This causes problems for low field values.'
            write(*,*)'  To get a stable result dr migth have to be'
            write(*,*)'  reduced by a factor 100. It is not done.'
         endif
      endif
      else
         if (flagpu.lt.5) then
            write(*,*)'flagpu',flagpu
c           write(*,*)' 2nd time dz is not changed '
            write(*,*)' User specified range and depth steps'
         endif
         if (delz.lt.(wl/2)) then
            write(*,*)'**************************'
            Write(*,*)'Warning deltaZ less than lambda/2'
            write(*,*)'delz,lambda/2=',delz,wl/2 
            pause
        endif
        n=2**ln 
      endif ! dr-loop

c-added by pg to obtain an identical steps
      delz=dzout/(int(dzout/delz+0.99))
      zmax = float(n) * delz 
      dr=drout/dfloat(int(drout/dr+0.99))
      if (flagpu.lt.0) then
         write(*,*)' dr,delz,n,zmax,zlim,wl',dr,delz,n,zmax,zlim,wl
         write(*,*)' drout/dr',drout/dr
      endif

      dr2 = .5 * dr

c path loss constant, add to TL to get PL:

      plcnst=20.*alog10(2.*fko)

c Setup output height array
      
      do i = 1, vnp%nzout
         zout(i) = vnp%hmin+float(i-1) * dzout
         izout(i)=nint(zout(i)/delz)
      end do
      
c Initialize variables for free-space propagator phase calculations.

      delp = pi/zmax
      FNorm = 2. / N
      cnst = delp / fko
      nm1 = n - 1
      
c Initialize variables and set-up filter array.

      no4 = n/4
      n34 = 3.* no4
      cn75 = 4.* pi / N
      do i = 0, no4
         fj= cn75 * float(i)
         filt(i) = .5 + .5 * cos(fj)
      end do
                       
c Initialize starter field.
                       
      call xyinit( sv )
                           
c Transform to z-space.
       
      call sinfft( ln, X )
      call sinfft( ln, Y )
c      do i=1,n-1
c         write(10,*)I*delz,20*log10(sqrt(x(i)**2+y(i)**2))
c      enddo
      ylast = 0.
      if( fter ) ylast = tr%tery(1) 
      iylast=nint(ylast/delz)            
      ycur=0
      if( fter ) ycur = tr%tery(1) 
      iycur=nint(ycur/delz)            
      kbin = 0
      ycurm = 0.
      rout = 0.

c Define mesh array in height

      do i=0,n
         ht(i)= float(i)*delz
      end do

c Trace THETALAUNCH ray and store all heights at each output range step in
c array HLIM().
      
c      if( fter ) then
c         if(( thetalaunch .gt. angu ) .and. ( angu .gt. 0. ) .and.  
c     +       (rangu .le. 10000.)) thetalaunch = -thetalaunch
c      end if
c      call traceh( vnp.nrout, rf.lvlep, tr )        
      
c Determine the free-space propagator (p-space) arrays.

      call phase1    

c If smooth surface and range-independent case then initialize all refractivity
c and z-space propagator arrays now.
                       
      if(( .not. fter ) .and. (rf%nprof .eq. 1 ) 
     1     .and. (itpem_opt(6).eq.0)) then
         if (itpem_opt(5).eq.1)  then
          call paramprofile(rp%base(1),rp%thick(1),
     1           rp%offset(1),rp%mdef(1),rp%coef(1,1),sv%antht,iperror)
        else
          call profref( hminter, RF%LVLEP, HREF, REFREF, 0 ) 
          call intprof( REFREF, HREF, RF%LVLEP )
         endif
         call phase2 
      end if
c         do  i=0,N
c           write(7,*)i,profint(i),profint(i)/con
c        enddo
c
c
c      write(*,*)'exit init'
      end 

c ************************** SUBROUTINE GETFFTSZ *************************** 

c Purpose:  Determines the FFT size needed for the given problem.  

c Called from: PESTEP

c Routines called: NONE
       
      subroutine getfftsz( fter, thetamax, wl, mxnfft, zlim, dzout,
     +             DELZ, ZMAX, LN, N, IERROR )
      real*8  dzout 
      logical fter
      
      ierror = 0 
      sthetamax = sin( thetamax )
       delz =  wl * .5 / sthetamax              
c changed by pg           
       delz=dzout/(int(dzout/delz)+1)
                   
                 
c Set lower FFT limit to 2**9 for smooth surface cases, if terrain case then
c set lower FFT limit to 2**10.
                 
      ln = 9
      if( fter ) ln=10
c      ln=11
      N=2**LN

      zmax=delz*float(n)         

c Determine transform size needed to perform calculations to a height of ZLIM, 
c up to the maximum FFT size allowed.
      
      do while( .75*zmax .lt. zlim ) 
         ln = ln + 1
         n = 2**ln
         zmax = delz * float(n)
      end do
c If the transform size needed is too large then return error code.
c      write(*,*)'ln,mxnfft,zlim,zmax,delz',ln,mxnfft,zlim,zmax,delz
      if( ln .gt. mxnfft ) then 
         write(*,*)' FFT array not large enough'
         write(*,*)' requested size (ln)=', ln
         write(*,*)' max size (mxnfft)=',mxnfft
         ierror = -10
      endif
      end
       
c ************************** SUBROUTINE INTPROF ********************************

c Purpose:  Performs a linear interpolation vertically with height on the 
c           refractivity profile.  Stores interpolated profile in PROFINT().

c Called from: REFINTER, PESTEP

c Routines called: NONE

c Parameter list:

c REF() = Current NLEV refractivity profile
c HIT() = Heights of NLEV refractivity profile
c NLEV = Number of levels in profile
      
      SUBROUTINE intprof( REF, HIT, NLEV )

      include 'tpem.inc'

      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is

      dimension ref(*), hit(*)
      
      J=2

      DO I=0,N
         height = ht(i)
   40    IF((height .LE. hit(J)) .OR. (J .GE. Nlev)) then
            k = j - 1
            FRAC = (height - hit(k)) / (hit(J) - hit(k))
            profint(I) = (ref(k) + FRAC * (ref(J) - ref(k))) * con
         else
            J=J+1
            GO TO 40
         end if 
      end do

      END
      
c ***************************** SUBROUTINE PHASE1 ****************************

c Purpose: Initialize free-space propagator factors SK and CK for the high
c          angle propagator.

c Called from: PEINIT

c Routines called: NONE

      SUBROUTINE PHASE1

      include 'tpem.inc'

      common / arrays / X(0:maxpts), Y(0:maxpts), filt(0:maxn4),     
     +                  SK(0:maxpts), CK(0:maxpts), CM(0:maxpts),
     +                  SM(0:maxpts), xlst(0:maxpts), ylst(0:maxpts)
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1,
     +                   antref, rpe, hlim(mxrout), slp(mxter), fter
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
                                                                    
      logical fter
                                                                    
      double precision cak
      
      drfk = dr * fko
      DO I=0,N
         ak = float(i) * cnst
         aksq=ak * ak
         aksq = amin1( 1., aksq )
         cak = sqrt(1. - aksq)
         ang = drfk * ( 1.d0 - cak )
         SK(I) = SIN(ANG) *fnorm
         CK(I)=  COS(ANG) *fnorm
      end do

c Filter the upper 1/4 of the propagator arrays.

      do  i = n34, n
         attn = filt(i-n34)
         sk(i)=attn*sk(i)
         ck(i)=attn*ck(i)
      end do

      END  

c *************************** SUBROUTINE PHASE2 ******************************

c Purpose: Calculates the environmental phase factors in z-space.

c Called from: PEINIT, PESTEP 

c Routines called: NONE
           
      SUBROUTINE PHASE2
 
      include 'tpem.inc'

      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      common / arrays / X(0:maxpts), Y(0:maxpts), filt(0:maxn4),          
     +                  SK(0:maxpts), CK(0:maxpts), CM(0:maxpts),
     +                  SM(0:maxpts), xlst(0:maxpts), ylst(0:maxpts)
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is

      do i = 0, n                   
         ang = dr * profint(i)
         SM(I)= SIN(ANG)
         CM(I) = COS(ANG)
      end do

c Filter upper 1/4 of the arrays.

      do i = n34, n
         attn = filt(i-n34)
         sm(i) = attn * sm(i)
         cm(i) = attn * cm(i)
      end do

      END   
      
c ****************************** SUBROUTINE PROFREF *************************
            
c Purpose: This subroutine determines the refractivity profile with respect 
c          to the reference height YREF which, depending on the value of IFLAG,
c          can be HMINTER or the local ground height above HMINTER.

c Called from:  PESTEP, REFINTER

c Routines called: NONE

c Parameter list:

c YREF = Reference height in meters. 
c NLVL = Number of levels in profile
c H() = Heights of refractivity profile
c RM() = Refractivity array
c IFLAG = 0: Profile arrays RM() and H() will be referenced to height HMINTER, 
c            and will also be used to initialize REFDUM() and HTDUM().
c       = 1: Profile arrays RM() and H() will be referenced to the local ground
c            height.
            
      subroutine profref( yref, NLVL, H, RM, iflag )
      
      include 'tpem.inc' 
      
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is

      dimension h(mxlvls), rm(mxlvls) 
             
      if( yref .gt. 0. ) then       
         do i = 1, mxlvls
            h(i) = 0.
            rm(i) = 0.
         end do                                           
         
c Get refractivity profile level at which the height of the ground is 
c just above. This level is JS.

         js = 0 
         nlvlm1 = nlvl - 1
         do i = 1, nlvlm1
            if(( yref .le. htdum(i+1) ) .and. ( yref .gt. htdum(i) ))
     +         js = i
         end do
                        
c Determine the refractivity value at the ground and fill arrays H and RM
c with refractivity profile where height 0. now refers to the ground reference,
c either local ground height or HMINTER.
                        
         if( js .ne. 0 ) then
            frac = (yref - htdum(js))/(htdum(js+1) - htdum(js))
            rmu = refdum(js) + frac * (refdum(js+1) - refdum(js))
            if( int( frac ) .eq. 1 ) js = js + 1
            newl = nlvl - js +1
            rm(1) = rmu
            h(1) = 0.
            k = js + 1
            do jk = 2, newl
               rm(jk) = refdum(k)
               h(jk) = htdum(k) - yref
               k = k + 1
            end do 
            nlvl = newl                           
            if( iflag .eq. 0 ) then
               do i = 1, mxlvls
                  refdum(i) = rm(i)
                  htdum(i) = h(i)
               end do
            end if
         end if
      else  

c If the reference height is 0. then height/refractivity profile is unchanged.
      
         do i = 1, nlvl
            h(i) = htdum(i)
            rm(i) = refdum(i)
         end do
      end if   
                  
      end
      
c ************************* SUBROUTINE REFINIT *****************************

c Purpose: Initializes refractivity arrays used for subsequent PE calculations.

c Called from: PEINIT

c Routines called: REMDUP
      subroutine refinit( elerr12, vrmax, rf, IERROR )
      
      include 'tpem.inc'
      
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is
      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      
      type(refractivity) rf
      
      logical elerr12
      
      data hlarge/ 1.e4 /
      data rlarge / 1.e10 /
      
      ierror = 0
                            
c Test to see if last profile entered ( for range dependent case ) meets or
c exceeds RMAX, otherwise, return error (unless error trapping is turned off - 
c EF.LERR12 = .FALSE.).
       
      if( rf%nprof .gt. 1 ) then
         if(( rf%rngprof(rf%nprof) .lt. vrmax ) .and. ( elerr12 )) then
            ierror = -12 
            return
         end if
      end if   
                  
c  Add extra level to tabulated prof with extrapolated gradient.  Test on HDIF
c greater than 0 for profiles that contain multiple height/M-unit 
c values that are equal. 
      
      rf%lvlep = rf%lvlep + 1
      do i = 1,rf%nprof 
         hdif = 0. 
         lvlm1 = rf%lvlep
         lvlm2 = rf%lvlep
         do while( hdif .eq. 0. )          
            lvlm1 = lvlm1 - 1
            lvlm2 = lvlm1 - 1
            hdif = rf%hmsl(lvlm1,i) - rf%hmsl(lvlm2,i)
         end do
         grad = (rf%refmsl(lvlm1,i)-rf%refmsl(lvlm2,i)) / hdif

c If last gradient in refractivity profile is negative then return error.
     
         if( grad .lt. 0 ) then
            write(*,*) 'STOPPING: Last gradient in any refract', 
     +                 'profile entered is negative'
            write(*,*) 'Range',i
            write(*,*) 'points: lvlm1,lvlm2',lvlm1,lvlm2  
            write(*,*)'refractivity values:',
     +           rf%refmsl(lvlm1,i),rf%refmsl(lvlm2,i)
            write(*,*)'coresponding heigth values:',
     +           rf%hmsl(lvlm1,i),rf%hmsl(lvlm2,i)
            ierror = -14
            return
         end if                                        

         rf%hmsl(rf%lvlep, i) = hlarge
         rf%refmsl( rf%lvlep, i ) = (hlarge-rf%hmsl(lvlm1,i)) * grad +
     +                   rf%refmsl( lvlm1, i )
      end do
      
      is = 1
      rv2=rf%rngprof(is)

      do i = 1, rf%lvlep
         refdum(i) = rf%refmsl( i, is )
         htdum(i) = rf%hmsl( i, is )
      end do
      
      if( rf%nprof .eq. 1 ) call remdup( REFDUM, HTDUM, RF%LVLEP )   

      np = rf%nprof + 1
      rf%rngprof(np) = rlarge
      do i = 1, rf%lvlep  
         npm1 = np - 1
         rf%hmsl( i, np ) = rf%hmsl( i, npm1 )
         rf%refmsl( i, np ) = rf%refmsl( i, npm1 )
      end do

      end
      
c **************************** SUBROUTINE REFINTER ************************
                
c Purpose: Interpolates vertically and horizontally on the refractivity profs.

c Called from: PESTEP

c Routines called: REMDUP, PROFREF, INTPROF

c Parameter list:

c RANGE: Range for profile interpolation
c HMINTER: Reference height in meters

      subroutine refinter( rf, range, hminter )

      include 'tpem.inc'

      common / htvar / ylast, ycur, ycurm, kbin,iylast,iycur
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is
       
      type(refractivity) rf
       
      dimension href(mxlvls), refref(mxlvls)

      save j, rv1
      data j, rv1 / 0, 0. /

c One-line interpolation function

      pint( p1, p2 ) = p1 + fv * ( p2 - p1 ) 
      
      nlvl = rf%lvlep
      
c If there is a range-dependent refractivity prof then interpolate horizontally
c using the two surrounding profiles at range RANGE with all duplicate levels.
  
      if( rf%nprof .gt. 1 ) then
         IF( range .gt. rv2 ) then
            j = is
            IS=IS+1
            rv1=rv2
            rv2=rf%rngprof(IS)
         end if
      
         FV=(range-rv1)/(rv2-rv1)

         do i = 1, nlvl
            refdum(i) = pint( rf%refmsl(i,j), rf%refmsl(i,is) )
            htdum(i) = pint( rf%hmsl(i,j), rf%hmsl(i,is) )
         end do                                          
            
c Now remove all duplicate levels with NLVL now being the # of points in the
c profile at range RANGE.
            
         call remdup( REFDUM, HTDUM, NLVL )
         call profref( hminter, NLVL, HREF, REFREF, 0 )
         
c At this point REFDUM() and HTDUM(), also HREF() and REFREF(),  are referenced
c to HMINTER.
         
      end if   
                          
c Using BS method must
c determine height and M-unit profiles relative to ground, where YCURM is now 
c the height of the local ground above the reference height HMINTER.
 
      call profref( ycurm, NLVL, HREF, REFREF, 1 )
      
c Interpolate vertically with height.  PROFINT is now an N-point (N=2**NFFT)
c array containing the interpolated M-unit values for the refractivity at 
c range RANGE.
 
      call intprof( REFREF, HREF, NLVL )
      
      end    

c ************************* SUBROUTINE REMDUP *****************************

c Purpose: Removes duplicate refractivity levels in profile.

c Called from: REFINIT, REFINTER

c Routines called:  NONE
 
      subroutine remdup( REF, HYT, NLVL )
 
      dimension ref(*), hyt(*)

c  Remove all duplicate levels in interpolated profile
      
      i = 1
      do while( i .lt. nlvl )
         ht1 = hyt(i)
         ht2 = hyt(i+1)  
         if( abs(ht1-ht2) .le. 1.e-3 ) then
            nlvl = nlvl - 1
            do j = i, nlvl 
               jp1 = j + 1
               hyt(j) = hyt(jp1)
               ref(j) = ref(jp1)
            end do
            i = i - 1
         end if
         i = i + 1
      end do

      end  
      
c *************************** SUBROUTINE TRACEA **************************** 

c Purpose:  This routine performs a ray trace to determine the minimum angle 
c           required (based on the reflected ray) to obtain a PE solution for 
c           all heights up to HMAX and all ranges beyond 90% of RMAX. THETAMAX 
c           is then determined from this angle.  This is done only for smooth
c           surface.  For terrain cases, the minimum angle to reach HMAX and 
c           90% RMAX is also based on the reflected ray. However, the reflected
c           ray angle is based on a smooth surface bounce.  A true terrain ray
c           trace is not done - the reflected ray based on a smooth surface
c           bounce is just to get "near the ball park" of  THETAMAX.
c           From here THETAMAX is 
c           maximized within the FFT size determined from the "minimum" direct 
c           ray.  This is done in order to minimize DELZ so subsequent shifting
c           of the X and Y arrays will not produce crude coverage diagrams 
c           (note: this doesn't always work).

c           If PRANG is not equal to 0, then the user has overriden the default
c         calculation and THETAMAX is then determined based on PRANG.  However 
c           a ray trace must still be done in order to determine the initial
c           launch angle such that the local angle of the ray remains less than
c           PRANG. The initial launch angle is used in subroutine TRACEH.

c Called from: PEINIT

c Routines called: NONE

c Parameter list:
 
c PRANG:  Problem angle (rad) - input by the user          

c Returns:
c THETAMAX      (in common block MISCVAR)
c RPE           (in common block MISCVAR)
c THETALAUNCH   (in common block TRVAR)

      subroutine tracea( prang, rlvlep, tr )
      
      include 'tpem.inc'
      
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1,
     +                  antref, rpe, hlim(mxrout), slp(mxter), fter
      common / trvar / dmdh(mxlvls), zlim, jls, thetalaunch, rlim
      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is
      
      type(terrain) tr
      
      logical fter, loop
       
      integer rlvlep
              
      data deg15 / .2617994 /

c All heights and ranges are in meters, gradients are in M-unit/meter * 1.e-6 
c and angles are in radians.

c Define in line ray trace functions:

      rada1( a, b ) = a**2 + 2. * grad * b             !a=a0, b=h1-h0    
      rp( a, b ) = a + b / grad                        !a=r0, b=a1-a0
      ap( a, b ) = a + b * grad                        !a=a0, b=r1-r0
      hp( a, b, c ) = a + ( b**2 - c**2 ) / 2. / grad  !a=h0, b=a1, c=a0

c HGRND = function to determine height of the ground at a given range.
      
      hgrnd( h, qm, rdif ) = amax1( 0., h + qm * rdif ) 
      
c AS = Starting launch angle in radians.
c ASL = Last (or previous) starting launch angle.
c AMXCUR = Maximum of local angle along ray.
c AMXCURL = Last (or previous) AMXCUR
c ISET = flag to test whether or not to stop loop

      as = -thetamax
      asl = 0.
      amxcurl = 0.    
      idn = -1
      iset = 0  

      do while( iset .eq. 0 ) 
                                            
c Increase or decrease angle by 1 mrad depending on the value of IDN and AS.
c Initialize ray trace variables.

         as = as + idn * 1.e-3         
         h0 = antref
         r0 = 0. 
         rpe = 0.        
         a0 = as
         jl = jls 
         amxcur = 0.
         it = 1 
         loop = .true.

c Perform ray trace until ray has reached ZLIM and/or RLIM where  
c ZLIM = maximum of HMAX or ANTREF
c RLIM = .9 * RMAX  
       
         do while( loop ) 
    
            grad = dmdh(jl) 
            if( .not. fter ) then
               if( a0 .lt. 0. ) h1 = htdum(jl)
               if( a0 .gt. 0. ) h1 = htdum(jl + 1)
               if( a0 .eq. 0. ) then
                  if( grad .lt. 0. ) h1 = htdum(jl)
                  if( grad .gt. 0. ) h1 = htdum(jl+1)
               end if
               if( h1 .gt. zlim ) h1 = zlim
               rad = rada1( a0, h1-h0 ) 
               if( rad .gt. 0 ) then
                  a1 = sign( 1., a0 ) * sqrt( rad )
               else
                  a1 = 0.   
                  h1 = hp( h0, a1, a0 )
               end if

               r1 = rp( r0, a1-a0 ) 
           
            else
               r1 = tr%terx(it+1)
               a1 = ap( a0, r1-r0 )
               h1 = hp( h0, a1, a0 )
               hter = hgrnd( tr%tery(it), slp(it), r1-tr%terx(it) )
               if(( h1 .lt. hter ) .and. ( rpe .gt. 0. )) goto 948 !exit
               IF( h1 .lt. hter ) THEN
                   
c Determine slope of line segment going into ground

                  sl2 = (h1 - h0) / (r1 - r0)
         
c Determine range at which ray hits ground and the height of the ground
c at that range.

                  slpdif = sl2 - slp(it)
       
                  r1 = (tr%tery(it)-h0+sl2*r0-slp(it)*tr%terx(it)) / 
     +                  slpdif
                  hter = hgrnd( tr%tery(it), slp(it), r1-tr%terx(it) )
                  h1 = hter
                  a1 = ap( a0, r1-r0 )
                  ihg = 1
               END IF
            end if
               
            if(( a1 .le. 0. ) .and. ( h1 .le. htdum(jl) )) then 
               h1 = htdum(jl) 
               if( h1 .gt. hter ) ihg = 0
               rad = rada1( a0, h1-h0 ) 
               a1 = -sqrt( rad )
               r1 = rp( r0, a1-a0 )    
               jl = jl - 1
               if( jl .eq. 0 ) jl = 1    
            elseif(( a1 .ge. 0. ).and.( h1 .ge. htdum(jl+1) ))then
               h1 = htdum(jl+1)
               rad = rada1( a0, h1-h0 ) 
               a1 = sqrt( rad )
               r1 = rp( r0, a1-a0 )
               jl = jl + 1
               if( jl .gt. rlvlep ) jl = rlvlep
            end if 
 
            if( h1 .gt. zlim ) then
               h1 = zlim
               rad = rada1( a0, h1-h0 ) 
               a1 = sqrt( rad )
               r1 = rp( r0, a1-a0 )
            end if
                     
            h0 = h1
            r0 = r1
            a0 = a1  

c If smooth surf case set RPE to range at which reflected ray hits the ground.
            
            if( .not. fter ) then
               if( h0 .eq. 0. ) then 
                  a0 = -a0 
                  rpe = r0      
               end if
            else
               if( ihg .eq. 1 ) then 
                  ihg = 0
                  rpe = r0
                  a0 = -a0
               end if 
               if( r1 .ge. tr%terx(it+1) ) it = it + 1 
               
            end if
            if( a0 .ge. 1.57079 ) goto  948 !exit
            amxcur = amax1( amxcur, a0 )   
            if(( h0 .ge. zlim ) .and. ( a0 .gt. 0.)) loop = .false.
            if( r0 .gt. rlim ) loop = .false.
         end do 
 948    continue

c
c
c   pg:  we set rpe=0

        rpe=0

             

c Test to see if the current ray traced from launch angle AS meets criteria.  
c If ray traced does not reach ZLIM AND is not within RLIM the initial launch
c angle AS is increased (smooth surface) by 1 mrad and ray trace is repeated.
         
         if( fter ) then
            if(( h0 .ge. zlim ).and.( r0 .le. rlim )) iset = 1
         else
            if(( r0 .le. rlim ) .and. ( rpe .gt. 0. )) iset = 1   
         end if   

c If criteria is met then (if user specified an angle) make sure the local
c maximum angle is just within PRANG.  If not then repeat ray trace until this 
c occurs.
         
         if(( prang .gt. 0. ) .and. ( iset .eq. 1 )) then
            a = amax1( abs(as), amxcur )
            if( a .lt. prang ) then
               iset = 0                
            elseif( asl .ne. 0. ) then
               as = asl
               amxcur = amxcurl
            end if
         end if   

c Just as a safeguard - set absolute maximum of launch angle to 15 degrees.
                      
         if( as .le. -deg15 ) then
            iset = 1
            as = -deg15
            amxcur = deg15 
         end if
         
         asl = as
         amxcurl = amxcur
         
      end do 

      thetamax = amax1( abs(as), amxcur ) 
      thetalaunch = abs(as)
      
      end 
      
c ************************ SUBROUTINE TRACEH ********************************

c Purpose: Computes ray trace for a single ray with launch angle -THETALAUNCH.
c          Upon reflection the heights of this ray at each output range point 
c       RO is then stored for subsequent output of loss values in array MLOSS. 
c       This is done so that only loss values that fall within the valid PE
c       solution region are output or passed back in MLOSS.  For terrain cases,
c      a true terrain ray trace is not performed (the reflected angle is based 
c      on a smooth surface reflection) but usually THETALAUNCH for this case  
c      is so steep that one can get away with it for surface based transmitter 
c      heights.

c Called from: PEINIT 

c Routines called: NONE
 
c Returns:
c HLIM =array contains heights of the ray traced to each output range point RO.
c          For ranges less than RPE, HLIM() = 0. (In common block MISVAR)

      subroutine traceh( vnrout, rlvlep, tr ) 
      
      include 'tpem.inc'
      
      common / trvar / dmdh(mxlvls), zlim, jls, thetalaunch, rlim
      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2, zout(mxzout) 
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1, 
     +                   antref, rpe, hlim(mxrout), slp(mxter), fter
      common / parinit / rv2, refdum(mxlvls), htdum(mxlvls),
     +  profint(0:maxpts), ht(0:maxpts), is

      type(terrain) tr
      
      logical fter
      
      integer vnrout, rlvlep
      
c Define one-line ray trace functions:
      
      rada1( a, b ) = a**2 + 2. * grad * b      
      rp( a, b ) = a + b / grad 
      ap( a, b ) = a + b * grad 
      hp( a, b, c ) = a + ( b**2 - c**2 ) / 2. / grad
      
c HGRND = function to determine height of the ground at a given range.
      
      hgrnd( h, qm, rdif ) = amax1( 0., h + qm * rdif )
      
      a0 = -thetalaunch 
      h0 = antref
      jl = jls
      ro = drout
      ihu = 0  
      ihl = 0
      r0 = 0.
      rpe = 0.   
      ihg = 0
      if( fter ) it = 1
      
c Ray is traced through NROUT output range points.
      
      do i = 1, vnrout
 
c Trace until ray reaches output range point RO.
                              
         do while( r0 .lt. ro )
         
            r1 = ro    
            if( fter ) r1 = amin1( ro, tr%terx(it+1) )
            
            grad = dmdh(jl)
            a1 = ap( a0, r1-r0 )
            
            if( sign(1.,a0) .ne. sign(1.,a1) ) then
               a1 = 0.
               r1 = rp( r0, a1-a0 )
            end if
            h1 = hp( h0, a1, a0 )

c Check to see if ray has hit ground.

            if( fter ) then
               hter = hgrnd( tr%tery(it), slp(it), r1-tr%terx(it) )
               
               IF(( h1 .lt. hter ) .and. (rpe .eq. 0.)) THEN
                   
c Determine slope of line segment going into ground

                  sl2 = (h1 - h0) / (r1 - r0)
         
c Determine range at which ray hits ground and the height of the ground
c at that range.

                  slpdif = sl2 - slp(it)
       
                  r1 = (tr%tery(it)-h0+sl2*r0-slp(it)*tr%terx(it)) / 
     +                  slpdif
                  hter = hgrnd( tr%tery(it), slp(it), r1-tr%terx(it) )
                  h1 = hter
                  a1 = ap( a0, r1-r0 )
                  ihg = 1
               END IF
            end if 
            
            if(( a1 .le. 0. ) .and. ( h1 .le. htdum(jl) ))  then
               h1 = htdum(jl)
               if( h1 .gt. hter ) ihg = 0
               rad = rada1( a0, h1-h0 ) 
               a1 = -sqrt( rad )
               r1 = rp( r0, a1-a0 )    
               jl = jl - 1
               if( jl .eq. 0 ) jl = 1    
            elseif(( a1 .ge. 0. ) .and. ( h1 .ge. htdum(jl+1) )) then
               h1 = htdum(jl+1)
               rad = rada1( a0, h1-h0 ) 
               a1 = sqrt( rad )
               r1 = rp( r0, a1-a0 )
               jl = jl + 1
               if( jl .gt. rlvlep ) jl = rlvlep
            end if
       
            if( r1 .gt. ro ) then 
               r1 = ro
               a1 = ap( a0, r1-r0 )
               h1 = hp( h0, a1, a0 )
            end if   
            
            h0 = h1
            r0 = r1
            a0 = a1 
            
            if( fter ) then
               if( ihg .eq. 1 ) then 
                  ihg = 0
                  rpe = r0
                  a0 = -a0
               end if 
               if( r1 .ge. tr%terx(it+1) ) it = it + 1 
            else
               if( h0 .eq. 0. ) then
                  a0 = -a0
                  rpe = r0
               end if
            end if  

c If ray has reached ZLIM (maximum output height region) then all 
c heights for subsequent output range points will also be at ZLIM 
c - so can exit loop.
            
            if( h0 .gt. zlim ) then
               ihu = i
c               stop                         
               goto 999 ! exit
c     pg This exit is different !
            end if  
         end do 
 999     continue
         if( ihu .gt. 0 ) goto 958           !exit
         
         if( a0 .lt. 0. ) hlim(i) = 0.
         if( a0 .ge. 0. ) hlim(i) = h0
         
         ro = ro + drout
        
      end do   
958   continue
      
      if( ihu .gt. 0 ) then  
         do i = ihu, vnrout
            hlim(i) = zlim
         end do
      end if       
      
      end
      
c *************************** SUBROUTINE XYINIT *************************      

c Purpose:  Determines the initial starter field.

c Called from:  PEINIT

c Routines called: ANTPAT

      SUBROUTINE xyinit( sv )

      include 'tpem.inc'
      
      common / arrays / X(0:maxpts), Y(0:maxpts), filt(0:maxn4),   
     +                  SK(0:maxpts), CK(0:maxpts), CM(0:maxpts),
     +                  SM(0:maxpts), xlst(0:maxpts), ylst(0:maxpts)
      common / miscvar / pi2, fnorm, cnst, delp, thetamax, plcnst, nm1,  
     +                   antref, rpe, hlim(mxrout), slp(mxter), fter
      common / pevar / wl, fko, delz, n, ln, zmax, n34, con
      
      type(systemvar) sv
      
      logical fter
      
      complex refcoef, rterm, dterm, utotal
      
      sgain= sqrt( wl ) / zmax

      x(0)=0.
      y(0)=0.

      dtheta = delp / fko
      antko = fko * sv%antht
      
      DO I=1,N-1

         pk = float(i) * dtheta
         zpk = pk * antko

c Get antenna pattern factors for the direct and reflected rays.

         call antpat( sv%ipat, pk, FACD )
         call antpat( sv%ipat, -pk, FACR )

c Reflection coefficient (for now) is fixed at -1 for horizontal polarization.

         refcoef = cmplx( -1., 0. )

         rterm = cmplx( cos( zpk ), sin( zpk ) )
         dterm = conjg( rterm )

         utotal = sgain * ( facd * dterm + refcoef * facr * rterm )

         x(i) = real( utotal )
         y(i) = imag( utotal )

      end do  

      x(n) = 0
      y(n) = 0

c Filter upper 1/4 of the field.

      do i = n34, n
         attn = filt(i-n34)
         x(i) = attn*x(i)
         y(i) = attn*y(i)
      end do

      END
      
