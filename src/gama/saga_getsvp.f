      subroutine saga_getsvp(nline)
c
c: this subroutine reads in the svp file.
c
      implicit integer*4(i-n)
      include 'common/pii'
      include 'common/svp'
      include 'common/depth'
      include 'common/bottom' 
      include 'common/attcon' 
      include 'common/cbotcom'
      include 'common/charcom'
      include 'common/gamaoptions'
      include 'common/parms_conv'
      real zcarr(110),blug(2),zdep(0:100)
      integer*4 iblug(0:4)
      integer*4 kext(10),kadd(0:40),ksub(0:40,2)
      data iblug/0,0,2,2,2/
c
c: read the title of the svp; the parameters for material above ocean;
c: read the number of points in the ocean for which the sound velocity
c: will be given, nsvp.
      nline=0
cpln      call star(2,nline)
cpln      read(2,99,end=500,err=500) ktitle 
cpln 99    format(a)
      ktitle='                                   ,
     .                           '
      ktitle(1:12)='GAMA IN SAGA'
      call mname(ktitle,lktit)
cpln      call star(2,nline)
cpln      read(2,*,end=500,err=500) cp1(-1),cs1(-1),rho1(-1),akp1(-1),
cpln     .      aks1(-1)
      cp1(-1)=343.
      cs1(-1)=0.
      rho1(-1)=.00121
      akp1(-1)=0.0
      aks1(-1)=0.0
      if((cp1(-1) .le. 0.) .or. (cs1(-1) .lt. 0.) .or.
     .      (rho1(-1) .le. 0.) .or. (akp1(-1) .lt. 0.) .or. 
     .      (aks1(-1) .lt. 0.)) then
         print *,'illegal cp1, cs1, rho1, or k above ocean: ',cp1(-1),
     .      cs1(-1),rho1(-1),akp1(-1),aks1(-1)
         print *,'Check line number ',nline,' in .svp file...'
         stop
      endif
c
cpln      call star(2,nline)
cpln      read(2,*,end=500,err=500) nsvp
cpln      nsvp=3
      nsvp=R_ND0(NSCTOR)
      if((nsvp .lt. 2) .or. (nsvp .gt. 50)) then
         print *,'illegal # of svp points: ',nsvp 
         print *,'Check line number ',nline,' in .svp file...'
         stop
      endif
c: svp arrays start with 0, so subtract one from nsvp: 
      nsvp=nsvp-1
c: if a total of two points, make segments linear: 
      if(nsvp .eq. 1) kseg=1
c
c: read the depth (always=0.), s.v., and density at top of ocean
cpln      call star(2,nline)
cpln      read(2,*,end=500,err=500) zsvp(0),csvp(0),rho2(0)
cpln      backspace(2)
cpln      read(2,*,end=480,err=480) zsvp(0),csvp(0),rho2(0),akp2(0)
cpln      zsvp(0)=0.
cpln      csvp(0)=1508.0000
      zsvp(0)=R_Z0(1,NSCTOR)
      csvp(0)=R_C0(1,NSCTOR)
      rho2(0)=1.050
      akp2(0)=0.
      goto 482
480   akp2(0)=0.
482   if((zsvp(0) .ne. 0.) .or. (csvp(0) .le. 0.) .or. (rho2(0) .le. 0.)
     .      .or. (akp2(0) .lt. 0.)) then
         print *,'illegal svp point at ocean top: ',zsvp(0),csvp(0),
     .      rho2(0),akp2(0)
         print *,'Check line number ',nline,' in .svp file...'
         stop
      endif
      nfsvp=nsvp
      zfsvp(0)=zsvp(0)
      cfsvp(0)=csvp(0)
      cmin=csvp(0)
c
cvv   read(2,*,end=500,err=500) zsvp(1),csvp(1)
cpln      read(2,*,end=500,err=500) (zcarr(jcan),jcan=1,2*nsvp)
cpln      zcarr(1)=90.
cpln      zcarr(2)=1510.50000
cpln      zcarr(3)=129.984
cpln      zcarr(4)=1512.50000

      do 82 j=1,nsvp
         zcarr(2*j-1)=R_Z0(j+1,NSCTOR)
         zcarr(2*j)=R_C0(j+1,NSCTOR)
 82   continue

      do 83 j=1,nsvp
         zsvp(j)=zcarr(2*j-1)
         csvp(j)=zcarr(2*j)
83    continue
      if((zsvp(1) .le. 0.) .or. (csvp(1) .le. 0.)) then
         print *,'illegal svp point in ocean: ',1,zsvp(1),csvp(1)
         print *,'Check line number ',nline,' in .svp file...'
         stop
      endif
      zfsvp(1)=zsvp(1)
      cfsvp(1)=csvp(1)
      cmin=min(cmin,csvp(1))
c
      if((csvp(1) .eq. csvp(0)) .and. (kseg .ne. 1)) csvp(1)=
     .      csvp(1) + 1.e-4
      delc1=csvp(1) - csvp(0)
      delz1=zsvp(1)
      slope1=atan2(delc1,delz1)
      next=1
      kext(1)=0
c
c: read the ocean svp points.  if curved segmentation is used (kseg=3),
c: don't allow successive points to be the exactly the same.
      do 10 k=2,nsvp
cvv      read(2,*,end=500,err=500) zsvp(k),csvp(k)
         if((zsvp(k) .le. zsvp(k-1)) .or. (csvp(k) .le. 0.)) then
            print *,'illegal svp point in ocean: ',k,zsvp(k),csvp(k)
            print *,'Check line number ',nline,' in .svp file...'
            stop
         endif
         if((csvp(k) .eq. csvp(k-1)) .and. (kseg .ne. 1)) csvp(k)=
     .      csvp(k) + sign(1.e-4,csvp(k-1)-csvp(k-2))
         cmin=min(cmin,csvp(k))
         zfsvp(k)=zsvp(k)
         cfsvp(k)=csvp(k)
c: find extrema, and for the curved profile compute an averaged slope 
c: at each point.
         delc2=csvp(k)-csvp(k-1)
         delz2=zsvp(k)-zsvp(k-1)
         slope2=atan2(delc2,delz2)
         if(slope1*slope2 .lt. 0.) then 
            bsvp(1,k-1)=.75*tan(slope1)
            bsvp(2,k-1)=.75*tan(slope2)
            next=next+1
            if(next .gt. 9) then
               print *,'too many extrema in svp: ',10
               print *,'Check line number ',nline,' in .svp file...'
               stop 
            endif
            kext(next)=k-1
         else
            fac1=delz2/(delz1 + delz2)
            slopavg=tan(fac1*slope1 + (1.-fac1)*slope2)
            bsvp(1,k-1)=slopavg
            bsvp(2,k-1)=slopavg
         endif
         slope1=slope2
         delz1=delz2
10    continue
c
      next=next+1
      kext(next)=nsvp
c: decide on s.v. gradient at top and bottom: 
      if(kseg .eq. 3) then
c: find appropriate slope at ocean top and bottom:
         call bsvpfit(zsvp(1),zsvp(0),csvp(1),csvp(0),bsvp(1,1),
     .      bsvp(2,0))
         bsvp(1,0)=bsvp(2,0)
         call bsvpfit(zsvp(nsvp-1),zsvp(nsvp),csvp(nsvp-1),csvp(nsvp),
     .      bsvp(2,nsvp-1),bsvp(1,nsvp))
         bsvp(2,nsvp)=bsvp(1,nsvp)
      endif
c
      nlsub=0
      nsub=0
      nadd=0
      kadd(0)=-1
      ksub(0,1)=-1
c: kskip is the max # of svp points we try to fit with each segment:
      kskip=5
      if(svtol .eq. 0.) kskip=1
c: for each interval between extrema in the ocean svp, fit the data
c: points and try to eliminate layers if data point is within svtol of
c: the fitted profile.
      do 12 jext=2,next
         k1=kext(jext-1)
         k2=min0(kext(jext),k1+kskip)
14       call proflay(k1,k2,jok,bkold,iich)
c     print *,'called proflay: ',k1,k2,csvp(k1),csvp(k2),bsvp(2,k1),
c    .   bsvp(1,k2),jok,zsvp(k1),zsvp(k2)
         if(jok .eq. 0) then
c: if curve could not be fit, we must add a layer and then continue.
            if(k2 .eq. k1+1) then
               nadd=nadd+1
               if(nadd .gt. 40) then
                  print *,'too many layers added to svp: ',nadd
                  stop
               endif
               kadd(nadd)=k2
c: if k2 is not last point in monotonic interval, reset k1 and k2:
               if(k2 .ne. kext(jext)) then
                  k1=k2
                  k2=min0(kext(jext),k1+kskip)
                  goto 14
               endif
            else
               k2=k2-1
               goto 14
            endif
         else
c: if curve was fit, then check if the intermediate data points are
c: fit sufficiently closely.
            jfit=1
c: check if the profile fit the data points within svtol tolerance.
            do 16 kfit=k1+1,k2-1
               hfit=zsvp(kfit)-zsvp(k1) 
               if(kseg .eq. 1) then
                  cfit=csvp(k1) + g0(k2,2)*hfit
               elseif(kseg .eq. 3) then 
                  call cdcwine(csvp(k1),g0(k2,2),g1(k2,2),g2(k2,2),hfit,
     .               cfit,dcfit,dcfit2,0)
               endif
c     print *,'kfit loop: ',k1,k2,zsvp(k1),zsvp(k2),csvp(k1),
c    .  csvp(k2),kfit,hfit,cfit,csvp(kfit),abs(cfit-csvp(kfit))
               if(abs(cfit-csvp(kfit)) .gt. svtol) then
                  jfit=0
                  if(iich .eq. 1) then
                     bsvp(1,k)=bkold
                     bsvp(2,k)=bkold
                  endif
                  goto 18
               endif
16          continue
c
c: if curve did not fit the points, then skip one less layer and try
c: again. 
18          if(jfit .eq. 0) then
               k2=k2-1
               goto 14
c: if curve did fit the layer, then check if any layers have to be
c: deleted and then go on to next section of profile to fit.
            else
               if(k2 .ne. k1+1) then
                  nsub=nsub+1 
                  nlsub=nlsub + k2 - (k1+1)
                  if(nsub .gt. 40) then 
                     print *,'too many layers subtracted: ',nsub
                     stop
                  endif
                  ksub(nsub,1)=k1+1
                  ksub(nsub,2)=k2-k1-1
               endif
               if(k2 .ne. kext(jext)) then
                  k1=k2
                  k2=min0(kext(jext),k1+kskip)
                  goto 14
               endif
            endif
         endif
12    continue
c
c: add or subtract layers as necessary from the bottom up.
24    if((nsub .ne. 0) .or. (nadd .ne. 0)) then
         if(ksub(nsub,1) .gt. kadd(nadd)) then
            call sublay(ksub(nsub,1),ksub(nsub,2))
            nsub=nsub-1
         else
            call addlay(kadd(nadd))
            nadd=nadd-1
         endif
         goto 24
      endif
c
c: read the number of bottom layers, ntot=nlev
cpln      call star(2,nline)
cpln      read(2,*,end=500,err=500) ntot
      if (R_ND1(NSCTOR).eq.0) then
         ntot=0
      else
         ntot=R_ND1(NSCTOR)-1
      endif
c
      if((ntot .lt. 0) .or. (ntot .gt. 42)) then
         print *,'illegal number of layers: ',ntot
         print *,'Check line number ',nline,' in .svp file...'
         stop
      endif
c
      nlev=ntot
c: read the bottom layer acoustic parameters.
cpln      call star(2,nline)
      zdep(0)=0.
      do 50 j=1,ntot
cpln         read(2,*,end=500,err=500) kprof(j),z(j),cp1(j),cp2(j),cs1(j),
cpln     .      cs2(j),rho1(j),rho2(j),akp1(j),akp2(j),aks1(j),aks2(j),
cpln     .      (blug(jj),jj=1,iblug(kprof(j)))
         z(j)=R_Z1(j+1,NSCTOR)-R_Z1(j,NSCTOR)
         if(NTYPE.eq.1) then
            kprof(j)=1
            cp1(j)=R_C1(j,NSCTOR) ! P sound speed
            cp2(j)=R_C1(j+1,NSCTOR) !P sound speed
            rho1(j)=R_R1(1,NSCTOR)
            rho2(j)=R_R1(1,NSCTOR)
            akp1(j)=abs(R_BETA(1,NSCTOR))
            akp2(j)=abs(R_BETA(1,NSCTOR))
            blug(1)=0
         elseif(NTYPE.eq.2) then
            kprof(j)=2
            cp1(j)=R_C1(j,NSCTOR) ! P sound speed
            cp2(j)=R_BLUG1(j,NSCTOR) !P sound speed
            rho1(j)=R_R1(j,NSCTOR)
            rho2(j)=R_R1(j+1,NSCTOR)
            akp1(j)=abs(R_BETA(j,NSCTOR))
            akp2(j)=abs(R_BETA(j+1,NSCTOR))
            blug(1)=R_BLUG2(j,NSCTOR)
            blug(2)=0
         elseif(NTYPE.eq.10) then
            kprof(j)=1
            cp1(j)=R_C1(j,NSCTOR) ! P sound speed
            cp2(j)=R_C1(j+1,NSCTOR) !P sound speed
            rho1(j)=R_R1(j,NSCTOR)
            rho2(j)=R_R1(j+1,NSCTOR)
            akp1(j)=abs(R_BETA(j,NSCTOR))
            akp2(j)=abs(R_BETA(j+1,NSCTOR))
            blug(1)=0
         end if
         cs1(j)=0.
         cs2(j)=0.
         aks1(j)=0.
         aks2(j)=0.
c: 7-17-91: make negative z(j) mean absolute depth:
         if(z(j) .lt. 0.) z(j)=abs(z(j)) - zdep(j-1)
         if((z(j) .le. 0.) .or. (cp1(j) .le. 0.)
     .    .or. (cs1(j) .lt. 0.) .or. (cs2(j) .lt. 0.) .or. (rho1(j)
     .    .le. 0.) .or. (rho2(j) .le. 0.) .or. (akp1(j) .lt. 0.) .or. 
     .      (akp2(j) .lt. 0.) .or. (kprof(j) .lt. 0) .or.
     .      (kprof(j) .gt. 4) .or. (aks1(j) .lt. 0.) .or.
     .      (aks2(j) .lt. 0.)) then
            print *,'illegal profile, z, cp, cs, rho, akp, or aks: ', 
     .         j,kprof(j),z(j),cp1(j),cp2(j),cs1(j),cs2(j), 
     .         rho1(j),rho2(j),akp1(j),akp2(j),aks1(j),aks2(j)
            print *,'Check line number ',nline,' in .svp file...'
            stop
         endif
         if(kprof(j) .eq. 0) then
            bp(j)=cp2(j)
            bp2(j)=bp(j)
            cp2(j)=cp1(j) + bp(j)*z(j)
            if(cp2(j) .le. 0.) then
               print *,'illegal cp2: j,cp1,g,cp2 = ',j,cp1(j),
     .            bp(j),cp2(j)
               print *,'Check line number ',nline,' in .svp file...'
               stop 
            endif
         elseif(kprof(j) .eq. 1) then
            if(cp2(j) .le. 0.) then
               print *,'illegal cp2: ',j,cp2(j)
               print *,'Check line number ',nline,' in .svp file...'
               stop 
            endif
            bp(j)=(cp2(j)-cp1(j))/z(j)
            bp2(j)=bp(j)
         elseif(kprof(j) .eq. 2) then
            bet(j)=blug(1)
            rcon(j)=blug(2)
            bp(j)=cp2(j)
            fac=cp1(j)*(1.+bet(j))
            rad=fac**2 + 2.*bp(j)*fac*z(j)
            if(rad .lt. 0.) then
               print *,'illegal blug profile: ',j,cp1(j),bp(j),bet(j) 
               print *,'Check line number ',nline,' in .svp file...'
               stop 
            endif
            if(bp(j) .eq. 0.) then
               cp2(j)=cp1(j)
               bp2(j)=bp(j)
               kprof(j)=1
            else
               cp2(j)=sign(1.,fac)*sqrt(rad) - bet(j)*cp1(j)
               bp2(j)=bp(j)*(1.+bet(j))/(cp2(j)/cp1(j) + bet(j))
            endif
         elseif(kprof(j) .eq. 3) then
c: kprof=3 means that cp2 is given instead of gradient g for blug: 
            kprof(j)=2
            bet(j)=blug(1)
            rcon(j)=blug(2)
            delc=cp2(j) - cp1(j)
c: beta=-999. means to set beta such that gradient g is continuous: 
            if(abs(bet(j) + 999.) .lt. 1.e-5) then
               if(j .eq. 1) then
            print *,'kprof=3 with beta=-999 cannot be for layer 1'
            print *,'Check line number ',nline,' in .svp file...'
                  stop
               endif
               if(abs(cp2(j-1)-cp1(j)) .gt. .05) then
                  print *,'blug layer with cont g specified, 
     .                    but c not cont ',j
                  print *,'Check line number ',nline,' in .svp file...'
                  stop
               endif
               bp(j)=bp2(j-1) 
               bet(j)=delc**2/(2.*cp1(j)*(bp(j)*z(j) - delc)) - 1.
               bp2(j)=bp(j)*(1.+bet(j))/(cp2(j)/cp1(j) + bet(j))
               write(8,410) j,cp1(j),cp2(j),bp(j),bp2(j),bet(j)
410   format(/'BLUG LAYER WITH CONTINUOUS GRADIENT FOR LAYER ',i2,':'/
     .   'CP1,CP2,G1,G2,BETA = ',2(f8.3,2x),3(f7.3,2x))
            else
               bp(j)=delc*(1. + delc/(2.*cp1(j)*(1.+bet(j))))/z(j)
               bp2(j)=bp(j)*(1.+bet(j))/(cp2(j)/cp1(j) + bet(j))
            endif
         elseif(kprof(j) .eq. 4) then
c: kprof=4 means that cp2 is given 2nd, g instead of beta at end:
            kprof(j)=2
            bp(j)=blug(1)
            rcon(j)=blug(2)
            delc=cp2(j) - cp1(j)
            bden=bp(j)*z(j) - delc
            if(bden .eq. 0.) then
               print *,'blug layer ',j,'illegal: cp1,cp2,g gives ',
     .            'a linear profile.'
               stop
            endif
            bet(j)=delc**2/(2.*cp1(j)*bden) - 1.
            bp2(j)=bp(j)*(1.+bet(j))/(cp2(j)/cp1(j) + bet(j))
         endif
c
         dakp(j)=(akp2(j) - akp1(j))/z(j)
         zdep(j)=zdep(j-1) + z(j)
50    continue
c
c: read the sound speeds and density at the top of the substrate.
cpln      call star(2,nline)
cpln      read(2,*,end=500,err=500) cp1(ntot+1),cs1(ntot+1),rho1(ntot+1), 
cpln     .      akp1(ntot+1),aks1(ntot+1)
cpln      cp1(ntot+1)=1652.362
      cp1(ntot+1)=R_C2(NSCTOR)
      cs1(ntot+1)=0.
cpln      rho1(ntot+1)=1.8
      rho1(ntot+1)=R_R2(NSCTOR)
cpln      akp1(ntot+1)=.1437
      akp1(ntot+1)=abs(R_BETA(2,NSCTOR))
      aks1(ntot+1)=0.
      if((cp1(ntot+1) .le. 0.) .or. (cs1(ntot+1) .lt. 0.) .or.
     .   (rho1(ntot+1) .le. 0.) .or. (akp1(ntot+1) .lt. 0.) .or.
     .   (aks1(ntot+1) .lt. 0.)) then
         print *,'illegal cp1, cs1, rho1, or k in substrate: ',
     .      cp1(ntot+1),cs1(ntot+1),rho1(ntot+1),akp1(ntot+1),
     .      aks1(ntot+1)
         print *,'Check line number ',nline,' in .svp file...'
         stop
      endif
      bp(ntot+1)=0. 
c
c: the following quantities are used for convenience.
c: ocean surface: 
      cp2(-2)=csvp(0)
      cs2(-2)=0.
      rho2(-2)=rho2(0)
c: ocean bottom: 
      cp2(0)=csvp(nsvp)
      cs2(0)=0.
c
      aks2(-2)=0.
      aks2(0)=0.
c
c: set mismatch to 0 if cp and rho are continuous across interface.
c: set mismatch to -1 if gradients are also continuous.
      do 62 j=0,ntot
         if((cp2(j) .eq. cp1(j+1)) .and. (rho2(j) .eq. rho1(j+1)))
     .         then 
c: NEW 10-10-90: DON'T ALLOW MISMCH=-1 AT WATER-SEDIMENT INTERFACE:
            if(j .ne. 0 .and. abs(bp2(j) - bp(j+1)) .lt. .05) then
               mismch(j)=-1
            else
               mismch(j)=0
            endif
         else
            mismch(j)=1
         endif
62    continue
c
c: check to see if attenuations and velocities fit the criterion that 
c: the imaginary part of the sound velocity be much smaller than the
c: real part when accounting for attenuation in the reflection coeff. 
      cp2(ntot+1)=0.
      cs2(ntot+1)=0.
      akp2(ntot+1)=0.
      aks2(ntot+1)=0.
c
      prl=2728.
      do 95 j=1,ntot+1
         if((akp1(j)*cp1(j) .gt. prl) .or. (akp2(j)*cp2(j) .gt. prl)) 
     .         then 
         write(8,405) j
405   format(/'WARNING: COMPRESSIONAL WAVE ATTENUATION IN ',
     .   'BOTTOM LAYER ',i2/'IS TOO LARGE FOR COMPLEX SOUND ',
     .   'SPEED APPROXIMATION.')
c           akp1(j)=0.
c           akp2(j)=0.
         endif
         if((aks1(j)*cs1(j) .gt. prl) .or. (aks2(j)*cs2(j) .gt. prl)) 
     .         then 
         write(8,406) j
406   format(/'WARNING: SHEAR WAVE ATTENUATION IN ',
     .   'BOTTOM LAYER ',i2/'IS TOO LARGE FOR COMPLEX SOUND ',
     .   'SPEED APPROXIMATION.')
c           aks1(j)=0.
c           aks2(j)=0.
         endif
95    continue
c: EKW 5-3-93: compute layer parameters needed for rtc_calc, which
c: handles attenuation correctly:
c: ocean surface/"air" interface:
      call layparm(csvp(0),cp1(-1),0.,cs1(-1),rho2(0),rho1(-1),
     .   akp2(0),akp1(-1),0.,aks1(-1),fax1,nx_p2(1),nx_s1(1),
     .   nx_s2(1),mratx(1),cp1x(1),cp2x(1),cs1x(1),cs2x(1))
      call layparm(cp1(-1),csvp(0),cs1(-1),0.,rho1(-1),rho2(0),
     .   akp1(-1),akp2(0),aks1(-1),0.,fax1,ny_p1(1),ny_s2(1),
     .   ny_s1(1),mraty(1),dum,dum,dum,dum)
c: ocean bottom/sediment interface:
      call layparm(csvp(nsvp),cp1(1),0.,cs1(1),rho2(0),rho1(1),
     .   akp2(0),akp1(1),0.,aks1(1),fax1,nx_p2(2),nx_s1(2),
     .   nx_s2(2),mratx(2),cp1x(2),cp2x(2),cs1x(2),cs2x(2))
      call layparm(cp1(1),csvp(nsvp),cs1(1),0.,rho1(1),rho2(0),
     .   akp1(1),akp2(0),aks1(1),0.,fax1,ny_p1(2),ny_s2(2),
     .   ny_s1(2),mraty(2),dum,dum,dum,dum)
c: bottom layer interfaces:
      do 395 j=1,ntot
         k=j+2
         call layparm(cp2(j),cp1(j+1),cs2(j),cs1(j+1),rho2(j),
     .      rho1(j+1),akp2(j),akp1(j+1),aks2(j),aks1(j+1),fax1,
     .      nx_p2(k),nx_s1(k),nx_s2(k),mratx(k),
     .      cp1x(k),cp2x(k),cs1x(k),cs2x(k))
         call layparm(cp1(j+1),cp2(j),cs1(j+1),cs2(j),rho1(j+1),
     .      rho2(j),akp1(j+1),akp2(j),aks1(j+1),aks2(j),fax1,
     .      ny_p1(k),ny_s2(k),ny_s1(k),mraty(k),dum,dum,dum,dum)
         aturn(j+2)=1./max(cp1(j),cp2(j))
395   continue
      nd=ntot+2
cxx   print *,nx_p2(1:nd),ny_p1(1:nd),nx_s1(1:nd),ny_s1(1:nd),
cxx  .   nx_s2(1:nd),ny_s2(1:nd),mratx(1:nd),mraty(1:nd),
cxx  .   cp1x(1:nd),cp2x(1:nd),cs1x(1:nd),cs2x(1:nd)
c
c: for bottom loss calculations, make svp isospeed:
      if(iiblt .eq. 1) then
c        zsvp(1)=max(10.,100.*z(ntot))
c: keep user-chosen water depth, zs, and zr:
         zsvp(1)=zsvp(nsvp)
         csvp(1)=csvp(nsvp)
         csvp(0)=csvp(nsvp)
         nsvp=1
         g0(1,1)=0.
         g0(1,2)=0.
         kseg=1
      endif
      zlev(4)=zsvp(nsvp)
      do 60 j=1,ntot
         zlev(j+4)=zlev(j+3) + z(j)
60    continue
c
c: write out the svp to a scratch file for use in newdep.
      open(unit=56,status='scratch',form='unformatted')
      rewind(56)
      write(56) nsvp,zsvp,csvp,bsvp,g0,g1,g2

      return
      end 
