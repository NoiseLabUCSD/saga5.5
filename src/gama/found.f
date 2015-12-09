      subroutine found(al,ah,rl,rh,drl,drh,af,r,drf,wk1,wk,wk2,kr1,
     .   kr,kr2,iiterp,nwk,krwk,nlom,ndun,iiex,xlist,klist,nlist,
     .   range,kkk,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
c
c: this subroutine is called when rf is found to lie between
c: the ranges rl and rh, corresponding to si values al and ah.
c
      implicit integer*4(i-n)
      include 'common/dublcom'
      include 'common/srlox'
      include 'common/gamaoptions'
      include 'common/caustix'
      include 'common/pathway'
      include 'common/svp'
      include 'common/pii'
      include 'common/scatarr'
      include 'common/scatcom'
      include 'common/freqcom'
      include 'common/bdcom'
      include 'common/saga_gama_fail'
      integer*4 nlist(2,nrtot),kkk(nrtot),krwk(15)
      real range(nrtot),tmin(nrtot),dbmax(nrtot)
      real*4 angs(0:nasl),frsl(0:nfsl),dbsl(0:nfsl,0:nasl),
     .   phsl(0:nfsl,0:nasl),angb(0:nabl),frbl(0:nfbl),dbbl(0:nfbl,
     .   0:nabl),phbl(0:nfbl,0:nabl)
      logical wk1,wk,wk2,qweak,qwindow,qbdbad
      data klbad/0/
c
c      print *
c      print *,'found: al,ah,rl,rh,drl,drh = ',al,ah,rl,rh,drl,drh
c      print *,'r = ',r
c      print *
      kbad=0
      kcub=1
      a1=al
      a2=ah
      r1=rl
      r2=rh
      dr1=drl
      dr2=drh
      kloop=0
10    kloop=kloop + 1
      if(kloop .ge. 6) kcub=0 
      if(kloop .gt. 50) then
         if(abs(diff) .lt. 100.*rtol) goto 77
         print *,'kloop > 50 in found '
         print *,'af,1/af,rf,range = ',af,1./af,rf,r
         print *,'1/al,1/ah,rl,rh,drl,drh = ',1./al,1./ah,rl,rh,drl,drh
         iifail=1
         return
c         print *,'divide by zero to make dump: ',1./zzzz
c         stop
c        wk=.true.
c        nwk=nwk + 1
c        if(nwk .gt. 50) print *,'warning: nwk > 50 in found, rcheck'
c        krwk(nwk)=kr
c        goto 88
      endif
c: bisect to get next guess for af: 
      if(kcub .eq. 1) then
c      print *,'kloop = : ',kloop,1./max(.0001,a1),1./max(.0001,a2),
c     .   r1,r2,dr1,dr2
         call polfit(a1,a2,r1,r2,dr1,dr2)
         call cubroot(r,a1,a2,af,kcub)
      else
         af=(a1+a2)/2.
      endif
      call rdrcalc(af,rf,drf,ibd,rbd,drbd)
      if(drf*drl .lt. 0.) then
c        print *,'double extremum detected in found: ',1./a1,1./af,
c    .      1./a2,r1,rf,r2,dr1,drf,dr2
         kbad=1
         abad(1)=af
         abad(2)=rf
         abad(3)=drf
         do 15 kdun=1,ndun
            krr=kkk(kdun)
            nlist(1,krr)=nlist(1,krr) - 1
15       continue
         return
      endif
      diff=rf - r
c: update a1 or a2 if not close enough: 
      if(abs(diff) .gt. rtol) then
         if(diff .gt. 0.) then
            a2=af
            r2=rf
            dr2=drf 
         else
            a1=af
            r1=rf
            dr1=drf 
         endif
         goto 10
      endif
c: calculate eigenray characteristics:  
77    call timatt(af,e,t,ibd) 
c: EKW 5-3-93: rtc_calc handles attenuation correctly:
ccc   call trcalc(af,trm,trp,nps)
      call rtc_multi(af,trm,trp,nps)
      call sscalc(ntop,iisl,af,atop,amtop,ctop,ssigma,swind,nfmin,nfmin,
     .   frq,ss,angs,frsl,dbsl,phsl,nfsl,nasl)
      call sscalc(nbas,iibl,af,abas,ambas,cbas,bsigma,bwind,nfmin,nfmin,
     .   frq,bs,angb,frbl,dbbl,phbl,nfbl,nabl)
      ps=1. - (af*cs)**2
      pr=1. - (af*cr)**2
      gsl=sqrt(af*cs*cr/(sqrt(ps*pr)*abs(rf*drf)))
      raymag0=expo(fmin,e)*gsl*trm
      raymag=raymag0*ss(1,nfmin)*bs(1,nfmin)
      call record(kr,af,sign(gsl,drf),t,e,trm,trp,raymag,nps,ibd,
     .   xlist,klist,nlist,nda)
      ndun=ndun + 1 
      kkk(ndun)=kr
      raycut=dbmax(kr)*cutmag 
      tcut=tmin(kr) + twin
c: check for weakness or arrival time outside of window: 
      qbdbad=(rbd/rf .lt. -.5)
c     qbdbad=(rbd .lt. -1.*bdlam(kbdf))
      qweak=(raymag .lt. raycut)
      qwindow=(t .gt. tcut)
      if(qweak .or. qwindow .or. qbdbad) then
         if(qwindow .and. (.not. qweak)) ntcut=ntcut + 1
         if(qbdbad) nbdbad=nbdbad+1
         wk=.true.
         nwk=nwk + 1
         if(nwk .gt. 50) print *,'warning: nwk > 50 in found, rcheck'
         krwk(nwk)=kr
c: mrmin set to mmr2+1 means that a weak ray was found: 
         mrmin(nbotx1)=min0(mrmin(nbotx1),m21)
         krweak=max0(krweak,kr)
c: if ray is within 6 dB of cutoff, update krmin for safety:
         if(2.*raymag .gt. raycut) krmin=min0(krmin,kr)
      else
c: update smallest range that receives a strong ray: 
         mrmin(nbotx1)=min0(mrmin(nbotx1),kr)
         krmin=min0(krmin,kr) 
         if(ibd .eq. 1) mmbd=mmbd + 1
      endif
88    iiterp=0
      if(iiex .eq. 1) goto 99 
c: interpolate in range if delr<rterp and neither endpoint an extremum
c: or if delr<rterp/10: 
c: (interpolation of travel times near caustics does not work well.)
      delr=max(range(kr2) - r,r - range(kr1))
      if(((nlom .eq. 0) .and. (delr .lt. rterp)) .or.
     .   (delr .lt. rterp/5.)) then
         iiterp=1
         call linefit(range,kr1,kr,kr2,xlist,klist,nlist,
     .      wk1,wk,wk2,nda)
      endif
c
99    return
      end 
