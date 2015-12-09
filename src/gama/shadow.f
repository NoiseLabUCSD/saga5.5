      subroutine shadow(al,ah,kf,a1,a2,iifc,
     .   xlist,klist,nlist,range,wk,ap,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
c
c: this subroutine is called when a caustic point is found by
c: loma.  it checks whether any of the receiver ranges lie in
c: the shadow region around this point.  if so, an eigenray with
c: jcaus=1 is entered into the list.  for ranges just barely in
c: the insonified region (0<xairy<.05), the shadow zone formulas are
c: also used.
c
      implicit integer*4(i-n)
      include 'common/srlox'
      include 'common/caustix'
      include 'common/pathway'
      include 'common/svp'
      include 'common/gamaoptions'
      include 'common/freqcom'
      include 'common/pii'
      include 'common/scatarr'
      include 'common/scatcom'
      include 'common/causbad'
      include 'common/inivar'
      include 'common/discint'
      include 'common/bottom'
      include 'common/hccom'
      include 'common/bdcom'
      real al(3),ah(3),a1(3),a2(3)
      real range(nrtot),tmin(nrtot),dbmax(nrtot)
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot)
      complex csf
      logical wk(nrtot),ap(nrtot)
      integer*4 kr2next(15)
      data bfac/1.5/,shadext/100/
c
      if(iicaus .eq. 0) return
      iifc=0
c: if exact caustic was not found (kf=0) skip shadow check:
      if(kf .ne. 1) then
         if(al(3) .lt. 0.) then
            krd1=mmr2
            krd2=mmr0
            krinc=-1
         else
            krd1=mmr0
            krd2=mmr2
            krinc=1
         endif
         goto 22
      endif
c
      ac=(al(1)+ah(1))/2. 
c: incorp of r''' seems to hurt as often as it helps, so don't do
c     call dr3calc(ac,rc,drc,dr2c,dr3c,ibd,rbd,drbd,dr2bd)
      call dr2calc(ac,rc,drc,dr2c,ibd,rbd,drbd,dr2bd)
      sg=sign(1.,dr2c)
      dr2p=abs(dr2c/2.)**(-.33333)
      xfac=sg*dr2p*bdom23(kbdf)
c
c: calculate quantities that don't depend on range: 
      call timatt(ac,e,t,ibd) 
c: EKW 5-3-93: rtc_calc handles attenuation correctly:
ccc   call trcalc(ac,trm,trp,nturn)
      call rtc_multi(ac,trm,trp,nturn)
c     print *,'shadow: ',1./ac,rc,drc,dr2c,nturn
c: corrections due to r''' do not seem to help, so don't include: 
c     tcorr=dr3c/(6.*dr2c**2) 
      fdep=expo(bdf(kbdf),e)
      ps=1. - (ac*cs)**2
      pr=1. - (ac*cr)**2
c     gg2=(ps + pr - ps*pr)/(ps*pr*ac) - dr3c/(3.*dr2c)
c     g2corr=-1.*dr3c/(3.*dr2c)
      gg2=(ps + pr - ps*pr)/(ps*pr*ac)
      gsl=sqrt(ac*cs*cr/(sqrt(ps*pr)*abs(rc)))
      call sscalc(ntop,iisl,ac,atop,amtop,ctop,ssigma,swind,nfmin,nfmin,
     .   frq,ss,angs,frsl,dbsl,phsl,nfsl,nasl)
      call sscalc(nbas,iibl,ac,abas,ambas,cbas,bsigma,bwind,nfmin,nfmin,
     .   frq,bs,angb,frbl,dbbl,phbl,nfbl,nabl)
      gsltrm=gsl*trm
c
      sfac9=fdep*ss(1,nfmin)*bs(1,nfmin)
      if(sfac9 .eq. 0.) return
c: call smcaus to find s1 and s2 slightly on insonified side where
c: xairy should be about .10:
      rsm=rc + .10/xfac
      call smcaus(a1,al,ah,a2,rsm,s1,s2,
     .   xairy,nogo,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      if(nogo .eq. 1) then
         s1f=sq2pie*dr2p*gsltrm
         s2f=.5*sg*sq2pie*gg2*dr2p**2*gsltrm
         s1=s1f*sfac9*bdom16(kbdf)
         s2=s2f*sfac9/bdom16(kbdf)
c     print *,'shadow from dr2: 1/a,rc = ',1./ac,rc
      else
         xfac=xairy/(rsm - rc)
         s1f=s1/(sfac9*bdom16(kbdf))
         s2f=s2*bdom16(kbdf)/sfac9
c     print *,'shadow from smcaus: 1/a,rc = ',1./ac,rc
      endif
c
      delta=-3.33/xfac
      dellam=delta/bdlam(kbdf)
c: check if 40-dB shadow zone is > 500 lambda, in which case the caustic
c: is an irregular kink in the r-a plot and shadow zone is ignored:
      if(abs(dellam) .gt. shadext) then
         ncbad=ncbad + 1
         if(ncbad .le. 20) then
            ctbad(ncbad)=1./ac
            rcbad(ncbad)=rc
            delbad(ncbad)=dellam
         endif
         iifc=1
         r07=rc + sign(5.*bdlam(kbdf),xfac)
c     print *,'unreal caustic: 1./ac,rc,icpth = ',1./ac,rc,r07,icpth
      else
c: range at which xairy=.07 is dividing line for doublet/shadow zones:
         r07=rc + .07/xfac
      endif
c: krd1,krd2 point to the first,last ranges to check for dublets:
      if(xfac .gt. 0.) then
c: caustic is local minimum:
         krd2=mmr2
         krs2=mmr0
         krinc=1
         krd1=mmr2+1
         do 10 kr=mmr0,mmr2
            if(range(kr) .gt. r07) then
               krd1=kr
               goto 11
            endif
10       continue
11       krs1=krd1-1
      else
c: caustic is local maximum:
         krd2=mmr0
         krs2=mmr2
         krinc=-1
         krd1=mmr0-1
         do 12 kr=mmr2,mmr0,-1
            if(range(kr) .lt. r07) then
               krd1=kr
               goto 13
            endif
12       continue
13       krs1=krd1+1
      endif
c     print *,'shadow: r07,krd1,kdr2,krs1,krs2 = ',r07,krd1,krd2,
c    .   krs1,krs2
c
c: if this is a false caustic, delete the dublet rays within 5
c: wavelengths of the caustic:
      if(iifc .eq. 1) then
         krs22=krs2
         if(iihc .eq. 1) krs2=krs1 + krinc
         do 18 kr=krs1,krs22,-krinc
            xairy=xfac*(range(kr) - rc)
            if(xairy .gt. 0.) then
               nr=nlist(1,kr)
               nrst=nlist(2,kr)
               if(nr .le. nrst) goto 18
               jc=mod(klist(1,nr,kr)/100,10)
               if(jc .eq. 0) then
                  nr=nr-1
                  nlist(1,kr)=nr
                  jc=mod(klist(1,nr,kr)/100,10)
                  if((nr .gt. nrst) .and. (jc .eq. 0)) then
                     nlist(1,kr)=nr-1
                  elseif(iihc .eq. 1) then
c: for horizontal caustics, set krs2 for use in horcaus:
                     krs2=kr
                  endif
               endif
            else
               goto 22
            endif
18       continue
         goto 22
      endif
c
c: check if this is a beam displacement caustic:
      iibdc=0
      if(ibd .eq. 1) then
         a1cs=a1(1)*cs
         if(asin(ac*cs) - asin(a1cs) .lt. .035) then
            iibdc=1
            cturn=1./a1(1)
            ceps=cturn/500
            cdif1=min(abs(cturn-cp1(-1)),abs(cturn-cs1(-1)))
            cdif2=min(abs(cturn-cp1(ndp+1)),abs(cturn-cs1(ndp+1)))
            if(cdif2 .lt. cdif1) then
               nbnc=nttot(ndp)/2
            else
               nbnc=ntop
            endif
            thcr=acos(a1cs)
            ztot=cs*t*sqrt(1. - a1cs**2)
            sindth=bfac/(2.*(bdom(kbdf)/cs)*ztot**1.4*
     .         float(nbnc)**(-1.58))
            if(abs(sindth) .gt. 1.) then
               dth=thcr/20.
            else
               dth=2.*asin(sindth)
            endif
            aps=cos(thcr - dth)/cs
            apmin=.97*adlop(1) + .03*adhip(1)
            apmax=.03*adlop(1) + .97*adhip(1)
            aps=max(apmin,min(apmax,aps))
            call rdrcalc(aps,rps,drps,ibdps,rbdps,drbdps)
            xlim=(rps-rc)*xfac
c     print *,'xlim = ',xlim
            call timatt(aps,eps,tps,ibdps)
c: EKW 5-3-93: rtc_calc handles attenuation correctly:
ccc         call trcalc(aps,trmps,trpps,npsps)
            call rtc_multi(aps,trmps,trpps,npsps)
c     print *,'iibdc=1: thcr,1./aps,rps,cdif1,cdif2 = ',thcr*piedeg,
c    .   1./aps,rps,cdif1,cdif2
            ps=1. - (aps*cs)**2
            pr=1. - (aps*cr)**2
            gps=sqrt(aps*cs*cr/(sqrt(ps*pr)*abs(rps*drps)))
            psmag=expo(bdf(kbdf),eps)*gps*trmps
            psph=trpps
            csf=cmplx(.3550*s1,-.2588*s2)
            csmag=abs(csf)
            csph=trp + piequ + atan2(aimag(csf),real(csf))
         endif
      endif
c
      iiwk=0
      do 20 kr=krs1,krs2,-krinc
         call shadcal(iiwk,xfac,s1,s2,rc,range(kr),dbmax(kr),iibdc,
     .      ac,gsl,t,e,trm,trp,nturn,s1f,s2f,dr2c,kr,tcorr,tmin(kr),
     .      ibd,a1(1),xlist,klist,nlist,nda,psmag,csmag,psph,csph,
     .      xlim)
         if(iiwk .eq. 1) goto 22
20    continue
c     print *,'shadows done: kr = ',kr
c: if caustic is local maximum and all ranges>rc, don't allow krmin
c: to indicate that all rays are weak.  Those with more turnings
c: might make caustic move out in range and be stronger there.
      if((xfac .lt. 0.) .and. (krs1 .eq. mmr0)) then
         krmin=min0(krmin,mmr0)
      endif
22    continue
c
c: at ranges in insonified zone, check if insonified zone pair found:
      iino=0
c     print *,'insonifed checks: ',krd1,krd2,krinc
      do 30 kr=krd1,krd2,krinc
         call sonify(iino,kr,xlist,klist,nlist)
         if(iino .eq. 1) goto 33
30    continue
33    continue
c     print *,'insonifieds done: kr = ',kr
c: save info for horcaus if this could be a horizontal caustic:
      if(iihc .eq. 1) then
         krsh1(kp0)=krs1
         krsh2(kp0)=krs2
         krish(kp0)=krinc
         krdb1(kp0)=krd1
         krdb2(kp0)=krd2
      endif
c
      return
      end 
