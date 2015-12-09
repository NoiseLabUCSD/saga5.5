      subroutine shadcal(iiwk,xfac,s1,s2,rc,r,dbmax,iibdc,ac,gsl,t,e,
     .   trm,trp,nturn,s1f,s2f,dr2,kr,tcorr,tmin,ibd,alo,
     .   xlist,klist,nlist,nda,psmag,csmag,psph,csph,xlim)
c
      implicit integer*4(i-n)
      include 'common/caustix'
      include 'common/gamaoptions'
      include 'common/srlox'
      include 'common/pathway'
      include 'common/freqcom'
      include 'common/pii'
      include 'common/discint'
      include 'common/bdcom'
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot)
      real xlist(7,mray,nrtot)
      complex plpar,cc
      logical qweak,qwindow
      common /info/ienv_info,lu_env,iiwrite
      integer ienv_info,lu_env,iiwrite
c
      delr=r - rc
      xairy=xfac*delr
c     tc=t + delr*(ac - tcorr*delr)
      tc=t + delr*ac
      call airy(-1.*xairy,ai,aip)
      cc=cmplx(s1*ai,s2*aip)
      ccmag=abs(cc)
      raycut=dbmax*cutmag
      tcut=tmin + twin
      qweak=(ccmag .le. raycut)
      qwindow=(tc .gt. tcut)
      if(qweak .or. qwindow) then
         iiwk=1
         if(qwindow .and. (.not. qweak)) ntcut=ntcut+1
         goto 99
      endif
c
      nr=nlist(1,kr)
      nrst=nlist(2,kr)
      s1fx=s1f
      s2fx=s2f
      trpx=trp
c: check if a doublet pair is being treated as a shadow ray: 
      if(xairy .gt. 0.) then
c: delete second doublet ray if it was found and if not a shadow ray: 
         jc=mod(klist(1,nr,kr)/100,10)
cxx   print *,'2nd dub ch: jc,kr = ',jc,kr,'; nr,nrst = ',nr,nrst
         if((nr .gt. nrst) .and. (jc .eq. 0)) then 
            nr=nr-1
            nlist(1,kr)=nr
            jc=mod(klist(1,nr,kr)/100,10)
c: delete first doublet ray if it was found and if not a shadow ray: 
cxx   print *,'1st dub ch: jc,kr = ',jc,kr,'; nr,nrst = ',nr,nrst
            if((nr .gt. nrst) .and. (jc .eq. 0)) nlist(1,kr)=nr-1
         endif
      elseif(iibdc .eq. 1) then
         if(xairy .lt. xlim) then
            iiwk=1
            goto 99
         endif
         fac1=xairy/xlim
         fac2=1. - fac1
         fac=psmag**fac1*csmag**fac2
         phfac=fac1*psph + fac2*csph
c: check if we should correct partially transmitted ray:
cxx      print *,'calling bdparsh: kr,nr = ',kr,nr
         call bdparsh(nr,kr,nrp,xlist,klist,nlist,nda)
         if(nrp .ne. 0) then
            fdep=exp(max(-656.,fax2*bdf(kbdf)*xlist(4,nrp,kr)))
            trmold=xlist(5,nrp,kr)
            xlist(5,nrp,kr)=fac/(fdep*abs(xlist(2,nrp,kr)))
            xlist(7,nrp,kr)=phfac
            if(iiwrite.gt.0)
     .       print *,'partial ray modified: ',trmold,xlist(5,nrp,kr)
            goto 99
         else
c: correct shadow ray reflection/transmission coefficient:
            s1fx=fac*s1f/ccmag
            s2fx=fac*s2f/ccmag
cxx   print *,'cc,s1,s2,ai,aip = ',cc,s1,s2,ai,aip
            trpx=phfac - piequ - atan2(imag(cc),real(cc))
         endif
      endif
      nps=nturn+100
      if(ibd .eq. 1) mmbd=mmbd + 1
      call record(kr,ac,gsl,tc,e,trm,trpx,ccmag,nps,ibd,
     .   xlist,klist,nlist,nda)
      xfackr=xairy/bdom23(kbdf)
      call record(kr,xfackr,dr2,s1fx,s2fx,rc,rc,ccmag,nps,ibd,
     .   xlist,klist,nlist,nda)
      tmin=min(tmin,tc)
      dbmax=max(dbmax,ccmag)
c: update smallest range that receives a strong ray: 
      mrmin(nbotx1)=min0(mrmin(nbotx1),kr)
      krmin=min0(krmin,kr) 
c
99    return
      end 
