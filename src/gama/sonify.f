      subroutine sonify(iino,kr,xlist,klist,nlist)
c
c: this subroutine checks if a pair of rays should be treated with
c: caustic correction factors for the insonified zone of the caustic.
c
      implicit integer*4(i-n)
      include 'common/caustix'
      include 'common/srlox'
      include 'common/bdcom'
      include 'common/pii'
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot)
      real xlist(7,mray,nrtot)
c
      nr=nlist(1,kr)
      nrst=nlist(2,kr)
c: no insonified zone pair if two rays have not been found:
      if(nr-nrst .lt. 2) goto 99
      nrm1=nr - 1
      jc=mod(klist(1,nr,kr)/100,10)
      jcm1=mod(klist(1,nrm1,kr)/100,10)
c: no insonified zone pair if either of two are shadow zone rays:
      if((jc .eq. 1) .or. (jcm1 .eq. 1)) goto 99
c
      tdif=.75*(xlist(3,nr,kr) - xlist(3,nrm1,kr))
      phdif=.75*(xlist(7,nr,kr) - xlist(7,nrm1,kr))
      xairy=abs(bdom(kbdf)*tdif + phdif)**(.66667)
c: don't count if they are too far from caustic (xairy=0 at caustic):
      if(xairy .gt. 2.4) goto 99
c
c: if first ray is already an insonified zone ray:
      if(jcm1 .eq. 2) then
         nrm2=nrm1 - 1
c     print *,'double caustic pair found: xairys= ',xlist(6,nrm2,kr),
c    .   xairy
c: if xairy of this pair is less than last one, only include this one:
         if(xairy .lt. xlist(6,nrm2,kr)) then
            klist(1,nrm2,kr)=klist(1,nrm2,kr) - 200
            klist(1,nrm1,kr)=klist(1,nrm1,kr) - 200
            xlist(6,nrm2,kr)=xlist(6,nrm1,kr)
         else
c: else, include last one and skip this one:
            goto 99
         endif
      endif
c
      xaip=xairy**(.25)
      call airy(-1.*xairy,ai,aip)
      fac=xlist(6,nr,kr)
      fac2=xlist(6,nrm1,kr)
      s1=sqpie*(fac2 + fac)*xaip
      s2=sqpie*(fac2 - fac)/xaip
      ccmag=abs(cmplx(s1*ai,s2*aip))
c: for doublets use xlist(6) of first ray for xairy:
      xlist(6,nrm1,kr)=xairy
      xlist(6,nr,kr)=ccmag
      klist(1,nrm1,kr)=klist(1,nrm1,kr) + 200
      klist(1,nr,kr)=klist(1,nr,kr) + 200
      return
c
99    iino=1
      return
      end 
