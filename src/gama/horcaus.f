      subroutine horcaus(mbp,xlist,klist,nlist)
c
c: this subroutine checks for caustic rays that have different ocean
c: paths due to the fact that one travels up from/to src/rec and one
c: travels down from/to src/rec.  mbp1 is the path counter for the
c: first of the rays, the current mbp is that for the second.
c
      implicit integer*4(i-n)
      include 'common/srlox'
      include 'common/hccom'
      include 'common/caustix'
      include 'common/bdcom'
      include 'common/pii'
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot),ktmp(2)
      real xlist(7,mray,nrtot),xtmp(7)
c
c: loop through the ranges:
      if(krish(kp0) .eq. 0) return
      mbp1=mbphc(kp0)
      do 10 kr=krdb1(kp0),krdb2(kp0),krish(kp0)
         nr=nlist(1,kr)
c: skip if no rays found in previous disc. interval:
      print *,'horcaus: kr,nrst,nrnow = ',kr,nlist(2,kr),nr
         if(nlist(2,kr) .ge. nr) goto 10
         jcaus=mod(klist(1,nr,kr)/100,10)
      print *,'kr,jcaus = ',kr,jcaus
c: skip if ray was a shadow or insonified zone ray already:
         if(jcaus .ne. 0) goto 10
c: go through previous rays found at current range:
      print *,'checking others: kr,nr = ',kr,nr
         do 20 nr2=nr-1,1,-1
            jcaus=mod(klist(1,nr2,kr)/100,10)
c     print *,'kr,nr2,jcaus = ',kr,nr2,jcaus
c: skip if previous ray was a shadow or insonified zone ray:
            if(jcaus .ne. 0) goto 20
            mbp2=klist(2,nr2,kr)/1000
      print *,'nr2,mbp1,mbp2 = ',nr2,mbp1,mbp2 
c: if ray has correct mbp, check if the two rays are a caustic pair:
            if(mbp2 .eq. mbp1) then
c: add a -pi/2 phase shift for extra turning of 1st ray:
               xlist(7,nr2,kr)=xlist(7,nr2,kr) - pieh
               tdif=.75*(xlist(3,nr,kr) - xlist(3,nr2,kr))
               phdif=.75*(xlist(7,nr,kr) - xlist(7,nr2,kr))
      print *,'corrected phases: ',xlist(7,nr,kr)*piedeg,
     .   xlist(7,nr2,kr)*piedeg,phdif*piedeg
               xairy=abs(bdom(kbdf)*tdif + phdif)**(.66667)
      print *,'xairy = ',xairy
c: don't count if they are too far from caustic (xairy=0 at caustic):
               if(xairy .gt. 2.4) goto 99
               xaip=xairy**(.25)
               call airy(-1.*xairy,ai,aip)
               fac=xlist(6,nr,kr)
               fac2=xlist(6,nr2,kr)
               s1=sqpie*(fac2 + fac)*xaip
               s2=sqpie*(fac2 - fac)/xaip
               ccmag=abs(cmplx(s1*ai,s2*aip))
c: for doublets use xlist(6) of first ray for xairy:
               xlist(6,nr2,kr)=xairy
               xlist(6,nr,kr)=ccmag
               klist(1,nr2,kr)=klist(1,nr2,kr) + 200
               klist(1,nr,kr)=klist(1,nr,kr) + 200
c: slip current ray in just after its caustic partner and shift others:
               ktmp(1)=klist(1,nr,kr)
               ktmp(2)=klist(2,nr,kr)
               do 940 jcan=1,7
                  xtmp(jcan)=xlist(jcan,nr,kr)
940            continue
               do 30 nrr=nr,nr2+2,-1
                  nrrm1=nrr-1
                  klist(1,nrr,kr)=klist(1,nrrm1,kr)
                  klist(2,nrr,kr)=klist(2,nrrm1,kr)
                  do 942 jcan=1,7
                     xlist(jcan,nrr,kr)=xlist(jcan,nrrm1,kr)
942               continue
30             continue
               klist(1,nr2+1,kr)=ktmp(1)
               klist(2,nr2+1,kr)=ktmp(2)
               do 944 jcan=1,7
                  xlist(jcan,nr2+1,kr)=xtmp(jcan)
944            continue
      print *,'found dub and shifted: kr,nr = ',kr,nr
               goto 10
c: skip to next range if we have passed the mbp1 we are looking for:
            elseif(mbp2 .lt. mbp1) then
               goto 10
            endif
20       continue
10    continue
c
99    continue
      iifc=0
c: check if false caustic found, and if so, which mbp to delete:
      if(iifcs(kp0,1) .eq. 1) then
         mbpdel=mbp
         iifc=1
      elseif(iifcs(kp0,2) .eq. 1) then
         mbpdel=mbp1
         iifc=1
      endif
      print *,'horcaus: iifc,iifcs(kp0,1:2),mbpdel = ',
     .      iifc,(iifcs(kp0,jcan),jcan=1,2),mbpdel
      if(iifc .eq. 1) then
      print *,'krsh1,krsh2,krish = ',krsh1(kp0),krsh2(kp0),-krish(kp0),
     .   mbpdel
c: go through previous rays found at current range:
         do 40 kr=krsh1(kp0),krsh2(kp0),-krish(kp0)
            nr=nlist(1,kr)
            do 50 nr2=nr,1,-1
               mbp2=klist(2,nr2,kr)/1000
               if(mbp2 .eq. mbpdel) then
                  jcaus=mod(klist(1,nr2,kr)/100,10)
                  if(jcaus .ne. 0) goto 40
c: delete ray by shifting all subsequent rays downward in list:
      print *,'deleting loner ray from false caustic: kr,nr2 = ',kr,nr2
                  do 55 nrr=nr2,nr-1
                     klist(1,nrr,kr)=klist(1,nrr+1,kr)
                     klist(2,nrr,kr)=klist(2,nrr+1,kr)
                     do 948 jcan=1,7
                        xlist(jcan,nrr,kr)=xlist(jcan,nrr+1,kr)
948                  continue
55                continue
                  nlist(1,kr)=nr - 1
      print *,'nlist upd: kr,nr = ',kr,nr-1
                  goto 40
               elseif(mbp2 .lt. mbpdel) then
                  goto 40
               endif
50          continue
40       continue
      endif
      krish(kp0)=0
c
      return
      end
