      subroutine rcheck(a1,a2,nlom,xlist,klist,nlist,
     .   range,wk,ap,kkk,xxx,yyy,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
c
c: this subroutine is called to check if a receiver range is
c: bracketted by the monotonic range interval from a1 to a2.
c: if so, the subroutine "found" is called.
c
      implicit integer*4(i-n)
      include 'common/dublcom'
      include 'common/srlox'
      include 'common/gamaoptions'
      include 'common/saga_gama_fail'
      integer*4 nlist(2,nrtot),kr2next(50),krwk(50)
      real a1(3),a2(3)
      real range(nrtot),xxx(nrnf),yyy(nrnf)
      logical wk(nrtot),ap(nrtot),wkd
      data wkd/.false./
c
c: ndun keeps track of which kr's have been done, in case they have
c: to be undone: 
      if(iidiag .ne. 0) print *,'rcheck: 1/a1,a/a2,r1,r2,dr1,dr2,nlom=',
     .   1./max(1.e-8,a1(1)),1./a2(1),a1(2),a2(2),a1(3),a2(3),nlom
      ndun=0
c: set al, rl, drl to values corr to lower range: 
      if(a1(2) .lt. a2(2)) then
         al=a1(1)
         rl=a1(2)
         drl=a1(3)
         ah=a2(1)
         rh=a2(2)
         drh=a2(3)
      else
         ah=a1(1)
         rh=a1(2)
         drh=a1(3)
         al=a2(1)
         rl=a2(2)
         drl=a2(3)
      endif
c
c: go through the range intervals over which r is increasing: 
      if((range(mmr0) .gt. rh) .or. (range(mmr2) .lt. rl)) goto 199
      nwk=0
c: find interval in current range set between rl and rh: 
c: [start at mmr0, the smallest range for which strong eigenrays were 
c: found for the previous number of bottom bounces.]
      do 12 kr=mmr0,mmr2
         if(range(kr) .ge. rl) then
            nr1=kr
            goto 13 
         endif
12    continue
13    continue
      do 14 kr=mmr2,mmr0,-1
         if(range(kr) .le. rh) then
            nr2=kr
            goto 15 
         endif
14    continue
15    continue
c
      if(nr1 .gt. nr2) goto 199
      do 16 kr=nr1,nr2
         wk(kr)=.false.
16    continue
c
c: find rays at end points
c      print *,'first found call: nr1 = ',nr1
cpln      write(6,*)'rcheck: found #1'
      call found(al,ah,rl,rh,drl,drh,xxx(nr1),range(nr1),yyy(nr1),
     .   wkd,wk(nr1),wkd,nr1,nr1,nr1,iitp,nwk,krwk,nlom,ndun,1,
     .   xlist,klist,nlist,range,kkk,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      if(iifail.eq.1) return
      if(kbad .eq. 1) return
      if(nr1 .eq. nr2) goto 99
c      print *,'second found call: nr2 = ',nr2
cpln      write(6,*)'rcheck: found #2',nr1,nr2
cpln      write(6,*)xxx(nr1),range(nr1)
cpln      write(6,*)xxx(nr2),range(nr2)
cpln      write(6,*)yyy(nr1),yyy(nr2)
cpln      write(6,*)wk(nr1),wk(nr2)
cpln      write(6,*)
      call found(xxx(nr1),ah,range(nr1),rh,yyy(nr1),drh,xxx(nr2),
     .   range(nr2),yyy(nr2),wkd,wk(nr2),wkd,nr2,nr2,nr2,iitp,nwk,
     .   krwk,nlom,ndun,1,xlist,klist,nlist,range,kkk,nda,tmin,dbmax,
     .   angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      if(iifail.eq.1)return
      if(kbad .eq. 1) return
      if((nr1+1 .eq. nr2) .or. (wk(nr1) .and. wk(nr2))) goto 99
c
c: bisect range interval and find eigenray using nearest found values: 
      inc=0
      kr1=nr1
      kr2=nr2
19    continue
      idif=kr2-kr1
      kr=kr1 + idif/2
c      print *,'call found for kr1,kr,kr2 = ',kr1,kr,kr2
c      print *,'ranges(kr1,kr,kr2) = ',range(kr1),range(kr),range(kr2)
c      write(6,*)'rcheck: found #3'
      call found(xxx(kr1),xxx(kr2),range(kr1),range(kr2),yyy(kr1),
     .   yyy(kr2),xxx(kr),range(kr),yyy(kr),wk(kr1),wk(kr),wk(kr2),
     .   kr1,kr,kr2,iitp,nwk,krwk,nlom,ndun,0,xlist,klist,nlist,range,
     .   kkk,nda,tmin,dbmax,angs,frsl,dbsl,phsl,angb,frbl,dbbl,phbl)
      if(iifail.eq.1) return
      if(kbad .eq. 1) return
      if((idif .eq. 2) .or. (iitp .eq. 1)) then
c: jump to right if all kr's done or if interpolation done: 
         if(kr2 .eq. nr2) goto 99
         kr1=kr2
         kr2=kr2next(inc)
         inc=inc-1
      elseif(idif .eq. 3) then
c: slide to right:  
         kr1=kr
      else
c: slide to left: 
         inc=inc + 1
         kr2next(inc)=kr2
         kr2=kr
      endif
25    if(wk(kr1) .and. wk(kr2)) then
c: jump right if end points are weak: 
         if(kr2 .eq. nr2) goto 99
         kr1=kr2
         kr2=kr2next(inc)
         inc=inc-1
         goto 25
      endif
      goto 19
99    continue
c: delete rays that were found to be weak: 
      do 45 kwk=1,nwk
         kr=krwk(kwk)
         nlist(1,kr)=nlist(1,kr) - 1
45    continue
c
199   return
      end 
