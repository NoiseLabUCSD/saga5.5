      subroutine geom(rangx,tyme,range,jjnr,kkk,xxx,yyy,phi)
c
c: this subroutine computes the source track and the ranges to the
c: receiver array.
c
      implicit integer*4(i-n)
      include 'common/pii'
      include 'common/srlox'
      include 'common/gamaoptions'
      include 'common/depth'
      include 'common/freqcom'
      include 'common/traxcom'
      real*8 rangx(nrangx),tyme(nrangx),range(nrtot)
      integer*4 jjnr(nrtot),kkk(nrtot)
      real*8 xxx(nrnf),yyy(nrnf)
      real*4 phi(nrtot)
c
c: if bottom loss table requested, ignore the range spec in opt file:
      if(iiblt .eq. 1) then
c        h2=2.*zlev(4)
         h2=2.*zlev(4) - srloc(1) - xyz(1,3)
         if(nang .gt. 1) then
            dthg=(ang2 - ang1)/float(nang-1)
         else
            dthg=0.
         endif
         do 300 jang=1,nang
            thg=max(.05,min(89.95,ang1 + float(jang-1)*dthg))*pierad
            rangx(jang)=h2/tan(thg)
            tyme(jang)=rangx(jang)
            xxx(jang)=rangx(jang)
            yyy(jang)=0.
300      continue
      elseif(iirr .eq. 1) then
c: fill rangx for simple range specification:
         do 120 j=1,nrangx
            rangx(j)=1000.*rtmp(j)
            tyme(j)=rangx(j)
            xxx(j)=rangx(j)
            yyy(j)=0.
120      continue
      elseif(iirr .eq. -1) then
         do 121 j=1,nrangx
            rangx(j)=1000.*(rtmp(1) + float(j-1)*drtmp)
            tyme(j)=rangx(j)
            xxx(j)=rangx(j)
            yyy(j)=0.
121      continue 
c         write(6,*)'geom: '
c         write(6,*)'rtmp(1), drtmp: ',rtmp(1),drtmp
c         write(6,*)'rng: ',(rangx(j),j=1,nrangx)
c         write(6,*)
      elseif(iirr .eq. 0) then
c: convert source track specification:
         r0=0.
         tlast=0.
         krr=0
         do 10 nl=1,nleg
            if(iitype(nl) .eq. 1) then
               cosphi=cos(phid(nl)*pierad)
               sinphi=sin(phid(nl)*pierad)
c: for V=0, times are distances; for V.ne.0, convert to min to km
               fac=1.
               if(vs(nl) .ne. 0.) fac=vs(nl)*60./1000.
               if(iicont(nl) .eq. 0) then
                  sg=sign(1.,abs(phid(nl))-90.)
c: xxx0,yyy0 are (x,y) source coordinates at cpa: 
                  xxx0=sg*cpa(nl)*sinphi
                  yyy0=-sg*cpa(nl)*cosphi
                  x1x=xxx0 + fac*t1(nl)*cosphi
                  y1x=yyy0 + fac*t1(nl)*sinphi
                  x2x=xxx0 + fac*t2(nl)*cosphi
                  y2x=yyy0 + fac*t2(nl)*sinphi
                  time0=t1(nl)
               else
c: begin at kr=2 for continuous legs so there is no overlap of points: 
                  x1x=x2x
                  y1x=y2x
                  x2x=x1x + fac*(t2(nl)-t1(nl))*cosphi
                  y2x=y1x + fac*(t2(nl)-t1(nl))*sinphi
                  time0=tlast
               endif
            else
               if(iicont(nl) .eq. 0) then
                  x1x=x1(nl)
                  y1x=y1(nl)
                  time0=t1(nl)
               else
                  x1x=x2x
                  y1x=y2x
                  time0=tlast
               endif
               x2x=x2(nl)
               y2x=y2(nl)
            endif
            nrden=max0(1,nrleg(nl) - 1)
            delx=(x2x - x1x)/nrden
            dely=(y2x - y1x)/nrden
            delr=sqrt(delx**2 + dely**2)
            delt=(t2(nl) - t1(nl))/nrden
            do 12 kr=iicont(nl)+1,nrleg(nl)
               krr=krr + 1
               xxx(krr)=1000.*(x1x + (kr-1)*delx)
               yyy(krr)=1000.*(y1x + (kr-1)*dely)
c: compute source track distance: 
               rangx(krr)=r0 + 1000.*(kr-1)*delr
c: tyme is the current time (or distance along track for vs=0): 
               tyme(krr)=time0 + (kr-1)*delt
12          continue
            tlast=tyme(krr)
            r0=rangx(krr)
10       continue
      else
         print *,'bug in getdat2 or geom: iirr illegal: ',iirr
         stop
      endif
c     print *,'rangx = ',nrangx,(rangx(j),j=1,nrangx)
c     print *,'tyme = ',(tyme(j),j=1,nrangx)
c
299   continue
      nr=0
c: calculate ranges at each depth: 
      do 21 jzr=1,nzr
c: for the jzr'th receiver depth, the ranges are in range(mr1:mr2).
         mr1(jzr)=nr + 1
         do 22 kr=1,nrangx
            do 24 jxy=1,nxy(jzr)
c: compute s/r ranges for receivers (jx,jy) at each source pos kr
               nr=nr+1
               delx=xxx(kr) - xyz(kxy(jzr)+jxy,1)
               dely=yyy(kr) - xyz(kxy(jzr)+jxy,2)
               range(nr)=sqrt(delx**2 + dely**2)
               if(range(nr) .ne. 0.) then
                  phi(nr)=atan2(dely,delx)
               else
                  phi(nr)=0.
               endif
               if(range(nr) .lt. 1.e-3) range(nr)=.01
c: given jzr,kr,jxy, set mr=jjnr(mr1(jzr) + nxy(jzr)*(kr-1) + jxy - 1)
c: and range(mr) is the range corr to that source and rec pos.
               kkk(nr)=nr
24          continue
22       continue
         mr2(jzr)=nr
21    continue
      if(iitrx .ne. 0) call traxpl(xxx,yyy,tyme,rangx,nrangx)
c
c: bubble sort range(), keeping jjxy up to date also: 
      do 30 jzr=1,nzr
         do 32 nr1=1,mr2(jzr)-mr1(jzr)
            do 35 nr2=mr1(jzr),mr2(jzr) - nr1
               nr2p1=nr2+1
               if(range(nr2) .gt. range(nr2p1)) then
                  rtemp=range(nr2p1)
                  range(nr2p1)=range(nr2)
                  range(nr2)=rtemp
                  jtemp=kkk(nr2p1)
                  kkk(nr2p1)=kkk(nr2)
                  kkk(nr2)=jtemp
               endif
35          continue
32       continue
c: make sure no ranges are equal, (otherwise quadratic fit can bomb): 
c: BUG 4-10-90: add eps to ranges so that they can't go negative:
         eps0=.00001
         eps=eps0
c: 5-13-90: BUG: check for ranges being within eps/2 of each other:
cxx   print *,'ranges = ',range(mr1(jzr):mr2(jzr))
         do 40 kr=mr2(jzr)-1,mr1(jzr),-1
            if(abs(range(kr) - range(kr+1)) .lt. eps/2.) then
               range(kr+1)=range(kr+1) + eps
cxx   print *,'range increased: ',kr,kr+1,range(kr),range(kr+1)
               eps=.75*eps
            else
               eps=eps0
            endif
40       continue
30    continue
      do 50 nr=1,nrtot
         jj=kkk(nr) 
         jjnr(jj)=nr
50    continue
c     print *,'nzr,mr1,mr2 = ',nzr,mr1(1:nzr),mr2(1:nzr)
c     print *,'range: ',nrtot,(range(j),j=1,nrtot)
c     print *,'jjnr: ',(jjnr(j),j=1,nrtot)
c     do 99 jzr=1,nzr
c        do 98 jxy=1,nxy(jzr) 
c           print *,'r(jjnr) =  ',jzr,jxy,(range(jjnr(mr1(jzr) +
c    .         nxy(jzr)*(kr-1) + jxy - 1)),kr=1,nrangx)
98       continue
99    continue
c: nsr is the total number of source/receiver pairs: 
      nsr=nzs*nrtot 
c
      return
      end 

