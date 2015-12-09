      subroutine geomtilt(rangx,tyme,range,jjnr,kkk,xxx,yyy,phi)
c
c: this subroutine computes the source track and the ranges to the
c: receiver array.
c
      implicit integer*4(i-n)
      include 'common/pii'
      include 'common/srlox'
      include 'common/options'
      include 'common/depth'
      include 'common/freqcom'
      include 'common/traxcom'
      real rangx(nrangx),tyme(nrangx),range(nrtot)
      integer*4 jjnr(nrtot),kkk(nrtot),jzr1
      real xxx(nrnf),yyy(nrnf)
      real*4 phi(nrtot)
c
      if(iirr .eq. 1) then
c: fill rangx for simple range specification:
         do 120 j=1,nrangx
            rangx(j)=1000.*rtmp(j)
            tyme(j)=rangx(j)
            xxx(j)=rangx(j)
            yyy(j)=0.
 120     continue
      elseif(iirr .eq. -1) then
         do 121 j=1,nrangx
            rangx(j)=1000.*(rtmp(1) + float(j-1)*drtmp)
            tyme(j)=rangx(j)
            xxx(j)=rangx(j)
            yyy(j)=0.
 121     continue 
      end if
c     print *,'tyme = ',(tyme(j),j=1,nrangx)
c
299   continue
      nr=0
c: calculate ranges at each depth: 
      do 21 jzr1=1,nzr
c: for the jzr'th receiver depth, the ranges are in range(mr1:mr2).
         mr1(jzr1)=nr + 1
         do 22 kr=1,nrangx
            do 24 jxy=1,nxy(jzr1)
c: compute s/r ranges for receivers (jx,jy) at each source pos kr
               nr=nr+1
               delx=xxx(kr) - xyz(kxy(jzr1)+jxy,1)
               dely=yyy(kr) - xyz(kxy(jzr1)+jxy,2)
               range(nr)=sqrt(delx**2 + dely**2)
               if(range(nr) .ne. 0.) then
                  phi(nr)=atan2(dely,delx)
               else
                  phi(nr)=0.
               endif
c               write(6,*)'range(nr): ',range(nr)
c               write(6,*)'delx,dely: ',delx,dely
               if(range(nr) .lt. 1.e-3) range(nr)=.01
c: given jzr,kr,jxy, set mr=jjnr(mr1(jzr) + nxy(jzr)*(kr-1) + jxy - 1)
c: and range(mr) is the range corr to that source and rec pos.
               kkk(nr)=nr
24          continue
22       continue
         mr2(jzr1)=nr
21    continue
      if(iitrx .ne. 0) call traxpl(xxx,yyy,tyme,rangx,nrangx)
c
c: bubble sort range(), keeping jjxy up to date also: 
      do 30 jzr1=1,nzr
         do 32 nr1=1,mr2(jzr1)-mr1(jzr1)
            do 35 nr2=mr1(jzr1),mr2(jzr1) - nr1
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
         do 40 kr=mr2(jzr1)-1,mr1(jzr1),-1
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
