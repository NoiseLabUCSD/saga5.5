      subroutine discpts(kax)
c
c: this subroutine fills the array aaxis with the snell
c: invariant points at which the r(a) and r'(a) functions are to
c: be computed.  the discontinuity points are found from the svp
c: and put in adisc.  for each discontinuity point, two points
c: which incrementally bracket it are put into aaxis.  additional
c: points at conveniently spaced theta values are added between the
c: lowest and highest disc pts.
c
      implicit integer*4(i-n)
      include 'common/discpt' 
      include 'common/svp'
      include 'common/inivar' 
      include 'common/laydata'
      include 'common/bottom' 
      include 'common/pii'
c
      real adisc(76)
      integer*4 npt(3),nsample(76)
      data npt/4,40,13/
c
      jbset=0
      jtset=0
      jpset=0
      knext=1
c: adir is the maximum value of the si for which a ray can
c: travel from the source level to the receiver level.
      cbig=cmax(2)
      adir=1./cbig
      jair=1
      adisc(jair)=adir
      nsample(jair)=kseg
      jbd(-1)=jair
      jtrav(-1)=jair
      jpen(-1)=jair 
c
c: adisc(1) to adisc(jair) are the si discontinuity points in the
c: upward direction.
      jpen(-2)=jair 
      do 10 k=klev(1)-1,klev(0),-1
         if(kseg .eq. 1) then 
            call chdisc(adisc,csvp(k),cbig,jair,npt,nsample,kseg)
         elseif(kseg .eq. 3) then
            call chpeak(adisc,csvp(k),k,cbig,jair,npt,nsample,2)
         endif
10    continue
      jtrav(-2)=jair
c: check for beam displacement at ocean surface:  
      call chdisc(adisc,cs1(-1),cbig,jair,npt,nsample,3)
      call chdisc(adisc,cp1(-1),cbig,jair,npt,nsample,3)
      jbd(-2)=jair
      cbig=cmax(2)
      jsub=jair+1
      adisc(jsub)=adir
      nsample(jsub)=kseg
      jpen(0)=jsub
c
c: adisc(jair+1) to adisc(jsub) are the si disc pts in the
c: downward direction.
      do 20 k=klev(2)+1,klev(3)
         if(kseg .eq. 1) then 
            call chdisc(adisc,csvp(k),cbig,jsub,npt,nsample,kseg)
         elseif(kseg .eq. 3) then
            call chpeak(adisc,csvp(k),k,cbig,jsub,npt,nsample,2)
         endif
20    continue
c
      jtrav(0)=jsub 
c: check for beam displacement disc points at ocean bottom: 
      call chdisc(adisc,cs1(1),cbig,jsub,npt,nsample,3)
      call chdisc(adisc,cp1(1),cbig,jsub,npt,nsample,3)
      jbd(0)=jsub
c
c: adisc(jpen(k)) is the si for which a ray barely penetrates the
c: bottom layer k.  adisc(jtrav(k)) is the si for which a ray
c: barely traverses the bottom layer k. 
      do 30 k=1,ntot
c: if no mismatch at top of layer, set jpen to last real interface: 
c: 12/11/91:         if(mismch(k-1) .eq. -1) then
c        if(mismch(k-1) .ne. 1) then
c: BUG: 9/8/92: Looks like it was right to begin with:
c: See setbp: no-mismatch layers are treated as one only if gradients 
c: are cont
         if(mismch(k-1) .eq. -1) then
            do 32 kk=k-1,1,-1 
               jpen(k)=jpen(kk)
c: 12/11/91:               if(mismch(kk-1) .ne. -1) goto 33
c              if(mismch(kk-1) .eq. 1) goto 33
c: BUG: 9/8/92: Looks like it was right to begin with:
               if(mismch(kk-1) .ne. -1) goto 33
32          continue
33          continue
c     print *,'discpts mismch=-1 for layer k= ',k,1./adisc(jpen(k))
         else
            jpen(k)=jsub
         endif
c: if no mismatch at bottom of layer, don't put in disc point: 
         if(mismch(k) .eq. -1) then
            jtrav(k)=jsub
         else
            call chdisc(adisc,cp2(k),cbig,jsub,npt,nsample,1)
            jtrav(k)=jsub
c: check for beam disp discontinuity points at bottom of layer: 
            call chdisc(adisc,cs1(k+1),cbig,jsub,npt,nsample,3)
            call chdisc(adisc,cp1(k+1),cbig,jsub,npt,nsample,3)
         endif
         jbd(k)=jsub
c     print *,'discpts: ',k,1./adisc(jbd(k)),1./adisc(jtrav(k)),
c    .   1./adisc(jpen(k))
30    continue
c
c: theps is the theta increment used around the discontinuity points: 
      theps=1.e-9
c
c: make sure that the number of sample points will fit into the array: 
79    nstot=3
      do 80 j=1,jsub
         nstot=nstot + npt(nsample(j)) + 2
80    continue
      if(nstot .gt. 156) then 
         print *,'# sample point in snell inv axis being reduced.'
         do 85 j=1,3
            if(npt(j) .gt. 1) npt(j)=npt(j)-1
85       continue
         goto 79
      endif
c
c: for disc pts not caused by bd, space theta values evenly.
      do 89 j=1,2
         do 90 k=1,npt(j)
            space(j,k)=float(k)/float(npt(j)+1)
90       continue
89    continue
c: for disc pts caused by bd, abrupt bumps in the r(a) function
c: can occur, so bias the sample points toward the extremes of
c: the interval.
      den=.5*(float(npt(3)) + 1.)**2
      do 91 k=1,(npt(3)+1)/2
         space(3,k)=float(k)**2/den
91    continue
      do 92 k=(npt(3)+1)/2 + 1,npt(3)
         space(3,k)=1. - space(3,npt(3)-k+1)
92    continue
c
c: aaxis(1) is close to zero.  kax is the counter for aaxis.
      aaxis(1)=sin(theps)/cs
      kax=1
      jfall=jsub
c: kpt is the counter for kdisc.
      kpt=0
c: nbd, ntrav, and npen are counters for jbd, jtrav, and jpen, resp.
      nbd=ntot
      ntrav=ntot
      npen=ntot
c
c: we now merge the two lists (the upward and downward lists) of
c: disc pts.  for each disc pt, two points that incrementally
c: bracket it are entered into aaxis.  additional points are
c: entered between disc pts at convenient theta spacings.
c: kdisc holds information about which points in aaxis are disc
c: pts and which direction they come from.  the kpt'th disc pt
c: is stored in aaxis(int(kdisc(kpt)/10)).  let kdir=
c: mod(kdisc(kpt,10)).  if kdir=2 then the disc pt is from the
c: downward direction; if kdir=3 then it's from the upward; and
c: if kdir=1 then it's from both directions.
      do 50 jrise=jair,1,-1
51       if(adisc(jfall) .lt. adisc(jrise)) then
            call interp(adisc(jfall),npt,nsample(jfall),kax,cs) 
            kpt=kpt+1
            kdisc(kpt)=kax*10 + 2
            call laycrit(nbd,ntrav,npen,jfall)
            jfall=jfall-1
            goto 51 
         elseif(adisc(jrise) .lt. adisc(jfall)) then
            call interp(adisc(jrise),npt,nsample(jrise),kax,cs) 
            kpt=kpt+1
            kdisc(kpt)=kax*10 + 3
            call laycrit(nbd,ntrav,npen,jrise)
         elseif(adisc(jfall) .eq. adisc(jrise)) then
            nsam=nsample(jrise)
            if(npt(nsample(jfall)) .gt. npt(nsample(jrise))) nsam=
     .            nsample(jfall)
            call interp(adisc(jrise),npt,nsam,kax,cs) 
            kpt=kpt+1
            kdisc(kpt)=kax*10 + 1
            call laycrit(nbd,ntrav,npen,jfall)
            call laycrit(nbd,ntrav,npen,jrise)
            jfall=jfall-1
         endif
50    continue
c
      inx=kax-1
      jbd(-1)=kpt
      jtrav(-1)=kpt 
      jpen(-1)=kpt
cxx   print *,'1./aaxis: ',(1./aaxis(k),k=1,inx)
c    print *,(asin(aaxis(k)*cs)/pierad,k=1,inx),inx
c    print *,'  '
c     print *,'kdisc, etc: ',(kdisc(k),k=1,kpt),(jbd(k),k=-2,ntot),
c    .c     (jtrav(k),k=-2,ntot),(jpen(k),k=-2,ntot)
c     print *,'nstot = ',nstot
c
      return
      end 
