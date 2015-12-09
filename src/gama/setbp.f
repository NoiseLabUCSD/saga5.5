      subroutine setbp(kst,kstt,kend,kendd,norefl,cz,kdone,mmr0b,
     .   ntlay,NLAYM)
c
      implicit integer*4(i-n)
      include 'common/bottom' 
      include 'common/pathway'
      include 'common/paths'
      include 'common/discpt' 
      include 'common/weakcom'
      include 'common/srlox'
      include 'common/gamaoptions'
cpln      integer*4 ntlay(-1:NLAYM+1,5000)
c added by pln 11/02/00, otherwise array index 
c exceeded at line 179
      integer*4 ntlay(-1:NLAYM+2,50000)
      character*1 cz
c
cpln      write(6,*)'INSIDE SETBP NLAYM: ',NLAYM
55    kdone=0
      kmult=0
      norefl=0
      mmr0=mmr0b
      if(icpth .ne. icxx) then
c: if a new ocean path, reset ndp to 0 and return:
         ndp=0
         nttot(1)=0
c: if no mismatch at deepest layer, flag for use in rays:
         if(mismch(0) .le. 0) norefl=1
         goto 99
      elseif(ndp .eq. 0) then
c: if the second bottom path for the current ocean path (the first was
c: the ndp=0 one), set ndp=1
         ndp=1
88       if((ndp .gt. ntot) .or. (cz .eq. 'w') .or. (cz .eq. 'r')) then
cxx   print *,'kdone=1 in setbp: ndp,ntot,cz = ',ndp,ntot,cz
            kdone=1
            return
         endif
         if(mismch(ndp) .eq. -1) then
c: if no impedance mismatch at deepest layer, increment ndp:
            ndp=ndp + 1
            goto 88
         endif
         do 940 jcan=1,ndp-1
            nttot(jcan)=2
940      continue
         do 942 jcan=ndp,ntot+1
            nttot(jcan)=0
942      continue
      endif
      jptr=ndp
10    nttot(jptr)=nttot(jptr) + 2
      if(nttot(jptr) .gt. nbotx*maxtrav(jptr)) then
         if(jptr .eq. 1) then
c: all combinations for ndp are done, so increment ndp:
12          if(ndp .eq. ntot) then
c: all combinations for the current icpth are done, so return:
cxx   print *,'kdone=1: ndp,ntot,jptr,nttot,maxtrax = ',
cxx  .   ndp,ntot,jptr,nttot,maxtrax
               kdone=1
               return
            endif
            ndp=ndp + 1
c: ignore bottom path if there is no mismatch and sv gradient is
c: continuous (mismch=-1) at the bottom of the deepest layer.
c: the rays that refract in these layers are found along with those
c: that refract in the next deepest layer that does have mismatch
c: (see subroutine discpt):
            if(mismch(ndp) .eq. -1) goto 12
            jptr=ndp
            do 944 jcan=1,ndp-1
               nttot(jcan)=2
944         continue
            do 946 jcan=ndp,ntot+1
               nttot(jcan)=0
946         continue
         else
            do 948 jcan=jptr,ndp
               nttot(jcan)=2
948         continue
            jptr=jptr - 1
         endif
         goto 10
      endif
c
c: check if current bottom path has a reflection off a 
c: no-mismatch layer:
cxx   do 15 j=1,ndp-1
c: NEW 10-10-90: DETECT NO-MISMATCH AT WATER-SEDIMENT INTERFACE ALSO:
      do 15 j=0,ndp-1
         if(nttot(j) .gt. 2) kmult=1
         if(mismch(j) .le. 0) then
c: if reflections must occur at a no-mismatch layer, ignore bottom path:
            if(nttot(j) .ne. nttot(j+1)) then
               goto 55
            endif
         endif
15    continue
      if(nttot(ndp) .gt. 2) kmult=1
c: if no mismatch at deepest layer, flag for use in rays:
      if(mismch(ndp) .le. 0) norefl=1
c
c: check if current bottom path must be weaker than an already weak one:
      do 20 jwk=nbwk,1,-1
         mw=mweak(jwk,1)
c: to compare, require ndp to be the same:
         if(ndp .ne. ntlay(-1,mw)) goto 20
c: skip over bottom path if any Mj<Wj:
         do 30 j=1,ndp-1
            if(nttot(j) .lt. ntlay(j,mw)) goto 20
            if(nttot(j) .gt. ntlay(j,mw)) then
               if((ntlay(j,mw) .lt. ntlay(j-1,mw)) .or.
     .            (ntlay(j,mw) .lt. ntlay(j+1,mw))) goto 20
            endif
30       continue
c: if fewer bottom interactions here than in weak, discard:
         if(nttot(0) .lt. ntlay(0,mw)) goto 20
c: if more bottom interactions here and double bounces in weak path,
c: discard:
         if(nttot(0) .gt. ntlay(0,mw)) then
            if(ntlay(0,mw) .lt. ntlay(1,mw)) goto 20
         endif
         if(ndp .ne. 0) then
c: if fewer bounces off last layer here, discard:
            if(nttot(ndp) .lt. ntlay(ndp,mw)) goto 20
c: if more bounces off last layer here and double bounces in weak path,
c: discard:
            if(nttot(ndp) .gt. ntlay(ndp,mw)) then
               if(ntlay(ndp,mw) .lt. ntlay(ndp-1,mw)) goto 20
            endif
         endif
c: if a weak bottom path found, update mmr0:
         continue
c        print *,'mmr0 updated: jwk,mmr0,krm= ',jwk,mmr0,mweak(jwk,2)
         mmr0=max0(mmr0,mweak(jwk,2))
c     print *,'current nttot: ',nttot(0:ndp)/2,ntop,nbotx
c     print *,'weak nttot:    ',ntlay(0:ndp,mw)/2
         if(mmr0 .gt. mmr2) then
c           print *,'ignoring bottom path '
            goto 55
         endif
20    continue
c
c: nbas is the number of times the ray hits the basement (substrate):
99    nbas=nttot(ndp)/2
c
c: set kstt and kendd for the entire path.
      if(ndp .gt. 0) then
         kstt=min0(kst,jbd(ndp))
         kendd=min0(kend,jpen(ndp))
c: BUG 10-19-90: Don't allow rays to reflect off no-mismatch sed int:
cqq   elseif(norefl .eq. 1) then
cqq      kstt=max0(kst,jtrav(0))
cqq      kendd=kend
      else
         kstt=kst
         kendd=kend 
      endif
c: BUG 01-17-91: Don't allow rays to reflect off no-mismatch int:
      if(norefl .eq. 1) then
         kstt=max0(kst,jtrav(ndp))
cxx      print *,'norefl=1: ',ntop,nbotx,nttot(0:ndp),kstt
      endif
c
      icxx=icpth
      if(mbp .ge. mbpmax) then
         print *,'# bottom path combinations exceeded allocated space.',
     .      'mbpmax = ',mbpmax
         print *,'mbpmax = ',mbpmax,' must be increased.'
         kdone=1
         return
      endif
      mbp=mbp + 1
      ntlay(-1,mbp)=ndp
cxx   if(iidiag .ne. 0) print *,'setbp set mbp: mbp,ndp,ntlay = ',
cxx  .   mbp,ndp,ntlay(-1,1:mbp)
cpln      write(6,*)'NDP: ',ndp
      do 80 j=0,ndp+1
         ntlay(j,mbp)=nttot(j)
80    continue
      if(iidiag .ne. 0) print *,'setbp done: mbp,ntop,nbotx,nttot = ',
     .   mbp,ntop,nbotx,(nttot(jcan)/2,jcan=1,ndp)
c
      return
      end
