      subroutine getdat3tilt
c
c: this subroutine reads the option file for the src/rec geometry.
c
      implicit integer*4(i-n)
      include 'common/rthcom' 
      include 'common/options'
      include 'common/pcodes' 
      include 'common/pcode2' 
      include 'common/srlox'
      include 'common/svp'
      include 'common/depth'
      include 'common/paths'
      include 'common/caustix'
      include 'common/pii'
      include 'common/freqcom'
      include 'common/bottom' 
      include 'common/scatcom'
      include 'common/charcom' 
      character*1 fd,cz,cz2
      common /info/ienv_info,lu_env,iiwrite
      integer ienv_info,lu_env,iiwrite
c
c: no limit now on number of ranges:
c     if(nrtot .gt. 2000) then
c        print *,'total # ranges too high (nrtot<=2000): ',nrtot
c        stop
c     endif
c: check source/receiver depths: 
      do 10 j=1,nzs 
cxx      zrzs=srloc(j,1)
         zrzs=srloc(j)
         if((zrzs .lt. 0.) .or. (zrzs .gt. zsvp(nsvp))) then
            print *,'illegal source depth: j,zs = ',j,zrzs
            print *,'water depth from .svp file = ',zsvp(nsvp)
            print *,'Check line number 6 in .opt file...'
            stop
         endif
10    continue
      do 15 j=1,nzr 
cxx      zrzs=srloc(j,2)
         zrzs=xyz(kxy(j)+1,3)
         if((zrzs .lt. 0.) .or. (zrzs .gt. zsvp(nsvp))) then
            print *,'illegal receiver depth: ',j,zrzs
            print *,'water depth from .svp file = ',zsvp(nsvp)
            print *,'Check line numbers 10-12 in .opt file...'
            stop
         endif
15    continue
c
c: call pathfil to fill ncpth: 
      zs=srloc(1) 
      zr=xyz(1,3)
      kksr2=int(sign(1.,real(zs-zr)))
      call pathfil
c
c: mbpmax is the maximum number of bottom paths requested:
      nlay=ntot
      mbpmax=2
      do 52 icpth=1,ncp
         nbotx=ncpth(icpth,3)/2
         jdpsum=0
         do 54 jdp=1,nlay
            jprod=1 
            do 56 j=1,jdp
               jprod=jprod*nbotx*maxtrav(j)/2
56          continue
            jdpsum=jdpsum + jprod
54       continue
         mbpmax=mbpmax + 1 + jdpsum
52    continue
      mbpmax=mbpmax*nzr
      if(mbpmax .gt. 100000) then
         if(iiwrite.gt.0)
     .    print *,'limiting mbpmax = ',mbpmax,' to 100000.'
         mbpmax=100000
      endif
c
c: mray is the number of rays allocated in xlist,klist,etc.
      mray=50
c
      return
      end 
