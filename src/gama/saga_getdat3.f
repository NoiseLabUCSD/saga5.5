      subroutine saga_getdat3(nline)
c
c: this subroutine reads the option file for the src/rec geometry.
c
      implicit integer*4(i-n)
c
      include 'common/rthcom' 
      include 'common/gamaoptions'
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
      include 'common/saga_gama_fail' 

      character*1 fd,cz,cz2
c
      common /info/ienv_info,lu_env,iiwrite
      common/iterpar/iter,iforwpop
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
            iifail=1
            return
cpln            stop
         endif
10    continue
      do 15 j=1,nzr 
cxx      zrzs=srloc(j,2)
         zrzs=xyz(kxy(j)+1,3)
         if((zrzs .lt. 0.) .or. (zrzs .gt. zsvp(nsvp))) then
            print *,'illegal receiver depth: ',j,zrzs
            print *,'water depth from .svp file = ',zsvp(nsvp)
            print *,'Check line numbers 10-12 in .opt file...'
            iifail=1
            return
cpln            stop
         endif
15    continue
c
c: read # of standard sets and # of individual ocean paths. 
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) nst,ncp0
      nst=1
      ncp0=0
      if((nst .lt. 0) .or. (ncp0 .lt. 0) .or. (nst .gt. 10) .or.
     .      (ncp0 .gt. 200)) then
         print *,'illegal number of standard or individual ocean ',
     .      'path lines specified: ',nst,ncp0
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
c
c: each standard set of ocean paths gives a range of bottom interactions
c: all paths with #b in that range are entered in ncpth by calling 
c: traduc in newdep.
c: kb1 through kb2 is the range of bottom interactions. type can be: 
c: 'r'=refracting; 'w'=waterborne (at least totally reflected at ocean
c: bottom); 'p'=bottom penetrating; 'a'=all.
c
cpln      call star(1,nline)
      ncp=0
      do 20 k=1,nst 
cpln         read(1,*,end=500,err=500) kb1(k),kb2(k),cz(k)
cpln         kb1(k)=1
cpln         kb2(k)=48
         kb1(k)=0
cpln         kb2(k)=0
cpln         kb1(k)=1
cpln         kb2(k)=5
cpln         kb2(k)=48
         kb2(k)=40
         cz(k)='a'
c: check for errors: 
         if((kb1(k) .lt. 0) .or. (kb2(k) .lt. kb1(k))) then 
            print *,'illegal standard ocean path: ',kb1(k),kb2(k),cz(k)
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
         ncp=ncp + 4*(kb2(k)-kb1(k)+1)
         if(kb1(k) .eq. 0) ncp=ncp - 2
         if(ncp .gt. 195) then
            print *,'too many ocean paths due to standard set ',
     .         'specification: ',ncp,'. Limit is 195.'
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
20    continue
c
c: read the individually specified ocean paths.  each path is char-
c: acterized by its initial direction, fd='u' or 'd', and the number
c: of times it enters the top (ktop) and bottom (kbot) sections of
c: of the ocean.
cpln      call star(1,nline)
      do 25 k=1,ncp0
cpln         read(1,*,end=500,err=500) fd(k),ktop(k),kbot(k),cz2(k)
         fd(k)='d'
         ktop(k)=1
         kbot(k)=3
         cz2(k)='a'
25    continue
      ncp=ncp + ncp0
c: for bottom loss tables, only include single bottom reflection:
      if(iiblt .eq. 1) then
         ncp=1
         ncp0=1
         nst=0
         fd(1)='d'
         ktop(1)=0
         kbot(1)=1
         cz2(1)='a'
      endif
c
c: arrays in paths and pathchr require length ncp+3: 
c     ncp2=ncp+3
c     allocate(ncpth(ncp2,3),fdir(ncp2),czprop(ncp2))
c
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) (maxtrav(jcan),jcan=1,ntot)
      do jcan=1,ntot
cpln         maxtrav(jcan)=5
         maxtrav(jcan)=2
      end do
      iihuh=0
      do 40 j=1,ntot
         if(maxtrav(j) .gt. 3) iihuh=1
         if(maxtrav(j) .le. 0) then
            print *,'maxtrav(j) cannot be <= 0 for bottom layer j=',j
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
40    continue
      if(iihuh .eq. 1) then
         write(8,140) (maxtrav(jcan),jcan=1,ntot)
140   format('warning: maxtrav is the max # of complete down-up ',
     .   'traversals'/'in a layer on EACH bottom interaction. maxtrav ',
     .   'has been set to a number > 3,'/'which may cause a large ',
     .   'number of insignificant rays to be considered.'/'maxtrav = ',
     .   42(i2,1x))
      endif
      do 940 jcan=1,ntot
         maxtrav(jcan)=2*maxtrav(jcan)
940   continue
c
c: check for loads of plots asked for:  
      if((iipl .gt. 0) .and. (nfr .gt. 15)) then
         print *,'warning: ',nfr,' prop loss plots requested. '
cc       print *,'ok:  (1=yes, 0=stop program)'
cc       read *,ii
cc       if(ii .eq. 0) stop
      elseif((iitf .gt. 0) .and. (nsr .gt. 15)) then
         if(iiwrite.gt.0)
     .    print *,'warning: ',nsr,' transfer functions requested. '
cc       print *,'ok?  (1=yes, 0=stop program)'
cc       read *,ii
cc       if(ii .eq. 0) stop
      endif
      if((iipic .ne. 0) .and. ((nzs .gt. 15) .or. (nrtot .gt. 25)))
     .   then
         print *,'warning: nzs = ',nzs,', nrtot = ',nrtot,
     .      '; many ray plots will be generated.' 
cc       print *,'ok? (1=yes, 0=stop program)'
cc       read *,ii
cc       if(ii .eq. 0) stop
      endif
      if((iibar .ne. 0) .and. (nsr .gt. 15)) then 
         print *,'warning: ',nsr,' bar plots will be generated.'
cc       print *,'ok? (1=yes, 0=stop program)'
cc       read *,ii
cc       if(ii .eq. 0) stop
      endif
      if((iieig .ne. 0) .and. (nsr .gt. 25)) then 
         print *,'warning: ',nsr,' eigenray lists will be generated. '
cc       print *,'ok? (1=yes, 0=stop program)'
cc       read *,ii
cc       if(ii .eq. 0) stop
      endif
c
c: read surface loss model to be used:  
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) iisl
      iisl=0
cpln      backspace(1)
      if((iisl .eq. 1) .or. (iisl .eq. 2)) then
cpln         read(1,*,end=500,err=500) iisl,swind
         iisl=0
         swind=0
c: EKW fixed bug (6-21-93) sigma->ssigma in two spots:
         ssigma=1.4102686e-3*swind**2
      elseif(iisl .eq. 4) then
cpln         read(1,*,end=500,err=500) iisl,swind
         iisl=0
         swind=0
         ssigma=10.**(-1.*abs(swind)/20.) 
      elseif(iisl .eq. 0) then
cpln         read(1,*,end=500,err=500) iisl,swind
         iisl=0
         swind=0.
      elseif(iisl .eq. 5) then
cpln         read(1,*,end=500,err=500) iisl,slfile
         call mname(slfile,lslf)
         slname=svpn(1:lsvp)//slfile(1:lslf)
         open(53,file=slname,status='unknown')
         nline2=0
cpln         call star(53,nline2)
cpln         read(53,*,end=600,err=600) nasl,nfsl
         close(53)
         nasl=nasl + 1
         nfsl=nfsl + 1
         swind=0.
      elseif(iisl .ne. 3) then
         print *,'illegal surface scattering model iisl: ',iisl
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      ctop=csvp(0)
      atop=1./ctop
      if(ntot .gt. 0) then
         cbas=cp2(ntot)
         ambas=1./max(cp1(ntot),cp2(ntot))
      else
         cbas=csvp(nsvp)
      endif
      abas=1./cbas 
c
c: read bottom/basement loss model to be used:  
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) iibl
      iibl=0
cpln      backspace(1)
      if((iibl .eq. 1) .or. (iibl .eq. 2)) then
cpln         read(1,*,end=500,err=500) iibl,bwind
         iibl=0
         bwind=0
         bsigma=1.4102686e-3*bwind**2
      elseif(iibl .eq. 4) then
cpln         read(1,*,end=500,err=500) iibl,bwind
         iibl=0
         bwind=0
         bsigma=10.**(-1.*abs(bwind)/20.) 
      elseif(iibl .eq. 0) then
cpln         read(1,*,end=500,err=500) iibl,bwind
         iibl=0
         bwind=0.
      elseif(iibl .eq. 5) then
cpln         read(1,*,end=500,err=500) iibl,blfile
         call mname(blfile,lblf)
         blname=blfile(1:lblf)
         open(53,file=blname,status='unknown')
         nline3=0
cpln         call star(53,nline3)
cpln         read(53,*,end=700,err=700) nabl,nfbl
         close(53)
         nabl=nabl + 1
         nfbl=nfbl + 1
         bwind=0.
      elseif(iibl .ne. 3) then
         print *,'illegal surface scattering model iibl: ',iibl
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
c
cpln      call star(1,nline)
cpln      read(1,*,end=500,err=500) iirth,thmn,thmx,rmnn,rmxx
      iirth=0
      thmn=0.
      thmx=90.
      rmnn=0.
      rmxx=35000.
cpln      backspace(1)
cpln      read(1,*,end=498,err=498) iirth,thmn,thmx,rmnn,rmxx,nrth
      nrth=0
      goto 499
cpln 498   backspace(1)
      if(iirth .eq. 3) then
         print *,'for iirth=3 (r-t matlab file output), need nrth'
         stop
      endif
499   if((iirth .lt. 0) .or. (iirth .gt. 100)) then 
         print *,'illegal iirth option in line ',nline,' of opt file' 
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      if(iirth .eq. 1 .or. iirth .eq. 2) then
         if((thmn .lt. 0.) .or. (thmx .gt. 90.) .or. (rmxx .le. rmnn)
     .         .or. (rmnn .lt. 0.)) then
            print *,'error in r-theta spec: ',thmn,thmx,rmnn,rmxx
            print *,'theta is src/rec grazing angle (0-90 deg).  ',
     .         'range cannot be negative.'
            stop
         endif
      endif
c
      if(iipic+iibar+iieig+iirth+iipltf+iiffi+iidat .eq. 0) then
         print *,'all options were chosen to be 0 (no).'
         stop
      endif
c
c: call pathfil to fill ncpth: 
      zs=srloc(1) 
      zr=xyz(1,3)
      kksr2=int(sign(1.,real(zs-zr)))
      
      call pathfil
c
c: mbpmax is the maximum number of bottom paths requested:
      nlay=ntot
cpln      mbpmax=5
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
