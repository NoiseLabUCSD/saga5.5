      subroutine newrdep(range)
c
c: this subroutine sets a new receiver depth.  the svp arrays are
c: adjusted and the scratch files for eigenray data are reinitialized.
c
      implicit integer*4(i-n)
      include 'common/depth'
      include 'common/srlox'
      include 'common/svp'
      include 'common/caustix'
      include 'common/scatcom'
      include 'common/bottom'
      include 'common/gamaoptions'
      include 'common/hccom'
      include 'common/bdcom'
      include 'common/cbotcom'
      real range(nrtot)
      logical qmax,qupper,qlower,qsrup,qsrlow
c
      zr=xyz(kxy(jzr)+1,3)
cpln      write(6,*)'Depth: ',jzr,nzr,zr
c
      mmr1=mr1(jzr) 
      mmr2=mr2(jzr) 
      m21=mmr2+1
      m22=mmr2+2
      nnxy=nxy(jzr) 
      
c
      rewind(56)
      read(56) nsvp,zsvp,csvp,bsvp,g0,g1,g2
c
      ztol=(csvp(1)/fmax)/100.
      if((zs .eq. zr) .and. (iiblt .ne. 1)) then
c: NEW 4-24-90: zs=zr can inhibit detection of horizontal caustics,
c: so shift zr a tiny amount:
         if(zr .gt. ztol) then
            zr=zr - ztol
         else
            zr=zr + ztol
         endif
         do 490 jcan=kxy(jzr)+1,kxy(jzr)+nnxy
            xyz(jcan,3)=zr
490      continue
         print *,' '
         print *,'rec depth ',jzr,' changed by .01*(min lambda) =',ztol,
     .      ' so that zs not equal to zr: zs,zr = ',zs,zr
      endif
c
c: refill ncpth() if relative positions of source and rec have changed.
      kksr1=int(sign(1.,real(zs-zr)))
cxx   print *,'start of newrdep: ',jzr,zr,kksr1,kksr2
      if(kksr1 .ne. kksr2) call pathfil 
      kksr2=kksr1
c: add extra layers for source and receiver depths if necessary
      ks=-1
      kr=-1
      do 22 k=0,nsvp
         if(ks .eq. -1) then
            if(zs .le. zsvp(k)) ks=k
         endif
         if(kr .eq. -1) then
            if(zr .le. zsvp(k)) kr=k
         endif
22    continue
c
      if(zs .ge. zr) then
cxx      if(zs .ne. zsvp(ks)) call srlev(zs,ks)
         if(abs(zs - zsvp(ks)) .gt. .5*ztol) then
            call srlev(zs,ks)
         else
            zs=zsvp(ks)
         endif
         if(abs(zr - zsvp(kr)) .gt. .5*ztol) then
            call srlev(zr,kr) 
            ks=ks+1 
         else
            zr=zsvp(kr)
         endif
      else
         if(abs(zr - zsvp(kr)) .gt. .5*ztol) then
            call srlev(zr,kr)
         else
            zr=zsvp(kr)
         endif
         if(abs(zs - zsvp(ks)) .gt. .5*ztol) then
            call srlev(zs,ks) 
            kr=kr+1 
         else
            zs=zsvp(ks)
         endif
      endif
c
      klev(0)=0
      klev(1)=min0(ks,kr)
      klev(2)=max0(ks,kr)
      klev(3)=nsvp
      zlev(1)=0.
      zlev(2)=min(zs,zr)
      zlev(3)=max(zs,zr)
      cs=csvp(ks)
      cr=csvp(kr)
      csarr(jzs)=cs
      crarr(jzr)=cr
      do 23 jbdf=1,nbdf
         bdlam(jbdf)=cr/bdf(jbdf)
23    continue
      atol=(1./cs)*1.e-6
      rtol=(cr/fmax)/1500.
c
      do 25 j=1,3
         cmax(j)=0. 
         do 30 k=klev(j-1),klev(j)
            cmax(j)=max(cmax(j),csvp(k))
30       continue
25    continue
      artop=1./cmax(1)
      arbot=1./cmax(3)
c: EKW 5-3-93: needed for rtc_calc:
      aturn(1)=artop
      aturn(2)=arbot
      amtop=artop     
      if(ntot .eq. 0) ambas=arbot
c
c: fill kseq(2,3,3), which tells the correct sequence of ocean layers 
c: for each section of the ocean when tracing a ray path.
      do 70 leg=1,3 
         kseq(2,leg,1)=klev(leg-1)
         kseq(2,leg,2)=klev(leg)-1
         kseq(2,leg,3)=1
         kseq(1,leg,1)=klev(leg)
         kseq(1,leg,2)=klev(leg-1)+1
         kseq(1,leg,3)=-1
70    continue
c
      qmax=((cmax(2) .gt. cs) .and. (cmax(2) .gt. cr))
      qupper=(csvp(klev(1)) .lt. cmax(1))
      qlower=(csvp(klev(2)) .lt. cmax(3))
      qsrup=(csvp(klev(1)) .ge. csvp(klev(2)))
      qsrlow=(csvp(klev(2)) .ge. csvp(klev(1)))
c: variables needed for horizontal caustic check in rays:
      if(qmax) then
         nldif=0
      else
         if(qupper .and. qsrup) then
            nldif=1
         elseif(qlower .and. qsrlow) then
            nldif=3
         else
            nldif=0
         endif
      endif
      nlsam=4 - nldif
      krish(1)=0
      krish(2)=0
c
      return
      end 
