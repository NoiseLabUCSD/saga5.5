      subroutine bdsort(kl,xl,nray,nr,jcaus,mbp,kbdf1,kbdf2,nreig,
     .   nreig2)
c
c: this subroutine sorts through the beam displaced rays at various
c: frequencies and finds the ones that go together.
c
      implicit integer*4(i-n)
      include 'common/bdcom'
      include 'common/bdchar'
      integer*4 kl(2,2500)
      real xl(7,2500),another(10)
c
c: nreig points to the ray at the frequency for picture, list, etc:
      nreig=0
      nreig2=0
      if(kbdf .eq. nfeig) nreig=nr
c: fill the eigenray characteristics for frequency bdf(kbdf):
      call bdfill(kl,xl,nr,kbdf,gbd,tbd,ebd,trmbd,trpbd)
      if(jcaus .eq. 1) then
c: for shadow rays, fill additional characteristics:
         nr2=nr+1
         xaibd(kbdf)=xl(1,nr2)
         s1fbd(kbdf)=xl(3,nr2)
         s2fbd(kbdf)=xl(4,nr2)
         kl(1,nr2)=kl(1,nr2) - kbdf*1000000
      elseif(jcaus .eq. 2) then
c: for insonified zone pairs, fill second of two rays:
         call bdfill(kl,xl,nr+1,kbdf,g2bd,t2bd,e2bd,trm2bd,trp2bd)
      endif
c
      kbdf1=kbdf
      kbdf2=kbdf
      kk=kbdf2 + 1
      jc(kbdf)=jcaus
      jc2=jcaus
      jc1=jc2
      nother=0
      a0=xl(1,nr)
c: loop through other rays found at other frequencies:
      do 110 nrr=nr+1,nray
         if(jc2 .ne. 0) then
            jc2=0
            goto 110
         endif
c: make sure ray has same path:
         mbp2=kl(2,nrr)/1000
         if(mbp2 .ne. mbp) goto 112
         jc2=mod(kl(1,nrr)/100,10)
         kbdf=kl(1,nrr)/1000000
         if(kbdf .eq. 0) then
c: kbdf=0 means the ray has already been processed:
            goto 110
         elseif(kbdf .eq. kbdf2) then
c: keep track of other rays at same frequency to make sure assoc ok.
            nother=nother + 1
            another(nother)=xl(1,nrr)
            goto 110
         elseif(kbdf .gt. kk) then
c: if kbdf > kk, then the ray we are looking for isn't there:
            goto 112
         endif
c: make sure another ray doesn't match up better than current one:
         acur=xl(1,nrr)
c: for doublet pair, take acur to be that of stronger ray:
         if(jc2 .eq. 2) then
            if(abs(xl(2,nrr+1)) .gt. abs(xl(2,nrr))) acur=xl(1,nrr+1)
         endif
         adif0=abs(acur - a0)
         do 15 kother=1,nother
            adif=abs(acur - another(kother))
            if(adif .lt. adif0) goto 112
15       continue
c
c: correct ray found at frequency kbdf:
         if(kbdf .eq. kk) then
            kbdf2=kbdf
            kk=kbdf2 + 1
            if(kbdf .eq. nfeig) nreig=nrr
            nother=0
            a0=acur
            jc(kbdf)=jc2
c: fill in eigenray characteristics as before:
            call bdfill(kl,xl,nrr,kbdf,gbd,tbd,ebd,trmbd,trpbd)
c     print *,'kbdf2 incr: kbdf2,jc1,jc2 = ',kbdf2,jc1,jc2
            if(jc2 .eq. 1) then
               nrr2=nrr+1
               xaibd(kbdf)=xl(1,nrr2)
               s1fbd(kbdf)=xl(3,nrr2)
               s2fbd(kbdf)=xl(4,nrr2)
               kl(1,nrr2)=kl(1,nrr2) - kbdf*1000000
            elseif(jc2 .eq. 2) then
               call bdfill(kl,xl,nrr+1,kbdf,g2bd,t2bd,e2bd,trm2bd,
     .            trp2bd)
            elseif((jc1 .gt. 0) .and. (jc2 .eq. 0)) then
c: check for shadow or dublet changing to two normal rays:
               nrr2=nrr+1
               if(kbdf .eq. nfeig) nreig2=nrr2
               jc0=mod(kl(1,nrr2)/100,10)
               mbp2=kl(2,nrr2)/1000
               kbdf=kl(1,nrr2)/1000000
c     print *,'checking for dub-norm: nrr22,jc0,mbp,mbp2,kbdf= ',
c    .   nrr22,jc0,mbp,mbp2,kbdf,nrr2
               if((mbp2 .eq. mbp) .and. (kbdf .eq. kbdf2) .and.
     .               (jc0 .eq. 0)) then
                  call bdfill(kl,xl,nrr2,kbdf,g2bd,t2bd,e2bd,trm2bd,
     .               trp2bd)
c: treat two normal rays as a caustic doublet pair:
                  jc(kbdf)=2
                  jc2=2
c     print *,'detected shadow or dublet to two normal transition'
               endif
            endif
         endif
108      jc1=jc2
110   continue
112   continue
c     do 200 k=kbdf1,kbdf2
c        write(60,300) bdf(k),log10(bdf(k)),tbd(k),trpbd(k),gbd(k),
c    .      trmbd(k),ebd(k)
300      format(f8.2,1x,f6.3,1x,5(e15.8,1x))
200   continue
c     do 202 k=kbdf1,kbdf2
c        write(60,300) bdf(k),log10(bdf(k)),t2bd(k),trp2bd(k),g2bd(k),
c    .      trm2bd(k),e2bd(k)
202   continue
c
      return
      end
