      subroutine botlos(range,plcoh,jjnr,ang,blph,blmag)
c
c: this subroutine computes the bottom loss for the ranges r and
c: frequencies f according to the complex field values in plcoh.
c: a bottom loss table written out to blt file.
c
      implicit integer*4(i-n)
      include 'common/svp'
      include 'common/depth'
      include 'common/pii'
      include 'common/charcom'
      include 'common/srlox'
      include 'common/freqcom'
      complex*8 plcoh(nrtot,nfr2)
      complex cbl,eikr1,eyem
cmay  real*4 blmag(:),blph(:,:),ang(:)
      real*4 blmag(nfr),blph(nfr,nang),ang(nang)
      real range(nrtot)
      integer*4 jjnr(nrtot)
      data eyem/(0.,-1.)/
c
cmay  allocate(ang(nang),blmag(nfr),blph(nfr,nang))
      if(nang .gt. 1001 .or. nfr .gt. 1001) then
         print *,'Increase array sizes in botlos: nang,nfr = ',
     .      nang,nfr
         stop
      endif
      h2=2.*zsvp(nsvp) - zs - zr
      h2sq=h2**2
      do 10 kr=1,nang
         mr=jjnr(mmr1 + kr - 1)
         ang(kr)=atan2(h2,range(mr))*piedeg
10    continue
      open(14,file=outn(1:loutn)//'.blt',status='unknown',
     .     form='formatted')
      write(14,100)
100   format('*(1) # angles, na; # frequencies, nf')
      write(14,110) nang,nfr
110   format(13x,i4,15x,i4)
      write(14,120)
120   format('*(2) list of frequencies')
      write(14,130) (frq(jcan),jcan=1,nfr)
130   format(107(f9.2,1x))
      write(14,140)
140   format('*(3) angle  dbloss,phase  dbloss,phase...[nf pairs]',
     .   ' [na lines like this]')
      do 15 kr=1,nang
         mr=jjnr(mmr1 + kr - 1)
         r1=sqrt(h2sq + range(mr)**2)
         do 20 kfr=1,nfr
            xk=twpie*frq(kfr)/csvp(nsvp)
            eikr1=exp(eyem*xk*r1)
            cbl=plcoh(mr,kfr)*eikr1*r1
            if(cbl .ne. cmplx(0.,0.)) then
               blmag(kfr)=-20.*log10(abs(cbl))
               blph(kfr,kr)=atan2(aimag(cbl),real(cbl))*piedeg
            else
               blmag(kfr)=99.99 
               blph(kfr,kr)=0.
               print *,'warning: bottom loss infinite for kr,kf = ',
     .            kr,kfr,'; f,ang = ',frq(kfr),ang(kr)
            endif
20       continue
         write(14,150) ang(kr),ctab,(blmag(kfr),ctab,blph(kfr,kr),
     .      kfr=1,nfr)
150      format(f6.3,a1,107(f6.2,a1,f6.1,2x))
15    continue
c
c: check for phase of bottom loss varying too much:
      do 40 kfr=1,nfr
         do 50 kr=1,nang-1
            phdif=blph(kfr,kr+1) - blph(kfr,kr)
            if(abs(phdif) .gt. 180.) phdif=phdif - sign(360.,phdif)
            if(abs(phdif) .gt. 180.) then
               write(8,200) frq(kfr),ang(kr),ang(kr+1)
200   format('warning: phase of bottom loss changed by > 90 deg for '/
     .   '   freq = ',f9.2,'  angles = ',f5.2,2x,f5.2)
            endif
50       continue
         if(kfr .eq. nfr) goto 40
         do 60 kr=1,nang
            phdif=blph(kfr+1,kr) - blph(kfr,kr)
            if(abs(phdif) .gt. 180.) phdif=phdif - sign(360.,phdif)
            if(abs(phdif) .gt. 180.) then
               write(8,210) ang(kr),frq(kfr),frq(kfr+1)   
210   format('warning: phase of bottom loss changed by > 90 deg for '/
     .   '   angle = ',f5.2,'  frequencies = ',f9.2,2x,f9.2)
            endif
60       continue
40    continue
c
      return
      end
