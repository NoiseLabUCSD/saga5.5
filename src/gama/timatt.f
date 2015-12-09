      subroutine timatt(a,exptot,travtot,ibd)
c
c: this subroutine calculates the travel time travtot and the total
c: frequency-dependent attenuation factor (in db/khz) exptot for the ray
c: with snell invariant "a" and the path specified by variables in
c: pathway.
c
      implicit integer*4(i-n)
      include 'common/svp'
      include 'common/bottom' 
      include 'common/attcon' 
      include 'common/paths'
      include 'common/depth'
      include 'common/pathway'
      include 'common/gamaoptions'
c
      ibd=0
      travtot=0.
      exptot=0.
      ainv=1./a
c
c: add up the components from each of the ocean layers in each of the 
c: three ocean sections.  see also comments for subroutine rdrcalc.
c
      do 10 j=1,3
         if(ncpth(icpth,j) .ne. 0) then 
            tleg=0. 
            expleg=0.
            iup=-1
            if(j .eq. 1) iup=1
            jup=(3-iup)/2
c: allow concurrency on tmatlay calls inside loop:
            do 20 k=kseq(jup,j,1),kseq(jup,j,2),kseq(jup,j,3)
               kk=k+jup-1
c: the ocean is modeled to have constant frequency-dep attenuation:
               call tmatlay(kseg,a,ainv,csvp(k),csvp(k-iup),
     .              zsvp(kk)-zsvp(kk-1),g0(kk,jup),g1(kk,jup),
     .              g2(kk,jup),akp2(0),0.,expfac,travt,kturn)
               expleg=expleg + expfac
               tleg=tleg + travt
               if(kturn .ne. -1) goto 15
20          continue
            if(((j .eq. 3) .and. (ainv .le. cp1(1))) .or.
     .         ((j .eq. 1) .and. (ainv .lt. cp1(-1)))) then 
               jj=j-3
               jjp1=jj+1
               call delta(1,a,rho2(jj),cp2(jj),cs2(jj),
     .            rho1(jjp1),cp1(jjp1),cs1(jjp1),rb,drb,dr2b,ibd)
               arb=a*rb
               tleg=tleg + arb/2.
               expleg=expleg + akp1(jjp1)*arb/2
            endif
15          travtot=travtot + tleg*ncpth(icpth,j)
            exptot=exptot + expleg*ncpth(icpth,j)
         endif
10    continue
c
c: add up the components from each bottom layer.
      do 30 j=1,ndp 
         call tmatlay(kprof(j),a,ainv,cp1(j),cp2(j),z(j),bp(j),
     .      bet(j),0.,akp1(j),dakp(j),expfac,travt,kturn)
         if(kturn .eq. 0) goto 99
         if((kturn .eq. -1) .and. (ainv .le. cp1(j+1))) then
            jp1=j+1
            call delta(1,a,rho2(j),cp2(j),cs2(j),rho1(jp1), 
     .         cp1(jp1),cs1(jp1),rb,drb,dr2b,ibd) 
            arb=a*rb
            travt=travt + arb/2.
            expfac=expfac + akp1(jp1)*arb/2.
         endif
         travtot=travtot + travt*nttot(j)
         exptot=exptot + expfac*nttot(j)
30    continue
99    continue
c
      return
      end 
