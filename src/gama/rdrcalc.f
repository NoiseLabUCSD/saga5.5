      subroutine rdrcalc(a,r,dr,ibd,rbd,drbd)
c
c: this subroutine calculates r and r' for the ray with the 
c: snell invariant "a" and the path given by icpth, nttot, and nbot.
c
      implicit integer*4(i-n)
      include 'common/svp'
      include 'common/bottom' 
      include 'common/paths'
      include 'common/pathway'
c
      ainv=1./a
c: add up r and r' in all of the ocean layers in each of the
c: three ocean sections.
      r=0.
      dr=0.
      ibd=0
      rbd=0.
      drbd=0.
      do 10 j=1,3
         if(ncpth(icpth,j) .eq. 0) goto 10
         rleg=0.
         drleg=0.
         iup=-1
         if(j .eq. 1) iup=1
         jup=(3-iup)/2
c
c: the counter for the ocean layers, k, is found from kseq: 
         do 20 k=kseq(jup,j,1),kseq(jup,j,2),kseq(jup,j,3)
            kk=k+jup-1
            call rdrlay(kseg,a,ainv,csvp(k),
     .         csvp(k-iup),zsvp(kk)-zsvp(kk-1),g0(kk,jup),
     .         g1(kk,jup),g2(kk,jup),rr,drdr,kturn)
            rleg=rleg + rr
            drleg=drleg + drdr
c: if ray turns in this layer, exit.
            if(kturn .ne. -1) goto 15
20       continue
c: check for bd at ocean surface or bottom: 
         if(((j .eq. 3) .and. (ainv .lt. cp1(1))) .or.
     .         ((j .eq. 1) .and. (ainv .lt. cp1(-1)))) then 
            jj=j-3
            call delta(2,a,rho2(jj),cp2(jj),cs2(jj),
     .         rho1(jj+1),cp1(jj+1),cs1(jj+1),rb,drb,dr2b,ibd)
            rleg=rleg + rb/2. 
            drleg=drleg + drb/2.
            rbd=rbd + rb*ncpth(icpth,j)/2
            drbd=drbd + drb*ncpth(icpth,j)/2
         endif
15       r=r + rleg*ncpth(icpth,j)
         dr=dr + drleg*ncpth(icpth,j)
10    continue
c
c: add up the r-a and r'-a diagrams in the bottom layers
c: allow concurrency on rdrlay calls inside loop:
      do 30 j=1,ndp 
         call rdrlay(kprof(j),a,ainv,cp1(j),cp2(j),z(j),bp(j),
     .      bet(j),0.,rr,drdr,kturn)
c: check for bd at bottom of deepest layer: 
         if(j .eq. ndp) then
            if((kturn .eq. -1) .and. (ainv .le. cp1(j+1))) then
               call delta(2,a,rho2(j),cp2(j),cs2(j),rho1(j+1),
     .            cp1(j+1),cs1(j+1),rb,drb,dr2b,ibd)
               rr=rr + rb/2.
               drdr=drdr + drb/2.
               rbd=rbd + rb*nttot(j)/2. 
               drbd=drbd + drb*nttot(j)/2.
            endif
         endif
         r=r + rr*nttot(j)
         dr=dr + drdr*nttot(j)
c: if the ray turns or is totally reflected, the summation stops.
         if(kturn .ne. -1) goto 50
30    continue
c
50    return
      end 
