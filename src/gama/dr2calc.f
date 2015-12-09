      subroutine dr2calc(a,r,dr,dr2,ibd,rbd,drbd,dr2bd)
c
c: this subroutine calculates r, r', and r" for the ray with
c: the snell inv "a" and the path given by icpth, nttot, and nbot.
c: see rdrcalc for additional comments. 
c
      implicit integer*4(i-n)
      include 'common/svp'
      include 'common/bottom' 
      include 'common/paths'
      include 'common/pathway'
      include 'common/attcon' 
      include 'common/gamaoptions'
c
      ainv=1./a
c: add up r and r' in all of the ocean layers in each of the
c: three ocean sections.
      r=0.
      dr=0.
      dr2=0.
      ibd=0
c     rbd=0.
c     drbd=0.
      dr2bd=0.
      do 10 j=1,3
         if(ncpth(icpth,j) .ne. 0) then 
            rleg=0. 
            drleg=0.
            dr2leg=0.
            iup=-1
            if(j .eq. 1) iup=1
            jup=(3-iup)/2
            do 20 k=kseq(jup,j,1),kseq(jup,j,2),kseq(jup,j,3)
               kk=k+jup-1
               call dr2lay(kseg,a,ainv,csvp(k),
     .            csvp(k-iup),zsvp(kk)-zsvp(kk-1),g0(kk,jup),
     .            g1(kk,jup),g2(kk,jup),rr,drdr,dr2dr2,kturn)
               rleg=rleg + rr 
               drleg=drleg + drdr
               dr2leg=dr2leg + dr2dr2
               if(kturn .ne. -1) goto 15
20          continue
            if(((j .eq. 3) .and. (ainv .lt. cp1(1))) .or.
     .         ((j .eq. 1) .and. (ainv .lt. cp1(-1)))) then 
               jj=j-3
               call bd2calc(a,rho2(jj),cp2(jj),cs2(jj),rho1(jj+1),
     .            cp1(jj+1),cs1(jj+1),rb,drb,drb2,ibd)
               rleg=rleg + rb/2.
               drleg=drleg + drb/2.
               dr2leg=dr2leg + drb2/2.
c              rbd=rbd + rb*ncpth(icpth,j)/2
c              drbd=drbd + drb*ncpth(icpth,j)/2
c              dr2bd=dr2bd + drb2*ncpth(icpth,j)/2
            endif
15          r=r + rleg*ncpth(icpth,j)
            dr=dr + drleg*ncpth(icpth,j) 
            dr2=dr2 + dr2leg*ncpth(icpth,j)
         endif
10    continue
c
c: add up the r-a and r'-a diagrams in the bottom layers
      do 30 j=1,ndp 
         call dr2lay(kprof(j),a,ainv,cp1(j),cp2(j),z(j),bp(j),
     .      bet(j),0.,rr,drdr,dr2dr2,kturn)
         if(j .eq. ndp) then
         if((kturn .eq. -1) .and. (ainv .le. cp1(j+1))) then
            call bd2calc(a,rho2(j),cp2(j),cs2(j),rho1(j+1),cp1(j+1),
     .         cs1(j+1),rb,drb,drb2,ibd)
            rr=rr + rb/2.
            drdr=drdr + drb/2.
            dr2dr2=dr2dr2 + drb2/2.
c           rbd=rbd + rb*nttot(j)/2.
c           drbd=drbd + drb*nttot(j)/2. 
c           dr2bd=dr2bd + drb2*nttot(j)/2.
         endif
         endif
         r=r + rr*nttot(j)
         dr=dr + drdr*nttot(j)
         dr2=dr2 + dr2dr2*nttot(j)
c: if the ray turns or is totally reflected, the summation stops.
         if(kturn .ne. -1) goto 50
30    continue
c
50    return
      end 
