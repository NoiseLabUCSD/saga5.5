      subroutine srlev(zr,kr) 
c
c: this subroutine inserts the source or receiver depth, zr, in the
c: ocean svp between levels kr-1 and kr.
c
      implicit integer*4(i-n)
      include 'common/svp'
c
      call morelay(kr)
      zsvp(kr)=zr
      kr0=kr - 1
      kr1=kr + 1
      if(kseg .eq. 1) then
         csvp(kr)=csvp(kr0) + g0(kr,2)*(zr - zsvp(kr0))
         bsvp(1,kr)=bsvp(1,kr0)
         bsvp(2,kr)=bsvp(2,kr0)
         call proflay(kr,kr1,jok,bkold,iich)
      elseif(kseg .eq. 3) then
c: get slopes at end points. If zero, set as would be done in bsvpfit:
         h1=zsvp(kr) - zsvp(kr0)
         call cdcwine(csvp(kr0),g0(kr1,2),g1(kr1,2),g2(kr1,2),h1,
     .      csvp(kr),bsvp(1,kr),dc2,1)
         bsvp(2,kr)=bsvp(1,kr)
         call winefit(kr,csvp(kr0),csvp(kr),bsvp(2,kr0),bsvp(1,kr),h1,
     .      dc2max,jok)
         if(jok .eq. 0) then
            print *,'bad winefit call from srlev: ',1./zz
            stop
         endif
         h1=zsvp(kr1) - zsvp(kr)
         call winefit(kr1,csvp(kr),csvp(kr1),bsvp(2,kr),bsvp(1,kr1),h1,
     .      dc2max,jok)
         if(jok .eq. 0) then
            print *,'bad winefit call from srlev: ',1./zz
            stop
         endif
      endif
c
      return
      end 
