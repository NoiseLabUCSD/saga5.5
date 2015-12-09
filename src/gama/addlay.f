      subroutine addlay(k)
c
c: this subroutine adds a curvilinear layer between two layers for which
c: a curve could not be fit.
c
      implicit integer*4(i-n)
      include 'common/svp'
c
      k1=k-1
      k3=k+1
      dztot=zsvp(k) - zsvp(k1)
      dctot=csvp(k) - csvp(k1)
c: make space for new layer:
      call morelay(k)
c: compute c for which slope criterion can be met on both sides:
      iiok=1
      cbfac=sqrt(csvp(k1)*bsvp(2,k1)/(csvp(k3)*bsvp(1,k3)))
      zsvp(k)=.5*(zsvp(k1) + zsvp(k3))
      csvp(k)=sqrt((cbfac*csvp(k3)**2 + csvp(k1)**2)/(1. + cbfac))
      iiok=1
      if((csvp(k)-csvp(k1))*(csvp(k)-csvp(k3)) .ge. 0.) then
         cbfac=-1.*cbfac
      print *,'trying -cbfac: '
         csvp(k)=sqrt((cbfac*csvp(k3)**2 + csvp(k1)**2)/(1. + cbfac))
         if((csvp(k)-csvp(k1))*(csvp(k)-csvp(k3)) .ge. 0.) then
            csvp(k)=.5*(csvp(k1) + csvp(k3))
      print *,'both failed in addlay: ',k1,k,k3
            iiok=0
         endif
      endif
c: find slope at the new point from both sides:
      call bsvpfit(zsvp(k1),zsvp(k),csvp(k1),csvp(k),
     .   bsvp(2,k1),bsvp(1,k))
      bsvp(2,k)=bsvp(1,k)
      call bsvpfit(zsvp(k3),zsvp(k),csvp(k3),csvp(k),bsvp(1,k3),btest)
c: take shallower of two slopes if not same:
      if(iiok .eq. 0) then
         bsvp(1,k)=sign(1.,btest)*min(abs(bsvp(1,k)),abs(btest))
         bsvp(2,k)=bsvp(1,k)
      endif
c: fit to weinberg profile:
      h2=.5*dztot
      call winefit(k,csvp(k1),csvp(k),bsvp(2,k1),bsvp(1,k),
     .   h2,dc2max,jok)
      if(jok .eq. 0) then
         print *,'bad winefit call from addlay: ',1./zz
         stop
      endif
      call winefit(k3,csvp(k),csvp(k3),bsvp(2,k),bsvp(1,k3),
     .   h2,dc2max,jok)
      if(jok .eq. 0) then
         print *,'bad winefit call from addlay: ',1./zz
         stop
      endif
c
      return
      end 
