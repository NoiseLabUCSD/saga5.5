      subroutine proflay(j,k,jok,bkold,iich)
c
c: this subroutine fits either a linear segment or a curved segment
c: to the data points given in locations j and k in the ocean svp.  if
c: the curved profile cannot be fit to the points and slopes at the
c: two points, then jok is set to zero. 
c
      implicit integer*4(i-n)
      include 'common/svp'
c
      h=zsvp(k)-zsvp(j)
      jok=1
      if(kseg .eq. 1) then
c: linear segmentation:
         g0(k,2)=(csvp(k)-csvp(j))/h
         g0(k,1)=-1.*g0(k,2)
         return
      endif
c: curved segmentation:  
      iich=0
      bkold=bsvp(1,k)
      call bsvpfit(zsvp(j),zsvp(k),csvp(j),csvp(k),bsvp(2,j),btest)
      bsrat=bsvp(1,k)/btest
      if(bsvp(1,j)*bsvp(2,j) .lt. 0.) then
         bbfac=bsvp(2,j)*btest
         bfit=bbfac/bsvp(1,k)
         blin=(csvp(j+1)-csvp(j))/(zsvp(j+1)-zsvp(j))
         if(bfit/blin .lt. .5) goto 99
c     print *,'b at ext: ',j,zsvp(j),bfit
         bsvp(2,j)=bfit
      elseif(bsvp(1,k)*bsvp(2,k) .lt. 0.) then
         blin=(csvp(k-1)-csvp(k))/(zsvp(k-1)-zsvp(k))
         if(btest/blin .lt. .5) goto 99
c     print *,'b at ext: ',k,zsvp(k),btest
         bsvp(1,k)=btest
      elseif(bsrat .gt. .75 .and. bsrat .lt. 1.25) then
         bsvp(1,k)=btest
         bsvp(2,k)=btest
         iich=1
      else
         goto 99
      endif
      dz=zsvp(k) - zsvp(j)
      dc2ok=abs(5.*(bsvp(1,k) - bsvp(2,j))/h)
      call winefit(k,csvp(j),csvp(k),bsvp(2,j),bsvp(1,k),h,dc2max,jok)
      if(dc2max .gt. dc2ok) then
c        print *,'c" big: k,dc2max,dc2ok = ',j,k,zsvp(j),zsvp(k),
c    .      dc2max/dc2ok
         jok=0
         if(iich .eq. 1) then
            bsvp(1,k)=bkold
            bsvp(2,k)=bkold
         endif
      endif
c
      return
99    jok=0
      return
      end 
