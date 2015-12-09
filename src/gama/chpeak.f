      subroutine chpeak(adisc,c,k,cbig,j,npt,nsample,kseg2) 
c
c: this subroutine checks if the sound velocity point is a local
c: maximum.  for curved profiles this generates a discontinuity point.
      implicit integer*4(i-n)
      include 'common/svp'
      real adisc(76)
      integer*4 nsample(76),npt(3)
c
      if(c .gt. cbig) then
         if((c .ge. csvp(max0(k-1,0))) .and. (c .ge. csvp(min0(
     .      k+1,nsvp)))) then 
            j=j+1
            if(j .gt. 76) then
               print *,'discontinuity point limit reached in ocean ', 
     .                 'svp: ',j
               stop 
            endif
            adisc(j)=1./c
            cbig=c
            nsample(j)=kseg2
         endif
      endif
c
      return
      end 
