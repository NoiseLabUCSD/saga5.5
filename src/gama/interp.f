      subroutine interp(anext,npt,kint,kax,cs)
c
c: this subroutine enters values in aaxis from the previous disc
c: point to the next one, anext.  points in between are filled in
c: at appropriate intervals found from the array space.
c
      implicit integer*4(i-n)
      include 'common/laydata'
      include 'common/inivar' 
      integer*4 npt(3)
c
      theps1=theps
      thnext=asin(min(1.,anext*cs))
      if(kax .eq. 1) then
c: do not fill in any points from a=0 to the first disc pt. 
         kax=kax+2
      else
c: place sample points at appropriatly spaced theta values: 
         thlast=asin(min(1.,aaxis(kax)*cs))
         thdelt=thnext-thlast 
         do 10 k=1,npt(knext)-1
            aaxis(kax+k)=sin(thlast + space(knext,k)*thdelt)/cs
10       continue
         kax=kax + npt(knext) + 1
      endif
c
c: fill in the two points which bracket the next disc pt.
20    aaxis(kax-1)=sin(thnext-theps1)/cs 
cxx   print *,anext-aaxis(kax-1),1.d0/anext
      if(abs(anext-aaxis(kax-1)) .le. 1.e-14) then
c     if(abs(anext-aaxis(kax-1)) .le. 1.e-15) then
         theps1=theps1*10.
         goto 20
      endif
      if(kax .le. 150) aaxis(kax)=sin(thnext+theps1)/cs
c
c: npt(knext) is the number of sample points to put between the current
c: and the next discontinuity point.
      knext=kint
c
      return
      end 
