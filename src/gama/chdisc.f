      subroutine chdisc(adisc,c,cbig,j,npt,nsample,kseg)
c
c: this subroutine checks if the sound speed c is larger than
c: cbig, in which case it corresponds to a discontinuity point,
c: a=1/c, and is entered in the array adisc.
c
      implicit integer*4(i-n)
      include 'common/gamaoptions'
      real adisc(76)
      integer*4 npt(3),nsample(76)
c
c     if((kseg .eq. 3) .and. (iibd .eq. 0)) return
      if(c .gt. cbig) then
         j=j+1
         if(j .gt. 76) then
            print *,'program cannot handle the number of ', 
     .      'discontinuity points in this svp.'
            stop
         endif
         adisc(j)=1.d0/c
         cbig=c
         nsample(j)=kseg
      endif
c
      return
      end 
