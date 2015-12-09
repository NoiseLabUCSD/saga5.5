      subroutine ocean(r0,dr0,ibd0,rbd0,drbd0)
c
c
c: this subroutine calculates the r(a) and r'(a) functions at the
c: snell invariant points in the array aaxis for the three ocean
c: sections and places them in the arrays r0 and dr0.  the 196th,
c: 197th, and 198th ocean paths are used in calling rdrcalc to
c: isolate each of the three ocean sections: top, middle, and
c: bottom.
c
      implicit integer*4(i-n)
      include 'common/laydata'
      include 'common/pathway'
      include 'common/bdcom'
      real r0(3,inx),dr0(3,inx),rbd0(3,inx),drbd0(3,inx)
      integer*4 ibd0(3,inx)
c
c: initialize frequency at which beam displacement will be calculated:
      bdfreq=bdf(1)
      bdomega=bdom(1)
      ndp=0
      do 10 k=1,inx 
         a=aaxis(k) 
         do 20 icleg=1,3
            icpth=icleg+195
            call rdrcalc(a,r0(icleg,k),dr0(icleg,k),ibd0(icleg,k),
     .         rbd0(icleg,k),drbd0(icleg,k))
20       continue
10    continue
c
      return
      end 
