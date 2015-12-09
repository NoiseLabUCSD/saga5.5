      subroutine rdrsum(kax,r,dr,rc,drc,ibdc,rbdc,drbdc,
     .                  rlay,drlay,ibdlay,rbdlay,drbdlay)
c
c: this subroutine finds r and r' at the kax'th sample point in the
c: snell invariant axis.  it starts with r and r' in the ocean
c: part of the ray path (rc and drc) and adds r and r' in the bottom
c: layers.
c
      implicit integer*4(i-n)
      include 'common/laydata'
      include 'common/pathway'
      include 'common/bdcom'
      real rlay(nlayer,inx),drlay(nlayer,inx),rbdlay(nlayer,inx),
     .   drbdlay(nlayer,inx),rbdc(inx),drbdc(inx),rc(inx),drc(inx)
      integer*4 ibdc(inx),ibdlay(nlayer,inx)
c
      r=rc(kax)
      dr=drc(kax)
      if(ibdc(kax) .eq. 1) then
         r=r + rbdc(kax)*bdfac
         dr=dr + drbdc(kax)*bdfac
      endif
      do 10 j=1,ndp 
         r=r + rlay(j,kax)*nttot(j)
         dr=dr + drlay(j,kax)*nttot(j)
10    continue
      if(ndp .gt. 0) then
         if(ibdlay(ndp,kax) .eq. 1) then
            r=r + rbdlay(ndp,kax)*bdfac
            dr=dr + drbdlay(ndp,kax)*bdfac
         endif
      endif
c
      return
      end 
