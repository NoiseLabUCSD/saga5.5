      subroutine layers(rlay,drlay,ibdlay,rbdlay,drbdlay)
c
c: this subroutine calculates the r(a) and r'(a) functions at the points
c: in the array aaxis and puts them in the arrays rlay and drlay.
c
      implicit integer*4(i-n)
      include 'common/bottom' 
      include 'common/laydata'
      include 'common/discpt' 
      include 'common/bdcom' 
      real rlay(nlayer,inx),drlay(nlayer,inx),rbdlay(nlayer,inx),
     .   drbdlay(nlayer,inx)
      integer*4 ibdlay(nlayer,inx)
c
      do 20 j=1,ntot
        do 30 k=1,int(kdisc(jpen(j))/10) - 1
            a=aaxis(k)
            ainv=1./a
            call rdrlay(kprof(j),a,ainv,cp1(j),cp2(j),z(j),bp(j),
     .            bet(j),0.,rlay(j,k),drlay(j,k),kturn)
            if((kturn .eq. -1) .and. (ainv .lt. cp1(j+1))) then
               call delta(2,a,rho2(j),cp2(j),cs2(j),rho1(j+1),
     .            cp1(j+1),cs1(j+1),rbdlay(j,k),drbdlay(j,k),dr2b,
     .            ibdlay(j,k))
               rlay(j,k)=rlay(j,k) + rbdlay(j,k)/2.
               drlay(j,k)=drlay(j,k) + drbdlay(j,k)/2.
            endif
30       continue
20    continue
c
      return
      end 
