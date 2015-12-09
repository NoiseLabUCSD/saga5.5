      subroutine sublay(kk,nk)
c
c: this subroutine deletes nk layers from the ocean svp starting at
c: layer kk.
c
      implicit integer*4(i-n)
      include 'common/svp'
c
      do 10 k=kk+nk,nsvp
         j=k-nk
         zsvp(j)=zsvp(k)
         csvp(j)=csvp(k)
         bsvp(1,j)=bsvp(1,k)
         bsvp(2,j)=bsvp(2,k)
         g0(j,1)=g0(k,1)
         g1(j,1)=g1(k,1)
         g2(j,1)=g2(k,1)
         g0(j,2)=g0(k,2)
         g1(j,2)=g1(k,2)
         g2(j,2)=g2(k,2)
10    continue
      nsvp=nsvp-nk
c
      return
      end 
