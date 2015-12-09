      subroutine morelay(kk)
c
c: this subroutine makes a space in the ocean svp between layers kk-1 
c: and kk.
c
      implicit integer*4(i-n)
      include 'common/svp'
c
      if(nsvp .ge. 51) then
         print *,'tried to add too many layers to ocean svp: ',51
         stop
      endif
      do 10 k=nsvp,kk,-1
         zsvp(k+1)=zsvp(k)
         csvp(k+1)=csvp(k)
         bsvp(1,k+1)=bsvp(1,k)
         bsvp(2,k+1)=bsvp(2,k)
         g0(k+1,1)=g0(k,1)
         g1(k+1,1)=g1(k,1)
         g2(k+1,1)=g2(k,1)
         g0(k+1,2)=g0(k,2)
         g1(k+1,2)=g1(k,2)
         g2(k+1,2)=g2(k,2)
10    continue
      nsvp=nsvp+1
c
      return
      end 
