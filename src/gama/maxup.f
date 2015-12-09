      subroutine maxup(xlist,klist,nlist,tmin,dbmax)
c
c: this subroutine updates tmin and dbmax for the rays found for the
c: current ray path.
c
      implicit integer*4(i-n)
      include 'common/caustix'
      include 'common/srlox'
      include 'common/pii'
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot)
      real xlist(7,mray,nrtot),tmin(nrtot),dbmax(nrtot)
c
      do 10 kr=mmr0,mmr2
         jcaus=0
         do 20 nr=nlist(2,kr)+1,nlist(1,kr)
            if(jcaus .gt. 0) then
               jcaus=0 
               goto 20 
            endif
            jcaus=mod(klist(1,nr,kr)/100,10)
c: jcaus=0 indicates a normal eigenray: 
            if(jcaus .eq. 0) then
               tmin(kr)=min(tmin(kr),xlist(3,nr,kr))
               dbmax(kr)=max(dbmax(kr),xlist(6,nr,kr))
c: jcaus=2 indicates a pair of caustic eigenrays (insonified zone): 
            elseif(jcaus .eq. 2) then
               nr2=nr+1
               dbmax(kr)=max(dbmax(kr),xlist(6,nr2,kr))
               tmin(kr)=min(tmin(kr),xlist(3,nr,kr),xlist(3,nr2,kr))
            endif
20       continue
10    continue
c
      return
      end 
