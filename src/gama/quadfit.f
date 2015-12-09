      subroutine quadfit(x,kr1,kr2,kr3,xlist,klist,nlist,
     .   wk1,wk2,wk3,nda)
c
c: this subroutine quadratically fits the eigenray characteristics.
c
      implicit integer*4(i-n)
      include 'common/pii'
      include 'common/srlox'
      logical wk1,wk2,wk3
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot),nda(nrtot,0:50) 
      real xlist(7,mray,nrtot),x(nrtot) 
c
c: initialize interpolation constants:  
      x1sq=x(kr1)**2
      dx12=x(kr1) - x(kr2)
      dx13=x(kr1) - x(kr3)
      dxs12=x1sq - x(kr2)**2
      dxs13=x1sq - x(kr3)**2
      aden=dxs12*dx13 - dxs13*dx12
      brat=dxs12/dx12
c: increment the number of rays found for the interpolated ranges
c: and interpolate klist: 
      if(.not. wk1 .or. .not. wk2) then 
         do 10 kr=kr1+1,kr2-1 
            if(nlist(1,kr) .ge. mray) then
               call dawrite(kr,xlist,klist,nlist,nda)
            endif
            nlist(1,kr)=nlist(1,kr) + 1 
            klist(1,nlist(1,kr),kr)=klist(1,nlist(1,kr1),kr1)
            klist(2,nlist(1,kr),kr)=klist(2,nlist(1,kr1),kr1)
10       continue
      endif
      if(.not. wk2 .or. .not. wk3) then 
         do 20 kr=kr2+1,kr3-1 
            if(nlist(1,kr) .ge. mray) then
               call dawrite(kr,xlist,klist,nlist,nda)
            endif
            nlist(1,kr)=nlist(1,kr) + 1 
            klist(1,nlist(1,kr),kr)=klist(1,nlist(1,kr2),kr2)
            klist(2,nlist(1,kr),kr)=klist(2,nlist(1,kr2),kr2)
20       continue
      endif
c
c: do interpolation for a,gsl,t,e,trm,trp:  
      nr1=nlist(1,kr1)
      nr2=nlist(1,kr2)
      nr3=nlist(1,kr3)
      do 30 nc=1,7
         xl1=xlist(nc,nr1,kr1)
         dy12=xl1 - xlist(nc,nr2,kr2)
         dy13=xl1 - xlist(nc,nr3,kr3)
         aq=(dy12*dx13 - dy13*dx12)/aden
         bq=dy12/dx12 - aq*brat
         cq=xl1 - aq*x1sq - bq*x(kr1)
      if(nc .eq. 1) print *,'aq,bq,cq,dy12,dy13 = ',aq,bq,cq,dy12,dy13
         if(.not. wk1 .or. .not. wk2) then
            do 40 kr=kr1+1,kr2-1
               xlist(nc,nlist(1,kr),kr)=aq*x(kr)**2 + bq*x(kr) + cq
40          continue
         endif
         if(.not. wk2 .or. .not. wk3) then
            do 45 kr=kr2+1,kr3-1
               xlist(nc,nlist(1,kr),kr)=aq*x(kr)**2 + bq*x(kr) + cq
45          continue
         endif
30    continue
c
c     print *,'quadph done: ',(xlist(7,nlist(1,kr),kr),kr=kr1,kr3)
c
      return
      end 
