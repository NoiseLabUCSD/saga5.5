      subroutine linefit(x,kr1,kr2,kr3,xlist,klist,nlist,
     .   wk1,wk2,wk3,nda)
c
c: this subroutine linearly fits the eigenray characteristics.
c
      implicit integer*4(i-n)
      include 'common/pii'
      include 'common/srlox'
      logical wk1,wk2,wk3
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot),nda(nrtot,0:50) 
      real xlist(7,mray,nrtot),x(nrtot) 
c
      nr1=nlist(1,kr1)
      nr2=nlist(1,kr2)
      nr3=nlist(1,kr3)
      if((kr2-kr1 .gt. 1) .and. (.not. wk1 .or. .not. wk2)) then 
         do 10 kr=kr1+1,kr2-1 
            if(nlist(1,kr) .ge. mray) then
               call dawrite(kr,xlist,klist,nlist,nda)
            endif
            nlist(1,kr)=nlist(1,kr) + 1 
            klist(1,nlist(1,kr),kr)=klist(1,nr2,kr2)
            klist(2,nlist(1,kr),kr)=klist(2,nr2,kr2)
10       continue
         dx12=x(kr2) - x(kr1)
         do 20 nc=1,7
            afac=(xlist(nc,nr2,kr2) - xlist(nc,nr1,kr1))/dx12
            bfac=xlist(nc,nr1,kr1) - afac*x(kr1)
            do 15 kr=kr1+1,kr2-1 
               xlist(nc,nlist(1,kr),kr)=afac*x(kr) + bfac
15          continue
20       continue
      endif
c
      if((kr3-kr2 .gt. 1) .and. (.not. wk2 .or. .not. wk3)) then 
         do 30 kr=kr2+1,kr3-1 
            if(nlist(1,kr) .ge. mray) then
               call dawrite(kr,xlist,klist,nlist,nda)
            endif
            nlist(1,kr)=nlist(1,kr) + 1 
            klist(1,nlist(1,kr),kr)=klist(1,nr2,kr2)
            klist(2,nlist(1,kr),kr)=klist(2,nr2,kr2)
30       continue
         dx12=x(kr3) - x(kr2)
         do 40 nc=1,7
            afac=(xlist(nc,nr3,kr3) - xlist(nc,nr2,kr2))/dx12
            bfac=xlist(nc,nr2,kr2) - afac*x(kr2)
            do 35 kr=kr2+1,kr3-1 
               xlist(nc,nlist(1,kr),kr)=afac*x(kr) + bfac
35          continue
40       continue
      endif
c     if(kr1 .lt. 78 .and. kr3 .gt. 78) then
c        print *,'linefit: kr1,kr2,kr3,wk1,wk2,wk3,nr1,nr2,nr3,mray = ',
c    .      kr1,kr2,kr3,wk1,wk2,wk3,nr1,nr2,nr3,mray
c        do 300 kr=kr1,kr3
c           print *,'kr,klist(1),(2)= ',kr,klist(1,nlist(1,kr),kr),
c    .   klist(2,nlist(1,kr),kr)
300      continue
c     endif
c
      return
      end 
