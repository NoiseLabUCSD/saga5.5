      subroutine record(kr,a,gsl,t,e,trm,trp,raymag,nps,ibd,
     .   xlist,klist,nlist,nda)
c
      implicit integer*4(i-n)
      include 'common/srlox'
      include 'common/pathway'
      include 'common/bdcom'
      include 'common/gamaoptions'
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot),nda(nrtot,0:50) 
      real xlist(7,mray,nrtot)
c
      if(nlist(1,kr) .ge. mray) then
cxx   print *,'calling dawrite: kr,nray = ',kr,nlist(1,kr)
         call dawrite(kr,xlist,klist,nlist,nda)
      endif
10    nlist(1,kr)=nlist(1,kr) + 1
      nr=nlist(1,kr)
      xlist(1,nr,kr)=a
      xlist(2,nr,kr)=gsl
      xlist(3,nr,kr)=t
      xlist(4,nr,kr)=e
      xlist(5,nr,kr)=trm
      xlist(6,nr,kr)=raymag
      xlist(7,nr,kr)=trp
      klist(1,nr,kr)=nps + ibd*1000 + nbas*10000 + kbdf*1000000
      klist(2,nr,kr)=icpth + mbp*1000
      if(iidiag .ne. 0) print *,'record ray: kr,nr,1/a,t = ',kr,nr,
     .   1/a,t,raymag 
c
      return
      end 
