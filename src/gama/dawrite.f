      subroutine dawrite(kr,xlist,klist,nlist,nda)
c
      implicit integer*4(i-n)
      include 'common/srlox'
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot),nda(nrtot,0:50) 
      real xlist(7,mray,nrtot)
c
c: write out all eigenrays found up to last disc point to file    
c: [leave others there for dublet check and weakness check in rcheck]
      nogo=0
      nray=max0(2,nlist(2,kr))
      nleft=nlist(1,kr) - nray
      if(nda(kr,0) .ge. 50) then
         print *,'more than 2500 rays found for receiver kr = ',
     .      kr,'; some rays ignored.'
         nda(kr,0)=nda(kr,0) - 1
      endif
c: NEW 4-13-90: nnda is counter for direct access writes. nda(kr,0)
c: holds the # of writes for kr'th receiver. nda(kr,j) is the record
c: number for the j'th direct access write.
      nnda=nnda + 1
      write(54,rec=nnda,err=500) nray,((xlist(jcan,jcan2,kr),
     .   jcan=1,7),jcan2=1,nray),((klist(jcan,jcan2,kr),
     .   jcan2=1,nray),jcan=1,2)
      nda(kr,0)=nda(kr,0) + 1
      nda(kr,nda(kr,0))=nnda
      do 940 jcan=1,nleft
         do 942 jcan2=1,7
            xlist(jcan2,jcan,kr)=xlist(jcan2,nray+jcan,kr)
942      continue
         klist(1,jcan,kr)=klist(1,nray+jcan,kr)
         klist(2,jcan,kr)=klist(2,nray+jcan,kr)
940   continue
cmay  xlist(1:7,1:nleft,kr)=xlist(1:7,nray+1:nlist(1,kr),kr)
cmay  klist(1:2,1:nleft,kr)=klist(1:2,nray+1:nlist(1,kr),kr)
c     print *,'kr,nray,nleft = ',kr,nray,nleft
      nlist(1,kr)=nleft
      nlist(2,kr)=0
c
      return
c
500   print *,'error writing to temp direct access file for eig char.'
      print *,'probably ran out of space on /tmp.  some eigs ignored.'
      nnda=nnda - 1
c
      return
      end 
