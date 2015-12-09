      subroutine newsdep(nlist,plcoh,plinc,nda,tmin,dbmax,range,NDAW)
c
      implicit integer*4(i-n)
      include 'common/freqcom'
      include 'common/depth'
      include 'common/srlox'
      include 'common/pathway'
      real*4 plinc(nrtot,nfr2)
      complex*8 plcoh(nrtot,nfr2)
      integer*4 nlist(2,nrtot),nda(nrtot,0:50) 
      real tmin(nrtot),dbmax(nrtot)
      real range(nrtot)
c
c: nnda is the direct access file record counter. nda tells the record
c: numbers written for each receiver.
      nnda=0
      do 5 kr=1,nrtot
         nda(kr,0)=0
5     continue
c: open file in case eigenrays overflow from xlist,klist arrays: 
c: record length is 4*(# integer vars) + 8*(# real vars): 
      lrec=NDAW*((1+2*mray) + 8*(7*mray))
      open(54,access='direct',recl=lrec,form='unformatted',
     .   status='scratch')
c     open(54,file='eigtemp',access='direct',recl=lrec,
c    .   form='unformatted')
      zs=srloc(jzs)
c: BUG discovered 4-2-90: Initialize mbp=0 in newsdep only:
      mbp=0
c
c: initialize variables for new source depth
      do 10 kr=1,nrtot
         nlist(1,kr)=0
         tmin(kr)=1.e13
         dbmax(kr)=0.
         do 20 kfr=1,nfr2
            plcoh(kr,kfr)=cmplx(0.,0.)
            plinc(kr,kfr)=0.
20       continue
10    continue
c
      return
      end 
