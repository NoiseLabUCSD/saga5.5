      subroutine bdparsh(nr,kr,nrp,xlist,klist,nlist,nda)
c
c: this subroutine is called from shadcal when a beam displacement
c: caustic is found (iibdc=1).  a partially reflected ray is searched
c: for, and, if found corrected using an ad hoc formula.
c
      implicit integer*4(i-n)
      include 'common/srlox'
      include 'common/pii'
      include 'common/pathway'
      include 'common/discint'
      include 'common/bdcom'
      integer*4 klist(2,mray,nrtot),nlist(2,nrtot),ktemp(2,20)
      real xlist(7,mray,nrtot),xtemp(7,20)
c
c: try to find previous partially reflected eigenray at the 
c: same range and frequency:
12    nrp=0
      do 10 nrr=nr,1,-1
         mbp2=klist(2,nrr,kr)/1000
c: stop looking if we've gone to a different path:
         if(mbp2 .ne. mbp) goto 40
         apar=xlist(1,nrr,kr)
         jcpar=mod(klist(1,nrr,kr)/100,10)
c: check if a normal ray in correct discontinuity interval:
         if((jcpar .eq. 0) .and. (apar .lt. adhip(1)) .and.
     .      (apar .gt. adlop(1))) then
            ibd2=mod(klist(1,nrr,kr)/1000,10)
c: if no bd in ray, make copies of it so that f-dep can be included:
            if(ibd2 .eq. 0) then
c: store rays after nrr in temporary array:
               nrtemp=0
               do 15 nrc=nrr+1,nr
                  nrtemp=nrtemp + 1
                  do 940 jcan=1,7
                     xtemp(jcan,nrtemp)=xlist(jcan,nrc,kr)
940               continue
cmay              xtemp(1:7,nrtemp)=xlist(1:7,nrc,kr)
                  ktemp(1,nrtemp)=klist(1,nrc,kr)
                  ktemp(2,nrtemp)=klist(2,nrc,kr)
15             continue
               nlist(1,kr)=nlist(1,kr) - nrtemp
c: put in nbdf rays with ibd=1 and kbdf=1,nbdf:
               jbdf=klist(1,nrr,kr)/1000000
               klist(1,nrr,kr)=klist(1,nrr,kr)+1000+(jbdf-1)*1000000
               kcode1=klist(1,nrr,kr)
               kcode2=klist(2,nrr,kr)
               do 30 nrc=1,nbdf-1
                  if(nlist(1,kr) .ge. mray) then
                     call dawrite(kr,xlist,klist,nlist,nda)
                  endif
                  nlist(1,kr)=nlist(1,kr) + 1
                  nra=nlist(1,kr)
                  do 942 jcan=1,7
                     xlist(jcan,nra,kr)=xlist(jcan,nra-1,kr)
942               continue
cmay              xlist(1:7,nra,kr)=xlist(1:7,nra-1,kr)
                  klist(1,nra,kr)=kcode1 + nrc*1000000
                  klist(2,nra,kr)=kcode2
c                 jbdf=klist(1,nra,kr)/1000000
c                 ibd2=mod(klist(1,nra,kr)/1000,10)
c                 jc=mod(klist(1,nra,kr)/100,10)
c       print *,'partials added: nra,jbdf,ibd2,jc = ',nra,jbdf,ibd2,jc
30             continue
c: put the nrtemp eigenrays back into eigenray list:
               do 35 nrc=1,nrtemp
                  if(nlist(1,kr) .ge. mray) then
                     call dawrite(kr,xlist,klist,nlist,nda)
                  endif
                  nlist(1,kr)=nlist(1,kr) + 1
                  do 944 jcan=1,7
                     xlist(jcan,nlist(1,kr),kr)=xtemp(jcan,nrc)
944               continue
cmay              xlist(1:7,nlist(1,kr),kr)=xtemp(1:7,nrc)
                  klist(1,nlist(1,kr),kr)=ktemp(1,nrc)
                  klist(2,nlist(1,kr),kr)=ktemp(2,nrc)
35             continue
c: reset nr and go back to look for ray at correct frequency kbdf:
               nr=nlist(1,kr)
               goto 12
            endif
c: check if partial ray found at correct frequency:
            kbdf2=klist(1,nrr,kr)/1000000
            if(kbdf2 .eq. kbdf) then
               nrp=nrr
               goto 40
            endif
         endif
10    continue
40    continue
c
      return
      end 
