      subroutine pathfil
c
c: this subroutine calls traduc to translate the opt file data into
c: the info found in ncpth().
c
      implicit integer*4(i-n)
      include 'common/pcodes' 
      include 'common/pcode2' 
      include 'common/depth'
      character*1 cz,cz2,fd
c
c: translate ocean paths to format required in array ncpth:  
      ncp1=0
      do 41 k=1,nst 
         if(kb1(k) .eq. 0) then
            if(zs .lt. zr) then
               call traduc(ncp1+1,'d',0,0,'w')
            else
               call traduc(ncp1+1,'u',0,0,'w')
            endif
            call traduc(ncp1+2,'u',1,0,'w')
            ncp1=ncp1 + 2
         endif
         do 42 kb=max0(1,kb1(k)),kb2(k) 
            call traduc(ncp1+1,'d',kb-1,kb,cz(k)) 
c: enter ocean paths in order of increasing path length: 
            if(zs .ge. zr) then
               call traduc(ncp1+2,'d',kb,kb,cz(k))
               call traduc(ncp1+3,'u',kb,kb,cz(k))
            else
               call traduc(ncp1+2,'u',kb,kb,cz(k))
               call traduc(ncp1+3,'d',kb,kb,cz(k))
            endif
            call traduc(ncp1+4,'u',kb+1,kb,cz(k)) 
            ncp1=ncp1 + 4
42       continue
41    continue
c
      do 20 k=1,ncp0
         call traduc(ncp1+k,fd(k),ktop(k),kbot(k),cz2(k))
20    continue
c: check for duplicate ocean and bottom paths: 
      call ckdup
c
      return
      end 
