      subroutine laycrit(nbd,ntrav,npen,jfall)
c
c: this subroutine changes jpen and jtrav so that kdisc(jpen(k))
c: is the point in aaxis for which the ray barely penetrates the
c: kth layer.  likewise, jtrav is for the ray which barely traverses
c: the kth layer.
c
      implicit integer*4(i-n)
      include 'common/discpt' 
      include 'common/inivar' 
c
5     if(nbd .ge. 0) then
         if(jbd(nbd) .eq. jfall) then
            jbd(nbd)=kpt
            nbd=nbd-1
            goto 5
         endif
      endif
      if(jbset .eq. 0) then
         if(jbd(-2) .eq. jfall) then
            jbd(-2)=kpt
            jbset=1 
         endif
      endif
c
10    if(ntrav .ge. 0) then
         if(jtrav(ntrav) .eq. jfall) then
            jtrav(ntrav)=kpt
            ntrav=ntrav-1
            goto 10 
         endif
      endif
      if((jtset .eq. 0) .and. (jtrav(-2) .eq. jfall)) then
         jtrav(-2)=kpt
         jtset=1
      endif
c
20    if(npen .ge. 0) then
         if(jpen(npen) .eq. jfall) then 
            jpen(npen)=kpt
            npen=npen-1
            goto 20 
         endif
      endif
      if((jpset .eq. 0) .and. (jpen(-2) .eq. jfall)) then
         jpen(-2)=kpt
         jpset=1
      endif
c
      return
      end 
