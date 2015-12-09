      subroutine bd2calc(a,rho1,cp1,cs1,rho2,cp2,cs2,r,dr,dr2,ibd)
c
c: this subroutine calculates r, r', and r" due to beam displacement
c: at an interface.  r" is calculated numerically.
      implicit integer*4(i-n)
      include 'common/gamaoptions'
c
      if(iibd .eq. 0) then
         ibd=0
         r=0.
         dr=0.
         dr2=0.
         return
      endif
      ibd=1
c
      if((cs1 .eq. 0.) .and. (cs2 .eq. 0.)) then
         call delta(3,a,rho1,cp1,cs1,rho2,cp2,cs2,r,dr,dr2,ibd)
         return
      endif
c
      ksm=0
      kbig=0
      kdis=0
      aeps=min(5.e-4/cp1,a - 1./cp2, 1./cp1 - a)
      if(aeps .ne. 5.e-4/cp1) kdis=1
10    call delta(2,a-aeps/2.,rho1,cp1,cs1,rho2,cp2,cs2,rlo,drlo,d,ibd)
      call delta(2,a+aeps/2.,rho1,cp1,cs1,rho2,cp2,cs2,rhi,drhi,d,ibd)
      ratio=abs(rhi-rlo)/max(abs(rhi),abs(rlo)) 
      if((ratio .gt. .1) .and. (kbig .eq. 0)) then
         ksm=1
         aeps=.2*aeps
         goto 10
      elseif((ratio .lt. .0005) .and. (ksm .eq. 0) .and. (kdis .eq. 0))
     .      then
         kbig=1
         aept=5.*aeps
         aeps=min(aept,a - 1./cp2, 1./cp1 - a)
         if(aeps .ne. aept) kdis=1
         goto 10
      endif
      dr2=(drhi-drlo)/aeps
      call delta(2,a,rho1,cp1,cs1,rho2,cp2,cs2,r,dr,d,ibd)
c
      return
      end 
