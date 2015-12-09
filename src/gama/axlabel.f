      subroutine axlabel(rmin,rmax,rlo,rhi,delt)
c
c: this subroutine finds round number limits (rlo,rhi) and an increment
c: size delt for an axis that includes (rmin,rmax).
c
      implicit integer*4(i-n)
      r1=rmin
      r2=rmax
      if(r2 .eq. r1) then
         r1=r1 - 1. 
         r2=r2 + 1. 
      endif
c
      rrr=r2-r1
      fac=1.
10    if(rrr .le. 100.) goto 20
      rrr=rrr/10.
      fac=fac*10.
      goto 10
20    if(rrr .gt. 10.) go to 30
      rrr=rrr*10.
      fac=fac/10.
      goto 20
30    delt=100.
      if(rrr .le. 50.) delt=50.
      if(rrr .le. 20.) delt=20.
      delt=delt*fac/10.
      rlo=float(int(r1/delt))*delt
      rhi=float(int(r2/delt))*delt
      if(rhi .lt. r2) rhi=rhi + delt
      if(rlo .gt. r1) rlo=rlo - delt
c
      return
      end 
