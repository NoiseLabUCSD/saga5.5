      subroutine mem_lim(n,nlim,eline,le,vname,lv,iibad)
c
      implicit none
c
c      implicit integer*4(i-n)
      integer*4 n,nlim,le,lv,iibad
      character(len=*) eline,vname
c
      if(n .gt. nlim) then
         print *,' '
         print *,eline(1:le)
         print *,'VARIABLE NAME = ',vname(1:lv),'; LIMIT = ',nlim
         print *,'ENTERED OR COMPUTED VALUE FOR THIS RUN = ',n
         iibad=1
      endif
c
      return
      end
