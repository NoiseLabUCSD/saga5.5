      subroutine setop(kst,kend,kdir,cz)
c
c: this subroutine sets up the next ocean path.
c
      implicit integer*4(i-n)
      include 'common/discpt' 
      include 'common/paths'
      include 'common/pathway'
      character*1 cz
c
c: set kdir, kst and kend according to which sections of the
c: ocean the current ocean path includes.  kst points to the
c: lowest disc pt in kdisc; kend to the highest.
      if(ntop + nbotx .eq. 0) then
         kdir=2
         kst=jbd(-1)
         kend=jpen(-1)
      elseif(nbotx .eq. 0) then
         kdir=3
         kst=jbd(-2)
         kend=jpen(-2)
      elseif(ntop .eq. 0) then
         kdir=2
         kst=jbd(0) 
         kend=jpen(0)
      else
         kdir=6
         kst=min0(jbd(-2),jbd(0))
         kend=jpen(0)
      endif
c
c: if only bottom section refracting eigenrays are desired: 
      if((cz .eq. 'r') .and. (nbotx .ne. 0)) kst=jtrav(0)
c: if only bottom penetrating eignrays desired: 
      if((cz .eq. 'p') .and. (nbotx .ne. 0)) kend=jtrav(0)
c
      return
      end 
