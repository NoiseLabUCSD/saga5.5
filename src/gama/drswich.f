      subroutine drswich(al,ah)
c
      implicit integer*4 (i-n)
      include './common/saga_gama_fail'

      common /info/ienv_info,lu_env,iiwrite
      integer ienv_info,lu_env,iiwrite
c
c: this subroutine finds the point abad between al and ah for which
c: drbad is of opposite sign as drl and drh.
      include 'common/dublcom'
      real al(3),ah(3)
c
      a1=al(1)
      a2=ah(1)
      kloop=0
10    kloop=kloop + 1
      if(kloop .gt. 50) then
         if(iiwrite.gt.0)
     .     print *,'kloop>50 in drswich.' 
         dummy=0.
         print *,'Almost building a core ... '
         iifail=1
         return
c         stop
      endif
      abad(1)=(a1 + a2)/2.
      call rdrcalc(abad(1),abad(2),abad(3),ibd,rbd,drbd)
      if(al(3)*abad(3) .le. 0.) return
      if(al(3) .gt. 0.) then
         if(abad(2) .gt. al(2)) then
            a1=abad(1) 
         else
            a2=abad(1)
         endif
      else
         if(abad(2) .lt. al(2)) then
            a1=abad(1)
         else
            a2=abad(1)
         endif
      endif
      goto 10
      end 
