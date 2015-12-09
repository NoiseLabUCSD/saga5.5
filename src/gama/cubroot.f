      subroutine cubroot(rr,x1,x2,xhit,kcub)
c
c: this subroutine finds the root, xhit, between x1 and x2 of the
c: polynomial f(x)=a*x**3 + b*x**2 + c*x + d - rr, where the coeffi-
c: cients a,b,c,d are found in the common block poly.
c
      implicit integer*4(i-n)
      include 'common/poly'
      include 'common/pii'
c
      a1=b/a
      a2=c/a
      a3=(d-rr)/a
      q=(3.*a2 - a1**2)/9.
      r=(9.*a1*a2 - 27.*a3 - 2.*a1**3)/54.
      dd=q**3 + r**2
      if(dd .ge. 0.) then
         rdd=sqrt(dd)
         xhit=sign(1.,r+rdd)*abs(r+rdd)**(.33333) +
     .      sign(1.,r-rdd)*abs(r-rdd)**(.33333) - a1/3.
      else
         rq=sqrt(-1.*q)
         arg=r/rq**3
         if((arg .lt. -1.d0) .or. (arg .gt. 1.d0)) then
            xhit=.5*(x1+x2)
            kcub=0
            return
         else
            th=acos(arg)/3.
         endif
         xhit=2.*rq*cos(th) - a1/3.
         if((xhit .lt. 0.d0) .or. (xhit .gt. 1.d0)) then
            xhit=2.*rq*cos(th + pietth) - a1/3.
            if((xhit .lt. 0.d0) .or. (xhit .gt. 1.d0)) then 
               xhit=2.*rq*cos(th + piefth) - a1/3.
            endif
         endif
      endif
      if((xhit .le. 0.d0) .or. (xhit .ge. 1.d0)) then
         xhit=0.5d0
         kcub=0
      endif
      xhit=x1 + xhit*(x2-x1)
c     if(kcub .eq. 0) print *,'kcub=0 in cubroot' 
c
      return
      end 
