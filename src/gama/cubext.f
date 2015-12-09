      subroutine cubext(a1,a2,amid,dr2,kbad)
c
c: this subroutine finds the extremum, amid, between a1 and a2 of the 
c: the function f(x)=a*x**3 + b*x**2 + c*x + d.  the coefficients a,b,c,
c: are given in the common block poly.
c
      implicit integer*4(i-n)
      include 'common/poly'
c
      radcand=b**2 - 3.*a*c
      if(radcand .lt. 0.) then
         kbad=1
         amid=.5d0*(a1+a2)
         return
      endif
      rad=sqrt(radcand)
      amid=(-1.*b + sign(1.,dr2)*rad)/(3.*a)
      if((amid .le. 0.d0) .or. (amid .ge. 1.d0)) then
         kbad=1
         amid=0.d0
      endif
      amid=a1 + amid*(a2-a1)
c
      return
      end 
