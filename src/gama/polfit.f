      subroutine polfit(x1,x2,fx1,fx2,fpx1,fpx2)
c
c: this subroutine fits a cubic polynomial "f(x)= a*x**3 + b*x**2 + c*x
c: + d" to the points x1 and x2, their functional values fx1
c: and fx2, and their derivatives fpx1 and fpx2.
c
      include 'common/poly'
      real*8 x1,x2,fx1,fx2,fpx1,fpx2
c
      d=fx1
      c=fpx1
      b=3.d0*(fx2 - c - d) - (fpx2 - c)
      a=fx2 - c - d - b
c
      return
      end 
