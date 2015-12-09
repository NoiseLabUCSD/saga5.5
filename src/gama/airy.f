      subroutine airy(z,ai,aip)
c
c: this subroutine computes the airy function ai and its derivative aip
c: of the real number z.  ekw added asymptotic forms 07-19-89.
c
      implicit integer*4(i-n)
      include 'common/pii'
      data cc1/.069444444444/ 
      data dd1/-.097222222222/
      data c1/.35502805388782/
      data c2/.258819403792807/
c
c: check if we can use asymptotic expansion (abs(z) > 4): 
c: this results in errors of .1% or less.
      if(abs(z) .gt. 4.) then 
c: neg z: 
         if(z .lt. 0.) then
            zqu=(-1.*z)**(.25)
            zeta=(.66667)*(-1.*z)**(1.5) 
            zetainv=1./zeta
            ai=(sqpinv/zqu)*(sin(zeta + piequ) - cos(zeta + 
     .         piequ)*cc1*zetainv)
            aip=-1.*sqpinv*zqu*(cos(zeta + piequ) + sin(zeta +
     .         piequ)*dd1*zetainv)
c: z < 0. : 
         else
            zqu=z**(.25)
            zeta=(.66667)*z**(1.5)
            zetainv=1./zeta
            expzet=exp(-1.*zeta)
            ai=(.5*sqpinv*expzet/zqu)*(1. - cc1*zetainv)
            aip=(-.5*sqpinv*expzet*zqu)*(1. - dd1*zetainv)
         endif
c: else do a series expansion: 
      else
         eps=.000001
         z2=z**2
         z3=z**3
c
         f=1.
         fp=0.
         df=f
c        kloop=0
         do 10 n=3,100,3
c           kloop=kloop + 1
            dfp=df*z2/float(n-1)
            df=df*z3/float(n*(n-1))
            f=f + df
            fp=fp + dfp
            if(abs(df) .lt. eps) goto 15
10       continue
c
15       g=z
         gp=1.
         dg=g
         do 20 n=4,100,3
c           kloop=kloop + 1
            dgp=dg*z2/float(n-1)
            dg=dg*z3/float(n*(n-1))
            g=g + dg
            gp=gp + dgp
            if(abs(dg) .lt. eps) goto 25
20       continue
c
25       ai=c1*f - c2*g
         aip=c1*fp - c2*gp
      endif
c
      return
      end 
