      subroutine cdcwine(c1,g0,g1,g2,z,c,dc,dc2,nd)
c
c: this subroutine calculates the sound speed and its first two
c: derivatives for the Weinberg profile.
c
      implicit integer*4(i-n)
      g2fac=1. + g2*z
      g2facsq=g2fac**2
      c=(1./c1**2 + z*(g0+g1*z)/g2facsq)**(-.5)
      if(nd .eq. 0) return
      ggg=2.*g1 - g0*g2
      fac1=(g0 + z*ggg)/(g2fac*g2facsq)
      dc=-.5*c**3*fac1
      if(nd .eq. 1) return
      dfac1=(g2fac*ggg - 3.*g2*(g0 + z*ggg))/g2facsq**2
      dc2=-.5*c**2*(3.*fac1*dc + c*dfac1)
c
      return
      end
