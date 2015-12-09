      subroutine layparm(cp1x,cp2x,cs1x,cs2x,rho1x,rho2x,
     .   alp1x,alp2x,als1x,als2x,fax1,nx_p2,nx_s1,nx_s2,xmrat,
     .   cp1o,cp2o,cs1o,cs2o)
c
c: Computes ratios of sound speeds and densities required for new
c: version of reflection coefficient subroutine that works with
c: attenuation (rtc_calc).
c
      implicit integer*4(i-n)
      complex nx_p2,nx_s1,nx_s2
c
      fac1=alp1x*cp1x
      ximfac=(fac1 - alp2x*cp2x)*fax1
      nx_p2=((cp1x/cp2x)*cmplx(1.,ximfac))**2
      ximfac=(fac1 - als1x*cs1x)*fax1
      nx_s1=((cp1x/max(1.e-6,cs1x))*cmplx(1.,ximfac))**2
      ximfac=(fac1 - als2x*cs2x)*fax1
      nx_s2=((cp1x/max(1.e-6,cs2x))*cmplx(1.,ximfac))**2
      xmrat=rho2x/rho1x
      cp1o=cp1x
      cp2o=cp2x
      cs1o=cs1x
      cs2o=cs2x
c
      return
      end
