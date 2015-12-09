      subroutine rtpres(a,rho1,ccp1,ccs1,rho2,ccp2,
     .   ccs2,rc12,tc12)
c
c: this subroutine calculates the p-p reflection coefficient, rc12,
c: and the p-p transmission coefficient, tc12, for the ray incident
c: from layer 1 into layer 2 with snell invariant "a".
c
      implicit integer*4(i-n)
      include 'common/pii'
      complex rc12,tc12
      complex l1,l2,l3,l4,m1,m2,m3,m4,d,phip1,phip2,phis1,phis2
      complex l1l3,m2m4,l2l4,m1m3,q,q2
      complex ccp1,ccs1,ccp2,ccs2,ccs1ssq,ccs2ssq 
c
      assq=a**2
      phip1=sqrt(cmplx(1.,0.) - assq*ccp1**2)/ccp1
      phip2=sqrt(cmplx(1.,0.) - assq*ccp2**2)/ccp2
      rhorat=rho2/rho1
c
      cs1=real(ccs1)
      cs2=real(ccs2)
c: liquid into liquid: 
      if((cs1 .eq. 0.) .and. (cs2 .eq. 0.)) then
         l1=phip1*rhorat
         d=l1 + phip2
         rc12=(l1 - phip2)/d
         tc12=2.*l1/d
c
c: liquid into solid: 
      elseif(cs1 .eq. 0.) then
         ccs2ssq=ccs2**2
         phis2=sqrt(cmplx(1.,0.) - assq*ccs2ssq)/ccs2
         m1=phis2**2 - cmplx(assq,0.)
         l1=phip1*m1**2
         l2=4.*assq*phip1*phip2*phis2
         l3=phip2/(ccs2ssq**2*rhorat)
         d=l1 + l2 + l3
         rc12=(l1+l2-l3)/d
         tc12=2.*m1*phip1/(d*ccs2ssq)
c
c: solid into liquid: 
      elseif(cs2 .eq. 0.) then
         ccs1ssq=ccs1**2
         phis1=sqrt(cmplx(1.,0.) - assq*ccs1ssq)/ccs1
         m1=phis1**2 - cmplx(assq,0.)
         l1=4.*assq*phip1*phip2*phis1
         l2=rhorat*phip1/ccs1ssq**2
         l3=phip2*m1**2
         d=l1+l2+l3 
c
         rc12=(l1+l2-l3)/d
         tc12=2.*l2*ccs1ssq*m1/d
c
c: solid into solid: 
      elseif((cs1 .ne. 0.) .and. (cs2 .ne. 0.)) then
         ccs1ssq=ccs1**2
         ccs2ssq=ccs2**2
         phis1=sqrt(cmplx(1.,0.) - assq*ccs1ssq)/ccs1
         phis2=sqrt(cmplx(1.,0.) - assq*ccs2ssq)/ccs2
         q2=2.*a*(rho1*ccs1ssq - rho2*ccs2ssq)
         q=a*q2
         h=rho1-rho2
c
         l1=q + cmplx(rho2,0.)
         l2=(cmplx(h,0.) - q)*a
         l3=(cmplx(rho1,0.) - q)*phip2
         l4=-1.*phip2*q2
         m1=phis2*q2
         m2=(cmplx(rho1,0.) - q)*phis2
         m3=-1.*l2
         m4=l1
c
         l1l3=phip1*l1 + l3
         m2m4=m2 + phis1*m4
         l2l4=l2 + phis1*l4
         m1m3=phip1*m1 + m3
         d=l1l3*m2m4 - l2l4*m1m3
c
         rc12=((phip1*l1-l3)*m2m4 - l2l4*(phip1*m1-m3))/d
         tc12=2.*phip1*m2m4*rho2/d
      endif
c
      return
      end 
