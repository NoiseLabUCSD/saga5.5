      subroutine pwreal(a,rho1,cp1,cs1,rho2,cp2,
     .   cs2,rc12,tc12,rc21,tc21) 
c
c: this subroutine calculates the p-p reflection coefficients, rc12
c: and rc21, and the transmission coeffs, tc12 and tc21 for the plane 
c: wave with snell invariant a=cos(theta)/c
c
      implicit integer*4(i-n)
      complex rc12,tc12,rc21,tc21
      complex l1,l2,l3,l4,m1,m2,m3,m4,d,phip1,phip2,phis1,phis2
      complex l1l3,m2m4,l2l4,m1m3,q,q2
c
      assq=a**2
      phip1=sqrt(dcmplx(1.d0,0.d0) - assq*cp1**2)/cp1
      phip2=sqrt(dcmplx(1.d0,0.d0) - assq*cp2**2)/cp2
      rhorat=rho2/rho1
c
c: liquid into liquid: 
      if((cs1 .eq. 0.) .and. (cs2 .eq. 0.)) then
         l1=phip1*rhorat
         d=l1 + phip2
         rc12=(l1 - phip2)/d
         tc12=dcmplx(1.d0,0.d0) + rc12
         rc21=-1.*rc12
         tc21=dcmplx(1.d0,0.d0) + rc21
c
c: liquid into solid: 
      elseif(cs1 .eq. 0.) then
         cs2ssq=cs2**2
         phis2=sqrt(dcmplx(1.d0,0.d0) - assq*cs2ssq)/cs2
         m1=phis2**2 - dcmplx(assq,0.0e0)
         l1=phip1*m1**2
         l2=4.*assq*phip1*phip2*phis2
         l3=phip2/(cs2ssq**2*rhorat)
         d=l1 + l2 + l3
         rc12=(l1+l2-l3)/d
         tc12=2.*m1*phip1/(d*cs2ssq)
c
         rc21=(l2+l3-l1)/d
         tc21=2.*l3*cs2ssq*m1/d
c
c: solid into liquid: 
      elseif(cs2 .eq. 0.) then
         cs1ssq=cs1**2
         phis1=sqrt(dcmplx(1.d0,0.d0) - assq*cs1ssq)/cs1
         m1=phis1**2 - dcmplx(assq,0.e0)
         l1=4.*assq*phip1*phip2*phis1
         l2=rhorat*phip1/cs1ssq**2
         l3=phip2*m1**2
         d=l1+l2+l3 
         rc12=(l1+l2-l3)/d
         tc12=2.*l2*cs1ssq*m1/d
c
         rc21=(l3+l1-l2)/d
         tc21=2.*m1*phip2/(d*cs1ssq)
c
c: solid into solid: 
      elseif((cs1 .ne. 0.) .and. (cs2 .ne. 0.)) then
         cs1ssq=cs1**2
         cs2ssq=cs2**2
         phis1=sqrt(dcmplx(1.d0,0.d0) - assq*cs1ssq)/cs1
         phis2=sqrt(dcmplx(1.d0,0.d0) - assq*cs2ssq)/cs2
         q2=2.*a*(rho1*cs1ssq - rho2*cs2ssq)
         q=a*q2
         h=rho1-rho2
c
         l1=q + dcmplx(rho2,0.e0)
         l2=(dcmplx(h,0.e0) - q)*a
         l3=(dcmplx(rho1,0.e0) - q)*phip2
         l4=-1.*phip2*q2
         m1=phis2*q2
         m2=(dcmplx(rho1,0.e0) - q)*phis2
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
c
         q2=-1.*q2
         q=a*q2
         h=-1.*h
         l1=q + dcmplx(rho1,0.e0)
         l2=-1.*l2
         l3=(dcmplx(rho2,0.e0) - q)*phip1
         l4=-1.*phip1*q2
         m1=phis1*q2
         m2=(dcmplx(rho2,0.e0) - q)*phis1
         m3=-1.*l2
         m4=l1
c
         l1l3=phip2*l1 + l3
         m2m4=m2 + phis2*m4
         l2l4=l2 + phis2*l4
         m1m3=phip2*m1 + m3
         d=l1l3*m2m4 - l2l4*m1m3
c
         rc21=((phip2*l1-l3)*m2m4 - l2l4*(phip2*m1-m3))/d
         tc21=2.*phip2*m2m4*rho1/d
      endif
c
      return
      end 
