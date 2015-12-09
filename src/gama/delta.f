      subroutine delta(kind,a,rho1,cp1,cs1,rho2,cp2,cs2,
     .   r,dr,dr2,ibd) 
c
c: this subroutine calculates the beam displacement r and its
c: derivative dr from analytic expressions.
c
      implicit integer*4(i-n)
      include 'common/caustix'
      include 'common/gamaoptions'
      include 'common/bdcom'
c
      real l1,l2,l3,l4,lsum,lpsum,m1,m2,l1p,l2p,l3p,l4p,m1p,m2p
      real l4m1p,l1l3p,l2l4p,l3m2p,l2m1p,l1m2p
      real n1,n2,n3,n4,n1p,n2p,n3p,n4p
      real l1pp,l2pp,l3pp,l4pp,m1pp,m2pp,n1pp,n2pp,n3pp,n4pp
c
      if(iibd .eq. 0) then
         r=0.
         dr=0.
         ibd=0
         return
      endif
c
      ibd=1
      asq=a**2
      rh=rho2/rho1
      rhsq=rh**2
      sp1=(1. - (a*cp1)**2)/cp1**2
      sp2=((a*cp2)**2 - 1.)/cp2**2
      phip1=sqrt(sp1)
      phip2=sqrt(sp2)
      pp12=phip1*phip2
      if((cs1 .eq. 0.) .and. (cs2 .eq. 0.)) then
         u=2.*a*rh*(sp1 + sp2)
         v=pp12*(sp2 + rhsq*sp1)
         r=u/(bdomega*v)
         if(kind .eq. 1) goto 99
         up=u/a
         vp=2.*a*pp12*(1.-rhsq) + a*v*(1./sp2 - 1./sp1)
         dr=(v*up - u*vp)/(bdomega*v**2) 
         if(kind .eq. 3) then 
            vpp=2.*(1.-rhsq)*pp12*(asq*(1./sp2 - 1./sp1) + 1.)
     .         - 2.*asq*v*(1./sp2**2 + 1./sp1**2) + (a*vp+v)*
     .         (1./sp2 - 1./sp1)
            dr2=-1.*(u*vpp/bdomega + 2.*dr*v*vp)/v**2
         endif
      elseif(cs1 .eq. 0.) then
         sg=sign(1.,1. - a*cs2)
         cs2sq=cs2**2
         ss2=sg*(1. - asq*cs2sq)/cs2sq
         phis2=sqrt(ss2)
c
         m1=(1. - 2.*asq*cs2sq)/cs2sq
         l1=phip1*m1**2
         l2=sg*4.*asq*pp12*phis2
         l3=phip2/(cs2**4*rh) 
         l1p=-1.*a*(8.*phip1*m1 + l1/sp1)
         l2p=a*l2*(2./asq - 1./sp1 + 1./sp2 - sg/ss2)
         l3p=a*l3/sp2
c
         if(sg .eq. 1.) then
            u=2.*l1*l3
            v=l3**2 - l1**2 - l2**2
            up=2.*(l1*l3p+l3*l1p)
            vp=2.*(l3*l3p - l1*l1p - l2*l2p)
         else
            lsum=l1 + l2
            lpsum=l1p + l2p
            u=2.*l3*lsum
            v=l3**2 - lsum**2 
            up=2.*(l3*lpsum + l3p*lsum) 
            vp=2.*(l3*l3p - lsum*lpsum) 
         endif
c
         uvsq=u**2 + v**2
         uvp=u*vp - v*up
         r=uvp/(bdomega*uvsq)
c
         if(kind .eq. 1) goto 99
         l1pp=a*(8.*a*(4.*phip1 + m1/phip1) - (sp1*l1p+2.*a*l1)
     .      /sp1**2) + l1p/a
         l2pp=-2.*asq*l2*(2./asq**2 + 1./sp1**2 + 1./sp2**2 
     .      + 1./ss2**2) + l2p*(l2p/l2 + 1./a)
         l3pp=(sp2*(a*l3p+l3) - 2.*asq*l3)/sp2**2 
         if(sg .eq. 1.) then
            upp=2.*(l1*l3pp + 2.*l1p*l3p + l3*l1pp)
            vpp=2.*(l3*l3pp + l3p**2 - l1*l1pp - l1p**2 - l2*l2pp
     .         - l2p**2)
         else
            upp=2.*(l3*(l1pp+l2pp) + 2.*l3p*lpsum + lsum*l3pp)
            vpp=2.*(l3*l3pp + l3p**2 - lsum*(l1pp+l2pp) - lpsum**2)
         endif
         dr=(uvsq*(u*vpp-v*upp)-2.*uvp*(u*up+v*vp))/(bdomega*uvsq**2)
c
      elseif(cs2 .eq. 0.) then
         cs1sq=cs1**2
         ss1=(1. - asq*cs1sq)/cs1sq
         phis1=sqrt(ss1)
         m1=ss1 - asq
         l1=4.*asq*pp12*phis1 
         l1p=a*l1*(2./asq - 1./ss1 - 1./sp1 + 1./sp2)
         l2=rh*phip1/cs1**4
         l2p=-1.*a*l2/sp1
         l3=phip2*m1**2
         l3p=a*(l3/sp2 - 8.*phip2*m1)
         u=2.*l2*l3 
         up=2.*(l2*l3p+l3*l2p)
         v=l3**2 - l2**2 - l1**2
         vp=2.*(l3*l3p - l2*l2p - l1*l1p)
c
         uvsq=u**2 + v**2
         uvp=u*vp - v*up
         r=uvp/(bdomega*uvsq)
         if(kind .eq. 1) goto 99
         l1pp=-2.*asq*l1*(2./asq**2 + 1./ss1**2 + 1./sp1**2 +
     .      1./sp2**2) + l1p*(l1p/l1 + 1./a)
         l2pp=(-2.*asq*l2 - sp1*(a*l2p+l2))/sp1**2
         l3pp=a*((sp2*l3p-2.*a*l3)/sp2**2 + 8.*a*(4.*phip2- 
     .      m1/phip2)) + l3p/a
         upp=2.*(l2*l3pp+2.*l2p*l3p+l3*l2pp)
         vpp=2.*(l3*l3pp+l3p**2 - l2*l2pp-l2p**2 - l1*l1pp-l1p**2)
         dr=(uvsq*(u*vpp-v*upp)-2.*uvp*(u*up+v*vp))/(bdomega*uvsq**2)
c
      elseif(1 .eq. 1) then
         sg=sign(1.,1.-a*cs2) 
         cs1sq=cs1**2
         ss1=(1. - asq*cs1sq)/cs1sq
         phis1=sqrt(ss1)
         cs2sq=cs2**2
         ss2=sg*(1. - asq*cs2sq)/cs2sq
         phis2=sqrt(ss2)
         g=rho1*cs1sq - rho2*cs2sq
         q2p=2.*g
         q2=q2p*a
         q=a*q2
         qp=2.*q2
         h=rho1 - rho2
c
         l1=q + rho2
         l1p=qp
         l2=(h-q)*a 
         l2p=(h-q) - a*qp
         l3=(rho1-q)*phip2
         l3p=(rho1-q)*a/phip2 - phip2*qp
         l4=-1.*phip2*q2
         l4p=-1.*(phip2*q2p + q2*a/phip2)
         m1=phis2*q2
         m1p=phis2*q2p - sg*q2*a/phis2
         m2=(rho1-q)*phis2
         m2p=-1.*(sg*(rho1-q)*a/phis2 + phis2*qp) 
         l4m1p=l4*m1p + m1*l4p
         l1l3p=l1*l3p + l3*l1p
         l2l4p=l2*l4p + l4*l2p
         l3m2p=l3*m2p + m2*l3p
         l2m1p=l2*m1p + m1*l2p
         l1m2p=l1*m2p + m2*l1p
c
         if(sg .eq. 1.) then
            n1=phis1*phip1*l4*m1
            n1p=phis1*phip1*(l4m1p) - a*l4*m1*(phis1/phip1
     .         + phip1/phis1) 
            n2=phis1*(l1*l3 + l2*l4) + l3*m2
            n2p=phis1*((l1l3p) + (l2l4p)) - (l1*l3
     .         + l2*l4)*a/phis1 + (l3m2p)
            n3=phip1*(l1*m2 - l2*m1 + phis1*l1**2)
            n3p=phip1*((l1m2p) - (l2m1p) + 2.*phis1*
     .         l1*l1p - a*l1**2/phis1) - a*n3/sp1 
            n4=l2**2
            n4p=2.*l2*l2p
         else
            n1=phip1*(l2*m1 - l1*m2)
            n1p=phip1*((l2m1p) - (l1m2p)) - a*n1
     .         /sp1 
            n2=phis1*(l1*l3 + l2*l4)
            n2p=phis1*((l1l3p) + (l2l4p)) - a*n2
     .         /ss1 
            n3=phip1*phis1*(l1**2 + l4*m1)
            n3p=phip1*phis1*(2.*l1*l1p + (l4m1p)) - a*
     .         (l1**2 + l4*m1)*(phip1/phis1 + phis1/phip1)
            n4=l2**2 - m2*l3
            n4p=2.*l2*l2p - (l3m2p)
         endif
c
         u=-2.*(n1*n4 + n2*n3)
         up=-2.*((n1*n4p+n4*n1p) + (n2*n3p+n3*n2p))
         v=n1**2 + n3**2 - n2**2 - n4**2
         vp=2.*(n1*n1p + n3*n3p - n2*n2p - n4*n4p)
c
         uvsq=u**2 + v**2
         uvp=u*vp - v*up
         r=uvp/(bdomega*uvsq)
c
         if(kind .eq. 1) goto 99
         qpp=4.*g
         l1pp=qpp
         l2pp=-2.*qp - a*qpp
         l3pp=((q-rho1)/(cp2**2*sp2) - 2.*a*qp)/phip2 - phip2*qpp
         l4pp=(q2/(cp2**2*sp2) - 2.*a*q2p)/phip2
         m1pp=-1.*(sg*2.*a*q2p + q2/(cs2sq*ss2))/phis2
         m2pp=((q-rho1)/(cs2sq*ss2) + sg*2.*a*qp)/phis2 - phis2*qpp
         psrat=phip1/phis1
         sprat=1./psrat
         psm=phip1*phis1
         psmp=-1.*a*(psrat + sprat)
         psmpp=asq*((psrat-sprat)/sp1+(sprat-psrat)/ss1)-(psrat+sprat)
         if(sg .eq. 1.) then
            n1pp=psm*(l4*m1pp+2.*l4p*m1p+m1*l4pp) + 2.*psmp*l4m1p
     .         + l4*m1*psmpp
            n2pp=phis1*(l1*l3pp+2.*l1p*l3p+l3*l1pp + l2*l4pp+2.*l2p*l4p
     .         +l4*l2pp) - (2.*(l1l3p + l2l4p)*a + (l1*l3+l2*l4)/(cs1sq
     .         *ss1))/phis1 + (l3*m2pp+2.*l3p*m2p+m2*l3pp)
            n3pp=phip1*(l1*m2pp+2.*l1p*m2p+m2*l1pp - l2*m1pp-2.*l2p*m1p
     .         -m1*l2pp) - (2.*(l1m2p - l2m1p)*a + (l1*m2-l2*m1)/(cp1**2
     .         *sp1))/phip1 + psm*2.*(l1*l1pp+l1p**2) + 4.*psmp*l1*l1p
     .         + l1**2*psmpp
            n4pp=2.*(l2*l2pp + l2p**2)
c
         else
            n1pp=phip1*(l2*m1pp+2.*l2p*m1p+m1*l2pp - l1*m2pp-2.*l1p*m2p
     .         -m2*l1pp) - (2.*(l2m1p - l1m2p)*a + (l2*m1-l1*m2)/(cp1**2
     .         *sp1))/phip1
            n2pp=phis1*(l1*l3pp+2.*l1p*l3p+l3*l1pp + l2*l4pp+2.*l2p*l4p
     .         +l4*l2pp) - (2.*(l1l3p + l2l4p)*a + (l1*l3+l2*l4)/(cs1sq
     .         *ss1))/phis1
            n3pp=psm*(2.*(l1*l1pp+l1p**2) + l4*m1pp+2.*l4p*m1p+m1*l4pp)
     .         + 2.*psmp*(2.*l1*l1p + l4m1p) + (l1**2 + l4*m1)*psmpp
            n4pp=2.*(l2*l2pp + l2p**2) - m2*l3pp-2.*m2p*l3p-l3*m2pp
         endif
         upp=-2.*(n1*n4pp+2.*n1p*n4p+n4*n1pp + n2*n3pp+2.*n2p*n3p+
     .      n3*n2pp)
         vpp=2.*(n1*n1pp+n1p**2 + n3*n3pp+n3p**2 - n2*n2pp-n2p**2 -
     .      n4*n4pp-n4p**2)
         dr=(uvsq*(u*vpp-v*upp)-2.*uvp*(u*up+v*vp))/(bdomega*uvsq**2)
c
      endif
99    continue
c
      return
      end 
