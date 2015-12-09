      subroutine dr2lay(kprof,a,ainv,c,cp,z,g0,g1,g2,r,dr,dr2,kturn)
c
c: this subroutine calculates r, r', and r" for a ray with si
c: "a".  see rdrlay for additional comments.
c
      implicit integer*4(i-n)
      include 'common/pii'
c: ray reflects at entrance to layer.
      if(ainv .le. c) then
         r=0.
         dr=0.
         dr2=0.
         kturn=0
         return
      endif
      ac=a*c
      acp=a*cp
      psi=sqrt(1.-(ac)**2)
c: linear sound velocity profile
      if(kprof .eq. 1) then
c: ray turns in layer.
         if(ainv .le. cp) then
            r=psi/(a*g0)
            dr=-1./(g0*psi*a**2)
            dr2=(3.*psi**2 - 1.)/(a**3*g0*psi**3) 
            kturn=1 
c: ray passes though layer.
         else
            kturn=-1
            if(g0 .ne. 0.) then
               psipr=sqrt(1.-(acp)**2)
               r=(psi-psipr)/(a*g0)
               dr=r/(psi*psipr*a)
               pd=psi-psipr
               pp=psi*psipr
               dr2=pd*(3.*pp*(1.-pp) + pd**2)/(g0*(a*pp)**3)
            else
               r=ac*z/psi
               dr=c*z/psi**3
               dr2=3.*a*c**3*z/psi**5
            endif
         endif
c
c: blug sound velocity profile
      elseif (kprof .eq. 2) then
         fac=1.+g1
         den=a*g0*fac
c: ray turns in layer
         if(ainv .le. cp) then
            kturn=1 
            s=(asin(ac) - pieh)/(ac)
            r=(psi*(1.+2.*g1) - s)/(2.*den)
            dr=(s - fac/psi)/(a*den)
            dr2=((3./psi - 1./psi**3)*fac + 1./psi - 3.*s)/(a**2*den) 
c: ray traverses layer
         else
            kturn=-1
            psipr=sqrt(1.-(acp)**2)
            fac2=cp/c + g1
            s=(asin(ac) - asin(acp))/(ac)
            r=(psi*(1.+2.*g1) - psipr*(cp/c + 2.*g1) - s)/(2.*den)
            dr=(fac2/psipr - fac/psi + s)/(a*den) 
            dr2=((3./psi - 1./psi**3)*fac - (3./psipr -
     .         1./psipr**3)*fac2 + 1./psi - cp/(c*psipr) -
     .         3.*s)/(a**2*den)
         endif
c: curvilinear profile: 
      elseif(kprof .eq. 3) then
         phi1=psi/c 
         f1=g1 + g2**2*phi1**2
         f1p=-2.*a*g2**2
         f1pp=f1p/a 
         rabf1=sqrt(abs(f1))
         f3=g1 - g0*g2/2.
         f4=g1 - g0*g2
         if(ainv .le. cp) then
            kturn=1 
            h=zwine(sg,z,g0,g1,g2,phi1,f1,f4)
            f8=((g0/2.)**2 - f4*phi1**2)**(.5)
            f8p=a*f4/f8
            f8pp=f4*(f8 - a*f8p)/f8**2
            hp=(2.*a*g2 - h*f1p + sg*f8p)/f1
            hpp=(2.*(g2 - hp*f1p) - h*f1pp + sg*f8pp)/f1
            d=h/phi1
            dp=hp/phi1 + a*h/phi1**3
            dpp=hpp/phi1 + (a*(2.*hp + 3.*a*d**2/h) + h)/phi1**3
         else
            kturn=-1
            phi2=sqrt(1.-(acp)**2)/cp
            h=z
            f2=1.+g2*h
            f9=1./phi1 + f2/phi2
            f9p=a*(1./phi1**3 + f2/phi2**3)
            d=h/(phi1 + f2*phi2)
            dp=a*f9*d**2/h
            dpp=(2.*a*f9*d*dp + d**2*(a*f9p+f9))/h
            hp=0.
            hpp=0.
         endif
         u=rabf1*d
         if(u**2 .le. 1.e-5) then
            s=d**2/3.
            sp=2.*d*dp/3.
            spp=2.*(d*dpp + dp**2)/3.
         else
            if(f1 .lt. 0.) then
               hh=2.*atan(u)/rabf1
            else
               archt=log((1.+u)/(1.-u))/2.
               hh=2.*archt/rabf1
            endif
            s=(hh/(2.*d) - 1.)/f1
            f1d2=1. - f1*d**2 
            f1dp=2.*f1*dp + f1p*d
            hhp=(f1dp/f1d2 - hh*f1p/2.)/f1
            sp=((d*hhp-hh*dp)/(2.*d**2) - s*f1p)/f1
            hhpp=(2./f1d2)*(dpp+divrat(d,f1p,f1,dp,f1pp,f1p)/2.) +
     .         (d/f1)*(f1dp/f1d2)**2 - divrat(hh,f1p,f1,hhp,f1pp,f1p)/2.
            spp=((d*hhpp-hh*dpp) - (d*hhp-hh*dp)*f1dp/(f1*d))/(2.*
     .         f1*d**2) - divrat(s,f1p,f1,sp,f1pp,f1p)
         endif
         atd=a*d
         adpd=a*dp + d
         r=2.*atd*(1. + g2*h/2. + f3*s) 
         dr=2.*atd*(g2*hp/2. + f3*sp) + r*adpd/(atd)
         dr2=2.*atd*(g2*hpp/2. + f3*spp) + adpd*(g2*hp + 2.*f3*sp
     .      + (atd*dr - r*adpd)/(atd)**2) + r*(a*dpp + 2.*dp)/(atd)
c
      endif
c
      return
      end 
