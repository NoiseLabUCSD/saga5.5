      subroutine rdrlay(kprof,a,ainv,c,cp,z,g0,g1,g2,
     .   r,dr,kturn)
c
c: this subroutine calculates r and r' for a ray with
c: snell inv, "a", in a layer with entry point sound speed c,
c: exit point sound velocity cp, and positive
c: thickness z.  kturn is set to 1 if the ray turns in the layer
c: and to 0 if the ray is totally reflected before entering the layer.
c: kprof=1: linear svp with g0=gradient 
c: kprof=2: blug svp with g0=initial gradient; g1=beta.
c: kprof=3: weinberg curved profile.
c
      implicit integer*4(i-n)
      include 'common/pii'
c: ray reflects at entrance to layer.
      if(ainv .le. c) then
         r=0.
         dr=0.
         kturn=0
         return
      endif
      ac=a*c
      acp=a*cp
      psi=sqrt(1.-ac**2)
c: linear sound velocity profile
      if(kprof .eq. 1) then
c: ray turns in layer.
         if(ainv .le. cp) then
            r=psi/(a*g0)
            dr=-1./(g0*psi*a**2)
            kturn=1 
c: ray passes though layer.
         else
            kturn=-1
            if(g0 .ne. 0.) then
               psipr=sqrt(1.-acp**2)
               r=(psi-psipr)/(a*g0)
               dr=r/(psi*psipr*a)
            else
               temp=c*z/psi
               r=a*temp
               dr=temp/psi**2 
            endif
         endif
c
c: blug sound velocity profile
      elseif(kprof .eq. 2) then
         fac=1.+g1
         den=a*g0*fac
c: ray turns in layer
         if(ainv .le. cp) then
            kturn=1 
            s=(asin(ac) - pieh)/ac
            r=(psi*(1.+2.*g1) - s)/(2.*den)
            dr=(s - fac/psi)/(a*den)
c: ray traverses layer
         else
            kturn=-1
            psipr=sqrt(1.-acp**2)
            fac2=cp/c + g1
            s=(asin(ac) - asin(acp))/ac 
            r=(psi*(1.+2.*g1) - psipr*(cp/c + 2.*g1) - s)/(2.*den)
            dr=(fac2/psipr - fac/psi + s)/(a*den) 
         endif
c: weinberg curved, continuous-gradient profile:  
      elseif(kprof .eq. 3) then
         phi1=psi/c 
         phi1sq=phi1**2
         f1=g1 + g2**2*phi1sq
         f1p=-2.*a*g2**2
         rabf1=sqrt(abs(f1))
         g0tg2=g0*g2
         f3=g1 - g0tg2/2.
         f4=g1 - g0tg2
         if(ainv .le. cp) then
            kturn=1 
            h=zwine(sg,z,g0,g1,g2,phi1,f1,f4)
            f8=((g0/2.)**2 - f4*phi1sq)**(.5)
            f8p=a*f4/f8
            hp=(2.*a*g2 - h*f1p + sg*f8p)/f1
            d=h/phi1
            dsq=d**2
c           dp=hp/phi1 + a*h/phi1**3
            dp=(hp + a*h/phi1sq)/phi1
         else
            kturn=-1
            phi2=sqrt(1.-acp**2)/cp
            h=z
            f2=1. + g2*h
            hp=0.
            d=h/(phi1 + f2*phi2)
            dsq=d**2
            dp=a*(1./phi1 + f2/phi2)*dsq/h
         endif
         u=rabf1*d
c     print *,'rdrlay: c,f1,f1p: ',1./a,f1,f1p
c     print *,'   h,hp,phi2: ',h,hp,phi2
c     print *,'   d,dp: ',d,dp
c     print *,'   u,u**2: ',u,u**2
c        if(u**2 .le. 1.e-5) then
         if(abs(u) .le. .003162) then
c           print *,'power series used in rdrlay '
            s=dsq/3.
            sp=2.*d*dp/3.
         else
            if(f1 .lt. 0.) then
               hh=2.*atan(u)/rabf1
            else
               archt=log((1.+u)/(1.-u))/2.
               hh=2.*archt/rabf1
            endif
            s=(hh/(2.*d) - 1.)/f1
            hhp=(2./(1.-f1*dsq))*(dp + d*f1p/(2.*f1)) - hh*f1p/(2.*f1)
            sp=divrat(hh,1.,d,hhp,0.,dp)/(2.*f1) - s*f1p/f1 
         endif
         ad2=2.*a*d
         r=ad2*(1. + g2*h/2. + f3*s) 
         dr=ad2*(g2*hp/2. + f3*sp) + r*(a*dp+d)/(a*d)
c     print *,'   s,sp: ',s,sp
c     print *,'   hh,hhp: ',hh,hhp
      endif
c
      return
      end 
