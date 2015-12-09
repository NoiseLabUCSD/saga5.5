      subroutine tmatlay(kprof,a,ainv,c,cp,z,g0,g1,g2,
     .   akp,dakp,expfac,travt,kturn)
c
c: this subroutine calculates the travel time and frequency-dependent 
c: attenuation factor for a ray with snell invariant a in a layer
c: with initial s.v. c, final s.v. cp, thickness z, attenuation akp
c: (db/m-khz), and attenuation gradient dakp.
c
      implicit integer*4(i-n)
      include 'common/pii'
      logical attn
c: ray does not enter the layer.
      if(ainv .le. c) then
         expfac=0.
         travt=0.
         kturn=0
         return
      endif
c
c: flag for ocean layers where attenuation is zero: 
      if((akp .eq. 0.) .and. (dakp .eq. 0.)) then
         expfac=0.
         attn=.false.
      else
         attn=.true.
      endif
      ac=a*c
      acp=a*cp
      psi=sqrt(1.-(ac)**2)
c: linear profile:  
      if(kprof .eq. 1) then
c: ray turns in the layer.
         ag0=a*g0
         if(ainv .le. cp) then
            if(attn) then
               fac1=(akp-dakp*c/g0)/(ag0)
               fac2=dakp/(ag0)**2
               expfac=fac1*acos(ac) + fac2*psi
            endif
            travt=log((1.+psi)/(ac))/g0
            kturn=1 
c: ray passes through the layer.
         else
            kturn=-1
            if(g0 .ne. 0.) then
               psipr=sqrt(1.-(acp)**2)
               if(attn) then
                  fac1=(akp-dakp*c/g0)/(ag0)
                  fac2=dakp/(ag0)**2
                  expfac=fac1*(asin(acp)-asin(ac)) + fac2*(psi-psipr) 
               endif
               travt=log(cp*(1.+psi)/(c*(1.+psipr)))/g0
            else
               if(attn) then
                  expfac=(akp*z + dakp*(z**2)/2.)/psi
               endif
               travt=z/(c*psi)
            endif
         endif
c: blug profile: 
      elseif(kprof .eq. 2) then
         fac=1.+g1
         if(attn) then
            f1=g1*c 
            f12=f1**2
            f2=c*fac
            f22=f2**2
            f3=1./a**2
            f4=akp/(a*g0*f2)
            q=dakp/(2.*a*g0**2*f22)
            f9=f1*(q*(f12 - f22 + 1.5*f3) + f4)
            f8=(q*(f3 + 3.*f12 - f22) + f4)/a
            f7=q*f3/(3.*a)
            f6=3.*q*f1/(2.*a) 
         endif
c
         if(ainv .le. cp) then
            kturn=1 
            s=pieh - asin(ac) 
            if(attn) then
               expfac=s*f9 + psi*f8 - psi**3*f7 + c*psi*f6
            endif
            travt=(s/(ac) + g1*log((1.+psi)/(ac)))/(g0*fac)
         else
            kturn=-1
            psipr=sqrt(1.-(acp)**2)
            s=asin(acp) - asin(ac)
c: 8-14-89: expfac goes to infinity for snell inv very small
            if(attn) then
               if(acp .lt. .017) then
                  expfac=akp*z + .5*dakp*z**2
               else 
                  expfac=s*f9+(psi-psipr)*f8-(psi**3-psipr**3)*f7
     .               + (c*psi-cp*psipr)*f6
               endif
            endif
            travt=(s/(ac) + g1*log(cp*(1.+psi)/(c*(1.+psipr))))
     .         /(g0*fac)
         endif
c: curvilinear profile: 
      elseif(kprof .eq. 3) then
         if(attn) then
c: approx curved profile in ocean with linear. dakp=0 in ocean:
            g00=(cp - c)/z
            if(ainv .le. cp) then
               expfac=akp*acos(ac)/(a*g00)
            else
               if(g00 .ne. 0.) then
                  expfac=akp*(asin(acp)-asin(ac))/(a*g00)
               else
                  expfac=akp*z/psi
               endif
            endif
         endif
         phi1=psi/c 
         f1=g1 + g2**2*phi1**2
         rabf1=sqrt(abs(f1))
         f3=g1 - g0*g2/2.
         f4=g1 - g0*g2
         rabf4=sqrt(abs(f4))
         if(ainv .le. cp) then
            kturn=1 
            phi2=0. 
            h=zwine(sg,z,g0,g1,g2,phi1,f1,f4)
            f2=1.+g2*h
            d=h/phi1
         else
            kturn=-1
            phi2=sqrt(1.-(acp)**2)/cp
            h=z
            f2=1.+g2*h
            d=h/(phi1+f2*phi2)
         endif
         u=rabf1*d
         if(u**2 .le. 1.e-5) then
            s=d**2/3.
            hh=2.*d*(1. + f1*d**2/3.)
         else
            if(f1 .lt. 0.) then
               hh=2.*atan(u)/rabf1
            else
               archt=log((1.+u)/(1.-u))/2.
               hh=2.*archt/rabf1
            endif
            s=(hh/(2.*d) - 1.)/f1
         endif
         r=2.*a*d*(1. + g2*h/2. + f3*s) 
         f5=g0/2. + f3*h
         f9=rabf4*(g0*f2*phi2/2. - f5*phi1)/(g0*f5/2. - f2*f4*
     .      phi1*phi2)
         if(f4 .gt. 0.) then
            archt=log((1.+f9)/(1.-f9))/2.
            e=rabf4*archt
         else
            e=rabf4*atan(f9)
         endif
         travt=r/(a*c**2) + (g1*r/a - 2.*f3*hh + e)/g2**2
      endif
c
      return
      end 
