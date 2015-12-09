      function rnglay(kprof,a,c,g0,g1,g2,z)
c
c: this function returns the horizontal progress of a ray in a layer
c: of depth z, initial sound speed cn, and gradient g0.
      implicit integer*4(i-n)
      include 'common/kwine'
c
      if(kprof .eq. 1) then
         if(g0 .eq. 0.) then
            rnglay=a*c*z/psi
         else
            cp=c + g0*z
            psipr=sqrt(max(1.-(a*cp)**2,0.))
            rnglay=(psi-psipr)/(a*g0)
         endif
      elseif(kprof .eq. 2) then
         ac=a*c
         fac=c*(1.+g1)
         cp=sign(1.,fac)*sqrt(fac**2 + 2.*g0*fac*z) - g1*c
         acp=a*cp
         psipr=sqrt(max(1.-(acp)**2,0.))
         s=(asin(ac) - asin(min(acp,1.)))/(ac)
         rnglay=(psi*(1.+2.*g1) - psipr*(cp/c+2.*g1) - s)
     .         /(2.*a*g0*(1.+g1))
      elseif(kprof .eq. 3) then
         call cdcwine(c,g0,g1,g2,z,cp,dcp,dcp2,0)
         phi2=sqrt(max(0.,1.-(a*cp)**2))/cp
         f2=1. + g2*z
         d=z/(phi1 + f2*phi2) 
         if(f1*d**2 .lt. 1.e-5) then
            s=d**2/3.
         else
            if(f1 .lt. 0.) then
               hh=2.*atan(rabf1*d)/rabf1
            else
               arg=rabf1*d
               archt=log((1.+arg)/(1.-arg))/2.
               hh=2.*archt/rabf1
            endif
            s=(hh/(2.*d) - 1.)/f1
         endif
         rnglay=2.*a*d*(1. + g2*z/2. + f3*s)
      endif
c
      return
      end 
