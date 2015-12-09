       subroutine rp_slay(Rmat,p,q,u,v,e11,e12,e21,e22,
     .   ze,s11,s12,s21,s22,zs,gami1,beti1,gami2,beti2,
     .   isoe,isos,Vmat,ndv,Wmat,iiw,rhorat)
c
c: Finds the four reflection coefficients referenced to the top of
c: layer 1 given the four reflection coefficients at the top of layer
c: layer 2, the (p,q,u,v) constants for the 1-2 layer interface, and
c: the compressional and shear wave propagator matrices e and s.
c: gam1,gam2 and bet1,bet2 are the vertical comp and shear wavenumbers 
c: at the top of layers 1 & 2.
c: Note that ze and zs are exponential factors common to eij and sij.
c: vps and vsp need to account for them.
c
      implicit complex*16(a-z)
      integer*4 isoe,isos,ndv,j,iiw
      complex*16 e11(3),e12(3),e21(3),e22(3),s11(3),s12(3),s21(3),
     .   s22(3),Rmat(3,5),p(3),q(3),u(3),v(3),Vmat(3,5),Wmat(6)
      complex*16 gami1(3),beti1(3),gami2(3),beti2(3),zzero
      real*8 rhorat,magsq,deps
      data zzero/(0.d0,0.d0)/,deps/1.d-50/
c
cxx   if(Rmat(1,4) .eq. zzero .or. Rmat(1,2) .eq. zzero) then
      if(magsq(Rmat(1,2)) .lt. deps .or.
     .   magsq(Rmat(1,4)) .lt. deps) then
         call rp_slay0(Rmat,p,q,u,v,e11,e12,e21,e22,ze,s11,s12,
     .      s21,s22,zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,
     .      ndv,Wmat,iiw,rhorat)
         return
      endif
c
      q1=1.d0 + Rmat(1,1)
      q2diff=1.d0 - Rmat(1,1)
      q2=gami2(1)*q2diff
      q3=Rmat(1,2)
      q4=-beti2(1)*Rmat(1,2)
      r1=Rmat(1,4)
      r2=-gami2(1)*Rmat(1,4)
      r3=1.d0 + Rmat(1,3)
      r4diff=1.d0 - Rmat(1,3)
      r4=beti2(1)*r4diff
c
c: rps and rsp are multiplied by exp(zsp) (see r1,r2,q3,q4):
      spfac=cdexp(Rmat(1,5))
      r1sp=spfac*r1
      r2sp=spfac*r2
      q3sp=spfac*q3
      q4sp=spfac*q4
c: Multiply out to get simpler form:
cxx1  delta1=q2*r1 - q1*r2
cxx1  delta2=q4*r3 - q3*r4
      delta1=-2.d0*r2
      delta2=2.d0*q4
      delta3=q4sp*r1sp - q1*r4
      delta4=q3sp*r2sp - q2*r3
c
      f1=q4sp*p(1) + q1*v(1)
      f2=q4sp*u(1) - q1*q(1)
      f3=q3sp*q(1) + q2*u(1)
      f4=q3sp*v(1) - q2*p(1)
c: Note: I switched x2 and x3 from notes:
      x1=delta4*f1
      x2=delta3*f3
      x3=delta4*f2
      x4=delta3*f4
c
      del2q1=delta2*q1
      del2q2=delta2*q2
      del1q3=delta1*q3sp
      del1q4=delta1*q4sp
      w1=del2q2*p(1) - del1q3*v(1)
      w2=del2q1*q(1) - del1q4*u(1)
      w3=del2q2*u(1) + del1q3*q(1)
      w4=del2q1*v(1) + del1q4*p(1)
c
      if(isoe .eq. 0) then
         e12g=gami1(1)*e12(1)
         e22g=gami1(1)*e22(1)
         e1p=e11(1) + e12g
         e1m=e11(1) - e12g
         e2p=e21(1) + e22g
         e2m=e21(1) - e22g
      else
         e1p=e11(1)
         e1m=e12(1)
         e2p=e21(1)
         e2m=e22(1)
      endif
      n1=w1*e1p - w2*e2p
      n2=x1*e1p - x2*e2p
      k1=w1*e1m - w2*e2m
      k2=x1*e1m - x2*e2m
c
      if(isos .eq. 0) then
         s12b=beti1(1)*s12(1)
         s22b=beti1(1)*s22(1)
         s1p=s11(1) + s12b
         s1m=s11(1) - s12b
         s2p=s21(1) + s22b
         s2m=s21(1) - s22b
      else
         s1p=s11(1)
         s1m=s12(1)
         s2p=s21(1)
         s2m=s22(1)
      endif
      m1=w3*s2p - w4*s1p
      m2=x3*s2p - x4*s1p
      y1=w3*s2m - w4*s1m
      y2=x3*s2m - x4*s1m
c
      den=(y2*k1 - y1*k2)
      deninv=1.d0/den
      num1=-(y2*n1 - y1*n2)
      Vmat(1,1)=num1*deninv
      num2=(w1*x2 - w2*x1)
      gamnum2=gami1(1)*num2
      Vmat(1,2)=2.d0*gamnum2*deninv
      num3=(k2*m1 - k1*m2)
      Vmat(1,3)=num3*deninv
      num4=(w4*x3 - w3*x4)
      betnum4=beti1(1)*num4
      Vmat(1,4)=-2.d0*betnum4*deninv
c
c: Exponential factor [exp(zexp_sp)] to be multiplied by vps and vsp:
      Vmat(1,5)=-(ze+zs)
c
      if(iiw .eq. 1) then
         a7=e1m*f4 - e2m*f2
         a8=e2m*(r1sp*q(1) - r4*u(1)) - e1m*(r2sp*p(1) - r3*v(1))
         b7=s1m*f1 - s2m*f3
         b8=s1m*(r1sp*v(1) + r4*p(1)) - s2m*(r2sp*u(1) + r3*q(1))
         wden=a7*b8 - b7*a8
         fac=-2.d0*rhorat/wden
cxx      a9_den=fac*gami1(1)*cdexp(-ze)
         a9_den=fac*gami1(1)
         Wmat(1)=a9_den*b8
         Wmat(2)=-a9_den*b7
c: Include separate exponential factor for Wpp and Wps:
         Wmat(5)=-ze
cxx      a9bar_den=fac*beti1(1)*cdexp(-zs)
         a9bar_den=fac*beti1(1)
         Wmat(3)=-a9bar_den*a7
         Wmat(4)=a9bar_den*a8
c: Include separate exponential factor for Wss and Wsp:
         Wmat(6)=-zs
cxx   print *,'prop 1',q1*Wmat(1)+r1sp*Wmat(2),cdexp(ze)*p(1)*
cxx  .   (e1p+e1m*Vmat(1,1))+cdexp(-ze)*u(1)*s2m*Vmat(1,2)
cxx   print *,'prop 2',q2*Wmat(1)+r2sp*Wmat(2),cdexp(ze)*q(1)*
cxx  .   (e2p+e2m*Vmat(1,1))+cdexp(-ze)*v(1)*s1m*Vmat(1,2)
cxx   print *,'prop 3',q3sp*Wmat(1)+r3*Wmat(2),-cdexp(ze)*u(1)*
cxx  .   (e2p+e2m*Vmat(1,1))+cdexp(-ze)*p(1)*s1m*Vmat(1,2)
cxx   print *,'prop 4',q4sp*Wmat(1)+r4*Wmat(2),-cdexp(ze)*v(1)*
cxx  .   (e1p+e1m*Vmat(1,1))+cdexp(-ze)*q(1)*s2m*Vmat(1,2)
cxx   print *,'prop 5',q3sp*Wmat(4)+r3*Wmat(3),cdexp(zs)*p(1)*
cxx  .   (s1p+s1m*Vmat(1,3))-cdexp(-zs)*u(1)*e2m*Vmat(1,4)
cxx   print *,'prop 6',q4sp*Wmat(4)+r4*Wmat(3),cdexp(zs)*q(1)*
cxx  .   (s2p+s2m*Vmat(1,3))-cdexp(-zs)*v(1)*e1m*Vmat(1,4)
cxx   print *,'prop 7',q1*Wmat(4)+r1sp*Wmat(3),cdexp(zs)*u(1)*
cxx  .   (s2p+s2m*Vmat(1,3))+cdexp(-zs)*p(1)*e1m*Vmat(1,4)
cxx   print *,'prop 8',q2*Wmat(4)+r2sp*Wmat(3),cdexp(zs)*v(1)*
cxx  .   (s1p+s1m*Vmat(1,3))+cdexp(-zs)*q(1)*e2m*Vmat(1,4)
      endif
c
c: If first derivatives desired, compute them now:
      do j=2,ndv
         q1p=Rmat(j,1)
         q2p=-gami2(1)*Rmat(j,1) + gami2(j)*q2diff
         q3p=Rmat(j,2)
         q4p=-beti2(1)*Rmat(j,2) - beti2(j)*Rmat(1,2)
         r1p=Rmat(j,4)
         r2p=-gami2(1)*Rmat(j,4) - gami2(j)*Rmat(1,4)
         r3p=Rmat(j,3)
         r4p=-beti2(1)*Rmat(j,3) + beti2(j)*r4diff
         r1spp=spfac*r1p
         r2spp=spfac*r2p
         q3spp=spfac*q3p
         q4spp=spfac*q4p
c
         delta1p=-2.d0*r2p
         delta2p=2.d0*q4p
         delta3p=q4sp*r1spp + q4spp*r1sp - q1*r4p - q1p*r4
         delta4p=q3sp*r2spp + q3spp*r2sp - q2*r3p - q2p*r3
c
         f1p=q4sp*p(j) + q4spp*p(1) + q1*v(j) + q1p*v(1)
         f2p=q4sp*u(j) + q4spp*u(1) - q1*q(j) - q1p*q(1)
         f3p=q3sp*q(j) + q3spp*q(1) + q2*u(j) + q2p*u(1)
         f4p=q3sp*v(j) + q3spp*v(1) - q2*p(j) - q2p*p(1)
         x1p=delta4*f1p + delta4p*f1
         x2p=delta3*f3p + delta3p*f3
         x3p=delta4*f2p + delta4p*f2
         x4p=delta3*f4p + delta3p*f4
c
         del2q1p=delta2*q1p + delta2p*q1
         del2q2p=delta2*q2p + delta2p*q2
         del1q3p=delta1*q3spp + delta1p*q3sp
         del1q4p=delta1*q4spp + delta1p*q4sp
         w1p=del2q2*p(j) + del2q2p*p(1) - del1q3*v(j) - del1q3p*v(1)
         w2p=del2q1*q(j) + del2q1p*q(1) - del1q4*u(j) - del1q4p*u(1)
         w3p=del2q2*u(j) + del2q2p*u(1) + del1q3*q(j) + del1q3p*q(1)
         w4p=del2q1*v(j) + del2q1p*v(1) + del1q4*p(j) + del1q4p*p(1)
         if(isoe .eq. 0) then
            e12gp=gami1(1)*e12(j) + gami1(j)*e12(1)
            e22gp=gami1(1)*e22(j) + gami1(j)*e22(1)
            e1pp=e11(j) + e12gp
            e1mp=e11(j) - e12gp
            e2pp=e21(j) + e22gp
            e2mp=e21(j) - e22gp
         else
            e1pp=e11(j)
            e1mp=e12(j)
            e2pp=e21(j)
            e2mp=e22(j)
         endif
         n1p=w1*e1pp + w1p*e1p - w2*e2pp - w2p*e2p
         n2p=x1*e1pp + x1p*e1p - x2*e2pp - x2p*e2p
         k1p=w1*e1mp + w1p*e1m - w2*e2mp - w2p*e2m
         k2p=x1*e1mp + x1p*e1m - x2*e2mp - x2p*e2m
c
         if(isos .eq. 0) then
            s12bp=beti1(1)*s12(j) + beti1(j)*s12(1)
            s22bp=beti1(1)*s22(j) + beti1(j)*s22(1)
            s1pp=s11(j) + s12bp
            s1mp=s11(j) - s12bp
            s2pp=s21(j) + s22bp
            s2mp=s21(j) - s22bp
         else
            s1pp=s11(j)
            s1mp=s12(j)
            s2pp=s21(j)
            s2mp=s22(j)
         endif
         m1p=w3*s2pp + w3p*s2p - w4*s1pp - w4p*s1p
         m2p=x3*s2pp + x3p*s2p - x4*s1pp - x4p*s1p
         y1p=w3*s2mp + w3p*s2m - w4*s1mp - w4p*s1m
         y2p=x3*s2mp + x3p*s2m - x4*s1mp - x4p*s1m
c
         denp=(y2*k1p + y2p*k1 - y1*k2p - y1p*k2)
         deninvp=-denp*deninv*deninv
         num1p=-(y2*n1p + y2p*n1 - y1*n2p - y1p*n2)
         Vmat(j,1)=num1*deninvp + num1p*deninv
         num2p=(w1*x2p + w1p*x2 - w2*x1p - w2p*x1)
         gamnum2p=gami1(1)*num2p + gami1(j)*num2
         Vmat(j,2)=2.d0*(gamnum2*deninvp + gamnum2p*deninv)
         num3p=(k2*m1p + k2p*m1 - k1*m2p - k1p*m2)
         Vmat(j,3)=num3*deninvp + num3p*deninv
         num4p=(w4*x3p + w4p*x3 - w3*x4p - w3p*x4)
         betnum4p=beti1(1)*num4p + beti1(j)*num4
         Vmat(j,4)=-2.d0*(betnum4*deninvp + betnum4p*deninv)
      enddo
c
      return
      end
