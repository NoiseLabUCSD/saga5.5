      subroutine rp_slay0(Rmat,p,q,u,v,e11,e12,e21,e22,ze,s11,s12,
     .   s21,s22,zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,
     .   Wmat,iiw,rhorat)
c
c: Finds the four reflection coefficients and their derivatives
c: referenced to the top of layer 1 given the (p,q,u,v) constants 
c: for the 1-2 layer interface, and the compressional and shear 
c: wave propagator matrices e and s.  Same as rp_slay, except that
c: Rps=Rsp=0.  When Rpp=Rss=0 also, layer 2 is a halfspace, and slightly
c: shorter formulas are used.
c: gam1,gam2 and bet1,bet2 are the vertical comp and shear wavenumbers 
c: at the top of layers 1 & 2.
c
      implicit complex*16 (a-z)
      integer*4 isoe,isos,ndv,j,iiw,iiz
      complex*16 Rmat(3,5),e11(3),e12(3),e21(3),e22(3),s11(3),
     .   s12(3),s21(3),s22(3),p(3),q(3),u(3),v(3),Vmat(3,5),
     .   gami1(3),gami2(3),beti1(3),beti2(3),Wmat(6),zzero
      real*8 rhorat
      data zzero/(0.d0,0.d0)/
c
      if(Rmat(1,1) .eq. zzero .and. Rmat(1,3) .eq. zzero) then
         iiz=1
         q1_q=q(1)
         q1_v=v(1)
         r3_q=q(1)
         r3_v=v(1)
         q2_p=gami2(1)*p(1)
         q2_u=gami2(1)*u(1)
         r4_p=beti2(1)*p(1)
         r4_u=beti2(1)*u(1)
      else
         iiz=0
         q1=1.d0 + Rmat(1,1)
         q2diff=1.d0 - Rmat(1,1)
         q2=gami2(1)*q2diff
         r3=1.d0 + Rmat(1,3)
         r4diff=1.d0 - Rmat(1,3)
         r4=beti2(1)*r4diff
         q1_q=q1*q(1)
         q1_v=q1*v(1)
         r3_q=r3*q(1)
         r3_v=r3*v(1)
         q2_p=q2*p(1)
         q2_u=q2*u(1)
         r4_p=r4*p(1)
         r4_u=r4*u(1)
      endif
c
      L=2.d0*(r4_u*q2_p - q1_q*r3_v)
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
      m1=q2_p*e1p - q1_q*e2p
      f1=q2_p*e1m - q1_q*e2m
      m2=r4_u*e2p - r3_v*e1p
      f2=r4_u*e2m - r3_v*e1m
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
      t1=q2_u*s2p - q1_v*s1p
      c1=q2_u*s2m - q1_v*s1m
      t2=r4_p*s1p - r3_q*s2p
      c2=r4_p*s1m - r3_q*s2m
c
      den=(c2*f1 + c1*f2)
      deninv=1./den
      num1=-(c2*m1 + c1*m2)
      Vmat(1,1)=num1*deninv
      gamL=gami1(1)*L
      Vmat(1,2)=gamL*deninv
      num2=-(t1*f2 + t2*f1)
      Vmat(1,3)=num2*deninv
      betL=beti1(1)*L
      Vmat(1,4)=-betL*deninv
c
c: vps(1:2) and vsp(1:2) are to be multiplied by exp(zexp_sp):
      Vmat(1,5)=-(ze+zs)
c
      if(iiw .eq. 1) then
c: Take from rp_prop, but set r1=r2=q3=q4=0:
         a7=e2m*q1_q - e1m*q2_p
         a8=e1m*r3_v - e2m*r4_u
         b7=s1m*q1_v - s2m*q2_u
         b8=s1m*r4_p - s2m*r3_q
         wden=a7*b8 - b7*a8
cxx      fac=-2.d0*(p(1)*q(1) + u(1)*v(1))/wden
         fac=-2.d0*rhorat/wden
         a9_den=fac*gami1(1)
         Wmat(1)=a9_den*b8
         Wmat(2)=-a9_den*b7
         Wmat(5)=-ze
         a9bar_den=fac*beti1(1)
         Wmat(3)=-a9bar_den*a7
         Wmat(4)=a9bar_den*a8
         Wmat(6)=-zs
cxx   print *,'Wmat 1',q1*Wmat(1),cdexp(ze)*p(1)*
cxx  .   (e1p+e1m*Vmat(1,1))+cdexp(-ze)*u(1)*s2m*Vmat(1,2)
cxx   print *,'Wmat 2',q2*Wmat(1),cdexp(ze)*q(1)*
cxx  .   (e2p+e2m*Vmat(1,1))+cdexp(-ze)*v(1)*s1m*Vmat(1,2)
cxx   print *,'Wmat 3',r3*Wmat(2),-cdexp(ze)*u(1)*(e2p+e2m*
cxx  .   Vmat(1,1))+cdexp(-ze)*p(1)*s1m*Vmat(1,2)
cxx   print *,'Wmat 4',r4*Wmat(2),-cdexp(ze)*v(1)*(e1p+e1m*
cxx  .   Vmat(1,1))+cdexp(-ze)*q(1)*s2m*Vmat(1,2)
cxx   print *,'Wmat 5 ',q1*Wmat(4),cdexp(zs)*u(1)*(s2p + s2m*
cxx  .   Vmat(1,3)) + cdexp(-zs)*p(1)*e1m*Vmat(1,4)
cxx   print *,'Wmat 6 ',q2*Wmat(4),cdexp(zs)*v(1)*(s1p + s1m*
cxx  .   Vmat(1,3)) + cdexp(-zs)*q(1)*e2m*Vmat(1,4)
cxx   print *,'Wmat 7 ',r3*Wmat(3),cdexp(zs)*p(1)*(s1p + s1m*
cxx  .   Vmat(1,3)) - cdexp(-zs)*u(1)*e2m*Vmat(1,4)
cxx   print *,'Wmat 8 ',r4*Wmat(3),cdexp(zs)*q(1)*(s2p + s2m*
cxx  .   Vmat(1,3)) - cdexp(-zs)*v(1)*e1m*Vmat(1,4)
      endif
c
      do j=2,ndv
         if(iiz .eq. 1) then
            q1_qp=q(j)
            q1_vp=v(j)
            r3_qp=q(j)
            r3_vp=v(j)
            q2_pp=gami2(1)*p(j) + gami2(j)*p(1)
            q2_up=gami2(1)*u(j) + gami2(j)*u(1)
            r4_pp=beti2(1)*p(j) + beti2(j)*p(1)
            r4_up=beti2(1)*u(j) + beti2(j)*u(1)
         else
            q1p=Rmat(j,1)
            q2diffp=-Rmat(j,1)
            q2p=gami2(1)*q2diffp + gami2(j)*q2diff
            r3p=Rmat(j,3)
            r4diffp=-Rmat(j,3)
            r4p=beti2(1)*r4diffp + beti2(j)*r4diff
            q1_qp=q1*q(j) + q1p*q(1)
            q1_vp=q1*v(j) + q1p*v(1)
            r3_qp=r3*q(j) + r3p*q(1)
            r3_vp=r3*v(j) + r3p*v(1)
            q2_pp=q2*p(j) + q2p*p(1)
            q2_up=q2*u(j) + q2p*u(1)
            r4_pp=r4*p(j) + r4p*p(1)
            r4_up=r4*u(j) + r4p*u(1)
         endif
c
         Lp=2.d0*(r4_u*q2_pp + r4_up*q2_p - q1_q*r3_vp - q1_qp*r3_v)
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
         m1p=q2_p*e1pp + q2_pp*e1p - q1_q*e2pp - q1_qp*e2p
         m2p=r4_u*e2pp + r4_up*e2p - r3_v*e1pp - r3_vp*e1p
         f1p=q2_p*e1mp + q2_pp*e1m - q1_q*e2mp - q1_qp*e2m
         f2p=r4_u*e2mp + r4_up*e2m - r3_v*e1mp - r3_vp*e1m
c
         if(isos .eq. 0) then
            bets12p=beti1(1)*s12(j) + beti1(j)*s12(1)
            bets22p=beti1(1)*s22(j) + beti1(j)*s22(1)
            s1pp=s11(j) + bets12p
            s1mp=s11(j) - bets12p
            s2pp=s21(j) + bets22p
            s2mp=s21(j) - bets22p
         else
            s1pp=s11(j)
            s1mp=s12(j)
            s2pp=s21(j)
            s2mp=s22(j)
         endif
         t1p=q2_u*s2pp + q2_up*s2p - q1_v*s1pp - q1_vp*s1p
         t2p=r4_p*s1pp + r4_pp*s1p - r3_q*s2pp - r3_qp*s2p
         c1p=q2_u*s2mp + q2_up*s2m - q1_v*s1mp - q1_vp*s1m
         c2p=r4_p*s1mp + r4_pp*s1m - r3_q*s2mp - r3_qp*s2m
c
         denp=(c2*f1p + c2p*f1 + c1*f2p + c1p*f2)
         deninvp=-denp*deninv*deninv
         num1p=-(c2*m1p + c2p*m1 + c1*m2p + c1p*m2)
         Vmat(j,1)=num1*deninvp + num1p*deninv
         gamLp=gami1(1)*Lp + gami1(j)*L
         Vmat(j,2)=gamL*deninvp + gamLp*deninv
         num2p=-(t1*f2p + t1p*f2 + t2*f1p + t2p*f1)
         Vmat(j,3)=num2*deninvp + num2p*deninv
         betLp=beti1(1)*Lp + beti1(j)*L
         Vmat(j,4)=-betL*deninvp - betLp*deninv
      enddo
c
      return
      end
