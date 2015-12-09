       subroutine rp_nomm(Rmat,e11,e12,e21,e22,ze,s11,s12,s21,s22,
     .   zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,Wmat,iiw)
c
c: Same as rp_slay, except the interface does not exist.  Formulas
c: derived from rp_slay by setting P=Q=1, U=V=0.
c
      implicit complex*16(a-z)
      integer*4 isoe,isos,ndv,j,iiw
      complex*16 e11(3),e12(3),e21(3),e22(3),s11(3),s12(3),s21(3),
     .   s22(3),Rmat(3,5),Vmat(3,5),Wmat(6)
      complex*16 gami1(3),beti1(3),gami2(3),beti2(3),zzero
      data zzero/(0.d0,0.d0)/
c
      if(Rmat(1,4) .eq. zzero .or. Rmat(1,2) .eq. zzero) then
         call rp_nomm0(Rmat,e11,e12,e21,e22,ze,s11,s12,s21,s22,zs,
     .      gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,Wmat,iiw)
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
      delta1=-2.d0*r2
      delta2=2.d0*q4
      delta3=q4sp*r1sp - q1*r4
      delta4=q3sp*r2sp - q2*r3
c
      x1=delta4*q4sp
      x2=delta3*q3sp
      x3=-delta4*q1
      x4=-delta3*q2
c
      w1=delta2*q2
      w2=delta2*q1
      w3=delta1*q3sp
      w4=delta1*q4sp
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
         a7=e2m*q1 - e1m*q2
         a8=e2m*r1sp - e1m*r2sp
         b7=s1m*q4sp - s2m*q3sp
         b8=s1m*r4 - s2m*r3
         wden=a7*b8 - b7*a8
         fac=-2.d0/wden
         a9_den=fac*gami1(1)
         Wmat(1)=a9_den*b8
         Wmat(2)=-a9_den*b7
         Wmat(5)=-ze
         a9bar_den=fac*beti1(1)
         Wmat(3)=-a9bar_den*a7
         Wmat(4)=a9bar_den*a8
         Wmat(6)=-zs
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
         x1p=delta4*q4spp + delta4p*q4sp
         x2p=delta3*q3spp + delta3p*q3sp
         x3p=-delta4*q1p - delta4p*q1
         x4p=-delta3*q2p - delta3p*q2
c
         w1p=delta2*q2p + delta2p*q2
         w2p=delta2*q1p + delta2p*q1
         w3p=delta1*q3spp + delta1p*q3sp
         w4p=delta1*q4spp + delta1p*q4sp
c
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
