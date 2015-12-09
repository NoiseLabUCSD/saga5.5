      subroutine rp_nomm0(Rmat,e11,e12,e21,e22,ze,s11,s12,s21,s22,
     .   zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,Wmat,iiw)
c
c: Same as rp_slay0, except no interface exists at bottom.  Formulas
c: derived from rp_slay0 by setting P=Q=1 and U=V=0.
c
      implicit complex*16 (a-z)
      integer*4 isoe,isos,ndv,j,iiw,iiz
      complex*16 Rmat(3,5),e11(3),e12(3),e21(3),e22(3),s11(3),
     .   s12(3),s21(3),s22(3),Vmat(3,5),gami1(3),gami2(3),beti1(3),
     .   beti2(3),Wmat(6),zzero
      data zzero/(0.d0,0.d0)/
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
c
      if(Rmat(1,1) .eq. zzero .and. Rmat(1,3) .eq. zzero) then
         iiz=1
         q2=gami2(1)
         r4=beti2(1)
         m1=q2*e1p - e2p
         f1=q2*e1m - e2m
         t2=r4*s1p - s2p
         c2=r4*s1m - s2m
      else
         iiz=0
         q1=1.d0 + Rmat(1,1)
         q2diff=1.d0 - Rmat(1,1)
         q2=gami2(1)*q2diff
         r3=1.d0 + Rmat(1,3)
         r4diff=1.d0 - Rmat(1,3)
         r4=beti2(1)*r4diff
         m1=q2*e1p - q1*e2p
         f1=q2*e1m - q1*e2m
         t2=r4*s1p - r3*s2p
         c2=r4*s1m - r3*s2m
      endif
c
      Vmat(1,1)=-m1/f1
      Vmat(1,2)=zzero
      Vmat(1,3)=-t2/c2
      Vmat(1,4)=zzero
c
c: vps(1:2) and vsp(1:2) are to be multiplied by exp(zexp_sp):
      Vmat(1,5)=zzero
c
      if(iiw .eq. 1) then
         Wmat(1)=2.d0*gami1(1)/f1
         Wmat(2)=zzero
         Wmat(5)=-ze
         Wmat(3)=2.d0*beti1(1)/c2
         Wmat(4)=zzero
         Wmat(6)=-zs
      endif
c
      do j=2,ndv
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
c
         if(iiz .eq. 1) then
            q2p=gami2(j)
            r4p=beti2(j)
            m1p=q2*e1pp + q2p*e1p - e2pp
            f1p=q2*e1mp + q2p*e1m - e2mp
            t2p=r4*s1pp + r4p*s1p - s2pp
            c2p=r4*s1mp + r4p*s1m - s2mp
         else
            q1p=Rmat(j,1)
            q2diffp=-Rmat(j,1)
            q2p=gami2(1)*q2diffp + gami2(j)*q2diff
            r3p=Rmat(j,3)
            r4diffp=-Rmat(j,3)
            r4p=beti2(1)*r4diffp + beti2(j)*r4diff
            m1p=q2*e1pp + q2p*e1p - q1*e2pp - q1p*e2p
            f1p=q2*e1mp + q2p*e1m - q1*e2mp - q1p*e2m
            t2p=r4*s1pp + r4p*s1p - r3*s2pp - r3p*s2p
            c2p=r4*s1mp + r4p*s1m - r3*s2mp - r3p*s2m
         endif
c
         Vmat(j,1)=-(m1p/f1 + Vmat(1,1)*f1p/f1)
         Vmat(j,2)=zzero
         Vmat(j,3)=-(t2p/c2 + Vmat(1,3)*c2p/c2)
         Vmat(j,4)=zzero
      enddo
c
      return
      end
