      subroutine rp_sfint(Rmat,Acon1,Bcon1,ikcon,rhorat,gami1,
     .   gami2,beti2,R,ndv,Wmat,iiw)
c
c: Finds the p-p reflection coefficient just above a fluid-solid interface
c: given p-p, p-s, s-s, and s-p reflection coefficients just below the
c: fluid-solid interface.  The interface is characterized by P,Q,U,V
c: gami2,beti2 are the p- and w-wave vertical wavenumbers in the solid,
c: and gami1 is the p-wave vertical wavenumber in the fluid.
c:
      implicit complex*16(a-z)
      integer*4 ndv,j,iiw,iiz,iiz0,iir
      complex*16 Rmat(3,5),Acon1(3),Bcon1(3),gami1(3),
     .   gami2(3),beti2(3),R(3),Wmat(6),ikcon(3),zzero
      real*8 rhorat,magsq,deps,rat
      data zzero/(0.d0,0.d0)/,deps/1.d-50/
c
      if(Rmat(1,1) .eq. zzero .and. Rmat(1,3) .eq. zzero) then
         iiz0=1
         q1=(1.d0,0.d0)
         q2=gami2(1)
         r3=(1.d0,0.d0)
         r4=beti2(1)
      else
         iiz0=0
         q1=1.d0 + Rmat(1,1)
         q2diff=1.d0 - Rmat(1,1)
         q2=gami2(1)*q2diff
         r3=1.d0 + Rmat(1,3)
         r4diff=1.d0 - Rmat(1,3)
         r4=beti2(1)*r4diff
      endif
c
      q2_rho=q2*rhorat
      ik_rho=ikcon(1)*rhorat
      a2gam=Acon1(1)*gami1(1)
      b2gam=Bcon1(1)*gami1(1)
      q1_a2gam=q1*a2gam
c
      c1=q2_rho + q1_a2gam
      d1=q2_rho - q1_a2gam
      f1=q2*Bcon1(1)
c: Case where vps and vsp below are zero:
cxx   if(Rmat(1,2) .eq. zzero .and. Rmat(1,4) .eq. zzero) then
      if(magsq(Rmat(1,2)) .lt. deps .or. 
     .   magsq(Rmat(1,4)) .lt. deps) then
         iiz=1
         delta3=-q1*r4
         delta4=-q2*r3
         d3=-q1*ik_rho
         f3=-q1*Acon1(1)
         d4=q2*b2gam
c
         j1=c1
         k1=d1
         m1=f1
      else
         iiz=0
c: rps and rsp are multiplied by exp(zsp) (see r1,r2,q3,q4):
         spfac=cdexp(Rmat(1,5))
         r1sp=Rmat(1,4)*spfac
         r2sp=-gami2(1)*r1sp
         q3sp=Rmat(1,2)*spfac
         q4sp=-beti2(1)*q3sp
c: exp(zsp) factor must be accounted for in delta3,delta4, etc:
         delta1=2.d0*gami2(1)*Rmat(1,4)
         delta2=-2.d0*beti2(1)*Rmat(1,2)
         delta3=q4sp*r1sp - q1*r4
         delta4=q3sp*r2sp - q2*r3
c
         q4_b2gam=q4sp*b2gam
         q3_ikrho=q3sp*ik_rho
         c2=q4_b2gam - q3_ikrho
         d2=q4_b2gam + q3_ikrho
         f2=q3sp*Acon1(1)
c
         d3=q4sp*rhorat - q1*ik_rho
         f3=q4sp*Bcon1(1) - q1*Acon1(1)
         d4=q3sp*a2gam + q2*b2gam
c
         j1=delta2*c1 - delta1*c2
         k1=delta2*d1 + delta1*d2
         m1=delta2*f1 + delta1*f2
      endif
c
      del4_d3=delta4*d3
      del3_d4=delta3*d4
      j2=del4_d3 + del3_d4
      k2=del4_d3 - del3_d4
      m2=delta4*f3
c
      den=m2*k1 - m1*k2
      num=m1*j2 - m2*j1
      R(1)=num/den
      if(iiw .eq. 1) then
         h2=r3*ik_rho + r4*b2gam
         if(iiz .eq. 1) then
            h4=r3*Acon1(1)
            h5=-h2
            g4=q2*Bcon1(1)
            g5=-d1
         else
            h1=r1sp*a2gam - r2sp*rhorat
            h4=r2sp*Bcon1(1) + r3*Acon1(1)
            h5=h1 - h2
            g4=q2*Bcon1(1) + q3sp*Acon1(1)
            g5=-(d1 + d2)
         endif
         L5=2.d0*gami1(1)*rhorat*(Acon1(1) - ikcon(1)*Bcon1(1))
c: OLD:
cxx      wden=h4*g5 - h5*g4
cxx      fac=L5/wden
cxx      Wmat(1)=h4*fac
cxx      Wmat(2)=-g4*fac
c: NEW:
         wden2=round_check(h4*g5,-h5*g4,1.d0/L5,iir,rat)
         if(iir .eq. 1) then
            Wmat(1)=zzero
            Wmat(2)=zzero
cxx         print *,'iir=1 in rp_sfint: '
         else
            Wmat(1)=h4/wden2
            Wmat(2)=-g4/wden2
         endif
         Wmat(3)=(0.d0,0.d0)
         Wmat(4)=(0.d0,0.d0)
         Wmat(5)=(0.d0,0.d0)
         Wmat(6)=(0.d0,0.d0)
      endif
c
c: Derivatives:
      do j=2,ndv
         if(iiz0 .eq. 1) then
            q1p=(0.d0,0.d0)
            q2p=gami2(j)
            r3p=(0.d0,0.d0)
            r4p=beti2(j)
         else
            q1p=Rmat(j,1)
            q2p=-gami2(1)*Rmat(j,1) + gami2(j)*q2diff
            r3p=Rmat(j,3)
            r4p=-beti2(1)*Rmat(j,3) + beti2(j)*r4diff
         endif
c
         q2_rhop=q2p*rhorat
         ik_rhop=ikcon(j)*rhorat
         a2gamp=Acon1(1)*gami1(j) + Acon1(j)*gami1(1)
         b2gamp=Bcon1(1)*gami1(j) + Bcon1(j)*gami1(1)
         q1_a2gamp=q1*a2gamp + q1p*a2gam
c
         c1p=q2_rhop + q1_a2gamp
         d1p=q2_rhop - q1_a2gamp
         f1p=q2*Bcon1(j) + q2p*Bcon1(1)
c
         if(iiz .eq. 1) then
            delta3p=-q1*r4p - q1p*r4
            delta4p=-q2*r3p - q2p*r3
            d3p=-q1*ik_rhop - q1p*ik_rho
            f3p=-q1*Acon1(j) - q1p*Acon1(1)
            d4p=q2*b2gamp + q2p*b2gam
c
            j1p=c1p
            k1p=d1p
            m1p=f1p
         else
            r1spp=Rmat(j,4)*spfac
            r2spp=(-gami2(1)*Rmat(j,4) - gami2(j)*Rmat(1,4))*spfac
            q3spp=Rmat(j,2)*spfac
            q4spp=(-beti2(1)*Rmat(j,2) - beti2(j)*Rmat(1,2))*spfac
c
            delta1p=2.d0*(gami2(1)*Rmat(j,4) + gami2(j)*Rmat(1,4))
            delta2p=-2.d0*(beti2(1)*Rmat(j,2) + beti2(j)*Rmat(1,2))
            delta3p=q4sp*r1spp + q4spp*r1sp - q1*r4p - q1p*r4
            delta4p=q3sp*r2spp + q3spp*r2sp - q2*r3p - q2p*r3
c
            q4_b2gamp=q4sp*b2gamp + q4spp*b2gam
            q3_ikrhop=q3sp*ik_rhop + q3spp*ik_rho
            c2p=q4_b2gamp - q3_ikrhop
            d2p=q4_b2gamp + q3_ikrhop
            f2p=q3sp*Acon1(j) + q3spp*Acon1(1)
c
            d3p=q4spp*rhorat - q1*ik_rhop - q1p*ik_rho
            f3p=q4sp*Bcon1(j) + q4spp*Bcon1(1) - q1*Acon1(j) - 
     .         q1p*Acon1(1)
            d4p=q3sp*a2gamp + q3spp*a2gam + q2*b2gamp + q2p*b2gam
c
            j1p=delta2*c1p + delta2p*c1 - delta1*c2p - delta1p*c2
            k1p=delta2*d1p + delta2p*d1 + delta1*d2p + delta1p*d2
            m1p=delta2*f1p + delta2p*f1 + delta1*f2p + delta1p*f2
         endif
c
         del4_d3p=delta4*d3p + delta4p*d3
         del3_d4p=delta3*d4p + delta3p*d4
         j2p=del4_d3p + del3_d4p
         k2p=del4_d3p - del3_d4p
         m2p=delta4*f3p + delta4p*f3
c
         denp=m2*k1p + m2p*k1 - m1*k2p - m1p*k2
         nump=m1*j2p + m1p*j2 - m2*j1p - m2p*j1
cxx      R(j)=(den*nump - denp*num)/(den*den)
         R(j)=(nump-denp*R(1))/den
      enddo
c
      return
      end
