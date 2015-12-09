      subroutine rx_mode_fun(jm,phiz,dphiz,jlay_ref,ii_ref,vg)
c
c: Finds mode functions for the current eigenvalue.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jlay_ref,ii_ref,jlay,j,ii1,ii2,jz,jx0,jx,jm,jsr,jzx
      real*8 gamma,igamma,gamsq,
     .   gamz,cosgamz,singamz,phi_exp1,phi_exp2,A1,B1,A1_exp,
     .   B1_exp,A1tot,B1tot,xi0,xi1,ai1,bi1,aip1,bip1,zzexp1(3),
     .   sg,sg_rho,expmin,rzero,phsum,
     .   dphsum,expfaca,expfacb,expA,expB
      real*8 vg,phiq,dphiq
      complex*8 phiz(nzsr,jm),dphiz(nzsr,jm)
      complex*16 zzero
c: Set expmin to ln of min value of phi,dphi before being set to zero 
c: (expmin = -23 for phi = 1.d-10 here):
      data expmin/-23.d0/,zzero/(0.d0,0.d0)/,rzero/0.d0/
c
c: Find mode group velocity and attenuation:
      call rx_vg_atten
c
c: Find mode function values at desired depths:
      do jlay=1,nlwz
         j=jlwz(jlay)
         jx0=jzmx(j)-1
         if(j .ge. jlay_ref+ii_ref-1) then
c: Use top of layer to compute A1,B1:
            ii1=1
            ii2=2
         else
c: Use bottom of layer to compute A1,B1:
            ii1=2
            ii2=1
         endif
         phi_exp1=Aplay(ii1,j)
         phi_exp2=Aplay(ii2,j)
c: Check if mode function is sufficiently evanescent to set to zero:
         if(phi_exp1 .lt. expmin .and. phi_exp2 .lt. expmin) then
            do jz=1,nzmx(j)
               jzx=jx0 + jz
               phix(jzx)=zzero
               dphix(jzx)=zzero
            enddo
         elseif(iabs(isp(j)) .eq. 1) then
            gamsq=dreal(xksq(1,j)) - dreal(xkhsq)
            if(gamsq .gt. 0.d0) then
c: Isospeed propagating: phi=A cos(gamma*z) + B*sin(gamma*z) and set z=0
c: (zmx=0 at top of layer in rx_zmx_init):
               gamma=sqrt(gamsq)
               A1_exp=dexp(phi_exp1)
               A1=dreal(philay(1,j))*A1_exp
               B1=dreal(dphilay(1,j))*A1_exp/gamma
               do jz=1,nzmx(j)
                  jzx=jx0 + jz
                  gamz=gamma*zmx(jzx)
                  cosgamz=cos(gamz)
                  singamz=sin(gamz)
                  phsum=A1*cosgamz + B1*singamz
                  dphsum=-A1*singamz + B1*cosgamz
                  phix(jzx)=phsum
                  dphix(jzx)=gamma*dphsum
               enddo
            elseif(j .ne. nlay) then
c: Isospeed evanescent (see p. 41, ORCA II):
               igamma=sqrt(-gamsq)
               call rx_a1b1_iso(ii1,ii2,philay(1,j),dphilay(1,j),
     .            Aplay(1,j),zetalay(1,1,j),igamma,A1,B1,expA,expB)
               do jz=1,nzmx(j)
                  jzx=jx0 + jz
                  gamz=igamma*zmx(jzx)
                  expfaca=dexp(gamz + expA)
                  expfacb=dexp(-gamz + expB)
                  A1_exp=A1*expfaca
                  B1_exp=B1*expfacb
                  phix(jzx)=A1_exp + B1_exp
                  dphix(jzx)=igamma*(A1_exp - B1_exp)
               enddo
            else
c: Isospeed evanescent for halfspace, where A1=0 (see p. 41, ORCA II):
               igamma=sqrt(-gamsq)
               B1=dreal(philay(1,j))
               do jz=1,nzmx(j)
                  jzx=jx0 + jz
                  gamz=igamma*zmx(jzx)
                  expfacb=dexp(-gamz + phi_exp1)
                  phix(jzx)=B1*expfacb
                  dphix(jzx)=-igamma*phix(jzx)
               enddo
            endif
         else
c: Airy layer (see p. 42, ORCA II):
            xi0=(xkhsq - xksq(1,j))/etasq(j)
c
co            call rx_a1b1_airy(ii1,ii2,philay(1,j),dphilay(1,j),
co     .         Aplay(1,j),eta(j),ailay(1,1,1,j),bilay(1,1,1,j),
co     .         zetalay(1,1,j),A1,B1,expA,expB,pie)
            call rx_a1b1_airy(ii1,ii2,philay(1,j),dphilay(1,j),
     .         Aplay(1,j),dreal(eta(j)),ailay(1,1,1,j),bilay(1,1,1,j),
     .         zetalay(1,1,j),A1,B1,expA,expB,pie)
c
            do jz=1,nzmx(j)
               jzx=jx0 + jz
               xi1=xi0 - eta(j)*zmx(jzx)
               call rx_airy_sm(xi1,ai1,bi1,aip1,bip1,zzexp1)
               expfaca=dexp(expA - zzexp1(1))
               expfacb=dexp(expB + zzexp1(1))
               A1tot=A1*expfaca
               B1tot=B1*expfacb
               phsum=A1tot*aip1 + B1tot*bip1
               phix(jzx)=A1tot*ai1 + B1tot*bi1
               dphix(jzx)=-eta(j)*phsum
            enddo
         endif
      enddo
c
c: Make dphi at surface positive as a convention:
      sg=1.d0
      if(dreal(dphilay(1,2)) .lt. 0.d0) sg=-1.d0
c: Convert from phix to phi using mx_m pointers computed in zsr_init:
      do jsr=1,nzsr
         jx=mx_m(jsr)
         sg_rho=sg*rho_sr(jsr)
         phiq=sg_rho*dreal(phix(jx))
         dphiq=sg_rho*dreal(dphix(jx))
c: INCLUDE RHO in PHI,DPHI,PSI,DPSI:
         phiz(jsr,jm)=sg_rho*phix(jx)
         dphiz(jsr,jm)=sg_rho*dphix(jx)
      enddo
c      print *,nzsr,mx_m(1:nzsr),phiz(:,jm)
c
      return
      end
ccc
      subroutine rx_vg_atten
c
c: Finds modal group velocity and attenuation.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jlay,j,ii
      real*8 phi_duct,dphi_duct,gamsq_duct,mfun_mag,mfun_max,Apfac,
     .   N_n,vg,dn,phi_t,phi_b,dphi_t,dphi_b,I1,I2
c
c: Find duct in which mode has its maximum value:
      mfun_max=-1.d0
      nzref(nmode)=1
      do j=1,ndrx
         jlay=jval(1,j)
         ii=jval(2,j)
         phi_duct=philay(ii,jlay)
         dphi_duct=dphilay(ii,jlay)
         gamsq_duct=max(0.d0,dreal(xksq(ii,jlay)) - dreal(xkhsq))
         Apfac=dexp(Aplay(ii,jlay))
         mfun_mag=Apfac*(gamsq_duct*phi_duct*phi_duct + 
     .      dphi_duct*dphi_duct)
         if(mfun_mag .gt. mfun_max) then
            nzref(nmode)=j
            mfun_max=mfun_mag
         endif
      enddo
c
c: Compute mode attenuation and group velocity (see ORCA II, p.50, Lev etal):
      N_n=0.d0
      vg=0.d0
      dn=0.d0
cze   nzero=0
cii   do j=1,ndrx
cii      I2int(j)=0.d0
cii   enddo
cp    print *,'nmode,kn = ',nmode,kn(nmode)
      do j=2,nlay
         phi_t=philay(1,j)
         phi_b=philay(2,j)
         dphi_t=dphilay(1,j)
         dphi_b=dphilay(2,j)
co         call rx_int(phi_t,phi_b,dphi_t,dphi_b,Aplay(1,j),Aplay(2,j),
co     .      xilay(1,j),xilay(2,j),geo(1,3,j),geo(2,3,j),eta(j),etasq(j),
co     .      h(j),isp(j),gm(j),I1,I2)
         call rx_int(phi_t,phi_b,dphi_t,dphi_b,Aplay(1,j),Aplay(2,j),
     .      xilay(1,j),xilay(2,j),geo(1,3,j),geo(2,3,j),
     .      dreal(eta(j)),dreal(etasq(j)),
     .      h(j),isp(j),gm(j),I1,I2)
         N_n=N_n + I1
         vg=vg + (am(j)*I1 + bm(j)*I2)
         dn=dn + (betm(j)*I1 + gm(j)*I2)
cp    if(betm(j) .ne. 0.d0) print *,'j = ',j,betm(j),
cp   .   I1*rhom(j),gm(j),I2*rhom(j),rhom(j)
cii      I2int(jl2jd(j))=I2int(jl2jd(j)) + I2
c
c: Add up number of zeros of mode function:
cze      nzz=int(dimag(zetlay(1,1,j))/pie)
cze      if(mod(nzz,2) .eq. 0) then
cze         if((dreal(phi_t) .gt. 0.d0 .and. dreal(phi_b) .le. 0.d0) 
cze  .         .or. (dreal(phi_t) .lt. 0.d0 .and. 
cze  .          dreal(phi_b) .ge. 0.d0)) nzz=nzz+1
cze      else
cze         if((dreal(phi_t) .gt. 0. .and. dreal(phi_b) .ge. 0.) 
cze  .         .or. (dreal(phi_t) .lt. 0. .and. 
cze  .         dreal(phi_b) .le. 0.)) nzz=nzz+1
cze      endif
cze      nzero=nzero + nzz
      enddo
cii   I2=0.d0
cii   do j=1,ndrx
cii      if(I2int(j) .gt. I2) then
cii         nzref(nmode)=j
cii         I2=I2int(j)
cii      endif
cii   enddo
cc    vg=kn(nmode)/(w*(am(nlay)*N_inf/N_n + vg))
cc    dn=(w/kn(nmode))*(betm(nlay)*N_inf/N_n + dn)
      if(vg .eq. 0.d0 .or. w .eq. 0.d0) then
         print *,'vg bad'
      endif
      vg=kn(nmode)/(w*vg)
      dn=(w/kn(nmode))*dn
cze   mode_phz(2,nmode)=nzero
cc    print *,'j,N_n,vg,dn = ',nmode,N_n,vg,dn
ccc   print *,'nmode,dn = ',nmode,dn
c: Enter attenuation as imaginary part of kn:
      kn(nmode)=dcmplx(dreal(kn(nmode)),dn)
c
      return
      end
ccc
      subroutine rx_int(phi1,phi2,dphi1,dphi2,Ap1,Ap2,xi1,xi2,
     .   rho1,rho2,eta,etasq,h,iso,gm,I1,I2)
c
      implicit none
      integer*4 iso
      real*8 phi1,phi2,dphi1,dphi2,Ap1,Ap2,xi1,xi2,rho1,rho2,eta,
     .   etasq,h,gm,expmin,I1,I2,hsq2,term1,term2,term3,term4,
     .   dzero,ex1,ex2,phi1sq,phi2sq,dphi1sq,dphi2sq,mix2,done,
     .   rho_avg
      data dzero/0.d0/,expmin/-15.d0/,done/1.d0/
c
c: 1-16-96: Multiply phi1^2, etc by rho1, phi2^2 etc by rho2, rather than
c: whole integrals by rho1 or rho_avg.  This helps make attenuations more
c: accurate for layers with density gradients.
c
      I1=dzero
      I2=dzero
      ex1=done
      ex2=done
c: 7-8-99: Use of different rhos (rho1 at surface, rho2 at bottom) can
c: result in *negative* attenuations for higher-order modes near cutoff.
c: Go back to rho_avg to avoid this calamity:
      rho_avg=.5*(rho1 + rho2)
      if(iso .eq. 0) then
cc       I1=((dphi2*dphi2 - dphi1*dphi1)/eta - 
cc   .      (xi2*phi2*phi2 - xi1*phi1*phi1))/etasq
cc       I2=(3.d0*h + 2.d0*(xi2*dphi2*dphi2 - xi1*dphi1*dphi1)/eta - 
cc   .      (phi2*dphi2 - phi1*dphi1) - 3.d0*etasq*xi2*phi2sq*h -
cc   .      2.d0*eta*(xi2*xi2*phi2sq - xi1*xi1*phi1*phi1))/
cc   .      (3.d0*eta*etasq)
         if(Ap2 .gt. expmin .or. abs(phi2) .gt. 10.d0) then
            if(Ap2 .ne. 0.d0) ex2=dexp(2.d0*Ap2)
            phi2sq=phi2*phi2
            dphi2sq=dphi2*dphi2
            term2=rho_avg*(dphi2sq/etasq - xi2*phi2sq)
            term4=rho_avg*(dphi2sq*(3.d0*h + 2.d0*xi2/eta) -
     .         phi2*dphi2 - phi2sq*xi2*(3.d0*etasq*h + 
     .         2.d0*eta*xi2))
            if(Ap1 .gt. expmin .or. abs(phi1) .gt. 10.d0) then
               if(Ap1 .ne. 0.d0) ex1=dexp(2.d0*Ap1)
               phi1sq=phi1*phi1
               dphi1sq=dphi1*dphi1
               term1=rho_avg*(dphi1sq/etasq - xi1*phi1sq)
               I1=(term2*ex2 - term1*ex1)/eta
               term3=rho_avg*(2.d0*xi1*dphi1sq/eta -
     .            phi1*dphi1 - 2.d0*eta*xi1*xi1*phi1sq)
               I2=(term4*ex2 - term3*ex1)/(3.d0*eta*etasq)
            else
               I1=term2*ex2/eta
               I2=term4*ex2/(3.d0*eta*etasq)
            endif
         elseif(Ap1 .gt. expmin .or. abs(phi1) .gt. 10.d0) then
            if(Ap1 .ne. 0.d0) ex1=dexp(2.d0*Ap1)
            phi1sq=phi1*phi1
            dphi1sq=dphi1*dphi1
            term1=rho_avg*(dphi1sq/etasq - xi1*phi1sq)
            I1=-term1*ex1/eta
            term3=rho_avg*(2.d0*xi1*dphi1sq/eta -
     .         phi1*dphi1 - 2.d0*eta*xi1*xi1*phi1sq)
            I2=-term3*ex1/(3.d0*eta*etasq)
         endif
      else
         if(Ap2 .gt. expmin .or. abs(phi2) .gt. 10.d0) then
            if(Ap2 .ne. 0.d0) ex2=dexp(2.d0*Ap2)
            phi2sq=phi2*phi2
            dphi2sq=dphi2*dphi2
            mix2=phi2*dphi2
            term2=rho_avg*(phi2sq*h + (dphi2sq*h - mix2)/xi1)
            hsq2=2.d0*h*h
            term4=rho_avg*(dphi2sq*(hsq2 - 1.d0/xi1) + 
     .         phi2sq*(hsq2*xi1 + 1.d0) - 4.d0*mix2*h)
            if(Ap1 .gt. expmin .or. abs(phi1) .gt. 10.d0) then
               if(Ap1 .ne. dzero) ex1=dexp(2.d0*Ap1)
               term1=rho_avg*(-phi1*dphi1/xi1)
               I1=0.5d0*(term2*ex2 - term1*ex1)
               term3=rho_avg*(-dphi1*dphi1/xi1 + phi1*phi1)
               I2=(term4*ex2 - term3*ex1)/(8.d0*xi1)
            else
               I1=0.5d0*term2*ex2
               I2=term4*ex2/(8.d0*xi1)
            endif
         elseif(Ap1 .gt. expmin .or. abs(phi1) .gt. 10.d0) then
            if(Ap1 .ne. dzero) ex1=dexp(2.d0*Ap1)
            term1=rho_avg*(-phi1*dphi1/xi1)
            I1=-0.5d0*term1*ex1
            term3=rho_avg*(-dphi1*dphi1/xi1 + phi1*phi1)
            I2=-term3*ex1/(8.d0*xi1)
         endif
cc       I1=0.5d0*(phi2*phi2*h + (dphi2*dphi2*h - 
cc   .      (phi2*dphi2 - phi1*dphi1))/xi1)
cc       hsq2=2.d0*h*h
cc       I2=(dphi2*dphi2*(hsq2 - 1.d0/xi1) + dphi1*dphi1/xi1 +
cc   .      phi2*phi2*(hsq2*xi1 + 1.d0) - phi1*phi1 -
cc   .      4.d0*phi2*dphi2*h)/
      endif
c
      return
      end
ccc
      subroutine rx_a1b1_iso(ii1,ii2,philay,dphilay,Aplay,zetalay,
     .   igamma,A1,B1,expA,expB)
c
      implicit none
      integer*4 ii1,ii2,ii
      complex*16 philay(2),dphilay(2),zetalay(2)
      real*8 Aplay(2),igamma,A1,B1,expA,expB,term1,term2,
     .   rx_round,rata,ratb
c
      ii=ii1
      term1=dreal(philay(ii))
      term2=dreal(dphilay(ii))/igamma
      A1=rx_round(term1,term2,0.5d0,rata)
cc    if(rata .lt. 1.d-3) then
cc       term1h=dreal(philay(ii2))
cc       term2h=dreal(dphilay(ii2))/igamma
cc       A1x=rx_round(term1h,term2h,0.5d0,ratax)
cc       if(ratax .gt. rata) then
cc          A1=A1x
cc          ii=ii2
cc       endif
cc    endif
      expA=Aplay(ii) - dreal(zetalay(ii))
c
      ii=ii2
      term1=dreal(philay(ii))
      term2=dreal(dphilay(ii))/igamma
      B1=rx_round(term1,-term2,0.5d0,ratb)
cc    if(ratb .lt. 1.d-3) then
cc       term1h=dreal(philay(ii1))
cc       term2h=-dreal(dphilay(ii1))/igamma
cc       B1x=rx_round(term1h,term2h,0.5d0,ratbx)
cc       if(ratbx .gt. ratb) then
cc          B1=B1x
cc          ii=ii1
cc       endif
cc    endif
      expB=Aplay(ii) + dreal(zetalay(ii))
c
      return
      end
ccc
      subroutine rx_a1b1_airyx(ii1,ii2,philay,dphilay,Aplay,
     .   eta,ailay,bilay,zetalay,A1,B1,expA,expB,pie)
c
      implicit none
      integer*4 ii1,ii2,ii
      complex*16 philay(2),dphilay(2),ailay(2,2),bilay(2,2),zetalay(2)
      real*8 Aplay(2),eta,A1,B1,expA,expB,pie,term1,term2,
     .   rx_round,dphi_eta,rata,ratb
c
c: See ORCA II, p.65:
      ii=ii1
      dphi_eta=dreal(dphilay(ii1))/eta
      term1=bilay(2,ii1)*dreal(philay(ii1))
      term2=bilay(1,ii1)*dphi_eta
      A1=rx_round(term1,term2,pie,rata)
cc    if(rata .lt. 1.d-3) then
cc       dphi_etah=dreal(dphilay(ii2))/eta
cc       term1h=bilay(2,ii2)*dreal(philay(ii2))
cc       term2h=bilay(1,ii2)*dphi_etah
cc       A1x=rx_round(term1h,term2h,pie,ratax)
cc       if(ratax .gt. rata) then
cc          A1=A1x
cc          ii=ii2
cc       endif
cc    endif
      expA=Aplay(ii) + dreal(zetalay(ii))
c
      ii=ii2
      term1=ailay(2,ii)*dreal(philay(ii))
      term2=ailay(1,ii)*dphi_eta
      B1=rx_round(term1,term2,-pie,ratb)
cc    if(ratb .lt. 1.d-3) then
cc       dphi_etah=dreal(dphilay(ii1))/eta
cc       term1h=ailay(2,ii1)*dreal(philay(ii1))
cc       term2h=ailay(1,ii1)*dphi_etah
cc       B1x=rx_round(term1h,term2h,-pie,ratbx)
cc       if(ratbx .gt. ratb) then
cc          B1=B1x
cc          ii=ii1
cc       endif
cc    endif
      expB=Aplay(ii) - dreal(zetalay(ii))
c
      print *,'Airy phi1 check: ',philay(1)*exp(Aplay(1)),
     .A1*ailay(1,1)*dexp(expA-dreal(zetalay(1)))+
     .   B1*bilay(1,1)*dexp(expB+dreal(zetalay(1)))
      print *,'Airy phi2 check: ',philay(2)*exp(Aplay(2)),
     .A1*ailay(1,2)*dexp(expA-dreal(zetalay(2)))+
     .   B1*bilay(1,2)*dexp(expB+dreal(zetalay(2)))
      print *,'Airy dphi1 check: ',dphilay(1)*exp(Aplay(1)),
     .(A1*ailay(2,1)*dexp(expA-dreal(zetalay(1)))+
     .   B1*bilay(2,1)*dexp(expB+dreal(zetalay(1))))*(-eta)
      print *,'Airy dphi2 check: ',dphilay(2)*exp(Aplay(2)),
     .(A1*ailay(2,2)*dexp(expA-dreal(zetalay(2)))+
     .   B1*bilay(2,2)*dexp(expB+dreal(zetalay(2))))*(-eta)
c
      return
      end
ccc
      subroutine rx_a1b1_airy(ii1,ii2,philay,dphilay,Aplay,
     .   eta,ailay,bilay,zetalay,A1,B1,expA,expB,pie)
c
      implicit none
      integer*4 ii1,ii2,ii
      complex*16 philay(2),dphilay(2),ailay(2,2),bilay(2,2),zetalay(2)
      real*8 Aplay(2),eta,A1,B1,expA,expB,pie,term1,term2,
     .   rx_round,dphi_eta,rata,ratb
c
c: See ORCA II, p.65:
      ii=ii1
      dphi_eta=dreal(dphilay(ii))/eta
      term1=bilay(2,ii)*dreal(philay(ii))
      term2=bilay(1,ii)*dphi_eta
      A1=rx_round(term1,term2,pie,rata)
cc    if(rata .lt. 1.d-3) then
cc       dphi_etah=dreal(dphilay(ii2))/eta
cc       term1h=bilay(2,ii2)*dreal(philay(ii2))
cc       term2h=bilay(1,ii2)*dphi_etah
cc       A1x=rx_round(term1h,term2h,pie,ratax)
cc       if(ratax .gt. rata) then
cc          A1=A1x
cc          ii=ii2
cc       endif
cc    endif
      expA=Aplay(ii) + dreal(zetalay(ii))
c
      ii=ii2
      dphi_eta=dreal(dphilay(ii))/eta
      term1=ailay(2,ii)*dreal(philay(ii))
      term2=ailay(1,ii)*dphi_eta
      B1=rx_round(term1,term2,-pie,ratb)
cc    if(ratb .lt. 1.d-3) then
cc       dphi_etah=dreal(dphilay(ii1))/eta
cc       term1h=ailay(2,ii1)*dreal(philay(ii1))
cc       term2h=ailay(1,ii1)*dphi_etah
cc       B1x=rx_round(term1h,term2h,-pie,ratbx)
cc       if(ratbx .gt. ratb) then
cc          B1=B1x
cc          ii=ii1
cc       endif
cc    endif
      expB=Aplay(ii) - dreal(zetalay(ii))
c
cc    print *,'Airy phi1 check: ',philay(1)*exp(Aplay(1)),
cc   .A1*ailay(1,1)*dexp(expA-dreal(zetalay(1)))+
cc   .   B1*bilay(1,1)*dexp(expB+dreal(zetalay(1)))
cc    print *,'Airy phi2 check: ',philay(2)*exp(Aplay(2)),
cc   .A1*ailay(1,2)*dexp(expA-dreal(zetalay(2)))+
cc   .   B1*bilay(1,2)*dexp(expB+dreal(zetalay(2)))
cc    print *,'Airy dphi1 check: ',dphilay(1)*exp(Aplay(1)),
cc   .(A1*ailay(2,1)*dexp(expA-dreal(zetalay(1)))+
cc   .   B1*bilay(2,1)*dexp(expB+dreal(zetalay(1))))*(-eta)
cc    print *,'Airy dphi2 check: ',dphilay(2)*exp(Aplay(2)),
cc   .(A1*ailay(2,2)*dexp(expA-dreal(zetalay(2)))+
cc   .   B1*bilay(2,2)*dexp(expB+dreal(zetalay(2))))*(-eta)
c
      return
      end
ccc
      function rx_round(term1,term2,deninv,ratio)
c
c: Checks to see if cancellation occurs in term1+term2.
c: Returns round_check=(term1 + term2)*deninv if not.
      implicit none
      real*8 rx_round,term1,term2,deninv,numerator,magsum,ratio
c
      numerator=term1 + term2
      magsum=max(dabs(term1),dabs(term2))
      ratio=dabs(numerator)/magsum
cc    if(ratio .gt. 1.d-8) then
c: EKW FIX 5-11-90: bad mode functions when ratio = 3e-7 (see ORCA II, p77):
cc    if(ratio .gt. 1.d-8) then
      if(ratio .gt. 1.d-5) then
         rx_round=numerator*deninv
      else
         rx_round=0.d0
cc       print *,'Round to zero: ',term1,term2,numerator,deninv
      endif
c
      return
      end
