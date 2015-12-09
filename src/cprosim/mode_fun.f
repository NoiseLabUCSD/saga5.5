      subroutine mode_fun(kmode,R1x,R2x,dW_dk,phiz,dphiz,psiz,dpsiz,
     .   expz_gbs,jm)
c
c: Finds the mode function whose eigenvalue is kn as a function of depth.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 jm,ii,jflu1,jflu2,inc,jsol1,jsol2,jx1,j,jx,jsr,
     .   ii1,ii2,jjhsp,joff,k,j1
      complex*8 phiz(nzsr,jm),dphiz(nzsr,jm),psiz(nzsr,jm),
     .   dpsiz(nzsr,jm)
      complex*16 kmode,R1x,R2x,dW_dk
      complex*16 phiref,dphiref,phi0,dphi0,psi0,dpsi0,mnorm(2),arg,
     .   Ap,As,phi2,dphi2,psi2,dpsi2
      real*8 expz_gbs(nzsr,jm),sg,Apexp,Asexp
c
c: Include minus sign here so that mode functions turn out real:
c     mnorm(1)=1.d0/cdsqrt(-dW_dk)
c: Include sqrt(k) in normalization & keep track of sheets for sqrt operation:
      arg=-2.d0*kmode/(dW_dk*rho_duct)
c: (-1) factor now determined by dphi at top (see sg):
      mnorm(1)=cdsqrt(arg)
      mnorm(2)=R1x*mnorm(1)
      phiref=mnorm(1)*(1.d0 + R1x)
      dphiref=gamiref*mnorm(1)*(1.d0 - R1x)
      do ii=1,2
         ii1=ii
         ii2=3 - ii
         jflu1=jflu(ii,1)
         jsol1=jsol(ii,1)
         inc=jflu(ii,3)
         joff=jflu(ii,4) - inc
         jjhsp=jhsp(ii)
         if(allf(ii) .eq. 1) then
c: If all fluid, include halfspace by adding inc:
            jflu2=jflu(ii,2) + inc
            jsol2=jsol(ii,2)
         else
            jflu2=jflu(ii,2)
c: If any solid layers, include halfspace by adding inc:
            jsol2=jsol(ii,2)+inc
         endif
c: Limit counters to layers in which receivers are located:
c: Don't limit counters now (for mode id and -1 norm fac):
cxx      if(ii .eq. 1) then
cxx         jflu2=min(jflu2,jlmax)
cxx         jsol2=min(jsol2,jlmax)
cxx      else
cxx         jflu2=max(jflu2,jlmin)
cxx         jsol2=max(jsol2,jlmin)
cxx      endif
c
c: Initialize phi,dphi at reference depth:
         phi0=phiref
c: Make derivative continuous by multiplying dphi0 by -1 for upgoing 
c: reflection coefficient (see p. 104-105):
         dphi0=inc*dphiref
c: Amplitude of plane wave downward is 1, upward is R1x [both norm by 
c: 1/sqrt(W')]
         Ap=mnorm(ii)
         Apexp=0.d0
c: Do ref depth layer if we are doing bottom half and ref depth
c: at top of layer or if we are doing top half and ref depth is at bottom:
         if(isvmin .eq. ii) then
            j=jflu1
            jx1=jzmx(j)
c: FIX 2-6-95: Give em_calc phi,dphi at bottom of current layer:
            call phi_prop(j,ii1,ii2,inc,phi0,dphi0,Ap,Apexp)
            call em_calc(xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .etasq(j),gami(1,ii1,j),isp(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .philay(1,j),dphilay(1,j),Aplay(1,j),phix(jx1),dphix(jx1),
     .expx_gbs(jx1),xi,ai,aip,bi,bip,zzexp,inc,ailay(1,1,1,j),
     .bilay(1,1,1,j),zetalay(1,1,j),aisoln(1,j),ii1,ii2,j,jjhsp,kw0)
         endif
c
         do j=jflu1+inc,jflu2,inc
c: Propagate phi across fluid-fluid interface:
            phi0=rhorat(j+joff)*phi0
            jx1=jzmx(j)
            if(j .eq. jjhsp) then
               philay(ii1,j)=phi0
               dphilay(ii1,j)=inc*dphi0
               Aplay(ii1,j)=Apexp
            else
c: FIX 2-6-95: Give em_calc phi,dphi at bottom of current layer:
c: Compute phi,dphi at bottom of fluid layer:
               call phi_prop(j,ii1,ii2,inc,phi0,dphi0,Ap,Apexp)
            endif
            call em_calc(xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .etasq(j),gami(1,ii1,j),isp(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .philay(1,j),dphilay(1,j),Aplay(1,j),phix(jx1),dphix(jx1),
     .expx_gbs(jx1),xi,ai,aip,bi,bip,zzexp,inc,ailay(1,1,1,j),
     .bilay(1,1,1,j),zetalay(1,1,j),aisoln(1,j),ii1,ii2,j,jjhsp,kw0)
         enddo
c
         As=dcmplx(0.d0,0.d0)
         Asexp=-1.d+200
         if(inc*(jsol2-jsol1) .ge. 0) then
            j=jsol1
c: Compute p-wave and s-wave potentials at top of first solid layer:
            call phipsi(Ap,As,Vmat(1,1,j,ii),Wmat(1,j,ii),phi2,dphi2,
     .         Apexp,psi2,dpsi2,Asexp,j,ii1,inc)
         endif
         do j=jsol1,jsol2,inc
            if(j .ne. jjhsp) then
c: For layers "above" halfspace:
               j1=j + inc
c: Compute p-wave and s-wave potentials at top of next layer:
               call phipsi(Ap,As,Vmat(1,1,j1,ii),Wmat(1,j1,ii),
     .            phi0,dphi0,Apexp,psi0,dpsi0,Asexp,j1,ii1,inc)
c: Save these quantities at top of next layer for when we do halfspace:
               phi2=phi0
               dphi2=dphi0
               psi2=psi0
               dpsi2=dpsi0
               k=j1 + joff
c: Convert p- and s-wave potentials to bottom of current layer:
               call phi_xlay(phi0,dphi0,Apexp,psi0,dpsi0,Asexp,
     .            Pcon(1,k),Qcon(1,k),Ucon(1,k),Vcon(1,k),ikcon,
     .            rhorat(k),j,ii2,inc,mm(k))
            else
c: When doing halfspace, set phi0,dphi0,psi0,dpsi0 to values at top
c: of halfspace (note that Apexp will be unchanged and correct):
               phi0=phi2
               dphi0=dphi2
               psi0=psi2
               dpsi0=dpsi2
            endif
            jx1=jzmx(j)
c: Find p-wave potential in layer:
            call em_calc(xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .etasq(j),gami(1,ii1,j),isp(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .philay(1,j),dphilay(1,j),Aplay(1,j),phix(jx1),dphix(jx1),
     .expx_gbs(jx1),xi,ai,aip,bi,bip,zzexp,inc,ailay(1,1,1,j),
     .bilay(1,1,1,j),zetalay(1,1,j),aisoln(1,j),ii1,ii2,j,jjhsp,kw0)
c: Find s-wave potential in layer:
            call em_calc(xkhsq,xbsq(ii1,j),xbsq(ii2,j),inc*etb(j),
     .etbsq(j),beti(1,ii1,j),iss(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .psilay(1,j),dpsilay(1,j),Aslay(1,j),psix(jx1),dpsix(jx1),
     .expx_gbs(jx1),xis,ais,aips,bis,bips,zzexps,inc,ailay(1,1,2,j),
     .bilay(1,1,2,j),zetalay(1,2,j),aisoln(2,j),ii1,ii2,j,jjhsp,kw0)
         enddo
         if(ii .eq. 2) sg=-dsign(1.d0,dreal(dphi0))
      enddo
c
c: Convert from phix to phi using mx_m pointers computed in zsr_init:
      do jsr=1,nzsr
         jx=mx_m(jsr)
c: INCLUDE RHO in PHI,DPHI,PSI,DPSI:
c:pln sometimes phix is huge (no call to em_calc for this run)
c:pln Engineering solution
         if(cdabs(phix(jx)).gt.1.d10) then
            jjfail=1
            return
         end if
         phiz(jsr,jm)=sg*rho_sr(jsr)*phix(jx)
         dphiz(jsr,jm)=sg*rho_sr(jsr)*dphix(jx)
         psiz(jsr,jm)=sg*rho_sr(jsr)*psix(jx)
         dpsiz(jsr,jm)=sg*rho_sr(jsr)*dpsix(jx)
c: Exponential for complex source points:
         expz_gbs(jsr,jm)=expx_gbs(jx)
      enddo
cc    print *,'jm,phi_src,phi_rec = ',jm,cabs(phiz(mzsrc(1),jm)),
cc   .   cabs(phiz(mzsrc(1)-1,jm)),cabs(phiz(mzsrc(1)+1,jm))
c
      return
      end
ccc
      subroutine phi_prop(j,ii1,ii2,inc,phi0,dphi0,Ap,Apexp)
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 j,ii1,ii2,inc
      complex*16 Vmatx,Wmatx,Wexp,gami2,phi0,dphi0,Ap
      real*8 Apexp,Wphase
c
      philay(ii1,j)=phi0
      dphilay(ii1,j)=inc*dphi0
      Aplay(ii1,j)=Apexp
c: Compute phi0,dphi0 [to be multiplied by exp(Apexp)] at bottom of layer
c: from phi0,dphi0 at top of layer:
      Wexp=Wmat(5,j,ii1)
      Apexp=Apexp + dreal(Wexp)
      Wphase=dimag(Wexp)
c: Multiply plane wave amplitude by transmission coefficient [except for
c: the exp(Apexp) part]:
      Wmatx=Wmat(1,j,ii1)
      Ap=Ap*Wmatx*dcmplx(dcos(Wphase),dsin(Wphase))
      Vmatx=Vmat(1,1,j,ii1)
      gami2=gami(1,ii2,j)
      phi0=Ap*(1.d0 + Vmatx)
      dphi0=Ap*gami2*(1.d0 - Vmatx)
      philay(ii2,j)=phi0
      dphilay(ii2,j)=inc*dphi0
      Aplay(ii2,j)=Apexp
c
      return
      end
ccc
      subroutine phipsi(Ap,As,Vmatx,Wmatx,phi0,dphi0,Apexp,psi0,dpsi0,
     .   Asexp,j,ii1,inc)
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 j,ii1,inc
      complex*16 Ap,As,Vmatx(3,5),Wmatx(6),phi0,dphi0,psi0,
     .   dpsi0,vpp,Ap_vps,vss,As_vsp,cfacp,cfacs,cfac,
     .   term1p,term2p,term1s,term2s
      real*8 Apexp,Asexp,delex,spreal,Wphase,Vphase,ex1,ex2
c
c: Do imaginary part of cdexp(Wmatx(5)) term for p-wave transmission coeffs:
      Wphase=dimag(Wmatx(5))
      cfacp=Ap*dcmplx(dcos(Wphase),dsin(Wphase))
c: Do imaginary part of cdexp(Wmatx(6)) term for s-wave transmission coeffs:
      Wphase=dimag(Wmatx(6))
      cfacs=As*dcmplx(dcos(Wphase),dsin(Wphase))
c: Ap is to be multiplied by exp(Apexp), As by exp(Asexp):
cxx   Apexp=Apexp + dreal(Wmatx(5))
cxx   Asexp=Asexp + dreal(Wmatx(6))
cxx   delex=Apexp - Asexp
cxx   efac=dexp(delex)
c: Compute new amplitudes of downward-going p-wave and s-wave plane waves:
      term1p=Wmatx(1)*cfacp
      term2p=Wmatx(4)*cfacs
      ex1=Apexp + dreal(Wmatx(5))
      term1s=Wmatx(3)*cfacs
      term2s=Wmatx(2)*cfacp
      ex2=Asexp + dreal(Wmatx(6))
c: Make sure that Ap,As have magnitudes on the order of unity:
      call norm_add(term1p,ex1,term2p,ex2,Ap,Apexp)
      call norm_add(term1s,ex2,term2s,ex1,As,Asexp)
cxx   Ap=Apx*Wmatx(1)*cfacp + Asx*Wmatx(4)*cfacs/efac
cxx   As=Asx*Wmatx(3)*cfacs + Apx*Wmatx(2)*cfacp*efac
c: Compute p-wave and s-wave potentials at top of layer:
      Vphase=dimag(Vmatx(1,5))
      cfac=dcmplx(dcos(Vphase),dsin(Vphase))
      spreal=dreal(Vmatx(1,5))
      delex=Apexp - Asexp
      vpp=Vmatx(1,1)
c: Include exponential factors due to Vps and Ap,As:
      Ap_vps=Ap*Vmatx(1,2)*cfac*dexp(spreal + delex)
      vss=Vmatx(1,3)
c: Include exponential factors due to Vsp and Ap,As:
      As_vsp=As*Vmatx(1,4)*cfac*dexp(spreal - delex)
c: phi0,dphi0 must be multiplied by exp(Apexp):
      phi0=Ap*(1.d0 + vpp) + As_vsp
      dphi0=gami(1,ii1,j)*(Ap*(1.d0 - vpp) - As_vsp)
c: psi0,dpsi0 must be multiplied by exp(Asexp):
      psi0=As*(1.d0 + vss) + Ap_vps
      dpsi0=beti(1,ii1,j)*(As*(1.d0 - vss) - Ap_vps)
c
      philay(ii1,j)=phi0
      dphilay(ii1,j)=inc*dphi0
      Aplay(ii1,j)=Apexp
      psilay(ii1,j)=psi0
      dpsilay(ii1,j)=inc*dpsi0
      Aslay(ii1,j)=Asexp
c
      return
      end
ccc
      subroutine norm_add(term1,ex1,term2,ex2,A,Aexp)
c
c: Computes A*exp(Aexp) = term1*exp(ex1) + term2*exp(ex2), such
c: that the magnitude of A is on the order of unity.
c
      implicit none
      complex*16 term1,term2,A
      real*8 ex1,ex2,Aexp,mag1,mag2,exx1,exx2
c
      if(term1 .ne. (0.d0,0.d0)) then
         mag1=cdabs(term1)
         term1=term1/mag1
         exx1=ex1 + dlog(mag1)
      else
         exx1=-1.d200
      endif
      if(term2 .ne. (0.d0,0.d0)) then
         mag2=cdabs(term2)
         term2=term2/mag2
         exx2=ex2 + dlog(mag2)
      else
         exx2=-1.d200
      endif
      if(exx1 .gt. exx2) then
c: term1 dominates:
         A=term1 + term2*dexp(exx2 - exx1)
         Aexp=exx1
      else
c: term2 dominates:
         A=term1*dexp(exx1 - exx2) + term2
         Aexp=exx2
      endif
c
      return
      end
ccc
      subroutine phi_xlay(phi0,dphi0,Apexp,psi0,dpsi0,Asexp,Pconx,
     .   Qconx,Uconx,Vconx,ikconx,rhoratx,j,ii2,inc,mmx)
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 j,ii2,inc,mmx
      complex*16 phi0,dphi0,psi0,dpsi0,Pconx,Qconx,Uconx,Vconx,ikconx,
     .   phib,dphib,psib,dpsib,Pbar,Qbar,Ubar,Vbar
      real*8 Apexp,Asexp,rhoratx,efac
c
cc    print *,'uz below: ',dphi0*dexp(Apexp) + ikconx*psi0*dexp(Asexp)
c: Find P,Q,U,V for going backward across interface (see p. 121):
c: If no mismatch, phi and psi are continuous:
      if(mmx .eq. 1) then
         Pbar=Qconx/rhoratx
         Qbar=Pconx/rhoratx
         Ubar=Uconx/(-rhoratx)
         Vbar=ikconx*(1.d0 - Pbar)
c
c: Propagate Phi backward across solid-solid interface (see p. 87):
         efac=exp(Apexp - Asexp)
         phib=Pbar*phi0 + Ubar*dpsi0/efac
         dphib=Qbar*dphi0 + Vbar*psi0/efac
         psib=-Ubar*dphi0*efac + Pbar*psi0
         dpsib=-Vbar*phi0*efac + Qbar*dpsi0
c
         phi0=phib
         dphi0=dphib
         psi0=psib
         dpsi0=dpsib
      endif
cc    print *,'uz above: ',dphi0*dexp(Apexp) + ikconx*psi0*dexp(Asexp)
c
      philay(ii2,j)=phi0
      dphilay(ii2,j)=inc*dphi0
      Aplay(ii2,j)=Apexp
      psilay(ii2,j)=psi0
      dpsilay(ii2,j)=inc*dpsi0
      Aslay(ii2,j)=Asexp
c
      return
      end
