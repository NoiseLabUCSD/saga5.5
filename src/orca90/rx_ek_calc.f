      subroutine rx_ek_calc(k,e_k,phm,ndv)
c
c: Computes the plane wave reflection coefficient for a series of
c: fluid and solid layers.
c
      use parms_com
      use i_o_com
      use gen_com
c: Local variables:
      integer*4 ndv,j,jm1,jdv
      real*8 k,e_k(6),phm,dzero,g12_k,g22_k,g12_w,g22_w,
     .   g12_kk,g22_kk,capx_k,capx_w,capx_fac,gamsq,g12x,g22x
      data dzero/0.d0/
c
      nctot=nctot + 1
      xkh=k
      xkhsq=dcmplx(k*k,0.d0)
c
c: Phi(nlay) = G Phi(j), where values are taken at BOTTOMS of layers [except
c: for the halfspace layer (layer nlay)]:
c: Set values for G at bottom of layer above lower halfspace:
      jm1=nlay-1
      g11(1,2,jm1)=rhofac(jm1)
      g12(1,2,jm1)=dzero
      g21(1,2,jm1)=dzero
      g22(1,2,jm1)=1.d0
      g_exp(1,jm1)=dzero
      do j=2,ndv
         g11(j,2,jm1)=dzero
         g12(j,2,jm1)=dzero
         g21(j,2,jm1)=dzero
         g22(j,2,jm1)=dzero
         g_exp(j,jm1)=dzero
      enddo
      phtot=0.d0
      do j=nlay-1,2,-1
         jm1=j-1
         call rx_ep_calc(xkh,w,wsq,xkhsq,xksq(1,j),xksq(2,j),eta(j),
     .      etasq(j),h(j),u11(1,j),u12(1,j),u21(1,j),u22(1,j),
     .      u_exp(1,j),isp(j),ndv,ailay(1,1,1,j),bilay(1,1,1,j),
     .      xilay(1,j),zetalay(1,1,j),1,2,phtot)
55       continue
         call rx_2x2(g11(1,2,j),g12(1,2,j),g21(1,2,j),g22(1,2,j),
     .      u11(1,j),u12(1,j),u21(1,j),u22(1,j),
     .      g11(1,1,j),g12(1,1,j),g21(1,1,j),g22(1,1,j),ndv)
         do jdv=1,ndv
            g_exp(jdv,jm1)=g_exp(jdv,j) + u_exp(jdv,j)
         enddo
c: Take G across interface from top of layer j to bottom of layer jm1:
         call rx_iface(g11(1,1,j),g12(1,1,j),g21(1,1,j),g22(1,1,j),
     .      g11(1,2,jm1),g12(1,2,jm1),g21(1,2,jm1),g22(1,2,jm1),
     .      rhofac(jm1),ndv)
      enddo
c
      gamsq=dreal(xksq(1,nlay)) - dreal(xkhsq)
      xilay(1,nlay)=gamsq
      xilay(2,nlay)=gamsq
      migam2(1)=sqrt(-gamsq)
      g12x=g12(1,2,1)
      g22x=g22(1,2,1)
      e_k(1)=g22x + migam2(1)*g12x
      migam2(2)=xkh/migam2(1)
c: Include derivative of exponential factor so that e_k(2) is the derivative
c: of the part without the exponential term:
      capx_k=g_exp(2,1)
      g22_k=g22(2,2,1) - g22x*capx_k
      g12_k=g12(2,2,1) - g12x*capx_k
      e_k(2)=g22_k + migam2(1)*g12_k + migam2(2)*g12x
      if(ndv .gt. 2) then
         migam2(3)=-xksq(1,nlay)/(migam2(1)*w)
         capx_w=g_exp(3,1)
cc       e_k(3)=g22(3) + migam2(1)*g12(3) + migam2(3)*g12(1)
         g22_w=g22(3,2,1) - g22x*capx_w
         g12_w=g12(3,2,1) - g12x*capx_w
         e_k(3)=g22_w + migam2(1)*g12_w + migam2(3)*g12x
         if(ndv .gt. 3) then
            capx_fac=capx_k*capx_k + g_exp(5,1)
            g22_kk=g22(5,2,1) - g22x*capx_fac - 2.d0*g22_k*capx_k
            g12_kk=g12(5,2,1) - g12x*capx_fac - 2.d0*g12_k*capx_k
            migam2(5)=xksq(1,nlay)/(migam2(1)*gamsq)
            e_k(5)=g22_kk + migam2(1)*g12_kk + migam2(5)*g12x +
     .         2.d0*migam2(2)*g12_k
         endif
      endif
      phm=phtot/pie
c
      return
      end
ccc
      subroutine rx_norm(jlay,ii,ndv,vg,D2k_Df)
c
c: Computes mode normalization and finds field at top of lower halfspace
c: (see ORCA II, p. 43).
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jlay,ndv,j,jp1,jdv,jlaymax,j1,j2,kd,ntry,iipol,ii
      real*8 rho_ref,dzero,phmx,err(25),errmin,vg,Dvg_w0,D2k_Df,
     .   knr,magsq,Lmag
      complex*16 igam0(6),igam1(6),igam2(6),N_n,DN_n,DR1_w,
     .   Aigam(6),DWrk_w,Digam1,Digam2,Digam0,DTC2,Aifac,vgc
      data dzero/0.d0/
c
c: Need w derivative for group velocity:
      iipol=0
      iidone=0
      kd=kduct0
10    continue
      ntry=0
      call rx_igam_calc(xksq(1,nlay),xkh,xkhsq,w,wsq,igam2,ndv)
      call rx_igam_calc(xksq(1,2),xkh,xkhsq,w,wsq,igam0,ndv)
c
      jlay=nsvmin
      ii=isvmin
      jlaymax=jlay
c: Match at bottom of layer SVP minimum:
c
c: Phi(j) = H Phi(1), where values are taken at BOTTOMS of layers [except
c: for the upper halfspace layer (layer 1)]:
      jp1=1
      h11(1,2,jp1)=1.d0
      h12(1,2,jp1)=dzero
      h21(1,2,jp1)=dzero
      h22(1,2,jp1)=1.d0
      h_exp(1,jp1)=dzero
      do j=2,ndv
         h11(j,2,jp1)=dzero
         h12(j,2,jp1)=dzero
         h21(j,2,jp1)=dzero
         h22(j,2,jp1)=dzero
         h_exp(j,jp1)=dzero
      enddo
c
      j1=1
20    continue
      j2=jlay
      if(ii .eq. 1) j2=jlay-1
      do j=j1,j2
         jp1=j + 1
c: Take H across interface from bottom of layer j to top of layer jp1:
         call rx_iface(h11(1,2,j),h12(1,2,j),h21(1,2,j),h22(1,2,j),
     .      h11(1,1,jp1),h12(1,1,jp1),h21(1,1,jp1),h22(1,1,jp1),
     .      1.d0/rhofac(j),ndv)
c: Take H from top of layer jp1 to bottom of layer jp1:
         call rx_2x2inv(h11(1,1,jp1),h12(1,1,jp1),h21(1,1,jp1),
     .      h22(1,1,jp1),u11(1,jp1),u12(1,jp1),u21(1,jp1),u22(1,jp1),
     .      h11(1,2,jp1),h12(1,2,jp1),h21(1,2,jp1),h22(1,2,jp1),ndv)
         do jdv=1,ndv
            h_exp(jdv,jp1)=h_exp(jdv,j) + u_exp(jdv,jp1)
         enddo
      enddo
      if(ii .eq. 1) then
         j=jlay-1
         jp1=jlay
c: Take H across interface from bottom of layer j to top of layer jp1:
         call rx_iface(h11(1,2,j),h12(1,2,j),h21(1,2,j),h22(1,2,j),
     .      h11(1,1,jp1),h12(1,1,jp1),h21(1,1,jp1),h22(1,1,jp1),
     .      1.d0/rhofac(j),ndv)
      endif
c
30    continue
c: Vertical wavenumbers at match point:
      call rx_igam_calc(xksq(ii,jlay),xkh,xkhsq,w,wsq,igam1,ndv)
c: Compute downward-looking reflection coefficient:
      call rx_rp_flay(g11(1,ii,jlay),g12(1,ii,jlay),g21(1,ii,jlay),
     .   g22(1,ii,jlay),igam1,igam2,rho_prod(1,jlay),ndv,R1,TC1)
c: Compute upward-looking reflection coefficient:
      call rx_rp_vac(h11(1,ii,jlay),h12(1,ii,jlay),h21(1,ii,jlay),
     .   h22(1,ii,jlay),igam1,igam0,rho_prod(2,jlay),ndv,R2,TC2,
     .   -1.d0,jlay,ii)
c
c: Check if R1*R2 close enough to 1:
      ntry=ntry + 1
      rx_r1r2(1)=R1(1)*R2(1) - 1.d0
      rx_r1r2(2)=R1(1)*R2(2) + R1(2)*R2(1)
      Lmag=cdabs(rx_r1r2(1)/rx_r1r2(2))
cc    if(rx_r1r2(2) .eq. (0.d0,0.d0)) then
cc       print *,'zero'
cc    call rx_rp_flay(g11(1,ii,jlay),g12(1,ii,jlay),g21(1,ii,jlay),
cc   .   g22(1,ii,jlay),igam1,igam2,rho_prod(1,jlay),ndv,R1,TC1)
cc    endif
      err(ntry)=magsq(rx_r1r2(1))
      jlayx(ntry)=jlay
cc    if(iidone .eq. 0 .and. err(ntry) .gt. 1.d-8) then
c: R1*R2 - 1 can be far from zero even when Lmag is small:
      if(iidone .eq. 0 .and. (Lmag .gt. dk_max .or. 
     .   err(ntry) .gt. 1.d-6)) then
40       continue
         if(ntry .lt. nduct) then
            kd=mod(kduct0+ntry-1,nduct) + 1
            jlay=jduct(1,kd)
            ii=jduct(2,kd)
            if(dreal(xkh) .gt. dreal(xk(ii,jlay))) then
c: Don't bother to check duct if mode is not propagating in that duct:
               ntry=ntry + 1
               err(ntry)=1.d100
               goto 40
            endif
            if(jlay .gt. jlaymax) then
               j1=jlaymax
               jlaymax=jlay
               goto 20
            endif
            goto 30
         elseif(iipol .eq. 0) then
c: Polish root to see if more accurate value of kn needed:
            kn(nmode)=kn(nmode) - ekn(1,nmode)/ekn(2,nmode)
            call rx_ek_calc(kn(nmode),ekn(1,nmode),phmx,ndv)
            mode_phz(1,nmode)=phmx
            ncalc(nmode)=ncalc(nmode) + 1
            iipol=1
            if(iidiag .ge. 2) print *,'Polished roots ...',
     .   (err(j),j=1,nduct)
            goto 10
         elseif(iipol .eq. 1) then
c: Polish root again to see if more accurate value of kn needed:
            kn(nmode)=kn(nmode) - ekn(1,nmode)/ekn(2,nmode)
            call rx_ek_calc(kn(nmode),ekn(1,nmode),phmx,ndv)
            mode_phz(1,nmode)=phmx
            ncalc(nmode)=ncalc(nmode) + 1
            iipol=2
            if(iidiag .ge. 2) print *,'Polished roots again ...',
     .   (err(j),j=1,nduct)
            goto 10
         else
c: If still error after polishing, take match point with smallest error:
            errmin=err(1)
            jlay=jlayx(1)
            do j=2,nduct
               if(err(j) .lt. errmin) then
                  errmin=err(j)
                  jlay=jlayx(j)
               endif
            enddo
c: Set flag so that this is the last try:
            iidone=1
      if(iidiag .ge. 2) print *,'Using best match point: ',nmode,jlay,
     .   (err(j),j=1,nduct)
c: FIX (4-1-96):
            goto 30
cc          goto 10
         endif
      endif
c: Form Wronskian from downward and upward looking reflection coeffs:
cc    print *,'args = ',datan2(dimag(R1(1)),dreal(R1(1)))*180./pie,
cc   .   datan2(dimag(R2(1)),dreal(R2(1)))*180./pie
      knr=dreal(kn(nmode))
      rx_r1r2(3)=R1(1)*R2(3) + R1(3)*R2(1)
c
c: Group velocity:
      vgc=-rx_r1r2(2)/rx_r1r2(3)
      vg=dreal(vgc)
c: Better derivatives are found when small error is not included:
      Aifac=-2.d0*R1(1)
      Aigam(1)=Aifac*igam1(1)
      Aigam(2)=Aifac*igam1(2)
      Wr(1)=Aigam(1)*rx_r1r2(1)
      Wr(2)=Aigam(1)*rx_r1r2(2) + Aigam(2)*rx_r1r2(1)
      rho_ref=geo(ii,3,jlay)
      N_n=cdsqrt(-2.d0*knr/(Wr(2)*rho_ref))
c
c: Compute potential and derivative at water surface:
      philay(1,2)=dcmplx(0.d0,0.d0)
c: Factor of 2 from 1-R=2, minus sign from looking upward:
      dphilay(1,2)=-2.d0*igam0(1)*R1(1)*TC2(1)*N_n
      Aplay(1,2)=-h_exp(1,jlay)
c: Compute potential and derivative at top of basement:
      philay(1,nlay)=TC1(1)*N_n
      dphilay(1,nlay)=igam2(1)*philay(1,nlay)
c: Use g_exp value in layer above if ref depth is at top of layer:
      Aplay(1,nlay)=-g_exp(1,jlay+ii-2)
c
c: Compute potential and derivative at reference depth and just below:
      philay(ii,jlay)=N_n*(1.d0 + R1(1))
      dphilay(ii,jlay)=igam1(1)*N_n*(1.d0 - R1(1))
      Aplay(ii,jlay)=0.d0
      if(ii .eq. 2) then
         philay(1,jlay+1)=philay(ii,jlay)*rhofac(jlay)
         dphilay(1,jlay+1)=dphilay(ii,jlay)
         Aplay(1,jlay+1)=0.d0
      else
         philay(2,jlay-1)=philay(ii,jlay)/rhofac(jlay-1)
         dphilay(2,jlay-1)=dphilay(ii,jlay)
         Aplay(2,jlay-1)=0.d0
      endif
c
      if(ndv .ge. 4) then
         rx_r1r2(4)=R1(1)*R2(4) + R1(4)*R2(1) +
     .              R1(2)*R2(3) + R1(3)*R2(2)
         rx_r1r2(5)=R1(1)*R2(5) + R1(5)*R2(1) + 2.d0*R1(2)*R2(2)
         rx_r1r2(6)=R1(1)*R2(6) + R1(6)*R2(1) + 2.d0*R1(3)*R2(3)
c: Total derivative of group velocity:
         Dvg_w0=dreal((vgc*vgc*rx_r1r2(6) + 2.d0*vgc*rx_r1r2(4) +
     .      rx_r1r2(5))/rx_r1r2(2))
c: Second total derivative of kn w.r.t. f (see ORCA II, p.63):
         D2k_Df=(-(twpie/vgc)**2)*Dvg_w0
cc       Aigam(2)=-2.d0*(R1(1)*igam1(2) + R1(2)*igam1(1))
cc       Aigam(3)=-2.d0*(R1(1)*igam1(3) + R1(3)*igam1(1))
cc       Aigam(4)=-2.d0*(R1(1)*igam1(4) + R1(4)*igam1(1) + 
cc   .      R1(2)*igam1(3) + R1(3)*igam1(3))
cc       Aigam(5)=-2.d0*(R1(1)*igam1(5) + R1(5)*igam1(1) + 
cc   .      2.d0*R1(2)*igam1(2))
cc       Aigam(6)=-2.d0*(R1(1)*igam1(6) + R1(6)*igam1(1) + 
cc   .      2.d0*R1(3)*igam1(3))
c: R1 is constant A here, so no derivatives should be used(?):
         Aigam(3)=Aifac*igam1(3)
         Aigam(4)=Aifac*igam1(4)
         Aigam(5)=Aifac*igam1(5)
         Aigam(6)=Aifac*igam1(6)
c: Not used:
cc       Wr(3)=Aigam(1)*rx_r1r2(3) + Aigam(3)*rx_r1r2(1)
         Wr(4)=Aigam(1)*rx_r1r2(4) + Aigam(3)*rx_r1r2(2) + 
     .      Aigam(2)*rx_r1r2(3) + Aigam(4)*rx_r1r2(1)
         Wr(5)=Aigam(1)*rx_r1r2(5) + 2.d0*Aigam(2)*rx_r1r2(2) +
     .      Aigam(5)*rx_r1r2(1)
         Wr(6)=Aigam(1)*rx_r1r2(6) + 2.d0*Aigam(3)*rx_r1r2(3) +
     .      Aigam(6)*rx_r1r2(1)
c: Total derivativ of N_n w.r.t. w (p. 63, ORCA II):
         DWrk_w=Wr(4) + Wr(5)/vgc
         DN_n=-(Wr(2)/vgc - DWrk_w*knr)/(N_n*rho_ref*Wr(2)*Wr(2))
c: Compute total derivative w.r.t. w of phi,dphi at water surface:
         Dphi_w(1,2)=dcmplx(0.d0,0.d0)
         DR1_w=R1(3) + R1(2)/vgc
         DTC2=TC2(3) + TC2(2)/vgc
         Digam0=igam0(3) + igam0(2)/vgc
         Ddphi_w(1,2)=-2.d0*(igam0(1)*R1(1)*(TC2(1)*DN_n + DTC2*N_n)
     .      + TC2(1)*N_n*(igam0(1)*DR1_w + Digam0*R1(1)))
c: Compute total derivative w.r.t. w of phi,dphi at top of basement:
         Dphi_w(1,nlay)=(TC1(3) + TC1(2)/vgc)*N_n + TC1(1)*DN_n
         Digam2=igam2(3)+igam2(2)/vgc
         Ddphi_w(1,nlay)=Digam2*philay(1,nlay) +
     .      igam2(1)*Dphi_w(1,nlay)
c: Compute total derivative w.r.t. w of phi,dphi at reference depth:
cc       philay(ii,jlay)=N_n*(1.d0 + R1(1))
         Dphi_w(ii,jlay)=N_n*DR1_w + DN_n*(1.d0 + R1(1))
         Digam1=igam1(3) + igam1(2)/vgc
         Ddphi_w(ii,jlay)=-igam1(1)*N_n*DR1_w + 
     .      (igam1(1)*DN_n + Digam1*N_n)*(1.d0 - R1(1))
         if(ii .eq. 2) then
            Dphi_w(1,jlay+1)=Dphi_w(ii,jlay)*rhofac(jlay)
            Ddphi_w(1,jlay+1)=Ddphi_w(ii,jlay)
         else
            Dphi_w(2,jlay-1)=Dphi_w(ii,jlay)/rhofac(jlay-1)
            Ddphi_w(2,jlay-1)=Ddphi_w(ii,jlay)
         endif
      endif
cc
      eig_char(1,nmode)=phmx
      eig_char(2,nmode)=rx_r1r2(2)
      eig_char(3,nmode)=rx_r1r2(3)
      eig_char(4,nmode)=dcmplx(vg,0.d0)
      eig_char(5,nmode)=dcmplx(D2k_Df,0.d0)
      nzref(nmode)=kd
c
      eig_char(1,nmode)=vgc
      eig_char(2,nmode)=Dvg_w0
cc    eig_char(1,nmode)=N_n
cc    eig_char(2,nmode)=DN_n
cc    eig_char(1,nmode)=TC2(1)*dexp(-h_exp(1,jlay))
cc    eig_char(2,nmode)=(TC2(3) + TC2(2)/vgc)*dexp(-h_exp(1,jlay))
cc    eig_char(1,nmode)=TC1(1)*dexp(-g_exp(1,jlay))
cc    eig_char(2,nmode)=(TC1(3) + TC1(2)/vgc)*dexp(-g_exp(1,jlay))
c
cc    eig_char(1,nmode)=dphilay(1,nlay)*dexp(-g_exp(1,jlay))
cc    eig_char(2,nmode)=Ddphi_w(1,nlay)*dexp(-g_exp(1,jlay))
c: 10-11-95: I've now checked total derivatives of phi,dphi at top and bottom
c
      return
      end
ccc
      subroutine rx_igam_calc(Ksq,xkh,xkhsq,w,wsq,igam0,ndv)
c
      implicit none
      integer*4 ndv
      real*8 Ksq,xkh,xkhsq,w,wsq,gamsq
      complex*16 igam0(6),gamma
c
      gamsq=Ksq - xkhsq
      gamma=cdsqrt(dcmplx(gamsq,0.d0))
      igam0(1)=dcmplx(-dimag(gamma),dreal(gamma))
      igam0(2)=xkh/igam0(1)
      if(ndv .gt. 2) then
         igam0(3)=-Ksq/(w*igam0(1))
         if(ndv .gt. 3) then
            igam0(4)=-igam0(2)*igam0(3)/igam0(1)
            igam0(5)=Ksq/(igam0(1)*gamsq)
            igam0(6)=Ksq*(Ksq/gamsq - 1.d0)/(wsq*igam0(1))
         endif
      endif
c
      return
      end
ccc
      subroutine rx_2x2(e11,e12,e21,e22,f11,f12,f21,f22,
     .   g11,g12,g21,g22,ndv)
c: Multiplies 2x2 matrices E*F=G and computes derivatives up to order ndv.
c
      implicit none
      real*8 e11(6),e12(6),e21(6),e22(6),f11(6),f12(6),f21(6),f22(6),
     .   g11(6),g12(6),g21(6),g22(6)
      integer*4 ndv,j
c
      g11(1)=e11(1)*f11(1) + e12(1)*f21(1)
      g21(1)=e21(1)*f11(1) + e22(1)*f21(1)
      g12(1)=e11(1)*f12(1) + e12(1)*f22(1)
      g22(1)=e21(1)*f12(1) + e22(1)*f22(1)
      if(ndv .ge. 2) then
         j=2
         g11(j)=e11(1)*f11(j)+e11(j)*f11(1)+e12(1)*f21(j)+e12(j)*f21(1)
         g21(j)=e21(1)*f11(j)+e21(j)*f11(1)+e22(1)*f21(j)+e22(j)*f21(1)
         g12(j)=e11(1)*f12(j)+e11(j)*f12(1)+e12(1)*f22(j)+e12(j)*f22(1)
         g22(j)=e21(1)*f12(j)+e21(j)*f12(1)+e22(1)*f22(j)+e22(j)*f22(1)
      endif
      if(ndv .ge. 3) then
         j=3
         g11(j)=e11(1)*f11(j)+e11(j)*f11(1)+e12(1)*f21(j)+e12(j)*f21(1)
         g21(j)=e21(1)*f11(j)+e21(j)*f11(1)+e22(1)*f21(j)+e22(j)*f21(1)
         g12(j)=e11(1)*f12(j)+e11(j)*f12(1)+e12(1)*f22(j)+e12(j)*f22(1)
         g22(j)=e21(1)*f12(j)+e21(j)*f12(1)+e22(1)*f22(j)+e22(j)*f22(1)
      endif
c: Do second derivative w.r.t. k and w:
      if(ndv .ge. 4) then
         j=4
         g11(j)=e11(1)*f11(j)+e11(j)*f11(1)+e12(1)*f21(j)+e12(j)*f21(1)
     .      + e11(2)*f11(3)+e11(3)*f11(2)+e12(2)*f21(3)+e12(3)*f21(2)
         g21(j)=e21(1)*f11(j)+e21(j)*f11(1)+e22(1)*f21(j)+e22(j)*f21(1)
     .      + e21(2)*f11(3)+e21(3)*f11(2)+e22(2)*f21(3)+e22(3)*f21(2)
         g12(j)=e11(1)*f12(j)+e11(j)*f12(1)+e12(1)*f22(j)+e12(j)*f22(1)
     .      + e11(2)*f12(3)+e11(3)*f12(2)+e12(2)*f22(3)+e12(3)*f22(2)
         g22(j)=e21(1)*f12(j)+e21(j)*f12(1)+e22(1)*f22(j)+e22(j)*f22(1)
     .      + e21(2)*f12(3)+e21(3)*f12(2)+e22(2)*f22(3)+e22(3)*f22(2)
c: Do second derivatives w.r.t. k and w:
         j=5
         g11(j)=e11(1)*f11(j)+e11(j)*f11(1)+e12(1)*f21(j)+
     .      e12(j)*f21(1) + 2.d0*(e11(2)*f11(2)+e12(2)*f21(2))
         g21(j)=e21(1)*f11(j)+e21(j)*f11(1)+e22(1)*f21(j)+
     .      e22(j)*f21(1) + 2.d0*(e21(2)*f11(2)+e22(2)*f21(2))
         g12(j)=e11(1)*f12(j)+e11(j)*f12(1)+e12(1)*f22(j)+
     .      e12(j)*f22(1) + 2.d0*(e11(2)*f12(2)+e12(2)*f22(2))
         g22(j)=e21(1)*f12(j)+e21(j)*f12(1)+e22(1)*f22(j)+
     .      e22(j)*f22(1) + 2.d0*(e21(2)*f12(2)+e22(2)*f22(2))
         j=6
         g11(j)=e11(1)*f11(j)+e11(j)*f11(1)+e12(1)*f21(j)+
     .      e12(j)*f21(1) + 2.d0*(e11(3)*f11(3)+e12(3)*f21(3))
         g21(j)=e21(1)*f11(j)+e21(j)*f11(1)+e22(1)*f21(j)+
     .      e22(j)*f21(1) + 2.d0*(e21(3)*f11(3)+e22(3)*f21(3))
         g12(j)=e11(1)*f12(j)+e11(j)*f12(1)+e12(1)*f22(j)+
     .      e12(j)*f22(1) + 2.d0*(e11(3)*f12(3)+e12(3)*f22(3))
         g22(j)=e21(1)*f12(j)+e21(j)*f12(1)+e22(1)*f22(j)+
     .      e22(j)*f22(1) + 2.d0*(e21(3)*f12(3)+e22(3)*f22(3))
      endif
c
      return
      end
ccc
      subroutine rx_2x2inv(e11,e12,e21,e22,h11,h12,h21,h22,
     .   g11,g12,g21,g22,ndv)
c: Multiplies 2x2 matrices E*H^(-1)=G and computes derivatives up to order ndv.
c: rx_2x2inv takes the inverse of the propagator matrix H by setting
c: f11 ==> h22, f22 ==> h11, f12 ==> -h12, f21 ==> -h21.  NO!!!
c: Must set
c: f11 ==> h22, f22 ==> h11, f12 ==> h12, f21 ==> h21.
c: in order to work.  Minus signs go away because slope eta changes sign also.
c
      implicit none
      real*8 e11(6),e12(6),e21(6),e22(6),h11(6),h12(6),h21(6),h22(6),
     .   g11(6),g12(6),g21(6),g22(6)
      integer*4 ndv,j
c
      g11(1)=e11(1)*h22(1) - e12(1)*h21(1)
      g21(1)=e21(1)*h22(1) - e22(1)*h21(1)
      g12(1)=-e11(1)*h12(1) + e12(1)*h11(1)
      g22(1)=-e21(1)*h12(1) + e22(1)*h11(1)
      if(ndv .ge. 2) then
         j=2
         g11(j)=e11(1)*h22(j)+e11(j)*h22(1)-e12(1)*h21(j)-e12(j)*h21(1)
         g21(j)=e21(1)*h22(j)+e21(j)*h22(1)-e22(1)*h21(j)-e22(j)*h21(1)
         g12(j)=-e11(1)*h12(j)-e11(j)*h12(1)+e12(1)*h11(j)+e12(j)*h11(1)
         g22(j)=-e21(1)*h12(j)-e21(j)*h12(1)+e22(1)*h11(j)+e22(j)*h11(1)
      endif
      if(ndv .ge. 3) then
         j=3
         g11(j)=e11(1)*h22(j)+e11(j)*h22(1)-e12(1)*h21(j)-e12(j)*h21(1)
         g21(j)=e21(1)*h22(j)+e21(j)*h22(1)-e22(1)*h21(j)-e22(j)*h21(1)
         g12(j)=-e11(1)*h12(j)-e11(j)*h12(1)+e12(1)*h11(j)+e12(j)*h11(1)
         g22(j)=-e21(1)*h12(j)-e21(j)*h12(1)+e22(1)*h11(j)+e22(j)*h11(1)
      endif
      if(ndv .ge. 4) then
         j=4
         g11(j)=e11(1)*h22(j)+e11(j)*h22(1)-e12(1)*h21(j)-e12(j)*h21(1)
     .      + e11(2)*h22(3)+e11(3)*h22(2)-e12(2)*h21(3)-e12(3)*h21(2)
         g21(j)=e21(1)*h22(j)+e21(j)*h22(1)-e22(1)*h21(j)-e22(j)*h21(1)
     .      + e21(2)*h22(3)+e21(3)*h22(2)-e22(2)*h21(3)-e22(3)*h21(2)
         g12(j)=-e11(1)*h12(j)-e11(j)*h12(1)+e12(1)*h11(j)+e12(j)*h11(1)
     .      - e11(2)*h12(3)-e11(3)*h12(2)+e12(2)*h11(3)+e12(3)*h11(2)
         g22(j)=-e21(1)*h12(j)-e21(j)*h12(1)+e22(1)*h11(j)+e22(j)*h11(1)
     .      - e21(2)*h12(3)-e21(3)*h12(2)+e22(2)*h11(3)+e22(3)*h11(2)
         j=5
         g11(j)=e11(1)*h22(j)+e11(j)*h22(1)-e12(1)*h21(j)-e12(j)*h21(1)
     .      + 2.d0*(e11(2)*h22(2)-e12(2)*h21(2))
         g21(j)=e21(1)*h22(j)+e21(j)*h22(1)-e22(1)*h21(j)-e22(j)*h21(1)
     .      + 2.d0*(e21(2)*h22(2)-e22(2)*h21(2))
         g12(j)=-e11(1)*h12(j)-e11(j)*h12(1)+e12(1)*h11(j)+e12(j)*h11(1)
     .      + 2.d0*(-e11(2)*h12(2)+e12(2)*h11(2))
         g22(j)=-e21(1)*h12(j)-e21(j)*h12(1)+e22(1)*h11(j)+e22(j)*h11(1)
     .      + 2.d0*(-e21(2)*h12(2)+e22(2)*h11(2))
         j=6
         g11(j)=e11(1)*h22(j)+e11(j)*h22(1)-e12(1)*h21(j)-e12(j)*h21(1)
     .      + 2.d0*(e11(3)*h22(3)-e12(3)*h21(3))
         g21(j)=e21(1)*h22(j)+e21(j)*h22(1)-e22(1)*h21(j)-e22(j)*h21(1)
     .      + 2.d0*(e21(3)*h22(3)-e22(3)*h21(3))
         g12(j)=-e11(1)*h12(j)-e11(j)*h12(1)+e12(1)*h11(j)+e12(j)*h11(1)
     .      + 2.d0*(-e11(3)*h12(3)+e12(3)*h11(3))
         g22(j)=-e21(1)*h12(j)-e21(j)*h12(1)+e22(1)*h11(j)+e22(j)*h11(1)
     .      + 2.d0*(-e21(3)*h12(3)+e22(3)*h11(3))
      endif
c
      return
      end
ccc
      subroutine rx_iface(h11t,h12t,h21t,h22t,h11b,h12b,h21b,h22b,
     .   rho1_rho2,ndv)
c
c: Takes H across interface from bottom of layer j to top of layer jp1:
c
      implicit none
      integer*4 ndv,jdv
      real*8 h11t(6),h12t(6),h21t(6),h22t(6),h11b(6),h12b(6),
     .   h21b(6),h22b(6),rho1_rho2
c
      if(rho1_rho2 .ne. 1.d0) then
         do jdv=1,ndv
            h11b(jdv)=h11t(jdv)*rho1_rho2
            h21b(jdv)=h21t(jdv)*rho1_rho2
            h12b(jdv)=h12t(jdv)
            h22b(jdv)=h22t(jdv)
         enddo
      else
         do jdv=1,ndv
            h11b(jdv)=h11t(jdv)
            h21b(jdv)=h21t(jdv)
            h12b(jdv)=h12t(jdv)
            h22b(jdv)=h22t(jdv)
         enddo
      endif
c
      return
      end
ccc
      subroutine rx_ep_calc(xkh,w,wsq,xkhsq,xk1sq,xk2sq,eta,etasq,h,
     .   e11,e12,e21,e22,zzexp,iso,ndv,ailay,bilay,xilay,zetalay,
     .   ii1,ii2,phtot)
c
c: Computes the propagator matrix and its derivative for a given 
c: 1/c^2 linear layer: see Notebook92, p. 69-70.
c: xkh=horizontal wavenumber; xkqsq=xkh^2; (xk1sq=(w/c1)^2; 
c: xk2sq=(w/c2)^2; eta=(w^2*beta)^(1/3); etasq=eta^2; beta=gradient;
c: h=layer thickness.
c: E and E'=dE/dk must be multiplied by exp(zzexp) to get true values.
c: For ndv=2, the derivative w.r.t. k is computed.
c: For ndv=3, the derivative w.r.t. w is computed.
c: For ndv=4, the second derivatives w.r.t. k,w k,k w,w are computed.
c: For iso=-1, the layer thickness is assumed to be inversely proportional
c: to frequency (used for faster broadband calculations with false bottoms).
c 
      use scairy_com
      integer*4 iso,ndv,ii1,ii2
      real*8 xkh,xkhsq,xk1sq,xk2sq,eta,etasq,e11(6),e12(6),e21(6),
     .   e22(6),zzexp(6),xilay(2),phlay,phtot,pie,twpie
      real*8 xi1,xi2,ai1,aip1,bi1,bip1,ai2,aip2,bi2,bip2,
     .   zzexp1(3),zzexp2(3),mgam1sq,mgam2sq,gam1cub,gam2cub
      real*8 gam,gamh,cosgamh,singamh,igam,igamh,ei2gamh,
     .   cosigamh,sinigamh,gamsq,h,hsq,hcub,w,wsq,dzero,khprod,
     .   khprodsq,Kfac,Kfac_w,Kfac_k,dh_dw,d2h_dw
      complex*16 ailay(2,2),bilay(2,2),zetalay(2)
      data dzero/0.d0/,pie/3.14159265358979/,twpie/6.28318530717959/
c
cc    if(h .eq. 0.d0) return
      if(iso .eq. 1 .or. iso .eq. -1) then
c: Compute propagator matrix for isospeed layer:
         gamsq=xk1sq - xkhsq
         xilay(1)=gamsq
         xilay(2)=gamsq
         khprod=xkh*h
         if(abs(gamsq) .lt. 1.d-11) then
c: Vertical wavenumber essentially zero:
            e11(1)=1.d0
            e22(1)=1.d0
            e12(1)=h
            e21(1)=0.d0
            zzexp(1)=dzero
            zzexp(2)=dzero
            zzexp(5)=dzero
            phlay=0.d0
            if(ndv .ge. 2) then
c: k derivatives (use L'Hopital's rule for propagating case) (see p. 91):
               e11(2)=khprod*h
               e22(2)=e11(2)
c: Should this be: e11(2)*h/2.d0?
cc             e12(2)=e11(2)*h/3.d0
               e12(2)=0.5d0*e11(2)*h
               e21(2)=2.*khprod
               if(iso .gt. 0) then
                  if(ndv .ge. 3) then
c: w derivatives (same as k, but replace -xksq/w for k):
                     Kfac=-xk1sq/(w*xkh)
                     e11(3)=Kfac*e11(2)
                     e22(3)=e11(3)
                     e12(3)=Kfac*e12(2)
                     e21(3)=Kfac*e21(2)
                  endif
                  if(ndv .ge. 4) then
c: k&w derivatives (see ORCA II, p.55):
                     hsq=h*h
                     hcub=h*hsq
                     e11(4)=khprod*e12(3)
                     e22(4)=e11(4)
                     e12(4)=0.2d0*h*e11(4)
cc                   e21(4)=Kfac*xkhsq*hcub/(0.75d0)
                     e21(4)=4.d0*e11(4)/h
                     khprodsq=khprod*khprod
                     e11(5)=hsq*(1.d0 + khprodsq/3.d0)
                     e22(5)=e11(5)
cc                   e12(5)=hcub*(1.d0 + .6d0*khprodsq)/3.d0
cc                   e21(5)=h*(2.d0 + khprodsq/(.75d0))
c: Usual form of these equations is faster:
                     e12(5)=e12(4)/Kfac + e12(2)/xkh
                     e21(5)=e21(4)/Kfac + e21(2)/xkh
c: Use regular form of these equations also:
                     Kfac_w=Kfac/w
                     e11(6)=Kfac*e11(4) + Kfac_w*e11(2)
                     e22(6)=e11(6)
                     e12(6)=Kfac*e12(4) + Kfac_w*e12(2)
                     e21(6)=Kfac*e21(4) + Kfac_w*e21(2)
                  endif
               else
c: Case for frequency-dependent thickness h:
c: Frequency dependent depth h(w)=n*lambda=2*pi*n/K:
                  if(ndv .ge. 3) then
                     dh_dw=-h/w
                     Kfac=-xk1sq/(w*xkh)
                     e11(3)=Kfac*e11(2)
                     e22(3)=e11(3)
c: Must keep track of "usual" expressions in order to use in subsequent
c: "usual" derivatives that involve those expressions, "usual" meaning
c: without the h(w) dependence:
                     e12(3)=Kfac*e12(2) + dh_dw
                     e21(3)=Kfac*e21(2)
                  endif
                  if(ndv .ge. 4) then
                     hsq=h*h
                     hcub=h*hsq
c: k&k derivatives (see ORCA II, p.55) (same):
                     khprodsq=khprod*khprod
                     e11(5)=hsq*(1.d0 + khprodsq/3.d0)
                     e22(5)=e11(5)
                     e12(5)=hcub*(1.d0 + .6d0*khprodsq)/3.d0
                     e21(5)=h*(2.d0 + khprodsq/(.75d0))
c: k&w derivatives (see ORCA II, p.73) (different):
                     Kfac_k=-Kfac/xkh
                     e11(4)=Kfac*e11(5) + Kfac_k*e11(2) + dh_dw*e21(2)
                     e22(4)=e11(4)
                     e12(4)=Kfac*e12(5) + Kfac_k*e12(2) + dh_dw*e11(2)
                     e21(4)=Kfac*e21(5) + Kfac_k*e21(2) - dh_dw*
     .                  (-2.d0*xkh*e11(1))
c
c: Use regular form of these equations also:
                     Kfac_w=Kfac/w
                     d2h_dw=2.d0*h/wsq
                     e11(6)=Kfac*e11(4) + Kfac_w*e11(2) + 
     .                  dh_dw*e21(3) + d2h_dw*e21(1)
                     e22(6)=e11(6)
                     e12(6)=Kfac*e12(4) + Kfac_w*e12(2) +
     .                  dh_dw*e11(3) + d2h_dw*e11(1)
                     e21(6)=Kfac*e21(4) + Kfac_w*e21(2) -
     .                  dh_dw*(2.d0*xk1sq*e11(1)/w)
                  endif
               endif
            endif
            zetalay(ii1)=dcmplx(0.d0,phlay)
            zetalay(ii2)=dcmplx(0.d0,phlay)
            return
         elseif(gamsq .gt. 0.d0) then
c: Case where plane wave is propagating:
cc    print *,'isospeed prop ',w/dreal(xkh),w/dreal(sqrt(xk1sq))
            gam=sqrt(gamsq)
            gamh=h*gam
            cosgamh=cos(gamh)
            singamh=sqrt(1.d0 - cosgamh*cosgamh)
            if(mod(gamh,twpie) .gt. pie) singamh=-singamh
            e11(1)=cosgamh
            e22(1)=cosgamh
            e12(1)=singamh/gam
            e21(1)=-gam*singamh
            zzexp(1)=dzero
            zzexp(2)=dzero
            zzexp(5)=dzero
            phlay=gamh
            phtot=phtot + phlay
            zetalay(ii1)=dcmplx(0.d0,phlay)
            zetalay(ii2)=dcmplx(gamh,phlay)
         elseif(gamsq .lt. 0.d0) then
cc    print *,'isospeed evan ',w/dreal(xkh),w/dreal(sqrt(xk1sq))
cc    print *,'determinant check: ',(e11(1)*e22(1)-e12(1)*e21(1))
cc   .   *exp(2.*zzexp)
c: Case where plane wave is evanescent:
            igam=sqrt(-gamsq)
            igamh=igam*h
            ei2gamh=0.5d0*exp(-2.*igamh)
            sinigamh=0.5d0 - ei2gamh
            cosigamh=0.5d0 + ei2gamh
            e11(1)=cosigamh
            e22(1)=cosigamh
            e12(1)=sinigamh/igam
            e21(1)=igam*sinigamh
            zzexp(1)=igamh
            zzexp(2)=khprod/igam
            phlay=0.d0
            zetalay(ii1)=dcmplx(0.d0,phlay)
            zetalay(ii2)=dcmplx(igamh,phlay)
         endif
c: Derivatives can be expressed in terms of functions themselves, 
c: independent of whether field is propagating or evanescent 
c: (see ORCA II, p.56-58):
         if(ndv .ge. 2) then
c: k derivatives:
            e11(2)=khprod*e12(1)
            e22(2)=e11(2)
            e12(2)=-xkh*(h*e11(1) - e12(1))/gamsq
            e21(2)=xkh*(h*e11(1) + e12(1))
            if(ndv .ge. 3) then
c: w derivatives:
               if(iso .eq. 1) then
                  Kfac=-xk1sq/(w*xkh)
                  e11(3)=Kfac*e11(2)
                  e22(3)=e11(3)
                  e12(3)=Kfac*e12(2)
                  e21(3)=Kfac*e21(2)
                  if(ndv .ge. 4) then
c: k&w derivatives:
                     e11(4)=khprod*e12(3)
                     e22(4)=e11(4)
cc                   e12(4)=(-xkh*(h*e11(3) - e12(3)) - 
cc   .                  2.d0*xk1sq*e12(2)/w)/gamsq
c: Simpler form:
                     e12(4)=-xkh*(h*e11(3) - 3.d0*e12(3))/gamsq
                     e21(4)=xkh*(h*e11(3) + e12(3))
c: k&k derivatives:
                     e11(5)=khprod*e12(2) + h*e12(1)
                     e22(5)=e11(5)
cc                   e12(5)=-(xkh/gamsq)*(h*e11(2) - 3.d0*e12(2)) + 
cc   .                  e12(2)/xkh
cc                   e21(5)=xkh*(h*e11(2) + e12(2)) + e21(2)/xkh
c: Simpler forms:
                     e12(5)=e12(4)/Kfac + e12(2)/xkh
                     e21(5)=e21(4)/Kfac + e21(2)/xkh
                     if(gamsq .lt. 0.d0) zzexp(5)=h*xk1sq/(gamsq*igam)
c: w&w derivatives:
                     Kfac_w=Kfac/w
                     e11(6)=Kfac*e11(4) + Kfac_w*e11(2)
                     e22(6)=e11(6)
                     e12(6)=Kfac*e12(4) + Kfac_w*e12(2)
                     e21(6)=Kfac*e21(4) + Kfac_w*e21(2)
                  endif
               else
c: Frequency dependent depth h(w)=n*lambda=2*pi*n/K (see ORCA II, p. 72-73):
                  dh_dw=-h/w
                  Kfac=-xk1sq/(w*xkh)
                  e11(3)=Kfac*e11(2) + dh_dw*e21(1)
                  e22(3)=e11(3)
                  e12(3)=Kfac*e12(2) + dh_dw*e11(1)
                  e21(3)=Kfac*e21(2) - dh_dw*gamsq*e11(1)
                  if(ndv .ge. 4) then
c: k&k derivatives (same, see p. 57, ORCA II):
                     e11(5)=khprod*e12(2) + h*e12(1)
                     e22(5)=e11(5)
                     e12(5)=-xkh*(h*e11(2) - 3.d0*e12(2))/gamsq + 
     .                  e12(2)/xkh
                     e21(5)=xkh*(h*e11(2) + e12(2)) + e21(2)/xkh
                     if(gamsq .lt. 0.d0) zzexp(5)=h*xk1sq/(gamsq*igam)
c: k&w derivatives (different):
                     Kfac_k=-Kfac/xkh
                     e11(4)=Kfac*e11(5) + Kfac_k*e11(2) + dh_dw*e21(2)
                     e22(4)=e11(4)
                     e12(4)=Kfac*e12(5) + Kfac_k*e12(2) + dh_dw*e11(2)
                     e21(4)=Kfac*e21(5) + Kfac_k*e21(2) - dh_dw*
     .                  (gamsq*e11(2) - 2.d0*xkh*e11(1))
c: w&w derivatives (different):
                     d2h_dw=2.d0*h/wsq
                     Kfac_w=Kfac/w
                     e11(6)=Kfac*e11(4) + Kfac_w*e11(2) + 
     .                  dh_dw*e21(3) + d2h_dw*e21(1)
                     e22(6)=e11(6)
                     e12(6)=Kfac*e12(4) + Kfac_w*e12(2) +
     .                  dh_dw*e11(3) + d2h_dw*e11(1)
                     e21(6)=Kfac*e21(4) + Kfac_w*e21(2) -
     .                  dh_dw*(gamsq*e11(3) + 2.d0*xk1sq*e11(1)/w) - 
     .                  d2h_dw*gamsq*e11(1)
                  endif
               endif
            endif
         endif
cc    print *,'determinant check: ',sngl(gamsq),sngl((e11(1)*e22(1)-
cc   .   e12(1)*e21(1))*exp(2.*zzexp)),e11(1),e12(1)
cc    print *,'e11,e12 = ',e11(1)*exp(zzexp),e12(1)*exp(zzexp),
cc   .   e21(1)*exp(zzexp),e22(1)*exp(zzexp)
      else
c: Compute propagator matrix for Airy layer:
c: Compute Airy function arguments at top and bottom of layer:
         mgam1sq=xkhsq - xk1sq
         mgam2sq=xkhsq - xk2sq
         xi1=mgam1sq/etasq
         xi2=mgam2sq/etasq
cc    print *,'Airy ...',xi1,xi2
cc       xi1=(xkhsq - xk1sq)/etasq
cc       xi2=(xkhsq - xk2sq)/etasq
c: Compute Airy functions of real arguments (ai and bi):
         call rx_airy_sm(xi1,ai1,bi1,aip1,bip1,zzexp1)
         call rx_airy_sm(xi2,ai2,bi2,aip2,bip2,zzexp2)
         if(xi1 .lt. 0.d0) then
            gam1cub=-sqrt(-mgam1sq)*mgam1sq
         else
            gam1cub=0.d0
         endif
         if(xi2 .lt. 0.d0) then
            gam2cub=-sqrt(-mgam2sq)*mgam2sq
         else
            gam2cub=0.d0
         endif
         phlay=abs(gam1cub - gam2cub)/(1.5*abs(eta)*etasq)
         phtot=phtot + phlay
c
         ailay(1,ii1)=ai1
         ailay(1,ii2)=ai2
         ailay(2,ii1)=aip1
         ailay(2,ii2)=aip2
         bilay(1,ii1)=bi1
         bilay(1,ii2)=bi2
         bilay(2,ii1)=bip1
         bilay(2,ii2)=bip2
         zetalay(ii1)=dcmplx(zzexp1(1),phlay)
         zetalay(ii2)=dcmplx(zzexp2(1),phlay)
         xilay(1)=xi1
         xilay(2)=xi2
c
c: Compute propagator matrix for Airy layer:
         call rx_airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,
     .      xi2,ai2,aip2,bi2,bip2,zzexp2,eta,e11,e12,e21,e22,
     .      zzexp,xkh,xkhsq,xk1sq,xk2sq,etasq,w,wsq,ndv)
c
cc    if(abs(w/sqrt(xk1sq) - w/sqrt(xk2sq)) .lt. 1.) print *,
cc   .   'Airy: ',e11(1)*exp(zzexp),e12(1)*exp(zzexp),
cc   .   e21(1)*exp(zzexp),e22(1)*exp(zzexp)
      endif
cc    print *,'determinant check: ',(e11(1)*e22(1)-e12(1)*e21(1))
cc   .   *exp(2.*zzexp)
c
      return
      end
ccc
      subroutine rx_airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,
     .   xi2,ai2,aip2,bi2,bip2,zzexp2,eta,e11,e12,e21,e22,
     .   zzexp,xkh,xkhsq,xk1sq,xk2sq,etasq,w,wsq,ndv)
c
      implicit none
      integer*4 ndv
      real*8 xi1,ai1,aip1,bi1,bip1,zzexp1(3),xi2,ai2,aip2,
     .   bi2,bip2,zzexp2(3),eta,e11(6),e12(6),e21(6),e22(6),
     .   zzexp(6),xkh,xkhsq,xk1sq,xk2sq,etasq,delex,zfac,
     .   e21_eta,eta_e12,dxi_dk,dxi_dw1,dxi_dw2,dxi_fac,edif,
     .   twok_eta,twok_eta3,e21fac,e5fac,xi_dxi1,xi_dxi2,
     .   eta_e12_w,d2xi_dw1,d2xi_dw2,xi_dxi1_w,xi_dxi2_w,e21_eta_w,
     .   e12num,e12num_w,e21num,e21num_w,d2den,delzp
      real*8 w,wsq,deta_eta,deta_eta_w,det,sg
      data det/0.31830988618379/
c
c: 1st terms of E must be multiplied by exp(delex); 2nd terms by exp(-delex):
      delex=zzexp1(1) - zzexp2(1)
      if(delex .gt. 0.d0) then
c: Exponential factor of first of two terms dominates:
         sg=1.d0
         zzexp(1)=delex
         zfac=dexp(-2.d0*delex)
c: Remove old factor from and insert new factor to second terms:
         e11(1)=(ai2*bip1 - zfac*bi2*aip1)/det
         e12(1)=(-ai2*bi1 + zfac*bi2*ai1)/(-det*eta)
         e21(1)=-eta*(aip2*bip1 - zfac*bip2*aip1)/det
         e22(1)=(-aip2*bi1 + zfac*bip2*ai1)/det
      elseif(delex .lt. 0.d0) then
c: Exponential factor of second of two terms dominates:
         sg=-1.d0
         zzexp(1)=-delex
         zfac=dexp(2.d0*delex)
c: Remove old factor from and insert new factor into first terms:
         e11(1)=(zfac*ai2*bip1 - bi2*aip1)/det
         e12(1)=(-zfac*ai2*bi1 + bi2*ai1)/(-det*eta)
         e21(1)=-eta*(zfac*aip2*bip1 - bip2*aip1)/det
         e22(1)=(-zfac*aip2*bi1 + bip2*ai1)/det
      else
c: Exponential factor is zero:
         sg=0.d0
c: Remove old factor from and insert new factor into first terms:
         e11(1)=(ai2*bip1 - bi2*aip1)/det
         e12(1)=(-ai2*bi1 + bi2*ai1)/(-det*eta)
         e21(1)=-eta*(aip2*bip1 - bip2*aip1)/det
         e22(1)=(-aip2*bi1 + bip2*ai1)/det
         zzexp(1)=0.d0
      endif
c
      if(ndv .ge. 2) then
         e21_eta=e21(1)/eta
         eta_e12=eta*e12(1)
c: Compute derivative w.r.t. k of exponential factor to be removed:
c: (see p.101-102 of Notebook):
         twok_eta=2.d0*xkh/eta
         dxi_dk=twok_eta/eta
         twok_eta3=dxi_dk/eta
         e11(2)=dxi_dk*(xi1*eta_e12 - e21_eta)
         e12(2)=twok_eta3*(e11(1) - e22(1))
         edif=xi1*e22(1) - xi2*e11(1)
         e21(2)=twok_eta*edif
         e22(2)=dxi_dk*(e21_eta - xi2*eta_e12)
         if(sg .ne. 0.d0) then
            delzp=zzexp1(2) - zzexp2(2)
            zzexp(2)=sg*dxi_dk*delzp
         else
            zzexp(2)=0.d0
         endif
         if(ndv .ge. 3) then
c: Compute derivative w.r.t. w (see p.108 of Notebook):
            deta_eta=1.d0/(1.5d0*w)
            dxi_fac=-deta_eta/etasq
            dxi_dw1=dxi_fac*(xk1sq + 2.d0*xkhsq)
            dxi_dw2=dxi_fac*(xk2sq + 2.d0*xkhsq)
            xi_dxi1=xi1*dxi_dw1
            xi_dxi2=xi2*dxi_dw2
            e11(3)=(xi_dxi1*eta_e12 - e21_eta*dxi_dw2)
            e12num=e11(1)*dxi_dw1 - e22(1)*dxi_dw2
            e12(3)=e12num/eta - e12(1)*deta_eta
            e21num=xi_dxi1*e22(1) - xi_dxi2*e11(1)
            e21(3)=e21num*eta + e21(1)*deta_eta
            e22(3)=(e21_eta*dxi_dw1 - xi_dxi2*eta_e12)
            if(sg .ne. 0.d0) then
               zzexp(3)=sg*(zzexp1(2)*dxi_dw1 - zzexp2(2)*dxi_dw2)
            else
               zzexp(3)=0.d0
            endif
         endif
         if(ndv .ge. 4) then
            e21fac=(e21(3) - 3.d0*e21(1)*deta_eta)/etasq
            e11(4)=twok_eta*(xi1*e12(3) + e12(1)*(dxi_dw1 - 
     .         xi1*deta_eta) - e21fac)
cc          e12(4)=(dxi_dk/eta)*(e11(3) - e22(3) - 3.d0*
cc   .         (e11(1) - e22(1))*deta_eta)
            e12(4)=twok_eta3*(e11(3) - e22(3)) - 
     .         3.d0*e12(2)*deta_eta
            e21(4)=twok_eta*(xi1*e22(3) + dxi_dw1*e22(1)
     .         - xi2*e11(3) - dxi_dw2*e11(1) - edif*deta_eta)
            e22(4)=-twok_eta*(xi2*e12(3) + e12(1)*(dxi_dw2 - 
     .         xi2*deta_eta) - e21fac)
            e5fac=dxi_dk*e12(1) - e21(2)/etasq
            e11(5)=twok_eta*(xi1*e12(2) + e5fac) + e11(2)/xkh
            e22(5)=-twok_eta*(xi2*e12(2) + e5fac) + e22(2)/xkh
            e12(5)=twok_eta3*(e11(2) - e22(2)) + e12(2)/xkh
cc          e21(5)=twok_eta*(xi1*e22(2) - xi2*e11(2)
cc   .         + dxi_dk*(e22(1) - e11(1))) + e21(2)/xkh
            e21(5)=twok_eta*(xi1*e22(2) - xi2*e11(2) - eta*e12(2))
     .         + e21(2)/xkh
            if(sg .ne. 0.d0) then
               zzexp(5)=sg*(2.d0*delzp/etasq + (zzexp1(3)-zzexp2(3))
     .            *dxi_dk*dxi_dk)
            else
               zzexp(5)=0.d0
            endif
            d2den=4.5d0*wsq*etasq
            d2xi_dw1=(xk1sq + 14.d0*xkhsq)/d2den
            d2xi_dw2=(xk2sq + 14.d0*xkhsq)/d2den
            xi_dxi1_w=xi1*d2xi_dw1 + dxi_dw1*dxi_dw1
            xi_dxi2_w=xi2*d2xi_dw2 + dxi_dw2*dxi_dw2
            eta_e12_w=eta*(e12(3) + deta_eta*e12(1))
            e21_eta_w=(-e21(1)*deta_eta + e21(3))/eta
            e11(6)=(xi_dxi1*eta_e12_w + xi_dxi1_w*eta_e12) - 
     .         (e21_eta*d2xi_dw2 + e21_eta_w*dxi_dw2)
            e22(6)=(e21_eta*d2xi_dw1 + e21_eta_w*dxi_dw1) - 
     .         (xi_dxi2*eta_e12_w + xi_dxi2_w*eta_e12)
            deta_eta_w=-deta_eta/w
            e12num_w=e11(1)*d2xi_dw1 + e11(3)*dxi_dw1 - 
     .         (e22(1)*d2xi_dw2 + e22(3)*dxi_dw2)
            e12(6)=(e12num_w - e12num*deta_eta)/eta - 
     .         (e12(1)*deta_eta_w + e12(3)*deta_eta)
            e21num_w=xi_dxi1*e22(3) + xi_dxi1_w*e22(1) - 
     .         (xi_dxi2*e11(3) + xi_dxi2_w*e11(1))
            e21(6)=e21num*eta*deta_eta + e21num_w*eta + 
     .         (e21(1)*deta_eta_w + e21(3)*deta_eta)
         endif
      endif
c
      return
      end
ccc
      subroutine rx_airy(argz,ai,bi,aip,bip,fexp)
c  arl normal mode - stickler   3/76, Modified by Levinson for ModeLab  5/93 
c: Modified by Westwood for ORCA 6/95 (derivatives of fexp).
c--      argz    arguement of the airy function
c--      ai      solution of the form  a(i) = ai * exp(-fexp)
c--      aip     ai prime of the form  a=(i) = aip * exp(-fexp)
c--      bi      solution of the form  b(i) = bi * exp(fexp)
c--      bip     bi prime of the form  b=(i) = bip * exp(fexp)
c--      fexp(1) power of e
c--      fexp(2) derivative of fexp w.r.t. argz
c--      fexp(3) 2nd derivative of fexp w.r.t. argz
c--  all coefficients were taken from nbs (big red)
      implicit none
      integer j,nstart
      real*8 ai,aip,argz,bi,bip,dzz,f,factor,fexp(3),g,h,l,n
     1   ,rspi,sq3,term1,term2,t1,t1i,t2,zs,zss,z2
     2   ,c1(6),c2(6),c3(6),c4(6),y1,y2,y3,y4,e2fx,fp
      data sq3/1.732050807568877d0/, rspi/.5641895835477563d0/
      data c1/ 1513.43953784973d0, 54.73651753365d0, 3.01523703552023d0,
     1 .283034838766718d0, .550974151234568d-01, .347222222222222d-01/
      data c2/ 274.026860098044d0,12.0840675164303d0,.851882583060333d0,
     1 .110951699674211d0, .355259773662552d-01, .694444444444444d-01/
      data c3/
     1 -1605.28176154639d0,-58.8223830649494d0,-3.30764827632364d0,
     2 -.321728616094393d0,-.670090663580247d-01,-.486111111111111d-01/
      data c4/
     1 -293.589241022865d0,-13.1945274712417d0,-.96065471819526d0,
     2 -.132692204646776d0,-.470357510288066d-01,-.972222222222222d-01/
c--
      fexp(1)=0.d0
      fexp(2)=0.d0
      fexp(3)=0.d0
      dzz  = abs(argz)
      if(dzz .le. 5.d0) then
c: power series expansion
         term1 = .355028053887817d0
         term2 = .258819403792807d0
         f     = term1
         l     = term2
         term2 = argz*term2
         g     = term2
         z2    = argz*argz
         h     = 0.d0
         n     = 1.d0
10       n     = n+1.d0
         term1 = term1*z2/n
         h     = h+term1
         n     = n+1.d0
         term1 = term1*argz/n
         fp    = f
         f     = f+term1
         term2 = term2*z2/n
         l     = l+term2
         n     = n+1.d0
         term2 = term2*argz/n
         g     = g+term2
         if(fp .ne. f) goto 10
         ai=f-g
         aip=h-l
         bi=sq3*(f+g)
         bip=sq3*(h+l)
c        if (argz .ge.-0.5069919730792688d0) then
c           factor=bi/3.775395359797091d-01
c           fexp=log(factor)
c           ai=ai*factor
c           aip=aip*factor
c           bip=bip/factor
c           bi= 3.775395359797091d-01
c        end if
         if(argz .gt. 0.d0) then
            zs=sqrt(argz)
            fexp(1)=zs*zs*zs/1.5d0
            fexp(2)=zs
            fexp(3)=0.5d0/zs
            factor=exp(fexp(1))
            ai=ai*factor
            aip=aip*factor
            bip=bip/factor
            bi=bi/factor
         endif
      else
c: phase amplitude expansions
         zs=sqrt(dzz)
         zss=sqrt(zs)
         t1=1.5d0/(zs*zs*zs)
         t2=t1*t1
         if(argz .lt. 0.d0) t2=-t2
         nstart=1
         if(dzz .gt. 12.d0) nstart=4
         y1=0.d0
         y2=0.d0
         y3=0.d0
         y4=0.d0
         do j=nstart,6
            y1=(y1+c1(j))*t2
            y2=(y2+c2(j))*t2
            y3=(y3+c3(j))*t2
            y4=(y4+c4(j))*t2
         enddo
         y1=y1+1.d0
         y2=y2+1.d0
         y3=y3+1.d0
         y4=y4+1.d0
         t1i=1.d0/t1
         y1=y1*rspi/zss
         y2=y2*t1i
         y3=y3*rspi*zss
         y4=y4*t1i
c
         if(argz .le. 0.d0) then
c: negative z
            y2=y2-0.7853981633974483d0
            y4=y4+3.9269908169872415d0
            ai=y1*cos(y2)
            bi=-y1*sin(y2)
            aip=y3*cos(y4)
            bip=-y3*sin(y4)
         else
c: positive z
            fexp(1)=y2
            y4=exp(y4-fexp(1))
c
            if (argz.lt.10.d0) then
               e2fx=exp(-2.d0*fexp(1))
               bi=y1-y1*e2fx
               bip=y3*(y4-e2fx/y4)
            else
               bi=y1
               bip=y3*y4
            endif
            ai=y1*.5d0
            aip=-y3/y4*.5d0
c: new: Change exponential factor to (2/3) z^(3/2):
            factor=exp(t1i - fexp(1))
            ai=ai*factor
            aip=aip*factor
            bip=bip/factor
            bi=bi/factor
            fexp(1)=t1i
            fexp(2)=zs
            fexp(3)=0.5d0/zs
         endif
      endif
c
      return
      end
ccc
      subroutine rx_airy_sm(argz,ai,bi,aip,bip,fexp)
c  arl normal mode - stickler   3/76, Modified by Levinson for ModeLab  5/93 
c: Modified by Westwood for ORCA 6/95 (derivatives of fexp).
c: Modified by Westwood 9/1/95 to make second derivative of f smooth at 0.
c--      argz    arguement of the airy function
c--      ai      solution of the form  a(i) = ai * exp(-fexp)
c--      aip     ai prime of the form  a=(i) = aip * exp(-fexp)
c--      bi      solution of the form  b(i) = bi * exp(fexp)
c--      bip     bi prime of the form  b=(i) = bip * exp(fexp)
c--      fexp(1) power of e
c--      fexp(2) derivative of fexp w.r.t. argz
c--      fexp(3) 2nd derivative of fexp w.r.t. argz
c--  all coefficients were taken from nbs (big red)
      implicit none
      integer j,nstart
      real*8 ai,aip,argz,bi,bip,dzz,f,factor,fexp(3),g,h,l,n
     1   ,rspi,sq3,term1,term2,t1,t1i,t2,zs,zss,z2
     2   ,c1(6),c2(6),c3(6),c4(6),y1,y2,y3,y4,e2fx,fp,
     .   c4x,c5x,c6x,c4xx,c5xx,c6xx,c4xxx,c5xxx,c6xxx,arg3,arg4,arg5
      data sq3/1.732050807568877d0/, rspi/.5641895835477563d0/
      data c1/ 1513.43953784973d0, 54.73651753365d0, 3.01523703552023d0,
     1 .283034838766718d0, .550974151234568d-01, .347222222222222d-01/
      data c2/ 274.026860098044d0,12.0840675164303d0,.851882583060333d0,
     1 .110951699674211d0, .355259773662552d-01, .694444444444444d-01/
      data c3/
     1 -1605.28176154639d0,-58.8223830649494d0,-3.30764827632364d0,
     2 -.321728616094393d0,-.670090663580247d-01,-.486111111111111d-01/
      data c4/
     1 -293.589241022865d0,-13.1945274712417d0,-.96065471819526d0,
     2 -.132692204646776d0,-.470357510288066d-01,-.972222222222222d-01/
      data c4x/1.03119738923038d0/,c5x/-0.61871843353823d0/,
     .   c6x/0.11048543456040d0/,c4xx/3.09359216769115d0/,
     .   c5xx/-2.47487373415292d0/,c6xx/0.55242717280199d0/,
     .   c4xxx/6.18718433538228d0/,c5xxx/-7.42462120245876d0/,
     .   c6xxx/2.20970869120800d0/
c--
      fexp(1)=0.d0
      fexp(2)=0.d0
      fexp(3)=0.d0
      dzz  = abs(argz)
      if(dzz .le. 5.d0) then
c: power series expansion
         term1 = .355028053887817d0
         term2 = .258819403792807d0
         f     = term1
         l     = term2
         term2 = argz*term2
         g     = term2
         z2    = argz*argz
         h     = 0.d0
         n     = 1.d0
10       n     = n+1.d0
         term1 = term1*z2/n
         h     = h+term1
         n     = n+1.d0
         term1 = term1*argz/n
         fp    = f
         f     = f+term1
         term2 = term2*z2/n
         l     = l+term2
         n     = n+1.d0
         term2 = term2*argz/n
         g     = g+term2
         if(fp .ne. f) goto 10
         ai=f-g
         aip=h-l
         bi=sq3*(f+g)
         bip=sq3*(h+l)
c        if (argz .ge.-0.5069919730792688d0) then
c           factor=bi/3.775395359797091d-01
c           fexp=log(factor)
c           ai=ai*factor
c           aip=aip*factor
c           bip=bip/factor
c           bi= 3.775395359797091d-01
c        end if
         if(argz .gt. 0.d0) then
            if(argz .gt. 2.d0) then
               zs=sqrt(argz)
               fexp(1)=zs*zs*zs/1.5d0
               fexp(2)=zs
               fexp(3)=0.5d0/zs
            else
               arg3=z2*argz
               arg4=arg3*argz
               arg5=arg4*argz
               fexp(1)=c4x*arg3 + c5x*arg4 + c6x*arg5
               fexp(2)=c4xx*z2 + c5xx*arg3 + c6xx*arg4
               fexp(3)=c4xxx*argz + c5xxx*z2 + c6xxx*arg3
            endif
            factor=exp(fexp(1))
            ai=ai*factor
            aip=aip*factor
            bip=bip/factor
            bi=bi/factor
         endif
      else
c: phase amplitude expansions
         zs=sqrt(dzz)
         zss=sqrt(zs)
         t1=1.5d0/(zs*zs*zs)
         t2=t1*t1
         if(argz .lt. 0.d0) t2=-t2
         nstart=1
         if(dzz .gt. 12.d0) nstart=4
         y1=0.d0
         y2=0.d0
         y3=0.d0
         y4=0.d0
         do j=nstart,6
            y1=(y1+c1(j))*t2
            y2=(y2+c2(j))*t2
            y3=(y3+c3(j))*t2
            y4=(y4+c4(j))*t2
         enddo
         y1=y1+1.d0
         y2=y2+1.d0
         y3=y3+1.d0
         y4=y4+1.d0
         t1i=1.d0/t1
         y1=y1*rspi/zss
         y2=y2*t1i
         y3=y3*rspi*zss
         y4=y4*t1i
c
         if(argz .le. 0.d0) then
c: negative z
            y2=y2-0.7853981633974483d0
            y4=y4+3.9269908169872415d0
            ai=y1*cos(y2)
            bi=-y1*sin(y2)
            aip=y3*cos(y4)
            bip=-y3*sin(y4)
         else
c: positive z
            fexp(1)=y2
            y4=exp(y4-fexp(1))
c
            if (argz.lt.10.d0) then
               e2fx=exp(-2.d0*fexp(1))
               bi=y1-y1*e2fx
               bip=y3*(y4-e2fx/y4)
            else
               bi=y1
               bip=y3*y4
            endif
            ai=y1*.5d0
            aip=-y3/y4*.5d0
c: new: Change exponential factor to (2/3) z^(3/2):
            factor=exp(t1i - fexp(1))
            ai=ai*factor
            aip=aip*factor
            bip=bip/factor
            bi=bi/factor
            fexp(1)=t1i
            fexp(2)=zs
            fexp(3)=0.5d0/zs
         endif
      endif
c
      return
      end
ccc
      subroutine rx_rp_flay(e11,e12,e21,e22,gami1,gami2,
     .   rhorat,ndv,V,TC)
c
c: Finds the reflection coefficient just above a fluid-fluid interface
c: (characterized by density ratio rho1/rho2=rhorat) given the 
c: reflection coefficient at the bottom of layer 2 below and the 
c: propagator matrix E for the fluid layer 2.  gami1 and gami2 are i*(the
c: vertical wavenumbers) at the bottom of layer 1 and bottom of layer 2
c: respectively.  The transmission coefficient W at the bottom
c: of layer 2 is also computed.  Cases where the thickness of layer 2 is
c: zero, the density ratio is one, and layer 2 is a halfspace (R=0) are 
c: treated.
c
      implicit none
      integer*4 ndv,j
      real*8 e11(6),e12(6),e21(6),e22(6),rhorat
      complex*16 V(6),TC(6),gami1(6),gami2(6),game11(6),game12(6),
     .   e1m(6),e2m(6),j1(6),k1(6),N(6),D(6),M(6),X(6),densq
c
      game11(1)=gami2(1)*e11(1)
      game12(1)=gami2(1)*e12(1)
      e1m(1)=game11(1) - e21(1)
      e2m(1)=game12(1) - e22(1)
      j1(1)=gami1(1)*e2m(1)
      k1(1)=e1m(1)
c
      N(1)=j1(1) + k1(1)
      D(1)=j1(1) - k1(1)
      if(D(1) .eq. (0.d0,0.d0)) then
         V(j)=(0.d0,0.d0) 
         TC(j)=(0.d0,0.d0) 
         do j=2,ndv
            V(j)=(1.d0,0.d0) 
            TC(j)=(1.d0,0.d0) 
         enddo
cc       if(iidiag .ge. 2) print *,'Info msg: D=0 in rx_rp_flay'
         return
      endif
      V(1)=N(1)/D(1)
c: Note that TC also has exp(-ze) factor associated with it:
      TC(1)=-2.d0*rhorat*gami1(1)/D(1)
c
c: If first derivatives desired, compute them now:
      densq=D(1)*D(1)
      if(ndv .ge. 2) then
         j=2
         game11(j)=gami2(1)*e11(j) + gami2(j)*e11(1)
         game12(j)=gami2(1)*e12(j) + gami2(j)*e12(1)
         e1m(j)=game11(j) - e21(j)
         e2m(j)=game12(j) - e22(j)
         j1(j)=gami1(1)*e2m(j) + gami1(j)*e2m(1)
         k1(j)=e1m(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=M(j)/densq
         X(j)=D(1)*gami1(j) - D(j)*gami1(1)
         TC(j)=-2.d0*rhorat*X(j)/densq
      endif
      if(ndv .ge. 3) then
         j=3
         game11(j)=gami2(1)*e11(j) + gami2(j)*e11(1)
         game12(j)=gami2(1)*e12(j) + gami2(j)*e12(1)
         e1m(j)=game11(j) - e21(j)
         e2m(j)=game12(j) - e22(j)
         j1(j)=gami1(1)*e2m(j) + gami1(j)*e2m(1)
         k1(j)=e1m(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=M(j)/densq
         X(j)=D(1)*gami1(j) - D(j)*gami1(1)
         TC(j)=-2.d0*rhorat*X(j)/densq
      endif
      if(ndv .gt. 3) then
         j=4
         game11(j)=gami2(1)*e11(j) + gami2(j)*e11(1)
     .      + gami2(2)*e11(3) + gami2(3)*e11(2)
         game12(j)=gami2(1)*e12(j) + gami2(j)*e12(1)
     .      + gami2(2)*e12(3) + gami2(3)*e12(2)
         e1m(j)=game11(j) - e21(j)
         e2m(j)=game12(j) - e22(j)
         j1(j)=gami1(1)*e2m(j) + gami1(j)*e2m(1)
     .      + gami1(2)*e2m(3) + gami1(3)*e2m(2)
         k1(j)=e1m(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) + D(3)*N(2) - D(j)*N(1) - D(2)*N(3)
         V(j)=(M(j) - 2.d0*M(2)*D(3)/D(1))/densq
c: Transmission coefficient:
         X(j)=D(1)*gami1(j) + D(3)*gami1(2) - D(j)*gami1(1) - 
     .      D(2)*gami1(3)
         TC(j)=-2.d0*rhorat*(X(j) - 2.d0*X(2)*D(3)/D(1))/densq
         j=5
         game11(j)=gami2(1)*e11(j) + gami2(j)*e11(1)
     .      + 2.d0*gami2(2)*e11(2)
         game12(j)=gami2(1)*e12(j) + gami2(j)*e12(1)
     .      + 2.d0*gami2(2)*e12(2)
         e1m(j)=game11(j) - e21(j)
         e2m(j)=game12(j) - e22(j)
         j1(j)=gami1(1)*e2m(j) + gami1(j)*e2m(1)
     .      + 2.d0*gami1(2)*e2m(2)
         k1(j)=e1m(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=(M(j) - 2.d0*M(2)*D(2)/D(1))/densq
c: Transmission coefficient:
         X(j)=D(1)*gami1(j) - D(j)*gami1(1)
         TC(j)=-2.d0*rhorat*(X(j) - 2.d0*X(2)*D(2)/D(1))/densq
         j=6
         game11(j)=gami2(1)*e11(j) + gami2(j)*e11(1)
     .      + 2.d0*gami2(3)*e11(3)
         game12(j)=gami2(1)*e12(j) + gami2(j)*e12(1)
     .      + 2.d0*gami2(3)*e12(3)
         e1m(j)=game11(j) - e21(j)
         e2m(j)=game12(j) - e22(j)
         j1(j)=gami1(1)*e2m(j) + gami1(j)*e2m(1)
     .      + 2.d0*gami1(3)*e2m(3)
         k1(j)=e1m(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=(M(j) - 2.d0*M(3)*D(3)/D(1))/densq
c: Transmission coefficient:
         X(j)=D(1)*gami1(j) - D(j)*gami1(1)
         TC(j)=-2.d0*rhorat*(X(j) - 2.d0*X(3)*D(3)/D(1))/densq
      endif
c
      return
      end
ccc
      subroutine rx_rp_vac(e11,e12,e21,e22,gami1,gami2,
     .   rho_prod,ndv,V,TC,sg,jlay,ii)
c
c: Finds the reflection coefficient given fluid propagator matrix E
c: and assuming a pressure release halfspace (see ORCA II, p.47).
c
      implicit none
      integer*4 ndv,j,jlay,ii
      real*8 e11(6),e12(6),e21(6),e22(6),rho_prod,sg
      complex*16 V(6),TC(6),gami1(6),gami2(6),j1(6),k1(6),N(6),
     .   D(6),M(6),X(6),Dw(6),densq,denwsq,zzero
      data zzero/(0.d0,0.d0)/
c
      if(jlay .eq. 2 .and. ii .eq. 1) then
c: Case where only pressure-release interface is above ref depth (V=-1, TC=1):
         V(1)=dcmplx(-1.d0,0.d0)
         TC(1)=dcmplx(1.d0,0.d0)
         do j=2,ndv
            V(j)=zzero
            TC(j)=zzero
         enddo
         return
      endif
      j1(1)=sg*gami1(1)*e12(1)
      k1(1)=e11(1)
      N(1)=j1(1) + k1(1)
      D(1)=j1(1) - k1(1)
      if(D(1) .eq. (0.d0,0.d0)) then
         print *,'D=0'
      endif
      V(1)=N(1)/D(1)
      Dw(1)=gami2(1)*D(1)
      TC(1)=-rho_prod*gami1(1)/Dw(1)
      densq=D(1)*D(1)
      denwsq=Dw(1)*Dw(1)
c
c: If first derivatives desired, compute them now:
      if(ndv .ge. 2) then
         j=2
         j1(j)=sg*(gami1(1)*e12(j) + gami1(j)*e12(1))
         k1(j)=e11(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=M(j)/densq
         Dw(j)=gami2(1)*D(j) + gami2(j)*D(1)
         X(j)=Dw(1)*gami1(j) - Dw(j)*gami1(1)
         TC(j)=-rho_prod*X(j)/denwsq
      endif
      if(ndv .ge. 3) then
         j=3
         j1(j)=sg*(gami1(1)*e12(j) + gami1(j)*e12(1))
         k1(j)=e11(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=M(j)/densq
         Dw(j)=gami2(1)*D(j) + gami2(j)*D(1)
         X(j)=Dw(1)*gami1(j) - Dw(j)*gami1(1)
         TC(j)=-rho_prod*X(j)/denwsq
      endif
c
c: If second derivative desired, compute them now:
      if(ndv .ge. 4) then
         j=4
         j1(j)=sg*(gami1(1)*e12(j) + gami1(j)*e12(1)
     .      + gami1(2)*e12(3) + gami1(3)*e12(2))
         k1(j)=e11(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) + D(3)*N(2) - D(j)*N(1) - D(2)*N(3)
         V(j)=(M(j) - 2.d0*M(2)*D(3)/D(1))/densq
c: Transmission coefficient:
         Dw(j)=gami2(1)*D(j) + gami2(j)*D(1) + gami2(2)*D(3) + 
     .      gami2(3)*D(2)
         X(j)=Dw(1)*gami1(j) + Dw(3)*gami1(2) - Dw(j)*gami1(1) - 
     .      Dw(2)*gami1(3)
         TC(j)=-rho_prod*(X(j) - 2.d0*X(2)*Dw(3)/Dw(1))/denwsq
         j=5
         j1(j)=sg*(gami1(1)*e12(j) + gami1(j)*e12(1)
     .      + 2.d0*gami1(2)*e12(2))
         k1(j)=e11(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=(M(j) - 2.d0*M(2)*D(2)/D(1))/densq
c: Transmission coefficient:
         Dw(j)=gami2(1)*D(j) + gami2(j)*D(1) + 2.d0*gami2(2)*D(2)
         X(j)=Dw(1)*gami1(j) - Dw(j)*gami1(1)
         TC(j)=-rho_prod*(X(j) - 2.d0*X(2)*Dw(2)/Dw(1))/denwsq
         j=6
         j1(j)=sg*(gami1(1)*e12(j) + gami1(j)*e12(1)
     .      + 2.d0*gami1(3)*e12(3))
         k1(j)=e11(j)
         N(j)=j1(j) + k1(j)
         D(j)=j1(j) - k1(j)
         M(j)=D(1)*N(j) - D(j)*N(1)
         V(j)=(M(j) - 2.d0*M(3)*D(3)/D(1))/densq
c: Transmission coefficient:
         Dw(j)=gami2(1)*D(j) + gami2(j)*D(1) + 2.d0*gami2(3)*D(3)
         X(j)=Dw(1)*gami1(j) - Dw(j)*gami1(1)
         TC(j)=-rho_prod*(X(j) - 2.d0*X(3)*Dw(3)/Dw(1))/denwsq
      endif
c
      return
      end
ccc
      subroutine num_deriv(f1,fp1,e1,f2,fp2,e2,dx,pct,var)
      implicit none
      character*3 var
      real*8 f1,fp1,e1,f2,fp2,e2,dx,pct,fp1x,fp2x,fpnx
      fp1x=fp1*exp(e1)
      fp2x=fp2*exp(e2)
      fpnx=(f2*exp(e2) - f1*exp(e1))/dx
      pct=(fpnx - fp1x)/(fp2x - fp1x)
      if(pct .lt. .48 .or. pct .gt. .52d0) then
         print *,'pct ',var,' ',sngl(100.*pct)
         print *,fp1x,fpnx,fp2x
      endif
      return
      end
