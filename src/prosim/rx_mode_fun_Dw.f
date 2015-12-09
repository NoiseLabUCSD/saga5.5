c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
      subroutine rx_mode_fun_Dw(jm,
     .   jlay_ref,ii_ref,vg,jdf)
c
c: Finds mode functions for the current eigenvalue.
      implicit none
      include 'Parms_com'
      include 'i_o_1b_com'
      common /ny_com/ narray(NVRMAX),act_sect
      integer*4 narray,act_sect
      common /out_com3/ nzsr, nlay,NSCTOR,dmmy3
      integer*4 nzsr,nlay,NSCTOR,dmmy3
      include 'i_o_opt_com'
      common /eig2_com1/ phi_bb(NVRMAX,NM_MAX,NFCMAX)
      real*4 phi_bb
      common /eig2_com3/ dphiz_bb(NVRMAX,NM_MAX,NFCMAX)
      real*4 dphiz_bb
      include 'gen1_com'
      include 'gen2_com'
      integer*4 jlay,j,jz,jx0,jx,jm,jsr,jzx,
     .          jlay_ref,ii1,ii2,ii_ref,jjx
      integer*4 kzsr,ncount,displcm,ny,jdf,dplcm,jd
      real*8 gamma,igamma,gamsq,
     .   gamz,cosgamz,singamz,A1,B1,A1_exp,B1_exp,
     .   A1tot,B1tot,xi1,ai1,bi1,aip1,bip1,zzexp1(3),
     .   sg,sg_rho,expmin,rzero,phsum,
     .   dphsum,Dgamma,Digamma,DA1,DB1,Dgamz,DA1_exp,DB1_exp,
     .   Dxi(2),Dxi1,Deta,phi_exp1,phi_exp2,xix(2),
     .   Dai1,Dbi1,Daip1,Dbip1,
     .   DA1tot,DB1tot,Dphsum_w,
     .   expfaca,expfacb,expA,expB
      real*8 phixs(NVRMAX),Dphixs_w(NVRMAX),
     .       vg,pie
      real*8 phixa(NVRMAX)

c Variables introduced for the Gaussian Beam source
      complex*16 zzero
c: Set expmin to ln of min value of phi,dphi before being set to zero 
c: (expmin = -23 for phi = 1.d-10 here):
      data expmin/-23.d0/,zzero/(0.d0,0.d0)/,rzero/0.d0/
c
c: Find mode group velocity and attenuation:
      call rx_vg_atten
c
      pie=3.14159265358979d0
c
      displcm=0
      dplcm=0
      if(act_sect .eq. 1) then
        ncount=narray(act_sect)
      else
        ncount=narray(act_sect)+narray(act_sect-1)
      endif
      if((ncount .eq. 0).or.(NSCTOR .eq. 1)) ncount=1
      do ny=1,ncount

c: Find mode function values at desired depths:
      do jlay=1,nlwz(ny)
         j=jlwz(jlay+displcm)
c fbv         j=jlwz(jlay) is this the problem?
c ORCA DOCET?         jx0=jzmx(j+displcm)-1
         jx0=jzmx(j+dplcm)-1
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
c ORCA DOCET?            do jz=1,nzmx(j+displcm)
            do jz=1,nzmx(j+dplcm)
               jzx=jx0 + jz
               phixs(jzx+displcm)=zzero
               Dphixs_w(jzx+displcm)=rzero
            enddo
         elseif(iabs(isp(j)) .eq. 1) then
c            gamsq=dble(xksq(1,j)) - dble(xkhsq)
            gamsq=xksq(1,j) - xkhsq
            if(gamsq .gt. 0.d0) then
c: Isospeed propagating (see p. 41, ORCA II):
c: (zmx=0 at top of layer in rx_zmx_init)
c: Correction by Westwood (1 instead of ii1 index) 16/10-97
               gamma=sqrt(gamsq)
               A1_exp=dexp(phi_exp1)
               A1=dble(philay(1,j))*A1_exp
               B1=dble(dphilay(1,j))*A1_exp/gamma
c               Dgamma=(dble(xksq(1,j))/w - dble(xkh)/vg)/gamma
               Dgamma=(xksq(1,j)/w - xkh/vg)/gamma
               DA1=Dphi_w(1,j)*A1_exp
               DB1=(gamma*Ddphi_w(1,j)-Dgamma*dble(dphilay(1,j)))
     .            *A1_exp/gamsq
c ORCA DOCET?               do jz=1,nzmx(j+displcm)
               do jz=1,nzmx(j+dplcm)
                  jzx=jx0 + jz + displcm
                  gamz=gamma*zmx(jzx)
                  Dgamz=Dgamma*zmx(jzx)
                  cosgamz=cos(gamz)
                  singamz=sin(gamz)
                  phsum=A1*cosgamz + B1*singamz
                  dphsum=-A1*singamz + B1*cosgamz
                  phixs(jzx)=phsum
                  Dphixs_w(jzx)=Dgamz*dphsum + DA1*cosgamz + 
     .              DB1*singamz
               enddo
            elseif(j .ne. nlay) then
c: Isospeed evanescent (see p. 41, ORCA II):
               igamma=sqrt(-gamsq)
c               Digamma=-(dble(xksq(1,j))/w - dble(xkh)/vg)/igamma
               Digamma=-(xksq(1,j)/w - xkh/vg)/igamma
               call rx_a1b1_iso_Dw(ii1,ii2,philay(1,j),dphilay(1,j),
     .            Aplay(1,j),zetalay(1,1,j),igamma,A1,B1,expA,expB,
     .            gamsq,Dphi_w(1,j),Ddphi_w(1,j),Digamma,DA1,DB1)
c ORCA DOCET?               do jz=1,nzmx(j+displcm)
               do jz=1,nzmx(j+dplcm)
                  jzx=jx0 + jz + displcm
                  gamz=igamma*zmx(jzx)
                  Dgamz=Digamma*zmx(jzx)
                  expfaca=dexp(gamz + expA)
                  expfacb=dexp(-gamz + expB)
                  A1_exp=A1*expfaca
                  B1_exp=B1*expfacb
                  DA1_exp=DA1*expfaca
                  DB1_exp=DB1*expfacb
                  phixs(jzx)=A1_exp + B1_exp
                  Dphixs_w(jzx)=DA1_exp + DB1_exp
               enddo
            else
c: Isospeed evanescent for halfspace, where A1=0 (see p. 41, ORCA II):
               igamma=sqrt(-gamsq)
c               Digamma=-(dble(xksq(1,j))/w - dble(xkh)/vg)/igamma
               Digamma=-(xksq(1,j)/w - xkh/vg)/igamma
               B1=dble(philay(ii1,j))
               DB1=Dphi_w(ii1,j)
c ORCA DOCET?               do jz=1,nzmx(j+displcm)
               do jz=1,nzmx(j+dplcm)
                  jzx=jx0 + jz + displcm
                  gamz=igamma*zmx(jzx)
                  Dgamz=Digamma*zmx(jzx)
                  expfacb=dexp(-gamz + phi_exp1)
                  phixs(jzx)=B1*expfacb
                  Dphixs_w(jzx)=DB1*expfacb
               enddo
            endif
         else
c: Airy layer (see p. 42, ORCA II):
            xix(1)=(xkhsq - xksq(1,j))/etasq(j)
            xix(2)=(xkhsq - xksq(2,j))/etasq(j)
c: See ORCA II, p.65:
            Dxi(1)=2.d0*((xksq(1,j)+2.d0*xkhsq)/(-3.d0*w) + 
     .         xkh/vg)/etasq(j)
            Dxi(2)=2.d0*((xksq(2,j)+2.d0*xkhsq)/(-3.d0*w) + 
     .         xkh/vg)/etasq(j)
            Deta=dble(eta(j))/(1.5d0*w)
c
            call rx_a1b1_airy_Dw(ii1,ii2,philay(1,j),dphilay(1,j),
     .         Aplay(1,j),eta(j),ailay(1,1,1,j),bilay(1,1,1,j),
     .         zetalay(1,1,j),A1,B1,expA,expB,pie,Deta,xix,Dxi,
     .         etasq(j),Dphi_w(1,j),Ddphi_w(1,j),DA1,DB1)
c     .         Aplay(1,j),eta(j),ailay(1,1,1,j),bilay(1,1,1,j),
c     .         Aplay(1,j),eta(j),ailay(1,j),bilay(1,j),
c
c ORCA DOCET?            do jz=1,nzmx(j+displcm)
            do jz=1,nzmx(j+dplcm)
               jzx=jx0 + jz + displcm
               xi1=xix(1) - eta(j)*zmx(jzx)
               Dxi1=Dxi(1) - Deta*zmx(jzx)
               call rx_airy_sm(xi1,ai1,bi1,aip1,bip1,zzexp1)
               Dai1=aip1*Dxi1
               Daip1=(xi1*ai1)*Dxi1
               Dbi1=bip1*Dxi1
               Dbip1=(xi1*bi1)*Dxi1
               expfaca=dexp(expA - zzexp1(1))
               expfacb=dexp(expB + zzexp1(1))
               A1tot=A1*expfaca
               B1tot=B1*expfacb
               DA1tot=DA1*expfaca
               DB1tot=DB1*expfacb
               phsum=A1tot*aip1 + B1tot*bip1
               phixs(jzx)=A1tot*ai1 + B1tot*bi1
               Dphixs_w(jzx)=A1tot*Dai1 + B1tot*Dbi1 + 
     .            DA1tot*ai1 + DB1tot*bi1
               Dphsum_w=A1tot*Daip1 + B1tot*Dbip1 + 
     .            DA1tot*aip1 + DB1tot*bip1 
            enddo
         endif
      enddo
      displcm=displcm+nrecusr+1
      dplcm=dplcm+nlay
      enddo
c
c: Make dphi at surface positive as a convention:
      sg=1.d0
      if(dble(dphilay(1,2)) .lt. 0.d0) sg=-1.d0
c: Convert from phixs to phi using mx_m pointers computed in zsr_init:
      displcm=0
      dplcm=0
      kzsr=ncount*(nrecusr+1)
      do jd=1,ncount
       do jsr=1,nrecusr+1
         jx=mx_m(jsr+displcm)
c Added for compliance with new ORCA
         jjx=zsr_indx(jx+displcm)
c End of addendum
         sg_rho=sg*rho_sr(jsr+displcm)
c: INCLUDE RHO in PHI,DPHI,PSI,DPSI:
         phixa(jjx+displcm)=phixs(jx+displcm)
         phi_bb(jjx+displcm,jm,jdf)=sg_rho*phixs(jx+displcm)
c fbv         phi_bb(jsr+displcm,jm,jdf)=sg_rho*phixs(jx+displcm)
c fbv         phi_bb(jsr+displcm,jm,jdf)=sg_rho*phixs(jjx+displcm)
         dphiz_bb(jjx+displcm,jm,jdf)=sg_rho*Dphixs_w(jx+displcm)
c fbv         dphiz_bb(jsr+displcm,jm,jdf)=sg_rho*Dphixs_w(jx+displcm)
c fbv         dphiz_bb(jsr+displcm,jm,jdf)=sg_rho*Dphixs_w(jjx+displcm)
       enddo
         displcm=displcm+nrecusr+1
         dplcm=dplcm+nrecusr
      enddo
c
      return
      end
c
      subroutine rx_vg_atten
c
c: Finds modal group velocity and attenuation.
c
      implicit none
      include 'Parms_com'
      include 'i_o_1b_com'
      common /out2_com/ kn(0:NM_MAX)
      complex*16 kn
      include 'i_o_2_com'
      include 'gen1_com'
      integer*4 jlay,j,ii
      real*8 phi_duct,dphi_duct,gamsq_duct,mfun_mag,mfun_max,Apfac,
     .   vg,dn,phi_t,phi_b,dphi_t,dphi_b,I1,I2
c
c: Find duct in which mode has its maximum value:
      mfun_max=-1.d0
      nzref(nmode)=1
      do j=1,ndrx
         jlay=jval(1,j)
         ii=jval(2,j)
         phi_duct=philay(ii,jlay)
         dphi_duct=dphilay(ii,jlay)
c          gamsq_duct=max(0.d0,dble(xksq(ii,jlay)) - dble(xkhsq))
         gamsq_duct=max(0.d0,xksq(ii,jlay) - xkhsq)
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
c fbv      N_n=0.d0
      vg=0.d0
      dn=0.d0
      do j=2,nlay
         phi_t=philay(1,j)
         phi_b=philay(2,j)
         dphi_t=dphilay(1,j)
         dphi_b=dphilay(2,j)
         call rx_int(phi_t,phi_b,dphi_t,dphi_b,Aplay(1,j),Aplay(2,j),
     .      xilay(1,j),xilay(2,j),geo(1,3,j),geo(2,3,j),eta(j),etasq(j),
     .      h(j),isp(j),I1,I2)
c fbv         N_n=N_n + I1
         vg=vg + (am(j)*I1 + bm(j)*I2)
         dn=dn + (betm(j)*I1 + gm(j)*I2)
      enddo
      vg=kn(nmode)/(w*vg)
      dn=(w/kn(nmode))*dn
c: Enter attenuation as imaginary part of kn:
      kn(nmode)=dcmplx(dble(kn(nmode)),dn)
c
      return
      end
ccc
      subroutine rx_int(phi1,phi2,dphi1,dphi2,Ap1,Ap2,xi1,xi2,
     .   rho1,rho2,eta,etasq,h,iso,I1,I2)
c
      implicit none
      integer*4 iso
      real*8 phi1,phi2,dphi1,dphi2,Ap1,Ap2,xi1,xi2,rho1,rho2,eta,
     .   etasq,h,expmin,I1,I2,hsq2,term1,term2,term3,term4,
     .   dzero,ex1,ex2,phi1sq,phi2sq,dphi1sq,dphi2sq,mix2,done
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
      if(iso .eq. 0) then
         if(Ap2 .gt. expmin) then
            if(Ap2 .ne. 0.d0) ex2=dexp(2.d0*Ap2)
            phi2sq=phi2*phi2
            dphi2sq=dphi2*dphi2
            term2=rho2*(dphi2sq/etasq - xi2*phi2sq)
               term4=rho2*(dphi2sq*(3.d0*h + 2.d0*xi2/eta) -
     .            phi2*dphi2 - phi2sq*xi2*(3.d0*etasq*h + 
     .            2.d0*eta*xi2))
            if(Ap1 .gt. expmin) then
               if(Ap1 .ne. 0.d0) ex1=dexp(2.d0*Ap1)
               phi1sq=phi1*phi1
               dphi1sq=dphi1*dphi1
               term1=rho1*(dphi1sq/etasq - xi1*phi1sq)
               I1=(term2*ex2 - term1*ex1)/eta
                  term3=rho1*(2.d0*xi1*dphi1sq/eta -
     .               phi1*dphi1 - 2.d0*eta*xi1*xi1*phi1sq)
                  I2=(term4*ex2 - term3*ex1)/(3.d0*eta*etasq)
            else
               I1=term2*ex2/eta
                  I2=term4*ex2/(3.d0*eta*etasq)
            endif
         elseif(Ap1 .gt. expmin) then
            if(Ap1 .ne. 0.d0) ex1=dexp(2.d0*Ap1)
            phi1sq=phi1*phi1
            dphi1sq=dphi1*dphi1
            term1=rho1*(dphi1sq/etasq - xi1*phi1sq)
            I1=-term1*ex1/eta
               term3=rho1*(2.d0*xi1*dphi1sq/eta -
     .            phi1*dphi1 - 2.d0*eta*xi1*xi1*phi1sq)
               I2=-term3*ex1/(3.d0*eta*etasq)
         endif
      else
         if(Ap2 .gt. expmin) then
            if(Ap2 .ne. 0.d0) ex2=dexp(2.d0*Ap2)
            phi2sq=phi2*phi2
            dphi2sq=dphi2*dphi2
            mix2=phi2*dphi2
            term2=rho2*(phi2sq*h + (dphi2sq*h - mix2)/xi1)
               hsq2=2.d0*h*h
               term4=rho2*(dphi2sq*(hsq2 - 1.d0/xi1) + 
     .            phi2sq*(hsq2*xi1 + 1.d0) - 4.d0*mix2*h)
            if(Ap1 .gt. expmin) then
               if(Ap1 .ne. dzero) ex1=dexp(2.d0*Ap1)
               term1=rho1*(-phi1*dphi1/xi1)
               I1=0.5d0*(term2*ex2 - term1*ex1)
                  term3=rho1*(-dphi1*dphi1/xi1 + phi1*phi1)
                  I2=(term4*ex2 - term3*ex1)/(8.d0*xi1)
            else
               I1=0.5d0*term2*ex2
                  I2=term4*ex2/(8.d0*xi1)
            endif
         elseif(Ap1 .gt. expmin) then
            if(Ap1 .ne. dzero) ex1=dexp(2.d0*Ap1)
            term1=rho1*(-phi1*dphi1/xi1)
            I1=-0.5d0*term1*ex1
               term3=rho1*(-dphi1*dphi1/xi1 + phi1*phi1)
               I2=-term3*ex1/(8.d0*xi1)
         endif
      endif
c
      return
      end
ccc
      function rx_round(term1,term2,deninv,iir)
c
c: Checks to see if cancellation occurs in term1+term2.
c: Returns round_check=(term1 + term2)*deninv if not.
      implicit none
      integer*4 iir
      real*8 rx_round,term1,term2,deninv,numerator
      real*8 magsum,ratio
c
      numerator=term1 + term2
      magsum=max(dabs(term1),dabs(term2))
      ratio=dabs(numerator)/magsum
      if(ratio .gt. 1.d-10) then
         rx_round=numerator*deninv
         iir=0
      else
         iir=1
         rx_round=0.d0
cc       print *,'Round to zero: ',term1,term2,numerator,deninv
      endif
c
      return
      end
ccc
      subroutine rx_a1b1_iso_Dw(ii1,ii2,philay,dphilay,Aplay,zetalay,
     .   igamma,A1,B1,expA,expB,gamsq,Dphi_w,Ddphi_w,Digamma,DA1,DB1)
     .            
c
      implicit none
      integer*4 ii1,ii2,ii
      complex*16 philay(2),dphilay(2),zetalay(2)
      real*8 Aplay(2),igamma,A1,B1,expA,expB,term1,term2,rx_round,
     .   rata,ratb,gamsq,Dphi_w(2),Ddphi_w(2),Digamma,DA1,DB1,Dterm2
c
      ii=ii1
      term1=dble(philay(ii))
      term2=dble(dphilay(ii))/igamma
      A1=rx_round(term1,term2,0.5d0,rata)
      expA=Aplay(ii) - dble(zetalay(ii))
      Dterm2=(igamma*Ddphi_w(ii) - Digamma*dble(dphilay(ii)))/gamsq
      DA1=(Dphi_w(ii) + Dterm2)*0.5d0
c
      ii=ii2
      term1=dble(philay(ii))
      term2=dble(dphilay(ii))/igamma
c: Correction by Westwood (-term2 instead of term2) 16/10-97
      B1=rx_round(term1,-term2,0.5d0,ratb)
      expB=Aplay(ii) + dble(zetalay(ii))
      Dterm2=(igamma*Ddphi_w(ii) - Digamma*dble(dphilay(ii)))/gamsq
      DB1=(Dphi_w(ii) - Dterm2)*0.5d0
c
      return
      end
ccc
      subroutine rx_a1b1_airy_Dw(ii1,ii2,philay,dphilay,Aplay,
     .   eta,ailay,bilay,zetalay,A1,B1,expA,expB,pie,Deta,xi,Dxi,
     .   etasq,Dphi_w,Ddphi_w,DA1,DB1)
c
      implicit none
      integer*4 ii1,ii2,ii
      complex*16 philay(2),dphilay(2),zetalay(2)
      real*8 ailay(2,2),bilay(2,2)
      real*8 Aplay(2),eta,A1,B1,expA,expB,pie,term1,term2,rx_round,
     .   dphi_eta,rata,ratb,Deta,xi(2),Dxi(2),etasq,Dphi_w(2),
     .   Ddphi_w(2),DA1,DB1,Ai0,Ai0p,Bi0,Bi0p,DAi0,DAi0p,Dbi0,DBi0p,
     .   Ddphi_eta
c
c: ALWAYS FIT AI COEFFICIENT A1 AT INTERFACE CLOSEST TO DUCT IN WHICH MODE
c: IS TRAPPED, AND BI COEFFICIENT B1 AT INTERFACE FARTHEST FROM DUCT.
c
      ii=ii1
c      Bi0=dble(bilay(1,ii))
c      Bi0p=dble(bilay(2,ii))
      Bi0=bilay(1,ii)
      Bi0p=bilay(2,ii)
      dphi_eta=dble(dphilay(ii))/eta
      term1=Bi0p*dble(philay(ii))
      term2=Bi0*dphi_eta
      A1=rx_round(term1,term2,pie,rata)
      expA=Aplay(ii) + dble(zetalay(ii))
      Ddphi_eta=(eta*Ddphi_w(ii) - Deta*dble(dphilay(ii)))/etasq
      DBi0=Bi0p*Dxi(ii)
      DBi0p=(xi(ii)*Bi0)*Dxi(ii)
      DA1=pie*(Bi0p*Dphi_w(ii) + DBi0p*
     .   dble(philay(ii)) + Bi0*Ddphi_eta + DBi0*dphi_eta)
c
      ii=ii2
      dphi_eta=dble(dphilay(ii))/eta
c      Ai0=dble(ailay(1,ii))
c      Ai0p=dble(ailay(2,ii))
      Ai0=ailay(1,ii)
      Ai0p=ailay(2,ii)
      term1=Ai0p*dble(philay(ii))
      term2=Ai0*dphi_eta
      B1=rx_round(term1,term2,-pie,ratb)
      expB=Aplay(ii) - dble(zetalay(ii))
c
      Ddphi_eta=(eta*Ddphi_w(ii) - Deta*dble(dphilay(ii)))/etasq
      DAi0=Ai0p*Dxi(ii)
      DAi0p=(xi(ii)*Ai0)*Dxi(ii)
      DB1=-pie*(Ai0p*Dphi_w(ii) + DAi0p*dble(philay(ii)) + 
     .   Ai0*Ddphi_eta + DAi0*dphi_eta)
c
      return
      end
