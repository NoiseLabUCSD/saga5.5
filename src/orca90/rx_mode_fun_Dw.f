      subroutine rx_mode_fun_Dw(jm,phiz,dphiz,Dphiz_wx,Ddphiz_wx,
     .   jlay_ref,ii_ref,vg)
c
c: Finds mode functions for the current eigenvalue.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jlay_ref,ii_ref,jlay,j,ii1,ii2,jz,jx0,jx,jm,jsr,jzx
      real*8 gamma,igamma,gamsq,gamz,
     .   cosgamz,singamz,phi_exp1,phi_exp2,A1,B1,A1_exp,B1_exp,A1tot,
     .   B1tot,xix(2),xi1,ai1,bi1,aip1,bip1,zzexp1(3),sg,
     .   sg_rho,expmin,rzero,phsum,dphsum,Dgamma,
     .   Digamma,DA1,DB1,Dgamz,DA1_exp,DB1_exp,Dxi(2),Dxi1,Deta,
     .   DAi1,DAip1,DBi1,DBip1,DA1tot,DB1tot,Dphsum_w,
     .   expfaca,expfacb,expA,expB
      real*8 Dphiz_wx(nzsr,jm),Ddphiz_wx(nzsr,jm),vg,phiq,dphiq
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
               Dphix_w(jzx)=rzero
               Ddphix_w(jzx)=rzero
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
               Dgamma=(dreal(xksq(1,j))/w - dreal(xkh)/vg)/gamma
               DA1=Dphi_w(1,j)*A1_exp
               DB1=(gamma*Ddphi_w(1,j)-Dgamma*dreal(dphilay(1,j)))
     .            *A1_exp/gamsq
               do jz=1,nzmx(j)
                  jzx=jx0 + jz
                  gamz=gamma*zmx(jzx)
                  Dgamz=Dgamma*zmx(jzx)
                  cosgamz=cos(gamz)
                  singamz=sin(gamz)
                  phsum=A1*cosgamz + B1*singamz
                  dphsum=-A1*singamz + B1*cosgamz
                  phix(jzx)=phsum
                  dphix(jzx)=gamma*dphsum
                  Dphix_w(jzx)=Dgamz*dphsum + DA1*cosgamz + DB1*singamz
                  Ddphix_w(jzx)=-gamma*(phsum*Dgamz + DA1*singamz -
     .               DB1*cosgamz) + Dgamma*dphsum
               enddo
            elseif(j .ne. nlay) then
c: Isospeed evanescent (see p. 41, ORCA II):
               igamma=sqrt(-gamsq)
               Digamma=-(dreal(xksq(1,j))/w - dreal(xkh)/vg)/igamma
               call rx_a1b1_iso_Dw(ii1,ii2,philay(1,j),dphilay(1,j),
     .            Aplay(1,j),zetalay(1,1,j),igamma,A1,B1,expA,expB,
     .            gamsq,Dphi_w(1,j),Ddphi_w(1,j),Digamma,DA1,DB1)
               do jz=1,nzmx(j)
                  jzx=jx0 + jz
                  gamz=igamma*zmx(jzx)
                  Dgamz=Digamma*zmx(jzx)
                  expfaca=dexp(gamz + expA)
                  expfacb=dexp(-gamz + expB)
                  A1_exp=A1*expfaca
                  B1_exp=B1*expfacb
                  DA1_exp=DA1*expfaca
                  DB1_exp=DB1*expfacb
                  phix(jzx)=A1_exp + B1_exp
                  dphix(jzx)=igamma*(A1_exp - B1_exp)
                  Dphix_w(jzx)=DA1_exp + DB1_exp
                  Ddphix_w(jzx)=igamma*(DA1_exp - DB1_exp) + 
     .               Digamma*(A1_exp - B1_exp)
               enddo
            else
c: Isospeed evanescent for halfspace, where A1=0 (see p. 41, ORCA II):
               igamma=sqrt(-gamsq)
               Digamma=-(dreal(xksq(1,j))/w - dreal(xkh)/vg)/igamma
               B1=dreal(philay(1,j))
               DB1=Dphi_w(1,j)
               do jz=1,nzmx(j)
                  jzx=jx0 + jz
                  gamz=igamma*zmx(jzx)
                  Dgamz=Digamma*zmx(jzx)
                  expfacb=dexp(-gamz + phi_exp1)
                  phix(jzx)=B1*expfacb
                  dphix(jzx)=-igamma*phix(jzx)
                  Dphix_w(jzx)=DB1*expfacb
                  Ddphix_w(jzx)=-igamma*Dphix_w(jzx) - Digamma*phix(jzx)
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
            Deta=dreal(eta(j))/(1.5d0*w)
c
            call rx_a1b1_airy_Dw(ii1,ii2,philay(1,j),dphilay(1,j),
     .         Aplay(1,j),eta(j),ailay(1,1,1,j),bilay(1,1,1,j),
     .         zetalay(1,1,j),A1,B1,expA,expB,pie,Deta,xix,Dxi,etasq(j),
     .         Dphi_w(1,j),Ddphi_w(1,j),DA1,DB1)
co            call rx_a1b1_airy_Dw(ii1,ii2,philay(1,j),dphilay(1,j),
cp     .         Aplay(1,j),dreal(eta(j)),ailay(1,1,1,j),bilay(1,1,1,j),
co     .         zetalay(1,1,j),A1,B1,expA,expB,pie,Deta,xix,Dxi,
co     .		 dreal(etasq(j)),
co     .         Dphi_w(1,j),Ddphi_w(1,j),DA1,DB1)


c
            do jz=1,nzmx(j)
               jzx=jx0 + jz
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
               phix(jzx)=A1tot*ai1 + B1tot*bi1
               dphix(jzx)=-eta(j)*phsum
               Dphix_w(jzx)=A1tot*Dai1 + B1tot*Dbi1 + 
     .            DA1tot*ai1 + DB1tot*bi1
               Dphsum_w=A1tot*Daip1 + B1tot*Dbip1 + 
     .            DA1tot*aip1 + DB1tot*bip1 
               Ddphix_w(jzx)=-eta(j)*Dphsum_w - Deta*phsum
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
         Dphiz_wx(jsr,jm)=sg_rho*Dphix_w(jx)
         Ddphiz_wx(jsr,jm)=sg_rho*Ddphix_w(jx)
cc       call Dphz_calc(w,wsq,vg,cp_sr(jsr),xkh,xkhsq,phiq,dphiq,
cc   .      Dphix_w(jx),Ddphix_w(jx),phimagx(jsr,jm),phiphzx(jsr,jm),
cc   .      Dphimagx(jsr,jm),Dphiphzx(jsr,jm),iiev)
cc    print *,jm,jsr,sngl(phimagx(jsr,jm)),
cc   .   sngl(phiphzx(jsr,jm)/pie)
cc    if(jsr .eq. 2) then
cc       if(mod(nmode,2) .eq. 0) then
cc          eig_char(1,nmode)=phimagx(jsr,jm)
cc          eig_char(2,nmode)=Dphimagx(jsr,jm)
cc       else
cc          eig_char(1,nmode)=phiphzx(jsr,jm)
cc          eig_char(2,nmode)=Dphiphzx(jsr,jm)
cc       endif
cc    endif
      enddo
c
      return
      end
ccc
      subroutine Dphz_calc(w,wsq,vg,cp,xkh,xkhsq,phi,dphi,Dphi_w,
     .   Ddphi_w,phimag,phiphz,Dphimag,Dphiphz,iiev)
c
      implicit none
      integer*4 iiev
      real*8 w,wsq,vg,cp,xkh,xkhsq,phi,dphi,Dphi_w,Ddphi_w,Ksq,gamsq,
     .   gamma,Dgamma,gamphi,Dgamphi,Mgamsq,phimag,phiphz,Dphimag,
     .   Dphiphz
c
      Ksq=wsq/(cp*cp)
      if(Ksq .gt. xkhsq) then
         iiev=0
         gamsq=Ksq - xkhsq
         gamma=sqrt(gamsq)
         Dgamma=(Ksq/w - xkh/vg)/gamma
         gamphi=gamma*phi
         Dgamphi=gamma*Dphi_w + Dgamma*phi
         Mgamsq=gamphi*gamphi + dphi*dphi
         phimag=sqrt(Mgamsq/gamsq)
cc       if(gamphi .eq. 0.d0 .and. dphi .eq. 0.d0) print *,
cc   .      'datan2 in Dphz_calc: '
         phiphz=datan2(gamphi,dphi)
         Dphimag=phimag*((dphi*Ddphi_w + gamphi*Dgamphi)/Mgamsq 
     .      - Dgamma/gamma)
         Dphiphz=(dphi*Dgamphi - gamphi*Ddphi_w)/Mgamsq
      else
         iiev=1
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
      Dterm2=(igamma*Ddphi_w(ii) - Digamma*dreal(dphilay(ii)))/gamsq
      DA1=(Dphi_w(ii) + Dterm2)*0.5d0
c
      ii=ii2
      term1=dreal(philay(ii))
      term2=dreal(dphilay(ii))/igamma
      B1=rx_round(term1,-term2,0.5d0,ratb)
cc    if(ratb .lt. 1.d-3) then
cc       term1h=dreal(philay(ii2))
cc       term2h=dreal(dphilay(ii2))/igamma
cc       B1x=rx_round(term1h,-term2h,0.5d0,ratbx)
cc       if(ratbx .gt. ratb) then
cc          B1=B1x
cc          ii=ii2
cc          Dterm2=(igamma*Ddphi_w(ii) - Digamma*dreal(dphilay(ii)))
cc   .         /gamsq
cc       endif
cc    endif
      expB=Aplay(ii) + dreal(zetalay(ii))
      Dterm2=(igamma*Ddphi_w(ii) - Digamma*dreal(dphilay(ii)))/gamsq
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
      complex*16 philay(2),dphilay(2),ailay(2,2),bilay(2,2),zetalay(2)
      real*8 Aplay(2),eta,A1,B1,expA,expB,pie,term1,term2,rx_round,
     .   dphi_eta,rata,ratb,Deta,xi(2),Dxi(2),etasq,Dphi_w(2),
     .   Ddphi_w(2),DA1,DB1,Ai0,Ai0p,Bi0,Bi0p,DAi0,DAi0p,Dbi0,DBi0p,
     .   Ddphi_eta
c
c: ALWAYS FIT AI COEFFICIENT A1 AT INTERFACE CLOSEST TO DUCT IN WHICH MODE
c: IS TRAPPED, AND BI COEFFICIENT B1 AT INTERFACE FARTHEST FROM DUCT.
c
      ii=ii1
      Bi0=dreal(bilay(1,ii))
      Bi0p=dreal(bilay(2,ii))
      dphi_eta=dreal(dphilay(ii))/eta
      term1=Bi0p*dreal(philay(ii))
      term2=Bi0*dphi_eta
      A1=rx_round(term1,term2,pie,rata)
cc    if(rata .lt. 1.d-3) then
cc       dphi_etah=dreal(dphilay(ii2))/eta
cc       Bi0=dreal(bilay(1,ii2))
cc       Bi0p=dreal(bilay(2,ii2))
cc       term1h=Bi0p*dreal(philay(ii2))
cc       term2h=Bi0*dphi_etah
cc       A1x=rx_round(term1h,term2h,pie,ratax)
cc       if(ratax .gt. rata) then
cc          A1=A1x
cc          ii=ii2
cc          dphi_eta=dphi_etah
cc       endif
cc    endif
      expA=Aplay(ii) + dreal(zetalay(ii))
      Ddphi_eta=(eta*Ddphi_w(ii) - Deta*dreal(dphilay(ii)))/etasq
      DBi0=Bi0p*Dxi(ii)
      DBi0p=(xi(ii)*Bi0)*Dxi(ii)
      DA1=pie*(Bi0p*Dphi_w(ii) + DBi0p*
     .   dreal(philay(ii)) + Bi0*Ddphi_eta + DBi0*dphi_eta)
cc    DBi0=bilay(2,ii)*Dxi(ii)
cc    DBi0p=(xi(ii)*bilay(1,ii))*Dxi(ii)
cc    DA1=pie*(dreal(bilay(2,ii))*Dphi_w(ii) + DBi0p*
cc   .   dreal(philay(ii)) + dreal(bilay(1,ii))*Ddphi_eta + 
cc   .   DBi0*dphi_eta)
c
      ii=ii2
      dphi_eta=dreal(dphilay(ii))/eta
      Ai0=dreal(ailay(1,ii))
      Ai0p=dreal(ailay(2,ii))
      term1=Ai0p*dreal(philay(ii))
      term2=Ai0*dphi_eta
      B1=rx_round(term1,term2,-pie,ratb)
cc    if(ratb .lt. 1.d-3) then
cc       dphi_etah=dreal(dphilay(ii1))/eta
cc       Ai0=dreal(ailay(1,ii1))
cc       Ai0p=dreal(ailay(2,ii1))
cc       term1=Ai0p*dreal(philay(ii1))
cc       term2=Ai0*dphi_etah
cc       B1x=rx_round(term1h,term2h,-pie,ratbx)
cc       if(ratbx .gt. ratb) then
cc          B1=B1x
cc          ii=ii1
cc          dphi_eta=dphi_etah
cc       endif
cc    endif
      expB=Aplay(ii) - dreal(zetalay(ii))
c
      Ddphi_eta=(eta*Ddphi_w(ii) - Deta*dreal(dphilay(ii)))/etasq
      DAi0=Ai0p*Dxi(ii)
      DAi0p=(xi(ii)*Ai0)*Dxi(ii)
      DB1=-pie*(Ai0p*Dphi_w(ii) + DAi0p*dreal(philay(ii)) + 
     .   Ai0*Ddphi_eta + DAi0*dphi_eta)
cc    DAi0=ailay(2,ii)*Dxi(ii)
cc    DAi0p=(xi(ii)*ailay(1,ii))*Dxi(ii)
cc    DB1=-pie*(ailay(2,ii)*Dphi_w(ii) + DAi0p*dreal(
cc   .   philay(ii)) + ailay(1,ii)*Ddphi_eta + DAi0*dphi_eta)
c
      return
      end
