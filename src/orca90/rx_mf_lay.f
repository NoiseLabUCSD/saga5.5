      subroutine rx_mf_lay(jlay_ref,ii_ref,vg,nzero)
c
c: Finds mode functions for the current eigenvalue.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jlay_ref,ii_ref,jlay,j,jp1,nm1,nzero,j2
      real*8 vg,phi_s,dphi_s,phi_s_exp,gamsq,rzero,
     .   phztop,phzbot,gamphi,phz_est,phi_re,dphi_re
c
c: Compute mode functions at layer interfaces above reference depth:
      phi_s=philay(1,2)
      dphi_s=dphilay(1,2)
      phi_s_exp=Aplay(1,2)
      j2=jlay_ref-1
      if(ii_ref .eq. 1) j2=j2-1
      do j=2,j2
         jp1=j+1
c: Compute phi,dphi at bottom of layer j using inverse of
c: cumulative downward propagator matrix (see p. 44, ORCA II):
         call rx_philay(philay(2,j),dphilay(2,j),Aplay(2,j),
     .      philay(1,jp1),dphilay(1,jp1),Aplay(1,jp1),phi_s,dphi_s,
     .      phi_s_exp,h11(1,2,j),h12(1,2,j),h21(1,2,j),h22(1,2,j),
     .      h_exp(1,j),rho_prod(2,j),rhofac(j))
      enddo
c
c: Compute mode functions at layer interfaces below reference depth:
      phi_s=dreal(philay(1,nlay))
      dphi_s=dreal(dphilay(1,nlay))
      phi_s_exp=Aplay(1,nlay)
c: Compute phi,dphi just above basement interface:
      nm1=nlay-1
      philay(2,nm1)=dcmplx(phi_s/rhofac(nm1),rzero)
      dphilay(2,nm1)=dcmplx(dphi_s,rzero)
      Aplay(2,nm1)=phi_s_exp
c
      j2=jlay_ref+1
      if(ii_ref .eq. 1) j2=j2-1
      do j=nlay-2,j2,-1
         jp1=j+1
c: Compute phi,dphi at bottom of layer j using inverse of
c: cumulative upward propagator matrix (see p. 44, ORCA II):
         call rx_philay(philay(2,j),dphilay(2,j),Aplay(2,j),
     .      philay(1,jp1),dphilay(1,jp1),Aplay(1,jp1),phi_s,dphi_s,
     .      phi_s_exp,g11(1,2,j),g12(1,2,j),g21(1,2,j),g22(1,2,j),
     .      g_exp(1,j),rho_prod(1,j),rhofac(j))
      enddo
c
c: Compute number of zeros in mode function:
      if(dreal(dphilay(1,2)) .gt. 0.d0) then
         phztop=2.d0
      else
         phztop=1.d0
      endif
c: Count zero at ocean surface
      nzero=1
      do jlay=2,nlay-1
         gamsq=xksq(2,jlay) - xkhsq
         phi_re=dreal(philay(2,jlay))
         dphi_re=dreal(dphilay(2,jlay))
         if(gamsq .gt. 0.d0) then
            if(phi_re .eq. 0.d0) then
               phzbot=0.d0
            else
               gamphi=sqrt(gamsq)*phi_re
               phzbot=datan2(gamphi,dphi_re)/pie
            endif
         elseif(phi_re .gt. 0.d0) then
            if(dphi_re .gt. 0.d0) then
               phzbot=0.25d0
            else
               phzbot=0.75d0
            endif
         else
            if(dphi_re .gt. 0.d0) then
               phzbot=1.75d0
            else
               phzbot=1.25d0
            endif
         endif
         phz_est=phztop + dimag(zetalay(1,1,jlay))/pie
         phzbot=phzbot + 2.d0*nint((phz_est - phzbot)/2.d0)
         if(phzbot - int(phzbot) .eq. 0.d0) phzbot=phzbot - 1.d-8
         nzero=nzero + int(phzbot) - int(phztop)
         phztop=phzbot
      enddo
cc    if(nmode .ne. nzero) print *,'nmode,nzero = ',nmode,nzero
      mode_phz(2,nmode)=nzero
c
      return
      end
ccc
      subroutine rx_mf_lay_Dw(jlay_ref,ii_ref,vg,nzero)
c
c: Finds mode functions for the current eigenvalue.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jlay_ref,ii_ref,jlay,j,jp1,nm1,nzero,j2
      real*8 vg,phi_s,dphi_s,phi_s_exp,Dphi_sw,Ddphi_sw,gamsq,rzero,
     .   phztop,phzbot,gamphi,phz_est,phi_re,dphi_re
c
c: Compute mode functions at layer interfaces above reference depth:
      phi_s=philay(1,2)
      dphi_s=dphilay(1,2)
      phi_s_exp=Aplay(1,2)
      Dphi_sw=Dphi_w(1,2)
      Ddphi_sw=Ddphi_w(1,2)
      j2=jlay_ref-1
      if(ii_ref .eq. 1) j2=j2-1
      do j=2,j2
         jp1=j+1
c: Compute phi,dphi at bottom of layer j using inverse of
c: cumulative downward propagator matrix (see p. 44, ORCA II):
         call rx_philay_Dw(philay(2,j),dphilay(2,j),Aplay(2,j),
     .      philay(1,jp1),dphilay(1,jp1),Aplay(1,jp1),phi_s,dphi_s,
     .      phi_s_exp,h11(1,2,j),h12(1,2,j),h21(1,2,j),h22(1,2,j),
     .      h_exp(1,j),rho_prod(2,j),rhofac(j),Dphi_w(2,j),
     .      Ddphi_w(2,j),Dphi_w(1,jp1),Ddphi_w(1,jp1),Dphi_sw,
     .      Ddphi_sw,vg)
      enddo
c
c: Compute mode functions at layer interfaces below reference depth:
      phi_s=dreal(philay(1,nlay))
      dphi_s=dreal(dphilay(1,nlay))
      phi_s_exp=Aplay(1,nlay)
      Dphi_sw=Dphi_w(1,nlay)
      Ddphi_sw=Ddphi_w(1,nlay)
c: Compute phi,dphi just above basement interface:
      nm1=nlay-1
      philay(2,nm1)=dcmplx(phi_s/rhofac(nm1),rzero)
      dphilay(2,nm1)=dcmplx(dphi_s,rzero)
      Aplay(2,nm1)=phi_s_exp
      Dphi_w(2,nm1)=Dphi_sw/rhofac(nm1)
      Ddphi_w(2,nm1)=Ddphi_sw
cc    eig_char(1,nmode)=philay(2,nm1)*dexp(Aplay(2,nm1))
cc    eig_char(2,nmode)=Dphi_w(2,nm1)*dexp(Aplay(2,nm1))
c
      j2=jlay_ref+1
      if(ii_ref .eq. 1) j2=j2-1
      do j=nlay-2,j2,-1
         jp1=j+1
c: Compute phi,dphi at bottom of layer j using inverse of
c: cumulative upward propagator matrix (see p. 44, ORCA II):
         call rx_philay_Dw(philay(2,j),dphilay(2,j),Aplay(2,j),
     .      philay(1,jp1),dphilay(1,jp1),Aplay(1,jp1),phi_s,dphi_s,
     .      phi_s_exp,g11(1,2,j),g12(1,2,j),g21(1,2,j),g22(1,2,j),
     .      g_exp(1,j),rho_prod(1,j),rhofac(j),Dphi_w(2,j),
     .      Ddphi_w(2,j),Dphi_w(1,jp1),Ddphi_w(1,jp1),Dphi_sw,
     .      Ddphi_sw,vg)
      enddo
cc    eig_char(1,nmode)=dphilay(2,nm1-2)*dexp(Aplay(2,nm1-2))
cc    eig_char(2,nmode)=Ddphi_w(2,nm1-2)*dexp(Aplay(2,nm1-2))
c
c: Compute number of zeros in mode function:
      if(dreal(dphilay(1,2)) .gt. 0.d0) then
         phztop=2.d0
      else
         phztop=1.d0
      endif
c: Count zero at ocean surface
      nzero=1
      do jlay=2,nlay-1
         gamsq=xksq(2,jlay) - xkhsq
         phi_re=dreal(philay(2,jlay))
         dphi_re=dreal(dphilay(2,jlay))
         if(gamsq .gt. 0.d0) then
            if(phi_re .eq. 0.d0) then
               phzbot=0.d0
            else
               gamphi=sqrt(gamsq)*phi_re
               phzbot=datan2(gamphi,dphi_re)/pie
            endif
         elseif(phi_re .gt. 0.d0) then
            if(dphi_re .gt. 0.d0) then
               phzbot=0.25d0
            else
               phzbot=0.75d0
            endif
         else
            if(dphi_re .gt. 0.d0) then
               phzbot=1.75d0
            else
               phzbot=1.25d0
            endif
         endif
         phz_est=phztop + dimag(zetalay(1,1,jlay))/pie
         phzbot=phzbot + 2.d0*nint((phz_est - phzbot)/2.d0)
         if(phzbot - int(phzbot) .eq. 0.d0) phzbot=phzbot - 1.d-8
         nzero=nzero + int(phzbot) - int(phztop)
         phztop=phzbot
      enddo
cc    if(nmode .ne. nzero) print *,'nmode,nzero = ',nmode,nzero
      mode_phz(2,nmode)=nzero
c
      return
      end
ccc
      subroutine rx_philay(phi1,dphi1,Ap1,phi2,dphi2,Ap2,phi_s,
     .   dphi_s,phi_s_exp,g11,g12,g21,g22,g_exp,rho_prod,rhofac)
c
      implicit none
      real*8 phi1,dphi1,Ap1,phi2,dphi2,Ap2,phi_s,dphi_s,
     .   phi_s_exp,g11,g12,g21,g22,g_exp,rho_prod,rhofac
c
c: Compute phi,dphi at bottom of layer j using inverse of
c: cumulative downward propagator matrix (see p. 44, ORCA II):
      phi1=(g22*phi_s - g12*dphi_s)/rho_prod
      dphi1=(g11*dphi_s - g21*phi_s)/rho_prod
      Ap1=phi_s_exp + g_exp
c: Go downward across interface to top of layer j+1:
      phi2=rhofac*phi1
      dphi2=dphi1
      Ap2=Ap1
c
      return
      end
ccc
      subroutine rx_philay_Dw(phi1,dphi1,Ap1,phi2,dphi2,Ap2,phi_s,
     .   dphi_s,phi_s_exp,g11,g12,g21,g22,g_exp,rho_prod,rhofac,
     .   Dphi_w1,Ddphi_w1,Dphi_w2,Ddphi_w2,Dphi_sw,Ddphi_sw,vg)
c
      implicit none
      real*8 phi1,dphi1,Ap1,phi2,dphi2,Ap2,phi_s,dphi_s,
     .   phi_s_exp,g11,g12,g21,g22,g_exp,rho_prod,rhofac,
     .   Dphi_w1,Ddphi_w1,Dphi_w2,Ddphi_w2,Dphi_sw,Ddphi_sw,Df_Dw,vg
c
c: Compute phi,dphi at bottom of layer j using inverse of
c: cumulative downward propagator matrix (see p. 44, ORCA II):
      phi1=(g22*phi_s - g12*dphi_s)/rho_prod
      dphi1=(g11*dphi_s - g21*phi_s)/rho_prod
      Ap1=phi_s_exp + g_exp
      Dphi_w1=(g22*Dphi_sw - g12*Ddphi_sw + Df_Dw(g22,vg)*phi_s - 
     .   Df_Dw(g12,vg)*dphi_s)/rho_prod
      Ddphi_w1=(g11*Ddphi_sw - g21*Dphi_sw + Df_Dw(g11,vg)*dphi_s - 
     .   Df_Dw(g21,vg)*phi_s)/rho_prod
c: Go downward across interface to top of layer j+1:
      phi2=rhofac*phi1
      dphi2=dphi1
      Ap2=Ap1
      Dphi_w2=rhofac*Dphi_w1
      Ddphi_w2=Ddphi_w1
c
      return
      end
ccc
      function Df_Dw(f,vg)
c: Computes the total derivative of the function f w.r.t. w, where the
c: f(2) is the k-derivative and f(3) is the w derivative, and 1/vg is
c: dk/dw.
      implicit none
      real*8 f(3),vg,Df_Dw
c
      Df_Dw=f(3) + f(2)/vg
c
      return
      end
