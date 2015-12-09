      subroutine phz_calc
c
c: Computes number of zeros in mode functions.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 j,iiev_top,iiev_bot
      real*4 phz_top,phz_bot,phz_lay(NLMAX),phz_phi,
     .   twpiex,piex,phz_bl,phz_psi,phz_lays
      data twpiex/6.283185307/,piex/3.14159265358979/
c
      phz_phi=0.
      phz_psi=0.
      phz_bl=0.
      mode_no(1,nmode)=0
      mode_no(2,nmode)=0
      do j=2,nlay-1
         call phi_phase(philay(1,j),dphilay(1,j),gami(1,1,j),
     .      xkh,xk(1,j),phz_top,iiev_top)
         call phi_phase(philay(2,j),dphilay(2,j),gami(1,2,j),
     .      xkh,xk(2,j),phz_bot,iiev_bot)
         call phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .      isp(j),h(j),gami(1,1,j),zetalay(1,1,j),zetalay(2,1,j),
     .      piex,twpiex,phz_lay(j))
c
         phz_phi=phz_phi + phz_lay(j)
         phz_bl=phz_bot
         if(iisol(j) .eq. 1) then
            call phi_phase(psilay(1,j),dpsilay(1,j),beti(1,1,j),
     .         xkh,xb(1,j),phz_top,iiev_top)
            call phi_phase(psilay(2,j),dpsilay(2,j),beti(1,2,j),
     .         xkh,xb(2,j),phz_bot,iiev_bot)
            call phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .         iss(j),h(j),beti(1,1,j),zetalay(1,2,j),zetalay(2,2,j),
     .         piex,twpiex,phz_lays)
            phz_psi=phz_psi + phz_lays
         endif
      enddo
c
      mode_phz(1,nmode)=phz_phi/piex
      mode_phz(2,nmode)=phz_psi/piex
      mode_phz(3,nmode)=mode_phz(1,nmode) + mode_phz(2,nmode)
      print *,'nmode = ',nmode,(mode_phz(j,nmode),j=1,3)
c
      return
      end
ccc
      subroutine phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .   iso,h,gami,zeta_top,zeta_bot,piex,twpie,phz_lay)
c
c: Computes phase change in a layer given phase at top and bottom and
c: iso or airy parameters.
c
      implicit none
      integer*4 iso,iiev_top,iiev_bot
      real*4 phz_top,phz_bot,phz_mod,phz_dif,dphz,twpie,piex,phz_lay
      real*8 h
      complex*16 gami,zeta_top,zeta_bot
c
c: Force layer phase difference to be zero when evanescent:
      if(iiev_top .eq. 1 .and. iiev_bot .eq. 1) then
         phz_lay=0.
         return
      endif
      phz_dif=phz_bot - phz_top
      if(phz_dif .lt. 0.) phz_dif=phz_dif + twpie
      if(iso .eq. 1) then
         phz_lay=dabs(h*dimag(gami))
      else
         phz_lay=abs(dimag(zeta_top)-dimag(zeta_bot))
      endif
      phz_mod=mod(phz_lay,twpie)
      dphz=phz_dif - phz_mod
      if(abs(dphz) .gt. piex) dphz=dphz - sign(twpie,dphz)
      phz_lay=phz_lay + dphz
c
      return
      end
ccc
      subroutine phi_phase(phi,dphi,gami,xkh,xk,phz,iiev)
c
      implicit none
      integer*4 iiev
      complex*16 gami,xkh,xk,phi,dphi
      real*4 phz,gamma,num_re,num_im,mag_re,mag_im
c
      gamma=dimag(gami)
      num_re=dreal(phi)*gamma
      num_im=dimag(phi)*gamma
      mag_re=num_re*num_re+dreal(dphi)*dreal(dphi)
      mag_im=num_im*num_im+dimag(dphi)*dimag(dphi)
      if(mag_re .gt. mag_im) then
         phz=atan2(num_re,sngl(dreal(dphi)))
      else
         phz=atan2(num_im,sngl(dimag(dphi)))
      endif
      iiev=0
      if(dreal(xkh) .gt. dreal(xk)) iiev=1
c
      return
      end
ccc
