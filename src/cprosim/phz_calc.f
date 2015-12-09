      subroutine phz_calc
c
c: Computes number of zeros in mode functions.
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
      subroutine nzero_lay(phz_lay,piex,phi_t,phi_b,dphi_t,dphi_b,
     .   iiev_t,iiev_b,nz)
c
      implicit none
      integer*4 nz,iiev_t,iiev_b,iwave_t,iwave_b
      real*4 phz_lay,piex
      complex*8 phi_t,phi_b,dphi_t,dphi_b
c
      if(iiev_t .eq. 1 .and. iiev_b .eq. 1) then
         nz=0
         iwave_t=0
         if((real(phi_t) .gt. 0. .and. real(dphi_t) .lt. 0.) .or.
     .      (real(phi_t) .lt. 0. .and. real(dphi_t) .gt. 0.)) iwave_t=1
         iwave_b=0
         if((real(phi_b) .gt. 0. .and. real(dphi_b) .gt. 0.) .or.
     .      (real(phi_b) .lt. 0. .and. real(dphi_b) .lt. 0.)) iwave_b=1
c: If interface wave at top and bottom, keep nz=0 and return:
         if(iwave_t .eq. 1 .and. iwave_b .eq. 1) return
      else
         nz=int(phz_lay/piex)
      endif
      if(mod(nz,2) .eq. 0) then
         if((real(phi_t) .gt. 0. .and. real(phi_b) .le. 0.) .or.
     .      (real(phi_t) .lt. 0. .and. real(phi_b) .ge. 0.)) nz=nz+1
      else
         if((real(phi_t) .gt. 0. .and. real(phi_b) .ge. 0.) .or.
     .      (real(phi_t) .lt. 0. .and. real(phi_b) .le. 0.)) nz=nz+1
      endif
c
      return
      end
ccc
      subroutine phz_calc_old
c
c: Computes total phase of mode function.
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
      integer*4 j,iiev(2,NLMAX),iievs
      real*4 phz_top,phz_bot,phz_lay(NLMAX),phz_mod,
     .   phz_phi,phz_dif,dphz,twpiex,piex,phz_bl,
     .   phz_duct,phz_psi,phz_lays
      data twpiex/6.28318530717959/,piex/3.14159265358979/
c
      phz_phi=0.
      phz_psi=0.
      phz_bl=0.
      do j=2,nlay-1
         call phi_phase(philay(1,j),dphilay(1,j),gami(1,1,j),
     .      xkh,xk(1,j),phz_top,iiev(1,j))
         call phi_phase(philay(2,j),dphilay(2,j),gami(1,2,j),
     .      xkh,xk(2,j),phz_bot,iiev(2,j))
         phz_dif=phz_bot - phz_top
c: Force layer phase difference to be zero when evanescent:
         if(iiev(1,j) .eq. 1 .and. iiev(2,j) .eq. 1) phz_dif=0
         if(phz_dif .lt. 0.) phz_dif=phz_dif + twpiex
         if(isp(j) .eq. 1) then
            phz_lay(j)=dabs(h(j)*dimag(gami(1,1,j)))
         else
            phz_lay(j)=abs(dimag(zetalay(2,1,j))-dimag(zetalay(1,1,j)))
         endif
         phz_mod=mod(phz_lay(j),twpiex)
         dphz=phz_dif - phz_mod
         if(abs(dphz) .gt. piex) dphz=dphz - sign(twpiex,dphz)
         phz_lay(j)=phz_lay(j) + dphz
         phz_phi=phz_phi + phz_lay(j)
      print *,nmode,j,iiev(1,j),iiev(2,j),phz_lay(j)/piex,phz_top/piex,
     .   phz_bot/piex,(phz_top-phz_bl)/piex
         phz_bl=phz_bot
         if(iisol(j) .eq. 1) then
            call phi_phase(psilay(1,j),dpsilay(1,j),beti(1,1,j),
     .         xkh,xb(1,j),phz_top,iievs)
            call phi_phase(psilay(2,j),dpsilay(2,j),beti(1,2,j),
     .         xkh,xb(2,j),phz_bot,iievs)
            phz_dif=phz_bot - phz_top
            if(phz_dif .lt. 0.) phz_dif=phz_dif + twpiex
            if(iss(j) .eq. 1) then
               phz_lays=dabs(h(j)*dimag(beti(1,1,j)))
            else
               phz_lays=abs(dimag(zetalay(2,2,j))-dimag(zetalay(1,2,j)))
            endif
            phz_mod=mod(phz_lay(j),twpiex)
            dphz=phz_dif - phz_mod
            if(abs(dphz) .gt. piex) dphz=dphz - sign(twpiex,dphz)
            phz_lays=phz_lays + dphz
            phz_psi=phz_psi + phz_lays
         endif
      enddo
c
c: Compute phase only in current duct (until mode goes evanescent):
      phz_duct=0.
      do j=nsvmin,2,-1
         phz_duct=phz_duct + phz_lay(j)
         if(iiev(1,j) .eq. 1 .or. iiev(2,j) .eq. 1) goto 10
      enddo
10    continue
      do j=nsvmin+1,nlay-1
         phz_duct=phz_duct + phz_lay(j)
         if(iiev(1,j) .eq. 1 .or. iiev(2,j) .eq. 1) goto 20
      enddo
20    continue
      mode_no(1,nmode)=nint(phz_phi/piex)
      mode_no(2,nmode)=nint(phz_psi/piex)
      mode_no(3,nmode)=mode_no(1,nmode) + mode_no(2,nmode)
      print *,'nmode = ',nmode,(mode_no(j,nmode),j=1,3),
     .   phz_phi/piex,phz_psi/piex,phz_duct/piex
c
      return
      end
ccc
      function deg(z)
c
      implicit none
      complex*8 z
      real*4 deg
c
      deg=atan2(aimag(z),real(z))*180./3.14159265
      return
      end
ccc
      function mag(z)
c
      implicit none
      complex*8 z
      real*4 mag
c
      mag=sqrt(real(z)*real(z) + aimag(z)*aimag(z))
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
      subroutine phi_phase_old(phi,dphi,gami,xkh,xk,phz,iiev)
c
      implicit none
      integer*4 iiev
      complex*16 gami,xkh,xk,phi,dphi
      real*4 phz,gamma,num
c
      gamma=dimag(gami)
      num=dreal(phi)*dimag(gami)
      phz=atan2(num,sngl(dreal(dphi)))
      iiev=0
      if(dreal(xkh) .gt. dreal(xk)) iiev=1
c
      return
      end
ccc
      subroutine phz_calc_next
c
c: Computes number of zeros in mode functions.
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
      integer*4 j,iiev_top,iiev_bot,nz
      real*4 phz_top,phz_bot,phz_lay(NLMAX),phz_phi,
     .   twpiex,piex,phz_bl,phz_psi,phz_lays
      data twpiex/6.28318530717959/,piex/3.14159265358979/
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
         call nzero_lay(phz_lay(j),piex,philay(1,j),philay(2,j),
     .      dphilay(1,j),dphilay(2,j),iiev_top,iiev_bot,nz)
         mode_no(1,nmode)=mode_no(1,nmode) + nz
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
c
            call nzero_lay(phz_lays,piex,psilay(1,j),psilay(2,j),
     .         dpsilay(1,j),dpsilay(2,j),iiev_top,iiev_bot,nz)
            mode_no(2,nmode)=mode_no(2,nmode) + nz
         endif
      enddo
c
      mode_no(3,nmode)=mode_no(1,nmode) + mode_no(2,nmode)
      print *,'nmode = ',nmode,(mode_no(j,nmode),j=1,3)
c
      return
      end
