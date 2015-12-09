      subroutine mode_ortho(phi_,dphi_,psi_,dpsi_)
c
c: Checks for mode orthogonality and normality.  See ORCA II pp. 4-5.
c: Note: In order to get normality, set iimf=1 and sample modes very 
c: finely from just inside the upper halfspace to just inside the lower
c: halfspace.  Orthogonality holds over any subinterval of the waveguide
c: as long as the modes are well sampled over the region.  See ORCA II,p36.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jm,jm2,jzmf,jsr,j,j1,ji
      complex*8 phi_(nzsr,nmode),dphi_(nzsr,nmode),
     .   psi_(nzsr,nmode),dpsi_(nzsr,nmode),Qn,Pn,kn_ks,vfac,dvfac
      complex*16 igral,hsp_lower,hsp_upper,dksq,int_terms,above,below,
     .   min_ik
c
      if(zmf(1) .ge. zdep(1) .or. zmf(nzmf) .le. zdep(nlay-1)) then
         print *,'WARNING: MODE NORMALITY WILL ONLY HOLD IF FIRST ' 
         print *,'   AND LAST MODE FUNCTION DEPTHS ARE IN UPPER AND '
         print *,'   LOWER HALFSPACES, RESPECTIVELY.'
      endif
      wsq=w*w
      do jm=1,nmode
         igral=(0.d0,0.d0)
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
c: NOTE 8/31/94: Conjugation does not work (Frisk wrong)!
cxx         igral=igral+conjg(phi_(jsr,jm))*phi_(jsr,jm)
            if(kksh(jsr) .eq. 0) then
               igral=igral + phi_(jsr,jm)*phi_(jsr,jm)/rho_sr(jsr)
            else
c: Try using v=-psi/(i*k) instead of v=psi:
               min_ik=dcmplx(dimag(kn(jm)),-dreal(kn(jm)))
               vfac=psi_(jsr,jm)/min_ik
               dvfac=dpsi_(jsr,jm)/min_ik
               kn_ks=2.*kn(jm)*kn(jm) - 1./ksm2_sr(jsr)
               Qn=phi_(jsr,jm) + dvfac
               Pn=2.*dphi_(jsr,jm) + kn_ks*vfac
               igral=igral + (Qn*Qn + Pn*vfac)/rho_sr(jsr)
            endif
         enddo
         igral=igral*(zmf(2)-zmf(1))
c: Use actual upper and lower limits from phi_,dphi_:
         jsr=mzmf(nzmf)
         call hsp_terms(phi_(jsr,jm),dphi_(jsr,jm),psi_(jsr,jm),
     .      kn(jm),wsq,cp_sr(jsr),cs_sr(jsr),rho_sr(jsr),kksh(jsr),
     .      hsp_lower)
         jsr=mzmf(1)
         call hsp_terms(phi_(jsr,jm),dphi_(jsr,jm),psi_(jsr,jm),
     .      kn(jm),wsq,cp_sr(jsr),cs_sr(jsr),rho_sr(jsr),kksh(jsr),
     .      hsp_upper)
         igral=igral - (hsp_lower - hsp_upper)
         print *,'jm=',jm,sngl(real(igral)),sngl(dimag(igral)),
     .      sngl(abs(hsp_lower)),sngl(abs(hsp_upper))
      enddo
c
      do jm=1,nmode
         do jm2=jm+1,nmode
            igral=(0.,0.)
            do jzmf=1,nzmf
               jsr=mzmf(jzmf)
c: NOTE 8/31/94: Conjugation does not work!
cxx            igral=igral+conjg(phi_(jsr,jm))*phi_(jsr,jm2)
               igral=igral+phi_(jsr,jm)*phi_(jsr,jm2)/rho_sr(jsr)
            enddo
            dksq=kn(jm)*kn(jm) - kn(jm2)*kn(jm2)
            igral=dksq*igral*(zmf(2)-zmf(1))
c: Account for fluid-solid and solid-solid interfaces:
            int_terms=cmplx(0.d0,0.d0)
            do ji=1,n_int
               j=mzint(ji)
               j1=j+1
               above=(phi_(j,jm2)*dphi_(j,jm)
     .             - phi_(j,jm)*dphi_(j,jm2))/rho_sr(j)
               below=(phi_(j1,jm2)*dphi_(j1,jm)
     .             - phi_(j1,jm)*dphi_(j1,jm2))/rho_sr(j1)
               int_terms=int_terms + (above - below)
            enddo
            jsr=mzmf(nzmf)
            hsp_lower=(phi_(jsr,jm2)*dphi_(jsr,jm) - 
     .         phi_(jsr,jm)*dphi_(jsr,jm2))/rho_sr(jsr)
            jsr=mzmf(1)
            hsp_upper=(phi_(jsr,jm2)*dphi_(jsr,jm) - 
     .         phi_(jsr,jm)*dphi_(jsr,jm2))/rho_sr(jsr)
            igral=igral - int_terms - (hsp_lower - hsp_upper)
         print *,'jm,jm2=',jm,jm2,sngl(abs(igral)),
     .      sngl(abs(hsp_lower)),sngl(abs(hsp_upper)),
     .      sngl(abs(int_terms))
         enddo
      enddo
c      
      return
      end
ccc
      subroutine hsp_terms(phi,dphi,psi,kn,wsq,cp,cs,rho,kksh,
     .   hsp_term)
c
c: Computes halfspace terms of mode normalization.
c
      implicit none
      complex*16 kn,cp,cs,knsq,gami,beti,ikn,hsp_term
      complex*8 phi,dphi,psi
      real*8 wsq,rho
      integer*4 kksh
c
      if(kksh .eq. 0) then
         if(dphi .ne. (0.,0.)) then
            hsp_term=phi*phi*phi/(2.*dphi*rho)
         else
            hsp_term=(0.d0,0.d0)
         endif
      else
         knsq=kn*kn
         gami=cdsqrt(knsq - wsq/(cp*cp))
         beti=cdsqrt(knsq - wsq/(cs*cs))
         ikn=dcmplx(-dimag(kn),dreal(kn))
         hsp_term=(phi*phi/(2.*gami) + 2.*phi*psi/ikn + psi*psi*
     .      (2.*beti*beti + knsq)/(2.*ikn*ikn*beti))/rho
      endif
c
      return
      end
ccc
      subroutine mode_ortho_analytic
c
c: Checks for mode orthogonality and normality using analytic formulas.  
c: See Levinson et al, Eqs. (28) and (A6).
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 j
      complex*16 igral,xx
c
      igral=(0.d0,0.d0)
      do j=2,nlay-1
         call mode_ortho_lay(philay(1,j),dphilay(1,j),Aplay(1,j),
     .      xkhsq,xksq(1,j),geo(1,3,j),eta(j),etasq(j),h(j),isp(j),xx)
         igral=igral + xx
      enddo
      j=nlay
      philay(2,j)=0.d0
      dphilay(2,j)=0.d0
      call mode_ortho_lay(philay(1,j),dphilay(1,j),Aplay(1,j),
     .   xkhsq,xksq(1,j),geo(1,3,j),eta(j),etasq(j),h(j),isp(j),xx)
      igral=igral + xx
c
      print *,'nmode = ',nmode,kduct,dreal(xkh)/kw0,dimag(xkh)*8685.9,
     .   cdabs(igral)
c
      return
      end
ccc
      subroutine mode_ortho_lay(philay,dphilay,Aplay,xkhsq,xksq,rho,
     .   eta,etasq,h,isp,xx)
c
      implicit none
      complex*16 philay(2),dphilay(2),xkhsq,xksq(2),eta,etasq,xx,
     .   phi1sq,dphi1sq,phi2sq,dphi2sq,mgam1sq,mgam2sq,
     .   term1,term2,pdp1,pdp2
      real*8 Aplay(2),rho(2),h,rhoexp1,rhoexp2
      integer*4 isp
c
      phi1sq=philay(1)*philay(1)
      dphi1sq=dphilay(1)*dphilay(1)
      phi2sq=philay(2)*philay(2)
      dphi2sq=dphilay(2)*dphilay(2)
      mgam1sq=xkhsq - xksq(1)
      mgam2sq=xkhsq - xksq(2)
      rhoexp1=rho(1)*dexp(2.d0*Aplay(1))
      rhoexp2=rho(2)*dexp(2.d0*Aplay(2))
      if(isp .eq. 0) then
         term1=rhoexp1*(dphi1sq - mgam1sq*phi1sq)
         term2=rhoexp2*(dphi2sq - mgam2sq*phi2sq)
         xx=(term2 - term1)/(eta*etasq)
      else
         pdp1=philay(1)*dphilay(1)
         pdp2=philay(2)*dphilay(2)
         term1=rhoexp1*pdp1/mgam1sq
         term2=rhoexp2*(phi2sq*h - (dphi2sq*h - pdp2)/mgam2sq)
         xx=0.5d0*(term2 - term1)
      endif
c
      return
      end
