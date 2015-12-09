      subroutine sigma_fill(phi_,dphi_,psi_,dpsi_,uzx,sig,nmm,
     .   ir,jrun,nrunx,jfcw)
c
c: Fills the array uzx and sig with the magnitude of the vertical and
c: horizontal displacements (uzx) and the normal and shear stress (sig)
c: and outputs them.  See p.67 and p.117.
c
      use parms_com
      use i_o_com
      use gen_com
      implicit none
      integer*4 nmm,ir,jrun,nrunx,jfcw,jmf,jsr,jm
      complex*8 phi_(nzsr,nmode),dphi_(nzsr,nmode),psi_(nzsr,nmode),
     .   dpsi_(nzsr,nmode),ux,uz,sigzz,sigzx
      complex*16 ik_,Acon_,Bcon_
      complex*8 uzx(nzmf,nmm,2),sig(nzmf,nmm,2)
      real*4 re_im_ax(2)
c
c: NOTE: WE DIVIDE PHI,DPHI,PSI,DPSI BY RHO BECAUSE RHO HAS BEEN FACTORED
C: INTO THEM IN MODE_FUN.
c
      do jm=1,nmm
         xmode(jm)=float(jm)
      enddo
      wsq=w*w
      do jmf=1,nzmf
         jsr=mzmf(jmf)
         if(kksh(jsr) .eq. 0) then
            do jm=1,nmode
               ik_=dcmplx(-dimag(kn(jm)),dreal(kn(jm)))
               uzx(jmf,jm,1)=dphi_(jsr,jm)/rho_sr(jsr)
               uzx(jmf,jm,2)=ik_*phi_(jsr,jm)/rho_sr(jsr)
               sig(jmf,jm,1)=-wsq*phi_(jsr,jm)
               sig(jmf,jm,2)=(0.,0.)
            enddo
         else
            do jm=1,nmode
               ik_=dcmplx(-dimag(kn(jm)),dreal(kn(jm)))
               Acon_=wsq*(1.d0 - 2.d0*ksm2_sr(jsr)*kn(jm)*kn(jm))
               Bcon_=ik_*wsq*2.d0*ksm2_sr(jsr)
               uz=(dphi_(jsr,jm) + ik_*psi_(jsr,jm))/rho_sr(jsr)
               ux=(ik_*phi_(jsr,jm) - dpsi_(jsr,jm))/rho_sr(jsr)
               sigzz=-Acon_*phi_(jsr,jm) + Bcon_*dpsi_(jsr,jm)
               sigzx=Bcon_*dphi_(jsr,jm) + Acon_*psi_(jsr,jm)
               uzx(jmf,jm,1)=uz
               uzx(jmf,jm,2)=ux
               sig(jmf,jm,1)=sigzz
               sig(jmf,jm,2)=sigzx
            enddo
         endif
         do jm=nmode+1,nmm
            uzx(jmf,jm,1)=NaN
            uzx(jmf,jm,2)=NaN
            sig(jmf,jm,1)=NaN
            sig(jmf,jm,2)=NaN
         enddo
      enddo
      re_im_ax(1)=1
      re_im_ax(2)=2
cpg      call hdf_write_gen(4+ir,outroot,lout,re_im_ax,2,zmf,nzmf,
cpg     .   xmode,nmm,fcw,nfcw,var_ax,nrunx,uzx(1,1,1),1,2,1,nzmf,
cpg     .   1,nmode,jfcw,jfcw,jrun,jrun,'Re/Im',5,'Depth - m',
cpg     .   9,'Mode',4,'Frequency - Hz',14,'Parameter',9,'uz',
cpgcpg     .   2,20,1,2)
cpg      call hdf_write_gen(4+ir,outroot,lout,re_im_ax,2,zmf,nzmf,
cpg     .   xmode,nmm,fcw,nfcw,var_ax,nrunx,uzx(1,1,2),1,2,1,nzmf,
cpg     .   1,nmode,jfcw,jfcw,jrun,jrun,'Re/Im',5,'Depth - m',
cpg     .   9,'Mode',4,'Frequency - Hz',14,'Parameter',9,'ux',
cpg     .   2,21,1,2)
cpg      call hdf_write_gen(4+ir,outroot,lout,re_im_ax,2,zmf,nzmf,
cpg     .   xmode,nmm,fcw,nfcw,var_ax,nrunx,sig(1,1,1),1,2,1,nzmf,
cpg     .   1,nmode,jfcw,jfcw,jrun,jrun,'Re/Im',5,'Depth - m',
cpg     .   9,'Mode',4,'Frequency - Hz',14,'Parameter',9,'sigzz',
cpg     .   5,22,1,2)
cpg      call hdf_write_gen(4+ir,outroot,lout,re_im_ax,2,zmf,nzmf,
cpg     .   xmode,nmm,fcw,nfcw,var_ax,nrunx,sig(1,1,2),1,2,1,nzmf,
cpg     .   1,nmode,jfcw,jfcw,jrun,jrun,'Re/Im',5,'Depth - m',
cpg     .   9,'Mode',4,'Frequency - Hz',14,'Parameter',9,'sigzx',
cpg     .   5,23,1,2)
      return
      end
