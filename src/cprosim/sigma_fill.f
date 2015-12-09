      subroutine sigma_fill(phi,dphi,psi,dpsi,nm,nzsr,mzmf,nzmf,uzx,
     .   sig,rho_sr,ksm2_sr,kksh,kn,zmf,outroot,lout,hdf_suf,lsuf)
c
c: Fills the array uzx and sig with the magnitude of the vertical and
c: horizontal displacements (uzx) and the normal and shear stress (sig)
c: and outputs them.  See p.67 and p.117.
c
      implicit none
      include 'Parms_com'
      include 'gen_com'
      include 'lab_com'
      integer*4 nm,nzsr,nzmf,mzmf(nzmf),kksh(nzsr),jm,jmf,jmf2,jsr,
     .   lout,lsuf
      complex*8 phi(nzsr,nm),dphi(nzsr,nm),psi(nzsr,nm),dpsi(nzsr,nm),
     .   ux,uz,sigzz,sigzx
      complex*16 ksm2_sr(nzsr),kn(0:nm),ik,Acon,Bcon
      real*4 uzx(nm,2*nzmf,2),sig(nm,2*nzmf,2),zmf(nzmf)
      real*8 rho_sr(nzsr)
      character*64 outroot,hdf_suf
c
c: NOTE: WE DIVIDE PHI,DPHI,PSI,DPSI BY RHO BECAUSE RHO HAS BEEN FACTORED
C: INTO THEM IN MODE_FUN.
c
      do jm=1,nm
         xmode(jm)=float(jm)
      enddo
      wsq=w*w
      do jmf=1,nzmf
         jsr=mzmf(jmf)
         jmf2=jmf + nzmf
         if(kksh(jsr) .eq. 0) then
            do jm=1,nm
               ik=dcmplx(-dimag(kn(jm)),dreal(kn(jm)))
               uzx(jm,jmf,1)=real(dphi(jsr,jm)/rho_sr(jsr))
               uzx(jm,jmf,2)=imag(dphi(jsr,jm)/rho_sr(jsr))
               uzx(jm,jmf2,1)=real(ik*phi(jsr,jm)/rho_sr(jsr))
               uzx(jm,jmf2,2)=imag(ik*phi(jsr,jm)/rho_sr(jsr))
               sig(jm,jmf,1)=real(-wsq*phi(jsr,jm))
               sig(jm,jmf,2)=imag(-wsq*phi(jsr,jm))
               sig(jm,jmf2,1)=0.
               sig(jm,jmf2,2)=0.
            enddo
         else
            do jm=1,nm
               ik=dcmplx(-dimag(kn(jm)),dreal(kn(jm)))
               Acon=wsq*(1.d0 - 2.d0*ksm2_sr(jsr)*kn(jm)*kn(jm))
               Bcon=ik*wsq*2.d0*ksm2_sr(jsr)
               uz=(dphi(jsr,jm) + ik*psi(jsr,jm))/rho_sr(jsr)
               uzx(jm,jmf,1)=real(uz)
               uzx(jm,jmf,2)=imag(uz)
               ux=(ik*phi(jsr,jm) - dpsi(jsr,jm))/rho_sr(jsr)
               uzx(jm,jmf2,1)=real(ux)
               uzx(jm,jmf2,2)=imag(ux)
               sigzz=-Acon*phi(jsr,jm) + Bcon*dpsi(jsr,jm)
               sig(jm,jmf,1)=real(sigzz)
               sig(jm,jmf,2)=imag(sigzz)
               sigzx=Bcon*dphi(jsr,jm) + Acon*psi(jsr,jm)
               sig(jm,jmf2,1)=real(sigzx)
               sig(jm,jmf2,2)=imag(sigzx)
            enddo
         endif
      enddo
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'uzre',
     .   lsuf+5,uzx(1,1,1),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,z4,2,
     .   'uzre','m',' ',' ','f6.2','f5.2','f7.1',ncall)
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'uzim',
     .   lsuf+5,uzx(1,1,2),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,z4,2,
     .   'uzim','m',' ',' ','f6.2','f5.2','f7.1',ncall)
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'uxre',
     .   lsuf+5,uzx(1,nzmf+1,1),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,
     .   z4,2,'uxre','m',' ',' ','f6.2','f5.2','f7.1',ncall)
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'uxim',
     .   lsuf+5,uzx(1,nzmf+1,2),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,
     .   z4,2,'uxim','m',' ',' ','f6.2','f5.2','f7.1',ncall)
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'sigzzre',
     .   lsuf+8,sig(1,1,1),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,z4,2,
     .   'sigzzre','m',' ',' ','f6.2','f5.2','f7.1',ncall)
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'sigzzim',
     .   lsuf+8,sig(1,1,2),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,z4,2,
     .   'sigzzim','m',' ',' ','f6.2','f5.2','f7.1',ncall)
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'sigzxre',
     .   lsuf+8,sig(1,nzmf+1,1),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,
     .   z4,2,'sigzxre','m',' ',' ','f6.2','f5.2','f7.1',ncall)
      call out_writex(outroot,lout,hdf_suf(1:lsuf)//SUFX//'sigzxim',
     .   lsuf+8,sig(1,nzmf+1,2),zmf,xmode,nzmf,nm,dlab,mnlab,z4,z4,z4,
     .   z4,2,'sigzxim','m',' ',' ','f6.2','f5.2','f7.1',ncall)
c
      return
      end
