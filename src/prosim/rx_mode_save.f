       subroutine rx_mode_save(numsectr,nmode,nzsr,ncx,eig_charsv,
     . dphizw, ddphizw,phisv,dphisv,knsv,dim1,dim2,dim3,jdf)
c
c: Saves Modes amplitudes & eigenvalues evaluated at fmin.
c
      implicit none
      include 'Parms_com'
      include 'modrng_com'
c
      integer*4 numsectr,jm,jk,nmode,nzsr,ncx
      integer*4 dim1,dim2,dim3,jj,jdf
      real*4 dphizw(dim1,dim2,dim3), ddphizw(dim1,dim2,dim3)
      real*4 phisv(dim1,dim2,dim3),dphisv(dim1,dim2,dim3)
      complex*16 eig_charsv(5,dim2)
      complex*16 knsv(0:dim2)
c
      nmode_rng(numsectr)=nmode
      ncx_rng(numsectr)=ncx
c
      do jm=1,nmode
        kn_rng(jm,numsectr)=knsv(jm)
        eig4_rng(jm,numsectr)=eig_charsv(4,jm)
        eig5_rng(jm,numsectr)=eig_charsv(5,jm)
      enddo

      do jj=1,nmode
       do jk=1,nzsr
        phi_rng(jk,jj,numsectr)=phisv(jk,jj,jdf)
        dphi_rng(jk,jj,numsectr)=dphisv(jk,jj,jdf)
        dphiz_rng(jk,jj,numsectr)=dphizw(jk,jj,jdf)
        ddphiz_rng(jk,jj,numsectr)=ddphizw(jk,jj,jdf)
       enddo
      enddo

      return
      end

ccc

       subroutine rx_mode_rcl(numsectr,nmode,nzsr,ncx,eig_charsv,
     . dphizw, ddphizw,phisv,dphisv,knsv,dim1,dim2,dim3,jdf)
c
c: Recall modes amplitudes & eigenvalues evaluated at fmin.
c
      implicit none
      include 'Parms_com'
      include 'modrng_com'
c
      integer*4 numsectr,jm,jk,nmode,nzsr,ncx
      integer*4 dim1,dim2,dim3,jj,jdf
      real*4 dphizw(dim1,dim2,dim3), ddphizw(dim1,dim2,dim3)
      real*4 phisv(dim1,dim2,dim3),dphisv(dim1,dim2,dim3)
      complex*16 eig_charsv(5,dim2)
      complex*16 knsv(0:dim2)
c
      nmode=nmode_rng(numsectr)
      ncx=ncx_rng(numsectr)
c
      do jm=1,nmode
        knsv(jm)=kn_rng(jm,numsectr)
        eig_charsv(4,jm)=eig4_rng(jm,numsectr)
        eig_charsv(5,jm)=eig5_rng(jm,numsectr)
      enddo

      do jj=1,nmode
       do jk=1,nzsr
        phisv(jk,jj,jdf)=phi_rng(jk,jj,numsectr)
        dphisv(jk,jj,jdf)=dphi_rng(jk,jj,numsectr)
        dphizw(jk,jj,jdf)=dphiz_rng(jk,jj,numsectr)
        ddphizw(jk,jj,jdf)=ddphiz_rng(jk,jj,numsectr)
       enddo
      enddo

      return
      end
