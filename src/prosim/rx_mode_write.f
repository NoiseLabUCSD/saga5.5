       subroutine rx_mode_write(jj,nmode,nzsr,phisv,knsv,
     .                          dim1,dim2,dim3,jdf)
c
c: Write modes amplitudes & eigenvalues evaluated.
c
      implicit none
      include 'Parms_com'
c
      integer*4 jm,jk,nmode,nzsr
      integer*4 dim1,dim2,dim3,jj,jdf
      real*4 phisv(dim1,dim2,dim3)
      complex*16 knsv(0:dim2)
c
      do jm=1,nmode
c        write(70+jj,*) dreal(knsv(jm)),dimag(knsv(jm)),jm

       do jk=1,nzsr
        write(71+jj,*) phisv(jk,jm,jdf),jm
       enddo
      enddo

      return
      end
