      subroutine mode_traj(k,r1r2,ii)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 ii
      complex*16 k,r1r2(3,4)
c
cxx   write(21,100) dreal(k)/kw0,dimag(k)/kw0,ii,iish(1,1),iish(1,2),
      write(21,100) dreal(k)/kw0,8685.9*dimag(k),ii,iish(1,1),
     .   iish(1,2),iish(2,1),iish(2,2),dreal(r1r2(1,4)),
     .   dimag(r1r2(1,4)),kduct
100   format(e14.8,1x,e14.8,1x,i2,4(1x,i2),1x,e9.3,1x,f7.2,1x,i1)
c
      return
      end
ccc
      subroutine mode_traj_bp(rr0)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 ii,jj
      complex*16 rr0(3,4)
c
      do ii=1,2
         do jj=1,2
            if(real(xkbp(ii,jj)) .ge. kremin .and. 
     .         real(xkbp(ii,jj)) .le. w/cphlo .and.
     .         xkbp(ii,jj) .ne. (0.d0,0.d0)) then
               call mode_traj(xkbp(ii,jj),rr0,-5)
            endif
         enddo
      enddo
c
      return
      end
