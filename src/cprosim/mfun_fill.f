      subroutine mfun_fill(phi_,mfun,mfunph,suf,lsuf)
c
c: Outputs Re/Im and/or Mag/Phase of mode function phi_.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm,jzmf,jsr,lsuf
      complex*8 phi_(nzsr,nmode)
      real*4 mfun(nmode,nzmf),mfunph(nmode,nzmf),piedeg
      character*64 suf
c
      piedeg=acos(-1.)/180.
      do jm=1,nmode
         xmode(jm)=jm
      enddo
c
      if(iiri .ne. 0) then
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
            do jm=1,nmode
               mfun(jm,jzmf)=real(phi_(jsr,jm))
               mfunph(jm,jzmf)=aimag(phi_(jsr,jm))
            enddo
         enddo
c      
         if(iiri .eq. 1 .or. iiri .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'re',lsuf+2,
     .         mfun,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         mrlab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
         if(iiri .eq. 2 .or. iiri .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'im',lsuf+2,
     .         mfunph,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         milab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
      endif
      if(iimp .ne. 0) then
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
            do jm=1,nmode
               mfun(jm,jzmf)=abs(phi_(jsr,jm))
               mfunph(jm,jzmf)=atan2(aimag(phi_(jsr,jm)),
     .            real(phi_(jsr,jm)))/piedeg
            enddo
         enddo
         if(iimp .eq. 1 .or. iimp .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'mag',lsuf+3,
     .         mfun,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         malab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
         if(iimp .eq. 2 .or. iimp .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'ph',lsuf+2,
     .         mfunph,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         mplab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
      endif
c
      return
      end
