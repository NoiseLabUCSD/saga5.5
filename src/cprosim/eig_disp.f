      subroutine eig_disp(xmode,nmode,nfreq,kn,iiopt,kw0,
     .   eig_re,eig_im)
c
      implicit none
      integer*4 iiopt,jm,jf,nmode,nfreq
      complex*16 kn(nfreq,nmode)
      real*8 kw0,twpie,piedeg
      real*4 xmode(nmode),eig_re(nmode,nfreq),eig_im(nmode,nfreq)
      data piedeg/57.29577951308232/,twpie/6.28318530717959/
c
      do jm=1,nmode
         xmode(jm)=jm
      enddo
      if(iiopt .eq. 1) then
         do jf=1,nfreq
            do jm=1,nmode
               eig_re(jm,jf)=real(kn(jf,jm))
c: Compute attenuation of mode as a function of range in dB/km:
c: [20*log10(exp(-Im(kn)*r))=-Im(kn)*r*(20*log10(e))=-Im(kn)*r*8.6859]
cxx            eig_im(jm,jf)=8685.9*dimag(kn(jf,jm))
               eig_im(jm,jf)=dimag(kn(jf,jm))
            enddo
         enddo
      else
         do jf=1,nfreq
            do jm=1,nmode
               eig_re(jm,jf)=real(kn(jf,jm))/kw0
c: Compute attenuation of mode as a function of range in dB/km:
c: [20*log10(exp(-Im(kn)*r))=-Im(kn)*r*(20*log10(e))=-Im(kn)*r*8.6859]
cxx            eig_im(jm,jf)=8685.9*dimag(kn(jf,jm))
               eig_im(jm,jf)=dimag(kn(jf,jm))
            enddo
         enddo
      endif
c
      if(iiopt .eq. 3) then
         do jf=1,nfreq
            do jm=1,nmode
               eig_re(jm,jf)=acos(eig_re(jm,jf))*piedeg
            enddo
         enddo
      endif
c
      return
      end
