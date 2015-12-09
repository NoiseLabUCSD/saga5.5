      subroutine eig_disp(kn,nf_kn,nmode,nfreq,jf1,jf2,xmode,iiopt,
     .   kw0,r4mat1,r4mat2)
c
      implicit none
      integer*4 nf_kn,jf1,jf2,iiopt,jm,jf,nmode,nfreq
      complex*16 kn(nf_kn,nmode)
      real*8 kw0,twpie,piedeg
      real*4 xmode(nmode),r4mat1(nfreq,nmode),r4mat2(nfreq,nmode)
      data piedeg/57.29577951308232/,twpie/6.28318530717959/
c
      do jm=1,nmode
         xmode(jm)=jm
      enddo
      if(iiopt .eq. 1) then
         do jf=jf1,jf2
            do jm=1,nmode
               r4mat1(jf,jm)=real(kn(jf,jm))
c: Compute attenuation of mode as a function of range in dB/km:
c: [20*log10(exp(-Im(kn)*r))=-Im(kn)*r*(20*log10(e))=-Im(kn)*r*8.6859]
cxx            r4mat2(jf,jm)=8685.9*dimag(kn(jf,jm))
               r4mat2(jf,jm)=dimag(kn(jf,jm))
            enddo
         enddo
      else
         do jf=jf1,jf2
            do jm=1,nmode
               r4mat1(jf,jm)=real(kn(jf,jm))/kw0
c: Compute attenuation of mode as a function of range in dB/km:
c: [20*log10(exp(-Im(kn)*r))=-Im(kn)*r*(20*log10(e))=-Im(kn)*r*8.6859]
cxx            r4mat2(jf,jm)=8685.9*dimag(kn(jf,jm))
               r4mat2(jf,jm)=dimag(kn(jf,jm))
            enddo
         enddo
      endif
c
      if(iiopt .eq. 3) then
         do jf=jf1,jf2
            do jm=1,nmode
               r4mat1(jf,jm)=acos(r4mat1(jf,jm))*piedeg
            enddo
         enddo
      endif
c
      return
      end
