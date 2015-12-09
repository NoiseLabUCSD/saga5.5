      subroutine eigend(plin,plco)
c
c: this subroutine writes the final line to the eigenray list file.
      implicit integer*4(i-n)
      include 'common/gamaoptions'
      include 'common/timing' 
      include 'common/vhfacs' 
      include 'common/pii' 
      complex*8 plco
      complex plco2
      real*4 plin,c60
      data c60/60./
c
      if(iieig .eq. 1) then
         plin2=9999.99
         if(plin .gt. 1.2e-38) plin2=10.*log10(plin)
         plc=abs(plco)
         if(plc .gt. 1.2e-38) then
            plc=20.*log10(plc)
            plco2=plco
            plph=atan2(aimag(plco2),real(plco2))*piedeg
            iplph=int(plph + sign(.49,plph))
         else
            plc=9999.99
            iplph=0 
         endif
         cptim=etime(cpsec)
         write(59,175) plc,iplph,plin2,dtau,int(cptim/c60),
     .      int(mod(cptim,c60))
175      format('COHERENT PL = ',f7.2,' dB; PHASE = ',i4,' deg',
     .      '; INCOH PL = ',f7.2,' dB','; DTAU = ',f7.3,
     .      ' s; CP MIN = ',i4,':',i2.2)
         write(59,160)
160      format('EOD -----------------------------------------',
     .'-------------------------------------------------------------')
      endif
cxx   if(iibar .eq. 1) call pltend(0.0) 
c
      return
      end 
