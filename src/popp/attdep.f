      FUNCTION ATTDEP (Z, F)

C Attenuation in water as a function of depth and frequency - 
C   Thorpe version

C AUTHOR:	Dale D Ellis
C DATE WRITTEN: 29-Jan-91
C LAST EDIT:    20-Oct-94

C MODIFICATIONS:
C  August-94: Modifications to attenuation handling, including 
C	        the addition of a depth dependence. J Theriault.

C REMARKS:
C  This is to return the version in VMODES that was changed by 
C    Terry Deveau. Note: dB/nepers is now required. 

C  Z - depth in meters
C  F - frequency in Hz
C  ATTDEP - loss in Nepers/m

C  THORP - dB/kiloyard
C  CONV - Conversion dB/kyd to nepers/m
C  CONV1 - Conversion dB/km to nepers/m
C  URICK - See p.108 of Urick's book (3rd ed.) 
C  ALPWC - attenuation in dB/(km-Hz) as a function of equally spaced 
C	   depth points.

      implicit none
      REAL CONV, RLN10,conv1
      PARAMETER (RLN10 = 2.3025851, CONV = 0.001/0.9144 * RLN10/20.) 
      PARAMETER ( CONV1 = 0.001 * RLN10/20.)
      REAL	ATTDEP, Z, F
      REAL	FKSQ, THORP, URICK,tabinp
      real   attflg,depthinc,alpwc(100)
      common /vol_atten/attflg,depthinc,alpwc
      integer iz

      if (attflg.eq.-1.) then
         FKSQ = (F/1000.)**2
         THORP = 0.1*FKSQ/(1.+FKSQ) + 40.*FKSQ/(4100.+FKSQ) 
         ATTDEP = CONV * THORP
      elseif (attflg.eq.-2.) then
         FKSQ = (F/1000.)**2
         THORP = 0.1*FKSQ/(1.+FKSQ) + 40.*FKSQ/(4100.+FKSQ) 
         URICK = THORP + 0.003 + 2.75E-4*FKSQ
         ATTDEP = CONV * URICK
      elseif (attflg.eq.1) then
         attdep =conv1*alpwc(1)*f
      else
         iz = int(z/depthinc)+1
         if (iz.ge.200) then
            iz = 199
            print*, 'WARNING: from ATTDEP'
            print*, 'Depth below those specified for attenuations' 
         endif
         tabinp = (alpwc(iz)-alpwc(iz+1))/depthinc* 
     *             (z-depthinc*(iz+1)) + alpwc(iz+1)
         attdep =conv1*tabinp*f
      endif
C
      RETURN
      END
