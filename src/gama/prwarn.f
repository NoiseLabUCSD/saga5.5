      subroutine prwarn
c
c: this subroutine prints warning messages before program termination.
      implicit integer*4(i-n)
      include 'common/svp'
      include 'common/pii'
      include 'common/causbad'
      include 'common/paths'
      include 'common/timing'
      include 'common/freqcom'
      include 'common/srlox'
      include 'common/bdcom'
c
      if(kph .gt. 0) then
         write(8,100) kph
100   format(/'WARNING: REFL/TRANS PHASE VARIED BY > PI/2 ',i4,
     .   ' TIMES DURING RANGE INTERPOLATION.'/'RTERP MAY HAVE',
     .   ' BEEN TOO LARGE FOR THIS RUN.')
      endif
      if(ntcut .gt. 0) then
         write(8,200) ntcut,twin
200   format(/'WARNING: ',i4,' RAYS WITHIN CUTOFF WERE CUT DUE TO',
     . ' TRAVEL TIMES OUTSIDE WINDOW'/'OF LENGTH TWIN=NFFT/FS = ',f7.4,
     . '; INCREASE NFFT, DECREASE FS, OR CHOOSE IIIR=2 '/
     .   'OR IIFFT=2 or 3 TO INCLUDE THESE RAYS')
      endif
      if(ncbad .gt. 0) then
         write(8,110) ncbad
110      format(/'WARNING: ',i6,' UNREALISTIC CAUSTICS DETECTED',
     .   ' DURING THE EIGENRAY SEARCH.'/'SUCH CAUSTICS ARE FORMED',
     .   ' BY RAYS THAT TURN AT DEPTHS WHERE THE SVP IS'/          
     .   'IRREGULAR.  LIST OF (SOURCE INCIDENT ANGLE, SV WHERE RAY ',
     .   'TURNS, CAUSTIC RANGE):'/'( THETAS ,   SVTURN ,    CAUSTIC R,',
     .   '   EXTENT-WVL,     EXTENT-m)')
         do 80 jcbad=1,min0(ncbad,20)
            thc=acos(cs/ctbad(jcbad))*piedeg
            write(8,120) thc,ctbad(jcbad),rcbad(jcbad),delbad(jcbad),
     .         delbad(jcbad)*bdlam(nfeig)
120   format('(',f7.3,' , ',f8.3,' , ',f12.2,' , ',f11.2,' , ',f11.2')')
80       continue
      endif
      if(nbdbad .gt. 0) then
         write(8,140) nbdbad
140      format(/'WARNING: ',i6,' RAYS WITH LARGE NEGATIVE BEAM ',
     .      'DISPLACEMENT WERE IGNORED.'/'SUCH BEAM DISPLACEMENT ',
     .      'OCCURS AT FLUID-SOLID AND SOLID-SOLID INTERFACES.'/
     .      'THE PROGRAM AUTHOR DOES NOT KNOW IF THE RAYS ARE VALID.')
      endif
c
      return
      end 
