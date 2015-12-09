      SUBROUTINE TLRAN(NP,*,F,MSP,RNG,
     & US,UROLD,URNEW,XT1,XT2,xsou,ISP1,NSECT,pressr,nptot,nr_start,
     & xranges)
 
      real xranges(*)
      REAL*8 EKFST,EKOLD,EKNEW,EKAVR
      PARAMETER (ICF=6,MODEN=200,NRANGE=2001,KSRD=500)
c
C   FLAG IS THE MATRIX THAT CONTAINS ALL FLAG INFORMATION
      COMMON/H/FLAG(ICF)
C
C   RD  CONTAINS THE RECEIVER   DEPTH 
      COMMON/M/RD(KSRD),SD,nrd
      COMMON /I/ SECD(3)
      DIMENSION XT1(moden,MSP),Xt2(moden,MSP),xsou(moden,msp)
      DIMENSION RNG(NRANGE)
      DIMENSION US(MODEN),UROLD(2,MODEN),URNEW(2,MODEN)
      COMMON/B/ALFNEW(MODEN)
      COMMON/C/EKNEW(MODEN)
      COMMON/N/MINMOD,MAXMOD,HORIG,NFREQ
      COMMON/D/EKFST(MODEN),EKOLD(MODEN),EKAVR(MODEN),
     *ALFOLD(MODEN),ALFAVR(MODEN)
      COMMON/NRNGD/HNEW,HOLD,RSTART,REND,SLNGTH,IND1
      COMMON/LUNIT/MSOURC,MODOLD,MODNEW,LUPRT
      COMMON/FLAGPULSE/FLAGPU
      complex pressr(*)
      real ranref
      logical tilt,incoh
      real dtilt
      common /tiltparm/tilt,dtilt,incoh
C
  300 FORMAT(1H ,F10.2,2X,3(e11.4,2x))
  400 FORMAT(1H ,/,'   TRANSMISSION LOSS VERSUS RANGE BY COHERENT ADDI
     *TION OF MODES',/,'   FREQUENCY      =',F9.1,' Hz',
     *            /,'   SOURCE DEPTH   =',F9.1,' m',/,'   RECEIVER DEPTH
     * =',F9.1,' m',/,'   RANGE(km)      TL(dB)')
C
c      write(*,*)'tlran0:',nptot,k,nr_start,np1,np2,np
      IF(NP.EQ.0)RETURN
C
      DELTA=REND-RSTART
C

      MODQTY=MAXMOD-MINMOD+1
      DELR=SECD(3)
      rngref=rng(1)
       IF(MAXMOD.GT.0) CALL SOURCE(SD,MODQTY,
     &  *4001,Xsou,MSP,US,moden)

      DO 4000   K=1,nrd
        NP1=1
        NPMAX=0
        NPMIN=0
 1900   CONTINUE
        if (tilt) then
          rng(1)=rngref+dtilt/(nrd-1)*(k-1)
        endif
c        write(*,*)'rng(1),tilt',rng(1),tilt
        CALL RECEIV(RD(K),MODQTY,NP,NP1,NP2,*2400,RNG,
     &  NRANGE,MSP,UROLD,URNEW,XT1,XT2,MODEN,DELR,ISP1,NSECT)
        IF(NPMIN.EQ.0)NPMIN=MAX(1,NP1)
        NPMAX=NP2
C   COHERENT ADDITION OF MODES STARTS HERE
c        write(*,*)'tlran1:',nptot,k,nr_start,np1,np2,np
        CALL COHRNT(RD(K),NP1,NP2,RNG,
     &    pressr(nptot*(k-1)+nr_start),US,
     &   UROLD,URNEW,MODEN,DELTA,MSP,ALFOLD,ALFNEW,ALFAVR,
     &   EKOLD,EKNEW,EKAVR)
        NP1=NP2+1
        IF(NP2.NE.NP)GOTO 1900
2400  CONTINUE

c      IF(flagpu.lt.0)  then
c        WRITE(LUPRT,400)F,SD,rd(K)
c        DO  2600   L=NPMIN,NPMAX
c         RKM=RNG(L)*1.0E-3
c 2600    WRITE(LUPRT,300)RKM,pressr(nptot*(k-1)+nr_start-1+l),
c     1             -20*log10(abs(pressr(nptot*(k-1)+nr_start-1+l)))
c      endif

4000  CONTINUE
 4001 CONTINUE
      nr_start=nr_start+np
      if (np.ne.(npmax-npmin+1)) then
           write(*,*)'***** tlran: warning a receiver is in', 
     1                ' the sediment !'
       stop
      endif
      RETURN 1
      END








