      SUBROUTINE TLRANPG(NP,*,F,MSP,RNG,
     & US,UROLD,URNEW,XT1,XT2,xsou,ISP1,NSECT,pressr,nptot,nr_start,
     & xranges,isect)
 
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
      parameter(Mdep=1000)
      real rdep(Mdep)
      common /receiverdepth/rdep

      common /tiltparm/tilt,dtilt,incoh

C
  300 FORMAT(1H ,F10.2,2X,3(e11.4,2x))
  400 FORMAT(1H ,/,'   TRANSMISSION LOSS VERSUS RANGE BY COHERENT ADDI
     *TION OF MODES',/,'   FREQUENCY      =',F9.1,' Hz',
     *            /,'   SOURCE DEPTH   =',F9.1,' m',/,'   RECEIVER DEPTH
     * =',F9.1,' m',/,'   RANGE(km)      TL(dB)')
C
c      write(*,*)'tlran0:',nptot,nr_start,np,nrd
c      IF(NP.EQ.0)RETURN
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
c        NPMAX=0
c        NPMIN=0
        NP2=1
        NP=1
c        write(*,*)'isect,nsect,rng(k),rd(k),rstart,rend',
c     &             isect,nsect,rng(k),rdep(k),rstart,rend
        if   (((isect.eq.nsect) .and. (rng(K).gt.Rstart))
     &   .or. ((rng(k).le.rend) .and. (rng(K).gt.Rstart))) then
c        write(*,*)' computing', k
c 1900   CONTINUE
c        if (tilt) then
c          rng(1)=rngref+dtilt/(nrd-1)*(k-1)
c        endif
c        write(*,*)'rd(k),rng(k)',rdep(k),MODQTY,NP,NP1,NP2,RNG(k),
c     &  NRANGE,MSP,MODEN,DELR,ISP1,NSECT
        CALL RECEIV(rdep(K),MODQTY,NP,NP1,NP2,*2400,RNG(k),
     &  NRANGE,MSP,UROLD,URNEW,XT1,XT2,MODEN,DELR,ISP1,NSECT)
c        IF(NPMIN.EQ.0)NPMIN=MAX(1,NP1)
c        NPMAX=NP2
        NP2=1
C   COHERENT ADDITION OF MODES STARTS HERE
c       write(*,*)'tlran1:',nptot,k,nr_start,np1,np2,np
c        CALL COHRNT(RD(K),NP1,NP2,RNG(k),
c     &    pressr(nptot*(k-1)+nr_start),US,
        CALL COHRNT(rdep(K),NP1,NP2,RNG(k),
     &    pressr(k),US,
     &   UROLD,URNEW,MODEN,DELTA,MSP,ALFOLD,ALFNEW,ALFAVR,
     &   EKOLD,EKNEW,EKAVR)
c         write(*,*)'pressr(k)',k,pressr(k)
c        NP1=NP2+1
c        IF(NP2.NE.NP)GOTO 1900
2400  CONTINUE

c      IF(flagpu.lt.0)  then
c        WRITE(LUPRT,400)F,SD,rd(K)
c        DO  2600   L=NPMIN,NPMAX
c         RKM=RNG(L)*1.0E-3
c 2600    WRITE(LUPRT,300)RKM,pressr(nptot*(k-1)+nr_start-1+l),
c     1             -20*log10(abs(pressr(nptot*(k-1)+nr_start-1+l)))
c      endif
      else
c        write(*,*)' skipping'

      endif !range K
4000  CONTINUE
 4001 CONTINUE
      nr_start=nr_start+np
c      if (np.ne.(npmax-npmin+1)) then
c           write(*,*)'***** tlran: warning a receiver is in', 
c     1                ' the sediment !'
c          write(*,*)'nr_start,np,npmax,npmin'
c          write(*,*)nr_start,np,npmax,npmin
c       stop
c     endif
c      write(*,*)'exit tlranpg'
      RETURN 1
      END








