C ***
C *** SUBROUTINE MODE CALCULATION **************************************
C **********************************************************************
      SUBROUTINE MODE(ISECT,IFREQ,F,MODQTY,*,XTS,MSP,MODPLT,EPSINP,
     & JUMP,JF2,MODAVR,ALFA,EK,MODEN,IPARAM,NMES,DELTAF,SS,FLAGP,
     1      cc0,cc1)
C
C   NDEP IS ALSO DEFINED IN SUB DRIVER AND INPUT

      PARAMETER (NDEP=201)
      PARAMETER (NMODES=1420,NBEG=5,NP555=44500+5*NBEG*NMODES+2,
     & LTOT=NP555-2)

      CHARACTER*3 MODPLT
      DOUBLE PRECISION  CC0, CC1
      REAL MODAVR

      DIMENSION ALFA(MODEN), EK(MODEN), MODAVR(MODEN)
      DIMENSION XTS(moden,MSP), SS(MSP)

      DOUBLE PRECISION DELTAF, EK, TWOPI, F
      DOUBLE PRECISION ADA(NP555), SPEED(NP555), EIGF(NP555),
     & ISO(NMODES), MY(8,NMODES), slow(np555)

      COMMON /A/ C0(NDEP),C1(NDEP),Z0(NDEP),Z1(NDEP),ND0,ND1
      COMMON /AB/ BETA(3),SCATT(2),C2S
      COMMON /G/ R0,R1,R2,C2,H0,H1,FACT0,FACT1,C11
      COMMON /LUNIT/ MSOURC,MODOLD,MODNEW,LUPRT
      COMMON /N/ MINMOD,MAXMOD,HORIG,NFREQ
      COMMON /NRNGD/ HNEW,HOLD,RSTART,REND,SLNGTH,IND1
      COMMON/FLAGPULSE/FLAGPU
      logical tilt
      real dtilt
      common /tiltparm/tilt,dtilt
C
 220  FORMAT(1H ,'  SEGMENT NO.   ',I4,/,'**************',/,
     $ '  SEGMENT LENGTH',2X,'=',F7.1,' KM',/)
 240  FORMAT(1H ,' SAMPLE POINTS:   MSP ',12X,'=',I6,/)
 250  FORMAT(1H ,' DEPTHS (M):',/,'  WATER',11X,'=',F8.2,
     $'     SEDIMENT',8X,'=',F8.2,/)
 260  FORMAT(1H ,' DENSITIES (G/CM**3):',/,'  WATER',11X,'=',F8.2,
     $/,'  SEDIMENT',8X,'=',F8.2,/,'  SUBBOTTOM',7X,'=',F8.2,/)
 270  FORMAT(3X,F8.2,4X,F8.2)
 280  FORMAT(1H ,/,5X,'SOUND SPEED PROFILE',/,11X,' WATER     ',/,
     * '   DEPTH (M)  SPEED (M/S)  ')
 290  FORMAT(1H ,/,'          SEDIMENT ',
     $/,'    DEPTH (M)  SPEED (M/S)')
 500  FORMAT(' MODE CUTOFF AT SOURCE FREQUENCY F=',F10.4,' HZ',/,
     & ' EXECUTION IS TERMINATED FOR THIS SOURCE FREQUENCY')
 600  FORMAT(1H ,' BOTTOM SOUND SPEED   =',F8.2,' M/S',/
     *,' SHEAR SPEED',10X,'=',F8.2,' M/S',/)
 700  FORMAT(1H ,' ATTENUATION COEFFICIENTS (DB/WL): ',
     $        /,'  SEDIMENT',3X,'=',F8.2,'  SUBBOTTOM',3X,'=',F8.2,
     $ '  SHEAR  ',3X,'=',F8.2,/,'  RMS ROUGHNESSES (M):',5X,/,
     $'  SEA SURFACE',5X,'=',F8.2,'  SEA FLOOR',7X,'=',F8.2,/)
C
      MAXMOD=MIN(MAXMOD,NMODES-2)
C
c      DO 1000   I=1,NP555
c        ADA(I)=0.0D0
c        SPEED(I)=0.0D0
c        EIGF(I)=0.0D0
c 1000 CONTINUE
C
C   INPUT PARAMETERS ARE READ FROM FILE 10
C
C   R0 IS DENSITY OF  WATER LAYER.
      R0=1.
      HNEW=H0
 1150 CONTINUE
      NNP=MSP
      IF((flagpu.lT.0) .and.(ifreq.eq.1.))   then
        write(luprt,*)' FREQUENCY NO',ifreq
        SLKM=SLNGTH*1.0E-3
        WRITE(LUPRT,220)ISECT,SLKM
        WRITE(LUPRT,240) NNP
        WRITE(LUPRT,250) H0,H1
        WRITE(LUPRT,260) R0,R1,R2
        WRITE(LUPRT,700) BETA,SCATT
        WRITE(LUPRT,280)
        WRITE(LUPRT,270) (Z0(I)*H0,C0(I),I=1,ND0)
        IF(H1.LE.0.0)   GO TO 1350
          WRITE(LUPRT,290)
          WRITE(LUPRT,270) (Z1(I)*H0,C1(I),I=1,ND1)
 1350   CONTINUE
        WRITE(LUPRT,600) C2,C2S
        if (tilt) write(luprt,*)' Tilt (m)=',dtilt 
      endif
      C11=C1(1)
C
C
      CALL PORTER(F,MINMOD,MAXMOD,NMES,JF2,DELTAF,
     &  MODPLT,EK,XTS,CC0,CC1,ALFA,MODAVR,MODEN,ADA,SPEED,
     &  EIGF,ISO,LTOT,NP555,SS,MSP,NBEG,NMODES,MY,slow)
      MODQTY=MAXMOD-MINMOD+1
      IF(MODQTY.LE.0)   THEN
        WRITE(LUPRT,500) F
        RETURN 1
      END IF
c      REWIND(MODNEW)
      END
C
C
C***********************************************************
C
C
      SUBROUTINE PORTER(F,MMIN,MMAX,NMES1,JF2,
     & DELTAF,MODPLT,EK,XTS,CC0,CC1,ALFA,MODAVR,MODEN,ADA,
     & SPEED,EIGF,ISO,LTOT,NP555,SS,MSP,NBEG,NMODES,MY,slow)
C__________________________________________________________
C                                                          |
C     This routine is acting as main for                   |
C     the routines: ISOINT, STURM, BRENT, NEWTON, CHARAC,  |
C     RICH, SPEED and EIGVEC                               |
C     The subject is to find an approximate solution for   |
C     the eigenvalues of a continuous diff. equation by    |
C     finite difference and extrapolation                  |
C     and to calculate the eigenfunctions and losses.      |
C__________________________________________________________|
C

      PARAMETER (NDEP=201,MAXMSH=8)

      INTEGER N(8),NS(8),MREF(10)

      CHARACTER*3 MODPLT

      REAL NRAT,ALFA(MODEN),MODAVR(MODEN)
      REAL XTS(moden,MSP), SS(MSP)

      DOUBLE PRECISION MY(8,NMODES), ISO(NMODES), H(8), ADA(NP555),
     &                 HS(8), DH0SQ(9)
      DOUBLE PRECISION DELTAF, RKOP, TEMP, BASE, CC0, OMEGA, MYBR
      DOUBLE PRECISION  CC1
      DOUBLE PRECISION DSN
      DOUBLE PRECISION SPEED(NP555),EIGF(NP555),slow(np555)
      DOUBLE PRECISION F,PI,W,K,A,B,BOT,HS1,HRAT,STIFF,Z,DRAT
      DOUBLE PRECISION CON1,CON2,CON3,CON4,CON5,EK(MODEN)
      real*8  twopi
      common /twopie/twopi

      COMMON /A/ C0(NDEP),C1(NDEP),Z0(NDEP),Z1(NDEP),ND0,ND1
      COMMON /AB/ BETA(3),SCATT(2),C2S
      COMMON /ATTEN/ ALF0,ALF1,ALF2,ALF2S,ALFOS,ALFOB
      COMMON /CONST/ CON1,CON2,CON3,CON4,CON5
      COMMON /EXPMAX/ TRESH
      COMMON /FLAGPULSE/ FLAGPU
      COMMON /G/ R0,R1,R2,C2,H0,DS,FACT0,FACT1,C11
      COMMON /GEN/ PI,K,ROB,ROS,BOT,STIFF,CB
      COMMON /LUNIT/ MSOURC,MODOLD,MODNEW,LUPRT
      COMMON /XREFL/ CINTFC,RINTFC
      DATA MREF/50,64,80,100,128,160,200,256,320,400/
      DATA DH0SQ(9)/0.0D0/
C
C     FORMATS
C
  320 FORMAT(1X,5(F10.7,2X))
 360  FORMAT(1H ,//,14X,'WAVE NUMBER',6X,'ALPHA',7X,'A0 ',7X,'A1 ',7X
     *,'A2 ',7X,'A2S',7X,'A0S',7X,'A0B',7X,' RB',2X,' SS',2X,' SB')
  420 FORMAT(I6,3X,' MININUM ORDER MODE',/,
     & I6,3X,' NUMBER OF MESH POINTS IN THE WATER LAYER ',/,
     & I6,3X,' NUMBER OF MESH POINTS IN THE SEDIMENT LAYER ')
c      DO 1234 I12=1,NMODES
c      DO 1234 I13=1,8
c 1234 MY(I13,I12)=xnul
C
C
C     DATA INTERFACE
C
      NMES=NMES1
      ROS=R1/R0
      ROB=R2/R0
      CB=C2/CC0
      PI=twopi*0.5
      K=twopi*F*H0/CC0
      OMEGA=twopi*F
      DSN=DBLE(DS)/DBLE(H0)
      IF(DSN.LT.1.0E-6) ROS=1.0
C
C     INITIALIZATION
C
C     DEFINITION OF MESHES
C
      RKOP=2.0*(F+(2-JF2)*DELTAF)*H0/CC0
      RKOP=MAX(1.0D0,RKOP)
      TEMP=NBEG*RKOP
      TMH0=FACT0*TEMP
      MRH0=MAX( (TMH0/MREF(1) + 0.5) , 1.)
      IF(DS .GT. 0.0)   THEN
       TMSED=FACT1*TEMP*(DS*CC0)/(H0*CC1)
       MRSED=MAX( (TMSED/MREF(1) + 0.5) , 1.)
      ELSE
       TMSED=0.
       MRSED=0
      END IF
      DO 1200 I=1,8
        N(I)= MREF(I)*MRH0
c        write(*,*)'sspnap2,N(I)= MREF(I)*MRH0',N(I),MREF(I),MRH0
 1200 CONTINUE
C
      IF(DS.GT.0.0)   THEN
       DO 1300 I=1,8
       NS(I)= MREF(I)*MRSED
 1300  CONTINUE
      ELSE
       DO 1400 I=1,8
       NS(I)=0
 1400  CONTINUE
      END IF
 2100 CONTINUE
      IF(FLAGPU .LT.-1.0)
     1  WRITE(LUPRT,601) (JJ,N(JJ),NS(JJ),JJ=1,NMES)
 601  FORMAT(1H0,'CALCULATION MESHES',
     1      /1H ,'STEP   WATER   SEDIMENT',
     2     (/1H ,I3,I9,I9))
C
        DO 67 I=1,8
        H(I)=1.0D0/DFLOAT(N(I))
        DH0SQ(I)=1.0D0/DFLOAT(N(I)**2)
        IF(DSN.NE.0.0)   THEN
          HS(I)=DSN/DFLOAT(NS(I))
        ELSE 
          HS(I)=0.0
        END IF
   67 CONTINUE
      DRAT=1.D0/(1.D0+DSN)
C
C     INITIALIZE REFLECTION COEFFICIENT CALCULATION
C
      IF (C2S.GT.0.0) THEN
       IF (DSN.GT.0.0) THEN
        CINTFC=C1(ND1)
        RINTFC=R1
       ELSE
        CINTFC=C0(ND0)
        RINTFC=R0
       END IF
       CALL REFL1(C2,C2S,BETA(2),BETA(3))
      END IF
C
C     FIND THE FIRST APPROXIMATION
C
c      WRITE(LUPRT,*) ' FINDING FIRST APPROXIMATION '
      N1=N(1)
      NS1=NS(1)
      HS1=HS(1)
      IF(DSN.LT.1E-6) THEN
        HRAT=1.
        STIFF=1.
        ELSE
        HRAT=H(1)/HS1
        STIFF=(C0(ND0)/C1(1))**2/ROS
      ENDIF
      BOT=((OMEGA*H0)/(C2*N1))**2
      CON1=2.*(STIFF-HRAT**2/ROS)/(STIFF+HRAT)
      CON2=2.*HRAT/(STIFF+HRAT)
      CON3=1.D0/HRAT**2
      CON4=2./ROS*HRAT**2/(STIFF+HRAT)
      CON5=ROS/(ROB*HRAT)
c      write(*,*) 'ssnap2, n1' ,n1
      CALL ISOINT(ISO,MMAX,MMIN,H(1),N1,NS1,HS1,HRAT,CC0,F,SPEED,
     & ADA,LTOT,NP555,NMODES,OMEGA)
      IF(MMAX.LT.MMIN)   RETURN
      NMOD=MMAX-MMIN+1
      DO 1 M=1,NMOD
      MN=M+MMIN-1 
      A=ISO(M+1)
      B=ISO(M)
      CALL BRENT(A,B,N1,MYBR,NS1,HRAT,ADA,NP555,
     & OMEGA,H(1),HS1)
C
C
      MY(1,M)=MYBR/H(1)**2
      IF(DS.LE.0.0)  THEN
       MY(1,M)=MY(1,M)-((K**2
     & -(DSIN((DFLOAT(MN)-.5D0)*PI*H(1)*0.5D0)/H(1)*2.D0*DRAT)**2)
     &  -(K**2-((DFLOAT(MN)-.5D0)*PI*DRAT)**2))
      END IF
      MY(2,M)=MY(1,M)
   1  CONTINUE
C
C     FIND THE SUBSEQUENT APPROXIMATIONS
C
      DO 9 I=2,NMES
C
C     INITIAL GUESS
C
c      WRITE(LUPRT,*) ' FINDING APPROXIMATION FOR MESH NO. : ',I
C
      DO 7 M=1,NMOD
      MN=M+MMIN-1
      IF(DS.LE.0.0)  THEN
      MY(I,M)=MY(I,M)+((K**2
     &        -(DSIN((DFLOAT(MN)-.5D0)*PI*H(I)*.5D0)/H(I)*2.D0*DRAT)**2)
     &        -(K**2-((DFLOAT(MN)-.5D0)*PI*DRAT)**2))
      END IF
      MY(I,M)=MY(I,M)*H(I)**2
   7  CONTINUE
C
      NI=N(I)
      NSI=NS(I)
C
      CALL NEWTON(MY,NMOD,H,NI,I,HS,NSI,ADA,SPEED,NP555,NMODES,
     & CC0,OMEGA)
   
C
C
      IF(I.GT.2)   THEN
      CALL REVISE(MY,NMODES,NMOD,NMES,I,N,NS,*2100)
      END IF
C
C     LAST TIME FIND EIGENVECTORS AND LOSSES
C
      IF(I.EQ.NMES) THEN
C
C
c      WRITE(LUPRT,*) ' FINDING FINAL EIGENVECTORS AND LOSSES '
      MAXMOD= MMIN + NMOD - 1
c      WRITE(LUPRT,*)'MAX ORDER COMPUTED MODE : ',MAXMOD
c      IF( ((FLAGPU.le.0.0) )   WRITE(LUPRT,360)
      MMAX=0
        if (nsi.gt.0)then
          nsp=nsi+ni+2
        else
          nsp=ni+1
        endif
        do isp=1,nsp
          slow(isp)=1.0/speed(isp)
        enddo
        DO 15 M=1,NMOD
        CALL EIGVEC(MY,M,NI,NSI,xtS,slow,CC0,CC1,H(I),HS(I),
     &     F,ALFA(M),MSP,NMES,INVSTT,DRAT,MMIN,MMAX,EK,
     & EIGF,MODAVR,H,JF2,ADA,MODEN,NP555,DH0SQ)

      IF(MMAX.GT.0)   RETURN
  15    CONTINUE
      ELSE
C
C     CORRECTION
C
      DO 8 M=1,NMOD
      MN=M+MMIN-1
      MY(I,M)=MY(I,M)/H(I)**2
      IF(DS.LE.0.0)  THEN
      MY(I,M)=MY(I,M)-((K**2
     &        -(DSIN((DFLOAT(MN)-.5D0)*PI*H(I)*.5D0)/H(I)*2.D0*DRAT)**2)
     &        -(K**2-((DFLOAT(MN)-.5D0)*PI*DRAT)**2))
      END IF
   8  CONTINUE
C
C***      CALL RICH(MY,NMODES,NMOD,H,I,NMES)
      CALL LAGRANGE(1,NMOD,I,MY,NMODES,DH0SQ,MAXMSH,1)
C
      END IF
   9  CONTINUE
C
        MMAX=MMIN+NMOD-1
      RETURN
      END
C
C
C***********************************************************
C
      SUBROUTINE ISOINT(ISO,MMAX,MMIN,H,NMAX,NSMAX,HS,HRAT,CC0,
     & F,SPEED,ADA,LTOT,NP555,NMODES,OMEGA)
C
C___________________________________________________________
C                                                          |
C     This routine isolates the eigenvalues for the        |
C     first mesh by counting sign changes in the Sturm     |
C     sequence. The counting is done by the routine STURM  |
C__________________________________________________________|
C
      DOUBLE PRECISION OMEGA
      DOUBLE PRECISION K,EIGMAX,EIGVAL,ADA(NP555),PI,H,Z,
     &                 CW,CMIN,C0, CC0
      DOUBLE PRECISION ISO(NMODES),S,SPEED(NP555),BOT,HS,
     &                 CS,HRAT,STIFF, F
      DOUBLE PRECISION STEP

      COMMON /GEN/ PI,K,ROB,ROS,BOT,STIFF,CB
      COMMON /LUNIT/ MSOURC,MODOLD,MODNEW,LUPRT
C
  300 FORMAT(1H ,/    /,' SOURCE FREQUENCY     =',F8.2,' HZ',/,
     $ ' ESTIMATED MAXIMUM NO. OF MODES =',I5,/)
C
      DO 3 I=1,NMODES
   3  ISO(I)=0.D0
      CMIN=1.0D10
C
C     DEFINE MATRIX DIAGONAL
C
      CALL VELOCITY(NMAX,NSMAX,CC0,SPEED,NP555,H,HS)
      DO 1 N=1,NMAX
        CW=SPEED(N+1)
        IF(CW.LT.CMIN) CMIN=CW
        ADA(N)=(K*H/CW)**2-2.
   1  CONTINUE
      DO 4 N=NMAX+1,NSMAX+NMAX
        CS=SPEED(N+2)
        IF(CS.LT.CMIN) CMIN=CS
        ADA(N)=(K*HS/CS)**2-2.
   4  CONTINUE
C
C     FIND EIGMAX 
C
      EIGMAX=(K*H/CMIN)**2
C
C     FIND THE ESTIMATED NUMBER OF PROPAGATING MODES, MMAX
C

      EIGVAL=BOT
      CALL STURM(EIGVAL,NMAX,S,M,NSMAX,ADA,NP555,OMEGA,H,HS)
C
c      WRITE(LUPRT,300) F, M
      MMAX=MIN0(MMAX,M)
      IF((MMAX-MMIN).LT.0)   RETURN
C
        INDEX=MIN(NMODES,M-MMIN+2)
        ISO(INDEX)=EIGVAL
C
C     ISOLATE THE EIGENVALUES
C
      IF(MMIN.EQ.1) THEN
        ISO(1)=EIGMAX
        M1=1
        ELSE
        M1=MMIN-1
      ENDIF
  10  CONTINUE
      EIGVAL=(EIGMAX+BOT)*.5D0
      STEP=(EIGMAX-BOT)*.5D0
  11  CONTINUE
      CALL STURM(EIGVAL,NMAX,S,M,NSMAX,ADA,NP555,OMEGA,H,HS)
      STEP=STEP*.5D0
      INDEX=MIN(NMODES,M-MMIN+2)
      INDEX=MAX(1,INDEX)
      IF(INDEX.GT.1)   THEN
      IF(ISO(INDEX).EQ.0.D0.AND.M+1.GE.MMIN.AND.M.LE.MMAX) 
     &  ISO(INDEX)=EIGVAL
      ELSE
      ISO(1)=EIGVAL
      END IF
      IF(M.GT.M1) THEN
        EIGVAL=EIGVAL+STEP
        GO TO 11
        ELSEIF(M.LT.M1) THEN
        EIGVAL=EIGVAL-STEP
        GO TO 11
      ENDIF
  12  CONTINUE
      M1=M1+1
      IF(M1.EQ.MMAX+1) RETURN
      IF(ISO(M1-MMIN+2).EQ.0.D0) THEN
        EIGMAX=ISO(M1-MMIN+1)
        GO TO 10
        ELSE
        GO TO 12
      ENDIF
      END
C
C
C***********************************************************
C
      SUBROUTINE VELOCITY(NMAX,NSMAX,CC0,SPEED,NP555,H,HS)
C
C___________________________________________________________
C                                                          |
C     This routine calculates the speed of sound.          |
C__________________________________________________________|
C

      PARAMETER (NDEP=201)

      DOUBLE PRECISION SPEED(NP555),H,HS,CNORM,Z,DIV,CC0

      COMMON /A/ C0(NDEP),C1(NDEP),Z0(NDEP),Z1(NDEP),ND0,ND1

C
C     WATER
C
      CNORM=1D0/CC0
      I=2
      DIV=(C0(2)-C0(1))/(Z0(2)-Z0(1))
      SPEED(1)=C0(1)*CNORM
      DO 2000 N=2,NMAX
      Z=H*(N-1)
      IF(Z.GT.Z0(I)) THEN
 1000  CONTINUE  
       I=I+1
       IF(Z.GT.Z0(I))   GO TO 1000
       DIV=(C0(I)-C0(I-1))/(Z0(I)-Z0(I-1))
      END IF
      SPEED(N)=(DIV*(Z-Z0(I-1))+C0(I-1))*CNORM
 2000 CONTINUE
c      write(*,*)'ND0,NMAX',ND0,NMAX
      SPEED(NMAX+1)=C0(ND0)*CNORM
C
C     SEDIMENT
C
      IF(NSMAX.GT.0)   THEN
      I=2
      DIV=(C1(2)-C1(1))/(Z1(2)-Z1(1))
      SPEED(NMAX+2)=C1(1)*CNORM
      DO 4000 N=2,NSMAX
      Z=HS*(N-1)+1.D0
      IF(Z.GT.Z1(I)) THEN 
 3000  CONTINUE
       I=I+1
       IF(Z.GT.Z1(I))    GO TO 3000
       DIV=(C1(I)-C1(I-1))/(Z1(I)-Z1(I-1))
      END IF
      SPEED(N+NMAX+1)=(DIV*(Z-Z1(I-1))+C1(I-1))*CNORM
 4000 CONTINUE
      SPEED(NSMAX+NMAX+2)=C1(ND1)*CNORM
      END IF
      RETURN
      END
C
C
C***********************************************************
C
C
C
C**********************************************************
C
      SUBROUTINE BRENT(A,B,MH0I,X,MSEDI,HRAT,ADA,NP555,
     & OMEGA,DH0I,DSEDI)
C
C MICHAEL B. PORTER  8/84
C
C FORTRAN CONVERSION OF ALGOL PROGRAM PUBLISHED IN
C    THE COMPUTER JOURNAL 14(4):422-425 (1971)
C    BY R. P. BRENT
C
C RETURNS A ZERO X OF THE FUNCTION F IN THE GIVEN INTERVAL
C    [A,B], TO WITHIN A TOLERANCE 6*MACHEP*ABS(X)+2*T, WHERE
C    MACHEP IS THE RELATIVE MACHINE PRECISION AND T IS A POSITIVE
C    TOLERANCE.  THE PROCEDURE ASSUMES THAT FUNCT(A) AND FUNCT(B) HAVE
C    DIFFERENT SIGNS.
C
C    THIS IS THE EXTENDED RANGE VERSION WHICH EXPECT THE FUNCTION
C    TO HAVE THE FORM
C         SUBROUTINE FUNCT(X,G,IPOW)
C    WHERE G*10**IPOW GIVES FUNCT(X)
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DOUBLE PRECISION DH0I, DSEDI, HRAT, OMEGA
      DOUBLE PRECISION MACHEP,M,T,TREF
      DOUBLE PRECISION ADA(NP555)
      CHARACTER*80 ERRMSG
      COMMON /LUNIT/ MSOURC,MODOLD,MODNEW,LUPRT
      DATA TREF/1.0D-13/

      ERRMSG = ' '
      MACHEP = 1.0E-16
      TEN = 10.0
      T=TREF*A
C      CALL FUNCT(A,FA,IEXPA)
C      CALL FUNCT(B,FB,IEXPB)
      CALL CHARAC(A,MH0I,FA,MSEDI,NOVFLA,ADA,NP555,
     & OMEGA,DH0I,DSEDI)
      CALL CHARAC(B,MH0I,FB,MSEDI,NOVFLB,ADA,NP555,
     & OMEGA,DH0I,DSEDI)
      IEXPA=NOVFLA*20
      IEXPB=NOVFLB*20
      IF ( ( (FA .GT. 0.0) .AND. (FB .GT. 0.0) ) .OR.
     &     ( (FA .LT. 0.0) .AND. (FB .LT. 0.0) ) ) THEN
         ERRMSG = ' *** ZBRENT ERROR: FUNCT SGN SAME AT INTRVL ENDPTS'
         WRITE(6,*) ERRMSG
         STOP
C         RETURN
      ENDIF
C
C     INTERNAL ROOT
C
 2000 C = A
      FC = FA
      IEXPC = IEXPA
      E = B-A
      D = E
C
C     EXTERNAL ROOT
C
      IF (IEXPA .LT. IEXPB) THEN
         F1 = FC*TEN**(IEXPC-IEXPB)
         F2 = FB
      ELSE
         F1 = FC
         F2 = FB*TEN**(IEXPB-IEXPC)
      ENDIF
 3000 IF (ABS(F1) .LT. ABS(F2)) THEN
         A = B
         B = C
         C = A
         FA = FB
         IEXPA = IEXPB
         FB = FC
         IEXPB = IEXPC
         FC = FA
         IEXPC = IEXPA
      ENDIF
      TOL = 2.0*MACHEP*ABS(B)+T
      M = 0.5*(C-B)
      IF ((ABS(M) .GT. TOL) .AND. (FB .NE. 0.0)) THEN
C     ------ SEE IF A BISECTION IS FORCED
      IF (IEXPA .LT. IEXPB) THEN
         F1 = FA*TEN**(IEXPA-IEXPB)
         F2 = FB
      ELSE
         F1 = FA
         F2 = FB*TEN**(IEXPB-IEXPA)
      ENDIF
         IF ( (ABS(E) .LT. TOL) .OR. 
     &       (ABS( F1 ) .LE. ABS(F2)) ) THEN
            E = M
            D = E
         ELSE
            S = FB/FA*TEN**(IEXPB-IEXPA)
            IF (A .EQ. C) THEN
C              ------ LINEAR INTERPOLATION
               P = 2.0*M*S
               Q = 1.0-S
            ELSE
C              ------ INVERSE QUADRATIC INTERPOLATION
               Q = FA/FC*TEN**(IEXPA-IEXPC)
               R = FB/FC*TEN**(IEXPB-IEXPC)
               P = S*(2.0*M*Q*(Q-R)-(B-A)*(R-1.0))
               Q = (Q-1.0)*(R-1.0)*(S-1.0)
            ENDIF
            IF (P .GT. 0.0) THEN
               Q = -Q
            ELSE
               P = -P
            ENDIF
            S = E
            E = D
            IF ((2.0*P .LT. 3.0*M*Q-ABS(TOL*Q)) .AND.
     &          (P .LT. ABS(0.5*S*Q))) THEN
               D = P/Q
            ELSE
               E = M
               D = E
            ENDIF
         ENDIF
         A = B
         FA = FB
         IEXPA = IEXPB
         IF (ABS(D) .GT. TOL) THEN
            B = B+D
         ELSE
            IF (M .GT. 0.0) THEN
               B = B+TOL
            ELSE
               B = B-TOL
            ENDIF
         ENDIF
      CALL CHARAC(B,MH0I,FB,MSEDI,NOVFLB,ADA,NP555,
     & OMEGA,DH0I,DSEDI)
      IEXPB=NOVFLB*20
C         CALL FUNCT(B,FB,IEXPB)
         IF ((FB .GT. 0.0) .EQV. (FC .GT. 0.0)) GOTO 2000
         GOTO 3000
      ENDIF
      X = B
      RETURN
      END

      SUBROUTINE NEWTON(MY,MMAX,H,NMAX,I,HS,NSMAX,ADA,SPEED,NP555,
     & NMODES,CC0,OMEGA)
C
C____________________________________________________________
C                                                            |
C     This routine finds the zeros in the characteristic     |
C     equation for the second and subsequent meshes.         |
C____________________________________________________________|
C

      DOUBLE PRECISION MY(8,NMODES),K,H(8),PI,EPS,Z,CW,X,S,X0,X1,X2
      DOUBLE PRECISION F0,F1,HS(8),STIFF,CS,BOT,HRAT,C0,HI,HSI, CC0
      DOUBLE PRECISION F00,F11
      DOUBLE PRECISION OMEGA
      DOUBLE PRECISION MYM
      DOUBLE PRECISION CON1,CON2,CON3,CON4,CON5,ADA(NP555),
     &                 SPEED(NP555)

      COMMON /G/ R0,R1,R2,C2,H0,H1,FACT0,FACT1,C11
      COMMON /GEN/ PI,K,ROB,ROS,BOT,STIFF,CB
      COMMON /CONST/ CON1,CON2,CON3,CON4,CON5
      COMMON /LUNIT/ MSOURC,MODOLD,MODNEW,LUPRT
      EPS=1.D-11
C
C     DEFINE MATRIX DIAGONAL
C
      HI=H(I)
      HSI=HS(I)
      CALL VELOCITY(NMAX,NSMAX,CC0,SPEED,NP555,HI,HSI)
      DO 1 N=1,NMAX
        CW=SPEED(N+1)
        ADA(N)=(K*H(I)/CW)**2-2.D0
   1  CONTINUE
      NFROM=NMAX+1
      NTO=NSMAX+NMAX
      DO 5 N=NFROM,NTO
        CS=SPEED(N+2)
        ADA(N)=(K*HS(I)/CS)**2-2.
   5  CONTINUE
      BOT=((OMEGA*H0*H(I))/C2)**2
      IF(NSMAX.EQ.0) THEN
        HRAT=1.
      ELSE
        HRAT=H(I)/HS(I)
      ENDIF
      CON1=2.*(STIFF-HRAT**2/ROS)/(STIFF+HRAT)
      CON2=2.*HRAT/(STIFF+HRAT)
      CON3=1.D0/HRAT**2
      CON4=2./ROS*HRAT**2/(STIFF+HRAT)
      CON5=ROS/(ROB*HRAT)
C
C     FIND THE MMAX EIGENVALUES
C
      DO 2 M=1,MMAX
        X0=MY(I,M)
        X1=X0+10.D0*EPS*X0
        CALL CHARAC(X0,NMAX,F0,NSMAX,NOVFL0,ADA,NP555,
     &   OMEGA,HI,HSI)
        NITER=0
   3  CONTINUE
        CALL CHARAC(X1,NMAX,F1,NSMAX,NOVFL1,ADA,NP555,
     &      OMEGA,HI,HSI)
        NODIF=NOVFL1-NOVFL0
        F00=F0
        F11=F1
        IF (NODIF.NE.0) THEN
          IF (NODIF.GT.0) THEN
            F00=F0*(1E-20)**NODIF
          ELSE
            F11=F1*(1E-20)**(-NODIF)
          END IF
        END IF
        IF(F11-F00.EQ.0.D0) THEN
          X2=(X0+X1)*0.5d0
        ELSE
          X2=(X0*F11-X1*F00)/(F11-F00)
        ENDIF
        NITER=NITER+1
       IF (NITER.GT.10) THEN
c        WRITE(LUPRT,*) '*** MAX ITERATIONS EXCEEDED IN NEWTON***'
        GO TO 31
      END IF
      IF(ABS(X1-X2).GT.ABS(EPS*X2)) THEN
        X0=X1
        F0=F1
        NOVFL0=NOVFL1
        X1=X2
        GO TO 3
      ENDIF
 31   MY(I,M)=X2
      IF(X2.LT.BOT) THEN
        NMOD=M-1
        GO TO 4
      ENDIF
   2  CONTINUE
      NMOD=MMAX
   4  CONTINUE
      MMAX=NMOD
      RETURN
      END















