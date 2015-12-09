C**********************************************************
C
      SUBROUTINE STURM(EIGVAL,NMAX,S,M,NSMAX,ADA,NP555,
     $ OMEGA,DH0I,DSEDI )
C
C__________________________________________________________
C                                                         |
C     This routine counts the number of sign changes in   |
C     the Sturm sequence and calculates the value of the  |
C     characteristic equation.                            |
C_________________________________________________________|  
C
      REAL*8 U1H0,U1SED,DH0I,DSEDI
      REAL*8 SEDK,OMEGA
      REAL*8 EIGVAL,ADA(NP555),S0,S1,S2,S,BOT,STIFF,K,PI
      REAL*8 CON1,CON2,CON3,CON4,CON5,TEMP
      COMMON/CONST/ CON1,CON2,CON3,CON4,CON5
      COMMON/EXPMAX/TRESH
      COMMON/G/ R0,R1,R2,C2,H0,H1,FACT0,FACT1,C11
      COMMON/GEN/ PI,K,ROB,ROS,BOT,STIFF,CB
      DATA SCFAC/1.0E-20/
C
C     'M' COUNTS SIGN CHANGES
C
C     WATER LAYERS
C
      S0=1.D0
      S1=(EIGVAL-ADA(1))*S0
      M=0
      DO 1000   N=2,NMAX-1
      S2=(EIGVAL-ADA(N))*S1-S0
      IF( ((S2 .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S2 .LT. 0.0))     )   M=M+1
      IF( ABS(S2) .GE. TRESH )   THEN
       S0=S1*SCFAC
       S1=S2*SCFAC
      ELSE 
       S0=S1
       S1=S2
      END IF
 1000 CONTINUE
C
C     SEDIMENT LAYERS
C
      IF(NSMAX.GT.0) THEN
C
C
C
C   FICTITIOUS POINT FROM WATER TO SEDIMENT 
       S2=(EIGVAL - ADA(NMAX))*S1 - S0
C   DERIVATIVE AT S1
       U1H0=(S2-S0)/(2.0*DH0I)
       U1SED=U1H0
C
C  INTERFACE POINT AS SEEN FROM THE SEDIMENT
       S1=S1/ROS
C
       SEDK=((OMEGA*DBLE(H1))/(DBLE(C11)*DFLOAT(NSMAX)))**2
C      S2=(U1SED*2*DSEDI+2*S1-(SEDK-EIGVAL*CON3)*S1)/2.0
       S2=U1SED*DSEDI - 0.5*(SEDK-EIGVAL*CON3-2.0D0)*S1
C
      IF( ((S2 .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S2 .LT. 0.0))     )   M=M+1
       IF( ABS(S2) .GE. TRESH )   THEN 
        NOVFL=NOVFL + 1
        S0=S1*SCFAC
        S1=S2*SCFAC
       ELSE 
        S0=S1
        S1=S2
       END IF
C
C
       DO 2000   N=NMAX+1,NMAX+NSMAX-1
      S2=(EIGVAL*CON3-ADA(N))*S1-S0
      IF( ((S2 .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S2 .LT. 0.0))     )   M=M+1
      IF( ABS(S2) .GE. TRESH )   THEN 
       S0=S1*SCFAC
       S1=S2*SCFAC
      ELSE 
       S0=S1
       S1=S2
      END IF
 2000 CONTINUE
C
      END IF
C
C     BOTTOM
C
      IF(BOT.GT.EIGVAL) THEN
        S=(EIGVAL*CON3-ADA(NMAX+NSMAX))*S1-2.*S0
      ELSE
        S=(EIGVAL*CON3-ADA(NMAX+NSMAX)
     &         +2.D0*SQRT(EIGVAL-BOT)*CON5)*S1-2.*S0
      ENDIF
      IF( ((S .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S .LT. 0.0))     )   M=M+1
      RETURN
      END
C
C
C************************************************************ 	
C
      SUBROUTINE CHARAC(EIGVAL,NMAX,S,NSMAX,NOVFL,ADA,NP555,
     $ OMEGA,DH0I,DSEDI)
C
C_____________________________________________________________
C                                                             |
C     This routine calculates the value of the characteristic |
C     equation of the matrix.                                 |
C_____________________________________________________________|
C
      REAL*8 U1H0,U1SED,DH0I,DSEDI
      REAL*8 SEDK,OMEGA,ztemp
      REAL*8 EIGVAL,ADA(NP555),S1,S2,S0,S,s3,BOT,STIFF,K,PI
      REAL*8 CON1,CON2,CON3,CON4,CON5
      COMMON/CONST/ CON1,CON2,CON3,CON4,CON5
      COMMON/EXPMAX/TRESH
      COMMON/G/ R0,R1,R2,C2,H0,H1,FACT0,FACT1,C11
      COMMON/GEN/ PI,K,ROB,ROS,BOT,STIFF,CB
      DATA SCFAC/1.0E-20/ 
C
C     WATER LAYERS
C
      NOVFL=0
      S0=1.D0
      S1=(EIGVAL-ADA(1))
      DO 1000   N=2,NMAX-1,2
        S2=(EIGVAL-ADA(N))*S1-S0
c          S0=S1
c          S1=S2
        S3=(EIGVAL-ADA(N+1))*S2-S1
C
        IF( ABS(S2) .LE. TRESH )   THEN 
          S0=S2
          S1=S3
        ELSE 
          NOVFL=NOVFL + 1
          S0=S2*SCFAC
          S1=S3*SCFAC
        END IF
C
 1000 CONTINUE
C
C    
C
C     SEDIMENT LAYERS
C
      IF(NSMAX.GT.0) THEN
C
C   FICTITIOUS POINT FROM WATER TO SEDIMENT 
       S2=(EIGVAL - ADA(NMAX))*S1 - S0
C   DERIVATIVE AT S1
       U1H0=(S2-S0)/(2.0*DH0I)
       U1SED=U1H0
C
C  INTERFACE POINT AS SEEN FROM THE SEDIMENT
       S1=S1/ROS
C
       SEDK=((OMEGA*DBLE(H1))/(DBLE(C11)*DFLOAT(NSMAX)))**2
       S2=U1SED*DSEDI - 0.5*(SEDK-EIGVAL*CON3-2.0D0)*S1
C
       IF( ABS(S2) .GE. TRESH )   THEN 
         NOVFL=NOVFL + 1
         S0=S1*SCFAC
         S1=S2*SCFAC
       ELSE 
         S0=S1
         S1=S2
       END IF
C
C
       ztemp=eigval*con3
       ntemp=nmax+nsmax-1
       DO 2000   N=NMAX+1,ntemp
       S2=(ztemp-ADA(N)  )*S1-S0
c          S0=S1
c          S1=S2
c       S2=(ztemp-ADA(N+1))*S1-S0
C
      IF( ABS(S2) .GE. TRESH )   THEN 
       NOVFL=NOVFL + 1
       S0=S1*SCFAC
       S1=S2*SCFAC
      ELSE 
       S0=S1
       S1=S2
      END IF
C
 2000 CONTINUE
      ENDIF
C
C     BOTTOM
C
      IF(BOT.GT.EIGVAL) THEN
        S=(EIGVAL*CON3-ADA(NMAX+NSMAX))*S1-2.*S0
      ELSE
        S=(EIGVAL*CON3-ADA(NMAX+NSMAX)
     &        +2.D0*SQRT(EIGVAL-BOT)*CON5)*S1-2.D0*S0
      ENDIF
      RETURN
      END
C
      SUBROUTINE CHARACcopy(EIGVAL,NMAX,S,NSMAX,NOVFL,ADA,NP555,
     $ OMEGA,DH0I,DSEDI)
C
C_____________________________________________________________
C                                                             |
C     This routine calculates the value of the characteristic |
C     equation of the matrix.                                 |
C_____________________________________________________________|
C
      REAL*8 U1H0,U1SED,DH0I,DSEDI
      REAL*8 SEDK,OMEGA
      REAL*8 EIGVAL,ADA(NP555),S1,S2,S0,S,BOT,STIFF,K,PI
      REAL*8 CON1,CON2,CON3,CON4,CON5
      COMMON/CONST/ CON1,CON2,CON3,CON4,CON5
      COMMON/EXPMAX/TRESH
      COMMON/G/ R0,R1,R2,C2,H0,H1,FACT0,FACT1,C11
      COMMON/GEN/ PI,K,ROB,ROS,BOT,STIFF,CB
      DATA SCFAC/1.0E-20/ 
C
C     WATER LAYERS
C
      NOVFL=0
      S0=1.D0
      S1=(EIGVAL-ADA(1))*S0
      DO 1000   N=2,NMAX-1
        S2=(EIGVAL-ADA(N))*S1-S0
C
        IF( ABS(S2) .GE. TRESH )   THEN 
          NOVFL=NOVFL + 1
          S0=S1*SCFAC
          S1=S2*SCFAC
        ELSE 
          S0=S1
          S1=S2
        END IF
C
 1000 CONTINUE
C
C    
C
C     SEDIMENT LAYERS
C
      IF(NSMAX.GT.0) THEN
C
C   FICTITIOUS POINT FROM WATER TO SEDIMENT 
       S2=(EIGVAL - ADA(NMAX))*S1 - S0
C   DERIVATIVE AT S1
       U1H0=(S2-S0)/(2.0*DH0I)
       U1SED=U1H0
C
C  INTERFACE POINT AS SEEN FROM THE SEDIMENT
       S1=S1/ROS
C
       SEDK=((OMEGA*DBLE(H1))/(DBLE(C11)*DFLOAT(NSMAX)))**2
       S2=U1SED*DSEDI - 0.5*(SEDK-EIGVAL*CON3-2.0D0)*S1
C
       IF( ABS(S2) .GE. TRESH )   THEN 
         NOVFL=NOVFL + 1
         S0=S1*SCFAC
         S1=S2*SCFAC
       ELSE 
         S0=S1
         S1=S2
       END IF
C
C
       DO 2000   N=NMAX+1,NMAX+NSMAX-1
       S2=(EIGVAL*CON3-ADA(N))*S1-S0
C
      IF( ABS(S2) .GE. TRESH )   THEN 
       NOVFL=NOVFL + 1
       S0=S1*SCFAC
       S1=S2*SCFAC
      ELSE 
       S0=S1
       S1=S2
      END IF
C
 2000 CONTINUE
      ENDIF
C
C     BOTTOM
C
      IF(BOT.GT.EIGVAL) THEN
        S=(EIGVAL*CON3-ADA(NMAX+NSMAX))*S1-2.*S0
      ELSE
        S=(EIGVAL*CON3-ADA(NMAX+NSMAX)
     &        +2.D0*SQRT(EIGVAL-BOT)*CON5)*S1-2.D0*S0
      ENDIF
      RETURN
      END
C
C
C**********************************************************
C
      SUBROUTINE EIGVEC(MY,MEIG,NMAX,NSMAX,xtS,Slow,CC0,CC1,
     &  H,HS,F,ALFA,MSP,NMES,INVSTT,DRAT,MMIN,MMAX,EK,
     &  EIGF,MODAVR,HH,JF2,ADA,MODEN,NP555,DH0SQ)
C__________________________________________________________
C                                                         |
C     This routine calculates the eigenvectors by         |
C     means of inverse iteration. The eigenfunctions      |
C     are normalized, and the attenuation coefficients    |
C     are determined.                                     |
C_________________________________________________________|
C
      PARAMETER (NMODES=1420,NBEG=5,NPOINT=44500+5*NBEG*NMODES+2)
      PARAMETER (MAXMSH=8)
      REAL*8 EIGVAL,DEL,P,ZZMAX, ABS1, ABS2, CC0,zzinv
      REAL*8 NORMW, NORMS, NORMB, NORM, SQNORM
      REAL*8 ATTW, ATTS, ATTB
      REAL*8 MY(8,NMODES), WORK, DRAT
      REAL*8 HH(8), ADA(NP555), DH0SQ(9)
      DOUBLE PRECISION  CC1
C     LOCAL ARRAYS USED FOR INVERSE ITERATION
C
      REAL*8 Q,A3(NPOINT),B3(NPOINT),C3(NPOINT),DD,OMEGAC
      REAL*8 EE(NPOINT),ZZ(NPOINT)
      REAL*8 SSOLD(NPOINT),EIGF(NP555),DIFF
      LOGICAL EXCH(NPOINT)
C     COMMON /INVIT/ Q,A3,B3,C3,DD,EE,ZZ,SSOLD,S,EXCH
C
C
      REAL*8 F,f2,PI,K,KM,CIN,COST,BOT,STIFF,H,HS,Slow(NP555)
      REAL*8 CON1,CON2,CON3,CON4,CON5,EK(MODEN)
      REAL K0,K1,K2,KH0, MODAVR(MODEN),      xtS(moden,MSP)
      COMMON /ATTEN/ ALF0,ALF1,ALF2,ALF2S,ALFOS,ALFOB
      COMMON/AB/ BETA(3),SCATT(2),CBS
      COMMON/CONST/ CON1,CON2,CON3,CON4,CON5
      COMMON/FLAGPULSE/FLAGPU
      COMMON/G/ R0,R1,R2,C2,D,DS,FACT0,FACT1,C11
      COMMON/GEN/ PI,K,ROB,ROS,BOT,STIFF,CB
      COMMON/LUNIT/MSOURC,MODOLD,MODNEW,LUPRT
      COMMON /XREFL/ CINTFC,RINTFC
      real xnul
      real*8 dnul
      data xnul/0.e0/,dnul/0.d0/

C     EIGVAL=SQRT(MY(NMES,MEIG))/(H*D)
C
      EIGVAL=MY(NMES,MEIG)
      KM=SQRT(EIGVAL)/H
C
      MN=MEIG+MMIN-1
      MY(NMES,MEIG)=MY(NMES,MEIG)/H**2
      IF(DS.LE.0.0) THEN
      MY(NMES,MEIG)=MY(NMES,MEIG)-((K**2
     &        -(DSIN((DFLOAT(MN)-.5D0)*PI*H*.5D0)/H*2.D0*DRAT)**2)
     &        -(K**2-((DFLOAT(MN)-.5D0)*PI*DRAT)**2))
      END IF

      DH0SQ(NMES+1)=0.0D0
      CALL LAGRANGE(MEIG,MEIG,NMES,MY,NMODES,DH0SQ,MAXMSH,1-NMES)
      WORK=MY(1,MEIG)
C
      IF(WORK.LT.(K/CB)**2)   THEN
      MMAX=MN-1
      RETURN
      END IF
C
      MY(NMES,MEIG)=WORK
      EK(MEIG)=SQRT(WORK)/D
C
c      IF(JF2.GT.0)   RETURN
      OMEGAC=2.0*PI*F/CC0
C
C     SET UP OF COEFFICIENT MATRIX
C
C     WATER LAYERS
C
 2020 CONTINUE
      A3(1)=ADA(1)-EIGVAL
      B3(1)=1.0D0
      DO 1 N=2,NMAX-1
        A3(N)=ADA(N)-EIGVAL
        B3(N)=1.0D0
        C3(N-1)=1.0D0
 1    CONTINUE
C
C     SEDIMENT LAYERS
C
      IF (NSMAX.GT.0) THEN
        A3(NMAX)=ADA(NMAX)-EIGVAL + (H/(HS*ROS))*
     &  ( (((OMEGAC*D*HS)*slow(NMAX+2))**2-2.) - EIGVAL*CON3 )
        B3(NMAX)=2.0D0*(H/HS)
        C3(NMAX-1)=2.0D0
        A3(NMAX+1)=ADA(NMAX+1)-EIGVAL*CON3
        B3(NMAX+1)=1.0D0
        C3(NMAX)=1.0D0/ROS
        DO 2 N=NMAX+2,NMAX+NSMAX-1
          A3(N)=ADA(N)-EIGVAL*CON3
          B3(N)=1.0D0
          C3(N-1)=1.0D0      
   2    CONTINUE
      END IF
      IF (BOT.GT.EIGVAL) THEN
        A3(NMAX+NSMAX)=ADA(NMAX+NSMAX)-EIGVAL*CON3
      ELSE
        A3(NMAX+NSMAX)=ADA(NMAX+NSMAX)-EIGVAL*CON3-2.0D0*SQRT(
     &  EIGVAL-BOT)*CON5
      END IF
      C3(NMAX+NSMAX-1)=2.0D0
      B3(NMAX+NSMAX)=dnul
      C3(NMAX+NSMAX)=dnul
C
C     ELIMINATION IN COEFFICIENT MATRIX
C
      NTOT=NMAX+NSMAX
      DO 20 I=1,NTOT-1
        i1=i+1
        IF (ABS(C3(I)).GT.ABS(A3(I))) THEN
          EXCH(I)=.TRUE.
          Q=C3(I)
          C3(I)=A3(I)
          A3(I)=Q
          Q=A3(I1)
          A3(I1)=B3(I)
          B3(I)=Q
          DD=B3(I1)
          A3(I)=1.0D0/A3(I)
          Q=-C3(I)*A3(I)
          EE(I)=Q
          A3(I1)=A3(I1)+Q*B3(I)
          B3(I1)=Q*DD
          C3(I) =DD
        ELSE
          EXCH(I)=.FALSE.
          A3(I)=1.0D0/A3(I)
          Q=-C3(I)*A3(I)
          EE(I)=Q
          A3(I1)=A3(I1)+Q*B3(I)
          C3(I) =Dnul
        END IF
 20   CONTINUE
      IF(A3(NTOT).EQ.0.0D0)   THEN
c       WRITE(LUPRT,*)'**** INVERSE ITERATION RESTARTED WITH ',
c     & 'OFFSET EIGENVALUE'
        EIGVAL=EIGVAL*(1.0D0+1.0D-12)
        GO TO 2020
      END IF
      A3(NTOT)=1.0D0/A3(NTOT)
C
C     ELIMINATION OF RIGHT HAND SIDE
C
      DO 25 I=1,NTOT
        ZZ(I)=1.0D0
c        SSOLD(I)=1.0D0
 25   CONTINUE
      NIT=0
 30   IF (NIT.GT.0) THEN
        DO 35 I=1,NTOT-1
        IF (EXCH(I)) THEN
          Q=ZZ(I+1)
          ZZ(I+1)=ZZ(I)
          ZZ(I)=Q
        END IF
        ZZ(I+1)=ZZ(I+1)+EE(I)*ZZ(I)
 35     CONTINUE
      END IF
C
C     BACK SUBSTITUTION
C
c       xdum=10e10
c       xdummy=xdum**100
      ZZMAX=dnul
      EIGF(NTOT)=ZZ(NTOT)*A3(NTOT)
      EIGF(NTOT-1)=(ZZ(NTOT-1)-B3(NTOT-1)*EIGF(NTOT))*A3(NTOT-1)
      ABS1=ABS(EIGF(NTOT))
      ABS2=ABS(EIGF(NTOT-1))
      ZZMAX=MAX(ABS1,ABS2)
      DO 40 I=NTOT-2,1,-1
        EIGF(I)=(ZZ(I)-B3(I)*EIGF(I+1)-C3(I)*EIGF(I+2))*A3(I)
        ZZMAX=MAX(ZZMAX,ABS(eigf(i)))
 40   CONTINUE
C
C     SCALE AND CHECK ACCURACY
C
      NIT=NIT+1
      Q=dnul
      zzinv=1.0d0/zzmax
      DO 50 I=1,NTOT
        EIGF(I)=EIGF(I)*ZZinv
        q=max(q,ABS(ABS(EIGF(I))-SSOLD(I)))
 50   CONTINUE
      IF (Q.LT.1.0D-10) GO TO 100
      IF (NIT.GT.10) THEN 
         WRITE(LUPRT,*) 
     &  '*** MAX NUMBER OF ITERATIONS EXCEEDED IN EIGVEC ***'
        WRITE(LUPRT,*)'*** MAX DIFFERENCE: ',Q
      GO TO 100
      END IF
      DO 55 I=1,NTOT
        ZZ(I)=EIGF(I)
        SSOLD(I)=abs(EIGF(I))
 55   CONTINUE
      GO TO 30
 100  CONTINUE
C
C     NORMALIZATION
C
      ATTW=0.D0
      NORMW=0.D0
      DO 3030 N=1,NMAX-1,2
        ATTW=ATTW + 2.0*EIGF(N)**2*slow(N+1) +
     &                  EIGF(N+1)**2*slow(N+2)
        NORMW=NORMW + 2.0*EIGF(N)**2 + EIGF(N+1)**2
 3030 CONTINUE
      ATTW=2.*ATTW - EIGF(NMAX)**2*slow(NMAX+1)
      NORMW=2.*NORMW - EIGF(NMAX)**2
      ATTW=ATTW/3.0*H
      NORMW=NORMW/3.0*H
      IF(NSMAX.GT.0)   THEN
        ATTS=EIGF(NMAX)**2/(ROS*ROS)*Slow(NMAX+2)
        NORMS=EIGF(NMAX)**2/(ROS*ROS)
        DO 3040 N=NMAX+1,NMAX+NSMAX-1,2
          ATTS=ATTS + 4.0*EIGF(N)**2*Slow(N+2) +
     &                2.0*EIGF(N+1)**2*Slow(N+3)
          NORMS=NORMS + 4.0*EIGF(N)**2 + 2.0*EIGF(N+1)**2
 3040   CONTINUE
        ATTS=ATTS - EIGF(NMAX+NSMAX)**2*slow(NSMAX+NMAX+2)
        NORMS=NORMS - EIGF(NMAX+NSMAX)**2
        ATTS=ATTS/3.0*HS*ROS
        NORMS=NORMS/3.0*HS*ROS
      ELSE
        ATTS=dnul
        NORMS=dnul
      END IF
      NORMB=H*ROB*(ROS/ROB*EIGF(NMAX+NSMAX))**2*.5/SQRT(EIGVAL-BOT)
      ATTB=NORMB/CB
      NORM=NORMW+NORMS+NORMB
      SQNORM=1.0/SQRT(D*NORM)
C
      zzinv=ek(meig)/omegac
      DO 3044   N=1,NTOT
        IF(slow(N).LE.zzinv)   GO TO 3044
        IF(MOD(MEIG+MMIN-1,2).EQ.0)   THEN
C         ODD ORDER MODES
          SQNORM=-SIGN(SQNORM,EIGF(N))
        ELSE
C         EVEN ORDER MODES
        SQNORM=SIGN(SQNORM,EIGF(N))
      END IF
      GO TO 3046
 3044 CONTINUE
 3046 CONTINUE
C
      DO 3050 N=1,NMAX+NSMAX
        EIGF(N)=EIGF(N)*SQNORM 
 3050 CONTINUE
      MODAVR(MEIG)=SQRT(NORMW)*ABS(SQNORM)
C
C     SAMPLING THE EIGENFUNCTIONS
C
      xts(meig,1)=0.
      DEL=1.0D0/DFLOAT(MSP-1)
      DO 10 I=2,MSP
        P=(I-1)*DEL*DFLOAT(NMAX)
        N1=INT(P)
        N2=N1+1
        IF(N1.EQ.0) THEN
         xtS(meig,I)=(P-DFLOAT(N1))*EIGF(N2)
        ELSE
         xtS(meig,I)=(P-DFLOAT(N1))*(EIGF(N2)-EIGF(N1))+EIGF(N1)
        ENDIF
  10  CONTINUE
      xtS(meig,MSP)=EIGF(NMAX)
C
C     CALCULATE THE MODAL ATTENUATION
C
C     COMPRESSIONAL ATTENUATION
C
      ATTW=ATTW/NORM*K/KM
      ATTS=ATTS/NORM*K/KM
      ATTB=ATTB/NORM*K/KM
      F2= F * F * 1.0D-6
      BETA0= F2 * (0.11 / (1.0 + F2) + 44.0 / (4100.0 + F2) +
     &       3.0D-4)
      BETA0= BETA0 * 1.0D-3 * CC0 / F
      ALF0= BETA0 * F * ATTW / (CC0 * 8.68589)
      IF(NSMAX.GT.0)   THEN
      ALF1=BETA(1)*F*ATTS/(CC1*8.68589)
      ELSE
      ALF1=0.0
      END IF
      ALF2=BETA(2)*F*ATTB/(CB*CC0*8.68589)
C
C     BOTT0M LOSSES
C
      IF(CBS.GT.0.) THEN
C
C     REFLECTION COEFFICIENT
C
        IF (NSMAX.GT.0) THEN
          CIN=1.0d0/Slow(NSMAX+NMAX+2)
        ELSE
          CIN=1.0d0/Slow(NMAX+1)
        END IF
        COST=KM/K*CIN
        CALL REFL2(CINTFC,RINTFC,R2,COST,Q)
C
C
C     THE FOLLOWING IS TAKEN FROM SNAP, BUT I DOUBT
C     IT'S CORRECTNESS FOR (K/CINTFC)<KM, I.E. THE
C     CASE WHERE THE MODE IS EVANESCENT AT THE
C     SUBBOTTOM INTERFACE
C
C     H. SCHMIDT, 841002
C
        K1=ABS((K/CIN)**2-KM*KM)
        SQRK1=SQRT(K1)
        K2=KM*KM-(K/CB)**2
        ALF2S=RINTFC*EIGF(NMAX+NSMAX)**2*SQRK1/(8.0*KM)
        ALF2S=ALF2S*(1.0+(RINTFC/R2)**2*K2/K1)*Q
        ALFA=ALF0+ALF1+ALF2S
      ELSE
        ALF2S=0.0
        ALFA=ALF0+ALF1+ALF2
      ENDIF
C
C     SCATTER LOSSES
C
C
      DER0=EIGF(1)/(D*H)
      IF (NSMAX.GT.0) THEN
      DER1=(CON1*EIGF(NMAX)-(CON1+CON4)*EIGF(NMAX-1)+
     & CON4*ROS*EIGF(NMAX+1))/(2.0*H*D)
      ELSE
      K1=KM*KM-(K/CB)**2
      SQK1=SQRT(K1)
      DER1=-EIGF(NMAX)*SQK1/(ROB*D)
      END IF
      SQKM=KM*KM
      RATS=(K*slow(1))**2
      RATB=(K*slow(NMAX+1))**2
      IF (SQKM.LT.RATS) THEN
      K0=SQRT(RATS-SQKM)
      ELSE
      K0=0.0
      END IF
      IF (SQKM.LT.RATB) THEN
      KH0=SQRT(RATB-SQKM)
      ELSE
      KH0=0.0
      END IF
      ALFOS=SCATT(1)**2*K0*DER0**2/(2.*KM)
      ALFOB=SCATT(2)**2*KH0/(2.*KM)*(DER1**2+(EIGF(NMAX)*KH0/D)**2)
      ALFA=ALFA+ALFOB+ALFOS
C
      CHECKS=2.0*((RATS-SQKM)/(D*D))*SCATT(1)**2
      CHECKB=2.0*((RATB-SQKM)/(D*D))*SCATT(2)**2
C
      IF(FLAGPU .LT.-1.0)   WRITE(LUPRT,376)MN,EK(MEIG),ALFA,ALF0,
     & ALF1,ALF2,ALF2S,ALFOS,ALFOB,Q,CHECKS,CHECKB
C
  376 FORMAT(1X,'MODE=',I4, 1X,D17.10,1X,E10.4,6(1X,E9.3),
     $2X,2(F5.2,1X),F5.2)
      RETURN
      END
C
C
C*************************************************************
C
      SUBROUTINE LAGRANGE(M1,M2,NMESH,MY,MODEN,DH0SQ,MAXMSH,IOFF)
C
      REAL*8 MY(MAXMSH,MODEN)
      REAL*8 DH0SQ(8), ZL, X, Y
      COMMON /LUNIT/ MSOURC,MODOLD,MODNEW,LUPRT
C
  100 FORMAT(1X,'*** EXTRAP NOT PERFORMED BECAUSE OF ARRAY SIZE',
     & ' LIMITATION ***')

      IF(NMESH+1 .GT. MAXMSH)   THEN
       WRITE(LUPRT,100)
       RETURN
      END IF
      IF(NMESH .GE. 2)   THEN
       X=DH0SQ(NMESH+1)
       DO 3000   M=M1,M2
       Y=0.0
       DO 2000  J=1,NMESH
       ZL=1.0
       DO 1000  K=1,NMESH
       IF(K.NE.J)  ZL=ZL*(X-DH0SQ(K))/(DH0SQ(J)-DH0SQ(K))
 1000  CONTINUE
       Y=Y+ZL*MY(J,M)
 2000  CONTINUE
       MY(NMESH+IOFF,M)= Y 
 3000  CONTINUE
      ELSE
       DO 4000   M=M1,M2
       MY(2,M)=MY(1,M)
 4000  CONTINUE
      END IF
      RETURN
      END





