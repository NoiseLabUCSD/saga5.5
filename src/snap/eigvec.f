      SUBROUTINE EIGVEC(MEIG,MH0I,MSEDI,XTS,
     &    MMM,Slow,CREF,CC0,
     &  CC1,DH0I,DSEDI,FRQ,ALFA,NMES,MINMOD,MAXMOD,MY,
     &  EIGF,MODAVR,ADA,A3,B3,C3,EE,ZZ,SSOLD,EXCH,EKM,EIGVAL)
C__________________________________________________________
C                                                         |
C     This routine calculates the eigenvectors by         |
C     means of inverse iteration. The eigenfunctions      |
C     are normalized, and the attenuation coefficients    |
C     are determined.                                     |
C_________________________________________________________|
C
      INTEGER EXTPOL
      DOUBLE PRECISION CREF, CC0, CC1, H0, H1
      DOUBLE PRECISION EKM, EIGVAL
      DOUBLE PRECISION DEL, PZ, ZZMAX, ABS1, ABS2,zzinv
      DOUBLE PRECISION NORMW, NORMS, NORMB, NORM, SQNORM
      DOUBLE PRECISION ATTW, ATTS, ATTB
      DOUBLE PRECISION ADA(NPOINT)
C
C     LOCAL ARRAYS USED FOR INVERSE ITERATION
C
      LOGICAL EXCH(NPOINT)

      REAL MODAVR(MODEN), XTS(moden,MSP)
      REAL K0, K1, K2, KH0

      DOUBLE PRECISION Q, DD, OMEGAC, DIFF
      DOUBLE PRECISION A3(NPOINT), B3(NPOINT), C3(NPOINT), 
     &                 EE(NPOINT), ZZ(NPOINT)
      DOUBLE PRECISION SSOLD(NPOINT), EIGF(NPOINT)
      DOUBLE PRECISION slow(NPOINT)
      DOUBLE PRECISION ROB, ROS
      DOUBLE PRECISION KM, CIN, COST, STIFF, DH0I, DSEDI
      DOUBLE PRECISION FRQ, TWOPI, PI, OMEGA
      DOUBLE PRECISION EIGREF, EIGMIN, EIGMAX
      DOUBLE PRECISION CON1, CON2, CON3, CON4, CON5, SEDK
      DOUBLE PRECISION MY(MAXMSH,MODEN)

      COMMON /ATTEN/ ALF0, ALF1, ALF2, ALF2S, ALFOS, ALFOB
      COMMON /AB/ BETA(-1:3), SCATT(2), CBS, C2
      COMMON /CONST/ CON1, CON2, CON3, CON4, CON5, SEDK
      COMMON /DENS/ R0, R1, R2
      COMMON /DENS8/ ROB, ROS
      integer flagpu
      COMMON /FLAGPU/ FLAGPU
       COMMON /FLAGPULSE/  EXTPOL, CORREC
      COMMON /GSNAP/ H0, H1, TWOPI, PI, OMEGA
      COMMON /GEN/ EIGREF, EIGMIN, EIGMAX, STIFF
      COMMON /LUNIT/ LUPLP, LUPLT, LUPRT
      COMMON /PARA1/ NFF, MSP, NDEP, NOPT, ICF, NDP, KSRD, MODEN
      COMMON /PARA2/ NPOINT, MAXMSH
      COMMON /XREFL/ CINTFC, RINTFC
      common /print/iprint
      real xnul
      real*8 dnul
      data xnul/0.e0/,dnul/0.d0/

  200 FORMAT(1X,'MODE=',I4, 1X,D17.10,1X,E10.4,6(1X,E9.3),
     &2X,2(F5.2,1X),F5.2)
  300 FORMAT(1X,' *** WARNING FOR MODE NO.',I4,' :',/,
     & '      INVERSE ITERATION RESTARTED WITH OFFSET EIGENVALUE.')
  400 FORMAT(1X,' *** MAX NUMBER OF ITERATIONS EXCEEDED IN EIGVEC',
     &' ***',/,' *** MAX DIFFERENCE: ',D15.8)
C
      KM=SQRT(EIGVAL)/DH0I
      NTOT=MH0I+MSEDI
      MN=MEIG+MINMOD-1
C
C
C     SET UP OF COEFFICIENT MATRIX
C
C     WATER LAYERS
C
 1000 CONTINUE

      A3(1)=ADA(1)-EIGVAL
      B3(1)=1.0D0
      DO 1200   N=2,MH0I-1
      A3(N)=ADA(N)-EIGVAL
      B3(N)=1.0D0
      C3(N-1)=1.0D0
 1200 CONTINUE
C
C     SEDIMENT LAYERS
C
      IF (MSEDI.GT.0) THEN
       A3(MH0I)=ADA(MH0I)-EIGVAL + (DH0I/(DSEDI*ROS))*
     &      ( (((EIGREF*DSEDI)*slow(MH0I+2))**2-2.) - EIGVAL*CON3 )
       B3(MH0I)=2.0D0*(DH0I/DSEDI)
       C3(MH0I-1)=2.0D0
       A3(MH0I+1)=ADA(MH0I+1)-EIGVAL*CON3
       B3(MH0I+1)=1.0D0
       C3(MH0I)=1.0D0/ROS
       DO 1400   N=MH0I+2,NTOT-1
       A3(N)=ADA(N)-EIGVAL*CON3
       B3(N)=1.0D0
       C3(N-1)=1.0D0      
 1400  CONTINUE
      END IF
      IF (EIGMIN.GT.EIGVAL) THEN
       A3(NTOT)=ADA(NTOT)-EIGVAL*CON3
      ELSE
       A3(NTOT)=ADA(NTOT)-EIGVAL*CON3-2.0D0*SQRT(EIGVAL-EIGMIN)*CON5
      END IF
      C3(NTOT-1)=2.0D0
      B3(NTOT)=dnul
      C3(NTOT)=dnul
C
C     ELIMINATION IN COEFFICIENT MATRIX
C
      DO 1600   I=1,NTOT-1
        I1=I+1
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
          C3(I) =DNUL
         endif
 1600 CONTINUE
      IF(A3(NTOT).EQ.dnul)   THEN
CF8    WRITE(08,300) MEIG
       EIGVAL=EIGVAL*(1.0D0+1.0D-12)
       GO TO 1000
      END IF
      A3(NTOT)=1.0D0/A3(NTOT)
C
C     ELIMINATION OF RIGHT HAND SIDE
C
      DO 2000   I=1,NTOT
      ZZ(I)=1.0D0
c      SSOLD(I)=1.0D0
 2000 CONTINUE
      NIT=0
 2200 IF (NIT.GT.0) THEN
        DO 2400   I=1,NTOT-1
          IF (EXCH(I)) THEN
            Q=ZZ(I+1)
            ZZ(I+1)=ZZ(I)
            ZZ(I)=Q
          END IF
          ZZ(I+1)=ZZ(I+1)+EE(I)*ZZ(I)
 2400   CONTINUE
      END IF
C
C     BACK SUBSTITUTION
C
      ZZMAX=dnul
      EIGF(NTOT)=ZZ(NTOT)*A3(NTOT)
      EIGF(NTOT-1)=(ZZ(NTOT-1)-B3(NTOT-1)*EIGF(NTOT))*A3(NTOT-1)
      ABS1=ABS(EIGF(NTOT))
      ABS2=ABS(EIGF(NTOT-1))
      ZZMAX=MAX(ABS1,ABS2)
C     ZZMAX=AMAX1(ABS(EIGF(NTOT)),ABS(EIGF(NTOT-1)))
      DO 3000   I=NTOT-2,1,-1
        EIGF(I)=(ZZ(I)-B3(I)*EIGF(I+1)-C3(I)*EIGF(I+2))*A3(I)
        ZZMAX=MAX(ZZMAX,ABS(EIGF(I)))
 3000 CONTINUE
C
C     SCALE AND CHECK ACCURACY
C
C     ZZMAX=SIGN(ZZMAX,EIGF(NTOT))
      NIT=NIT+1
      Q=dnul
      ZZINV=1.0D0/ZZMAX
      DO 4000   I=1,NTOT
        EIGF(I)=EIGF(I)*ZZINV
        Q=MAX(Q,ABS(ABS(EIGF(I))-SSOLD(I)))
 4000 CONTINUE
      IF (Q.LT.1.0D-10) GO TO 4400
      IF (NIT.GT.10) THEN 
        WRITE(luprt,400) Q
        GO TO 4400
      END IF
      DO 4200   I=1,NTOT
        ZZ(I)=EIGF(I)
        SSOLD(I)=abs(EIGF(I))
 4200 CONTINUE
      GO TO 2200
 4400 CONTINUE

C
C     NORMALIZATION
C
      ATTW=0.D0
      NORMW=0.D0
      DO 5000   N=1,MH0I-1,2
        ATTW=ATTW + 2.0*EIGF(N)**2*slow(N+1) +
     &            EIGF(N+1)**2*slow(N+2)
        NORMW=NORMW + 2.0*EIGF(N)**2 + EIGF(N+1)**2
 5000 CONTINUE
      ATTW=2*ATTW - EIGF(MH0I)**2*slow(MH0I+1)
      NORMW=2*NORMW - EIGF(MH0I)**2
      ATTW=ATTW/3.0*DH0I
      NORMW=NORMW/3.0*DH0I
      IF(MSEDI.GT.0)   THEN
        ATTS=EIGF(MH0I)**2/(ROS*ROS*slow(MH0I+2))
        NORMS=EIGF(MH0I)**2/(ROS*ROS)
        DO 5200   N=MH0I+1,NTOT-1,2
          ATTS=ATTS + 4.0*EIGF(N)**2*slow(N+2) +
     &                2.0*EIGF(N+1)**2*Slow(N+3)
          NORMS=NORMS + 4.0*EIGF(N)**2 + 2.0*EIGF(N+1)**2
 5200   CONTINUE
        ATTS=ATTS - EIGF(NTOT)**2*slow(NTOT+2)
        NORMS=NORMS - EIGF(NTOT)**2
        ATTS=ATTS/3.0*DSEDI*ROS
        NORMS=NORMS/3.0*DSEDI*ROS
      ELSE
        ATTS=dnul
        NORMS=dnul
      END IF
      NORMB=DH0I*ROB*((ROS/ROB)*EIGF(NTOT))**2*.5/SQRT(EIGVAL-EIGMIN)
      ATTB=(NORMB*CREF)/DBLE(C2)
      NORM=NORMW+NORMS+NORMB
      SQNORM=1.0/SQRT(H0*NORM)
C
c      OMEGAC=OMEGA/CREF
      zzinv=ekm/omega*cref
      DO 5400   N=1,NTOT
c       IF(OMEGAC*slow(N).LE.EKM)   GO TO 5400
        IF(slow(N).LE.zzinv)   GO TO 5400
        IF(MOD(MEIG+MINMOD-1,2).EQ.0)   THEN
C ODD ORDER MODES
          SQNORM=-SIGN(SQNORM,EIGF(N))
        ELSE
C EVEN ORDER MODES
          SQNORM=SIGN(SQNORM,EIGF(N))
        END IF
      GO TO 5600
 5400 CONTINUE
 5600 CONTINUE
C 
c the eigenfunctions are sampled at all points in the xts array
c  pg 930322
      DO 5800   N=1,NTOT
        EIGF(N)=EIGF(N)*SQNORM 
 5800 CONTINUE

      XTS(MMM,1)=dnul
      
      MODAVR(MEIG)=SQRT(NORMW)*ABS(SQNORM)

c     SAMPLING THE EIGENFUNCTIONS
C
      DEL=DFLOAT(MH0I)/DFLOAT(MSP-1)
      DO 6000   I=2,MSP
      PZ=DEL*DFLOAT(I-1)
      N1=INT(PZ)
      N2=N1+1
      IF(N1.EQ.0) THEN
        XTS(MMM,I)=(PZ-DFLOAT(N1))*EIGF(N2)
      ELSE
        XTS(MMM,I)=(PZ-DFLOAT(N1))*(EIGF(N2)-EIGF(N1))+EIGF(N1)
      END IF
 6000 CONTINUE
      XTS(MMM,MSP)=EIGF(MH0I)


C     MODE ATTENUATION

      ATTW=ATTW/NORM*EIGREF/KM
      ATTS=ATTS/NORM*EIGREF/KM
      ATTB=ATTB/NORM*EIGREF/KM

c
c--- added by pg 7/1 94
c
      if (beta(0).eq.0) then
      F2= FRQ*FRQ*1.0D-6
      BETAm1= F2 *
     &          (0.11 / (1.0 + F2) + 44.0 / (4100.0 + F2) + 3.0D-4)
       BETAm0= BETAm1 * 1.0D-3 * CC0/FRQ
       ALF0= BETAm0 * FRQ * ATTW / (CC0 * 8.68589)
c       write(*,*)' ALF0,BETAm0, F2, BETAm1', ALF0,BETAm0, F2, BETAm1
      else
        ALF0= BETA(0) * FRQ * ATTW / (CC0 * 8.68589)
      endif
      IF(MSEDI.GT.0)   THEN
      ALF1=BETA(1)*FRQ*ATTS/(CC1*8.68589)
      ELSE
      ALF1=0.0
      END IF
      ALF2=BETA(2)*FRQ*ATTB/(C2*8.68589)
C
C     BOTTOM LOSSES
C
      IF(CBS.GT.0.) THEN
C
C     REFLECTION COEFFICIENT
C
      IF (MSEDI.GT.0) THEN
      CIN=1.0d0/Slow(NTOT+2)
      ELSE
      CIN=1.0d0/Slow(MH0I+1)
      END IF
      COST=KM/EIGREF*CIN
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
      K1=ABS((EIGREF/CIN)**2-KM*KM)
      SQRK1=SQRT(K1)
      K2=KM*KM-((OMEGA*H0)/DBLE(C2))**2
      ALF2S=RINTFC*EIGF(NTOT)**2*SQRK1/(8.0*KM)
      ALF2S=ALF2S*(1.0+(RINTFC/R2)**2*K2/K1)*Q
      ALFA=ALF0+ALF1+ALF2S
      ELSE
      ALF2S=0.0
      ALFA=ALF0+ALF1+ALF2
      END IF
C
C     SCATTER LOSSES
C
C
      DER0=EIGF(1)/(H0*DH0I)
      IF (MSEDI.GT.0) THEN
      DER1=(CON1*EIGF(MH0I)-(CON1+CON4)*EIGF(MH0I-1)+
     & CON4*ROS*EIGF(MH0I+1))/(2.0*DH0I*H0)
      ELSE
      K1=KM*KM-((OMEGA*H0)/DBLE(C2))**2
      SQK1=SQRT(K1)
      DER1=-EIGF(MH0I)*SQK1/(ROB*H0)
      END IF
      SQKM=KM*KM
      RATS=(EIGREF*slow(1))**2
      RATB=(EIGREF*slow(MH0I+1))**2
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
      ALFOB=SCATT(2)**2*KH0/(2.*KM)*(DER1**2+(EIGF(MH0I)*KH0/H0)**2)
      ALFA=ALFA+ALFOB+ALFOS
C
      CHECKS=2.0*((RATS-SQKM)/(H0*H0))*SCATT(1)**2
      CHECKB=2.0*((RATB-SQKM)/(H0*H0))*SCATT(2)**2
C
      IF(FLAGPU .LT.1.0)   WRITE(LUPRT,200)MN,EKM,ALFA,ALF0,
     & ALF1,ALF2,ALF2S,ALFOS,ALFOB,Q,CHECKS,CHECKB
C
C
C
      RETURN
      END
