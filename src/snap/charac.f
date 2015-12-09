      SUBROUTINE CHARAC(EIGVAL,MH0I,MSEDI,DH0I,DSEDI,S2,NOVFL,
     & ADA,NPOINT)
C
C_____________________________________________________________
C                                                             |
C     This routine calculates the value of the characteristic |
C     equation of the matrix.                                 |
C_____________________________________________________________|
C
      DOUBLE PRECISION H0, H1, ROB, ROS
      DOUBLE PRECISION TWOPI, PI, OMEGA
      DOUBLE PRECISION DH0I, DSEDI, U1H0, U1SED, SEDK
      DOUBLE PRECISION EIGVAL, S0, S1, S2,s3, STIFF,ztemp
      DOUBLE PRECISION ADA(NPOINT)
      DOUBLE PRECISION EIGREF, EIGMIN, EIGMAX, EIGX
      DOUBLE PRECISION CON1, CON2, CON3, CON4, CON5

      COMMON /CONST/ CON1, CON2, CON3, CON4, CON5, SEDK
      COMMON /DENS/ R0, R1, R2
      COMMON /DENS8/ ROB, ROS
      COMMON /EXPMAX/ TRESH, EPS, RRMAX
      COMMON /GSNAP/ H0, H1, TWOPI, PI, OMEGA
      COMMON /GEN/ EIGREF, EIGMIN, EIGMAX, STIFF

      DATA SCFAC/1.0E-20/ 
C
C     WATER LAYERS
C
      EIGX=EIGVAL
      NOVFL=0
      S0=1.D0
      S1=(EIGX-ADA(1))
C
      DO 1000   N=2,MH0I-1,2
        S2=(EIGVAL-ADA(N))*S1-S0
        S3=(EIGVAL-ADA(N+1))*S2-S1
C
        S0=S2
        S1=S3
        IF( ABS(S2) .GE. TRESH )   THEN 
          NOVFL=NOVFL + 1
          S0=S0*SCFAC
          S1=S1*SCFAC
        END IF
C
 1000 CONTINUE
C
C     SEDIMENT LAYERS
C
      IF(MSEDI.GT.0) THEN
C
C   FICTITIOUS POINT FROM WATER TO SEDIMENT 
       S2=(EIGX-ADA(MH0I))*S1 - S0
C   DERIVATIVE AT S1
       U1H0=(S2-S0)/(2.0*DH0I)
       U1SED=U1H0
C
C  INTERFACE POINT AS SEEN FROM THE SEDIMENT
       S1=S1/ROS
       EIGX=EIGVAL*CON3
C
      S2=U1SED*DSEDI - 0.5*(SEDK-EIGX-2.0D0)*S1
C
        S0=S1
        S1=S2
        IF( ABS(S2) .GE. TRESH )   THEN 
          NOVFL=NOVFL + 1
          S0=S0*SCFAC
          S1=S1*SCFAC
        END IF
C
C
       DO 2000   N=MH0I+1,MH0I+MSEDI-1
        S2=(EIGX-ADA(N))*S1-S0
        S0=S1
        S1=S2
        IF( ABS(S2) .GE. TRESH )   THEN 
          NOVFL=NOVFL + 1
          S0=S0*SCFAC
          S1=S1*SCFAC
        END IF
 2000  CONTINUE
      ENDIF
C
C     BOTTOM
C
      IF(EIGMIN.GT.EIGVAL) THEN
        S2=(EIGX-ADA(MH0I+MSEDI))*S1-2.D0*S0
      ELSE
        S2=(EIGX-ADA(MH0I+MSEDI) +
     &     2.D0*SQRT(EIGVAL-EIGMIN)*CON5)*S1-2.D0*S0
      ENDIF
      RETURN
      END
C
C
C**********************************************************
C
