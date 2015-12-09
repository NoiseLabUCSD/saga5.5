C**********************************************************
C
      SUBROUTINE STURM(EIGVAL,MH0I,MSEDI,DH0I,DSEDI,S2,M,ADA,
     & NPOINT)
C
C__________________________________________________________
C                                                         |
C     This routine counts the number of sign changes in   |
C     the Sturm sequence and calculates the value of the  |
C     characteristic equation.                            |
C_________________________________________________________|  
C

      DOUBLE PRECISION H0, H1, ROB, ROS
      DOUBLE PRECISION TWOPI, PI, OMEGA
      DOUBLE PRECISION DH0I, DSEDI, U1H0, U1SED, SEDK
      DOUBLE PRECISION ADA(NPOINT)
      DOUBLE PRECISION EIGVAL, S0, S1, S2, STIFF
      DOUBLE PRECISION EIGREF,EIGMIN,EIGMAX,EIGX
      DOUBLE PRECISION CON1, CON2, CON3, CON4, CON5

      COMMON /CONST/ CON1, CON2, CON3, CON4, CON5, SEDK
      COMMON /DENS/ R0, R1, R2
      COMMON /DENS8/ ROB, ROS
      COMMON /EXPMAX/ TRESH, EPS, RRMAX
      COMMON /GSNAP/ H0, H1, TWOPI, PI, OMEGA
      COMMON /GEN/ EIGREF, EIGMIN, EIGMAX, STIFF

      DATA SCFAC/1.0E-20/
C
C     'M' COUNTS SIGN CHANGES
C
C     WATER LAYERS
C
      EIGX=EIGVAL
      S0=1.D0
      S1=(EIGX-ADA(1))*S0
      M=0
c      write(*,*)'tresh,s0,s1',tresh,s0,s1
      DO 1000   N=2,MH0I-1
      S2=(EIGX-ADA(N))*S1-S0
c      write(*,*)'s2:',n,s2
      IF( ((S2 .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S2 .LT. 0.0))     )   M=M+1
       S0=S1
       S1=S2
      IF( ABS(S2) .GE. TRESH )   THEN
        S0=S0*SCFAC
        S1=S1*SCFAC
      END IF
c      write(*,*)'Stum N,s2,m 1000:',n,s2,m
 1000 CONTINUE
c      write(*,*)'Stum s2,m 1000:',s2,m
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
C
C      SEDK=((OMEGA*H1)/(C11*DFLOAT(MSEDI)))**2
C      S2=(U1SED*2*DSEDI+2*S1-(SEDK-EIGX)*S1)/2.0
      S2=U1SED*DSEDI - 0.5*(SEDK-EIGX-2.0D0)*S1
C
      IF( ((S2 .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S2 .LT. 0.0))     )   M=M+1
       S0=S1
       S1=S2
      IF( ABS(S2) .GE. TRESH )   THEN
        S0=S0*SCFAC
        S1=S1*SCFAC
      END IF
c      IF( ABS(S2) .GE. TRESH )   THEN 
c       S0=S1*SCFAC
c       S1=S2*SCFAC
c      ELSE 
c       S0=S1
c       S1=S2
c      END IF
C
      DO 2000   N=MH0I+1,MH0I+MSEDI-1
      S2=(EIGX-ADA(N))*S1-S0
      IF( ((S2 .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S2 .LT. 0.0))     )   M=M+1
       S0=S1
       S1=S2
      IF( ABS(S2) .GE. TRESH )   THEN
        S0=S0*SCFAC
        S1=S1*SCFAC
      END IF
c      IF( ABS(S2) .GE. TRESH )   THEN 
c       S0=S1*SCFAC
c       S1=S2*SCFAC
c      ELSE 
c       S0=S1
c       S1=S2
c      END IF
 2000 CONTINUE
c      write(*,*)'Stum s2,m 2000:',s2,m
C
      END IF
C
C     BOTTOM
C
      IF(EIGMIN.GT.EIGVAL) THEN
        S2=(EIGX-ADA(MH0I+MSEDI))*S1-2.D0*S0
      ELSE
        S2=(EIGX-ADA(MH0I+MSEDI) +
     &     2.D0*SQRT(EIGVAL-EIGMIN)*CON5)*S1-2.D0*S0
      ENDIF
      IF( ((S2 .GT.0.0) .AND. (S1 .LT. 0.0)) .OR.
     &    ((S1 .GT.0.0) .AND. (S2 .LT. 0.0))     )   M=M+1
c      write(*,*)'Stum s2,m bot :',s2,m
      RETURN
      END 









