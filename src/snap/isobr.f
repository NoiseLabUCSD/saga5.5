      SUBROUTINE ISOBR(I,H1N,DRAT,FRQ,MINMOD,MAXMOD,MODQTY,
     & DH0,DSED,MH0,MSED,
     & CREF,ADA,SPEED,ISO,NBEG,MY,C0,Z0,C1,Z1,FIRST)

C__________________________________________________________
C                                                          |
C     The subject is to find an approximate solution for   |
C     the eigenvalues of a continuous diff. equation by    |
C     finite difference and extrapolation                  |
C__________________________________________________________|

      INTEGER MH0(8),MSED(8)

      DOUBLE PRECISION CREF, CMIN, H0, H1, H1N, ROB, ROS
      DOUBLE PRECISION MY(MAXMSH, MODEN), MYBR, ISO(MODEN),
     &                 ADA(NPOINT)
      DOUBLE PRECISION DH0(8), DSED(8)
      DOUBLE PRECISION SPEED(NPOINT)
      DOUBLE PRECISION FRQ, A, B, STIFF, DRAT
      DOUBLE PRECISION TWOPI, PI, OMEGA
      DOUBLE PRECISION EIGREF, EIGMIN, EIGMAX
      DOUBLE PRECISION HRAT, HRATSQ
      DOUBLE PRECISION CON1, CON2, CON3, CON4, CON5, SEDK
      DOUBLE PRECISION C0(NDEP), Z0(NDEP), C1(NDEP), Z1(NDEP)

      COMMON /ATTEN/ ALF0, ALF1, ALF2, ALF2S, ALFOS, ALFOB
      COMMON /CONST/ CON1, CON2, CON3, CON4, CON5, SEDK
      COMMON /DENS/ R0, R1, R2
      COMMON /DENS8/ ROB, ROS
      COMMON /GSNAP/ H0, H1, TWOPI, PI, OMEGA
      COMMON /GEN/ EIGREF, EIGMIN, EIGMAX, STIFF
      COMMON /NA/ ND0, ND1, CMIN
      COMMON /PARA1/ NFF, MSP, NDEP, NOPT, ICF, NDP, KSRD, MODEN
      COMMON /PARA2/ NPOINT, MAXMSH
C
C     INITIALIZATION
C
      IF(H1N.LT.1.0D-6) THEN
       HRAT=1.
       HRATSQ=1.
       CON3=1.
      ELSE
       HRAT=(MSED(I)*H0)/(MH0(I)*H1)
       HRATSQ=HRAT**2
       CON3=(MH0(I)*H1)**2/(MSED(I)*H0)**2
      END IF
      CON1=2.*(STIFF-HRATSQ/ROS)/(STIFF+HRAT)
      CON4=2./ROS*HRATSQ/(STIFF+HRAT)
      CON5=ROS/(ROB*HRAT)
c      write(*,*)' isobr 141'
      CALL ISOINT(ISO,MAXMOD,MINMOD,DH0(I),MH0(I),MSED(I),DSED(I),
     & CREF,FRQ,SPEED,ADA,C0,Z0,C1,Z1,FIRST)
      MODQTY=MAXMOD-MINMOD+1
      IF(MODQTY .LE. 0)   RETURN


      DO 1000 M=1,MODQTY
      A=ISO(M+1)
      B=ISO(M)
      CALL BRENT(A,B,MYBR,MH0(I),MSED(I),DH0(I),DSED(I),ADA,
     & NPOINT)
      MY(I,M)=MYBR
 1000 CONTINUE
      RETURN
      END
C
C
C***********************************************************
C





