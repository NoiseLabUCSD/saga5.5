C**********************************************************
C
      SUBROUTINE BRENT(A,B,X,MH0I,MSEDI,DH0I,DSEDI,ADA,
     & NPOINT)
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
 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      CHARACTER*80 ERRMSG

      DOUBLE PRECISION DH0I, DSEDI
      DOUBLE PRECISION CON1, CON2, CON3, CON4, CON5, SEDK
      DOUBLE PRECISION MACHEP, M, T, TREF
      DOUBLE PRECISION ADA(NPOINT)

      COMMON /CONST/ CON1, CON2, CON3, CON4, CON5, SEDK
      COMMON /LUNIT/ LUPLP, LUPLT, LUPRT

      DATA TREF/1.0D-13/


      ERRMSG = ' '
      MACHEP = 1.0E-16
      TEN = 10.0
      T=TREF*A
C      CALL FUNCT(A,FA,IEXPA)
C      CALL FUNCT(B,FB,IEXPB)
      CALL CHARAC(A,MH0I,MSEDI,DH0I,DSEDI,FA,NOVFLA,ADA,
     & NPOINT)
      CALL CHARAC(B,MH0I,MSEDI,DH0I,DSEDI,FB,NOVFLB,ADA,
     & NPOINT)
c       WRITE(*,*) ' SUB BRENT : ',A,FA,B,FB
c        WRITE(*,*) MH0I, ADA(1), ADA(2),ADA(3)
      IEXPA=NOVFLA*20
      IEXPB=NOVFLB*20
      IF ( ( (FA .GT. 0.0) .AND. (FB .GT. 0.0) ) .OR.
     &     ( (FA .LT. 0.0) .AND. (FB .LT. 0.0) ) ) THEN
         ERRMSG = ' *** ZBRENT ERROR: FUNCT SGN SAME AT INTRVL ENDPTS'
         WRITE(6,*) ERRMSG,'FA,FB=',fa,fb
         IF(LUPRT .NE. 6)   WRITE(LUPRT,*) ERRMSG
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
      CALL CHARAC(B,MH0I,MSEDI,DH0I,DSEDI,FB,NOVFLB,ADA,
     & NPOINT)
      IEXPB=NOVFLB*20
C         CALL FUNCT(B,FB,IEXPB)
         IF ((FB .GT. 0.0) .EQV. (FC .GT. 0.0)) GOTO 2000
         GOTO 3000
      ENDIF
      X = B
      RETURN
      END
