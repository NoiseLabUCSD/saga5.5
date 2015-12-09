      SUBROUTINE SCAIRY (ZS,AIS,BIS,AIPS,BIPS,ZETAS)
      IMPLICIT COMPLEX*16 (A,B,Z)
      REAL*8 ZPR(2)
      EQUIVALENCE (ZPR(1),ZPP)
      COMPLEX*8 ZS,AIS,BIS,AIPS,BIPS,ZETAS
      COMPLEX*16 PHSE13,PHSE23,PHSE16,PHSE56
      COMPLEX*16 PHSC13,PHSC23,PHSC16,PHSC56
      DATA     PHSE13/( 0.5D0, 0.8660254037844387D0)/,
     1         PHSE23/(-0.5D0, 0.8660254037844387D0)/,
     2         PHSE16/( 0.8660254037844387D0, 0.5D0)/,
     3         PHSE56/(-0.8660254037844387D0, 0.5D0)/,
     4         PHSC13/( 0.5D0,-0.8660254037844387D0)/,
     5         PHSC23/(-0.5D0,-0.8660254037844387D0)/,
     6         PHSC16/( 0.8660254037844387D0,-0.5D0)/,
     7         PHSC56/(-0.8660254037844387D0,-0.5D0)/
      IRET=1
      Z=ZS
      ZPP=Z
      CALL CLAIRY(Z,IRET,AI,BI,AIP,BIP,ZETA)
      IF (ZPR(2).GE.0.D0)  THEN
        BI  = 2.D0*PHSC16  * BI
        BIP = 2.D0*PHSC56  * BIP
      ELSE
        BI  = 2.D0*PHSE16  * BI
        BIP = 2.D0*PHSE56  * BIP
      END IF
      AIS=AI
      BIS=BI
      AIPS=AIP
      BIPS=BIP
      ZETAS=ZETA
      RETURN
      END
      SUBROUTINE CLAIRY (Z,IRET, AI,BI,AIP,BIP,ZETA)
      IMPLICIT COMPLEX*16 (A,B,E,X,Z)
      PARAMETER (ONEOOF=1D0/1.5D0)

C     COMPLEX AIRY FUNCTION SUBROUTINE
C     Z; THE COMPLEX ARGUMENT
C     IRET; RETURN CONTROL PARAMETER FOR SCALING
C           < 0 RETURN SCALED AIRY FUNCTIONS OF THE FIRST AND SECOND TYPES.
C           = 0 RETURN UNSCALED AIRY FUNCTIONS OF THE FIRST AND SECOND TYPES.
C           > 0 RETURN SCALED RESULTS OF THE NUMERICALLY STABLE INDEPENDENT 
C               TYPE.  THESE ARE AI(Z) AND AI(Z*EXP([+,-]2IPI/3)).
C     AI AND AIP; THE COMPLEX AIRY FUNCTION OF THE FIRST TYPE AND ITS
C                   DERIVATIVE DEPENDING ON IRET.
C               IRET < 0, SCALED AIRY FUNCTION OF THE FIRST TYPE.  AI AND AIP
C                   ARE TO BE DIVIDED BY EXP(ZETA) TO OBTAIN THE TRUE AIRY 
C                   FUNCTIONS.
C               IRET = 0, UNSCALED AIRY FUNCTION OF THE FIRST TYPE
C               IRET > 0, SCALED AIRY FUNCTION OF THE FIRST TYPE.  AI AND AIP
C                   ARE TO BE DIVIDED BY EXP(ZETA) TO OBTAIN THE TRUE AIRY 
C                   FUNCTIONS.
C     BI AND BIP; THE COMPLEX AIRY FUNCTION OF THE SECOND TYPE AND ITS
C               DERIVATIVE DEPENDING ON IRET AND THE ANGLE OF Z.
C               IRET < 0, SCALED AIRY FUNCTION OF THE SECOND TYPE.  BI AND BIP 
C                   ARE TO BE MULTIPLIED BY EXP(ZETA), FOR ARG(Z) < 60 AND 
C                   DIVIDED BY EXP(ZETA) FOR ARG(Z) > 60.
C               IRET = 0, UNSCALED AIRY FUNCTION OF THE SECOND TYPE
C               IRET > 0, SCALED AIRY FUNCTION OF OF THE FIRST TYPE FOR THE
C     ARGUMENT Z*EXP(-2IPI/3) WHEN ARG(Z) > 0 AND FOR
C     Z*EXP(2IPI/3) WHEN ARG(Z) < 0. THE RESULTS ARE TO
C     BE MULTIPLIED BY EXP(ZETA). 
C      ZETA; THE EXPONTIAL SCALING FACTOR
C
C     THIS ROUTINE USES THE ALGORITHMS DESCRIBED IN "AN ALGORITHM FOR 
C     THE EVALUATION OF THE COMPLEX AIRY FUNCTIONS,"  Z. SCHULTEN, 
C     D.G.M. ANDERSON, AND R.G. GORDEN, IN JOUR. COMP. PHYSICS 31, 
C     PP. 60-75, (1979).
C     A BETTER DESCRIPTION OF THE SERIES EXPANSION IS AVAILABLE IN
C     HANDBOOK OF MATHEMATICAL FUNCTIONS, APPLIED MATHEMATICS SERIES #55,
C     M. ABRAMOWITZ AND I. A. STEGAN, P. 446 (1964).

C      PARAMETER (TPIO3=2.094395102393195D0, I=(0.0D0,1.0D0))
      REAL*8 TPIO3
      REAL*8 ZPR(2)
      EQUIVALENCE (ZPR(1),ZPP)
      COMPLEX*16 PHSE13,PHSE23,PHSE16,PHSE56
      COMPLEX*16 PHSC13,PHSC23,PHSC16,PHSC56
      DATA     PHSE13/( 0.5D0, 0.8660254037844387D0)/,
     1         PHSE23/(-0.5D0, 0.8660254037844387D0)/,
     2         PHSE16/( 0.8660254037844387D0, 0.5D0)/,
     3         PHSE56/(-0.8660254037844387D0, 0.5D0)/,
     4         PHSC13/( 0.5D0,-0.8660254037844387D0)/,
     5         PHSC23/(-0.5D0,-0.8660254037844387D0)/,
     6         PHSC16/( 0.8660254037844387D0,-0.5D0)/,
     7         PHSC56/(-0.8660254037844387D0,-0.5D0)/

      DATA TPIO3 /2.094395102393195D0/
C   * STATEMENT FUNCTION DEFINITION OF THE PHSE *

C      THETA (Z) = DATAN2 (DIMAG (Z), REAL (Z))

C   * Z OR CONJUGATE Z IN POSITIVE IMAGINARY PLANE *

      ZPP=Z

      IF (ZPR(2) .LT. 0.0D0)  THEN
          ZP = CONJG (Z)
          ICONJG = 1
      ELSE
          ZP = Z
          ICONJG = 0
      END IF


C   *********************************************************
C   *               *
C   *   SERIES SOLUTIONS FOR AI & BI WITHIN 3.5 OF ORIGIN   *
C   *               *
C   *********************************************************

      IF (ABS (Z) .LE. 3.5D0)  THEN
          CALL ABSERIES (Z, AI,AIP,BI,BIP)
          ZETA = (0.0D0, 0.0D0)
          IF (IRET.NE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,0,IRET)


C   ******************************************************************
C   *    *
C   *   ANGLES LESS THAN 120 DEGS. - AI(Z) AND AI(Z*EXP(-2IPI/3))    *
C   *   ARE CALCULATED FOR Z IN THE UPPER HALF PLANE, CONJGATE OF    *
C   *   Z IF IN THE LOWER.  SERIES OR GAUSSIAN SOLUTION AS NEEDED.   *
C   *    *
C   ******************************************************************

      ELSE IF (ATAN2(ZPR(2),ZPR(1)) .LE. TPIO3)  THEN
          IF (ABS(ZP).GT.7.9111D0)  THEN          
              CALL ABGAUSS (ZP, AI, AIP, BI, BIP, ZETA)
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  BI  = CONJG (BI)
                  AIP = CONJG (AIP)
                  BIP = CONJG (BIP)
              END IF
          ELSE
              IF (ABS (ZP-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
                  CALL ASERIES (ZP, AI,AIP, ZETA)
              ELSE
                  CALL AGAUSS (ZP, AI, AIP, ZETA)
              END IF
              ZP2 = CONJG (ZP) * PHSE23
              IF (ABS (ZP2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
                  CALL ASERIES (ZP2, BI, BIP, ZETA2)
              ELSE
                  CALL AGAUSS (ZP2, BI, BIP, ZETA2)
              END IF
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  AIP = CONJG (AIP)
              ELSE
                  BI  = CONJG (BI)
                  BIP = CONJG (BIP)
              END IF
          END IF
          IF (IRET.LE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,1,IRET)


C   **********************************************************
C   *                *
C   *   CONNECTION REGION FOR ANGLE GREATER THAN 120 DEGS.   *
C   *   CALCULATE AI(Z*EXP(-2IPI/3)) AND AI(Z*EXP(2IPI/3))   *
C   *   SOLUTIONS, USE CONNECTION FORMULAS FOR RESULTS.      *
C   *                *
C   **********************************************************

      ELSE 
          Z1 = Z * PHSC23    
          Z2 = CONJG (Z) * PHSC23                
          
          IF (ABS (Z1-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES (Z1, AI1, AIP1, ZETA1)
          ELSE
              CALL AGAUSS (Z1, AI1, AIP1, ZETA1)
          END IF
          IF (ABS (Z2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES (Z2, AI2, AIP2, ZETA2)
          ELSE
              CALL AGAUSS (Z2, AI2, AIP2, ZETA2)
          END IF
          AI2  = CONJG (AI2)  
          AIP2 = CONJG (AIP2)
          ZETA2 = CONJG (ZETA2)

          IF (IRET.EQ.0)  THEN
              E1 = EXP (-ZETA1)
              E2 = EXP (-ZETA2)
              AI  = PHSE13 * AI1*E1  + PHSC13 * AI2*E2
              AIP = PHSC13 * AIP1*E1 + PHSE13 * AIP2*E2
              BI  = PHSE16 * AI2*E2  + PHSC16 * AI1*E1
              BIP = PHSE56 * AIP2*E2 + PHSC56 * AIP1*E1
              ZETA = (0.0D0, 0.0D0)

          ELSE IF (IRET.LT.0)  THEN               
              ZETA  = Z * SQRT(Z) * ONEOOF
              E1  = EXP (ZETA-ZETA1)
              E2  = EXP (ZETA-ZETA2)
              AI  = PHSE13 * AI1*E1  + PHSC13 * AI2*E2
              AIP = PHSC13 * AIP1*E1 + PHSE13 * AIP2*E2
              BI  = PHSC16 * AI1*E1  + PHSE16 * AI2*E2
              BIP = PHSC56 * AIP1*E1 + PHSE56 * AIP2*E2

          ELSE IF (IRET.GT.0)  THEN               
              ZETA  = Z * SQRT(Z) * ONEOOF
              E1  = EXP (ZETA-ZETA1)
              E2  = EXP (ZETA-ZETA2)
              AI  = PHSE13 * AI1*E1  + PHSC13 * AI2*E2
              AIP = PHSC13 * AIP1*E1 + PHSE13 * AIP2*E2
              IF (ZPR(2).GE.0.0D0)  THEN
                  E1  = EXP (-ZETA-ZETA1)
                  BI  = AI1  * E1
                  BIP = AIP1 * E1
              ELSE
                  E2  = EXP (-ZETA-ZETA2)
                  BI  = AI2  * E2
                  BIP = AIP2 * E2
              END IF
          END IF
      END IF

      RETURN
      END
      SUBROUTINE ABSERIES (Z, AI,AIP,BI,BIP)
      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      REAL*8 C1,C2,SQR3,FAC1,FAC2,FAC3
      DATA C1 /0.355028053887817D0/, C2 /0.258819403792807D0/,
     1          SQR3 /1.732050807568877D0/

      FC = (1.0D0, 0.0D0)
      GC = Z
      Z2 = Z * Z
      AI = C1*FC - C2*GC
      BI = C1*FC + C2*GC
      AIP = -C2
      BIP =  C2

      DO 10 K = 3,600,3
      FAC1=1.0D0/(K-1)
      FAC2=1.0D0/K
      FAC3=1.0D0/(K+1)
      FC = FC * Z2 * FAC1    
      GC = GC * Z2 * FAC2
      AIP = AIP + (C1*FC - C2*GC)
      BIP = BIP + (C1*FC + C2*GC)
      FC = FC * Z * FAC2       
      GC = GC * Z * FAC3
      AI = AI + (C1*FC - C2*GC)
      BI = BI + (C1*FC + C2*GC)
      IF (ABS(FC).LT.1.D-17*ABS(AI) .AND. 
     &    ABS(GC).LT.1.D-17*ABS(AI))  GO TO 12
   10 CONTINUE
      WRITE (6,*) ' SERIES APPROXIMATION FAILED TO CONVERGE'

   12 BI = SQR3 * BI
      BIP = SQR3 * BIP
      RETURN
      END
      SUBROUTINE ABGAUSS (Z, AI,AIP,BI,BIP,ZETA)

      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      PARAMETER (ONEOOF=1D0/1.5D0)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (ONEOPI=1.0D0/PI) 
      COMPLEX*16 SQRTZ,SQRPIZ,CFAC,CFAC2,ONEOZ
      REAL*8 ZPR(2)
      EQUIVALENCE (ZPR(1),ZPP)
      COMPLEX*16 EYE
      REAL*8 X4(4), W4(4), X6(6), W6(6)
      COMPLEX*16 PHSE13,PHSE23,PHSE16,PHSE56
      COMPLEX*16 PHSC13,PHSC23,PHSC16,PHSC56
      DATA     PHSE13/( 0.5D0, 0.8660254037844387D0)/,
     1         PHSE23/(-0.5D0, 0.8660254037844387D0)/,
     2         PHSE16/( 0.8660254037844387D0, 0.5D0)/,
     3         PHSE56/(-0.8660254037844387D0, 0.5D0)/,
     4         PHSC13/( 0.5D0,-0.8660254037844387D0)/,
     5         PHSC23/(-0.5D0,-0.8660254037844387D0)/,
     6         PHSC16/( 0.8660254037844387D0,-0.5D0)/,
     7         PHSC56/(-0.8660254037844387D0,-0.5D0)/

      DATA EYE/(0.0D0,1.0D0)/

      DATA X4/3.9198329554455091D0,  1.6915619004823504D0,
     1        5.0275532467263918D-1, 1.9247060562015692D-2/,
     2     W4/4.7763903057577263D-5, 4.9914306432910959D-3,
     3        8.6169846993840312D-2, 9.0879095845981102D-1/,
     4     X6/7.1620871339075440D0,  4.2311006706187214D0,
     5        2.3361772245064852D0,  1.0856431202004936D0,
     6        3.3391648924379639D-1, 1.3115888501576988D-2/,
     7     W6/4.9954496303045166D-8, 1.8066384626280827D-5,
     8        9.5530673977919037D-4, 1.5715675321710695D-2,
     9        1.1588902608004444D-1, 8.6742187551934309D-1/

      ZPP=Z

      SQRTZ = SQRT (Z)
      ZETA   = Z * SQRTZ * ONEOOF
      SQRPIZ = SQRT (PI * SQRTZ)
      AI = (0.0D0, 0.0D0)
      BI = (0.0D0, 0.0D0)
      AIP = (0.0D0, 0.0D0)
      BIP = (0.0D0, 0.0D0)

      IF (ABS(Z).GT.25.D0)  THEN
          DO 10 I = 1,4
          CFAC=1.0D0/(ZETA+X4(I))
          CFAC2=1.0D0/(ZETA-X4(I))
          AI = AI + W4(I) * ZETA * CFAC
          AIP = AIP + W4(I) * X4(I) * SQRTZ * CFAC**2
          BIP = BIP + W4(I) * X4(I) * SQRTZ * CFAC2**2
   10     BI = BI + W4(I) * ZETA * CFAC2
      ELSE
          DO 20 I = 1,6
          CFAC=1.0D0/(ZETA+X6(I))
          CFAC2=1.0D0/(ZETA-X6(I))
          AI = AI + W6(I) * ZETA * CFAC
          AIP = AIP + W6(I) * X6(I) * SQRTZ * CFAC**2
          BIP = BIP + W6(I) * X6(I) * SQRTZ * CFAC2**2
   20     BI = BI + W6(I) * ZETA * CFAC2
      END IF
      CFAC=1.0D0/SQRPIZ
      ONEOZ=1.0D0/Z
      AIP = - AI *0.125D0 * CFAC * ONEOZ
     1      - AI * SQRPIZ * 0.5D0 * ONEOPI
     1      + AIP * 0.5D0 * CFAC
      BIP = - BI * 0.25D0 * CFAC * ONEOZ
     1      + BI * SQRPIZ * ONEOPI
     1      - BIP * CFAC
      AI = AI * 0.5D0 * CFAC
      BI = BI * CFAC
      IF (ZPR(2).GE.0)  THEN                
          BI  = PHSE16 * BI * 0.5D0
          BIP = PHSE56 * BIP * 0.5D0
      ELSE
          BI  = PHSC16 * BI * 0.5D0
          BIP = PHSC56 * BI * 0.5D0
      END IF
      RETURN
      END
      SUBROUTINE ASERIES (Z, AI, AIP, ZETA)
      IMPLICIT COMPLEX*16 (A-H,O-Z) 
      PARAMETER (ONEOOF=1D0/1.5D0)
      REAL*8 C1,C2,FAC1,FAC2,FAC3
      DATA C1/0.355028053887817D0/, C2/0.258819403792807D0/

      FC = (1.0D0, 0.0D0)
      GC = Z
      Z2 = Z * Z
      AI = C1*FC - C2*GC
      AIP= -C2

      DO 10 K = 3,600,3
      FAC1=1.0D0/(K-1)
      FAC2=1.0D0/K
      FAC3=1.0D0/(K+1)
      FC = FC * Z2 * FAC1                   
      GC = GC * Z2 * FAC2
      AIP = AIP + (C1*FC - C2*GC)
      FC = FC * Z * FAC2    
      GC = GC * Z * FAC3
      AI = AI + (C1*FC - C2*GC)
      IF (ABS(FC).LT.1.D-17*ABS(AI) .AND. 
     &    ABS(GC).LT.1.D-17*ABS(AI))  GO TO 12
   10 CONTINUE
      WRITE (6,*) ' SERIES APPROXIMATION FAILED TO CONVERGE'

   12 ZETA = Z * SQRT(Z) * ONEOOF
      EXPZETA = EXP (ZETA)
      AI = AI * EXPZETA
      AIP = AIP * EXPZETA
      RETURN
      END
      SUBROUTINE AGAUSS (Z, AI,AIP,ZETA)
      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      PARAMETER (ONEOOF=1D0/1.5D0)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (ONEOPI=1.0D0/PI) 
      REAL*8 X4(4), W4(4), X6(6), W6(6)
     
      DATA X4/3.9198329554455091D0,  1.6915619004823504D0,
     1        5.0275532467263918D-1, 1.9247060562015692D-2/,
     2     W4/4.7763903057577263D-5, 4.9914306432910959D-3,
     3        8.6169846993840312D-2, 9.0879095845981102D-1/,
     4     X6/7.1620871339075440D0,  4.2311006706187214D0,
     5        2.3361772245064852D0,  1.0856431202004936D0,
     6        3.3391648924379639D-1, 1.3115888501576988D-2/,
     7     W6/4.9954496303045166D-8, 1.8066384626280827D-5,
     8        9.5530673977919037D-4, 1.5715675321710695D-2,
     9        1.1588902608004444D-1, 8.6742187551934309D-1/

      SQRTZ = SQRT (Z)
      ZETA = Z * SQRTZ * ONEOOF
      SQRPIZ = SQRT (PI * SQRTZ)
      AI  = (0.0D0, 0.0D0)
      AIP = (0.0D0, 0.0D0)

      IF (ABS(Z).GT.25.D0)  THEN
          DO 10 I = 1,4
          CFAC=1.0D0/(ZETA+X4(I))
          AI  = AI + W4(I) * ZETA * CFAC
   10     AIP = AIP + W4(I) * X4(I) * SQRTZ * CFAC**2
      ELSE
          DO 20 I = 1,6
          CFAC=1.0D0/(ZETA+X6(I))
          AI  = AI + W6(I) * ZETA * CFAC
   20     AIP = AIP + W6(I) * X6(I) * SQRTZ * CFAC**2
      END IF
      CFAC=1.0D0/SQRPIZ
      ONEOZ=1.0D0/Z
      AIP = - AI * 0.125D0 * CFAC * ONEOZ
     1      - AI * SQRPIZ * 0.5D0 * ONEOPI
     1      + AIP * 0.5D0 * CFAC
      AI  = AI * 0.5D0 * CFAC
      RETURN
      END
      SUBROUTINE FORMAIRY (AI, BI, AIP, BIP, Z, ZETA, IRETIN, IRETOUT)
      IMPLICIT COMPLEX*16 (A-E,X-Z)

C     THIS ROUTINE IS DESIGNED TO CONVERT FROM ONE TYPE OF AIRY FUNCTIONS
C     TO ANOTHER.  THE THREE REPRESENTATIONS CORRESPOND TO 
C     IRET = 0, AI AND BI UNSCALED
C     IRET < 0, AI AND BI SCALED, AI ALWAYS BY DIVIDING BY EXP(ZETA), BI BY
C               MULTIPLYING BY EXP(ZETA) FOR REAL(ZETA) > 0 (OR ARG(Z) < 60) AND
C               BY DIVIDING BY EXP(ZETA) FOR REAL(ZETA) < 0 (OR ARG(Z) > 60).
C     IRET > 0  AI(Z) AND AI(Z*EXP(2IPI/3)), THE FIRST BY DIVIDING BY EXP(ZETA)
C               AND THE SECOND BY MULTIPLYING BY EXP(ZETA).
C     FOR A TOTAL OF NINE CASES, THREE OF WHICH ARE TRIVAL.

      REAL*8 ZPR(2)
      EQUIVALENCE (ZPR(1),ZPP)
      COMPLEX*16 PHSE13,PHSE23,PHSE16,PHSE56
      COMPLEX*16 PHSC13,PHSC23,PHSC16,PHSC56
      COMPLEX*16 COMPI
      DATA     PHSE13/( 0.5D0, 0.8660254037844387D0)/,
     1         PHSE23/(-0.5D0, 0.8660254037844387D0)/,
     2         PHSE16/( 0.8660254037844387D0, 0.5D0)/,
     3         PHSE56/(-0.8660254037844387D0, 0.5D0)/,
     4         PHSC13/( 0.5D0,-0.8660254037844387D0)/,
     5         PHSC23/(-0.5D0,-0.8660254037844387D0)/,
     6         PHSC16/( 0.8660254037844387D0,-0.5D0)/,
     7         PHSC56/(-0.8660254037844387D0,-0.5D0)/,
     8          COMPI/(0.0D0,1.0D0)/         



      ZPP=Z

C   * PRELIMINARIES *

      IF (IRETIN.EQ.IRETOUT)  RETURN              
      ZETA = Z * SQRT(Z) / 1.5D0

C   * SCALED AI & BI OUPUT REQUESTED *

      IF (IRETOUT.LT.0)  THEN

          IF (IRETIN.EQ.0)  THEN                  
              EXPZETA = EXP (ZETA)
              AI  = AI * EXPZETA
              AIP = AIP * EXPZETA
              IF (DREAL(ZETA).GE.0.D0)  THEN
                  EXPMZETA = EXP (-ZETA)
                  BI  = BI * EXPMZETA
                  BIP = BIP * EXPMZETA
              ELSE
                  BI  = BI * EXPZETA
                  BIP = BIP * EXPZETA
              END IF

          ELSE IF (IRETIN.GT.0)  THEN             
              IF (ZPR(2).GE.0.D0)  THEN         
                  IF (DREAL(ZETA).GE.0.D0)  THEN              
                     EXPTZETA = EXP (-2.D0*ZETA)
                     BI  = 2.D0*PHSC16*BI  + COMPI*AI*EXPTZETA
                     BIP = 2.D0*PHSC56*BIP + COMPI*AIP*EXPTZETA
                  ELSE                   
                     EXPTZETA = EXP (2.D0*ZETA)
                     BI  = 2.D0*PHSC16*BI*EXPTZETA  + COMPI*AI
                     BIP = 2.D0*PHSC56*BIP*EXPTZETA + COMPI*AIP
                  END IF
              ELSE
                  IF (DREAL(ZETA).GE.0.D0)  THEN              
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = 2.D0*PHSE16*BI  - COMPI*AI*EXPTZETA
                    BIP = 2.D0*PHSE56*BIP - COMPI*AIP*EXPTZETA
                  ELSE                   
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = 2.D0*PHSE16*BI*EXPTZETA  - COMPI*AI
                    BIP = 2.D0*PHSE56*BIP*EXPTZETA - COMPI*AIP
                  END IF
              END IF
          END IF

C   * UNSCALED AI-BI OUTPUT REQUESTED *

        ELSE IF (IRETOUT.EQ.0)  THEN
          EXPMZETA = EXP (-ZETA)
          AI  = AI * EXPMZETA
          AIP = AIP * EXPMZETA

          IF (IRETIN.LT.0)  THEN                  
              IF (DREAL(ZETA).GE.0.D0) THEN
                  EXPZETA = EXP (ZETA)
                  BI  = BI * EXPZETA
                  BIP = BIP * EXPZETA
              ELSE
                  BI  = BI * EXPMZETA
                  BIP = BIP * EXPMZETA
              END IF

          ELSE IF (IRETIN.GT.0) THEN              
              EXPZETA = EXP (ZETA)
              IF (ZPR(2).GE.0.D0)  THEN
                  BI  = 2.D0*PHSC16  * BI*EXPZETA  + COMPI*AI
                  BIP = 2.D0*PHSC56  * BIP*EXPZETA + COMPI*AIP
              ELSE
                  BI  = 2.D0*PHSE16  * BI*EXPZETA  - COMPI*AI
                  BIP = 2.D0*PHSE56  * BIP*EXPZETA - COMPI*AIP
              END IF
          END IF
          ZETA = (0.0D0, 0.0D0)

C   * NUMERICALLY STABLE SCALED AI OUTPUT REQUESTED *
C     (OBVIOUSLY IF YOU ARE WORKING FROM UNSCALED OR SCALED AI/BI
C      SOME PRECISION LOSS MUST OCCUR, BE FOREWARNED)

      ELSE IF (IRETOUT.GT.0)  THEN

          IF (IRETIN.LT.0)  THEN                  
              IF (ZPR(2).GE.0.D0)  THEN         
                  IF (DREAL(ZETA).GE.0.D0)  THEN             
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = PHSE16*(BI  - COMPI*AI*EXPTZETA)*0.5D0
                    BIP = PHSE56*(BIP - COMPI*AIP*EXPTZETA)*0.5D0
                  ELSE                  
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = PHSE16*(BI*EXPTZETA  - COMPI*AI)*0.5D0
                    BIP = PHSE56*(BIP*EXPTZETA - COMPI*AIP)*0.5D0
                  END IF
              ELSE
                  IF (DREAL(ZETA).GE.0.D0)  THEN             
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = PHSC16*(BI  + COMPI*AI*EXPTZETA)*0.5D0
                    BIP = PHSC56*(BIP + COMPI*AIP*EXPTZETA)*0.5D0
                  ELSE                  
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = PHSC16*(BI*EXPTZETA  + COMPI*AI)*0.5D0
                    BIP = PHSC56*(BIP*EXPTZETA + COMPI*AIP)*0.5D0
                  END IF
              END IF

          ELSE IF (IRETIN.EQ.0)  THEN             
              EXPZETA = EXP (ZETA)
              EXPMZETA = EXP (-ZETA)
              IF (ZPR(2).GE.0.D0)  THEN
                  BI  = PHSE16*(BI  - COMPI*AI)*EXPMZETA*0.5D0
                  BIP = PHSE56*(BIP - COMPI*AIP)*EXPMZETA*0.5D0
              ELSE
                  BI  = PHSC16*(BI  + COMPI*AI)*EXPMZETA*0.5D0
                  BIP = PHSC56*(BIP + COMPI*AIP)*EXPMZETA*0.5D0
              END IF
              AI  = AI*EXPZETA
              AIP = AIP*EXPZETA
          END IF
      END IF
      RETURN
      END

