      subroutine scairy2(z1,ai1,aip1,bi1,bip1,zeta1,det1,z2,ai2,aip2,
     .   bi2,bip2,zeta2,det2,aisoln,iibad)
c
c: This subroutine computes the two numerically stable solutions to the
c: Airy differential equation for the two arguments z1 and z2.
c: The derivatives with respect to z1 and z2 are also computed.
c: For abs(arg(z)) < pi/3, the two solutions are Ai and Bi.
c: For pi/3 <= arg(x) <= pi, the two sol'n are Ai(z) and Ai(z*exp(-2pi/3)).
c: For -pi <= arg(x) <= pi/3, the two sol'n are Ai(z) and Ai(z*exp(2pi/3)).
c: Uses the Airy function routines from SAFARI, where one bug has been 
c: fixed in CLAIRY (see EKW note).  Note also that the derivatives of
c: Ai(z*exp(2pi/3)) are with respect to z, not the total argument.
c: Special care is taken that when z1 and z2 straddle the negative real
c: axis, the same pair of independent solutions are computed.
c
      IMPLICIT COMPLEX*16 (A,B,D,Z)
      include 'scairy_com'
      integer*4 aisoln,iibad
c
      iibad=0
      call clairy(z1,1,ai1,bi1,aip1,bip1,zeta1)
c: EKW changed: For imag(z)>0, we want (d/dz)(Ai(z*eim23)), whereas 
c: CLAIRY is giving BIP = (d/d(z*eim23))(Ai(z*eim23)), 
c: so multiply BIP by eim23:
      IF (dimag(z1).GE.0.D0)  THEN
         bip1=eim23*bip1
         det1=det_pos
         aisoln=1
      ELSE
         bip1=ei23*bip1
         det1=det_neg
         aisoln=-1
      END IF
c
      call clairy(z2,1,ai2,bi2,aip2,bip2,zeta2)
c: EKW changed: For imag(z)>0, we want (d/dz)(Ai(z*eim23)), whereas 
c: CLAIRY is giving BIP = (d/d(z*eim23))(Ai(z*eim23)), 
c: so multiply BIP  by eim23:
      IF (dimag(z2).GE.0.D0)  THEN
         bip2=eim23*bip2
         det2=det_pos
      ELSE
         bip2=ei23*bip2
         det2=det_neg
      END IF
c
c: Check if z1 and z2 straddle the real axis:
      if(dimag(z1) .ge. 0. .and. dimag(z2) .lt. 0.) then
         if(dreal(zeta2) .gt. dreal(zeta1)) then
c: z2 on neg side, so change Ai(z2*ei23) to Ai(z2*eim23):
            call ai_flip(zeta2,ai2,aip2,bi2,bip2,ei23,eim23,
     .         det2,det_pos,iibad)
         else
c: z1 on pos side, so change Ai(z1*eim23) to Ai(z1*ei23):
            call ai_flip(zeta1,ai1,aip1,bi1,bip1,eim23,ei23,
     .         det1,det_neg,iibad)
            aisoln=-1
         endif
      elseif(dimag(z2) .ge. 0. .and. dimag(z1) .lt. 0.) then
         if(dreal(zeta1) .gt. dreal(zeta2)) then
c: z1 on neg side, so change Ai(z1*ei23) to Ai(z1*eim23):
            call ai_flip(zeta1,ai1,aip1,bi1,bip1,ei23,eim23,
     .         det1,det_pos,iibad)
            aisoln=1
         else
c: z2 on pos side, so change Ai(z2*eim23) to Ai(z2*ei23):
            call ai_flip(zeta2,ai2,aip2,bi2,bip2,eim23,ei23,
     .         det2,det_neg,iibad)
         endif
      endif
c
      RETURN
      END
ccc
      subroutine ai_flip(zeta2,ai2,aip2,bi2,bip2,ei23,eim23,
     .   det2,det_pos,iibad)
c
      implicit none
      integer*4 iibad
      complex*16 zeta2,ai2,aip2,bi2,bip2,ei23,eim23,det2,det_pos,zfac,
     .   diff
c
      if(dreal(zeta2) .gt. 0.d0) then
         zfac=cdexp(-2.d0*zeta2)
         bi2=-(ei23*ai2*zfac + eim23*bi2)
         bip2=-(ei23*aip2*zfac + eim23*bip2)
      elseif(dreal(zeta2) .lt. -50.d0) then
c: Need to split into two layers in ai_strad later:
         iibad=1
      else
cxx      print *,'ai_flip for re(zeta)<0: ',zeta2,ai2,aip2,bi2,bip2
         zfac=cdexp(-zeta2)
         bi2=-(ei23*ai2*zfac + eim23*bi2/zfac)
         bip2=-(ei23*aip2*zfac + eim23*bip2/zfac)
         ai2=ai2*zfac
         aip2=aip2*zfac
         det2=ai2*bip2 - aip2*bi2
         diff=det_pos - det2
         if(dreal(diff)*dreal(diff)+dimag(diff)*dimag(diff) .gt. 
     .      1.e-12) iibad=1
         zeta2=(0.d0,0.d0)
      endif
      det2=det_pos
cxx   call wronsk_ch(ai2,aip2,bi2,bip2,det2)
c
      return
      end
ccc
      subroutine wronsk_ch(ai2,aip2,bi2,bip2,det2)
c
      implicit none
      complex*16 ai2,aip2,bi2,bip2,det2,det
c
      det=ai2*bip2 - aip2*bi2
      if(abs(det-det2) .gt. 1.e-4) then
         print *,'W = ',det,det2
      endif
c
      return
      end
ccc
      subroutine scairy3(nz,z,z0,ai,aip,bi,bip,zeta,det,aisoln)
c
c: Identical to scairy2, except computes a set of nz airy functions
c: at z(1:nz), assuring that all are of the same solution.
c
      implicit none
      integer*4 nz,jz,aisoln,iibad
      complex*16 z(nz),z0,ai(nz),aip(nz),bi(nz),bip(nz),zeta(nz),det
      include 'scairy_com'
c
cc    if(dimag(z0).GE.0.D0) then
      if(aisoln .eq. 1) then
         do jz=1,nz
            call clairy(z(jz),1,ai(jz),bi(jz),aip(jz),bip(jz),zeta(jz))
            bip(jz)=eim23*bip(jz)
c: z on neg side, so change Ai(z*ei23) to Ai(z*eim23):
            if(dimag(z(jz)) .lt. 0.d0) then
               call ai_flip(zeta(jz),ai(jz),aip(jz),bi(jz),bip(jz),
     .            ei23,eim23,det,det_pos,iibad)
            endif
         enddo
         det=det_pos
      else
         do jz=1,nz
            call clairy(z(jz),1,ai(jz),bi(jz),aip(jz),bip(jz),zeta(jz))
            bip(jz)=ei23*bip(jz)
c: z on pos side, so change Ai(z*eim23) to Ai(z*ei23):
            if(dimag(z(jz)) .gt. 0.d0) then
               call ai_flip(zeta(jz),ai(jz),aip(jz),bi(jz),bip(jz),
     .            eim23,ei23,det,det_neg,iibad)
            endif
         enddo
         det=det_neg
      endif
c
      return
      end
ccc
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
      include 'scairy_com'
      EQUIVALENCE (ZPR(1),ZPP)

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
          CALL ABSERIES2 (Z, AI,AIP,BI,BIP)
          ZETA = (0.0D0, 0.0D0)
          IF (IRET.NE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,0,IRET)


C   ******************************************************************
C   *    *
C   *   ANGLES LESS THAN 120 DEGS. - AI(Z) AND AI(Z*EXP(-2IPI/3))    *
C   *   ARE CALCULATED FOR Z IN THE UPPER HALF PLANE, CONJGATE OF    *
C   *   Z IF IN THE LOWER.  SERIES OR GAUSSIAN SOLUTION AS NEEDED.   *
C   *    *
C   ******************************************************************

c: EKW: This was a bug originally because ZPR is equivalenced to original
c: input value of Z, not ZP, which is always in upper half plane:
cxx   ELSE IF (ATAN2(ZPR(2),ZPR(1)) .LE. TPIO3)  THEN
cxx   ELSE IF (abs(ATAN2(ZPR(2),ZPR(1))) .LE. TPIO3)  THEN
      ELSE IF (dimag(zp) .gt. -1.73205080756888d0*dreal(zp)) then
          IF (ABS(ZP).GT.7.9111D0)  THEN          
              CALL ABGAUSS2 (ZP, AI, AIP, BI, BIP, ZETA)
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  BI  = CONJG (BI)
                  AIP = CONJG (AIP)
                  BIP = CONJG (BIP)
              END IF
          ELSE
              IF (ABS (ZP-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
                  CALL ASERIES2 (ZP, AI,AIP, ZETA)
              ELSE
                  CALL AGAUSS2 (ZP, AI, AIP, ZETA)
              END IF
              ZP2 = CONJG (ZP) * ei23
              IF (ABS (ZP2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
                  CALL ASERIES2 (ZP2, BI, BIP, ZETA2)
              ELSE
                  CALL AGAUSS2 (ZP2, BI, BIP, ZETA2)
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
c: z1=z*exp(-i2pi/3); z2=conjg[z*exp(i2pi/3] (to make z2 in upper half plane;
c: take conjugate of ai,aip,zeta afterwards):
          Z1 = Z * eim23    
          Z2 = CONJG (Z) * eim23                
          
          IF (ABS (Z1-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z1, AI1, AIP1, ZETA1)
          ELSE
              CALL AGAUSS2 (Z1, AI1, AIP1, ZETA1)
          END IF
          IF (ABS (Z2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z2, AI2, AIP2, ZETA2)
          ELSE
              CALL AGAUSS2 (Z2, AI2, AIP2, ZETA2)
          END IF
          AI2  = CONJG (AI2)  
          AIP2 = CONJG (AIP2)
          ZETA2 = CONJG (ZETA2)

c: Use connection formula: 
c: Ai(z)=exp(ipi/3)*Ai(z*exp(-i2pi/3)) + exp(-ipi/3)*z*exp(i2pi/3)
          IF (IRET.EQ.0)  THEN
              E1 = EXP (-ZETA1)
              E2 = EXP (-ZETA2)
              AI  = ei13 * AI1*E1  + eim13 * AI2*E2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
              BI  = ei16 * AI2*E2  + eim16 * AI1*E1
              BIP = ei56 * AIP2*E2 + eim56 * AIP1*E1
              ZETA = (0.0D0, 0.0D0)

          ELSE IF (IRET.LT.0)  THEN               
              ZETA  = Z * SQRT(Z) * ONEOOF
              E1  = EXP (ZETA-ZETA1)
              E2  = EXP (ZETA-ZETA2)
              AI  = ei13 * AI1*E1  + eim13 * AI2*E2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
              BI  = eim16 * AI1*E1  + ei16 * AI2*E2
              BIP = eim56 * AIP1*E1 + ei56 * AIP2*E2

          ELSE IF (IRET.GT.0)  THEN               
              ZETA  = Z * SQRT(Z) * ONEOOF
              E1  = EXP (ZETA-ZETA1)
              E2  = EXP (ZETA-ZETA2)
              AI  = ei13 * AI1*E1  + eim13 * AI2*E2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
cekw          IF (ZPR(2).GE.0.0D0)  THEN
              IF (ICONJG .eq. 0) then
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
      SUBROUTINE ABSERIES2 (Z, AI,AIP,BI,BIP)
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
      SUBROUTINE ABGAUSS2 (Z, AI,AIP,BI,BIP,ZETA)

      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      PARAMETER (ONEOOF=1D0/1.5D0)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (ONEOPI=1.0D0/PI) 
      COMPLEX*16 SQRTZ,SQRPIZ,CFAC,CFAC2,ONEOZ
      REAL*8 ZPR(2)
      EQUIVALENCE (ZPR(1),ZPP)
      REAL*8 X4(4), W4(4), X6(6), W6(6)
      include 'scairy_com'

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
          BI  = ei16 * BI * 0.5D0
          BIP = ei56 * BIP * 0.5D0
      ELSE
          BI  = eim16 * BI * 0.5D0
          BIP = eim56 * BI * 0.5D0
      END IF
      RETURN
      END
      SUBROUTINE ASERIES2 (Z, AI, AIP, ZETA)
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
ccc
      SUBROUTINE AGAUSS2 (Z,AI,AIP,ZETA)
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
      COMPLEX*16 COMPI
      include 'scairy_com'
      DATA     COMPI/(0.0D0,1.0D0)/         

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
              IF (REAL(ZETA).GE.0.D0)  THEN
                  EXPMZETA = EXP (-ZETA)
                  BI  = BI * EXPMZETA
                  BIP = BIP * EXPMZETA
              ELSE
                  BI  = BI * EXPZETA
                  BIP = BIP * EXPZETA
              END IF

          ELSE IF (IRETIN.GT.0)  THEN             
              IF (ZPR(2).GE.0.D0)  THEN         
                  IF (REAL(ZETA).GE.0.D0)  THEN              
                     EXPTZETA = EXP (-2.D0*ZETA)
                     BI  = 2.D0*eim16*BI  + COMPI*AI*EXPTZETA
                     BIP = 2.D0*eim56*BIP + COMPI*AIP*EXPTZETA
                  ELSE                   
                     EXPTZETA = EXP (2.D0*ZETA)
                     BI  = 2.D0*eim16*BI*EXPTZETA  + COMPI*AI
                     BIP = 2.D0*eim56*BIP*EXPTZETA + COMPI*AIP
                  END IF
              ELSE
                  IF (REAL(ZETA).GE.0.D0)  THEN              
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = 2.D0*ei16*BI  - COMPI*AI*EXPTZETA
                    BIP = 2.D0*ei56*BIP - COMPI*AIP*EXPTZETA
                  ELSE                   
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = 2.D0*ei16*BI*EXPTZETA  - COMPI*AI
                    BIP = 2.D0*ei56*BIP*EXPTZETA - COMPI*AIP
                  END IF
              END IF
          END IF

C   * UNSCALED AI-BI OUTPUT REQUESTED *

        ELSE IF (IRETOUT.EQ.0)  THEN
          EXPMZETA = EXP (-ZETA)
          AI  = AI * EXPMZETA
          AIP = AIP * EXPMZETA

          IF (IRETIN.LT.0)  THEN                  
              IF (REAL(ZETA).GE.0.D0) THEN
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
                  BI  = 2.D0*eim16  * BI*EXPZETA  + COMPI*AI
                  BIP = 2.D0*eim56  * BIP*EXPZETA + COMPI*AIP
              ELSE
                  BI  = 2.D0*ei16  * BI*EXPZETA  - COMPI*AI
                  BIP = 2.D0*ei56  * BIP*EXPZETA - COMPI*AIP
              END IF
          END IF
          ZETA = (0.0D0, 0.0D0)

C   * NUMERICALLY STABLE SCALED AI OUTPUT REQUESTED *
C     (OBVIOUSLY IF YOU ARE WORKING FROM UNSCALED OR SCALED AI/BI
C      SOME PRECISION LOSS MUST OCCUR, BE FOREWARNED)

      ELSE IF (IRETOUT.GT.0)  THEN

          IF (IRETIN.LT.0)  THEN                  
              IF (ZPR(2).GE.0.D0)  THEN         
                  IF (REAL(ZETA).GE.0.D0)  THEN             
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = ei16*(BI  - COMPI*AI*EXPTZETA)*0.5D0
                    BIP = ei56*(BIP - COMPI*AIP*EXPTZETA)*0.5D0
                  ELSE                  
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = ei16*(BI*EXPTZETA  - COMPI*AI)*0.5D0
                    BIP = ei56*(BIP*EXPTZETA - COMPI*AIP)*0.5D0
                  END IF
              ELSE
                  IF (REAL(ZETA).GE.0.D0)  THEN             
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = eim16*(BI  + COMPI*AI*EXPTZETA)*0.5D0
                    BIP = eim56*(BIP + COMPI*AIP*EXPTZETA)*0.5D0
                  ELSE                  
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = eim16*(BI*EXPTZETA  + COMPI*AI)*0.5D0
                    BIP = eim56*(BIP*EXPTZETA + COMPI*AIP)*0.5D0
                  END IF
              END IF

          ELSE IF (IRETIN.EQ.0)  THEN             
              EXPZETA = EXP (ZETA)
              EXPMZETA = EXP (-ZETA)
              IF (ZPR(2).GE.0.D0)  THEN
                  BI  = ei16*(BI  - COMPI*AI)*EXPMZETA*0.5D0
                  BIP = ei56*(BIP - COMPI*AIP)*EXPMZETA*0.5D0
              ELSE
                  BI  = eim16*(BI  + COMPI*AI)*EXPMZETA*0.5D0
                  BIP = eim56*(BIP + COMPI*AIP)*EXPMZETA*0.5D0
              END IF
              AI  = AI*EXPZETA
              AIP = AIP*EXPZETA
          END IF
      END IF
      RETURN
      END
ccc
      SUBROUTINE airy_only (Z,AI,AIP,ZETA)
      IMPLICIT COMPLEX*16 (A,B,E,X,Z)
      PARAMETER (ONEOOF=1D0/1.5D0)

c: Modified by EKW to only return Ai(z) assuming IRET=1 below
c: for use in h_space, where only one solution is required.
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
      include 'scairy_com'
      EQUIVALENCE (ZPR(1),ZPP)

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
cekw      CALL ABSERIES2 (Z, AI,AIP,BI,BIP)
          call airy_only_series (Z, AI,AIP,ZETA)
cekw      IF (IRET.NE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,0,IRET)


C   ******************************************************************
C   *    *
C   *   ANGLES LESS THAN 120 DEGS. - AI(Z) AND AI(Z*EXP(-2IPI/3))    *
C   *   ARE CALCULATED FOR Z IN THE UPPER HALF PLANE, CONJGATE OF    *
C   *   Z IF IN THE LOWER.  SERIES OR GAUSSIAN SOLUTION AS NEEDED.   *
C   *    *
C   ******************************************************************

c: EKW: This was a bug originally because ZPR is equivalenced to original
c: input value of Z, not ZP, which is always in upper half plane:
cxx   ELSE IF (ATAN2(ZPR(2),ZPR(1)) .LE. TPIO3)  THEN
cxx   ELSE IF (abs(ATAN2(ZPR(2),ZPR(1))) .LE. TPIO3)  THEN
      ELSE IF (dimag(zp) .gt. -1.73205080756888d0*dreal(zp)) then
          IF (ABS(ZP).GT.7.9111D0)  THEN          
cekw          CALL ABGAUSS2 (ZP, AI, AIP, BI, BIP, ZETA)
              call AGAUSS2 (ZP,AI,AIP,ZETA)
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  AIP = CONJG (AIP)
              END IF
          ELSE
              IF (ABS (ZP-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
cekw              CALL ASERIES2 (ZP, AI,AIP, ZETA)
                  call airy_only_series(ZP, AI,AIP,ZETA)
              ELSE
                  CALL AGAUSS2 (ZP, AI, AIP, ZETA)
              END IF
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  AIP = CONJG (AIP)
              END IF
          END IF
cekw      IF (IRET.LE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,1,IRET)

C   **********************************************************
C   *                *
C   *   CONNECTION REGION FOR ANGLE GREATER THAN 120 DEGS.   *
C   *   CALCULATE AI(Z*EXP(-2IPI/3)) AND AI(Z*EXP(2IPI/3))   *
C   *   SOLUTIONS, USE CONNECTION FORMULAS FOR RESULTS.      *
C   *                *
C   **********************************************************

      ELSE 
c: z1=z*exp(-i2pi/3); z2=conjg[z*exp(i2pi/3] (to make z2 in upper half plane;
c: take conjugate of ai,aip,zeta afterwards):
          Z1 = Z * eim23    
          Z2 = CONJG (Z) * eim23                
          
          IF (ABS (Z1-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z1, AI1, AIP1, ZETA1)
          ELSE
              CALL AGAUSS2 (Z1, AI1, AIP1, ZETA1)
          END IF
          IF (ABS (Z2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z2, AI2, AIP2, ZETA2)
          ELSE
              CALL AGAUSS2 (Z2, AI2, AIP2, ZETA2)
          END IF
          AI2  = CONJG (AI2)  
          AIP2 = CONJG (AIP2)
          ZETA2 = CONJG (ZETA2)

c: Use connection formula: 
c: Ai(z)=exp(ipi/3)*Ai(z*exp(-i2pi/3)) + exp(-ipi/3)*z*exp(i2pi/3)
cekw          ZETA  = Z * SQRT(Z) * ONEOOF
cekw          E1  = EXP (ZETA-ZETA1)
cekw          E2  = EXP (ZETA-ZETA2)
cekw          AI  = ei13 * AI1*E1  + eim13 * AI2*E2
cekw          AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
cekw: Take out exp(zeta1) or exp(zeta2), whichever has smaller real part:
           if(dreal(zeta1) .lt. dreal(zeta2)) then
              zeta=zeta1
              E2=exp(zeta-zeta2)
              AI  = ei13 * AI1  + eim13 * AI2*E2
              AIP = eim13 * AIP1 + ei13 * AIP2*E2
           else
              zeta=zeta2
              E1=exp(zeta-zeta1)
              AI  = ei13 * AI1*E1  + eim13 * AI2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2
           endif
      END IF

      RETURN
      END
ccc
      SUBROUTINE airy_only_series (Z, AI,AIP,ZETA)
      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      REAL*8 C1,C2,FAC1,FAC2,FAC3
      DATA C1 /0.355028053887817D0/, C2 /0.258819403792807D0/

      ZETA = (0.0D0, 0.0D0)
c
      FC = (1.0D0, 0.0D0)
      GC = Z
      Z2 = Z * Z
      AI = C1*FC - C2*GC
      AIP = -C2

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
     &      ABS(GC).LT.1.D-17*ABS(AI))  GO TO 12
10    CONTINUE
      WRITE (6,*) ' SERIES APPROXIMATION FAILED TO CONVERGE'

12    continue
c
      RETURN
      END
