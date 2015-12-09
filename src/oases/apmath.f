C*********************************************
C
C     SIMMULATION OF APMATH64 LIBRARY ROUTINES
C     AND RELATED APAL64 ROUTINES.
C
C     H. SCHMIDT    SACLANTCEN    830217
C
C*********************************************
      SUBROUTINE CFFT(DATA,N,IFORW)
      DIMENSION DATA(2*N)
      ISIGN=-IFORW
      IP0=2                                                             FF1  23 
      IP3=IP0*N                                                         FF1  24 
      I3REV=1                                                           FF1  25 
      DO 50 I3=1,IP3,IP0                                                FF1  26 
      IF(I3-I3REV)10,20,20                                              FF1  27 
10    TEMPR=DATA(I3)                                                    FF1  28 
      TEMPI=DATA(I3+1)                                                  FF1  29 
      DATA(I3)=DATA(I3REV)                                              FF1  30 
      DATA(I3+1)=DATA(I3REV+1)                                          FF1  31 
      DATA(I3REV)=TEMPR                                                 FF1  32 
      DATA(I3REV+1)=TEMPI                                               FF1  33 
20    IP1=IP3/2                                                         FF1  34 
30    IF(I3REV-IP1)50,50,40                                             FF1  35 
40    I3REV=I3REV-IP1                                                   FF1  36 
      IP1=IP1/2                                                         FF1  37 
      IF(IP1-IP0)50,30,30                                               FF1  38 
50    I3REV=I3REV+IP1                                                   FF1  39 
      IP1=IP0                                                           FF1  40 
60    IF(IP1-IP3)70,100,100                                             FF1  41 
70    IP2=IP1*2                                                         FF1  42 
      THETA=6.283185307/FLOAT(ISIGN*IP2/IP0)                                    
      SINTH=SIN(THETA/2.)                                               FF1  44 
      WSTPR=-2.*SINTH*SINTH                                             FF1  45 
      WSTPI=SIN(THETA)                                                  FF1  46 
      WR=1.                                                             FF1  47 
      WI=0.                                                             FF1  48 
      DO 90 I1=1,IP1,IP0                                                FF1  49 
      DO 80 I3=I1,IP3,IP2                                               FF1  50 
      I2A=I3                                                            FF1  51 
      I2B=I2A+IP1                                                       FF1  52 
      TEMPR=WR*DATA(I2B)-WI*DATA(I2B+1)                                 FF1  53 
      TEMPI=WR*DATA(I2B+1)+WI*DATA(I2B)                                 FF1  54 
      DATA(I2B)=DATA(I2A)-TEMPR                                         FF1  55 
      DATA(I2B+1)=DATA(I2A+1)-TEMPI                                     FF1  56 
      DATA(I2A)=DATA(I2A)+TEMPR                                         FF1  57 
80    DATA(I2A+1)=DATA(I2A+1)+TEMPI                                     FF1  58 
      TEMPR=WR                                                          FF1  59 
      WR=WR*WSTPR-WI*WSTPI+WR                                           FF1  60 
90    WI=WI*WSTPR+TEMPR*WSTPI+WI                                        FF1  61 
      IP1=IP2                                                           FF1  62 
      GO TO 60                                                          FF1  63 
100   RETURN                                                            FF1  64 
      END                                                               FF1  65 
      SUBROUTINE RFT(A,B,N,ISN)                                                 
C    C,10/11/80.V1 RTD <PN> :FFT LIBRARY                                        
C    IF ISN=1,THIS SUBROUTINE COMPLETES THE FOURIER TRANSFORM                   
C     OF 2*N REAL DATA VALUES,WHERE THE ORIGINAL DATA VALUES ARE                
C    STORED ALTERNATELY IN ARRAYS A AND B, AND ARE FIRST                        
C    TRANSFORMED BY A COMPLEX FOURIER TRANSFORM OF DIMENSION N.                 
C    THE COSINE COEFFICIENTS ARE IN A(1),A(2),...A(N+1) AND                     
C     THE SINE COEFFICIENTS ARE IN B(1),B(2),...B(N+1).                         
C    A TYPICAL CALLING SEQUENCE IS                                              
C         CALL FFT(A,B,N,N,N,1)                                                 
C         CALL RFT(A,B,N,1)                                                     
C    THE RESULTS SHOULD BE MULTIPLIED BY 0.5 /N TO GIVE THE                     
C    USUAL SCALING OF COEFFICIENTS.                                             
C    IF ISN=-1, THE INVERSE TRANSFORMATION IS DONE, THE FIRST                   
C    STEP IN EVALUATING A REAL FOURIER SERIES.                                  
C    A TYPICAL CALLING SEQUENCE IS                                              
C         CALL RFT(A,B,N,-1)                                                    
C         CALL FFT(A,B,N,N,N,-1)                                                
C    THE RESULTS SHOULD BE MULTIPLIED BY 0.5 TO GIVE THE USUAL                  
C    SCALING , AND THE TIME DOMAIN RESULTS ALTERNATE IN ARRAYS                  
C     A AND B, I.E.,A(1),B(1),A(2),B(2),...A(N),B(N).                           
C    THE DATA MAY ALTERNATIVELY BE STORED IN A SINGLE COMPLEX                   
C     ARRAY A, THEN THE MAGNITUDE OF ISN CHANGED TO TWO TO                      
C    GIVE THE CORRECT INDEXING INCREMENT AND A(2) USED TO PASS                  
C    THE INITIAL ADDRESS FOR THE SEQUENCE OF IMAGINARY                          
C    VALUES,E.G.                                                                
C          CALL FFT(A,A(2),N,N,N,2)                                             
C          CALL RFT(A,A(2),N,2)                                                 
C    IN THIS CASE, THE COSINE AND SINE COEFFICIENTS ALTERNATE IN                
C    A,  BY R.C. SINGLETON, STANFORD RESEARCH INSTITUTE, OCT. 1968              
C                                                                               
C NOTE FROM THE USERS:                                                          
C                                                                               
C     BUGS WERE FOUND IN THIS VERSION OF RFT, WHICH HAS BEEN USED IN            
C    BASIC I,II,AND III.   THESE ARE:                                           
C                                                                               
C      1)  TWO VALUES OUTSIDE THE USER'S BUFFER, X(N+1) AND X(N+2) FOR          
C     A TRANSFORM OF N DATA POINTS, ARE CLOBBERED IN AN ATTEMPT AT              
C     SCRATCH STORAGE.  THIS HAS BEEN CORRECTED IN THE LATER VERSION, AN        
C     THE OLD PROGRAM, IDENTICAL IN EXECUTION, BUT WHICH DID NOT CONTAIN        
C     CORREKT THINKING IN ITS NONUSE OF COMMENTS, HAS BEEN DESTROYED.           
C                                                                               
C     2)  THE INVERSE TRANSFORM STATEMENTS SHOULD BE CHANGED, SO THAT           
C     A(K),A(K+1) (BETWEEN STATEMENTS 5 & 10) IS DEFINED.                       
C                                                                               
C     3)   ROUNDED ARITHMETIC SHOULD BE USED.                                   
C                                                                               
       DIMENSION A(1),B(1)                                                      
       INC=IABS(ISN)                                                            
       NK=N*INC+2                                                               
       NH=NK/2                                                                  
       SD=-2.0*ATAN(1.0)/FLOAT(N)                                               
       CD=2.0*SIN(SD)**2                                                        
       SD=SIN(SD+SD)                                                            
      IF(ISN.LT.0)GOTO 30                                                       
      CN=1.                                                                     
C DO FIRST LOOP                                                                
C  STORING NYQUIST REAL IN IMAGINARY 0 SAMPLE & VICEVERSA                       
10    AA  =2.*(A(1)+B(1))                                                       
      B(1)=2.*(A(1)-B(1))                                                       
C                                                                               
12    A(1)=AA                                                                   
C INCREMENT SIN & COSINE                                                        
      SN=SD*CN                                                                  
      CN=CN*(1.-CD)                                                             
      L=1+INC                                                                   
C MAIN LOOP                                                                     
      DO 20 J=L,NH,INC                                                          
C OK, FROM HERE ON IT'S ALMOST THE SAME.                                        
        K=NK-J                                                                  
       AA=A(J)+A(K)                                                             
       AB=A(J)-A(K)                                                             
       BA=B(J)+B(K)                                                             
       BB=B(J)-B(K)                                                             
        RE=CN*BA+SN*AB                                                          
       FIM=SN*BA-CN*AB                                                          
       B(K)=FIM-BB                                                              
        B(J)=FIM+BB                                                             
       A(K)=AA-RE                                                               
       A(J)=AA+RE                                                               
       AA=CN-(CD*CN+SD*SN)                                                      
       SN=(SD*CN-CD*SN)+SN                                                      
C    THE FOLLOWING THREE STATEMENTS COMPENSATE FOR TRUNCATION                   
C     ERROR.  IF ROUNDED ARITHMETIC IS USED, SUBSTITUTE                         
   20  CN=AA                                                                    
C       CN=0.5/(AA**2+SN**2)+0.5                                                
C      SN=CN*SN                                                                 
C20    CN=CN*AA                                                                 
       RETURN                                                                   
C INVERSE TRANSFORM                                                             
 30    CN=-1.0                                                                  
       SD=-SD                                                                   
      AA  = A(1)+B(1)                                                           
      B(1)= A(1)-B(1)                                                           
       GO TO 12                                                                 
       END                                                                      
      SUBROUTINE RFFT(A,NX,IFORW)
      DIMENSION A(NX)
      IF (IFORW.GT.0) THEN
      CALL CFFT(A,NX/2,1)
      CALL RFT(A(1),A(2),NX/2,2)
      ELSE
      CALL RFT(A(1),A(2),NX/2,-2)
      CALL CFFT(A,NX/2,-1)
      END IF
      RETURN
      END
      SUBROUTINE RFFTSC(A,NX,IPACK,IFAC)
      DIMENSION A(*)
      IF (IFAC.EQ.1) THEN
      FAC=.5/NX
      ELSE IF (IFAC.EQ.-1) THEN
      FAC=0.25/NX
      ELSE
      FAC=1E0
      END IF
      IF (IPACK.EQ.2) THEN
      A(2)=0E0
      ELSE IF (IPACK.EQ.3) THEN
      A(NX+1)=A(2)
      A(NX+2)=0E0
      A(2)=0E0
      ELSE IF (IPACK.EQ.-3) THEN
      A(2)=A(NX+1)
      ELSE IF (IPACK.EQ.-2) THEN
      A(2)=0E0
      ELSE
      END IF
      IF (IFAC.NE.0) THEN
      DO 10 I=1,NX
 10   A(I)=A(I)*FAC
      END IF
      RETURN
      END
C
C
      SUBROUTINE VSMUL(A,I,C,D,J,N)
      DIMENSION A(1),D(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      D(J1)=C*A(I1)
      I1=I1+I
      J1=J1+J
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VMAX(A,I,C,N)
      DIMENSION A(1),D(1)
      I1=1
      C=-1E30
      DO 10 L=1,N
      C=MAX(C,A(I1))
      I1=I1+I
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VMIN(A,I,C,N)
      DIMENSION A(1),D(1)
      I1=1
      C=1E30
      DO 10 L=1,N
      C=MIN(C,A(I1))
      I1=I1+I
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VNEG(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
CDEC$ INIT_DEP_FWD
cvd$  nodepchk
      DO 10 L=1,N
      C(J1)=-A(I1)
      I1=I1+I
 10   J1=J1+J
      RETURN
      END
C
C
      SUBROUTINE VMOV(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=A(I1)
      I1=I1+I
 10   J1=J1+J
      RETURN
      END
C
C
      SUBROUTINE VSUM(A,I,C,J,N,H)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
      SUM=0E0
      HH=H
      DO 10 L=1,N
      SUM=SUM+A(I1)*HH
      C(J1)=SUM
      I1=I1+I
 10   J1=J1+J
      RETURN
      END
C
C
      SUBROUTINE VTRAPZ(A,I,C,J,N,H)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
      SUM=0E0
      HH=.5*H
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N-1
      SUM=SUM+HH*(A(I1)+A(I1+I))
      C(J1)=SUM
      I1=I1+I
 10   J1=J1+J
      c(1+(n-1)*j)=sum
      RETURN
      END
C
C
      SUBROUTINE VSQ(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=A(I1)*A(I1)
      I1=I1+I
 10   J1=J1+J
      RETURN
      END
C
C
      SUBROUTINE VSQRT(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=SQRT(A(I1))
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VEXP(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=EXP(A(I1))
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVEXP(A,I,C,J,N)
      DIMENSION A(1),C(1)
      REAL*8 TWOPI,ONOTPI
      DATA TWOPI,ONOTPI /6.28318530717959D0,0.159154943091895D0/
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
C      RR=A(I1)
      RR=A(I1)-INT(A(I1)*ONOTPI)*TWOPI
      C(J1)=COS(RR)
      C(J1+1)=SIN(RR)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVSQRT(A,I,C,J,N)
      COMPLEX CC
      REAL A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      CC=CMPLX(A(I1),A(I1+1))
      CC=CSQRT(CC)
      C(J1)=REAL(CC)
      C(J1+1)=AIMAG(CC)
      I1=I1+I
      J1=J1+J
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVMOV(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=A(I1)
      C(J1+1)=A(I1+1)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVNEG(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=-A(I1)
      C(J1+1)=-A(I1+1)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVCONJ(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=A(I1)
      C(J1+1)=-A(I1+1)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVMAGS(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=A(I1)*A(I1)+A(I1+1)*A(I1+1)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVMAX(A,I,C,N)
      DIMENSION A(1)
      I1=1
      C=0
      DO 10 L=1,N
      C=MAX(A(I1)*A(I1)+A(I1+1)*A(I1+1),C)
      I1=I1+I
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VCLR(C,J,N)
      DIMENSION C(1)
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=0.
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VFILL(A,C,J,N)
      DIMENSION C(1)
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=A
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VRAMP(A,B,C,J,N)
      DIMENSION C(1)
      J1=1
      AA=A
      BB=B
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=AA+BB*(L-1)
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVFILL(A,C,J,N)
      DIMENSION A(2),C(1)
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=A(1)
      C(J1+1)=A(2)
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVSMUL(A,I,C,D,J,N)
      DIMENSION A(1),D(1)
      I1=1
      J1=1
      CC=C
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      D(J1)=CC*A(I1)
      D(J1+1)=CC*A(I1+1)
      I1=I1+I
      J1=J1+J
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VMUL(A,I,C,J,D,K,N)
      DIMENSION A(1),C(1),D(1)
      I1=1
      J1=1
      K1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      D(K1)=C(J1)*A(I1)
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
      SUBROUTINE VADD(A,I,C,J,D,K,N)
      DIMENSION A(1),C(1),D(1)
      I1=1
      J1=1
      K1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      D(K1)=C(J1)+A(I1)
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVMEXP(A,I,C,J,D,K,N)
      DIMENSION A(1),C(1),D(1)
      REAL*8 TWOPI,ONOTPI
      DATA TWOPI,ONOTPI /6.28318530717959D0,0.159154943091895D0/
      I1=1
      J1=1
      K1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      RA=A(I1)-INT(A(I1)*ONOTPI)*TWOPI
      RC=C(J1)
      D(K1)=RC*COS(RA)
      D(K1+1)=RC*SIN(RA)
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVMUL(A,I,C,J,D,K,N,IFL)
      DIMENSION A(1),C(1),D(1)
      I1=1
      J1=1
      K1=1
      IF (IFL.GE.0) THEN
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      RR=A(I1)*C(J1)-A(I1+1)*C(J1+1)
      RI=A(I1+1)*C(J1)+A(I1)*C(J1+1)
      D(K1)=RR
      D(K1+1)=RI
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      ELSE
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 20 L=1,N
      RR=A(I1)*C(J1)+A(I1+1)*C(J1+1)
      RI=-A(I1+1)*C(J1)+A(I1)*C(J1+1)
      D(K1)=RR
      D(K1+1)=RI
      I1=I1+I
      J1=J1+J
      K1=K1+K
 20   CONTINUE
      END IF
      RETURN
      END
C
C
      SUBROUTINE VCLIP(A,I,B,C,D,K,N)
      DIMENSION A(1),D(1)
      I1=1
      K1=1
      DO 10 L=1,N
      IF (A(I1).LT.B) THEN
      D(K1)=B
      ELSE IF (A(I1).GT.C) THEN
           D(K1)=C
      ELSE
      D(K1)=A(I1)
      END IF
      I1=I1+I
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVSUB(A,I,C,J,D,K,N)
      DIMENSION A(1),C(1),D(1)
      I1=1
      J1=1
      K1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      D(K1)=A(I1)-C(J1)
      D(K1+1)=A(I1+1)-C(J1+1)
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVADD(A,I,C,J,D,K,N)
      DIMENSION A(1),C(1),D(1)
      I1=1
      J1=1
      K1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      D(K1)=A(I1)+C(J1)
      D(K1+1)=A(I1+1)+C(J1+1)
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CRVMUL(A,I,C,J,D,K,N)
      DIMENSION A(1),C(1),D(1)
      I1=1
      J1=1
      K1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      RR=C(J1)
      D(K1)=RR*A(I1)
      D(K1+1)=RR*A(I1+1)
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVRCIP(A,I,C,J,N)
      DIMENSION A(1),C(1)
C     REAL*8 RR,RR1,RR2
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      RR1=ABS(A(I1))+ABS(A(I1+1))
      RR1=1E0/RR1
      RR=A(I1)*RR1
      RI=A(I1+1)*RR1
      RS=1E0/(RR*RR+RI*RI)
      C(J1)=RR1*(RR*RS)
      C(J1+1)=RR1*(-RI*RS)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVIMVI(A,INDA,I,INDC,J,C,N)
      DIMENSION A(1),C(1)
      INTEGER INDA(1),INDC(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      II=INDA(I1)
      JJ=INDC(J1)
      C(JJ)=A(II)
      C(JJ+1)=A(II+1)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
      SUBROUTINE VSMA(A,I,B,C,J,D,K,N)
      DIMENSION A(1),C(1),D(1)
      I1=1
      J1=1
      K1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      D(K1)=B*A(I1)+C(J1)
      I1=I1+I
      J1=J1+J
      K1=K1+K
 10   CONTINUE
      RETURN
      END
C
C
      SUBROUTINE VALG10(A,I,C,J,N)
      DIMENSION A(1),C(1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      C(J1)=ALOG10(A(I1))
      I1=I1+I
 10   J1=J1+J
      RETURN
      END
      SUBROUTINE PKVAL(A,B,N,C,NP,R,MODE)
      DIMENSION A(1),B(1),C(1)
      RVAL=R
      RPEAK=1E0/R
      R=0.
      ITOG=MODE
      IF (ITOG.GT.0) THEN
      B(1)=0.
      ELSE
      B(1)=1E30
      END IF
      IP=1
      DO 10 I=1,N
      IF (ITOG.GT.0) THEN
        IF (A(I).GT.B(IP)) THEN
          B(IP)=A(I)
          C(IP)=I
          CHECK=RPEAK*B(IP)
        ELSE IF (A(I).LE.CHECK) THEN
          R=IP
          IP=IP+1
          ITOG=-ITOG
          IF (IP.GT.NP) RETURN
          B(IP)=A(I)
          C(IP)=I
        ELSE
        END IF
      ELSE
        IF (A(I).LT.B(IP)) THEN
          B(IP)=A(I)
          C(IP)=I
          CHECK=RVAL*B(IP)
        ELSE IF (A(I).GE.CHECK) THEN
          R=IP
          IP=IP+1
          ITOG=-ITOG
          IF (IP.GT.NP) RETURN
          B(IP)=A(I)
          C(IP)=I
        ELSE
        END IF
      END IF
 10   CONTINUE
      RETURN
      END
      SUBROUTINE CVIMOV(A,INDA,I,C,J,N)
      DIMENSION A(1),C(1)
      INTEGER INDA(1)
cvd$r relation(I.gt.0)
cvd$r relation(J.gt.1)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      II=INDA(I1)
      C(J1)=A(II)
      C(J1+1)=A(II+1)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CVMOVI(A,I,INDC,J,C,N)
      DIMENSION A(1),C(1)
      INTEGER INDC(1)
cvd$r relation(I.gt.1)
cvd$r relation(J.gt.0)
cvd$r permutation(INDC)
      I1=1
      J1=1
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 L=1,N
      JJ=INDC(J1)
      C(JJ)=A(I1)
      C(JJ+1)=A(I1+1)
      I1=I1+I
      J1=J1+J
  10  CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CBGEBS(A,AR,NR,N,IBW,EPS)     
      COMPLEX V,XC,CC1,A(1),AR(NR)      
      REAL X,ANORM
C     SPECIAL SOLVER MEANT FOR FIP SYSTEM       
C     SOLVES A SET OF LINEAR EQUATIONS, WHERE THE       
C     COEFFICIENT MATRIX IS OF BAND FORM WITH THE       
C     BAND-WIDTH IBW (NUMBER OF SUB- OR SUPER DIAGONALS)
C     BY GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING.    
C       
C     THE COEFFICIENT MATRIX SHOULD BE STORED IN BAND FORM  
C     FOR THIS SOLVER!!!
C *** REVISED 14.1.83 TO WORK WITH ONE-DIMENSIONAL ARRAYS
C *** CLEANED UP 8.1.89 TO IMPROVE ALLIANT OPTIMIZATION
C
C *** STATEMENT FUNCTION TO CALCULATE INDEX IN BANDED FORM
      IB(IR,IC) = (IBW+IC-IR)*NR+IR
cvd$r relation(NR.gt.0)
cvd$r relation(IBW.gt.0)
      EPPS=EPS
      EPS=0
      IBM=IBW   
      NRM1=NR-1
      DO 91 I=1,N       
C     --------------------------------- 
C     SEARCH FOR PIVOT ROW      
C     --------------------------------- 

      L1=MIN0(N,I+IBW)  
      ANORM=-1 
      II=IB(I,I)
      JI=II
cvd$  shortloop
      DO 30 J=I,L1      
       V=A(JI)
       JI=JI-NRM1
       X=ABS(REAL(V))+ABS(AIMAG(V))
       IF (X.GT.ANORM) THEN
        ANORM=X   
        K=J       
       END IF
 30   CONTINUE  
      IF (K.NE.I) THEN
C     -------------------------------   
C     CALCULATE NEW BAND-WIDTH  
C     -------------------------------   
       IBM=MAX0(IBM,K-I+IBW)     
C     ------------------------------    
C     INTERCHANGE ROWS I AND K  
C     ------------------------------    
       L2=MIN0(N,I+IBM)  
       IJ=II
       KJ=IB(K,I)
C$DIR NO_RECURRENCE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
cvd$  shortloop
       DO 50 J=I,L2      
        XC=A(IJ) 
        A(IJ)=A(KJ)     
        A(KJ)=XC 
        IJ=IJ+NR
        KJ=KJ+NR
 50    CONTINUE  
       XC=AR(I) 
       AR(I)=AR(K)     
       AR(K)=XC 
      END IF
C     -------------------------------   
C     TEST FOR SINGULARITY      
C     -------------------------------   
      IF (ANORM.LT.EPPS) THEN  
C     --------------------------------  
C     MATRIX IS SINGULAR
C     --------------------------------  
       EPS=1
       RETURN
      END IF  
C     --------------------------------  
C     REDUCTION OF REMAINING ROWS       
C     --------------------------------  
      V=1E0/A(II)
      A(II)=V
      L2=MIN(I+IBM,N)
cvd$  shortloop
C$DIR NO_RECURRENCE
cvd$  nodepchk
cvd$  select(concurrent)
      DO 85 K=I+1,L1    
       KI=IB(K,I)
       CC1=-A(KI)
       RR=REAL(CC1)
       RI=AIMAG(CC1)
       IF ((RR.NE.0.0) .OR. (RI.NE.0.0)) THEN
        XC=CC1*V
        A(KI)=CMPLX(0.,0.)    
        IJ=II
        KJ=KI
cvd$  shortloop
C$DIR NO_RECURRENCE
cvd$  nodepchk
cvd$  select(vector)
        DO 80 J=I+1,L2    
         IJ=IJ+NR
         KJ=KJ+NR
         A(KJ)=A(KJ)+XC*A(IJ)
 80     CONTINUE
        AR(K)=AR(K)+XC*AR(I)
       END IF
 85   CONTINUE  
 91   CONTINUE  
C     -------------------------------   
C     BACK SUBSTITUTION 
C     -------------------------------   
      AR(N)=AR(N)*A(IB(N,N))      
      DO 100 I=N-1,1,-1 
       L2=MIN(N,I+IBM)
       II=IB(I,I)
       IJ=II
       V=(0.,0.)
cvd$  shortloop
cvd$  select(vector,concurrent)
       DO 95 J=I+1,L2    
        IJ=IJ+NR
        V=V+A(IJ)*AR(J)
 95    CONTINUE
       AR(I)=(AR(I)-V)*A(II)  
 100  CONTINUE
      RETURN    
      END       

