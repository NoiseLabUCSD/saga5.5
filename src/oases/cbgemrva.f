      SUBROUTINE CBGEMR(A,AR,NR,N,NRHS,IBW,EPS)     
C     SPECIAL SOLVER MEANT FOR SAFARI PACKAGE       
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
      COMPLEX V,XC,CC1,A(1),AR(NR,NRHS)      
      REAL X,ANORM
C *** DYNAMIC BUFFER ARRAY ONLY FOR ALLIANT
C      COMPLEX VR(NRHS)
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
cvd$  nodepchk
cvd$  shortloop
       DO 50 J=I,L2      
        XC=A(IJ) 
        A(IJ)=A(KJ)     
        A(KJ)=XC 
        IJ=IJ+NR
        KJ=KJ+NR
 50    CONTINUE  
       DO 51 J=1,NRHS
        XC=AR(I,J) 
        AR(I,J)=AR(K,J)     
        AR(K,J)=XC
 51    CONTINUE 
      END IF
C     -------------------------------   
C     TEST FOR SINGULARITY      
C     -------------------------------   
c      IF (ANORM.LT.EPPS) THEN  
C     --------------------------------  
C     MATRIX IS SINGULAR
C     --------------------------------  
c       EPS=1
c       RETURN
c      END IF  
C     --------------------------------  
C     REDUCTION OF REMAINING ROWS       
C     --------------------------------  
      V=1E0/A(II)
      A(II)=V
      L2=MIN(I+IBM,N)
cvd$  shortloop
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
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
cvd$  select(vector)
        DO 80 J=I+1,L2    
         IJ=IJ+NR
         KJ=KJ+NR
         A(KJ)=A(KJ)+XC*A(IJ)
 80     CONTINUE
        DO 81 J=1,NRHS
         AR(K,J)=AR(K,J)+XC*AR(I,J)
 81     CONTINUE
       END IF
 85   CONTINUE  
 91   CONTINUE  
C     -------------------------------   
C     BACK SUBSTITUTION 
C     -------------------------------   
      DO 100 K=1,NRHS
c      DO 94 K=1,NRHS
      AR(N,K)=AR(N,K)*A(IB(N,N))      
c 94   CONTINUE
      DO 100 I=N-1,1,-1 
       L2=MIN(N,I+IBM)
       II=IB(I,I)
cvd$  select(concurrent)
c       DO 98 K=1,NRHS
c        VR(K)=(0.,0.)
        V=(0.,0.)
        IJ=II
cvd$  shortloop
cvd$  select(vector)
        DO 95 J=I+1,L2    
         IJ=IJ+NR
c         VR(K)=VR(K)+A(IJ)*AR(J,K)
         V=V+A(IJ)*AR(J,K)
 95     CONTINUE
c        AR(I,K)=(AR(I,K)-VR(K))*A(II)  
        AR(I,K)=(AR(I,K)-V)*A(II)  
c 98    CONTINUE
 100  CONTINUE
      RETURN    
      END       
