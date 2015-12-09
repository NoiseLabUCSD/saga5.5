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
c     SEARCH FOR PIVOT ROW      
C     --------------------------------- 
c
       L1=MIN0(N,I+IBW)  
c      ANORM=-1 
      II=IB(I,I)
      JI=II
cvd$  shortloop
c      DO 30 J=I,L1      
c       V=A(JI)
c       JI=JI-NRM1
c       X=ABS(REAL(V))+ABS(AIMAG(V))
c       IF (X.GT.ANORM) THEN
c        ANORM=X   
c        K=J       
c       END IF
c 30   CONTINUE  
c      IF (K.NE.I) THEN
C     -------------------------------   
C     CALCULATE NEW BAND-WIDTH  
C     -------------------------------   
c       IBM=MAX0(IBM,K-I+IBW)     
C     ------------------------------    
C     INTERCHANGE ROWS I AND K  
C     ------------------------------    
c       L2=MIN0(N,I+IBM)  
c       IJ=II
c       KJ=IB(K,I)
ccvd$  nodepchk
ccvd$  shortloop
c       DO 50 J=I,L2      
c        XC=A(IJ) 
c        A(IJ)=A(KJ)     
c        A(KJ)=XC 
c        IJ=IJ+NR
c        KJ=KJ+NR
c 50    CONTINUE  
c       DO 51 J=1,NRHS
c        XC=AR(I,J) 
c        AR(I,J)=AR(K,J)     
c        AR(K,J)=XC
c 51    CONTINUE 
c      END IF
cC     -------------------------------   
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
