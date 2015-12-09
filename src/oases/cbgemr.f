      SUBROUTINE CBGEMR(A,AR,NR,N,NRHS,IBW,EPS)     
c     special solver meant for safari package       
c     solves a set of linear equations, where the       
c     coefficient matrix is of band form with the       
c     band-width ibw (number of sub- or super diagonals)
c     by gaussian elimination with partial pivoting.    
c       
c     the coefficient matrix should be stored in band form  
c     for this solver!!!
c *** revised 14.1.83 to work with one-dimensional arrays
c *** cleaned up 8.1.89 to improve alliant optimization
c
      COMPLEX V,XC,CC1,A(1),AR(NR,NRHS)      
      REAL X,ANORM
c *** dynamic buffer array only for alliant
      COMPLEX VR(NRHS)
c *** statement function to calculate index in banded form
      IB(IR,IC) = (IBW+IC-IR)*NR+IR
cvd$r relation(nr.gt.0)
cvd$r relation(ibw.gt.0)
      EPPS=EPS
      EPS=0
      IBM=IBW   
      NRM1=NR-1
      DO 91 I=1,N       
c     --------------------------------- 
c     search for pivot row      
c     --------------------------------- 

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
c     -------------------------------   
c     calculate new band-width  
c     -------------------------------   
       IBM=MAX0(IBM,K-I+IBW)     
c     ------------------------------    
c     interchange rows i and k  
c     ------------------------------    
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
c     -------------------------------   
c     test for singularity      
c     -------------------------------   
      IF (ANORM.LT.EPPS) THEN  
c     --------------------------------  
c     matrix is singular
c     --------------------------------  
       EPS=1
       RETURN
      END IF  
c     --------------------------------  
c     reduction of remaining rows       
c     --------------------------------  
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
c     -------------------------------   
c     back substitution 
c     -------------------------------   
      DO 94 K=1,NRHS
      AR(N,K)=AR(N,K)*A(IB(N,N))      
 94   CONTINUE
      DO 100 I=N-1,1,-1 
       L2=MIN(N,I+IBM)
       II=IB(I,I)
cvd$  select(concurrent)
       DO 98 K=1,NRHS
        VR(K)=(0.,0.)
        IJ=II
cvd$  shortloop
cvd$  select(vector)
        DO 95 J=I+1,L2    
         IJ=IJ+NR
         VR(K)=VR(K)+A(IJ)*AR(J,K)
 95     CONTINUE
        AR(I,K)=(AR(I,K)-VR(K))*A(II)  
 98    CONTINUE
 100  CONTINUE
      RETURN    
      END       
