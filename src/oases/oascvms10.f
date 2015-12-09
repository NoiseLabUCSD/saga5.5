      COMPLEX FUNCTION FILON(CFUN,CEX,ARG,N,FAC)
      COMPLEX AI,CSUM,CFUN(1),CEX(1),DFUN,DEX,DARG,DPROD
      DIMENSION ARG(1)
      AI=CMPLX(0.,1E0)
      CSUM=0E0
      DO 10 I=1,N-1
      DARG=-AI/(ARG(I+1)-ARG(I))
      DFUN=CFUN(I+1)-CFUN(I)
      DEX=CEX(I+1)-CEX(I)
      DPROD=CFUN(I+1)*CEX(I+1)-CFUN(I)*CEX(I)
      CSUM=CSUM+DARG*(DPROD-DARG*DFUN*DEX)
  10  CONTINUE
      FILON=FAC*CSUM
      RETURN
      END
      SUBROUTINE CMMUL(A,NRA,NCA,B,NCB,C)
      COMPLEX A(NRA,NCA),B(NCA,NCB),C(NRA,NCB)
      COMPLEX SUM
      DO 30 I=1,NRA
      DO 20 J=1,NCB
      SUM=CMPLX(0E0,0E0)
      DO 10 K=1,NCA
      SUM=SUM+A(I,K)*B(K,J)
  10  CONTINUE
      C(I,J)=SUM
  20  CONTINUE
  30  CONTINUE
      RETURN
      END        

      SUBROUTINE CTMMUL(A,NRA,NCA,B,NCB,C)
c    
c     Subroutine for calculating the matrix product
c     of the complex transpose of A(NRA,NCA) and
c     the matrix B(NRA,NCB). The result is placed in
c     matrix C(NCA,NCB).
c
      COMPLEX A(NRA,NCA),B(NRA,NCB),C(NCA,NCB)
      COMPLEX SUM
      DO 30 I=1,NCA
      DO 20 J=1,NCB
      SUM=CMPLX(0E0,0E0)
      DO 10 K=1,NRA
      SUM=SUM+CONJG(A(K,I))*B(K,J)
  10  CONTINUE
      C(I,J)=SUM
  20  CONTINUE
  30  CONTINUE
      RETURN
      END        

        SUBROUTINE CMATIN(NM,A,AINV)
        PARAMETER (NMAX=100)
        COMPLEX A(NM,NM),AINV(NM,NM)
C        DIMENSION AR(NMAX,NMAX),AI(NMAX,NMAX)
        DIMENSION BR(NMAX),BI(NMAX),CR(NMAX),
     1            CI(NMAX),IP(NMAX),IQ(NMAX)
C       DOUBLE PRECISION AR,AI,CR,CI,BR,BI,ZR,ZI,PR,PI
C    
C      REAL AND IMAGINARY PARTS OF MATRIX ADDRESSED BY 
C      STATEMENT FUNCTIONS
C
       AR(II,JJ)=REAL(AINV(II,JJ))
       AI(II,JJ)=AIMAG(AINV(II,JJ))
       IF (NM.GT.NMAX) THEN
         STOP '*** ORDER OF MATRIX TOO HIGH IN CMATIN ***'
       END IF
       DO 51 I=1,NM
        DO 51 J=1,NM
 51     AINV(I,J)=A(I,J)
       DO 1 K=1,NM
       PR=0.0D0
       PI=PR
       DO 100 I=K,NM
       DO 100 J=K,NM
       ZR=AR(I,J)
       ZI=AI(I,J)
       IF(ZR*ZR+ZI*ZI-PR*PR-PI*PI)100,100,101
 101    PR=AR(I,J)
       PI=AI(I,J)
        IP(K)=I
       IQ(K)=J
 100   CONTINUE
       IF(ABS(PR)+ABS(PI).LT.1.0D-37)GO TO 35
       IPK=IP(K)
       IQK=IQ(K)
       IF(IPK-K)200,299,200
 200   DO 201 J=1,NM
       ZR=AR(IPK,J)
       ZI=AI(IPK,J)
       AINV(IPK,J)=AINV(K,J)
       AINV(K,J)=CMPLX(ZR,ZI)
 201   CONTINUE
 299   CONTINUE
       IF(IQK-K)300,399,300
 300   DO 301 I=1,NM
       ZR=AR(I,IQK)
       ZI=AI(I,IQK)
       AINV(I,IQK)=AINV(I,K)
       AINV(I,K)=CMPLX(ZR,ZI)
 301   CONTINUE
 399   CONTINUE
        ZR=1.0D0/(PR*PR+PI*PI)
       DO 400 J=1,NM
       IF(J-K)403,402,403
 402   BR(J)=PR*ZR
       BI(J)=-PI*ZR
       CR(J)=1.0D0
       CI(J)=0.0D0
       GO TO 404
 403   BR(J)=-(AR(K,J)*PR+AI(K,J)*PI)*ZR
       BI(J)=-(AI(K,J)*PR-AR(K,J)*PI)*ZR
       CR(J)=AR(J,K)
       CI(J)=AI(J,K)
 404   AINV(K,J)=CMPLX(0E0,0E0)
       AINV(J,K)=CMPLX(0E0,0E0)
  400  CONTINUE
       DO 405 I=1,NM
       DO 405 J=1,NM
       AINV(I,J)=AINV(I,J)+CMPLX(CR(I)*BR(J)-CI(I)*BI(J),
     1       CR(I)*BI(J)+CI(I)*BR(J))
 405   CONTINUE
  1    CONTINUE
       K=NM
       DO 500 KM=1,NM
       IPK=IP(K)
       IQK=IQ(K)
       IF(IPK-K)501,502,501
 501   DO 503 I=1,NM
       ZR=AR(I,IPK)
       ZI=AI(I,IPK)
       AINV(I,IPK)=AINV(I,K)
       AINV(I,K)=CMPLX(ZR,ZI)
 503   CONTINUE
 502   IF(IQK-K)504,500,504
 504   DO 506 J=1,NM
       ZR=AR(IQK,J)
       ZI=AI(IQK,J)
       AINV(IQK,J)=AINV(K,J)
       AINV(K,J)=CMPLX(ZR,ZI)
 506   CONTINUE
 500   K=K-1
          RETURN
  35      PRINT 37
  37      FORMAT(' ***** MATRIX SINGULAR PROGRAM STOPS ***** ')
           STOP
          END
        SUBROUTINE CMINV(NM,A,AINV,IERR)
        PARAMETER (NMAX=150)
        COMPLEX A(NM,NM),AINV(NM,NM)
C        DIMENSION AR(NMAX,NMAX),AI(NMAX,NMAX)
        DIMENSION BR(NMAX),BI(NMAX),CR(NMAX),
     1            CI(NMAX),IP(NMAX),IQ(NMAX)
C       DOUBLE PRECISION AR,AI,CR,CI,BR,BI,ZR,ZI,PR,PI
C    
C      REAL AND IMAGINARY PARTS OF MATRIX ADDRESSED BY 
C      STATEMENT FUNCTIONS
C
       AR(II,JJ)=REAL(AINV(II,JJ))
       AI(II,JJ)=AIMAG(AINV(II,JJ))

       IERR=0
       IF (NM.GT.NMAX) THEN
         STOP '*** ORDER OF MATRIX TOO HIGH IN CMINV ***'
       END IF
       DO 51 I=1,NM
        DO 51 J=1,NM
 51     AINV(I,J)=A(I,J)
       DO 1 K=1,NM
       PR=0.0D0
       PI=PR
       DO 100 I=K,NM
       DO 100 J=K,NM
       ZR=AR(I,J)
       ZI=AI(I,J)
       IF(ZR*ZR+ZI*ZI-PR*PR-PI*PI)100,100,101
 101    PR=AR(I,J)
       PI=AI(I,J)
        IP(K)=I
       IQ(K)=J
 100   CONTINUE
       IF(ABS(PR)+ABS(PI).LT.1.0D-37)GO TO 35
       IPK=IP(K)
       IQK=IQ(K)
       IF(IPK-K)200,299,200
 200   DO 201 J=1,NM
       ZR=AR(IPK,J)
       ZI=AI(IPK,J)
       AINV(IPK,J)=AINV(K,J)
       AINV(K,J)=CMPLX(ZR,ZI)
 201   CONTINUE
 299   CONTINUE
       IF(IQK-K)300,399,300
 300   DO 301 I=1,NM
       ZR=AR(I,IQK)
       ZI=AI(I,IQK)
       AINV(I,IQK)=AINV(I,K)
       AINV(I,K)=CMPLX(ZR,ZI)
 301   CONTINUE
 399   CONTINUE
        ZR=1.0D0/(PR*PR+PI*PI)
       DO 400 J=1,NM
       IF(J-K)403,402,403
 402   BR(J)=PR*ZR
       BI(J)=-PI*ZR
       CR(J)=1.0D0
       CI(J)=0.0D0
       GO TO 404
 403   BR(J)=-(AR(K,J)*PR+AI(K,J)*PI)*ZR
       BI(J)=-(AI(K,J)*PR-AR(K,J)*PI)*ZR
       CR(J)=AR(J,K)
       CI(J)=AI(J,K)
 404   AINV(K,J)=CMPLX(0E0,0E0)
       AINV(J,K)=CMPLX(0E0,0E0)
  400  CONTINUE
       DO 405 I=1,NM
       DO 405 J=1,NM
       AINV(I,J)=AINV(I,J)+CMPLX(CR(I)*BR(J)-CI(I)*BI(J),
     1                           CR(I)*BI(J)+CI(I)*BR(J))
 405   CONTINUE
  1    CONTINUE
       K=NM
       DO 500 KM=1,NM
       IPK=IP(K)
       IQK=IQ(K)
       IF(IPK-K)501,502,501
 501   DO 503 I=1,NM
       ZR=AR(I,IPK)
       ZI=AI(I,IPK)
       AINV(I,IPK)=AINV(I,K)
       AINV(I,K)=CMPLX(ZR,ZI)
 503   CONTINUE
 502   IF(IQK-K)504,500,504
 504   DO 506 J=1,NM
       ZR=AR(IQK,J)
       ZI=AI(IQK,J)
       AINV(IQK,J)=AINV(K,J)
       AINV(K,J)=CMPLX(ZR,ZI)
 506   CONTINUE
 500   K=K-1
          RETURN
  35      IERR=1
          RETURN
          END
      subroutine vdecim(a,i,b,j,nold,ndec,nnew)
      parameter (pi=3.14159)
c
c     subroutine for decimating a vector by running averaging
c     Output vector will be of dimension 1+(nold-1)/ndec 
C
      dimension a(1),b(1)
      dimension w(51)
      np=2*ndec + 1
      nnew=1+(nold-1)/ndec
      if (nnew.gt.nold) then
        write(6,*) '>>>> ERROR in VDECIM: nnew > nold <<<<'
        stop
      else if (np.gt.51) then
        stop '>>>> ERROR in VDECIM: Cannot decimate by more than 25'
      end if
c
c     determine window parameters
c
      faca=pi/(np-1)
      wsum=0E0
      do 10 l=1,np          
       w(l)=sin((l-1)*faca)**2
       wsum=wsum+w(l)
 10   continue
      fac=1E0/wsum 
      INEW=I*NDEC
      i1=1
      j1=1
      b(1)=a(1)
      i1=i1+INEW
      j1=j1+j
      do 20 l=2,nnew-1
        sum=0
        iof=-ndec*i
        do 15 k=-ndec,ndec
          sum=sum+w(ndec+1+k)*a(i1+iof)
          iof=iof+i
 15     continue
        b(j1)=fac*sum
        i1=i1+INEW
        j1=j1+j
 20   continue
      b(j1)=A(I1)
      return
      end
      COMPLEX FUNCTION CINTPL(C,XMIN,ONODX,N,X)
C *** INTERPOLATES IN A COMPLEX ARRAY C TO DETERMINE VALUE
C *** FOR ARGUMENT X
      COMPLEX C(N)
      RINDEX=(X-XMIN)*ONODX+1
      INDEX=IFIX(RINDEX)
      IF (INDEX.GE.N) THEN
        WRITE(6,*) '>>> WARNING: ARRAY OVERFLOW IN CINTPL, ARG=',X
      END IF
      INDEX=MIN(INDEX,N-1)
      INDEX=MAX(INDEX,1) 
      REM=RINDEX-INDEX
      CINTPL=(1E0-REM)*C(INDEX) + REM*C(INDEX+1)
      RETURN
      END
      REAL FUNCTION RINTPL(C,XMIN,ONODX,N,X)
C *** INTERPOLATES IN A REAL ARRAY C TO DETERMINE VALUE
C *** FOR ARGUMENT X
      REAL C(N)
      RINDEX=(X-XMIN)*ONODX+1
      INDEX=IFIX(RINDEX)
      IF (INDEX.GE.N) THEN
        WRITE(6,*) '>>> WARNING: ARRAY OVERFLOW IN RINTPL, ARG=',X
      END IF
      INDEX=MIN(INDEX,N-1)
      INDEX=MAX(INDEX,1) 
      REM=RINDEX-INDEX
      RINTPL=(1E0-REM)*C(INDEX) + REM*C(INDEX+1)
      RETURN
      END
      SUBROUTINE PREPBF(RMAX,WKMAX)
      PARAMETER (NBDIV=40)
      INCLUDE 'compar.f'
      INCLUDE 'combes.f'
      RKMAX=RMAX*WKMAX
      NRKMAX=MAX(NINT(NBDIV*RKMAX/(2*PI)),2)
      IF (NRKMAX.GT.NP) THEN
       WRITE(6,*) '>>> WARNING: NUMBER OF BESSEL INTERPOLATION'
       WRITE(6,*) '>>>          POINTS TRUNCATED TO'
       WRITE(6,*) '>>>         ', NBDIV*FLOAT(NP)/NRKMAX,'PER LAMBDA'
       NRKMAX=NP
      END IF
      DRK=MAX(1E-6,RKMAX)/(NRKMAX-1.001)
      ONODRK=1E0/DRK
      DO 10 II=1,NRKMAX
       RARG=DRK*(II-1)
       BF0(II)=BESSJ0(RARG)
       BF1(II)=BESSJ1(RARG)
C       CALL BESJ01(RARG,BF0(II),BF1(II))
 10   CONTINUE
      RETURN
      END
      SUBROUTINE BESJ01(ARG,B0,B1)
      REAL ARG,B0,B1
      B0=BESSJ0(ARG)
      B1=BESSJ1(ARG)
      RETURN
      END
      REAL FUNCTION bessj0(x)
      REAL X,AX,Z
      REAL*8 XX,Y,ANS,ANS1,ANS2
      AX=ABS(X)
      if (ax.LT. 8.0) THEN
       y=x*x
       ans1=57568490574.0D0+y*(-13362590354.0D0+y*(651619640.7D0
     &      +y*(-11214424.18D0+y*(77392.33017D0+y*(-184.9052456D0)))))
       ans2=57568490411.0D0+y*(1029532985.0D0+y*(9494680.718D0
     &      +y*(59272.64853D0+y*(267.8532712D0+y*1.0D0))))
       ans=ans1/ans2
      else 
       z=8.0/ax
       y=z*z
       xx=ax-0.785398164
       ans1=1.0D0+y*(-0.1098628627D-2+y*(0.2734510407D-4
     &      +y*(-0.2073370639D-5+y*0.2093887211D-6)))
       ans2 = -0.1562499995D-1+y*(0.1430488765D-3
     &        +y*(-0.6911147651D-5+y*(0.7621095161D-6
     &        -y*0.934935152D-7)))
       ans=sqrt(0.636619772D0/ax)*(cos(xx)*ans1-z*sin(xx)*ans2)
      END IF
      BESSJ0=ANS
      RETURN
      END
      REAL FUNCTION bessj1(x)
      REAL X,ax,z
      REAL*8 xx,y,ans,ans1,ans2
      AX=ABS(X)
      if (ax.LT. 8.0) THEN
       y=x*x
       ans1=x*(72362614232.0D0+y*(-7895059235.0D0+y*(242396853.1D0
     &      +y*(-2972611.439D0+y*(15704.48260D0+y*(-30.16036606D0))))))
       ans2=144725228442.0D0+y*(2300535178.0D0+y*(18583304.74D0
     &      +y*(99447.43394D0+y*(376.9991397D0+y*1.0D0))))
       ans=ans1/ans2
      else 
       z=8.0/ax
       y=z*z
       xx=ax-2.356194491D0
       ans1=1.0D0+y*(0.183105D-2+y*(-0.3516396496D-4
     &      +y*(0.2457520174D-5+y*(-0.240337019D-6))))
       ans2=0.04687499995D0+y*(-0.2002690873D-3
     &      +y*(0.8449199096D-5+y*(-0.88228987D-6
     &      +y*0.105787412D-6)))
       ans=sqrt(0.636619772D0/ax)*(cos(xx)*ans1-z*sin(xx)*ans2)
       if (x .LT. 0.0) ans = -ans
      END IF
      BESSJ1=ANS
      RETURN
      END


