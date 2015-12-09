      SUBROUTINE CALELC
C
C     CALCULATES ELASTIC PARAMETERS FOR TRANS. ISOTR. LAYERS
C
C     CHANGED TO INCLUDE ATTENUATION 860124
C
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comti.f'
C
      DUM=0E0
      DO 10 I=1,NTISOL
      LAYN=LAYNTI(I)
      SLNORM(LAYN)=1E10
      DO 5 J=1,NSL(LAYN)
 5    SLNORM(LAYN)=MIN(SLNORM(LAYN),REAL(BSP(J,LAYN)))
      CALL FINELAY(NSL(LAYN),ASP(1,LAYN),BSP(1,LAYN),ARO(1,LAYN),
     &             AH(1,LAYN))
      CALL VMOV(A,1,ELAG(1,LAYN),1,24)
      if (debug) then
       WRITE(26,*)'A,B1,B0,C2,C1,C0,C44,C66,C33,C13,C11,RH:'
       WRITE(26,*) (ELAG(J,LAYN),J=1,12)
      end if
 10   CONTINUE
      RETURN
      END
      SUBROUTINE CALSLN(SLOW)
C
C     CALCULATES VERTICAL SLOWNESSES, DISPL. AND STRESS COEFF. FOR
C     HORIZONTAL SLOWNESS SLOW
C
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comti.f'
      COMPLEX SLOW
C
      DO 10 I=1,NTISOL
      LAYN=LAYNTI(I)
      CALL VMOV(ELAG(1,LAYN),1,A,1,24)
      CALL TIEIG(SLOW)
      CALL WRBUF(41,S3UP,52)
 10   CONTINUE
      RETURN
      END
      SUBROUTINE FINELAY(NLAY,ALPHA,BETA,RHO,H)
C NLAY media are required (NLAY=2 or 3) consisisting of
C           ALPHA(I),BETA(I),RHO(I),H(I), I=1 to NLAY
C               where ALPHA(I) is compresional speed,
C                     BETA(I) is shear speed,
C                     RHO(I) is density,
C                     H(I) is volume fraction         of Ith layer.
C                          H(1)+H(2)+H(3)=1
C Outputs are the 5 moduli and the density of the long wavelength
C  equivalent transversely isotropic medium.
C  Note: Usual case has 2 constituents so H(2)=1-H(1)
C  CHANGED TO COVER LOSSY CASE (COMPLEX VEL.) 860124
C
C     IMPLICIT NONE
      INCLUDE 'compar.f'
      INCLUDE 'comti.f'	

      INTEGER I,NLAY
      COMPLEX ALPHA(3),BETA(3),CRHO(3),X
      REAL*4  RHO(3),H(3)
      COMPLEX XMU(3),GAM(3),XMTG(3),XMDG(3)

      DO 3 I=1,NLAY
 3    CRHO(I)=CMPLX(RHO(I),0E0)
      DO 5 I=1,NLAY
      XMU(I)=RHO(I)*BETA(I)**2
      GAM(I)=(BETA(I)/ALPHA(I))**2
      XMTG(I)=XMU(I)*GAM(I)
 5      XMDG(I)=XMU(I)/GAM(I)
      CALL AV(NLAY,XMU,H,-1,C44)
      CALL AV(NLAY,XMDG,H,-1,C33)
      CALL AV(NLAY,GAM,H,1,X)
      C13=(1.-2.*X)*C33
      CALL AV(NLAY,XMU,H,1,C66)
      CALL AV(NLAY,XMTG,H,1,X)
      C11=4.*(C66-X)+C13*C13/C33
      CALL AV(NLAY,CRHO,H,1,RH)
      A=C13+C44
      B1=2.*C13/C33+C13*C13/C33/C44-C11/C44
      B0=RH*(1./C44+1./C33)
      C2=C11/C33
      C1=-RH*(1+C11/C44)/C33
      C0=RH*RH/C33/C44
      RETURN
      END
C
      SUBROUTINE AV(N,X,H,J,XA)      
C     IMPLICIT NONE
      INTEGER I,J,N
      COMPLEX X(3),XA
      REAL*4 H(3),V
      XA=0.
      DO 1 I=1,N
 1      XA=XA+H(I)*X(I)**J
      V=1./J
C  SMALL IMAGINARY PART
      R=REAL(XA)
      ANG=AIMAG(XA)*V/REAL(XA)
      XA=(R**V)*EXP(CMPLX(0E0,ANG))
      RETURN
      END
C
      SUBROUTINE TIEIG(S1)
      INCLUDE 'compar.f'
      INCLUDE 'comti.f'
C     IMPLICIT NONE
      INTEGER J
      COMPLEX S1,SQ1
      COMPLEX SQ3(2),B,C,SCFAC(2)
C needed for eigenvalues & eigenvectors from SUBROUTINE FINELAY
      SQ1=S1**2
      B=(SQ1*B1+B0)/2.
      C=SQ1*SQ1*C2+SQ1*C1+C0
      C=B*B-C
      IF (REAL(C).GE.0) THEN
        SQ3(1)=B-CSQRT(C)
      ELSE
        SQ3(1)=B-AI*CSQRT(-C)
      END IF  
      SQ3(2)=2.*B-SQ3(1)
      IF(REAL(SQ3(1)).GE.0.)THEN
       S3DN(2)=SQRT(SQ3(2))
       S3UP(2)=-S3DN(2)
       S3DN(1)=SQRT(SQ3(1))
       S3UP(1)=-S3DN(1)
      END IF
C 2 real, 2 imaginary roots
      IF(REAL(SQ3(1)).LT.0.AND.REAL(SQ3(2)).GE.0.)THEN
       S3DN(2)=SQRT(SQ3(2))
       S3UP(2)=-S3DN(2)
       S3DN(1)=-AI*SQRT(-SQ3(1))
       S3UP(1)=-S3DN(1)
      END IF
C 4 imaginary roots
      IF(REAL(SQ3(2)).LT.0.)THEN
       S3DN(2)=-AI*SQRT(-SQ3(2))
       S3UP(2)=-S3DN(2)
       S3DN(1)=-AI*SQRT(-SQ3(1))
       S3UP(1)=-S3DN(1)
      END IF
C To find eigenvectors, stress (here called SIG) is stress/-i omega
      IF(ABS(SQ3(1)).NE.0.)GO TO 20
       U3DN(1)=0.
       U1DN(1)=A/ABS(A)
       U3UP(1)=U3DN(1)
       U1UP(1)=U1DN(1)
      GO TO 21
 20     U3DN(1)=(RH-C11*SQ1-C44*SQ3(1))
      U1DN(1)=A*S1*S3DN(1)
      U3UP(1)=-U3DN(1)
      U1UP(1)=-A*S1*S3UP(1)
 21     IF(ABS(SQ3(2)).NE.0.)GO TO 30
       U1DN(2)=0.
       U3DN(2)=-A/ABS(A)
       U1UP(2)=U1DN(2)
       U3UP(2)=U3DN(2)
      GO TO 31
 30     U1DN(2)=-(RH-C44*SQ1-C33*SQ3(2))
      U3DN(2)=-A*S1*S3DN(2)
      U1UP(2)=-U1DN(2)
      U3UP(2)=A*S1*S3UP(2)
 31    DO 2 J=1,2
      SIG3DN(J)=C13*S1*U1DN(J)+C33*S3DN(J)*U3DN(J)
      SIG13DN(J)=C44*(S1*U3DN(J)+S3DN(J)*U1DN(J))
      SIG3UP(J)=C13*S1*U1UP(J)+C33*S3UP(J)*U3UP(J)
 2      SIG13UP(J)=C44*(S1*U3UP(J)+S3UP(J)*U1UP(J))
C
C     SCALING FOR SAFARI COMPATIBILITY
C
      U3DN(1)=-U3DN(1)
      U3UP(1)=-U3UP(1)
      U1DN(1)=AI*U1DN(1)
      U1UP(1)=AI*U1UP(1)
      SIG3DN(1)=-SIG3DN(1)
      SIG3UP(1)=-SIG3UP(1)
      SIG13DN(1)=AI*SIG13DN(1)
      SIG13UP(1)=AI*SIG13UP(1)
      U1DN(2)=-AI*U1DN(2)
      U1UP(2)=-AI*U1UP(2)
      SIG13DN(2)=-AI*SIG13DN(2)
      SIG13UP(2)=-AI*SIG13UP(2)
      SCFAC(1)=-AI/(A*S3DN(1))
      SCFAC(2)=1E0/(A*S3DN(2))
      DO 10 I=1,2
      U1DN(I)=U1DN(I)*SCFAC(I)
      U1UP(I)=U1UP(I)*SCFAC(I)
      U3DN(I)=U3DN(I)*SCFAC(I)
      U3UP(I)=U3UP(I)*SCFAC(I)
      SIG13DN(I)=SIG13DN(I)*SCFAC(I)
      SIG13UP(I)=SIG13UP(I)*SCFAC(I)
      SIG3DN(I)=SIG3DN(I)*SCFAC(I)
      SIG3UP(I)=SIG3UP(I)*SCFAC(I)
      S3DN(I)=AI*S3DN(I)
      S3UP(I)=AI*S3UP(I)
 10   CONTINUE
C
C    SH WAVES ADDED 880128. HS
C
      S3SHUP=CSQRT((C66*SQ1-RH)/C44)
      S3SHDN=-S3SHUP
      U2UP=S1
      U2DN=S1
      SIG23DN=S3SHDN*S1*C44
      SIG23UP=S3SHUP*S1*C44
      RETURN
      END
C
      SUBROUTINE GETSLN(IL,NSLOW)
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comti.f'

      DIMENSION X(NP2,3)
      equivalence (X(1,1),CFF(1,1))

      CALL RWDBUF(41)
      DO 10 ISL=1,NSLOW
      DO 5 I=1,NTISOL
      LAYN=LAYNTI(I)
      CALL RDBUF(41,S3UP,52)
      IF (LAYN.EQ.IL) THEN
      X(ISL,1)=ABS(AIMAG(S3UP(1)))*SLNORM(IL)
      X(ISL,2)=ABS(AIMAG(S3UP(2)))*SLNORM(IL)
      X(ISL,3)=ABS(AIMAG(S3SHUP))*SLNORM(IL)
      END IF
 5    CONTINUE
 10   CONTINUE
      RETURN
      END
      SUBROUTINE PLSLOW(NSLOW,TITLE,
     1     XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      YUP,YDOWN,YINC)
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comti.f'
      INCLUDE 'complo.f'
      DIMENSION X(NP2,3)
      equivalence (X(1,1),CFF(1,1))
      DIMENSION XS(NP2)
      EQUIVALENCE (XS(1),CFFS(1))
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      OPTION(1)=PROGNM
      OPTION(2)='SLOWNS'
      DO 10 I=1,NTISOL
      IL=LAYNTI(I)
      CALL GETSLN(IL,NSLOW)
      PTIT='VERTICAL SLOWNESS'
 811  FORMAT('Layer:',I4,'$')
 812  FORMAT('Cr: ',F6.1,' m/s$')
      NLAB=2
      WRITE(LAB(1),811) IL
      WRITE(LAB(2),812) SLNORM(IL)
      XTXT='Norm. horizontal slowness$'
      YTXT='Norm. vertical slowness$'
      XTYP='LIN'
      YTYP='LIN'
      XDIV=1
      YDIV=1
      IGRID=1
      NC=3
C *** WRITE PLP FILE
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
      CALL VSMUL(XS,1,SLNORM(IL),ARG,1,NSLOW)
      CALL PLTWRI(NSLOW,0.,0.,0.,0.,ARG,1,X(1,1),1)
      CALL PLTWRI(NSLOW,0.,0.,0.,0.,ARG,1,X(1,2),1)
      CALL PLTWRI(NSLOW,0.,0.,0.,0.,ARG,1,X(1,3),1)
 10   CONTINUE
      RETURN
      END
