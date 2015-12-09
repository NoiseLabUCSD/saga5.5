      SUBROUTINE INITS
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      COMPLEX S
      COMPLEX TS2,CCC,TS
c
C     RESET DEPTH-DERIVATIVE FLAG
C
      IDERIV=0
      S=WVNO
      S2=S*S
      THICK(1)=0E0
      THICK(NUML)=0E0
      DO 2 JJJ=1,4
      DO 2 IL=1,NUML
      ROIN(IL,JJJ)=CMPLX(0E0,0E0)
      RUIN(IL,JJJ)=CMPLX(0E0,0E0)
 2    CONTINUE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 I=1,NUMT(1)
      LL=LAYT(I,1)
      ALFA(LL)=CSQRT(S2-AK2(LL,1))
      EZALFM(LL)=CEXP(-ALFA(LL)*THICK(LL))
 10   CONTINUE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 30 I=1,NUMT(3)
      LL=LAYT(I,3)
      ALFA(LL)=CSQRT(S2-AK2(LL,1))
      BETA(LL)=CSQRT(S2-AK2(LL,2))
 30   CONTINUE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 31 I=1,NUMT(3)
      LL=LAYT(I,3)
      CON2(LL)=CON1(LL)*(2E0*S2-AK2(LL,2))
      CON3(LL)=CON1(LL)*2E0*S*ALFA(LL)
      CON4(LL)=CON1(LL)*2E0*S*BETA(LL)
 31   CONTINUE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 32 I=1,NUMT(3)
      LL=LAYT(I,3)
      EZALFM(LL)=CEXP(-ALFA(LL)*THICK(LL))
      EZBETM(LL)=CEXP(-BETA(LL)*THICK(LL))
 32   CONTINUE
C
C     READ IN SLOWNESSES, DISP. AND STRESSES FOR ANISOTRIPIC LAYERS
C
        DO 40 I=1,NUMT(4)
        LAYN=LAYT(I,4)
        CALL RDBUF(41,ANSTD(1,LAYN),52)
 40     CONTINUE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
        DO 41 I=1,NUMT(4)
        LAYN=LAYT(I,4)
        ALFA(LAYN)=DSQ*ANSTD(3,LAYN)
        BETA(LAYN)=DSQ*ANSTD(4,LAYN)
 41     CONTINUE
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
        DO 42 I=1,NUMT(4)
        LAYN=LAYT(I,4)
        EZALFM(LAYN)=CEXP(-THICK(LAYN)*ALFA(LAYN))
        EZBETM(LAYN)=CEXP(-THICK(LAYN)*BETA(LAYN))
 42     CONTINUE
      RETURN
      END
      SUBROUTINE BUILD               
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      COMPLEX S
      S=WVNO
C *** build local coefficient matrices
       IF (NUMT(1).GT.0) CALL LIQLAY
       IF (NUMT(2).GT.0) CALL AIRLAY
       IF (NUMT(3).GT.0) CALL SOLLAY
       IF (NUMT(4).GT.0) CALL TISOLL
C *** SCATTERING PERTURBATION
c      CALL SCATT
C *** MERGE RIGHT HAND SIDES
      DO 420 I=1,NUMI
      DO 410 J=1,4
      R(I,J)=ROIN(I,J)+RUIN(I,J)
 410  CONTINUE
 420  CONTINUE
C *** NORMALIZE DISPLACEMENT EQUATIONS
      DO 530 L=1,NUML
      DO 520 I=1,2
      R(L,I)=R(L,I)*DISNRM
      DO 510 J=1,4
      AUP(L,I,J)=AUP(L,I,J)*DISNRM
      ALO(L,I,J)=ALO(L,I,J)*DISNRM
 510  CONTINUE
 520  CONTINUE
 530  CONTINUE
      RETURN        
      END           

      SUBROUTINE SOLLAY     
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      COMPLEX S
      LOGICAL IFS
      COMPLEX CC,CC0,CC1,CC2,CC3
C     FIRST ROW IS VERTICAL DISPLACEMENT
cvd$r permutation(LAYT)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 I=1,NUMT(3)
       LN=LAYT(I,3)              
       AUP(LN,1,1)=ALFA(LN)  
       AUP(LN,1,2)=-WVNO        
       AUP(LN,1,3)=-ALFA(LN) 
       AUP(LN,1,4)=-WVNO        
C     SECOND ROW IS HORIZONTAL DISPLACEMENT           
       AUP(LN,2,1)=WVNO         
       AUP(LN,2,2)=-BETA(LN) 
       AUP(LN,2,3)=WVNO         
       AUP(LN,2,4)=BETA(LN)  
C     THIRD ROW IS NORMAL STRESS  
       AUP(LN,3,1)=-CON2(LN)
       AUP(LN,3,2)=CON4(LN)
       AUP(LN,3,3)=-CON2(LN)
       AUP(LN,3,4)=-CON4(LN)
C     LAST ROW IS SHEAR STRESS    
       AUP(LN,4,1)=-CON3(LN)
       AUP(LN,4,2)=CON2(LN)
       AUP(LN,4,3)=CON3(LN)
       AUP(LN,4,4)=CON2(LN)
 10   CONTINUE
C     THE LOWER INTERFACE MATRIX FOLLOWS BY A CHANGE OF SIGN
C     AND MULTLIPICATION WITH THE EXPONENTIALS        
      DO 15 I=1,4           
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
       DO 15 J=1,NUMT(3)
        LN=LAYT(J,3)
        ALO(LN,I,1)=-AUP(LN,I,1)*EZALFM(LN)                 
        ALO(LN,I,2)=-AUP(LN,I,2)*EZBETM(LN)                 
        ALO(LN,I,3)=-AUP(LN,I,3)                 
        ALO(LN,I,4)=-AUP(LN,I,4)
        AUP(LN,I,3)=AUP(LN,I,3)*EZALFM(LN)
        AUP(LN,I,4)=AUP(LN,I,4)*EZBETM(LN)                 
 15   CONTINUE
C *** SOURCE TERMS
cvd$  novector
cvd$  noconcur
      DO 40 J=1,NUMTS(3)
       I=NSPNT(J,3)
       LN=LAYS(I)
       IF (.NOT.SHEAR) THEN
        CC0=1E0/ALFA(LN)
        CC1=CON2(LN)
        CC2=-CON3(LN)
        IF (LN.GT.1) THEN
         CC=CPHFAC(I)*CEXP(-ZUS(I)*ALFA(LN))
         CC3=CC*CC0
         RUIN(LN-1,1)=CC+RUIN(LN-1,1)
         RUIN(LN-1,2)=-WVNO*CC3+RUIN(LN-1,2)
         RUIN(LN-1,3)=CC1*CC3+RUIN(LN-1,3)
         RUIN(LN-1,4)=CC2*CC3+RUIN(LN-1,4)
        END IF
        IF (LN.LT.NUML) THEN
         CC=CPHFAC(I)*CEXP(-ZLS(I)*ALFA(LN))
         CC3=CC*CC0
         ROIN(LN,1)=CC+ROIN(LN,1)
         ROIN(LN,2)=WVNO*CC3+ROIN(LN,2)
         ROIN(LN,3)=-CC1*CC3+ROIN(LN,3)
         ROIN(LN,4)=CC2*CC3+ROIN(LN,4)
        END IF
      ELSE
        CC0=1E0/BETA(LN)
        IF (LN.GT.1) THEN
         CC=CPHFAC(I)*CEXP(-ZUS(I)*ALFA(LN))
         CC1=-ZUS(I)*BETA(LN)
         CC1=CPHFAC(I)*CEXP(CC1)
         RUIN(LN-1,1)=-CC*ALFA(LN)+S2*CC0*CC1+RUIN(LN-1,1)
         RUIN(LN-1,2)=WVNO*(CC-CC1)+RUIN(LN-1,2)
         RUIN(LN-1,3)=-CON2(LN)*CC+2E0*CON1(LN)*S2*CC1+RUIN(LN-1,3)
         RUIN(LN-1,4)=CON3(LN)*CC
     &               -CON1(LN)*WVNO*(S2*CC0+BETA(LN))*CC1+RUIN(LN-1,4)
        END IF
        IF (LN.LT.NUML) THEN
         CC=CPHFAC(I)*CEXP(-ZLS(I)*ALFA(LN))
         CC1=-ZLS(I)*BETA(LN)
         CC1=CPHFAC(I)*CEXP(CC1)
         ROIN(LN,1)=ALFA(LN)*CC-S2*CC0*CC1+ROIN(LN,1)
         ROIN(LN,2)=WVNO*(CC-CC1)+ROIN(LN,2)
         ROIN(LN,3)=-CON2(LN)*CC+2E0*CON1(LN)*S2*CC1+ROIN(LN,3)
         ROIN(LN,4)=-CON3(LN)*CC
     &              +CON1(LN)*WVNO*(S2*CC0+BETA(LN))*CC1+ROIN(LN,4)
        END IF
       END IF
 40   CONTINUE
      RETURN                
      END                   
      SUBROUTINE TISOLL     
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      COMPLEX CC,CC0,CC1,CC2,CC3
cvd$r permutation(LAYT)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 I=1,NUMT(4)
       LN=LAYT(I,4)
C     FIRST ROW IS VERTICAL DISPLACEMENT              
       AUP(LN,1,1)=DSQ*ANSTD(7,LN)
       AUP(LN,1,2)=DSQ*ANSTD(8,LN)            
       AUP(LN,1,3)=DSQ*ANSTD(11,LN)
       AUP(LN,1,4)=DSQ*ANSTD(12,LN)
C     SECOND ROW IS HORIZONTAL DISPLACEMENT           
       AUP(LN,2,1)=DSQ*ANSTD(5,LN)
       AUP(LN,2,2)=DSQ*ANSTD(6,LN)
       AUP(LN,2,3)=DSQ*ANSTD(9,LN)
       AUP(LN,2,4)=DSQ*ANSTD(10,LN)
C     THIRD ROW IS NORMAL STRESS  
       AUP(LN,3,1)=CON1(LN)*ANSTD(13,LN)
       AUP(LN,3,2)=CON1(LN)*ANSTD(14,LN)
       AUP(LN,3,3)=CON1(LN)*ANSTD(17,LN)
       AUP(LN,3,4)=CON1(LN)*ANSTD(18,LN)
C     LAST ROW IS SHEAR STRESS    
       AUP(LN,4,1)=CON1(LN)*ANSTD(15,LN)
       AUP(LN,4,2)=CON1(LN)*ANSTD(16,LN)
       AUP(LN,4,3)=CON1(LN)*ANSTD(19,LN)
       AUP(LN,4,4)=CON1(LN)*ANSTD(20,LN)
 10   CONTINUE
C     THE LOWER INTERFACE MATRIX FOLLOWS BY A CHANGE OF SIGN
C     AND MULTLIPICATION WITH THE EXPONENTIALS        
      DO 15 I=1,4
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
       DO 15 J=1,NUMT(4)
        LN=LAYT(J,4)           
        ALO(LN,I,1)=-AUP(LN,I,1)*EZALFM(LN)                 
        ALO(LN,I,2)=-AUP(LN,I,2)*EZBETM(LN)                 
        ALO(LN,I,3)=-AUP(LN,I,3)                 
        ALO(LN,I,4)=-AUP(LN,I,4)                 
        AUP(LN,I,3)=AUP(LN,I,3)*EZALFM(LN)                 
        AUP(LN,I,4)=AUP(LN,I,4)*EZBETM(LN)                 
 15   CONTINUE
C 
C     NO SOURCES ALLOWED IN TRANSVERSE ISOTROPIC LAYERS
C
      RETURN                
      END                   
      SUBROUTINE LIQLAY             
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      COMPLEX CC
      COMPLEX CC1,CC2
C
C     ISOVELOCITY FLUID LAYERS
C
cvd$r permutation(LAYT)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 I=1,NUMT(1)
       LN=LAYT(I,1)
C     FIRST ROW IS VERTICAL DISPLACEMENT  
       AUP(LN,1,1)=ALFA(LN)                
       AUP(LN,1,3)=-ALFA(LN)               
C     SECOND ROW IS HORIZONTAL DISPLACEMENT                  
       AUP(LN,2,1)=WVNO    
       AUP(LN,2,3)=WVNO    
C     THIRD ROW IS NORMAL STRESS          
       AUP(LN,3,1)=CON1(LN)
       AUP(LN,3,3)=CON1(LN)
C     VANISHING SHEAR STRESS
       AUP(LN,4,1)=0E0
       AUP(LN,4,3)=0E0
 10   CONTINUE
C     THE LOWER INTERFACE MATRIX FOLLOWS BY A CHANGE OF SIGN 
C     AND MULTLIPICATION WITH THE EXPONENTIALS               
      DO 15 I=1,4      
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
       DO 15 J=1,NUMT(1)
        LN=LAYT(J,1)
        ALO(LN,I,1)=-AUP(LN,I,1)*EZALFM(LN)     
        ALO(LN,I,3)=-AUP(LN,I,3)     
        AUP(LN,I,3)=AUP(LN,I,3)*EZALFM(LN)     
 15   CONTINUE
C *** SOURCE TERMS
cvd$  novector
cvd$  noconcur
      DO 40 J=1,NUMTS(1)
       I=NSPNT(J,1)
       LN=LAYS(I)
       CC1=1E0/ALFA(LN)
        IF (LN.GT.1) THEN
         CC=CPHFAC(I)*CEXP(-ZUS(I)*ALFA(LN))
c         write(*,*)'before',ruin(ln-1,3)
         RUIN(LN-1,1)=CC+RUIN(LN-1,1)
         RUIN(LN-1,2)=-WVNO*CC*CC1+RUIN(LN-1,2) 
         RUIN(LN-1,3)=-CON1(LN)*CC*CC1+RUIN(LN-1,3)          
        END IF
        IF (LN.LT.NUML) THEN
         CC=CPHFAC(I)*CEXP(-ZLS(I)*ALFA(LN))
         ROIN(LN,1)=CC+ROIN(LN,1)   
         ROIN(LN,2)=WVNO*CC*CC1+ROIN(LN,2)
         ROIN(LN,3)=CON1(LN)*CC*CC1+ROIN(LN,3)             
c         write(*,*)' CC,PHFAC(I),ZUS(I),ALFA(LN)'
c         write(*,*) CC,cPHFAC(I),ZlS(I),ALFA(LN),cc1
c         write(*,*)' press at lower interface',roin(ln,3)
c         write(*,*)' vel lower interface',roin(ln,1), ROIN(LN,2)
        END IF
 40   CONTINUE
      RETURN           
      END              
      SUBROUTINE AIRLAY             
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      COMPLEX CC
C
C     AIRY SOLUTION ADDED 840907
C
C     NOTE THE GRADIENT DEPENDENT COLUMN INTERCHANGE
C
      COMPLEX CC1,CC2
      COMPLEX ZETAU(NLA),AIRYU(NLA),BIRYU(NLA),
     &        AIRYDU(NLA),BIRYDU(NLA),ZTAMU(NLA)
      COMPLEX ZETAL(NLA),AIRYL(NLA),BIRYL(NLA),
     &        AIRYDL(NLA),BIRYDL(NLA),ZTAML(NLA)
      COMPLEX      AIRU,BIRU,AIRDU,BIRDU
      COMPLEX      AIRL,BIRL,AIRDL,BIRDL
      COMPLEX ZETAS,AIRYS,BIRYS,AIRYDS,BIRYDS,ZTAMS
cvd$r permutation(LAYT)
cvd$  cncall
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
cvd$  select(concur)
      DO 10 J=1,NUMT(2)
       LN=LAYT(J,2)
        ZETAU(LN)=CCO(LN)*S2-BCO(LN)
        CALL SCAIRY(ZETAU(LN),AIRYU(LN),BIRYU(LN),
     &              AIRYDU(LN),BIRYDU(LN),ZTAMU(LN))
        ZETAL(LN)=ZETAU(LN)-ACO(LN)*THICK(LN)
        CALL SCAIRY(ZETAL(LN),AIRYL(LN),BIRYL(LN),
     &              AIRYDL(LN),BIRYDL(LN),ZTAML(LN))
        IF (REAL(ACO(LN)).LT.0) THEN
          AISC(LN)=ZTAMU(LN)
          BISC(LN)=ZTAML(LN)
        ELSE
          AISC(LN)=ZTAML(LN)
          BISC(LN)=ZTAMU(LN)
        END IF
        CC1=CEXP(AISC(LN)-ZTAMU(LN))
        CC2=CEXP(ZTAMU(LN)-BISC(LN))
        AIRU=AIRYU(LN)*CC1
        AIRDU=AIRYDU(LN)*CC1
        BIRU=BIRYU(LN)*CC2
        BIRDU=BIRYDU(LN)*CC2
        CC1=CEXP(AISC(LN)-ZTAML(LN))
        CC2=CEXP(ZTAML(LN)-BISC(LN))
        AIRL=AIRYL(LN)*CC1
        AIRDL=AIRYDL(LN)*CC1
        BIRL=BIRYL(LN)*CC2
        BIRDL=BIRYDL(LN)*CC2
        IF (REAL(ACO(LN)).LT.0) THEN
          AUP(LN,1,1)=AIRDU*ACO(LN)
          AUP(LN,1,3)=BIRDU*ACO(LN)
          AUP(LN,2,1)=WVNO
          AUP(LN,3,1)=CON1(LN)*AIRU
          AUP(LN,3,3)=CON1(LN)*BIRU
        ELSE
          AUP(LN,1,1)=BIRDU*ACO(LN)
          AUP(LN,1,3)=AIRDU*ACO(LN)
          AUP(LN,3,1)=CON1(LN)*BIRU
          AUP(LN,3,3)=CON1(LN)*AIRU
        END IF
          AUP(LN,4,1)=0E0
          AUP(LN,4,3)=0E0
        IF (REAL(ACO(LN)).LT.0) THEN
          ALO(LN,1,1)=-AIRDL*ACO(LN)
          ALO(LN,1,3)=-BIRDL*ACO(LN)
          ALO(LN,3,1)=-CON1(LN)*AIRL
          ALO(LN,3,3)=-CON1(LN)*BIRL
        ELSE
          ALO(LN,1,1)=-BIRDL*ACO(LN)
          ALO(LN,1,3)=-AIRDL*ACO(LN)
          ALO(LN,3,1)=-CON1(LN)*BIRL
          ALO(LN,3,3)=-CON1(LN)*AIRL
        END IF
        ALO(LN,4,1)=0E0
        ALO(LN,4,3)=0E0
 10   CONTINUE
C
C     SOURCES IN AIRY LAYERS INCLUDED 850214
C
cvd$  noconcur     
cvd$  novector
       DO 15 J=1,NUMTS(2)
        I=NSPNT(J,2)
        LN=LAYS(I)
        ZETAS=CCO(LN)*S2-ACO(LN)*ZUS(I)-BCO(LN)
        CALL SCAIRY(ZETAS,AIRYS,BIRYS,AIRYDS,BIRYDS,ZTAMS)
        CC=2E0*CPHFAC(I)/((-AIRYS*BIRYDS
     1       +BIRYS*AIRYDS)*(-ACO(LN)))
        IF (IDERIV.NE.1) THEN
          CSAIR1(I)=CC*BIRYS
          CSAIR2(I)=CC*AIRYS
        ELSE
          CSAIR1(I)=-ACO(LN)*CC*BIRYDS
          CSAIR2(I)=-ACO(LN)*CC*AIRYDS
        END IF
        CSAIR3(I)=ZTAMS
        CC1=CSAIR1(I)*CEXP(ZTAMS-BISC(LN))
        CC2=CSAIR2(I)*CEXP(AISC(LN)-ZTAMS)
        IF (REAL(ACO(LN)).GE.0) THEN
C       NEGATIVE VELOCITY GRADIENT
          RUIN(LN-1,1)=RUIN(LN-1,1)-ACO(LN)*CC1*AIRYDU(LN)
          RUIN(LN-1,3)=RUIN(LN-1,3)-CON1(LN)*CC1*AIRYU(LN)
          ROIN(LN,1)=ROIN(LN,1)+ACO(LN)*CC2*BIRYDL(LN)
          ROIN(LN,3)=ROIN(LN,3)+CON1(LN)*CC2*BIRYL(LN)
        ELSE
          RUIN(LN-1,1)=RUIN(LN-1,1)-ACO(LN)*CC2*BIRYDU(LN)
          RUIN(LN-1,3)=RUIN(LN-1,3)-CON1(LN)*CC2*BIRYU(LN)
          ROIN(LN,1)=ROIN(LN,1)+ACO(LN)*CC1*AIRYDL(LN)
          ROIN(LN,3)=ROIN(LN,3)+CON1(LN)*CC1*AIRYL(LN)
        END IF
 15    CONTINUE
      RETURN
      END              
      SUBROUTINE KERNEL(CKERN,NRCV)        
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      COMPLEX CKERN(NRCV,3)
      COMPLEX ERALFA,ERBETA,ERALFM,ERBETM 
      COMPLEX CC,CC1,CC2,CC3,CC4,CWUFAC
      COMPLEX ZETA(NRD),AIRY(NRD),BIRY(NRD),
     &        AIRYD(NRD),BIRYD(NRD),ZTAM(NRD)
C
      CWUFAC=AI*DSQ*PCORR
C *** RECEIVERS IN ISOVELOCITY FLUID LAYERS
cvd$  permutation(NUMTR)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 J=1,NUMTR(1)
       INT=NRPNT(J,1)
       LL=LAY(INT)      
       ZZ=Z(INT)        
       IF (LL.NE.1) THEN
         ERALFM=CEXP(-ZZ*ALFA(LL))
       ELSE
         ERALFM=0E0
       END IF                 
       IF (LL.NE.NUML) THEN
         ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
       ELSE
         ERALFA=0E0
       END IF
       CC1=SS(LL,1)*ERALFM
       CC3=SS(LL,3)*ERALFA
       CKERN(INT,1)=-CON1(LL)*(CC1+CC3)
       CKERN(INT,2)=ALFA(LL)*(-CC1+CC3)
       CKERN(INT,3)=-WVNO*(CC1+CC3)
C *** SOURCE CONTRIBUTION
       IF (NOSOU(LL).GT.0) THEN
         CC=1E0/ALFA(LL)
         CC1=-CON1(LL)*CC
         CC3=-WVNO*CC
         DO 120 I=IFSOU(LL),ILSOU(LL)
           CC=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*ALFA(LL))
           CKERN(INT,1)=CKERN(INT,1)+CC*CC1
           CKERN(INT,2)=CKERN(INT,2)-SIGN(1.0,ZZ-ZUS(I))*CC
           CKERN(INT,3)=CKERN(INT,3)+CC*CC3
 120     CONTINUE
       END IF
 10    CONTINUE
c          write(6,*)'ckern(1,1)',ckern(1,1),int
C
C     AIRY SOLUTION IMPLEMENTED 840907
C
C *** RECEIVERS IN non-ISOVELOCITY FLUID LAYERS
cvd$  permutation(NRPNT)
cvd$  cncall
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
cvd$  select(concur)
      DO 20 J=1,NUMTR(2)
       INT=NRPNT(J,2)
       LL=LAY(INT)      
       ZZ=Z(INT)        
       ZETA(INT)=CCO(LL)*S2-ZZ*ACO(LL)-BCO(LL)
       CALL SCAIRY(ZETA(INT),AIRY(INT),BIRY(INT),
     &             AIRYD(INT),BIRYD(INT),ZTAM(INT))
        CC1=CEXP(AISC(LL)-ZTAM(INT))
        CC2=CEXP(ZTAM(INT)-BISC(LL))
        IF ((REAL(ACO(LL))).LT.0) THEN
         CC3=SS(LL,1)*CC1*AIRY(INT)+SS(LL,3)*CC2*BIRY(INT)
         CC4=SS(LL,1)*CC1*AIRYD(INT)+SS(LL,3)*CC2*BIRYD(INT)
        ELSE
         CC3=SS(LL,1)*CC2*BIRY(INT)+SS(LL,3)*CC1*AIRY(INT)
         CC4=SS(LL,1)*CC2*BIRYD(INT)+SS(LL,3)*CC1*AIRYD(INT)
        END IF
        CKERN(INT,1)=-CON1(LL)*CC3
        CKERN(INT,2)=-ACO(LL)*CC4
        CKERN(INT,3)=-WVNO*CC3          
C
C      SOURCES INCLUDED IN AIRY LAYERS 840214
C
        IF (NOSOU(LL).GT.0) THEN
         DO 5 I=IFSOU(LL),ILSOU(LL)
         IF (REAL(ACO(LL))*(ZZ-ZUS(I)).LE.0.) THEN
          CC3=CEXP(CSAIR3(I)-ZTAM(INT))
          CC1=AIRY(INT)*CSAIR1(I)*CC3
          CC2=AIRYD(INT)*CSAIR1(I)*CC3
         ELSE
          CC3=CEXP(ZTAM(INT)-CSAIR3(I))
          CC1=BIRY(INT)*CSAIR2(I)*CC3
          CC2=BIRYD(INT)*CSAIR2(I)*CC3
         END IF
         CKERN(INT,1)=CKERN(INT,1)-CON1(LL)*CC1
         CKERN(INT,2)=CKERN(INT,2)-ACO(LL)*CC2
         CKERN(INT,3)=CKERN(INT,3)-WVNO*CC1
 5       CONTINUE
        END IF
 20   CONTINUE
C *** RECEIVERS IN SOLID LAYERS
cvd$  permutation(NRPNT)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 30 J=1,NUMTR(3)
       INT=NRPNT(J,3)
       LL=LAY(INT)      
       ZZ=Z(INT)        
        IF (LL.NE.1) THEN
         ERALFM=CEXP(-ZZ*ALFA(LL))                 
         ERBETM=CEXP(-ZZ*BETA(LL))            
        ELSE
         ERALFM=0E0
         ERBETM=0E0
        END IF
        IF (LL.NE.NUML) THEN
         ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
         ERBETA=CEXP((ZZ-THICK(LL))*BETA(LL))
        ELSE
         ERALFA=0E0
         ERBETA=0E0
        END IF
        CC1=SS(LL,1)*ERALFM
        CC2=SS(LL,2)*ERBETM
        CC3=SS(LL,3)*ERALFA
        CC4=SS(LL,4)*ERBETA
        CKERN(INT,1)=CON2(LL)*(CC1+CC3)+CON4(LL)*(CC4-CC2)
        CKERN(INT,2)=ALFA(LL)*(CC3-CC1)+WVNO*(CC2+CC4)
        CKERN(INT,3)=-WVNO*(CC1+CC3)+BETA(LL)*(CC2-CC4)
 30   CONTINUE
C *** SOURCE CONTRIBUTION
cvd$  permutation(NRPNT)
cvd$  concur
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 35 J=1,NUMTR(3)
       INT=NRPNT(J,3)
       LL=LAY(INT)      
       ZZ=Z(INT)        
        IF (NOSOU(LL).GT.0) THEN
         IF (.NOT.SHEAR) THEN
          CC=1E0/ALFA(LL)
          CC1=CON2(LL)*CC
          CC3=-WVNO*CC
cvd$  novector
cvd$  noconcur
          DO 320 I=IFSOU(LL),ILSOU(LL)
           CC=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*ALFA(LL))
           CKERN(INT,1)=CKERN(INT,1)+CC*CC1
           CKERN(INT,2)=CKERN(INT,2)-SIGN(1.0,ZZ-ZUS(I))*CC
           CKERN(INT,3)=CKERN(INT,3)+CC*CC3
 320      CONTINUE
         ELSE
          IF (IOUT(2).GT.0) CC=S2/BETA(LL)
cvd$  novector
cvd$  noconcur
          DO 330 I=IFSOU(LL),ILSOU(LL)
           CC1=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*ALFA(LL))
           CC2=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*BETA(LL))
           CKERN(INT,1)=CKERN(INT,1)+SIGN(1.0,ZZ-ZUS(I))*
     &       ( CON2(LL)*CC1-CON1(LL)*2E0*S2*CC2  )
           CKERN(INT,2)=CKERN(INT,2)-ALFA(LL)*CC1+CC*CC2
           CKERN(INT,3)=CKERN(INT,3)
     &                 -SIGN(1.0,ZZ-ZUS(I))*WVNO*(CC1-CC2)
 330      CONTINUE
         END IF
        END IF
 35   continue
C *** RECEIVERS IN TRANSVERSILY ISOTROPIC LAYERS
cvd$  permutation(NRPNT)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 40 J=1,NUMTR(4)
       INT=NRPNT(J,4)
       LL=LAY(INT)      
       ZZ=Z(INT)        
        IF (LL.NE.1) THEN
         ERALFM=CEXP(-ZZ*ALFA(LL))                 
         ERBETM=CEXP(-ZZ*BETA(LL))            
        ELSE
         ERALFM=0E0
         ERBETM=0E0
        END IF
        IF (LL.NE.NUML) THEN
         ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
         ERBETA=CEXP((ZZ-THICK(LL))*BETA(LL))
        ELSE
         ERALFA=0E0
         ERBETA=0E0
        END IF
        CC1=SS(LL,1)*ERALFM
        CC2=SS(LL,2)*ERBETM
        CC3=SS(LL,3)*ERALFA
        CC4=SS(LL,4)*ERBETA
        IF (IOUT(1).GT.0) CKERN(INT,1)=-CON1(LL)*
     &       (ANSTD(13,LL)*CC1+ANSTD(14,LL)*CC2+
     &        ANSTD(17,LL)*CC3+ANSTD(18,LL)*CC4)
        IF (IOUT(2).GT.0) CKERN(INT,2)=-DSQ*
     &       (ANSTD(7,LL)*CC1+ANSTD(8,LL)*CC2+
     &        ANSTD(11,LL)*CC3+ANSTD(12,LL)*CC4)
        IF (IOUT(3).GT.0) CKERN(INT,3)=-DSQ*
     &       (ANSTD(5,LL)*CC1+ANSTD(6,LL)*CC2+
     &        ANSTD(9,LL)*CC3+ANSTD(10,LL)*CC4)
 40   continue
C *** CONVERT DISPLACEMENTS TO VELOCITIES
      DO 500 J=1,IR
       CKERN(J,2)=CKERN(J,2)*CWUFAC
       CKERN(J,3)=-DSQ*PCORR*CKERN(J,3)
 500  CONTINUE
      RETURN           
      END              

      SUBROUTINE KERDEC(CKERN,NRCV)        
C 
C     DETERMINES THE DECOMPOSED WAVEFIELD KERNELS.
C     ON EXIT, THE PARAMETERS ARE PLACED IN ARRAY CKERN(I,J,K)
C     AS FOLLOWS:
C
C     CKERN(I,1,K)	NORMAL STRESS, RECEIVER I.
C     CKERN(I,2,K)      VERTICAL PARTICLE VELOCITY, RECEIVER I.
C     CKERN(I,3,K)	HORIZONTAL PARTICLE VELOCITY, RECEIVER I.
C
C     CKERN(I,J,1)	TOTAL KERNEL
C     CKERN(I,J,2)	DOWN-GOING COMPRESSIONAL WAVES ONLY
C     CKERN(I,J,3) 	DOWN-GOING SHEAR WAVES ONLY
C     CKERN(I,J,4)	UP-GOING COMPRESSIONAL WAVES ONLY
C     CKERN(I,J,5)	UP-GOING SHEAR WAVES ONLY
C 
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
      COMPLEX CKERN(NRCV,3,5)
      COMPLEX ERALFA,ERBETA,ERALFM,ERBETM 
      COMPLEX CC,CC1,CC2,CC3,CC4,CC5,CC6,CWUFAC
      COMPLEX ZETA,AIRY,BIRY,AIRYD,BIRYD,ZTAM
c
      CALL VCLR(CKERN,1,30*NRCV)
      CWUFAC=AI*DSQ*PCORR
C *** RECEIVERS IN ISOVELOCITY FLUID LAYERS
cvd$  permutation(NUMTR)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 10 J=1,NUMTR(1)
       INT=NRPNT(J,1)
       LL=LAY(INT)      
       ZZ=Z(INT)        
       IF (LL.NE.1) THEN
         ERALFM=CEXP(-ZZ*ALFA(LL))
       ELSE
         ERALFM=0E0
       END IF                 
       IF (LL.NE.NUML) THEN
         ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
       ELSE
         ERALFA=0E0
       END IF
        CC1=SS(LL,1)*ERALFM
        CC3=SS(LL,3)*ERALFA
        CKERN(INT,1,2)=-CON1(LL)*CC1
        CKERN(INT,1,4)=-CON1(LL)*CC3
        CKERN(INT,2,2)=-ALFA(LL)*CC1
        CKERN(INT,2,4)= ALFA(LL)*CC3
        CKERN(INT,3,2)=-WVNO*CC1
        CKERN(INT,3,4)=-WVNO*CC3
C *** SOURCE CONTRIBUTION
        IF (NOSOU(LL).GT.0) THEN
          CC=1E0/ALFA(LL)
          CC1=-CON1(LL)*CC
          CC3=-WVNO*CC
          DO 120 I=IFSOU(LL),ILSOU(LL)
           ISS=2+NINT(1.0+SIGN(1.0,ZUS(I)-ZZ))
           CC=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*ALFA(LL))
           CKERN(INT,1,ISS)=CKERN(INT,1,ISS)+CC*CC1
           CKERN(INT,2,ISS)=CKERN(INT,2,ISS)
     &                       -SIGN(1.0,ZZ-ZUS(I))*CC
           CKERN(INT,3,ISS)=CKERN(INT,3,ISS)+CC*CC3
 120      CONTINUE
        END IF
 10   CONTINUE
C
C     AIRY SOLUTION IMPLEMENTED 840907
C
C *** RECEIVERS IN non-ISOVELOCITY FLUID LAYERS
cvd$  permutation(NUMTR)
      DO 20 J=1,NUMTR(2)
       INT=NRPNT(J,2)
       LL=LAY(INT)      
       ZZ=Z(INT)        
        ZETA=CCO(LL)*S2-ZZ*ACO(LL)-BCO(LL)
        CALL SCAIRY(ZETA,AIRY,BIRY,AIRYD,BIRYD,ZTAM)
        CC5=CEXP(AISC(LL)-ZTAM)
        CC6=CEXP(ZTAM-BISC(LL))
        IF ((REAL(ACO(LL))).LT.0) THEN
         CC1=SS(LL,1)*AIRY*CC5
         CC2=SS(LL,3)*BIRY*CC6
         CC3=SS(LL,1)*AIRYD*CC5
         CC4=SS(LL,3)*BIRYD*CC6
        ELSE
         CC1=SS(LL,1)*BIRY*CC6
         CC2=SS(LL,3)*AIRY*CC5
         CC3=SS(LL,1)*BIRYD*CC6
         CC4=SS(LL,3)*AIRYD*CC5
        END IF
        CKERN(INT,1,2)=-CON1(LL)*CC1
        CKERN(INT,1,4)=-CON1(LL)*CC2
        CKERN(INT,2,2)=-ACO(LL)*CC3
        CKERN(INT,2,4)=-ACO(LL)*CC4
        CKERN(INT,3,2)=-WVNO*CC1
        CKERN(INT,3,4)=-WVNO*CC2          
C
C      SOURCES INCLUDED IN AIRY LAYERS 840214
C
        IF (NOSOU(LL).GT.0) THEN
         CC1=-CON1(LL)*AIRY
         CC2=-CON1(LL)*BIRY
cvd$  novector
         DO 5 I=IFSOU(LL),ILSOU(LL)
         ISS=2+NINT(1.0+SIGN(1.0,ZUS(I)-ZZ))
         IF (REAL(ACO(LL))*(ZZ-ZUS(I)).LE.0.) THEN
          CC3=CEXP(CSAIR3(I)-ZTAM)
          CC1=AIRY*CSAIR1(I)*CC3
          CC2=AIRYD*CSAIR1(I)*CC3
         ELSE
          CC3=CEXP(ZTAM-CSAIR3(I))
          CC1=BIRY*CSAIR2(I)*CC3
          CC2=BIRYD*CSAIR2(I)*CC3
         END IF
         CKERN(INT,1,ISS)=CKERN(INT,1,ISS)-CON1(LL)*CC1
         CKERN(INT,2,ISS)=CKERN(INT,2,ISS)-ACO(LL)*CC2
         CKERN(INT,3,ISS)=CKERN(INT,3,ISS)-WVNO*CC1
 5       CONTINUE
        END IF
 20   CONTINUE
C *** RECEIVERS IN SOLID LAYERS
cvd$  permutation(NRPNT)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 30 J=1,NUMTR(3)
       INT=NRPNT(J,3)
       LL=LAY(INT)      
       ZZ=Z(INT)        
        IF (LL.NE.1) THEN
         ERALFM=CEXP(-ZZ*ALFA(LL))                 
         ERBETM=CEXP(-ZZ*BETA(LL))            
        ELSE
         ERALFM=0E0
         ERBETM=0E0
        END IF
        IF (LL.NE.NUML) THEN
         ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
         ERBETA=CEXP((ZZ-THICK(LL))*BETA(LL))
        ELSE
         ERALFA=0E0
         ERBETA=0E0
        END IF
        CC1=SS(LL,1)*ERALFM
        CC2=SS(LL,2)*ERBETM
        CC3=SS(LL,3)*ERALFA
        CC4=SS(LL,4)*ERBETA
        CKERN(INT,1,2)=CON2(LL)*CC1
        CKERN(INT,1,4)=CON2(LL)*CC3
        CKERN(INT,1,3)=-CON4(LL)*CC2
        CKERN(INT,1,5)=CON4(LL)*CC4
        CKERN(INT,2,2)=-ALFA(LL)*CC1
        CKERN(INT,2,4)=ALFA(LL)*CC3
        CKERN(INT,2,3)=WVNO*CC2
        CKERN(INT,2,5)=WVNO*CC4
        CKERN(INT,3,2)=-WVNO*CC1
        CKERN(INT,3,4)=-WVNO*CC3
        CKERN(INT,3,3)=BETA(LL)*CC2
        CKERN(INT,3,5)=-BETA(LL)*CC4
 30   CONTINUE
C *** SOURCE CONTRIBUTION
cvd$  permutation(NRPNT)
cvd$  concur
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 35 J=1,NUMTR(3)
       INT=NRPNT(J,3)
       LL=LAY(INT)      
       ZZ=Z(INT)        
        IF (NOSOU(LL).GT.0) THEN
         IF (.NOT.SHEAR) THEN
          CC=1E0/ALFA(LL)
          CC1=CON2(LL)*CC
          CC3=-WVNO*CC
cvd$  novector
cvd$  noconcur
          DO 320 I=IFSOU(LL),ILSOU(LL)
           ISS=2+NINT(1.0+SIGN(1.0,ZUS(I)-ZZ))
           CC=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*ALFA(LL))
           CKERN(INT,1,ISS)=CKERN(INT,1,ISS)+CC*CC1
           CKERN(INT,2,ISS)=CKERN(INT,2,ISS)-SIGN(1.0,ZZ-ZUS(I))*CC
           CKERN(INT,3,ISS)=CKERN(INT,3,ISS)+CC*CC3
 320      CONTINUE
         ELSE
          CC=S2/BETA(LL)
cvd$  novector
cvd$  noconcur
          DO 330 I=IFSOU(LL),ILSOU(LL)
           ISS=2+NINT(1.0+SIGN(1.0,ZUS(I)-ZZ))
           CC1=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*ALFA(LL))
           CC2=CPHFAC(I)*CEXP(-ABS(ZUS(I)-ZZ)*BETA(LL))
           CKERN(INT,1,ISS)=CKERN(INT,1,ISS)+SIGN(1.0,ZZ-ZUS(I))*
     &                      CON2(LL)*CC1
           CKERN(INT,1,ISS+1)=CKERN(INT,1,ISS+1)+SIGN(1.0,ZZ-ZUS(I))*
     &                        (-CON1(LL)*2E0*S2*CC2)
           CKERN(INT,2,ISS)=CKERN(INT,2,ISS)-ALFA(LL)*CC1
           CKERN(INT,2,ISS+1)=CKERN(INT,2,ISS+1)+CC*CC2
           CKERN(INT,3,ISS)=CKERN(INT,3,ISS)
     &                 -SIGN(1.0,ZZ-ZUS(I))*WVNO*CC1
           CKERN(INT,3,ISS+1)=CKERN(INT,3,ISS+1)
     &                 +SIGN(1.0,ZZ-ZUS(I))*WVNO*(CC2)
 330      CONTINUE
         END IF
        END IF
 35    CONTINUE
C *** RECEIVERS IN TRANSVERSILY ISOTROPIC LAYERS
cvd$  permutation(NRPNT)
cvd$  nodepchk
CDEC$ INIT_DEP_FWD
      DO 40 J=1,NUMTR(4)
       INT=NRPNT(J,4)
       LL=LAY(INT)      
       ZZ=Z(INT)        
        IF (LL.NE.1) THEN
         ERALFM=CEXP(-ZZ*ALFA(LL))                 
         ERBETM=CEXP(-ZZ*BETA(LL))            
        ELSE
         ERALFM=0E0
         ERBETM=0E0
        END IF
        IF (LL.NE.NUML) THEN
         ERALFA=CEXP((ZZ-THICK(LL))*ALFA(LL))
         ERBETA=CEXP((ZZ-THICK(LL))*BETA(LL))
        ELSE
         ERALFA=0E0
         ERBETA=0E0
        END IF
        CC1=SS(LL,1)*ERALFM
        CC2=SS(LL,2)*ERBETM
        CC3=SS(LL,3)*ERALFA
        CC4=SS(LL,4)*ERBETA
          CKERN(INT,1,2)=-CON1(LL)*ANSTD(13,LL)*CC1
          CKERN(INT,1,3)=-CON1(LL)*ANSTD(14,LL)*CC2
          CKERN(INT,1,4)=-CON1(LL)*ANSTD(17,LL)*CC3
          CKERN(INT,1,5)=-CON1(LL)*ANSTD(18,LL)*CC4
          CKERN(INT,2,2)=-DSQ*ANSTD(7,LL)*CC1
          CKERN(INT,2,3)=-DSQ*ANSTD(8,LL)*CC2
          CKERN(INT,2,4)=-DSQ*ANSTD(11,LL)*CC3
          CKERN(INT,2,5)=-DSQ*ANSTD(12,LL)*CC4
          CKERN(INT,3,2)=-DSQ*ANSTD(5,LL)*CC1
          CKERN(INT,3,3)=-DSQ*ANSTD(6,LL)*CC2
          CKERN(INT,3,4)=-DSQ*ANSTD(9,LL)*CC3
          CKERN(INT,3,5)=-DSQ*ANSTD(10,LL)*CC4
 40   CONTINUE
C *** CONVERT DISPLACEMENTS TO VELOCITIES
      DO 490 IIS=2,5
       DO 490 J=1,IR
       CKERN(J,2,IIS)=CKERN(J,2,IIS)*CWUFAC
       CKERN(J,3,IIS)=-DSQ*PCORR*CKERN(J,3,IIS)
 490  CONTINUE
C *** TOTAL FIELD
      DO 500 IIS=2,5
       DO 500 IIP=1,3
        DO 500 J=1,IR
         CKERN(J,IIP,1)=CKERN(J,IIP,1)+CKERN(J,IIP,IIS)
 500   CONTINUE
      RETURN           
      END              
      SUBROUTINE SOLVE                  
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      CALL CVIMOV(ALO,INDA,1,WORK1,2,NNA)
      CALL CVFILL(CNUL,WORK2,2,NNB)
      CALL CVMOVI(WORK1,2,INDB,1,WORK2,NNA)
      CALL CVIMOV(R,INDR,1,RHS,2,NEQ)
c      IF (DEBUG) THEN
c        DO 10 IS=1,1        !LSTOT
c        WRITE(*,*) 'right hand side for source no',IS
c        DO 10 I=1,NEQ
c          WRITE(*,*)I,RHS(I+(IS-1)*NEQ)
c10      CONTINUE
c      ENDIF
      EPS=1E-10
      CALL CBGEBS(WORK2,RHS,NEQ,NEQ,IBW,EPS)
      IERR=EPS
c      IF (DEBUG) THEN
c        DO 20 IS=1,1        !LSTOT
c        WRITE(*,*) 'solution vector for source no',IS
c        DO 20 I=1,NEQ
c          WRITE(*,*)I,RHS(I+(IS-1)*NEQ)
c20      CONTINUE
c      ENDIF
        IF (IERR.NE.0) RETURN
 999  CALL CVMOVI(RHS,2,INDS,1,SS,NEQ)
      RETURN        
      END           
