      SUBROUTINE INENVI
C
C     SUBROUTINE FOR READING IN ENVIRONMENTAL DATA
C
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'
C           
C ********************FORMATS*******************               
C           
C           
 350  FORMAT('     DEPTH        ALPHA       BETA      ATTENA      ',
     1       ' ATTENB         RHO      ROUGHNESS')
 351  FORMAT(3F12.5,2F12.8,2E12.5)              
C           
C           
C **********************************************               
C           
C     DEFAULT: ISOTROPIC LAYERS
C
      NTISOL=0
C           
C     NUML IS THE NUMBER OF LAYERS INCLUDING THE HALF SPACES.  
C     NUMI IS THE NUMBER OF INTERFACES        
      DO 1 I=1,NLTYP
 1      NUMT(I)=0
      IF (NUML.GT.NLA) THEN
      WRITE(6,*) '*** TOO MANY LAYERS ***'
      STOP
      END IF
      NUMI=NUML-1
c      write(*,*) 'flagpu',flagpu
      if (flagpu.lt.2) then
        WRITE (7,*)'     DEPTH        ALPHA       BETA',
     1       '      ATTENA       ATTENB         RHO      ROUGHNESS'
      endif
C
C     READ IN ENVIRONMENTAL DATA
      DO 110 M=1,NUML            
C     IF (ROUGH(M).LT.-1E-10) THEN
C
C     NON-KIRCHHOFF SCATTERING
C     BACKSPACE AND READ ALSO CORRELLATION LENGTH
C        BACKSPACE(1)
C        READ(1,*) (V(M,N),N=1,6),ROUGH(M),CLEN(M)
C        WRITE(6,351) (V(M,N),N=1,6),ROUGH(M)
C        WRITE(6,*) 'ROUGHNESS CORRELATION LENGTH (m):',CLEN(M)
C      ELSE
         CLEN(M)=1E6
C      END IF
c--- change of slope for preceeding layer
        IF (M.GT.2) THEN
           IF ((V(M-1,3).gt.-1000.000000) .and.
     &          (V(M-1,3).lt.-999.998)) V(M-1,3)=-V(M,2)
c          write(6,*)'celerite v(M,3)',M-1,v(m-1,3)
        END IF
C     TYPE OF LAYER
      IF (ABS(V(M,2)).LT.1E-10) THEN
         LAYTYP(M)=-1
      ELSE IF (V(M,2).GT.0E0.AND.ABS(V(M,3)).LT.1E-10) THEN
         LAYTYP(M)=1
         NUMT(1)=NUMT(1)+1
         LAYT(NUMT(1),1)=M
      ELSE IF (V(M,2).GT.0E0.AND.V(M,3).LT.0) THEN
         LAYTYP(M)=2
         NUMT(2)=NUMT(2)+1
         LAYT(NUMT(2),2)=M
      ELSE IF (V(M,2).GT.0E0.AND.V(M,3).GT.0) THEN
         LAYTYP(M)=3
         NUMT(3)=NUMT(3)+1
         LAYT(NUMT(3),3)=M
C
C     CHECK WHETHER PARAMETERS ARE PHYSICALLY MEANINGFUL
C
         GAMMA2=(V(M,3)/V(M,2))**2
         IF (GAMMA2.GT.0.75) THEN
         WRITE(6,*) '>>>>>WARNING: UNPHYSICAL SPEED RATIO, Cs/Cp>0.75'
         END IF
         IF ((GAMMA2*V(M,5)).GT.(0.75*V(M,4))) THEN
            WRITE(6,*) 
     &   '>>>>>WARNING: UNPHYSICAL ATTENUATION, (As/Ap)*(Cs/Cp)**2>0.75'
         END IF
      ELSE IF (V(M,2).LT.0E0) THEN
         write(*,*)' *** layertype not supported'
c     LAYTYP(M)=4
c     NUMT(4)=NUMT(4)+1
c     LAYT(NUMT(4),4)=M
c     NTISOL=NTISOL+1
c     LAYNTI(NTISOL)=M
c     c        READ(1,*) NSL(M)
c     AHSUM=0E0
c     WRITE(6,*) 'LAYER',M,' IS TRANSV. ISOTROPIC'
c     WRITE(6,*) 'DEPTH: ',V(M,1),' M'
c     WRITE(6,*) 'LAYERS:',NSL(M)
c     WRITE(6,352)
c 352  FORMAT(1H ,8X,'   ALPHA         BETA         ATTNA        ATTNB'
c     &       ,'        RO       FRACTION')       
c     AMEAN=0E0
c     BMEAN=0E0
c     DO 109 NNN=1,NSL(M)
c     READ(1,*) RASP,RBSP,ATTA,ATTB,
c     &            ARO(NNN,M),AH(NNN,M)
c     AMEAN=AMEAN+RASP
c     BMEAN=BMEAN+RBSP
c     WRITE(6,353) RASP,RBSP,ATTA,ATTB,ARO(NNN,M),AH(NNN,M)
c     353    FORMAT((/1H ,8X,2(3X,F9.2),2(3X,F9.6),2(3X,F9.3)))
C
C     CONVERT ATTENUATIONS TO IMAGINARY PART OF VELOCITY
C     
c     CIMA=RASP*ATTA/54.57512
c     CIMB=RBSP*ATTB/54.57512
c     ASP(NNN,M)=CMPLX(RASP,CIMA)
c     BSP(NNN,M)=CMPLX(RBSP,CIMB)
c     AHSUM=AHSUM+AH(NNN,M)
c     109    CONTINUE
c     V(M,2)=AMEAN/NSL(M)
c     V(M,3)=BMEAN/NSL(M)
c     IF (ABS(AHSUM-1E0).GT.1E-4) 
c     &     STOP '*** SUM OF FRACTIONS IS NOT 1.0 ***'
      ELSE
         STOP '*** UNKNOWN LAYER TYPE ***'
      END IF
 110  CONTINUE
      ROUGH(1)=ROUGH(2)
      DO 111 M=1,NUML
         ROUGH2(M)=ROUGH(M)**2
        if ((flagpu.lt.2).or.(debug))
     1      WRITE(7,351) (V(M,N),N=1,6),ROUGH(M)
 111  CONTINUE   
c      DO 1111 M=2,NUML
c      IF (ROUGH2(M).GT.1E-10) THEN
c        IF (LAYTYP(M-1).EQ.2.OR.LAYTYP(M).EQ.2) THEN
c          WRITE(6,*) '**** ROUGHNESS NOT ALLOWED BETWEEN LAYERS ****'
c          WRITE(6,*) '**** WITH SOUND SPEED GRADIENT. INSERT    ****'
c          WRITE(6,*) '**** A DUMMY ISOVELOCITY LAYER.           ****'
c          STOP '>>>> EXECUTION TERMINATED <<<<'
c        ELSE IF (LAYTYP(M-1).EQ.4.OR.LAYTYP(M).EQ.4) THEN
c          WRITE(6,*) '**** NON-KIRCHHOFF SCATTERING NOT INCLUDED ***'
c          WRITE(6,*) '**** FOR ANISOTROPIC LAYERS. KIRCHHOFF     ***'
c          WRITE(6,*) '**** APPROXIMATION WILL BE APPLIED.        ***'
c          CLEN(M)=1E6
c        ELSE
c        END IF
c      END IF
c 1111 CONTINUE
C
C     THE DENSITIES ARE CONVERTED FROM G/CM**3 TO KG/M**3      
C           
c      DO 112 M=1,NUML
c      IF (LAYTYP(M).NE.4) THEN
c        V(M,6)=V(M,6)*1E3      
c      ELSE 
c        DO 212 N=1,NSL(M)
c 212    ARO(N,M)=1E3*ARO(N,M)
c      END IF
c 112  CONTINUE
C
C     CALCULATE ELASTIC CONSTANTS FOR TRANS. ISOTR. LAYERS
C
      IF (NTISOL.GT.0) THEN
c         CALL CALELC          
         STOP '*** ANISOTROPIC MEDIA NOT ALLOWED IN THIS VERSION ***'
      END IF
C           
      RETURN
      END


      SUBROUTINE RECEIV(VL,NUMLL,RD,LR,Z)                    
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      V(1,1)=V(2,1)                   
      DO 10 I=2,NUML                  
      IF (V(I,1)-RD) 10,10,20         
 10   CONTINUE   
      LR=NUML    
      Z=RD-V(NUML,1)                  
      RETURN     
 20   LR=I-1     
      Z=RD-V(LR,1)                    
      RETURN     
      END        
      SUBROUTINE SOURCE(VL,NUMLL,SD,LSL,ZUP,ZLO)              
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      DO 10 I=2,NUML                  
      IF (V(I,1)-SD) 10,20,20         
 10   CONTINUE   
      LSL=NUML    
      ZUP=SD-V(NUML,1)                
      ZLO=-ZUP     
      RETURN     
 20   LSL=I-1     
      ZUP=SD-V(LSL,1)                  
      ZLO=V(I,1)-SD                   
      RETURN     
      END        
      SUBROUTINE PINIT1         
C     INITIALIZATION OF VARIABLES               
C               
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'
c
      IERR=0
      IOERR=0
      IFN=0
      PI=4.0*ATAN(1E0)
      CNUL=CMPLX(0.,0.)         
C
C     LAYER CLASSIFICATION:
C
C     LAYTYP=1:   ISOVELOCITY FLUID: V(I,3) = 0
C
C     LAYTYP=2:   1/C**2 LINEAR:     V(I,3) < 0 
C               -V(I,3) = VEL. AT LOWER BOUNDARY
C
C     LAYTYP=3:   ISOVELOCITY SOLID: V(I,3) > 0
C
C     LAYTYP=4:   TRANSVERSE ISOTROPIC SOLID: V(I,3) > 0
C
      LSOLF=NUML+1
      LSOLL=0
C
C     UPPER AND LOWER HALF SPACES MUST BE ISOVELOCITY
C
      IF (LAYTYP(1).EQ.2) THEN
        LAYTYP(1)=1
        V(1,3)=0
      END IF
      IF (LAYTYP(NUML).EQ.2) THEN
        LAYTYP(NUML)=1
        V(NUML,3)=0
      END IF
      DO 10 I=1,NUML
C
C     FOR SMALL GRADIENTS LAYER TREATED AS ISOVELOCITY
C   pg 14/7/93 but the laytyp was not changed, only v(i,3)-- this must be
c   an error 
c      IF (ABS(V(I,2)+V(I,3)).LT.1E-3) V(I,3)=0E0
      IF ((laytyp(i).eq.2).and.(ABS(V(I,2)+V(I,3)).LT.1E-3)) then
        write(*,*)'changing for',i,v(i,3)
        V(I,3)=0E0
        laytyp(i)=1
      endif
      IF (LSOLL.EQ.0.AND.LAYTYP(I).EQ.3) THEN
      LSOLF=I
      END IF
      IF (LAYTYP(I).EQ.3) THEN
      LSOLL=I
      END IF
 10   CONTINUE
C
C     NUMBER OF SOLID LAYERS
C
      NSOL=LSOLL-LSOLF+1
      DO 20 I=1,NUML            
        DO 20 J=1,4               
          IPS(I,J)=0
          SS(I,J)=CNUL              
          R(I,J)=CNUL               
          DO 20 K=1,4               
            AUP(I,J,K)=CNUL           
            ALO(I,J,K)=CNUL           
 20   CONTINUE  
      NEQ=0     
      IF (LAYTYP(1).LT.0) GO TO 100            
      IF (LAYTYP(1).EQ.3.OR.LAYTYP(1).EQ.4) GO TO 50             
C     LIQUID HALF SPACE         
      IPS(1,3)=1
      NEQ=1     
      GO TO 100 
C     SOLID HALF SPACE          
 50   IPS(1,3)=1
      IPS(1,4)=2
      NEQ=2     
 100  IF (NUMI.EQ.1) GO TO 200  
      DO 190 I=2,NUMI           
        IF (LAYTYP(I).EQ.3.OR.LAYTYP(I).EQ.4) GO TO 150            
C       LIQUID LAYER              
        IPS(I,1)=NEQ+1            
        IPS(I,3)=NEQ+2            
        NEQ=NEQ+2 
        GO TO 190 
 150    CONTINUE  
C     SOLID LAYER               
        DO 160 J=1,4              
 160      IPS(I,J)=NEQ+J            
        NEQ=NEQ+4 
 190  CONTINUE  
 200  IF (LAYTYP(NUML).LT.0) GO TO 300         
      IF (LAYTYP(NUML).EQ.3.OR.LAYTYP(NUML).EQ.4) GO TO 250         
      IPS(NUML,1)=NEQ+1         
      NEQ=NEQ+1 
      GO TO 300 
 250  IPS(NUML,1)=NEQ+1         
      IPS(NUML,2)=NEQ+2         
      NEQ=NEQ+2 
 300  CALL DETBW
      IBW=MIN0(NEQ-1,IBW)
      CALL DETPNT
      RCC=1.0/V(LAYS((LS-1)/2+1),6)
c      write(*,*)'ls',ls
c      write(*,*)'V(LAYS((LS-1)/2+1),6)',V(LAYS((LS-1)/2+1),6)
c      write(*,*)'numl',numl
c      write(*,*)'v(i,6):',(v(i,6),i=1,5)
c      write(*,*)'rcc',rcc

      CALL VSMUL(V(1,6),1,RCC,RCON1,1,NUML)

      DO 340 I=1,LS
 340    RLIND(I)=FLOAT(LAYS(I)*2-1)
      CALL VNEG(ZUS,1,RZUBUF,1,LS)
      CALL VNEG(ZLS,1,RZLBUF,1,LS)
C
C     DETERMINE SOURCE AND RECEIVER POINTERS 
C
      DO 350 I=1,NUML
        NOSOU(I)=0
        IFSOU(I)=NRD
        ILSOU(I)=1
 350  CONTINUE
      do 355 I=1,4
        NUMTS(I)=0
        NUMTR(I)=0
 355  CONTINUE
      DO 360 I=1,LS
        LL=LAYS(I)
        NOSOU(LL)=NOSOU(LL)+1
        IFSOU(LL)=MIN(IFSOU(LL),I)
        ILSOU(LL)=MAX(ILSOU(LL),I)
        IF (NOSOU(LL).NE.(ILSOU(LL)-IFSOU(LL)+1)) THEN
          WRITE(6,*) '*** PINIT1: Source no.',I,' out of order ***'
          STOP
        END IF
        LT=LAYTYP(LL)
        NUMTS(LT)=NUMTS(LT)+1
        NSPNT(NUMTS(LT),LT)=I
 360  CONTINUE
      DO 365 I=1,IR
        LL=LAY(I) 
        LT=LAYTYP(LL)
        NUMTR(LT)=NUMTR(LT)+1
        NRPNT(NUMTR(LT),LT)=I
 365  CONTINUE
      RETURN    
      END       
      SUBROUTINE PINIT2         
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'
C
C     AIRY FUNCTION SOLUTION ADDED 840907
C     SPEED AT LOWER INTERFACE OF LAYER I IS
C     GIVEN BY -V(I,3) (REPLACING SHEAR).
C
C     COMPLEX CONTOUR INTEGRATION ADDED AS OPTION 'J'  850321
C
      OFFIMA=FREQ*OFFDB/(8.68588964*V(LAYS((LS-1)/2+1),2))
C
      DISNRM=1E0/DSQ
      DO 10 I=1,NUML
        THICK(I)=0E0            
        IF (I.GT.1.AND.I.LT.NUML) THICK(I)=V(I+1,1)-V(I,1)        
        IF (LAYTYP(I).LT.0.OR.LAYTYP(I).EQ.4) GO TO 10             
        IF (LAYTYP(I).LE.2.AND.V(I,4).LE.0) THEN
          FK=FREQ/1000.
          VV=V(I,2)*FK*1E-6*(.007+.264/(2.89+FK**2))
        ELSE
          VV=V(I,4)
        END IF
        VI4=FREQ*VV/(8.68588964*V(I,2))           
c        write(*,*)'v(i,2)',v(i,2),freq,vi4
        AK(I,1)=DSQ/V(I,2)+CMPLX(0E0,-VI4)            
        AK2(I,1)=AK(I,1)*AK(I,1)  
c        write(*,*)'before for laytyp2',ak2(i,1)
        ALAME(I,1)=CSQ*V(I,6)/AK2(I,1)            
C
C     AIRY FUNCTION SOLUTION ADDED 840907
C     SPEED AT LOWER INTERFACE OF LAYER I IS
C     GIVEN BY -V(I,3) (REPLACING SHEAR).
C
       IF (LAYTYP(I).EQ.2) THEN
c        write(*,*)'entering for laytyp2',ak(i,1)
        DELTAK=AIMAG(AK(I,1))/REAL(AK(I,1))
        AKU=REAL(AK(I,1))
        AKL=-DSQ/V(I,3)
        AKU2=AKU*AKU*(1E0-DELTAK*DELTAK)
        AKL2=AKL*AKL*(1E0-DELTAK*DELTAK)
        GRAD=(AKL2-AKU2)/THICK(I)
        AA=SIGN(ABS(GRAD)**0.333333333333333,GRAD)
        ACO(I)=AA*CMPLX(1E0,DELTAK*0.666666666666667)
        AAM2=1E0/(AA*AA)
        CCO(I)=AAM2*CMPLX(1E0,-DELTAK*1.33333333333333)
        BCO(I)=CCO(I)*AKU2*CMPLX(1E0,2*DELTAK)
C     WRITE(6,876) I,ACO(I),BCO(I),CCO(I)
C876  FORMAT(1H ,I3,3(/1H ,2G20.8))
      END IF
      IF (V(I,3).LT.1E-10) GO TO 10             
        VI5=FREQ*V(I,5)/(8.68588964*V(I,3))           
      AK(I,2)=DSQ/V(I,3)+CMPLX(0E0,-VI5)            
      AK2(I,2)=AK(I,2)*AK(I,2)  
      ALAME(I,2)=CSQ*V(I,6)/AK2(I,2)            
      ALAME(I,1)=ALAME(I,1)-2*ALAME(I,2)        
 10   CONTINUE  
C
C     DETERMINE INTEGRATION SAMPLING FOR N-K SCATTERING
      DO 12 I=2,NUML
      IF (LAYTYP(I-1).LE.0) THEN
        DLWNK(I)=REAL(AK(I,1))/NSIP
      ELSE IF (LAYTYP(I).LE.0) THEN
        DLWNK(I)=REAL(AK(I-1,1))/NSIP
      ELSE IF (LAYTYP(I).NE.4.AND.LAYTYP(I-1).NE.4) THEN
       IF (LAYS((LS-1)/2+1).LT.I) THEN
        DLWNK(I)=REAL(AK(I-1,1))/NSIP
       ELSE
        DLWNK(I)=REAL(AK(I,1))/NSIP
       END IF
      ELSE
        DLWNK(I)=1E10
      END IF
      IMX(I)=8E0/(CLEN(I)*DLWNK(I))
      IF (IMX(I).GE.1.AND.IMX(I).LT.IMXMIN) THEN
        IMX(I)=IMXMIN
        DLWNK(I)=8E0/(CLEN(I)*IMX(I))
      END IF
C      WRITE(6,*) I,DLWNK(I),IMX(I)
 12   CONTINUE
      PCORR=CSQ*V(LAYS((LS-1)/2+1),6)         
      PCORR=1E0/PCORR
      DO 13 I=1,NUMT(1)
        LL=LAYT(I,1)
        CON1(LL)=PCORR*V(LL,6)*CSQ
 13   CONTINUE
      DO 14 I=1,NUMT(2)
        LL=LAYT(I,2)
        CON1(LL)=PCORR*V(LL,6)*CSQ
 14   CONTINUE
      DO 15 I=1,NUMT(3)
        LL=LAYT(I,3)
        CON1(LL)=PCORR*ALAME(LL,2)
 15     CONTINUE
      IF (NTISOL.GT.0) THEN
        CALL RWDBUF(41)
        DO 20 I=1,NUMT(4)
        LL=LAYT(I,4)
        RCON1(LL)=DSQ
        CON1(LL)=-AI*CSQ*PCORR
 20     CONTINUE
      END IF
      RETURN    
      END       
      SUBROUTINE LINARR(LSL,SD,DELTA,THETA,LTY)
      INCLUDE 'compar.f'
      INCLUDE 'comnrd.f'
C
      DO 10 I=1,LSL  
 10   SDC(I)=SD+(I-1-(LSL-1)/2.0)*ABS(DELTA)  
      RETURN        
      END           
      SUBROUTINE PHASES(LSL,FREQL,V,DELTA,THETA,LTYP,FOCDEP)  
      INCLUDE 'compar.f'
      INCLUDE 'comnrd.f'
      DIMENSION V(NLA,6)  
      PI=4.*ATAN(1.)
      VMEAN=0       
      IF (SHEAR.AND.V(LAYS((LS-1)/2+1),3).GT.1E-10) THEN
        INDV=3
      ELSE 
        INDV=2
      END IF

      DO 10 I=1,LS  
 10   VMEAN=VMEAN+V(LAYS(I),INDV)/LS       
      ANG=2*PI*FREQ*ABS(DELTA)*SIN(THETA)/VMEAN                  
c      write(*,*)'phases,THETA,DELTA,VMEAN,FREQ,PI,ang'
c      write(*,*)'phases',THETA,DELTA,VMEAN,FREQ,PI,ang
      HMEAN=.5*(SDC(LS)+SDC(1))
      ALEN=ABS(SDC(LS)-SDC(1))
      IF (LS.EQ.1) GO TO 15
      GO TO (15,25,35,45,55),LTYP                
 15   DO 20 I=1,LS  
 20   CPHFAC(I)=CEXP(CMPLX(0.,-(I-1-(LS-1)/2.0)*ANG))/LS       
c      write(*,*)'phases',  CPHFAC(1),ls,ang
      RETURN        
 25   DO 30 I=1,LS  
      CPHFAC(I)=CEXP(CMPLX(0.,-(I-1-(LS-1)/2.0)*ANG))       
      FAC=1.0-COS((I-1)*2.0*PI/FLOAT(LS-1))                 
      CPHFAC(I)=CPHFAC(I)*FAC/LS           
 30   CONTINUE      
      RETURN        
 35   CONTINUE
      IF (ABS(THETA).GT.1E-8) THEN
      HH=FOCDEP-HMEAN
      DM=HH/SIN(THETA)
      R2=(HH/TAN(THETA))**2
      ELSE
      DM=FOCDEP
      R2=DM*DM
      END IF
      ALS=1E0/LS
      DO 40 I=1,LS
      IF (ABS(THETA).GT.1E-8) THEN
      H=FOCDEP-SDC(I)
      ELSE
      H=HMEAN-SDC(I)
      END IF
      D=SQRT(H*H+R2)
      ANG=(D-DM)*2*PI*FREQ/VMEAN
      CPHFAC(I)=CEXP(CMPLX(0.,ANG))
      FAC=1.0-COS((I-1)*2.0*PI/FLOAT(LS-1))
      CPHFAC(I)=CPHFAC(I)*FAC*ALS
 40   CONTINUE
      RETURN
 45   CONTINUE
      FSUM=0.
      DO 48 I=1,LS  
      CPHFAC(I)=CEXP(CMPLX(0.,-(I-1-(LS-1)/2.0)*ANG))       
      FAC=EXP(-(4.*(SDC(I)-HMEAN)/ALEN)**2)
      FSUM=FSUM+FAC
      CPHFAC(I)=CPHFAC(I)*FAC           
 48   CONTINUE      
      DO 49 I=1,LS
 49   CPHFAC(I)=CPHFAC(I)/FSUM
      RETURN        
 55   CONTINUE
      IF (ABS(THETA).GT.1E-8) THEN
      HH=FOCDEP-HMEAN
      DM=HH/SIN(THETA)
      R2=(HH/TAN(THETA))**2
      ELSE
      DM=FOCDEP
      R2=DM*DM
      END IF
      FSUM=0.
      DO 58 I=1,LS
      IF (ABS(THETA).GT.1E-8) THEN
      H=FOCDEP-SDC(I)
      ELSE
      H=HMEAN-SDC(I)
      END IF
      D=SQRT(H*H+R2)
      ANG=(D-DM)*2*PI*FREQ/VMEAN
      CPHFAC(I)=CEXP(CMPLX(0.,ANG))
      FAC=EXP(-(4.*(SDC(I)-HMEAN)/ALEN)**2)
      FSUM=FSUM+FAC
      CPHFAC(I)=CPHFAC(I)*FAC
 58   CONTINUE
      DO 59 I=1,LS
 59   CPHFAC(I)=CPHFAC(I)/FSUM
      RETURN
      END 
c****************************************************************          
      SUBROUTINE CALINT    
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'
c      INCLUDE '[.-]cominv.f'
      COMPLEX FACSQ
      NS=ICUT2-ICUT1+1
C *** OPEN SCRATCH FILE FOR KERNELS
       ierr=0
c       if ((iopt(1).eq.2)) then
c      IF (DECOMP) THEN
c       NGFILS=5
c       IF (ISIZE.LT.IR*15) STOP '>>> PARAMETER ISIZE TOO SMALL <<<'
c      ELSE
c       NGFILS=1
c       IF (ISIZE.LT.IR*3) STOP '>>> PARAMETER ISIZE TOO SMALL <<<'
c      END IF
c      LUOFF=LUGRN-1 
c      DO 5 IFC=1,NGFILS
c       CALL OPNBUF(LUOFF+IFC,2*IR,NS*NOUT,8000/NGFILS)
c 5    CONTINUE
c      IF (SCTOUT) THEN
C *** OPEN FILE AND WRITE HEADER
c        CALL OPFILW(45,IOER)
c        WRITE(45,*) FREQ,LAYS(1)
c      END IF
c      endif

C *** WAVENUMBER RAMP
      CALL VRAMP(WK0,DLWVNO,FAC,1,NWVNO)
      NGVALS=2*IR*3
C *** WAVENUMBER LOOP
      DO 20 II=ICUT1,ICUT2         
       WVNO=CMPLX(FAC(II),OFFIMA)
       IF (ICDR.EQ.0) THEN
         FACSQ=CSQRT(WVNO)
       else
         facsq=1  
       END IF
c      write(*,*)'calint, waveno',wvno,icdr,wk0,dlwvno
         facsqrt(ii)=facsq
       CALL INITS
       CALL BUILD                  
       CALL SOLVE    
       IF (IERR.GT.0) RETURN
C       IF (DEBUG) WRITE(6,*) 'CALLING KERNEL'
c       IF (DECOMP) THEN
c        CALL KERDEC(CFILE,IR)
c       ELSE
        CALL KERNEL(CFILE,IR)
c       END IF
c        write(*,*)'calint:cfile',cfile(1)
cpg       if (iopt(1).eq.1) then

c *** tapering
      IF (II.LT.ICW1) THEN
       TFAC=(0.5*(1E0+COS((II-ICW1)*PI/(ICUT1-ICW1-1))))**4
       CALL VSMUL(CFILE(1),1,TFAC,CFILE(1),1,NGVALS)
      ELSE IF (II.GT.ICW2) THEN
       TFAC=(0.5*(1E0+COS((II-ICW2)*PI/(ICUT2-ICW2+1))))**4
       CALL VSMUL(CFILE(1),1,TFAC,CFILE(1),1,NGVALS)
      else
       tfac=1e0
      END IF
c      if (mod(ii,100).eq.0) write(*,*)ii,tfac,icw1,icw2,ngvals
        i0=0
       disp_to_vel= dsq*1000*1500*1e-6*1000  ! dsq*dens*Cref*um* (1000 ?)
c      write(*,*)'iout',iout(1) ,iout(2) ,iout(3), disp_to_vel,dsq
      IF (IOUT(1).GT.0) THEN
          INDXCF=IR*3*(0)
          do i=1,ir
             wavenoint(ii,i+i0)=cfile(indxcf+i)
          enddo
          i0=i0+ir
       endif 
       IF (IOUT(2).GT.0) THEN
          INDXCF=IR*(1)
          do i=1,ir
             wavenoint(ii,i+i0)=cfile(indxcf+i)* disp_to_vel
c          write(*,*)'cfile(indxcf+i)',cfile(indxcf+i)
          enddo
          i0=i0+ir
       endif 
       IF (IOUT(3).GT.0) THEN
          INDXCF=IR*(2)
          do i=1,ir
             wavenoint(ii,i+i0)=cfile(indxcf+i)* disp_to_vel
c          write(*,*)'cfile(indxcf+i)',cfile(indxcf+i)
          enddo
       endif 
c       else
c          write(*,*)'calint:unknown parameter'
c       write(*,*)'i0',i0
             
cpg         endif
c       endif
c         write(*,*)'cfile',ir,cfile(1)
c         write(*,*)'cfile,facsq',ir,cfile(1)*facsq,facsq
c
c       if ((iopt(1).eq.2)) then
c         DO 10 IFC=1,NGFILS
c          IF (IOUT(1).GT.0) THEN
c          INDXCF=1+IR*3*(IFC-1)
c          IF (ICDR.EQ.0) THEN
c            CALL CVMUL(CFILE(INDXCF),2,FACSQ,0,CFILE(INDXCF),2,IR,1)
c          END IF
c          CALL WRBUF(LUOFF+IFC,CFILE(INDXCF),2*IR)
c        END IF
c        IF (IOUT(2).GT.0) THEN
c         INDXCF=1+IR*(1+3*(IFC-1))
c         IF (ICDR.EQ.0) THEN
c           CALL CVMUL(CFILE(INDXCF),2,FACSQ,0,CFILE(INDXCF),2,IR,1)
c         END IF
c         CALL WRBUF(LUOFF+IFC,CFILE(INDXCF),2*IR)
c        END IF
c        IF (IOUT(3).GT.0) THEN
c         INDXCF=1+IR*(2+3*(IFC-1))
c         IF (ICDR.EQ.0) THEN
c          CALL CVMUL(CFILE(INDXCF),2,FACSQ,0,CFILE(INDXCF),2,IR,1)
c         END IF
c         CALL WRBUF(LUOFF+IFC,CFILE(INDXCF),2*IR)
c        END IF 
c 10    CONTINUE
c       endif
C
C *** OUTPUT ROUGH SURFACE DISCONTINUITIES FOR CALCULATING SCATTERED
C     FIELD
C
c      IF (SCTOUT) CALL SCTRHS(CMPLX(1E0,0E0))
  20   CONTINUE     
c      DO 30 IFC=1,NGFILS
c      CALL ENFBUF(LUOFF+IFC)
c 30   CONTINUE
C                   
      RETURN        
C                   
      END           
      SUBROUTINE CHKSOL
      INCLUDE 'compar.f'

          IF (IERR.NE.0) THEN
C            CALL PREQV(NUML,NUMI)
            WRITE(6,9987) IERR
9987        FORMAT(//1H ,'**** EXECUTION TERMINATED ***',
     -            //1H ,'**** ERROR NUMBER : ',I3)
            STOP
         END IF
         IF (IOERR.NE.0) THEN
            WRITE(6,9887) IOERR
9887        FORMAT(//1H ,'**** EXECUTION TERMINATED ****',
     -            //1H ,'**** IO- ERROR NUMBER : ',I3)
            STOP
         END IF
       RETURN
       END
      SUBROUTINE GETKNL(NREC)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      IF (DEBUG) WRITE(6,*) 'ENTER GETKNL'
C
C    READ KERNELS FROM SCRATCH FILE
C
      IF (DEBUG) WRITE(6,*) 'RWDBUF'
      CALL RWDBUF(LUTGRN)
      IF (DEBUG) WRITE(6,*) 'EXIT RWDBUF'
      NN=2*NP
      DO 2 I=1,3
       IF (IOUT(I).GT.0) THEN
        DO 1 JR=1,ICUT1-1
 1      CFF(JR,I)=CNUL
        DO 11 JR=ICUT2+1,NWVNO
 11     CFF(JR,I)=CNUL
       END IF
 2    CONTINUE
      LREC=2*IR
      DO 4 JR=ICUT1,ICUT2
      DO 3 I=1,3
      IF (IOUT(I).GT.0) THEN
       CALL RDBUF(LUTGRN,CFILE,LREC)
       IF (DEBUG) THEN 
        WRITE(6,*) JR
        WRITE(6,*) (CFILE(LL),LL=1,IR) 
       END IF
       CFF(JR,I)=CFILE(NREC)
      END IF
 3    CONTINUE
 4    CONTINUE
C
C
      RETURN
      END
      SUBROUTINE SCTRHS(cfac)              
C
C     WRITES SCATTERING RIGHT HAND SIDES TO FILE 45
c
      PARAMETER (SQ2PI=2.506628,ONEO2PI=1.591549E-1)
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'
C
C     LOCAL COEFFICIENT MATRICES RE-DEFINED
C
      COMPLEX DBDZ(4),BLC(4)
      COMPLEX S,cfac
      S=WVNO
      DO 200 I=1,NUMI
       DO 200 J=1,4
       ROIN(I,J) = CMPLX(0E0,0E0)
       RUIN(I,J) = CMPLX(0E0,0E0)
 200  CONTINUE
       IF (NUMT(1).GT.0) CALL LIQLAY
       IF (NUMT(2).GT.0) CALL AIRLAY
       IF (NUMT(3).GT.0) CALL SOLLAY
       IF (NUMT(4).GT.0) CALL TISOLL
      DO 10 IN=1,NUMI
      IN1=IN+1      
      IF (ROUGH2(IN1).LT.1E-10) GO TO 10

C     HOMEGENEOUS PART

      DO 252 JJ=1,4
      DBDZ(JJ)   =   -(-ALFA(IN))*ALO(IN,JJ,1)*SS(IN,1)
     1              -(-BETA(IN))*ALO(IN,JJ,2)*SS(IN,2)
     2              -ALFA(IN)*ALO(IN,JJ,3)*SS(IN,3)
     3              -BETA(IN)*ALO(IN,JJ,4)*SS(IN,4)
      DBDZ(JJ)   = DBDZ(JJ) - (-ALFA(IN1))*AUP(IN1,JJ,1)*SS(IN1,1)
     2           - (-BETA(IN1))*AUP(IN1,JJ,2)*SS(IN1,2)
     3           - ALFA(IN1)*AUP(IN1,JJ,3)*SS(IN1,3)
     4           - BETA(IN1)*AUP(IN1,JJ,4)*SS(IN1,4)
 252  CONTINUE

C     SOURCE TERMS SIMMILARLY (ONLY PROP TOWARDS INTERFACE)

      DO 256 JJ=1,4
      DBDZ(JJ) = DBDZ(JJ) + (-ALFA(IN))*ROIN(IN,JJ)
      DBDZ(JJ) = DBDZ(JJ) + ALFA(IN1)*RUIN(IN,JJ)
 256  CONTINUE     
C
C     ROTATION TERMS
C
        DO 2502 JJ=1,4
 2502     BLC(JJ)=CMPLX(0E0,0E0)
        DO 2504 JJ=1,4
          BLC(1)= BLC(1)-(-AI*ALO(IN,2,JJ))*SS(IN,JJ)
          BLC(2)= BLC(2)-(-AI*ALO(IN,1,JJ))*SS(IN,JJ)
          BLC(3)= BLC(3)-(-2E0*AI*ALO(IN,4,JJ))*SS(IN,JJ)
          BLC(1)= BLC(1)-(-AI*AUP(IN1,2,JJ))*SS(IN1,JJ)
          BLC(2)= BLC(2)-(-AI*AUP(IN1,1,JJ))*SS(IN1,JJ)
          BLC(3)= BLC(3)-(-2E0*AI*AUP(IN1,4,JJ))*SS(IN1,JJ)
 2504   CONTINUE

C     rotation for shear stress

        IF (LAYTYP(IN).EQ.3) THEN 
         BLC(4)=-2E0*(-AI)*CON1(IN)*SS(IN,1)*
     &          (-ALFA(IN)*ALO(IN,1,1)-WVNO*ALO(IN,2,1))
         BLC(4)=BLC(4)-2E0*(-AI)*CON1(IN)*SS(IN,2)*
     &          (-BETA(IN)*ALO(IN,1,2)-WVNO*ALO(IN,2,2))
         BLC(4)=BLC(4)-2E0*(-AI)*CON1(IN)*SS(IN,3)*
     &          (ALFA(IN)*ALO(IN,1,3)-WVNO*ALO(IN,2,3))
         BLC(4)=BLC(4)-2E0*(-AI)*CON1(IN)*SS(IN,4)*
     &          (BETA(IN)*ALO(IN,1,4)-WVNO*ALO(IN,2,4))
        END IF
        IF (LAYTYP(IN1).EQ.3) THEN
         BLC(4)=BLC(4)-2E0*(-AI)*CON1(IN1)*SS(IN1,1)*
     &          (-ALFA(IN1)*AUP(IN1,1,1)-WVNO*AUP(IN1,2,1))
         BLC(4)=BLC(4)-2E0*(-AI)*CON1(IN1)*SS(IN1,2)*
     &          (-BETA(IN1)*AUP(IN1,1,2)-WVNO*AUP(IN1,2,2))
         BLC(4)=BLC(4)-2E0*(-AI)*CON1(IN1)*SS(IN1,3)*
     &          (ALFA(IN1)*AUP(IN1,1,3)-WVNO*AUP(IN1,2,3))
         BLC(4)=BLC(4)-2E0*(-AI)*CON1(IN1)*SS(IN1,4)*
     &          (BETA(IN1)*AUP(IN1,1,4)-WVNO*AUP(IN1,2,4))
        END IF
C       SOURCE TERMS

         BLC(1)=BLC(1)+(-AI*ROIN(IN,2))
         BLC(2)=BLC(2)+(-AI*ROIN(IN,1))
         BLC(3)=BLC(3)+(-2E0*AI*ROIN(IN,4))
        IF (LAYTYP(IN).EQ.3) THEN
         BLC(4)=BLC(4)+2*(-AI)*CON1(IN)*
     &          (-ALFA(IN)*ROIN(IN,1)-WVNO*ROIN(IN,2))
        END IF
         BLC(1)=BLC(1)+(-AI*RUIN(IN,2))
         BLC(2)=BLC(2)+(-AI*RUIN(IN,1))
         BLC(3)=BLC(3)+(-2*AI*RUIN(IN,4))
        IF (LAYTYP(IN1).EQ.3) THEN
         BLC(4)=BLC(4)+2*(-AI)*CON1(IN1)*
     &          (ALFA(IN1)*RUIN(IN,1)-WVNO*RUIN(IN,2))
        END IF
       WRITE(45,*) WVNO,IN1,(DBDZ(JJ)*cfac,JJ=1,4),(BLC(JJ)*cfac,JJ=1,4)
 10   CONTINUE
      RETURN        
      END           
      SUBROUTINE DETPNT      
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
C
C     DETERMINE NUMBER OF ROWS FOR EACH INTERFACE
C
      ICNT=1
      DO 510 IN=1,NUMI      
      IN1=IN+1              
      IF (V(IN,2).GT.1E-10) GO TO 100                 
C 
C     VACUUM UPPER HALFSPACE
C     **********************
      IF (V(IN1,3).GT.1E-10) GO TO 50                 
C     NEXT LAYER OR HALFSPACE IS LIQUID               
C     ---------------------------------               
      NRI(IN)=1
      GO TO 500             
 50   CONTINUE              
C     NEXT LAYER OR HALFSPACE IS SOLID                
C     --------------------------------                
      NRI(IN)=2
      GO TO 500             
 100  IF (V(IN,3).GT.1E-10) GO TO 200                 
C 
C     LAYER OR HALFSPACE OVER INTERFACE IS LIQUID     
C     *******************************************     
      IF (V(IN1,2).GT.1E-10) GO TO 125                
C     LOWER HALFSPACE IS VACUUM   
C     -------------------------   
      NRI(IN)=1
      GO TO 500             
 125  IF (V(IN1,3).GT.1E-10) GO TO 150                
C     NEXT LAYER OR HALFSPACE IS LIQUID               
C     ---------------------------------               
      NRI(IN)=2
      GO TO 500             
C 
C     NEXT LAYER OR HALFSPACE IS SOLID                
C     --------------------------------                
 150  CONTINUE              
      NRI(IN)=3
      GO TO 500             
 200  CONTINUE              
C     LAYER OR HALFSPACE OVER INTERFACE IS SOLID      
C     *******************************************     
C     LOWER HALFSPACE VACUUM
C     ----------------------
      IF (V(IN1,2).GT.1E-10) GO TO 225                
      NRI(IN)=2
      GO TO 500             
 225  IF (V(IN1,3).GT.1E-10) GO TO 250                
C     NEXT LAYER OR HALFSPACE IS LIQUID               
C     ---------------------------------               
      NRI(IN)=3
      GO TO 500             
 250  CONTINUE              
C     NEXT LAYER OR HALFSPACE IS SOLID                
C     --------------------------------                
      NRI(IN)=4
 500  IRST(IN)=ICNT
      ICNT=ICNT+NRI(IN)
 510  CONTINUE              
      NRI(NUML)=0
      ICNT=1
      DO 10 IL=1,NUML
      NCL(IL)=0
      IF (V(IL,2).LT.1E-10) GO TO 10
      IF (V(IL,3).GE.1E-10) GO TO 5
      NCL(IL)=1
      GO TO 8
 5    NCL(IL)=2
 8    IF (IL.NE.1.AND.IL.NE.NUML) NCL(IL)=2*NCL(IL)
      ICST(IL)=ICNT
      ICNT=ICNT+NCL(IL)
 10   CONTINUE
C
C     DETERMINE FIRST POINTER FOR EACH LAYER
C
      ISTART(1)=1
      ISTART(2)=ISTART(1)+NRI(1)*NCL(1)
      IF (NUML.LE.2) GO TO 20
      DO 19 IL=3,NUML
 19   ISTART(IL)=ISTART(IL-1)+(NRI(IL-2)+NRI(IL-1))*NCL(IL-1)
 20   CONTINUE
C
C     DETERMINE POINTERS
C
      IF (NCL(1).EQ.0) GO TO 30
      DO 25 J=1,NCL(1)
      DO 22 L=1,NRI(1)
      IRN(ISTART(1)+(J-1)*NRI(1)+L-1)=IRST(1)+L-1
 22   CONTINUE
      ICP(ICST(1)+J-1)=ISTART(1)+(J-1)*NRI(1)
      IDP(ICST(1)+J-1)=ISTART(1)+(J-1)*(NRI(1)+1)
 25   CONTINUE
 30   DO 40 IL=2,NUML
      INA=IL-1
      INB=IL
      IF (NCL(IL).EQ.0) GO TO 40
      DO 35 J=1,NCL(IL)
      NR=NRI(INA)+NRI(INB)
      DO 32 L=1,NR
      IRN(ISTART(IL)+(J-1)*NR+L-1)=IRST(INA)+L-1
 32   CONTINUE
      ICP(ICST(IL)+J-1)=ISTART(IL)+(J-1)*NR
      IDP(ICST(IL)+J-1)=ISTART(IL)+(ICST(IL)-IRST(INA))+(J-1)*(NR+1)
 35   CONTINUE
 40   CONTINUE
C
C     DETERMINE NUMBER OF SPARSE ELEMENTS
C
      NNA=ISTART(NUML)+NRI(NUMI)*NCL(NUML)-1
      ICP(NEQ+1)=NNA+1
      DO 710 I=1,NEQ
      DO 705 J=ICP(I),ICP(I+1)-1
      IIRR=IRN(J)
      IICC=I-IIRR+IBW+1
      IICC=2*(IIRR+(IICC-1)*NEQ)-1
      INDB(J)=IICC
 705  CONTINUE
 710  CONTINUE
      IRHCOL=3*IBW+2
      NNB=(IRHCOL-1)*NEQ
      CALL DETIND
      RETURN                
      END                   
      SUBROUTINE DETIND              
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
C
      DO 10 IN=1,NUMI
      IN1=IN+1      
      IF (V(IN,2).GT.1E-10) GO TO 100   
C                   
C     VACUUM UPPER HALFSPACE            
C     **********************            
      IF (V(IN1,3).GT.1E-10) GO TO 50   
C     NEXT LAYER OR HALFSPACE IS LIQUID 
C     --------------------------------- 
      I=IPS(IN1,1)  
      J=IPS(IN1,3)  
      INDR(I)=INDX2(IN,3)
      IST=ISTART(IN1)
      INDA(IST)=INDXUP(IN1,3,1)
      IF (IN.GE.NUMI) GO TO 10
      IST=IST+NRI(IN)+NRI(IN1)
      INDA(IST)=INDXUP(IN1,3,3)
      GO TO 10        
 50   CONTINUE      
C     NEXT LAYER OR HALFSPACE IS SOLID  
C     --------------------------------  
      I=IPS(IN1,1)  
      J=IPS(IN1,2)  
      INDR(I)=INDX2(IN,3)
      INDR(J)=INDX2(IN,4)
      IST=ISTART(IN1)
      IF (IN.GE.NUMI) THEN
      INC=NRI(IN)
      J2=2
      ELSE
      INC=NRI(IN)+NRI(IN1)
      J2=4
      END IF
      DO 60 J=1,J2
      DO 55 I=0,1
 55   INDA(IST+I)=INDXUP(IN1,3+I,J)
 60   IST=IST+INC
      GO TO 10        
 100  IF (V(IN,3).GT.1E-10) GO TO 200   
C                   
C     LAYER OR HALFSPACE OVER INTERFACE IS LIQUID           
C     *******************************************           
      IF (V(IN1,2).GT.1E-10) GO TO 125  
C     LOWER HALFSPACE IS VACUUM         
C     -------------------------         
      I=IPS(IN,3)   
      INDR(I)=INDX2(IN,3)
      IF (IN.EQ.1) THEN
      J1=3
      IST=ISTART(1)
      INC=NRI(1)
      ELSE
      J1=1
      IST=ISTART(IN)+NRI(IN-1)
      INC=NRI(IN-1)+NRI(IN)
      END IF
      DO 110 J=J1,3,2
      INDA(IST)=INDXLO(IN,3,J)
 110  IST=IST+INC
      GO TO 10        
 125  IF (V(IN1,3).GT.1E-10) GO TO 150  
C     NEXT LAYER OR HALFSPACE IS LIQUID 
C     --------------------------------- 
      I=IPS(IN,3)   
      J=IPS(IN1,1)  
      INDR(I)=INDX2(IN,1)
      INDR(J)=INDX2(IN,3)
      IF (IN.EQ.1) THEN
      J1=3
      IST=ISTART(1)
      INC=NRI(1)
      ELSE
      J1=1
      IST=ISTART(IN)+NRI(IN-1)
      INC=NRI(IN-1)+NRI(IN)
      END IF
      DO 130 J=J1,3,2
      INDA(IST)=INDXLO(IN,1,J)
      INDA(IST+1)=INDXLO(IN,3,J)
 130  IST=IST+INC
      IF (IN.GE.NUMI) THEN
      J2=1
      IST=ISTART(IN1)
      INC=NRI(IN)
      ELSE
      J2=3
      IST=ISTART(IN1)
      INC=NRI(IN)+NRI(IN1)
      END IF
      DO 140 J=1,J2,2
      INDA(IST)=INDXUP(IN1,1,J)
      INDA(IST+1)=INDXUP(IN1,3,J)
 140  IST=IST+INC
      GO TO 10        
C                   
C     NEXT LAYER OR HALFSPACE IS SOLID  
C     --------------------------------  
 150  CONTINUE      
      I=IPS(IN,3)   
      J=IPS(IN1,1)  
      K=IPS(IN1,2)  
      INDR(I)=INDX2(IN,3)
      INDR(J)=INDX2(IN,1)
      INDR(K)=INDX2(IN,4)
      IF (IN.EQ.1) THEN
      J1=3
      IST=ISTART(1)
      INC=NRI(1)
      ELSE
      J1=1
      IST=ISTART(IN)+NRI(IN-1)
      INC=NRI(IN-1)+NRI(IN)
      END IF
      DO 170 J=J1,3,2
      INDA(IST)=INDXLO(IN,3,J)
      INDA(IST+1)=INDXLO(IN,1,J)
      INDA(IST+2)=INDXLO(IN,4,J)
 170  IST=IST+INC
      IF (IN.GE.NUMI) THEN
      J2=2
      IST=ISTART(IN1)
      INC=NRI(IN)
      ELSE
      J2=4
      IST=ISTART(IN1)
      INC=NRI(IN)+NRI(IN1)
      END IF
      DO 180 J=1,J2
      INDA(IST)=INDXUP(IN1,3,J)
      INDA(IST+1)=INDXUP(IN1,1,J)
      INDA(IST+2)=INDXUP(IN1,4,J)
 180  IST=IST+INC
      GO TO 10        
 200  CONTINUE      
C     LAYER OR HALFSPACE OVER INTERFACE IS SOLID            
C     *******************************************           
C     LOWER HALFSPACE VACUUM            
C     ----------------------            
      IF (V(IN1,2).GT.1E-10) GO TO 225  
      I=IPS(IN,3)   
      J=IPS(IN,4)   
      INDR(I)=INDX2(IN,3)
      INDR(J)=INDX2(IN,4)
      IF (IN.LE.1) THEN
      J1=3
      IST=ISTART(1)
      INC=NRI(1)
      ELSE
      J1=1
      IST=ISTART(IN)+NRI(IN-1)
      INC=NRI(IN-1)+NRI(IN)
      END IF
      DO 220 J=J1,4
      INDA(IST)=INDXLO(IN,3,J)
      INDA(IST+1)=INDXLO(IN,4,J)
 220  IST=IST+INC
      GO TO 10        
 225  IF (V(IN1,3).GT.1E-10) GO TO 250  
C     NEXT LAYER OR HALFSPACE IS LIQUID 
C     --------------------------------- 
      I=IPS(IN,3)   
      J=IPS(IN,4)   
      K=IPS(IN1,1)  
      INDR(I)=INDX2(IN,1)
      INDR(J)=INDX2(IN,3)
      INDR(K)=INDX2(IN,4)
      IF (IN.LE.1) THEN
      J1=3
      IST=ISTART(1)
      INC=NRI(1)
      ELSE
      J1=1
      IST=ISTART(IN)+NRI(IN-1)
      INC=NRI(IN-1)+NRI(IN)
      END IF
      DO 230 J=J1,4
      INDA(IST)=INDXLO(IN,1,J)
      INDA(IST+1)=INDXLO(IN,3,J)
      INDA(IST+2)=INDXLO(IN,4,J)
 230  IST=IST+INC
      IF (IN.GE.NUMI) THEN
      J2=1
      IST=ISTART(IN1)
      INC=NRI(IN)
      ELSE
      J2=3
      IST=ISTART(IN1)
      INC=NRI(IN)+NRI(IN1)
      END IF
      DO 240 J=1,J2,2
      INDA(IST)=INDXUP(IN1,1,J)
      INDA(IST+1)=INDXUP(IN1,3,J)
      INDA(IST+2)=INDXUP(IN1,4,J)
 240  IST=IST+INC
      GO TO 10        
 250  CONTINUE      
C     NEXT LAYER OR HALFSPACE IS SOLID  
C     --------------------------------  
      I=IPS(IN,3)   
      J=IPS(IN,4)   
      K=IPS(IN1,1)  
      L=IPS(IN1,2)  
      INDR(I)=INDX2(IN,1)
      INDR(J)=INDX2(IN,2)
      INDR(K)=INDX2(IN,3)
      INDR(L)=INDX2(IN,4)
      IF (IN.LE.1) THEN
      J1=3
      IST=ISTART(1)
      INC=NRI(1)
      ELSE
      J1=1
      IST=ISTART(IN)+NRI(IN-1)
      INC=NRI(IN-1)+NRI(IN)
      END IF
      DO 260 J=J1,4
      DO 255 I=0,3
 255  INDA(IST+I)=INDXLO(IN,I+1,J)
 260  IST=IST+INC
      IST=ISTART(IN1)
      IF (IN.GE.NUMI) THEN
      J2=2
      INC=NRI(IN)
      ELSE
      J2=4
      INC=NRI(IN)+NRI(IN1)
      END IF
      DO 280 J=1,J2
      DO 275 I=0,3
 275  INDA(IST+I)=INDXUP(IN1,I+1,J)
 280  IST=IST+INC
 10   CONTINUE
      DO 20 I=1,NUML
      DO 20 J=1,4
 20   IF (IPS(I,J).GT.0) INDS(IPS(I,J))=INDX2(I,J)
      RETURN        
      END           
C                   
      INTEGER FUNCTION INDX2(I,J)
      INCLUDE 'compar.f'
      II=I+(J-1)*NLA
      INDX2=2*II-1
      RETURN
      END
      INTEGER FUNCTION INDXLO(I,J,K)
      INCLUDE 'compar.f'
      II=I+(J-1)*NLA+(K-1)*NLA*4
      INDXLO=2*II-1
      RETURN
      END
      INTEGER FUNCTION INDXUP(I,J,K)
      INCLUDE 'compar.f'
      II=I+(J-1)*NLA+(K-1)*NLA*4+NLA*16
      INDXUP=2*II-1
      RETURN
      END
      SUBROUTINE HERMIT(I,FF,NP,LS,ICUT1,ICUT2,DLWVNO,PI,WK0,                   
     *IPLOT1,IPLOT2)
      DIMENSION FF(2,1)
      IPLOT1=ICUT1                
      NCUT=100
      IF(ICUT1.LE.2)   GO TO 4000       
      IPL1=MAX0(2,ICUT1-NCUT)           
      IPL2=ICUT1-1  
C                   
C   A AND B ARE THE LOWEST AND THE HIGHEST WAVE NUMBERS OF INTEREST             
      A=(WK0+(IPL1-1)*DLWVNO)*1.0E3     
      B=1.0E3*(WK0+(ICUT1-1)*DLWVNO)    
      FA=0.0        
      FB=FF(I,ICUT1)
      SIGNFB=SIGN(1.0,FB)               
      F1A=0.0       
      F1B=1.0E-3*(FF(I,ICUT1+1)-FF(I,ICUT1))/DLWVNO         
C   F1B IS NOT ALLOWED TO BE NEGATIVE WHEN FB IS POSITIVE, AND VICEVERSA        
      IF(SIGNFB.GT.0.0)   F1B=AMAX1(F1B,0.)                 
      IF(SIGNFB.LT.0.0)   F1B=AMIN1(F1B,0.)                 
C                   
C   THE FOLLOWING STATEMENT CHECKS THE CONDITION THAT THE HERMITE               
C   FUNCTION DOES NOT CHANGE SIGN IN THE INTERVAL A-B       
      IF(((3.0*FB*SIGNFB)/(B-A)).GE.F1B*SIGNFB)   GO TO 1000
C   WE INCREASE A TO SATISFY THE ABOVE CONDITION            
      A=B-3*FB/F1B  
      IPL1=1.5+(A*1.0E-3-WK0)/DLWVNO    
 1000 CONTINUE      
C                   
      DO 3000   II=IPL1,IPL2            
      X = (WK0 + FLOAT(II-1) * DLWVNO) * 1.0E3              
C   IN THIS CASE P1 IS ALWAYS ZERO (FA=0,F1A=0)             
C     P1=((X-B)/(B-A))**2*((1.0+2.0*(X-A)/(B-A))*FA+(X-A)*F1A)                  
      P2=((X-A)/(B-A))**2*((1.0-2.0*(X-B)/(B-A))*FB+(X-B)*F1B)                  
C     FF(I,II)=P1+P2
      FF(I,II)=P2   
 3000 CONTINUE      
       IPLOT1=IPL1  
C                   
C                   
 4000 CONTINUE      
      IPLOT2=ICUT2  
      IF(ICUT2.EQ.LS)   RETURN          
      IPL1=ICUT2+1  
      IPL2=MIN0(LS,ICUT2+NCUT)          
C                   
C   A AND B ARE THE LOWEST AND THE HIGHEST WAVE NUMBERS OF INTEREST             
      A = (WK0 + FLOAT(ICUT2-1) * DLWVNO) * 1.0E3           
      B=1.0E3*(WK0+(IPL2-1)*DLWVNO)     
      FA=FF(I,ICUT2)
      SIGNFA=SIGN(1.0,FA)               
      FB=0.0        
      F1A=1.0E-3*(FF(I,ICUT2)-FF(I,ICUT2-1))/DLWVNO         
C   F1A IS NOT ALLOWED TO BE POSITIVE   
      IF(SIGNFA.GT.0.0)F1A=AMIN1(F1A,0.)
      IF(SIGNFA.LT.0.0)F1A=AMAX1(F1A,0.)
      F1B=0.0       
C                   
C   THE FOLLOWING STATEMENT CHECKS THE CONDITION THAT THE HERMITE               
C   FUNCTION DOES NOT BECOME NEGATIVE IN THE INTERVAL A-B   
      IF(((3.0*FA*SIGNFA)/(B-A)).GE.-F1A*SIGNFA)   GO TO 5000                   
C   WE DECREASE B TO SATISFY THE ABOVE CONDITION            
      B=(-3.0*FA)/F1A+A                 
      IPL2=1+IFIX((B*1.0E-3-WK0)/DLWVNO)
 5000 CONTINUE      
C                   
      DO 7000   II=IPL1,IPL2            
      X = (WK0 + FLOAT(II-1) * DLWVNO) * 1.0E3              
      P1=((X-B)/(B-A))**2*((1.0+2.0*(X-A)/(B-A))*FA+(X-A)*F1A)                  
C   IN THIS CASE P2 IS ALWAYS ZERO (FB=0,F1B=0)             
C     P2=((X-A)/(B-A))**2*((1.0-2.0*(X-B)/(B-A))*FB+(X-B)*F1B)                  
C     FF(I,II)=P1+P2
      FF(I,II)=P1   
 7000 CONTINUE      
      IPLOT2=IPL2   
C     WRITE(6,100) ICUT1,ICUT2,IPLOT1,IPLOT2                
C     WRITE(6,300)  
      RETURN        
      END           
      SUBROUTINE CHERMIT(CFF,LS,ICUT1,ICUT2,DLWVNO,WK0,   
     *IPLOT1,IPLOT2)
C
C *** COMPLEX HERMITE EXTRAPOLATOR
C     HS 89/10/3
      COMPLEX CFF(LS)
      IPLOT1=ICUT1                
      NCUT=min(100,NINT(LS*.05))
      IF(ICUT1.LE.2)   GO TO 4000       
      IPL1=MAX0(2,ICUT1-NCUT)           
      IPL2=ICUT1-1  
C                   
C   A AND B ARE THE LOWEST AND THE HIGHEST WAVE NUMBERS OF INTEREST             
      A=(WK0+(IPL1-1)*DLWVNO)     
      B=(WK0+(ICUT1-1)*DLWVNO)    
      FA=0.0        
      FB=ABS(CFF(ICUT1))
      F1A=0.0       
      F1B=(ABS(CFF(ICUT1+1))-ABS(CFF(ICUT1)))/DLWVNO         
C *** PHASES
      RR=REAL(CFF(ICUT1))
      RI=AIMAG(CFF(ICUT1))
      PB=ATAN2(RI,RR)
      RR=REAL(CFF(ICUT1+1))
      RI=AIMAG(CFF(ICUT1+1))
      PB2=ATAN2(RI,RR)
      DP=PB2-PB
c *** f1b must be positive
      F1B=AMAX1(0E0,F1B)
C                   
C   THE FOLLOWING STATEMENT CHECKS THE CONDITION THAT THE HERMITE               
C   FUNCTION DOES NOT CHANGE SIGN IN THE INTERVAL A-B       
      IF(((3.0*FB)/(B-A)).GE.F1B)   GO TO 1000
C   WE INCREASE A TO SATISFY THE ABOVE CONDITION            
      A=B-3*FB/F1B  
      IPL1=1.5+(A-WK0)/DLWVNO    
 1000 CONTINUE      
C *** changed 89/10/3 HS
      C0=0
      C1=0
      C2=FB
      C3=(F1B*(B-A)-2E0*FB)                 
C **********************
      DO 3000   II=IPL1,IPL2            
      X = (WK0 + FLOAT(II-1) * DLWVNO)               
C   IN THIS CASE P1 IS ALWAYS ZERO (FA=0,F1A=0)             
C     P1=((X-B)/(B-A))**2*((1.0+2.0*(X-A)/(B-A))*FA+(X-A)*F1A)             
C      P2=((X-A)/(B-A))**2*((1.0-2.0*(X-B)/(B-A))*FB+(X-B)*F1B)            
C     FF(I,II)=P1+P2
      XB=(X-B)/(B-A)
      XA=(X-A)/(B-A)
      P2=XA*XA*(C2+C3*XB)
      CFF(II)=P2*CMPLX(COS(PB+(II-ICUT1)*DP),SIN(PB+(II-ICUT1)*DP))   
 3000 CONTINUE      
       IPLOT1=IPL1  
C                   
C                   
 4000 CONTINUE      
      IPLOT2=ICUT2  
      IF(ICUT2.EQ.LS)   RETURN          
      IPL1=ICUT2+1  
      IPL2=MIN0(LS,ICUT2+NCUT)          
C                   
C   A AND B ARE THE LOWEST AND THE HIGHEST WAVE NUMBERS OF INTEREST             
      A = (WK0 + FLOAT(ICUT2-1) * DLWVNO) 
      B = (WK0+(IPL2-1)*DLWVNO)     
      FA=abs(CFF(ICUT2))
      F1A=(ABS(CFF(ICUT2))-ABS(CFF(ICUT2-1)))/DLWVNO         
      FB=0.0        
      F1B=0
C *** PHASES
      RR=REAL(CFF(ICUT2))
      RI=AIMAG(CFF(ICUT2))
      PA=ATAN2(RI,RR)
      RR=REAL(CFF(ICUT2-1))
      RI=AIMAG(CFF(ICUT2-1))
      PA2=ATAN2(RI,RR)
      DP=PA-PA2
C *** F1A MUST BE NEGATIVE
      F1A=AMIN1(F1A,0E0)
C                   
C   THE FOLLOWING STATEMENT CHECKS THE CONDITION THAT THE HERMITE               
C   FUNCTION DOES NOT BECOME NEGATIVE IN THE INTERVAL A-B   
      IF(((3.0*FA)/(B-A)).GE.-F1A)   GO TO 5000                   
C   WE DECREASE B TO SATISFY THE ABOVE CONDITION            
      B=(-3.0*FA)/F1A+A                 
      IPL2=1+IFIX((B-WK0)/DLWVNO)
 5000 CONTINUE      
C                   
C *** changed 89/10/3 HS
      C0=FA
      C1=F1A*(B-A)
      C2=(-F1A*(B-A)-FA)
      C3=(F1A*(B-A)+2E0*FA)                 
C **********************
      DO 7000   II=IPL1,IPL2            
      X = (WK0 + FLOAT(II-1) * DLWVNO)
C      P1=((X-B)/(B-A))**2*((1.0+2.0*(X-A)/(B-A))*FA+(X-A)*F1A)    
C   IN THIS CASE P2 IS ALWAYS ZERO (FB=0,F1B=0)             
C     P2=((X-A)/(B-A))**2*((1.0-2.0*(X-B)/(B-A))*FB+(X-B)*F1B) 
C     FF(I,II)=P1+P2
      XB=(X-B)/(B-A)
      XA=(X-A)/(B-A)
      P1=C0+XA*(C1+XA*(C2+C3*XB))
      CFF(II)=P1*CMPLX(COS(PA+(II-ICUT2)*DP),SIN(PA+(II-ICUT2)*DP))      
 7000 CONTINUE      
      IPLOT2=IPL2   
C     WRITE(6,100) ICUT1,ICUT2,IPLOT1,IPLOT2                
C     WRITE(6,300)  
      RETURN        
      END           
      SUBROUTINE DETBW                  
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      IBW=0         
      DO 500 IN=1,NUMI                  
      IN1=IN+1      
      IF (V(IN,2).GT.1E-10) GO TO 100   
C                   
C     VACUUM UPPER HALFSPACE            
C     **********************            
      IF (V(IN1,3).GT.1E-10) GO TO 50   
C     NEXT LAYER OR HALFSPACE IS LIQUID 
C     --------------------------------- 
      I=IPS(IN1,1)  
      J=IPS(IN1,3)  
      IF (IN1.LT.NUML) IBW=MAX0(IBW,IABS(I-J))              
      GO TO 500     
 50   CONTINUE      
C     NEXT LAYER OR HALFSPACE IS SOLID  
C     --------------------------------  
      I=IPS(IN1,1)  
      J=IPS(IN1,2)  
      IBW=MAX0(IBW,IABS(I-J))           
      IF (IN1.GE.NUML) GO TO 500        
      K=IPS(IN1,3)  
      L=IPS(IN1,4)  
      IBW=MAX0(IBW,IABS(I-K),IABS(I-L),IABS(J-K),IABS(J-L)) 
      GO TO 500     
 100  IF (V(IN,3).GT.1E-10) GO TO 200   
C                   
C     LAYER OR HALFSPACE OVER INTERFACE IS LIQUID           
C     *******************************************           
      IF (V(IN1,2).GT.1E-10) GO TO 125  
C     LOWER HALFSPACE IS VACUUM         
C     -------------------------         
      I=IPS(IN,3)   
      J=IPS(IN,1)   
      IF (IN.GT.1) IBW=MAX0(IBW,IABS(I-J))                  
      GO TO 500     
 125  IF (V(IN1,3).GT.1E-10) GO TO 150  
C     NEXT LAYER OR HALFSPACE IS LIQUID 
C     --------------------------------- 
      I=IPS(IN,3)   
      J=IPS(IN1,1)  
      IBW=MAX0(IBW,IABS(I-J))           
      IF (IN.GE.NUMI) GO TO 130         
      K=IPS(IN1,3)  
      IBW=MAX0(IBW,IABS(I-K),IABS(J-K)) 
 130  IF (IN.LE.1) GO TO 500            
      K=IPS(IN,1)   
      IBW=MAX0(IBW,IABS(I-K),IABS(J-K)) 
      GO TO 500     
C                   
C     NEXT LAYER OR HALFSPACE IS SOLID  
C     --------------------------------  
 150  CONTINUE      
      I=IPS(IN,3)   
      J=IPS(IN1,1)  
      K=IPS(IN1,2)  
      IBW=MAX0(IBW,IABS(I-J),IABS(I-K),IABS(J-K))           
      IF (IN.GE.NUMI) GO TO 180         
      M=IPS(IN1,3)  
      L=IPS(IN1,4)  
      IBW=MAX0(IBW,IABS(I-M),IABS(I-L),IABS(J-M),IABS(J-L), 
     1   IABS(K-M),IABS(K-L))           
 180  IF (IN.LE.1) GO TO 500            
      M=IPS(IN,1)   
      IBW=MAX0(IBW,IABS(I-M),IABS(J-M),IABS(K-M))           
      GO TO 500     
 200  CONTINUE      
C     LAYER OR HALFSPACE OVER INTERFACE IS SOLID            
C     *******************************************           
C     LOWER HALFSPACE VACUUM            
C     ----------------------            
      IF (V(IN1,2).GT.1E-10) GO TO 225  
      I=IPS(IN,3)   
      J=IPS(IN,4)   
      IBW=MAX0(IBW,IABS(I-J))           
      IF (IN.LE.1) GO TO 500            
      K=IPS(IN,1)   
      L=IPS(IN,2)   
      IBW=MAX0(IBW,IABS(I-K),IABS(I-L),IABS(J-K),IABS(J-L)) 
      GO TO 500     
 225  IF (V(IN1,3).GT.1E-10) GO TO 250  
C     NEXT LAYER OR HALFSPACE IS LIQUID 
C     --------------------------------- 
      I=IPS(IN,3)   
      J=IPS(IN,4)   
      K=IPS(IN1,1)  
      IBW=MAX0(IBW,IABS(I-J),IABS(I-K),IABS(J-K))           
      IF (IN.GE.NUMI) GO TO 230         
      M=IPS(IN1,3)  
      IBW=MAX0(IBW,IABS(I-M),IABS(J-M),IABS(K-M))           
 230  IF (IN.LE.1) GO TO 500            
      M=IPS(IN,1)   
      N=IPS(IN,2)   
      IBW=MAX0(IBW,IABS(I-M),IABS(I-N),IABS(J-M),IABS(J-N), 
     1         IABS(K-M),IABS(K-N))     
      GO TO 500     
 250  CONTINUE      
C     NEXT LAYER OR HALFSPACE IS SOLID  
C     --------------------------------  
      I=IPS(IN,3)   
      J=IPS(IN,4)   
      K=IPS(IN1,1)  
      L=IPS(IN1,2)  
      IBW=MAX0(IBW,IABS(I-J),IABS(I-K),IABS(I-L),IABS(J-K), 
     1         IABS(J-L),IABS(K-L))     
      IF (IN.GE.NUMI) GO TO 280         
      M=IPS(IN1,3)  
      N=IPS(IN1,4)  
      IBW=MAX0(IBW,IABS(I-M),IABS(I-N),IABS(J-M),IABS(J-N), 
     1         IABS(K-M),IABS(K-N),IABS(L-M),IABS(L-N))     
 280  IF (IN.LE.1) GO TO 500            
      M=IPS(IN,1)   
      N=IPS(IN,2)   
      IBW=MAX0(IBW,IABS(I-M),IABS(I-N),IABS(J-M),IABS(J-N), 
     1         IABS(K-M),IABS(K-N),IABS(L-M),IABS(L-N))     
 500  CONTINUE      
      RETURN        
      END           
      SUBROUTINE FLLNEG
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
C
      DO 10 I=1,3
      IF (IOUT(I).NE.1) GO TO 10
      IF (I.LT.3) THEN
      CALL CVMOV(CFF(ICUT1,I),2,CFF(ICUT1-1,I),-2,NWVNO/2)
      ELSE
      CALL CVNEG(CFF(ICUT1,I),2,CFF(ICUT1-1,I),-2,NWVNO/2)
      END IF
 10   CONTINUE
      RETURN
      END

