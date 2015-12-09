      SUBROUTINE INTGRN_PG(xranges,nrec,nflag)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'combes.f'
      COMPLEX FILON
      COMPLEX CSFAC
      LOGICAL NFLAG
      REAL CBR(2,NP),xranges(*)
      EQUIVALENCE (CBUF(1),CBR(1,1))
      REAL FFS(2,NP)
      EQUIVALENCE (CFFS(1),FFS(1,1))
      integer nrec
c
c      write(*,*)' Entering intgrn,inttyp=',inttyp
c
      NS=ICUT2-ICUT1+1
      IF (ICDR.EQ.1) THEN
       TERM=0.
      ELSE
       TERM=PI/4E0
      END IF
      IFN=IFN+1
c      DO 30 JD=1,IR
c      CALL RWDBUF(LUTGRN)
      DO 4 JR=ICUT1,ICUT2
c      DO 3 I=1,3
       i=1
c      IF (IOUT(I).EQ.0) GO TO 3
c      CALL RDBUF(LUTGRN,CFILE,2*IR)
c      CFF(JR,I)=CFILE(JD)
       CFF(JR,I)=wavenoint(jr,NREC)
c       write(*,*)'integrand,jr',jr, CFF(JR,I)
c 3    CONTINUE
 4    CONTINUE
C
C
      IF (NFLAG) THEN
       CALL FLLNEG
       IC1=NWVNO-ICUT2+1
       IC2=ICUT2
      ELSE
       IC1=ICUT1
       IC2=ICUT2
      END IF
C
C
c      write(*,*)'nflag,ic1,ic2',nflag,ic1,ic2,nwvno
      IF (IC1.GT.1.OR.IC2.LT.NWVNO) THEN
       IPL1=NWVNO
       IPL2=0
c       DO 6 I=1,3
        i=1
c       IF (IOUT(I).GT.0) THEN
        if (IC1.GT.1) CALL VCLR(CFF(1,I),1,2*(IC1-1))
        if (IC2.LT.NWVNO) CALL VCLR(CFF(IC2+1,I),1,2*(NWVNO-IC2))
        CALL CHERMIT(CFF(1,I),NWVNO,IC1,IC2,DLWVNO,
     1              WK0,IPL1,IPL2)
c       END IF
 6     CONTINUE
      END IF

C

      WKM=WK0+(NWVNO-1)*DLWVNO
c      write(*,*)' loop 10',nplots,wkm
      DO 10 J=1,NPLOTS
c     RANGEM=R1+(J-1)*DLRAN
      rangem=xranges(j)
c      write(*,*)' ranges,wk0',rangeM,wk0
      RK=RANGEM*WKM
      IF (INTTYP.EQ.2) THEN
c$$$C *** FULL BESSEL FUNCTION INTEGRATION
c$$$        DO 11 II=1,NWVNO
c$$$         WN=WK0+(II-1)*DLWVNO
c$$$         RKL=WN*RANGEM
c$$$         ARG(II)=WN*RINTPL(BF0,0E0,ONODRK,NRKMAX,RKL)
c$$$         FAC(II)=WN*RINTPL(BF1,0E0,ONODRK,NRKMAX,RKL)
c$$$ 11     CONTINUE
c$$$        FCI=DLWVNO
      ELSE IF (RK.GT.1E-3.OR.ICDR.EQ.1) THEN
c         write(*,*)' #1'
        RST=TERM-WK0*RANGEM
        RSTP=-DLWVNO*RANGEM
        CALL VRAMP(RST,RSTP,ARG,1,NWVNO)
        CALL CVEXP(ARG,1,CBUF,2,NWVNO)
        RST=EXP(OFFIMA*RANGEM)
        CALL VSMUL(CBUF,1,RST,CBUF,1,2*NWVNO)
        IF (ICDR.EQ.1) THEN
          FCI=FNI5
        ELSE
          FCI=FNI5/SQRT(RANGEM)
        END IF
      ELSE
c        write(*,*)' #2'
        DO 18 I=1,NWVNO
        CBUF(I)=CSQRT(CMPLX(WK0+(I-1)*DLWVNO,OFFIMA))
 18     CONTINUE
        FCI=DLWVNO
      END IF
c      write(*,*)' intgrn #2'
      ICNT=0
c      DO 20 I=1,3
       i=1
c      IF (IOUT(I).EQ.0) GO TO 20
c$$$      IF (INTTYP.EQ.2) THEN
c$$$        CALL VCLR(CBUF,1,2*NWVNO)
c$$$        IF (I.EQ.3) THEN
c$$$          CALL VNEG(FAC,1,CBR(2,1),2,NWVNO)
c$$$        ELSE
c$$$          CALL VMOV(ARG,1,CBR(1,1),2,NWVNO)
c$$$        END IF
c$$$      END IF
      IF (NWVNO.EQ.1) THEN
        CFILE(J+ICNT*NPLOTS)=CFF(1,I)*CBUF(1)*FCI
c      write(*,*)' CFILE(J+ICNT*NPLOTS) one w#', CFILE(J+ICNT*NPLOTS)

c$$$      ELSE IF (INTTYP.EQ.1.AND.RK.GT.1E-3) THEN
c$$$        CFILE(J+ICNT*NPLOTS)=FILON(CFF(1,I),CBUF(1),ARG(1),
c$$$     &         NWVNO,FCI)
      ELSE IF (I.NE.3.OR.RK.GT.1E-3) THEN
        CALL CVMUL(CFF(1,I),2,CBUF(1),2,CFFS(1),2,NWVNO,1)
        CALL VTRAPZ(FFS(1,1),2,RR,0,NWVNO,FCI)
        CALL VTRAPZ(FFS(2,1),2,RI,0,NWVNO,FCI)
        CFILE(J+ICNT*NPLOTS)=CMPLX(RR,RI)
c        write(*,*)' CFILE(J+ICNT*NPLOTS)', CFILE(J+ICNT*NPLOTS)
      ELSE
        CFILE(J+ICNT*NPLOTS)=CNUL
      END IF
      ICNT=ICNT+1
 20   CONTINUE
 10   CONTINUE
c      CALL WRBUF(LUTTRF,CFILE,2*NOUT*NPLOTS)
      ICNT=0
c      DO 30 I=1,3
          i=1
c         IF (IOUT(I).EQ.0) GO TO 30
         DO  J=1,NPLOTS
c         write(*,*)'i=,=CFILE(J+ICNT*NPLOTS)',i,j,CFILE(J+ICNT*NPLOTS)
            CFF(j,i)=CFILE(J+ICNT*NPLOTS)
         enddo
         ICNT=ICNT+1
 30   CONTINUE
c      write(*,*)' Exit intgrn'
c      CALL CLSBUF(LUTGRN)
      RETURN
      END

      SUBROUTINE TLOSS2      
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comfip.f'
      REAL TEN
      EQUIVALENCE (LF,NUMFR),(NREC,ISPACE)
      DATA TEN /10./
      IF (DEBUG) WRITE(6,*) 'ENTER TLOSS'
      DO 15 I=1,3
       IF (IOUT(I).EQ.0) GO TO 15
c       do j=1,1000,10
c          write(6,*)'j,ccf bef fft',j,cff(j,i)
c       enddo
       CALL CFFT(CFF(1,I),NWVNO,1)
c       do j=1,1000,10
c          write(6,*)'j,aft fft',j,cff(j,i)
c       enddo
 15   CONTINUE
      DO 10 J=1,LF
       RANGEM=R1+(J-1)*DLRAN
       IF (ICDR.EQ.1) THEN
        FAC(J)=FNI5
       ELSE
        FAC(J)=FNI5/SQRT(RANGEM)
       END IF
       FAC(J)=FAC(J)*EXP(OFFIMA*RANGEM)
       ARG(J)=-RANGEM*WK0
 10   CONTINUE
      CALL CVMEXP(ARG,1,FAC,1,CFFS,2,LF)
      DO 20 I=1,3
       IF (IOUT(I).EQ.0) GO TO 20
c       write(*,*)'lf',lf,numfr
c       write(*,*)'cff(1,i):2',cff(1,i),cff(2,i),cff(3,i),cff(4,i) 
       CALL CVMUL(CFF(1,I),2,CFFS,2,CFF(1,I),2,LF,1)
c       do j=1,1000,10
c          write(6,*)'j,ccf tloss',j,cff(j,i)
c       enddo
c       write(*,*)'cff(1,i):1',cff(1,i),cff(2,i),cff(3,i),cff(4,i)
c**** the to next lines are for getting real tl-loss
c       CALL CVMAGS(CFF(1,I),2,CFF(1,I),1,LF)
c       call vsqrt(cff(1,i),1,cff(1,i),1,lf)
c**** end for getting tlloss
c       IF (DEPTAV) THEN
c        LOGNUM=33+I
c        CALL WRBUF(LOGNUM,CFF(1,I),LF)
c       END IF
c       CALL VCLIP(CFF(1,I),1,1E-30,1E30,CFF(1,I),1,LF)
c       CALL VALG10(CFF(1,I),1,CFF(1,I),1,LF)
c       CALL VSMUL(CFF(1,I),1,TEN,CFF(1,I),1,LF)
c       CALL VNEG(CFF(1,I),1,CFF(1,I),1,LF)
 20   CONTINUE
      RETURN
      END
c********************
      SUBROUTINE TLOSS      
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comfip.f'
      REAL TEN
      EQUIVALENCE (LF,NUMFR),(NREC,ISPACE)
      DATA TEN /10./
      IF (DEBUG) WRITE(6,*) 'ENTER TLOSS'
      DO 15 I=1,3
       IF (IOUT(I).EQ.0) GO TO 15
       CALL CFFT(CFF(1,I),NWVNO,1)
 15   CONTINUE
      DO 10 J=1,LF
       RANGEM=R1+(J-1)*DLRAN
       IF (ICDR.EQ.1) THEN
        FAC(J)=FNI5
       ELSE
        FAC(J)=FNI5/SQRT(RANGEM)
       END IF
       FAC(J)=FAC(J)*EXP(OFFIMA*RANGEM)
       ARG(J)=-RANGEM*WK0
 10   CONTINUE
      CALL CVMEXP(ARG,1,FAC,1,CFFS,2,LF)
      DO 20 I=1,3
       IF (IOUT(I).EQ.0) GO TO 20
       CALL CVMUL(CFF(1,I),2,CFFS,2,CFF(1,I),2,LF,1)
       CALL CVMAGS(CFF(1,I),2,CFF(1,I),1,LF)
       IF (DEPTAV) THEN
        LOGNUM=33+I
        CALL WRBUF(LOGNUM,CFF(1,I),LF)
       END IF
       CALL VCLIP(CFF(1,I),1,1E-30,1E30,CFF(1,I),1,LF)
       CALL VALG10(CFF(1,I),1,CFF(1,I),1,LF)
       CALL VSMUL(CFF(1,I),1,TEN,CFF(1,I),1,LF)
       CALL VNEG(CFF(1,I),1,CFF(1,I),1,LF)
 20   CONTINUE
      RETURN
      END
      SUBROUTINE GETTLAV        
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      REAL TEN
      DIMENSION TLAV(NP2,3)
      EQUIVALENCE (TLAV(1,1),CFF(1,1))
      EQUIVALENCE (LF,NUMFR),(NREC,ISPACE)
      DATA TEN /10./
      IF (DEBUG) WRITE(6,*) 'ENTER GETTLAV'
      SFAC=1E0/IR
      DO 20 I=1,3
      IF (IOUT(I).EQ.0) GO TO 20
        CALL VCLR(TLAV(1,I),1,LF)
        LOGNUM=33+I
        CALL RWDBUF(LOGNUM)
        DO 10 JR=1,IR
          CALL RDBUF(LOGNUM,ARG,LF)
          CALL VSMA(ARG,1,SFAC,TLAV(1,I),1,TLAV(1,I),1,LF)
 10     CONTINUE
      CALL VALG10(TLAV(1,I),1,TLAV(1,I),1,LF)
      CALL VSMUL(TLAV(1,I),1,TEN,TLAV(1,I),1,LF)
      CALL VNEG(TLAV(1,I),1,TLAV(1,I),1,LF)
 20   CONTINUE
      RETURN
      END
      SUBROUTINE PHINT(NFLAG)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      LOGICAL NFLAG
      DIMENSION FFS(2,NP)
      EQUIVALENCE (CFFS(1),FFS(1,1))
      EQUIVALENCE (LF,NUMFR),(NREC,ISPACE)
      IF (DEBUG) WRITE(6,*) 'ENTER PHINT'
C
C    READ KERNELS FROM SCRATCH FILE 30
C
      DO 2 I=1,3
       IF (IOUT(I).GT.0) THEN
         CALL VCLR(CFF(1,I),1,2*(ICUT1-1))
         CALL VCLR(CFF(ICUT2+1,I),1,2*(NWVNO-ICUT2))
cpg       END IF
cpg 2    CONTINUE
cpg      IF (DEBUG) WRITE(6,*) 'RWDBUF'
cpg      CALL RWDBUF(30)
cpg      IF (DEBUG) WRITE(6,*) 'EXIT RWDBUF'
      DO 4 JR=ICUT1,ICUT2
cpg      DO 3 I=1,3
cpg      IF (IOUT(I).GT.0) THEN
cpg       CALL RDBUF(30,CFILE,2*IR)
cpg       CFF(JR,I)=CFILE(NREC)
       CFF(JR,I)=wavenoint(jr,NREC)
 4    CONTINUE
      END IF
 2    continue
 3    CONTINUE
C
C
c       do jr=1,10
c       write(6,*)'phint:cff(jr,1)',jr,cff(jr,1)
c       enddo
      IF (NFLAG) THEN
       CALL FLLNEG
       IC1=NWVNO-ICUT2+1
       IC2=ICUT2
      ELSE
       IC1=ICUT1
       IC2=ICUT2
      END IF
C
C
      IF (DEBUG) WRITE(6,*) 'ENTER HERMIT'
      IF (IC1.GT.2.OR.IC2.LT.NWVNO) THEN
      IPL1=NWVNO
      IPL2=0
      DO 6 I=1,3
      IF (IOUT(I).GT.0) THEN
c      do j=1,200
c        write(6,*)'j,ccf',j,cff(j,i)
c      enddo
c      write(6,*)'nflg,ic1,ic2',nflag,ic1,ic2,nwvno,wk0,dlwvno,ipl1,ipl2      
        CALL CHERMIT(CFF(1,I),NWVNO,IC1,IC2,DLWVNO,
     1              WK0,IPL1,IPL2)
c      do j=1,200
c        write(6,*)'j,ccf after',j,cff(j,i)
c      enddo
      END IF
 6    CONTINUE
      END IF
C
C    CALCULATE WAVENUMBER FACTORS (ONLY FOR FIRST DEPTH)
C
      IF (NREC.EQ.1) THEN       
       IF (ICDR.EQ.1) THEN
        TERM=0.
       ELSE
        TERM=PI/4E0
       END IF
       RSTP=-DLWVNO*R1
       CALL VRAMP(TERM,RSTP,ARG,1,NWVNO)
       CALL CVEXP(ARG,1,CBUF,2,NWVNO)
      END IF
      DO 20 I=1,3
       IF (IOUT(I).EQ.0) GO TO 20
       CALL CVMUL(CFF(1,I),2,CBUF(1),2,CFF(1,I),2,NWVNO,1)
 20   CONTINUE
c      do j=1,20
c        write(6,*)'phint exit: j,ccf after',j,cff(j,1)
c      enddo

      RETURN
      END

      SUBROUTINE PLTLOS0(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      YUP,YDOWN,YINC,SD,RD)
C *** TRANSMISSION LOSS VS RANGE     
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      DIMENSION TLAV(NP2,3)
      EQUIVALENCE (TLAV(1,1),CFF(1,1))
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPTTL(3)
      DATA OPTTL /'STLRAN','WTLRAN','UTLRAN'/

      OPTION(1)=PROGNM
      OPTION(2)=OPTTL(INR)
c      write(*,*)'prognm',prognm
c      write(*,*)'option(1)',option(1)
      IF (DEBUG) WRITE(6,*) 'ENTERING PLTLOS'
      PTIT='TRANSMISSION LOSS'
      I=MAX(INR,1)
      NLAB=3
      WRITE(LAB(1),810) FREQ
      WRITE(LAB(2),811) SD
      WRITE(LAB(3),812) RD
      XTXT='Range (km)$'
      GOTO (901,902,903),INR
 901  YTXT='Normal stress (dB//1Pa)$'
      GO TO 904
 902  YTXT='Vert. particle velocity (dB//1m/s)$'
      go to 904
 903  YTXT='Hor. particle velocity (dB//1m/s)$'
 904  CONTINUE
      XTYP='LIN'
      YTYP='LIN'
      XDIV=1
      YDIV=1
      IGRID=0
      NC=1
C *** WRITE PLP FILE
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
      CALL PLTWRI(LF,RMIN,RSTEP,0.,0.,TLAV(1,I),1,TLAV(1,I),1)
C *** FORMATS
 810  FORMAT('Freq:',F7.1,' Hz$')
 811  FORMAT('SD:',F9.1,' m$')
 812  FORMAT('RD:',F9.1,' m$')
      RETURN
      END
      SUBROUTINE PLDAV(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      YUP,YDOWN,YINC,SD)
C *** DEPTH-AVERAGED LOSS VS RANGE
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      DIMENSION TLAV(NP2,3)
      EQUIVALENCE (TLAV(1,1),CFF(1,1))
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPTAV(3)
      DATA OPTAV /'STLDAV','WTLDAV','UTLDAV'/
      IF (DEBUG) WRITE(6,*) 'ENTERING PLDAV'
      OPTION(1)=PROGNM
      OPTION(2)=OPTAV(INR)
      PTIT='DEPTH AVERAGED LOSS'
      I=MAX(INR,1)
      NLAB=2
      WRITE(LAB(1),810) FREQ
      WRITE(LAB(2),811) SD
      XTXT='Range (km)$'
      GOTO (901,902,903),INR
 901  YTXT='Normal stress (dB//1Pa)$'
      GO TO 904
 902  YTXT='Vert. particle velocity (dB//1m/s)$'
      go to 904
 903  YTXT='Hor. particle velocity (dB//1m/s)$'
 904  CONTINUE
      XTYP='LIN'
      YTYP='LIN'
      XDIV=1
      YDIV=1
      IGRID=0
      NC=1
C *** WRITE PLP FILE
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
      CALL PLTWRI(LF,RMIN,RSTEP,0.,0.,TLAV(1,I),1,TLAV(1,I),1)
C *** FORMATS
 810  FORMAT('Freq:',F7.1,' Hz$')
 811  FORMAT('SD:',F9.1,' m$')
      RETURN
      END
      SUBROUTINE PTLDEP(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      YUP,YDOWN,YINC,SD,RD)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      DATA OPT2 /'STLDEP','WTLDEP','UTLDEP'/
      OPTION(1)=PROGNM
      OPTION(2)=OPT2(INR)
      I=MAX(1,INR)
      PTIT='TRANSMISSION LOSS'
 810  FORMAT('Freq:',F7.1,' Hz$')
 811  FORMAT('SD:',F9.1,' m$')
 812  FORMAT('Range:',F6.1,' km$')
      NLAB=3
      WRITE(LAB(1),810) FREQ
      WRITE(LAB(2),811) SD
      WRITE(LAB(3),812) RD
      GOTO (901,902,903),INR
 901  XTXT='Normal stress (dB//1Pa)$'
      GO TO 904
 902  XTXT='Vert. particle velocity (dB//1m/s)$'
      go to 904
 903  XTXT='Hor. particle velocity (dB//1m/s)$'
 904  continue
      YTXT='Depth (m)$'
      XTYP='LIN'
      YTYP='LIN'
      XDIV=1
      YDIV=1
      IGRID=0
      NC=1
C *** WRITE PLP FILE
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
      CALL PLTWRI(LF,0.,0.,RMIN,RSTEP,CFF(1,I),1,CFF(1,I),1)
      RETURN
      END
c****************
      SUBROUTINE CONDRW(TITLE,NPX,NPY,NX,NY,XLEFT,XRIGHT   
     $,XSCALE,XINC,YUP,YDOWN,YSCALE,YINC,ZMIN              
     $,ZMAX,ZSTEP,FREQ,SD,RECUP,RECLO,X1,XL,PX)                 
      DIMENSION SECTOR(28),PX(1)      
      CHARACTER*50 FILENM
      CHARACTER*4 TITLE(20),TITLEX(20),TITLEY(20)      
      DATA X1PL,Y1PL/2.,2.0/,HGTPT,HGTC,LABPT,NDIV,        
     *NARC/0.1,0.14,-3,1,5/,LABC,LWGT/-1,-1/,NSM/0/         
      DATA DUMMY /0./
      DATA TITLEX /'Rang','e (k','m)  ',17*'    '/
      DATA TITLEY /'Dept','h (m',')   ',17*'    '/
C                
C   FORMATS      
 401  FORMAT(1H ,F15.4,3X,'  NUMBER OF DATA POINTS ALONG THE X AXIS')           
 402  FORMAT(1H ,F15.4,3X,'  NUMBER OF DATA POINTS ALONG THE Y AXIS')           
403   FORMAT(1H ,F15.4,3X,'  DIVX ' )
404   FORMAT(1H ,F15.4,3X,'  DIVY ' )
405   FORMAT(1H ,F15.4,3X,'  FLAGRC ' )
406   FORMAT(1H ,F15.4,3X,'  RDUP ' )
407   FORMAT(1H ,F15.4,3X,'  RDLO ' )           
408   FORMAT(1H ,F15.4,3X,'  SOURCE DEPTH (M) ' )          
 409  FORMAT(1H ,F15.4,3X,'  NUMBER OF GRID POINTS ALONG THE X AXIS ' )         
 410  FORMAT(1H ,F15.4,3X,'  NUMBER OF GRID POINTS ALONG THE Y AXIS ' )         
  411 FORMAT(1H ,F15.4,3X,'  FREQUENCY (HZ)' )             
  412 FORMAT(1H ,F15.4,3X,'  DUMMY ' )
  413 FORMAT(1H ,F15.4,3X,'  CAY ' )
  414 FORMAT(1H ,F15.4,3X,'  NRNG ' )
  415 FORMAT(1H ,F15.4,3X,'  ZMIN ' ) 
  416 FORMAT(1H ,F15.4,3X,'  ZMAX ' ) 
  417 FORMAT(1H ,F15.4,3X,'  ZINC ' ) 
  418 FORMAT(1H ,F15.4,3X,'  X ORIGIN OF PLOT IN INCHES ' )
  419 FORMAT(1H ,F15.4,3X,'  DUMMY ' )
  420 FORMAT(1H ,F15.4,3X,'  Y ORIGIN OF PLOT IN INCHES ' )
  421 FORMAT(1H ,F15.4,3X,'  NSM   ' )
  422 FORMAT(1H ,F15.4,3X,'  HGTPT ' )
  423 FORMAT(1H ,F15.4,3X,'  HGTC ' ) 
  424 FORMAT(1H ,F15.4,3X,'  LABPT ' )
  425 FORMAT(1H ,F15.4,3X,'  NDIV ' ) 
  426 FORMAT(1H ,F15.4,3X,'  NARC ' ) 
  427 FORMAT(1H ,F15.4,3X,'  LABC ' ) 
  428 FORMAT(1H ,F15.4,3X,'  LWGT ' ) 
  800 FORMAT('CONDR,FIP,FMT,CPX              ')             
 801  FORMAT(A50)                
  850 FORMAT(20A4)                    
  900 FORMAT(1X,F15.4,3X,'  XLEFT',/,F15.4,4X,'  XRIGHT',/,F15.4,3X,            
     *'   XSCALE',/,F15.4,4X,'  XINC')
  901 FORMAT(1X,F15.4,3X,'  YUP',/,F15.4,4X,'  YDOWN',/,F15.4,3X,               
     *'   YSCALE',/,F15.4,4X,'  YINC')
  950 FORMAT(1H ,F15.4,1X,'    RMIN',/,F15.4,2X,'    RMAX')
      WRITE(28,800)             
      WRITE(28,850)TITLE              
      CALL VCLR(SECTOR,1,28)
      SECTOR(1)=NPX                   
C                
C   SECTOR(4) IS A FLAG WHICH IS SET TO ZERO IN THE RANGE  
C   DEPENDENT VERSION OF SNAP FOR ALL SECTORS EXCEPT THE LAST                   
C   ONE. HERE IS USED TO INDICATE THAT THIS IS THE LAST SECTOR                  
       SECTOR(4)=1.0                  
      WRITE(29,444) (SECTOR(L),L=1,28)
 444  FORMAT(1H ,6G13.5)
      INQUIRE(UNIT=29,NAME=FILENM)
      WRITE(28,801) FILENM             
      IF (ABS(XL-X1).LT.1.0) THEN
      TITLEX(2)='e (m'
      TITLEX(3)=')   '
      DIVX=1E0
      ELSE
      TITLEX(2)='e (k'
      TITLEX(3)='m)  '
      DIVX=1E-3
      END IF
      DIVY=1E0
      CAY=5.
      NRNG=5
      FLAGRC=0.
      WRITE(28,850)TITLEX             
      R1=X1*1.0E3
      R2=XL*1.0E3
      WRITE(28,950)R1,R2              
      AX1=XLEFT*1.0E3                 
      AX2=XRIGHT*1.0E3                
      AX3=XSCALE*1.0E3                
      AX4=XINC*1.0E3                  
      WRITE(28,900)AX1,AX2,AX3,AX4    
      WRITE(28,850)TITLEY             
      WRITE(28,901)YUP,YDOWN,YSCALE,YINC                   
      WRITE(28,401)FLOAT(NPX)
      WRITE(28,402)FLOAT(NPY)
      WRITE(28,403)DIVX              
      WRITE(28,404)DIVY              
      WRITE(28,405)FLAGRC              
      WRITE(28,406)RECUP              
      WRITE(28,407)RECLO                   
      WRITE(28,408)SD     
C   NUMBER OF GRID POINTS ALONG THE X AXIS                 
      WRITE(28,409)FLOAT(NX)          
C   NUMBER OF GRID POINTS ALONG THE Y AXIS                 
      WRITE(28,410)FLOAT(NY)          
      WRITE(28,411)FREQ      
      WRITE(28,412)DUMMY              
      WRITE(28,413)CAY              
      WRITE(28,414)FLOAT(NRNG)              
      WRITE(28,415)ZMIN               
      WRITE(28,416)ZMAX               
      WRITE(28,417)ZSTEP              
C X ORIGIN  OF PLOT IN INCHES         
      WRITE(28,418)X1PL               
      WRITE(28,419)DUMMY              
C Y ORIGIN  OF PLOT IN INCHES         
      WRITE(28,420)Y1PL               
      WRITE(28,421)FLOAT(NSM)                
      WRITE(28,422)HGTPT              
      WRITE(28,423)HGTC               
      WRITE(28,424)FLOAT(LABPT)       
      WRITE(28,425)FLOAT(NDIV)        
      WRITE(28,426)FLOAT(NARC)        
      WRITE(28,427)FLOAT(LABC)        
      WRITE(28,428)FLOAT(LWGT)        
      RETURN     
      END       

c*************************
      SUBROUTINE PLSPECT(DLWVNL,WK0L,SD,RD,TITLE      
     *,XLEN,YLEN,LAY)  
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnla.f'
      INCLUDE 'complo.f'
      LOGICAL NEGASP
      DIMENSION SPEC(NP,2)  
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      equivalence (SPEC(1,1),CFFS(1))
      DATA OPT2 /'SSPECT','WSPECT','USPECT'/
C   
      NEGASP=.FALSE.
      OMR=180./PI
      IPLOT1=ICUT2
      IPLOT2=ICUT1
      DO 2401 II=ICUT1,ICUT2           
      WVNO=WK0L+(II-1)*DLWVNL
      SSS=WVNO/REAL(AK(LAY,1))
C      IF (SSS.LE.1.AND.SSS.GE.0) THEN
      IF (ABS(SSS).LE.1) THEN
       IF (SSS.LE.0E0) NEGASP=.TRUE.
       IPLOT1=MIN0(IPLOT1,II)
       IPLOT2=MAX0(IPLOT2,II)
       SPEC(II,1)=ACOS(SSS)
       CC=SIN(SPEC(II,1))
       SPEC(II,2)=REAL(AK(LAY,1))*CC**2
       SPEC(II,1)=OMR*SPEC(II,1)
      END IF
 2401  CONTINUE
C   
C XAXIS DEFINITION     
C   
      NN=IPLOT2-IPLOT1+1
      XLEFT=0.
      IF (NEGASP) THEN
       XRIGHT=180.
       XINC=30.
      ELSE
       XRIGHT=90.
       XINC=15.
      END IF
      XDIV=1.
      NXDIF=0
C   
      DO 2701 I=1,3    
      IF (IOUT(I).EQ.0) GO TO 2701        
      OPTION(1)=PROGNM
      OPTION(2)=OPT2(I)
C   
C  YAXIS DEFINITION    
C   
      CALL CVMAGS(CFF(IPLOT1,I),2,ARG,1,NN)
      CALL VMUL(ARG,1,SPEC(IPLOT1,2),1,ARG,1,NN)
C      CALL VSQRT(ARG,1,ARG,1,NN)
      CALL VMAX(ARG,1,YMAX,NN)
      YMIN=0.0         
      CALL AUTOAX(YMIN,YMAX,YDOWN,YUP,YINC,YDIV,NYDIF)
      PTIT='ANGULAR SPECTRUM'
 810  FORMAT('Freq:',F7.1,' Hz$')
 811  FORMAT('SD:',F9.1,' m$')
 812  FORMAT('RD:',F9.1,' m$')
      NLAB=3
      WRITE(LAB(1),810) FREQ
      WRITE(LAB(2),811) SD
      WRITE(LAB(3),812) RD
      XTXT='Grazing angle (degrees)$'
      WRITE(YTXT,821) NYDIF
 821  FORMAT('Power (10**',I3,')$')
      XTYP='LIN'
      YTYP='LIN'
      XDIV=1
      IGRID=0
      NC=1
C *** WRITE PLP FILE
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
      CALL PLTWRI(NN,0.,0.,0.,0.,SPEC(IPLOT1,1),1,ARG,1)
 2701 CONTINUE
      RETURN           
      END              

