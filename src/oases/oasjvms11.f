      SUBROUTINE REFLEC(ANGLE1,ANGLE2,IPARM) 
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnp.f'
      INCLUDE 'comnrd.f'
      COMPLEX RFCOEF
      EXTERNAL RFCOEF
      DIMENSION X(NP2,3)
      EQUIVALENCE (CFF(1,1),X(1,1))
C       
c      IF (SCTOUT) THEN
c       CALL OPFILW(45,IOER)
c       WRITE(45,*) FREQ,LAYS(1)
c      END IF
      ANGRA=ANGLE1*PI/180.0
      IF (NWVNO.GT.1) THEN
        DLANGLE=(ANGLE2-ANGLE1)/(NWVNO-1)
        DLANGLE=DLANGLE*PI/180.
      ELSE
        DLANGLE=1
      END IF
      OMR=180./PI

      DO 2400 II=1,NWVNO
        ANG=ANGRA+(II-1)*DLANGLE
        WVNO=REAL(AK(1,1))*COS(ANG)
        CALL INITS
        CALL BUILD  
        CALL SOLVE    
        IF (IERR.GT.0) then
            write(*,*) 'error in reflec'
            stop
        endif
        CFF(II,1)=RFCOEF(IPARM)
c
c *** write scattered field to file 45
c      IF (SCTOUT) CALL SCTRHS(ALFA(1))
       IF (DEBUG) WRITE(6,*) FREQ,ANG,CABS(CFF(II,1))
 2400  CONTINUE     
      CALL CVMAGS(CFF(1,1),2,CFF(1,2),1,NWVNO)
      CALL VALG10(CFF(1,2),1,CFF(1,2),1,NWVNO)
      CALL VNEG(CFF(1,2),1,CFF(1,2),1,NWVNO)
      CALL VSMUL(CFF(1,2),1,10.0,CFF(1,2),1,NWVNO)
      DO 2500 II=1,NWVNO
        X(II,3)=OMR*ATAN2(AIMAG(CFF(II,1)),REAL(CFF(II,1)))  ! angle
 2500 CONTINUE
c      WRITE(30) (X(JJ,2),JJ=1,NWVNO)    ! db scale
c      WRITE(31) (X(JJ,3),JJ=1,NWVNO)    ! angular
      RETURN
      END   
c
c**********************************8
c
      COMPLEX FUNCTION RFCOEF(INT) 
      INCLUDE 'compar.f'
      INCLUDE 'comnla.f'
      INCLUDE 'comnrd.f'
        RFCOEF=ALFA(1)*SS(1,2+INT)/(CPHFAC(1)*EXP(-ZLS(1)*ALFA(1)))
      RETURN
      END   
c
c----
c 
      SUBROUTINE PLRC(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      YUP,YDOWN,YINC,FREQL)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      DIMENSION X(NP2,3)
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      EQUIVALENCE (X(1,1),CFF(1,1))
      DATA OPT2 /' P-P  ',' P-SV ','      '/
      OPTION(1)=PROGNM
      OPTION(2)=OPT2(INR)
      I=MAX(1,INR)
      PTIT='REFLECTION COEFFICIENT'
      NLAB=1
      WRITE(LAB(1),801) FREQL
 801  FORMAT('FREQ:',F7.1,' Hz$')
      XTXT='Grazing angle (deg)$'
      YTXT='Reflection coefficient |R|**2$'
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
      CALL CVMAGS(CFF(1,1),2,ARG,1,LF)
      CALL PLTWRI(LF,RMIN,RSTEP,0.,0.,ARG(1),1,ARG(1),1)
      RETURN
      END
      SUBROUTINE PLREFL(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      YUP,YDOWN,YINC,FREQL)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      DIMENSION X(NP2,3)
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      EQUIVALENCE (X(1,1),CFF(1,1))
      DATA OPT2 /' P-P  ',' P-SV ','      '/
      OPTION(1)=PROGNM
      OPTION(2)=OPT2(INR)
      I=MAX(1,INR)
      PTIT='REFLECTION COEFFICIENT'
      NLAB=1
      WRITE(LAB(1),801) FREQL
 801  FORMAT('FREQ:',F7.1,' Hz$')
      XTXT='Grazing angle (deg)$'
      YTXT='Reflection loss (dB)$'
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
      CALL PLTWRI(LF,RMIN,RSTEP,0.,0.,X(1,2),1,X(1,2),1)
      RETURN
      END
      SUBROUTINE PLREFP(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,FREQL)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      DIMENSION X(NP2,3)
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      EQUIVALENCE (X(1,1),CFF(1,1))
      DATA OPT2 /' P-P  ',' P-SV ','      '/
      OPTION(1)=PROGNM
      OPTION(2)=OPT2(INR)
      I=MAX(1,INR)
      PTIT='REFLECTION COEFFICIENT'
      NLAB=1
      WRITE(LAB(1),801) FREQL
 801  FORMAT('FREQ:',F7.1,' Hz$')
      XTXT='Grazing angle (deg)$'
      YTXT='Phase angle (deg)$'
      XTYP='LIN'
      YTYP='LIN'
      XDIV=1
      ydown=-180
      YUP=180
      YINC=90
      YDIV=1
      IGRID=0
      NC=1
C *** WRITE PLP FILE
      CALL PLPWRI(OPTION,PTIT,TITLE,NLAB,LAB,XLEN,YLEN,
     &                  IGRID,XLEFT,XRIGHT,XINC,XDIV,XTXT,XTYP,
     &                  YDOWN,YUP,YINC,YDIV,YTXT,YTYP,NC)
      CALL PLTWRI(LF,RMIN,RSTEP,0.,0.,X(1,3),1,X(1,3),1)
      RETURN
      END
      SUBROUTINE PLREFF(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      YUP,YDOWN,YINC,ANGLE,IANGLE)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      DIMENSION X(NP2,3)
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      EQUIVALENCE (X(1,1),CFF(1,1))
      DATA OPT2 /' P-P  ',' P-SV ','      '/
      OPTION(1)=PROGNM
      OPTION(2)=OPT2(INR)
      I=MAX(1,INR)
      PTIT='REFLECTION COEFFICIENT'
      NLAB=1
      WRITE(LAB,801) ANGLE
 801  FORMAT('ANGLE:',F5.1,' deg$')
      XTXT='Frequency (Hz)$'
      YTXT='Reflection loss (dB)$'
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
      REWIND(30)
      DO 1350 J=1,LF     
      READ(30) (X(JJ,2),JJ=1,NWVNO)
      X(J,1)=X(IANGLE,2)       
 1350 CONTINUE  
      CALL PLTWRI(LF,RMIN,RSTEP,0.,0.,X(1,1),1,X(1,1),1)
      RETURN
      END
      SUBROUTINE PLRFPH(LF,RMIN,RSTEP,TITLE,INR,
     1      XLEN,YLEN,XLEFT,XRIGHT,XINC,
     2      ANGLE,IANGLE)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      INCLUDE 'complo.f'
      DIMENSION X(NP2,3)
      CHARACTER*80 TITLE
      CHARACTER*6 OPTION(2),OPT2(3)
      EQUIVALENCE (X(1,1),CFF(1,1))
      DATA OPT2 /' P-P  ',' P-SV ','      '/
      OPTION(1)=PROGNM
      OPTION(2)=OPT2(INR)
      I=MAX(1,INR)
      PTIT='REFLECTION COEFFICIENT'
      NLAB=1
      WRITE(LAB,801) ANGLE
 801  FORMAT('ANGLE:',F5.1,' deg$')
      XTXT='Frequency (Hz)$'
      YTXT='Phase angle (deg)$'
      YDOWN=-180.
      YUP=180.
      YINC=90.
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
      REWIND(31)
      DO 1350 J=1,LF     
      READ(31) (X(JJ,3),JJ=1,NWVNO)
      X(J,1)=X(IANGLE,3)       
 1350 CONTINUE  
      CALL PLTWRI(LF,RMIN,RSTEP,0.,0.,X(1,1),1,X(1,1),1)
      RETURN
      END
      SUBROUTINE CONFAW(INR,NPX,NPY,AMN,AMX,XLEFT,XRIGHT,XINC,
     1                  XLEN,FMIN,FMAX,YDOWN,YUP,YINC,YLEN,
     2                  ZMIN,ZMAX,ZINC,TITLE,DLANGLE)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      DIMENSION X(NP2,3)
      CHARACTER*80 TITLE,TITLEX,TITLEY,HEADING
      CHARACTER*3 XTYP,YTYP
      CHARACTER*40 LAB
      CHARACTER*6 OPTION(2),OPT2(3)
      EQUIVALENCE (X(1,1),CFF(1,1))
      DATA OPT2 /'P-P   ','PSV   ','      '/
      DATA XTYP,YTYP /'LIN','LOG'/
      OPTION(1)='CONFR,'
      OPTION(2)=OPT2(INR)
      I=INR
      IF (INR.EQ.0) I=1
      WRITE(28,777) OPTION
      WRITE(28,778) TITLE
 777  FORMAT(2A6,'                                           ')
 778  FORMAT(A80)
 779  FORMAT(A40)
 780  FORMAT(A3)
      WRITE(28,6010) NPX,'NPX'
      WRITE(28,6010) NPY,'NPY'
      WRITE(28,6010) NPX,'NX'
      WRITE(28,6010) NPY,'NY'
      WRITE(28,6030) 1000.*AMN,'MINIMUM ANGLE'
      WRITE(28,6030) 1000.*AMX,'MAXIMUM ANGLE'
      WRITE(28,6030) 1000.*XLEFT,'XLEFT'
      WRITE(28,6030) 1000.*XRIGHT,'XRIGHT'
      XSCALE=ABS(XRIGHT-XLEFT)/XLEN
      WRITE(28,6030) 1000.*XSCALE,'XSCALE'
      WRITE(28,6030) 1000.*XINC,'XINC'
      TITLEX='Grazing angle (degrees)'
      WRITE(28,778) TITLEX
      WRITE(28,780) XTYP
      WRITE(28,6030) FMIN,'MINIMUM FREQUENCY'
      WRITE(28,6030) FMAX,'MAXIMUM FREQUENCY'
      WRITE(28,6030) YDOWN,'YDOWN'
      WRITE(28,6030) YUP,'YUP'
C     FCC=ALOG10(ABS(YUP/YDOWN))/ALOG10(2.0)
C     YSCALE=ABS(YLEN/FCC)
      WRITE(28,6030) YLEN,'NUMBER OF CM PR OCTAVE'
      WRITE(28,6030) YINC,'YINC'
      TITLEY='Frequency (Hz)'
      WRITE(28,778) TITLEY
      WRITE(28,780) YTYP
      WRITE(28,6010) 5,'CAY'
      WRITE(28,6010) 5,'NRNG'
      WRITE(28,6010) 0,'NSM'
      WRITE(28,6030) ZMAX,'ZMIN'
      WRITE(28,6030) ZMIN,'ZMAX'
      WRITE(28,6030) ZINC,'ZINC'
      WRITE(28,6030) 2.0,'X1PL'
      WRITE(28,6030) 2.0,'Y1PL'
      WRITE(28,6030) 0.1,'HGTPT'
      WRITE(28,6030) 0.14,'HGTC'
      WRITE(28,6010) -3,'LABPT'
      WRITE(28,6010) 1,'NDIV'
      WRITE(28,6010) 5,'NARC'
      LABC=-1
      IF (ZINC.LT.1.0) LABC=1
      WRITE(28,6010) LABC,'LABC'
      WRITE(28,6010) -1,'LWGT'
      WRITE(28,6030) 0.,'DUMMY'
      WRITE(28,6030) 0.,'DUMMY'
      WRITE(28,811)
 811  FORMAT('FILENAME')
      WRITE(28,6030) 1000.*DLANGLE
 6010 FORMAT(1H ,I8,8X,A40)
 6020 FORMAT(1H ,F15.6,1X,A40)
 6030 FORMAT(1H ,G15.6,1X,A40)
      RETURN
      END
      SUBROUTINE CONFAB(NPX,FREQL)
      INCLUDE 'compar.f'
      INCLUDE 'comnp.f'
      DIMENSION X(NP2,3)
      EQUIVALENCE (X(1,1),CFF(1,1))
      DIMENSION SECTOR(28)
      CALL VCLR(SECTOR,1,28)
      SECTOR(1)=NPX
      SECTOR(2)=FREQL
      SECTOR(3)=1.
      WRITE(29,444) (SECTOR(JJ),JJ=1,28)
      WRITE(29,444) (X(JJ,2),JJ=1,NPX)
 444  FORMAT(1H ,6G13.5)
      RETURN
      END
