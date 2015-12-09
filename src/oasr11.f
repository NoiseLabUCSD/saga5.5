c     234567890123456789012345678901234567890123456789012345678901234567890
c     1         2         3         4         5         6         7
      subroutine forwardmodel(iopt,mopt)
      integer   mopt,i,iopt(mopt)
      DO i=1,mopt
         iopt(i)=0.
      ENDDO
      iopt(1)=3
      iopt(30)=4
      end
c*****************************
      SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization program
c     PETER GERSTOFT, 1992
c     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INTEGER m,i,j,iq,jn,ierr,jj
      REAL delf,rdstep
      CHARACTER*40 optread
      EQUIVALENCE (dumch2,opt)
      CHARACTER*80 dumch,dumch2
      ierr=0
c---  read the GA-parameters
      call readinputstart
        write(*,*)' any0'

      READ(1,'(a)') dumch2
      WRITE(prtfil,'(a)') 'options for oases ',dumch2
      WRITE(*,'(a)') 'options for oases ',dumch2
       READ(1,*) nlay       
      if (nlay.gt.NLA) then
         stop 'oast11: To many layers; increase mlay'
      endif       
        write(*,*)' any2'
c     WRITE(prtfil,*) 'nb of layers ',nlay
      DO 110 m=1,nlay            
         READ(1,*)v(m,1),v(m,2),v(m,3),v(m,4),v(m,5),v(m,6)
         WRITE(prtfil,'(6e12.3)')
     &        v(m,1),v(m,2),v(m,3),v(m,4),v(m,5),v(m,6)
         v(m,6)=v(m,6)*1000.    ! convert units of density
 110  CONTINUE
      write(*,*)' Read in layers'
c     
c     receiver data
c     

      READ(1,*) FRQ(1),FRQ(2),NFRQ
      READ(1,*) ANGLE1,ANGLE2,NANG
      WRITE(prtfil,*)'nb of frequencies',nfrq, frq(1)             
c**** 
      freqmax = FRQ(2)
      if (nfrq.ge.2) then
         do i=2,nfrq
            frq(i)=frq(1)+(freqmax-frq(1))*(i-1)/(nfrq-1)
         enddo      
      endif

      nx=nang
      iopt(1)=3
c     
c---  number of curves to optimize on
      ncurv=nfrq
      ndep=1
c     
c     

      if (iopt(13).gt.0) then   ! for covariance matrix
         stop 'covariance not allowed for oasr'
      endif
      if (nang*nfrq.gt.mobs) then
         write(*,*) ' nx*ndep >mobs, nx,ndep,mobs',nx,ndep,mobs
         ierr=1
      ENDIF
      
c***  
      if (iopt(17).eq.1) then
         call eofinit
      endif

c***  create text strings for physical parameters
      phystxt(1) ='Depth $'
      phystxt(2) ='P-velocity $'
      phystxt(3) ='S-velocity $'
      phystxt(4) ='P-attenuation $'
      phystxt(5) ='S-attenuation $'
      phystxt(6) ='Density $'
      phystxt(7) ='Thickness $'
      phystxt(8) ='Source depth $'
      phystxt(9) ='Source range $'
      phystxt(11)='Shape coefficient $'
c***  

      do 8 j=1,mphys
         phystxt2(j)='                                               '
 6       DO 7 I=40,1,-1
            IF(phystxt(j)(I:I).NE.'$') GO TO 7
            phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c     write(*,*)phystxt2(j)
            GO TO 8
 7       CONTINUE
    8 CONTINUE

c---- read the optimization parameter
      call readoptparm

      DO i=1,nparm
c     write(*,*)'parm no ',i
         IF (par2phy(i).eq.6)THEN
            fmin(i)=fmin(i)*1e3
            fmax(i)=fmax(i)*1e3
            df(i)= df(i)*1.e3
c     DO j=0,ndigit(i)-1
c     fval(j,i)=fval(j,i)*1e3
c     ENDDO
         ENDIF
c     write(*,*)'parm no ',i
         
         IF ((fmin(i).gt.fmax(i)) .and. (par2phy(i).ne.3))THEN
            WRITE(*,*)' *** fmin > fmax for parm',i
            ierr=1
         ENDIF
c     write(*,*)'parm no ',i
         IF (par2phy(i).eq.11 )THEN
            if (par2lay(i).gt.neof) then
               write(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** The shapecoffient number is not defined'
               WRITE(*,*)' par2lay(i)', par2lay(i)
               ierr=1
            endif
         ENDIF
         IF ((par2phy(i).ne.11) .and.
     &        (par2lay(i).lt.1 .or. par2lay(i).gt.nlay))THEN
            WRITE(*,*)' *** par2lay not correct, parm',i
            WRITE(*,*)' it is larger than the number of layers'
            ierr=1
         ENDIF
         IF (par2phy(i).lt.1 .or. par2phy(i).gt.11)THEN
            WRITE(*,*)' *** par2phy not correct, parm',i
            ierr=1
         ENDIF
      ENDDO

c***  errors ?
      IF (ierr.eq.1)stop 'stopped in input'
c***  close input unit
      CLOSE(1)                                
      END

c****************************************
      subroutine forwinit
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     C
C     OASES V 1.1                       C
C     C
C     Ocean Acoustic and Seismic Exploration    	C
C     Synthesis         		C
C     C
C     (C)  Henrik Schmidt                   C
C     Department of Ocean Engineering		C
C     MIT				C
C     C
C     1988				C
C     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     REFLECTION COEFFICIENT VERSION
C     
      USE global
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
      INCLUDE './oases/comnrd.f'
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
      COMPLEX SLOW
c     CHARACTER*50 FILENM
c     LOGICAL NFLAG
      DIMENSION X(NP2,3),FF(2,NP3)  
      DIMENSION FFS(2,NP),XS(NP2)
c     CHARACTER*6  OPTION(2)
      LOGICAL PHPLOT,CONTUR
      EQUIVALENCE (X(1,1),CFF(1,1)),(XS(1),CFFS(1))              
      EQUIVALENCE (FF(1,1),CFF(1,1)),(FFS(1,1),CFFS(1))              
c     COMPLEX CCSQ,CDSQ
      integer iparm, LTYP, LINA
      real DELTA,theta,FOCDEP, OFFDB
      common /iparmsaved/ DELTA,theta,FOCDEP, 
     1     iparm, LTYP,PHPLOT,CONTUR
C     
C     ********************FORMATS*******************               
C     
C     
 350  FORMAT(//1H ,'     DEPTH        ALPHA       BETA      ATTENA   ',
     1     '    AT        TENB         RHO     ',
     1     ' ROUGHNESS'//(1H ,3F12.5,2F12.8,2F12.5))              
C     
C     **********************************************               
      DEBUG=.FALSE.
c     DEBUG=.TRUE.
      write(*,*)' Using OASR  as forward model'
      WRITE(prtfil,*)' Using OASR  as forward model'
      NUML=Nlay
      FREQ=frq(1)     
      do 5 m=1,nlay
         ROUGH(m)=0.
 5    CONTINUE   
      
      FREQ1=frq(1)
      FREQ2=frq(nfrq)
      NFREQ=nfrq
      WRITE(prtfil,*)'nb of frequencies',nfreq, freq1,freq2             
C     
C     
      ITXX=NP2
      PI=4E0*ATAN(1E0)      
      AI=CMPLX(0.,1.)        
      CNUL=CMPLX(0.,0.)
      IR=1  
      LS=1
      DELTA=1.
      THETA=0.
      FOCDEP=0.
      LTYP=1
      LINA=0
      OFFDB=0E0
C     
      CALL GETOPT(IPARM,PHPLOT,CONTUR,IPROF,opt,prtfil)
      if ((iout(1)+iout(2)+iout(3)).gt.1) then
         write(*,*)' only one inversion whith one physical variable'
         stop
      endif
C     
      SD=V(2,1)-1.0
      SDC(1)=SD
      rdc(1)=sd
C     
C     IF LOSSLESS WATER, DISABLE SCRETTING-LEROY ATTENUATION
C     
      IF (ABS(V(1,3)).LT.1E-10.AND.V(1,4).LT.1E-10) THEN
         V(1,4)=1E-8
      END IF
C     
C     DETERMINATION OF SOURCE LAYERS
C     
      WRITE(prtfil,908)
 908  FORMAT(//1H ,'SOURCE DATA:',//1H ,'  N0. ','   DEPTH  ',
     1     'LAYER','      ZU        ZL')
      DO 906 I=1,LS
 906     CALL SOURCE(V,NUML,SDC(I),LAYS(I),ZUS(I),ZLS(I))
         WRITE(prtfil,907) 1,SDC(1),LAYS(1),ZUS(1),ZLS(1)
 907     FORMAT(1H ,I6,F10.1,I6,2F10.1)
         CALL RECEIV(V,NUML,RDC(1),LAY(1),Z(1))
         NUMI=NUML-1            
C     
         IF (FREQ1*FREQ2.EQ.0.0) THEN
            STOP '*** FREQUENCIES MUST BE NON-ZERO, ABORTING ***'
         ELSE IF (CONTUR) THEN
            IF (NFREQ.LE.1) THEN
               STOP '*** CONTOURS REQUIRE NRFR>1 YOU STUPID FOOL ***'
            END IF
            F1LOG=ALOG(FREQ1)
            F2LOG=ALOG(FREQ2)
            DFLOG=(F2LOG-F1LOG)/(NFREQ-1)
         ELSE
         END IF
         WRITE(prtfil,1401) FREQ1,FREQ2,NFREQ
 1401    FORMAT(//1H ,'REFLECTION COEFFICIENTS',
     1        /1H ,'***********************',
     2        //1H ,'MINIMUM FREQUENCY:     ',G14.6,' Hz',
     3        /1H ,'MAXIMUM FREQUENCY:     ',G14.6,' Hz',
     4        /1H ,'NUMBER OF FREQUENCIES: ',I14)
         WRITE(prtfil,1402) ANGLE1,ANGLE2,NANG
 1402    FORMAT(//1H ,'MINIMUM ANGLE:         ',G14.6,' deg.',
     1        /1H ,'MAXIMUM ANGLE:         ',G14.6,' deg.',
     2        /1H ,'NUMBER OF ANGLES:      ',I14)
C     
         NWVNO=NANG
         DLFREQ=0.
         IF (NFREQ.GT.1) THEN
            DLFREQ=(FREQ2-FREQ1)/(NFREQ-1)
         END IF
         ICUT1=1
         ICUT2=NWVNO
         NWVNO=MIN0(NWVNO,NP)
         IF (NANG.GT.1) THEN
            DLANGLE=(ANGLE2-ANGLE1)/(NANG-1)
         ELSE
            DLANGLE=1
         END IF
 6010    FORMAT(1H ,I8,10X,A40)
c     
C     FREQUENCY LOOP         
C     --------------         
         flagpu=-2
         if (iopt(6).eq.1) then
            flagpu=-10000
         endif
c     
c**************************
c     
         entry forw2( )

         flagpu=flagpu+1
         lwascomp=1
c     
c---- Use of EOF
c     
         if (iopt(17).eq.1) then  
            call eofval
            if (lwascomp.eq.-1) then
               return
            endif
         endif

         CALL INENVI(prtfil)

         CALL PINIT1
         DO 15 JJ=1,NFREQ
            IF (CONTUR) THEN
               FREQ=EXP(F1LOG+(JJ-1)*DFLOG)
            ELSE
               FREQ=FREQ1+(JJ-1)*DLFREQ
            END IF
c     write(*,*)' freq',freq,FREQ1,DLFREQ,jj
            DSQ=2E0*PI*FREQ        
            CSQ=DSQ*DSQ            
            CALL PINIT2            
            CALL PHASES(LS,FREQ,V,DELTA,THETA,LTYP,FOCDEP)
            CALL REFLEC(ANGLE1,ANGLE2,IPARM)
            
            do ii=1,nang
c     write(*,*)'x(ii,2)',ii,x(ii,2)
               resp(ii+(jj-1)*nang)=x(ii,2) ! transferring magnitude
            enddo
 15      CONTINUE               

         lwascomp=1
         return
         END   

c     
C     ************************************
c     
      SUBROUTINE GETOPT(IPARM,PHPLOT,CONTUR,IPROF,opt,prtfil)
C     
C     INPUT OF OPTIONS       
C     
      integer prtfil      
      INCLUDE './oases/compar.f'
      LOGICAL PHPLOT,CONTUR
      CHARACTER*1 OPT(40)
      WRITE(prtfil,300)           
 300  FORMAT(//1H ,'OPTIONS:',/)                
      NOUT=0
      IREF=0
      ISTYP=0                
      ICDR=0
      ISTACK=0
      IBODY=0
      ISDEP=0
      IGRP=0
      PHPLOT=.FALSE.
      CONTUR=.FALSE.
      SHEAR=.FALSE.
      IPROF=0
      DO 10 I=1,3            
 10      IOUT(I)=0              
c     READ(1,200) OPT        
c     200  FORMAT(40A1)           
         DO 50 I=1,40           
            IF (OPT(I).EQ.'N') THEN             
               IF (IOUT(1).GT.0) GO TO 50              
               NOUT=NOUT+1            
               IOUT(1)=1              
               IPARM=1
               WRITE(prtfil,301)           
 301           FORMAT(1H ,'P-P REFLECTION COEFFICIENT')             
               GO TO 50               
            ELSE IF (OPT(I).EQ.'S') THEN             
               IF (IOUT(2).GT.0) GO TO 50              
               NOUT=NOUT+1            
               IOUT(2)=1              
               IPARM=2
               WRITE(prtfil,302)           
 302           FORMAT(1H ,'P-SV REFLECTION COEFFICIENT')        
               GO TO 50               
            ELSE IF (OPT(I).EQ.'P') THEN
               IF (PHPLOT) GO TO 50
               PHPLOT=.TRUE.
               WRITE(prtfil,304)
 304           FORMAT(1H ,'PHASE PLOTS')
               GO TO 50
            ELSE IF (OPT(I).EQ.'C') THEN
               IF (CONTUR) GO TO 50
               CONTUR=.TRUE.
               WRITE(prtfil,305)
 305           FORMAT(1H ,'CONTURS PLOTTED IN ANGLE/FREQUENCY')
               GO TO 50
            ELSE IF (OPT(I).EQ.'Z') THEN
               IF (IPROF.GT.0) GO TO 50
               IPROF=1
               WRITE(prtfil,314)
 314           FORMAT(1H ,'PLOT OF VELOCITY PROFILES')
               GO TO 50
            ELSE IF (OPT(I).EQ.'s') THEN
               IF (SCTOUT) GO TO 50
               SCTOUT=.TRUE.
               WRITE(prtfil,315)
 315           FORMAT(1H ,'OUTPUT OF SCATTERING DISCONTINUITIES')
               GO TO 50
            ELSEif (opt(i).ne.' ') then
               WRITE(prtfil,320) OPT(I)
 320           FORMAT(1H ,'>>> UNKNOWN OPTION: ',A,' <<<')
            END IF
 50      CONTINUE               
         IF (NOUT.EQ.1) RETURN  
         IF (NOUT.GT.1) STOP '*** ONLY ONE FIELD PARAMETER ALLOWED***'
         IOUT(1)=1              
         NOUT=1
         IPARM=1
         WRITE(prtfil,301)           
         RETURN
         END

      BLOCK DATA OASRBLK
      INCLUDE './oases/compar.f'
C
      DATA OMEGIM /0.0/
      DATA PROGNM /'OASRC '/
      DATA LUGRN,LUTRF,LUTGRN,LUTTRF /30,35,30,35/
      DATA SHEAR,DECOMP,SCTOUT /.FALSE.,.FALSE.,.FALSE./
      END
