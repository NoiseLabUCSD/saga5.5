c     234567890123456789012345678901234567890123456789012345678901234567890
c     1         2         3         4         5         6         7
c*********************************************
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
      INTEGER m,i,j,iq,jn,ierr,jj,irflag
      integer  ifrflag
      REAL delf,rdstep
      EQUIVALENCE (dumch2,opt)
      CHARACTER*80 dumch,dumch2
      ierr=0
      write(*,*)'Reading input file'
c---  read the GA-parameters
      call readinputstart
      if (iopt(13).eq.1)then    ! we are in range domain
         iopt(1)=2
      endif

      READ(1,'(a)') dumch2
      WRITE(prtfil,'(a)') 'options for oases ',dumch2
      WRITE(*,'(a)') 'options for oases ',dumch2
      READ(1,*)frq(1),freqmax,nfrq
      ifrflag=nfrq
      nfrq=abs(nfrq)
      if (nfrq.gt.Mfreq) stop' increase mfreq'
      frq(2)=freqmax
      if (ifrflag.gt.0) then
         if (nfrq.ge.2) then
            do i=2,nfrq
               frq(i)=frq(1)+(freqmax-frq(1))*(i-1)/(nfrq-1)
            enddo      
         endif
      else                      ! read in individual frequencies  
         read(1,*)(frq(i),i=1,nfrq)
      endif
      
      WRITE(prtfil,*)' number of frequencies',nfrq, frq(1)             
      READ(1,*) nlay       
      if (nlay.gt.NLA) then
         stop ' oast11: To many layers; increase mlay'
      endif       
      WRITE(prtfil,*) 'number of layers ',nlay
      DO 110 m=1,nlay            
         READ(1,*)v(m,1),v(m,2),v(m,3),v(m,4),v(m,5),v(m,6)
         WRITE(prtfil,'(6e12.3)')
     &        v(m,1),v(m,2),v(m,3),v(m,4),v(m,5),v(m,6)
         v(m,6)=v(m,6)*1000.    ! convert units of density
 110  CONTINUE
      write(*,*)' Read in layers'
      READ(1,*) sd         
c     
c     receiver data
c     
      read(1,'(a80)')dumch
      write(*,*)'Peter',dumch
      READ(dumch,*,err=22,end=22) rd,rdlow,ndep,dtilt
      goto 23
 22   READ(dumch,*) rd,rdlow,ndep
 23   continue
      IF ((ndep.gt.mdep)) THEN
         write(*,*) 'ndep must grater than mdep',ndep,mdep
         ierr=1
      ENDIF
      IF ((ndep.gt.NRD)) THEN
         write(*,*) 'ndep must grater than NRD in oases/compar.f'
         write(*,*)'ndep,NRD',ndep,NRD
         ierr=1
      ENDIF
      IF (rd.eq.rdlow) THEN
         rdlow=rdlow*1.00001   
      ENDIF
      irflag=ndep
      ndep=abs(ndep)
      IF (ndep.gt.1) THEN
         rdstep=(rdlow-rd)/float(ndep-1)
      ELSE
         rdstep=1.
      ENDIF
      IF (irflag.gt.0) THEN
c     equidistant receiver depths
         DO 920 jj=1,ndep

 920        rdep(jj)=(jj-1)*rdstep+rd
         ELSE
c     read in individual receiver depths
            READ(1,*) (rdep(jj),jj=1,ndep)
         ENDIF

         READ(1,*)cmin,cmax        
         READ(1,*)nwave,ic1,ic2       
         write(*,*)'nwave,ic1,ic2',nwave,ic1,ic2
c     read(1,*)
c**** 
c---- for wavenumber domain
         if (iopt(1).eq.1) then
            nx=nwave
         endif
         
c***  

         write(*,*)'calling range block'
c---  read the range-block
         call  read_range_inp(ierr,nwave)

c---  should the  EOF be read ? 
         if (iopt(17).eq.1) then
            write(*,*)' Reading EOFs ...  '
            call eofinit
         endif
         write(*,*)'done'

c***  create text strings for physical parameters
         phystxt(1)='Depth of Layer (m)$'
         phystxt(2)='P sound speed (m/s)$'
         phystxt(3)='S sound speed (m/s)$'
         phystxt(4)='P-attenuation (dB/lambda)$'
         phystxt(5)='S-attenuation (dB/lambda)$'
         phystxt(6)='Density (g/cm3)$'
         phystxt(7)='Thickness (m)$'
         phystxt(8)='Source depth (m)$'
         phystxt(9)='Source range (m)$'
         phystxt(11)='Shape coefficient $'
         phystxt(15)='Receiver depth (m)$'
c***  

         do 8 j=1,mphys
         phystxt2(j)= '                                                '
 6          DO 7 I=40,1,-1
               IF(phystxt(j)(I:I).NE.'$') GO TO 7
               phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c     write(*,*)phystxt2(j)
               GO TO 8
 7          CONTINUE
 8       CONTINUE
c     phystxt2(9)='Source range (m)'

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
            IF ((par2phy(i).ne.11) .and.
     &           (par2lay(i).lt.1 .or. par2lay(i).gt.nlay))THEN
               WRITE(*,*)' *** par2lay not correct, parm',i
               WRITE(*,*)' it is larger than the number of layers'
               ierr=1
            ENDIF
            IF (par2phy(i).eq.11 )THEN
               if (par2lay(i).gt.neof) then
                  write(*,*)' *** Optimzation variable #:',i
                  WRITE(*,*)' *** The shapecoffient number',
     &                 ' is not defined'
                  WRITE(*,*)' par2lay(i)', par2lay(i)
                  ierr=1
               endif
            ENDIF
            IF ((par2phy(i).lt.1 .or. par2phy(i).gt.11).and.
     &           par2phy(i).ne.15) THEN
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
C     OASES V-1.1                       C
C     C
C     Ocean Acoustic and Seismic Exploration      	C
C     Synthesis        		C
C     C
C     (C) Henrik Schmidt			C
C     Massachusetts Institute of Technology		C
C     1988				C
C     C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     
C     MS-FORTRAN - TRANSVERSE ISOTROPY VERSION
C     VERSION 1.1, UPDATE 12-FEB-88.     
C     
      USE global
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnp.f'
      INCLUDE './oases/comnrd.f'
      INCLUDE './oases/comfip.f'
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comoas.h'
c     prognm
      COMPLEX SLOW
      CHARACTER*50 FILENM
c     logical lstop
      real theta1,theta2,dtheta,xrange_ref
      complex cpoint2
      DIMENSION X(NP2,3),FF(2,NP3)           
      DIMENSION FFS(2,NP),XS(NP2),AKM(100),AKM1(100),AKM2(100)     
      LOGICAL NFLAG,ICONTU,lincoh,nofft,ausamp
      EQUIVALENCE (NREC,ISPACE),(LF,NUMFR)
      EQUIVALENCE (FF(1,1),CFF(1,1)),(FFS(1,1),CFFS(1))
      EQUIVALENCE (X(1,1),CFF(1,1)),(XS(1),CFFS(1))        
      integer ibody      
      common /oasesvar/ IPOW,DELTA,THETA,FOCDEP,LTYP,NFLAG,ICNTIN,
     .     ICONTU,lincoh,nofft,ausamp,
     .     ibody
      real              freqm,cminin,cmaxin,cref,ranref,rm,rmi,rr
      common /oasesgard/freqm,cminin,cmaxin,cref,ranref,rm,rmi,rr
C     
C     ********************FORMATS*******************         
C     
 200  FORMAT(20A4)                
 350  FORMAT(//1H ,'    DEPTH        ALPHA       BETA',
     1     '      ATTENA       ATTENB         RHO       ROUGHNESS',
     2     //(1H ,3F12.5,2F12.8,2F12.5))  
 500  FORMAT(//1H ,'SOURCE FREQUENCY:',F10.2,' HZ ')       
 550  FORMAT(//1H ,3X,'NW',5X,'  IC1',5X,'  IC2',/(1X,I5,2I10))
 600  FORMAT(//,'  CMIN = ',G15.6,' M/S ',/,'  CMAX = ',G15.6,' M/S ')
C     
C     
C     **********************************************         
C     

      write(*,*)' Using OAST  as forward model'
      WRITE(prtfil,*)' Using OAST  as forward model'
      DEBUG=.FALSE. 
      DECOMP=.FALSE.
      LUGRN=30
      cminin=cmin
      cmaxin=cmax
C     
C     CHECK THAT NP IS INTEGER POWER OF 2
C     
      IPOW=0
 992  IPOW=IPOW+1
      II=2**IPOW
      IF (II.LT.NP) GO TO 992
      IF (II.NE.NP) THEN
         WRITE(prtfil,*)
         WRITE(prtfil,*) 
     &        '>>>> FATAL: PARAMETER NP MUST BE',
     &        ' INTEGER POWER OF 2 <<<<'
         WRITE(prtfil,*)
     &        '>>>>        CHANGE IN FILE compar.f AND RELINK. ',
     &        '    <<<<'      
      END IF
      NUML=Nlay
      FREQ=frq(1)     
      freqm=frq(nfrq)
      DO 5 M=1,NUML
         ROUGH(m)=0.
 5    CONTINUE   
      
      ir=ndep
      do i=1,ir
         rdc(i)=rdep(i)
      enddo
      nwvno=nwave
      write(*,*)' nwvno=nwave', nwvno,nwave
      icw1=ic1
      icw2=ic2
      ITXX=NP2
      PI=4E0*ATAN(1E0)           
      AI=CMPLX(0.,1.)             
      CNUL=CMPLX(0.,0.)
      LS=1
      DELTA=1.
      THETA=0.
      FOCDEP=0.
      LTYP=1
      LINA=0
      NFLAG=.FALSE.
      ICNTIN=0
      nplots=nx
c     >>> Hanning tapering
      icut1=1
      icut2=nwvno
C     *** FOR TL TAPERING BY CHERMIT 
c     ICUT1=ICW1
c     ICUT2=ICW2
C     
c     
c     we use slowness integration for non-automatic integration
      ibody=1
C     
C     DEFAULT: ISOTROPIC LAYERS
C     
      NTISOL=0
C     
      CALL GETOPT(IPROF,ICNTIN,ICONTU,OPt,prtfil,lincoh,nofft,tilt)

      nsensor=nout
      if (nsensor.gt.1) then
c     write(*,*)' only one inversion whith one physical variable'
c     stop
         write(*,*)' ***** Vector sensor processing with',nsensor,
     &        ' variables'
         if (mod(nfrq,nsensor).ne.0) then
            write(*,*)'nfrq,nsenso=r',nfrq,nsensor
            write(*,*)' number of "frequencies" must be a multiplum of'
            write(*,*)' number of sensors and number of frequencies' 
            stop
         endif
         if (nsensor*ir.gt.nrd) then
            write(*,*)' Maximum Number of receiver depths  (nrd) in',
     &           ' file compar.f is too small '
            stop
         endif
      endif
      OFFDB=0E0
C     
      SDC(1)=SD
      IF (IR.GT.NRD) THEN
         WRITE(prtfil,*) '*** TOO MANY RECEIVERS ***',ir,nrd
         STOP
      END IF
C     LAYER NUMBERS OF SOURCE AND RECEIVER ARE DETERMINED
C     
C     DETERMINATION OF SOURCE LAYERS
C     
      WRITE(prtfil,910) NUML,LS,IR
 910  FORMAT(//1H ,'NUMBER OF LAYERS:      ',I3,           
     1     /1H ,'NUMBER OF SOURCES:     ',I3,           
     2     /1H ,'NUMBER OF RECEIVERS:   ',I3)           
      DO  I=1,LS
         CALL SOURCE(V,NUML,SDC(I),LAYS(I),ZUS(I),ZLS(I))
      enddo

      NUMI=NUML-1                 
C     
C     CHECK THAT NWVNO IS INTEGER POWER OF 2
C     
c     IPOW=0
c     991  IPOW=IPOW+1
c     II=2**IPOW
c     IF (II.LT.NWVNO) GO TO 991
c     IF ((II.NE.NWVNO).and.(iopt(1).ne.1)) THEN
c     NWVNO=II
c     WRITE(prtfil,*)
c     WRITE(prtfil,*)'>> NW MUST BE INTEGER POWER OF 2, CHANGED TO',
c     &             NWVNO,' <<<<'      
c     END IF

      NWVNO=MIN0(NWVNO,NP)
      ICw2=MIN0(NWVNO,ICw2)
      ICw1=MAX0(1,ICw1)
      ICUT2=MIN0(NWVNO,ICUT2)
      ICUT1=MAX0(1,ICUT1)
      IF (CMIN.EQ.0) then
         write(*,*)' cmin',cmin
         STOP '*** CMIN MUST BE NON-ZERO ***,cmin='
      endif
      IF (((CMAX.GT.0).AND.(CMIN.GT.0)).OR.
     1     ((CMAX.LT.0).AND.(CMIN.LT.0))) THEN
         IF (CMAX.LT.0.AND.ICDR.NE.1) THEN
            STOP '*** NEGATIVE SPECTRUM ONLY ALLOWED FOR
     1           PLANE GEOMETRY'
         END IF
         NFLAG=.FALSE.
         WK0 = 2*PI*FREQ / CMAX
         WKMAX = 2*PI*FREQ / CMIN               
      ELSE
         IF (ICDR.NE.1) THEN
            STOP '*** NEGATIVE SPECTRUM ONLY ALLOWED FOR
     1           PLANE GEOMETRY ***'
         END IF
         IF (CMIN.LE.0) STOP '*** CMIN/CMAX CONFLICT ***'
         WKMAX=2*PI*FREQ/CMIN
         WK0=-WKMAX
         CMAX=-CMIN
         NFLAG=.TRUE.
         ICUT2=NWVNO
         ICUT1=NWVNO/2+1
      END IF

      WRITE(prtfil,600)CMIN,CMAX       
      WRITE(prtfil,550)NWVNO,ICw1,ICw2                      
      IF (NFLAG) WRITE(prtfil,*) '*** NEGATIVE SPECTRUM BY',
     1     ' SYMMETRY ***'
C     
C     RANGE DATA
C     
      IF (NWVNO.GT.1) THEN
         DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )            
      ELSE
         DLWVNO=1.0
      END IF
      write(*,*)' wavenumber', DLWVNO,wkmax,wk0,icw1,icw2
      if (icw1.ne.1)
     1     write(*,*)' cmax_icw1', 2*PI*FREQ/(wk0+dlwvno*icw1) 
      if (icw2.ne.nwvno)
     2     write(*,*)' cmin_icw2', 2*PI*FREQ/(wk0+dlwvno*icw2) 
      DLRAN=2E0*PI/(NWVNO*DLWVNO) 
      R1=dlran
      LF=NWVNO                    
      RANMAX=NWVNO*DLRAN
      WRITE(prtfil,360) DLRAN,RANMAX*1E-3
 360  FORMAT(1H ,' ',/1H ,'RANGE STEP:   ',F12.3,' m',
     &     ' MAX FFT RANGE:',F12.3,' km')
C     
C     IF OPTION 'J', SET DEFAULT OFFDB
C     
      IF (ICNTIN.GT.0) THEN
         WRITE(prtfil,*)
         IF (NFLAG) THEN
            STOP '*** CONTOUR OFFSET NOT ALLOWED FOR
     &           NEGATIVE SPECTRUM ***'
         ELSE IF (OFFDB.LT.1E-10) THEN
c     OFFDB=60.0*V(LAYS((LS-1)/2+1),2)/(FREQ*RANMAX)
            OFFDB=60.0*V(LAYS((LS-1)/2+1),2)*(1E0/CMIN-1E0/CMAX)/NWVNO
            WRITE(prtfil,*) 'THE DEFAULT CONTOUR OFFSET IS APPLIED'
         ELSE
            WRITE(prtfil,*) 'THE USER DEFINED CONTOUR ', 
     &           'OFFSET IS APPLIED'
         END IF
         ATTMAX=OFFDB*FREQ*RANMAX/V(LAYS((LS-1)/2+1),2)
c     WRITE(prtfil,361) OFFDB,ATTMAX,ATTMAX*XRIGHT*1E3/RANMAX
 361     FORMAT(1H ,'CONTOUR OFFSET:         ',F12.6,' dB/wavelength',
     &        /1H ,'AT MAX FFT RANGE:       ',F12.2,' dB',
     &        /1H ,'AT MAX WINDOW RANGE:    ',F12.2,' dB',
     &        /1H ,'>> NOTE THAT COMPENSATION IS AUTOMATIC <<')
      END IF

      FNIFAC =  SQRT( 2.0/PI ) * 0.5               
      IF (ICDR.EQ.1) THEN
         FNI5=DLWVNO*SQRT(FREQ/V(LAYS((LS-1)/2+1),2))
      ELSE
         FNI5=DLWVNO*FNIFAC
      END IF
 6010 FORMAT(1H ,I8,10X,A40)
      flagpu=-3
      AUSAMP=(NWVNO.LT.1)
c     write(*,*)' Automatic sampling',ausamp,prtfil,nwvno,(NWVNO.LT.1)
c     
c**************************
      entry forw2( )
      do i=1,ir
         rdc(i)=rdep(i)
      enddo
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

      CALL INENVI



      if (flagpu.lt.0) WRITE(prtfil,918)
 918  FORMAT(//1H ,'RECEIVER DATA:',//1H ,'  NO. ','        DEPTH  ',
     1     'LAYER','           Z')
C     
C     EQUIDISTANT RECEIVER DEPTHS
C     
      DO 921 JJ=1,IR
         CALL RECEIV(V,NUML,RDC(JJ),LAY(JJ),Z(JJ))                 
         if (flagpu.lt.0) WRITE(prtfil,907) JJ,RDC(JJ),LAY(JJ),Z(JJ)
 921  CONTINUE
c---  for a new source reset the depths
      if (flagpu.lt.0)WRITE(prtfil,908)
 908  FORMAT(//1H ,'SOURCE DATA:',//1H ,'  N0. ','        DEPTH  ',
     1     'LAYER','           ZU             ZL')
      DO 906 I=1,LS
         CALL SOURCE(V,NUML,SDC(I),LAYS(I),ZUS(I),ZLS(I))
 906  CONTINUE
      if (flagpu.lt.0) then
         WRITE(prtfil,907) 1,SDC(1),LAYS(1),ZUS(1),ZLS(1)
         IF (LS.GT.1) THEN
            WRITE(prtfil,907) LS,SDC(LS),LAYS(LS),ZUS(LS),ZLS(LS)
         END IF
 907     FORMAT(1H ,I6,G15.6,I6,2G15.6)
      endif
c---  for analytic gradient we uses pointers between input and forward model.
c     w/o gradient it is a dummy routine in oastsub
      call map2forw(1)

      if (tilt) then
         if (nplots.gt.1) stop 'tilt is only allowed for one rage'
         xrange_ref=xranges(1)
      endif


      IF (AUSAMP) THEN
         AUSAMP=.TRUE.
         cref=0e0
         do 981 l=1,numl
            if (v(L,2).lt.2E4) cref=max(cref,v(L,2))
 981     continue
         RANREF=0
c     >>> max and min ranges
         rm =0e0
         rmi=1e20
         do 988 ii=1,nplots
c     rr=abs(1E3*(r0+(ii-1)*rspace))
            rr=abs(xranges(ii))
            if (tilt) then
               rr=rr+dtilt
            endif
            rm=max(rm,rr)
            if (rr.gt.0e0) rmi=min(rmi,rr)
 988     continue
         if (rmi.gt.1e19) rmi=0e0
c     RANREF=RANREF+RM
         ranref=min(ranref+2*rm,6*rm)
         OFFDBIN=0E0
         if (flagpu.lt.0) then
            WRITE(6,*)
            WRITE(6,*) '>>> AUTOMATIC SAMPLING '
            write(6,*) '    REFERENCE SPEED:',CREF
            WRITE(6,*) '    REFERENCE RANGE:',RANREF
         endif
         FREQ=FREQM
         ibody=0
         CALL AUTSAM(CMININ,CMAXIN,rmi,RANREF,CMIN,CMAX,
     &        NWVNO,ICW1,ICW2,ibody)
c     WRITE(6,*) '    MAX NO. OF WAVENUMBERS:',NWVNO
         ICUT1=1
         ICUT2=NWVNO
      ELSE
         CMIN=CMININ
         CMAX=CMAXIN
      END IF                    ! Atomatic sampling loop

      CALL PINIT1

c---------------------Frequency loop
      do 15 jj=1,nfrq,nsensor
         freq=frq(jj)

         DSQ=2E0*PI*FREQ             
         CSQ=DSQ*DSQ                 
         RDSQ=DSQ

         IF (AUSAMP) THEN
C     ***  AUTOMATIC SAMPLING
            
            IF (IBODY.GT.0) THEN
C     ***   SLOWNESS INTEGRATION
               CALL AUTSAM(CMININ,CMAXIN,rmi,RANREF,CMIN,CMAX,
     &              NWVNO,ICW1,ICW2,ibody)
            ELSE
C     ***   WAVENUMBER INTEGRATION
               CALL AUTSAM(CMININ*FREQ/FREQM,CMAXIN,rmi,RANREF,
     &              CMIN,CMAX,NWVNO,ICW1,ICW2,ibody)
            END IF
            if (flagpu.lt.2) then
               WRITE(6,*)
               WRITE(6,*) '>>> FREQUENCY:',FREQ
               WRITE(6,*) '    CMIN,CMAX:  ',CMIN,CMAX
               WRITE(6,*) '    WAVENUMBERS:',NWVNO,ICW1,ICW2
            endif
            ICUT1=1
            ICUT2=NWVNO
            WK0=RDSQ/CMAX           
            WKMAX=RDSQ/CMIN
            IF (NWVNO.GT.1) THEN
               DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )         
            ELSE
               DLWVNO=1.0
            END IF
            IF (ICNTIN.GT.0) THEN
               OFFDB=60.0*V(LAYS((LS-1)/2+1),2)*
     &              (1E0/CMIN-1E0/CMAX)/NWVNO
               if (flagpu.lt.2) then
                  WRITE(6,*) 'DEFAULT CONTOUR OFFSET APPLIED,',OFFDB,
     &                 ' dB/wavelength'
               endif
            END IF
            IF (NFLAG) THEN
               WRITE(6,*) 'NEGATIVE SPECTRUM BY SYMMETRY'
            END IF
         ELSE
            IF (.NOT.NFLAG) THEN
               IF (CMAX.GT.0) THEN
                  WK0=RDSQ/CMAX           
               ELSE
                  WK0=RDSQ/CMIN-(NWVNO-1)*DLWVNO
               END IF
            END IF
            IF (IBODY.GT.0) THEN
               WK0=RDSQ/CMAX           
               WKMAX=RDSQ/CMIN
               IF (NWVNO.GT.1) THEN
                  DLWVNO = ( WKMAX - WK0 ) / ( FLOAT(NWVNO-1) )            
               ELSE
                  DLWVNO=1.0
               END IF
            END IF
         END IF

c     WK0=RDSQ/CMAX           
         DLRAN=2E0*PI/(NWVNO*DLWVNO) 
         R1=dlran
         RANMAX=NWVNO*DLRAN
         
         FNIFAC =  SQRT( 2.0/PI ) * 0.5               
         IF (ICDR.EQ.1) THEN
            FNI5=DLWVNO*SQRT(FREQ/V(LAYS((LS-1)/2+1),2))
         ELSE
            FNI5=DLWVNO*FNIFAC
         END IF 
c     
c***  and now set the discrete sampling points
c     
         if (iopt(1).eq.2) then
            do i=1,nx
               xhelp=xranges(i)/dlran
               ixpoints(i)=xhelp          
               xweight(i)=xhelp-ixpoints(i)
c     write(*,*)' from oast',i,xranges(i),xweight(i),ixpoints(i)
               if (xranges(i).gt.0.66*ranmax) then
                  write(*,'(//,a,f12.3)')' >>> Error For frequency',
     &                 freq
                  WRITE(*,362) DLRAN,RANMAX*1E-3
 362              FORMAT(' RANGE STEP:   ',F12.3,' m',
     &                 '  MAX FFT RANGE:',F12.3,' km')
c     write(*,*)frq(jj),ibody,dlwvno
                  write(*,*)'** ERROR ** range is greater than 0.66',
     &                 ' times the max fft-range, for range',
     &                 xranges(i)  
                  stop
               endif
               if (ixpoints(i).eq.0) then
                  write(*,*)'range to small relative to sampling'
                  write(*,*)'for range',xranges(i)       
                  stop
               endif
            enddo

         endif

         CALL PINIT2                 
         CALL PHASES(LS,FREQ,V,DELTA,THETA,LTYP,FOCDEP)
c     if (flagpu.lt.0)  write(*,*)'frq:',jj,freq
         CALL CALINT
c     if (flagpu.lt.0) write(*,*)'exited calint'
C******************
         DO 120 NREC=1,IR
            do 120 isensor=1,nsensor ! number of vector measurements 
               if (iopt(1).eq.1) then
                  index=(nrec-1+(jj-1)*ir)*(icut2-icut1+1)
               else
                  index=(nrec-1+(jj-1)*ir)*(nx) 
     &                 + (isensor-1)*ir*nx
               endif
               ir0=ir*(isensor-1)

               if (iopt(1).eq.1) then !inversion in wavenumber domain
                  do ii=icut1,icut2
                     resp(ii+index)=wavenoint(ii,nrec)
                  enddo
               elseif (iopt(1).eq.2) then
c     if ((flagpu.lt.-1)) write(*,*)'entering tloss loop'
c---  multiply with sqrt (wavenumber) (moved from calint)
                  do ii=icut1,icut2
                     wavenoint(ii,nrec+ir0)=
     &                    wavenoint(ii,nrec+ir0)*facsqrt(ii)
                  enddo
c     WRITE(prtfil,*)'msoast11 wavenoint:',nrec+ir0, 
c     &           wavenoint(2,nrec+ir0)
                  if (nofft) then
c     if((flagpu.lt.0).and.(nrec.eq.1))
c     &              write(*,*)' Using trapezoidal integration'
                     if (tilt) then
                        xranges(1)=xrange_ref+dtilt/(ir-1)*(nrec-1)
                        if (flagpu.lt.0) 
     &                       write(prtfil,*)'response for range',
     &                       xranges(1)
                     endif
                     call intgrn_pg(xranges,nrec+ir0,nflag)
c     write(*,*)' exited intgrn',nrec
                  else
                     CALL PHINT(NFLAG)
                     CALL TLOSS2
                  endif
                  if (lincoh) then
c     Johns incoherent averaging
                     delr=dlran
                     do 90 jr=1,nx
                        ix=(xranges(jr)-.1)/delr
                        iave=xranges(jr)*.116/delr
                        rpowave=0.
                        j=ix
                        jlow=j-iave
                        if(j-iave.lt.1) jlow=1
                        jup=j+iave
                        if(jup.gt.nwvno) jup=nwvno
                        do 70 k=jlow,jup
                           rpowave=rpowave+abs(cff(k,1))
 70                     continue
c     write(*,*)rpowave,jlow,jup,nwave,jr
                        resp(jr+index)=(rpowave)/(jup-jlow+1)
c     write(*,*)' oast', resp(jr+index)
 90                  continue
                     
                  elseif (nofft) then
c     write(*,*)'writing resp'
                     do i=1,nx
                        resp(i+index)=cff(i  ,1)
                     enddo
                     if (iopt(6).eq.1)
     .                    write(*,*)'resp(ran=1,dep,frq):'
     .                    ,nrec,jj,resp(1+index),index
                  else
                     dtheta= dlran*2*pi*freq/1500
                     do i=1,nx  ! ranges
                        theta1=
     &                       ATAN2(AIMAG(cff(ixpoints(i)  ,1)),
     &                       REAL(cff(ixpoints(i)  ,1)))
                        cpoint2=cff(ixpoints(i)+1,1)
                        theta2=
     &                       ATAN2(AIMAG(cpoint2),REAL(cpoint2))+dtheta
                        if (theta2.gt.2*pi) theta2=theta2-2*pi
                        if (abs(theta1-theta2).gt.pi) then
                           if ((theta1.lt.0).and.(theta2.gt.0))then
                              theta1=theta1+2.*pi
                           elseif ((theta2.lt.0).and.(theta1.gt.0))then
                              theta2=theta2+2.*pi
                           endif
                        endif
                        theta=theta1*(1.-xweight(i))
     &                       +(theta2-dtheta)*(xweight(i))
                        resp(i+index)=
     &                       (abs(cff(ixpoints(i)  ,1))*(1.-xweight(i))
     &                       +abs(cff(ixpoints(i)+1,1))*xweight(i))
     &                       *exp(cmplx(0.,theta))
                     enddo
                  endif
               endif
 120        continue            ! receiver depth
            if (tilt) then
               xranges(1)=xrange_ref
            endif
 15      continue               ! frequency
         lwascomp=1

c     c: Write transfer function in a file:
         if ((iWriteTrf.gt.0)) then
            write(6,*)' Calling  WRITEoasTRF for best model'
            write(6,*)' This is not yet working'
            CALL WRITEoasTRF(LUN,FILENAME,TITLE,RD,RDLOW,R0,RSPACE,
     &           NX,LX,MX,DT,FREQS,SD)


         endif
         
         return
         END  
C     
      SUBROUTINE GETOPT(IPROF,ICNTIN,ICONTU,opt,prtfil,lincoh,nofft,
     &     tilt)
C     
C     INPUT OF OPTIONS            
C     
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comfip.f'
      CHARACTER OPT(40)
      LOGICAL ICONTU,lincoh,nofft,tilt
      integer prtfil
      write(*,*)'from getopt',opt
      WRITE(prtfil,300)                
 300  FORMAT(//,' OPTIONS:',/)    
      DEBUG=.FALSE.
      DEPTAV=.FALSE.
      DRCONT=.FALSE.
      PLTL=.FALSE.
      PLKERN=.FALSE.
      ANSPEC=.FALSE.
      TLDEP=.FALSE.
      SHEAR=.FALSE.
      ICONTU=.FALSE.
      SCTOUT=.FALSE.
      lincoh=.FALSE.
      nofft=.true.
      inttyp=0
      NOUT=0                      
      IREF=0                      
      ISTYP=0                     
      ICDR=0
      IPROF=0
      ICNTIN=0
      DO 10 I=1,3                 
 10      IOUT(I)=0                   
         DO 50 I=1,40                
            IF (OPT(I).EQ.'N') THEN 
               IF (IOUT(1).GT.0) GO TO 50  
               NOUT=NOUT+1                 
               IOUT(1)=1                   
               WRITE(prtfil,301)                
 301           FORMAT(1H ,'NORMAL STRESS') 
            ELSE IF (OPT(I).EQ.'V') THEN 
               IF (IOUT(2).GT.0) GO TO 50  
               NOUT=NOUT+1                 
               IOUT(2)=1                   
               WRITE(prtfil,302)                
 302           FORMAT(1H ,'VERTICAL VELOCITY')                   
            ELSE IF (OPT(I).EQ.'H') THEN 
               IF (IOUT(3).GT.0) GO TO 50  
               NOUT=NOUT+1                 
               IOUT(3)=1                   
               WRITE(prtfil,303)                
 303           FORMAT(1H ,'HORIZONTAL VELOCITY')                 
            ELSE IF (OPT(I).EQ.'A') THEN 
               IF (DEPTAV) GO TO 50     
               DEPTAV=.TRUE.
               WRITE(prtfil,304)                
 304           FORMAT(1H ,'DEPTH AVERAGE')
            ELSE IF (OPT(I).EQ.'C') THEN 
               IF (DRCONT) GO TO 50    
               DRCONT=.TRUE.
               WRITE(prtfil,305)           
 305           FORMAT(1H ,'DEPTH-RANGE CONTOURS')
            ELSE IF (OPT(I).EQ.'c') THEN 
               IF (ICONTU) GO TO 50    
               ICONTU=.TRUE.
               WRITE(prtfil,306)           
 306           FORMAT(1H ,'DEPTH INTEGRAND CONTOURS')
            ELSE IF (OPT(I).EQ.'T') THEN 
               IF (PLTL) GO TO 50    
               PLTL=.TRUE.
               WRITE(prtfil,307)           
 307           FORMAT(1H ,'TRANSMISSION LOSS')
            ELSE IF (OPT(I).EQ.'I') THEN
               IF (PLKERN) GO TO 50
               PLKERN=.TRUE.
               WRITE(prtfil,309)
 309           FORMAT(1H ,'HANKEL TRANSFORM INTEGRANDS')
            ELSE IF (OPT(I).EQ.'P') THEN
               IF (ICDR.GT.0) GO TO 50
               ICDR=1
               WRITE(prtfil,313)
 313           FORMAT(1H ,'PLANE GEOMETRY')
            ELSE IF (OPT(I).EQ.'L') THEN
               IF (LINA.GT.0) GO TO 50
               LINA=1
               WRITE(prtfil,308)
 308           FORMAT(1H ,'LINEAR VERTICAL SOURCE ARRAY')
            ELSE IF (OPT(I).EQ.'S') THEN
               IF (ANSPEC) GO TO 50
               ANSPEC=.TRUE.
               WRITE(prtfil,3091)
 3091          FORMAT(1H ,'ANGULAR SPECTRA')
            ELSE IF (OPT(I).EQ.'s') THEN
               IF (SCTOUT) GO TO 50
               SCTOUT=.TRUE.
               WRITE(prtfil,3092)
 3092          FORMAT(1H ,'OUTPUT OF SCATTERING DISCONTINUITIES')
            ELSE IF (OPT(I).EQ.'Z') THEN
               IF (IPROF.GT.0) GO TO 50
               IPROF=1
               WRITE(prtfil,314)
 314           FORMAT(1H ,'PLOT OF VELOCITY PROFILES')
            ELSE IF (OPT(I).EQ.'J') THEN
               IF (ICNTIN.GT.0) GO TO 50
               ICNTIN=1
               WRITE(prtfil,315)
 315           FORMAT(1H ,'COMPLEX INTEGRATION CONTOUR')
            ELSE IF (OPT(I).EQ.'D') THEN
               IF (TLDEP) GO TO 50
               TLDEP=.TRUE.
               WRITE(prtfil,316)
 316           FORMAT(1H ,'TRANSMISSION LOSS AS FUNCTION OF DEPTH')
            ELSE IF (OPT(I).EQ.'X') THEN
               IF (SHEAR) GO TO 50
               SHEAR=.TRUE.
               WRITE(prtfil,317)
 317           FORMAT(1H ,'SHEAR SOURCES IN SOLID MEDIA')
            ELSE IF (OPT(I).EQ.'Q') THEN
               IF (DEBUG) GO TO 50
               DEBUG=.TRUE.
               WRITE(prtfil,318)
 318           FORMAT(1H ,'***** DEBUGGING *****')
            ELSE IF (OPT(I).EQ.'i') THEN
               IF (lincoh) GO TO 50
               lincoh=.TRUE.
               nofft=.false.
               WRITE(prtfil,319)
 319           FORMAT(1H ,'JOHNs INCOHERENT AVERAGING ! ')
            ELSE IF (OPT(I).EQ.'p') THEN
               IF (nofft) GO TO 50
               nofft=.TRUE.
               WRITE(prtfil,320)
 320           FORMAT(1H ,'Trapezodal integration ! ')
            ELSE if (opt(i).eq.'t') then
               tilt=.true.
               write(*,*)' tilt of array is included'
               write(prtfil,*)' tilt of array is included'
            ELSE IF (opt(I).EQ.'!') THEN
               goto 60
            ELSE IF (OPT(I).NE.' ') THEN
               WRITE(prtfil,399) OPT(I)
 399           FORMAT(1H ,'>>>> UNKNOWN OPTION: ',A1,' <<<<')
            ELSE
            END IF
 50      CONTINUE                    
 60      CONTINUE                    
         IF (NOUT.NE.0) RETURN       
         IOUT(1)=1                   
         NOUT=1                      
         WRITE(prtfil,301)                
         RETURN                      
         END  



      SUBROUTINE AUTSAM(C1,C2,rmin,RMAX,CMIN,CMAX,NW,IC1,IC2,ibody)
c     PARAMETER (RFAC=1.5)
      PARAMETER (NR=150)
c     >>> Minimum # of k, # of periods of exponential over taper interval
      PARAMETER (NKMIN=200,WTAPER=10.0)
      INCLUDE './oases/compar.f'
      INCLUDE './oases/comnla.f'
      INCLUDE './oases/comnrd.f'
      if (inttyp.eq.1) then
c     filon:
         RFAC=1.e0
      else
         RFAC=1.5e0
      end if
C     *** DETERMINE WAVENUMBER SAMPLING INTERVAL
      DK=2*PI/(RFAC*RMAX)
c     >>> Determine tapering interval from 
c     minimum source receiver separations
c     >>> Depth-separation
      zsep=1E20
      do 200 isou=1,ls
         do 200 ircv=1,ir
            zsep=min(zsep,abs(rdc(ircv)-sdc(isou)))
c     write(*,*)rdc(ircv),sdc(isou)
 200     continue
         rref=rmin+wtaper*zsep
c     >>> if zero set to  one  wavelength
         if (rref.eq.0e0) rref= v(lays(1),2)/freq
         taperk=2*pi*wtaper/rref
C     *** INTEGRATION LIMITS
         WN1=2*PI*FREQ/C2
         if (ibody.eq.2) then
            WN2=2*PI*(FREQ+0.2*freq2)/C1
         else
            WN2=2*PI*FREQ/C1
         end if
c     if (ibody.gt.0) then
c     WNMAX=(1.1+(freq2-freq)/freq2*0.9)*(WN2-WN1)+WN1
c     else
c     WNMAX=1.1*(WN2-WN1)+WN1
c     end if
c     WNMIN=WN1-0.1*(WN2-WN1)
         wnmax=wn2+taperk
         wnmin=wn1-taperk
         WNMIN=MAX(WNMIN,DK*1e-3)
         DK=MIN(DK,(WNmax-WNmin)/(NKMIN-1))
C     *** NUMBER OF WAVENUMBER SAMPLING POINTS
         NW=(WNMAX-WNMIN)/DK+1
         WNMAX=WNMIN+(NW-1)*DK
         IC1=(WN1-WNMIN)/(WNMAX-WNMIN)*(NW-1)+1
         IC1=MAX(IC1,1)
         IC2=(WN2-WNMIN)/(WNMAX-WNMIN)*(NW-1)+1
         IC2=MIN(IC2,NWVNO)
         CMIN=2*PI*FREQ/WNMAX
         CMAX=2*PI*FREQ/WNMIN
         if (flagpu.lt.0) then
            write(6,*) '>>> Reference depth for tapering:', zsep,' m'
            write(6,*) '>>> Reference range for sampling:',rmax,' m'
            write(6,*) '>>> Reference range for tapering:', rref,' m'
         endif
         if (nw.gt.np) then
            write(*,*)' Maximum wavenumber dimension ',NP
            write(*,*)' Requested wavenumber dimension ',NW
            Write(*,*)' Stop increase dimension'
            stop
         endif
         RETURN
         END





      BLOCK DATA SAFBK1
      INCLUDE './oases/compar.f'      
C**** DEFINITION OF MAX REAL ARGUMENT TO THE EXPONENTIAL FUNCTION
      COMMON /ARGMAX/ AM
C**** THE FOLLOWING DEFINITION SHOULD BE USED FOR THE FPS164
C     FPS  DATA AM /300./
C**** THE FOLLOWING DEFINITION SHOULD BE USED FOR THE VAX
      DATA AM /65./     
      DATA PROGNM /'OASTL '/
      DATA OMEGIM /0.0/
      DATA LUGRN,LUTRF,LUTGRN,LUTTRF /30,35,30,35/
      DATA SHEAR,DECOMP,SCTOUT /.FALSE.,.FALSE.,.FALSE./
      END














