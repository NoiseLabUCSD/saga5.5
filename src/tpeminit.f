c234567890123456789012345678901234567890123456789012345678901234567890
c     1         2         3         4         5         6         7
      SUBROUTINE forwardmodel(iopt,mopt)
      INTEGER  mopt,i,iopt(mopt)
      DO i=1,mopt
         iopt(i)=0.
      ENDDO
      iopt(12)=1                ! using three indexes for adressing variable
      iopt(30)=6                ! 6 is for tpem.
      iopt(1)=2
      END
c******************************
      SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization PROGRAM
c     PETER GERSTOFT, 1992
c     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'    
      REAL wl, fko,delz,zmax,con
      INTEGER n,ln, n34 
      COMMON / pevar / wl, fko, delz, n, ln, zmax, n34, con
c     REAL*8  dr, drout,dzout, dr2
c     COMMON / irhstps / izout(mxzout)
c     COMMON / rhstps / dr, drout, dzout, dr2, zout(mxzout)
      INTEGER j,i,ierr,ifrflag
      REAL freqmax
      CHARACTER*80 dumch2, dumch
      EQUIVALENCE (dumch2,opt)
c     CHARACTER*80 dumtitle
      ierr=0
      nfrq=1
c---  READ the GA-parameters
      CALL readinputstart
      READ(1,'(a)') dumch2       ! This line has no effect in snap
      CALL tpemoption(opt)
      WRITE(prtfil,'(a)') 'options for TPEM ',dumch2
      WRITE(*,*)' reading frequencies...'
      READ(1,*,err=201)frq(1),freqmax,nfrq
      GOTO 202
      
 201  WRITE(*,*)' PLEASE USE NEW FORMAT FOR READING FREQUENCIES'
      WRITE(*,*)' first frequency read',frq(1)
      nfrq=1
      freqmax=frq(1)
      PAUSE
 202  CONTINUE
      ifrflag=nfrq
      nfrq=ABS(nfrq)
      IF (itpem_opt(7).EQ.1 .AND. nfrq.GT.1) THEN
         STOP 'multiple beam inversion works only with one freq'
      ENDIF   
      IF (nfrq.GT.Mfreq) STOP' increase mfreq'
      frq(2)=freqmax
      IF (ifrflag.GT.0) THEN
         IF (nfrq.GE.2) THEN
            DO i=2,nfrq
               frq(i)=frq(1)+(freqmax-frq(1))*(i-1)/(nfrq-1)
            ENDDO      
         ENDIF
      ELSE                      ! read in individual frequencies  
         READ(1,*)(frq(i),i=1,nfrq)
      ENDIF
      sv%freq=frq(1)
      
c     READ(1,*)fmin, fmax, nfreq
c     READ(1,*)sv.freq          !Frequency 
      WRITE(*,*)'Frequency ',sv%freq !Frequency 
      dr=0. 
c     IF (sv.freq.LT.0) THEN
      IF (itpem_opt(1).EQ.1) THEN
c     sv.freq =- sv.freq  
         WRITE(*,*)' Reading dr, delz, ln (the fft-power)'
         READ(1,*)  dr, delz, ln
         WRITE(*,*)' dr, delz, ln =', dr, delz, ln
      ENDIF
c     frq(1) = sv.freq
      IF (itpem_opt(2).EQ.1) THEN
         WRITE(*,*)' Reading Teds coefficients....'
         READ(1,*) ncap
         READ(1,*) cap_coef(1),( cap_coef(i+1),i=1,ncap)
         WRITE(*,*)' ....done'
      ENDIF

      READ(1,*)fieldtype        ! 1:complex field, 0:magnitude    
      IF ( fieldtype.EQ.0) WRITE(*,*) ' PROPAGATION LOSS'                  
      IF ( fieldtype.EQ.1) WRITE(*,*) ' COMPLEX FIELD'                  

      READ(1,*)sv%antht         !antenna height 
      WRITE(*,*)'antenna height ',sv%antht
      READ(1,*)sv%ipat          !antenna pattern
      READ(1,*)sv%bwidth        !beamwidth. This is ignored for Omni ant.
      WRITE(*,*)'antenna bwidth ',sv%bwidth
      WRITE(*,*)'antenna pattern ',sv%ipat
      IF (itpem_opt(7).EQ.1) THEN
         READ(1,*)nbeam 
         IF (nbeam.GT.mbeam) THEN
            WRITE(*,*)'nbeam,mbeam=',nbeam,mbeam
            STOP 'increase mbeam'
         ENDIF
         READ(1,*)(beamelev(i),i=1,nbeam)
         nfrq=nbeam 
         sv%elev= beamelev(1)
      ELSE
         READ(1,*)sv%elev       !This value is ignored for Omni antenna
      ENDIF
      WRITE(*,*)'antenna elevation ',sv%elev
c---- Get min and max height
      READ(1,*,err=101)vnp%hmin, vnp%hmax !Coverage up to 1000. m.
      GOTO 102      
 101  WRITE(*,*)' >>>>>>>>>>> Warning using old format, '
      WRITE(*,*)' please specify both hmin and hmax'
      vnp%hmin=0
      BACKSPACE(1)
      READ(1,*)vnp%hmax   
 102  WRITE(*,*)'min & max height ',vnp%hmin, vnp%hmax

c---- Get min and max range
      READ(1,*,err=103)vnp%rmin, vnp%rmax !Range up to 100 km.      
      GOTO 104      
 103  WRITE(*,*)' >>>>>>>>>>> Warning using old format, '
      WRITE(*,*)' please specify both rmin and rmax'
      vnp%rmin=0
      BACKSPACE(1)
      READ(1,*)vnp%rmax   
 104  CONTINUE


      READ(1,*)vnp%nrout,vnp%nzout !number of  height points to output.
      WRITE(*,*)'Output height and range points ',vnp%nzout,vnp%nrout
      ndep=vnp%nzout            !number of  height points to output.

      IF (vnp%rmin.LE.0) THEN
         vnp%dr=(vnp%rmax-vnp%rmin)/vnp%nrout
         vnp%rmin= vnp%dr
      ENDIF    

      IF (vnp%nrout.NE.1) THEN
         vnp%dr = (vnp%rmax-vnp%rmin)/(vnp%nrout-1)
      ELSE
         vnp%dr =  vnp%rmax
      ENDIF
c     pg 8 aug 00 this does not make sense
c     pg 8 aug     vnp.rmin =  vnp.rmin/(INT(vnp.dr/vnp.rmin)+1)
      vnp%rmin =  vnp%dr*(INT(vnp%rmin/vnp%dr+0.99))
      vnp%rmin = MAX(vnp%dr,vnp%rmin)
      vnp%rmax = vnp%rmin+vnp%dr*(vnp%nrout-1)
      
      WRITE(*,*)'minimum, maximum rec range and range step',
     1     vnp%rmin,vnp%rmax,vnp%dr

      IF (vnp%hmin.LE.0) THEN
         vnp%dz  =(vnp%hmax-vnp%hmin)/vnp%nzout
         vnp%hmin= vnp%dz
      ELSE    
         IF (vnp%nzout.NE.1) THEN
            vnp%dz = (vnp%hmax-vnp%hmin)/(vnp%nzout-1)
         ELSE
            vnp%dz = vnp%hmax
         ENDIF
         vnp%hmin = MAX(vnp%dz, vnp%hmin/(INT(vnp%dz/vnp%hmin)+1))
      ENDIF
      WRITE(*,*)'minimum receiver height ',vnp%hmin

      IF (vnp%nzout.GT.mxzout) THEN
         WRITE(*,*)'vnp.nzout,mxzout',vnp%nzout,mxzout
         STOP 'Too many points in height !!!'
      ENDIF
      IF (vnp%nrout.GT.mxrout) THEN
         WRITE(*,*)'vnp.nzout,mxzout',vnp%nrout,mxrout
         STOP 'Too many points in range !!!'
      ENDIF


c---  initialize receiver depth (ONLY for SAGA)
      WRITE(*,*)' First five Receiver heights:'
      IF (vnp%nzout.EQ.1) THEN
         rdep(1)=vnp%hmin
      ELSE 
         DO i=1,vnp%nzout
            rdep(i)=vnp%hmin+(vnp%hmax-vnp%hmin)/(vnp%nzout-1)*(i-1)
         ENDDO
      ENDIF
      DO  i=1,MIN(vnp%nzout,5)
         WRITE(*,*)i,rdep(i)
      ENDDO
      WRITE(*,*)'Output height and range points ',
     1     vnp%nzout,vnp%nrout
      
      IF (itpem_opt(4).EQ.1) THEN
         WRITE(*,*)' Reading noise floor and clutter cross sections...'
         READ(1,*)znoise  
         WRITE(*,*) ' Noise floor: ',znoise 
c---  clutter cross-section reading:
         READ(1,*) Nno
         IF (nno.GT.Mno) STOP 'Too many clutter-cross-sect nodes'
         READ(1,*) (elenod(i),i=0,nno+1)
         READ(1,*) (ccs(i),i=1,nno)
         WRITE(*,*)' x-coor, clutter CS'
         DO i=1,nno
            WRITE(*,*)elenod(i),ccs(i)
            IF (elenod(i).LT.elenod(i-1) ) THEN
               WRITE(*,*) 'x coordinates for node',i,'must increase'
               STOP
            ENDIF
         ENDDO
         IF (elenod(nno+1).LT.elenod(nno) ) THEN
            WRITE(*,*) 'x coordinates for node',nno,'must increase'
            STOP
         ENDIF
         
      ENDIF

      READ(1,*) rf%nprof,rf%lvlep !This is a range-independent case 
      WRITE(*,*)'refractivity profile points (H,R): ',
     &     rf%nprof,rf%lvlep  
      lvlep_start=rf%lvlep
      IF (rf%lvlep.GT.mxlvls) THEN
         WRITE(*,*) 'STOPPING: maximun number of m-profile',
     1        'points exceeded'
         WRITE(*,*)'rf%lvlep,mxlvls',rf%lvlep,mxlvls
         STOP
      ENDIF
c     
c---  parametric refractivity profile
c     
      IF (itpem_opt(5).EQ.1)  THEN
         DO j=1, rf%nprof       !loop over profiles
            READ(1,*)rf%rngprof(j) !Range of profile j   
            WRITE(*,*)'reading profile at range:',rf%rngprof(j) !Range of profile j   
            READ(1,*) rp%base(j),rp%thick(j),rp%offset(j),
     1           rp%mdef(j),rp%maxref(j)
            READ(1,*)(rp%coef(i,j),i=1,mpcoef)

            WRITE(*,*)' base height:', rp%base(j)
            WRITE(*,*)' thickness:', rp%thick(j)
            WRITE(*,*)' offset', rp%offset(j)
            WRITE(*,*)' M-defecit:', rp%mdef(j)
            WRITE(*,*)' Maximum height:', rp%maxref(j)
            WRITE(*,*)' Coefficients:',(rp%coef(i,j),i=1,mpcoef)
         ENDDO
         IF (itpem_opt(6).EQ.1)  THEN
            READ(1,*)rp%npoly   ! number of polynomials
            WRITE(*,*)' Number of polynomials', rp%npoly
            READ(1,*)(rp%factor(j,1),j=1,rp%npoly) ! number of polynomials
            WRITE(*,*)' The factor to the polynomials',
     1           (rp%factor(j,1),j=1,rp%npoly)
         ENDIF
c     
c---  reading refractivity profile
c     
      ELSE
         DO j=1, rf%nprof       !loop over profiles
            READ(1,*)rf%rngprof(j) !Range of profile j   
            WRITE(*,*)'reading profile at range:',rf%rngprof(j) !Range of profile j   
            DO i=1,rf%lvlep
               READ(1,*)rf%hmsl(i,j),rf%refmsl(i,j)
c     WRITE(*,*)rf.hmsl(i,j),rf.refmsl(i,j)
            ENDDO
            WRITE(*,*)'Last refractivity point:'
            WRITE(*,*)rf%hmsl(rf%lvlep,j),rf%refmsl(rf%lvlep,j)
         ENDDO
      ENDIF
C     
c     check the environment
c     
      IF (rf%nprof .GT. 1) THEN
         DO j=1, rf%nprof
            IF (rf%rngprof(j).LT.0) THEN
               WRITE(*,*)'Rangeprofile', j, 'less than zero range' 
               STOP
            ENDIF
         ENDDO
         DO j=2, rf%nprof
            IF (rf%rngprof(j).LT.rf%rngprof(j-1)) THEN
               WRITE(*,*)'Rangeprofiles must be increasing' 
               STOP
            ENDIF
         ENDDO
         IF (rf%rngprof(1).GT.0) THEN
            WRITE(*,*)'First rangeprofile less than zero range' 
            STOP
         ENDIF
         IF (rf%rngprof(rf%nprof).GT.vnp%rmax) THEN
            WRITE(*,*)' Last Rangeprofile greater than maximum-range' 
            STOP
         ENDIF
      ENDIF 
C     
      READ(1,*)tr%itp           !number of pairs in terrain profile 
      WRITE(*,*)'terrain points ',tr%itp 
      tr%itpmax=tr%itp 
      DO i=1,tr%itp
         READ(1,*)tr%terx(i),tr%tery(i)
         WRITE(*,*)tr%terx(i),tr%tery(i)
      ENDDO

      nx=vnp%nrout
c     ncurv=nx
c     nbart=1

c*****************Copied from SUBROUTINE "read_ran_in" **************
c     WRITE(prtfil,*)' reading ranges in....'
c     IF (iopt(2).EQ.1 .OR. iopt(2).EQ.3
c     1    .OR.iopt(13).EQ.1) THEN           ! for using readdat2
c     201     READ(1,*)nx
c     READ(1,*,err=203)(xranges(i),i=1,nx)
c     WRITE(prtfil,*) 'number of ranges',nx
c     ELSE   ! for using readdata and for syntetic
c     pg aug00       IF (iopt(2) .EQ. 2) THEN 
c     202     READ(1,*,err=200) rng1,rng2
      rng1=vnp%rmin
      rng2=vnp%rmax
      WRITE(prtfil,*)'rng1,rng2', rng1,rng2
      IF (rng1.GT.rng2) STOP 'rng1> rng2'
      nx=vnp%nrout              !=min(100,mx/nfrq)
c     IF ((iopt(2).EQ.0) .OR.(iopt(2).EQ.2)) THEN  ! for synthetic data
      drng= vnp%dr              !(rng2-rng1)/nx
c     nx=(rng2-rng1)/drng+1
c     nx=MIN(nx,mx/nfrq)
      DO i=1,nx
         xranges(i)=rng1+drng*(i-1)
      ENDDO
c     pg aug00        ENDIF
      ndep=vnp%nzout
c     ncurv=nx
c     nbart=1

c      WRITE(*,*)' RANGES',(xranges(i),i=1,nx)
      ncurv=nfrq*ndep
c---   for a vertical array of sensors.
       if (iopt(2)==3 .or. iopt(13)==1)       ncurv=nfrq*nx

c     pg aug00      ENDIF
      GOTO 205
 200  WRITE(*,*)' stop: range input section had wrong format'
      STOP
c     BACKSPACE(1)
c     GOTO 201
 203  WRITE(*,*)' stop: range input section had wrong format'
      STOP
c     BACKSPACE(1)
c     BACKSPACE(1)
c     GOTO 202
 205  CONTINUE
c*****************Copied from SUBROUTINE "read_ran_in" **************
      IF (nx.GT.mx) THEN
         WRITE(*,*) 'nx, mx=',nx,mx
         STOP 'nx is greater than mx'
      ENDIF

c---  should the  EOF be READ ? 

      IF (iopt(17).EQ.1) THEN
         CALL eofinit
      ENDIF

c***  create text strings for physical parameters
c     123456789012345678901234567890
      phystxt(1) = 'Refractivity $'
      phystxt(2) = 'Heights of refractivity point (m)$'
      phystxt(3) = 'Terrain x-coordinates (m)$'
      phystxt(4) = 'Terrain y-coordinates (m)$'
      phystxt(5) = 'Range (m)$'
      phystxt(6) = 'Height (m)$'
      phystxt(9) = 'Range (km)$'
      phystxt(11)= 'Shape coefficient$'
      phystxt(12)= 'Baseheight (m)$'
      phystxt(13)= 'Thickness (m)$'
      phystxt(14)= 'Offset (M-units)$'
      phystxt(15)= 'M-defecit (M-units)$'
      phystxt(16)= 'Noise-floor$'
      phystxt(17)= 'Slopes ......$'
      phystxt(18)= 'Clutter cross section (dB)$'
      phystxt(19)= 'Range factors   baseheight$'
c***  

      DO 8 j=1,mphys
         phystxt2(j)='                                               '
 6       DO 7 I=40,1,-1
            IF(phystxt(j)(I:I).NE.'$') GO TO 7
            phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c     WRITE(*,*)phystxt2(j)
            GO TO 8
 7       CONTINUE
    8 CONTINUE
      
c---- READ the optimization PARAMETER
      write(*,*)'entering readoptparm'
      CALL readoptparm
      write(*,*)'end readoptparm'
      

      DO i=1,nparm
         IF (fmin(i).GT.fmax(i))THEN
            WRITE(*,*)' *** fmin > fmax for parm',i
            ierr=1
         ENDIF
         IF (par2phy(i).LT.1 .OR. par2phy(i).GT.19 .OR.
     &        par2phy(i).EQ.10 )THEN
            WRITE(*,*)' *** par2phy not correct, parm',i
            ierr=1
         ENDIF
         IF (par2phy(i).EQ.11 )THEN
            IF (par2lay(i).GT.neof) THEN
               WRITE(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** The shapecoffient number is not defined'
               WRITE(*,*)' par2lay(i)', par2lay(i)
               ierr=1
            ENDIF
         ENDIF
c     IF (par2phy(i).EQ.1 )THEN
c     IF (fmin(i).LT.z0(nd0-1)) THEN
c     ENDIF
c     ENDIF
c     IF (par2phy(i).EQ.17)THEN
c     ENDIF

      ENDDO

      write(*,*)'... end of input subroutine'
c***  errors ?
      IF (ierr.EQ.1)STOP 'stopped in input'
c***  CLOSE input unit
      CLOSE(1)                                    

      END
c**********************************************
      SUBROUTINE tpemoption(options)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      CHARACTER options(40)
      INTEGER i
c     tilt=.FALSE.
      DO i=1,40
         itpem_opt(i)=0
      ENDDO 
c     itpem_opt(6)=1 range dependent steps.
      DO 50 i=1,40
         IF (options(i).EQ.'u') THEN
            itpem_opt(1)=1
            WRITE(*,*)'User defined stepsizes'
            WRITE(prtfil,*)'User defined stepsizes'
         ELSE IF (options(i).EQ.'E') THEN
            itpem_opt(2)=1
            WRITE(*,*)' Using TEDS SEOF'
            WRITE(prtfil,*)'Using TEDS SEOF '
         ELSEIF (options(i).EQ.'m') THEN
            itpem_opt(3)=1
            WRITE(*,*)'magnitude of field'
            WRITE(prtfil,*)'magnitude of field'
         ELSE IF (options(i).EQ.'c') THEN
            itpem_opt(4)=1
            WRITE(*,*)' Using clutter power'
            WRITE(prtfil,*)'Using clutter power'
         ELSEIF (options(i).EQ.'p') THEN
            itpem_opt(5)=1
            WRITE(*,*)'parametric  refrac profile'
            WRITE(prtfil,*)'parametric refrac. profile'
         ELSEIF (options(i).EQ.'a') THEN
            itpem_opt(6)=1
            WRITE(*,*)'range dependent polynomial'
            WRITE(prtfil,*)'range dependent polynomial'
         ELSEIF (options(i).EQ.'d') THEN
            itpem_opt(7)=1
            WRITE(*,*)' Multiple beam inversion'
            WRITE(prtfil,*) ' Multiple beam inversion'
         ELSEIF (options(i).EQ.'r') THEN
            itpem_opt(8)=1
            WRITE(*,*)' maximum M-deffecit'
            WRITE(prtfil,*) '  maximum M-deffecit '
         ELSEIF (options(i).EQ.'e') THEN
            itpem_opt(9)=1
            WRITE(*,*)' Meteorologic constrain inversion'
            WRITE(prtfil,*) ' Meteorologic constrain inversion'

         ELSE IF (options(I).EQ.'!') THEN
            GOTO 60
         ELSE IF (options(I).NE.' ') THEN
            WRITE(prtfil,399) options(I)
 399        FORMAT(1H ,' >>>> unknown TPEM OPTION: ',A1,' <<<<')
         END IF
 50   CONTINUE
 60   CONTINUE
      WRITE(*     ,*)'Tpem options: ',(itpem_opt(i),i=1,10)
      WRITE(prtfil,*)'Tpem options: ',(itpem_opt(i),i=1,10)
      END
c*******************
      SUBROUTINE bugfind(aaa)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'     
      CHARACTER*4 aaa
      
      WRITE(*,*)aaa
      WRITE(*,*)'ndep,frq(1),frq(2)',ndep,frq(1),frq(2)

      END
      
c**********************************************************************
      SUBROUTINE forwinit
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'      
      INTEGER i,j,id,jid,ierror,index
      INTEGER jstart,jend, i_valid
      REAL hminter
      INTEGER mloss(mxzout) 
      COMPLEX c_field(mxzout,mxrout)
      REAL*8 rout
      LOGICAL user_spec_dr
      INTEGER irejt,iaccept,ifrq,iperror
      SAVE irejt,iaccept
      INTEGER luttrf        ! logical unit for trf
      data luttrf/16/
      REAL offset
      REAL base_par(4)          ! used with shifted eof
      REAL realresp(1000),xmove !refmin
      INTEGER iresp(1000),ii,ir,idum !,iref
      INTEGER iprof
      REAL mexcess,metcon,basehvar(10000)
      COMMON/ mexc/ mexcess,metcon,basehvar,iprof
      iperror=0
      WRITE(*,*)' uses TPEM as forward model'                     
      vnp%propang = 0.          !Automatic internal calculation.       
c     
c---  pass number of ranges back to tpem
c     
      IF (iopt(2).EQ.2) THEN
         vnp%nrout=nx
         IF (vnp%nrout.GT.mxrout) THEN
            WRITE(*,*)'vnp%nzout,mxzout',vnp%nrout,mxrout
            STOP 'Too many points in range !!!'
         ENDIF
      ENDIF

      IF (itpem_opt(4).EQ.1) THEN
         id=1                   ! only one depth
         ifrq=1                 ! only one frequency
         IF (ierrinfile.EQ.0) THEN
            data_mean=0
            DO ifrq=1,nfrq
               index=(id +((ifrq-1))*ndep-1)*nx
               DO i = 1,  vnp%nrout !ranges
                  data_mean=data_mean+DATA(i+index)
               ENDDO
            ENDDO
            data_mean = data_mean/ vnp%nrout/nfrq
         ELSE
            data_mean=100
         ENDIF
      ENDIF
c     Set LOGICAL flags to trap for errors.  LERR6=.TRUE.-PEINIT returns error
c     IF last point in terrain profile is less than maximum plot range.
c     LERR12=.TRUE.-PEINIT RETURN error IF last refractivity profile entered
c     (for range-dependent environment) is less than maximum plot range.
      
      ef%lerr6 = .TRUE.
      ef%lerr12 = .TRUE.
      flagpu=-2
      IF (iopt(6).EQ.1) flagpu=-10
      IF (dr.GT.0) THEN
         user_spec_dr=.TRUE.
      ELSE
         user_spec_dr=.FALSE.
      ENDIF

c****************************************************
      ENTRY forw2 
c***************************************************

      flagpu=1 +flagpu
      iprof=0
      mexcess=0
      metcon=1000
      IF (iWriteTrf.EQ.1)THEN
         vnp%nzout= 100 
         ndep= vnp%nzout
         vnp%hmin = 0 
         vnp%hmax = 200
         vnp%dz  =(vnp%hmax-vnp%hmin)/vnp%nzout
         vnp%hmin= vnp%dz           
         vnp%nrout=75
         vnp%rmin=0
         vnp%rmax=100000        !Range up to 100 km.      
         vnp%nrout=100
         vnp%rmin=00
         vnp%rmax=60000         !Range up to 100 km.      
         nfrq=1
         vnp%dr=(vnp%rmax-vnp%rmin)/vnp%nrout
         vnp%rmin= vnp%dr
         CALL TRFHEAD('TRF   ',TITLE,REAL(vnp%hmax),REAL(vnp%hmin),
     .        REAL(vnp%rmin),REAL(vnp%dr),sv%antht)                
      ENDIF


      DO 100 ifrq=1,nfrq 
         IF (itpem_opt(2).EQ.1) THEN
            IF (flagpu.LT.2) THEN
               WRITE(*,*)' calling TEDs SEOF- generation' 
               WRITE(prtfil,*)'cap_coef',cap_coef(1),
     1              (cap_coef(1+i),i=1,ncap)
            ENDIF
            CALL seofgen(cap_coef,Ncap, rf%hmsl(1,1),
     1           rf%refmsl(1,1),lvlep_start,i_valid,base_par)
            IF (i_valid.EQ.0) THEN ! this is not a valid profile
               irejt=irejt+1
               WRITE(*,*)'profile rejected, rejt/accpt=',irejt,iaccept
c     WRITE(99,*)' irejt',irejt
c     DO j=1, rf.nprof     !loop over profiles
c     DO i=1,rf.lvlep
c     WRITE(99,*)i,j,rf.hmsl(i,j),rf.refmsl(i,j)
c     ENDDO
c     ENDDO
               lwascomp=0
               RETURN
            ELSE
               iaccept=iaccept+1
            ENDIF
         ENDIF
         IF (itpem_opt(7).EQ.1) THEN
            sv%elev=beamelev(ifrq)
         ELSE
            sv%freq=frq(ifrq)
         ENDIF
         IF (flagpu.LT.0) WRITE(*,*)'computing for freq',sv%freq
         IF (.NOT. user_spec_dr ) THEN
            dr=0       
         ELSE
c     WRITE(*,*)'dr is not changed'
         ENDIF
         lwascomp=1
c     WRITE(*,*) 'has entered tpem'
c---  
         tr%itp=tr%itpmax

c     
c---- USE of EOF
c     
         IF (iopt(17).EQ.1) THEN
            CALL eofval
            IF (lwascomp.EQ.-1) THEN
               RETURN
            ENDIF
         ENDIF

         IF (vnp%rmin.LE.0) THEN
            vnp%dr=(vnp%rmax-vnp%rmin)/vnp%nrout
            vnp%rmin= vnp%dr
         ELSE    
            IF (vnp%nrout.NE.1) THEN
               vnp%dr = (vnp%rmax-vnp%rmin)/(vnp%nrout-1)
            ELSE
               vnp%dr   =  vnp%rmax
               vnp%rmin =  vnp%rmax 
            ENDIF
            vnp%rmin =  vnp%rmin/(INT(vnp%dr/vnp%rmin)+1)
            vnp%rmin = MAX(vnp%dr,vnp%rmin)
         ENDIF
         rf%lvlep=lvlep_start
         

         IF (flagpu.LT.2) THEN
            CALL writetpem
         ENDIF                  ! debugging output
c         write(*,*)'hello'

         IF (itpem_opt(5).EQ.0) THEN
            DO j=1, rf%nprof    !loop over profiles
               DO i=1,rf%lvlep
                  IF ((rf%hmsl(i,j).LT.0).
     &                 OR.(rf%refmsl(i,j).LT.0)) THEN
                     WRITE(*,*)'point',i,j,' heigth and refractivity',
     1                    rf%hmsl(i,j),rf%refmsl(i,j)
                     STOP 
     1                    'tpeminit:  heigth and
     1                    refractivity cannot be negative'
                  ENDIF
               ENDDO
            ENDDO
         ENDIF

         
c     
c     Variables in CAPS are returned.
c      WRITE(*,*)'* Enter PEINIT *'
         
         CALL peinit( ef, vnp, rf, sv, tr,rp,
     1        HMINTER, ROUT, IERROR,itpem_opt,iperror)  

c      WRITE(*,*)' Exit PEINIT '

         IF( ierror .NE. 0 ) THEN
            WRITE(*,*)'*********** ERROR IN PEINIT **************'
            WRITE(*,*)'******** IERROR = ', ierror,' **********'
            STOP
         END IF
         IF (iperror.LT.0 ) THEN
            WRITE(*,*)'on return peinit:iperror',iperror
            lwascomp=-1
            RETURN
         ENDIF
c     WRITE(15,*)'Reference height in m = ', hminter
c*********************************************


         DO i = 1,  vnp%nrout   !ranges
            IF (iopt(6).EQ.1 .AND. flagpu.LT.2) THEN
               WRITE(*,*)'calling pestep, range=',i,rout,ndep   
            ENDIF   
            IF (iopt(2).EQ.2) THEN
               vnp%nrout=nx
            ENDIF
c     CALL bugfind('bbbb')
            CALL pestep( hminter, vnp, rf, tr,rp,sv,
     +           ROUT, MLOSS,c_field(1,i),
     +           JSTART, JEND, fieldtype,itpem_opt,iperror)
c     ,xranges(i) ) 
c     CALL bugfind('bbb2')
c     WRITE(*,*)'called pestep, range=',i,iprof  
            IF (iperror.LT.0 ) THEN
c     WRITE(*,*)'on return pestep:iperror',iperror
               lwascomp=-1
               RETURN
            ENDIF

c     Recall that MLOSS is the propagation loss in centibels, i.e., 
c     MLOSS() = NINT( propagation loss in dB * 10. ).  JSTART = start of
c     valid loss points, JEND = END of valid loss points.
            IF ((flagpu.LT.0) .AND. (ifrq.EQ.1)) THEN
               WRITE(15,*)
               WRITE(15,*)'range in km = ', rout*1.e-3
c     WRITE( *,*)'range in km = ', rout*1.e-3  !Output to screen
               IF (fieldtype.EQ.1) THEN ! complex field  
                  DO j = jstart, jend
                     WRITE(15,*)'Height (m)= ',j*vnp%dz+vnp%hmin,
     1                    ' c_field ',c_field(j,i)
                  END DO                               
               ELSE
                  DO j = jstart, jend
                     WRITE(15,*)'Height (m)= ',j**vnp%dz+vnp%hmin,
     1                    '  loss in dB = ',mloss(j)*.1
                  END DO                               
               ENDIF
            ENDIF
            IF ((jend-jstart+1).NE.ndep)THEN
               WRITE(*,*)'************** ERROR *********************'  
               WRITE(*,*)' Number of calculated depth responses'
               WRITE(*,*)' do not correspond to number of requested'
               WRITE(*,*)' message from tpeminit'
               WRITE(*,*)'jstart,jend,ndep=',jstart,jend,ndep
            ENDIF
c     
c---- writing TRF file
c     
            IF ((iWriteTrf.EQ.1) ) THEN
               IF (fieldtype.EQ.0) THEN
                  DO j = jend,1,-1
                     WRITE(LUTTRF,*) 10**(-mloss(j)*0.005),
     &                    0
c     WRITE(*,*)'writing to unit',luttrf
                  ENDDO 
               ELSE
                  DO j = jend,1,-1
                     c_field(j,i)=CONJG(c_field(j,i))
                     WRITE(LUTTRF,*) REAL(c_field(j,i)),
     &                    imag(c_field(j,i))
c     WRITE(*,*)'writing to unit',luttrf
                  ENDDO 
               ENDIF
            ENDIF

            DO id=1,ndep        ! irin
               jid=id+jstart-1
               index=(id +((ifrq-1))*ndep-1)*nx
               IF (fieldtype.EQ.0) THEN ! centibel field      
                  resp(i+index)=-mloss(jid)*.1
c     WRITE(98,*)'orig',i,id,ifrq,resp(i+index)
c     WRITE(*,*)'fieldtype',fieldtype,ndep,-mloss(jid)
c     WRITE(*,*)'resp', jid,i,index, resp(i+index)
                  mloss(jid)=0
               ELSE             ! complex field     
                  resp(i+index)=c_field(jid,i)
               ENDIF
c     
c---  Clutter model.
c     
               IF (itpem_opt(4).EQ.1) THEN
                  IF (fieldtype.EQ.0) THEN  
                     resp(i+index)=2*resp(i+index)
     1                    +10*LOG10(rout/1000)
                  ELSE          ! complex field     
                     resp(i+index)=40*LOG10(ABS(resp(i+index)))
     1                    +10*LOG10(rout/1000)
                  ENDIF
c      WRITE(*,*)'respclutter', i,index, resp(i+index),rout
c     WRITE(*,*)'resp,range', resp(i+index),rout
c---- clutter power
                  CALL addclutterpower(i+index,rout)            
                  IF (REAL(resp(i+index)).LE.znoise)THEN
                     resp(i+index)=znoise 
                  ENDIF       
c     WRITE(*,*)'resp clutter', resp(i+index),rout
               ENDIF            ! end clutter

               IF (flagpu.LT.3) THEN
                  WRITE(prtfil,*)rdep(id),rout,resp(i+index),id,i+index
               ENDIF
            ENDDO               ! ndep loop.                              
         END DO                 ! range_loop
 100  CONTINUE                  !frequency loop

      IF (itpem_opt(4).EQ.1) THEN
         id=1                   !only one depth
         jid=id+jstart-1
         resp_mean=0
         DO ifrq=1,nfrq
            index=(id +((ifrq-1))*ndep-1)*nx
            DO i = 1,  vnp%nrout !ranges
               resp_mean=resp_mean+resp(i+index)
            ENDDO
         ENDDO
         resp_mean = resp_mean/ vnp%nrout/nfrq
         offset    = data_mean - resp_mean
c      WRITE(*,*)'resp_mean,data_mean, offset'
c       WRITE(*,*) resp_mean,data_mean,offset
         DO ifrq=1,nfrq
            index=(id +((ifrq-1))*ndep-1)*nx
            DO i = 1,  vnp%nrout !ranges
               resp(i+index)=resp(i+index)+offset
               realresp(i+vnp%nrout*(ifrq-1))=REAL( resp(i+index))
c     WRITE(*,*) 'realresp(i),i',realresp(i),i
            ENDDO
         ENDDO
c     
c-----truncating....
c     
         idum=1
         IF (idum.EQ.1) THEN

c     WRITE(*,*)' I should never be here'
c     pointers to sorted array iresp 
            CALL indexx( vnp%nrout*nfrq, realresp,iresp)
c     iresp is an array that points to the realresp array in the order 
c     in increasing order 
c     
            xmove=0
c     truncate points below zero
            DO i = 1,vnp%nrout*nfrq,1 !samples
c     WRITE(*,*) 'iresp(i),i',iresp(i),i
               IF ((realresp(iresp(i))-xmove).LT. 0) THEN
                  xmove=-(realresp(iresp(i))-xmove)
     1                 /(vnp%nrout*nfrq-i)+xmove
                  realresp(iresp(i))=0
               ELSE
                  ii=i
                  GOTO 79
               ENDIF
            ENDDO
            
 79         DO i = ii,vnp%nrout*nfrq,1 !ranges
c     WRITE(*,*)'i,ii,index,iresp(i)',i,iresp(i),xmove
               realresp(iresp(i))= realresp(iresp(i))-xmove
            ENDDO

            resp_mean=0
            ii=0
            DO ifrq=1,nfrq
               index=(id +((ifrq-1))*ndep-1)*nx
               DO i = 1,  vnp%nrout !ranges
                  ii=ii+1
c     resp(i+index)= realresp(iresp(ii))
                  resp((ii))= realresp(ii)
                  resp_mean=resp_mean+resp(iresp(ii))
c     WRITE(98,*)'resp',i,id,ifrq,resp(ii)
c     WRITE(*,*)'iresp(ii),i+index',iresp(ii),i+index
               ENDDO
            ENDDO
            resp_mean = resp_mean/ vnp%nrout/nfrq
c     WRITE(*,*)'new resp mean', resp_mean
c     DO i = 1,vnp.nrout,1   !ranges
c     WRITE(*,*)'clutter',resp(i+index),id,i+index
c     ENDDO
         ENDIF
c     
c-----truncating ended....
c     
      ENDIF



c     
c---  From post writeout refractivity profile
c     
      IF ((iWriteTrf.EQ.1) ) THEN
         CALL writeMfile
      ENDIF
      
      END                                            

      SUBROUTINE plotstddev(expfit,totobs)
c     computes standard deviation... genetic algorithm
c     optimization PROGRAM
c     PETER GERSTOFT, 1992
c     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
c     REAL fval              ! function for computation of real values
      INTEGER i,j,iq
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      REAL xsed(10,mq_post),xmeansed(10),svarsed(10)
      REAL*8 xsumsed(10),xsumsqrsed(10),zsedi(10) ! param for statis.
      REAL*8 expfit(mq_post),totobs,xhelp
c     REAL   xbestbest(meofvar),xmeanmean(meofvar)
      INTEGER npoint,nbase_par,i_valid
      REAL base_par(4)

      WRITE(*,*)' std dev for mean model...'

      nbase_par=4
      npoint=4
      DO i=1,npoint
         xsumsed(i)=0.
         xsumsqrsed(i)=0.
      ENDDO
      
c      OPEN(unit=8,file='velppd.m',access='append',
c     &     status='unknown')
      WRITE(8,*)'velsed=['
      
      DO iq=1,q 
         CALL setmodel(iq)
         IF (itpem_opt(2).EQ.1) THEN
            CALL seofgen(cap_coef,Ncap, rf%hmsl(1,1),
     1           rf%refmsl(1,1),lvlep_start,i_valid,base_par)
            IF (i_valid.EQ.0) THEN ! this is not a valid profile
               WRITE(*,*)'profile not valid, iq=',iq
            ENDIF
         ELSE
            IF (iopt(17).EQ.1) THEN
               CALL eofval
            ENDIF
            base_par(1) = rf%hmsl(2,1) ! base height
            base_par(2) = rf%refmsl(2,1)-rf%refmsl(3,1) ! mdeffecit
            base_par(3) = rf%hmsl(3,1)-   rf%hmsl(2,1) ! 
            base_par(4) = rf%refmsl(4,1) ! offset
         ENDIF
         DO i=1,nbase_par
            j=i
            xsed(i,iq)=base_par(i)
            xhelp = xsed(j,iq)*expfit(iq)
            xsumsed(j) = xsumsed(j)+xhelp
            xsumsqrsed(j)=xsumsqrsed(j)+xhelp*xsed(j,iq)
         ENDDO
         WRITE(8,'(i3,f8.4,10f10.2)')iq,expfit(iq), 
     1        (xsed(j,iq),j=1,nbase_par)
      ENDDO                     !number of populations, iq
      WRITE(8,*)'];'
      CLOSE(8)
c     c
c---- for standard deviation
c     
      DO i=1,npoint         
         xmeansed(i)=xsumsed(i)/totobs ! /2 
         svarsed(i)
c     abs in formula below is to avoid problems...
     &        =SQRT(ABS(xsumsqrsed(i)/totobs-(xsumsed(i)/totobs)**2))
c     IF (xmean(i).NE.0)svar(i)=svar(i)/ABS(xmean(i))
         WRITE(*,*)'mean,std.dev', xmeansed(i), svarsed(i)
      ENDDO
      zsedi(1)=0
      zsedi(2)=5
      zsedi(3)=10
      zsedi(4)=20
c      OPEN(unit=8,file='vel.m',access='append',
c     &     status='unknown')
      
      WRITE(8,*)' vel=['
      DO i=1,npoint
         WRITE(8,*)zsedi(i),xmeansed(i)-svarsed(i),xmeansed(i),
     &        xmeansed(i)+svarsed(i)
      ENDDO
      WRITE(8,*)' ];'
      END
c***************************************
      SUBROUTINE addclutterpower(ipoint,rout)
c     -    This adds the clutter power  of the response.
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'
      INTEGER i,ipoint
      REAL   xtemp,xtemp2,xtemp3
      REAL*8 rout

      DO i=1,nno
         xtemp=(rout-elenod(i-1))
     1        *(rout-elenod(i+1))
c      WRITE(*,*)'temp1:',  elenod(i),rout,
c     1                xtemp,ccs(i)
         IF (xtemp.LT.0) THEN
            xtemp3=rout-elenod(i)
            IF (xtemp3.LT.0) THEN !left  element
               xtemp2=1-xtemp3/(elenod(i-1)-elenod(i))
            ELSE                !right element
               xtemp2=1-xtemp3/(elenod(i+1)-elenod(i))
            ENDIF
            resp(ipoint) =resp(ipoint)+ccs(i)*xtemp2
         ENDIF
      ENDDO
c     WRITE(*,*)' from cluttercs: resp', resp(ix+index)

      END

      SUBROUTINE writeMfile
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'      
      REAL wl,fko,delz,zmax,con
      INTEGER  n, ln, n34
      COMMON / pevar / wl, fko, delz, n, ln, zmax, n34, con
c     ONLY used WITH generating profile
      REAL rv2, refdum(mxlvls), htdum(mxlvls),
     +     profint(0:maxpts), ht(0:maxpts)
      INTEGER  is
      COMMON / parinit / rv2, refdum, htdum,
     +     profint, ht, is
c     local vars
      INTEGER i
c     OPEN(unit=9,file='res.m',status='unknown',access='append')
      WRITE(13,'(a)') 'vel=['
      DO i=1,N
         WRITE(13,'(2f12.3)')  ht(i),profint(i)/con 
      ENDDO
      WRITE(13,'(a)') '];'
      IF (itpem_opt(4).EQ.1) THEN
         WRITE(13,'(a)') 'clut=['
         WRITE(13,'(2f12.3)')elenod(0)/1000,0
         DO i=1, nno
            WRITE(13,'(2f12.3)')elenod(i)/1000,ccs(i)
         ENDDO
         WRITE(13,'(2f12.3)')elenod(nno+1)/1000,0
         WRITE(13,'(a)') '];'
      ENDIF
      WRITE(*,*)'... the best found refractivity', 
     1     ' profile has been written to file res.m'
      END
c**********************

      SUBROUTINE writetpem     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './tpem/tpem.inc'
      INCLUDE 'comtpem.h'      
      INTEGER i,j
      INTEGER luttrf       ! logical unit for trf
      data luttrf/16/
      WRITE(prtfil,*)
      WRITE(prtfil,*)' Calling TPEM' 
      WRITE(prtfil,*)' flagpu',flagpu
      IF ( fieldtype.EQ.0) THEN
         WRITE(*,*) ' PROPAGATION LOSS'                  
      ELSE IF ( fieldtype.EQ.1) THEN
         WRITE(*,*) ' COMPLEX FIELD'                  
      ELSE
         WRITE(*,*) ' FIELDTYPE NOT DEFINED'
         STOP        
      ENDIF
      WRITE(prtfil,*)' Frequency ',sv%freq !Frequency 
      WRITE(prtfil,*)' antenna height ',sv%antht
      WRITE(prtfil,*)' antenna pattern ',sv%ipat
      WRITE(prtfil,*)' antenna bwidth ',sv%bwidth
      WRITE(prtfil,*)' antenna elevation ',sv%elev
      WRITE(prtfil,*)' min height ',vnp%hmin
      WRITE(prtfil,*)' max height ',vnp%hmax
      WRITE(prtfil,*)' max range ',vnp%rmax 
c     
c---  parametric refractivity profile
c     
      IF (itpem_opt(5).EQ.1)  THEN
         DO j=1, rf%nprof       !loop over profiles
            WRITE(*,*)'profile at range:',rf%rngprof(j) !Range of profile j   
            
            WRITE(*,*)'  base height:   ', rp%base(j)
            WRITE(*,*)'  thickness:     ', rp%thick(j)
            WRITE(*,*)'  offset         ', rp%offset(j)
            WRITE(*,*)'  M-defecit:    ', rp%mdef(j)
            WRITE(*,*)'  Maximum height:', rp%maxref(j)
            WRITE(*,*)' Coefficients:',(rp%coef(i,j),i=1,mpcoef)
         ENDDO
         IF (itpem_opt(5).EQ.1)  THEN
            WRITE(*,*)' The factor to the polynomials',
     1           (rp%factor(j,1),j=1,rp%npoly)
            
         ENDIF     
c     
c---  reading refractivity profile
c     
      ELSE
         DO j=1, rf%nprof       !loop over profiles
            WRITE(prtfil,*)'profile at range:',rf%rngprof(j) !Range of profile j   
            DO i=1,rf%lvlep
               WRITE(prtfil,*)rf%hmsl(i,j),rf%refmsl(i,j)
            ENDDO
            WRITE(prtfil,*)'Finished with profile',rf%lvlep
         ENDDO
      ENDIF
      IF (tr%itp .GT. 0) THEN
         WRITE(prtfil,*)'terrain points ',tr%itp 
         WRITE(prtfil,*)'   x    y'
         DO i=1,tr%itp
            WRITE(prtfil,*)tr%terx(i),tr%tery(i)
         ENDDO
      ENDIF
      IF (itpem_opt(4).EQ.1) THEN
         WRITE(*,*)'clutter '
         DO i=1, nno
            WRITE(*,*)elenod(i),ccs(i)
         ENDDO
      ENDIF
      END


