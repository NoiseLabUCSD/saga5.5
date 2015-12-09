      PROGRAM ga
c****************************************************************************
c     PROGRAM for inversion of seismo-acoustic DATA using genetic algorithms.
c     The RESULT is presented as the most likely model, the mean model and the 
c     standard variation.
c
c     PETER GERSTOFT, 1992 
c****************************************************************************  

                                ! Common block definitions are
      USE global
      INCLUDE 'comopt.h'        ! INCLUDED in these files.
      INCLUDE 'comforw.h'
      REAL fval                 ! function for computation of real values
C**   local variables
      INTEGER iter,iq,iiq,i,J,  ! iterative indices.
     >     ifreq                ! frequency index.
      INTEGER ipar(mq),iforwloc
      INTEGER ihelp,lun,index
      INTEGER modelball(mpar)   ! best-of-all model
      INTEGER modelbsf(mpar)    ! best-so-far model
c     counters forward-models
      INTEGER iforwtot,iforwpop,iforwtotal,iforwreject 
      REAL    objbsf,objball
      data    objbsf /1e10/,objball /1e10/
      REAL t1,ahelp,fit_teo
      real fitpar(mq),expfit(mq),
     >     pmodel(mq),pcum(mq),
     >     expobj,temp
      REAL Xstart(mpar)
      REAL ran2,xtheta(mpar)
      REAL nsearch,xdum
      INTEGER ix,itrn,fork,ispawn,wait ,retval,status,nspawn,nspawn0
      INTEGER ifrq,idep
      EXTERNAL fork,wait
      INTEGER getpid,id,iterpowell
      CHARACTER*20 COMFILE
      INTEGER qpar,irep,qtemp,iparm,maxiter,
c     error flags for read_cov and read_hp
     >     ierrdat,ierrcov,ierrhp,ierraun,
     >     ierrinp,ioer,ierr,ierrD
      data qtemp/1/,maxiter/100/,
c     error flags for read_cov and read_hp
     >     ierrdat/0/,ierrcov/0/,ierrhp/0/,ierraun/0/, 
     >     ierrinp/0/,ioer/0/,ierr/0/,ierrD/0/
c     INTEGER iter,iforwpop
      COMMON /iterpar/iter,iforwpop 
      INTEGER flagpu
      COMMON /flagpu/flagpu
      INTEGER jparm,M,itmax,mq_post1
c     PARAMETER(M=50) 
      REAL, DIMENSION(:,:), ALLOCATABLE :: xiun
      real fret,tol,xstep(mpar)
      REAL xlocobj,xlocmod(mpar)
      REAL, DIMENSION(:), ALLOCATABLE ::  fit
      WRITE(*,*)'Running SAGA version 5.5, September 2010' 
      xlocobj=1e10
      IF (iopt(34).EQ.1) THEN
         iseed=-getpid()-ipop   !  random seed
      ELSE
         iseed=-6511
      ENDIF
      xdum=ran2(iseed)          ! initialize random generator
c     
c---- initialization
c
      iWriteTrf=0
c      OPEN(6, STATUS='NEW',CARRIAGECONTROL='LIST')
c
      CALL opfilr(1,ierrinp)
      IF (ierrinp.NE.0) THEN
         WRITE(*,*)' The input file *.dat does not exist'
         STOP
      ENDIF
      CALL opfilw(7,ioer)
c
c---  Set some forward model specific flags
c
      CALL forwardmodel(iopt,mopt)
c
c---  READ input DATA
c
      CALL input 
      CALL flush(7)
      nsearch=0
      DO iparm=1,nparm
         nsearch=nsearch+LOG10(1.0*ndigit(Iparm))
      ENDDO
      WRITE(*,'(a,f6.2)')' Size of the search space: 10**',nsearch
      WRITE(prtfil,'(a,f6.2)')' Size of the search space: 10**',nsearch

c---  allocate
      mq_post1=(niter+q)*npop
      allocate(model(nparm+2,mq_post1))
      allocate(allmodel(nparm+2,mq_post1))
      allocate(fit(mq_post1))
      allocate(allfit(mq_post1))
      allocate(parents(nparm,q))
      if (itrans(4).gt.0) then 
         write(*,*) 'allocating weight matrix'
         allocate(weight(mobs))   
            weight=1.0
      endif
      if (iopt(13).gt.0) then
         write(*,*) 'allocating cov matrix'
         allocate(cov(ndep*ndep*nfrq*nx))
      else
         write(*,*)'ndep,nfrq,mx',ndep,nfrq,mx,nx
         allocate(data(ndep*nfrq*nx))
      endif
C .cfh.
      if (iopt(35).gt.0) then
         write(*,*) 'allocating unc matrix'
         allocate(AUN(ndep*ndep*nfrq*nx))
      endif

c--- allocate for the powell method.
      M=nparm
      allocate(xiun(M,M))

C     Reading DATA from COV or IN file (unit 3 or 2). This can be in
C     covariance matrix form (option C), or it can be in hydrophone
C     spectra vectors (option d) or ? (D) or Saplot FORMAT (e).

      IF (iopt(2).GE.1) THEN    ! Options c, d, D, or e.
         CALL opfilr(2,ierrinfile)
         IF (ierrinfile.NE.0) THEN
            WRITE(*,*) 'No in-file  with observed read'      
         ELSE         
            IF (iopt(13).EQ.1) THEN ! Option c: Read covariance matrix
               WRITE(*,*)' Reading covariance matrix' ! from COV file.
               CALL READ_COV(ierrcov)
               
            ELSEIF (iopt(2).EQ.1) THEN ! Option d: Read h/p data vectors
               CALL readdat2(ierrdat) ! file in ? format.
            
            ELSEIF (iopt(2).EQ.2) THEN ! Option D: Read data from IN
               CALL readdata(ierrD) ! file in ? format.
               
            ELSEIF (iopt(2).EQ.3) THEN ! Option e: Read data from IN
               WRITE(*,*)' Reading hydrophone data' ! from IN file.
               CALL READ_HP(ierrhp)
c     CALL readdatiso           ! file in 'rev' format.
            ELSEIF (iopt(2).EQ.4) THEN ! Option T: Read data from IN
               WRITE(*,*)' Reading HA data' ! from IN file.
               CALL READ_HA(ierrhp)
c     CALL readdatiso           ! file in 
            ENDIF
            CLOSE(2)
         END IF
      END IF
      write(*,*)' ... in-file read'

C .cfh.****************************
      IF (iopt(35) .EQ. 1) THEN
         WRITE(*,*)' Reading unc. cov. matrix' ! from UNC file.
         CALL opfilr(4,ierraun)
         CALL READ_AUN(ierraun)
         CALL ROTMAT
         CLOSE(4)
      ENDIF

C**********************************

      IF (itrans(4).GE.1) THEN
         CALL opfilr(3,ioer)
         CALL Read_weight
         CLOSE(3)
      ENDIF
      if ((iopt(23)==1)) then
c     reading prior info
         CALL opfilr(31,ierrinfile)
         write(*,*)' reading prior information'
         call readprior
         close(31)
      endif
C     Run the forward model using the parameters specified in the run
C     (.DAT) file.  The results of this run can later be written to the
C     OBS file in either covariance matrix or hydrophone DATA FORMAT.
      CALL forwinit
       
c
c---- should DATA be weighted
c
      IF (iopt(3).GE.1) THEN
c     WRITE(*,*)' the DATA is being weighted...'
         CALL norm
      ENDIF   

      iforwtot=1

      CALL getmodelreal(xstart)
      CALL modeldecode(xstart,1,ierr)
      DO iparm=1,nparm
         xstar(iparm)=xstart(iparm)
         xprior(iparm)=xstart(iparm)
      ENDDO

C     Option W: Writes forward model output to the OBS file. DATA can
C     be in Covariance Matrix (option c), or h/p DATA vector (d) FORMAT.
      IF (iopt(21).EQ.1) THEN   ! option W.
         CALL opfilw(30,ioer)
         IF (iopt(13).EQ.1) THEN ! Option c (Covariance matrix)
            WRITE(30,'(a1,a31)') '!',' File format: Covariance matrix' 
            CALL header(30)
            WRITE(*,*) ' The generated cov-matrix will be',
     >           ' written to the OBS file'
            DO ifreq = 1, nfrq
               DO ix = 1, nx
                  CALL WRITE_COV(ifreq,ix)
               ENDDO
            ENDDO
c---  reread covariance matrix
            CALL opfilr(2,ierrinp)
            IF (ierrinp.EQ.0) THEN
               WRITE(*,*)' Reading covariance matrix' ! from COV file.
               CALL READ_COV(ierrcov)
               CLOSE(2)
            ENDIF
         ELSE
            IF (iopt(2).EQ.3) THEN ! Option e (hydrophone data vectors)
               WRITE(30,'(a1,a32)') 
     >              '!',' File format: Hydrophone vectors' 
               CALL header(30)
               WRITE(*,*) ' The generated h/p data will be',
     >              ' written to the OBS file'       
               DO ix = 1, nx
                  DO ifreq = 1, nfrq
                     CALL WRITE_HP(ifreq,ix)
                  ENDDO
               ENDDO
            ELSEIF (iopt(2).EQ.4) THEN ! Option T (HA data vectors)
               WRITE(30,'(a1,a32)') 
     >              '!',' File format: HA array' 
               CALL header(30)
               WRITE(*,*) ' The HA data will be',
     >              ' written to the OBS file'       
               DO idep = 1, ndep
                  DO ifreq = 1, nfrq
                     CALL WRITE_HA(ifreq,idep)
                  ENDDO
               ENDDO
            ELSEIF (iopt(2).EQ.1) THEN ! option d
               WRITE(30,'(a1,a)')'!',title
               CALL header(30)
               
               WRITE(*,*)'ndep,nx',ndep,nx,iopt(2)
               DO ifrq=1,nfrq
                  DO idep=1,ndep
c     DO j=1,ncurv
                     j=(ifrq-1)*ndep +idep
                     index=(j-1)*nx
                     DO i=1,nx
                        WRITE(30,*)ifrq,idep,i,resp(i+index),index+i
                     ENDDO      ! nx
                  ENDDO         ! ndep
               ENDDO            ! nfrq
            ELSEIF (iopt(2).EQ.2) THEN
c     DO j=1,ndep
               WRITE(30,'(a1,a)')'!',title
               CALL header(30)
               
               index=(j-1)*nx
               DO i=1,nx
                  WRITE(30,'(100G13.6)')
     &                 xranges(i),(REAL(resp(i+(j-1)*nx)),j=1,ncurv)
               ENDDO
c     ENDDO
               
            ELSE
               WRITE(*,*)'Error: Output format for option W',
     &              'not defined.'
               STOP
            ENDIF
         ENDIF
         CLOSE(30)
         WRITE(*,*)' The synthetic generated data has been written',
     &        ' to the *.obs file'
      ENDIF
      
c     IF (iopt(2).EQ.0.AND.iopt(13).EQ.0) THEN
c---  move DATA from resp to DATA  for synthetic inversions
c     DO j=1,ncurv
c     index=(j-1)*nx
c     DO i=1,nx
c     DATA(i+index)=resp(i+index)
c     c          ENDDO
c     ENDDO
c     ENDIF
c     
c---- should DATA be weighted
c
      IF (iopt(3).GE.1) THEN
         WRITE(*,*)' Data is being weighted...'
         CALL normdata
      ENDIF   
c--- for response
      IF (iopt(25).EQ.1) THEN
         WRITE(*,*)' transforming data: only magnitude is used'
         DO j=1,ncurv
            index=(j-1)*nx
            DO i=1,nx
c     WRITE(prtfil,*)' data',DATA(j,idep),ABS(DATA(j,idep))
               DATA(i+index)=ABS(DATA(i+index))
            ENDDO
         ENDDO
      ENDIF
      IF (iopt(27).EQ.1) THEN
         WRITE(*,*)' transforming data: to unit norm'
         ahelp=0
         DO j=1,ncurv
            index=(j-1)*nx
            DO i=1,nx
               ahelp=ahelp+ABS(DATA(i+index))**2
            ENDDO
         ENDDO
         ahelp=1.0/SQRT(ahelp)
         DO j=1,ncurv
            index=(j-1)*nx
            DO i=1,nx
               DATA(i+index)=(DATA(i+index))*ahelp
            ENDDO
         ENDDO
      ENDIF
c     
c---  energy normalization
c   
c      enorm=0.
c      DO j=1,ncurv
c        index=(j-1)*nx
c         DO i=1,nx
c            enorm=enorm+ABS(DATA(i+index))**2
c     WRITE(*,*)'from ga data',j,i,DATA(i+index)
c         ENDDO
c      ENDDO
c
c---  should noise be added ??? 
c
c      IF (iopt(11).EQ.1) THEN
c        WRITE(*,*)' ***** warning noise is being added ********'  
c        xnoise=SQRT(enorm)/5./SQRT(2.)
c        xnoise=0.
c       DO j=1,ncurv
c          index=(j-1)*nx
c          DO i=1,nx 
c           DATA(i+index)=DATA(i+index)+xnoise*CMPLX(1,1)*ran2(1)
c          ENDDO
c        ENDDO
c        enorm=0.
c        DO j=1,ncurv
c          index=(j-1)*nx
c          DO i=1,nx
c            enorm=enorm+ABS(DATA(i+index))**2
c          ENDDO
c        ENDDO
c     ENDIF
c
c---  check IF the initialization was successful:
      IF (ierrinfile.GT.0) STOP '>> IN file does not exist'
      IF (ierrcov.GT.0) STOP 'error reading IN file (cov-format)'
      IF (ierraun.GT.0) STOP 'error reading unc file (cov-format)'
      IF (ierrdat.GT.0) STOP 'error reading IN file(d-format)'
      IF (ierrD.GT.0) STOP 'error reading IN file (D-format)'
      IF (ierrhp.GT.0 ) STOP 'error reading IN file (e-format)'
      WRITE(*,*)
      WRITE(prtfil,*)
c
c---  finished ?
c
      IF (iopt(3).GE.1) CALL norm ! normalize the response 
      CALL cost(fit(1))
      WRITE(*,*)' initial fitness',fit(1)
      WRITE(prtfil,*)' initial fitness',fit(1)
c     IF (iopt(13).EQ.1) THEN
c     WRITE(*,900)(1-fit(1)),10*LOG10(1-fit(1))
c     WRITE(prtfil,900)(1-fit(1)),10*LOG10(1-fit(1))
c     900       FORMAT('  best of bartlett-power',f10.4,
c     &            ' lin',f10.4,' dB')
c     ENDIF
      IF (iforwtot.GE.  niter) THEN
         STOP ' >>>>>>>>> Only one forward model'
      ENDIF

c
c---  OPEN files
c
      IF   ((iopt(8).GE.1).OR.iopt(4).EQ.2.OR.iopt(32).GE.1 ) THEN
         CALL opfilw(13,ioer)   ! the matlab plotting file
c-----Used options to matlab file
         WRITE(13,'(a,100i3,a)')' iopt=[',(iopt(i),i=1,mopt)
         WRITE(13,'(a)') '];'   
      ENDIF
      IF ((iopt(8).GE.1).OR.iopt(4).EQ.2) THEN ! contour plot files
         CALL opfilw(28,ioer)
         CALL opfilw(29,ioer)
      ELSE
         CALL opfilw(60,ioer)
         CALL opfilw(10,ioer)
      ENDIF
      IF ((iopt(22).EQ.1)) THEN ! fipplot plot files
         CALL opfilw(19,ioer)
         WRITE(19,*)'      1024          MODU' ! This is for compatability
         CALL opfilw(20,ioer)
      ENDIF
      IF ((iopt(26).EQ.1)) THEN ! environmental file
         CALL opfilw(90,ioer)
      ENDIF
c
c---  size of parental distribution
c
      qpar=2*INT(0.5*pu*q)

c---  start the clock
      CALL cltime
      IF (iopt(4).EQ.1)THEN
c**** CALL SA
         CALL sa()
         GOTO 999
      ELSEIF (iopt(4).EQ.2)THEN
c**** CALL Gauss Newton
c         CALL gaunew(xstart,maxiter)

c powell
         DO iparm=1,nparm
            DO jparm=1,nparm
               xiun( iparm, jparm)=0.e0
            ENDDO 
            xiun( iparm, iparm)=1.e0
            xstep(iparm)=df(iparm)
         ENDDO
         tol=0.01
c         xiun(1,1)=1
c         xiun(2,1)=1 
c         xiun(1,2)=-1
c        xiun(2,2)=1 
c         xiun(1,1)=1
c         xiun(2,1)=.6/.8
c         xiun(1,2)=-.6/.8
c        xiun(2,2)=1 
c         CALL powell(xstar,xiun,xstep,
c     1        nparm,mpar,tol,iterpowell,fret)
         WRITE(*,*)'start point',(xstar(iparm),iparm=1,nparm),fit(1)
         fret=fit(1)
         itmax=20
         write(*,*)'call to powell:',xstar(1),xstar(2)
         CALL powell(xstar,xiun,xstep,fmin,
     1        nparm,mpar,M,tol,iterpowell,fret,iforwloc,itmax)
         WRITE(*,*)'retuned point',(xstar(iparm),iparm=1,nparm)
         WRITE(*,*)'retuned value',fret
         
         
         GOTO 999
      ELSEIF (iopt(4).EQ.3)THEN
c**** CALL VFSA
         WRITE(*,*)'calling VFSA module...'
         CALL vfsa()
         GOTO 999
      ELSEIF (iopt(4).EQ.4)THEN
c**** CALL MCMC
         CALL opfilw(40,ioer)
         WRITE(*,*)'calling Metropolis-Hastings module...'
         CALL mcmc
         GOTO 999
      ELSEIF (iopt(8).GE.1) THEN
c**** plot ambiguity
         WRITE(*,*)'calling contour plot'
         IF (nparm.EQ.1) THEN 
            WRITE(*,*)' at least two inversion parameters must'
            WRITE(*,*)' be specified for a contour plot'
            STOP
         ENDIF
         CALL conplot()
         GOTO 999
      ELSEIF (iopt(32).GE.1) THEN
c**** plot ambiguity
         WRITE(*,*)'calling line plot'
         CALL lineplot()
         GOTO 999
      ENDIF

c      
c---  WRITE out
c
      WRITE(prtfil,*)
c      WRITE(prtfil,*)' norm of the data         ',enorm
      WRITE(prtfil,*)' crossover probability    ',px
      WRITE(prtfil,*)' mutation  probability    ',pm
      WRITE(prtfil,*)' update    probability    ',pu
      WRITE(prtfil,*)' number of iterations     ',niter
      WRITE(prtfil,*)' size   of populations    ',q
      WRITE(prtfil,*)' number of populations    ',npop
      WRITE(prtfil,*)
c
c***  loop over number of generations
c
      iforwtotal  = 0

c*********** Number of spawning processes
      nspawn=1
c************ For DSO we USE 
c      nspawn=1
c      IF (iopt(30).EQ.7)nspawn=1
    
      
      iallmodel=0 
      iallstart= 1
      ispawn=0
      nspawn0=nspawn
      DO 100 ipop=1,npop
c      DO 100 ipop=5,npop
         WRITE(*,*)' ipop',ipop
c         WRITE(80,*)' ipop',ipop
         IF (iopt(34).EQ.1) THEN
            iseed=-getpid()*101-ipop !  "random" random seed
         ELSE
            iseed=-1001*ipop +1 ! "deterministic" random seed 
         ENDIF
         xdum=ran2(iseed)   
         itrn=0
         IF ((MOD(ipop,nspawn).NE.0) .AND. (ipop.NE.npop) ) THEN
c-- population one to nspawn-1 is spawned
            itrn=fork()      
            IF (itrn.NE.0)   THEN ! this is for the parent
               IF (ispawn.LE.nspawn) THEN
                  ispawn= ispawn+1
               ELSE
                  ispawn=1
               ENDIF
               WRITE(*,*)'ipop,itrn',ipop,itrn
               IF (itrn.GE.0) THEN
                  GOTO 99
               ENDIF
            ELSE
               ispawn=-1
               WRITE(*,*)' child process'
c     so that child processes does not WRITE to file. 
c     This causes problems for the intel compiler.
               flagpu=200  
            ENDIF
         ELSE
            IF (ipop.EQ.npop) then
               nspawn=npop-(npop/nspawn0)*nspawn0
               if (nspawn==0) nspawn=nspawn0
            endif
         ENDIF
         WRITE(*,*)'ispawn,flagpu',ispawn,flagpu,nspawn
         objbsf=1e10
         iforwpop=0
         iforwreject = 0

         IF (MOD(ipop,10).EQ.0) WRITE(prtfil,*)' Population',ipop
c     WRITE(*,*)' Population:',ipop
c---  obtain initial models and their objective FUNCTION
         CALL initga
         IF ((iopt(24).EQ.1).AND.(ipop.LE.temp0) ) THEN
            WRITE(prtfil,*)' Initial model included for 
     .           population',ipop
            CALL modeldecode(xstart,1,ierr)
            IF (ierr.NE.0) STOP ' A parameter is out of bounds'
         ENDIF
c     IF (ipop.EQ.1 .AND. iopt(26).EQ.1) THEN
c         WRITE the bartlett power out to the env file
c               IF (iopt(5).EQ.4) THEN
c                  WRITE(90,'(i4, 1000f8.1)')
c     1                 iforwpop, (frq(ifrq), ifrq=1,nfrq)
c               ELSEIF (iopt(5).EQ.5) THEN
c                  WRITE(90,'(i4, 1000f8.1)')
c     1                 iforwpop, (rdep(idep), idep=1,ndep)
c               ENDIF
c            ENDIF
c
c---    Initialize the initial population
c
         DO iq=1,q
c---      set the model and CALL the modelling routine
            CALL setmodel(iq)
c     WRITE(*,*)' calling forward model',iforwpop+1,ipop
            CALL forw2
c            WRITE(*,*)' finished forward model',iforwpop+1,ipop
            IF (lwascomp.EQ.1) THEN
               iforwtotal=iforwtotal+1
               iforwpop=iforwpop+1
               IF (iopt(3).GE.1) CALL norm ! normalize the response 
               CALL cost(fit(iq))
            ELSE
               iforwreject= iforwreject+1
               fit(iq)=10
            ENDIF 
         ENDDO
c--     sort according to fitness
         CALL sortfit(fit,fitpar)
c        STOP
c WRITE extended info for initial population.
         IF ((iopt(5).EQ.4).AND.(iopt(28).EQ.1).AND.(nfrq.GE.2)) THEN
            DO iq=1,q
               WRITE(60,'(2i5,100G13.5)')0,iq,fit(iq),
     1              (fval(model(iparm,iq),iparm),iparm=1,nparm),
     1              (fit_bart_ind(ifrq),ifrq=1,nfrq)
c     WRITE(10,'(i4,i4,g12.5,100i5)')
c     &        iq,ipop,fit(iq),(model(iparm,iq),iparm=1,nparm)
c     CALL flush(10)
               CALL flush(60)
            ENDDO
         ENDIF
         
         DO iq=1,q
            iallmodel=iallmodel+1
            allfit(iallmodel)=fit(iq)
            allmodel(nparm+1,iallmodel)=iq
            allmodel(nparm+2,iallmodel)=ipop
            DO j=1,nparm
               allmodel(j,iallmodel)= model(j,iq) ! move
            ENDDO
         ENDDO
c     ... debugging
         If (iopt(6).eq.1) then
            DO iq=1,q
               write(81,*)
     1              (fval(model(iparm,iq),iparm),iparm=1,nparm)
            enddo
         endif
c     ... debugging
c
c----  estimate run time
c 
         IF ((ipop.EQ.nspawn) .OR. 
     1        ((ipop.EQ.npop).AND.(npop.LT.nspawn))) THEN
            CALL rdtime(t1)
            WRITE(prtfil,320) t1*(niter*npop)/q
            WRITE(*,320) t1*(niter*npop)/q
 320        FORMAT('  >>> Estimated serial Inversion time: ',
     1           F12.1,' secs.')
c     WRITE(*,320) t1*(niter*NINT(npop/nspawn))/q

            
         ENDIF
c
c---  finished ?
c
         IF (iforwpop.GE.niter)
     &        STOP ' >>>>>>>>> Only one generation'
c     
c%%%%%%loop over iterations
c     
         
         DO iter=1,niter
c     WRITE(prtfil,*)'iteration',iter
            IF (iopt(6) .EQ. 1) WRITE(*,*)'iteration',iter
            IF (MOD(iter,10).EQ.0) THEN
               WRITE(*,901)iter,ipop,iforwpop      
 901           FORMAT(' Iteration:',i5,' population:',I5,
     &              ' forward model:',i9)
            ENDIF
c     
c---- obtain temperature 
c     
            temp=fit(qtemp)
            IF (temp.LE.0)  THEN
               WRITE(*,*)' temp < 0, has been corrected....'
               WRITE(*,*)' probaly the model is matched exactly'
               WRITE(*,903)qtemp,fit(qtemp)
 903           FORMAT( 'fit(',i2,')=',e12.4)
               temp=10e-3
            ENDIF
 11         CONTINUE
            expobj=0.
c***  loop for each individual
            DO iq=1,q
               expfit(iq)=EXP(-DBLE(fit(iq)/Temp) )
               expobj=expobj+expfit(iq)
            ENDDO
            IF (expobj.LE.(1.e-6)) THEN
               WRITE(*,*)'temperature to low relative to obj-function'
               WRITE(*,*)'iter,temp,expobj:',iter,temp,expobj
               temp=temp*2
               GOTO 11
            ENDIF
c     
c---  obtain probabilities  
c     
            DO iq=1,q
               pmodel(iq)=expfit(iq)/expobj
            ENDDO
c---  obtain accumulated probabilities
            pcum(1)=pmodel(1)    
            DO iq=2,q
               pcum(iq)=pcum(iq-1)+pmodel(iq)    
            ENDDO
c     
c---  generate parental distribution
c     
            DO iq=1,qpar
               DO iiq=1,q
                  IF (ran2(1).LT.pcum(iiq)) THEN
                     ipar(iq)=iiq
                     DO iparm=1,nparm
                        parents(iparm,iq)=model(iparm,iiq)
                     ENDDO
                     GOTO 101
                  ENDIF
            ENDDO
 101        CONTINUE
         ENDDO
c     
c---      gray incoding
c
         IF (iopt(33).EQ.1) THEN 
            CALL grayincode(qpar)      
         ENDIF
c
c---  carry out crossover & mutations 
c     
         CALL cross(qpar)
         CALL mutation(qpar)
         
         IF (iopt(33).EQ.1) THEN 
            CALL graydecode(qpar)      
         ENDIF
         
c     inserts the children into the model vector and 
c     perform uniqueness test
         CALL insertmodel(qpar,irep)

         DO iq=q-irep+1,q
c---  set the model and CALL the modelling routine
c     ... debugging
            If (iopt(6).eq.1) then
               write(81,*)
     1             (fval(model(iparm,iq),iparm),iparm=1,nparm)
            endif

            CALL setmodel(iq)
            CALL forw2()

c     WRITE(*,*)' finished forward model',iforwpop+1,ipop
            IF (lwascomp.EQ.1) THEN
               iforwtotal=iforwtotal+1
               iforwpop=iforwpop+1
               IF (iopt(3).GE.1) CALL norm ! normalize the response 
               CALL cost(fit(iq))
               IF ((iopt(5).EQ.4).AND.(iopt(28).EQ.1).AND.(nfrq.GE.2))
     .              THEN
                  WRITE(60,'(2i5,100G13.5)')0,iq,fit(iq),
     1                 (fval(model(iparm,iq),iparm),iparm=1,nparm),
     1                 (fit_bart_ind(ifrq),ifrq=1,nfrq)
               ENDIF
               CALL flush(60)
               

c              WRITE(10,'(i4,i4,g12.5,100i5)')
c     &        iq,ipop,fit(iq),(model(iparm,iq),iparm=1,nparm)
c              CALL flush(10)
            
               IF (iopt(15).EQ.1) THEN
c-----for optimization using hybrid scheme
                  CALL getmodelreal(xtheta)
c     gn                CALL gaunew(xtheta,5)            !maxiter)
c     gn                iforwtotal=iforwtotal+5
c     gn                 iforwpop = iforwpop+5
                  DO iparm=1,nparm
                     xstar(iparm)=xtheta(iparm)
                     xstep(iparm)=df(iparm)
                  ENDDO
                  fret=fit(iq)
                  CALL covmat(allmodel,mpar,mq_post,Nparm,
     .                 iallmodel,xiun)
                  IF (iopt(6).EQ.1) THEN
                     DO  j=1,nparm
                        WRITE(2,*)'ev',j,(xiun(j,i),i=1,nparm)
                     ENDDO
                     tol=0.01
                     WRITE(prtfil,*)'start point',
     .                    (xstar(iparm),iparm=1,nparm)
                     WRITE(prtfil,*)'start value',fit(iq)
                  ENDIF
                  WRITE(*,*)'start point',(xstar(iparm),
     .                 iparm=1,nparm),fit(iq)
c     WRITE(70,*)(xstar(iparm),iparm=1,nparm),objbsf
                  itmax=1
                  CALL powell(xstar,xiun,xstep,fmin,
     1                 nparm,mpar,M,tol,iterpowell,fret,iforwloc,itmax)
                  CALL setmodelreal(xstar)
         
                  CALL forw2
                  CALL modeldecode(xstar,iq,ierr)
                  IF (iopt(3).GE.1) CALL norm ! normalize the response 
                  CALL cost(fit(iq))
                  WRITE(*,*)'ret point',(xstar(iparm),iparm=1,nparm),
     .                 fit(iq),fret
                  iforwpop=iforwpop+iforwloc
                  iforwtot=iforwtot+iforwpop
                  IF (ierr.EQ.1) STOP 'Decoding of Powell did not work'
               ENDIF
            ELSE
               iforwreject= iforwreject+1
               fit(iq)=10
               WRITE(*,*)'reject/accept',iforwreject,iforwpop
            ENDIF 
            iallmodel=iallmodel+1
            allfit(iallmodel)=fit(iq)
            allmodel(nparm+1,iallmodel)=iq
            allmodel(nparm+2,iallmodel)=ipop
            DO j=1,nparm
               allmodel(j,iallmodel)= model(j,iq) ! move
            ENDDO
            
         ENDDO                  !iq
c--   sort according to fitness
         CALL sortfit(fit,fitpar)
  
         
         IF (fit(1).LT.objbsf) THEN
            objbsf=fit(1)
            DO i=1,nparm
               modelbsf(i)=model(i,1)
            ENDDO
            WRITE(*,902)objbsf,' at iteration',iter
            WRITE(*,*)' new best energy',objbsf,
     &           ' at iteration',iter
 902        FORMAT(' new best energy',E12.4,A,I4)
            IF (objbsf.LT.objball) THEN
c     WRITE(60,'(1x,I7,I7,G12.5,100I5)')
c     &          ipop,iforwtotal,fit(1),(model(iparm,1),iparm=1,nparm)
            ENDIF
         ENDIF
         
         IF (iforwpop.GT.niter)THEN
c---  maximum number of generation has been carried out
            GOTO 90
         ENDIF
      ENDDO                     ! iter
 90   CONTINUE

c
c---    find best of all population
c
      IF (objbsf.LT.objball) THEN
         objball=objbsf
         DO i=1,nparm
            modelball(i)=modelbsf(i)
         ENDDO
         WRITE(*,902)objbsf,'at population',ipop
c          WRITE(prtfil,902)objbsf,'at population',ipop
      ENDIF    
      WRITE(*,*)'  best energy in population',objbsf,ipop
c        WRITE(prtfil,*)'  best energy in population',objbsf,ipop

c
c---    sum up statistics
c
c        DO i=iallstart,iallmodel
c            WRITE(10,'(i4,i4,g12.5,100i5)')allmodel(nparm+1,i),
c     &          allmodel(nparm+2,i),
c     &          allfit(i),(allmodel(iparm,i),iparm=1,nparm)
c        ENDDO
c        CALL flush(10)
c
c DO a local run
c
      DO iparm=1,nparm
         xstar(iparm)=fval(modelbsf(iparm),iparm)
         xstep(iparm)=df(iparm)
      ENDDO
      fret=objbsf
c         DO iparm=1,nparm
c            DO jparm=1,nparm
c               xiun( iparm, jparm)=0.e0
c            ENDDO 
c            xiun( iparm, iparm)=1.e0
c         ENDDO
      CALL covmat(allmodel,nparm+2,mq_post1,Nparm,iallmodel,xiun,M)
      IF (iopt(6).EQ.1) THEN
         DO  j=1,nparm
            WRITE(2,*)'ev',j,(xiun(j,i),i=1,nparm)
         ENDDO
         tol=0.01
         WRITE(prtfil,*)'start point',(xstar(iparm),iparm=1,nparm)
         WRITE(prtfil,*)'start value',objbsf
         CALL flush(prtfil)
      ENDIF
      WRITE(*,*)'start point',(xstar(iparm),iparm=1,nparm)
c     WRITE(70,*)(xstar(iparm),iparm=1,nparm),objbsf
      itmax=5
      CALL powell(xstar,xiun,xstep,fmin,
     1     nparm,mpar,M,tol,iterpowell,fret,iforwloc,itmax)
c     IF (iopt(6).EQ.1) THEN
      WRITE(prtfil,*)'retuned point',(xstar(iparm),iparm=1,nparm)
      WRITE(prtfil,*)'retuned value',fret
c     ENDIF
c     WRITE(71,*)(xstar(iparm),iparm=1,nparm),fret
      IF (fret.LT. xlocobj) THEN
         xlocobj=fret
         DO iparm=1,nparm
            xlocmod(iparm)=xstar(iparm)
         ENDDO
      ENDIF
      
      
      iforwpop=iforwpop+iforwloc
      iforwtot=iforwtot+iforwpop
      WRITE(*,*)'total calls for population',ipop,iforwpop
      IF (iforwreject.GT.0) THEN
         WRITE(*,*)'rejected forward models for population'
     &        ,ipop,iforwreject
      ENDIF
      write(*,*)'gapeter ispawn,ipop,npop',ispawn,ipop,npop
      IF ((ispawn==-1).AND. (ipop.NE.npop)) THEN
         WRITE(*,*)'finnished population',ipop
         ID=getpid()
         WRITE(comfile,'(A3,I5.5)') 'dum',ID
         OPEN(UNIT=49,FILE=comfile,FORM='UNFORMATTED',STATUS='NEW')
c     OPEN(UNIT=49,FILE=comfile,FORM='FORMATTED',STATUS='NEW')         
         WRITE(*,*)' Start writing unit 49'
         WRITE(49)(xstar(iparm),iparm=1,nparm),fret
         WRITE(49) iforwpop
         DO i=iallstart,iallmodel
c             WRITE(49,'(i4,i4,g12.5,100i5)')
             WRITE(49)
     &            allmodel(nparm+1,i),
     &            allmodel(nparm+2,i),
     &            allfit(i),(allmodel(iparm,i),iparm=1,nparm)
          ENDDO
          WRITE(*,*)'...unit 49 written'
          STOP
        ELSE IF ((nspawn.GT.1).AND.(ipop.GT.1)) THEN
           DO i=1,ispawn
              retval=wait(status)
              WRITE(*,*)'retval,status',retval,status
              IF (retval.GE.0) THEN 
                 WRITE(comfile,'(A3,I5.5)') 'dum',retval
                 OPEN(UNIT=49,FILE=comfile,FORM='UNFORMATTED')         
c     OPEN(UNIT=49,FILE=comfile,FORM='FORMATTED')         
                 READ(49)(xstar(iparm),iparm=1,nparm),fret
                 IF (fret.LT. xlocobj) THEN
                    xlocobj=fret
                    DO iparm=1,nparm
                       xlocmod(iparm)=xstar(iparm)
                    ENDDO
                 ENDIF
                 READ(49) iforwpop
                 iforwtot=iforwtot+iforwpop
                 DO iq=1,mq_post
                    iallmodel=iallmodel+1
                    
c     READ(49,'(i4,i4,g12.5,100i5)',END=9001)
                    READ(49,END=9001)
     &                   allmodel(nparm+1,iallmodel),
     &                   allmodel(nparm+2,iallmodel),
     &                   allfit(iallmodel),
     &                   (allmodel(iparm,iallmodel),iparm=1,nparm)
        IF (allfit(iallmodel).LT.objball) THEN
          objball=allfit(iallmodel)
          DO iparm=1,nparm
            modelball(iparm)=allmodel(iparm,iallmodel)
          ENDDO
c          WRITE(*,902)objbsf,'at population',ipop
c          WRITE(prtfil,902)objbsf,'at population',ipop
        ENDIF    
                 ENDDO
 9001            CONTINUE
                 CLOSE(49)
                 iallmodel=iallmodel-1
                 IF (iopt(6).EQ.0) THEN
                    CALL system('/bin/rm -f '//comfile//' &')
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
        WRITE(*,*)'writing iallstart iallmodel=',
     &       iallstart,iallmodel
        
        DO iq=iallstart,iallmodel
c     WRITE(10,'(i5,i5,g14.7,50i5)')
           WRITE(10,*)
     &          allmodel(nparm+1,iq),allmodel(nparm+2,iq),
     &          allfit(iq),(allmodel(iparm,iq),iparm=1,nparm)
        ENDDO
        WRITE(*,*)' ... writing done'
        
        CALL flush(10)
        iallstart=iallmodel+1
        IF (LOG10(1.0*iallmodel).GT.nsearch) THEN
           WRITE(*,*)'All posible combinations have been explored'
           GOTO 102
        ENDIF
 99     CONTINUE                ! ipop  for spawned process
 100  CONTINUE                  ! ipop  number of populations 
 102  CONTINUE                  ! outside all loops
c     
c---  WRITE unit 60 
c
      IF ((iopt(5).NE.4).OR.(iopt(28).NE.1).OR.(nfrq.EQ.1)) THEN
         DO iq=1,iallmodel
c     WRITE(60,'(2i5,100G13.5)')0,iq,allfit(iq),
c     1           (fval(allmodel(iparm,iq),iparm),iparm=1,nparm)
c     1           ,(allfit_bart_ind(ifrq,iq),ifrq=1,nfrq)
            WRITE(60,'(2i7,100G13.5)')0,iq,allfit(iq),
     1           (fval(allmodel(iparm,iq),iparm),iparm=1,nparm)
         ENDDO
      ENDIF
c
c--------- 
c
      CALL cost_max(fit_teo)
      WRITE(lun,*)' best obtainable theoretical fit', fit_teo
c
c---- WRITE out the RESULT
c
      DO 111 ihelp=1,2
         IF (ihelp.EQ.1) lun=6
         IF (ihelp.EQ.2) lun=prtfil
         WRITE(lun,*)
         
         WRITE(lun,*)'  best of all energy',objball
         IF (iopt(13).EQ.1) THEN
            IF (iopt(18).EQ.0) THEN
c     WRITE(lun,900)(1-objball),10*LOG10(1-objball)
c     ELSEIF IF (iopt(20).EQ.0) THEN
c     WRITE(lun,*)objball,objball
            ENDIF
         ENDIF
         WRITE(lun,*)' best-of all      deviation  '       
         DO i=1,nparm
            WRITE(lun,'((i4,f10.3,7x),3f10.3)')i,
     &           fval(modelball(i),i),fval(modelball(i),i)-xstart(i)
         ENDDO

         WRITE(lun,*)'local fit',xlocobj
         DO i=1,nparm
            WRITE(lun,'((i4,f10.3,7x),3f10.3)')i,
     &           xlocmod(i),xlocmod(i)-xstart(i)
         ENDDO
         
         
         WRITE(lun,'(/a,i8)')' Forward  modelling calls:', iforwtot
 111  CONTINUE
      
      
      WRITE(*,*)' writing to results.m !!!'
c
c---  WRITE the best RESULT to RESULTS.M
cf90: position='append', f77: access='append',
      OPEN(unit=8,file='results.m',position='append',
     &     status='unknown')
      WRITE(8,*)' %%%%%%%%%%%%%%% '
      WRITE(8,*)' fit=[ fit; ',objball,'];'
      WRITE(8,*)' res=[ res; ['
      DO iparm=1,nparm      
         WRITE(8,'(100f15.3)')
     &        fval(modelball(iparm),iparm)
      ENDDO 
      WRITE(8,*)' ]];'
      WRITE(8,*)' pointer=[ pointer; ['
      DO iparm=1,nparm      
         WRITE(8,'(3i6)')par2phy(iparm),par2lay(iparm),par3(iparm)
      ENDDO 
      WRITE(8,*)' ]];'
      WRITE(8,*)' fitloc=[ fitloc; ',xlocobj,'];'
      WRITE(8,*)' resloc=[ resloc; ['
      DO         iparm=1,nparm      
         WRITE(8,'(100f15.3)')
     &        xlocmod(iparm)
      ENDDO 
      WRITE(8,*)' ]];'
      CLOSE(8)
      WRITE(*,*)'...results.m is written !'

c--- write best results from powell to snap_powell.dat
      iWriteTrf=1      
      CALL setmodelreal(xlocmod)
      CALL forw2
      write(8,*)'Finished best forward model'
      if (iopt(30).eq.1) then
         call system('/bin/mv -f snap.dat snap_powell.dat')
      endif
 999  CALL rdtime(t1)
      WRITE(prtfil,310) t1
      WRITE(*,310) t1
 310  FORMAT(' Inversion, time: ',F12.3,' secs.')
      
c     CALL lib$show_timer()
      END
