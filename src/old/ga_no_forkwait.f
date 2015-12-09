      PROGRAM ga
c****************************************************************************
c   PROGRAM for inversion of seismo-acoustic DATA using genetic algorithms.
c   The RESULT is presented as the most likely model, the mean model and the 
c   standard variation.
c
c   PETER GERSTOFT, 1992 
c****************************************************************************  

                          ! Common block definitions are
      INCLUDE 'comopt.h'  ! INCLUDEd in these files.
      INCLUDE 'comforw.h'
      REAL fval              ! function for computation of real values
C**   local variables
      INTEGER iter,iq,iiq,i,J,  ! iterative indices.
     >     ifreq                ! frequency index.
      INTEGER ipar(mq)
      INTEGER ihelp,ipop,lun,index
      INTEGER modelball(mpar)   ! best of all model
      INTEGER modelbsf(mpar)    ! best so far model
      INTEGER iforwtot,iforwpop,iforwtotal ! counters forward-models
      REAL    objbsf /1e10/,objball /1e10/
      REAL t1,ahelp
      REAL fit(mq),fitpar(mq),expfit(mq),
     >     pmodel(mq),pcum(mq),
     >     expobj,temp
      REAL Xstart(mpar)
      REAL ran2,xtheta(mpar)
      REAL nsearch,xdum
      INTEGER ix,itrn, fork,ispawn,wait ,retval,status,nspawn,ifrq,idep
      EXTERNAL fork,wait
      INTEGER qpar,irep,qtemp/1/,iparm,maxiter/100/,
     >     ierrdat/0/,ierrcov/0/,ierrhp/0/, !error flags for read_cov and read_hp
     >     ierrinp/0/,ioer/0/,ierr/0/,ierrD/0/
c      INTEGER iter,iforwpop
      COMMON/iterpar/iter,iforwpop 
      WRITE(*,*)'Running SAGA version 3.1, August 1999' 

      iseed=6511

      xdum=ran2(iseed)          ! initialize random generator
c
c---- initialization
c
      DO i=1,mobs
         weight(i)=1.0
      ENDDO
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

      nsearch=0
      DO iparm=1,nparm
         nsearch=nsearch+LOG10(1.0*ndigit(Iparm))
      ENDDO
      WRITE(*,'(a,f6.2)')' Size of the search space: 10**',nsearch
      WRITE(prtfil,'(a,f6.2)')' Size of the search space: 10**',nsearch

      IF (itrans(4).GE.1) THEN
         CALL opfilr(3,ioer)
         CALL Read_weight
         CLOSE(3)
      ENDIF

C     Reading DATA from COV or IN file (unit 3 or 2). This can be in
C     covariance matrix form (option C), or it can be in hydrophone
C     spectra vectors (option d) or ? (D) or Saplot FORMAT (e).
      IF (iopt(2).GE.1) THEN    ! Options c, d, D, or e.
         CALL opfilr(2,ierrinp)
         IF (ierrinp.NE.0) GOTO 21 
         
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
            
         ENDIF
         CLOSE(2)
 21   END IF

C     Run the forward model using the parameters specified in the run
C     (.DAT) file.  The results of this run can later be written to the
C     OBS file in either covariance matrix or hydrophone DATA FORMAT.
      CALL forwinit

c
c---- should DATA be weighted
c
      IF (iopt(3).GE.1) THEN
c     WRITE(*,*)' the  is being weighted...'
         CALL norm
      ENDIF   

      iforwtot=1

      CALL getmodelreal(xstart)
      CALL modeldecode(xstart,1,ierr)
      DO iparm=1,nparm
         xstar(iparm)=xstart(iparm)
      ENDDO

C     Option W: Writes forward model output to the OBS file. DATA can
C     be in Covariance Matrix (option c), or h/p DATA vector (d) FORMAT.
      IF (iopt(21).EQ.1) THEN ! option W.
         CALL opfilw(30,ioer)
         IF (iopt(13).EQ.1) THEN ! Option c (Covariance matrix)
            WRITE(30,'(a1,a31)') '!',' File format: Covariance matrix' 
            CALL header(30)
            WRITE(*,*) ' The generated cov-matrix will be',
     >           ' written to the OBS file'
            DO ix = 1, nx
               DO ifreq = 1, nfrq
                  CALL WRITE_COV(ifreq,ix)
               ENDDO
            ENDDO
c--- reread covariance matrix
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
            ELSEIF (iopt(2).EQ.1) THEN
               WRITE(30,'(a1,a)')'!',title
               CALL header(30)
               DO iparm=1,nparm
                  WRITE(30,'(a1,g)')'!',xstart(iparm)
               ENDDO
               WRITE(*,*)'ndep,nx',ndep,nx,iopt(2)
               DO j=1,ncurv
                  index=(j-1)*nx
                  DO i=1,nx
                     WRITE(30,*)j,i,resp(i+index),index+i
                  ENDDO
               ENDDO
            ELSEIF (iopt(2).EQ.2) THEN
c     DO j=1,ndep
               index=(j-1)*nx
               DO i=1,nx
                  WRITE(30,'(100G13.6)')
     &                 xranges(i),(ABS(resp(i+(j-1)*nx)),j=1,ncurv)
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

c      IF (iopt(2).EQ.0.AND.iopt(13).EQ.0) THEN
c---    move DATA from resp to DATA  for synthetic inversions
c        DO j=1,ncurv
c          index=(j-1)*nx
c          DO i=1,nx
c            DATA(i+index)=resp(i+index)
cc          ENDDO
c       ENDDO
c      ENDIF
c
c---- should DATA be weighted
c
      IF (iopt(3).GE.1) THEN
         WRITE(*,*)' Data is being weighted...'
         CALL normdata
      ENDIF   
c--- for response
      IF (iopt(25).EQ.1) THEN
         WRITE(*,*)' transforming data; Only magnitude is used'
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
      IF (ierrinp.GT.0) STOP '>> IN file does not exist'
      IF (ierrcov.GT.0) STOP 'error reading IN file'
      IF (ierrdat.GT.0) STOP 'error reading IN file'
      IF (ierrhp.GT.0 ) STOP 'error reading IN file'
      WRITE(*,*)
      WRITE(prtfil,*)
c
c---  finished ?
c
      IF (iforwtot.GE.  niter) THEN
         IF (iopt(3).GE.1) CALL norm ! normalize the response 
         CALL cost(fit(1))
         WRITE(*,*)' initial fitness',fit(1)
         WRITE(prtfil,*)' initial fitness',fit(1)
         IF (iopt(13).EQ.1) THEN
            WRITE(*,900)(1-fit(1)),10*LOG10(1-fit(1))
            WRITE(prtfil,900)(1-fit(1)),10*LOG10(1-fit(1))
 900        FORMAT('  best of all bartlett- power',f10.4,
     &           ' lin',f10.4,' dB')
         ENDIF
         STOP ' >>>>>>>>> Only one forward model'
      ENDIF

c
c---  OPEN files
c
      IF ((iopt(8).GE.1).OR.iopt(4).EQ.2) THEN ! contour plot files
         CALL opfilw(28,ioer)
         CALL opfilw(29,ioer)
         CALL opfilw(13,ioer)   ! the matlab plotting file
c-----Used options to matlab file
         WRITE(13,'(a,30i3,a)')' iopt=[',(iopt(i),i=1,30),'];'   
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
c      OPEN(60,STATUS='NEW',CARRIAGECONTROL='LIST',recl=32000)
c      OPEN(10,STATUS='NEW',CARRIAGECONTROL='LIST',recl=32000)

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
      ELSEIF (iopt(4).EQ.3)THEN
c**** CALL VFSA
         CALL vfsa()
         GOTO 999
      ELSEIF (iopt(4).EQ.2)THEN
c**** CALL Gauss Newton
         CALL gaunew(xstart,maxiter)
         GOTO 999
      ELSEIF (iopt(8).GE.1) THEN
c**** plot ambiguity
         WRITE(*,*)'calling contour plot'
         CALL conplot()
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
      iforwtotal=0
c*********** Number of spawning processes
c      nspawn=2
c************ For DSO we USE 
      nspawn=1
c      IF (iopt(30).EQ.7)nspawn=1
    

      DO 100 ipop=1,npop
c       DO 100 ipop=8,npop
         WRITE(*,*)' ipop',ipop
         iseed=-1001*ipop +1
         xdum=ran2(iseed)       ! 
         itrn=0
         ispawn=0
         IF ((MOD(ipop,nspawn).NE.0) .AND. (ipop.NE.npop) ) THEN
c-- population one to nspawn-1 is spawned
c            itrn=fork()      
            IF (itrn.NE.0)   THEN ! this is for the parent
               WRITE(*,*)'ipop,itrn',ipop,itrn
               IF (itrn.GE.0) THEN
                  GOTO 99
               ENDIF
            ELSE
               ispawn=1
            ENDIF
         ENDIF
         objbsf=1e10
         iforwpop=0

         IF (MOD(ipop,10).EQ.0) WRITE(prtfil,*)' Population',ipop
         WRITE(*,*)' Population:',ipop
c---    obtain initial models and their objective FUNCTION
         CALL initga
         IF ( (iopt(24).EQ.1).AND.(ipop.LE.temp0) ) THEN
            WRITE(prtfil,*)' Initial model included for population',
     1           ipop
            CALL modeldecode(xstart,1,ierr)
            IF (ierr.NE.0) STOP ' A parameter is out of bounds'
         ENDIF
         IF (ipop.EQ.1 .AND. iopt(26).EQ.1) THEN
c     WRITE the bartlett power out to the env file
            IF (iopt(5).EQ.4) THEN
               WRITE(90,'(i4, 1000f8.1)')
     1              iforwpop, (frq(ifrq), ifrq=1,nfrq)
            ELSEIF (iopt(5).EQ.5) THEN
               WRITE(90,'(i4, 1000f8.1)')
     1              iforwpop, (rdep(idep), idep=1,ndep)
            ENDIF
         ENDIF
c     
c---    Initialize the initial population
c     
         DO iq=1,q
c---      set the model and CALL the modelling routine
            CALL setmodel(iq)
c            WRITE(*,*)' calling forward model',iforwpop+1,ipop
            CALL forw2
c            WRITE(*,*)' finished forward model',iforwpop+1,ipop
c          WRITE(*,*)' lwascomp',lwascomp
            IF (lwascomp.EQ.1) THEN
               iforwtotal=iforwtotal+1
               iforwpop=iforwpop+1
               IF (iopt(3).GE.1) CALL norm ! normalize the response 
               CALL cost(fit(iq))
               IF (ipop.EQ.1 .AND. iopt(26).EQ.1) THEN
c     WRITE the bartlett power out to the env file
                  IF (iopt(5).EQ.4) THEN
                     WRITE(90,'(i4, 1000f8.5)')
     1                    iforwpop, (fit_bart_ind(ifrq), ifrq=1,nfrq)
                  ELSEIF (iopt(5).EQ.5) THEN
                     WRITE(90,'(i4, 1000f8.5)')
     1                    iforwpop, (fit_bart_ind(idep), idep=1,ndep)
                  ENDIF
               ENDIF
            ELSE
c     GOTO 20
               fit(iq)=10
            ENDIF 
         ENDDO
c--   sort according to fitness
         CALL sortfit(fit,fitpar)
         
c     DO iq=1,q
c     WRITE(99,*)0,iq,fit(iq),
c     1                     (fval(model(iparm,iq),iparm),iparm=1,2)
c     ENDDO
         
c     
c----  estimate run time
c 
         IF ((ipop.EQ.nspawn) .OR. 
     1        ((ipop.EQ.npop).AND.(npop.LT.nspawn))) THEN
            CALL rdtime(t1)
            WRITE(prtfil,320) t1*(niter*npop)/q
            WRITE(*,320) t1*(niter*npop)/q
 320        FORMAT(/1H ,' >>> Estimated serial Inversion time: ',
     1           F12.1,' secs.')
c     WRITE(*,320) t1*(niter*NINT(npop/nspawn))/q


         ENDIF
c
c---  finished ?
c
         IF (iforwpop.GE.niter)
     &        STOP ' >>>>>>>>> Only one generation'
c
c%%%%%%    loop over iterations
c

         DO iter=1,niter
            WRITE(prtfil,*)'iteration',iter
            IF (iopt(6).EQ.1) WRITE(*,*)'iteration',iter
            IF (MOD(iter,10).EQ.0) THEN
               WRITE(*,901)iter,ipop,iforwpop      
 901           FORMAT(' Iteration:',i5,' population:',I5,
     &              ' forward model:',i9)
            ENDIF
c     DO iq=1,q
c     WRITE(60,'(1x,I6,I4,G12.5,100I5)')
c     &          iforwtotal,iq,fit(iq),(model(iparm,iq),iparm=1,nparm)
c     ENDDO
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
c***      loop for each individual
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
c---      obtain probabilities  
c
            DO iq=1,q
               pmodel(iq)=expfit(iq)/expobj
            ENDDO
c---      obtain accumulated probabilities
            pcum(1)=pmodel(1)    
            DO iq=2,q
               pcum(iq)=pcum(iq-1)+pmodel(iq)    
            ENDDO
c
c---      generate parental distribution
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
 101           CONTINUE
            ENDDO
c
c---      carry out crossover & mutations 
c
            CALL cross(qpar)
            CALL mutation(qpar)
c         inserts the children into the model vector and 
c         perform uniqueness test
            CALL insertmodel(qpar,irep)

            DO iq=q-irep+1,q
c---  set the model and CALL the modelling routine
               CALL setmodel(iq)
               CALL forw2()
c     WRITE(*,*)' finished forward model',iforwpop+1,ipop
               IF (lwascomp.EQ.1) THEN
                  iforwtotal=iforwtotal+1
                  iforwpop=iforwpop+1
                  IF (iopt(3).GE.1) CALL norm ! normalize the response 
                  CALL cost(fit(iq))
                  
                  IF (iopt(15).EQ.1) THEN
c-----          for optimization using hybrid scheme
                     CALL getmodelreal(xtheta)
                     CALL gaunew(xtheta(1),5) !maxiter)
                     iforwtotal=iforwtotal+5
                     iforwpop = iforwpop+5
                     CALL modeldecode(xtheta,iq,ierr)
                     IF (iopt(3).GE.1) CALL norm ! normalize the response 
                     CALL cost(fit(iq))
                  ENDIF
               ELSE
                  fit(iq)=10
               ENDIF 
            ENDDO               !iq
c--       sort according to fitness
            CALL sortfit(fit,fitpar)
  
c      DO iq=1,q
c       WRITE(99,*)iter,iq,fit(iq),
c     1                     (fval(model(iparm,iq),iparm),iparm=1,2)
c      ENDDO
c
c---      find bsf population
c     
            WRITE(60,'(1x,I7,I7,G12.5,100I5)')
     &           ipop,iforwpop,fit(1),(model(iparm,1),iparm=1,nparm)
            CALL flush(60)
            IF (fit(1).LT.objbsf) THEN
               objbsf=fit(1)
               DO i=1,nparm
                  modelbsf(i)=model(i,1)
               ENDDO
               WRITE(*,902)objbsf,' at iteration',iter
               WRITE(prtfil,902)objbsf,' at iteration',iter
 902           FORMAT(' new best energy',E12.4,A,I4)
               IF (objbsf.LT.objball) THEN
c               WRITE(60,'(1x,I7,I7,G12.5,100I5)')
c     &          ipop,iforwtotal,fit(1),(model(iparm,1),iparm=1,nparm)
               ENDIF
            ENDIF
            
            IF (iforwpop.GT.niter)THEN
c---        maximum number of generation has been carried out
               GOTO 90
            ENDIF
         ENDDO                  ! iter
 90      CONTINUE

c
c---    find best of all population
c
         IF (objbsf.LT.objball) THEN
            objball=objbsf
            DO i=1,nparm
               modelball(i)=modelbsf(i)
            ENDDO
            WRITE(*,902)objbsf,'at population',ipop
            WRITE(prtfil,902)objbsf,'at population',ipop
         ENDIF    
         WRITE(*,*)'  best energy in population',objbsf,ipop
         WRITE(prtfil,*)'  best energy in population',objbsf,ipop

c
c---    sum up statistics
c
         DO i=1,q
            WRITE(10,'(i4,i4,g12.5,100i5)')
     &           i,ipop,fit(i),(model(iparm,i),iparm=1,nparm)
         ENDDO
         CALL flush(10)
         iforwtot=iforwtot+iforwpop
         WRITE(*,*)'total calls for population',ipop,iforwpop

         IF (ispawn.EQ.1) THEN
            WRITE(*,*)'finnished population',ipop
            STOP
         ELSE IF ((nspawn.GT.1).AND.(ipop.GE.nspawn)) THEN
            DO i=1,nspawn
c            retval=wait(status)
c            WRITE(*,*)'retval,staus',retval,status
            ENDDO
         ENDIF
 99      CONTINUE               ! ipop  number of populations 
 100  CONTINUE                  ! ipop  number of populations 

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

         WRITE(lun,'(/a,i8)')' Forward  modelling calls:', iforwtot
 111  CONTINUE

 999  CALL rdtime(t1)
      WRITE(prtfil,310) t1
      WRITE(*,310) t1
 310  FORMAT(/1H ,' Inversion, time: ',F12.3,' secs.')

c     CALL lib$show_timer()
      END
