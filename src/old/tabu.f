      PROGRAM ga
c****************************************************************************
c   program for inversion of seismo-acoustic data using genetic algorithms.
c   The result is presented as the most likely model, the mean model and the 
c   standard variation.
c
c   PETER GERSTOFT, 1992 
c****************************************************************************  

                          ! Common block definitions are
      INCLUDE 'comopt.h'  ! INCLUDEd in these files.
      INCLUDE 'comforw.h'
C**   local variables
      INTEGER iter,iq,iiq,i,J,  ! iterative indices.
     >   ifreq    ! frequency index.
      INTEGER ipar(mq)
      INTEGER ihelp,ipop,lun,index
      INTEGER modelball(mpar) ! best of all model
      INTEGER modelbsf(mpar)  ! best so far model
      LOGICAL lcomp(mq)
      INTEGER iforwtot,iforwpop,iforwtotal  ! counters forward-models
      REAL    objbsf /1e10/,objball /1e10/
      REAL t1,ahelp
      REAL fit(mq),fitpar(mq),expfit(mq),
     >     pmodel(mq),pcum(mq),
     >     expobj,temp
      REAL Xstart(mpar)
      REAL ran2,xtheta(mpar)
      REAL xnoise,nsearch,xdum
      integer itrn, fork,ispawn,wait ,retval,status,nspawn
      external fork,wait
      INTEGER qpar,irep,qtemp/1/,iparm,maxiter/100/,
     >    idep,                           ! depth index.
     >    ierrdat/0/,ierrcov/0/,ierrhp/0/,!error flags for read_cov and read_hp.
     >    ierrinp/0/,ioer/0/,ierr/0/

c---  for tabu serach
      integer tabu(mpar,2),modelcur(Mpar),ntabu/6/,delta,iparm_delt
      real fitcur
      integer ii
      integer i1,i2,i3,i4,idum
      iseed=6511
      xdum=ran2(iseed)   ! initialize random generator
c
c---- initialization
c
      do i=1,mobs
         weight(i)=1.0
      enddo
c      OPEN(6, STATUS='NEW',CARRIAGECONTROL='LIST')
c
      call opfilr(1,ioer)
      call opfilw(7,ioer)

c
c---  Set some forward model specific flags
c
      CALL forwardmodel(iopt,mopt)
c
c---  read input data
c
      CALL input 

      nsearch=0
      DO iparm=1,nparm
        nsearch=nsearch+log10(1.0*ndigit(Iparm))
      ENDDO
      WRITE(*,'(a,f6.2)')' Size of the search space: 10**',nsearch
      WRITE(prtfil,'(a,f6.2)')' Size of the search space: 10**',nsearch

      if (iopt(3).eq.4) then
        call opfilr(3,ioer)
        call Read_weight
        close(3)
      endif

C     Reading data from COV or IN file (unit 3 or 2). This can be in
C     covariance matrix form (option C), or it can be in hydrophone
C     spectra vectors (option d) or ? (D) or Saplot format (e).
      IF (iopt(2).ge.1) THEN ! Options c, d, D, or e.
        call opfilr(2,ierrinp)
        if (ierrinp.ne.0) goto 21 

        IF (iopt(13).eq.1) THEN     ! Option c: Read covariance matrix
          WRITE(*,*)' Reading covariance matrix' ! from COV file.
          CALL READ_COV(ierrcov)

        ELSEIF (iopt(2).eq.1) THEN ! Option d: Read h/p data vectors
          CALL readdat2(ierrdat)             ! file in ? format.

        ELSEIF (iopt(2).eq.2) THEN ! Option D: Read data from IN
          CALL readdata             ! file in ? format.

        ELSEIF (iopt(2).eq.3) THEN ! Option e: Read data from IN
          WRITE(*,*)' Reading hydrophone data' ! from IN file.
          CALL READ_HP(ierrhp)
c         CALL readdatiso           ! file in 'rev' format.
        
        ENDIF
        close(2)
21    END IF

C     Run the forward model using the parameters specified in the run
C     (.DAT) file.  The results of this run can later be written to the
C     OBS file in either covariance matrix or hydrophone data format.
      CALL forwinit
      iforwtot=1
      CALL getmodelreal(xstart)
      CALL modeldecode(xstart,1,ierr)
      DO iparm=1,nparm
         xstar(iparm)=xstart(iparm)
      ENDDO

C     Option W: Writes forward model output to the OBS file. Data can
C     be in Covariance Matrix (option c), or h/p data vector (d) format.
      IF (iopt(21).eq.1) THEN ! option W.
        call opfilw(30,ioer)
        IF (iopt(13).eq.1) THEN ! Option c (Covariance matrix)
          WRITE(30,'(a1,a31)') '!',' File format: Covariance matrix' 
          WRITE(*,*)
     >    ' The generated cov-matrix will be written to the OBS file'
          WRITE(*,*)
     >    ' and it will be used as the baseline model in the inversion' 
          DO ifreq = 1, nfrq
            CALL WRITE_COV(ifreq)
          ENDDO

         ELSE
          IF (iopt(2).eq.3) THEN ! Option e (hydrophone data vectors)
            WRITE(30,'(a1,a32)') '!',' File format: Hydrophone vectors' 
            WRITE(*,*)
     >    ' The generated h/p data will be written to the OBS file'       
            DO ifreq = 1, nfrq
               CALL WRITE_HP(ifreq)
            ENDDO
          ELSEIF (iopt(2).eq.1) THEN
            WRITE(30,'(a1,a)')'!',title
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
          ELSEIF (iopt(2).eq.2) THEN
            DO j=1,ndep
              index=(j-1)*nx
              DO i=1,nx
                WRITE(30,*)xranges(i),abs(resp(i+index))
  	      ENDDO
            ENDDO

          ELSE
            WRITE(*,*)'Error: Output format for option W not defined.'
            STOP
          ENDIF
        ENDIF
	CLOSE(30)
      ENDIF

      IF (iopt(2).eq.0.and.iopt(13).eq.0) THEN
c---    move data from resp to data  for synthetic inversions
        DO j=1,ncurv
          index=(j-1)*nx
          DO i=1,nx
            data(i+index)=resp(i+index)
          ENDDO
        ENDDO
      ENDIF
c
c---- should data be weighted
c
      IF (iopt(3).ge.1) THEN
         write(*,*)' Data is being weighted...'
         CALL normdata
      ENDIF   
c--- for response
      if (iopt(25).eq.1) then
       write(*,*)' transforming data; Only magnitude is used'
        do j=1,ncurv
          index=(j-1)*nx
          do i=1,nx
c           WRITE(prtfil,*)' data',data(j,idep),abs(data(j,idep))
           data(i+index)=abs(data(i+index))
         enddo
        enddo
      endif
      if (iopt(27).eq.1) then
        write(*,*)' transforming data: to unit norm'
        ahelp=0
        do j=1,ncurv
          index=(j-1)*nx
          do i=1,nx
           ahelp=ahelp+abs(data(i+index))**2
          enddo
        enddo
        ahelp=1.0/sqrt(ahelp)
        do j=1,ncurv
          index=(j-1)*nx
          do i=1,nx
           data(i+index)=(data(i+index))*ahelp
         enddo
        enddo
      endif
cc
c---  energy normalization
c   
      enorm=0.
      DO j=1,ncurv
        index=(j-1)*nx
        DO i=1,nx
          enorm=enorm+abs(data(i+index))**2
c           WRITE(*,*)'from ga data',j,i,data(i+index)
        ENDDO
      ENDDO
c
c---  should noise be added ??? 
c
      IF (iopt(11).eq.1) THEN
        WRITE(*,*)' ***** warning noise is being added ********'  
        xnoise=sqrt(enorm)/5./sqrt(2.)
        xnoise=0.
        DO j=1,ncurv
          index=(j-1)*nx
          DO i=1,nx 
           data(i+index)=data(i+index)+xnoise*cmplx(1,1)*ran2(1)
          ENDDO
        ENDDO
        enorm=0.
        DO j=1,ncurv
          index=(j-1)*nx
          DO i=1,nx
            enorm=enorm+abs(data(i+index))**2
          ENDDO
        ENDDO
      ENDIF
c
c---  check if the initialization was successful:
      IF (ierrinp.gt.0) STOP '>> IN file does not exist'
      IF (ierrcov.gt.0) STOP 'error reading IN file'
      IF (ierrdat.gt.0) STOP 'error reading IN file'
      IF (ierrhp.gt.0 ) STOP 'error reading IN file'
c
c---  finished ?
c
      IF (iforwtot.ge.  niter) THEN
            IF (iopt(3).ge.1) CALL norm   ! normalize the response 
            CALL cost(fit(1))
        WRITE(*,*)' initial fitness',fit(1)
        WRITE(prtfil,*)' initial fitness',fit(1)
        IF (iopt(13).eq.1) THEN
           WRITE(*,900)(1-fit(1)),10*log10(1-fit(1))
           WRITE(prtfil,900)(1-fit(1)),10*log10(1-fit(1))
 900       FORMAT('  best of all bartlett- power',f10.4,
     &            ' lin',f10.4,' dB')
        ENDIF
        STOP ' >>>>>>>>> Only one forward model'
      ENDIF

c
c---  open files
c
      if ((iopt(8).eq.1)) then          ! contour plot files
        call opfilw(28,ioer)
        call opfilw(29,ioer)
      else
        call opfilw(60,ioer)
        call opfilw(10,ioer)
      endif
      if ((iopt(22).eq.1)) then          ! fipplot plot files
        call opfilw(19,ioer)
        call opfilw(20,ioer)
      endif
c      OPEN(60,STATUS='NEW',CARRIAGECONTROL='LIST',recl=32000)
c      OPEN(10,STATUS='NEW',CARRIAGECONTROL='LIST',recl=32000)

c
c---  size of parental distribution
c
      qpar=2*int(0.5*pu*q)

c---  start the clock
      CALL cltime

      IF (iopt(4).eq.1)THEN
c**** call SA
         CALL sa()
         goto 999
      ELSEIF (iopt(4).eq.2)THEN
c**** call Gauss Newton
         CALL gaunew(xstart,maxiter)
         goto 999
      ELSEIF (iopt(8).eq.1) THEN
c**** plot ambiguity
        WRITE(*,*)'calling contour plot'
        CALL conplot()
        goto 999
      ENDIF
c
c This is for grid search
c
      idum=1
      if (idum.eq.1) then
        iforwtotal=0
        do i1=0,ndigit(1)-1
          call setmodelx(par2lay(1),par2phy(1),fval(i1,1))
          do i2=0,ndigit(2)-1
            call setmodelx(par2lay(2),par2phy(2),fval(i2,2))
            if (nparm.eq.2) then
              call forw2
              call norm
              call cost(fit(1))
              iforwtotal=iforwtotal+1
              WRITE(10,'(1x,I8,I4,G15.8,100I5)')
     &          iforwtotal,1,fit(1),i1,i2
            else if (nparm.ge.3) then
              do i3=0,ndigit(3)-1
                call setmodelx(par2lay(3),par2phy(3),fval(i3,3))
                if (nparm.eq.3) then
                 call forw2
                 call norm
                 call cost(fit(1))
                 iforwtotal=iforwtotal+1
                 WRITE(10,'(1x,I6,I4,G15.8,100I5)')
     &           iforwtotal,1,fit(1),i1,i2,i3
                else if (nparm.eq.4) then
                 do i4=0,ndigit(4)-1
                  call setmodelx(par2lay(4),par2phy(4),fval(i4,4))
                  call forw2
                  call norm
                  call cost(fit(1))
                  iforwtotal=iforwtotal+1
                  WRITE(10,'(1x,I6,I4,G15.8,100I5)')
     &              iforwtotal,1,fit(1),i1,i2,i3,i4
                 enddo
                endif
               enddo
            endif
          enddo
        enddo
        stop
      endif

c      
c---  write out
c
      WRITE(prtfil,*)
      WRITE(prtfil,*)' norm of the data         ',enorm
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
      nspawn=1

      DO 100 ipop=1,npop
          write(*,*)' ipop',ipop
          iseed=-1001*ipop +1
          xdum=ran2(iseed)   ! 
        itrn=0
        ispawn=0
        if ((mod(ipop,nspawn).ne.0) .and. (ipop.ne.npop) ) then
           itrn=fork()      
            if (itrn.ne.0)   then  ! this is for the parent
               write(*,*)'ipop,itrn',ipop,itrn
              if (itrn.ge.0) then
                goto 99
              endif
            else
               ispawn=1
            endif
        endif
        objbsf=1e10
        iforwpop=0

c---    reset tabu list
        do iparm=1,nparm
           tabu(iparm,1)=0
           tabu(iparm,2)=0
        enddo
        q=1

        IF (mod(ipop,10).eq.0) WRITE(prtfil,*)' Population',ipop
        WRITE(*,*)' Population:',ipop
c---    obtain initial models and their objective function
        CALL initga
        IF ( (iopt(24).eq.1).and.(ipop.le.temp0) ) THEN
          WRITE(prtfil,*)' Initial model included for population',ipop
          CALL modeldecode(xstart,1,ierr)
          if (ierr.ne.0) stop ' A parameter is out of bounds'
        ENDIF
c
c---    Initialize the initial population
c
 	DO iq=1,q
c---      set the model and call the modelling routine
          CALL setmodel(iq)
          CALL forw2
          IF (lwascomp.eq.1) THEN
            iforwtotal=iforwtotal+1
            iforwpop=iforwpop+1
            IF (iopt(3).ge.1) CALL norm   ! normalize the response 
            CALL cost(fit(iq))
          ELSE
            fit(iq)=10
          ENDIF 
        ENDDO
c--     sort according to fitness
c        CALL sortfit(fit,fitpar)
        do iparm=1,nparm
           modelcur(iparm)=model(iparm,1)
        enddo
        fitcur=fit(1)

c
c---    finished ?
c
        IF (iforwpop.ge.niter)
     &                STOP ' >>>>>>>>> Only one generation'
c
c%%%%%%    loop over iterations
c

        DO iter=1,niter
c          WRITE(prtfil,*)'iteration',iter
          IF (iopt(6).eq.1) WRITE(*,*)'iteration',iter
          IF (mod(iter,10).eq.0)WRITE(*,901)iter,ipop,iforwpop      
 901      Format(' Iteration:',i5,' population:',I5,
     &               ' forward model:',i9)

      write(*,'(a,f8.4,100i4)')
     &          ' curent model',fitcur,(modelcur(iparm),iparm=1,nparm)
c---  fill up posible models
          do iparm=1,nparm
            ihelp=2*iparm
            do ii=1,nparm
              model(ii,ihelp-1)=modelcur(ii)
              model(ii,ihelp)=modelcur(ii)
            enddo
          enddo 
c---  put step length 
c         do iq=1,4
c            write(*,*)'before',model(1,iq),model(2,iq)
c         enddo
          q=0 
          do iparm=1,nparm
            if (tabu(iparm,1).eq.1) then
              q=q+1
              model(iparm,q)= model(iparm,q)-1
            elseif(tabu(iparm,1).eq.-1) then
              q=q+1
              model(iparm,q)= model(iparm,q)+1
            else
              q=q+1
              model(iparm,q)= model(iparm,q)-1
              q=q+1
              model(iparm,q)= model(iparm,q)+1
            endif
          enddo

          iq=0
 502      iq=iq+1
 500      if (iq.gt.q) goto 503
            do iparm=1,nparm
              if (( model(iparm,iq).le.-1) .or. 
     &           model(iparm,iq).ge.ndigit(iparm)) then
                write(*,*)'outside of limits',iq,model(iparm,iq)
                goto 501
              endif
            enddo
            goto 502
 501        do iiq=iq+1,q
              do iparm=1,nparm
                model(iparm,iiq-1)=model(iparm,iiq)   
              enddo
            enddo
            q=q-1
            write(*,*)' one model eliminated'
            goto 500
 503        continue

c
c---- try only using two models
c
         if (mod(iter,10).ne.1) then
c         write(*,*)'iter',iter
         iq=q*ran2(iseed)+0.999
              do iparm=1,nparm
                model(iparm,1)=model(iparm,iq)   
              enddo

              iparm_delt=iparm
              if (model(iparm_delt,1)-modelcur(iparm_delt).eq.delta)then
                 q=1
              else
                 do iparm=1,nparm
                   model(iparm,2)=modelcur(iparm)   
                 enddo
                 model(iparm_delt,2)=model(iparm_delt,2)+delta
                 q=2
              endif
          endif

c         do iq=1,4
c            write(*,*)'after',model(1,iq),model(2,iq)
c         enddo

c          write(*,*)' using',q,' models in iteration',iter 

          DO iq=1,q
c---          set the model and call the modelling routine
            CALL setmodel(iq)
            CALL forw2()
            IF (lwascomp.eq.1) THEN
              iforwtotal=iforwtotal+1
              iforwpop=iforwpop+1
              IF (iopt(3).ge.1) CALL norm   ! normalize the response 
              CALL cost(fit(iq))
c              write(*,*)model(1,iq),model(2,iq),fit(iq)
              if (iopt(15).eq.1) then
c-----          for optimization using hybrid scheme
                call getmodelreal(xtheta)
                CALL gaunew(xtheta(1),5)            !maxiter)
                iforwtotal=iforwtotal+5
                iforwpop = iforwpop+5
                call modeldecode(xtheta,iq)
                IF (iopt(3).ge.1) CALL norm   ! normalize the response 
                CALL cost(fit(iq))
               endif
             ELSE
              fit(iq)=10
            ENDIF 
          ENDDO !iq
c--       sort according to fitness
          CALL sortfit(fit,fitpar)
c
c---      updata tabu list
c
          do iparm=1,nparm
            tabu(iparm,2)=tabu(iparm,2)+1
            if ((model(iparm,1)-modelcur(iparm)).ne.0) then
              tabu(iparm,1)=-(model(iparm,1)-modelcur(iparm))
              tabu(iparm,2)=0
              iparm_delt=iparm
              delta=(model(iparm,1)-modelcur(iparm))
            endif
            if (tabu(iparm,2).gt.ntabu) then
              tabu(iparm,1)=0
            endif
c            write(*,*)'tabu',tabu(iparm,1),tabu(iparm,2)
          enddo
c
c----     make best fit model current model
c
          do iparm=1,nparm
            modelcur(iparm)=model(iparm,1)
          enddo
          fitcur=fit(1)
c
c---      find bsf population
c
          IF (fit(1).lt.objbsf) THEN
            objbsf=fit(1)
            DO i=1,nparm
              modelbsf(i)=model(i,1)
            ENDDO
            WRITE(*,902)objbsf,' at iteration',iter
            WRITE(prtfil,902)objbsf,' at iteration',iter
 902        FORMAT(' new best energy',E12.4,A,I4)
            IF (objbsf.lt.objball) THEN
               WRITE(60,'(1x,I6,I4,G12.5,100I5)')
     &          iforwtotal,1,fit(1),(model(iparm,1),iparm=1,nparm)
            ENDIF
          ENDIF

          IF (iforwpop.gt.niter)THEN
c---        maximum number of generation has been carried out
            GOTO 90
          ENDIF
        ENDDO         ! iter
90      CONTINUE
c
c---    find best of all population
c
        IF (objbsf.lt.objball) THEN
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
            write(10,'(i4,i4,g12.5,100i5)')
     &        i,ipop,fit(i),(model(iparm,i),iparm=1,nparm)
        ENDDO

        iforwtot=iforwtot+iforwpop
        WRITE(*,*)'total calls for population',ipop,iforwpop

        if (ispawn.eq.1) then
          write(*,*)'finnished population',ipop
          stop
        else if (nspawn.gt.1) then
          do i=1,nspawn
            retval=wait(status)
            write(*,*)'retval,staus',retval,status
          enddo
        endif
 99     CONTINUE                ! ipop  number of populations 
100   CONTINUE    ! ipop  number of populations 

c
c---- write out the result
c
      DO 111 ihelp=1,2
        IF (ihelp.eq.1) lun=6
        IF (ihelp.eq.2) lun=prtfil
        WRITE(lun,*)

        WRITE(lun,*)'  best of all energy',objball
        IF (iopt(13).eq.1) THEN
          IF (iopt(18).eq.0) THEN
c             WRITE(lun,900)(1-objball),10*log10(1-objball)
c          ELSEIF IF (iopt(20).eq.0) THEN
c             WRITE(lun,*)objball,objball
          ENDIF
        Endif
        WRITE(lun,*)' best-of all      deviation  '       
        DO i=1,nparm
          WRITE(lun,'((i4,f10.3,7x),3f10.3)')i,
     &      fval(modelball(i),i),fval(modelball(i),i)-xstart(i)
        ENDDO

        WRITE(lun,'(/a,i8)')' Forward  modelling calls:', iforwtot
111   CONTINUE

999      CALL rdtime(t1)
      WRITE(prtfil,310) t1
      WRITE(*,310) t1
 310  FORMAT(/1H ,' Inversion, time: ',F12.3,' secs.')

c      CALL lib$show_timer()
      END




















