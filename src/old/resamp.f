      PROGRAM resamp
c**********************************************************************
c
c**********************************************************************
                                ! Common block definitions are
      USE global
      INCLUDE 'comopt.h'        ! INCLUDEd in these files.
      INCLUDE 'comforw.h'
      REAL fval                 ! function for computation of real values
C**   local variables
      INTEGER iter,iq,i,J,ii,  ! iterative indices.
     >     ifreq                ! frequency index.
      INTEGER ipar(mq),iforwloc
      INTEGER ihelp,lun,index
      INTEGER modelball(mpar)   ! best of all model
      INTEGER modelbsf(mpar)    ! best so far model
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
      INTEGER ix,itrn, fork,ispawn,wait ,retval,status,nspawn,nspawn0
      INTEGER ifrq,idep
      EXTERNAL fork,wait
      INTEGER getpid,id,iterpowell
      CHARACTER*20 COMFILE
      INTEGER qpar,irep,qtemp,iparm,maxiter,
c     error flags for read_cov and read_hp
     >     ierrdat,ierrcov,ierrhp, 
     >     ierrinp,ioer,ierr,ierrD
      data qtemp/1/,maxiter/100/,
c     error flags for read_cov and read_hp
     >     ierrdat/0/,ierrcov/0/,ierrhp/0/, 
     >     ierrinp/0/,ioer/0/,ierr/0/,ierrD/0/
c     INTEGER iter,iforwpop
      COMMON/iterpar/iter,iforwpop 
      INTEGER flagpu
      COMMON /flagpu/flagpu
      INTEGER jparm,M,itmax,mq_post1
c     PARAMETER(M=50) 
      real fret,tol,xstep(mpar)
      REAL xlocobj,xlocmod(mpar)
      REAL,DIMENSION(:), ALLOCATABLE ::  fit
      WRITE(*,*) 'Running SAGA version 5.0, June 2006' 
      xlocobj=1e10
c     
c---- initialization
c
      iWriteTrf = 0
c     OPEN(6, STATUS='NEW',CARRIAGECONTROL='LIST')
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
      write(*,*)' ... exit input'
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

C     Run the forward model using the parameters specified in the run
C     (.DAT) file.  The results of this run can later be written to the
C     OBS file in either covariance matrix or hydrophone DATA FORMAT.
        write(*,*)' enter forwinit ...'

      CALL forwinit
c    
        write(*,*)' ... exit forwinit'
      iforwtot=1

      CALL getmodelreal(xstart)
      CALL modeldecode(xstart,1,ierr)
      DO iparm=1,nparm
         xstar(iparm) = xstart(iparm)
      ENDDO

c     
c
c---  OPEN files
c
c
c---  size of parental distribution
c
      qpar=2*INT(0.5*pu*q)

c---  start the clock
      CALL cltime
      iq = 1
c sep 08 iWriteTrf = 2 was changed to 0  

      iWriteTrf = 0 

      CALL opfilw(16,ioer)
      iforwt = 0

C     Thing starts here .cfh.
      DO ii = 1,1000000
         iforwt = iforwt+1
c     cfh 11/08/2005
         READ(81,*,err=999,end=999)(xtheta(iparm),iparm=1,nparm)
         WRITE(*,*)'ii = ',ii
c     WRITE(*,*)'xtheta = ', (xtheta(iparm),iparm=1,nparm)
        
c---  set the model and CALL the modelling routine
         CALL setmodelreal(xtheta)
         CALL forw2()
c     WRITE(*,*)' finished forward model',iforwpop+1,ipop
         if (iopt(30)==1 .or. iopt(30)==2 .or. 
     *         iopt(30)==9 .or. iopt(30)==11) then  !snap or snaprd
            if (iopt(2) == 3) then ! option e 
               call writecomplex 
            elseif (iopt(13) == 1) then ! option c 
               call writecomplex 
             elseif (iopt(2) == 4) then ! option k 
       write(*,*)' Peter horizontal'
               call writecomplex 
            else
               call writeTL
            endif
         endif
      ENDDO                     !iq
      
 999  CALL rdtime(t1)
      WRITE(prtfil,310) t1
      WRITE(*,310) t1
 310  FORMAT(' Run, time: ',F12.3,' secs.')
      
      END

