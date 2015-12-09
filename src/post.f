      PROGRAM post
c**************************************************************************
c     post processor PROGRAM for inversion of seismo-acoustic DATA 
c     using genetic algorithms.
c     The RESULT is presented as the most likely model, the mean model and the 
c     standard variation.
c     
c     PETER GERSTOFT, 1992 
c**************************************************************************

      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL fval                 ! function for computation of real values
C     
C**   local variables
      INTEGER iq,i,J,jn,jj,k
      INTEGER ihelp,ii,lun,index
      INTEGER maxpdp(mpar)      ! best  model in population
      REAL t1,xbestfit,fit_teo
      REAL fit(mq_post),fitpar(mq_post),fitmin,fit_best(mfreq)
      REAL xmeanmean(meofvar),xmean(mpar),xmeancost,xpdpcost
      REAL xbestbest(meofvar),xstart(mpar),xbest(mpar)
      REAL xmaxpdp(meofvar), fitmin_freq(mfreq),xdum
      REAL*8  svar(mpar),xsum,xsumsqr ! param for statis.
      REAL*8 pdp(mpar,0:(mdig-1)),xhelp
      REAL*8 covar(mpar,mpar)
      REAL*8 objmean,expfit(mq_post),nsearch,totobs     
      INTEGER qpar,jdum,ifreq,ifitmin
      INTEGER npdf(mpar),npdfplot,nmean,iparm,npdfvar,      mq_post1
      REAL   ppd2d(0:(mdig-1),0:(mdig-1)),xdiv,zsum,xsum1
c     CHARACTER*80 char
      INTEGER 
     >     ierrdat,ierrcov,ierrhp, ierraun,
     >     ierrinp,ioer
c     error flags for read_cov and read_hp
      data
     >     ierrdat/0/,ierrcov/0/,ierrhp/0/, 
     >     ierrinp/0/,ioer/0/, ierraun/0/
c     error flags for read_cov and read_hp
      
c---  for a priori knowledge
      REAL a1,a2,a3
c     
c---- initialization
c     
      WRITE(*,*)' You are now running the SAGA postprocessor '
      iWriteTrf=0

      DO i=1,mpar
         DO j=0,mdig-1
            pdp(i,j)=0
         ENDDO
         DO j=1,mpar
            covar(i,j)=0.
         ENDDO
      ENDDO
      DO i=0,mdig-1
         DO j=0,mdig-1
            ppd2d(i,j)=0
         ENDDO
      ENDDO
c     
c---  Set some forward model specific flags
c     
      CALL forwardmodel(iopt,mopt)
c     
c---  READ input DATA
c     
c     OPEN(7, STATUS='NEW',CARRIAGECONTROL='LIST')
      CALL opfilr(1,ioer)
      CALL opfilw(7,ioer)
      CALL input 
c     OPEN(11,STATUS='NEW',CARRIAGECONTROL='LIST',recl=32000)

c     
c---  size of search space
c     
      nsearch=0
      DO iparm=1,nparm
         nsearch=nsearch+LOG10(1.0*ndigit(Iparm))
      ENDDO
      WRITE(*,'(a,f6.2)')' Size of the search space: 10**',nsearch
      WRITE(prtfil,'(a,f6.2)')' Size of the search space: 10**',nsearch
      if (itrans(4).gt.0) then 
         write(*,*) 'allocating weight matrix'
         allocate(weight(mobs))   
         weight=1.0
      endif
      if (iopt(13).gt.0) then
         write(*,*) 'allocating cov matrix'
         allocate(cov(ndep*ndep*nfrq*mx))
      else
         allocate(data(ndep*nfrq*mx))
      endif
C     .cfh.
      write(*,*) 'nx = ',nx
      if (iopt(35) .eq. 1) then
         write(*,*) 'allocating unc matrix'
         allocate(AUN(ndep*ndep*nfrq*nx))
      endif

      IF (iopt(2).GE.1) THEN
c---  reading observed DATA from file 2
c     CALL readdata
         CALL opfilr(2,ierrinp)
         IF (ierrinp.NE.0) STOP ' >>> In file does not exist'    
         
         IF (iopt(13).EQ.1) THEN
            WRITE(*,*)' Reading covariance matrix'
            CALL read_cov(ierrcov)
         ELSEIF (iopt(2).EQ.1) THEN
            CALL readdat2(ierrdat)
         ELSEIF (iopt(2).EQ.2) THEN
            CALL readdata
         ELSEIF (iopt(2).EQ.3) THEN
            CALL read_hp(ierrhp)
         ELSEIF (iopt(2).EQ.4) THEN
            CALL read_ha(ierrhp)
c     CALL readdatiso
         ENDIF
         CLOSE(2)
      ENDIF

C     .cfh.****************************
      IF (iopt(35) .EQ. 1) THEN
         WRITE(*,*)' Reading unc. cov. matrix' ! from UNC file.
         CALL opfilr(4,ierraun)
         CALL READ_AUN(ierraun)
         CALL ROTMAT
         CLOSE(4)
      ENDIF
      IF (ierrinp.GT.0) STOP '>> IN file does not exist'
      IF (ierraun.GT.0) STOP 'error reading UNC file (cov-format)'
      IF (ierrcov.GT.0) STOP 'error reading IN file'
      IF (ierrdat.GT.0) STOP 'error reading IN file'
      IF (ierrhp.GT.0 ) STOP 'error reading IN file'
c     
c---  OPEN files
c     

      CALL opfilw(11,ioer)     
      IF ((iopt(16).EQ.1)) THEN ! contour plot files
         CALL opfilw(28,ioer)
         CALL opfilw(29,ioer)
      ENDIF
      CALL opfilw(19,ioer)      ! fipplot plot files
      WRITE(19,*)'      1024          MODU' ! This is for compatability
      CALL opfilw(20,ioer)
      CALL opfilr(10,ioer)      ! the inversion result
      CALL opfilr(60,ioer)      ! the inversion result
      CALL opfilw(13,ioer)      ! the matlab plotting file
c-----Used options to matlab file
      WRITE(13,'(a,30i3,a)')' iopt=[',(iopt(i),i=1,30),'];'   
      WRITE(13,'(a)')' freq=['
      do i = 1,nfrq
         WRITE(13,'(f15.4)') frq(i)
      enddo
      WRITE(13,'(a)')' ];'
      CALL opfilw(15,ioer)      ! the correlation coefficient file
c     
c---  sort to find how many plots to DO
c     
      npdfplot=1
      ii=1
      ihelp=1
      DO i=2,nparm
c     pg 29/5        IF (par2phy(I).NE.par2phy(i-1)) THEN
         IF ((par2phy(I).NE.par2phy(i-1)).OR.(fmin(i).NE.fmin(i-1))
     1        .OR.(fmax(i).NE.fmax(i-1))) THEN
            npdfplot=npdfplot+1
            npdf(ihelp)=ii
            ii   = 0
            ihelp=ihelp+1
         ENDIF
         ii=ii+1
      ENDDO
      npdf(ihelp)=ii
      npdfvar=npdfplot
C*****We always plot them individually, to avoid problems WITH different 
c*****scaling
c---- This is for plotting single trace pdp
c     DO i=2,nparm
c     npdf(i)=1
c     ENDDO
c     npdfplot=nparm
c************
c     WRITE(*,*)npdfplot,(npdf(i),i= 1,npdfplot)
c     
c---  initialize forward model
c     
c     
C**********************************
      IF (itrans(4).GE.1) THEN  !weight for each data point
         CALL opfilr(3,ioer)
         CALL Read_weight
         CLOSE(3)
      ENDIF
C     
      write(*,*) 'before forw'
      CALL forwinit
      write(*,*) 'after forw'
      
      CALL getmodelreal(xstart)

      DO iparm = 1,nparm
         xprior(iparm)=xstart(iparm)
      ENDDO

      IF (iopt(2).EQ.0) THEN
c---  move DATA from resp to DATA
         DO j=1,ndep
            index=(j-1)*nx
            DO i=1,nx
               DATA(i+index)=resp(i+index)
            ENDDO
         ENDDO
      ENDIF
c     
c     DO i=1,7
c     WRITE(*,*)'data 0',i,DATA(i)
c     ENDDO
      IF (iopt(3).GE.1) CALL normdata
      
      IF (iopt(3).GE.1) CALL norm ! normalize the response 

      write(*,*)' before cost'
      CALL cost(fit(1))
      WRITE(*,*)' initial fitness',fit(1)
      WRITE(prtfil,*)' initial fitness',fit(1)
      
c     DO i=1,7
c     WRITE(*,*)'data 1',i,DATA(i)
c     ENDDO
c     
c---  energy normalization
c     
c     
c---  writeout
      WRITE(prtfil,*)
      WRITE(prtfil,*)' crossover probability    ',px
      WRITE(prtfil,*)' mutation  probability    ',pm
      WRITE(prtfil,*)' update    probability    ',pu
      WRITE(prtfil,*)' number  of iterations    ',niter
      WRITE(prtfil,*)' size   of populations    ',q
      WRITE(prtfil,*)' number of populations    ',npop
      WRITE(prtfil,*)
      
      IF (iopt(8).EQ.1) THEN
c**** plot ambiguity
         WRITE(*,*)' contour plot, option C,' 
         WRITE(*,*)' is plotted directly from SAGA'
         WRITE(*,*)' Stopping; no need for running POST'
         STOP
      ENDIF
      
      CALL cltime

      qpar=2*INT(0.5*pu*q)
c*****warning q is redefined
c     q=npop*q
c     WRITE(*,*)' Warning q is redefined to',q
c     IF (q.GT.mq_post) THEN
c     STOP 'increase mq'
c     ENDIF

c     
c---  READ from *.mat file
c     
      WRITE(*,*)' >>>> reading model vectors'
c---  allocate
      mq_post1=(niter+q)*npop
      allocate(model(nparm+2,mq_post1))
 21   DO  i=1,mq_post
c     DO  i=1,q
         READ(10,*,END=999,err=998)
     &        model(nparm+1,i),model(nparm+2,i),fit(i),
     &        (model(jn,i),jn=1,nparm)
c     &        idum,jdum,fit(i),(model(jn,i),jn=1,nparm)
c     READ(10,'(a80)',END=999,err=998)char
c     WRITE(*,*)char
c     READ(CHAR(10:80),*,END=999,err=998)
c     &        jdum,fit(i),(model(jn,i),jn=1,nparm)
c     &        CHAR(1:9),jdum,fit(i),(model(jn,i),jn=1,nparm)
c     WRITE(60,*)fit(i)
       if (i.eq.mq_post1) then
          mq_post1=  mq_post1+ mq_post1
          WRITE(*,*)'MODEL SIZE INCREASED!'
          deallocate(model)
          allocate(model(nparm+2,mq_post1))
          rewind(10)
          goto 21
       endif
         IF (i.EQ.mq_post) THEN
            WRITE(*,*)' post: parameter mq_post not large ',
     &           'enough to obtain'
            WRITE(*,*)' all observations, increase mq_post=',mq_post
            STOP
         ENDIF
         IF (MOD(i,100000).EQ.1)
     &        WRITE(*,*) i,jdum,fit(i),(model(jn,i),jn=1,nparm)
      ENDDO
      
      
 998  WRITE(*,*)'problems reading matfile: ',i,'lines read'
      WRITE(*,*) jdum,fit(i),(model(jn,i),jn=1,nparm)
      STOP

 999  CONTINUE
      IF ((i-1).EQ.0) THEN
         STOP' No indivduals is read from mat-file'
      ELSEIF ((i-1).LE.npop*niter) THEN
         WRITE(*,*)' WARNING: number of read individuals',i-1
         WRITE(*,*)' is less than npop*niter=',npop*niter
      ELSEIF ((i-1).GE.(npop*(q+niter))) THEN
         WRITE(*,*)' WARNING: number of read individuals',i-1
         WRITE(*,*)' is too much larger than npop*niter=',npop*niter
      ENDIF
      q=i-1
      WRITE(*,*)'number of read model vectors',i-1

c---  sort the models according to the value of the object FUNCTION
      IF (q.LT.200000) THEN

c     nparm=nparm+2
         IF (iopt(28).EQ.1) CALL elimiatemodels(q,fit)
         CALL sortfitPost(fit,fitpar)
c---  WRITE sorted output to *b.mat
c     
         DO  i=1,q
            WRITE(11,'(i6,i4,g12.5,100i5)')
     &           i,1,fit(i),(model(jn,i),jn=1,nparm)
         ENDDO
c     nparm=nparm-2
         CLOSE(11)
         WRITE(*,*)' The models have been sorted and written to *b.mat'
         fitmin=fit(1)
         ifitmin=1
c     
c-----for NEW post proc.
c     
         IF ((iopt(28).EQ.2)) THEN
            q=q*0.5
            WRITE(*,*)' For unscaled ppd q is reduced, q=',q
         ENDIF
         IF (iopt(28).EQ.1) THEN
            CALL setmodel(ifitmin)
            CALL forw2
            IF (iopt(3).GE.1) CALL norm ! normalize the response 
            CALL cost(xdum)
            xsum=0
            DO ifreq=1,nfrq
               IF (iopt(5).EQ.3) THEN
                  fitmin_freq(ifreq)=  fit_bart_ind(ifreq)/nx*5
               ELSEIF (iopt(5).EQ.4) THEN
                  fit_best(ifreq)=   fit_bart_ind(ifreq)
                  fitmin_freq(ifreq)=  fit_bart_ind(ifreq)/nummodes !ndep 
                  xsum=xsum+fitmin_freq(ifreq)
               ELSE             ! covariance matrix
                  fit_best(ifreq)=   fit_bart_ind(ifreq)
                  fitmin_freq(ifreq)=  fit_bart_ind(ifreq)/nummodes !ndep 
                  xsum=xsum+fitmin_freq(ifreq)
c     WRITE(*,*) 'post:fitmin',ifreq,  zcov_noise(ifreq)
c     WRITE(*,*) 'post:fitmin_cov',ifreq, fitmin_freq(ifreq)
               ENDIF
            ENDDO 
            xsum=xsum/nfrq
         ENDIF
      ELSE
         WRITE(*,*)' The model vectors has NOT been sorted'
         fitmin=10e6
         DO i=1,q
            IF (fitmin.GT.fit(i)) THEN
               fitmin=fit(i)
               ifitmin=i
            ENDIF
         ENDDO
      ENDIF
c*****warning q is redefined, so that we ONLY USE the best half 
c*****of the selected DATA
c     q=0.5*q
c     WRITE(*,*)' Q is not changed'
      npop=q
      objmean=0.
      nmean=MIN(q,50)
      DO  i=1,nmean
         objmean=objmean+fit(i)
      ENDDO
      objmean=objmean/nmean
c     objmean=objmean-fit(ifitmin)

c     
c     
c---  weigthing for the probability distributions
c     
      xdiv=objmean-fitmin
      IF (xdiv.LE.0) THEN
         WRITE(*,*)' All observations have same energy' 
         xdiv=1
      ENDIF
c     xdiv=1
      DO  ipop=1,q
c     expfit(ipop)=EXP(-fit(ipop)/objmean)  
         IF((iopt(28).EQ.0) .AND. (q.LT.200000)) THEN
c---- classical post processing
            expfit(ipop)=EXP(-(fit(ipop)-fit(ifitmin))/xdiv )  
         ELSEIF ((iopt(28).EQ.2)) THEN
            expfit(ipop)=1
         ELSE
            IF (nfrq.EQ.1) THEN
               expfit(ipop)=EXP(-(fit(ipop)-fit(1))/fitmin_freq(1))
c     IF (MOD(ipop,200).EQ.0)
c     &           WRITE(*,*)' fitmin',ipop, expfit(ipop), fitmin_freq(1)
            ELSE
               IF (iopt(5).EQ.4) THEN
                  
                  READ(60,*,END=999,err=998)
     &                 xdum,xdum,xdum,(xdum,jn=1,nparm),
     &                 (fit_bart_ind(ifreq),ifreq=1,nfrq)
                  WRITE(*,*)'forwardmodel',ipop,
     &                 (fit_bart_ind(ifreq),ifreq=1,nfrq)
               ELSE
                  CALL setmodel(ipop)
                  CALL forw2
                  IF (iopt(3).GE.1) CALL norm ! normalize the response 
                  CALL cost(xdum)
                  IF (ipop.LT.100)
     &                 WRITE(*,*)'forwardmodel',ipop,xdum,fit(ipop)
               ENDIF
               expfit(ipop)=1d0 
               xsum1=0
               DO ifreq=1,nfrq
                  expfit(ipop)= expfit(ipop)
     &                 *dexp(DBLE(-(fit_bart_ind(ifreq)-
     &                 fit_best(ifreq))
     &                 /fitmin_freq(ifreq))) 
                  xsum1=xsum1+fit_bart_ind(ifreq)-fit_best(ifreq)
c     WRITE(71,*)ipop,expfit(ipop),fit_best(ifreq),
c     &                 fit_bart_ind(ifreq),fitmin_freq(ifreq),xsum1,xsum
               ENDDO 
c     expfit(ipop)= EXP(-xsum1/xsum)
            ENDIF
         ENDIF
c     expfit(ipop)=1
      ENDDO    

      IF (iopt(23).EQ.1) THEN
c     
c---- adding a priori distributions
c     
         DO iparm=1,nparm
            a1=fmin(iparm)-xstart(iparm)
            a2=fmax(iparm)-xstart(iparm)
            WRITE(*,*)'a1,a2',a1,a2
            IF (iapri(iparm).GT.0) THEN
               WRITE(*,*)' Prior knowledge is for paramater',iparm
               DO iq=1,q
                  a3=fval(model(iparm,iq),iparm)
                  WRITE(*,*)'a3',a3
                  IF ((a1*(a3-xstart(iparm))).GE.0) THEN ! it belongs to lower part
                     expfit(iq)= expfit(iq)*(a3-fmin(iparm))/(-a1)
                  ELSE          ! then it is lower part
                     expfit(iq)= expfit(iq)*(a3-fmax(iparm))/(-a2)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDIF


c     
c---  for probability of a model set
c     
      WRITE(*,*)' >>> estimating probability...' 
      DO 100 ipop=1,npop 
         DO jn=1,nparm 
            pdp(jn,model(jn,ipop))=
     &           pdp(jn,model(jn,ipop))+expfit(ipop)
         ENDDO                  !jn
 100  CONTINUE                  ! ipop  number of populations 
c     
c---  find marginal probability
c     
      WRITE(*,*)'normalizing distributions'
      DO i=1,nparm
c---  normalize the pdp
         xhelp=0.
         DO j=0,ndigit(i)-1
            xhelp=xhelp+pdp(i,j)
         ENDDO       
         totobs=xhelp
c     WRITE(*,*)i,xhelp,ndigit(i)
c     xhelp=1
         DO j=0,ndigit(i)-1
            pdp(i,j)=pdp(i,j)/xhelp
         ENDDO       
      ENDDO       

      WRITE(*,*)'marginal 1d probabilites found'
c     
c---  Find marginal 2d-Ppd
c     
      IF (iopt(16).EQ.1) THEN
         DO  ipop=1,npop 
            ppd2d(model(ippd1,ipop),model(ippd2,ipop))= 
     &           ppd2d(model(ippd1,ipop),model(ippd2,ipop)) 
     &           +expfit(ipop)
         ENDDO	
c     +      normalize
         zsum=0
         DO j=0,ndigit(ippd1)-1
            DO  i=0,ndigit(ippd2)-1
               ppd2d(j,i)=ppd2d(j,i)/xhelp
               zsum=MAX(zsum,ppd2d(j,i))
            ENDDO
         ENDDO
         DO j=0,ndigit(ippd1)-1
            DO  i=0,ndigit(ippd2)-1
               ppd2d(j,i)=ppd2d(j,i)/zsum
            ENDDO
         ENDDO
         CALL conppd2d(ppd2d) 
      ENDIF                     !End Find marginal 2d-Ppd
c     
c---- for standard deviation
c     
      
      DO i=1,nparm
         xsum=0.
         xsumsqr=0.
         DO j=0,ndigit(i)-1
            xhelp=fval(j,i)*pdp(i,j)
            xsum=xsum+xhelp
            xsumsqr=xsumsqr+xhelp*fval(j,i)
         ENDDO
         xmean(i)=xsum
         svar(i)=SQRT(ABS(1.00000000000001*xsumsqr-xsum**2)) ! not correct for non-norm
         IF (xmean(i).NE.0)svar(i)=svar(i)/ABS(xmean(i))
      ENDDO

c     
c---  most likely values
c     
      WRITE(*,*)' Finding most likely values'
      DO i=1,nparm
         xhelp=0.
         DO j=0,ndigit(i)-1
            IF (xhelp.LT.pdp(i,j)) THEN
               xhelp=pdp(i,j)
               maxpdp(i)=j
            ENDIF 
         ENDDO      
      ENDDO      

c     
c---  forward model for most likely vector
c     
      WRITE(*,*)'setting model for maxppd parameters'
      WRITE(prtfil,*)'setting model for maxppd parameters'
      CALL setmodelbest(maxpdp)
      CALL forw2
      IF (iopt(3).GE.1) CALL norm ! normalize the response 
      CALL cost(xpdpcost)
c     
c---  forward model mean vector
c     
      WRITE(*,*)'setting model for mean parameters'
      WRITE(prtfil,*)'setting model for mean parameters'
      CALL setmodelreal(xmean)
      CALL forw2
      IF (iopt(3).GE.1) CALL norm ! normalize the response 
      CALL cost(xmeancost)
c     
c---  construct the covariance matrix
c     
      WRITE(*,*)'parameter covariance is estimated...'
      DO 101 ipop=npop,1,-1 
         DO jn=1,nparm
            DO jj=1,nparm
               covar(jn,jj)=covar(jn,jj)+
     &              1.D0*(fval(model(jn,ipop),jn)/xmean(jn)-1d0)
     &              *(fval(model(jj,ipop),jj)/xmean(jj)-1d0)
     &              *(expfit(ipop)/totobs)
c     &        -1.D0/npop        
            ENDDO
         ENDDO                  !jn
 101  CONTINUE                  ! ipop  number of populations 
      DO i=1,nparm
         svar(i)=SQRT(ABS(covar(i,i)))
      ENDDO
      DO i=1,nparm
         DO ii=1,nparm
            IF (svar(i).NE.0 .AND. svar(ii).NE.0) THEN
c               covar(i,ii)=ABS(covar(i,ii))/svar(i)/svar(ii) 
               covar(i,ii) = covar(i,ii)/svar(i)/svar(ii) 
            ELSEIF (svar(i).EQ.0 .AND. svar(ii).EQ.0) THEN
               covar(i,ii)=1.0
            ENDIF
         ENDDO
      ENDDO

c---  renormalize standard deviation
      DO i=1,nparm
         svar(i)=svar(i)*ABS(xmean(i)/(fmax(i)-fmin(i)))
      ENDDO
c     
c---- WRITE out coovariance      (DISABLED)
c     
      DO i=1,nparm
         WRITE(15,'(20f15.6)')(covar(i,ii),ii=1,nparm)
      ENDDO

c     
c---- map the shape FUNCTION to the model parameters
c     for most likely, mean and best
c     
      DO i=1,nparm
         IF (par2phy(i).EQ.11) THEN
            aeof(par2lay(i))=fval(maxpdp(i),i)
         ENDIF
      ENDDO
      CALL eofvalpoint(xmaxpdp)

      DO i=1,nparm
         IF (par2phy(i).EQ.11) THEN
            aeof(par2lay(i))=xmean(i)
         ENDIF
      ENDDO
      CALL eofvalpoint(xmeanmean)

      DO i=1,nparm
         IF (par2phy(i).EQ.11) THEN
            aeof(par2lay(i))=Fval(model(i,1),i)
         ENDIF
      ENDDO
      CALL eofvalpoint(xbestbest)
c---- plot the bedst RESULT

      iwrite_soundsp=1  
      WRITE(prtfil,*)' Calling forwardmodel with best model' 
      WRITE(*,*)' Calling forward model with best model ...'
      CALL setmodel(ifitmin)
      CALL forw2
      iwrite_soundsp=0
      IF (iopt(3).GE.1) CALL norm ! normalize the response 
      CALL cost(xbestfit)
      WRITE(*,*)'the best energy from post', xbestfit
      WRITE(*,*)'the best energy from inv ', fit(ifitmin)
      IF ( (xbestfit.GT. 1.001*fit(ifitmin)) .OR. 
     &     (xbestfit.LT. 0.999*fit(ifitmin))) THEN
         WRITE(*,*)'**********'
         WRITE(*,*)' The fitness from post for the best model'
         WRITE(*,*)' does not correspond to that from the inversion'
         PAUSE                  !'enter "continue" to proceed'
      ENDIF
      WRITE(*,*)'iopt(10)',iopt(10)
      write(*,*)'number of ranges = ',nx
      IF (iopt(10).GE.1 ) THEN
         IF ((iopt(13).EQ.1 ) .OR. (iopt(5).EQ.4) 
     &        .OR. (iopt(5).EQ.6) )  THEN
            CALL pltphase
         ELSEIF (iopt(5).EQ.5)  THEN
            CALL pltphaseF
         ELSEIF (iopt(5).EQ.7)  THEN
            CALL pltphasek
         ELSEIF ((iopt(5).EQ.1).OR. (iopt(5).EQ.3)  
     &           .OR. (iopt(5).EQ.0) .OR. (iopt(5).EQ.8)
     &           .OR. (iopt(5).EQ.10)) THEN
            IF (nx.LE.1) THEN
c     WRITE(*,*)' >> Only one range no plot vs range'
               CALL pltlos_dep
            ELSE
               WRITE(*,*)' calling pltlos'
               WRITE(prtfil,*)' calling pltlos'
               CALL pltlos
            ENDIF
         ELSEIF ((iopt(5).EQ.5)) THEN 
            CALL plt_tim
         ENDIF
      ELSE
         WRITE(13,*)'nplots=0;'
      ENDIF
c---  forward model  mean vector
c     
c     WRITE(*,*)' ***** PLOTTING WITH GA-MEAN PROFILE *****'
c     WRITE(prtfil,*)'setting model for mean parameters'
      CALL setmodelreal(xmean)
      CALL forw2
      IF (iopt(3).GE.1) CALL norm ! normalize the response 
      CALL cost(xmeancost)

c     
c---  plot of the pdp-functions
c     
c     WRITE(*,*)' ***********************************************'
c     WRITE(*,*)' **** scalepdp has been commented out **********'
      CALL scalepdp(pdp,maxpdp)
c     
c>>>> this is for plotting a pdp for each PARAMETER-TYPE
c     
      IF (iopt(14).EQ.1) THEN
         ii=1
c     individual
         IF (npdfvar.EQ.1) THEN
            npdfvar=1
            iopt(14)=0
            WRITE(13,*)'unitplot=1;'
         ELSE
            npdfplot=nparm
         ENDIF
         DO i=1,npdfplot        ! this loops over each parameter type
            ihelp=ii+npdf(i)
c     individual
            ihelp=ii+1
c---- plot plp
            CALL PLtpdp(ii,ihelp,pdp)
            IF (iopt(19).EQ.1) THEN !   starting values are added to the pdp 
               CALL pltmod(-1,ii,ihelp,xstart)
            ENDIF 
            ii=ihelp
         ENDDO
      ELSE
c>>>> Plotting all parameters on the same plot
c-----pltmod(# of best curves, first parm,last parm, reference solution)
c     # of best curves =-1 => ONLY xstart is plottet         
c     CALL pltmod(20,1,nparm,xstart) ! plotting the 20 best profiles
         WRITE(*,*)' calling pltpdp'
         CALL PLtpdp(1,nparm,pdp)
         WRITE(*,*)' pltpdp was called'
         
         IF ((iopt(2).GT.0).AND.(iopt(19).EQ.0)) THEN
c---  The best values are plotted on each ppd plot
            DO iparm=1,nparm
               xbest(iparm)=fval(model(iparm,1),iparm)
            ENDDO
            CALL pltmod(-1,1,nparm,xbest)
         ELSE
c---  The inupt values  is plotted on each to ppd plot
            CALL pltmod(-1,1,nparm,xstart)
         ENDIF
      ENDIF

c     --- WRITE best results to *.m file
      WRITE(13,*)' res=[ '
      DO k=1,nparm      
         WRITE(13,'(100f15.3)') fval(model(k,1),k),xmean(k),svar(k)
      ENDDO
      WRITE(13,*)' ]; '
      WRITE(13,*)' bestfit=',fit(ifitmin),';'
c---  .cfh. 2007     
      WRITE(13,*)' meanfit=', xmeancost,';'
      WRITE(13,*)' pcov = ['
      DO i=1,nparm
         WRITE(13,'(20f10.5)')(covar(i,ii),ii=1,nparm) 
         WRITE(13,*)';'
      ENDDO
      WRITE(13,*)'];'
c---
      WRITE(*,*)' writing to results.m !!!'
c     
c---  WRITE the best RESULT to RESULTS.M
c     
c     this has been moved to Ga.f     
c     
c---  WRITE the best RESULT to RESULTS.M
c     
c     OPEN(unit=8,file='results_mean.m',access='append',
c     &           carriagecontrol='list',recl=32000,status='unknown')
c     WRITE(8,*)' %%%%%%%%%%%%%%% '
c     WRITE(8,*)' res(:,  )=['
c     DO k=1,nparm      
c     WRITE(8,'(100f10.2)')
c     &      (xmean(k),k=1,nparm)
c     ENDDO 
c     WRITE(8,*)' ]'';'
c     CLOSE(8)
c     
c---- Writeout to log-file
c     
c     WRITE(*,*)'...results_mean.m is written !'


c     OPEN(unit=9,file='res.m',status='unknown')
c     WRITE(9,'(a)') 'vel(1,:)=['
c     DO i=1,neofvar
c     WRITE(9,'(f12.1)')   xbestbest(i)
c     ENDDO
c     WRITE(9,'(a)') '];'

c     &                   xbestbest(neofvar),'];'
c     CLOSE(9)
c     WRITE(*,*)'...res.m is written !'
c%%%%%%OPEN(unit=9,file='../fit.m',access='append',
c     &           carriagecontrol='list',status='unknown')
c     WRITE(9,'(f8.4)') fit(ifitmin)
c     CLOSE(9)
c     OPEN(unit=9,file='../ran.m',access='append',
c     &           carriagecontrol='list',status='unknown')
c     DO i=1,nparm
c     IF (par2phy(i).EQ.9) THEN
c     WRITE(9,'(f8.1)')  fval(model(i,1),i)
c     ENDIF
c     ENDDO
c     CLOSE(9)
c     WRITE(*,*)'...ran.m is written !'

cf90: position='append', f77: access='append',
      OPEN(unit=8,file='results',position='append',
     &     status='unknown')

c     WRITE(8,*)
      WRITE(8,'(a80)')title
      WRITE(8,*)'***************'

      DO 111 ihelp=1,3
c     USE lun=5 on (DEC alpha)
C     and  lun=6 on (SUN, Linux, SGI)
         IF (ihelp.EQ.1) lun=6
         IF (ihelp.EQ.2) lun=prtfil
         IF (ihelp.EQ.3) lun=8
         IF ((iopt(13).EQ.1) .OR. (iopt(5).EQ.4)) THEN
            IF (iopt(18).EQ.0) THEN
c     WRITE(lun,*)' best of all barlett- power',10*LOG10(1-fit(1))
            ENDIF
         ENDIF
         CALL cost_max(fit_teo)
         WRITE(lun,*)' best obtainable theoretical fit', fit_teo
         WRITE(lun,*)' best fit (best, ppd, mean)',fit(ifitmin),
     &        xpdpcost, xmeancost
         WRITE(lun,'(/,31x,a)')
     &        'best of all  most likely      mean   std-dev'       
         IF (iopt(12).EQ.0) THEN ! two index for optimization variable
            DO i=1,nparm
               WRITE(lun,'((i3,x,a26,i2,f10.3,3x),3f10.3)')
     &              i,phystxt2(par2phy(i)),par2lay(i),
     &              fval(model(i,ifitmin),i),
     &              fval(maxpdp(i),i),xmean(i),svar(i)
            ENDDO
         ELSE                   !
            DO i=1,nparm
               WRITE(lun,'((i3,x,a26,2i3,f10.3,3x),3f10.3)')
     &              i,phystxt2(par2phy(i)),par2lay(i),par3(I),
     &              fval(model(i,ifitmin),i),
     &              fval(maxpdp(i),i),xmean(i),svar(i)
            ENDDO
         ENDIF

         IF (iopt(17).EQ.1) THEN
 1000       FORMAT(/,1x,'- Values for the ',I2,
     &           ' eof''s points from the',
     &           i2,' eof''s functions')
c            WRITE(lun,*)'                            ',
c     &           'best-of all  most likely    mean '    
         WRITE(lun,'(/,31x,a)')
     &        'best of all  most likely      mean'          
            DO i=1,neofvar
               WRITE(lun,'((i3,x,a26,i2,f10.3,3x),2f10.3)')i,
     &              phystxt2(par2phy_eof(i)),par2lay_eof(i),
     &              xbestbest(i),xmaxpdp(i),xmeanmean(i)
            ENDDO
         ENDIF

c---  WRITE out deviation
         IF ((iopt(2).GT.0) .AND.(lun.NE.8)) THEN
            WRITE(lun,*)
            WRITE(lun,*)' Deviation from initial model:'
            DO i=1,nparm
               WRITE(lun,'((i3,x,a26,f10.3,3x),3f10.3)')
     &              i,phystxt2(par2phy(i)),
     &              fval(model(i,1),i)-xstart(i),
     &              fval(maxpdp(i),i)-xstart(i),
     &              xmean(i)-xstart(i),svar(i)
            ENDDO
         ENDIF
 111  CONTINUE
      

      WRITE(8,*)
      CLOSE(8)

      IF ((iopt(13).EQ.1) .OR. (iopt(5).EQ.4)) THEN
         CALL costbart
      ENDIF

c%%%%%%%This is for computing the fit for the best models
c     DO iq=1, 1        !q
c     WRITE(*,*)' Calling forward model with  model ...',iq
c     CALL setmodel(iq)
c     CALL forw2
c     IF (iopt(3).GE.1) CALL norm   ! normalize the response 
c     CALL cost(fit(iq))
c     WRITE(*,*)' fitness', fit(iq)
c     ENDDO
c%%%%%%%
c     WRITE(12,*)' *****',iq,fit(iq)
c     DO ifreq=1,nfrq
c     WRITE(12,*)' **',iq,ifreq, frq(ifreq)
c     DO idep = 1, Ndep
c     i = (ifreq-1)*Ndep + idep
c     WRITE( 12, * ) idep, resp(1+(i-1)*1)
c     END DO
c     END DO
c     WRITE(12,*)
c     
c---- for chi-sqared test
c     
c     IF (iopt(13).NE.1) THEN
c     CALL chisqr(xmean,xmeancost)
c     ENDIF
c     
c     
c-----for plotting velocity profile
c     
c     WRITE(*,*)'setting mean model...'
c     CALL setmodelreal(xmean)
c     CALL plotvelprof
c     
c     
c---- map the shape FUNCTION to the physical values
c     for most likely, mean and best
c     
      IF (iopt(30).EQ.6) THEN
         WRITE(*,*)' calling plotstddev...'
c     CALL   plotstddev(expfit,totobs)     
      ENDIF
c     
c     
      IF (iopt(30).EQ.7 .OR. iopt(30).EQ.1 
c     oases      1            .OR. iopt(30).EQ.3 
     1     .OR. iopt(30).EQ.6  .OR. iopt(30).EQ.9 ) THEN
      CALL opfilw(16,ioer)
      WRITE(prtfil,*)' Calling forwardmodel with best model' 
      WRITE(*,*)' Calling forward model with best model'
      WRITE(*,*)' to write trf file'
      iWriteTrf=1
      CALL setmodel(ifitmin)
      CALL forw2
      ENDIF

c     iwrite_soundsp=0
c     IF (iopt(3).GE.1) CALL norm ! normalize the response 
c     CALL cost(xbestfit)
c     WRITE(*,*)'the best energy from post',xbestfit
c     WRITE(*,*)'the best energy from inv ', fit(ifitmin)
c     IF ( (xbestfit.GT. 1.001*fit(ifitmin)) .OR. 
c     &        (xbestfit.LT. 0.999*fit(ifitmin))) THEN
c     WRITE(*,*)'**********'
c     WRITE(*,*)' The fitness from post for the best model'
c     WRITE(*,*)' does not correspond to that from the  inversion'
c     PAUSE
c     ENDIF


      CALL rdtime(t1)
c     CALL lib$show_timer()
      WRITE(prtfil,310) t1
      WRITE(*,310) t1
 310  FORMAT(' Postprossing, time: ',F12.3,' secs.')
      END
