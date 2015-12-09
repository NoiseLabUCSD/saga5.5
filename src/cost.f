c***************************
      SUBROUTINE cost(fit)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL fit,fit2
      REAL*8 fitlong,e1,e2,e3,fitlag,fittemp
      INTEGER i,j,ix,ifrq,ibart,j1start,jdep,index,index2
      INTEGER index3
      COMPLEX*16 cc,cx,cc_diag,fitlagcmplx
c---  variables for regularization
      REAL regsum,xcurr(mpar)
      INTEGER iparm,ilag,iminlag,ilow,iup,mprof
      INTEGER iprof
      REAL mexcess,metcon,basehvar(10000)
      COMMON/ mexc/ mexcess,metcon,basehvar,iprof
      REAL xonel,xoneu,xu,xl
      REAL amp,tau
      PARAMETER(xonel=0.999,xoneu=1.00001)

C     .cfh.****************************
      IF (iopt(35) .EQ. 1 .and. iopt(8) .EQ. 3) THEN
         DO 200 ifrq=1,nfrq
            DO 201 ix=1,nx 
               ibart  = ix+((ifrq-1))*nx ! which observation
               index3 = (ibart-1)*ndep-1 
               e1     = 0.d0
               j1start= ((ifrq-1))*ndep ! where the response starts
               DO j=j1start,ndep+j1start-1
                  index = ix+(j)*nx
                  e1 = e1+(ABS(resp(index)))**2 
               ENDDO
               DO j=j1start,ndep+j1start-1
                  index = ix+(j)*nx
                  resp(index) = resp(index)/sqrt(e1)
               ENDDO
 201        CONTINUE            !nx
 200     CONTINUE               !nfrq 
      ENDIF
      IF (iopt(35) .Eq. 1) THEN
         CALL ROTVEC
      ENDIF


C**********************************
c     
c      WRITE(*,*)' Entering cost'
C      WRITE(*,*)'lwascomp', lwascomp
      IF (lwascomp.NE.1) THEN
         fit=10
         RETURN
      ENDIF
      fitlong=0.
      IF (iopt(18).EQ.1) THEN   ! Change to dB scale
         fitlong=1d0
      ENDIF
      fit=0.
      e1=0
      e2=0
      e3=0
      cc=(0.d0,0.d0) 
c     WRITE(*,*) 'calling costmbeam iopt,isubopt',
c     1    iopt(30),isubopt(30)
      
      IF ((iopt(30).EQ.1).AND.( isubopt(30).EQ.1)) THEN
c     WRITE(*,*) 'calling costmbeam'
         CALL costmbeam(fit)
         GOTO 900
      ENDIF

      IF (iopt(13).EQ.1) THEN
c     
c-------BARTLETT estimator
c     
         DO 100 ifrq=1,nfrq
            DO 101 ix=1,nx 
               ibart  =ix+((ifrq-1))*nx ! which observation
               index3 =(ibart-1)*ndep-1 
c>>>> now for each covariance matrix
               cc     =(0.d0,0.d0) 
               cc_diag=(0.d0,0.d0) 
               e1     =0.d0
               j1start=((ifrq-1))*ndep ! where the response starts
               DO j=j1start,ndep+j1start-1
                  index =ix+(j)*nx
                  e1 =e1+(ABS(resp(index)))**2 
               ENDDO
               xrepl_sum(ibart)=e1
               DO i=1,ndep
                  cx=(0.d0,0.d0)
                  DO j=i+1,ndep
                     index =ix+(j+j1start-1)*nx
c          WRITE(*,*)'pg ix,j,j1start,nx,index',ix,j,j1start,nx,index
                     index2=i+(j+index3)*ndep 
                     cx    =cx+cov(index2)*resp(index)     
c                WRITE(*,*)'cx',cx,cov(index2),resp(index),index2,index
                  ENDDO
                  index  =ix+(i+j1start-1)*nx
                  index2 =i+(i+index3)*ndep 
                  cc_diag=cc_diag+cov(index2)*ABS(resp(index))**2     
                  cc     =cc+CONJG(resp(index))*cx
               ENDDO
               cc             =2.*cc +cc_diag
               xcor_sum(ibart)=cc
               IF (iopt(35) .EQ. 1 .and. iopt(8) .EQ. 3) THEN
                  write(*,*) 'e1',e1
                  cc =cc/(e1**2) ! normalize with replica
               ELSE
                  cc =cc/e1 ! normalize with replica
               ENDIF
                  fittemp=xcov_trace(ifrq)-dreal(cc)
                  if (fittemp.le.0) fittemp=0
               fit_bart_ind(ibart)= fittemp
               IF (iopt(18).EQ.0) THEN
                  fitlong = fitlong+fit_bart_ind(ibart)
               ELSE IF (iopt(18).EQ.1) THEN
                  fitlong = fitlong*fit_bart_ind(ibart)
c     cfh 2005/1/23   fitlong=fitlong*(fit_bart_ind(ibart)**(1./nbart))
               ELSE IF  (iopt(18).EQ.2) THEN
                  STOP 'cost: iopt(18) not defined'
               ENDIF
               
 101        CONTINUE            !nx
 100     CONTINUE               !nfrq 
         
         IF (iopt(18).EQ.1) THEN 
            fit = fitlong**(1./nbart) ! cfh geometric mean
c     cfh 2005/1/23    fit = fitlong
         ELSE
            fit = fitlong/nbart ! cfh arithmetic mean
         ENDIF
         GOTO 900
      ENDIF                     ! FINISHED WITH COVARIANCE
C     
c     
c---- Here starts the evaluation based on REAL observations
c     WITH PHASE information
c     
      IF (iopt(5).EQ.4) THEN    ! option f
c---- Incoherent addition over frequencies, PHASES information
         DO ifrq=1,nfrq
            DO ix=1,nx            
               e1=0.
               e2=0.
               cc=0.
               DO jdep=1,ndep
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  e1=e1+ABS(resp(ix+index))**2 ! resp
                  e2=e2+ABS(DATA(ix+index))**2 ! data
                  cc=cc+CONJG(DATA(ix+index))*(resp(ix+index))
c     WRITE(99,*)REAL(resp(ix+index)),imag(resp(ix+index))
c     WRITE(*,*)'cost,data',i,
c     1         (DATA(i+index)),(resp(i+index))
               ENDDO            ! ndep--loop
c     WRITE(98,*)ABS(cc)
               IF (e1.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: All calculated data are zero '
                  WRITE(*,*)' cost: All calculated data are zero '
                  IF (iopt(31).EQ.0) THEN
                     WRITE(*,*)'...Stopping from cost'
                     STOP
                  ELSE
                     fit =1e10
                     RETURN
                  ENDIF
               ELSEIF (e2.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: observed data are zero'
                  WRITE(*,*)' cost: observed data are zero'
               ENDIF
               
               IF (isubopt(5).EQ.0) THEN              
c     ---        each freq is weigthed according to received power
                  fit_bart_ind(ifrq)=e2-ABS(cc)**2/(e1)
                  fit =fit+e2-ABS(cc)**2/(e1    )
               ELSEIF (isubopt(5).EQ.1) THEN !f1
c     ----    each frequency is weighted the same.
                  fit_bart_ind(ifrq)=1-ABS(cc)**2/(e1*e2)
                  fit =fit+ABS(cc)**2/(e1*e2)
               ENDIF
c     WRITE(97,*)e2-ABS(cc)**2/(e1)
            ENDDO               ! nx-loop
         ENDDO                  ! freq loop
c     WRITE(*,*)'isubopt(5)',isubopt(5)
         IF (isubopt(5).EQ.0) THEN
c     ---        each freq is weigthed according to received power
            fit=fit/nfrq/nx
         ELSEIF (isubopt(5).EQ.1) THEN !f1
c     ----    each frequency is weighted the same.
            fit=1-fit/nfrq/nx
         ENDIF
         GOTO 900
      ENDIF                     ! iopt(5)=4
c     
c     
c     
      IF (iopt(5).EQ.5) THEN    !OPTION F
c---  Matched frequency (?); Incoherent addition over DEPTH, PHASES information
               IF ((isubopt(5).EQ.1) .or. (isubopt(5).EQ.2)) THEN              
                  fitlong =0
               ELSEIF ((isubopt(5).EQ.3) .or. (isubopt(5).EQ.4)) THEN
                  fitlong =1
               ENDIF
c               WRITE(*,*)'isubopt(5)',isubopt(5),fitlong
c       pause
       DO jdep=1,ndep
            DO i=1,nx            
               e1=0.
               e2=0.
               cc=0.
               DO ifrq=1,nfrq
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  e1=e1+ABS(resp(i+index))**2 ! resp
                  e2=e2+ABS(DATA(i+index))**2 ! data
                  cc=cc+CONJG(DATA(i+index))*(resp(i+index))
c     WRITE(*,*)'cost,data, jdep,i,ifrq,',jdep,i,ifrq,
c     1         REAL(DATA(i+index)),REAL(resp(i+index)),e1,e2,cc
               ENDDO            ! freq loop
               IF (e1.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: All calculated data are zero '
                  WRITE(*,*)' cost: All calculated data are zero '
                  IF (iopt(31).EQ.0) THEN
                     WRITE(*,*)'...Stopping from cost'
                     STOP
                  ELSE
                     fit =1e10
                     RETURN
                  ENDIF
               ELSEIF (e2.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: observed data are zero'
                  WRITE(*,*)' cost: observed data are zero'
                  WRITE(*,*)'...Stopping from cost'
                  STOP
               ENDIF
               IF ((isubopt(5).EQ.1) .or. (isubopt(5).EQ.3)) THEN              
c     ---        each freq is weigthed according to received power
                 fittemp=e2-ABS(cc)**2/(e1)
                 if (fittemp.le.0) fittemp=0
                  fit_bart_ind(jdep)=fittemp
               ELSEIF ((isubopt(5).EQ.2).or. (isubopt(5).EQ.4)) THEN !f1
c     ----    each frequency is weighted the same.
                  fittemp=1-ABS(cc)**2/(e1*e2)
                  if (fittemp.le.0) fittemp=0
                  fit_bart_ind(jdep)=fittemp
               endif 
c         write(*,*)'fitlong2',fitlong,fit_bart_ind(jdep),ndep,nx
              IF ((isubopt(5).EQ.1) .or. (isubopt(5).EQ.2)) THEN              
                  fitlong =fitlong+fit_bart_ind(jdep)
               ELSEIF ((isubopt(5).EQ.3) .or. (isubopt(5).EQ.4)) THEN
                  fitlong =fitlong*fit_bart_ind(jdep)**(1.0/ndep/nx)
               ENDIF
             ENDDO               ! nx-loop
         ENDDO                  ! ndep--loop
         
         IF ((isubopt(5).EQ.1) .or. (isubopt(5).EQ.2)) THEN              
            fit=fitlong/ndep/nx   ! nfrq
         ELSEIF ((isubopt(5).EQ.3) .or. (isubopt(5).EQ.4)) THEN
            fit=real(fitlong)+1e-11
         ENDIF
c      write(*,*)'fitlong', fitlong,fit,ndep,nx
c      write(*,*) fit_bart_ind(1),fit_bart_ind(128)
c               write(*,*)'fitlong', fitlong
         GOTO 900
      ENDIF                     ! iopt(5)=5
      IF (iopt(5).EQ.6) THEN    ! option j
c---- Incoherent addition over frequencies, PHASES information
c     WRITE(*,*)' Using magnitude and phase'
         DO ifrq=1,nfrq
            DO ix=1,nx            
               e1=0.
               e2=0.
               cc=0.
               DO jdep=1,ndep
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  e1=e1+ABS(resp(ix+index))**2 ! resp
                  e2=e2+ABS(DATA(ix+index))**2 ! data
                  cc=cc+CONJG(DATA(ix+index))*(resp(ix+index))
c     WRITE(99,*)'data,resp',ix,
c     1         (DATA(ix+index)),(resp(ix+index)),
c     1   20*LOG10(ABS(DATA(ix+index))),20*LOG10(ABS(resp(ix+index)))
               ENDDO            ! ndep--loop
               IF (e1.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: All calculated data are zero '
                  WRITE(*,*)' cost: All calculated data are zero '
                  IF (iopt(31).EQ.0) THEN
                     WRITE(*,*)'...Stopping from cost'
                     STOP
                  ELSE
                     fit =1e10
                     RETURN
                  ENDIF
               ELSEIF (e2.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: observed data are zero'
                  WRITE(*,*)' cost: observed data are zero'
                  WRITE(*,*)'...Stopping from cost'
                  STOP
               ENDIF
               fit_bart_ind(ifrq)=e1+e2-2*ABS(cc)
               fit =fit+ fit_bart_ind(ifrq)
c     WRITE(*,*)e1,e2,cc,fit
            ENDDO               ! nx-loop
         ENDDO                  ! freq loop
c     fit=1-fit/nfrq/nx
         GOTO 900
      ENDIF                     ! iopt(5)=6
C     summing over a horizontal array.
      IF (iopt(5).EQ.7) THEN    ! option k
c---- Incoherent addition over frequencies, PHASES information
c     WRITE(*,*)' Using magnitude and phase'
         DO ifrq=1,nfrq
            DO jdep=1,ndep
               e1=0.
               e2=0.
               cc=0.
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO ix=1,nx            
                  e1=e1+ABS(resp(ix+index))**2 ! resp
                  e2=e2+ABS(DATA(ix+index))**2 ! data
                  cc=cc+CONJG(DATA(ix+index))*(resp(ix+index))
c     WRITE(*,*)'cost,data',i,
c     1         (DATA(i+index)),(resp(i+index))
               ENDDO            ! ndep--loop
               IF (e1.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: All calculated data are zero '
                  WRITE(*,*)' cost: All calculated data are zero '
                  IF (iopt(31).EQ.0) THEN
                     WRITE(*,*)'...Stopping from cost'
                     STOP
                  ELSE
                     fit =1e10
                     RETURN
                  ENDIF
               ELSEIF (e2.EQ.0) THEN 
                  WRITE(prtfil,*)' cost: observed data are zero'
                  WRITE(*,*)' cost: observed data are zero'
                  WRITE(*,*)'...Stopping from cost'
                  STOP
               ENDIF
               IF (isubopt(5).EQ.0) THEN
c     ---        each freq is weigthed according to received power
                  fit_bart_ind(ifrq)=e2-ABS(cc)**2/(e1)
                  fit =fit+e2-ABS(cc)**2/(e1    )
               ELSEIF (isubopt(5).EQ.1) THEN
c     ----    each frequency is weighted the same.
                  fit_bart_ind(ifrq)=1-ABS(cc)**2/(e1*e2)
                  fit =fit+ABS(cc)**2/(e1*e2)
               ENDIF
            ENDDO               ! ndep-loop
         ENDDO                  ! freq loop
c     WRITE(*,*)'isubopt(5)',isubopt(5)
         IF (isubopt(5).EQ.0) THEN
c     ---        each freq is weigthed according to received power
            fit=fit/nfrq/ndep
         ELSEIF (isubopt(5).EQ.1) THEN
c     ----    each frequency is weighted the same.
            fit=1-fit/nfrq/ndep
         ENDIF
         GOTO 900
      ENDIF                     ! iopt(5)=7
c     
c     
c---- optimization w/o USE of phase
c     
c     
      DO jdep=1,ndep 
         DO ifrq=1,nfrq
            j=(ifrq-1)*ndep+jdep
            index=(j-1)*nx
            DO i=1,nx
               e1=e1+ABS(resp(i+index))**2 ! resp
               e2=e2+ABS(DATA(i+index))**2 ! data
               e3=e3+ABS(DATA(i+index))*ABS(resp(i+index))
c     WRITE(prtfil,*)' data-loop',(DATA(i+index))**2,e2 ! data
c     WRITE(prtfil,*)' RESP-loop',(RESP(i+index))**2,e1 ! data
            ENDDO
         ENDDO
      ENDDO
      IF (e1.EQ.0) THEN 
         WRITE(prtfil,*)' cost: All calculated data are zero '
         WRITE(*,*)' cost: All calculated data are zero '
         IF (iopt(31).EQ.0) THEN
            DO jdep=1,ndep 
               DO ifrq=1,nfrq
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  DO i=1,nx
                     e2=e2+ABS(DATA(i+index))**2 ! data
                     e1=e1+ABS(resp(i+index))**2 ! resp
                     WRITE(prtfil,*)' RESP-loop',(RESP(i+index))**2,e1 ! data
                     WRITE(prtfil,*)' data',(DATA(i+index))**2,e2 ! data
                  ENDDO
               ENDDO
            ENDDO

            WRITE(*,*)'...Stopping from cost'
            STOP
         ELSE
            fit =1e10
            RETURN
         ENDIF
      ELSEIF (e2.EQ.0) THEN 
         WRITE(prtfil,*)' cost: observed data are zero'
         WRITE(*,*)' cost: observed data are zero'
         WRITE(*,*)'...Stopping from cost'
         STOP
      ENDIF
c     
      amp=SQRT(e2/e1)           ! (data/resp)
      
      IF (iopt(5).EQ.0) THEN    ! This is based on PHASE information...
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
                  cc=cc+CONJG(DATA(i+index))*resp(i+index)
               ENDDO
            ENDDO
         ENDDO
         fit =1-ABS(cc)**2/e1/e2
         
      ELSE IF(iopt(5).EQ.1) THEN !   ! option n  
c--   this is based ONLY on AMPLITUDES
       amp=e3/e1       
c        e2=1.0/SQRT(e2)   ! data
c         e1=1.0/SQRT(e1)   ! resp
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
                  fitlong=fitlong+(REAL(DATA(i+index))
     1                 -REAL(resp(i+index))*amp)**2
c                  fitlong=fitlong+(REAL(DATA(i+index))*e2
c     1                 -REAL(resp(i+index))*e1)**2
c     resp(i+index)=amp*ABS((resp(i+index)))
                  
c      WRITE(*,*)'cost',(i+index),resp(i+index),DATA(i+index),amp,fitlong
               ENDDO
c     WRITE(*,*)'cost tpem', resp(nx+index),nx,amp
            ENDDO
         ENDDO
         fit =fitlong
         
      ELSE IF(iopt(5).EQ.3) THEN ! this is based on AMPLITUDES---   opt X
                                ! WITHOUT correction for offset
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               fitlong=0
               DO i=1,nx
                  fitlong=fitlong
     1                 +(REAL(DATA(i+index))-REAL(resp(i+index)))**2
c     WRITE(*,*)'cost,data',i,
c     1         REAL(DATA(i+index)),REAL(resp(i+index))
               ENDDO
               fit_bart_ind(ifrq)=fitlong
               fit=fit+fitlong
            ENDDO
         ENDDO
         fit =fit/e2
         
      ELSE IF(iopt(5).EQ.8) THEN ! opt Y searching for lag 
c--   this is MBMF  JPH IEEE 1999. 1-R^2, Eq (13)
         e2=1.0/SQRT(e2)
         e1=1.0/SQRT(e1)
         fitlong=1e10
         iup=100
         DO ilag=1,iup
            fitlag=0
           fitlagcmplx=0
            tau=1/(frq(nfrq)-frq(1))*(ilag-1)/iup
            DO i=1,nx
               DO jdep=1,ndep 
                  DO ifrq=1,nfrq
                     j=(ifrq-1)*ndep+jdep
                     index=(j-1)*nx
                     fitlagcmplx=fitlagcmplx+CONJG(DATA(index+i))*e2
     1                     *resp(index+i)*e1
     2                     *cexp(cmplx(0d0,frq(ifrq)*8*atan(1.)*tau))
                  ENDDO
               ENDDO
            ENDDO
            fitlag=1-abs(fitlagcmplx)**2
                  IF (fitlag.LT.fitlong) THEN
                     fitlong=fitlag
                     iminlag=ilag
c     WRITE(*,*)' fitlong,iminlag', fitlong,iminlag,fitlagcmplx
                  ENDIF
         ENDDO
         WRITE(*,*)'iminlag=',iminlag
         fit =fitlong
         
      ELSE IF(iopt(5).EQ.10) THEN ! this is based on AMPLITUDES 
c     WITH uncertainty in both DATA and range
         
         DO j=1,ncurv
            index=(j-1)*nx
            DO i=1,nx
               fit=100000
               DO ix=1,nx
                  fit2=(DATA(i+index)-resp(ix+index))**2/sigmad**2
     &                 +(xranges(i)-xranges(ix))**2/sigmar**2
                  IF (fit2.LE.fit) THEN
                     fit=fit2
                     ilag=ix
                  ENDIF
               ENDDO
               IF (i.NE.ilag)       WRITE(*,*)i,ilag,        
     &              xranges(i),xranges(ilag),fit,fitlong
               fitlong=fitlong+fit
            ENDDO
         ENDDO
         fit =ABS(fitlong)
      ENDIF
c     
c     
c---  Here we INCLUDE regularizing terms in the object FUNCTION
c     
 900  CONTINUE

      CALL getmodelreal(xcurr)
      DO iparm=1,nparm
c         if (fmin(iparm).gt.0) then
c            xl=fmin(iparm)*xonel
c         else
c            xl=fmin(iparm)*xoneu
c         endif
c         if (fmax(iparm).gt.0) then
c            xu=fmax(iparm)*xoneu
c         else
c            xu=fmax(iparm)*xonel
c         endif
            xl= fmin(iparm)-df(iparm)
            xu= fmax(iparm)+df(iparm)

         IF  ((xcurr(iparm).LT.xl).OR.
     1        (xcurr(iparm).GT.xu)) THEN
            WRITE(*,*)'cost: regularizing,iparm,xcurr(iparm):',
     1           iparm,xcurr(iparm),xl,xu
            fit=fit+10
         ENDIF
      ENDDO

CC    cfh 2004,9,24
c      write(*,*)'isubopt(36)',isubopt(36)
      IF (isubopt(36) .EQ. 1) THEN ! optimizing for nu
c     ndep added by chenfen 27 Sep 04
c     ndep is replaced by rankCd + 1 (Jeffrey's prior 1/nu is used)  
c         write(*,*) 'old fit',fit
         fit=fit/numh+(rankCd+1)*LOG(numh)
c         write(*,*) 'new fit',fit
      ENDIF

c  Enumerative integration; fit is divided by numh; 
c        To fix this, here fit is multiplied by numh  
c---pg 23 sep 07 I dont understand this; maybe it should be inside the ifg then condition above 
c      IF  (iopt(4) .EQ. 4 .and. isubopt(4) .EQ. 2) THEN
c         fit = numh*fit
c      ENDIF


c ................................................      
      IF (iopt(9).EQ.1) THEN
         IF (iopt(30).EQ.6) THEN
            mprof=iprof
            IF ((metcon+40).LT.0) THEN
               fit= fit-0.01*(metcon+40)/20
            ELSE
            ENDIF
         ELSE
            CALL getmodelreal(xcurr)
            regsum=0
            DO iparm=1,nparm
               IF (iapri(iparm).GE.0) THEN
                  regsum = regsum+ ABS((xcurr(iparm)-xprior(iparm))/
     1                 (Fmax(iparm)-Fmin(iparm)))*
     1                 ABS((xcurr(iparm)-xprior(iparm))/
     1                 (Fmax(iparm)-Fmin(iparm)))*
     1                 iapri(iparm)
               ENDIF
            ENDDO
            regsum=regsum/nprior
            fit=fit+regsum
         ENDIF
      ENDIF
      


c     IF (.NOT.( ABS(fit).LE.10)) THEN  
c     WRITE(*,*)' '
c     WRITE(*,*)'fit=',fit 
c     WRITE(*,*)' The objective function is larger than 10'
c     WRITE(*,*)' This indicates that the specified environment or'
c     WRITE(*,*)' forward model is wrong....'
c     WRITE(*,*)' The current model will be rejected'  
c     WRITE(*,*)' '
c     fit=10 
c     ENDIF
c     WRITE(*,*)'fit of obj',fit
      END
c     c***************************
      SUBROUTINE cost_max(fit)
c     computes the best obtainable fit for a given DATA
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL fit
      REAL*8 fitlong
      INTEGER i,ix,ifrq,ibart,jdep
      fit=0
c     
c     WRITE(*,*)' Entering cost'
      IF (iopt(13).EQ.1) THEN
         fitlong = 0
c     
c-------BARTLETT estimator
c     
         DO 100 ifrq=1,nfrq
            DO 101 ix=1,nx 
               ibart=ix+((ifrq-1))*nx ! which observation
c>>>> now for each covariance matrix
               fit_bart_ind(ibart)=xcov_trace(ibart)-xcov_lamb(ibart)
               IF (iopt(18).EQ.0) THEN
                  fitlong=fitlong+fit_bart_ind(ibart)
               ELSE IF (iopt(18).EQ.1) THEN
                  fitlong=fitlong+10.*LOG10(fit_bart_ind(ibart))
               ELSE IF  (iopt(18).EQ.2) THEN

               ENDIF
 101        CONTINUE            !nx
 100     CONTINUE               !nfrq 
         fit=fitlong/nbart
         IF (iopt(18).EQ.1) THEN ! Change to dB scale
            fit=10**(fit/10.)
         ENDIF
         
         GOTO 900
      ENDIF                     ! FINISHED WITH COVARIANCE
c     
c---- Here starts the evaluation based on REAL observations
c     WITH PHASE information
c     
      IF (iopt(5).EQ.4) THEN    ! option f
c---- Incoherent addition over frequencies, PHASES information
         DO ifrq=1,nfrq
            DO ix=1,nx            
               fit_bart_ind(ifrq)=0
               fit =fit+1
            ENDDO               ! nx-loop
         ENDDO                  ! freq loop
         fit=1-fit/nfrq/nx
         GOTO 900
      ENDIF                     ! iopt(5)=4
c     
      IF (iopt(5).EQ.5) THEN
c---  Matched frequency (?); Incoherent addition over DEPTH, PHASES information
         DO jdep=1,ndep
            DO i=1,nx            
               fit_bart_ind(jdep)=1-1
               fit =fit+1
            ENDDO               ! nx-loop
         ENDDO                  ! ndep--loop
         fit=1-fit/ndep/nx      ! nfrq
         GOTO 900
      ENDIF                     ! iopt(5)=5
c     
c     
c---- optimization w/o USE of phase
c     
      IF (iopt(5).EQ.0) THEN    ! This is based on PHASE information...
         fit =1-1
      ELSE IF(iopt(5).EQ.1) THEN ! opt n  
c--   this is based ONLY on AMPLITUDES
         fit =0
      ELSE IF(iopt(5).EQ.3) THEN ! this is based on AMPLITUDES---   opt X
                                ! WITHOUT correction for offset
         fit =0
      ENDIF
c     
 900  CONTINUE
      END
c************************************************************
      SUBROUTINE norm()
c     Normalize the response
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER i,j,index,jdep,ifrq
      
      IF (itrans(1).EQ.1) THEN  ! MULTIPLYING with SQRT(R)
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
                  resp(i+index)=resp(i+index)*SQRT(xranges(i))
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (itrans(2).EQ.3) THEN ! DIVIDING by SQRT R'
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
                  resp(i+index)=resp(i+index)/SQRT(xranges(i))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (MOD(itrans(3),10).EQ.2) THEN ! transforming  resp to 'DB SCALE'
c     PRINT *,'DB-loop'
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
c     WRITE(*,*)'index,j,i,resp(i+index)'
c     WRITE(7,*)'dB loop',index,j,i,resp(i+index)
c     1                 -20*LOG10(ABS(resp(i+index)))
                  resp(i+index)=-20*LOG10(ABS(resp(i+index)))
c     WRITE(7,*)'dB loop',index,j,i,resp(i+index)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF ((itrans(4)-10*INT(itrans(4)/10)).EQ.4) THEN 
c     ! Multiplying with weight
c     WRITE(*,*) 'resp is multiplied'
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
c     WRITE(*,*)'jdep,ifrq,i,resp,weight',
c     &     jdep,ifrq,i,resp(i+index),weight(i+index)
                  resp(i+index)=resp(i+index)*weight(i+index)
c     WRITE(98,*)'wei resp',i,jdep,ifrq,resp(i+index),weight(i+index)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (itrans(5).EQ.5) THEN  ! taking magnitude
         DO jdep=1,ndep
            DO i=1,nx            
               DO ifrq=1,nfrq
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  resp(i+index)=ABS(resp(i+index))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      END
c************************************************************
      SUBROUTINE normdata()
c     Normalize the DATA
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER i,j,index,ix,ifrq,ibart,jdep,
     1     index2,index3,indexw1,indexw2,j1start,jtemp
      IF (iopt(13).EQ.1) THEN
c---  For covariance matrix ONLY option M (iopt(3)=4) is allowed
         IF (itrans(3).EQ.4) THEN ! Multiplying with weight
            DO 100 ifrq=1,nfrq
               DO 101 ix=1,nx 
                  ibart=ix+((ifrq-1))*nx ! which observation
                  index3=(ibart-1)*ndep-1 
c>>>> now for each covariance matrix
                  j1start=((ifrq-1))*ndep ! where the response starts
                  DO i=1,ndep
                     jtemp=(ifrq-1)*ndep+i
                     indexw2=(jtemp-1)*nx
                     DO j=1,ndep
                        jtemp=(ifrq-1)*ndep+j
                        indexw1=(jtemp-1)*nx
                        index=ix+(j+j1start-1)*nx
                        index2=i+(j+index3)*ndep 
                        cov(index2)=cov(index2)
     1                       *weight(indexw1+ix)*weight(indexw2+ix)
                     ENDDO
                  ENDDO
 101           CONTINUE         !nx
 100        CONTINUE            !nfrq 
         ELSE
            WRITE(*,*)'>>>> For covariance matrix only option M ' 
            WRITE(*,*)'>>>> is allowed for weighting the data '        
            STOP
         ENDIF      
      ENDIF
      IF (itrans(1).EQ.1) THEN  ! MULTIPLYING with SQRT(R)
         
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
                  DATA(i+index)=DATA(i+index)*SQRT(xranges(i))
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (itrans(2).EQ.3) THEN ! DIVIDING by SQRT R'
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
                  DATA(i+index)=DATA(i+index)/SQRT(xranges(i))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (itrans(3).EQ.12) THEN ! transforming data to 'DB SCALE'
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
c     WRITE(*,*)'index,j,i,resp(i+index)'
c     WRITE(*,*)index,j,i,resp(i+index)
                  DATA(i+index)=-20*LOG10(ABS(DATA(i+index)))
c     WRITE(80,*) i, REAL( DATA(i+index))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
c     pg 8 apr 02: changed iopt(4) to itrans(4) 
      WRITE(*,*)'Itrans from cost',itrans(4)
      IF ((itrans(4)/10).EQ.4) THEN ! Multiplying with weight
         WRITE(*,*) 'data is weighted'
         DO jdep=1,ndep 
            DO ifrq=1,nfrq
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx
                  DATA(i+index)=DATA(i+index)*weight(i+index)
c     WRITE(98,*)'wei data',i,jdep,ifrq,DATA(i+index),weight(i+index)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      IF (itrans(3).EQ.5) THEN  ! taking magnitude
         DO jdep=1,ndep
            DO i=1,nx            
               DO ifrq=1,nfrq
                  j=(ifrq-1)*ndep+jdep
                  index=(j-1)*nx
                  DATA(i+index)=ABS(DATA(i+index))
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      END
c*************************************************************
c     c***************************
      SUBROUTINE costmbeam(fit)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL fit
      REAL*8 fitlong,e1,e2
      REAL amp
      INTEGER i,j,ix,ifrq,ibart,j1start,jdep,index,index2
      INTEGER index3
      COMPLEX*16 cc,cx,cc_diag, cc2,cx2,cc_diag2
c---  variables for regularization
c     REAL regsum,xcurr(mpar)
c     INTEGER iparm,ilag,iminlag,ilow,iup,mprof
c     INTEGER iprof
c     REAL mexcess,metcon,basehvar(1000)
c     COMMON/ mexc/ mexcess,metcon,basehvar,iprof
c-----to DO svd
      INTEGER Msour,info
      PARAMETER(msour=2)
      COMPLEX zx(Mdep,Msour)    ! matrix to find SVD of 
      COMPLEX ssvd(msour+1),esvd(msour),usvd(mdep,msour ),
     1     vsvd(msour,msour)
      COMPLEX work2(mdep),wsvd(msour)             
      INTEGER nsvd              ! number of observations

c     
c     WRITE(*,*)' Entering cost'
      fitlong=0.
      IF (iopt(18).EQ.1) THEN   ! Change to dB scale
         fitlong=1d0
      ENDIF
      fit=0.
      e1=0
      e2=0
      cc=(0.d0,0.d0) 
      IF (iopt(13).EQ.1) THEN
c     
c-------BARTLETT estimator
c     
         DO 100 ifrq=1,nfrq
            DO 101 ix=1,nx 
               ibart=ix+((ifrq-1))*nx ! which observation
               index3=(ibart-1)*ndep-1 
c>>>> now for each covariance matrix
               cc      =(0.d0,0.d0) 
               cc_diag =(0.d0,0.d0) 
               cc2     =(0.d0,0.d0) 
               cc_diag2=(0.d0,0.d0) 
c     e1=0.d0
               j1start =((ifrq-1))*ndep ! where the response starts
c     DO j=j1start,ndep+j1start-1
c     index=ix+(j)*nx
c     e1=e1+(ABS(resp(index)))**2 
c     ENDDO
c     WRITE(*,*)'e1',e1
               xrepl_sum(ibart)=1
               DO j=1,ndep
                  index=ix+(j+j1start-1)*nx
                  zx(j,1)=resp(index) 
                  zx(j,2)=resp1(index)
c     WRITE(79,*) j,REAL(zx(j,1)),imag(zx(j,1)),
c     &                 REAL(zx(j,2)),imag(zx(j,2))
               ENDDO
               nsvd=2 
               CALL csvdc(zx,mdep,ndep,nsvd,ssvd,esvd,usvd,mdep,
     &              vsvd,msour,work2,21,info)
               IF (info.NE.0) THEN
                  WRITE(*,*)'********* info from zsvdc 
     &                 ************',info
                  STOP
               ENDIF
               xrepl_sum(ibart)=e1
               DO i=1,ndep
                  cx =(0.d0,0.d0)
                  cx2=(0.d0,0.d0)
                  DO j=i+1,ndep
                     index =ix+(j+j1start-1)*nx
                     index2=i+(j+index3)*ndep 
                     cx    =cx+cov(index2)*usvd(j,1)     
                     cx2   =cx2+cov(index2)*usvd(j,2)     
c     WRITE(*,*)'cx',cx,cov(index2),resp(index),index2,index
                  ENDDO
                  index   =ix+(i+j1start-1)*nx
                  index2  =i+(i+index3)*ndep 
                  cc_diag =cc_diag+cov(index2)*ABS(usvd(i,1))**2     
                  cc      =cc +CONJG(usvd(i,1))*cx
                  cc_diag2=cc_diag2+cov(index2)*ABS(usvd(i,2))**2     
                  cc2     =cc2+CONJG(usvd(i,2))*cx2
               ENDDO
               cc             =2.*cc +cc_diag
               cc2            =2.*cc2 +cc_diag2
               xcor_sum(ibart)=cc
c     cc=cc/e1         ! normalize with replica
c     WRITE(*,*)'cc,cc2',cc,cc2
               fit_bart_ind(ibart)=xcov_trace(ifrq)-dreal(cc+cc2)
               IF (iopt(18).EQ.0) THEN
                  fitlong=fitlong+fit_bart_ind(ibart)
               ELSE IF (iopt(18).EQ.1) THEN
c     pg 4 feb03     fitlong=fitlong+10.*LOG10(fit_bart_ind(ibart))
                  fitlong=fitlong*(fit_bart_ind(ibart))**(1./nbart)
               ELSE IF  (iopt(18).EQ.2) THEN
                  STOP 'cost: opt(18) not defined'
               ENDIF
               
 101        CONTINUE            !nx
            
 100     CONTINUE               !nfrq 
         
         IF (iopt(18).EQ.1) THEN ! Change to dB scale
c     pg 4 feb03       fit=10**(fit/10.)
            fit= fitlong
         ELSE
            fit=fitlong/nbart
         ENDIF
         
      ENDIF                     ! FINISHED WITH COVARIANCE
C     
c     
      END
