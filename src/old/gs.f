c     Metropolis-Hastings sampler ...
c     Peter Gerstoft 02
      SUBROUTINE mcmc
      USE global
      INCLUDE 'comopt.h'
      REAL ppd1(mpar,mdig), ppd2(mpar,mdig), f1(mpar), f2(mpar)
      REAL bound(mpar,2),deltaf(mpar)
      REAL eps,nye,tot,xdum,epssumsq,epssum  ,tot1,tot2
      DOUBLE PRECISION cov1(mpar,mpar), cov2(mpar,mpar),
     &     total1(mpar),total2(mpar)
      REAL u(mpar,mpar)
      CHARACTER*20 COMFILE
      INTEGER itrn, fork,wait ,retval,status
      REAL ran2,diff1,epsold(10000)  ,dfchi
      EXTERNAL fork,wait
      INTEGER getpid,id
      INTEGER iloop,i,jp,ii,irot,lunfit

      epssum=0.
      epssumsq=0.
      nppd=50
      nye = nugibs

      irot=0
      iforwt=0
c---  initialize 
      IF  (isubopt(4).EQ.2) THEN
         CALL OPFILB(80,ierr)
         CALL enumer
         STOP 'Enumerative integration done'
      ENDIF
c---  
      DO i=1,nparm
         DO jp=1,nppd
            ppd1(i,jp)=0
            ppd2(i,jp)=0
         ENDDO
      ENDDO
      DO i=1,nparm
         DO jp=1,nparm
            cov1(i,jp)=0
            cov2(i,jp)=0
         ENDDO
         total1(i)=0
         total2(i)=0
         deltaf(i)=(fmax(i)- fmin(i))/(nppd-1)
      ENDDO

      CALL getmodelreal(f1)
      CALL getmodelreal(f2)

      CALL scalevec(f1)
      CALL scalevec(f2)
      
      WRITE(*,*)' .... SAGA Metropolis-Hastings sampling module ....'
      WRITE(*,*)' Stop criteria, eps', epsstop
      WRITE(*,*)' Rotation is done for eps', epscov
      WRITE(*,*)
      WRITE(*,*)' first call', f1(1)

      DO iloop = 1,10000
c     population one to nspawn-1 is spawned
         WRITE(*,*)' Loop number', iloop
         itrn=fork()      

         IF (itrn.NE.0)   THEN  ! this is for the parent
c     WRITE(*,*)'parent, child process itrn',ipop,itrn
         ELSE
c     WRITE(*,*)'calling 1st sampling', iloop
            iseed=-1001*iloop +1
            xdum=ran2(iseed)    !
            IF (iloop.EQ.1) THEN
               CALL system('/bin/rm -f samp1.mat') 
               OPEN(UNIT=85,FILE='samp1.mat',FORM='UNFORMATTED',
     &              STATUS='new')
               CALL writematlabhead(85)         
c     stop
            ELSE
               OPEN(UNIT=85,FILE='samp1.mat',FORM='UNFORMATTED', 
     1              STATUS='old',position='Append')
            ENDIF
c     CALL opfilr(85,ioer)
            lunfit=85
            CALL gs(f1,ppd1,nye,cov1,total1,u,irot,bound,lunfit)
            CALL flush(lunfit)
            ID=getpid()
c     WRITE(*,*)'From child ID=', ID
            WRITE(comfile,'(A3,I5.5)') 'dum',ID
            OPEN(UNIT=49,FILE=comfile,FORM='UNFORMATTED',STATUS='NEW') 
c     OPEN(UNIT=49,FILE=comfile,FORM='FORMATTED',STATUS='NEW')
c     WRITE(49,*)(f1(i),i=1,nparm)
            WRITE(49)(f1(i),i=1,nparm)
            WRITE(49)((ppd1(i,jp),i=1,nparm),jp=1,nppd)              
            WRITE(49)((cov1(i,jp),i=1,nparm),jp=1,nparm) 
            WRITE(49)(total1(i),i=1,nparm)
            CLOSE(49)
            STOP 'Child finished'
         ENDIF
         IF  (iloop.EQ.1) THEN
c     CALL opfilw(86,ioer)
            OPEN(UNIT=86,FILE='samp2.mat',FORM='unFORMATTED')
            CALL writematlabhead(86)         
         ELSE
            OPEN(UNIT=86,FILE='samp2.mat',FORM='unFORMATTED',
     1           position='Append')

         ENDIF
c     WRITE(*,*)'calling 2nd sampling', iloop
         iseed=-1001*iloop +2
         xdum=ran2(iseed)       ! 
         lunfit=86
         CALL gs(f2,ppd2,nye,cov2,total2,u,irot,bound,lunfit)
         CLOSE(86)
c     
         retval=wait(status) 
c         WRITE(*,*)'From wait retval=', retval
         WRITE(comfile,'(A3,I5.5)') 'dum',retval
c     OPEN(UNIT=49,FILE=comfile,FORM='FORMATTED')         
         OPEN(UNIT=49,FILE=comfile,FORM='UNFORMATTED')         
c     READ(49,*)(f1(i),i=1,nparm)
         READ(49) (f1(i),i=1,nparm)
         READ(49) ((ppd1(i,jp),i=1,nparm),jp=1,nppd)                 
         READ(49) ((cov1(i,jp),i=1,nparm),jp=1,nparm)                 
         READ(49) (total1(i),i=1,nparm)
         CLOSE(49)
         IF (iopt(6).EQ.0) THEN
            CALL system('/bin/rm -f '//comfile//' &')
         ENDIF
         WRITE(*,*)'iforwt',iforwt

         IF (nparm.GE.2) THEN
            CALL rotconv(cov1,cov2,total1,total2,u,bound,
     &           iforwt,diff1,epscov)
            IF ((diff1.LT.epscov).AND. (irot.EQ.0)) THEN
               irot=1
               WRITE(*,*)'Rotating is OK'
            ENDIF
         ENDIF

         CALL gsconv(ppd1,ppd2,eps) 
         IF (eps.LE. epsstop) THEN 
            GOTO 1000
         ELSE 
            epsold(iloop)=eps
            call convergence(epsstop,epsold,iloop)
         ENDIF
      ENDDO
 1000 CONTINUE

      DO i=1,nparm
         tot=0 
         DO jp=1,nppd
            ppd1(i,jp)=ppd1(i,jp)+ppd2(i,jp)
            tot= tot +ppd1(i,jp)
         ENDDO
         DO jp=1,nppd
            ppd1(i,jp)=ppd1(i,jp)/tot
         ENDDO
      ENDDO
C     WRITE out
      WRITE(40,'(a,30i3,a)')' iopt=[',(iopt(i),i=1,30),'];'
      WRITE(40,'(a,i4,a)') 'nparm=',nparm,'; % number of parameters'
      WRITE(40,'(a,i4,a)') 'nppd=',nppd,'; % number of values/parm '
      WRITE(40,'(a,100g14.6,a)') 'e_stop=' , epsstop,';'
      
      DO ii=1,nparm
         IF ((iopt(14).EQ.1) ) THEN
c            IF (par2phy(ii) .EQ. 9 .and. iopt(30) .ne. 11) THEN ! source range
c               WRITE(40,810)' f_min(',ii,')=',fmin(ii)/1000,';'
c               WRITE(40,810)' f_max(',ii,')=',fmax(ii)/1000,';'
c            ELSE
               WRITE(40,810)' f_min(',ii,')=',fmin(ii),     ';'
               WRITE(40,810)' f_max(',ii,')=',fmax(ii),     ';'
c            ENDIF
            WRITE(40,'(a,i2,a,i3,a)')'par2phy(',ii,')=',par2phy(ii),';'
         ELSE
            WRITE(40,'(a,i2,a,i3,a)')'par2phy(',1 ,')=',par2phy(1 ),';'
            WRITE(40,810)' f_min(',1,')=',fmin(1),     ';'
            WRITE(40,810)' f_max(',1,')=',fmax(1),     ';'
 810        FORMAT(a,i2,a,e15.6,a)
         ENDIF
      ENDDO

      WRITE(40,'(a)')'ppd=['
      DO jp=1,nppd
         WRITE(40,'(100g14.6)') (ppd1(i,jp),i=1,nparm)
      ENDDO
      WRITE(40,'(a)')'];'
      WRITE(40,'(a)')'x=['
      DO jp=1,nppd
         WRITE(40,'(100g14.6)') (fmin(i)+deltaf(i)*(jp-0.5),i=1,nparm)
      ENDDO
      WRITE(40,'(a)')'];'
      DO i=1,nparm
         WRITE(40,'(a,i3,a,i3,a)')'par2phy(',i,')= ',par2phy(i),';'
      ENDDO
      WRITE(*,*)'number of forward modelling calls',iforwt
      WRITE(*,*)' End of Metropolis-Hastings '

      CALL system('/bin/mv -f samp1.mat ${filename}.gs1') 
      CALL system('/bin/mv -f samp2.mat ${filename}.gs2') 
      CLOSE(90) 
      CALL system('/bin/mv -f fort.90 ${filename}.conv') 
      END

c***********************************************************      
      SUBROUTINE  convergence(epsstop,eps,nobs)
c***********************************************************      
      integer nobs,i,n2,n3
      real eps(nobs),x(nobs),y(nobs)
      real epsstop,sxx,sxy,sx,sy,a,b,a1,b1
      if (nobs==1) return
      do i=1,nobs
         x(i)=i
      enddo
      y=eps
      sxx=sum(x*x)
      sxy=sum(y*x)
      sx=sum(x);
      sy=sum(y);
      b=(sxy-sx*sy/nobs)/(sxx-sx*sx/nobs)
      a=(sy-sx*b)/nobs
      if (nobs.gt.10) then
         n3=int(nobs/2)
         n2=nobs-n3+1
         x=x-nobs
         sxx=sum(x(n2:nobs)*x(n2:nobs))
         sxy=sum(y(n2:nobs)*x(n2:nobs))
         sx =sum(x(n2:nobs));
         sy =sum(y(n2:nobs));
         b1 =(sxy-sx*sy/(n3))/(sxx-sx*sx/(n3))
         a1 =(sy-sx*b1)/(n3)
         n2=nobs
      else
         b1 = b
         a1 = a
         n2=0
      endif
      WRITE(*,*)'ESTIMATED additional number of runs',
     &     INT((epsstop-eps(nobs))/b)+1, INT((epsstop-eps(nobs))/b1)+1
      END

c***********************************************************      
      SUBROUTINE gs(fopt,ppd,nye,covar,total,u,irot,bound,lunfit)
c***********************************************************      
      USE global
      INCLUDE 'comopt.h'
      REAL ppd(mpar,mdig),fopt(mpar),nye,u(mpar,mpar)
      DOUBLE PRECISION covar(mpar,mpar),total(mpar)
      INTEGER irot,lunfit
      REAL bound(mpar,2)
      INTEGER i,j,jmov,jtemp,jj,maxtemps,nmov
      REAL temp,tmp(mpar)
      REAL e,ebest,arg,eold 
      REAL fold(mpar),fp(mpar),fpold(mpar)
      REAL fr(mpar),frold(mpar),fbest(mpar)
      REAL fmin0(mpar),fmax0(mpar),df0(mpar),dfmax(mpar)
      REAL pde,decay
      INTEGER ibest,id,iaccept,ireject,ipar,iskip
c     functions
      REAL gswalk,ran2,boxwalk,boxwalksmall
c     fopt the dimesionless  model (which could be rotated)
c     fr   dimensionless in unrotated space
c     fp   physical model used for calling forward model
c     irot=1
c     u(1,1)= - SQRT(2.)/2
c     u(1,2)=   SQRT(2.)/2
c     u(2,1)= - SQRT(2.)/2
c     u(2,2)= - SQRT(2.)/2
c     bound(1,1)=-SQRT(2.)/2
c     bound(1,2)=-1.4
c     bound(2,1)=-0.7
c     bound(2,2)=0
c     
      decay = 1.0
      temp0=1
      temp   = temp0            !*exp(-decay*(jtemp-1)**(1./nparm))
      iaccept=0
      ireject=0
      nmov=nparm 
      maxtemps=niter
      IF (isubopt(4).EQ.0)      maxtemps=niter/nmov

c---  initialize
      CALL checkmodel(fopt)
      CALL scaleinv(fopt,fp)
      CALL setmodelreal(fp)
      CALL forw2
      iforwt=iforwt+1
      CALL cost(e)
      IF (isubopt(36) == 0) THEN
         e=e/nye
      ELSE  ! isubopt(36) = 1 is used for optimizing nu
         e=e
      ENDIF
      If (iopt(6).eq.1) then
         write(81,*)(fp(ipar),ipar=1,nparm),e
      endif
c     WRITE(lunfit,'(i6,20g14.8)')iforwt,(fp(ipar),ipar=1,nparm),e
      EBEST=E
      DO JJ=1,Nparm
         FBEST(JJ)=fopt(JJ)
      ENDDO
      DO ipar=1,nparm
         fr(ipar)= fopt(ipar)   !fr is just temporary
      ENDDO
      IF (irot.EQ.1) THEN
         CALL Rotparm(fr,fopt,u,nparm,mpar,iskip,+1) ! fr ->fopt
c     WRITE(*,*)'fr,foptb',(fr(ipar),fopt(ipar),ipar=1,nparm)
      ENDIF

c     Initialize tempratures
      DO i=1,nparm 
         tmp(i) = 1             ! temp0 !*exp(-decay*(jtemp-1)**(1./nparm))
         fmin0(i)=0
         fmax0(i)=1
         dfmax(i)=1
         df0(i)=0.01
         IF (irot.EQ.1) THEN
            fmin0(i)=bound(i,1)
            fmax0(i)=bound(i,2)
            dfmax(i)= (fmax0(i)-fmin0(i))
            df0(i)= dfmax(i) /1000
         ENDIF
      ENDDO
      
      DO jtemp=1,maxtemps
         DO jmov=1,nmov         ! metropolis       
c     try a discrete random walk for xmod & ymod
c     & get etrial, the NEW error
            eold=e
            DO i=1,nparm
               fold(i)=fopt(i)
               frold(i)=fr(i)
               fpold(i)=fp(i)
            ENDDO
c     WRITE(*,*)'fr,foptb',(fr(ipar),fopt(ipar),ipar=1,nparm)
 22         IF (isubopt(4).EQ.0) THEN ! one move at a time
               IF (isubopt(35).EQ.0) THEN
c     fopt(jmov)=     gswalk(fold(jmov),df0(jmov),
c     1                     fmin0(jmov),fmax0(jmov),tmp(jmov))
                  fopt(jmov)=  boxwalk(fold(jmov),df0(jmov),
     &                 fmin0(jmov),fmax0(jmov),tmp(jmov))
               ELSE
                  fopt(jmov)= boxwalksmall(fold(jmov),df0(jmov),
     &                 fmin0(jmov),fmax0(jmov),tmp(jmov))
               ENDIF
            ELSE                ! all at once
               DO i=1,nparm
                  IF (isubopt(35).EQ.0) THEN
                     fopt(i)=  boxwalk(fold(i),df0(i),
     1                    fmin0(i),fmax0(i),tmp(i))
                  ELSE
                     fopt(i)=  boxwalksmall(fold(i),df0(i),
     1                    fmin0(i),fmax0(i),tmp(i))
                  ENDIF
               ENDDO
            ENDIF
c     WRITE(*,*)'parameters',(fopt(i),i=1,nparm)
c     WRITE(*,*)fopt(1),fold(1),df(1),fmin(1),fmax(1),tmp(1)
            IF (irot.EQ.1) THEN !with rotation
               CALL Rotparm(fr,fopt,u,nparm,mpar,iskip,-1) ! fr<-fopt
               IF (iskip.EQ.1) THEN
c     WRITE(72,*)'move outside bounds',
c     2                 iforwt,(fr(ipar),ipar=1,nparm)
                  GOTO 22
               ENDIF 
            ELSE
               DO ipar=1,nparm
                  fr(ipar)= fopt(ipar)
               ENDDO
            ENDIF
            CALL scaleinv(fr,fp)
            CALL setmodelreal(fp)
            CALL forw2
            CALL cost(e)
            iforwt=iforwt+1
c     .cfh.
c            write(*,*) 'isubopt(36) = ', isubopt(36)
            IF (isubopt(36) .eq. 0) THEN
               e = e/nye
               pde = EXP(-(e-eold)/temp)       
            ELSEIF  (isubopt(36) .eq. 1) THEN 
c     isubopt(36) = 1 is used for optimizing nu
               pde = EXP(-(e-eold)/temp)       
            ELSEIF  (isubopt(36) .eq. 2) THEN
               pde = (eold/e)**rankCd  
            ENDIF

            IF (pde.LT.ran2(1)) THEN
c     move rejected;    revert to old solution
c---  write out
               if (iopt(6)==1) 
     &              WRITE(lunfit+10,'(30f13.3)')
     &              (fp(iparm),iparm=1,nparm),e
               ireject=1+ireject
               e=eold
               DO i=1,nparm
                  fopt(i) = fOLD(i)
                  fp(i) = fpOLD(i)
                  Fr(i) = frOLD(i)
               ENDDO
            ELSE
c---- move is accepted!
               if (iopt(6)==1)
     &              WRITE(lunfit+12,'(30f13.3)')(fp(iparm),
     &              iparm=1,nparm),e
               iaccept=1+iaccept
               IF (isubopt(4).EQ.0) THEN
                  df0(jmov)=MAX(df0(jmov),
     &                 kgrow*ABS(tmp(jmov)))*.99
                  df0(jmov)=MIN(df0(jmov),dfmax(jmov))
c     IF (jmov.EQ.1) WRITE(70,*)iforwt,df0(jmov)
               ELSE
                  DO i=1,nparm
                     df0(i)=MAX(df0(i),kgrow*ABS(tmp(i)))*.99
                     df0(i)=MIN(df0(i),dfmax(i))
                  ENDDO
               ENDIF   
            ENDIF
c     WRITE(lunfit,'(i8,20g14.6)')
c cfh            WRITE(lunfit)  (fp(ipar),ipar=1,nparm),e
c     CALL flush(lunfit)
            
            DO i=1,nparm
               DO j=1,nparm
                  covar(i,j)=
     &                 covar(i,j)+(fr(i)-0.5d0)*(fr(j)-0.5d0)
               ENDDO
               total(i)=total(i)+fr(i) -0.5d0
c     WRITE(*,*)' covar',(covar(i,j),j=1,nparm)
            ENDDO
            
            
c     id=1+(fopt(jmov)-fmin(jmov))/(fmax(jmov)-fmin(jmov))*nppd
            DO i=1,nparm
               id=1+fr(i)*nppd
               ppd(i,id)= ppd(i,id)+1 
            ENDDO
c     WRITE(*,*) id,fr(jmov),irot,jmov,fp(jmov)
            
            IF (E.LT.EBEST) THEN
               WRITE(60,'(1x,I7,I7,G12.5,100f12.2)')
     &              1,iforwt,e,(fopt(jj),jj=1,nparm)
c     WRITE(*,*)'   **** NEW BEST ENERGY ****'
               EBEST=E
               IBEST=jtemp
               DO JJ=1,Nparm
                  FBEST(JJ)=fopt(JJ)
               ENDDO
               WRITE(prtfil,*)jtemp,Ebest
               WRITE(prtfil,*)jtemp,(fbest(jj),jj=1,nparm)
            ENDIF
            
         ENDDO                  ! metropolis step
         WRITE(lunfit)  (fp(ipar),ipar=1,nparm),e
      ENDDO                     ! temperature loop
      
      
      IF (irot.EQ.1) THEN
c---  RETURN the unrotated model
         DO ipar=1,nparm
            fopt(ipar)= fr(ipar)
         ENDDO
      ENDIF
      WRITE(*,*)'accept/reject', iaccept,  ireject
c     WRITE(*,'(a,10e14.4)')'end one gs',(fr(jJ),jj=1,Nparm)
      
      END

c***********************************************************
      SUBROUTINE rotconv(cov1,cov2,total1,total2,u, bound,nsamp,
     &     diff1,epsmax)
c***********************************************************
      USE global
      INCLUDE 'comopt.h'
      INTEGER nsamp
      DOUBLE PRECISION cov1(mpar,mpar), cov2(mpar,mpar),
     &     total1(mpar),total2(mpar)
      REAL cross1(mpar,mpar),cross2(mpar,mpar),mbar1(mpar),mbar2(mpar)
      REAL corr1(mpar,mpar),corr2(mpar,mpar)
      REAL sinv1(mpar),sinv2(mpar),m(mpar),mp(mpar),bound(mpar,2)
      REAL u(mpar,mpar)
      REAL diff1,tot,epsmax,ran2
      INTEGER i, jp, imax,ipar,jpar,iskip,ii
      REAL xx,xlow,xup
c--   estimate covariance matrix
      DO ipar=1,nparm
         DO jpar=1,nparm
            cross1(ipar,jpar)=
     &           (cov1(ipar,jpar)-(total1(ipar)/nsamp*
     &           total1(jpar)))/nsamp
            cross2(ipar,jpar)=
     &           (cov2(ipar,jpar)-(total2(ipar)/nsamp*
     &           total2(jpar)))/nsamp
         ENDDO
      ENDDO
c------
      DO ipar=1,nparm
         sinv1(ipar)=SQRT (cross1(ipar,ipar))
         sinv2(ipar)=SQRT (cross2(ipar,ipar))
         IF  (sinv1(ipar).NE.0) THEN
            sinv1(ipar)=1./ sinv1(ipar)
         ELSE
            sinv1(ipar)=0.
         ENDIF
         IF  (sinv2(ipar).NE.0) THEN
            sinv2(ipar)=1./ sinv2(ipar)
         ELSE
            sinv2(ipar)=0.
         ENDIF
      ENDDO

c--   estimate correlation coefficient matrix
      DO ipar=1,nparm
         DO jpar=1,nparm
            corr1(ipar,jpar)=
     &           cross1(ipar,jpar)*sinv1(ipar)*sinv1(jpar)
            corr2(ipar,jpar)=
     &           cross2(ipar,jpar)*sinv2(ipar)*sinv2(jpar)
         ENDDO
      ENDDO


c---  check for convergence   
      diff1=-1
      DO ipar=1,nparm
         DO jpar=1,nparm
            diff1=MAX(diff1,ABS(corr1(ipar,jpar)-corr2(ipar,jpar)))
         ENDDO
         WRITE(7,*)' cor',
     1        (0.5*(corr1(ipar,jpar)+corr2(ipar,jpar)),jpar=1,nparm)
      ENDDO
      WRITE(7,*) ' correlation eps ', diff1
      WRITE(*,*) ' correlation eps ', diff1
c     IF (diff1.GE.epsmax) THEN
c     STOP ' gs: coverged'
c     RETURN
c     ENDIF

c-----start the rotation
      DO ipar=1,nparm
         DO jpar=1,nparm
            u(ipar,jpar)= 0.5*(cross1(ipar,jpar)+cross2(ipar,jpar))
         ENDDO
      ENDDO

c     find rotation 
      CALL svdcmp(u,nparm,nparm,mpar,mpar,cross2,sinv1)  

      DO ipar=1,nparm
         WRITE(7,*)' EV', (u(ipar,jpar),jpar=1,nparm)
      ENDDO
c     cross1 CONTAINS the svd.     

c     generate sampling bounds
C     DO ipar=1,nparm
C     bound(ipar,1)=10000
C     bound(ipar,2)=-10000
C     ENDDO
C     DO ii=1,10000
C     DO ipar=1,nparm
C     m(ipar)=ran2(1)
C     ENDDO
C     CALL  Rotparm (m,mp,u,nparm,mpar,iskip,+1)
C     DO ipar=1,nparm
C     bound(ipar,1)=MIN(mp(ipar),bound(ipar,1))
C     bound(ipar,2)=MAX(mp(ipar),bound(ipar,2))
C     ENDDO
C     ENDDO

C     DO ipar=1,nparm
C     WRITE(*,*)'Bounds, low up',bound(ipar,1), bound(ipar,2)
C     ENDDO
      
      DO ipar=1,nparm
         xup=0
         xlow=0
         DO jpar=1,nparm
            IF (u(jpar,ipar).GT.0) THEN
               xup=xup+ u(jpar,ipar)
            ELSE
               xlow=xlow+ u(jpar,ipar)
            ENDIF
         ENDDO
         bound(ipar,1)=xlow
         bound(ipar,2)=xup
c     WRITE(*,*)'Bounds2, low up',bound(ipar,1), bound(ipar,2)
      ENDDO
      END

c***********************************************************
      SUBROUTINE Rotparm (m,mp,U,nparm,mpar,iskip,irot)
c***********************************************************
c--   Transform back/forth between physical and rotated parameters.
c--   IF irot=+1 transform  m  -> mp (physical to rotated),
c--   IF irot=-1 transform  mp -> m  (rotated to physical). 

      INTEGER nparm,mpar,iskip,irot
      REAL m(mpar),mp(mpar),U(mpar,mpar)
c     REAL minlim(ndim),maxlim(ndim),maxl(ndim)
      REAL sum
      INTEGER id,jd 
c     pg     maxl=maxlim-minlim
      
      IF (irot .EQ. +1) THEN    ! physical -> rotated   
         DO id=1,nparm
            sum=0.
            DO jd=1,nparm
c     pg            sum=sum+U(jd,id)*((m(jd)-minlim(jd))/maxl(jd))
               sum=sum+U(jd,id)*((m(jd)))
            ENDDO
            mp(id)=sum
         ENDDO
      ENDIF

      IF (irot .EQ. -1) THEN    ! rotated -> physical   
         iskip=0
         DO id=1,nparm
            sum=0.
            DO jd=1,nparm
               sum=sum+U(id,jd)*mp(jd)
            ENDDO
            m(id)=sum           !*maxl(id)+minlim(id)    
            IF ((m(id).LT.0.) .OR. (m(id).GT.1.)) iskip=1
         ENDDO
      ENDIF
      RETURN      
      END
 
c***********************************************************
      SUBROUTINE gsconv(ppd1,ppd2,epsmax)
c***********************************************************      
      USE global
      INCLUDE 'comopt.h'
      REAL ppd1(mpar,mdig), ppd2(mpar,mdig),epsmax,dfchi
      REAL diff1,tot,eps(mpar), tot1,tot2,propconv
      INTEGER i, jp, imax
      epsmax=0
      IF (1 .eq. 2) THEN
         tot1 = 0
         tot2 = 0
         DO jp=1,nppd
            tot1 = tot1 +ppd1(1,jp)
            tot2 = tot2 +ppd2(1,jp)
         ENDDO
C         write(*,*) 'tot1,tot2',tot1,tot2
      DO i=1,nparm
         dfchi = nppd-1
         diff1=0
         propconv = 0
         DO jp=1,nppd
            IF (ppd1(i,jp) .eq. 0 .and. ppd2(i,jp) .eq. 0) THEN
               dfchi = dfchi - 1
            ELSE
               diff1 = diff1 + 
     1           (sqrt(tot2/tot1)*ppd1(i,jp)-
     1              sqrt(tot1/tot2)*ppd2(i,jp))**2
     2              /(ppd1(i,jp) + ppd2(i,jp))
            ENDIF
            WRITE(91,*)jp,ppd1(i,jp), ppd2(i,jp)
         ENDDO
         IF (diff1 .EQ. 0) THEN
            STOP ' gs: gaconv all pdp is zero'
         ENDIF

         WRITE(*,*) 'dfchi, diff1',dfchi,diff1 
         propconv = gammq(0.5*dfchi,0.5*diff1)
         WRITE(*,*) 'dfchi,diff1,tot', dfchi,diff1,propconv 
         eps(i) = propconv
         IF (eps(i) .LE. epsmax) THEN
            epsmax=eps(i)
            imax=i
         ENDIF
      ENDDO
      epsmax = 1 - epsmax

      else

      DO i=1,nparm
         diff1=0
         tot=0
         DO jp=1,nppd
            diff1=diff1+ABS(ppd1(i,jp)-ppd2(i,jp))
            tot=tot+ppd1(i,jp)+ppd2(i,jp)
            WRITE(91,*)jp,ppd1(i,jp), ppd2(i,jp)
         ENDDO
         IF (tot.EQ.0) THEN
            STOP ' gs: gaconv all pdp is zero'
         ENDIF
         eps(i)=(diff1/tot)*2
         IF (eps(i).GE.epsmax) THEN
            epsmax=eps(i)
            imax=i
         ENDIF
      ENDDO

      endif
      WRITE(90,'(i8,i6,20g14.6)')iforwt,imax,(eps(i),i=1,nparm)
      CALL flush(91)
      CALL flush(90)
      WRITE(*,*) ' conv eps   ', epsmax,imax
      END
      
c***********************************************************      
      REAL FUNCTION gswalk(xmod,dx,xmin,xmax,tmp)
      REAL xmod,dx,xmin,xmax,tmp
c---  local variables
      REAL ran2
      REAL arand,ayy,dif,xmod1,yy,pwr
      INTEGER ntry
      
c     generate a NEW state following Ingber ----
      
      ntry = 1
 123  CONTINUE
      arand = ran2(1)
      ayy = 0.0
      dif = arand - 0.5
      IF (dif.LT.0.0) ayy = -1.0
      IF (dif.GE.0.0) ayy = +1.0
      
      pwr = ABS(2*arand-1.)
      yy = ayy*tmp*( (1+1/tmp)**pwr - 1.)
      xmod1 = xmod + yy*(xmax-xmin)
      
      IF (xmod1.LT.xmin.OR.xmod1.GT.xmax) THEN
         ntry = ntry + 1
         IF( ntry.LT.100) THEN
            go to 123
         ELSE
            WRITE(*,*)'ingber ntry xmin,xmax,  xmod, xmod1=',
     1           ntry,xmin,xmax,  xmod, xmod1
         END IF
      END IF
      
      gswalk = xmod1
      RETURN
      END

c***********************************************************      
      REAL FUNCTION boxwalk(xmod,dx,xmin,xmax,tmp)
      REAL xmod,dx,xmin,xmax,tmp
c---  local variables
      REAL ran2,xmod1
c     REAL arand,ayy,dif,xmod1,yy,pwr
      INTEGER ntry

      ntry = 1
 123  CONTINUE

      xmod1 = xmin + ran2(1)*(xmax-xmin)
      IF (xmod1.LT.xmin.OR.xmod1.GT.xmax) THEN
         ntry = ntry + 1
         IF( ntry.LT.100) THEN
            go to 123
         ELSE
            WRITE(*,*)'ntry, xmin, xmax, xmod1 = ',
     1           ntry, xmin, xmax, xmod1
         END IF
      END IF
      boxwalk = xmod1
      RETURN
      END

c***********************************************************           
      REAL FUNCTION boxwalksmall(xmod,dx,xmin,xmax,tmp)
      REAL xmod,dx,xmin,xmax,tmp
c---  local variables
      REAL ran2,xmod1
c     REAL arand,ayy,dif,xmod1,yy,pwr
      INTEGER ntry

      ntry = 1
 123  CONTINUE
      tmp=(ran2(1)-0.5)*dx*2
      xmod1 = xmod + tmp
      IF (xmod1.LT.xmin.OR.xmod1.GT.xmax) THEN
         ntry = ntry + 1
         IF( ntry.LT.100) THEN
            go to 123
         ELSE
            WRITE(*,*)'ntry,xmin,xmax, xmod1,xmod,dx=',
     1           ntry,xmin,xmax, xmod1,xmod,dx
         END IF
      END IF
      boxwalksmall = xmod1
      RETURN
      END

c***********************************************************         
      SUBROUTINE scalevec(xstart)
c***********************************************************      
      USE global
      INCLUDE 'comopt.h'
      REAL xstart(*)
      INTEGER i
      
      DO i=1,nparm
         xstart(i)=(xstart(i)-fmin(i))/(fmax(i)-fmin(i))
      ENDDO
      END

c***********************************************************         
      SUBROUTINE scaleinv(xstart,x)
c***********************************************************      
      USE global
      INCLUDE 'comopt.h'
      REAL xstart(*),x(*)
      INTEGER i
      
      DO i=1,nparm
         x(i)=(xstart(i))*(fmax(i)-fmin(i))+fmin(i)
      ENDDO
      END

c***********************************************************         
      SUBROUTINE checkmodel(x)
c***********************************************************      
      USE global
      INCLUDE 'comopt.h'
      REAL x(*)
      INTEGER i
      
      DO i=1,nparm
         IF (x(i).LT.0) WRITE(*,*)'*** WARNING ***',
     1        'Model below search bounds for parameter',i,
     1        'iforwt=',iforwt
         IF (x(i).GT.1) WRITE(*,*)'*** WARNING ***',
     1        'Model below search bounds for parameter',i,
     1        'iforwt=',iforwt
      ENDDO
      END

c***********************************************************       
      SUBROUTINE writematlabhead(lunmat)
c***********************************************************      
      USE global
      INCLUDE 'comopt.h'
      INTEGER iparm,lunmat
      WRITE(lunmat)nparm
      WRITE(lunmat)nppd
      WRITE(lunmat)(ndigit(iparm),iparm=1,nparm)
      WRITE(lunmat)(fmin(iparm),iparm=1,nparm)
      WRITE(lunmat)(fmax(iparm),iparm=1,nparm)
      WRITE(lunmat)(df(iparm),iparm=1,nparm)
      WRITE(lunmat)(par2phy(iparm),iparm=1,nparm)
      WRITE(lunmat)(iopt(iparm),iparm=1,40)
      END

c***********************************************************      
      SUBROUTINE enumer
c***********************************************************      
      USE global
      INCLUDE 'comopt.h'
      REAL t1,xtt(Mpar),funcSWE,pri
      INTEGER i1,i2,i3,i4,i5, iparm

      Write(*,*)' ... Performing enumerative integration'
c      OPEN(UNIT=80,FILE='enum.mat',FORM='UNFORMATTED')

      call writematlabhead(80)
      pri=0

      IF  (nparm.EQ.1) THEN
c---  1 PARAMETER
         DO i1=1,ndigit(1)
            IF (i1.EQ.2) THEN
               CALL rdtime(t1)
               WRITE(*,320) t1*ndigit(1)
 320           FORMAT(' >>> Estimated enumerative cpu time:',
     &              F12.3,' secs.')
            ENDIF
            xtt(1)=fmin(1)+(i1-1)*df(1)
            CALL setmodelreal(xtt)  
            CALL forw2
            CALL cost(func) 
            iforwt=iforwt+1
            if (iopt(23)==1) then
               pri=prior(xtt)
               pri=max(-1e30,log(pri))
            endif      
            WRITE(80)(xtt(iparm),iparm=1,nparm),(func/nugibs-pri)
            if (iopt(6)==1)    
     &           WRITE(81,*)(xtt(iparm),iparm=1,nparm),
     &           (func/nugibs-pri)
         ENDDO
      ELSEIF (nparm.EQ.2) THEN
c---  2 PARAMETERS
         DO i1=1,ndigit(1)
            IF (i1.EQ.2) THEN
               CALL rdtime(t1)
               WRITE(*,320) t1*ndigit(1)
            ENDIF
            xtt(1)=fmin(1)+(i1-1)*df(1)
            DO i2=1,ndigit(2)
               xtt(2)=fmin(2)+(i2-1)*df(2)
               CALL setmodelreal(xtt)  
               CALL forw2
               CALL cost(func) 
               iforwt=iforwt+1
               if (iopt(23)==1) then
                  pri=prior(xtt)
                  pri=max(-1e30,log(pri))
               endif
c               WRITE(*,*)'nugibs, func',nugibs,func
               WRITE(80)(xtt(iparm),iparm=1,nparm),(func/nugibs-pri)
               if (iopt(6)==1)
     &              WRITE(81,*)(xtt(iparm),iparm=1,nparm),
     &              (func/nugibs-pri)
            ENDDO
         ENDDO
      ELSEIF (nparm.EQ.3) THEN
c---  3 PARAMETERS
         DO i1=1,ndigit(1)
            IF (i1.EQ.2) THEN
               CALL rdtime(t1)
               write(*,*)'classic!'
               WRITE(*,320) t1*ndigit(1)
            ENDIF
            xtt(1)=fmin(1)+(i1-1)*df(1)
            DO i2=1,ndigit(2)
               IF ((i1.EQ.1).and. (i2==2)) THEN
                  CALL rdtime(t1)
                  WRITE(*,320) t1*ndigit(1)*ndigit(2)
               ENDIF
               xtt(2)=fmin(2)+(i2-1)*df(2)
               DO i3=1,ndigit(3)
                  xtt(3)=fmin(3)+(i3-1)*df(3)
                  CALL setmodelreal(xtt) 
                  CALL forw2
                  CALL cost(func) 
                  iforwt=iforwt+1
                  if (iopt(23)==1) then
                     pri=prior(xtt)
                     pri=max(-1e30,log(pri))
                  endif
                  WRITE(80)(xtt(iparm),iparm=1,nparm),(func/nugibs-pri)
                  if (iopt(6)==1)
     1                 WRITE(81,*)(xtt(iparm),iparm=1,nparm),
     2                 (func/nugibs-pri)
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (nparm.EQ.4) THEN
c---  4 PARAMETERS
         DO i1=1,ndigit(1)
            IF (i1.EQ.2) THEN
               CALL rdtime(t1)
               WRITE(*,320) t1*ndigit(1)
            ENDIF
            xtt(1)=fmin(1)+(i1-1)*df(1)
            DO i2=1,ndigit(2)
               xtt(2)=fmin(2)+(i2-1)*df(2)
               DO i3=1,ndigit(3)
                  IF ((i1*i2.EQ.1) .and. (i3==2) ) THEN
                     CALL rdtime(t1)
                     WRITE(*,320) t1*ndigit(1)*ndigit(2)*ndigit(3)
                  ENDIF
                  xtt(3)=fmin(3)+(i3-1)*df(3)
                  DO i4=1,ndigit(4)
                     xtt(4)=fmin(4)+(i4-1)*df(4)
                     CALL setmodelreal(xtt)  
                     CALL forw2
                     CALL cost(func) 
                     iforwt=iforwt+1
                     if (iopt(23)==1) then
                        pri=prior(xtt)
                        pri=max(-1e30,log(pri))
                     endif
                     WRITE(80)(xtt(iparm),iparm=1,nparm),
     .                    (func/nugibs-pri)
                     if (iopt(6)==1)
     1                    WRITE(81,*)(xtt(iparm),iparm=1,nparm),
     2                    (func/nugibs-pri)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ELSEIF (nparm.EQ.5) THEN
c---  5 PARAMETERS
         DO i1=1,ndigit(1)
c     IF (i1.EQ.2) THEN
c     CALL rdtime(t1)
c     WRITE(*,320) t1*ndigit(1)
c     ENDIF
            xtt(1)=fmin(1)+(i1-1)*df(1)
            DO i2=1,ndigit(2)
               xtt(2)=fmin(2)+(i2-1)*df(2)
               DO i3=1,ndigit(3)
                  xtt(3)=fmin(3)+(i3-1)*df(3)
                  DO i4=1,ndigit(4)
                     IF ((i1*i2*i3==1).and. (i4==2) ) THEN
                        CALL rdtime(t1)
                        WRITE(*,320) 
     .                       t1*ndigit(1)*ndigit(2)*ndigit(3)*ndigit(4)
                     ENDIF
                     xtt(4)=fmin(4)+(i4-1)*df(4)
                     DO i5=1,ndigit(5)
                        xtt(5)=fmin(5)+(i5-1)*df(5)
                        CALL setmodelreal(xtt)  
                        CALL forw2
                        CALL cost(func) 
                        iforwt=iforwt+1
                        if (iopt(23)==1) then
                           pri=prior(xtt)
                           pri=max(-1e30,log(pri))
                        endif
                        WRITE(80)(xtt(iparm),iparm=1,nparm),
     .                       (func/nugibs-pri)
                        if (iopt(6)==1)
     1                       WRITE(81,*)(xtt(iparm),iparm=1,nparm),
     2                       (func/nugibs-pri)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ELSE
         WRITE(*,*)'*** Number of parameters too large',
     1        ' to enumerative integration, nparm = ',nparm
         STOP
      ENDIF
      END

      FUNCTION gammq(a,x)
      REAL a,gammq,x
CU    USES gcf,gser
      REAL gammcf,gamser,gln
      if(x.lt.0..or.a.le.0.)pause 'bad arguments in gammq'
      if(x.lt.a+1.)then
        call gser(gamser,a,x,gln)
        gammq=1.-gamser
      else
        call gcf(gammcf,a,x,gln)
        gammq=gammcf
      endif
      return
      END

      SUBROUTINE gcf(gammcf,a,x,gln)
      INTEGER ITMAX
      REAL a,gammcf,gln,x,EPS,FPMIN
      PARAMETER (ITMAX=100,EPS=3.e-7,FPMIN=1.e-30)
CU    USES gammln
      INTEGER i
      REAL an,b,c,d,del,h,gammln
      gln=gammln(a)
      b=x+1.-a
      c=1./FPMIN
      d=1./b
      h=d
      do 11 i=1,ITMAX
        an=-i*(i-a)
        b=b+2.
        d=an*d+b
        if(abs(d).lt.FPMIN)d=FPMIN
        c=b+an/c
        if(abs(c).lt.FPMIN)c=FPMIN
        d=1./d
        del=d*c
        h=h*del
        if(abs(del-1.).lt.EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gcf'
1     gammcf=exp(-x+a*log(x)-gln)*h
      return
      END

      SUBROUTINE gser(gamser,a,x,gln)
      INTEGER ITMAX
      REAL a,gamser,gln,x,EPS
      PARAMETER (ITMAX=100,EPS=3.e-7)
CU    USES gammln
      INTEGER n
      REAL ap,del,sum,gammln
      gln=gammln(a)
      if(x.le.0.)then
        if(x.lt.0.)pause 'x < 0 in gser'
        gamser=0.
        return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,ITMAX
        ap=ap+1.
        del=del*x/ap
        sum=sum+del
        if(abs(del).lt.abs(sum)*EPS)goto 1
11    continue
      pause 'a too large, ITMAX too small in gser'
1     gamser=sum*exp(-x+a*log(x)-gln)
      return
      END

      FUNCTION gammln(xx)
      REAL gammln,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,
     *24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
     *-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gammln=tmp+log(stp*ser/x)
      return
      END
