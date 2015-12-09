ccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     This file CONTAINS the subprograms for the genetic algorithm
c     optimization PROGRAM
c     
c     Peter Gerstoft
c     Saclant Undersea Research Centre
c     
c*************************************************
c     
      SUBROUTINE modeldecode(xstart,iq,ierr)
c     The model vector given in xstart is 'decoded' into a INTEGER POINTER 

      USE global
      INCLUDE 'comopt.h'
      REAL xstart(*)
      INTEGER iq,i,ierr
 99   FORMAT(' >>>Warning, The start parm is not in the search int',
     &     ' For parameter',i3,' with value',e12.4)
c     
      ierr=0
      DO i=1,nparm
         model(i,iq)=(xstart(i)-fmin(i))/df(i)+0.5
c     IF (model(i,iq).LT.0 .OR. model(i,iq).GT.ndigit(i)) THEN
         IF (((xstart(i)-fmin(i))*(xstart(i)-fmax(i))).GT.0) THEN
            WRITE(prtfil,99) i,xstart(i)
            WRITE(*,99)      i,xstart(i)
            ierr=1
         ENDIF
      ENDDO
      END
c     
c*****************************************************
c     
      SUBROUTINE header(iun)
c     The model vector given in xstart is 'decoded' into a INTEGER POINTER 

      USE global
      INCLUDE 'comopt.h'
      INTEGER iparm,iun
      DO iparm=1,nparm 
c     This does not work for g77/linux   WRITE(iun,'(a,a,f)')'! ',phystxt2(par2phy(iparm)),xstar(iparm)
c     WRITE(iun,'(a2,a,f)')'! ',phystxt2(par2phy(iparm)),xstar(iparm)
         WRITE(iun,99)'! ',phystxt2(par2phy(iparm)),xstar(iparm)
 99      FORMAT(a2,a40,f15.4)
      ENDDO
      END
c     
c******************************************
c     
      SUBROUTINE readinputstart
c     
c     read the optimization part of the input file
      USE global
      INCLUDE 'comopt.h'
      INTEGER i
      CHARACTER*80 dumch
c      CHARACTER  optread(mopt)
      EQUIVALENCE (dumch,optinv)

 200  FORMAT(A80)
      READ(1,200)title            
      WRITE(prtfil,200)title
      
      READ(1,'(a)')  dumch      !optread
      WRITE(prtfil,*)optinv
c     
c---  decode the option string
c     
      CALL gaoptions
      WRITE(prtfil,'(''                   '',a)')
     1     ' 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9',
     1     ' 0 1 2 3 4 5 6 7 8 9 0'
      WRITE(prtfil,'('' Options:          '',40i2)')(iopt(i),i=1,mopt)
      WRITE(prtfil,'('' Options (norm):   '',40i2)')
     1     (itrans(i),i=1,mopt)
      WRITE(prtfil,'('' Options (isubopt):'',40i2)')
     1     (isubopt(i),i=1,mopt)

      IF ((iopt(4).EQ.1).OR.(iopt(4).EQ.3)) THEN
c-----for simulated annealing
         READ(1,*,err=30) niter,q ,npop,temp0,mter
         GOTO 50
 30      WRITE(*,*)' **** For simulated annealing supply 5 parameters'
         WRITE(*,*)' **** in the niter line'
         WRITE(*,*)' **** the aditional parameters are T0 and Mter'
         WRITE(*,*)' **** usual values are 1, 10000'
         STOP 'Error'
      ELSEIF (iopt(24).EQ.1) THEN
c-----for the initial population
         READ(1,*,err=40)niter,q ,npop, temp0
         GOTO 50
 40      WRITE(*,*)'  **** For option x supply 4 parameters'
         WRITE(*,*)'  **** in the niter line'
         STOP 'Error'
      ELSE
         READ(1,*)niter,q ,npop
      ENDIF
      qin=q
      npopin=npop
 50   CONTINUE
      WRITE(prtfil,*)' Number of forward modelling runs:',npop*niter
      IF ((niter.LT.0) .OR. (npop.LT.0)) THEN
         WRITE(*,*)'wrong input',niter,npop
         STOP
      ENDIF
      IF (q.GT.Mq) THEN
         WRITE(*,*)'q > mq'
         STOP
      ENDIF
      IF (q.GT.Mq_post) THEN
         WRITE(*,*)'q > mq_post'
         STOP
      ENDIF
      IF (MOD(q,2).EQ.1) THEN
         WRITE(*,*) 'q must be even, q=',q
         STOP
      ENDIF
c     
c     READ next BLOCK
c     
      IF (iopt(11).EQ.1) THEN  
         READ(1,*,err=60)px,pu,pm,snr_db
         GOTO 70
 60      WRITE(*,*)' The SNR was not supplied; default SNR assumed.'
         snr_db=40
 70      CONTINUE
         WRITE(prtfil,*)' Gaussian noise added to',
     &        'the data in the *.obs file,  SNR=',snr_db
      ELSEIF (iopt(5).EQ.10) THEN
         READ(1,*,err=61)px,pu,pm, sigmar,sigmad
         WRITE(prtfil,*)' nomalizing standard deviations for r,y  ',
     &        'sigma_r, sigma_d',sigmar, sigmad
         WRITE(*,*)' nomalizing standard deviations for r,y  ',
     &        'sigma_r, sigma_d',sigmar, sigmad
         
         
       ELSE
         READ(1,*)px,pu,pm
       ENDIF
      IF (pm.LT.0 .OR. pm.GT.1) THEN
         WRITE(*,*)'0< pm < 1 !!, pm=',pm
         STOP
      ENDIF
      IF (pu.LT.0 .OR. pu.GT.1) THEN
         WRITE(*,*)'0< pu < 1 !!, pm=',pu
         STOP
      ENDIF
      IF (px.LT.0 .OR. px.GT.1) THEN
         WRITE(*,*)'0< px < 1 !!, px=',px
         STOP
      ENDIF
c---  Metropolis-Hastings sampling 
      IF (iopt(4).EQ.4) THEN
         WRITE(*,*)' Reading input line for Metropolis-Hastings ',
     &        'sampling'
c     .cfh.
         READ(1,*)numh,epsstop,epscov,kgrow,rankCd 
         WRITE(*,*)'numh,epsstop,epscov,kgrow,rankCd'
         WRITE(*,*) numh,epsstop,epscov,kgrow,rankCd 
      ENDIF
       write(*,*)' Read system parameters'
      RETURN
 61   STOP 'sigmax and sigmad was not supplied'
      END
c     
c***************************************************************
c     Reads the range BLOCK
c**** 
      SUBROUTINE read_range_inp(ierr,nwave)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER i,ierr,nwave

      WRITE(prtfil,*)' reading ranges in....'
      IF (iopt(2).EQ.1 .OR. iopt(2).EQ.3 .OR. iopt(2).EQ.4
     1     .OR.iopt(13).EQ.1) THEN ! for using readdat2
 201     READ(1,*)nx
         IF(nx .GT. 0) THEN
            READ(1,*,err=203)(xranges(i),i=1,nx)
         ELSE
            nx=ABS(nx)
            READ(1,*,err=203)rng1,rng2
            IF (rng1.GT.rng2) STOP 'rng1> rng2'
            drng=(rng2-rng1)/(nx-1)
            DO i=1,nx
               xranges(i)=rng1+drng*(i-1)
            ENDDO
         END IF
         WRITE(prtfil,*) 'number of ranges',nx
      ELSE                      ! for using readdata and for syntetic
 202     READ(1,*,err=200) rng1,rng2
         WRITE(prtfil,*)'rng1,rng2', rng1,rng2
         IF (rng1.GT.rng2) STOP 'rng1> rng2'
         nx=MIN(100,mx/nfrq)
         IF ((iopt(2).EQ.0) .OR.(iopt(2).EQ.2)) THEN ! for synthetic data
            drng=(rng2-rng1)/nx
c     nx=(rng2-rng1)/drng+1
c     nx=MIN(nx,mx/nfrq)
            DO i=1,nx
               xranges(i)=rng1+drng*(i-1)
            ENDDO
         ENDIF
      ENDIF
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
c---- for wavenumber domain
      IF (iopt(1).EQ.1) THEN
         nx=nwave
      ENDIF

      ncurv=nfrq*ndep           ! iopt(5)=1,3,7
      IF ((iopt(5).EQ.4) .OR. (iopt(5).EQ.6).OR. (iopt(13).EQ.1)) THEN
         ncurv=nfrq*nx
      ENDIF
      IF (iopt(13).GT.0) THEN   ! for covariance matrix
         nbart=nfrq*nx
         WRITE(prtfil,*)'number of bartlett estimators',nbart
         IF (nbart.GT.mfreq) THEN
            WRITE(*,*)'nbart > mfreq'
            ierr=1
         ENDIF
         ncov_siz=nbart*ndep*ndep
         WRITE(*,*)' Elements in cov-matrix:',ncov_siz
         IF (ncov_siz.GT.mobs) THEN
            WRITE(*,*)'ncov_siz > mobs,  mobs=',mobs
            ierr=1
         ENDIF
      ELSEIF (iopt(5).EQ.5) THEN
         ncurv=nx*ndep          ! broadband matched filter
      ENDIF

      WRITE(prtfil,*)' ncurv=',ncurv,'nbart=',nbart
      IF ((nx.EQ.1).AND.(ndep.EQ.1).AND.(nfrq.EQ.1)) THEN
         WRITE(*,*) ' nx, ndep and nfrq must not all be 1'
         ierr=1
      ENDIF
      IF ((iopt(5).EQ.7).AND. (nx.EQ.1)) THEN
         WRITE(*,*)' *****   Error   ******'
         WRITE(*,*)' Only one phone specified in range '
         WRITE(*,*)' at least two phones required for this inversion'
         ierr=1
      ENDIF
      IF ((iopt(5).EQ.4).AND. (ndep.EQ.1)) THEN
         WRITE(*,*)' *****   Error   ******'
         WRITE(*,*)' Only one phone specified in depth'
         WRITE(*,*)' at least two phones required for this inversion'
         ierr=1
      ENDIF
      IF ((iopt(5).EQ.5).AND. (nfrq.EQ.1)) THEN
         WRITE(*,*)' *****   Error   ******'
         WRITE(*,*)' Only one frequency specified'
         WRITE(*,*)'at least two frequencies required for this',
     &        ' inversion'
         ierr=1
      ENDIF
      IF (nx*ndep.GT.mobs) THEN
         WRITE(*,*) ' nx*ndep >mobs, nx,ndep,mobs',nx,ndep,mobs
         ierr=1
      ENDIF
      END
c****************************************************************
c     Reads the parameters to be optimized
c     
      SUBROUTINE readoptparm
      USE global
      INCLUDE 'comopt.h'
      REAL delf
      INTEGER i,j,ierr
      WRITE(*,*)' '
      WRITE(*,*)'  Reading search intervals...'
      READ(1,*) nparm           ! parameters to be optimized 
      ierr=0
      IF (nparm.GT.Mpar)THEN
         WRITE(*,*)'nparm > mpar'
         WRITE(*,*)' nparm,mpar=',nparm,mpar
         ierr=1
      ENDIF
c     
c---  READ each PARAMETER
c     
      nprior=0
      DO i=1,nparm
         IF (iopt(12).EQ.0) THEN
            READ (1,*) par2phy(i),par2lay(i), fmin(i),fmax(i),ndigit(i)
            IF ((iopt(9).EQ.1) .or. (iopt(23).EQ.1)) THEN
c  03/18/2005
               BACKSPACE(1)
               READ (1,*) par2phy(i),par2lay(i),
     &              fmin(i),fmax(i),ndigit(i),iapri(i)
               WRITE(*,*)' iapri(i)',i,iapri(i)
               nprior=nprior+1
            ENDIF
         ELSE
            READ (1,*) par2phy(i),par2lay(i), par3(i),
     1           fmin(i),fmax(i),ndigit(i)
            IF ((iopt(9).EQ.1) .or. (iopt(23).EQ.1)) THEN
c  03/18/2005
               BACKSPACE(1)
               READ (1,*) par2phy(i),par2lay(i), 
     &              fmin(i),fmax(i),ndigit(i),iapri(i)
               WRITE(*,*)' iapri(i)',i,iapri(i)
               nprior=nprior+1
            ENDIF
         ENDIF
         nbit(i)=NINT(LOG10(1.*ndigit(i))/LOG10(2.)+0.4999)
         IF (ndigit(i).GT.mdig) THEN
            IF (fmin(i).EQ.1)  THEN
               WRITE(*,*)' Probably 3 indices was specified '
               WRITE(*,*)' for the parameter; and only ',
     &              'two were required'
            ENDIF
            WRITE(*,*)' Problems for parameter No.',i
            WRITE(*,*)' Number of requested discretization',  ndigit(i)
            WRITE(*,*)' But maximum number of discretizations are',mdig
            ierr=1
         ENDIF
         IF (par2phy(i).GT.mphys) THEN
            PRINT*,'par2phy(i)>mphys;par2phy(i),mphys=',par2phy(i),mphys 
            ierr=1
         ENDIF
c     
c>>>  control 2 parameters WITH same par2phy have different par2lay
c     
         DO j=1,i-1
            IF (par2phy(i).EQ.par2phy(j)) THEN
               IF (par2lay(i).EQ.par2lay(j)) THEN
                  IF (iopt(12).EQ.0) THEN
                     WRITE(prtfil,*)' 2 param. points to',
     &                    ' same variable:',i,j
                     WRITE(*,*)' 2 param. points to same variable:',i,j
                     STOP
                  ELSEIF (par3(i).EQ.par3(j)) THEN
                     WRITE(prtfil,*)' 2 param. points to',
     &                    ' same variable:',i,j
                     WRITE(*,*)' 2 parameters points to ',
     &                    'same variable:',i,j
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDDO                  ! loop j
      ENDDO
c---  generate the discrete values
      WRITE(prtfil,*)
      WRITE(prtfil,*)' parameters to be optimized:'
      WRITE(prtfil,*)'                  physparm     layer  min.val  ',
     &     'max.val   del.val   nbit  ndig'
      DO i=1,nparm
         IF (ndigit(i).GT.1) THEN
            delf=(fmax(i)-fmin(i))/(ndigit(i)-1)
         ELSE
            WRITE(*,*) 'There is only one step for parameter',i
            delf=fmax(i)-fmin(i)
         ENDIF
         IF (delf.EQ.0) THEN 
            WRITE(*,*)' Stopping for parameter number',i
            WRITE(*,*)'Min and Max bounds are identical'
            WRITE(prtfil,*)' Stopping for parameter number',i
            WRITE(prtfil,*)'Min and Max bounds are identical'
            STOP
         ENDIF
         df(i)=delf
         ilin(i)=0
         IF (iopt(12).EQ.0) THEN ! two index for optimization variable
            WRITE(prtfil,'(i3,x,a22,2i4,3f10.2,2i5)')
     &           i,phystxt2(par2phy(i)),par2phy(i),
     &           par2lay(i), fmin(i),fmax(i),delf,nbit(i),ndigit(i)
         ELSE                   ! tree indexes
            WRITE(prtfil,'(i3,x,a22,3i4,3f10.2,2i5)')
     &           i,phystxt2(par2phy(i)),par2phy(i),par2lay(i),par3(I),
     &           fmin(i),fmax(i),delf,nbit(i),ndigit(i)
         ENDIF
c     DO j=0,ndigit(i)-1
c     fval(j,i)=fmin(i)+delf*j
c     ENDDO
      ENDDO

      IF (ierr.EQ.1) STOP 'from readoptparm'
      END
c***********************
      REAL FUNCTION fval(j,i)
      USE global
      INCLUDE 'comopt.h'
      INTEGER j,i
      fval=fmin(i)+df(i)*j
      END

c****************************************************************
      SUBROUTINE eofinit
c     initialize the EOF. The SUBROUTINE is called from input. 
c     
      USE global
      INCLUDE 'comopt.h'
      INTEGER i,j,k,newblockflag
      INTEGER n1,n2,n1sum,n2sum,neofsec,ieor,jstart
      CHARACTER*250 dumch
c--   This SUBROUTINE initializes the eof and reads them from file 9.
c--   
c---  openfile
      CALL opfilr(9,ieor)
 1    READ(9,'(a80)',END=2)dumch
      IF (dumch(1:1).NE.'!') THEN
         READ(dumch,*,err=99) Neof, Neofvar, Neofsec
      ELSE 
         GOTO 1
      ENDIF
      
      IF (Neof.GT.Meof)THEN
         WRITE(*,*)' Neof, Meof= ',Neof, Meof
         WRITE(*,*)' ***Stopping: Neof > Meof; '
         STOP
      ENDIF
      IF (Neofvar.GT.Meofvar)THEN
         WRITE(*,*)' *** Neofvar > Meofvar'
         STOP
      ENDIF
      
      DO i=1,Neofvar
 3       READ(9,'(a80)',END=2)dumch
         IF (dumch(1:1).NE.'!') THEN
            IF (iopt(12).EQ.0) THEN
               READ(dumch,*,err=99)par2phy_eof(I),par2lay_eof(i)
            ELSE
               READ(dumch,*,err=99)par2phy_eof(I),par2lay_eof(i),
     &              par3_eof(i)
            ENDIF
            WRITE(*,'(a,3i8)')
     &           ' Read EOF point',i,par2phy_eof(I),par2lay_eof(i)
         ELSE 
c     IF (i.EQ.1) THEN
            GOTO 3
c     ELSE
c     WRITE(*,*)' *********EOF file error**********'
c     WRITE(*,*)'Not all EOF points is read in'
c     WRITE(*,*)'actual line, NEOFvar',i,Neofvar
c     WRITE(*,*)' last line read:'
c     WRITE(*,*)dumch
c     ENDIF
         ENDIF
      ENDDO
      WRITE(*,*)Neofvar,' eof points read in'
      DO i=1,Neofvar
         DO j=1,Neof
            eofcoef(i,j)=0
         ENDDO
      ENDDO
      n1sum=0
      n2sum=0
      DO k=1,Neofsec
         newblockflag=0
 4       READ(9,'(a80)',END=2)dumch
c     WRITE(*,*)'pg',dumch
         IF (dumch(1:1).NE.'!') THEN
            READ(dumch,*,err=99)n2,n1 ! n eof, nvar in block
         ELSE 
            newblockflag=1
            GOTO 4
         ENDIF
         IF (newblockflag.NE.1) THEN
            WRITE(*,*)' ****** EOF file error **********'
            WRITE(*,*)' error reading block',k
            WRITE(*,*)' first line should begin with  "!"'
            WRITE(*,*)dumch
         ENDIF
         DO i=1,n1
            READ(9,*,err=199)(eofcoef(i+n1sum,j+n2sum),j=1,N2)
            IF (iopt(6).EQ.1) THEN
               WRITE(prtfil,*)k,i,(eofcoef(i+n1sum,j+n2sum),j=1,N2)
            ENDIF
         ENDDO
         n1sum=n1sum+n1
         n2sum=n2sum+n2
         WRITE(*,*)' Block',k,'is read in. Total blocks:',Neofsec
      ENDDO
      IF (n2sum.NE.neof) THEN
         WRITE(*,*)' Neof ne sum of neofs in blocks'
         WRITE(*,*)'neof,n2sum=',neof,n2sum
         STOP
      ENDIF
      IF (n1sum.NE.neofvar) THEN
         WRITE(*,*)' Neofvar ne sum of neofvar in blocks'
         WRITE(*,*)'neofvar,n1sum=',neofvar,n1sum
         STOP
      ENDIF

      WRITE(prtfil,'(/,a)')' The EOF are'
      WRITE(prtfil,*)'      physparm layer  ... shape functions'
      DO i=1,Neofvar
         IF (iopt(12).EQ.0) THEN
            WRITE(prtfil,'(x,i3,x,2i3,2x,100f5.1)')
     &           i,par2phy_eof(i),
     &           par2lay_eof(i),(eofcoef(i,j),j=1,Neof)
         ELSE
            WRITE(prtfil,'(x,i3,x,3i3,2x,100f6.2)')
     &           i,par2phy_eof(i),
     &           par2lay_eof(i),par3_eof(i),
     &           (eofcoef(i,j),j=1,Neof)
         ENDIF
c>>>  control 2 parameters WITH same par2phy have different par2lay
         DO j=1,i-1
            IF (par2phy_eof(i).EQ.par2phy_eof(j)) THEN
               IF (par2lay_eof(i).EQ.par2lay_eof(j)) THEN
                  IF (iopt(12).EQ.0) THEN
                     WRITE(prtfil,*)' 2 param-eof points to',
     &                    ' same variable:',i,j
                     WRITE(*,*)' 2 param. eof points to',
     &                    ' same variable:',i,j
                     STOP
                  ELSEIF (par3_eof(i).EQ.par3_eof(j)) THEN
                     WRITE(prtfil,*)' 2 param-eof points to',
     &                    ' same variable:',i,j
                     WRITE(*,*)' 2 parameter eof points to',
     &                    ' same variable:',i,j
                     STOP
                  ENDIF
               ENDIF
            ENDIF
         ENDDO                  ! loop j
      ENDDO
c     
c---  starting coefficients
      WRITE(*,*)' Reading starting coefficints....'
      jstart=1
      GOTO 5
 50   jstart=j
 5    READ(9,'(a)',END=2)dumch
      IF (dumch(1:1).NE.'!') THEN
         BACKSPACE(9)
         READ(9,*,err=50)(aeof(j),j=jstart,neof)     
c     READ(dumch,'(200g)',err=50)(aeof(j),j=jstart,neof)     
      ELSE 
         GOTO 5
      ENDIF
      WRITE(prtfil,*)' The starting values of the EOF amplitudes are:' 
      WRITE(prtfil,'(10f10.2)')(aeof(j),j=1,neof)
c--   limits on each PARAMETER in the EOF
      nlimits=0
c     READ(9,*)nlimits
c     WRITE(prtfil,*)' There are limits on the following points'
c     WRITE(prtfil,*)' Point  fmin  fmax'
c     DO i=1,Nlimits
c     READ(9,*)point_eof(I),fmin_lim(i),fmax_lim(i)
c     WRITE(prtfil,'(i5,2f10.2)')point_eof(I),fmin_lim(i),fmax_lim(i)
c     ENDDO 

      WRITE(prtfil,*)
      CLOSE(9)
      RETURN
 99   WRITE(*,*)' Problems reading EOF file, last line read:'
      WRITE(*,*)dumch
      STOP
 199  WRITE(*,*)' Problems reading EOF file'
      WRITE(*,*)' read ',i,'lines of the EOF-coefficient file'
      
 2    STOP ' Reached End of file of EOF-file'
      END
c************************************************
c     Computes the eof value and check IF the results is within the limits
      SUBROUTINE eofval
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL xeof(500)
      INTEGER i
      CALL  eofvalpoint(xeof)
      DO i=1,nlimits
         IF ((xeof(point_eof(i)).LT.fmin_lim(i)) .OR.
     &        (xeof(point_eof(i)).GT.fmax_lim(i))) THEN
            lwascomp=-1
            RETURN
         ENDIF
      ENDDO
      CALL setmodelEOF(xeof)      
      END
c*******************************************
      
      SUBROUTINE eofval_nocheck
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL xeof(100)
      CALL  eofvalpoint(xeof)
c     DO i=1,nlimits
c     IF ((xeof(point_eof(i)).LT.fmin_lim(i)) .OR.
c     &      (xeof(point_eof(i)).GT.fmax_lim(i))) THEN
c     WRITE(*,*)' EOF funtion',i,xeof(point_eof(i))
c     lwascomp=-1
c     RETURN
c     ENDIF
c     ENDDO
      CALL setmodelEOF(xeof)      
      RETURN
      END
c*******************************************
      SUBROUTINE eofvalpoint(xeof)
c---  computes the parameters corresponding to each of the EOF functions
      USE global
      INCLUDE 'comopt.h'
      REAL xeof(*)
      INTEGER i,j
      DO i=1,Neofvar
         xeof(i)=0.
         DO j=1,neof
            xeof(i)=xeof(i)+aeof(j)*eofcoef(i,j)
         ENDDO
      ENDDO
      END
c******************************************
c*****************************************
c     insert the initial model into the model vector
      SUBROUTINE initga
      USE global
      INCLUDE 'comopt.h'
      INTEGER jn,iq
      REAL ran2
c---  initial models
      DO jn=1,nparm
         DO iq=1,q
            model(jn,iq)=INT(ran2(1)*ndigit(jn))
         ENDDO
      ENDDO          
c     copy all models into one matrix to keep statistics of all attempted models.
c     DO jn=1,nparm
c     DO iq=1,q
c     allmodel(jn,iq)=model(jn,iq)
c     ENDDO
c     ENDDO
c     DO iq=1,q
c     allmodel(nparm+1,iq)=iq
c     allmodel(nparm+2,iq)=ipop
c     allfit(iq)=fit(iq)
c     ENDDO
      
c     iallmodel=iallmodel+q               !counter        
      

c     WRITE(prtfil,*)' initial model'
c     DO jn=1,nparm
c     WRITE(prtfil,'(26i3)')(model(jn,iq),iq=1,q)
c     WRITE(prtfil,*)(fval(model(jn,iq),jn),iq=1,q)   
c     ENDDO   
      END
c*******************************************************
      SUBROUTINE costbart
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL fit,e1,e2
      DOUBLE PRECISION  fitlong
      INTEGER i,j,ix,ifrq,ibart,j1start,jdep,index,index2
      INTEGER index3
      COMPLEX*16 cc,cx


c     
      fitlong=0.
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
               cc=(0.d0,0.d0) 
               e1=0.
               j1start=((ifrq-1))*ndep ! where the response starts
               DO j=j1start,ndep+j1start-1
                  index=ix+(j)*nx
                  e1=e1+REAL(resp(index))**2+AIMAG(resp(index))**2 
               ENDDO

               DO i=1,ndep
                  cx=(0.d0,0.d0)
                  DO j=1,ndep
                     index=ix+(j+j1start-1)*nx
                     index2=i+(j+index3)*ndep 
                     cx=cx+cov(index2)*resp(index)     
c     WRITE(*,*)'cx',cx,cov(index2),resp(index),index2,index
                  ENDDO
                  cc=cc+CONJG(resp(ix+(i+j1start-1)*nx))*cx
               ENDDO
               cc=cc/e1         ! normalize with replica
               fit=(1.0-cc)
               WRITE(prtfil,'(a,i4,f9.2)')
     &              ' Bartlett power (dB) for freq',
     &              ifrq,10*LOG10(1-fit) 
 101        CONTINUE            !nx
 100     CONTINUE               !nfrq 
         GOTO 900
      ENDIF                     ! FINISHED WITH COVARIANCE
C     
c     
c---- Here starts the evaluation based on REAL observations
c     
c     
      IF (iopt(5).EQ.4) THEN
c---- Icoherent addition over frequencies, PHASES information
         DO ifrq=1,nfrq
            e1=0.
            e2=0.
            cc=0.
            DO jdep=1,ndep
               j=(ifrq-1)*ndep+jdep
               index=(j-1)*nx
               DO i=1,nx            
                  e1=e1+ABS(resp(i+index))**2 ! resp
                  e2=e2+ABS(DATA(i+index))**2 ! data
                  cc=cc+CONJG(DATA(i+index))*(resp(i+index))
               ENDDO            ! nx-loop
            ENDDO               ! ndep--loop
            IF ((e1.EQ.0).OR.  (e1.EQ.0)) THEN
               WRITE(*,*)' cost: calculated or observed data is '
               WRITE(*,*)'       all zero'
               STOP
            ENDIF
            fit =ABS(cc)**2/(e1*e2)
            WRITE(prtfil,'(a,i4,f9.2)')
     &           'Bartlett power (dB) for freq',ifrq,10*LOG10(1-fit) 
         ENDDO                  ! freq loop
         GOTO 900
      ENDIF
c     
 900  CONTINUE
      
      END

c************************************************************
      SUBROUTINE chisqr(xmean,xmeancost)
      USE global
      REAL xmean(*),xmeancost
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      REAL fit,fit2,e1,e2
      DOUBLE PRECISION fitlong
      INTEGER i,j,ix,ifrq,ibart,j1start,ndof,index,index2,index3
      INTEGER mdof,info
      DOUBLE COMPLEX cc,cx
      DOUBLE COMPLEX zx(mdep,mdep) ! matrix to find SVD of 
      DOUBLE COMPLEX ssvd(mdep+1),esvd(mdep),usvd(mdep,mdep)
      DOUBLE COMPLEX vsvd(mdep,mdep),work2(mdep)            
      WRITE(*,*)'xmean(1)',xmean(1)
      WRITE(prtfil,*)'entering chisqr...'
      CALL setmodelreal(xmean)
c     CALL setmodel(1)
      CALL forw2
      IF (iopt(3).GE.1) CALL norm ! normalize the response 

c     
      fitlong=0.
      fit=0.
      fit2=0.
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
               cc=(0.d0,0.d0) 
               e1=0.
               j1start=((ifrq-1))*ndep ! where the response starts
               DO j=j1start,ndep+j1start-1
                  index=ix+(j)*nx
                  e1=e1+REAL(resp(index))**2+AIMAG(resp(index))**2 
               ENDDO
c---  find the first eigenvector
               DO j=1,ndep
                  DO i=1,ndep
                     index2=i+(j+index3)*ndep 
                     zx(i,j)=cov(index2)             
                  ENDDO 
               ENDDO 
               
               CALL zsvdc(zx,mdep,ndep,ndep,ssvd,esvd,
     &              usvd,mdep,vsvd,mdep,work2,21,info)
               IF (info.NE.0) THEN
                  WRITE(*,*)'********* info from zsvdc',
     &                 ' ************',info
                  STOP
               ENDIF

               cx=(0.d0,0.d0)
               e2=(0.d0,0.d0)
               DO j=1,ndep
                  index=ix+(j+j1start-1)*nx
                  cx=cx+CONJG(usvd(j,1))*resp(index)     
                  e2=e2+CONJG(usvd(j,1))*usvd(j,1)    
c     WRITE(*,*)'cx',cx,cov(index2),resp(index),index2,index
               ENDDO
               cc=cc/e1         ! normalize with replica
               fit=fit+(1.0-cx/SQRT(e1))
               WRITE(*,*)'fit,e1,e2',fit,e1,e2
               WRITE(*,*)'cx',cx
c     WRITE(*,*)'Bartlett for ifreq',ifrq,10*LOG10(REAL(cc))
c     fit2=fit2+1-cc
 101        CONTINUE            !nx
 100     CONTINUE               !nfrq 
c     fit2=fit2/nbart
      ENDIF                     ! FINISHED WITH COVARIANCE
C     

      mdof=(2*nfrq*ndep)
      ndof=(2*nfrq*ndep-nparm)

      fit=fit/ndof

c      OPEN(unit=8,file='modelorder', status='unknown')
c      WRITE(*,*)
c      WRITE(*,*)' Chi squared test ',fit,fit*ndof,fit2
c      WRITE(*,'(3(a,i4))')
c     1     ' Nparm =',nparm,' Nobs =',nx*ncurv,' NDOF =',ndof
c      WRITE(8,'((a15,3i4,5f10.3))')
c     1     title,nparm,mdof,ndof,xmeancost,
c     1     xmeancost+1.0*nparm/mdof,5.96890970E-03/fit !3.1356984E-03/fit 
c     1     ,        4.1089943E-03/fit !4.1055516E-03/fit 
c      CLOSE(8)
      END


c**********************************************************************
      SUBROUTINE  svdmax(ibart,smax,noise)
C     *** find the largest singular value of the covariance matrix   
      USE global
      DOUBLE PRECISION  smax,noise
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER idep,jdep,ibart
      INTEGER info,index2
      COMPLEX zx(mdep,mdep)     ! matrix to find SVD of 
      COMPLEX ssvd(mdep+1),esvd(mdep),usvd(mdep,mdep),vsvd(mdep,mdep)
      COMPLEX work2(mdep)         

c      WRITE(*,*) ' entering into svdmax....'
c     
c---  find the first eigenvector
      DO jdep=1,ndep
         index2= (jdep+(ibart-1)*ndep-1)*ndep 
         DO idep=1,ndep
            zx(idep,jdep)=cov(idep+index2)             
         ENDDO 
      ENDDO 
      
      CALL csvdc(zx,mdep,ndep,ndep,ssvd,esvd,usvd,mdep,vsvd,mdep,
     &     work2,21,info)
      IF (iopt(6) .EQ. 1) THEN
         print *, "Matrix cov usvd follows"
         call PRINTMAT(usvd,ndep)
      ENDIF
      noise = 0
      DO jdep=2,ndep
C     write(*,*) 'jdep, ssvd',jdep, ssvd(jdep)
         noise=noise+ssvd(jdep)
      ENDDO
      noise=noise/(ndep-1)
      WRITE(*,*)' Largest SVD:', ssvd(1)
      WRITE(*,*)' Second SVD:', ssvd(2)
      WRITE(*,*)' Noise SVD:', noise
c     noise=ssvd(2)
      smax=ssvd(1) 
      IF (info.NE.0) THEN
         WRITE(*,*)'********* info from csvdc ************',info
         STOP
      ENDIF
 1000 CONTINUE
      END


C     .cfh.
c**********************************************************************
      subroutine PRINTMAT(x,ndep)
      INTEGER idep,jdep,ndep
      COMPLEX x(ndep,ndep)         
      do idep = 1, ndep
         write(*,*) idep ,' column vector'
         write(*,*) (x(jdep,idep), jdep=1,ndep)
      enddo
      end subroutine PRINTMAT
