c234567890123456789012345678901234567890123456789012345678901234567890
c        1         2         3         4         5         6         7
      subroutine forwardmodel(iopt,mopt)
      integer  mopt,i,iopt(mopt)
      DO i=1,mopt
        iopt(i)=0.
      ENDDO
        iopt(30)=7   ! slave program
      end
c******************************
      SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization program
c     PETER GERSTOFT, 1992
c
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'
      INTEGER iii,m,i,j,iq,jn,ierr,jj,irflag,ncov_siz,ipar,ioer,idum
      REAL delf,rdstep,zloc
      real*8  ztemp
      ierr=0
      nfrq=1
c---  open slave files:
      call opfilw(31,ioer)
      call opfilw(32,ioer)
c---  read the GA-parameters
      call readinputstart

c---- read the frequencies
         READ(1,*,err=92) nfrq
92        write(*,'(a,i3,a,i4)') ' nfrq=', nfrq
        if (nfrq.gt.mfreq) stop ' nfreq > mfreq'    ! test
        READ(1,*)(frq(I),i=1,nfrq)
        write(*,'(a,10f10.3)')' frequencies:',(frq(I),i=1,nfrq)
     
      read(1,*)nx
      read(1,*)ndep
      write(*,*)' nx, ndep=', nx, ndep
      read(1,*)nslavepar
      write(*,*)' nslavepar',nslavepar

      do ipar=1,nslavepar
         read(1,*)idum,slavepar(ipar)
         write(*,*)ipar,slavepar(ipar)
      enddo

      IF (ndep.gt.mdep) THEN
         write(*,*) ' ndep must grater than mdep',ndep,mdep
         ierr=1
      ENDIF

c****
      nfrq=1
      ncurv=nfrq*ndep
      if ((iopt(2).eq.3).or.(iopt(13).eq.1)) then
        ncurv=nfrq*nx
      endif
      if (iopt(13).gt.0) then  ! for covariance matrix
        nbart=nfrq*nx
         WRITE(prtfil,*)'number of bartlett estimators',nbart
         if (nbart.gt.mfreq) then
           WRITE(*,*)'nbart > mfreq'
           ierr=1
         ENDIF
         ncov_siz=nbart*ndep*ndep
         write(*,*)' Elements in cov-matrix:',ncov_siz
         if (ncov_siz.gt.mobs) then
           WRITE(*,*)'ncov_siz > mobs,  mobs=',mobs
           ierr=1
         ENDIF
      endif
      WRITE(prtfil,*)' ncurv=',ncurv,'nbart=',nbart
      if ((nx.eq.1).and.(ndep.eq.1)) then
         write(*,*) ' nx and ndep must not both be 1'
         ierr=1
      ENDIF
       if (nx*ndep.gt.mobs) then
         write(*,*) ' nx*ndep >mobs, nx,ndep,mobs',nx,ndep,mobs
         ierr=1
      ENDIF

      if (iopt(17).eq.1) then
        call eofinit
      endif

c*** create text strings for physical parameters
c                   123456789012345678901234567890
      phystxt(1) = 'Parameter $'
      phystxt(11)= 'Shape coefficient $'
c***  

      do 8 j=1,mphys
       phystxt2(j)='                                               '
    6 DO 7 I=40,1,-1
      IF(phystxt(j)(I:I).NE.'$') GO TO 7
      phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c      write(*,*)phystxt2(j)
      GO TO 8
    7 CONTINUE
    8 CONTINUE
          
c---- read the optimization parameter
      call readoptparm

      DO i=1,nparm
        IF (fmin(i).gt.fmax(i))THEN
          WRITE(*,*)' *** fmin > fmax for parm',i
          ierr=1
        ENDIF
        IF (par2phy(i).ne.1 .and. par2phy(i).ne.11)THEN
          WRITE(*,*)' *** par2phy not correct, parm',i
          write(*,*)' par2phy:',par2phy(i)
          ierr=1
        ENDIF
        IF (par2phy(i).eq.11 )THEN
          if (par2lay(i).gt.neof) then
          write(*,*)' *** Optimzation variable #:',i
          WRITE(*,*)' *** The shapecoffient number is not defined'
          WRITE(*,*)' par2lay(i)', par2lay(i)
          ierr=1
          endif
        ENDIF
      ENDDO

c*** errors ?
      IF (ierr.eq.1)stop 'stopped in input'
c***  close input unit
      CLOSE(1)                                    
      END

c**************************************************************
      subroutine forwinit
c SNAP subroutine to interface to snap- initialization of parameters
      character*4 titdum(20)
      integer ipar,i,id,index 
      integer infomode
      common /infomode/infomode
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comslave.h'
      write(*,*)' uses slave as forward model'
      WRITE(prtfil,*)' uses slave as forward model'

      flagpu=-0
      if (iopt(6).eq.1) flagpu=-10
c***************************************************
      entry forw2()
c***************************************************
       lwascomp=1
c      write(*,*) 'has entered snap'
      flagpu=1 +flagpu
c
c---- Use of EOF
c
      if (iopt(17).eq.1) then
        call eofval
        if (lwascomp.eq.-1) then
        return
        endif
      endif
c
      rewind(31)
      do ipar=1,nslavepar
        write(31,*)slavepar(ipar)
      enddo
      rewind(31)

      call system('slave')

      rewind(32)
c      if ((flagpu.le.2)) then
c        WRITE(prtfil,*)'freq,number of modes',ifreq,modqty
c      endif 

c         write(*,*)'snapinit:np,srd(15,1,1)',np,srd(15,1,1)     
      do id=1,ndep     ! irin
        index=(id -1)*nx
        do i=1,nx         ! ranges
           read(32,*)  resp(i+index)
c            write(*,*)index+i,fldpr(i+(id-1)*np)
        enddo  ! ranges
        if (flagpu.lt.0) then
          do i=1,nx         ! ranges
            WRITE(prtfil,*) resp(i+index),i+index
          enddo  ! ranges
        endif
      enddo   ! ir

      lwascomp=1
      return
      end       
