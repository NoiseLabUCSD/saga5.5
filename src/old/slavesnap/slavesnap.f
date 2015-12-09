      PROGRAM slave
      INCLUDE '../comopt.h'
      INCLUDE '../comforw.h'
c
      open(unit=1,STATUS= 'OLD',file='slave.input')
      call input
      close(1)
      open(unit=31,STATUS= 'OLD',file='slave.par')
      read(31,*)xranges(1)  !source range  
      read(31,*)sd          !source depth   
      call forwinit
      Write(*,*)'Slave finnished'
c
      end

c******************************
      SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization program
c     PETER GERSTOFT, 1992
c
      INCLUDE '../comopt.h'
      INCLUDE '../comforw.h'
      INTEGER iii,m,i,j,iq,jn,ierr,jj,irflag,ncov_siz
      REAL delf,rdstep,zloc
      real*8  ztemp
      INCLUDE '../comsnap.h'
      ierr=0
      nfrq=1
c---  read the GA-parameters
c      call readinputstart
      
      read(1,*)  ! this is the title line
      READ(1,'(a)') opt      ! This line has no effect in snap
      call snapoption(opt)
      iopt(1)=2

        write(*,*)'reading nfrq..'
        READ(1,*,err=91) nfrq, maxnomode
        goto 92
91      backspace(1)
        READ(1,*,err=91) nfrq
        maxnomode=500
92        write(*,'(a,i3,a,i4)') ' nfrq=', nfrq,' maxmode=' ,maxnomode
        if (nfrq.gt.mfreq) stop ' nfreq > mfreq'    ! test
        READ(1,*)(frq(I),i=1,nfrq)
        write(*,'(a,10f10.3)')' frequencies:',(frq(I),i=1,nfrq)

c---  now we are reading in a snap format
        read(1,*)ztemp,scatt(1),scatt(2),beta(0) 
        write(*,*)' hw=',ztemp
        hw=ztemp
        write(*,*)' Water speed profile:'
        do m=1,maxdep
          read(1,*)z0(m),c0(m)
          write(*,*)m,z0(m),c0(m)
          if (z0(m).eq.ztemp) goto 101
        enddo
        stop 'maxdep points read in for z0'

101     nwater=m
        nd0=m
        WRITE(prtfil,*)' Points in water:',nd0
        read(1,*)h1,r1,beta(1) 
        write(*,*)' h_sed =',h1
        IF(H1.GT.0.0)   then
         write(*,*)' Sediment speed profile:'
         do m=1,maxdep
           read(1,*)z1(m),c1(m)
           write(*,*)m,z1(m),c1(m)
           if (z1(m).eq.h1) goto 102
         enddo
         stop 'maxdep points read in for z1'
102      nd1=m
        else      ! no sediment
         ND1=0
         BETA(1)=0.0
         R1=0.0
        endif
         write(*,*)' Points in sediment:',nd1
         write(*,*)' Bottom parameters:'
         read(1,*)r2,beta(2),c2
         write(*,'(a,f7.2,a,f8.4,a,f10.2)')
     1            '  Dens:',r2,' P-att:',beta(2),' P-vel:',c2
         read(1,*)beta(3),c2s 
         write(*,'(a,f8.4,a,f9.2)')
     1            '  S-att:',beta(3),' S-vel:',c2s
        write(*,*)
        READ(1,*) sd         
c
c     receiver data
c
      if (tilt) then
        READ(1,*) rd,rdlow,ndep,dtilt
      else
        READ(1,*) rd,rdlow,ndep
      endif
      IF (ndep.gt.mdep) THEN
         write(*,*) ' ndep must grater than mdep',ndep,mdep
         ierr=1
      ENDIF
      irflag=ndep
      ndep=abs(ndep)
      IF (ndep.gt.1) THEN
        rdstep=(rdlow-rd)/float(ndep-1)
      ELSE
        rdstep=1.
      ENDIF
      IF (irflag.gt.0) THEN
c       equidistant receiver depths
        DO 920 jj=1,ndep
920       rdep(jj)=(jj-1)*rdstep+rd
      ELSE
c       read in individual receiver depths
        READ(1,*) (rdep(jj),jj=1,ndep)
      ENDIF
      WRITE(prtfil,*)' number of receiver depths read',ndep

c****
      WRITE(prtfil,*)' reading ranges in....'
201     READ(1,*)nx
        READ(1,*,err=203)(xranges(i),i=1,nx)
        WRITE(prtfil,*) 'number of ranges',nx
 203    continue
205   continue




c*** errors ?
      IF (ierr.eq.1)stop 'stopped in input'
c***  close input unit
      CLOSE(1)                                    

      END
      subroutine snapoption(options)
      INCLUDE '../comopt.h'
      INCLUDE '../comforw.h'
      INCLUDE '../comsnap.h'
      INCLUDE '../snap/common.f'
      character options(40)
      integer i
      tilt=.false.
      do 50 i=1,40
        if (options(i).eq.'i') then
          flag(1,6)=1
          write(*,*)'incoherent addition of modes'
        ELSE if (options(i).eq.'t') then
          tilt=.true.
          write(*,*)' tilt of array is included'
          write(prtfil,*)' tilt of array is included'
        ELSE IF (options(I).EQ.'!') THEN
          goto 60
        ELSE IF (options(I).NE.' ') THEN
          WRITE(prtfil,399) options(I)
 399      FORMAT(1H ,' >>>> unknown snap OPTION: ',A1,' <<<<')
        END IF
50    CONTINUE
60    CONTINUE
      end

c**************************************************************
      subroutine forwinit
c SNAP subroutine to interface to snap- initialization of parameters
      character*4 titdum(20)
      integer i,ifreqy,isour,isub,id,modqty,index 
      integer infomode
      common /infomode/infomode
      INCLUDE '../comopt.h'
      INCLUDE '../comforw.h'
      INCLUDE '../snap/common.f'
      INCLUDE '../comsnap.h'
      INCLUDE '../snap/a.f'
      write(*,*)' uses snap as forward model'
      WRITE(prtfil,*)' uses snap as forward model'

       if (maxnomode.gt.moden) then
         write(*,*)' Warning maxnumber of modes reduced to',moden
         write(*,*)' Snap does not accommodate more modes'
         maxnomode=moden
       endif

c---  transfer water properties
      
      r0=1.0				  ! density in water
c---  transfer receiver parametes
      fldrd(1)= rd		      ! minimum receiver depth
      fldrd(2)= rdlow		      ! maximum receiver depth
      if (ndep.eq.1) then
        fldrd(3)= (rdlow-rd)          ! step in receiver depth
      else
        fldrd(3)= (rdlow-rd)/(ndep-1)   ! step in receiver depth
      endif
      srd(15,1,1)=sd               ! srd(isub, up to 20 sd, not app. to field)
c---  range parameters
      np= nx                           ! number of range steps
      do i=1,np
        rng(i)= xranges(i)               ! range points
      enddo
      hstart=h0				  ! also water depth !!!  
      infomode=0
      flagpu=-2
      CORREC=-1.0
c***************************************************
      entry forw2()
c***************************************************
       lwascomp=1
c      write(*,*) 'has entered snap'
      flagpu=1 +flagpu
c
      h0=hw
      hstart=h0				  ! also water depth !!!  
      z0(nd0)=h0	      ! last point is the water depth
c
        z0(nwater)=hw
        if (nd1.gt.0) z1(nd1)=h1
C***********
C     NORMALIZATION OF DEPTHS
C
      CMINsnap=1.0D38
      cc0=0
      DO 1450   I=1,ND0
        CMINsnap=MIN(CMINsnap,C0(I))
        Zz0(I)=Z0(I)/H0
 1450 CONTINUE
      DO    I=1,ND0-1
        cc0= cc0 + (c0(i+1) + c0(i))*.5*(z0(i+1) - z0(i))
      enddo
      cc0=cc0/h0

      cc1=0
      if (h1.gt.0) then
        DO 1500   I=1,ND1
          CMINsnap=MIN(CMINsnap,C1(I))
          Zz1(I)=(H0+Z1(I))/H0
 1500   CONTINUE
        DO I=1,ND1-1
          cc1= cc1 + (c1(i+1) + c1(i))*.5*(z1(i+1) - z1(i))
        enddo
        cc1=cc1/h1
      endif
      if (flagpu.lt.2) then
       WRITE(prtfil,*)
       WRITE(prtfil,*)'Calling SNAP' 
       WRITE(prtfil,*)'flagpu',flagpu
c      WRITE(prtfil,*)'max no of modes:',maxnomode
       WRITE(prtfil,*)' Frequency:',(frq(ifreq),ifreq=1,nfrq) 
       WRITE(prtfil,*)' R1,R2',R1,r2
       WRITE(prtfil,*)' SCATT(1),SCATT(2)',SCATT(1),SCATT(2)
       WRITE(prtfil,*)' H0,H1',H0,H1
       WRITE(prtfil,*)' ND0,ND1',ND0,ND1
       WRITE(prtfil,*)' BETA(1),BETA(2),BETA(3)',
     &                 BETA(1),BETA(2),BETA(3)
c       WRITE(prtfil,*)' cC0,cC1' ,cC0,cC1
       WRITE(prtfil,*)' c2,c2s' ,c2,c2s
       WRITE(prtfil,'(a)')'zz0,cc0'
       WRITE(prtfil,'(a,/,100(i3,3F10.3/))')
     &     ' z0,zz0,c0',(i,z0(i),zz0(i),c0(I),i=1,nd0)
       IF (ND1.GT.0) then
         WRITE(prtfil,'(a,/,100(i3,3F10.3/))')
     &    ' z1,zz1,c1',(i,z1(i),zz1(i),c1(I),i=1,nd1)
       endif
       if (tilt) write(luprt,*)' Tilt (m)=',dtilt
         WRITE(prtfil,*)' np,srd(15,1,1)',np,srd(15,1,1)     
         WRITE(prtfil,*)' fldrd(1-3)',fldrd(1),fldrd(2),fldrd(3)
         WRITE(prtfil,*)' ranges',(rng(i),i=1,min(np,20))
       endif
     
      do ifreq=1,nfrq
       if (ifreq.ge.2)flagpu=flagpu+3
       freqy=dble(frq(ifreq))
       OMEGA=TWOPI*FReqy
c       WRITE(prtfil,*)'calling snap;freqy ',freqy
c       write(*,*)'calling snap;freqy ',freqy
      PHVMAX=C2
      PHVMIN=0.
      MINMOD=1
      MAXMOD=maxnomode
      CALL PORTER(FReQy,MINMOD,MAXMOD,NMES,JF2,DELFRQ,
     &  MODPLT,EK,EGV,XTS,CC0,CC1,ALFA,MODAVR,ADA,SPEED,
     &  EIGF,ISO,NBEG,MY,C0,Zz0,C1,Zz1,
     &  A3,B3,C3,EE,ZZ,SSOLD,EXCH,slow)

      MODQTY=MAXMOD-MINMOD+1
      IF(MODQTY.EQ.1)  then
        write(*,*)' Warning: snapinit found only one mode'
      elseIF((MODQTY.EQ.maxnomode).and.(infomode.le.100))  then
        write(*,*)' Warning: snapinit max number of modes reached'
        infomode=infomode+1
      elseIF(MODQTY.EQ.0)  then
        write(*,*)' Warning: snapinit found no  modes...returning'
        lwascomp=-1
        return
      endif
      if ((flagpu.le.20)) then
        WRITE(prtfil,*)'freq, number of modes',ifreq,modqty
      endif 

       isub=15                      ! we are selecting field option     
c---  transfer source parametes
c       maxmod=1 ! only one mode
       ifreqy=1

       
       call FIELD(NP,frq(ifreq),TITdum,ISUB,IFREQy,MSP,XTS,
     & TEMPOR,RNG,FLDPR,US,UR,SECD,ALFA,EK,TLUS,SRD,C0,Z0,
     & ND0,ISD,MODEN,NOPT,KSRD)

1000   continue  ! normal exit from field  
       open(unit=32,STATUS= 'UNKNOWN',file='slave.obs')
       
c         write(*,*)'snapinit:np,srd(15,1,1)',np,srd(15,1,1)     
	do id=1,ndep     ! irin
           index=(id +((ifreq-1))*ndep-1)*np
          do i=1,np         ! ranges
c           resp(i+index)=fldpr(i+(id-1)*np)

            write(32,*)fldpr(i+(id-1)*np)
          enddo  ! ranges
          if (flagpu.lt.0) then
          do i=1,np         ! ranges
            WRITE(prtfil,*)fldpr(i+(id-1)*np) ,i+index
          enddo  ! ranges

          endif
          if (abs(fldpr(np)).le.0) then 
            Write(*,*)'the pressure is zero; a dump follows'
            WRITE(prtfil,*)'the pressure is zero;  a dump follows'
            do i=1,np         ! ranges
              WRITE(prtfil,*)rng(i),fldpr(i+(id-1)*np)
            enddo 
            WRITE(prtfil,*)'number of modes',modqty
            WRITE(prtfil,*)' Frequency:',frq(1)  
            WRITE(prtfil,*)' R1,R2',R1,r2
            WRITE(prtfil,*)' SCATT(1),SCATT(2)',SCATT(1),SCATT(2)
            WRITE(prtfil,*)' H0,H1',H0,H1
            WRITE(prtfil,*)' ND0,ND1',ND0,ND1
            WRITE(prtfil,*)' BETA(1),BETA(2),BETA(3)',
     &                      BETA(1),BETA(2),BETA(3)
            WRITE(prtfil,*)' c2,c2s' ,c2,c2s
              WRITE(prtfil,'(a,/,100(i3,3F10.3/))')
     &      ' z0,zz0,c0',(i,z0(i),zz0(i),c0(I),i=1,nd0)
            IF (ND1.GT.0) then
             WRITE(prtfil,'(a,/,100(i3,3F10.3/))')
     &        ' z1,zz1,c1',(i,z1(i),zz1(i),c1(I),i=1,nd1)
            endif
          endif
        enddo   ! ir
        if (ifreq.ge.2)flagpu=flagpu-3
       enddo  ! loop over frequencies

       lwascomp=1
      return
      end       


