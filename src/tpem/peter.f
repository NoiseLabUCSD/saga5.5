C *** WRITE TRANSFER FUNCTION FILE
c       DO 300 IS=1,ISROW
c       DO 300 M=1,MSUFT
c       DO 300 JRH=1,NPLOTS
c       DO 300 JRV=1,IR
c        do i=1,nout
c        cdummy(i)=CFFX(I,JRV,JRH,M,IS)
c        end do
c        if (bintrf) then
c         WRITE(LUTTRF) (cdummy(i),I=1,NOUT)
c        else
c         WRITE(LUTTRF,*) (real(cdummy(i)),
c     &                   imag(cdummy(i)),I=1,NOUT)
c        end if
c 300   CONTINUE
c      RETURN



      SUBROUTINE TRFHEAD(trfext,TITLE,RD,RDLOW,R0,RSPACE,SD)
c      INCLUDE 'compar.f'
c      INCLUDE 'comnp.f'
c      INCLUDE 'comnrd.f'
      include 'tpem.inc'
      include 'tpem.h'      
      CHARACTER*(*) trfext
      CHARACTER*8 FILEID
      character*16 bufch
      CHARACTER*80 TITLE,filenm
      CHARACTER SIGNN
      LOGICAL bintrf
      CHARACTER*6 PROGNM
      INTEGER IPARM(12)
      integer idummy,nout,inttyp,ir,icdr
c >>> Dummy real*4 to force single prec trf file.
      real*4 dummy,adummy(10)
      bufch=trfext
c      inquire(1,name=filenm)
c      ii=indexs(filenm,'.')
      nout=1
      msuft=1
      PROGNM='TPEM  '
      
      ir=vnp.nzout
      nplots=vnp.nrout
      icdr=1
      omegim=0
      inttyp=0
      bintrf=.false.
      if (bintrf) then
c        OPEN(LUTTRF,FILE=Filenm(1:ii)//bufch,
c     &          STATUS='UNKNOWN',FORM='UNFORMATTED')
cunix        OPEN(LUTTRF,FILE=FILEnm(1:ii)//trfext,
cunix     &       STATUS='UNKNOWN',FORM='UNFORMATTED')
c        FILEID='PULSETRF'
c        WRITE(LUTTRF) FILEID
c        WRITE(LUTTRF) PROGNM
c        WRITE(LUTTRF) NOUT
c        ICNT=1
c        DO 10 I=1,NPAR
c         IF (IOUT(I).NE.0) THEN
c           IPARM(ICNT)=I
c           ICNT=ICNT+1
c         END IF
c 10     CONTINUE
c        WRITE(LUTTRF) (IPARM(J),J=1,NOUT)
c        WRITE(LUTTRF) TITLE
c        SIGNN='+'
c        WRITE(LUTTRF) SIGNN
c        dummy=freqs
c        WRITE(LUTTRF) dummy
c        dummy=sd
c        WRITE(LUTTRF) dummy
c       adummy(1)=rd
c        adummy(2)=rdlow
c        WRITE(LUTTRF) adummy(1),adummy(2),IR
cc        IF (IR.LT.0) THEN
c         do L = 1,abs(ir)
c          dummy=rdc(L)
c          WRITE(LUTTRF) dummy
c         end do
c         write(luttrf) (rdc(l),l=1,abs(ir))
c        END IF
c        adummy(1)=r0
c        adummy(2)=rspace
c        WRITE(LUTTRF) adummy(1),adummy(2),NPLOTS
c        dummy=dt
c        WRITE(LUTTRF) NX,LX,MX,dummy
c        WRITE(LUTTRF) ICDR
c        dummy=omegim
c        WRITE(LUTTRF) dummy
C ***  EXTRA FIELDS ADDED 891211 HS
c        WRITE(LUTTRF) MSUFT
c        write(6,*) 'trfhead: msuft=',msuft
c        WRITE(LUTTRF) ISROW
c        write(LUTTRF) inttyp
c        idummy=0
c        DO 300 I=1,2
c         WRITE(LUTTRF) IDUMMY
c 300    CONTINUE
c        dummy=0
c        DO 400 I=1,5
c         WRITE(LUTTRF) DUMMY
c 400    CONTINUE
        else
c        OPEN(LUTTRF,FILE=Filenm//bufch,
c pg        OPEN(LUTTRF,FILE='tpem.trf',
c pg    &          STATUS='UNKNOWN',FORM='FORMATTED')
cunix        OPEN(LUTTRF,FILE=FILEnm(1:ii)//trfext,
cunix     &       STATUS='UNKNOWN',FORM='UNFORMATTED')
        FILEID='PULSETRF'
        WRITE(LUTTRF,'(1x,a)') FILEID
        WRITE(LUTTRF,'(1x,a)') PROGNM
        WRITE(LUTTRF,*) NOUT
c        ICNT=1
c        DO 11 I=1,NPAR
c         IF (IOUT(I).NE.0) THEN
c          IPARM(ICNT)=I
c           ICNT=ICNT+1
c         END IF
 11     CONTINUE
        WRITE(LUTTRF,*) 1         !(IPARM(J),J=1,NOUT)
c        title='                                                    '
        WRITE(LUTTRF,'(1x,a80)') TITLE
        SIGNN='+'
        WRITE(LUTTRF,'(1x,a)') SIGNN
        WRITE(LUTTRF,*) FREQS
        WRITE(LUTTRF,*) SD
        WRITE(LUTTRF,*) RD,RDLOW,IR
c        IF (IR.LT.0) THEN
c         WRITE(LUTTRF,*) (RDC(L),L=1,ABS(IR))
c        END IF
        WRITE(LUTTRF,*) R0/1000,RSPACE/1000,NPLOTS   ! transformed to km

        WRITE(LUTTRF,*) NX,LX,MX,DT
        WRITE(LUTTRF,*) ICDR
        WRITE(LUTTRF,*) OMEGIM
C ***  EXTRA FIELDS ADDED 891211 HS
        WRITE(LUTTRF,*) MSUFT
c        write(6,*) 'trfhead: msuft=',msuft
        WRITE(LUTTRF,*) ISROW
        write(LUTTRF,*) inttyp
        idummy=0
        DO 301 I=1,2
           WRITE(LUTTRF,*) IDUMMY
 301    CONTINUE
        dummy=0e0
        DO 401 I=1,5
           WRITE(LUTTRF,*) DUMMY
 401    CONTINUE
        end if
        RETURN
        
        END
c************************************TIMIMG**************************

      SUBROUTINE CLTIME
      real*4 etime
      DIMENSION SECS(2)
      COMMON /TIMECT/ OLDSECS
      OLDSECS = ETIME(SECS)
C      OLDSECS=SECS(1)
C      CALLSTAT=LIB$INIT_TIMER()
      RETURN
      END
      SUBROUTINE RDTIME(T1)
      real*4 etime
c
c Returns CPU time in secs.
c
C *** UNIX SYSTEMS
      DIMENSION SECS(2)
      COMMON /TIMECT/ OLDSECS
      T1=ETIME(SECS)-OLDSECS
C *** VMS SYSTEMS
CVMS  INTEGER*4 IT4
CVMS  REAL*8 IT8
CVMS  EQUIVALENCE (IT4,IT8)
CVMS  CALLSTAT=LIB$STAT_TIMER(2,IT4)
CVMS  T1=1E-2*IT4
      RETURN
      END
c
c******************* END TIMIMG***************************
c
      Integer Function Indexs(Str,Substr)
      Character*(*) Str,Substr
      L1=Len(Str)
      L2=Len(Substr)
      Do 10 I=L1-L2+1,1,-1
       If (Substr.Eq.Str(I:I+L2-1)) Then
        Indexs=I
        Return
       end if
 10   continue
      indexs=0
      return
      end



      SUBROUTINE OPFILW(IUN,IOER)
      CHARACTER*40 FILENM
      character*6 ENVVAR
      character*3 exten
      IOER=0
C
C     OPENS SEQUENCIAL ASCII FILE FOR WRITE
C
C***  VMS
C     OPEN(UNIT=IUN,STATUS='NEW',FORM='FORMATTED',ERR=300)
C***  ultrix 
      WRITE(ENVVAR,'(A3,I3.3)') 'FOR',IUN
      CALL GETENV(ENVVAR,FILENM)
      IF (FILENM.EQ.' ') THEN
       inquire(1,name=filenm)
       ii=indexs(filenm,'.')
       write(exten,'(i3.3)') iun
       OPEN(UNIT=IUN,file=filenm(1:ii)//exten,
     &      STATUS='UNKNOWN',FORM='FORMATTED',ERR=300)
      ELSE
       OPEN(UNIT=IUN,FILE=FILENM,STATUS='UNKNOWN',
     &      FORM='FORMATTED',ERR=300)
      END IF
C***  MS-DOS         
c      IF (IUN.LT.10) THEN
c        WRITE(FILENM,100) IUN
c      ELSE
c        WRITE(FILENM,200) IUN
c      END IF
c 100  FORMAT('FOR00',I1,'.DAT')
c 200  FORMAT('FOR0',I2,'.DAT')
c      OPEN(UNIT=IUN,FILE=FILENM,STATUS='UNKNOWN',FORM='FORMATTED',
c     1     ERR=300)
      RETURN
 300  IOER=1
      RETURN
C
      ENTRY OPFILR(IUN,IOER)
      IOER=0
C
C     OPENS SEQUENCIAL ASCII FILE FOR READ
C
C***  VMS
C     OPEN(UNIT=IUN,STATUS='OLD',ERR=400)
C***  ultrix 
      WRITE(ENVVAR,'(A3,I3.3)') 'FOR',IUN
      CALL GETENV(ENVVAR,FILENM)
      IF (FILENM.EQ.' ') THEN
       OPEN(UNIT=IUN,STATUS='OLD',FORM='FORMATTED',ERR=400)
      ELSE
       OPEN(UNIT=IUN,FILE=FILENM,STATUS='OLD',
     &      FORM='FORMATTED',ERR=400)
      END IF
C***  MS-DOS         
c      IF (IUN.LT.10) THEN
c        WRITE(FILENM,100) IUN
c      ELSE
c        WRITE(FILENM,200) IUN
c      END IF
c      OPEN(UNIT=IUN,FILE=FILENM,STATUS='OLD',ERR=400)
      RETURN
 400  IOER=2
      RETURN
      END
      SUBROUTINE OPFILB(IUN,IOER)
      CHARACTER*40 FILENM
      character*6 ENVVAR
      character*3 exten
      IOER=0
C
C     OPENS SEQUENCIAL ASCII FILE FOR WRITE
C
C***  VMS
C     OPEN(UNIT=IUN,STATUS='NEW',FORM='FORMATTED',ERR=300)
C***  ultrix 
      WRITE(ENVVAR,'(A3,I3.3)') 'FOR',IUN
      CALL GETENV(ENVVAR,FILENM)
      IF (FILENM.EQ.' ') THEN
       inquire(1,name=filenm)
       ii=indexs(filenm,'.')
       write(exten,'(i3.3)') iun
       OPEN(UNIT=IUN,file=filenm(1:ii)//exten,
     &      STATUS='UNKNOWN',FORM='UNFORMATTED',ERR=300)
      ELSE
       OPEN(UNIT=IUN,FILE=FILENM,STATUS='UNKNOWN',
     &      FORM='UNFORMATTED',ERR=300)
      END IF
C***  MS-DOS         
c      IF (IUN.LT.10) THEN
c        WRITE(FILENM,100) IUN
c      ELSE
c        WRITE(FILENM,200) IUN
c      END IF
c 100  FORMAT('FOR00',I1,'.DAT')
c 200  FORMAT('FOR0',I2,'.DAT')
c      OPEN(UNIT=IUN,FILE=FILENM,STATUS='UNKNOWN',FORM='UNFORMATTED',
c     1     ERR=300)
      RETURN
 300  IOER=1
      RETURN
      end




