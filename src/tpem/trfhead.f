
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
      
      ir=vnp%nzout
      nplots=vnp%nrout
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


