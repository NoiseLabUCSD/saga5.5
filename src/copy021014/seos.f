      Program seos
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER nloop/1000/
      REAL tempr(100),tempi(100)
      INTEGER ifreq,idep,i
      CHARACTER*2 ch
      CHARACTER*80 filein,fileout
      ndep=100
      WRITE(*,*)' This program translate a filename in seos format'
      WRITE(*,*)' into SAGA hydrophone pressure vector format'
      WRITE(*,*)' '
c     >>> input file names
 1    WRITE(6,*) 'Input file?'
      READ(5,'(A)') filein
      OPEN(9,FILE=FILEin,STATUS='OLD',FORM='FORMATTED',ERR=10)
      GOTO 20 
 10   WRITE(*,*)' ***** input file does not exist, try again'
      GOTO 1
 20   WRITE(6,*) 'Output file?'
      READ(5,'(A)') fileout
      OPEN(30,FILE=FILEout,STATUS='unknown',FORM='FORMATTED')
      title='seos file: '//filein
      DO idep=1,ndep
         rdep(idep)=idep
      ENDDO
c     OPEN(file='respr.m',unit=98,STATUS= 'UNKNOWN')
      Do ifreq=1,nloop
         READ(9,'(a2,f6.2)',err=100)ch,frq(ifreq)
         WRITE(*,*)ch,frq(ifreq)
         READ(9,'(a2,100f8.3)')ch,(tempr(idep),idep=1,ndep)
         READ(9,'(a2,100f8.3)')ch,(tempi(idep),idep=1,ndep)
         DO idep=1,ndep
            i = (ifreq-1)*Ndep + idep
            resp(1+(i-1)*1)=CMPLX(tempr(idep),tempi(idep))/1000
         ENDDO
      ENDDO
 100  CONTINUE

      nfrq=ifreq-1
      WRITE(*,*)' Number of frequencies:',nfrq

      WRITE(30,'(a1,a32)') '!',' File format: Hydrophone vectors' 
      DO ifreq=1,nfrq
         CALL WRITE_HP(ifreq)
      ENDDO
      END

c********************************************************************
      SUBROUTINE WRITE_HP(ifreq)

C     Writes the generated hydrophone fourier DATA to the OBS file
C     in the "hydrophone data" FORMAT. This consists of COMPLEX vectors
C     of the hydrophone output at each frequency.  Note this routine
C     is called once for each required frequency.    

      INTEGER OBSFIL
      INTEGER i, idep, ifreq

      PARAMETER ( OBSFIL = 30 )

      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      WRITE( OBSFIL, * ) TITLE(1:60) ! Title
      WRITE( OBSFIL, * ) FRQ(ifreq), Nfrq ! Current frequency,(optional: NFRQS)

      WRITE( OBSFIL, * ) Ndep
      WRITE( OBSFIL, 100 )( rdep(i), i=1,Ndep ) ! The actual depths.
 100  FORMAT( 1F10.2 )

C     *** WRITE pressure for current frequency to the OBS file. ***
C     Note resp() CONTAINS the Nfrq fourier components at each hydrophone
C     depth.  Frequency varies slowest.  The first index represents range.
      DO idep = 1, Ndep
         i = (ifreq-1)*Ndep + idep
         WRITE( OBSFIL, * ) idep, resp(1+(i-1)*1)
     .        , ABS(resp(1+(i-1)*1))
      END DO
      
      END
