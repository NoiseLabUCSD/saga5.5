      SUBROUTINE writetrf()
      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'

      INTEGER ifreq,i,id,index,mx1
c     pg      INTEGER*4 jj0,jsrc,jrec,j,mx1
      REAL*4 rspace,freqs,dt1
      LOGICAL bintrf
      INTEGER luttrf
      data luttrf/16/
      CHARACTER*30 outroot      !not used
      bintrf=.FALSE.
c     bintrf=.TRUE.
c     
c     rspace = range increment for receiver position
c     pln      rkm(1)=rkm(1)*1.0e-3
      IF (np.EQ.1) THEN
         rspace=1
      ELSE
         rspace=rng(2)-rng(1)
      END IF
      mx1=nf1+nfrq-1
      dt1=1.E0/fsbb
      freqs=(fmaxdum-fmindum)/2.+fmindum
c     sd=zsr(mzsrc(1))
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     rng(1),rspace,nfftbb,nf1,mx1,dt1,freqs,sd,
     &     bintrf,ndep,1,np)
c     
      IF(iifft .NE. 0) THEN
c     
c     pg         IF (bintrf) THEN
c     pg            DO j=1,nfbb
c     pg               DO jsrc=1,nsrc
c     pg                  jj0=(jsrc-1)*nfbb*nrecusr
c     pg                  DO jrec=1,nrecusr
c     c fbv                     WRITE(95, *) -REAL(tf(j+jj0)), faxbb(j)
c     Output for Normal Stress
c     WRITE(luttrf)-REAL(tf(j+jj0)), AIMAG(tf(j+jj0))
c     WRITE(luttrf)-dreal(tf(j+jj0))
c     the previous line for the output on linux workstations
c     Output for Pressure
c     pg                     WRITE(luttrf)REAL(tf(j+jj0)), -AIMAG(tf(j+jj0))
c     pg                     jj0=jj0 + nfbb
c     pg                  END DO
c     pg               END DO
c     pg            END DO
c     pg         ELSE
c     pg            DO j=1,nfbb
c     pg               DO jsrc=1,nsrc
c     pg                  jj0=(jsrc-1)*nfbb*nrecusr
c     pg                  DO jrec=1,nrecusr
c     c fbv                     WRITE(95, *) REAL(tf(j+jj0)),
c     c fbv     &                            faxbb(j),jrec
c     Output for Normal Stress
c     WRITE(luttrf,*)-REAL(tf(j+jj0)),AIMAG(tf(j+jj0))
c     WRITE(luttrf,*)-dreal(tf(j+jj0))
c     the previous line for the output on linux workstations
c     Output for Pressure
c     pg                     WRITE(luttrf,*)REAL(tf(j+jj0)),-AIMAG(tf(j+jj0))
c     pg                     jj0=jj0 + nfbb
c     pg                  END DO
c     pg               END DO
c     pg            END DO
         DO ifreq=1,nfrq
            DO i=1,np           ! ranges
               DO id=1,ndep     ! irin
                  index=(id +((ifreq-1))*ndep-1)*np
                  WRITE(luttrf,*)REAL(resp(i+index)),
     1                 +imag(resp(i+index))
c     pg 1 may 2000 sign for prosim      -imag(resp(i+index))
               ENDDO            ! depth
            ENDDO               ! ranges
         ENDDO                  ! freq
c     pg         END IF
c     
c         CLOSE(luttrf)
      ENDIF
c     
      RETURN
      END
c     
      SUBROUTINE trfhead(filenm,title,rd,rdlow,r0,rspace,
     &     nx,lx,mx,dt1,freqs,sd,bintrf,ir,is,nplots)
c     
      IMPLICIT NONE
c     
c     
      LOGICAL bintrf
      INTEGER luttrf
      DATA luttrf/16/
      INTEGER*4 i,idummy,inttyp,isrow,icdr,msuft,nx,lx,mx,
     &     nplots,ir,is,icnt
      REAL*4 dummy,adummy(10),omegim,r0,rd,rspace,rdlow,dt1,
     &     sd,freqs
      CHARACTER*6 prognm
      CHARACTER*8 fileid
      CHARACTER*64 filenm
      CHARACTER*80 title
      CHARACTER signn
c     >>> Dummy REAL*4 to force single prec trf file.
c     
c     
c      write(*,*)' Entering trfhead...'

      omegim=0.
c     title(65:80)=' '
      IF (bintrf) THEN
         IF(is.GT.1) STOP'in bb_modes'
c     OPEN(luttrf, file='prosim_out.dat',
c     &        status='unknown',form='unformatted')
         fileid='PULSETRF'
         WRITE(luttrf) fileid
         prognm='PROSIM'
         WRITE(luttrf) prognm
c     nout=1
c     WRITE(luttrf) nout
         WRITE(luttrf) INT(1.)
c     icnt=1
c     WRITE(luttrf) icnt
         WRITE(luttrf) INT(1.)
         WRITE(luttrf) title(1:80)
         signn='+'
         WRITE(luttrf) signn
c     CENTER FREQUENCY
         dummy=freqs
         WRITE(luttrf) dummy
c     SOURCE DEPTH
         dummy=sd
         WRITE(luttrf) dummy
c     UPPER MOST RECEIVER DEPTH
         adummy(1)=rd
c     LOWER MOST RECEIVER DEPTH
         adummy(2)=rdlow
c     IR=NO OF RECEIVERS BETWEEN R0 AND RDLOW
         WRITE(luttrf) adummy(1),adummy(2),ir
C     
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
c     DO L = 1,ABS(ir)
c     dummy=rdc(L)
c     WRITE(LUTTRF) dummy
c     END DO
c     RDC(L) CONTAINS ALL THE RECEIVER DEPTHS
C     WRITE(luttrf) (rdc(l),l=1,ABS(ir))
C     END IF
C     
C     R0= THE FIRST RANGE TO PLOT IN RANGE STACKED PLOT
         adummy(1)=r0
C     RSPACE= THE RANGE INCREMENT TO PLOT IN RANGE STACKED PLOT
         adummy(2)=rspace
c     WRITE R0, RSPACE AND THE NO OF PLOTS ASSOCIATED WITH IR
         WRITE(LUTTRF) adummy(1),adummy(2),nplots
C     DT=TIME SAMPLING
         dummy=dt1
C     NX=NO OF TIME SAMPLES (DENOTED NT IN OASES MANUAL) MAYBE?
C     LX=INDEX OF FIRST FREQUENCY COMPONENT (INT(FR1*DT))
C     MX=INDEX OF LAST FREQUENCY COMPONENT (INT(FR2*DT))
         WRITE(luttrf) nx,lx,mx,dt1
         icdr=0
         WRITE(luttrf) icdr
         dummy=omegim
         WRITE(luttrf) dummy
C     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         WRITE(luttrf) msuft
c     WRITE(6,*) 'trfhead: msuft=',msuft
         isrow=1
         WRITE(luttrf) isrow
         inttyp=1
         WRITE(luttrf) inttyp
         idummy=0
         DO 300 i=1,2
            WRITE(luttrf) idummy
 300     CONTINUE
         dummy=0
         DO 400 i=1,5
            WRITE(luttrf) dummy
 400     CONTINUE
      ELSE
c      WRITE(6,*)'No of char in filenm: ',lout
c      WRITE(6,'(a)')'Filename: ',filenm(1:lout)
c     PAUSE
c     OPEN(luttrf, file='prosim_out.dat',
c     &        status='unknown',form='formatted')
c     OPEN(luttrf,file='pek.asc',
c     &        status='unknown',form='formatted')
         fileid='PULSETRF'
         prognm='SNAP '
         WRITE(luttrf,'(1x,a)') fileid
         WRITE(luttrf,'(1x,a)') prognm
c     WRITE(LUTTRF,*) NOUT
         WRITE(luttrf,*) INT(1.)
         icnt=1
         WRITE(luttrf,*) INT(1.)
         WRITE(luttrf,'(1x,a)') title(1:64)
         signn='+'
         WRITE(luttrf,'(1x,a)') signn
         WRITE(luttrf,*) freqs
         WRITE(luttrf,*) sd
c     IR=1
         WRITE(luttrf,*) rd,rdlow,ir
C     
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
C     WRITE(LUTTRF,*) (RDC(L),L=1,ABS(IR))
C     END IF
C     NPLOTS=1
C     
         WRITE(luttrf,*) r0,rspace,nplots
         WRITE(luttrf,'(3I6,E20.10)') nx,lx,mx,dt1
         icdr=0
         WRITE(luttrf,*) icdr
         WRITE(luttrf,*) omegim
c     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         WRITE(luttrf,*) msuft
c      WRITE(6,*) 'trfhead: msuft=',msuft
         isrow=1
         WRITE(luttrf,*) isrow
         inttyp=1
         WRITE(luttrf,*) inttyp
         idummy=0
         DO 301 i=1,2
            WRITE(luttrf,*) idummy
 301     CONTINUE
         dummy=0e0
         DO 401 i=1,5
            WRITE(luttrf,*) dummy
 401     CONTINUE
      END IF
      RETURN
      END


      SUBROUTINE trfhead2(filenm,title,rd,rdlow,r0,rspace,
     &     nx,lx,mx,dt1,freqs,sd,bintrf,ir,is,nplots)
c     
      IMPLICIT NONE     
      LOGICAL bintrf
      INTEGER luttrf
      DATA luttrf/16/
      INTEGER*4 i,idummy,inttyp,isrow,icdr,msuft,nx,lx,mx,
     &     nplots,ir,is,icnt
      REAL*4 dummy,adummy(10),omegim,r0,rd,rspace,rdlow,dt1,
     &     sd,freqs
      CHARACTER*6 prognm
      CHARACTER*8 fileid
      CHARACTER*64 filenm
      CHARACTER*80 title
      CHARACTER signn
c     >>> Dummy REAL*4 to force single prec trf file.
c     
      omegim=0.
c     title(65:80)=' '
      IF (bintrf) THEN
         IF(is.GT.1) STOP'in bb_modes'
c     OPEN(luttrf, file='prosim_out.dat',
c     &        status='unknown',form='unformatted')
         fileid='PULSETRF'
         WRITE(luttrf) fileid
         prognm='PROSIM'
         WRITE(luttrf) prognm
c     nout=1
c     WRITE(luttrf) nout
         WRITE(luttrf) INT(1.)
c     icnt=1
c     WRITE(luttrf) icnt
         WRITE(luttrf) INT(1.)
         WRITE(luttrf) title(1:80)
         signn='+'
         WRITE(luttrf) signn
c     CENTER FREQUENCY
         dummy=freqs
         WRITE(luttrf) dummy
c     SOURCE DEPTH
         dummy=sd
         WRITE(luttrf) dummy
c     UPPER MOST RECEIVER DEPTH
         adummy(1)=rd
c     LOWER MOST RECEIVER DEPTH
         adummy(2)=rdlow
c     IR=NO OF RECEIVERS BETWEEN R0 AND RDLOW
         WRITE(luttrf) adummy(1),adummy(2),ir
C     
C     R0= THE FIRST RANGE TO PLOT IN RANGE STACKED PLOT
         adummy(1)=r0
C     RSPACE= THE RANGE INCREMENT TO PLOT IN RANGE STACKED PLOT
         adummy(2)=rspace
c     WRITE R0, RSPACE AND THE NO OF PLOTS ASSOCIATED WITH IR
         WRITE(LUTTRF) adummy(1),adummy(2),nplots
C     DT=TIME SAMPLING
         dummy=dt1
C     NX=NO OF TIME SAMPLES (DENOTED NT IN OASES MANUAL) MAYBE?
C     LX=INDEX OF FIRST FREQUENCY COMPONENT (INT(FR1*DT))
C     MX=INDEX OF LAST FREQUENCY COMPONENT (INT(FR2*DT))
         WRITE(luttrf) nx,lx,mx,dt1
         icdr=0
         WRITE(luttrf) icdr
         dummy=omegim
         WRITE(luttrf) dummy
C     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         WRITE(luttrf) msuft
c     WRITE(6,*) 'trfhead: msuft=',msuft
         isrow=1
         WRITE(luttrf) isrow
         inttyp=1
         WRITE(luttrf) inttyp
         idummy=0
         DO 300 i=1,2
            WRITE(luttrf) idummy
 300     CONTINUE
         dummy=0
         DO 400 i=1,5
            WRITE(luttrf) dummy
 400     CONTINUE
      ELSE
c     WRITE(6,*)'No of char in filenm: ',lout
c     WRITE(6,'(a)')'Filename: ',filenm(1:lout)
c     PAUSE
c     OPEN(luttrf, file='prosim_out.dat',
c     &        status='unknown',form='formatted')
c     OPEN(luttrf,file='pek.asc',
c     &        status='unknown',form='formatted')
         fileid='PULSETRF'
         prognm='SNAP '
         WRITE(luttrf,'(1x,a)') fileid
         WRITE(luttrf,'(1x,a)') prognm
c     WRITE(LUTTRF,*) NOUT
         WRITE(luttrf,*) INT(1.)
         icnt=1
         WRITE(luttrf,*) INT(1.)
         WRITE(luttrf,'(1x,a)') title(1:64)
         signn='+'
         WRITE(luttrf,'(1x,a)') signn
         WRITE(luttrf,*) freqs
         WRITE(luttrf,*) sd
c     IR=1
         WRITE(luttrf,*) rd,rdlow,ir
         WRITE(luttrf,*) r0,rspace,nplots
         WRITE(luttrf,'(3I6,E20.10)') nx,lx,mx,dt1
         icdr=0
         WRITE(luttrf,*) icdr
         WRITE(luttrf,*) omegim
c     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         WRITE(luttrf,*) msuft
c     WRITE(6,*) 'trfhead: msuft=',msuft
         isrow=1
         WRITE(luttrf,*) isrow
         inttyp=1
         WRITE(luttrf,*) inttyp
         idummy=0
         DO 301 i=1,2
            WRITE(luttrf,*) idummy
 301     CONTINUE
         dummy=0e0
         DO 401 i=1,5
            WRITE(luttrf,*) dummy
 401     CONTINUE
      END IF
      RETURN
      END
