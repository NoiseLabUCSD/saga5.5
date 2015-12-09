c
      subroutine trfhead(filenm,title,rd,rdlow,r0,rspace,
     &  nx,lx,mx,dt1,freqs,sd,bintrf,ir,is,nplots)
c
      implicit none
c
      logical bintrf
      integer*4 lsvp,lopt,out,lout
      integer*4 i,idummy,inttyp,isrow,icdr,msuft,nx,lx,mx,
     &          nplots,ir,is,icnt
      real*8 dummy,adummy(10),omegim,rd,rdlow,dt1,
     &       sd,freqs
      real*8 r0,rspace
      character*7 prognm
      character*8 fileid
      character*64 filenm
      character*64 title
      character signn
c >>> Dummy real*4 to force single prec trf file.
c
      omegim=0.
      if (bintrf) then
         if(is.gt.1) stop'in trfhead'
         fileid='PULSETRF'
         write(16) fileid
         prognm='GAMARAY'
         write(16) prognm
c         nout=1
c         write(16) nout
         write(16) int(1.)
c         icnt=1
c         write(16) icnt
         write(16) int(1.)
         write(16) title(1:64)
         signn='+'
         write(16) signn
c     CENTER FREQUENCY
         dummy=freqs
         write(16) dummy
c     SOURCE DEPTH
         dummy=sd
         write(16) dummy
c     UPPER MOST RECEIVER DEPTH
         adummy(1)=rd
c     LOWER MOST RECEIVER DEPTH
         adummy(2)=rdlow
c     IR=NO OF RECEIVERS BETWEEN R0 AND RDLOW
         write(16) adummy(1),adummy(2),ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
c     do L = 1,abs(ir)
c     dummy=rdc(L)
c     WRITE(16) dummy
c     end do
c     RDC(L) CONTAINS ALL THE RECEIVER DEPTHS
C     write(16) (rdc(l),l=1,abs(ir))
C     END IF
C
C     R0= THE FIRST RANGE TO PLOT IN RANGE STACKED PLOT
         adummy(1)=r0
C     RSPACE= THE RANGE INCREMENT TO PLOT IN RANGE STACKED PLOT
         adummy(2)=rspace
c     WRITE R0, RSPACE AND THE NO OF PLOTS ASSOCIATED WITH IR
         WRITE(16) adummy(1),adummy(2),nplots
C     DT=TIME SAMPLING
         dummy=dt1
C     NX=NO OF TIME SAMPLES (DENOTED NT IN OASES MANUAL) MAYBE?
C     LX=INDEX OF FIRST FREQUENCY COMPONENT (INT(FR1*DT))
C     MX=INDEX OF LAST FREQUENCY COMPONENT (INT(FR2*DT))
         write(16) nx,lx,mx,dt1
         icdr=0
         write(16) icdr
         dummy=omegim
         write(16) dummy
C     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(16) msuft
c         write(6,*) 'trfhead: msuft=',msuft
         isrow=1
         write(16) isrow
         inttyp=1
         write(16) inttyp
         idummy=0
         do 300 i=1,2
            write(16) idummy
 300     continue
         dummy=0
         do 400 i=1,5
            write(16) dummy
 400     continue
      else
         fileid='PULSETRF'
         prognm='GAMARAY'
         write(16,'(1x,a)') fileid
         write(16,'(1x,a)') prognm
c     WRITE(16,*) NOUT
         write(16,*) int(1.)
         icnt=1
         write(16,*) int(1.)
         write(16,'(1x,a)') title(1:64)
         signn='+'
         write(16,'(1x,a)') signn
         write(16,*) freqs
         write(16,*) sd
c     IR=1
         write(16,*) rd,rdlow,ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
C     WRITE(16,*) (RDC(L),L=1,ABS(IR))
C     END IF
C     NPLOTS=1
C
         write(16,*) r0,rspace,nplots
         write(16,'(3I6,E20.10)') nx,lx,mx,dt1
         icdr=0
         write(16,*) icdr
         write(16,*) omegim
c     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(16,*) msuft
c         write(6,*) 'trfhead: msuft=',msuft
         isrow=1
         write(16,*) isrow
         inttyp=1
         write(16,*) inttyp
         idummy=0
         do 301 i=1,2
            write(16,*) idummy
 301     continue
         dummy=0e0
         do 401 i=1,5
            write(16,*) dummy
 401     continue
      end if
      return
      end
