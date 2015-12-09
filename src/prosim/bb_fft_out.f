c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************

      subroutine bb_fft_out()
      implicit none
      include 'Parms_com'
      include 'i_o_svp_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      include 'i_o_1tf_com'
      include 'i_o_2_com'
      include 'gen3_com'
      integer*4 jj0,jsrc,jrec,j,mx
      real*4 rspace,freqs,dt1,sd
      logical bintrf
      bintrf=.false.
c      bintrf=.true.
c
c     rspace = range increment for receiver position
cpln      rkm(1)=rkm(1)*1.0e-3
      if (nsrc.eq.1) then
         rspace=rkm(1)*1.0e-3
      else
         rspace=((rkm(nsrc)-rkm(1))*1.0e-3)/float(nsrc-1)
      end if
      mx=nf1+nfbb-1
      dt1=1.E0/fsbb
      freqs=(fmax-fmin)/2.+fmin
c      sd=zsr(mzsrc(1))
      sd=zsrc(1)
       call trfhead(outroot,svp_title,zrecusr(1),zrecusr(nrecusr),
     &   rkm(1)*1.0e-3,rspace,nfftbb,nf1,mx,dt1,freqs,sd,
     &   bintrf,nrecusr,nzs,nsrc)
c
      if(iifft .ne. 0) then
c
         if (bintrf) then
            do j=1,nfbb
               do jsrc=1,nsrc
                  jj0=(jsrc-1)*nfbb*nrecusr
                  do jrec=1,nrecusr
cc fbv                     write(95, *) -real(tf(j+jj0)), faxbb(j)
c Output for Normal Stress
c                     write(luttrf)-real(tf(j+jj0)), aimag(tf(j+jj0))
c                     write(luttrf)-dreal(tf(j+jj0))
c the previous line for the output on linux workstations
c Output for Pressure
                     write(luttrf)real(tf(j+jj0)), -aimag(tf(j+jj0))
                     jj0=jj0 + nfbb
                  end do
               end do
            end do
         else
            do j=1,nfbb
               do jsrc=1,nsrc
                  jj0=(jsrc-1)*nfbb*nrecusr
                  do jrec=1,nrecusr
cc fbv                     write(95, *) real(tf(j+jj0)),
cc fbv     &                            faxbb(j),jrec
c Output for Normal Stress
c                     write(luttrf,*)-real(tf(j+jj0)),aimag(tf(j+jj0))
c                     write(luttrf,*)-dreal(tf(j+jj0))
c the previous line for the output on linux workstations
c Output for Pressure
                     write(luttrf,*)real(tf(j+jj0)),-aimag(tf(j+jj0))
                     jj0=jj0 + nfbb
                  end do
               end do
            end do
         end if
c     
         close(luttrf)
      endif
c
      return
      end
c
      subroutine trfhead(filenm,title,rd,rdlow,r0,rspace,
     &  nx,lx,mx,dt1,freqs,sd,bintrf,ir,is,nplots)
c
      implicit none
c
      include 'Parms_com'
      include 'i_o_svp_com'
      common /out_com2/ lsvp,lopt,out,lout
c
      logical bintrf
      integer*4 lsvp,lopt,out,lout
      integer*4 i,idummy,inttyp,isrow,icdr,msuft,nx,lx,mx,
     &          nplots,ir,is,icnt
      real*4 dummy,adummy(10),omegim,r0,rd,rspace,rdlow,dt1,
     &       sd,freqs
      character*6 prognm
      character*8 fileid
      character*64 filenm
      character*80 title
      character signn
c >>> Dummy real*4 to force single prec trf file.
c
      omegim=0.
      title(65:80)=' '
      if (bintrf) then
         if(is.gt.1) stop'in bb_modes'
c         open(luttrf, file='prosim_out.dat',
c     &        status='unknown',form='unformatted')
         fileid='PULSETRF'
         write(luttrf) fileid
         prognm='PROSIM'
         write(luttrf) prognm
c         nout=1
c         write(luttrf) nout
         write(luttrf) int(1.)
c         icnt=1
c         write(luttrf) icnt
         write(luttrf) int(1.)
         write(luttrf) title(1:80)
         signn='+'
         write(luttrf) signn
c     CENTER FREQUENCY
         dummy=freqs
         write(luttrf) dummy
c     SOURCE DEPTH
         dummy=sd
         write(luttrf) dummy
c     UPPER MOST RECEIVER DEPTH
         adummy(1)=rd
c     LOWER MOST RECEIVER DEPTH
         adummy(2)=rdlow
c     IR=NO OF RECEIVERS BETWEEN R0 AND RDLOW
         write(luttrf) adummy(1),adummy(2),ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
c     do L = 1,abs(ir)
c     dummy=rdc(L)
c     WRITE(LUTTRF) dummy
c     end do
c     RDC(L) CONTAINS ALL THE RECEIVER DEPTHS
C     write(luttrf) (rdc(l),l=1,abs(ir))
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
         write(luttrf) nx,lx,mx,dt1
         icdr=0
         write(luttrf) icdr
         dummy=omegim
         write(luttrf) dummy
C     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(luttrf) msuft
c         write(6,*) 'trfhead: msuft=',msuft
         isrow=1
         write(luttrf) isrow
         inttyp=1
         write(luttrf) inttyp
         idummy=0
         do 300 i=1,2
            write(luttrf) idummy
 300     continue
         dummy=0
         do 400 i=1,5
            write(luttrf) dummy
 400     continue
      else
c         write(6,*)'No of char in filenm: ',lout
c         write(6,'(a)')'Filename: ',filenm(1:lout)
c         pause
c         open(luttrf, file='prosim_out.dat',
c     &        status='unknown',form='formatted')
c         open(luttrf,file='pek.asc',
c     &        status='unknown',form='formatted')
         fileid='PULSETRF'
         prognm='PROSIM'
         write(luttrf,'(1x,a)') fileid
         write(luttrf,'(1x,a)') prognm
c     WRITE(LUTTRF,*) NOUT
         write(luttrf,*) int(1.)
         icnt=1
         write(luttrf,*) int(1.)
         write(luttrf,'(1x,a)') title(1:64)
         signn='+'
         write(luttrf,'(1x,a)') signn
         write(luttrf,*) freqs
         write(luttrf,*) sd
c     IR=1
         write(luttrf,*) rd,rdlow,ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
C     WRITE(LUTTRF,*) (RDC(L),L=1,ABS(IR))
C     END IF
C     NPLOTS=1
C
         write(luttrf,*) r0,rspace,nplots
         write(luttrf,'(3I6,E20.10)') nx,lx,mx,dt1
         icdr=0
         write(luttrf,*) icdr
         write(luttrf,*) omegim
c     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(luttrf,*) msuft
c         write(6,*) 'trfhead: msuft=',msuft
         isrow=1
         write(luttrf,*) isrow
         inttyp=1
         write(luttrf,*) inttyp
         idummy=0
         do 301 i=1,2
            write(luttrf,*) idummy
 301     continue
         dummy=0e0
         do 401 i=1,5
            write(luttrf,*) dummy
 401     continue
      end if
      return
      end
