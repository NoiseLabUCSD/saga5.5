      SUBROUTINE forwardmodel(iopt,mopt)
      INTEGER  mopt,i,iopt(mopt)
      DO i=1,mopt
        iopt(i)=0
      ENDDO
      iopt(1)=2
      iopt(12)=1            ! using three indexes for adressing variable
      iopt(30)=10            ! 10 is GAMA
      END 
c
       SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization PROGRAM
c     PETER GERSTOFT, 1992
c
      USE global
      IMPLICIT INTEGER (i-n)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './gama/common/parms_conv'
      INCLUDE './gama/common/srlox'
      INCLUDE './gama/common/saga_freqcom'
      INCLUDE './gama/common/saga_caustix'
      INCLUDE './gama/common/depth'
      INCLUDE './gama/common/tilt'
c
      INTEGER i,j,jj,ndum,ierr,irflag,iicw
      REAL*8 rdstep,rderr
      CHARACTER*80 dumch, dumch2
      EQUIVALENCE (dumch2,opt)

      ierr=0
c---  READ the GA-parameters
      CALL readinputstart
      
      READ(1,'(a)') dumch2         ! This line has no effect in snap
      CALL gamaoption(opt)
      
      WRITE(*,*)'reading nfrq..'
      READ(1,*) nfrq
 92   WRITE(*,'(a,i3,a,i4)') ' nfrq=', nfrq
      IF (iopt(5).EQ.5) THEN
         BACKSPACE(1)
         WRITE(*,*)'for option F reading nfrq and time delay'
         READ(1,*) nfrq,del_time
      ENDIF
      IF ((nfrq).GT.mfreq) STOP ' nfreq > mfreq' ! test
      IF (nfrq.GT.0) THEN
         nfr=nfrq
         READ(1,*)(frqsaga(I),i=1,nfrq)
         WRITE(*,'(a,10f10.3)')' frequencies:',(frqsaga(i),i=1,nfrq)
         DO i=1,nfrq
            frq(i)=frqsaga(i)
         END DO
         iicw=1
         fminsaga=frqsaga(1)
         fmaxsaga=frqsaga(nfr)
         flo=fminsaga
         fhi=fmaxsaga
         flotmp=flo
         fhitmp=fhi
         fs=4*fmaxsaga
      ELSE
         iicw=2
         nfr=1
         frqsaga(1)=100.
 10      WRITE(*,*)'reading frequencies...'
         READ(1,'(a80)')dumch
         WRITE(*,*)dumch
         READ(dumch,*,err=11,END=11)
     &    dfrq, fminsaga, fmaxsaga,dflast,fminlast,fmaxlast
c         write(6,*)nfr,nfrq,(frqsaga(I),i=1,nfrq),dfrq, fminsaga, 
c     &        fmaxsaga,dflast,fminlast,fmaxlast
c         stop
         GOTO 12
 11    CONTINUE
c      BACKSPACE(1)
         READ(dumch,*,err=10)dfrq, fminsaga, fmaxsaga

 12      CONTINUE
         nf1=MAX(1,NINT(fminsaga/dfrq))+1
         nf2=NINT(fmaxsaga/dfrq)+1
         fminsaga=(nf1-1)*dfrq
         fmaxsaga=(nf2-1)*dfrq
         nfbb=nf2 - nf1 + 1
         nfrq=nfbb
         nfftbb=2**(NINT(LOG10(float(nf2))/LOG10(2.) + .9))*2
         nfft=nfftbb
         fs= dfloat(nfftbb)*dfrq
         WRITE(*,*)' FS, Fmin, Fmax, DF, nfrq: ', 
     .        fs, Fminsaga, Fmaxsaga, dfrq,nf1,nf2,nfftbb
         IF (nfrq.GT.mfreq) THEN
           WRITE(*,*)' mfreq,nfrq=',mfreq,nfrq
           STOP 'increase mfreq! '
         ENDIF
         flo=fminsaga
         fhi=fmaxsaga
         flotmp=flo
         fhitmp=fhi
         DO i=1,nfbb
            frq(i)=fminsaga+(i-1)*dfrq
         END DO
      ENDIF
c     BACKSPACE(1)
c     BACKSPACE(1)
      
      CALL svp_read_gama
c

      write(6,*)'PETER IN'
      write(6,*)freq,fminsaga,fmaxsaga
      write(6,*)'PETER OUT'

      nsctor=1
c
c---  for several source locations
      READ(1,*) sd         
      WRITE(*,*)' Source depth',sd
c     
c     receiver DATA
c     
      IF ((tiltv).OR.(tilth)) THEN
         READ(1,*) rd,rdlow,ndep,dtiltv
         dtilth=dtiltv
      ELSE
         READ(1,*) rd,rdlow,ndep
         dtiltv=0.0
         dtilth=0.0
      ENDIF
c
      IF (ndep.GT.mdep) THEN
         WRITE(*,*) 'ndep must grater than mdep',ndep,mdep
         ierr=1
      ENDIF
c     
c---  READ the range-BLOCK
      CALL  read_range_inp(ierr,ndum)
c---  transfer receiver parameters
      nrangx=nx
      DO i=1,nx
         rtmp(i)=DBLE(xranges(i)*0.001)
      ENDDO
c
      IF(ABS(nx).GT.1) THEN
         drtmp=rtmp(2)-rtmp(1)
      ELSE
         drtmp=1.
      END IF
      drtmp0=drtmp
c
c---  source info
      nzs=1
      srloc(1)=DBLE(sd)
c---  passing receiver depth
c
      irflag=ndep
      ndep=ABS(ndep)
      nzr=ndep
      zr1=rd
      IF (ndep.GT.1) THEN
         rdstep=DBLE(rdlow-rd)/float(ndep-1)
         dz0(1)=rdstep
         IF(tiltv) THEN
c            nrangx=ndep
c            nrec=nrangx

c            DO i=1,nrangx
c               rtmp(i)=(DBLE(xranges(1))+
c     >              DBLE((i-1)/(nrangx-1))*dtiltv)*0.001
c               xranges(i)=rtmp(i)
c            END DO

            rtmp(1)=DBLE(xranges(1))*0.001
            nrec=ndep
            nrangx=1
            drtmp=1.E0/DBLE(nrangx-1)*dtiltv*0.001           
            rdstep=dsqrt((dz0(1))**2
     >           -(ABS(dtiltv)/DBLE(ndep-1))**2)
         END IF
      ELSE
         IF(tilth) THEN
            nrec = nrangx
            nzr=nrangx
            rdstep = dtilth/float(nrec-1)
            drtmp=dsqrt((drtmp0)**2
     >            -(0.001*ABS(dtilth)/DBLE(nrec-1))**2)
            nrangx=1
         ELSE
            rdstep=1.
         END IF
      ENDIF

      IF (irflag.GT.0) THEN
c     equidistant receiver depths
         DO  jj=1,ndep
            dz(jj)=rdstep
            rdep(jj)=zr1+DBLE(jj-1)*dz(jj)
            WRITE(6,*)'RDEP: ',rdep(jj)
         ENDDO
      ELSE
         WRITE(*,*)'Only equidistant receivers supported'
         STOP
c     READ in individual receiver depths
c     READ(1,*) (rdep(jj),jj=1,ndep)
      ENDIF
      WRITE(prtfil,*)'number of receiver depths read',ndep
c     
      rderr=MAX(rdep(1),rdep(ndep))
      DO jj=1,nsctor
         IF(R_H0(jj).LE.rderr) THEN
            WRITE(6,*)'   ******Receiver in bottom*****'
            WRITE(6,*)'Revise receiver/water depth in sector:',jj
            STOP
         END IF
      END DO
c
cpln      nrec=ndep
      WRITE(*,*)' Number of depth', ndep
      WRITE(*,*)' first & last depth',rd,rdlow
      WRITE(*,*)' first & last depth',zr1,(ndep-1)*dz(1)+zr1
c---  range parameters
      
      WRITE(*,*)' Number of ranges', nx
      WRITE(*,*)' first & last range',xranges(1),xranges(nx)
      
         
c---  should the  EOF be READ ? 
      
      IF (iopt(17).EQ.1) THEN
         CALL eofinit
      ENDIF
      
c***  create text strings for physical parameters
      phystxt(1) = 'Water depth (m) $'
      phystxt(2) = 'Water speed (m/s) $'
      phystxt(3) = 'Sediment sound speed (m/s) $'
      phystxt(4) = 'Attenuation (dB/lambda) $'
      phystxt(5) = 'Surface roughness (m) $'
      phystxt(6) = 'Sediment Density (g/cm3) $'
      phystxt(8) = 'Source depth (m) $'
      phystxt(9) = 'Source range (km) $'
      phystxt(11)= 'Shape coefficient $'
      phystxt(12)= 'Bottom speed (m/s) $'
      phystxt(13)= 'Bottom S-speed (m/s) $'
      phystxt(14)= 'Bottom density (g/cm3) $'
      phystxt(15)= 'Receiver depth (m) $'
      phystxt(16)= 'Depth of water speed pt.(m) $'
      phystxt(17)= 'Sediment depth (m) $'
      phystxt(18)= 'Depth of sed. speed pt.(m) $'
      phystxt(19)= 'Vertical array tilt (m)$'
      phystxt(20)= 'Time delay $'
      phystxt(21)= 'Horizontal array tilt (m)$'
      phystxt(22)= 'SSP grd g (1/s) $'
      phystxt(23)= 'SSP curvature BETA $'
      phystxt(24)= 'Receiver separation in range (m) $'
c***  

      DO 8 j=1,mphys
         phystxt2(j)='                                               '
 6       DO 7 I=40,1,-1
            IF(phystxt(j)(I:I).NE.'$') GO TO 7
            phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c      WRITE(*,*)phystxt2(j)
            GO TO 8
 7       CONTINUE
    8 CONTINUE
      phystxt2(9)='Source range (m)'
      
c---- READ the optimization PARAMETER
      CALL readoptparm
      

      DO i=1,nparm
         IF (fmin(i).GT.fmax(i))THEN
            WRITE(*,*)' *** fmin > fmax for parm',i
            ierr=1
         ENDIF
         IF (par2phy(i).LT.1 .OR. par2phy(i).GT.24 .OR.
     &       par2phy(i).EQ.5 .OR.
     &       par2phy(i).EQ.7 .OR.
     &       par2phy(i).EQ.10 .OR.
     &       par2phy(i).EQ.13) THEN
cpln     &       par2phy(i).EQ.13 .OR.
cpln     &       par2phy(i).EQ.19 )THEN
            WRITE(*,*)' *** par2phy not correct, parm',i
            ierr=1
        ENDIF
        IF (par2phy(i).EQ.1 )THEN
c     IF (fmin(i).LT.R_z0(R_nd0(par3(i))-1,par3(i))) THEN
c          WRITE(*,*)' *** Optimzation variable #:',i
c     WRITE(*,*)' *** Durring optimization'
c     WRITE(*,*)' *** the water-depth can get above the second-last'
c          WRITE(*,*)' *** point in the water'
c          ierr=1
c          ENDIF
        ENDIF
        IF (par2phy(i).EQ.11 )THEN
           IF (par2lay(i).GT.neof) THEN
              WRITE(*,*)' *** Optimzation variable #:',i
              WRITE(*,*)' *** The shapecoffient number is not defined'
              WRITE(*,*)' par2lay(i)', par2lay(i)
              ierr=1
           ENDIF
        ENDIF
        IF (par2phy(i).EQ.17)THEN
c     IF (fmin(i).LT.R_z1(R_nd1(par3(i))-1,par3(i))) THEN
c          WRITE(*,*)' *** Optimzation variable #:',i
c          WRITE(*,*)' *** During optimization'
c          WRITE(*,*)' *** the sediment-depth can get above the'
c          WRITE(*,*)' *** second-last point in the sediment'
c          ierr=1
c          ENDIF
        ENDIF
      ENDDO
      WRITE(*,*)' Number of data points',nfrq*nx*ndep
      WRITE(prtfil,*)' Number of data points',nfrq*nx*ndep
      IF (nfrq*nx*ndep.GT. mobs) THEN
        WRITE(*,*)' Number of data points',nfrq*nx*ndep
        WRITE(*,*)' Number of data point is too large to be contained'
        WRITE(*,*)' in data matrix, with dimension mobs=',mobs
        ierr=1
      ENDIF

c*** errors ?
      IF (ierr.EQ.1)STOP 'stopped in input'
c***  CLOSE input unit
      CLOSE(1)                                
      END
c**************************************************************
      SUBROUTINE gamaoption(options)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './gama/common/tilt'
c      INCLUDE 'comsnaprd.h'
c      INCLUDE 'comprosim.h'

      CHARACTER options(40)
      INTEGER i
      INTEGER ipos,iposV
      COMMON /div/ipos,iposV
c
      ipos=0
      tilth=.FALSE.
      tiltv=.FALSE.
      hilbt=.FALSE.
c      incoh=.FALSE.
      DO 50 i=1,40
         IF (options(I).EQ.'z') THEN
            ipos=2
         END IF
         IF (options(I).EQ.'V') THEN
            iposV=3
         END IF
         IF (options(I).EQ.'p') THEN
            ipos=1
            WRITE(prtfil,*)
     &           ' Only positive gradients in sediment and bottom'
         ELSE IF (options(i).EQ.'t') THEN
            tiltv=.TRUE.
cpln            tiltv=.FALSE.
            WRITE(*,*)' tilt of vertical array is included'
            WRITE(prtfil,*)' tilt of vertical array is included'
         ELSE IF (options(i).EQ.'T') THEN
            tilth=.TRUE.
            WRITE(*,*)' tilt of horizontal array is included'
            WRITE(prtfil,*)' tilt of horizontal array is included'
         ELSE IF (options(i).EQ.'h') THEN
            hilbt=.TRUE.
            WRITE(*,*)' Envelope of impulse response'
            WRITE(prtfil,*)' Envelope of impulse response'
         ELSE IF (options(I).EQ.'!') THEN
            GOTO 60
         ELSE IF (options(I).NE.' ') THEN
            WRITE(prtfil,399) options(I)
 399        FORMAT(1H ,' >>>> unknown GAMMARAY OPTION: ',A1,' <<<<')
         END IF
 50   CONTINUE
 60   CONTINUE
      IF((tiltv).AND.(tilth)) THEN
         WRITE(6,*)'Horizontal and vertical tilt together ',
     .             ' is not allowed'
         STOP
      END IF
      END
c
      SUBROUTINE bb_fft_out(rangx,tfz)
c
      INCLUDE './gama/common/depth'
      INCLUDE './gama/common/srlox'
      INCLUDE './gama/common/saga_freqcom'
c
      INTEGER*4 luttrf
      REAL*8 rangx(nrangx),dt1,freqs,zrnz
      COMPLEX*8 tfz(nrangx*nrec*nfbb)
      CHARACTER*64 outfile,outroot
       CHARACTER*80 svp_title
      LOGICAL bintrf
      integer j,jj0,jrec,jsrc,lout,mx1
      write(*,*)'enter gama bbfftout'
c
c      bintrf=.TRUE.
      bintrf=.FALSE.
      luttrf=16
c
c      WRITE(6,*)'fh,fl: ',fhitmp,flotmp
      freqs=(fhitmp+flotmp)/2
      bintrf=.FALSE.
      dt1=1/fs
cpln            nf1=MAX(1,NINT(flotmp/dfrq)) + 1
cpln            nf2=NINT(fhitmp/dfrq) + 1
      mx1=nf1+nffth1-1
      zrnz=zr1+float(nzr-1)*dz(1)
c      WRITE(6,*)'freqs,fs,fmin,zr1,dz: ',
c     &     freqs,fs,flotmp,zr1,dz(1)
c      WRITE(6,*)'rangx(1),dr: ',
c     &     rangx(1)/1000.,(rangx(2)-rangx(1))/1000.
      lout=11
      outroot(1:lout)='gamaout'
      outfile=outroot(1:lout)//'.trf'
      svp_title='TRF OUTPUT FROM SAGA'
c      WRITE(6,*)'1st rec depth: ',zrec(1)
c      WRITE(6,*)'last rec depth: ',zrec(nrec)
c      WRITE(6,*)'Delta time: ',dt1
c      WRITE(6,*)'Freq components: ',mx
c      WRITE(6,*)'Center frequency: ',freqs
c      WRITE(6,*)'Source depth: ',sd
c      WRITE(6,*)'Range space: ',drtmp
c      WRITE(6,*)'Logical bintrf ',bintrf
c      WRITE(6,*)'Before call to trfhead'
c      PAUSE
c        write(*,*)'Input',outroot(1:lout),svp_title,zr1,zrnz,
c     &     rangx(1)/1000.,drtmp,
c     &     nfft,nf1,nf2,dt1,freqs,zs,
cpln     & nfft,1,nffth1,dt1,freqs,zs,
c     &     bintrf,nrec,nzs,nrangx
c       write(*,*)'before trfhead'
c      CALL trfhead(outroot(1:lout),svp_title,zr1,zrnz,
      CALL trfhead(outroot,svp_title,zr1,zrnz,
     &     rangx(1)/1000.,drtmp,
     &     nfft,nf1,nf2,dt1,freqs,zs,
cpln     & nfft,1,nffth1,dt1,freqs,zs,
     &     bintrf,nrec,nzs,nrangx)
c       Write(*,*)'after trfhead'
      IF(iifft .NE. 0) THEN
c
         IF (bintrf) THEN
            DO j=1,nfbb
               DO jsrc=1,nrangx
                  jj0=(jsrc-1)*nfbb*nrec
                  DO jrec=1,nrec
                     WRITE(luttrf)REAL(tfz(j+jj0)),AIMAG(tfz(j+jj0))
                     jj0=jj0 + nfbb
                  END DO
               END DO
            END DO
         ELSE
            DO j=1,nfbb
               DO jsrc=1,nrangx
                  jj0=(jsrc-1)*nfbb*nrec
                  DO jrec=1,nrec
                     WRITE(luttrf,*)REAL(tfz(j+jj0)),AIMAG(tfz(j+jj0))
                     jj0=jj0 + nfbb
                  END DO
               END DO
            END DO
         END IF
c     
          CLOSE(luttrf)
      ENDIF
c
      RETURN
      END
c
      SUBROUTINE saga_out(tfz)
c
c: *******************************
c: *     REVISION (1996):  
c************************************
c Developed exclusively for the USE by SAGA 
c      *
c: *     P R O S I M  GROUP      *
c: *     S A C L A N T C E N     *
c: *******************************
c
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './gama/common/depth'
      INCLUDE './gama/common/srlox'
      INCLUDE './gama/common/saga_freqcom'
      INCLUDE './gama/common/tilt'
c
      COMPLEX*8 tfz(nrangx*nrec*nfbb),tfh1(nfft),tfh2(nfft),
     .          zzero
      CHARACTER*64 svp_title,outfile,outroot
      LOGICAL bintrf, dohilb
      INTEGER jj0,jj1,jsrc,jrec,j,index, ntap1
      REAL pie, fact, ptap, dpi, cshade(nfbb),summ
      DATA zzero/(0.,0.)/
      integer n1,n2
c
      pie=4.0*ATAN(1.0)
      
c      WRITE(6,*)'#####################################'
c      WRITE(6,*)
c      WRITE(6,*)'SAGA_OUT Link from GAMA to SAGA'
c      WRITE(6,*)'NFFT: ',nfft
c      WRITE(6,*)'nf1,nf2: ',nf1,nf2
c      WRITE(6,*)'nrangx,nrec,nfbb: ',nrangx,nrec,nfbb
c      WRITE(6,*)
c      WRITE(6,*)'#####################################'
c
      ptap=50.0
      ntap1=INT(nfbb*ptap/100)
      dpi=pie/(ntap1-1)
      n1=ntap1+1
      n2=nfbb-ntap1
      DO j=1,nfbb
         cshade(j)=1.0
      END DO
      DO j=1,n1-1
         cshade(j)=0.5*(1-COS((j-1)*dpi))
         cshade(nfbb-j+1)=cshade(j)
      END DO
c
      IF(hilbt) THEN
          DO jsrc=1,nrangx       !  range
            DO jrec=1,nrec      !depth
               jj0=(jsrc-1)*nfbb*nrec
               DO j=1,nf1-1
                  tfh1(j)=zzero
               END DO
               DO j=1,nfbb
c                  tfh1(nf1-1+j)=-2.0*(CONJG(tfz(j+jj0)))
                  tfh1(nf1-1+j)=-2.0*(CONJG(tfz(j+jj0)))*cshade(j)
c
c                  WRITE(74,*)REAL(-2.0*(CONJG(tfz(j+jj0)))),
c     >            imag(-2.0*(CONJG(tfz(j+jj0))))
c                  WRITE(74,*)REAL(tfh1(nf1-1+j)),
c     >            imag(tfh1(nf1-1+j))
c                  WRITE(99,*)cshade(j)
c
               END DO
c               CLOSE(74)
c               CLOSE(99)
c               PAUSE

               DO j=nfbb+1,nfft-nf1+1
                  tfh1(nf1-1+j)=zzero
               END DO

c               DO j=1,nfft
c                  WRITE(75,*)REAL(tfh1(j)),imag(tfh1(j))
c               END DO

               CALL four1(tfh1,nfft,-1)

c               DO j=1,nfft
c                  WRITE(76,*)REAL(tfh1(j)),imag(tfh1(j))
c               END DO

               fact=1.0/REAL(nfft)
               summ=0.
               DO j=1,nfft
c                  tfh2(j)=REAL(fact*tfh1(j))**2+imag(fact*tfh1(j))**2
                  tfh2(j)=fact*SQRT(REAL(tfh1(j))**2+imag(tfh1(j))**2)
                  summ=summ+tfh2(j)
               END DO
               summ=summ*fact
c               WRITE(6,*)'AVG: ',summ

               DO j=1,nfft
                  tfh2(j)=tfh2(j)-summ
               END DO

c               DO j=1,nfft
c                  WRITE(76,*)REAL(tfh2(j)),imag(tfh2(j))
c               END DO

               CALL four1(tfh2,nfft,1)

c
c               DO j=1,nfft
c                  WRITE(77,*)REAL(tfh2(j)),imag(tfh2(j))
c               END DO
c               CLOSE(74)
c               CLOSE(75)
c               CLOSE(76)
c               CLOSE(77)

               DO j=1,nfbb      ! frequency
                  index=(jrec +((j-1))*nrec-1)*nrangx
                  resp(jsrc+index)= tfh2(j)
c                  WRITE(77,*)j,jsrc,jrec,jj0,index
                  IF (iopt(5).EQ.5) THEN
                     resp(jsrc+index)=resp(jsrc+index)
     .                    *cexp(CMPLX(0.,frq(j)*2.*pie*del_time))
                  ENDIF
               END DO
               jj0=jj0 + nfbb
            END DO
         END DO
      ELSE
         DO j=1,nfbb            ! frequency
            DO jsrc=1,nrangx    !  range
               jj0=(jsrc-1)*nfbb*nrec
               DO jrec=1,nrec   !depth
                  index=(jrec +((j-1))*nrec-1)*nrangx
                  resp(jsrc+index)= -(CONJG(tfz(j+jj0)))
c     pln               resp(jsrc+index)= -(tfz(j+jj0))
c                  WRITE(77,*)j,jsrc,jrec,jj0,index
                  IF (iopt(5).EQ.5) THEN
                     resp(jsrc+index)=resp(jsrc+index)
     .                    *cexp(CMPLX(0.,frq(j)*2.*pie*del_time))
                  ENDIF
                  jj0=jj0 + nfbb
               END DO
            END DO
         END DO
      END IF
c     
c Test on Hilbert transform



      lwascomp=1
      RETURN
      END
