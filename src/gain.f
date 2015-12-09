c******************************************************
      SUBROUTINE readdata(ierr) 
c     Reads DATA in Range FORMAT
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER i,j,index,ierr
      CHARACTER*250 dumch
      REAL rng,unr(1000) 
c
c---- First READ the file to find number of points (nx).
c
      i=0
 3    READ(2,'(a250)',END=2)dumch
      IF (dumch(1:1).NE.'!') THEN
         BACKSPACE(2)
 1       READ(2,*,END=2)rng,(unr(j),j=1,ncurv)
      ELSE 
         GOTO 3
      ENDIF
c      WRITE(*,*)'rng,rng1,rng2',rng,rng1,rng2
      IF(rng.LT.rng1)GOTO 3
      IF(rng.GT.rng2)GOTO 3
      i=i+1
      GOTO 3
    2 CONTINUE
      nx=i
      IF (nx*ncurv.GT.mobs) THEN
         WRITE(*,*)' nx,ncurv,mobs=',nx,ncurv,mobs
         WRITE(*,*)' Mobs (in comopth.h) is not large enough'
         ierr=1
      ENDIF
      IF (nx.EQ.0) THEN
         WRITE(*,*)' ******** Nx is zero'
         ierr=1
      ENDIF
      WRITE(*,*)' reading ',nx, ' data points'
      deallocate(data)
      allocate(data(ndep*nfrq*nx))
c
c--- Now READ the DATA
c
      REWIND(2)
      i=0
 13   READ(2,'(a250)',END=12)dumch
      IF (dumch(1:1).NE.'!') THEN
         BACKSPACE(2)
         READ(2,*,END=12)rng,(unr(j),j=1,ncurv)
      ELSE 
         GOTO 13
      ENDIF
      IF(rng.LT.rng1)GOTO 13
      IF(rng.GT.rng2)GOTO 13
c     rplus=rng+drng 
      i=i+1
      write(*,*)ndep,nfrq,nx,i
      IF (i.gt.(ndep*nfrq*nx))  then
         write(*,*)'Realdata: size of data-array exceeded.'
         write(*,*)'Realdata: increse nx in the input data file'
         stop
c     deallocate(data)
c     mobs=mobs+1
c     allocate(data(ndep*nfrq*nx))
      endif
      IF (i.GT.mx) THEN 
         WRITE(*,*)' Error durring reading of the real data'
         WRITE(*,*)' Readdata: Mx is not large enough for all ranges'
         WRITE(*,*)' increase Mx, mx=',mx
         STOP
      ENDIF
      xranges(i)=rng
      DO j=1,ncurv
         index=(j-1)*nx
         DATA(i+index)=unr(j)
         WRITE(*,*)'gain,data',index,i, DATA(i+index)
      ENDDO
      GOTO 13
 12   CONTINUE
      
      WRITE(*,*)' last data points:'
      WRITE(*,*)rng,(unr(j),j=1,ncurv)
      WRITE(*,*)'read ',nx,'data points'
      END

c*************************************************************
      SUBROUTINE readdat2(ierrdat2)

C     Reads hydrophone fourier vector DATA from the IN file.

      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      INTEGER  ierrdat2,i,j,index,ifrq,idep
      REAL dum,xhelp
      CHARACTER*80 dumch
      
      DO ifrq=1,nfrq
         DO idep=1,ndep
c     DO j=1,ncurv
            j=(ifrq-1)*ndep +idep
            WRITE(*,'(a,i3,a,i3)')
     >           ' Reading data for freq ',ifrq,' depth ',idep
            DO i=1,nx
 3             READ(2,'(a80)',END=2,err=2)dumch
               IF (iopt(6).EQ.1) WRITE(*,*)'reading line:',dumch
               IF (dumch(1:1).NE.'!') THEN
                  index=(j-1)*nx
                  BACKSPACE(2)
                  READ(2,*,err=31)dum,dum,dum,DATA(i+index)
                  GOTO 32
 31               CONTINUE
                  BACKSPACE(2)
                  READ(2,*,err=2)dum,dum,dum,xhelp
                  DATA(i+index)=xhelp
 32               CONTINUE
               ELSE 
                  GOTO 3
               ENDIF
            ENDDO               !nx loop
            WRITE(*,*)' data for point',index+nx,DATA(nx+index)
         ENDDO                  !ndep
      ENDDO                     !nfrq
    2 CONTINUE

c this works for reverberation DATA
      WRITE(*,*)'read ',(i-1)*(ifrq-1)*(idep-1),'data points'
 4    CONTINUE    
      IF ((ifrq-1)*(idep-1).NE.nfrq*ndep) THEN
c     IF ((j-1).NE.ncurv) THEN
         WRITE(*,'(a,i3,a,i3,a,i3,a,i3)')' ifrq',ifrq-1,
     >        ' idep',idep-1,' nfrq',nfrq, ' ndep',ndep
         WRITE(*,*) 'readdat2: not enough data points'
         ierrdat2=1
      ENDIF
      IF ((i-1).NE.nx) THEN
         WRITE(*,*)' i',i,' nx',nx
         WRITE(*,*) 'readdat2: not enough data points'
         WRITE(*,*)nx,' points should be  read in but only'
         WRITE(*,*)i,' points was read'
         ierrdat2=2
      ENDIF

      RETURN
      END
c**************************************************************
      SUBROUTINE readdatiso 
c       Dummy routine the REAL SUBROUTINE is in REV     
      END        
c
c**************************************************************
      SUBROUTINE READ_COV(ierrcov)
C     Reads DATA from the COV file which CONTAINS either REAL or
C     synthetic covariance matrices.

      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      INTEGER idep,idum,idep1,ibart,i,ndep_ch,ifrq,irng,ierrcov
      INTEGER iwarn,index2
      REAL rdum,xfreq
      REAL*8 xsum,xnorm,znoise
      CHARACTER*80 dumch
c     WRITE(*,*)'entering  read_cov'
      xcov_scale=0
      xcov_sum=0
      iwarn=0
      ibart=0
      DO ifrq=1,nfrq
         DO irng=1,nx
            ibart=ibart+1
 999        DO i=1,3 ! there should be 3 title lines. (Not including 'Format:...')
 2             READ(2,'(a80)',END=991) dumch
               IF (dumch(1:1).NE.'!') THEN
                  IF (i.EQ.2) THEN
                     BACKSPACE(2)
                     READ(2,*,END=991) xfreq
                  ELSEIF (i.EQ.3) THEN
                     BACKSPACE(2)
                     READ(2,*,END=991) ndep_ch
                     IF (ndep_ch.NE.ndep) THEN
                        WRITE(*,*)' *****Error reading cov mat******'
                        WRITE(*,*)' read ndep:', ndep_ch
                        WRITE(*,*)' req. ndep:', ndep
                        WRITE(*,*)' *****Error reading cov mat******'
                        Ierrcov=1
                        RETURN
                     ENDIF
                  ENDIF
               ELSE 
                  GOTO 2
               ENDIF
            ENDDO
c
c       Reading receiver depths
            DO idep=1,ndep
 3             READ(2,'(a80)',END=991)dumch
               IF (dumch(1:1).NE.'!') THEN
                  BACKSPACE(2)
                  READ(2,*) rdum !rd
               ELSE 
                  GOTO 3
               ENDIF
               IF (iwarn.EQ.0) THEN
                  IF (idep.EQ.1) THEN
                     IF (rdum.NE.rdep(1)) THEN
                        WRITE(*,*)'*** WARNING*** 1st rec 
     .                       depth in covmat',rdum
                        WRITE(*,*)' is different from 
     .                       the input specified',rdep(1)
                        iwarn=1
                     ENDIF
                  ELSEIF (idep.EQ.ndep) THEN
c     WRITE(*,*)'ndep=',ndep
                     IF (rdum.NE.rdep(ndep)) THEN
                        WRITE(*,*)'*** WARNING*** last rec
     .                       depth in covmat',rdum
                        WRITE(*,*)' is different from the
     .                       input specified', rdep(ndep)
                        iwarn=1
                     ENDIF  
                  ENDIF  
               ENDIF
            ENDDO
            
c >>>>> Here comes the covariance matrix
            DO idep1=1,ndep
               DO idep=1,ndep
 4                READ(2,'(a80)',END=992)dumch
                  IF (dumch(1:1).NE.'!') THEN
                     BACKSPACE(2)
                     index2=(idep+(ibart-1)*ndep-1)*ndep 
                     READ(2,*)idum,idum,cov(idep1+index2) ! cov
                  ELSE 
                     GOTO 4
                  ENDIF
               ENDDO
            ENDDO               ! loop over ndep
c-----
            IF (xfreq.LE.frq(ifrq)+0.05 .AND. 
     &           xfreq.GE.frq(ifrq)-0.05) THEN
               WRITE(*,*)' using    cov mat for freq',xfreq
            ELSE
               WRITE(*,*)' skipping cov mat for freq:',xfreq
c     !               ' Req. freq: ', frq(ifrq)
               GOTO 999
            ENDIF
c----   find the trace of covariance matrix
            xsum=0
            DO idep1=1,ndep
               index2=(idep1+(ibart-1)*ndep-1)*ndep 
               idep=idep1
               xsum=xsum+cov(idep+index2) 
            ENDDO
            xcov_trace(ibart)=xsum
c-----   Find largest eigenvalue
            IF (ndep .GT. 1000) THEN
               WRITE(*,*)' The covariance matrix is large,',
     &              ' thus the computation of the SNR has been skipped'
            ELSE
               CALL svdmax(ibart,xsum,znoise)
               zcov_noise(ibart)=znoise
               xcov_lamb(ibart)=xsum
c          WRITE(*,*)' trace and the lagest eigenvalue of cov-matrix:'
c          WRITE(*,*) xcov_trace(ibart),xcov_lamb(ibart)
               WRITE(*,'(a,f10.2)')' estimated SNR from cov matrix:', 
     &              10*LOG10((xcov_lamb(ibart)-znoise)/(znoise))
               IF (xcov_trace(ibart).LT.0) THEN 
                  WRITE(*,*)' *** Error in covariance matrix: ',
     &                 'The trace of the covariance matrix 
     &                 is negative! '
                  STOP
               ENDIF
            ENDIF
c
c-----  normalize covarince matrix WITH diagonal or max eigenvalue
c

            IF (iopt(20).GE. 1) THEN
               IF (iopt(20).EQ.1   ) THEN 
                  xnorm=xcov_trace(ibart)
                  xcov_lamb(ibart)=xcov_lamb(ibart)/xnorm
                  xcov_trace(ibart)=1
               ELSEIF (iopt(20).EQ.2   ) THEN 
                  xnorm=xcov_lamb(ibart)
                  xcov_trace(ibart)=xcov_trace(ibart)/xnorm
                  xcov_lamb(ibart)=1
               ELSE
                  STOP ' iopt(20)>= 2 not defined'
               ENDIF
C     .cfh.
               IF (iopt(35) .NE. 1) THEN
                  WRITE(*,*)' The covariance matrix is renormalized'
                  DO idep1=1,ndep
                     index2=(idep1+(ibart-1)*ndep-1)*ndep 
                     DO idep=1,ndep
                        cov(idep+index2)=cov(idep+index2)/xnorm
                     ENDDO
                  ENDDO
                  WRITE(*,*)'...renormalized with xnorm ', xnorm
               ENDIF
            ENDIF
c---  Find a scaling factor to be used when multipling bartlett power.
c     WRITE(*,*) xcov_trace(ibart),xcov_lamb(ibart)
            xcov_scale = xcov_scale+
     &           10.*LOG10(xcov_lamb(ibart))
c     &       10.*LOG10(xcov_trace(ibart)-xcov_lamb(ibart))
c     WRITE(*,*)'xcov_scale',xcov_scale
            xcov_sum=xcov_sum+xcov_trace(ibart)
         ENDDO                  ! loop over bartlett      
      ENDDO                     ! loop over bartlett      
      WRITE(*,*)' readcov, cov', cov(1), cov(nbart*ndep*ndep )

      RETURN
991   WRITE(*,*)' ***Problems reading COV FILE *************'  
      ierrcov=2
      RETURN
992   WRITE(*,*)' *** EOF encountered on COV FILE ****'
      ierrcov=3
      RETURN
      END


c********************************************************************
      SUBROUTINE WRITE_COV(ibart,ix)

C     Writes a covariance matrix to the OBS file. Writes the matrix
C     for ONLY a single frequency. IF multiple frequencies are to be
C     processed, THEN USE the hydrophone vector FORMAT instead of the
C     covariance matrix FORMAT.
C     The covariance matrices are normalized before writing to file.
     
      USE global
      INTEGER   OBSFIL,ibart
      INTEGER I, J,index2,index,index3,ix,j1start

      PARAMETER (  OBSFIL = 30)
      REAL      avgsignal,noise
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

C     --- Zero out covariance matrix ---
c     DO I = 1, Ndep
c     index2=(i+(ibart-1)*ndep-1)*ndep 
c     DO J = 1, Ndep
c     COV(  J +index2) = 0.0
c     END DO
c     END DO
      j1start=((ibart-1))*ndep

C     --- Compute the normalizing value ---
      AVGSIGNAL = 0.0
      DO I = 1, Ndep            ! loop over receiver depths
         index=ix+(i+j1start-1)*nx
         AVGSIGNAL = AVGSIGNAL + ABS( resp(index ) ) ** 2
c     AVGSIGNAL = AVGSIGNAL + ABS( resp(1+ (I-1)*1 ) ) ** 2
      END DO                    ! next receiver depth

      WRITE(*,*)'Sum sqr signal',avgsignal 

C     --- Compute cross-sensor correlation matrix R
      DO I = 1, Ndep
         DO J = 1, Ndep
            index=ix+(i+j1start-1)*nx
            index3=ix+(j+j1start-1)*nx
            index2=(j+(ibart-1)*ndep-1)*ndep 
            COV(i+index2 ) =    !COV( i+index2 ) +
     &           resp(index ) * CONJG( resp(index3 ) )
         END DO
      END DO

      IF (iopt(11).EQ.1) THEN
         WRITE(*,*)' Adding white noise to the data in the *.obs file'
         WRITE(*,'(a,f10.3,a)')'   with a  SNR =',snr_db,' dB'
         WRITE(*,*)
         noise=avgsignal*10**(-snr_db/10)
         WRITE(*,*)'Noise per hydrophone: ',  noise
         DO J = 1, Ndep
            i=j
            index2=(j+(ibart-1)*ndep-1)*ndep 
            COV(i+index2 ) =  COV( i+index2 ) + noise
         END DO
      ENDIF
      
C     *** WRITE normalized cross-sensor correlation matrix column-wise

      WRITE( OBSFIL, * ) TITLE(1:60)
      WRITE( OBSFIL, * ) FRQ(ibart)
      WRITE( OBSFIL, * ) Ndep, Nfrq ! Note: Nfrq is optional.
      DO I = 1, Ndep
         WRITE( OBSFIL, * ) rdep(i)
      END DO

      DO I = 1, Ndep
         DO J = 1, Ndep
            index2=(j+(ibart-1)*ndep-1)*ndep 
            WRITE( OBSFIL, 1000 ) I, J, COV( i+index2 )
 1000       FORMAT( I5, I5, ' (', E16.8, ',', E16.8, ' )' )
         END DO
      END DO
      END


c**************************************************************
      SUBROUTINE READ_HP(ierr)
      
C     Reads DATA from the IN file which CONTAINS either REAL or
C     synthetic hydrophone fourier vectors.
      
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER m_iphone
      PARAMETER(m_iphone=mdep*100)
      INTEGER idep, idum, i, ndep_ch, ifreq, ifflag, ierr
      REAL rdum, xfreq
      CHARACTER*80 dumch ! Used to pre-read each file record.
      INTEGER iphone_arr(m_iphone),iphone,idep_ch,index1,iran,index2
      INTEGER izero
      PARAMETER (izero=0)
C     READ the first line which should indicate the DATA FORMAT.

      READ(2,100) dumch
      WRITE(*,100)dumch
 100  FORMAT( a60 )             ! Format to read one record from the IN file.

C     Check to see that the FORMAT is "Hydrophone"
      IF ( INDEX(dumch,'Hydrophone') .EQ. 0 ) THEN
c         WRITE(*,*)' *** Warning *** The IN file format line does',
c     >        ' not contain "Hydrophone".'
c         WRITE(*,*)'Continuing anyway, with no error code introduced.'
         BACKSPACE(2)
      ENDIF

c
c
c
      DO 1100 iran=1,nx      
      WRITE(*,*)' ******* Range number: ',iran,' *****'
      ifreq = 1
      ifflag = 1
      DO WHILE (.TRUE.)         ! Cycle through the frequencies.
C     There should normally be 3 important lines following the FORMAT line,
C     or preceeding the DATA for the next frequency.
C     They should be: 1) Title, 2) Frequency, 3) Number of depths.
C     Note we don't do anything with the title here.
      
         DO i=1,3
 2          READ(2,100,end=991)dumch ! Move current record into dumch.
            if (iopt(6).eq.1 ) then
               write(*,*)dumch
            endif 
            IF (dumch(1:1).ne.'!') then ! check to see if it is important.
               IF (i.EQ.2) THEN ! dumch should contain the current frequency.
                  BACKSPACE(2)
                  READ(2,*,END=991) xfreq
C     Check to see that the frequency matches the specified one:
c       WRITE(*,*)'freqs', xfreq,frq(ifrq)
                  IF ( ABS( xfreq-frq(ifreq) ) .LT. 0.05 ) THEN
                     IFFLAG = 1
                  ELSE
                     ifflag = 0
                  ENDIF
               ELSEIF (i.EQ.3) THEN ! dumch should contain # of depths.
                  BACKSPACE(2)
                  READ(2,*,END=991) ndep_ch
                  IF (iopt(6).EQ.1 ) THEN
                     WRITE(*,*)' Number of depths in IN-file:',ndep_ch
                  ENDIF 
                  IF (m_iphone.LT.ndep_ch) STOP 
     1                 'increase dimension of ndep'
                  IF ((ndep_ch.NE.ndep).AND.
     1                 (ifflag.EQ.1).AND.(ifreq.EQ.1)) THEN
                     WRITE(*,*)'*** Warning ***.  Number of
     1                    depths mismatch.'
                     WRITE(*,*)'Number of depths reported in IN file:',
     1                    ndep_ch
                     WRITE(*,*)'Number of depths specified in DAT 
     1                    file:',ndep 
c     ierr=1
c     WRITE(*,*)'Error code set. Returning to main program.' 
c     RETURN
                  ENDIF
               ENDIF
            ELSE 
               GOTO 2           ! The line read was just a comment so get next line.
            ENDIF
         ENDDO                  ! For reading the 3 lines at the top.
c      WRITE(*,*) 'the first 3 lines read'
      
C     The next few lines should contain the actual channel depths.
C     Note there are 1 depths in each record.
         iphone=1
         DO i = 1, ndep_ch
c     WRITE(*,*)'idep,mdep',i,mdep
 3          READ(2,100,END=991)dumch ! move current record into dumch.
            IF (dumch(1:1).NE.'!') THEN ! check to see if it was just a comment.
c     IF ((ifreq.EQ.1).AND.(ifflag.EQ.1)) THEN ! write out depths once.
c     WRITE(*,100) dumch ! write the records with the depths to screen.
c     ENDIF
c     BACKSPACE(2)
               READ(dumch,*) rdum ! store the first depth in the rec.
c     READ(2,*) rdum      ! store the first depth in the rec.

c
c--- pg 000503 for now all are READ in !
c
c            IF ((iphone.LE.ndep).AND.
c     &           (ABS(rdum-rdep(iphone)).LT.0.005)) THEN
               iphone_arr(i)=1
               iphone=iphone+1
c            ELSE
c               iphone_arr(i)=izero
c            ENDIF
c--- pg 000503 for now all are READ in !


               IF ( (i.EQ.1).AND.(ifflag.EQ.1).AND.
     1              (ifreq.EQ.1).AND.(rdum.NE.rdep(1)) ) THEN
                  WRITE(*,*)'*** WARNING*** 1st rec depth in IN file:',
     1                 rdum
                  WRITE(*,*)'is different than the depth in DAT file:',
     1                 rdep(1)
                  WRITE(*,*)'No error code is set. Continuing.'
               ENDIF
            ELSE 
               GOTO 3           ! record was commented out so get next one.
            ENDIF
         ENDDO                  ! i -loop
c      WRITE(*,*)'finished reading depths'
         IF ((iphone-1).NE.ndep) THEN
            WRITE(*,*)' *****'
            WRITE(*,*) 'number of read phones from in file',
     1           iphone-1
            WRITE(*,*) 'number of requested  phones from dat file',
     1           ndep
            STOP 'phone mismatch between *.dat and *.in'
         ELSE 
            IF (ifreq.EQ.1) WRITE(*,*)' Found', iphone-1,
     1           ' phones that are matching in depth'
         ENDIF
         
C     Now we should be positioned at top of hydrophone DATA. These DATA
C     records should consist of the channel number and the COMPLEX
C     pressure for the current frequency.

         IF ((nfrq.LT.10).OR. (iopt(6).EQ.1)) THEN
            IF (ifflag.EQ.1) THEN
               WRITE(*,*)' Reading data for frequency:  ',xfreq
            ELSE
               WRITE(*,*)' Skipping data for frequency: ',xfreq
            ENDIF
         ELSE
            IF (ifflag.EQ.1) THEN
               WRITE(prtfil,*)' Reading data for frequency:  ',xfreq
            ELSE
               WRITE(prtfil,*)' Skipping data for frequency: ',xfreq
            ENDIF
         ENDIF
c            WRITE(*,*)' Reading data for frequency:  ',xfreq
         idep=0
         DO idep_ch=1,ndep_ch   ! cycle through depths at current freq.
 4          READ(2,100,END=991)dumch
            IF (dumch(1:1).NE.'!') THEN
               IF (ifflag.EQ.1) THEN ! Move stuff into data if valid freq.
                  IF (iphone_arr(idep_ch).EQ.1 ) THEN
                     idep=idep+1
c     BACKSPACE(2)
                     index2=(ifreq-1)*ndep+idep
                     index1=(index2-1)*nx
                     READ(dumch,*)idum, DATA( iran + index1 )
c     READ(2,*)idum, DATA( iran + index1 )
c     READ(2,*)idum, DATA( 1 +((ifreq-1)*Ndep+idep-1) )
c     WRITE(*,*)'reading data, idep,idep_ch',idep,idep_ch
                     IF (idum.NE.idep_ch) THEN
                        WRITE(*,*)' *** Warning ***',
     >                       ' Depth index mismatch at idep_ch = ',
     >                       idep_ch
                        WRITE(*,*)
     >                       'Error  set. Returning to main program'
                        ierr=2
                        RETURN
                     ENDIF
                  ENDIF
               ENDIF
            ELSE 
               GOTO 4           ! Record was commented out so try next one.
            ENDIF
         ENDDO
c      WRITE(*,*)' the data was read'  
         IF (ifflag.EQ.1) THEN  ! Check to see if ifreq should be incremented.
            ifreq = ifreq + 1
            IF (ifreq.GT.nfrq) GOTO 900 ! All required freqs read, so finish.
         ENDIF
         
      ENDDO                     ! for main do-while in frequency
      
 900  CONTINUE
C     Should have finished reading DATA for the current frequency.
      WRITE(*,*)' All',ifreq-1,' frequencies was read from In-file'
 1100 CONTINUE                  ! range loop

      RETURN
      
 991  WRITE(*,*)' *** EOF reached in IN file before all data read ***'
      WRITE(*,*)' *** was attempting to find frequency',frq(ifreq)
      WRITE(*,*)
      ierr=3
      END

c********************************************************************
      SUBROUTINE WRITE_HP(ifreq,ix)

C     Writes the generated hydrophone fourier DATA to the OBS file
C     in the "hydrophone data" FORMAT. This consists of COMPLEX vectors
C     of the hydrophone output at each frequency.  Note this routine
C     is called once for each required frequency.    

      USE global
      INTEGER OBSFIL
      INTEGER index, idep, ifreq,ix

      PARAMETER ( OBSFIL = 30 )

      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      WRITE( OBSFIL, * ) TITLE(1:60)  ! Title
      WRITE( OBSFIL, * ) FRQ(ifreq), Nfrq, ix ! Current frequency,(optional: NFRQS)

      WRITE( OBSFIL, * ) Ndep
      WRITE( OBSFIL, 100 )( rdep(idep), idep=1,Ndep ) ! The actual depths.
100   FORMAT( 1F10.2 )
      
C     *** WRITE pressure for current frequency to the OBS file. ***
C     Note resp() CONTAINS the Nfrq fourier components at each hydrophone
C     depth.  Frequency varies slowest.  The first index represents range.
      DO idep = 1, Ndep
         index = (ifreq-1)*Ndep + idep
         WRITE( OBSFIL, * ) idep, resp(ix+(index-1)*nx)
     .        , ABS(resp(ix+(index-1)*nx))
      END DO
 
      END

c**************************************************************
      SUBROUTINE READ_HA(ierr)
      
C     Reads DATA from the IN file which CONTAINS either REAL or
C     synthetic hydrophone fourier vectors.
      
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER m_iphone
      PARAMETER(m_iphone=mx*100)
      INTEGER idep, idum, i, nran_ch, ifreq, ifflag, ierr
      REAL rdum, xfreq
      CHARACTER*80 dumch ! Used to pre-read each file record.
      INTEGER iphone_arr(m_iphone),iphone,iran_ch,index1,iran,index2
      INTEGER izero
      PARAMETER (izero=0)
C     READ the first line which should indicate the DATA FORMAT.

      READ(2,100) dumch
      WRITE(*,100)dumch
 100  FORMAT( a60 )             ! Format to read one record from the IN file.

C     Check to see that the FORMAT is "Hydrophone"
      IF ( INDEX(dumch,'Hydrophone') .EQ. 0 ) THEN
c         WRITE(*,*)' *** Warning *** The IN file format line does',
c     >        ' not contain "Hydrophone".'
c         WRITE(*,*)'Continuing anyway, with no error code introduced.'
         BACKSPACE(2)
      ENDIF

c
c
c
      DO 1100 idep=1,ndep      
         WRITE(*,*)' ******* Range number: ',idep,' *****'
         ifreq = 1
         ifflag = 1
         DO WHILE (.TRUE.)      ! Cycle through the frequencies.
C     There should normally be 3 important lines following the FORMAT line,
C     or preceeding the DATA for the next frequency.
C     They should be: 1) Title, 2) Frequency, 3) Number of depths.
C     Note we don't do anything with the title here.
      
            DO i=1,3
 2             READ(2,100,end=991)dumch ! Move current record into dumch.
               IF (iopt(6).eq.1 ) then
                  write(*,*)dumch
               endif 
               IF (dumch(1:1).ne.'!') then ! check to see if it is important.
                  IF (i.EQ.2) THEN ! dumch should contain the current frequency.
                     BACKSPACE(2)
                     READ(2,*,END=991) xfreq
C     Check to see that the frequency matches the specified one:
c       WRITE(*,*)'freqs', xfreq,frq(ifrq)
                     IF ( ABS( xfreq-frq(ifreq) ) .LT. 0.05 ) THEN
                        IFFLAG = 1
                     ELSE
                        ifflag = 0
                     ENDIF
                  ELSEIF (i.EQ.3) THEN ! dumch should contain # of depths.
                     BACKSPACE(2)
                     READ(2,*,END=991) nran_ch
                     IF (iopt(6).EQ.1 ) THEN
                        WRITE(*,*)' Number of ranges in IN-file:',
     .                       nran_ch
                     ENDIF 
                     IF (m_iphone.LT.nran_ch) STOP 
     1                    'increase dimension of mx'
                     IF ((nran_ch.NE.nx).AND.
     1                    (ifflag.EQ.1).AND.(ifreq.EQ.1)) THEN
                        WRITE(*,*)
     1                 '*** Warning ***.  Number of ranges mismatch.'
                        WRITE(*,*)
     1                 'Number of ranges reported in IN file:',nran_ch
                        WRITE(*,*)
     1                 'Number of ranges specified in DAT file:',nx 
c               ierr=1
c               WRITE(*,*)'Error code set. Returning to main program.' 
c               RETURN
                     ENDIF
                  ENDIF
               ELSE 
                  GOTO 2        ! The line read was just a comment so get next line.
               ENDIF
            ENDDO               ! For reading the 3 lines at the top.
c      WRITE(*,*) 'the first 3 lines read'
      
C     The next few lines should contain the actual channel depths.
C     Note there are 1 depths in each record.
            iphone=1
            DO i = 1, nran_ch
c     WRITE(*,*)'idep,mdep',i,mdep
 3             READ(2,100,END=991)dumch ! move current record into dumch.
               IF (dumch(1:1).NE.'!') THEN ! check to see if it was just a comment.
c     IF ((ifreq.EQ.1).AND.(ifflag.EQ.1)) THEN ! write out depths once.
c     WRITE(*,100) dumch ! write the records with the depths to screen.
c     ENDIF
c     BACKSPACE(2)
                  READ(dumch,*) rdum ! store the first depth in the rec.
c     READ(2,*) rdum      ! store the first depth in the rec.

c
c--- pg 000503 for now all are READ in !
c
c     IF ((iphone.LE.ndep).AND.
c     &           (ABS(rdum-rdep(iphone)).LT.0.005)) THEN
                  iphone_arr(i)=1
                  iphone=iphone+1
c            ELSE
c               iphone_arr(i)=izero
c            ENDIF
c--- pg 000503 for now all are READ in !


                  IF ( (i.EQ.1).AND.(ifflag.EQ.1).AND.
     1                 (ifreq.EQ.1).AND.(rdum.NE.xranges(1)) ) THEN
                     WRITE(*,*)
     1                    '*** WARNING*** 1st rec range in IN file:',
     2                    rdum
                     WRITE(*,*)
     1                    'is different than the range in DAT file:',
     2                    xranges(1)
                     WRITE(*,*)'No error code is set. Continuing.'
                  ENDIF
               ELSE 
                  GOTO 3        ! record was commented out so get next one.
               ENDIF
            ENDDO               ! i -loop
c      WRITE(*,*)'finished reading depths'
            IF ((iphone-1).NE.nx) THEN
               WRITE(*,*)' *****'
               WRITE(*,*) 'number of read phones from in file',
     1              iphone-1
               WRITE(*,*) 'number of requested  phones from dat file',
     1              nx
               STOP 'phone mismatch between *.dat and *.in'
            ELSE 
               IF (ifreq.EQ.1) WRITE(*,*)' Found', iphone-1,
     1              ' phones that are matching in depth'
            ENDIF

C     Now we should be positioned at top of hydrophone DATA. These DATA
C     records should consist of the channel number and the COMPLEX
C     pressure for the current frequency.

            IF ((nfrq.LT.10).OR. (iopt(6).EQ.1)) THEN
               IF (ifflag.EQ.1) THEN
                  WRITE(*,*)' Reading data for frequency:  ',xfreq
               ELSE
                  WRITE(*,*)' Skipping data for frequency: ',xfreq
               ENDIF
            ELSE
               IF (ifflag.EQ.1) THEN
                  WRITE(prtfil,*)' Reading data for frequency:  ',xfreq
               ELSE
                  WRITE(prtfil,*)' Skipping data for frequency: ',xfreq
               ENDIF
            ENDIF
c            WRITE(*,*)' Reading data for frequency:  ',xfreq
            iran=0
            DO iran_ch=1,nran_ch ! cycle through depths at current freq.
 4             READ(2,100,END=991)dumch
               IF (dumch(1:1).NE.'!') THEN
                  IF (ifflag.EQ.1) THEN ! Move stuff into data if valid freq.
                     IF (iphone_arr(iran_ch).EQ.1 ) THEN
                        iran=iran+1
c     BACKSPACE(2)
                        index2=(ifreq-1)*ndep+idep
                        index1=(index2-1)*nx
                        READ(dumch,*)idum, DATA( iran + index1 )
c     READ(2,*)idum, DATA( iran + index1 )
c     READ(2,*)idum, DATA( 1 +((ifreq-1)*Ndep+idep-1) )
c     WRITE(*,*)'reading data, idep,idep_ch',idep,idep_ch
                        IF (idum.NE.iran_ch) THEN
                           WRITE(*,*)' *** Warning ***',
     1                          ' Depth index mismatch at idep_ch = ',
     2                          iran_ch
                           WRITE(*,*)
     1                          'Error  set. Returning to main program'
                           ierr=2
                     RETURN
                  ENDIF
               ENDIF
            ENDIF
         ELSE 
            GOTO 4              ! Record was commented out so try next one.
         ENDIF
      ENDDO
c      WRITE(*,*)' the data was read'  
      IF (ifflag.EQ.1) THEN     ! Check to see if ifreq should be incremented.
         ifreq = ifreq + 1
         IF (ifreq.GT.nfrq) GOTO 900 ! All required freqs read, so finish.
      ENDIF
        
      ENDDO                     ! for main do-while in frequency

 900  CONTINUE
C     Should have finished reading DATA for the current frequency.
      WRITE(*,*)' All',ifreq-1,' frequencies was read from In-file'
 1100 CONTINUE                  ! range loop

      RETURN
      
 991  WRITE(*,*)' *** EOF reached in IN file before all data read ***'
      WRITE(*,*)' *** was attempting to find frequency',frq(ifreq)
      WRITE(*,*)
      ierr=3
      END
      
c********************************************************************
      SUBROUTINE WRITE_HA(ifreq,idep)

C     Writes the generated hydrophone fourier DATA to the OBS file
C     in the "HA-array data" FORMAT. This consists of COMPLEX vectors
C     of the hydrophone output at each frequency.  Note this routine
C     is called once for each required frequency.    

      USE global
      INTEGER OBSFIL
      INTEGER index, idep, ifreq,ix

      PARAMETER ( OBSFIL = 30 )

      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      WRITE( OBSFIL, * ) TITLE(1:60) ! Title
      WRITE( OBSFIL, * ) FRQ(ifreq), Nfrq, idep ! Current frequency,(optional: NFRQS)

      WRITE( OBSFIL, * ) Nx
      WRITE( OBSFIL, 100 )( xranges(ix), ix=1,Nx ) ! The actual ranges.
 100  FORMAT( 1F10.2 )
      
C     *** WRITE pressure for current frequency to the OBS file. ***
C     Note resp() CONTAINS the Nfrq fourier components at each hydrophone
C     depth.  Frequency varies slowest.  The first index represents range.
      DO ix = 1, Nx
         index = (ifreq-1)*Ndep + idep
         WRITE( OBSFIL, * ) ix, resp(ix+(index-1)*nx)
     .        , ABS(resp(ix+(index-1)*nx))
      END DO
      
      END

c**************************************************************
      SUBROUTINE READ_weight   
c(ierr)
      
C     Reads DATA from the WEI file which CONTAINS either REAL or
C     synthetic hydrophone fourier vectors.
      
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER m_iphone
      PARAMETER(m_iphone=mdep*100)
      INTEGER idep, idum, i, ndep_ch, ifreq, ifflag, ierr
      REAL rdum, xfreq
      CHARACTER*80 dumch ! Used to pre-read each file record.
      INTEGER iphone_arr(m_iphone),iphone,idep_ch,index1,iran,index2
      INTEGER izero,ifrq
      PARAMETER (izero=0)
C     READ the first line which should indicate the DATA FORMAT.
      READ(3,100) dumch
      WRITE(*,100)dumch
 100  FORMAT( a60 )             ! Format to read one record from the WEI file.
C     Check to see that the FORMAT is "Hydrophone"
      IF ( INDEX(dumch,'Hydrophone') .EQ. 0 ) THEN
c     pg         WRITE(*,*)' *** Warning *** The WEI file format line does',
c     pg     >        ' not contain "Hydrophone".'
c     pg         WRITE(*,*)'Continuing anyway, with no error code introduced.'
         BACKSPACE(3)
      ENDIF

c
c
c
      DO 1100 iran=1,nx      
         WRITE(*,*)' WEI: ******* Range number: ',iran,' *****'
         ifreq = 1
         ifflag = 1
         DO ifrq=1,nfrq         !WHILE (.true.)         ! Cycle through the frequencies.
C     There should normally be 3 important lines following the FORMAT line,
C     or preceeding the DATA for the next frequency.
C     They should be: 1) Title, 2) Frequency, 3) Number of depths.
C     Note we don't do anything with the title here.
      
            DO i=1,3
 2             READ(3,100,end=991)dumch ! Move current record into dumch.
               IF (dumch(1:1).ne.'!') then ! check to see if it is important.
                  IF (i.EQ.2) THEN ! dumch should contain the current frequency.
                     BACKSPACE(3)
                     READ(3,*,END=991) xfreq
C     Check to see that the frequency matches the specified one:
c     WRITE(*,*)'freqs', xfreq,frq(ifrq)
c     IF ( abs( xfreq-frq(ifreq) ) .lt. 0.05 ) then
                     IFFLAG = 1
c     ELSE
c     ifflag = 0
c     ENDIF
                  ELSEIF (i.EQ.3) THEN ! dumch should contain # of depths.
                     BACKSPACE(3)
                     READ(3,*,END=991) ndep_ch
                     IF (m_iphone.LT.ndep_ch) STOP 
     1                    'increase dimension of ndep'
                     IF ((ndep_ch.NE.ndep).AND.
     1                    (ifflag.EQ.1).AND.(ifreq.EQ.1)) THEN
                        WRITE(*,*)
     1                 '*** Warning ***.  Number of depths mismatch.'
                        WRITE(*,*)
     1                 'Number of depths reported in WEI file:',ndep_ch
                        WRITE(*,*)
     1                 'Number of depths specified in DAT file:',ndep 
c     ierr=1
c     WRITE(*,*)'Error code set. Returning to main program.' 
c     RETURN
                     ENDIF
                  ENDIF
               ELSE 
                  GOTO 2 ! The line read was just a comment so get next line.
               ENDIF
            ENDDO               ! For reading the 3 lines at the top.
c     WRITE(*,*) 'the first 3 lines read'
            
C     The next few lines should contain the actual channel depths.
C     Note there are 1 depths in each record.
            iphone=1
            DO i = 1, ndep_ch
c     WRITE(*,*)'idep,mdep',i,mdep
 3             READ(3,100,END=991)dumch ! move current record into dumch.
               IF (dumch(1:1).NE.'!') THEN ! check to see if it was just a comment.
c     IF ((ifreq.EQ.1).AND.(ifflag.EQ.1)) THEN ! write out depths once.
c     WRITE(*,100) dumch ! write the records with the depths to screen.
c     ENDIF
c     BACKSPACE(3)
c     READ(3,*) rdum      ! store the first depth in the rec.
                  READ(dumch,*) rdum ! store the first depth in the rec.
c     READ(3,*) rdum      ! store the first depth in the rec.
                  IF ((iphone.LE.ndep).AND.
     &                 (ABS(rdum-rdep(iphone)).LT.0.005)) THEN
                     iphone_arr(i)=1
c     WRITE(*,*)'iphone,rdum',iphone,rdum
                     iphone=iphone+1
                  ELSE
                     iphone_arr(i)=izero
                  ENDIF


                  IF ( (i.EQ.1).AND.(ifflag.EQ.1).AND.
     1                 (ifreq.EQ.1).AND.(rdum.NE.rdep(1)) ) THEN
                     WRITE(*,*)
     1                    '*** WARNING*** 1st rec depth in WEI file:',
     2                    rdum
                     WRITE(*,*)
     1                    'is different than the depth in DAT file:',
     2                    rdep(1)
                     WRITE(*,*)'No error code is set. Continuing.'
                  ENDIF
               ELSE 
                  GOTO 3        ! record was commented out so get next one.
               ENDIF
            ENDDO               ! i -loop
c     WRITE(*,*)'finished reading depths'
            IF ((iphone-1).NE.ndep) THEN
               
               WRITE(*,*)' *****'
               WRITE(*,*) 'number of read phones from *.WEI file',
     1              iphone-1
               WRITE(*,*) 'number of requested  phones from dat file',
     1              ndep
               STOP 'phone mismatch between *.dat and *.wei'
            ELSE 
               IF (ifreq.EQ.1) WRITE(*,*)' Found', iphone-1,
     1              ' phones that are matching in depth'
            ENDIF

C     Now we should be positioned at top of hydrophone DATA. These DATA
C     records should consist of the channel number and the COMPLEX
C     pressure for the current frequency.

            IF ((nfrq.LT.10).OR. (iopt(6).EQ.1)) THEN
               IF (ifflag.EQ.1) THEN
                  WRITE(*,*)' Reading data for frequency:  ',xfreq
               ELSE
                  WRITE(*,*)' Skipping data for frequency: ',xfreq
               ENDIF
            ELSE
               IF (ifflag.EQ.1) THEN
                  WRITE(prtfil,*)' Reading data for frequency:  ',xfreq
               ELSE
                  WRITE(prtfil,*)' Skipping data for frequency: ',xfreq
               ENDIF
            ENDIF
c     WRITE(*,*)' Reading data for frequency:  ',xfreq
            idep=0
            DO idep_ch=1,ndep_ch ! cycle through depths at current freq.
 4             READ(3,100,END=991)dumch
               IF (dumch(1:1).NE.'!') THEN
                  IF (ifflag.EQ.1) THEN ! Move stuff into data if valid freq.
                     IF (iphone_arr(idep_ch).EQ.1 ) THEN
                        idep=idep+1
c     BACKSPACE(3)
                        index2=(ifreq-1)*ndep+idep
                        index1=(index2-1)*nx
                        READ(dumch,*)idum, weight( iran + index1 )
c     READ(3,*)idum, DATA( iran + index1 )
c     READ(3,*)idum, DATA( 1 +((ifreq-1)*Ndep+idep-1) )
                        WRITE(*,*)'reading data, idep,idep_ch',
     >                       iran,idep_ch,ifreq,weight( iran + index1 )
                        IF (idum.NE.idep_ch) THEN
                           WRITE(*,*)' *** Warning ***',
     1                          ' Depth index mismatch at idep_ch = ',
     2                          idep_ch
                           WRITE(*,*)
     1                          'Error  set. Returning to main program'
                           ierr=2
                           RETURN
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE 
                  GOTO 4        ! Record was commented out so try next one.
               ENDIF
            ENDDO
            IF (ifflag.EQ.1) THEN ! Check to see if ifreq should be incremented.
               ifreq = ifreq + 1
               IF (ifreq.GT.nfrq) GOTO 900 ! All required freqs read, so finish.
            ENDIF
            
         ENDDO                  ! for main do-while in frequency
         
 900     CONTINUE
C     Should have finished reading DATA for the current frequency.
         WRITE(*,*)' All',ifreq-1,' frequencies was read from WEI-file'
 1100 CONTINUE                  ! range loop

      RETURN
      
 991  WRITE(*,*)' *** EOF reached in WEI file before all data read ***'
      WRITE(*,*)' *** was attempting to find frequency',frq(ifreq)
      WRITE(*,*)
      PAUSE
      ierr=3
      END

C .cfh.
c*************************************************************
      SUBROUTINE READ_AUN(ierraun)
C     Reads AUN from the unc file which CONTAINS either REAL or
C     synthetic covariance matrices.

      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      INTEGER idep,idum,idep1,ibart,i,ndep_ch,ifrq,irng,ierraun
      INTEGER iwarn,index2
      REAL rdum,xfreq
      REAL*8 xsum,xnorm,znoise
      CHARACTER*80 dumch
C      WRITE(*,*) 'entering read_aun'
      xcov_scale=0
      xcov_sum=0
      iwarn=0
      ibart=0
      DO ifrq=1,nfrq            ! Frequency loop
         DO irng=1,nx           ! Range loop
            ibart=ibart+1
 999        DO i=1,3            ! there should be 3 title lines. (Not including 'Format:...')
 2             READ(4,'(a80)',END=991) dumch
               IF (dumch(1:1).NE.'!') THEN
                  IF (i.EQ.2) THEN
                     BACKSPACE(4)
                     READ(4,*,END=991) xfreq
                  ELSEIF (i.EQ.3) THEN
                     BACKSPACE(4)
                     READ(4,*,END=991) ndep_ch
c                     write(*,*) 'ndep_ch', ndep_ch
                     IF (ndep_ch.NE.ndep) THEN
                        WRITE(*,*)' *****Error reading AUN mat******'
                        WRITE(*,*)' read ndep:', ndep_ch
                        WRITE(*,*)' req. ndep:', ndep
                        WRITE(*,*)' *****Error reading AUN mat******'
                        ierraun=1
                        RETURN
                     ENDIF
                  ENDIF
               ELSE 
                  GOTO 2
               ENDIF
            ENDDO
            
c     
c     Reading receiver depths
            DO idep=1,ndep
 3             READ(4,'(a80)',END=991) dumch
               IF (dumch(1:1) .NE. '!') THEN
                  BACKSPACE(4)
                  READ(4,*) rdum !rd
               ELSE 
                  GOTO 3
               ENDIF
               IF (iwarn.EQ.0) THEN
                  IF (idep.EQ.1) THEN
                     IF (rdum.NE.rdep(1)) THEN
                        WRITE(*,*)'*** WARNING*** 1st rec 
     .                       depth in covmat',rdum
                        WRITE(*,*)' is different from 
     .                       the input specified',rdep(1)
                        iwarn=1
                     ENDIF
                  ELSEIF (idep.EQ.ndep) THEN
c     WRITE(*,*)'ndep=',ndep
                     IF (rdum.NE.rdep(ndep)) THEN
                        WRITE(*,*)'*** WARNING*** last rec
     .                       depth in covmat',rdum
                        WRITE(*,*)' is different from the
     .                       input specified', rdep(ndep)
                        iwarn=1
                     ENDIF  
                  ENDIF  
               ENDIF
            ENDDO
            
c     >>>>> Here comes the squared root of data error covariance matrix
            DO idep1=1,ndep
               DO idep=1,ndep
 4                READ(4,'(a80)',END=992) dumch
                  IF (dumch(1:1).NE.'!') THEN
                     BACKSPACE(4)
                     index2=(idep+(ibart-1)*ndep-1)*ndep 
                     READ(4,*) idum,idum, Aun(idep1+index2) ! AUN
                  ELSE 
                     GOTO 4
                  ENDIF
               ENDDO
            ENDDO               ! loop over ndep
c-----
            IF (xfreq.LE.frq(ifrq)+0.05 .AND. 
     &           xfreq.GE.frq(ifrq)-0.05) THEN
               WRITE(*,*)' using Aun mat for freq',xfreq
            ELSE
               WRITE(*,*)' skipping Aun mat for freq:',xfreq
c     !               ' Req. freq: ', frq(ifrq)
               GOTO 999
            ENDIF
         ENDDO                  ! loop over bartlett    
      ENDDO                     ! loop over Frequency   
      WRITE(*,*)' readAun, Aun', Aun(1), Aun(nbart*ndep*ndep)
      RETURN
991   WRITE(*,*)' ***Problems reading Aun FILE *************'  
      ierraun=2
      RETURN
992   WRITE(*,*)' *** EOF encountered on Aun FILE ****'
      ierraun=3
      RETURN

      END

C .cfh.
C*******************************************************
      SUBROUTINE ROTMAT
C     Rotate the data covariance matrix: A*COV*A' 
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER ifrq,irng,index2
      INTEGER idep,jdep,kdep,ibart,idep1
      COMPLEX zx(ndep,ndep), zx1(ndep,ndep), zx2(ndep,ndep) 
      REAL*8  xsum,xnorm,znoise
      xcov_scale=0
      xcov_sum=0
      ibart=0
      DO ifrq=1,nfrq            ! Frequency loop
         DO irng=1,nx           ! Range loop
            ibart=ibart+1
C     put COV and A in the matrix form            
            DO jdep=1,ndep
               index2= (jdep+(ibart-1)*ndep-1)*ndep 
               DO idep=1,ndep
                  zx1(idep,jdep) = cov(idep+index2)
                  zx2(idep,jdep) = Aun(idep+index2)             
               ENDDO 
            ENDDO
            IF (iopt(6) .EQ. 1) THEN
               print *, "Matrix total cov follows"
               call PRINTMAT(zx1,ndep)
               print *, "Matrix Aun follows"
               call PRINTMAT(zx2,ndep)
            ENDIF
           
C     compute COV*A'
            Do idep = 1, ndep
               do jdep = 1, ndep
                  zx(idep,jdep) = (0.0e0,0.0e0)
                  do kdep = 1, ndep
                     zx(idep,jdep) = zx(idep,jdep) + 
     &                    zx1(idep,kdep) * conjg(zx2(jdep,kdep))
                  enddo
               enddo
            enddo
C     compute A*(COV*A')
            Do idep = 1, ndep
               do jdep = 1, ndep
                  zx1(idep,jdep) = (0.0e0,0.0e0)
                  do kdep = 1, ndep
                     zx1(idep,jdep) = zx1(idep,jdep) + 
     &                    zx2(idep,kdep) * zx(kdep,jdep)
                  enddo
               enddo
            enddo
C---- put A*COV*A' in the vector form            
            DO jdep=1,ndep
               index2= (jdep+(ibart-1)*ndep-1)*ndep 
               DO idep=1,ndep
                  cov(idep+index2) = zx1(idep,jdep)
               ENDDO 
            ENDDO 
C---- find the trace of covariance matrix
            xsum=0
            DO idep1=1,ndep
               index2=(idep1+(ibart-1)*ndep-1)*ndep 
               idep=idep1
               xsum=xsum+cov(idep+index2) 
            ENDDO
            xcov_trace(ibart)=xsum
c---- Find largest eigenvalue
            IF (ndep .GT. 1000) THEN
               WRITE(*,*)' The new covariance matrix is large,',
     &              ' thus the computation of the SNR has been skipped'
            ELSE
               CALL svdmax(ibart,xsum,znoise)
               zcov_noise(ibart)=znoise
               xcov_lamb(ibart)=xsum
c          WRITE(*,*)' trace and the lagest eigenvalue of cov-matrix:'
c          WRITE(*,*) xcov_trace(ibart),xcov_lamb(ibart)
               WRITE(*,'(a,f10.2)')' estimated SNR from cov matrix:', 
     &              10*LOG10((xcov_lamb(ibart)-znoise)/(znoise))
               IF (xcov_trace(ibart).LT.0) THEN 
                  WRITE(*,*)' *** Error in covariance matrix: ',
     &                 'The trace of the covariance matrix 
     &                 is negative! '
                  STOP
               ENDIF
            ENDIF
c
c-----  normalize covarince matrix with diagonal or max eigenvalue
c
            IF (iopt(20) .GE. 1) THEN
               IF (iopt(20) .EQ. 1) THEN 
                  xnorm = xcov_trace(ibart)
                  xcov_lamb(ibart) = xcov_lamb(ibart)/xnorm
                  xcov_trace(ibart) = 1
               ELSEIF (iopt(20) .EQ. 2) THEN 
                  xnorm=xcov_lamb(ibart)
                  xcov_trace(ibart)=xcov_trace(ibart)/xnorm
                  xcov_lamb(ibart)=1
               ELSE
                  STOP ' iopt(20)>= 2 not defined'
               ENDIF
               WRITE(*,*)' The covariance matrix is renormalized'
               DO idep1=1,ndep
                  index2=(idep1+(ibart-1)*ndep-1)*ndep 
                  DO idep=1,ndep
                     cov(idep+index2)=cov(idep+index2)/xnorm
                  ENDDO
               ENDDO
               WRITE(*,*)'...renormalized with xnorm ', xnorm
            ENDIF
c---  Find a scaling factor to be used when multipling bartlett power.
c     WRITE(*,*) xcov_trace(ibart),xcov_lamb(ibart)
            xcov_scale = xcov_scale+
     &           10.*LOG10(xcov_lamb(ibart))
c     &       10.*LOG10(xcov_trace(ibart)-xcov_lamb(ibart))
c     WRITE(*,*)'xcov_scale',xcov_scale
            xcov_sum=xcov_sum+xcov_trace(ibart)
         ENDDO
      ENDDO
      IF (iopt(6) .EQ. 1) THEN
         print *, "column new cov follows"
         write(*,*) 'new cov', cov
      ENDIF
      END

C .cfh. 
C*******************************************************
      SUBROUTINE ROTVEC
C     Rotate the replica vector: A*resp
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER ifrq,ix, index, index2, index3,j1start
      INTEGER idep,jdep,kdep,ibart,i,j
      COMPLEX cx,cxt(mdep)

c      write(*,*) 'nfrq,nx,ndep',nfrq,nx,ndep
      IF (iopt(6) .EQ. 1) THEN
         print *, "column Aun follows"
         write(*,*) Aun
      ENDIF

      DO ifrq=1,nfrq            ! Frequency loop
         DO ix=1,nx             ! Range loop
            ibart=ix+((ifrq-1))*nx ! which observation 
            index3=(ibart-1)*ndep-1 
c     now for each covariance matrix 
            j1start=((ifrq-1))*ndep ! where the response starts
            DO i=1,ndep
               cx=(0.d0,0.d0)
               DO j=1,ndep
C     write(*,*) 'ix+(j+j1start-1)*nx',ix,j,j1start,nx
                  index=ix+(j+j1start-1)*nx
                  index2=i+(j+index3)*ndep 
C                  write(*,*) 'index2, index',index2, index
                  cx=cx+Aun(index2)*resp(index)     
c     WRITE(*,*)'cx',cx,cov(index2),resp(index),index2,index
               ENDDO
               index=ix+(i+j1start-1)*nx
               IF (iopt(6) .EQ. 1) THEN 
                  write(*,*) 'index,resp,cx', index,resp(index),cx
               ENDIF
               cxt(index) = cx
            ENDDO 
         ENDDO
      ENDDO

      ix = nfrq*nx*ndep
      DO i = 1,ix
         resp(i) = cxt(i)
      ENDDO
      END

