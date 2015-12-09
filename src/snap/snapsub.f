

      SUBROUTINE READ_SP(ierr)
      
C     Reads the source spectrum 
c     from the SOU file which contains either real or
C     synthetic hydrophone fourier vectors.
      
      INCLUDE '../comopt.h'
      INCLUDE '../comforw.h'
      include '../comsnap.h'
      integer m_iphone
      parameter(m_iphone=mdep*100)
      INTEGER idep, idum, i, ndep_ch, ifreq, ifflag, ierr
      REAL rdum, xfreq
      CHARACTER*80 dumch ! Used to pre-read each file record.
      integer iphone_arr(m_iphone),iphone,idep_ch,index1,index2
      integer izero
      parameter (izero=0)
C     Read the first line which should indicate the data format.

      READ(42,100) dumch
      WRITE(*,100)dumch
 100  FORMAT( a60 )             ! Format to read one record from the IN file.

C     Check to see that the format is "Hydrophone"
      IF ( index(dumch,'Hydrophone') .eq. 0 ) then
         WRITE(*,*)' *** Warning *** The IN file format line does',
     >        ' not contain "Hydrophone".'
         WRITE(*,*)'Continuing anyway, with no error code introduced.'
         backspace(42)
      ENDIF

c
c
      ifreq = 1
      ifflag = 1
      DO WHILE (.true.)         ! Cycle through the frequencies.
C     There should normally be 3 important lines following the format line,
C     or preceeding the data for the next frequency.
C     They should be: 1) Title, 2) Frequency, 3) Number of depths.
C     Note we don't do anything with the title here.
      
      DO i=1,3
 2       READ(42,100,end=991)dumch ! Move current record into dumch.
         IF (dumch(1:1).ne.'!') then ! check to see if it is important.
            IF (i.eq.2) then    ! dumch should contain the current frequency.
               backspace(42)
               READ(42,*,end=991) xfreq
C     Check to see that the frequency matches the specified one:
c       write(*,*)'freqs', xfreq,frq(ifrq)
               IF ( abs( xfreq-frq(ifreq) ) .lt. 0.05 ) then
                  IFFLAG = 1
               ELSE
                  ifflag = 0
               ENDIF
            ELSEIF (i.eq.3) then ! dumch should contain # of depths.
               backspace(42)
               READ(42,*,end=991) ndep_ch
               if (msrd.lt.ndep_ch) stop 
     1                              'increase dimension of mdep'
               IF ((ndep_ch.ne.nsrd).and.
     1            (ifflag.eq.1).and.(ifreq.eq.1)) then
               WRITE(*,*)'*** Warning ***.  Number of depths mismatch.'
               WRITE(*,*)'Number of depths reported in IN file:',ndep_ch
               WRITE(*,*)'Number of depths specified in DAT file:',nsrd 
c               ierr=1
c               WRITE(*,*)'Error code set. Returning to main program.' 
c               RETURN
               ENDIF
            ENDIF
         ELSE 
            GOTO 2  ! The line read was just a comment so get next line.
         ENDIF
      ENDDO         ! For reading the 3 lines at the top.
c      write(*,*) 'the first 3 lines read'
      
C     The next few lines should contain the actual channel depths.
C     Note there are 1 depths in each record.
      iphone=1
      DO i = 1, ndep_ch
c          write(*,*)'idep,mdep',i,mdep
 3       READ(42,100,end=991)dumch ! move current record into dumch.
         IF (dumch(1:1).ne.'!') then ! check  if it was just a comment.
c            IF ((ifreq.eq.1).and.(ifflag.eq.1)) then ! write out dep once.
c             WRITE(*,100) dumch ! write the records with the dep to screen.
c            ENDIF
c            backspace(42)
c            READ(42,*) rdum      ! store the first depth in the rec.
            READ(dumch,*) rdum      ! store the first depth in the rec.
c            READ(42,*) rdum      ! store the first depth in the rec.
            if ((iphone.le.nsrd).and.
     &           (abs(rdum-sdep(iphone)).lt.0.005)) then
               iphone_arr(i)=1
c     write(*,*)'iphone,rdum',iphone,rdum
               iphone=iphone+1
            else
               iphone_arr(i)=izero
            endif


            IF ( (i.eq.1).and.(ifflag.eq.1).and.
     1           (ifreq.eq.1).and.(rdum.ne.sdep(1)) ) then
               WRITE(*,*)'*** WARNING***1st rec depth in SOU file:',rdum
               WRITE(*,*)'is different than the depth in DAT file:',
     1              sdep(1)
               WRITE(*,*)'No error code is set. Continuing.'
            ENDif
         ELSE 
            GOTO 3              ! record was commented out so get next one.
         ENDIF
      ENDDO                     ! i -loop
c      write(*,*)'finished reading depths'
      if ((iphone-1).ne.nsrd) then
          stop 'phone mismatch between *.dat and *.in'
      else 
          if (ifreq.eq.1) write(*,*)' Found', iphone-1,
     1              ' phones that are matching in depth'
      endif

C     Now we should be positioned at top of hydrophone data. These data
C     records should consist of the channel number and the complex
C     pressure for the current frequency.

      if ((nfrq.lt.10).or. (iopt(6).eq.1)) then
         IF (ifflag.eq.1) then
            WRITE(*,*)' Reading data for frequency:  ',xfreq
         ELSE
            WRITE(*,*)' Skipping data for frequency: ',xfreq
         ENDIF
      else
         IF (ifflag.eq.1) then
            WRITE(prtfil,*)' Reading data for frequency:  ',xfreq
         ELSE
            WRITE(prtfil,*)' Skipping data for frequency: ',xfreq
         ENDIF
      endif
c            WRITE(*,*)' Reading data for frequency:  ',xfreq
      idep=0
      DO idep_ch=1,ndep_ch      ! cycle through depths at current freq.
 4       READ(42,100,end=991)dumch
         IF (dumch(1:1).ne.'!') then
            IF (ifflag.eq.1) then ! Move stuff into data if valid freq.
               if (iphone_arr(idep_ch).eq.1 ) then
                  idep=idep+1
c     backspace(42)
c                  index2=(ifreq-1)*ndep+idep
c                  index1=(index2-1)*nx
                  READ(dumch,*)idum, sourcesp( idep,ifreq)
c     READ(42,*)idum, data( iran + index1 )
c     READ(42,*)idum, data( 1 +((ifreq-1)*Ndep+idep-1) )
c     write(*,*)'reading data, idep,idep_ch',idep,idep_ch
                  IF (idum.ne.idep_ch) then
                     WRITE(*,*)' *** Warning ***',
     >                    ' Depth index mismatch at idep_ch = ',idep_ch
                     WRITE(*,*)'Error  set. Returning to main program'
                     ierr=2
                     RETURN
                  ENDIF
               endif
            ENDIF
         ELSE 
            GOTO 4              ! Record was commented out so try next one.
         ENDIF
      ENDDO
c      write(*,*)' the data was read'  
        IF (ifflag.eq.1) then ! Check to see if ifreq should be incremented.
           ifreq = ifreq + 1
           IF (ifreq.gt.nfrq) goto 900 ! All required freqs read, so finish.
        ENDIF
        
      ENDDO                     ! for main do-while in frequency

 900  CONTINUE
C     Should have finished reading data for the current frequency.
      write(*,*)' All',ifreq-1,' frequencies was read from In-file'
 1100 continue ! range loop

      RETURN
      
991   write(*,*)' *** EOF reached in SOU file before all data read ***'
      write(*,*)' *** was attempting to find frequency',frq(ifreq)
      write(*,*)
      ierr=3
      END


