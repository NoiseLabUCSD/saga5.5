      subroutine forwardmodel(iopt,mopt)
      integer  mopt,i,iopt(mopt)
      DO i=1,mopt
         iopt(i)=0
      ENDDO
      iopt(1)=2
      iopt(12)=1                ! using three indexes for adressing variable
      iopt(30)=11               ! 11 is ORCA
      end 
c==================================================================
      Subroutine input
c==================================================================

      USE global
      use fgs_com      
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'
c  Local use variables
      character*120 line
      character*1 lineopt(120)
      equivalence (line,lineopt)
      real delf,f_min,f_max
c     pg      real sdep 
      real rdup,rdlo,rdstep   
      real rrlo,rrhi,rrstep   
      integer j1,j2,nline,nsou,nelem,idum(3),kd
      
      integer ndhs,nstep,nave,idhs,t0,tfact,ntemp,npert   
      real rdum(3)
      integer itemp,icount,ierr,jj,ndum,inu
      inu=0                     !cfh
      ierr=0
c
c---	Title, options
c
c---  read the GA-parameters
      call readinputstart

c       read(1,'(A)') line
c  Be careful not to read in user comments in this line,
c  look for the index of '!' else use 50

c      j1=index(line,'!')
c      if(j1==0) j1=50

c      if( index(line(1:j1),'V')  ) arraytype=v_arr
c      if( index(line(1:j1),'H')  ) arraytype=h_arr

c ASSA or FGS parameters
c      nline=3
c      read(1,*,err=510,end=510) ntemp,npert,idum(1)

c      nline=4
c      read(1,*,err=510,end=510) t0,tfact,snr

c  Forward model options

       read(1,'(A)') line
       write(*,*)line
       j1=index(line,'!')
       if(j1==0) j1=50
       write(*,*)'index ..',  index(line(1:j1),'V')
c       stop
       do j=1,j1  
          if( lineopt(j).eq.'V' ) arraytype=v_arr
          if( lineopt(j).eq.'H' ) arraytype=h_arr
          if( lineopt(j).eq.'t' ) tilt=.TRUE.
          if( lineopt(j).eq.'C' ) optc=.TRUE.
          if( lineopt(j).eq.'p' ) optp=.TRUE.
       enddo
       if( tilt ) write(*,*) ' Tilt (straight-line) of Vertical Array '
       if( optc ) write(*,*) ' Using ORCA complex root finder '
       if( optp ) write(*,*) ' Allow only positive velocity gradients '
       if( arraytype.eq.v_arr ) then 
          write(*,*) ' Vertical array'
       else if ( arraytype.eq.h_arr ) then 
          write(*,*) ' Horizotal array'
       else
          stop 'array type not defined'
       endif
c------------------------------------------------------------------------
c  Frequencies:
c------------------------------------------------------------------------

       write(*,*)'reading nfrq..'
       nline=6
      read(1,*,err=510,end=510) nfrq

c      allocate( frq(iabs(nfrq)) )

      nline=7
      if (nfrq.gt.0) then
         READ(1,*,err=510,end=510)(frq(i),i=1,nfrq)
         write(*,'(A,10f10.3)')' Frequencies:',(frq(i),i=1,nfrq)
      else
         nfrq=-nfrq
         read(1,'(A)') line
         READ(line,*,err=510,end=510) delf,f_min,f_max
         do j=1,nfrq
            frq(j)=f_min+float(j-1)*delf
         enddo
      endif
      fminorca=frq(1)
      fmaxorca=frq(nfrq)

c------------------------------------------------------------------------
c  Water depth and SSP (max 500 pts)
c------------------------------------------------------------------------

      nline=8
      read(1,*,err=510,end=510) nsvp
      nline=nline+1      
      write(*,*)'Number of sound speed points:', nsvp		
      do j=1,nsvp
         read(1,*,end=510,err=510)zsvp(j),csvp(j)
         nline=nline+1
      enddo
      
c------------------------------------------------------------------------
c  Bottom layers (ORCA style, max 100 pts)
c------------------------------------------------------------------------

      read(1,*,err=510,end=510) nlayb
      nline=nline+1        
      write(*,*) 'Number of bottom layers:', nlayb                
      if (nlayb .gt. 100 .or. nlayb .lt. 0) then
         write(*,*) 'Error in number of sediment layers, nlayb =',nlayb
         stop
      endif
      do j=1,nlayb
         read(1,'(A)') line
         read(line,*,err=510,end=510) 
     &      hb(j),((geob(j1,j2,j),j1=1,2),j2=1,5)
         nline=nline+1
      enddo
                        
c------------------------------------------------------------------------
c  Subbottom (ORCA style)
c------------------------------------------------------------------------

      read(1,'(A)') line
      read(line,*,end=510,err=510) (geob(1,j,nlayb+1),j=1,5)
      nline=nline+1
      read(1,*)      ! Skip blank line
       
c------------------------------------------------------------------------ 
c  Read source depth  
c------------------------------------------------------------------------

      nsou=1
      read(1,*,err=510,end=510) sdep(1)
      nline=nline+1
      write(*,*)'Source depth (m) :',sdep(1)
     
c------------------------------------------------------------------------ 
c  Read receiver depths  
c------------------------------------------------------------------------

      read(1,'(A)',err=510,end=510) line
      nline=nline+1

      if (tilt) then
         read(line,*,err=510,end=510) rdup,rdlo,ndep,dtilt
      else
         read(line,*,err=510,end=510) rdup,rdlo,ndep
      endif
      
      if( tilt .and. ndep < 0 ) then
         write(*,*) ' Tilt requires equi-distant spacing '
         stop
      endif

      allocate( rdarr(iabs(ndep)+1) )         ! ???
	
      if( ndep .lt. 0 ) then
         read(1,'(A)',err=510,end=510) line
         read(line,*,err=510,end=510) (rdarr(j),j=1,abs(ndep))
         nline=nline+1
      else if( ndep .eq. 1 ) then
         rdarr(1)=rdup
      else                                    ! ( ndep > 1 )		
         rdstep=(rdlo-rdup)/float(ndep-1)
         do  j=1,ndep
            rdarr(j)=float(j-1)*rdstep+rdup
         enddo
       endif
       	
      ndep = iabs(ndep)

       do  j=1,ndep
          rdep(j)=rdarr(j)
       enddo

       write(*,*) 'Receiver depths (m):',(rdarr(j),j=1,iabs(ndep))
        
c------------------------------------------------------------------------ 
c  Read receiver ranges.  
c------------------------------------------------------------------------
c---  read the range-block
      call  read_range_inp(ierr,ndum)

c      read(1,*,err=510,end=510) rrlo,rrhi,nran
c      nline=nline+1
      nran=nx

      allocate( rrarr(nran+1) )         ! ???
      do i=1,nx
         rrarr(i)= xranges(i)
      enddo

	
      write(*,*) 'Receiver ranges (km) :',(rrarr(j),j=1,nran)

      read(1,*)   ! Skip one line
 

c------------------------------------------------------------------------
c  Map to ORCA parameters:
c------------------------------------------------------------------------

      nfcw = abs(nfrq)                   ! Frequencies
      fcw(1:nfrq) = frq(1:abs(nfrq))

      nzs = abs(nsou)                    ! Source depth
      zsrc(1) = sdep(1)

      nrec = abs(ndep)                   ! Receiver depths
      zrec(1:nrec) = rdarr(1:abs(ndep))

      nsrc = abs(nran)                   ! Receiver ranges
      rkm(1:nsrc) = rrarr(1:abs(nran))

      select case( arraytype )

      case( v_arr )                        
         nHyd = ndep                             
         print *, ' Vertical Array ',nHyd,' hydrophones '
         valength=rdarr(ndep)-rdarr(1)               
         write(*,*) ' Array length [m]: ',valength         	
         if( tilt ) then   ! Def range pt for every h/ph
            nsrc=ndep
            rkm(1:nsrc)=rrarr(1)
         endif
		
      case( H_arr )      
         nHyd = nran   
         print *, ' Horizontal Array ',nHyd,' hydrophones '

      case default  !  Unspecified array
         print *, ' Error: Array type not specified '
         stop	
             
      end select

c--- should the  EOF be read ? 
      if (iopt(17).eq.1) then
        call eofinit
      endif


c*** create text strings for physical parameters
c                   123456789012345678901234567890
      phystxt(1) = 'Water depth (m) $'
      phystxt(2) = 'Sediment sound speed (m/s)$'
      phystxt(3) = 'Sediment shear sound speed (m/s)$'
      phystxt(4) = 'Attenuation (dB/lambda)$'
      phystxt(5) = 'S-attenuation (dB/lambda)$'
      phystxt(6) = 'Sediment density (g/cm3)$'
      phystxt(7) = 'Layer thickness (m)$'
      phystxt(8) = 'Source depth (m) $'
      phystxt(9) = 'Source range (km) $'
      phystxt(11)= 'Shape coefficient $'
      phystxt(14)= 'Bottom density (g/cm3)$'
      phystxt(15)= 'Receiver depth (m)$'
      phystxt(17)= 'Sediment depth (m)$'
      phystxt(19)= 'Vertical array tilt$'
      phystxt(20)= 'Sound speed in water $'
      phystxt(29)= 'Error variance$'

      do 8 j=1,mphys
         phystxt2(j)='                                               '
 6       DO 7 I=40,1,-1
            IF (phystxt(j)(I:I).NE.'$') GO TO 7
            phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c     write(*,*)phystxt2(j)
            GO TO 8
 7       CONTINUE
    8 CONTINUE
      phystxt2(9)='Source range (m)'
          
c---- read the optimization parameter
      call readoptparm
      DO i=1,nparm
         IF (fmin(i).gt.fmax(i)) THEN
            WRITE(*,*)' *** fmin > fmax for parm',i
            ierr=1
         ENDIF
c     .cfh.   
         IF ( par2phy(i).lt.1  .or. par2phy(i).gt.29 .or.
     &        par2phy(i).eq.21 .or. par2phy(i).eq.22 .or.
     &        par2phy(i).eq.23 .or. par2phy(i).eq.24 .or.
     &        par2phy(i).eq.25 .or. par2phy(i).eq.26 .or.
     &        par2phy(i).eq.27 .or. par2phy(i).eq.28 .or.
     &        par2phy(i).eq.10 .or. par2phy(i).eq.12 .or.
     &        par2phy(i).eq.13 .or. par2phy(i).eq.14 .or.
     &        par2phy(i).eq.16 .or. par2phy(i).eq.17 .or.
     &        par2phy(i).eq.18 ) THEN
            WRITE(*,*)' *** par2phy not correct, parm',i
            ierr=1
         ENDIF
         IF (par2phy(i).eq.1 ) THEN
            if (fmin(i).lt.zsvp(nsvp-1)) then
               write(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** Durring optimization'
               WRITE(*,*)' the water-depth can get above the',
     &              ' second-last point in the water'
               ierr=1
            endif
         ENDIF
         IF (par2phy(i).eq.11 ) THEN
            if (par2lay(i).gt.neof) then
               write(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** The shapecoffient number is not defined'
               WRITE(*,*)' par2lay(i)', par2lay(i)
               ierr=1
            endif
         ENDIF
         IF (par2phy(i).eq.11 ) THEN
            if (par2lay(i).gt.neof) then
               write(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** the second pointer must be less than',
     &              ' the number of shapefunctions, Neof=',neof
               WRITE(*,*)' par2lay(i)', par2lay(i)
               ierr=1
            endif
         ENDIF
         IF (par2phy(i) .eq. 29 ) inu = 1
        
         IF ((par2phy(i).eq.1) .or. (par2phy(i).eq.8)  .or. 
     1        (par2phy(i).eq.9).or. (par2phy(i).eq.13) .or. 
     1        (par2phy(i).eq.14).or. (par2phy(i).eq.15) .or.
     1        (par2phy(i).eq.17).or. (par2phy(i).eq.19)) THEN
            itemp=par2phy(i) 
            icount=0 
            do jj=1,nparm
               if (itemp.eq.par2phy(jj)) icount=icount+1 
            enddo
            if (icount.ne.1) then
               write(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** The parameter is defined twice as',
     1              ' optimization parameter'
               ierr=1
            endif
         ENDIF
      enddo
 50   continue
c     .cfh.

      IF (inu .eq. 1 .and. isubopt(36) .eq. 0) THEN
         WRITE(*,*)' *** Optimzation error variance #',i
         WRITE(*,*)' *** Need to specify the saga option *' 
         ierr = 1
      ELSEIF (inu .eq. 0 .and. isubopt(36) .eq. 1) THEN
         isubopt(36) = 2
         IF (iopt(4) .EQ. 4 .and. isubopt(4) .EQ. 2) ierr = 1
      ENDIF
c***  errors ?
      IF (ierr .eq. 1) stop 'stopped in input'
      return
510   print *,'End or error reading input file at line ',nline
      stop
      END
      
c======================================================================
      Subroutine forwinit 
c======================================================================
      USE global
      use fgs_com
      use i_o_com
      implicit none
      INCLUDE 'comopt.h'
      include 'comforw.h'
      integer kd
      real E_out,rref,rad
      complex p_out(nfrq,nhyd)	
      integer id,ifreq,index
      integer   flagpu,J1,J2
      common /flagpu/flagpu
c------------------------------------------------------------------------
c     Adjust parameters of i_o_com of ORCA v2.01 
c     Only required at first pass.
c------------------------------------------------------------------------
      ocount=0
      if ( ocount .eq. 0 ) then       
         iiAih(1)=-1.0          ! homogeneous halfspace
         iiAih(2)=-1.0          ! homogeneous halfspace
         kkblm=1                ! search for branch line modes
         iirx=1                 ! real axis search                          
         if(optc) then          ! complex-plane search
            iirx=0
            iiAih(1)=90.005
            iiAih(2)=-1.0
            kkblm=0
         endif          
      endif

c------------------------------------------------------------------------
c     Map inversion parameters to ORCA. Avoid for first call due to
c     precision problem (real to double conversion),
c-----------------------------------------------------------------------
      flagpu=0
      if (iopt(6).eq.1) flagpu=-30000
      entry forw2

      lwascomp=1
      flagpu=1 +flagpu
      
c     
c---- Use of EOF
c     
      if (iopt(17).eq.1) then
         call eofval
         if (lwascomp.eq.-1) then
            return
         endif
      endif
c------------------------------------------------------------------------
c     Vertical Array. Elements numbered from #1 at top. 
c     Assume equidistant spaced array. Ignore adjustment of depths.
c------------------------------------------------------------------------

      rad=3.14159265/180.
      if( V_arr==arraytype  ) then
        rref=rkm(1)
         if( tilt ) then
            rkm(1)=rref
            do j=2,nrec
               rkm(j)=rref+
     &              1.e-3*sin(dtilt*rad)*valength*(j-1)/(nrec-1)
            enddo
         endif
      endif
      
c------------------------------------------------------------------------
c     Horizontal Array. Receiver Ranges [km] must translate the whole array
c     Assuming rkm(1) is closest element of array, holds the perturbed 
c     starting point of the array
c------------------------------------------------------------------------

      if( H_arr==arraytype  ) then
         rref=rkm(1)
         do j=2,nran
            rkm(j)=rref+(rrarr(j)-rrarr(1))
         enddo
      endif

c------------------------------------------------------------------------
c     Must adjust rmin,rmax if array has been moved.
C     Call ORCA
c------------------------------------------------------------------------

      rmin=min(rkm(1),rkm(nsrc))
      rmax=max(rkm(1),rkm(nsrc))

      ocount=ocount+1
      if (flagpu.lt.3) then
         write(*,*) 'Calling orca ...'
         write(*,*) 'frequencies',(frq(j),j=1,nfrq)
         write(*,*)'Number of sound speed points:', nsvp
         do j=1,nsvp
            write(*,*)zsvp(j),csvp(j)
         enddo
         write(*,*) 'Number of bottom layers:', nlayb                
         do j=1,nlayb
            write(*,*)hb(j),((geob(j1,j2,j),j1=1,2),j2=1,5)
         enddo
         
         write(*,*)'Bottom', (geob(1,j,nlayb+1),j=1,5)
         write(*,*)'Source depth (m) :',sdep(1)
         write(*,*) 'Receiver depths (m):',(rdarr(j),j=1,iabs(ndep))
         write(*,*) 'Receiver ranges (km) :',(rrarr(j),j=1,iabs(nran))
      endif


      call orca( ocount )
      
c------------------------------------------------------------------------
c     Extract complex pressures across array.
c     Stored in complex*8 tf(1:nfrq,1:nrec,1:nsrc).
c     where indexes are: (frequency:rec_depth:rec_range).
c------------------------------------------------------------------------
      select case( arraytype )
      
      case( V_arr )
         i=1                    ! just 1 range
         if( tilt ) then
            do ifreq=1,nfrq
               do id=1,ndep	! irin
                  index=(id +((ifreq-1))*ndep-1)*nran
                  resp(i+index)=-conjg(tf(ifreq,id,id))
               enddo
            enddo
         else
	    do ifreq=1,nfrq
	       do id=1,ndep	! irin
		  index=(id +((ifreq-1))*ndep-1)*nran
		  resp(i+index)=-conjg(tf(ifreq,id,1))
c     write(*,*)'resp',resp(i+index)
	       enddo
	    enddo
	 endif

         
      case( H_arr )
         do ifreq=1,nfrq
            id=1                ! irin
            do i=1,nran		! just 1 range
               index=(id +((ifreq-1))*ndep-1)*nran
               resp(i+index)=-conjg(tf(ifreq,1,i))
            enddo
         enddo
c     p_out = tf(:,1,:)
         
      end select
 
      ocount=ocount+1
      lwascomp=1
      return
      end

c*********************************************
      SUBROUTINE writeTL()
      USE global
      use fgs_com
      use i_o_com
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
 
      INTEGER ifreq,ii,id,index,mx1
c     pg      INTEGER*4 jj0,jsrc,jrec,j,mx1
      REAL*4 rspace,freqs,dt1
      LOGICAL bintrf
      INTEGER luttrf,np,nf1
      data luttrf/16/
c      CHARACTER*30 outroot      !not used
      bintrf=.FALSE.
      np=nran

c     bintrf=.TRUE.


c     
c     rspace = range increment for receiver position
c     pln      rkm(1)=rkm(1)*1.0e-3
      IF (nran.EQ.1) THEN
         rspace=1
      ELSE
         rspace=rkm(2)-rkm(1)
      END IF
      mx1=nf1+nfrq-1
      if (fsbb==0) then
         dt1=1
      else
         dt1=1.E0/fsbb
      endif
      freqs=(fcw(1)-fcw(nfcw))/2.+fcw(1)
      if (nfrq>1) then
         nf1=MAX(1,NINT(fcw(1)/(fcw(2)-fcw(1)))) + 1
      else
         nf1=1
      endif
c     sd=zsr(mzsrc(1))
      if (iforwt==1) then
      write(*,*)'Calling trfhead'
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     rkm(1),rspace,nfftbb,nf1,mx1,dt1,freqs,sd,
     &     bintrf,ndep,1,nran)
      endif
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
c         DO ifreq=1,nfrq
c            DO i=1,np           ! ranges
c               DO id=1,ndep     ! irin
c                  index=(id +((ifreq-1))*ndep-1)*np
c                  write(*,*)'writing trf..'
                  WRITE(luttrf,*)
                  WRITE(luttrf,'(1000000f8.3)')
     &    ((((20*log10(abs(resp(ii+(id +((ifreq-1))*ndep-1)*np))))
     &        ,id=1,ndep),ii=1,nran), ifreq=1,nfrq)
c     pg 1 may 2000 sign for prosim      -imag(resp(i+index))
c               ENDDO            ! depth
c            ENDDO               ! ranges
c         ENDDO                  ! freq
c     pg         END IF
c     
c         CLOSE(luttrf)
c      ENDIF
c     
      RETURN
      END  !wrtetrf

ccc
 
      SUBROUTINE writecomplex
      USE global
      use fgs_com
      use i_o_com
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'

      INTEGER ifreq,ii,id,index,mx1
c     pg      INTEGER*4 jj0,jsrc,jrec,j,mx1
      REAL*4 rspace,freqs,dt1
      LOGICAL bintrf
      INTEGER luttrf,np,nf1
      data luttrf/16/
c     CHARACTER*30 outroot      !not used
      bintrf=.FALSE.

c     bintrf=.TRUE.
      np=nran
c     
c     rspace = range increment for receiver position
c     pln      rkm(1)=rkm(1)*1.0e-3
      IF (nran.EQ.1) THEN
         rspace=1
      ELSE
         rspace=rkm(2)-rkm(1)
      END IF
      dt1=1.E0 !/fsbb
c      freqs=(fmaxdum-fmindum)/2.+fmindum
      freqs=(fcw(1)-fcw(nfcw))/2.+fcw(1)
c     sd=zsr(mzsrc(1))
      if (nfrq>1) then
         nf1=MAX(1,NINT(fcw(1)/(fcw(2)-fcw(1)))) + 1
      else
         nf1=1
      endif
      mx1=nf1+nfrq-1

      if (iforwt==1) then
      write(*,*)'Calling trfhead'
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     rkm(1),rspace,nfftbb,nf1,mx1,dt1,freqs,sd,
     &     bintrf,ndep,1,nran)
      endif
c      write(*,*) np,ndep,nfrq
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
c         DO ifreq=1,nfrq
c            DO i=1,np           ! ranges
c               DO id=1,ndep     ! irin
c                  index=(id +((ifreq-1))*ndep-1)*np
c                  write(*,*)'writing trf..'
                  WRITE(luttrf,*)
                  WRITE(luttrf,'(1000000f15.8)')
     &    ((( REAL(resp(ii+(id +((ifreq-1))*ndep-1)*np)), 
     &        imag(resp(ii+(id +((ifreq-1))*ndep-1)*np))
     &        ,id=1,ndep),ii=1,np), ifreq=1,nfrq)
                  WRITE(34)
     &    ((( REAL(resp(ii+(id +((ifreq-1))*ndep-1)*np)), 
     &        imag(resp(ii+(id +((ifreq-1))*ndep-1)*np))
     &        ,id=1,ndep),ii=1,np), ifreq=1,nfrq)

      RETURN
      END  !writecomplex
ccc
c
