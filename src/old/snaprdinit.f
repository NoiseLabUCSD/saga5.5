c234567890123456789012345678901234567890123456789012345678901234567890
c        1         2         3         4         5         6         7
      SUBROUTINE forwardmodel(iopt,mopt)
      INTEGER  mopt,i,iopt(mopt)
      DO i=1,mopt
        iopt(i)=0
      ENDDO
      iopt(1)=2
      iopt(12)=1    ! using three indexes for adressing variable
      iopt(30)=2
      END
c*****************************************************
      SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization PROGRAM
c     PETER GERSTOFT, 1992
c
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER m,i,j,ierr,jj,irflag,ndum
      REAL rdstep
      REAL*8  ztemp
      INTEGER isect
      REAL  fmindum,fmaxdum,df_temp
      INTEGER nf1,nf2
      CHARACTER*80 dumch2, dumch
      EQUIVALENCE (dumch2,opt)
      INCLUDE 'comsnaprd.h'
      ierr=0
      nfrq=1
c---  READ the GA-parameters
      CALL readinputstart
     
      READ(1,'(a)') dumch2      ! This line has no effect in snap
      CALL snapoption(opt)

      WRITE(*,*)'reading nfrq..'
      READ(1,*,err=91) nfrq, maxnomode
      GOTO 92
91    BACKSPACE(1)
      WRITE(*,*)' maxmode not defined, reading just nfreq'
      READ(1,*) nfrq
      maxnomode=100
92    WRITE(*,'(a,i3,a,i4)') ' nfrq=', nfrq,' maxmode=' ,maxnomode
      IF (nfrq.GT.mfreq) STOP ' nfreq > mfreq' ! test

      IF (nfrq.GT.0) THEN
         READ(1,*)(frq(I),i=1,nfrq)
         WRITE(*,'(a,10f10.3)')' frequencies:',(frq(I),i=1,nfrq)
      ELSE
         READ(1,*)df_temp, fmindum, fmaxdum
         nf1=MAX(1,NINT(fmindum/df_temp)) + 1
         nf2=NINT(fmaxdum/df_temp) + 1
         fmindum=(nf1 - 1)*df_temp
         fmaxdum=(nf2 - 1)*df_temp
         nfrq=nf2 - nf1 + 1
         WRITE(*,*)' Fmin, Fmax, DF, nfrq: ', 
     .        Fmindum, Fmaxdum, DF_temp,nfrq
         IF (nfrq.GT.mfreq) STOP 'increase mfreq! '
         
         DO j=1,nfrq
            frq(j)=fmindum + (j-1)*df_temp
         ENDDO           
      ENDIF
      
      READ(1,*)nsect
      WRITE(*,*)' number of sectors',nsect
      IF (nsect.GT.msec) THEN
         WRITE(*,*) ' Max number of sectors exceed'
         WRITE(*,*) 'nsect,msec=',nsect,msec
      ENDIF
c---  now we are reading in a snap FORMAT
c
c---  loop over sectors
c     
      DO isect=1,nsect
         WRITE(*,*) ' SECTOR No ',isect
         READ(1,*) R_SLNGTH(isect)
         R_SLNGTH(isect) = 1000.*R_SLNGTH(isect)
         WRITE(*,*) ' sector length (m):',R_slngth(isect) 
         READ(1,*) ztemp,R_scatt(1,isect),
     1        R_scatt(2,isect),R_beta(0,isect) 
         WRITE(*,*)' hw=',ztemp
         R_h0(isect)=ztemp
         WRITE(*,*)' Water speed profile:'
         DO m = 1,maxdep
            READ(1,*) R_z0(m,isect),R_c0(m,isect)
            WRITE(*,*) m,R_z0(m,isect),R_c0(m,isect)
            IF (R_z0(m,isect).EQ.ztemp) GOTO 101
         ENDDO
         STOP 'maxdep points read in for z0'
         
 101     nwater=m
         R_nd0(isect)=m
         WRITE(prtfil,*)' Points in water:',R_nd0(isect)
         READ(1,*)R_h1(isect),R_r1(isect),R_beta(1,isect) 
         WRITE(*,*)' h_sed =',R_h1(isect)
         IF(R_H1(isect).GT.0.0)   THEN
            WRITE(*,*)' Sediment speed profile:'
            DO m=1,maxdep
               READ(1,*)R_z1(m,isect),R_c1(m,isect)
               WRITE(*,*)m,R_z1(m,isect),R_c1(m,isect)
               IF (R_z1(m,isect).EQ.R_h1(isect)) GOTO 102
            ENDDO
            STOP 'maxdep points read in for z1'
 102        R_nd1(isect)=m
         ELSE                   ! no sediment
            R_ND1(isect)=0
            R_BETA(1,isect)=0.0
            R_R1(isect)=0.0
         ENDIF
         WRITE(*,*)' Points in sediment:',R_nd1(isect)
         WRITE(*,*)' Bottom parameters:'
         READ(1,*)R_r2(isect),R_beta(2,isect),R_c2(isect)
         WRITE(*,'(a,f7.2,a,f8.4,a,f10.2)')'  Dens:',
     1        R_r2(isect),' P-att:',R_beta(2,isect),' P-vel:',
     2        R_c2(isect)
         READ(1,*)R_beta(3,isect),R_c2s(isect) 
         WRITE(*,'(a,f8.4,a,f9.2)')
     1        '  S-att:',R_beta(3,isect),' S-vel:',R_c2s(isect)
         WRITE(*,*)
      ENDDO                     ! nsect
c---  for several source locations
      READ(1,*) sd         
c     
c     receiver DATA
c     
      IF (tilt) THEN
         READ(1,*) rd,rdlow,ndep,dtilt
      ELSE
         READ(1,*) rd,rdlow,ndep
      ENDIF
      IF (ndep.GT.mdep) THEN
         WRITE(*,*) 'ndep must grater than mdep',ndep,mdep
         ierr=1
      ENDIF
      IF (rd.EQ.rdlow) THEN
         rdlow=rdlow*1.00001   
      ENDIF
      irflag=ndep
      ndep=ABS(ndep)
      IF (ndep.GT.1) THEN
         rdstep=(rdlow-rd)/float(ndep-1)
      ELSE
         rdstep=1.
      ENDIF
      IF (irflag.GT.0) THEN
c     equidistant receiver depths
         DO 920 jj=1,ndep
 920        rdep(jj)=(jj-1)*rdstep+rd
         ELSE
            WRITE(*,*)'Only equidistant receivers supported'
            STOP
c     READ in individual receiver depths
c     READ(1,*) (rdep(jj),jj=1,ndep)
         ENDIF
         WRITE(prtfil,*)'number of receiver depths read',ndep


c---  READ the range-BLOCK
         CALL  read_range_inp(ierr,ndum)

c--- should the  EOF be READ ? 

         IF (iopt(17).EQ.1) THEN
            CALL eofinit
         ENDIF

c*** create text strings for physical parameters
         phystxt(1) = 'Water depth (m) $'
         phystxt(2) = 'Water speed (m/s)$'
         phystxt(3) = 'Sediment sound speed (m/s)$ '
         phystxt(4) = 'Attenuation (dB/lambda)$'
         phystxt(5) = 'Surface roughness (m)'
         phystxt(6) = 'Sediment Density (g/cm3) $'
         phystxt(8) = 'Source depth (m) $'
         phystxt(9) = 'Source range (km) $'
         phystxt(11)= 'Shape coefficient $'
         phystxt(12)= 'Bottom speed (m/s)$'
         phystxt(13)= 'Bottom S-speed (m/s)$'
         phystxt(14)= 'Bottom density (g/cm3)$'
         phystxt(15)= 'Receiver depth (m)$'
         phystxt(16)= 'Depth of water speed pt.(m)$'
         phystxt(17)= 'Sediment depth (m)$'
         phystxt(18)= 'Depth of sed. speed pt.(m)$'
         phystxt(19)= 'Array tilt (m)$'
         phystxt(27)= 'Sector lenght (km)$'
c***  

         DO 8 j=1,mphys
       phystxt2(j)='                                               '
    6    DO 7 I=40,1,-1
            IF(phystxt(j)(I:I).NE.'$') GO TO 7
            phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c     WRITE(*,*)phystxt2(j)
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
         IF (par2phy(i).LT.1 .OR. par2phy(i).GT.19 .OR.
     &        par2phy(i).EQ.10 )THEN
            WRITE(*,*)' *** par2phy not correct, parm',i
            ierr=1
         ENDIF
         IF (par2phy(i).EQ.1 )THEN
            IF (fmin(i).LT.R_z0(R_nd0(par3(i))-1,par3(i))) THEN
               WRITE(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** During optimization'
               WRITE(*,*)' *** the water-depth can get ',
     &              'above the second-last'
               WRITE(*,*)' *** point in the water'
               ierr=1
            ENDIF
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
            IF (fmin(i).LT.R_z1(R_nd1(par3(i))-1,par3(i))) THEN
               WRITE(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** Durring optimization'
               WRITE(*,*)' *** the sediment-depth can get above the'
               WRITE(*,*)' *** second-last point in the sediment'
               ierr=1
            ENDIF
         ENDIF
      ENDDO

c*** errors ?
      IF (ierr.EQ.1)STOP 'stopped in input'
c***  CLOSE input unit
      CLOSE(1)                                
      END
      SUBROUTINE snapoption(options)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnaprd.h'
      CHARACTER options(40)
      INTEGER i
      tilt=.FALSE.
      incoh=.FALSE.
      DO 50 i=1,40
         IF (options(i).EQ.'i') THEN
c     flag(6)=1
c     WRITE(*,*)'NO incoherent addition of modes  for SNAPRD !'
            WRITE(*,*)' incoherent addition of modes  for SNAPRD'
            incoh=.TRUE.
c     STOP
c     IF (iopt(13).EQ.1) THEN
c     WRITE(*,*)'************* WARNING******************'
c     WRITE(*,*)' incoherent adition of modes does not give'
c     WRITE(*,*)' the phase of the pressure'
c     WRITE(prtfil,*)'************* WARNING******************'
c     WRITE(prtfil,*)' incoherent adition of modes does not give'
c     WRITE(prtfil,*)' the phase of the pressure'
c     ENDIF
         ELSE IF (options(i).EQ.'t') THEN
            tilt=.TRUE.
            WRITE(*,*)' tilt of array is included'
            WRITE(prtfil,*)' tilt of array is included'
            IF (iopt(13).EQ.0 .OR. nx.GT.1) THEN
               WRITE(*,*)' the tilt can only be computed for one range'
            ENDIF
         ELSE IF (options(I).EQ.'!') THEN
            GOTO 60
         ELSE IF (options(I).NE.' ') THEN
            WRITE(prtfil,399) options(I)
 399        FORMAT(1H ,' >>>> unknown snap OPTION: ',A1,' <<<<')
         END IF
 50   CONTINUE
 60   CONTINUE
      END

c**************************************************************
      SUBROUTINE forwinit
c SNAP SUBROUTINE to INTERFACE to snap- initialization of parameters
      USE global
      INTEGER i,ifreqy,id,index 
      INTEGER isect
      DOUBLE PRECISION ccmean
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnaprd.h'
      DOUBLE PRECISION  R_CC0(msec), R_CC1(msec)
      WRITE(*,*)' uses snapRD as forward model'
      WRITE(prtfil,*)' uses snapRD as forward model'

c       IF (maxnomode.GT.moden) THEN
c         WRITE(*,*)' Warning maxnumber of modes reduced to',moden
c         WRITE(*,*)' Snap does not accommodate more modes'
c         maxnomode=moden
c       ENDIF

c---  transfer water properties
      DO isect=1,nsect
         R_r0(isect)=1.0        ! density in water
      ENDDO	
c---  transfer receiver parametes
      fldrd(1)= rd		      ! minimum receiver depth
      fldrd(2)= rdlow		      ! maximum receiver depth
      IF (ndep.EQ.1) THEN
        fldrd(3)= (rdlow-rd)          ! step in receiver depth
      ELSE
        fldrd(3)= (rdlow-rd)/(ndep-1)   ! step in receiver depth
      ENDIF
      sddum=sd
      DO i=1,ndep
c        srd(i,1)=sd
        rddum(i)=rd+fldrd(3)*(i-1)
      ENDDO 
      nrd=ndep
c---  range parameters
      np= nx                           ! number of range steps
      secd(1)=xranges(1)
      secd(2)=xranges(np)
      IF (np.GT.1) THEN
        secd(3)=(xranges(np)-xranges(1))/(np-1)
      ELSE
        secd(3)=0.
      ENDIF   
      secd(2)=secd(2)+0.1*secd(3)
      

       WRITE(*,*)' min and max range from snaprdinit',secd
c      WRITE(*,*)'fldrd from snaprdinit',rddum(1),rddum(2)

      flagpu=-2

      IF (iopt(6).EQ.1) flagpu=-200
      CALL file(prtfil)

c***************************************************
      ENTRY forw2()
c***************************************************
       lwascomp=1
c      WRITE(*,*) 'has entered snaprd'
      flagpu=1 +flagpu
c
c---- USE of EOF
c
      IF (iopt(17).EQ.1) THEN
        CALL eofval
        IF (lwascomp.EQ.-1) THEN
        RETURN
        ENDIF
      ENDIF
c
      IF (np.EQ.1) secd(2)=secd(1)
      DO isect=1,nsect
        R_z0(R_nd0(isect),isect)=R_h0(isect)      ! last point is the water depth
        IF (R_nd1(isect).GT.0) R_z1(R_nd1(isect),isect)=R_h1(isect)
      ENDDO
C***********

C     NORMALIZATION OF DEPTHS
C
      DO 3000 isect=1,nsect
         ccmean=0
         DO I=1,R_ND0(isect)-1
            ccmean= ccmean + 
     1           (R_c0(i+1,isect)+R_c0(i,isect))
     2           *.5*(R_z0(i+1,isect)-R_z0(i,isect))
         ENDDO
         R_cc0(isect)=ccmean/R_h0(isect) ! mean velocity
         
         ccmean=0
         IF (R_h1(isect).GT.0) THEN
            DO 1500   I=1,R_ND1(isect)
 1500       CONTINUE
            DO I=1,R_ND1(isect)-1
               ccmean= ccmean + 
     1              (R_c1(i+1,isect)+R_c1(i,isect))
     2              *.5*(R_z1(i+1,isect)-R_z1(i,isect))
            ENDDO
            R_cc1(isect)=ccmean/R_h1(isect) ! mean velocity sed
         ENDIF
 3000 CONTINUE
      R_beta(1,isect)=ABS( R_beta(1,isect))
      R_beta(2,isect)=ABS( R_beta(2,isect))

      DO ifreqy=1,nfrq
         IF (ifreqy.GE.2)flagpu=flagpu+3
         MINMOD=1
         MAXMOD=maxnomode
         CALL snaprd(nsect,frq(ifreqy),maxmod,minmod,R_slngth, 
     &        R_R0, R_R1, R_R2, R_BETA, R_SCATT, R_C2S, R_C2,
     &        R_C0, R_Z0, R_C1, R_Z1, 
     &        R_CC0, R_CC1, R_H0, R_H1, R_ND0, R_ND1,
     &        msec,maxdep,pressr,np,xranges)
         
         DO id=1,ndep           ! irin
            index=(id +((ifreqy-1))*ndep-1)*np
            DO i=1,np           ! ranges
               resp(i+index)=pressr(i+(id-1)*np)
            ENDDO               ! ranges
         ENDDO                  ! depth
c     
c---- WRITE out
         IF (flagpu.LT.-2)THEN
            WRITE(prtfil,*)'range,depth',sddum,secd(1)
            WRITE(prtfil,*)' The pressures are:'
            DO id=1,ndep        ! irin
               index=(id +((ifreqy-1))*ndep-1)*np
               DO i=1,np        ! ranges
                  WRITE(prtfil,*)i+index,pressr(i+(id-1)*np),i+(id-1)*np
               ENDDO            ! ranges
            ENDDO               ! depth
         ENDIF
         IF (ifreqy.GE.2)flagpu=flagpu-3
      ENDDO                     ! ifreq 
      
      IF (ABS(resp(np+index)).LE. 1.e-20) THEN 
         WRITE(*,*)' the pressure of the last point is less than'
         WRITE(*,*)'  -400dB; '
         WRITE(*,*)' This is probably due to a too low sound '
         WRITE(*,*)'speed in the'
         WRITE(*,*)' sediment or bottom.'
c     WRITE(*,*)' a dump follows '
         WRITE(prtfil,*)' the pressure of the last point is' 
         WRITE(prtfil,*)' less than -400dB; '
         WRITE(prtfil,*)' This is probably due to a too low sound '
         WRITE(prtfil,*)'speed in the'
         WRITE(prtfil,*)' sediment or bottom.'
c     WRITE(prtfil,*)' a dump follows: '
         
      ENDIF

c     IF (flagpu.EQ.1) CLOSE(0)
      lwascomp=1
      RETURN
      END       

c*********************************************
      SUBROUTINE writeTL()
      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnaprd.h'

      INTEGER ifreq,i,id,index,mx1,nf1
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
         rspace= drng ! rng(2)-rng(1)
      END IF
      mx1=1 !nf1+nfrq-1
      dt1=1 !1.E0/fsbb
      freqs= 1             !(fmaxdum-fmindum)/2.+fmindum
c     sd=zsr(mzsrc(1))
      if (iforwt==1) then
      write(*,*)'Calling trfhead'
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     xranges(1),rspace,nfrq,nf1,mx1,dt1,freqs,sd,
     &     bintrf,ndep,1,np)
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
     &    ((((20*log10(abs(resp(i+(id +((ifreq-1))*ndep-1)*np))))
     &        ,id=1,ndep),i=1,np), ifreq=1,nfrq)
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
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnaprd.h'

      INTEGER ifreq,i,id,index,mx1,nf1
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
         rspace=drng
      END IF
      mx1=1 !nf1+nfrq-1
      dt1=1 !1.E0/fsbb
      freqs=1 !(fmaxdum-fmindum)/2.+fmindum
c     sd=zsr(mzsrc(1))
      if (iforwt==1) then
      write(*,*)'Calling trfhead'
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     xranges(1),rspace,nfrq,nf1,mx1,dt1,freqs,sd,
     &     bintrf,ndep,1,np)
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
                  WRITE(luttrf,'(1000000f15.8)')
     &    ((( REAL(resp(i+(id +((ifreq-1))*ndep-1)*np)), 
     &        imag(resp(i+(id +((ifreq-1))*ndep-1)*np))
     &        ,id=1,ndep),i=1,np), ifreq=1,nfrq)
                  WRITE(34)
     &    ((( REAL(resp(i+(id +((ifreq-1))*ndep-1)*np)), 
     &        imag(resp(i+(id +((ifreq-1))*ndep-1)*np))
     &        ,id=1,ndep),i=1,np), ifreq=1,nfrq)

      RETURN
      END  !writecomplex
ccc
c     




