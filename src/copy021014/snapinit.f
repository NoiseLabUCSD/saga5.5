c     234567890123456789012345678901234567890123456789012345678901234567890
c     1         2         3         4         5         6         7
      SUBROUTINE forwardmodel(iopt,mopt)
      INTEGER  mopt,i,iopt(mopt)
      DO i=1,mopt
         iopt(i)=0.
      ENDDO
      iopt(1)=2
      iopt(30)=1                ! 1 is for snap.
      END
c******************************
      SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization PROGRAM
c     PETER GERSTOFT, 1992
c     
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INTEGER iii,m,i,j,ierr,jj,irflag,ndum,inu
      REAL rdstep
      REAL*8 ztemp
      REAL sdlow
      INTEGER itemp,icount,ierrinp
      CHARACTER*80 dumch2, dumch
      EQUIVALENCE (dumch2,opt)
      INCLUDE 'comsnap.h' 
      inu=0                     !cfh
      ierr=0
      nfrq=1
c---  READ the GA-parameters
      CALL readinputstart
      
      READ(1,'(a)') dumch2      ! This line has no effect in snap
      CALL snapoption(opt)

      WRITE(*,*) 'reading nfrq...'
      READ(1,*,err=91) nfrq, maxnomode
      IF (maxnomode.LE.0) THEN
         maxnomode=500
         WRITE(*,*)' Max number of modes < 0; changed to 500'
      ENDIF
      GOTO 92
 91   BACKSPACE(1)
      READ(1,*,err=91) nfrq
      maxnomode=500
 92   WRITE(*,'(a,i3,a,i4)') ' nfrq=', nfrq,' maxmode=' ,maxnomode
c     .cfh.
      IF (iopt(5).EQ.5) THEN
         del_time = 0.0 
         WRITE(*,'(a,i3,a,i4,a,10f10.3)') ' nfrq=', nfrq,
     &        ' maxmode=',maxnomode,' time delay=',del_time
      ENDIF


      IF (ABS(nfrq).GT.mfreq) STOP ' nfreq > mfreq' ! test
c      WRITE(*,*)'hello snapinit0'
      
      IF (nfrq.GT.0) THEN
         READ(1,*)(frq(I),i=1,nfrq)
         WRITE(*,'(a,10f10.3)')' frequencies:',(frq(I),i=1,nfrq)
      ELSE
         iifft =1
 10      WRITE(*,*)'reading frequencies...'
         READ(1,'(a80)')dumch
         READ(dumch,*,err=11)
     &        df_temp, fmindum, fmaxdum,dflast,fminlast,fmaxlast
         GOTO 12
 11      CONTINUE
c     BACKSPACE(1)
         READ(dumch,*,err=10)df_temp, fmindum, fmaxdum
         dflast   = df_temp
         fminlast = fmindum
         fmaxlast = fmaxdum
 12      CONTINUE
c     WRITE(*,*)'hello snapinit'
         nf1=MAX(1,NINT(fmindum/df_temp)) + 1
         nf2=NINT(fmaxdum/df_temp) + 1
         fmindum=(nf1 - 1)*df_temp
         fmaxdum=(nf2 - 1)*df_temp
c     old         nfbb=nf2 - nf1 + 1
         nfrq=nf2 - nf1 + 1
         nfftbb=2**(INT(LOG10(float(nf2))/LOG10(2.) + .9))*2
         fsbb= dfloat(nfftbb)*df_temp
c     old         READ(1,*)df_temp, fmindum, fmaxdum
c     old         nf1=MAX(1,NINT(fmindum/df_temp)) + 1
c     old         nf2=NINT(fmaxdum/df_temp) + 1
c     old         fmindum=(nf1 - 1)*df_temp
c     old         fmaxdum=(nf2 - 1)*df_temp
c     old         nfrq=nf2 - nf1 + 1
         WRITE(*,*)' Fmin, Fmax, DF, nfrq: ', 
     .        Fmindum, Fmaxdum, DF_temp,nfrq
         IF (nfrq.GT.mfreq) STOP 'increase mfreq! '
         
         DO j=1,nfrq
            frq(j)=fmindum + (j-1)*df_temp
         ENDDO           
      ENDIF

c---  now we are reading in a snap FORMAT
      READ(1,*)ztemp,scatt(1),scatt(2),beta(0) 
      WRITE(*,*)' hw=',ztemp
      hw=ztemp
      WRITE(*,*)' Water speed profile:'
      IF (beta(0).GE.0.1) THEN
         beta(0)=0
         WRITE(*,*)'att in water >0.1; put to zero !!'
      ENDIF
      DO m=1,maxdep
         READ(1,*)z0(m),c0(m)
         WRITE(*,*)m,z0(m),c0(m)
         IF (z0(m).EQ.ztemp) GOTO 101
      ENDDO
      STOP 'maxdep points read in for z0'
      
 101  nwater=m
      nd0=m
      IF ((nd0.EQ.1) .OR. (z0(1).NE.0) ) THEN
         WRITE(*,*)'*** Snap requires a sound speed',
     1        ' point at start of ocean  (z=0)'
         STOP
      ENDIF 
      WRITE(prtfil,*)' Points in water:',nd0
      READ(1,*)h1,r1,beta(1) 
c     IF (beta(1).LE.0
      WRITE(*,*)' h_sed =',h1
      IF(H1.GT.0.0)   THEN
         WRITE(*,*)' Sediment speed profile:'
         DO m=1,maxdep
            READ(1,*)z1(m),c1(m)
            WRITE(*,*)m,z1(m),c1(m)
            IF (z1(m).EQ.h1) GOTO 102
         ENDDO
         STOP 'maxdep points read in for z1'
 102     nd1=m
         IF ((nd1.EQ.1) .OR. (z1(1).NE.0) ) THEN
            WRITE(*,*)'*** Snap requires a sound speed',
     1           ' point at start of sediment (z=0)'
            STOP
         ENDIF
      ELSE                      ! no sediment
         ND1=0
         BETA(1)=0.0
         R1=0.0
      ENDIF
      WRITE(*,*)' Points in sediment:',nd1
      
      WRITE(*,*)' Bottom parameters:'
      READ(1,*)r2,beta(2),c2
      WRITE(*,'(a,f7.2,a,f8.4,a,f10.2)')
     1     '  Dens:',r2,' P-att:',beta(2),' P-vel:',c2
      READ(1,*)beta(3),c2s 
      WRITE(*,'(a,f8.4,a,f9.2)')
     1     '  S-att:',beta(3),' S-vel:',c2s
      WRITE(*,*)
      
      IF (igrain) THEN
         READ(1,*)  phim        !,  max_dep 
         WRITE(*,*)'read  phim ', phim
      ENDIF
      WRITE(*,*)'MSOURCE',msource
      IF (msource) THEN
         READ(1,*) sd,sdlow,nsrd
         IF (nsrd.GT.msrd) THEN
            WRITE(*,*) ' ndep must greater than mdep',ndep,mdep
            ierr=1
         ENDIF
         IF (sd.EQ.sdlow) THEN
            sdlow=sdlow*1.00001   
         ENDIF
         irflag=nsrd
         nsrd=ABS(nsrd)
         IF (nsrd.GT.1) THEN
            rdstep=(sdlow-sd)/float(nsrd-1)
         ELSE
            rdstep=1.
         ENDIF
         IF (irflag.GT.0) THEN
c     equidistant receiver depths
            DO 921 jj=1,nsrd
 921           sdep(jj)=(jj-1)*rdstep+sd
            ELSE
c     READ in individual receiver depths
               READ(1,*) (sdep(jj),jj=1,nsrd)
            ENDIF
            WRITE(prtfil,*)' number of source depths read',nsrd
            IF (msourcesp) THEN
               WRITE(*,*)' ... reading source spectrum'
               CALL opfilr(42,ierrinp)
               IF (ierrinp.NE.0) THEN
                  WRITE(*,*)' The source file *.sou does not exist'
                  STOP
               ENDIF 
               CALL read_sp(ierr)
            ENDIF
         ELSE
            READ(1,*) sd         
            sdep(1)=sd
            nsrd=1
         ENDIF
         WRITE(*,*)'source depths:',(sdep(jj),jj=1,nsrd)
c     
c---  for multistatic
c     
         IF (multistatic) THEN
            WRITE(*,*)'reading multistatic repeater range and depth' 
            READ(1,*) multigeo(1),multigeo(2)
            WRITE(*,*)'repeater range,depth', multigeo(1),multigeo(2)
         ENDIF
c     
c---  for multiple beams
c     
         IF (mbeam) THEN
            WRITE(*,*)
     1           'reading 2nd source range and depth,angle, rel-power' 
            READ(1,*) multigeo(1),multigeo(2),multigeo(3),multigeo(4)
            WRITE(*,*)'2nd source range,depth ',multigeo(1),multigeo(2)
            WRITE(*,*)'2nd source angle,power ',multigeo(3),multigeo(4)
            WRITE(prtfil,*)
     1           '2nd source range,depth ',multigeo(1),multigeo(2)
            WRITE(prtfil,*)
     2           '2nd source angle,power ',multigeo(3),multigeo(4)
         ENDIF
c     
c     receiver DATA
c     
         IF (tilt) THEN
            READ(1,*) rd,rdlow,ndep,dtilt
         ELSE
            READ(1,*) rd,rdlow,ndep
         ENDIF
c AARON
         if (offset) then
            IF (abs(ndep).GT.256) THEN 
               WRITE(*,*)'dr_ind cannot > 256 elements'
               ierr=1
            ENDIF 
            Write(*,*)' Reading range offsets'
            READ(1,*) (dr_ind(jj),jj=1,(abs(ndep)))
            WRITE(prtfil,*)' number of range offsets read',abs(ndep)
         endif
         if (timeoffset) then
            IF (abs(ndep).GT.256) THEN 
               WRITE(*,*)'time_ind cannot > 256 elements'
               ierr=1
            ENDIF 
            Write(*,*)' Reading range offsets'
            READ(1,*) (time_ind(jj),jj=1,(abs(ndep)))
            WRITE(prtfil,*)' number of range offsets read',abs(ndep)
         endif
  
         irflagg=ndep
         ndep=ABS(ndep)
         IF (ndep.GT.mdep) THEN
            WRITE(*,*) ' ndep must grater than mdep',ndep,mdep
            ierr=1
         ENDIF
         IF (rd.EQ.rdlow) THEN
            rdlow=rdlow*1.00001   
         ENDIF
         IF (ndep.GT.1) THEN
            rdstep=(rdlow-rd)/float(ndep-1)
         ELSE
            rdstep=1.
         ENDIF
         IF (irflagg.GT.0 ) THEN
c     equidistant receiver depths
            DO 920 jj=1,ndep
 920           rdep(jj)=(jj-1)*rdstep+rd
            ELSE
c     READ in individual receiver depths
               WRITE(*,*)' Reading individual receiver depths'
               READ(1,*) (rdep(jj),jj=1,ndep)
            ENDIF
            WRITE(prtfil,*)' number of receiver depths read',ndep
            DO  jj=1,ndep
               rdref(jj)=rdep(jj)
            ENDDO
            IF (out_plane) THEN        
               WRITE(*,*)' Reading array coordinates'
               IF (xplane.GE.1) THEN
                  DO jj=1,ndep 
                     zhor(jj)=rdep(jj)
                     READ(1,*)xhor(jj),yhor(jj)
                  ENDDO
               ENDIF            !xplane=1
               READ(1,*)(arrayshape(jj),jj=1,2)
               READ(1,*)(arrayshape(jj),jj=3,5)
               WRITE(*,*)' initial array param:',
     .              (arrayshape(jj),jj=1,5)
               IF (xplane.EQ.0) THEN
                  IF (ndep.GE.2) THEN
                     rdstep=arrayshape(1)/float(ndep-1)
                  ELSE
                     rdstep=1.
                  ENDIF
c     equidistant receiver array
                  DO  jj=1,ndep
                     rdep(jj)=(jj-1)*rdstep+rd
                  ENDDO
               ENDIF            !xplane=0
               DO  jj=1,ndep
                  zhor(jj)=rdep(jj)
               ENDDO
            ENDIF               ! out_plane
            DO 931 jj=1,ndep
 931           rdref(jj)=rdep(jj)
c---  READ the range-BLOCK
               CALL read_range_inp(ierr,ndum)
               
c---  should the EOF be READ ? 
               IF (iopt(17).EQ.1) THEN
                  CALL eofinit
               ENDIF
               
c***  create text strings for physical parameters
c     123456789012345678901234567890
               phystxt(1) = 'Water depth (m) $'
               phystxt(2) = 'Water sound speed (m/s)$'
               phystxt(3) = 'Sediment sound speed (m/s)$'
               phystxt(4) = 'Attenuation (dB/\lambda)$'
               phystxt(5) = 'Surface roughness (m)'
               phystxt(6) = 'Sediment density (g/cm^3)$'
               phystxt(8) = 'Source depth (m) $'
               phystxt(9) = 'Source range (km) $'
               phystxt(11)= 'Shape coefficient $'
               phystxt(12)= 'Bottom sound speed (m/s)$'
               phystxt(13)= 'Bottom S-speed (m/s)$'
               phystxt(14)= 'Bottom density (g/cm^3)$'
               phystxt(15)= 'Receiver depth (m)$'
               phystxt(16)= 'Depth of water speed pt.(m)$'
               phystxt(17)= 'Sediment depth (m)$'
               phystxt(18)= 'Depth of sed. speed pt.(m)$'
               phystxt(19)= 'Array tilt (m)$'
               phystxt(20)= 'Sound speed profile # $'
               phystxt(21)= 'Array out of plane param$'
               phystxt(22)= 'Source strength $'
               phystxt(23)= 'Repeater geometry$'
               phystxt(24)= 'Mean grain size (\phi)$'
c AARON
               phystxt(25)= 'Range offset (m)$'
               phystxt(26)= 'Time offset (s)$'
c     .cfh. for matched filter
               phystxt(28)= 'Time delay $'
               phystxt(29)= 'Error variance$'
c***  
               
               DO 8 j = 1,mphys
                  phystxt2(j)='                                       '
 6                DO 7 I= 40, 1, -1
                     IF (phystxt(j)(I:I).NE.'$') GO TO 7
                     phystxt2(j)(1:I-1) = phystxt(j)(1:I-1)
c     WRITE(*,*) phystxt2(j)
                     GO TO 8
 7                CONTINUE
 8             CONTINUE
               phystxt2(9)='Source range (m)'
               
c---- READ the optimization PARAMETER
               CALL readoptparm
               ibeta=0

               DO i=1,nparm
                  IF (fmin(i).GT.fmax(i))THEN
                     WRITE(*,*)' *** fmin > fmax for parm',i
                     ierr=1
                  ENDIF
c     AARON .cfh.
                  IF (par2phy(i).LT.1 .OR. par2phy(i) .GT. 29 .OR.
     &                 par2phy(i).EQ.10 ) THEN
                     WRITE(*,*)' *** par2phy not correct, parm',i
                     ierr=1
                  ENDIF
                  IF (par2phy(i).EQ.1 )THEN
                     IF (fmin(i).LT.z0(nd0-1)) THEN
                        WRITE(*,*)' *** Optimzation variable #:',i
                        WRITE(*,*)' *** Durring optimization'
                        WRITE(*,*)' the water-depth ',
     .                       'can get above the second-last'
                        WRITE(*,*)' point in the water'
                        ierr=1
                     ENDIF
                  ENDIF
                  IF (par2phy(i).EQ.8 )THEN
                     IF (fmin(i).LE.0) THEN
                        WRITE(*,*)' *** Optimzation variable #:',i
                        WRITE(*,*)' *** Durring optimization'
                        WRITE(*,*)' the source depth ',
     .                       'can not be above the water (z=0)'
                        ierr=1
                     ENDIF
                  ENDIF
                  IF (par2phy(i).EQ.11 )THEN
                     IF (par2lay(i).GT.neof) THEN
                        WRITE(*,*)' *** Optimzation variable #:',i
                        WRITE(*,*)' *** The shapecoffient',
     .                       ' number is not defined'
                        WRITE(*,*)' par2lay(i)', par2lay(i)
                        ierr=1
                     ENDIF
                  ENDIF
                  IF (par2phy(i).EQ.11 )THEN
                     IF (par2lay(i).GT.neof) THEN
                        WRITE(*,*)' *** Optimzation variable #:',i
                        WRITE(*,*)' *** the second pointer',
     .                       ' must be less than',
     .                       ' the number of shapefunctions, 
     .                       Neof=',neof
                        WRITE(*,*)' par2lay(i)', par2lay(i)
                        ierr=1
                     ENDIF
                  ENDIF
                  IF (par2phy(i) .eq. 29 ) inu = 1
                  IF ((par2phy(i).EQ.6) .OR. (par2phy(i).EQ.8)  .OR. 
     1                 (par2phy(i).EQ.12).OR. (par2phy(i).EQ.13) .OR. 
     1                 (par2phy(i).EQ.14).OR. (par2phy(i).EQ.15) .OR.
     1                 (par2phy(i).EQ.17).OR. (par2phy(i).EQ.19)) THEN
                     itemp=par2phy(i) 
                     icount=0 
                     DO jj=1,nparm
                        IF (itemp.EQ.par2phy(jj)) icount=icount+1 
                     ENDDO
                     IF (icount.NE.1) THEN
                        WRITE(*,*)' *** Optimzation variable #:',i
                        WRITE(*,*)' *** The parameter is ',
     1                       ' defined twice as optimization parameter'
                        ierr=1
                     ENDIF
                  ENDIF
c     
                  
                  IF (par2phy(i).EQ.22)THEN
                     IF (iscale.NE.1) iscale=1
                  ENDIF
                  IF (par2phy(i).EQ.4)THEN
                     ibeta=1
                  ENDIF

                  
                  IF (par2phy(i).EQ.17)THEN
                     IF (fmin(i).LT.z1(nd1-1)) THEN
                        WRITE(*,*)' *** Optimzation variable #:',i
                        WRITE(*,*)' *** Durring optimization'
                        WRITE(*,*)' the sediment-depth can get',
     .                       ' above the second-last'
                        WRITE(*,*)' point in the sediment'
                        ierr=1
                     ENDIF
                  ENDIF
                  
                  IF (par2phy(i).EQ.20)THEN
                     OPEN(unit=44,file='soundspeed.obs',status='old')
                     READ(44,*) nobsssp, nobspts
                     WRITE(*,*) ' nobsssp=', nobsssp, ' nobspts=',
     1                    nobspts
                     DO iii=1,nobsssp
                        READ(44,*)(xobsssp(iii,j),j=1,nobspts)
                     ENDDO
                     CLOSE(44)
                     IF (nobsssp.GT.maxprof)THEN
                        ierr=1
                        WRITE(*,*)' The number of observed ssp is',
     1                       nobsssp
                        WRITE(*,*)' max dimension of ssp-matrix is',
     1                       maxprof
                     ENDIF
c     IF (fmax.NE.nobsssp) THEN
c     WRITE(*,*)' The number of observed ssp is',nobsssp
c     WRITE(*,*)' The requested ssp profiles is',par2lay(i)
c     WRITE(*,*)' These should be identical'
c     ierr=1
c     ENDIF
                     IF (nwater.LT.nobspts) THEN
                        WRITE(*,*)' sound speed points in ',
     1                       'the actual profile is'
                        WRITE(*,*)' too small relative ',
     1                       'to the observed profiles'
                        ierr=1
                     ENDIF
                  ENDIF
                  
               ENDDO
c AARON added equation 25 statement
c     .cfh. added 28 & 29 statements
               DO i=1,nparm
                  IF ((par2phy(i).EQ. 8).OR.
     1                 (par2phy(i).EQ. 9).OR.
     1                 (par2phy(i).EQ.15).OR.
     1                 (par2phy(i).EQ.19).OR.
     1                 (par2phy(i).EQ.21).OR.
     1                 (par2phy(i).EQ.22).OR.
     1                 (par2phy(i).EQ.23).OR.
     1                 (par2phy(i).EQ.25).OR.
     1                 (par2phy(i).EQ.26).OR.
     1                 (par2phy(i).EQ.28).OR.
     1                 (par2phy(i).EQ.29).OR.
     1                 (par2phy(i).EQ.1*ibathy)) THEN
                     i_geom=0   ! we DONT have to compute the modes
                  ELSE
                     i_geom=1   ! we have to compute the modes
                     GOTO 50
                  ENDIF
               ENDDO
 50            CONTINUE
               
c     i_geom=1  
               WRITE(*,*)
               WRITE(*,*)' Should modes be calculated (1=Yes,0=NO)? ',
     1              i_geom
               WRITE(*,*)
      IF (inu .eq. 1 .and. isubopt(36) .eq. 0) THEN
         WRITE(*,*)' *** Optimzation error variance #',i
         WRITE(*,*)' *** Need to specify the saga option *' 
         ierr = 1
      ELSEIF (inu .eq. 0 .and. isubopt(36) .eq. 1) THEN
         isubopt(36) = 2
         IF (iopt(4) .EQ. 4 .and. isubopt(4) .EQ. 2) ierr = 1
      ENDIF
               
c***  errors ?
               IF (ierr.EQ.1)STOP 'stopped in input'
c***  CLOSE inroptput unit
               CLOSE(1)                                    
               END
      SUBROUTINE snapoption(options)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'comsnap.h'
      REAL   FLAGopt(6)
      COMMON /flagarray/flagopt
      CHARACTER options(40)
      INTEGER i,j
      tilt      = .FALSE.
c     AARON
      offset    = .FALSE.
      timeoffset= .FALSE.
      out_plane = .FALSE.
      xplane    =  0
      msource   = .FALSE.
      mbeam     = .FALSE.
      msourcesp = .FALSE.
      multistatic= .FALSE.      
      two_way   = .FALSE.
      ibathy    = 0
      iscale    = 0
      igrain= .FALSE.
      i=0
      DO 50 j=1,40
         i=i+1
         IF (options(i).EQ.'i') THEN
            flagopt(6)=1
            WRITE(*,*)'incoherent addition of modes'
            IF (iopt(13).EQ.1) THEN
               WRITE(*,*)'************* WARNING******************'
               WRITE(*,*)' incoherent adition of modes does not give'
               WRITE(*,*)' the phase of the pressure'
               WRITE(prtfil,*)'************* WARNING******************'
               WRITE(prtfil,*)' incoherent adition of',
     .              ' modes does not give'
               WRITE(prtfil,*)' the phase of the pressure'
               PAUSE
            ENDIF
         ELSE IF (options(i).EQ.'t') THEN
            tilt=.TRUE.
            WRITE(*,*)' tilt of array is included'
            WRITE(prtfil,*)' tilt of array is included'
            IF (iopt(13).EQ.0 .OR. nx.GT.1) THEN
               WRITE(*,*)' the tilt can only be computed for one range'
            ENDIF
c     AARON
         ELSE IF (options(i).EQ.'d') THEN
            offset=.TRUE.
            WRITE(*,*)'individual receiver offsets are included'
            WRITE(prtfil,*)' individual receiver offsets are included'
            IF (iopt(13).EQ.0 .OR. nx.GT.1) THEN
               WRITE(*,*)' the offsets only compute for 1 range'
            ENDIF
          ELSE IF (options(i).EQ.'D') THEN
            timeoffset=.TRUE.
            WRITE(*,*)'individual receiver time offsets are included'
            WRITE(prtfil,*)
     1           ' individual receiver time offsets are included'
            IF (iopt(13).EQ.0 .OR. nx.GT.1) THEN
               WRITE(*,*)' the offsets only compute for 1 range'
            ENDIF
         ELSE IF (options(i).EQ.'o') THEN
            out_plane=.TRUE.
            WRITE(*,*)' out of plane array'
            WRITE(prtfil,*)' out of plane array'
         ELSE IF (options(i).EQ.'x') THEN
            out_plane=.TRUE.
            xplane=1
            WRITE(*,*)' out of plane array'
            WRITE(prtfil,*)' out of plane array'
            IF (ICHAR(options(i+1)).EQ.50) THEN ! 50 =2
               i=i+1
               xplane=2
               WRITE(*,*)' The x & z  array coordinate ',
     .              'is read from file'
               WRITE(prtfil,*)
     1              ' The x & z  array coordinate is read from file'
            ELSE IF (ICHAR(options(i+1)).EQ.51) THEN ! 51 =3
               i=i+1
               xplane=3
               WRITE(*,*)
     1              ' The x & y & z  array coordinate ',
     2              'is read from file'
               WRITE(prtfil,*)
     1              'The x & y & z  array coordinate is read from file'
            ELSE
               WRITE(*,*)' The x array coordinate is read from file'
               WRITE(prtfil,*)' The x array coordinate',
     1              ' is read from file'
            ENDIF
         ELSE IF (options(i).EQ.'m') THEN
            WRITE(*,*)' multiple sources'
            WRITE(prtfil,*)' multiple sources'
            msource=.TRUE.
         ELSE IF (options(i).EQ.'b') THEN
            WRITE(*,*)' MFP with two sources'
            WRITE(prtfil,*)' MFP with two sources'
            mbeam=.TRUE.
            isubopt(30)=1
         ELSE IF (options(i).EQ.'M') THEN
            WRITE(*,*)' multiple sources with spectrum'
            WRITE(prtfil,*)' multiple sources with spectrum'
            msource=.TRUE.
            msourcesp=.TRUE.
         ELSE IF (options(i).EQ.'h') THEN
            WRITE(*,*)' modes not recalculated for new Bathymetry '
            WRITE(prtfil,*)' modes not recalculated for new Bathymetry'
            ibathy=1
         ELSE IF (options(i).EQ.'T') THEN
            WRITE(*,*)' Two-way propagation'
            WRITE(prtfil,*)' Two-way propagation'
c     two_way=.TRUE.
            multistatic=.TRUE.
         ELSE IF (options(i).EQ.'s') THEN
            WRITE(*,*)' scaling factor multiplied on pressure '
            WRITE(prtfil,*)' scaling factor multiplied on pressure '
            iscale=1
         ELSE IF (options(i).EQ.'p') THEN
            WRITE(*,*)' Mean grain size inversion '
            WRITE(prtfil,*)'  Mean grain size inversion '
            igrain=.TRUE.
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
c     SNAP SUBROUTINE to INTERFACE to snap- initialization of parameters
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './snap/common.f'
      INCLUDE 'comsnap.h'
      INTEGER infomode
      COMMON /infomode/infomode
      DOUBLE PRECISION CC0, CC1, ZZ0(maxDEP), ZZ1(maxDEP)
      INTEGER i,ifreqy,isub,id,modqty(mfreq),index 
      DOUBLE PRECISION omeg(mfreq)
      COMMON /modes/omeg,modqty
      REAL rdstep,xhelp,pie
      INTEGER jj,ifreq
      REAL temp2,temp3,temp4,temparr(4)
      INTEGER ntemp2,ntemp3,ntemp4
      INTEGER j,ii
      DOUBLE PRECISION cw_bot 

      INCLUDE './snap/a.f'

      WRITE(*,*)' uses snap as forward model'
      WRITE(prtfil,*)' uses snap as forward model'
      ierror= 0
      IF (nfrq.GT.mfreq_sn) THEN 
         WRITE(*,*) ' Max number of frequencies in snap',mfreq_sn
         WRITE(*,*) ' is smaller than number of frequencies', nfrq
         STOP
      ENDIF
      WRITE(*,*)'maxnomode,moden',maxnomode,moden
      IF (maxnomode.GT.moden-2) THEN
         WRITE(*,*)' Warning maxnumber of modes reduced to',moden
         WRITE(*,*)' Snap does not accommodate more modes'
         maxnomode=moden-2
      ENDIF
      WRITE(*,*)'maxnomode,moden',maxnomode,moden
c---  check dimensions
c     WRITE(*,*)'nrange,ndep*nx',nrange,ndep,nx
      IF (nrange.LT.ndep*nx) THEN
         WRITE(*,*)' variable "nrange" in file saga/src/snap/common.f' 
         WRITE(*,*)' is not large enough to contain all computed '
         WRITE(*,*)' pressure for all ranges and depths'
         WRITE(*,*)' presently nrange=',nrange
         WRITE(*,*)' requested size',ndep*nx
      ENDIF
c---  transfer water properties
      
      r0=1.0                    ! density in water
      
      IF (iscale.EQ.1) THEN
         DO ifreq=1,nfrq
            pfact(ifreq)=0
         ENDDO
      ENDIF
      IF (msource) THEN
         IF (nsrd.GT.ksrd) THEN
            WRITE(*,*)'  nsrd,ksrd',nsrd,ksrd
            STOP 'source array ksrd in snap/common.f not large enough'
         ENDIF
         DO i=1,nsrd
            srd(i) = sdep(i)    ! srd(isub, up to 20 sd, not app. to field)
         ENDDO
         
      ELSE
         srd(1)=sd              ! srd(isub, up to 20 sd, not app. to field)
      ENDIF
c---  range parameters
      np= nx                    ! number of range steps
      DO i=1,np
         rng(i)= xranges(i)     ! range points
         IF (rng(i).LT.0) THEN
            WRITE(*,*)' Range less than zero for range number',i
            WRITE(prtfil,*)' Range less than zero for range number',i
         ENDIF
      ENDDO
      hstart=h0                 ! also water depth !!!  
      infomode=0
      flagpu=-0
      CORREC=-1.0
      IF (iopt(6).EQ.1) flagpu=-30
      i_call_porter=1
c     IF (i_geom.EQ.0) THEN
c     OPEN(UNIT=40,STATUS='SCRATCH',FORM='UNFORMATTED')
c     OPEN(UNIT=40,STATUS='UNKNOWN',FORM='UNFORMATTED')
c     ENDIF
c***************************************************
      ENTRY forw2()
c***************************************************
c     WRITE(*,*)'flagpu',flagpu,iWriteTrf
c     IF ((i_geom.EQ.0.) .AND. (i_call_porter.EQ.0)) THEN
c     REWIND(40)
c     ENDIF
      lwascomp=1
c     WRITE(*,*) 'has entered snap'
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
      IF  (out_plane .AND. (xplane.EQ.1)) THEN
         IF (ndep.GE.1) THEN
            rdstep=arrayshape(1)/float(ndep-1)
         ELSE
            rdstep=1.
         ENDIF
c     equidistant receiver array
         DO  jj=2,ndep
            rdep(jj)=(jj-1)*rdstep+rdep(1)
         ENDDO
      ENDIF
      

c     
      h0=hw
      hstart=h0                 ! also water depth !!!  
      z0(nd0)=h0                ! last point is the water depth
c     
      z0(nwater)=hw
      IF (nd1.GT.0) z1(nd1)=h1
      
c     
c---- optimizing for grain size
c     
      IF (igrain) THEN
         IF (iopt(6).EQ.1) WRITE(*,*)'calling PHI2PARM ', phim
         cw_bot = c0(nd0)
         z1(1) = 0.0
         DO i = 2, nd1-1
            z1(i) = 10 ** ( -1.301 + (i-2.) * 
     1           (LOG10(z1(nd1)) + 1.301) / (nd1-1.) )
         END DO
         
         CALL PHI2PARM(cw_bot, nd1, z1, phim, c1, c1s, r1_phi,
     1        alp1p, alp1s) 
         IF (iopt(6) .EQ. 1) THEN 
            DO jj = 1, nd1
               WRITE(*,'(6g12.5)') z1(jj), c1(jj), c1s(jj), 
     1              r1_phi(jj), alp1p(jj), alp1s(jj)  
            ENDDO
         ENDIF
         r1 = SUM(r1_phi)/nd1
         if (ibeta.eq.0) then
c            write(*,*)' changing beta'
            BETA(1) = SUM(alp1p)/nd1
         endif
      ENDIF
      

c     
c---- IF writing out trf FUNCTION
c     
      if (iWriteTrf.EQ.1) call writesnap
      IF (iWriteTrf.GE.1 .AND. 
     &     dflast.GT.0 .AND. fmaxlast.GT.fminlast) THEN
         flagpu=0
         df_temp=dflast
         fmindum=fminlast
         fmaxdum=fmaxlast
         nf1=MAX(1,NINT(fmindum/df_temp)) + 1
         nf2=NINT(fmaxdum/df_temp) + 1
         fmindum=(nf1 - 1)*df_temp
         fmaxdum=(nf2 - 1)*df_temp
         nfrq=nf2 - nf1 + 1
         nfftbb=2**(INT(LOG10(float(nf2))/LOG10(2.) + .9))*2
         fsbb= dfloat(nfftbb)*df_temp
         WRITE(*,*)' For optimised env Fmin, Fmax, DF, nfrq: ', 
     .        Fmindum, Fmaxdum, DF_temp,nfrq


         IF (nfrq.GT.mfreq) THEN
            WRITE(*,*)' mfreq,nfrq=',mfreq,nfrq
            STOP 'increase mfreq! '
         ENDIF 
         IF (nfrq.GT.mfreq_sn) THEN 
            WRITE(*,*)' Max number of frequencies in snap',mfreq_sn
            WRITE(*,*)' is smaller than number of frequencies', nfrq
            STOP
         ENDIF
         
         DO j=1,nfrq
            frq(j)=fmindum + (j-1)*df_temp
         ENDDO           
         i_call_porter=1
      ELSE
c     flagpu=1
      ENDIF

C***********
C     NORMALIZATION OF DEPTHS
C     
      CMINsnap=1.0D38
      cc0=0
      DO 1450   I=1,ND0
         CMINsnap=MIN(CMINsnap,C0(I))
         Zz0(I)=Z0(I)/H0
 1450 CONTINUE
      DO    I=1,ND0-1
         cc0= cc0 + (c0(i+1) + c0(i))*.5*(z0(i+1) - z0(i))
      ENDDO
      cc0=cc0/h0
      
      cc1=0
      IF (h1.GT.0) THEN
         DO 1500   I=1,ND1
            CMINsnap=MIN(CMINsnap,C1(I))
            Zz1(I)=(H0+Z1(I))/H0
 1500    CONTINUE
         DO I=1,ND1-1
            cc1= cc1 + (c1(i+1) + c1(i))*.5*(z1(i+1) - z1(i))
         ENDDO
         cc1=cc1/h1
      ENDIF
c---- positive attenuation
      DO i=1,3
         beta(i)=ABS(beta(i))
      ENDDO
      IF (flagpu.LT.2) THEN
         WRITE(prtfil,*)
         WRITE(prtfil,*)'Calling SNAP' 
         WRITE(prtfil,*)'flagpu',flagpu
         CALL writemodel
      ENDIF

c      write(*,*) 'beforeoutplane: rdep = ',(rdep(ii),ii=1,ndep)
c      write(*,*) 'beforeoutplane: rng = ',(rng(ii),ii=1,ndep)
      IF  (out_plane) THEN
         DO ii=1,ndep 
            zhor(ii)=  rdep(1)-rdref(1)+rdref(ii)
         ENDDO 
c         write(*,*) 'zhor',(zhor(ii),ii=1,ndep)
c         do irang=1,np 
            CALL geom(rng,rdep,ndep)
c         enddo
      ENDIF
c      write(*,*) 'afteroutplane: rdep = ',(rdep(ii),ii=1,ndep)
c      write(*,*) 'afteroutplane: rng = ',(rng(ii),ii=1,ndep)

      DO ifreq=1,nfrq
         IF(ifreq.GE.2)flagpu=flagpu+2
         freqy=DBLE(frq(ifreq))
         OMEGA=TWOPI*FReqy
c     WRITE(prtfil,*)'calling snap;freqy ',freqy
c     WRITE(*,*)'calling snap;freqy ',freqy
         PHVMAX=C2
         PHVMIN=0.
         MINMOD=1
         MAXMOD=maxnomode
c     i_call_porter=1  ! calling porter
c     i_geom=1         ! geometry & writing to unit 40 
c     WRITE(*,*)'to enter porter'
c     WRITE(*,*) 'i_call_porter, i_geom',i_call_porter, i_geom
         IF (i_call_porter.EQ.1) THEN
c     WRITE(*,*)'entering Porter, nmes=',nmes
            CALL PORTER(FReQy,MINMOD,MAXMOD,NMES,JF2,DELFRQ,
     &           MODPLT,EK(1,ifreq),EGV,XTS(1,1,ifreq),
     &           CC0,CC1,ALFA(1,ifreq),MODAVR,ADA,SPEED,
     &           EIGF,ISO,NBEG,MY,C0,ZZ0,C1,ZZ1,
     &           A3,B3,C3,EE,ZZ,SSOLD,EXCH,slow)
c     WRITE(*,*)'exit porter',nmes
            MODQTY(ifreq)=MAXMOD-MINMOD+1
            omeg(ifreq) =omega
            IF(MODQTY(ifreq).EQ.1)  THEN
               WRITE(*,*)' Warning: snapinit found only one mode'
            ELSEIF((MODQTY(ifreq).EQ.maxnomode)
     &              .AND.(infomode.LE.30))THEN
               WRITE(*,*)' Warning: snapinit max number', 
     &              ' of modes reached'
               infomode=infomode+1
            ELSEIF(MODQTY(ifreq).EQ.0)  THEN
               WRITE(*,*)' Warning: snapinit found no ',  
     &              'modes...returning'
               lwascomp=-1
               RETURN
            ENDIF
            IF ((flagpu.LE.20)) THEN
               WRITE(prtfil,*)'freq,number of modes',
     &              ifreq,modqty(ifreq)
               WRITE(*,*)'freq,number of modes',ifreq,modqty(ifreq)
            ENDIF 
c     IF (i_geom.EQ.0) THEN
c     WRITE(*,*)'Writing modes'
c     WRITE(40)freqy,modqty,minmod,maxmod,msp,omega
c     WRITE(40)(ek(i),i=1,modqty)
c     WRITE(40)(alfa(i),i=1,modqty)
c     WRITE(40)((xts(i,j),i=1,modqty),j=1,msp)
c     ENDIF  
c     ELSE
c     WRITE(*,*)'Reading modes'
c     READ(40)freqy,modqty,minmod,maxmod,msp,omega
c     WRITE(*,*)    'freqy,modqty,minmod,maxmod,msp,omega'      
c     WRITE(*,*)    freqy,modqty,minmod,maxmod,msp,omega     
c     READ(40)(ek(i),i=1,modqty)
c     READ(40)(alfa(i),i=1,modqty)
c     READ(40)((xts(i,j),i=1,modqty),j=1,msp)
         ENDIF

         isub=15                ! we are selecting field option     
c     
c---- multistatic
         IF (multistatic) THEN
            temp3=srd(1)
            srd(1)=multigeo(1)
            ntemp3=nsrd
            nsrd=1
         ENDIF
c---  transfer source parametes
         ifreqy=1
         minmod=1
         maxmod=modqty(ifreq)
         omega =omeg(ifreq)
         

c     WRITE(*,*)'entering field'
c         write(*,*) 'beforefield: rng = ',(rng(irang),irang=1,np)
         CALL FIELD(NP,frq(ifreq),IFREQy,MSP,XTS(1,1,ifreq),
     &        TEMPOR,RNG,FLDPR,FLDPR1,US,UR,ALFA(1,ifreq),EK(1,ifreq),
     &        TLUS,SRD, ND0,MODEN,KSRD,rdep,ndep,nsrd,nrange,
     &        sourcesp(1,ifreq),*4000)
c     WRITE(*,*)'Field was called'
         
c         write(*,*) 'afterfield: rng = ',(rng(irang),irang=1,np)

         DO id=1,ndep           ! irin
            index=(id +((ifreq-1))*ndep-1)*np
            DO i=1,np           ! ranges
c     IF(two_way) THEN
c     WRITE(6,*)'Two_way: ',two_way
c     resp(i+index)=(fldpr(i+(id-1)*np))**2
c     ELSE
               resp(i+index)=fldpr(i+(id-1)*np)
c     END IF
c     WRITE(*,*)index+i,fldpr(i+(id-1)*np)
            ENDDO               ! ranges
            
            
            IF (flagpu.LT.0) THEN
c     WRITE(prtfil,*)'Pressure for each range:'
               DO i=1,np        ! ranges
                  WRITE(prtfil,*) resp(i+index),i+index
               ENDDO            ! ranges         
            ENDIF
            IF (ABS(fldpr(np)).LE. 1.e-20) THEN 
               ierror= ierror+1
               CALL mumble
               CALL writemodel
               WRITE(prtfil,*)'number of modes',modqty(ifreq)
               DO i=1,np        ! ranges
                  WRITE(prtfil,*)rng(i),fldpr(i+(id-1)*np)
               ENDDO 
               IF (ierror.GT.10) STOP 'stoping due to too many errors'
            ENDIF               ! dump loop
         ENDDO                  ! ir
         
c     IF (ifreq.GE.2)flagpu=flagpu-3
c     
c----------multistatic
c     
         IF (multistatic ) THEN
c--   move REAL sd back
            nsrd=ntemp3
            srd(1)= temp3
c---  USE first range and 
            temp2=rng(1)
            ntemp2=np
            rng(1)=multigeo(2)
            np=1
c---  USE just one depth
            temp4=rdep(1)
            ntemp4=ndep
            rdep(1)=multigeo(1)
            ndep=1
c     WRITE(*,*)'srd(1),nsrd',srd(1),nsrd
c     WRITE(*,*)'rng(1),ndep,np',rng(1),ndep,np
c     WRITE(*,*)'Calling Field'
            
            CALL FIELD(NP,frq(ifreq),IFREQy,MSP,XTS(1,1,ifreq),
     &           TEMPOR,RNG,FLDPR,FLDPR1,US,UR,ALFA(1,ifreq),
     &           EK(1,ifreq),TLUS,SRD,ND0,MODEN,KSRD,rdep,ndep,
     &           nsrd,nrange,sourcesp(1,ifreq),*4000)
c---  move back
            rdep(1)=temp4
            rng(1)= temp2
            ndep=ntemp4
            np=ntemp2
            DO id=1,ndep        ! irin
               index=(id +((ifreq-1))*ndep-1)*np
               DO i=1,np        ! ranges
c     WRITE(*,*)index+i,resp(I+index),fldpr(1)
                  resp(i+index)=resp(i+index)*fldpr(1)
c     WRITE(*,*)'after',index+i,resp(I+index),fldpr(1)
               ENDDO            ! ranges
            ENDDO               ! depth
         ENDIF                  ! multistatic
         
         IF (mbeam ) THEN
c--   move REAL sd back
            temparr(1)=rng(1)
            temparr(2)=srd(1)
            temparr(3)=arrayshape(5)
            rng(1)=multigeo(1)
            srd(1)= multigeo(2)
            arrayshape(5)=multigeo(3)
            IF  (out_plane ) THEN
               DO ii=1,ndep
                  zhor(ii)=  rdep(1)-rdref(1)+rdref(ii)
               ENDDO
               CALL geom(rng,rdep,ndep)
            ENDIF
            
c     WRITE(*,*)'srd(1),nsrd',srd(1),nsrd
c     WRITE(*,*)'rng(1),ndep,np',rng(1),ndep,np
c     WRITE(*,*)'Calling Field'
            
            CALL FIELD(NP,frq(ifreq),IFREQy,MSP,XTS(1,1,ifreq),
     &           TEMPOR,RNG,FLDPR,FLDPR1,US,UR,ALFA(1,ifreq),
     &           EK(1,ifreq),TLUS,SRD,ND0,MODEN,KSRD,
     &           rdep,ndep,nsrd,nrange,sourcesp(1,ifreq),*4000)
c---  move back
            rng(1) =  temparr(1)
            srd(1) = temparr(2)
            arrayshape(5) = temparr(3)
            IF  (out_plane ) THEN
               DO ii=1,ndep
                  zhor(ii)=  rdep(1)-rdref(1)+rdref(ii)
               ENDDO
               CALL geom(rng,rdep,ndep)
            ENDIF
c---  get field 
            xhelp = 10.**(multigeo(4)/20.)
c     WRITE(*,*)'xhelp',xhelp
            DO id=1,ndep        ! irin
               index=(id +((ifreq-1))*ndep-1)*np
               DO i=1,np        ! ranges
c     WRITE(*,*)index+i,resp(I+index),fldpr(1)
                  resp1(i+index)=fldpr(i+(id-1)*np)
               ENDDO            ! ranges
            ENDDO               ! depth
         ENDIF                  ! mbeam
         
         
      ENDDO                     ! loop over frequencies
c     
c---- multiply WITH scaling factor
c     
      IF (iscale.EQ.1) THEN
c     WRITE(*,*)' scaling of pressure',iscale,pfact
         DO ifreq=1,nfrq
            xhelp=10.**(pfact(ifreq)/20.)
            DO id=1,ndep        ! irin
               index=(id +((ifreq-1))*ndep-1)*np
               IF(two_way) THEN
c     WRITE(6,*)'Two_way: ',two_way
                  DO i=1,np     ! ranges
                     resp(i+index)=(xhelp*fldpr(i+(id-1)*np))**2
                  ENDDO         ! ranges
               ELSE
                  DO i=1,np     ! ranges
                     resp(i+index)=xhelp*fldpr(i+(id-1)*np)
                  ENDDO         ! ranges
               END IF
            ENDDO               ! depth
         ENDDO                  ! freq
      ENDIF                     ! iscale loop
c     .cfh.
c---- multiply WITH delta time
c     
      IF (iopt(5).EQ.5) THEN
         pie = 4.*ATAN(-1.)
         DO ifreq=1,nfrq
c            WRITE(6,*)'frq', frq(ifreq)
            DO id=1,ndep        ! irin
               index=(id +((ifreq-1))*ndep-1)*np
               DO i=1,np        ! ranges
                  resp(i+index) = resp(i+index)
     .                 *cexp(CMPLX(0.,frq(ifreq)*2.*pie*del_time))
c     write(6,*) cexp(CMPLX(0.,frq(ifreq)*2.*pie*del_time))
               ENDDO            ! ranges
            ENDDO               ! depth
         ENDDO                  ! freq
      ENDIF                     ! del_time loop
     
      
      
      lwascomp=1
c     IF ((i_geom.EQ.0).AND.(i_call_porter.EQ.1)) THEN 
      IF ((i_geom.EQ.0)) THEN 
         i_call_porter=0
c     CLOSE(40)
c     OPEN(UNIT=40,STATUS='OLD',readonly,FORM='UNFORMATTED')
      ENDIF
      
c     c: WRITE transfer FUNCTION in a file:
c      write(*,*)'for snapinitit iifft,iWriteTrf=',iifft,iWriteTrf
      IF ((iifft .NE. 0).AND.(iWriteTrf.GT.0)) THEN
         WRITE(6,*)'Calling bb_fft_out for best model'
         CALL writetrf()
      ENDIF
      
      IF ((flagpu.LE.20)) THEN
c     WRITE(prtfil,*)' Exiting snap...'
      ENDIF 
      RETURN
c---- EXIT WITH error
 4000 lwascomp=-1
      RETURN
      END       
      
      SUBROUTINE mumble
      USE global
      INCLUDE 'comopt.h'
      WRITE(*,*)' the pressure of the last point is more than'
      WRITE(*,*)'  400dB; '
      WRITE(*,*)' This is probably due to a too low sound '
      WRITE(*,*)'speed in the'
      WRITE(*,*)' sediment or bottom.'
      WRITE(*,*)' a dump follows  in the *.out file'
      WRITE(prtfil,*)' the pressure of the last point is' 
      WRITE(prtfil,*)' less than -400dB; '
      WRITE(prtfil,*)' This is probably due to a too low sound '
      WRITE(prtfil,*)' speed in the'
      WRITE(prtfil,*)' sediment or bottom.'
      WRITE(prtfil,*)' a dump follows: '
      END
      
      SUBROUTINE writemodel
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
c     INCLUDE './snap/common.f'
      INTEGER i,jj,ifreq
      INCLUDE 'comsnap.h'
      
      WRITE(prtfil,*)' Frequency:',(frq(ifreq),ifreq=1,nfrq) 
      WRITE(prtfil,*)' R1,R2',R1,r2
      WRITE(prtfil,*)' SCATT(1),SCATT(2)',SCATT(1),SCATT(2)
      WRITE(prtfil,*)' H0,H1',H0,H1
      WRITE(prtfil,*)' ND0,ND1',ND0,ND1
      WRITE(prtfil,*)' BETA(1),BETA(2),BETA(3)',
     &     BETA(1),BETA(2),BETA(3)
      WRITE(prtfil,*)' c2,c2s' ,c2,c2s
      WRITE(prtfil,'(a,/,100(i3,2F10.3/))')
     &     ' z0,c0',(i,z0(i),c0(I),i=1,nd0)
      IF (ND1.GT.0) THEN
         WRITE(prtfil,'(a,/,100(i3,2F10.3/))')
     &        ' z1,c1',(i,z1(i),c1(I),i=1,nd1)
      ENDIF
      IF (tilt) WRITE(prtfil,*)' Tilt (m)=',dtilt
c AARON
      IF (offset) then
         WRITE(prtfil,*)'range offsets: ',(dr_ind(i),i=1,ndep)
      ENDIF
      IF (timeoffset) then
         WRITE(prtfil,*)'range offsets: ',(time_ind(i),i=1,ndep)
      ENDIF
      WRITE(prtfil,*)' np',np     
      WRITE(prtfil,*)' first reciever depth',rdep(1)
      WRITE(prtfil,*)' ranges',(rng(i),i=1,MIN(np,20))
      WRITE(prtfil,*)'number of sources',nsrd
      DO i=1,nsrd
         WRITE(prtfil,*)' source depth,i',srd(i),i  
      ENDDO
      WRITE(prtfil,*)' ndep,rd,rdlow: ',ndep,rd,rdlow        
      WRITE(prtfil,*)'arrayshape: ',(arrayshape(jj),jj=1,5)
      WRITE(prtfil,*)' multistatic: ',multistatic
      WRITE(prtfil,*)' tilt: ',tilt
      WRITE(prtfil,*)' mbeam:   ',mbeam
      WRITE(prtfil,*)' msource: ',msource
      WRITE(prtfil,*)' iscale:  ',  iscale
      WRITE(prtfil,*)' xplane:  ',  xplane
      
      CALL flush(prtfil)
      END

      
      SUBROUTINE writesnap
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
c     INCLUDE './snap/common.f'
      INTEGER i,jj,ifreq
      INCLUDE 'comsnap.h'
       CHARACTER*40 dumch
c      CHARACTER  optread(mopt)
      EQUIVALENCE (dumch,optinv)
     
      open(unit=40,file='snap.dat',form='formatted')
      write(40,'(a60)')title
      write(40,'(60a)')optinv,'!options'
      write(40,*)niter,qin,npopin
      write(40,*)px,pu,pm
      write(40,'(60a)')opt,' !snap options'
      write(40,*)nfrq,maxnomode
      WRITE(40,*)(frq(ifreq),ifreq=1,nfrq),' !frequency' 
      WRITE(40,*)H0,SCATT(1),SCATT(2),beta(0),' ocean ssp'
      do i=1,nd0 
        WRITE(40,*)z0(i),c0(I)
      enddo

      WRITE(40,*)H1,r1,beta(1),' sediment'
      do i=1,nd1 
        WRITE(40,*)z1(i),c1(I)
      enddo
      WRITE(40,*)r2,beta(2),c2,' bottom'
      WRITE(40,*)beta(3),c2s    
      WRITE(40,*)
      WRITE(40,*)srd(1), ' Source depth'
      WRITE(40,*)
      IF (tilt) then
c        WRITE(40,*)rdep(1),rdep(np),ndep, ' receivers'
         WRITE(40,*)rdep(1),rdep(ndep),irflagg,dtilt, ' receivers+tilt'
      else
c          WRITE(40,*)rdep(1),rdep(np),ndep,dtilt, ' receivers+tilt'
          WRITE(40,*)rdep(1),rdep(ndep),irflagg,dtilt,
     1 ' receivers '
      endif
c AARON added
      if (offset) then
	  write(40,*)(dr_ind(jj),jj=1,ndep)
      endif
      if (timeoffset) then
	  write(40,*)(time_ind(jj),jj=1,ndep)
      endif
      WRITE(40,*)
      if (irflagg.lt.0) then
	write(40,*)(rdep(jj),jj=1,ndep)
      endif
      IF  (out_plane) THEN 
         IF (xplane.GE.1) THEN
            DO jj=1,ndep 
               write(40,*)xhor(jj),yhor(jj)
            ENDDO
         ENDIF                  !xplane=1
         write(40,*)(arrayshape(jj),jj=1,2),' !array length and bow'
         write(40,*)(arrayshape(jj),jj=3,5),' !array rotations'
      ENDIF               ! out_plane

      WRITE(40,*)np,' number of ranges'
      WRITE(40,*)(rng(i),i=1,np)
      WRITE(40,*)
      
      WRITE(40,*)nparm
      if (iopt(9).eq.1.or.iopt(23).eq.1) then
         do i=1,nparm
            WRITE(40,99)par2phy(i),par2lay(i),fmin(i),fmax(i),
     1           ndigit(i),iapri(i)
         enddo
 99      format(i3,i3,g15.3,g15.3,i5,g15.3)
      else
         do i=1,nparm
            WRITE(40,99)par2phy(i),par2lay(i),fmin(i),fmax(i),
     1        ndigit(i)
         enddo
      endif
      WRITE(prtfil,*)' ndep,rd,rdlow: ',ndep,rd,rdlow        
      WRITE(prtfil,*)' arrayshape: ',(arrayshape(jj),jj=1,5)
      WRITE(prtfil,*)' multistatic: ',multistatic
      WRITE(prtfil,*)' tilt: ',tilt
      WRITE(prtfil,*)' mbeam:   ',mbeam
      WRITE(prtfil,*)' msource: ',msource
      WRITE(prtfil,*)' iscale:  ',  iscale
      WRITE(prtfil,*)' xplane:  ',  xplane
      close(40)
      CALL flush(prtfil)
      END

c*********************************************
      SUBROUTINE writeTL()
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
      if (fsbb==0) then
         dt1=1
      else
         dt1=1.E0/fsbb
      endif
      freqs=(fmaxdum-fmindum)/2.+fmindum
c     sd=zsr(mzsrc(1))
      if (iforwt==1) then
      write(*,*)'Calling trfhead'
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     rng(1),rspace,nfftbb,nf1,mx1,dt1,freqs,sd,
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
      if (iforwt==1) then
      write(*,*)'Calling trfhead'
      CALL trfhead(outroot,title,rdep(1),rdep(ndep),
     &     rng(1),rspace,nfftbb,nf1,mx1,dt1,freqs,sd,
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



