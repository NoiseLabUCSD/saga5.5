      SUBROUTINE forwardmodel(iopt,mopt)
      INTEGER  mopt,i,iopt(mopt)
      DO i=1,mopt
        iopt(i)=0
      ENDDO
      iopt(1)=2
      iopt(12)=1            ! using three indexes for adressing variable
      iopt(30)=7
      END 
c
       SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization PROGRAM
c     PETER GERSTOFT, 1992
c
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './cprosim/Parms_com'
      INCLUDE 'gen_i_o_saga.h'
      INCLUDE './cprosim/sector_env_com'
      INTEGER i,j,jj,ndum,ierr,irflag,nsctor
      REAL rdstep,rderr
      CHARACTER*80 dumch
      EQUIVALENCE (dumch,opt)
      REAL zreclo,zrechi,zsrclo,zsrchi,
     .     zwatlo,zwathi,zhtiltlo,zhtilthi,dzrec,dzsrc,
     .     dzwat,dzhtilt,dzlimit
      COMMON /parlimits/ zreclo,zrechi,zsrclo,zsrchi,
     .                   zwatlo,zwathi,zhtiltlo,zhtilthi,
     .                   dzrec,dzsrc,
     .                   dzwat,dzhtilt,dzlimit

      ierr=0


c---  READ the GA-parameters
      CALL readinputstart
      
      READ(1,'(a)') dumch         ! This line has no effect in snap
      CALL cprosimoption(opt)
      
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
         nfcw=nfrq
         nfbb=nfcw
         READ(1,*)(frq(I),i=1,nfrq)
         WRITE(*,'(a,10f10.3)')' frequencies:',(frq(i),i=1,nfrq)
         iicw=1
         DO i=1,nfrq
            fcw(i)=frq(i)
         ENDDO
         fmindum=fcw(1)
         fmaxdum=fcw(nfrq)
         fsbb=4*fmaxdum
      ELSE
         iicw=2
 10      WRITE(*,*)'reading frequencies...'
         READ(1,'(a80)')dumch
         WRITE(*,*)dumch
         READ(dumch,*,err=11)
     &    df_temp, fmindum, fmaxdum,dflast,fminlast,fmaxlast
         GOTO 12
 11    CONTINUE
c      BACKSPACE(1)
         READ(dumch,*,err=10)df_temp, fmindum, fmaxdum

 12      CONTINUE
         nf1=MAX(1,NINT(fmindum/df_temp)) + 1
         nf2=NINT(fmaxdum/df_temp) + 1
         fmindum=(nf1 - 1)*df_temp
         fmaxdum=(nf2 - 1)*df_temp
c--- pln inorder to SELECT a specific start point
c         nf1=4
c         nf2=9
c         fmindum=225.0
c         fmaxdum=725.0
         nfbb=nf2 - nf1 + 1
         nfrq=nfbb
         nfftbb=2**(INT(LOG10(float(nf2))/LOG10(2.) + .9))*2
         fsbb= dfloat(nfftbb)*df_temp
         WRITE(*,*)' Fmin, Fmax, DF, nfrq: ', 
     .        Fmindum, Fmaxdum, DF_temp,nfrq,nf2,nfftbb
         IF (nfrq.GT.mfreq) THEN
           WRITE(*,*)' mfreq,nfrq=',mfreq,nfrq
           STOP 'increase mfreq! '
         ENDIF 
         
         DO j=1,nfrq
            faxbb(j)=fmindum + (j-1)*(df_temp)
            frq(j)=faxbb(j)
            WRITE(99,*)j,frq(j)
         ENDDO 
      ENDIF
c     BACKSPACE(1)
c     BACKSPACE(1)

      
      CALL svp_read_pro
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
      IF (ndep.GT.mdep) THEN
         WRITE(*,*) 'ndep must grater than mdep',ndep,mdep
         ierr=1
      ENDIF
c     
c---  READ the range-BLOCK
      CALL  read_range_inp(ierr,ndum)
c---  transfer receiver parameters
      DO i=1,nx
         rkm(i)=xranges(i)
      ENDDO
      rmin=rkm(1)/1000.e0
      rmin=1.e-20
      nsrc=nx
c---  source info
      nzs=1
      zsrc(1)=sd
c---  passing receiver depth
c
      irflag=ndep
      ndep=ABS(ndep)
      nrec=ndep
      IF (ndep.GT.1) THEN
         rdstep=(rdlow-rd)/float(ndep-1)
      ELSE
         IF(tilth) THEN
            nrec = nsrc
            rdstep = dtilth/float(nrec-1)
         ELSE
            rdstep=1.
         END IF
      ENDIF

      IF (irflag.GT.0) THEN
c     equidistant receiver depths
         DO  jj=1,ndep
            rdep(jj)=(jj-1)*rdstep+rd
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
      DO i=1,nrec
         zrec(i)=(i-1)*rdstep+rd
      ENDDO 
c     
cpln      nrec=ndep
      WRITE(*,*)' Number of depth', ndep
      WRITE(*,*)' Number of ranges,nx,nrec,nsrc', nx,nrec,nsrc
      WRITE(*,*)' first & last depth',rdep(1),rdep(ndep)
      WRITE(*,*)' first & last depth',zrec(1),zrec(nrec)
c---  range parameters
      
c     nrng=nx                        ! number of range steps
c     DO i=1,nx
c     xsrc(i)=xranges(i)
c     ENDDO
      WRITE(*,*)' Number of ranges', nx
      WRITE(*,*)' first & last range',xranges(1),xranges(nx)
      
         
c---  should the  EOF be READ ? 
      
      IF (iopt(17).EQ.1) THEN
         CALL eofinit
      ENDIF
      
c***  create text strings for physical parameters
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
      phystxt(19)= 'Vertical array tilt (m)$'
      phystxt(20)= 'Time delay $'
      phystxt(21)= 'Horizontal array tilt (m)$'
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
         IF (par2phy(i).LT.1 .OR. par2phy(i).GT.21 .OR.
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
      DO i=1,nparm
         IF ((par2phy(i).EQ. 1).OR.
     1        (par2phy(i).EQ. 8).OR.
     1        (par2phy(i).EQ. 9).OR.
     1        (par2phy(i).EQ.15).OR.
     1        (par2phy(i).EQ.19).OR.
     1        (par2phy(i).EQ.20).OR.
     1        (par2phy(i).EQ.21)) THEN
            i_geom=0            ! we DONT have to compute the modes
         ELSE
            i_geom=1            ! we have to compute the modes
            GOTO 50
         ENDIF
      ENDDO
 50   CONTINUE

      i_geom=1  
      WRITE(*,*)
      WRITE(*,*)' Should modes be calculated (1=Yes,0=NO)? ',i_geom
      WRITE(*,*)

      
c
c     Put limits of receiver depths in zreclo and zrechi to
c     be used IF ONLY geometry search WITH CPROSIM
        IF(par2phy(i).EQ.1) THEN
           zwatlo=fmin(i)
           zwathi=fmax(i)
        END IF
        IF(par2phy(i).EQ.8) THEN
           zsrclo=fmin(i)
           zsrchi=fmax(i)
        END IF
        IF(par2phy(i).EQ.15) THEN
           zreclo=fmin(i)
           zrechi=fmax(i)
        END IF
        IF(par2phy(i).EQ.21) THEN
           zhtiltlo=fmin(i)
           zhtilthi=fmax(i)
        END IF



c*** errors ?
      IF (ierr.EQ.1)STOP 'stopped in input'
c***  CLOSE input unit
      CLOSE(1)                                
      END
c**************************************************************
      SUBROUTINE cprosimoption(options)
      USE global
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
c      INCLUDE 'comsnaprd.h'
c      INCLUDE 'comprosim.h'

      CHARACTER options(40)
      INTEGER i
      INTEGER ipos,iposV
      COMMON /div/ipos,iposV
c
      COMMON /tiltparm/tiltv,tilth,dtiltv,dtilth
      LOGICAL tiltv,tilth
      REAL dtiltv,dtilth
      ipos=0
      tiltv=.FALSE.
      tilth=.FALSE.
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
            WRITE(*,*)' tilt of vertical array is included'
            WRITE(prtfil,*)' tilt of vertical array is included'
         ELSE IF (options(i).EQ.'T') THEN
            tilth=.TRUE.
            WRITE(*,*)' tilt of horizontal array is included'
            WRITE(prtfil,*)' tilt of horizontal array is included'
         ELSE IF (options(I).EQ.'!') THEN
            GOTO 60
         ELSE IF (options(I).NE.' ') THEN
            WRITE(prtfil,399) options(I)
 399        FORMAT(1H ,' >>>> unknown prosim OPTION: ',A1,' <<<<')
         END IF
 50   CONTINUE
 60   CONTINUE
      IF((tiltv).AND.(tilth)) THEN
         WRITE(6,*)'Horizontal and vertical tilt together ',
     .             ' is not allowed'
         STOP
      END IF
      END
c******************************************************************
       SUBROUTINE forwinit
c
c: Normal mode model for acoustic propagation in range dependent ocean
c: environments WITH fluid and/or solid layered structures above and below.  
c
c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: References for the code as of 3/22/95 are:
c:    E. K. Westwood, C. T. Tindle, and N. R. Chapman, "A normal mode
c:       model for multilayered acoustoelastic ocean environments based on
c:       an analytic reflection coefficient method," J. Acoust. Soc. Am.,
c:       95, No. 5, Pt. 2, 2908 (1994).
c:    E. K. Westwood, "An efficient broadband normal-mode model for
c:       acoustoelastic ocean environments," J. Acoust. Soc. Am., 96,
c:       No. 5, Pt. 2, 3352 (1994).
c
c: *******************************
c: *     REVISION (1996):        *
c: *     P R O S I M  GROUP      *
c: *     S A C L A N T C E N     *
c: *******************************

      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE './cprosim/Parms_com'
      INCLUDE 'gen_i_o_saga.h'
      INCLUDE './cprosim/depth_com'
      INCLUDE './cprosim/sector_env_com'

      INTEGER*4 isectr,intrp(150),isubb,nsctor
      REAL*4  prsus(1),dt2,tmax,suswgh,rngval,fcntr
      REAL*4 rdstep, zrectemp(NSRMAX)
      INTEGER ioer, izero, iilay
      INTEGER noff,ndt
      INTEGER*4 nrectemp,nxtemp
      COMMON /info/ienv_info,lu_env
      COMMON/iterpar/iter,iforwpop
      INTEGER ienv_info,lu_env,iter,iforwpop
      INTEGER ipos,j
      COMMON /div/ipos
      LOGICAL dosus,dohilb
      EXTERNAL time
      INTEGER*4 time,isec1,kdat(3),ktim(3)

      REAL zreclo,zrechi,zsrclo,zsrchi,
     .     zwatlo,zwathi,zhtiltlo,zhtilthi,dzrec,dzsrc,
     .     dzwat,dzhtilt,dzlimit
      COMMON /parlimits/ zreclo,zrechi,zsrclo,zsrchi,
     .                   zwatlo,zwathi,zhtiltlo,zhtilthi,
     .                   dzrec,dzsrc,
     .                   dzwat,dzhtilt,dzlimit



c
      iiwrite=3

      IF(iiwrite.GT.0)
     .  PRINT *,'Running COMPLEX PLANE PROSIM ...'

c      IF (Mseg .NE. Nsegmax) THEN 
c         WRITE(*,*)' Mseg in file comsnaprd is not equal to  Nsegmax'
c         WRITE(*,*)' in prosim/Params_com. please set them equal.'
c         STOP
c      ENDIF
      IF (iopt(6).EQ.1) THEN
         iiwrite=10
         ienv_info=1         ! write the environment for each fm to file 90
         lu_env=90
         CALL opfilw(lu_env,ioer)
      ENDIF
      i_call_porter=1
      tot_jjfail=0
C***********************************************
c
c
      ENTRY forw2
c
c---- IF writing out trf FUNCTION
c
        IF (iWriteTrf.EQ.1 .AND. 
     &     dflast.GT.0 .AND. fmaxlast.GT.fminlast) THEN
         df_temp=dflast
         fmindum=fminlast
         fmaxdum=fmaxlast
         nf1=MAX(1,NINT(fmindum/df_temp)) + 1
         nf2=NINT(fmaxdum/df_temp) + 1
         fmindum=(nf1 - 1)*df_temp
         fmaxdum=(nf2 - 1)*df_temp
         nfbb=nf2 - nf1 + 1
         nfrq=nfbb
         nfftbb=2**(INT(LOG10(float(nf2))/LOG10(2.) + .9))*2
         fsbb= dfloat(nfftbb)*df_temp
         WRITE(*,*)' For optimized env Fmin, Fmax, DF, nfrq: ', 
     .        Fmindum, Fmaxdum, DF_temp,nfrq,nf2,nfftbb
         IF (nfrq.GT.mfreq) STOP 'increase mfreq! '
         
         DO j=1,nfrq
            faxbb(j)=fmindum + (j-1)*df_temp
            frq(j)=faxbb(j)
         ENDDO           
         i_call_porter=1
         i_geom=1
      ENDIF

c
c SAVE receiver depths in temp array until modes and eigenvalues 
c have been calculated
      IF((i_geom.EQ.0) .AND. (i_call_porter.EQ.1)) THEN
         nrectemp=nrec
         dzrec=zrechi-zreclo
         dzsrc=zsrchi-zsrclo
         dzwat=zwathi-zwatlo
         dzhtilt=zhtilthi-zhtiltlo
         dzlimit=(dzrec+dzsrc+dzwat+dzhtilt)/4.e0
         dzlimit=0.1e0*r_h0(1)
         DO j=1,nrec
            zrectemp(j)=zrec(j)
         END DO
         WRITE(*,*)'nrecsamp,nrec',nrecsamp,nrec

         rdstep=(MAX(zrec(nrec),zsrc(1))-MIN(zrec(1),zsrc(1))
     .          +2.e0*dzlimit)/float(NRECSAMP-1)
         nrec=NRECSAMP
      WRITE(*,*)'PETER ISECTR LOOP, nrec,nsrd',nrec,nsrc
         DO j=1,NRECSAMP
            zrec(j)=(j-1)*rdstep+zrectemp(1)-dzlimit
         END DO
         nxtemp=nx
         nx=1
         nsrc=1
         iwat0=r_h0(1)
      END IF
c
      CALL clean_array
      CALL vbl_init
      ver_cur=2.0
      IF(iiwrite.GT.0)
     .   PRINT *,'Running CPROSIM v2.0b ...'
cpln      CALL flush(6)

c      isec1=time()
c
c: For SGI:
cxx   CALL idate(kdat(2),kdat(1),kdat(3))
c: For Sun:
c      CALL idate(kdat(1),kdat(2),kdat(3))
c      kdat(3)=MOD(kdat(3),100)
c      CALL itime(ktim)
c
      IF(ienv_info .NE. 1)
     &      iiwrite= iiwrite-1
      lwascomp=1
c
c---- USE of EOF
c
      IF (iopt(17).EQ.1) THEN
         CALL eofval
         IF (lwascomp.EQ.-1) THEN
            RETURN
         ENDIF
      ENDIF
      nsctor=1
c      PRINT *, 'nap time'
c
      IF (ipos.EQ.1) THEN
         DO isectr=1,nsctor
c the impedance contrast should be positive
            IF(R_H1(isectr) .NE. 0.0) THEN
               IF (r_r1(isectr).GT.r_r2(isectr)) THEN
                  lwascomp=-1
                  WRITE(*,*)' Model',iforwpop,' rejected. dens-contr'
                  RETURN
               ELSEIF (R_c1(R_ND1(ISECTr),isectr).GT.
     .                (R_c2(isectr))) THEN
                  lwascomp=-1
                  WRITE(*,*)' Model',iforwpop,' rejected. Csed > Cbot'
                  RETURN
               ENDIF
c     for sediment profiles WITH negative slope of ssp
               DO iilay=2,R_ND1(ISECTR)
                  IF(R_c1(iilay-1,isectr).GT. 
     .                 R_c1(iilay,isectr)) THEN
                     lwascomp=-1
                     WRITE(*,*)' Model',iforwpop,' rejected. neg slope'
                     WRITE(*,*)' between layer',iilay-1,' and ',iilay
                     RETURN
                  END IF
               END DO
            ENDIF
         ENDDO
      ENDIF
      IF (ipos.EQ.2) THEN
         DO isectr=1,nsctor
            IF(R_c2(isectr).LT.1530)THEN
               WRITE(*,*)' Model',iforwpop,' rejected. Cbot<1530'
               lwascomp=-1
               RETURN
            ENDIF
         ENDDO
      ENDIF
      IF(zsrc(1).GE.r_h0(1)) THEN
         WRITE(*,*)' Model',iforwpop,'rejected. Source in the bottom'
         lwascomp=-1
         RETURN
      END IF
c
      iifail=0
      jjfail=0
      isubb = 1
c: Set iiwrite=1 so that output files and messages are sent:
c      
c      PRINT *, 'before opt_read'
      CALL opt_read_pro
c      PRINT *, 'after opt_read'
c     
      DO isectr=1,nsctor
         R_z0(R_nd0(isectr),isectr)=R_h0(isectr) ! last point is the water depth
         IF (R_nd1(isectr).GT.0) THEN
            R_z1(R_nd1(isectr),isectr)=R_h1(isectr)
c            WRITE(6,*)'R_z1,R_h1,R_nd1: ',R_z1(R_nd1(isectr),isectr),
c     .           R_h1(isectr),R_nd1(isectr)
c            PAUSE
         END IF
      ENDDO
      
      iirx=0

      DO isectr=1,nsctor 
c         WRITE(*,*)'sector',isectr,iiwrite
c         WRITE(*,*)'nap time'
         CALL tenv_read_pro(isectr,intrp,isubb)
         CALL svp_check(1)
         CALL svp_check2
         IF(jjfail.GT.1)
     .      GOTO 1001
         IF (iiwrite.GE.2) THEN
            PRINT *, '# OF SECTOR: ', isectr
            PRINT *, 'LENGHT OF SECT: ', secleng(isectr)
         ENDIF
c
         IF(iicw .EQ. 1) THEN
c     : CW mode computations:
c pg           CALL cw_modes(0)
               CALL bb_brute
          ELSEIF(iicw .EQ. 2) THEN
c     : BB mode computations:
            IF(iirx .EQ. 0) THEN
               CALL bb_brute
            ELSEIF(iirx.EQ.-1) THEN
               CALL bb_modes
            END IF
         ENDIF
         
         isubb=isubb+1          !subband counter
         
      ENDDO                     !FINISH SUBBAND CALCULATIONS
c
      IF ((i_geom.EQ.0) .AND. (i_call_porter.EQ.1)) THEN 
         i_call_porter=0
         nrec=nrectemp
         DO j=1,nrec
            zrec(j)=zrectemp(j)
         END DO
         nx=nxtemp
         nsrc=nx
         WRITE(6,*)'nrec,nsrc: ',nrec,nsrc
c        CLOSE(40)
c        OPEN(UNIT=40,STATUS='OLD',readonly,FORM='UNFORMATTED')
      ENDIF

      IF ((iWriteTrf.GT.0)) THEN
         WRITE(6,*)'Calling bb_fft_out for best model'
         CALL bb_fft_out()
      ENDIF
      CALL bb_fft_out()
c     PAUSE
c     
c     DO Hilbert transform (remove carrier
c     dohilb=.TRUE.
      dohilb=.FALSE.
c
      fcntr=(fmaxdum+fmindum)/2.
c     Generate source replica of SUS charge
c      dosus=.TRUE.
      dosus=.FALSE.
      IF(dosus) THEN
         tmax=1/df_temp
c         dt=1/(8*fmaxdum)
c         IF(tmax .LT. .5) THEN
c            tmax=2/df_temp
c         ENDIF
c         ndt=INT(tmax/dt)
c         npow2=ndt
c         DO j=1,10000
c            npow2=2**j
c            IF(npow2 .GT. ndt) GOTO 111
c         END DO
c 111     ndt=npow2
         ndt=32768
         dt2=tmax/(ndt-1)
c
         noff=0
         suswgh=820.e0
         rngval=1.e0

         IF(iiwrite .GT. 0) THEN
            WRITE(6,*)
            WRITE(6,*)'======================================='
            WRITE(6,*)'Source depth      = ',zsrc(1)
            WRITE(6,*)'Source weight     = ',suswgh
            WRITE(6,*)'Source range      = ',rngval
            WRITE(6,*)'Frequency sampling= ',df_temp
            WRITE(6,*)'Min frequency     = ',fmindum
            WRITE(6,*)'Max frequency     = ',fmaxdum
            WRITE(6,*)'No frequency      = ',nfbb
            WRITE(6,*)'First frequency   = ',nf1
            WRITE(6,*)'Last frequency    = ',nf2
            WRITE(6,*)'Time sampling     = ',dt2
            WRITE(6,*)'No time sample    = ',ndt
            WRITE(6,*)'Total time window = ',dt2*(ndt-1)
            WRITE(6,*)'======================================='
            WRITE(6,*)
         ENDIF
c
         CALL susequ(noff,prsus,dt2,ndt,suswgh,zsrc(1),rngval)
c
c         DO j=1,ndt
c            WRITE(78,*)(j-1)*dt,prsus(2*j-1)
c         END DO
c
         CALL four1(prsus,ndt,1)
c
c         DO j=1,ndt/2
c            WRITE(79,*)(j-1)*df_temp,prsus(2*j-1)*dt,
c     .           prsus(2*j)*dt
c         END DO
c         CLOSE(78)
c         CLOSE(79)
c         STOP
      END IF
c
      CALL saga_out(nsrc,nrec,nfbb,prsus,dosus,dohilb,nf1,
     .     fcntr)
      
      isubb=1
      
c     c: Check IF failure occurred in mode_find or eig_find:

cpln      WRITE(6,*)'JJFAIL returned from CPROSIM: ',jjfail
 1001 IF(jjfail .GT. 0) THEN
         tot_jjfail=tot_jjfail+1
         lwascomp=-1
         WRITE(*,*)
cpln         WRITE(*,*)'FAILURE OCCURRED FOR THIS CPROSIM RUN.'
         WRITE(*,*)' Model',iforwpop,' rejected.'
         WRITE(*,*)' Failure for this CPROSIM run'
         WRITE(*,*)
         WRITE(21,*)' Model',iforwpop,' rejected.'
         WRITE(21,*)' Failure for this CPROSIM run'
         WRITE(21,*)' Number of failures for this inversion: ',
     .                tot_jjfail
         WRITE(21,*)
c         STOP
         RETURN
      ENDIF
      RETURN
      
      END
      
      SUBROUTINE saga_out(nsrc,nrec,nfbb,prsus,dosus,dohilb,
     .           nf1,fcntr)
c
c: *******************************
c: *     REVISION (1996):  
c************************************
c Developed exclusively for the USE by SAGA 
c      *
c: *     P R O S I M  GROUP      *
c: *     S A C L A N T C E N     *
c: *******************************
      USE global

      IMPLICIT NONE
      INCLUDE './cprosim/Parms_com'
c     INCLUDE 'i_o_com'
c pg      INCLUDE 'i_o_1a_com'
c pg      INCLUDE 'i_o_1b_com'
c pln      INCLUDE './prosim/i_o_1tf_com'
c pg      INCLUDE 'i_o_2_com'
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      COMMON /out_com1/ tf(NTFMAX)
      INTEGER   nsrc,nrec,nfbb
      INTEGER nf1,nrec_htilt
      INTEGER jj0,jsrc,jrec,j,index
      REAL prsus(65536),fcntr,pie
      COMPLEX*8 tf
      LOGICAL dosus,dohilb
      COMMON /tiltparm/tiltv,tilth,dtiltv,dtilth
      LOGICAL tiltv,tilth
      REAL dtiltv,dtilth
c
c      DO jj=nf1,nfbb+nf1-1
c         WRITE(80,*)(jj-1)*df_temp,SQRT(prsus(2*jj-1)**2
c     .           +prsus(2*jj)**2)
c      END DO
c      CLOSE(80)
c      STOP
c
      nrec_htilt=1
      pie=4.0*ATAN(1.0)
      IF(tilth) THEN
         DO j=1,nfbb            ! frequency
            DO jsrc=1,nsrc      !  range
               jj0=(jsrc-1)*nfbb*nrec+(jsrc-1)*nfbb
               DO jrec=1,nrec_htilt   !depth
                  index=(jrec +((j-1))*nrec_htilt-1)*nsrc
                  IF(dohilb) THEN
                     pie=4.*ATAN(-1.)
                     resp(jsrc+index)= CONJG(2.*tf(j+jj0)*
     .                    cexp(-CMPLX(0.,2.*pie*fcntr)))
                     
                  ENDIF
                  IF(dosus) THEN
                     resp(jsrc+index)= -(CONJG(tf(j+jj0)*
     .                    CMPLX(prsus(2*(j-1+nf1)-1),
     .                    prsus(2*(j+nf1-1)))))
                  ELSE
                     resp(jsrc+index)= -(CONJG(tf(j+jj0)))
c                     WRITE(6,*)resp(jsrc+index)
c                     WRITE(6,*)-(CONJG(tf(j+jj0)))
                  END IF
                  IF (iopt(5).EQ.5) THEN
                     resp(jsrc+index)=resp(jsrc+index)
     .                    *cexp(CMPLX(0.,frq(j)*2.*pie*del_time))
                  ENDIF
               END DO
c               PAUSE
            END DO
         END DO
      ELSE
         DO j=1,nfbb            ! frequency
            DO jsrc=1,nsrc      !  range
               jj0=(jsrc-1)*nfbb*nrec
               DO jrec=1,nrec   !depth
                  index=(jrec +((j-1))*nrec-1)*nsrc
                  IF(dohilb) THEN
                     pie=4.*ATAN(-1.)
c                  WRITE(6,*)'======================'
c                  WRITE(6,*)'pie,fcntr,dohilb',pie,fcntr,dohilb
c                  WRITE(6,*)'======================'
                     resp(jsrc+index)= CONJG(2.*tf(j+jj0)*
     .                    cexp(-CMPLX(0.,2.*pie*fcntr)))
                     
                  ENDIF
                  IF(dosus) THEN
                     resp(jsrc+index)= -(CONJG(tf(j+jj0)*
     .                    CMPLX(prsus(2*(j-1+nf1)-1),
     .                    prsus(2*(j+nf1-1)))))
c     WRITE(77,*)frq(j),-REAL(resp(jsrc+index)),
c     .                       -imag(resp(jsrc+index))
                  ELSE
                     resp(jsrc+index)= -(CONJG(tf(j+jj0)))
c     WRITE(77,*)frq(j),REAL(resp(jsrc+index)),
c     .                       imag(resp(jsrc+index))
                  END IF
c     WRITE(*,*)'resp,tf',index,resp(jsrc+index),tf(j+jj0)
c     tf(j+jj0)=zzero
                  jj0=jj0 + nfbb
                  IF (iopt(5).EQ.5) THEN
                     resp(jsrc+index)=resp(jsrc+index)
     .                    *cexp(CMPLX(0.,frq(j)*2.*pie*del_time))
                  ENDIF
               END DO
            END DO
         END DO
      END IF
c
c      CLOSE(77)
c     
      lwascomp=1
      RETURN
      END
c
      SUBROUTINE SUSEQU(NOFF, PRES, DT2, NDT, W, D, R)
C*----******************************************************
C*----* THIS SUBROUTINE FILL PRES WITH THE SOURCE WAVEFORM *
C*----* FOR A SUS CHARGE.                                  *
C*----******************************************************
c
      INTEGER NDT,K,J,NOFF
      REAL PIE,WROOT,WRRAT,PRES,W,D,R,TST,DT2,T,TB,P,PT
      REAL E0,E1,PP
      DIMENSION P(0:4), PP(4), PRES(65536), T(0:4), TB(0:4)

c      WRITE(6,*)'inside susequ'
c      WRITE(6,*)'noff= ',noff
c      WRITE(6,*)'dt  = ',dt
c      WRITE(6,*)'ndt=  ',ndt
c      WRITE(6,*)'w=    ',w
c      WRITE(6,*)'d=    ',d
c      WRITE(6,*)'r= ',r

      PIE=3.141592654
      WROOT=W**(1./3.)
      WRRAT=WROOT/R

      P(0)=3.74E12*(WRRAT**1.13)
      P(1)=9.02E11*WRRAT
      P(2)=.22*P(1)
      P(3)=.10*P(1)
      P(4)=.03*P(1)

      T(0)=1.75E-5*WROOT*(WRRAT**(-.26))
      T(1)=1.48E-4*WROOT/((D+10.1)**(1./6.))
      T(2)=1.91*T(1)
      T(3)=2.79*T(1)
      T(4)=2.79*T(1)

      TB(0)=0.
      TB(1)=.21*WROOT/((D+10.1)**(5./6.))
      TB(2)=1.71*TB(1)
      TB(3)=2.28*TB(1)
      TB(4)=2.81*TB(1)

      DO 5 J=1,4
         PP(J)=(PIE/2.)*(P(J-1)*T(J-1) + P(J)*T(J))
     &    /(TB(J-1)-TB(J))
5     CONTINUE

      DO 10 K=1,NDT
         TST=FLOAT(K-1)*DT2
         PT=0.
         E0=-1*ABS(TST/T(0))
         IF((E0 .GT. -656.) .AND. (TST .GE. 0.)) PT=PT+P(0)*EXP(E0)
         DO 20 J=1,4
           E1=-1*ABS((TST-TB(J))/T(J))
           IF(E1 .GT. -656.) PT=PT + P(J)*EXP(E1)
           IF((TST-TB(J-1))*(TST-TB(J)) .LE. 0.)
     &         PT=PT + PP(J)*SIN(PIE*(TST-TB(J-1))/(TB(J)-TB(J-1)))
20       CONTINUE
         PRES(2*K)=0.0
         PRES(2*K-1)=PT
10    CONTINUE

      NOFF=0

      RETURN
      END
C
      SUBROUTINE four1(DATA,nn,isign)
      INTEGER isign,nn
      REAL DATA(2*nn)
      INTEGER i,istep,j,m,mmax,n
      REAL tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
      DO 11 i=1,n,2
        IF(j.GT.i)THEN
          tempr=DATA(j)
          tempi=DATA(j+1)
          DATA(j)=DATA(i)
          DATA(j+1)=DATA(i+1)
          DATA(i)=tempr
          DATA(i+1)=tempi
        ENDIF
        m=n/2
1       IF ((m.GE.2).AND.(j.GT.m)) THEN
          j=j-m
          m=m/2
        GOTO 1
        ENDIF
        j=j+m
11    CONTINUE
      mmax=2
2     IF (n.GT.mmax) THEN
        istep=2*mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=-2.d0*SIN(0.5d0*theta)**2
        wpi=SIN(theta)
        wr=1.d0
        wi=0.d0
        DO 13 m=1,mmax,2
          DO 12 i=m,n,istep
            j=i+mmax
            tempr=sngl(wr)*DATA(j)-sngl(wi)*DATA(j+1)
            tempi=sngl(wr)*DATA(j+1)+sngl(wi)*DATA(j)
            DATA(j)=DATA(i)-tempr
            DATA(j+1)=DATA(i+1)-tempi
            DATA(i)=DATA(i)+tempr
            DATA(i+1)=DATA(i+1)+tempi
12        CONTINUE
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      CONTINUE
        mmax=istep
      GOTO 2
      ENDIF
      RETURN
      END
