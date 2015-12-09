      SUBROUTINE forwardmodel(iopt,mopt)
      INTEGER  mopt,i,iopt(mopt)
      DO i=1,mopt
         iopt(i)=0.
      ENDDO
      iopt(1)=4
      iopt(30)=5
      END
c******************************
      SUBROUTINE input
c     reads and interpretates the input file to the genetic algorithm 
c     optimization PROGRAM
c     PETER GERSTOFT, 1995

C     Input for POPP CMODES/OMODES: READ from DATA file 
C     DATE WRITTEN: 7-Oct-91
C     LAST EDIT:   20-Oct-94

C     MODIFICATIONS:
C     7-Oct-91: Edits to READ things from a file 
C     August 94: READ depth dependent attenuation 

      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'
      INTEGER i,ierr,j,ifrq
      REAL xhelp
      CHARACTER*80 dumch
      luinp=1

      IERR = 0
c---  READ the GA-parameters
      CALL readinputstart

C     Sound speed profile
      READ (LUINP,*) H, RHOWC, ATTFLG
      IF (ATTFLG.GT.100) THEN
         PRINT*,' TOO MANY ATTENUATION POINTS'
         STOP
      ENDIF

      IF ((1.0*INT(ATTFLG).NE.attflg).OR.(attflg.EQ.0.)) THEN 
         depthinc = 1.0
         alpwc(1) = attflg
         attflg = 1
      ELSEIF (ATTFLG.GT.0.) THEN
         READ(LUINP,*) depthinc, (alpwc(i),i=1,attflg) 
      ENDIF
      I = 0
 80   I = I + 1
      READ (LUINP,*) ZP(I), CP(I)
      WRITE (*,*) ZP(I), CP(I)
      IF (CP(I).GT.0.) GO TO 80
      NWL = I-1
      WRITE(*,*)' Number of water layers',nwl

C     Bottom properties
      READ (LUINP,*) C2S, ALP2S
      I = 0
 18   I = I + 1
      READ (LUINP,*) HBOT(I), RHOBOT(I), CBOT(I), ALPBOT(I) 
      WRITE    (*,*) HBOT(I), RHOBOT(I), CBOT(I), ALPBOT(I) 
      IF (HBOT(I).GT.0.) GO TO 18
      NBL = I

      READ (LUINP,*) ROUGHS, ROUGHB(1)

      READ (LUINP,*) NFRQ, (FRQ(ifrq),ifrq=1,nfrq)
      WRITE(*,*)' Number of frequencies',nfrq
      WRITE(*,*)' Frequencies:',(frq(ifrq),ifrq=1,nfrq)
      READ (LUINP,*) NS, (ZS(I), I= 1, NS)
      READ (LUINP,*) NR, (ZR(I), I= 1, NR)
      READ (LUINP,*) NL, NMAX
      IF (NL.LE.0) NL = 251
      IF (NMAX.LE.0) NMAX = MAXNM
      READ (luinp,*) isb
      READ (luinp,*) TAU0
      READ  (luinp,*) (sour_lev(ifrq),ifrq=1,nfrq)
      READ  (luinp,*) (nois_lev(ifrq),ifrq=1,nfrq)
      READ (luinp,*) (DMU(ifrq),ifrq=1,nfrq)
      READ (luinp,*) lampow
      WRITE(*,*)' tau0  =',tau0
      WRITE(*,*)' source strenght=',(sour_lev(ifrq),ifrq=1,nfrq)
      WRITE(*,*)' noise level=    ',(nois_lev(ifrq),ifrq=1,nfrq)
      WRITE(*,*)' dmu = ',(dmu(ifrq),ifrq=1,nfrq)
      WRITE(*,*)' lampow =',lampow
c     ri0= 10.**(xhelp/10.)
      DO ifrq=1,nfrq
         sour_lev(ifrq)=10.**(sour_lev(ifrq)/10)
         nois_lev(ifrq)=10.**(nois_lev(ifrq)/10)
      ENDDO
c     WRITE(*,*)'source level',xhelp,'dB, lin=',ri0
      READ (luinp,*) TMIN,TMAX,TINC
      nx=(tmax-tmin)/tinc +0.5
      xincre=tinc
      xminimum=tmin
      xmaxmum=tmax
      WRITE(*,*)'xminimum,xmaxmum,xincre'
      WRITE(*,*)xminimum,xmaxmum,xincre
      ncurv=nfrq*nr
      ndep=1
c     nfrq=1
c***  create text strings for physical parameters
c     123456789012345678901234567890
      phystxt(1) = 'Water depth (m) $'
      phystxt(2) = 'Water sound speed (m/s)$'
      phystxt(3) = 'Sediment sound speed (m/s)$'
      phystxt(4) = 'Sediment attenuation ()$'
      phystxt(5) = 'sediment thickness (m)'
      phystxt(6) = 'Sediment density (g/cm3)$'
      phystxt(8) = 'Source depth (m) $'
c     phystxt(9) = 'Source range (km) $'
      phystxt(11)= 'Shape coefficient $'
      phystxt(15)= 'Receiver depth (m)$'
      phystxt(16)= 'Depth of water speed pt.(m)$'
      phystxt(17)= 'Lambert coefficient$'
      phystxt(18)= 'Lambert power$'
c***  

      DO 8 j=1,mphys
         phystxt2(j)='                                               '
 6       DO 7 I=40,1,-1
            IF(phystxt(j)(I:I).NE.'$') GO TO 7
            phystxt2(j)(1:I-1)=phystxt(j)(1:I-1)
c     WRITE(*,*)phystxt2(j)
            GO TO 8
 7       CONTINUE
    8 CONTINUE
      
      IF (iopt(17).EQ.1) THEN
         CALL eofinit
      ENDIF

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
            IF (fmin(i).LT.zp(nwl-1)) THEN
               WRITE(*,*)' *** Optimzation variable #:',i
               WRITE(*,*)' *** Durring optimization'
               WRITE(*,*)' the water-depth can get above',
     &              ' the second-last'
               WRITE(*,*)' point in the water'
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

      ENDDO


c***  errors ?
      IF (ierr.EQ.1)STOP 'stopped in input'
c***  CLOSE input unit
      CLOSE(1)                                    
      
      END

      SUBROUTINE forwinit

C     Produces normal mode inputs for reverberation model OGOPOGO/REVERB;
C     also transmission loss. Basically calls the PROLOS/MODES subprogram 
C     OMODES and transmission loss SUBROUTINE PRPLOS. 

C     AUTHOR:       Dale D. Ellis
C     DATE WRITTEN: February 1991
C     LAST EDIT:    14-Dec-94

C     USAGE:
C     LINK  POPP, ORCMOD, nogrp
C     Input  file: POPP.INP  (See also files INPMOD.*) 
C     (See POPP.HLP for description of input parameters)
C     Output files:

C     MODIFICATION HISTORY:
C     May 1985:   PROGRAM MUNK
C     March 1988: PROGRAM GMODE
C     30-Jan-91:  Quick and dirty conversion to CALL OMODES and PRPLOS 
C     13-Sep-91:  OGOBIN added. Renamed to POPP, Inp/Out... formatting. 
C     7-Oct-91:  Added INPMOD to READ from a file 
C     11-Nov-92:  Edited initial comments
C     5-Jun-94:  Removed GMODES references
C     1993-94:  Converted to HP - Joe Farrell? 
C     OPEN statements and user supplied file names
C     More modes.
C     22-Aug-94:  J Theriault: Modifications to attenuation 
C     Modified INPMOD to READ attenuation array
C     Modified FUNCTION ATTDEP to interpolate input array
C     20-Oct-94:  Merge of recent Centre and HP modifications [DDE] 
C     The HP code (Version 4.1e) was the starting point
C     14-Dec-94:  Cosmetic edits and revised file POPP.HLP 


C     REFERENCES: (for the PROLOS/MODES package) 
C     D. D. Ellis, "A two-ended shooting technique for calculating 
C     normal modes in underwater acoustic propagation,"
C     DREA Report 85/105, September 1985.
C     T. J. Deveau, "Enhancements to the PROLOS normal-mode acoustic 
C     propagation model," Oceanroutes Canada Inc., DREA Contractor
C     Report CR/89/442, December 1989. Limited distribution.
C     B. A. Leverman, "User's guide to the normal mode propagation loss 
C     package PROLOS," DREA Research Note AM/82/3, June 1982. Informal
C     communication.
C     B. A. Leverman and D. D. Ellis, "Software documentation for the 
C     normal mode subprogram MODES," DREA Research Note AM/82/4,
C     June 1982. Informal communication.


C     COMMENTS and wish list:
C     This is a quick and dirty propagation loss model modified to 
C     produce a binary file for the reverberation model OGOPOGO. 
C     The calling SEQUENCE of OMODES should be modified: 
C     - C2S, ALP2S, IOPT, ROPT, ALPWC
C     - the amplitude and phase arrays are horrendously large 
C     There are lots of unused variables in the code 
C     The PRPLOS SUBROUTINE should give a scintillation index too 
C     This PROGRAM should eventually produce output for energy pulse 
C     propagation and reverberation calculations. 
C     The peak energy and time spreading for 3dB or 6dB down points of 
C     the (shaded?) pulse would be useful sonar values. 
C     Thick sub-bottom layers cause a MODES bug (underflow); this can be 
C     patched here by changing the PARAMETER MAXNBL, but MODES should be 
C     corrected.
C     Re HP version: User cannot supply blank file names as defaults 


C------------------------------------------------------------------------------

      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      INCLUDE 'compopp.h'
      INTEGER  NMCAL, LOWER, IERR,ifrq 
      REAL     UN(MAXNM,MAXNSR),
     3     RKN(MAXNM), RKNI(MAXNM), GRPVEL(MAXNM)
      DOUBLE PRECISION WORK(MAXWRK)

C     Additional declarations for OMODES:
      INTEGER I, NMODES
      REAL    RHO(0:MAXNBL), ALPBT(MAXNBL), ROUGH(0:MAXNBL), BH,
     *     AMPL(MAXNL,MAXNM), PHSE(MAXNL,MAXNM),
     *     PRTKNI(MAXNM), PRTKNB(MAXNM), PRTWAT(MAXNM),
     *     PRTSED(MAXNM), PRTSUB(MAXNM),
     *     DBNEP
      DOUBLE PRECISION ZI(MAXNL), CI(MAXNL), KND(MAXNM) 

C     Additional declarations for tolerances: 
      REAL   tolk, tolu, tolgam
      COMMON /tol/ tolk, tolu, tolgam

C     Additional declarations for OGOBIN
      DOUBLE PRECISION PI, HOMEGA
      PARAMETER (pi=3.141592653589793238D0)
c---  peters 
      INTEGER modmax,j
      REAL ssp(maxnl), depth(maxnl),
     *     UNS(MAXnM,MAXNSR), UNR(MAXnM,MAXNSR)

      LUOUT = 6
      tolk = 1.e-15
      tolu = 1.e-10
      tolgam = 1.e-10

C-----------------------------------------------------------------------
C     Preprocessing for calling OMODES

      flagpu=-1
      IF (iopt(6).EQ.1) flagpu=-10

c***************************************************
      ENTRY forw2()
c***************************************************
      lwascomp=1
c     WRITE(*,*) 'has entered popp'
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
c---  last water layer should be waterdepht   
      zp(nwl)=h

      DBNEP = 20000./ALOG(10.)
      NMODES = NMax
      BH = H
      RHO(0) = RHOWC
      ROUGH(0) = ROUGHS
c     
c-----frequency loop
      DO 100  ifrq=1,nfrq
         freq=frq(ifrq)

         DO I = 1, NBL
            RHO(I) = RHOBOT(I)
            ALPBT(I) = ALPBOT(I)*FREQ/DBNEP
            ROUGH(I) = ROUGHB(I)
         ENDDO
c     
c---  WRITE initial model
c     

         IF (flagpu.LE.0) THEN
            WRITE(*,*)'calling popp for freq',freq
            WRITE(prtfil,*)
            WRITE(prtfil,*)'Calling POPP' 
            WRITE(prtfil,*)'flagpu',flagpu
            WRITE(prtfil,*)' Frequency:',freq !     (frq(ifreq),ifreq=1,nfrq) 
            WRITE(prtfil,*)' H,RHOWC', H,RHOWC
c     WRITE(prtfil,*)' SCATT(1),SCATT(2)',SCATT(1),SCATT(2)
c     WRITE(prtfil,*)' H0,H1',H0,H1
            WRITE(prtfil,*)' Nwl,Nbl',Nwl,Nbl
c     WRITE(prtfil,*)' cC0,cC1' ,cC0,cC1
c     WRITE(prtfil,*)' c2,c2s' ,c2,c2s
            WRITE(prtfil,'(a)')'water profile'
            WRITE(prtfil,'(a,/,100(i3,2F10.3/))')
     &           ' ',(i,zp(i),cp(I),i=1,nwl)
            IF (Nwl.GT.0) THEN
               WRITE(prtfil,*)' sediment'
               WRITE(prtfil,'(a,/,100(i3,4F10.3/))')
     &              ' ',(i,hbot(i),rhobot(i),cbot(i),alpbot(I),i=1,nbl)
            ENDIF
            WRITE(prtfil,*)' sd    ',zs(1)     
            WRITE(prtfil,*)' zr    ',(zr(i),i=1,nr)
            WRITE(prtfil,*)' dmu   ',dmu(ifrq)
            WRITE(prtfil,*)' lampow',lampow
         ENDIF


         CALL OMODES (NWL, ZP, CP, BH, H, NBL, HBOT, CBOT, ALPBT, RHO, 
     1        ROUGH, NS, NR, ZS, ZR, NL, ZI, CI, FREQ, MAXNM, NMODES,
     2        LOWER, NMCAL, WORK, KND, RKNI, UN, MAXNL, AMPL, PHSE,
     3        PRTKNI, PRTKNB, PRTWAT, PRTSED, PRTSUB, GRPVEL, IERR)
         IF(nmCAL.EQ.0)  THEN
            WRITE(*,*)' Warning: snapinit found no  modes'
            lwascomp=-1
            RETURN
         ENDIF
         IF (flagpu.LT.20) THEN
            WRITE(prtfil,*)'number of modes',nmcal
         ENDIF 
         
         homega=2d0*pi*freq*h
         DO i=1,nl
            depth(i)= H*ZI(I)
            ssp(i)  = SNGL(HOMEGA/SQRT(CI(I)))
         ENDDO
         modmax = LOWer+NMCAL-1

         DO i=1,nmcal
            DO j=1,ns
               uns(i,j)=un(i,j)
            ENDDO
            DO j=1,nr
               unr(i,j)=un(i,j+ns)
            ENDDO
         ENDDO

         CALL ddeNOGRP (LUINP,MAXNL,MAXNM,maxzs,maxzr,maxnsr,
     *        RHOwc, FREQ, 
     *	      NL, NL-NBL, NS, ZS, NR, ZR,
     *	      LOwer, modmax, KND, rKNI, grpvel, uns,unr, 
     *        AMPL, PHSE,depth,ssp,
     *        isb,tmin,tmax,tinc,tau0,sour_lev(ifrq),nois_lev(ifrq),
     *        dmu(ifrq),flagpu,lampow,ifrq)

 100  CONTINUE
      lwascomp=1
      END


C**************************************
      SUBROUTINE ddenogrp(LUINP, maxlyr, maxmod,maxsrc,maxrcv,maxnsr, 
     *     rhowc, freq,  
     *     nlayer, nlaywc,  nsrc, srcdep, nrcv, rcvdep,
     *     modmin, modmax, rkn, attn, group, uns,unr, 
     *     ampl, phse,depth,ssp,
     *     isb,tmin,tmax,tinc,tau0,ri0,rnoise,
     *     dmu,flagpu,lampow,ifrq)

C     References: For MODBIN: DREA CR/89/444, p.5.
C     For NOGRP:  Notes of Mar 94; also JASA resubmission.

C     AUTHOR:       Dale D Ellis
C     DATE WRITTEN: 28-Jan-91	! As MODBIN
C     LAST EDIT:     5-Jun-94

C     MODIFICATIONS (for NOGRP):

      USE global
      IMPLICIT NONE
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      
      INTEGER  MAXSRC, MAXRCV, MAXNSR, MAXLYR, MAXMOD
      CHARACTER *(*) FMT4G
      PARAMETER (FMT4G = '(1X,I4,2G15.7,I10, 2G15.7)' )

      INTEGER  NLAYER, NLAYWC, NSRC, NRCV,  MODMIN, MODMAX
      REAL     RHOWC, FREQ
      INTEGER  J, NMODES
      REAL   SRCDEP(MAXSRC), RCVDEP(MAXRCV), DEPTH(MAXLYR), 
     *     SSP(MAXLYR),
     *     ATTN(MAXMOD), GROUP(MAXMOD), 
     *     UNS(MAXMOD,MAXNSR), UNR(MAXMOD,MAXNSR), 
     *     AMPL(MAXLYR,MAXMOD), PHSE(MAXLYR,MAXMOD)
      DOUBLE PRECISION  RKN(MAXMOD)
      REAL lampow,rnoise,rgvl_noise
      INTEGER flagpu,ifrq
C------------------------------------------------------------------------------
C     Added declarations for reverberation calculations
c     COMPLEX  cinc, csca, rcoh, pr
      INTEGER  m, n, isb, ism, isrc, ircv, it, ntime
      REAL     tol, pi, cgrp, tau0, dmu, rmu, tmin, tmax, tinc, omega
      REAL     r, area, rkm, cb, vphase !pinc, psca, rinc
      REAL     t, rgvl, sm, sn, rmn, vmn, gmn, cphs, ri0
      REAL     cosm, sinm, gamm, cosn, sinn, gamn ! phasm,, phasn, tl, tlc
      REAL     gm0, gm1, phsem(1000) !phsem(maxmod)
C------------------------------------------------------------------------------
      INTEGER luinp
      
c     DATA NFREQ/1/

c     WRITE (*,*)  NLAYER, NLAYWC, RHOWC, NSRC, NRCV,
c     *        (SRCDEP(I),I=1,NSRC), (RCVDEP(I),I=1,NRCV)
c     WRITE (*,FMT4G)  (I,DEPTH(I),SSP(I), I=1, NLAYER)

c     DO 100 IFREQ = 1, NFREQ
c     WRITE (*,*)  FREQ, MODMIN, MODMAX
c     WRITE(*,*)' nfreq',nfreq
      NMODES = MODMAX-MODMIN+1
      DO 99 J = 1, NMODES
c     WRITE (*,*)  J, RKN(J), GROUP(J), ATTN(J), 
c     *                   (UNS(J,I),I=1,NSRC), (UNR(J,I),I=1,NRCV)

c     WRITE (*,FMT4G)  (I,AMPL(I,J), PHSE(I,J),I=1,NLAYER)
 99   CONTINUE

C------------------------------------------------------------------------------
c     OPEN (20, FORM='FORMATTED', CARRIAGECONTROL='LIST', STATUS='NEW')
c     OPEN (23, FORM='FORMATTED', CARRIAGECONTROL='LIST', STATUS='NEW')
      PI = 4.*ATAN(1.0D0)
      TOL = 1.E-19
      CGRP = 1500.
      CPHS = 1500.
      OMEGA = 2.*PI*FREQ

      IF (ISB.EQ.1)  THEN
         CB = SSP(1)
         ISM = 1
         DO 55 J = 1, NMODES
            PHSEM(J) = -PHSE(ISM,J)*PI/180.
            IF (flagpu.LT.0)
     >           WRITE (prtfil,*)  j, ampl(ism,j), phse(ism,j),
     >           phsem(j)
 55      CONTINUE
      ELSE
         ISM = NLAYWC
         CB = SSP(NLAYWC)
         DO 56 J = 1, NMODES
c     VPHASE = OMEGA/RKN(J)
c     IF (VPHASE.LE.CB)  THEN
C     Possibly a density correction needed; not used anyway.
c     PHSEM(J) = PHSE(ISM+1,J)*PI/180.
c     ELSE
c     PHSEM(J) = PHSE(ISM,J)*PI/180. + DEPTH(ISM+1)/NLAYWC
c     >                            * OMEGA*SQRT(1./CB**2-1./VPHASE**2)
c     ENDIF
            IF (flagpu.LT.0)
     >           WRITE(prtfil,*) j, ampl(ism,j), phse(ism,j), phsem(j)
 56      CONTINUE
      ENDIF

      NTIME = (TMAX-TMIN)/TINC + 0.5
      RMU = SQRT(10.**(DMU/10.))

c     WRITE(*,*)'Number of time samples:',ntime
      ISRC = 1
      IRCV = 1

      DO 90 IT = 1, NTIME
         T = TMIN + (IT-1)*TINC
         R = 0.5*CGRP*(T-0.5*TAU0)
         AREA = PI*R*CGRP*TAU0

c     RINC = TOL**2
c     RCOH = TOL
         RGVL = TOL**2
c     pr = tol
c     tl = tol**2

         DO 80 M =1, NMODES
c     PHASM = R*RKN(M) + PHSEM(M)
c     pr = pr + uns(m,isrc)*(-1)**m*ampl(ism,m)*EXP(-r*attn(m))
c     >              * CMPLX(COS(phasm),SIN(phasm)) /SQRT(rkn(m))
c     tl = tl + (uns(m,isrc)*unr(m,ircv)*EXP(-r*attn(m)))**2/rkn(m)
            VPHASE = OMEGA/RKN(M)
            IF (VPHASE.LE.CB)  GO TO 80
            GAMM = SQRT ( (OMEGA/CB)**2 - RKN(M)**2 )
            COSM = CB/VPHASE
            SINM = SQRT (1. - COSM**2)
            SM  = UNS(M,ISRC)/SQRT(RKN(M))*AMPL(ISM,M)
c     PINC = (-1)**M * SM * EXP(-R*ATTN(M))
c     CINC = PINC * CMPLX (COS(PHASM), SIN(PHASM))
c--   angle 
            DO 70 N = 1, NMODES
               VPHASE = OMEGA/RKN(N)
               IF (VPHASE.LE.CB)  GO TO 70
               GAMN = SQRT((OMEGA/CB)**2-RKN(N)**2)
               COSN = CB/VPHASE
               SINN = SQRT (1. - COSN**2)
c     PHASN = R*RKN(N) + PHSEM(N)
               SN   = UNR(N,IRCV)/SQRT(RKN(N))*AMPL(ISM,N)
c     PSCA = SN * (-1)**N * EXP(-R*ATTN(N))
c     CSCA = PSCA * CMPLX (COS(PHASN), SIN(PHASN))

               RMN = (T-TAU0/2.) / (1./GROUP(M) + 1./GROUP(N) )
               VMN = OMEGA / (RKN(M)+RKN(n))
c---- pg 15/june 1995 lampow=1 gives laberts law 
               GMN = (RMU)*(SINM*SINN)**(Lampow/2)
c     Don't comment these out
               gm0 = gmn
               gm1 = gmn
c     !  Feature scattering
c     if (rmn.gt.20000 .and. rmn.lt.20200)  gm0 = gm0*10.
c     if (r  .gt.20000 .and. r  .lt.20200)  gm1 = gm1*10.
c     *  Slope effect [grpvel only]
*     if (rmn.gt.25000 .and. rmn.lt.25100)  then
*     gm0 = sqrt(rmu*(sinm+0.2)*(sinn+0.2))
*     elseif (rmn.ge.25100 .and. rmn.lt.25200)  then
*     gm0 = sqrt(rmu*max(0.,sinm-0.2)*max(0.,sinn-0.2))
*     endif
c------------eq (16) in ellis, jasa95
c     RINC = RINC + (PINC*PSCA)**2 * gm1**2
c------------eq (14) in ellis, jasa95
c     RCOH = RCOH +  CINC*CSCA * gm1
c------------eq (27) in ellis, jasa95
               RGVL = RGVL + (SM*SN*gm0)**2 / RMN * VMN
     >              * EXP(-2.*RMN*(ATTN(M)+ATTN(N)))

 70         CONTINUE

 80      CONTINUE

c     rinc = ri0* 1./16.*rinc * area*(2.*pi/r)**2
c     rcoh = ri0* 1./16.*rcoh * conjg(rcoh)*area*(2.*pi/r)**2
         rgvl = ri0* 1./16.*rgvl * (2.*pi)**3 * tau0

c     tl  = 10.*alog10(2.*pi/r * tl)
c     tlc = 10.*alog10(2.*pi/r * (cabs(pr))**2)

c     write (20,1010)  t, r/1000., 10.*alog10(rinc), 
c     >              10.*alog10(abs(rcoh)), 10.*alog10(rgvl), tl, tlc
c     1010 FORMAT (2F7.2,5G13.5)
 1011    FORMAT (F7.2, F10.4)
c     write (22,1011)  t, 10.*alog10(rinc)
c     
c     pg 19/11 96 modify rgvl so that it includes the noise level 
c     
         rgvl_noise=rgvl +rnoise
         resp(it+(ifrq-1)*ntime)=10*alog10(rgvl_noise)
c     
 90   CONTINUE


 100  CONTINUE
      if (flagpu.le.0) then
         write (prtfil,*)' revebration level'
         do it=1,ntime
            T = TMIN + (IT-1)*TINC
            write(prtfil,*)t,it+(ifrq-1)*ntime,
     1           real(resp(it+(ifrq-1)*ntime))
         enddo
      endif
      END






