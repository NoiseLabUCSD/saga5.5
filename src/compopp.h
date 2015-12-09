      INTEGER  MAXNWL, MAXNBL, MAXNM, MAXNL, MAXWRK, 
     1         MAXZS, MAXZR, MAXNSR
      PARAMETER (MAXNWL= 50, MAXZS = 4, MAXZR = 10, MAXNSR = 14, 
     1	         MAXNBL= 10, MAXNM = 500, MAXNL = 1001,
     2	         MAXWRK= 4*MAXNL + MAXNM + MAXNBL+1 )

c      CHARACTER*80 TITLE
      INTEGER LUINP, LUOUT, NWL, NBL, NS, NR, NMAX, NL,NFREQ
      REAL             ZP(MAXNWL), CP(MAXNWL), RHOWC, H, ROUGHS,
     1	               HBOT(MAXNBL), CBOT(MAXNBL), RHOBOT(MAXNBL),
     2	               ALPBOT(MAXNBL), ROUGHB(MAXNBL),
     3                 FREQ, ZS(MAXZS), ZR(MAXZR), 
     4	               ALP2S, C2S
      common /INPMOD/  LUINP, LUOUT,  
     *    NWL, ZP, CP, RHOWC, H, 
     1    ROUGHS, NBL, HBOT, CBOT, RHOBOT, ALPBOT, ROUGHB, C2S, ALP2S, 
     2    FREQ, NS, ZS, NR, ZR, NMAX, NL 

C volume attenuation declaration
      real   attflg,depthinc,alpwc(100)
      common /vol_atten/attflg,depthinc,alpwc
c-- source and noise level
      real sour_lev(mfreq),nois_lev(mfreq)
      common /levels/sour_lev,nois_lev
c--- nogrp input
      integer isb
      real tmin,tmax,tinc,tau0,ri0,dmu(mfreq)
      common /nogrp/isb,tmin,tmax,tinc,tau0,ri0,dmu

c---  various
      integer flagpu
      common /flagpu/flagpu
      real lampow
      common /var/lampow  





