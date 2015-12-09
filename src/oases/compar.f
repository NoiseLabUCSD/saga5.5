      PARAMETER (    NLTYP = 4    )
      PARAMETER (    NLA   = 30,
     &               NLA1  = (NLA-1)*32-16,
     &               NLA2  = 4*NLA+1,
     &               NLA3  = 4*NLA,
     &               NLA4  = 9*NLA+3*NLA3+NLA2+NLA1,
     &               NLA5  = 2*NLA1+2*NLA3          )
      PARAMETER (    IBSI  = 64*NLA)
      PARAMETER (NDV=10,NLA10=(NDV+1)*NLA)
      PARAMETER (    NRD   = 201   )
      PARAMETER (    NMOD  = 10   )
      PARAMETER (    NPEXP=13,
     &               NP    = 2**NPEXP,
     &               NP2   = 2*NP,
     &               ITWONP= 2*NP,
     &               NP3   = 3*NP,
     &               NP4   = 4*NP,
     &               NP6   = 6*NP,
     &               NPHP1 = NP/2+1 )
      PARAMETER (    MODULO=  128 )
      PARAMETER (    ISIZE =  2000 )
      PARAMETER (    NSIP  = 101,
     &               IMXMIN= 20   )
      PARAMETER (    NSMAX = 101   ,
     &               NRMAX = 101   ,
     &               NPMAX = 5    ,
     &               NRNP  = NRMAX*NPMAX,
     &               NRNR  = NRMAX*NRMAX,
     &               NPNP  = NPMAX*NPMAX) 
      PARAMETER (    NCONM = 150  )
C
C     COMMON CONSTANTS
C
      COMPLEX AI,CNUL
      REAL PI
      COMMON /CONST/ PI,CNUL,AI
      COMMON /FREQS/ NFREQ,FREQ1,FREQ2,DLFREQ
      COMPLEX DSQ,CSQ
      COMMON /VARS1/ FREQ,DSQ,CSQ,NUMI,NUML,LS,IR,ZU,ZL,PCORR
      COMMON /VARS2/ NWVNO,ICUT1,ICUT2,ICW1,ICW2
      COMMON /VARS3/ DLWVNO,DLRAN,R1,FNI5,FNIFAC,LMAX,ISPACE,
     &               NPLOTS,NUMFR,IFN
      COMMON /VARS4/ WK0
      COMPLEX WVNO,S2
      COMMON /WAVENO/ WVNO,S2
      COMMON /ERRCOM/ IERR,IOERR
      COMMON /CONTUR/ OFFDB,OFFIMA,OFFDBS
      COMMON /POPT/   KPLOT,NOUT,IOUT(3),IREF,ISTYP,LINA,ICDR,
     &                NDEC,NCON
      COMMON /INTG/   INTTYP
      COMMON /DFLAG/ IDERIV
      LOGICAL SHEAR,DECOMP,SCTOUT
      COMMON /LOGFLG/ SHEAR,DECOMP,SCTOUT
      COMMON /LOGUNI/ LUGRN,LUTRF,LUTGRN,LUTTRF
      CHARACTER*6 PROGNM
      COMMON /SAFANM/ PROGNM
      COMMON /DISNRM/ DISNRM
      logical DEBUG
      COMMON /DBG/ DEBUG,IPRINT,TIMF
      COMMON /OMEGIM/ OMEGIM
      integer flagpu
      common /flagpu/flagpu
      integer iofrhs(nla),lstot
      complex rgrad(nla,4,4*nla+1),ssgrad(nla,4,4*nla+1),ssr(nla,4)
      common /oas_grad_var/iofrhs,rgrad,ssgrad,lstot,ssr
