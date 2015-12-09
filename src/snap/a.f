
      INTEGER EXTPOL

      CHARACTER*3 modplt
c     CHARACTER*7 STATUS
c     CHARACTER*80 INPFILE, PULSEFILE


      DOUBLE PRECISION  DELFRQ
      DOUBLE PRECISION ROB, ROS

      COMMON /DENS8/ ROB, ROS

      real    TRESH, EPS, RRMAX
      COMMON /EXPMAX/ TRESH, EPS, RRMAX
      real     FACT0, FACT1, FLAGF0, FLAGF1
      COMMON /FACTORS/ FACT0, FACT1, FLAGF0, FLAGF1
      integer flagpu
      COMMON /FLAGPU/ FLAGPU
      real     CORREC
      COMMON /FLAGPULSE/  EXTPOL, CORREC
      real    PLANE, FFP, EK0, SQEK0
      COMMON /FLAGSnap/ PLANE, FFP, EK0, SQEK0

      integer  LUPLP, LUPLT, LUPRT
      COMMON /LUNIT/ LUPLP, LUPLT, LUPRT
      integer MH0(8), MSED(8), ICOUNT 
      real   USEPAST
      COMMON /MESHIST/ MH0, MSED, ICOUNT, USEPAST
      integer   MINMOD, MAXMOD,NFREQ
      real hstart
      COMMON /N/ MINMOD, MAXMOD, HSTART, NFREQ

      integer NFFZ,MSP,NDEPZ,NOPTZ,ICFZ,NDPZ,KSRDZ,
     &        MODENZ,NPOINTZ,MAXMSHZ
      COMMON /PARA1/ NFFZ, MSP, NDEPZ, NOPTZ, ICFZ, NDPZ, KSRDZ,
     &               MODENZ
      COMMON /PARA2/ NPOINTZ, MAXMSHZ
      real    PHVMIN, PHVMAX
      COMMON /PHVEL/ PHVMIN, PHVMAX



      integer JF2,NMES
c      DATA PLANE/0.0/
c      DATA MODPLT/'   '/
c      DATA JF2/ 0/
c      DATA NMES/ 4 /
c      DATA ICOUNT, USEPAST/ 0, 0 /
c      DATA FACT0, FACT1, FLAGF0, FLAGF1/ 1., 1., 0., 0./
c      DATA FLAGPU/0/
c      DATA EXTPOL/ 2 /
c     DATA STATUS/'NEW    '/
c      DATA TRESH/ 1.0E17 /
c      DATA RRMAX/ 0 /
      
       common /snapconst/  JF2, NMES 
       PLANE=0.0
       MODPLT='   '
       JF2= 0
       NMES= 4 
       ICOUNT=0
       USEPAST= 0
       FACT0=1.
       FACT1=1.
       FLAGF0=0.
       FLAGF1=0.
       FLAGPU=0
       EXTPOL= 2 
c     DATA STATUS/'NEW    '/
       TRESH= 1.0E17 
       RRMAX= 0 
 
       PI=DACOS(-1.0D0)
       TWOPI=2.0*PI

C      NFFZ=NFF

c      MSP=MSPMAX
      MSP=301
      NDEPZ=maxdep
      NOPTZ=NOPT
      ICFZ=ICF
      NDPZ=NDP
      KSRDZ=KSRD
      MODENZ=MODEN
      NPOINTZ=NPOINT
      MAXMSHZ=MAXMSH/2
c      NRANGEZ=nrange
c      NR1Z=NR1
c      ISDZ=ISD
c      ISD2Z=ISD2
c      NRPZ=NRP
c      MODNPZ=MODNP
C
      NFREQ= 1
      luprt=prtfil
C     END
