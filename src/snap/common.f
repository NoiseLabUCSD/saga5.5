COMMON PROC
C
      integer  NOPT,ICF,MODEN,KSRD,nrange,
     &     NDP,
     &     IFLD,MODNP,
     &     NFF,NRP,NR1,NBEG,NPOINT,
     &     mspmax,MAXMSH,mfreq_sn

       PARAMETER ( NOPT=18,ICF=6,MODEN=112,KSRD=1,nrange=15000,
c pg     PARAMETER ( NOPT=18,ICF=6,MODEN=500,KSRD=10,nrange=mx*mdep,
     &     NDP=20010,mfreq_sn=201,
     &     IFLD=20000 )
      PARAMETER (  MODNP=5, NFF=100, NRP=5, NR1=100)
c     p may98      PARAMETER ( NBEG=5, NPOINT=10000+5*NBEG*MODEN+2)
      PARAMETER ( NBEG=5, NPOINT=10000+5*NBEG*MODEN+2)
c     parameter (mspmax=npoint) ! redefined from 3001 by pg 930322
      parameter (mspmax=2001) ! redefined from 3001 by pg 930322
      PARAMETER ( MAXMSH=8 )
C
      DOUBLE PRECISION ADA(NPOINT), SPEED(NPOINT), EIGF(NPOINT),
     &     A3(NPOINT), B3(NPOINT), C3(NPOINT), slow(npoint),
     &     EE(NPOINT), ZZ(NPOINT), SSOLD(NPOINT)
      DOUBLE PRECISION ISO(MODEN), MY(MAXMSH,MODEN)
C
C
      LOGICAL EXCH(NPOINT)
C
      REAL XTS(moden,MSPMAX,mfreq_sn), ALFA(MODEN,mfreq_sn), 
     &     MODAVR(MODEN)      
      DOUBLE PRECISION EK(MODEN,mfreq_sn), EGV(MODEN)
C
      real   FLAGopt(ICF)
      common /flagarray/flagopt
C
C ****
      REAL US(MODEN), UR(MODEN)
      DOUBLE PRECISION TEMPOR(MODEN)
C
C ****
      REAL TLUS(kSrD,MODEN)
C ****      
C FIELD ONLY
      COMPLEX FLDPR(nrange)
      COMPLEX FLDPR1(nrange)
C

c Just to avoid Warnings "variable is unused"

      common /varius/slow,egv,ek,xts,alfa,modavr,exch










