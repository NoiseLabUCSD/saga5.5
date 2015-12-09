      integer   maxdep,maxrange,msec
      parameter(maxdep=mlay,maxrange=mx,MSEC=25)
c
c---  environmental common
c
      integer nsect
      real R_slngth(msec)
      common /slnght/R_slngth,nsect
      real            R_R0(msec), R_R1(msec), R_R2(msec)
      COMMON /R_DENS/ R_R0, R_R1, R_R2
      real R_BETA(-1:3,msec), R_SCATT(2,msec), R_C2S(msec), R_C2(msec)
      COMMON /R_AB/  R_BETA, R_SCATT, R_C2S, R_C2
      DOUBLE PRECISION R_C0(maxDEP,msec), R_Z0(maxDEP,msec), 
     1                 R_C1(maxDEP,msec),R_Z1(maxDEP,msec)
     2                ,cminsnap,
     1            R_H0(msec), R_H1(msec)
      integer  R_ND0(msec), R_ND1(msec)
      COMMON /R_NA/ R_ND0, R_ND1, CMINsnap
      COMMON /R_G/ R_H0, R_H1, TWOPI, PI, OMEGA
c
      real     RNG(maxRANGE)
c      REAL SRD(mdep,2)
       common /R_rec_sou/rng
c,srd
      double precision freqy
      common /matprop/  R_c0,R_z0,R_c1,R_z1,freqy
      DOUBLE PRECISION TWOPI, PI, OMEGA
      integer np,ndepth,maxnomode
      common /divparm/np,ndepth,maxnomode
c
      real FLDRD,PULRD
      COMMON/CFIELD/FLDRD(3),PULRD(3)
      integer minmod,maxmod,nfreq,nrd
      real hstart
      common /n/ minmod,maxmod,hstart,nfreq
      real rddum(500),sddum,secd(3)
      COMMON/M/RDdum,SDdum,nrd
      integer FLAGPU
      COMMON /FLAGPU/ FLAGPU
      complex pressr(mdep*mx)
      common/pressr/pressr
      COMMON /I/ SECD
      logical tilt,incoh
      real dtilt
      common /tiltparm/tilt,dtilt,incoh

      integer xplane, irflagg
      logical out_plane
      real arrayshape(10)
      common /arrayparm/out_plane,arrayshape,xplane,irflagg

      real xhor(Mdep), yhor(Mdep), zhor(Mdep)
      common /hparmx/xhor
      common /hparmy/yhor
      common /hparmz/zhor
