c      include 'fftsiz.inc':
c MXNFFT:  Maximum power of 2 for transform size
c MAXPTS:  Maximum size of arrays for the real and imaginary fields
      
C REMEMBERER ALSO TO UPDATE THESE VARIABLE IN FFTSIZ.INC
      integer*4 maxpts, mxnfft
      parameter ( mxnfft = 12)
      parameter ( maxpts = 2**mxnfft ) 
C REMEMBERER ALSO TO UPDATE THESE VARIABLE IN FFTSIZ.INC
      
      integer*4 maxn4, mxzout, mxrout, mxlvls, mxnprof, mxter    
      
      parameter ( maxn4 = maxpts/4 )  !used for filter array - filters 
                                      !upper 1/4 of field.
      real pi
      parameter ( pi = 3.1415926 )  !Self-explanatory
      parameter ( mxzout = 500 ) !Maximum number of output height points
      parameter ( mxrout = 640 )  !Maximum number of output range points
      parameter ( mxlvls = 600 )  !Maximum number of height/M-unit levels
      parameter ( mxnprof = 40 )  !Maximum number of profiles allowed for
                                  !range-dependent environment.
      parameter ( mxter = 102 )  !Maximum number of height/range points
                                  !allowed for terrain profile
      
c ERRORFLAG:
c  LERR6 =   Logical flag that allows for greater flexibility in allowing error 
c            -6 to be bypassed.  If set to .TRUE. then trapping for this error
c            occurs, otherwise it can be totally ignored by main driver program.
c            (Within TPEM program it is handled as a warning).  If this error is 
c            bypassed (LERR6 = .FALSE.) terrain profile is extended to RMAX with
c            same elevation height of last valid terrain profile point. 
c  LERR12 =  Same as LERR6 - allows for trapping of this error.  If LERR12 = 
c            .FALSE., then (for range-dependent case) if range of last refractivity 
c            profile entered is less than RMAX, the environment is treated as
c            homogeneous from the last profile entered to RMAX.

      TYPE,BIND(C) :: errorflag
         logical lerr6
         logical lerr12
      end TYPE
      
c INPUTVAR:
c   HMAX = maximum output height in meters
c   RMAX = maximum output range in meters  
c   NZOUT = integer number of output height points desired
c   NROUT = integer number of output range points desired
c   PROPANG = Maximum problem (propagation) angle in degrees desired for solution.
c             If set to 0., then TPEM will determine it's own.

      TYPE,BIND(C) :: inputvar
         real*8 hmax
         real*8 rmax
         integer nzout
         integer nrout
         real*8 propang
         real*8 rmin
         real*8 dr
         real*8 hmin
         real*8 dz
      end TYPE
      
c REFRACTIVITY:
c   LVLEP = number of levels in refractivity profile (for range dependent case
c           all profiles must have same number of levels)
c   REFMSL() = 2-dimensional array containing refractivity with respect to mean
c              sea level of each profile.  Array format must be REFMSL(I,J) = 
c              M-unit value at Ith level of Jth profile.  J = 1 for range-
c              independent cases.
c   HMSL() = 2-dimensional array containing heights in meters with respect to
c            mean sea level of each profile.  Array format must be HMSL(I,J) =
c            height of Ith level of Jth profile.  J = 1 for range-independent
c            cases.
c   RNGPROF() = ranges of each profile in meters, i.e., RNGPROF(I) = range of
c               Ith profile.  RNGPROF(1) should always be equal to 0.
c   NPROF = number of profiles.  Equals 1 for range-independent cases.
      
      TYPE,BIND(C) :: refractivity
         integer*4 lvlep
         real refmsl(mxlvls, mxnprof)
         real hmsl(mxlvls, mxnprof)
         real rngprof(mxnprof)
         integer*4 nprof
      end TYPE
      
c SYSTEMVAR:
c   FREQ = frequency in MHz
c   ANTHT = transmitting antenna height above local ground in meters.
c   BWIDTH = half-power (3 dB) antenna pattern beamwidth in degrees (.5 to 45.)
c   ELEV = antenna pattern elevation angle in degrees. (-10 to 10)
c   IPAT = integer value indicating type of antenna pattern desired
c          IPAT = 0 -> omni
c          IPAT = 1 -> gaussian
c          IPAT = 2 -> sinc x
c          IPAT = 3 -> csc**2 x
c          IPAT = 4 -> generic height-finder      
       
      TYPE,BIND(C) :: systemvar
         real freq
         real antht
         real bwidth
         real elev
         integer*4 ipat
      end TYPE
       
c TERRAIN:
c  TERX() = range points of terrain profile in meters
c  TERY() = height points of terrain profile in meters
c  ITP = number of height/range pairs in profile
      
      TYPE,BIND(C) ::  terrain 
         real terx(mxter)
         real tery(mxter)
         integer*4 itp
         integer*4 itpmax
      end TYPE
      
c
c--- for modelling of refraction profile
c
      integer mpcoef,mfactor
      parameter (mpcoef=4,mfactor=10)
      TYPE,BIND(C) :: refparam
        real base(mxnprof)
        real thick(mxnprof)
        real offset(mxnprof)
        real mdef(mxnprof)
        real maxref(mxnprof)
        real coef(mpcoef,mxnprof)
        real factor(	mfactor,1)
        integer npoly
      end TYPE

 


