c
c---- oases variables, as used in input files
c
      INTEGER nlay             !number of layers
      REAL cmin,cmax	       ! min and max phase velocity
      real angle1,angle2
      integer nang
      INTEGER Nwave,IC1,IC2 ! # of wavenumbers, lowerst and highest cut       
      COMMON /oaspar/nlay,cmin,cmax,
     &              nwave,IC1,IC2,ANGLE1,ANGLE2,NANG
      real dtilt
      LOGICAL tilt
      common /tiltparm/tilt,dtilt
