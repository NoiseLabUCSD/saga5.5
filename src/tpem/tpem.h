      TYPE(errorflag) ef
      real     flagpu
      common /flag/ flagpu
      TYPE(inputvar) vnp
      TYPE(refractivity) rf
      TYPE(systemvar) sv 
      TYPE(terrain) tr 
      TYPE(refparam) rp
      integer  fieldtype
      common / tpemcom/ fieldtype,ef,rf,sv,tr,rp
      common / tpemcom2/ vnp
      real     freqs, FR1,FR2,DT,tmin,tmax
      integer  nx,mx,lx
      character*80 dumtitle

      common / tf_var/ freqs,NX,FR1,FR2,DT,mx,lx,tmin,tmax,dumtitle
      integer luttrf   ! logical unit for trf
      data luttrf/16/
c
      integer itpem_opt(40)
      real frq(100)
      common /tpembb/itpem_opt,frq
c
      integer prtfil   
      data prtfil /7/  
      real znoise
      common /noisefl/znoise
