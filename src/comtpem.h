
      TYPE(inputvar) vnp
      TYPE(refractivity) rf
      TYPE(systemvar) sv 
      TYPE(terrain) tr 
      TYPE(errorflag) ef
      TYPE(refparam) rp
      integer flagpu
      common /flagpu/flagpu
      integer fieldtype
       common / tpemcom/ fieldtype,ef,rf,sv,tr,rp
      common /tpemcom2/ vnp
  
      real*8  dr, drout,dzout, dr2
      common / rhstps / dr, drout, dzout, dr2
 
      integer  lvlep_start,itpem_opt(40)
      common  /sagatpem/lvlep_start,itpem_opt
      
      real cap_coef(100)
      integer ncap
      common /capping/cap_coef, ncap

      real znoise
      common /noisefl/znoise

c****** clutter      
      integer Mno  ! max number of nodes
      parameter (Mno=15)
      real elenod(0:Mno+1), ccs(Mno)
      integer Nno  ! number of nodes
      real data_mean,resp_mean
      common /clutterCS/Nno,elenod,ccs,data_mean,resp_mean
c***** beam
      integer mbeam
      parameter (mbeam=10)
      integer nbeam
      real beamelev(mbeam) 
      common /beamdata/ beamelev, nbeam

