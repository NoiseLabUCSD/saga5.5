c
c---- This is for the slave program...
c      
      integer  maxslavepar
      parameter (maxslavepar=100)
      integer  nslavepar
      real     slavepar(maxslavepar)
      common /slave_par/nslavepar,slavepar
      integer flagpu
      common /flagpu/flagpu


