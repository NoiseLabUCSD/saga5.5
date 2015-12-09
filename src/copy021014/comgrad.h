      integer ngrad
      real jacscale(mpar)
      complex grad(Mobs,Mpar)    ! the gradient
      INTEGER    par2phy_forw(mpar),par2lay_forw(mpar) !pointers from parm2phy and lay
      integer nforw_eof,nforw_eofvar
      integer forw2eof(meof),forw2parm(mpar),forw2opt(mpar),nthet
      common  /gradient/grad,ngrad,
     &                par2phy_forw,par2lay_forw,
     &       nforw_eof,nforw_eofvar,forw2eof,forw2parm,forw2opt,
     &       nthet,jacscale
c**** for use in oases grad
      integer irecderiv
      common/oasesgrad/irecderiv
      integer    lay2parm(100,10),Nparmlay(100) !100 is max no of layers
      common /pointergrad/ lay2parm,Nparmlay
