      CHARACTER*40 opt
      REAL frq(mfreq),
     &     rd,rdlow,rdep(mdep),rdref(mdep)         ! upper and lower receiver
      REAL sd,sdep(msrd)                  ! source depth
      integer ndep
      integer nsrd
      common /forw_parm/frq,ndep,rd,rdlow,rdep,sd,opt,rdref,sdep,nsrd
      complex resp(Mobs)
      integer nx,nfrq,nbart
      real    xranges(mx)
      integer ixpoints(mx)
      real    xweight(mx)
      real    zsed 
      integer nwater  ! number of points in the water
      common  /zsed/zsed,nwater               !a snap common block
      common /invres/resp
      common /invflags/nx,nfrq,nbart
      common /range/xranges,ixpoints,xweight
      integer lwascomp,iWriteTrf,ierrinfile
      common /logcomp/lwascomp,iWriteTrf,ierrinfile
      real xminimum,xmaxmum,xincre
      common /xvalues/xminimum,xmaxmum,xincre
c time delay
      real del_time
      common /deltime/del_time
c trf
      real dflast,fminlast,fmaxlast
      common /trf_out/dflast,fminlast,fmaxlast



