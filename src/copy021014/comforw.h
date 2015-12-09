c      CHARACTER(len=40) opt
      CHARACTER opt(mopt)
      REAL frq(mfreq),
     &     rd,rdlow,rdep(mdep),rdref(mdep) ! upper/lower rec
      REAL sd,sdep(msrd)                  ! source depth
      integer ndep
      integer nsrd
      common /forw_parm/frq,ndep,rd,rdlow,sd,opt,rdref,sdep,nsrd
      common /receiverdepth/rdep
      Complex resp(Mobs)
      integer nx,nfrq,nbart
      real    xranges(mx)
      integer ixpoints(mx)
      real    xweight(mx)
      real    zsed 
      integer nwater  ! number of points in the water
      common  /zsed/zsed,nwater   !a snap common block
      complex resp1(Mobs)
      common /invres/resp, resp1
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



