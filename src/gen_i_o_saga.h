c
c THIS IS A PART OF THE I_O_COM 
      common /out_com1/ tf(NTFMAX),phi_reuse(NFBBMAX,NSR_NM_MAX),
     .       zsrgeom(NSRMAX),iwat0,nzsrgeom
      complex*8 tf,phi_reuse
      real*8 zsrgeom,iwat0
      integer*4 nzsrgeom
c
      common /out_com3/ nzs,nrec,nsrc,nfcw,iicw,iiwrite,iirx,
     .       zsrc(NSRMAX),zrec(NSRMAX),rkm(NSRMAX),fcw(NFBBMAX),
     .       fsbb,fmindum,fmaxdum,rmin,ver_cur
      integer*4 nzs,nrec,nsrc,nfcw,iicw,iiwrite,iirx
      real*4 zsrc,zrec,rkm,fcw,fsbb,fmindum,fmaxdum,rmin,
     .       ver_cur

      common /out_com4/ iifail,jjfail,nfbb,nfftbb,
     .       faxbb(NFBBMAX)
      integer*4 iifail,jjfail,nfbb,nfftbb
      real*4 faxbb
c
c THIS IS A PART OF THE GEN_COM 
      common /bb_com1/ df_temp, nf1, nf2, nmbb(NFBBMAX)
      integer*4 nf1,nf2, nmbb
      real*4 df_temp
c
      common /tiltparm/tiltv,tilth,dtiltv,dtilth
      logical tiltv,tilth
      real dtiltv,dtilth
c
c     logics for determine if modes should be calculated.
      integer i_call_porter,i_geom
      common /log_mode/i_call_porter,i_geom
c
c     Keep track of number of jjfails
      integer tot_jjfail
      common /fails/tot_jjfail
