C This common block is a copy of i_o_1a_com
c where Fmin, fmax nparm is replaced by dummy variables
c
c: i_o_svp_com
c: VARIABLES FROM SVP FILE (read by svp_read):
      common /svp_com1/
     . nsvp,nlayb,nlayt,svp_ver,
     . rend,rstart,f_max
      integer*4 nsvp,nlayb,nlayt
      real*4 svp_ver
      real*8 rend,rstart,f_max
      common /svp_com1a/ ktb(NLMAX),ktt(NLMAX)
      integer*4 ktb,ktt
      common /svp_com2/ 
     .   svp_title,ctol,zsvp(NLMAX),csvp(NLMAX),rho_svp,alpha_svp,
     .   hb(NLMAX),geob(2,5,NLMAX),bpb(2,NLMAX),
     .   ht(NLMAX),geot(2,5,NLMAX),bpt(2,NLMAX)
      real*8 ctol,zsvp,csvp,rho_svp,alpha_svp,hb,geob,bpb,ht,geot,bpt
      character*64 svp_title
c
c: i_o_opt_com
c: VARIABLES FROM OPTION FILE (read by opt_read):
      common /opt_com/ 
c: Lines 5 or 8:
     .   nzs,zsrc(NVRMAX),nrec,zrec(NVRMAX),nsrc,rkm(RGMAX),
     .   zrecusr(NSRMAX),
c: Line 3 
     .   fcw(NFMAX),
     .   nfcw,nrecusr,ncountr,
c: Line 1: (set iiwrite=1 to output files, 0 to suppress output)
     .   ver_no,iicw,iikpl,iirc,n_env,iiwrite,
c: Line 2:
     .   cphmin,cphmax,rmin,rmax,phfac,db_cut,iidiag,iirx,iifb,
c: Line 4:
     .   iimt,iidc,
c: Line 7:
     .   fsbb,Tw,fmindum,fmaxdum,iifft,iiout,iift
      integer*4 iicw,iikpl,iirc,n_env,iiwrite,iidiag,iirx,iifb,
     .   nfcw,nzs,nrec,nsrc,
     .   iifft,iiout,iift,iidc,iimt,
     .   nrecusr,ncountr
      real*4 ver_no,cphmin,cphmax,rmin,rmax,phfac,db_cut,fcw,zsrc,
     .   zrecusr,zrec,rkm,fsbb,Tw,fmindum,fmaxdum
      common /lgc_com/ multcw
      logical multcw 
c
c: gen3_com
      common /bb1_com/ nf1,nf2,df_temp
      integer*4 nf1,nf2
      real*8 df_temp
c
      logical tilt
      real dtilt
      common /tiltparm/tilt,dtilt

      real rdstep
      common /recsep/rdstep

