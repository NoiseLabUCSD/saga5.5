c: VARIABLES FROM SVP FILE (read by svp_read):
      common /svp_com0/ nsvp,nlayb,nlayt,ktb(NLMAX),ktt(NLMAX),svp_ver
      integer*4 nsvp,nlayb,nlayt,ktb,ktt
      real*4 svp_ver
      common /svp_com/ 
     .   svp_title,ctol,zsvp(5*NLMAX),csvp(5*NLMAX),rho_svp,alpha_svp,
     .   hb(NLMAX),geob(2,5,NLMAX),bpb(2,NLMAX),
     .   ht(NLMAX),geot(2,5,NLMAX),bpt(2,NLMAX)
      real*8 ctol,zsvp,csvp,rho_svp,alpha_svp,hb,geob,bpb,ht,geot,bpt
      character*64 svp_title
c: VARIABLES FROM OPTION FILE (read by opt_read):
      common /opt_com/ 
c: Line 1: (set iiwrite=1 to output files, 0 to suppress output)
     .   ver_cur,ver_no,iicw,iikpl,iirc,iiparm,n_env,iifmt,iiwrite,
c: Line 2:
     .   iirx,cphmin,cphmax,rmin,rmax,phfac,db_cut,iifb,iigbs,iidiag,
c: Line 3: 
     .   nfcw,fcw(NFBBMAX),
c: Line 4:
     .   iitl,iimf,iisig,iimt,iidc,iikn,iilist,iikrak,iioas,iifepe,
     .   iimlab,iitsp,
c: Lines 5 or 8:
     .   nzs,zsrc(NSRMAX),nrec,zrec(NSRMAX),nsrc,rkm(NSRMAX),
     .   nth_gbs,nb_gbs,th_gbs(NSRMAX),b_gbs(NSRMAX),
c: Line 6:
     .   iiri,iimp,nzmf,zmf(NSRMAX),
c: Line 7 (deleted):
     .   nr_tsp,r1_tsp,r2_tsp,nt_tsp,Tw_tsp,iimex,mr_tsp,pct_tsp,
     .   nrm_tsp,
c: Line 7:
     .   fsbb,Tw,fmindum,fmaxdum,iifft,iiout,iift,
c: Line 9:
     .   nrun,nparmdum,rseed,kvar(4,NPMAX),xvar(2,NPMAX),
c: Line 10 and 11:
     .   nseg,geom_file,lgeom,iitype(NSEGMAX),iicont(NSEGMAX),
     .   vs(NSEGMAX),t1(NSEGMAX),t2(NSEGMAX),cpa(NSEGMAX),
     .   phid(NSEGMAX),x2(NSEGMAX),y2(NSEGMAX),
c: Line 12:
     .   fkpl,iivar,iiform,xkhr1,xkhr2,nreal,xkhi1,xkhi2,nimag,kduc,
     .      iiwr,iishp,iishs,
c: Line 13 and 14:
     .   freq1,freq2,nfreq,iilog,th1,th2,nang,fsrc,nfft
      integer*4 iicw,iikpl,iirc,iiparm,n_env,iifmt,iiwrite,iigbs,
     .   iidiag,iirx,iifb,iitl,iimf,iisig,iimt,nfcw,nzs,nrec,nsrc,
     .   nth_gbs,nb_gbs,iiri,iimp,nzmf,nr_tsp,nt_tsp,iimex,
     .   mr_tsp,nrm_tsp,iifft,iiout,iift,
     .   iidc,iikn,iilist,iikrak,iioas,iifepe,iimlab,iitsp,
     .   nrun,nparmdum,rseed,kvar,nseg,lgeom,iitype,iicont,iivar,
     .   iiform,nreal,nimag,kduc,iiwr,iishp,iishs,nfreq,iilog,nang,
     .   nfft
      real*4 ver_cur,ver_no,cphmin,cphmax,rmin,rmax,phfac,db_cut,fcw,
     .   zsrc,zrec,rkm,th_gbs,b_gbs,zmf,r1_tsp,r2_tsp,Tw_tsp,pct_tsp,
     .   fsbb,Tw,fmindum,fmaxdum,xvar,vs,t1,t2,cpa,phid,x2,y2,fkpl,
     .   xkhr1,xkhr2,xkhi1,xkhi2,freq1,freq2,th1,th2,fsrc
      character*64 geom_file
c
c: OUTPUT VARIABLES:
      common /out_com/ svp_file,opt_file,outroot,outfile,
c: Geoacoustic profile (geo index 1:1=top, 2=bottom;
c: index 2: 1=cp, 2=cs, 3=rho, 4=alphap, 5=alphas; index 3=layer #):
     .   geo(2,5,NLMAX),h(NLMAX),zdep(NLMAX),
c: Complex field array [refer to as plc(1:nzs,1:nsrc,1:nrec)]:
     .   plc(NTLMAX),
c: Mode eigenvalues [refer to kn(1:nmode); kn(0) is dummy]:
     .   kn(0:NM_MAX),
c: Mode functions [refer to phi,dphi,psi,dpsi as phi(1:nzsr,1:nmode)]:
     .   phi(NSR_NM_MAX),dphi(NSR_NM_MAX),psi(NSR_NM_MAX),
     .   dpsi(NSR_NM_MAX),exp_gbs(NSR_NM_MAX),
c: Mode characteristics (row 1=L=ln(R1*R2), 2=dL/dk, 3=dL/dw, 4=vG, 5=R1):
     .   eig_char(5,NM_MAX),
c: Miscellaneous arrays:
     .   zsr(NSRMAX),zsr_im_gbs(NSRMAX),rho_sr(NSRMAX),cp_sr(NSRMAX),
     .   cs_sr(NSRMAX),range(NSNRMAX),
     .   sq2pir(RGMAX,NSNRMAX),ksm2_sr(NSRMAX),
cpln     .   sq2pir(NSNRMAX),ksm2_sr(NSRMAX),
c: Broadband eigenvalues and characteristics:
     .   knbb(NM_NF_MAX),eig_bb(5*NM_NF_MAX),
     .   dtiltp(NSRMAX)
c: Broadband transfer functions [refer to as tf(nfbb,nrec)]:
      common /out_com1/ tf(NTFMAX)
      complex*16 kn,eig_char,knbb,eig_bb,ksm2_sr,cp_sr,cs_sr
      character*64 svp_file,opt_file,outroot,outfile
      real*8 geo,h,zdep,zsr,zsr_im_gbs,rho_sr,range,exp_gbs,
     .       dtiltp
      complex*8 plc,phi,dphi,psi,dpsi,tf,sq2pir
c: Integer*4 and real*4 arrays:
      common /out_com2/
c: nlay=# layers in geo,h,zdep.
     .   nlay,lsvp,lopt,out,lout,loutf,nzref(NM_MAX),kksh(NSRMAX),
c: Number of modes, mode numbers, max # modes for dispersion curves:
     .   nmode,nm_lim,mode_no(3,NM_MAX),nm_miss,mode_phz(3,0:NM_MAX),
     .   nm_cw_max,iishn(NM_MAX),
c: iifail=1 when failure occurs in cw_modes or bb_modes:
     .   iifail,xlam_fb1,xlam_fb2,
c: Arrays for mapping zsr(1:nzsr) to zsrc,zrec,zmf (see opt_com for nsrc, etc):
     .   nzsr,mzsrc(NSRMAX),mzrec(NSRMAX),mzmf(NSRMAX),rng_sr(NSNRMAX),
c: Arrays for mapping sorted ranges to plc(1:nsrc,1:nrec):
     .   nrng,nrec_jr(NSNRMAX),krec_jr(NSNRMAX),jrec_jr(2,NSNRMAX),
c: Broadband nfft, # frequencies, frequency axis:
     .   nfftbb,nfbb,faxbb(NFBBMAX),iish_bb(NM_NF_MAX),
c: TL array [refer to as tl(1:nzs,1:nsrc,1:nrec)]:
     .   tl(NTLMAX),
c: Source/receiver geometry arrays (see opt_com for zsrc,zrec,xsrc):
     .   xsrc(NRNGMAX),ysrc(NRNGMAX),xrec(NSRMAX),yrec(NSRMAX),iigeom,
     .   kn_indx(NM_MAX),zsr_indx(NSNRMAX)
      integer*4 nlay,lsvp,lopt,out,lout,loutf,nmode,nm_lim,mode_no,
     .   nm_miss,nm_cw_max,iishn,iifail,nzref,nzsr,kksh,mzsrc,
     .   mzrec,mzmf,nrng,nrec_jr,krec_jr,jrec_jr,nfftbb,nfbb,iish_bb,
     .   iigeom,kn_indx,zsr_indx
      real*4 tl,faxbb,xsrc,ysrc,xrec,yrec,rng_sr,mode_phz,
     .   xlam_fb1,xlam_fb2
c
      common /bb_com/ wbb(NFBBMAX),nffth1bb,xhbb(20),nf1,nf2,iifull,
     .   df_temp,phibb(NZ_NF_MAX),dpsibb(NZ_NF_MAX),nmbb(NFBBMAX),
     .   kim_bb(NFBBMAX)
      integer*4 nffth1bb,nf1,nf2,iifull,nmbb
      real*8 wbb,df_temp,kim_bb
      real*4 xhbb
      complex*8 phibb,dpsibb
ccc
      common /geo_com/ 
     .   Pcon(3,NLMAX),Qcon(3,NLMAX),Ucon(3,NLMAX),Vcon(3,NLMAX),
     .   Alay(3,2),Blay(3,2),ikcon(3),rholay(2),f_hz,w,wsq,xkh,xkhsq,
     .   xkhrat,xkhratp,xkref,kw,kw0,cref,pie,phtot,gami(3,2,NLMAX),
     .   beti(3,2,NLMAX),xk(2,NLMAX),xb(2,NLMAX),xksq(2,NLMAX),
     .   xbsq(2,NLMAX),eta(NLMAX),etb(NLMAX),etasq(NLMAX),
     .   etbsq(NLMAX),rhorat(NLMAX),gamiref,xbsqinv(2,NLMAX),
     .   xkbp(2,2),xkrat(2,2),xkrat_ref(2),cfmin,csmin,cpfake(2),crmax,
     .   kcrmin,cphlo,cphhi,zduct(NLMAX),rho_duct,twpie,kcut,lncut,
     .   phcut,chspmax,khspmin,f_min,f_max
      complex*16 Pcon,Qcon,Ucon,Vcon,Alay,Blay,ikcon,xkh,xkhsq,xkhrat,
     .   xkhratp,xkref,gami,beti,xk,xb,xksq,xbsq,eta,etb,
     .   etasq,etbsq,xbsqinv,xkbp,xkrat,xkrat_ref,gamiref,kcut,lncut
      real*8 rholay,f_hz,w,wsq,kw,kw0,cref,pie,phtot,rhorat,cfmin,
     .   csmin,cpfake,crmax,kcrmin,cphlo,cphhi,zduct,rho_duct,twpie,
     .   phcut,chspmax,khspmin,f_min,f_max
      common /geo_com2/ isp(NLMAX),iss(NLMAX),iiww(NLMAX),
     .   mm(NLMAX),jflu(2,5),jsol(2,2),nduct,jduct(5,NDMAX),
     .   mzduct(NDMAX),nsvmin,nsvmin0,isvmin,kduct,kduct0,allf(2),
     .   jsurf,jobot,iich,iich_ref,iish(2,2),iish_ref(2),iish0(2,2),
     .   iishr0(2),iisol(NLMAX),jlmin,jlmax,jhsp(2),nhigh,jlfake(2),
     .   iimst,n_int,j_int(NLMAX),mzint(NLMAX),iilk,iiccw,iicut
      integer*4 isp,iss,iiww,mm,jflu,jsol,nduct,jduct,mzduct,nsvmin,
     .   nsvmin0,isvmin,kduct,kduct0,allf,jsurf,jobot,iich,iich_ref,
     .   iish,iish_ref,iish0,iishr0,iisol,jlmin,jlmax,jhsp,nhigh,
     .   jlfake,iimst,n_int,j_int,mzint,iilk,iiccw,iicut
ccc
      common /var_com/ theta(NSRMAX),xh(20),xkhr(NSRMAX),xkhi(NSRMAX),
     .   fftfile
      real*4 theta,xh,xkhr,xkhi
      character*64 fftfile
ccc
      common /vw_com/ Vmat(3,5,NLMAX,2),Wmat(6,NLMAX,2)
      complex*16 Vmat,Wmat
ccc
      common /zsr_com/ jsr2j(NSRMAX),nzmx(NLMAX),nzmxtot,
     .   jzmx(NLMAX),jsrmx(NSRMAX),mx_m(NSRMAX),xmode(NM_MAX),
     .   nctot,ncalc(NM_MAX),nclast,ncall,ncmat,nm_max2,
     .   nrleg(NSEGMAX),t_src(NRNGMAX),rec_lab(NSRMAX),
     .   r4mat1(NHDFMAX),r4mat2(NHDFMAX),aisoln(2,NLMAX)
      integer*4 nzmx,nzmxtot,jzmx,jsrmx,mx_m,nctot,ncalc,nclast,
     .   nrleg,ncall,ncmat,nm_max2,jsr2j,aisoln
      real*4 xmode,r4mat1,r4mat2,t_src,rec_lab
      common /zsr_com2/ zmx(NSRMAX),zmx_im_gbs(NSRMAX),phisr(NSRMAX),
     .   dpsisr(NSRMAX),dkim,kim_min,kim_max,kim_fac,errdkms,errdk2,
     .   errdk100,kremin,phfac0,magfac,magfac0,ph_step,mag_step,
     .   phix(NSRMAX),dphix(NSRMAX),psix(NSRMAX),dpsix(NSRMAX),
     .   expx_gbs(NSRMAX),
     .   xi(NSRMAX),ai(NSRMAX),aip(NSRMAX),bi(NSRMAX),bip(NSRMAX),
     .   xis(NSRMAX),ais(NSRMAX),aips(NSRMAX),bis(NSRMAX),
     .   bips(NSRMAX),zzexp(NSRMAX),zzexps(NSRMAX),r8mat(NHDFMAX)
      real*8 zmx,zmx_im_gbs,dkim,kim_min,kim_max,kim_fac,errdkms,
     .   errdk2,errdk100,kremin,phfac0,magfac,magfac0,
     .   ph_step,mag_step,r8mat,expx_gbs
      complex*8 phisr,dpsisr
      complex*16 phix,dphix,psix,dpsix,xi,ai,aip,bi,bip,xis,ais,aips,
     .   bis,bips,zzexp,zzexps
      common /phz_com/ ailay(2,2,2,NLMAX),bilay(2,2,2,NLMAX),
     .   zetalay(2,2,NLMAX),philay(2,NLMAX),dphilay(2,NLMAX),
     .   psilay(2,NLMAX),dpsilay(2,NLMAX),Aplay(2,NLMAX),
     .   Aslay(2,NLMAX)
      complex*16 ailay,bilay,zetalay,philay,dphilay,psilay,dpsilay
      real*8 Aplay,Aslay
      common /traj_com/ k_sdp,ln_sdp,dln_sdp,k_spt,ln_spt,dln_spt,
     .   k_cont(150),ln_cont(150),dln_cont(150)
      complex*16 k_sdp,ln_sdp,dln_sdp,k_spt,ln_spt,dln_spt,
     .   k_cont,ln_cont,dln_cont
      common /traj_com2/ npt
      integer*4 npt
ccc
      common /rx_com/ dk_max,dk_max0,g11(6,2,NLMAX),g12(6,2,NLMAX),
     .   g21(6,2,NLMAX),g22(6,2,NLMAX),g_exp(6,NLMAX),h11(6,2,NLMAX),
     .   h12(6,2,NLMAX),h21(6,2,NLMAX),h22(6,2,NLMAX),h_exp(6,NLMAX),
     .   u11(6,NLMAX),u12(6,NLMAX),u21(6,NLMAX),u22(6,NLMAX),
     .   u_exp(6,NLMAX),migam2(6),ekn(6,0:NM_MAX),rhofac(NLMAX),
     .   rho_prod(2,NLMAX),Dphix_w(NSRMAX),Ddphix_w(NSRMAX),
     .   Dphi_w(2,NLMAX),Ddphi_w(2,NLMAX),Dphiz_w(NSR_NM_MAX),
     .   Ddphiz_w(NSR_NM_MAX),am(NLMAX),bm(NLMAX),betm(NLMAX),
     .   gm(NLMAX),rhom(NLMAX),xilay(2,NLMAX),cspan(0:NDMAX),
     .   cduct(NDMAX),dz_duct(NDMAX),ccr_lo,ccr_hi,
     .   iia1b1(NLMAX),indx_duct(NDMAX),ndrx,jl2jd(NLMAX),
     .   jpeak(NDMAX),jval(2,NDMAX),jmlink(NM_MAX),nlwz,jlwz(NLMAX),
     .   nzoff,jlayx(NDMAX)
      integer*4 iia1b1,indx_duct,ndrx,jl2jd,jpeak,jmlink,jval,nlwz,
     .   jlwz,nzoff,jlayx
      real*8 dk_max,dk_max0,g11,g12,g21,g22,g_exp,h11,h12,h21,h22,
     .   h_exp,u11,u12,u21,u22,u_exp,migam2,ekn,rhofac,rho_prod,
     .   Dphix_w,Ddphix_w,Dphi_w,Ddphi_w,Dphiz_w,Ddphiz_w,
     .   am,bm,betm,gm,rhom,xilay,cspan,cduct,dz_duct,ccr_lo,ccr_hi
      common /rx_com2/ rx_r1r2(6),R1(6),R2(6),TC1(6),TC2(6),Wr(6)
      complex*16 rx_r1r2,R1,R2,TC1,TC2,Wr
ccc
      common /tiltparm/tiltv,tilth,dtiltv,dtilth
      logical tiltv,tilth
      real dtiltv,dtilth


ccccc these are used for saving modal info

      real zreclo,zrechi,zsrclo,zsrchi,
     .     zwatlo,zwathi,zhtiltlo,zhtilthi,dzrec,dzsrc,
     .     dzwat,dzhtilt,dzlimit
      common /parlimits/ zreclo,zrechi,zsrclo,zsrchi,
     .                   zwatlo,zwathi,zhtiltlo,zhtilthi,
     .                   dzrec,dzsrc,
     .                   dzwat,dzhtilt,dzlimit






