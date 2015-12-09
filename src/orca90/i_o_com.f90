!---- pg: notice fmin and fmax has been taken out of this common,
!---- so all subroutines taht depends on these do not work!
!---- 
!---- 
MODULE i_o_com
USE parms_com
IMPLICIT NONE
SAVE
!: VARIABLES FROM SVP FILE (read by svp_read):
INTEGER*4 :: nsvp,nlayb,nlayt,ktb(NLMAX),ktt(NLMAX)
REAL*4 :: svp_ver
REAL*8 :: ctol,zsvp(5*NLMAX),csvp(5*NLMAX),rho_svp,alpha_svp,&
        hb(NLMAX),geob(2,5,NLMAX),bpb(4,NLMAX),ht(NLMAX),&
        geot(2,5,NLMAX),bpt(4,NLMAX),htop,Htot,cs_cw_rat,&
        tau_2way,gbar
CHARACTER*64 :: svp_title
!: VARIABLES FROM OPTION FILE (read by opt_read):
!     common /opt_com/ 
!: Line 1: (set iiwrite=1 to output files, 0 to suppress output)
!    .   ver_cur,ver_no,iicw,iikpl,iirc,iiparm,n_env,iifmt,iiwrite,
!: Line 2:
!    .   iirx,cphmin,cphmax,rmin,rmax,phfac,db_cut,iiAih(2),kkblm,
!    .   iigbs,iidiag,
!: Line 3: 
!    .   nfcw,fcw(NFCWMAX),
!: Line 4:
!    .   iitl,iimf,iisig,iimt,iidc,iikn,iilist,iikrak,iioas,iifepe,
!    .   iimlab,iitsp,
!: Lines 5 or 8:
!    .   nzs,zsrc(NSRMAX),nrec,zrec(NSRMAX),nsrc,rkm(NSRMAX),
!    .   nth_gbs,nb_gbs,th_gbs(NSRMAX),b_gbs(NSRMAX),
!: Line 6:
!    .   iiri,iimp,nzmf,zmf(NSRMAX),
!: Line 7 (deleted):
!    .   nr_tsp,r1_tsp,r2_tsp,nt_tsp,Tw_tsp,iimex,mr_tsp,pct_tsp,
!    .   nrm_tsp,
!: Line 7:
!    .   fsbb,Tw,fmin,fmax,iifft,iiout,iift,
!: Line 9:
!    .   nrun,nparm,rseed,kvar(4,NPMAX),xvar(2,NPMAX),
!: Line 10 and 11:
!    .   nseg,geom_file,lgeom,iitype(NSEGMAX),iicont(NSEGMAX),
!    .   vs(NSEGMAX),t1(NSEGMAX),t2(NSEGMAX),dt(NSEGMAX),cpa(NSEGMAX),
!    .   phid(NSEGMAX),x2(NSEGMAX),y2(NSEGMAX),
!: Line 12:
!    .   fkpl,iivar,iiform,xkhr1,xkhr2,nreal,xkhi1,xkhi2,nimag,kduc,
!    .      iiwr,iishp,iishs,
!: Line 13 and 14:
!    .   freq1,freq2,nfreq,iilog,th1,th2,nang,fsrc,nfft
!     integer*4 iicw,iikpl,iirc,iiparm,n_env,iifmt,iiwrite,kkblm,iigbs,
!    .   iidiag,iirx,iitl,iimf,iisig,iimt,nfcw,nzs,nrec,nsrc,
!    .   nth_gbs,nb_gbs,iiri,iimp,nzmf,nr_tsp,nt_tsp,iimex,
!    .   mr_tsp,nrm_tsp,iifft,iiout,iift,
!    .   iidc,iikn,iilist,iikrak,iioas,iifepe,iimlab,iitsp,
!    .   nrun,nparm,rseed,kvar,nseg,lgeom,iitype,iicont,iivar,
!    .   iiform,nreal,nimag,kduc,iiwr,iishp,iishs,nfreq,iilog,nang,nfft
!     real*4 ver_cur,ver_no,cphmin,cphmax,rmin,rmax,phfac,db_cut,iiAih,
!    .   fcw,zsrc,zrec,rkm,th_gbs,b_gbs,zmf,r1_tsp,r2_tsp,Tw_tsp,
!    .   pct_tsp,fsbb,Tw,fmin,fmax,xvar,vs,t1,t2,dt,cpa,phid,x2,y2,fkpl,
!    .   xkhr1,xkhr2,xkhi1,xkhi2,freq1,freq2,th1,th2,fsrc
!     character*64 geom_file
INTEGER :: iicw,iikpl,iirc,iiparm,n_env,iifmt,iiwrite,kkblm,iigbs,&
	iidiag,iirx,iitl,iimf,iisig,iimt,nfcw,nzs,nrec,nsrc,&
	nth_gbs,nb_gbs,iiri,iimp,nzmf,nr_tsp,nt_tsp,iimex,&
	mr_tsp,nrm_tsp,iifft,iiout,iift,iidc,iikn,&
	iilist,iikrak,iioas,iifepe,iimlab,iitsp,nrun,rseed,&
	kvar(4,NPMAX),nseg,lgeom,iitype(NSEGMAX),iicont(NSEGMAX),iivar,&
	iiform,nreal,nimag,kduc,iiwr,iishp,iishs,nfreq,iilog,nang,nfft
! pg	iilist,iikrak,iioas,iifepe,iimlab,iitsp,nrun,nparm,rseed,
!pg fmin&fmax deleted	pct_tsp,fsbb,Tw,xvar(4,NPMAX),vs(NSEGMAX),t1(NSEGMAX),
REAL :: ver_cur,ver_no,cphmin,cphmax,rmin,rmax,phfac,db_cut,iiAih(2),&
	fcw(NFCWMAX),r1_tsp,r2_tsp,Tw_tsp,fminorca,fmaxorca,&
	pct_tsp,fsbb,Tw,xvar(4,NPMAX),vs(NSEGMAX),t1(NSEGMAX),&
	t2(NSEGMAX),dt(NSEGMAX),cpa(NSEGMAX),phid(NSEGMAX),x2(NSEGMAX),&
	y2(NSEGMAX),fkpl,xkhr1,xkhr2,xkhi1,xkhi2,freq1,freq2,&
	th1,th2,fsrc
REAL zsrc(NSRMAX),zrec(NSRMAX),rkm(NRNGMAX),th_gbs(NSRMAX),b_gbs(NSRMAX),zmf(NSRMAX)	! NSRMAX
CHARACTER(len=64) :: geom_file
!
!: OUTPUT VARIABLES:
!     common /out_com/ svp_file,opt_file,outroot,outfile,
!: Geoacoustic profile (geo index 1:1=top, 2=bottom;
!: index 2: 1=cp, 2=cs, 3=rho, 4=alphap, 5=alphas; index 3=layer #):
!    .   geo(2,5,NLMAX),h(NLMAX),zdep(NLMAX),fexp(2,NLMAX),
!: Coherent and incoherent TL [refer to as tlc(1:nzs,1:nsrc,1:nrec)]:
!    .   tlc(NTLMAX),tli(NTLMAX),
!: Mode eigenvalues [refer to kn(1:nmode); kn(0) is dummy]:
!    .   kn(0:NM_MAX),
!: Mode functions [refer to phi,dphi,psi,dpsi as phi(1:nzsr,1:nmode)]:
!    .   phi(NSR_NM_MAX),dphi(NSR_NM_MAX),psi(NSR_NM_MAX),
!    .   dpsi(NSR_NM_MAX),exp_gbs(NSR_NM_MAX),
!: Mode characteristics (row 1=L=ln(R1*R2), 2=dL/dk, 3=dL/dw, 4=vG, 5=R1):
!    .   eig_char(5,NM_MAX),
!: Miscellaneous arrays:
!    .   zsr(NSRMAX),zsr_im_gbs(NSRMAX),rho_sr(NSRMAX),cp_sr(NSRMAX),
!    .   cs_sr(NSRMAX),range(NSNRMAX),
!    .   sq2pir(NSNRMAX),ksm2_sr(NSRMAX),
!: Broadband eigenvalues and characteristics:
!    .   knbb(NM_NF_MAX),eig_bb(5*NM_NF_MAX),
!: Broadband transfer functions [refer to as tf(nfbb,nrec)]:
!    .   tf(NTFMAX)
!     complex*16 kn,eig_char,knbb,eig_bb,ksm2_sr
!     character*64 svp_file,opt_file,outroot,outfile
!     real*8 geo,h,zdep,fexp,zsr,zsr_im_gbs,rho_sr,range,exp_gbs,
!    .   cp_sr,cs_sr
!     complex*8 tlc,phi,dphi,psi,dpsi,tf,sq2pir
!     real*4 tli
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: eig_char
COMPLEX*16,	DIMENSION(:),	ALLOCATABLE :: kn,ksm2_sr
COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: knbb
COMPLEX*16, DIMENSION(:,:,:), ALLOCATABLE :: eig_bb
INTEGER,    DIMENSION(:,:), ALLOCATABLE :: iish_bb
REAL*8,	   DIMENSION(:,:), ALLOCATABLE :: vg_bb,Dvg_w	! for rx_bb
REAL,	   DIMENSION(:,:,:,:), ALLOCATABLE :: phi_bb	        ! for rx_bb
CHARACTER(len=64) :: svp_file,opt_file,outroot,outfile
REAL*8 :: geo(2,5,NLMAX),h(NLMAX),zdep(NLMAX),fexp(2,NLMAX)
REAL*8,	   DIMENSION(:),	ALLOCATABLE :: zsr,zsr_im_gbs,rho_sr,cp_sr,cs_sr	! NSRMAX
REAL*8,	   DIMENSION(:),	ALLOCATABLE :: range	! NSNRMAX
COMPLEX*8, DIMENSION(:,:), ALLOCATABLE :: phi,dphi,psi,dpsi	
REAL*8,	   DIMENSION(:,:), ALLOCATABLE :: exp_gbs
COMPLEX*8, DIMENSION(:,:), ALLOCATABLE :: tlc	! NTLMAX
COMPLEX*8, DIMENSION(:,:,:), ALLOCATABLE :: tf	! NTFMAX
COMPLEX*8, DIMENSION(:),	ALLOCATABLE :: sq2pir	! NSNRMAX
REAL,	   DIMENSION(:,:), ALLOCATABLE :: tli	! NTLMAX


INTEGER*4 :: nlay,lsvp,lopt,out,lout,loutf,nmode,nm_put,nm_tot,&
	nm_lim,nm_miss,nm_cw_max,&
	iifail,iidone,nzsr,nrng,iigeom,&
	nfftbb,nfbb
INTEGER*4,	DIMENSION(:),	ALLOCATABLE :: mzsrc,mzrec,mzmf,nzref,iishn,kn_indx,kksh
INTEGER*4,	DIMENSION(:,:), ALLOCATABLE :: mode_no
INTEGER*4,	DIMENSION(:),	ALLOCATABLE :: nrec_jr,krec_jr,zsr_indx ! NSNRMAX 
INTEGER*4,	DIMENSION(:,:), ALLOCATABLE :: jrec_jr	! (2,NSNRMAX)
REAL*4,		DIMENSION(:,:), ALLOCATABLE :: tl			! NTLMAX
REAL*4,		DIMENSION(:),	ALLOCATABLE :: faxbb		! NFBBMAX
REAL*4,		DIMENSION(:,:), ALLOCATABLE :: mode_phz
REAL*4 :: xsrc(NRNGMAX),ysrc(NRNGMAX)
REAL*4,		DIMENSION(:),	ALLOCATABLE :: xrec,yrec	! NSRMAX
REAL*4,		DIMENSION(:,:), ALLOCATABLE :: rng_sr	! NSNRMAX (nsrc,nrec)
REAL*4 :: xlam_fb1,xlam_fb2
!
!c added by pg 7 march 08
!integer np
!common /divparm/np,ndepth,maxnomode
END MODULE i_o_com
