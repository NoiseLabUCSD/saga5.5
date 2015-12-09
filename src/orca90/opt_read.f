      subroutine opt_read
c
c: Reads option file.
c
c	Converted for use of ORCA as a subroutine
c	Variables are either hard-coded herein or passed in
c	call to subroutine orca_sub
c	Unused lines have been commented out by 'co'
c	Note that several input file options will not be recognized,
c	a stop command will be issued.

      use i_o_com
      use gen_com
c-----
      integer*4 nline,j,iierr,nrd,nffth1
      real*4 fac,f1log
      character*64 eline
      data eline/'INVALID INPUT IN OPT FILE: '/

c---------------------------------------------------------
c	Variables now set here 
c---------------------------------------------------------

c	Line 1 - Type of Computations
	ver_no	= 2.01
      iicw	= 1
	iikpl	= 0
	iirc	= 0
	iiparm	= 0
      n_env	= 0
      iifmt	= 0

c	Line 2 - Mode computations
co      iirx	= 0	! 1=Real axis option
	cphmin	= 0
	cphmax	= 0
!	rmin	= 0
!	rmax	= 0
      phfac	= 0
      db_cut	= 0
!	iifb	= 0	! 1=false bottom *** NOT used
!	iiAih(1)= -1.0 ! note: set to -1.0 below
!	iiAih(2)= -1.0 ! note: set to -1.0 below
!	kkblm   = 1    ! note: set to 1 below
	iigbs	= 0
	iidiag	= 0

c	Line 3 - Frequencies
c	nfcw,(fcw(j),j=1,nrd(nfcw))


c	Line 4 - Mode Output Options
      iitl	= 1
	iimf	= 0
	iisig	= 0
	iimt	= 0
	iidc	= 0
      iikn	= 0
      iilist	= 0
	iikrak	= 0
	iioas	= 0
	iifepe	= 0
	iimlab	= 0
c	Line 5 - Source depths, receiver depths, ranges
co               nzs,(zsrc(j),j=1,nrd(nzs))
co     .         nrec,(zrec(j),j=1,nrd(nrec))
co     .         nsrc,(rkm(j),j=1,nrd(nsrc))


c	Line 6 - Mode function depths
co            call star2(46,nline,2,opt_file,1,iiwrite)
co            read(46,*,end=110,err=110) iiri,iimp,nzmf,
co     .         (zmf(j),j=1,nrd(nzmf))




c	----------------------------------------------------------
co      write(2,95)
95    format(/'### OPTION FILE INFORMATION ###')
co      open(46,file=opt_file(1:lopt),form='formatted',err=110)
      nline=0
      iierr=0
      iigeom=0
co      call star2(46,nline,2,opt_file,1,iiwrite)
co      read(46,*,end=110,err=110) ver_no
co      call check_val_r4(ver_no,0.90e0,2.10e0,eline,27,nline,'ver_no',
co     .   6,iierr)
co      backspace(46)
      if(ver_no .lt. 1.10) then
co         read(46,*,end=110,err=110) ver_no,iicw,iikpl,iirc,n_env
         iifmt=1
         iiparm=0
         if(iicw .eq. 3) then
            iiparm=1
            iicw=1
         endif
      else
co         read(46,*,end=110,err=110) ver_no,iicw,iikpl,iirc,iiparm,
co     .      n_env,iifmt
         call check_val_i(iifmt,0,3,eline,27,nline,'iifmt',5,iierr)
         call check_val_i(iiparm,0,1,eline,27,nline,'iiparm',6,iierr)
      endif
      call check_val_i(iicw,0,2,eline,27,nline,'iicw',4,iierr)
      call check_val_i(iikpl,-4,5,eline,27,nline,'iikpl',5,iierr)
      call check_val_i(iirc,0,2,eline,27,nline,'iirc',4,iierr)
      call check_val_i(n_env,0,10000,eline,27,nline,'n_env',5,iierr)
c
co      iiAih(1)=-1.
co      iiAih(2)=-1.
co      kkblm=1
co      call star2(46,nline,2,opt_file,1,iiwrite)
co      if(ver_no .ge. 2.01e0) then
co         read(46,*,end=110,err=110) iirx,cphmin,cphmax,rmin,rmax,
co     .      phfac,db_cut,iiAih(1),iiAih(2),kkblm,iigbs,iidiag
co      elseif(ver_no .ge. 2.0e0) then
co         read(46,*,end=110,err=110) iirx,cphmin,cphmax,rmin,rmax,
co     .      phfac,db_cut,iiAih(1),iiAih(2),iigbs,iidiag
co      elseif(ver_no .ge. 1.6e0) then
co         read(46,*,end=110,err=110) iirx,cphmin,cphmax,rmin,rmax,
co     .      phfac,db_cut,iiAih(1),iigbs,iidiag
co      elseif(ver_no .ge. 1.e0) then
co         read(46,*,end=110,err=110) iirx,cphmin,cphmax,rmin,rmax,
co     .      phfac,db_cut,iiAih(1),iidiag
co      elseif(ver_no .ge. 0.95e0) then
co         read(46,*,end=110,err=110) cphmin,cphmax,rmin,rmax,phfac,
co     .      db_cut,iidiag
co         iirx=0
co         iigbs=0
co      else
co         read(46,*,end=110,err=110) cphmin,cphmax,rmin,rmax,phfac,
co     .      db_cut
co         iidiag=0
co         iirx=0
co         iigbs=0
co      endif
      call check_val_r4(cphmin,-1.e0,1.e10,eline,27,nline,
     .   'cphmin',6,iierr)
      call check_val_r4(cphmax,-89.99e0,1.e10,eline,27,nline,
     .   'cphmax',6,iierr)
      call check_val_r4(rmin,-10000.,10000.,eline,27,nline,
     .   'rmin',4,iierr)
      call check_val_r4(rmax,0.e0,1.e10,eline,27,nline,
     .   'rmax',4,iierr)
      call check_val_r4(phfac,0.e0,512.e0,eline,27,nline,
     .   'phfac',5,iierr)
      call check_val_r4(db_cut,0.e0,120.e0,eline,27,nline,
     .   'db_cut',6,iierr)
      call check_val_i(iirx,-1,2,eline,27,nline,'iirx',4,iierr)
      call check_val_r4(iiAih(1),-1.,361.,eline,27,nline,'iiAih(1)',
     .   7,iierr)
      call check_val_r4(iiAih(2),-1.,361.,eline,27,nline,'iiAih(2)',
     .   7,iierr)
      call check_val_i(kkblm,0,1,eline,27,nline,'kkblm',5,iierr)
      call check_val_i(iigbs,0,1,eline,27,nline,'iigbs',5,iierr)
      if(cphmax .ne. 0. .and. rmin .lt. 999.) then
         iierr=1
         print *,'Limit # modes using cphmax or rmin, but not both!'
         print *,'Set cphmax=0 to disable phase speed limit OR '
         print *,'Set rmin>=999 to disable minimum range limit.'
      elseif(cphmax .eq. 0. .and. rmin .ge. 999.) then
         iierr=1
         print *,'Must use cphmax or rmin to limit # modes!'
         print *,'Set cphmax nonzero to enable phase speed limit OR '
         print *,'Set rmin<999 to enable minimum range limit.'
      endif
c
      f_min=1.d20
      f_max=0.d0
c
      if(iicw .eq. 1) then
         iikn=0
co         call star2(46,nline,2,opt_file,1,iiwrite)
co         read(46,*,end=110,err=110) nfcw,(fcw(j),j=1,nrd(nfcw))
co         call star2(46,nline,2,opt_file,1,iiwrite)
         if(ver_no .ge. 1.e0) then
co            read(46,*,end=110,err=110) iitl,iimf,iisig,iimt,iidc,
co     .         iikn,iilist,iikrak,iioas,iifepe,iimlab
            call check_val_i(iikn,0,1,eline,27,nline,'iikn',4,iierr)
            call check_val_i(iikrak,0,100,eline,27,nline,'iikrak',6,
     .         iierr)
            call check_val_i(iioas,0,5,eline,27,nline,'iioas',5,
     .         iierr)
            call check_val_i(iifepe,0,2,eline,27,nline,'iifepe',6,
     .         iierr)
            call check_val_i(iimlab,0,1,eline,27,nline,'iimlab',6,
     .         iierr)
            iitsp=0
         elseif(ver_no .ge. 0.96e0) then
co            read(46,*,end=110,err=110) iitl,iimf,iisig,iimt,iidc,
co     .         iilist,iitsp,iikrak,iikn
            call check_val_i(iikn,0,1,eline,27,nline,'iikn',4,iierr)
            call check_val_i(iikrak,-3,100,eline,27,nline,'iikrak',6,
     .         iierr)
         elseif(ver_no .ge. 0.95e0) then
co            read(46,*,end=110,err=110) iitl,iimf,iisig,iimt,iidc,
co     .         iilist,iitsp,iikrak
            call check_val_i(iikrak,-1,100,eline,27,nline,'iikrak',6,
     .         iierr)
         elseif(ver_no .ge. 0.92e0) then
co            read(46,*,end=110,err=110) iitl,iimf,iisig,iimt,iidc,
co     .         iilist,iitsp,iidiag,iikrak
            call check_val_i(iikrak,-1,100,eline,27,nline,'iikrak',6,
     .         iierr)
         elseif(ver_no .ge. 0.91e0) then
co            read(46,*,end=110,err=110) iitl,iimf,iisig,iimt,iidc,
co     .         iilist,iitsp,iidiag
            iikrak=0
         else
            print *,'ver_no<0.91 no longer supported'
            stop
         endif
c
         call mem_lim(nfcw,NFCWMAX,MLINE,LML,'nfcw',4,'NFCWMAX',7,
     .      1,1)
         call check_val_i(iitl,-3,3,eline,27,nline,'iitl',4,iierr)
         call check_val_i(iimf,0,4,eline,27,nline,'iimf',4,iierr)
         call check_val_i(iisig,0,1,eline,27,nline,'iisig',5,iierr)
         call check_val_i(iimt,0,1,eline,27,nline,'iimt',4,iierr)
         call check_val_i(iidc,0,3,eline,27,nline,'iidc',4,iierr)
         call check_val_i(iilist,0,1,eline,27,nline,'iilist',6,iierr)
         call check_val_i(iitsp,0,1,eline,27,nline,'iitsp',5,iierr)
         call check_val_i(iidiag,-10,3,eline,27,nline,'iidiag',6,
     .      iierr)
         do j=1,nrd(nfcw)
            call check_val_r4(fcw(j),0.e0,1.e10,eline,27,nline,
     .         'fcw(j)',6,iierr)
         enddo
         f_min=fcw(1)
         f_max=fcw(1)
         do j=2,nrd(nfcw)
            f_min=amin1(fcw(j),sngl(f_min))
            f_max=amax1(fcw(j),sngl(f_max))
         enddo
c
         if(iabs(iitl) .eq. 1 .or. iabs(iitl) .eq. 3) then
            iigeom=1
co            call star2(46,nline,2,opt_file,0,iiwrite)
            if(iigbs .eq. 0) then
co               read(46,*,end=110,err=110) nzs,(zsrc(j),
co     .            j=1,nrd(nzs)),nrec,(zrec(j),j=1,nrd(nrec)),
co     .            nsrc,(rkm(j),j=1,nrd(nsrc))
co               write(2,202) nzs,(zsrc(j),j=1,nrd(nzs))
co               write(2,202) nrec,(zrec(j),j=1,nrd(nrec))
co               write(2,202) nsrc,(rkm(j),j=1,nrd(nsrc))
            else
			 stop ' input for iigbs=1 not supported '
               read(46,*,end=110,err=110) nzs,(zsrc(j),
     .            j=1,nrd(nzs)),nrec,(zrec(j),j=1,nrd(nrec)),
     .            nsrc,(rkm(j),j=1,nrd(nsrc)),nth_gbs,(th_gbs(j),
     .            j=1,nrd(nth_gbs)),nb_gbs,(b_gbs(j),j=1,nrd(nb_gbs))
               write(2,202) nzs,(zsrc(j),j=1,nrd(nzs))
               write(2,202) nrec,(zrec(j),j=1,nrd(nrec))
               write(2,202) nsrc,(rkm(j),j=1,nrd(nsrc))
               write(2,202) nth_gbs,(th_gbs(j),j=1,nrd(nth_gbs))
               write(2,202) nb_gbs,(b_gbs(j),j=1,nrd(nb_gbs))
               call check_val_i(iabs(nth_gbs),iabs(nzs),iabs(nzs),
     .            eline,27,nline,'nzs for Gaussian beam source angle',
     .            34,iierr)
               call check_val_i(iabs(nb_gbs),iabs(nzs),iabs(nzs),
     .            eline,27,nline,'nzs for Gaussian beam beamwidth',
     .            31,iierr)
            endif
202         format(i4,2x,7(f9.3,1x)/1000(6x,7(f9.3,1x)/))
            call check_val_i(iabs(nzs),1,NSRMAX,eline,27,nline,
     .         'nzs',7,iierr)
            call check_val_i(iabs(nrec),1,NSRMAX,eline,27,nline,
     .         'nrec',7,iierr)
            call check_val_i(iabs(nsrc),1,NRNGMAX,eline,27,nline,
     .         'nsrc',4,iierr)
         elseif(iabs(iitl) .eq. 2) then
co           call star2(46,nline,2,opt_file,0,iiwrite)
            iigeom=2
            if(iigbs .eq. 1) then
               print *,'Gaussian Beam Source not implemented for ',
     .            'array geometry file. Setting to zero'
               iigbs=0
            endif
         else
co            call star2(46,nline,2,opt_file,0,iiwrite)
            iigbs=0
         endif

         if(iimf .ne. 0 .or. iisig .ne. 0) then
co            call star2(46,nline,2,opt_file,1,iiwrite)
co            read(46,*,end=110,err=110) iiri,iimp,nzmf,
co     .         (zmf(j),j=1,nrd(nzmf))
            call check_val_i(iiri,0,3,eline,27,nline,'iiri',4,iierr)
            call check_val_i(iimp,0,3,eline,27,nline,'iimp',4,iierr)
            call check_val_i(iabs(nzmf),0,NSRMAX,eline,27,nline,
     .         'nzmf',4,iierr)
         else
co            call star2(46,nline,2,opt_file,0,iiwrite)
         endif
c
         if(iitsp .eq. 1) then
			 stop ' input for iitsp=1 not supported '
            call star2(46,nline,2,opt_file,1,iiwrite)
            read(46,*,end=110,err=110) nr_tsp,r1_tsp,r2_tsp,
     .         nt_tsp,Tw_tsp,iimex,mr_tsp,pct_tsp,nrm_tsp
            call check_val_i(nr_tsp,1,10000,eline,27,nline,
     .         'nr_tsp',6,iierr)
            call check_val_r4(r1_tsp,0.e0,r2_tsp,eline,27,nline,
     .         'r1_tsp',6,iierr)
            call check_val_r4(r2_tsp,r1_tsp,1.e10,eline,27,nline,
     .         'r2_tsp',6,iierr)
            call check_val_r4(Tw_tsp,1.e-6,1.e10,eline,27,nline,
     .         'Tw_tsp',6,iierr)
            call check_val_i(iimex,0,1,eline,27,nline,'iimex',5,
     .         iierr)
            call check_val_i(mr_tsp,1,3000,eline,27,nline,
     .         'mr_tsp',6,iierr)
            call check_val_r4(pct_tsp,1.e0,100.e0,eline,27,nline,
     .         'pct_tsp',7,iierr)
            call check_val_i(nrm_tsp,0,1,eline,27,nline,
     .         'nrm_tsp',7,iierr)
         elseif(ver_no .lt. 1.e0) then
c: For ver_no < 1.0, skip line of inputs for time-spread plot:
co            call star2(46,nline,2,opt_file,0,iiwrite)
         endif
      else
			 stop ' input for iicw NE 1 not supported '
         call star2(46,nline,2,opt_file,1,iiwrite)
cc       read(46,*,end=110,err=110) nfcw,(fcw(j),j=1,nrd(nfcw))
         call star2(46,nline,2,opt_file,1,iiwrite)
         if(ver_no .ge. 1.e0) then
            read(46,*,end=110,err=110) iitl,iimf,iisig,iimt,iidc,
     .         iikn,iilist,iikrak,iioas,iifepe,iimlab
         endif
         do j=1,2
            call star2(46,nline,2,opt_file,0,iiwrite)
         enddo
         if(ver_no .lt. 1.e0) then
c: For ver_no < 1.0, skip line of inputs for time-spread plot:
            call star2(46,nline,2,opt_file,0,iiwrite)
         endif
      endif
c
      if(iicw .eq. 2) then
			 stop ' input for iicw NE 1 not supported '
         call star2(46,nline,2,opt_file,1,iiwrite)
         if(ver_no .ge. 1.01) then
            read(46,*,end=110,err=110) fsbb,Tw,fminorca,fmaxorca,iifft,
     .         iiout,iift,iimt,iidc,iimf
            call check_val_i(iimf,0,4,eline,27,nline,'iimf',4,iierr)
         else
            read(46,*,end=110,err=110) fsbb,Tw,fminorca,fmaxorca,iifft,
     .         iiout,iift,iimt,iidc
            iimf=0
         endif
         call check_val_r4(fsbb,0.e0,1.e10,eline,27,nline,
     .      'fsbb',4,iierr)
         call check_val_r4(Tw,-1.e10,131072e0,eline,27,nline,
     .      'nfft/Tw',7,iierr)
         if(fmaxorca .lt. 0.d0) then
            read(46,*,end=110,err=110) nfcw,(fcw(j),j=1,nrd(nfcw))
            do j=1,nrd(nfcw)
               call check_val_r4(fcw(j),0.e0,1.e10,eline,27,nline,
     .            'fcw(j)',6,iierr)
            enddo
            f_min=fcw(1)
            f_max=fcw(1)
            do j=2,nrd(nfcw)
               f_min=amin1(fcw(j),sngl(f_min))
               f_max=amax1(fcw(j),sngl(f_max))
            enddo
            call uni_space(nfcw,fcw,1.e0)
c: Make sure we do bb_brute or rx_bb_brute at list of frequencies:
            if(iirx .ne. 0) iirx=2
c: Make sure we do not try to output an FFT for a list of frequencies:
            if(iifft .ne. 0) then
               print *,'iifft set to 0 for list of frequencies'
               write(2,*) 'iifft set to 0 for list of frequencies'
               iifft=0
            endif
         else
            call check_val_r4(fminorca,1.e-3,fmaxorca,
     .         eline,27,nline,'fmin',4,iierr)
            call check_val_r4(fmaxorca,fminorca,fsbb/2.e0,
     .         eline,27,nline,'fmax',4,iierr)
            f_min=fminorca
            f_max=fmaxorca
         endif
         call check_val_i(iifft,0,2,eline,27,nline,'iifft',5,iierr)
         call check_val_i(iiout,0,2,eline,27,nline,'iiout',5,iierr)
         call check_val_i(iift,0,5,eline,27,nline,'iift',4,iierr)
         call check_val_i(iimt,0,1,eline,27,nline,'iimt',4,iierr)
         call check_val_i(iidc,0,3,eline,27,nline,'iidc',4,iierr)
         if(iifft*iiout .ne. 0 .and. iifft .ne. iiout) then
            print *,'Require iifft=iiout when iifft>0 and iiout>0.'
            print *,'Assuming zs,zr,r to be read according to iifft.'
            iiout=iifft
         endif
c
         call star2(46,nline,2,opt_file,0,iiwrite)
         if(iifft .eq. 1 .or. iiout .eq. 1 .or. iimf .ne. 0) then
            iigeom=1
            if(iigbs .eq. 0) then
               read(46,*,end=110,err=110) nzs,(zsrc(j),
     .            j=1,nrd(nzs)),nrec,(zrec(j),j=1,nrd(nrec)),
     .            nsrc,(rkm(j),j=1,nrd(nsrc))
               write(2,202) nzs,(zsrc(j),j=1,nrd(nzs))
               write(2,202) nrec,(zrec(j),j=1,nrd(nrec))
               write(2,202) nsrc,(rkm(j),j=1,nrd(nsrc))
            else
               read(46,*,end=110,err=110) nzs,(zsrc(j),
     .            j=1,nrd(nzs)),nrec,(zrec(j),j=1,nrd(nrec)),
     .            nsrc,(rkm(j),j=1,nrd(nsrc)),nth_gbs,(th_gbs(j),
     .            j=1,nrd(nth_gbs)),nb_gbs,(b_gbs(j),j=1,nrd(nb_gbs))
               write(2,202) nzs,(zsrc(j),j=1,nrd(nzs))
               write(2,202) nrec,(zrec(j),j=1,nrd(nrec))
               write(2,202) nsrc,(rkm(j),j=1,nrd(nsrc))
               write(2,202) nth_gbs,(th_gbs(j),j=1,nrd(nth_gbs))
               write(2,202) nb_gbs,(b_gbs(j),j=1,nrd(nb_gbs))
               call check_val_i(iabs(nth_gbs),iabs(nzs),iabs(nzs),
     .            eline,27,nline,'nzs for Gaussian beam source angle',
     .            34,iierr)
               call check_val_i(iabs(nb_gbs),iabs(nzs),iabs(nzs),
     .            eline,27,nline,'nzs for Gaussian beam beamwidth',
     .            31,iierr)
            endif
            call check_val_i(iabs(nzs),1,NSRMAX,eline,27,nline,
     .         'nzs',7,iierr)
            call check_val_i(iabs(nrec),1,NSRMAX,eline,27,nline,
     .         'nrec',7,iierr)
            call check_val_i(iabs(nsrc),1,NRNGMAX,eline,27,nline,
     .         'nsrc',4,iierr)
         elseif(iifft .eq. 2 .or. iiout .eq. 2) then
            iigeom=2
            if(iigbs .eq. 1) then
               print *,'Gaussian Beam Source not implemented for ',
     .            'array geometry file. Setting to zero'
               iigbs=0
            endif
         endif
      else
co         call star2(46,nline,2,opt_file,0,iiwrite)
co         call star2(46,nline,2,opt_file,0,iiwrite)
      endif

      if(iiparm .eq. 1) then
			 stop ' input for iiparm = 1 not supported '
cpg         call star2(46,nline,2,opt_file,0,iiwrite)
c: Make sure we write out nparm lines to output file:
cpg          read(46,*,end=110,err=110) nrun,nparm
cpg          backspace(46)
cpg          backspace(46)
cpg         backspace(2)
cpg         call star2(46,nline,2,opt_file,nparm,iiwrite)
cpg          if(ver_no .lt. 1.5) then
cpg             read(46,*,end=110,err=110) nrun,nparm,(kvar(1,j),
cpg      .         kvar(2,j),kvar(3,j),kvar(4,j),xvar(1,j),xvar(2,j),
cpg      .         j=1,min(NPMAX,iabs(nparm)))
cpg          else
cpg             read(46,*,end=110,err=110) nrun,nparm,rseed,(kvar(1,j),
cpg      .         kvar(2,j),kvar(3,j),kvar(4,j),xvar(1,j),xvar(2,j),
cpg      .         j=1,min(NPMAX,iabs(nparm)))
cpg          endif
cpg         call check_val_i(nrun,-1000,1000,eline,27,nline,
cpg      .      'nrun',4,iierr)
cpg          call check_val_i(nparm,1,NPMAX,eline,27,nline,
cpg      .      'nparm',5,iierr)
cpg       else
co         call star2(46,nline,2,opt_file,0,iiwrite)
      endif
      if(iigeom .eq. 2) then
	 stop ' input for iigeom=2 not supported '
cpg          call star2(46,nline,2,opt_file,1,iiwrite)
cpg          nzs=1
cpg          read(46,*,end=110,err=110) zsrc(1),nseg,geom_file
cpg          call lname64(geom_file,lgeom)
cpg          call check_val_i(nseg,1,NSEGMAX,eline,27,nline,
cpg      .      'nseg',4,iierr)
cpg          call src_track(nline,eline,iierr)
      else
co         call star2(46,nline,2,opt_file,0,iiwrite)
co         call star2(46,nline,2,opt_file,0,iiwrite)
      endif
      if(iikpl .ne. 0) then
			 stop ' input for iikpl > 0 not supported '
         call star2(46,nline,2,opt_file,1,iiwrite)
         read(46,*,end=110,err=110) fkpl,iivar,iiform,xkhr1,xkhr2,
     .      nreal,xkhi1,xkhi2,nimag,kduc,iiwr,iishp,iishs
         call check_val_r4(fkpl,1.e-3,1.e10,eline,27,nline,
     .      'f for iikpl=1',11,iierr)
         call check_val_i(iivar,1,3,eline,27,nline,'iivar',5,iierr)
         call check_val_i(iiform,1,3,eline,27,nline,'iikf',6,iierr)
         call check_val_r4(xkhr1,-1.e10,xkhr2,eline,27,nline,
     .      'kr1',5,iierr)
         call check_val_r4(xkhr2,xkhr1,1.e10,eline,27,nline,
     .      'kr2',5,iierr)
         call check_val_i(nreal,1,100000,eline,27,nline,'nkr',3,iierr)
         call check_val_r4(xkhi1,-1.e10,1.e10,eline,27,nline,
     .      'ki1',5,iierr)
         call check_val_r4(xkhi2,-1.e10,1.e10,eline,27,nline,
     .      'ki2',5,iierr)
         call check_val_i(nimag,1,100000,eline,27,nline,'nki',3,iierr)
         call check_val_i(kduc,0,25,eline,27,nline,'nduct',5,iierr)
         call check_val_i(iiwr,1,2,eline,27,nline,'iiph',4,iierr)
         call check_val_i(iishp,-1,1,eline,27,nline,'iishp',5,iierr)
         call check_val_i(iishs,-1,1,eline,27,nline,'iishs',5,iierr)
         f_max=amax1(sngl(f_max),fkpl)
         f_min=amin1(sngl(f_min),fkpl)
      else
         nreal=0
         nimag=0
co         call star2(46,nline,2,opt_file,0,iiwrite)
      endif
      if(iirc .eq. 1) then
			 stop ' input for iirc EQ 1 not supported '
         call star2(46,nline,2,opt_file,1,iiwrite)
         if(ver_no .lt. 2.1) then
            read(46,*,end=110,err=110) freq1,freq2,nfreq,iilog,
     .         th1,th2,nang
            call mem_lim(nfcw,NFCWMAX,MLINE,LML,'nfreq',5,'NFCWMAX',7,
     .         1,1)
            call mem_lim(nang,NSRMAX,MLINE,LML,'nang',4,'NSRMAX',6,
     .         1,1)
            call check_val_r4(freq1,1.e-3,freq2,eline,27,nline,
     .         'freq1',5,iierr)
            call check_val_r4(freq2,freq1,1.e10,eline,27,nline,
     .         'freq2',5,iierr)
            call check_val_i(iilog,0,1,eline,27,nline,
     .         'iilog',5,iierr)
            call check_val_r4(th1,0.e0,90.e0,eline,27,nline,
     .         'theta1',5,iierr)
            call check_val_r4(th2,0.e0,90.e0,eline,27,nline,
     .         'theta2',5,iierr)
            if(iilog .eq. 1) then
               f1log=log10(freq1)
               fac=(log10(freq2)-f1log)/max(nfreq-1,1)
               do j=1,nfreq
                  fcw(j)=10.**(f1log + (j-1)*fac)
               enddo
            else
               fac=(freq2-freq1)/max(nfreq-1,1)
               do j=1,nfreq
                  fcw(j)=freq1 + (j-1)*fac
               enddo
            endif
            fac=(th2-th1)/max(nang-1,1)
            do j=1,nang
               theta(j)=th1 + (j-1)*fac
            enddo
         else
            read(46,*,end=110,err=110) nfreq,(fcw(j),j=1,
     .         nrd(nfreq)),nang,(theta(j),j=1,nrd(nang))
            call mem_lim(iabs(nfreq),NFCWMAX,MLINE,LML,'nfreq',5,
     .         'NF_RC',5,1,1)
            call mem_lim(iabs(nang),NSRMAX,MLINE,LML,'nang',4,
     .         'NSRMAX',6,1,1)
            do j=1,nrd(nfreq)
               call check_val_r4(fcw(j),0.e0,1.e10,eline,27,nline,
     .            'f_rc(j)',7,iierr)
               f_min=amin1(fcw(j),sngl(f_min))
               f_max=amax1(fcw(j),sngl(f_max))
            enddo
            call uni_space(nfreq,fcw,1.e0)
            call uni_space(nang,theta,1.e0)
            iilog=0
         endif
         call check_val_i(nfreq,1,NFCWMAX,eline,27,nline,
     .      'nfreq',5,iierr)
         call check_val_i(nang,1,NSRMAX,eline,27,nline,
     .      'ntheta',6,iierr)
         f_max=amax1(sngl(f_max),freq2)
         f_min=amin1(sngl(f_min),freq1)
      else
co        call star2(46,nline,2,opt_file,0,iiwrite)
      endif
      if(iirc .eq. 2) then
			 stop ' input for iirc EQ 2 not supported '
         call star2(46,nline,2,opt_file,1,iiwrite)
         read(46,*,end=110,err=110) freq1,freq2,fsrc,nfft,
     .      th1,th2,nang
         call check_val_r4(freq1,1.e-3,freq2,eline,27,nline,
     .      'freq1',5,iierr)
         call check_val_r4(freq2,freq1,fsrc/2.e0,eline,27,nline,
     .      'freq2',5,iierr)
         call check_val_r4(fsrc,1.e-3,1.e10,eline,27,nline,
     .      'fsrc',4,iierr)
         call check_val_i(nfft,32,131072,eline,27,nline,
     .      'nfft',5,iierr)
         call check_val_r4(th1,0.e0,90.e0,eline,27,nline,
     .      'theta1',5,iierr)
         call check_val_r4(th2,0.e0,90.e0,eline,27,nline,
     .      'theta2',5,iierr)
         call check_val_i(nang,1,10000,eline,27,nline,
     .      'ntheta',6,iierr)
         f_max=amax1(sngl(f_max),freq2)
         f_min=amin1(sngl(f_min),freq1)
c
         df=fsrc/nfft
         nffth1=nfft/2 + 1
         nf1=max(2,nint(freq1/df + 1))
         freq1=df*(nf1-1)
         nf2=min(nffth1,int(freq2/df) + 1)
         freq2=df*(nf2-1)
         nfreq=nf2-nf1+1
      else
co         call star2(46,nline,2,opt_file,0,iiwrite)
      endif
c
      if(iigeom .eq. 2) call arr_geom(iierr)
c
      if(iierr .eq. 1) then
         print *,' '
         print *,'Execution terminating.  Check input option file '//
     .      'for error(s).'
         stop
      endif

      RETURN
c
110   print *,'Error opening or reading option file ',opt_file(1:lopt)
      print *,'*() line number = ',nline
      stop
c
      end
ccc
      subroutine src_track(nline,eline,iierr)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 nline,iierr,nrden,kr,j
      real*4 r0,tlast,cosphi,sinphi,pierad,fac,xxx0,yyy0,x1x,y1x,
     .   x2x,y2x,delx,dely,delr,delt,time0
      character*64 eline
      data pierad/0.01745329251994/
c
      call star2(46,nline,2,opt_file,nseg,iiwrite)
      nsrc=0
      do j=1,nseg
         read(46,*,end=110,err=110) iitype(j),iicont(j),
     .      vs(j),t1(j),t2(j),dt(j),cpa(j),phid(j),x2(j),y2(j)
         call check_val_i(iitype(j),1,2,eline,27,nline,
     .      'iitype',6,iierr)
         call check_val_i(iicont(j),0,1,eline,27,nline,
     .      'iicont',6,iierr)
cxx      if(iimet .eq. 0) then
cxx         vs(j)=.51480*vs(j)
cxx         cpa(j)=1.8520*cpa(j)
cxx         if(iitype(j) .eq. 2) phid(j)=1.8520*phid(j)
cxx         x2(j)=1.8520*x2(j)
cxx         y2(j)=1.8520*y2(j)
cxx      endif
c: For type=2 and vs not zero, set t2 according to dist traveled and velocity:
         if(vs(j) .ne. 0 .and. iitype(j) .eq. 2) then
            t2(j)=t1(j) + 1000.*sqrt((x2(j)-cpa(j))**2 +
     .         (y2(j)-phid(j))**2)/(vs(j)*60.)
         endif
         if(dt(j) .lt. 0.) then
            nrleg(j)=iabs(nint(dt(j)))
            if(nrleg(j) - abs(dt(j)) .ne. 0.) then
               print *,'NPT on source track must be integer: ',j,dt(j)
               stop
            endif
         else
            nrleg(j)=nint((t2(j) - t1(j))/dt(j)) + 1
         endif
         if((iicont(j) .eq. 0 .and. nrleg(j) .lt. 1) .or.
     .      (iicont(j) .eq. 1 .and. nrleg(j) .lt. 2)) then
            print *,'NPT for iicont=0 leg must be > 0, ',
     .         'NPT for iicont=1 leg must be > 1: ',j,dt(j),nrleg(j)
            stop
         endif
         nsrc=nsrc + nrleg(j)
         if(iicont(j) .eq. 1) nsrc=nsrc - 1
      enddo
      iicont(1)=0
      call mem_lim(nsrc,NRNGMAX,MLINE,LML,'nsrc',4,'NRNGMAX',7,1,1)
c
c: Compute source track positions:
      r0=0.
      tlast=0.
      nsrc=0
      do j=1,nseg
         if(iitype(j) .eq. 1) then
            cosphi=cos(phid(j)*pierad)
            sinphi=sin(phid(j)*pierad)
c: for vs=0, times are distances; for vs not 0, convert to min to km
            fac=1.
            if(vs(j) .ne. 0.) fac=vs(j)*60./1000.
            if(iicont(j) .eq. 0) then
c: xxx0,yyy0 are (x,y) source coordinates at cpa:
c: Change phid to be in deg E of N, rather than from x- to y- axis,
c: +cpa means a clockwise track, -cpa means a ccw track:
               xxx0=-cpa(j)*cosphi
               yyy0=cpa(j)*sinphi
               x1x=xxx0 + fac*t1(j)*sinphi
               y1x=yyy0 + fac*t1(j)*cosphi
               x2x=xxx0 + fac*t2(j)*sinphi
               y2x=yyy0 + fac*t2(j)*cosphi
               time0=t1(j)
            else
c: begin at kr=2 for continuous legs so there is no overlap of points:
               x1x=x2x
               y1x=y2x
c: Change phid to be in deg E of N, rather than from x- to y- axis:
               x2x=x1x + fac*(t2(j)-t1(j))*sinphi
               y2x=y1x + fac*(t2(j)-t1(j))*cosphi
               time0=tlast
            endif
         else
            if(iicont(j) .eq. 0) then
               x1x=cpa(j)
               y1x=phid(j)
               time0=t1(j)
            else
               x1x=x2x
               y1x=y2x
               time0=tlast
            endif
            x2x=x2(j)
            y2x=y2(j)
         endif
         nrden=max(1,nrleg(j) - 1)
         delx=(x2x - x1x)/nrden
         dely=(y2x - y1x)/nrden
         delr=sqrt(delx**2 + dely**2)
         delt=(t2(j) - t1(j))/nrden
         do kr=iicont(j)+1,nrleg(j)
            nsrc=nsrc + 1
            xsrc(nsrc)=1000.*(x1x + (kr-1)*delx)
            ysrc(nsrc)=1000.*(y1x + (kr-1)*dely)
c: compute source track distance:
cxx         rangx(nsrc)=r0 + 1000.*(kr-1)*delr
c: t_src is the current time (or distance along track for vs=0):
            t_src(nsrc)=time0 + (kr-1)*delt
         enddo
         tlast=t_src(nsrc)
cxx      r0=rangx(nsrc)
      enddo
cxx   write(2,218)
218   format('### SOURCE TRACK')
cxx   do j=1,nsrc
cxx      write(2,220) .001*xsrc(j),.001*ysrc(j)
220      format(f8.3,2x,f8.3)
cxx   enddo
c
      return
110   print *,'Error opening or reading option file ',opt_file(1:lopt)
      print *,'*() line number = ',nline
      stop
c
      end
ccc
      subroutine arr_geom(iierr)
c
c: Reads in the receiver array geometry.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 iierr,nline,jz,jy,jx,jrec,iiarr,nx,ny,nz,j,nxy,nrd
      real*4 zrj,zr1,yr,xr,dx,dy,dz
      character*64 eline
      data eline/'INVALID INPUT IN ARRAY FILE: '/
c
      write(2,95)
95    format(/'### RECEIVER ARRAY GEOMETRY FILE INFORMATION ###')
      open(62,file=geom_file(1:lgeom),form='formatted',
     .   status='old',err=99)
      nline=0
      call star2(62,nline,2,geom_file,1,iiwrite)
      read(62,*) iiarr
      call check_val_i(iiarr,1,2,eline,29,nline,'iiarr',5,iierr)
      nrec=0
      if(iiarr .eq. 1) then
         call star2(62,nline,2,geom_file,1,iiwrite)
         read(62,*) zr1,nx,ny,nz,dx,dy,dz
         do 110 jz=1,nz
            zrj=zr1 + (jz-1)*dz
            do 110 jy=1,ny
               yr=(jy-1)*dy
               do 110 jx=1,nx
                  xr=(jx-1)*dx
                  nrec=nrec + 1
                  xrec(nrec)=xr
                  yrec(nrec)=yr
                  zrec(nrec)=zrj
110      continue
         call star2(62,nline,2,geom_file,0,iiwrite)
         call star2(62,nline,2,geom_file,0,iiwrite)
      else
         call star2(62,nline,2,geom_file,0,iiwrite)
         call star2(62,nline,2,geom_file,1,iiwrite)
         read(62,*) nz,(zsr(j),j=1,nrd(nz))
         call uni_space(nz,zsr,1.e0)
         call star2(62,nline,2,geom_file,0,iiwrite)
         do 112 jz=1,iabs(nz)
            read(62,*) nxy,(xrec(j),yrec(j),j=nrec+1,nrec+nxy)
            write(2,206) nxy,(xrec(j),yrec(j),j=nrec+1,nrec+1)
206         format(i8,f8.1,',',f8.1)
            do j=nrec+2,nrec+nxy
               write(2,207) xrec(j),yrec(j)
207            format(8x,f8.1,',',f8.1)
            enddo
            do jrec=nrec+1,nrec+nxy
               zrec(jrec)=zsr(jz)
            enddo
            nrec=nrec + nxy
112      continue
      endif
c: Receiver label for TL plots for source track option:
      do j=1,nrec
         rec_lab(j)=j
      enddo
      close(62)
c
      return
99    print *,'Error opening array geometry file ',geom_file
      stop
      end
ccc
      subroutine star2(nf,nline,nfout,fname,necho,iiwrite)
c
c: Searches for the next line of an input file that begins
c: with the key symbol '*', which indicates that the next
c: line contains data.
c
      implicit none
      integer*4 nf,nline,nfout,necho,iiwrite,j,lchfs
      character*1 chstar
      character*64 fname
      character*80 chfs
      data chstar/'*'/
c
      nline=nline + 1
10    read(nf,100,end=400,err=500) chfs(1:80)
100   format(a80)
      if(chfs(1:1) .ne. chstar) goto 10
      call lname80(chfs,lchfs)
      if(iiwrite .ne. 0) then
         write(nfout,110) chfs(1:lchfs)
110      format(a)
         do j=1,necho
            read(nf,102,end=390,err=390) chfs(1:80)
102         format(a80)
            call lname80(chfs,lchfs)
            write(nfout,110) chfs(1:lchfs)
         enddo
         do j=1,necho
            backspace(nf)
         enddo
      endif
c
390   return
400   continue
      print *,'search for * reached end of file for line # ',nline,
     .   '; file ',fname
      stop
500   continue
      print *,'search for * encountered end of file or error for '//
     .   'line # ',nline,' in input file ',fname
      stop
      end
ccc
      subroutine check_val_i(val,val_lo,val_hi,eline,le,nline,
     .   vname,lv,iierr)
c
c: Checks integer input val to see if it is in the range of allowable 
c: values, val_lo to val_hi.
c
      implicit none
      integer*4 val,val_lo,val_hi,nline,le,lv,iierr
      character*64 eline,vname
c
      if(val .lt. val_lo .or. val .gt. val_hi) then
         iierr=1
         print *,' '
         print *,eline(1:le),' LINE # ',nline
         print *,'VARIABLE NAME = ',vname(1:lv),'; VALID RANGE = ',
     .      val_lo,' ,',val_hi
         print *,'ENTERED VALUE = ',val
      endif
c
      return
      end
ccc
      subroutine check_val_r4(val,val_lo,val_hi,eline,le,nline,
     .   vname,lv,iierr)
c
c: Checks real*8 input val to see if it is in the range of allowable 
c: values, val_lo to val_hi.
c
      implicit none
      real*4 val,val_lo,val_hi
      integer*4 nline,le,lv,iierr
      character*64 eline,vname
c
      if(val .lt. val_lo .or. val .gt. val_hi) then
         iierr=1
         print *,' '
         print *,eline(1:le),' LINE # ',nline
         print *,'VARIABLE NAME = ',vname(1:lv),'; VALID RANGE = ',
     .      val_lo,' ,',val_hi
         print *,'ENTERED VALUE = ',val
      endif
c
      return
      end
ccc
      subroutine check_val_r8(val,val_lo,val_hi,eline,le,nline,
     .   vname,lv,iierr)
c
c: Checks real*8 input val to see if it is in the range of allowable 
c: values, val_lo to val_hi.
c
      implicit none
      real*8 val,val_lo,val_hi
      integer*4 nline,le,lv,iierr
      character*64 eline,vname
c
      if(val .lt. val_lo .or. val .gt. val_hi) then
         iierr=1
         print *,' '
         print *,eline(1:le),' LINE # ',nline
         print *,'VARIABLE NAME = ',vname(1:lv),'; VALID RANGE = ',
     .      val_lo,' ,',val_hi
         print *,'ENTERED VALUE = ',val
      endif
c
      return
      end
ccc
      function nrd(n)
c
      implicit none
      integer*4 nrd,n
c
      if(n .ge. 0) then
         nrd=n
      else
         nrd=2
      endif
c
      return
      end
