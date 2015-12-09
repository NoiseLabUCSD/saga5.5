      subroutine bb_init
c
c: Checks and initializes variables for broadband mode calculations.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 n,nfft0,j,iibad,ii1st,jr,jzs,j1,j2,j3
      real*8 rng_im
      complex*16 cfac
      data cfac/(1.77245385090552d0,1.77245385090552d0)/
      data ii1st/1/
c
      if(ii1st .eq. 0) goto 10
c
c: Initialize source and receiver geometry arrays:
      call sr_geom(1)
      if(iigbs .eq. 0) then
         do jr=1,nrng
            sq2pir(jr)=cfac/dsqrt(range(jr))
         enddo
      else
         jzs=1
         rng_im=-b_gbs(jzs)*cos(th_gbs(jzs)*pie/180.d0)
         do jr=1,nrng
            sq2pir(jr)=cfac/cdsqrt(dcmplx(range(jr),rng_im))
         enddo
      endif
c
      if(Tw .gt. 0.e0) then
         nfftbb=nint(Tw)
         Tw=nfftbb/fsbb
      else
         nfftbb=-Tw*fsbb
      endif
      n=int(log10(float(nfftbb))/log10(2.) + .9)
      nfft0=nfftbb
      nfftbb=2**n
      nffth1bb=nfftbb/2 + 1
      Tw=float(nfftbb)/fsbb
      if(nfftbb .ne. nfft0) then
         print *,'NFFT set to a power of 2: ',nfftbb,'; Tw = ',Tw
         write(2,*) 'NFFT set to a power of 2: ',nfftbb,'; Tw = ',Tw
      endif
      df=dble(fsbb)/dfloat(nfftbb)
      iifull=0
      if(fmaxorca .gt. 0.d0) then
         nf1=max(1,nint(fminorca/df)) + 1
         nf2=nint(fmaxorca/df) + 1
         fminorca=(nf1 - 1)*df
         fmaxorca=(nf2 - 1)*df
         nfbb=nf2 - nf1 + 1
         if(nf2 .eq. nffth1bb) iifull=1
      else
         nfbb=nfcw
      endif
      iibad=0
      allocate(faxbb(nfbb))
      allocate(wbb(nfbb))
      allocate(nmbb(nfbb))
      allocate(kim_bb(nfbb))
c
      allocate(tf(nfbb,nrec,nrng))
      allocate(phibb(nfbb,nzsr))
      allocate(dpsibb(nfbb,nzsr))
      if(iibad .eq. 1) stop
c
10    continue
c
      do j=1,nfbb
         nmbb(j)=0
      enddo
c
c      do j=1,nfbb*nrec*nsrc
c         tf(j)=(0.,0.)
c      enddo
	do j1=1,nfbb
	do j2=1,nrec
	do j3=1,nrng
		tf(j1,j2,j3)=(0.,0.)
	enddo
	enddo
	enddo
      if(fmaxorca .gt. 0.d0) then
         do j=1,nfbb
            faxbb(j)=fminorca + (j-1)*df
            wbb(j)=twpie*faxbb(j)
            kim_bb(j)=1.d100
cc          phim_bb(j)=0.
         enddo
      else
         do j=1,nfbb
            fmaxorca=max(fmaxorca,fcw(j))
            faxbb(j)=fcw(j)
            wbb(j)=twpie*faxbb(j)
            kim_bb(j)=1.d100
cc          phim_bb(j)=0.
         enddo
      endif
c
      if(iimt .ne. 0) then
         open(21,file=outroot(1:lout)//'_mtraj',form='formatted')
      endif
c
      ii1st=0
c
      return
      end
ccc
      subroutine bb_out_init(nm_fmax)
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 nm_fmax
c
c: Open direct access for broadband mode characteristics, if desired:
cc    if(iirx .eq. -1) then
cc       lrecw=(4+4*nzsr)*nfbb
cc    elseif(iirx .eq. 0) then
c: (Complex*16 eigenvalue + complex*8 mode function at nzsr depth) at
c: each frequency bin:
cc       lrecw=4*nm_fmax + 2*nzsr*nm_fmax
cc    else
c: (Complex*16 eigenvalue + real*4 mode function at nzsr depth) at
c: each frequency bin:
cc       lrecw=(4+nzsr)*nm_fmax
cc    endif
c: 7-17-99: Simplify by outputting complex*8 mode functions in all cases:
      lrecw=4*nm_fmax + 2*nzsr*nm_fmax
      open(33,file=outroot(1:lout)//'_bbeig',access='direct',
     .   recl=NRECL*lrecw)
      lheadw=11 + 3*nzsr + 2*nfbb
      nh_off=(lheadw - 1)/lrecw + 1
c
      return
      end
