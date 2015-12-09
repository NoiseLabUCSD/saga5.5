      subroutine pw_refco(ir,jrun,nrunx)
c
      use parms_com
      use i_o_com
      use gen_com
c
      implicit none
      integer*4 ir,jrun,nrunx,jang,jf,j,nffth1,iibad,lfft
      complex*8 zzero
      complex*16 rc1(3)
      data zzero/(0.,0.)/
c
      iibad=0
      call mem_lim(nang,NSRMAX,MLINE,LML,'nang',4,'NSRMAX',6,iibad,0)
      call mem_lim(nfreq,NFCWMAX,MLINE,LML,'nfreq',5,'NFCWMAX',7,
     .   iibad,0)
cc    call mem_lim(nang*nfreq,NTLMAX,MLINE,LML,'nang*nfreq',10,
cc   .   'NTLMAX',6,iibad,1)
      if(iibad .eq. 1) stop
c
c
c: Don't allow sheet changes at branch cuts:
      iich=0
      iich_ref=0
c: Loop over frequencies:
      allocate(tlc(nfreq,nang))
      allocate(r4mat1(nfreq,nang))
      allocate(r4mat2(nfreq,nang))
      do jf=1,nfreq
         f_hz=fcw(jf)
         call freq_init
c: Make reference layer at bottom of ocean so that BL is done there:
         kduct=1
         jduct(1,kduct)=jobot
         jduct(2,kduct)=2
         call zref_chng
         xkh=xkref - dcmplx(1.d-7*dreal(xkref),0.d0)
         call sheet_init(xkh,1,iish,iish_ref)
c: Loop over angles:
         do jang=1,nang
            xkh=xkref*dcos(theta(jang)*pie/180.)
            call xkh_init(1)
            call rp_calc(1,rc1,1,0)
            tlc(jf,jang)=rc1(1)
            r4mat1(jf,jang)=20.*dlog10(abs(rc1(1)))
            r4mat2(jf,jang)=datan2(dimag(rc1(1)),dreal(rc1(1)))
     .         *180./pie
         enddo
      enddo
c
      if(iirc .eq. 1) then
cpg         call hdf_write_gen(2+ir,outroot,lout,fcw,nfreq,theta,nang,
cpg     .      var_ax,nrunx,theta,1,theta,1,r4mat1,1,nfreq,1,nang,
cpg     .      jrun,jrun,1,1,1,1,'Frequency - Hz',14,'Angle - deg',11,
cpg     .      'Parameter',9,' ',1,' ',1,'Reflection Coeff - dB',31,4,1,2)
cpg        call hdf_write_gen(2+ir,outroot,lout,fcw,nfreq,theta,nang,
cpg    .      var_ax,nrunx,theta,1,theta,1,r4mat2,1,nfreq,1,nang,
cpg     .      jrun,jrun,1,1,1,1,'Frequency - Hz',14,'Angle - deg',11,
cpg     .      'Parameter',9,' ',1,' ',1,'Refl Coeff Phase - deg',22,5,1,2)
      elseif(iirc .eq. 2) then
         xh(1)=0.
         xh(11)=nang
         call openfftout2(10,outroot,lout,fftfile,lfft,nfft,
     .      fsrc,freq1,freq2,xh,NRECL)
         do jang=1,nang
            xh(4)=jang
            write(10,rec=jang) (xh(j),j=1,20),(tlc(j,jang),j=1,nfreq)
         enddo
         close(10)
      endif
      deallocate(tlc)
      deallocate(r4mat1)
      deallocate(r4mat2)
c
      return
      end
