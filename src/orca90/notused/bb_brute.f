       subroutine bb_brute
c
c: Computes the broadband field from fmin to fmax in steps of df
c: using brute force calls of the usual CW routines
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jm,jfbb,nmbbtot,nm_fmax,nctot0,ncjf,jzjm
      real*4 ncperm
c
      call bb_init
      call zmx_init
      f_hz=faxbb(nfbb)
      call freq_init
      call bb_fft_init
c
      nmbbtot=0
      do jfbb=nfbb,1,-1
         nctot0=nctot
         f_hz=faxbb(jfbb)
         call freq_chng
         call mode_find(0)
         if(jfbb .eq. nfbb) then
            if(iiAih(1) .eq. -1 .and. iiAih(2) .eq. -1) then
               nm_fmax=nm_tot
            else
               nm_fmax=max(nm_tot+int(nm_tot/10),nm_tot + 12)
            endif
            if(iiout .ne. 0) call bb_out_init(nm_fmax)
cc          NM_NF_MAX=nm_fmax*nfbb
            allocate(knbb(nfbb,nm_fmax))
            allocate(eig_bb(5,nfbb,nm_fmax))
            allocate(iish_bb(nfbb,nm_fmax))
            if(iikn .gt. 0) then
               allocate(r4mat1(nfbb,nm_fmax))
               allocate(r4mat2(nfbb,nm_fmax))
               r4mat1=NaN
               r4mat2=NaN
            endif
            if(iimf .ne. 0) then
               allocate(r4mat3(nrec,nm_fmax))
               allocate(r4mat4(nrec,nm_fmax))
            endif
            if(iidc .ne. 0) then
               allocate(r4mat5(nfbb,nm_fmax))
               allocate(r4mat6(nfbb,nm_fmax))
               r4mat5=NaN
               r4mat6=NaN
            endif
            do jm=1,nm_fmax
               xmode(jm)=jm
            enddo
         else
            nm_tot=min(nm_fmax,nm_tot)
         endif
         nmbb(jfbb)=min(nm_fmax,nm_tot)
         do jm=1,nm_tot
            call bb_field(kn(jm),phi,dphi,dpsi,exp_gbs,jm,tf,jfbb,jm)
            knbb(jfbb,jm)=kn(jm)
            eig_bb(1,jfbb,jm)=eig_char(1,jm)
            eig_bb(2,jfbb,jm)=eig_char(2,jm)
            eig_bb(3,jfbb,jm)=eig_char(3,jm)
            eig_bb(4,jfbb,jm)=eig_char(4,jm)
            eig_bb(5,jfbb,jm)=eig_char(5,jm)
            call iish_code(iish,iish_ref,iish_bb(jfbb,jm),1)
            if(iift .eq. 1) call bb_write(jm,jfbb,faxbb,kn(jm),kw0,iish)
         enddo
c
         if(iikn .ne. 0) then
            do jm=1,nm_tot
               r4mat1(jfbb,jm)=dreal(kn(jm))
               r4mat2(jfbb,jm)=dimag(kn(jm))
               if(iikn .gt. 1) r4mat2(jfbb,jm)=r4mat2(jfbb,jm)/kw0
            enddo
         endif
         if(iimf .ne. 0) then
            do jrec=1,nrec
               jz=mzrec(jrec)
               do jm=1,nm_tot
                  jj=(jm-1)*nzsr + jz
                  r4mat3(jrec,jm)=abs(phi(jz,jm))
                  r4mat4(jrec,jm)=atan2(aimag(phi(jz,jm)),
     .               real(phi(jz,jm)))/piedeg
               enddo
               do jm=nm_tot+1,nm_fmax
                  r4mat3(jrec,jm)=NaN
                  r4mat4(jrec,jm)=NaN
               enddo
            enddo
            call hdf_write_gen(3,outroot,lout,zrec,nrec,xmode,
     .         nm_fmax,faxbb,nfbb,1.,1,1.,1,r4mat3,1,nrec,
     .         1,nm_fmax,jfbb,jfbb,1,1,1,1,'Depth - m',9,'Mode',4,
     .         'Frequency - Hz',14,' ',1,' ',1,'mag(phi)',8,12,1,2)
            call hdf_write_gen(3,outroot,lout,zrec,nrec,xmode,
     .         nm_fmax,faxbb,nfbb,1.,1,1.,1,r4mat4,1,nrec,
     .         1,nm_fmax,jfbb,jfbb,1,1,1,1,'Depth - m',9,'Mode',4,
     .         'Frequency - Hz',14,' ',1,' ',1,'arg(phi)',8,13,1,2)
         endif
         if(iimf .eq. 4) then
            do jrec=1,nrec
               jz=mzrec(jrec)
               do jm=1,nm_tot
                  jj=(jm-1)*nzsr + jz
                  r4mat3(jrec,jm)=real(dphi(jz,jm))
                  r4mat4(jrec,jm)=aimag(dphi(jz,jm))
               enddo
               do jm=nm_tot+1,nm_fmax
                  r4mat3(jrec,jm)=NaN
                  r4mat4(jrec,jm)=NaN
               enddo
            enddo
            call hdf_write_gen(3,outroot,lout,zrec,nrec,xmode,
     .         nm_fmax,faxbb,nfbb,1.,1,1.,1,r4mat3,1,nrec,
     .         1,nm_tot,jfbb,jfbb,1,1,1,1,'Depth - m',9,'Mode',4,
     .         'Frequency - Hz',14,' ',1,' ',1,'Re(dphi)',8,28,1,2)
            call hdf_write_gen(3,outroot,lout,zrec,nrec,xmode,
     .         nm_fmax,faxbb,nfbb,1.,1,1.,1,r4mat4,1,nrec,
     .         1,nm_tot,jfbb,jfbb,1,1,1,1,'Depth - m',9,'Mode',4,
     .         'Frequency - Hz',14,' ',1,' ',1,'Im(dphi)',8,29,1,2)
         endif
         if(iidc .ne. 0) then
            do jm=1,nm_tot
               r4mat5(jfbb,jm)=real(eig_char(4,jm))
               r4mat6(jfbb,jm)=w/dreal(kn(jm))
            enddo
         endif
c
         if(iiout .ne. 0) then
            if (nm_tot .gt. nm_fmax) print *,'nm_tot>nm_fmax',
     .         f_hz,nm_tot,nm_fmax
            write(33,rec=nh_off+jfbb) (kn(jm),jm=1,nmbb(jfbb)),
     .         ((phi(jz,jm),jz=1,nzsr),jm=1,nmbb(jfbb))
         endif
         nmbbtot=nmbbtot + nm_tot
         ncjf=nctot - nctot0
         print *,'Done f,nm_tot,#R1R2/mode = ',sngl(f_hz),nm_tot,
     .      float(ncjf)/max(1,nm_tot)
      enddo
c
      if(iikn .gt. 0.) then
         call hdf_write_gen(2,outroot,lout,faxbb,nfbb,xmode,
     .      nm_fmax,1.,1,1.,1,1.,1,r4mat1,1,nfbb,1,nm_fmax,
     .      1,1,1,1,1,1,'Frequency - Hz',14,'Mode',4,
     .      ' ',1,' ',1,' ',1,'Re(kn)',6,10,1,2)
         call hdf_write_gen(2,outroot,lout,faxbb,nfbb,xmode,
     .      nm_fmax,1.,1,1.,1,1.,1,r4mat2,1,nfbb,1,nm_fmax,
     .      1,1,1,1,1,1,'Frequency - Hz',14,'Mode',4,
     .      ' ',1,' ',1,' ',1,'Im(kn)',6,11,1,2)
         deallocate(r4mat1)
         deallocate(r4mat2)
      endif
      if(iimf .gt. 0) then
         deallocate(r4mat3)
         deallocate(r4mat4)
      endif
c
      if(iidc .eq. 1 .or. iidc .eq. 3) then
         call hdf_write_gen(2,outroot,lout,faxbb,nfbb,xmode,
     .      nm_fmax,1.,1,1.,1,1.,1,r4mat5,1,nfbb,1,nm_fmax,
     .      1,1,1,1,1,1,'Frequency - Hz',14,'Mode',4,
     .      ' ',1,' ',1,' ',1,'Group Speed - m/s',
     .      17,24,1,2)
      endif
      if(iidc .eq. 2 .or. iidc .eq. 3) then
         call hdf_write_gen(2,outroot,lout,faxbb,nfbb,xmode,
     .      nm_fmax,1.,1,1.,1,1.,1,r4mat6,1,nfbb,1,nm_fmax,
     .      1,1,1,1,1,1,'Frequency - Hz',14,'Mode',4,
     .      ' ',1,' ',1,' ',1,'Phase Speed - m/s',
     .      17,25,1,2)
      endif
      if(iidc .ne. 0) then
         deallocate(r4mat5)
         deallocate(r4mat6)
      endif
c
c: Close files:
      if(iimt .ne. 0) close(21)
      if(iift .ne. 0) close(14)
      if(iiout .ne. 0) then
         close(33)
         open(33,file=outroot(1:lout)//'_bbeig',access='direct',
     .      recl=NRECL*lheadw)
         write(33,rec=1) fsbb,nfftbb,fmin,fmax,nfbb,nzsr,nm_fmax,
     .      rmin,rmax,sngl(cfmin),
     .      (zsr(jsr),jsr=1,nzsr),(nmbb(jf),jf=1,nfbb),iirx,
     .      (sngl(rho_sr(jsr)),jsr=1,nzsr),(faxbb(jf),jf=1,nfbb)
cc       print *,'nmbb at end= ',(nmbb(jf),jf=1,nfbb)
cc       print *,'lheadw,lrecw = ',lheadw,lrecw,nh_off
         close(33)
      endif
c
c: Output FFT file:
      call bb_fft_out
c
      ncperm=float(nctot)/float(max(1,nmbbtot))
      write(2,120) nmbbtot,nctot,ncperm
120   format('CW LOOP # MODES = ',i8,'; # R1R2 CALCS = ',i8,
     .   '; #CALCS/MODE = ',f5.2)
c
      return
      end
ccc
      subroutine bb_write(jmo,jfbb,faxbb,k,kw0,iish)
c
      implicit none
      integer*4 jmo,jfbb,iish(2,2)
      real*4 faxbb(jfbb)
      real*8 kw0,twpie
      complex*16 k
      data twpie/6.28318530717959/
c
ccc   write(14,100) jmo,jfbb,faxbb(jfbb),dreal(k)/kw0,dimag(k)/kw0,
ccc  .   iish(1,1),iish(2,1),iish(1,2),iish(2,2)
c: Output IM(k) in terms of dB/km:
      write(14,100) jmo,jfbb,faxbb(jfbb),dreal(k)/kw0,8685.9*dimag(k),
     .   iish(1,1),iish(2,1),iish(1,2),iish(2,2)
100   format(i3,1x,i4,1x,f8.3,1x,e14.8,1x,e14.8,4(1x,i2))
c
      return
      end
ccc
      subroutine bb_fft_out
c
      use parms_com
      use i_o_com
      use gen_com
      implicit none
      integer*4 irec,jsrc,jrec,j,lfft
      real*8 vgref,tbuf,tstart
c
      if(iifft .ne. 0) then
c: Open output FFT file:
         call openfftout2(10,outroot,lout,fftfile,lfft,nfftbb,
     .      fsbb,fmin,fmax,xhbb,NRECL)
cc       open(10,file=outroot(1:lout)//'_fft',access='direct',
cc   .      recl=NRECL*(22+nfftbb))
         irec=0
c: Set up FFT header:
c: Flag that FFT is a simulation, not data:
c         xhbb(1)=-1.
         xhbb(2)=nfftbb
         xhbb(3)=1
         xhbb(10)=0.
         xhbb(11)=nrec
         xhbb(12)=1
c: Make ref time for impulse resp 10% of the way through total time window:
         tbuf=.10*Tw
c: Output transfer functions to the FFT file:
c: Make reference sound speed, used to shift impulse responses as a 
c: function of range, the minimum speed in the profile.
         vgref=cfmin
         do jsrc=1,nsrc
            tstart=rng_sr(jsrc,1)/vgref - tbuf
c: Try letting MOY be km so that Time in min is really km:
            xhbb(1)=rng_sr(jsrc,1)/1000.
c: tf(nfbb,nrec,nrng):
            call bb_align(tf(1,1,jsrc),tstart,wbb,nfbb,nrec,iifull)
            xhbb(14)=tstart
            do jrec=1,nrec
               xhbb(4)=jsrc
c: Not compatible with new FFT format:
cc             xhbb(15)=zsr(mzrec(jrec))
cc             xhbb(17)=cp_sr(mzrec(jrec))
               xhbb(18)=rng_sr(jsrc,jrec)
               xhbb(19)=zsr(mzrec(jrec))
               irec=irec + 1
c: Note that CONJUGATE of transfer function is output to FFT file:
               write(10,rec=irec) (xhbb(j),j=1,20),
     .            (conjg(tf(j,jrec,jsrc)),j=1,nfbb)
            enddo
         enddo
         close(10)
      endif
c
      return
      end
ccc
      subroutine bb_fft_init
c
      use parms_com
      use i_o_com
      use gen_com
c
c: Open file for the frequency trajectories of the modes:
      if(iift .eq. 1) then
         open(14,file=outroot(1:lout)//'_ftraj',form='formatted')
         write(14,100) -1,faxbb(1),faxbb(nfbb),fsbb,nfftbb,kw0,0,0,0
100      format(i3,2x,3(f8.3,2x),i5,2x,e10.4,3(1x,i2))
c: Write branch point(s) to frequency trajectory file:
         call bb_write(0,nfbb,faxbb,xkbp(1,1),kw0,iish)
         call bb_write(0,nfbb,faxbb,xkbp(1,2),kw0,iish)
         if(geo(2,1,1) .gt. 500.) then
            call bb_write(0,nfbb,faxbb,xkbp(2,1),kw0,iish)
            call bb_write(0,nfbb,faxbb,xkbp(2,2),kw0,iish)
         endif
      endif
c
      return
      end
ccc
      subroutine bb_align(tf,tstart,wbb,nfbb,nrec,iifull)
c
      implicit none
      integer*4 nfbb,nrec,iifull,jzr,jf
      complex*8 tf(nfbb,nrec),eiph
      real*8 tstart,phase,wbb(nfbb)
c
      do jf=1,nfbb
         phase=-wbb(jf)*tstart
         eiph=cmplx(cos(phase),sin(phase))
         do jzr=1,nrec
            tf(jf,jzr)=tf(jf,jzr)*eiph
         enddo
      enddo
c: If highest frequency is Nyquist, make sure imaginary part is zero:
      if(iifull .eq. 1) then
         do jzr=1,nrec
            tf(nfbb,jzr)=cmplx(real(tf(nfbb,jzr)),0.e0)
         enddo
      endif
c
      return
      end
