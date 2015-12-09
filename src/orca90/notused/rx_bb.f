      subroutine rx_bb
c
c: Does broadband calculations for modes on the real axis.
c
      use parms_com
      use i_o_com
      use gen_com
c
      implicit none
      real*8 kmax,kmin,f1,f2,fint_max
      integer*4, dimension(:), allocatable :: nct_f
      integer*4 ndv,nm1,nm2,jflo,jfhi,ndfmax,nm_fmax,nmex,nmint,nmtot,
     .   ncx,nback,nbackm,nminok,iibad,nmlo,ndf,
     .   j1,j2,jmlo,jmhi,iick,jsr,jf,jfbb,jm,jrec,jz
c
      iiwrite=0
      fint_max=15.d0
      iick=0
      if(iidiag .eq. -10) iick=1
      ndv=6
      nback=0
      nbackm=0
      nminok=0
c: Check input variables, sizes of arrays, and initialize variables:
      call bb_init
      allocate(nct_f(nfbb))
c: Make max frequency jump close to fint_max:
      ndfmax=min(nfbb,nint(fint_max/df)+1)
      if(ndfmax .gt. 200) then
         print *,'Increase fm array size in rx_bb_mterp: ',200,ndfmax
         stop
      endif
c
      jflo=nfbb
      f_hz=faxbb(jflo)
      call rx_prep
      call rx_zmx_init
      nctot=0
      nclast=0
      nzoff=0
c
      call rx_freq_init(kmax,kmin)
      nmode=0
      ncx=nctot
      kim_min=1.d100
c: Find modes at highest frequency in whole band:
      call rx_mode_int(kmax,kmin,nm_lim,ndv,0,1,1,0)
      nct_f(nfbb)=nctot - ncx
      nmlo=nmode
      nm_fmax=nmode
      do jm=1,nmode
         xmode(jm)=jm
      enddo
c
      j1=nfbb*nm_fmax/2 + 2
      j2=j1 + j1
      allocate(knbb(nfbb,nm_fmax))
      knbb=NaN
      allocate(vg_bb(nfbb,nm_fmax))
      vg_bb=NaN
      allocate(Dvg_w(nfbb,nm_fmax))
      allocate(phi_bb(4,nzsr,nfbb,nm_fmax))
      allocate(iish_bb(nfbb,nm_fmax))
      allocate(pmat(4,nzsr))
      allocate(phi_f(nzsr))
      if(iimf .ne. 0) then
         allocate(r4mat_3d(nrec,nm_fmax,nfbb))
      endif
c
c: Enter results in lower f of block (to be shifted to upper in loop):
      call rx_bb_enter(1,nmode,jflo,1,ndfmax,phi,dphi,
     .   Dphiz_w,Ddphiz_w,nmex)
c
c: Open mode characteristics file if desired:
      if(iiout .ne. 0) then
         call bb_out_init(nm_fmax)
         call rx_bb_out(ndfmax,nm_fmax,jflo,jflo+1,nh_off)
      endif
c
      do while(jflo .ne. 1)
         jfhi=jflo
         nm2=nmlo
         jflo=max(jfhi - ndfmax + 1,1)
         f1=faxbb(jflo)
         f2=faxbb(jfhi)
c: Transfer mode quantities from former lower f to current upper f in block:
         ndf=jfhi - jflo + 1
cc         print *,'block shift: '
         call rx_block_shift(nzsr,ndfmax,nmlo,1,ndf,phi_bb,NaN)
c
c: Find modes at lower frequency in block:
         f_hz=f1
         call rx_freq_init(kmax,kmin)
         ncx=nctot
         nmode=0
         kim_min=1.d100
c: Don't allow to find more modes than at f2:
         call rx_mode_int(kmax,kmin,nm2,ndv,0,1,1,0)
         nct_f(jflo)=nctot - ncx
c
         nm1=nmode
         nmlo=nmode
cc         print *,'enter: '
         call rx_bb_enter(1,nmode,jflo,1,ndfmax,phi,dphi,
     .      Dphiz_w,Ddphiz_w,nmex)
c
         call rx_mcross(nm1,nm2,jflo,jfhi,jmlo,jmhi,ndrx,f1,f2,
     .      ccr_lo,ccr_hi,pie)
c
         call rx_bb_mterp(ndfmax,nm1,nm2,jflo,jfhi,jmlo,jmhi,
     .      nmex,nmint,nct_f,ndv)
c
         if(iick .eq. 1) then
            call interp_check(ndfmax,nm2,phi,jflo,jfhi)
         endif
c: Output mode characteristics if desired:
         if(iiout .ne. 0) then
            call rx_bb_out(ndfmax,nm_fmax,jflo,jfhi,nh_off)
         endif
         if(iimf .ne. 0) then
            do jfbb=jflo,jfhi
               jf=jfbb-jflo+1
               do jm=1,nm_fmax
                  do jrec=1,nrec
                     jz=mzrec(jrec)
                     r4mat_3d(jrec,jm,jfbb)=phi_bb(1,jz,jf,jm)
                  enddo
               enddo
            enddo
         endif
      enddo
c
      iiwrite=1
cc    print *,'nct_f = ',(nct_f(ncx),ncx=1,nfbb)
cc    print *,'nback,nbackm,nminok = ',nback,nbackm,nminok
c
      if(nzoff .gt. 0) then
         print *,'Info msg: # of zeros of mode function off by one ',
     .      nzoff,' times.  Can happen with multiple duct environments.'
         write(2,*) 'Info msg: # of zeros of mode function off by one ',
     .      nzoff,' times.  Can happen with multiple duct environments.'
      endif
c
      nmtot=nmex + nmint
      write(2,120) nmtot,nmex,nmint,nctot,nctot/float(nmtot)
120   format('TOTAL # MODES = ',i8,'; # EXACT =',i8,'; # INTERP =',i8/
     .   '   # R1R2 CALCS = ',i8,'; # CALCS/MODE = ',f5.2)
c
      if(iimf .ne. 0) then
         call hdf_write_gen(3,outroot,lout,zrec,nrec,xmode,
     .      nm_fmax,faxbb,nfbb,1.,1,1.,1,r4mat_3d,1,nrec,1,nm_fmax,
     .      1,nfbb,1,1,1,1,'Depth - m',9,'Mode',4,
     .      'Frequency - Hz',14,' ',1,' ',1,'Re(phi)',7,12,1,2)
         deallocate(r4mat_3d)
      endif
      if(iidc .ne. 0) then
         allocate(r4mat5(nfbb,nm_fmax))
         allocate(r4mat6(nfbb,nm_fmax))
         do jfbb=1,nfbb
            do jm=1,nm_fmax
               r4mat5(jfbb,jm)=vg_bb(jfbb,jm)
               r4mat6(jfbb,jm)=wbb(jfbb)/dreal(knbb(jfbb,jm))
            enddo
         enddo
      endif
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
      if(iikn .ne. 0) then
         allocate(r4mat1(nfbb,nm_fmax))
         allocate(r4mat2(nfbb,nm_fmax))
         do jfbb=1,nfbb
            do jm=1,nm_fmax
               r4mat1(jfbb,jm)=real(knbb(jfbb,jm))
               r4mat2(jfbb,jm)=dimag(knbb(jfbb,jm))
               if(iikn .gt. 1) r4mat2(jfbb,jm)=r4mat2(jfbb,jm)/kw0
            enddo
         enddo
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
c
      if(iifft .ne. 0) call bb_fft_out
c
      if(iiout .ne. 0) then
         close(33)
         open(33,file=outroot(1:lout)//'_bbeig',access='direct',
     .      recl=NRECL*lheadw)
         write(33,rec=1) fsbb,nfftbb,fmin,fmax,nfbb,nzsr,nm_fmax,
     .      rmin,rmax,sngl(cfmin),
     .      (zsr(jsr),jsr=1,nzsr),(nmbb(jf),jf=1,nfbb),iirx,
     .      (sngl(rho_sr(jsr)),jsr=1,nzsr),(faxbb(jf),jf=1,nfbb)
         close(33)
      endif
c
      return
      end
ccc
      subroutine rx_bb_mterp(ndfmax,nm1,nm2,jflo,jfhi,jmlo,jmhi,
     .   nmex,nmint,nct_f,ndv)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 ndfmax,nm1,nm2,jf1,jflo,jfhi,jmlo,jmhi,nmex,nmint,ndv,
     .   ndvm,nct_f(nfbb),jm1,nm1x,nm2x,nvgbad,nvgok,iibad,n2k,nm0
      real*8 kmin
c
      nvgbad=0
      nvgok=0
      n2k=0
c
c: Find modes, if any, that are trapped so deep in deepest duct 
c: that they cannot cross:
      do jm1=1,jmlo-1
         call rx_bb_fint(ndfmax,jm1,jflo,jfhi,jflo,jfhi,
     .      nmex,nmint,nct_f,ndv,nvgok,nvgbad,n2k)
      enddo
c
      if(jmlo .le. jmhi) then
c: Find modes from jmlo to jmhi, which can cross due to multiple ducts, 
c: using safe method:
         iibad=0
c: Min number of modes that should be found: jmhi (when nm2>jmhi) or
c: nm1 (when all modes are in ducts, nm1<=jmhi):
         nm0=min(nm1,jmhi)
c: Only compute derivatives needed for dispersion curves:
         ndvm=min(ndv,3)
         if(iidc .eq. 0) ndvm=2
         do jf1=jflo+1,jfhi-1
c: Use found mode at slightly lower freq as kmin since k increases with f:
            kmin=knbb(jf1-1,jmhi)
cc      print *,'mint from mterp: ',jmlo,jmhi,jf1,jflo
            call rx_bb_mint(ndfmax,jmlo,jmhi,jf1,jflo,kmin,
     .         nmex,nct_f,ndvm)
c: FIX (4-2-96):
            if(nmode .lt. nm0) then
cc          if(nmode .ne. jmhi) then
c: This should never happen:
cc               print *,'rx_bb_mint did not find all 1st try: ',
cc     .            jf1,jflo,jfhi,f_hz,jmlo,jmhi,nm1,nm2,nmode
c: Try one more time, setting kmin to a lower value:
               if(jmhi .lt. nm1) then
                  kmin=knbb(jflo,jmhi+1)
               else
                  kmin=0.d0
               endif
      print *,'This should never happen mint: ',nmode+1,jmhi,jf1,jflo
               call rx_bb_mint(ndfmax,nmode+1,jmhi,jf1,jflo,kmin,
     .            nmex,nct_f,ndvm)
               if(nmode .lt. nm0) then
                  iibad=jf1
                  print *,'rx_bb_mint did not find all 1st try: ',
     .               jf1,jflo,jfhi,f_hz,jmlo,jmhi,nm1,nm2,nmode
               else
                  print *,'rx_bb_mint found 2nd try'
               endif
            endif
         enddo
         if(iibad .ne. 0) then
            do jf1=jflo,jfhi
               write(99,100) (dreal(knbb(jf1,jm1)),jm1=jmlo,jmhi)
100            format(100(e15.8,1x))
            enddo
            jf1=iibad
            kmin=knbb(jf1-1,jmhi)
cc      print *,'iibad mint call: ',iibad,jmlo,jmhi,jf1,jflo
            call rx_bb_mint(ndfmax,jmlo,jmhi,jf1,jflo,kmin,
     .         nmex,nct_f,ndv)
         endif
      endif
      do jm1=jmhi+1,nm1
         call rx_bb_fint(ndfmax,jm1,jflo,jfhi,jflo,jfhi,
     .      nmex,nmint,nct_f,ndv,nvgok,nvgbad,n2k)
      enddo
c
c: Do triangle from jflo to jfhi, where # modes increases from nm1 to nm2
      nm1x=nm1+1
      nm2x=min(nm1x+3,nm2)
      jf1=jflo+1
      kmin=0.d0
      do while(jf1 .lt. jfhi .and. nm1x .le. nm2x)
cc      print *,'mint triangle: ',nm1x,nm2,jf1,jflo
         call rx_bb_mint(ndfmax,nm1x,nm2,jf1,jflo,kmin,
     .      nmex,nct_f,ndv)
         if(nmode .ge. nm2x .and. jf1 .lt. jfhi-1) then
            do jm1=nm1x,nmode
               call rx_bb_fint(ndfmax,jm1,jf1,jfhi,jflo,jfhi,
     .            nmex,nmint,nct_f,ndv,nvgok,nvgbad,n2k)
            enddo
            nm1x=nmode+1
            nm2x=min(nm1x+3,nm2)
         endif
         jf1=jf1 + 1
      enddo
cc    print *,'nvgok,nvgbad,n2k final = ',jflo,jfhi,nvgok,nvgbad,n2k
      print *,'Finished frequency block ',faxbb(jflo),faxbb(jfhi),nmex,
     .   nmint
c
      return
      end
c
      subroutine rx_bb_fint(ndfmax,jm1,jf1x,jf2x,jflo,jfhi,
     .   nmex,nmint,nct_f,ndv,nvgok,nvgbad,n2k)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 ndfmax,jf1x,jf2x,jflo,jfhi,nmex,nmint,ndv,
     .   nct_f(nfbb),jm1,jf1,jf2,jfmid,jdf1,jdf2,jdf,ncx,
     .   ndf,jf,nstack,stack(2,50),ntry,jm,jlay_ref,ii_ref,nvgbad,
     .   nvgok,joff,nzero,j0,jzsr,n2k,iim0
      real*4 phi1,phi2,Dphi1_f,Dphi2_f,
     .   f1x,f2x,twpiex,delfx,att1,att2,delatt,attx
      real*8 f1,f2,k1,k2,k1p,k2p,
     .   k1pp,k2pp,kmin,kmax,klast,fm,fmsq,fmcu,kminx,kmid,ekmid(6),
     .   phmid,fac,dk_err,c1,c2,c3,c4,c5,c6,delf,
     .   vgmid,vgex,vg_err,vg,D2k_Df,D2k_err,kmidp,kmidpp
cc    real*8 err_w,err_k
      complex*16 knx
      data twpiex/6.28318530718/
c
cc	print *,'fint: jm1,jf1,jf2 = ',jm1,jf1x,jf2x
      jf1=jf1x
      jf2=jf2x
10    continue
c: Find mid frequency between f1 and f2:
      jfmid=(jf1 + jf2)/2
      f1=faxbb(jf1)
      f2=faxbb(jf2)
c: jdf's are frequency counters relative to current frequency block:
      jdf1=jf1 - jflo + 1
      jdf2=jf2 - jflo + 1
      jdf=jfmid - jflo + 1
      ndf=jf2 - jf1 + 1
      fac=dfloat(ndf-1)
      fm=(jdf - jdf1)/fac
c: Initialize frequency dependent mode eigenvalue variables:
      f_hz=faxbb(jfmid)
      call rx_freq_init(kmax,kminx)
c
      k1=dreal(knbb(jf1,jm1))
      k2=dreal(knbb(jf2,jm1))
      k1p=twpie/vg_bb(jf1,jm1)
      k2p=twpie/vg_bb(jf2,jm1)
      k1pp=Dvg_w(jf1,jm1)
      k2pp=Dvg_w(jf2,jm1)
cc    print *,'f1,f_hz,f2 = ',f1,f_hz,f2
cc    print *,'kp = ',k1p,(k2-k1)/(f2-f1),k2p
cc    print *,'kpp = ',k1pp,(k2p-k1p)/(f2-f1),k2pp
      call poly5_fit(f1,f2,k1,k2,k1p,k2p,k1pp,k2pp,
     .   c1,c2,c3,c4,c5,c6,delf)
      call poly5_eval(c1,c2,c3,c4,c5,c6,delf,fm,kmid,kmidp,kmidpp,
     .   fmsq,fmcu)
      call rx_ek_calc(kmid,ekmid,phmid,ndv)
      nct_f(jfmid)=nct_f(jfmid) + 1
      dk_err=ekmid(1)/ekmid(2)
      vgmid=twpie/kmidp
c
      if(abs(dk_err) .le. dk_max0) then
c: This routine no longer needed:
ccc      call rx_bb_norm_Dvg(kmid,ekmid,vgex,D2k_Df)
         call rx_bb_Dvg(g12(1,2,1),g22(1,2,1),xkh,xkhsq,xksq(1,nlay),
     .      w,wsq,twpie,vgex,D2k_Df)
cc       if(abs(vgex0-vgex)/vgex .gt. .001) print *,'vgex bad: ',
cc   .      vgex0,vgex
cc       if(abs(D2k_Df0-D2k_Df)/D2k_Df .gt. .001) print *,
cc   .      'D2k_Df bad: ',D2k_Df0,D2k_Df
         vg_err=dabs(vgex-vgmid)/vgex
         D2k_err=dabs((D2k_Df-kmidpp)/D2k_Df)
         if(vg_err .lt. 1.e-4 .and. D2k_err .ge. .1d0) n2k=n2k+1
         if(vg_err .lt. 1.d-4 .and. D2k_err .lt. .1d0) then
cc    print *,'D2k_err = ',jfmid,jm1,D2k_df,kmidpp,D2k_err
c: Fit good: interpolate from jf1 to jf2:
            nvgok=nvgok + 1
            f1x=sngl(f1)
            f2x=sngl(f2)
c: Do cubic fit of mode function at each depth:
            do jzsr=1,nzsr
               phi1=phi_bb(1,jzsr,jdf1,jm1)
               phi2=phi_bb(1,jzsr,jdf2,jm1)
               Dphi1_f=twpiex*phi_bb(3,jzsr,jdf1,jm1)
               Dphi2_f=twpiex*phi_bb(3,jzsr,jdf2,jm1)
               call cub_fit_r4(f1x,f2x,phi1,phi2,Dphi1_f,Dphi2_f,
     .            pmat(1,jzsr),pmat(2,jzsr),pmat(3,jzsr),pmat(4,jzsr),
     .            delfx)
            enddo
c: Interpolate attenuation linearly:
            att1=dimag(knbb(jf1,jm1))
            att2=dimag(knbb(jf2,jm1))
            delatt=att2 - att1
cc            print *,'fint interp: ',jf1+1,jf2-1
            do jf=jf1+1,jf2-1
               fm=(jf - jf1)/fac
               call poly5_eval(c1,c2,c3,c4,c5,c6,delf,fm,kmid,
     .            kmidp,kmidpp,fmsq,fmcu)
               attx=att1 + delatt*fm
               knx=dcmplx(kmid,dble(attx))
               knbb(jf,jm1)=knx
               vg_bb(jf,jm1)=twpie/kmidp
               jdf=jf-jflo + 1
               do jzsr=1,nzsr
                  phi_f(jzsr)=pmat(1,jzsr) + pmat(2,jzsr)*fm + 
     .               pmat(3,jzsr)*fmsq + pmat(4,jzsr)*fmcu
                  phi_bb(1,jzsr,jdf,jm1)=phi_f(jzsr)
               enddo
cc    print *,'k interp: ',jm1,jf,knbb(jf1,jm1),knbb(jf2,jm1),knx
cxx   print *,'phi interp: nzsr ',phi1,phi2,phi_f(nzsr)
cxx   print *,'phi interp: 1 ',phi_bb(1,jdf1,jm1),
cxx  .   phi_bb(1,jdf2,jm1),phi_f(1)
               call rx_bb_field(knx,phi_f,jm1,jf)
ccx            if(iift .eq. 5) write(88,108) faxbb(jf),jm1,1
108            format(f9.2,2x,i4,2x,i1)
            enddo
            nmint=nmint + (jf2-jf1-1)
c
            if(nstack .gt. 0) then
               jf1=stack(1,nstack)
               jf2=stack(2,nstack)
               nstack=nstack-1
               goto 10
            else
               goto 99
            endif
         else
            nvgbad=nvgbad + 1
         endif
      else
c: Compute error in group velocity for next comparison:
         vgex=-ekmid(2)/ekmid(3)
         vg_err=dabs(vgex-vgmid)/vgex
      endif
c
c: Fit bad: find mode exactly and enter into knbb,vg_bb,Dvg_w,phi_bb:
c: Allow more than 100* error???
      if(abs(dk_err) .le. 200.d0*dk_max0 .and. vg_err .lt. 1.d-2) then
c: Guess close enough to get eigenvalue by polishing root:
         ntry=0
35       if(abs(dk_err) .gt. dk_max0) then
            ntry=ntry+1
            if(ntry .gt. 2) then
               if(iidiag .ge. 2) print *,'polish ntry>2: ',ntry,jm1,
     .            jf1,jf2
            endif
c: If trouble finding mode, give up and use mode_int:
            if(ntry .gt. 4) goto 45
            kmid=kmid - dk_err
            call rx_ek_calc(kmid,ekmid,phmid,ndv)
            nct_f(jfmid)=nct_f(jfmid) + 1
            dk_err=ekmid(1)/ekmid(2)
            goto 35
         endif
c
         nmode=jm1
         call rx_enter(nmode,ncalc,nclast,nctot,kn,kmid,
     .      ekmid,ekn,phmid,mode_phz,w,iidiag)
         call rx_norm(jlay_ref,ii_ref,ndv,vg,D2k_Df)
         call rx_mf_lay_Dw(jlay_ref,ii_ref,vg,nzero)
         if(nzero .eq. nmode) then
            call rx_mode_fun_Dw(nmode,phi,dphi,Dphiz_w,Ddphiz_w,
     .         jlay_ref,ii_ref,vg)
c: Success using kmid guess.  Jump to rx_bb_enter:
            goto 50
         else
c: Found wrong mode, so use mode_int below:
c           print *,'Informative Message: nzero~=nmode in rx_bb_fint ',
c    .         nzero,nmode,f_hz
         endif
      endif
c
45    continue
      nmode=jm1-1
      joff=1
      kmin=max(kminx,kmid - 100.d0*dk_max0)
      if(nmode .gt. 0) then
         klast=dreal(knbb(jfmid,nmode))
         if(klast .lt. kmin) then
c: This is nothing to worry about, just a bad guess for kmid:
cc          print *,'kmin > klast: ',kmid,kmin,kminx,k1,k2
            kmin=max(k1,kminx)
         endif
         kmax=klast
         iim0=1
      else
         klast=kmax + 2.d0*dk_max0
         iim0=0
      endif
      ntry=0
c: Find next mode only:
30    continue
      ncx=nctot
      ntry=ntry + 1
cc      print *,'mode_int 30: ',ntry,kmax,kmin,knbb(jfmid,nmode),
cc     .   kminx,kmid,k1
      call rx_mode_int(kmax,kmin,jm1,ndv,iim0,0,1,0)
      nct_f(jfmid)=nct_f(jfmid) + nctot - ncx
      if(nmode .lt. jm1) then
c: No mode found over that interval, so increase interval:
         if(ntry .gt. 1) then
            print *,'trouble: ',jm1,nmode,kmax,kmin
            nmode=0
            call rx_freq_init(kmax,kmin)
            call rx_mode_int(kmax,kmin,jm1,ndv,0,0,1,0)
            f_hz=faxbb(jf1)
            nmode=0
            call rx_freq_init(kmax,kmin)
            call rx_mode_int(kmax,kmin,jm1,ndv,0,0,1,0)
cc    print *,'at f1: nmode = ',nmode,(kn(jm),jm=jm1-4,nmode)
            kmin=kminx
            stop 
         else
c: Leave kmax same so that pairs of closely spaced modes can be found better:
cc          kmax=kmin
cc          iim0=0
            kmin=max(k1,kminx)
         endif
cc    print *,'trying kminx = ',kminx,kmax
         goto 30
      endif
      if(abs(klast - dreal(kn(nmode))) .lt. 2.d0*dk_max0) then
         if(iidiag .ge. 1) print *,'Info msg: dup mode??',
     .      klast,kn(nmode)
         nmode=jm1-1
         if(klast .lt. dreal(kn(nmode))) then
            kmax=klast - 2.d0*dk_max0
            iim0=0
         else
            kmax=dreal(kn(nmode))
            iim0=1
         endif
         kmin=max(k1,kminx)
         ntry=0
         goto 30
      endif
c
cc    print *,'success by mode_int: ',k1,kn(nmode),k2
c
50    continue
cc	print *,'bb_enter call from fint: '
      call rx_bb_enter(jm1,jm1,jfmid,jdf,ndfmax,
     .   phi,dphi,Dphiz_w,Ddphiz_w,nmex)
c
c: Push right interval onto stack if it exists:
      if(jf2 .gt. jfmid+1) then
         nstack=nstack+1
         stack(1,nstack)=jfmid
         stack(2,nstack)=jf2
      endif
      if(jfmid .eq. jf1+1) then
c: Left interval empty, check stack:
         if(nstack .gt. 0) then
            jf1=stack(1,nstack)
            jf2=stack(2,nstack)
            nstack=nstack-1
            goto 10
         endif
c: Fall through to mode loop when no intervals left on stack.
      else
c: Prepare to do left interval:
         jf2=jfmid
         goto 10
      endif
99    continue
c
      return
      end
ccc
      subroutine rx_bb_mint(ndfmax,jm1,jm2,jf1,jflo,kminx,nmex,
     .   nct_f,ndv)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 ndfmax,jm1,jm2,jf1,jflo,nmex,ndv,nct_f(nfbb),jdf1,ncx,
     .   iim0
      real*8 kmin,kmax,kminx
c
      f_hz=faxbb(jf1)
      call rx_freq_init(kmax,kmin)
      if(kminx .gt. kmin) kmin=kminx
      nmode=jm1-1
      if(nmode .gt. 0) then
         kmax=dreal(knbb(jf1,nmode))
         iim0=1
      else
         kim_min=1.d100
         iim0=0
      endif
      ncx=nctot
      call rx_mode_int(kmax,kmin,jm2,ndv,iim0,0,1,0)
      nct_f(jf1)=nct_f(jf1) + nctot - ncx
cc    print *,'triangle: ',jf1,nm1x,nmode,nctot-ncx
      jdf1=jf1 - jflo + 1
cc      print *,'call enter from mint: '
      call rx_bb_enter(jm1,nmode,jf1,jdf1,ndfmax,
     .   phi,dphi,Dphiz_w,Ddphiz_w,nmex)
c
      return
      end
ccc
      subroutine interp_check(ndfmax,nm,phiz,jflo,jfhi)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 ndfmax,nm,jflo,jfhi,jf,jm,jdf,j0,jzsr
      real*8 dk_err,kmax,kmin
      real*4 ph_err,phi_ex
      complex*8 phi_(nzsr,nm)
c
c: TEMP: check how well interpolation worked:
      do jf=jflo+1,jfhi-1
         nmode=0
         f_hz=faxbb(jf)
         call rx_freq_init(kmax,kmin)
         kim_min=1.d100
         call rx_mode_int(kmax,kmin,nm_lim,3,0,1,1,0)
         do jm=1,nmode
            dk_err=dabs(dreal(kn(jm))-dreal(knbb(jf,jm)))
            if(dk_err .gt. dk_max0) then
               print *,jf,jm,dk_err/dk_max0,
     .            sngl(dreal(kn(jm))),sngl(dreal(knbb(jf,jm)))
            endif
            jdf=jf-jflo + 1
            do jzsr=2,2
               j0=1 + (jzsr-1)*4
               phi_ex=real(phi_(jzsr,jm))
               ph_err=100.*abs(phi_ex - phi_bb(1,jzsr,jdf,jm))/
     .            abs(phi_ex)
               if(ph_err .gt. 5. .and. abs(phi_ex) .gt. 1.e-4) then
                  print *,'jf,jm,j0 = ',sngl(f_hz),jf,jm,jzsr,ph_err,
     .               phi_ex,phi_bb(1,jzsr,jdf,jm)
               endif
            enddo
         enddo
      enddo
c
      return
      end
ccc
      subroutine quad_fit(k1,k1p,k1pp,c1,c2,c3,delf,kmid)
c
      implicit none
      real*8 k1,k1p,k1pp,c1,c2,c3,delf,kmid
c
      c3=0.5d0*k1pp
      c2=k1p
      c1=k1
      kmid=c1 + c2*delf + c3*delf*delf
c
      return
      end
ccc
      subroutine rx_bb_enter(jm1,jm2,jfbb,jdf,ndfmax,
     .   phi_,dphi_,Dphiz_w,Ddphiz_w,nmex)
c
      use i_o_com
      implicit none
      integer*4 jm1,jm2,jfbb,jdf,ndfmax,nmex,jm,jzsr
      complex*8 phi_(nzsr,nmode),dphi_(nzsr,nmode)
      real*8 Dphiz_w(nzsr,nmode),Ddphiz_w(nzsr,nmode)
c
cc      print *,'enter: ',jm1,jm2,jfbb,nfbb
      do jm=jm1,jm2
         knbb(jfbb,jm)=kn(jm)
         vg_bb(jfbb,jm)=eig_char(4,jm)
         Dvg_w(jfbb,jm)=eig_char(5,jm)
         iish_bb(jfbb,jm)=nzref(jm)
cc       call rx_bb_enterx(nzsr,phi_(1,jm),dphi_(1,jm),Dphiz_w(1,jm),
cc   .      Ddphiz_w(1,jm),phi_bb(1,1,jdf,jm))
         do jzsr=1,nzsr
            phi_bb(1,jzsr,jdf,jm)=real(phi_(jzsr,jm))
            phi_bb(2,jzsr,jdf,jm)=real(dphi_(jzsr,jm))
            phi_bb(3,jzsr,jdf,jm)=Dphiz_w(jzsr,jm)
            phi_bb(4,jzsr,jdf,jm)=Ddphiz_w(jzsr,jm)
         enddo
c: Compute field for this mode at this frequency:
         call rx_bb_field(kn(jm),phi_(1,jm),jm,jfbb)
         if(iift .eq. 5) write(88,108) faxbb(jfbb),jm
108      format(f9.2,2x,i4,2x,i1)
      enddo
      nmex=nmex + jm2-jm1+1
c
      return
      end
ccc
      subroutine rx_block_shift(nzsr,ndfmax,nmode,jflo,jfhi,
     .   phi_bbx,NaN)
c
c: Transfers mode info at jflo to jfhi.
      implicit none
      integer*4 nzsr,ndfmax,nmode,jflo,jfhi,j0,jm,n0,jf
      real*4 phi_bbx(4*nzsr,ndfmax,nmode),NaN
c
      n0=4*nzsr
      do jm=1,nmode
         do j0=1,n0
            phi_bbx(j0,jfhi,jm)=phi_bbx(j0,jflo,jm)
            do jf=jflo,jfhi-1
               phi_bbx(j0,jf,jm)=NaN
            enddo
         enddo
      enddo
c
      return
      end
ccc
      subroutine rx_mcross(nm1,nm2,jf1,jf2,jmlo,jmhi,
     .  ndrx,f1,f2,ccr_lo,ccr_hi,pie)
c
c: Finds mode interval over which modes might cross between frequencies
c: jf1 and jf2.
c
      use i_o_com
      implicit none
      integer*4 nm1,nm2,jf1,jf2,jmlo,jmhi,ndrx,jm,nm_min
      real*8 f1,f2,ccr_lo,ccr_hi,pie,kcr_lo,kcr_hi
c
      if(ndrx .le. 1) then
c: Easy case of zero or one duct (one-to-one correspondence):
         jmlo=1
         jmhi=0
      else
         kcr_hi=2.d0*pie*f1/ccr_lo
c: EKW FIX 6-2-97: Do not let jmhi exceed # modes at high or low freqs:
         nm_min=min(nm1,nm2)
c: Find first mode at f1 that could be propagating at ss ccr_lo:
c: EKW FIX 5-9-97 (in case all modes are well trapped):
         jmlo=nm_min+1
         do jm=1,nm_min
            if(dreal(knbb(jf1,jm)) .le. kcr_hi) then
               jmlo=jm
               goto 10
            endif
         enddo
10       continue
         kcr_lo=2.d0*pie*f2/ccr_hi
c: Find highest order mode at f2 that could be trapped in duct with
c: highest sound speed between ducts of ccr_hi:
         do jm=nm_min,1,-1
            if(dreal(knbb(jf2,jm)) .ge. kcr_lo) then
               jmhi=jm
               goto 20
            endif
         enddo
20       continue
c: To find modes more safely (avoid possibility of missing two at end of
c: interval), increase jmhi by one if possible:
         jmhi=min(jmhi+1,nm_min)
      endif
c
      return
      end
ccc
      subroutine rx_bb_brute
c
c: Does broadband calculations for modes on the real axis in the brute
c: force manner.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 jf,nmex,nmtot,nmint,jm,nm_fmax,jzjm,jsr,jj
      real*8 kmax,kmin
c
      if(iift .ne. 0) then
         call rx_bb_image
         print *,'Done rx_bb_image'
      endif
      nmint=0
c
      call bb_init
      f_hz=faxbb(nfbb)
      call rx_prep
      call rx_zmx_init
      nmex=0
      nmint=0
      do jfbb=nfbb,1,-1
         f_hz=faxbb(jfbb)
         call rx_freq_init(kmax,kmin)
         nmode=0
         kim_min=1.d100
         call rx_mode_int(kmax,kmin,nm_lim,3,0,1,0,0)
         if(jfbb .eq. nfbb) then
            nm_fmax=nmode
            do jm=1,nmode
               xmode(jm)=jm
            enddo
            allocate(knbb(nfbb,nm_fmax))
            allocate(eig_bb(5,nfbb,nm_fmax))
            allocate(iish_bb(nfbb,nm_fmax))
            if(iikn .gt. 0) then
               allocate(r4mat1(nfbb,nm_fmax))
               allocate(r4mat2(nfbb,nm_fmax))
            endif
            if(iimf .ne. 0) then
               allocate(r4mat_3d(nrec,nm_fmax,nfbb))
            endif
            if(iidc .ne. 0) then
               allocate(r4mat5(nfbb,nm_fmax))
               allocate(r4mat6(nfbb,nm_fmax))
            endif
c: Open mode characteristics file if desired:
            if(iiout .ne. 0) call bb_out_init(nm_fmax)
         else
            nmode=min(nmode,nm_fmax)
         endif
         if(iiout .ne. 0) then
            nmbb(jfbb)=min(nm_fmax,nm_tot)
            if(nm_tot .gt. nm_fmax) print *,'nm_tot>nm_fmax',
     .         f_hz,nm_tot,nm_fmax
            write(33,rec=nh_off+jfbb) (kn(jm),jm=1,nmbb(jfbb)),
     .         ((phi(jz,jm),jz=1,nzsr),jm=1,nmbb(jfbb))
         endif
c
         do jm=1,nmode
            knbb(jfbb,jm)=kn(jm)
            eig_bb(1,jfbb,jm)=eig_char(4,jm)
         enddo
c
         if(iimf .ne. 0) then
            do jrec=1,nrec
               jz=mzrec(jrec)
               do jm=1,nmode
                  r4mat_3d(jrec,jm,jfbb)=real(phi(jz,jm))
               enddo
               do jm=nmode+1,nm_fmax
                  r4mat_3d(jrec,jm,jfbb)=NaN
               enddo
            enddo
            call hdf_write_gen(3,outroot,lout,zrec,nrec,xmode,
     .         nm_fmax,faxbb,nfbb,1.,1,1.,1,r4mat3,1,nrec,
     .         1,nm_fmax,jfbb,jfbb,1,1,1,1,'Depth - m',9,'Mode',4,
     .         'Frequency - Hz',14,' ',1,' ',1,'Re(phi)',7,12,1,2)
         endif
         if(iidc .ne. 0) then
            do jm=1,nmode
               r4mat5(jfbb,jm)=real(eig_char(4,jm))
               r4mat6(jfbb,jm)=wbb(jfbb)/dreal(kn(jm))
            enddo
            do jm=nmode+1,nm_fmax
               r4mat5(jfbb,jm)=NaN
               r4mat6(jfbb,jm)=NaN
            enddo
         endif
         if(iikn .ne. 0) then
            do jm=1,nmode
               r4mat1(jfbb,jm)=dreal(kn(jm))
               r4mat2(jfbb,jm)=dimag(kn(jm))
               if(iikn .gt. 1) r4mat2(jfbb,jm)=r4mat2(jfbb,jm)/kw0
            enddo
            do jm=nmode+1,nm_fmax
               r4mat1(jfbb,jm)=NaN
               r4mat2(jfbb,jm)=NaN
            enddo
         endif
c
         do jm=1,nmode
c: Compute field for this mode at this frequency:
            call rx_bb_field(kn(jm),phi(1,jm),jm,jfbb)
         enddo
         nmtot=nmtot + nmode
         nmex=nmex + nmode
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
      endif
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
      if(iikn .ne. 0) then
         deallocate(r4mat1)
         deallocate(r4mat2)
      endif
      if(iimf .ne. 0) then
         deallocate(r4mat3)
      endif
      if(iidc .ne. 0) then
         deallocate(r4mat5)
         deallocate(r4mat6)
      endif
c
      if(iifft .ne. 0) call bb_fft_out
c
      write(2,120) nmtot,nmex,nmint,nctot,nctot/float(nmtot)
120   format('TOTAL # MODES = ',i8,'; # EXACT =',i8,'; # INTERP =',i8/
     .   '   # R1R2 CALCS = ',i8,'; # CALCS/MODE = ',f5.2)
c
      if(iiout .ne. 0) then
         close(33)
         open(33,file=outroot(1:lout)//'_bbeig',access='direct',
     .      recl=NRECL*lheadw)
         write(33,rec=1) fsbb,nfftbb,fmin,fmax,nfbb,nzsr,nm_fmax,
     .      rmin,rmax,sngl(cfmin),
     .      (zsr(jsr),jsr=1,nzsr),(nmbb(jfbb),jfbb=1,nfbb),iirx,
     .      (sngl(rho_sr(jsr)),jsr=1,nzsr),(faxbb(jfbb),jfbb=1,nfbb)
cc       print *,'nmbb at end= ',(nmbb(jfbb),jfbb=1,nfbb)
cc       print *,'lheadw,lrecw = ',lheadw,lrecw,nh_off
         close(33)
      endif
c
      return
      end
ccc
      subroutine rx_bb_image
c
c: Calculates a 2-D image of e(k) and e'(k) vs. k and f.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 ncos,nf,jcos,jf,jj0
      real*8 kmin,kmax,cos1,cos2,f1,f2,cosfac,ffac,k,ek(6),phm
c
      print *,'Enter cos1,cos2,ncos,f1,f2,nf: (cos1=0 ==> kmin/kmax)'
      read(5,*) cos1,cos2,ncos,f1,f2,nf
      f_hz=f1
      call rx_prep
      call rx_zmx_init
      call rx_freq_init(kmax,kmin)
      if(cos1 .eq. 0.d0) cos1=1.0001d0*kmin/kmax
c
      cosfac=(cos2 - cos1)/(max(1,ncos-1))
      ffac=(f2 - f1)/(max(1,nf-1))
      do jf=1,nf
         faxbb(jf)=f1 + (jf-1)*ffac
      enddo
      do jcos=1,ncos
         xmode(jcos)=cos1 + (jcos-1)*cosfac
      enddo
c
      allocate(r4mat5(nf,ncos))
      allocate(r4mat6(nf,ncos))
      do jf=1,nf
         f_hz=faxbb(jf)
         call rx_freq_init(kmax,kmin)
         jj0=(jf-1)*ncos
         do jcos=1,ncos
            k=max(kmin,xmode(jcos)*kmax)
            call rx_ek_calc(k,ek,phm,2)
            r4mat1(jf,jcos)=ek(1)
            r4mat2(jf,jcos)=phm
         enddo
      enddo
      call hdf_write_gen(2,outroot,lout,faxbb,nf,xmode,
     .   ncos,1.,1,1.,1,1.,1,r4mat5,1,nf,1,ncos,
     .   1,1,1,1,1,1,'Frequency - Hz',14,'cos(k/kw)',9,
     .   ' ',1,' ',1,' ',1,'e(k)',4,26,1,2)
      call hdf_write_gen(2,outroot,lout,faxbb,nf,xmode,
     .   ncos,1.,1,1.,1,1.,1,r4mat6,1,nf,1,ncos,
     .   1,1,1,1,1,1,'Frequency - Hz',14,'cos(k/kw)',9,
     .   ' ',1,' ',1,' ',1,'ep(k)',5,27,1,2)
c
      deallocate(r4mat5)
      deallocate(r4mat6)
c
      return
      end
ccc
      subroutine rx_bb_field(knx,phiz,jmx,jf)
c
c: Computes field tf for a single mode characterized by eigenvalue k and mode 
c: functions phi and dpsi.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jmx,jf,jr,jsr,kk,kk0,jsrc,jrec,jzs
      complex*8 phiz(nzsr)
      complex*16 knx,sqkn,cfac2,exp_sqk,iknx
      real*4 phi_src
c
c: Find mode function at source:
c: Normalize by density at source (9-1-94):
      jzs=mzsrc(1)
c: FIX 3-23-95:
      phi_src=phiz(jzs)/rho_sr(jzs)
      iknx=dcmplx(-dimag(knx),dreal(knx))
      sqkn=cdsqrt(knx)
      do jr=1,nrng
         exp_sqk=cdexp(iknx*range(jr))/sqkn
         cfac2=sq2pir(jr)*phi_src*exp_sqk
         kk0=krec_jr(jr)
         do kk=kk0+1,kk0+nrec_jr(jr)
            jsrc=jrec_jr(1,kk)
            jrec=jrec_jr(2,kk)
            jsr=mzrec(jrec)
            tf(jf,jrec,jsrc)=tf(jf,jrec,jsrc) + cfac2*phiz(jsr)
         enddo
      enddo
c
c: Total number of modes at frequency jf:
      nmbb(jf)=max(nmbb(jf),jmx)
c
      return
      end
ccc
      subroutine cub_fit_r4(k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,delk)
c
c: This subroutine fits a cubic polynomial y(x)=c1 + c2*x + c3*x**2 +
c: c4*x**3, where x=(k-k1)/(k2-k1), to the points (k1,y1) and (k2,y2)
c: and the derivatives y1p and y2p at those points.
c: Note that k is normalized from (k1,k2) to (0,1), so the fit values
c: must be computed accordingly: k=k1 + kx*delk.
c
      implicit none
      real*4 k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,b1,b2,delk,y1px,y2px
c
      delk=k2 - k1
      y1px=y1p*delk
      y2px=y2p*delk
      c1=y1
      c2=y1px
      b1=y2 - y1 - y1px
      b2=y2px - y1px
      c3=3.d0*b1 - b2
      c4=b2 - 2.d0*b1
c
      return
      end
ccc
      subroutine rx_bb_out(ndfmax,nm_fmax,jflo,jfhi,nh_off)
c
      use i_o_com
      implicit none
c
      integer*4 ndfmax,nm_fmax,jflo,jfhi,nh_off,j,jf,jm,jz
c
      j=0
      do jf=jflo,jfhi-1
         j=j+1
         write(33,rec=nh_off+jf) (knbb(jf,jm),jm=1,nm_fmax),
     .      ((phi_bb(1,jz,j,jm),0.,jm=1,nm_fmax),jz=1,nzsr)
      enddo
c
      return
      end
ccc
      subroutine rx_bb_Dvg(g12,g22,xkh,xkhsq,Ksq,w,wsq,twpie,
     .   vg,D2k_Df)
c
      implicit none
      real*8 g12(6),g22(6),xkh,xkhsq,Ksq,w,wsq,twpie,vg,D2k_Df,e_k(6),
     .   migam2(6),gamsq,Dvg_w
c
      gamsq=Ksq - xkhsq
      migam2(1)=sqrt(-gamsq)
      e_k(1)=g22(1) + migam2(1)*g12(1)
      migam2(2)=xkh/migam2(1)
      migam2(3)=-Ksq/(migam2(1)*w)
      e_k(2)=g22(2) + migam2(1)*g12(2) + migam2(2)*g12(1)
      e_k(3)=g22(3) + migam2(1)*g12(3) + migam2(3)*g12(1)
      migam2(4)=-migam2(2)*migam2(3)/migam2(1)
      migam2(5)=Ksq/(migam2(1)*gamsq)
      migam2(6)=Ksq*(Ksq/gamsq - 1.d0)/(wsq*migam2(1))
      e_k(4)=g22(4) + migam2(1)*g12(4) + migam2(4)*g12(1) + 
     .   migam2(2)*g12(3) + migam2(3)*g12(2)
      e_k(5)=g22(5) + migam2(1)*g12(5) + migam2(5)*g12(1) + 
     .   2.d0*migam2(2)*g12(2)
      e_k(6)=g22(6) + migam2(1)*g12(6) + migam2(6)*g12(1) + 
     .   2.d0*migam2(3)*g12(3)
c
      vg=-e_k(2)/e_k(3)
c: Total derivative of group velocity:
      Dvg_w=(vg*vg*e_k(6) + 2.d0*vg*e_k(4) + e_k(5))/e_k(2)
      D2k_Df=(-(twpie/vg)**2)*Dvg_w
c
      return
      end
