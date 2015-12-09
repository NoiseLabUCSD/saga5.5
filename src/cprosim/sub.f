      subroutine airy_hsp(ii1,ii2,iiAih,geolay,Aih_mag,Aih_dir,ihf,
     .   ii_xi_mv,hlay,f_3lam,pie,cphlo)
c
c: Introduces Airy halfspace instead of homogeneous halfspace.
c
      implicit none
c
      integer*4 ii1,ii2,ihf,ii_xi_mv
      real*4 iiAih,Aih_mag,Aih_dir,Aih_dc,Aih_da,lrat,Aihx
      real*8 geolay(2,5),hlay,f_3lam,pie,cphlo
c
      ii_xi_mv=1
      if(geolay(ii1,1) .le. cphlo) then
c: Don't insert gradient in air, for example:
         ii_xi_mv=0
         return
      elseif(geolay(1,1) .ne. geolay(2,1) .or.
     .   geolay(1,4) .ne. geolay(2,4)) then
c: Don't insert gradient in halfspace for which user has defined gradient 
c: explicitly
         print *,'Using user-defined Airy halfspace gradient ',
     .      geolay(ii1,1) - geolay(ii2,1),geolay(ii1,4) - geolay(ii2,4)
         ii_xi_mv=0
         return
      elseif(iiAih .eq. -1.) then
c: Flag to keep homogeneous halfspace
         ii_xi_mv=0
         return
      elseif(iiAih .lt. 0.) then
         ihf=1
      else
         ihf=0
         hlay=3.d0*geolay(ii1,1)/f_3lam
      endif
c
      Aihx=abs(iiAih)
      if(Aihx .eq. 0.) then
c: Default gradient such that modes move straight up in k-plane:
         Aih_dir=135.
         Aih_mag=0.001
      elseif(Aihx .lt. 1.) then
c: Default gradient such that modes move straight up in k-plane:
         Aih_dir=135.
         Aih_mag=Aihx
      elseif(Aihx .lt. 2.) then
c: 1.005 means you want gradient in c only:
         Aih_dir=0.
         Aih_mag=Aihx - int(Aihx)
      else
         Aih_dir=int(Aihx)
         Aih_mag=Aihx - int(Aihx)
      endif
      if(Aih_dir .lt. 0. .or. Aih_dir .gt. 180.) then
         print *,'Error: Airy gradient direction ',
     .      '[atan2(da,dc)] must be in [0,180] interval.'
         stop
      endif
      if(Aih_mag .lt. 0.0005) then
         print *,'Error: Airy gradient magnitude must ',
     .      'be > 0.0005 (recommended interval = [.001,.01]).'
         stop
      endif
      if(Aih_mag .gt. 0.5) then
         print *,'Warning: Airy gradient should be ',
     .      'less than 0.5 (recommended interval = [.001,.01]).'
      endif
      Aih_dc=Aih_mag*cos(Aih_dir*pie/180.)
      Aih_da=Aih_mag*sin(Aih_dir*pie/180.)
      geolay(ii2,1)=geolay(ii1,1) + Aih_dc*geolay(ii1,1)
c: Convert attenuation to dB/m-kHz:
      geolay(ii2,4)=geolay(ii1,4) + Aih_da*54.575/(.001*geolay(ii1,1))
c: Check if same gradient needs to be applied to shear wave profile:
      if(geolay(ii1,2) .gt. 0.d0) then
c: Multiply by ratio of wavelengths to convert thickness of layer
c: (hlay) to wavelength of shear waves:
         lrat=geolay(ii1,1)/geolay(ii1,2)
         geolay(ii2,2)=geolay(ii1,2) + lrat*Aih_dc*geolay(ii1,2)
         geolay(ii2,5)=geolay(ii1,5)+lrat*Aih_da*54.575/
     .      (.001*geolay(ii1,2))
      endif
c
      return
      end
ccc
      subroutine phimag_check(phiz,dphiz,phim_max,k,jm)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 jm,nsv,isv,jd,jzsr,ndone,ndone_max,iiask
      real*4 phim_max,magsq_c8,ampl,abar,ki_eff,phiw_max,
     .   xi_re,xi_rat,xi_res,xi_rats,xi_rat_max,phiw_last,
     .   many_fac
      complex*8 phiz(nzsr,jm),dphiz(nzsr,jm),dphi_z
      complex*16 k
c      data ndone/0/,iiask/0/,xi_rat_max/20./
c
      ndone=0
      iiask=0
      xi_rat_max=20.
c
cc    if(iiask .eq. 0) then
cc       print *,'Enter max # weak BLM modes to include'
cc       read(5,*) ndone_max
cc       ndone_max=0
cc       iiask=1
cc    endif
      ndone_max=0
      if(nmode .eq. 1) then
         ndone=0
      endif
      iiblm=0
c
csi   phim_max=max(phim_max,phiw_max)
c: Simple estimate of maximum amplitude in waveguide:
      phim_max=dsqrt(2.d0/Htot)
      if(nsvmin .ne. nlay) then
         xi_re=dreal(xi_hsp(1,1))
         xi_rat=abs(xi_re/dimag(xi_hsp(1,1)))
         if(xi_re .lt. 0. .and. xi_rat .gt. xi_rat_max) then
            if(iisol(nlay) .eq. 0) then
c: Only set flag if BLMs will be found separately:
               iiblm=1
            else
               iiblm=-1
            endif
            if(iidiag .ge. 1) print *,
     .         'BLM non-hspace: ',sngl(dreal(k)/kw0),
     .         sngl(dimag(k)*8685.9),xi_re,xi_rat,sngl(f_hz)
         elseif(iisol(nlay) .eq. 1) then
            xi_res=dreal(xi_hsp(1,2))
            xi_rats=abs(xi_res/dimag(xi_hsp(1,2)))
            if(xi_res .lt. 0. .and. xi_rats .gt. xi_rat_max) then
               if(iisol(nlay) .eq. 0) then
c: Only set flag if BLMs will be found separately:
                  iiblm=2
               else
                  iiblm=-2
               endif
            endif
         endif
      endif
c: EKW FIX: only increment ndone when nsvmin=nlay, otherwise
c: fewer than ndone_max BLM will be kept:
cc    if(nsvmin .eq. nlay .or. iiblm .gt. 0) then
cc0   if(nsvmin .eq. nlay .or. (ii_xi_ax .gt. 0 .and.
cc0  .   iiblm .eq. 0)) then
cc    if(nsvmin .eq. nlay) then
c
      if(iiblm .ne. 0) then
c: Find maximum amplitude of mode inside waveguide:
         phiw_max=0.
         do jd=1,nduct
            jzsr=mzduct(jd)
            nsv=jduct(1,jd)
            isv=jduct(2,jd)
            dphi_z=dphiz(jzsr,jm)/gami(1,isv,nsv)
            ampl=magsq_c8(phiz(jzsr,jm)) + magsq_c8(dphi_z)
            if(nsv .ne. nlay) phiw_max=max(phiw_max,ampl)
         enddo
c
         if(phiw_max .lt. phiw_last) then
c: Account for potential of linear decrease to zero of mode amplitudes:
            many_fac=0.5*phiw_max/(phiw_last - phiw_max)
         else
c: If mode amplitudes are rising, make factor very safe:
            many_fac=100.
         endif
c: Check if mode is weak because of weakness in mode function inside
c: waveguide:
cc       abar=phiw_max/phim_max
         abar=many_fac*phiw_max/phim_max
         ki_eff=dimag(k) - log(abar)/(1000.*rmin)
cpln         print *,jm,phiw_last,phiw_max,many_fac,ki_eff,kim_max
         if(iidiag .ge. 2) print *,'ki_eff = ',jm,
     .      8.6859*log(abar)/rmin,sngl(8685.9*kim_max),abar
         if(ki_eff .gt. kim_max) then
            ndone=ndone + 1
cc          if(ndone .gt. ndone_max) then
               iidone=1
               nhigh=max(1,nhigh)
cc          endif
            if(iidiag .ge. 1) then
               print *,'Airy contour done due to phimag weakness',
     .            jm,(8.6859*log(abar)/rmin),xi_re,xi_rat
            endif
         endif
         phiw_last=phiw_max
      else
         phiw_last=0.
      endif
c
      return
      end
      subroutine bb_field(knx,phiz,dphiz,dpsiz,expz_gbs,jmk,tfz,
     .   jfbb,jmo)
c
c: Computes field tf for a single mode characterized by eigenvalue k and mode 
c: functions phi and dpsi.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jmk,jfbb,jmo,jr,jsr,kk,kk0,jsrc,jrec,j,jzs,msrc,nrmax,
     .   jd,nsv,isv,jzsr
      complex*8 phiz(nzsr,jmk),dphiz(nzsr,jmk),dpsiz(nzsr,jmk),
     .   tfz(nfbb,nrec,nrng),dphi_z
      complex*16 knx,sqkn,cfac2,h_arg,H0,iknx,phi_src
      real*8 expz_gbs(nzsr,jmk),xn_beam,rsig_max,magsq
      real*4 phiw_max,magsq_c8,ampl
c
c: Find mode function at source:
c: Normalize by density at source (9-1-94):
      jzs=1
      msrc=mzsrc(jzs)
      phi_src=phiz(msrc,jmk)/rho_sr(msrc)
c
cpln      write(6,*)zsr(msrc),phi_src
c
      xn_beam=(w/dreal(cp_sr(msrc)))*b_gbs(jzs) - expz_gbs(msrc,jmk)
      iknx=dcmplx(-dimag(knx),dreal(knx))
      sqkn=cdsqrt(knx)
      do jrec=1,nrec
         jsr=mzrec(jrec)
         phisr(jsr)=phiz(jsr,jmk)
         if(kksh(jsr) .eq. 1) then
c: Add shear wave potential contributions to recs in shear layers (see p. 117):
            phisr(jsr)=phisr(jsr) - 2.d0*ksm2_sr(jsr)*
     .         (knx*knx*phisr(jsr) + iknx*dpsiz(jsr,jmk))
         endif
      enddo
c
c: Find range beyond which mode is insignificant (ranges have been sorted):
      rsig_max=kim_fac/dmax1(1.d-20,dimag(knx)-kim_bb(jfbb))
      nrmax=0
      call hunt(range,nrng,rsig_max,nrmax)
c
      do jr=1,nrmax
         kk0=krec_jr(jr)
         if(tilth) then
            kk=kk0+jr
            jsrc=jrec_jr(1,kk)
            jrec=jrec_jr(2,kk)
            h_arg=knx*range(jr)
            if(magsq(h_arg) .gt. 25.d0) then
czs: Include normalization by exp(-xn_beam) here:
               cfac2=sq2pir(jr,jrec)*phi_src*
     .              cdexp(iknx*range(jr) - xn_beam)/sqkn
            else
               call cdhankel(h_arg,1.d-6,H0)
               cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
            endif
            jsr=mzrec(jrec)
            tfz(jfbb,jrec,jsrc)=tfz(jfbb,jrec,jsrc) + 
     .           cfac2*phisr(jsr)
         else
            do kk=kk0+1,kk0+nrec_jr(jr)
               jsrc=jrec_jr(1,kk)
               jrec=jrec_jr(2,kk)
               h_arg=knx*(range(jr)+dtiltvp(jrec))
               if(magsq(h_arg) .gt. 25.d0) then
czs: Include normalization by exp(-xn_beam) here:
                  cfac2=sq2pir(jr,jrec)*phi_src*
     .                 cdexp(iknx*(range(jr)+dtiltvp(jrec))
     .                  - xn_beam)/sqkn
               else
                  call cdhankel(h_arg,1.d-6,H0)
                  cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
               endif
               jsr=mzrec(jrec)
               tfz(jfbb,jrec,jsrc)=tfz(jfbb,jrec,jsrc) + 
     .         cfac2*phisr(jsr)
            enddo
         end if
      enddo
c
c: Keep track of strongest mode amplitude outside of halfspace:
      phiw_max=0.
      do jd=1,nduct
         nsv=jduct(1,jd)
         if(nsv .ne. nlay) then
            jzsr=mzduct(jd)
            isv=jduct(2,jd)
            dphi_z=dphiz(jzsr,jmk)/gami(1,isv,nsv)
            ampl=magsq_c8(phiz(jzsr,jmk)) + magsq_c8(dphi_z)
            phiw_max=max(phiw_max,ampl)
         endif
      enddo
      phim_bb(jfbb)=max(phim_bb(jfbb),phiw_max)
c
      if(iiout .ne. 0 .or. iimf .ne. 0) then
c: Fill phibb(nfbb,nzsr),dpsibb(nfbb,nzsr) with mode functions:
         do jsr=1,nzsr
            j=(jsr-1)*nfbb + jfbb
            phibb(j)=phiz(jsr,jmk)
            dpsibb(j)=dpsiz(jsr,jmk)
         enddo
      endif
c
      return
      end
      subroutine bb_init
c
c: Checks and initializes variables for broadband mode calculations.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 n,nfft0,j,iibad,ii1st,jr,jzs,jd
      real*8 rng_im
      complex*16 cfac
      data cfac/(1.77245385090552d0,1.77245385090552d0)/
c      data ii1st/1/

      ii1st=1
c
      if(ii1st .eq. 0) goto 10
c
c: Initialize source and receiver geometry arrays:
      call sr_geom(rng_sr,iabs(nsrc),iabs(nrec))
      if(iigbs .eq. 0) then
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/dsqrt(range(jr)+dtiltvp(jd))
            end do
cpln            sq2pir(jr)=cfac/dsqrt(range(jr))
         enddo
      else
         jzs=1
         rng_im=-b_gbs(jzs)*cos(th_gbs(jzs)*pie/180.d0)
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .              +dtiltvp(jd),rng_im+dtiltvp(jd)))
            end do
cpln            sq2pir(jr)=cfac/cdsqrt(dcmplx(range(jr),rng_im))
         enddo
      endif
10    continue
c
      if(iicw .eq. 2) then
         if(iiwrite.ge.1) then
            write(6,*)'Delta frequency= ',df
            write(6,*)'Maximum time window= ',1./df
c     fmc
cpln            print *,' from subroutine bb_init.f :'
            print *,' fs,   nfft, Tw ',  fsbb, nfftbb, Tw
            print *,' fmin, fmax, df :', fmin, fmax, df
            print *,' nf1,  nf2,  nfreq :', nf1, nf2, nfbb
c     fmc
            iibad=0
            call mem_lim(nfbb,NFBBMAX,MLINE,LML,
     .           'nfbb',4,'NFBBMAX',7,iibad,0)
            call mem_lim(nfbb*nrec*nsrc,NTFMAX,
     .           MLINE,LML,'nfbb*nrec*nsrc',14,
     .           'NTFMAX',6,iibad,0)
            if(iibad .eq. 1) then
              write(*,*)'   nfbb,NFBBMAX,nfbb,nrec,nsrc'           
              write(*,*)   nfbb,NFBBMAX,nfbb,nrec,nsrc          
               stop 'from bb_init'
c     
            endif
         endif
      endif
c
      if(i_geom .ne. 0) then
         do j=1,nfbb
            nmbb(j)=0
         enddo
      end if
c
      do j=1,nfbb*nrec*nsrc
         tf(j)=(0.,0.)
      enddo
      do j=1,nfbb
         faxbb(j)=fmin + (j-1)*df
         wbb(j)=twpie*faxbb(j)
         kim_bb(j)=1.d100
         phim_bb(j)=0.
      enddo
c
      if(iimt .ne. 0) then
         open(21,file=outroot(1:lout)//'_mtraj',status='unknown',
     .        form='formatted')
      endif
c
      ii1st=0
c
      return
      end
ccc
      subroutine bb_out_init(nm_fmax)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 nm_fmax
c
c: Open direct access for broadband mode characteristics, if desired:
      if(iirx .eq. -1) then
         lrecw=(4+4*nzsr)*nfbb
      elseif(iirx .eq. 0) then
c: (Complex*16 eigenvalue + complex*8 mode function at nzsr depth) at
c: each frequency bin:
         lrecw=4*nm_fmax + 2*nzsr*nm_fmax
      else
c: (Complex*16 eigenvalue + real*4 mode function at nzsr depth) at
c: each frequency bin:
         lrecw=(4+nzsr)*nm_fmax
      endif
      open(33,file=outroot(1:lout)//'_bbeig',status='unknown',
     .   access='direct',recl=NRECL*lrecw)
      lheadw=11 + 3*nzsr + nfbb
      nh_off=(lheadw - 1)/lrecw + 1
c
      return
      end
      subroutine bb_fmarch(jmk,jmo,k,r1r2,kp0,Lp0,dL_dwp0,dL_dkp0,
     .   iisg,dfinc,dfincx,nm2,jfbb2,iixi_cut,kdone)
c
c: Marches mode to new frequency.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jmk,jmo,iisg,ntry,nm2,jfbb2,iixi_cut,kdone
      complex*16 k,r1r2(3,4),k0,rr0(3,4),kp0,Lp0,dL_dwp0,dL_dkp0,
     .   kp,Lp,dL_dwp,dL_dkp,lnr1r2,branch,delk,dW_dk
      real*8 f_hz0,delf,dfinc,dfincx,df_fac,Lmagsq,Lerrok,magsq,
     .   frac_done,f_final
      real*8 xi_real,xi_rat
c
      frac_done=0.
      f_hz0=f_hz
      f_final=f_hz + iisg*df
c: Save input values of kp, Lp, etc. in case merged modes found:
      kp=kp0
      Lp=Lp0
      dL_dwp=dL_dwp0
      dL_dkp=dL_dkp0
c
10    ntry=0
20    continue
      ntry=ntry + 1
      delf=iisg*dfinc*df
      f_hz=f_hz0 + delf
      if(f_hz .lt. fmin) then
         f_hz=fmin
         delf=f_hz - f_hz0
         dfinc=delf/(iisg*df)
      elseif(f_hz .gt. fmax) then
         f_hz=fmax
         delf=f_hz - f_hz0
         dfinc=delf/(iisg*df)
      endif
c: Obtain guess for eigenvalue at new frequency using derivatives with
c: respect to k and w:
cxx   print *,'delk dir = ',atan2(dimag(delk),dreal(delk))*180./pie
      call freq_chng
      delk=-(Lp + dL_dwp*delf*twpie)/dL_dkp
      k=kp + delk
      call r1r2_calc(k,r1r2,3,1,jjfail)
      if(jjfail .gt. 0) return
      nctot=nctot + 1
      lnr1r2=r1r2(1,4)
      if(ntry .gt. 5) then
         print *,'jm,f,ntry = ',jmo,f_hz,ntry,dfinc
         print *,'Re(k): ',dreal(kp0)/kw0,dreal(kp)/kw0,dreal(k)/kw0
         print *,'Im(k): ',dimag(kp0)*8686.,dimag(kp)*8686.,
     .      dimag(k)*8686.
         print *,'Lp0,Lp,L = ',Lp0,Lp,lnr1r2,Lmagsq,Lerrok
         print *,'delk = ',dreal(delk)/kw0,dimag(delk)*8686
         print *,'iish,iish_ref,iishn = ',iish,iish_ref,iishn(jmk)
c        call r1r2_calc(kp,r1r2,3,1)
c        print *,'kp again: ',Lp,r1r2(1,4),dL_dkp,r1r2(2,4),
c    .      dL_dwp,r1r2(3,4)
cc       iidiag=1
         if(ntry .gt. 8) then
            print *,'Stopping ...'
            stop
         endif
      endif
c
      Lmagsq=magsq(lnr1r2)
      Lerrok=magsq(r1r2(2,4))*errdkms
      dfincx=dfinc
c: FIX 4-12-93: Check magnitude squared of ln(R1*R2) also:
      if(Lmagsq .le. Lerrok .and. Lmagsq .lt. .01d0) then
c: Guess perfect:
         if(abs(dimag(lnr1r2)) .gt. .157) then
            print *,'arg(r1r2) > .157 !!',nmode,f_hz,lnr1r2
         endif
         if(abs(dreal(lnr1r2)) .gt. .1) then
            print *,'mag(r1r2) > .1 !!',nmode,f_hz,lnr1r2
         endif
         if(Lerrok/Lmagsq .gt. 6.d0) dfincx=2.*dfinc
cc       if(Lerrok/Lmagsq .gt. 6.d0) dfinc=2.*dfinc
      elseif(abs(dreal(lnr1r2)) .lt. .1 .and.
     .   abs(dimag(lnr1r2)) .lt. .3927) then
c: Guess ok, but go back if trying to jump more than df since linear 
c: interpolation won't be so good here:
         if(dfinc .gt. 1.01d0) then
            df_fac=(Lerrok/Lmagsq)**(0.25)
ctt   print *,'Guess OK, but dfinc > 1: ',dfinc,df_fac,dfinc*df_fac
            dfinc=dmax1(1.d0,dmin1(dfinc-1.d0,
     .         dfloat(nint(dfinc*df_fac))))
ctt   print *,'dfinc set to ',dfinc
            call xkh_backup
            goto 20
         endif
c: Use eig_final to zero in on eigenvalue:
         call eig_final(k,r1r2,errdkms,0,kw0,iish,iifail,jjfail)
         if(iifail .eq. 1) then
            print *,'iifail=1 in eig_final from bb_fmarch'
            return
         endif
      else
c: Guess bad:
         dfinc=.5d0*dfinc
         if(dfinc .gt. 1.d0) dfinc=dfloat(nint(dfinc))
         call xkh_backup
         goto 20
      endif
c
      if(dimag(k) .lt. 0.d0 .and. dimag(k)*dimag(k) .gt. errdkms) 
     .   then
cxx      if(iidiag .ne. 0) print *,'MODE FOUND IN LOWER HP: ',
cxx  .      nmode,k,k/kw
cxx      iilhp=1
cxx      goto 888
         if(iidiag .ge. 1) print *,'neg im(k): ',k,k/kw0,f_hz,iish
      endif
c: Check for negative sheets, but don't change sheet to pos if found:
      if(isp(nlay) .eq. 1) then
         if(iish(1,1) + iish(2,1) .lt. 2) then
            kdone=1
            print *,'Mode crossed Pekeris cut: ',jmk,
     .         dreal(k)/kw0,dimag(k)*8685.9
            return
         endif
         call sheet_look(0,branch)
         if(branch .ne. (0.,0.)) then
            if(dfinc .gt. 1.001d0) then
               dfinc=1.d0
               call xkh_backup
               goto 20
            endif
            if(dreal(k) .lt. dreal(branch) .or. dimag(k) .lt. 
     .         dimag(branch)) then
               call fix_path(k0,rr0,1)
               if(iidone .eq. 1 .or. iifail .eq. 1) then
                  print *,'bb fix_path trouble',iidone,iifail
                  return
               endif
c: EKW FIX: make phcut=2*pi so that dph_left is close to zero in dph_calc:
               phcut=twpie
               call eig_find0(k0,rr0,k,r1r2,nmode+10)
               if(iidone .eq. 1 .or. iifail .eq. 1) then
                  print *,'bb eig_find0 trouble',iidone,iifail
                  return
               endif
            else
               if(iidiag .ge. 1) print *,'OK(?) MODE ON -1 SHEET ',
     .            'FOUND: ',nmode,sngl(f_hz),iish
            endif
         endif
cc    elseif(dfinc .lt. .999d0) then
cc    elseif(iixi_cut .eq. 0) then
      else
         if(-iisg*dimag(xi_hsp(1,1)) .lt. 0.d0) then
            kdone=1
            print *,'fmarch crossed gradient BL: ',jmk,xi_hsp(1,1),
     .         dreal(k)/kw0,dimag(k)*8685.9
            return
         elseif(iisol(nlay) .eq. 1) then
            if(-iisg*dimag(xi_hsp(1,2)) .lt. 0.d0) then
               kdone=1
               print *,'fmarch crossed gradient BL: ',jmk,xi_hsp(1,1),
     .            dreal(k)/kw0,dimag(k)*8685.9
               return
            endif
         endif
         xi_real=dreal(xi_hsp(1,1))
         xi_rat=dabs(xi_real/dimag(xi_hsp(1,1)))
         if(nsvmin .ne. nlay) then
            xi_real=dreal(xi_hsp(1,1))
            xi_rat=dabs(xi_real/dimag(xi_hsp(1,1)))
            if(xi_rat .gt. 20. .and. xi_real .lt. -1.) then
ccc         dW_dk=-2.*r1r2(1,1)*gamiref*r1r2(2,3)
ccc         call mode_fun(k,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
ccc  .         psi,dpsi,exp_gbs,nm2)
c: Compute mode amplitudes at duct reference depths:
ccc         call phimag_check(phi,dphi,phim_bb(jfbb2),k,nm2)
ccc         if(iidone .eq. 1) then
      print *,'BLM CROSS: ',jmo,sngl(f_hz),sngl(dfinc),xi_rat,xi_real
ccc            iidone=0
               f_hz=f_final
               call freq_chng
               call bb_blm_exit2(k,r1r2)
               iixi_cut=1
               dfinc=1.
               dfincx=1.
            endif
         endif
      endif
      nsave=nsave + max(0,nint(dfinc-1.d0))
cpln      write(6,*)'nsave,max0: ',nsave,max(0,nint(dfinc-1.d0)),dfinc
cpln      pause
c     
c: Check if entire freq bin has been crossed:
      if(dfinc .lt. .999d0) then
c: If only jumped by a fraction of df, check if a whole df done:
         frac_done=frac_done + dfinc
c: Allow dfinc to increase again:
         dfinc=dfincx
         if(frac_done .lt. .999) then
c: If total fraction less than 1, update f_hz0 and jump more in f:
            if(iift .eq. 1) call bb_write(jmo,jfbb2,faxbb,k,kw0,iish)
            f_hz0=f_hz
            kp=k
            Lp=r1r2(1,4)
            dL_dkp=r1r2(2,4)
            dL_dwp=r1r2(3,4)
            dfinc=min(dfinc,1.-frac_done)
            goto 10
         endif
      endif
c
      return
      end
ccc
      subroutine bb_blm_exit(k,r1r2)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      real*8 I0_re,del_xi_zero
      complex*16 k,r1r2(3,4),k0,rr0(3,4),beta,delk,lnr1r2,I_k,del_I,
     .   xi_lower
c
      k0=k
      beta=(etasq(nlay)/k0)*dcmplx(0.d0,-1.d0)
      delk=(.01d0*pie/Htot)*beta/cdabs(beta)
      k0=k0 + delk
      xi_lower=xi_hsp(1,1) - 
     .   dcmplx(0.d0,min(2.d0,-.05d0*dreal(xi_hsp(1,1))))
      k0=cdsqrt(xksq(1,nlay) + xi_lower*etasq(nlay))

      del_xi_zero=pie/dsqrt(dabs(dreal(xi_hsp(1,1))))
      xi_lower=xi_hsp(1,1) - dcmplx(0.d0,4.d0*del_xi_zero)
      k0=cdsqrt(xksq(1,nlay) + xi_lower*etasq(nlay))
      call r1r2_calc(k0,rr0,3,1,jjfail)
      if(jjfail .gt. 0) return
c
10    continue
      lnr1r2=rr0(1,4)
      I0_re=dreal(lnr1r2)
      if(dabs(I0_re) .gt. .05d0) then
         I_k=rr0(2,4)
         beta=(-I_k*etasq(nlay)/k0)
         del_I=dcmplx(-I0_re,-I0_re*dimag(beta)/dreal(beta))
         delk=del_I/I_k
         k0=k0 + delk
c
         call r1r2_calc(k0,rr0,3,1,jjfail)
         if(jjfail .gt. 0) return
         print *,'I = ',dreal(k0)/kw0,dimag(k0)*8685.9,lnr1r2,
     .      rr0(1,4),xi_hsp(1,1)
         goto 10
      endif
c
c: Cal traj_mag if not close enough to mode:
      iimst=-1
      phcut=2.d0*pie
      call eig_find0(k0,rr0,k,r1r2,0)
cc    dW_dk=-2.*r1r2(1,1)*gamiref*r1r2(2,3)
cc    call mode_fun(k,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
cc   .   psi,dpsi,exp_gbs,nm2)
c
      return
      end
ccc
      subroutine bb_blm_exit2(k,r1r2)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      real*8 del_xi_zero,ln_re1,ln_re2,ln_re0
      complex*16 k,r1r2(3,4),k0,rr0(3,4),xi1,xi2,xi0
      integer*4 ntry
c
      xi1=xi_hsp(1,1)
      call i_xi_calc(xi1,k0,rr0,ln_re1)
c
      ntry=0
      if(ln_re1 .lt. 0.d0) then
10       continue
         xi2=dcmplx(2.d0*dreal(xi1),dimag(xi1))
         ntry=ntry + 1
         if(ntry .gt. 10) print *,'ntry trouble in exit'
         call i_xi_calc(xi2,k0,rr0,ln_re2)
         if(ln_re2 .lt. 0.d0) then
            xi1=xi2
            goto 10
         endif
      else
         xi2=xi1
15       continue
         xi1=dcmplx(0.50*dreal(xi2),dimag(xi2))
         ntry=ntry + 1
         if(ntry .gt. 10) print *,'ntry trouble in exit'
         call i_xi_calc(xi1,k0,rr0,ln_re1)
         if(ln_re1 .gt. 0.d0) then
            xi2=xi1
            goto 15
         endif
      endif
c
c: We now have |R1R2|=1 bracketed on xi axis:
20    xi0=0.5d0*(xi1 + xi2)
      call i_xi_calc(xi0,k0,rr0,ln_re0)
      ntry=ntry + 1
      if(ntry .gt. 20) print *,'ntry trouble in exit'
      if(dabs(ln_re0) .gt. 0.2d0) then
         if(ln_re0 .lt. 0.d0) then
            xi1=xi0
         else
            xi2=xi0
         endif
         goto 20
      endif
      print *,'ntry final = ',ntry
c
c: Call traj_mag if not close enough to mode:
      iimst=-1
      phcut=2.d0*pie
cc    iimst=0
      call eig_find0(k0,rr0,k,r1r2,0)
c
      return
      end
ccc
      subroutine i_xi_calc(xi0,k0,rr0,i0)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      complex*16 xi0,k0,rr0(3,4)
      real*8 i0,del_xi_zero
c
      del_xi_zero=pie/dsqrt(dabs(dreal(xi0)))
      xi0=dcmplx(dreal(xi0),dimag(xi0)-4.d0*del_xi_zero)
      k0=cdsqrt(xksq(1,nlay) + xi0*etasq(nlay))
      call r1r2_calc(k0,rr0,3,1,jjfail)
      if(jjfail .gt. 0) return
      i0=dreal(rr0(1,4))
c
      return
      end
      subroutine bb_mloop(jm0,nmode_bb,nmbbtot,iisg,f_st,jf_st,
     .   ndied)
c
c: Loops over nmode_bb modes, whose indices for kn,phi,eig_char lie in 
c: kn_indx, and calls routine for following them in frequency.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm0,nmode_bb,nmbbtot,iisg,jf_st,jmx,jmk,
     .   jmo,jfbb,nfinc,jfbb2,nm1,nm2,nm3,jj0,
     .   nsave0,ii,ncmode,nm_fmax,iimrg,jm_mrg,nfout,kductp,iichng,
     .   kdone,nf_done,ndied,nf_blm,iixi_cut
      real*8 dfinc,dfincx
      real*4 f_st
      complex*16 k,r1r2(3,4),Lp,dL_dkp,dL_dwp,dW_dk,kp,vg,k_mrg
c
      nm_fmax=kn_indx(nmode_bb)
      nsave=0
      nskip=0
      ndied=0
      do jmx=1,nmode_bb
         jmk=kn_indx(jmx)
         jmo=jm0 + jmx
         nhigh=0
         iidone=0
         nrise=0
         iilk=0
         iixi_cut=0
         kim_max=0
         nsave0=nsave
c: Initialize frequency to that where modes were found:
         f_hz=f_st
c: Re-initialize frequency so that roundoff errors don't cause trouble:
         call freq_init
c: Initialize sheets and set xkhratp to kn(jmk)/kw (at f_st):
cc       call sheet_init(kn(jmk),1,iish,iish_ref)
c: Initialize sheet variables as found for mode jmk at f_st:
         xkhrat=kn(jmk)/kw0
         call iish_code(iish,iish_ref,iishn(jmk),-1)
c: Compute field at all receivers for this mode at f_st:
         jfbb=jf_st
         call bb_field(kn(jmk),phi,dphi,dpsi,exp_gbs,jmk,tf,jfbb,jmo)
         nm1=jmk
         nm2=nm_fmax + 2
         nm3=nm_fmax + 3
         ii=-1
         call bb_enter(kn(jmk),eig_char(1,jmk),eig_char(4,jmk),
     .      eig_char(5,jmk),knbb,eig_bb,jfbb,jmo,nfbb,nm_fmax,
     .      iish,iish_ref,iish_bb,nmbb)
         if(iift .eq. 1) call bb_write(jmo,jfbb,faxbb,kn(jmk),kw0,iish)
c
         dfinc=1.d0
         kp=kn(jmk)
         Lp=eig_char(1,jmk)
         dL_dkp=eig_char(2,jmk)
         dL_dwp=eig_char(3,jmk)
c: Check for change in reference depth:
         if(nzref(jmk) .ne. kduct) then
            kduct=nzref(jmk)
            call zref_chng
            if(iidiag .ge. 2) print *,'changed reference '//
     .         'depth for mode ',jmo
         endif
         ncmode=nctot
c: For each mode found at f_st, track it as a function of frequency:
         kdone=0
         nf_done=0
         nf_blm=0
         do while(kdone .eq. 0)
15          call bb_fmarch(jmk,jmo,k,r1r2,kp,Lp,dL_dwp,dL_dkp,
     .         iisg,dfinc,dfincx,nm2,jfbb2,iixi_cut,kdone)
            if(kdone .eq. 1) goto 50
            nfinc=max(1,nint(dfinc))
            jfbb2=jfbb + iisg*nfinc
            phim_bb(jfbb2)=.5d0*phim_bb(jfbb)
            if(iifail .eq. 1) then
               print *,'Mode ',jmk,' could not be tracked past f = ',
     .            f_hz
               iifail=0
               goto 50
            endif
c
c: Check for mode merging:
            call bb_merge(k,knbb,nfbb,nm_fmax,jfbb2,jmo,errdk2,iimrg,
     .         jm_mrg,k_mrg,nmerge,faxbb,iidiag)
            if(iimrg .eq. 1 .or. iimrg .eq. 2) then
c: If modes merged, first try going back to prev freq with current mode:
c: iimrg=1 ==> dfinc=.125, iimrg=2 ==> dfinc=.05:
               dfinc=.125/((iimrg-1)*3 + 1)
               f_hz=faxbb(jfbb)
               goto 15
            elseif(iimrg .eq. 3 .or. iimrg .eq. 4) then
c: If modes merged, next try going back to prev freq with duplicate mode:
               f_hz=faxbb(jfbb)
               call freq_chng
               jj0=(jm_mrg - 1)*nfbb + jfbb
               call iish_code(iish,iish_ref,iish_bb(jj0),-1)
               if(nzref(jmk) .ne. nzref(jm_mrg)) then
                  kduct=nzref(jm_mrg)
                  call zref_chng
                  if(iidiag .ge. 2) then
                     print *,'Changed ref depth for dup mode ',jm_mrg
                  endif
               endif
               kp=k_mrg
cc             call sheet_init(kp,0,iish,iish_ref)
               xkhrat=kp/kw0
               call r1r2_calc(kp,r1r2,3,1,jjfail)
               if(jjfail .gt. 0) return
               Lp=r1r2(1,4)
               dL_dkp=r1r2(2,4)
               dL_dwp=r1r2(3,4)
c: iimrg=3 ==> dfinc=.125, iimrg=4 ==> dfinc=.05:
               dfinc=.125/((iimrg-3)*3 + 1)
               goto 15
            elseif(iimrg .gt. 4) then
               print *,'Modes found to merge. Unable to resolve, '//
     .            'skipping: ',jmo,jm_mrg,faxbb(jfbb)
               nskip=nskip + 1
               goto 50
            endif
c
c: Keep track of minimum Im(k) in order to stop mode search due to rmin:
            kim_bb(jfbb2)=dmin1(max(0.d0,dimag(k)),kim_bb(jfbb2))
c: Flag for leaky modes:
            if(dreal(k) .lt. kcrmin .or. nsvmin .eq. nlay) iilk=1
c: Check if k to left of kremin (related to cphmax):
            if(dreal(k) .lt. kremin) then
               if(iidiag .ge. 1) print *,'BB cutoff due to kremin: ',
     .            jmo,k/kw0,f_hz
               goto 50
            endif
c: Check if Im(k) so high that mode is weak at shortest range of interest:
            if(iilk .eq. 1 .or. nrise .gt. 5) then
               kim_max=kim_bb(jfbb2) + dkim
            else
c: Safety factor for modes not leaky yet:
               kim_max=kim_bb(jfbb2) + 5.*dkim
            endif
            if(dimag(k) .gt. kim_max) then
               nhigh=nhigh + 1
c: Check if 2 modes in a row have been high & traj heading higher (see p.150):
c: Ignore this criterion since island modes can have any k-derivative:
ccx            if(nhigh .ge. 2 .and. dreal(r1r2(2,4)) .gt. 0.d0) then
               if(nhigh .ge. 2) then
                  if(iidiag .ge. 1) print *,'BB cutoff due to nhigh: ',
     .               jmo,jfbb,k/kw0,f_hz
                  goto 50
               endif
            else
               nhigh=0
            endif
c
            if(dimag(k) .gt. dimag(kn(nmode-1))) then
               nrise=nrise + 1
            else
               nrise=0
            endif
c
c: Numerical derivative test for dL/dw:
ctemp f_hz=f_hz + .01
ctemp call freq_chng
ctemp call r1r2_calc(k,r1r2x,3,1)
ctemp print *,'An,Nu = ',r1r2(3,4),(r1r2x(1,4)-r1r2(1,4))/(.01*2.*pie)
ctemp f_hz=f_hz - .01
ctemp call freq_chng
c
c: Valid eigenvalue found at new frequency:
            dW_dk=-2.*r1r2(1,1)*gamiref*r1r2(2,3)
c: Place mode functions in spare location at nm_fmax+1:
            call mode_fun(k,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
     .         psi,dpsi,exp_gbs,nm2)
c: Enter k into knbb, etc:
            vg=-r1r2(2,4)/r1r2(3,4)
            call bb_enter(k,r1r2(1,4),vg,r1r2(1,1),knbb,eig_bb,
     .         jfbb2,jmo,nfbb,nm_fmax,iish,iish_ref,iish_bb,nmbb)
c: Interpolate mode functions and compute field for interpolated freqs:
            if(nfinc .gt. 1) call bb_interp(phi,dphi,psi,dpsi,exp_gbs,
     .         phix,dphix,psix,dpsix,expx_gbs,nzsr,nm1,nm2,nm3,iisg,
     .         nfinc,knbb,tf,nfbb,jfbb,nm_fmax,jmo,faxbb,cfmin,iish,
     .         eig_bb,iift,kim_bb,phim_bb,iish_bb)
c
c: Compute field at all receivers at f_hz:
            jfbb=jfbb2
            call bb_field(k,phi,dphi,dpsi,exp_gbs,nm2,tf,jfbb,jmo)
            if(iift .eq. 1) call bb_write(jmo,jfbb,faxbb,k,kw0,iish)
c
c: Check if best duct in which to find mode has changed:
cc          if(nduct .gt. 1) then
cc             kductp=kduct
cc             call zduct_chng(phi,dphi,jmk,nm2,iichng,k)
cc             if(iichng .eq. 1) then
cc                call zref_chng
cc                call r1r2_calc(k,r1r2,3,1)
cc                if(iidiag .ge. 1) 
cc   .               print *,'changed duct: ',kductp,kduct,jmo,f_hz
cc             endif
cc          endif
c
            nm1=nm2
            nm2=nm2 + ii
            ii=-ii
cxx   print *,'nm1,nm2,ii = ',nm1,nm2,ii
c: Update "previous values" for next jump in frequency:
            dfinc=dfincx
            kp=k
            Lp=r1r2(1,4)
            dL_dkp=r1r2(2,4)
            dL_dwp=r1r2(3,4)
            if(dfinc .gt. 1.d0) dfinc=nint(dfinc)
            if(iisg .lt. 0) then
               if(jfbb .le. 1) kdone=1
            else
               if(jfbb .ge. nfbb) kdone=1
            endif
            nf_done=nf_done + 1
         enddo
50       continue
         ncmode=nctot - ncmode
         print *,'jm,nsave,fcut,nm = ',jmo,nsave-nsave0,faxbb(jfbb),
     .      nf_done,float(ncmode)/max(1,nf_done)
c         pause
         nmbbtot=nmbbtot + nf_done
c
c: Output mode characteristics if desired:
         if(iiout .ne. 0) then
            call bb_out(nfbb,nm_fmax,nzsr,jmo,knbb,phibb,dpsibb,nh_off)
         endif
         if(iimf .eq. 1) then
            nfout=nfbb-jfbb+1
            call bb_mfout(jmo,jfbb,nfout,phibb,r4mat1,r4mat2)
         elseif(iimf .eq. 2) then
cpln            call bb_mfout2(jmo,jfbb,nfout,nm_fmax,phibb,r4mat1,r4mat2)
         endif
      enddo
c
      return
      end
ccc
      subroutine bb_blmodes
c
c: Find the branch line modes from fmax to fmin in a brute force manner.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm,jfbb,nmbbtot,nm_fmax,ndup,ndup_max,jmo
c
      kduct=indx_duct(1)
      call zref_chng
      if(nsvmin .ne. nlay) then
         print *,'Bad nsvmin from bb_blm'
         return
      endif
      print *,'Enter nblm_max: '
      read(5,*) nblm_max

      nmbbtot=0
      do jfbb=nfbb,1,-1
         f_hz=faxbb(jfbb)
         call freq_chng
         nmode=0
         nm_put=0
c
         nblm=0
         write(6,*)'#4 mode_branch'
         call mode_branch(1,xkref,1,0,0,0,ndup,ndup_max)
         jmo=nmbb(jfbb)
c
         nm_fmax=max(nm_fmax,nmbb(jfbb) + nm_put)
         do jm=1,nm_put
            jmo=jmo + 1
            call bb_field(kn(jm),phi,dphi,dpsi,exp_gbs,jm,tf,jfbb,jm)
            call bb_enter(kn(jm),eig_char(1,jm),eig_char(4,jm),
     .         eig_char(5,jm),knbb,eig_bb,jfbb,jmo,nfbb,nm_put,
     .         iish,iish_ref,iish_bb,nmbb)
            if(iift .eq. 1) call bb_write(jmo,jfbb,faxbb,kn(jm),
     .         kw0,iish)
         enddo
         nmbbtot=nmbbtot + nm_put
         print *,'Done f = ',sngl(f_hz),nm_put,nmbb(jfbb),
     .      sngl(dimag(kn(nm_put))*8685.9)
      enddo
c
      print *,'nm_fmax = ',nm_fmax
      write(lusvp,120) nmbbtot
120   format('BRUTE FORCE BLM TOTAL # MODES = ',i8)
c
      return
      end
       subroutine bb_modes
c
c: Computes the broadband field from fmin to fmax in steps of df
c: using mode theory.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jm,nm_fmax,nm_fmin,nmbbtot,iibad,nmiss,ndied,nrose
      real*4 ncperm
c
      nctot=0
      nclast=0
      nmerge=0
      ncperm=0.e0
      nctot=0
      nclast=0
      nmode=0
      nm_fmax=0
      nsave=0
c
c: Check input variables, sizes of arrays, and initialize variables:
      call bb_init
c
c
      nmode=0
      call zmx_init
      f_hz=fmax
      call mode_find(1)
      if(jjfail.gt.0) return
      kim_bb(nfbb)=kim_min
      phim_bb(nfbb)=phi_mag_max
      nm_fmax=nmode
      print *,'fmax,nmode = ',fmax,nmode
c
c: Check if knbb and eig_bb are large enough:
      iibad=0
      call mem_lim(nmode*nfbb,NM_NF_MAX,MLINE,LML,'nmode*nfbb',10,
     .   'NM_NF_MAX',9,iibad,1)
      if(iibad .eq. 1) stop
c
c: Open mode characteristics file if desired:
      if(iiout .ne. 0) call bb_out_init(nm_fmax)
c
      call bb_ftinit
c
      do jm=1,nm_fmax
         kn_indx(jm)=jm
      enddo
      nmbbtot=nm_fmax
      call bb_mloop(0,nm_fmax,nmbbtot,-1,fmax,nfbb,ndied)
      if(isp(nlay) .eq. 0) then
         call bb_blmodes
      endif
c
c: Find modes at lowest frequency:
      f_hz=fmin
      call mode_find(0)
      if(jjfail.gt.0)return
      nm_fmin=nmode
c
      call bb_miss(knbb,nm_fmax,kn(1),nm_fmin,kn(nm_fmin+3),nmiss)
      if(nmiss .gt. 0) then
         call bb_mloop(nm_fmax,nmiss,nmbbtot,1,fmin,1,nrose)
         nm_fmax=nm_fmax + nmiss
      else
         nrose=0
      endif
      if(ndied .gt. nrose) then
         print *,'ndied>nrose: ',ndied,nrose
      endif
c
      call bb_done(nm_fmax)
c
      ncperm=float(nctot)/float(max(1,nmbbtot))
      print *,'nmbbtot,nctot,nsave= ',nmbbtot,nctot,nsave,ncperm
      if(nmerge .gt. 0) then
         print *,'# modes found to merge: ',nmerge
         print *,'# modes skipped due to merging = ',nskip
      endif
      write(lusvp,120) nmbbtot,nctot,ncperm,nsave
120   format('TOTAL # MODES = ',i8,'; # R1R2 CALCS = ',i8,
     .   '; #CALCS/MODE = ',f5.2,'; NSAVE = ',i6)
c
      return
      end
ccc
      subroutine bb_interp(phi,dphi,psi,dpsi,exp_gbs,phix,dphix,
     .   psix,dpsix,expx_gbs,nzsr,nm1,nm2,nm3,iisg,nfinc,knbb,tf,
     .   nfbb,jfbb,nmode,jmo,faxbb,cfmin,iish,eig_bb,iift,kim_bb,
     .   phim_bb,iish_bb)
c
      implicit none
      integer*4 nzsr,nm1,nm2,nm3,iisg,nfinc,nfbb,jfbb,nmode,jmo,
     .   iish(2,2)
      integer*4 jf1,jf2,jf,jz,jchar,iift,iish_bb(nfbb,nmode),jfx
      real*4 faxbb(nfbb),phim_bb(nfbb)
      real*8 kim_bb(nfbb),cfmin,kw0x,twpie
      complex*8 phi(nzsr,nm3),dphi(nzsr,nm3),psi(nzsr,nm3),
     .   dpsi(nzsr,nm3),exp_gbs(nzsr,nm3),phix(nzsr),dphix(nzsr),
     .   psix(nzsr),dpsix(nzsr),expx_gbs(nzsr),tf
      complex*16 knbb(nfbb,nmode),eig_bb(5,nfbb,nmode),fac
      data twpie/6.28318530717959/
c
c: Interpolate eigenvalues and mode characteristics from jf1 to jf2:
      jf1=jfbb
      jf2=jfbb + iisg*nfinc
      fac=(knbb(jf2,jmo) - knbb(jf1,jmo))/nfinc
      do jf=1,nfinc-1
         jfx=jf1+iisg*jf
         knbb(jfx,jmo)=knbb(jf1,jmo) + jf*fac
         iish_bb(jfx,jmo)=iish_bb(jf1,jmo)
      enddo
      do jchar=1,5
         fac=(eig_bb(jchar,jf2,jmo) - eig_bb(jchar,jf1,jmo))/nfinc
         do jf=1,nfinc-1
            jfx=jf1+iisg*jf
            eig_bb(jchar,jfx,jmo)=eig_bb(jchar,jf1,jmo) + jf*fac
         enddo
      enddo
c
c: Use nm3 column as temporary storage for interpolation factors:
      do jz=1,nzsr
         phi(jz,nm3)=(phi(jz,nm2)-phi(jz,nm1))/nfinc
         dphi(jz,nm3)=(dphi(jz,nm2)-dphi(jz,nm1))/nfinc
         psi(jz,nm3)=(psi(jz,nm2)-psi(jz,nm1))/nfinc
         dpsi(jz,nm3)=(dpsi(jz,nm2)-dpsi(jz,nm1))/nfinc
         exp_gbs(jz,nm3)=(exp_gbs(jz,nm2)-exp_gbs(jz,nm1))/nfinc
      enddo
      do jf=1,nfinc-1
         do jz=1,nzsr
            phix(jz)=phi(jz,nm1)+jf*phi(jz,nm3)
            dphix(jz)=dphi(jz,nm1)+jf*dphi(jz,nm3)
            psix(jz)=psi(jz,nm1)+jf*psi(jz,nm3)
            dpsix(jz)=dpsi(jz,nm1)+jf*dpsi(jz,nm3)
            expx_gbs(jz)=exp_gbs(jz,nm1)+jf*exp_gbs(jz,nm3)
         enddo
c: Compute field at all receivers at interpolated freq between f_hz0 & f_hz:
         jfx=jf1+iisg*jf
         call bb_field(knbb(jfx,jmo),phix,dphix,dpsix,expx_gbs,1,tf,
     .      jfx,jmo)
         if(iift .eq. 1) then
            kw0x=twpie*faxbb(jfx)/cfmin
            call bb_write(jmo,jfx,faxbb,knbb(jfx,jmo),kw0x,iish)
         endif
         kim_bb(jfx)=dmin1(kim_bb(jf1),kim_bb(jf2))
         phim_bb(jfx)=min(phim_bb(jf1),phim_bb(jf2))
      enddo
c
      return
      end
ccc
      subroutine bb_enter(k,lnrr,vg,R1,knbb,eig_bb,jfbb,jmo,nfbb,
     .   nmode,iish,iish_ref,iish_bb,nmbb)
c
c: Fills knbb and eig_bb with eigenvalue and mode characteristics.
c
      implicit none
      integer jfbb,nfbb,jmo,nmode,iish(2,2),iish_ref(2),
     .   iish_bb(nfbb,nmode),nmbb(nfbb)
      complex*16 k,lnrr(3),vg,R1,knbb(nfbb,nmode),eig_bb(5,nfbb,nmode)
c
      knbb(jfbb,jmo)=k
      eig_bb(1,jfbb,jmo)=lnrr(1)
      eig_bb(2,jfbb,jmo)=lnrr(2)
      eig_bb(3,jfbb,jmo)=lnrr(3)
      eig_bb(4,jfbb,jmo)=vg
      eig_bb(5,jfbb,jmo)=R1
      nmbb(jfbb)=jmo
      call iish_code(iish,iish_ref,iish_bb(jfbb,jmo),1)
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
ccc
      subroutine bb_out(nfbb,nmode,nzsr,jmo,knbb,phibb,dpsibb,nh_off)
c
      implicit none
      integer*4 nfbb,nmode,nzsr,jmo,nh_off,jf,jsr
      complex*16 knbb(nfbb,nmode)
      complex*8 phibb(nfbb,nzsr),dpsibb(nfbb,nzsr)
c
      write(33,rec=nh_off+jmo) (knbb(jf,jmo),jf=1,nfbb),
     .   ((phibb(jf,jsr),jf=1,nfbb),jsr=1,nzsr),
     .   ((dpsibb(jf,jsr),jf=1,nfbb),jsr=1,nzsr)
c
      return
      end
ccc
      subroutine bb_merge(k,knbb,nfbb,nmode,jfbb2,jmo,errdk2,iimrg,
     .   jm_mrg,k_mrg,nmerge,faxbb,iidiag)
c
c: Checks if this mode jmo has merged with any of previous modes found.
c
      implicit none
      integer*4 nfbb,nmode,jfbb2,jmo,jm_mrg,iimrg,nmerge,iidiag
      complex*16 knbb(nfbb,nmode),k,k_mrg
      real*8 magsq,kdif,errdk2
      real*4 faxbb(nfbb)
c
      do jm_mrg=jmo-1,1,-1
         kdif=magsq(knbb(jfbb2,jm_mrg) - k)
         if(kdif .le. errdk2) then
            if(iidiag .ge. 2) then
               print *,'Merged modes found: ',faxbb(jfbb2),jm_mrg,
     .            jmo,knbb(jfbb2,jm_mrg),k
               print *,'k = ',knbb(jfbb2+1,jm_mrg),knbb(jfbb2+1,jmo)
            endif
            iimrg=iimrg + 1
c: Mode found at previous (higher) frequency that was valid for sure:
            k_mrg=knbb(jfbb2+1,jm_mrg)
            if(iimrg .eq. 1) nmerge=nmerge + 1
            return
         endif
      enddo
      iimrg=0
c
      return
      end
c
      subroutine bb_fft_out
c
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 irec,jj0,jsrc,jrec,j,mx,nrd
      real*8 vgmax,tbuf
      real*4 rspace,freqs,dt1,sd
      complex*8 zzero
      character*16 trfext
      logical bintrf
      data zzero/(0.,0.)/
c
c      bintrf=.true.
      bintrf=.false.
c
c     rspace = range increment for receiver position
      if (nsrc.eq.1) then
         rspace=rkm(1)
      else
         rspace=(rkm(nsrc)-rkm(1))/float(nsrc-1)
      end if
      mx=nf1+nfbb-1
      dt1=1./fsbb
      freqs=(fmax-fmin)/2.+fmin
      sd=zsr(mzsrc(1))
      outroot(1:lout)='cprosim_out'
      outfile=outroot(1:lout)//'.dat'
      svp_title='TRF OUTPUT FROM SAGA'
c      write(6,*)'1st rec depth: ',zrec(1)
c      write(6,*)'last rec depth: ',zrec(nrec)
c      write(6,*)'Delta time: ',dt1
c      write(6,*)'Freq components: ',mx
c      write(6,*)'Center frequency: ',freqs
c      write(6,*)'Source depth: ',sd
c      write(6,*)'Range space: ',rspace
c      write(6,*)'Logical bintrf ',bintrf
c      write(6,*)'Before call to trfhead'
c      pause
      if(tilth) then
        call trfhead(outroot,svp_title,zrec(1),zrec(1),
     &    rkm(1),rspace,nfftbb,nf1,mx,dt1,freqs,sd,
     &    bintrf,1,nzs,nsrc)
      else
c        write(*,*)' calling trfhead nrec=',nrec
        call trfhead(outroot,svp_title,zrec(1),zrec(nrec),
     &    rkm(1),rspace,nfftbb,nf1,mx,dt1,freqs,sd,
     &    bintrf,nrec,nzs,nsrc)
      end if
c      write(6,*)'After call to trfhead'
c
      if(iifft .ne. 0) then
         tbuf=0
         vgmax= cfmin
c         do jsrc=1,nsrc
c            jj0=(jsrc-1)*nfbb*nrec
c            call bb_align(tf(jj0+1),nfbb,nrec,iifull)
c         enddo
c
         if (bintrf) then
            if(tilth) then
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec+(jsrc-1)*nfbb
                     write(luttrf)real(tf(j+jj0)),
     .                    -imag(tf(j+jj0))
                  end do
               end do
            else
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec
c     write(6,*)'No of nsrc: ',nsrc
                     do jrec=1,nrec
c     write(6,*)'No of nrec: ',nrec
c     write(6,*)'range: ',range(jsrc*jrec)
                        write(luttrf)real(tf(j+jj0)),-aimag(tf(j+jj0))
                        jj0=jj0 + nfbb
                     end do
                  end do
               end do
            end if
         else
            if(tilth) then
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec+(jsrc-1)*nfbb
                     write(luttrf,*)real(tf(j+jj0)),
     .                    -imag(tf(j+jj0))
                  end do
               end do
            else
               do j=1,nfbb
                  do jsrc=1,nsrc
                     jj0=(jsrc-1)*nfbb*nrec
                     do jrec=1,nrec
                        write(luttrf,*)real(tf(j+jj0)),-aimag(tf(j+jj0))
                        jj0=jj0 + nfbb
                     end do
                  end do
               end do
            end if
         end if
c     
         close(luttrf)
      endif
c
      return
      end
c
      subroutine trfhead(filenm,title,rd,rdlow,r0,rspace,
     &  nx,lx,mx,dt1,freqs,sd,bintrf,ir,is,
     &  nplots)
c
      implicit none
c
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
c
      logical bintrf
      integer*4 i,idummy,inttyp,isrow,icdr,msuft,nx,lx,mx,
     &          nplots,ir,is,nout,icnt
      real*4 dummy,adummy(10),omegim,r0,rd,rspace,rdlow,dt1,
     &       sd,freqs
      character*6 prognm
      character*8 fileid
      character*64 filenm
      character*64 title
      character signn
c >>> Dummy real*4 to force single prec trf file.
c
      omegim=0.
      if (bintrf) then
         if(is.gt.1) stop'in bb_fft_out'
c         open(luttrf,file=filenm(1:lout)//'.dat',
c     &        status='unknown',form='unformatted')
         fileid='PULSETRF'
         write(luttrf) fileid
         prognm='ORCA'
         write(luttrf) prognm
c         nout=1
c         write(luttrf) nout
         write(luttrf) int(1.)
c         icnt=1
c         write(luttrf) icnt
         write(luttrf) int(1.)
         write(luttrf) title(1:64)
         signn='+'
         write(luttrf) signn
c     CENTER FREQUENCY
         dummy=freqs
         write(luttrf) dummy
c     SOURCE DEPTH
         dummy=sd
         write(luttrf) dummy
c     UPPER MOST RECEIVER DEPTH
         adummy(1)=rd
c     LOWER MOST RECEIVER DEPTH
         adummy(2)=rdlow
c     IR=NO OF RECEIVERS BETWEEN R0 AND RDLOW
         write(luttrf) adummy(1),adummy(2),ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
c     do L = 1,abs(ir)
c     dummy=rdc(L)
c     WRITE(LUTTRF) dummy
c     end do
c     RDC(L) CONTAINS ALL THE RECEIVER DEPTHS
C     write(luttrf) (rdc(l),l=1,abs(ir))
C     END IF
C
C     R0= THE FIRST RANGE TO PLOT IN RANGE STACKED PLOT
         adummy(1)=r0
C     RSPACE= THE RANGE INCREMENT TO PLOT IN RANGE STACKED PLOT
         adummy(2)=rspace
c     WRITE R0, RSPACE AND THE NO OF PLOTS ASSOCIATED WITH IR
         WRITE(LUTTRF) adummy(1),adummy(2),nplots
C     DT=TIME SAMPLING
         dummy=dt1
C     NX=NO OF TIME SAMPLES (DENOTED NT IN OASES MANUAL) MAYBE?
C     LX=INDEX OF FIRST FREQUENCY COMPONENT (INT(FR1*DT))
C     MX=INDEX OF LAST FREQUENCY COMPONENT (INT(FR2*DT))
         write(luttrf) nx,lx,mx,dummy
         icdr=0
         write(luttrf) icdr
         dummy=omegim
         write(luttrf) dummy
C     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(luttrf) msuft
c         write(6,*) 'trfhead: msuft=',msuft
         isrow=1
         write(luttrf) isrow
         inttyp=1
         write(luttrf) inttyp
         idummy=0
         do 300 i=1,2
            write(luttrf) idummy
 300     continue
         dummy=0
         do 400 i=1,5
            write(luttrf) dummy
 400     continue
      else
c         write(6,*)'No of char in filenm: ',lout
c         write(6,'(a)')'Filename: ',filenm(1:lout)
c         pause
c         open(luttrf,file=filenm(1:lout)//'.dat',
c     &        status='unknown',form='formatted')
c         open(luttrf,file='pek.asc',
c     &        status='unknown',form='formatted')
         fileid='PULSETRF'
         prognm='ORCA'
         write(luttrf,'(1x,a)') fileid
         write(luttrf,'(1x,a)') prognm
c     WRITE(LUTTRF,*) NOUT
         write(luttrf,*) int(1.)
         icnt=0
         write(luttrf,*) int(1.)
         write(luttrf,'(1x,a)') title(1:64)
         signn='+'
         write(luttrf,'(1x,a)') signn
         write(luttrf,*) freqs
         write(luttrf,*) sd
c     IR=1
         write(luttrf,*) rd,rdlow,ir
C
C     MAY BE ADDED IN ORCA
C     IF (IR.LT.0) THEN
C     WRITE(LUTTRF,*) (RDC(L),L=1,ABS(IR))
C     END IF
C     NPLOTS=1
C
         write(luttrf,*) r0,rspace,nplots
         write(luttrf,*) nx,lx,mx,dt1
         icdr=0
         write(luttrf,*) icdr
         write(luttrf,*) omegim
c     ***  EXTRA FIELDS ADDED 891211 HS
         msuft=1
         write(luttrf,*) msuft
         isrow=1
         write(luttrf,*) isrow
         inttyp=1
         write(luttrf,*) inttyp
         idummy=0
         do 301 i=1,2
            write(luttrf,*) idummy
 301     continue
         dummy=0e0
         do 401 i=1,5
            write(luttrf,*) dummy
 401     continue
      end if
      return
      end
ccc
      subroutine bb_fft_out_original
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 irec,jj0,jsrc,jrec,j,lfft
      real*8 vgref,tbuf,tstart
      complex*8 zzero
      data zzero/(0.,0.)/
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
         xhbb(1)=-1.
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
            tstart=rng_sr(jsrc)/vgref - tbuf
            jj0=(jsrc-1)*nfbb*nrec
            call bb_align(tf(jj0+1),tstart,wbb,nfbb,nrec,iifull)
            xhbb(14)=tstart
            do jrec=1,nrec
               xhbb(4)=jsrc
c: Not compatible with new FFT format:
cc             xhbb(15)=zsr(mzrec(jrec))
cc             xhbb(17)=cp_sr(mzrec(jrec))
               xhbb(18)=rng_sr((jrec-1)*nsrc + jsrc)
               xhbb(19)=zsr(mzrec(jrec))
               irec=irec + 1
c: Note that CONJUGATE of transfer function is output to FFT file:
               write(10,rec=irec) (xhbb(j),j=1,20),
     .            (conjg(tf(j)),j=jj0+1,jj0+nfbb)
               jj0=jj0 + nfbb
            enddo
         enddo
         close(10)
      endif
c
      return
      end
ccc
      subroutine bb_mfout(jmo,jfbb,nfout,phibbx,phi_re,phi_im)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
c
      integer*4 jmo,jfbb,nfout,jrec,jz,jf,jfx,jfile
      complex*8 phibbx(nfbb,nzsr)
      real*4 phi_re(nfout,nrec),phi_im(nfout,nrec)
      character*10 mrsuf,misuf
      data mrsuf/'phire_m000'/,misuf/'phiim_m000'/
c
      do jrec=1,nrec
         jz=mzrec(jrec)
         do jf=jfbb,nfbb
            jfx=jf-jfbb+1
            phi_re(jfx,jrec)=real(phibbx(jf,jz))
            phi_im(jfx,jrec)=aimag(phibbx(jf,jz))
         enddo
      enddo
c
      write(mrsuf(8:10),'(i3.3)') jmo
      write(misuf(8:10),'(i3.3)') jmo
      jfile=0
      if(mod(jmo,10) .eq. 1) jfile=2
      call out_writex(outroot,lout,SUFX//mrsuf,11,phi_re,zrec,
     .   faxbb(jfbb),nrec,nfout,dlab,flab,z4,z4,z4,z4,jfile,
     .   'Re(phi) vs z,f',' ',' ',' ','f7.2','f7.2','f8.5',ncall)
      call out_writex(outroot,lout,SUFX//misuf,11,phi_im,zrec,
     .   faxbb(jfbb),nrec,nfout,dlab,flab,z4,z4,z4,z4,jfile,
     .   'Im(phi) vs z,f',' ',' ',' ','f7.2','f7.2','f8.5',ncall)
c
      return
      end
ccc
      subroutine bb_ftinit
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
c: Open file for the frequency trajectories of the modes:
      if(iift .eq. 1) then
         open(14,file=outroot(1:lout)//'_ftraj',status='unknown',
     .        form='formatted')
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
      subroutine bb_done(nm_fmax)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 nm_fmax,jsr,jf
c
c: Close files:
      if(iimt .ne. 0) close(21)
      if(iift .ne. 0) close(14)
      if(iiout .ne. 0) then
         close(33)
         open(33,file=outroot(1:lout)//'_bbeig',status='unknown',
     .      access='direct',
     .      recl=NRECL*lheadw)
         write(33,rec=1) fsbb,nfftbb,fmin,fmax,nfbb,nzsr,nm_fmax,
     .      rmin,rmax,sngl(cfmin),
     .      (zsr(jsr),jsr=1,nzsr),(nmbb(jf),jf=1,nfbb),iirx,
     .      (sngl(rho_sr(jsr)),jsr=1,nzsr)
cc       print *,'nmbb at end= ',(nmbb(jf),jf=1,nfbb)
cc       print *,'lheadw,lrecw = ',lheadw,lrecw,nh_off
         close(33)
      endif
c
c: Output FFT file:
c pln 020500       call bb_fft_out
c
c: Only output kni if not already done in disp_curv:
c pln 020500     if(iidc .ne. 2 .and. iidc .ne. 3) then
c         call out_writex(outroot,lout,SUFX//'kni',4,r4mat2,faxbb,xmode,
c     .        nfbb,nm_fmax,flab,mnlab,z4,z4,z4,z4,2,'Im(kn) vs Mode No',
c     .        ' ',' ','dB/km','f5.0','f6.1','f7.1',ncall)
c      endif
c
      return
      end
ccc
      subroutine bb_miss(knbbx,nm_fmax,knx,nm_fmin,kn_trk,nmiss)
c
c: Compares modes at fmin that were tracked from fmax with those that
c: were found at fmin using mode_find.  Finds any missing modes.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 nm_fmax,nm_fmin,nmiss,
     .   nm_trk,jm,jm_trk,nm1,jmf
      complex*16 knbbx(nfbb,nm_fmax),knx(nm_fmin),kn_trk(nm_fmax),
     .   k_found
      real*8 magsq
c
      nm_trk=0
      nm1=nm_fmax+1
      do jm=1,nm_fmax
c: Copy all non-zero eigenvalues to new list in knbb(:,nm1)
         if(knbbx(1,jm) .ne. (0.d0,0.d0)) then
            nm_trk=nm_trk + 1
            kn_trk(nm_trk)=knbbx(1,jm)
         endif
      enddo
c: Sort eigenvalues in knbb(:,nm1) for comparison with knx(1:nm_fmin:
      call hpsort_indx_c16(nm_trk,kn_trk,kn_indx)
      print *,'Sorted kn_trk: ',nm_trk,(kn_trk(jm),jm=1,nm_trk)
      print *,'kn at fmin: ',nm_fmin,(knx(jm),jm=1,nm_fmin)
c
      nmiss=0
      if(nm_trk .eq. nm_fmin) then
         print *,'BB tracking found same # modes as CW at fmin: ',
     .      nm_trk,nm_fmin
      elseif(nm_trk .gt. nm_fmin) then
         print *,'BB tracking found more modes than CW at fmin: ',
     .      nm_trk,nm_fmin
      else
         print *,'BB tracking found fewer modes than CW at fmin: ',
     .      nm_trk,nm_fmin
      endif
c
      jm_trk=1
      kn_trk(nm_trk+1)=(1.d100,0.d0)
      do jmf=1,nm_fmin
         k_found=knx(jmf)
10       continue
         if(magsq(k_found-kn_trk(jm_trk)) .gt. errdk2) then
            if(dreal(k_found) .gt. dreal(kn_trk(jm_trk))) then
               nmiss=nmiss + 1
               kn_indx(nmiss)=jmf
            else
               jm_trk=jm_trk + 1
               goto 10
            endif
         else
            jm_trk=jm_trk + 1
         endif
      enddo
      print *,'# missing modes: ',nmiss
      print *,'Missing modes: ',(knx(kn_indx(jm)),jm=1,nmiss)
c
      return
      end
      subroutine cub_fit_new(k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,delk)
c
c: This subroutine fits a cubic polynomial y(x)=c1 + c2*x + c3*x**2 + 
c: c4*x**3, where x=(k-k1)/(k2-k1), to the points (k1,y1) and (k2,y2) 
c: and the derivatives y1p and y2p at those points.
c
      implicit none
      real*8 k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,b1,b2,delk,y1px,y2px
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
c
      return
      end 
ccc
      subroutine cub_root_new(c1,c2,c3,c4,xhit,kcub)
c
c: Finds the root, xhit, between 0 and 1 of the
c: polynomial f(x)=c1 + c2*x + c3*x**2 + c4*x**4.
c
      implicit none
      integer*4 kcub
      real*8 c1,c2,c3,c4,xhit,a1,a2,a3,a1sq,q,r,
     .   dd,rdd,rq,arg,th,pietth,piefth,one_third
      data pietth/2.09439510239320/,piefth/4.18879020478639/,
     .   one_third/0.33333333333333/
c
      kcub=1
      a1=c3/c4
      a2=c2/c4
      a3=c1/c4
      a1sq=a1*a1
      q=(3.d0*a2 - a1sq)/9.
      r=(9.d0*a1*a2 - 27.d0*a3 - 2.d0*a1*a1sq)/54.d0
      dd=q**3 + r**2
      if(dd .ge. 0.d0) then
         rdd=sqrt(dd)
         xhit=dsign(1.d0,r+rdd)*abs(r+rdd)**(one_third) +
     .      dsign(1.d0,r-rdd)*abs(r-rdd)**(one_third) - a1*one_third
      else
         rq=sqrt(-q)
         arg=r/rq**3
         if((arg .lt. -1.d0) .or. (arg .gt. 1.d0)) then
            xhit=0.5d0
            kcub=0
            return
         endif
         th=acos(arg)/3.d0
         xhit=2.d0*rq*cos(th) - a1*one_third
         if((xhit .lt. 0.d0) .or. (xhit .gt. 1.d0)) then
            xhit=2.d0*rq*cos(th + pietth) - a1*one_third
            if((xhit .lt. 0.d0) .or. (xhit .gt. 1.d0)) then 
               xhit=2.d0*rq*cos(th + piefth) - a1*one_third
            endif
         endif
      endif
      if((xhit .le. 0.d0) .or. (xhit .ge. 1.d0)) then
         xhit=.5d0
         kcub=0
      endif
c     if(kcub .eq. 0) print *,'kcub=0 in cubroot' 
c
      return
      end 
ccc
      subroutine blug(geo,h,kt,bp,klay,nadd)
c
      implicit none
      include 'Parms_com'
      integer*4 kt,j,klay,nadd,kdo
      real*8 geo(2,5,NLMAX),h(NLMAX),bp(2),
     .   fac,rad,delc,beta,g0,bden,zmax,cmax,cdifmax
c
      if(kt .eq. 2) then
         g0=geo(2,1,klay)
         beta=bp(1)
         fac=geo(1,1,klay)*(1. + beta)
         rad=fac**2 + 2.*g0*fac*h(klay)
         geo(2,1,klay)=dsign(1.d0,fac)*dsqrt(rad) - beta*geo(1,1,klay)
      elseif(kt .eq. 3) then
c: kt=3 means that cp2 is given instead of gradient g for blug:
         beta=bp(1)
         delc=geo(2,1,klay) - geo(1,1,klay)
         g0=delc*(1. + delc/(2.*geo(1,1,klay)*(1. + beta)))/h(klay)
      elseif(kt .eq. 4) then
c: kt=4 means that cp2 is given 2nd, g instead of beta at end:
         g0=bp(1)
         delc=geo(2,1,klay) - geo(1,1,klay)
         bden=g0*h(klay) - delc
         if(bden .eq. 0.) then
            print *,'blug layer ',j,'illegal: cp1,cp2,g gives ',
     .         'a linear profile.'
            stop
         endif
         beta=delc**2/(2.*geo(1,1,klay)*bden) - 1.
      endif
      if(g0 .eq. 0. .or. rad .lt. 0.) then
         print *,'illegal blug profile: j,g0,rad = ',j,g0,rad,
     .      geo(1,1,klay),geo(2,1,klay),beta
         stop
      endif
c
      kdo=klay
      nadd=0
99    call c_dif_max(geo(1,1,kdo),geo(2,1,kdo),h(kdo),beta,
     .   cdifmax,cmax,zmax)
      if(cdifmax .gt. bp(2)) then
c: Create a new layer by splitting current layer:
         call add_layer(kdo,kdo+1,klay,zmax,cmax,geo,h)
         nadd=nadd + 1
      else
         kdo=kdo + 1
      endif
      if(kdo .le. klay) goto 99
c
      return
      end
ccc
      subroutine c_dif_max(c1,c2,h,beta,cdifmax,cmax,zmax)
c
      implicit none
      include 'Parms_com'
      integer*4 j,nvalx
      real*8 c1,c2,h,beta,g0,delc,fac,sg,g0fac2,fac3,cdifmax,zpt,cai,
     .   cblug,rad,cmax,zmax,cdif,facsq,alph
c
      delc=c2 - c1
      g0=delc*(1. + delc/(2.*c1*(1. + beta)))/h
      fac=c1*(1. + beta)
      sg=dsign(1.d0,fac)
      facsq=fac**2
      g0fac2=2.*g0*fac
      fac3=beta*c1
      alph=(1./c2**2 - 1./c1**2)/h
c: Find maximum difference between BLUG profile and 1/c**2 profile:
      nvalx=19
      cdifmax=-1.
      do j=1,nvalx
         zpt=j*h/(nvalx+1)
         cai=c1/dsqrt(1.d0 + alph*c1*c1*zpt)
         rad=facsq + g0fac2*zpt
         cblug=sg*dsqrt(rad) - fac3
         cdif=abs(cai - cblug)
         if(cdif .gt. cdifmax) then
            zmax=zpt
            cmax=cblug
            cdifmax=cdif
         endif
      enddo
c
      return
      end
ccc
      subroutine add_layer(k,kp1,ktot,zpt,cpt,geo,h)
c
      implicit none
      integer*4 k,kp1,ktot,j,j1,j2,ji
      real*8 zpt,cpt,geo(2,5,kp1),h(kp1),fac,hold
c
c: Copy any layers below to one layer down:
      do j=ktot,kp1,-1
         j1=j+1
         h(j1)=h(j)
         do ji=1,2
            do j2=1,5
               geo(ji,j2,j1)=geo(ji,j2,j)
            enddo
         enddo
      enddo
c
      hold=h(k)
c: Make thickness of new layer the leftover amount:
      h(kp1)=h(k) - zpt
c: Make thickness of old layer the depth of max sound speed diff:
      h(k)=zpt
c: Make bottom of new layer same as bottom of old layer:
      do j2=1,5
         geo(2,j2,kp1)=geo(2,j2,k)
      enddo
c
c: Enter blug value of sound speed into profile:
      geo(2,1,k)=cpt
c: Linearly interpolate (except cp) to get bottom of old layer:
ccc   fac=zpt/hold
c: Interpolate other parameters in same non-linear way as BLUG for cp:
      fac=(cpt - geo(1,1,k))/(geo(2,1,kp1) - geo(1,1,k))
      do j2=2,5
         geo(2,j2,k)=geo(1,j2,k) + fac*(geo(2,j2,k)-geo(1,j2,k))
      enddo
c
c: Make top of new layer same as bottom of old layer:
      do j2=1,5
         geo(1,j2,kp1)=geo(2,j2,k)
      enddo
c
c: Increment # of layers:
      ktot=ktot + 1
c
      return
      end
      subroutine contour_find(k_jump,iduct,k0,rr0,iidiff)
c
c: Checks if previously found mode at k_jump is now an island mode
c: relative to current reference depth.  If so, finds |R1R2|=1 contour 
c: by following path of steepest descent to saddle point and then
c: the constant arg(R1R2) contour to the |R1R2|=1 contour.
c: iidiff=1 is flag that a different main contour has been found.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iduct,iidiff,ndup,nc00,nm1,ntry1
      complex*16 k_jump,k0,rr0(3,4),r1r2(3,4),knd,kstart,dk_rad
      real*8 mag_des0,emagmax0,ephmax0,magsq
c
      ntry1=0
      nm1=0
      mag_des0=0.d0
      emagmax0=0.1d0
      ephmax0=0.19635d0
      iimst=1
      k0=k_jump
      call sheet_init(k0,1,iish,iish_ref)
c: Need to call r1r2_calc again since derivatives will change in
c: different duct:
      call r1r2_calc(k0,rr0,2,0,jjfail)
      if(jjfail .gt. 0) return
c
10    continue
      if(iimt .eq. 1) call mode_traj(k0,rr0,0)
      nm1=nmode + 1
      nc00=nctot
c: Use iiccw=-2 to go to right in k-plane:
      iiccw=-2
      nhigh=0
      call eig_findm(k0,rr0,r1r2,nm1)
      if(jjfail.eq.1) return
      ndup=0
      knd=kn(nmode) - k_jump
      if(magsq(knd) .lt. errdk2) then
c: Last mode found in first duct is now an island mode in this duct.
c: This means the mode branches are separate, and we now need to find the 
c: main branch for this duct:
         call sdp_find(rr0,0,dk_rad)
25       if(iimt .eq. 1) call mode_traj(k_sdp,rr0,0)
         kstart=k_sdp
c: Find saddle point, move beyond and find unity magnitude contour:
         call traj_sdp(k_sdp,ln_sdp,dln_sdp,k0,rr0,kstart,dk_rad)
         if(jjfail.gt.0) return
         if(iifail .ne. 0) then
            call sdp_find(rr0,1,dk_rad)
            if(deepest .ne. 0.) then
               iifail=0
               goto 25
            else
               if(iiwrite .gt. 0)
     .          print *,'traj_sdp unsuccessful '
               return
            endif
         endif
         iidiff=1
      else
         if(dimag(k0) .gt. kim_min + dkim) then
            iidiff=0
         else
            if(iiblm .ne. 0) then
               iidiff=0
               if(iidiag .ge. 1) print *,'BLM found from contour_find'
            else
c: Check if new mode has already been found:
cpln               write(6,*)'Enter from contour_find'
cpln               write(6,*)'nm_put,ndup,iduct: ',nm_put,ndup,iduct
               call duct_dupl(nm_put,ndup,r1r2,phi,dphi,iduct)
               if(jjfail.eq.1) return
               if(ndup .eq. 1) then
c: New mode has already been found, which means the main branch for this
c: duct is the same as the one for the first duct.  We do not need to
c: search this branch any more.
ccp            print *,'Main branches found to be the same'
                  iidiff=0
               else
c: New mode has not been found, so we should continue searching the 
c: branch:
                  if(iiwrite .gt. 0)
     .            print *,'Main branches different???  Tell EKW.'
                  iidiff=0
c: EKW FIX 6-2-98. Don't include mode that was found.
                  nm_put=nm_put - 1
               endif
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine sdp_find(rr0,ii,dk_rad)
c
c: For ii=0, finds local minima in magsq(dln) on previous island-mode 
c: contour.  Uses deepest minimum as first guess.
c: For ii=1, uses other local minima found.
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 ii,j,jmin,jopp,jdl
      complex*16 rr0(3,4),dk_rad
      real*8 magsq,dln_mag,dln_magp,dln_magn
c
      deepest=0.
      if(ii .eq. 0) then
         nmin_dln=0
         dln_magp=magsq(dln_cont(npt))
         dln_magn=magsq(dln_cont(1))
         k_cont(npt+1)=k_cont(1)
         ln_cont(npt+1)=ln_cont(1)
         dln_cont(npt+1)=dln_cont(1)
         do j=1,npt
            dln_mag=dln_magn
            dln_magn=magsq(dln_cont(j+1))
ccp         print *,'dln_mag = ',j,dln_mag
            if(dln_mag .lt. dln_magp .and. dln_mag .lt. dln_magn) then
               nmin_dln=nmin_dln + 1
               jmin_dln(nmin_dln)=j
c: Measure of how deep the minimum was:
               dln_deep(nmin_dln)=dsqrt((dln_magp-dln_mag)*
     .            (dln_magn-dln_mag))/dln_mag
               if(dln_deep(nmin_dln) .gt. deepest) then
                  jdl=nmin_dln
                  deepest=dln_deep(nmin_dln)
               endif
c              print *,'dln_mag = ',j,nmin_dln,dln_mag
            endif
            dln_magp=dln_mag
         enddo
         if(nmin_dln .eq. 1) then
c: Add another point to try (on opposite side of island) in case first one 
c: did not work:
            nmin_dln=2
            jmin_dln(2)=jmin_dln(1) + npt/2
            if(jmin_dln(2) .gt. npt) jmin_dln(2)=jmin_dln(2)-npt
            dln_deep(2)=0.
         endif
      else
         do j=1,nmin_dln
            if(dln_deep(j) .gt. deepest) then
               jdl=j
               deepest=dln_deep(j)
            endif
         enddo
         if(deepest .eq. 0.) return
      endif
c
      jmin=jmin_dln(jdl)
      k_sdp=k_cont(jmin)
      ln_sdp=ln_cont(jmin)
      dln_sdp=dln_cont(jmin)
      rr0(1,4)=ln_sdp
      rr0(2,4)=dln_sdp
c: Find direction radially away from center of island contour:
      jopp=jmin + npt/2
      if(jopp .gt. npt) jopp=jopp-npt
      dk_rad=k_sdp - k_cont(jopp)
c: Set deep negative so that this point not chosen again:
      dln_deep(jdl)=-1.
c
      return
      end
ccc
      subroutine traj_sdp(klast,lnlast,dlnlast,k,r1r2,kstart,dk_rad)
c
c: Moves to saddle point by fitting two points and first derivatives
c: to a cubic.  Then continues in the same direction along path of 
c: steepest descent or ascent to the unity magnitude contour.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 ntry1,ntry2,kinc,kdec,ii_usek,j,k_count
      complex*16 klast,lnlast,dlnlast,k,r1r2(3,4),k2,ln2,dln2,kstart,
     .   dk_rad,dk_spt,lnrrdes,delk
      real*8 magsq,dot,maglast,mag_spt,ph_des0,lnmag,dmag_des,
     .   emagmax0,dmag_left0,dk_lin,drat,mag0,delk_magsq,mag_max
c
c: Try not moving along phase contour:
caa   ph_des0=dimag(lnlast)
c: Compute correct sign of delta mag based on radial from island center:
caa   dmag_left0=dsign(20.d0,dreal(dk_rad*dlnlast))
caa   emagmax0=0.1d0
c: Move on constant phase contour until derivative reaches local min:
cc    call traj_ph_sdp(ph_des0,dmag_left0,klast,lnlast,dlnlast,k,r1r2,
cc   .   emagmax0,0.19635d0)
caa   emagmax0=0.02d0
caa   call traj_ph_sdp(ph_des0,dmag_left0,klast,lnlast,dlnlast,k,r1r2,
caa  .   emagmax0,0.019635d0)
caa   klast=k
caa   lnlast=r1r2(1,4)
caa   dlnlast=r1r2(2,4)
c
      iifail=0
      ntry1=0
      maglast=magsq(dlnlast)
      mag0=maglast
      mag_max=1.d-8*maglast
      dk_lin=0.05d0/dsqrt(maglast)
      ii_usek=0
c
      if(iiwrite .gt. 0)
     .  print *,'Starting point from circle: ',dreal(klast)/kw0,
     .   dimag(klast)*8685.9
10    continue
c: Choose neighboring point with which to compute cubic fit:
      kinc=0
      kdec=0
      k_count=0
      if(ii_usek .eq. 0) then
5        k2=klast + dk_lin
         call r1r2_calc(k2,r1r2,2,0,jjfail)
         k_count=k_count+1
         if(k_count .gt. 100)
     .     jjfail=1
         if(jjfail .gt. 0) return
         ln2=r1r2(1,4)
         dln2=r1r2(2,4)
c: Make sure k2 not too close to k_spt:
         if(kinc .ne. 1 .or. kdec .ne. 1) then
            drat=magsq(dln2-dlnlast)/mag0
            if(drat .lt. 1.d-7) then
               dk_lin=10.d0*dk_lin
               kinc=1
               goto 5
            elseif(drat .gt. 1.d-2) then
               dk_lin=0.1d0*dk_lin
               kdec=1
               goto 5
            endif
         endif
      endif
c
      ntry1=ntry1 + 1
      if(ntry1 .gt. 50) then
         if(iiwrite .gt. 0)
     .    print *,'Failure in traj_sdp: ntry1 = ',ntry1
         iifail=1
         return
      endif
c: Find delk to get to point where derivative of L is zero (saddle point):
      call lnk_cub_ext(klast,k2,lnlast,ln2,dlnlast,dln2,delk)
      k_spt=klast + delk
      ntry2=0
c
20    continue
      ntry2=ntry2 + 1
      if(ntry2 .gt. 4) then
         if(iiwrite .gt. 0) then
            print *,'Failure in traj_sdp: ntry2 = ',ntry2,
     .           maglast,mag_spt,
     .           dreal(k_spt)/kw0,dimag(k_spt)*8685.9
            print *,'Num: ',(ln_spt-lnlast)/(2.*delk),
     .           dln_spt,dlnlast
         end if
         iifail=1
         return
      endif
      call r1r2_calc(k_spt,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      ln_spt=r1r2(1,4)
      dln_spt=r1r2(2,4)
      mag_spt=magsq(dln_spt)
c
      delk_magsq=magsq(delk)
c: If delk very small or derivative very small, we are close enough
c: to saddle point:
      if(delk_magsq .le. errdk2 .or. mag_spt .lt. mag_max) goto 50
c
      if(mag_spt .lt. maglast) then
         if(iiwrite .gt. 0)
     .   print *,'Lower mag of deriv: ',dreal(klast)/kw0,
     .      dimag(klast)*8685.9
c: Magnitude of derivative decreased:
         ii_usek=0
         if(delk_magsq .lt. dk_lin*dk_lin) then
            ii_usek=1
            k2=klast
            ln2=lnlast
            dln2=dlnlast
         endif
         klast=k_spt
         lnlast=ln_spt
         dlnlast=dln_spt
         maglast=mag_spt
         if(iimt .ne. 0) call mode_traj(k_spt,r1r2,0)
         goto 10
      else
c: If derivative increased, decrease delk:
         if(iimt .ne. 0) call mode_traj(k_spt,r1r2,-1)
         delk=0.5d0*delk
         k_spt=klast + delk
         goto 20
      endif
c
50    continue
c: Now we are at the saddle point k_spt:
      if(magsq(k_spt) .lt. errdk2) then
c: Check if at origin, where derivative is zero:
         iifail=1
         return
      endif
c
c: Fit cubic at saddle point k_spt and nearby point k2:
c
      if(iiwrite .gt. 0)
     .  print *,'Saddle pt: ',dreal(k_spt)/kw0,dimag(k_spt)*8685.9
      lnmag=dreal(ln_spt)
c: Determine sign of desired step in magnitude (go toward |I|=1):
      dmag_des=sign(min(0.1d0,dabs(lnmag)),-lnmag)
      ph_des0=dimag(ln_spt)
      lnrrdes=dcmplx(lnmag + dmag_des,ph_des0)
      dk_spt=k_spt - kstart
      dk_lin=1.d-5*cdabs(dk_spt)
c: Choose k2 near k_spt to make cubic fit from two neighboring points:
25    k2=k_spt + dk_lin
      call r1r2_calc(k2,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      ln2=r1r2(1,4)
      dln2=r1r2(2,4)
c: Make sure k2 not too close to k_spt:
      if(magsq(ln2-ln_spt)/magsq(ln_spt) .lt. 1.d-7) then
         dk_lin=10.d0*dk_lin
         goto 25
      endif
c
c: Fit ln and ln at k_spt and k2 to a cubic and solve for delk that
c: makes function equal to lnrrdes:
      call lnk_cub(k_spt,k2,ln_spt,ln2,dln_spt,dln2,lnrrdes,delk)
c: Near zero of cubic, we can go +delk or -delk.  Determine which one:
      if(dot(delk,dk_spt) .lt. 0.d0) then
         k_spt=k_spt - delk
      else
         k_spt=k_spt + delk
      endif
c: Now step past saddle point and follow constant phase contour to 
c: the unity magnitude contour:
      call r1r2_calc(k_spt,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      ln_spt=r1r2(1,4)
      dln_spt=r1r2(2,4)
      if(iimt .ne. 0) call mode_traj(k_spt,r1r2,0)
c
      emagmax0=min(0.1d0*dabs(dreal(ln_spt)),0.1d0)
      dmag_left0=-dreal(ln_spt)
      ph_des0=dimag(ln_spt)
      call traj_phase(ph_des0,dmag_left0,k_spt,ln_spt,dln_spt,k,
     .   r1r2,emagmax0,0.19635d0)
      if(jjfail.gt.0) return
      if(iiwrite .gt. 0)
     .  print *,'After traj_phase: ',dreal(k)/kw0,dimag(k)*8685.9
c
      return
      end
ccc
      subroutine lnk_cub_ext(k1,k2,lnk1,lnk2,lnk1p,lnk2p,delk)
c
c: Fits a cubic to the function ln(k), given two points (k1,lnk1) and 
c: (k2,lnk2) and the derivatives lnk1p and lnk2p.  Solves the cubic to 
c: find the point kext=k1 + delk nearest k1 where (d/dk) ln(k) = 0
c
      implicit none
      complex*16 k1,k2,lnk1,lnk2,lnk1p,lnk2p,delk,c1,c2,c3,c4,
     .   dk,ax,bx,cx,q,kext1,kext2,kext
cce  .   a2,b2,fpp1,fpp2,dk_re_min1,dk_re_min2,kext1_norm,kext2_norm
      real*8 magsq
c
c: Fit cubic:
      call lnk_cub_fit(k1,k2,lnk1,lnk2,lnk1p,lnk2p,c1,c2,c3,c4,dk)
c: Take derivative to get quadratic ax*x^2 + bx*x + c:
      ax=3.*c4
      bx=2.*c3
      cx=c2
c: Find roots of quadratic, which are extrema of cubic:
      q=-.5d0*(bx + sign(1.d0,dreal(bx))*cdsqrt(bx*bx - 4.d0*ax*cx))
      kext1=q/ax
      kext2=cx/q
c: Take second derivative (a2*x + bx) and evaluate at extrema to find if 
c: local max or min:
cce   a2=2.*ax
cce   b2=bx
cce   fpp1=a2*kext1+b2
cce   fpp2=a2*kext2+b2
cce   dk_re_min1=cdsqrt(cdabs(fpp1)/fpp1)
cce   dk_re_min2=cdsqrt(cdabs(fpp2)/fpp2)
c: Normalize vectors from midpoint of k1,k2 to extrema:
cce   kext1_norm=dcmplx(real(kext1) - 0.5d0,dimag(kext1))
cce   kext1_norm=kext1_norm/cdabs(kext1_norm)
cce   kext2_norm=dcmplx(real(kext2) - 0.5d0,dimag(kext2))
cce   kext2_norm=kext2_norm/cdabs(kext2_norm)
cce   dot1=dot(dk_re_min1,kext1_norm)
cce   dot2=dot(dk_re_min2,kext2_norm)
cce   print *,'cubic extrema: ',kext1,kext2,dot1,dot2
c     if(abs(dot1) .gt. abs(dot2)) then
c        kext=kext1
c     print *,'kext1 chosen: ',magsq(kext1),magsq(kext2),dot1,dot2
c     else
c        kext=kext2
c     print *,'kext2 chosen: ',magsq(kext2),magsq(kext1),dot2,dot1
c     endif
ccp   print *,'magsq(kext) = ',magsq(kext1-.5),magsq(kext2-.5)
      if(magsq(kext1-.5d0) .lt. magsq(kext2-.5d0)) then
         kext=kext1
      else
         kext=kext2
      endif
c: Transform from (0,1) to (k1,k2):
      delk=kext*dk
cc    kext=k1 + kext*dk
c
      return
      end
ccc
      subroutine lnk_cub_spt(k1,k2,lnk1,lnk2,lnk1p,lnk2p,delk,dln_re)
c
c: Fits a cubic to the function ln(k), given two points (k1,lnk1) and 
c: (k2,lnk2) and the derivatives lnk1p and lnk2p.  Solves the cubic to 
c: find the point kext=k1 + delk nearest k1 where (d/dk) ln(k) = 0
c
      implicit none
      complex*16 k1,k2,lnk1,lnk2,lnk1p,lnk2p,delk,c1,c2,c3,c4,
     .   dk,ax,bx,cx,q,kext1,kext2,kext,a2,b2,
     .   fpp1,fpp2,dk_re_min1,dk_re_min2,kext1_norm,kext2_norm
      real*8 dln_re,magsq,dot,dot1,dot2
c
c: Fit cubic:
      call lnk_cub_fit(k1,k2,lnk1,lnk2,lnk1p,lnk2p,c1,c2,c3,c4,dk)
c: Take derivative to get quadratic ax*x^2 + bx*x + c:
      ax=3.*c4
      bx=2.*c3
      cx=c2
c: Find roots of quadratic, which are extrema of cubic:
      q=-.5d0*(bx + sign(1.d0,dreal(bx))*cdsqrt(bx*bx - 4.d0*ax*cx))
      kext1=q/ax
      kext2=cx/q
c: Take second derivative (a2*x + bx) and evaluate at extrema to find if 
c: local max or min:
      a2=2.*ax
      b2=bx
      fpp1=a2*kext1+b2
      fpp2=a2*kext2+b2
      dk_re_min1=cdsqrt(cdabs(fpp1)/fpp1)
      dk_re_min2=cdsqrt(cdabs(fpp2)/fpp2)
c: Normalize vectors from midpoint of k1,k2 to extrema:
      kext1_norm=dcmplx(real(kext1) - 0.5d0,dimag(kext1))
      kext1_norm=kext1_norm/cdabs(kext1_norm)
      kext2_norm=dcmplx(real(kext2) - 0.5d0,dimag(kext2))
      kext2_norm=kext2_norm/cdabs(kext2_norm)
      dot1=dot(dk_re_min1,kext1_norm)
      dot2=dot(dk_re_min2,kext2_norm)
ccp   print *,'cubic extrema: ',kext1,kext2,dot1,dot2
c     if(abs(dot1) .gt. abs(dot2)) then
c        kext=kext1
c     print *,'kext1 chosen: ',magsq(kext1),magsq(kext2),dot1,dot2
c     else
c        kext=kext2
c     print *,'kext2 chosen: ',magsq(kext2),magsq(kext1),dot2,dot1
c     endif
ccp   print *,'magsq(kext) = ',magsq(kext1-.5),magsq(kext2-.5)
      if(magsq(kext1-.5d0) .lt. magsq(kext2-.5d0)) then
         kext=kext1
      else
         kext=kext2
      endif
c: Transform from (0,1) to (k1,k2):
      delk=kext*dk
cc    kext=k1 + kext*dk
c
      return
      end
ccc
      subroutine traj_ph_sdp(ph_des0,dmag_left0,klast,lnlast,
     .   dlnlast,k,r1r2,emagmax0,ephmax0)
c
c: Moves along line of constant phase of the function ln(r1*r2) in the
c: complex k-plane.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 ntry
      complex*16 klast,lnlast,dlnlast,lnrrdes,k,delk,r1r2(3,4),
     .   lnr1r2,dlnr1r2
      real*8 ph_des0,dmag_left0,emagmax0,ephmax0,emagmax,
     .   dmag_des,mag_des,ph_des,dph_des,mag_err,dmag_diff,
     .   eph_rat,emag_rat,mag_step_max,ph_err,dmag_left,ph_fac,
     .   mag_dlnp,mag_dln,magsq
c
      ph_fac=1.d0
      dmag_left=dmag_left0
      mag_step_max=1.d0
      mag_dlnp=magsq(dlnlast)
10    dmag_des=dmag_left/max(1,nint(dabs(dmag_left)/mag_step))
      ph_fac=min(2.d0*ph_fac,1.d0)
      dph_des=ph_fac*(ph_des0 - dimag(lnlast))
      ntry=0
      ph_des=ph_des0
20    mag_des=dreal(lnlast) + dmag_des
      ph_des=dimag(lnlast) + dph_des
      lnrrdes=dcmplx(mag_des,ph_des)
c
      if(ntry .le. 3) then
c: Compute next guess for k based on last value and first derivative:
         delk=(lnrrdes-lnlast)/dlnlast
      else
         call lnk_cub(klast,k,lnlast,lnr1r2,dlnlast,dlnr1r2,
     .      lnrrdes,delk)
cc       print *,'cubic fit in traj_phase: ',klast,k,lnlast,lnr1r2,
cc   .      dlnlast,dlnr1r2,delk,(lnrrdes-lnlast)/dlnlast
      endif
      k=klast + delk
      ntry=ntry + 1
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnr1r2=r1r2(1,4)
      dlnr1r2=r1r2(2,4)
c
      if(ntry .gt. 10) then
         if(iiwrite .gt. 0) then
            print *,'ntry>10 in traj_phase: ,ntry,k,k/kw,iish = ',
     .      ntry,k,k/kw,iish
            print *,'k,klast,delk = ',k/kw,klast/kw,delk/kw,lnrrdes,
     .      lnlast,dlnlast
            print *,'xkbp,iish = ',xkbp(1,1)/kw,iish
            print *,'num deriv = ',(r1r2(1,4)-lnlast)/(k-klast),
     .      dlnlast,r1r2(2,4)
         end if
         if(ntry .gt. 11) then
            iifail=1
cc          call stop_run
            return
         endif
      endif
c
      ph_err=ph_des - dimag(lnr1r2)
      eph_rat=dabs(ph_err/ephmax0)
c: Check for bad phase:
      if(eph_rat .gt. 1.d0) then
         dmag_des=dmag_des/2.d0
         dph_des=dph_des/2.d0
         mag_step=mag_step/2.d0
         ph_fac=ph_fac/2.d0
         if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
         call xkh_backup
         goto 20
      endif
c
c: Make sure max mag error emagmax is less than desired mag jump:
      emagmax=min(emagmax0,.75d0*dabs(dmag_des))
c: Check for phase overshoot:
      mag_err=dabs(mag_des - dreal(lnr1r2))
c: Don't allow phase to be off in other direction either:
      emag_rat=mag_err/emagmax
      if(emag_rat .gt. 1.d0) then
         dmag_des=dmag_des/2.d0
         dph_des=dph_des/2.d0
         mag_step=mag_step/2.d0
         ph_fac=ph_fac/2.d0
         if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
         call xkh_backup
         goto 20
      else
         if(emag_rat .lt. 0.25d0 .and. eph_rat .lt. 0.25d0) then
            mag_step=dmin1(2.0d0*mag_step,mag_step_max)
         elseif(emag_rat .lt. 0.5d0 .and. eph_rat .lt. 0.5d0) then
            mag_step=dmin1(1.5d0*mag_step,mag_step_max)
         endif
      endif
c
      dmag_diff=dreal(lnr1r2) - dreal(lnlast)
c: Update how much phase left to go:
      dmag_left=dmag_left - dmag_diff
c
      mag_dln=magsq(dlnr1r2)
      if(mag_dln .gt. mag_dlnp) then
         k=klast
         r1r2(1,4)=lnlast
         r1r2(2,4)=dlnlast
ccp   print *,'traj_ph_sdp ended for dln: ',mag_dlnp,mag_dln
         return
      endif
      mag_dlnp=mag_dln
c
      klast=k
      lnlast=lnr1r2
      dlnlast=dlnr1r2
      if(iimt .ne. 0) call mode_traj(k,r1r2,0)
c
c: Check if close enough to call eig_final:
      if(dabs(dmag_left) .gt. emagmax0) goto 10
c
ccp   print *,'traj_ph_sdp ended for dmag_left: ',mag_dlnp,mag_dln
      return
      end
ccc
      subroutine contour_find2(k_st,rr_st,iidiff)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iidiff
      complex*16 k_st,rr_st(3,4),k,r1r2(3,4),lnr1r2,dlnr1r2,kmid
      real*8 ph_des0,dmag_left0
c
      iidiff=1
      if(iabs(iiblm) .eq. 1) then
         kmid=0.5d0*xk(1,nlay)
         if(rmin .lt. 999.) then
cc          k=dcmplx(0.d0,dkim)
            k=kmid
            k=kmid + dcmplx(0.d0,.5d0*dkim)
         else
            k=dcmplx(kremin,0.d0)
         endif
      elseif(iabs(iiblm) .eq. 2) then
c: If contour found BLM associated with shear waves, find |R1R2|=1 contour
c: *between* p- and s-wave branch points:
         kmid=0.5d0*(xk(1,nlay) + xb(1,nlay))
         if(rmin .lt. 999.) then
            k=kmid
         else
            if(dreal(kmid) .gt. kremin) then
               k=kmid
            else
               k=dcmplx(kremin,0.d0)
            endif
         endif
      endif
c
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnr1r2=r1r2(1,4)
      dlnr1r2=r1r2(2,4)
      ph_des0=dimag(lnr1r2)
      dmag_left0=-dreal(lnr1r2)
      if(iimt .eq. 1) call mode_traj(k,r1r2,0)
c
c: Follow constant arg(R1R2) to |R1R2|=1 contour:
      call traj_phase(ph_des0,dmag_left0,k,lnr1r2,dlnr1r2,k_st,
     .   rr_st,0.1d0,0.19635d0)
      if(jjfail.gt.0) return
      if(iifail .ne. 0) then
         print *,'Failure in traj_phase from contour_find2'
         iidiff=0
      endif
c
      return
      end
      subroutine cw_modes(iiwrt)
c
c: Performs CW mode computations
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 iiwrt,jfcw,jm,lsuf,lsuf0,lsuf2,jzs,
     .   jj,iibad,j
      real*8 vph,vg,db_km
      real*4 cpsec0,cpsec,etime,time(2),pierad
      complex*8 pl0
      character*64 hdf_suf
      data pierad/0.0174532925/
c
      nctot=0
      nclast=0
      nmode=0
      iibad=0
      call uni_space(nfcw,fcw,1.e0)
      call sr_geom(rng_sr,iabs(nsrc),iabs(nrec))
      if(iitl .ne. 0) then
         call mem_lim(nrec*nsrc,NTLMAX,MLINE,LML,'nrec*nsrc',
     .      9,'NTLMAX',7,iibad,1)
      endif
c
      if(iiwrt .eq. 1) then
         open(7,file=outroot(1:lout)//'_modes',status='unknown',
     .        form='formatted')
         if(iimt .ne. 0) then
            open(21,file=outroot(1:lout)//'_mtraj',status='unknown',
     .        form='formatted')
         endif
      endif
c
      call zmx_init
c
      do jfcw=1,nfcw
         call suffix(' ',0,jfcw,nfcw,SUFX//'f',2,hdf_suf,lsuf)
         lsuf0=lsuf
         nctot=0
         nclast=0
         f_hz=fcw(jfcw)
         cpsec0=etime(time)
         if(iirx .le. 0) then
            call mode_find(iiwrt)
            if(jjfail.gt.0) return
         else
c            call rx_modes
            write(6,*)'No real-axis version'
         endif
         cpsec=etime(time) - cpsec0
         if(iiwrt .eq. 1) then
            write(6,220) 'Time taken to find modes = ',cpsec
            write(lusvp,220) 'Time taken to find modes = ',cpsec
220         format(a,f7.2)
         endif
c
         if(iitl .ne. 0 .and. iiwrt .eq. 1) then
            do jzs=1,nzs
               call suffix(hdf_suf,lsuf0,jzs,nzs,SUFX//'zs',3,hdf_suf,
     .            lsuf)
ccjj           jj=(jzs-1)*nrec*nsrc + 1
               jj=1
               call mode_field(phi,dpsi,exp_gbs,plc(jj),tl(jj),jzs)
               lsuf2=lsuf + 3
               hdf_suf(lsuf+1:lsuf2)=SUFX//'tl'
               call out_writex(outroot,lout,hdf_suf,lsuf2,tl(jj),
     .            rec_lab,t_src,nrec,nsrc,dlab,rlab,z4,z4,z4,z4,2,
     .            pllab,'m','km',dblab,'f7.2','f9.3','f7.3',ncall)
               if(iitl .lt. 0) then
                  do j=1,nrec*nsrc
                     pl0=plc(jj+j-1)
                     r4mat1(j)=atan2(aimag(pl0),real(pl0))/pierad
                  enddo
                  hdf_suf(lsuf+1:lsuf2)=SUFX//'ph'
                  call out_writex(outroot,lout,hdf_suf,lsuf2,r4mat1,
     .               rec_lab,t_src,nrec,nsrc,dlab,rlab,z4,z4,z4,z4,2,
     .               'phase','m','km','deg','f7.2','f9.3','f7.3',
     .               ncall)
               endif
            enddo
         endif

         if(iifft .ne. 0) then
            do jzs=1,nzs
               jj=1
               call mode_cmplx_field(tf(jj),phi,dpsi,exp_gbs,jzs,
     .                               jfcw,nfcw)
            enddo
         endif

         if((iimf .eq. 1 .or. iimf .eq. 3) .and. iiwrt .eq. 1) then
            lsuf=lsuf0 + 4
            hdf_suf(lsuf0+1:lsuf)=SUFX//'phi'
            call mfun_fill(phi,r4mat1,r4mat2,hdf_suf,lsuf)
            if(iidiag .eq. -2) call mode_ortho(phi,dphi,psi,dpsi)
c: Do this if you want to output d(phi)/dz also:
cxx         lsuf=lsuf0 + 5
cxx         hdf_suf(lsuf0+1:lsuf)=SUFX//'dphi'
cxx         call mfun_fill(dphi,r4mat1,r4mat2,hdf_suf,lsuf)
         endif
         if((iimf .eq. 2 .or. iimf .eq. 3) .and. iiwrt .eq. 1) then
            lsuf=lsuf0 + 4
            hdf_suf(lsuf0+1:lsuf)=SUFX//'psi'
            call mfun_fill(psi,r4mat1,r4mat2,hdf_suf,lsuf)
            if(iidiag .eq. -2) call mode_ortho(psi,dpsi,psi,dpsi)
         endif
c
c
         if(iiwrt .eq. 1) then
            write(7,202) f_hz,nmode,nctot,kw0
            write(lusvp,203) f_hz,nmode,nctot
202         format('f = ',f8.2,'; # of Modes =',i5,
     .         '; # R1*R2 Calcs =',i6,'; Kw = ',e14.8/
     .         'Mode#        Re(k)/Kw    Im(k)-dB/km     Phase Vel   ',
     .         'Group Vel  Duct#  #Calc')
203         format('f = ',f8.2,'; # of Modes =',i5,
     .         '; # R1*R2 Calcs =',i6)
            do jm=1,nmode
               vph=w/dreal(kn(jm))
               vg=dreal(eig_char(4,jm))
               db_km=8685.889638*dimag(kn(jm))
               if(abs(vph) .lt. 100000.d0) then
                  write(7,201) jm,dreal(kn(jm))/kw0,
     .               db_km,vph,vg,nzref(jm),ncalc(jm)
               else
                  write(7,204) jm,dreal(kn(jm))/kw0,
     .               db_km,vph,vg,nzref(jm),ncalc(jm)
               endif
201   format(i5,2x,e14.8,4x,g11.5,1x,f13.6,1x,f11.5,2x,i5,2x,i5)
204   format(i5,2x,e14.8,4x,g11.5,1x,e13.7,1x,f11.5,2x,i5,2x,i5)
            enddo
         endif
      end do
c
      if(iiwrite .eq. 1) then
c: Check if leaky modes computed in upper or lower isospeed halfspace:
         if(dreal(kn(nmode)) .lt. dreal(xkbp(1,1)) .and.
     .      zsr(nzsr) .gt. zdep(nlay-1) .and. isp(nlay) .eq. 1) 
     .      call print_warn('lower')
         if(dreal(kn(nmode)) .lt. dreal(xkbp(2,1)) .and.
     .      zsr(1) .lt. zdep(1)) call print_warn('upper')
      endif
c
      nfbb=nfcw
c pln 020500      call bb_fft_out(cfmin)
c
      return
      end
ccc
      subroutine print_warn(ch_hsp)
c
      implicit none
      include 'Parms_com'
      character*5 ch_hsp
c
      print *,'WARNING: Field or mode functions requested in ',
     .   'lower halfspace and leaky mode(s) found.'
      write(lusvp,*) 'WARNING: Field or mode functions requested in ',
     .   ch_hsp,' halfspace and leaky mode(s) found.'
      print *,'Be aware that leaky modes are not valid in ',
     .   'isospeed halfspaces!!'
      write(lusvp,*) 'Be aware that leaky modes are not valid in ',
     .   'isospeed halfspaces!!'
      print *,'One solution is to insert a false layer.'
      write(lusvp,*) 'One solution is to insert a false layer.'
c
      return
      end
      subroutine zref_chng
c
c: Changes reference depth to kduct'th duct.
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 j,nsvold,isvold,jx,jx1
c
c: Redo all computations involving reference depth zsvp(nsvmin):
c: From svp_check: Reset ref depth for reflection coefficient calcs:
      nsvold=nsvmin
      isvold=isvmin
      nsvmin=jduct(1,kduct)
      isvmin=jduct(2,kduct)
      rho_duct=geo(isvmin,3,nsvmin)
      jflu(1,1)=nsvmin
      jflu(2,1)=nsvmin
      jflu(1,5)=isvmin-1
      jflu(2,5)=isvmin-2
c: EKW FIX 12/13/95: Account for fact that a fake layer with different 
c: density may be there:
      if(jflu(1,2) .ne. jduct(4,kduct)) then
         jflu(1,2)=jduct(4,kduct)
         rholay(1)=geo(2,3,jflu(1,2))/geo(1,3,jsol(1,1))
      endif
      if(jflu(2,2) .ne. jduct(5,kduct)) then
         jflu(2,2)=jduct(5,kduct)
         rholay(2)=geo(1,3,jflu(2,2))/geo(2,3,jsol(2,1))
      endif
c: Check if this is a fake layer used for seismic modes:
c: Not necessary to go slower in seismic layers? 
cc    if(jduct(3,kduct) .gt. 0) then
cc       phfac0=max(8.e0,phfac)
cc    else
cc       phfac0=phfac
cc    endif
c
c: Invert density ratios at interfaces between old and new reference depths:
      do j=min(nsvold,nsvmin),max(nsvold,nsvmin)-1
         rhorat(j)=1.d0/rhorat(j)
      enddo
c: Change relative depths for mode function depths:
c: FIX: 8-27-97.  For case where nsvold=nsvmin, but isvold .ne. isvmin:
cc    if(nsvold .lt. nsvmin) then
      if(nsvold .le. nsvmin) then
         do j=nsvold+isvold-1,nsvmin+isvmin-2
            jx1=jzmx(j)-1
            do jx=1,nzmx(j)
               zmx(jx1+jx)=h(j) - zmx(jx1+jx)
            enddo
         enddo
      else
         do j=nsvmin+isvmin-1,nsvold+isvold-2
            jx1=jzmx(j)-1
            do jx=1,nzmx(j)
               zmx(jx1+jx)=h(j) - zmx(jx1+jx)
            enddo
         enddo
      endif
      jlmin=min(nsvmin,jsr2j(1))
      jlmax=max(nsvmin,jsr2j(nzsr))
      do j=1,nlay
         iiww(j)=0
         if(j .ge. jlmin .and. j .le. jlmax) iiww(j)=1
c: Temp:
         iiww(j)=1
      enddo
c
      xkref=xk(isvmin,nsvmin)
      xkrat_ref(1)=xkref/kw0
      kw=dreal(xkref)
      cref=w/kw0
c
c: kduct0 holds current duct:
      kduct0=kduct
c
      return
      end
ccc
      subroutine duct_dupl(nm_ok,ndup,r1r2,phiz,dphiz,iduct)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer nm_ok,ndup,iduct,jm0,jm1,nm_temp,jm2,inc,
     .        maxmode
      complex*8 phiz(nzsr,nm_ok),dphiz(nzsr,nm_ok)
      complex*16 kn_found,knd,r1r2(3,4),rr0(3,4),dW_dk
      real*8 magsq
      real*4 pct_diff
c
      if(iiccw .gt. 0) then
         jm1=1
         jm2=nm_ok
         inc=1
      else
         jm1=nm_ok
         jm2=1
         inc=-1
      endif
      kn_found=kn(nmode)
      do jm0=jm1,jm2,inc
         knd=kn(jm0)-kn_found
c: Check if eigenvalues extremely close:
         if(magsq(knd) .le. errdk2) then
c: Check if duct minima are nearly identical. If so, possibility that modes
c: might still be distinct:
c
c: Find current mode more precisely:
            call eig_final(kn(nmode),r1r2,.01d0*errdkms,0,kw0,
     .         iish,iifail,jjfail)
            if(iifail .eq. 1) return
            call eig_enter(r1r2)
            dW_dk=-2.d0*r1r2(1,1)*gamiref*r1r2(2,3)
            call mode_fun(kn(nmode),r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
     .         psi,dpsi,exp_gbs,nmode)
            if(jjfail.eq.1) return
c
c: Find duplicate mode more precisely:
            kduct=nzref(jm0)
c: Change ducts:
            call zref_chng
            nm_temp=nmode
            nmode=jm0
            rr0(1,4)=eig_char(1,nmode)
            rr0(2,4)=eig_char(2,nmode)
            call eig_final(kn(nmode),rr0,.01d0*errdkms,0,kw0,
     .         iish,iifail,jjfail)
            if(iifail .eq. 1) return
            call eig_enter(rr0)
            dW_dk=-2.d0*rr0(1,1)*gamiref*rr0(2,3)
            call mode_fun(kn(nmode),rr0(1,1),rr0(1,2),dW_dk,phi,dphi,
     .         psi,dpsi,exp_gbs,nmode)
            if(jjfail.eq.1) return
c: Change reference depth back:
            nmode=nm_temp
            kduct=nzref(nmode)
            call zref_chng
c
cpln            write(6,*)
cpln            write(6,*)'nmode,ndup,nm_ok,jm0: ',nmode,ndup,nm_ok,jm0
cpln            write(6,*)'jm1,jm2,inc,iiccw: ',jm1,jm2,inc,iiccw
            maxmode=max(jm0,nmode)
            call mfun_comp(nzsr,nmode,jm0,maxmode,phiz,dphiz,nduct,
     .         mzduct,pct_diff)
            if(iidiag .ge. 2) print *,'Dup eigenvalues in different ',
     .         'ducts: pct_diff = ',pct_diff
c
            if(pct_diff .lt. 1.e-2) then
               if(iidiag .ge. 2) then
                  print *,'Duplicate mode found for duct: ',iduct,
     .               kn(jm0),kn(nmode)
               endif
               ndup=ndup + 1
               return
            endif
            if(iidiag .ge. 2) print *,'Close eigenvalues found: ',
     .         kn(jm0),kn(nmode)
         endif
c: EKW FIX. Check imaginary parts also to accommodate evanescent contours
c: that oscillate sideways:
         if(iiccw .gt. 0) then
            if(dreal(knd).lt.0.d0 .and. dimag(knd).gt.0.d0) goto 99
         else
            if(dreal(knd).gt.0.d0 .and. dimag(knd).lt.0.d0) goto 99
         endif
      enddo
c     
99    continue
c: Place new mode into nm_put spot:
      nm_put=nm_put + 1
      kn(nm_put)=kn(nmode)
      if(nm_put .gt. nmode) then
         if(iiwrite .gt. 0)
     .     write(6,*)'nm_put > nmode: ',nm_put,nmode
      end if
c
      call eig_insert(nmode,nm_put,eig_char,nzref,ncalc,mode_phz,
     .   phi,dphi,psi,dpsi,exp_gbs,iishn,max(nm_put,nmode),nzsr)
c
c: New: reset ndup so that two dups in a row must be found to stop looking:
      ndup=0
c
      return
      end
ccc
      subroutine eig_insert(jfrom,jto,eig_char,nzref,ncalc,mode_phz,
     .   phi,dphi,psi,dpsi,exp_gbs,iishn,nmode,nzsr)
c
c: Inserts mode into correct place and moves other up:
      implicit none
      integer*4 jfrom,jto,nmode,nzref(nmode),ncalc(nmode),
     .   iishn(nmode),
     .   nzsr,ji,jzsr
      real*4 mode_phz(3,nmode)
      complex*16 eig_char(5,nmode)
      complex*8 phi(nzsr,nmode),dphi(nzsr,nmode),psi(nzsr,nmode),
     .   dpsi(nzsr,nmode)
      real*8 exp_gbs(nzsr,nmode)
c
      ncalc(jto)=ncalc(jfrom)
      mode_phz(1,jto)=mode_phz(1,jfrom)
      mode_phz(2,jto)=mode_phz(2,jfrom)
      mode_phz(3,jto)=mode_phz(3,jfrom)
      nzref(jto)=nzref(jfrom)
      iishn(jto)=iishn(jfrom)
      do ji=1,5
         eig_char(ji,jto)=eig_char(ji,jfrom)
      enddo
      do jzsr=1,nzsr
         phi(jzsr,jto)=phi(jzsr,jfrom)
         dphi(jzsr,jto)=dphi(jzsr,jfrom)
         psi(jzsr,jto)=psi(jzsr,jfrom)
         dpsi(jzsr,jto)=dpsi(jzsr,jfrom)
         exp_gbs(jzsr,jto)=exp_gbs(jzsr,jfrom)
      enddo
c
      return
      end
ccc
      subroutine eig_sort(nmode,nzsr,kn_indx,eig_char,nzref,
     .   ncalc,iishn,mode_phz,phi,dphi,psi,dpsi,exp_gbs)
c
c: Inserts mode into correct place and moves other up:
      implicit none
      integer*4 nmode,kn_indx(nmode),nzsr,nzref(nmode),ncalc(nmode),
     .   iishn(nmode),jsave,jto,jfrom,iic,j,nm1
      real*4 mode_phz(3,nmode)
      complex*16 eig_char(5,nmode)
      complex*8 phi(nzsr,nmode),dphi(nzsr,nmode),psi(nzsr,nmode),
     .   dpsi(nzsr,nmode)
      real*8 exp_gbs(nzsr,nmode)
c
c: Save first mode:
      jfrom=1
      jsave=1
      nm1=nmode + 1
      call eig_insert(jsave,nmode+1,eig_char,nzref,ncalc,
     .   mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
      do j=1,nmode
         jto=jfrom
         jfrom=kn_indx(jto)
         if(jto .eq. jfrom) then
            iic=1
         elseif(jsave .eq. jfrom) then
            call eig_insert(nmode+1,jto,eig_char,nzref,ncalc,
     .         mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
            iic=1
         else
            call eig_insert(jfrom,jto,eig_char,nzref,ncalc,
     .         mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
            iic=0
         endif
         kn_indx(jto)=0
         if(iic .eq. 1) then
            jsave=jsave+1
            do while(jsave .le. nmode .and. kn_indx(jsave) .eq. 0)
               jsave=jsave+1
            enddo
            call eig_insert(jsave,nmode+1,eig_char,nzref,ncalc,
     .         mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
            jfrom=jsave
         endif
      enddo
c
      return
      end
ccc
      subroutine mfun_comp(nzsr,nmode,jm0,maxmode,phiz,dphiz,nduct,
     .   mzduct,pct_diff)
c
c: Compares mode functions at the duct depths to see if eigenvalue is same.
c
      implicit none
      integer*4 nzsr,nmode,jm0,nduct,mzduct(nduct),jd,jzsr,
     .          maxmode
      real*4 psum,dpsum,pct_diff,psum1,psum2,dpsum1,dpsum2
      complex*8 phiz(nzsr,maxmode),dphiz(nzsr,maxmode)
c
      psum1=0.0
      psum2=0.0
      dpsum1=0.0
      dpsum2=0.0
      psum=0.0
      dpsum=0.0
      do jd=1,nduct
         jzsr=mzduct(jd)
c: Take into account possibility that mode function differ by -1:
         call sum_sq(phiz(jzsr,jm0),phiz(jzsr,nmode),1.0,psum1)
         call sum_sq(phiz(jzsr,jm0),phiz(jzsr,nmode),-1.0,psum2)
         call sum_sq(phiz(jzsr,jm0),phiz(jzsr,jm0),0.0,psum)
         call sum_sq(dphiz(jzsr,jm0),dphiz(jzsr,nmode),1.0,dpsum1)
         call sum_sq(dphiz(jzsr,jm0),dphiz(jzsr,nmode),-1.0,dpsum2)
         call sum_sq(dphiz(jzsr,jm0),dphiz(jzsr,jm0),0.0,dpsum)
      enddo
cxx   pct_diff=min(psum1,psum2)/max(1.e-200,psum)
cxx   pct_diff=min(psum1,psum2)/max(1.e-200,psum)
      pct_diff=max(min(psum1,psum2)/max(1.e-37,psum),
     .   min(dpsum1,dpsum2)/max(1.e-35,dpsum))
c
      return
      end
ccc
      subroutine sum_sq(z1,z2,sg,psum)
c
      implicit none
      complex*8 z1,z2,diff
      real*4 sg,psum
c
      diff=z1 - sg*z2
      psum=psum + real(diff)*real(diff) + aimag(diff)*aimag(diff)
c
      return
      end
ccc
      subroutine zduct_chng(phiz,dphiz,jm,jm2,iichng,k)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jm,jm2,iichng,jzsr,nsv,isv,kdmax,jd
      complex*8 phiz(nzsr,jm2),dphiz(nzsr,jm2),dphi_fac
      real*4 magsq_c8,phimag(25),fac
      complex*16 k
c
      kdmax=1
      do jd=1,nduct
         jzsr=mzduct(jd)
         nsv=jduct(1,jd)
         isv=jduct(2,jd)
         dphi_fac=dphiz(jzsr,jm2)/gami(1,isv,nsv)
         phimag(jd)=magsq_c8(phiz(jzsr,jm2)) + magsq_c8(dphi_fac)
         if(phimag(jd) .gt. phimag(kdmax)) kdmax=jd
      enddo
c
      iichng=0
      isv=jduct(2,kdmax)
      nsv=jduct(1,kdmax)
      fac=2.
      if(nsv .eq. nlay) fac=5.
      if(phimag(kdmax)/phimag(kduct) .gt. fac .and.
     .   dreal(k) .lt. dreal(xk(isv,nsv))) then
c: Change reference duct larger magnitude in other duct exists:
         kduct=kdmax
         iichng=1
      elseif(nzref(jm) .ne. kduct) then
         if(phimag(nzref(jm))/phimag(kduct) .gt. .5) then
c: Change reference duct back to standard one if magnitudes are not
c: clearly different:
            kduct=nzref(jm)
            iichng=1
         endif
      endif
c
      return
      end
      subroutine eig_findm(k0,rr0,r1r2,nm1)
c
c: Finds eigenvalues by moving along constant ln|R1*R2|=1 in complex 
c: k-plane.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iibad,nm1
      complex*16 k0,rr0(3,4),r1r2(3,4),dW_dk
      real*8 xi_real,xi_rat
c
      nmode=nmode + 1
c
      iibad=0
      call mem_lim(nmode,NM_MAX,MLINE,LML,'nmode',5,'NM_MAX',6,iibad,0)
      call mem_lim(nmode,nm_max2,MLINE,LML,'nmode*nzsr',10,
     .   'NSR_NM_MAX',10,iibad,0)
      if(iibad .ne. 0) then
         iifail=1
         nmode=nmode - 1
         return
      elseif(nmode .gt. nm_lim) then
         iidone=1
         return
      endif
      ncalc(nmode)=0
c
      call eig_find0(k0,rr0,kn(nmode),r1r2,nm1)
      k0=kn(nmode)
      rr0(1,4)=r1r2(1,4)
      rr0(2,4)=r1r2(2,4)
c
      if(iidone .eq. 1 .or. iifail .eq. 1) return
c
c: Save eigenvalue characteristics that may be needed later:
      call eig_enter(r1r2)
c
c: Compute mode function at desired depths:
c: Compute analytic derivative of Wronskian with respect to k and w
c: (see p. 105):
c: dW/dk = 2*R1*i*gamref*[-d(R1*R2)/dk]:
c: dW/dw = 2*R1*i*gamref*[-d(R1*R2)/dw]:
      dW_dk=-2.d0*r1r2(1,1)*gamiref*r1r2(2,3)
cpln: Sometimes phix becomes and dW_dk is huge
      call mode_fun(k0,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,psi,dpsi,
     .   exp_gbs,nmode)
      if(jjfail.eq.1) return
      if(iidiag .eq. -2) then
         call phz_calc
      endif
c
c pln
      if(abs(dimag(xi_hsp(1,1))).gt.0.d0) then
         xi_real=dreal(xi_hsp(1,1))
         xi_rat=dabs(xi_real/dimag(xi_hsp(1,1)))
      else
         xi_real=0.d0
         xi_rat=0.d0
      end if
      if(xi_rat .gt. 20. .and. xi_real .lt. 1.) nblm=nblm + 1
c
c: Compute mode amplitudes at duct reference depths:
cxx   if(nblm .gt. nblm_max .and. isp(nlay) .eq. 0) then
      if(isp(nlay) .eq. 0) then
cc       if(iicw .le. 1) then
            call phimag_check(phi,dphi,phi_mag_max,k0,nmode)
cc       else
cc          call phimag_check(phi,dphi,phim_bb(jfbb),k0,nmode)
cc       endif
      endif
c
      return
      end
ccc
      subroutine eig_find0(k0,rr0,k,r1r2,nm1)
c
c: Finds eigenvalues by moving along constant ln|R1*R2|=1 in complex 
c: k-plane.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 nm1,iilhp,nlhp,nfix,nc_ref,iicirc,nmp,ii999
      complex*16 k0,rr0(3,4),k,r1r2(3,4),klast,lnlast,dlnlast,
     .   branch,dk_rad,kstart
      real*8 mag_des0,emagmax0,ephmax0,dph_left,magsq,magsq_knd,
     .   del_r1r2
c
c: iimst=1 means starting at a found mode or at xkref: go 2*pi.
c: iimst=0 means starting anywhere (in a duct, e.g.): always go forward.
c: iimst=-1 means starting near branch cut: go forward if phase > pi/2,
c:          backward if phase < pi/2.
c
      iicirc=0
      nlhp=0
      nc_ref=0
      nfix=0
c: EKW FIX??
cc    phcut=dimag(lncut)
c: Parameters for following mode contour:
      mag_des0=0.d0
      emagmax0=0.1d0
      ephmax0=0.19635d0
999   continue
      iilhp=0
      klast=k0
      lnlast=rr0(1,4)
      dlnlast=rr0(2,4)
      call dph_calc(lnlast,iiccw,iimst,phcut,dph_left,pie,twpie,iidiag)
      iimst=1
c
c: Move along constant magnitude contour:
      call traj_mag(mag_des0,dph_left,klast,lnlast,dlnlast,k,r1r2,
     .   emagmax0,ephmax0,iilhp)
      if(iidone .eq. 1 .or. iifail .eq. 1) return
c
      if(iilhp .eq. 1) then
         call fix_path(k0,rr0,1)
         if(iidone .eq. 1 .or. iifail .eq. 1) return
         nlhp=nlhp + 1
         if(nlhp .gt. 1) then
c: EKW FIX 6/11/97:
            iifail=1
            return
         endif
         goto 999
      endif
c
      call eig_final(k,r1r2,errdkms,iimt,kw0,iish,iifail,jjfail)
      if(iifail .eq. 1) then
         print *,'iifail=1 in eig_final: nmode = ',nmode
         call stop_run
         return
      endif
c
c: Check if we've circled back on xkref:
      magsq_knd=magsq(k - xkref)
      if(magsq_knd .lt. errdk2) then
         if(nsvmin .eq. nlay .and. dreal(xi_hsp(1,1)) .lt. -.2d0) then
         else
            if(iiwrite .gt. 0)
     .       print *,'FALSE MODE AT XKREF IGNORED'
            kcut=xkref + dcmplx(0.d0,.1d0/8686.d0)
            call fix_path(k0,rr0,1)
            phcut=dimag(rr0(1,4))
            if(iidone .eq. 1 .or. iifail .eq. 1) return
            nc_ref=nc_ref + 1
            if(nc_ref .gt. 1) then
               if(iiwrite .gt. 0)
     .          print *,'XKREF island mode was not resolved.'
               iifail=1
               return
            endif
            goto 999
         endif
      endif
c
c: Check if mode trajectory has inadvertently gone in a circle 
c: around an almost-pinched-off island mode:
      magsq_knd=magsq(k - kn(nmode-1))
      if(nmode .gt. nm1 .and. magsq_knd .lt. errdk2) then
         del_r1r2=magsq(r1r2(2,3))*magsq_knd
         if(iidiag .ge. 1) print *,'iicirc incr: nmode = ',nmode,
     .      iicirc,dreal(k)/kw0,dimag(k)*8685.9d0,del_r1r2
         if(del_r1r2 .lt. .1d0) then
            iicirc=iicirc + 1
            if(iicirc .eq. 1) then
c: Set up k0 and rr0 from last mode on true contour found:
c: Check for case where current mode is actually an island mode:
               nmp=nmode-2
               k0=kn(nmp)
               rr0(1,4)=eig_char(1,nmp)
               rr0(2,4)=eig_char(2,nmp)
               emagmax0=emagmax0/20.d0
               ephmax0=ephmax0/20.d0
               if(iimt .ne. 0) call mode_traj(k0,rr0,0)
               goto 999
            elseif(iicirc .eq. 2) then
c: Set up k0 and rr0 from island mode contour. Make sure it is indeed an 
c: island mode:
               nmp=nmode-1
               k0=kn(nmp)
               rr0(1,4)=eig_char(1,nmp)
               rr0(2,4)=eig_char(2,nmp)
               if(iimt .ne. 0) call mode_traj(k0,rr0,0)
               goto 999
            elseif(iicirc .eq. 3) then
               call sdp_find(rr0,0,dk_rad)
25             if(iimt .eq. 1) call mode_traj(k_sdp,rr0,0)
               kstart=k_sdp
               call traj_sdp(k_sdp,ln_sdp,dln_sdp,k0,rr0,kstart,dk_rad)
               if(iifail .ne. 0) then
                  call sdp_find(rr0,1,dk_rad)
                  if(deepest .ne. 0.) then
                     iifail=0
                     goto 25
                  else
                     if(iiwrite.gt.0) then
                        print *,'traj_sdp unsuccessful '
                        print *,'iicirc>4: Tell EKW ',nmode
                     end if
                     nmode=nmode - 1
                     return
                  endif
               endif
               iimst=0
               goto 999
            else
               if(iiwrite.gt.0)
     .          print *,'iicirc>4: Tell EKW ',nmode
               nmode=nmode - 1
               iifail=1
               return
            endif
         endif
      endif
c
      if(iicirc .gt. 0) then
c: Island mode successfully bypassed, so reset emagmax0,ephmax0:
         if(iidiag .ge. 1) print *,'Island mode successfully '//
     .      'bypassed: ',nmode,iicirc,dreal(k)/kw0,dimag(k)*8685.9d0
         iicirc=0
         emagmax0=0.1d0
         ephmax0=0.19635d0
      endif
c  
c: Flag for leaky modes (with respect to p-waves):
      if(dreal(k) .lt. khspmin .or. nsvmin .eq. nlay) then
         iilk=1
         nhigh_max=1
      endif
c
c: Once Re(k) to left of w/crmax (so that mode is not evanescent), keep
c: track of minimum Im(k) in order to stop mode search due to rmin:
      kim_min=dmin1(dimag(k),kim_min)
      if(iilk .eq. 1 .or. nrise .ge. 3) then
         kim_max=kim_min + dkim
      else
c: Safety factor when k not to left of w/crmax:
         kim_max=kim_min + 5.*dkim
      endif
c
c: Check for negative sheets:
      call sheet_look(0,branch)
      if(dimag(k) .lt. 0.d0 .and. dimag(k)*dimag(k) .gt. errdkms) then
         if(iimt .ne. 0) call mode_traj(k,r1r2,-2)
         if(branch .ne. (0.d0,0.d0)) then
c: Mode in lower half plane because we are on wrong sheet:
            if(iidiag .ge. 1) then
               print *,'Informative message: MODE FOUND IN LOWER HP ',
     .            'ON BAD SHEET: ',nmode,k,k/kw0
               print *,'   MODE IGNORED.  ATTEMPTING TO PICK UP PATH',
     .            ' ALONG BRANCH CUT ...'
            endif
            call fix_path(k0,rr0,1)
            if(iidone .eq. 1 .or. iifail .eq. 1) return
         else
c: Continue on trajectory, but skip mode found in lower half plane:
            if(iidiag .ge. 1) then
               print *,'Informative message: MODE FOUND IN LOWER HP ',
     .            'ON OK SHEET: ',nmode,k,k/kw0
               print *,'  MODE IGNORED, SEARCH CONTINUING ON PATH ...'
            endif
            k0=k
            rr0(1,4)=r1r2(1,4)
            rr0(2,4)=r1r2(2,4)
            nlhp=nlhp + 1
            if(nlhp .gt. 2) then
               iifail=1
               return
            endif
         endif
         goto 999
      endif
      if(branch .ne. (0.d0,0.d0)) then
         call bad_sheet(branch,k,r1r2,nfix,k0,rr0,ii999)
         if(iifail .eq. 1 .or. iidone .eq. 1) return
         if(ii999 .eq. 1) goto 999
      endif
c
c: Check if Im(k) so high that mode is weak at shortest range of interest:
cc    if(dimag(k) .gt. kim_max) then
c: EKW FIX 6-2-98:
      if(dimag(k) .gt. kim_max .and. iiccw .gt. 0) then
         nhigh=nhigh + 1
c: Check if nhigh_max modes in a row have been high & traj heading 
c: higher (see p.150):
         if(nhigh .ge. nhigh_max) then
            iidone=1
            return
         endif
      else
         nhigh=0
      endif
      if(dimag(k) .gt. dimag(kn(nmode-1))) then
         nrise=nrise + 1
      else
         nrise=0
      endif
c
      if(iimt .ne. 0) call mode_traj(k,r1r2,1)
c
cxx   if(iidiag .eq. 1) then
c: Compute d ln(R1*R2) / dk numerically as a check:
cxx      kder=k + 1.d-8*kw
cxx      call r1r2_calc(kder,r1r2der,1,0)
cxx      print *,'dL/dk ana: ',nmode,r1r2(2,4)
cxx      print *,'dL/dk num: ',nmode,(r1r2der(1,4)-r1r2(1,4))/(kder-k)
cxx   endif
c
c: Set up k0 and rr0 for next eig_find call:
      k0=k
      rr0(1,4)=r1r2(1,4)
      rr0(2,4)=r1r2(2,4)
c
      return
      end
ccc
      subroutine eig_enter(r1r2)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      complex*16 r1r2(3,4)
c
      eig_char(1,nmode)=r1r2(1,4)
      eig_char(2,nmode)=r1r2(2,4)
      eig_char(3,nmode)=r1r2(3,4)
c: Compute group velocity of mode dw/dk=-d ln(R1*R2)/dk / d ln(R1*R2)/dw
c: (see p. 106):
      eig_char(4,nmode)=-r1r2(2,4)/r1r2(3,4)
      eig_char(5,nmode)=r1r2(1,1)
      nzref(nmode)=kduct
      ncalc(nmode)=ncalc(nmode) + nctot - nclast
      nclast=nctot
c: Keep a record of the sheets as well:
      call iish_code(iish,iish_ref,iishn(nmode),1)
c
      return
      end
ccc
      subroutine eig_final(k,r1r2,errdkms,iimt,kw0,iish,iifail,jjfail)
c
c: Finds the mode eigenvalue once we are in neighborhood of it.
c
      implicit none
      integer*4 iimt,iish(2,2),ntry,ntry2,iifail,jjfail
      complex*16 k,r1r2(3,4),klast,delk,lnr1r2
      real*8 Lmagsq,errdkms,Lmagsq0,magsq,kw0,Lmag_ok
c
      lnr1r2=r1r2(1,4)
      Lmagsq=magsq(lnr1r2)
      ntry=0
30    klast=k
      Lmagsq0=Lmagsq
      delk=-lnr1r2/r1r2(2,4)
      k=klast + delk
      if(k .eq. klast) return
      call r1r2_calc(k,r1r2,3,1,jjfail)
      if(jjfail .gt. 0) return
      if(iimt .ne. 0) call mode_traj(k,r1r2,-2)
      lnr1r2=r1r2(1,4)
      Lmagsq=magsq(lnr1r2)
      Lmag_ok=magsq(r1r2(2,4))*errdkms
c: FIX 4-12-93: Check magnitude squared of ln(R1*R2) also:
      if(Lmagsq .le. Lmag_ok .and. Lmagsq .lt. .01d0) then
         return
      endif
c
      ntry2=0
      do while(Lmagsq .ge. Lmagsq0)
         delk=.5*delk
         k=klast + delk
         if(k .eq. klast) return
         call r1r2_calc(k,r1r2,3,1,jjfail)
         if(jjfail .gt. 0) return
         if(iimt .ne. 0) call mode_traj(k,r1r2,-2)
         lnr1r2=r1r2(1,4)
         Lmagsq=magsq(lnr1r2)
         ntry2=ntry2 + 1
         if(ntry2 .gt. 45) then
            print *,'ntry2 trouble in eig_final: ',dreal(k)/kw0,
     .         dimag(k)*8685.9,delk
            if(ntry2 .gt. 50) then
               iifail=1
               return
            endif
         endif
      enddo
      ntry=ntry + 1
      if(ntry .gt. 45) then
         print *,'ntry trouble in eig_final: ',dreal(k)/kw0,
     .      dimag(k)*8685.9,delk
         if(ntry .gt. 50) then
            iifail=1
            return
         endif
      endif
      goto 30
c
      end
ccc
      subroutine dph_calc(lnlast,iiccw,iimst,phcut,dph_left,
     .   pie,twpie,iidiag)
c
      implicit none
      integer*4 iiccw,iimst,iidiag
      complex*16 lnlast
      real*8 phcut,dph_left,pie,twpie,lnph,dph_left_g,dph_diff
c
      lnph=dimag(lnlast)
c: Compute amount of phase to go to get to next mode:
      if(iiccw .gt. 0) then
         if(iimst .eq. 1) then
c: If starting from a found mode, always go about 2*pi:
            dph_left=twpie - lnph
         elseif(iimst .eq. 0) then
c: If starting from anywhere, always go forward (dph_left positive):
            dph_left=pie + sign(pie,lnph) - lnph
         elseif(iimst .eq. -1) then
c: If starting near branch cut, go to right if mode is within a quarter
c: circle:
            dph_left=pie + sign(pie,lnph) - lnph
            dph_left_g=pie + sign(pie,phcut) - phcut
            dph_diff=dph_left - dph_left_g
cc    print *,'dph_left,_g,diff = ',sngl(dph_left),
cc   .   sngl(dph_left_g),sngl(dph_diff)
            if(dph_diff .gt. pie) then
               dph_left=dph_left - twpie
            elseif(dph_diff .lt. -pie) then
               dph_left=dph_left + twpie
            endif
         endif
      else
         if(iimst .eq. 1) then
c: If starting from a found mode, always go about 2*pi:
            dph_left=-twpie - lnph
c: For iiccw=1, always go to right, or else bad modes on negative sheets
c: can be found:
         elseif(iimst .eq. 0) then
c: If starting from anywhere, always go CW (dph_left negative):
            dph_left=-pie + sign(pie,lnph) - lnph
         elseif(iimst .eq. -1) then
c: If starting near branch cut, go to right if mode is within a quarter
c: circle:
            dph_left=-pie + sign(pie,lnph) - lnph
            dph_left_g=-pie + sign(pie,phcut) - phcut
            dph_diff=dph_left - dph_left_g
cc    print *,'dph_left,_g,diff = ',sngl(dph_left),
cc   .   sngl(dph_left_g),sngl(dph_diff)
            if(dph_diff .gt. pie) then
               dph_left=dph_left - twpie
            elseif(dph_diff .lt. -pie) then
               dph_left=dph_left + twpie
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine stop_run
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
      integer*4 j,jj
c
      print *,'Unable to follow |R1R2|=1 contour.  Giving up.'
      print *,'Will continue with modes already found... '
      if(iidiag .ge. 1) then
         print *,'PROFILE TO FOLLOW. nmode = ',nmode
         do j=1,nlay
            print *,j,1,h(j),(geo(1,jj,j),jj=1,5)
            print *,j,2,h(j),(geo(2,jj,j),jj=1,5)
         enddo
      endif
      nmode=nmode - 1
c
      return
      end
ccc
      subroutine bad_sheet(branch,k,r1r2,nfix,k0,rr0,ii999)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      complex*16 branch,k,r1r2(3,4),k0,rr0(3,4)
      integer*4 nfix,ii999
      real*8 delk_re
c
      ii999=0
      delk_re=(dreal(branch)-dreal(k))*sign(1,iiccw)
      if(delk_re .gt. 0.d0) then
         if(nfix .eq. 0) then
            if(iidiag .ge. 2) then
         print *,'Informative message: MODE ON BAD SHEET FOUND: ',
     .      k/kw0,branch/kw0
         print *,'   MODE IGNORED.  ATTEMPTING TO PICK UP PATH',
     .      ' ALONG BRANCH CUT ...'
            endif
            if(iimt .ne. 0) call mode_traj(k,r1r2,-3)
            call fix_path(k0,rr0,1)
            if(iidone .eq. 1 .or. iifail .eq. 1) return
            nfix=nfix + 1
            ii999=1
         else
      print *,'Mode on wrong sheet near branch point kept',nmode,
     .   dreal(k)/kw0,dimag(k)*8685.9
         endif
      else
         if(iidiag .ge. 1) print *,'OK(?) MODE ON -1 SHEET FOUND: ',
     .      k/kw0,branch/kw0,iish
      endif
c
      return
      end
      subroutine em_calc(xkhsq,xk1sq,xk2sq,eta,etasq,gami,iso,zmr,zmi,
     .   nz,philay,dphilay,Aplay,phix,dphix,exp_gbs,xi,ai,aip,bi,bip,
     .   zzexp,inc,ailay,bilay,zetalay,aisoln,ii1,ii2,jlay,jjhsp,kw0)
c
c: Computes the propagator matrix and its derivative for a given 
c: 1/c^2 linear layer: see Notebook92, p. 69-70.
c: xkhsq=k^2=(w*cos(th)/c0)^2; xk1sq=(w/c1)^2; 
c: etasq=(w^2*beta)^(2/3).
c 
      implicit none
ccx   include 'scairy_com'
      integer*4 iso,nz,inc,aisoln,ii1,ii2,jlay,jjhsp,jz,iir1,iir2,jsign
      complex*16 xkhsq,xk1sq,xk2sq,eta,etasq,gami,philay(2),
     .   dphilay(2),phix(nz),dphix(nz),edif,zeta1,zetb1,A1exp,B1exp
      complex*16 xi(nz),ai(nz),aip(nz),bi(nz),
     .   bip(nz),zzexp(nz),det,A1,B1,A1fac,B1fac,xi0
      complex*16 ailay(2,2),bilay(2,2),zetalay(2)
      real*8 zmr(nz),zmi(nz),Aplay(2),exp_gbs(nz),Aprefa,Aprefb,
     .   exlim,exmax,kw0
      logical stradflag
      data exlim/-50.d0/
c
c: Return if no depths required for this layer:
      if(nz .eq. 0) return
c         write(6,*)'No receivers in this layer!! ',nz
c         return
c      end if
c
      if(iso .eq. 1) then
c: For isovelocity layers:
         if(jlay .eq. jjhsp) then
c: For isospeed halfspace, only allow to match at top of layer (set ii2=ii1):
            call a1b1_iso(philay,dphilay,Aplay,inc,ii1,ii1,det,gami,
     .         zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .         jlay,jjhsp)
         else
            call a1b1_iso(philay,dphilay,Aplay,inc,ii1,ii2,det,gami,
     .         zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .         jlay,jjhsp)
         endif
         do jz=1,nz
            exp_gbs(jz)=0.d0
            edif=gami*dcmplx(zmr(jz),zmi(jz))
            A1exp=Aprefa + (edif - zeta1)
            B1exp=Aprefb - (edif - zetb1)
czs: Take out exponential factor from A1fac and B1fac if necessary:
            if(zmi(jz) .ne. 0.d0) then
               exmax=max(dreal(A1exp),dreal(B1exp))
               A1exp=A1exp - exmax
               B1exp=B1exp - exmax
               exp_gbs(jz)=exmax
            endif
            if(dreal(A1exp) .gt. exlim) then
               if(dreal(B1exp) .gt. exlim) then
                  A1fac=A1*cdexp(A1exp)
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=A1fac + B1fac
                  dphix(jz)=inc*gami*(A1fac - B1fac)
               else
                  A1fac=A1*cdexp(A1exp)
                  phix(jz)=A1fac
                  dphix(jz)=inc*gami*A1fac
               endif
            else
               if(dreal(B1exp) .gt. exlim) then
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=B1fac
                  dphix(jz)=-inc*gami*B1fac
               else
                  phix(jz)=cmplx(0.,0.)
                  dphix(jz)=cmplx(0.,0.)
               endif
            endif
cxx         phix(jz)=A1fac + B1fac
cxx         dphix(jz)=inc*gami*(A1fac - B1fac)
         enddo
      else
c: For linear 1/c^2 profile, use Airy function solutions.
c: Compute Airy function arguments at top and bottom of layer:
         xi0=(xkhsq - xk1sq)/etasq
         do jz=1,nz
            xi(jz)=xi0 - eta*dcmplx(zmr(jz),zmi(jz))
         enddo
cxx      call scairy3(nz,xi,xi0,ai,aip,bi,bip,zzexp,det,-aisoln)
         call scairy3(nz,xi,xi0,ai,aip,bi,bip,zzexp,det,aisoln)
c: Compute A1,B1 coeffs of Ai and Aihat for phi in layer (see p. 87,82):
         if(jlay .eq. jjhsp) then
c: For Airy halfspace, only allow to match at top of layer (set ii2=ii1):
            B1=(0.d0,0.d0)
            A1=philay(ii1)/ailay(1,ii1)
            zeta1=zetalay(ii1)
            Aprefa=Aplay(ii1)
            Aprefb=-1.d300
c
c           call a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii1,eta,det,
c    .         ailay,bilay,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,
c    .         Aprefa,Aprefb,jlay,jjhsp)
cc          print *,'-eta = ',-eta,
cc   .         atan2(-dimag(eta),-dreal(eta))*180/acos(-1.)
c           write(66,100) dreal(cdsqrt(xkhsq))/kw0,
c    .         dimag(cdsqrt(xkhsq))*8685.9,xi0,xi(nz)
100         format(6(e12.5,2x))
         else
            call a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii2,eta,det,
     .         ailay,bilay,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,
     .         Aprefa,Aprefb,jlay,jjhsp)
         endif
cxx   print *,'bot: a1,b1 = ',a1*cdexp(zeta1),b1*cdexp(-zetb1),
cxx  .   iir1,iir2
         do jz=1,nz
c: EKW FIX:
            if(zmi(jz) .ne. 0.d0) then
cj3: check to see if xi(source) needs different solution for stability
               stradflag=(jsign(dimag(xi(jz))) .ne. jsign(dimag(xi0)))
               if(stradflag) then
                  A1exp=Aprefa+zeta1
                  B1exp=Aprefb-zetb1
                  call src_recalc(A1,B1,xi(jz),eta,xkhsq,xk1sq,xk2sq,
     .               etasq,A1exp,B1exp,ii1,ii2,exlim,phix(jz),dphix(jz),
     .               ai(jz),aip(jz),bi(jz),bip(jz),zzexp(jz),
     .               exp_gbs(jz),inc,jlay,jjhsp)
cj3: needed values calculated already by src_recalc, so exit and proceed
                  goto 101
               endif
            endif
c: Make ai,aip get factor of exp(-zeta1); bi,bip get exp(zetb1):
            A1exp=Aprefa + (zeta1 - zzexp(jz))
            B1exp=Aprefb - (zetb1 - zzexp(jz))
            if(zmi(jz) .ne. 0.d0) then
               exmax=max(dreal(A1exp),dreal(B1exp))
               A1exp=A1exp - exmax
               B1exp=B1exp - exmax
               exp_gbs(jz)=exmax
            endif
            if(dreal(A1exp) .gt. exlim) then
               if(dreal(B1exp) .gt. exlim) then
                  A1fac=A1*cdexp(A1exp)
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=A1fac*ai(jz) + B1fac*bi(jz)
                  dphix(jz)=-inc*eta*(A1fac*aip(jz) + B1fac*bip(jz))
               else
                  A1fac=A1*cdexp(A1exp)
                  phix(jz)=A1fac*ai(jz)
                  dphix(jz)=-inc*eta*(A1fac*aip(jz))
               endif
            else
               if(dreal(B1exp) .gt. exlim) then
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=B1fac*bi(jz)
                  dphix(jz)=-inc*eta*(B1fac*bip(jz))
               else
                  phix(jz)=cmplx(0.,0.)
                  dphix(jz)=cmplx(0.,0.)
               endif
            endif
cxx         phix(jz)=A1fac*ai(jz) + B1fac*bi(jz)
cxx         dphix(jz)=-inc*eta*(A1fac*aip(jz) + B1fac*bip(jz))
101      enddo
      endif
99    continue
c
      return
      end
ccc
      subroutine a1b1_iso(philay,dphilay,Aplay,inc,ii1,ii2,det,
     .   gami,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .   jlay,jjhsp)
c
      implicit none
      integer*4 inc,ii1,ii2,iir1,iir2,jlay,jjhsp,ii
      complex*16 philay(2),dphilay(2),det,gami,zetalay(2),A1,B1,
     .   zeta1,zetb1,detinv,term1,term2,term1x,term2x,round_check,
     .   A1x,B1x
      real*8 Aplay(2),Aprefa,Aprefb,rat1,rat2
c
      detinv=0.5d0/gami
      term1=gami*philay(ii2)
      term2=inc*dphilay(ii2)
      A1=round_check(term1,term2,detinv,iir1,rat1)
      ii=ii2
      if(iir1 .eq. 1) then
         A1x=A1
         term1x=gami*philay(ii1)
         term2x=inc*dphilay(ii1)
         A1=round_check(term1x,term2x,detinv,iir1,rat2)
         ii=ii1
         if(iir1 .eq. 1) then
cc          print *,'A1 iso unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               A1=A1x
               ii=ii2
            endif
         endif
      endif
      zeta1=zetalay(ii)
      Aprefa=Aplay(ii)
c
      if(jlay .eq. jjhsp) then
c: For halfspace, B1=0:
         B1=dcmplx(0.d0,0.d0)
         zetb1=zetalay(ii)
         Aprefb=-1.d100
         return
      endif
c
      term2=-term2
      B1=round_check(term1,term2,detinv,iir2,rat1)
      ii=ii2
      if(iir2 .eq. 1) then
         B1x=B1
         term1x=gami*philay(ii1)
         term2x=-inc*dphilay(ii1)
         B1=round_check(term1x,term2x,detinv,iir2,rat2)
         ii=ii1
         if(iir2 .eq. 1) then
cc          print *,'B1 iso unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               B1=B1x
               ii=ii2
            endif
         endif
      endif
      zetb1=zetalay(ii)
      Aprefb=Aplay(ii)
c
      return
      end
ccc
      subroutine a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii2,eta,det,
     .   ailay,bilay,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .   jlay,jjhsp)
c
      implicit none
      integer*4 inc,ii1,ii2,iir1,iir2,ii,jlay,jjhsp
      complex*16 philay(2),dphilay(2),eta,det,ailay(2,2),bilay(2,2),
     .   zetalay(2),A1,B1,zeta1,zetb1,dphfac,dphfac2,detinv,
     .   term1,term2,round_check,A1x,B1x
cc    complex*16 A1exp,B1exp,A1fac,B1fac,ph,dph
      real*8 Aplay(2),Aprefa,Aprefb,rat1,rat2
c
      detinv=1.d0/det
      dphfac=inc*dphilay(ii2)/eta
      term1=bilay(2,ii2)*philay(ii2)
      term2=bilay(1,ii2)*dphfac
      A1=round_check(term1,term2,detinv,iir1,rat1)
cxx   A1=detinv*(bilay(2,ii2)*philay(ii2) + bilay(1,ii2)*dphfac)
      ii=ii2
      if(iir1 .eq. 1) then
         A1x=A1
         dphfac2=inc*dphilay(ii1)/eta
         term1=bilay(2,ii1)*philay(ii1)
         term2=bilay(1,ii1)*dphfac2
         A1=round_check(term1,term2,detinv,iir1,rat2)
         ii=ii1
         if(iir1 .eq. 1) then
cc          print *,'A1 unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               A1=A1x
               ii=ii2
            endif
         endif
      endif
      Aprefa=Aplay(ii)
      zeta1=zetalay(ii)
c
      if(jlay .eq. jjhsp) then
         B1=dcmplx(0.d0,0.d0)
         Aprefb=-1.d100
         zetb1=zetalay(ii)
         return
      endif
c
      term1=-ailay(2,ii2)*philay(ii2)
      term2=-ailay(1,ii2)*dphfac
      B1=round_check(term1,term2,detinv,iir2,rat1)
cxx   B1=detinv*(-ailay(2,ii2)*philay(ii2) - ailay(1,ii2)*dphfac)
      ii=ii2
      if(iir2 .eq. 1) then
         B1x=B1
         dphfac2=inc*dphilay(ii1)/eta
         term1=-ailay(2,ii1)*philay(ii1)
         term2=-ailay(1,ii1)*dphfac2
         B1=round_check(term1,term2,detinv,iir2,rat2)
         ii=ii1
         if(iir2 .eq. 1) then
cc          print *,'B1 unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               B1=B1x
               ii=ii2
            endif
         endif
      endif
      Aprefb=Aplay(ii)
      zetb1=zetalay(ii)
c
cc       do ii=1,2
cc          A1exp=Aprefa + (zeta1 - zetalay(ii))
cc          B1exp=Aprefb - (zetb1 - zetalay(ii))
cc          A1fac=A1*cdexp(A1exp)
cc          B1fac=B1*cdexp(B1exp)
cc          ph=A1fac*ailay(1,ii) + B1fac*bilay(1,ii)
cc          dph=-inc*eta*(A1fac*ailay(2,ii) + B1fac*bilay(2,ii))
cc    print *,'a1b1_airy: ii,phi = ',ii,ph,philay(ii)*dexp(Aplay(ii))
cc    print *,'   dphi = ',ii,dph,dphilay(ii)*dexp(Aplay(ii))
cc       enddo
c
      return
      end
ccc
      function round_check(term1,term2,deninv,iir,ratio)
c
c: Checks to see if cancellation occurs in term1+term2.
c: Returns round_check=(term1 + term2)*deninv if not.
      implicit none
      integer*4 iir
      complex*16 round_check,term1,term2,deninv,numerator
      real*8 magsq,magsum,ratio
c
      numerator=term1 + term2
      magsum=magsq(term1) + magsq(term2)
      ratio=magsq(numerator)/magsum
cxx   if(ratio .gt. 1.d-18) then
      if(ratio .gt. 1.d-14) then
         round_check=numerator*deninv
         iir=0
      else
         iir=1
         round_check=(0.,0.)
cxx      print *,'Round to zero: ',term1,term2,numerator,deninv
      endif
c
      return
      end
ccc
      function jsign(number)
c
      implicit none
      real*8 number
      integer*4 jsign
c
      if (number .gt. 0d0) then
         jsign=1
         return
      else
         if(number .lt. 0d0) then
            jsign=-1
            return
         else
            jsign=0
            return
         endif
      endif
c
      return
      end
ccc
      subroutine src_recalc(A1,B1,xisrc,eta,xkhsq,xk1sq,xk2sq,etasq,
     . A1exp,B1exp,ii1,ii2,exlim,phisrc,dphisrc,aisrc,aipsrc,bisrc,
     . bipsrc,zzsrc,srcexp,inc,jlay,jjhsp)
c
c: Subroutine to recalculate source excitations if there a problem with    ::
c: xi(source,receivers) on opposite sides of real axis.  (part of patch j3)::
c
      implicit none
ccx   include 'scairy_com'
c: INPUT                ::
      complex*16 A1,B1,xisrc,eta,xkhsq,xk1sq,xk2sq,etasq,A1exp,
     .   B1exp
      integer*4  ii1,ii2,inc,jlay,jjhsp
      real*8 exlim
c: OUTPUT       ::
      complex*16 phisrc,dphisrc,aisrc,aipsrc,bisrc,bipsrc,zzsrc
      real*8  srcexp
c: OTHER                ::
      complex*16 Gii1,Gii2,ailay(2,2),bilay(2,2),A1f,B1f,philay(2),
     .   dphilay(2),A2,B2,A2exp,B2exp,A2fac,B2fac,xii1,xii2,det,
     .   detii1,detii2,zetalay(2),zeta2,zetb2,A1fac,B1fac
      real*8 eps,Aplay(2),exmax,Aprefa2,Aprefb2
      integer*4 iibad,jsign,iir1,iir2,ii,aisoln
c
c:  Calculate old xis, and zero out imag. part::
      Gii1=(xkhsq-xk1sq)/etasq
      Gii2=(xkhsq-xk2sq)/etasq
      eps=1d-20 
      xii1=cmplx(dreal(Gii1),-eps*jsign(dimag(xisrc)))
      xii2=cmplx(dreal(Gii2),-eps*jsign(dimag(xisrc)))
c:  First calculate phi, dphi on the same side of real axis as orig layer depths::
      call scairy2(xii1,ailay(1,ii1),ailay(2,ii1),bilay(1,ii1),
     .   bilay(2,ii1),zetalay(ii1),detii1,xii2,ailay(1,ii2),
     .   ailay(2,ii2),bilay(1,ii2),bilay(2,ii2),
     .   zetalay(ii2),detii2,aisoln,iibad)
      A1fac=A1*cdexp(A1exp) 
      B1fac=B1*cdexp(B1exp)
      do ii=1,2
         A1f=A1fac*cdexp(-zetalay(ii))
         B1f=B1fac*cdexp(zetalay(ii))
         philay(ii)=(A1f*ailay(1,ii)+B1f*bilay(1,ii))
         dphilay(ii)=-inc*eta*(A1f*ailay(2,ii)+B1f*bilay(2,ii))
         Aplay(ii)=0.d0
      enddo
c:  Now get ailay, etc... on other side of line ::
      xii1=cmplx(dreal(xii1),eps*jsign(dimag(xisrc)))
      xii2=cmplx(dreal(xii2),eps*jsign(dimag(xisrc)))
      call scairy2(xii1,ailay(1,ii1),ailay(2,ii1),bilay(1,ii1),
     .   bilay(2,ii1),zetalay(ii1),detii1,xii2,ailay(1,ii2),
     .   ailay(2,ii2),bilay(1,ii2),bilay(2,ii2),
     .   zetalay(ii2),detii2,aisoln,iibad)
      call scairy3(1,xisrc,xisrc,aisrc,aipsrc,bisrc,bipsrc,zzsrc,det,
     .   aisoln)
      call a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii2,eta,det,ailay,
     .     bilay,zetalay,A2,B2,iir1,iir2,zeta2,zetb2,Aprefa2,
     .     Aprefb2,jlay,jjhsp)
      A2exp=Aprefa2+(zeta2-zzsrc)
      B2exp=Aprefb2-(zetb2-zzsrc)
      exmax=max(dreal(A2exp),dreal(B2exp))
      A2exp=A2exp-exmax
      B2exp=B2exp-exmax
      srcexp=exmax
      if(dreal(A2exp) .gt. exlim) then
         if(dreal(B2exp) .gt. exlim) then
            A2fac=A2*cdexp(A2exp)
            B2fac=B2*cdexp(B2exp)
            phisrc=A2fac*aisrc + B2fac*bisrc
            dphisrc=-inc*eta*(A2fac*aipsrc + B2fac*bipsrc)
         else
            A2fac=A2*cdexp(A2exp)
            phisrc=A2fac*aisrc
            dphisrc=-inc*eta*(A2fac*aipsrc)
         endif
      else
         if(dreal(B2exp) .gt. exlim) then
            B2fac=B2*cdexp(B2exp)
            phisrc=B2fac*bisrc
            dphisrc=-inc*eta*(B2fac*bipsrc)
         else
            phisrc=cmplx(0.,0.)
            dphisrc=cmplx(0.,0.)
        endif
      endif
c
      return
      end
      subroutine ep_calc(xkh,w,xkhsq,xk1sq,xk2sq,eta,etasq,gami,h,
     .   e11,e12,e21,e22,zzexp,iso,ihf,isx,ndv,iiw,ailay,bilay,
     .   zetalay,aisoln,ii1,ii2)
c
c: Computes the propagator matrix and its derivative for a given 
c: 1/c^2 linear layer: see Notebook92, p. 69-70.
c: xkh=horizontal wavenumber; xkqsq=xkh^2; (xk1sq=(w/c1)^2; 
c: xk2sq=(w/c2)^2; eta=(w^2*beta)^(1/3); etasq=eta^2; beta=gradient;
c: gami=i*sqrt(xk1sq - xkhsq); h=layer thickness.
c: E and E'=dE/dk must be multiplied by exp(zzexp) to get true values.
c: For ndv=2, the derivative w.r.t. k is computed.
c: For ndv=3, the derivative w.r.t. w is computed.
c: ihf=1 means h is inversely proportional to frequency.
c 
      implicit none
cxx   include 'scairy_com'
      integer*4 iso,ihf,isx,ndv,j,iiw,iibad,aisoln,ii1,ii2
      complex*16 xkh,xkhsq,xk1sq,xk2sq,eta,etasq,gami(3),
     .   e11(3),e12(3),e21(3),e22(3),zzexp
      complex*16 xi1,xi2,ai1,aip1,bi1,bip1,ai2,aip2,bi2,bip2,
     .   det1,det2,zfac,zzexp1,zzexp2,igamhp
      complex*16 ailay(2,2),bilay(2,2),zetalay(2),gamih
      real*8 h,w,pie
      data pie/3.14159265358979/
c
      if(h .eq. 0.d0) return
      if(iso .eq. 1) then
         gamih=gami(1)*h
c: For isovelocity layers, combine terms that will be needed in
c: rp_slay, rp_slay0, or rp_flay:
         if(isx .eq. 1) then
c: For isx=1, set e11=e11 + gami*e12, e12=e11 - gami*e12,
c: e21=e21 + gami*e22, e22=e21 - gami*e22 (see rp_slay,rp_slay0):
            if(dreal(gami(1)) .gt. 0.d0) then
c: Remove zzexp=exp(gami*h):
               zzexp=gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=(1.d0,0.d0)
               e12(1)=zfac
               e21(1)=gami(1)
               e22(1)=-gami(1)*zfac
               if(ndv .ge. 2) then
                  j=2
                  e11(j)=gami(j)*h
                  e12(j)=-e11(j)*zfac
                  e21(j)=gami(j)*(zzexp + 1.d0)
                  e22(j)=gami(j)*(zzexp - 1.d0)*zfac
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e11(j)=gami(j)*h
                     e12(j)=-e11(j)*zfac
                     e21(j)=gami(j)*(zzexp + 1.d0)
                     e22(j)=gami(j)*(zzexp - 1.d0)*zfac
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=igamhp
                     e12(j)=-igamhp*zfac
                     e21(j)=gami(1)*igamhp + gami(j)
                     e22(j)=(gami(1)*igamhp - gami(j))*zfac
                  endif
               endif
            else
c: Remove zzexp=exp(-gami*h):
               zzexp=-gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=zfac
               e12(1)=(1.d0,0.d0)
               e21(1)=gami(1)*zfac
               e22(1)=-gami(1)
               if(ndv .ge. 2) then
                  j=2
                  e12(j)=-gami(j)*h
                  e11(j)=-e12(j)*zfac
                  e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                  e22(j)=-gami(j)*(zzexp + 1.d0)
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e12(j)=-gami(j)*h
                     e11(j)=-e12(j)*zfac
                     e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                     e22(j)=-gami(j)*(zzexp + 1.d0)
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=igamhp*zfac
                     e12(j)=-igamhp
                     e21(j)=(gami(1)*igamhp + gami(j))*zfac
                     e22(j)=gami(1)*igamhp - gami(j)
                  endif
               endif
            endif
         else
c: For isx=2, set e11=gami*e11 + e21, e21=gami*e11 - e21, 
c: e12=gami*e12 + e22, e22=gami*e12 - e22 (see rp_flay): 
            if(dreal(gami(1)) .gt. 0.d0) then
c: Remove zzexp=exp(gami*h):
               zzexp=gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=gami(1)
               e21(1)=gami(1)*zfac
               e12(1)=(1.d0,0.d0)
               e22(1)=-zfac
               if(ndv .ge. 2) then
                  j=2
                  e11(j)=gami(j)*(zzexp + 1.d0)
                  e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                  e12(j)=gami(j)*h
                  e22(j)=e12(j)*zfac
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e11(j)=gami(j)*(zzexp + 1.d0)
                     e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                     e12(j)=gami(j)*h
                     e22(j)=e12(j)*zfac
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=gami(1)*igamhp + gami(j)
                     e21(j)=(-gami(1)*igamhp + gami(j))*zfac
                     e12(j)=igamhp
                     e22(j)=igamhp*zfac
                  endif
               endif
            else
c: Remove zzexp=exp(-gami*h):
               zzexp=-gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=gami(1)*zfac
               e21(1)=gami(1)
               e12(1)=zfac
               e22(1)=(-1.d0,0.d0)
               if(ndv .ge. 2) then
                  j=2
                  e11(j)=gami(j)*(-zzexp + 1.d0)*zfac
                  e21(j)=gami(j)*(zzexp + 1.d0)
                  e22(j)=gami(j)*h
                  e12(j)=e22(j)*zfac
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e11(j)=gami(j)*(-zzexp + 1.d0)*zfac
                     e21(j)=gami(j)*(zzexp + 1.d0)
                     e22(j)=gami(j)*h
                     e12(j)=e22(j)*zfac
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=(gami(1)*igamhp + gami(j))*zfac
                     e21(j)=-gami(1)*igamhp + gami(j)
                     e12(j)=igamhp*zfac
                     e22(j)=igamhp
                  endif
               endif
            endif
         endif
         if(iiw .eq. 1) then
            zetalay(ii1)=0.d0
            zetalay(ii2)=gamih
         endif
c
c: Compute propagator matrix for Airy layer:
      else
c: For linear 1/c^2 profile, use Airy function solutions.
c: Compute Airy function arguments at top and bottom of layer:
         xi1=(xkhsq - xk1sq)/etasq
         xi2=(xkhsq - xk2sq)/etasq
c: Compute numerically stable Airy function solutions and derivatives:
c: Note: bi,bip are actually Ai(xi*eim23) and (d/dxi)(Ai(xi*eim23)) for
c: Im(xi)>=0, and Ai(xi*ei23) and (d/dxi)(Ai(xi*ei23)) for Im(xi)<0.
         call scairy2(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,ai2,aip2,
     .      bi2,bip2,zzexp2,det2,aisoln,iibad)
c
         if(iiw .eq. 1) then
            ailay(1,ii1)=ai1
            ailay(1,ii2)=ai2
            ailay(2,ii1)=aip1
            ailay(2,ii2)=aip2
            bilay(1,ii1)=bi1
            bilay(1,ii2)=bi2
            bilay(2,ii1)=bip1
            bilay(2,ii2)=bip2
            zetalay(ii1)=zzexp1
            zetalay(ii2)=zzexp2
         endif
c
         if(iibad .eq. 0) then
c: Compute propagator matrix for Airy layer:
            call airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,ai2,
     .         aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,zzexp,
     .         xkh,xkhsq,xk1sq,xk2sq,etasq,w,ndv)
         else
            call ai_strad(xi1,xk1sq,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,
     .         xk2sq,ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,
     .         zzexp,xkh,xkhsq,etasq,w,ndv)
         endif
c
      endif
c
      return
      end
ccc
      subroutine airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,
     .   ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,zzexp,
     .   xkh,xkhsq,xk1sq,xk2sq,etasq,w,ndv)
c
      implicit none
      integer*4 ndv
      complex*16 xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,ai2,aip2,
     .   bi2,bip2,zzexp2,det2,eta,e11(3),e12(3),e21(3),e22(3),
     .   zzexp,xkh,xkhsq,xk1sq,xk2sq,etasq,dexp,zfac,e21_eta,eta_e12,
     .   dxi_dk,dxi_dw1,dxi_dw2,dxi_fac
      real*8 w,deta_eta
c
c: 1st terms of E must be multiplied by exp(dexp); 2nd terms by exp(-dexp):
      dexp=zzexp1 - zzexp2
      if(dreal(dexp) .gt. 0.d0) then
c: Exponential factor of first of two terms dominates:
         zzexp=dexp
         zfac=cdexp(-2.d0*dexp)
c: Remove old factor from and insert new factor to second terms:
         e11(1)=(ai2*bip1 - zfac*bi2*aip1)/det1
         e12(1)=(-ai2*bi1 + zfac*bi2*ai1)/(-det1*eta)
         e21(1)=-eta*(aip2*bip1 - zfac*bip2*aip1)/det1
         e22(1)=(-aip2*bi1 + zfac*bip2*ai1)/det1
      else
c: Exponential factor of second of two terms dominates:
         zzexp=-dexp
         zfac=cdexp(2.d0*dexp)
c: Remove old factor from and insert new factor to first terms:
         e11(1)=(zfac*ai2*bip1 - bi2*aip1)/det1
         e12(1)=(-zfac*ai2*bi1 + bi2*ai1)/(-det1*eta)
         e21(1)=-eta*(zfac*aip2*bip1 - bip2*aip1)/det1
         e22(1)=(-zfac*aip2*bi1 + bip2*ai1)/det1
      endif
c
      if(ndv .ge. 2) then
         e21_eta=e21(1)/eta
         eta_e12=eta*e12(1)
c: Compute derivative w.r.t. k of exponential factor to do removed:
c: (see p.101-102 of Notebook):
         dxi_dk=2.d0*xkh/etasq
         e11(2)=dxi_dk*(xi1*eta_e12 - e21_eta)
         e12(2)=dxi_dk*(e11(1) - e22(1))/eta
         e21(2)=dxi_dk*(xi1*e22(1) - xi2*e11(1))*eta
         e22(2)=dxi_dk*(e21_eta - xi2*eta_e12)
         if(ndv .ge. 3) then
c: Compute derivative w.r.t. w (see p.108 of Notebook):
            deta_eta=1.d0/(1.5d0*w)
cxx         dxi_fac=-2.d0/(3.d0*w*etasq)
            dxi_fac=-deta_eta/etasq
            dxi_dw1=dxi_fac*(xk1sq + 2.d0*xkhsq)
            dxi_dw2=dxi_fac*(xk2sq + 2.d0*xkhsq)
            e11(3)=(xi1*eta_e12*dxi_dw1 - e21_eta*dxi_dw2)
            e12(3)=(e11(1)*dxi_dw1 - e22(1)*dxi_dw2)/eta - 
     .         e12(1)*deta_eta
            e21(3)=(xi1*e22(1)*dxi_dw1 - xi2*e11(1)*dxi_dw2)*eta +
     .         e21(1)*deta_eta
            e22(3)=(e21_eta*dxi_dw1 - xi2*eta_e12*dxi_dw2)
         endif
      endif
c
      return
      end
c
      subroutine ai_strad(xi1,xk1sq,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,
     .   xk2sq,ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,zzexp,
     .   xkh,xkhsq,etasq,w,ndv)
c
c: Computes the propagator matrix for a layer when xi1 and xi2 straddle
c: the real axis so that different pairs of solutions were originally
c: computed in scairy2.  When the determinant is bad in scairy2,
c: ai_strad is called. 
c: The intersection from xi1 to xi2 with the real axis is found,
c: and the layer is split into two layers, each of which now has the
c: same solution in scairy2.  The propagator matrices are then multiplied
c: together.
c
      implicit none
      integer*4 ndv,iibad,j,aisoln
      real*8 w,rfac
      complex*16 xi1,xk1sq,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,xk2sq,
     .   ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11(3),e12(3),e21(3),e22(3),
     .   zzexp,xkh,xkhsq,etasq
      complex*16 xi01,xk01sq,ai01,aip01,bi01,bip01,zzexp01,det01,
     .   xi02,xk02sq,ai02,aip02,bi02,bip02,zzexp02,det02,f11(3),f12(3),
     .   f21(3),f22(3),g11(3),g12(3),g21(3),g22(3),zzexpf,zzexpg
c
c: Compute factor for depth at which xi crosses real axis:
      rfac=-dimag(xi1)/(dimag(xi2) - dimag(xi1))
c: Let xi01 be slightly on xi1's side of real axis at intersection point:
      xi01=dcmplx(dreal(xi1) + (dreal(xi2)-dreal(xi1))*rfac,
     .   dsign(1.d-20,dimag(xi1)))
c: Compute interpolated value for xksq:
      xk01sq=xk1sq + rfac*(xk2sq - xk1sq)
c: Get Airy functions for first layer from xi1 to xi01:
      call scairy2(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi01,ai01,aip01,
     .   bi01,bip01,zzexp01,det01,aisoln,iibad)
      if(iibad .eq. 1) print *,'iibad=1 in ai_strad!!'
c: Compute propagator matrix for first layer:
      call airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi01,
     .   ai01,aip01,bi01,bip01,zzexp01,det01,eta,f11,f12,f21,
     .   f22,zzexpf,xkh,xkhsq,xk1sq,xk01sq,etasq,w,ndv)
c
c: Let xi02 be slightly on xi2's side of real axis at intersection point:
      xi02=dconjg(xi01)
      xk02sq=xk01sq
c: Get Airy functions for first layer from xi02 to xi2:
      call scairy2(xi02,ai02,aip02,bi02,bip02,zzexp02,det02,xi2,
     .   ai2,aip2,bi2,bip2,zzexp2,det2,aisoln,iibad)
      if(iibad .eq. 1) print *,'iibad=1 in ai_strad!!'
c: Compute propagator matrix for first layer:
      call airy_prop(xi02,ai02,aip02,bi02,bip02,zzexp02,det02,
     .   xi2,ai2,aip2,bi2,bip2,zzexp2,det2,eta,g11,g12,g21,g22,
     .   zzexpg,xkh,xkhsq,xk02sq,xk2sq,etasq,w,ndv)
c
c: Multiply propagator matrices together (G*F):
      e11(1)=g11(1)*f11(1) + g12(1)*f21(1)
      e21(1)=g21(1)*f11(1) + g22(1)*f21(1)
      e12(1)=g11(1)*f12(1) + g12(1)*f22(1)
      e22(1)=g21(1)*f12(1) + g22(1)*f22(1)
      zzexp=zzexpf + zzexpg
c
c: Compute derivatives using the chain rule:
      do j=2,ndv
         e11(j)=g11(1)*f11(j) + g11(j)*f11(1) + g12(1)*f21(j) +
     .      g12(j)*f21(1)
         e21(j)=g21(1)*f11(j) + g21(j)*f11(1) + g22(1)*f21(j) +
     .      g22(j)*f21(1)
         e12(j)=g11(1)*f12(j) + g11(j)*f12(1) + g12(1)*f22(j) +
     .      g12(j)*f22(1)
         e22(j)=g21(1)*f12(j) + g21(j)*f12(1) + g22(1)*f22(j) +
     .      g22(j)*f22(1)
      enddo
c
      return
      end
       subroutine bb_brute
c
c: Computes the broadband field from fmin to fmax in steps of df
c: using brute force calls of the usual CW routines
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      include 'sector_env_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm,jfbb,nmbbtot,nm_fmax,nctot0,ncjf,jzjm,j
      real*8 watscale 
      real*4 ncperm
c
c      write(*,*)'hello'
      nm_put=0
      nmode=0
      nctot0=0
      nctot=0
      ncperm=0
      nm_tot=0
      nm_fmax=0
c
      call bb_init
      call zmx_init
      f_hz=faxbb(nfbb)
      call freq_init
      call bb_ftinit
c
      nmbbtot=0
      if(i_call_porter .eq. 1) then
         if(iiwrite. gt. 0)
     .        write(6,*)'Full-blown calculation'
         do jfbb=nfbb,1,-1
            nctot0=nctot
            f_hz=faxbb(jfbb)
            call freq_chng
            call mode_find(0)
            if(jjfail.gt. 0) return
            if(jfbb .eq. nfbb) then
               if(iiAih(1) .eq. -1 .and. iiAih(2) .eq. -1) then
                  nm_fmax=nm_tot
               else
                  nm_fmax=max(1.1*nm_tot,float(nm_tot) + 12)
               endif
               if(iiout .ne. 0) call bb_out_init(nm_fmax)
            endif
            do jm=1,nm_tot
               call bb_field(kn(jm),phi,dphi,dpsi,exp_gbs,
     .                            jm,tf,jfbb,jm)
               call bb_enter(kn(jm),eig_char(1,jm),
     .              eig_char(4,jm),eig_char(5,jm),knbb,eig_bb,
     .              jfbb,jm,nfbb,nm_tot,iish,iish_ref,iish_bb,nmbb)
               if(iift .eq. 1) call bb_write(jm,jfbb,faxbb,kn(jm),
     .                                       kw0,iish)
            enddo
            nmbb(jfbb)=min(nm_fmax,nm_tot)
            if(iiout .ne. 0) then
               if (nm_tot .gt. nm_fmax) print *,'nm_tot>nm_fmax',
     .              f_hz,nm_tot,nm_fmax
               write(33,rec=nh_off+jfbb) (kn(jm),jm=1,nmbb(jfbb)),
     .              (phi(jzjm),jzjm=1,nzsr*nmbb(jfbb))
            endif
            nmbbtot=nmbbtot + nm_tot
            ncjf=nctot - nctot0
            if(iiwrite .gt. 0)
     .           print *,'Done f,nm_tot,#R1R2/mode = ',
     .           sngl(f_hz),nm_tot,float(ncjf)/max(1,nm_tot)
            if(i_geom.eq.0) then
               do jzjm=1,nzsr*nmbb(jfbb)
                  phi_reuse(jfbb,jzjm)=phi(jzjm)
               end do
            end if
         enddo
      else
         if(iiwrite .gt. 0)
     .        write(6,*)'Reuse eigenvalues and mode functions'
c: pln
c: scale mode functions when searching for water depth
         if(r_h0(1).ne.iwat0) then
            watscale=r_h0(1)/iwat0
            do j=1,nzsrgeom
               zsrgeom(j)=zsrgeom(j)*watscale
            end do
            if(iiwrite.gt.0) then
               write(6,*)'New water depth:  ',r_h0(1)
               write(6,*)'Ref. water depth: ',iwat0
               write(6,*)
            end if
         end if
         do jfbb=nfbb,1,-1
            do jzjm=1,nzsrgeom*nmbb(jfbb)
               phi(jzjm)=phi_reuse(jfbb,jzjm)
            end do
            do jm=1,nmbb(jfbb)
               call bb_field_reuse(knbb((jm-1)*nfbb+jfbb),
     .              phi,exp_gbs,jm,tf,jfbb)
            end do
            nmbbtot=nmbbtot + nm_tot
            ncjf=nctot - nctot0
            if(iiwrite .gt. 0)
     .           print *,'Done f,nm_tot,#R1R2/mode = ',
     .           faxbb(jfbb),nm_tot,float(ncjf)/max(1,nm_tot)
         end do
c: pln
c: scale mode functions back to original for water depth
         if(r_h0(1).ne.iwat0) then
            watscale=iwat0/r_h0(1)
            do j=1,nzsrgeom
               zsrgeom(j)=zsrgeom(j)*watscale
            end do
         end if
      end if
c
      call bb_done(nm_fmax)
c
c: Output FFT file:
      
c no time for this Peter gerstoft      
c      call bb_fft_out
c
      ncperm=float(nctot)/float(max(1,nmbbtot))
      write(lusvp,120) nmbbtot,nctot,ncperm
120   format('CW LOOP # MODES = ',i8,'; # R1R2 CALCS = ',i8,
     .   '; #CALCS/MODE = ',f5.2)
c
      return
      end
      subroutine fix_path(k0,rr0,iileft)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iileft,j
      complex*16 k0,rr0(3,4),klast,lnlast,dlnlast,branch,delk,kcut0
      real*8 ph_des0,dmag_left
c
      k0=kcut
      kcut0=kcut
      call sheet_init(k0,1,iish,iish_ref)
c
c: phcut is phase along |R1R2|=1 contour at point where contour crossed
c: branch cut:
      phcut=dimag(lncut)
c
      if(iimt .ne. 0) call mode_traj(k0,rr0,0)
      call r1r2_calc(k0,rr0,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnlast=rr0(1,4)
      if(dabs(dreal(lnlast)) .gt. 0.05d0) then
c: Move along phase contour to |R1R2|=1 contour:
         klast=k0
         dlnlast=rr0(2,4)
         dmag_left=-dreal(lnlast)
         ph_des0=dimag(lnlast)
         call traj_phase(ph_des0,dmag_left,klast,lnlast,
     .      dlnlast,k0,rr0,.05d0,0.19635d0)
         if(jjfail.gt.0) return
         if(iifail .eq. 1) then
            print *,'Failure in traj_phase from fix_path: ',nmode,k0,rr0
            call stop_run
            return
         endif
         if(iidiag .ge. 2) print *,'traj_phase from fix_path: ',
     .      k0,(rr0(j,4),j=1,2)
         call sheet_look(0,branch)
         if(branch .ne. (0.d0,0.d0)) then
            delk=k0 - branch
            if(dreal(delk) .lt. 0.d0 .or. dimag(delk) .lt. 0.d0) then
               k0=kcut0
               call sheet_init(k0,1,iish,iish_ref)
      print *,'Calling cut_cross2 after bad traj_phase: ',k0/kw0
               call r1r2_calc(k0,rr0,2,0,jjfail)
               if(jjfail .gt. 0) return
               call cut_cross2(k0,rr0)
            endif
         endif
      endif
      if(iileft .eq. 1 .and. dimag(rr0(2,4)) .ge. 0.d0) then
c: Found contour, but heading to right in complex k-plane.  Must find
c: contour crossing farther up:
         call cut_cross(k0,rr0)
cc       print *,'cut_cross: ',k0,(rr0(j,4),j=1,2)
         if(iidone .eq. 1) then
            nmode=nmode - 1
            print *,'no crossing found '
            return
         elseif(iifail .eq. 1) then
            nmode=nmode - 1
            print *,'Failure in cut_cross from fix_path: ',
     .         nmode,k0,rr0
            call stop_run
            return
         endif
      endif
      if(iimt .ne. 0) call mode_traj(k0,rr0,0)
c: Set iimst so that we go to right (backing up) if phase is < pi/2:
      iimst=-1
c
      return
      end
ccc
      subroutine cut_cross(kmid,rrmid)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 kcub,ntry,nwt,jpt,jpt2,iilo
      complex*16 kmid,rrmid(3,4),khi,rrhi(3,4)
      real*8 kilo,kihi,maglo,maghi,dmaglo,dmaghi,c1,c2,c3,c4,delk,kx,
     .   k_real,lnmid,kimid,magmid,dmagmid,wtlo(9),wthi(9),kimid2
      data nwt/9/
      data wtlo/.5,.75,.25,.875,.125,.9375,.0625,.01,.99/
      data wthi/.5,.25,.75,.125,.875,.0625,.9375,.99,.01/
c
c: Found contour, but heading to right in complex k-plane.  Must find
c: contour crossing farther up:
      k_real=dreal(kmid)
      call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
c: Make sure maglo is not slightly positive (cubic root won't work):
      if(maglo .gt. 0.d0) maglo=-1.d-10
c: Find maximum imaginary value for k of interest:
      if(rmin .ge. 999.) then
c: Use .1 km as arbitrary rmin here:
         khi=dcmplx(k_real,dble(db_cut)/(8685.9*(.1)))
      elseif(kim_max .lt. 1.d100) then
         khi=dcmplx(k_real,kim_max)
      else
         khi=dcmplx(k_real,dkim)
      endif
      if(dimag(khi) .le. kilo) then
c: Return if higher value of ki not larger than starting value:
         iidone=1
         return
      endif
c
      call r1r2_calc(khi,rrhi,2,0,jjfail)
      if(jjfail .gt. 0) return
      call ki_copy(khi,rrhi,kihi,maghi,dmaghi)
c
      if(maghi .lt. 0.d0) then
c: Value is still negative, so no crossing
         if(dmaghi .lt. 0.d0) then
c: But derivative indicates there could be a maximum with positive value:
c: Sample at several points between, searching for positive values
c: or positive derivatives:
            do jpt=1,nwt
               kimid=wtlo(jpt)*kilo + wthi(jpt)*kihi
               kmid=dcmplx(k_real,kimid)
               call r1r2_calc(kmid,rrmid,2,0,jjfail)
               if(jjfail .gt. 0) return
               call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
               if(magmid .gt. 0.d0) then
                  call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value '
                  goto 20
               elseif(dmagmid .gt. 0.d0) then
c: Found positive derivative between kimid and kihi, but still 
c: need positive value:
                  do jpt2=1,nwt
                     kimid2=wtlo(jpt2)*kimid + wthi(jpt2)*kihi
                     kmid=dcmplx(k_real,kimid2)
                     call r1r2_calc(kmid,rrmid,2,0,jjfail)
                     if(jjfail .gt. 0) return
                     call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
                     if(magmid .gt. 0.d0) then
                        call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value after pos deriv'
                        goto 20
                     endif
                  enddo
cc    print *,'cut_cross found pos deriv, but not pos value'
                  iidone=1
                  return
               endif
            enddo
            iidone=1
            return
         else
            iidone=1
            return
         endif
      endif
c
20    continue
      ntry=0
cccc  iilo=0
      iilo=1
10    continue
      ntry=ntry + 1
      if(iilo .eq. 0) then
c: On first try, try to get away from branch point:
         kmid=dcmplx(k_real,0.999*kilo+.001*kihi)
      elseif(ntry .gt. 30) then
         print *,'Failure in cut_cross'
         iifail=1
         return
      elseif(ntry .lt. 18 .and. mod(ntry,3) .ne. 0) then
         call cub_fit_new(kilo,kihi,maglo,maghi,dmaglo,dmaghi,
     .      c1,c2,c3,c4,delk)
         call cub_root_new(c1,c2,c3,c4,kx,kcub)
         if(kcub .eq. 0) print *,'kcub = 0 from cut_cross'
         kmid=dcmplx(k_real,kilo + kx*delk)
      else
         kmid=dcmplx(k_real,0.5d0*(kilo+kihi))
      endif
      call r1r2_calc(kmid,rrmid,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnmid=dreal(rrmid(1,4))
      if(dabs(lnmid) .lt. .05d0) then
cc       print *,'cut_cross succcessful: ',ntry
         return
      elseif(lnmid .lt. 0.d0) then
         call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
         iilo=1
      else
         call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
      endif
      goto 10
c
      end
ccc
      subroutine ki_copy(klo,rrlo,kilo,maglo,dmaglo)
c
c: Copies real values from complex values for use in cut_cross.
c
      implicit none
      complex*16 klo,rrlo(3,4)
      real*8 kilo,maglo,dmaglo
c
      kilo=dimag(klo)
      maglo=dreal(rrlo(1,4))
      dmaglo=-dimag(rrlo(2,4))
c
      return
      end
ccc
      subroutine move_up_cut2(k0,rr0,kw,emagmax,iifail,nctot,iidiag,
     .           jjfail)
c
      implicit none
      integer*4 iifail,nctot,ntry,ntry2,iidiag,jjfail
      complex*16 k0,rr0(3,4),kneg,kpos,rrneg(3,4),rrpos(3,4),kneg2
      real*8 kw,emagmax,delki,delk
c
c: Move up in complex plane until |R1*R2| has pos gradient:
cxx   delki=.05d0*kw
      delki=.01d0*kw
      kneg=k0
5     call r1r2_calc(kneg,rrneg,2,0,jjfail)
      if(jjfail .gt. 0) return
      if(dreal(rrneg(1,4)) .gt. 0.d0) then
         delk=.1*delki
         ntry2=0
54       ntry=0
55       ntry=ntry+1
         kneg2=dcmplx(dreal(kneg),dimag(kneg)+delk)
         call r1r2_calc(kneg2,rrneg,2,0,jjfail)
         if(jjfail .gt. 0) return
         if(ntry .gt. 8) then
            kneg=dcmplx(dreal(kneg),dimag(kneg)+.1*delki)
            ntry2=ntry2 + 1
            if(ntry2 .gt. 20) then
               print *,'Failure to pick up path in move_up_cut2: ',
     .            k0/kw,kneg/kw
               iifail=1
               return
            endif
            goto 54
         elseif(dreal(rrneg(1,4)) .gt. 0.d0) then
            delk=0.5d0*delk
            goto 55
         endif
         kneg=kneg2
         if(iidiag .ge. 2) print *,'Pos re(ln) done: ',k0/kw,
     .      kneg/kw,rrneg(1,4),rrneg(2,4)
      endif
c
      kpos=k0
      ntry=0
10    kpos=dcmplx(dreal(kpos),dimag(kpos) + delki)
      ntry=ntry + 1
      if(ntry .gt. 50) then
         print *,'Failure to find kpos in move_up_cut2: ',k0/kw,kpos/kw
         iifail=1
         return
      endif
      call r1r2_calc(kpos,rrpos,2,0,jjfail)
      if(jjfail .gt. 0) return
      if(dreal(rrpos(1,4)) .lt. 0.d0) then
         delki=1.25*delki
         goto 10
      endif
cxx   print *,'move_up2 start: ',kneg/kw,rrneg(1,4),rrneg(2,4),
cxx  .   kpos/kw,rrpos(1,4),rrpos(2,4)
      ntry=0
c
c: We now have kneg with neg Re[ln(r1r2)], kpos with pos Re[ln(r1r2)]:
20    k0=dcmplx(dreal(k0),.5d0*(dimag(kneg)+dimag(kpos)))
      ntry=ntry + 1
      if(ntry .gt. 20) then
         print *,'Unsuccessful move_up_cut2: ',
     .         kneg/kw,rrneg(2,4),kpos/kw,rrpos(2,4)
         iifail=1
         return
      endif
      call r1r2_calc(k0,rr0,2,0,jjfail)
      if(jjfail .gt. 0) return
c: Note that change in Re[f(k)] along imaginary k axis is -IM[df/dk]:
      if(abs(dreal(rr0(1,4))) .lt. emagmax) then
         if(dimag(rrpos(2,4)) .lt. 0.d0 .and.
     .      dimag(rrneg(2,4)) .lt. 0.d0) then
cxx   print *,'cut2 done: ',rrpos(1,4),rrneg(1,4),
cxx  .   kpos/kw,kneg/kw,k0/kw,rr0(1,4),emagmax
            return
         else
            if(iidiag .ge. 2) then
               print *,'emagmax OK, but derivatives wrong sign: ',
     .            kneg/kw,rrneg(2,4),kpos/kw,rrpos(2,4)
            endif
         endif
      endif
      if(dreal(rr0(1,4)) .gt. 0.d0) then
         kpos=k0
         rrpos(1,4)=rr0(1,4)
         rrpos(2,4)=rr0(2,4)
      else
         kneg=k0
         rrneg(1,4)=rr0(1,4)
         rrneg(2,4)=rr0(2,4)
      endif
cxx   print *,'move_up2 loop: ',kneg/kw,rrneg(1,4),rrpos(2,4),
cxx  .   kpos/kw,rrpos(1,4),rrpos(2,4)
      goto 20
c
      end
ccc
      subroutine cut_cross2(kmid,rrmid)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 kcub,ntry,nwt,jpt,jpt2
      complex*16 kmid,rrmid(3,4),khi,rrhi(3,4)
      real*8 kilo,kihi,maglo,maghi,dmaglo,dmaghi,c1,c2,c3,c4,delk,kx,
     .   k_real,lnmid,kimid,magmid,dmagmid,wtlo(9),wthi(9),kimid2
      data nwt/9/
      data wtlo/.5,.75,.25,.875,.125,.9375,.0625,.01,.99/
      data wthi/.5,.25,.75,.125,.875,.0625,.9375,.99,.01/
c
c: Found contour, but heading to right in complex k-plane.  Must find
c: contour crossing farther up:
      k_real=dreal(kmid)
      call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
c: Make sure maglo is not slightly positive (cubic root won't work):
ccc   if(maglo .gt. 0.d0) maglo=-1.d-10
c: Find maximum imaginary value for k of interest:
      if(rmin .ge. 999.) then
c: Use .1 km as arbitrary rmin here:
         khi=dcmplx(k_real,dble(db_cut)/(8685.9*(.1)))
      elseif(kim_max .lt. 1.d100) then
         khi=dcmplx(k_real,kim_max)
      else
         khi=dcmplx(k_real,dkim)
      endif
      if(dimag(khi) .le. kilo) then
c: Return if higher value of ki not larger than starting value:
         iidone=1
         return
      endif
c
      call r1r2_calc(khi,rrhi,2,0,jjfail)
      if(jjfail .gt. 0) return
      call ki_copy(khi,rrhi,kihi,maghi,dmaghi)
c
      if(maghi .lt. 0.d0) then
c: Value is still negative, so no crossing
         if(dmaghi .lt. 0.d0) then
c: But derivative indicates there could be a maximum with positive value:
c: Sample at several points between, searching for positive values
c: or positive derivatives:
            do jpt=1,nwt
               kimid=wtlo(jpt)*kilo + wthi(jpt)*kihi
               kmid=dcmplx(k_real,kimid)
               call r1r2_calc(kmid,rrmid,2,0,jjfail)
               if(jjfail .gt. 0) return
               call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
               if(magmid .gt. 0.d0) then
                  call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value '
                  goto 20
               elseif(dmagmid .gt. 0.d0) then
c: Found positive derivative between kimid and kihi, but still 
c: need positive value:
                  do jpt2=1,nwt
                     kimid2=wtlo(jpt2)*kimid + wthi(jpt2)*kihi
                     kmid=dcmplx(k_real,kimid2)
                     call r1r2_calc(kmid,rrmid,2,0,jjfail)
                     if(jjfail .gt. 0) return
                     call ki_copy(kmid,rrmid,kimid,magmid,dmagmid)
                     if(magmid .gt. 0.d0) then
                        call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
cc    print *,'cut_cross found pos value after pos deriv'
                        goto 20
                     endif
                  enddo
cc    print *,'cut_cross found pos deriv, but not pos value'
                  iidone=1
                  return
               endif
            enddo
            iidone=1
            return
         else
            iidone=1
            return
         endif
      endif
c
20    continue
      ntry=0
10    continue
      ntry=ntry + 1
      if(ntry .eq. 1) then
c: On first try, try to get away from branch point:
         kmid=dcmplx(k_real,0.98*kilo+.02*kihi)
      elseif(ntry .gt. 25) then
         print *,'Failure in cut_cross'
         iifail=1
         return
      elseif(mod(ntry,3) .ne. 0) then
         call cub_fit_new(kilo,kihi,maglo,maghi,dmaglo,dmaghi,
     .      c1,c2,c3,c4,delk)
         call cub_root_new(c1,c2,c3,c4,kx,kcub)
         if(kcub .eq. 0) print *,'kcub = 0 from cut_cross'
         kmid=dcmplx(k_real,kilo + kx*delk)
      else
         kmid=dcmplx(k_real,0.5d0*(kilo+kihi))
      endif
      call r1r2_calc(kmid,rrmid,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnmid=dreal(rrmid(1,4))
      if(dabs(lnmid) .lt. .05d0) then
cc       print *,'cut_cross succcessful: ',ntry
         return
      elseif(lnmid .lt. 0.d0) then
         call ki_copy(kmid,rrmid,kilo,maglo,dmaglo)
      else
         call ki_copy(kmid,rrmid,kihi,maghi,dmaghi)
      endif
      goto 10
c
      end
      subroutine freq_init
c
c: Initializes variables associated with a new frequency.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 j,ii,jsr,inc,ii1,ii2
      data pie/3.14159265358979d0/,twpie/6.28318530717959d0/
c
      w=twpie*f_hz
      wsq=w*w
c: Compute minimum Re(k) to search for eigenvalues:
      if(cphmax .gt. 0.) then
         kremin=w/cphmax
      else
         kremin=-1.d100
      endif
      kcrmin=w/crmax
      khspmin=w/chspmax
c
      if(nthorpe .gt. 0) call thorpe_attn
c
c: Set Airy halfspace thicknesses to 3 wavelengths:
      if(ihf(1) .eq. 1) h(1)=3.d0*geo(2,1,1)/f_hz
      if(ihf(nlay) .eq. 1) h(nlay)=3.d0*geo(1,1,nlay)/f_hz
c
      do j=1,nlay
         call k_attn(w,geo(1,1,j),geo(1,4,j),xk(1,j),xksq(1,j))
         call k_attn(w,geo(2,1,j),geo(2,4,j),xk(2,j),xksq(2,j))
         call eta_calc(xksq(1,j),xksq(2,j),h(j),eta(j),etasq(j),
     .      isp(j),j,nlay)
         if(iisol(j) .eq. 1) then
            call k_attn2(w,geo(1,2,j),geo(1,5,j),xb(1,j),xbsq(1,j),
     .         xbsqinv(1,j))
            call k_attn2(w,geo(2,2,j),geo(2,5,j),xb(2,j),xbsq(2,j),
     .         xbsqinv(2,j))
            call eta_calc(xbsq(1,j),xbsq(2,j),h(j),etb(j),etbsq(j),
     .         iss(j),j,nlay)
         endif
      enddo
c
      kw0=w/cfmin
      xkref=xk(isvmin,nsvmin)
c: Set branch point ratio for reference depth:
      xkrat_ref(1)=xkref/kw0
      kw=dreal(xkref)
      cref=w/kw0
      errdkms=dmin1(errdkms,(1.d-5*kw)**2)
      errdk100=-100.d0*dsqrt(errdkms)
c
      do ii=1,2 
c: xkbp(iihs,iips) holds the branch points for halfspace iihs (1=lower,
c: 2=upper) and wave type iips (1=p-wave,2=s-wave):
         ii1=ii
         ii2=3 - ii
         xkbp(ii,1)=xk(ii1,jsol(ii,2)+jflu(ii,3))
         xkbp(ii,2)=xb(ii1,jsol(ii,2)+jflu(ii,3))
         xkrat(ii,1)=xkbp(ii,1)/kw0
         xkrat(ii,2)=xkbp(ii,2)/kw0
         if(iidiag .ge. 2) then
            print *,'BP = ',ii,xkbp(ii,1),xkbp(ii,2)
            print *,'k/kw = ',ii,xkrat(ii,1),xkrat(ii,2)
         endif
         inc=jflu(ii,3)
         if(allf(ii) .eq. 0) then
            do jsr=1,nzsr
               ksm2_sr(jsr)=(0.d0,0.d0)
               if(kksh(jsr) .eq. 1) ksm2_sr(jsr)=(cs_sr(jsr)**2)/wsq
            enddo
         endif
      enddo
c
      return
      end
ccc
      subroutine k_attn(w,c0,alpha,xk,xksq)
c
c: Computes complex wavenumber for medium with sound speed c0 and
c: attenuation alpha (in dB/(m-kHz)).
      implicit none
      complex*16 xk,xksq
      real*8 w,c0,alpha,fac
c: fac=2*pi*1000*20*log(e):
      data fac/5.457505415367364d+04/
c
      if(c0 .ne. 0.) then
c: Convert alpha (db/m-kHz) to alphap (db/m):
         xk=(w/c0)*dcmplx(1.d0,alpha*c0/fac)
         xksq=xk*xk
      endif
c
      return
      end
ccc
      subroutine k_attn2(w,c0,alpha,xk,xksq,xksqinv)
c
c: Computes complex wavenumber for medium with sound speed c0 and
c: attenuation alpha (in dB/(m-kHz)).
      implicit none
      complex*16 xk,xksq,xksqinv
      real*8 w,c0,alpha,fac
c: fac=2*pi*1000*20*log(e):
      data fac/5.457505415367364d+04/
c
      if(c0 .ne. 0.) then
c: Convert alpha (db/m-kHz) to alphap (db/m):
         xk=(w/c0)*dcmplx(1.d0,alpha*c0/fac)
         xksq=xk*xk
         xksqinv=1.d0/xksq
      else
         xksqinv=(0.d0,0.d0)
      endif
c
      return
      end
ccc
      subroutine eta_calc(xk1sq,xk2sq,h,eta,etasq,iso,jlay,nlay)
c
      implicit none
      integer*4 iso,jlay,nlay
      complex*16 xk1sq,xk2sq,eta,etasq,dk,eip23
      real*8 h,third,sqrt3,pie23
      data third/0.333333333333333d0/,sqrt3/1.73205080756888/
     .   eip23/(-0.5d0,0.86602540378444d0)/,
     .   pie23/2.09439510239320d0/
c
      if(iso .eq. 1 .or. h .le. 0.d0) then
         eta=(0.d0,0.d0)
         etasq=(0.d0,0.d0)
         return
      endif
      if(jlay .ne. 1 .and. jlay .ne. nlay) then
         dk=xk2sq - xk1sq
         if(dreal(dk) .gt. 0.d0) then
            eta=(dk/h)**third
         else
c: For dk/h in left half plane, make sure eta lies close to neg
c: real axis:
            eta=-((-dk/h)**third)
         endif
      else
c: For Airy halfspaces, make sure line in -eta direction in xi-plane 
c: points between +-pi/3 direction:
         dk=xk2sq - xk1sq
         eta=-((-dk/h)**third)
      endif
      etasq=eta*eta
c
      return
      end
ccc
      subroutine freq_chng
c
c: Changes variables associated with a new frequency after freq_init has 
c: already been called.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 j,ii,jx1,jx,jsr,inc
      real*8 w_old,w_rat,w_ratsq,w_23,w_43,rat23
c
c: If Thorpe layers exist, might as well call freq_init:
      if(nthorpe .gt. 0) then
         call freq_init
         return
      endif
c
      w_old=w
      w=twpie*f_hz
      w_rat=w/w_old
      w_ratsq=w_rat*w_rat
      rat23=2.d0/3.d0
      w_23=w_rat**rat23
      w_43=w_23*w_23
c: Compute minimum Re(k) to search for eigenvalues:
      if(cphmax .gt. 0.) then
         kremin=w/cphmax
      else
         kremin=-1.d100
      endif
      kcrmin=w/crmax
c
      do j=1,nlay
         call xk_chng(xk(1,j),xksq(1,j),w_rat,w_ratsq,eta(j),
     .      etasq(j),w_23,w_43,isp(j),h(j),ihf(j))
         if(iisol(j) .eq. 1) then
            call xk_chng2(xb(1,j),xbsq(1,j),w_rat,w_ratsq,etb(j),
     .         etbsq(j),w_23,w_43,iss(j),xbsqinv(1,j),ihf(j))
         endif
      enddo
c
      xkref=xk(isvmin,nsvmin)
      kw=dreal(xkref)
      kw0=w/cfmin
      errdkms=dmin1(errdkms,(1.d-5*kw)**2)
      errdk100=-100.d0*dsqrt(errdkms)
c
      do ii=1,2 
c: xkbp(iihs,iips) holds the branch points for halfspace iihs (1=lower,
c: 2=upper) and wave type iips (1=p-wave,2=s-wave):
         xkbp(ii,1)=xk(ii,jsol(ii,2)+jflu(ii,3))
         xkbp(ii,2)=xb(ii,jsol(ii,2)+jflu(ii,3))
         xkrat(ii,1)=xkbp(ii,1)/kw0
         xkrat(ii,2)=xkbp(ii,2)/kw0
         inc=jflu(ii,3)
         if(allf(ii) .eq. 0) then
            do j=jsol(ii,1),jsol(ii,2)+inc,inc
               jx1=jzmx(j)
               do jx=jx1,jx1+nzmx(j)-1
                  jsr=jsrmx(jx)
                  ksm2_sr(jsr)=ksm2_sr(jsr)/w_ratsq
               enddo
            enddo
         endif
      enddo
c
      return
      end
ccc
      subroutine xk_chng(xk,xksq,w_rat,w_ratsq,eta,etasq,w_23,
     .   w_43,isp,h,ihf)
c
      implicit none
      integer*4 isp,ihf
      complex*16 xk(2),xksq(2),eta,etasq
      real*8 w_rat,w_ratsq,w_23,w_43,h
c
      xk(1)=xk(1)*w_rat
      xksq(1)=xksq(1)*w_ratsq
      if(isp .eq. 0) then
         xk(2)=xk(2)*w_rat
         xksq(2)=xksq(2)*w_ratsq
         if(ihf .eq. 0) then
            eta=eta*w_23
            etasq=etasq*w_43
         else
            eta=eta*w_rat
            etasq=etasq*w_ratsq
c: Layer thickness h inversely proportional to frequency:
            h=h/w_rat
         endif
      else
         xk(2)=xk(1)
         xksq(2)=xksq(1)
c: Layer thickness h inversely proportional to frequency:
         if(ihf .ne. 0) h=h/w_rat
      endif
c
      return
      end
ccc
      subroutine xk_chng2(xk,xksq,w_rat,w_ratsq,eta,etasq,w_23,
     .   w_43,isp,xbsqinv,ihf)
c
      implicit none
      integer*4 isp,ihf
      complex*16 xk(2),xksq(2),eta,etasq,xbsqinv(2)
      real*8 w_rat,w_ratsq,w_23,w_43
c
      xk(1)=xk(1)*w_rat
      xksq(1)=xksq(1)*w_ratsq
      xbsqinv(1)=xbsqinv(1)/w_ratsq
      if(isp .eq. 0) then
         xk(2)=xk(2)*w_rat
         xksq(2)=xksq(2)*w_ratsq
         xbsqinv(2)=xbsqinv(2)/w_ratsq
         if(ihf .eq. 0) then
            eta=eta*w_23
            etasq=etasq*w_43
         else
            eta=eta*w_rat
            etasq=etasq*w_ratsq
         endif
      else
         xk(2)=xk(1)
         xksq(2)=xksq(1)
         xbsqinv(2)=xbsqinv(1)
      endif
c
      return
      end
ccc
      subroutine thorpe_attn
c
c: Set attenuation in Thorpe attenuation layers:
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      real*8 fsq,attn_thorpe,a_nu_db_kyd
      integer*4 j,jlay
c
      fsq=(1.d-3*f_hz)**2
      do j=1,nthorpe
         jlay=jthorpe(j)
	 a_nu_db_kyd=.1*fsq/(1.+fsq)+40.*fsq/(4100.+fsq)+2.75e-4*fsq
	 attn_thorpe=a_nu_db_kyd*(1./914.4)/(.001*f_hz)
         geo(1,4,jlay)=attn_thorpe
         geo(2,4,jlay)=attn_thorpe
      enddo
c
      return
      end
      subroutine h_space(ii,isp,xkh,xkhsq,xksq,eta,etasq,gami,
     .   w,Vmatx,ailay,zetalay,ndv,iiw,jj1,jj2,ihf,xi)
c
c: Sets the plane wave reflection coefficients for a halfspace to zero.
c
      implicit none
      include 'scairy_com'
      integer ii,isp,ndv,iiw,jj1,jj2,ihf,j
      complex*16 xkh,xkhsq,xksq,eta,etasq,gami(3),Vmatx(3,5),ailay(2),
     .   zetalay,zzero,xi,ai1,aip1,zzexp1,
     .   F,G,num,den,xi_k,xi_w,F_dot,G_dot,hdensq,eta_w
      real*8 w,magsq
      data zzero/(0.d0,0.d0)/
c
      if(isp .eq. 0) then
c: Airy halfspace:
         xi=(xkhsq - xksq)/etasq
         call airy_only (xi,ai1,aip1,zzexp1)
cc       call clairy(xi,1,ai1,bi1,aip1,bip1,zzexp1)
cc    if(magsq(ai1x*exp(-zzexp1x)-ai1*exp(-zzexp1)) .gt. 1.d-16) then
cc       print *,'airy_only problem: ',ai1x*exp(zzexp1x),
cc   .      ai1*exp(zzexp1)
cc    endif
c     print *,'clairy   : ',xi,ai1,aip1,zzexp1
c
         if(iiw .eq. 1) then
            ailay(1)=ai1
            ailay(2)=aip1
            zetalay=zzexp1
         endif
c
         F=gami(1)*ai1
         G=eta*aip1
         num=F + G
         den=F - G
         Vmatx(1,jj1)=num/den
         Vmatx(1,jj2)=zzero
         if(ndv .gt. 1) then
            xi_k=2.d0*xkh/etasq
            F_dot=gami(1)*aip1*xi_k + gami(2)*ai1
            G_dot=eta*xi*ai1*xi_k
            hdensq=0.5d0*den*den
            Vmatx(2,jj1)=(F*G_dot - G*F_dot)/hdensq
            Vmatx(2,jj2)=zzero
            if(ndv .gt. 2) then
               if(ihf .eq. 0) then
                  eta_w=eta/(1.5d0*w)
                  xi_w=(xksq + 2.d0*xkhsq)/(-1.5d0*w*etasq)
               else
c: Frequency-dependent gradient in halfspace:
                  eta_w=eta/w
                  xi_w=-2.d0*xkhsq/(w*etasq)
               endif
c: More general xi_w for different eta_w:
cc             xi_w=-2.d0*(xksq/w + eta*eta_w*xi)/etasq
               F_dot=gami(1)*aip1*xi_w + gami(3)*ai1
               G_dot=eta*xi*ai1*xi_w + eta_w*aip1
               Vmatx(3,jj1)=(F*G_dot - G*F_dot)/hdensq
               Vmatx(3,jj2)=zzero
            endif
         endif
      else
c: Homogeneous halfspace:
         do j=1,ndv
            Vmatx(j,jj1)=zzero
            Vmatx(j,jj2)=zzero
         enddo
      endif
      Vmatx(1,5)=zzero
c
      return
      end
      subroutine k_plane(nre,nim,mat1,mat2,mat3,mat4)
c
c: Computes and outputs HDF files for complex k-plane plots of the real
c: and imaginary parts of ln(R1*R2) and/or ln(1-R1*R2).
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
c
      integer*4 nre,nim,jr,ji,iibad
      real*4 mat1(nre,nim),mat2(nre,nim),mat3(nre,nim),mat4(nre,nim),
     .   facr,faci,fac_kr,fac_ki
      complex*16 rr0(3,4),lnr1r2,lnW,psi_fac,arg,H0,kernel,psi_fac_neg
      real*8 pienorm,r_pbli,pierad
c
      call mem_lim(nre,NSRMAX,MLINE,LML,'nreal',5,'NSRMAX',6,iibad,0)
      call mem_lim(nim,NSRMAX,MLINE,LML,'nimag',5,'NSRMAX',6,iibad,0)
      call mem_lim(nre*nim,NTLMAX,MLINE,LML,'nre*nim',7,'NTLMAX',6,
     .   iibad,1)
c
      f_hz=fkpl
      call freq_init
c: Don't allow sheet changes at branch cuts:
      iich=0
      iich_ref=0
c
      iish0(1,1)=iishp
      iish0(2,1)=iishp
      iish0(1,2)=iishs
      iish0(2,2)=iishs
      iishr0(1)=1
      iishr0(2)=1
cc    print *,'enter iishref:'
cc    read(5,*) iishr0(1),iishr0(2)
      call sheet_init(xkh,0,iish0,iishr0)
      if(kduc .eq. 0) then
         kduct=indx_duct(nduct)
      else
         if(kduc .gt. nduct .or. kduc .lt. 1) then
            print *,'Illegal duct number: kduc,nduct = ',kduc,nduct
            return
         endif
         kduct=kduc
      endif
      call zref_chng
      print *,'For k-plane plots, kduct,z,c = ',kduct,
     .   zduct(kduct),geo(isvmin,1,nsvmin)
c
      pierad=pie/180.
      if(iivar .eq. 1) then
         pienorm=pierad
      else
         pienorm=1.d0
      endif
c: For k/kw0, use absolute minimum in SVP, rather than other duct if requested:
      if(iiform .eq. 1) t hen
         fac_kr=kw0
         fac_ki=kw0
      elseif(iiform .eq. 2) then
         fac_kr=1.
         fac_ki=1.
      elseif(iiform .eq. 3) then
         fac_kr=kw0
         fac_ki=1./8685.9
      endif
c
      facr=(xkhr2-xkhr1)/max(nre-1,1)
      do jr=1,nre
         xkhr(jr)=xkhr1 + (jr-1)*facr
      enddo
      faci=(xkhi2-xkhi1)/max(nim-1,1)
      do ji=1,nim
         xkhi(ji)=xkhi1 + (ji-1)*faci
      enddo
      if(iikpl .eq. 4 .or. iikpl .eq. 5) then
         if(rmin .eq. 0.) then
            print *,'rmin = 0 in PBLI for k-plane.  Enter rmin in km: '
            read *,rmin
         endif
         r_pbli=1000.*rmin
      endif
c
      do ji=1,nim
         do jr=1,nre
            xkh=dcmplx(fac_kr*xkhr(jr),fac_ki*xkhi(ji))
            call r1r2_calc(xkh,rr0,iivar,1,jjfail)
            if(jjfail .gt. 0) return
            if(iikpl .eq. 1 .or. iikpl .eq. 3) then
               lnr1r2=rr0(iivar,4)
               mat1(jr,ji)=real(lnr1r2)
               mat2(jr,ji)=dimag(lnr1r2)/pienorm
            elseif(iikpl .eq. -1 .or. iikpl .eq. -3) then
               lnr1r2=rr0(iivar,3)
               mat1(jr,ji)=real(lnr1r2)
               mat2(jr,ji)=dimag(lnr1r2)
            endif
            if(iikpl .eq. 2 .or. iikpl .eq. 3) then
               lnW=cdlog(1.d0 - rr0(iivar,3))
               mat3(jr,ji)=real(lnW)
               mat4(jr,ji)=dimag(lnW)/pierad
            elseif(iikpl .eq. -2 .or. iikpl .eq. -3) then
               lnW=1.d0 - rr0(iivar,3)
               mat3(jr,ji)=real(lnW)
               mat4(jr,ji)=dimag(lnW)
            endif
c: Plot integrand kernel for looking at Pekeris branch line integral
            if(iikpl .eq. 4) then
               psi_fac=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
     .            (gamiref*(1.d0 - rr0(1,3)))
               arg=xkh*r_pbli
               if(cdabs(arg) .gt. 5.) then
                  H0=cdsqrt(dcmplx(0.d0,-2.d0)/(pie*xkh*r_pbli))*
     .               cdexp(dcmplx(0.,1.)*xkh*r_pbli)
               else
                  call cdhankel(arg,1.d-6,H0)
               endif
               kernel=cdlog(0.5*psi_fac*H0*xkh)
               mat1(jr,ji)=real(kernel)
               mat2(jr,ji)=dimag(kernel)/pierad
            elseif(iikpl .eq. 5) then
               psi_fac=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
     .            (gamiref*(1.d0 - rr0(1,3)))
               arg=xkh*r_pbli
               if(cdabs(arg) .gt. 5.) then
                  H0=cdsqrt(dcmplx(0.d0,-2.d0)/(pie*xkh*r_pbli))*
     .               cdexp(dcmplx(0.,1.)*xkh*r_pbli)
               else
                  call cdhankel(arg,1.d-6,H0)
               endif
               kernel=0.5*psi_fac*H0*xkh
               iish0(1,1)=-iish0(1,1)
               call sheet_init(xkh,0,iish0,iishr0)
               call r1r2_calc(xkh,rr0,iivar,1,jjfail)
               if(jjfail .gt. 0) return
               psi_fac_neg=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
     .            (gamiref*(1.d0 - rr0(1,3)))
               kernel=cdlog(0.5*H0*xkh*(psi_fac_neg - psi_fac))
               mat1(jr,ji)=real(kernel)
               mat2(jr,ji)=dimag(kernel)/pierad
               iish0(1,1)=-iish0(1,1)
               call sheet_init(xkh,0,iish0,iishr0)
            endif
         enddo
      enddo
      if(iikpl .gt. 0 .and. iiwr .eq. 2) then
         if(iikpl .eq. 1 .or. iikpl .ge. 3) then
            do ji=1,nim
               do jr=1,nre
                  if(mat2(jr,ji) .lt. 0.) mat2(jr,ji)=mat2(jr,ji) + 360.
               enddo
            enddo
         endif
         if(iikpl .eq. 3) then
            do ji=1,nim
               do jr=1,nre
                  if(mat4(jr,ji) .lt. 0.) mat4(jr,ji)=mat4(jr,ji) + 360.
               enddo
            enddo
         endif
      endif
      if(abs(iikpl) .eq. 1 .or. abs(iikpl) .eq. 3) then
         call out_writex(outroot,lout,SUFX//'re',3,mat1,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'im',3,mat2,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Im[ln(R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
      endif
      if(abs(iikpl) .eq. 2 .or. abs(iikpl) .eq. 3) then
         call out_writex(outroot,lout,SUFX//'reW',4,mat3,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(1-R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'imW',4,mat4,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Im[ln(1-R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
      endif
      if(abs(iikpl) .eq. 4) then
         call out_writex(outroot,lout,SUFX//'lnker',6,mat1,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(kernel)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'argker',7,mat2,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Im[ln(kernel)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
      endif
      if(abs(iikpl) .eq. 5) then
         call out_writex(outroot,lout,SUFX//'lnkerdif',9,mat1,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(kernel)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'argkerdif',10,mat2,xkhi,
     .      xkhr,nim,nre,kilab,krlab,0.,0.,0.,0.,2,
     .      'Im[ln(kernel)] vs. k',' ',' ',' ','f7.4','f7.4','f7.2',
     .      ncall)
      endif
c
      return
      end
c: *******************************
c: *   AUTHOR (1996):            *
c: *     P R O S I M  GROUP      *
c: *     S A C L A N T C E N     *
c: *******************************
      subroutine tenv_read_pro(psect,isubx)
c     
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
      include 'sector_env_com'
      include 'debug_com'
c:    Local variables
      integer*4 i,j,j1,psect,isubx,iierr,nline
      character*64 eline
      common /info/ienv_info,lu_env
      common/iterpar/iter,iforwpop,ipop
      integer iter, iforwpop, ipop
      integer lu_env,ienv_info,ii
      real*8 hhb
      common /tiltparm/tiltv,tilth,dtiltv,dtilth
      logical tiltv,tilth
      real dtiltv,dtilth
c
      data eline/'INVALID INPUT IN TENV FILE: '/
c
      nline=0
      j1=1
      geot(2,1,j1)= 343.0
      geot(2,2,j1)= 0.0
      geot(2,3,j1)= 1.21e-3
      geot(2,4,j1)= 0.0
      geot(2,5,j1)= 0.0
      ht(j1)=1.d+20
      ktt(j1)= 1
      do j=1,5
         geot(1,j,j1)= geot(2,j,j1)
      enddo
c
      call svp_check_val_lay(2,0.d0,geot(1,1,j1),'upper h-space',13,
     .      eline,nline,iierr)
c
      nsvp=R_ND0(psect)
c
      if(nrec .gt. 1) then
         if(tiltv) then
            do i=1,nrec
               dtiltvp(i)=dtiltv/(nrec-1)*(i-1)
            end do
         else
            do i=1,nrec
               dtiltvp(i)=0.0e0
            end do
         end if
c      else
c         tiltv=.false.
      end if
c
      if(nsrc .gt. 1) then
         if(tilth) then
            dtilthp=dtilth/float(nsrc-1)
            do i=2,nsrc
               zrec(i)=zrec(1)+dtilthp*(i-1)
            end do
c            write(6,*)
c            write(6,*)'Horiz tilt: ',dtilth,zrec(1),zrec(nsrc)
         else
            do i=1,nsrc
               dtilthp=0.0e0
            end do
         end if
c      else
c         tilth=.false.
      end if
c
      if(iidebug .eq. 1)
     . write(6,*)'No of layers in water: ',nsvp

      do j=1,nsvp
         zsvp(j)=R_Z0(j,psect)
         csvp(j)=R_C0(j,psect)
         if(iidebug .eq. 1)
     .        write(6,*)zsvp(j),csvp(j)
      end do
c
      do j=2,nsvp
         call check_val_r8(zsvp(j),zsvp(j-1),1.d+10,eline,27,nline,
     .      'zsvp(j)',7,iierr)
         call check_val_r8(csvp(j),1.d-10,1.d+10,eline,27,nline,
     .      'csvp(j)',7,iierr)
      enddo
c
c This is ONLY for IVWK
cpln      zsvp(nsvp)=r_h0(psect)
cpln      csvp(nsvp)=1495.-0.04*zsvp(nsvp)
cpln      R_C0(nsvp,psect)=csvp(nsvp)
c
      rho_svp=1.0
      alpha_svp=0.0
c
      if (R_ND1(psect).eq.0) then
         nlayb=0
      else
         nlayb=R_ND1(psect)-1
      endif
      if(iidebug .eq. 1)
     . write(6,*)'No of layers in bottom: ',nlayb
c
      do j=1,nlayb
         ktb(j)= 1
         hb(j)=R_Z1(j+1,psect)-R_Z1(j,psect)
cpln     SAGA RETURNS LAYER THICKNESS AND NOT ABSOLUTE
cpln     SEDIMENT DEPTH RE WATER-SEDIMENT INTERFACE
cpln         hb(j)=R_Z1(j+1,psect)
         geob(1,1,j)= R_C1(j,psect) ! P sound speed


cpln Original
         geob(2,1,j)= R_C1(j+1,psect) !P sound speed
cpln ONLY ISO SPEED LAYERING IN BOTTOM:
cpln         geob(2,1,j)= R_C1(j,psect) !P sound speed


c Only for IVWK workshop
cpln         geob(2,1,j)= R_C1(j,psect) !P sound speed
c
cpln         write(6,*)hb(j),geob(1,1,j),geob(2,1,j)
cpln         write(6,*)hb(j),R_C1(j,psect),R_C1(j+1,psect)
cpln ONLY ASCOT01
cpln         if(nlayb.gt.2) then
cpln            geob(2,1,j)= R_C1(j+1,psect) !P sound speed
cpln         else
cpln            geob(2,1,j)= R_C1(j,psect) !P sound speed
cpln         end if
         geob(1,2,j)= 0.0       ! shear speed
         geob(2,2,j)= 0.0       ! shear speed
         geob(1,3,j)= R_R1(psect)
         geob(2,3,j)= R_R1(psect)
C ONLY FOR ITWK
cpln         geob(1,4,j)= R_BETA(j,psect)
cpln         geob(2,4,j)= R_BETA(j,psect)
         geob(1,4,j)= R_BETA(1,psect)
         geob(2,4,j)= R_BETA(1,psect)
         geob(1,5,j)= 0.0       ! shear att
         geob(2,5,j)= 0.0       ! shear att
         if(iidebug .eq. 1) then
            write(6,*)hb(j),geob(1,1,j),geob(2,1,j)
            write(6,*)geob(1,2,j),geob(2,2,j)
            write(6,*)geob(1,3,j),geob(2,3,j)
            write(6,*)geob(1,4,j),geob(2,4,j)
            write(6,*)geob(1,5,j),geob(2,5,j)
         end if
         if(j .eq. 1) then
c: Check for negative h, meaning two-way travel time, and negative c,
c: meaning csed/cwater ratio:
            if(geob(1,1,1) .lt. 0.d0) then
               cs_cw_rat=-geob(1,1,1)
               geob(1,1,1)=cs_cw_rat*csvp(nsvp)
      print *,'cb(1) = ',geob(1,1,1)
            endif
            if(hb(1) .lt. 0.d0) then
c: Nominal average gradient (since we don't know it for sure):
               gbar=0.75
               tau_2way=-hb(1)
               hb(1)=geob(1,1,1)*(exp(gbar*tau_2way/2.) - 1.)/gbar
      print *,'hb(1) = ',hb(1)
            endif
         endif
         call svp_check_val_lay(1,hb(j),geob(1,1,j),'bottom layer',12,
     .      eline,nline,iierr)
         call svp_check_val_lay(2,hb(j),geob(1,1,j),'bottom layer',12,
     .      eline,nline,iierr)
         call zero_sh(geob(1,2,j),geob(1,5,j),j,'top   ','bottom')
         call zero_sh(geob(2,2,j),geob(2,5,j),j,'bottom','bottom')
      end do
cpln      pause
c
c     read subbottom parameters from temp file LUSVP
      j1=nlayb+1
      
      hb(j1)=1.d+20
      ktb(j1)=1
      geob(1,1,j1)= R_C2(psect) ! P sound speed
      geob(2,1,j1)= R_C2(psect) !P sound speed
      geob(1,2,j1)= 0.0         ! shear speed
      geob(2,2,j1)= 0.0         ! shear speed
      geob(1,3,j1)= R_R2(psect)
      geob(2,3,j1)= R_R2(psect)
C ONLY ITWK
cpln      geob(1,4,j1)= R_BETA(j1,psect)
cpln      geob(2,4,j1)= R_BETA(j1,psect)
      geob(1,4,j1)= R_BETA(2,psect)
      geob(2,4,j1)= R_BETA(2,psect)
      geob(1,5,j1)= 0.0         ! shear att
      geob(2,5,j1)= 0.0         ! shear att

cpln      do ii=1,j1
cpln         write(6,*)'R_BETA: ',R_BETA(ii,psect)
cpln      end do
cpln      pause
c
c: c_hsp < 0 means use previous layer as halfspace:
c: Set thickness of halfspaces to large numbers (for use in zmx_init):
      call svp_check_val_lay(1,0.d0,geob(1,1,j1),'lower h-space',
     .     13,eline,nline,iierr)
      call zero_sh(geob(1,2,j1),geob(1,5,j1),j1,'top   ','bottom')
c
      if(iidebug .eq. 1) then
         write(6,*)'Sub-bottom'
         write(6,*)geob(1,1,j1),geob(2,1,j1)
         write(6,*)geob(1,2,j1),geob(2,2,j1)
         write(6,*)geob(1,3,j1),geob(2,3,j1)
         write(6,*)geob(1,4,j1),geob(2,4,j1)
         write(6,*)geob(1,5,j1),geob(2,5,j1)
      endif
c
      if(iiwrite .gt. 0) then
c     
C WRITE OUT ENVIRONMENT TO INPUT FILE FOR STAND ALONE PROSIM
         write(lualon,*)'OUTPUT FROM SAGA'
         write(lualon,'(I5)')0
         write(lualon,'(2I5)')0,0
         write(lualon,'(3I5)')1,0,0
         write(lualon,'(2I5)')0,0
         write(lualon,'(F33.26,2I5,F5.1,I5)')zsvp(nsvp),
     .                 0,0,10.0,0
         do j=1,nsvp
            write(lualon,'(2F33.26)')zsvp(j),csvp(j)
         end do
         hhb=0
         if (R_ND1(psect).gt.0) then
            do j1=1,nlayb
               hhb=hhb+hb(j1)
            end do
            write(lualon,'(3F33.26)')hhb,geob(1,3,nlayb),
     .            geob(1,4,nlayb)
            hhb=0
            write(lualon,'(2F33.26)')hhb,geob(1,1,1)
            do j1=1,nlayb
               hhb=hhb+hb(j1)
               write(lualon,'(2F33.26)')hhb,geob(2,1,j1)
            end do
         else
            write(lualon,'(3F33.26)')hhb,geob(1,3,nlayb+1),
     .            geob(1,4,nlayb+1)
         end if
         j1=nlayb+1
         write(lualon,'(3F33.26)')geob(1,3,j1),geob(1,4,j1),
     .        geob(1,1,j1)
         write(lualon,'(2I5)')0,0
         write(lualon,*)'PULSE'
         write(lualon,'(1I5,3F33.26)')nfftbb,fmin,fmax,-fsbb
         write(lualon,'(2F33.26,i5)')rkm(1)*1.D-3,rkm(nsrc)*1.D-3,nsrc
         write(lualon,'(2F33.26,i5)')zrec(1),zrec(nrec),
     .        nrec
         write(lualon,'(F33.26)')zsrc(1)
         close(lualon)

      end if
c
      nlayt=0
c
      return
      end
c
      subroutine zero_sh(cs,as,j,ch_tb1,ch_tb2)
c
      implicit none
      include 'Parms_com'
      real*8 cs,as
      integer*4 j
      character*6 ch_tb1,ch_tb2
c
      if(cs .eq. 0.d0 .and. as .ne. 0.d0) then
         write(lusvp,125) ch_tb1,ch_tb2,j
125      format('/*** SHEAR ATTENUATION SET TO 0 for 0 shear speed at ',
     .         a6,' of ',a6,' layer # ',i3,'***'/)
         as=0.d0
      endif
c
      return
      end
c
      subroutine svp_check_val_lay(ii,h,geo,char_lay,nch,eline,nline,
     .   iierr)
c
      implicit none
      real*8 h,geo(2,5)
      integer*4 ii,nch,nline,iierr
      character*64 char_lay,eline
c
c: Make negative h mean two-way travel time:
      call check_val_r8(h,-10.0d0,1.5d4,eline,27,nline,
     .   'h '//char_lay(1:nch),2+nch,iierr)
c: Make negative c mean sound speed ratio:
      call check_val_r8(geo(ii,1),-10.d0,1.d+10,eline,27,nline,
     .   'cp '//char_lay(1:nch),3+nch,iierr)
      call check_val_r8(geo(ii,2),0.0d0,1.d+10,eline,27,nline,
     .   'cs '//char_lay(1:nch),3+nch,iierr)
      call check_val_r8(geo(ii,3),1.d-50,1.d+50,eline,27,nline,
     .   'rho '//char_lay(1:nch),4+nch,iierr)
      call check_val_r8(geo(ii,4),-200.d0,999.d0,eline,27,nline,
     .   'ap '//char_lay(1:nch),3+nch,iierr)
      call check_val_r8(geo(ii,5),-200.d0,200.d0,eline,27,nline,
     .   'as '//char_lay(1:nch),3+nch,iierr)
c
      return
      end
      subroutine mode_branch(iiref,k_st,iiccwx,iimstx,iilk0,
     .   iidup,ndup,ndup_max)
c
c: Finds modes along a given branch (|R1R2|=1 contour) starting at k_st
c: in the direction given by iiccwx (1=to left in k plane, -1=to right).
c: For iidup=1, check found modes for duplicated with modes kn(1:nm_ok).
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 nm_ok,iiref,iiccwx,iimstx,iilk0,iidup,iduct,ndup,
     .   ndup_max,nm1
      complex*16 k_st,k0,rr0(3,4),r1r2(3,4)
c
c pln 6/4/00 (problems found for certain environments
c solved by following or by increaseing db_cut to cut weak modes
      if(nmode.gt.NM_MAX) then
         jjfail = 1
         return
      end if
c
      kn(nmode)=k_st
      iiccw=iiccwx
      iimst=iimstx
      iilk=iilk0
      nhigh_max=3
      if(iilk .gt. 0) nhigh_max=1
c
      iifail=0
      iidone=0
      nhigh=0
      nrise=0
      ndup=0
      nm_ok=nm_put
      nm1=0
c
c: Start at k of reference depth and find |R1R2|=1 contour:
cc    kcut=dcmplx((1.d0-1.d-8)*dreal(xkref),dimag(xkref))
      if(iiref .eq. 1) then
c: Starting at reference depth, move slightly to left:
cc       k0=dcmplx((1.d0-1.d-8)*dreal(k_st),dimag(k_st))
         if(nsvmin .eq. nlay) then
            k0=cdsqrt(k_st*k_st - .5d0*etasq(nlay))
            call sheet_init(k0,1,iish,iish_ref)
            if(dreal(k0) .gt. dreal(k_st)) iish_ref(1)=-1
cc          print *,'HSP k0 = ',dreal(k0)/kw0,dimag(k0)*8685.9
         else
c: Move to left by 1/1000 of nominal mode spacing:
            k0=k_st - min(.001d0*pie/Htot,.1d0*dreal(k_st))
            call sheet_init(k0,1,iish,iish_ref)
         endif
      else
c: Starting at point on contour, start on k_st:
         k0=k_st
      endif
      call r1r2_calc(k0,rr0,2,0,jjfail)
      if(jjfail .gt. 0) return
cc    call fix_path(k0,rr0,0)
cc    if(iidone .eq. 1) goto 88
      if(iimt .eq. 1) call mode_traj(k0,rr0,0)
c
      nm1=nmode+1
10    continue
         call eig_findm(k0,rr0,r1r2,nm1)
         if(jjfail.eq.1) return
         if(iifail .eq. 1) then
            if(iiwrite .gt. 0) then
               print *,'Mode finding failure in duct ',kduct,
     .              ' at depth ',zduct(kduct)
               print *,'Continuing with modes found ...'
            end if
            iifail=0
            return
         endif
c: Return if branch line mode found on main contour and BLMs will
c: be found later (for all-fluid when ref depth can be placed in 
c: halfspace):
         if(iiblm .gt. 0) return
         if(iidone .eq. 1) then
c: Subtract off any modes that were too weak at end of mode search:
            if(nhigh .gt. 1) nm_put=nm_put + 1 - nhigh
            return
         endif
         if(iidup .eq. 1) then
cpln               write(6,*)'Enter from mode_branch'
            call duct_dupl(nm_ok,ndup,r1r2,phi,dphi,iduct)
            if(jjfail.eq.1) return
            if(ndup .ge. ndup_max) return
         else
            nm_put=nm_put + 1
         endif
      goto 10
c
      end
      subroutine mfun_fill(phi_,mfun,mfunph,suf,lsuf)
c
c: Outputs Re/Im and/or Mag/Phase of mode function phi_.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm,jzmf,jsr,lsuf
      complex*8 phi_(nzsr,nmode)
      real*4 mfun(nmode,nzmf),mfunph(nmode,nzmf),piedeg
      character*64 suf
c
      piedeg=acos(-1.)/180.
      do jm=1,nmode
         xmode(jm)=jm
      enddo
c
      if(iiri .ne. 0) then
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
            do jm=1,nmode
               mfun(jm,jzmf)=real(phi_(jsr,jm))
               mfunph(jm,jzmf)=aimag(phi_(jsr,jm))
            enddo
         enddo
c      
         if(iiri .eq. 1 .or. iiri .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'re',lsuf+2,
     .         mfun,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         mrlab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
         if(iiri .eq. 2 .or. iiri .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'im',lsuf+2,
     .         mfunph,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         milab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
      endif
      if(iimp .ne. 0) then
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
            do jm=1,nmode
               mfun(jm,jzmf)=abs(phi_(jsr,jm))
               mfunph(jm,jzmf)=atan2(aimag(phi_(jsr,jm)),
     .            real(phi_(jsr,jm)))/piedeg
            enddo
         enddo
         if(iimp .eq. 1 .or. iimp .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'mag',lsuf+3,
     .         mfun,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         malab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
         if(iimp .eq. 2 .or. iimp .eq. 3) then
            call out_writex(outroot,lout,suf(1:lsuf)//'ph',lsuf+2,
     .         mfunph,zmf,xmode,nzmf,nmode,dlab,mnlab,z4,z4,z4,z4,2,
     .         mplab,'m',' ',' ','f5.0','f5.1','f7.2',ncall)
         endif
      endif
c
      return
      end
      subroutine mode_field(phiz,dpsiz,expz_gbs,plcx,tlx,jzs)
c
c: Computes field given mode eigenvalues and mode function values.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jzs,jm,jr,jd,jsr,kk,kk0,jsrc,jrec,msrc,nrmax
      complex*16 knx,iknx,iknr,sqkn,cfac2,phi_src,H0,h_arg,
     .   cfac,exp_sqk
      complex*8 phiz(nzsr,nmode),dpsiz(nzsr,nmode),plcx(nsrc,nrec),tlz
      real*8 expz_gbs(nzsr,nmode),rsig_max,hmsq,xn_b,xn_beam,
     .   magsq,rng_im
      real*4 tlx(nsrc,nrec),tlmagsq
c: cfac = sqrt(2*pi)*exp(i*pi/4) (see p. 7, ORCA II):
      data cfac/(1.77245385090552d0,1.77245385090552d0)/
c
      msrc=mzsrc(jzs)
      if(iigbs .ne. 0) then
         rng_im=-b_gbs(jzs)*cos(th_gbs(jzs)*pie/180.d0)
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .         +dtiltvp(jd),rng_im+dtiltvp(jd)))
            end do
cpln            sq2pir(jr)=cfac/cdsqrt(dcmplx(range(jr),rng_im))
         enddo
      else
         rng_im=0.d0
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .         +dtiltvp(jd),rng_im+dtiltvp(jd)))
            end do
cpln            sq2pir(jr)=cfac/dsqrt(range(jr))
         enddo
      endif
c
      hmsq=25.d0
c: Initialize complex propagation loss array:
      do jsrc=1,nsrc
         do jrec=1,nrec
            plcx(jsrc,jrec)=cmplx(0.e0,0.e0)
         enddo
      enddo
c
      xn_b=(w/dreal(cp_sr(msrc)))*b_gbs(jzs)
      do jm=1,nmode
         knx=kn(jm)
         iknx=dcmplx(-dimag(knx),dreal(knx))
         sqkn=cdsqrt(knx)
c: rsig_max is maximum range at which this mode is significant (50 dB down 
c: from strongest mode):
         rsig_max=kim_fac/dmax1(1.d-20,dimag(knx)-kim_min)
c: Obtain p-wave mode excitation at source depth:
         phi_src=phiz(msrc,jm)/rho_sr(msrc)
         xn_beam=xn_b - expz_gbs(msrc,jm)
c: Obtain p-wave mode excitations at receiver depths:
         do jrec=1,nrec
            jsr=mzrec(jrec)
c: rho included in mode_fun:
            phisr(jsr)=phiz(jsr,jm)
            if(kksh(jsr) .eq. 1) then
c: Add shear wave potential contributions to recs in shear layers (see p. 117):
               phisr(jsr)=phisr(jsr) - 2.d0*ksm2_sr(jsr)*
     .            (knx*knx*phisr(jsr) + iknx*dpsiz(jsr,jm))
            endif
         enddo
c: Find range beyond which mode is insignificant (ranges have been sorted):
         nrmax=0
         call hunt(range,nrng,rsig_max,nrmax)
         do jr=1,nrmax
            kk0=krec_jr(jr)
            do kk=kk0+1,kk0+nrec_jr(jr)
               jsrc=jrec_jr(1,kk)
               jrec=jrec_jr(2,kk)
               h_arg=knx*dcmplx(range(jr)+dtiltvp(jrec),
     .               rng_im+dtiltvp(jrec))
cpln               h_arg=knx*dcmplx(range(jr),rng_im)
               if(magsq(h_arg) .gt. hmsq) then
c     zs: Include normalization by exp(-xn_beam) here:
                  iknr=dcmplx(-dimag(h_arg),dreal(h_arg))
                  exp_sqk=cdexp(dcmplx(dreal(iknr)-xn_beam,dimag(iknr)))
     .                 /sqkn
                  cfac2=sq2pir(jr,jrec)*phi_src*exp_sqk
cpln               cfac2=sq2pir(jr)*phi_src*exp_sqk
               else
                  call cdhankel(h_arg,1.d-6,H0)
c     c             cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
                  cfac2=dcmplx(-pie*dimag(H0),pie*dreal(H0))*phi_src*
     .                 dexp(-xn_beam)
               endif
               jsr=mzrec(jrec)
               plcx(jsrc,jrec)=plcx(jsrc,jrec) + cfac2*phisr(jsr)
            enddo
         enddo
99       continue
      enddo
      do jsrc=1,nsrc
         do jrec=1,nrec
            tlz=plcx(jsrc,jrec)
            tlmagsq=real(tlz)*real(tlz)+aimag(tlz)*aimag(tlz)
            tlx(jsrc,jrec)=10.*alog10(amax1(1.e-37,tlmagsq))
         enddo
      enddo
c
      return
      end
      subroutine mode_find(iiwrt)
c
c: Finds eigenvalues and computes mode functions.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iiwrt,iduct,nm_duct,nm_found,iidup,ndup,ndup_max,
     .   nm_put0,nm_tot0,iidiff,ntry1
      complex*16 k_st,rr_st(3,4),k_jump
c
      ntry1=0
      nmode=0
      nblm=0
      nblm_max=0
      phi_mag_max=0.
c: nm_put is index in kn to put valid modes (excluding duplicates):
      nm_put=nmode
c: nm_tot is total number of modes at this frequency:
      nm_tot=0
c: Flag to change sheets when crossing Pekeris branch cuts:
      iich=1
      iich_ref=1
c
      call freq_init
c
      if(iimt .ne. 0) call mode_traj_bp(rr_st)
c: Initialize kim_min (the min modal attenuation in neper/m), and
c: kim_max (the max attenuation a mode can have to be kept):
      kim_min=.25d0*dkim
      kim_max=kim_min + 5.d0*dkim
c
      phfac0=dmax1(2.d0,dble(phfac))
c: Set ph_step, the phase step initially used along the |R1R2|=1 contour:
      ph_step=twpie/nint(phfac0)
c: Set mag_step, the magnitude step initially used along the phase contour:
      mag_step=0.1d0
c
      iidup=0
c: Loop over ducts in which to set reference depth:
      do iduct=nduct,1,-1
c
         kduct=indx_duct(iduct)
         if(kduct .ne. kduct0) call zref_chng
c
         if(iduct .eq. nduct) then
            call xi_cut_set
            if(jjfail .gt. 0) then
               if(iiwrite .gt. 0)
     .            write(6,*)'Return from xi_cut_set'
               return
            end if
         endif
c: Make space for copying new modes on end of existing list:
         if(iduct .lt. nduct) nmode=nmode + 1
c
         ndup_max=2
         if(nsvmin .eq. nlay) ndup_max=10
c
         nm_tot0=nm_tot
         nm_put0=nm_put
c: Find modes on contour associated with current ref depth:
cpln         write(6,*)'#1 mode_branch'
         call mode_branch(1,xkref,1,0,0,iidup,ndup,ndup_max)
            if(jjfail .gt. 0) then
               if(iiwrite .gt. 0)
     .            write(6,*)'Return from mode_branch'
               return
            end if
         nm_found=nm_put - nm_put0
         nm_duct=nm_found
c
         nm_tot=nm_tot + nm_found
         nmode=nm_tot
c: Merge new modes found in this duct into a sorted list (but leave 
c: branch line modes on end of list:
         if(nm_tot0 .gt. 0 .and. nm_found .gt. 0 .and. 
     .      nsvmin .ne. nlay) then
            call hpsort_indx_c16(nmode,kn(1),kn_indx)
            call eig_sort(nmode,nzsr,kn_indx,eig_char,nzref,
     .         ncalc,iishn,mode_phz,phi,dphi,psi,dpsi,exp_gbs)
         endif
c
         if(iifail .ne. 0) goto 20
c
15       continue
         iidiff=0
         if(iidone .eq. 0 .and. ndup .gt. 1) then
c: When mode search encountered duplicate modes before getting to end of
c: search (Im(kn) limit, e.g.), find contour starting at upper left end
c: and proceed back down and to right on it:
cpln            write(6,*)'Enter from mode_find #1',iduct
            call contour_find(k_jump,iduct,k_st,rr_st,iidiff)
cpln            write(6,*)'Exit contour_find #1',jjfail
            if(jjfail .gt. 0) return
         elseif(iiblm .ne. 0) then
c: Mode search ended by finding a branch line mode, so look for contour
c: in upper left corner of k plane:
            if(iduct .lt. nduct) then
cpln               write(6,*)'Enter contour_find #2',iifail,jjfail,
cpln     .                    iduct,nduct
               call contour_find(k_jump,iduct,k_st,rr_st,iidiff)
cpln               write(6,*)'Exit contour_find #2',iifail,jjfail,
cpln     .                    iduct,nduct
cpln               write(6,*)k_jump,k_st,rr_st,iidiff
               if(jjfail .gt. 0) return
            else
               call contour_find2(k_st,rr_st,iidiff)
               if(jjfail .gt. 0) return
            endif
c: EKW FIX 6/2/98:
            iiblm=0
         endif
         if(iidiff .eq. 1) then
c: Different contours found, so look for modes on new contour:
            nm_put0=nm_put
            nmode=nmode + 1
            call iish_xfer(iish,iish_ref,iishx,iish_refx)
c: Find modes to right of k_st:
cpln            write(6,*)'#2 mode_branch'
            call mode_branch(0,k_st,-1,0,1,1,ndup,1)
            if(jjfail.eq.1) return
            if(ndup .eq. 1) then
               if(iiwrite .gt. 0)
     .         print *,'Contour_find looped back on same branch'
               nmode=nm_tot
cpln               write(6,*)'Im going to 20'
               goto 20
            endif
            nmode=nmode + 1
            call sheet_init(k_st,0,iishx,iish_refx)
c: Find modes to left of k_st:
            iiblm=0
cpln            write(6,*)'#3 mode_branch'
            call mode_branch(0,k_st,1,0,1,1,ndup,1)
            if(jjfail.eq.1) return
            nm_found=nm_put - nm_put0
            nm_duct=nm_duct + nm_found
            nm_tot=nm_tot + nm_found
c
            nmode=nm_tot
            if(nm_found .gt. 0) then
c: Merge new modes found in this duct into a sorted list:
               call hpsort_indx_c16(nmode,kn(1),kn_indx)
               call eig_sort(nmode,nzsr,kn_indx,eig_char,nzref,
     .            ncalc,iishn,mode_phz,phi,dphi,psi,dpsi,exp_gbs)
            endif
c: Go back and check for iiblm again if we just did region between p- and
c: s-wave branch points:
cpln 14/4/2000 Can't find modes->continue to goto 15
            ntry1=ntry1+1
            if(ntry1 .gt. 200) then
               jjfail=1
               write(6,*)'Exhausted trying to find mode'
               return
            end if
            goto 15
         else
            nmode=nm_tot
         endif
20       continue
c
         if(iduct .eq. nduct .and. nduct .gt. 1) then
c: After finding modes at first reference depth, choose mode at which 
c: to start looking for different mode branches for other ducts:
            call k_jump_find(nmode,kn,k_jump)
         endif
c
         if(iiwrt .eq. 1) then
            print *,'Informative message: Finished checking duct '//
     .         'at depth ',zduct(kduct)
            print *,'   # modes found in duct = ',nm_duct
            write(lusvp,'(a,f7.2)') 'Informative message: Finished '//
     .         'checking duct at depth ',zduct(kduct)
            write(lusvp,'(a,i4)') '   # modes found in duct = ',nm_duct
         endif
c
199      continue
c
         iidup=1
      enddo
c
      return
      end
ccc
      subroutine k_jump_find(nmode,kn,k_jump)
c
      implicit none
      integer*4 nmode,jm
      complex*16 kn(0:nmode),k_jump,delk
c
c: After first of several ducts, find mode at bottom of vertical section 
c: (first "evanescent" mode) to use to jump to other mode branches:
      k_jump=kn(nmode)
      do jm=nmode,2,-1
         delk=kn(jm) - kn(jm-1)
         if(-real(delk) .gt. dimag(delk)) then
            k_jump=kn(jm)
            return
         endif
      enddo
c
      return
      end
      subroutine mode_fun(kmode,R1x,R2x,dW_dk,phiz,dphiz,psiz,dpsiz,
     .   expz_gbs,jm)
c
c: Finds the mode function whose eigenvalue is kn as a function of depth.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 jm,ii,jflu1,jflu2,inc,jsol1,jsol2,jx1,j,jx,jsr,
     .   ii1,ii2,jjhsp,joff,k,j1
      complex*8 phiz(nzsr,jm),dphiz(nzsr,jm),psiz(nzsr,jm),
     .   dpsiz(nzsr,jm)
      complex*16 kmode,R1x,R2x,dW_dk
      complex*16 phiref,dphiref,phi0,dphi0,psi0,dpsi0,mnorm(2),arg,
     .   Ap,As,phi2,dphi2,psi2,dpsi2
      real*8 expz_gbs(nzsr,jm),sg,Apexp,Asexp
c
c: Include minus sign here so that mode functions turn out real:
c     mnorm(1)=1.d0/cdsqrt(-dW_dk)
c: Include sqrt(k) in normalization & keep track of sheets for sqrt operation:
      arg=-2.d0*kmode/(dW_dk*rho_duct)
c: (-1) factor now determined by dphi at top (see sg):
      mnorm(1)=cdsqrt(arg)
      mnorm(2)=R1x*mnorm(1)
      phiref=mnorm(1)*(1.d0 + R1x)
      dphiref=gamiref*mnorm(1)*(1.d0 - R1x)
      do ii=1,2
         ii1=ii
         ii2=3 - ii
         jflu1=jflu(ii,1)
         jsol1=jsol(ii,1)
         inc=jflu(ii,3)
         joff=jflu(ii,4) - inc
         jjhsp=jhsp(ii)
         if(allf(ii) .eq. 1) then
c: If all fluid, include halfspace by adding inc:
            jflu2=jflu(ii,2) + inc
            jsol2=jsol(ii,2)
         else
            jflu2=jflu(ii,2)
c: If any solid layers, include halfspace by adding inc:
            jsol2=jsol(ii,2)+inc
         endif
c: Limit counters to layers in which receivers are located:
c: Don't limit counters now (for mode id and -1 norm fac):
cxx      if(ii .eq. 1) then
cxx         jflu2=min(jflu2,jlmax)
cxx         jsol2=min(jsol2,jlmax)
cxx      else
cxx         jflu2=max(jflu2,jlmin)
cxx         jsol2=max(jsol2,jlmin)
cxx      endif
c
c: Initialize phi,dphi at reference depth:
         phi0=phiref
c: Make derivative continuous by multiplying dphi0 by -1 for upgoing 
c: reflection coefficient (see p. 104-105):
         dphi0=inc*dphiref
c: Amplitude of plane wave downward is 1, upward is R1x [both norm by 
c: 1/sqrt(W')]
         Ap=mnorm(ii)
         Apexp=0.d0
c: Do ref depth layer if we are doing bottom half and ref depth
c: at top of layer or if we are doing top half and ref depth is at bottom:
         if(isvmin .eq. ii) then
            j=jflu1
            jx1=jzmx(j)
c: FIX 2-6-95: Give em_calc phi,dphi at bottom of current layer:
            call phi_prop(j,ii1,ii2,inc,phi0,dphi0,Ap,Apexp)
            call em_calc(xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .etasq(j),gami(1,ii1,j),isp(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .philay(1,j),dphilay(1,j),Aplay(1,j),phix(jx1),dphix(jx1),
     .expx_gbs(jx1),xi,ai,aip,bi,bip,zzexp,inc,ailay(1,1,1,j),
     .bilay(1,1,1,j),zetalay(1,1,j),aisoln(1,j),ii1,ii2,j,jjhsp,kw0)
         endif
c
         do j=jflu1+inc,jflu2,inc
c: Propagate phi across fluid-fluid interface:
            phi0=rhorat(j+joff)*phi0
            jx1=jzmx(j)
            if(j .eq. jjhsp) then
               philay(ii1,j)=phi0
               dphilay(ii1,j)=inc*dphi0
               Aplay(ii1,j)=Apexp
            else
c: FIX 2-6-95: Give em_calc phi,dphi at bottom of current layer:
c: Compute phi,dphi at bottom of fluid layer:
               call phi_prop(j,ii1,ii2,inc,phi0,dphi0,Ap,Apexp)
            endif
            call em_calc(xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .etasq(j),gami(1,ii1,j),isp(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .philay(1,j),dphilay(1,j),Aplay(1,j),phix(jx1),dphix(jx1),
     .expx_gbs(jx1),xi,ai,aip,bi,bip,zzexp,inc,ailay(1,1,1,j),
     .bilay(1,1,1,j),zetalay(1,1,j),aisoln(1,j),ii1,ii2,j,jjhsp,kw0)
         enddo
c
         As=dcmplx(0.d0,0.d0)
         Asexp=-1.d+200
         if(inc*(jsol2-jsol1) .ge. 0) then
            j=jsol1
c: Compute p-wave and s-wave potentials at top of first solid layer:
            call phipsi(Ap,As,Vmat(1,1,j,ii),Wmat(1,j,ii),phi2,dphi2,
     .         Apexp,psi2,dpsi2,Asexp,j,ii1,inc)
         endif
         do j=jsol1,jsol2,inc
            if(j .ne. jjhsp) then
c: For layers "above" halfspace:
               j1=j + inc
c: Compute p-wave and s-wave potentials at top of next layer:
               call phipsi(Ap,As,Vmat(1,1,j1,ii),Wmat(1,j1,ii),
     .            phi0,dphi0,Apexp,psi0,dpsi0,Asexp,j1,ii1,inc)
c: Save these quantities at top of next layer for when we do halfspace:
               phi2=phi0
               dphi2=dphi0
               psi2=psi0
               dpsi2=dpsi0
               k=j1 + joff
c: Convert p- and s-wave potentials to bottom of current layer:
               call phi_xlay(phi0,dphi0,Apexp,psi0,dpsi0,Asexp,
     .            Pcon(1,k),Qcon(1,k),Ucon(1,k),Vcon(1,k),ikcon,
     .            rhorat(k),j,ii2,inc,mm(k))
            else
c: When doing halfspace, set phi0,dphi0,psi0,dpsi0 to values at top
c: of halfspace (note that Apexp will be unchanged and correct):
               phi0=phi2
               dphi0=dphi2
               psi0=psi2
               dpsi0=dpsi2
            endif
            jx1=jzmx(j)
c: Find p-wave potential in layer:
            call em_calc(xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .etasq(j),gami(1,ii1,j),isp(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .philay(1,j),dphilay(1,j),Aplay(1,j),phix(jx1),dphix(jx1),
     .expx_gbs(jx1),xi,ai,aip,bi,bip,zzexp,inc,ailay(1,1,1,j),
     .bilay(1,1,1,j),zetalay(1,1,j),aisoln(1,j),ii1,ii2,j,jjhsp,kw0)
c: Find s-wave potential in layer:
            call em_calc(xkhsq,xbsq(ii1,j),xbsq(ii2,j),inc*etb(j),
     .etbsq(j),beti(1,ii1,j),iss(j),zmx(jx1),zmx_im_gbs(jx1),nzmx(j),
     .psilay(1,j),dpsilay(1,j),Aslay(1,j),psix(jx1),dpsix(jx1),
     .expx_gbs(jx1),xis,ais,aips,bis,bips,zzexps,inc,ailay(1,1,2,j),
     .bilay(1,1,2,j),zetalay(1,2,j),aisoln(2,j),ii1,ii2,j,jjhsp,kw0)
         enddo
         if(ii .eq. 2) sg=-dsign(1.d0,dreal(dphi0))
      enddo
c
c: Convert from phix to phi using mx_m pointers computed in zsr_init:
      do jsr=1,nzsr
         jx=mx_m(jsr)
c: INCLUDE RHO in PHI,DPHI,PSI,DPSI:
c:pln sometimes phix is huge (no call to em_calc for this run)
c:pln Engineering solution
         if(cdabs(phix(jx)).gt.1.d10) then
            jjfail=1
            return
         end if
         phiz(jsr,jm)=sg*rho_sr(jsr)*phix(jx)
         dphiz(jsr,jm)=sg*rho_sr(jsr)*dphix(jx)
         psiz(jsr,jm)=sg*rho_sr(jsr)*psix(jx)
         dpsiz(jsr,jm)=sg*rho_sr(jsr)*dpsix(jx)
c: Exponential for complex source points:
         expz_gbs(jsr,jm)=expx_gbs(jx)
      enddo
cc    print *,'jm,phi_src,phi_rec = ',jm,cabs(phiz(mzsrc(1),jm)),
cc   .   cabs(phiz(mzsrc(1)-1,jm)),cabs(phiz(mzsrc(1)+1,jm))
c
      return
      end
ccc
      subroutine phi_prop(j,ii1,ii2,inc,phi0,dphi0,Ap,Apexp)
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 j,ii1,ii2,inc
      complex*16 Vmatx,Wmatx,Wexp,gami2,phi0,dphi0,Ap
      real*8 Apexp,Wphase
c
      philay(ii1,j)=phi0
      dphilay(ii1,j)=inc*dphi0
      Aplay(ii1,j)=Apexp
c: Compute phi0,dphi0 [to be multiplied by exp(Apexp)] at bottom of layer
c: from phi0,dphi0 at top of layer:
      Wexp=Wmat(5,j,ii1)
      Apexp=Apexp + dreal(Wexp)
      Wphase=dimag(Wexp)
c: Multiply plane wave amplitude by transmission coefficient [except for
c: the exp(Apexp) part]:
      Wmatx=Wmat(1,j,ii1)
      Ap=Ap*Wmatx*dcmplx(dcos(Wphase),dsin(Wphase))
      Vmatx=Vmat(1,1,j,ii1)
      gami2=gami(1,ii2,j)
      phi0=Ap*(1.d0 + Vmatx)
      dphi0=Ap*gami2*(1.d0 - Vmatx)
      philay(ii2,j)=phi0
      dphilay(ii2,j)=inc*dphi0
      Aplay(ii2,j)=Apexp
c
      return
      end
ccc
      subroutine phipsi(Ap,As,Vmatx,Wmatx,phi0,dphi0,Apexp,psi0,dpsi0,
     .   Asexp,j,ii1,inc)
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 j,ii1,inc
      complex*16 Ap,As,Vmatx(3,5),Wmatx(6),phi0,dphi0,psi0,
     .   dpsi0,vpp,Ap_vps,vss,As_vsp,cfacp,cfacs,cfac,
     .   term1p,term2p,term1s,term2s
      real*8 Apexp,Asexp,delex,spreal,Wphase,Vphase,ex1,ex2
c
c: Do imaginary part of cdexp(Wmatx(5)) term for p-wave transmission coeffs:
      Wphase=dimag(Wmatx(5))
      cfacp=Ap*dcmplx(dcos(Wphase),dsin(Wphase))
c: Do imaginary part of cdexp(Wmatx(6)) term for s-wave transmission coeffs:
      Wphase=dimag(Wmatx(6))
      cfacs=As*dcmplx(dcos(Wphase),dsin(Wphase))
c: Ap is to be multiplied by exp(Apexp), As by exp(Asexp):
cxx   Apexp=Apexp + dreal(Wmatx(5))
cxx   Asexp=Asexp + dreal(Wmatx(6))
cxx   delex=Apexp - Asexp
cxx   efac=dexp(delex)
c: Compute new amplitudes of downward-going p-wave and s-wave plane waves:
      term1p=Wmatx(1)*cfacp
      term2p=Wmatx(4)*cfacs
      ex1=Apexp + dreal(Wmatx(5))
      term1s=Wmatx(3)*cfacs
      term2s=Wmatx(2)*cfacp
      ex2=Asexp + dreal(Wmatx(6))
c: Make sure that Ap,As have magnitudes on the order of unity:
      call norm_add(term1p,ex1,term2p,ex2,Ap,Apexp)
      call norm_add(term1s,ex2,term2s,ex1,As,Asexp)
cxx   Ap=Apx*Wmatx(1)*cfacp + Asx*Wmatx(4)*cfacs/efac
cxx   As=Asx*Wmatx(3)*cfacs + Apx*Wmatx(2)*cfacp*efac
c: Compute p-wave and s-wave potentials at top of layer:
      Vphase=dimag(Vmatx(1,5))
      cfac=dcmplx(dcos(Vphase),dsin(Vphase))
      spreal=dreal(Vmatx(1,5))
      delex=Apexp - Asexp
      vpp=Vmatx(1,1)
c: Include exponential factors due to Vps and Ap,As:
      Ap_vps=Ap*Vmatx(1,2)*cfac*dexp(spreal + delex)
      vss=Vmatx(1,3)
c: Include exponential factors due to Vsp and Ap,As:
      As_vsp=As*Vmatx(1,4)*cfac*dexp(spreal - delex)
c: phi0,dphi0 must be multiplied by exp(Apexp):
      phi0=Ap*(1.d0 + vpp) + As_vsp
      dphi0=gami(1,ii1,j)*(Ap*(1.d0 - vpp) - As_vsp)
c: psi0,dpsi0 must be multiplied by exp(Asexp):
      psi0=As*(1.d0 + vss) + Ap_vps
      dpsi0=beti(1,ii1,j)*(As*(1.d0 - vss) - Ap_vps)
c
      philay(ii1,j)=phi0
      dphilay(ii1,j)=inc*dphi0
      Aplay(ii1,j)=Apexp
      psilay(ii1,j)=psi0
      dpsilay(ii1,j)=inc*dpsi0
      Aslay(ii1,j)=Asexp
c
      return
      end
ccc
      subroutine norm_add(term1,ex1,term2,ex2,A,Aexp)
c
c: Computes A*exp(Aexp) = term1*exp(ex1) + term2*exp(ex2), such
c: that the magnitude of A is on the order of unity.
c
      implicit none
      complex*16 term1,term2,A
      real*8 ex1,ex2,Aexp,mag1,mag2,exx1,exx2
c
      if(term1 .ne. (0.d0,0.d0)) then
         mag1=cdabs(term1)
         term1=term1/mag1
         exx1=ex1 + dlog(mag1)
      else
         exx1=-1.d200
      endif
      if(term2 .ne. (0.d0,0.d0)) then
         mag2=cdabs(term2)
         term2=term2/mag2
         exx2=ex2 + dlog(mag2)
      else
         exx2=-1.d200
      endif
      if(exx1 .gt. exx2) then
c: term1 dominates:
         A=term1 + term2*dexp(exx2 - exx1)
         Aexp=exx1
      else
c: term2 dominates:
         A=term1*dexp(exx1 - exx2) + term2
         Aexp=exx2
      endif
c
      return
      end
ccc
      subroutine phi_xlay(phi0,dphi0,Apexp,psi0,dpsi0,Asexp,Pconx,
     .   Qconx,Uconx,Vconx,ikconx,rhoratx,j,ii2,inc,mmx)
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 j,ii2,inc,mmx
      complex*16 phi0,dphi0,psi0,dpsi0,Pconx,Qconx,Uconx,Vconx,ikconx,
     .   phib,dphib,psib,dpsib,Pbar,Qbar,Ubar,Vbar
      real*8 Apexp,Asexp,rhoratx,efac
c
cc    print *,'uz below: ',dphi0*dexp(Apexp) + ikconx*psi0*dexp(Asexp)
c: Find P,Q,U,V for going backward across interface (see p. 121):
c: If no mismatch, phi and psi are continuous:
      if(mmx .eq. 1) then
         Pbar=Qconx/rhoratx
         Qbar=Pconx/rhoratx
         Ubar=Uconx/(-rhoratx)
         Vbar=ikconx*(1.d0 - Pbar)
c
c: Propagate Phi backward across solid-solid interface (see p. 87):
         efac=exp(Apexp - Asexp)
         phib=Pbar*phi0 + Ubar*dpsi0/efac
         dphib=Qbar*dphi0 + Vbar*psi0/efac
         psib=-Ubar*dphi0*efac + Pbar*psi0
         dpsib=-Vbar*phi0*efac + Qbar*dpsi0
c
         phi0=phib
         dphi0=dphib
         psi0=psib
         dpsi0=dpsib
      endif
cc    print *,'uz above: ',dphi0*dexp(Apexp) + ikconx*psi0*dexp(Asexp)
c
      philay(ii2,j)=phi0
      dphilay(ii2,j)=inc*dphi0
      Aplay(ii2,j)=Apexp
      psilay(ii2,j)=psi0
      dpsilay(ii2,j)=inc*dpsi0
      Aslay(ii2,j)=Asexp
c
      return
      end
      subroutine mode_ortho(phi_,dphi_,psi_,dpsi_)
c
c: Checks for mode orthogonality and normality.  See ORCA II pp. 4-5.
c: Note: In order to get normality, set iimf=1 and sample modes very 
c: finely from just inside the upper halfspace to just inside the lower
c: halfspace.  Orthogonality holds over any subinterval of the waveguide
c: as long as the modes are well sampled over the region.  See ORCA II,p36.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jm,jm2,jzmf,jsr,j,j1,ji
      complex*8 phi_(nzsr,nmode),dphi_(nzsr,nmode),
     .   psi_(nzsr,nmode),dpsi_(nzsr,nmode),Qn,Pn,kn_ks,vfac,dvfac
      complex*16 igral,hsp_lower,hsp_upper,dksq,int_terms,above,below,
     .   min_ik
c
      if(zmf(1) .ge. zdep(1) .or. zmf(nzmf) .le. zdep(nlay-1)) then
         print *,'WARNING: MODE NORMALITY WILL ONLY HOLD IF FIRST ' 
         print *,'   AND LAST MODE FUNCTION DEPTHS ARE IN UPPER AND '
         print *,'   LOWER HALFSPACES, RESPECTIVELY.'
      endif
      wsq=w*w
      do jm=1,nmode
         igral=(0.d0,0.d0)
         do jzmf=1,nzmf
            jsr=mzmf(jzmf)
c: NOTE 8/31/94: Conjugation does not work (Frisk wrong)!
cxx         igral=igral+conjg(phi_(jsr,jm))*phi_(jsr,jm)
            if(kksh(jsr) .eq. 0) then
               igral=igral + phi_(jsr,jm)*phi_(jsr,jm)/rho_sr(jsr)
            else
c: Try using v=-psi/(i*k) instead of v=psi:
               min_ik=dcmplx(dimag(kn(jm)),-dreal(kn(jm)))
               vfac=psi_(jsr,jm)/min_ik
               dvfac=dpsi_(jsr,jm)/min_ik
               kn_ks=2.*kn(jm)*kn(jm) - 1./ksm2_sr(jsr)
               Qn=phi_(jsr,jm) + dvfac
               Pn=2.*dphi_(jsr,jm) + kn_ks*vfac
               igral=igral + (Qn*Qn + Pn*vfac)/rho_sr(jsr)
            endif
         enddo
         igral=igral*(zmf(2)-zmf(1))
c: Use actual upper and lower limits from phi_,dphi_:
         jsr=mzmf(nzmf)
         call hsp_terms(phi_(jsr,jm),dphi_(jsr,jm),psi_(jsr,jm),
     .      kn(jm),wsq,cp_sr(jsr),cs_sr(jsr),rho_sr(jsr),kksh(jsr),
     .      hsp_lower)
         jsr=mzmf(1)
         call hsp_terms(phi_(jsr,jm),dphi_(jsr,jm),psi_(jsr,jm),
     .      kn(jm),wsq,cp_sr(jsr),cs_sr(jsr),rho_sr(jsr),kksh(jsr),
     .      hsp_upper)
         igral=igral - (hsp_lower - hsp_upper)
         print *,'jm=',jm,sngl(real(igral)),sngl(dimag(igral)),
     .      sngl(abs(hsp_lower)),sngl(abs(hsp_upper))
      enddo
c
      do jm=1,nmode
         do jm2=jm+1,nmode
            igral=(0.,0.)
            do jzmf=1,nzmf
               jsr=mzmf(jzmf)
c: NOTE 8/31/94: Conjugation does not work!
cxx            igral=igral+conjg(phi_(jsr,jm))*phi_(jsr,jm2)
               igral=igral+phi_(jsr,jm)*phi_(jsr,jm2)/rho_sr(jsr)
            enddo
            dksq=kn(jm)*kn(jm) - kn(jm2)*kn(jm2)
            igral=dksq*igral*(zmf(2)-zmf(1))
c: Account for fluid-solid and solid-solid interfaces:
            int_terms=cmplx(0.d0,0.d0)
            do ji=1,n_int
               j=mzint(ji)
               j1=j+1
               above=(phi_(j,jm2)*dphi_(j,jm)
     .             - phi_(j,jm)*dphi_(j,jm2))/rho_sr(j)
               below=(phi_(j1,jm2)*dphi_(j1,jm)
     .             - phi_(j1,jm)*dphi_(j1,jm2))/rho_sr(j1)
               int_terms=int_terms + (above - below)
            enddo
            jsr=mzmf(nzmf)
            hsp_lower=(phi_(jsr,jm2)*dphi_(jsr,jm) - 
     .         phi_(jsr,jm)*dphi_(jsr,jm2))/rho_sr(jsr)
            jsr=mzmf(1)
            hsp_upper=(phi_(jsr,jm2)*dphi_(jsr,jm) - 
     .         phi_(jsr,jm)*dphi_(jsr,jm2))/rho_sr(jsr)
            igral=igral - int_terms - (hsp_lower - hsp_upper)
         print *,'jm,jm2=',jm,jm2,sngl(abs(igral)),
     .      sngl(abs(hsp_lower)),sngl(abs(hsp_upper)),
     .      sngl(abs(int_terms))
         enddo
      enddo
c      
      return
      end
ccc
      subroutine hsp_terms(phi,dphi,psi,kn,wsq,cp,cs,rho,kksh,
     .   hsp_term)
c
c: Computes halfspace terms of mode normalization.
c
      implicit none
      complex*16 kn,cp,cs,knsq,gami,beti,ikn,hsp_term
      complex*8 phi,dphi,psi
      real*8 wsq,rho
      integer*4 kksh
c
      if(kksh .eq. 0) then
         if(dphi .ne. (0.,0.)) then
            hsp_term=phi*phi*phi/(2.*dphi*rho)
         else
            hsp_term=(0.d0,0.d0)
         endif
      else
         knsq=kn*kn
         gami=cdsqrt(knsq - wsq/(cp*cp))
         beti=cdsqrt(knsq - wsq/(cs*cs))
         ikn=dcmplx(-dimag(kn),dreal(kn))
         hsp_term=(phi*phi/(2.*gami) + 2.*phi*psi/ikn + psi*psi*
     .      (2.*beti*beti + knsq)/(2.*ikn*ikn*beti))/rho
      endif
c
      return
      end
ccc
      subroutine mode_ortho_analytic
c
c: Checks for mode orthogonality and normality using analytic formulas.  
c: See Levinson et al, Eqs. (28) and (A6).
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 j
      complex*16 igral,xx
c
      igral=(0.d0,0.d0)
      do j=2,nlay-1
         call mode_ortho_lay(philay(1,j),dphilay(1,j),Aplay(1,j),
     .      xkhsq,xksq(1,j),geo(1,3,j),eta(j),etasq(j),h(j),isp(j),xx)
         igral=igral + xx
      enddo
      j=nlay
      philay(2,j)=0.d0
      dphilay(2,j)=0.d0
      call mode_ortho_lay(philay(1,j),dphilay(1,j),Aplay(1,j),
     .   xkhsq,xksq(1,j),geo(1,3,j),eta(j),etasq(j),h(j),isp(j),xx)
      igral=igral + xx
c
      print *,'nmode = ',nmode,kduct,dreal(xkh)/kw0,dimag(xkh)*8685.9,
     .   cdabs(igral)
c
      return
      end
ccc
      subroutine mode_ortho_lay(philay,dphilay,Aplay,xkhsq,xksq,rho,
     .   eta,etasq,h,isp,xx)
c
      implicit none
      complex*16 philay(2),dphilay(2),xkhsq,xksq(2),eta,etasq,xx,
     .   phi1sq,dphi1sq,phi2sq,dphi2sq,mgam1sq,mgam2sq,
     .   term1,term2,pdp1,pdp2
      real*8 Aplay(2),rho(2),h,rhoexp1,rhoexp2
      integer*4 isp
c
      phi1sq=philay(1)*philay(1)
      dphi1sq=dphilay(1)*dphilay(1)
      phi2sq=philay(2)*philay(2)
      dphi2sq=dphilay(2)*dphilay(2)
      mgam1sq=xkhsq - xksq(1)
      mgam2sq=xkhsq - xksq(2)
      rhoexp1=rho(1)*dexp(2.d0*Aplay(1))
      rhoexp2=rho(2)*dexp(2.d0*Aplay(2))
      if(isp .eq. 0) then
         term1=rhoexp1*(dphi1sq - mgam1sq*phi1sq)
         term2=rhoexp2*(dphi2sq - mgam2sq*phi2sq)
         xx=(term2 - term1)/(eta*etasq)
      else
         pdp1=philay(1)*dphilay(1)
         pdp2=philay(2)*dphilay(2)
         term1=rhoexp1*pdp1/mgam1sq
         term2=rhoexp2*(phi2sq*h - (dphi2sq*h - pdp2)/mgam2sq)
         xx=0.5d0*(term2 - term1)
      endif
c
      return
      end
c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************
                                                 
      subroutine opt_read_pro
c
c: Reads option file.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
c
      integer*4 nline,iierr
      integer*4 IOP

      character*64 eline

c: Local variables:
      integer*4 jj

      data eline/'INVALID INPUT IN OPT FILE: '/
c
      iierr=0
      nline=0
      svp_ver= 2.0
      rho_svp= 1.0
      alpha_svp= 0.0
c
      nline=0
      iierr=0
      iigeom=1
      ver_no= 2.0
      IOP= 16
      iikpl= 0
      iirc= 0
      n_env= 0
c      iirx= -1
      cphmin= 0.0
      cphmax= 0.0
      rmin=0.0
      rmax= 0.0
      phfac= 4
      db_cut= 48
      iidiag= 0
      lout=11
      iiaih(1)=0.0
      iiaih(2)=0.0
      iigbs=0
c
cfmc  iifft   = Output FFT File (0=no;1=zs,zr,r on Line 9; 2=read file on Line 11);
      iifft= 1
cfmc  iiout   = Output BB Eigenvalues and Functions (same options as iifft above);
      iiout= 0
cfmc   iift    = Freq Traj(ASCII); iimt = Mode Traj(ASCII); 
      iift= 0
cfmc    iimt= ????
      iimt= 0
cfmc  iidc    = Disp Curves (0=no,1=vg,2=vph,3=both); 
      iidc= 0

      call check_val_r4(fsbb,0.e0,1.e10,eline,27,nline,
     .         'fsbb',4,iierr)
      call check_val_r4(Tw,-1.e10,131072e0,eline,27,nline,
     .         'nfft/Tw',7,iierr)
      call check_val_r4(fmin,1.e-3,fmax,eline,27,nline,
     .         'fmin',4,iierr)
      call check_val_r4(fmax,fmin,fsbb/2.e0,eline,27,nline,
     .         'fmax',4,iierr)
c
      if(iierr .eq. 1) then
         print *,' '
         print *,'Execution terminating.  Check input option file '//
     .      'for error(s).'
         stop
      endif
      return
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
      subroutine out_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
     .   xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt,ncall)
c
c: This subroutine chooses the subroutine to call in order to output 
c: the array cmat and the axes x and y.
c
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
      integer*4 ncall,lout0,lsuf,nx,ny,nf,ncasc,mat_fid,iilong
      character*64 outroot0,ascsuf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt
      character*30 xlab,ylab
      real*4 cmat(ny,nx),x(nx),y(ny),xmin,xmax,ymin,ymax
      data iilong/0/
c
      if(iiwrite .eq. 0) then
         return
      endif
c: ncall<0 means to close Matlab file:
      if(ncall .lt. 0) then
cpln         stat = matClose(mat_fid)
         if(iilong .gt. 0) then
            print *,'Warning: Variable names in .mat file '//
     .         'shortened'
            write(nf,*) 'Warning: Variable names in .mat file '//
     .         'shortened'
         endif
         return
      endif
c
      if(iifmt .eq. 1 .or. iifmt .eq. 0) then
cpln         call hdf_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
cpln     .      xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
cpln     .      xfmt,yfmt,dfmt,ncall)
      endif
      if(iifmt .eq. 2 .or. iifmt .eq. 0) then
cpln         call mat_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
cpln     .      xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
cpln     .      xfmt,yfmt,dfmt,ncmat,mat_fid,iilong)
      endif
      if(iifmt .eq. 3 .or. iifmt .eq. 0) then
         call asc_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
     .      xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
     .      xfmt,yfmt,dfmt,ncasc)
      endif
c
      return
      end
ccc
      subroutine asc_writex(outroot,lout,ascsuf,lsuf,cmat,x,y,nx,ny,
     .   xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt,ncall)
c
c: This subroutine outputs the array cmat and the axes x and y
c: to an ASCII file named outroot(1:lout)//ascsuf(1:lsuf)//.asc
c: If xmin not equal to xmax, x(1:nx) is filled from xmin to xmax.
c: Likewise with ymin,ymax and y(1:ny).  If nf is not equal to 
c: zero, a message is written to that file number.
c
c: This subroutine writes out an ascii file of the data in cmat.
c
      implicit none
      include 'Parms_com'
      integer*4 ncall,lout,lsuf,nx,ny,nf,j,jx,ldat
      character*64 outroot,ascsuf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt,dataname
      character*30 xlab,ylab
      real*4 cmat(ny,nx),x(nx),y(ny),xmin,xmax,ymin,ymax,delx,dely
      real*8 rc_code
c
      dataname=outroot(1:lout)//'_asc'//ascsuf(1:lsuf)
      ldat=lout + 4 + lsuf
c
      if(nx*ny .le. 0) then
         write(nf,110) dataname(1:ldat),nx,ny
110      format('OUTPUT ASCII FILE = ',a,'; # ROWS =',i6,'; # COLS =',
     .      i6/'  NOT WRITTEN OUT SINCE ZERO SIZE')
         return
      endif
c
      if(xmin .ne. xmax .and. xmin .ne. -999.) then
         delx=(xmax-xmin)/float(max0(1,nx-1))
         do j=1,nx
            x(j)=xmin + (j-1)*delx
         enddo
      endif
      if(ymin .ne. ymax .and. ymin .ne. -999.) then
         dely=(ymax-ymin)/float(max0(1,ny-1))
         do j=1,ny
            y(j)=ymin + (j-1)*dely
         enddo
      endif
c
      open(83,file=dataname(1:ldat),status='unknown',form='formatted')
c: Include number of rows (ny) and number of columns (nx) in dummy word in
c: corner.  Decode as ny=int(rc_code); nx=10000*(rc_code-ny).
      rc_code=ny + nx/10000.d0
      write(83,200) rc_code,(y(j),j=1,ny)
200   format(2000(e14.8,1x))
      do jx=1,nx
         write(83,200) x(jx),(cmat(j,jx),j=1,ny)
      enddo
      close(83)
c
      if(nf .ne. 0) then
         write(nf,100) dataname(1:ldat),nx,ny
100      format('ASCII FILE = ',a,'; # ROWS =',i6,'; # COLS =',i6)
      endif
      ncall=ncall + 1
c
      return
      end
ccc
      subroutine phz_calc
c
c: Computes number of zeros in mode functions.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 j,iiev_top,iiev_bot
      real*4 phz_top,phz_bot,phz_lay(NLMAX),phz_phi,
     .   twpiex,piex,phz_bl,phz_psi,phz_lays
      data twpiex/6.283185307/,piex/3.14159265358979/
c
      phz_phi=0.
      phz_psi=0.
      phz_bl=0.
      mode_no(1,nmode)=0
      mode_no(2,nmode)=0
      do j=2,nlay-1
         call phi_phase(philay(1,j),dphilay(1,j),gami(1,1,j),
     .      xkh,xk(1,j),phz_top,iiev_top)
         call phi_phase(philay(2,j),dphilay(2,j),gami(1,2,j),
     .      xkh,xk(2,j),phz_bot,iiev_bot)
         call phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .      isp(j),h(j),gami(1,1,j),zetalay(1,1,j),zetalay(2,1,j),
     .      piex,twpiex,phz_lay(j))
c
         phz_phi=phz_phi + phz_lay(j)
         phz_bl=phz_bot
         if(iisol(j) .eq. 1) then
            call phi_phase(psilay(1,j),dpsilay(1,j),beti(1,1,j),
     .         xkh,xb(1,j),phz_top,iiev_top)
            call phi_phase(psilay(2,j),dpsilay(2,j),beti(1,2,j),
     .         xkh,xb(2,j),phz_bot,iiev_bot)
            call phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .         iss(j),h(j),beti(1,1,j),zetalay(1,2,j),zetalay(2,2,j),
     .         piex,twpiex,phz_lays)
            phz_psi=phz_psi + phz_lays
         endif
      enddo
c
      mode_phz(1,nmode)=phz_phi/piex
      mode_phz(2,nmode)=phz_psi/piex
      mode_phz(3,nmode)=mode_phz(1,nmode) + mode_phz(2,nmode)
      print *,'nmode = ',nmode,(mode_phz(j,nmode),j=1,3)
c
      return
      end
ccc
      subroutine phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .   iso,h,gami,zeta_top,zeta_bot,piex,twpie,phz_lay)
c
c: Computes phase change in a layer given phase at top and bottom and
c: iso or airy parameters.
c
      implicit none
      integer*4 iso,iiev_top,iiev_bot
      real*4 phz_top,phz_bot,phz_mod,phz_dif,dphz,twpie,piex,phz_lay
      real*8 h
      complex*16 gami,zeta_top,zeta_bot
c
c: Force layer phase difference to be zero when evanescent:
      if(iiev_top .eq. 1 .and. iiev_bot .eq. 1) then
         phz_lay=0.
         return
      endif
      phz_dif=phz_bot - phz_top
      if(phz_dif .lt. 0.) phz_dif=phz_dif + twpie
      if(iso .eq. 1) then
         phz_lay=dabs(h*dimag(gami))
      else
         phz_lay=abs(dimag(zeta_top)-dimag(zeta_bot))
      endif
      phz_mod=mod(phz_lay,twpie)
      dphz=phz_dif - phz_mod
      if(abs(dphz) .gt. piex) dphz=dphz - sign(twpie,dphz)
      phz_lay=phz_lay + dphz
c
      return
      end
ccc
      subroutine nzero_lay(phz_lay,piex,phi_t,phi_b,dphi_t,dphi_b,
     .   iiev_t,iiev_b,nz)
c
      implicit none
      integer*4 nz,iiev_t,iiev_b,iwave_t,iwave_b
      real*4 phz_lay,piex
      complex*8 phi_t,phi_b,dphi_t,dphi_b
c
      if(iiev_t .eq. 1 .and. iiev_b .eq. 1) then
         nz=0
         iwave_t=0
         if((real(phi_t) .gt. 0. .and. real(dphi_t) .lt. 0.) .or.
     .      (real(phi_t) .lt. 0. .and. real(dphi_t) .gt. 0.)) iwave_t=1
         iwave_b=0
         if((real(phi_b) .gt. 0. .and. real(dphi_b) .gt. 0.) .or.
     .      (real(phi_b) .lt. 0. .and. real(dphi_b) .lt. 0.)) iwave_b=1
c: If interface wave at top and bottom, keep nz=0 and return:
         if(iwave_t .eq. 1 .and. iwave_b .eq. 1) return
      else
         nz=int(phz_lay/piex)
      endif
      if(mod(nz,2) .eq. 0) then
         if((real(phi_t) .gt. 0. .and. real(phi_b) .le. 0.) .or.
     .      (real(phi_t) .lt. 0. .and. real(phi_b) .ge. 0.)) nz=nz+1
      else
         if((real(phi_t) .gt. 0. .and. real(phi_b) .ge. 0.) .or.
     .      (real(phi_t) .lt. 0. .and. real(phi_b) .le. 0.)) nz=nz+1
      endif
c
      return
      end
ccc
      subroutine phz_calc_old
c
c: Computes total phase of mode function.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 j,iiev(2,NLMAX),iievs
      real*4 phz_top,phz_bot,phz_lay(NLMAX),phz_mod,
     .   phz_phi,phz_dif,dphz,twpiex,piex,phz_bl,
     .   phz_duct,phz_psi,phz_lays
      data twpiex/6.28318530717959/,piex/3.14159265358979/
c
      phz_phi=0.
      phz_psi=0.
      phz_bl=0.
      do j=2,nlay-1
         call phi_phase(philay(1,j),dphilay(1,j),gami(1,1,j),
     .      xkh,xk(1,j),phz_top,iiev(1,j))
         call phi_phase(philay(2,j),dphilay(2,j),gami(1,2,j),
     .      xkh,xk(2,j),phz_bot,iiev(2,j))
         phz_dif=phz_bot - phz_top
c: Force layer phase difference to be zero when evanescent:
         if(iiev(1,j) .eq. 1 .and. iiev(2,j) .eq. 1) phz_dif=0
         if(phz_dif .lt. 0.) phz_dif=phz_dif + twpiex
         if(isp(j) .eq. 1) then
            phz_lay(j)=dabs(h(j)*dimag(gami(1,1,j)))
         else
            phz_lay(j)=abs(dimag(zetalay(2,1,j))-dimag(zetalay(1,1,j)))
         endif
         phz_mod=mod(phz_lay(j),twpiex)
         dphz=phz_dif - phz_mod
         if(abs(dphz) .gt. piex) dphz=dphz - sign(twpiex,dphz)
         phz_lay(j)=phz_lay(j) + dphz
         phz_phi=phz_phi + phz_lay(j)
      print *,nmode,j,iiev(1,j),iiev(2,j),phz_lay(j)/piex,phz_top/piex,
     .   phz_bot/piex,(phz_top-phz_bl)/piex
         phz_bl=phz_bot
         if(iisol(j) .eq. 1) then
            call phi_phase(psilay(1,j),dpsilay(1,j),beti(1,1,j),
     .         xkh,xb(1,j),phz_top,iievs)
            call phi_phase(psilay(2,j),dpsilay(2,j),beti(1,2,j),
     .         xkh,xb(2,j),phz_bot,iievs)
            phz_dif=phz_bot - phz_top
            if(phz_dif .lt. 0.) phz_dif=phz_dif + twpiex
            if(iss(j) .eq. 1) then
               phz_lays=dabs(h(j)*dimag(beti(1,1,j)))
            else
               phz_lays=abs(dimag(zetalay(2,2,j))-dimag(zetalay(1,2,j)))
            endif
            phz_mod=mod(phz_lay(j),twpiex)
            dphz=phz_dif - phz_mod
            if(abs(dphz) .gt. piex) dphz=dphz - sign(twpiex,dphz)
            phz_lays=phz_lays + dphz
            phz_psi=phz_psi + phz_lays
         endif
      enddo
c
c: Compute phase only in current duct (until mode goes evanescent):
      phz_duct=0.
      do j=nsvmin,2,-1
         phz_duct=phz_duct + phz_lay(j)
         if(iiev(1,j) .eq. 1 .or. iiev(2,j) .eq. 1) goto 10
      enddo
10    continue
      do j=nsvmin+1,nlay-1
         phz_duct=phz_duct + phz_lay(j)
         if(iiev(1,j) .eq. 1 .or. iiev(2,j) .eq. 1) goto 20
      enddo
20    continue
      mode_no(1,nmode)=nint(phz_phi/piex)
      mode_no(2,nmode)=nint(phz_psi/piex)
      mode_no(3,nmode)=mode_no(1,nmode) + mode_no(2,nmode)
      print *,'nmode = ',nmode,(mode_no(j,nmode),j=1,3),
     .   phz_phi/piex,phz_psi/piex,phz_duct/piex
c
      return
      end
ccc
      function deg(z)
c
      implicit none
      complex*8 z
      real*4 deg
c
      deg=atan2(aimag(z),real(z))*180./3.14159265
      return
      end
ccc
      function mag(z)
c
      implicit none
      complex*8 z
      real*4 mag
c
      mag=sqrt(real(z)*real(z) + aimag(z)*aimag(z))
      return
      end
ccc
      subroutine phi_phase(phi,dphi,gami,xkh,xk,phz,iiev)
c
      implicit none
      integer*4 iiev
      complex*16 gami,xkh,xk,phi,dphi
      real*4 phz,gamma,num_re,num_im,mag_re,mag_im
c
      gamma=dimag(gami)
      num_re=dreal(phi)*gamma
      num_im=dimag(phi)*gamma
      mag_re=num_re*num_re+dreal(dphi)*dreal(dphi)
      mag_im=num_im*num_im+dimag(dphi)*dimag(dphi)
      if(mag_re .gt. mag_im) then
         phz=atan2(num_re,sngl(dreal(dphi)))
      else
         phz=atan2(num_im,sngl(dimag(dphi)))
      endif
      iiev=0
      if(dreal(xkh) .gt. dreal(xk)) iiev=1
c
      return
      end
ccc
      subroutine phi_phase_old(phi,dphi,gami,xkh,xk,phz,iiev)
c
      implicit none
      integer*4 iiev
      complex*16 gami,xkh,xk,phi,dphi
      real*4 phz,gamma,num
c
      gamma=dimag(gami)
      num=dreal(phi)*dimag(gami)
      phz=atan2(num,sngl(dreal(dphi)))
      iiev=0
      if(dreal(xkh) .gt. dreal(xk)) iiev=1
c
      return
      end
ccc
      subroutine phz_calc_next
c
c: Computes number of zeros in mode functions.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 j,iiev_top,iiev_bot,nz
      real*4 phz_top,phz_bot,phz_lay(NLMAX),phz_phi,
     .   twpiex,piex,phz_bl,phz_psi,phz_lays
      data twpiex/6.28318530717959/,piex/3.14159265358979/
c
      phz_phi=0.
      phz_psi=0.
      phz_bl=0.
      mode_no(1,nmode)=0
      mode_no(2,nmode)=0
      do j=2,nlay-1
         call phi_phase(philay(1,j),dphilay(1,j),gami(1,1,j),
     .      xkh,xk(1,j),phz_top,iiev_top)
         call phi_phase(philay(2,j),dphilay(2,j),gami(1,2,j),
     .      xkh,xk(2,j),phz_bot,iiev_bot)
         call phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .      isp(j),h(j),gami(1,1,j),zetalay(1,1,j),zetalay(2,1,j),
     .      piex,twpiex,phz_lay(j))
c
         call nzero_lay(phz_lay(j),piex,philay(1,j),philay(2,j),
     .      dphilay(1,j),dphilay(2,j),iiev_top,iiev_bot,nz)
         mode_no(1,nmode)=mode_no(1,nmode) + nz
         phz_bl=phz_bot
         if(iisol(j) .eq. 1) then
            call phi_phase(psilay(1,j),dpsilay(1,j),beti(1,1,j),
     .         xkh,xb(1,j),phz_top,iiev_top)
            call phi_phase(psilay(2,j),dpsilay(2,j),beti(1,2,j),
     .         xkh,xb(2,j),phz_bot,iiev_bot)
            call phase_lay(phz_top,phz_bot,iiev_top,iiev_bot,
     .         iss(j),h(j),beti(1,1,j),zetalay(1,2,j),zetalay(2,2,j),
     .         piex,twpiex,phz_lays)
            phz_psi=phz_psi + phz_lays
c
            call nzero_lay(phz_lays,piex,psilay(1,j),psilay(2,j),
     .         dpsilay(1,j),dpsilay(2,j),iiev_top,iiev_bot,nz)
            mode_no(2,nmode)=mode_no(2,nmode) + nz
         endif
      enddo
c
      mode_no(3,nmode)=mode_no(1,nmode) + mode_no(2,nmode)
      print *,'nmode = ',nmode,(mode_no(j,nmode),j=1,3)
c
      return
      end
      subroutine r1r2_calc(k,r1r2,ndv,iiw,jjfail)
c
c: Computes the product of the downlooking (rc1)  and uplooking (rc2)
c: reflection coefficients at the horizontal wavenumber k, and (if ndv>0) 
c: Also computes deriviates with respect to k (for ndv>=2) and with respect
c: to w (for ndv=3).
c: Output: r1r2(3,4), where 
c:    r1r2(1:3,1)=R1 and its two derivatives,
c:    r1r2(1:3,2)=R2 and its two derivatives,
c:    r1r2(1:3,3)=R1*R2 and its two derivatives,
c:    r1r2(1:3,4)=ln(R1*R2) and its two derivatives,
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer ndv,iiw,j,jjfail
      complex*16 k,r1r2(3,4)
c
      nctot=nctot + 1
      xkh=k
      phtot=0.
      iicut=0
      call xkh_init(ndv)
      call rp_calc(1,r1r2(1,1),ndv,iiw)
      if(jjfail .gt. 0) return
      call rp_calc(2,r1r2(1,2),ndv,iiw)
      if(jjfail .gt. 0) return
c      if(cdabs(r1r2(1,3)).lt.1.d-306) then
c         jjfail = 1
c         return
c      endif
c
      r1r2(1,3)=r1r2(1,1)*r1r2(1,2)
      r1r2(1,4)=cdlog(r1r2(1,3))
      do j=2,ndv
         r1r2(j,3)=r1r2(1,1)*r1r2(j,2) + r1r2(j,1)*r1r2(1,2)
         r1r2(j,4)=r1r2(j,3)/r1r2(1,3)
      enddo
      if(iicut .eq. 1) lncut=r1r2(1,4)
c
      return
      end
      subroutine rp_calc(ii,R,ndv,iiw)
c
c: Computes the plane wave reflection coefficient for a series of
c: fluid and solid layers.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c: Local variables:
      integer ii,ns,jflu1,jflu2,inc,j,ndv,iiw,k,ii1,ii2,
     .   jsol1,jsol2,iiwww,j1
      complex*16 e11(3),e12(3),e21(3),e22(3),s11(3),s12(3),s21(3),
     .   s22(3),zexpe,zexps,R(3)
c
      ii1=ii
      ii2=3 - ii
      inc=jflu(ii,3)
      jflu1=jflu(ii,1)
      jflu2=jflu(ii,2)
      jsol1=jsol(ii,1)
      jsol2=jsol(ii,2)
      ns=jhsp(ii)
c
c: Compute reflection coefficient at top of halfspace (0 for homogeneous):
      call h_space(ii,isp(ns),xkh,xkhsq,xksq(ii1,ns),eta(ns),etasq(ns),
     .   gami(1,ii1,ns),w,Vmat(1,1,ns,ii),ailay(1,ii,1,ns),
     .   zetalay(ii,1,ns),ndv,iiw,1,2,ihf(ns),xi_hsp(ii,1))
      if(allf(ii) .eq. 0) then
         call h_space(ii,iss(ns),xkh,xkhsq,xbsq(ii1,ns),etb(ns),
     .      etbsq(ns),beti(1,ii1,ns),w,Vmat(1,1,ns,ii),
     .      ailay(1,ii,2,ns),zetalay(ii,2,ns),ndv,iiw,3,4,
     .      ihf(ns),xi_hsp(ii,2))
      endif
c: Test for duct at top of Airy halfspace:
      if(jflu1 .eq. ns) then
         do j=1,ndv
            R(j)=Vmat(j,1,ns,ii)
         enddo
         return
      endif
c
c: Loop over solid layers:
      do j=jsol2,jsol1,-inc
c: j1 points to layer to be propagated from:
         j1=j + inc
c: k points to interface:
         k=j + jflu(ii,4)
         iiwww=iiw*iiww(j1)
c: Compute 2x2 propagator matrices for compressional and shear profiles:
         call ep_calc(xkh,w,xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .      etasq(j),gami(1,ii2,j),h(j),e11,e12,e21,e22,zexpe,isp(j),
     .      ihf(j),1,ndv,iiwww,ailay(1,1,1,j),bilay(1,1,1,j),
     .      zetalay(1,1,j),aisoln(1,j),ii1,ii2)
         call ep_calc(xkh,w,xkhsq,xbsq(ii1,j),xbsq(ii2,j),inc*etb(j),
     .      etbsq(j),beti(1,ii2,j),h(j),s11,s12,s21,s22,zexps,iss(j),
     .      ihf(j),1,ndv,iiwww,ailay(1,1,2,j),bilay(1,1,2,j),
     .      zetalay(1,2,j),aisoln(2,j),ii1,ii2)
         if(mm(k) .eq. 1) then
c: Mismatch at bottom of layer.  Call rp_slay as usual:
            call rp_slay(Vmat(1,1,j1,ii),Pcon(1,k),Qcon(1,k),
     .         Ucon(1,k),Vcon(1,k),e11,e12,e21,e22,zexpe,s11,s12,
     .         s21,s22,zexps,gami(1,ii1,j),beti(1,ii1,j),
     .         gami(1,ii1,j1),beti(1,ii1,j1),isp(j),iss(j),
     .         Vmat(1,1,j,ii),ndv,Wmat(1,j1,ii),iiwww,rhorat(k))
         else
c: No mismatch at bottom of layer.  Call rp_nomm:
            call rp_nomm(Vmat(1,1,j1,ii),e11,e12,e21,e22,zexpe,
     .         s11,s12,s21,s22,zexps,gami(1,ii1,j),beti(1,ii1,j),
     .         gami(1,ii1,j1),beti(1,ii1,j1),isp(j),iss(j),
     .         Vmat(1,1,j,ii),ndv,Wmat(1,j1,ii),iiwww)
         endif
      enddo
c
      if(allf(ii) .eq. 1) then
c: Propagate across fluid-fluid halfspace interface to bottom of last 
c: (fluid) layer:
         iiwww=iiw*iiww(ns)
         j1=jsol2
         k=j1 + jflu(ii,4)
         call rp_flay(Vmat(1,1,ns,ii),e11,e12,e21,e22,zexpe,
     .      gami(1,ii2,jflu2),gami(1,ii1,ns),1,0.d0,rhorat(k),
     .      mm(k),iiwww,ndv,Vmat(1,1,j1,ii),Wmat(1,ns,ii),jjfail)
         if(jjfail .gt. 0) return
      else
c: Propagate across solid-fluid interface to bottom of first fluid layer:
         j1=jflu2
         iiwww=iiw*iiww(jsol1)
         call rp_sfint(Vmat(1,1,jsol1,ii),Alay(1,ii),Blay(1,ii),
     .      ikcon,rholay(ii),gami(1,ii2,j1),gami(1,ii1,jsol1),
     .      beti(1,ii1,jsol1),Vmat(1,1,j1,ii),ndv,Wmat(1,jsol1,ii),
     .      iiwww)
      endif
c
      do j=jflu2,jflu1+inc,-inc
         j1=j - inc
c: Propagate through fluid layers from bottom of j'th layer to top of
c: j1'th layer (the next layer toward the svp minimum from the j'th):
         k=j + jflu(ii,4) - inc
         iiwww=iiw*iiww(j)
c: Compute 2x2 propagator matrices for compressional and shear profiles:
         call ep_calc(xkh,w,xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .      etasq(j),gami(1,ii2,j),h(j),e11,e12,e21,e22,zexpe,isp(j),
     .      ihf(j),2,ndv,iiw,ailay(1,1,1,j),bilay(1,1,1,j),
     .      zetalay(1,1,j),aisoln(1,j),ii1,ii2)
         call rp_flay(Vmat(1,1,j,ii),e11,e12,e21,e22,zexpe,
     .      gami(1,ii2,j1),gami(1,ii2,j),isp(j),h(j),rhorat(k),
     .      mm(k),iiwww,ndv,Vmat(1,1,j1,ii),Wmat(1,j,ii),jjfail)
         if(jjfail .gt. 0) return
      enddo
c
c: Do last layer (with no interface) if we are doing bottom half and ref depth
c: at top of layer or if we are doing top half and ref depth is at bottom:
      if(isvmin .eq. ii) then
         j=jflu1
         j1=j - inc
         iiwww=iiw*iiww(j)
         call ep_calc(xkh,w,xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .      etasq(j),gami(1,ii2,j),h(j),e11,e12,e21,e22,zexpe,isp(j),
     .      ihf(j),2,ndv,iiw,ailay(1,1,1,j),bilay(1,1,1,j),
     .      zetalay(1,1,j),aisoln(1,j),ii1,ii2)
         call rp_flay(Vmat(1,1,j,ii),e11,e12,e21,e22,zexpe,
     .      gami(1,ii1,j),gami(1,ii2,j),isp(j),h(j),1.d0,0,iiwww,ndv,
     .      Vmat(1,1,j1,ii),Wmat(1,j,ii),jjfail)
         if(jjfail .gt. 0) return
      endif
c
c: Copy desired reflection coefficient and its derivatives from last Vmat:
      do j=1,ndv
         R(j)=Vmat(j,1,j1,ii)
      enddo
c
      return
      end
      subroutine mode_traj(k,r1r2,ii)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 ii
      complex*16 k,r1r2(3,4)
c
cxx   write(21,100) dreal(k)/kw0,dimag(k)/kw0,ii,iish(1,1),iish(1,2),
      write(21,100) dreal(k)/kw0,8685.9*dimag(k),ii,iish(1,1),
     .   iish(1,2),iish(2,1),iish(2,2),dreal(r1r2(1,4)),
     .   dimag(r1r2(1,4)),kduct
100   format(e14.8,1x,e14.8,1x,i2,4(1x,i2),1x,e9.3,1x,f7.2,1x,i1)
c
      return
      end
ccc
      subroutine mode_traj_bp(rr0)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 ii,jj
      complex*16 rr0(3,4)
c
      do ii=1,2
         do jj=1,2
            if(real(xkbp(ii,jj)) .ge. kremin .and. 
     .         real(xkbp(ii,jj)) .le. w/cphlo .and.
     .         xkbp(ii,jj) .ne. (0.d0,0.d0)) then
               call mode_traj(xkbp(ii,jj),rr0,-5)
            endif
         enddo
      enddo
c
      return
      end
      subroutine rp_flay(R,e11,e12,e21,e22,ze,gami1,gami2,isoe,h,
     .   rhorat,iimm,iiw,ndv,V,Wflu,jjfail)
c
c: Finds the reflection coefficient just above a fluid-fluid interface
c: (characterized by density ratio rho1/rho2=rhorat) given the 
c: reflection coefficient at the bottom of layer 2 below and the 
c: propagator matrix E for the fluid layer 2.  gami1 and gami2 are i*(the
c: vertical wavenumbers) at the bottom of layer 1 and bottom of layer 2
c: respectively.  For iiw=1, the transmission coefficient W at the bottom
c: of layer 2 is also computed.  Cases where the thickness of later 2 is
c: zero, the density ratio is one, and layer 2 is a halfspace (R=0) are 
c: treated.
c: Note that ze is exponential factor common to E.
c
      implicit complex*16(a-z)
      integer*4 isoe,ndv,j,iiw,iimm,jjfail
      complex*16 e11(3),e12(3),e21(3),e22(3),R(3),V(3),Wflu(6),
     .   gami1(3),gami2(3),zzero
      real*8 h,rhorat
      data zzero/(0.d0,0.d0)/
c
      if(h .eq. 0.d0) then
c: Interface only:
cxx      if(gami1(1) .eq. gami2(1) .and. rhorat .eq. 1.d0) then
         if(iimm .eq. 0) then
c: Zero thickness, no-mismatch case: V=R, W=1:
            V(1)=R(1)
            V(2)=R(2)
            V(3)=R(3)
            Wflu(1)=dcmplx(1.d0,0.d0)
            Wflu(5)=dcmplx(0.d0,0.d0)
            return
         elseif(R(1) .ne. zzero) then
            j1=-gami1(1)*(1.d0 + R(1))
            k1=gami2(1)*(1.d0 - R(1))
         else
c: Halfspace below:
            j1=-gami1(1)
            k1=gami2(1)
         endif
      else
         if(isoe .eq. 0) then
            game11=gami2(1)*e11(1)
            game12=gami2(1)*e12(1)
            e1m=game11 - e21(1)
            e1p=game11 + e21(1)
            e2m=game12 - e22(1)
            e2p=game12 + e22(1)
         else
c: Special expressions for isospeed layer:
            if(rhorat .eq. 1.d0 .and. gami1(1) .eq. gami2(1)) then
c: Special case for isospeed layer in continuous medium (see. p.124):
               num=R(1)*e11(1)
               if(cdabs(e21(1)).eq.0.d0) then
                  jjfail=1
                  return
               end if
               V(1)=num/e21(1)
               if(iiw .eq. 1) then
                  Wflu(1)=gami1(1)/e21(1)
                  Wflu(5)=-ze
               endif
               do j=2,ndv
                  nump=R(1)*e11(j) + R(j)*e11(1)
                  V(j)=(nump - V(1)*e21(j))/e21(1)
               enddo
               return
            endif
            e1p=e11(1)
            e1m=e21(1)
            e2p=e12(1)
            e2m=e22(1)
         endif
         if(R(1) .ne. zzero) then
            term2=e2m - R(1)*e2p
            j1=gami1(1)*term2
            k1=e1m - R(1)*e1p
         else
            j1=gami1(1)*e2m
            k1=e1m
         endif
      endif
c
      if(rhorat .ne. 1.d0) k1=rhorat*k1
      num=j1 + k1
      den=j1 - k1
cpln Added 03/02-2000 for soft layers
cpln taken from PROSIM in rx_rp_flay
      if(cdabs(den). eq. 0.d0) then
         V(1)=(0.d0,0.d0) 
         Wflu(1)=(0.d0,0.d0) 
         do j=2,ndv
            V(j)=(1.d0,0.d0) 
            Wflu(j)=(1.d0,0.d0) 
         enddo
         print *,'Info msg: D=0 in rp_flay'
         jjfail=1
         return
      endif
      V(1)=num/den
      if(iiw .eq. 1) then
         Wflu(1)=-2.d0*gami1(1)/den
         Wflu(5)=-ze
         if(rhorat .ne. 1.d0) Wflu(1)=rhorat*Wflu(1)
      endif
c
c: If first derivatives desired, compute them now:
      do j=2,ndv
         if(h .eq. 0.d0) then
c: Interface only:
            if(R(1) .ne. zzero) then
               j1p=-gami1(1)*R(j) - gami1(j)*(1.d0 + R(1))
               k1p=-gami2(1)*R(j) + gami2(j)*(1.d0 - R(1))
            else
c: Halfspace below:
               j1p=-gami1(j)
               k1p=gami2(j)
            endif
         else
            if(isoe .eq. 0) then
               game11p=gami2(1)*e11(j) + gami2(j)*e11(1)
               game12p=gami2(1)*e12(j) + gami2(j)*e12(1)
               e1mp=game11p - e21(j)
               e1pp=game11p + e21(j)
               e2mp=game12p - e22(j)
               e2pp=game12p + e22(j)
            else
c: Special expressions for isospeed layer:
               e1pp=e11(j)
               e1mp=e21(j)
               e2pp=e12(j)
               e2mp=e22(j)
            endif
            if(R(1) .ne. zzero) then
               term2p=e2mp - R(1)*e2pp - R(j)*e2p
               j1p=gami1(1)*term2p + gami1(j)*term2
               k1p=e1mp - R(1)*e1pp - R(j)*e1p
            else
               j1p=gami1(1)*e2mp + gami1(j)*e2m
               k1p=e1mp
            endif
         endif
c
         if(rhorat .ne. 1.d0) k1p=rhorat*k1p
         nump=j1p + k1p
         denp=j1p - k1p
         V(j)=(den*nump - denp*num)/(den*den)
      enddo
c
      return
      end
       subroutine rp_nomm(Rmat,e11,e12,e21,e22,ze,s11,s12,s21,s22,
     .   zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,Wmat,iiw)
c
c: Same as rp_slay, except the interface does not exist.  Formulas
c: derived from rp_slay by setting P=Q=1, U=V=0.
c
      implicit complex*16(a-z)
      integer*4 isoe,isos,ndv,j,iiw
      complex*16 e11(3),e12(3),e21(3),e22(3),s11(3),s12(3),s21(3),
     .   s22(3),Rmat(3,5),Vmat(3,5),Wmat(6)
      complex*16 gami1(3),beti1(3),gami2(3),beti2(3),zzero
      data zzero/(0.d0,0.d0)/
c
      if(Rmat(1,4) .eq. zzero .or. Rmat(1,2) .eq. zzero) then
         call rp_nomm0(Rmat,e11,e12,e21,e22,ze,s11,s12,s21,s22,zs,
     .      gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,Wmat,iiw)
         return
      endif
c
      q1=1.d0 + Rmat(1,1)
      q2diff=1.d0 - Rmat(1,1)
      q2=gami2(1)*q2diff
      q3=Rmat(1,2)
      q4=-beti2(1)*Rmat(1,2)
      r1=Rmat(1,4)
      r2=-gami2(1)*Rmat(1,4)
      r3=1.d0 + Rmat(1,3)
      r4diff=1.d0 - Rmat(1,3)
      r4=beti2(1)*r4diff
c
c: rps and rsp are multiplied by exp(zsp) (see r1,r2,q3,q4):
      spfac=cdexp(Rmat(1,5))
      r1sp=spfac*r1
      r2sp=spfac*r2
      q3sp=spfac*q3
      q4sp=spfac*q4
      delta1=-2.d0*r2
      delta2=2.d0*q4
      delta3=q4sp*r1sp - q1*r4
      delta4=q3sp*r2sp - q2*r3
c
      x1=delta4*q4sp
      x2=delta3*q3sp
      x3=-delta4*q1
      x4=-delta3*q2
c
      w1=delta2*q2
      w2=delta2*q1
      w3=delta1*q3sp
      w4=delta1*q4sp
c
      if(isoe .eq. 0) then
         e12g=gami1(1)*e12(1)
         e22g=gami1(1)*e22(1)
         e1p=e11(1) + e12g
         e1m=e11(1) - e12g
         e2p=e21(1) + e22g
         e2m=e21(1) - e22g
      else
         e1p=e11(1)
         e1m=e12(1)
         e2p=e21(1)
         e2m=e22(1)
      endif
      n1=w1*e1p - w2*e2p
      n2=x1*e1p - x2*e2p
      k1=w1*e1m - w2*e2m
      k2=x1*e1m - x2*e2m
c
      if(isos .eq. 0) then
         s12b=beti1(1)*s12(1)
         s22b=beti1(1)*s22(1)
         s1p=s11(1) + s12b
         s1m=s11(1) - s12b
         s2p=s21(1) + s22b
         s2m=s21(1) - s22b
      else
         s1p=s11(1)
         s1m=s12(1)
         s2p=s21(1)
         s2m=s22(1)
      endif
      m1=w3*s2p - w4*s1p
      m2=x3*s2p - x4*s1p
      y1=w3*s2m - w4*s1m
      y2=x3*s2m - x4*s1m
c
      den=(y2*k1 - y1*k2)
      deninv=1.d0/den
      num1=-(y2*n1 - y1*n2)
      Vmat(1,1)=num1*deninv
      num2=(w1*x2 - w2*x1)
      gamnum2=gami1(1)*num2
      Vmat(1,2)=2.d0*gamnum2*deninv
      num3=(k2*m1 - k1*m2)
      Vmat(1,3)=num3*deninv
      num4=(w4*x3 - w3*x4)
      betnum4=beti1(1)*num4
      Vmat(1,4)=-2.d0*betnum4*deninv
c
c: Exponential factor [exp(zexp_sp)] to be multiplied by vps and vsp:
      Vmat(1,5)=-(ze+zs)
c
      if(iiw .eq. 1) then
         a7=e2m*q1 - e1m*q2
         a8=e2m*r1sp - e1m*r2sp
         b7=s1m*q4sp - s2m*q3sp
         b8=s1m*r4 - s2m*r3
         wden=a7*b8 - b7*a8
         fac=-2.d0/wden
         a9_den=fac*gami1(1)
         Wmat(1)=a9_den*b8
         Wmat(2)=-a9_den*b7
         Wmat(5)=-ze
         a9bar_den=fac*beti1(1)
         Wmat(3)=-a9bar_den*a7
         Wmat(4)=a9bar_den*a8
         Wmat(6)=-zs
      endif
c
c: If first derivatives desired, compute them now:
      do j=2,ndv
         q1p=Rmat(j,1)
         q2p=-gami2(1)*Rmat(j,1) + gami2(j)*q2diff
         q3p=Rmat(j,2)
         q4p=-beti2(1)*Rmat(j,2) - beti2(j)*Rmat(1,2)
         r1p=Rmat(j,4)
         r2p=-gami2(1)*Rmat(j,4) - gami2(j)*Rmat(1,4)
         r3p=Rmat(j,3)
         r4p=-beti2(1)*Rmat(j,3) + beti2(j)*r4diff
         r1spp=spfac*r1p
         r2spp=spfac*r2p
         q3spp=spfac*q3p
         q4spp=spfac*q4p
c
         delta1p=-2.d0*r2p
         delta2p=2.d0*q4p
         delta3p=q4sp*r1spp + q4spp*r1sp - q1*r4p - q1p*r4
         delta4p=q3sp*r2spp + q3spp*r2sp - q2*r3p - q2p*r3
c
         x1p=delta4*q4spp + delta4p*q4sp
         x2p=delta3*q3spp + delta3p*q3sp
         x3p=-delta4*q1p - delta4p*q1
         x4p=-delta3*q2p - delta3p*q2
c
         w1p=delta2*q2p + delta2p*q2
         w2p=delta2*q1p + delta2p*q1
         w3p=delta1*q3spp + delta1p*q3sp
         w4p=delta1*q4spp + delta1p*q4sp
c
         if(isoe .eq. 0) then
            e12gp=gami1(1)*e12(j) + gami1(j)*e12(1)
            e22gp=gami1(1)*e22(j) + gami1(j)*e22(1)
            e1pp=e11(j) + e12gp
            e1mp=e11(j) - e12gp
            e2pp=e21(j) + e22gp
            e2mp=e21(j) - e22gp
         else
            e1pp=e11(j)
            e1mp=e12(j)
            e2pp=e21(j)
            e2mp=e22(j)
         endif
         n1p=w1*e1pp + w1p*e1p - w2*e2pp - w2p*e2p
         n2p=x1*e1pp + x1p*e1p - x2*e2pp - x2p*e2p
         k1p=w1*e1mp + w1p*e1m - w2*e2mp - w2p*e2m
         k2p=x1*e1mp + x1p*e1m - x2*e2mp - x2p*e2m
c
         if(isos .eq. 0) then
            s12bp=beti1(1)*s12(j) + beti1(j)*s12(1)
            s22bp=beti1(1)*s22(j) + beti1(j)*s22(1)
            s1pp=s11(j) + s12bp
            s1mp=s11(j) - s12bp
            s2pp=s21(j) + s22bp
            s2mp=s21(j) - s22bp
         else
            s1pp=s11(j)
            s1mp=s12(j)
            s2pp=s21(j)
            s2mp=s22(j)
         endif
         m1p=w3*s2pp + w3p*s2p - w4*s1pp - w4p*s1p
         m2p=x3*s2pp + x3p*s2p - x4*s1pp - x4p*s1p
         y1p=w3*s2mp + w3p*s2m - w4*s1mp - w4p*s1m
         y2p=x3*s2mp + x3p*s2m - x4*s1mp - x4p*s1m
c
         denp=(y2*k1p + y2p*k1 - y1*k2p - y1p*k2)
         deninvp=-denp*deninv*deninv
         num1p=-(y2*n1p + y2p*n1 - y1*n2p - y1p*n2)
         Vmat(j,1)=num1*deninvp + num1p*deninv
         num2p=(w1*x2p + w1p*x2 - w2*x1p - w2p*x1)
         gamnum2p=gami1(1)*num2p + gami1(j)*num2
         Vmat(j,2)=2.d0*(gamnum2*deninvp + gamnum2p*deninv)
         num3p=(k2*m1p + k2p*m1 - k1*m2p - k1p*m2)
         Vmat(j,3)=num3*deninvp + num3p*deninv
         num4p=(w4*x3p + w4p*x3 - w3*x4p - w3p*x4)
         betnum4p=beti1(1)*num4p + beti1(j)*num4
         Vmat(j,4)=-2.d0*(betnum4*deninvp + betnum4p*deninv)
      enddo
c
      return
      end
      subroutine rp_nomm0(Rmat,e11,e12,e21,e22,ze,s11,s12,s21,s22,
     .   zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,Wmat,iiw)
c
c: Same as rp_slay0, except no interface exists at bottom.  Formulas
c: derived from rp_slay0 by setting P=Q=1 and U=V=0.
c
      implicit complex*16 (a-z)
      integer*4 isoe,isos,ndv,j,iiw,iiz
      complex*16 Rmat(3,5),e11(3),e12(3),e21(3),e22(3),s11(3),
     .   s12(3),s21(3),s22(3),Vmat(3,5),gami1(3),gami2(3),beti1(3),
     .   beti2(3),Wmat(6),zzero
      data zzero/(0.d0,0.d0)/
c
      if(isoe .eq. 0) then
         e12g=gami1(1)*e12(1)
         e22g=gami1(1)*e22(1)
         e1p=e11(1) + e12g
         e1m=e11(1) - e12g
         e2p=e21(1) + e22g
         e2m=e21(1) - e22g
      else
         e1p=e11(1)
         e1m=e12(1)
         e2p=e21(1)
         e2m=e22(1)
      endif
c
      if(isos .eq. 0) then
         s12b=beti1(1)*s12(1)
         s22b=beti1(1)*s22(1)
         s1p=s11(1) + s12b
         s1m=s11(1) - s12b
         s2p=s21(1) + s22b
         s2m=s21(1) - s22b
      else
         s1p=s11(1)
         s1m=s12(1)
         s2p=s21(1)
         s2m=s22(1)
      endif
c
      if(Rmat(1,1) .eq. zzero .and. Rmat(1,3) .eq. zzero) then
         iiz=1
         q2=gami2(1)
         r4=beti2(1)
         m1=q2*e1p - e2p
         f1=q2*e1m - e2m
         t2=r4*s1p - s2p
         c2=r4*s1m - s2m
      else
         iiz=0
         q1=1.d0 + Rmat(1,1)
         q2diff=1.d0 - Rmat(1,1)
         q2=gami2(1)*q2diff
         r3=1.d0 + Rmat(1,3)
         r4diff=1.d0 - Rmat(1,3)
         r4=beti2(1)*r4diff
         m1=q2*e1p - q1*e2p
         f1=q2*e1m - q1*e2m
         t2=r4*s1p - r3*s2p
         c2=r4*s1m - r3*s2m
      endif
c
      Vmat(1,1)=-m1/f1
      Vmat(1,2)=zzero
      Vmat(1,3)=-t2/c2
      Vmat(1,4)=zzero
c
c: vps(1:2) and vsp(1:2) are to be multiplied by exp(zexp_sp):
      Vmat(1,5)=zzero
c
      if(iiw .eq. 1) then
         Wmat(1)=2.d0*gami1(1)/f1
         Wmat(2)=zzero
         Wmat(5)=-ze
         Wmat(3)=2.d0*beti1(1)/c2
         Wmat(4)=zzero
         Wmat(6)=-zs
      endif
c
      do j=2,ndv
         if(isoe .eq. 0) then
            e12gp=gami1(1)*e12(j) + gami1(j)*e12(1)
            e22gp=gami1(1)*e22(j) + gami1(j)*e22(1)
            e1pp=e11(j) + e12gp
            e1mp=e11(j) - e12gp
            e2pp=e21(j) + e22gp
            e2mp=e21(j) - e22gp
         else
            e1pp=e11(j)
            e1mp=e12(j)
            e2pp=e21(j)
            e2mp=e22(j)
         endif
c
         if(isos .eq. 0) then
            bets12p=beti1(1)*s12(j) + beti1(j)*s12(1)
            bets22p=beti1(1)*s22(j) + beti1(j)*s22(1)
            s1pp=s11(j) + bets12p
            s1mp=s11(j) - bets12p
            s2pp=s21(j) + bets22p
            s2mp=s21(j) - bets22p
         else
            s1pp=s11(j)
            s1mp=s12(j)
            s2pp=s21(j)
            s2mp=s22(j)
         endif
c
         if(iiz .eq. 1) then
            q2p=gami2(j)
            r4p=beti2(j)
            m1p=q2*e1pp + q2p*e1p - e2pp
            f1p=q2*e1mp + q2p*e1m - e2mp
            t2p=r4*s1pp + r4p*s1p - s2pp
            c2p=r4*s1mp + r4p*s1m - s2mp
         else
            q1p=Rmat(j,1)
            q2diffp=-Rmat(j,1)
            q2p=gami2(1)*q2diffp + gami2(j)*q2diff
            r3p=Rmat(j,3)
            r4diffp=-Rmat(j,3)
            r4p=beti2(1)*r4diffp + beti2(j)*r4diff
            m1p=q2*e1pp + q2p*e1p - q1*e2pp - q1p*e2p
            f1p=q2*e1mp + q2p*e1m - q1*e2mp - q1p*e2m
            t2p=r4*s1pp + r4p*s1p - r3*s2pp - r3p*s2p
            c2p=r4*s1mp + r4p*s1m - r3*s2mp - r3p*s2m
         endif
c
         Vmat(j,1)=-(m1p/f1 + Vmat(1,1)*f1p/f1)
         Vmat(j,2)=zzero
         Vmat(j,3)=-(t2p/c2 + Vmat(1,3)*c2p/c2)
         Vmat(j,4)=zzero
      enddo
c
      return
      end
      subroutine rp_sfint(Rmat,Acon1,Bcon1,ikcon,rhorat,gami1,
     .   gami2,beti2,R,ndv,Wmat,iiw)
c
c: Finds the p-p reflection coefficient just above a fluid-solid interface
c: given p-p, p-s, s-s, and s-p reflection coefficients just below the
c: fluid-solid interface.  The interface is characterized by P,Q,U,V
c: gami2,beti2 are the p- and w-wave vertical wavenumbers in the solid,
c: and gami1 is the p-wave vertical wavenumber in the fluid.
c:
      implicit complex*16(a-z)
      integer*4 ndv,j,iiw,iiz,iiz0,iir
      complex*16 Rmat(3,5),Acon1(3),Bcon1(3),gami1(3),
     .   gami2(3),beti2(3),R(3),Wmat(6),ikcon(3),zzero
      real*8 rhorat,magsq,deps,rat
      data zzero/(0.d0,0.d0)/,deps/1.d-50/
c
      if(Rmat(1,1) .eq. zzero .and. Rmat(1,3) .eq. zzero) then
         iiz0=1
         q1=(1.d0,0.d0)
         q2=gami2(1)
         r3=(1.d0,0.d0)
         r4=beti2(1)
      else
         iiz0=0
         q1=1.d0 + Rmat(1,1)
         q2diff=1.d0 - Rmat(1,1)
         q2=gami2(1)*q2diff
         r3=1.d0 + Rmat(1,3)
         r4diff=1.d0 - Rmat(1,3)
         r4=beti2(1)*r4diff
      endif
c
      q2_rho=q2*rhorat
      ik_rho=ikcon(1)*rhorat
      a2gam=Acon1(1)*gami1(1)
      b2gam=Bcon1(1)*gami1(1)
      q1_a2gam=q1*a2gam
c
      c1=q2_rho + q1_a2gam
      d1=q2_rho - q1_a2gam
      f1=q2*Bcon1(1)
c: Case where vps and vsp below are zero:
cxx   if(Rmat(1,2) .eq. zzero .and. Rmat(1,4) .eq. zzero) then
      if(magsq(Rmat(1,2)) .lt. deps .or. 
     .   magsq(Rmat(1,4)) .lt. deps) then
         iiz=1
         delta3=-q1*r4
         delta4=-q2*r3
         d3=-q1*ik_rho
         f3=-q1*Acon1(1)
         d4=q2*b2gam
c
         j1=c1
         k1=d1
         m1=f1
      else
         iiz=0
c: rps and rsp are multiplied by exp(zsp) (see r1,r2,q3,q4):
         spfac=cdexp(Rmat(1,5))
         r1sp=Rmat(1,4)*spfac
         r2sp=-gami2(1)*r1sp
         q3sp=Rmat(1,2)*spfac
         q4sp=-beti2(1)*q3sp
c: exp(zsp) factor must be accounted for in delta3,delta4, etc:
         delta1=2.d0*gami2(1)*Rmat(1,4)
         delta2=-2.d0*beti2(1)*Rmat(1,2)
         delta3=q4sp*r1sp - q1*r4
         delta4=q3sp*r2sp - q2*r3
c
         q4_b2gam=q4sp*b2gam
         q3_ikrho=q3sp*ik_rho
         c2=q4_b2gam - q3_ikrho
         d2=q4_b2gam + q3_ikrho
         f2=q3sp*Acon1(1)
c
         d3=q4sp*rhorat - q1*ik_rho
         f3=q4sp*Bcon1(1) - q1*Acon1(1)
         d4=q3sp*a2gam + q2*b2gam
c
         j1=delta2*c1 - delta1*c2
         k1=delta2*d1 + delta1*d2
         m1=delta2*f1 + delta1*f2
      endif
c
      del4_d3=delta4*d3
      del3_d4=delta3*d4
      j2=del4_d3 + del3_d4
      k2=del4_d3 - del3_d4
      m2=delta4*f3
c
      den=m2*k1 - m1*k2
      num=m1*j2 - m2*j1
      R(1)=num/den
      if(iiw .eq. 1) then
         h2=r3*ik_rho + r4*b2gam
         if(iiz .eq. 1) then
            h4=r3*Acon1(1)
            h5=-h2
            g4=q2*Bcon1(1)
            g5=-d1
         else
            h1=r1sp*a2gam - r2sp*rhorat
            h4=r2sp*Bcon1(1) + r3*Acon1(1)
            h5=h1 - h2
            g4=q2*Bcon1(1) + q3sp*Acon1(1)
            g5=-(d1 + d2)
         endif
         L5=2.d0*gami1(1)*rhorat*(Acon1(1) - ikcon(1)*Bcon1(1))
c: OLD:
cxx      wden=h4*g5 - h5*g4
cxx      fac=L5/wden
cxx      Wmat(1)=h4*fac
cxx      Wmat(2)=-g4*fac
c: NEW:
         wden2=round_check(h4*g5,-h5*g4,1.d0/L5,iir,rat)
         if(iir .eq. 1) then
            Wmat(1)=zzero
            Wmat(2)=zzero
cxx         print *,'iir=1 in rp_sfint: '
         else
            Wmat(1)=h4/wden2
            Wmat(2)=-g4/wden2
         endif
         Wmat(3)=(0.d0,0.d0)
         Wmat(4)=(0.d0,0.d0)
         Wmat(5)=(0.d0,0.d0)
         Wmat(6)=(0.d0,0.d0)
      endif
c
c: Derivatives:
      do j=2,ndv
         if(iiz0 .eq. 1) then
            q1p=(0.d0,0.d0)
            q2p=gami2(j)
            r3p=(0.d0,0.d0)
            r4p=beti2(j)
         else
            q1p=Rmat(j,1)
            q2p=-gami2(1)*Rmat(j,1) + gami2(j)*q2diff
            r3p=Rmat(j,3)
            r4p=-beti2(1)*Rmat(j,3) + beti2(j)*r4diff
         endif
c
         q2_rhop=q2p*rhorat
         ik_rhop=ikcon(j)*rhorat
         a2gamp=Acon1(1)*gami1(j) + Acon1(j)*gami1(1)
         b2gamp=Bcon1(1)*gami1(j) + Bcon1(j)*gami1(1)
         q1_a2gamp=q1*a2gamp + q1p*a2gam
c
         c1p=q2_rhop + q1_a2gamp
         d1p=q2_rhop - q1_a2gamp
         f1p=q2*Bcon1(j) + q2p*Bcon1(1)
c
         if(iiz .eq. 1) then
            delta3p=-q1*r4p - q1p*r4
            delta4p=-q2*r3p - q2p*r3
            d3p=-q1*ik_rhop - q1p*ik_rho
            f3p=-q1*Acon1(j) - q1p*Acon1(1)
            d4p=q2*b2gamp + q2p*b2gam
c
            j1p=c1p
            k1p=d1p
            m1p=f1p
         else
            r1spp=Rmat(j,4)*spfac
            r2spp=(-gami2(1)*Rmat(j,4) - gami2(j)*Rmat(1,4))*spfac
            q3spp=Rmat(j,2)*spfac
            q4spp=(-beti2(1)*Rmat(j,2) - beti2(j)*Rmat(1,2))*spfac
c
            delta1p=2.d0*(gami2(1)*Rmat(j,4) + gami2(j)*Rmat(1,4))
            delta2p=-2.d0*(beti2(1)*Rmat(j,2) + beti2(j)*Rmat(1,2))
            delta3p=q4sp*r1spp + q4spp*r1sp - q1*r4p - q1p*r4
            delta4p=q3sp*r2spp + q3spp*r2sp - q2*r3p - q2p*r3
c
            q4_b2gamp=q4sp*b2gamp + q4spp*b2gam
            q3_ikrhop=q3sp*ik_rhop + q3spp*ik_rho
            c2p=q4_b2gamp - q3_ikrhop
            d2p=q4_b2gamp + q3_ikrhop
            f2p=q3sp*Acon1(j) + q3spp*Acon1(1)
c
            d3p=q4spp*rhorat - q1*ik_rhop - q1p*ik_rho
            f3p=q4sp*Bcon1(j) + q4spp*Bcon1(1) - q1*Acon1(j) - 
     .         q1p*Acon1(1)
            d4p=q3sp*a2gamp + q3spp*a2gam + q2*b2gamp + q2p*b2gam
c
            j1p=delta2*c1p + delta2p*c1 - delta1*c2p - delta1p*c2
            k1p=delta2*d1p + delta2p*d1 + delta1*d2p + delta1p*d2
            m1p=delta2*f1p + delta2p*f1 + delta1*f2p + delta1p*f2
         endif
c
         del4_d3p=delta4*d3p + delta4p*d3
         del3_d4p=delta3*d4p + delta3p*d4
         j2p=del4_d3p + del3_d4p
         k2p=del4_d3p - del3_d4p
         m2p=delta4*f3p + delta4p*f3
c
         denp=m2*k1p + m2p*k1 - m1*k2p - m1p*k2
         nump=m1*j2p + m1p*j2 - m2*j1p - m2p*j1
cxx      R(j)=(den*nump - denp*num)/(den*den)
         R(j)=(nump-denp*R(1))/den
      enddo
c
      return
      end
       subroutine rp_slay(Rmat,p,q,u,v,e11,e12,e21,e22,
     .   ze,s11,s12,s21,s22,zs,gami1,beti1,gami2,beti2,
     .   isoe,isos,Vmat,ndv,Wmat,iiw,rhorat)
c
c: Finds the four reflection coefficients referenced to the top of
c: layer 1 given the four reflection coefficients at the top of layer
c: layer 2, the (p,q,u,v) constants for the 1-2 layer interface, and
c: the compressional and shear wave propagator matrices e and s.
c: gam1,gam2 and bet1,bet2 are the vertical comp and shear wavenumbers 
c: at the top of layers 1 & 2.
c: Note that ze and zs are exponential factors common to eij and sij.
c: vps and vsp need to account for them.
c
      implicit complex*16(a-z)
      integer*4 isoe,isos,ndv,j,iiw
      complex*16 e11(3),e12(3),e21(3),e22(3),s11(3),s12(3),s21(3),
     .   s22(3),Rmat(3,5),p(3),q(3),u(3),v(3),Vmat(3,5),Wmat(6)
      complex*16 gami1(3),beti1(3),gami2(3),beti2(3),zzero
      real*8 rhorat,magsq,deps
      data zzero/(0.d0,0.d0)/,deps/1.d-50/
c
cxx   if(Rmat(1,4) .eq. zzero .or. Rmat(1,2) .eq. zzero) then
      if(magsq(Rmat(1,2)) .lt. deps .or.
     .   magsq(Rmat(1,4)) .lt. deps) then
         call rp_slay0(Rmat,p,q,u,v,e11,e12,e21,e22,ze,s11,s12,
     .      s21,s22,zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,
     .      ndv,Wmat,iiw,rhorat)
         return
      endif
c
      q1=1.d0 + Rmat(1,1)
      q2diff=1.d0 - Rmat(1,1)
      q2=gami2(1)*q2diff
      q3=Rmat(1,2)
      q4=-beti2(1)*Rmat(1,2)
      r1=Rmat(1,4)
      r2=-gami2(1)*Rmat(1,4)
      r3=1.d0 + Rmat(1,3)
      r4diff=1.d0 - Rmat(1,3)
      r4=beti2(1)*r4diff
c
c: rps and rsp are multiplied by exp(zsp) (see r1,r2,q3,q4):
      spfac=cdexp(Rmat(1,5))
      r1sp=spfac*r1
      r2sp=spfac*r2
      q3sp=spfac*q3
      q4sp=spfac*q4
c: Multiply out to get simpler form:
cxx1  delta1=q2*r1 - q1*r2
cxx1  delta2=q4*r3 - q3*r4
      delta1=-2.d0*r2
      delta2=2.d0*q4
      delta3=q4sp*r1sp - q1*r4
      delta4=q3sp*r2sp - q2*r3
c
      f1=q4sp*p(1) + q1*v(1)
      f2=q4sp*u(1) - q1*q(1)
      f3=q3sp*q(1) + q2*u(1)
      f4=q3sp*v(1) - q2*p(1)
c: Note: I switched x2 and x3 from notes:
      x1=delta4*f1
      x2=delta3*f3
      x3=delta4*f2
      x4=delta3*f4
c
      del2q1=delta2*q1
      del2q2=delta2*q2
      del1q3=delta1*q3sp
      del1q4=delta1*q4sp
      w1=del2q2*p(1) - del1q3*v(1)
      w2=del2q1*q(1) - del1q4*u(1)
      w3=del2q2*u(1) + del1q3*q(1)
      w4=del2q1*v(1) + del1q4*p(1)
c
      if(isoe .eq. 0) then
         e12g=gami1(1)*e12(1)
         e22g=gami1(1)*e22(1)
         e1p=e11(1) + e12g
         e1m=e11(1) - e12g
         e2p=e21(1) + e22g
         e2m=e21(1) - e22g
      else
         e1p=e11(1)
         e1m=e12(1)
         e2p=e21(1)
         e2m=e22(1)
      endif
      n1=w1*e1p - w2*e2p
      n2=x1*e1p - x2*e2p
      k1=w1*e1m - w2*e2m
      k2=x1*e1m - x2*e2m
c
      if(isos .eq. 0) then
         s12b=beti1(1)*s12(1)
         s22b=beti1(1)*s22(1)
         s1p=s11(1) + s12b
         s1m=s11(1) - s12b
         s2p=s21(1) + s22b
         s2m=s21(1) - s22b
      else
         s1p=s11(1)
         s1m=s12(1)
         s2p=s21(1)
         s2m=s22(1)
      endif
      m1=w3*s2p - w4*s1p
      m2=x3*s2p - x4*s1p
      y1=w3*s2m - w4*s1m
      y2=x3*s2m - x4*s1m
c
      den=(y2*k1 - y1*k2)
      deninv=1.d0/den
      num1=-(y2*n1 - y1*n2)
      Vmat(1,1)=num1*deninv
      num2=(w1*x2 - w2*x1)
      gamnum2=gami1(1)*num2
      Vmat(1,2)=2.d0*gamnum2*deninv
      num3=(k2*m1 - k1*m2)
      Vmat(1,3)=num3*deninv
      num4=(w4*x3 - w3*x4)
      betnum4=beti1(1)*num4
      Vmat(1,4)=-2.d0*betnum4*deninv
c
c: Exponential factor [exp(zexp_sp)] to be multiplied by vps and vsp:
      Vmat(1,5)=-(ze+zs)
c
      if(iiw .eq. 1) then
         a7=e1m*f4 - e2m*f2
         a8=e2m*(r1sp*q(1) - r4*u(1)) - e1m*(r2sp*p(1) - r3*v(1))
         b7=s1m*f1 - s2m*f3
         b8=s1m*(r1sp*v(1) + r4*p(1)) - s2m*(r2sp*u(1) + r3*q(1))
         wden=a7*b8 - b7*a8
         fac=-2.d0*rhorat/wden
cxx      a9_den=fac*gami1(1)*cdexp(-ze)
         a9_den=fac*gami1(1)
         Wmat(1)=a9_den*b8
         Wmat(2)=-a9_den*b7
c: Include separate exponential factor for Wpp and Wps:
         Wmat(5)=-ze
cxx      a9bar_den=fac*beti1(1)*cdexp(-zs)
         a9bar_den=fac*beti1(1)
         Wmat(3)=-a9bar_den*a7
         Wmat(4)=a9bar_den*a8
c: Include separate exponential factor for Wss and Wsp:
         Wmat(6)=-zs
cxx   print *,'prop 1',q1*Wmat(1)+r1sp*Wmat(2),cdexp(ze)*p(1)*
cxx  .   (e1p+e1m*Vmat(1,1))+cdexp(-ze)*u(1)*s2m*Vmat(1,2)
cxx   print *,'prop 2',q2*Wmat(1)+r2sp*Wmat(2),cdexp(ze)*q(1)*
cxx  .   (e2p+e2m*Vmat(1,1))+cdexp(-ze)*v(1)*s1m*Vmat(1,2)
cxx   print *,'prop 3',q3sp*Wmat(1)+r3*Wmat(2),-cdexp(ze)*u(1)*
cxx  .   (e2p+e2m*Vmat(1,1))+cdexp(-ze)*p(1)*s1m*Vmat(1,2)
cxx   print *,'prop 4',q4sp*Wmat(1)+r4*Wmat(2),-cdexp(ze)*v(1)*
cxx  .   (e1p+e1m*Vmat(1,1))+cdexp(-ze)*q(1)*s2m*Vmat(1,2)
cxx   print *,'prop 5',q3sp*Wmat(4)+r3*Wmat(3),cdexp(zs)*p(1)*
cxx  .   (s1p+s1m*Vmat(1,3))-cdexp(-zs)*u(1)*e2m*Vmat(1,4)
cxx   print *,'prop 6',q4sp*Wmat(4)+r4*Wmat(3),cdexp(zs)*q(1)*
cxx  .   (s2p+s2m*Vmat(1,3))-cdexp(-zs)*v(1)*e1m*Vmat(1,4)
cxx   print *,'prop 7',q1*Wmat(4)+r1sp*Wmat(3),cdexp(zs)*u(1)*
cxx  .   (s2p+s2m*Vmat(1,3))+cdexp(-zs)*p(1)*e1m*Vmat(1,4)
cxx   print *,'prop 8',q2*Wmat(4)+r2sp*Wmat(3),cdexp(zs)*v(1)*
cxx  .   (s1p+s1m*Vmat(1,3))+cdexp(-zs)*q(1)*e2m*Vmat(1,4)
      endif
c
c: If first derivatives desired, compute them now:
      do j=2,ndv
         q1p=Rmat(j,1)
         q2p=-gami2(1)*Rmat(j,1) + gami2(j)*q2diff
         q3p=Rmat(j,2)
         q4p=-beti2(1)*Rmat(j,2) - beti2(j)*Rmat(1,2)
         r1p=Rmat(j,4)
         r2p=-gami2(1)*Rmat(j,4) - gami2(j)*Rmat(1,4)
         r3p=Rmat(j,3)
         r4p=-beti2(1)*Rmat(j,3) + beti2(j)*r4diff
         r1spp=spfac*r1p
         r2spp=spfac*r2p
         q3spp=spfac*q3p
         q4spp=spfac*q4p
c
         delta1p=-2.d0*r2p
         delta2p=2.d0*q4p
         delta3p=q4sp*r1spp + q4spp*r1sp - q1*r4p - q1p*r4
         delta4p=q3sp*r2spp + q3spp*r2sp - q2*r3p - q2p*r3
c
         f1p=q4sp*p(j) + q4spp*p(1) + q1*v(j) + q1p*v(1)
         f2p=q4sp*u(j) + q4spp*u(1) - q1*q(j) - q1p*q(1)
         f3p=q3sp*q(j) + q3spp*q(1) + q2*u(j) + q2p*u(1)
         f4p=q3sp*v(j) + q3spp*v(1) - q2*p(j) - q2p*p(1)
         x1p=delta4*f1p + delta4p*f1
         x2p=delta3*f3p + delta3p*f3
         x3p=delta4*f2p + delta4p*f2
         x4p=delta3*f4p + delta3p*f4
c
         del2q1p=delta2*q1p + delta2p*q1
         del2q2p=delta2*q2p + delta2p*q2
         del1q3p=delta1*q3spp + delta1p*q3sp
         del1q4p=delta1*q4spp + delta1p*q4sp
         w1p=del2q2*p(j) + del2q2p*p(1) - del1q3*v(j) - del1q3p*v(1)
         w2p=del2q1*q(j) + del2q1p*q(1) - del1q4*u(j) - del1q4p*u(1)
         w3p=del2q2*u(j) + del2q2p*u(1) + del1q3*q(j) + del1q3p*q(1)
         w4p=del2q1*v(j) + del2q1p*v(1) + del1q4*p(j) + del1q4p*p(1)
         if(isoe .eq. 0) then
            e12gp=gami1(1)*e12(j) + gami1(j)*e12(1)
            e22gp=gami1(1)*e22(j) + gami1(j)*e22(1)
            e1pp=e11(j) + e12gp
            e1mp=e11(j) - e12gp
            e2pp=e21(j) + e22gp
            e2mp=e21(j) - e22gp
         else
            e1pp=e11(j)
            e1mp=e12(j)
            e2pp=e21(j)
            e2mp=e22(j)
         endif
         n1p=w1*e1pp + w1p*e1p - w2*e2pp - w2p*e2p
         n2p=x1*e1pp + x1p*e1p - x2*e2pp - x2p*e2p
         k1p=w1*e1mp + w1p*e1m - w2*e2mp - w2p*e2m
         k2p=x1*e1mp + x1p*e1m - x2*e2mp - x2p*e2m
c
         if(isos .eq. 0) then
            s12bp=beti1(1)*s12(j) + beti1(j)*s12(1)
            s22bp=beti1(1)*s22(j) + beti1(j)*s22(1)
            s1pp=s11(j) + s12bp
            s1mp=s11(j) - s12bp
            s2pp=s21(j) + s22bp
            s2mp=s21(j) - s22bp
         else
            s1pp=s11(j)
            s1mp=s12(j)
            s2pp=s21(j)
            s2mp=s22(j)
         endif
         m1p=w3*s2pp + w3p*s2p - w4*s1pp - w4p*s1p
         m2p=x3*s2pp + x3p*s2p - x4*s1pp - x4p*s1p
         y1p=w3*s2mp + w3p*s2m - w4*s1mp - w4p*s1m
         y2p=x3*s2mp + x3p*s2m - x4*s1mp - x4p*s1m
c
         denp=(y2*k1p + y2p*k1 - y1*k2p - y1p*k2)
         deninvp=-denp*deninv*deninv
         num1p=-(y2*n1p + y2p*n1 - y1*n2p - y1p*n2)
         Vmat(j,1)=num1*deninvp + num1p*deninv
         num2p=(w1*x2p + w1p*x2 - w2*x1p - w2p*x1)
         gamnum2p=gami1(1)*num2p + gami1(j)*num2
         Vmat(j,2)=2.d0*(gamnum2*deninvp + gamnum2p*deninv)
         num3p=(k2*m1p + k2p*m1 - k1*m2p - k1p*m2)
         Vmat(j,3)=num3*deninvp + num3p*deninv
         num4p=(w4*x3p + w4p*x3 - w3*x4p - w3p*x4)
         betnum4p=beti1(1)*num4p + beti1(j)*num4
         Vmat(j,4)=-2.d0*(betnum4*deninvp + betnum4p*deninv)
      enddo
c
      return
      end
      subroutine rp_slay0(Rmat,p,q,u,v,e11,e12,e21,e22,ze,s11,s12,
     .   s21,s22,zs,gami1,beti1,gami2,beti2,isoe,isos,Vmat,ndv,
     .   Wmat,iiw,rhorat)
c
c: Finds the four reflection coefficients and their derivatives
c: referenced to the top of layer 1 given the (p,q,u,v) constants 
c: for the 1-2 layer interface, and the compressional and shear 
c: wave propagator matrices e and s.  Same as rp_slay, except that
c: Rps=Rsp=0.  When Rpp=Rss=0 also, layer 2 is a halfspace, and slightly
c: shorter formulas are used.
c: gam1,gam2 and bet1,bet2 are the vertical comp and shear wavenumbers 
c: at the top of layers 1 & 2.
c
      implicit complex*16 (a-z)
      integer*4 isoe,isos,ndv,j,iiw,iiz
      complex*16 Rmat(3,5),e11(3),e12(3),e21(3),e22(3),s11(3),
     .   s12(3),s21(3),s22(3),p(3),q(3),u(3),v(3),Vmat(3,5),
     .   gami1(3),gami2(3),beti1(3),beti2(3),Wmat(6),zzero
      real*8 rhorat
      data zzero/(0.d0,0.d0)/
c
      if(Rmat(1,1) .eq. zzero .and. Rmat(1,3) .eq. zzero) then
         iiz=1
         q1_q=q(1)
         q1_v=v(1)
         r3_q=q(1)
         r3_v=v(1)
         q2_p=gami2(1)*p(1)
         q2_u=gami2(1)*u(1)
         r4_p=beti2(1)*p(1)
         r4_u=beti2(1)*u(1)
      else
         iiz=0
         q1=1.d0 + Rmat(1,1)
         q2diff=1.d0 - Rmat(1,1)
         q2=gami2(1)*q2diff
         r3=1.d0 + Rmat(1,3)
         r4diff=1.d0 - Rmat(1,3)
         r4=beti2(1)*r4diff
         q1_q=q1*q(1)
         q1_v=q1*v(1)
         r3_q=r3*q(1)
         r3_v=r3*v(1)
         q2_p=q2*p(1)
         q2_u=q2*u(1)
         r4_p=r4*p(1)
         r4_u=r4*u(1)
      endif
c
      L=2.d0*(r4_u*q2_p - q1_q*r3_v)
      if(isoe .eq. 0) then
         e12g=gami1(1)*e12(1)
         e22g=gami1(1)*e22(1)
         e1p=e11(1) + e12g
         e1m=e11(1) - e12g
         e2p=e21(1) + e22g
         e2m=e21(1) - e22g
      else
         e1p=e11(1)
         e1m=e12(1)
         e2p=e21(1)
         e2m=e22(1)
      endif
      m1=q2_p*e1p - q1_q*e2p
      f1=q2_p*e1m - q1_q*e2m
      m2=r4_u*e2p - r3_v*e1p
      f2=r4_u*e2m - r3_v*e1m
c
      if(isos .eq. 0) then
         s12b=beti1(1)*s12(1)
         s22b=beti1(1)*s22(1)
         s1p=s11(1) + s12b
         s1m=s11(1) - s12b
         s2p=s21(1) + s22b
         s2m=s21(1) - s22b
      else
         s1p=s11(1)
         s1m=s12(1)
         s2p=s21(1)
         s2m=s22(1)
      endif
      t1=q2_u*s2p - q1_v*s1p
      c1=q2_u*s2m - q1_v*s1m
      t2=r4_p*s1p - r3_q*s2p
      c2=r4_p*s1m - r3_q*s2m
c
      den=(c2*f1 + c1*f2)
      deninv=1./den
      num1=-(c2*m1 + c1*m2)
      Vmat(1,1)=num1*deninv
      gamL=gami1(1)*L
      Vmat(1,2)=gamL*deninv
      num2=-(t1*f2 + t2*f1)
      Vmat(1,3)=num2*deninv
      betL=beti1(1)*L
      Vmat(1,4)=-betL*deninv
c
c: vps(1:2) and vsp(1:2) are to be multiplied by exp(zexp_sp):
      Vmat(1,5)=-(ze+zs)
c
      if(iiw .eq. 1) then
c: Take from rp_prop, but set r1=r2=q3=q4=0:
         a7=e2m*q1_q - e1m*q2_p
         a8=e1m*r3_v - e2m*r4_u
         b7=s1m*q1_v - s2m*q2_u
         b8=s1m*r4_p - s2m*r3_q
         wden=a7*b8 - b7*a8
cxx      fac=-2.d0*(p(1)*q(1) + u(1)*v(1))/wden
         fac=-2.d0*rhorat/wden
         a9_den=fac*gami1(1)
         Wmat(1)=a9_den*b8
         Wmat(2)=-a9_den*b7
         Wmat(5)=-ze
         a9bar_den=fac*beti1(1)
         Wmat(3)=-a9bar_den*a7
         Wmat(4)=a9bar_den*a8
         Wmat(6)=-zs
cxx   print *,'Wmat 1',q1*Wmat(1),cdexp(ze)*p(1)*
cxx  .   (e1p+e1m*Vmat(1,1))+cdexp(-ze)*u(1)*s2m*Vmat(1,2)
cxx   print *,'Wmat 2',q2*Wmat(1),cdexp(ze)*q(1)*
cxx  .   (e2p+e2m*Vmat(1,1))+cdexp(-ze)*v(1)*s1m*Vmat(1,2)
cxx   print *,'Wmat 3',r3*Wmat(2),-cdexp(ze)*u(1)*(e2p+e2m*
cxx  .   Vmat(1,1))+cdexp(-ze)*p(1)*s1m*Vmat(1,2)
cxx   print *,'Wmat 4',r4*Wmat(2),-cdexp(ze)*v(1)*(e1p+e1m*
cxx  .   Vmat(1,1))+cdexp(-ze)*q(1)*s2m*Vmat(1,2)
cxx   print *,'Wmat 5 ',q1*Wmat(4),cdexp(zs)*u(1)*(s2p + s2m*
cxx  .   Vmat(1,3)) + cdexp(-zs)*p(1)*e1m*Vmat(1,4)
cxx   print *,'Wmat 6 ',q2*Wmat(4),cdexp(zs)*v(1)*(s1p + s1m*
cxx  .   Vmat(1,3)) + cdexp(-zs)*q(1)*e2m*Vmat(1,4)
cxx   print *,'Wmat 7 ',r3*Wmat(3),cdexp(zs)*p(1)*(s1p + s1m*
cxx  .   Vmat(1,3)) - cdexp(-zs)*u(1)*e2m*Vmat(1,4)
cxx   print *,'Wmat 8 ',r4*Wmat(3),cdexp(zs)*q(1)*(s2p + s2m*
cxx  .   Vmat(1,3)) - cdexp(-zs)*v(1)*e1m*Vmat(1,4)
      endif
c
      do j=2,ndv
         if(iiz .eq. 1) then
            q1_qp=q(j)
            q1_vp=v(j)
            r3_qp=q(j)
            r3_vp=v(j)
            q2_pp=gami2(1)*p(j) + gami2(j)*p(1)
            q2_up=gami2(1)*u(j) + gami2(j)*u(1)
            r4_pp=beti2(1)*p(j) + beti2(j)*p(1)
            r4_up=beti2(1)*u(j) + beti2(j)*u(1)
         else
            q1p=Rmat(j,1)
            q2diffp=-Rmat(j,1)
            q2p=gami2(1)*q2diffp + gami2(j)*q2diff
            r3p=Rmat(j,3)
            r4diffp=-Rmat(j,3)
            r4p=beti2(1)*r4diffp + beti2(j)*r4diff
            q1_qp=q1*q(j) + q1p*q(1)
            q1_vp=q1*v(j) + q1p*v(1)
            r3_qp=r3*q(j) + r3p*q(1)
            r3_vp=r3*v(j) + r3p*v(1)
            q2_pp=q2*p(j) + q2p*p(1)
            q2_up=q2*u(j) + q2p*u(1)
            r4_pp=r4*p(j) + r4p*p(1)
            r4_up=r4*u(j) + r4p*u(1)
         endif
c
         Lp=2.d0*(r4_u*q2_pp + r4_up*q2_p - q1_q*r3_vp - q1_qp*r3_v)
         if(isoe .eq. 0) then
            e12gp=gami1(1)*e12(j) + gami1(j)*e12(1)
            e22gp=gami1(1)*e22(j) + gami1(j)*e22(1)
            e1pp=e11(j) + e12gp
            e1mp=e11(j) - e12gp
            e2pp=e21(j) + e22gp
            e2mp=e21(j) - e22gp
         else
            e1pp=e11(j)
            e1mp=e12(j)
            e2pp=e21(j)
            e2mp=e22(j)
         endif
         m1p=q2_p*e1pp + q2_pp*e1p - q1_q*e2pp - q1_qp*e2p
         m2p=r4_u*e2pp + r4_up*e2p - r3_v*e1pp - r3_vp*e1p
         f1p=q2_p*e1mp + q2_pp*e1m - q1_q*e2mp - q1_qp*e2m
         f2p=r4_u*e2mp + r4_up*e2m - r3_v*e1mp - r3_vp*e1m
c
         if(isos .eq. 0) then
            bets12p=beti1(1)*s12(j) + beti1(j)*s12(1)
            bets22p=beti1(1)*s22(j) + beti1(j)*s22(1)
            s1pp=s11(j) + bets12p
            s1mp=s11(j) - bets12p
            s2pp=s21(j) + bets22p
            s2mp=s21(j) - bets22p
         else
            s1pp=s11(j)
            s1mp=s12(j)
            s2pp=s21(j)
            s2mp=s22(j)
         endif
         t1p=q2_u*s2pp + q2_up*s2p - q1_v*s1pp - q1_vp*s1p
         t2p=r4_p*s1pp + r4_pp*s1p - r3_q*s2pp - r3_qp*s2p
         c1p=q2_u*s2mp + q2_up*s2m - q1_v*s1mp - q1_vp*s1m
         c2p=r4_p*s1mp + r4_pp*s1m - r3_q*s2mp - r3_qp*s2m
c
         denp=(c2*f1p + c2p*f1 + c1*f2p + c1p*f2)
         deninvp=-denp*deninv*deninv
         num1p=-(c2*m1p + c2p*m1 + c1*m2p + c1p*m2)
         Vmat(j,1)=num1*deninvp + num1p*deninv
         gamLp=gami1(1)*Lp + gami1(j)*L
         Vmat(j,2)=gamL*deninvp + gamLp*deninv
         num2p=-(t1*f2p + t1p*f2 + t2*f1p + t2p*f1)
         Vmat(j,3)=num2*deninvp + num2p*deninv
         betLp=beti1(1)*Lp + beti1(j)*L
         Vmat(j,4)=-betL*deninvp - betLp*deninv
      enddo
c
      return
      end
      subroutine mode_cmplx_field(tfz,phiz,dpsiz,expz_gbs,jzs,jf,nf)
c
c: Computes field given mode eigenvalues and mode function values.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      integer*4 jf,jzs,jm,jr,jd,jsr,kk,kk0,jsrc,jrec,msrc,
     .          nrmax,nf,jj0,jj1,jj2
      complex*16 knx,iknx,iknr,sqkn,cfac2,phi_src,H0,h_arg,
     .   cfac,exp_sqk
      complex*8 phiz(nzsr,nmode),dpsiz(nzsr,nmode),
     .          tfz(nf*nzsr*nrec)
      real*8 expz_gbs(nzsr,nmode),rsig_max,hmsq,xn_b,xn_beam,
     .   magsq,rng_im
c: cfac = sqrt(2*pi)*exp(i*pi/4) (see p. 7, ORCA II):
c      data cfac/(1.77245385090552d0,1.77245385090552d0)/
c

      cfac=(1.77245385090552d0,1.77245385090552d0)
      msrc=mzsrc(jzs)
      if(iigbs .ne. 0) then
         rng_im=-b_gbs(jzs)*cos(th_gbs(jzs)*pie/180.d0)
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .              +dtiltvp(jd),rng_im)+dtiltvp(jd))
            end do
         enddo
      else
         rng_im=0.d0
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/dsqrt(range(jr)+dtiltvp(jd))
            end do
         enddo
      endif
c
      hmsq=25.d0
c: Initialize complex propagation loss array:
c      do jsrc=1,nsrc
c         do jrec=1,nrec
c            tfz(jrec,jsrc)=cmplx(0.e0,0.e0)
c         enddo
c      enddo
c
      xn_b=(w/dreal(cp_sr(msrc)))*b_gbs(jzs)
      do jm=1,nmode
         knx=kn(jm)
         iknx=dcmplx(-dimag(knx),dreal(knx))
         sqkn=cdsqrt(knx)
c: rsig_max is maximum range at which this mode is significant (50 dB down 
c: from strongest mode):
         rsig_max=kim_fac/dmax1(1.d-20,dimag(knx)-kim_min)
c: Obtain p-wave mode excitation at source depth:
         phi_src=phiz(msrc,jm)/rho_sr(msrc)
         xn_beam=xn_b - expz_gbs(msrc,jm)
c: Obtain p-wave mode excitations at receiver depths:
         do jrec=1,nrec
            jsr=mzrec(jrec)
c: rho included in mode_fun:
            phisr(jsr)=phiz(jsr,jm)
            if(kksh(jsr) .eq. 1) then
c: Add shear wave potential contributions to recs in shear layers (see p. 117):
               phisr(jsr)=phisr(jsr) - 2.d0*ksm2_sr(jsr)*
     .            (knx*knx*phisr(jsr) + iknx*dpsiz(jsr,jm))
            endif
         enddo
c: Find range beyond which mode is insignificant (ranges have been sorted):
         nrmax=0
         call hunt(range,nrng,rsig_max,nrmax)
         do jr=1,nrmax
cpln            h_arg=knx*dcmplx(range(jr),rng_im)
            kk0=krec_jr(jr)
            do kk=kk0+1,kk0+nrec_jr(jr)
               jsrc=jrec_jr(1,kk)
               jrec=jrec_jr(2,kk)
               h_arg=knx*dcmplx(range(jr)+dtiltvp(jrec),
     .               rng_im+dtiltvp(jrec))
               if(magsq(h_arg) .gt. hmsq) then
czs: Include normalization by exp(-xn_beam) here:
                  iknr=dcmplx(-dimag(h_arg),dreal(h_arg))
                  exp_sqk=cdexp(dcmplx(dreal(iknr)-xn_beam,dimag(iknr)))
     .                 /sqkn
                  cfac2=sq2pir(jr,jrec)*phi_src*exp_sqk
               else
                  call cdhankel(h_arg,1.d-6,H0)
cc             cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
                  cfac2=dcmplx(-pie*dimag(H0),pie*dreal(H0))*phi_src*
     .                 dexp(-xn_beam)
               endif
               jj0=(jsrc-1)*nf*nrec
               jj1=(jrec-1)*nf
               jsr=mzrec(jrec)
               jj2=jj0+jj1
               tfz(jf+jj2)=tf(jf+jj2) + cfac2*phisr(jsr)
cpln               tfz(jsrc,jrec)=tfz(jsrc,jrec)+cfac2*phisr(jsr)
            enddo
         enddo
99       continue
      enddo
c
c               write(6,*)'jrec,jsrc,jf: ',jrec,jsrc,jf
c               write(6,*)'jj0,jj1,jj2: ',jj0,jj1,jj2
c               tlz(1,jrec,jsrc)=tlz(1,jrec,jsrc) + cfac2*phisr(jsr)
c               write(6,*)tf(jf+jj2)
c               pause
      return
      end


      subroutine scairy2(z1,ai1,aip1,bi1,bip1,zeta1,det1,z2,ai2,aip2,
     .   bi2,bip2,zeta2,det2,aisoln,iibad)
c
c: This subroutine computes the two numerically stable solutions to the
c: Airy differential equation for the two arguments z1 and z2.
c: The derivatives with respect to z1 and z2 are also computed.
c: For abs(arg(z)) < pi/3, the two solutions are Ai and Bi.
c: For pi/3 <= arg(x) <= pi, the two sol'n are Ai(z) and Ai(z*exp(-2pi/3)).
c: For -pi <= arg(x) <= pi/3, the two sol'n are Ai(z) and Ai(z*exp(2pi/3)).
c: Uses the Airy function routines from SAFARI, where one bug has been 
c: fixed in CLAIRY (see EKW note).  Note also that the derivatives of
c: Ai(z*exp(2pi/3)) are with respect to z, not the total argument.
c: Special care is taken that when z1 and z2 straddle the negative real
c: axis, the same pair of independent solutions are computed.
c
      IMPLICIT COMPLEX*16 (A,B,D,Z)
      include 'scairy_com'
      integer*4 aisoln,iibad
c
      iibad=0
      call clairy(z1,1,ai1,bi1,aip1,bip1,zeta1)
c: EKW changed: For imag(z)>0, we want (d/dz)(Ai(z*eim23)), whereas 
c: CLAIRY is giving BIP = (d/d(z*eim23))(Ai(z*eim23)), 
c: so multiply BIP by eim23:
      IF (dimag(z1).GE.0.D0)  THEN
         bip1=eim23*bip1
         det1=det_pos
         aisoln=1
      ELSE
         bip1=ei23*bip1
         det1=det_neg
         aisoln=-1
      END IF
c
      call clairy(z2,1,ai2,bi2,aip2,bip2,zeta2)
c: EKW changed: For imag(z)>0, we want (d/dz)(Ai(z*eim23)), whereas 
c: CLAIRY is giving BIP = (d/d(z*eim23))(Ai(z*eim23)), 
c: so multiply BIP  by eim23:
      IF (dimag(z2).GE.0.D0)  THEN
         bip2=eim23*bip2
         det2=det_pos
      ELSE
         bip2=ei23*bip2
         det2=det_neg
      END IF
c
c: Check if z1 and z2 straddle the real axis:
      if(dimag(z1) .ge. 0. .and. dimag(z2) .lt. 0.) then
         if(dreal(zeta2) .gt. dreal(zeta1)) then
c: z2 on neg side, so change Ai(z2*ei23) to Ai(z2*eim23):
            call ai_flip(zeta2,ai2,aip2,bi2,bip2,ei23,eim23,
     .         det2,det_pos,iibad)
         else
c: z1 on pos side, so change Ai(z1*eim23) to Ai(z1*ei23):
            call ai_flip(zeta1,ai1,aip1,bi1,bip1,eim23,ei23,
     .         det1,det_neg,iibad)
            aisoln=-1
         endif
      elseif(dimag(z2) .ge. 0. .and. dimag(z1) .lt. 0.) then
         if(dreal(zeta1) .gt. dreal(zeta2)) then
c: z1 on neg side, so change Ai(z1*ei23) to Ai(z1*eim23):
            call ai_flip(zeta1,ai1,aip1,bi1,bip1,ei23,eim23,
     .         det1,det_pos,iibad)
            aisoln=1
         else
c: z2 on pos side, so change Ai(z2*eim23) to Ai(z2*ei23):
            call ai_flip(zeta2,ai2,aip2,bi2,bip2,eim23,ei23,
     .         det2,det_neg,iibad)
         endif
      endif
c
      RETURN
      END
ccc
      subroutine ai_flip(zeta2,ai2,aip2,bi2,bip2,ei23,eim23,
     .   det2,det_pos,iibad)
c
      implicit none
      integer*4 iibad
      complex*16 zeta2,ai2,aip2,bi2,bip2,ei23,eim23,det2,det_pos,zfac,
     .   diff
c
      if(dreal(zeta2) .gt. 0.d0) then
         zfac=cdexp(-2.d0*zeta2)
         bi2=-(ei23*ai2*zfac + eim23*bi2)
         bip2=-(ei23*aip2*zfac + eim23*bip2)
      elseif(dreal(zeta2) .lt. -50.d0) then
c: Need to split into two layers in ai_strad later:
         iibad=1
      else
cxx      print *,'ai_flip for re(zeta)<0: ',zeta2,ai2,aip2,bi2,bip2
         zfac=cdexp(-zeta2)
         bi2=-(ei23*ai2*zfac + eim23*bi2/zfac)
         bip2=-(ei23*aip2*zfac + eim23*bip2/zfac)
         ai2=ai2*zfac
         aip2=aip2*zfac
         det2=ai2*bip2 - aip2*bi2
         diff=det_pos - det2
         if(dreal(diff)*dreal(diff)+dimag(diff)*dimag(diff) .gt. 
     .      1.e-12) iibad=1
         zeta2=(0.d0,0.d0)
      endif
      det2=det_pos
cxx   call wronsk_ch(ai2,aip2,bi2,bip2,det2)
c
      return
      end
ccc
      subroutine wronsk_ch(ai2,aip2,bi2,bip2,det2)
c
      implicit none
      complex*16 ai2,aip2,bi2,bip2,det2,det
c
      det=ai2*bip2 - aip2*bi2
      if(abs(det-det2) .gt. 1.e-4) then
         print *,'W = ',det,det2
      endif
c
      return
      end
ccc
      subroutine scairy3(nz,z,z0,ai,aip,bi,bip,zeta,det,aisoln)
c
c: Identical to scairy2, except computes a set of nz airy functions
c: at z(1:nz), assuring that all are of the same solution.
c
      implicit none
      integer*4 nz,jz,aisoln,iibad
      complex*16 z(nz),z0,ai(nz),aip(nz),bi(nz),bip(nz),zeta(nz),det
      include 'scairy_com'
c
cc    if(dimag(z0).GE.0.D0) then
      if(aisoln .eq. 1) then
         do jz=1,nz
            call clairy(z(jz),1,ai(jz),bi(jz),aip(jz),bip(jz),zeta(jz))
            bip(jz)=eim23*bip(jz)
c: z on neg side, so change Ai(z*ei23) to Ai(z*eim23):
            if(dimag(z(jz)) .lt. 0.d0) then
               call ai_flip(zeta(jz),ai(jz),aip(jz),bi(jz),bip(jz),
     .            ei23,eim23,det,det_pos,iibad)
            endif
         enddo
         det=det_pos
      else
         do jz=1,nz
            call clairy(z(jz),1,ai(jz),bi(jz),aip(jz),bip(jz),zeta(jz))
            bip(jz)=ei23*bip(jz)
c: z on pos side, so change Ai(z*eim23) to Ai(z*ei23):
            if(dimag(z(jz)) .gt. 0.d0) then
               call ai_flip(zeta(jz),ai(jz),aip(jz),bi(jz),bip(jz),
     .            eim23,ei23,det,det_neg,iibad)
            endif
         enddo
         det=det_neg
      endif
c
      return
      end
ccc
      SUBROUTINE CLAIRY (Z,IRET, AI,BI,AIP,BIP,ZETA)
      IMPLICIT COMPLEX*16 (A,B,E,X,Z)
      PARAMETER (ONEOOF=1D0/1.5D0)

C     COMPLEX AIRY FUNCTION SUBROUTINE
C     Z; THE COMPLEX ARGUMENT
C     IRET; RETURN CONTROL PARAMETER FOR SCALING
C           < 0 RETURN SCALED AIRY FUNCTIONS OF THE FIRST AND SECOND TYPES.
C           = 0 RETURN UNSCALED AIRY FUNCTIONS OF THE FIRST AND SECOND TYPES.
C           > 0 RETURN SCALED RESULTS OF THE NUMERICALLY STABLE INDEPENDENT 
C               TYPE.  THESE ARE AI(Z) AND AI(Z*EXP([+,-]2IPI/3)).
C     AI AND AIP; THE COMPLEX AIRY FUNCTION OF THE FIRST TYPE AND ITS
C                   DERIVATIVE DEPENDING ON IRET.
C               IRET < 0, SCALED AIRY FUNCTION OF THE FIRST TYPE.  AI AND AIP
C                   ARE TO BE DIVIDED BY EXP(ZETA) TO OBTAIN THE TRUE AIRY 
C                   FUNCTIONS.
C               IRET = 0, UNSCALED AIRY FUNCTION OF THE FIRST TYPE
C               IRET > 0, SCALED AIRY FUNCTION OF THE FIRST TYPE.  AI AND AIP
C                   ARE TO BE DIVIDED BY EXP(ZETA) TO OBTAIN THE TRUE AIRY 
C                   FUNCTIONS.
C     BI AND BIP; THE COMPLEX AIRY FUNCTION OF THE SECOND TYPE AND ITS
C               DERIVATIVE DEPENDING ON IRET AND THE ANGLE OF Z.
C               IRET < 0, SCALED AIRY FUNCTION OF THE SECOND TYPE.  BI AND BIP 
C                   ARE TO BE MULTIPLIED BY EXP(ZETA), FOR ARG(Z) < 60 AND 
C                   DIVIDED BY EXP(ZETA) FOR ARG(Z) > 60.
C               IRET = 0, UNSCALED AIRY FUNCTION OF THE SECOND TYPE
C               IRET > 0, SCALED AIRY FUNCTION OF OF THE FIRST TYPE FOR THE
C     ARGUMENT Z*EXP(-2IPI/3) WHEN ARG(Z) > 0 AND FOR
C     Z*EXP(2IPI/3) WHEN ARG(Z) < 0. THE RESULTS ARE TO
C     BE MULTIPLIED BY EXP(ZETA). 
C      ZETA; THE EXPONTIAL SCALING FACTOR
C
C     THIS ROUTINE USES THE ALGORITHMS DESCRIBED IN "AN ALGORITHM FOR 
C     THE EVALUATION OF THE COMPLEX AIRY FUNCTIONS,"  Z. SCHULTEN, 
C     D.G.M. ANDERSON, AND R.G. GORDEN, IN JOUR. COMP. PHYSICS 31, 
C     PP. 60-75, (1979).
C     A BETTER DESCRIPTION OF THE SERIES EXPANSION IS AVAILABLE IN
C     HANDBOOK OF MATHEMATICAL FUNCTIONS, APPLIED MATHEMATICS SERIES #55,
C     M. ABRAMOWITZ AND I. A. STEGAN, P. 446 (1964).

C      PARAMETER (TPIO3=2.094395102393195D0, I=(0.0D0,1.0D0))
      REAL*8 TPIO3
      REAL*8 ZPR(2)
      include 'scairy_com'
      EQUIVALENCE (ZPR(1),ZPP)

      DATA TPIO3 /2.094395102393195D0/
C   * STATEMENT FUNCTION DEFINITION OF THE PHSE *

C      THETA (Z) = DATAN2 (DIMAG (Z), REAL (Z))

C   * Z OR CONJUGATE Z IN POSITIVE IMAGINARY PLANE *

      ZPP=Z

      IF (ZPR(2) .LT. 0.0D0)  THEN
          ZP = CONJG (Z)
          ICONJG = 1
      ELSE
          ZP = Z
          ICONJG = 0
      END IF


C   *********************************************************
C   *               *
C   *   SERIES SOLUTIONS FOR AI & BI WITHIN 3.5 OF ORIGIN   *
C   *               *
C   *********************************************************

      IF (ABS (Z) .LE. 3.5D0)  THEN
          CALL ABSERIES2 (Z, AI,AIP,BI,BIP)
          ZETA = (0.0D0, 0.0D0)
          IF (IRET.NE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,0,IRET)


C   ******************************************************************
C   *    *
C   *   ANGLES LESS THAN 120 DEGS. - AI(Z) AND AI(Z*EXP(-2IPI/3))    *
C   *   ARE CALCULATED FOR Z IN THE UPPER HALF PLANE, CONJGATE OF    *
C   *   Z IF IN THE LOWER.  SERIES OR GAUSSIAN SOLUTION AS NEEDED.   *
C   *    *
C   ******************************************************************

c: EKW: This was a bug originally because ZPR is equivalenced to original
c: input value of Z, not ZP, which is always in upper half plane:
cxx   ELSE IF (ATAN2(ZPR(2),ZPR(1)) .LE. TPIO3)  THEN
cxx   ELSE IF (abs(ATAN2(ZPR(2),ZPR(1))) .LE. TPIO3)  THEN
      ELSE IF (dimag(zp) .gt. -1.73205080756888d0*dreal(zp)) then
          IF (ABS(ZP).GT.7.9111D0)  THEN          
              CALL ABGAUSS2 (ZP, AI, AIP, BI, BIP, ZETA)
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  BI  = CONJG (BI)
                  AIP = CONJG (AIP)
                  BIP = CONJG (BIP)
              END IF
          ELSE
              IF (ABS (ZP-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
                  CALL ASERIES2 (ZP, AI,AIP, ZETA)
              ELSE
                  CALL AGAUSS2 (ZP, AI, AIP, ZETA)
              END IF
              ZP2 = CONJG (ZP) * ei23
              IF (ABS (ZP2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
                  CALL ASERIES2 (ZP2, BI, BIP, ZETA2)
              ELSE
                  CALL AGAUSS2 (ZP2, BI, BIP, ZETA2)
              END IF
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  AIP = CONJG (AIP)
              ELSE
                  BI  = CONJG (BI)
                  BIP = CONJG (BIP)
              END IF
          END IF
          IF (IRET.LE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,1,IRET)


C   **********************************************************
C   *                *
C   *   CONNECTION REGION FOR ANGLE GREATER THAN 120 DEGS.   *
C   *   CALCULATE AI(Z*EXP(-2IPI/3)) AND AI(Z*EXP(2IPI/3))   *
C   *   SOLUTIONS, USE CONNECTION FORMULAS FOR RESULTS.      *
C   *                *
C   **********************************************************

      ELSE 
c: z1=z*exp(-i2pi/3); z2=conjg[z*exp(i2pi/3] (to make z2 in upper half plane;
c: take conjugate of ai,aip,zeta afterwards):
          Z1 = Z * eim23    
          Z2 = CONJG (Z) * eim23                
          
          IF (ABS (Z1-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z1, AI1, AIP1, ZETA1)
          ELSE
              CALL AGAUSS2 (Z1, AI1, AIP1, ZETA1)
          END IF
          IF (ABS (Z2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z2, AI2, AIP2, ZETA2)
          ELSE
              CALL AGAUSS2 (Z2, AI2, AIP2, ZETA2)
          END IF
          AI2  = CONJG (AI2)  
          AIP2 = CONJG (AIP2)
          ZETA2 = CONJG (ZETA2)

c: Use connection formula: 
c: Ai(z)=exp(ipi/3)*Ai(z*exp(-i2pi/3)) + exp(-ipi/3)*z*exp(i2pi/3)
          IF (IRET.EQ.0)  THEN
              E1 = EXP (-ZETA1)
              E2 = EXP (-ZETA2)
              AI  = ei13 * AI1*E1  + eim13 * AI2*E2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
              BI  = ei16 * AI2*E2  + eim16 * AI1*E1
              BIP = ei56 * AIP2*E2 + eim56 * AIP1*E1
              ZETA = (0.0D0, 0.0D0)

          ELSE IF (IRET.LT.0)  THEN               
              ZETA  = Z * SQRT(Z) * ONEOOF
              E1  = EXP (ZETA-ZETA1)
              E2  = EXP (ZETA-ZETA2)
              AI  = ei13 * AI1*E1  + eim13 * AI2*E2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
              BI  = eim16 * AI1*E1  + ei16 * AI2*E2
              BIP = eim56 * AIP1*E1 + ei56 * AIP2*E2

          ELSE IF (IRET.GT.0)  THEN               
              ZETA  = Z * SQRT(Z) * ONEOOF
              E1  = EXP (ZETA-ZETA1)
              E2  = EXP (ZETA-ZETA2)
              AI  = ei13 * AI1*E1  + eim13 * AI2*E2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
cekw          IF (ZPR(2).GE.0.0D0)  THEN
              IF (ICONJG .eq. 0) then
                  E1  = EXP (-ZETA-ZETA1)
                  BI  = AI1  * E1
                  BIP = AIP1 * E1
              ELSE
                  E2  = EXP (-ZETA-ZETA2)
                  BI  = AI2  * E2
                  BIP = AIP2 * E2
              END IF
          END IF
      END IF

      RETURN
      END
      SUBROUTINE ABSERIES2 (Z, AI,AIP,BI,BIP)
      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      REAL*8 C1,C2,SQR3,FAC1,FAC2,FAC3
      DATA C1 /0.355028053887817D0/, C2 /0.258819403792807D0/,
     1          SQR3 /1.732050807568877D0/

      FC = (1.0D0, 0.0D0)
      GC = Z
      Z2 = Z * Z
      AI = C1*FC - C2*GC
      BI = C1*FC + C2*GC
      AIP = -C2
      BIP =  C2

      DO 10 K = 3,600,3
      FAC1=1.0D0/(K-1)
      FAC2=1.0D0/K
      FAC3=1.0D0/(K+1)
      FC = FC * Z2 * FAC1    
      GC = GC * Z2 * FAC2
      AIP = AIP + (C1*FC - C2*GC)
      BIP = BIP + (C1*FC + C2*GC)
      FC = FC * Z * FAC2       
      GC = GC * Z * FAC3
      AI = AI + (C1*FC - C2*GC)
      BI = BI + (C1*FC + C2*GC)
      IF (ABS(FC).LT.1.D-17*ABS(AI) .AND. 
     &    ABS(GC).LT.1.D-17*ABS(AI))  GO TO 12
   10 CONTINUE
      WRITE (6,*) ' SERIES APPROXIMATION FAILED TO CONVERGE'

   12 BI = SQR3 * BI
      BIP = SQR3 * BIP
      RETURN
      END
      SUBROUTINE ABGAUSS2 (Z, AI,AIP,BI,BIP,ZETA)

      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      PARAMETER (ONEOOF=1D0/1.5D0)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (ONEOPI=1.0D0/PI) 
      COMPLEX*16 SQRTZ,SQRPIZ,CFAC,CFAC2,ONEOZ
      REAL*8 ZPR(2)
      EQUIVALENCE (ZPR(1),ZPP)
      REAL*8 X4(4), W4(4), X6(6), W6(6)
      include 'scairy_com'

      DATA X4/3.9198329554455091D0,  1.6915619004823504D0,
     1        5.0275532467263918D-1, 1.9247060562015692D-2/,
     2     W4/4.7763903057577263D-5, 4.9914306432910959D-3,
     3        8.6169846993840312D-2, 9.0879095845981102D-1/,
     4     X6/7.1620871339075440D0,  4.2311006706187214D0,
     5        2.3361772245064852D0,  1.0856431202004936D0,
     6        3.3391648924379639D-1, 1.3115888501576988D-2/,
     7     W6/4.9954496303045166D-8, 1.8066384626280827D-5,
     8        9.5530673977919037D-4, 1.5715675321710695D-2,
     9        1.1588902608004444D-1, 8.6742187551934309D-1/

      ZPP=Z

      SQRTZ = SQRT (Z)
      ZETA   = Z * SQRTZ * ONEOOF
      SQRPIZ = SQRT (PI * SQRTZ)
      AI = (0.0D0, 0.0D0)
      BI = (0.0D0, 0.0D0)
      AIP = (0.0D0, 0.0D0)
      BIP = (0.0D0, 0.0D0)

      IF (ABS(Z).GT.25.D0)  THEN
          DO 10 I = 1,4
          CFAC=1.0D0/(ZETA+X4(I))
          CFAC2=1.0D0/(ZETA-X4(I))
          AI = AI + W4(I) * ZETA * CFAC
          AIP = AIP + W4(I) * X4(I) * SQRTZ * CFAC**2
          BIP = BIP + W4(I) * X4(I) * SQRTZ * CFAC2**2
   10     BI = BI + W4(I) * ZETA * CFAC2
      ELSE
          DO 20 I = 1,6
          CFAC=1.0D0/(ZETA+X6(I))
          CFAC2=1.0D0/(ZETA-X6(I))
          AI = AI + W6(I) * ZETA * CFAC
          AIP = AIP + W6(I) * X6(I) * SQRTZ * CFAC**2
          BIP = BIP + W6(I) * X6(I) * SQRTZ * CFAC2**2
   20     BI = BI + W6(I) * ZETA * CFAC2
      END IF
      CFAC=1.0D0/SQRPIZ
      ONEOZ=1.0D0/Z
      AIP = - AI *0.125D0 * CFAC * ONEOZ
     1      - AI * SQRPIZ * 0.5D0 * ONEOPI
     1      + AIP * 0.5D0 * CFAC
      BIP = - BI * 0.25D0 * CFAC * ONEOZ
     1      + BI * SQRPIZ * ONEOPI
     1      - BIP * CFAC
      AI = AI * 0.5D0 * CFAC
      BI = BI * CFAC
      IF (ZPR(2).GE.0)  THEN                
          BI  = ei16 * BI * 0.5D0
          BIP = ei56 * BIP * 0.5D0
      ELSE
          BI  = eim16 * BI * 0.5D0
          BIP = eim56 * BI * 0.5D0
      END IF
      RETURN
      END
      SUBROUTINE ASERIES2 (Z, AI, AIP, ZETA)
      IMPLICIT COMPLEX*16 (A-H,O-Z) 
      PARAMETER (ONEOOF=1D0/1.5D0)
      REAL*8 C1,C2,FAC1,FAC2,FAC3
      DATA C1/0.355028053887817D0/, C2/0.258819403792807D0/

      FC = (1.0D0, 0.0D0)
      GC = Z
      Z2 = Z * Z
      AI = C1*FC - C2*GC
      AIP= -C2

      DO 10 K = 3,600,3
      FAC1=1.0D0/(K-1)
      FAC2=1.0D0/K
      FAC3=1.0D0/(K+1)
      FC = FC * Z2 * FAC1                   
      GC = GC * Z2 * FAC2
      AIP = AIP + (C1*FC - C2*GC)
      FC = FC * Z * FAC2    
      GC = GC * Z * FAC3
      AI = AI + (C1*FC - C2*GC)
      IF (ABS(FC).LT.1.D-17*ABS(AI) .AND. 
     &    ABS(GC).LT.1.D-17*ABS(AI))  GO TO 12
   10 CONTINUE
      WRITE (6,*) ' SERIES APPROXIMATION FAILED TO CONVERGE'

   12 ZETA = Z * SQRT(Z) * ONEOOF
      EXPZETA = EXP (ZETA)
      AI = AI * EXPZETA
      AIP = AIP * EXPZETA
      RETURN
      END
ccc
      SUBROUTINE AGAUSS2 (Z,AI,AIP,ZETA)
      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      PARAMETER (ONEOOF=1D0/1.5D0)
      PARAMETER (PI=3.1415926535897932D0)
      PARAMETER (ONEOPI=1.0D0/PI) 
      REAL*8 X4(4), W4(4), X6(6), W6(6)
     
      DATA X4/3.9198329554455091D0,  1.6915619004823504D0,
     1        5.0275532467263918D-1, 1.9247060562015692D-2/,
     2     W4/4.7763903057577263D-5, 4.9914306432910959D-3,
     3        8.6169846993840312D-2, 9.0879095845981102D-1/,
     4     X6/7.1620871339075440D0,  4.2311006706187214D0,
     5        2.3361772245064852D0,  1.0856431202004936D0,
     6        3.3391648924379639D-1, 1.3115888501576988D-2/,
     7     W6/4.9954496303045166D-8, 1.8066384626280827D-5,
     8        9.5530673977919037D-4, 1.5715675321710695D-2,
     9        1.1588902608004444D-1, 8.6742187551934309D-1/

      SQRTZ = SQRT (Z)
      ZETA = Z * SQRTZ * ONEOOF
      SQRPIZ = SQRT (PI * SQRTZ)
      AI  = (0.0D0, 0.0D0)
      AIP = (0.0D0, 0.0D0)

      IF (ABS(Z).GT.25.D0)  THEN
          DO 10 I = 1,4
          CFAC=1.0D0/(ZETA+X4(I))
          AI  = AI + W4(I) * ZETA * CFAC
   10     AIP = AIP + W4(I) * X4(I) * SQRTZ * CFAC**2
      ELSE
          DO 20 I = 1,6
          CFAC=1.0D0/(ZETA+X6(I))
          AI  = AI + W6(I) * ZETA * CFAC
   20     AIP = AIP + W6(I) * X6(I) * SQRTZ * CFAC**2
      END IF
      CFAC=1.0D0/SQRPIZ
      ONEOZ=1.0D0/Z
      AIP = - AI * 0.125D0 * CFAC * ONEOZ
     1      - AI * SQRPIZ * 0.5D0 * ONEOPI
     1      + AIP * 0.5D0 * CFAC
      AI  = AI * 0.5D0 * CFAC
      RETURN
      END
      SUBROUTINE FORMAIRY (AI, BI, AIP, BIP, Z, ZETA, IRETIN, IRETOUT)
      IMPLICIT COMPLEX*16 (A-E,X-Z)

C     THIS ROUTINE IS DESIGNED TO CONVERT FROM ONE TYPE OF AIRY FUNCTIONS
C     TO ANOTHER.  THE THREE REPRESENTATIONS CORRESPOND TO 
C     IRET = 0, AI AND BI UNSCALED
C     IRET < 0, AI AND BI SCALED, AI ALWAYS BY DIVIDING BY EXP(ZETA), BI BY
C               MULTIPLYING BY EXP(ZETA) FOR REAL(ZETA) > 0 (OR ARG(Z) < 60) AND
C               BY DIVIDING BY EXP(ZETA) FOR REAL(ZETA) < 0 (OR ARG(Z) > 60).
C     IRET > 0  AI(Z) AND AI(Z*EXP(2IPI/3)), THE FIRST BY DIVIDING BY EXP(ZETA)
C               AND THE SECOND BY MULTIPLYING BY EXP(ZETA).
C     FOR A TOTAL OF NINE CASES, THREE OF WHICH ARE TRIVAL.

      REAL*8 ZPR(2)
      EQUIVALENCE (ZPR(1),ZPP)
      COMPLEX*16 COMPI
      include 'scairy_com'
      DATA     COMPI/(0.0D0,1.0D0)/         

      ZPP=Z

C   * PRELIMINARIES *

      IF (IRETIN.EQ.IRETOUT)  RETURN              
      ZETA = Z * SQRT(Z) / 1.5D0

C   * SCALED AI & BI OUPUT REQUESTED *

      IF (IRETOUT.LT.0)  THEN

          IF (IRETIN.EQ.0)  THEN                  
              EXPZETA = EXP (ZETA)
              AI  = AI * EXPZETA
              AIP = AIP * EXPZETA
              IF (REAL(ZETA).GE.0.D0)  THEN
                  EXPMZETA = EXP (-ZETA)
                  BI  = BI * EXPMZETA
                  BIP = BIP * EXPMZETA
              ELSE
                  BI  = BI * EXPZETA
                  BIP = BIP * EXPZETA
              END IF

          ELSE IF (IRETIN.GT.0)  THEN             
              IF (ZPR(2).GE.0.D0)  THEN         
                  IF (REAL(ZETA).GE.0.D0)  THEN              
                     EXPTZETA = EXP (-2.D0*ZETA)
                     BI  = 2.D0*eim16*BI  + COMPI*AI*EXPTZETA
                     BIP = 2.D0*eim56*BIP + COMPI*AIP*EXPTZETA
                  ELSE                   
                     EXPTZETA = EXP (2.D0*ZETA)
                     BI  = 2.D0*eim16*BI*EXPTZETA  + COMPI*AI
                     BIP = 2.D0*eim56*BIP*EXPTZETA + COMPI*AIP
                  END IF
              ELSE
                  IF (REAL(ZETA).GE.0.D0)  THEN              
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = 2.D0*ei16*BI  - COMPI*AI*EXPTZETA
                    BIP = 2.D0*ei56*BIP - COMPI*AIP*EXPTZETA
                  ELSE                   
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = 2.D0*ei16*BI*EXPTZETA  - COMPI*AI
                    BIP = 2.D0*ei56*BIP*EXPTZETA - COMPI*AIP
                  END IF
              END IF
          END IF

C   * UNSCALED AI-BI OUTPUT REQUESTED *

        ELSE IF (IRETOUT.EQ.0)  THEN
          EXPMZETA = EXP (-ZETA)
          AI  = AI * EXPMZETA
          AIP = AIP * EXPMZETA

          IF (IRETIN.LT.0)  THEN                  
              IF (REAL(ZETA).GE.0.D0) THEN
                  EXPZETA = EXP (ZETA)
                  BI  = BI * EXPZETA
                  BIP = BIP * EXPZETA
              ELSE
                  BI  = BI * EXPMZETA
                  BIP = BIP * EXPMZETA
              END IF

          ELSE IF (IRETIN.GT.0) THEN              
              EXPZETA = EXP (ZETA)
              IF (ZPR(2).GE.0.D0)  THEN
                  BI  = 2.D0*eim16  * BI*EXPZETA  + COMPI*AI
                  BIP = 2.D0*eim56  * BIP*EXPZETA + COMPI*AIP
              ELSE
                  BI  = 2.D0*ei16  * BI*EXPZETA  - COMPI*AI
                  BIP = 2.D0*ei56  * BIP*EXPZETA - COMPI*AIP
              END IF
          END IF
          ZETA = (0.0D0, 0.0D0)

C   * NUMERICALLY STABLE SCALED AI OUTPUT REQUESTED *
C     (OBVIOUSLY IF YOU ARE WORKING FROM UNSCALED OR SCALED AI/BI
C      SOME PRECISION LOSS MUST OCCUR, BE FOREWARNED)

      ELSE IF (IRETOUT.GT.0)  THEN

          IF (IRETIN.LT.0)  THEN                  
              IF (ZPR(2).GE.0.D0)  THEN         
                  IF (REAL(ZETA).GE.0.D0)  THEN             
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = ei16*(BI  - COMPI*AI*EXPTZETA)*0.5D0
                    BIP = ei56*(BIP - COMPI*AIP*EXPTZETA)*0.5D0
                  ELSE                  
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = ei16*(BI*EXPTZETA  - COMPI*AI)*0.5D0
                    BIP = ei56*(BIP*EXPTZETA - COMPI*AIP)*0.5D0
                  END IF
              ELSE
                  IF (REAL(ZETA).GE.0.D0)  THEN             
                    EXPTZETA = EXP (-2.D0*ZETA)
                    BI  = eim16*(BI  + COMPI*AI*EXPTZETA)*0.5D0
                    BIP = eim56*(BIP + COMPI*AIP*EXPTZETA)*0.5D0
                  ELSE                  
                    EXPTZETA = EXP (2.D0*ZETA)
                    BI  = eim16*(BI*EXPTZETA  + COMPI*AI)*0.5D0
                    BIP = eim56*(BIP*EXPTZETA + COMPI*AIP)*0.5D0
                  END IF
              END IF

          ELSE IF (IRETIN.EQ.0)  THEN             
              EXPZETA = EXP (ZETA)
              EXPMZETA = EXP (-ZETA)
              IF (ZPR(2).GE.0.D0)  THEN
                  BI  = ei16*(BI  - COMPI*AI)*EXPMZETA*0.5D0
                  BIP = ei56*(BIP - COMPI*AIP)*EXPMZETA*0.5D0
              ELSE
                  BI  = eim16*(BI  + COMPI*AI)*EXPMZETA*0.5D0
                  BIP = eim56*(BIP + COMPI*AIP)*EXPMZETA*0.5D0
              END IF
              AI  = AI*EXPZETA
              AIP = AIP*EXPZETA
          END IF
      END IF
      RETURN
      END
ccc
      SUBROUTINE airy_only (Z,AI,AIP,ZETA)
      IMPLICIT COMPLEX*16 (A,B,E,X,Z)
      PARAMETER (ONEOOF=1D0/1.5D0)

c: Modified by EKW to only return Ai(z) assuming IRET=1 below
c: for use in h_space, where only one solution is required.
C     COMPLEX AIRY FUNCTION SUBROUTINE
C     Z; THE COMPLEX ARGUMENT
C     IRET; RETURN CONTROL PARAMETER FOR SCALING
C           < 0 RETURN SCALED AIRY FUNCTIONS OF THE FIRST AND SECOND TYPES.
C           = 0 RETURN UNSCALED AIRY FUNCTIONS OF THE FIRST AND SECOND TYPES.
C           > 0 RETURN SCALED RESULTS OF THE NUMERICALLY STABLE INDEPENDENT 
C               TYPE.  THESE ARE AI(Z) AND AI(Z*EXP([+,-]2IPI/3)).
C     AI AND AIP; THE COMPLEX AIRY FUNCTION OF THE FIRST TYPE AND ITS
C                   DERIVATIVE DEPENDING ON IRET.
C               IRET < 0, SCALED AIRY FUNCTION OF THE FIRST TYPE.  AI AND AIP
C                   ARE TO BE DIVIDED BY EXP(ZETA) TO OBTAIN THE TRUE AIRY 
C                   FUNCTIONS.
C               IRET = 0, UNSCALED AIRY FUNCTION OF THE FIRST TYPE
C               IRET > 0, SCALED AIRY FUNCTION OF THE FIRST TYPE.  AI AND AIP
C                   ARE TO BE DIVIDED BY EXP(ZETA) TO OBTAIN THE TRUE AIRY 
C                   FUNCTIONS.
C     BI AND BIP; THE COMPLEX AIRY FUNCTION OF THE SECOND TYPE AND ITS
C               DERIVATIVE DEPENDING ON IRET AND THE ANGLE OF Z.
C               IRET < 0, SCALED AIRY FUNCTION OF THE SECOND TYPE.  BI AND BIP 
C                   ARE TO BE MULTIPLIED BY EXP(ZETA), FOR ARG(Z) < 60 AND 
C                   DIVIDED BY EXP(ZETA) FOR ARG(Z) > 60.
C               IRET = 0, UNSCALED AIRY FUNCTION OF THE SECOND TYPE
C               IRET > 0, SCALED AIRY FUNCTION OF OF THE FIRST TYPE FOR THE
C     ARGUMENT Z*EXP(-2IPI/3) WHEN ARG(Z) > 0 AND FOR
C     Z*EXP(2IPI/3) WHEN ARG(Z) < 0. THE RESULTS ARE TO
C     BE MULTIPLIED BY EXP(ZETA). 
C      ZETA; THE EXPONTIAL SCALING FACTOR
C
C     THIS ROUTINE USES THE ALGORITHMS DESCRIBED IN "AN ALGORITHM FOR 
C     THE EVALUATION OF THE COMPLEX AIRY FUNCTIONS,"  Z. SCHULTEN, 
C     D.G.M. ANDERSON, AND R.G. GORDEN, IN JOUR. COMP. PHYSICS 31, 
C     PP. 60-75, (1979).
C     A BETTER DESCRIPTION OF THE SERIES EXPANSION IS AVAILABLE IN
C     HANDBOOK OF MATHEMATICAL FUNCTIONS, APPLIED MATHEMATICS SERIES #55,
C     M. ABRAMOWITZ AND I. A. STEGAN, P. 446 (1964).

C      PARAMETER (TPIO3=2.094395102393195D0, I=(0.0D0,1.0D0))
      REAL*8 TPIO3
      REAL*8 ZPR(2)
      include 'scairy_com'
      EQUIVALENCE (ZPR(1),ZPP)

      DATA TPIO3 /2.094395102393195D0/
C   * STATEMENT FUNCTION DEFINITION OF THE PHSE *

C      THETA (Z) = DATAN2 (DIMAG (Z), REAL (Z))

C   * Z OR CONJUGATE Z IN POSITIVE IMAGINARY PLANE *

      ZPP=Z

      IF (ZPR(2) .LT. 0.0D0)  THEN
          ZP = CONJG (Z)
          ICONJG = 1
      ELSE
          ZP = Z
          ICONJG = 0
      END IF


C   *********************************************************
C   *               *
C   *   SERIES SOLUTIONS FOR AI & BI WITHIN 3.5 OF ORIGIN   *
C   *               *
C   *********************************************************

      IF (ABS (Z) .LE. 3.5D0)  THEN
cekw      CALL ABSERIES2 (Z, AI,AIP,BI,BIP)
          call airy_only_series (Z, AI,AIP,ZETA)
cekw      IF (IRET.NE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,0,IRET)


C   ******************************************************************
C   *    *
C   *   ANGLES LESS THAN 120 DEGS. - AI(Z) AND AI(Z*EXP(-2IPI/3))    *
C   *   ARE CALCULATED FOR Z IN THE UPPER HALF PLANE, CONJGATE OF    *
C   *   Z IF IN THE LOWER.  SERIES OR GAUSSIAN SOLUTION AS NEEDED.   *
C   *    *
C   ******************************************************************

c: EKW: This was a bug originally because ZPR is equivalenced to original
c: input value of Z, not ZP, which is always in upper half plane:
cxx   ELSE IF (ATAN2(ZPR(2),ZPR(1)) .LE. TPIO3)  THEN
cxx   ELSE IF (abs(ATAN2(ZPR(2),ZPR(1))) .LE. TPIO3)  THEN
      ELSE IF (dimag(zp) .gt. -1.73205080756888d0*dreal(zp)) then
          IF (ABS(ZP).GT.7.9111D0)  THEN          
cekw          CALL ABGAUSS2 (ZP, AI, AIP, BI, BIP, ZETA)
              call AGAUSS2 (ZP,AI,AIP,ZETA)
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  AIP = CONJG (AIP)
              END IF
          ELSE
              IF (ABS (ZP-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
cekw              CALL ASERIES2 (ZP, AI,AIP, ZETA)
                  call airy_only_series(ZP, AI,AIP,ZETA)
              ELSE
                  CALL AGAUSS2 (ZP, AI, AIP, ZETA)
              END IF
              IF (ICONJG.EQ.1)  THEN              
                  ZETA = CONJG (ZETA)
                  AI  = CONJG (AI)
                  AIP = CONJG (AIP)
              END IF
          END IF
cekw      IF (IRET.LE.0)  CALL FORMAIRY (AI,BI,AIP,BIP,Z,ZETA,1,IRET)

C   **********************************************************
C   *                *
C   *   CONNECTION REGION FOR ANGLE GREATER THAN 120 DEGS.   *
C   *   CALCULATE AI(Z*EXP(-2IPI/3)) AND AI(Z*EXP(2IPI/3))   *
C   *   SOLUTIONS, USE CONNECTION FORMULAS FOR RESULTS.      *
C   *                *
C   **********************************************************

      ELSE 
c: z1=z*exp(-i2pi/3); z2=conjg[z*exp(i2pi/3] (to make z2 in upper half plane;
c: take conjugate of ai,aip,zeta afterwards):
          Z1 = Z * eim23    
          Z2 = CONJG (Z) * eim23                
          
          IF (ABS (Z1-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z1, AI1, AIP1, ZETA1)
          ELSE
              CALL AGAUSS2 (Z1, AI1, AIP1, ZETA1)
          END IF
          IF (ABS (Z2-(-0.9D0,2.8D0)) .LT. 4.97D0)  THEN
              CALL ASERIES2 (Z2, AI2, AIP2, ZETA2)
          ELSE
              CALL AGAUSS2 (Z2, AI2, AIP2, ZETA2)
          END IF
          AI2  = CONJG (AI2)  
          AIP2 = CONJG (AIP2)
          ZETA2 = CONJG (ZETA2)

c: Use connection formula: 
c: Ai(z)=exp(ipi/3)*Ai(z*exp(-i2pi/3)) + exp(-ipi/3)*z*exp(i2pi/3)
cekw          ZETA  = Z * SQRT(Z) * ONEOOF
cekw          E1  = EXP (ZETA-ZETA1)
cekw          E2  = EXP (ZETA-ZETA2)
cekw          AI  = ei13 * AI1*E1  + eim13 * AI2*E2
cekw          AIP = eim13 * AIP1*E1 + ei13 * AIP2*E2
cekw: Take out exp(zeta1) or exp(zeta2), whichever has smaller real part:
           if(dreal(zeta1) .lt. dreal(zeta2)) then
              zeta=zeta1
              E2=exp(zeta-zeta2)
              AI  = ei13 * AI1  + eim13 * AI2*E2
              AIP = eim13 * AIP1 + ei13 * AIP2*E2
           else
              zeta=zeta2
              E1=exp(zeta-zeta1)
              AI  = ei13 * AI1*E1  + eim13 * AI2
              AIP = eim13 * AIP1*E1 + ei13 * AIP2
           endif
      END IF

      RETURN
      END
ccc
      SUBROUTINE airy_only_series (Z, AI,AIP,ZETA)
      IMPLICIT COMPLEX*16 (A-H,O-Z)    
      REAL*8 C1,C2,FAC1,FAC2,FAC3
      DATA C1 /0.355028053887817D0/, C2 /0.258819403792807D0/

      ZETA = (0.0D0, 0.0D0)
c
      FC = (1.0D0, 0.0D0)
      GC = Z
      Z2 = Z * Z
      AI = C1*FC - C2*GC
      AIP = -C2

      DO 10 K = 3,600,3
         FAC1=1.0D0/(K-1)
         FAC2=1.0D0/K
         FAC3=1.0D0/(K+1)
         FC = FC * Z2 * FAC1    
         GC = GC * Z2 * FAC2
         AIP = AIP + (C1*FC - C2*GC)
         FC = FC * Z * FAC2       
         GC = GC * Z * FAC3
         AI = AI + (C1*FC - C2*GC)
         IF (ABS(FC).LT.1.D-17*ABS(AI) .AND. 
     &      ABS(GC).LT.1.D-17*ABS(AI))  GO TO 12
10    CONTINUE
      WRITE (6,*) ' SERIES APPROXIMATION FAILED TO CONVERGE'

12    continue
c
      RETURN
      END
      subroutine sr_geom(rng_srx,nsrcx,nrecx)
c
c: Computes horizontal ranges between nsrc sources at (xsrc,ysrc,zsrc)
c: and nrec receivers at (xrec,yrec,zrec) and places them in the
c: array rng_sr(1:nsrc,1:nrec).  Then places sorted ranges in 
c: range(1:nrng), where nrng=nsrc*nrec.
c  
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 nsrcx,nrecx,jsrc,jrec,jrng,j,jj,kk,jr,
     .   iibad,kk0,nzsr0,jd
      real*4 rng_srx(nsrcx,nrecx),rngcur,zp1,zp2
      real*8 pierad,rlast,delz
      data pierad/0.01745329251994/
c
      if(iimf .eq. 0 .and. iisig .eq. 0) nzmf=0
      call uni_space(nzmf,zmf,1.e0)
      call uni_space(nzs,zsrc,1.e0)
      call uni_space(nrec,zrec,1.e0)
      if(iigbs .ne. 0) then
         call uni_space(nth_gbs,th_gbs,1.e0)
         call uni_space(nb_gbs,b_gbs,1.e0)
      endif

      if(iigeom .eq. 1) call uni_space(nsrc,rkm,1.e0)
c
c: Find interfaces across which to compute terms for mode orthogonality:
      n_int=0
      if(iidiag .eq. -2 .and. (iimf .ne. 0 .or. iisig .ne. 0)) then
         do j=1,nlay-1
            jj=j+1
            if(zdep(j) .ge. zmf(1) .and. zdep(j) .le. zmf(nzmf) .and.
     .         mm(j) .eq. 1 .and. (iisol(j) .eq. 1 .or.
     .         iisol(jj) .eq. 1)) then
               n_int=n_int + 1
               j_int(n_int)=j
            endif
         enddo
         if(iicw .eq. 1) then
            delz=.01*cfmin/fcw(1)
         else
            delz=.01*cfmin/fmin
         endif
      endif
c
      nzsr=nzs + nrec + nzmf + 2*n_int
c: Don't include duct depths in zsr if doing real axis version:
      if(nduct .gt. 1 .and. iirx .le. 0) nzsr=nzsr + nduct
      nm_max2=NSR_NM_MAX/max(1,nzsr)
      iibad=0
      call mem_lim(nzsr,NSRMAX,MLINE,LML,'nzsr',4,'NSRMAX',6,iibad,0)
      call mem_lim(nsrc*nrec,NSNRMAX,MLINE,LML,'nsrc*nrec',9,
     .   'NSNRMAX',7,iibad,0)
      if(iibad .eq. 1) then
         write(*,*)'nzsr,NSRMAX,nsrc,nrec,NSNRMAX'
         write(*,*)nzsr,NSRMAX,nsrc,nrec,NSNRMAX
         stop 'stopped from sr_geom'
      endif
c
      if(iigeom .eq. 1) then
         do j=1,nsrc
            xsrc(j)=rkm(j)
cpln            xsrc(j)=rkm(j)*1000.e0
            ysrc(j)=0.
cpln            t_src(j)=rkm(j)
            t_src(j)=rkm(j)/1000.e0
            range(j)=rkm(j)
         enddo
c: Compute receiver depths:
         do j=1,nrec
            xrec(j)=0.
            yrec(j)=0.
            rec_lab(j)=zrec(j)
         enddo
      endif
c
c: Check which ducts have sources and/or receivers in them:
      do jd=1,nduct
         jd_ch(jd)=0
         if(jduct(1,jd) .eq. nlay .or. jduct(1,jd) .eq. 1) goto 15
         zp1=zpeak(jd)
         zp2=zpeak(jd+1)
         do jrec=1,nrec
            if(zrec(jrec) .ge. zp1 .and. zrec(jrec) .le. zp2) then
               jd_ch(jd)=1
               goto 15
            endif
         enddo
         do j=1,nzs
            if(zsrc(j) .ge. zp1 .and. zsrc(j) .le. zp2) then
               jd_ch(jd)=1
               goto 15
            endif
         enddo
         do j=1,nzmf
            if(zmf(j) .ge. zp1 .and. zmf(j) .le. zp2) then
               jd_ch(jd)=1
               goto 15
            endif
         enddo
15       continue
      enddo
c
c: Create an array zsr of src depths zsrc, rec depths zrec, and mode function
c: depths zmf:
      do jrec=1,nrec
         zsr(jrec)=zrec(jrec)
         zsr_im_gbs(nrec+j)=0.d0
         zsrmin=min(zsrmin,zrec(jrec))
         zsrmax=max(zsrmax,zrec(jrec))
      enddo
      do j=1,nzs
         zsr(nrec+j)=zsrc(j)
         zsr_im_gbs(nrec+j)=0.d0
         zsrmin=min(zsrmin,zsrc(j))
         zsrmax=max(zsrmax,zsrc(j))
      enddo
      if(iigbs .ne. 0) then
         do j=1,nzs
            zsr_im_gbs(nrec+j)=-b_gbs(j)*sin(th_gbs(j)*pierad)
         enddo
      endif
      nzsr0=nrec+nzs
      do j=1,nzmf
         zsr(nzsr0+j)=zmf(j)
         zsr_im_gbs(nzsr0+j)=0.d0
      enddo
      nzsr0=nzsr0+nzmf
      if(nduct .gt. 1 .and. iirx .le. 0) then
         do j=1,nduct
            zsr(nzsr0+j)=zduct(j)
            zsr_im_gbs(nzsr0+j)=0.d0
         enddo
         nzsr0=nzsr0 + nduct
      endif
      do j=1,n_int
         zsr(nzsr0+2*j-1)=zdep(j_int(j))
         zsr_im_gbs(nzsr0+2*j-1)=0.d0
         zsr(nzsr0+2*j)=zdep(j_int(j)) + delz
         zsr_im_gbs(nzsr0+2*j)=0.d0
      enddo
c: Sort real depths in zsr, keeping imaginary depths in zsr_im_gbs:
      nzsr0=nzsr+1
      call hpsort_re_im(nzsr,zsr,zsr_im_gbs,zsr_indx(nzsr0))
c: Sort resulting indices so that zsr_indx(i)=original index in zsr of 
c: current i'th element:
      call hpsort_i4_indx(nzsr,zsr_indx(nzsr0),zsr_indx)
c
c: Discard duplicates:
      if(nzsr .gt. 1) then
         j=1
         nzsr0=nzsr
10       if(zsr(j) .eq. zsr(j+1) .and. 
     .      zsr_im_gbs(j) .eq. zsr_im_gbs(j+1)) then
            nzsr=nzsr-1
            do jj=j+1,nzsr
               zsr(jj)=zsr(jj+1)
               zsr_im_gbs(jj)=zsr_im_gbs(jj+1)
            enddo
c: Keep indices correct by decrementing those above j:
            do jj=1,nzsr0
               if(zsr_indx(jj) .gt. j) zsr_indx(jj)=zsr_indx(jj)-1
            enddo
            j=j-1
         endif
         j=j+1
         if(j .lt. nzsr) goto 10
      endif
c
c: Copy zsr into zsrgeom if only geometry
c: (reuse eigenvalues and mode functions)
c
      if((i_geom.eq.0).and.(i_call_porter.eq.1)) then
         do jj=1,nzsr
            zsrgeom(jj)=zsr(jj)
         end do
         nzsrgeom=nzsr
      end if
c
c: Find indices of zrec,zsrc,zmf in zsr:
      do jrec=1,nrec
         mzrec(jrec)=zsr_indx(jrec)
      enddo
      do j=nrec+1,nrec+nzs
         mzsrc(j-nrec)=zsr_indx(j)
      enddo
      nzsr0=nrec+nzs
      do j=nzsr0+1,nzsr0+nzmf
         mzmf(j-nzsr0)=zsr_indx(j)
      enddo
      nzsr0=nzsr0+nzmf
      if(nduct .gt. 1 .and. iirx .le. 0) then
         do j=nzsr0+1,nzsr0+nduct
            mzduct(j-nzsr0)=zsr_indx(j)
         enddo
         nzsr0=nzsr0 + nduct
      endif
      do j=nzsr0+1,nzsr0+2*n_int
         mzint(j-nzsr0)=zsr_indx(j)
      enddo
c
      if(iigeom .eq. 1) then
         nrng=nsrc
         call mem_lim(nrng,NRNGMAX,MLINE,LML,'nrng',4,'NRNGMAX',7,
     .      iibad,1)
         do jrng=1,nrng
            range(jrng)=xsrc(jrng)
         enddo
c
         call hpsort_indx(nrng,range,zsr_indx)
c
         do jrng=1,nrng
            rlast=range(jrng)
            jsrc=zsr_indx(jrng)
            nrec_jr(jrng)=nrec
            kk0=(jrng-1)*nrec
            krec_jr(jrng)=kk0
            do jrec=1,nrec
               jrec_jr(1,kk0+jrec)=jsrc
               jrec_jr(2,kk0+jrec)=jrec
               rng_srx(jsrc,jrec)=rlast
            enddo
         enddo
c
      else
c: Compute ranges between source positions and receivers:
         nrng=0
         do jsrc=1,nsrc
            do jrec=1,nrec
               rngcur=sqrt((xsrc(jsrc)-xrec(jrec))**2 + 
     .            (ysrc(jsrc)-yrec(jrec))**2)
               rng_srx(jsrc,jrec)=rngcur
               do jr=nrng,1,-1
c: Skip duplicate ranges, but keep track of # s/r pairs had that range:
                  if(rngcur .eq. range(jr)) then
                     nrec_jr(jr)=nrec_jr(jr) + 1
                     goto 45
                  endif
               enddo
               nrng=nrng + 1
               call mem_lim(nrng,NRNGMAX,MLINE,LML,'nrng',4,'NRNGMAX',
     .            7,iibad,1)
               range(nrng)=rngcur
               nrec_jr(nrng)=1
45          enddo
         enddo
c: Sort ranges:
         call hpsort(nrng,range)
c: Set krec_jr(1:nrng), the starting index in jrec_jr(1:2,1:nsrc*nrec):
         krec_jr(1)=0
         do jrng=2,nrng
            krec_jr(jrng)=krec_jr(jrng-1) + nrec_jr(jrng-1)
         enddo
c: Reset nrec_jr to 0 for use in setting jrec_jr next:
         do jrng=1,nrng
            nrec_jr(jrng)=0
         enddo
         do jsrc=1,nsrc
            do jrec=1,nrec
c: Find index of rng_srx(jsrc,jrec) in range(1:nrng):
               rlast=rng_srx(jsrc,jrec)
               jrng=1
               call hunt(range,nrng,rlast,jrng)
               if(range(jrng) .ne. rlast) jrng=jrng + 1
c: Increment # s/r pairs that have this range:
               nrec_jr(jrng)=nrec_jr(jrng) + 1
               kk=krec_jr(jrng) + nrec_jr(jrng)
c: Keep track of source and rec index that had this range:
               jrec_jr(1,kk)=jsrc
               jrec_jr(2,kk)=jrec
            enddo
         enddo
      endif
c
      iibad=0
      if(rmax .eq. 0.e0 .and. nrng .gt. 0) rmax=range(nrng)/1000.e0
      if(iirx .le. 0 .and. rmax .le. 0.e0) then
         print *,'rmax = 0 in sr_geom!  Set rmax in _opt file.'
         iibad=1
      endif
      if(rmin .eq. 0.e0) then
         if(nrng .gt. 0) then
            rmin=min(998.d0,range(1)/1000.d0)
cpln            write(6,*)'rmin: ',rmin
        else
            print *,'rmin=0, but no S/R geometry given to compute'//
     .         ' rmin automaticallly.'
            iibad=1
         endif
      endif
c ONLY FOR ITWK
cpln      rmin=0.5
      if(rmin .lt. 0.e0) then
         nm_lim=nint(-rmin)
         rmin=999.
      else
         nm_lim=NM_MAX
      endif
      if(nsrc .gt. 0 .and. range(1) .le. 0.d0) then
         print *,'All ranges must be > 0!',(range(j),j=1,nrng)
         iibad=1
      endif
      if(iibad .eq. 1) stop
c: Compute maximum error in eigenvalues so that phase at maximum range
c: is within 1 degree of exact:
c: Be sure rmax is at least 50 water depths here:
cpln      write(6,*)
cpln          write(6,*)'Before RMAX: ',rmax
cpln      write(6,*)
      rmax=amax1(rmax,.05*sngl(zdep(nlay-1)))
c
cpln      write(6,*)
cpln        write(6,*)'rmin,rmax,zdep: ',rmin,rmax,0.05*sngl(zdep(nlay-1))
cpln      write(6,*)
cpln      do j=1,nsrc
cpln      write(6,*)'rkm =',(rkm(j),j=1,nsrc)
cpln      write(6,*)'range: ',(range(j),j=1,nsrc)
cpln      end do
cpln      write(6,*)
cpln      write(6,*)'Vertical tilt: ',dtiltvp(1),dtiltvp(nrec)
cpln      write(6,*)
cpln      write(6,*)'Horizontal tilt: ',dtilthp
cpln      write(6,*)
cpln      write(6,*)'fmin,fmax,df: ',fmin,fmax,df
cpln      write(6,*)faxbb(1),faxbb(nfbb)
cpln  write(6,*)
c      do jsrc=1,nzs
cpln       write(6,*)'sd= ',(zsrc(jsrc),jsrc=1,nzs)
c      end do
cpln  write(6,*)
cpln      do jrec=1,nrec
cpln      write(6,*)'rd= ',(zrec(jrec),jrec=1,nrec)
c      end do
c      write(6,*)
c      pause
c
c: Compute maximum error in eigenvalues so that phase at maximum range
c: is within 1 degree of exact:
c: Be sure rmax is at least 50 water depths here:
      rmax=amax1(rmax,.05*sngl(zdep(nlay-1)))
      errdkms=(twpie/(360.*1000.*rmax))**2
      errdk2=4.d0*errdkms
      if(phfac .lt. 2.e0) phfac=2.e0
      if(db_cut .eq. 0.e0) then
         db_cut=50.e0
      elseif(db_cut .lt. 30.e0) then
         print *,'Warning: db_cut low (>=30 recommended)...',db_cut
      endif
c
c: kim_fac=-ln[10**(-dB_down/20)] is used to find max range mode is
c: significant for (used in mode_field) see p. 150,131:
ccc   kim_fac=-dlog(10.**(-db_cut/20.d0))
ccc   dkim=kim_fac/rmin
c: Set max IM[kn] to allow for min range of interest rmin(km) (see pp.131,150):
      if(rmin .lt. 999.) then
         dkim=db_cut/(8685.9*rmin)
      else
         dkim=1.d100
      endif
c: kim_fac will be used with range in m in mode_field:
      kim_fac=db_cut/8.6859
c
      if(isp(nlay) .eq. 0 .and. allf(1) .eq. 1 .and. rmin .ge. 999.
     .   .and. cphmax .gt. geo(1,1,nlay)) then
         write(6,120) ' '
         write(6,120) 'WARNING: # branch line modes not '//
     .      'limited by cphmax: ',cphmax,sngl(geo(1,1,nlay))
         write(6,120) '   Use of RMIN instead highly recommended.'//
     .      '  Control-c to terminate now.'
         write(6,120) ' '
         write(lusvp,120) 'WARNING: # branch line modes not '//
     .      'limited by cphmax: ',cphmax,sngl(geo(1,1,nlay))
         write(lusvp,120) '   Use of RMIN instead highly recommended.'
120      format(a,f9.2,2x,f9.2)
      endif
c
      return
      end
      subroutine svp_check(iiwrt)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c: Local variables:
      integer*4 iiwrt,j,ii,ii2,jflubot,jflutop,nadd,ji,j1,j2,jj,i,
     .   jtop1,jbot2,inc
c
c: Set first top layer, last bottom layer depending on whether halfspaces
c: are homogeneous or given by first top layer or last bottom layer:
      jtop1=1
      if(geot(2,1,1) .lt. 0.d0) jtop1=2
      jbot2=nlayb+1
      if(geob(1,1,jbot2) .lt. 0.d0) jbot2=nlayb
c: Find the number of fluid layers in bottom and top layering:
      call last_flu(nlayt+1,jtop1,-1,geot,jflutop,allf(2))
      call last_flu(1,jbot2,1,geob,jflubot,allf(1))
c
c: NOTE: TOP LAYERS ASSUMED TO BE GIVEN FROM UPPER HALFSPACE TO OCEAN,
c: WHILE BOTTOM LAYERS ARE GIVEN FROM OCEAN BOTTOM TO LOWER HALFSPACE.
c
      jlfake(1)=0
      jlfake(2)=0
c: Copy top layering to master arrays:
      nlay=0
c: Make type of profile linear for halfspaces:
      ktt(1)=1
      ktb(nlayb+1)=1
      htop=0.d0
      do j=jtop1,nlayt+1
         if(j .gt. jtop1) htop=htop + ht(j)
c: Skip a layer between fluid and solid layers in case seismic modes
c: are desijred (see duct_check):
         if(j .eq. jflutop) then
            nlay=nlay+1
            h(nlay)=0.d0
            jlfake(2)=nlay
         endif
c: Copy geoacoustic parameters of top layers to master array:
         nlay=nlay+1
         h(nlay)=ht(j)
         do ji=1,2
            do j2=1,5
               geo(ji,j2,nlay)=geot(ji,j2,j)
            enddo
         enddo
         if(ktt(j) .ne. 1) then
            call blug(geo,h,ktt(j),bpt(1,j),nlay,nadd)
            if(iiwrt .eq. 1 .and. nadd .gt. 0) then
               write(lusvp,209)
209            format('   BLUG layer subdivided: ')
               do i=nlay-nadd,nlay
                  write(lusvp,208) h(i),((geo(ji,ii,i),ji=1,2),ii=1,5)
208               format(11(f10.4,1x))
               enddo
            endif
         endif
      enddo
      if(jflutop .eq. nlayt+2) then
         nlay=nlay+1
         h(nlay)=0.
         jlfake(2)=nlay
      endif
c: Copy SVP layers to master arrays:
      jsurf=nlay+1
      do j=1,nsvp-1
         nlay=nlay+1
         h(nlay)=zsvp(j+1) - zsvp(j)
         geo(1,1,nlay)=csvp(j)
         geo(2,1,nlay)=csvp(j+1)
         geo(1,3,nlay)=rho_svp
         geo(2,3,nlay)=rho_svp
         geo(1,4,nlay)=alpha_svp
         geo(2,4,nlay)=alpha_svp
      enddo
      jobot=nlay
c
      if(jflubot .eq. 0) then
         nlay=nlay+1
         h(nlay)=0.
         jlfake(1)=nlay
      endif
c: Copy layers in bottom layering to master arrays:
      do j=1,jbot2
         nlay=nlay+1
c: Copy geoacoustic parameters of bottom layers to master array:
         h(nlay)=hb(j)
         do ji=1,2
            do j2=1,5
               geo(ji,j2,nlay)=geob(ji,j2,j)
            enddo
         enddo
         if(ktb(j) .ne. 1) then
            call blug(geo,h,ktb(j),bpb(1,j),nlay,nadd)
            if(iiwrt .eq. 1 .and. nadd .gt. 0) then
               write(lusvp,209)
               do i=nlay-nadd,nlay
                  write(lusvp,208) h(i),((geo(ji,ii,i),ji=1,2),ii=1,5)
               enddo
            endif
         endif
c: If we just put in the last fluid layer, skip a layer in case seismic 
c: modes are desired (see duct_check):
         if(j .eq. jflubot) then
            nlay=nlay+1
            h(nlay)=0.
            jlfake(1)=nlay
         endif
      enddo
c
c: Make upper and lower halfspaces homogeneous:
cxx   do j=1,5
cxx      geo(1,j,1)=geo(2,j,1)
cxx      geo(2,j,nlay)=geo(1,j,nlay)
cxx   enddo
c
c: Make fake layers transparent for now:
c: Make fake layers have density and sound speeds close to what we
c: will set them to later so that mismatch variable mm set correctly:
      do ii=1,2
         j1=jlfake(ii)
         if(j1 .ne. 0) then
            ii2=3-ii
            inc=2*ii - 3
cc          j2=j1 + inc
            j2=j1 - inc
            geo(ii,1,j1)=geo(ii2,2,j2)
            geo(ii,2,j1)=0.d0
            geo(ii,3,j1)=geo(ii2,3,j2)
            geo(ii,4,j1)=geo(ii2,5,j2)
            geo(ii,5,j1)=0.d0
            do jj=1,5
               geo(ii2,jj,j1)=geo(ii,jj,j1)
            enddo
         endif
      enddo
c
      return
      end
ccc
      subroutine svp_check2
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c: Local variables:
      integer*4 j,ii,j1,iibad,jj,mm_dn,mm_up,mm_upx,jmin,jjmin,
     .   jflu1,jflu2,nd_next,jlast
      real*8 mmx(5),c3,c4,c_up,c_dn,c_upx,cpint,rhoratx,rhomin(2),
     .   alpmin(2)
      data mmx/.05,.05,.001,.0005,.0005/
c
      nthorpe=0
      do j=1,nlay
c: Frequency-dependent thickness flag:
         ihf(j)=0
c: Convert p-wave attenuations:
         if(geo(1,4,j) .eq. 999.d0) then
            nthorpe=nthorpe + 1
            jthorpe(nthorpe)=j
         else
            call alpha_conv(geo(1,4,j),geo(1,4,j),geo(1,1,j))
            call alpha_conv(geo(2,4,j),geo(2,4,j),geo(2,1,j))
         endif
c: Convert s-wave attenuations:
         call alpha_conv(geo(1,5,j),geo(1,5,j),geo(1,2,j))
         call alpha_conv(geo(2,5,j),geo(2,5,j),geo(2,2,j))
      enddo
c
      if(iirx .eq. 1) then
c: Don't allow gradients if real axis option chosen:
         iiAih(1)=-1.
         iiAih(2)=-1.
      endif
c: Put gradient(s) into lower halfspace if desired:
      call airy_hsp(1,2,iiAih(1),geo(1,1,nlay),Aih_mag(1),Aih_dir(1),
     .   ihf(nlay),ii_xi_mv(1),h(nlay),f_max,pie,cphlo)
c: Put gradient(s) into upper halfspace if desired:
      call airy_hsp(2,1,iiAih(2),geo(1,1,1),Aih_mag(2),Aih_dir(2),
     .   ihf(1),ii_xi_mv(2),h(1),f_max,pie,cphlo)
c
c: Compute absolute depths of layer interfaces (ocean surface is at z=0):
c: Note that zdep(j) is depth of BOTTOM of layer j.
      zdep(1)=-htop
      do j=2,nlay
         zdep(j)=zdep(j-1) + h(j)
      enddo
      Htot=zdep(nlay-1) - zdep(1)
c
c
      do j=1,nlay
c: Isospeed in cp flag:
         isp(j)=1
         if(geo(1,1,j) .ne. geo(2,1,j) .or.
     .      geo(1,4,j) .ne. geo(2,4,j)) isp(j)=0
c: Solid layer flag:
         iisol(j)=0
         iss(j)=0
         if(geo(1,2,j) .gt. 0. .or. geo(2,2,j) .gt. 0.) then
            iisol(j)=1
c: Check for restrictions on cs/cp and as/ap:
            call crat_check(geo(1,1,j),geo(1,2,j),geo(1,4,j),
     .         geo(1,5,j),j,iicw,iirx)
            call crat_check(geo(2,1,j),geo(2,2,j),geo(2,4,j),
     .         geo(2,5,j),j,iicw,iirx)
c: Isospeed in cs flag:
            if(geo(1,2,j) .eq. geo(2,2,j) .and.
     .         geo(1,5,j) .eq. geo(2,5,j)) iss(j)=1
         endif
      enddo
c
c: Check for mismatch across all layer interfaces:
      do j=1,nlay-1
c: Set shear attenuations to zero if shear speeds are zero:
         if(geo(1,2,j) .eq. 0.d0) geo(1,5,j)=0.d0
         if(geo(2,2,j) .eq. 0.d0) geo(2,5,j)=0.d0
         j1=j+1
c: Check for mismatch from layer j to layer j+1:
         mm(j)=1
         do jj=1,5
           if(abs(geo(1,jj,j1)-geo(2,jj,j)) .gt. mmx(jj)) goto 99
         enddo
c: Parameters across interface close enough. Make exactly equal:
         mm(j)=0
         do jj=1,5
           geo(2,jj,j)=geo(1,jj,j1)
         enddo
99       continue
      enddo
c
c: Last solid layer in top layering:
      jsol(2,2)=2
c: Last solid layer in bottom layering:
      jsol(1,2)=nlay-1
c
c: jflu(ii,1) gives the first fluid layer, jflu(ii,2) gives the last fluid 
c: layer, jflu(ii,3) gives increment, jflu(ii,4) gives constant needed
c: in rp_calc for direction ii (1=down,2=up):
      if(allf(1) .eq. 1) then
         jsol(1,1)=nlay
         jflu(1,2)=nlay-1
      else
         do j=nlay,1,-1
            if(geo(1,2,j) .eq. 0.d0 .and. geo(2,2,j) .eq. 0.d0) then
c: If first fluid layer going up, set jsol1 to layer below, jflu2 to 
c: layer above (taking into account room for fake layer):
               jsol(1,1)=j+1
               jflu(1,2)=j-1
               goto 60
            endif
         enddo
60       continue
      endif
      if(allf(2) .eq. 1) then
         jsol(2,1)=1
         jflu(2,2)=2
      else
         do j=1,nlay
            if(geo(1,2,j) .eq. 0.d0 .and. geo(2,2,j) .eq. 0.d0) then
c: If first fluid layer going down, set jsol1 to layer above, jflu2 to 
c: layer below (taking into account room for fake layer):
               jsol(2,1)=j-1
               jflu(2,2)=j+1
               goto 62
            endif
         enddo
62       continue
      endif
c
c: Find local minima (or ducts) in fluid layers:
      nduct=0
      nd_next=1
      cfmin=1.d+100
c
      jflu1=jflu(2,2)
      jflu2=jflu(1,2)
      jlast=1
c
c: Always put duct at surface in SVP is increasing:
      mm_up=1
cc    rhoratx=geo(2,3,jflu1-1)/geo(1,3,jflu1)
cc    if(rhoratx .gt. 1.5d0 .or. rhoratx .lt. .667d0) mm_up=1
      mm_upx=mm_up
      c_upx=geo(2,1,jflu1-1)
c
      zpeak(1)=zdep(1)
      do j=jflu1,jflu2
c: c3,c4 are sound speeds at top and bottom of current layer:
         c3=geo(1,1,j)
         c4=geo(2,1,j)
c
         c_up=geo(2,1,j-1)
         if(c3 .eq. geo(2,1,j-1)) c_up=geo(1,1,j-1)
         c_dn=geo(1,1,j+1)
         if(c4 .eq. geo(1,1,j+1)) c_dn=geo(2,1,j+1)
c
         rhoratx=geo(2,3,j)/geo(1,3,j+1)
         mm_dn=0
cxx      if(rhoratx .gt. 1.5d0 .or. rhoratx .lt. .667d0 .or.
         if(rhoratx .gt. 1.1d0 .or. rhoratx .lt. .909d0 .or.
     .      j .eq. jflu2) mm_dn=1
c
         if(c3 .eq. c4) then
c: Isospeed: c_up,c_dn are speeds to compare to:
            if((c3 .lt. c_up .or. mm_up .eq. 1 .or.
c: Case of two isospeeds together:
     .         (c3 .eq. c_up .and. (c3 .lt. c_upx .or. mm_upx .eq. 1)))
     .         .and. (c4 .lt. c_dn .or. mm_dn .eq. 1)) then
               call duct_enter(j,1,c3,0,nd_next)
            endif
            if((c3 .gt. c_up .and. c4 .gt. c_dn) .or. 
     .         mm_dn .eq. 1) then
               call peak_enter(nd_next,j,jlast,c4,c_dn,mm_dn)
            endif
            if(c3 .ne. c_up) then
               c_upx=c_up
               mm_upx=mm_up
            endif
         elseif(c3 .lt. c4) then
c: Positive gradient with depth, check for duct at top of layer:
            if(mm_up .eq. 1) then
               call duct_enter(j,1,c3,0,nd_next)
            endif
c: Check if peak at bottom of layer:
            if(c4 .gt. c_dn .or. mm_dn .eq. 1) then
               call peak_enter(nd_next,j,jlast,c4,c_dn,mm_dn)
            endif
            c_upx=c_up
            mm_upx=mm_up
         elseif(c4 .lt. c3) then
c: Check if peak at top of layer:
            if(c3 .gt. c_up) then
               call peak_enter(nd_next,j-1,jlast,c4,c_dn,mm_up)
            endif
c: Negative gradient with depth, check for duct at bottom of layer:
            if(c4 .lt. c_dn .or. mm_dn .eq. 1) then
               call duct_enter(j,2,c4,0,nd_next)
            endif
c: Check if peak at bottom of layer:
            if(mm_dn .eq. 1) then
               call peak_enter(nd_next,j,jlast,c4,c_dn,mm_dn)
            endif
            c_upx=c_up
            mm_upx=mm_up
         endif
         mm_up=mm_dn
      enddo
c: Add on bottom layers thickness to last duct width:
      nd_next=nd_next + 1
      zpeak(nduct+1)=zdep(nlay-1)
      dz_duct(nduct)=dz_duct(nduct) + zdep(nlay-1)-zdep(jlast)
c
c: Place duct at top of lower halfspace if gradient and all fluid:
      if(isp(nlay) .eq. 0 .and. allf(1) .eq. 1) then 
         call duct_enter(nlay,1,geo(1,1,nlay),0,nd_next)
         zpeak(nduct+1)=zdep(nlay-1)
      endif
c: Place duct at bottom of upper halfspace if gradient and all fluid:
      if(isp(1) .eq. 0 .and. allf(1) .eq. 1) then 
         call duct_enter(1,1,geo(2,1,1),0,nd_next)
         zpeak(nduct+1)=zdep(nlay-1)
      endif
c
c: Temp to place ducts at user's request:
      if(iidiag .eq. -8) then
55       continue
         print *,'Enter layer#, ii(1=top,2=bot) to place ref depth: '
         read(5,*) j,ii
         if(j .ne. 0) then
            c4=geo(ii,1,j)
            call peak_enter(nduct,zdep,j,jlast,c4,c4,1)
            call duct_enter(j,ii,geo(ii,1,j),0,nd_next)
            goto 55
         endif
      endif
c
      if(nduct .eq. 0 .or. cfmin .ge. 1.d100) then
         print *,'nduct=0 or cfmin bad: BUG - tell EKW ',
     .      (csvp(j),j=1,nsvp)
         stop
      endif
c
c: Sort dz_duct:
c: Sort ducts by (decreasing) attenuation rather than by increasing thickness:
      do j=1,nduct
         jj=jduct(1,j)
         ii=jduct(2,j)
         if(jj .ne. nlay) then
            dz_duct(j)=1.d0/max(1.d-20,geo(ii,4,jj))
         else
c: Make sure branch line duct is last duct checked:
            dz_duct(j)=-1.d0
         endif
      enddo
      call hpsort_indx(nduct,dz_duct,indx_duct)
c: TEMP:
cc    print *,'Change order of ducts?'
cc    read(5,*) isvmin
cc    if(isvmin .eq. 1) then
cc       isvmin=indx_duct(nduct)
cc       indx_duct(nduct)=indx_duct(1)
cc       indx_duct(1)=isvmin
cc    endif
c
      kduct0=indx_duct(nduct)
c
      nsvmin0=jduct(1,kduct0)
      nsvmin=nsvmin0
      isvmin=jduct(2,kduct0)
      kduct=kduct0
      rho_duct=geo(isvmin,3,nsvmin)
      jflu(1,1)=nsvmin
      jflu(1,3)=1
      jflu(1,4)=0
      jflu(1,5)=isvmin-1
      jflu(2,1)=nsvmin
      jflu(2,3)=-1
      jflu(2,4)=-1
      jflu(2,5)=isvmin-2
      jhsp(1)=nlay
      jhsp(2)=1
c
c: Set cpfake for seismic mode checking in duct_check:
      cphlo=cfmin
      if(cphmin .eq. 0.e0) then
         cpfake(1)=0.d0
         cpfake(2)=0.d0
      else
         do ii=1,2
            if(allf(ii) .eq. 1) then
               cpfake(ii)=0.d0
            else
               jjmin=3-ii
               jmin=jflu(ii,2)
               cpint=geo(jjmin,1,jmin)
c: Start csmin at fluid sound speed at interface since c_sch=.92*min(cs,cw):
               csmin=cpint
               rhomin(ii)=geo(jjmin,3,jmin)
               alpmin(ii)=geo(jjmin,4,jmin)
               do j=jsol(ii,1),jsol(ii,2)+jflu(ii,3),jflu(ii,3)
                  do jj=1,2
                     if(geo(jj,2,j) .gt. 10.) then
                        if(geo(jj,2,j) .lt. csmin) then
                           csmin=geo(jj,2,j)
                           rhomin(ii)=geo(jj,3,j)
                           alpmin(ii)=geo(jj,4,j)
                           jjmin=jj
                           jmin=j
                        endif
                     endif
                  enddo
               enddo
               if(cphmin .lt. 0.e0) then
                  cpfake(ii)=0.85*csmin
               elseif(cphmin .ge. cpint) then
                  cpfake(ii)=0.d0
               else
                  cpfake(ii)=dmax1(dble(cphmin),0.85d0*csmin)
               endif
               cphlo=min(cphlo,cpfake(ii))
            endif
         enddo
      endif
c
c: Include ducts for seismic mode checking:
      do ii=1,2
         if(cpfake(ii) .ne. 0.d0) then
ccf         j=jflu(ii,2) + jflu(ii,3)
            j=jlfake(ii)
            geo(1,1,j)=cpfake(ii)
            geo(2,1,j)=cpfake(ii)
ccf         geo(1,3,j)=geo(jjmin,3,jmin)
ccf         geo(2,3,j)=geo(1,3,j)
            geo(1,3,j)=rhomin(ii)
            geo(2,3,j)=rhomin(ii)
            geo(1,4,j)=alpmin(ii)
            geo(2,4,j)=alpmin(ii)
            call duct_enter(j,3-ii,1.d100,1,nd_next)
c: Make this duct be done last:
            do j=nduct,2,-1
               indx_duct(j)=indx_duct(j-1)
            enddo
            indx_duct(1)=nduct
         endif
      enddo
c
c: Check for isospeed layers next to no-mismatch interfaces (causes
c: problems in rp_nomm when on negative sheets):
      do j=1,nlay-1
         if(mm(j) .eq. 0 .and. h(j+1) .ne. 0.d0) then
cc       if(geo(2,1,j) .eq. geo(1,1,j+1) .and. 
cc   .      geo(2,4,j) .eq. geo(1,4,j+1) .and. h(j+1) .ne. 0.d0) then
            if(isp(j) .eq. 1 .and. isp(j+1) .eq. 1) then
               print *,'ORCA CANNOT HANDLE ISOSPEED P-WAVE LAYERS '//
     .   'ON BOTH SIDES OF NO-MISMATCH-IN-P INTERFACE!! SORRY!'
               print *,'Layer # ',j,j+1
               print *,'Put in slight discontinuity or gradient'
               iibad=1
               jjfail=1
               return
            endif
         endif
         if(iisol(j) .eq. 1 .and. iisol(j+1) .eq. 1) then
            if(geo(2,2,j) .eq. geo(1,2,j+1) .and. 
     .         geo(2,5,j) .eq. geo(1,5,j+1) .and. h(j+1) .ne. 0.d0 .and.
     .         iss(j) .eq. 1 .and. iss(j+1) .eq. 1) then
               print *,'ORCA CANNOT HANDLE ISOSPEED S-WAVE LAYERS '//
     .   'ON BOTH SIDES OF NO-MISMATCH-IN-S INTERFACE!! SORRY!'
               print *,'Layer # ',j,j+1
               print *,'Put in slight discontinuity or gradient'
               iibad=1
               jjfail=1
               return
            endif
         endif
      enddo
cpln      if(iibad .eq. 1) stop
c
      do j=1,nlay-1
         j1=j + 1
c: Set density ratio:
         if(mm(j) .eq. 0) then
            rhorat(j)=1.d0
         elseif(j .ge. nsvmin) then
            rhorat(j)=geo(2,3,j)/geo(1,3,j1)
         else
            rhorat(j)=geo(1,3,j1)/geo(2,3,j)
         endif
      enddo
c: EKW FIX 12/13/95: Account for fact that a fake layer with different 
c: density may be there:
      if(allf(1) .eq. 0) then
cbb      rholay(1)=geo(2,3,jsol(1,1)-1)/geo(1,3,jsol(1,1))
         rholay(1)=geo(2,3,jflu(1,2))/geo(1,3,jsol(1,1))
      endif
      if(allf(2) .eq. 0) then
cbb      rholay(2)=geo(1,3,jsol(2,1)+1)/geo(2,3,jsol(2,1))
         rholay(2)=geo(1,3,jflu(2,2))/geo(2,3,jsol(2,1))
      endif
c
c: Convert cphmax from grazing angle to phase velocity if necessary:
      if(cphmax .lt. 0.) then
c: Interpret -cphmax as maximum grazing angle at min sound speed desired:
         cphmax=cfmin/cos(-cphmax*twpie/360.)
      endif
c
      crmax=cfmin
c: Find maximum halfspace sound speed:
      chspmax=.999d0*geo(1,1,nlay)
      if(geo(2,1,1) .gt. cfmin .and. geo(2,1,1) .lt. chspmax) then
c: If upper halfspace p-wave speed is larger than cfmin, but smaller
c: than lower halfspace p-wave speed, then is should be used to decide
c: when modes get leaky (usually does not happen when air above):
         chspmax=.999d0*geo(2,1,1)
      endif
c
      return
500   print *,'Error opening SVP file ',svp_file
      stop
      end
ccc
      subroutine last_flu(jlay1,jlay2,inc,geo,nlay_fl,allf)
c
      implicit none
      integer*4 jlay1,jlay2,inc,nlay_fl,j,ii_sol,allf
      real*8 geo(2,5,jlay2)
c
      allf=0
      ii_sol=0
      nlay_fl=jlay1-inc
      do j=jlay1,jlay2,inc
         if(geo(1,2,j) .gt. .0 .or. geo(2,2,j) .gt. .0) then
c: ii_sol set to 1 when first solid layer encountered: 
            ii_sol=1
c: Don't allow shear speed to start or end at zero:
            geo(1,2,j)=dmax1(.5d0,geo(1,2,j))
            geo(2,2,j)=dmax1(.5d0,geo(2,2,j))
         else
            if(ii_sol .eq. 0) then
c: All fluid layers up to here:
               nlay_fl=j
            else
c: If fluid layer detected below solid layer, make it slightly solid:
               geo(1,2,j)=.5
               geo(2,2,j)=.5
            endif
         endif
      enddo
c: If all fluid, set nlay_fl to one beyond:
      if(nlay_fl .eq. jlay2) then
         allf=1
         nlay_fl=nlay_fl + inc
      endif
c
      return
      end
ccc
      subroutine alpha_conv(alpha1,alpha2,c)
c
      implicit none
      real*8 alpha1,alpha2,c
c
      if(dabs(c) .lt. 1.d-100) return
      if(alpha1 .eq. 999.d0) then
c: For Thorpe attenuation set alpha later in freq_init:
         alpha2=0.d0
      elseif(alpha1 .ge. 0.d0) then
c: Attenuation given in dB/m-kHz:
         alpha2=alpha1
      else
c: Attenuation given in dB/lambda:
         alpha2=-alpha1/(.001d0*c)
      endif
c
      return
      end
ccc
      subroutine crat_check(cp,cs,ap,as,jlay,iicw,iirx)
c
      implicit none
      include 'Parms_com'
      real*8 cp,cs,as,ap,aratlim,cratlim,csmax,asmax,safe_fac
      integer*4 jlay,iicw,iirx
      data cratlim/0.86602540378444d0/
c
      csmax=cp*cratlim
      if(cs .gt. csmax) then
         print *,'WARNING: Positive compressibility requires ',
     .      'cs/cp < sqrt(3/4).'
         print *,'         Layer # ',jlay,'; cs,csmax = ', cs,csmax
         write(lusvp,100) jlay,cs,csmax
         cs=csmax
      endif
100   format('WARNING: Positive compressibility requires ',
     .      'cs/cp < sqrt(3/4).'/
     .      '          Layer # = ',i3,'; cs,csmax = ',f7.1,2x,f7.1)
c
c: See SAFARI manual, Sec. 2, p. 5 and ORCA II, p. 11 (note ap,as in
c: db/m-kHz, not db/lambda):
      aratlim=.75*(cp/cs)**3
      asmax=ap*aratlim
      if(iicw .eq. 2 .and. iirx .eq. -1) then
c: For mode-following in complex k plane, don't allow maximum shear
c: wave attenuation because modes can become very difficult to find:
         safe_fac=.75d0
      else
         safe_fac=1.0d0
      endif
      if(as .gt. safe_fac*asmax) then
         print *,'WARNING: Conservation of energy requires ',
     .      'as/ap < (3/4)*(cp/cs)**3.'
         print *,'         Layer # ',jlay,'; as,asmax = ',as,asmax
         write(lusvp,110) jlay,as,asmax
         as=safe_fac*asmax
      endif
110   format('WARNING: Conservation of energy requires ',
     .      'as/ap < (3/4)*(cp/cs)**3.'/
     .      '          Layer # = ',i3,'; as,asmax = ',f9.4,2x,f9.4)
c
      return
      end
ccc
      subroutine duct_enter(jlay,ii,cmin,iifake,nd_next)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jlay,ii,iifake,nd_next
      real*8 cmin
c
      if(nduct .ge. nd_next) then
         if(iidiag .ge. 1) then
            print *,'Duct skipped (two in a row): ',jlay
         endif
         return
      endif
c
c: Make sure not to put duct at no-mismatch interface to a halfspace:
      if(jlay .eq. nlay-1 .and. mm(nlay-1) .eq. 0 .and.
     .   isp(nlay-1) .eq. 1) then
         print *,'Duct skipped (no-mismatch lower h-space)'
         return
      elseif(jlay .eq. 2 .and. mm(1) .eq. 0 .and.
     .   isp(2) .eq. 1) then
         print *,'Duct skipped (no-mismatch upper h-space)'
         return
cpln
      elseif(jlay+ii-2 .lt. 1) then
         if(iiwrite .gt. 0)
     .        print *,'Zero or negative index in zdep'
         return
      end if
      nduct=nduct + 1
      jduct(1,nduct)=jlay
      jduct(2,nduct)=ii
      jduct(3,nduct)=iifake
      jduct(4,nduct)=jflu(1,2)
      jduct(5,nduct)=jflu(2,2)
      zduct(nduct)=zdep(jlay+ii-2)
      if(iifake .eq. 1 .or. iifake .eq. 2) then
         jduct(3+iifake,nduct)=jflu(iifake,2) + jflu(iifake,3)
      endif
c: FIX 2-27-95.  Min sound speed need not be in water column:
cxx   if(cmin .lt. cfmin .and. jlay .ge. jsurf .and. 
cxx  .   jlay .le. jobot) then
      if(cmin .lt. cfmin) then
         cfmin=cmin
      endif
c
      return
      end
ccc
      subroutine peak_enter(nd_next,jlay,jp,c4,c_dn,mm_dn)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 nd_next,jlay,jp,mm_dn
      real*8 c4,c_dn,dwid,lwid
c
      if(nd_next .gt. nduct) then
         if(iidiag .ge. 1) then
            print *,'Peak skipped (two in a row): ',jlay
         endif
         return
      endif
c
      dwid=zdep(jlay) - zdep(jp)
      lwid=dwid/(min(c4,c_dn)/f_max)
c: Do not include duct if too thin in terms of wavelengths (unless
c: a significant portion of entire waveguide thickness):
      if((mm_dn .eq. 1 .or. lwid .gt. .2d0) .and.
     .   (lwid .gt. .05d0 .or. dwid/Htot .gt. .5d0)) then
c: Temp to keep fake ducts for testing:
cc   .   lwid .lt. .0001)) then
         dz_duct(nd_next)=dwid
         cspan(nd_next)=max(c4,c_dn)
         jp=jlay
         nd_next=nd_next + 1
         zpeak(nd_next)=zdep(jlay)
      else
c: If previous duct less than lambda/5, delete it:
         nduct=nduct-1
         if(iidiag .ge. 1) then
            print *,'Duct with width < lambda/5 deleted.'
         endif
      endif
c
      return
      end
c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************
      subroutine svp_read_pro
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'debug_com'
      include 'sector_env_com'
      include 'depth_com'
      integer*4 ISECT,NSCTOR,NSECT
c
c: Local variables:
c      integer*2 jc
c      integer*4 jd,nshft,templuinp
      integer*4 templuinp
c      integer*4 nline,j,j1,j2,iiblug(-1:4),iierr,ndel,nr
      integer*4 nline,j,iierr,nfcr
c      real*8 zdel(NLMAX),hbb(NLMAX),geobb(2,5,NLMAX),bpbb(2,NLMAX)
c fbv
      integer  MINMOD, MAXMOD, MODCUT
      real*8 dummy1, dummy2, seclen
      character*64 eline

cc       COMMON /NN/  H0BEG, H1BEG, NFREQ_FMC, MSP_FMC

      data eline/'INVALID INPUT IN SVP FILE: '/
c      data iiblug/0,0,0,2,2,2/
c
c
      write(lusvp,95)
95    format(/'### SVP FILE INFORMATION ###')
c
c      open(luinp, err=500,status= 'OLD')
c      open(luinp, file='cprosim_in',err=500,status= 'OLD')
      templuinp = luinp
      iierr=0
      nline=0
      svp_ver= 2.0
      rho_svp= 1.0
      alpha_svp= 0.0
      iidebug=0
c
c

      nfcr = nfcw
      READ(luinp,*)NSECT
      READ(luinp,*,ERR=1525) MINMOD, MAXMOD, MODCUT
c      READ(luinp,*,ERR=1525) MINMOD, cphmax, MODCUT
      write(6,*)'IICW: ',iicw
      write(6,*)'NSECT: ',NSECT
      write(6,*)'MINMOD, cphmax, MODCUT: ',MINMOD, cphmax, MODCUT
c
C REVIEW 1525 CONTINUE
 1525 CONTINUE

 490  CONTINUE

c: Sound speed profile in ocean:
      DO ISECT=1,NSECT
         READ(luinp,*)   r_h0(isect), dummy1, dummy2, seclen
c
c
c The distance from the source is stored in a common variable:
         secleng(ISECT) = seclen*1000.0
c
         if(iiwrite .eq.0) then
            write(6,*)'Water depth: ',r_h0(isect)
            write(6,*)'Sector length: ',secleng(isect)
         end if
c
cc      H0BEG= h0_fmc
         nline = 1
         read(luinp,*,end=510,err=510)R_Z0(1,ISECT),R_C0(1,ISECT)
         call check_val_r8(r_z0(1,isect),0.d0,0.d0,eline,27,
     .        nline,'zsvp(1) MUST BE 0',17,iierr)
         call check_val_r8(r_c0(1,isect),1.d-10,1.d+10,eline,27,
     .        nline,'csvp(1)',7,iierr)
         if(iidebug .eq. 1) then
            write(6,*)'Sound Speed in water at sector: ',isect
            write(6,*)R_Z0(1,ISECT),R_C0(1,ISECT)
         end if
         do j=2,nlmax
            read(luinp,*,end=510,err=510)R_Z0(j,ISECT),R_C0(j,ISECT)
c
            if(iidebug .eq. 1)
     .           write(6,*)R_Z0(j,ISECT),R_C0(j,ISECT)
c     
c     nline= nline+1
            call check_val_r8(r_z0(j,isect),r_z0(j-1,isect),1.d10,
     .           eline,27,nline,'zsvp(1) MUST BE 0',17,iierr)
            call check_val_r8(r_c0(j,isect),1.d-10,1.d+10,eline,27,
     .           nline,'csvp(1)',7,iierr)
c     
            if( r_z0(j,isect) .eq. r_h0(isect) )   go to 1200
         end do
         print *, ' sub svp_read, too many points in svp '
         stop
 1200    continue
c     
         nsvp= j
c     pln      write(*,*)'Number of sound speed points', nsvp
         R_ND0(ISECT)=nsvp
         secz(isect)=r_h0(isect)
c     
c     : Bottom layering:
         read(luinp,*)  R_H1(isect), R_R1(isect),R_BETA(1,isect)
         if(iidebug .eq. 1) then
            write(6,*)' Sediment values'
            write(6,*) R_H1(isect), R_R1(isect),R_BETA(1,isect)
         end if
         if( R_H1(isect) .gt. 0.0 ) then
            do j= 1, NLMAX
cpln Read attenuation in bottom:
cpln               read(luinp,*,end=510,err=510) R_Z1(j,isect),
cpln     .              R_C1(j,isect),R_BETA(j,isect)
cpln     .              write(6,*) R_Z1(j,isect), R_C1(j,isect), 
cpln     .                         R_BETA(j,isect)
               read(luinp,*,end=510,err=510) R_Z1(j,isect),
     .              R_C1(j,isect)
               if(iidebug .eq. 1)
cpln: Write attenuation in bottom
cpln     .              write(6,*) R_Z1(j,isect), R_C1(j,isect),
cpln     .                         R_BETA(j,isect)
     .              write(6,*) R_Z1(j,isect), R_C1(j,isect)
               if( r_z1(j,isect) .eq. r_h1(isect) )   go to 1400
            enddo
            print *,' sub opt_read, e^ successo un ...... '
            stop
 1400       continue
            nlayb=j
         else
            nlayb=0
         end if
cpln ONLY ASCOT01
cpln         if(nlayb.eq.2)
cpln     .     R_C1(2,isect)=R_C1(1,isect)
         R_ND1(ISECT)=nlayb

c: Lower halfspace:
cpln: Attenuation in sediment
cpln         read(luinp,*)  R_R2(isect), R_BETA(nlayb+1,isect), R_C2(isect)
cpln         read(luinp,*)  R_BETA(nlayb+2,isect),  R_C2S(isect)
         read(luinp,*)  R_R2(isect), R_BETA(2,isect), R_C2(isect)
         read(luinp,*)  R_BETA(3,isect),  R_C2S(isect)
c
         if(iidebug .eq. 1) then
            write(6,*)'Sub-bottom values ' 
cpln: Attenuation in sediment
cpln            write(6,*)R_R2(isect), R_BETA(nlayb+1,isect), R_C2(isect)
cpln            write(6,*)  R_BETA(nlayb+2,isect),  R_C2S(isect)
            write(6,*)R_R2(isect), R_BETA(2,isect), R_C2(isect)
            write(6,*)  R_BETA(3,isect),  R_C2S(isect)
         end if
      end do
c
c Use this line for DEC/HP/INTEL (unix) workstations
 480  continue
      ISECT=ISECT-1
      write(6,*)'ISECT: ',ISECT
c
c Use this line for SUN (sunos) workstations
c 480  continue
c
      NSCTOR = ISECT
c
      if(iierr .ne. 0) stop
c
      do ISECT = NSCTOR, 2, -1
       secleng(ISECT) = secleng(ISECT) - secleng(ISECT-1)
      enddo
      secz(NSCTOR+1)=secz(NSCTOR)
c
      f_max=fmax
      if (iicw .eq. 1) then
         f_max=fcw(1)
         do j=1,nfcw
            f_max=dmax1(dble(fcw(j)),f_max)
         enddo
         call hpsort_r4(nfcr,fcw)
         if(iiwrite .eq. 0)
     &      write(6,*)'Frequency band: ',fcw(1),fcw(nfcr)
      end if
c
      return
 500  print *,'Error opening input file'
      stop
 510  print *,'Endo or error reading input file at line ',nline
      stop
      end
      subroutine lname64(chstr,leng)
c
c: Subroutine to determine the length of a character
c: string without blanks
c: chstr -- the string (length 64)
c: leng -- the length without blanks
c
      character*64 chstr,blank
      character*1 char1
      data blank/'                                                  
     .              '/
c
      do 40 i=1,64
         char1 = chstr(i:i)
       	 if(char1 .eq. ' ') then
            leng = i-1
            chstr(i:64)=blank(i:64)
            return
       	 endif
40    continue
      leng=64
c
      return
      end
ccc
      subroutine lname80(chstr,leng)
c
c: Subroutine to determine the length of a character
c: string without blanks.  Starts at END of string.
c: chstr -- the string (length 80)
c: leng -- the length without blanks
c
      character*80 chstr
c
      do 40 i=80,1,-1
       	 if(chstr(i:i) .ne. ' ') then
            leng = i
            return
       	 endif
40    continue
      leng=0
c
      return
      end
ccc
      subroutine usage(progname,lprog)
c
c: This subroutine updates a file containing a list of the times
c: certain programs have been run and the users that ran them.
c
      integer kdat(3),ktim(3)
      character*20 progname
      character*64 home
c
      open(83,file='/usr/people/westwood/.usage',status='unknown',
     .     form='formatted',
     .   access='append')
      call getenv('HOME',home)
      call idate(kdat(1),kdat(2),kdat(3))
      call itime(ktim)
      write(83,100) progname(1:20),home(1:20),
     .   kdat(2),kdat(1),kdat(3),(ktim(j),j=1,3)
100   format(a20,'; ',a20,'; ',2(i2.2,'/'),i4,'; ',i2.2,2(':',i2.2))
      close(83)
c
      return
      end
ccc
      subroutine hpsort_re_im(n,ra_re,ra_im,indx)
c
      implicit none
      integer*4 n,indx(n)
      integer*4 i,ir,j,l,iia
      real*8 ra_re(n),ra_im(n)
c: Sorts an array ra_re(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra_re is replaced on output by its sorted 
c: rearrangement. ra_im is rearranged in the same manner as ra_re is.
c
      real*8 rra,rrb
c
      do j=1,n
         indx(j)=j
      enddo
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra_re(l)
            rrb=ra_im(l)
            iia=indx(l)
         else
            rra=ra_re(ir)
            rrb=ra_im(ir)
            iia=indx(ir)
            ra_re(ir)=ra_re(1)
            ra_im(ir)=ra_im(1)
            indx(ir)=indx(1)
            ir=ir-1
            if(ir .eq. 1) then
               ra_re(1)=rra
               ra_im(1)=rrb
               indx(1)=iia
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra_re(j) .lt. ra_re(j+1)) j=j+1
            endif
            if(rra .lt. ra_re(j)) then
               ra_re(i)=ra_re(j)
               ra_im(i)=ra_im(j)
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra_re(i)=rra
         ra_im(i)=rrb
         indx(i)=iia
      goto 10
c
      end
ccc
      subroutine hpsort(n,ra)
c
      implicit none
      integer*4 n
      integer*4 i,ir,j,l
      real*8 ra(n)
c: Sorts an array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted 
c: rearrangement.
c
      real*8 rra
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
         else
            rra=ra(ir)
            ra(ir)=ra(1)
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
      goto 10
c
      end
ccc
      subroutine hpsort_r4(n,ra)
c
      implicit none
      integer*4 n
      integer*4 i,ir,j,l
      real*4 ra(n)
c: Sorts an array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted 
c: rearrangement.
c
      real*4 rra
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      do l=1,n
         write(6,*)l,ra(l)
      end do
c
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
         else
            rra=ra(ir)
            ra(ir)=ra(1)
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
      goto 10
c
      end
ccc
      subroutine hpsort_indx(n,ra,indx)
c
      implicit none
      integer*4 n,indx(n)
      real*8 ra(n)
      integer*4 i,ir,j,l,iia
c: Sorts an array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      real*8 rra
c
      do j=1,n
         indx(j)=j
      enddo
c
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
            iia=indx(l)
         else
            rra=ra(ir)
            iia=indx(ir)
            ra(ir)=ra(1)
            indx(ir)=indx(1)
c
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               indx(1)=iia
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
         indx(i)=iia
      goto 10
c
      end
ccc
      subroutine hpsort_indx_c16(n,ra,indx)
c
      implicit none
      integer*4 n,indx(n)
      complex*16 ra(n)
      integer*4 i,ir,j,l,iia
c: Sorts real part of complex*16 array ra(1:n) into DESCENDING order 
c: using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      complex*16 rra
c
      do j=1,n
         indx(j)=j
      enddo
c
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
            iia=indx(l)
         else
            rra=ra(ir)
            iia=indx(ir)
            ra(ir)=ra(1)
            indx(ir)=indx(1)
c
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               indx(1)=iia
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(dble(ra(j)) .gt. dble(ra(j+1))) j=j+1
            endif
            if(dble(rra) .gt. dble(ra(j))) then
               ra(i)=ra(j)
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
         indx(i)=iia
      goto 10
c
      end
ccc
      subroutine hpsort_i4(n,ra)
c
      implicit none
      integer*4 n
      integer*4 ra(n)
      integer*4 i,ir,j,l
c: Sorts an integer*4 array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      integer*4 rra
c
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
         else
            rra=ra(ir)
            ra(ir)=ra(1)
c
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
      goto 10
c
      end
ccc
      subroutine hpsort_i4_indx(n,ra,indx)
c
      implicit none
      integer*4 n,indx(n)
      integer*4 ra(n)
      integer*4 i,ir,j,l
c: Sorts an integer*4 array ra(1:n) into ascending order using the Heapsort
c: algorithm. n is input; ra is replaced on output by its sorted
c: rearrangement.
c
      integer*4 rra,iia
c
      do j=1,n
         indx(j)=j
      enddo
c
      if(n .lt. 2) return
c: The index l will be decremented from its initial value down to 1
c: during the "hiring" (heap creation) phase.  Once it reaches 1, the
c: index ir will be decremented from its initial value down to 1
c: during the "retirement and promotion" (heap selection) phase.
      l=n/2 + 1
      ir=n
10    continue
         if(l .gt. 1) then
            l=l-1
            rra=ra(l)
            iia=indx(l)
         else
            rra=ra(ir)
            iia=indx(ir)
            ra(ir)=ra(1)
            indx(ir)=indx(1)
c
            ir=ir-1
            if(ir .eq. 1) then
               ra(1)=rra
               indx(1)=iia
               return
            endif
         endif
         i=l
         j=l+l
20       if(j .le. ir) then
            if(j .lt. ir) then
               if(ra(j) .lt. ra(j+1)) j=j+1
            endif
            if(rra .lt. ra(j)) then
               ra(i)=ra(j)
               indx(i)=indx(j)
               i=j
               j=j+j
            else
               j=ir+1
            endif
            goto 20
         endif
         ra(i)=rra
         indx(i)=iia
      goto 10
c
      end
ccc
      subroutine openfftout(nfile,fftroot,lfft,fftfile,nfft)
c
c: This subroutine opens an FFT file with root name fftroot and 
c: suffix _fft for output.  Opens as direct access file.
c
      character*64 fftroot,fftfile
      logical qopen
c
c: open input _fft file:
      fftfile(1:lfft+4)=fftroot(1:lfft)//'_fft'
      lfft=lfft + 4
      inquire(nfile,opened=qopen)
      if(qopen) close(nfile)
c: reopen as direct access file with record length = 22+nfft words:
can   lenrec=4*(20+(nfft+2))
c: IRIS uses words rather than bytes:
      lenrec=(20+(nfft+2))
      open(nfile,file=fftfile(1:lfft),status='unknown',
     .   form='unformatted',
     .   access='direct',recl=lenrec)
c
      return
      end
ccc
      subroutine openfftout2(nfile,fftroot,lroot,fftfile,lfft,nfft,
     .   fs,fmin,fmax,xh,NRECL)
c
c: This subroutine opens an FFT file with root name fftroot and 
c: suffix _fft for output.  Opens as direct access file.
c: For fmin=-999., it uses the old FFT format, where all nff2/2 + 1
c: frequency bins are output. Otherwise, only bins from fmin to fmax 
c: are output.  xh(7)=1 means one frequency band from xh(15)*df to 
c: xh(16)*df is included in the file, where df=fs/nfft.  The total 
c: number of bins present is xh(16)-xh(15)+1.
c
      implicit none
      integer*4 nfile,lroot,lfft,nfft,NRECL,nf1,nf2,nffth1,lenrec,
     .   total_bins
      real*4 fs,fmin,fmax,xh(20),df
      character*64 fftroot,fftfile
c
      xh(5)=nfft
      xh(6)=fs
c: open input _fft file:
      fftfile(1:lroot+4)=fftroot(1:lroot)//'_fft'
      lfft=lroot + 4
c: Open as direct access file with record length = 22+nfft words:
c: SUN, Alliant has record length in bytes:
      df=fs/nfft
      nf1=nint(fmin/df) + 1
      nf2=nint(fmax/df) + 1
      nffth1=nfft/2 + 1
c: Set fmin=-999. to output all bins in the old way:
      if(fmin .eq. -999. .or. (nf1 .eq. 1 .and. nf2 .eq. nffth1)) then
         total_bins=nffth1
         xh(7)=0.
      else
         if(nf1 .lt. 0 .or. nf2 .gt. nffth1 .or. nf1 .gt. nf2) then
            print *,'Bad fmin,fmax,fs in openfftout2: ',fmin,fmax,fs
            stop
         endif
         total_bins=nf2-nf1+1
         xh(7)=1.
         xh(15)=nf1-1
         xh(16)=nf2-1
      endif
c: For direct acces files, NRECL should be 4 for lenrec in bytes
c: (SUN), 1 for lenrec in words (IRIS).  Set in Parms_com.
      lenrec=NRECL*(20 + 2*total_bins)
      open(nfile,file=fftfile(1:lfft),status='unknown',
     .   form='unformatted',
     .   access='direct',recl=lenrec)
c
      return
      end
ccc
      subroutine mem_lim(n,nlim,eline,le,vname,lv,lname,ll,iibad,
     .   iistop)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
c
      integer n,nlim,le,lv,ll,iibad,iistop
      character*128 eline,vname,lname
c
      if(n .gt. nlim) then
         print *,' '
         print *,eline(1:le)
         print *,'VARIABLE NAME = ',vname(1:lv),'; LIMIT = ',nlim,
     .      '; LIMIT NAME = ',lname(1:ll)
         print *,'ENTERED OR COMPUTED VALUE FOR THIS RUN = ',n
cc       do j=1,nlay+1
cc          print *,j,1,h(j),(geo(1,jj,j),jj=1,5)
cc          print *,j,2,h(j),(geo(2,jj,j),jj=1,5)
cc       enddo
         if(iistop .eq. 1) stop
         iibad=1
      endif
c
      return
      end
ccc
      subroutine uni_space(n,x,fac)
c
      implicit none
      integer*4 n,j
      real*4 x(10000),fac,xfac
c
      if(abs(n) .gt. 10000) then
         write(6,*)'Stop in uni_space; array index exceeded'
      end if
c
      if(n .lt. 0) then
         n=iabs(n)
         xfac=(x(2) - x(1))/max(1,n-1)
         do j=2,n
            x(j)=x(1) + float(j-1)*xfac
         enddo
         if(fac .ne. 1.e0) then
            do j=1,n
               x(j)=fac*x(j)
            enddo
         endif
      endif
c
      return
      end
ccc
      subroutine suffix(root,lroot,j,n,suf,lsuf,root2,lroot2)
c
      implicit none
      integer*4 lroot,j,n,lsuf,lroot2,ln
      character*128 root,suf,root2,nsuf
c
      if(n .eq. 1) then
         lroot2=lroot
         root2(1:lroot)=root(1:lroot)
      else
         if(n .lt. 10) then
            write(nsuf(1:1),'(i1)') j
            ln=1
         elseif(n .lt. 100) then
            write(nsuf(1:2),'(i2.2)') j
            ln=2
         elseif(n .lt. 1000) then
            write(nsuf(1:3),'(i3.3)') j
            ln=3
         elseif(n .lt. 10000) then
            write(nsuf(1:4),'(i4.4)') j
            ln=4
         else
            print *,'error in suffix: ',j,n
         endif
         lroot2=lroot + lsuf + ln
         root2(1:lroot2)=root(1:lroot)//suf(1:lsuf)//nsuf(1:ln)
      endif
c
      return
      end
ccc
      subroutine ludcmp_r8(a,n,np,indx,d)
c
c: Modified for real*8 matrices and vectors.
c: To solve A x = b the first time:
c: call ludcmp(a,n,np,indx,d)  [a destroyed here]
c: call lubksb(a,n,np,indx,b)  [b replaced by answer x]
c: Answer x will be returned in b.
c: To solve A x = b for the same A, but different b, just call
c: call lubksb(a,n,np,indx,b) 
c: with a and indx as computed by ludcmp the first time. Again, b is replaced
c: by the answer x.
c
      INTEGER n,np,indx(n),NMAX
      real*8 a(np,np),sum,cdum,TINY
      REAL*8 d
      PARAMETER (NMAX=4,TINY=(1.0d-20,0.d0))
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (dabs(a(i,j)).gt.aamax) aamax=dabs(a(i,j))
11      continue
c: Take out check for zeros:
cc      if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*dabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            cdum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=cdum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j) .eq. 0.d0) a(j,j)=TINY
        if(j.ne.n)then
          cdum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*cdum
18        continue
        endif
19    continue
      return
      END
ccc
      subroutine lubksb_r8(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      real*8 a(np,np),b(n),sum
      INTEGER i,ii,j,ll
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
ccc
      subroutine ludcmp(a,n,np,indx,d)
c
c: Modified for complex*16 matrices and vectors.
c: To solve A x = b the first time:
c: call ludcmp(a,n,np,indx,d)  [a destroyed here]
c: call lubksb(a,n,np,indx,b)  [b replaced by answer x]
c: Answer x will be returned in b.
c: To solve A x = b for the same A, but different b, just call
c: call lubksb(a,n,np,indx,b) 
c: with a and indx as computed by ludcmp the first time. Again, b is replaced
c: by the answer x.
c
      INTEGER n,np,indx(n),NMAX
      complex*16 a(np,np),sum,cdum,TINY
      REAL*8 d
      PARAMETER (NMAX=4,TINY=(1.0d-20,0.d0))
      INTEGER i,imax,j,k
      REAL*8 aamax,dum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (cdabs(a(i,j)).gt.aamax) aamax=cdabs(a(i,j))
11      continue
c: Take out check for zeros:
cc      if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*cdabs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            cdum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=cdum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j) .eq. dcmplx(0.d0,0.d0)) a(j,j)=TINY
        if(j.ne.n)then
          cdum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*cdum
18        continue
        endif
19    continue
      return
      END
ccc
      subroutine lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      complex*16 a(np,np),b(n),sum
      INTEGER i,ii,j,ll
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
ccc
      SUBROUTINE zroots(a,m,roots,polish)
      INTEGER m,MAXM
      REAL EPS
      COMPLEX a(m+1),roots(m)
      LOGICAL polish
      PARAMETER (EPS=1.e-6,MAXM=101)
CU    USES laguer
      INTEGER i,j,jj,its
      COMPLEX ad(MAXM),x,b,c
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      do 13 j=m,1,-1
        x=cmplx(0.,0.)
        call laguer(ad,j,x,its)
        if(abs(aimag(x)).le.2.*EPS**2*abs(real(x))) x=cmplx(real(x),0.)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
13    continue
      if (polish) then
        do 14 j=1,m
          call laguer(a,m,roots(j),its)
14      continue
      endif
      do 16 j=2,m
        x=roots(j)
        do 15 i=j-1,1,-1
          if(real(roots(i)).le.real(x))goto 10
          roots(i+1)=roots(i)
15      continue
        i=0
10      roots(i+1)=x
16    continue
      return
      END
ccc
      SUBROUTINE laguer(a,m,x,its)
      INTEGER m,its,MAXIT,MR,MT
      REAL EPSS
      COMPLEX a(m+1),x
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      REAL abx,abp,abm,err,frac(MR)
      COMPLEX dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=cmplx(0.,0.)
        f=cmplx(0.,0.)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
11      continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=sqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.) then
            dx=m/gp
          else
            dx=exp(cmplx(log(1.+abx),real(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1)return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      pause 'too many iterations in laguer'
      return
      END
ccc
      SUBROUTINE zroot8_int(a,m,xroot,xlo,xhi)
      implicit none
      INTEGER m,MAXM
      REAL*8 EPS,xlo,xhi
      COMPLEX*16 a(m+1),xroot,roots(20)
      PARAMETER (EPS=1.d-6,MAXM=101)
CU    USES laguer8
      INTEGER j,jj,its,nf
      COMPLEX*16 ad(MAXM),x,b,c
      do 11 j=1,m+1
        ad(j)=a(j)
11    continue
      nf=0
      x=xroot
      do 13 j=m,1,-1
cc      x=cmplx(0.,0.)
        call laguer8(ad,j,x,its)
c: Check if root found is in interval and essentially real:
        if(dble(x) .ge. xlo .and. dble(x) .le. xhi .and.
     .     dabs(dimag(x)) .lt. EPS) then
c: Polish if necessary
           if(j .lt. m) then
              call laguer8(a,m,x,its)
           endif
           xroot=dcmplx(dble(x),0.d0)
      nf=nf + 1
           if(j .lt. m) print *,'root not first: ',m-j+1
c          return
        endif
cc      if(abs(dimag(x)).le.2.*EPS**2*abs(dble(x))) 
cc   .     x=dcmplx(real(x),0.d0)
        roots(j)=x
        b=ad(j+1)
        do 12 jj=j,1,-1
          c=ad(jj)
          ad(jj)=b
          b=x*b+c
12      continue
        x=dcmplx(0.5d0,0.d0)
13    continue
cc    print *,'no roots found between xlo and xhi: ',xlo,xhi,roots
      if(nf .gt. 1) print *,'>1 root found: ',roots,xroot
      return
      END
ccc
      SUBROUTINE laguer8(a,m,x,its)
      INTEGER m,its,MAXIT,MR,MT
      REAL*8 EPSS
      COMPLEX*16 a(m+1),x
      PARAMETER (EPSS=2.e-7,MR=8,MT=10,MAXIT=MT*MR)
      INTEGER iter,j
      REAL*8 abx,abp,abm,err,frac(MR)
      COMPLEX*16 dx,x1,b,d,f,g,h,sq,gp,gm,g2
      SAVE frac
      DATA frac /.5,.25,.75,.13,.38,.62,.88,1./
      do 12 iter=1,MAXIT
        its=iter
        b=a(m+1)
        err=abs(b)
        d=dcmplx(0.d0,0.d0)
        f=dcmplx(0.d0,0.d0)
        abx=abs(x)
        do 11 j=m,1,-1
          f=x*f+d
          d=x*d+b
          b=x*b+a(j)
          err=abs(b)+abx*err
11      continue
        err=EPSS*err
        if(abs(b).le.err) then
          return
        else
          g=d/b
          g2=g*g
          h=g2-2.*f/b
          sq=cdsqrt((m-1)*(m*h-g2))
          gp=g+sq
          gm=g-sq
          abp=abs(gp)
          abm=abs(gm)
          if(abp.lt.abm) gp=gm
          if (max(abp,abm).gt.0.) then
            dx=m/gp
          else
            dx=cdexp(dcmplx(dlog(1.d0+abx),dble(iter)))
          endif
        endif
        x1=x-dx
        if(x.eq.x1) return
        if (mod(iter,MT).ne.0) then
          x=x1
        else
          x=x-dx*frac(iter/MT)
        endif
12    continue
      pause 'too many iterations in laguer'
      return
      END
ccc
      subroutine cdhankel(z,tol,H0)
c
c: Computes the zero'th order hankel function of the first kind of the
c: complex argument z using a tolerance of tol.
c
      implicit none
      complex*16 z,H0,zsq,J0,Y0,Y0sum,term,termx
      real*8 tol,gamma,ser,tw_o_pie
      integer*4 j,k,sg
      data gamma/.5772156659015d0/,tw_o_pie/0.63661977236758d0/
c
      zsq=z*z
      k=2 
      ser=1.d0
      term=zsq/4.d0
      Y0sum=term
      sg=1
      J0=1.d0 - term
c
      do j=2,100
         k=k+2
         sg=-sg
         ser=ser + 1.d0/dfloat(j)
         term=term*zsq/dfloat(k*k)
         termx=term*ser
         Y0sum=Y0sum + sg*termx
         J0=J0 - sg*term
         if(cdabs(termx/Y0sum) .lt. tol) then
            Y0=tw_o_pie*((cdlog(0.5d0*z) + gamma)*J0 + Y0sum)
            H0=J0 + dcmplx(-dimag(Y0),dble(Y0))
cc          print *,'j = ',j,z,H0
            return
         endif
      enddo
      print *,'cdhankel failed to converge: ',z,termx,Y0sum
c
      return
      end
ccc
      function magsq(z)
c
      implicit none
      complex*16 z
      real*8 magsq
c
      magsq=dble(z)*dble(z) + dimag(z)*dimag(z)
c
      return
      end
ccc
      function magsq_c8(z)
c
      implicit none
      complex*8 z
      real*4 magsq_c8
c
      magsq_c8=real(z)*real(z) + imag(z)*imag(z)
c
      return
      end
ccc
      function dot(z1,z2)
c
      implicit none
      complex*16 z1,z2
      real*8 dot
      dot=dble(z1)*dble(z2) + dimag(z1)*dimag(z2)
c
      return
      end
ccc
      subroutine hunt(xx,n,x,jlo)
      integer*4 jlo,n
      real*8 x,xx(n)
      integer*4 inc,jhi,jm
      logical ascnd
c
c: EKW FIX (always ascending for us, but tricked when n=1):
cc    ascnd=xx(n) .gt. xx(1)
      ascnd=xx(n) .ge. xx(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x .ge. xx(jlo) .eqv. ascnd) then
1        jhi=jlo+inc
         if(jhi .gt. n)then
            jhi=n+1
         elseif(x .ge. xx(jhi) .eqv. ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=jlo
2        jlo=jhi-inc
         if(jlo .lt. 1) then
            jlo=0
         elseif(x .lt. xx(jlo) .eqv. ascnd) then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
3     if(jhi-jlo .eq. 1) return
      jm=(jhi+jlo)/2
      if(x .gt. xx(jm) .eqv. ascnd) then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
c
      end
ccc
      subroutine hunt_re_im(xx_re,n,x,jlo)
      integer*4 jlo,n
      real*8 x,xx_re(n)
      integer*4 inc,jhi,jm
      logical ascnd
c
      ascnd=xx_re(n) .gt. xx_re(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x .ge. xx_re(jlo) .eqv. ascnd) then
1        jhi=jlo+inc
         if(jhi .gt. n)then
            jhi=n+1
         elseif(x .ge. xx_re(jhi) .eqv. ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=jlo
2        jlo=jhi-inc
         if(jlo .lt. 1) then
            jlo=0
         elseif(x .lt. xx_re(jlo) .eqv. ascnd) then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
3     if(jhi-jlo .eq. 1) return
      jm=(jhi+jlo)/2
      if(x .gt. xx_re(jm) .eqv. ascnd) then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
c
      end
      subroutine traj_mag(mag_des0,dph_left0,klast,lnlast,dlnlast,
     .   k,r1r2,emagmax0,ephmax0,iilhp)
c
c: Moves along line of constant magnitude of the function ln(r1*r2) in the
c: complex k-plane.
c: EKW FIX areas done on 4/25/97 to make routine more robust.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iilhp,ntry,ntryx,iimode1,n10
      complex*16 klast,lnlast,dlnlast,lnrrdes,k,delk,r1r2(3,4),
     .   lnr1r2,dlnr1r2,branch,kp,lnp,dlnp
      real*8 mag_des0,dph_left0,emagmax0,ephmax0,ephmax,
     .   mag_des,ph_des,dmag_des,dph_des,ph_err,ph_step_max,dph_diff,
     .   mag_err,eph_rat,mag_rat,krat,dph_left,mag_fac,delph
c
      npt=0
      mag_fac=1.d0
      dph_left=dph_left0
      ph_step_max=twpie/phfac0
c
      iimode1=0
      krat=dreal(klast)/dreal(xkref)
      if(dabs(krat - 1.d0) .lt. 1.d-6) iimode1=1
c
      ntryx=0
      n10=0
10    dph_des=dph_left/max(1,nint(dabs(dph_left)/ph_step))
      delph=dabs(dph_left)/ph_step
      if(delph .lt. 1000000.d0) then
         dph_des=dph_left/max(1,nint(delph))
      else
         dph_des=dph_left/delph
         if(iiwrite .gt. 0)
     .    print *,'delph used: ',dreal(k)/kw0,dimag(k)*8685.9d0
      endif
      mag_fac=min(2.d0*mag_fac,1.d0)
      dmag_des=mag_fac*(mag_des0 - dreal(lnlast))
      ntry=0
20    ph_des=dimag(lnlast) + dph_des
      mag_des=dreal(lnlast) + dmag_des
      lnrrdes=dcmplx(mag_des,ph_des)
c: Compute next guess for k based on last value and first derivative:
      if(ntry .le. 3 .and. dabs(dph_des) .gt. 0.1257d0) then
c: Compute next guess for k based on last value and first derivative:
         if(dlnlast.eq.0.0) then
            jjfail=0
            return
         end if
         delk=(lnrrdes-lnlast)/dlnlast
      else
cc       print *,'cubic: ',dreal(klast)/kw0,dimag(klast)*8685.9
         if(klast .eq. k) then
            call lnk_cub(klast,kp,lnlast,lnp,dlnlast,dlnp,
     .         lnrrdes,delk)
         else
            call lnk_cub(klast,k,lnlast,lnr1r2,dlnlast,dlnr1r2,
     .         lnrrdes,delk)
         endif
      endif
      k=klast + delk
      ntry=ntry + 1
      ntryx=ntryx + 1
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnr1r2=r1r2(1,4)
      dlnr1r2=r1r2(2,4)
c
      if(ntry .gt. 12) then
         if(iiwrite .gt. 0) then
         print *,'ntry>10: ,ntry,k,k/kw0,iish = ',ntry,k,k/kw0
         print *,'k,klast,delk = ',k/kw0,klast/kw0,delk/kw0
         print *,'lnrr,lnlast,lnrrdes = ',lnr1r2,lnlast,lnrrdes
         print *,'xkbp,iish = ',xkbp(1,1)/kw0,iish
         print *,'xkref,iish_ref = ',xkref/kw0,iish_ref
         end if
         iifail=1
cc          call stop_run
         return
         print *,'num deriv = ',(lnr1r2-lnlast)/(k-klast),dlnlast,
     .      dlnr1r2
cc       call r1r2_calc(k,rr,2,0)
cc       call r1r2_calc(klast,r1r2,2,0)
cc       print *,'num deriv: R1 = ',(rr(1,1)-r1r2(1,1))/(k-klast),
cc   .      rr(2,1),r1r2(2,1)
cc       print *,'num deriv: R2 = ',(rr(1,2)-r1r2(1,2))/(k-klast),
cc   .      rr(2,2),r1r2(2,2)
      endif
      mag_err=mag_des - dreal(lnr1r2)
      mag_rat=dabs(mag_err/emagmax0)
c: Check for bad magnitude (off |R1R2|=1 contour too far):
      if(mag_rat .gt. 1.d0) then
         dph_des=dph_des/2.d0
c: Decrease jump in magnitude as well in case this is what is holding
c: us back:
         dmag_des=dmag_des/2.d0
         ph_step=ph_step/2.d0
         mag_fac=mag_fac/2.d0
         if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
         call xkh_backup
         goto 20
      endif
c
      ph_err=dimag(lnr1r2) - ph_des
      if(dabs(ph_err) .gt. pie) ph_err=ph_err - dsign(twpie,ph_err)
      dph_diff=dph_des + ph_err
c
c: Make sure max phase error ephmax is less than desired phase jump:
      ephmax=min(ephmax0,.75d0*dabs(dph_des))
c: Check for bad phase:
c: EKW FIX: take abs value here:
c     eph_rat=ph_err/ephmax
      eph_rat=dabs(ph_err/ephmax)
      if(eph_rat .gt. 1.d0) then
c: Reduce phase step if phase was not close enough:
         dph_des=dph_des/2.d0
         dmag_des=dmag_des/2.d0
         ph_step=ph_step/2.d0
         mag_fac=mag_fac/2.d0
c: Accept guess if it's mode one and phase jump was too small:
         if(iimode1 .eq. 1 .and. n10 .eq. 0) then
cap      if(iimode1 .eq. 1 .and. dph_diff .gt. 0.d0 .and. 
cap  .      ph_err .lt. 0.d0) then
         else
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
            call xkh_backup
            if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
            goto 20
         endif
      else
         if(eph_rat .lt. 0.25d0 .and. mag_rat .lt. 0.25d0) then
            ph_step=dmin1(2.0d0*ph_step,ph_step_max)
         elseif(eph_rat .lt. 0.5d0 .and. mag_rat .lt. 0.5d0) then
            ph_step=dmin1(1.5d0*ph_step,ph_step_max)
         endif
      endif
c
      dph_left=dph_left - dph_diff
c
      kp=klast
      lnp=lnlast
      dlnp=dlnlast
c
      klast=k
      lnlast=lnr1r2
      dlnlast=dlnr1r2
c
c: Save values on contour for use in contour_find:
      npt=npt+1
      k_cont(npt)=klast
      ln_cont(npt)=lnlast
      dln_cont(npt)=dlnlast
      if(iimt .ne. 0) call mode_traj(k,r1r2,0)
c
cxx   if(dimag(k) .lt. 0.d0 .and. dreal(delk) .gt. 0.d0 .and. 
cxx  .   dimag(delk) .lt. 0.d0) then
c: FIX 2-2-95: Don't require path to heading to right in LHP:
      if(dimag(k) .lt. 0.d0 .and. dimag(delk) .lt. 0.d0) then
         call sheet_look(0,branch)
         if(branch .ne. (0.d0,0.d0)) then
            if(iidiag .ge. 2) print *,'PATH HEADING TOWARD LOWER HP: ',
     .         nmode,k,k/kw
            iilhp=1
            if(iimt .ne. 0) call mode_traj(k,r1r2,-2)
            return
         endif
      endif
cc    if(dreal(k) .lt. kremin .and. real(delk) .lt. 0.d0) then
      if(dreal(k) .lt. kremin .and. iiccw .gt. 0) then
         iidone=1
cc       if(iiwrite .eq. 1) then
cc          print *,'Informative message: Mode search ',
cc   .         'terminating due to small Re(k) (see cphmax).'
cc          write(lusvp,'(a)') 'Informative message: Mode search '//
cc   .         'terminating due to small Re(k) (see cphmax).'
cc       endif
         return
      endif
c: Check if close enough to call eig_final:
      if(dabs(dph_left) .gt. ephmax) then
         if(ntryx .gt. 100) then
            if(iiwrite .gt. 0)
     .       print *,'ntryx>100'
         endif
         if(ntryx .gt. 250) then
            print *,'Total # tries exceeded 250 in traj_mag: ',
     .         k/kw0,iish,dph_left
            iifail=1
cc          call stop_run
            return
         endif
         n10=1
         goto 10
      endif
c
      return
      end
ccc
      subroutine traj_phase(ph_des0,dmag_left0,klast,lnlast,
     .   dlnlast,k,r1r2,emagmax0,ephmax0)
c
c: Moves along line of constant phase of the function ln(r1*r2) in the
c: complex k-plane.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 ntry
      complex*16 klast,lnlast,dlnlast,lnrrdes,k,delk,r1r2(3,4),lnr1r2
      real*8 ph_des0,dmag_left0,emagmax0,ephmax0,emagmax,
     .   dmag_des,mag_des,ph_des,dph_des,mag_err,dmag_diff,
     .   eph_rat,emag_rat,mag_step_max,ph_err,dmag_left,ph_fac
c
      ph_fac=1.d0
      dmag_left=dmag_left0
      mag_step_max=1.d0
 10   if(mag_step.eq.0.d0) then
         jjfail=1
         return
      end if
      dmag_des=dmag_left/max(1,nint(dabs(dmag_left)/mag_step))
      ph_fac=min(2.d0*ph_fac,1.d0)
      dph_des=ph_fac*(ph_des0 - dimag(lnlast))
      ntry=0
      ph_des=ph_des0
20    mag_des=dreal(lnlast) + dmag_des
      ph_des=dimag(lnlast) + dph_des
      lnrrdes=dcmplx(mag_des,ph_des)
c
      if(ntry .le. 3) then
c: Compute next guess for k based on last value and first derivative:
         delk=(lnrrdes-lnlast)/dlnlast
      else
c PLN 070600
         if(klast .eq. k) then
            iifail = 1
            return
         else
            call lnk_cub(klast,k,lnlast,lnr1r2,dlnlast,r1r2(2,4),
     .           lnrrdes,delk)
         end if
cc       print *,'cubic fit in traj_phase: ',klast,k,lnlast,lnr1r2,
cc   .      dlnlast,r1r2(2,4),delk,(lnrrdes-lnlast)/dlnlast
      endif
      k=klast + delk
      ntry=ntry + 1
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
c
      if(ntry .gt. 10) then
         if(iiwrite .gt. 0) then
         print *,'ntry>10 in traj_phase: ,ntry,k,k/kw,iish = ',
     .      ntry,k,k/kw,iish
         print *,'k,klast,delk = ',k/kw,klast/kw,delk/kw,lnrrdes,
     .      lnlast,dlnlast
         print *,'xkbp,iish = ',xkbp(1,1)/kw,iish
         end if
         iifail = 1
         return
c         print *,'num deriv = ',(r1r2(1,4)-lnlast)/(k-klast),dlnlast,
c     .      r1r2(2,4)

         if(ntry .gt. 11) then
            iifail=1
cc          call stop_run
            return
         endif
      endif
c
      lnr1r2=r1r2(1,4)
      ph_err=ph_des - dimag(lnr1r2)
      eph_rat=dabs(ph_err/ephmax0)
c: Check for bad phase:
      if(eph_rat .gt. 1.d0) then
         dmag_des=dmag_des/2.d0
         dph_des=dph_des/2.d0
         mag_step=mag_step/2.d0
         ph_fac=ph_fac/2.d0
         if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
         call xkh_backup
         goto 20
      endif
c
c: Make sure max mag error emagmax is less than desired mag jump:
      emagmax=min(emagmax0,.75d0*dabs(dmag_des))
c: Check for phase overshoot:
      mag_err=dabs(mag_des - dreal(lnr1r2))
c: Don't allow phase to be off in other direction either:
      emag_rat=mag_err/emagmax
      if(emag_rat .gt. 1.d0) then
         dmag_des=dmag_des/2.d0
         dph_des=dph_des/2.d0
         mag_step=mag_step/2.d0
         ph_fac=ph_fac/2.d0
         if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
         call xkh_backup
         goto 20
      else
         if(emag_rat .lt. 0.25d0 .and. eph_rat .lt. 0.25d0) then
            mag_step=dmin1(2.0d0*mag_step,mag_step_max)
         elseif(emag_rat .lt. 0.5d0 .and. eph_rat .lt. 0.5d0) then
            mag_step=dmin1(1.5d0*mag_step,mag_step_max)
         endif
      endif
c
      dmag_diff=dreal(lnr1r2) - dreal(lnlast)
c: Update how much phase left to go:
      dmag_left=dmag_left - dmag_diff
c
      klast=k
      lnlast=lnr1r2
      dlnlast=r1r2(2,4)
c
      if(iimt .ne. 0) call mode_traj(k,r1r2,0)
c
c: Check if close enough to call eig_final:
      if(dabs(dmag_left) .gt. emagmax0) goto 10
c
      return
      end
ccc
      subroutine traj_hor(k,r1r2,emagmax0,mag_step0,ii_ext)
c
c: Moves along horizontal line in complex k-plane until |R1R2|=1 
c: contour is crossed.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 ii_ext,ntry,ntryx,iisg0,iisg
      complex*16 k,r1r2(3,4)
      real*8 emagmax0,mag_step0,emagmax,dmag_des,mag_des,mag_err,
     .   dmag_diff,emag_rat,mag_step_max,dmag_left,ln_re,dln_re,
     .   ki,delk_re,klast_re,lnlast_re,dlnlast_re
      data mag_step_max/2.d0/
c
      ii_ext=0
      mag_step=mag_step0
      ki=dimag(k)
      klast_re=dreal(k)
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnlast_re=dreal(r1r2(1,4))
      dlnlast_re=dreal(r1r2(2,4))
      iisg0=nint(sign(1.d0,dlnlast_re))
      dmag_left=-lnlast_re
      ntryx=0
10    dmag_des=dmag_left/max(1,nint(dabs(dmag_left)/mag_step))
      ntry=0
20    mag_des=lnlast_re + dmag_des
c
      delk_re=dmag_des/dlnlast_re
      k=dcmplx(klast_re+delk_re,ki)
      ntry=ntry + 1
      ntryx=ntryx + 1
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      ln_re=dreal(r1r2(1,4))
      dln_re=dreal(r1r2(2,4))
      iisg=nint(sign(1.d0,dln_re))
      if(iisg*iisg0 .lt. 0) then
         print *,'ii_ext=1 in traj_hor: ',k/kw,ntryx
         ii_ext=1
         return
      endif
c
      if(ntry .gt. 10) then
         iifail=1
cc       call stop_run
      print *,'ntryx = ',ntryx
         return
      endif
c
c: Make sure max mag error emagmax is less than desired mag jump:
      emagmax=min(emagmax0,.75d0*dabs(dmag_des))
      mag_err=dabs(mag_des - ln_re)
c: Don't allow phase to be off in other direction either:
      emag_rat=mag_err/emagmax
      if(emag_rat .gt. 1.d0) then
         dmag_des=dmag_des/2.d0
         mag_step=mag_step/2.d0
         if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
         call xkh_backup
         goto 20
      else
         if(emag_rat .lt. 0.25d0) then
            mag_step=dmin1(2.0d0*mag_step,mag_step_max)
         elseif(emag_rat .lt. 0.5d0) then
            mag_step=dmin1(1.5d0*mag_step,mag_step_max)
         endif
      endif
c
      dmag_diff=ln_re - lnlast_re
c: Update how much phase left to go:
      dmag_left=dmag_left - dmag_diff
c
      klast_re=dreal(k)
      lnlast_re=ln_re
      dlnlast_re=dln_re
c
      if(iimt .ne. 0) call mode_traj(k,r1r2,0)
c
c: Check if close enough to call eig_final:
      if(dabs(dmag_left) .gt. emagmax0) goto 10
c
      print *,'ntryx = ',ntryx
      return
      end
ccc
      subroutine traj_ver(k,r1r2,emagmax0,mag_step0,ii_ext)
c
c: Moves along vertical line in complex k-plane until |R1R2|=1 
c: contour is crossed.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 ii_ext,ntry,ntryx,iisg0,iisg
      complex*16 k,r1r2(3,4)
      real*8 emagmax0,mag_step0,emagmax,dmag_des,mag_des,mag_err,
     .   dmag_diff,emag_rat,mag_step_max,dmag_left,ln_re,dln_im,
     .   kr,delk_im,klast_im,lnlast_re,dlnlast_im
      data mag_step_max/2.d0/
c
      ii_ext=0
      mag_step=mag_step0
      kr=dreal(k)
      klast_im=dimag(k)
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      lnlast_re=dreal(r1r2(1,4))
      dlnlast_im=dimag(r1r2(2,4))
      iisg0=nint(sign(1.d0,dlnlast_im))
      dmag_left=-lnlast_re
      ntryx=0
10    dmag_des=dmag_left/max(1,nint(dabs(dmag_left)/mag_step))
      ntry=0
20    mag_des=lnlast_re + dmag_des
c
      delk_im=-dmag_des/dlnlast_im
      k=dcmplx(kr,klast_im+delk_im)
      ntry=ntry + 1
      ntryx=ntryx + 1
      call r1r2_calc(k,r1r2,2,0,jjfail)
      if(jjfail .gt. 0) return
      ln_re=dreal(r1r2(1,4))
      dln_im=dimag(r1r2(2,4))
      iisg=nint(sign(1.d0,dln_im))
      if(iisg*iisg0 .lt. 0) then
         print *,'ii_ext=1 in traj_ver: ',k/kw,ntryx
         ii_ext=1
         return
      endif
c
      if(ntry .gt. 10) then
         iifail=1
cc       call stop_run
      print *,'ntryx = ',ntryx
         return
      endif
c
c: Make sure max mag error emagmax is less than desired mag jump:
      emagmax=min(emagmax0,.75d0*dabs(dmag_des))
      mag_err=dabs(mag_des - ln_re)
c: Don't allow phase to be off in other direction either:
      emag_rat=mag_err/emagmax
      if(emag_rat .gt. 1.d0) then
         dmag_des=dmag_des/2.d0
         mag_step=mag_step/2.d0
         if(iimt .ne. 0) call mode_traj(k,r1r2,-1)
c: For purposes of sheet changing, go back to previous point so that
c: cut is not crossed in a "triangle":
         call xkh_backup
         goto 20
      else
         if(emag_rat .lt. 0.25d0) then
            mag_step=dmin1(2.0d0*mag_step,mag_step_max)
         elseif(emag_rat .lt. 0.5d0) then
            mag_step=dmin1(1.5d0*mag_step,mag_step_max)
         endif
      endif
c
      dmag_diff=ln_re - lnlast_re
c: Update how much phase left to go:
      dmag_left=dmag_left - dmag_diff
c
      klast_im=dimag(k)
      lnlast_re=ln_re
      dlnlast_im=dln_im
c
      if(iimt .ne. 0) call mode_traj(k,r1r2,0)
c
c: Check if close enough to call eig_final:
      if(dabs(dmag_left) .gt. emagmax0) goto 10
c
      print *,'ntryx = ',ntryx
      return
      end
ccc
      subroutine lnk_cub(k1,k2,lnk1,lnk2,lnk1p,lnk2p,lndes,delk)
c
c: Fits a cubic to the function e(k), given two points (k1,lnk1) and (k2,ek2)
c: and the derivatives lnk1p and lnk2p.  Solves the cubic to find the point
c: kx where e(k1+delk)=0.
c
      implicit none
      complex*16 k1,k2,lnk1,lnk2,lnk1p,lnk2p,lndes,delk,c1,c2,c3,c4,
     .   dk,kx
c
c: Fit cubic:
      call lnk_cub_fit(k1,k2,lnk1,lnk2,lnk1p,lnk2p,c1,c2,c3,c4,dk)
      call lnk_cub_root(c1,c2,c3,c4,lndes,kx)
      delk=kx*dk
c
      return
      end
ccc
      subroutine lnk_cub_fit(k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,dk)
c
c: This subroutine fits a cubic polynomial y(x)=c1 + c2*x + c3*x**2 + 
c: c4*x**3, where x=(k-k1)/(k2-k1), to the points (k1,y1) and (k2,y2) 
c: and the derivatives y1p and y2p at those points.
c
      implicit none
      complex*16 k1,k2,y1,y2,y1p,y2p,c1,c2,c3,c4,b1,b2,dk,y1px,y2px
c
      dk=k2 - k1
c: Convert derivatives from d/dk to d/dx by mult by dk/dx:
      y1px=y1p*dk
      y2px=y2p*dk
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
      subroutine lnk_cub_root(c1,c2,c3,c4,ydes,xhit)
c
c: Finds the roots, x1,x2,x3, of the complex-valued
c: polynomial y(x)-ydes=0: c1-ydes + c2*x + c3*x^2 + c4*x^3 = 0.
c
      implicit none
      complex*16 c1,c2,c3,c4,ydes,xhit,x1,x2,x3,a1,a2,a3,a1sq,Q,R,A,B,
     .   dd,rdd,a1_3,AB_sum,ifac,c1x,x2_term1,x2_term2
      real*8 one_third,sq3_2,x1magsq,x2magsq,x3magsq,re_r_rdd
      data one_third/0.33333333333333d0/,sq3_2/0.86602540378444d0/
c
      c1x=c1 - ydes
      a1=c3/c4
      a2=c2/c4
      a3=c1x/c4
c: Use Numerical Recipes, p. 179:
      a1sq=a1*a1
      Q=(a1sq - 3.d0*a2)/9.d0
      R=(2.d0*a1*a1sq  - 9.d0*a1*a2 + 27.d0*a3)/54.d0
      dd=R*R - Q*Q*Q
      rdd=cdsqrt(dd)
      re_r_rdd=dreal(r)*dreal(rdd) + dimag(r)*dimag(rdd)
      if(re_r_rdd .lt. 0.d0) rdd=-rdd
c
      A=-(r + rdd)**one_third
      if(A .ne. dcmplx(0.d0,0.d0)) then
         B=Q/A
      else
         B=dcmplx(0.d0,0.d0)
      endif
c
      AB_sum=A + B
      a1_3=one_third*a1
      x1=AB_sum - a1_3
cc    print *,'f(x1): ',c1x + c2*x1 + c3*x1**2 + c4*x1**3.,c4
      x2_term1=-0.5d0*AB_sum - a1_3
      ifac=sq3_2*(A - B)
      x2_term2=dcmplx(-dimag(ifac),dreal(ifac))
      x2=x2_term1 + x2_term2
      x3=x2_term1 - x2_term2
cc    print *,'f(x2): ',c1x + c2*x2 + c3*x2**2 + c4*x2**3.
cc    print *,'f(x3): ',c1x + c2*x3 + c3*x3**2 + c4*x3**3.
c
      x1magsq=dreal(x1)*dreal(x1) + dimag(x1)*dimag(x1)
      x2magsq=dreal(x2)*dreal(x2) + dimag(x2)*dimag(x2)
      x3magsq=dreal(x3)*dreal(x3) + dimag(x3)*dimag(x3)
      if(x1magsq .lt. x2magsq) then
         if(x1magsq .lt. x3magsq) then
            xhit=x1
         else
            xhit=x3
         endif
      else
         if(x2magsq .lt. x3magsq) then
            xhit=x2
         else
            xhit=x3
         endif
      endif
c
      return
      end 
      subroutine vbl_init
c
c: Initializes variables in scairy_com and lab_com.
c
      implicit none
      include 'scairy_com'
      include 'lab_com'
c
      ei13=( 0.5D0, 0.8660254037844387D0)
      ei23=(-0.5D0, 0.8660254037844387D0)
      ei16=( 0.8660254037844387D0, 0.5D0)
      ei56=(-0.8660254037844387D0, 0.5D0)
      eim13=( 0.5D0, -0.8660254037844387D0)
      eim23=(-0.5D0, -0.8660254037844387D0)
      eim16=( 0.8660254037844387D0, -0.5D0)
      eim56=(-0.8660254037844387D0, -0.5D0)
      eye=(0.D0,1.D0)
cc    zpt=(-0.9D0,2.8D0)
      pie23=2.094395102393195D0
      pie_x=3.1415926535897932D0
      pie_inv=0.31830988618379D0
      det_bi=(0.31830988618379,0.)
      det_pos=(0.13783222385545,0.07957747154595)
      det_neg=(0.13783222385545,-0.07957747154595)
      sqrt3=1.73205080756888
cc
      flab='Frequency  - Hz'
      rlab='Range  - km'
      mnlab='Mode Number'
      tlab='Time  - s'
      krlab='RE[k]'
      kilab='IM[k]'
      thlab='Grazing Angle  - deg'
      mlab='|R| - dB'
      phlab='Phase of R  - deg'
      dlab='Depth  - m'
      pllab='Propagation Loss - dB'
      mrlab='RE[Mode Amplitude]'
      milab='IM[Mode Amplitude]'
      malab='Mode Amplitude'
      mplab='Mode Phase'
      dblab='dB'
      mtlab='m'
      kmlab='km'
      z4=0. 
c
      return
      end
      subroutine xkh_init(ndv)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer ndv,j,ii,jj,jp1,jm1,ii2
      complex*16 ik,xbx
c
c: Keep track of previous ratio in case we need to back up in eig_find:
      xkhratp=xkhrat
      call iish_xfer(iish,iish_ref,iish0,iishr0)
c: Variable used to check for sheet changes (use kw0, which does not change
c: when the duct changes):
      xkhrat=xkh/kw0
      xkhsq=xkh*xkh
c
      ik=dcmplx(-dimag(xkh),dreal(xkh))
      ikcon(1)=ik
      ikcon(2)=(0.d0,1.d0)
      ikcon(3)=(0.d0,0.d0)
c
      if(isp(nlay) .eq. 0) then
c: For Airy halfspaces, work outwards from reference depth:
c: Reference layer:
         j=nsvmin
         ii=isvmin
         ii2=3-isvmin
         call sheet(xksq(ii,j),xkhsq,xk(ii,j),xkh,w,gami(1,ii,j),
     .      ndv,iich_ref,iish_ref(1),xkhratp,xkhrat,xkrat_ref(1))
         if(isp(j) .eq. 0) then
            call gami_calc(xkh,xkhsq,xk(ii2,j),xksq(ii2,j),w,
     .         gami(1,ii2,j),ndv)
         else
            do jj=1,ndv
               gami(jj,ii2,j)=gami(jj,ii,j)
            enddo
         endif
c pln: Added 140900
         if(abs(xksq(ii,j)).eq.abs(xkhsq)) then
            jjfail=1
            return
         end if
c
c: Above reference layer:
         do j=nsvmin-1,2,-1
            jp1=j+1
            call gami_lay(mm(j),isp(j),xkh,xkhsq,gami(1,2,j),
     .         gami(1,1,jp1),gami(1,1,j),xk(2,j),xksq(2,j),
     .         xk(1,j),xksq(1,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(j),iss(j),xkh,xkhsq,beti(1,2,j),
     .            beti(1,1,jp1),beti(1,1,j),xb(2,j),xbsq(2,j),
     .            xb(1,j),xbsq(1,j),w,ndv)
            endif
         enddo
c: Below reference layer:
         do j=nsvmin+1,nlay-1
            jm1=j-1
            call gami_lay(mm(jm1),isp(j),xkh,xkhsq,gami(1,1,j),
     .         gami(1,2,jm1),gami(1,2,j),xk(1,j),xksq(1,j),
     .         xk(2,j),xksq(2,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(jm1),iss(j),xkh,xkhsq,beti(1,1,j),
     .            beti(1,2,jm1),beti(1,2,j),xb(1,j),xbsq(1,j),
     .            xb(2,j),xbsq(2,j),w,ndv)
            endif
         enddo
c
c: Lower halfspace:
         call gami_hsp(nlay,1,ndv,mm(nlay-1),gami(1,2,nlay-1))
c: Upper halfspace:
         call gami_hsp(1,2,ndv,mm(1),gami(1,1,2))
      else
c: For homogeneous halfspaces, work inwards from halfspaces:
c: Lower halfspace:
         call gami_hsp(nlay,1,ndv,mm(nlay-1),gami(1,2,nlay-1))
c: Upper halfspace:
         call gami_hsp(1,2,ndv,mm(1),gami(1,1,2))
c: Above reference layer:
         do j=2,nsvmin-1
            jm1=j-1
            call gami_lay(mm(jm1),isp(j),xkh,xkhsq,gami(1,1,j),
     .         gami(1,2,jm1),gami(1,2,j),xk(1,j),xksq(1,j),xk(2,j),
     .         xksq(2,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(jm1),iss(j),xkh,xkhsq,beti(1,1,j),
     .            beti(1,2,jm1),beti(1,2,j),xb(1,j),xbsq(1,j),
     .            xb(2,j),xbsq(2,j),w,ndv)
            endif
         enddo
c: Below reference layer:
         do j=nlay-1,nsvmin+1,-1
            jp1=j+1
            call gami_lay(mm(j),isp(j),xkh,xkhsq,gami(1,2,j),
     .         gami(1,1,jp1),gami(1,1,j),xk(2,j),xksq(2,j),
     .         xk(1,j),xksq(1,j),w,ndv)
            if(iisol(j) .eq. 1) then
               call gami_lay(mm(j),iss(j),xkh,xkhsq,beti(1,2,j),
     .            beti(1,1,jp1),beti(1,1,j),xb(2,j),xbsq(2,j),
     .            xb(1,j),xbsq(1,j),w,ndv)
            endif
         enddo
c: Reference layer:
         if(nsvmin .ne. nlay .and. nsvmin .ne. 1) then
            j=nsvmin
            ii=isvmin
            ii2=3-isvmin
            call sheet_ref(xksq(ii,j),xkhsq,xk(ii,j),xkh,w,
     .         gami(1,ii,j),ndv,iich_ref,iish_ref(1),xkhratp,xkhrat,
     .         xkrat_ref(1))
            if(isp(j) .eq. 0) then
               call gami_calc(xkh,xkhsq,xk(ii2,j),xksq(ii2,j),w,
     .            gami(1,ii2,j),ndv)
            else
               do jj=1,ndv
                  gami(jj,ii2,j)=gami(jj,ii,j)
               enddo
            endif
         endif
      endif
c
c: Vertical wavenumber at reference depth:
      gamiref=gami(1,isvmin,nsvmin)
c
      do j=1,jsol(2,1)-1
         if(mm(j) .eq. 1) then
            jp1=j + 1
            call pquv_calc(rhorat(j),w,ik,xbsqinv(2,j),xbsqinv(1,jp1),
     .         xkh,Pcon(1,j),Qcon(1,j),Ucon(1,j),Vcon(1,j),ndv)
         endif
      enddo
      do j=jsol(1,1),nlay-1
         if(mm(j) .eq. 1) then
            jp1=j + 1
            call pquv_calc(rhorat(j),w,ik,xbsqinv(1,jp1),xbsqinv(2,j),
     .         xkh,Pcon(1,j),Qcon(1,j),Ucon(1,j),Vcon(1,j),ndv)
         endif
      enddo
c
      do ii=1,2
c: Alay,Blay(1:2,ii),ikcon are elements and derivatives 
c: of the T^(-1) matrix at the top of the first solid layer 
c: in bottom (ii=1) and top (ii=2) layering.
         if(allf(ii) .eq. 0) then
            xbx=xbsqinv(ii,jsol(ii,1))
            Alay(1,ii)=-1.d0 + 2.d0*xkhsq*xbx
            Blay(1,ii)=-2.d0*ik*xbx
            if(ndv .ge. 2) then
               Alay(2,ii)=4.d0*xkh*xbx
               Blay(2,ii)=(0.,-2.d0)*xbx
               if(ndv .ge. 3) then
                  Alay(3,ii)=-4.d0*xkhsq*xbx/w
                  Blay(3,ii)=-2.d0*Blay(1,ii)/w
               endif
            endif
         endif
      enddo
c
      return
      end
ccc
      subroutine gami_hsp(j,ii,ndv,mmx,gamip)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer j,ii,ndv,mmx,jj
      complex*16 gamip(3)
c
      if(isp(j) .eq. 1) then
c: Homogeneous halfspace:
         call sheet(xksq(ii,j),xkhsq,xk(ii,j),xkh,w,gami(1,ii,j),
     .      ndv,iich,iish(ii,1),xkhratp,xkhrat,xkrat(ii,1))
      else
c: Airy halfspace:
         if(j .ne. nsvmin) then
            if(mmx .eq. 0) then
               do jj=1,ndv
                  gami(jj,ii,j)=gamip(jj)
               enddo
            else
               call gami_calc(xkh,xkhsq,xk(ii,j),xksq(ii,j),w,
     .            gami(1,ii,j),ndv)
            endif
         endif
      endif
      if(iisol(j) .eq. 1) then
         if(iss(j) .eq. 1) then
            call sheet(xbsq(ii,j),xkhsq,xb(ii,j),xkh,w,beti(1,ii,j),
     .         ndv,iich,iish(ii,2),xkhratp,xkhrat,xkrat(ii,2))
         else
            call gami_calc(xkh,xkhsq,xb(ii,j),xbsq(ii,j),w,
     .         beti(1,ii,j),ndv)
         endif
      endif
c
      return
      end
ccc
      subroutine gami_lay(mmx,ispx,xkh,xkhsq,gami1,gamip,gami2,
     .   xk1,xk1sq,xk2,xk2sq,w,ndv)
c
      implicit none
      integer*4 mmx,ispx,ndv,jj
      complex*16 xkh,xkhsq,gami1(3),gamip(3),gami2(3),xk1,xk1sq,
     .   xk2,xk2sq
      real*8 w
c
      if(mmx .eq. 0) then
         do jj=1,ndv
            gami1(jj)=gamip(jj)
         enddo
      else
         call gami_calc(xkh,xkhsq,xk1,xk1sq,w,gami1,ndv)
      endif
c
      if(ispx .eq. 0) then
         call gami_calc(xkh,xkhsq,xk2,xk2sq,w,gami2,ndv)
      else
         do jj=1,ndv
            gami2(jj)=gami1(jj)
         enddo
      endif
c
      return
      end
ccc
      subroutine sheet_init(k0,ii_ones,iish_in,iish_ref_in)
c
      implicit none
      include 'Parms_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'gen_com'
      integer*4 ii_ones,iish_in(4),iish_ref_in(2)
      complex*16 k0
c
      xkhratp=k0/kw0
      xkhrat=k0/kw0
      if(ii_ones .eq. 1) then
         iish(1,1)=1
         iish(1,2)=1
         iish(2,1)=1
         iish(2,2)=1
         iish_ref(1)=1
         iish_ref(2)=1
      else
         call iish_xfer(iish_in,iish_ref_in,iish,iish_ref)
      endif
c
      return
      end
ccc
      subroutine gami_calc(xkh,xkhsq,xk,xksq,omega,gami,ndv)
c
c: Computes vertical wavenumber and its derivative w.r.t k using 
c: the Pekeris cut.
      implicit none
      integer ndv
      complex*16 xkh,xkhsq,xk,xksq,gami(3),gamma,gamsq
      real*8 omega
c
      gamsq=xksq - xkhsq
      gamma=cdsqrt(gamsq)
c: Pekeris cut:
      if(dreal(xkh) .gt. dreal(xk)) then
         if(dimag(gamma) .lt. 0.d0) then
            gamma=-gamma
         endif
      endif
c
      gami(1)=dcmplx(-dimag(gamma),dreal(gamma))
      if(ndv .ge. 2) then
         if(gami(1) .ne. (0.d0,0.d0)) then
            gami(2)=xkh/gami(1)
            if(ndv .ge. 3) gami(3)=-xksq/(omega*gami(1))
         else
            gami(2)=1.d100*xkh
            if(ndv .ge. 3) gami(3)=1.d100*xkh
         endif
      endif
c
      return
      end
ccc
      subroutine sheet(xksq,xkhsq,xk,xkh,omega,gami,ndv,iich,iish,
     .   xkhratp,xkhrat,xkrat)
c
c: Computes vertical wavenumber and its derivative w.r.t k using 
c: Pekeris cut.  Sheets are changed if iich=1 and present and previous
c: k-values straddle the branch cut. Positive sheet always taken when
c: iich=-1
      implicit none
      integer ndv,iich,iish
      complex*16 xksq,xkhsq,xk,xkh,gamma,gami(3),gamsq,
     .   xkhratp,xkhrat,xkrat
      real*8 omega
c
      gamsq=xksq - xkhsq
      gamma=cdsqrt(gamsq)
c: Check if present and previous points straddle branch cut:
      if(iich .eq. 1) then
c: Check sheet changing in cos(theta)=k/kw0 space instead of k-space
c: so that frequency changes won't affect the process:
         call cross_cut_pek(xkhratp,xkhrat,xkrat,iish)
      endif
c: Get on positive sheet of Pekeris branch cut (do EJP only when 
c: k to right of branch point):
      if(dreal(xkhrat) .gt. dreal(xkrat)) then
         if(dimag(gamma) .lt. 0.) then
            gamma=-gamma
         endif
      endif
c
c: Apply current sheet to gamma:
      gamma=iish*gamma
      gami(1)=dcmplx(-dimag(gamma),dreal(gamma))
      if(ndv .ge. 2) then
         if(gami(1) .ne. (0.,0.)) then
            gami(2)=xkh/gami(1)
         else
            gami(2)=1.d100*xkh
         endif
         if(ndv .ge. 3) then
            if(gami(1) .ne. (0.d0,0.d0)) then
               gami(3)=-xksq/(omega*gami(1))
            else
               gami(3)=1.d100*xkh
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine sheet_ref(xksq,xkhsq,xk,xkh,omega,gami,ndv,iich_ref,
     .   iish_ref,xkhratp,xkhrat,xkrat)
c
c: Computes vertical wavenumber and its derivative w.r.t k using 
c: Pekeris cut.  Sheets are changed if iich_ref=1 and present and previous
c: k-values straddle the branch cut.
      implicit none
      integer ndv,iich_ref,iish_ref
      complex*16 xksq,xkhsq,xk,xkh,gamma,gami(3),gamsq,
     .   xkhratp,xkhrat,xkrat
      real*8 omega
c
      gamsq=xksq - xkhsq
      gamma=cdsqrt(gamsq)
c: Check if present and previous points straddle branch cut:
      if(iich_ref .eq. 1) then
         call cross_cut_pek(xkhratp,xkhrat,xkrat,iish_ref)
      endif
c: Get on positive sheet of Pekeris branch cut (do EJP only when 
c: k to right of branch point):
      if(dreal(xkhrat) .gt. dreal(xkrat)) then
         if(dimag(gamma) .lt. 0.d0) then
            gamma=-gamma
         endif
      endif
c
c: Apply current sheet to gamma:
      gamma=iish_ref*gamma
      gami(1)=dcmplx(-dimag(gamma),dreal(gamma))
      if(ndv .ge. 2) then
         if(gami(1) .ne. (0.,0.)) then
            gami(2)=xkh/gami(1)
         else
            gami(2)=1.d100*xkh
         endif
         if(ndv .ge. 3) then
            if(gami(1) .ne. (0.d0,0.d0)) then
               gami(3)=-xksq/(omega*gami(1))
            else
               gami(3)=1.d100*xkh
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine cross_cut_pek(kp,k,kx,iishq)
c
c: Checks to see if the Pekeris branch cut has been crossed.
c: kp,k are previous and current values of xkh, kx is branch point.
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iishq
      complex*16 kp,k,kx,delk
      real*8 ypt
c
      if((dreal(kp) .le. dreal(kx) .and. dreal(k) .gt. dreal(kx)) .or.
     .   (dreal(kp) .gt. dreal(kx) .and. dreal(k) .le. dreal(kx))) then
         delk=k - kp
         ypt=dimag(kp) + (dreal(kx)-dreal(kp))*dimag(delk)/dreal(delk)
         if(ypt .gt. dimag(kx)) then
            iishq=-iishq
            iicut=1
            if(iiccw*dreal(delk) .lt. 0.d0) then
               kcut=k*kw0
            else
               kcut=kp*kw0
            endif
            if(iidiag .ge. 2) then 
               print *,'pek sheet changed: ',-iishq,' to ',iishq,
     .            dreal(kp),dreal(k),dreal(kx),
     .            dimag(kp)*8685.9*kw0,dimag(k)*8685.9*kw0,
     .            dimag(kx)*8685.9*kw0
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine pquv_calc(rhorat,w,ik,xbsqinv1,xbsqinv2,xkh,
     .   P,Q,U,V,ndv)
c
      implicit none
      integer*4 ndv
      complex*16 ik,xbsqinv1,xbsqinv2,xkh,P(3),Q(3),U(3),V(3),F,
     .   i1mp,k_w
      real*8 rhorat,w
c
      F=2.d0*(xbsqinv1 - rhorat*xbsqinv2)
      U(2)=dcmplx(-dimag(F),dreal(F))
      U(1)=xkh*U(2)
      Q(2)=-2.d0*xkh*F
      Q(1)=0.5d0*xkh*Q(2) + 1.d0
      P(2)=-Q(2)
      P(1)=0.5d0*xkh*P(2) + rhorat
      i1mp=dcmplx(0.d0,1.d0)*(1.d0 - P(1))
      V(2)=-ik*P(2) + i1mp
      V(1)=xkh*i1mp
      if(ndv .ge. 3) then
         k_w=-xkh/w
         P(3)=k_w*P(2)
cxx      Q(3)=k_w*Q(2)
         Q(3)=-P(3)
         U(3)=2.d0*k_w*U(2)
         V(3)=-ik*P(3)
      endif
c
      return
      end
ccc
      subroutine iish_xfer(iish_in,iish_ref_in,iish_out,iish_ref_out)
c
      implicit none
      integer*4 iish_in(2,2),iish_out(2,2),iish_ref_in(2),
     .   iish_ref_out(2)
c
      iish_out(1,1)=iish_in(1,1)
      iish_out(2,1)=iish_in(2,1)
      iish_out(1,2)=iish_in(1,2)
      iish_out(2,2)=iish_in(2,2)
      iish_ref_out(1)=iish_ref_in(1)
      iish_ref_out(2)=iish_ref_in(2)
c
      return
      end
ccc
      subroutine iish_code(iish,iish_ref,iicode,ii)
c
c: Codes iish(2,2) into a single integer iicode for storage (ii=1).
c: or decodes iicode into iish(2,2) (ii=-1).
c
      implicit none
      integer*4 iish(4),iish_ref(2),iicode,ii,ipow(4),ipow2(2),j,iic
      data ipow/-1,-2,-4,-8/,ipow2/-16,-32/
c
      if(ii .eq. 1) then
         iicode=min(iish_ref(2),0)*32 + min(iish_ref(1),0)*16 + 
     .      min(iish(4),0)*8 + min(iish(3),0)*4 + 
     .      min(iish(2),0)*2 + min(iish(1),0)
      else
         iish_ref(1)=1
         iish_ref(2)=1
         iish(1)=1
         iish(2)=1
         iish(3)=1
         iish(4)=1
c: Short cut:
         if(iicode .eq. 0) return
         iic=iicode
c: Decode the sheets:
         do j=2,1,-1
            if(iic .le. ipow2(j)) then
               iish_ref(j)=-1
               iic=iic - ipow2(j)
            endif
         enddo
         do j=4,1,-1
            if(iic .le. ipow(j)) then
               iish(j)=-1
               iic=iic - ipow(j)
            endif
         enddo
      endif
c
      return
      end
ccc
      subroutine sheet_look(kkch,branch)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer kkch,ii0
      complex*16 branch
c
      branch=(0.d0,0.d0)
      ii0=0
      if(isp(1) .eq. 1) call sheet_ch(iish(1,1),xkbp(1,1),
     .   branch,kkch,ii0)
      if(iss(1) .eq. 1) call sheet_ch(iish(1,2),xkbp(1,2),
     .   branch,kkch,ii0)
      if(isp(nlay) .eq. 1) call sheet_ch(iish(2,1),xkbp(2,1),
     .   branch,kkch,ii0)
      if(iss(nlay) .eq. 1) call sheet_ch(iish(2,2),xkbp(2,2),
     .   branch,kkch,ii0)
c
      return
      end
ccc
      subroutine sheet_ch(iish,xkbp,branch,kkch,ii0)
c
      implicit none
      integer iish,kkch,ii,ii0,jj
      complex*16 xkbp,branch
c
      if(iish .lt. 0) then
c: If more than one negative sheet, choose cut with branch point farthest
c: to the right:
         if(dreal(xkbp) .gt. dreal(branch)) then
            branch=xkbp
            if(kkch .eq. 1) then
               iish=1
c: If more than on -1 sheet, change left ones back to -1:
               if(ii0 .gt. 0) iish=-1
           ii0=ii
            endif
         endif
      endif
c
      return
      end
ccc
      subroutine xkh_backup
c
c: Called from k-plane search routines when we want to back up to previously
c: called value of k.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      iish(1,1)=iish0(1,1)
      iish(2,1)=iish0(2,1)
      iish(1,2)=iish0(1,2)
      iish(2,2)=iish0(2,2)
      iish_ref(1)=iishr0(1)
      iish_ref(2)=iishr0(2)
      xkhrat=xkhratp
c
      return
      end
      subroutine zmx_init
c
c: Initializes depth arrays for computation of mode functions.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jptr,jj,j,jj1,jj2,jsr,ii,jx1
      real*8 z1,z2,rhfac,fac
      complex*16 c1inv,c1invsq,c2inv,c2invsq,cfac
c: fac=2*pi*1000*20*log(e):
      data fac/5.457505415367364d+04/
c
      jptr=0
      jj=1 
      z1=-1.d+20
      do j=1,jflu(2,1)+jflu(2,5)
         z2=zdep(j)
c: For receiver depths equal to interface depths, put receiver in layer
c: closest to top ocean layer:
         if(j .ge. jsurf) z2=z2 + 1.d-10*abs(z2)
         do while(zsr(jj) .lt. z1 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj1=jj
         do while(zsr(jj) .lt. z2 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj2=jj-1
         jzmx(j)=jptr+1
         nzmx(j)=jj2-jj1+1
         rhfac=(geo(1,3,j) - geo(2,3,j))/dmax1(1.d-20,h(j))
cxx      c1inv=1./geo(2,1,j)**2
cxx      cfac=(1./geo(1,1,j)**2 - c1inv)/dmax1(1.d-20,h(j))
c: see p. 66 (from SAFARI manual for attenuation in k=w/c):
         c1inv=dcmplx(1.d0,geo(2,4,j)*geo(2,1,j)/fac)/geo(2,1,j)
         c1invsq=c1inv*c1inv
         c2inv=dcmplx(1.d0,geo(1,4,j)*geo(1,1,j)/fac)/geo(1,1,j)
         c2invsq=c2inv*c2inv
         cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
         do jsr=jj2,jj1,-1
            jptr=jptr + 1
            zmx(jptr)=zdep(j)-zsr(jsr)
            zmx_im_gbs(jptr)=-zsr_im_gbs(jsr)
c: jsrmx(jptr)=jsr, the index of depth zmx(jptr) in the array zsr(1:nzsr):
            jsrmx(jptr)=jsr
            rho_sr(jsr)=geo(2,3,j) + zmx(jptr)*rhfac
            cp_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            cs_sr(jsr)=(0.d0,0.d0)
c: Pointer from zmx to zm [find depth zsr(jsr) in zmx(mx_m(jsr))]:
            mx_m(jsr)=jptr
c: Pointer from zsr to geo [find s/r depth zsr(jsr) in geo(:,:,j)]:
            jsr2j(jsr)=j
         enddo
         if(iisol(j) .eq. 1 .and. jj2 .ge. jj1) then
            c1inv=dcmplx(1.d0,geo(2,5,j)*geo(2,2,j)/fac)/geo(2,2,j)
            c1invsq=c1inv*c1inv
            c2inv=dcmplx(1.d0,geo(1,5,j)*geo(1,2,j)/fac)/geo(1,2,j)
            c2invsq=c2inv*c2inv
            cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
            do jsr=jj2,jj1,-1
               kksh(jsr)=1
               cs_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            enddo
         endif
cxx   print *,'j,z1,z2,h,zmx = ',j,z1,z2,h(j),(zmx(jx1),jx1=
cxx  .   jzmx(j),jzmx(j)+nzmx(j)-1)
         z1=z2
      enddo
c
c: For layers below reference depth:
      do j=jflu(1,1)+jflu(1,5),nlay
         if(j .ne. nlay) then
            z2=zdep(j)
            if(j .ge. jsurf) z2=z2 + 1.d-10*abs(z2)
         else
            z2=1.d20
         endif
         do while(zsr(jj) .lt. z1 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj1=jj
         do while(zsr(jj) .lt. z2 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj2=jj-1
         jzmx(j)=jptr+1
         nzmx(j)=jj2-jj1+1
cxx      c1inv=1./geo(1,1,j)**2
cxx      cfac=(1./geo(2,1,j)**2 - c1inv)/dmax1(1.d-20,h(j))
         c1inv=dcmplx(1.d0,geo(1,4,j)*geo(1,1,j)/fac)/geo(1,1,j)
         c1invsq=c1inv*c1inv
         c2inv=dcmplx(1.d0,geo(2,4,j)*geo(2,1,j)/fac)/geo(2,1,j)
         c2invsq=c2inv*c2inv
         cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
         rhfac=(geo(2,3,j) - geo(1,3,j))/dmax1(1.d-20,h(j))
         do jsr=jj1,jj2
            jptr=jptr + 1
            zmx(jptr)=zsr(jsr)-zdep(j-1)
            zmx_im_gbs(jptr)=zsr_im_gbs(jsr)
c: jsrmx(jptr)=jsr, the index of depth zmx(jptr) in the array zsr(1:nzsr):
            jsrmx(jptr)=jsr
            rho_sr(jsr)=geo(1,3,j) + zmx(jptr)*rhfac
            cp_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            cs_sr(jsr)=(0.d0,0.d0)
c: Pointer from zmx to zm [find depth zsr(jsr) in zmx(mx_m(jsr))]:
            mx_m(jsr)=jptr
c: Pointer from zsr to geo [find s/r depth zsr(jsr) in geo(:,:,j)]:
            jsr2j(jsr)=j
         enddo
         if(iisol(j) .eq. 1 .and. jj1 .le. jj2) then
            c1inv=dcmplx(1.d0,geo(1,5,j)*geo(1,2,j)/fac)/geo(1,2,j)
            c1invsq=c1inv*c1inv
            c2inv=dcmplx(1.d0,geo(2,5,j)*geo(2,2,j)/fac)/geo(2,2,j)
            c2invsq=c2inv*c2inv
            cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
            do jsr=jj1,jj2
               kksh(jsr)=1
               cs_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            enddo
         endif
cxx   print *,'j,z1,z2,h,zmx = ',j,z1,z2,h(j),(zmx(jx1),jx1=
cxx  .   jzmx(j),jzmx(j)+nzmx(j)-1)
         z1=z2
      enddo
c
c: Total number of depths in zmx, which includes source, receiver, and
c: layer interface depths:
      nzmxtot=jptr
c
c: Initialize s-wave potentials to zero in water and fluid layers:
      do ii=1,2
         do j=jflu(ii,1),jflu(ii,2),jflu(ii,3)
            do jx1=jzmx(j),jzmx(j)+nzmx(j)-1
               psix(jx1)=(0.,0.)
               dpsix(jx1)=(0.,0.)
            enddo
         enddo
      enddo
c: iiww(j) is a flag that transmission coefficients need to be computed in
c: the layer j:
      if(nzsr .gt. 0) then
         jlmin=min(nsvmin,jsr2j(1))
         jlmax=max(nsvmin,jsr2j(nzsr))
      else
         jlmin=nsvmin
         jlmax=nsvmin
      endif
      do j=1,nlay
         iiww(j)=0
         if(j .ge. jlmin .and. j .le. jlmax) iiww(j)=1
c: Temp (see also duct_check):
         iiww(j)=1
      enddo
c: Find maximum receiver sound speed for use in mode_field:
      crmax=cfmin
      do jsr=1,nzsr
         crmax=dmax1(crmax,dreal(cp_sr(jsr)))
      enddo
c
c: DON'T DO THIS ANY MORE SINCE WE CAN HAVE MULTIPLE SOURCES:
c: Normalize densities to density at source:
cxx   rho_src=rho_sr(mzsrc(1))
cxx   do j=1,nzsr
cxx      rho_sr(j)=rho_sr(j)/rho_src
cxx   enddo
c
      return
      end
      subroutine xi_cut_set
c
c: Sets gradient in Airy halfspaces that pass through |R1R2|=1 contours
c: as far from modes as possible in order to minimize the contribution of
c: the branch line modes.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      complex*16 xksqx,etasqx
      integer*4 iibt,jlay,ii2
c
c: Make main contour cross cut with phase within pi/16 of +/-pi:
      jlay=nlay
      do iibt=1,2
         if(iibt .eq. 2) jlay=1
         ii2=3-iibt
         if(ii_xi_mv(iibt) .eq. 1 .and. 
     .      dreal(xk(iibt,jlay)) .lt. dreal(xkref)) then
            xksqx=xksq(iibt,jlay)
            etasqx=etasq(jlay)
            call sheet_init(xk(iibt,jlay),1,iish,iish_ref)
            call xi_cross_calc(xksqx,etasqx,h(jlay),eta(jlay),
     .         etasq(jlay),xk(ii2,jlay),xksq(ii2,jlay),pie,twpie,kw0)
c
            if(iisol(jlay) .eq. 1 .and. iss(jlay) .eq. 0 .and. 
     .         dreal(xb(iibt,jlay)) .lt. dreal(xkref)) then
               xksqx=xbsq(iibt,jlay)
               etasqx=etbsq(jlay)
               call sheet_init(xk(iibt,jlay),1,iish,iish_ref)
               call xi_cross_calc(xksqx,etasqx,h(jlay),etb(jlay),
     .            etbsq(jlay),xb(ii2,jlay),xbsq(ii2,jlay),pie,twpie,kw0)
            endif
         endif
      enddo
c
      return
      end
ccc
      subroutine xi_cross_calc(xksqx,etasqx,hx,eta_new,etasq_new,
     .   xk_new,xksq_new,pie,twpie,kw0)
c
c: Finds point in complex k plane where negative xi axis is crossed
c: by |R1R2|=1 contour.
c
      implicit none
      complex*16 xksqx,etasqx,eta_new,etasq_new,xk_new,xksq_new,
     .   xi_st,k1,ln1,dln1,k2,ln2,dln2,xi1,xi2,kpie,r1r2(3,4),d_ksq,
     .   eimpi3,eimpi6,kt
      real*8 hx,pie,twpie,kw0,ln_im1,ln_im2,ln_lim,dph_left,
     .   etasq_mag,xi_im,meta_dir,ln_re_eps,ln_mid,ln_dif
      integer*4 iilhp,ii2nd
      data eimpi3/(0.5d0,-0.866025403784439d0)/
      data eimpi6/(0.866025403784439d0,-0.5d0)/
c
c: Make arg(R1R2) be within pi/8 of pi (ln_lim=7*pi/16):
      ln_lim=2.74889
      xi_im=0.1d0
      xi_im=0.5d0
      xi_im=1.0d0
      ln_re_eps=0.1d0
cc    if(ii2nd .eq. 0) then
cc       print *,'Enter xi_im,ln_re_eps: ',xi_im,ln_re_eps
cc       read(5,*) xi_im,ln_re_eps
cc    endif
      ii2nd=1
c
c: Find point just to right of Im(xi)=0 cut for which |R1R2|=1:
      xi_st=dcmplx(-4.d0,xi_im)
      call xi_cross(xksqx,etasqx,xi_st,k1,ln1,dln1,xi1,ln_im1,
     .   ln_re_eps)
c
c: If contour crosses xi axis too close to origin, do not try to move cut:
      if(dreal(xi1) .gt. -0.5d0) then
cc       print *,'xi1 too close to origin: ',xi1
         return
      endif
c
c: Find point just to left of Im(xi)=0 cut for which |R1R2|=1:
      xi_st=dcmplx(dreal(xi1),-xi_im)
      call xi_cross(xksqx,etasqx,xi_st,k2,ln2,dln2,xi2,ln_im2,
     .   ln_re_eps)
c
cc    print *,'ln_im1,ln_im2 = ',ln_im1*180./pie,ln_im2*180./pie
c
      ln_dif=ln_im2 - ln_im1
      if(abs(ln_dif) .gt. pie) ln_dif=ln_dif-sign(twpie,ln_dif)
      ln_mid=ln_im1 + 0.5d0*ln_dif
      if(abs(ln_mid) .gt. pie) ln_mid=ln_mid-sign(twpie,ln_mid)
c
      if(abs(ln_mid) .gt. ln_lim) then
cc       print *,'Did not have to move: ',ln_im1*180./pie,
cc   .      ln_im2*180./pie
         return
      endif
c
      if(ln_mid .le. 0.d0) then
         dph_left=-pie - ln_im1
         kt=k1
         call traj_mag(0.d0,dph_left,k1,ln1,dln1,kpie,r1r2,
     .      0.1d0,0.19635d0,iilhp)
cc       print *,'Moved to right: k,arg(R1R2) = ',dreal(kt)/kw0,
cc   .      dreal(kpie)/kw0,ln_im1*180/pie,dimag(ln1)*180./pie
      else
         dph_left=pie - ln_im2
         kt=k2
         call traj_mag(0.d0,dph_left,k2,ln2,dln2,kpie,r1r2,
     .      0.1d0,0.19635d0,iilhp)
cc       print *,'Moved to left: k,arg(R1R2) = ',dreal(kt)/kw0,
cc   .      dreal(kpie)/kw0,ln_im2*180/pie,dimag(ln2)*180./pie
      endif
c
c: Now adjust phase of etasq so that cut passes through kpie:
      etasq_mag=cdabs(etasqx)
      d_ksq=kpie*kpie - xksqx
c: Make new gradient have same magnitude as old, but with same phase
c: as kpie^2 - K^2 so that Im(xi)=0 cut passes through kpie:
      etasq_new=dsign(etasq_mag,dreal(xi1))*d_ksq/cdabs(d_ksq)
      meta_dir=datan2(-dimag(eta_new),-dreal(eta_new))*180.d0/pie
cc    print *,'old eta: ',meta_dir
c: Place eta_new somewhere in left half plane so that arg(-eta) has chance
c: at being between -pi/6 and -pi/3:
      eta_new=-cdsqrt(etasq_new)
      meta_dir=datan2(-dimag(eta_new),-dreal(eta_new))*180.d0/pie
c: Make sure arg(-eta) between -60 and 0 degrees:
      if(meta_dir .lt. -60.d0 .or. meta_dir .gt. -30.d0) then
cc       print *,'eta placed outside of valid range: ',meta_dir
         if(meta_dir .lt. -60.d0) then
            eta_new=-dsqrt(etasq_mag)*eimpi3
         else
            eta_new=-dsqrt(etasq_mag)*eimpi6
         endif
         etasq_new=eta_new*eta_new
      endif
c
      meta_dir=datan2(-dimag(eta_new),-dreal(eta_new))*180.d0/pie
cc    print *,'new eta: ',meta_dir
      if(abs(meta_dir-45.d0) .lt. 5.d0) then
         print *,'WARNING: mode_traj may have paralleled cut '
      endif
c
      xksq_new=xksqx + hx*eta_new*etasq_new
      xk_new=cdsqrt(xksq_new)
c
         xi_st=dcmplx(-2.d0,xi_im)
         call xi_cross(xksq_new,etasq_new,xi_st,k1,ln1,dln1,xi1,ln_im1,
     .      ln_re_eps)
cc       print *,'New phase to right of xi cut: ',ln_im1*180./pie
         xi_st=dcmplx(-2.d0,xi_im)
         call xi_cross(xksq_new,etasq_new,xi_st,k2,ln2,dln2,xi2,ln_im2,
     .      ln_re_eps)
cc       print *,'New phase to left of xi cut: ',ln_im2*180./pie
c
      return
      end
ccc
      subroutine xi_cross(xksqx,etasqx,xi_st,k0,ln0,dln0,xi0,ln_im0,
     .   ln_re_eps)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      real*8 ln_im0,ln_re_eps,ln_re1,ln_re2,ln_re0,xi_im
      complex*16 xi_st,k0,ln0,dln0,xi0,xi1,xi2,xksqx,
     .   etasqx
      integer*4 ntry
c
      xi_im=dimag(xi_st)
      xi1=xi_st
      call i_xi_calc2(xi1,xksqx,etasqx,1,k0,ln0,dln0,ln_re1,ln_im0,
     .     jjfail)
      if(jjfail .gt. 0) return
c
      ntry=0
      if(ln_re1 .lt. 0.d0) then
10       continue
         xi2=dcmplx(2.d0*dreal(xi1),xi_im)
         ntry=ntry + 1
         if(ntry .gt. 20) print *,'ntry trouble in exit'
         call i_xi_calc2(xi2,xksqx,etasqx,1,k0,ln0,dln0,ln_re2,ln_im0,
     .        jjfail)
         if(jjfail .gt. 0) return
         if(ln_re2 .lt. 0.d0) then
            xi1=xi2
            goto 10
         endif
      else
         xi2=xi1
         ln_re2=ln_re1
15       continue
         if(dreal(xi2) .gt. -2.0d0) then
c: Allow to go across zero if necessary:
            xi1=dcmplx(dreal(xi2)+2.d0,xi_im)
         else
            xi1=dcmplx(0.5d0*dreal(xi2),xi_im)
         endif
         ntry=ntry + 1
         if(ntry .gt. 20) print *,'ntry trouble in exit'
         call i_xi_calc2(xi1,xksqx,etasqx,1,k0,ln0,dln0,ln_re1,ln_im0,
     .        jjfail)
         if(jjfail .gt. 0) return
         if(ln_re1 .gt. 0.d0) then
            if(dreal(xi1) .gt. -0.5d0) then
               xi0=xi1
               return
            endif
            xi2=xi1
            goto 15
         endif
      endif
c
c: We now have |R1R2|=1 bracketed on xi axis:
20    xi0=0.5d0*(xi1 + xi2)
      call i_xi_calc2(xi0,xksqx,etasqx,2,k0,ln0,dln0,ln_re0,ln_im0,
     .     jjfail)
      if(jjfail .gt. 0) return
      ntry=ntry + 1
      if(ntry .gt. 25) print *,'ntry trouble in exit'
      if(dabs(ln_re0) .gt. ln_re_eps) then
         if(ln_re0 .lt. 0.d0) then
            xi1=xi0
         else
            xi2=xi0
         endif
         goto 20
      endif
c
      return
      end
ccc
      subroutine i_xi_calc2(xi0,xksqx,etasqx,ndv,k0,ln0,dln0,
     .   ln_re0,ln_im0,jjfail)
c
      implicit none
      complex*16 xi0,xksqx,etasqx,k0,ln0,dln0,rr0(3,4)
      real*8 ln_re0,ln_im0
      integer*4 ndv,jjfail
c
      k0=cdsqrt(xksqx + xi0*etasqx)
      call r1r2_calc(k0,rr0,ndv,1,jjfail)
      if(jjfail .gt. 0) return
      ln0=rr0(1,4)
      dln0=rr0(2,4)
      ln_re0=dreal(ln0)
      ln_im0=dimag(ln0)
c
      return
      end
      subroutine bb_field_reuse(knx,phiz,expz_gbs,jmk,tfz,jfbb)
c
c: Computes field tf for a single mode characterized by eigenvalue k and mode 
c: functions phi and dpsi.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jmk,jfbb,jmo,jr,jsr,kk,kk0,jsrc,jrec,j,jzs,msrc,nrmax,
     .   jd,nsv,isv,jzsr,msrcmax,nsrcmax
      complex*8 phiz(nzsrgeom,jmk),tfz(nfbb,nrec,nrng),dphi_z
      complex*16 knx,sqkn,cfac2,h_arg,H0,iknx,phi_src
      real*8 expz_gbs(nzsrgeom,jmk),xn_beam,rsig_max,magsq
      real*8 phi_src_re,phi_src_im,dphi_src_re,dphi_src_im,
     .   dzsrgeom
      real*4 phiw_max,magsq_c8,ampl
c
c: Find mode function at source:
c: Normalize by density at source (9-1-94):
      jzs=1
      msrc=mzsrc(jzs)
      msrcmax=0
c: pln
c: hunt finds the index before zsrgeom becomes greater than zsr
      call hunt(zsrgeom,nzsrgeom,zsr(msrc),msrcmax)

      dzsrgeom = (zsr(msrc)-zsrgeom(msrcmax))
     .     /(zsrgeom(msrcmax+1)-zsrgeom(msrcmax))

      dphi_src_re=(real(phiz(msrcmax+1,jmk))
     .    -real(phiz(max(msrcmax,1),jmk)))*dzsrgeom
      dphi_src_im=(imag(phiz(msrcmax+1,jmk))
     .    -imag(phiz(max(msrcmax,1),jmk)))*dzsrgeom

      phi_src=dcmplx(dphi_src_re+real(phiz(msrcmax,jmk)),
     .     dphi_src_im+imag(phiz(msrcmax,jmk)))
     .     /rho_sr(msrcmax)

cpln      write(6,*)zsrgeom(msrcmax),zsrgeom(msrcmax+1)
cpln      write(6,*)zsr(msrc),phi_src

      xn_beam=(w/dreal(cp_sr(msrcmax)))*b_gbs(jzs) - 
     .     expz_gbs(msrcmax,jmk)
      iknx=dcmplx(-dimag(knx),dreal(knx))
      sqkn=cdsqrt(knx)
c
      do jrec=1,nrec
         jsr=mzrec(jrec)
         nsrcmax=0
         call hunt(zsrgeom,nzsrgeom,zsr(jsr),nsrcmax)
         dzsrgeom = (zsr(jsr)-zsrgeom(nsrcmax))
     .        /(zsrgeom(nsrcmax+1)-zsrgeom(nsrcmax))
         dphi_src_re=real(phiz(nsrcmax+1,jmk)
     .        -phiz(max(nsrcmax,1),jmk))*dzsrgeom
         dphi_src_im=imag(phiz(nsrcmax+1,jmk)
     .        -phiz(max(nsrcmax,1),jmk))*dzsrgeom
         phisr(jsr)=dcmplx(dphi_src_re+real(phiz(nsrcmax,jmk)),
     .        dphi_src_im+imag(phiz(nsrcmax,jmk)))
cpln         write(6,*)zsrgeom(nsrcmax+1),zsrgeom(nsrcmax)
cpln         write(6,*)jmk,zsr(jsr),phisr(jsr)
cpln         phisr(jsr)=phiz(jsr,jmk)
      enddo
c
c: Find range beyond which mode is insignificant (ranges have been sorted):
      rsig_max=kim_fac/dmax1(1.d-20,dimag(knx)-kim_bb(jfbb))
      nrmax=0
      call hunt(range,nrng,rsig_max,nrmax)
c
      do jr=1,nrmax
         kk0=krec_jr(jr)
         if(tilth) then
            kk=kk0+jr
            jsrc=jrec_jr(1,kk)
            jrec=jrec_jr(2,kk)
            h_arg=knx*range(jr)
            if(magsq(h_arg) .gt. 25.d0) then
czs: Include normalization by exp(-xn_beam) here:
               cfac2=sq2pir(jr,jrec)*phi_src*
     .              cdexp(iknx*range(jr) - xn_beam)/sqkn
            else
               call cdhankel(h_arg,1.d-6,H0)
               cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
            endif
            jsr=mzrec(jrec)
            tfz(jfbb,jrec,jsrc)=tfz(jfbb,jrec,jsrc) + 
     .           cfac2*phisr(jsr)
         else
            do kk=kk0+1,kk0+nrec_jr(jr)
               jsrc=jrec_jr(1,kk)
               jrec=jrec_jr(2,kk)
               h_arg=knx*(range(jr)+dtiltvp(jrec))
               if(magsq(h_arg) .gt. 25.d0) then
czs: Include normalization by exp(-xn_beam) here:
                  cfac2=sq2pir(jr,jrec)*phi_src*
     .                 cdexp(iknx*(range(jr)+dtiltvp(jrec))
     .                  - xn_beam)/sqkn
               else
                  call cdhankel(h_arg,1.d-6,H0)
                  cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
               endif
               jsr=mzrec(jrec)
               tfz(jfbb,jrec,jsrc)=tfz(jfbb,jrec,jsrc) + 
     .         cfac2*phisr(jsr)
            enddo
         end if
      enddo
c
      return
      end
      subroutine clean_array
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      include 'lab_com'
      include 'debug_com'

      integer j,jj,izero
      real*4 zero
      real*8 dzero
      complex*16 dczero
      data dczero/(0.d0,0.d0)/
      data dzero/0.d0/
      data zero/0.0/
      data izero/0/

cpln      nzs=0
cpln      nrec=0
cpln      nzmf=0
cpln      nsrc=0
      nzsr=0
      nlay=izero
      nlayb=izero
      nlayt=izero
      nduct=izero
      nzmxtot=izero
      npt=0
      nmin_dln=izero
      kduct=izero
      kduct0=izero
      isvmin=izero
      nsvmin=izero
      jsurf=izero
      jobot=izero
      ii_xi_mv(1)=izero
      ii_xi_mv(2)=izero
c This has caused me a lot of problems
c Has to be reset for multible calls of bb_brute
      cphlo=dzero
c
c Variables used in traj_sdp
c
      k_sdp=dczero
      ln_sdp=dczero
      dln_sdp=dczero
      k_spt=dczero
      ln_spt=dczero
      dln_spt=dczero
c
      do j=1,NLMAX
         zduct(j)=dzero
         isp(j)=izero
         iss(j)=izero
         iiww(j)=izero
         zdep(j)=dzero
cpln         nzmx(j)=izero
cpln         jzmx(j)=izero
         do jj=1,2
            aisoln(jj,j)=izero
         end do
      end do

      do j=1,NSRMAX
         jsrmx(j)=izero
         jsr2j(j)=izero
         mx_m(j)=izero
cpln         zmx(j)=dzero
cpln         zrec(j)=dzero
cpln         zsrc(j)=dzero
         zsr(j)=dzero
cpln         phix(j)=dczero
cpln         dphix(j)=dczero
      end do
      
      do j=1,5
         do jj=1,2
            jflu(jj,j)=izero
         end do
      end do
      
      do j=1,NDMAX
         mzduct(j)=izero
         zpeak(j)=dzero
         do jj=1,5
            jduct(jj,j)=izero
         end do
      end do
      zpeak(NDMAX+1)=dzero

c      do j=1,32
c         jmin_dln(j)=izero
c         dln_deep(j)=zero
c      end do
c
c      do j=1,150
c         k_cont(j)=dczero
c         ln_cont(j)=dczero
c         dln_cont(j)=dczero
c      end do
c
      return
      end

