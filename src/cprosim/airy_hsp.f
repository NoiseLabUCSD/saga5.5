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
