      subroutine eig_findm(k0,rr0,r1r2,nm1)
c
c: Finds eigenvalues by moving along constant ln|R1*R2|=1 in complex 
c: k-plane.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 iibad,nm1
      complex*16 k0,rr0(3,4),r1r2(3,4),dW_dk
      real*8 xi_real,xi_rat
c
      nmode=nmode + 1
      ncalc(nmode)=0
c
      iibad=0
      call mem_lim(nmode,nm_max2,MLINE,LML,'nmode',5,'nm_max2',7,
     .   iibad,0)
      if(iibad .ne. 0) then
         iifail=1
         nmode=nmode - 1
         print *,'Mode search stopped because of nmode*nzsr limit.'  
         write(2,*) 'Mode search stopped because of nmode*nzsr limit.'  
         return
      elseif(nmode .gt. nm_lim) then
         if(iiwrite .eq. 1) then
            print *,'Mode search stopped due to mode number.'
            write(2,*) 'Mode search stopped due to mode number.'
         endif
         iidone=1
         return
      endif
c
      call eig_find0(k0,rr0,kn(nmode),r1r2,nm1)
      k0=kn(nmode)
      rr0(1,4)=r1r2(1,4)
      rr0(2,4)=r1r2(2,4)
c
      if(iidone .eq. 1 .or. iifail .gt. 0) return
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
cc      print *,'Calling mode_fun from eig_findm: '
      call mode_fun(k0,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,psi,dpsi,
     .   exp_gbs,nmode)
cc      print *,'done'
      if(iidiag .eq. -2) then
         call phz_calc
      endif
c
c: Compute mode amplitudes at duct reference depths:
cxx   if(nblm .gt. nblm_max .and. isp(nlay) .eq. 0) then
      if(isp(nlay) .eq. 0) then
         xi_real=dreal(xi_hsp(1,1))
         xi_rat=dabs(xi_real/dimag(xi_hsp(1,1)))
         if(xi_rat .gt. 20. .and. xi_real .lt. 1.) nblm=nblm + 1
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
      use parms_com
      use i_o_com
      use gen_com
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
      if(iidone .eq. 1 .or. iifail .gt. 0) return
c
      if(iilhp .eq. 1) then
         call fix_path(k0,rr0,1,1)
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
      call eig_final(k,r1r2,errdkms,iimt,kw0,iish,iifail)
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
c: OK if we found a branch line mode very close to branch point.
         else
            nc_ref=nc_ref + 1
            if(nc_ref .gt. 1) then
               if(iidiag .ge. 1) print *,'Informative message: ',
     .            'XKREF island mode not resolved.'
               iifail=2
               return
            endif
            if(iidiag .ge. 1) print *,'FALSE MODE AT XKREF IGNORED: ',
     .         dreal(xkref)/kw0,dimag(xkref)*8685.9
            kcut=xkref + dcmplx(0.d0,.1d0/8686.d0)
            call fix_path(k0,rr0,1,0)
            phcut=dimag(rr0(1,4))
            if(iidone .eq. 1 .or. iifail .eq. 1) return
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
                     print *,'traj_sdp unsuccessful '
                     print *,'iicirc>4: Tell EKW ',nmode
                     nmode=nmode - 1
                     return
                  endif
               endif
               iimst=0
               goto 999
            else
               print *,'iicirc>4: Tell EKW ',nmode
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
            call fix_path(k0,rr0,1,0)
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
            if(iidiag .ne. 0 .and. iiwrite .eq. 1) then
               print *,'Mode search stopped due to high mode '//
     .            'attenuations'
               write(2,*) 'Mode search stopped due to high mode '//
     .            'attenuations'
            endif
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
      use parms_com
      use i_o_com
      use gen_com
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
      subroutine eig_final(k,r1r2,errdkms,iimt,kw0,iish,iifail)
c
c: Finds the mode eigenvalue once we are in neighborhood of it.
c
      implicit none
      integer*4 iimt,iish(2,2),ntry,ntry2,iifail
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
      call r1r2_calc(k,r1r2,3,1)
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
         call r1r2_calc(k,r1r2,3,1)
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
      use parms_com
      use i_o_com
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
      use parms_com
      use i_o_com
      use gen_com
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
            call fix_path(k0,rr0,1,0)
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
