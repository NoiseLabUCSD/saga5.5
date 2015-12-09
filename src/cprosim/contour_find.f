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
