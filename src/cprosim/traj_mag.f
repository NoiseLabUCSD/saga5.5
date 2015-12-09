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
