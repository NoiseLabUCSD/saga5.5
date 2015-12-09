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
