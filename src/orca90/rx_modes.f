      subroutine rx_modes(iiwrt)
co	Modified for use as subroutine
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 iiwrt
      real*8 kmax,kmin
      real*4 cpu1,cpu2
c
co      call cpu_time(cpu1)
c
      call rx_prep
      call rx_zmx_init
      call rx_freq_init(kmax,kmin)
c
      nmode=0
      nctot=0
      nclast=0
      nzoff=0
      kim_min=1.d100
c
      call rx_mode_int(kmax,kmin,nm_lim,3,0,1,0,iiwrt)
c
      if(nmode .ge. nm_max2) then
         print *,'MAX # MODES REACHED. INCREASE nm_max2 CALC',
     .      ' IN sr_geom ',nm_max2
         print *,'Solution will include only ',nmode,' modes ...'
         write(2,*) 'MAX # MODES REACHED. INCREASE nm_max2 CALC',
     .      ' IN sr_geom ',nm_max2
         write(2,*) 'Solution will include only ',nmode,' modes ...'
      endif
c
      if(nzoff .gt. 0) then      
co         print *,'Info msg: # of zeros of mode function off by one ',
co     .      nzoff,' times.  Can happen with multiple duct environments.'
co         write(2,*) 'Info msg: # of zeros of mode function off by one ',
co     .      nzoff,' times.  Can happen with multiple duct environments.'
      endif
c
      if(nm_miss .gt. 0) then
         print *,'WARNING: ORCA missed ',nm_miss,' modes. Tell EKW.'
         write(2,*) 'WARNING: ORCA missed ',nm_miss,' modes. Tell EKW.'
      endif
      if(iiwrt .eq. 1) then
co         call cpu_time(cpu2)
co         print *,'nmode,nctot = ',nmode,nctot,float(nctot)/nmode
co         print *,'Time taken to find modes and mode functions = ',
co     .      cpu2-cpu1
      endif
c
      return
      end
ccc
      subroutine rx_mode_int(kmax,kmin,nmax,ndv,iim0,iiwk_check,iibb,
     .   iiwrt)
c
c: Finds modes on real axis from kmax down to kmin.
c: iim0=1 means we are starting at or near an already found mode, 
c: so do not duplicate it.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 nmax,ndv,iim0,iiwk_check,iibb,iiwrt,iipol5,iiek,iiekp,
     .   iiekpp,nstack,iiroots,iival,iimok,ndv0,ndvr,nloopm,ii25,
     .   mode
      real*8 k1,k2,ek1(6),ek2(6),phm1,phm2,kmin,kmax,dphm,kx,
     .   stack(6,100),kext,ekext(6),phmext,ekx(6),phmx,
     .   kmid,ekmid(6),phmid,phmn1,ek_rat
      real*8 kn_imag

c
      nm_miss=0
cc    print *,'enter iipol5 (0=cubic,1=poly5)'
cc    read(5,*),iipol5
c: TEMP check of e'':
cpp   print *,'enter factor (1.d-10): '
cpp   read(5,*),kx
cpp   kx=(kmax-kmin)*kx
cpp   do iifnd=1,50
cpp      k1=kmax - iifnd*(kmax-kmin)/50
cpp      call rx_ek_calc(k1,ek1,phm1,6)
cpp      k2=k1 + kx
cpp      call rx_ek_calc(k2,ek2,phm2,6)
cpp      kmid=(ek2(1)-ek1(1))/kx
cpp      phmid=(ek2(2)-ek1(2))/kx
cpp      print *,'j: ',iifnd,sngl((kmid-ek1(2))/(ek2(2)-ek1(2))),
cpp  .      sngl((phmid-ek1(5))/(ek2(5)-ek1(5))),sngl(w/k1)
cpp   enddo
c
      iipol5=0
      ndv0=2
      ndvr=ndv
      if(kmax .le. kmin) return
      ncalc(nmode+1)=0
c
      k1=kmax
      call rx_ek_calc(k1,ek1,phm1,ndv0)
c
c: EKW FIX: 4-10-97 Do this before possibly changing k1 and ek1:
c: Copy mode quantities at kmax to kn(0),ekn(:,0), etc:
	call ek_copy(k1,kn(mode),ek1,ekn(1,nmode),phm1,phm1)

      mode_phz(1,nmode)=phm1
c
c: Check if starting from an already-found mode:
      if(iim0 .eq. 1) then
c: EKW FIX: 4-10-97:
c: Check if current mode might be found twice:
         ek_rat=ek1(1)/ek1(2)
         if(ek_rat .gt. 0.d0) then
c: If so, move left across mode using first derivative fit:
            ek1(1)=-ek1(1)
            k1=k1 - 2.d0*ek_rat
         endif
      endif
c
      k2=kmin
      call rx_ek_calc(k2,ek2,phm2,ndv0)
c
      nstack=0
c: iiek,iiekp=-1 when ek,ekp have changed sign:
25    continue
      iiekp=sign(1.d0,ek1(2)*ek2(2))
      if(iiekp .eq. 1) then
         iiekpp=sign(1.d0,ek1(2)*(ek1(1)-ek2(1)))
c: Check for case where e derivatives same sign, but change in e opposite:
         if(iiekpp .eq. -1) then
            call ek_stack(nstack,k2,ek2,phm2,stack,1,0)
            k2=.5d0*(k1 + k2)
            call rx_ek_calc(k2,ek2,phm2,ndv0)
            goto 25
         endif
      endif
      iiek=sign(1.d0,ek1(1)*ek2(1))
      dphm=phm2 - phm1
      if(dphm .gt. 1.5d0) then
         call ek_stack(nstack,k2,ek2,phm2,stack,1,0)
         k2=.5d0*(k1 + k2)
         call rx_ek_calc(k2,ek2,phm2,ndv0)
         goto 25
      elseif(iiek .eq. -1) then
c: e(k) has changed sign.  Find root between:
         call ek_root(k1,k2,ek1,ek2,kx,ekx,phmx,dk_max,ndvr,iipol5)
         nmode=nmode + 1
c
         call rx_nz_check(kx,ekx,phmx,ndvr,ndv0,
     .      k1,ek1,phm1,k2,ek2,phm2,nstack,stack,ii25,iibb)
         if(ii25 .eq. 1) goto 25
c
c: Check for weak modes at rmin:
         if(iiwk_check .eq. 1) then
            kn_imag=dimag(kn(nmode))
            kim_min=dmin1(kn_imag,kim_min)
            kim_max=kim_min + dkim
            if(kn_imag .gt. kim_max) then
               nhigh=nhigh + 1
c: Modes in different ducts may have different attenuations (ndrx = #ducts):
               if(nhigh .ge. ndrx+2) then
                  nmode=nmode - 1
                  if(iidiag .ne. 0 .and. iiwrt .eq. 1) then
                     print *,'Informative message: Mode search ',
     .                  'terminating due to large Im(k) (see rmin).'
                     write(2,'(a)') 'Informative message: Mode search '
     .                  //'terminating due to large Im(k) (see rmin).'
                  endif
                  return
               endif
            else
               nhigh=0
            endif
         endif
c
cc    print *,'nmode,phmx = ',nmode,dreal(kx),sngl(phmx),ncalc(nmode),
cc   .   nzref(nmode),sngl(eig_char(4,nmode))
c: Return if number of modes reaches nmax:
         if(nmode .ge. nmax) return
      elseif(iiekp .eq. -1 .and. ek1(2)*ek2(1) .gt. 0.) then
c: Derivative of e(k) changed sign and extremum is toward zero:
c: Find extremum:
         call ek_extreme(k1,k2,ek1,ek2,kext,ekext,phmext,dk_max,
     .      iiroots,ndv0)
         if(iidiag .ge. 2 .and. dphm .le. .5d0) 
     .      print *,'ep changed sign: ',dphm,iiroots
         if(iiroots .eq. 1) then
c: If extremum on other side of zero, push k2 on stack, set k2 to 
c: extremum, and go back:
            call ek_stack(nstack,k2,ek2,phm2,stack,1,0)
            call ek_copy(kext,k2,ekext,ek2,phmext,phm2)
            goto 25
cc       elseif(iidiag .ge. 2) then
cc          print *,'IIROOTS=0 from ek_extreme'
         endif
      endif
38    continue
      if(nstack .gt. 0) then
         call ek_copy(k2,k1,ek2,ek1,phm2,phm1)
40       call ek_stack(nstack,k2,ek2,phm2,stack,-1,iival)
         if(iival .gt. 0) then
c: If k value on stack was an eigenvalue (from ek_miss2), enter it and 
c: get next k value from stack:
            nmode=nmode + 1
c
c: Must make new call to rx_ek_calc since we've lost G info:
            kx=k2
            call rx_ek_calc(kx,ekx,phmx,ndvr)
c
c: Pull off k to left of found mode to serve as k2:
            call ek_stack(nstack,k2,ek2,phm2,stack,-1,iival)
c
            call rx_nz_check(kx,ekx,phmx,ndvr,ndv0,
     .         k1,ek1,phm1,k2,ek2,phm2,nstack,stack,ii25,iibb)
            if(ii25 .eq. 1) then
               if(iidiag .ge. 2) print *,
     .            'Info msg: ii25=1 after iival=1'
               goto 25
            endif
c
c: Return if number of modes reaches nmax:
            if(nmode .ge. nmax) return
c
cc    print *,'nmode,phmx = ',nmode,dreal(kx),sngl(phmx),ncalc(nmode),
cc   .   nzref(nmode),sngl(eig_char(4,nmode))
c: Set value of ek just slightly on correct side so that mode not found twice:
cc??        ek2(1)=-sign(1.d-10,ek2(2))
            goto 38
         endif
         goto 25
      endif
c: Before exiting, check for large phase change from last mode found
c: to kmin in case two modes have been missed at end:
      phmn1=mode_phz(1,nmode)
      dphm=phm2 - phmn1
      nloopm=0
      if(dphm .gt. 2.d0) then
c: Certain that two modes are missing:
         nloopm=10
      elseif(dphm .gt. 1.25) then
c: Not certain that two modes are missing:
         nloopm=3
      endif
      if(nloopm .gt. 0) then
         call ek_copy(k2,kx,ek2,ekx,phm2,phmx)
         call ek_miss2_safe(kn(nmode),kx,ekn(1,nmode),ekx,
     .      phmn1,phmx,k1,k2,ek1,ek2,phm1,phm2,
     .      kmid,ekmid,phmid,dk_max,ndv0,nloopm,iimok,iidiag)

         if(iidiag .ge. 2) print *,'ek_miss2 at end: ',iimok,
     .      kn(nmode),kx,ekn(1,nmode),ekx(1),ekn(2,nmode),ekx(2),
     .      phmn1,phmx,kmid,ekmid(1),ekmid(2),phmid
         if(iimok .eq. 0) then
            if(nloopm .eq. 10) then
c: Two modes not missed after all (this should never happen):
               print *,'Two modes missed at end? nmode,phmn1,phmx,f= ',
     .            nmode,phmn1,phmx,f_hz
            endif
         else
cc    print *,'Two modes missed ',nmode,phmx
c: Push k value to left of missed extremum onto stack:
            call ek_stack(nstack,k2,ek2,phm2,stack,1,0)
            call ek_copy(kmid,k2,ekmid,ek2,phmid,phm2)
c: If k1 is at a mode, adjust e value so that it is not found twice:
            if(k1 .eq. dreal(kn(nmode))) 
     .         ek1(1)=-sign(1.d-10,ek1(2))
            goto 25
         endif
      endif
c
      return
      end
ccc
      subroutine rx_nz_check(kx,ekx,phmx,ndvr,ndv0,
     .   k1,ek1,phm1,k2,ek2,phm2,nstack,stack,ii25,iibb)
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 iibb,jlay_ref,nstack,
     .   iival,iimok,ndv0,ndvr,nzero,ntry,nz_err,j,ii25,ii_ref
      real*8 k1,k2,ek1(6),ek2(6),phm1,phm2,stack(6,100),
     .   kx,ekx(6),phmx,phmx2,kmid,ekmid(6),phmid,phmn1,vg,D2k_Df
c
      ii25=0
      ntry=1
12    call rx_enter(nmode,ncalc,nclast,nctot,kn,kx,ekx,ekn,phmx,
     .   mode_phz,w,iidiag)
      call rx_norm(jlay_ref,ii_ref,ndvr,vg,D2k_Df)
      if(iibb .eq. 0) then
         call rx_mf_lay(jlay_ref,ii_ref,vg,nzero)
      else
         call rx_mf_lay_Dw(jlay_ref,ii_ref,vg,nzero)
      endif
c: Account for possibly missing nm_miss modes:
      nz_err=iabs(nmode+nm_miss-nzero)
      if(nz_err .eq. 1) then
c: nzero is off by one, indicating mode may need to be found more precisely:
         if(iidiag .ge. 1) print *,'Info msg: nzero off by 1: ',
     .      nmode,nzero,ntry,f_hz
         if(ntry .le. 2) then
c: Polish root two times more:
            do j=1,2
               kx=kx - ekx(1)/ekx(2)
               call rx_ek_calc(kx,ekx,phmx,ndvr)
            enddo
            ntry=ntry + 1
            goto 12
         else
            nzoff=nzoff + 1
            if(iidiag .ge. 1) print *,'Warning msg: nzero off ',
     .         'by 1 (tell EKW): ',nmode,nzero,ntry,f_hz
         endif
      endif
c: Check if two (or more) modes missed between current and previous mode:
c: (Nearly degenerate modes can cause nzero to be off by one, which 
c: does not mean we have missed a mode.)
      if(nz_err .gt. 1) then
c: Push leftmost k value onto stack:
         call ek_stack(nstack,k2,ek2,phm2,stack,1,0)
c: Push current mode onto stack:
         call ek_stack(nstack,kx,ekx,phmx,stack,1,1)
c: Find value of k between modes (kmid) that has opposite sign:
         nmode=nmode - 1
         phmn1=mode_phz(1,nmode)
         if(iidiag .ge. 2) print *,'CALLING ek_miss2_safe: ',
     .      k2,kx,k1,kn(nmode),f_hz
         call ek_miss2_safe(kn(nmode),kx,ekn(1,nmode),ekx,
     .      phmn1,phmx2,k1,k2,ek1,ek2,phm1,phm2,
     .      kmid,ekmid,phmid,dk_max,ndv0,10,iimok,iidiag)

         if(iimok .eq. 0) then
c: Two modes not missed after all (this should never happen):
            print *,'Two modes missed: nmode,nzero,phase,f = ',
     .         nmode,nzero,phmx,f_hz
            nm_miss=nm_miss + 2
            call ek_stack(nstack,kx,ekx,phmx,stack,-1,iival)
            call ek_stack(nstack,k2,ek2,phm2,stack,-1,iival)
c: Must make new call to rx_ek_calc since we've lost G info:
            nmode=nmode + 1
            call rx_ek_calc(kx,ekx,phmx,ndvr)
            call rx_norm(jlay_ref,ii_ref,ndvr,vg,D2k_Df)
            if(iibb .eq. 0) then
               call rx_mf_lay(jlay_ref,ii_ref,vg,nzero)
            else
               call rx_mf_lay_Dw(jlay_ref,ii_ref,vg,nzero)
            endif
         else
cc    print *,'Two modes missed ',nmode,phmx
c: If left value of k is the mode, set sign of e to prevent double finding:
            if(k2 .eq. kx) ek2(1)=sign(1.d-10,ek2(2))
c: Push k value to left of missed extremum onto stack:
            call ek_stack(nstack,k2,ek2,phm2,stack,1,0)
            call ek_copy(kmid,k2,ekmid,ek2,phmid,phm2)
c: If k1 is at a mode, adjust e value so that it is not found twice:
            if(k1 .eq. dreal(kn(nmode))) 
     .         ek1(1)=-sign(1.d-10,ek1(2))
            ii25=1
            return
         endif
      endif
      if(iibb .eq. 0) then
         call rx_mode_fun(nmode,phi,dphi,jlay_ref,ii_ref,vg)
      else
         call rx_mode_fun_Dw(nmode,phi,dphi,Dphiz_w,Ddphiz_w,
     .      jlay_ref,ii_ref,vg)
      endif
cc    print *,'nmode,phmx = ',nmode,dreal(kx),sngl(phmx),ncalc(nmode),
cc   .   nzref(nmode),sngl(eig_char(4,nmode))
c
      return
      end
ccc
      subroutine rx_enter(nmode,ncalc,nclast,nctot,kn,kx,ekx,ekn,
     .   phmx,mode_phz,w,iidiag)
c
c: Enter mode into arrays.
c
      implicit none
      integer*4 nmode,ncalc(nmode+1),nclast,nctot,iidiag
      complex*16 kn(0:nmode)
      real*8 kx,ekx(6),ekn(6,0:nmode),phmx,w
      real*4 mode_phz(3,0:nmode)
c
      ncalc(nmode)=ncalc(nmode) + nctot - nclast
      nclast=nctot
c: Polish root one last time based on first derivative:
cc    kn(nmode)=kx - ekx(1)/ekx(2)
      kn(nmode)=dcmplx(kx,0.d0)
      ekn(1,nmode)=ekx(1)
      ekn(2,nmode)=ekx(2)
      ekn(3,nmode)=ekx(3)
      mode_phz(1,nmode)=phmx
      if(iidiag .ge. 2) then
         print *,nmode,kx,w/kx,phmx,nint(phmx+.25),ncalc(nmode),
     .      ekn(2,nmode)
      endif
      ncalc(nmode+1)=0
c
      return
      end
ccc
      subroutine ek_root(k1,k2,ek1,ek2,kx,ekx,phmx,dk_max,ndv,iipol5)
c
c: Finds roots bracketted by k1 and k2 (k1>k2).  Calls rx_ek_calc
c: asking for ndv derivatives (3 for CW, 6 for BB).
      implicit none
      integer*4 ndv,iipol5,ninv,ncub,nbis,nloop,iibis,nloopm
      real*8 k1,k2,ek1(6),ek2(6),kx,ekx(6),phmx,dk_max,
     .   k1x,k2x,ek1x(6),ek2x(6),dk_err
c
      nloop=0
      nloopm=4
c: Flag to revert to bisection:
      iibis=0
      ninv=0
      ncub=0
      nbis=0
      call ek_copy(k1,k1x,ek1,ek1x,phmx,phmx)
      call ek_copy(k2,k2x,ek2,ek2x,phmx,phmx)
10    continue
      nloop=nloop + 1
      if(nloop .gt. 200) then
         print *,'Failure in ek_root: see Evan ...',k1,k2,ek1,ek2
         print *,'k1x... = ',k1x,k2x,ek1x,ek2x
         stop
      endif
      if(nloop .gt. nloopm) then
cc       print *,'nloop>nloopm in ek_root: ',k1,k2,ek1,ek2
cc       print *,'  k1x,k2x = ',nbis,ncub,k1x,k2x,ek1x,ek2x
         iibis=1
c: Set nloopm so that cubic and bisection are alternated:
         nloopm=nloop + 1
      endif
      if(iibis .eq. 1) then
         kx=0.5d0*(k1x + k2x)
         nbis=nbis + 1
         iibis=0
      elseif(iipol5 .eq. 0) then
         call ek_cub(k1x,k2x,ek1x,ek2x,kx)
         ncub=ncub + 1
      else
         call ek_poly5(k1x,k2x,ek1x,ek2x,kx)
         ncub=ncub + 1
      endif
      call rx_ek_calc(kx,ekx,phmx,ndv)
      dk_err=abs(ekx(1)/ekx(2))
      if(dk_err .gt. dk_max) then
         if(ekx(1)*ek1x(1) .gt. 0.d0) then
            call ek_copy(kx,k1x,ekx,ek1x,phmx,phmx)
         else
            call ek_copy(kx,k2x,ekx,ek2x,phmx,phmx)
         endif
         goto 10
      endif
cc    print *,'ek_root: nloop = ',nloop
cc    print *,'nloop = ',nloop,ncub,nbis,sngl(dk_err)
c
      return
      end
ccc
      subroutine ek_extreme(k1,k2,ek1,ek2,kext,ekext,phmext,
     .   dk_max,iiroots,ndv)
c
      implicit none
      integer*4 iiroots,ndv,nloop,iibis,nloopm
      real*8 k1,k2,ek1(6),ek2(6),kext,ekext(6),dk_max,
     .   phmext,k1x,k2x,ek1x(6),ek2x(6)
cx    real*8 ecross,kcross
c
      iiroots=0
      nloopm=3
      iibis=0
      nloop=0
      call ek_copy(k1,k1x,ek1,ek1x,phmext,phmext)
      call ek_copy(k2,k2x,ek2,ek2x,phmext,phmext)
c
10    continue
      nloop=nloop + 1
      if(nloop .gt. nloopm) then
cc       print *,'nloop > nloopm in ek_ext: ',k1,k2,ek1,ek2
cc       print *,'k1x ... = ',k1x,k2x,ek1x,ek2x,phmext
c: THIS TURNED OUT TO BE UNRELIABLE:
c: Find intersection of tangents to function on either side of extremum:
cx       kcross=(ek2(1)-ek1(1) + ek1(2)*k1 - ek2(2)*k2)/
cx   .      (ek1(2)-ek2(2))
cx       ecross=ek1(1) + ek1(2)*(kcross-k1)
cx       if(ecross*ek1(1) .gt. 0.d0) then
cc          print *,'ek_extreme returning since ecross wrong sign:'
cc          print *,ecross,ek2(1) + ek2(2)*(kcross-k2)
cx          return
cx       endif
         iibis=1
         nloopm=nloop+1
      endif
      if(abs(k2x-k1x) .lt. dk_max) then
cc       print *,'ek_extreme found ext did not cross: ',nloop,
cc   .      k1x,k2x,ek1x,ek2x
         return
      endif
      if(iibis .eq. 0) then
         call ek_cub_ext(k1x,k2x,ek1x,ek2x,kext)
      else
         kext=0.5d0*(k1x + k2x)
         iibis=0
      endif
      call rx_ek_calc(kext,ekext,phmext,ndv)
      if(ekext(1)*ek1(1) .gt. 0.d0) then
         if(ekext(2)*ek1x(2) .gt. 0.d0) then
            call ek_copy(kext,k1x,ekext,ek1x,phmext,phmext)
         else
            call ek_copy(kext,k2x,ekext,ek2x,phmext,phmext)
         endif
         goto 10
      endif
      iiroots=1
c
      return
      end
ccc
      subroutine ek_cub_inv(k1,k2,ek1,ek2,kx)
c
c: Fits a cubic to the function k(e), given two points (k1,ek1) and (k2,ek2)
c: and the derivatives ek1p and ek2p.  Solves the cubic to find the point
c: kx where e(kx)=0.
c
      implicit none
      real*8 k1,k2,ek1(6),ek2(6),kx,k1p,k2p,
     .   a,b,c,d,ec,ecp,kxp,ek0,ek0sq,delek
c
c: Compute derivatives of k(e), rather than e(k):
      k1p=1./ek1(2)
      k2p=1./ek2(2)
c: Fit cubic:
cc    call cub_fit(ek1(1),ek2(1),k1,k2,k1p,k2p,d,c,b,a)
cc    kx=d
      call cub_fit_new(ek1(1),ek2(1),k1,k2,k1p,k2p,d,c,b,a,delek)
      ek0=-ek1(1)/delek
      ek0sq=ek0*ek0
      kx=d + c*ek0 + b*ek0sq + a*ek0*ek0sq
c: Check if cubic is monotonic between y1 and y2:
      ec=-b/(3.*a)
c: If point of inflection of cubic is outside range (y1sh,y2sh), then OK:
c: Assume k1>k2:
      if(ec .lt. 0.d0 .or. ec .gt. 1.d0) return
cc    if(ek1(1) .gt. ek2(1)) then
cc       if(ec .gt. max(0.,ek1(1)) .or. ec .lt. min(0.,ek2(1))) return
cc    else
cc       if(ec .lt. min(0.,ek1(1)) .or. ec .gt. max(0.,ek2(1))) return
cc    endif
c: If derivative at point of inflection same as e1p,e2p, then OK:
      ecp=b*ec + c
cc    kxp=c
      kxp=(c + 2.d0*b*ek0 + 3.d0*a*ek0sq)/delek
      if(ecp*ek1(2) .gt. 0. .and. kxp*ek1(2) .gt. 0) return
c: If not OK, then fit with linear:
      kx=k1 - ek1(1)*(k2-k1)/delek
      print *,'cub_inv_root fit with linear: ',kx
c
      return
      end
ccc
      subroutine ek_cub(k1,k2,ek1,ek2,kx)
c
c: Fits a cubic to the function e(k), given two points (k1,ek1) and (k2,ek2)
c: and the derivatives ek1p and ek2p.  Solves the cubic to find the point
c: kx where e(kx)=0.
c
      implicit none
      integer*4 kcub
      real*8 k1,k2,ek1(6),ek2(6),kx,c1,c2,c3,c4,delk
c
c: Fit cubic:
      call cub_fit_new(k1,k2,ek1(1),ek2(1),ek1(2),ek2(2),
     .   c1,c2,c3,c4,delk)
      call cub_root_new(c1,c2,c3,c4,kx,kcub)
      kx=k1 + kx*delk
cc    if(kcub .eq. 0) print *,'kcub = 0',k1,k2,ek1,ek2,c1,c2,c3,c4,kx
c
      return
      end
ccc
      subroutine ek_poly5(k1,k2,ek1,ek2,kx)
c
c: Fits a cubic to the function e(k), given two points (k1,ek1) and (k2,ek2)
c: and the derivatives ek1p and ek2p.  Solves the cubic to find the point
c: kx where e(kx)=0.
c
      implicit none
cc    integer*4 its
      real*8 k1,k2,ek1(6),ek2(6),kx,delk
      complex*16 c(6),kc
c
c: Fit fifth-order polynomial:
      call poly5_fit(k1,k2,ek1(1),ek2(1),ek1(2),ek2(2),ek1(5),ek2(5),
     .   c(1),c(2),c(3),c(4),c(5),c(6),delk)
c: Make a linear guess for root between 0 and 1:
      kc=dcmplx(-ek1(1)/(ek2(1) - ek1(1)),0.d0)
10    continue
cc    call laguer8(c,5,kc,its)
      call zroot8_int(c,5,kc,0.d0,1.d0)
      kx=k1 + dreal(kc)*delk
c
      return
      end
ccc
      subroutine ek_cub_ext(k1,k2,ek1,ek2,kext)
c
      implicit none
      real*8 k1,k2,ek1(6),ek2(6),kext,c1,c2,c3,c4,ax,bx,cx,q,delk
c
c: Fit cubic:
cc    call cub_fit(k1,k2,ek1(1),ek2(1),ek1(2),ek2(2),a,b,c,d)
      call cub_fit_new(k1,k2,ek1(1),ek2(1),ek1(2),ek2(2),
     .   c1,c2,c3,c4,delk)
      ax=3.*c4
      bx=2.*c3
      cx=c2
      q=-.5d0*(bx + sign(sqrt(bx*bx - 4.d0*ax*cx),bx))
      kext=q/ax
c: Assume k2<k1:
      if(kext .lt. 0.d0 .or. kext .gt. 1.d0) then
         kext=cx/q
      endif
c: Transform from (0,1) to (k1,k2):
      kext=k1 + kext*delk
c
      return
      end
ccc
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
      subroutine cub_fit_old(x1,x2,fx1,fx2,fpx1,fpx2,a,b,c,d)
c
c: This subroutine fits a cubic polynomial "f(x)= a*x**3 + b*x**2 + c*x
c: + d" to the points x1 and x2, their functional values fx1
c: and fx2, and their derivatives fpx1 and fpx2.
c
      implicit none
      real*8 x1,x2,fx1,fx2,fpx1,fpx2,a,b,c,d,xd,xsd,xcd,xmix,
     .   fpd,fd,xfpmix,x1sq,x1cub,x2sq
c
      xd=x1-x2
      x1sq=x1*x1
      x1cub=x1*x1sq
      x2sq=x2*x2
      xsd=x1sq - x2sq
      xcd=x1cub - x2*x2sq
      xmix=x1sq*x2 - x1*x2sq
      fpd=fpx1 - fpx2
      fd=fx1 - fx2
      xfpmix=x2*fpx1 - x1*fpx2
c
      a=(fd-(fpd*xsd)/(2.*xd)+xfpmix)/(xcd-3.*xsd*xsd/(2.*xd)+3.*xmix)
      b=(fpd-3.*a*xsd)/(2.*xd) 
      c=(3.*a*xmix - xfpmix)/xd
      d=fx1-a*x1cub-b*x1sq-c*x1
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
      subroutine cub_root(a,b,c,d,x1,x2,xhit,kcub)
c
c: this subroutine finds the root, xhit, between x1 and x2 of the
c: polynomial f(x)=a*x**3 + b*x**2 + c*x + d.
c
      implicit none
      integer*4 kcub
      real*8 a,b,c,d,x1,x2,xhit,xmin,xmax,a1,a2,a3,a1sq,q,r,
     .   dd,rdd,rq,arg,th,pietth,piefth,one_third
      data pietth/2.09439510239320/,piefth/4.18879020478639/,
     .   one_third/0.33333333333333/
c
      kcub=1
      xmin=min(x1,x2)
      xmax=max(x1,x2)
      a1=b/a
      a2=c/a
      a3=d/a
      a1sq=a1*a1
      q=(3.*a2 - a1sq)/9.
      r=(9.*a1*a2 - 27.*a3 - 2.*a1*a1sq)/54.
      dd=q**3 + r**2
      if(dd .ge. 0.) then
         rdd=sqrt(dd)
         xhit=sign(1.d0,r+rdd)*abs(r+rdd)**(one_third) +
     .      sign(1.d0,r-rdd)*abs(r-rdd)**(one_third) - a1*one_third
      else
         rq=sqrt(-q)
         arg=r/rq**3
         if((arg .lt. -1.) .or. (arg .gt. 1.)) then
            xhit=(x1+x2)/2.
            kcub=0
            return
         endif
         th=acos(arg)/3.
         xhit=2.*rq*cos(th) - a1*one_third
         if((xhit .lt. xmin) .or. (xhit .gt. xmax)) then
            xhit=2.*rq*cos(th + pietth) - a1*one_third
            if((xhit .lt. xmin) .or. (xhit .gt. xmax)) then 
               xhit=2.*rq*cos(th + piefth) - a1*one_third
            endif
         endif
      endif
      if((xhit .le. xmin) .or. (xhit .ge. xmax)) then
         xhit=.5d0*(x1 + x2)
         kcub=0
      endif
c     if(kcub .eq. 0) print *,'kcub=0 in cubroot' 
c
      return
      end 
ccc
      subroutine ek_stack(nstack,k2,ek2,phm2,stack,ii,iival)
c
      implicit none
      integer*4 nstack,ii,iival
      real*8 k2,ek2(6),phm2,stack(6,100)
c
      if(ii .eq. 1) then
cc       if(stack(1,nstack) .eq. k2) then
c: Do not put duplicate values onto stack:
cc          print *,'DUP STACK NOT PUT'
cc          if(iival .ne. 0) stack(6,nstack)=iival
cc       else
            nstack=nstack+1
            stack(1,nstack)=k2
            stack(2,nstack)=ek2(1)
            stack(3,nstack)=ek2(2)
            stack(4,nstack)=ek2(5)
            stack(5,nstack)=phm2
            stack(6,nstack)=iival
cc       endif
      else
         k2=stack(1,nstack)
         ek2(1)=stack(2,nstack)
         ek2(2)=stack(3,nstack)
         ek2(5)=stack(4,nstack)
         phm2=stack(5,nstack)
         iival=nint(stack(6,nstack))
         nstack=nstack-1
      endif
c
      return
      end
ccc
      subroutine rx_poly6(x1,x2,x3,y1,y2,y3,xmm,xpp)
c
c: Solves for the coefficients of a 6th order polynomial that fits
c: y(x) and y'(x) at three points.  x2 is NOT assumed to be the midpoint of
c: x1 and x3.
      implicit none
      integer*4 j,indx(4)
      real*8 x1,x2,x3,y1(2),y2(2),y3(2),xmm(4),a(4,4),d,b(4),
     .   xext,yext,yextp,y1p,y2p,y3p,delx,xsq,xpp(4),xjp1,delx21
      complex*8 c(5),roots(4)
c
c: Form A matrix:
      delx21=x2 - x1
      a(1,1)=delx21*delx21
      a(2,1)=1.d0
      a(3,1)=2.d0*delx21
      a(4,1)=2.d0
      do j=2,4
         a(1,j)=a(1,j-1)*delx21
         a(2,j)=1.d0
         xjp1=dfloat(j+1)
         a(3,j)=xjp1*a(1,j-1)
         a(4,j)=xjp1
      enddo
c: Perform LU decomposition of A
      call ludcmp_r8(a,4,4,indx,d)
c: Normalize derivatives to accout to shifting x1,x3 to 0,1:
      delx=x3 - x1
      y1p=y1(2)*delx
      y2p=y2(2)*delx
      y3p=y3(2)*delx
      b(1)=y2(1) - y1(1) - delx21*y1p
      b(2)=y3(1) - y1(1) - y1p
      b(3)=y2p - y1p
      b(4)=y3p - y1p
c: Do back substitution:
c: Coefficients c3,c4,c5,c6 will be in b(1:4) (c1=y1, c2=y1'):
      call lubksb_r8(a,4,4,indx,b)
cc    xext=1.
cc    yext=y1(1)+y1p*xext+b(1)*xext**2+b(2)*xext**3+b(3)*xext**4+
cc   .   b(4)*xext**5
cc    yextp=y1p+2.*b(1)*xext+3.*b(2)*xext**2+4.*b(3)*xext**3+
cc   .   5.*b(4)*xext**4
cc    print *,'xext,yext,yextp = ',xext,yext,yextp,y3(1),y3p
      c(1)=cmplx(y1p,0.d0)
      c(2)=cmplx(2.*b(1),0.d0)
      c(3)=cmplx(3.*b(2),0.d0)
      c(4)=cmplx(4.*b(3),0.d0)
      c(5)=cmplx(5.*b(4),0.d0)
      call zroots(c,4,roots,.TRUE.)
cc    print *,'roots = ',roots
      do j=1,201
         xext=float(j-1)/200.
         yext=y1(1)+y1p*xext+b(1)*xext**2+b(2)*xext**3+
     .      b(3)*xext**4+b(4)*xext**5
         yextp=y1p+2.*b(1)*xext+3.*b(2)*xext**2+4.*b(3)*xext**3+
     .      5.*b(4)*xext**4
         xext=x1 + (j-1)*delx/200
         write(99,100) xext,yext,yextp,0.
100      format(e14.7,1x,e14.7,1x,e14.7,1x,f10.3)
      enddo
      write(99,100) 0.,0.,0.,0.
      do j=1,4
         xmm(j)=x1 + real(roots(j))*delx
         xsq=roots(j)*roots(j)
         xpp(j)=(2.*b(1) + 6.*b(2)*roots(j) + 12.*b(3)*xsq + 
     .      20.*b(4)*roots(j)*xsq)/(delx*delx)
      enddo
c
      return
      end
ccc
      subroutine rx_poly6_mid(x1,x2,x3,y1,y2,y3,xmm,xpp)
c
c: Solves for the coefficients of a 6th order polynomial that fits
c: y(x) and y'(x) at three points.  x2 is assumed to be the midpoint of 
c: x1 and x3.
      implicit none
      integer*4 j,indx(4)
      real*8 x1,x2,x3,y1(2),y2(2),y3(2),xmm(4),a(4,4),d,b(4),
     .   xext,yext,yextp,y1p,y2p,y3p,delx,xsq,xpp(4)
      complex*8 c(5),roots(4)
c
c: Form A matrix:
      a(1,1)=.25d0
      a(1,2)=.125d0
      a(1,3)=.0625d0
      a(1,4)=.03125d0
      a(2,1)=1.d0
      a(2,2)=1.d0
      a(2,3)=1.d0
      a(2,4)=1.d0
      a(3,1)=1.d0
      a(3,2)=.75d0
      a(3,3)=.5d0
      a(3,4)=.3125d0
      a(4,1)=2.d0
      a(4,2)=3.d0
      a(4,3)=4.d0
      a(4,4)=5.d0
c: Perform LU decomposition of A
      call ludcmp_r8(a,4,4,indx,d)
c: Normalize derivatives to accout to shifting x1,x3 to 0,1:
      delx=x3 - x1
      y1p=y1(2)*delx
      y2p=y2(2)*delx
      y3p=y3(2)*delx
      b(1)=y2(1) - y1(1) - .5d0*y1p
      b(2)=y3(1) - y1(1) - y1p
      b(3)=y2p - y1p
      b(4)=y3p - y1p
c: Do back substitution:
c: Coefficients c3,c4,c5,c6 will be in b(1:4) (c1=y1, c2=y1'):
      call lubksb_r8(a,4,4,indx,b)
cc    xext=1.
cc    yext=y1(1)+y1p*xext+b(1)*xext**2+b(2)*xext**3+b(3)*xext**4+
cc   .   b(4)*xext**5
cc    yextp=y1p+2.*b(1)*xext+3.*b(2)*xext**2+4.*b(3)*xext**3+
cc   .   5.*b(4)*xext**4
cc    print *,'xext,yext,yextp = ',xext,yext,yextp,y3(1),y3p
      c(1)=cmplx(y1p,0.d0)
      c(2)=cmplx(2.*b(1),0.d0)
      c(3)=cmplx(3.*b(2),0.d0)
      c(4)=cmplx(4.*b(3),0.d0)
      c(5)=cmplx(5.*b(4),0.d0)
      call zroots(c,4,roots,.TRUE.)
cc    print *,'roots = ',roots
      do j=1,201
         xext=float(j-1)/200.
         yext=y1(1)+y1p*xext+b(1)*xext**2+b(2)*xext**3+
     .      b(3)*xext**4+b(4)*xext**5
         yextp=y1p+2.*b(1)*xext+3.*b(2)*xext**2+4.*b(3)*xext**3+
     .      5.*b(4)*xext**4
         xext=x1 + (j-1)*(x3-x1)/200
         write(99,100) xext,yext,yextp,0.
100      format(e14.7,1x,e14.7,1x,e14.7,1x,f10.3)
      enddo
      write(99,100) 0.,0.,0.,0.
      do j=1,4
         xmm(j)=x1 + real(roots(j))*(x3 - x1)
         xsq=roots(j)*roots(j)
         xpp(j)=(2.*b(1) + 6.*b(2)*roots(j) + 12.*b(3)*xsq + 
     .      20.*b(4)*roots(j)*xsq)/(delx*delx)
      enddo
c
      return
      end
ccc
      subroutine rat_fun_6fit(x1,x2,y1,y2,y1px,y2px,y1ppx,y2ppx,
     .   p0,p1,p2,p3,q1,q2,delx)
c
c: Solves for coefficients of rational function y=P(x)/Q(x), where 
c: P(x)=p0 + p1*x + p2*x^2 + p3*x^3, Q(x)=1 + q1*x + q2*x^2,
c: given the function values y1,y2 at x1,x2, the derivatives y1p,y2p,
c: and the second derivatives y1pp,y2pp. The interval (x1,x2) is
c: normalized to (0,1), so yp computed from the coeffs must be divided
c: by delx, and ypp must be divided by delx*delx.
      implicit none
      integer*4 indx(5)
      real*8 x1,x2,y1,y2,y1px,y2px,y1ppx,y2ppx,p0,p1,p2,p3,q1,q2,
     .   a(5,5),delx,delxsq,y1p,y2p,y1pp,y2pp,b(5),d
c
c: Constants:
      a(1,1)=1.d0
      a(1,2)=0.d0
      a(1,3)=0.d0
      a(1,5)=0.d0
      a(2,1)=0.d0
      a(2,2)=2.d0
      a(2,3)=0.d0
      a(3,1)=1.d0
      a(3,2)=1.d0
      a(3,3)=1.d0
      a(4,1)=1.d0
      a(4,2)=2.d0
      a(4,3)=3.d0
      a(5,1)=0.d0
      a(5,2)=2.d0
      a(5,3)=6.d0
c
      delx=x2 - x1
      delxsq=delx*delx
      y1p=y1px*delx
      y2p=y2px*delx
      y1pp=y1ppx*delxsq
      y2pp=y2ppx*delxsq
      a(1,4)=-y1
      a(2,4)=-2.d0*y1p
      a(3,4)=-y2
      a(4,4)=-(y2 + y2p)
      a(5,4)=-(2.d0*y2p + y2pp)
      a(2,5)=-2.d0*y1
      a(3,5)=-y2
      a(4,5)=-(2.d0*y2 + y2p)
      a(5,5)=-(2.d0*y2  + 4.d0*y2p + y2pp)
c: Perform LU decomposition of A
      call ludcmp_r8(a,5,5,indx,d)
c: Do back substitution:
c: Coefficients p1,p2,p3,q1,q2 will be in b(1:5) (p0=y1):
      b(1)=y1p
      b(2)=y1pp
      b(3)=y2 - y1
      b(4)=y2p
      b(5)=y2pp
      call lubksb_r8(a,5,5,indx,b)
c: Pull out coefficients:
      p0=y1
      p1=b(1)
      p2=b(2)
      p3=b(3)
      q1=b(4)
      q2=b(5)
c
      return
      end
ccc
      subroutine rat_fun_6eval(p0,p1,p2,p3,q1,q2,delx,x,y,yp,ndv)
c
c: x should be between 0 and 1 for interpolation.
      implicit none
      integer*4 ndv
      real*8 p0,p1,p2,p3,q1,q2,delx,x,y,yp,xsq,P,Q,Pp,Qp
c
      xsq=x*x
      P=p0 + p1*x + p2*xsq + p3*xsq*x
      Q=1.d0 + q1*x + q2*xsq
      y=P/Q
      if(ndv .gt. 1) then
         Pp=p1 + 2.d0*p2*x + 3.d0*p3*xsq
         Qp=q1 + 2.d0*q2*x
         yp=(Pp - Qp*y)/(Q*delx)
      endif
c
      return
      end
ccc
      subroutine poly5_fit(x1,x2,y1,y2,y1px,y2px,y1ppx,y2ppx,
     .   c1,c2,c3,c4,c5,c6,delx)
c
c: Solves for coefficients of polynomial 
c: y=c1 + c2*x + c3*x^2 + c4*x^3 + c5*x^4 + c6*x^5, 
c: given the function values y1,y2 at x1,x2, the derivatives y1p,y2p,
c: and the second derivatives y1pp,y2pp. The interval (x1,x2) is
c: normalized to (0,1), so yp computed from the coeffs must be divided
c: by delx, and ypp must be divided by delx*delx.
      implicit none
      real*8 x1,x2,y1,y2,y1px,y2px,y1ppx,y2ppx,c1,c2,c3,c4,c5,c6,
     .   delx,delxsq,y1p,y2p,y1pp,y2pp,b1,b2,b3
c
      delx=x2 - x1
      delxsq=delx*delx
      y1p=y1px*delx
      y2p=y2px*delx
      y1pp=y1ppx*delxsq
      y2pp=y2ppx*delxsq
c
      c1=y1
      c2=y1p
      c3=0.5d0*y1pp
      b1=y2 - c1 - c2 - c3
      b2=y2p - c2 - 2.d0*c3
      b3=y2pp - 2.d0*c3
      c4=10.d0*b1 - 4.d0*b2 + 0.5d0*b3
      c5=-15.d0*b1 + 7.d0*b2 - b3
      c6=6.d0*b1 - 3.d0*b2 + 0.5d0*b3
c
      return
      end
ccc
      subroutine poly5_eval(c1,c2,c3,c4,c5,c6,delx,x,y,yp,ypp,x2,x3)
c
c: Evaluates the polynomial 
c: y=c1 + c2*x + c3*x^2 + c4*x^3 + c5*x^4 + c6*x^5 and its derivative.
      implicit none
      real*8 c1,c2,c3,c4,c5,c6,delx,x,y,yp,ypp,x2,x3,x4,x5
c
      x2=x*x
      x3=x*x2
      x4=x*x3
      x5=x*x4
      y=c1 + c2*x + c3*x2 + c4*x3 + c5*x4 + c6*x5
      yp=(c2 + 2.d0*c3*x + 3.d0*c4*x2 + 4.d0*c5*x3 + 5.d0*c6*x4)/delx
      ypp=(2.d0*c3 + 6.d0*c4*x + 12.d0*c5*x2 + 20.d0*c6*x3)/(delx*delx)
c
      return
      end
ccc
      subroutine ek_miss2_safe(k1n,k2n,ek1n,ek2n,phm1n,phm2n,
     .   k1,k2,ek1,ek2,phm1,phm2,kx,ekx,phmx,dk_max,ndv0,
     .   nloopm,iimok,iidiag)
c
      implicit none
      integer*4 ndv0,iidiag,kord(2049),nloop,nax,jax,jprev,
     .   j1,j2,iiroots,nloopm,iimok,nprobe
      real*8 k1n,k2n,ek1n(6),ek2n(6),phm1n,phm2n,k1,k2,ek1(6),ek2(6),
     .   phm1,phm2,kx,ekx(6),phmx,dk_max,kax(2049),ekax(3,2049),
     .   phax(2049),k0,delk,dk_lim
c
      dk_lim=100.d0*dk_max
      iimok=1
      call ek_copy(k2n,kax(1),ek2n,ekax(1,1),phm2n,phax(1))
      call ek_copy(k1n,kax(2),ek1n,ekax(1,2),phm1n,phax(2))
c: kord(j) holds the index of the next value of k in kax:
      kord(1)=2
      kord(2)=0
      nax=2
      delk=k1n - k2n
      do nloop=0,nloopm
         nprobe=2**nloop
         jprev=1
         k0=k2n - .5d0*delk
c: Sample interval more finely:
         do jax=nax+1,nax+nprobe
            k0=k0 + delk
            kax(jax)=k0
            call rx_ek_calc(k0,ekax(1,jax),phax(jax),ndv0)
            if(ekax(1,jax)*ek1n(2) .gt. 0.d0) then
c: If sign of e is as desired, copy to output and return:
               call ek_copy(k0,kx,ekax(1,jax),ekx,phax(jax),phmx)
               j1=jprev
               j2=kord(jprev)
               call ek_copy(kax(j2),k1,ekax(1,j2),ek1,phax(j2),phm1)
               call ek_copy(kax(j1),k2,ekax(1,j1),ek2,phax(j1),phm2)
               return
            endif
            kord(jax)=kord(jprev)
            kord(jprev)=jax
            jprev=kord(jax)
         enddo
         delk=.5*delk
         nax=nax + nprobe
c: Check for correct change in sign of derivatives between neighboring pts:
         if(ek1n(2) .gt. 0.d0) then
            j1=1
            do jax=1,nax-1
               j2=kord(j1)
               if(ekax(2,j1).gt.0.d0 .and. ekax(2,j2).lt.0.d0) then
c: Now find extr between kax(j1) and kax(j2) and check if on other side of 0:
                  call ek_extreme(kax(j2),kax(j1),ekax(1,j2),ekax(1,j1),
     .               kx,ekx,phmx,dk_max,iiroots,ndv0)
                  if(iiroots .ne. 0) then
                     goto 10
                  else
                     if(iidiag .ge. 2) print *,
     .                  'Info msg: iiroots=0 from ek_miss2_safe '
                  endif
               endif
               j1=j2
            enddo
         else
            j1=1
            do jax=1,nax-1
               j2=kord(j1)
               if(ekax(2,j1).lt.0.d0 .and. ekax(2,j2).gt.0.d0) then
c: Now find extr between kax(j1) and kax(j2) and check if on other side of 0:
                  call ek_extreme(kax(j2),kax(j1),ekax(1,j2),ekax(1,j1),
     .               kx,ekx,phmx,dk_max,iiroots,ndv0)
                  if(iiroots .ne. 0) then
                     goto 10
                  else
                     if(iidiag .ge. 2) print *,
     .                  'Info msg: iiroots=0 from ek_miss2_safe '
                  endif
               endif
               j1=j2
            enddo
         endif
c: If delk gets very small, exit early:
         if(delk .lt. dk_lim) goto 20
      enddo
20    continue
      if(iidiag .ge. 1) then
         print *,'ek_miss2_safe could not find point: tell EKW.'
      endif
      iimok=0
      return
10    continue
c
c: Now find extremum between kax(j1) and kax(j2):
      call ek_copy(kax(j2),k1,ekax(1,j2),ek1,phax(j2),phm1)
      call ek_copy(kax(j1),k2,ekax(1,j1),ek2,phax(j1),phm2)
c
      return
      end
ccc
      subroutine ek_copy(kin,kout,ekin,ekout,phmin,phmout)
c
      implicit none
      real*8 kin,kout,ekin(6),ekout(6),phmin,phmout
c
      kout=kin
      ekout(1)=ekin(1)
      ekout(2)=ekin(2)
cc    ekout(5)=ekin(5)
      phmout=phmin
c
      return
      end
ccc
      subroutine rx_freq_init(kmax,kmin)
c
c: Initializes variables associated with a new frequency.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 j
      real*8 kmin,kmax,nepcon,nepfac,hj,dk_est,kminx,w0
      data nepcon/54575.05415d0/,w0/6.283185307179586d+03/
c
      w=2.d0*pie*f_hz
      wsq=w*w
      kcrmin=w/crmax
c
      if(nthorpe .gt. 0) call thorpe_attn
c
      do j=2,nlay-1
         call rx_k_attn(w,geo(1,1,j),geo(1,4,j),xk(1,j),xksq(1,j))
         call rx_k_attn(w,geo(2,1,j),geo(2,4,j),xk(2,j),xksq(2,j))
         call rx_eta_calc(xksq(1,j),xksq(2,j),h(j),eta(j),etasq(j),
     .      isp(j))
co		modified to bypass compiler warning
co         call rx_eta_calc(dreal(xksq(1,j)),dreal(xksq(2,j)),
co     .		h(j),eta(j),etasq(j),isp(j))				

      enddo
      call rx_k_attn(w,geo(1,1,nlay),geo(1,4,nlay),xk(1,nlay),
     .   xksq(1,nlay))
c: Set xk in vacuum to xk in first layer for convenience:
      xk(1,1)=xk(1,2)
      xk(2,1)=xk(1,2)
      xksq(1,1)=xksq(1,2)
      xksq(2,1)=xksq(1,2)
c
      kw0=w/cfmin
cc    xkref=xk(isvmin,nsvmin)
c: Set branch point ratio for reference depth:
cc    xkrat_ref(1)=xkref/kw0
cc    kw=dreal(xkref)
cc    cref=w/kw0
cc    errdkms=dmin1(errdkms,(1.d-5*kw)**2)
cc    errdk100=-100.d0*dsqrt(errdkms)
c
cc    nepfac=w/nepcon
c: New 3-26-99. Incorporate power law for attenuation vs freq such that
c: when fexp=1, we get alpha_w=alpha * w (as before), and when w=w0 (1000 Hz),
c: alpha_w=alpha * w0 (as before):
c
c: Set up additional arrays for perturbation calc of vg, attenuation:
cc    nlatt=0
      do j=2,nlay
         hj=dmax1(1.d-100,h(j))
c: Not needed since vg is computed from derivatives:
         am(j)=1.d0/(geo(1,1,j)*geo(1,1,j))
         bm(j)=(1.d0/(geo(2,1,j)*geo(2,1,j)) - am(j))/hj
c: Bottom line is: [w]  ==> [(w/w0)^fexp * w0]
         nepfac=((w/w0)**fexp(1,j)) * w0 / nepcon
         betm(j)=(geo(1,4,j)*nepfac)/geo(1,1,j)
         gm(j)=((geo(2,4,j)*nepfac)/geo(2,1,j) - betm(j))/hj
c: 1-16-96: Better agreement when surface density used in attenuation calcs:
ccx      rhom(j)=0.5d0*(geo(1,3,j) + geo(2,3,j))
         rhom(j)=geo(1,3,j)
      enddo
c
c: Find horizontal wavenumber k interval over which to find modes:
      dk_est=.0001/zdep(nlay-1)
      kmax=kw0 - dk_est
      kmin=w/(0.99999d0*geo(1,1,nlay))
      if(cphmax .gt. 0.) then
         kminx=w/cphmax
         if(kminx .lt. kmin) then
            print *,'Requested CPHMAX > cp in lower h-space. ',
     .         'Consider using false bottom (iifb).'
            write(2,*) 'Requested CPHMAX > cp in lower h-space. ',
     .         'Consider using false bottom (iifb).'
         endif
         kmin=max(kminx,kmin)
      endif
c
      return
      end
ccc
      subroutine rx_k_attn(w,c0,alpha,xk,xksq)
c
c: Computes complex wavenumber for medium with sound speed c0 and
c: attenuation alpha (in dB/(m-kHz)).
      implicit none
      complex*16 xk,xksq
      real*8 w,c0,xkr,alpha,fac
c: fac=2*pi*1000*20*log(e):
      data fac/5.457505415367364d+04/
c
      xkr=w/c0
      xk=dcmplx(xkr,0.d0)
      xksq=dcmplx(xkr*xkr,0.d0)
c
      return
      end
ccc
      subroutine rx_eta_calc(xk1sq,xk2sq,h,eta,etasq,iso)
c
      implicit none
      integer*4 iso
      complex*16 eta,etasq
      real*8 xk1sq,xk2sq,dk,h,third,etar
      data third/0.333333333333333d0/
c
      if(iso .ne. 0 .or. h .le. 0.d0) then
         eta=(0.d0,0.d0)
         etasq=(0.d0,0.d0)
         return
      endif
      dk=xk2sq - xk1sq
      if(dk .gt. 0.d0) then
         etar=(dk/h)**third
      else
c: For dk/h in left half plane, make sure eta lies close to neg
c: real axis:
         etar=-((-dk/h)**third)
      endif
      eta=dcmplx(etar,0.d0)
      etasq=dcmplx(etar*etar,0.d0)
c
      return
      end
ccc
      subroutine rx_prep
c
c: Makes a few additional computations required by rx routines.
c
      use parms_com
      use i_o_com
      use gen_com
c
      integer*4 j,jd
      real*8 c3,c4,ccr_min,c_up,c_upx,c_dn
      complex*16 zzero
      data zzero/(0.d0,0.d0)/
c
      if(allf(1) .ne. 1 .or. allf(2) .ne. 1) then
         print *,'NOTE: Real axis option chosen for environment '//
     .      'with elastic properties.'
         print *,'      Shear wave speeds will be ignored.'
         write(2,*) 'NOTE: Real axis option chosen for environment '//
     .      'with elastic properties.'
         write(2,*) '      Shear wave speeds will be ignored.'
      endif
c: Redo isospeed calculation, not taking into account changes in attenuation:
      do j=2,nlay
         isp(j)=0
         if(geo(1,1,j) .eq. geo(2,1,j)) isp(j)=1
      enddo
c
c: Make rmax at least 100 km since not finding k exactly causes problems
c: in mode function finding:
      if(iicw .eq. 2) then
         dk_max=2.d0*pie/(200.d0*1000.d0*amax1(10000.e0,rmax))
      else
         dk_max=2.d0*pie/(200.d0*1000.d0*amax1(100.e0,rmax))
      endif
      dk_max0=2.d0*pie/(200.d0*1000.d0*rmax)
c
c: Set rhofac at surface to 1 for convenience:
      rhofac(1)=1.d0
      rho_prod(1,nlay)=1.d0
      do j=nlay-1,2,-1
         rhofac(j)=geo(2,3,j)/geo(1,3,j+1)
         rho_prod(1,j)=rho_prod(1,j+1)*rhofac(j)
      enddo
c: For convenience, make density of upper halfspace (assumed vaccuum) 
c: same as water:
      geo(2,3,1)=geo(1,3,2)
      rho_prod(1,1)=rho_prod(1,2)
c: In downward direction:
      rho_prod(2,1)=1.d0
      do j=2,nlay
         rho_prod(2,j)=rho_prod(2,j-1)/rhofac(j-1)
      enddo
c: Set phi,dphi at bottom of lower halfspace to zero:
      philay(2,nlay)=zzero
      dphilay(2,nlay)=zzero
      Aplay(2,nlay)=-300.d0
c
c: Set up mapping from layer number to duct number:
      ndrx=0
      c_upx=geo(2,1,1)
      do j=2,nlay-1
c: c3,c4 are sound speeds at top and bottom of current layer:
         c3=geo(1,1,j)
         c4=geo(2,1,j)
c
c: c_up,c_dn are sound speed above and below that matter:
         c_up=geo(2,1,j-1)
         if(c3 .eq. geo(2,1,j-1)) c_up=geo(1,1,j-1)
         if(j .eq. 2) c_up=1.d10
         c_dn=geo(1,1,j+1)
         if(c4 .eq. geo(1,1,j+1)) c_dn=geo(2,1,j+1)
         if(j .eq. nlay-1) c_dn=1.d10
c
         if(c3 .eq. c4) then
c: Isospeed: c_up,c_dn are speeds to compare to:
            if((c3 .lt. c_up .or. 
c: Case of two isospeeds together:
     .         (c3 .eq. c_up .and. (c3 .lt. c_upx)))
     .         .and. (c4 .lt. c_dn)) then
               call rx_duct(ndrx,j,1,jval,c3,cduct)
            endif
            if(c3 .ne. c_up) then
               c_upx=c_up
            endif
         elseif(c3 .lt. c4) then
c: Positive gradient with depth, check for duct at top of layer:
            if(c3 .lt. c_up .or. (c3 .eq. c_up .and. c3 .lt. c_upx)) 
     .         then
               call rx_duct(ndrx,j,1,jval,c3,cduct)
            endif
            c_upx=c_up
         elseif(c4 .lt. c3) then
c: Negative gradient with depth, check for duct at bottom of layer:
            if(c4 .lt. c_dn) then
               call rx_duct(ndrx,j,2,jval,c4,cduct)
            endif
            c_upx=c_up
         endif
      enddo
c: Find phase speed interval over which modes may cross in multiple duct envs:
      if(ndrx .gt. 1) then
         ccr_hi=0.d0
c: ccr_hi is largest sound speed between ducts:
         do jd=jval(1,1),jval(1,ndrx)-1
            ccr_hi=max(ccr_hi,geo(2,1,jd))
            ccr_hi=max(ccr_hi,geo(1,1,jd+1))
         enddo
c: ccr_lo is second highest duct sound speed:
         ccr_min=cduct(1)
         do jd=2,ndrx
            ccr_min=min(ccr_min,cduct(jd))
         enddo
         ccr_lo=ccr_hi
         do jd=1,ndrx
            if(cduct(jd) .ne. ccr_min) ccr_lo=min(ccr_lo,cduct(jd))
         enddo
      endif
      if (ccr_hi .gt. geo(1,1,nlay)) then
         print *,'Lower halfspace sound speed must be largest '//
     .        'speed in profile!! ', geo(1,1,nlay), ccr_hi
         stop
      endif
cc    print *,'jl2jd = ',(jl2jd(j),j=2,nlay)
cc    print *,'jpeak = ',(jpeak(j),j=1,ndrx)
c
      return
      end
ccc
      subroutine rx_duct(ndrx,jlay,ii,jval,cd,cduct)
c
      implicit none
      integer*4 ndrx,jlay,ii,jval(2,ndrx+1)
      real*8 cd,cduct(ndrx+1)
c
      ndrx=ndrx + 1
      jval(1,ndrx)=jlay
      jval(2,ndrx)=ii
      cduct(ndrx)=cd
c
      return
      end
c
      subroutine rx_zmx_init
c
c: Initializes depth arrays for computation of mode functions.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jptr,jj,j,jj1,jj2,jsr
      real*8 z1,z2,rhfac,c1inv,c1invsq,c2inv,c2invsq,cfac,zpar
c
      jptr=0
      jj=1 
      z1=-1.d+20
c
c: For layers below reference depth (force reference depth to ocean surface):
c: nlwz = number of layers with zs or zr in them. jlwz holds the layer #s:
      nlwz=0
      do j=1,nlay
         z2=zdep(j)
         if(j .ge. jsurf) z2=z2 + 1.d-10*abs(z2)
         do while(zsr(jj) .lt. z1 .and. jj .lt. nzsr)
            jj=jj+1
         enddo
         if(jj .eq. nzsr) jj=jj+1
         jj1=jj
         do while(zsr(jj) .lt. z2 .and. jj .lt. nzsr)
            jj=jj+1
         enddo
         if(jj .eq. nzsr) jj=jj+1
         jj2=jj-1
c        print *,'jlay,jj1,jj2,nzsr = ',jlay,jj1,jj2,nzsr
         jzmx(j)=jptr+1
         nzmx(j)=jj2-jj1+1
         if(nzmx(j) .gt. 0) then
            nlwz=nlwz + 1
            jlwz(nlwz)=j
         endif
         c1inv=1.d0/geo(1,1,j)
         c1invsq=c1inv*c1inv
         c2inv=1.d0/geo(2,1,j)
         c2invsq=c2inv*c2inv
         cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
         rhfac=(geo(2,3,j) - geo(1,3,j))/dmax1(1.d-100,h(j))
c: iia1b1=1 means compute A1,B1 in mode_fun at top:
         iia1b1(j)=1
c: iia1b2=2 means compute A1,B1 at bottom (force ii=1 for lower h-space):
         if(geo(1,1,j) .le. geo(2,1,j) .and. j .ne. nlay) iia1b1(j)=2
         do jsr=jj1,jj2
            jptr=jptr + 1
c: zmx measured from top of layer for c1>c2, from bottom for c1<c2:
c: NEW: zmx measured from top of layer always:
            zmx(jptr)=zsr(jsr)-zdep(j-1)
c: jsrmx(jptr)=jsr, the index of depth zmx(jptr) in the array zsr(1:nzsr):
            jsrmx(jptr)=jsr
            zpar=zsr(jsr) - zdep(j-1)
            rho_sr(jsr)=geo(1,3,j) + zpar*rhfac
            cp_sr(jsr)=1./sqrt(c1invsq + zpar*cfac)
            cs_sr(jsr)=0.
c: Pointer from zmx to zm [find depth zsr(jsr) in zmx(mx_m(jsr))]:
            mx_m(jsr)=jptr
c: Pointer from zsr to geo [find s/r depth zsr(jsr) in geo(:,:,j)]:
            jsr2j(jsr)=j
         enddo
         z1=z2
      enddo
c
c: Total number of depths in zmx, which includes source, receiver, and
c: layer interface depths:
      nzmxtot=jptr
c
c: Find maximum receiver sound speed for use in mode_field:
      crmax=cfmin
      do jsr=1,nzsr
         crmax=dmax1(crmax,cp_sr(jsr))
      enddo
c
      return
      end
