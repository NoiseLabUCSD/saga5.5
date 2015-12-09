      subroutine vary_env
c
c: Calls mode program varying environmental parameters.
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jrun,jparm,kobt,kparm,klay,k_tb,kgrad,iiw,
     .   nrunx,joff,nrow,jf,idum
      real*8 dxparm,xparm,fac
      real*4 ran2
c
      nrunx=iabs(nrun)
      joff=0
      iiwrite=1
c
      if(rseed .ne. 0) idum=-iabs(rseed)
      ncall=0
      allocate(var_ax(nrunx))
      do jrun=1,nrunx
         iiw=iiwrite
         iiwrite=0
         if(jrun .gt. 1) then
            call svp_read
         endif
         iiwrite=iiw
         do jparm=1,nparm
            kobt=kvar(1,jparm)
            klay=kvar(2,jparm)
            k_tb=iabs(kvar(3,jparm))
            kgrad=isign(1,kvar(3,jparm))
            kparm=kvar(4,jparm)
            if(rseed .ne. 0) then
               fac=ran2(idum)
            else
               fac=float(jrun-1)/float(max(1,nrunx-1))
            endif
            dxparm=xvar(2,jparm) - xvar(1,jparm)
            xparm=xvar(1,jparm) + dxparm*fac
            if(iidiag .ge. 1) then
               print *,'jrun,jparm,xparm = ',jrun,jparm,xparm
            endif
            call vary_parm(kobt,klay,k_tb,kgrad,kparm,xparm)
         enddo
         var_ax(jrun)=xparm
         call svp_check(1)
         call svp_check2
         if(iicw .eq. 1) then
            call cw_modes(1,1,jrun,nrunx)
         elseif(iicw .eq. 2) then
            if(iirx .le. 0) then
               call bb_brute
cc          elseif(iirx .eq. -1) then
cc             call bb_modes
            elseif(iirx .eq. 1) then
               call rx_bb
            elseif(iirx .eq. 2) then
               call rx_bb_brute
            endif
         endif
         if(iikpl .ne. 0) then
c: Complex k-plane computations:
            call k_plane(1,jrun,nrunx)
         endif
c
         if(iirc .ne. 0) then
c: Plane wave reflection coefficient computations:
            call pw_refco(1,jrun,nrunx)
         endif
      enddo
c
      deallocate(var_ax)
c
      return
      end
ccc
      subroutine vary_parm(kobt,klay,k_tb,kgrad,kparm,xparm)
c
c: Translates kobt,klay,k_tb,kgrad,kparm and changes appropriate parameter in
c: input geoacoustic parameter arrays.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 kobt,klay,k_tb,kgrad,kparm
      real*8 xparm,delc
c
c: Convert to positive layer if a top layer:
      if(kobt .eq. 0) then
         if(kparm .eq. 0) then
c: For ocean (kobt=0) and thickness (kparm=0), vary ocean DEPTH, regardless
c: of layer number:
            call vary_ocean_depth(xparm)
         elseif(kparm .eq. 1) then
            if(kgrad .eq. -1) then
c: Change other part of layer (top/bottom) first:
               delc=xparm - csvp(klay+k_tb-1)
               csvp(klay+2-k_tb)=csvp(klay+2-k_tb) + delc
            endif
            csvp(klay+k_tb-1)=xparm
         else
            print *,'Illegal pc for obt=0 in vary_env (0 or 1): ',
     .         kobt,kparm
            stop
         endif
      elseif(kobt .eq. 1) then
         call vary_lay(klay,k_tb,kgrad,kparm,xparm,hb,geob,bpb)
      elseif(kobt .eq. 2) then
c: For top layers, add 1 to klay since upper halfspace is layer 1:
         call vary_lay(klay+1,k_tb,kgrad,kparm,xparm,ht,geot,bpt)
      elseif(kobt .eq. 3) then
         rkm(1)=xparm
      endif
      if(iidiag .eq. 1) print *,'parm = ',kparm,klay,k_tb,kgrad,xparm
c
      return
      end
ccc
      subroutine vary_lay(klay,k_tb,kgrad,kparm,xparm,hb,geob,bpb)
c
      implicit none
      integer*4 klay,k_tb,kgrad,kparm,jp
      real*8 xparm,hb(klay),geob(2,5,klay),grad,hold,bpb(4,klay)
c
      if(kparm .eq. 0) then
         hold=hb(klay)
         hb(klay)=xparm
c: When changing thickness of layer, kgrad=-1 means keep gradient same:
         if(kgrad .eq. -1) then
            do jp=1,5
               grad=(geob(2,jp,klay)-geob(1,jp,klay))/max(1.d-100,hold)
               geob(2,jp,klay)=geob(1,jp,klay) + grad*hb(klay)
            enddo
         endif
c: kparm=6,7 means to change fexpp,fexps:
      elseif(kparm .eq. 6) then
         bpb(3,klay)=xparm
      elseif(kparm .eq. 7) then
         bpb(4,klay)=xparm
      else
         if(kgrad .eq. -1) then
            grad=(geob(3-k_tb,kparm,klay)-geob(k_tb,kparm,klay))/
     .         max(1.d-100,hb(klay))
            geob(k_tb,kparm,klay)=xparm
            geob(3-k_tb,kparm,klay)=xparm + grad*hb(klay)
         else
            geob(k_tb,kparm,klay)=xparm
         endif
      endif
c
      return
      end
ccc
      subroutine vary_xtract(kobt,klay,k_tb,kparm,xparm)
c
c: Obtains current value in the input geoacoustic profile corresponding
c: to the variable given by kobt,klay,k_tb,kparm.
c
      use parms_com
      use i_o_com
      integer*4 kobt,klay,k_tb,kparm
      real*8 xparm
c
      if(kobt .eq. 0) then
         if(kparm .eq. 0) then
            xparm=zsvp(nsvp)
         elseif(kparm .eq. 1) then
            xparm=csvp(klay)
         else
            print *,'Illegal kparm for obt=0: ',kparm
            stop
         endif
      elseif(kobt .eq. 1) then
         if(kparm .eq. 0) then
            xparm=hb(klay)
         else
            xparm=geob(k_tb,kparm,klay)
         endif
      elseif(kobt .eq. 2) then
         if(kparm .eq. 0) then
            xparm=ht(klay)
         else
            xparm=geot(k_tb,kparm,klay)
         endif
      elseif(kobt .eq. 3) then
         xparm=rkm(1)
      endif
c
      return
      end
ccc
      subroutine vary_ocean_depth(bathy)
c
c: Changes ocean depth from zsvp(nsvp) to bathy.
c
      use parms_com
      use i_o_com
      integer*4 jj,jsvp
      real*8 bathy,cpt,alppt,xppt,rhopt
c
      do jj=nsvp,1,-1
         jsvp=jj
         if(zsvp(jsvp) .lt. bathy) goto 99
      enddo
99    if(jsvp .eq. nsvp) jsvp=nsvp-1
      xppt=bathy
      call c_interp(zsvp(jsvp),zsvp(jsvp+1),csvp(jsvp),
     .   csvp(jsvp+1),alpha_svp,alpha_svp,rho_svp,rho_svp,
     .   xppt,cpt,alppt,rhopt,1)
      csvp(jsvp+1)=cpt
      zsvp(jsvp+1)=bathy
      nsvp=jsvp+1
c
      return
      end
