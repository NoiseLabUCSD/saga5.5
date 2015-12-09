c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
      subroutine svp_check
c
      implicit none
      include 'Parms_com'
      include 'i_o_svp_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      common /out_com3/ nzsr, nlay, NSCTOR,wclmn
      integer*4 nzsr, nlay, NSCTOR,wclmn
      include 'gen1_com'
c: Local variables:
      integer*4 j,ii,jflubot,jflutop,nadd,ji,j1,j2,jj,jtop1,
     .          inc,ii2
      real*8 mmx(5),htop
      data mmx/.05,.05,.001,.0005,.0005/
c
c: Find the number of fluid layers in bottom and top layering:
      call last_flu(nlayt+1,1,-1,geot,jflutop,allf(2))
      call last_flu(1,nlayb+1,1,geob,jflubot,allf(1))
c
c: NOTE: TOP LAYERS ASSUMED TO BE GIVEN FROM UPPER HALFSPACE TO OCEAN,
c: WHILE BOTTOM LAYERS ARE GIVEN FROM OCEAN BOTTOM TO LOWER HALFSPACE.
c
      jlfake(1)=0
      jlfake(2)=0
c: Copy top layering to master arrays:
      nlay=0
c: Set thickness of halfspaces to large numbers (for use in zmx_init):
      ht(1)=1.d+20
      hb(nlayb+1)=1.d+20
c: Make type of profile linear for halfspaces:
      ktt(1)=1
      ktb(nlayb+1)=1
      htop=0.d0
      jtop1=1
      do j=1,nlayt+1
         if(j .gt. jtop1) htop=htop + ht(j)
c: Skip a layer between fluid and solid layers in case seismic modes
c: are desired (see duct_check):
         if(j .eq. jflutop) then
            nlay=nlay+1
            h(nlay)=0.
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
      do j=1,nlayb+1
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
      do j=1,5
         geo(1,j,1)=geo(2,j,1)
         geo(2,j,nlay)=geo(1,j,nlay)
      enddo
c
c: Compute absolute depths of layer interfaces (ocean surface is at z=0):
c: Note that zdep(j) is depth of BOTTOM of layer j.
      zdep(1)=-htop
      do j=2,nlay
         zdep(j)=zdep(j-1) + h(j)
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
c: Make fake layers transparent for now:
      do ii=1,2
         j1=jlfake(ii)
         if(j1 .ne. 0) then
            ii2=3-ii
            inc=2*ii - 3
            j2=j1 + inc
            do jj=1,5
               geo(1,jj,j1)=geo(ii2,jj,j2)
               geo(2,jj,j1)=geo(1,jj,j1)
            enddo
         endif
      enddo
c
      return
      end
ccc
      subroutine svp_check2(iifail)
c
      implicit none
      include 'Parms_com'
      include 'i_o_svp_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      common /out_com3/ nzsr, nlay, NSCTOR,wclmn
      integer*4 nzsr, nlay, NSCTOR,wclmn,iifail
      include 'gen1_com'
c: Local variables:
      integer*4 j,ii,ii2,nhsp,j1,iibad,jj,mm_dn,mm_up,mm_upx,
     .          jmin,jjmin,jflu1,jflu2,nd_next,jlast
      real*8 c1,c2,c3,c4,cpint,c_up,c_dn,c_upx,rhoratx,rhomin(2),
     .   alpmin(2)
c
c
      cphmax = 0.0
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
      do j=1,25
        dz_duct(j)=0.
        jduct(1,j)=0
        jduct(2,j)=0
        jduct(3,j)=0
        jduct(4,j)=0
        jduct(5,j)=0
        zduct(j)=0.
      end do
c
      nduct=0
      nd_next=1
      cfmin=1.d+100
      c1=1.d100
      c2=2.d100
      jflu1=jflu(2,2)
      jflu2=jflu(1,2)
      jlast=1
c
      mm_up=0
      rhoratx=geo(2,3,jflu1-1)/geo(1,3,jflu1)
      if(rhoratx .gt. 1.5d0 .or. rhoratx .lt. .667d0) mm_up=1
      mm_upx=mm_up
      c_upx=geo(2,1,jflu1-1)
c
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
         if(rhoratx .gt. 1.5d0 .or. rhoratx .lt. .667d0) mm_dn=1
c
         if(c3 .eq. c4) then
c: Isospeed: c_up,c_dn are speeds to compare to:
            if((c3 .lt. c_up .or. mm_up .eq. 1 .or.
c: Case of two isospeeds together:
     .         (c3 .eq. c_up .and. (c3 .lt. c_upx .or. mm_upx .eq. 1)))
     .         .and. (c4 .lt. c_dn .or. mm_dn .eq. 1)) then
               call duct_enter(j,1,c3,0,nd_next)
            endif
            if(mm_dn .eq. 1) then
               call peak_enter(nduct,nd_next,dz_duct,zdep,j,jlast,
     .            c4,c_dn,cspan,iidiag,f_max,mm_dn)
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
               call peak_enter(nduct,nd_next,dz_duct,zdep,j,jlast,
     .            c4,c_dn,cspan,iidiag,f_max,mm_dn)
            endif
            c_upx=c_up
            mm_upx=mm_up
         elseif(c4 .lt. c3) then
c: Check if peak at top of layer:
            if(c3 .gt. c_up) then
               call peak_enter(nduct,nd_next,dz_duct,zdep,j-1,jlast,
     .            c4,c_dn,cspan,iidiag,f_max,mm_up)
            endif
c: Negative gradient with depth, check for duct at bottom of layer:
            if(c4 .lt. c_dn .or. mm_dn .eq. 1) then
               call duct_enter(j,2,c4,0,nd_next)
            endif
c: Check if peak at bottom of layer:
            if(mm_dn .eq. 1) then
               call peak_enter(nduct,nd_next,dz_duct,zdep,j,jlast,
     .            c4,c_dn,cspan,iidiag,f_max,mm_dn)
            endif
            c_upx=c_up
            mm_upx=mm_up
         endif
         mm_up=mm_dn
      enddo
      nd_next=nd_next + 1
      dz_duct(nduct)=dz_duct(nduct) + zdep(nlay-1)-zdep(jlast)
c
      if(iifb .ge. 999) then
         call duct_enter(nlay,1,geo(1,1,nlay),0,nd_next)
      endif
c
      if(nduct .eq. 0 .or. cfmin .ge. 1.d100) then
         print *,'nduct=0 or cfmin bad: BUG - tell EKW ',
     .      (csvp(j),j=1,nsvp)
         stop
      endif
c
c: Sort dz_duct:
      call hpsort_indx(nduct,dz_duct,indx_duct)
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
      do ii=1,2
         nhsp=jhsp(ii)
         ii2=3 - ii
c: Homogeneous halfspace:
         geo(ii2,1,nhsp)=geo(ii,1,nhsp)
         geo(ii2,4,nhsp)=geo(ii,4,nhsp)
         geo(ii2,2,nhsp)=geo(ii,2,nhsp)
         geo(ii2,3,nhsp)=geo(ii,3,nhsp)
         geo(ii2,5,nhsp)=geo(ii,5,nhsp)
      enddo
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
            j=jflu(ii,2) + jflu(ii,3)
            geo(1,1,j)=cpfake(ii)
            geo(2,1,j)=cpfake(ii)
            geo(1,3,j)=geo(jjmin,3,jmin)
            geo(2,3,j)=geo(1,3,j)
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
      do j=1,nlay
c: Convert p-wave attenuations:
         call alpha_conv(geo(1,4,j),geo(1,4,j),geo(1,1,j))
         call alpha_conv(geo(2,4,j),geo(2,4,j),geo(2,1,j))
c: Isospeed in cp flag:
         isp(j)=1
         if(geo(1,1,j) .ne. geo(2,1,j) .or.
     .      geo(1,4,j) .ne. geo(2,4,j)) isp(j)=0
c: Solid layer flag:
         iisol(j)=0
         iss(j)=1
      enddo
c
c: Check for isospeed layers next to no-mismatch interfaces (causes
c: problems in rp_nomm when on negative sheets):
      do j=1,nlay-1
         if(mm(j) .eq. 0 .and. h(j+1) .ne. 0.d0) then
            if(isp(j) .eq. 1 .and. isp(j+1) .eq. 1) then
               print *,'CANNOT HAVE ISOSPEED P-WAVE LAYERS ON BOTH '//
     .            'SIDES OF NO-MISMATCH INTERFACE!!'
               print *,'Layer # ',j,j+1
cpln ONLY FOR ISO SPEED LAYER
               geo(1,1,j+1)=geo(1,1,j)+0.10
               geo(2,1,j+1)=geo(2,1,j)+0.10
cpln
cpln               iibad=1
cpln               iifail=1
cpln               return
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
         rholay(1)=geo(2,3,jflu(1,2))/geo(1,3,jsol(1,1))
      endif
      if(allf(2) .eq. 0) then
         rholay(2)=geo(1,3,jflu(2,2))/geo(2,3,jsol(2,1))
      endif
c
c: Convert cphmax from grazing angle to phase velocity if necessary:
      if(cphmax .lt. 0.) then
         if(abs(cphmax) .gt. 90.) then
            print *,'Max angle cannot be >= 90 deg'
            stop
         endif
c: Make maximum angle 89.9 degrees to avoid infinite cphmax:
         cphmax=cfmin/cos(amin1(89.9e0,abs(cphmax))*acos(-1.)/180.)
      elseif(cphmax .eq. 0.) then
         cphmax=-89.9d0
      endif
c
      return
c 500   print *,'Error opening SVP file ',svp_file
c      stop
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
c: Correction by Westwood (1 instead of ii1 index) 16/10-97
      if(dabs(c) .lt. 1.d-100) return
      if(alpha1 .ge. 0.) then
         alpha2=alpha1
      else
         alpha2=-alpha1/(.001*c)
      endif
c
      return
      end
ccc
      subroutine duct_enter(jlay,ii,cmin,iifake,nd_next)
c
      implicit none
      include 'Parms_com'
      include 'i_o_svp_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      common /out_com3/ nzsr, nlay, NSCTOR,wclmn
      integer*4 nzsr, nlay, NSCTOR,wclmn
      include 'gen1_com'
c
      integer*4 jlay,ii,iifake,nd_next
      real*8 cmin
c
      if(nduct .ge. nd_next) then
         if(iidiag .ne. 0) then
            print *,'Duct skipped (two in a row): ',jlay
         endif
         return
      endif
c
c: Make sure not to put duct at no-mismatch interface to a halfspace:
      if(jlay .eq. nlay-1 .and. ii .eq. 2 .and. 
     .   mm(nlay-1) .eq. 0) then
         print *,'Duct skipped (no-mismatch lower h-space)'
         return
      elseif(jlay .eq. 2 .and. ii .eq. 1 .and. mm(1) .eq. 0) then
         print *,'Duct skipped (no-mismatch upper h-space)'
         return
      endif
      nduct=nduct + 1
      jduct(1,nduct)=jlay
      jduct(2,nduct)=ii
      jduct(3,nduct)=iifake
      jduct(4,nduct)=jflu(1,2)
      jduct(5,nduct)=jflu(2,2)
      zduct(nduct)=zdep(jlay+ii-2)
      if(iifake .gt. 0) then
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
      subroutine peak_enter(nduct,nd_next,dz_duct,zdep,jlay,jlast,
     .   c4,c_dn,cspan,iidiag,f_max,mm_dn)
c
      implicit none
      integer*4 nduct,nd_next,jlay,jlast,iidiag,mm_dn
      real*8 dz_duct(nduct),zdep(jlay),c4,c_dn,cspan(nduct),f_max,
     .   duct_width,lambda_width
c
      if(nd_next .gt. nduct) then
         if(iidiag .ne. 0) then
            print *,'Peak skipped (two in a row): ',jlay
         endif
         return
      endif
c
      duct_width=zdep(jlay) - zdep(jlast)
      lambda_width=duct_width/(min(c4,c_dn)/f_max)
      if(mm_dn .eq. 1 .or. lambda_width .gt. .2d0) then
         dz_duct(nd_next)=duct_width
         cspan(nd_next)=max(c4,c_dn)
         jlast=jlay
         nd_next=nd_next + 1
      else
c: If previous duct less than lambda/5, delete it:
         nduct=nduct-1
         if(iidiag .ne. 0) then
            print *,'Duct with width < lambda/5 deleted.'
         endif
      endif
c
      return
      end
