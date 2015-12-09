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
