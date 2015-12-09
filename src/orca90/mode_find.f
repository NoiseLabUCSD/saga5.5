      subroutine mode_find(iiwrt)
c
c: Finds eigenvalues and computes mode functions.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 iiwrt,iduct,nm_duct,nm_found,iidup,ndup,ndup_max,
     .   nm_put0,nm_tot0,iidiff,k_jump_done
      complex*16 k_st,rr_st(3,4),k_jump
c
      k_jump_done=0
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
c: Flag that we've found a mode to jump to in contour_find:
      k_jump_done=0
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
         call mode_branch(1,xkref,1,0,0,iidup,ndup,ndup_max)
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
         if(k_jump_done .eq. 0) goto 20
         iidiff=0
         if(iidone .eq. 0 .and. ndup .gt. 1) then
c: When mode search encountered duplicate modes before getting to end of
c: search (Im(kn) limit, e.g.), find contour starting at upper left end
c: and proceed back down and to right on it:
            call contour_find(k_jump,iduct,k_st,rr_st,iidiff)
         elseif(iiblm .ne. 0) then
c: Mode search ended by finding a branch line mode, so look for contour
c: in upper left corner of k plane:
            if(iduct .lt. nduct) then
               call contour_find(k_jump,iduct,k_st,rr_st,iidiff)
            else
               call contour_find2(k_st,rr_st,iidiff)
            endif
c: EKW FIX 6/2/98:
         iiblm=0
         endif
         if(iifail .ne. 0) goto 20
         if(iidiff .eq. 1) then
c: Different contours found, so look for modes on new contour:
            nm_put0=nm_put
            nmode=nmode + 1
            call iish_xfer(iish,iish_ref,iishx,iish_refx)
c: Find modes to right of k_st:
            call mode_branch(0,k_st,-1,0,1,1,ndup,1)
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
cc          if(ndup .eq. 1) then
cc             print *,'Contour_find looped back on same branch'
cc             nmode=nm_tot
cc             nm_put=nm_put0
cc             goto 20
cc          endif
            nm_put0=nm_put
            nmode=nmode + 1
            call sheet_init(k_st,0,iishx,iish_refx)
c: Find modes to left of k_st:
            iiblm=0
            call mode_branch(0,k_st,1,0,1,1,ndup,1)
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
            goto 15
         else
            nmode=nm_tot
         endif
20       continue
c
c: After finding modes at first reference depth, choose mode at which 
c: to start looking for different mode branches for other ducts:
         call k_jump_find(nmode,kn,k_jump,k_jump_done)
c
         if(iiwrt .eq. 1) then
co            print *,'Informative message: Finished checking duct '//
co     .         'at depth ',zduct(kduct)
co            print *,'   # modes found in duct = ',nm_duct
co            write(2,'(a,f7.2)') 'Informative message: Finished '//
co     .         'checking duct at depth ',zduct(kduct)
co            write(2,'(a,i4)') '   # modes found in duct = ',nm_duct
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
      subroutine k_jump_find(nmode,kn,k_jump,k_jump_done)
c
      implicit none
      integer*4 nmode,k_jump_done,jm
      complex*16 kn(0:nmode),k_jump,k_jumpx,delk
c
c: After first of several ducts, find mode at bottom of vertical section 
c: (first "evanescent" mode) to use to jump to other mode branches:
      if(nmode .gt. 1) then
         k_jumpx=kn(nmode)
         do jm=nmode,2,-1
            delk=kn(jm) - kn(jm-1)
            if(-real(delk) .gt. dimag(delk)) then
               k_jumpx=kn(jm)
               goto 10
            endif
         enddo
10       k_jump=k_jumpx
         k_jump_done=1
      endif
c
      return
      end
