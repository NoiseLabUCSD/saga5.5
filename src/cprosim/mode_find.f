      subroutine mode_find(iiwrt)
c
c: Finds eigenvalues and computes mode functions.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 iiwrt,iduct,nm_duct,nm_found,iidup,ndup,ndup_max,
     .   nm_put0,nm_tot0,iidiff,ntry1
      complex*16 k_st,rr_st(3,4),k_jump
c
      ntry1=0
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
            if(jjfail .gt. 0) then
               if(iiwrite .gt. 0)
     .            write(6,*)'Return from xi_cut_set'
               return
            end if
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
cpln         write(6,*)'#1 mode_branch'
         call mode_branch(1,xkref,1,0,0,iidup,ndup,ndup_max)
            if(jjfail .gt. 0) then
               if(iiwrite .gt. 0)
     .            write(6,*)'Return from mode_branch'
               return
            end if
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
         iidiff=0
         if(iidone .eq. 0 .and. ndup .gt. 1) then
c: When mode search encountered duplicate modes before getting to end of
c: search (Im(kn) limit, e.g.), find contour starting at upper left end
c: and proceed back down and to right on it:
cpln            write(6,*)'Enter from mode_find #1',iduct
            call contour_find(k_jump,iduct,k_st,rr_st,iidiff)
cpln            write(6,*)'Exit contour_find #1',jjfail
            if(jjfail .gt. 0) return
         elseif(iiblm .ne. 0) then
c: Mode search ended by finding a branch line mode, so look for contour
c: in upper left corner of k plane:
            if(iduct .lt. nduct) then
cpln               write(6,*)'Enter contour_find #2',iifail,jjfail,
cpln     .                    iduct,nduct
               call contour_find(k_jump,iduct,k_st,rr_st,iidiff)
cpln               write(6,*)'Exit contour_find #2',iifail,jjfail,
cpln     .                    iduct,nduct
cpln               write(6,*)k_jump,k_st,rr_st,iidiff
               if(jjfail .gt. 0) return
            else
               call contour_find2(k_st,rr_st,iidiff)
               if(jjfail .gt. 0) return
            endif
c: EKW FIX 6/2/98:
            iiblm=0
         endif
         if(iidiff .eq. 1) then
c: Different contours found, so look for modes on new contour:
            nm_put0=nm_put
            nmode=nmode + 1
            call iish_xfer(iish,iish_ref,iishx,iish_refx)
c: Find modes to right of k_st:
cpln            write(6,*)'#2 mode_branch'
            call mode_branch(0,k_st,-1,0,1,1,ndup,1)
            if(jjfail.eq.1) return
            if(ndup .eq. 1) then
               if(iiwrite .gt. 0)
     .         print *,'Contour_find looped back on same branch'
               nmode=nm_tot
cpln               write(6,*)'Im going to 20'
               goto 20
            endif
            nmode=nmode + 1
            call sheet_init(k_st,0,iishx,iish_refx)
c: Find modes to left of k_st:
            iiblm=0
cpln            write(6,*)'#3 mode_branch'
            call mode_branch(0,k_st,1,0,1,1,ndup,1)
            if(jjfail.eq.1) return
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
cpln 14/4/2000 Can't find modes->continue to goto 15
            ntry1=ntry1+1
            if(ntry1 .gt. 200) then
               jjfail=1
               write(6,*)'Exhausted trying to find mode'
               return
            end if
            goto 15
         else
            nmode=nm_tot
         endif
20       continue
c
         if(iduct .eq. nduct .and. nduct .gt. 1) then
c: After finding modes at first reference depth, choose mode at which 
c: to start looking for different mode branches for other ducts:
            call k_jump_find(nmode,kn,k_jump)
         endif
c
         if(iiwrt .eq. 1) then
            print *,'Informative message: Finished checking duct '//
     .         'at depth ',zduct(kduct)
            print *,'   # modes found in duct = ',nm_duct
            write(lusvp,'(a,f7.2)') 'Informative message: Finished '//
     .         'checking duct at depth ',zduct(kduct)
            write(lusvp,'(a,i4)') '   # modes found in duct = ',nm_duct
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
      subroutine k_jump_find(nmode,kn,k_jump)
c
      implicit none
      integer*4 nmode,jm
      complex*16 kn(0:nmode),k_jump,delk
c
c: After first of several ducts, find mode at bottom of vertical section 
c: (first "evanescent" mode) to use to jump to other mode branches:
      k_jump=kn(nmode)
      do jm=nmode,2,-1
         delk=kn(jm) - kn(jm-1)
         if(-real(delk) .gt. dimag(delk)) then
            k_jump=kn(jm)
            return
         endif
      enddo
c
      return
      end
