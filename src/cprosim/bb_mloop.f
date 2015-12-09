      subroutine bb_mloop(jm0,nmode_bb,nmbbtot,iisg,f_st,jf_st,
     .   ndied)
c
c: Loops over nmode_bb modes, whose indices for kn,phi,eig_char lie in 
c: kn_indx, and calls routine for following them in frequency.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm0,nmode_bb,nmbbtot,iisg,jf_st,jmx,jmk,
     .   jmo,jfbb,nfinc,jfbb2,nm1,nm2,nm3,jj0,
     .   nsave0,ii,ncmode,nm_fmax,iimrg,jm_mrg,nfout,kductp,iichng,
     .   kdone,nf_done,ndied,nf_blm,iixi_cut
      real*8 dfinc,dfincx
      real*4 f_st
      complex*16 k,r1r2(3,4),Lp,dL_dkp,dL_dwp,dW_dk,kp,vg,k_mrg
c
      nm_fmax=kn_indx(nmode_bb)
      nsave=0
      nskip=0
      ndied=0
      do jmx=1,nmode_bb
         jmk=kn_indx(jmx)
         jmo=jm0 + jmx
         nhigh=0
         iidone=0
         nrise=0
         iilk=0
         iixi_cut=0
         kim_max=0
         nsave0=nsave
c: Initialize frequency to that where modes were found:
         f_hz=f_st
c: Re-initialize frequency so that roundoff errors don't cause trouble:
         call freq_init
c: Initialize sheets and set xkhratp to kn(jmk)/kw (at f_st):
cc       call sheet_init(kn(jmk),1,iish,iish_ref)
c: Initialize sheet variables as found for mode jmk at f_st:
         xkhrat=kn(jmk)/kw0
         call iish_code(iish,iish_ref,iishn(jmk),-1)
c: Compute field at all receivers for this mode at f_st:
         jfbb=jf_st
         call bb_field(kn(jmk),phi,dphi,dpsi,exp_gbs,jmk,tf,jfbb,jmo)
         nm1=jmk
         nm2=nm_fmax + 2
         nm3=nm_fmax + 3
         ii=-1
         call bb_enter(kn(jmk),eig_char(1,jmk),eig_char(4,jmk),
     .      eig_char(5,jmk),knbb,eig_bb,jfbb,jmo,nfbb,nm_fmax,
     .      iish,iish_ref,iish_bb,nmbb)
         if(iift .eq. 1) call bb_write(jmo,jfbb,faxbb,kn(jmk),kw0,iish)
c
         dfinc=1.d0
         kp=kn(jmk)
         Lp=eig_char(1,jmk)
         dL_dkp=eig_char(2,jmk)
         dL_dwp=eig_char(3,jmk)
c: Check for change in reference depth:
         if(nzref(jmk) .ne. kduct) then
            kduct=nzref(jmk)
            call zref_chng
            if(iidiag .ge. 2) print *,'changed reference '//
     .         'depth for mode ',jmo
         endif
         ncmode=nctot
c: For each mode found at f_st, track it as a function of frequency:
         kdone=0
         nf_done=0
         nf_blm=0
         do while(kdone .eq. 0)
15          call bb_fmarch(jmk,jmo,k,r1r2,kp,Lp,dL_dwp,dL_dkp,
     .         iisg,dfinc,dfincx,nm2,jfbb2,iixi_cut,kdone)
            if(kdone .eq. 1) goto 50
            nfinc=max(1,nint(dfinc))
            jfbb2=jfbb + iisg*nfinc
            phim_bb(jfbb2)=.5d0*phim_bb(jfbb)
            if(iifail .eq. 1) then
               print *,'Mode ',jmk,' could not be tracked past f = ',
     .            f_hz
               iifail=0
               goto 50
            endif
c
c: Check for mode merging:
            call bb_merge(k,knbb,nfbb,nm_fmax,jfbb2,jmo,errdk2,iimrg,
     .         jm_mrg,k_mrg,nmerge,faxbb,iidiag)
            if(iimrg .eq. 1 .or. iimrg .eq. 2) then
c: If modes merged, first try going back to prev freq with current mode:
c: iimrg=1 ==> dfinc=.125, iimrg=2 ==> dfinc=.05:
               dfinc=.125/((iimrg-1)*3 + 1)
               f_hz=faxbb(jfbb)
               goto 15
            elseif(iimrg .eq. 3 .or. iimrg .eq. 4) then
c: If modes merged, next try going back to prev freq with duplicate mode:
               f_hz=faxbb(jfbb)
               call freq_chng
               jj0=(jm_mrg - 1)*nfbb + jfbb
               call iish_code(iish,iish_ref,iish_bb(jj0),-1)
               if(nzref(jmk) .ne. nzref(jm_mrg)) then
                  kduct=nzref(jm_mrg)
                  call zref_chng
                  if(iidiag .ge. 2) then
                     print *,'Changed ref depth for dup mode ',jm_mrg
                  endif
               endif
               kp=k_mrg
cc             call sheet_init(kp,0,iish,iish_ref)
               xkhrat=kp/kw0
               call r1r2_calc(kp,r1r2,3,1,jjfail)
               if(jjfail .gt. 0) return
               Lp=r1r2(1,4)
               dL_dkp=r1r2(2,4)
               dL_dwp=r1r2(3,4)
c: iimrg=3 ==> dfinc=.125, iimrg=4 ==> dfinc=.05:
               dfinc=.125/((iimrg-3)*3 + 1)
               goto 15
            elseif(iimrg .gt. 4) then
               print *,'Modes found to merge. Unable to resolve, '//
     .            'skipping: ',jmo,jm_mrg,faxbb(jfbb)
               nskip=nskip + 1
               goto 50
            endif
c
c: Keep track of minimum Im(k) in order to stop mode search due to rmin:
            kim_bb(jfbb2)=dmin1(max(0.d0,dimag(k)),kim_bb(jfbb2))
c: Flag for leaky modes:
            if(dreal(k) .lt. kcrmin .or. nsvmin .eq. nlay) iilk=1
c: Check if k to left of kremin (related to cphmax):
            if(dreal(k) .lt. kremin) then
               if(iidiag .ge. 1) print *,'BB cutoff due to kremin: ',
     .            jmo,k/kw0,f_hz
               goto 50
            endif
c: Check if Im(k) so high that mode is weak at shortest range of interest:
            if(iilk .eq. 1 .or. nrise .gt. 5) then
               kim_max=kim_bb(jfbb2) + dkim
            else
c: Safety factor for modes not leaky yet:
               kim_max=kim_bb(jfbb2) + 5.*dkim
            endif
            if(dimag(k) .gt. kim_max) then
               nhigh=nhigh + 1
c: Check if 2 modes in a row have been high & traj heading higher (see p.150):
c: Ignore this criterion since island modes can have any k-derivative:
ccx            if(nhigh .ge. 2 .and. dreal(r1r2(2,4)) .gt. 0.d0) then
               if(nhigh .ge. 2) then
                  if(iidiag .ge. 1) print *,'BB cutoff due to nhigh: ',
     .               jmo,jfbb,k/kw0,f_hz
                  goto 50
               endif
            else
               nhigh=0
            endif
c
            if(dimag(k) .gt. dimag(kn(nmode-1))) then
               nrise=nrise + 1
            else
               nrise=0
            endif
c
c: Numerical derivative test for dL/dw:
ctemp f_hz=f_hz + .01
ctemp call freq_chng
ctemp call r1r2_calc(k,r1r2x,3,1)
ctemp print *,'An,Nu = ',r1r2(3,4),(r1r2x(1,4)-r1r2(1,4))/(.01*2.*pie)
ctemp f_hz=f_hz - .01
ctemp call freq_chng
c
c: Valid eigenvalue found at new frequency:
            dW_dk=-2.*r1r2(1,1)*gamiref*r1r2(2,3)
c: Place mode functions in spare location at nm_fmax+1:
            call mode_fun(k,r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
     .         psi,dpsi,exp_gbs,nm2)
c: Enter k into knbb, etc:
            vg=-r1r2(2,4)/r1r2(3,4)
            call bb_enter(k,r1r2(1,4),vg,r1r2(1,1),knbb,eig_bb,
     .         jfbb2,jmo,nfbb,nm_fmax,iish,iish_ref,iish_bb,nmbb)
c: Interpolate mode functions and compute field for interpolated freqs:
            if(nfinc .gt. 1) call bb_interp(phi,dphi,psi,dpsi,exp_gbs,
     .         phix,dphix,psix,dpsix,expx_gbs,nzsr,nm1,nm2,nm3,iisg,
     .         nfinc,knbb,tf,nfbb,jfbb,nm_fmax,jmo,faxbb,cfmin,iish,
     .         eig_bb,iift,kim_bb,phim_bb,iish_bb)
c
c: Compute field at all receivers at f_hz:
            jfbb=jfbb2
            call bb_field(k,phi,dphi,dpsi,exp_gbs,nm2,tf,jfbb,jmo)
            if(iift .eq. 1) call bb_write(jmo,jfbb,faxbb,k,kw0,iish)
c
c: Check if best duct in which to find mode has changed:
cc          if(nduct .gt. 1) then
cc             kductp=kduct
cc             call zduct_chng(phi,dphi,jmk,nm2,iichng,k)
cc             if(iichng .eq. 1) then
cc                call zref_chng
cc                call r1r2_calc(k,r1r2,3,1)
cc                if(iidiag .ge. 1) 
cc   .               print *,'changed duct: ',kductp,kduct,jmo,f_hz
cc             endif
cc          endif
c
            nm1=nm2
            nm2=nm2 + ii
            ii=-ii
cxx   print *,'nm1,nm2,ii = ',nm1,nm2,ii
c: Update "previous values" for next jump in frequency:
            dfinc=dfincx
            kp=k
            Lp=r1r2(1,4)
            dL_dkp=r1r2(2,4)
            dL_dwp=r1r2(3,4)
            if(dfinc .gt. 1.d0) dfinc=nint(dfinc)
            if(iisg .lt. 0) then
               if(jfbb .le. 1) kdone=1
            else
               if(jfbb .ge. nfbb) kdone=1
            endif
            nf_done=nf_done + 1
         enddo
50       continue
         ncmode=nctot - ncmode
         print *,'jm,nsave,fcut,nm = ',jmo,nsave-nsave0,faxbb(jfbb),
     .      nf_done,float(ncmode)/max(1,nf_done)
c         pause
         nmbbtot=nmbbtot + nf_done
c
c: Output mode characteristics if desired:
         if(iiout .ne. 0) then
            call bb_out(nfbb,nm_fmax,nzsr,jmo,knbb,phibb,dpsibb,nh_off)
         endif
         if(iimf .eq. 1) then
            nfout=nfbb-jfbb+1
            call bb_mfout(jmo,jfbb,nfout,phibb,r4mat1,r4mat2)
         elseif(iimf .eq. 2) then
cpln            call bb_mfout2(jmo,jfbb,nfout,nm_fmax,phibb,r4mat1,r4mat2)
         endif
      enddo
c
      return
      end
ccc
      subroutine bb_blmodes
c
c: Find the branch line modes from fmax to fmin in a brute force manner.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      include 'lab_com'
      integer*4 jm,jfbb,nmbbtot,nm_fmax,ndup,ndup_max,jmo
c
      kduct=indx_duct(1)
      call zref_chng
      if(nsvmin .ne. nlay) then
         print *,'Bad nsvmin from bb_blm'
         return
      endif
      print *,'Enter nblm_max: '
      read(5,*) nblm_max

      nmbbtot=0
      do jfbb=nfbb,1,-1
         f_hz=faxbb(jfbb)
         call freq_chng
         nmode=0
         nm_put=0
c
         nblm=0
         write(6,*)'#4 mode_branch'
         call mode_branch(1,xkref,1,0,0,0,ndup,ndup_max)
         jmo=nmbb(jfbb)
c
         nm_fmax=max(nm_fmax,nmbb(jfbb) + nm_put)
         do jm=1,nm_put
            jmo=jmo + 1
            call bb_field(kn(jm),phi,dphi,dpsi,exp_gbs,jm,tf,jfbb,jm)
            call bb_enter(kn(jm),eig_char(1,jm),eig_char(4,jm),
     .         eig_char(5,jm),knbb,eig_bb,jfbb,jmo,nfbb,nm_put,
     .         iish,iish_ref,iish_bb,nmbb)
            if(iift .eq. 1) call bb_write(jmo,jfbb,faxbb,kn(jm),
     .         kw0,iish)
         enddo
         nmbbtot=nmbbtot + nm_put
         print *,'Done f = ',sngl(f_hz),nm_put,nmbb(jfbb),
     .      sngl(dimag(kn(nm_put))*8685.9)
      enddo
c
      print *,'nm_fmax = ',nm_fmax
      write(lusvp,120) nmbbtot
120   format('BRUTE FORCE BLM TOTAL # MODES = ',i8)
c
      return
      end
