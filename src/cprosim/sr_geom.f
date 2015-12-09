      subroutine sr_geom(rng_srx,nsrcx,nrecx)
c
c: Computes horizontal ranges between nsrc sources at (xsrc,ysrc,zsrc)
c: and nrec receivers at (xrec,yrec,zrec) and places them in the
c: array rng_sr(1:nsrc,1:nrec).  Then places sorted ranges in 
c: range(1:nrng), where nrng=nsrc*nrec.
c  
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 nsrcx,nrecx,jsrc,jrec,jrng,j,jj,kk,jr,
     .   iibad,kk0,nzsr0,jd
      real*4 rng_srx(nsrcx,nrecx),rngcur,zp1,zp2
      real*8 pierad,rlast,delz
      data pierad/0.01745329251994/
c
      if(iimf .eq. 0 .and. iisig .eq. 0) nzmf=0
      call uni_space(nzmf,zmf,1.e0)
      call uni_space(nzs,zsrc,1.e0)
      call uni_space(nrec,zrec,1.e0)
      if(iigbs .ne. 0) then
         call uni_space(nth_gbs,th_gbs,1.e0)
         call uni_space(nb_gbs,b_gbs,1.e0)
      endif

      if(iigeom .eq. 1) call uni_space(nsrc,rkm,1.e0)
c
c: Find interfaces across which to compute terms for mode orthogonality:
      n_int=0
      if(iidiag .eq. -2 .and. (iimf .ne. 0 .or. iisig .ne. 0)) then
         do j=1,nlay-1
            jj=j+1
            if(zdep(j) .ge. zmf(1) .and. zdep(j) .le. zmf(nzmf) .and.
     .         mm(j) .eq. 1 .and. (iisol(j) .eq. 1 .or.
     .         iisol(jj) .eq. 1)) then
               n_int=n_int + 1
               j_int(n_int)=j
            endif
         enddo
         if(iicw .eq. 1) then
            delz=.01*cfmin/fcw(1)
         else
            delz=.01*cfmin/fmin
         endif
      endif
c
      nzsr=nzs + nrec + nzmf + 2*n_int
c: Don't include duct depths in zsr if doing real axis version:
      if(nduct .gt. 1 .and. iirx .le. 0) nzsr=nzsr + nduct
      nm_max2=NSR_NM_MAX/max(1,nzsr)
      iibad=0
      call mem_lim(nzsr,NSRMAX,MLINE,LML,'nzsr',4,'NSRMAX',6,iibad,0)
      call mem_lim(nsrc*nrec,NSNRMAX,MLINE,LML,'nsrc*nrec',9,
     .   'NSNRMAX',7,iibad,0)
      if(iibad .eq. 1) then
         write(*,*)'nzsr,NSRMAX,nsrc,nrec,NSNRMAX'
         write(*,*)nzsr,NSRMAX,nsrc,nrec,NSNRMAX
         stop 'stopped from sr_geom'
      endif
c
      if(iigeom .eq. 1) then
         do j=1,nsrc
            xsrc(j)=rkm(j)
cpln            xsrc(j)=rkm(j)*1000.e0
            ysrc(j)=0.
cpln            t_src(j)=rkm(j)
            t_src(j)=rkm(j)/1000.e0
            range(j)=rkm(j)
         enddo
c: Compute receiver depths:
         do j=1,nrec
            xrec(j)=0.
            yrec(j)=0.
            rec_lab(j)=zrec(j)
         enddo
      endif
c
c: Check which ducts have sources and/or receivers in them:
      do jd=1,nduct
         jd_ch(jd)=0
         if(jduct(1,jd) .eq. nlay .or. jduct(1,jd) .eq. 1) goto 15
         zp1=zpeak(jd)
         zp2=zpeak(jd+1)
         do jrec=1,nrec
            if(zrec(jrec) .ge. zp1 .and. zrec(jrec) .le. zp2) then
               jd_ch(jd)=1
               goto 15
            endif
         enddo
         do j=1,nzs
            if(zsrc(j) .ge. zp1 .and. zsrc(j) .le. zp2) then
               jd_ch(jd)=1
               goto 15
            endif
         enddo
         do j=1,nzmf
            if(zmf(j) .ge. zp1 .and. zmf(j) .le. zp2) then
               jd_ch(jd)=1
               goto 15
            endif
         enddo
15       continue
      enddo
c
c: Create an array zsr of src depths zsrc, rec depths zrec, and mode function
c: depths zmf:
      do jrec=1,nrec
         zsr(jrec)=zrec(jrec)
         zsr_im_gbs(nrec+j)=0.d0
         zsrmin=min(zsrmin,zrec(jrec))
         zsrmax=max(zsrmax,zrec(jrec))
      enddo
      do j=1,nzs
         zsr(nrec+j)=zsrc(j)
         zsr_im_gbs(nrec+j)=0.d0
         zsrmin=min(zsrmin,zsrc(j))
         zsrmax=max(zsrmax,zsrc(j))
      enddo
      if(iigbs .ne. 0) then
         do j=1,nzs
            zsr_im_gbs(nrec+j)=-b_gbs(j)*sin(th_gbs(j)*pierad)
         enddo
      endif
      nzsr0=nrec+nzs
      do j=1,nzmf
         zsr(nzsr0+j)=zmf(j)
         zsr_im_gbs(nzsr0+j)=0.d0
      enddo
      nzsr0=nzsr0+nzmf
      if(nduct .gt. 1 .and. iirx .le. 0) then
         do j=1,nduct
            zsr(nzsr0+j)=zduct(j)
            zsr_im_gbs(nzsr0+j)=0.d0
         enddo
         nzsr0=nzsr0 + nduct
      endif
      do j=1,n_int
         zsr(nzsr0+2*j-1)=zdep(j_int(j))
         zsr_im_gbs(nzsr0+2*j-1)=0.d0
         zsr(nzsr0+2*j)=zdep(j_int(j)) + delz
         zsr_im_gbs(nzsr0+2*j)=0.d0
      enddo
c: Sort real depths in zsr, keeping imaginary depths in zsr_im_gbs:
      nzsr0=nzsr+1
      call hpsort_re_im(nzsr,zsr,zsr_im_gbs,zsr_indx(nzsr0))
c: Sort resulting indices so that zsr_indx(i)=original index in zsr of 
c: current i'th element:
      call hpsort_i4_indx(nzsr,zsr_indx(nzsr0),zsr_indx)
c
c: Discard duplicates:
      if(nzsr .gt. 1) then
         j=1
         nzsr0=nzsr
10       if(zsr(j) .eq. zsr(j+1) .and. 
     .      zsr_im_gbs(j) .eq. zsr_im_gbs(j+1)) then
            nzsr=nzsr-1
            do jj=j+1,nzsr
               zsr(jj)=zsr(jj+1)
               zsr_im_gbs(jj)=zsr_im_gbs(jj+1)
            enddo
c: Keep indices correct by decrementing those above j:
            do jj=1,nzsr0
               if(zsr_indx(jj) .gt. j) zsr_indx(jj)=zsr_indx(jj)-1
            enddo
            j=j-1
         endif
         j=j+1
         if(j .lt. nzsr) goto 10
      endif
c
c: Copy zsr into zsrgeom if only geometry
c: (reuse eigenvalues and mode functions)
c
      if((i_geom.eq.0).and.(i_call_porter.eq.1)) then
         do jj=1,nzsr
            zsrgeom(jj)=zsr(jj)
         end do
         nzsrgeom=nzsr
      end if
c
c: Find indices of zrec,zsrc,zmf in zsr:
      do jrec=1,nrec
         mzrec(jrec)=zsr_indx(jrec)
      enddo
      do j=nrec+1,nrec+nzs
         mzsrc(j-nrec)=zsr_indx(j)
      enddo
      nzsr0=nrec+nzs
      do j=nzsr0+1,nzsr0+nzmf
         mzmf(j-nzsr0)=zsr_indx(j)
      enddo
      nzsr0=nzsr0+nzmf
      if(nduct .gt. 1 .and. iirx .le. 0) then
         do j=nzsr0+1,nzsr0+nduct
            mzduct(j-nzsr0)=zsr_indx(j)
         enddo
         nzsr0=nzsr0 + nduct
      endif
      do j=nzsr0+1,nzsr0+2*n_int
         mzint(j-nzsr0)=zsr_indx(j)
      enddo
c
      if(iigeom .eq. 1) then
         nrng=nsrc
         call mem_lim(nrng,NRNGMAX,MLINE,LML,'nrng',4,'NRNGMAX',7,
     .      iibad,1)
         do jrng=1,nrng
            range(jrng)=xsrc(jrng)
         enddo
c
         call hpsort_indx(nrng,range,zsr_indx)
c
         do jrng=1,nrng
            rlast=range(jrng)
            jsrc=zsr_indx(jrng)
            nrec_jr(jrng)=nrec
            kk0=(jrng-1)*nrec
            krec_jr(jrng)=kk0
            do jrec=1,nrec
               jrec_jr(1,kk0+jrec)=jsrc
               jrec_jr(2,kk0+jrec)=jrec
               rng_srx(jsrc,jrec)=rlast
            enddo
         enddo
c
      else
c: Compute ranges between source positions and receivers:
         nrng=0
         do jsrc=1,nsrc
            do jrec=1,nrec
               rngcur=sqrt((xsrc(jsrc)-xrec(jrec))**2 + 
     .            (ysrc(jsrc)-yrec(jrec))**2)
               rng_srx(jsrc,jrec)=rngcur
               do jr=nrng,1,-1
c: Skip duplicate ranges, but keep track of # s/r pairs had that range:
                  if(rngcur .eq. range(jr)) then
                     nrec_jr(jr)=nrec_jr(jr) + 1
                     goto 45
                  endif
               enddo
               nrng=nrng + 1
               call mem_lim(nrng,NRNGMAX,MLINE,LML,'nrng',4,'NRNGMAX',
     .            7,iibad,1)
               range(nrng)=rngcur
               nrec_jr(nrng)=1
45          enddo
         enddo
c: Sort ranges:
         call hpsort(nrng,range)
c: Set krec_jr(1:nrng), the starting index in jrec_jr(1:2,1:nsrc*nrec):
         krec_jr(1)=0
         do jrng=2,nrng
            krec_jr(jrng)=krec_jr(jrng-1) + nrec_jr(jrng-1)
         enddo
c: Reset nrec_jr to 0 for use in setting jrec_jr next:
         do jrng=1,nrng
            nrec_jr(jrng)=0
         enddo
         do jsrc=1,nsrc
            do jrec=1,nrec
c: Find index of rng_srx(jsrc,jrec) in range(1:nrng):
               rlast=rng_srx(jsrc,jrec)
               jrng=1
               call hunt(range,nrng,rlast,jrng)
               if(range(jrng) .ne. rlast) jrng=jrng + 1
c: Increment # s/r pairs that have this range:
               nrec_jr(jrng)=nrec_jr(jrng) + 1
               kk=krec_jr(jrng) + nrec_jr(jrng)
c: Keep track of source and rec index that had this range:
               jrec_jr(1,kk)=jsrc
               jrec_jr(2,kk)=jrec
            enddo
         enddo
      endif
c
      iibad=0
      if(rmax .eq. 0.e0 .and. nrng .gt. 0) rmax=range(nrng)/1000.e0
      if(iirx .le. 0 .and. rmax .le. 0.e0) then
         print *,'rmax = 0 in sr_geom!  Set rmax in _opt file.'
         iibad=1
      endif
      if(rmin .eq. 0.e0) then
         if(nrng .gt. 0) then
            rmin=min(998.d0,range(1)/1000.d0)
cpln            write(6,*)'rmin: ',rmin
        else
            print *,'rmin=0, but no S/R geometry given to compute'//
     .         ' rmin automaticallly.'
            iibad=1
         endif
      endif
c ONLY FOR ITWK
cpln      rmin=0.5
      if(rmin .lt. 0.e0) then
         nm_lim=nint(-rmin)
         rmin=999.
      else
         nm_lim=NM_MAX
      endif
      if(nsrc .gt. 0 .and. range(1) .le. 0.d0) then
         print *,'All ranges must be > 0!',(range(j),j=1,nrng)
         iibad=1
      endif
      if(iibad .eq. 1) stop
c: Compute maximum error in eigenvalues so that phase at maximum range
c: is within 1 degree of exact:
c: Be sure rmax is at least 50 water depths here:
cpln      write(6,*)
cpln          write(6,*)'Before RMAX: ',rmax
cpln      write(6,*)
      rmax=amax1(rmax,.05*sngl(zdep(nlay-1)))
c
cpln      write(6,*)
cpln        write(6,*)'rmin,rmax,zdep: ',rmin,rmax,0.05*sngl(zdep(nlay-1))
cpln      write(6,*)
cpln      do j=1,nsrc
cpln      write(6,*)'rkm =',(rkm(j),j=1,nsrc)
cpln      write(6,*)'range: ',(range(j),j=1,nsrc)
cpln      end do
cpln      write(6,*)
cpln      write(6,*)'Vertical tilt: ',dtiltvp(1),dtiltvp(nrec)
cpln      write(6,*)
cpln      write(6,*)'Horizontal tilt: ',dtilthp
cpln      write(6,*)
cpln      write(6,*)'fmin,fmax,df: ',fmin,fmax,df
cpln      write(6,*)faxbb(1),faxbb(nfbb)
cpln  write(6,*)
c      do jsrc=1,nzs
cpln       write(6,*)'sd= ',(zsrc(jsrc),jsrc=1,nzs)
c      end do
cpln  write(6,*)
cpln      do jrec=1,nrec
cpln      write(6,*)'rd= ',(zrec(jrec),jrec=1,nrec)
c      end do
c      write(6,*)
c      pause
c
c: Compute maximum error in eigenvalues so that phase at maximum range
c: is within 1 degree of exact:
c: Be sure rmax is at least 50 water depths here:
      rmax=amax1(rmax,.05*sngl(zdep(nlay-1)))
      errdkms=(twpie/(360.*1000.*rmax))**2
      errdk2=4.d0*errdkms
      if(phfac .lt. 2.e0) phfac=2.e0
      if(db_cut .eq. 0.e0) then
         db_cut=50.e0
      elseif(db_cut .lt. 30.e0) then
         print *,'Warning: db_cut low (>=30 recommended)...',db_cut
      endif
c
c: kim_fac=-ln[10**(-dB_down/20)] is used to find max range mode is
c: significant for (used in mode_field) see p. 150,131:
ccc   kim_fac=-dlog(10.**(-db_cut/20.d0))
ccc   dkim=kim_fac/rmin
c: Set max IM[kn] to allow for min range of interest rmin(km) (see pp.131,150):
      if(rmin .lt. 999.) then
         dkim=db_cut/(8685.9*rmin)
      else
         dkim=1.d100
      endif
c: kim_fac will be used with range in m in mode_field:
      kim_fac=db_cut/8.6859
c
      if(isp(nlay) .eq. 0 .and. allf(1) .eq. 1 .and. rmin .ge. 999.
     .   .and. cphmax .gt. geo(1,1,nlay)) then
         write(6,120) ' '
         write(6,120) 'WARNING: # branch line modes not '//
     .      'limited by cphmax: ',cphmax,sngl(geo(1,1,nlay))
         write(6,120) '   Use of RMIN instead highly recommended.'//
     .      '  Control-c to terminate now.'
         write(6,120) ' '
         write(lusvp,120) 'WARNING: # branch line modes not '//
     .      'limited by cphmax: ',cphmax,sngl(geo(1,1,nlay))
         write(lusvp,120) '   Use of RMIN instead highly recommended.'
120      format(a,f9.2,2x,f9.2)
      endif
c
      return
      end
