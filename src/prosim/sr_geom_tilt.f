C:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************
      subroutine sr_geom(isum)
c
c: Computes horizontal ranges between nsrc sources at (xsrc,ysrc,zsrc)
c: and nrec receivers at (xrec,yrec,zrec) and places them in the
c: array rng_sr(1:nsrc,1:nrec). It places no more sorted ranges in 
c: range(1:nrng), where nrng=nsrc*nrec.
c  
      implicit none
      include 'Parms_com'
      include 'i_o_opt_com'
      include 'i_o_1b_com'
      include 'i_o_2_com'
      include 'gen1_com'
      common /svp_gaus/ theta1, theta2, f_max
      real*4 theta1,theta2
      real*8 f_max
  
      integer*4 nsrcx,jrec,jrng,j,jsr,jzs,jj,
     .   iibad,nzsr0,iisort,isum,idx,jdx,displcm,dplcm
      integer*4 z_indx(NVRMAX)
      real*4 zsrx(NVRMAX)
      real*4 c0,twpie
      real*4 zapp(NVRMAX)
      real*8 rlast
      complex*16 cfac
c
      data twpie/6.28318530717959/,c0/1500.0/
c: cfac = sqrt(2*pi)*exp(i*pi/4) (see p. 7, ORCA II):
      data cfac/(1.77245385090552d0,1.77245385090552d0)/
c
      nsrcx = iabs(nsrc)
c
c: Find interfaces across which to compute terms for mode orthogonality:
c
      idx = nrec/nrecusr
c
      nzsr=idx*(nzs + nrecusr)
c: Don't include duct depths in zsr if doing real axis version:
      if(nduct .gt. 1 .and. iirx .eq. 0) nzsr=nzsr + nduct
      call mem_lim(nzsr,NVRMAX,MLINE,LML,'nzsr',4,'NVRMAX',6,iibad,0)
c
c: Create an array zsr of src depths zsrc, rec depths zrec, and mode function
c: depths zmf:
      displcm = 0
      dplcm = 0
      nzsr0=nrecusr+nzs
      if(nduct .gt. 1 .and. iirx .eq. 0) nzsr0=nzsr0 + nduct
      do jdx=1,idx
       do jrec=1,nrecusr
          zsrx(jrec+dplcm)=zrec(jrec+displcm)
       enddo
       do j=1,nzs
          zsrx(dplcm+nrecusr+j)=zsrc(j)
       enddo
       if(nduct .gt. 1 .and. iirx .eq. 0) then
          do j=1,nduct
             zsrx(dplcm+nzsr0+j)=zduct(j)
          enddo
c          nzsr0=nzsr0 + nduct
       endif

       displcm = displcm + nrecusr
       dplcm = dplcm + nzsr0

      enddo


c       call hpsort(nzsr0,zsrx)


      displcm = 0
      dplcm = 0
      do jdx=1,idx
c
       do jrec=1,nzsr0
          zapp(jrec)=zsrx(jrec+dplcm)
       enddo
       call r_hpsort_indx(nzsr0,zapp,z_indx)
c
c: Find indices of zrec,zsrc,zmf in zsr:
cc fbv(!)      if(isum .eq. 2) then
       do jzs=1,nzs
c: mzsrc(jzs) (jzs=1:nzs) gives the index of zsrc(jzs) in zsr:
         jsr=1
c fbv         do while (zsrx(jsr+dplcm) .ne. zsrc(jzs))
         do while (zapp(jsr) .ne. zsrc(jzs))
           jsr=jsr+1
         enddo
cc fbv         call hunt(zsrx,nzsr,dble(zsrc(jzs)),jsr)
cc fbv         if(zsrx(jsr) .ne. zsrc(jzs)) jsr=jsr + 1
cc!         mzsrc(jzs+displcm)=jsr+dplcm
         mzsrc(jzs+dplcm)=jsr+dplcm
       enddo
cc!cc fbv(!)      endif

cc! mzrec is not used any more so do not waste time!
cc!c fbv       do jrec=1,nrecusr Probably nzsr0 instead of nrecusr
cc!       do jrec=1,nzsr0
cc!         jsr=1
cc!c: mzrec(jrec) (jrec=1:nrec) gives the index of zrec(jrec) in zsr:
cc!c fbv         do while (zsrx(jsr+dplcm) .ne. zrec(jrec+displcm)) !Pbly
cc!         do while (zapp(jsr) .ne. zrec(jrec+displcm))
cc!           jsr=jsr+1
cc!         enddo
cc!cc fbv         call hunt(zsrx,nzsr,dble(zrec(jrec+displcm)),jsr)
cc!cc fbv         if(zsrx(jsr) .ne. zrec(jrec+displcm)) jsr=jsr + 1
cc!c fbv         mzrec(jrec+displcm)=jsr+dplcm
cc!         mzrec(jrec+dplcm)=jsr+dplcm
cc!       enddo
ccc
       if(nduct .gt. 1) then
         do j=1,nduct
            jsr=1
            call hunt(zsrx,nzsr,sngl(zduct(j)),jsr)
c fbvv            if(zsrx(jsr) .ne. zduct(j)) jsr=jsr + 1
c fbv            mzduct(j)=jsr
         enddo
       endif

c fbv       do jrec=1,nzsr0
c fbv         zapp(jrec)=zsrx(jrec+dplcm)
c fbv       enddo
c fbv       call hpsort_indx(nzsr0,zapp,z_indx)

       do jrec=1,nzsr0
         zsr(jrec+dplcm)=zapp(jrec)
         zsr_indx(jrec+dplcm)=z_indx(jrec)
       enddo
       displcm = displcm + nrecusr
       dplcm = dplcm + nzsr0
       jrec=0
c       nzsr0=0
      enddo
c
c fbv      call hpsort_indx(nzsr0,zsr,zsr_indx)
c
      if(isum .eq. 2) then
       iisort=1
       if(iigeom .eq. 1) then
c         nrng=nsrc
         nrng=nsrcx
         call mem_lim(nrng,RGMAX,MLINE,LML,'nrng',4,'RGMAX',
     .      7,iibad,1)
         do jrng=1,nrng
            range(jrng)=rkm(jrng)
         enddo
c
c: Check if ranges in increasing order:
         rlast=-1.d0
         jrng=1
         do while(range(jrng) .ge. rlast)
            rlast=range(jrng)
            jrng=jrng + 1
            if(jrng .gt. nsrc) goto 15
         enddo
15       continue
       endif
c
       do j=1,nrng
          do jj=1,nrecusr
             sq2pir(j,jj)=cfac/dsqrt(range(j)+dtiltp(jj))
          end do
       enddo
      endif
c
      if(rmax .eq. 0.e0 .and. nrng .gt. 0) rmax=range(nrng)/1000.e0
c      if(rmin .eq. 0.e0 .and. nrng .gt. 0) rmin=range(1)/1000.e0
      if(rmin .gt. range(1)/1000.e0) rmin=range(1)/1000.e0
      if(rmin .eq. 0.e0) then
         if(cphmax .lt. 1.e10) then
            rmin=1.e-20
         else
            print *,'rmin = 0 in sr_geom!  Set rmin in _opt file.'
            stop
         endif
      endif
c      if(rmin .lt. 0.e0) then
c         nm_lim=nint(-rmin)
c         rmin=1.e-20
c      else
         nm_lim=NM_MAX
c      endif
      if(nsrc .gt. 0 .and. range(1) .le. 0.d0) then
         print *,'All ranges must be > 0!',(range(j),j=1,nrng)
         stop
      endif
c: Compute maximum error in eigenvalues so that phase at maximum range
c: is within 1 degree of exact:
c fbv      errdkms=(twpie/(360.*1000.*rmax))**2
c fbv      errdk2=4.d0*errdkms
c: Be sure rmax is at least 50 water depths here:
      rmax=dmax1(dble(rmax),.05d0*zdep(nlay-1))
      if(phfac .lt. 4.e0) phfac=4.e0
      if(db_cut .eq. 0.e0) then
         db_cut=50.e0
      elseif(db_cut .lt. 30.e0) then
         print *,'Warning: db_cut low (>=30 recommended)...',db_cut
      endif
c
c: Initialize minimum IM[kn]:
      kim_min=1.d100
      kim_max=1.d100
c: kim_fac=-ln[10**(-dB_down/20)] is used to find max range mode is
c: significant for (used in mode_field) see p. 150,131:
ccc   kim_fac=-dlog(10.**(-db_cut/20.d0))
ccc   dkim=kim_fac/rmin
c: Set max IM[kn] to allow for min range of interest rmin(km) (see pp.131,150):
      dkim=db_cut/(8685.9*rmin)
c: kim_fac will be used with range in m in mode_field:
      kim_fac=db_cut/8.6859
c
      return
      end
ccc
      subroutine hunt(xx,n,x,jlo)
      integer*4 jlo,n
      real*4 x,xx(n)
c fbvfbv      real*8 x,xx(n)
      integer*4 inc,jhi,jm
      logical ascnd
c
c: EKW FIX (always ascending for us, but tricked when n=1):
cc    ascnd=xx(n) .gt. xx(1)
      ascnd=xx(n) .ge. xx(1)
      if(jlo .le. 0 .or. jlo .gt. n) then
         jlo=0
         jhi=n+1
         goto 3
      endif
      inc=1
      if(x .ge. xx(jlo) .eqv. ascnd) then
1        jhi=jlo+inc
         if(jhi .gt. n)then
            jhi=n+1
         elseif(x .ge. xx(jhi) .eqv. ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
         endif
      else
         jhi=jlo
2        jlo=jhi-inc
         if(jlo .lt. 1) then
            jlo=0
         elseif(x .lt. xx(jlo) .eqv. ascnd) then
            jhi=jlo
            inc=inc+inc
            goto 2
         endif
      endif
3     if(jhi-jlo .eq. 1) return
      jm=(jhi+jlo)/2
      if(x .gt. xx(jm) .eqv. ascnd) then
         jlo=jm
      else
         jhi=jm
      endif
      goto 3
c
      end
