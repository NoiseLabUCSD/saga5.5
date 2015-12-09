      subroutine zref_chng
c
c: Changes reference depth to kduct'th duct.
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 j,nsvold,isvold,jx,jx1
c
c: Redo all computations involving reference depth zsvp(nsvmin):
c: From svp_check: Reset ref depth for reflection coefficient calcs:
      nsvold=nsvmin
      isvold=isvmin
      nsvmin=jduct(1,kduct)
      isvmin=jduct(2,kduct)
      rho_duct=geo(isvmin,3,nsvmin)
      jflu(1,1)=nsvmin
      jflu(2,1)=nsvmin
      jflu(1,5)=isvmin-1
      jflu(2,5)=isvmin-2
c: EKW FIX 12/13/95: Account for fact that a fake layer with different 
c: density may be there:
      if(jflu(1,2) .ne. jduct(4,kduct)) then
         jflu(1,2)=jduct(4,kduct)
         rholay(1)=geo(2,3,jflu(1,2))/geo(1,3,jsol(1,1))
      endif
      if(jflu(2,2) .ne. jduct(5,kduct)) then
         jflu(2,2)=jduct(5,kduct)
         rholay(2)=geo(1,3,jflu(2,2))/geo(2,3,jsol(2,1))
      endif
c: Check if this is a fake layer used for seismic modes:
c: Not necessary to go slower in seismic layers? 
cc    if(jduct(3,kduct) .gt. 0) then
cc       phfac0=max(8.e0,phfac)
cc    else
cc       phfac0=phfac
cc    endif
c
c: Invert density ratios at interfaces between old and new reference depths:
      do j=min(nsvold,nsvmin),max(nsvold,nsvmin)-1
         rhorat(j)=1.d0/rhorat(j)
      enddo
c: Change relative depths for mode function depths:
c: FIX: 8-27-97.  For case where nsvold=nsvmin, but isvold .ne. isvmin:
cc    if(nsvold .lt. nsvmin) then
      if(nsvold .le. nsvmin) then
         do j=nsvold+isvold-1,nsvmin+isvmin-2
            jx1=jzmx(j)-1
            do jx=1,nzmx(j)
               zmx(jx1+jx)=h(j) - zmx(jx1+jx)
            enddo
         enddo
      else
         do j=nsvmin+isvmin-1,nsvold+isvold-2
            jx1=jzmx(j)-1
            do jx=1,nzmx(j)
               zmx(jx1+jx)=h(j) - zmx(jx1+jx)
            enddo
         enddo
      endif
      jlmin=min(nsvmin,jsr2j(1))
      jlmax=max(nsvmin,jsr2j(nzsr))
      do j=1,nlay
         iiww(j)=0
         if(j .ge. jlmin .and. j .le. jlmax) iiww(j)=1
c: Temp:
         iiww(j)=1
      enddo
c
      xkref=xk(isvmin,nsvmin)
      xkrat_ref(1)=xkref/kw0
      kw=dreal(xkref)
      cref=w/kw0
c
c: kduct0 holds current duct:
      kduct0=kduct
c
      return
      end
ccc
      subroutine duct_dupl(nm_ok,ndup,r1r2,phiz,dphiz,iduct)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer nm_ok,ndup,iduct,jm0,jm1,nm_temp,jm2,inc,
     .        maxmode
      complex*8 phiz(nzsr,nm_ok),dphiz(nzsr,nm_ok)
      complex*16 kn_found,knd,r1r2(3,4),rr0(3,4),dW_dk
      real*8 magsq
      real*4 pct_diff
c
      if(iiccw .gt. 0) then
         jm1=1
         jm2=nm_ok
         inc=1
      else
         jm1=nm_ok
         jm2=1
         inc=-1
      endif
      kn_found=kn(nmode)
      do jm0=jm1,jm2,inc
         knd=kn(jm0)-kn_found
c: Check if eigenvalues extremely close:
         if(magsq(knd) .le. errdk2) then
c: Check if duct minima are nearly identical. If so, possibility that modes
c: might still be distinct:
c
c: Find current mode more precisely:
            call eig_final(kn(nmode),r1r2,.01d0*errdkms,0,kw0,
     .         iish,iifail,jjfail)
            if(iifail .eq. 1) return
            call eig_enter(r1r2)
            dW_dk=-2.d0*r1r2(1,1)*gamiref*r1r2(2,3)
            call mode_fun(kn(nmode),r1r2(1,1),r1r2(1,2),dW_dk,phi,dphi,
     .         psi,dpsi,exp_gbs,nmode)
            if(jjfail.eq.1) return
c
c: Find duplicate mode more precisely:
            kduct=nzref(jm0)
c: Change ducts:
            call zref_chng
            nm_temp=nmode
            nmode=jm0
            rr0(1,4)=eig_char(1,nmode)
            rr0(2,4)=eig_char(2,nmode)
            call eig_final(kn(nmode),rr0,.01d0*errdkms,0,kw0,
     .         iish,iifail,jjfail)
            if(iifail .eq. 1) return
            call eig_enter(rr0)
            dW_dk=-2.d0*rr0(1,1)*gamiref*rr0(2,3)
            call mode_fun(kn(nmode),rr0(1,1),rr0(1,2),dW_dk,phi,dphi,
     .         psi,dpsi,exp_gbs,nmode)
            if(jjfail.eq.1) return
c: Change reference depth back:
            nmode=nm_temp
            kduct=nzref(nmode)
            call zref_chng
c
cpln            write(6,*)
cpln            write(6,*)'nmode,ndup,nm_ok,jm0: ',nmode,ndup,nm_ok,jm0
cpln            write(6,*)'jm1,jm2,inc,iiccw: ',jm1,jm2,inc,iiccw
            maxmode=max(jm0,nmode)
            call mfun_comp(nzsr,nmode,jm0,maxmode,phiz,dphiz,nduct,
     .         mzduct,pct_diff)
            if(iidiag .ge. 2) print *,'Dup eigenvalues in different ',
     .         'ducts: pct_diff = ',pct_diff
c
            if(pct_diff .lt. 1.e-2) then
               if(iidiag .ge. 2) then
                  print *,'Duplicate mode found for duct: ',iduct,
     .               kn(jm0),kn(nmode)
               endif
               ndup=ndup + 1
               return
            endif
            if(iidiag .ge. 2) print *,'Close eigenvalues found: ',
     .         kn(jm0),kn(nmode)
         endif
c: EKW FIX. Check imaginary parts also to accommodate evanescent contours
c: that oscillate sideways:
         if(iiccw .gt. 0) then
            if(dreal(knd).lt.0.d0 .and. dimag(knd).gt.0.d0) goto 99
         else
            if(dreal(knd).gt.0.d0 .and. dimag(knd).lt.0.d0) goto 99
         endif
      enddo
c     
99    continue
c: Place new mode into nm_put spot:
      nm_put=nm_put + 1
      kn(nm_put)=kn(nmode)
      if(nm_put .gt. nmode) then
         if(iiwrite .gt. 0)
     .     write(6,*)'nm_put > nmode: ',nm_put,nmode
      end if
c
      call eig_insert(nmode,nm_put,eig_char,nzref,ncalc,mode_phz,
     .   phi,dphi,psi,dpsi,exp_gbs,iishn,max(nm_put,nmode),nzsr)
c
c: New: reset ndup so that two dups in a row must be found to stop looking:
      ndup=0
c
      return
      end
ccc
      subroutine eig_insert(jfrom,jto,eig_char,nzref,ncalc,mode_phz,
     .   phi,dphi,psi,dpsi,exp_gbs,iishn,nmode,nzsr)
c
c: Inserts mode into correct place and moves other up:
      implicit none
      integer*4 jfrom,jto,nmode,nzref(nmode),ncalc(nmode),
     .   iishn(nmode),
     .   nzsr,ji,jzsr
      real*4 mode_phz(3,nmode)
      complex*16 eig_char(5,nmode)
      complex*8 phi(nzsr,nmode),dphi(nzsr,nmode),psi(nzsr,nmode),
     .   dpsi(nzsr,nmode)
      real*8 exp_gbs(nzsr,nmode)
c
      ncalc(jto)=ncalc(jfrom)
      mode_phz(1,jto)=mode_phz(1,jfrom)
      mode_phz(2,jto)=mode_phz(2,jfrom)
      mode_phz(3,jto)=mode_phz(3,jfrom)
      nzref(jto)=nzref(jfrom)
      iishn(jto)=iishn(jfrom)
      do ji=1,5
         eig_char(ji,jto)=eig_char(ji,jfrom)
      enddo
      do jzsr=1,nzsr
         phi(jzsr,jto)=phi(jzsr,jfrom)
         dphi(jzsr,jto)=dphi(jzsr,jfrom)
         psi(jzsr,jto)=psi(jzsr,jfrom)
         dpsi(jzsr,jto)=dpsi(jzsr,jfrom)
         exp_gbs(jzsr,jto)=exp_gbs(jzsr,jfrom)
      enddo
c
      return
      end
ccc
      subroutine eig_sort(nmode,nzsr,kn_indx,eig_char,nzref,
     .   ncalc,iishn,mode_phz,phi,dphi,psi,dpsi,exp_gbs)
c
c: Inserts mode into correct place and moves other up:
      implicit none
      integer*4 nmode,kn_indx(nmode),nzsr,nzref(nmode),ncalc(nmode),
     .   iishn(nmode),jsave,jto,jfrom,iic,j,nm1
      real*4 mode_phz(3,nmode)
      complex*16 eig_char(5,nmode)
      complex*8 phi(nzsr,nmode),dphi(nzsr,nmode),psi(nzsr,nmode),
     .   dpsi(nzsr,nmode)
      real*8 exp_gbs(nzsr,nmode)
c
c: Save first mode:
      jfrom=1
      jsave=1
      nm1=nmode + 1
      call eig_insert(jsave,nmode+1,eig_char,nzref,ncalc,
     .   mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
      do j=1,nmode
         jto=jfrom
         jfrom=kn_indx(jto)
         if(jto .eq. jfrom) then
            iic=1
         elseif(jsave .eq. jfrom) then
            call eig_insert(nmode+1,jto,eig_char,nzref,ncalc,
     .         mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
            iic=1
         else
            call eig_insert(jfrom,jto,eig_char,nzref,ncalc,
     .         mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
            iic=0
         endif
         kn_indx(jto)=0
         if(iic .eq. 1) then
            jsave=jsave+1
            do while(jsave .le. nmode .and. kn_indx(jsave) .eq. 0)
               jsave=jsave+1
            enddo
            call eig_insert(jsave,nmode+1,eig_char,nzref,ncalc,
     .         mode_phz,phi,dphi,psi,dpsi,exp_gbs,iishn,nm1,nzsr)
            jfrom=jsave
         endif
      enddo
c
      return
      end
ccc
      subroutine mfun_comp(nzsr,nmode,jm0,maxmode,phiz,dphiz,nduct,
     .   mzduct,pct_diff)
c
c: Compares mode functions at the duct depths to see if eigenvalue is same.
c
      implicit none
      integer*4 nzsr,nmode,jm0,nduct,mzduct(nduct),jd,jzsr,
     .          maxmode
      real*4 psum,dpsum,pct_diff,psum1,psum2,dpsum1,dpsum2
      complex*8 phiz(nzsr,maxmode),dphiz(nzsr,maxmode)
c
      psum1=0.0
      psum2=0.0
      dpsum1=0.0
      dpsum2=0.0
      psum=0.0
      dpsum=0.0
      do jd=1,nduct
         jzsr=mzduct(jd)
c: Take into account possibility that mode function differ by -1:
         call sum_sq(phiz(jzsr,jm0),phiz(jzsr,nmode),1.0,psum1)
         call sum_sq(phiz(jzsr,jm0),phiz(jzsr,nmode),-1.0,psum2)
         call sum_sq(phiz(jzsr,jm0),phiz(jzsr,jm0),0.0,psum)
         call sum_sq(dphiz(jzsr,jm0),dphiz(jzsr,nmode),1.0,dpsum1)
         call sum_sq(dphiz(jzsr,jm0),dphiz(jzsr,nmode),-1.0,dpsum2)
         call sum_sq(dphiz(jzsr,jm0),dphiz(jzsr,jm0),0.0,dpsum)
      enddo
cxx   pct_diff=min(psum1,psum2)/max(1.e-200,psum)
cxx   pct_diff=min(psum1,psum2)/max(1.e-200,psum)
      pct_diff=max(min(psum1,psum2)/max(1.e-37,psum),
     .   min(dpsum1,dpsum2)/max(1.e-35,dpsum))
c
      return
      end
ccc
      subroutine sum_sq(z1,z2,sg,psum)
c
      implicit none
      complex*8 z1,z2,diff
      real*4 sg,psum
c
      diff=z1 - sg*z2
      psum=psum + real(diff)*real(diff) + aimag(diff)*aimag(diff)
c
      return
      end
ccc
      subroutine zduct_chng(phiz,dphiz,jm,jm2,iichng,k)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jm,jm2,iichng,jzsr,nsv,isv,kdmax,jd
      complex*8 phiz(nzsr,jm2),dphiz(nzsr,jm2),dphi_fac
      real*4 magsq_c8,phimag(25),fac
      complex*16 k
c
      kdmax=1
      do jd=1,nduct
         jzsr=mzduct(jd)
         nsv=jduct(1,jd)
         isv=jduct(2,jd)
         dphi_fac=dphiz(jzsr,jm2)/gami(1,isv,nsv)
         phimag(jd)=magsq_c8(phiz(jzsr,jm2)) + magsq_c8(dphi_fac)
         if(phimag(jd) .gt. phimag(kdmax)) kdmax=jd
      enddo
c
      iichng=0
      isv=jduct(2,kdmax)
      nsv=jduct(1,kdmax)
      fac=2.
      if(nsv .eq. nlay) fac=5.
      if(phimag(kdmax)/phimag(kduct) .gt. fac .and.
     .   dreal(k) .lt. dreal(xk(isv,nsv))) then
c: Change reference duct larger magnitude in other duct exists:
         kduct=kdmax
         iichng=1
      elseif(nzref(jm) .ne. kduct) then
         if(phimag(nzref(jm))/phimag(kduct) .gt. .5) then
c: Change reference duct back to standard one if magnitudes are not
c: clearly different:
            kduct=nzref(jm)
            iichng=1
         endif
      endif
c
      return
      end
