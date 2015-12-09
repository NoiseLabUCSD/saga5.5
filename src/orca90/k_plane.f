      subroutine k_plane(ir,jrun,nrunx)
c
c: Computes and outputs HDF files for complex k-plane plots of the real
c: and imaginary parts of ln(R1*R2) and/or ln(1-R1*R2).
c
      use parms_com
      use i_o_com
      use gen_com
c
      implicit none
      integer*4 ir,jrun,nrunx,ji,jr,iibad
      real*4 facr,faci,fac_kr,fac_ki
      complex*16 rr0(3,4),lnr1r2,lnW,psi_fac,arg,H0,kernel,psi_fac_neg
      real*8 pienorm,r_pbli,pierad
c
      call mem_lim(nreal,NSRMAX,MLINE,LML,'nreal',5,'NSRMAX',6,iibad,0)
      call mem_lim(nimag,NSRMAX,MLINE,LML,'nimag',5,'NSRMAX',6,iibad,0)
c
      f_hz=fkpl
      call freq_init
c: Don't allow sheet changes at branch cuts:
      iich=0
      iich_ref=0
c
      iish0(1,1)=iishp
      iish0(2,1)=iishp
      iish0(1,2)=iishs
      iish0(2,2)=iishs
      iishr0(1)=1
      iishr0(2)=1
cc    print *,'enter iishref:'
cc    read(5,*) iishr0(1),iishr0(2)
      call sheet_init(xkh,0,iish0,iishr0)
      if(kduc .eq. 0) then
         kduct=indx_duct(nduct)
      else
         if(kduc .gt. nduct .or. kduc .lt. 1) then
            print *,'Illegal duct number: kduc,nduct = ',kduc,nduct
            return
         endif
         kduct=kduc
      endif
      call zref_chng
      print *,'For k-plane plots, kduct,z,c = ',kduct,
     .   zduct(kduct),geo(isvmin,1,nsvmin)
c
      pierad=pie/180.
      if(iivar .eq. 1) then
         pienorm=pierad
      else
         pienorm=1.d0
      endif
c: For k/kw0, use absolute minimum in SVP, rather than other duct if requested:
      if(iiform .eq. 1) t hen
         fac_kr=kw0
         fac_ki=kw0
      elseif(iiform .eq. 2) then
         fac_kr=1.
         fac_ki=1.
      elseif(iiform .eq. 3) then
         fac_kr=kw0
         fac_ki=1./8685.9
      endif
c
      facr=(xkhr2-xkhr1)/max(nreal-1,1)
      do jr=1,nreal
         xkhr(jr)=xkhr1 + (jr-1)*facr
      enddo
      faci=(xkhi2-xkhi1)/max(nimag-1,1)
      do ji=1,nimag
         xkhi(ji)=xkhi1 + (ji-1)*faci
      enddo
      if(iikpl .eq. 4 .or. iikpl .eq. 5) then
         if(rmin .eq. 0.) then
            print *,'rmin = 0 in PBLI for k-plane.  Enter rmin in km: '
            read *,rmin
         endif
         r_pbli=1000.*rmin
      endif
c
      allocate(r4mat1(nimag,nreal))
      allocate(r4mat2(nimag,nreal))
      allocate(r4mat3(nimag,nreal))
      allocate(r4mat4(nimag,nreal))
      do ji=1,nimag
         do jr=1,nreal
            xkh=dcmplx(fac_kr*xkhr(jr),fac_ki*xkhi(ji))
            call r1r2_calc(xkh,rr0,iivar,1)
            if(iikpl .eq. 1 .or. iikpl .eq. 3) then
               lnr1r2=rr0(iivar,4)
               r4mat1(ji,jr)=real(lnr1r2)
               r4mat2(ji,jr)=dimag(lnr1r2)/pienorm
            elseif(iikpl .eq. -1 .or. iikpl .eq. -3) then
               lnr1r2=rr0(iivar,3)
               r4mat1(ji,jr)=real(lnr1r2)
               r4mat2(ji,jr)=dimag(lnr1r2)
            endif
            if(iikpl .eq. 2 .or. iikpl .eq. 3) then
               lnW=cdlog(1.d0 - rr0(iivar,3))
               r4mat3(ji,jr)=real(lnW)
               r4mat4(ji,jr)=dimag(lnW)/pierad
            elseif(iikpl .eq. -2 .or. iikpl .eq. -3) then
               lnW=1.d0 - rr0(iivar,3)
               r4mat3(ji,jr)=real(lnW)
               r4mat4(ji,jr)=dimag(lnW)
            endif
c: Plot integrand kernel for looking at Pekeris branch line integral
cc          if(iikpl .eq. 4) then
cc             psi_fac=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
cc   .            (gamiref*(1.d0 - rr0(1,3)))
cc             arg=xkh*r_pbli
cc             if(cdabs(arg) .gt. 5.) then
cc                H0=cdsqrt(dcmplx(0.d0,-2.d0)/(pie*xkh*r_pbli))*
cc   .               cdexp(dcmplx(0.,1.)*xkh*r_pbli)
cc             else
cc                call cdhankel(arg,1.d-6,H0)
cc             endif
cc             kernel=cdlog(0.5*psi_fac*H0*xkh)
cc             r4mat1(ji,jr)=real(kernel)
cc             r4mat2(ji,jr)=dimag(kernel)/pierad
cc          elseif(iikpl .eq. 5) then
cc             psi_fac=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
cc   .            (gamiref*(1.d0 - rr0(1,3)))
cc             arg=xkh*r_pbli
cc             if(cdabs(arg) .gt. 5.) then
cc                H0=cdsqrt(dcmplx(0.d0,-2.d0)/(pie*xkh*r_pbli))*
cc   .               cdexp(dcmplx(0.,1.)*xkh*r_pbli)
cc             else
cc                call cdhankel(arg,1.d-6,H0)
cc             endif
cc             kernel=0.5*psi_fac*H0*xkh
cc             iish0(1,1)=-iish0(1,1)
cc             call sheet_init(xkh,0,iish0,iishr0)
cc             call r1r2_calc(xkh,rr0,iivar,1)
cc             psi_fac_neg=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
cc   .            (gamiref*(1.d0 - rr0(1,3)))
cc             kernel=cdlog(0.5*H0*xkh*(psi_fac_neg - psi_fac))
cc             r4mat1(ji,jr)=real(kernel)
cc             r4mat2(ji,jr)=dimag(kernel)/pierad
cc             iish0(1,1)=-iish0(1,1)
cc             call sheet_init(xkh,0,iish0,iishr0)
cc          endif
         enddo
      enddo
      if(iikpl .gt. 0 .and. iiwr .eq. 2) then
         if(iikpl .eq. 1 .or. iikpl .ge. 3) then
            do ji=1,nimag
               do jr=1,nreal
                  if(r4mat2(ji,jr) .lt. 0.) r4mat2(ji,jr)=
     .               r4mat2(ji,jr) + 360.
               enddo
            enddo
         endif
         if(iikpl .eq. 3) then
            do ji=1,nimag
               do jr=1,nreal
                  if(r4mat4(ji,jr) .lt. 0.) r4mat4(ji,jr)=
     .               r4mat4(ji,jr) + 360.
               enddo
            enddo
         endif
      endif
cpg      if(abs(iikpl) .eq. 1 .or. abs(iikpl) .eq. 3) then
cpg        call hdf_write_gen(2+ir,outroot,lout,fcw,nfcw,xmode,
cpg          call hdf_write_gen(2+ir,outroot,lout,fcw,nfcw,xmode,
cpg        call hdf_write_gen(2+ir,outroot,lout,xkhi,nimag,xkhr,nreal,
cpg     .      var_ax,nrunx,xkhr,1,xkhr,1,r4mat1,1,nimag,1,nreal,jrun,jrun,
cpg     .      1,1,1,1,'Im(k)',5,'Re(k)',5,'Parameter',9,' ',1,' ',1,
cpg     .      'Re[ln(R1R2)]',12,6,1,2)
cpg         call hdf_write_gen(2+ir,outroot,lout,xkhi,nimag,xkhr,nreal,
cpg     .      var_ax,nrunx,xkhr,1,xkhr,1,r4mat2,1,nimag,1,nreal,jrun,jrun,
cpg     .      1,1,1,1,'Im(k)',5,'Re(k)',5,'Parameter',9,' ',1,' ',1,
cpg    .      'Im[ln(R1R2)]',12,7,1,2)
cpg      endif
cpg      if(abs(iikpl) .eq. 2 .or. abs(iikpl) .eq. 3) then
cpg         call hdf_write_gen(2+ir,outroot,lout,xkhi,nimag,xkhr,nreal,
cpg     .      var_ax,nrunx,xkhr,1,xkhr,1,r4mat3,1,nimag,1,nreal,
cpg     .      jrun,jrun,1,1,1,1,'Im(k)',5,'Re(k)',5,'Parameter',9,' ',1,
cpg     .      ' ',1,'Re[ln(1-R1R2)]',14,8,1,2)
cpg         call hdf_write_gen(2+ir,outroot,lout,xkhi,nimag,xkhr,nreal,
cpg     .      var_ax,nrunx,xkhr,1,xkhr,1,r4mat4,1,nimag,1,nreal,
cpg     .      jrun,jrun,1,1,1,1,'Im(k)',5,'Re(k)',5,'Parameter',9,' ',1,
cpg    .      ' ',1,'Im[ln(1-R1R2)]',14,9,1,2)
cpg      endif
c
      deallocate(r4mat1)
      deallocate(r4mat2)
      deallocate(r4mat3)
      deallocate(r4mat4)
c
      return
      end
