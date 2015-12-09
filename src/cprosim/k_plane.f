      subroutine k_plane(nre,nim,mat1,mat2,mat3,mat4)
c
c: Computes and outputs HDF files for complex k-plane plots of the real
c: and imaginary parts of ln(R1*R2) and/or ln(1-R1*R2).
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
c
      integer*4 nre,nim,jr,ji,iibad
      real*4 mat1(nre,nim),mat2(nre,nim),mat3(nre,nim),mat4(nre,nim),
     .   facr,faci,fac_kr,fac_ki
      complex*16 rr0(3,4),lnr1r2,lnW,psi_fac,arg,H0,kernel,psi_fac_neg
      real*8 pienorm,r_pbli,pierad
c
      call mem_lim(nre,NSRMAX,MLINE,LML,'nreal',5,'NSRMAX',6,iibad,0)
      call mem_lim(nim,NSRMAX,MLINE,LML,'nimag',5,'NSRMAX',6,iibad,0)
      call mem_lim(nre*nim,NTLMAX,MLINE,LML,'nre*nim',7,'NTLMAX',6,
     .   iibad,1)
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
      facr=(xkhr2-xkhr1)/max(nre-1,1)
      do jr=1,nre
         xkhr(jr)=xkhr1 + (jr-1)*facr
      enddo
      faci=(xkhi2-xkhi1)/max(nim-1,1)
      do ji=1,nim
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
      do ji=1,nim
         do jr=1,nre
            xkh=dcmplx(fac_kr*xkhr(jr),fac_ki*xkhi(ji))
            call r1r2_calc(xkh,rr0,iivar,1,jjfail)
            if(jjfail .gt. 0) return
            if(iikpl .eq. 1 .or. iikpl .eq. 3) then
               lnr1r2=rr0(iivar,4)
               mat1(jr,ji)=real(lnr1r2)
               mat2(jr,ji)=dimag(lnr1r2)/pienorm
            elseif(iikpl .eq. -1 .or. iikpl .eq. -3) then
               lnr1r2=rr0(iivar,3)
               mat1(jr,ji)=real(lnr1r2)
               mat2(jr,ji)=dimag(lnr1r2)
            endif
            if(iikpl .eq. 2 .or. iikpl .eq. 3) then
               lnW=cdlog(1.d0 - rr0(iivar,3))
               mat3(jr,ji)=real(lnW)
               mat4(jr,ji)=dimag(lnW)/pierad
            elseif(iikpl .eq. -2 .or. iikpl .eq. -3) then
               lnW=1.d0 - rr0(iivar,3)
               mat3(jr,ji)=real(lnW)
               mat4(jr,ji)=dimag(lnW)
            endif
c: Plot integrand kernel for looking at Pekeris branch line integral
            if(iikpl .eq. 4) then
               psi_fac=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
     .            (gamiref*(1.d0 - rr0(1,3)))
               arg=xkh*r_pbli
               if(cdabs(arg) .gt. 5.) then
                  H0=cdsqrt(dcmplx(0.d0,-2.d0)/(pie*xkh*r_pbli))*
     .               cdexp(dcmplx(0.,1.)*xkh*r_pbli)
               else
                  call cdhankel(arg,1.d-6,H0)
               endif
               kernel=cdlog(0.5*psi_fac*H0*xkh)
               mat1(jr,ji)=real(kernel)
               mat2(jr,ji)=dimag(kernel)/pierad
            elseif(iikpl .eq. 5) then
               psi_fac=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
     .            (gamiref*(1.d0 - rr0(1,3)))
               arg=xkh*r_pbli
               if(cdabs(arg) .gt. 5.) then
                  H0=cdsqrt(dcmplx(0.d0,-2.d0)/(pie*xkh*r_pbli))*
     .               cdexp(dcmplx(0.,1.)*xkh*r_pbli)
               else
                  call cdhankel(arg,1.d-6,H0)
               endif
               kernel=0.5*psi_fac*H0*xkh
               iish0(1,1)=-iish0(1,1)
               call sheet_init(xkh,0,iish0,iishr0)
               call r1r2_calc(xkh,rr0,iivar,1,jjfail)
               if(jjfail .gt. 0) return
               psi_fac_neg=-(1.d0 + rr0(1,1))*(1.d0 + rr0(1,2))/
     .            (gamiref*(1.d0 - rr0(1,3)))
               kernel=cdlog(0.5*H0*xkh*(psi_fac_neg - psi_fac))
               mat1(jr,ji)=real(kernel)
               mat2(jr,ji)=dimag(kernel)/pierad
               iish0(1,1)=-iish0(1,1)
               call sheet_init(xkh,0,iish0,iishr0)
            endif
         enddo
      enddo
      if(iikpl .gt. 0 .and. iiwr .eq. 2) then
         if(iikpl .eq. 1 .or. iikpl .ge. 3) then
            do ji=1,nim
               do jr=1,nre
                  if(mat2(jr,ji) .lt. 0.) mat2(jr,ji)=mat2(jr,ji) + 360.
               enddo
            enddo
         endif
         if(iikpl .eq. 3) then
            do ji=1,nim
               do jr=1,nre
                  if(mat4(jr,ji) .lt. 0.) mat4(jr,ji)=mat4(jr,ji) + 360.
               enddo
            enddo
         endif
      endif
      if(abs(iikpl) .eq. 1 .or. abs(iikpl) .eq. 3) then
         call out_writex(outroot,lout,SUFX//'re',3,mat1,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'im',3,mat2,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Im[ln(R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
      endif
      if(abs(iikpl) .eq. 2 .or. abs(iikpl) .eq. 3) then
         call out_writex(outroot,lout,SUFX//'reW',4,mat3,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(1-R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'imW',4,mat4,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Im[ln(1-R1*R2)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
      endif
      if(abs(iikpl) .eq. 4) then
         call out_writex(outroot,lout,SUFX//'lnker',6,mat1,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(kernel)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'argker',7,mat2,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Im[ln(kernel)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
      endif
      if(abs(iikpl) .eq. 5) then
         call out_writex(outroot,lout,SUFX//'lnkerdif',9,mat1,xkhi,xkhr,
     .      nim,nre,kilab,krlab,0.,0.,0.,0.,2,'Re[ln(kernel)] vs. k',
     .      ' ',' ',' ','f7.4','f7.4','f7.2',ncall)
         call out_writex(outroot,lout,SUFX//'argkerdif',10,mat2,xkhi,
     .      xkhr,nim,nre,kilab,krlab,0.,0.,0.,0.,2,
     .      'Im[ln(kernel)] vs. k',' ',' ',' ','f7.4','f7.4','f7.2',
     .      ncall)
      endif
c
      return
      end
