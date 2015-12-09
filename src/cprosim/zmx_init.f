      subroutine zmx_init
c
c: Initializes depth arrays for computation of mode functions.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jptr,jj,j,jj1,jj2,jsr,ii,jx1
      real*8 z1,z2,rhfac,fac
      complex*16 c1inv,c1invsq,c2inv,c2invsq,cfac
c: fac=2*pi*1000*20*log(e):
      data fac/5.457505415367364d+04/
c
      jptr=0
      jj=1 
      z1=-1.d+20
      do j=1,jflu(2,1)+jflu(2,5)
         z2=zdep(j)
c: For receiver depths equal to interface depths, put receiver in layer
c: closest to top ocean layer:
         if(j .ge. jsurf) z2=z2 + 1.d-10*abs(z2)
         do while(zsr(jj) .lt. z1 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj1=jj
         do while(zsr(jj) .lt. z2 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj2=jj-1
         jzmx(j)=jptr+1
         nzmx(j)=jj2-jj1+1
         rhfac=(geo(1,3,j) - geo(2,3,j))/dmax1(1.d-20,h(j))
cxx      c1inv=1./geo(2,1,j)**2
cxx      cfac=(1./geo(1,1,j)**2 - c1inv)/dmax1(1.d-20,h(j))
c: see p. 66 (from SAFARI manual for attenuation in k=w/c):
         c1inv=dcmplx(1.d0,geo(2,4,j)*geo(2,1,j)/fac)/geo(2,1,j)
         c1invsq=c1inv*c1inv
         c2inv=dcmplx(1.d0,geo(1,4,j)*geo(1,1,j)/fac)/geo(1,1,j)
         c2invsq=c2inv*c2inv
         cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
         do jsr=jj2,jj1,-1
            jptr=jptr + 1
            zmx(jptr)=zdep(j)-zsr(jsr)
            zmx_im_gbs(jptr)=-zsr_im_gbs(jsr)
c: jsrmx(jptr)=jsr, the index of depth zmx(jptr) in the array zsr(1:nzsr):
            jsrmx(jptr)=jsr
            rho_sr(jsr)=geo(2,3,j) + zmx(jptr)*rhfac
            cp_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            cs_sr(jsr)=(0.d0,0.d0)
c: Pointer from zmx to zm [find depth zsr(jsr) in zmx(mx_m(jsr))]:
            mx_m(jsr)=jptr
c: Pointer from zsr to geo [find s/r depth zsr(jsr) in geo(:,:,j)]:
            jsr2j(jsr)=j
         enddo
         if(iisol(j) .eq. 1 .and. jj2 .ge. jj1) then
            c1inv=dcmplx(1.d0,geo(2,5,j)*geo(2,2,j)/fac)/geo(2,2,j)
            c1invsq=c1inv*c1inv
            c2inv=dcmplx(1.d0,geo(1,5,j)*geo(1,2,j)/fac)/geo(1,2,j)
            c2invsq=c2inv*c2inv
            cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
            do jsr=jj2,jj1,-1
               kksh(jsr)=1
               cs_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            enddo
         endif
cxx   print *,'j,z1,z2,h,zmx = ',j,z1,z2,h(j),(zmx(jx1),jx1=
cxx  .   jzmx(j),jzmx(j)+nzmx(j)-1)
         z1=z2
      enddo
c
c: For layers below reference depth:
      do j=jflu(1,1)+jflu(1,5),nlay
         if(j .ne. nlay) then
            z2=zdep(j)
            if(j .ge. jsurf) z2=z2 + 1.d-10*abs(z2)
         else
            z2=1.d20
         endif
         do while(zsr(jj) .lt. z1 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj1=jj
         do while(zsr(jj) .lt. z2 .and. jj .le. nzsr)
            jj=jj+1
         enddo
         jj2=jj-1
         jzmx(j)=jptr+1
         nzmx(j)=jj2-jj1+1
cxx      c1inv=1./geo(1,1,j)**2
cxx      cfac=(1./geo(2,1,j)**2 - c1inv)/dmax1(1.d-20,h(j))
         c1inv=dcmplx(1.d0,geo(1,4,j)*geo(1,1,j)/fac)/geo(1,1,j)
         c1invsq=c1inv*c1inv
         c2inv=dcmplx(1.d0,geo(2,4,j)*geo(2,1,j)/fac)/geo(2,1,j)
         c2invsq=c2inv*c2inv
         cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
         rhfac=(geo(2,3,j) - geo(1,3,j))/dmax1(1.d-20,h(j))
         do jsr=jj1,jj2
            jptr=jptr + 1
            zmx(jptr)=zsr(jsr)-zdep(j-1)
            zmx_im_gbs(jptr)=zsr_im_gbs(jsr)
c: jsrmx(jptr)=jsr, the index of depth zmx(jptr) in the array zsr(1:nzsr):
            jsrmx(jptr)=jsr
            rho_sr(jsr)=geo(1,3,j) + zmx(jptr)*rhfac
            cp_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            cs_sr(jsr)=(0.d0,0.d0)
c: Pointer from zmx to zm [find depth zsr(jsr) in zmx(mx_m(jsr))]:
            mx_m(jsr)=jptr
c: Pointer from zsr to geo [find s/r depth zsr(jsr) in geo(:,:,j)]:
            jsr2j(jsr)=j
         enddo
         if(iisol(j) .eq. 1 .and. jj1 .le. jj2) then
            c1inv=dcmplx(1.d0,geo(1,5,j)*geo(1,2,j)/fac)/geo(1,2,j)
            c1invsq=c1inv*c1inv
            c2inv=dcmplx(1.d0,geo(2,5,j)*geo(2,2,j)/fac)/geo(2,2,j)
            c2invsq=c2inv*c2inv
            cfac=(c2invsq - c1invsq)/dmax1(1.d-100,h(j))
            do jsr=jj1,jj2
               kksh(jsr)=1
               cs_sr(jsr)=1.d0/cdsqrt(c1invsq + zmx(jptr)*cfac)
            enddo
         endif
cxx   print *,'j,z1,z2,h,zmx = ',j,z1,z2,h(j),(zmx(jx1),jx1=
cxx  .   jzmx(j),jzmx(j)+nzmx(j)-1)
         z1=z2
      enddo
c
c: Total number of depths in zmx, which includes source, receiver, and
c: layer interface depths:
      nzmxtot=jptr
c
c: Initialize s-wave potentials to zero in water and fluid layers:
      do ii=1,2
         do j=jflu(ii,1),jflu(ii,2),jflu(ii,3)
            do jx1=jzmx(j),jzmx(j)+nzmx(j)-1
               psix(jx1)=(0.,0.)
               dpsix(jx1)=(0.,0.)
            enddo
         enddo
      enddo
c: iiww(j) is a flag that transmission coefficients need to be computed in
c: the layer j:
      if(nzsr .gt. 0) then
         jlmin=min(nsvmin,jsr2j(1))
         jlmax=max(nsvmin,jsr2j(nzsr))
      else
         jlmin=nsvmin
         jlmax=nsvmin
      endif
      do j=1,nlay
         iiww(j)=0
         if(j .ge. jlmin .and. j .le. jlmax) iiww(j)=1
c: Temp (see also duct_check):
         iiww(j)=1
      enddo
c: Find maximum receiver sound speed for use in mode_field:
      crmax=cfmin
      do jsr=1,nzsr
         crmax=dmax1(crmax,dreal(cp_sr(jsr)))
      enddo
c
c: DON'T DO THIS ANY MORE SINCE WE CAN HAVE MULTIPLE SOURCES:
c: Normalize densities to density at source:
cxx   rho_src=rho_sr(mzsrc(1))
cxx   do j=1,nzsr
cxx      rho_sr(j)=rho_sr(j)/rho_src
cxx   enddo
c
      return
      end
