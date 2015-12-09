      subroutine bb_field(knx,phiz,dphiz,dpsiz,expz_gbs,jmk,tfz,
     .   jfbb,jmo)
c
c: Computes field tf for a single mode characterized by eigenvalue k and mode 
c: functions phi and dpsi.
c
      use parms_com
      use i_o_com
      use gen_com
      integer*4 jmk,jfbb,jmo,jr,jsr,kk,kk0,jsrc,jrec,j,jzs,msrc,nrmax
      complex*8 phiz(nzsr,jmk),dphiz(nzsr,jmk),dpsiz(nzsr,jmk),
     .   tfz(nfbb,nrec,nrng),dphi_z
      complex*16 knx,sqkn,cfac2,h_arg,H0,iknx,phi_src
      real*8 expz_gbs(nzsr,jmk),xn_beam,rsig_max,magsq
cc    integer*4 isv,jzsr,nsv,jd
cc    real*4 phiw_max,magsq_c8,ampl
c
c: Find mode function at source:
c: Normalize by density at source (9-1-94):
      jzs=1
      msrc=mzsrc(jzs)
      phi_src=phiz(msrc,jmk)/rho_sr(msrc)
      xn_beam=(w/cp_sr(msrc))*b_gbs(jzs) - expz_gbs(msrc,jmk)
      iknx=dcmplx(-dimag(knx),dreal(knx))
      sqkn=cdsqrt(knx)
      do jrec=1,nrec
         jsr=mzrec(jrec)
         phisr(jsr)=phiz(jsr,jmk)
         if(kksh(jsr) .eq. 1) then
c: Add shear wave potential contributions to recs in shear layers (see p. 117):
            phisr(jsr)=phisr(jsr) - 2.d0*ksm2_sr(jsr)*
     .         (knx*knx*phisr(jsr) + iknx*dpsiz(jsr,jmk))
         endif
      enddo
c
c: Find range beyond which mode is insignificant (ranges have been sorted):
      rsig_max=kim_fac/dmax1(1.d-20,dimag(knx)-kim_bb(jfbb))
      nrmax=0
      call hunt(range,nrng,rsig_max,nrmax)
c
      do jr=1,nrmax
         h_arg=knx*range(jr)
         if(magsq(h_arg) .gt. 25.d0) then
czs: Include normalization by exp(-xn_beam) here:
            cfac2=sq2pir(jr)*phi_src*
     .         cdexp(iknx*range(jr) - xn_beam)/sqkn
         else
            call cdhankel(h_arg,1.d-6,H0)
            cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
         endif
         kk0=krec_jr(jr)
         do kk=kk0+1,kk0+nrec_jr(jr)
            jsrc=jrec_jr(1,kk)
            jrec=jrec_jr(2,kk)
            jsr=mzrec(jrec)
            tfz(jfbb,jrec,jsrc)=tfz(jfbb,jrec,jsrc) + cfac2*phisr(jsr)
         enddo
      enddo
c
c: Keep track of strongest mode amplitude outside of halfspace:
cc    phiw_max=0.
cc    do jd=1,nduct
cc       nsv=jduct(1,jd)
cc       if(nsv .ne. nlay) then
cc          jzsr=mzduct(jd)
cc          isv=jduct(2,jd)
cc          dphi_z=dphiz(jzsr,jmk)/gami(1,isv,nsv)
cc          ampl=magsq_c8(phiz(jzsr,jmk)) + magsq_c8(dphi_z)
cc          phiw_max=max(phiw_max,ampl)
cc       endif
cc    enddo
cc    phim_bb(jfbb)=max(phim_bb(jfbb),phiw_max)
c
      if(iiout .ne. 0 .or. iimf .ne. 0) then
c: Fill phibb(nfbb,nzsr),dpsibb(nfbb,nzsr) with mode functions:
         do jsr=1,nzsr
            j=(jsr-1)*nfbb + jfbb
            phibb(jfbb,jsr)=phiz(jsr,jmk)
            dpsibb(jfbb,jsr)=dpsiz(jsr,jmk)
         enddo
      endif
c
      return
      end
