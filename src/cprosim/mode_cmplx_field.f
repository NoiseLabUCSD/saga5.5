      subroutine mode_cmplx_field(tfz,phiz,dpsiz,expz_gbs,jzs,jf,nf)
c
c: Computes field given mode eigenvalues and mode function values.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      integer*4 jf,jzs,jm,jr,jd,jsr,kk,kk0,jsrc,jrec,msrc,
     .          nrmax,nf,jj0,jj1,jj2
      complex*16 knx,iknx,iknr,sqkn,cfac2,phi_src,H0,h_arg,
     .   cfac,exp_sqk
      complex*8 phiz(nzsr,nmode),dpsiz(nzsr,nmode),
     .          tfz(nf*nzsr*nrec)
      real*8 expz_gbs(nzsr,nmode),rsig_max,hmsq,xn_b,xn_beam,
     .   magsq,rng_im
c: cfac = sqrt(2*pi)*exp(i*pi/4) (see p. 7, ORCA II):
c      data cfac/(1.77245385090552d0,1.77245385090552d0)/
c

      cfac=(1.77245385090552d0,1.77245385090552d0)
      msrc=mzsrc(jzs)
      if(iigbs .ne. 0) then
         rng_im=-b_gbs(jzs)*cos(th_gbs(jzs)*pie/180.d0)
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .              +dtiltvp(jd),rng_im)+dtiltvp(jd))
            end do
         enddo
      else
         rng_im=0.d0
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/dsqrt(range(jr)+dtiltvp(jd))
            end do
         enddo
      endif
c
      hmsq=25.d0
c: Initialize complex propagation loss array:
c      do jsrc=1,nsrc
c         do jrec=1,nrec
c            tfz(jrec,jsrc)=cmplx(0.e0,0.e0)
c         enddo
c      enddo
c
      xn_b=(w/dreal(cp_sr(msrc)))*b_gbs(jzs)
      do jm=1,nmode
         knx=kn(jm)
         iknx=dcmplx(-dimag(knx),dreal(knx))
         sqkn=cdsqrt(knx)
c: rsig_max is maximum range at which this mode is significant (50 dB down 
c: from strongest mode):
         rsig_max=kim_fac/dmax1(1.d-20,dimag(knx)-kim_min)
c: Obtain p-wave mode excitation at source depth:
         phi_src=phiz(msrc,jm)/rho_sr(msrc)
         xn_beam=xn_b - expz_gbs(msrc,jm)
c: Obtain p-wave mode excitations at receiver depths:
         do jrec=1,nrec
            jsr=mzrec(jrec)
c: rho included in mode_fun:
            phisr(jsr)=phiz(jsr,jm)
            if(kksh(jsr) .eq. 1) then
c: Add shear wave potential contributions to recs in shear layers (see p. 117):
               phisr(jsr)=phisr(jsr) - 2.d0*ksm2_sr(jsr)*
     .            (knx*knx*phisr(jsr) + iknx*dpsiz(jsr,jm))
            endif
         enddo
c: Find range beyond which mode is insignificant (ranges have been sorted):
         nrmax=0
         call hunt(range,nrng,rsig_max,nrmax)
         do jr=1,nrmax
cpln            h_arg=knx*dcmplx(range(jr),rng_im)
            kk0=krec_jr(jr)
            do kk=kk0+1,kk0+nrec_jr(jr)
               jsrc=jrec_jr(1,kk)
               jrec=jrec_jr(2,kk)
               h_arg=knx*dcmplx(range(jr)+dtiltvp(jrec),
     .               rng_im+dtiltvp(jrec))
               if(magsq(h_arg) .gt. hmsq) then
czs: Include normalization by exp(-xn_beam) here:
                  iknr=dcmplx(-dimag(h_arg),dreal(h_arg))
                  exp_sqk=cdexp(dcmplx(dreal(iknr)-xn_beam,dimag(iknr)))
     .                 /sqkn
                  cfac2=sq2pir(jr,jrec)*phi_src*exp_sqk
               else
                  call cdhankel(h_arg,1.d-6,H0)
cc             cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
                  cfac2=dcmplx(-pie*dimag(H0),pie*dreal(H0))*phi_src*
     .                 dexp(-xn_beam)
               endif
               jj0=(jsrc-1)*nf*nrec
               jj1=(jrec-1)*nf
               jsr=mzrec(jrec)
               jj2=jj0+jj1
               tfz(jf+jj2)=tf(jf+jj2) + cfac2*phisr(jsr)
cpln               tfz(jsrc,jrec)=tfz(jsrc,jrec)+cfac2*phisr(jsr)
            enddo
         enddo
99       continue
      enddo
c
c               write(6,*)'jrec,jsrc,jf: ',jrec,jsrc,jf
c               write(6,*)'jj0,jj1,jj2: ',jj0,jj1,jj2
c               tlz(1,jrec,jsrc)=tlz(1,jrec,jsrc) + cfac2*phisr(jsr)
c               write(6,*)tf(jf+jj2)
c               pause
      return
      end


