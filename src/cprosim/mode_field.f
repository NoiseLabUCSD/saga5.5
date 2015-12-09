      subroutine mode_field(phiz,dpsiz,expz_gbs,plcx,tlx,jzs)
c
c: Computes field given mode eigenvalues and mode function values.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jzs,jm,jr,jd,jsr,kk,kk0,jsrc,jrec,msrc,nrmax
      complex*16 knx,iknx,iknr,sqkn,cfac2,phi_src,H0,h_arg,
     .   cfac,exp_sqk
      complex*8 phiz(nzsr,nmode),dpsiz(nzsr,nmode),plcx(nsrc,nrec),tlz
      real*8 expz_gbs(nzsr,nmode),rsig_max,hmsq,xn_b,xn_beam,
     .   magsq,rng_im
      real*4 tlx(nsrc,nrec),tlmagsq
c: cfac = sqrt(2*pi)*exp(i*pi/4) (see p. 7, ORCA II):
      data cfac/(1.77245385090552d0,1.77245385090552d0)/
c
      msrc=mzsrc(jzs)
      if(iigbs .ne. 0) then
         rng_im=-b_gbs(jzs)*cos(th_gbs(jzs)*pie/180.d0)
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .         +dtiltvp(jd),rng_im+dtiltvp(jd)))
            end do
cpln            sq2pir(jr)=cfac/cdsqrt(dcmplx(range(jr),rng_im))
         enddo
      else
         rng_im=0.d0
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .         +dtiltvp(jd),rng_im+dtiltvp(jd)))
            end do
cpln            sq2pir(jr)=cfac/dsqrt(range(jr))
         enddo
      endif
c
      hmsq=25.d0
c: Initialize complex propagation loss array:
      do jsrc=1,nsrc
         do jrec=1,nrec
            plcx(jsrc,jrec)=cmplx(0.e0,0.e0)
         enddo
      enddo
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
            kk0=krec_jr(jr)
            do kk=kk0+1,kk0+nrec_jr(jr)
               jsrc=jrec_jr(1,kk)
               jrec=jrec_jr(2,kk)
               h_arg=knx*dcmplx(range(jr)+dtiltvp(jrec),
     .               rng_im+dtiltvp(jrec))
cpln               h_arg=knx*dcmplx(range(jr),rng_im)
               if(magsq(h_arg) .gt. hmsq) then
c     zs: Include normalization by exp(-xn_beam) here:
                  iknr=dcmplx(-dimag(h_arg),dreal(h_arg))
                  exp_sqk=cdexp(dcmplx(dreal(iknr)-xn_beam,dimag(iknr)))
     .                 /sqkn
                  cfac2=sq2pir(jr,jrec)*phi_src*exp_sqk
cpln               cfac2=sq2pir(jr)*phi_src*exp_sqk
               else
                  call cdhankel(h_arg,1.d-6,H0)
c     c             cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
                  cfac2=dcmplx(-pie*dimag(H0),pie*dreal(H0))*phi_src*
     .                 dexp(-xn_beam)
               endif
               jsr=mzrec(jrec)
               plcx(jsrc,jrec)=plcx(jsrc,jrec) + cfac2*phisr(jsr)
            enddo
         enddo
99       continue
      enddo
      do jsrc=1,nsrc
         do jrec=1,nrec
            tlz=plcx(jsrc,jrec)
            tlmagsq=real(tlz)*real(tlz)+aimag(tlz)*aimag(tlz)
            tlx(jsrc,jrec)=10.*alog10(amax1(1.e-37,tlmagsq))
         enddo
      enddo
c
      return
      end
