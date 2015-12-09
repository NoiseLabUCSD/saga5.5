      subroutine bb_field_reuse(knx,phiz,expz_gbs,jmk,tfz,jfbb)
c
c: Computes field tf for a single mode characterized by eigenvalue k and mode 
c: functions phi and dpsi.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
      integer*4 jmk,jfbb,jmo,jr,jsr,kk,kk0,jsrc,jrec,j,jzs,msrc,nrmax,
     .   jd,nsv,isv,jzsr,msrcmax,nsrcmax
      complex*8 phiz(nzsrgeom,jmk),tfz(nfbb,nrec,nrng),dphi_z
      complex*16 knx,sqkn,cfac2,h_arg,H0,iknx,phi_src
      real*8 expz_gbs(nzsrgeom,jmk),xn_beam,rsig_max,magsq
      real*8 phi_src_re,phi_src_im,dphi_src_re,dphi_src_im,
     .   dzsrgeom
      real*4 phiw_max,magsq_c8,ampl
c
c: Find mode function at source:
c: Normalize by density at source (9-1-94):
      jzs=1
      msrc=mzsrc(jzs)
      msrcmax=0
c: pln
c: hunt finds the index before zsrgeom becomes greater than zsr
      call hunt(zsrgeom,nzsrgeom,zsr(msrc),msrcmax)

      dzsrgeom = (zsr(msrc)-zsrgeom(msrcmax))
     .     /(zsrgeom(msrcmax+1)-zsrgeom(msrcmax))

      dphi_src_re=(real(phiz(msrcmax+1,jmk))
     .    -real(phiz(max(msrcmax,1),jmk)))*dzsrgeom
      dphi_src_im=(imag(phiz(msrcmax+1,jmk))
     .    -imag(phiz(max(msrcmax,1),jmk)))*dzsrgeom

      phi_src=dcmplx(dphi_src_re+real(phiz(msrcmax,jmk)),
     .     dphi_src_im+imag(phiz(msrcmax,jmk)))
     .     /rho_sr(msrcmax)

cpln      write(6,*)zsrgeom(msrcmax),zsrgeom(msrcmax+1)
cpln      write(6,*)zsr(msrc),phi_src

      xn_beam=(w/dreal(cp_sr(msrcmax)))*b_gbs(jzs) - 
     .     expz_gbs(msrcmax,jmk)
      iknx=dcmplx(-dimag(knx),dreal(knx))
      sqkn=cdsqrt(knx)
c
      do jrec=1,nrec
         jsr=mzrec(jrec)
         nsrcmax=0
         call hunt(zsrgeom,nzsrgeom,zsr(jsr),nsrcmax)
         dzsrgeom = (zsr(jsr)-zsrgeom(nsrcmax))
     .        /(zsrgeom(nsrcmax+1)-zsrgeom(nsrcmax))
         dphi_src_re=real(phiz(nsrcmax+1,jmk)
     .        -phiz(max(nsrcmax,1),jmk))*dzsrgeom
         dphi_src_im=imag(phiz(nsrcmax+1,jmk)
     .        -phiz(max(nsrcmax,1),jmk))*dzsrgeom
         phisr(jsr)=dcmplx(dphi_src_re+real(phiz(nsrcmax,jmk)),
     .        dphi_src_im+imag(phiz(nsrcmax,jmk)))
cpln         write(6,*)zsrgeom(nsrcmax+1),zsrgeom(nsrcmax)
cpln         write(6,*)jmk,zsr(jsr),phisr(jsr)
cpln         phisr(jsr)=phiz(jsr,jmk)
      enddo
c
c: Find range beyond which mode is insignificant (ranges have been sorted):
      rsig_max=kim_fac/dmax1(1.d-20,dimag(knx)-kim_bb(jfbb))
      nrmax=0
      call hunt(range,nrng,rsig_max,nrmax)
c
      do jr=1,nrmax
         kk0=krec_jr(jr)
         if(tilth) then
            kk=kk0+jr
            jsrc=jrec_jr(1,kk)
            jrec=jrec_jr(2,kk)
            h_arg=knx*range(jr)
            if(magsq(h_arg) .gt. 25.d0) then
czs: Include normalization by exp(-xn_beam) here:
               cfac2=sq2pir(jr,jrec)*phi_src*
     .              cdexp(iknx*range(jr) - xn_beam)/sqkn
            else
               call cdhankel(h_arg,1.d-6,H0)
               cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
            endif
            jsr=mzrec(jrec)
            tfz(jfbb,jrec,jsrc)=tfz(jfbb,jrec,jsrc) + 
     .           cfac2*phisr(jsr)
         else
            do kk=kk0+1,kk0+nrec_jr(jr)
               jsrc=jrec_jr(1,kk)
               jrec=jrec_jr(2,kk)
               h_arg=knx*(range(jr)+dtiltvp(jrec))
               if(magsq(h_arg) .gt. 25.d0) then
czs: Include normalization by exp(-xn_beam) here:
                  cfac2=sq2pir(jr,jrec)*phi_src*
     .                 cdexp(iknx*(range(jr)+dtiltvp(jrec))
     .                  - xn_beam)/sqkn
               else
                  call cdhankel(h_arg,1.d-6,H0)
                  cfac2=dcmplx(0.d0,pie)*phi_src*H0*dexp(-xn_beam)
               endif
               jsr=mzrec(jrec)
               tfz(jfbb,jrec,jsrc)=tfz(jfbb,jrec,jsrc) + 
     .         cfac2*phisr(jsr)
            enddo
         end if
      enddo
c
      return
      end
