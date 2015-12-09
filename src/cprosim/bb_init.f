      subroutine bb_init
c
c: Checks and initializes variables for broadband mode calculations.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 n,nfft0,j,iibad,ii1st,jr,jzs,jd
      real*8 rng_im
      complex*16 cfac
      data cfac/(1.77245385090552d0,1.77245385090552d0)/
c      data ii1st/1/

      ii1st=1
c
      if(ii1st .eq. 0) goto 10
c
c: Initialize source and receiver geometry arrays:
      call sr_geom(rng_sr,iabs(nsrc),iabs(nrec))
      if(iigbs .eq. 0) then
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/dsqrt(range(jr)+dtiltvp(jd))
            end do
cpln            sq2pir(jr)=cfac/dsqrt(range(jr))
         enddo
      else
         jzs=1
         rng_im=-b_gbs(jzs)*cos(th_gbs(jzs)*pie/180.d0)
         do jr=1,nrng
            do jd=1,nrec
               sq2pir(jr,jd)=cfac/cdsqrt(dcmplx(range(jr)
     .              +dtiltvp(jd),rng_im+dtiltvp(jd)))
            end do
cpln            sq2pir(jr)=cfac/cdsqrt(dcmplx(range(jr),rng_im))
         enddo
      endif
10    continue
c
      if(iicw .eq. 2) then
         if(iiwrite.ge.1) then
            write(6,*)'Delta frequency= ',df
            write(6,*)'Maximum time window= ',1./df
c     fmc
cpln            print *,' from subroutine bb_init.f :'
            print *,' fs,   nfft, Tw ',  fsbb, nfftbb, Tw
            print *,' fmin, fmax, df :', fmin, fmax, df
            print *,' nf1,  nf2,  nfreq :', nf1, nf2, nfbb
c     fmc
            iibad=0
            call mem_lim(nfbb,NFBBMAX,MLINE,LML,
     .           'nfbb',4,'NFBBMAX',7,iibad,0)
            call mem_lim(nfbb*nrec*nsrc,NTFMAX,
     .           MLINE,LML,'nfbb*nrec*nsrc',14,
     .           'NTFMAX',6,iibad,0)
            if(iibad .eq. 1) then
              write(*,*)'   nfbb,NFBBMAX,nfbb,nrec,nsrc'           
              write(*,*)   nfbb,NFBBMAX,nfbb,nrec,nsrc          
               stop 'from bb_init'
c     
            endif
         endif
      endif
c
      if(i_geom .ne. 0) then
         do j=1,nfbb
            nmbb(j)=0
         enddo
      end if
c
      do j=1,nfbb*nrec*nsrc
         tf(j)=(0.,0.)
      enddo
      do j=1,nfbb
         faxbb(j)=fmin + (j-1)*df
         wbb(j)=twpie*faxbb(j)
         kim_bb(j)=1.d100
         phim_bb(j)=0.
      enddo
c
      if(iimt .ne. 0) then
         open(21,file=outroot(1:lout)//'_mtraj',status='unknown',
     .        form='formatted')
      endif
c
      ii1st=0
c
      return
      end
ccc
      subroutine bb_out_init(nm_fmax)
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c
      integer*4 nm_fmax
c
c: Open direct access for broadband mode characteristics, if desired:
      if(iirx .eq. -1) then
         lrecw=(4+4*nzsr)*nfbb
      elseif(iirx .eq. 0) then
c: (Complex*16 eigenvalue + complex*8 mode function at nzsr depth) at
c: each frequency bin:
         lrecw=4*nm_fmax + 2*nzsr*nm_fmax
      else
c: (Complex*16 eigenvalue + real*4 mode function at nzsr depth) at
c: each frequency bin:
         lrecw=(4+nzsr)*nm_fmax
      endif
      open(33,file=outroot(1:lout)//'_bbeig',status='unknown',
     .   access='direct',recl=NRECL*lrecw)
      lheadw=11 + 3*nzsr + nfbb
      nh_off=(lheadw - 1)/lrecw + 1
c
      return
      end
