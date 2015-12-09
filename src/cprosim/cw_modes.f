      subroutine cw_modes(iiwrt)
c
c: Performs CW mode computations
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
      integer*4 iiwrt,jfcw,jm,lsuf,lsuf0,lsuf2,jzs,
     .   jj,iibad,j
      real*8 vph,vg,db_km
      real*4 cpsec0,cpsec,etime,time(2),pierad
      complex*8 pl0
      character*64 hdf_suf
      data pierad/0.0174532925/
c
      nctot=0
      nclast=0
      nmode=0
      iibad=0
      call uni_space(nfcw,fcw,1.e0)
      call sr_geom(rng_sr,iabs(nsrc),iabs(nrec))
      if(iitl .ne. 0) then
         call mem_lim(nrec*nsrc,NTLMAX,MLINE,LML,'nrec*nsrc',
     .      9,'NTLMAX',7,iibad,1)
      endif
c
      if(iiwrt .eq. 1) then
         open(7,file=outroot(1:lout)//'_modes',status='unknown',
     .        form='formatted')
         if(iimt .ne. 0) then
            open(21,file=outroot(1:lout)//'_mtraj',status='unknown',
     .        form='formatted')
         endif
      endif
c
      call zmx_init
c
      do jfcw=1,nfcw
         call suffix(' ',0,jfcw,nfcw,SUFX//'f',2,hdf_suf,lsuf)
         lsuf0=lsuf
         nctot=0
         nclast=0
         f_hz=fcw(jfcw)
         cpsec0=etime(time)
         if(iirx .le. 0) then
            call mode_find(iiwrt)
            if(jjfail.gt.0) return
         else
c            call rx_modes
            write(6,*)'No real-axis version'
         endif
         cpsec=etime(time) - cpsec0
         if(iiwrt .eq. 1) then
            write(6,220) 'Time taken to find modes = ',cpsec
            write(lusvp,220) 'Time taken to find modes = ',cpsec
220         format(a,f7.2)
         endif
c
         if(iitl .ne. 0 .and. iiwrt .eq. 1) then
            do jzs=1,nzs
               call suffix(hdf_suf,lsuf0,jzs,nzs,SUFX//'zs',3,hdf_suf,
     .            lsuf)
ccjj           jj=(jzs-1)*nrec*nsrc + 1
               jj=1
               call mode_field(phi,dpsi,exp_gbs,plc(jj),tl(jj),jzs)
               lsuf2=lsuf + 3
               hdf_suf(lsuf+1:lsuf2)=SUFX//'tl'
               call out_writex(outroot,lout,hdf_suf,lsuf2,tl(jj),
     .            rec_lab,t_src,nrec,nsrc,dlab,rlab,z4,z4,z4,z4,2,
     .            pllab,'m','km',dblab,'f7.2','f9.3','f7.3',ncall)
               if(iitl .lt. 0) then
                  do j=1,nrec*nsrc
                     pl0=plc(jj+j-1)
                     r4mat1(j)=atan2(aimag(pl0),real(pl0))/pierad
                  enddo
                  hdf_suf(lsuf+1:lsuf2)=SUFX//'ph'
                  call out_writex(outroot,lout,hdf_suf,lsuf2,r4mat1,
     .               rec_lab,t_src,nrec,nsrc,dlab,rlab,z4,z4,z4,z4,2,
     .               'phase','m','km','deg','f7.2','f9.3','f7.3',
     .               ncall)
               endif
            enddo
         endif

         if(iifft .ne. 0) then
            do jzs=1,nzs
               jj=1
               call mode_cmplx_field(tf(jj),phi,dpsi,exp_gbs,jzs,
     .                               jfcw,nfcw)
            enddo
         endif

         if((iimf .eq. 1 .or. iimf .eq. 3) .and. iiwrt .eq. 1) then
            lsuf=lsuf0 + 4
            hdf_suf(lsuf0+1:lsuf)=SUFX//'phi'
            call mfun_fill(phi,r4mat1,r4mat2,hdf_suf,lsuf)
            if(iidiag .eq. -2) call mode_ortho(phi,dphi,psi,dpsi)
c: Do this if you want to output d(phi)/dz also:
cxx         lsuf=lsuf0 + 5
cxx         hdf_suf(lsuf0+1:lsuf)=SUFX//'dphi'
cxx         call mfun_fill(dphi,r4mat1,r4mat2,hdf_suf,lsuf)
         endif
         if((iimf .eq. 2 .or. iimf .eq. 3) .and. iiwrt .eq. 1) then
            lsuf=lsuf0 + 4
            hdf_suf(lsuf0+1:lsuf)=SUFX//'psi'
            call mfun_fill(psi,r4mat1,r4mat2,hdf_suf,lsuf)
            if(iidiag .eq. -2) call mode_ortho(psi,dpsi,psi,dpsi)
         endif
c
c
         if(iiwrt .eq. 1) then
            write(7,202) f_hz,nmode,nctot,kw0
            write(lusvp,203) f_hz,nmode,nctot
202         format('f = ',f8.2,'; # of Modes =',i5,
     .         '; # R1*R2 Calcs =',i6,'; Kw = ',e14.8/
     .         'Mode#        Re(k)/Kw    Im(k)-dB/km     Phase Vel   ',
     .         'Group Vel  Duct#  #Calc')
203         format('f = ',f8.2,'; # of Modes =',i5,
     .         '; # R1*R2 Calcs =',i6)
            do jm=1,nmode
               vph=w/dreal(kn(jm))
               vg=dreal(eig_char(4,jm))
               db_km=8685.889638*dimag(kn(jm))
               if(abs(vph) .lt. 100000.d0) then
                  write(7,201) jm,dreal(kn(jm))/kw0,
     .               db_km,vph,vg,nzref(jm),ncalc(jm)
               else
                  write(7,204) jm,dreal(kn(jm))/kw0,
     .               db_km,vph,vg,nzref(jm),ncalc(jm)
               endif
201   format(i5,2x,e14.8,4x,g11.5,1x,f13.6,1x,f11.5,2x,i5,2x,i5)
204   format(i5,2x,e14.8,4x,g11.5,1x,e13.7,1x,f11.5,2x,i5,2x,i5)
            enddo
         endif
      end do
c
      if(iiwrite .eq. 1) then
c: Check if leaky modes computed in upper or lower isospeed halfspace:
         if(dreal(kn(nmode)) .lt. dreal(xkbp(1,1)) .and.
     .      zsr(nzsr) .gt. zdep(nlay-1) .and. isp(nlay) .eq. 1) 
     .      call print_warn('lower')
         if(dreal(kn(nmode)) .lt. dreal(xkbp(2,1)) .and.
     .      zsr(1) .lt. zdep(1)) call print_warn('upper')
      endif
c
      nfbb=nfcw
c pln 020500      call bb_fft_out(cfmin)
c
      return
      end
ccc
      subroutine print_warn(ch_hsp)
c
      implicit none
      include 'Parms_com'
      character*5 ch_hsp
c
      print *,'WARNING: Field or mode functions requested in ',
     .   'lower halfspace and leaky mode(s) found.'
      write(lusvp,*) 'WARNING: Field or mode functions requested in ',
     .   ch_hsp,' halfspace and leaky mode(s) found.'
      print *,'Be aware that leaky modes are not valid in ',
     .   'isospeed halfspaces!!'
      write(lusvp,*) 'Be aware that leaky modes are not valid in ',
     .   'isospeed halfspaces!!'
      print *,'One solution is to insert a false layer.'
      write(lusvp,*) 'One solution is to insert a false layer.'
c
      return
      end
