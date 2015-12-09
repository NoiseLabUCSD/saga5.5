      subroutine cw_modes(iiwrt,ir,jrun,nrunx)
co	Modified for use as a subroutine with FGS
co	All verbose omitted.  Could still strip code
co	associated with unused options.  Note that deallocation
co	of memory is not done, see end of routine.

c
c: Performs CW mode computations
c
      use parms_com
      use i_o_com
      use gen_com
      implicit none
      integer*4 iiwrt,ir,jrun,nrunx,jfcw,jm,lsuf,lsuf0,
     .   lsuf2,jzs,iibad,j,jrec,jsrc,jf,ii_alloc
      real*8 vph,vg,db_km
      real*4 cpu1,cpu2,cpsec,time(2),pierad
      complex*8 pl0

      character*64 hdf_suf
      integer*4, dimension(:), allocatable:: indx
      real*4, dimension(:), allocatable:: fcw_sort
c      character*8 chnan
      data pierad/0.0174532925/

co	Make sure memory is allocated only once
co	See also deallocate statements at end of this routine
      if(iitl .ne. 0 .AND. jrun.NE.0 ) then
         allocate(tlc(iabs(nrec),iabs(nsrc)))
         allocate(tl(iabs(nrec),iabs(nsrc)))
         allocate(tli(iabs(nrec),iabs(nsrc)))
      endif
c
c: Create real*4 with NaN in it to initialize non-full HDF arrays:
c      chnan='7FC00000'
c      read(chnan,'(z8)') NaN
c
      iibad=0
      call uni_space(nfcw,fcw,1.e0)
      ii_alloc=0
      if(jrun .eq. 1) ii_alloc=1
      call sr_geom(ii_alloc)
c
co      if(iiwrt .eq. 1) then
co         open(7,file=outroot(1:lout)//'_modes',form='formatted')
co         if(iimt .ne. 0) then
co            open(21,file=outroot(1:lout)//'_mtraj',form='formatted')
co         endif
co      endif
c
      call zmx_init
c
       
		allocate(fcw_sort(nfcw))
		allocate(indx(nfcw))
c:		Sort in descending order:
		fcw_sort(1:nfcw)=-fcw(1:nfcw)
		call hpsort_r4_indx(nfcw,fcw_sort,indx)

c
co	Frequency loop
      do jf=1,nfcw
         jfcw=indx(jf)
         nctot=0
         nclast=0
         f_hz=fcw(jfcw)
co         call cpu_time(cpu1)


         if(iirx .le. 0) then
            call mode_find(iiwrt)
         else
            call rx_modes(iiwrt)
         endif


c         write(*,*)'number of modes',nmode
co         call cpu_time(cpu2)
co         cpsec=cpu2 - cpu1
co         if(iiwrt .eq. 1) then
co            write(6,220) 'Time taken to find modes = ',cpsec
co            write(2,220) 'Time taken to find modes = ',cpsec
co220         format(a,f7.2)
co         endif
c

Cdag280105         if(jf .eq. 1 .and. jrun .eq. 1) then
         if(jf .eq. 1 ) then
            nmode_hdf=nmode
c:			Take a conservative guess at the maximum # of modes at max CW frequency:
cc             nmode_hdf=nint((f_max/f_hz)*(nmode+2)) + 10
            do jm=1,nmode_hdf
               xmode(jm)=jm
            enddo
         else
            nmode=min(nmode_hdf,nmode)
         endif

         if(iikn .gt. 0 .and. iiwrt .eq. 1) then
            if(jf .eq. 1) then
               allocate(r4mat1(nfcw,nmode_hdf))
               allocate(r4mat2(nfcw,nmode_hdf))
               r4mat1(1:nfcw,1:nmode_hdf)=NaN
               r4mat2(1:nfcw,1:nmode_hdf)=NaN
            endif
            do jm=1,nmode
               r4mat1(jfcw,jm)=dreal(kn(jm))
               r4mat2(jfcw,jm)=dimag(kn(jm))
               if(iikn .gt. 1) r4mat2(jfcw,jm)=r4mat2(jfcw,jm)/kw0
            enddo
            do jm=nmode+1,nmode_hdf
               r4mat1(jfcw,jm)=NaN
               r4mat2(jfcw,jm)=NaN
            enddo
         endif
c
         if(iitl .ne. 0 .and. iiwrt .eq. 1) then

co		  Loop over sources.
            do jzs=1,nzs
               call mode_field(phi,dpsi,exp_gbs,tlc,tli,tl,jzs)
co               if(iitl .ne. 0) then
co                  call hdf_write_gen(4+ir,outroot,lout,zrec,nrec,t_src,
co     .               nsrc,zsrc,nzs,fcw,nfcw,var_ax,nrunx,tl,1,nrec,
co     .               1,nsrc,jzs,jzs,jfcw,jfcw,jrun,jrun,'Depth - m',9,
co     .               'Range - km',10,'Source Depth - m',16,
co     .               'Frequency - Hz',14,'Parameter',9,'TL - dB',
co     .               7,1,1,2)
co               endif
co               if(iabs(iitl) .eq. 3) then
co                  call hdf_write_gen(4+ir,outroot,lout,zrec,nrec,t_src,
co     .               nsrc,zsrc,nzs,fcw,nfcw,var_ax,nrunx,tli,1,nrec,
co     .               1,nsrc,jzs,jzs,jfcw,jfcw,jrun,jrun,'Depth - m',9,
co     .               'Range - km',10,'Source Depth - m',16,
co     .               'Frequency - Hz',14,'Parameter',9,'Inc TL - dB',
co     .               11,2,1,2)
co               endif
co               if(iitl .lt. 0) then
co                  do jrec=1,nrec
co                     do jsrc=1,nsrc
co                        pl0=tlc(jrec,jsrc)
co                        tli(jrec,jsrc)=atan2(aimag(pl0),real(pl0))/
co     .                     pierad
co                     enddo
co                  enddo
co                  call hdf_write_gen(4+ir,outroot,lout,zrec,nrec,t_src,
co     .               nsrc,zsrc,nzs,fcw,nfcw,var_ax,nrunx,tli,1,nrec,
co     .               1,nsrc,jzs,jzs,jfcw,jfcw,jrun,jrun,'Depth - m',9,
co     .               'Range - km',10,'Source Depth - m',16,
co     .               'Frequency - Hz',14,'Parameter',9,'TL Phase - deg',
co     .               14,3,1,2)
co               endif
            enddo
         endif

co		Transfer pressure field at this frequency to array
co		holding field at all frequencies.  Note that we are
co		'borrowing' use of an array actually intended for
co		use with bb computations
co		Field is summed over all sources
		tf(jfcw,1:nrec,1:nsrc)=tlc(1:nrec,1:nsrc)

co
co
         if(iisig .gt. 0 .and. iiwrt .eq. 1 .and. jf .eq. 1) then
c:		sigma_fill needs complex numbers for uz,ux,sigzz,sigzx:
            allocate(r4mat3(nzmf,4*nmode_hdf))
            allocate(r4mat4(nzmf,4*nmode_hdf))
            r4mat3(1:nzmf,1:4*nmode_hdf)=NaN
            r4mat4(1:nzmf,1:4*nmode_hdf)=NaN
         elseif(iimf .ne. 0 .and. iiwrt .eq. 1 .and. jf .eq. 1) then
c:		mfun_fill needs real numbers for re(phi),im(phi):
            allocate(r4mat3(nzmf,nmode_hdf))
            allocate(r4mat4(nzmf,nmode_hdf))
            r4mat3(1:nzmf,1:nmode_hdf)=NaN
            r4mat4(1:nzmf,1:nmode_hdf)=NaN
         endif
         if((iimf .eq. 1 .or. iimf .ge. 3) .and. iiwrt .eq. 1) then
            call mfun_fill(phi,r4mat3,r4mat4,nmode_hdf,ir,jrun,nrunx,
     .         12,iiri,iimp,jfcw,0)
            if(iidiag .eq. -2) call mode_ortho(phi,dphi,psi,dpsi)
         endif
         if((iimf .eq. 4) .and. iiwrt .eq. 1) then
c:		Always output Re(dphi), Im(dphi):
            call mfun_fill(dphi,r4mat3,r4mat4,nmode_hdf,ir,jrun,nrunx,
     .         28,3,0,jfcw,1)
         endif
         if((iimf .eq. 2 .or. iimf .eq. 3) .and. iiwrt .eq. 1) then
            call mfun_fill(psi,r4mat3,r4mat4,nmode_hdf,ir,jrun,nrunx,
     .         16,iiri,iimp,jfcw,2)
            if(iidiag .eq. -2) call mode_ortho(psi,dpsi,psi,dpsi)
         endif
c
         if(iisig .gt. 0) then
            call sigma_fill(phi,dphi,psi,dpsi,r4mat3,r4mat4,nmode_hdf,
     .         ir,jrun,nrunx,jfcw)
         endif
c
         if(iiwrt .eq. 1) then
co            write(7,202) f_hz,nmode,nctot,kw0
co            write(2,203) f_hz,nmode,nctot
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
co               if(abs(vph) .lt. 100000.d0) then
co                  write(7,201) jm,dreal(kn(jm))/kw0,
co     .               db_km,vph,vg,nzref(jm),ncalc(jm)
co               else
co                  write(7,204) jm,dreal(kn(jm))/kw0,
co     .               db_km,vph,vg,nzref(jm),ncalc(jm)
co               endif
co201			format(i5,2x,e14.8,4x,g11.5,1x,f13.6,1x,f11.5,2x,i5,2x,i5)
co204			format(i5,2x,e14.8,4x,g11.5,1x,e13.7,1x,f11.5,2x,i5,2x,i5)
            enddo
         endif
         if(iilist .ne. 0 .and. iiwrt .eq. 1) then
            call eig_list(jrun,nrunx,jf)
         endif
c
         if(iidc .ne. 0 .and. iiwrt .eq. 1) then
            if(jf .eq. 1) then
               allocate(r4mat5(nfcw,nmode_hdf))
               allocate(r4mat6(nfcw,nmode_hdf))
               r4mat5(1:nfcw,1:nmode_hdf)=NaN
               r4mat6(1:nfcw,1:nmode_hdf)=NaN
            endif
            do jm=1,nmode
               r4mat5(jfcw,jm)=dreal(eig_char(4,jm))
               r4mat6(jfcw,jm)=w/dreal(kn(jm))
            enddo
cc          call disp_fill(nmode,nfcw,jfcw,kn,eig_char,knbb,eig_bb,
cc   .         nm_cw_max,fcw,nmode_hdf)
         endif
      enddo ! jf

cpg      if(iiwrt .eq. 1) then
cpg         if(iidc .eq. 1 .or. iidc .eq. 3) then
cpg            call hdf_write_gen(2+ir,outroot,lout,fcw,nfcw,xmode,
cpg     .         nmode_hdf,var_ax,nrunx,1.,1,1.,1,r4mat5,1,nfcw,
cpg     .         1,nmode_hdf,jrun,jrun,1,1,1,1,'Frequency - Hz',14,
cpg     .         'Mode',4,'Parameter',9,' ',1,' ',1,'Group Speed - m/s',
cpg     .         17,24,1,2)
cpg            call hdf_write_gen(2+ir,outroot,lout,fcw,nfcw,xmode,
cpg     .         nmode_hdf,var_ax,nrunx,1.,1,1.,1,r4mat6,1,nfcw,1,
cpg     .         nmode_hdf,jrun,jrun,1,1,1,1,'Frequency - Hz',14,
cpg     .         'Mode',4,'Parameter',9,' ',1,' ',1,'Phase Speed - m/s',
cpg     .         17,25,1,2)
cpg         endif
cpg        if(iimt .eq. 1) close(21)
co         close(7)
c
cpg         if(iikn .gt. 0) then
cpg            call hdf_write_gen(2+ir,outroot,lout,fcw,nfcw,xmode,
cpg     .         nmode_hdf,var_ax,nrunx,1.,1,1.,1,r4mat1,1,nfcw,1,
cpg     .         nmode_hdf,jrun,jrun,1,1,1,1,'Frequency - Hz',14,'Mode',4,
cpg     .         'Parameter',9,' ',1,' ',1,'Re(kn)',6,10,1,2)
cpg            call hdf_write_gen(2+ir,outroot,lout,fcw,nfcw,xmode,
cpg     .         nmode_hdf,var_ax,nrunx,1.,1,1.,1,r4mat2,1,nfcw,1,
cpg     .         nmode_hdf,jrun,jrun,1,1,1,1,'Frequency - Hz',14,'Mode',4,
cpg     .         'Parameter',9,' ',1,' ',1,'Im(kn)',6,11,1,2)
cpg         endif
cpg      endif
c
      if(iiwrite .eq. 1) then
c:		Check if leaky modes computed in upper or lower isospeed halfspace:
         if(dreal(kn(nmode)) .lt. dreal(xkbp(1,1)) .and.
     .      zsr(nzsr) .gt. zdep(nlay-1) .and. isp(nlay) .eq. 1) 
     .      call print_warn('lower')
         if(dreal(kn(nmode)) .lt. dreal(xkbp(2,1)) .and.
     .      zsr(1) .lt. zdep(1)) call print_warn('upper')
      endif
c
co	Should deallocate at ultimate run only
co      if(iitl .ne. 0) then
co         deallocate(tlc)
co         deallocate(tli)
co         deallocate(tl)
co      endif
 
      if(iikn .ne. 0 .and. iiwrt .eq. 1) then
         deallocate(r4mat1)
         deallocate(r4mat2)
      endif
      if((iimf .ne. 0 .or. iisig .ne. 0) .and. iiwrt .eq. 1) then
         deallocate(r4mat3)
         deallocate(r4mat4)
      endif
      if(iidc .ne. 0 .and. iiwrt .eq. 1) then
         deallocate(r4mat5)
         deallocate(r4mat6)
      endif


	deallocate(fcw_sort)
	deallocate(indx)
c
      return
      end
ccc
      subroutine print_warn(ch_hsp)
c
      implicit none
      character*5 ch_hsp
c
      print *,'WARNING: Field or mode functions requested in ',
     .   'lower halfspace and leaky mode(s) found.'
      write(2,*) 'WARNING: Field or mode functions requested in ',
     .   ch_hsp,' halfspace and leaky mode(s) found.'
      print *,'Be aware that leaky modes are not valid in ',
     .   'isospeed halfspaces!!'
      write(2,*) 'Be aware that leaky modes are not valid in ',
     .   'isospeed halfspaces!!'
      print *,'One solution is to insert a false layer.'
      write(2,*) 'One solution is to insert a false layer.'
c
      return
      end
