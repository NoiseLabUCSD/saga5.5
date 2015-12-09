      subroutine orca( ocount )

c	ORCA in subroutine version (unauthorized)
c	lines not required commented out by 'co' or removed

c
c: Normal mode model for acoustic propagation in range independent ocean
c: environments with fluid and/or solid layered structures above and below.  
c: Computes leaky modes and seismic modes.   Users must construct an SVP
c: file (see the template file orca_svp) that specifies the ocean
c: environment and an OPTION file (see the template file orca_opt) that 
c: specifies the run parameters and outputs desired.
c:
c: Questions or comments should be referred to:
c:
c: Evan Westwood
c: Applied Research Laboratories
c: The University of Texas at Austin
c: P. O. Box 8029
c: Austin, TX  78713-8029
c
c:   References for the complex k-plane model at this time (05-21-96) are:
c:	E. K. Westwood, C. T. Tindle, and N. R. Chapman, "A normal mode
c:         model for acousto-elastic ocean environments," J. Acoust. Soc. Am.,
c:         100, 3631-3645 (1996).
c:	E. K. Westwood, C. T. Tindle, and N. R. Chapman, "A normal mode 
c:	   model for multilayered acoustoelastic ocean environments based on 
c:	   an analytic reflection coefficient method," J. Acoust. Soc. Am., 
c:	   95, No. 5, Pt. 2, 2908 (1994).
c:	E. K. Westwood, "An efficient broadband normal-mode model for 
c:	   acoustoelastic ocean environments," J. Acoust. Soc. Am., 96, 
c:	   No. 5, Pt. 2, 3352 (1994).
c:   References for the real k-axis option are:
c:	E. K. Westwood, "Improvements to narrow-band and broadband normal-mode
c:	   algorithms for fluid ocean environments," J. Acoust. Soc. Am.,
c:	   99, No. 4, Pt. 2, 2524 (1996).
c:	S. J. Levinson, E. K. Westwood, R. A. Koch, S. K. Mitchell, and 
c:	   C. V. Sheppard, "An efficient and robust method for underwater 
c:	   acoustic normal-mode computations," J. Acoust. Soc. Am., 97, 
c:	   1576-1585 (1995).
c
      use parms_com
      use i_o_com
      use gen_com
      implicit none
      include 'sd_com'

c	This variable is input
c	at ocount==1 storage is allocated 
	integer ocount

      integer*4 values(8),j
      real*4 cpu1,cpu2,nsec,cptim
      character*10 jday,jtime,jzone
c      external time,etime,clock,jdate,date_and_time
c
co      open(22,file='orca_status',form='formatted')
co      write(22,*) 'Starting ORCA ...'
co      close(22)
      iifirst(1:12)=1   
      call vbl_init
      ver_cur=2.0
co      open(44,file='orca_in',form='formatted',err=100)
co      read(44,440,end=100,err=120) svp_file
co      read(44,440,end=100,err=120) opt_file
co      read(44,440,end=100,err=120) outroot
440   format(a)
co      close(44)
co      call lname64(svp_file,lsvp)
co      call lname64(opt_file,lopt)
co      call lname64(outroot,lout)
c
c: For SGI:
cxx   call idate(kdat(2),kdat(1),kdat(3))
c: For Sun:
co      call cpu_time(cpu1)
co      call date_and_time(jday,jtime,jzone,values)
      outfile=outroot(1:lout)//'_out'
      loutf=lout + 4
co      open(2,file=outfile(1:loutf),form='formatted')
co      write(2,195) values(2),values(3),values(1),values(5),
co     .   values(6),values(7),values(8),
co     .   svp_file(1:lsvp),opt_file(1:lopt),outroot(1:lout)
195   format('ORCA NORMAL MODE MODEL RUN: DATE = ',2(i2.2,'/'),i4,
     .   '; TIME = ',i2.2,3(':',i2.2)/'SVP FILE = ',a
     .   /'OPT FILE = ',a/'OUTPUT FILE ROOT = ',a)
c
      print *,'Running ORCA v2.0 ...'
      ncall=0
      ncmat=0
c
c: Set iiwrite=1 so that output files and messages are sent:
co      iiwrite=1
	iiwrite=0
c: Read input svp file in file name svp_file(1:lsvp):
co	Should always be called for sanity checks on perturbed environment
      call svp_read
c: Read input option file in file name opt_file(1:lopt):
co	Need be called only once unless options are changed
	if( ocount == 1 ) then
	      call opt_read
	endif
c: iidiag=-5 means trap IEEE exceptions in subroutine hand:
cc    if(iidiag .eq. -5) j=ieee_handler( 'set', 'common', hand )
c
c: Check input parameters and set up master files:
      call svp_check(1)
      call svp_check2
c
co      if(n_env .gt. 0) then
c:		Write out environmental parameters:
co         allocate(r4mat1(2*(n_env+20),0:5))
co         call env_write(r4mat1)
co         deallocate(r4mat1)
co      endif
c
      if(iiparm .eq. 0) then
         allocate(var_ax(1))
         if(iicw .eq. 1) then
c: CW mode computations:
co		prevent new allocation of storage for repeated runs.
co		Must set first argument (iiwrt=1) to get mode fields
co		Note that memory is never deallocated...
			if( ocount == 1 ) then
co				borrowing use of an array otherwise
co				allocated in bb_init
				allocate(tf(iabs(nfcw),iabs(nrec),iabs(nsrc)))
				call cw_modes(1,0,1,1)
			else
				call cw_modes(1,0,0,1)
			endif
         elseif(iicw .eq. 2) then
c: BB mode computations:
            if(iirx .eq. 0) then
c: CW loop BB is safest and not always much slower, so make default
               call bb_brute
            elseif(iirx .eq. 1) then
               call rx_bb
c: bb_modes is too complicated to maintain - always do simpler brute force
c: in complex k-plane
cc          elseif(iirx .eq. -1) then
cc             call bb_modes
            elseif(iirx .eq. 2) then
               call rx_bb_brute
            endif
         endif
         if(iikpl .ne. 0) then
c: Complex k-plane computations:
            call k_plane(0,1,1)
         endif
         if(iirc .ne. 0) then
c: Plane wave reflection coefficient computations:
            call pw_refco(0,1,1)
         endif
         deallocate(var_ax)
      elseif(iiparm .eq. 1) then
c: Variation of environmental parameters:
         call vary_env
      endif
c
      if(iikrak .gt. 0) call krak_out
      if(iifepe .gt. 0) then
         f_hz=fcw(1)
         call fepe_out
      endif
      if(iioas .gt. 0) call oases_out
      if(iimlab .gt. 0) call modelab_out
c
c: Check if failure occurred in mode_find or eig_find:
      if(iifail .eq. 1) then
         print *,'FAILURE OCCURRED FOR THIS RUN.'
      endif
c
ctms
co      call hdf_write_end
c
co      call cpu_time(cpu2)
co      cptim=cpu2-cpu1
c      print *,'cpu_time = ',cpu1,cpu2,cptim
co      nsec=mod(cptim,60.)
co      write(2,330) int(cptim/60.e0),int(nsec),nint(100.*mod(nsec,1.))
330   format('TOTAL CP MINUTES TAKEN FOR RUN = ',i4.2,':',i2.2,'.',i2.2)
co      close(2)
co      print *,' '
co      print *,'ORCA finished. Output listing in file ',
co     .   outfile(1:loutf),':'
c      print *,' '
c      call system('more -10000 '//outfile(1:loutf))
      if(iicw .eq. 1) then
co         print *,' '
co         print *,'Mode list in file ',outroot(1:lout),'_modes.'
co         print *,' '
ctms     if(nmode .lt. 50) then
ctms         call system('more -10000 '//outroot(1:lout)//'_modes')
ctms     endif
      endif
      if(iifail .eq. 0) then
co         open(22,file='orca_status')
co         write(22,*) 'ORCA ended successfully'
co         close(22)
      endif
c
      write(6,*) 'LET ME OUT!'
      RETURN

100   print *,'Error opening or reading input file orca_in.'
      print *,'Place _svp, _opt, and output file names in file '//
     .   'named orca_in.'
      stop
120   print *,'End encoutered in input file orca_in.'
      print *,'Place _svp, _opt, and output file names in file '//
     .   'named orca_in.'
      stop
      end

