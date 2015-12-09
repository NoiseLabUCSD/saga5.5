      program cprosim
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
c:	   model for acousto-elastic ocean environments," submitted to 
c:	   J. Acoust. Soc. Am., Jan 1996.
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
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
      include 'lab_com'
      include 'debug_com'
c
      external time
      integer*4 time,isec1,isec2,no_sec,no_sec10,kdat(3),ktim(3),j
      real*4 cpsec,cptim,etime
cpln      external hand
      integer ieee_handler,intrp(RGMAX),isectr,isect,isubb
c
      call vbl_init
      ver_cur=2.0
c
c: For SGI:
cxx   call idate(kdat(2),kdat(1),kdat(3))
c: For Sun:
      call idate(kdat(1),kdat(2),kdat(3))
      kdat(3)=mod(kdat(3),100)
      call itime(ktim)


      lout= 11
      outroot(1:lout)= 'cprosim_out'
      outfile=outroot(1:lout)//'.dat'
      loutf=lout + 4
cpln      open(2,file=outfile(1:loutf),status='unknown',
cpln     .  form='formatted')
cpln      write(lusvp,195) kdat(2),kdat(1),kdat(3),(ktim(j),j=1,3),
cpln     .   svp_file(1:lsvp),opt_file(1:lopt),outroot(1:lout)
cpln195   format('ORCA NORMAL MODE MODEL RUN: DATE = ',2(i2.2,'/'),i2,
cpln     .   '; TIME = ',i2.2,2(':',i2.2)/'SVP FILE = ',a
cpln     .   /'OPT FILE = ',a/'OUTPUT FILE ROOT = ',a)
c
      print *,'Running CPROSIM v2.0b ...'
      call flush(6)
      ncall=0
      ncmat=0
      isec1=time()
c
c
c: Set iiwrite=1 so that output files and messages are sent:
      iiwrite=0
      iidebug=0
c: Read input svp file in file name svp_file(1:lsvp):
c      call svp_read
      call svp_read_pro
c: Read input option file in file name opt_file(1:lopt):
      call opt_read_pro
c: Read input option file in file name opt_file(1:lopt):
c
      isect=1
      isubb=1
      call tenv_read_pro(isect,intrp,isubb)
c
c Set constants from input run stream to force complex plane ORCA
      iirx=0
      iikpl=0
      iirc=0
      iiparm=0
      n_env=0
      iikrak=0
      iifepe=0
      iioas=0
      iimlab=0
      iifmt=0
c
c: iidiag=-5 means trap IEEE exceptions in subroutine hand:
c$$$      if(iidiag .eq. -5) j=ieee_handler( 'set', 'common', hand )
      if(iidiag .eq. -5) then
         write(6,*)'Should have called ieee_handler'
      end if
c
c: Check input parameters and set up master files:
      call svp_check(1)
      call svp_check2
c
      if(iiparm .eq. 0) then
         if(iicw .eq. 1) then
c: CW mode computations:
            call cw_modes(1)
         elseif(iicw .eq. 2) then
c: BB mode computations:
            call bb_modes
         endif
         if(iikpl .ne. 0) then
c: Complex k-plane computations:
            call k_plane(nreal,nimag,r4mat1,r4mat2,tl,plc)
         endif
      endif
c
c: Check if failure occurred in mode_find or eig_find:
      if((iifail .eq. 1) .or. (jjfail .gt. 0)) then
         print *,'FAILURE OCCURRED FOR THIS RUN.'
         print *,'IIFAIL, JJFAIL: ',iifail,jjfail
      endif
c
      if((iifmt .eq. 2 .or. iifmt .eq. 0) .and. ncmat .gt. 0) then
c: Close Matlab file:
         j=-1
         call out_writex(outroot,lout,' ',0,r4mat1,r4mat1,
     .      r4mat1,0,0,' ',' ',0.,0.,0.,0.,2,' ',
     .      ' ',' ',' ',' ',' ',' ',j)
      endif
c
      isec2=time() - isec1
      cptim=etime(cpsec)
      no_sec=int(mod(cptim,60.e0))
      no_sec10=nint(100.*(cptim - int(cptim)))
      write(lusvp,330) int(cptim/60.e0),no_sec,no_sec10,
     .   int(isec2/60),mod(isec2,60)
330   format('TOTAL CP MINUTES TAKEN FOR RUN = ',i4,':',i2.2,'.',i2.2,
     .   '; WALL CLOCK TIME = ',i4,':',i2.2)
      close(2)
      print *,'ORCA v1.1 finished. Listing on file ',outfile(1:loutf)
cpln      call system('more '//outfile(1:loutf))
cpln      if(iicw .eq. 1) call system('more '//outroot(1:lout)//'_modes')
c
      stop
      end
ccc
cpln      function hand(sig,sip,uap)
cpln      integer sig,address
cpln      structure /fault/
cpln           integer address
cpln      end structure
cpln
cpln      structure /siginfo/
cpln           integer si_signo
cpln           integer si_code
cpln           integer si_errno
cpln           record /fault/ fault
cpln      end structure
cpln      record /siginfo/ sip
cpln      address = sip.fault.address
cpln      write (*,10) address
cpln10    format('Exception at hex address ', z8 )
cpln      return
cpln      end

