c:**********************************************
c:*   AUTHOR:                                  *
c:*      Evan Westwood                         *
c:*      Applied Research Laboratories         *
c:*      The University of Texas at Austin     *
c:*      P. O. Box 8029                        *
c:*      Austin, TX  78713-8029                *
c:**********************************************
c: *******************************
c: *     REVISION (1996):        *
c: *         E M G  GROUP        *
c: *     S A C L A N T C E N     *
c: *******************************
                                                 
      subroutine opt_read_pro
c
c: Reads option file.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
c
      integer*4 nline,iierr
      integer*4 IOP

      character*64 eline

c: Local variables:
      integer*4 jj

      data eline/'INVALID INPUT IN OPT FILE: '/
c
      iierr=0
      nline=0
      svp_ver= 2.0
      rho_svp= 1.0
      alpha_svp= 0.0
c
      nline=0
      iierr=0
      iigeom=1
      ver_no= 2.0
      IOP= 16
      iikpl= 0
      iirc= 0
      n_env= 0
c      iirx= -1
      cphmin= 0.0
      cphmax= 0.0
      rmin=0.0
      rmax= 0.0
      phfac= 4
      db_cut= 48
      iidiag= 0
      lout=11
      iiaih(1)=0.0
      iiaih(2)=0.0
      iigbs=0
c
cfmc  iifft   = Output FFT File (0=no;1=zs,zr,r on Line 9; 2=read file on Line 11);
      iifft= 1
cfmc  iiout   = Output BB Eigenvalues and Functions (same options as iifft above);
      iiout= 0
cfmc   iift    = Freq Traj(ASCII); iimt = Mode Traj(ASCII); 
      iift= 0
cfmc    iimt= ????
      iimt= 0
cfmc  iidc    = Disp Curves (0=no,1=vg,2=vph,3=both); 
      iidc= 0

      call check_val_r4(fsbb,0.e0,1.e10,eline,27,nline,
     .         'fsbb',4,iierr)
      call check_val_r4(Tw,-1.e10,131072e0,eline,27,nline,
     .         'nfft/Tw',7,iierr)
      call check_val_r4(fmin,1.e-3,fmax,eline,27,nline,
     .         'fmin',4,iierr)
      call check_val_r4(fmax,fmin,fsbb/2.e0,eline,27,nline,
     .         'fmax',4,iierr)
c
      if(iierr .eq. 1) then
         print *,' '
         print *,'Execution terminating.  Check input option file '//
     .      'for error(s).'
         stop
      endif
      return
      end
ccc
      subroutine check_val_i(val,val_lo,val_hi,eline,le,nline,
     .   vname,lv,iierr)
c
c: Checks integer input val to see if it is in the range of allowable 
c: values, val_lo to val_hi.
c
      implicit none
      integer*4 val,val_lo,val_hi,nline,le,lv,iierr
      character*64 eline,vname
c
      if(val .lt. val_lo .or. val .gt. val_hi) then
         iierr=1
         print *,' '
         print *,eline(1:le),' LINE # ',nline
         print *,'VARIABLE NAME = ',vname(1:lv),'; VALID RANGE = ',
     .      val_lo,' ,',val_hi
         print *,'ENTERED VALUE = ',val
      endif
c
      return
      end
ccc
      subroutine check_val_r4(val,val_lo,val_hi,eline,le,nline,
     .   vname,lv,iierr)
c
c: Checks real*8 input val to see if it is in the range of allowable 
c: values, val_lo to val_hi.
c
      implicit none
      real*4 val,val_lo,val_hi
      integer*4 nline,le,lv,iierr
      character*64 eline,vname
c
      if(val .lt. val_lo .or. val .gt. val_hi) then
         iierr=1
         print *,' '
         print *,eline(1:le),' LINE # ',nline
         print *,'VARIABLE NAME = ',vname(1:lv),'; VALID RANGE = ',
     .      val_lo,' ,',val_hi
         print *,'ENTERED VALUE = ',val
      endif
c
      return
      end
ccc
      subroutine check_val_r8(val,val_lo,val_hi,eline,le,nline,
     .   vname,lv,iierr)
c
c: Checks real*8 input val to see if it is in the range of allowable 
c: values, val_lo to val_hi.
c
      implicit none
      real*8 val,val_lo,val_hi
      integer*4 nline,le,lv,iierr
      character*64 eline,vname
c
      if(val .lt. val_lo .or. val .gt. val_hi) then
         iierr=1
         print *,' '
         print *,eline(1:le),' LINE # ',nline
         print *,'VARIABLE NAME = ',vname(1:lv),'; VALID RANGE = ',
     .      val_lo,' ,',val_hi
         print *,'ENTERED VALUE = ',val
      endif
c
      return
      end
ccc
