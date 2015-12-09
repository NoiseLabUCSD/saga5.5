      subroutine out_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
     .   xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt,ncall)
c
c: This subroutine chooses the subroutine to call in order to output 
c: the array cmat and the axes x and y.
c
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
c      include 'i_o_com'
      integer*4 ncall,lout0,lsuf,nx,ny,nf,ncasc,mat_fid,iilong
      character*64 outroot0,ascsuf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt
      character*30 xlab,ylab
      real*4 cmat(ny,nx),x(nx),y(ny),xmin,xmax,ymin,ymax
      data iilong/0/
c
      if(iiwrite .eq. 0) then
         return
      endif
c: ncall<0 means to close Matlab file:
      if(ncall .lt. 0) then
cpln         stat = matClose(mat_fid)
         if(iilong .gt. 0) then
            print *,'Warning: Variable names in .mat file '//
     .         'shortened'
            write(nf,*) 'Warning: Variable names in .mat file '//
     .         'shortened'
         endif
         return
      endif
c
      if(iifmt .eq. 1 .or. iifmt .eq. 0) then
cpln         call hdf_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
cpln     .      xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
cpln     .      xfmt,yfmt,dfmt,ncall)
      endif
      if(iifmt .eq. 2 .or. iifmt .eq. 0) then
cpln         call mat_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
cpln     .      xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
cpln     .      xfmt,yfmt,dfmt,ncmat,mat_fid,iilong)
      endif
      if(iifmt .eq. 3 .or. iifmt .eq. 0) then
         call asc_writex(outroot0,lout0,ascsuf,lsuf,cmat,x,y,nx,ny,
     .      xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
     .      xfmt,yfmt,dfmt,ncasc)
      endif
c
      return
      end
ccc
      subroutine asc_writex(outroot,lout,ascsuf,lsuf,cmat,x,y,nx,ny,
     .   xlab,ylab,xmin,xmax,ymin,ymax,nf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt,ncall)
c
c: This subroutine outputs the array cmat and the axes x and y
c: to an ASCII file named outroot(1:lout)//ascsuf(1:lsuf)//.asc
c: If xmin not equal to xmax, x(1:nx) is filled from xmin to xmax.
c: Likewise with ymin,ymax and y(1:ny).  If nf is not equal to 
c: zero, a message is written to that file number.
c
c: This subroutine writes out an ascii file of the data in cmat.
c
      implicit none
      include 'Parms_com'
      integer*4 ncall,lout,lsuf,nx,ny,nf,j,jx,ldat
      character*64 outroot,ascsuf,dlab,xunit,yunit,dunit,
     .   xfmt,yfmt,dfmt,dataname
      character*30 xlab,ylab
      real*4 cmat(ny,nx),x(nx),y(ny),xmin,xmax,ymin,ymax,delx,dely
      real*8 rc_code
c
      dataname=outroot(1:lout)//'_asc'//ascsuf(1:lsuf)
      ldat=lout + 4 + lsuf
c
      if(nx*ny .le. 0) then
         write(nf,110) dataname(1:ldat),nx,ny
110      format('OUTPUT ASCII FILE = ',a,'; # ROWS =',i6,'; # COLS =',
     .      i6/'  NOT WRITTEN OUT SINCE ZERO SIZE')
         return
      endif
c
      if(xmin .ne. xmax .and. xmin .ne. -999.) then
         delx=(xmax-xmin)/float(max0(1,nx-1))
         do j=1,nx
            x(j)=xmin + (j-1)*delx
         enddo
      endif
      if(ymin .ne. ymax .and. ymin .ne. -999.) then
         dely=(ymax-ymin)/float(max0(1,ny-1))
         do j=1,ny
            y(j)=ymin + (j-1)*dely
         enddo
      endif
c
      open(83,file=dataname(1:ldat),status='unknown',form='formatted')
c: Include number of rows (ny) and number of columns (nx) in dummy word in
c: corner.  Decode as ny=int(rc_code); nx=10000*(rc_code-ny).
      rc_code=ny + nx/10000.d0
      write(83,200) rc_code,(y(j),j=1,ny)
200   format(2000(e14.8,1x))
      do jx=1,nx
         write(83,200) x(jx),(cmat(j,jx),j=1,ny)
      enddo
      close(83)
c
      if(nf .ne. 0) then
         write(nf,100) dataname(1:ldat),nx,ny
100      format('ASCII FILE = ',a,'; # ROWS =',i6,'; # COLS =',i6)
      endif
      ncall=ncall + 1
c
      return
      end
ccc
