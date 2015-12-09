      subroutine hdf_write_gen(irank,outroot,lout,yax,ny,xax,nx,
     .   zax,nz,aax,na,bax,nb,data,jy1,jy2,jx1,jx2,jz1,jz2,
     .   ja1,ja2,jb1,jb2,ylab,ly,xlab,lx,zlab,lz,alab,la,
     .   blab,lb,dlab,ld,jset,jhdf,nf)
c
c: Writes out multi-dimensional HDF dataset of rank irank to file
c: outroot(1:lout)//'.hdf'.  Up to 5 dimension scales can be given:
c: yax(1:ny),xax(1:nx),zax(1:nz),aax(1:na),bax(1:nb).  The entire
c: output data set is of dimension (ny,nx,nz,na,nb).  The data
c: being written out on this call is contained in data and is written
c: to dimension subsets (jy1:jy2,jx1:jx2,jz1:jz2,ja1:ja2,jb1:jb2).
c: Therefore, the input data array is of dimension:
c: data(jy2-jy1+1,jx2-jx1+1,jz2-jz1+1,ja2-ja1+1,jb2-jb1+1).
c: Dimension labels are given in ylab(1:ly), etc., and the data label
c: is given in dlab(1:ld).  The dataset number is jset.  Different
c: datasets can be written to the same file as long as the given jset
c: is different each time.  jhdf is the number of the HDF file being
c: written to.  For nf>0, a summary of the dataset are
c: written to file nf.
c: 
      implicit none
      integer*4 irank,lout,ny,nx,nz,na,nb,jy1,jy2,jx1,jx2,jz1,
     .   jz2,ja1,ja2,jb1,jb2,ly,lx,lz,la,lb,ld,
     .   jset,jhdf,nf,ival,lhdf,jdim,j,iibad
      real*4 yax(ny),xax(nx),zax(nz),aax(na),bax(nb),
     .   data(jy2-jy1+1,jx2-jx1+1,jz2-jz1+1,ja2-ja1+1,jb2-jb1+1)
      character*64 outroot,ylab,xlab,zlab,dlab,alab,blab,hdf_file
      character*2 jslab
      real*4, dimension(:,:,:,:,:), allocatable :: data2
      integer start(5),edges(5),stride(5),dim_id
cc    integer sfstart,sfcreate,sfsdtstr,sfdimid,sfsdmstr,sfsdmname,
cc   .   sfsdscale,sfwdata,sfendacc,sfend,sfsdmvc,sfselect
c     ** hdf attributes
      integer dfacc_create,dfnt_float32,dfnt_int32,dfacc_write,
     .   sd_unlimited
      data start,edges,stride/15*0/
      data dfacc_create/4/
      data dfnt_float32/5/
      data dfacc_write/2/
      data sd_unlimited/0/
      data dfnt_int32/24/
      include 'sd_com'
cc      data iifirst/12*1/
c
cc    print *,'starting hdf: ',jset,jhdf
cc    print *,'dims,iifirst = ',ny,nx,nz,na,nb,iifirst(jhdf)
cc    print *,jy1,jy2,jx1,jx2,jz1,jz2,ja1,ja2,jb1,jb2
cc    call flush(6)
c
	print *, 'Warning: HDF not supported (dto) '
	return

      if(jset .gt. 50) then
         print *,'Number of datasets currently limited to 50.'
         print *,'Change values in sd_com,hdf_write_gen.'
         return
      endif
      if(jhdf .gt. 12) then
         print *,'Number of HDF files currently limited to 12.'
         print *,'Change values in sd_com,hdf_write_gen.'
         return
      endif
c

      print*,'*** irank ', irank

      write(jslab,'(i2)') jset
      ival=0
      lhdf=lout + 4
      hdf_file(1:lhdf)=outroot(1:lout)//'.hdf'
      if(iifirst(jhdf) .eq. 1) then
         sd_id(jhdf)=sfstart(hdf_file(1:lhdf),dfacc_create)
cc       print *,'ststart: ',jhdf,sd_id(jhdf),hdf_file(1:lhdf)
         do j=1,50
            set_id(j,jhdf)=0
         enddo
         nset_open(jhdf)=0
      endif
      if(set_id(jset,jhdf) .eq. 0) then
         rank(jset,jhdf)=irank
         idims(1,jset,jhdf)=ny
         idims(2,jset,jhdf)=nx
         idims(3,jset,jhdf)=nz
         idims(4,jset,jhdf)=na
         idims(5,jset,jhdf)=nb
         set_id(jset,jhdf)=sfcreate(sd_id(jhdf),hdf_file(1:lhdf),
     .      dfnt_float32,irank,idims(1,jset,jhdf))  
         nset_open(jhdf)=nset_open(jhdf)+1
         jset_open(nset_open(jhdf),jhdf)=jset
         ival=ival+sfsdtstr(set_id(jset,jhdf),dlab(1:ld),' ','f8.3',
     .      'cartesian')
         dim_id=sfdimid(set_id(jset,jhdf),0)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
         ival=ival+sfsdmname(dim_id,ylab(1:ly)//jslab)
         ival=ival+sfsdmstr(dim_id,ylab(1:ly),' ','f9.4')
         ival=ival+sfsdscale(dim_id,ny,dfnt_float32,yax)
         dim_id=sfdimid(set_id(jset,jhdf),1)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
         ival=ival+sfsdmname(dim_id,xlab(1:lx)//jslab)
         ival=ival+sfsdmstr(dim_id,xlab(1:lx),' ','f9.4')
         ival=ival+sfsdscale(dim_id,nx,dfnt_float32,xax)
         if(irank .ge. 3) then
            dim_id=sfdimid(set_id(jset,jhdf),2)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
            ival=ival+sfsdmname(dim_id,zlab(1:lz)//jslab)
            ival=ival+sfsdmstr(dim_id,zlab(1:lz),' ','f9.4')
            ival=ival+sfsdscale(dim_id,nz,dfnt_float32,zax)
         endif
         if(irank .ge. 4) then
            dim_id=sfdimid(set_id(jset,jhdf),3)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
            ival=ival+sfsdmname(dim_id,alab(1:la)//jslab)
            ival=ival+sfsdmstr(dim_id,alab(1:la),' ','f9.4')
            ival=ival+sfsdscale(dim_id,na,dfnt_float32,aax)
         endif
         if(irank .ge. 5) then
            dim_id=sfdimid(set_id(jset,jhdf),4)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
            ival=ival+sfsdmname(dim_id,blab(1:lb)//jslab)
            ival=ival+sfsdmstr(dim_id,blab(1:lb),' ','f9.4')
            ival=ival+sfsdscale(dim_id,nb,dfnt_float32,bax)
         endif
         if(nf .ne. 0 .and. iifirst(jhdf) .eq. 1) then
            write(nf,99) hdf_file(1:lhdf)
99          format('HDF FILE = ',a)
         endif
         if(nf .ne. 0) then
            write(nf,100) dlab(1:ld),ny,nx,nz,na,nb
100         format('   DATASET = ',a22,'; DIMS = ',5(i6))
         endif
      else
         iibad=0
         if(rank(jset,jhdf) .ne. irank) then
            print *,'WARNING: rank on subsequent call not same: ',
     .         irank,rank(jset,jhdf)
            iibad=1
         endif
         if(jy2 .gt. idims(1,jset,jhdf)) then
      print *,'WARNING: jy2>ny (dataset cannot grow) ',jy2,ny,jset
      print *,idims(1:rank(jset,jhdf),jset,jhdf)
            iibad=1
         endif
         if(jx2 .gt. idims(2,jset,jhdf)) then
      print *,'WARNING: jx2>nx (dataset cannot grow) ',jx2,nx,jset
      print *,idims(1:rank(jset,jhdf),jset,jhdf)
            iibad=1
         endif
         if(irank .ge. 3 .and. jz2 .gt. idims(3,jset,jhdf)) then
      print *,'WARNING: jz2>nz (dataset cannot grow) ',jz2,nz,jset
      print *,idims(1:rank(jset,jhdf),jset,jhdf)
            iibad=1
         endif
         if(irank .ge. 4 .and. ja2 .gt. idims(4,jset,jhdf)) then
      print *,'WARNING: ja2>na (dataset cannot grow) ',ja2,na,jset
      print *,idims(1:rank(jset,jhdf),jset,jhdf)
            iibad=1
         endif
         if(irank .ge. 5 .and. jb2 .gt. idims(5,jset,jhdf)) then
      print *,'WARNING: jb2>nb (dataset cannot grow) ',jb2,nb,jset
      print *,idims(1:rank(jset,jhdf),jset,jhdf)
            iibad=1
         endif
         if(iibad .eq. 1) return
c: Select existing dataset:
cc       print *,'ival sel: ',ival
cc       ival=ival + sfselect(sd_id,set_id(jset))
cc       print *,'ival se2: ',ival
         dim_id=sfdimid(set_id(jset,jhdf),0)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
         ival=ival+sfsdscale(dim_id,ny,dfnt_float32,yax)
         dim_id=sfdimid(set_id(jset,jhdf),1)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
         ival=ival+sfsdscale(dim_id,nx,dfnt_float32,xax)
         if(irank .ge. 3) then
            dim_id=sfdimid(set_id(jset,jhdf),2)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
            ival=ival+sfsdscale(dim_id,nz,dfnt_float32,zax)
         endif
         if(irank .ge. 4) then
            dim_id=sfdimid(set_id(jset,jhdf),3)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
            ival=ival+sfsdscale(dim_id,na,dfnt_float32,aax)
         endif
         if(irank .ge. 5) then
            dim_id=sfdimid(set_id(jset,jhdf),4)
         
ctms     ** set backward compatible dims
         ival=sfsdmvc(dim_id,1)
         
            ival=ival+sfsdscale(dim_id,nb,dfnt_float32,bax)
         endif     
      endif
c
      start(1)=jy1-1
      stride(1)=1
      edges(1)=jy2-jy1+1
      start(2)=jx1-1
      stride(2)=1
      edges(2)=jx2-jx1+1
      start(3)=jz1-1
      stride(3)=1
      edges(3)=jz2-jz1+1
      start(4)=ja1-1
      stride(4)=1
      edges(4)=ja2-ja1+1
      start(5)=jb1-1
      stride(5)=1
      edges(5)=jb2-jb1+1
c
      ival=ival+sfwdata(set_id(jset,jhdf),start,stride,edges,data)
cc    ival=ival+sfendacc(set_id(jset,jhdf),jhdf)
      if(ival .ne. 0) then
         print *,'ival final bad in hdf_write_gen = ',ival,jset,irank
cc    else
cc       print *,'ival final OK  in hdf_write_gen = ',ival,jset,irank
      endif
c: Only do this once at end of program:
cc    ival=ival+sfend(sd_id)
c
      iifirst(jhdf)=0
c
      return
      end
