      program ramgeo
c
      open(unit=1,status='old',file='ramgeo.in')
      open(unit=32,status='unknown',file='tlgeo.line')
      open(unit=33,status='unknown',file='tlgeo.grid',form='unformatted')
c
      call input

      call forw
      write(2,*)'second call'
      write(*,*)'second call'
      call forw

      close(1)
      close(2)
      close(3)
      end
     
      subroutine forw
c
c     This version of ram handles multiple sediment layers that 
c     parallel the bathymetry.
c
c     ******************************************************************
c     ***** Range-dependent Acoustic Model, Version 1.4g, 05-Jul-00 ****
c     ******************************************************************
c     
c     This code was developed by Michael D. Collins at the Naval
c     Research Laboratory in Washington, DC. It solves range-dependent 
c     ocean acoustics problems with the split-step Pade algorithm 
c     [M. D. Collins, J. Acoust. Soc. Am. 93, 1736-1742 (1993)]. A 
c     user's guide and updates of the code are available via anonymous 
c     ftp from ram.nrl.navy.mil. 
c
c     Version 1.5g contains a correction to a bug in the dimension of
c     quantities passed to subroutines fndrt and guerre that Laurie
c     Fialkowski noticed. 
c
c     Version 1.4g contains a correction to a minor bug in subroutine
c     guerre that Dave King noticed (amp1 and amp2 were declared
c     twice) and a few other minor improvements. 
c
c     Version 1.3g contains a new root-finding subroutine.
c
c     Version 1.2g contains a minor modification. The output to tl.grid
c     is no longer zeroed out along the ocean bottom. This was done in
c     previous versions so that the ocean bottom would be highlighted
c     in graphical displays. The graphics codes ramclr, ramctr, and
c     ramcc read in the bathymetry from ram.in and plot the ocean
c     bottom directly. 
c
c     Version 1.1g contains two improvements:
c
c     (1) An improved self starter. Stability is improved by using the 
c     factor (1-X)**2 instead of (1+X)**2 to smooth the delta function. 
c     The factor (1+X)**2 is nearly singular for some problems involving
c     deep water and/or weak attenuation. Numerical problems associated 
c     with this singularity were detected by Eddie Scheer of Woods Hole 
c     Oceanographic Institute. 
c
c     (2) Elimination of underflow problems. A very small number is 
c     added to the solution in subroutine solve to prevent underflow,
c     which can adversely affect run time on some computers. This
c     improvement was suggested by Ed McDonald of the SACLANT Undersea
c     Research Centre. 
c
c
c     mr=bathymetry points, mz=depth grid, mp=pade terms.
c
c      parameter (mr=100,mz=20000,mp=10)
      INCLUDE 'comopt.h'
      INCLUDE 'comforw.h'
      include '../comramgeo.h'
      real cw(mz),cb(mz),rhob(mz),attn(mz),alpw(mz),
     >   alpb(mz),f1(mz),f2(mz),f3(mz),ksqw(mz),tlg(mz)
      complex ci,ksq(mz),ksqb(mz),u(mz),
     >   v(mz),r1(mz,mp),r2(mz,mp),r3(mz,mp),s1(mz,mp),
     >   s2(mz,mp),s3(mz,mp),pd1(mp),pd2(mp)
      real k0,r,eps,dir,eta,rp,rs
      integer ir,lz,mdr,ib
      integer iz


      call setup(mr,mz,nz,mp,np,ns,mdr,ndr,ndz,iz,nzplt,lz,ib,ir,dir,dr,
     >   dz,pi,eta,eps,omega,rmax,c0,k0,ci,r,rp,rs,rb,zb,cw,cb,rhob,
     >   attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,u,v,r1,r2,r3,s1,s2,s3,
     >   pd1,pd2,tlg)

c
c     March the acoustic field out in range.
c
    1 call updat(mr,mz,nz,mp,np,iz,ib,dr,dz,eta,omega,rmax,c0,k0,ci,r,
     >   rp,rs,rb,zb,cw,cb,rhob,attn,alpw,alpb,ksq,ksqw,ksqb,f1,f2,f3,
     >   r1,r2,r3,s1,s2,s3,pd1,pd2)
c
c-----
c
c      write(96,*)' range',r
c     do i=1,nz+2 
c         write(96,*)v(i)
c      enddo
      call solvegeo(mz,nz,mp,np,iz,u,v,r1,r2,r3,s1,s2,s3)
      r=r+dr
c      write(91,*)' range',r
c      do i=1,nz+2 
c         write(91,*)f1(i),f2(i),f3(i)
c      enddo
c      write(93,*)' range',r
c      do i=1,nz+2 
c         write(93,*)r1(i,1),r2(i,1),r3(i,1),r1(i,2),r2(i,2),r3(i,2)
c      enddo
c      write(94,*)' range',r
c      do i=1,nz+2 
c         write(94,*)s1(i,1),s2(i,1),s3(i,1),s1(i,2),s2(i,2),s3(i,2)
c      enddo
      call outpt(mz,mdr,ndr,ndz,iz,nzplt,lz,ir,dir,eps,r,f3,u,tlg)
c

      if(r.lt.rmax)go to 1
c      write(95,*)' range',r
c      do i=1,nz+2 
c         write(95,*)u(i)
c      enddo
c
c
      return
      end
c
c     Initialize the parameters, acoustic field, and matrices.
c
      subroutine input
      include  '../comramgeo.h'
      include 'ram.h'
      integer i,id,ic,irange,itype 
      real   zdep,zval,zvalold
       read(1,*)
      read(1,*)freq,zs,zr
      pi=4.0*atan(1.0)
      omega=2.0*pi*freq
      read(1,*)rmax,dr,ndr
      read(1,*)zmax,dz,ndz,zmplt
      read(1,*)c0,np,ns,rs0
      nz=zmax/dz-0.5
      nzplt=zmplt/dz-0.5

c---- read the bathymetry
      i=1
    1 read(1,*)rb(i),zb(i)
      if(rb(i).lt.0.0)go to 2
      i=i+1
      go to 1
    2 rb(i)=2.0*rmax
      zb(i)=zb(i-1)
c---- check dimensions
      if(nz+2.gt.mz)then
         write(*,*)'   Need to increase parameter mz to ',nz+2
         stop
      end if
      if(np.gt.mp)then
         write(*,*)'   Need to increase parameter mp to ',np
         stop
      end if
      if(i.gt.mr)then
         write(*,*)'   Need to increase parameter mr to ',i
         stop
      end if
      rangeprof(1)=0
      do irange=1,mrange
         do itype=1,4
            ic=1  
            zval=0
            do id=1,mpardep
             zvalold=zval
             read(1,*,err=20) zdep,zval
             geodep(ic,itype,irange)=zdep
             geopar(ic,itype,irange)=zval
             igeodep(ic,itype,irange)= zdep/dz+1.5

             if ((id.eq.1).and. (igeodep(id,itype,irange).ne.1)) then
                igeodep(ic+1,itype,irange)=igeodep(ic,itype,irange)
                geopar(ic+1,itype,irange)=geopar(ic,itype,irange)
                igeodep(ic,itype,irange)=1
                ic=ic+1
             endif

             if (geodep(ic,itype,irange).lt.0) then
                geodep(ic,itype,irange)=(nz+1)*dz
                geopar(ic,itype,irange)=zvalold
                igeodep(ic,itype,irange)=nz+2
                goto 22
             elseif (igeodep(ic,itype,irange).lt.(nz+2)) then
               if (ic.ge.2) then
               if (igeodep(ic,itype,irange).eq.
     #                 igeodep(ic-1,itype,irange)) then
                  igeodep(ic,itype,irange)=igeodep(ic-1,itype,irange)+1
               endif
               endif
               ic=ic+1
             endif
           enddo
 22        ndepprof(irange,itype)=ic
           if (ic.gt.Mpardep ) then
              write(*,*) 'Need to increase parameter',
     #                   ' Mpardep to ',(Nprof+1)
         stop
      endif
c234567890123456789012345678901234567890123456789012345678901234567890
        enddo
        read(1,*,err=20)rangeprof(irange+1)
      enddo
 20   Nprof=irange
      rangeprof(Nprof+1)=2*rmax
      if ((Nprof+1).gt.Mrange) then
         write(*,*) 'Need to increase parameter Mrange to ',(Nprof+1)
         stop
      endif

c
c----
c
      do irange=1,Nprof
        write(*,*)'Profile at range:',rangeprof(irange)
         do itype=1,4
           If (itype.eq.1) then
              write(*,*)' Ocean Sound Speed'
           elseIf (itype.eq.2) then
              write(*,*)' Bottom Sound Speed'
           elseIf (itype.eq.3) then
              write(*,*)' Bottom Density'
           elseIf (itype.eq.4) then
              write(*,*)' Bottom Attenuation'
           endif
           write(*,*)'Number of points:',ndepprof(irange,itype)

            do id=1,ndepprof(irange,itype)
c             write(*,*)'id,itype,irange',id,itype,irange
             write(*,*) 
     #              geodep(id,itype,irange),geopar(id,itype,irange),
     #              igeodep(id,itype,irange)
           enddo
        enddo
      enddo


      end    
