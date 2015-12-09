      subroutine getdat2tilt
c
c: this subroutine gets the second batch of input data.
c
      implicit integer*4(i-n)
      include 'common/options'
      include 'common/svp'
      include 'common/caustix'
      include 'common/pii'
      include 'common/freqcom'
      include 'common/srlox'
      include 'common/depth'
      include 'common/bdcom'
      include 'common/traxcom'
      include 'common/tilt'
      integer*4 jzr1
c: nbdfmax is the max # of freq at which to calc bd eigs, fsfac is the
c: factor that controls the spacing of the frequecies:
      data nbdfmax/19/,fsfac/2.5/
c
cpln      read(1,*,end=500,err=500) nrangx
      if(nrangx .gt. 0) then
         iirr=1
         if(nrangx .gt. 101) then
            print *,'# ranges limited to 101 when you list them.'
            print *,'Check line number ',nline,' in .opt file...'
            stop
         endif
cpln         read(1,*,end=500,err=500) nrangx,(rtmp(j),j=1,nrangx),
cpln     .        dtiltv,dtilth
      else
         iirr=-1
cpln         read(1,*,end=500,err=500) nrangx,rtmp(1),drtmp,dtiltv,dtilth
         nrangx=iabs(nrangx)
      endif
c
c: read in uniform array locations: 
cpln         read(1,*,end=500,err=500) zr1,nx,ny,nzr,dx,dy,dz(1)
c
      if(zr1 .lt. 0.) then
         print *,'negative rec depth ',zr1,' measured from ocean ',
     .        'bottom: ',zsvp(nsvp) + zr1
         zr1=zsvp(nsvp) + zr1
      endif
      dx=0
      dy=0
      nx=1
      ny=1
      nxny=nx*ny 
      if((nx .lt. 1) .or. (ny .lt. 1) .or. (nzr .lt. 1)
     .     .or. (nzr .gt. 128) .or. (nxny .gt. 32)) then
         print *,'illegal nx,ny,nz for iiarr=1 in opt file: ',
     .        nx,ny,nzr
         print *,'Check line number ',nline,' in .opt file...'
         stop
      endif
      dxtot=float(nx-1)*dx 
      dytot=float(ny-1)*dy 
      nrec=nxny*nzr
      jxy=0
      do 20 jzr1=1,nzr
c     : zr1 is depth of top of array; recs are numbered beginning from top:
         zr0=zr1 + float(jzr1-1)*dz(1)
         nxy(jzr1)=nxny
         kxy(jzr1)=jxy
c     : fill xyz(jxy,i) with (x,y) coordinates from center of array: 
         do 22 jy=1,ny
            do 24 jx=1,nx
               jxy=jxy + 1 
               xyz(jxy,1)=float(jx-1)*dx - dxtot/2.
               xyz(jxy,2)=float(jy-1)*dy - dytot/2.
               xyz(jxy,3)=zr0
c     pln                  write(6,*)'zr0,zr1,xyz: ',zr0,zr1,xyz(jxy,3)
 24         continue
 22      continue
 20   continue
c: read in nonuniform array locations:  
c
c: nrtot is the total number of ranges in range() for each zs.
      nrtot=nrec*nrangx
c     nrnf=max0(nrtot,nfr2)
c
      return
      end 
