      integer mr,mz,mp
      parameter (mr=100,mz=20000,mp=10)
        real freq,zs,zr,
     >    rmax,dr,dr0,
     >    zmax,dz,zmplt,
     >    c0,rs0,
     >    rb(mr),zb(mr),omega
      integer ndr,ndz,np,ns,nz,nzplt
      common /ramsaga/
     >    rmax,dr,dr0,ndr,
     >    dz,ndz,nz,nzplt,
     >    c0,np,ns,rs0,
     >    rb,zb,omega
      common /tosetup/freq,zs,zr, zmax,zmplt

      real rdmin,rdmax,rarray
      integer irdmin,nbathy
      common /ramrange/rdmin,rdmax,irdmin,nbathy,rarray
 
      real pi
      common /div/pi

      integer flagpu
      common /flagpu/flagpu

      integer  ichangedepth
      common /depth/ichangedepth

