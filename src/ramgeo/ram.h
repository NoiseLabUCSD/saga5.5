c common file for RAM
      integer  Mrange,Mpardep 
      parameter  (Mrange=10,Mpardep=200)
      real*8 geopar(Mpardep,4,Mrange)
      real     geodep(Mpardep,4,Mrange),
     >   rangeprof(Mrange)
      integer   igeodep(Mpardep,4,Mrange)
      integer   ndepprof(Mrange,4),iprofile,nprof

      common /geopar/ rangeprof,geodep,geopar,igeodep, 
     >                ndepprof, Nprof,iprofile
      real dr2,rmin
      common /rangestep/dr2,rmin







     
