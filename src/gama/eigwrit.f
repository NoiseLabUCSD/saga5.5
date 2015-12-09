      subroutine eigwrit(nr,jcaus,a,travtot,exptot,trm,trp, 
     .   gsl,dbattn,nps,r,dr2,delta,dellam,ssl,bsl,phr,ibd,
     .   ntlay,NLAYM)
c
c: this subroutine is called by eigbar and writes out eigenray
c: information to the eigenray list file.
c
      implicit integer*4(i-n)
      include 'common/bottom' 
      include 'common/pathchr'
      include 'common/paths'
      include 'common/pathway'
      include 'common/pii'
      include 'common/svp'
      include 'common/caustix'
      include 'common/gamaoptions'
c
      character*1 cq(0:3),cbd(0:1),cz 
      integer*4 zbuf(42)
      real*4 ssl,bsl
      integer*4 ntlay(-1:NLAYM+2,mbp)
      data cq/' ','$','&','!'/,cbd/' ','_'/
      data zbuf/42*0/
c
      if(iidiag .ne. 0) print *,'starting eigwrit: ',nr,jcaus,a,
     .   travtot,exptot,trm,trp,gsl,dbattn,nps,r,dr2,delta,dellam,
     .   ssl,bsl,phr,ibd
      cz=czprop(icpth)
      ndp=ntlay(-1,mbp)
cxx   print *,'eigwrit: mbp,ndp,ntlay = ',mbp,ndp,ntlay(0:ndp,mbp)
c
      ths=acos(min(1.,a*cs))*piedeg
      if(fdir(icpth) .eq. 'd') ths=-1.*ths
      thr=acos(min(1.,a*cr))*piedeg
      if(ldir(icpth) .eq. 'u') thr=-1.*thr
      trmdb=20.*log10(trm)
      gsldb=20.*log10(amax1(1.e-20,gsl))
      ssldb=20.*log10(ssl)
      bsldb=20.*log10(bsl)
      fdpdb=exptot*freq/(-1000.)
      iitrp=nint(trp)
      iiphr=nint(phr)
      iitrm=nint(trmdb)
      iigsl=nint(gsldb)
      iissl=nint(ssldb)
      iibsl=nint(bsldb)
      iifdp=nint(fdpdb)
      ainv=1./a
      svii=min(ainv,9999.9) 
      nlim=ndp
      do 44 j=1,ndp
         if(ainv .lt. cp1(j)) then
            nlim=j-1
            goto 45
         endif
44    continue
45    continue
      write(59,140) nr,cq(jcaus),cbd(ibd),fdir(icpth),
     .   ntop,nbotx,cz,(ntlay(j,mbp)/2,j=1,min0(nlim,5)),
     .   (zbuf(j),j=nlim+1,5),ths,thr,svii,travtot,iifdp,iitrm,iitrp,
     .   iigsl,iissl,iibsl,dbattn,iiphr,nps
140   format(i4,a1,a1,1x,a1,1x,i2,1x,i2,1x,a1,1x,i2,4(i2.0),
     .   2(1x,f6.2),1x,f6.1,1x,f9.5,1x,6(i4,1x),f7.2,1x,i4,1x,i2)
      if(nlim .gt. 5) then
         write(59,142) (ntlay(j,mbp)/2,j=6,nlim),(zbuf(j),j=nlim+1,42)
142      format(' other layers: ',i2,36(i2.0))
      endif
c
      if(jcaus .eq. 1) then
         write(59,180) r,delta,dellam,dr2
180      format(4x,'CAUSTIC R = ',f10.2,'; 40-dB SHADOW ZONE ',
     .      'EXTENT = ',f10.2,' M = ',f6.0,' LAMBDA; R" = ',e10.4)
      endif
c
      return
      end 
