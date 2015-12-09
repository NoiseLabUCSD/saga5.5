      subroutine traxpl(xs,ys,tyme,rangx,n)
c
c: this subroutine plots the source tracks as viewed from above.
c
      implicit integer*4(i-n)
      include 'common/srlox'
      include 'common/charcom'
      real xs(n),ys(n),tyme(n),rangx(n)
      character*32 labx,laby,labb
      data labx/'X (KM)    '/,laby/'Y (KM)    '/,labb/'          '/
c
      open(52,file=outn(1:loutn)//'_trx',status='unknown')
      xsmin=1.e13
      xsmax=-1.e13 
      ysmin=1.e13
      ysmax=-1.e13 
      do 10 kr=1,n
         xsmin=min(xsmin,xs(kr))
         ysmin=min(ysmin,ys(kr))
         xsmax=max(xsmax,xs(kr))
         ysmax=max(ysmax,ys(kr))
10    continue
      do 12 jzr=1,nzr
         do 14 jxy=1,nxy(jzr) 
            xsmin=min(xsmin,xyz(kxy(jzr)+jxy,1))
            xsmax=max(xsmax,xyz(kxy(jzr)+jxy,1))
            ysmin=min(ysmin,xyz(kxy(jzr)+jxy,2))
            ysmax=max(ysmax,xyz(kxy(jzr)+jxy,2))
14       continue
12    continue
      xsmin=xsmin - 1.
      xsmax=xsmax + 1.
      ysmin=ysmin - 1.
      ysmax=ysmax + 1.
c
      call axlabel(xsmin,xsmax,xslo,xshi,xdel)
      call axlabel(ysmin,ysmax,yslo,yshi,ydel)
c: make plot come out true scale: 
      hfac=min(10.75/(xshi-xslo),7.5/(yshi-yslo))
      vfac=hfac
      hdim=(xshi-xslo)*hfac
      vdim=(yshi-yslo)*vfac
c
cxx   call pltlfn(l"raypic")
cxx   call pltdim(11.75,9.0,-1)
cxx   call pltorg(.5,.5)
cxx   call pltaxis(0.,0.,hdim,0.,xslo,xshi,xdel,labx,-5,1)
cxx   call pltaxis(0.,0.,vdim,90.,yslo,yshi,ydel,laby,5,1)
cxx   call pltaxis(0.,vdim,hdim,0.,xslo,xshi,xdel,labb,1,0) 
cxx   call pltaxis(hdim,0.,vdim,90.,yslo,yshi,ydel,labb,-1,0)
      nl=3
      do 30 kr=1,nrangx
         xpt=(xs(kr) - xslo)*hfac
         ypt=(ys(kr) - yslo)*vfac
cxx      call plt(xpt,ypt,nl) 
cxx      call symbol(xpt,ypt,.12,3,0.,-1)
c: temporary write to file while plotting unavailable:
         write(52,200) .001*xs(kr),ctab,.001*ys(kr),ctab,
     .      tyme(kr),ctab,.001*rangx(kr)
200      format(3(f14.6,a1),f14.6)
         nl=2
30    continue
      do 40 jzr=1,nzr
         do 42 jxy=1,nxy(jzr) 
            xpt=(xyz(kxy(jzr)+jxy,1) - xslo)*hfac
            ypt=(xyz(kxy(jzr)+jxy,2) - yslo)*vfac
cxx         call symbol(xpt,ypt,.12,4,0.,-1)
c: temporary write to file while plotting unavailable:
            write(52,200) .001*xyz(kxy(jzr)+jxy,1),ctab,
     .         .001*xyz(kxy(jzr)+jxy,2),ctab,xyz(kxy(jzr)+jxy,3)
42       continue
40    continue
c
cxx   call pltline(hdim/2.,vdim+.400,-.16)
cxx   write(8,100) (srloc(j,1),j=1,min0(nzs,10))
100   format('zs = ',10(f9.4,2x))
cxx   if(nzs .gt. 10) then
cxx      call pltline(hdim/2.,vdim-.135,-.16)
cxx      write(8,100) (srloc(j,1),j=11,nzs)
cxx   endif
cxx   call pltline(hdim/2.,vdim+.135,-.16)
cxx   write(8,110) (srloc(j,2),j=1,min0(nzr,10))
110   format('zr = ',10(f9.2,4x))
cxx   if(nzr .gt. 10) then
cxx      call pltline(hdim/2.,vdim-.400,-.16)
cxx      write(8,100) (srloc(j,2),j=11,nzr)
cxx   endif
cxx   call pltend(0.)
      close(52)
c
      return
      end 
