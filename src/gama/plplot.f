      subroutine plplot(rax,tax,n,plco,plin,iipl,labx,ix,
     .   labx2,fr,labc,fr2,labc2,cpl)
c
c: this subroutine plots prop loss versus range.
c
      implicit integer*4(i-n)
      include 'common/timing' 
      include 'common/srlox'
      include 'common/depth'
      include 'common/pii'
      include 'common/charcom'
      real rax(n),tax(n),plin(n),amp(1000,2),ph(1000)
      real*4 cpl,c60
      complex plco(n)
      integer*4 alab(2),iline(2)
      character*10 labx,labx2,labc,labc2
cxx   data alab/"prop loss ","(db)"/,nlab/"          "/
      data iline/2,4/,c60/60./
c
      hdim=10.75
      vdim=7.5
      hmin=1.e38
      hmax=-1.e38
      do 5 k=1,n
         hmin=min(hmin,rax(k))
         hmax=max(hmax,rax(k))
5     continue
      call axlabel(hmin,hmax,rlo,rhi,rdel)
      hfac=hdim/(rhi - rlo)
      acmin=1.e38
      acmax=-1.e38 
      aimin=1.e38
      aimax=-1.e38 
      do 10 k=1,n
         if(plco(k) .ne. cmplx(0.,0.)) then
            amp(k,1)=20.*log10(abs(plco(k)))
            ph(k)=atan2(imag(plco(k)),real(plco(k)))*piedeg
         else
            amp(k,1)=-999.99
            ph(k)=0.
         endif
         acmin=min(acmin,amp(k,1))
         acmax=max(acmax,amp(k,1))
         if(plin(k) .ne. 0.) then
            amp(k,2)=10.*log10(plin(k))
         else
            amp(k,2)=-999.99
         endif
         aimin=min(aimin,amp(k,2))
         aimax=max(aimax,amp(k,2))
10    continue
      if(iipl .eq. 4) goto 88 
cxx   call pltlfn(l"raypic")
cxx   call pltdim(11.75,9.0,-1)
cxx   call pltorg(.5,.5)
c
cxx   call pltaxis(0.,0.,hdim,0.,rlo,rhi,rdel,labx,-1*ix,1) 
cxx   call pltaxis(0.,vdim,hdim,0.,rlo,rhi,rdel,nlab,1,0)
      if(iipl .eq. 3) then
         axmin=min(acmin,aimin)
         axmax=max(acmax,aimax)
         call axlabel(axmin,axmax,axlo,axhi,axdel)
         m1=1
         m2=2
      elseif(iipl .eq. 1) then
         m1=1
         m2=1
         call axlabel(acmin,acmax,axlo,axhi,axdel)
      elseif(iipl .eq. 2) then
         m1=2
         m2=2
         call axlabel(aimin,aimax,axlo,axhi,axdel)
      endif
      vfac=vdim/(axhi - axlo) 
cxx   call pltaxis(0.,0.,vdim,90.,axlo,axhi,axdel,alab,14,1)
cxx   call pltaxis(hdim,0.,vdim,90.,axlo,axhi,axdel,nlab,-1,0)
cxx   do 90 m=m1,m2 
cxx      call plt((rax(1)-rlo)*hfac,(amp(1,m)-axlo)*vfac,3) 
cxx      do 50 k=2,n
cxx         px=(rax(k) - rlo)*hfac
cxx         py=(amp(k,m) - axlo)*vfac
cxx         call plt(px,py,iline(m))
50       continue
90    continue
cxx   call pltline(hdim/2.,vdim+.135,-.16)
cxx   write(8,200) zs,zr,jxy,labc,fr,labc2,fr2
cxx200   format('zs =',f9.2,'; zr =',f9.2,'; jxy=',i2,'; ',
cxx  .   a10,' =',f9.2,'; ',a10,' =',f9.2)
cxx   call pltline(hdim/2.,vdim+.40,-.16)
cxx   write(8,202) ktitle(1:lktit),kdat(2),kdat(1),kdat(3),
cxx  .   ktim(1:3),int(cpl+.5),
cxx  .   int(rterp+.5)
cxx202   format(a,' DATE = ',2(i2,'/'),i4,' TIME = ',i2.2,
cxx  .   2(':',i2.2),' CPSEC = ',i5,' RTERP = ',i4)
cxx   call pltend(0.0)
88    continue
c
      write(12,99) labc,fr,labc2,fr2,int(cptim/c60),int(mod(cptim,c60))
99    format(/,a10,' =',f9.2,'; ',a10,' =',f9.2,'; CP MIN = ',i4,':',
     .   i2.2)
      write(12,100) n,zs,zr,jxy
100   format('N =',i5,'; ZS =',f9.2,' M; ZR =',f9.2,' M; JXY =',i2)
      write(12,104) labx,ctab,labx2,ctab,'COH PL',ctab,'COH PHASE',
     .   'INCOH PL'
104   format(a,a1,a,a1,a,a1,a,a1,a)
      write(12,110) '$$$ '
110   format(a4)
      do 120 k=1,n
         write(12,112) .001*rax(k),ctab,tax(k),ctab,amp(k,1),ctab,ph(k),
     .      ctab,amp(k,2) 
112      format(f12.3,a1,f12.3,a1,f8.3,a1,f7.2,a1,f8.3)
120   continue
cxx   write(12,110) 'EOD'
c
      return
      end 
