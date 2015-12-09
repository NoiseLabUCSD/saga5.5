      subroutine acouspl
c
c: this subroutine plots as a function of depth the sound speed,
c: attenuation, and density for the ocean/layered-layered-bottom
c: environment the plot appears to the right of the ray picture.
c
      implicit integer*4(i-n)
      include 'common/attcon' 
      include 'common/svp'
      include 'common/pic'
      include 'common/bottom' 
      include 'common/depth'
      character*32 svlab,akplab,rholab
      data svlab/'SOUND SPEED (M/S) +'/
      data rholab/'DENSITY (G/CC) 0'/
      data akplab/'ATTN (DB/M-KHZ) X'/
c
      nptot=1000
      svmin=cmin
      svmax=max(cmax(1),cmax(2),cmax(3))
      call axlabel(svmin,svmax,svlo,svhi,svdelt)
      fach=h2dim/(svhi-svlo)
c
cxx   call pltorg(0.,hdim)
cxx   call pltaxis(0.,0.,h2dim,90.,svlo,svhi,svdelt,svlab,19,2)
cxx   call symbol((zsvp(0)-zlo)*zfac,(csvp(0)-svlo)*fach,.12,3,90.,-1)
      do 15 k=1,klev(3)
         if(kseg .eq. 3) then 
            h=zsvp(k) - zsvp(k-1)
            klay=max0(int(nptot*h/zlev(4+nlev)),3) 
            fac=h/float(klay) 
            dc2ok=(bsvp(1,k) - bsvp(2,k-1))/h
c           do 17 kk=1,klay
            do 17 kk=0,klay
               zinc=float(kk)*fac
               zinc2=h-zinc
               call cdcwine(csvp(k-1),g0(k,2),g1(k,2),g2(k,2),zinc,
     .            cw,dcw,dcw2,2)
cxx            call plt((zsvp(k-1)+zinc-zlo)*zfac,
cxx  .                  (cinc-svlo)*fach,2)
c: temporary write to file while plotting unavailable:
               write(51,100) zsvp(k-1)+zinc,ctab,cw,ctab,dcw,ctab,
     .            dcw2,ctab,dc2ok,ctab,dcw2/dc2ok
100            format(f8.2,a1,f8.2,a1,e12.6,a1,f9.7,a1,f9.7,a1,f9.3)
17          continue
         endif
c: temporary write to file while plotting unavailable:
         write(51,100) zsvp(k),ctab,csvp(k),ctab,1.
cxx      call symbol((zsvp(k)-zlo)*zfac,(csvp(k)-svlo)*fach,.12,3,
cxx  .        90.,-2)
15    continue
c
cxx   call symbol((zfsvp(0)-zlo)*zfac,(cfsvp(0)-svlo)*fach,.12,15,
cxx  .      90.,-1) 
c: plot original s.v. data points (as opposed fitted ones)
      do 19 k=1,nfsvp
cxx      call symbol((zfsvp(k)-zlo)*zfac,(cfsvp(k)-svlo)*fach,.12,15, 
cxx  .      90.,-1) 
19    continue
c
c: don't plot attenuation if there are no bottom layers.
      if(ntot .eq. 0) goto 99 
      svmin=1.e38
      svmax=0.
      rhomin=1.e38 
      rhomax=0.
      akpmin=1.e38
      akpmax=0.
      do 13 j=1,ntot
         svmin=min(svmin,cp1(j),cp2(j))
         svmax=max(svmax,cp1(j),cp2(j))
         rhomin=min(rhomin,rho1(j),rho2(j))
         rhomax=max(rhomax,rho1(j),rho2(j))
         akpmin=min(akpmin,akp1(j),akp2(j))
         akpmax=max(akpmax,akp1(j),akp2(j))
13    continue
      call axlabel(svmin,svmax,svlo,svhi,svdelt)
      call axlabel(rhomin,rhomax,rholo,rhohi,rhodelt)
      call axlabel(akpmin,akpmax,akplo,akphi,aldelt)
      fach=h2dim/(svhi-svlo)
      fach2=h2dim/(rhohi-rholo)
      fach3=h2dim/(akphi-akplo)
      vax=(zlev(4) - zlo)*zfac
cxx   call pltaxis(vdim,0.,h2dim,90.,svlo,svhi,svdelt,svlab,-19,2)
cxx   call pltaxis(vax,0.,h2dim,90.,akplo,akphi,aldelt,
cxx  .        akplab,+17,2)
cxx   call pltaxis(vax-.6,0.,h2dim,90.,rholo,rhohi,
cxx  .     rhodelt,rholab,+16,2)
c
c: plot the three parameters in the bottom layers.
      do 20 j=1,ntot
         zst=(zlev(j+3)-zlo)*zfac
         zend=(zlev(j+4)-zlo)*zfac
cxx      call symbol(zst,(cp1(j)-svlo)*fach,.12,3,90.,-1)
c: temporary write to file while plotting unavailable:
         write(51,100) zlev(j+3),ctab,cp1(j),ctab,akp1(j)
         if(kprof(j) .eq. 2) then
            klay=max0(int(nptot*z(j)/zlev(4+nlev)),3)
            fac=z(j)/float(klay)
            do 30 k=1,klay-1
               zinc=float(k)*fac
               fac2=cp1(j)*(1.+bet(j))
               cinc=sign(1.,fac2)*sqrt(fac2**2 + 2.*bp(j)*zinc*fac2)
     .            - bet(j)*cp1(j)
               akpval=akp1(j) + zinc*dakp(j)
cxx            call plt((zlev(j+3)+zinc-zlo)*zfac,(cinc-svlo)*fach,2) 
c: temporary write to file while plotting unavailable:
               write(51,100) zlev(j+3)+zinc,ctab,cinc,ctab,akpval
30          continue
         endif
c: temporary write to file while plotting unavailable:
         write(51,100) zlev(j+4),ctab,cp2(j),ctab,akp2(j)
cxx      call symbol(zend,(cp2(j)-svlo)*fach,.12,3,90.,-2)
cxx      call symbol(zst,(akp1(j)-akplo)*fach3,.12,4,90.,-1)
cxx      call plt(zend,(akp2(j)-akplo)*fach3,5)
cxx      call symbol(zend,(akp2(j)-akplo)*fach3,.12,4,90.,-1)
cxx      call symbol(zst,(rho1(j)-rholo)*fach2,.12,15,90.,-1)
cxx      call plt(zend,(rho2(j)-rholo)*fach2,9)
cxx      call symbol(zend,(rho2(j)-rholo)*fach2,.12,15,90.,-1)
20    continue
99    continue
cxx   call pltend(0.0)
c
      return
      end 
