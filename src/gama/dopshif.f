      subroutine dopshif(a,kr,jcaus)
c
c: this subroutine shifts the frequencies frq such that due to
c: the doppler shift the frequency at the receiver would be frx.
c
      implicit integer*4(i-n)
      include 'common/freqcom'
      include 'common/srlox'
      include 'common/pii'
c
c     fac=1. - a*vs*cpsi(kr)
c     if(jcaus .gt. 0) then
c        do 10 kfr=1,nfr2
c           frq(kfr)=frx(kfr)*fac
c           wom(kfr)=twpie*frq(kfr)
c           w16(kfr)=wom(kfr)**(.16667) 
c           w23(kfr)=wom(kfr)**(.66667) 
c           frsq(kfr)=frq(kfr)**2
10       continue
c     else
c        do 20 kfr=1,nfr2
c           frq(kfr)=frx(kfr)*fac
c           wom(kfr)=twpie*frq(kfr)
c           frsq(kfr)=frq(kfr)**2
20       continue
c     endif
c     print *,'doppler: ',fac,1./a,vs,acos(cpsi(kr))/pierad 
c     print *,'dopshif frx= ',(frx(j),j=1,nfr2)
c     print *,'frq = ',(frq(j),j=1,nfr2)
c
      return
      end 
