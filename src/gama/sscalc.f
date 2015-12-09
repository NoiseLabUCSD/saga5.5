      subroutine sscalc(nb,iisl,a,aint,alim,cint,sigma,wind,nf1,nf2,
     .   freq,ss,ang,frbl,dbbl,phbl,nf,na)
c
c: this subroutine calculates the surface scattering for the ray
c: ray with snell invariant "a" at the frequencies in array freq: 
c
      implicit integer*4(i-n)
      include 'common/pii'
      real*4 ss(2,nf2),ffac,ffac1,afac,afac1,fdif,fa00,fa01,fa10,fa11,
     .   dbfit,phfit,ssmx,ssreal
      real*4 ang(0:na),frbl(0:nf),dbbl(0:nf,0:na),phbl(0:nf,0:na)
      real freq(nf2)
      data ssmx/.2818383/
c
      if(iisl .eq. 0) return
c     print *,'sscalc: iisl,nf1,nf2 = ',iisl,nf1,nf2,freq(nf1:nf2)
      if((a .gt. alim) .or. (nb .eq. 0)) then
         do 10 kfr=nf1,nf2
            ss(1,kfr)=1.
            ss(2,kfr)=0.
10       continue
      elseif(iisl .eq. 1) then
c: eckart scattering:
         sinsq=1. - (a*cint)**2
         arg=nb*(-2.*(sinsq*(twpie*sigma*aint)**2))
         do 20 kfr=nf1,nf2
            ss(1,kfr)=exp(arg*freq(kfr)**2)
            ss(2,kfr)=0.
20       continue
      elseif(iisl .eq. 2) then
c: modified eckart scattering:
         sinsq=1. - (a*cint)**2
         arg=-2.*(sinsq*(twpie*sigma*aint)**2)
         do 30 kfr=nf1,nf2
            arg2=arg*freq(kfr)**2
            ssreal=(besl0(-arg2)*exp(arg2))**nb
            ss(1,kfr)=max(ssreal,ssmx)
            ss(2,kfr)=0.
30       continue
      elseif(iisl .eq. 3) then
c: beckmann-spizzochino scattering:
         thg=acos(a*cint)
         sinth=sin(thg)
         smla=500./(3.0 + 2.6*wind)
         t=smla*thg**2/4.
         v=sinth - (sinth/thg)*exp(-t)/sqrt(smla*pie)
         sqvnb=sqrt(1. - min(0.99,max(v,sinth/2.)))**nb
         frqfac=1.16e-5*wind**2
         do 40 kfr=nf1,nf2
c: set frequency and wind speed parameters
            frqspd=frqfac*freq(kfr)
            cons=.3 + .7/(1. + .01*frqspd**2)
c: compute the scattering coefficient
            ss(1,kfr)=cons*sqvnb
            ss(2,kfr)=0.
40       continue
      elseif(iisl .eq. 4) then
c: uniform scattering:
         ek=sigma**nb
         do 50 kfr=nf1,nf2
            ss(1,kfr)=ek
            ss(2,kfr)=0.
50       continue
      elseif(iisl .eq. 5) then
c: surface/bottom loss scattering table:
         thg=piedeg*acos(a*cint)
         do 60 jabl=1,na
            if(thg .le. ang(jabl)) goto 62
60       continue
c: find interpolation factors in angle:
62       continue
         afac=(thg - ang(jabl-1))/(ang(jabl) - ang(jabl-1))
         afac1=1. - afac
         jfbl=0
         do 70 kfr=nf1,nf2
c: find interpolation factors in frequency:
            if(freq(kfr) .gt. frbl(jfbl)) then
68             jfbl=jfbl + 1
               if(freq(kfr) .gt. frbl(jfbl)) goto 68
               fdif=frbl(jfbl) - frbl(jfbl-1)
            endif
            if(freq(kfr) .lt. frbl(jfbl-1)) then
69             jfbl=jfbl-1
               if(freq(kfr) .lt. frbl(jfbl-1)) goto 69
               fdif=frbl(jfbl) - frbl(jfbl-1)
            endif
            ffac=(freq(kfr) - frbl(jfbl-1))/fdif
            ffac1=1. - ffac
c: use a four-point bivariate interpolation scheme:(A&S, p. 882)
            fa00=ffac*afac
            fa10=ffac1*afac
            fa01=ffac*afac1
            fa11=ffac1*afac1
            dbfit=fa11*dbbl(jfbl-1,jabl-1) + fa01*dbbl(jfbl,jabl-1) +
     .         fa10*dbbl(jfbl-1,jabl) + fa00*dbbl(jfbl,jabl)
            phfit=(fa11*phbl(jfbl-1,jabl-1) + fa01*phbl(jfbl,jabl-1) +
     .         fa10*phbl(jfbl-1,jabl) + fa00*phbl(jfbl,jabl))*pierad
            ss(1,kfr)=10.**(-nb*dbfit/20)
            ss(2,kfr)=nb*phfit 
c     print *,'angles: ',jabl,ang(jabl-1),thg,ang(jabl)
c     print *,'freqs: ',jfbl,frbl(jfbl-1),freq(kfr),frbl(jfbl)
c     print *,'dbfit = ',dbfit
c     print *,'db ang 1: ',dbbl(jfbl-1,jabl-1),dbbl(jfbl,jabl-1)
c     print *,'db ang 2: ',dbbl(jfbl-1,jabl),dbbl(jfbl,jabl)
c     print *,'phfit = ',phfit*piedeg,ss(2,kfr)*piedeg
c     print *,'ph ang 1: ',phbl(jfbl-1,jabl-1),phbl(jfbl,jabl-1)
c     print *,'ph ang 2: ',phbl(jfbl-1,jabl),phbl(jfbl,jabl)
70       continue
      endif
c
c     print *,'ss(1,:) = ',ss(1,nf1:min0(10,nf2))
c     print *,'ss(2,:) = ',ss(2,nf1:min0(10,nf2))
c     print *,'phbl(0:nfbl,0:nabl) = ',phbl(0:nf,0:na)
c
      return
      end
