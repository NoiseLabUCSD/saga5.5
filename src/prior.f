c     Peter Gerstoft 02
      SUBROUTINE readprior
      INCLUDE 'comopt.h'
      integer NMAX
      PARAMETER (NMAX=500)
      real xpri(nmax,mpar),ppri(nmax,mpar)
      integer npri(mpar)
      common /pri/xpri,ppri,npri
      integer eflag,ipar,ipri
      character*80 dummy

      eflag=0
      read(31,*)dummy
      do ipar=1,nparm
         read(31,*) npri(ipar)
         write(*,*) 'Prior for param',ipar
         do ipri=1,npri(ipar)
            read(31,*) xpri(ipri,ipar),ppri(ipri,ipar)
            write(*,*) xpri(ipri,ipar),ppri(ipri,ipar)
         enddo
c     check parameters
         if (xpri(1,ipar).gt.fmin(ipar)) then
            write(*,*)' ipar = ',ipar
            write(*,*)' fmin = ',fmin(ipar)
            write(*,*)' xpri = ',xpri(1,ipar)
            write(*,*)' Lower value of xpri greater than fmin'
            eflag=1
         endif      
         if (xpri(npri(ipar),ipar).lt.fmax(ipar)) then
            write(*,*)' ipar = ',ipar
            write(*,*)' fmin = ',fmax(ipar)
            write(*,*)' xpri = ',xpri(npri(ipar),ipar)
            write(*,*)' Upper value of xpri less than fmax'
            eflag=1
         endif
      enddo
      if (eflag==1) stop 'error reading prior file'
      end
      
      real function prior(fp)
      INCLUDE 'comopt.h'
      integer NMAX
      PARAMETER (NMAX=500)
      real xpri(nmax,mpar),ppri(nmax,mpar),pfac
      integer npri(mpar)
      common /pri/xpri,ppri,npri
      real fp(mpar)
      integer eflag,ipar,ipri,ip1,ip2
c     write(*,*)' entering prior'
      prior=1
      do ipar=1,nparm
         do ipri=2,npri(ipar)
            if (xpri(ipri,ipar)>fp(ipar)) goto 10
         enddo
 10      continue  
         ip1=ipri-1
         ip2=ipri
         pfac=ppri(ip1,ipar)+
     1        (ppri(ip2,ipar)-ppri(ip1,ipar))/
     2        (xpri(ip2,ipar)-xpri(ip1,ipar))*(fp(ipar)-xpri(ip1,ipar))
         prior=prior*pfac
         write(90,*)ipar,fp(ipar),pfac,prior,ip1,ip2
         write(90,*)xpri(ip1,ipar),xpri(ip2,ipar)
         write(90,*)ppri(ip1,ipar),ppri(ip2,ipar)
      enddo
      end
