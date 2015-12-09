      subroutine scfft1d(n,x,iido)
c
      implicit none
c: This subroutine is a replacement for the Alliant forward FFT
c: of a real time series.
c      implicit integer*4 (i-n)
      integer*4 n,iido,j
      real*4 x(n+2)
c
      if(iido .eq. 0) return
      call realfft(x,n,1)
c: Get real(Nyquist f) from second position and zero out the 
c: imaginary parts of zero f and Nyquist f:
      x(n+1)=x(2)
      x(2)=0.
      x(n+2)=0.
c: Time reverse to get to convention for Numerical Recipes routines:
      do j=2,n,2
         x(j)=-x(j)
      enddo
c 
      return
      end
ccc
      subroutine csfft1d(n,x,iido)
c
c: This subroutine is a replacement for the Alliant inverse FFT
c: of a Hermitian spectrum.  
      implicit none
      integer*4 n,iido,j
      real*4 x(n+2)
c
      if(iido .eq. 0) return
c: Time reverse to get to convention for Numerical Recipes routines:
      do j=2,n,2
         x(j)=-x(j)
      enddo
c: Place real(Nyquist f) in imag(zero frequency component):
      x(2)=x(n+1)
      call realfft(x,n,-1)
c
      return
      end
ccc
      subroutine cfft1d(n,x,iido)
c
c: This subroutine is a replacement for the Alliant inverse FFT
c: of a general complex spectrum.
      implicit none
      integer*4 n,iido,j
      real*4 fac
      complex*8 x(n)
c
      if(iido .eq. 0) return
c: Time reverse to get to convention for Numerical Recipes routines:
      do j=1,n
         x(j)=conjg(x(j))
      enddo
c: Call with +1 instead of -1 since Alliant FFT has +i in inverse FFT:
      call four1(x,n,-1)
c: Normalize inverse FFT since four1 does not:
      fac=1./float(n)
      do 10 j=1,n
         x(j)=fac*x(j)
10    continue
c 
      return
      end
ccc
      subroutine realfft(data,n,isign)
c
c: isign=1: Calculates the in-place FFT of a real time series data(1:n).
c: The real-valued first and last components of 
c: the complex transform are returned as data(1) and data(2), resp.
c: n MUST be a power of 2.  isign=-1: Calculates the inverse
c: transform of a complex data array if it is the transform of real 
c: data (result must be multiplied by 2/n, WHICH I DO).
C: TAKEN FROM NUMERICAL RECIPES, P. 507 (NAMED REALFT).  I ALSO MULTIPLY
C: BY 2/N FOR THE INVERSE FFT
c
      implicit none
      integer*4 n,isign,n2p3,i,i1,i2,i3,i4,j
      real*4 data(n),wrs,wis,c1,c2,h1r,h1i,h2r,h2i,fac
      real*8 theta,wi,wpi,wpr,wr,wtemp
c
      theta=3.141592653589793d0/dble(n/2)
      c1=0.5
      if(isign .eq. 1) then
         c2=-0.5
c: Call with -1 since Alliant FFT has -i in forward FFT:
c??      call four1(data,n/2,-1)
         call four1(data,n/2,+1)
      else
         c2=0.5
         theta=-theta
      endif
      wpr=-2.d0*sin(0.5d0*theta)**2
      wpi=sin(theta)
      wr=1.d0 + wpr
      wi=wpi
      n2p3=n + 3
      do 10 i=2,n/4
         i1=2*i - 1
         i2=i1 + 1
         i3=n2p3 - i2
         i4=i3 + 1
         wrs=sngl(wr)
         wis=sngl(wi)
         h1r=c1*(data(i1) + data(i3))
         h1i=c1*(data(i2) - data(i4))
         h2r=-c2*(data(i2) + data(i4))
         h2i=c2*(data(i1) - data(i3))
         data(i1)=h1r + wrs*h2r - wis*h2i
         data(i2)=h1i + wrs*h2i + wis*h2r
         data(i3)=h1r - wrs*h2r + wis*h2i
         data(i4)=-h1i + wrs*h2i + wis*h2r
         wtemp=wr
         wr=wr*wpr - wi*wpi + wr
         wi=wi*wpr + wtemp*wpi + wi
10    continue
c
      if(isign .eq. 1) then
         h1r=data(1)
         data(1)=h1r + data(2)
         data(2)=h1r - data(2)
      else
         h1r=data(1)
         data(1)=c1*(h1r + data(2))
         data(2)=c1*(h1r - data(2))
c: Call with +1 since Alliant FFT has +i in inverse FFT:
c??      call four1(data,n/2,+1)
         call four1(data,n/2,-1)
c: INCLUDE FACTOR OF 2/N NOW:
         fac=2./float(n)
         do 40 j=1,n
            data(j)=fac*data(j)
40       continue
      endif
c
      return
      end
ccc
      subroutine four1(data,nn,isign)
c
c: Replaces data(1:2*nn) by its discrete FFT, if isign=1, or replaces
c: data(1:2*nn) by nn times its inverse FFT, if isign=-1.  data is a
c: complex array of length nn or, equivalently, a real array of length
c: 2*nn.  nn MUST be an integer power of 2 (not checked).
c:
c: nn is the number of complex data points, data is the data array,
c: and isign (+1 or -1) is the sign of i in the exponential of:
c: H(n)=sum(k=0,N-1) {h(k) exp(2*pi*i*k*n/N)}.  This is Numerical Recipes'
c: definition of the forward FFT.  The inverse FFT is defined by:
c: h(k)=(1/N) sum(n=0,N-1) {H(n) exp(-2*pi*i*k*n/N)}.  When isign is set
c: to -1, the routine thus computes the inverse FFT, except that it does
c: not multiply by the normalizing factor (1/N).  
c:
c: In the time domain, data(1) is the real part of the h(0), data(2) is the
c: imaginary part of h(0), etc. (complex time series).  In the frequency
c: domain, the real and imaginary parts of the zero frequency component are
c: in data(1) and data(2), resp.; the smallest positive frequency 
c: component is in data(3) and data(4); the maximum frequency component 
c: (both pos and neg) is contained in data(nn+1) and data(nn+2); then
c: the frequencies travel positive on the negative frequency axis; and the
c: smallest (in abs value) negative frequency component is in data(2*nn-1)
c: and data(2*nn).
c
      implicit none
      integer*4 nn,n,i,j,isign,m,mmax,istep
      real*4 data(2*nn),tempr,tempi
      real*8 theta,wi,wpi,wpr,wr,wtemp
c
      n=2*nn
      j=1
      do 10 i=1,n,2
         if(j .gt. i) then
            tempr=data(j)
            tempi=data(j+1)
            data(j)=data(i)
            data(j+1)=data(i+1)
            data(i)=tempr
            data(i+1)=tempi
         endif
         m=n/2
5        if((m .ge. 2) .and. (j .gt. m)) then
            j=j-m
            m=m/2
            goto 5
         endif
         j=j+m
10    continue
      mmax=2
6     if(n .gt. mmax) then
         istep=2*mmax
         theta=6.28318530717959d0/(isign*mmax)
         wpr=-2.d0*sin(0.5d0*theta)**2
         wpi=sin(theta)
         wr=1.d0
         wi=0.d0
         do 20 m=1,mmax,2
            do 30 i=m,n,istep
               j=i+mmax
               tempr=sngl(wr)*data(j) - sngl(wi)*data(j+1)
               tempi=sngl(wr)*data(j+1) + sngl(wi)*data(j)
               data(j)=data(i) - tempr
               data(j+1)=data(i+1) - tempi
               data(i)=data(i) + tempr
               data(i+1)=data(i+1) + tempi
30          continue
            wtemp=wr
            wr=wr*wpr - wi*wpi + wr
            wi=wi*wpr + wtemp*wpi + wi
20       continue
         mmax=istep
         goto 6
      endif
c
      return
      end
