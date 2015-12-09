      subroutine rtc_calc(a,nx_p2,nx_s1,nx_s2,mrat,cp1,cs1,
     .   cs2,nd,aturn,rc12,tc12,nturn)
c
c: Calculates the p-p reflection coefficient rc12
c: and the p-p transmission coefficient tc12 for the 
c: plane wave incident at the 1-2 layer interface. 
c: Plane wave angle is given by a=cos(theta)/c, where theta
c: is the grazing angle.
c: Rewritten 3-19-93 to handle attenuation correctly.
      implicit integer*4(i-n)
      complex nx_p2(44),nx_s1(44),nx_s2(44)
      real a,cp1(44),cs1(44),cs2(44),aturn(44),mrat(44)
      complex rc12(44),tc12(44)
      integer*4 nturn(44)
      complex l1,l2,l3,l4,m1,m2,m3,m4,l2l4,m2m4,d,q,qpr2,r1mq
      complex qp1,qp2,qs1,qs2,psip1,psip2,psis1,psis2
c
      do 10 j=1,nd
         nturn(j)=0
         if(j .eq. 1 .or. j .eq. nd) then
            if(a .gt. aturn(j)) then
               rc12(j)=cmplx(0.,-1.)
               tc12(j)=cmplx(0.,0.)
               nturn(j)=1
               goto 10
            endif
         endif
c: FLUID OVER FLUID:
         if(cs1(j) .eq. 0. .and. cs2(j) .eq. 0.) then
               costh1=cp1(j)*a
               costh1sq=costh1**2
               qp1=sqrt(cmplx(1. - costh1sq,0.))
               qp2=sqrt(nx_p2(j) - costh1sq)
               qp1=dcmplx(real(qp1),abs(aimag(qp1)))
               qp2=dcmplx(real(qp2),abs(aimag(qp2)))
               psip1=qp1/costh1
               psip2=qp2/costh1
               l1=psip1*mrat(j)
               d=l1 + psip2
               rc12(j)=(l1 - psip2)/d
               tc12(j)=1. + rc12(j)
20          continue
c: FLUID OVER SOLID:
         elseif(cs1(j) .eq. 0. .and. cs2(j) .ne. 0.) then
            cs2sq=cs2(j)**2
            cs24=cs2sq**2
            costh1=cp1(j)*a
            costh1sq=costh1**2
            qp1=sqrt(cmplx(1. - costh1sq,0.))
            qp2=sqrt(nx_p2(j) - costh1sq)
            qp1=cmplx(real(qp1),abs(aimag(qp1)))
            qp2=cmplx(real(qp2),abs(aimag(qp2)))
            psip1=qp1/costh1
            psip2=qp2/costh1
            qs2=sqrt(nx_s2(j) - costh1sq)
            qs2=cmplx(real(qs2),abs(aimag(qs2)))
            psis2=qs2/costh1
            m1=(a**2)*(psis2**2 - 1.)
            l1=psip1*(m1**2)
            l2=4.*psip1*psip2*psis2*a**4
            l3=psip2/(cs24*mrat(j))
            d=l1 + l2 + l3
            rc12(j)=(l1+l2-l3)/d
            tc12(j)=2.*m1*psip1/(d*cs2sq)
c: SOLID OVER FLUID:
         elseif(cs1(j) .ne. 0. .and. cs2(j) .eq. 0.) then
            cs1sq=cs1(j)**2
            cs14=cs1sq**2
            costh1=cp1(j)*a
            costh1sq=costh1**2
            qp1=sqrt(cmplx(1. - costh1sq,0.))
            qp2=sqrt(nx_p2(j) - costh1sq)
            qp1=cmplx(real(qp1),abs(aimag(qp1)))
            qp2=cmplx(real(qp2),abs(aimag(qp2)))
            psip1=qp1/costh1
            psip2=qp2/costh1
            qs1=sqrt(nx_s1(j) - costh1sq)
            qs1=cmplx(real(qs1),abs(aimag(qs1)))
            psis1=qs1/costh1
            m1=(a**2)*(psis1**2 - 1.)
            l1=4.*psip1*psip2*psis1*a**4
            l2=mrat(j)*psip1/cs14
            l3=psip2*(m1**2)
            d=l1+l2 +l3
            rc12(j)=(l1+l2-l3)/d
            tc12(j)=2.*cs1sq*l2*m1/d
c: SOLID OVER SOLID:
         else
            cs1sq=cs1(j)**2
            cs2sq=cs2(j)**2
               costh1=cp1(j)*a
               costh1sq=costh1**2
               qp1=sqrt(cmplx(1. - costh1sq,0.))
               qp2=sqrt(nx_p2(j) - costh1sq)
               qp1=cmplx(real(qp1),abs(aimag(qp1)))
               qp2=cmplx(real(qp2),abs(aimag(qp2)))
               psip1=qp1/costh1
               psip2=qp2/costh1
               qs1=sqrt(nx_s1(j) - costh1sq)
               qs1=cmplx(real(qs1),abs(aimag(qs1)))
               psis1=qs1/costh1
               qs2=sqrt(nx_s2(j) - costh1sq)
               qs2=cmplx(real(qs2),abs(aimag(qs2)))
               psis2=qs2/costh1
c
               q=(2.*(cs1sq - mrat(j)*cs2sq))*a**2
               qpr2=q + mrat(j) 
               r1mq=1. - q
               l1=qpr2*psip1
               l2=1. - mrat(j) - q
               l3=r1mq*psip2
               l4=-q*psip2*psis1
               m1=q*psis2*psip1
               m2=r1mq*psis2
               m3=-l2
               m4=qpr2*psis1
               l2l4=l2 + l4
               m2m4=m2 + m4
               d=(l1+l3)*m2m4 - l2l4*(m1+m3)
               rc12(j)=((l1-l3)*m2m4 - l2l4*(m1-m3))/d
               tc12(j)=2.*mrat(j)*psip1*m2m4/d
35          continue
         endif
10    continue
cxx   print *,'end rtc_calc: ',
cxx  .   20.*log10(max(1.e-6,cabs(rc12(1:nd)))),
cxx  .   20.*log10(max(1.e-6,cabs(tc12(1:nd)))),nturn(1:nd)
c
      return
      end
