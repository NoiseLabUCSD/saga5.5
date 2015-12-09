      subroutine rtc_multi(a,rtcmag,rtcph,nt)
c: Calculates the total reflection/transmission coefficient
c: (rtc) for the ray with the snell invariant a and the path
c: specified by ntop, nbotx, nttot, and ndp.
c
      implicit integer*4(i-n)
      include 'common/pathway'
      include 'common/cbotcom'
      include 'common/pascal'
      complex rc12(44),tc12(44),rc21(44),tc21(44),rtc
      complex rr,tt,spr,spr_term,rrtt
      integer*4 nturn(44)
c
10    continue
      nd2=ndp+2
      call rtc_calc(a,nx_p2,nx_s1,nx_s2,mratx,cp1x,cs1x,
     .   cs2x,nd2,aturn,rc12,tc12,nturn)
      call rtc_calc(a,ny_p1,ny_s2,ny_s1,mraty,cp2x,cs2x,
     .   cs1x,nd2,aturn,rc21,tc21,nturn)
c
      nt=nturn(1)*ntop + nturn(nd2)*nttot(ndp)/2
cxx   print *,'rtc_multi: ',ntop,nbotx,nttot(1:ndp)/2
      if(ntop .gt. 0) then
         rtc=rc12(1)**ntop
      else
         rtc=cmplx(1.,0.)
      endif
      if(nbotx .gt. 0) then
         if(ndp .eq. 0) then
            rtc=rtc*rc12(2)**nbotx
ccc      elseif(kmult .eq. 0) then
ccc         do 25 j=2,nd2-1
ccc            rtc=rtc*tc12(j)*tc21(j)
25          continue
ccc         rtc=rtc*rc12(nd2)
         else
            nb1=nbotx
            do 30 j=2,nd2-1
               nb11=nb1+1
               nb2=nttot(j-1)/2
               nb=min(nb1,nb2)
               spr=cmplx(0.,0.)
               rr=rc12(j)*rc21(j)
               tt=tc12(j)*tc21(j)
c: Break off last term in case rr=0 (can't take 0.**0):
               do 40 jb=1,nb-1
c: EKW BUG 11-30-95: when the float is not included, the product here
c: turned out to negative, despite two positive integer factors.
c: Something crazy going on with integer*8 perhaps.
                  p_prod=float(p_mat(nb11,jb+1))*float(p_mat(nb2,jb))
cc    print *,'p_mat = ',jb,nb,p_mat(nb11,jb+1),p_mat(nb2,jb),p_prod
                  rrtt=(rr**(nb-jb))*(tt**jb)
                  spr_term=p_prod*rrtt
cc    print *,'rrtt = ',abs(rrtt),abs(spr_term)
                  spr=spr + spr_term
40             continue
               spr=spr + p_mat(nb11,nb+1)*p_mat(nb2,nb)*(tt**nb)
               if(nb1 .gt. nb2) then
                  spr=spr*rc12(j)**(nb1-nb2)
               elseif(nb2 .gt. nb1) then
                  spr=spr*rc21(j)**(nb2-nb1)
               endif
               rtc=rtc*spr
29             nb1=nb2
30          continue
            rtc=rtc*rc12(nd2)**nb1
         endif
      endif
      rtcmag=abs(rtc)
cc    if(rtcmag .gt. 1) then
cc       print *,'bad rtcmag: ',rtcmag
cc       goto 10
cc    endif
      if(rtcmag .gt. 0.) then
         rtcph=atan2(aimag(rtc),real(rtc))
      else
         rtcph=0.
      endif
cxx   print *,'rtc,rtcmag,rtcph = ',rtc,20.*log10(amax1
cxx  .   (1.e-6,rtcmag)),rtcph
c
      return
      end
