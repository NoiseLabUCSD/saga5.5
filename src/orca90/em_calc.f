      subroutine em_calc(xkhsq,xk1sq,xk2sq,eta,etasq,gami,iso,zmr,zmi,
     .   nz,philay,dphilay,Aplay,phix,dphix,exp_gbs,xi,ai,aip,bi,bip,
     .   zzexp,inc,ailay,bilay,zetalay,aisoln,ii1,ii2,jlay,jjhsp,kw0)
c
c: Computes the propagator matrix and its derivative for a given 
c: 1/c^2 linear layer: see Notebook92, p. 69-70.
c: xkhsq=k^2=(w*cos(th)/c0)^2; xk1sq=(w/c1)^2; 
c: etasq=(w^2*beta)^(2/3).
c 
      implicit none
ccx   use scairy_com
      integer*4 iso,nz,inc,aisoln,ii1,ii2,jlay,jjhsp,jz,iir1,iir2,jsign
      complex*16 xkhsq,xk1sq,xk2sq,eta,etasq,gami,philay(2),
     .   dphilay(2),phix(nz),dphix(nz),edif,zeta1,zetb1,A1exp,B1exp
      complex*16 xi(nz),ai(nz),aip(nz),bi(nz),
     .   bip(nz),zzexp(nz),det,A1,B1,A1fac,B1fac,xi0
      complex*16 ailay(2,2),bilay(2,2),zetalay(2)
      real*8 zmr(nz),zmi(nz),Aplay(2),exp_gbs(nz),Aprefa,Aprefb,
     .   exlim,exmax,kw0
      logical stradflag
      data exlim/-50.d0/
c
c: Return if no depths required for this layer:
      if(nz .eq. 0) return
c
      if(iso .eq. 1) then
c: For isovelocity layers:
         if(jlay .eq. jjhsp) then
c: For isospeed halfspace, only allow to match at top of layer (set ii2=ii1):
            call a1b1_iso(philay,dphilay,Aplay,inc,ii1,ii1,det,gami,
     .         zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .         jlay,jjhsp)
         else
            call a1b1_iso(philay,dphilay,Aplay,inc,ii1,ii2,det,gami,
     .         zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .         jlay,jjhsp)
         endif
         do jz=1,nz
            exp_gbs(jz)=0.d0
            edif=gami*dcmplx(zmr(jz),zmi(jz))
            A1exp=Aprefa + (edif - zeta1)
            B1exp=Aprefb - (edif - zetb1)
czs: Take out exponential factor from A1fac and B1fac if necessary:
            if(zmi(jz) .ne. 0.d0) then
               exmax=max(dreal(A1exp),dreal(B1exp))
               A1exp=A1exp - exmax
               B1exp=B1exp - exmax
               exp_gbs(jz)=exmax
            endif
            if(dreal(A1exp) .gt. exlim) then
               if(dreal(B1exp) .gt. exlim) then
                  A1fac=A1*cdexp(A1exp)
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=A1fac + B1fac
                  dphix(jz)=inc*gami*(A1fac - B1fac)
               else
                  A1fac=A1*cdexp(A1exp)
                  phix(jz)=A1fac
                  dphix(jz)=inc*gami*A1fac
               endif
            else
               if(dreal(B1exp) .gt. exlim) then
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=B1fac
                  dphix(jz)=-inc*gami*B1fac
               else
                  phix(jz)=cmplx(0.,0.)
                  dphix(jz)=cmplx(0.,0.)
               endif
            endif
cxx         phix(jz)=A1fac + B1fac
cxx         dphix(jz)=inc*gami*(A1fac - B1fac)
         enddo
      else
c: For linear 1/c^2 profile, use Airy function solutions.
c: Compute Airy function arguments at top and bottom of layer:
         xi0=(xkhsq - xk1sq)/etasq
         do jz=1,nz
            xi(jz)=xi0 - eta*dcmplx(zmr(jz),zmi(jz))
         enddo
cxx      call scairy3(nz,xi,xi0,ai,aip,bi,bip,zzexp,det,-aisoln)
         call scairy3(nz,xi,xi0,ai,aip,bi,bip,zzexp,det,aisoln)
c: Compute A1,B1 coeffs of Ai and Aihat for phi in layer (see p. 87,82):
         if(jlay .eq. jjhsp) then
c: For Airy halfspace, only allow to match at top of layer (set ii2=ii1):
            B1=(0.d0,0.d0)
            A1=philay(ii1)/ailay(1,ii1)
            zeta1=zetalay(ii1)
            Aprefa=Aplay(ii1)
            Aprefb=-1.d300
c
c           call a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii1,eta,det,
c    .         ailay,bilay,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,
c    .         Aprefa,Aprefb,jlay,jjhsp)
cc          print *,'-eta = ',-eta,
cc   .         atan2(-dimag(eta),-dreal(eta))*180/acos(-1.)
c           write(66,100) dreal(cdsqrt(xkhsq))/kw0,
c    .         dimag(cdsqrt(xkhsq))*8685.9,xi0,xi(nz)
100         format(6(e12.5,2x))
         else
            call a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii2,eta,det,
     .         ailay,bilay,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,
     .         Aprefa,Aprefb,jlay,jjhsp)
         endif
cxx   print *,'bot: a1,b1 = ',a1*cdexp(zeta1),b1*cdexp(-zetb1),
cxx  .   iir1,iir2
         do jz=1,nz
c: EKW FIX:
            if(zmi(jz) .ne. 0.d0) then
cj3: check to see if xi(source) needs different solution for stability
               stradflag=(jsign(dimag(xi(jz))) .ne. jsign(dimag(xi0)))
               if(stradflag) then
                  A1exp=Aprefa+zeta1
                  B1exp=Aprefb-zetb1
                  call src_recalc(A1,B1,xi(jz),eta,xkhsq,xk1sq,xk2sq,
     .               etasq,A1exp,B1exp,ii1,ii2,exlim,phix(jz),dphix(jz),
     .               ai(jz),aip(jz),bi(jz),bip(jz),zzexp(jz),
     .               exp_gbs(jz),inc,jlay,jjhsp)
cj3: needed values calculated already by src_recalc, so exit and proceed
                  goto 101
               endif
            endif
c: Make ai,aip get factor of exp(-zeta1); bi,bip get exp(zetb1):
            A1exp=Aprefa + (zeta1 - zzexp(jz))
            B1exp=Aprefb - (zetb1 - zzexp(jz))
            if(zmi(jz) .ne. 0.d0) then
               exmax=max(dreal(A1exp),dreal(B1exp))
               A1exp=A1exp - exmax
               B1exp=B1exp - exmax
               exp_gbs(jz)=exmax
            endif
            if(dreal(A1exp) .gt. exlim) then
               if(dreal(B1exp) .gt. exlim) then
                  A1fac=A1*cdexp(A1exp)
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=A1fac*ai(jz) + B1fac*bi(jz)
                  dphix(jz)=-inc*eta*(A1fac*aip(jz) + B1fac*bip(jz))
               else
                  A1fac=A1*cdexp(A1exp)
                  phix(jz)=A1fac*ai(jz)
                  dphix(jz)=-inc*eta*(A1fac*aip(jz))
               endif
            else
               if(dreal(B1exp) .gt. exlim) then
                  B1fac=B1*cdexp(B1exp)
                  phix(jz)=B1fac*bi(jz)
                  dphix(jz)=-inc*eta*(B1fac*bip(jz))
               else
                  phix(jz)=cmplx(0.,0.)
                  dphix(jz)=cmplx(0.,0.)
               endif
            endif
cxx         phix(jz)=A1fac*ai(jz) + B1fac*bi(jz)
cxx         dphix(jz)=-inc*eta*(A1fac*aip(jz) + B1fac*bip(jz))
101      enddo
      endif
99    continue
c
      return
      end
ccc
      subroutine a1b1_iso(philay,dphilay,Aplay,inc,ii1,ii2,det,
     .   gami,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .   jlay,jjhsp)
c
      implicit none
      integer*4 inc,ii1,ii2,iir1,iir2,jlay,jjhsp,ii
      complex*16 philay(2),dphilay(2),det,gami,zetalay(2),A1,B1,
     .   zeta1,zetb1,detinv,term1,term2,term1x,term2x,round_check,
     .   A1x,B1x
      real*8 Aplay(2),Aprefa,Aprefb,rat1,rat2
c
      detinv=0.5d0/gami
      term1=gami*philay(ii2)
      term2=inc*dphilay(ii2)
      A1=round_check(term1,term2,detinv,iir1,rat1)
      ii=ii2
      if(iir1 .eq. 1) then
         A1x=A1
         term1x=gami*philay(ii1)
         term2x=inc*dphilay(ii1)
         A1=round_check(term1x,term2x,detinv,iir1,rat2)
         ii=ii1
         if(iir1 .eq. 1) then
cc          print *,'A1 iso unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               A1=A1x
               ii=ii2
            endif
         endif
      endif
      zeta1=zetalay(ii)
      Aprefa=Aplay(ii)
c
      if(jlay .eq. jjhsp) then
c: For halfspace, B1=0:
         B1=dcmplx(0.d0,0.d0)
         zetb1=zetalay(ii)
         Aprefb=-1.d100
         return
      endif
c
      term2=-term2
      B1=round_check(term1,term2,detinv,iir2,rat1)
      ii=ii2
      if(iir2 .eq. 1) then
         B1x=B1
         term1x=gami*philay(ii1)
         term2x=-inc*dphilay(ii1)
         B1=round_check(term1x,term2x,detinv,iir2,rat2)
         ii=ii1
         if(iir2 .eq. 1) then
cc          print *,'B1 iso unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               B1=B1x
               ii=ii2
            endif
         endif
      endif
      zetb1=zetalay(ii)
      Aprefb=Aplay(ii)
c
      return
      end
ccc
      subroutine a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii2,eta,det,
     .   ailay,bilay,zetalay,A1,B1,iir1,iir2,zeta1,zetb1,Aprefa,Aprefb,
     .   jlay,jjhsp)
c
      implicit none
      integer*4 inc,ii1,ii2,iir1,iir2,ii,jlay,jjhsp
      complex*16 philay(2),dphilay(2),eta,det,ailay(2,2),bilay(2,2),
     .   zetalay(2),A1,B1,zeta1,zetb1,dphfac,dphfac2,detinv,
     .   term1,term2,round_check,A1x,B1x
cc    complex*16 A1exp,B1exp,A1fac,B1fac,ph,dph
      real*8 Aplay(2),Aprefa,Aprefb,rat1,rat2
c
      detinv=1.d0/det
      dphfac=inc*dphilay(ii2)/eta
      term1=bilay(2,ii2)*philay(ii2)
      term2=bilay(1,ii2)*dphfac
      A1=round_check(term1,term2,detinv,iir1,rat1)
cxx   A1=detinv*(bilay(2,ii2)*philay(ii2) + bilay(1,ii2)*dphfac)
      ii=ii2
      if(iir1 .eq. 1) then
         A1x=A1
         dphfac2=inc*dphilay(ii1)/eta
         term1=bilay(2,ii1)*philay(ii1)
         term2=bilay(1,ii1)*dphfac2
         A1=round_check(term1,term2,detinv,iir1,rat2)
         ii=ii1
         if(iir1 .eq. 1) then
cc          print *,'A1 unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               A1=A1x
               ii=ii2
            endif
         endif
      endif
      Aprefa=Aplay(ii)
      zeta1=zetalay(ii)
c
      if(jlay .eq. jjhsp) then
         B1=dcmplx(0.d0,0.d0)
         Aprefb=-1.d100
         zetb1=zetalay(ii)
         return
      endif
c
      term1=-ailay(2,ii2)*philay(ii2)
      term2=-ailay(1,ii2)*dphfac
      B1=round_check(term1,term2,detinv,iir2,rat1)
cxx   B1=detinv*(-ailay(2,ii2)*philay(ii2) - ailay(1,ii2)*dphfac)
      ii=ii2
      if(iir2 .eq. 1) then
         B1x=B1
         dphfac2=inc*dphilay(ii1)/eta
         term1=-ailay(2,ii1)*philay(ii1)
         term2=-ailay(1,ii1)*dphfac2
         B1=round_check(term1,term2,detinv,iir2,rat2)
         ii=ii1
         if(iir2 .eq. 1) then
cc          print *,'B1 unstable at top & bot: ',rat1,rat2
            if(rat1 .lt. rat2) then
               B1=B1x
               ii=ii2
            endif
         endif
      endif
      Aprefb=Aplay(ii)
      zetb1=zetalay(ii)
c
cc       do ii=1,2
cc          A1exp=Aprefa + (zeta1 - zetalay(ii))
cc          B1exp=Aprefb - (zetb1 - zetalay(ii))
cc          A1fac=A1*cdexp(A1exp)
cc          B1fac=B1*cdexp(B1exp)
cc          ph=A1fac*ailay(1,ii) + B1fac*bilay(1,ii)
cc          dph=-inc*eta*(A1fac*ailay(2,ii) + B1fac*bilay(2,ii))
cc    print *,'a1b1_airy: ii,phi = ',ii,ph,philay(ii)*dexp(Aplay(ii))
cc    print *,'   dphi = ',ii,dph,dphilay(ii)*dexp(Aplay(ii))
cc       enddo
c
      return
      end
ccc
      function round_check(term1,term2,deninv,iir,ratio)
c
c: Checks to see if cancellation occurs in term1+term2.
c: Returns round_check=(term1 + term2)*deninv if not.
      implicit none
      integer*4 iir
      complex*16 round_check,term1,term2,deninv,numerator
      real*8 magsq,magsum,ratio
c
      numerator=term1 + term2
      magsum=magsq(term1) + magsq(term2)
      ratio=magsq(numerator)/magsum
cxx   if(ratio .gt. 1.d-18) then
      if(ratio .gt. 1.d-14) then
         round_check=numerator*deninv
         iir=0
      else
         iir=1
         round_check=(0.,0.)
cxx      print *,'Round to zero: ',term1,term2,numerator,deninv
      endif
c
      return
      end
ccc
      function jsign(number)
c
      implicit none
      real*8 number
      integer*4 jsign
c
      if (number .gt. 0d0) then
         jsign=1
         return
      else
         if(number .lt. 0d0) then
            jsign=-1
            return
         else
            jsign=0
            return
         endif
      endif
c
      return
      end
ccc
      subroutine src_recalc(A1,B1,xisrc,eta,xkhsq,xk1sq,xk2sq,etasq,
     . A1exp,B1exp,ii1,ii2,exlim,phisrc,dphisrc,aisrc,aipsrc,bisrc,
     . bipsrc,zzsrc,srcexp,inc,jlay,jjhsp)
c
c: Subroutine to recalculate source excitations if there a problem with    ::
c: xi(source,receivers) on opposite sides of real axis.  (part of patch j3)::
c
      implicit none
ccx   use scairy_com
c: INPUT                ::
      complex*16 A1,B1,xisrc,eta,xkhsq,xk1sq,xk2sq,etasq,A1exp,
     .   B1exp
      integer*4  ii1,ii2,inc,jlay,jjhsp
      real*8 exlim
c: OUTPUT       ::
      complex*16 phisrc,dphisrc,aisrc,aipsrc,bisrc,bipsrc,zzsrc
      real*8  srcexp
c: OTHER                ::
      complex*16 Gii1,Gii2,ailay(2,2),bilay(2,2),A1f,B1f,philay(2),
     .   dphilay(2),A2,B2,A2exp,B2exp,A2fac,B2fac,xii1,xii2,det,
     .   detii1,detii2,zetalay(2),zeta2,zetb2,A1fac,B1fac
      real*8 eps,Aplay(2),exmax,Aprefa2,Aprefb2
      integer*4 iibad,jsign,iir1,iir2,ii,aisoln
c
c:  Calculate old xis, and zero out imag. part::
      Gii1=(xkhsq-xk1sq)/etasq
      Gii2=(xkhsq-xk2sq)/etasq
      eps=1d-20 
      xii1=cmplx(dreal(Gii1),-eps*jsign(dimag(xisrc)))
      xii2=cmplx(dreal(Gii2),-eps*jsign(dimag(xisrc)))
c:  First calculate phi, dphi on the same side of real axis as orig layer depths::
      call scairy2(xii1,ailay(1,ii1),ailay(2,ii1),bilay(1,ii1),
     .   bilay(2,ii1),zetalay(ii1),detii1,xii2,ailay(1,ii2),
     .   ailay(2,ii2),bilay(1,ii2),bilay(2,ii2),
     .   zetalay(ii2),detii2,aisoln,iibad)
      A1fac=A1*cdexp(A1exp) 
      B1fac=B1*cdexp(B1exp)
      do ii=1,2
         A1f=A1fac*cdexp(-zetalay(ii))
         B1f=B1fac*cdexp(zetalay(ii))
         philay(ii)=(A1f*ailay(1,ii)+B1f*bilay(1,ii))
         dphilay(ii)=-inc*eta*(A1f*ailay(2,ii)+B1f*bilay(2,ii))
         Aplay(ii)=0.d0
      enddo
c:  Now get ailay, etc... on other side of line ::
      xii1=cmplx(dreal(xii1),eps*jsign(dimag(xisrc)))
      xii2=cmplx(dreal(xii2),eps*jsign(dimag(xisrc)))
      call scairy2(xii1,ailay(1,ii1),ailay(2,ii1),bilay(1,ii1),
     .   bilay(2,ii1),zetalay(ii1),detii1,xii2,ailay(1,ii2),
     .   ailay(2,ii2),bilay(1,ii2),bilay(2,ii2),
     .   zetalay(ii2),detii2,aisoln,iibad)
      call scairy3(1,xisrc,xisrc,aisrc,aipsrc,bisrc,bipsrc,zzsrc,det,
     .   aisoln)
      call a1b1_airy(philay,dphilay,Aplay,inc,ii1,ii2,eta,det,ailay,
     .     bilay,zetalay,A2,B2,iir1,iir2,zeta2,zetb2,Aprefa2,
     .     Aprefb2,jlay,jjhsp)
      A2exp=Aprefa2+(zeta2-zzsrc)
      B2exp=Aprefb2-(zetb2-zzsrc)
      exmax=max(dreal(A2exp),dreal(B2exp))
      A2exp=A2exp-exmax
      B2exp=B2exp-exmax
      srcexp=exmax
      if(dreal(A2exp) .gt. exlim) then
         if(dreal(B2exp) .gt. exlim) then
            A2fac=A2*cdexp(A2exp)
            B2fac=B2*cdexp(B2exp)
            phisrc=A2fac*aisrc + B2fac*bisrc
            dphisrc=-inc*eta*(A2fac*aipsrc + B2fac*bipsrc)
         else
            A2fac=A2*cdexp(A2exp)
            phisrc=A2fac*aisrc
            dphisrc=-inc*eta*(A2fac*aipsrc)
         endif
      else
         if(dreal(B2exp) .gt. exlim) then
            B2fac=B2*cdexp(B2exp)
            phisrc=B2fac*bisrc
            dphisrc=-inc*eta*(B2fac*bipsrc)
         else
            phisrc=cmplx(0.,0.)
            dphisrc=cmplx(0.,0.)
        endif
      endif
c
      return
      end
