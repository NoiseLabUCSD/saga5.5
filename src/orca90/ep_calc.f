      subroutine ep_calc(xkh,w,xkhsq,xk1sq,xk2sq,eta,etasq,gami,h,
     .   e11,e12,e21,e22,zzexp,iso,ihf,isx,ndv,iiw,ailay,bilay,
     .   zetalay,aisoln,ii1,ii2)
c
c: Computes the propagator matrix and its derivative for a given 
c: 1/c^2 linear layer: see Notebook92, p. 69-70.
c: xkh=horizontal wavenumber; xkqsq=xkh^2; (xk1sq=(w/c1)^2; 
c: xk2sq=(w/c2)^2; eta=(w^2*beta)^(1/3); etasq=eta^2; beta=gradient;
c: gami=i*sqrt(xk1sq - xkhsq); h=layer thickness.
c: E and E'=dE/dk must be multiplied by exp(zzexp) to get true values.
c: For ndv=2, the derivative w.r.t. k is computed.
c: For ndv=3, the derivative w.r.t. w is computed.
c: ihf=1 means h is inversely proportional to frequency.
c 
      implicit none
cxx   use scairy_com
      integer*4 iso,ihf,isx,ndv,j,iiw,iibad,aisoln,ii1,ii2
      complex*16 xkh,xkhsq,xk1sq,xk2sq,eta,etasq,gami(3),
     .   e11(3),e12(3),e21(3),e22(3),zzexp
      complex*16 xi1,xi2,ai1,aip1,bi1,bip1,ai2,aip2,bi2,bip2,
     .   det1,det2,zfac,zzexp1,zzexp2,igamhp
      complex*16 ailay(2,2),bilay(2,2),zetalay(2),gamih
      real*8 h,w,pie
      data pie/3.14159265358979/
c
      if(h .eq. 0.d0) return
      if(iso .eq. 1) then
         gamih=gami(1)*h
c: For isovelocity layers, combine terms that will be needed in
c: rp_slay, rp_slay0, or rp_flay:
         if(isx .eq. 1) then
c: For isx=1, set e11=e11 + gami*e12, e12=e11 - gami*e12,
c: e21=e21 + gami*e22, e22=e21 - gami*e22 (see rp_slay,rp_slay0):
            if(dreal(gami(1)) .gt. 0.d0) then
c: Remove zzexp=exp(gami*h):
               zzexp=gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=(1.d0,0.d0)
               e12(1)=zfac
               e21(1)=gami(1)
               e22(1)=-gami(1)*zfac
               if(ndv .ge. 2) then
                  j=2
                  e11(j)=gami(j)*h
                  e12(j)=-e11(j)*zfac
                  e21(j)=gami(j)*(zzexp + 1.d0)
                  e22(j)=gami(j)*(zzexp - 1.d0)*zfac
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e11(j)=gami(j)*h
                     e12(j)=-e11(j)*zfac
                     e21(j)=gami(j)*(zzexp + 1.d0)
                     e22(j)=gami(j)*(zzexp - 1.d0)*zfac
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=igamhp
                     e12(j)=-igamhp*zfac
                     e21(j)=gami(1)*igamhp + gami(j)
                     e22(j)=(gami(1)*igamhp - gami(j))*zfac
                  endif
               endif
            else
c: Remove zzexp=exp(-gami*h):
               zzexp=-gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=zfac
               e12(1)=(1.d0,0.d0)
               e21(1)=gami(1)*zfac
               e22(1)=-gami(1)
               if(ndv .ge. 2) then
                  j=2
                  e12(j)=-gami(j)*h
                  e11(j)=-e12(j)*zfac
                  e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                  e22(j)=-gami(j)*(zzexp + 1.d0)
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e12(j)=-gami(j)*h
                     e11(j)=-e12(j)*zfac
                     e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                     e22(j)=-gami(j)*(zzexp + 1.d0)
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=igamhp*zfac
                     e12(j)=-igamhp
                     e21(j)=(gami(1)*igamhp + gami(j))*zfac
                     e22(j)=gami(1)*igamhp - gami(j)
                  endif
               endif
            endif
         else
c: For isx=2, set e11=gami*e11 + e21, e21=gami*e11 - e21, 
c: e12=gami*e12 + e22, e22=gami*e12 - e22 (see rp_flay): 
            if(dreal(gami(1)) .gt. 0.d0) then
c: Remove zzexp=exp(gami*h):
               zzexp=gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=gami(1)
               e21(1)=gami(1)*zfac
               e12(1)=(1.d0,0.d0)
               e22(1)=-zfac
               if(ndv .ge. 2) then
                  j=2
                  e11(j)=gami(j)*(zzexp + 1.d0)
                  e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                  e12(j)=gami(j)*h
                  e22(j)=e12(j)*zfac
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e11(j)=gami(j)*(zzexp + 1.d0)
                     e21(j)=gami(j)*(-zzexp + 1.d0)*zfac
                     e12(j)=gami(j)*h
                     e22(j)=e12(j)*zfac
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=gami(1)*igamhp + gami(j)
                     e21(j)=(-gami(1)*igamhp + gami(j))*zfac
                     e12(j)=igamhp
                     e22(j)=igamhp*zfac
                  endif
               endif
            else
c: Remove zzexp=exp(-gami*h):
               zzexp=-gamih
               zfac=cdexp(-2.d0*zzexp)
               e11(1)=gami(1)*zfac
               e21(1)=gami(1)
               e12(1)=zfac
               e22(1)=(-1.d0,0.d0)
               if(ndv .ge. 2) then
                  j=2
                  e11(j)=gami(j)*(-zzexp + 1.d0)*zfac
                  e21(j)=gami(j)*(zzexp + 1.d0)
                  e22(j)=gami(j)*h
                  e12(j)=e22(j)*zfac
               endif
               if(ndv .ge. 3) then
                  j=3
                  if(ihf .eq. 0) then
                     e11(j)=gami(j)*(-zzexp + 1.d0)*zfac
                     e21(j)=gami(j)*(zzexp + 1.d0)
                     e22(j)=gami(j)*h
                     e12(j)=e22(j)*zfac
                  else
c: See ORCA I, p. 124 bottom:
                     igamhp=-xkhsq*h/(gami(1)*w)
                     e11(j)=(gami(1)*igamhp + gami(j))*zfac
                     e21(j)=-gami(1)*igamhp + gami(j)
                     e12(j)=igamhp*zfac
                     e22(j)=igamhp
                  endif
               endif
            endif
         endif
         if(iiw .eq. 1) then
            zetalay(ii1)=0.d0
            zetalay(ii2)=gamih
         endif
c
c: Compute propagator matrix for Airy layer:
      else
c: For linear 1/c^2 profile, use Airy function solutions.
c: Compute Airy function arguments at top and bottom of layer:
         xi1=(xkhsq - xk1sq)/etasq
         xi2=(xkhsq - xk2sq)/etasq
c: Compute numerically stable Airy function solutions and derivatives:
c: Note: bi,bip are actually Ai(xi*eim23) and (d/dxi)(Ai(xi*eim23)) for
c: Im(xi)>=0, and Ai(xi*ei23) and (d/dxi)(Ai(xi*ei23)) for Im(xi)<0.
         call scairy2(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,ai2,aip2,
     .      bi2,bip2,zzexp2,det2,aisoln,iibad)
c
         if(iiw .eq. 1) then
            ailay(1,ii1)=ai1
            ailay(1,ii2)=ai2
            ailay(2,ii1)=aip1
            ailay(2,ii2)=aip2
            bilay(1,ii1)=bi1
            bilay(1,ii2)=bi2
            bilay(2,ii1)=bip1
            bilay(2,ii2)=bip2
            zetalay(ii1)=zzexp1
            zetalay(ii2)=zzexp2
         endif
c
         if(iibad .eq. 0) then
c: Compute propagator matrix for Airy layer:
            call airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,ai2,
     .         aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,zzexp,
     .         xkh,xkhsq,xk1sq,xk2sq,etasq,w,ndv)
         else
            call ai_strad(xi1,xk1sq,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,
     .         xk2sq,ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,
     .         zzexp,xkh,xkhsq,etasq,w,ndv)
         endif
c
      endif
c
      return
      end
ccc
      subroutine airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,
     .   ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,zzexp,
     .   xkh,xkhsq,xk1sq,xk2sq,etasq,w,ndv)
c
      implicit none
      integer*4 ndv
      complex*16 xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,ai2,aip2,
     .   bi2,bip2,zzexp2,det2,eta,e11(3),e12(3),e21(3),e22(3),
     .   zzexp,xkh,xkhsq,xk1sq,xk2sq,etasq,dexp,zfac,e21_eta,eta_e12,
     .   dxi_dk,dxi_dw1,dxi_dw2,dxi_fac
      real*8 w,deta_eta
c
c: 1st terms of E must be multiplied by exp(dexp); 2nd terms by exp(-dexp):
      dexp=zzexp1 - zzexp2
      if(dreal(dexp) .gt. 0.d0) then
c: Exponential factor of first of two terms dominates:
         zzexp=dexp
         zfac=cdexp(-2.d0*dexp)
c: Remove old factor from and insert new factor to second terms:
         e11(1)=(ai2*bip1 - zfac*bi2*aip1)/det1
         e12(1)=(-ai2*bi1 + zfac*bi2*ai1)/(-det1*eta)
         e21(1)=-eta*(aip2*bip1 - zfac*bip2*aip1)/det1
         e22(1)=(-aip2*bi1 + zfac*bip2*ai1)/det1
      else
c: Exponential factor of second of two terms dominates:
         zzexp=-dexp
         zfac=cdexp(2.d0*dexp)
c: Remove old factor from and insert new factor to first terms:
         e11(1)=(zfac*ai2*bip1 - bi2*aip1)/det1
         e12(1)=(-zfac*ai2*bi1 + bi2*ai1)/(-det1*eta)
         e21(1)=-eta*(zfac*aip2*bip1 - bip2*aip1)/det1
         e22(1)=(-zfac*aip2*bi1 + bip2*ai1)/det1
      endif
c
      if(ndv .ge. 2) then
         e21_eta=e21(1)/eta
         eta_e12=eta*e12(1)
c: Compute derivative w.r.t. k of exponential factor to do removed:
c: (see p.101-102 of Notebook):
         dxi_dk=2.d0*xkh/etasq
         e11(2)=dxi_dk*(xi1*eta_e12 - e21_eta)
         e12(2)=dxi_dk*(e11(1) - e22(1))/eta
         e21(2)=dxi_dk*(xi1*e22(1) - xi2*e11(1))*eta
         e22(2)=dxi_dk*(e21_eta - xi2*eta_e12)
         if(ndv .ge. 3) then
c: Compute derivative w.r.t. w (see p.108 of Notebook):
            deta_eta=1.d0/(1.5d0*w)
cxx         dxi_fac=-2.d0/(3.d0*w*etasq)
            dxi_fac=-deta_eta/etasq
            dxi_dw1=dxi_fac*(xk1sq + 2.d0*xkhsq)
            dxi_dw2=dxi_fac*(xk2sq + 2.d0*xkhsq)
            e11(3)=(xi1*eta_e12*dxi_dw1 - e21_eta*dxi_dw2)
            e12(3)=(e11(1)*dxi_dw1 - e22(1)*dxi_dw2)/eta - 
     .         e12(1)*deta_eta
            e21(3)=(xi1*e22(1)*dxi_dw1 - xi2*e11(1)*dxi_dw2)*eta +
     .         e21(1)*deta_eta
            e22(3)=(e21_eta*dxi_dw1 - xi2*eta_e12*dxi_dw2)
         endif
      endif
c
      return
      end
c
      subroutine ai_strad(xi1,xk1sq,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,
     .   xk2sq,ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11,e12,e21,e22,zzexp,
     .   xkh,xkhsq,etasq,w,ndv)
c
c: Computes the propagator matrix for a layer when xi1 and xi2 straddle
c: the real axis so that different pairs of solutions were originally
c: computed in scairy2.  When the determinant is bad in scairy2,
c: ai_strad is called. 
c: The intersection from xi1 to xi2 with the real axis is found,
c: and the layer is split into two layers, each of which now has the
c: same solution in scairy2.  The propagator matrices are then multiplied
c: together.
c
      implicit none
      integer*4 ndv,iibad,j,aisoln
      real*8 w,rfac
      complex*16 xi1,xk1sq,ai1,aip1,bi1,bip1,zzexp1,det1,xi2,xk2sq,
     .   ai2,aip2,bi2,bip2,zzexp2,det2,eta,e11(3),e12(3),e21(3),e22(3),
     .   zzexp,xkh,xkhsq,etasq
      complex*16 xi01,xk01sq,ai01,aip01,bi01,bip01,zzexp01,det01,
     .   xi02,xk02sq,ai02,aip02,bi02,bip02,zzexp02,det02,f11(3),f12(3),
     .   f21(3),f22(3),g11(3),g12(3),g21(3),g22(3),zzexpf,zzexpg
c
c: Compute factor for depth at which xi crosses real axis:
      rfac=-dimag(xi1)/(dimag(xi2) - dimag(xi1))
c: Let xi01 be slightly on xi1's side of real axis at intersection point:
      xi01=dcmplx(dreal(xi1) + (dreal(xi2)-dreal(xi1))*rfac,
     .   dsign(1.d-20,dimag(xi1)))
c: Compute interpolated value for xksq:
      xk01sq=xk1sq + rfac*(xk2sq - xk1sq)
c: Get Airy functions for first layer from xi1 to xi01:
      call scairy2(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi01,ai01,aip01,
     .   bi01,bip01,zzexp01,det01,aisoln,iibad)
      if(iibad .eq. 1) print *,'iibad=1 in ai_strad!!'
c: Compute propagator matrix for first layer:
      call airy_prop(xi1,ai1,aip1,bi1,bip1,zzexp1,det1,xi01,
     .   ai01,aip01,bi01,bip01,zzexp01,det01,eta,f11,f12,f21,
     .   f22,zzexpf,xkh,xkhsq,xk1sq,xk01sq,etasq,w,ndv)
c
c: Let xi02 be slightly on xi2's side of real axis at intersection point:
      xi02=dconjg(xi01)
      xk02sq=xk01sq
c: Get Airy functions for first layer from xi02 to xi2:
      call scairy2(xi02,ai02,aip02,bi02,bip02,zzexp02,det02,xi2,
     .   ai2,aip2,bi2,bip2,zzexp2,det2,aisoln,iibad)
      if(iibad .eq. 1) print *,'iibad=1 in ai_strad!!'
c: Compute propagator matrix for first layer:
      call airy_prop(xi02,ai02,aip02,bi02,bip02,zzexp02,det02,
     .   xi2,ai2,aip2,bi2,bip2,zzexp2,det2,eta,g11,g12,g21,g22,
     .   zzexpg,xkh,xkhsq,xk02sq,xk2sq,etasq,w,ndv)
c
c: Multiply propagator matrices together (G*F):
      e11(1)=g11(1)*f11(1) + g12(1)*f21(1)
      e21(1)=g21(1)*f11(1) + g22(1)*f21(1)
      e12(1)=g11(1)*f12(1) + g12(1)*f22(1)
      e22(1)=g21(1)*f12(1) + g22(1)*f22(1)
      zzexp=zzexpf + zzexpg
c
c: Compute derivatives using the chain rule:
      do j=2,ndv
         e11(j)=g11(1)*f11(j) + g11(j)*f11(1) + g12(1)*f21(j) +
     .      g12(j)*f21(1)
         e21(j)=g21(1)*f11(j) + g21(j)*f11(1) + g22(1)*f21(j) +
     .      g22(j)*f21(1)
         e12(j)=g11(1)*f12(j) + g11(j)*f12(1) + g12(1)*f22(j) +
     .      g12(j)*f22(1)
         e22(j)=g21(1)*f12(j) + g21(j)*f12(1) + g22(1)*f22(j) +
     .      g22(j)*f22(1)
      enddo
c
      return
      end
