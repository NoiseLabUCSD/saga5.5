      subroutine h_space(ii,isp,xkh,xkhsq,xksq,eta,etasq,gami,
     .   w,Vmatx,ailay,zetalay,ndv,iiw,jj1,jj2,ihf,xi)
c
c: Sets the plane wave reflection coefficients for a halfspace to zero.
c
      use scairy_com
      integer ii,isp,ndv,iiw,jj1,jj2,ihf,j
      complex*16 xkh,xkhsq,xksq,eta,etasq,gami(3),Vmatx(3,5),ailay(2),
     .   zetalay,zzero,xi,ai1,aip1,zzexp1,
     .   F,G,num,den,xi_k,xi_w,F_dot,G_dot,hdensq,eta_w
      real*8 w,magsq
      data zzero/(0.d0,0.d0)/
c
      if(isp .eq. 0) then
c: Airy halfspace:
         xi=(xkhsq - xksq)/etasq
         call airy_only (xi,ai1,aip1,zzexp1)
cc       call clairy(xi,1,ai1,bi1,aip1,bip1,zzexp1)
cc    if(magsq(ai1x*exp(-zzexp1x)-ai1*exp(-zzexp1)) .gt. 1.d-16) then
cc       print *,'airy_only problem: ',ai1x*exp(zzexp1x),
cc   .      ai1*exp(zzexp1)
cc    endif
c     print *,'clairy   : ',xi,ai1,aip1,zzexp1
c
         if(iiw .eq. 1) then
            ailay(1)=ai1
            ailay(2)=aip1
            zetalay=zzexp1
         endif
c
         F=gami(1)*ai1
         G=eta*aip1
         num=F + G
         den=F - G
         Vmatx(1,jj1)=num/den
         Vmatx(1,jj2)=zzero
         if(ndv .gt. 1) then
            xi_k=2.d0*xkh/etasq
            F_dot=gami(1)*aip1*xi_k + gami(2)*ai1
            G_dot=eta*xi*ai1*xi_k
            hdensq=0.5d0*den*den
            Vmatx(2,jj1)=(F*G_dot - G*F_dot)/hdensq
            Vmatx(2,jj2)=zzero
            if(ndv .gt. 2) then
               if(ihf .eq. 0) then
                  eta_w=eta/(1.5d0*w)
                  xi_w=(xksq + 2.d0*xkhsq)/(-1.5d0*w*etasq)
               else
c: Frequency-dependent gradient in halfspace:
                  eta_w=eta/w
                  xi_w=-2.d0*xkhsq/(w*etasq)
               endif
c: More general xi_w for different eta_w:
cc             xi_w=-2.d0*(xksq/w + eta*eta_w*xi)/etasq
               F_dot=gami(1)*aip1*xi_w + gami(3)*ai1
               G_dot=eta*xi*ai1*xi_w + eta_w*aip1
               Vmatx(3,jj1)=(F*G_dot - G*F_dot)/hdensq
               Vmatx(3,jj2)=zzero
            endif
         endif
      else
c: Homogeneous halfspace:
         do j=1,ndv
            Vmatx(j,jj1)=zzero
            Vmatx(j,jj2)=zzero
         enddo
      endif
      Vmatx(1,5)=zzero
c
      return
      end
