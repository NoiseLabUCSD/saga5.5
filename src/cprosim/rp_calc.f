      subroutine rp_calc(ii,R,ndv,iiw)
c
c: Computes the plane wave reflection coefficient for a series of
c: fluid and solid layers.
c
      implicit none
      include 'Parms_com'
      include 'i_o_cpro_com'
      include 'i_o_saga_com'
      include 'gen_cpro_com'
      include 'gen_saga_com'
c      include 'i_o_com'
c      include 'gen_com'
c: Local variables:
      integer ii,ns,jflu1,jflu2,inc,j,ndv,iiw,k,ii1,ii2,
     .   jsol1,jsol2,iiwww,j1
      complex*16 e11(3),e12(3),e21(3),e22(3),s11(3),s12(3),s21(3),
     .   s22(3),zexpe,zexps,R(3)
c
      ii1=ii
      ii2=3 - ii
      inc=jflu(ii,3)
      jflu1=jflu(ii,1)
      jflu2=jflu(ii,2)
      jsol1=jsol(ii,1)
      jsol2=jsol(ii,2)
      ns=jhsp(ii)
c
c: Compute reflection coefficient at top of halfspace (0 for homogeneous):
      call h_space(ii,isp(ns),xkh,xkhsq,xksq(ii1,ns),eta(ns),etasq(ns),
     .   gami(1,ii1,ns),w,Vmat(1,1,ns,ii),ailay(1,ii,1,ns),
     .   zetalay(ii,1,ns),ndv,iiw,1,2,ihf(ns),xi_hsp(ii,1))
      if(allf(ii) .eq. 0) then
         call h_space(ii,iss(ns),xkh,xkhsq,xbsq(ii1,ns),etb(ns),
     .      etbsq(ns),beti(1,ii1,ns),w,Vmat(1,1,ns,ii),
     .      ailay(1,ii,2,ns),zetalay(ii,2,ns),ndv,iiw,3,4,
     .      ihf(ns),xi_hsp(ii,2))
      endif
c: Test for duct at top of Airy halfspace:
      if(jflu1 .eq. ns) then
         do j=1,ndv
            R(j)=Vmat(j,1,ns,ii)
         enddo
         return
      endif
c
c: Loop over solid layers:
      do j=jsol2,jsol1,-inc
c: j1 points to layer to be propagated from:
         j1=j + inc
c: k points to interface:
         k=j + jflu(ii,4)
         iiwww=iiw*iiww(j1)
c: Compute 2x2 propagator matrices for compressional and shear profiles:
         call ep_calc(xkh,w,xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .      etasq(j),gami(1,ii2,j),h(j),e11,e12,e21,e22,zexpe,isp(j),
     .      ihf(j),1,ndv,iiwww,ailay(1,1,1,j),bilay(1,1,1,j),
     .      zetalay(1,1,j),aisoln(1,j),ii1,ii2)
         call ep_calc(xkh,w,xkhsq,xbsq(ii1,j),xbsq(ii2,j),inc*etb(j),
     .      etbsq(j),beti(1,ii2,j),h(j),s11,s12,s21,s22,zexps,iss(j),
     .      ihf(j),1,ndv,iiwww,ailay(1,1,2,j),bilay(1,1,2,j),
     .      zetalay(1,2,j),aisoln(2,j),ii1,ii2)
         if(mm(k) .eq. 1) then
c: Mismatch at bottom of layer.  Call rp_slay as usual:
            call rp_slay(Vmat(1,1,j1,ii),Pcon(1,k),Qcon(1,k),
     .         Ucon(1,k),Vcon(1,k),e11,e12,e21,e22,zexpe,s11,s12,
     .         s21,s22,zexps,gami(1,ii1,j),beti(1,ii1,j),
     .         gami(1,ii1,j1),beti(1,ii1,j1),isp(j),iss(j),
     .         Vmat(1,1,j,ii),ndv,Wmat(1,j1,ii),iiwww,rhorat(k))
         else
c: No mismatch at bottom of layer.  Call rp_nomm:
            call rp_nomm(Vmat(1,1,j1,ii),e11,e12,e21,e22,zexpe,
     .         s11,s12,s21,s22,zexps,gami(1,ii1,j),beti(1,ii1,j),
     .         gami(1,ii1,j1),beti(1,ii1,j1),isp(j),iss(j),
     .         Vmat(1,1,j,ii),ndv,Wmat(1,j1,ii),iiwww)
         endif
      enddo
c
      if(allf(ii) .eq. 1) then
c: Propagate across fluid-fluid halfspace interface to bottom of last 
c: (fluid) layer:
         iiwww=iiw*iiww(ns)
         j1=jsol2
         k=j1 + jflu(ii,4)
         call rp_flay(Vmat(1,1,ns,ii),e11,e12,e21,e22,zexpe,
     .      gami(1,ii2,jflu2),gami(1,ii1,ns),1,0.d0,rhorat(k),
     .      mm(k),iiwww,ndv,Vmat(1,1,j1,ii),Wmat(1,ns,ii),jjfail)
         if(jjfail .gt. 0) return
      else
c: Propagate across solid-fluid interface to bottom of first fluid layer:
         j1=jflu2
         iiwww=iiw*iiww(jsol1)
         call rp_sfint(Vmat(1,1,jsol1,ii),Alay(1,ii),Blay(1,ii),
     .      ikcon,rholay(ii),gami(1,ii2,j1),gami(1,ii1,jsol1),
     .      beti(1,ii1,jsol1),Vmat(1,1,j1,ii),ndv,Wmat(1,jsol1,ii),
     .      iiwww)
      endif
c
      do j=jflu2,jflu1+inc,-inc
         j1=j - inc
c: Propagate through fluid layers from bottom of j'th layer to top of
c: j1'th layer (the next layer toward the svp minimum from the j'th):
         k=j + jflu(ii,4) - inc
         iiwww=iiw*iiww(j)
c: Compute 2x2 propagator matrices for compressional and shear profiles:
         call ep_calc(xkh,w,xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .      etasq(j),gami(1,ii2,j),h(j),e11,e12,e21,e22,zexpe,isp(j),
     .      ihf(j),2,ndv,iiw,ailay(1,1,1,j),bilay(1,1,1,j),
     .      zetalay(1,1,j),aisoln(1,j),ii1,ii2)
         call rp_flay(Vmat(1,1,j,ii),e11,e12,e21,e22,zexpe,
     .      gami(1,ii2,j1),gami(1,ii2,j),isp(j),h(j),rhorat(k),
     .      mm(k),iiwww,ndv,Vmat(1,1,j1,ii),Wmat(1,j,ii),jjfail)
         if(jjfail .gt. 0) return
      enddo
c
c: Do last layer (with no interface) if we are doing bottom half and ref depth
c: at top of layer or if we are doing top half and ref depth is at bottom:
      if(isvmin .eq. ii) then
         j=jflu1
         j1=j - inc
         iiwww=iiw*iiww(j)
         call ep_calc(xkh,w,xkhsq,xksq(ii1,j),xksq(ii2,j),inc*eta(j),
     .      etasq(j),gami(1,ii2,j),h(j),e11,e12,e21,e22,zexpe,isp(j),
     .      ihf(j),2,ndv,iiw,ailay(1,1,1,j),bilay(1,1,1,j),
     .      zetalay(1,1,j),aisoln(1,j),ii1,ii2)
         call rp_flay(Vmat(1,1,j,ii),e11,e12,e21,e22,zexpe,
     .      gami(1,ii1,j),gami(1,ii2,j),isp(j),h(j),1.d0,0,iiwww,ndv,
     .      Vmat(1,1,j1,ii),Wmat(1,j,ii),jjfail)
         if(jjfail .gt. 0) return
      endif
c
c: Copy desired reflection coefficient and its derivatives from last Vmat:
      do j=1,ndv
         R(j)=Vmat(j,1,j1,ii)
      enddo
c
      return
      end
