      subroutine rp_flay(R,e11,e12,e21,e22,ze,gami1,gami2,isoe,h,
     .   rhorat,iimm,iiw,ndv,V,Wflu,jjfail)
c
c: Finds the reflection coefficient just above a fluid-fluid interface
c: (characterized by density ratio rho1/rho2=rhorat) given the 
c: reflection coefficient at the bottom of layer 2 below and the 
c: propagator matrix E for the fluid layer 2.  gami1 and gami2 are i*(the
c: vertical wavenumbers) at the bottom of layer 1 and bottom of layer 2
c: respectively.  For iiw=1, the transmission coefficient W at the bottom
c: of layer 2 is also computed.  Cases where the thickness of later 2 is
c: zero, the density ratio is one, and layer 2 is a halfspace (R=0) are 
c: treated.
c: Note that ze is exponential factor common to E.
c
      implicit complex*16(a-z)
      integer*4 isoe,ndv,j,iiw,iimm,jjfail
      complex*16 e11(3),e12(3),e21(3),e22(3),R(3),V(3),Wflu(6),
     .   gami1(3),gami2(3),zzero
      real*8 h,rhorat
      data zzero/(0.d0,0.d0)/
c
      if(h .eq. 0.d0) then
c: Interface only:
cxx      if(gami1(1) .eq. gami2(1) .and. rhorat .eq. 1.d0) then
         if(iimm .eq. 0) then
c: Zero thickness, no-mismatch case: V=R, W=1:
            V(1)=R(1)
            V(2)=R(2)
            V(3)=R(3)
            Wflu(1)=dcmplx(1.d0,0.d0)
            Wflu(5)=dcmplx(0.d0,0.d0)
            return
         elseif(R(1) .ne. zzero) then
            j1=-gami1(1)*(1.d0 + R(1))
            k1=gami2(1)*(1.d0 - R(1))
         else
c: Halfspace below:
            j1=-gami1(1)
            k1=gami2(1)
         endif
      else
         if(isoe .eq. 0) then
            game11=gami2(1)*e11(1)
            game12=gami2(1)*e12(1)
            e1m=game11 - e21(1)
            e1p=game11 + e21(1)
            e2m=game12 - e22(1)
            e2p=game12 + e22(1)
         else
c: Special expressions for isospeed layer:
            if(rhorat .eq. 1.d0 .and. gami1(1) .eq. gami2(1)) then
c: Special case for isospeed layer in continuous medium (see. p.124):
               num=R(1)*e11(1)
               if(cdabs(e21(1)).eq.0.d0) then
                  jjfail=1
                  return
               end if
               V(1)=num/e21(1)
               if(iiw .eq. 1) then
                  Wflu(1)=gami1(1)/e21(1)
                  Wflu(5)=-ze
               endif
               do j=2,ndv
                  nump=R(1)*e11(j) + R(j)*e11(1)
                  V(j)=(nump - V(1)*e21(j))/e21(1)
               enddo
               return
            endif
            e1p=e11(1)
            e1m=e21(1)
            e2p=e12(1)
            e2m=e22(1)
         endif
         if(R(1) .ne. zzero) then
            term2=e2m - R(1)*e2p
            j1=gami1(1)*term2
            k1=e1m - R(1)*e1p
         else
            j1=gami1(1)*e2m
            k1=e1m
         endif
      endif
c
      if(rhorat .ne. 1.d0) k1=rhorat*k1
      num=j1 + k1
      den=j1 - k1
cpln Added 03/02-2000 for soft layers
cpln taken from PROSIM in rx_rp_flay
      if(cdabs(den). eq. 0.d0) then
         V(1)=(0.d0,0.d0) 
         Wflu(1)=(0.d0,0.d0) 
         do j=2,ndv
            V(j)=(1.d0,0.d0) 
            Wflu(j)=(1.d0,0.d0) 
         enddo
         print *,'Info msg: D=0 in rp_flay'
         jjfail=1
         return
      endif
      V(1)=num/den
      if(iiw .eq. 1) then
         Wflu(1)=-2.d0*gami1(1)/den
         Wflu(5)=-ze
         if(rhorat .ne. 1.d0) Wflu(1)=rhorat*Wflu(1)
      endif
c
c: If first derivatives desired, compute them now:
      do j=2,ndv
         if(h .eq. 0.d0) then
c: Interface only:
            if(R(1) .ne. zzero) then
               j1p=-gami1(1)*R(j) - gami1(j)*(1.d0 + R(1))
               k1p=-gami2(1)*R(j) + gami2(j)*(1.d0 - R(1))
            else
c: Halfspace below:
               j1p=-gami1(j)
               k1p=gami2(j)
            endif
         else
            if(isoe .eq. 0) then
               game11p=gami2(1)*e11(j) + gami2(j)*e11(1)
               game12p=gami2(1)*e12(j) + gami2(j)*e12(1)
               e1mp=game11p - e21(j)
               e1pp=game11p + e21(j)
               e2mp=game12p - e22(j)
               e2pp=game12p + e22(j)
            else
c: Special expressions for isospeed layer:
               e1pp=e11(j)
               e1mp=e21(j)
               e2pp=e12(j)
               e2mp=e22(j)
            endif
            if(R(1) .ne. zzero) then
               term2p=e2mp - R(1)*e2pp - R(j)*e2p
               j1p=gami1(1)*term2p + gami1(j)*term2
               k1p=e1mp - R(1)*e1pp - R(j)*e1p
            else
               j1p=gami1(1)*e2mp + gami1(j)*e2m
               k1p=e1mp
            endif
         endif
c
         if(rhorat .ne. 1.d0) k1p=rhorat*k1p
         nump=j1p + k1p
         denp=j1p - k1p
         V(j)=(den*nump - denp*num)/(den*den)
      enddo
c
      return
      end
